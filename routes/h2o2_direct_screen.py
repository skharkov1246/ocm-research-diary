#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
MatterForge · h2o2-direct — Stage 2 ПРОДАКШЕН-СКРИН дешёвых 3d-центров (M–N₄ SAC).

Переносит «дескриптор Pd₁» на земные металлы: для M ∈ {Fe, Co, Ni, Pd} строит
металло-порфин (канонический M–N₄ карман), находит основное спиновое состояние,
оптимизирует геометрию, сажает *OOH, и считает дескриптор селективности 2e⁻→H₂O₂:
  • держит ли центр связь O–O (R(O–O) в *OOH; ~1.45 Å = пероксо/2e⁻, удлинение = к 4e⁻);
  • энергию связывания ΔE(*OOH);
  • мультиреференсность центра — N_u по NOON из AVAS(M 3d + O 2p)→CASCI (как в NiO/OCM).
Ключевой литературный контраст для проверки: Fe–N₄ тянет в 4e⁻ (→H₂O), Co–N₄ — в 2e⁻.

ВАЖНО (с оговоркой): открытооболочечный DFT металло-порфина (def2-SVP, ~420 базисных фн)
локально на CPU не сходится за разумное время (SCF-стена; проверено 3 попытки: plain,
density-fit, second-order Newton). Это AWS-нагрузка, как FeMoco / NiO-18q-VQE в др. вкладках.

env-ручки (для AWS-сплита по инстансам, см. routes/h2o2_screen_aws.py):
  H2O2_METALS        подмножество из Fe,Co,Ni,Pd через запятую (default все 4)
  H2O2_SUFFIX        суффикс results-JSON, напр. _feco → h2o2_direct_screen_results_feco.json
  H2O2_STAGE_TIMEOUT таймаут одной стадии, сек (default 5400; SIGALRM, честно пишется в JSON)
Чекпойнты: results-JSON дописывается после КАЖДОЙ стадии каждого металла; повторный
запуск резюмирует с места обрыва. RDKit нужен только для build_mp — при кешированной
в JSON геометрии не требуется.

  python3 routes/h2o2_direct_screen.py --selftest   # быстрая проверка ПЛУМБИНГА (мелкие системы)
  python3 routes/h2o2_direct_screen.py --screen      # полный скрин (для AWS/мощной машины)
"""
import json
import os
import signal
import sys
import time

import numpy as np

SMI = "c1cc2cc3ccc(cc4ccc(cc5ccc(cc1[nH]2)n5)[nH]4)n3"   # free-base porphine, C20H14N4
DLABEL = {"Fe": "Fe 3d", "Co": "Co 3d", "Ni": "Ni 3d", "Pd": "Pd 4d"}
# плауз. мультиплетности M(II)-порфина (2S+1) до *OOH (метрика — основное состояние)
SPINS = {"Fe": [1, 3, 5], "Co": [2, 4], "Ni": [1, 3], "Pd": [1, 3]}

SUFFIX = os.environ.get("H2O2_SUFFIX", "")
RESULTS_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                            f"h2o2_direct_screen_results{SUFFIX}.json")
STAGE_TIMEOUT = int(os.environ.get("H2O2_STAGE_TIMEOUT", "5400"))
# Функционал энергий/дескрипторов и ОТДЕЛЬНО функционал релаксации геометрии.
# Гибридные градиенты на открытооболочечном аддукте *OOH разваливались:
# Fe/Co — RuntimeError, Pd — таймаут (прогон 26.07). Тот же урок, что в
# spin_field_oxidation v2: релакс на PBE, гибрид — single-point (PBE0//PBE).
XC = os.environ.get("H2O2_XC", "pbe0")
RELAX_XC = os.environ.get("H2O2_RELAX_XC", "pbe")


def env_metals():
    raw = os.environ.get("H2O2_METALS", "Fe,Co,Ni,Pd")
    ms = [t.strip().capitalize() for t in raw.split(",") if t.strip()]
    bad = [m for m in ms if m not in SPINS]
    if bad:
        sys.exit(f"H2O2_METALS: неизвестные металлы {bad}; можно: {sorted(SPINS)}")
    return ms


# ---------- чекпойнты (поточечные: JSON дописывается после каждой стадии) ----------

def jat(atoms):
    """atoms → JSON-сериализуемое [[symbol, [x,y,z]], …]."""
    return [[s, [float(x) for x in c]] for s, c in atoms]


def unjat(rows):
    return [(s, tuple(c)) for s, c in rows]


def load_results(path=RESULTS_PATH):
    if os.path.exists(path):
        with open(path) as f:
            return json.load(f)
    return {"meta": {"script": "h2o2_direct_screen", "suffix": SUFFIX,
                     "stage_timeout_s": STAGE_TIMEOUT,
                     "basis": "def2-svp", "xc": XC, "relax_xc": RELAX_XC,
                     "protocol": (f"{XC}//{RELAX_XC}" if RELAX_XC != XC else XC),
                     "protocol_note": ("геометрия аддукта *OOH релаксируется на "
                                       "RELAX_XC, энергии/дескрипторы — на XC; "
                                       "гибридные градиенты открытой оболочки "
                                       "не переживали оптимизацию (прогон 26.07)")},
            "metals": {}}


def save_results(res, path=RESULTS_PATH):
    tmp = path + ".tmp"
    with open(tmp, "w") as f:
        json.dump(res, f, indent=1, ensure_ascii=False)
    os.replace(tmp, path)


# ---------- таймаут стадии (SIGALRM; статус честно пишется в JSON) ----------

class StageTimeout(Exception):
    pass


def with_timeout(fn, seconds):
    def _alrm(signum, frame):
        raise StageTimeout()
    old = signal.signal(signal.SIGALRM, _alrm)
    signal.alarm(int(seconds))
    try:
        return fn()
    finally:
        signal.alarm(0)
        signal.signal(signal.SIGALRM, old)


def run_stage(res, rec, name, fn):
    """одна стадия под SIGALRM-таймаутом; статус/время — в JSON сразу же.
    Возвращает результат fn() либо None (timeout/исключение — записаны честно)."""
    st = rec.setdefault("stages", {})
    t0 = time.time()
    try:
        out = with_timeout(fn, STAGE_TIMEOUT)
        st[name] = {"status": "ok", "wall_s": round(time.time() - t0, 1)}
        return out
    except StageTimeout:
        st[name] = {"status": "timeout", "timeout_s": STAGE_TIMEOUT,
                    "wall_s": round(time.time() - t0, 1)}
        print(f"  [stage {name}] TIMEOUT {STAGE_TIMEOUT}s — записано в JSON")
        return None
    except Exception as e:
        st[name] = {"status": f"fail:{type(e).__name__}", "error": str(e)[:300],
                    "wall_s": round(time.time() - t0, 1)}
        print(f"  [stage {name}] FAIL {type(e).__name__}: {e}")
        return None
    finally:
        save_results(res)


def stage_ok(rec, name):
    return rec.get("stages", {}).get(name, {}).get("status") == "ok"


# ---------- научное содержание (дескрипторы/уровни/спины — БЕЗ изменений) ----------

def Nu(occ):
    occ = np.asarray(occ)
    return float(np.sum(occ * (2.0 - occ)))


def build_mp(metal):
    """RDKit: порфин → 3D → металлирование (металл в центроид N₄, снять 2 пиррольных N–H).
    RDKit нужен ТОЛЬКО здесь; при кешированной в чекпойнте геометрии не вызывается."""
    from rdkit import Chem
    from rdkit.Chem import AllChem
    m = Chem.AddHs(Chem.MolFromSmiles(SMI))
    AllChem.EmbedMolecule(m, AllChem.ETKDGv3())
    AllChem.MMFFOptimizeMolecule(m, maxIters=2000)
    P = m.GetConformer().GetPositions()
    Ns = [a.GetIdx() for a in m.GetAtoms() if a.GetSymbol() == "N"]
    drop = {nb.GetIdx() for n in Ns for nb in m.GetAtomWithIdx(n).GetNeighbors()
            if nb.GetSymbol() == "H"}
    cen = P[Ns].mean(0)
    atoms = [(metal, tuple(cen))]
    atoms += [(a.GetSymbol(), tuple(P[a.GetIdx()]))
              for a in m.GetAtoms() if a.GetIdx() not in drop]
    return atoms


def n4_normal(atoms):
    """нормаль к плоскости N₄ (через SVD)."""
    N = np.array([c for s, c in atoms if s == "N"])
    c = N.mean(0)
    _, _, vt = np.linalg.svd(N - c)
    return vt[-1] / np.linalg.norm(vt[-1])


def place_ooh(atoms, m_o=1.85, o_o=1.40, o_h=0.97):
    """посадить *OOH аксиально на металл (геом-опт потом уточнит)."""
    metal = np.array(atoms[0][1])
    n = n4_normal(atoms)
    t = np.cross(n, [1.0, 0.0, 0.0]); t = t / (np.linalg.norm(t) + 1e-9)
    O1 = metal + m_o * n
    O2 = O1 + o_o * (0.5 * n + 0.866 * t)
    H = O2 + o_h * (0.5 * n - 0.866 * t)
    return atoms + [("O", tuple(O1)), ("O", tuple(O2)), ("H", tuple(H))]


def pyscf_atom(atoms):
    return "; ".join(f"{s} {x:.5f} {y:.5f} {z:.5f}" for s, (x, y, z) in atoms)


def make_mol(atoms, spin, charge=0, basis="def2-svp"):
    from pyscf import gto
    return gto.M(atom=pyscf_atom(atoms), basis=basis, ecp=basis, spin=spin,
                 charge=charge, verbose=0)


def robust_uks(mol, xc=None):
    """устойчивый открытооболочечный SCF: density-fit + level-shift + second-order Newton."""
    xc = XC if xc is None else xc
    from pyscf import dft
    mf = dft.UKS(mol, xc=xc).density_fit()
    mf.grids.level = 3
    mf.level_shift = 0.3
    mf.conv_tol = 1e-7
    mf = mf.newton()
    mf.max_cycle = 120
    mf.kernel()
    return mf


def avas_casci_Nu(mf, metal):
    """AVAS(M 3d + O 2p) → CASCI → N_u по заселённостям натуральных орбиталей."""
    from pyscf import mcscf
    from pyscf.mcscf import avas
    ncas, nelec, mo = avas.avas(mf, [DLABEL[metal], "O 2p"], canonicalize=True, verbose=0)
    mc = mcscf.CASCI(mf, ncas, nelec)
    mc.kernel(mo)
    dm1 = mc.fcisolver.make_rdm1(mc.ci, mc.ncas, mc.nelecas)
    occ = np.linalg.eigvalsh(dm1)[::-1]
    return ncas, nelec, Nu(occ), occ, mc.e_tot


def spins_ooh(metal):
    """мультиплетности для M-порфин+*OOH: аддукт с дублетным •OOH меняет ЧЁТНОСТЬ
    числа электронов — берём «голые ±1» (ферро/антиферро связка радикала с центром).
    Исходный скрипт передавал сюда голые SPINS -> pyscf мгновенно отбрасывал все
    спины (parity mismatch), стадия падала за 0.0 с — пойман AWS-прогоном 2026-07-25."""
    return sorted({s + d for s in SPINS[metal] for d in (-1, +1) if s + d >= 1})


def ground_spin(metal, atoms, spins=None):
    """скан мультиплетности → основное состояние (E, spin, mf)."""
    best = None
    for s in (spins if spins is not None else SPINS[metal]):
        try:
            mf = robust_uks(make_mol(atoms, spin=s - 1))   # spin = 2S = (2S+1)-1
            if mf.converged and (best is None or mf.e_tot < best[0]):
                best = (mf.e_tot, s, mf)
        except Exception as e:
            print(f"   spin {s}: FAIL {type(e).__name__}")
    return best


def ground_spin_or_raise(metal, atoms, spins=None):
    """ground_spin, но «ни один спин не сошёлся» — это fail стадии (честно в JSON)."""
    gs = ground_spin(metal, atoms, spins=spins)
    if gs is None:
        raise RuntimeError("SCF не сошёлся ни в одном спине")
    return gs


def uks_at(atoms, spin_mult, xc=None):
    """пересчёт SCF в уже известном (кешированном) спине — для resume между стадиями."""
    mf = robust_uks(make_mol(atoms, spin=spin_mult - 1), xc=xc)
    if not mf.converged:
        raise RuntimeError("SCF not converged on resume")
    return mf


def optimize_adsorbate(mf, n_atoms_total, maxsteps=80):
    """Релакс ТОЛЬКО адсорбата *OOH (последние 3 атома); каркас M–N₄ заморожен.

    Полная оптимизация разваливалась на «trust radius got too small» (Pd — и на
    PBE, и на PBE0): виноват не функционал, а оптимизатор на 40 степенях свободы
    жёсткого макроцикла. Дескриптору нужна только геометрия аддукта, поэтому
    каркас фиксируется — это и устойчивее, и дешевле. Фоллбэк на berny сохранён.
    """
    import tempfile
    from pyscf.geomopt import geometric_solver
    with tempfile.NamedTemporaryFile("w", suffix=".txt", delete=False) as fh:
        fh.write(f"$freeze\nxyz 1-{n_atoms_total - 3}\n")
        cons = fh.name
    return geometric_solver.optimize(mf, constraints=cons, maxsteps=maxsteps)


def mole_atoms(mol):
    """геометрия pyscf Mole → наш формат atoms (Å)."""
    return [(s, tuple(c)) for s, c in zip(mol.elements, mol.atom_coords() * 0.529177)]


# ---------- скрин с чекпойнтами и resume ----------

def screen_metal(M, res):
    from pyscf.geomopt.berny_solver import optimize
    rec = res["metals"].setdefault(M, {})
    # Флаг done ставится по завершении всех стадий. Если стадии были СБРОШЕНЫ
    # (смена протокола, пересчёт под единый уровень), флаг обязан уступить:
    # иначе металл молча пропускается и остаётся со старыми числами.
    if rec.get("done") and stage_ok(rec, "opt_ooh"):
        print(f"\n### {M}-порфин ### — уже сделан (resume), пропуск")
        return
    if rec.get("done"):
        print(f"\n### {M}-порфин ### — флаг done снят: стадии сброшены, пересчёт")
        rec.pop("done", None)
    print(f"\n### {M}-порфин ###")

    # 1) геометрия (RDKit только здесь; кеш в JSON снимает зависимость от rdkit)
    if rec.get("atoms_bare"):
        bare = unjat(rec["atoms_bare"])
        print(f"  геометрия из чекпойнта ({len(bare)} атомов), RDKit не нужен")
    else:
        bare = run_stage(res, rec, "build", lambda: build_mp(M))
        if bare is None:
            return
        rec["atoms_bare"] = jat(bare)
        save_results(res)

    # 2) основной спин голого центра
    mf = None
    if stage_ok(rec, "spin_bare"):
        spin = rec["bare"]["spin"]
        print(f"  spin_bare из чекпойнта: (2S+1)={spin} E={rec['bare']['E']:.4f}")
    else:
        gs = run_stage(res, rec, "spin_bare",
                       lambda: ground_spin_or_raise(M, bare))
        if gs is None:
            return
        E_bare, spin, mf = gs
        rec["bare"] = {"E": float(E_bare), "spin": int(spin)}
        save_results(res)
        print(f"  основной спин (2S+1)={spin}  E_bare={E_bare:.4f}")

    # 3) геом-опт голого центра
    if stage_ok(rec, "opt_bare"):
        bare_opt = unjat(rec["atoms_bare_opt"])
        print("  opt_bare из чекпойнта")
    else:
        if mf is None:   # resume: SCF заново в кешированном спине
            mf = run_stage(res, rec, "rescf_bare", lambda: uks_at(bare, spin))
            if mf is None:
                return
        mf_bare = mf
        mol_opt = run_stage(res, rec, "opt_bare",
                            lambda: optimize(mf_bare, maxsteps=60))
        if mol_opt is None:
            return
        bare_opt = mole_atoms(mol_opt)
        rec["atoms_bare_opt"] = jat(bare_opt)
        save_results(res)

    # 4) *OOH: посадка + основной спин
    ads = place_ooh(bare_opt)
    mf_a = None
    if stage_ok(rec, "spin_ooh"):
        spin_a = rec["ooh"]["spin"]
        print(f"  spin_ooh из чекпойнта: (2S+1)={spin_a} E={rec['ooh']['E']:.4f}")
    else:
        gsa = run_stage(res, rec, "spin_ooh",
                        lambda: ground_spin_or_raise(M, ads, spins=spins_ooh(M)))
        if gsa is None:
            return
        E_ads, spin_a, mf_a = gsa
        rec["ooh"] = {"E": float(E_ads), "spin": int(spin_a)}
        save_results(res)

    # 5) геом-опт *OOH + R(O–O)
    if stage_ok(rec, "opt_ooh"):
        roo = rec["ooh"]["roo"]
        print(f"  opt_ooh из чекпойнта: R(O–O)={roo:.3f}Å")
    else:
        if mf_a is None:
            mf_a = run_stage(res, rec, "rescf_ooh", lambda: uks_at(ads, spin_a))
            if mf_a is None:
                return
        # Релакс аддукта идёт на RELAX_XC (по умолчанию PBE): на гибриде
        # градиенты открытооболочечного *OOH не переживали оптимизацию.
        if RELAX_XC == XC:
            mfa = mf_a
        else:
            mfa = run_stage(res, rec, "rescf_ooh_relax",
                            lambda: uks_at(ads, spin_a, xc=RELAX_XC))
            if mfa is None:
                return
        def _opt_ooh():
            try:
                return optimize_adsorbate(mfa, len(ads), maxsteps=80)
            except Exception as e:                      # geomeTRIC нет/сорвался
                print(f"  [opt_ooh] констрейнд-опт не пошёл ({type(e).__name__}: "
                      f"{e}); фоллбэк на berny по всей геометрии")
                rec["opt_ooh_fallback"] = f"{type(e).__name__}: {e}"[:200]
                return optimize(mfa, maxsteps=80)

        mol_a = run_stage(res, rec, "opt_ooh", _opt_ooh)
        if mol_a is None:
            return
        rec["opt_ooh_mode"] = ("frozen-frame/geomeTRIC"
                               if "opt_ooh_fallback" not in rec else "berny-full")
        rec["relax_protocol"] = (f"{XC}//{RELAX_XC}" if RELAX_XC != XC else XC)
        ooh_opt = mole_atoms(mol_a)
        rec["atoms_ooh_opt"] = jat(ooh_opt)
        coords = np.array([c for _s, c in ooh_opt])
        # O–O расстояние (последние два O в списке)
        oidx = [i for i, (z, _c) in enumerate(ooh_opt) if z == "O"][-2:]
        roo = float(np.linalg.norm(coords[oidx[0]] - coords[oidx[1]]))
        rec["ooh"]["roo"] = roo
        save_results(res)

    # 6) AVAS→CASCI→N_u (как в оригинале — на mf основного спина *OOH, до-опт геометрия)
    if stage_ok(rec, "cas"):
        cas = rec["cas"]
        print(f"  cas из чекпойнта: CAS({cas['nelec']},{cas['ncas']}o) N_u={cas['nu']:.2f}")
    else:
        if mf_a is None:
            mf_a = run_stage(res, rec, "rescf_ooh", lambda: uks_at(ads, spin_a))
            if mf_a is None:
                return
        mfa2 = mf_a
        out = run_stage(res, rec, "cas", lambda: avas_casci_Nu(mfa2, M))
        if out is None:
            return
        ncas, nelec, nu, occ, e_cas = out
        rec["cas"] = {"ncas": int(ncas), "nelec": int(nelec), "nu": float(nu),
                      "noon": [float(o) for o in occ], "e_cas": float(e_cas)}
        save_results(res)
        print(f"  *OOH: спин={spin_a} R(O–O)={roo:.3f}Å  "
              f"AVAS->CAS({nelec},{ncas}o) N_u={nu:.2f}")

    rec["holds_oo_2e"] = bool(roo < 1.55)
    rec["done"] = True
    save_results(res)
    print(f"  => держит O–O (2e⁻): {'да' if roo < 1.55 else 'НЕТ (к 4e⁻)'}")


def screen(metals=None):
    metals = metals or env_metals()
    res = load_results()
    print(f"== Stage 2 продакшен-скрин M–N₄ SAC (тяжёлое; для AWS) ==")
    print(f"metals={','.join(metals)} suffix='{SUFFIX}' "
          f"stage_timeout={STAGE_TIMEOUT}s results={RESULTS_PATH}")
    for M in metals:
        screen_metal(M, res)
    save_results(res)
    print(f"\n[checkpoint] итог в {RESULTS_PATH}")


def selftest():
    """быстрая проверка плумбинга (без металло-порфиновой SCF-стены);
    rdkit локально НЕ обязателен — при его отсутствии build_mp подменяется мок-N₄."""
    import tempfile
    print("== SELFTEST: плумбинг скрина на мелких системах ==")
    try:
        fe = build_mp("Fe")
        print(f"[build] Fe-порфин: {len(fe)} атомов (ожид. 37)")
    except ImportError:
        fe = [("Fe", (0.0, 0.0, 0.0)), ("N", (2.0, 0.0, 0.0)), ("N", (-2.0, 0.0, 0.0)),
              ("N", (0.0, 2.0, 0.0)), ("N", (0.0, -2.0, 0.0))]
        print("[build] rdkit локально недоступен — build_mp пропущен, мок-N₄ "
              "(на AWS rdkit ставится pip-ом; при кеше геометрии в JSON не нужен)")
    ads = place_ooh(fe)
    roo = np.linalg.norm(np.array(ads[-2][1]) - np.array(ads[-3][1]))
    print(f"[place_ooh] стартовое R(O–O)={roo:.3f}Å (до опт)")
    # чекпойнт: json round-trip геометрии + save/load + resume-маркер
    with tempfile.TemporaryDirectory() as td:
        p = os.path.join(td, "ckpt.json")
        res = {"meta": {}, "metals": {"Fe": {"atoms_bare": jat(ads),
                                             "stages": {"build": {"status": "ok"}}}}}
        save_results(res, p)
        back = load_results(p)
        assert unjat(back["metals"]["Fe"]["atoms_bare"]) == [
            (s, tuple(c)) for s, c in ads]
        assert stage_ok(back["metals"]["Fe"], "build")
        print(f"[checkpoint] save/load/resume-маркер OK ({p})")
    # SIGALRM-таймаут стадии
    try:
        with_timeout(lambda: time.sleep(3), 1)
        print("[timeout] FAIL — SIGALRM не сработал"); sys.exit(1)
    except StageTimeout:
        print("[timeout] SIGALRM-стадийный таймаут работает (1s на sleep 3s)")
    # robust_uks на маленьком открытооболочечном O2 (триплет) — быстро
    mf = robust_uks(make_mol([("O", (0, 0, 0)), ("O", (0, 0, 1.207))], spin=2))
    print(f"[robust_uks] O2 триплет: conv={mf.converged} E={mf.e_tot:.4f}")
    # AVAS+CASCI+N_u код-путь на O2 (метка 'O 2p')
    from pyscf import mcscf
    from pyscf.mcscf import avas
    ncas, nelec, mo = avas.avas(mf, ["O 2p"], canonicalize=True, verbose=0)
    mc = mcscf.CASCI(mf, ncas, nelec); mc.kernel(mo)
    dm1 = mc.fcisolver.make_rdm1(mc.ci, mc.ncas, mc.nelecas)
    occ = np.linalg.eigvalsh(dm1)[::-1]
    print(f"[avas_casci_Nu] O2: CAS({nelec},{ncas}o) N_u={Nu(occ):.2f} (ожид. ~2.3 — бирадикал)")
    print("SELFTEST OK — плумбинг корректен; продакшен-скрин уезжает на AWS.")


if __name__ == "__main__":
    if "--screen" in sys.argv:
        screen()
    else:
        selftest()
