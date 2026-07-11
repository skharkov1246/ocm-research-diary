#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
routes/lao_cr_trimer.py — вариант 2 гребёнки (ЛАО/ПАО-скид): дескриптор
селективной ТРИмеризации этилена на Cr. Этилен со скида не возится (−104°C);
олигомеризация на месте в 1-гексен/ЛАО→ПАО ($2500–4000/т, жидкость) снимает
логистику. Квант-ниша: цикл Cr(I)/Cr(III) — метациклы со спин-кроссовером,
редокс-пара, где DFT исторически путается (известная проблема литературы
триммеризации Phillips/Sasol).

СМОК-СТАДИЯ (эта): жизнеспособность оснастки, без TS —
  1) spins: лестницы UKS-PBE0 для хромациклопентана и хромациклогептана
     [Cl2Cr(C4H8)] / [Cl2Cr(C6H12)] (модель-скелет, как bare-этапы OCM);
  2) int: релакс интермедиатов + CASSCF/NEVPT2 одноточки (AVAS [Cr 3d]) —
     сходимость и NOON-мультиреференсность (есть ли вообще наша ниша).
Дескриптор (следующая стадия, по смоку): ΔΔE‡ = барьер(внедрение C2H4 в
C7-цикл → рост/C8+) − барьер(β-H элиминирование → 1-гексен); >0 = селективно
к C6. DFT vs NEVPT2 — переворачивает ли корреляция C6/C8-чтение.

Запуск: python3 routes/lao_cr_trimer.py spins|int|sanity
Env: LAO_BASIS (def2-svp).
"""
import json
import math
import os
import sys
import time

import numpy as np
from pyscf import dft, fci, gto, mcscf, scf
from pyscf.mcscf import avas
from pyscf.mrpt import NEVPT

HARTREE_KCAL = 627.509474
DIR = os.path.dirname(os.path.abspath(__file__))
BASIS = os.environ.get("LAO_BASIS", "def2-svp")
XC = "pbe0"


def M(atoms, spin, charge=0):
    return gto.M(atom=atoms, basis=BASIS, ecp=BASIS, charge=charge, spin=spin,
                 verbose=0, max_memory=8000)


def uks(mol, dm0=None):
    mf = dft.UKS(mol).density_fit()
    mf.xc = XC
    mf.conv_tol = 1e-8
    mf.max_cycle = 200
    mf.kernel(dm0=dm0) if dm0 is not None else mf.kernel()
    if not mf.converged:
        mf.level_shift = 0.3
        mf.max_cycle = 300
        mf.kernel(mf.make_rdm1())
    return mf


def rohf(mol):
    mf = scf.ROHF(mol)
    mf.conv_tol = 1e-8
    mf.max_cycle = 200
    mf.kernel()
    if not mf.converged:
        mf.level_shift = 0.4
        mf.max_cycle = 300
        mf.kernel(mf.make_rdm1())
    return mf


# ---------------------------------------------------------------- геометрии
def _ring_carbons(n, r0=2.05, ccb=1.53):
    """Позиции C метацикла Cr(CH2)n: Cr в нуле, кольцо в плоскости xz дугой.
    Грубый старт — релакс доводит."""
    pts = []
    span = math.radians(100 + 12 * n)          # раствор дуги под размер кольца
    for k in range(n):
        th = -span / 2 + span * k / (n - 1)
        rad = r0 + (ccb - r0 * 0.35) * min(k, n - 1 - k) * 0.55
        pts.append((rad * math.sin(th), 0.0, rad * math.cos(th)))
    return pts


def chromacycle_atoms(n_c):
    """[Cl2Cr(CnH2n)] метацикл: Cr + 2 Cl (тетраэдрично сверху) + кольцо CH2."""
    a = [("Cr", (0.0, 0.0, 0.0)),
         ("Cl", (1.65, 1.45, -0.85)), ("Cl", (-1.65, 1.45, -0.85))]
    cs = _ring_carbons(n_c)
    for x, y, z in cs:
        a.append(("C", (x, y, z)))
    d = 1.09
    for i, (x, y, z) in enumerate(cs):
        # два H на каждый C, разнесены из плоскости кольца
        a.append(("H", (x * 1.12, y + d * 0.82, z * 1.12)))
        a.append(("H", (x * 1.12, y - d * 0.82, z * 1.12)))
    return a


def relax(atoms, spin, maxsteps=150):
    from pyscf.geomopt.geometric_solver import optimize
    mf = uks(M(atoms, spin))
    mol = optimize(mf, maxsteps=maxsteps, assert_convergence=False)
    na = [(mol.atom_symbol(i), tuple(c))
          for i, c in enumerate(mol.atom_coords(unit="Angstrom"))]
    mfx = uks(mol, dm0=mf.make_rdm1())
    return na, mfx


def nevpt2_point(atoms, spin, thr=0.5, labels=None):
    mol = M(atoms, spin)
    mol.max_memory = 24000     # c7 (211 AO) уронил NEVPT2 по MemoryError на 8 ГБ
    mf = rohf(mol)
    ncas, nelec, mo = avas.avas(mf, labels or ["Cr 3d"], threshold=thr)
    ss = spin / 2 * (spin / 2 + 1)
    mc = mcscf.CASSCF(mf, ncas, nelec)
    fci.addons.fix_spin_(mc.fcisolver, shift=0.2, ss=ss)
    mc.max_cycle_macro = 150
    mc.conv_tol = 1e-7
    mc.kernel(mo)
    dm1 = mc.fcisolver.make_rdm1(mc.ci, mc.ncas, mc.nelecas)
    noon = sorted(np.linalg.eigvalsh(dm1), reverse=True)
    return {"e_casscf": float(mc.e_tot),
            "e_nevpt2": float(mc.e_tot + NEVPT(mc).kernel()),
            "cas": [int(np.sum(nelec)), int(ncas)],
            "conv": bool(mc.converged),
            "noon": [round(float(x), 3) for x in noon]}


RINGS = {"c5": 4, "c7": 6}     # хромациклопентан (4 CH2), -гептан (6 CH2)


def stage_spins():
    """Лестницы спинов обоих циклов. Нейтральный Cl2Cr(кольцо) = формально
    Cr(IV) d2 (кольцо −2, 2 Cl −1) — чётный счёт: кандидаты по чётности."""
    out = {}
    path = os.path.join(DIR, "lao_cr_spins.json")
    for tag, n in RINGS.items():
        atoms = chromacycle_atoms(n)
        parity = sum(gto.charge(a[0]) for a in atoms) % 2
        cands = [s for s in range(parity, parity + 6, 2)][:3]
        out[tag] = {"ladder": [], "candidates_2S": cands}
        for s2 in cands:
            t0 = time.time()
            try:
                mf = uks(M(atoms, s2))
                out[tag]["ladder"].append(
                    {"spin_2S": s2, "e_h": float(mf.e_tot),
                     "converged": bool(mf.converged),
                     "spin_square": round(float(mf.spin_square()[0]), 3),
                     "wall_s": round(time.time() - t0, 1)})
            except Exception as e:
                out[tag]["ladder"].append({"spin_2S": s2, "error": str(e)[:150]})
            print(f"[spins][{tag}] 2S={s2}: {out[tag]['ladder'][-1]}", flush=True)
            with open(path, "w") as f:
                json.dump(out, f, indent=1)
        ok = [x for x in out[tag]["ladder"] if x.get("converged")]
        if ok:
            best = min(ok, key=lambda x: x["e_h"])
            out[tag]["ground_2S"] = best["spin_2S"]
            out[tag]["rel_kcal"] = {
                str(x["spin_2S"]): round((x["e_h"] - best["e_h"]) * HARTREE_KCAL, 2)
                for x in ok}
            print(f"[spins][{tag}] ground 2S={best['spin_2S']} "
                  f"ladder={out[tag]['rel_kcal']}", flush=True)
        with open(path, "w") as f:
            json.dump(out, f, indent=1)


def stage_int():
    """Релакс интермедиатов в основном спине + CASSCF/NEVPT2-одноточка:
    сходимость + NOON (мультиреференсность = есть ли ниша)."""
    sp = json.load(open(os.path.join(DIR, "lao_cr_spins.json")))
    path = os.path.join(DIR, "lao_cr_int.json")
    out = {}
    if os.path.exists(path):
        out = json.load(open(path))
    for tag, n in RINGS.items():
        if tag in out and out[tag].get("nevpt2"):
            print(f"[int][{tag}] resume: on disk", flush=True)
            continue
        s2 = sp[tag]["ground_2S"]
        t0 = time.time()
        atoms, mf = relax(chromacycle_atoms(n), s2)
        rec = {"spin_2S": s2, "e_dft": float(mf.e_tot),
               "dft_converged": bool(mf.converged),
               "atoms": [[s, [round(x, 6) for x in c]] for s, c in atoms]}
        try:
            rec.update({"nevpt2": nevpt2_point(atoms, s2)})
        except Exception as e:
            rec["nevpt2_error"] = str(e)[:200]
        rec["wall_s"] = round(time.time() - t0, 1)
        out[tag] = rec
        with open(path, "w") as f:
            json.dump(out, f, indent=1)
        print(f"[int][{tag}] 2S={s2} dft={rec['e_dft']:.6f} "
              f"nevpt2={rec.get('nevpt2')}", flush=True)
    # смок-вердикт: обе точки сошлись и NOON мультиреференсны?
    ver = {}
    for tag in RINGS:
        nv = out.get(tag, {}).get("nevpt2") or {}
        noon = nv.get("noon") or []
        frac = sum(1 for x in noon if 0.1 < x < 1.9)
        ver[tag] = {"casscf_conv": nv.get("conv"), "n_frac_noon": frac}
    out["_smoke_verdict"] = ver
    with open(path, "w") as f:
        json.dump(out, f, indent=1)
    print("[int] smoke verdict:", json.dumps(ver), flush=True)


def stage_sanity():
    for tag, n in RINGS.items():
        atoms = chromacycle_atoms(n)
        try:
            mol = M(atoms, sum(__import__("pyscf").gto.charge(a[0]) for a in atoms) % 2)
            print(f"[sanity] {tag}: natoms={len(atoms)} nao={mol.nao} "
                  f"nelec={mol.nelectron} build OK", flush=True)
        except Exception as e:
            print(f"[sanity] {tag} FAIL: {e}", flush=True)


if __name__ == "__main__":
    st = sys.argv[1] if len(sys.argv) > 1 else "sanity"
    {"spins": stage_spins, "int": stage_int, "sanity": stage_sanity,
     "bhe": stage_bhe, "ins": stage_ins, "desc": stage_desc}[st]()


# ---------------------------------------------------------- TS-стадии (v2)
# Дескриптор: ddE_LAO = E‡(внедрение C2H4 в c7) − E‡(β-H элиминирование из c7);
# >0 — элиминирование выигрывает => селективный 1-гексен (ЛАО/ПАО-скид);
# <0 — рост цепи (C8+/полимер). DFT vs NEVPT2 — флипает ли корреляция.
# Атомный порядок c7 (chromacycle_atoms(6), сохраняется релаксом):
# Cr=0, Cl=1,2, кольцевые C=3..8 (C3/C8 — альфа), H: C3->9,10; C4->11,12; ...
# β-H для элиминирования: H idx 11 (на C4, соседнем с альфа-C3).

def _c7_relaxed():
    d = json.load(open(os.path.join(DIR, "lao_cr_int.json")))
    rec = d["c7"]
    return [(s, tuple(c)) for s, c in rec["atoms"]], rec["spin_2S"]


def _scan(tag, build_start, pins, ref_e, spin, pin_atoms_1based):
    """Общий пин-скан с поточечным сохранением/resume (урок msa v6)."""
    from pyscf.geomopt.geometric_solver import optimize
    path = os.path.join(DIR, f"lao_cr_{tag}.json")
    pts = []
    if os.path.exists(path):
        try:
            pts = json.load(open(path))["pts"]
            print(f"[{tag}] resume: {len(pts)} точек", flush=True)
        except Exception:
            pts = []
    done = {round(p["pin"], 3) for p in pts}
    dm, start = None, None
    for pin in pins:
        if round(pin, 3) in done:
            saved = next(p for p in pts if round(p["pin"], 3) == round(pin, 3))
            start = [(s, tuple(c)) for s, c in saved["atoms"]]
            continue
        cf = os.path.join(DIR, f"_lao_{tag}_{pin:.2f}.txt")
        i, j = pin_atoms_1based
        open(cf, "w").write(f"$set\ndistance {i} {j} {pin:.4f}\n")
        try:
            atoms0 = start if start is not None else build_start(pin)
            mf0 = uks(M(atoms0, spin), dm0=dm)
            mol = optimize(mf0, constraints=cf, maxsteps=150,
                           assert_convergence=False)
            mfx = uks(mol, dm0=mf0.make_rdm1())
            na = [(mol.atom_symbol(k), tuple(c)) for k, c in
                  enumerate(mol.atom_coords(unit="Angstrom"))]
            dm = mfx.make_rdm1()
            start = na
            pts.append({"pin": pin, "e_h": float(mfx.e_tot),
                        "conv": bool(mfx.converged), "atoms": na})
            with open(path, "w") as f:
                json.dump({"pts": pts}, f, indent=1)
            print(f"  [{tag}] pin={pin} rel="
                  f"{(mfx.e_tot - ref_e) * HARTREE_KCAL:.2f} "
                  f"conv={mfx.converged}", flush=True)
        except Exception as ex:
            print(f"  [{tag}] pin={pin} failed: {ex}", flush=True)
        finally:
            if os.path.exists(cf):
                os.remove(cf)
    return pts, path


def _readout(pts, ref_e):
    """Интерьерный максимум + маска выбросов (>25 над обоими соседями)."""
    pts = sorted(pts, key=lambda p: -p["pin"])     # путь: пин уменьшается
    rel = [(p["e_h"] - ref_e) * HARTREE_KCAL for p in pts]
    keep = [i for i in range(len(pts)) if pts[i]["conv"] and not all(
        rel[i] - rel[j] > 25.0 for j in (i - 1, i + 1) if 0 <= j < len(pts))]
    kr = [rel[i] for i in keep]
    if len(kr) < 4:
        return None, None, rel, False
    imax = max(range(1, len(kr) - 1), key=lambda i: kr[i])
    interior = kr[imax] >= max(kr[0], kr[-1])
    return pts[keep[imax]], round(kr[imax], 2), [round(x, 2) for x in rel], interior


def stage_bhe():
    """β-H элиминирование из хромациклогептана: пин r(Cr–H_β) стягивается."""
    c7, spin = _c7_relaxed()
    mf_ref = uks(M(c7, spin))
    ref = float(mf_ref.e_tot)
    pts, path = _scan("bhe", lambda pin: c7, (2.40, 2.10, 1.90, 1.75, 1.62),
                      ref, spin, (1, 12))
    ts, bar, rel, interior = _readout(pts, ref)
    out = {"ref_e_h": ref, "rel_kcal": rel, "interior": bool(interior),
           "dft_barrier": bar}
    if ts is not None and interior and bar is not None and 0 < bar < 100:
        out["ts_atoms"] = [[s, [round(x, 6) for x in c]] for s, c in ts["atoms"]]
        try:
            n_ts = nevpt2_point(ts["atoms"], spin)
            n_r = nevpt2_point(c7, spin)
            for lv in ("e_casscf", "e_nevpt2"):
                out[lv.replace("e_", "") + "_barrier"] = round(
                    (n_ts[lv] - n_r[lv]) * HARTREE_KCAL, 2)
        except Exception as e:
            out["nevpt2_error"] = str(e)[:200]
    else:
        out["failed"] = "TS not interior/unphysical"
    with open(os.path.join(DIR, "lao_cr_bhe_result.json"), "w") as f:
        json.dump(out, f, indent=1)
    print("[bhe]", json.dumps({k: v for k, v in out.items()
                               if k != "ts_atoms"})[:400], flush=True)


def stage_ins():
    """Внедрение C2H4 в c7: этилен у Cr, пин r(C_alpha–C_eth) стягивается."""
    import numpy as _np
    c7, spin = _c7_relaxed()
    # этилен: релакс отдельно, референс = E(c7)+E(C2H4)
    eth0 = [("C", (0.0, 0.0, 0.0)), ("C", (1.33, 0.0, 0.0)),
            ("H", (-0.56, 0.92, 0.0)), ("H", (-0.56, -0.92, 0.0)),
            ("H", (1.89, 0.92, 0.0)), ("H", (1.89, -0.92, 0.0))]
    eth, mf_eth = relax(eth0, 0, maxsteps=60)
    mf_c7 = uks(M(c7, spin))
    ref = float(mf_c7.e_tot) + float(mf_eth.e_tot)

    def build(pin):
        # открытая грань: перебор направлений вокруг Cr, берём максимально
        # свободное (наивное «прочь от Cl-Cl» утыкалось в кольцо: 0.87 A клэш;
        # перебор даёт клиренс ~2.85 A)
        cr = _np.array(c7[0][1])
        others = _np.array([x[1] for x in c7[1:]])
        best_u, best_d = None, -1.0
        for i in range(200):
            th = _np.arccos(1 - 2 * ((i * 0.618034) % 1.0))
            ph = 2 * _np.pi * ((i * 0.381966) % 1.0)
            uu = _np.array([_np.sin(th) * _np.cos(ph),
                            _np.sin(th) * _np.sin(ph), _np.cos(th)])
            dmin = _np.min(_np.linalg.norm(others - (cr + uu * 2.25), axis=1))
            if dmin > best_d:
                best_d, best_u = dmin, uu
        u = best_u
        print(f"  [ins] placement clearance {best_d:.2f} A", flush=True)
        v = _np.cross(u, [0, 0, 1.0])
        if _np.linalg.norm(v) < 0.3:
            v = _np.cross(u, [0, 1.0, 0])
        v = v / _np.linalg.norm(v)
        w = _np.cross(u, v)
        c1 = cr + u * 2.35
        atoms = list(c7)
        ec = (_np.array(eth[0][1]) + _np.array(eth[1][1])) / 2
        for s, xyz in eth:
            rel_ = _np.array(xyz) - ec
            atoms.append((s, tuple(c1 + rel_[0] * v + rel_[1] * u
                                   + rel_[2] * w)))
        return atoms

    # пин: C_alpha(кольцевой C3, 1-based 4) – C_eth(первый C этилена, 1-based 22)
    pts, path = _scan("ins", build, (2.80, 2.45, 2.20, 2.00, 1.80),
                      ref, spin, (4, 22))
    ts, bar, rel, interior = _readout(pts, ref)
    out = {"ref_e_h": ref, "rel_kcal": rel, "interior": bool(interior),
           "dft_barrier": bar}
    if ts is not None and interior and bar is not None and 0 < bar < 100:
        out["ts_atoms"] = [[s, [round(x, 6) for x in c]] for s, c in ts["atoms"]]
        try:
            n_ts = nevpt2_point(ts["atoms"], spin,
                                labels=["Cr 3d", "21 C 2p", "22 C 2p"])
            n_c7 = nevpt2_point(c7, spin)
            n_eth = nevpt2_point(eth, 0, thr=0.4, labels=["C 2p"])
            for lv in ("e_casscf", "e_nevpt2"):
                out[lv.replace("e_", "") + "_barrier"] = round(
                    (n_ts[lv] - n_c7[lv] - n_eth[lv]) * HARTREE_KCAL, 2)
        except Exception as e:
            out["nevpt2_error"] = str(e)[:200]
    else:
        out["failed"] = "TS not interior/unphysical"
    with open(os.path.join(DIR, "lao_cr_ins_result.json"), "w") as f:
        json.dump(out, f, indent=1)
    print("[ins]", json.dumps({k: v for k, v in out.items()
                               if k != "ts_atoms"})[:400], flush=True)


def stage_desc():
    b = json.load(open(os.path.join(DIR, "lao_cr_bhe_result.json")))
    i = json.load(open(os.path.join(DIR, "lao_cr_ins_result.json")))
    res = {"descriptor": "ddE_LAO = E‡(insertion) − E‡(beta-H elim); "
                         ">0 = selective 1-hexene (LAO/PAO)",
           "bhe": {k: b.get(k) for k in ("dft_barrier", "casscf_barrier",
                                         "nevpt2_barrier", "interior", "failed")},
           "ins": {k: i.get(k) for k in ("dft_barrier", "casscf_barrier",
                                         "nevpt2_barrier", "interior", "failed")}}
    for lv in ("dft", "casscf", "nevpt2"):
        bb, ii = b.get(f"{lv}_barrier"), i.get(f"{lv}_barrier")
        if bb is not None and ii is not None:
            res[f"ddE_{lv}"] = round(ii - bb, 2)
    with open(os.path.join(DIR, "lao_cr_desc_results.json"), "w") as f:
        json.dump(res, f, indent=1)
    print(json.dumps(res, indent=1, ensure_ascii=False))
