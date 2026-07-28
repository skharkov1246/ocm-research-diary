#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
routes/msa_chain.py — высокомаржинальная мишень №1: прямой синтез
метансульфоновой кислоты (MSA) CH4 + SO3 -> CH3SO3H по радикальной цепи.
$2500-3500/т, жидкость, майнинг-иммунна (см. METHANE_HIGHVALUE_TARGETS.md).
Отличие от мёртвой газофазной уксусной: цепь MSA ЖИВЁТ (Grillo коммерциализировал),
значит считаем не «жива ли», а «где DFT промахивается в дизайне инициатора/
селективности» — прямой computable-win нашего NEVPT2-стека.

Радикальная цепь (скелет):
  (1) пропагация-добавление:  •CH3 + SO3        -> •CH3SO3   (метил к SO3)
  (2) β-распад (конкурент):    •CH3SO3           -> •CH3 + SO3 (обратно)
  (3) пропагация-HAT:          •CH3SO3 + CH4     -> CH3SO3H + •CH3 (несёт цепь)
Цепь живёт, если HAT (3) конкурентоспособен с β-распадом (2). Решают ТОЧНЫЕ
барьеры; RSO3•/сульфонил — открытая оболочка, спин-делокализация по 3 O →
DFT ненадёжен. Считаем UKS-PBE0 vs ROHF->AVAS->CASSCF->SC-NEVPT2.

Модель: газофазные молекулы, def2-SVP (S — всеэлектронный). TS(1) — пин r(C–S),
скан; TS(3) — коллинеарный [H3C···H···O–SO2CH3], пин r(C–H)+r(O–H). AVAS-метки
ИНДЕКСИРОВАНЫ по реагирующим атомам (урок antibep/acetic: жадные метки рвут CAS).

Запуск:
  python3 routes/msa_chain.py add     # стадии (1)+(2): скан C–S + концы
  python3 routes/msa_chain.py hat     # стадия (3): HAT CH4 -> метансульфонилокси
  python3 routes/msa_chain.py merge   # три числа + вердикт цепи
  python3 routes/msa_chain.py sanity  # дешёвая локальная проверка (build+SCF)
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
BASIS = os.environ.get("MSA_BASIS", "def2-svp")
MSA_SUF = "" if BASIS == "def2-svp" else "_" + BASIS.split("-")[-1]  # "" | "_tzvp"
XC = "pbe0"

# MSA_RADCAS=1 — симметричные активные пространства S-частиц: SOMO CH3SO3•
# делокализован по ВСЕМ трём O; кривобокий CAS (S 3p + один O) ломал CASSCF
# (β-scission NEVPT2 −5.9 vs DFT +19.6). Все три O в пространство; threshold
# повыше — AVAS сам подрежет слабые компоненты до NEVPT2-подъёмного размера.
RADCAS = os.environ.get("MSA_RADCAS", "0") == "1"
LAB_TS_ADD = (["0 S 3p", "1 O 2p", "2 O 2p", "3 O 2p", "4 C 2p"] if RADCAS
              else ["0 S 3p", "4 C 2p"])
LAB_SO3 = (["0 S 3p", "1 O 2p", "2 O 2p", "3 O 2p"] if RADCAS
           else ["0 S 3p", "1 O 2p"])
LAB_RAD = (["0 S 3p", "2 O 2p", "3 O 2p", "4 O 2p"] if RADCAS
           else ["0 S 3p", "2 O 2p"])
THR_S = 0.6 if RADCAS else 0.45



def M(atoms, spin, charge=0):
    return gto.M(atom=atoms, basis=BASIS, charge=charge, spin=spin,
                 verbose=0, max_memory=8000)


def uks(mol, dm0=None):
    mf = dft.UKS(mol); mf.xc = XC; mf.conv_tol = 1e-8; mf.max_cycle = 200
    mf.kernel(dm0=dm0) if dm0 is not None else mf.kernel()
    if not mf.converged:
        mf.level_shift = 0.3; mf.max_cycle = 300; mf.kernel(mf.make_rdm1())
    return mf


def rohf(mol):
    mf = scf.ROHF(mol); mf.conv_tol = 1e-8; mf.max_cycle = 200; mf.kernel()
    if not mf.converged:
        mf.level_shift = 0.4; mf.max_cycle = 300; mf.kernel(mf.make_rdm1())
    return mf


def relax(atoms, spin):
    from pyscf.geomopt.geometric_solver import optimize
    mf = uks(M(atoms, spin))
    mol = optimize(mf, maxsteps=100, assert_convergence=False)
    na = [(mol.atom_symbol(i), tuple(c))
          for i, c in enumerate(mol.atom_coords(unit="Angstrom"))]
    return na, uks(mol)


def nevpt2(atoms, spin, labels, thr=0.45):
    mol = M(atoms, spin); mf = rohf(mol)
    ncas, nelec, mo = avas.avas(mf, labels, threshold=thr)
    ss = spin / 2 * (spin / 2 + 1)
    mc = mcscf.CASSCF(mf, ncas, nelec)
    fci.addons.fix_spin_(mc.fcisolver, shift=0.2, ss=ss)
    mc.max_cycle_macro = 150; mc.conv_tol = 1e-7
    mc.kernel(mo)
    return {"e_casscf": float(mc.e_tot),
            "e_nevpt2": float(mc.e_tot + NEVPT(mc).kernel()),
            "cas": [int(np.sum(nelec)), int(ncas)], "conv": bool(mc.converged)}


# ---------------------------------------------------------------- геометрии
def ch3_atoms():
    d = 1.08
    return [("C", (0, 0, 0))] + [("H", (d*math.cos(math.radians(120*k)),
                                        d*math.sin(math.radians(120*k)), 0))
                                 for k in range(3)]


def ch4_atoms():
    d, th = 1.09, math.radians(109.47)
    a = [("C", (0, 0, 0)), ("H", (0, 0, d))]
    for k in range(3):
        phi = math.radians(120 * k)
        a.append(("H", (d*math.sin(th)*math.cos(phi), d*math.sin(th)*math.sin(phi),
                        d*math.cos(th))))
    return a


def so3_atoms():
    r = 1.43
    return [("S", (0, 0, 0))] + [("O", (r*math.cos(math.radians(120*k)),
                                        r*math.sin(math.radians(120*k)), 0))
                                 for k in range(3)]


def ch3so3_atoms(radical=True):
    """CH3–S(=O)2–O(•): S=0, CH3 сверху (+z), 3 O книзу тетраэдрично; H на one O если кислота."""
    a = [("S", (0.0, 0.0, 0.0)), ("C", (0.0, 0.0, 1.78))]
    for k in range(3):
        phi = math.radians(120 * k)
        a.append(("O", (1.35*math.cos(phi), 1.35*math.sin(phi), -0.55)))
    # CH3 водороды
    d, th = 1.09, math.radians(109.0)
    for k in range(3):
        phi = math.radians(120 * k + 60)
        a.append(("H", (d*math.sin(th)*math.cos(phi), d*math.sin(th)*math.sin(phi),
                        1.78 + d*abs(math.cos(th)))))
    return a


def ch3so3h_atoms():
    a = ch3so3_atoms()
    # H на первый O (индекс 2)
    ox = a[2][1]
    a.append(("H", (ox[0]*1.5, ox[1]*1.5, ox[2] - 0.6)))
    return a


def add_ts(rcs):
    """•CH3 + SO3, пин r(C–S): SO3 в плоскости, C подлетает по +z к S."""
    a = so3_atoms()                      # S=0, O=1,2,3
    zc = rcs
    a.append(("C", (0.0, 0.0, zc)))      # C = индекс 4
    d, th = 1.08, math.radians(100.0)
    for k in range(3):
        phi = math.radians(120 * k)
        a.append(("H", (d*math.sin(th)*math.cos(phi), d*math.sin(th)*math.sin(phi),
                        zc + d*abs(math.cos(th)))))
    return a


def hat_ts(rCH, rOH):
    """[H3C···H···O–SO2CH3]: Cα=0, H_t=1, акцепторный O=2, далее сульфонил-хвост."""
    a = [("C", (0.0, 0.0, 0.0)), ("H", (0.20, 0.0, rCH))]     # Cα, H_t
    zO = rCH + rOH
    a.append(("O", (0.40, 0.0, zO)))                          # акцепторный O = idx2
    a.append(("S", (0.40, 0.0, zO + 1.55)))                   # S
    a += [("O", (1.65, 0.0, zO + 2.05)), ("O", (-0.60, 1.25, zO + 2.05)),
          ("C", (0.40, -1.40, zO + 2.30))]                     # 2 =O + CH3-C
    d, th = 1.08, math.radians(109.0)
    for k in range(3):
        phi = math.radians(120 * k)
        a.append(("H", (0.40 + d*math.cos(phi), -1.40 + d*math.sin(phi),
                        zO + 2.30 + 0.4)))
    # метил-хвост метана (вниз)
    d2, th2 = 1.09, math.radians(105.0)
    for k in range(3):
        phi = math.radians(120 * k + 60)
        a.append(("H", (d2*math.sin(th2)*math.cos(phi), d2*math.sin(th2)*math.sin(phi),
                        -d2*abs(math.cos(th2)))))
    return a


def add_scan():
    from pyscf.geomopt.geometric_solver import optimize
    best = None; prof = []
    for rcs in (2.60, 2.35, 2.15, 1.95, 1.80):
        cf = os.path.join(DIR, f"_msa_cs_{rcs:.2f}.txt")
        open(cf, "w").write(f"$set\ndistance 1 5 {rcs:.4f}\n")   # S(1)–C(5) 1-based
        try:
            mf = uks(M(add_ts(rcs), 1))
            mol = optimize(mf, constraints=cf, maxsteps=70, assert_convergence=False)
            mfx = uks(mol); e = float(mfx.e_tot)
            na = [(mol.atom_symbol(i), tuple(c))
                  for i, c in enumerate(mol.atom_coords(unit="Angstrom"))]
            prof.append({"rcs": rcs, "e_h": e})
            if best is None or e > best[1]:
                best = (na, e, rcs)
        except Exception as ex:
            print(f"  [add] r={rcs} failed: {ex}", flush=True)
        finally:
            if os.path.exists(cf):
                os.remove(cf)
    if best is None:
        raise RuntimeError("add scan failed")
    return best, prof


def stage_add():
    t0 = time.time(); out = {"stage": "add", "model": f"gas-phase, {BASIS}"}
    ch3, m_ch3 = relax(ch3_atoms(), 1)
    so3, m_so3 = relax(so3_atoms(), 0)
    rad, m_rad = relax(ch3so3_atoms(), 1)
    (ts, e_ts, rcs), prof = add_scan()
    eR = m_ch3.e_tot + m_so3.e_tot
    out["dft"] = {"barrier_add_kcal": round((e_ts - eR) * HARTREE_KCAL, 2),
                  "dE_reac_kcal": round((m_rad.e_tot - eR) * HARTREE_KCAL, 2),
                  "barrier_bscission_kcal": round((e_ts - m_rad.e_tot) * HARTREE_KCAL, 2),
                  "ts_rcs": rcs, "profile": prof}
    n_ts = nevpt2(ts, 1, LAB_TS_ADD, thr=THR_S)
    n_ch3 = nevpt2(ch3, 1, ["0 C 2p"])
    n_so3 = nevpt2(so3, 0, LAB_SO3, thr=THR_S)
    n_rad = nevpt2(rad, 1, LAB_RAD, thr=THR_S)
    for lvl in ("casscf", "nevpt2"):
        eRc = n_ch3[f"e_{lvl}"] + n_so3[f"e_{lvl}"]
        out[lvl] = {
            "barrier_add_kcal": round((n_ts[f"e_{lvl}"] - eRc) * HARTREE_KCAL, 2),
            "barrier_bscission_kcal": round((n_ts[f"e_{lvl}"] - n_rad[f"e_{lvl}"])
                                            * HARTREE_KCAL, 2)}
    out["ts_cas"] = n_ts["cas"]; out["wall_s"] = round(time.time() - t0, 1)
    json.dump(out, open(os.path.join(DIR, "msa_add.json"), "w"), indent=1)
    print(f"[add] barrier add/bscission kcal: DFT {out['dft']['barrier_add_kcal']}/"
          f"{out['dft']['barrier_bscission_kcal']} | NEVPT2 "
          f"{out['nevpt2']['barrier_add_kcal']}/{out['nevpt2']['barrier_bscission_kcal']}"
          f" ({out['wall_s']}s)", flush=True)


def stage_hat():
    t0 = time.time()
    out = {"stage": "hat", "model": f"pinned-relax collinear HAT, {BASIS}"}
    ch4, m_ch4 = relax(ch4_atoms(), 0)
    rad, m_rad = relax(ch3so3_atoms(), 1)
    from pyscf.geomopt.geometric_solver import optimize
    # РЕЛАКС TS: пин r(C–H)=dist 1-2, r(H–O)=dist 2-3, хвост релаксируется.
    # Урок radcas-прогона (806 ккал): optimize с клэш-старта может вернуть
    # геометрию хуже стартовой, а свежий SCF на ней — чужую ветку. Поэтому:
    # warm-start плотностью по цепочке, проверка сходимости КАЖДОЙ точки,
    # весь скан в JSON, TS = максимум только по сошедшимся точкам,
    # sanity-гейт на нефизичный барьер.
    # v4: v3 дал плоский низкий профиль с максимумом на РАННЕМ краю ⇒ перевал
    # раньше rCH=1.25, а идеализированный коллинеарный старт там застревает
    # (v2: +82 ккал). Решение: ранние точки стартуют С РЕЛАКСИРОВАННОЙ
    # ГЕОМЕТРИИ якоря (1.25,1.15) — geomeTRIC сам дожимает пины; поздняя ветка
    # как в v3. Профиль собирается в порядке пути (по rCH).
    def _hat_point(start_atoms, rCH, rOH, dm):
        cf = os.path.join(DIR, f"_msa_hat_{rCH:.2f}.txt")
        open(cf, "w").write(
            f"$set\ndistance 1 2 {rCH:.4f}\ndistance 2 3 {rOH:.4f}\n")
        try:
            mf0 = uks(M(start_atoms, 1), dm0=dm)
            mol = optimize(mf0, constraints=cf, maxsteps=180,
                           assert_convergence=False)
            mfx = uks(mol, dm0=mf0.make_rdm1())
            na = [(mol.atom_symbol(i), tuple(c))
                  for i, c in enumerate(mol.atom_coords(unit="Angstrom"))]
            return (float(mfx.e_tot), bool(mfx.converged), na,
                    mfx.make_rdm1())
        except Exception as ex:
            print(f"  [hat] rCH={rCH} failed: {ex}", flush=True)
            return None, False, None, dm
        finally:
            if os.path.exists(cf):
                os.remove(cf)

    # v7: поточечное сохранение скана (v6 убит timeout за полчаса до конца —
    # 8 ч скана пропали, т.к. JSON писался только в финале). Resume по (rCH,rOH).
    scan_path = os.path.join(DIR, f"msa_hat_scan{MSA_SUF}.json")
    pts = []
    if os.path.exists(scan_path):
        try:
            pts = json.load(open(scan_path))["pts"]
            print(f"[hat] resume: {len(pts)} точек скана с диска", flush=True)
        except Exception:
            pts = []
    done_keys = {(p["rCH"], p["rOH"]) for p in pts}

    def _save_scan():
        with open(scan_path, "w") as f:
            json.dump({"pts": pts}, f, indent=1)

    # якорь + поздняя ветка (v3-проверенная): идеализированный старт, цепь dm
    anchor_atoms, anchor_dm = None, None
    dm = None
    for rCH, rOH in ((1.25, 1.15), (1.35, 1.05), (1.50, 0.98),
                     (1.65, 0.965), (1.80, 0.955)):
        if (rCH, rOH) in done_keys:
            saved = next(p for p in pts if (p["rCH"], p["rOH"]) == (rCH, rOH))
            if anchor_atoms is None:
                anchor_atoms = [(s, tuple(c)) for s, c in saved["atoms"]]
            continue
        e, conv, na, dm = _hat_point(hat_ts(rCH, rOH), rCH, rOH, dm)
        if e is not None and np.isfinite(e):
            pts.append({"rCH": rCH, "rOH": rOH, "e_h": e,
                        "scf_converged": conv, "atoms": na})
            _save_scan()
            print(f"  [hat] point rCH={rCH} rOH={rOH} E={e:.6f} conv={conv}",
                  flush=True)
            if anchor_atoms is None:
                anchor_atoms, anchor_dm = na, dm
    # ранняя ветка: старт с релакс-геометрии якоря, пины стягиваются к реагентам.
    # v5: +подходные точки (rCH≈равновесие C–H, O–H далеко) — v4 показал максимум
    # 18.8 на 1.12 (совпал со старым валидированным ~19.4), но без R-стороны
    # он краевой; подход делает его интерьерным.
    dm = anchor_dm
    start = anchor_atoms
    # TZVP-прогон 27.07: максимум упёрся в край (1.09,1.62) — TS на большом
    # базисе сдвинулся к реагентам (ранний TS). Продолжаем разведение, чтобы
    # увидеть спад к реагентному пределу и получить интерьерный перевал.
    for rCH, rOH in ((1.18, 1.24), (1.12, 1.36), (1.10, 1.48), (1.09, 1.62),
                     (1.08, 1.85), (1.07, 2.10), (1.06, 2.40)):
        if (rCH, rOH) in done_keys:
            saved = next(p for p in pts if (p["rCH"], p["rOH"]) == (rCH, rOH))
            start = [(s, tuple(c)) for s, c in saved["atoms"]]
            continue
        if start is None:
            break
        e, conv, na, dm = _hat_point(start, rCH, rOH, dm)
        if e is not None and np.isfinite(e):
            pts.append({"rCH": rCH, "rOH": rOH, "e_h": e,
                        "scf_converged": conv, "atoms": na})
            _save_scan()
            print(f"  [hat] early rCH={rCH} rOH={rOH} E={e:.6f} conv={conv}",
                  flush=True)
            start = na
    # v9: перепроверка SCF-ветки — интерьерная точка на >10 ккал выше ОБОИХ
    # соседей пере-релаксируется со старта из геометрии нижнего соседа
    # (warm-chain, как в OCM-протоколе). Это не маска: точка остаётся, берётся
    # МИНИМУМ, обе энергии документируются в чекпойнте. Ловит застрявший релакс
    # /чужую SCF-ветку идеализированного старта (кандидат: (1.50,0.98) v7).
    pts.sort(key=lambda p: p["rCH"])
    for i in range(1, len(pts) - 1):
        p = pts[i]
        if p.get("recheck_done"):
            continue
        gapL = (p["e_h"] - pts[i - 1]["e_h"]) * HARTREE_KCAL
        gapR = (p["e_h"] - pts[i + 1]["e_h"]) * HARTREE_KCAL
        if gapL > 10.0 and gapR > 10.0:
            nb = pts[i - 1] if pts[i - 1]["e_h"] < pts[i + 1]["e_h"] else pts[i + 1]
            st = [(s, tuple(c)) for s, c in nb["atoms"]]
            e, conv, na, _ = _hat_point(st, p["rCH"], p["rOH"], None)
            p["recheck_done"] = True
            p["e_h_first"] = p["e_h"]
            if (e is not None and np.isfinite(e) and conv
                    and e < p["e_h"] - 1e-4):
                print(f"  [hat] recheck rCH={p['rCH']}: нижняя ветка на "
                      f"{(p['e_h'] - e) * HARTREE_KCAL:.1f} ккал ниже — заменяю",
                      flush=True)
                p["e_h"], p["atoms"], p["scf_converged"] = e, na, conv
            else:
                print(f"  [hat] recheck rCH={p['rCH']}: исходная точка "
                      f"подтверждена", flush=True)
            _save_scan()
    eR = m_rad.e_tot + m_ch4.e_tot
    out["scan"] = [{k: p[k] for k in ("rCH", "rOH", "e_h", "scf_converged")}
                   for p in pts]
    out["scan_rel_kcal"] = [round((p["e_h"] - eR) * HARTREE_KCAL, 2) for p in pts]
    # порядок счёта (якорь-первым) != ход реакции: сортируем по rCH, иначе
    # соседство/интерьерность считаются по бессмысленным парам (баг вскрыт 27.07)
    ok = sorted([p for p in pts if p["scf_converged"]], key=lambda p: p["rCH"])
    if len(ok) < 3:
        raise RuntimeError(f"hat scan: only {len(ok)} converged points")
    # отбраковка выбросов: точка на >25 ккал выше ОБОИХ соседей (или
    # единственного соседа на краю) — застрявший релакс, не поверхность
    rel = [(p["e_h"] - eR) * HARTREE_KCAL for p in ok]
    keep = []
    for i, p in enumerate(ok):
        nb = [rel[j] for j in (i - 1, i + 1) if 0 <= j < len(ok)]
        if all(rel[i] - x > 25.0 for x in nb):
            print(f"  [hat] outlier dropped: rCH={p['rCH']} rel={rel[i]:.1f}",
                  flush=True)
            continue
        keep.append((p, rel[i]))
    if len(keep) < 3:
        raise RuntimeError("hat scan: too few points after outlier rejection")
    imax = max(range(len(keep)), key=lambda i: keep[i][1])
    best, dft_bar = keep[imax]
    interior = 0 < imax < len(keep) - 1
    ts, rCH, rOH = best["atoms"], best["rCH"], best["rOH"]
    out["dft"] = {"barrier_hat_kcal": round(dft_bar, 2),
                  "ts_rCH": rCH, "ts_rOH": rOH, "ts_interior": interior}
    if not interior or not (0.0 < dft_bar < 150.0):
        out["hat_stage_failed"] = (
            f"TS not bracketed/unphysical: barrier {dft_bar:.1f} kcal, "
            f"interior={interior} — extend scan, no verdict")
        out["wall_s"] = round(time.time() - t0, 1)
        json.dump(out, open(os.path.join(DIR, f"msa_hat{MSA_SUF}.json"), "w"), indent=1)
        print(f"[hat] FAILED sanity gate: {out['hat_stage_failed']}", flush=True)
        return
    ts_labels = (["0 C 2p", "1 H 1s", "2 O 2p", "3 S 3p", "4 O 2p", "5 O 2p"]
                 if RADCAS else ["0 C 2p", "1 H 1s", "2 O 2p"])
    out["ts_atoms"] = [[s, [round(x, 6) for x in c]] for s, c in ts]
    n_ts = nevpt2(ts, 1, ts_labels, thr=0.45)
    n_ch4 = nevpt2(ch4, 0, ["0 C 2p", "1 H 1s"], thr=0.6)
    n_rad = nevpt2(rad, 1, LAB_RAD, thr=THR_S)
    for lvl in ("casscf", "nevpt2"):
        eRc = n_rad[f"e_{lvl}"] + n_ch4[f"e_{lvl}"]
        out[lvl] = {"barrier_hat_kcal":
                    round((n_ts[f"e_{lvl}"] - eRc) * HARTREE_KCAL, 2)}
    out["ts_cas"] = n_ts["cas"]; out["wall_s"] = round(time.time() - t0, 1)
    json.dump(out, open(os.path.join(DIR, f"msa_hat{MSA_SUF}.json"), "w"), indent=1)
    print(f"[hat] barrier kcal: DFT {out['dft']['barrier_hat_kcal']} | NEVPT2 "
          f"{out['nevpt2']['barrier_hat_kcal']} ({out['wall_s']}s)", flush=True)


def stage_sanity():
    for name, atoms, spin in (("ch3so3•", ch3so3_atoms(), 1), ("SO3", so3_atoms(), 0),
                              ("MSA", ch3so3h_atoms(), 0)):
        try:
            mf = scf.ROHF(M(atoms, spin)); mf.max_cycle = 100; mf.kernel()
            print(f"[sanity] {name}: conv={mf.converged} E={mf.e_tot:.4f} nat={len(atoms)}",
                  flush=True)
        except Exception as e:
            print(f"[sanity] {name} FAIL: {e}", flush=True)


def stage_merge():
    add = json.load(open(os.path.join(DIR, "msa_add.json")))
    hat_path = os.path.join(DIR, "msa_hat.json")
    if os.path.exists(hat_path):
        hat = json.load(open(hat_path))
    else:  # v8: hat убит таймаутом до записи файла — merge не падает
        hat = {"hat_stage_failed": "msa_hat.json missing (stage killed?)"}
    if hat.get("hat_stage_failed") or "nevpt2" not in hat:
        res = {"route": "direct CH4 + SO3 -> CH3SO3H radical chain",
               "numbers_nevpt2": {
                   "add_barrier": add["nevpt2"]["barrier_add_kcal"],
                   "bscission_barrier": add["nevpt2"]["barrier_bscission_kcal"]},
               "chain_verdict": "NO VERDICT: HAT stage failed sanity gate — "
                                "needs rerun, add/bscission numbers stand",
               "hat_failure": hat.get("hat_stage_failed", "no nevpt2 block")}
        json.dump(res, open(os.path.join(DIR, "msa_results.json"), "w"), indent=1)
        print(json.dumps(res, ensure_ascii=False, indent=1))
        return
    nv = hat["nevpt2"]["barrier_hat_kcal"]
    if not (0.0 < nv < 100.0):
        res = {"route": "direct CH4 + SO3 -> CH3SO3H radical chain",
               "numbers_nevpt2": {
                   "add_barrier": add["nevpt2"]["barrier_add_kcal"],
                   "bscission_barrier": add["nevpt2"]["barrier_bscission_kcal"],
                   "hat_flagged": nv},
               "numbers_dft": {
                   "add_barrier": add["dft"]["barrier_add_kcal"],
                   "bscission_barrier": add["dft"]["barrier_bscission_kcal"],
                   "hat_propagation_barrier": hat["dft"]["barrier_hat_kcal"]},
               "chain_verdict": "DFT-ALIVE; NEVPT2-HAT flagged unphysical "
                                "(CAS imbalance) — no correlated verdict",
               }
        json.dump(res, open(os.path.join(DIR, "msa_results.json"), "w"),
                  indent=1)
        print(json.dumps(res, ensure_ascii=False, indent=1))
        return
    alive = nv <= add["nevpt2"]["barrier_bscission_kcal"] + 5
    res = {"route": "direct CH4 + SO3 -> CH3SO3H (methanesulfonic acid) radical chain",
           "numbers_nevpt2": {
               "add_barrier": add["nevpt2"]["barrier_add_kcal"],
               "bscission_barrier": add["nevpt2"]["barrier_bscission_kcal"],
               "hat_propagation_barrier": hat["nevpt2"]["barrier_hat_kcal"]},
           "numbers_dft": {
               "add_barrier": add["dft"]["barrier_add_kcal"],
               "bscission_barrier": add["dft"]["barrier_bscission_kcal"],
               "hat_propagation_barrier": hat["dft"]["barrier_hat_kcal"]},
           "chain_verdict": ("ALIVE: HAT propagation competitive with beta-scission — "
                             "radical MSA chain viable (consistent with commercial route)"
                             if alive else
                             "gas-phase skeleton marginal — real route uses initiator/"
                             "conditions to bias chain; DFT-NEVPT2 gap guides initiator"),
           "note": "gas-phase skeleton. DFT-vs-NEVPT2 gap on RSO3* open-shell states is "
                   "the de-risk / selectivity+initiator design lever. Screening grade. "
                   "$2500-3500/t liquid, mining-immune (value decoupled from energy)."}
    json.dump(res, open(os.path.join(DIR, "msa_results.json"), "w"), indent=1)
    print(json.dumps(res, ensure_ascii=False, indent=1))


if __name__ == "__main__":
    {"add": stage_add, "hat": stage_hat, "sanity": stage_sanity,
     "merge": stage_merge}[sys.argv[1]]()
