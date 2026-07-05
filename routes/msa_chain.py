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
XC = "pbe0"


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
    n_ts = nevpt2(ts, 1, ["0 S 3p", "4 C 2p"])
    n_ch3 = nevpt2(ch3, 1, ["0 C 2p"])
    n_so3 = nevpt2(so3, 0, ["0 S 3p", "1 O 2p"])
    n_rad = nevpt2(rad, 1, ["0 S 3p", "2 O 2p"])
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
    t0 = time.time(); out = {"stage": "hat", "model": f"rigid collinear HAT, {BASIS}"}
    ch4, m_ch4 = relax(ch4_atoms(), 0)
    rad, m_rad = relax(ch3so3_atoms(), 1)
    from pyscf.geomopt.geometric_solver import optimize
    best = None
    for rCH, rOH in ((1.15, 1.30), (1.25, 1.15), (1.35, 1.05), (1.50, 0.98)):
        # РЕЛАКС TS (не rigid!): пин r(C–H)=dist 1-2, r(H–O)=dist 2-3, остальное
        # (сульфонил-хвост) релаксируется → убирает раздутый rigid-барьер. Изогнутый
        # старт (H_t сдвинут) обходит линейную сингулярность. Fallback — одноточка.
        cf = os.path.join(DIR, f"_msa_hat_{rCH:.2f}.txt")
        open(cf, "w").write(f"$set\ndistance 1 2 {rCH:.4f}\ndistance 2 3 {rOH:.4f}\n")
        try:
            mf0 = uks(M(hat_ts(rCH, rOH), 1))
            mol = optimize(mf0, constraints=cf, maxsteps=60, assert_convergence=False)
            mfx = uks(mol); e = float(mfx.e_tot)
            na = [(mol.atom_symbol(i), tuple(c))
                  for i, c in enumerate(mol.atom_coords(unit="Angstrom"))]
        except Exception as ex:
            print(f"  [hat] rCH={rCH} relax fell back to SP: {ex}", flush=True)
            try:
                mfx = uks(M(hat_ts(rCH, rOH), 1)); e = float(mfx.e_tot)
                na = hat_ts(rCH, rOH)
            except Exception:
                continue
        finally:
            if os.path.exists(cf):
                os.remove(cf)
        if not np.isfinite(e):
            continue
        if best is None or e > best[1]:
            best = (na, e, rCH, rOH)
    if best is None:
        raise RuntimeError("hat scan failed")
    ts, e_ts, rCH, rOH = best
    eR = m_rad.e_tot + m_ch4.e_tot
    out["dft"] = {"barrier_hat_kcal": round((e_ts - eR) * HARTREE_KCAL, 2),
                  "ts_rCH": rCH, "ts_rOH": rOH}
    n_ts = nevpt2(ts, 1, ["0 C 2p", "1 H 1s", "2 O 2p"], thr=0.4)
    n_ch4 = nevpt2(ch4, 0, ["0 C 2p", "1 H 1s"], thr=0.6)
    n_rad = nevpt2(rad, 1, ["0 S 3p", "2 O 2p"])
    for lvl in ("casscf", "nevpt2"):
        eRc = n_rad[f"e_{lvl}"] + n_ch4[f"e_{lvl}"]
        out[lvl] = {"barrier_hat_kcal":
                    round((n_ts[f"e_{lvl}"] - eRc) * HARTREE_KCAL, 2)}
    out["ts_cas"] = n_ts["cas"]; out["wall_s"] = round(time.time() - t0, 1)
    json.dump(out, open(os.path.join(DIR, "msa_hat.json"), "w"), indent=1)
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
    hat = json.load(open(os.path.join(DIR, "msa_hat.json")))
    alive = hat["nevpt2"]["barrier_hat_kcal"] <= add["nevpt2"]["barrier_bscission_kcal"] + 5
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
