#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
routes/halogen_descriptor.py — дескриптор галоген-маршрута (маршрут 4 разведки):
барьер H-абстракции галоген-радикалом X• + CH4 → CH3• + HX для X = Cl, Br, I.

Почему это «наша ниша» (из METHANE_ROUTES_RESEARCH.md). TS реакции — коллинеарная
3-центровая-3-электронная система [CH3···H···X]‡ (открытая оболочка, дублет): рвётся
σ(C–H), образуется σ(H–X), плюс радикальная орбиталь галогена. Обычный DFT на таких
open-shell TS ненадёжен (лит.: ошибка C–X/C–H термохимии растёт с числом галогенов до
77 кДж/моль; барьеры считали QCI/коррелированными методами, не DFT). Лестница барьеров
Cl < Br < I управляет и селективностью, и тем, какой галоген тянет цикл. Референс из
литературы (для сверки порядка): Br+CH4 ≈ 18 ккал/моль, I+CH4 ≈ 34 ккал/моль.

Модель: газофазные молекулы, def2-SVP (+ECP для Br/I), дублет по всему пути.
Уровни: UKS-PBE0 (DFT) и ROHF→AVAS(C 2p + X np + H 1s)→CASSCF→SC-NEVPT2.
TS: 2-координатный пин r(C–H)+r(H–X), скан вдоль переносимого H, максимум = TS
(тот же приём, что HAT на Mn=O в Этапах 15/16, но здесь центр — сам галоген-радикал).

Барьер = E(TS) − [E(X•) + E(CH4)];  ΔE_reac = [E(CH3•)+E(HX)] − [E(X•)+E(CH4)].
Fan-out: один галоген на запуск → параллельно.

Запуск:
  python3 routes/halogen_descriptor.py run   Cl|Br|I
  python3 routes/halogen_descriptor.py merge
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
BASIS = "def2-svp"
XC = "pbe0"

# np-метка валентной p-оболочки галогена для AVAS
XCONF = {
    "F":  {"np": "F 2p",  "rHX": 0.92, "rHX_ts": 1.15},
    "Cl": {"np": "Cl 3p", "rHX": 1.27, "rHX_ts": 1.45},
    "Br": {"np": "Br 4p", "rHX": 1.41, "rHX_ts": 1.55},
    "I":  {"np": "I 5p",  "rHX": 1.61, "rHX_ts": 1.75},
}
RCH_TS = 1.30   # растянутая C–H в TS-области (скан вокруг)


def M(atoms, spin, charge=0):
    return gto.M(atom=atoms, basis=BASIS, ecp=BASIS, charge=charge, spin=spin,
                 verbose=0, max_memory=6000)


def uks(mol, dm0=None):
    mf = dft.UKS(mol).density_fit(); mf.xc = XC; mf.conv_tol = 1e-8
    mfn = scf.newton(mf)
    mfn.kernel(dm0=dm0) if dm0 is not None else mfn.kernel()
    if not mfn.converged:
        m2 = dft.UKS(mol).density_fit(); m2.xc = XC; m2.level_shift = 0.2
        m2.max_cycle = 300; m2.kernel(dm0=dm0)
        mfn = scf.newton(m2); mfn.kernel(m2.make_rdm1())
    return mfn


def rohf(mol):
    mf = scf.ROHF(mol); mf.conv_tol = 1e-8
    mfn = scf.newton(mf); mfn.kernel()
    if not mfn.converged:
        m2 = scf.ROHF(mol); m2.level_shift = 0.4; m2.max_cycle = 300; m2.kernel()
        mfn = scf.newton(m2); mfn.kernel(m2.make_rdm1())
    return mfn


def relax_uks(atoms, spin):
    from pyscf.geomopt.geometric_solver import optimize
    mf = uks(M(atoms, spin))
    mol = optimize(mf, maxsteps=100, assert_convergence=False)
    na = [(mol.atom_symbol(i), tuple(c))
          for i, c in enumerate(mol.atom_coords(unit="Angstrom"))]
    return na, uks(mol)


def nevpt2(atoms, spin, avas_labels, thr=0.5):
    """ROHF→AVAS→CASSCF→SC-NEVPT2 для дублет-радикалов и молекул."""
    mol = M(atoms, spin)
    mf = rohf(mol)
    ncas, nelec, mo = avas.avas(mf, avas_labels, threshold=thr)
    ss = spin / 2 * (spin / 2 + 1)
    mc = mcscf.CASSCF(mf, ncas, nelec)
    fci.addons.fix_spin_(mc.fcisolver, shift=0.2, ss=ss)
    mc.max_cycle_macro = 120; mc.conv_tol = 1e-7
    mc.kernel(mo)
    e_pt = float(mc.e_tot + NEVPT(mc).kernel())
    return {"e_casscf": float(mc.e_tot), "e_nevpt2": e_pt,
            "cas": [int(np.sum(nelec)), int(ncas)], "conv": bool(mc.converged)}


# ---------------------------------------------------------------- геометрии
def ch4_atoms():
    d, th = 1.09, math.radians(109.47)
    a = [("C", (0, 0, 0)), ("H", (0, 0, d))]
    for k in range(3):
        phi = math.radians(120 * k)
        a.append(("H", (d*math.sin(th)*math.cos(phi), d*math.sin(th)*math.sin(phi),
                        d*math.cos(th))))
    return a


def ch3_atoms():
    d = 1.08
    return [("C", (0, 0, 0))] + [("H", (d*math.cos(math.radians(120*k)),
                                        d*math.sin(math.radians(120*k)), 0))
                                 for k in range(3)]


def hx_atoms(X, r):
    return [("H", (0, 0, 0)), (X, (0, 0, r))]


def ts_atoms(X, rCH, rHX):
    """Коллинеарный [CH3···H···X]‡ по оси z: C(0) — H_t — X, метил зонтиком от C."""
    zC = 0.0; zH = rCH; zX = rCH + rHX
    d, th = 1.08, math.radians(105.0)
    a = [("C", (0, 0, zC)), ("H", (0, 0, zH)), (X, (0, 0, zX))]  # C, H_transfer, X
    for k in range(3):
        phi = math.radians(120 * k)
        a.append(("H", (d*math.sin(th)*math.cos(phi), d*math.sin(th)*math.sin(phi),
                        zC - d*math.cos(th))))
    return a


def ts_scan(X, rHX_ts):
    """Скан переносимого H: пин r(C–H)+r(H–X), релакс метила, максимум = TS."""
    from pyscf.geomopt.geometric_solver import optimize
    grid = [(1.20, rHX_ts+0.15), (1.30, rHX_ts), (1.40, rHX_ts-0.10),
            (1.50, rHX_ts-0.20)]
    best = None
    for rCH, rHX in grid:
        cf = os.path.join(DIR, f"_hx_{X}.txt")
        open(cf, "w").write(f"$set\ndistance 1 2 {rCH:.4f}\ndistance 2 3 {rHX:.4f}\n")
        mf = uks(M(ts_atoms(X, rCH, rHX), 1))
        mol = optimize(mf, constraints=cf, maxsteps=80, assert_convergence=False)
        os.remove(cf)
        mfx = uks(mol)
        e = float(mfx.e_tot)
        na = [(mol.atom_symbol(i), tuple(c))
              for i, c in enumerate(mol.atom_coords(unit="Angstrom"))]
        if best is None or e > best[1]:
            best = (na, e, rCH, rHX)
    return best   # (atoms, e_ts_dft, rCH, rHX)


def stage_run(X):
    cfg = XCONF[X]
    out = {"halogen": X, "model": "gas-phase X + CH4 -> CH3 + HX, def2-SVP+ECP"}
    path = os.path.join(DIR, f"halogen_{X}.json")
    t0 = time.time()
    # реагенты/продукты (релакс DFT, затем NEVPT2 на релакс-геометрии)
    ch4, mf_ch4 = relax_uks(ch4_atoms(), 0)
    ch3, mf_ch3 = relax_uks(ch3_atoms(), 1)
    hx, mf_hx = relax_uks(hx_atoms(X, cfg["rHX"]), 0)
    xr = M([(X, (0, 0, 0))], 1)                # X• радикал
    mf_x = uks(xr)
    # TS-скан (DFT)
    ts, e_ts_dft, rCH, rHX = ts_scan(X, cfg["rHX_ts"])
    # DFT-энергетика
    eR = mf_x.e_tot + mf_ch4.e_tot
    eP = mf_ch3.e_tot + mf_hx.e_tot
    out["dft"] = {
        "barrier_kcal": round((e_ts_dft - eR) * HARTREE_KCAL, 2),
        "dE_reac_kcal": round((eP - eR) * HARTREE_KCAL, 2),
        "ts_rCH": rCH, "ts_rHX": round(rHX, 3)}
    # NEVPT2-энергетика (те же геометрии)
    lab_ts = ["C 2p", cfg["np"], "H 1s"]        # 3c-3e {C,H,X}
    n_ts = nevpt2(ts, 1, lab_ts, thr=0.4)
    n_x = nevpt2([(X, (0, 0, 0))], 1, [cfg["np"]], thr=0.5)
    n_ch4 = nevpt2(ch4, 0, ["C 2p", "H 1s"], thr=0.6)
    n_ch3 = nevpt2(ch3, 1, ["C 2p"], thr=0.5)
    n_hx = nevpt2(hx, 0, [cfg["np"], "H 1s"], thr=0.6)
    for lvl in ("casscf", "nevpt2"):
        eR = n_x[f"e_{lvl}"] + n_ch4[f"e_{lvl}"]
        eP = n_ch3[f"e_{lvl}"] + n_hx[f"e_{lvl}"]
        out[lvl] = {"barrier_kcal": round((n_ts[f"e_{lvl}"] - eR) * HARTREE_KCAL, 2),
                    "dE_reac_kcal": round((eP - eR) * HARTREE_KCAL, 2)}
    out["ts_cas"] = n_ts["cas"]
    out["wall_s"] = round(time.time() - t0, 1)
    json.dump(out, open(path, "w"), indent=1)
    print(f"[{X}] barrier kcal/mol: DFT {out['dft']['barrier_kcal']} | "
          f"CASSCF {out['casscf']['barrier_kcal']} | NEVPT2 {out['nevpt2']['barrier_kcal']} "
          f"| ΔE_reac(NEVPT2) {out['nevpt2']['dE_reac_kcal']} | TS-CAS {n_ts['cas']} "
          f"({out['wall_s']}s)", flush=True)


def stage_merge():
    rows = []
    for X in ("F", "Cl", "Br", "I"):
        p = os.path.join(DIR, f"halogen_{X}.json")
        if os.path.exists(p):
            d = json.load(open(p))
            rows.append({"X": X,
                         "barrier_dft": d["dft"]["barrier_kcal"],
                         "barrier_casscf": d["casscf"]["barrier_kcal"],
                         "barrier_nevpt2": d["nevpt2"]["barrier_kcal"],
                         "dE_reac_nevpt2": d["nevpt2"]["dE_reac_kcal"]})
    res = {"descriptor": "H-abstraction barrier X. + CH4 -> CH3 + HX (Cl<Br<I ladder)",
           "note": "gas-phase, def2-SVP+ECP, open-shell doublet TS (3c-3e). NEVPT2 is "
                   "the multireference-corrected barrier; DFT-NEVPT2 gap = the de-risk. "
                   "Lit. reference order: Br~18, I~34 kcal/mol.",
           "halogens": rows}
    json.dump(res, open(os.path.join(DIR, "halogen_descriptor_results.json"), "w"),
              indent=1)
    print(json.dumps(res, ensure_ascii=False, indent=1))


if __name__ == "__main__":
    if sys.argv[1] == "run":
        stage_run(sys.argv[2])
    elif sys.argv[1] == "merge":
        stage_merge()
    else:
        raise SystemExit("usage: run Cl|Br|I|F, or merge")
