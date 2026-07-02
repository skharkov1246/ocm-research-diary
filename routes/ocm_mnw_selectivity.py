#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
routes/ocm_mnw_selectivity.py — Этап 15 OCM: закрыть двухсубстратный ΔΔE‡
на реальном радикально-кислородном центре Mn=O (Mn–Na₂WO₄/SiO₂, минимальная модель).

Чинит ровно то, на чём упал Этап 14:
  • ЕДИНОЕ курированное активное пространство CAS(9e,10o) для всех 4 структур
    (R/TS × CH₄/C₂H₆): AVAS(Mn 3d + O 2p + 1s переносимого H). Проверено:
    размер стабилен к порогу 0.05–0.2 (метки с C 2p «прыгают» 13o↔15o — не берём).
  • Настоящий HAT-профиль: 2-координатный пин r(C–H)+r(O–H) (geomeTRIC),
    каркас релаксируется; TS = внутренний максимум профиля (краевой — флагуется).
  • Оба субстрата до конца: CH₄ И C₂H₆, один протокол, одна сетка.

Протокол: UKS-PBE0/def2-SVP геометрия и DFT-барьеры; ROHF→AVAS→CASSCF(9e,10o)
(робастный: newton-SCF старт, натурбитали, fix_spin)→SC-NEVPT2 — энергетика.
Секстет (S=5/2) по всему пути (спин-сохраняющий HAT). Модель — диатомный MnO
(секстет ⁶Σ⁺ — известное основное состояние) как минимальный стенд-ин
рабочего Mn=O; НЕ периодический Mn–Na₂WO₄/SiO₂ (оговорка в дневнике).

Дескриптор (конвенция «карты потолка»): ΔΔE‡ = барьер(C₂H₆) − барьер(CH₄);
< 0 — центр доокисляет C₂-продукт быстрее метана (стена BEP), для 70% нужно ≳ +4.

Запуск (стадии независимы, результаты инкрементально в JSON):
  python3 -u routes/ocm_mnw_selectivity.py geom   ch4|c2h6   # скан-профиль
  python3 -u routes/ocm_mnw_selectivity.py energy ch4|c2h6   # CASSCF+NEVPT2 в R и TS
  python3 -u routes/ocm_mnw_selectivity.py merge             # итоговый results-JSON
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
SPIN = 5                    # секстет: MnO(⁶Σ⁺) + закрытооболочечный субстрат
XC = "pbe0"
CAS_TARGET = (9, 10)        # (ne, no) — единое для всех 4 структур
AVAS_LABELS = ["Mn 3d", "O 2p", "2 H 1s"]   # атом 2 (0-based) — переносимый H

# HAT-сетка: пары (r_CH, r_OH), Å. R-точка: пин только r(O–H)=2.2 (пре-комплекс).
R_POINT = (None, 2.20)
PATH = [(1.15, 1.70), (1.20, 1.50), (1.25, 1.35), (1.30, 1.25),
        (1.35, 1.17), (1.40, 1.10), (1.50, 1.02)]


# ---------------------------------------------------------------- геометрии
def start_atoms(sub):
    """Стартовая геометрия Mn=O ... H–R (коллинеарная HAT-ось по z).
    Порядок атомов фиксирован: Mn(0), O(1), H_перенос(2), C(3), ..."""
    zMn, zO = -1.65, 0.0
    rOH0, rCH0 = 2.20, 1.10
    zH = rOH0
    zC = rOH0 + rCH0
    d, th = 1.09, math.radians(108.0)
    atoms = [("Mn", (0.0, 0.0, zMn)), ("O", (0.0, 0.0, zO)),
             ("H", (0.0, 0.0, zH)), ("C", (0.0, 0.0, zC))]
    if sub == "ch4":
        for k in range(3):
            phi = math.radians(120.0 * k)
            atoms.append(("H", (d * math.sin(th) * math.cos(phi),
                                d * math.sin(th) * math.sin(phi),
                                zC - d * math.cos(th))))
    elif sub == "c2h6":
        rCC = 1.53
        zC2 = zC + rCC * math.cos(math.radians(35.0))
        xC2 = rCC * math.sin(math.radians(35.0))
        atoms.append(("C", (xC2, 0.0, zC2)))
        for k in (1, 2):   # два оставшихся H на C1
            phi = math.radians(120.0 * k + 60.0)
            atoms.append(("H", (d * math.sin(th) * math.cos(phi),
                                d * math.sin(th) * math.sin(phi),
                                zC - d * math.cos(th))))
        for k in range(3):  # H на C2
            phi = math.radians(120.0 * k)
            atoms.append(("H", (xC2 + d * math.sin(th) * math.cos(phi),
                                d * math.sin(th) * math.sin(phi),
                                zC2 + d * math.cos(th))))
    else:
        raise ValueError(sub)
    return atoms


def build_mol(atoms):
    return gto.M(atom=atoms, basis=BASIS, charge=0, spin=SPIN,
                 verbose=0, max_memory=6000)


def uks(mol, dm0=None):
    """UKS-PBE0 (RI/DF — кратно быстрее на градиентах), newton-primary
    (устойчив на Mn), DIIS-фоллбэк."""
    mf = dft.UKS(mol).density_fit()
    mf.xc = XC
    mf.conv_tol = 1e-8
    mfn = scf.newton(mf)
    mfn.kernel(dm0=dm0) if dm0 is not None else mfn.kernel()
    if not mfn.converged:
        mf2 = dft.UKS(mol).density_fit()
        mf2.xc = XC
        mf2.level_shift = 0.3
        mf2.max_cycle = 300
        mf2.kernel(dm0=dm0)
        mfn = scf.newton(mf2)
        mfn.kernel(mf2.make_rdm1())
    return mfn


def constrained_opt(atoms, rCH, rOH, tag):
    """Релакс-оптимизация с пином r(C–H) и/или r(O–H) (geomeTRIC, $set)."""
    from pyscf.geomopt.geometric_solver import optimize
    cfile = os.path.join(DIR, f"_constr_{tag}.txt")
    lines = ["$set"]
    if rOH is not None:
        lines.append(f"distance 2 3 {rOH:.4f}")   # O(2)–H(3), 1-based
    if rCH is not None:
        lines.append(f"distance 3 4 {rCH:.4f}")   # H(3)–C(4)
    with open(cfile, "w") as f:
        f.write("\n".join(lines) + "\n")
    mf = uks(build_mol(atoms))
    mol_eq = optimize(mf, constraints=cfile, maxsteps=120,
                      convergence_energy=1e-6, convergence_grms=3e-4,
                      convergence_gmax=6e-4, convergence_drms=1.5e-3,
                      convergence_dmax=2.5e-3)
    os.remove(cfile)
    new_atoms = [(mol_eq.atom_symbol(i), tuple(c))
                 for i, c in enumerate(mol_eq.atom_coords(unit="Angstrom"))]
    mf_final = uks(mol_eq)
    ss = mf_final.spin_square()[0]
    return new_atoms, float(mf_final.e_tot), bool(mf_final.converged), float(ss)


def stage_geom(sub):
    out = {"substrate": sub, "protocol": f"UKS-{XC.upper()}/{BASIS}, sextet, "
           "2D-pin r(C-H)+r(O-H) (geomeTRIC), frame relaxed", "points": []}
    path = os.path.join(DIR, f"ocm_mnw_{sub}_geom.json")
    atoms = start_atoms(sub)
    grid = [R_POINT] + PATH
    for i, (rCH, rOH) in enumerate(grid):
        t0 = time.time()
        atoms, e, conv, ss = constrained_opt(atoms, rCH, rOH, f"{sub}_{i}")
        out["points"].append({
            "rCH_pin": rCH, "rOH_pin": rOH, "e_h": e, "scf_converged": conv,
            "spin_square": round(ss, 3), "wall_s": round(time.time() - t0, 1),
            "atoms": [[s, [round(x, 6) for x in c]] for s, c in atoms]})
        with open(path, "w") as f:
            json.dump(out, f, indent=1)
        print(f"[{sub}] point {i}/{len(grid)-1} rCH={rCH} rOH={rOH} "
              f"E={e:.6f} conv={conv} <S2>={ss:.2f} ({out['points'][-1]['wall_s']}s)",
              flush=True)
    es = [p["e_h"] for p in out["points"]]
    rel = [(e - es[0]) * HARTREE_KCAL for e in es]
    imax = int(np.argmax(rel[1:]) + 1)
    out["rel_kcal"] = [round(x, 2) for x in rel]
    out["ts_index"] = imax
    out["ts_interior"] = bool(0 < imax < len(grid) - 1)
    out["dft_barrier_kcal"] = round(rel[imax], 2)
    with open(path, "w") as f:
        json.dump(out, f, indent=1)
    print(f"[{sub}] DFT barrier {rel[imax]:.2f} kcal/mol at point {imax} "
          f"(interior={out['ts_interior']})", flush=True)


def analyze_scan(g):
    """R = минимум профиля ДО максимума (истинный пре-комплекс: у CH₄ при
    сближении с O есть яма ниже точки 0), TS = внутренний максимум после него."""
    es = [p["e_h"] for p in g["points"]]
    rel = [(e - es[0]) * HARTREE_KCAL for e in es]
    imax = int(np.argmax(rel[1:]) + 1)
    imin = int(np.argmin(rel[:imax]))
    return rel, imin, imax, bool(0 < imax < len(rel) - 1)


# ---------------------------------------------------------------- энергетика
def rohf(mol):
    mf = scf.ROHF(mol)
    mf.conv_tol = 1e-8
    mfn = scf.newton(mf)
    mfn.kernel()
    if not mfn.converged:
        mf2 = scf.ROHF(mol)
        mf2.level_shift = 0.4
        mf2.max_cycle = 300
        mf2.kernel()
        mfn = scf.newton(mf2)
        mfn.kernel(mf2.make_rdm1())
    return mfn


def forced_avas(mf, target=CAS_TARGET):
    """AVAS с курированными метками; порог подбирается, чтобы (ne,no) совпали
    с целевыми у ВСЕХ структур. Если не выходит — честный отказ (не считаем)."""
    for thr in (0.20, 0.15, 0.12, 0.10, 0.08, 0.06, 0.05, 0.04, 0.03, 0.25, 0.30):
        ncas, nelec, mo = avas.avas(mf, AVAS_LABELS, threshold=thr)
        ne = nelec if isinstance(nelec, (int, np.integer)) else sum(nelec)
        if (int(ne), int(ncas)) == target:
            return mo, thr
    raise RuntimeError(f"AVAS cannot hit CAS{target} (last: {ne}e,{ncas}o)")


def cas_nevpt2(atoms, tag):
    """ROHF → AVAS(единое CAS(9e,10o)) → робастный CASSCF → SC-NEVPT2 + NOON."""
    mol = build_mol(atoms)
    mol.max_memory = 12000
    mf = rohf(mol)
    mo, thr = forced_avas(mf)
    ne, no = CAS_TARGET
    # шаг 1: CASCI + натурбитали — устойчивый старт для CASSCF (рецепт Этапа 14)
    mc0 = mcscf.CASCI(mf, no, ne)
    mc0.fcisolver = fci.direct_spin1.FCI(mol)
    fci.addons.fix_spin_(mc0.fcisolver, shift=0.2, ss=SPIN / 2 * (SPIN / 2 + 1))
    mc0.kernel(mo)
    mo_nat = mc0.cas_natorb()[0]
    # шаг 2: CASSCF с натурбитального старта
    mc = mcscf.CASSCF(mf, no, ne)
    fci.addons.fix_spin_(mc.fcisolver, shift=0.2, ss=SPIN / 2 * (SPIN / 2 + 1))
    mc.max_cycle_macro = 150
    mc.conv_tol = 1e-7
    mc.kernel(mo_nat)
    dm1 = mc.fcisolver.make_rdm1(mc.ci, mc.ncas, mc.nelecas)
    noon = sorted(np.linalg.eigvalsh(dm1), reverse=True)
    e_pt = NEVPT(mc).kernel()          # SC-NEVPT2 корреляционная поправка
    return {
        "tag": tag, "avas_threshold": thr, "cas": list(CAS_TARGET),
        "e_rohf_h": float(mf.e_tot), "rohf_converged": bool(mf.converged),
        "e_casci_h": float(mc0.e_tot),
        "e_casscf_h": float(mc.e_tot), "casscf_converged": bool(mc.converged),
        "e_nevpt2_h": float(mc.e_tot + e_pt), "nevpt2_corr_h": float(e_pt),
        "noon": [round(float(n), 3) for n in noon],
        "n_unpaired_beyond_spin": round(float(
            sum(min(n, 2 - n) for n in noon) - SPIN), 3),
    }


def stage_energy(sub):
    with open(os.path.join(DIR, f"ocm_mnw_{sub}_geom.json")) as f:
        g = json.load(f)
    rel, iR, iTS, interior = analyze_scan(g)
    out = {"substrate": sub, "cas": list(CAS_TARGET), "avas_labels": AVAS_LABELS,
           "r_index": iR, "ts_index": iTS, "ts_interior": interior,
           "rel_kcal": [round(x, 2) for x in rel],
           "dft_barrier_kcal": round(rel[iTS] - rel[iR], 2)}
    path = os.path.join(DIR, f"ocm_mnw_{sub}_energy.json")
    for name, idx in (("R", iR), ("TS", iTS)):
        atoms = [(s, tuple(c)) for s, c in g["points"][idx]["atoms"]]
        t0 = time.time()
        res = cas_nevpt2(atoms, f"{sub}_{name}")
        res["wall_s"] = round(time.time() - t0, 1)
        out[name] = res
        with open(path, "w") as f:
            json.dump(out, f, indent=1)
        print(f"[{sub}] {name}: CASSCF {res['e_casscf_h']:.6f} "
              f"(conv={res['casscf_converged']}) NEVPT2 {res['e_nevpt2_h']:.6f} "
              f"NOON={res['noon']} ({res['wall_s']}s)", flush=True)
    for lvl in ("casscf", "nevpt2"):
        out[f"{lvl}_barrier_kcal"] = round(
            (out["TS"][f"e_{lvl}_h"] - out["R"][f"e_{lvl}_h"]) * HARTREE_KCAL, 2)
    with open(path, "w") as f:
        json.dump(out, f, indent=1)
    print(f"[{sub}] barriers kcal/mol: DFT {out['dft_barrier_kcal']} | "
          f"CASSCF {out['casscf_barrier_kcal']} | NEVPT2 {out['nevpt2_barrier_kcal']}",
          flush=True)


def stage_merge():
    res = {"model": "MnO (sextet) minimal radical-oxygen HAT site, stand-in for "
                    "Mn-Na2WO4/SiO2; NOT the periodic catalyst",
           "protocol": {"geometry": f"UKS-{XC.upper()}/{BASIS}, 2D-pin r(C-H)+r(O-H) "
                                    "relaxed scan (geomeTRIC), sextet",
                        "energy": "ROHF -> AVAS(Mn 3d + O 2p + transferring H 1s) "
                                  "-> CASSCF(9e,10o) [fix_spin, natural-orbital start] "
                                  "-> SC-NEVPT2; SAME curated CAS for all 4 structures",
                        "descriptor": "ddE = barrier(C2H6) - barrier(CH4); "
                                      "<0 over-oxidizes C2 (BEP wall); +4 needed for 70%"},
           "substrates": {}}
    ok = True
    for sub in ("ch4", "c2h6"):
        p = os.path.join(DIR, f"ocm_mnw_{sub}_energy.json")
        if not os.path.exists(p):
            ok = False
            continue
        with open(p) as f:
            res["substrates"][sub] = json.load(f)
    if ok:
        res["ddE_kcal"] = {}
        for lvl in ("dft", "casscf", "nevpt2"):
            b = {s: res["substrates"][s][f"{lvl}_barrier_kcal"] for s in ("ch4", "c2h6")}
            res["ddE_kcal"][lvl] = round(b["c2h6"] - b["ch4"], 2)
        res["quantum_shift_kcal"] = round(
            res["ddE_kcal"]["nevpt2"] - res["ddE_kcal"]["dft"], 2)
        res["all_casscf_converged"] = all(
            res["substrates"][s][pt]["casscf_converged"]
            for s in ("ch4", "c2h6") for pt in ("R", "TS"))
        res["all_ts_interior"] = all(
            res["substrates"][s]["ts_interior"] for s in ("ch4", "c2h6"))
    with open(os.path.join(DIR, "ocm_mnw_selectivity_results.json"), "w") as f:
        json.dump(res, f, indent=1)
    print(json.dumps({k: v for k, v in res.items() if k != "substrates"}, indent=1))


if __name__ == "__main__":
    stage = sys.argv[1]
    if stage == "geom":
        stage_geom(sys.argv[2])
    elif stage == "energy":
        stage_energy(sys.argv[2])
    elif stage == "merge":
        stage_merge()
    else:
        raise SystemExit("usage: geom|energy ch4|c2h6, or merge")
