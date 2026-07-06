#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
routes/field_screen.py — скрин №3 off-label карты: пьезокерамика/внешнее поле как
«электрическое копьё» для C–H. Насколько электрическое поле вдоль оси реакции
снижает барьер H-абстракции?

Физика. По Хеллману–Фейнману dE/dF_z = −μ_z: наклон барьера по полю равен разности
дипольных моментов TS и реагентов вдоль оси. Реагенты (CH4 — ТД-симметричен, X• —
сферичен) почти бездипольны → весь рычаг в диполе TS [CH3···H···X]. Считаем μ_z(TS)
на жёстких TZVP-геометриях TS галогенов + finite-field валидация линейности
(±0.002 а.е.). Перевод: 1 а.е. поля = 514 В/нм; пьезо-горячие точки ~0.1–1 В/нм,
наноплазмоника/острия — до ~10 В/нм.

Выход: ккал/моль на 1 В/нм для Cl/Br/I TS → вердикт идеи №3.
"""
import json
import os

import numpy as np

os.environ.setdefault("HAL_BASIS", "def2-svp")
import importlib.util as ilu

DIR = os.path.dirname(os.path.abspath(__file__))
spec = ilu.spec_from_file_location("hal", os.path.join(DIR, "halogen_descriptor.py"))
hal = ilu.module_from_spec(spec)
spec.loader.exec_module(hal)

AU_FIELD_V_NM = 514.22          # 1 а.е. поля в В/нм
HARTREE_KCAL = 627.509474
AU_DIP_DEBYE = 2.541746


def uks_field(mol, Ez):
    """UKS в постоянном поле Ez (а.е.) вдоль z: стандартный finite-field рецепт."""
    from pyscf import dft, scf
    mf = dft.UKS(mol)
    mf.xc = "pbe0"
    mf.conv_tol = 1e-9
    if abs(Ez) > 0:
        field = np.array([0.0, 0.0, Ez])
        h0 = scf.hf.get_hcore(mol)
        r_ints = mol.intor("int1e_r", comp=3)
        mf.get_hcore = lambda *a: h0 + np.einsum("x,xij->ij", field, r_ints)
        # ядерный член: -Σ Z_A (F·R_A)
        chg = mol.atom_charges()
        crd = mol.atom_coords()
        enuc_shift = -float(np.einsum("i,ix,x->", chg, crd, field))
    else:
        enuc_shift = 0.0
    mf.kernel()
    return float(mf.e_tot) + enuc_shift, mf


def main():
    out = {"screen": "off-label #3: piezo/external E-field lowering of HAT barrier",
           "method": "Hellmann-Feynman slope dE/dFz = -mu_z on rigid TZVP TS "
                     "geometries + finite-field linearity check (+-0.002 au)",
           "ts": []}
    for X in ("Cl", "Br", "I"):
        p = os.path.join(DIR, f"halogen_tzvp_{X}.json")
        if not os.path.exists(p):
            continue
        d = json.load(open(p))
        rCH, rHX = d["dft"]["ts_rCH"], d["dft"]["ts_rHX"]
        atoms = hal.ts_atoms(X, rCH, rHX)
        mol = hal.M(atoms, 1)
        e0, mf0 = uks_field(mol, 0.0)
        mu = mf0.dip_moment(verbose=0)          # Debye, вектор
        mu_z_au = float(mu[2]) / AU_DIP_DEBYE
        # наклон барьера: реагенты бездипольны -> dBarrier/dF = -mu_z(TS)
        slope_kcal_per_Vnm = -mu_z_au * HARTREE_KCAL / AU_FIELD_V_NM
        # finite-field проверка линейности
        ep, _ = uks_field(mol, +0.002)
        em, _ = uks_field(mol, -0.002)
        slope_ff = (ep - em) / 0.004 * HARTREE_KCAL / AU_FIELD_V_NM
        out["ts"].append({
            "X": X, "rCH": rCH, "rHX": rHX,
            "mu_z_debye": round(float(mu[2]), 3),
            "dBarrier_kcal_per_Vnm_HF": round(slope_kcal_per_Vnm, 2),
            "dBarrier_kcal_per_Vnm_ff": round(slope_ff, 2),
            "at_1Vnm_kcal": round(slope_ff, 2),
            "at_0p1Vnm_kcal": round(slope_ff * 0.1, 3)})
        print(f"[{X}] mu_z={mu[2]:.2f} D | dBarrier/dF: HF {slope_kcal_per_Vnm:.2f} "
              f"vs FF {slope_ff:.2f} kcal/mol per V/nm", flush=True)
    strong = [t for t in out["ts"] if abs(t["dBarrier_kcal_per_Vnm_ff"]) > 3]
    out["verdict"] = {
        "reading": ("field lever is REAL at piezo hot-spot scale (>=3 kcal/mol per "
                    "V/nm on TS with aligned dipole)" if strong else
                    "TS dipole too small — piezo lever weak"),
        "honest_caveats": "rigid TS, orientation-averaged real reactor would see "
                          "~1/3 of aligned value; piezo hot-spot fields are local "
                          "and transient"}
    json.dump(out, open(os.path.join(DIR, "field_screen_results.json"), "w"), indent=1)
    print(json.dumps(out["verdict"], ensure_ascii=False))


if __name__ == "__main__":
    main()
