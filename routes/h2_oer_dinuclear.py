#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
routes/h2_oer_dinuclear.py — Этап 2b водородного трека: ΔG-лестница OER на
ДИНУКЛЕАРНОМ краевом кластере NiOOH — Ni₂(μ-OH)₂(OH)₄ (+X на апикальном Ni).

Зачем: дневник 2026-06-07 честно предупредил, что мононуклеарная модель может
гасить/искажать эффекты решётки: Ni(III) с мостиковыми O — вот где живёт и
сильная корреляция, и реальная OER-активность β/γ-NiOOH. Это минимальный шаг
к решётке: два Ni(III), два мостиковых μ-OH, край как в β-NiOOH.

Схема (паттерн co2_to_fuels_ts — замороженное ядро):
  1) полная оптимизация голого димера (pyberny), скан спина;
  2) адсорбат X = OH*/O*/OOH* на апикальной позиции Ni₁; релаксируются ТОЛЬКО
     DOF адсорбата (scipy L-BFGS по аналитическим градиентам), ядро заморожено
     в геометрии голого димера;
  3) CHE-лестница теми же формулами/поправками, что h2_oer_ladder.py.

ЧЕСТНО: замороженное ядро — жёсткое приближение (нет поверхностной релаксации
под адсорбатом → шаги завышаются); нет растворителя/потенциала; ZPE−TS
литературные. Это скрин тренда моно→ди, не количественный прогноз η.

Запуск:  python3 -u routes/h2_oer_dinuclear.py           # def2-SVP
Резюм:   повторный запуск досчитывает недостающее.
Результат: routes/h2_oer_dinuclear_results.json
"""
import os, sys, json, time, argparse
import numpy as np
from scipy.optimize import minimize
from pyscf import gto, dft
from pyscf.scf import addons

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from h2_oer_ladder import (HARTREE_EV, GCORR_EV, WATER_SPLIT_EV, U_EQ,
                           _scf as _uks_scf)

HERE = os.path.dirname(os.path.abspath(__file__))
RESULTS = os.path.join(HERE, "h2_oer_dinuclear_results.json")
BASIS_DEFAULT = "def2-svp"

# Ni₂(μ-OH)₂(OH)₄, старт (Å): Ni по x, мостики по y, терминальные OH наружу.
CORE = [
    ("Ni",  1.42,  0.00,  0.00),
    ("Ni", -1.42,  0.00,  0.00),
    ("O",   0.00,  1.25,  0.15), ("H", 0.00,  1.95,  0.85),
    ("O",   0.00, -1.25, -0.15), ("H", 0.00, -1.95, -0.85),
    ("O",   2.65,  1.15,  0.35), ("H", 3.35,  1.05,  1.00),
    ("O",   2.65, -1.15, -0.35), ("H", 3.35, -1.05, -1.00),
    ("O",  -2.65,  1.15, -0.35), ("H", -3.35,  1.05, -1.00),
    ("O",  -2.65, -1.15,  0.35), ("H", -3.35, -1.05,  1.00),
]
# Апикальный сайт на Ni₁ (положительный x), старт над Ni вдоль +z:
ADS = {
    "oh":  [("O", 1.42, 0.10, 1.85), ("H", 2.05, 0.45, 2.50)],
    "o":   [("O", 1.42, 0.05, 1.72)],
    "ooh": [("O", 1.42, 0.10, 1.85), ("O", 2.30, 0.75, 2.55),
            ("H", 2.00, 1.60, 2.80)],
}
SPINS = {"bare": (2, 0, 4), "oh": (1, 3), "o": (2, 0, 4), "ooh": (1, 3)}


def _mol(atoms, basis, spin):
    return gto.M(atom=[(s, (x, y, z)) for s, x, y, z in atoms],
                 basis=basis, charge=0, spin=spin, verbose=0, unit="Angstrom")


def opt_core(basis, spin, maxsteps=80):
    """Полная оптимизация голого димера."""
    from pyscf.geomopt.berny_solver import optimize as berny_opt
    mf = _uks_scf(_mol(CORE, basis, spin))
    mol_eq = berny_opt(mf, maxsteps=maxsteps, verbose=0)
    atoms = [(mol_eq.atom_symbol(i), *mol_eq.atom_coords(unit="Angstrom")[i])
             for i in range(mol_eq.natm)]
    mf2 = _uks_scf(_mol(atoms, basis, spin))
    return float(mf2.e_tot), atoms, float(mf2.spin_square()[0])


def relax_ads(core_atoms, ads_name, basis, spin, maxiter=40):
    """Замороженное ядро; L-BFGS только по координатам адсорбата."""
    nc = len(core_atoms)
    ads0 = ADS[ads_name]
    labels = [a[0] for a in ads0]
    x0 = np.array([a[1:] for a in ads0], float).ravel()
    cache = {"dm": None, "best": (np.inf, x0.copy())}

    def build(x):
        pos = x.reshape(-1, 3)
        return core_atoms + [(l, *pos[i]) for i, l in enumerate(labels)]

    def fun(x):
        mf = _uks_scf(_mol(build(x), basis, spin), dm0=cache["dm"])
        cache["dm"] = mf.make_rdm1()
        e = float(mf.e_tot)
        if e < cache["best"][0]:
            cache["best"] = (e, x.copy())
        # PySCF-градиент в Ha/Bohr → Ha/Å (координаты x в Å)
        gA = mf.nuc_grad_method().kernel()[nc:].ravel() / 0.52917721092
        return e, gA

    res = minimize(fun, x0, jac=True, method="L-BFGS-B",
                   options=dict(maxiter=maxiter, ftol=1e-8, gtol=3e-4))
    e_best, x_best = cache["best"]
    return float(min(res.fun, e_best)), build(x_best), bool(res.success)


def load():
    if os.path.exists(RESULTS):
        with open(RESULTS) as f:
            return json.load(f)
    return {"meta": {}, "runs": {}}


def save(d):
    tmp = RESULTS + ".tmp"
    with open(tmp, "w") as f:
        json.dump(d, f, indent=1, ensure_ascii=False)
    os.replace(tmp, RESULTS)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--basis", default=BASIS_DEFAULT)
    args = ap.parse_args()
    basis = args.basis

    data = load()
    data["meta"] = {"model": "Ni2(mu-OH)2(OH)4 edge dimer + apical X on Ni1",
                    "method": f"UKS-PBE(DF)/{basis}; core full-opt (bare), "
                              "adsorbate-only relax (frozen core, L-BFGS)",
                    "honest": "замороженное ядро под адсорбатом; ZPE-TS литературные; "
                              "нет растворителя/потенциала — скрин тренда моно→ди"}

    # 1) голый димер: полный opt по спинам
    for spin in SPINS["bare"]:
        key = f"{basis}:bare:s{spin}"
        if key in data["runs"] and "e_ha" in data["runs"][key]:
            print(f"[skip] {key}"); continue
        print(f"[core] {key} …", flush=True)
        t0 = time.time()
        try:
            e, atoms, ss = opt_core(basis, spin)
            data["runs"][key] = dict(species="bare", spin=spin, e_ha=e,
                                     ss=round(ss, 3), atoms=atoms,
                                     seconds=round(time.time() - t0, 1))
        except Exception as exc:
            data["runs"][key] = dict(species="bare", spin=spin,
                                     error=f"{type(exc).__name__}: {exc}",
                                     seconds=round(time.time() - t0, 1))
        save(data)
        r = data["runs"][key]
        print(f"       E={r.get('e_ha', float('nan'))} ⟨S²⟩={r.get('ss')} "
              f"{r['seconds']}s", flush=True)

    cores = {k: r for k, r in data["runs"].items()
             if r.get("species") == "bare" and "e_ha" in r and k.startswith(basis)}
    if not cores:
        print("[fatal] голый димер не сошёлся ни в одном спине"); return
    best_core = min(cores.values(), key=lambda r: r["e_ha"])
    core_atoms = [tuple(a) for a in best_core["atoms"]]
    print(f"[core] лучший спин 2S={best_core['spin']}  E={best_core['e_ha']:.6f}")

    # 2) адсорбаты на замороженном ядре
    for ads_name in ("oh", "o", "ooh"):
        for spin in SPINS[ads_name]:
            key = f"{basis}:{ads_name}:s{spin}"
            if key in data["runs"] and "e_ha" in data["runs"][key]:
                print(f"[skip] {key}"); continue
            print(f"[ads ] {key} …", flush=True)
            t0 = time.time()
            try:
                e, atoms, ok = relax_ads(core_atoms, ads_name, basis, spin)
                data["runs"][key] = dict(species=ads_name, spin=spin, e_ha=e,
                                         relax_ok=ok, atoms=atoms,
                                         seconds=round(time.time() - t0, 1))
            except Exception as exc:
                data["runs"][key] = dict(species=ads_name, spin=spin,
                                         error=f"{type(exc).__name__}: {exc}",
                                         seconds=round(time.time() - t0, 1))
            save(data)
            r = data["runs"][key]
            print(f"       E={r.get('e_ha', float('nan'))} ok={r.get('relax_ok')} "
                  f"{r['seconds']}s", flush=True)

    # 3) лестница: газовые H₂O/H₂ берём из results мономерного скрипта
    mono = os.path.join(HERE, "h2_oer_ladder_results.json")
    if not os.path.exists(mono):
        print("[warn] нет h2_oer_ladder_results.json — газовые опоры отсутствуют"); return
    with open(mono) as f:
        mruns = json.load(f)["runs"]
    gas = {}
    for k, r in mruns.items():
        b, sp = k.split(":")[0], k.split(":")[1]
        if b == basis and sp in ("h2o", "h2") and "e_ha" in r:
            if sp not in gas or r["e_ha"] < gas[sp]:
                gas[sp] = r["e_ha"]
    best = {}
    for k, r in data["runs"].items():
        if not k.startswith(basis) or "e_ha" not in r:
            continue
        sp = r["species"]
        if sp not in best or r["e_ha"] < best[sp]["e_ha"]:
            best[sp] = r
    need = {"bare", "oh", "o", "ooh"}
    if need <= set(best) and {"h2o", "h2"} <= set(gas):
        g = {sp: best[sp]["e_ha"] * HARTREE_EV + GCORR_EV[sp] for sp in need}
        g["h2o"] = gas["h2o"] * HARTREE_EV + GCORR_EV["h2o"]
        g["h2"] = gas["h2"] * HARTREE_EV + GCORR_EV["h2"]
        h = 0.5 * g["h2"]
        dg1 = g["oh"] + h - g["bare"] - g["h2o"]
        dg2 = g["o"] + h - g["oh"]
        dg3 = g["ooh"] + h - g["o"] - g["h2o"]
        dg4 = WATER_SPLIT_EV - (dg1 + dg2 + dg3)
        steps = {"dG1_bare_to_OH": dg1, "dG2_OH_to_O": dg2,
                 "dG3_O_to_OOH": dg3, "dG4_OOH_to_O2": dg4}
        lim = max(steps, key=steps.get)
        data["ladder"] = {"basis": basis,
                          "steps_eV": {k2: round(v, 3) for k2, v in steps.items()},
                          "limiting": lim, "eta_V": round(steps[lim] - U_EQ, 3),
                          "best_spins": {sp: best[sp]["spin"] for sp in best}}
        save(data)
        print(json.dumps(data["ladder"], indent=1, ensure_ascii=False))
    else:
        print("[warn] лестница не собрана:", sorted(need - set(best)))


if __name__ == "__main__":
    main()
