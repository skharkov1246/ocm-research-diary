#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
routes/h2_oer_ladder.py — Этап 2 водородного трека: ΔG-лестница OER
с ОПТИМИЗАЦИЕЙ ГЕОМЕТРИИ (то, чего не хватало предварительной версии 2026-06-08,
скрипт которой был потерян в эфемерном контейнере — восстановлено и закоммичено
ДО запуска, как требует CLAUDE.md).

Модель: мононуклеарный кластер [Ni(OH)₂(X)], X = * (пусто) / OH* / O* / OOH*.
Метод: UKS PBE/def2-SVP + density fitting, аналитические градиенты,
оптимизация геометрии (pyberny), скан спиновых состояний на каждом интермедиате
(берём низший по энергии ПОСЛЕ оптимизации). SCF-робастность: Newton(SOSCF)
первичен, fallback — Fermi-smearing → re-Newton (рецепт co2_to_fuels_ts).

Лестница — computational hydrogen electrode (Nørskov):
  ΔG1: * + H₂O → OH* + (H⁺+e⁻)
  ΔG2: OH*    → O*  + (H⁺+e⁻)
  ΔG3: O* + H₂O → OOH* + (H⁺+e⁻)
  ΔG4: OOH*   → * + O₂ + (H⁺+e⁻) = 4.92 эВ − (ΔG1+ΔG2+ΔG3)  [обход DFT-ошибки O₂]
  η = max(ΔGᵢ)/e − 1.23 В

ЧЕСТНО (границы применимости):
  • ZPE−TS поправки — ЛИТЕРАТУРНЫЕ стандартные (Valdés/Nørskov), не считаны
    для этого кластера (числа в GCORR_EV ниже); электронные энергии — реальные.
  • Мононуклеарный [Ni(OH)₂(X)] — минимальная модель; нет решётки NiOOH,
    растворителя, потенциала электрода. Динуклеарный шаг — h2_oer_dinuclear.py.
  • PBE (GGA) занижает барьеры/зазоры d-оболочки; мультиреференс-поправка
    лимитирующей стадии — Этап 3 (h2_oer_quantum.py).

Запуск:  python3 -u routes/h2_oer_ladder.py            # полный def2-SVP
         python3 -u routes/h2_oer_ladder.py --fast     # смок-тест (STO-3G, без opt)
Резюм:   повторный запуск досчитывает только недостающие (species,spin).
Результат: routes/h2_oer_ladder_results.json
"""
import os, sys, json, time, argparse
import numpy as np
from pyscf import gto, dft
from pyscf.scf import addons

HARTREE_EV = 27.211386245988
HERE = os.path.dirname(os.path.abspath(__file__))
RESULTS = os.path.join(HERE, "h2_oer_ladder_results.json")

# Литературные ZPE−TS поправки (эВ), Valdés et al. JPC C 2008 / Nørskov 2004:
# H₂O(l): ZPE 0.56 − TS 0.67 (0.035 бар ~ жидкость); H₂: 0.27 − 0.41;
# адсорбаты (TS≈0): OH* +0.35, O* +0.05, OOH* +0.41.
GCORR_EV = {"h2o": -0.11, "h2": -0.14, "bare": 0.0,
            "oh": +0.35, "o": +0.05, "ooh": +0.41}
WATER_SPLIT_EV = 4.92          # 2H₂O → O₂ + 2H₂ (эксперимент)
U_EQ = 1.23                    # В, равновесный потенциал OER

# Стартовые геометрии (Å). Слегка несимметричные, чтобы не застрять на седле.
SPECIES = {
    "bare": dict(atoms=[("Ni", 0.0, 0.0, 0.0),
                        ("O",  1.82, 0.15, 0.0), ("H", 2.45, 0.80, 0.15),
                        ("O", -1.82, 0.15, 0.0), ("H", -2.45, 0.80, -0.15)],
                 charge=0, spins=(2, 0)),
    "oh":   dict(atoms=[("Ni", 0.0, 0.0, 0.0),
                        ("O",  1.82, 0.15, 0.0), ("H", 2.45, 0.80, 0.15),
                        ("O", -1.82, 0.15, 0.0), ("H", -2.45, 0.80, -0.15),
                        ("O",  0.10, 1.85, 0.10), ("H", 0.75, 2.45, 0.35)],
                 charge=0, spins=(1, 3)),
    "o":    dict(atoms=[("Ni", 0.0, 0.0, 0.0),
                        ("O",  1.82, 0.15, 0.0), ("H", 2.45, 0.80, 0.15),
                        ("O", -1.82, 0.15, 0.0), ("H", -2.45, 0.80, -0.15),
                        ("O",  0.05, 1.70, 0.05)],
                 charge=0, spins=(0, 2, 4)),
    "ooh":  dict(atoms=[("Ni", 0.0, 0.0, 0.0),
                        ("O",  1.82, 0.15, 0.0), ("H", 2.45, 0.80, 0.15),
                        ("O", -1.82, 0.15, 0.0), ("H", -2.45, 0.80, -0.15),
                        ("O",  0.10, 1.85, 0.10), ("O", 1.05, 2.65, 0.45),
                        ("H",  0.75, 3.50, 0.60)],
                 charge=0, spins=(1, 3)),
    "h2o":  dict(atoms=[("O", 0.0, 0.0, 0.0),
                        ("H", 0.76, 0.59, 0.0), ("H", -0.76, 0.59, 0.0)],
                 charge=0, spins=(0,)),
    "h2":   dict(atoms=[("H", 0.0, 0.0, 0.0), ("H", 0.74, 0.0, 0.0)],
                 charge=0, spins=(0,)),
}


def _mol(species, basis, spin, coords=None):
    sp = SPECIES[species]
    atoms = sp["atoms"] if coords is None else [
        (a[0], *xyz) for a, xyz in zip(sp["atoms"], coords)]
    return gto.M(atom=[(s, (x, y, z)) for s, x, y, z in
                       ((a[0], a[1], a[2], a[3]) for a in atoms)],
                 basis=basis, charge=sp["charge"], spin=spin,
                 verbose=0, unit="Angstrom")


def _scf(mol, dm0=None):
    """UKS-PBE(DF), Newton первичен; fallback: level-shift DIIS → smearing → Newton."""
    mf = dft.UKS(mol).density_fit(); mf.xc = "pbe"; mf.verbose = 0
    mf.conv_tol = 1e-8
    mfn = mf.newton()
    try:
        mfn.kernel(dm0=dm0)
    except Exception:
        mfn.converged = False
    if not getattr(mfn, "converged", False):
        mf2 = dft.UKS(mol).density_fit(); mf2.xc = "pbe"; mf2.verbose = 0
        mf2.level_shift = 0.4; mf2.max_cycle = 300; mf2.conv_tol = 1e-6
        mf2 = addons.smearing_(mf2, sigma=0.01, method="fermi")
        mf2.kernel(dm0=dm0)
        mfn = mf2.newton(); mfn.kernel(mf2.make_rdm1())
    return mfn


def optimize_species(species, spin, basis, do_opt=True, maxsteps=60):
    """Оптимизация геометрии (pyberny) для (species, spin). Возвращает dict."""
    t0 = time.time()
    mol = _mol(species, basis, spin)
    mf = _scf(mol)
    e0 = float(mf.e_tot)
    coords = mol.atom_coords(unit="Angstrom").tolist()
    converged_opt = None
    if do_opt and mol.natm > 2:
        from pyscf.geomopt.berny_solver import optimize as berny_opt
        try:
            mol_eq = berny_opt(mf, maxsteps=maxsteps, verbose=0)
            mf = _scf(gto.M(atom=list(zip([mol.atom_symbol(i) for i in range(mol.natm)],
                                          mol_eq.atom_coords(unit="Angstrom"))),
                            basis=basis, charge=mol.charge, spin=spin,
                            verbose=0, unit="Angstrom"))
            coords = mol_eq.atom_coords(unit="Angstrom").tolist()
            converged_opt = True
        except Exception as exc:            # честно фиксируем неудачу оптимизации
            converged_opt = f"FAILED: {type(exc).__name__}: {exc}"
    elif do_opt and mol.natm == 2 and species == "h2":
        # H₂: явный 1D-скан связи (berny не любит двухатомники)
        best = (None, None)
        for d in np.arange(0.62, 0.92, 0.02):
            m = gto.M(atom=[("H", (0, 0, 0)), ("H", (float(d), 0, 0))],
                      basis=basis, spin=0, verbose=0, unit="Angstrom")
            e = float(_scf(m).e_tot)
            if best[0] is None or e < best[0]:
                best = (e, float(d))
        e0 = best[0]
        coords = [[0, 0, 0], [best[1], 0, 0]]
        converged_opt = True
        return dict(species=species, spin=spin, e_ha=e0, opt=converged_opt,
                    ss=0.0, coords=coords, seconds=round(time.time() - t0, 1))
    e_final = float(mf.e_tot)
    ss = float(mf.spin_square()[0]) if spin >= 0 else 0.0
    return dict(species=species, spin=spin,
                e_ha=e_final if converged_opt in (True, None) else e0,
                opt=converged_opt, ss=round(ss, 3), coords=coords,
                seconds=round(time.time() - t0, 1))


def load_results():
    if os.path.exists(RESULTS):
        with open(RESULTS) as f:
            return json.load(f)
    return {"meta": {}, "runs": {}}


def save_results(data):
    tmp = RESULTS + ".tmp"
    with open(tmp, "w") as f:
        json.dump(data, f, indent=1, ensure_ascii=False)
    os.replace(tmp, RESULTS)


def ladder(best):
    """CHE-лестница из лучших (низших) G каждой частицы. Энергии в эВ."""
    g = {k: best[k]["e_ha"] * HARTREE_EV + GCORR_EV[k] for k in best}
    h = 0.5 * g["h2"]
    dg1 = g["oh"] + h - g["bare"] - g["h2o"]
    dg2 = g["o"] + h - g["oh"]
    dg3 = g["ooh"] + h - g["o"] - g["h2o"]
    dg4 = WATER_SPLIT_EV - (dg1 + dg2 + dg3)
    steps = {"dG1_bare_to_OH": dg1, "dG2_OH_to_O": dg2,
             "dG3_O_to_OOH": dg3, "dG4_OOH_to_O2": dg4}
    lim = max(steps, key=steps.get)
    return {"steps_eV": {k: round(v, 3) for k, v in steps.items()},
            "limiting": lim, "eta_V": round(steps[lim] - U_EQ, 3),
            "sum_check_eV": round(dg1 + dg2 + dg3 + dg4, 3)}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fast", action="store_true",
                    help="смок: STO-3G, без оптимизации геометрии")
    ap.add_argument("--basis", default=None)
    ap.add_argument("--species", default=None, help="только одна частица")
    args = ap.parse_args()
    basis = args.basis or ("sto-3g" if args.fast else "def2-svp")
    do_opt = not args.fast

    data = load_results()
    data["meta"] = {"model": "[Ni(OH)2(X)] mononuclear",
                    "method": f"UKS-PBE(DF)/{basis}"
                              + (" + berny geom-opt" if do_opt else " single-point"),
                    "gcorr_eV": GCORR_EV,
                    "honest": "ZPE-TS поправки литературные; кластер минимальный; "
                              "нет решётки/растворителя/потенциала"}
    todo = [args.species] if args.species else list(SPECIES)
    for sp in todo:
        for spin in SPECIES[sp]["spins"]:
            key = f"{basis}:{sp}:s{spin}" + ("" if do_opt else ":sp")
            if key in data["runs"] and data["runs"][key].get("opt") in (True, None):
                print(f"[skip] {key} (готово)"); continue
            print(f"[run ] {key} …", flush=True)
            try:
                res = optimize_species(sp, spin, basis, do_opt=do_opt)
            except Exception as exc:
                res = dict(species=sp, spin=spin, error=f"{type(exc).__name__}: {exc}")
            data["runs"][key] = res
            save_results(data)
            print(f"       E={res.get('e_ha', float('nan')):.6f} Ha  "
                  f"opt={res.get('opt')}  ⟨S²⟩={res.get('ss')}  "
                  f"{res.get('seconds', '?')}s", flush=True)

    # лучшие (низшие) по спину — только успешные
    best = {}
    for key, r in data["runs"].items():
        b, sp = key.split(":")[0], key.split(":")[1]
        if b != basis or "e_ha" not in r:
            continue
        if do_opt and r.get("opt") not in (True, None):
            continue
        if sp not in best or r["e_ha"] < best[sp]["e_ha"]:
            best[sp] = r
    if all(k in best for k in ("bare", "oh", "o", "ooh", "h2o", "h2")):
        lad = ladder(best)
        data["ladder"] = {"basis": basis, "geom_opt": do_opt, **lad,
                          "best_spins": {k: best[k]["spin"] for k in best}}
        save_results(data)
        print(json.dumps(data["ladder"], indent=1, ensure_ascii=False))
    else:
        print("[warn] лестница не собрана — не все частицы посчитаны:",
              sorted(set(SPECIES) - set(best)))


if __name__ == "__main__":
    main()
