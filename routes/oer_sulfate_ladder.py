#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
routes/oer_sulfate_ladder.py — гейт G1 кейса N1 (анод электроэкстракции
Норникеля, ТЗ: неблагородное покрытие <= $390/м²).

Вопрос гейта: как СУЛЬФАТНАЯ координация (кислый H₂SO₄-электролит цеха
электроэкстракции) меняет OER-перенапряжение η на Ni-центре относительно
щелочного якоря routes/h2_oer_ladder.py (η = 0.435 В)?

Модель (мононуклеарная, как в образце h2_oer_ladder.py):
  • Сульфатная ветка:  [Ni(OH)(HSO4)(X)] — один OH-лиганд якорной модели
    заменён на бисульфат HSO4⁻ (O-монодентатный), X = * / OH* / O* / OOH*
    на апикальной позиции. Заряды/спины пересчитаны честно:
      sulf_bare 86 e⁻ (чёт)  → spins (2, 0)
      sulf_oh   95 e⁻ (нечёт)→ spins (1, 3)
      sulf_o    94 e⁻ (чёт)  → spins (0, 2, 4)
      sulf_ooh 103 e⁻ (нечёт)→ spins (1, 3)
    Все частицы нейтральны: Ni²⁺ + OH⁻ + HSO4⁻ = 0.
  • Контрольная ветка: [Ni(OH)₂(X)] — ТА ЖЕ частица, что в h2_oer_ladder
    (якорь старой лестницы), пересчитанная тем же кодом/базисом/XC, чтобы
    изолировать эффект сульфатного лиганда: Δη = η_sulf − η_ctrl.

Метод: UKS PBE(DF)/def2-SVP (env-ручки ниже), pyberny-оптимизация геометрии,
скан спиновых состояний, gradient-audit (gmax в JSON), SCF-каскад
Newton → level-shift DIIS → smearing → re-DIIS (рецепт h2_oer_ladder).

CHE-лестница (Nørskov) — в RHE-шкале НЕ зависит от pH, поэтому кислая среда
не меняет формулы шагов; эффект среды входит только через лиганд HSO4⁻:
  ΔG1: * + H₂O → OH* + (H⁺+e⁻)
  ΔG2: OH*     → O*  + (H⁺+e⁻)
  ΔG3: O* + H₂O → OOH* + (H⁺+e⁻)
  ΔG4: OOH*    → * + O₂ + (H⁺+e⁻) = 4.92 эВ − (ΔG1+ΔG2+ΔG3)
  η = max(ΔGᵢ)/e − 1.23 В

Честные оговорки (computed_here в meta):
  • электронные энергии и геометрии — реальный расчёт (computed_here: true);
  • ZPE−TS поправки — литературные (Valdés/Nørskov), НЕ считаны для этих
    кластеров (computed_here: false); для HSO4-лиганда допущение: поправки
    адсорбатов OH*/O*/OOH* не зависят от спектаторного лиганда;
  • замыкание 4.92 эВ — эксперимент (computed_here: false);
  • нет явного растворителя/потенциала/решётки; PBE (GGA) — известные
    ограничения по d-оболочке Ni.

env-ручки:
  OER_SPECIES        подмножество частиц через запятую (default: все)
  OER_BASIS          базис (default def2-svp)
  OER_XC             функционал pyscf (default "pbe,pbe")
  OER_STAGE_TIMEOUT  лимит секунд на (species,spin)-вариант (default 5400)
  OER_JSON_SUFFIX    суффикс файла результатов ("_sulf"/"_ctrl" при сплите)

Запуск:  python3 -u routes/oer_sulfate_ladder.py            # полный
         python3 -u routes/oer_sulfate_ladder.py --fast     # смок (STO-3G, без opt)
Резюм:   повторный запуск досчитывает только недостающие (species,spin) —
         каждая сошедшаяся частица немедленно пишется в results-JSON.
Результат: routes/oer_sulfate_ladder_results{OER_JSON_SUFFIX}.json
"""
import os, sys, json, time, signal, argparse

import numpy as np
try:
    from pyscf import gto, dft
    from pyscf.scf import addons
except ImportError as _exc:                       # в сандбоксе pyscf может отсутствовать
    sys.exit("[oer-sulfate] pyscf не установлен (%s). Скрипт предназначен для "
             "AWS-запуска через routes/oer_sulfate_aws.py; локально: "
             "pip install numpy scipy pyscf pyberny" % _exc)

HARTREE_EV = 27.211386245988
HERE = os.path.dirname(os.path.abspath(__file__))
SUFFIX = os.environ.get("OER_JSON_SUFFIX", "")
RESULTS = os.path.join(HERE, f"oer_sulfate_ladder_results{SUFFIX}.json")
STAGE_TIMEOUT = int(os.environ.get("OER_STAGE_TIMEOUT", "5400"))  # с на (species,spin)

# Литературные ZPE−TS поправки (эВ), Valdés et al. JPC C 2008 / Nørskov 2004.
# computed_here: false — НЕ считаны для этих кластеров. Допущение: поправки
# адсорбатов переносимы между [Ni(OH)2] и [Ni(OH)(HSO4)] (спектаторный лиганд).
GCORR_ADS_EV = {"bare": 0.0, "oh": +0.35, "o": +0.05, "ooh": +0.41}
GCORR_REF_EV = {"h2o": -0.11, "h2": -0.14}
WATER_SPLIT_EV = 4.92          # 2H₂O → O₂ + 2H₂, эксперимент
U_EQ = 1.23                    # В, равновесный потенциал OER

# --- Стартовые геометрии (Å), слегка несимметричные, чтобы не сесть на седло ---

# Контрольный остов [Ni(OH)2] — ровно как в h2_oer_ladder.py (якорь).
_CTRL_CORE = [("Ni", 0.0, 0.0, 0.0),
              ("O",  1.82, 0.15, 0.0), ("H", 2.45, 0.80, 0.15),
              ("O", -1.82, 0.15, 0.0), ("H", -2.45, 0.80, -0.15)]

# Сульфатный остов [Ni(OH)(HSO4)]: OH сохранён, второй OH заменён на
# O-монодентатный бисульфат (Ni–O ~1.85 Å, S–O(Ni) ~1.51, S=O ~1.45/1.52,
# S–OH ~1.61, O–H ~0.95).
_SULF_CORE = [("Ni", 0.0, 0.0, 0.0),
              ("O",  1.82, 0.15, 0.0), ("H", 2.45, 0.80, 0.15),   # OH⁻
              ("O", -1.85, 0.15, 0.0),                            # O(Ni)–SO3H
              ("S", -3.35, 0.25, 0.15),
              ("O", -4.05, -0.95, -0.30),                         # S=O
              ("O", -4.30, 1.35, 0.60),                           # S=O
              ("O", -2.85, 0.55, 1.65), ("H", -3.45, 0.35, 2.35)] # S–O–H

# Апикальные адсорбаты — те же стартовые позиции, что в образце.
_ADS = {"bare": [],
        "oh":   [("O", 0.10, 1.85, 0.10), ("H", 0.75, 2.45, 0.35)],
        "o":    [("O", 0.05, 1.70, 0.05)],
        "ooh":  [("O", 0.10, 1.85, 0.10), ("O", 1.05, 2.65, 0.45),
                 ("H", 0.75, 3.50, 0.60)]}

# Спины (2S = Nα−Nβ): чётность по числу электронов, набор — как в образце.
_SPINS_EVEN_BARE = (2, 0)      # Ni(II) d8: триплет обычно ниже
_SPINS_ODD = (1, 3)
_SPINS_EVEN_O = (0, 2, 4)

SPECIES = {}
for _ads, _spins in (("bare", _SPINS_EVEN_BARE), ("oh", _SPINS_ODD),
                     ("o", _SPINS_EVEN_O), ("ooh", _SPINS_ODD)):
    SPECIES[f"sulf_{_ads}"] = dict(atoms=_SULF_CORE + _ADS[_ads],
                                   charge=0, spins=_spins)
    SPECIES[f"ctrl_{_ads}"] = dict(atoms=_CTRL_CORE + _ADS[_ads],
                                   charge=0, spins=_spins)
SPECIES["h2o"] = dict(atoms=[("O", 0.0, 0.0, 0.0),
                             ("H", 0.76, 0.59, 0.0), ("H", -0.76, 0.59, 0.0)],
                      charge=0, spins=(0,))
SPECIES["h2"] = dict(atoms=[("H", 0.0, 0.0, 0.0), ("H", 0.74, 0.0, 0.0)],
                     charge=0, spins=(0,))

BRANCHES = {"sulf": "[Ni(OH)(HSO4)(X)] — сульфатная кислая среда",
            "ctrl": "[Ni(OH)2(X)] — контроль без сульфата (якорь h2_oer_ladder)"}


def _mol(species, basis, spin):
    sp = SPECIES[species]
    return gto.M(atom=[(a[0], (a[1], a[2], a[3])) for a in sp["atoms"]],
                 basis=basis, charge=sp["charge"], spin=spin,
                 verbose=0, unit="Angstrom")


def _scf(mol, xc, dm0=None):
    """UKS(DF), Newton первичен; fallback: level-shift DIIS →
    smearing-DIIS → warm-start level-shift DIIS. Незаcошедшийся SCF —
    ошибка, а не тихо возвращённая энергия (рецепт h2_oer_ladder)."""
    mf = dft.UKS(mol).density_fit(); mf.xc = xc; mf.verbose = 0
    mf.conv_tol = 1e-8
    mfn = mf.newton()
    try:
        mfn.kernel(dm0=dm0)
    except Exception:
        mfn.converged = False
    if getattr(mfn, "converged", False):
        return mfn
    mf2 = dft.UKS(mol).density_fit(); mf2.xc = xc; mf2.verbose = 0
    mf2.level_shift = 0.3; mf2.max_cycle = 400; mf2.conv_tol = 1e-8
    mf2.kernel(dm0=dm0)
    if mf2.converged:
        return mf2
    mfs = dft.UKS(mol).density_fit(); mfs.xc = xc; mfs.verbose = 0
    mfs = addons.smearing_(mfs, sigma=0.01, method="fermi")
    mfs.max_cycle = 400; mfs.conv_tol = 1e-6
    mfs.kernel(dm0=dm0)
    mf3 = dft.UKS(mol).density_fit(); mf3.xc = xc; mf3.verbose = 0
    mf3.level_shift = 0.3; mf3.max_cycle = 400; mf3.conv_tol = 1e-8
    mf3.kernel(dm0=mfs.make_rdm1())
    if not mf3.converged:
        raise RuntimeError("SCF unconverged after Newton/level-shift/smearing")
    return mf3


GMAX_TOL = 6e-4   # Ha/Bohr — gradient-audit: max-компонента градиента в финале


class StageTimeout(Exception):
    pass


def _alarm(_sig, _frm):
    raise StageTimeout(f"OER_STAGE_TIMEOUT={STAGE_TIMEOUT}s exceeded")


def optimize_species(species, spin, basis, xc, do_opt=True, maxsteps=60,
                     restarts=2):
    """Оптимизация геометрии (pyberny) для (species, spin). Сходимость
    проверяется ПРЯМО — нормой градиента в финальной точке (gradient-audit:
    gmax уходит в JSON); до `restarts` продолжений."""
    t0 = time.time()
    mol = _mol(species, basis, spin)
    mf = _scf(mol, xc)
    e0 = float(mf.e_tot)
    coords = mol.atom_coords(unit="Angstrom").tolist()
    converged_opt, gmax = None, None
    if do_opt and species == "h2":
        # H₂: 1D-скан + доводка (berny не любит двухатомники) — как в образце
        def e_at(d):
            m = gto.M(atom=[("H", (0, 0, 0)), ("H", (float(d), 0, 0))],
                      basis=basis, spin=0, verbose=0, unit="Angstrom")
            return float(_scf(m, xc).e_tot)
        coarse = {round(float(d), 3): e_at(d) for d in np.arange(0.62, 0.92, 0.02)}
        d0 = min(coarse, key=coarse.get)
        fine = {round(float(d), 4): e_at(d)
                for d in np.arange(d0 - 0.02, d0 + 0.021, 0.004)}
        d1 = min(fine, key=fine.get)
        return dict(species=species, spin=spin, e_ha=fine[d1], opt=True,
                    gmax=0.0, ss=0.0, coords=[[0, 0, 0], [d1, 0, 0]],
                    geom_opt=True, seconds=round(time.time() - t0, 1))
    if do_opt:
        from pyscf.geomopt.berny_solver import optimize as berny_opt
        try:
            cur_mf = mf
            for _attempt in range(restarts + 1):
                mol_eq = berny_opt(cur_mf, maxsteps=maxsteps, verbose=0)
                atoms = list(zip([mol.atom_symbol(i) for i in range(mol.natm)],
                                 mol_eq.atom_coords(unit="Angstrom")))
                cur_mf = _scf(gto.M(atom=atoms, basis=basis, charge=mol.charge,
                                    spin=spin, verbose=0, unit="Angstrom"), xc)
                gmax = float(np.abs(cur_mf.nuc_grad_method().kernel()).max())
                if gmax < GMAX_TOL:
                    break
            mf = cur_mf
            coords = mol_eq.atom_coords(unit="Angstrom").tolist()
            converged_opt = True if gmax < GMAX_TOL else \
                f"UNCONVERGED gmax={gmax:.5f} after {restarts+1} berny runs"
        except StageTimeout:
            raise                                # таймаут фиксируем выше, честно
        except Exception as exc:                 # неудачу оптимизации фиксируем
            converged_opt = f"FAILED: {type(exc).__name__}: {exc}"
    e_final = float(mf.e_tot)
    ss = float(mf.spin_square()[0]) if spin >= 0 else 0.0
    failed = isinstance(converged_opt, str) and converged_opt.startswith("FAILED")
    return dict(species=species, spin=spin,
                e_ha=e0 if failed else e_final,
                opt=converged_opt, gmax=gmax, ss=round(ss, 3), coords=coords,
                geom_opt=bool(do_opt),
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


def _gcorr(name):
    """ZPE−TS поправка (эВ) по имени частицы: sulf_oh/ctrl_oh → oh и т.д."""
    if name in GCORR_REF_EV:
        return GCORR_REF_EV[name]
    return GCORR_ADS_EV[name.split("_", 1)[1]]


def ladder(best, branch):
    """CHE-лестница ветки branch из лучших (низших) G. RHE-шкала: шаги
    не зависят от pH — кислая среда входит только через лиганд HSO4⁻."""
    need = [f"{branch}_bare", f"{branch}_oh", f"{branch}_o", f"{branch}_ooh",
            "h2o", "h2"]
    if not all(k in best for k in need):
        return None
    g = {k: best[k]["e_ha"] * HARTREE_EV + _gcorr(k) for k in need}
    h = 0.5 * g["h2"]
    dg1 = g[f"{branch}_oh"] + h - g[f"{branch}_bare"] - g["h2o"]
    dg2 = g[f"{branch}_o"] + h - g[f"{branch}_oh"]
    dg3 = g[f"{branch}_ooh"] + h - g[f"{branch}_o"] - g["h2o"]
    dg4 = WATER_SPLIT_EV - (dg1 + dg2 + dg3)
    steps = {"dG1_bare_to_OH": dg1, "dG2_OH_to_O": dg2,
             "dG3_O_to_OOH": dg3, "dG4_OOH_to_O2": dg4}
    lim = max(steps, key=steps.get)
    return {"branch": branch, "model": BRANCHES[branch],
            "steps_eV": {k: round(v, 3) for k, v in steps.items()},
            "limiting": lim, "eta_V": round(steps[lim] - U_EQ, 3),
            "sum_check_eV": round(dg1 + dg2 + dg3 + dg4, 3),
            "best_spins": {k: best[k]["spin"] for k in need}}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fast", action="store_true",
                    help="смок: STO-3G, без оптимизации геометрии")
    ap.add_argument("--basis", default=None)
    ap.add_argument("--species", default=None,
                    help="подмножество через запятую (приоритет над OER_SPECIES)")
    args = ap.parse_args()
    basis = args.basis or os.environ.get("OER_BASIS") or \
        ("sto-3g" if args.fast else "def2-svp")
    xc = os.environ.get("OER_XC", "pbe,pbe")
    do_opt = not args.fast

    sel = args.species or os.environ.get("OER_SPECIES") or ""
    todo = [s.strip() for s in sel.split(",") if s.strip()] or list(SPECIES)
    unknown = [s for s in todo if s not in SPECIES]
    if unknown:
        sys.exit(f"[oer-sulfate] неизвестные species: {unknown}; "
                 f"допустимые: {sorted(SPECIES)}")

    data = load_results()
    data["meta"] = {
        "purpose": "гейт G1 N1: влияние сульфатной координации на η; "
                   "ТЗ покрытие <=$390/м2",
        "case": "N1 анод электроэкстракции (H2SO4, кислая среда), "
                "неблагородное покрытие",
        "branches": BRANCHES,
        "method": f"UKS-{xc}(DF)/{basis}"
                  + (" + berny geom-opt" if do_opt else " single-point"),
        "che_note": "CHE в RHE-шкале не зависит от pH: формулы шагов те же, "
                    "что в щелочном якоре; эффект среды — через лиганд HSO4-",
        "gcorr_eV": {**GCORR_REF_EV, **GCORR_ADS_EV},
        "computed_here": {"electronic_energies": True, "geometries": True,
                          "gmax_gradient_audit": True,
                          "zpe_ts_gcorr": False,        # Valdés/Nørskov, лит.
                          "water_split_4.92eV": False,  # эксперимент
                          "anchor_eta_alkaline_0.435V": False},  # h2_oer_ladder
        "stage_timeout_s": STAGE_TIMEOUT,
        "honest": "ZPE-TS литературные и приняты переносимыми на HSO4-кластер; "
                  "мононуклеарная модель; нет решётки/растворителя/потенциала; "
                  "PBE занижает d-зазоры Ni"}
    save_results(data)

    have_alarm = hasattr(signal, "SIGALRM")
    if have_alarm:
        signal.signal(signal.SIGALRM, _alarm)
    for sp in todo:
        for spin in SPECIES[sp]["spins"]:
            key = f"{basis}:{xc}:{sp}:s{spin}" + ("" if do_opt else ":sp")
            prev = data["runs"].get(key, {})
            # резюм: пропускаем только УСПЕШНЫЕ записи (ошибки пересчитываем)
            if "e_ha" in prev and prev.get("opt") in (True, None):
                print(f"[skip] {key} (готово)"); continue
            print(f"[run ] {key} (timeout {STAGE_TIMEOUT}s) …", flush=True)
            if have_alarm:
                signal.alarm(STAGE_TIMEOUT)
            try:
                res = optimize_species(sp, spin, basis, xc, do_opt=do_opt)
            except StageTimeout as exc:
                res = dict(species=sp, spin=spin, timeout=True,
                           error=str(exc))
            except Exception as exc:
                res = dict(species=sp, spin=spin,
                           error=f"{type(exc).__name__}: {exc}")
            finally:
                if have_alarm:
                    signal.alarm(0)
            # поточечный чекпойнт: каждая частица пишется НЕМЕДЛЕННО
            data["runs"][key] = res
            save_results(data)
            print(f"       E={res.get('e_ha', float('nan')):.6f} Ha  "
                  f"opt={res.get('opt')}  gmax={res.get('gmax')}  "
                  f"⟨S²⟩={res.get('ss')}  {res.get('seconds', '?')}s",
                  flush=True)

    # лучшие (низшие) по спину — только успешные И того же режима (opt vs :sp)
    best = {}
    for key, r in data["runs"].items():
        parts = key.split(":")
        if len(parts) < 4:
            continue
        b, kxc, sp, is_sp = parts[0], parts[1], parts[2], key.endswith(":sp")
        if b != basis or kxc != xc or "e_ha" not in r or is_sp == do_opt:
            continue
        if do_opt and r.get("opt") is not True:   # UNCONVERGED/FAILED — вон
            continue
        if sp not in best or r["e_ha"] < best[sp]["e_ha"]:
            best[sp] = r

    ladders = {}
    for branch in BRANCHES:
        lad = ladder(best, branch)
        if lad is None:
            missing = sorted(k for k in
                             (f"{branch}_bare", f"{branch}_oh", f"{branch}_o",
                              f"{branch}_ooh", "h2o", "h2") if k not in best)
            print(f"[warn] ветка {branch}: лестница не собрана, "
                  f"не хватает {missing}")
            continue
        ladders[branch] = lad
        print(json.dumps(lad, indent=1, ensure_ascii=False))
    if ladders:
        entry = {"basis": basis, "xc": xc, "geom_opt": do_opt, **ladders}
        if "sulf" in ladders and "ctrl" in ladders:
            entry["delta_eta_sulf_minus_ctrl_V"] = round(
                ladders["sulf"]["eta_V"] - ladders["ctrl"]["eta_V"], 3)
        # изоляция режимов: --fast НЕ трогает продакшен-слот "ladder"
        data.setdefault("ladders", {})[f"{basis}:{xc}"
                                       + ("" if do_opt else ":sp")] = entry
        if do_opt:
            data["ladder"] = entry
        save_results(data)


if __name__ == "__main__":
    main()
