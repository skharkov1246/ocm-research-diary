#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""routes/nh3_nitride_ladder.py — ШОРТКАТ-ЛЕСТНИЦА (ammonia): химический LOOPING
вместо Габера-Боша. Вулкан-дескриптор двухшагового цикла на металл-мотивах.

ФЛИП технологии. Габер-Бош делает NH₃ в ОДНУ стадию — диссоциативная хемосорбция
N₂ на Fe/Ru при 150-300 атм и 400-500°C (весь бюджет — на разрыв N≡N в жёстких
условиях). Химический looping разбивает синтез на ДВА шага, каждый — мягче:
  (1) нитридизация:   M + ½N₂ → M–N   (реактивный металл рвёт N₂; можно вести
                                        электро-/механо-/при низкой T);
  (2) гидрирование:   M–N + 3⁄2 H₂ → M + NH₃  (нитрид отдаёт аммиак, металл
                                        регенерируется).
  Нетто:  ½N₂ + 3⁄2H₂ → NH₃  (то же, что ГБ), но T/P каждого шага РАЗВЯЗАНЫ.
Литературный якорь: chemical-looping ammonia synthesis (Li₃N-, Mn-нитрид-,
щёлочно-земельные imide/nitride циклы) — реально работает, вопрос в ВЫБОРЕ
металла: слишком стабильный нитрид (Li,Ca) жадно ловит N₂, но не отдаёт NH₃ без
высокой T; слишком слабый — плохо ловит N₂. Оптимум Сабатье — где ОБА шага
близки к термонейтральности (мягкие условия в обе стороны = дешевле всего).

Дескрипторы на металл M (мотив M₄ + μ-N):
  ΔE_nit = E(M₄N) − E(M₄) − ½E(N₂)              (фиксация N₂; хотим ≲0, но не −∞)
  ΔE_hyd = E(M₄) + E(NH₃) − E(M₄N) − 3⁄2 E(H₂)  (выпуск NH₃; хотим ≲0)
  balance = max(|ΔE_nit − C/2|, |ΔE_hyd − C/2|) → МИНИМИЗИРУЕМ (вулкан).

ВСТРОЕННЫЙ САМОТЕСТ (замыкание цикла): для ЛЮБОГО металла
  ΔE_nit + ΔE_hyd = E(NH₃) − ½E(N₂) − 3⁄2E(H₂) = C  (энергия образования NH₃,
НЕ зависит от металла). Считаем C напрямую из молекул-эталонов и требуем, чтобы
сумма двух шагов совпала с C (±0.05 эВ). Если не совпала — геометрия/спин
недосошлись, точку метим невалидной. Это делает вулкан честным: цена ошибки
в одном шаге видна как разъезд с C.

ЧЕСТНО: M₄-КЛАСТЕР, не решётка — числа суть ЭНЕРГЕТИКА МОТИВА (N на грани из 4
атомов), физический смысл в ТРЕНДЕ по металлам, не в абсолютной энтальпии
образования нитрида; газовая фаза, UKS-PBE/def2-SVP+DF, скан мультиплетности
(берём низшую по E); без ZPE/энтропии/растворителя. Проба отвечает на один
вопрос: КАКОЙ металл даёт наиболее сбалансированный (=мягкий, =дешёвый) цикл.

Параллелится по металлу: METAL=Fe python3 -u routes/nh3_nitride_ladder.py
Металлы:  Li, Mg, Ca, Mn, Fe (по env METAL; по умолчанию все по очереди).
Smoke:    --smoke (только Li + эталоны, грубо; проверка замыкания цикла, ~5 мин)
Выход:    routes/nh3_nitride_<metal>_results.json (atomic, incremental)
"""
import os, sys, json, time, argparse, tempfile
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
from pyscf import gto, dft, lib  # noqa: E402
from pyscf.geomopt.geometric_solver import optimize  # noqa: E402

EV = 27.211386245988

# металл-специфичные стартовые геометрии и сетки спинов (2S на весь кластер/мотив)
# рёбра M–M и M–N — грубые лит.-оценки; ГЕОМЕТРИЯ РЕЛАКСИРУЕТСЯ, старт лишь затравка.
METALS = {
    "Li": {"a_mm": 3.0, "a_mn": 1.9, "spins_m4": [0, 2], "spins_m4n": [1, 3]},
    "Mg": {"a_mm": 3.2, "a_mn": 2.1, "spins_m4": [0, 2], "spins_m4n": [1, 3]},
    "Ca": {"a_mm": 3.9, "a_mn": 2.4, "spins_m4": [0, 2], "spins_m4n": [1, 3]},
    "Mn": {"a_mm": 2.7, "a_mn": 2.0, "spins_m4": [0, 4, 8, 12, 16, 20],
           "spins_m4n": [1, 5, 9, 13, 17, 21]},
    "Fe": {"a_mm": 2.5, "a_mn": 1.9, "spins_m4": [0, 4, 8, 12, 16],
           "spins_m4n": [1, 5, 9, 13, 17]},
}

def say(m): print(f"[{time.strftime('%H:%M:%S')}] {m}", flush=True)

def outpath(metal): return os.path.join(HERE, f"nh3_nitride_{metal.lower()}_results.json")

def atomic_save(obj, path):
    fd, tmp = tempfile.mkstemp(dir=HERE, suffix=".tmp")
    with os.fdopen(fd, "w") as f:
        json.dump(obj, f, ensure_ascii=False, indent=1)
    os.replace(tmp, path)

def tetra(a):
    """4 вершины правильного тетраэдра с ребром a, центр в нуле."""
    v = np.array([[1, 1, 1], [1, -1, -1], [-1, 1, -1], [-1, -1, 1]], float)
    v *= a / np.sqrt(8.0)                      # ребро тетраэдра из этих вершин = a
    return v

def build_m4(metal, a):
    return [[metal, tuple(p)] for p in tetra(a)]

def build_m4n(metal, a_mm, a_mn):
    # N в центроиде тетраэдра (μ₄); релаксация уведёт его к грани если так глубже
    atoms = [[metal, tuple(p)] for p in tetra(a_mm)]
    atoms.append(["N", (0.0, 0.0, 0.0)])
    return atoms

def mk_mol(atoms, spin, charge=0):
    return gto.M(atom=atoms, basis="def2-svp", spin=spin, charge=charge, verbose=0)

def mkmf(mol):
    mf = dft.UKS(mol).density_fit()
    mf.xc = "pbe"; mf.conv_tol = 1e-9; mf.max_cycle = 200; mf.level_shift = 0.2
    return mf

def scf_e(mol):
    mf = mkmf(mol)
    e = mf.kernel()
    if not mf.converged:
        mfn = mf.newton(); mfn.max_cycle = 80
        e = mfn.kernel(mf.mo_coeff, mf.mo_occ)
        mf.converged = bool(mfn.converged); e = float(mfn.e_tot)
    return float(e), bool(mf.converged), mf

def relax_e(atoms, spin, charge=0, maxsteps=40):
    """Геом-оптимизация (geomeTRIC) на UKS-PBE; при провале — SP на затравке (флаг)."""
    mol = mk_mol(atoms, spin, charge)
    mf = mkmf(mol)
    try:
        mol_eq = optimize(mf, maxsteps=maxsteps)
        e, conv, _ = scf_e(mol_eq)
        return e, conv, True, [[mol_eq.atom_symbol(i),
                                tuple(float(x) for x in mol_eq.atom_coord(i, unit="Angstrom"))]
                               for i in range(mol_eq.natm)]
    except Exception as exc:
        say(f"      relax FAILED ({type(exc).__name__}: {str(exc)[:60]}) -> SP на затравке")
        e, conv, _ = scf_e(mol)
        return e, conv, False, atoms

def lowest_over_spins(tag, atoms_fn, spins, store, save):
    """Перебор мультиплетности; берём низшую по E релаксированную точку."""
    best = None
    for s in spins:
        skey = f"2S={s}"
        rec = store.setdefault(tag, {}).setdefault(skey, {})
        if rec.get("e") is not None:
            say(f"    [{tag}] {skey}: есть E={rec['e']:.5f} (рестарт)")
        else:
            t0 = time.time()
            try:
                e, conv, relaxed, geo = relax_e(atoms_fn(), s)
                rec.update(e=round(e, 6), converged=conv, relaxed=relaxed,
                           xyz=geo, t_s=round(time.time() - t0, 1))
                say(f"    [{tag}] {skey}: E={e:.5f} conv={conv} relaxed={relaxed} "
                    f"({rec['t_s']}s)")
            except Exception as exc:
                rec.update(e=None, converged=False,
                           error=f"{type(exc).__name__}: {str(exc)[:80]}")
                say(f"    [{tag}] {skey}: FAILED {exc}")
            save()
        if rec.get("e") is not None and rec.get("converged") and \
           (best is None or rec["e"] < best[1]):
            best = (skey, rec["e"])
    if best:
        store[tag]["lowest"] = {"spin": best[0], "e": best[1]}
    return best

def references(res, save):
    """Молекулы-эталоны N₂, H₂, NH₃ (релакс, синглеты) — метал-независимо."""
    refs = res.setdefault("refs", {})
    specs = {
        "N2": (("N", (0, 0, 0)), ("N", (0, 0, 1.10))),
        "H2": (("H", (0, 0, 0)), ("H", (0, 0, 0.74))),
        "NH3": (("N", (0, 0, 0.0)), ("H", (0.94, 0, -0.33)),
                ("H", (-0.47, 0.81, -0.33)), ("H", (-0.47, -0.81, -0.33))),
    }
    for name, atoms in specs.items():
        if refs.get(name, {}).get("e") is not None:
            say(f"  ref {name}: E={refs[name]['e']:.6f} (рестарт)"); continue
        t0 = time.time()
        a = [[el, tuple(map(float, xyz))] for el, xyz in atoms]
        e, conv, relaxed, geo = relax_e(a, 0)
        refs[name] = {"e": round(e, 6), "converged": conv, "relaxed": relaxed,
                      "t_s": round(time.time() - t0, 1)}
        say(f"  ref {name}: E={e:.6f} conv={conv} relaxed={relaxed} ({refs[name]['t_s']}s)")
        save()
    return refs

def descriptors(res):
    r = res.get("refs", {})
    if not all(r.get(k, {}).get("e") is not None for k in ("N2", "H2", "NH3")):
        return
    e_n2, e_h2, e_nh3 = r["N2"]["e"], r["H2"]["e"], r["NH3"]["e"]
    C = (e_nh3 - 0.5 * e_n2 - 1.5 * e_h2) * EV        # энергия образования NH₃, эВ
    res["C_nh3_formation_eV"] = round(C, 4)
    m4 = res.get("M4", {}).get("lowest", {}).get("e")
    m4n = res.get("M4N", {}).get("lowest", {}).get("e")
    if m4 is None or m4n is None:
        return
    dE_nit = (m4n - m4 - 0.5 * e_n2) * EV
    dE_hyd = (m4 + e_nh3 - m4n - 1.5 * e_h2) * EV
    closure = dE_nit + dE_hyd - C
    res["descriptors"] = {
        "dE_nit_eV": round(dE_nit, 4),
        "dE_hyd_eV": round(dE_hyd, 4),
        "sum_eV": round(dE_nit + dE_hyd, 4),
        "C_eV": round(C, 4),
        "cycle_closure_eV": round(closure, 4),
        "closure_ok": bool(abs(closure) < 0.05),
        "balance_eV": round(max(abs(dE_nit - C / 2), abs(dE_hyd - C / 2)), 4),
        "spin_M4": res["M4"]["lowest"]["spin"],
        "spin_M4N": res["M4N"]["lowest"]["spin"],
        "note": ("balance = max отклонение шага от C/2 (термонейтральной точки); "
                 "меньше = мягче обе стадии = дешевле цикл. closure_ok — сумма "
                 "шагов сошлась с энергией образования NH₃ (валидность точки)."),
    }

def run_metal(metal, smoke=False):
    if metal not in METALS:
        say(f"unknown metal {metal}"); return 1
    cfg = METALS[metal]
    OUT = outpath(metal)
    res = json.load(open(OUT)) if os.path.exists(OUT) else {}
    res.setdefault("meta", {
        "metal": metal,
        "what": f"chemical-looping ammonia motif {metal}4 + mu-N: dE_nit (N2 "
                "fixation) & dE_hyd (NH3 release), volcano balance descriptor; "
                "UKS-PBE/def2-SVP+DF, geomeTRIC relax, spin scan (lowest taken)",
        "honesty": "M4 CLUSTER not lattice -> motif energetics, physics is the "
                   "TREND across metals not the absolute enthalpy; gas phase; "
                   "PBE; no ZPE/entropy/solvent; def2-SVP",
        "selftest": "cycle closure: dE_nit + dE_hyd must equal C (NH3 formation "
                    "energy from refs) within 0.05 eV",
    })
    save = lambda: atomic_save(res, OUT)
    save()
    say(f"=== {metal} === threads={lib.num_threads()}")
    references(res, save)
    sm4 = [cfg["spins_m4"][0]] if smoke else cfg["spins_m4"]
    sm4n = [cfg["spins_m4n"][0]] if smoke else cfg["spins_m4n"]
    say(f"  M4 spin-scan {sm4}")
    lowest_over_spins("M4", lambda: build_m4(metal, cfg["a_mm"]), sm4, res, save)
    say(f"  M4N spin-scan {sm4n}")
    lowest_over_spins("M4N", lambda: build_m4n(metal, cfg["a_mm"], cfg["a_mn"]),
                      sm4n, res, save)
    descriptors(res)
    d = res.get("descriptors", {})
    res["status"] = "smoke" if smoke else (
        "ok" if d.get("dE_nit_eV") is not None else "incomplete")
    save()
    if d:
        say(f"  {metal}: ΔE_nit={d['dE_nit_eV']:+.3f}  ΔE_hyd={d['dE_hyd_eV']:+.3f}  "
            f"C={d['C_eV']:+.3f}  closure={d['cycle_closure_eV']:+.4f} "
            f"(ok={d['closure_ok']})  balance={d['balance_eV']:.3f} эВ")
    say(f"status={res['status']} -> {OUT}")
    return 0

def main(args):
    if args.smoke:
        return run_metal("Li", smoke=True)
    metal = os.environ.get("METAL")
    metals = [metal] if metal else list(METALS.keys())
    rc = 0
    for m in metals:
        rc |= run_metal(m, smoke=False)
    return rc

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--smoke", action="store_true")
    a = ap.parse_args()
    sys.exit(main(a))
