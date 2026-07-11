#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""routes/nh3_nitride_echem.py — ЭЛЕКТРОХИМИЧЕСКИЙ дескриптор мобильного устройства
для аммиака «из воздуха и воды» (медиированная eNRR). Надстройка над нитридной
лестницей: превращает «химический H₂-looping» в дескриптор КОРОБКИ-НА-ПОЛЕ.

Вопрос устройства (воздух+вода+электрон → NH₃, без метана, без завода): при
КАКОМ катодном потенциале нитрид M₄N отдаёт NH₃ протонированием, и обгоняет ли
это протонирование конкурирующее выделение водорода (HER)? Единственная реально
работающая ambient-версия — Li-опосредованная (Stanford/DTU, ~100% FE): Li⁰
рвёт N₂ химически, а H связывается на Li так плохо, что протон идёт на нитрид,
а не в H₂. Эту маржу считаем из первых принципов.

Дескрипторы на металл M (CHE = computational hydrogen electrode, U vs RHE):
  протон-сопряжённая лестница выпуска  M₄N → M₄NH → M₄NH₂ → M₄NH₃ (→ +NH₃):
    ΔG_n = E(прод) − E(исх) − ½E(H₂)             (электронно, U=0)
    U_L  = −max_n(ΔG_n)/e   — предельный потенциал (всё под гору) → МЯГЧЕ=лучше
    ΔG_des = E(M₄)+E(NH₃) − E(M₄NH₃)             (десорбция NH₃, химическая)
  HER-маржа:
    ΔG_H* = E(M₄H) − E(M₄) − ½E(H₂)  — H* на МЕТАЛЛЕ; ≫0 ⇒ HER подавлен (Li).
    her_margin = ΔG_H* − ΔG_1(N-протонирование) — >0 ⇒ протон идёт на N, не в H₂.

Источник геометрий/энергий: routes/nh3_nitride_<metal>_results.json (низший спин
M₄ и M₄N, релакс. эталоны N₂/H₂/NH₃). Интермедиаты строим добавлением H к N
(или к металлу для H*), релаксируем UKS-PBE/def2-SVP+DF, скан спина ±1.

ЧЕСТНО: ТОЛЬКО электронные энергии — без ZPE/энтропии/сольватации (для
протонирования это ±0.2–0.4 эВ, значимо; помечаем как v1-термодинамику, не
кинетику); кластерный мотив, не электрод; «кто выигрывает HER» — термодинам.
прокси, не барьер; CHE U vs RHE. Проба отвечает: при каком потенциале
термодинамика позволяет выпуск NH₃ и голодает ли водород.

Запуск:  METAL=Li python3 -u routes/nh3_nitride_echem.py
Smoke:   --smoke (Li, только 1-й протон + H*, грубо)
Выход:   routes/nh3_nitride_echem_<metal>_results.json (atomic, incremental)
"""
import os, sys, json, time, argparse, tempfile
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
from pyscf import gto, dft, lib  # noqa: E402
from pyscf.geomopt.geometric_solver import optimize  # noqa: E402

EV = 27.211386245988
# стандартные потенциалы осаждения M^n+/M (В vs SHE) — лит., для контекста размена
E_DEPOSIT_V = {"Li": -3.04, "Ca": -2.87, "Mg": -2.37, "Mn": -1.18, "Fe": -0.44,
               "Na": -2.71, "K": -2.93, "Sr": -2.89, "Ba": -2.91, "Al": -1.66,
               "V": -1.18, "Cr": -0.91, "Mo": -0.20}

def say(m): print(f"[{time.strftime('%H:%M:%S')}] {m}", flush=True)

def src_path(metal):  return os.path.join(HERE, f"nh3_nitride_{metal.lower()}_results.json")
def out_path(metal):  return os.path.join(HERE, f"nh3_nitride_echem_{metal.lower()}_results.json")

def atomic_save(obj, path):
    fd, tmp = tempfile.mkstemp(dir=HERE, suffix=".tmp")
    with os.fdopen(fd, "w") as f:
        json.dump(obj, f, ensure_ascii=False, indent=1)
    os.replace(tmp, path)

def mkmf(mol):
    mf = dft.UKS(mol).density_fit()
    mf.xc = os.environ.get("XC","pbe"); mf.conv_tol = 1e-9; mf.max_cycle = 200; mf.level_shift = 0.2
    return mf

def scf_e(mol):
    mf = mkmf(mol); e = mf.kernel()
    if not mf.converged:
        mfn = mf.newton(); mfn.max_cycle = 80
        e = mfn.kernel(mf.mo_coeff, mf.mo_occ)
        mf.converged = bool(mfn.converged); e = float(mfn.e_tot)
    return float(e), bool(mf.converged)

def relax_e(atoms, spin, maxsteps=40):
    mol = gto.M(atom=atoms, basis="def2-svp", ecp="def2-svp", spin=spin,
                charge=0, verbose=0)
    mf = mkmf(mol)
    try:
        meq = optimize(mf, maxsteps=maxsteps)
        e, conv = scf_e(meq)
        geo = [[meq.atom_symbol(i),
                tuple(float(x) for x in meq.atom_coord(i, unit="Angstrom"))]
               for i in range(meq.natm)]
        return e, conv, True, geo
    except Exception as exc:
        say(f"      relax FAILED ({type(exc).__name__}: {str(exc)[:50]}) -> SP")
        e, conv = scf_e(mol)
        return e, conv, False, atoms

def lowest_spin(atoms, spin_center, label, store, save, dspin=1):
    """Релакс со сканом мультиплетности вокруг spin_center; низшую по E берём."""
    if store.get(label, {}).get("e") is not None:
        say(f"    {label}: есть E={store[label]['e']:.5f} (рестарт)")
        return store[label]
    best = None
    nat = len(atoms)
    nel_parity = None  # спин должен совпадать по чётности с числом электронов
    cands = sorted(set(max(0, spin_center + d) for d in range(-dspin, dspin + 1)))
    for s in cands:
        t0 = time.time()
        try:
            e, conv, relaxed, geo = relax_e(atoms, s)
            if best is None or (conv and e < best["e"]):
                best = {"e": round(e, 6), "spin": s, "converged": conv,
                        "relaxed": relaxed, "xyz": geo, "t_s": round(time.time()-t0, 1)}
        except Exception as exc:
            say(f"    {label} 2S={s}: FAILED {exc}")
    if best:
        store[label] = best
        say(f"    {label}: E={best['e']:.5f} 2S={best['spin']} conv={best['converged']} "
            f"relaxed={best['relaxed']} ({best['t_s']}s)")
    save()
    return best or {}

def add_H_to(atoms, target_idx, existing_H_dirs, r=1.02):
    """Добавить H на 1.02 Å от target атома, в сторону от центроида, разводя от уже
    добавленных H."""
    coords = np.array([a[1] for a in atoms], float)
    centroid = coords.mean(0)
    t = coords[target_idx]
    base = t - centroid
    base = base / (np.linalg.norm(base) + 1e-9)
    # ортогональный джиттер, чтобы 2-й/3-й H не легли на 1-й
    jit = np.array([0.4 * np.cos(1.9 * len(existing_H_dirs)),
                    0.4 * np.sin(1.9 * len(existing_H_dirs)), 0.2 * len(existing_H_dirs)])
    d = base + jit; d = d / np.linalg.norm(d)
    return atoms + [["H", tuple(t + r * d)]], existing_H_dirs + [d]

def run_metal(metal, smoke=False):
    SRC, OUT = src_path(metal), out_path(metal)
    if not os.path.exists(SRC):
        say(f"нет {SRC} — сперва нужна нитридная лестница по {metal}"); return 2
    src = json.load(open(SRC))
    refs = src.get("refs", {})
    m4 = src.get("M4", {}).get("lowest", {})
    m4n = src.get("M4N", {}).get("lowest", {})
    if not (refs.get("H2", {}).get("e") and refs.get("NH3", {}).get("e")
            and m4.get("e") is not None and m4n.get("e") is not None):
        say(f"{metal}: неполная нитридная лестница (нет M4/M4N/эталонов)"); return 2
    # геометрии низшего спина
    m4_xyz = src["M4"][m4["spin"]]["xyz"]
    m4n_xyz = src["M4N"][m4n["spin"]]["xyz"]
    s_m4 = int(m4["spin"].split("=")[1]); s_m4n = int(m4n["spin"].split("=")[1])
    e_h2, e_nh3 = refs["H2"]["e"], refs["NH3"]["e"]

    res = json.load(open(OUT)) if os.path.exists(OUT) else {}
    res.setdefault("meta", {
        "metal": metal,
        "what": f"CHE proton ladder M4N->M4NH->M4NH2->M4NH3 (+NH3 desorb) and H* "
                f"HER descriptor on {metal}4 motif; U_L release potential & HER "
                "margin; UKS-PBE/def2-SVP+DF, relax, spin scan +-1",
        "honesty": "ELECTRONIC energies only (no ZPE/entropy/solvation, ~0.2-0.4 "
                   "eV for protonation); cluster motif not electrode; HER 'win' is "
                   "a thermodynamic proxy not a barrier; CHE U vs RHE; v1",
        "E_deposit_V_vs_SHE": E_DEPOSIT_V.get(metal),
    })
    res.setdefault("species", {})
    sp = res["species"]
    save = lambda: atomic_save(res, OUT)
    save()
    say(f"=== echem {metal} === threads={lib.num_threads()}  (base spins M4 {s_m4}, M4N {s_m4n})")

    # база из нитридной лестницы (не пересчитываем)
    sp.setdefault("M4", {"e": m4["e"], "spin": f"2S={s_m4}"})
    sp.setdefault("M4N", {"e": m4n["e"], "spin": f"2S={s_m4n}"})

    # H* на металле (HER-дескриптор): добавляем H к вершине M₄
    atoms_h, _ = add_H_to([[el, tuple(map(float, xyz))] for el, xyz in m4_xyz], 0, [])
    lowest_spin(atoms_h, max(0, s_m4 - 1) if s_m4 else 1, "M4H", sp, save)

    # протон-лестница на нитриде: N — последний атом M4N
    ladder = ["M4N"]
    cur_atoms = [[el, tuple(map(float, xyz))] for el, xyz in m4n_xyz]
    n_idx = len(cur_atoms) - 1
    dirs = []
    steps = 1 if smoke else 3
    cur_spin = s_m4n
    for n in range(1, steps + 1):
        cur_atoms, dirs = add_H_to(cur_atoms, n_idx, dirs)
        label = f"M4NH{n if n > 1 else ''}"
        cur_spin = max(0, cur_spin - 1)     # каждый H (e⁻) обычно спаривает
        best = lowest_spin(cur_atoms, cur_spin, label, sp, save)
        if best.get("xyz"):
            cur_atoms = [[el, tuple(map(float, xyz))] for el, xyz in best["xyz"]]
            cur_spin = int(best["spin"]) if isinstance(best.get("spin"), int) else cur_spin
        ladder.append(label)

    if not smoke:
        # CHE-лестница ΔG_n и U_L
        seq = [l for l in ladder if sp.get(l, {}).get("e") is not None]
        dG = []
        for a, b in zip(seq[:-1], seq[1:]):
            dG.append({"step": f"{a}->{b}",
                       "dG_eV": round((sp[b]["e"] - sp[a]["e"] - 0.5 * e_h2) * EV, 4)})
        u_l = round(-max(s["dG_eV"] for s in dG) / 1.0, 4) if dG else None  # эВ→В (e=1)
        m4nh3 = sp.get(f"M4NH{steps}" if steps > 1 else "M4NH", {}).get("e")
        dG_des = round((sp["M4"]["e"] + e_nh3 - m4nh3) * EV, 4) if m4nh3 is not None else None
        dG_h = (round((sp["M4H"]["e"] - sp["M4"]["e"] - 0.5 * e_h2) * EV, 4)
                if sp.get("M4H", {}).get("e") is not None else None)
        her_margin = (round(dG_h - dG[0]["dG_eV"], 4) if dG and dG_h is not None else None)
        res["descriptors"] = {
            "che_ladder": dG,
            "U_L_V_vs_RHE": u_l,
            "dG_desorb_NH3_eV": dG_des,
            "dG_Hstar_eV": dG_h,
            "her_margin_eV": her_margin,
            "E_deposit_V_vs_SHE": E_DEPOSIT_V.get(metal),
            "note": ("U_L: потенциал, при котором ВСЕ протонные шаги под гору "
                     "(мягче=положительнее=дешевле). her_margin>0 ⇒ протон "
                     "термодинамически идёт на N, а не в H₂ (эффект Li)."),
        }
        res["status"] = "ok" if (u_l is not None and dG_h is not None) else "incomplete"
        save()
        d = res["descriptors"]
        say(f"  {metal}: U_L={d['U_L_V_vs_RHE']} В  ΔG_H*={d['dG_Hstar_eV']} эВ  "
            f"HER-маржа={d['her_margin_eV']} эВ  (осаждение {d['E_deposit_V_vs_SHE']} В)")
    else:
        res["status"] = "smoke"; save()
        say(f"SMOKE {metal}: M4H={sp.get('M4H',{}).get('e')} M4NH={sp.get('M4NH',{}).get('e')}")
    say(f"status={res['status']} -> {OUT}")
    return 0

def main(args):
    if args.smoke:
        return run_metal("Li", smoke=True)
    metal = os.environ.get("METAL")
    metals = [metal] if metal else ["Li", "Mg", "Ca", "Mn", "Fe"]
    rc = 0
    for m in metals:
        rc |= run_metal(m, smoke=False)
    return rc

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--smoke", action="store_true")
    a = ap.parse_args()
    sys.exit(main(a))
