#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""routes/nh3_nitride_alloy.py — СПЛАВНЫЙ мотив петли аммиака: можно ли украсть
ОБА свойства сразу. Прямое продолжение Этапа 21.

Находка Этапа 21: ни один чистый металл не Парето-оптимален.
  Li  — единственный морит HER голодом (HER-маржа +2.68 эВ), но термодинамически
        перекошен (balance 1.64: жадно ловит N₂, реклюктантно отдаёт NH₃).
  Mg  — самая мягкая термодинамика (balance 0.67, лёгкий выпуск NH₃), но HER
        выигрывает (маржа −1.40).
Гипотеза: сплавной мотив LiₓMg₄₋ₓ+μ-N наследует H-голодание от Li-вершин (сайт
адсорбции H на Li остаётся плохим) И мягкий выпуск от Mg-химии нитрида. Если у
какого-то состава balance ниже Li ПРИ HER-марже >0 — это точка «почти
бесплатной» мобильной петли, к которой Этап 21 показал направление.

Считаем те же дескрипторы, что Этапы 20-21, но на смешанном мотиве:
  base:  ΔE_nit, ΔE_hyd, balance, замыкание цикла (самотест, тот же C эталонов);
  echem: CHE протон-лестница M₄N→M₄NH₃ → U_L; H* НА Li-ВЕРШИНЕ → ΔG_H*,
         HER-маржа = ΔG_H* − ΔG(1-й N-протон).
Составы: Li3Mg, Li2Mg2, LiMg3 (env COMPO выбирает один; по умолчанию все).

ЧЕСТНО: тот же кластерный уровень, что Этап 21 (мотив, не решётка; тренд, не
абсолют); H* сажаем на Li-вершину (гипотеза «Li-сайт морит HER» — проверяем её
же); только электронные энергии, UKS-PBE/def2-SVP+DF, скан спина; газовая фаза.

Запуск:  COMPO=Li3Mg python3 -u routes/nh3_nitride_alloy.py
Smoke:   --smoke (Li2Mg2, base + 1-й протон + H*, грубо)
Выход:   routes/nh3_nitride_alloy_<compo>_results.json
"""
import os, sys, json, time, argparse, tempfile
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
from pyscf import gto, dft, lib  # noqa: E402
from pyscf.data.elements import charge as Z  # noqa: E402
from pyscf.geomopt.geometric_solver import optimize  # noqa: E402

EV = 27.211386245988
COMPOS = {
    "Li3Mg": ["Li", "Li", "Li", "Mg"],
    "Li2Mg2": ["Li", "Li", "Mg", "Mg"],
    "LiMg3": ["Li", "Mg", "Mg", "Mg"],
}
# грубые ребра M-M (Å) для средней геометрии сплава; релаксируется
A_ELEM = {"Li": 3.0, "Na": 3.7, "K": 4.6, "Mg": 3.2, "Ca": 3.9, "Sr": 4.3,
          "Ba": 4.3, "Al": 2.8}

import re  # noqa: E402
def parse_compo(s):
    """'Li3Al' / 'Li2Ca2' / 'Mg3Al1' -> список из 4 символов вершин."""
    toks = re.findall(r"([A-Z][a-z]?)(\d*)", s)
    els = []
    for el, n in toks:
        if el: els += [el] * (int(n) if n else 1)
    if len(els) != 4:
        raise ValueError(f"состав {s} даёт {len(els)} вершин, нужно 4")
    return els

def a_mm(els):
    return sum(A_ELEM.get(e, 3.1) for e in els) / len(els)

def say(m): print(f"[{time.strftime('%H:%M:%S')}] {m}", flush=True)
TAG = os.environ.get("RUN_TAG", "")
def outpath(c): return os.path.join(HERE, f"nh3_nitride_alloy_{c.lower()}{TAG}_results.json")

def atomic_save(obj, path):
    fd, tmp = tempfile.mkstemp(dir=HERE, suffix=".tmp")
    with os.fdopen(fd, "w") as f:
        json.dump(obj, f, ensure_ascii=False, indent=1)
    os.replace(tmp, path)

def tetra(a):
    v = np.array([[1, 1, 1], [1, -1, -1], [-1, 1, -1], [-1, -1, 1]], float)
    return v * a / np.sqrt(8.0)

def parity_spins(atoms, charge, want):
    nel = sum(Z(el) for el, _ in atoms) - charge
    return [s for s in want if s % 2 == nel % 2] or [nel % 2]

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
    mol = gto.M(atom=atoms, basis="def2-svp", spin=spin, charge=0, verbose=0)
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

def lowest(atoms, want, label, store, save):
    if store.get(label, {}).get("e") is not None:
        say(f"    {label}: есть (рестарт)"); return store[label]
    best = None
    for s in parity_spins(atoms, 0, want):
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

def refs(store, save):
    specs = {"N2": (("N", (0, 0, 0)), ("N", (0, 0, 1.10))),
             "H2": (("H", (0, 0, 0)), ("H", (0, 0, 0.74))),
             "NH3": (("N", (0, 0, 0.0)), ("H", (0.94, 0, -0.33)),
                     ("H", (-0.47, 0.81, -0.33)), ("H", (-0.47, -0.81, -0.33)))}
    r = store.setdefault("refs", {})
    for name, atoms in specs.items():
        if r.get(name, {}).get("e") is not None: continue
        a = [[el, tuple(map(float, xyz))] for el, xyz in atoms]
        e, conv, relaxed, _ = relax_e(a, 0)
        r[name] = {"e": round(e, 6), "converged": conv}
        say(f"  ref {name}: E={e:.6f} conv={conv}"); save()
    return r

def add_H(atoms, target_idx, existing, r=1.02):
    coords = np.array([a[1] for a in atoms], float)
    t = coords[target_idx]; base = t - coords.mean(0)
    base = base / (np.linalg.norm(base) + 1e-9)
    jit = np.array([0.4*np.cos(1.9*len(existing)), 0.4*np.sin(1.9*len(existing)), 0.2*len(existing)])
    d = base + jit; d = d/np.linalg.norm(d)
    return atoms + [["H", tuple(t + r*d)]], existing + [d]

def run(compo, smoke=False):
    els = COMPOS[compo] if compo in COMPOS else parse_compo(compo)
    OUT = outpath(compo)
    res = json.load(open(OUT)) if os.path.exists(OUT) else {}
    res.setdefault("meta", {
        "compo": compo, "elements": els,
        "what": f"alloy motif {compo}4 + mu-N: base (dE_nit/dE_hyd/balance, cycle "
                "closure selftest) + echem (CHE U_L, H* on a Li vertex -> HER "
                "margin); UKS-PBE/def2-SVP+DF, relax, parity spin scan",
        "hypothesis": "Li vertices keep HER starved (H* bad on Li) while Mg "
                      "chemistry softens NH3 release -> lower balance at HER-margin>0",
        "honesty": "cluster motif not lattice (trend not absolute); H* placed on a "
                   "Li vertex (tests the very hypothesis); electronic energies only; "
                   "PBE; gas phase; def2-SVP",
        "ref_pure_Li": {"balance": 1.64, "her_margin": 2.68},
        "ref_pure_Mg": {"balance": 0.67, "her_margin": -1.40},
    })
    sp = res.setdefault("species", {})
    save = lambda: atomic_save(res, OUT)
    save()
    say(f"=== alloy {compo} {els} === threads={lib.num_threads()}")
    r = refs(res, save)
    e_n2, e_h2, e_nh3 = r["N2"]["e"], r["H2"]["e"], r["NH3"]["e"]

    coords = tetra(a_mm(els))
    m4 = [[els[i], tuple(coords[i])] for i in range(4)]
    m4n = m4 + [["N", (0.0, 0.0, 0.0)]]
    # широкие списки; parity_spins сам отберёт нужную чётность по числу электронов
    want = [0, 1, 2, 3, 4] if not smoke else [0, 1, 2]
    lowest(m4, want, "M4", sp, save)
    lowest(m4n, want, "M4N", sp, save)

    # H* на самой электроположительной вершине (Li если есть, иначе 1-й металл)
    li_idx = els.index("Li") if "Li" in els else 0
    m4_xyz = [[el, tuple(map(float, xyz))] for el, xyz in sp["M4"]["xyz"]]
    ah, _ = add_H(m4_xyz, li_idx, [])
    lowest(ah, [0, 1, 2, 3], "M4H", sp, save)

    # протон-лестница на нитриде (N последний)
    cur = [[el, tuple(map(float, xyz))] for el, xyz in sp["M4N"]["xyz"]]
    n_idx = len(cur) - 1; dirs = []
    steps = 1 if smoke else 3
    ladder = ["M4N"]
    for n in range(1, steps + 1):
        cur, dirs = add_H(cur, n_idx, dirs)
        label = f"M4NH{n if n > 1 else ''}"
        best = lowest(cur, [0, 1, 2, 3, 4], label, sp, save)
        if best.get("xyz"):
            cur = [[el, tuple(map(float, xyz))] for el, xyz in best["xyz"]]
        ladder.append(label)

    # дескрипторы
    m4e, m4ne = sp["M4"]["e"], sp["M4N"]["e"]
    C = (e_nh3 - 0.5*e_n2 - 1.5*e_h2) * EV
    dE_nit = (m4ne - m4e - 0.5*e_n2) * EV
    dE_hyd = (m4e + e_nh3 - m4ne - 1.5*e_h2) * EV
    d = {"dE_nit_eV": round(dE_nit, 4), "dE_hyd_eV": round(dE_hyd, 4),
         "C_eV": round(C, 4), "cycle_closure_eV": round(dE_nit+dE_hyd-C, 4),
         "closure_ok": bool(abs(dE_nit+dE_hyd-C) < 0.05),
         "balance_eV": round(max(abs(dE_nit-C/2), abs(dE_hyd-C/2)), 4)}
    if not smoke:
        seq = [l for l in ladder if sp.get(l, {}).get("e") is not None]
        dG = [{"step": f"{a}->{b}",
               "dG_eV": round((sp[b]["e"]-sp[a]["e"]-0.5*e_h2)*EV, 4)}
              for a, b in zip(seq[:-1], seq[1:])]
        d["U_L_V_vs_RHE"] = round(-max(x["dG_eV"] for x in dG), 4) if dG else None
        dG_h = (round((sp["M4H"]["e"]-sp["M4"]["e"]-0.5*e_h2)*EV, 4)
                if sp.get("M4H", {}).get("e") is not None else None)
        d["dG_Hstar_eV"] = dG_h
        d["her_margin_eV"] = round(dG_h - dG[0]["dG_eV"], 4) if (dG and dG_h is not None) else None
        d["che_ladder"] = dG
    res["descriptors"] = d
    res["status"] = "smoke" if smoke else ("ok" if d.get("her_margin_eV") is not None else "incomplete")
    save()
    say(f"  {compo}: balance={d['balance_eV']} эВ (Li 1.64 / Mg 0.67)  "
        f"HER-маржа={d.get('her_margin_eV')} эВ (Li +2.68 / Mg −1.40)  "
        f"closure_ok={d['closure_ok']}")
    say(f"status={res['status']} -> {OUT}")
    return 0

def main(a):
    if a.smoke: return run("Li2Mg2", smoke=True)
    c = os.environ.get("COMPO")
    for cc in ([c] if c else list(COMPOS)):
        run(cc, smoke=False)
    return 0

if __name__ == "__main__":
    ap = argparse.ArgumentParser(); ap.add_argument("--smoke", action="store_true")
    sys.exit(main(ap.parse_args()))
