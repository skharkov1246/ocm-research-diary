#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""routes/spin_field_oxidation.py — PAYOFF спинового транзистора: снижает ли
ЭЛЕКТРИЧЕСКОЕ поле барьер НАСТОЯЩЕЙ оксидации через спин-свитч O₂?

Этап 25 подтвердил: поле флипает основной спин активации O₂ на дешёвом Mn (в двух
функционалах). Вопрос-плата: превращается ли это в управление РЕАКЦИЕЙ? Считаем
барьер отрыва H от метана дистальным кислородом [M(NH₃)₄(O₂)] (первый и лимитирую-
щий шаг любой C–H оксидации, «rebound»-механизм) при полях −0.51 / 0 / +0.51 В/Å,
которые дают РАЗНЫЙ основной спин. Если барьер сильно зависит от поля — поле есть
электрическая ручка ВКЛ/ВЫКЛ мягкой оксидации метана, а на электроде поле бесплатно.

Метод: [M(NH₃)₄(O₂)] end-on + CH₄ над дистальным O (одна C–H смотрит на O).
Distinguished-coordinate: замораживаем d(O_dist–H_abs), скан 1.6→1.0 Å (H уходит с
C на O), релаксируем остальное (geomeTRIC, warm-chain); барьер = max(E)−E(старт).
Finite-field вдоль оси M–O–O (подмена get_hcore). Скан спина (низший берём).

ОГРАНИЧЕНИЯ: distinguished coordinate = верхняя оценка барьера (не истинное седло);
спиновые зазоры/барьеры функционал-зависимы (PBE0-проверка отдельно, env XC);
кластер, газовая фаза, замороженный каркас (честны СДВИГИ по полю); v1. Проба
отвечает: двигает ли поле барьер C–H оксидации через спиновую поверхность.

Запуск: METAL=Mn OMP_NUM_THREADS=8 python3 -u routes/spin_field_oxidation.py
Smoke:  --smoke (F=0, один спин, d=1.3 — проверка пайплайна)
Выход:  routes/spin_field_oxidation_<metal>[_xc]_results.json
"""
import os, sys, json, time, argparse, tempfile
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
from pyscf import gto, dft, lib  # noqa: E402
from pyscf.geomopt.geometric_solver import optimize  # noqa: E402
from pyscf.data.elements import charge as Z  # noqa: E402

METAL = os.environ.get("METAL", "Fe")
XC = os.environ.get("XC", "pbe")
CHARGE = int(os.environ.get("CHARGE", "1"))      # [MO]+ канонический катион (Shaik HAT)
RELAX = int(os.environ.get("RELAX", "0"))        # 1 = релаксированный барьер (честнее, дороже)
_tag = ("_" + METAL.lower()) + ("" if XC == "pbe" else "_" + XC) + ("_relax" if RELAX else "")
OUT = os.path.join(HERE, f"spin_field_oxidation{_tag}_results.json")
EV = 27.211386245988
VA = 51.42206747
FIELDS = [0.0, 0.010, -0.010]                    # 0, ±0.514 В/Å (разный осн. спин)
DGRID = [1.55, 1.40, 1.28, 1.18, 1.08, 1.00]     # d(O–H_abs), Å
SPINS = [0, 1, 2, 3, 4, 5, 6, 7]

def say(m): print(f"[{time.strftime('%H:%M:%S')}] {m}", flush=True)

def atomic_save(o):
    fd, tmp = tempfile.mkstemp(dir=HERE, suffix=".tmp")
    with os.fdopen(fd, "w") as f: json.dump(o, f, ensure_ascii=False, indent=1)
    os.replace(tmp, OUT)

def build():
    """[MO]+ + CH4: канонический оксо-HAT (Shaik). M@0, O@+z, CH4 над O, одна C–H к O."""
    m_o = 1.63                                    # M=O
    Oz = m_o
    Cz = Oz + 2.55                                # C над оксо-кислородом
    atoms = [[METAL, (0.0, 0.0, 0.0)], ["O", (0.0, 0.0, Oz)],
             ["C", (0.0, 0.0, Cz)], ["H", (0.0, 0.0, Cz - 1.09)]]   # abstractable H
    for a in (0, 2 * np.pi / 3, 4 * np.pi / 3):
        atoms.append(["H", (1.03 * np.cos(a), 1.03 * np.sin(a), Cz + 0.36)])
    return atoms                                  # M, O, C, H_abs, 3×H_methyl = 7

def idx():
    return 1, 3, 7                                # O=1, H_abs=3, nat=7

def mkmf(mol, f_au):
    mf = dft.UKS(mol).density_fit()
    mf.xc = XC; mf.conv_tol = 1e-8; mf.max_cycle = 200; mf.level_shift = 0.3
    if abs(f_au) > 0:
        mol.set_common_orig((0.0, 0.0, 0.0)); ao = mol.intor_symmetric("int1e_r", comp=3)
        h0 = mf.get_hcore(); E = np.array([0.0, 0.0, f_au])
        mf.get_hcore = lambda *a: h0 + np.einsum("x,xij->ij", E, ao)
    return mf

def scf(atoms, spin, f_au, dm0=None):
    mol = gto.M(atom=atoms, basis="def2-svp", spin=spin, charge=CHARGE, verbose=0)
    mf = mkmf(mol, f_au)
    e = mf.kernel(dm0=dm0) if dm0 is not None else mf.kernel()
    if not mf.converged:
        mfn = mf.newton(); mfn.max_cycle = 80
        e = mfn.kernel(mf.mo_coeff, mf.mo_occ); mf.converged = bool(mfn.converged); e = float(mfn.e_tot)
    return float(e), bool(mf.converged), mf.make_rdm1()

def crelax(atoms, spin, f_au, i, j, d, maxsteps=30):
    a = [[el, tuple(map(float, xyz))] for el, xyz in atoms]
    p_i = np.array(a[i][1]); p_j = np.array(a[j][1]); v = p_j - p_i
    v = v / (np.linalg.norm(v) + 1e-9); a[j] = [a[j][0], tuple(p_i + v * d)]
    mol = gto.M(atom=a, basis="def2-svp", spin=spin, charge=CHARGE, verbose=0)
    mf = mkmf(mol, f_au)
    fd, cf = tempfile.mkstemp(dir=HERE, suffix=".txt")
    with os.fdopen(fd, "w") as fh: fh.write(f"$freeze\ndistance {i+1} {j+1}\n")
    try:
        meq = optimize(mf, constraints=cf, maxsteps=maxsteps)
        geo = [[meq.atom_symbol(k), tuple(float(x) for x in meq.atom_coord(k, unit="Angstrom"))]
               for k in range(meq.natm)]
        e, conv, _ = scf(geo, spin, f_au); return e, conv, geo
    except Exception as exc:
        say(f"      relax d={d} FAILED ({type(exc).__name__}: {str(exc)[:45]})")
        try: e, conv, _ = scf(a, spin, f_au); return e, conv, a
        except Exception: return None, False, None
    finally:
        try: os.remove(cf)
        except Exception: pass

def place_H(base, o_i, h_j, d):
    """Жёсткая посадка переносимого H на расстоянии d от O вдоль оси O→C(H)."""
    a = [[el, tuple(map(float, xyz))] for el, xyz in base]
    p_o = np.array(a[o_i][1]); p_h = np.array(a[h_j][1])
    v = p_h - p_o; v = v / (np.linalg.norm(v) + 1e-9)
    a[h_j] = [a[h_j][0], tuple(p_o + v * d)]
    return a

def barrier_at_field(f_au, spin, o_i, h_j, store, save):
    """Жёсткий скан d(O–H) с warm-chain плотности (без релаксации): честны СДВИГИ."""
    blk = store.setdefault(f"{f_au:+.3f}", {}).setdefault(f"2S={spin}", {"points": {}})
    base = build(); dm = None; cur = base
    for d in DGRID:
        dk = f"{d:.2f}"
        if blk["points"].get(dk, {}).get("e") is not None:
            continue
        t0 = time.time()
        try:
            if RELAX:
                e, conv, geo = crelax(cur, spin, f_au, o_i, h_j, d)
                if geo is not None: cur = geo
            else:
                e, conv, dm = scf(place_H(base, o_i, h_j, d), spin, f_au, dm0=dm)
        except Exception as exc:
            e, conv = None, False
            say(f"      pt d={d} FAILED ({type(exc).__name__})")
        blk["points"][dk] = {"e": round(e, 6) if e is not None else None, "converged": conv}
        say(f"    F={f_au:+.3f} 2S={spin} d={dk}: E={e if e is None else round(e,5)} conv={conv} ({round(time.time()-t0,1)}s)")
        save()
    pts = sorted((float(k), v["e"]) for k, v in blk["points"].items() if v.get("converged") and v.get("e") is not None)
    if len(pts) >= 3:
        E = np.array([p[1] for p in pts]); blk["barrier_eV"] = round(float((E.max() - E[0]) * EV), 4)
    return blk.get("barrier_eV")

def main(a):
    o_i, h_j, nat = idx()
    nel = Z(METAL) + 8 + 6 + 4 - CHARGE           # [MO]+ + CH4
    spins = [s for s in SPINS if s % 2 == nel % 2] or [nel % 2]
    say(f"spin-field OXIDATION [{METAL}O]+ + CH4 (oxo HAT) — xc={XC} nat={nat} O_dist={o_i} H_abs={h_j} spins={spins} threads={lib.num_threads()}")
    res = json.load(open(OUT)) if os.path.exists(OUT) else {}
    res.setdefault("meta", {
        "what": f"C-H abstraction barrier of CH4 by [{METAL}O]+ oxo vs axial field "
                "(0,+-0.51 V/A); RIGID d(O-H) scan, density warm-chain (field-shifts "
                "honest, absolute barrier overestimated); does field move HAT barrier/spin?",
        "metal": METAL, "xc": XC, "fields_au": FIELDS,
        "limits": "distinguished-coord = upper-bound barrier (not true saddle); "
                  "functional-sensitive spins/barriers; cluster, gas phase, frozen "
                  "scaffold (field-shifts honest); v1",
    })
    res.setdefault("scan", {})
    fields = [0.0] if a.smoke else FIELDS
    sp = [spins[0]] if a.smoke else spins
    rows = []
    for f in fields:
        best = None
        for s in sp:
            b = barrier_at_field(f, s, o_i, h_j, res["scan"], lambda: atomic_save(res))
            if b is not None and (best is None or b < best[1]): best = (s, b)
        if best: rows.append({"field_V_per_A": round(f * VA, 3), "spin_2S": best[0], "barrier_eV": best[1]})
    res["summary_rows"] = rows
    if len(rows) >= 2:
        b0 = next((r["barrier_eV"] for r in rows if r["field_V_per_A"] == 0.0), None)
        res["barrier_zero_field_eV"] = b0
        res["field_shifts_barrier"] = {r["field_V_per_A"]: r["barrier_eV"] for r in rows}
        if len(rows) >= 3:
            F = np.array([r["field_V_per_A"] for r in rows]); B = np.array([r["barrier_eV"] for r in rows])
            res["dBarrier_dF_eV_per_VA"] = round(float(np.polyfit(F, B, 1)[0]), 4)
    res["status"] = "smoke" if a.smoke else ("ok" if len(rows) == len(FIELDS) else "incomplete")
    atomic_save(res)
    say(f"[{METAL}] status={res['status']} rows={rows} dBarrier/dF={res.get('dBarrier_dF_eV_per_VA')} -> {OUT}")
    return 0

if __name__ == "__main__":
    ap = argparse.ArgumentParser(); ap.add_argument("--smoke", action="store_true")
    sys.exit(main(ap.parse_args()))
