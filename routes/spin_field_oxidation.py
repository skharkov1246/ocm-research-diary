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
Distinguished-coordinate: замораживаем d(O_dist–H_abs), скан 2.55→1.0 Å (H стартует
НА углероде — истинный реагент — и уходит на O). Протокол v2 (после адверсариального
ревью 25.07): геометрии релаксируются ПРИ F=0 (geomeTRIC, warm-chain, кадр M@0,
ось M–O = z), поле прикладывается только в single-point — чистый электронный
Stark-эффект, энергия с ядерным членом −F·ΣZ_A·z_A; forward-барьер = max(E)−E(max d).
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
RELAX = int(os.environ.get("RELAX", "0"))        # 1 = релаксированный барьер (дороже, ближе к седлу)
SP_XC = os.environ.get("SP_XC")                  # гибрид-контроль: геометрия по XC (надёжный GGA),
                                                 # энергия барьера — single-point по SP_XC (напр. pbe0//pbe)
SUB = os.environ.get("SUBSTRATE", "ch4").lower() # ch4 | ch3oh (селективность: барьер перекисления продукта)
SOLV_EPS = float(os.environ["SOLV_EPS"]) if os.environ.get("SOLV_EPS") else None
                                                 # ddCOSMO single-point на релакс-геометрии (экранировка средой)
RUN_TAG = os.environ.get("RUN_TAG", "")          # per-instance тег: параллельные поля не затирают друг друга
_eff = SP_XC if SP_XC else XC                    # эффективный функционал для имени/меты
_tag = (("_" + METAL.lower()) + ("" if SUB == "ch4" else "_" + SUB)
        + ("" if _eff == "pbe" else "_" + _eff)
        + (f"_eps{int(SOLV_EPS)}" if SOLV_EPS else "")
        + ("_relax" if RELAX else "") + (("_" + RUN_TAG) if RUN_TAG else ""))
OUT = os.path.join(HERE, f"spin_field_oxidation{_tag}_results.json")
EV = 27.211386245988
VA = 51.42206747
if os.environ.get("FIELDS_AU"):                  # явный список полей (а.е.) — один инстанс = одно поле
    FIELDS = [float(x) for x in os.environ["FIELDS_AU"].split(",")]
elif os.environ.get("FINE") == "1":
    FIELDS = [-0.012,-0.008,-0.004,0.0,0.004,0.008,0.012]   # FINE=1 — кривая барьер/поле
else:
    FIELDS = [0.0, 0.010, -0.010]
DGRID = ([2.55, 2.15, 1.80, 1.50, 1.20, 1.00] if os.environ.get("COARSE")=="1"
         else [2.55, 2.30, 2.05, 1.85, 1.65, 1.45, 1.25, 1.05])  # d(O–H_abs), Å: реагент H-на-C → продукт O–H
SPINS = [int(x) for x in os.environ.get("SPINS","0,1,2,3,4,5,6,7").split(",")]

def say(m): print(f"[{time.strftime('%H:%M:%S')}] {m}", flush=True)

def atomic_save(o):
    fd, tmp = tempfile.mkstemp(dir=HERE, suffix=".tmp")
    with os.fdopen(fd, "w") as f: json.dump(o, f, ensure_ascii=False, indent=1)
    os.replace(tmp, OUT)

def build():
    """[MO]+ + субстрат: канонический оксо-HAT (Shaik). M@0, O@+z, C над O, одна C–H к O.
    ch4: метан. ch3oh: метанол, отрываем C–H (доминантный HAT-канал, BDE C–H < O–H) —
    путь перекисления продукта в формальдегид: мера селективности."""
    m_o = 1.63                                    # M=O
    Oz = m_o
    Cz = Oz + 3.69                                # C выше: d(O–H_abs) стартует ~2.6 Å (H ещё НА C, истинный реагент)
    atoms = [[METAL, (0.0, 0.0, 0.0)], ["O", (0.0, 0.0, Oz)],
             ["C", (0.0, 0.0, Cz)], ["H", (0.0, 0.0, Cz - 1.09)]]   # abstractable H, коллинеарно O–H–C
    angs = (0, 2 * np.pi / 3, 4 * np.pi / 3)
    if SUB == "ch4":
        for a in angs:
            atoms.append(["H", (1.03 * np.cos(a), 1.03 * np.sin(a), Cz + 0.36)])
    else:                                         # ch3oh: 2 H + OH-группа вместо третьего H (C–O 1.43, O–H 0.96)
        for a in angs[:2]:
            atoms.append(["H", (1.03 * np.cos(a), 1.03 * np.sin(a), Cz + 0.36)])
        a = angs[2]
        o_alc = np.array([1.352 * np.cos(a), 1.352 * np.sin(a), Cz + 0.472])
        atoms.append(["O", tuple(o_alc)])
        cpos = np.array([0.0, 0.0, Cz])           # H_alc под 108.5° к связи O–C (не коллинеарно!)
        u = (cpos - o_alc); u /= np.linalg.norm(u)
        zax = np.array([0.0, 0.0, 1.0]); perp = zax - (zax @ u) * u; perp /= np.linalg.norm(perp)
        th = np.deg2rad(108.5)
        atoms.append(["H", tuple(o_alc + 0.96 * (np.cos(th) * u + np.sin(th) * perp))])
    return atoms                                  # ch4: 7 атомов; ch3oh: 8

def idx():
    return 1, 3, (7 if SUB == "ch4" else 8)       # O_dist=1, H_abs=3, nat

def mkmf(mol, f_au, xc=None, eps=None):
    """SP-фабрика. Поле — ТОЛЬКО для single-point (лямбда замыкает интегралы текущей
    геометрии — в оптимизатор такой mf отдавать нельзя: интегралы протухнут, а градиент
    поля не увидит; релаксация всегда при F=0 через relax-mf в crelax)."""
    mf = dft.UKS(mol).density_fit()
    mf.xc = xc or XC; mf.conv_tol = 1e-8; mf.max_cycle = 200; mf.level_shift = 0.3
    if eps:                                       # солватацию вешаем ДО подмены hcore: ddCOSMO-обёртка
        from pyscf import solvent                 # создаёт новый объект и потеряла бы instance-лямбду
        mf = solvent.ddCOSMO(mf); mf.with_solvent.eps = float(eps)
    if abs(f_au) > 0:
        mol.set_common_orig((0.0, 0.0, 0.0)); ao = mol.intor_symmetric("int1e_r", comp=3)
        h0 = mf.get_hcore(); E = np.array([0.0, 0.0, f_au])
        mf.get_hcore = lambda *a: h0 + np.einsum("x,xij->ij", E, ao)
    return mf

def scf(atoms, spin, f_au, dm0=None, xc=None, eps=None):
    mol = gto.M(atom=atoms, basis="def2-svp", spin=spin, charge=CHARGE, verbose=0)
    mf = mkmf(mol, f_au, xc, eps)
    e = mf.kernel(dm0=dm0) if dm0 is not None else mf.kernel()
    if not mf.converged:
        try:
            mfn = mf.newton(); mfn.max_cycle = 80
            e = mfn.kernel(mf.mo_coeff, mf.mo_occ); mf.converged = bool(mfn.converged); e = float(mfn.e_tot)
        except Exception:                         # newton не для всех обёрток (ddCOSMO) — оставляем как есть
            pass
    if abs(f_au) > 0:                             # ядерный член поля −F·Σ Z_A z_A (тот же origin (0,0,0));
        e = float(e) - f_au * float(np.dot(mol.atom_charges(), mol.atom_coords()[:, 2]))
    return float(e), bool(mf.converged), mf.make_rdm1()

def align(geo):
    """Единый кадр после релаксации: M в нуле, ось M→O_dist вдоль +z (ось поля).
    Обязательно: заряженный кластер + поле → энергия зависит от положения/ориентации."""
    P = np.array([x[1] for x in geo], float); P = P - P[0]
    v = P[1] / (np.linalg.norm(P[1]) + 1e-12); z = np.array([0.0, 0.0, 1.0])
    ax = np.cross(v, z); s = np.linalg.norm(ax); c = float(v @ z)
    if s > 1e-8:
        ax = ax / s; K = np.array([[0, -ax[2], ax[1]], [ax[2], 0, -ax[0]], [-ax[1], ax[0], 0]])
        th = np.arctan2(s, c); R = np.eye(3) + np.sin(th) * K + (1 - np.cos(th)) * (K @ K)
        P = P @ R.T
    elif c < 0:                                   # антипараллельно: разворот на π вокруг x
        P = P * np.array([1.0, -1.0, -1.0])
    return [[geo[k][0], tuple(float(x) for x in P[k])] for k in range(len(geo))]

def crelax(atoms, spin, i, j, d, maxsteps=int(os.environ.get("RELAX_STEPS", "30"))):
    """Constrained-релаксация ПРИ F=0 (поле — только single-point на итоговой геометрии).
    Возвращает (geo, relaxed): relaxed=False — geomeTRIC упал, geo = жёсткая посадка."""
    a = [[el, tuple(map(float, xyz))] for el, xyz in atoms]
    p_i = np.array(a[i][1]); p_j = np.array(a[j][1]); v = p_j - p_i
    v = v / (np.linalg.norm(v) + 1e-9); a[j] = [a[j][0], tuple(p_i + v * d)]
    mol = gto.M(atom=a, basis="def2-svp", spin=spin, charge=CHARGE, verbose=0)
    mf = dft.UKS(mol).density_fit()               # relax-mf: БЕЗ level_shift (ронял сканер в волне 1)
    mf.xc = XC; mf.conv_tol = 1e-7; mf.max_cycle = 300
    fd, cf = tempfile.mkstemp(dir=HERE, suffix=".txt")
    with os.fdopen(fd, "w") as fh: fh.write(f"$freeze\ndistance {i+1} {j+1}\n")
    try:
        meq = optimize(mf, constraints=cf, maxsteps=maxsteps)
        geo = [[meq.atom_symbol(k), tuple(float(x) for x in meq.atom_coord(k, unit="Angstrom"))]
               for k in range(meq.natm)]
        return align(geo), True
    except Exception as exc:
        say(f"      relax d={d} FAILED ({type(exc).__name__}: {str(exc)[:60]})")
        return align(a), False
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

def relax_chain(spin, o_i, h_j, res, save):
    """Релакс-цепочка F=0 по DGRID (warm-chain геометрий), кэш в res['geoms'].
    Геометрии ОБЩИЕ для всех полей — поле прикладывается только в SP."""
    g = res.setdefault("geoms", {}).setdefault(f"2S={spin}", {})
    st = res.setdefault("relax_stats", {}).setdefault(f"2S={spin}", {"relaxed": 0, "fallback": 0})
    cur = build()
    for d in DGRID:
        dk = f"{d:.2f}"
        if g.get(dk, {}).get("xyz"):
            cur = [[el, tuple(xyz)] for el, xyz in g[dk]["xyz"]]; continue
        t0 = time.time()
        geo, relaxed = crelax(cur, spin, o_i, h_j, d)
        g[dk] = {"xyz": geo, "relaxed": relaxed}
        st["relaxed" if relaxed else "fallback"] += 1
        cur = geo
        say(f"  relax 2S={spin} d={dk}: relaxed={relaxed} ({round(time.time()-t0,1)}s)")
        save()
    return g

def barrier_at_field(f_au, spin, o_i, h_j, store, save, geoms=None):
    """Скан d(O–H): SP на релакс-геометриях (RELAX) или жёсткий скан; warm-chain плотности.
    Ядерный член поля включён в scf()."""
    blk = store.setdefault(f"{f_au:+.3f}", {}).setdefault(f"2S={spin}", {"points": {}})
    base = build(); dm = None
    for d in DGRID:
        dk = f"{d:.2f}"
        if blk["points"].get(dk, {}).get("e") is not None:
            continue
        t0 = time.time()
        try:
            atoms = ([[el, tuple(xyz)] for el, xyz in geoms[dk]["xyz"]] if geoms
                     else place_H(base, o_i, h_j, d))
            e, conv, dm = scf(atoms, spin, f_au, dm0=dm, xc=SP_XC, eps=SOLV_EPS)
        except Exception as exc:
            e, conv = None, False
            say(f"      pt d={d} FAILED ({type(exc).__name__})")
        blk["points"][dk] = {"e": round(e, 6) if e is not None else None, "converged": conv}
        say(f"    F={f_au:+.3f} 2S={spin} d={dk}: E={e if e is None else round(e,5)} conv={conv} ({round(time.time()-t0,1)}s)")
        save()
    pts = sorted((float(k), v["e"]) for k, v in blk["points"].items() if v.get("converged") and v.get("e") is not None)
    if len(pts) >= 3:
        E = np.array([p[1] for p in pts])  # pts по d возрастанию; реагент = max d = pts[-1]
        blk["barrier_eV"] = round(float((E.max() - E[-1]) * EV), 4)  # forward: TS − reagent
    return blk.get("barrier_eV")

def main(a):
    o_i, h_j, nat = idx()
    nel = Z(METAL) + 8 + (10 if SUB == "ch4" else 18) - CHARGE   # [MO]+ + CH4(10e) | CH3OH(18e)
    spins = [s for s in SPINS if s % 2 == nel % 2] or [nel % 2]
    say(f"spin-field OXIDATION [{METAL}O]+ + {SUB.upper()} (oxo HAT) — xc={XC} sp_xc={SP_XC} eps={SOLV_EPS} "
        f"nat={nat} O_dist={o_i} H_abs={h_j} spins={spins} threads={lib.num_threads()}")
    res = json.load(open(OUT)) if os.path.exists(OUT) else {}
    res.setdefault("meta", {
        "what": f"C-H abstraction barrier of {SUB.upper()} by [{METAL}O]+ oxo vs axial field; "
                + ("RELAXED distinguished-coordinate d(O-H) scan (geomeTRIC frozen-distance)"
                   if RELAX else "RIGID d(O-H) scan, density warm-chain")
                + "; forward barrier = E(TS)-E(reactant, H-on-C at 2.6 A); "
                  "does field move the HAT barrier, and does it move product over-oxidation equally?",
        "metal": METAL, "xc": XC, "sp_xc": SP_XC, "substrate": SUB, "solv_eps": SOLV_EPS,
        "method": (f"{SP_XC}//{XC}" if SP_XC else XC) + (f" + ddCOSMO(eps={SOLV_EPS})//vacuum-geo" if SOLV_EPS else ""),
        "fields_au": FIELDS,
        "limits": "distinguished-coord = upper-bound barrier (not true saddle); "
                  "functional-sensitive spins/barriers; cluster, gas phase; "
                  "geometries relaxed at F=0 (field on energies only — pure electronic "
                  "Stark effect; frame: M at origin, M-O axis = z); nuclear field term "
                  "-F*sum(Z_A z_A) included; v2",
    })
    res.setdefault("scan", {})
    fields = [0.0] if a.smoke else FIELDS
    sp = [spins[0]] if a.smoke else spins
    rows = []
    geoms_by_spin = {}
    if RELAX:
        for s in sp:                              # релакс-цепочка F=0 один раз на спин; SP-поля дальше
            geoms_by_spin[s] = relax_chain(s, o_i, h_j, res, lambda: atomic_save(res))
        st = res.get("relax_stats", {})
        say(f"  relax_stats: {st}")
    for f in fields:
        best = None
        for s in sp:
            b = barrier_at_field(f, s, o_i, h_j, res["scan"], lambda: atomic_save(res),
                                 geoms=geoms_by_spin.get(s))
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
