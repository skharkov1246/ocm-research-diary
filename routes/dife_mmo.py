#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
routes/dife_mmo.py — изобретение I3 / пробел карты L1 №12: антиферромагнитный
ди-железный бис-μ-оксо «алмаз» по мотиву интермедиата Q фермента метан-
монооксигеназы (sMMO). Природный эталон: CH4 -> CH3OH при КОМНАТНОЙ температуре
со 100% селективностью на ди-Fe(IV)=O центре. Это одновременно (карта) эталон,
к которому все стремятся, и (глубина) субстрат нашего динуклеарного изобретения.

Почему это КВИНТЭССЕНЦИЯ нашей ниши. Ядро [Fe2(μ-O)2]4+ — два Fe(IV) d^4 (HS S=2),
антиферромагнитно связанные в общий низкоспиновый терм. AFM-связка = сильнейшая
статическая корреляция (почти-вырожденные HS/BS-состояния): DFT/broken-symmetry
даёт неверную спин-щель и барьер HAT на десятки ккал/моль; правильно это описывает
только мультиреференс (CASSCF/NEVPT2 ловит AFM-синглет естественно, без BS-костыля).
Здесь наш стек не «полезен», а незаменим.

Минимальная обоснованная модель (стенд-ин, не белок): алмазное ядро [Fe2(μ-O)2], каждый
Fe кэпирован 2×OH до ~октаэдра, Fe(IV). Нейтрально: 2·Fe(IV)=+8, 2·μ-O2-=-4,
4·OH-=-4 → 0. HAT: пин r(C-H)+r(Oμ-H) 2D-скан к мостиковому O (тот же приём, что
Mn=O в Этапах 15/16). Спин: считаем FM (2S=8) и AFM-синглет (CASSCF S=0);
барьер и его DFT↔NEVPT2 расхождение — главный выход (де-риск + связь с I1/I3).

Уровни: UKS-PBE0 (геометрия + DFT-барьер, FM-поверхность — tractable);
ROHF→AVAS(2×Fe 3d + реагирующий Oμ 2p)→CASSCF(fix_spin S)→SC-NEVPT2 (AFM-с оговоркой).
Fan-out по точкам пути.

Запуск:
  python3 routes/dife_mmo.py geom            # 2D-пин HAT-профиль (UKS-PBE0, FM)
  python3 routes/dife_mmo.py profile <spin2S>  # CASSCF/NEVPT2 по точкам (0=AFM,8=FM)
  python3 routes/dife_mmo.py sanity          # дешёвая локальная проверка (SCF+AVAS)
  python3 routes/dife_mmo.py merge
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
BASIS = os.environ.get("DIFE_BASIS", "def2-svp")
XC = "pbe0"
# сетка HAT: r(C-H) растёт, r(Oμ-H) падает — переносимый H от метана к мостик. O
GRID = [(None, 2.4), (1.15, 1.6), (1.25, 1.35), (1.40, 1.15), (1.60, 1.02)]


def core_atoms(rCH=None, rOH=None):
    """[Fe2(μ-O)2(OH)4] алмаз по x; при заданных r — добавлен подлетающий CH4/CH3+H
    к мостиковому O_a (индекс 2). Индексы: Fe1=0 Fe2=1 O_a=2 O_b=3, далее OH-кэпы."""
    a = [("Fe", (-1.30, 0.0, 0.0)), ("Fe", (1.30, 0.0, 0.0)),
         ("O", (0.0, 1.05, 0.0)), ("O", (0.0, -1.05, 0.0))]          # μ-оксо
    # терминальные OH: по 2 на Fe, растянуты наружу
    for (fx, s) in ((-1.30, -1), (1.30, 1)):
        a += [("O", (fx + s*1.55, 0.75, 1.05)), ("H", (fx + s*2.05, 1.05, 1.45)),
              ("O", (fx + s*1.55, -0.75, -1.05)), ("H", (fx + s*2.05, -1.05, -1.45))]
    if rCH is not None:
        # метил + переносимый H над O_a (по +y), коллинеарно O_a···H···C
        zc = 1.05 + rOH + rCH
        a.append(("H", (0.0, 1.05 + rOH, 0.0)))                       # H_t
        d, th = 1.09, math.radians(107.0)
        a.append(("C", (0.0, zc, 0.0)))
        for k in range(3):
            phi = math.radians(120 * k)
            a.append(("C" if False else "H",
                      (d*math.sin(th)*math.cos(phi), zc + d*math.cos(th),
                       d*math.sin(th)*math.sin(phi))))
    return a


def build(atoms, spin, charge=0):
    return gto.M(atom=atoms, basis=BASIS, ecp=BASIS, charge=charge, spin=spin,
                 verbose=0, max_memory=16000)


def uks(atoms, spin, dm0=None):
    mol = build(atoms, spin); mol.max_memory = 16000
    mf = dft.UKS(mol); mf.xc = XC; mf.conv_tol = 1e-7; mf.max_cycle = 300
    mf.kernel(dm0=dm0)
    if not mf.converged:
        mf.level_shift = 0.4; mf.max_cycle = 500; mf.kernel(mf.make_rdm1())
    return mf


def rohf(mol):
    mf = scf.ROHF(mol); mf.conv_tol = 1e-7; mf.max_cycle = 300; mf.kernel()
    if not mf.converged:
        mf.level_shift = 0.5; mf.max_cycle = 500; mf.kernel(mf.make_rdm1())
    return mf


def constrained_opt(atoms, spin, rCH, rOH):
    from pyscf.geomopt.geometric_solver import optimize
    cf = os.path.join(DIR, f"_dife_{rCH}_{rOH}.txt")
    # пин r(O_a–H_t) = distance 3 (idx2+1=3-я? 1-based) и r(H_t–C); индексы 1-based:
    # O_a=3, H_t=13, C=14 (см. core_atoms порядок: 4 core +8 caps =12, +H_t=13,+C=14)
    open(cf, "w").write(f"$set\ndistance 13 14 {rCH:.4f}\ndistance 3 13 {rOH:.4f}\n")
    mf = uks(atoms, spin)
    try:
        mol = optimize(mf, constraints=cf, maxsteps=70, assert_convergence=False)
    finally:
        if os.path.exists(cf):
            os.remove(cf)
    na = [(mol.atom_symbol(i), tuple(c))
          for i, c in enumerate(mol.atom_coords(unit="Angstrom"))]
    return na, uks(na, spin)


def stage_geom():
    """FM-поверхность (2S=8), 2D-пин HAT-профиль."""
    spin = 8
    out = {"model": "[Fe2(mu-O)2(OH)4] diamond, FM 2S=8", "points": []}
    e0 = None
    for i, (rCH, rOH) in enumerate(GRID):
        t0 = time.time()
        if rCH is None:
            atoms = core_atoms(); mf = uks(atoms, spin)
            na = atoms
        else:
            na, mf = constrained_opt(core_atoms(rCH, rOH), spin, rCH, rOH)
        e = float(mf.e_tot)
        if e0 is None:
            e0 = e
        out["points"].append({"i": i, "rCH": rCH, "rOH": rOH, "e_h": e,
                              "rel_kcal": round((e - e0) * HARTREE_KCAL, 2),
                              "s2": round(float(mf.spin_square()[0]), 2),
                              "conv": bool(mf.converged), "atoms": na,
                              "wall_s": round(time.time() - t0, 1)})
        json.dump(out, open(os.path.join(DIR, "dife_geom.json"), "w"), indent=1)
        print(f"[geom] p{i} rCH={rCH} rOH={rOH} rel={out['points'][-1]['rel_kcal']} "
              f"kcal <S2>={out['points'][-1]['s2']} ({out['points'][-1]['wall_s']}s)",
              flush=True)
    es = [p["rel_kcal"] for p in out["points"]]
    im = es.index(max(es))
    out["dft_barrier_kcal"] = max(es)
    out["ts_index"] = im
    out["ts_interior"] = 0 < im < len(es) - 1
    json.dump(out, open(os.path.join(DIR, "dife_geom.json"), "w"), indent=1)
    print(f"[geom] DFT HAT barrier {max(es):.2f} kcal/mol at p{im} "
          f"(interior={out['ts_interior']})", flush=True)


def stage_profile(spin2s):
    spin = int(spin2s)
    g = json.load(open(os.path.join(DIR, "dife_geom.json")))
    out = {"spin_2S": spin, "points": []}
    for p in g["points"]:
        t0 = time.time()
        atoms = [(s, tuple(c)) for s, c in p["atoms"]]
        mol = build(atoms, spin); mol.max_memory = 24000
        mf = rohf(mol)
        # фикс-состав AVAS: ТОЛЬКО 2×Fe 3d. Полное «Fe 3d + O 2p» дало CAS(18,16) —
        # AFM-синглет там ~1e8 детерминантов, FCI-стена (ди-Fe MMO = DMRG-grade).
        # Fe-3d-only (~CAS(8,10)) ловит суть AFM Fe–Fe связи; корреляцию рвущейся
        # O–H связи NEVPT2 добирает пертурбативно. Честная screening-усечёнка.
        labels = ["Fe 3d"]
        ncas, nelec, mo = avas.avas(mf, labels, threshold=0.15)
        ss = spin / 2 * (spin / 2 + 1)
        mc = mcscf.CASSCF(mf, ncas, nelec)
        fci.addons.fix_spin_(mc.fcisolver, shift=0.3, ss=ss)
        mc.max_cycle_macro = 100; mc.conv_tol = 1e-6
        mc.kernel(mo)
        e_pt = float(mc.e_tot + NEVPT(mc).kernel())
        out["points"].append({"i": p["i"], "rCH": p["rCH"], "rOH": p["rOH"],
                              "cas": [int(np.sum(nelec)), int(ncas)],
                              "e_casscf_h": float(mc.e_tot), "e_nevpt2_h": e_pt,
                              "casscf_conv": bool(mc.converged),
                              "wall_s": round(time.time() - t0, 1)})
        json.dump(out, open(os.path.join(DIR, f"dife_profile_s{spin}.json"), "w"),
                  indent=1)
        print(f"[prof s{spin}] p{p['i']} CAS{out['points'][-1]['cas']} "
              f"CASSCF={mc.e_tot:.5f} NEVPT2={e_pt:.5f} "
              f"({out['points'][-1]['wall_s']}s)", flush=True)


def stage_sanity():
    """Дешёвая локальная проверка: строится ли ядро, сходится ли SCF, что даёт AVAS."""
    atoms = core_atoms()
    for spin in (8, 0):
        mol = build(atoms, spin); mol.max_memory = 8000
        mf = scf.ROHF(mol); mf.max_cycle = 80; mf.kernel()
        try:
            ncas, nelec, mo = avas.avas(mf, ["Fe 3d"], threshold=0.15)
            print(f"[sanity] 2S={spin} ROHF conv={mf.converged} E={mf.e_tot:.4f} "
                  f"AVAS CAS=({int(np.sum(nelec))},{int(ncas)})", flush=True)
        except Exception as e:
            print(f"[sanity] 2S={spin} AVAS FAIL: {e}", flush=True)


def stage_merge():
    g = json.load(open(os.path.join(DIR, "dife_geom.json")))
    res = {"invention": "I3 di-Fe MMO-Q antiferromagnetic diamond HAT",
           "dft_barrier_kcal": g.get("dft_barrier_kcal"),
           "dft_ts_interior": g.get("ts_interior"), "levels": {}}
    for spin in (0, 8):
        p = os.path.join(DIR, f"dife_profile_s{spin}.json")
        if os.path.exists(p):
            d = json.load(open(p)); pts = d["points"]
            e0c, e0n = pts[0]["e_casscf_h"], pts[0]["e_nevpt2_h"]
            bc = max((x["e_casscf_h"] - e0c) for x in pts) * HARTREE_KCAL
            bn = max((x["e_nevpt2_h"] - e0n) for x in pts) * HARTREE_KCAL
            res["levels"][f"2S={spin}"] = {"casscf_barrier": round(bc, 2),
                                           "nevpt2_barrier": round(bn, 2)}
    res["note"] = ("AFM (2S=0) vs FM (2S=8) HAT barrier + DFT-vs-NEVPT2 gap. The "
                   "correlated AFM description is where DFT/broken-symmetry fails; "
                   "our stack is not merely useful but necessary here. Cluster "
                   "stand-in, coarse grid — screening.")
    json.dump(res, open(os.path.join(DIR, "dife_mmo_results.json"), "w"), indent=1)
    print(json.dumps(res, ensure_ascii=False, indent=1))


if __name__ == "__main__":
    s = sys.argv[1]
    if s == "geom":
        stage_geom()
    elif s == "profile":
        stage_profile(sys.argv[2])
    elif s == "sanity":
        stage_sanity()
    elif s == "merge":
        stage_merge()
    else:
        raise SystemExit("usage: geom | profile <2S> | sanity | merge")
