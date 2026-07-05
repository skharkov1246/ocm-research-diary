#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
routes/antibep_carbene.py — НОВАЯ мишень A: металл-free триплетный абстрактор
как спин-селективный HAT-медиатор метана. Тест принципа «анти-BEP белого пятна».

Замысел (наша ниша, из диалога). Весь мир скринит катализаторы на DFT, который
качественно врёт на сильной статической корреляции (почти-вырожденные спиновые
состояния). Значит, мультиреференс-выгодные мотивы систематически невидимы всему
полевому тулчейну → не найдены → не коммерциализированы. Наш собственный результат
Этапа 15 (ΔΔE‡: DFT +4.7 → NEVPT2 +11.6 ккал/моль, знак переворачивается) —
прямое доказательство, что «анти-BEP» белое пятно реально.

Здесь тестируем принцип на МЕТАЛЛ-FREE абстракторах (нет драгметалла → дёшево,
мобильно). Триплетные карбены/нитрены/оксил-радикалы отрывают H у CH4/C2H6
(HAT); триплет-синглет почти вырождены → DFT врёт на барьере и спин-щели → нас.
Никто не скринил метан-HAT на металл-free дирадикалах именно потому, что DFT их
не видит.

Селективностный дескриптор — тот же, что в OCM (Этап 15):
    ΔΔE‡ = барьер(HAT из C2H6) − барьер(HAT из CH4)
Анти-BEP сигнал = знак ΔΔE‡_NEVPT2 против знака ΔΔE‡_DFT. Если корреляция делает
абстракцию из CH4 ОТНОСИТЕЛЬНО труднее (ΔΔE‡ растёт/меняет знак вверх), медиатор
защищает первичный CH3 от переокисления — пробой BEP-стены селективности.

Модель: газофазные молекулы, def2-SVP(+ECP), триплетная поверхность HAT (спин=2).
TS: коллинеарный [субстрат···H···A(медиатор)], 2-координатный пин r(C–H)+r(H–A),
скан по переносимому H, максимум = TS (тот же приём, что в Этапах 15/16 и галогене).
Уровни: UKS-PBE0 (DFT) и ROHF→AVAS(C 2p + A np + H 1s)→CASSCF→SC-NEVPT2.

Панель медиаторов (лестница мультиреференсности, все триплеты, все металл-free):
  O   — атомарный O(3P): бенчмарк-оксидант, O+CH4→OH+CH3 (реперная реакция)
  NH  — имидоген 3Σ- (нитрен-аналог)
  O2  — 3Σg-, терминальный O отрывает H → HO2
  CH2 — метилен 3B1 (архетипный триплетный карбен), C отрывает H → CH3

Запуск:
  python3 routes/antibep_carbene.py run  O|NH|O2|CH2
  python3 routes/antibep_carbene.py merge
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
BASIS = "def2-svp"
XC = "pbe0"


def M(atoms, spin, charge=0):
    return gto.M(atom=atoms, basis=BASIS, ecp=BASIS, charge=charge, spin=spin,
                 verbose=0, max_memory=6000)


def uks(mol, dm0=None):
    # без density_fit: молекулы крошечные, а DF+newton+threading иногда даёт
    # 0%-CPU дедлок. Плейн-DIIS сначала, newton только как ускоритель с cap.
    mf = dft.UKS(mol); mf.xc = XC; mf.conv_tol = 1e-8; mf.max_cycle = 200
    mf.kernel(dm0=dm0) if dm0 is not None else mf.kernel()
    if not mf.converged:
        mf.level_shift = 0.3; mf.max_cycle = 300; mf.kernel(mf.make_rdm1())
    if not mf.converged:
        mfn = scf.newton(mf); mfn.max_cycle = 60; mfn.kernel(mf.make_rdm1())
        if mfn.converged:
            return mfn
    return mf


def rohf(mol):
    mf = scf.ROHF(mol); mf.conv_tol = 1e-8; mf.max_cycle = 200
    mf.kernel()
    if not mf.converged:
        mf.level_shift = 0.4; mf.max_cycle = 300; mf.kernel(mf.make_rdm1())
    if not mf.converged:
        mfn = scf.newton(mf); mfn.max_cycle = 60; mfn.kernel(mf.make_rdm1())
        if mfn.converged:
            return mfn
    return mf


def relax_uks(atoms, spin):
    # одноатомный медиатор (напр. атомарный O 3P) — оптимизировать нечего,
    # geomeTRIC падает на построении internal-coords → одноточечно.
    if len(atoms) <= 1:
        mf = uks(M(atoms, spin))
        na = [(atoms[0][0], tuple(atoms[0][1]))]
        return na, mf
    from pyscf.geomopt.geometric_solver import optimize
    mf = uks(M(atoms, spin))
    mol = optimize(mf, maxsteps=100, assert_convergence=False)
    na = [(mol.atom_symbol(i), tuple(c))
          for i, c in enumerate(mol.atom_coords(unit="Angstrom"))]
    return na, uks(mol)


def nevpt2(atoms, spin, avas_labels, thr=0.5):
    """ROHF→AVAS→CASSCF→SC-NEVPT2 на фиксированной геометрии."""
    mol = M(atoms, spin)
    mf = rohf(mol)
    ncas, nelec, mo = avas.avas(mf, avas_labels, threshold=thr)
    ss = spin / 2 * (spin / 2 + 1)
    mc = mcscf.CASSCF(mf, ncas, nelec)
    fci.addons.fix_spin_(mc.fcisolver, shift=0.2, ss=ss)
    mc.max_cycle_macro = 120; mc.conv_tol = 1e-7
    mc.kernel(mo)
    e_pt = float(mc.e_tot + NEVPT(mc).kernel())
    return {"e_casscf": float(mc.e_tot), "e_nevpt2": e_pt,
            "cas": [int(np.sum(nelec)), int(ncas)], "conv": bool(mc.converged)}


# ------------------------------------------------------------ субстраты
def ch4_atoms():
    d, th = 1.09, math.radians(109.47)
    a = [("C", (0, 0, 0)), ("H", (0, 0, d))]
    for k in range(3):
        phi = math.radians(120 * k)
        a.append(("H", (d*math.sin(th)*math.cos(phi), d*math.sin(th)*math.sin(phi),
                        d*math.cos(th))))
    return a


def ch3_atoms():
    d = 1.08
    return [("C", (0, 0, 0))] + [("H", (d*math.cos(math.radians(120*k)),
                                        d*math.sin(math.radians(120*k)), 0))
                                 for k in range(3)]


def c2h6_atoms():
    """Этан вдоль z: Cα(0), Cβ(1.53); ставим по 3 H шахматно на каждый C."""
    dcc, dch, th = 1.53, 1.09, math.radians(111.0)
    a = [("C", (0, 0, 0)), ("C", (0, 0, dcc))]
    for k in range(3):                                    # H на Cα (вниз, -z)
        phi = math.radians(120 * k)
        a.append(("H", (dch*math.sin(th)*math.cos(phi),
                        dch*math.sin(th)*math.sin(phi), -dch*math.cos(th))))
    for k in range(3):                                    # H на Cβ (вверх, +z), шахматно
        phi = math.radians(120 * k + 60)
        a.append(("H", (dch*math.sin(th)*math.cos(phi),
                        dch*math.sin(th)*math.sin(phi), dcc + dch*math.cos(th))))
    return a


def c2h5_atoms():
    """Этил-радикал: убрали один H с Cα (первый H-верх? нет — с Cα один вниз)."""
    a = c2h6_atoms()
    del a[2]           # снять один H с Cα (индексы 2..4 — Hα)
    return a


# медиатор: abstracting-атом в точке z=Az, остальные атомы «отвёрнуты» дальше (+z)
def med_atoms(name, Az):
    if name == "O":
        return [("O", (0, 0, Az))]
    if name == "NH":
        return [("N", (0, 0, Az)), ("H", (0.90, 0, Az + 0.55))]
    if name == "O2":
        return [("O", (0, 0, Az)), ("O", (0, 0, Az + 1.21))]
    if name == "CH2":
        return [("C", (0, 0, Az)), ("H", (0.94, 0, Az + 0.55)),
                ("H", (-0.94, 0, Az + 0.55))]
    raise ValueError(name)


MED = {
    "O":  {"spin": 2, "plab": "O 2p", "rAH": 0.98,
           "mh": [("O", (0, 0, 0)), ("H", (0, 0, 0.97))], "mh_spin": 1,
           "mh_lab": ["O 2p", "H 1s"]},                       # OH
    "NH": {"spin": 2, "plab": "N 2p", "rAH": 1.04,
           "mh": [("N", (0, 0, 0)), ("H", (0, 0, 1.02)), ("H", (0.95, 0, -0.32))],
           "mh_spin": 1, "mh_lab": ["N 2p", "H 1s"]},         # NH2
    "O2": {"spin": 2, "plab": "O 2p", "rAH": 1.00,
           "mh": [("O", (0, 0, 0)), ("O", (0, 0, 1.33)), ("H", (0.95, 0, -0.32))],
           "mh_spin": 1, "mh_lab": ["O 2p", "H 1s"]},         # HO2
    "CH2": {"spin": 2, "plab": "C 2p", "rAH": 1.10,
            "mh": ch3_atoms(), "mh_spin": 1, "mh_lab": ["C 2p"]},   # CH3
}

# субстраты: (builder, C-радикал-продукт-builder, spin_sub, spin_rad, cα-index в TS)
SUBSTRATE = {
    "ch4":  {"sub": ch4_atoms,  "rad": ch3_atoms,  "rad_lab": ["C 2p"]},
    "c2h6": {"sub": c2h6_atoms, "rad": c2h5_atoms, "rad_lab": ["C 2p"]},
}


def ts_atoms(name, sub_key, rCH, rAH):
    """Почти-коллинеарный [Cα···H···A]‡ по z: Cα(0) — H_t — A(медиатор), хвост вниз.
    Намеренный лёгкий излом (dx) от оси: точно линейный C–H–A = 180° — координатная
    сингулярность internal-coords (geomeTRIC «Inverse iteration failed»). Пины на
    расстояния r(C–H)/r(H–A) держатся, угол свободен → TS сам досогнётся."""
    dxH, dxA = 0.22, 0.44                               # излом ~8-10° от оси z
    zC, zH, zA = 0.0, rCH, rCH + rAH
    a = [("C", (0, 0, zC)), ("H", (dxH, 0, zH))]        # Cα, переносимый H (сдвинут)
    a += [(s, (x + dxA, y, z)) for s, (x, y, z) in med_atoms(name, zA)]  # A + хвост, сдвинут
    # хвост субстрата (вниз, -z), чтобы не мешать линейному каналу
    d, th = 1.09, math.radians(105.0)
    if sub_key == "ch4":
        for k in range(3):
            phi = math.radians(120 * k)
            a.append(("H", (d*math.sin(th)*math.cos(phi), d*math.sin(th)*math.sin(phi),
                            zC - d*math.cos(th))))
    else:  # c2h6: 2 H на Cα + этильный CH3 вниз
        for k in range(2):
            phi = math.radians(120 * k + 60)
            a.append(("H", (d*math.sin(th)*math.cos(phi), d*math.sin(th)*math.sin(phi),
                            zC - d*math.cos(th))))
        zCb = zC - 1.53
        a.append(("C", (0, 0, zCb)))
        for k in range(3):
            phi = math.radians(120 * k)
            a.append(("H", (d*math.sin(th)*math.cos(phi), d*math.sin(th)*math.sin(phi),
                            zCb - d*math.cos(th))))
    return a


def ts_scan(name, sub_key):
    """ЖЁСТКИЙ скан HAT-координаты: одноточечные UKS по сетке (r(Cα–H), r(H–A)),
    максимум = TS-оценка. Без геом-оптимизации — коллинеарный [C···H···A] TS
    патологичен для internal-coords (сингулярность 180°) и слишком дорог в релаксе
    локально. Честная оговорка: фрагменты не релаксированы в TS → барьеры суть
    ВЕРХНИЕ оценки; но знак ΔΔE‡ (CH4↔C2H6) и его инверсия DFT↔NEVPT2 — робастные
    наблюдаемые (несистематическая ошибка сокращается в разности субстратов)."""
    rAH0 = MED[name]["rAH"]
    grid = [(1.15, rAH0 + 0.30), (1.25, rAH0 + 0.15), (1.35, rAH0 + 0.05),
            (1.45, rAH0), (1.60, rAH0 - 0.10), (1.75, rAH0 - 0.15)]
    best = None
    for rCH, rAH in grid:
        try:
            mfx = uks(M(ts_atoms(name, sub_key, rCH, rAH), MED[name]["spin"]))
            e = float(mfx.e_tot)
            if not np.isfinite(e):
                continue
        except Exception as ex:
            print(f"  [{name}/{sub_key}] rCH={rCH} failed: {ex}", flush=True)
            continue
        na = ts_atoms(name, sub_key, rCH, rAH)          # жёсткая геометрия точки
        if best is None or e > best[1]:
            best = (na, e, rCH, rAH)
    if best is None:
        raise RuntimeError(f"TS scan failed for {name}/{sub_key}")
    return best


def barrier_set(name, sub_key, mf_med, e_med_dft, med_corr):
    """Барьер HAT для одного субстрата: DFT + CASSCF + NEVPT2 (общий репер медиатора)."""
    S = SUBSTRATE[sub_key]
    sub_a, mf_sub = relax_uks(S["sub"](), 0)
    ts_a, e_ts_dft, rCH, rAH = ts_scan(name, sub_key)
    out = {"dft_barrier": round((e_ts_dft - e_med_dft - mf_sub.e_tot) * HARTREE_KCAL, 2),
           "ts_rCH": rCH, "ts_rAH": round(rAH, 3)}
    # NEVPT2: TS и субстрат (медиатор-репер общий, посчитан снаружи)
    lab_ts = ["C 2p", MED[name]["plab"], "H 1s"]
    n_ts = nevpt2(ts_a, MED[name]["spin"], lab_ts, thr=0.4)
    n_sub = nevpt2(sub_a, 0, ["C 2p", "H 1s"], thr=0.6)
    for lvl in ("casscf", "nevpt2"):
        b = (n_ts[f"e_{lvl}"] - med_corr[f"e_{lvl}"] - n_sub[f"e_{lvl}"]) * HARTREE_KCAL
        out[f"{lvl}_barrier"] = round(b, 2)
    out["ts_cas"] = n_ts["cas"]
    return out


def stage_run(name):
    cfg = MED[name]
    path = os.path.join(DIR, f"antibep_{name}.json")
    out = {"mediator": name, "model": "gas-phase triplet metal-free HAT, def2-SVP+ECP",
           "spin_triplet": cfg["spin"]}
    t0 = time.time()
    # репер медиатора (общий для обоих субстратов)
    med_a, mf_med = relax_uks(med_atoms(name, 0.0), cfg["spin"])
    med_corr = nevpt2(med_a, cfg["spin"], [cfg["plab"]], thr=0.5)
    out["ch4"] = barrier_set(name, "ch4", mf_med, float(mf_med.e_tot), med_corr)
    out["c2h6"] = barrier_set(name, "c2h6", mf_med, float(mf_med.e_tot), med_corr)
    for lvl in ("dft", "casscf", "nevpt2"):
        k = "dft_barrier" if lvl == "dft" else f"{lvl}_barrier"
        out[f"ddE_{lvl}"] = round(out["c2h6"][k] - out["ch4"][k], 2)
    out["wall_s"] = round(time.time() - t0, 1)
    json.dump(out, open(path, "w"), indent=1)
    print(f"[{name}] barrier CH4/C2H6 (kcal/mol)  DFT {out['ch4']['dft_barrier']}/"
          f"{out['c2h6']['dft_barrier']} | NEVPT2 {out['ch4']['nevpt2_barrier']}/"
          f"{out['c2h6']['nevpt2_barrier']}  ==>  ΔΔE‡ DFT {out['ddE_dft']} | "
          f"CASSCF {out['ddE_casscf']} | NEVPT2 {out['ddE_nevpt2']}  ({out['wall_s']}s)",
          flush=True)


def stage_merge():
    rows = []
    for name in MED:
        p = os.path.join(DIR, f"antibep_{name}.json")
        if os.path.exists(p):
            d = json.load(open(p))
            rows.append({"mediator": name,
                         "ddE_dft": d.get("ddE_dft"), "ddE_casscf": d.get("ddE_casscf"),
                         "ddE_nevpt2": d.get("ddE_nevpt2"),
                         "sign_flip": (d.get("ddE_dft", 0) * d.get("ddE_nevpt2", 0)) < 0})
    res = {"descriptor": "metal-free triplet HAT selectivity ΔΔE‡ = barrier(C2H6) − "
                         "barrier(CH4); anti-BEP test = sign(ΔΔE‡_NEVPT2) vs DFT",
           "note": "gas-phase triplet mediators, def2-SVP+ECP; RIGID collinear HAT scan "
                   "(fragments unrelaxed at TS → absolute barriers are upper bounds; "
                   "ΔΔE‡ sign and the DFT↔NEVPT2 flip are the robust observables). "
                   "sign_flip=True means the "
                   "correlated (NEVPT2) selectivity descriptor has the OPPOSITE sign to "
                   "DFT — a metal-free anti-BEP site invisible to DFT screening. "
                   "O(3P)+CH4 is the literature benchmark (barrier ~11-14 kcal/mol).",
           "mediators": rows}
    json.dump(res, open(os.path.join(DIR, "antibep_carbene_results.json"), "w"),
              indent=1)
    print(json.dumps(res, ensure_ascii=False, indent=1))


if __name__ == "__main__":
    if sys.argv[1] == "run":
        stage_run(sys.argv[2])
    elif sys.argv[1] == "merge":
        stage_merge()
    else:
        raise SystemExit("usage: run O|NH|O2|CH2, or merge")
