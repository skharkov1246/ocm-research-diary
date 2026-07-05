#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
routes/pyrolysis_descriptor.py — дескриптор пиролиза метана в расплаве металла
(маршрут 3 разведки): CH4 -> C(тв) + 2 H2 в жидком био-металлическом сплаве
(активный металл + низкоплавкий «растворитель»: Ni-Bi, Cu-Bi, Ni-In, Pt-Bi …).
Продукт — водород + ТВЁРДЫЙ углерод, CO2-free; самый bankable маршрут по
экономике (валидированный LCOH $1.3-3/кг H2, коммерческие скид-модули).

Почему это «наша ниша» (из METHANE_ROUTES_RESEARCH.md / _ECONOMICS.md). Селект-
ивность/жизнеспособность расплава управляется ОДНИМ балансом Сабатье на паре
активный-металл/растворитель:
  • слишком сильная связь M–C  → карбид/закоксовывание, расплав грязнет;
  • слишком слабая              → метан не активируется, нет конверсии.
Нужен ОПТИМУМ энергии связи углерода E_bind(C). Связь переходный-металл–C на
частично заполненных d-оболочках биметаллической пары — там, где DFT врёт со
спин-состояниями и прочностью связи M–адсорбат (функционал-чувствительно на
десятки ккал/моль). Прямой повод для CASSCF/NEVPT2 — как в Этапах 15/16 и в
вакансионном волкане (маршрут 2).

Минимальная честная модель (стенд-ин, НЕ периодический расплав — как и всюду
в этом дневнике): гетеро-биметаллический ДИМЕР M–Sol как «пара расплава», и на
нём — два дескриптора:
  1) E_bind(C)      = E(Sol–M–C) − E(M–Sol) − E(C•)      (ось коксования, главная)
  2) ΔE_diss(CH4)   = E(Sol–M(CH3)(H)) − E(M–Sol) − E(CH4) (ось активации, вторая)
Вулкан: хотим достаточную активацию, но НЕ настолько сильное связывание C, чтобы
образовался стабильный карбид. Наш репер — Ni-Bi (канонический рабочий сплав).

Уровни: UKS-PBE0/def2-SVP(+ECP для тяжёлых) — геометрия + DFT-энергетика (дёшево,
все пары); ROHF→AVAS(M nd + Sol вал. + C 2p)→CASSCF→SC-NEVPT2 — коррелированная
E_bind(C) (де-риск, на AWS). Спин каждого фрагмента — берём нижний из скана
допустимых по чётности мультиплетностей (честно, как в соседних треках).
Fan-out: одна пара на запуск → embarrassingly parallel.

Запуск:
  python3 routes/pyrolysis_descriptor.py dft   <PAIR>   # геометрии + DFT-дескрипторы
  python3 routes/pyrolysis_descriptor.py corr  <PAIR>   # + CASSCF/NEVPT2 E_bind(C) (AWS)
  python3 routes/pyrolysis_descriptor.py merge          # свести вулкан по всем парам
  PAIR из PANEL ниже, напр. Ni-Bi, Cu-Bi, Ni-In, Pt-Bi, Ni-Sn, Fe-Sn, Co-Ga, Pd-Bi
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

# Панель пар «активный металл – низкоплавкий растворитель» (лит. рабочие расплавы).
# orb — валентная метка для AVAS; d_MM/d_MC — стартовые длины связей (Å).
# nd/np — какую валентную оболочку растворителя тянуть в активное пространство.
ACTIVE = {
    "Ni": {"orb": "Ni 3d", "d": 3, "n": 8},
    "Cu": {"orb": "Cu 3d", "d": 3, "n": 10},
    "Pt": {"orb": "Pt 5d", "d": 5, "n": 9},
    "Pd": {"orb": "Pd 4d", "d": 4, "n": 10},
    "Fe": {"orb": "Fe 3d", "d": 3, "n": 6},
    "Co": {"orb": "Co 3d", "d": 3, "n": 7},
}
SOLVENT = {
    "Bi": {"orb": "Bi 6p"},
    "In": {"orb": "In 5p"},
    "Sn": {"orb": "Sn 5p"},
    "Ga": {"orb": "Ga 4p"},
}
PANEL = {
    "Ni-Bi": ("Ni", "Bi", 2.55, 1.75),   # репер: канонический рабочий сплав
    "Cu-Bi": ("Cu", "Bi", 2.55, 1.80),
    "Ni-In": ("Ni", "In", 2.55, 1.78),
    "Pt-Bi": ("Pt", "Bi", 2.65, 1.80),
    "Ni-Sn": ("Ni", "Sn", 2.55, 1.78),
    "Fe-Sn": ("Fe", "Sn", 2.55, 1.80),
    "Co-Ga": ("Co", "Ga", 2.45, 1.75),
    "Pd-Bi": ("Pd", "Bi", 2.65, 1.82),
    # волна 2 (Франкфурт): достроить решётку активный×растворитель
    "Pt-Sn": ("Pt", "Sn", 2.60, 1.80),
    "Pd-In": ("Pd", "In", 2.60, 1.80),
    "Cu-Sn": ("Cu", "Sn", 2.55, 1.80),
    "Fe-Bi": ("Fe", "Bi", 2.60, 1.80),
    "Ni-Ga": ("Ni", "Ga", 2.45, 1.75),
    "Co-Bi": ("Co", "Bi", 2.55, 1.78),
    "Cu-In": ("Cu", "In", 2.55, 1.80),
    "Pd-Sn": ("Pd", "Sn", 2.60, 1.80),
}
# ручки Франкфурт-прогонов: PYRO_NSPIN=3 — расширить лестницу спинов (пере-прогон
# Cu-Bi артефакта); PYRO_SUFFIX="_v2" — не затирать прежние JSON.
NSPIN = int(os.environ.get("PYRO_NSPIN", "2"))
SUF = os.environ.get("PYRO_SUFFIX", "")


def build(atoms, spin, charge=0):
    return gto.M(atom=atoms, basis=BASIS, ecp=BASIS, charge=charge, spin=spin,
                 verbose=0, max_memory=8000)


def valid_spins(atoms, charge=0, n=None):
    """Нижняя чётность-совместимая лестница мультиплетностей (2S = 0/1,+2,+4…).
    Скрин-компромисс: 2 мультиплетности (базовая + базовая+2) — обычная HS/LS
    конкуренция биметалл-димера; берём нижнюю по энергии (честная оговорка)."""
    base = None
    for s in (0, 1):
        try:
            build(atoms, s, charge)
            base = s
            break
        except Exception:
            continue
    if base is None:
        base = 0
    return [base + 2 * k for k in range(n or NSPIN)]


def uks(mol, dm0=None):
    mf = dft.UKS(mol).density_fit()
    mf.xc = XC
    mf.conv_tol = 1e-8
    mfn = scf.newton(mf)
    mfn.kernel(dm0=dm0) if dm0 is not None else mfn.kernel()
    if not mfn.converged:
        m2 = dft.UKS(mol).density_fit(); m2.xc = XC; m2.level_shift = 0.3
        m2.max_cycle = 400; m2.kernel(dm0=dm0)
        mfn = scf.newton(m2); mfn.kernel(m2.make_rdm1())
    return mfn


# PYRO_SP=1 — быстрый режим одноточечных энергий (без геом-оптимизации): для
# локальной валидации арифметики/I-O и как аварийный fallback. Голые металл-димеры
# патологичны для internal-coord оптимизатора → релакс идёт в декартовых (coordsys
# 'cart', устойчивее для крошечных систем), а при сбое — одноточечно на входной
# геометрии (честная оговорка: неполный релакс, скрин-оценка).
SINGLE_POINT = os.environ.get("PYRO_SP", "0") == "1"


def relax_lowest(atoms, charge=0, tag=""):
    """Релакс UKS по нескольким спинам, вернуть НИЖНИЙ по энергии (геометрия+E)."""
    from pyscf.geomopt.geometric_solver import optimize
    best = None
    for spin in valid_spins(atoms, charge):
        try:
            mf = uks(build(atoms, spin, charge))
            if SINGLE_POINT:
                mol = build(atoms, spin, charge); mfx = mf
            else:
                try:
                    mol = optimize(mf, maxsteps=80, assert_convergence=False,
                                   coordsys="cart")
                    mfx = uks(mol)
                except Exception as ox:
                    print(f"  [{tag}] spin={spin} opt fell back to SP: {ox}",
                          flush=True)
                    mol = build(atoms, spin, charge); mfx = mf
            e = float(mfx.e_tot)
            if not np.isfinite(e):
                raise ValueError("non-finite energy")
        except Exception as ex:
            print(f"  [{tag}] spin={spin} failed: {ex}", flush=True)
            continue
        na = [(mol.atom_symbol(i), tuple(c))
              for i, c in enumerate(mol.atom_coords(unit="Angstrom"))]
        s2 = float(mfx.spin_square()[0])
        if best is None or e < best[1]:
            best = (na, e, spin, bool(mfx.converged), s2)
    if best is None:
        raise RuntimeError(f"all spins failed for {tag}")
    return best   # (atoms, e, spin, conv, s2)


# ------------------------------------------------------------ геометрии-стартеры
def dimer(M, Sol, dMM):
    return [(M, (0.0, 0.0, 0.0)), (Sol, (0.0, 0.0, dMM))]


def dimer_C(M, Sol, dMM, dMC):
    """C атоп активного металла M, растворитель с другой стороны (линейный старт)."""
    return [(M, (0.0, 0.0, 0.0)), (Sol, (0.0, 0.0, dMM)),
            ("C", (0.0, 0.0, -dMC))]


def dimer_CH3_H(M, Sol, dMM, dMC):
    """Диссоц. адсорбция CH4: CH3 атоп M, H атоп Sol (продукт разрыва C–H)."""
    zc = -dMC
    d, th = 1.09, math.radians(110.0)
    a = [(M, (0.0, 0.0, 0.0)), (Sol, (0.0, 0.0, dMM)),
         ("C", (0.0, 0.0, zc))]
    for k in range(3):
        phi = math.radians(120 * k)
        a.append(("H", (d*math.sin(th)*math.cos(phi), d*math.sin(th)*math.sin(phi),
                        zc - d*math.cos(th))))
    a.append(("H", (0.0, 1.55, dMM + 0.55)))   # H атоп растворителя
    return a


def ch4_atoms():
    d, th = 1.09, math.radians(109.47)
    a = [("C", (0, 0, 0)), ("H", (0, 0, d))]
    for k in range(3):
        phi = math.radians(120 * k)
        a.append(("H", (d*math.sin(th)*math.cos(phi), d*math.sin(th)*math.sin(phi),
                        d*math.cos(th))))
    return a


def c_atom_energy(level="dft"):
    """E(C•) — атом углерода, триплет (3P), общий репер."""
    mol = gto.M(atom=[("C", (0, 0, 0))], basis=BASIS, spin=2, verbose=0)
    mf = uks(mol)
    if level == "dft":
        return float(mf.e_tot)
    no, ne, mo = avas.avas(mf, ["C 2p"], threshold=0.2)
    mc = mcscf.CASSCF(mf, no, ne)
    fci.addons.fix_spin_(mc.fcisolver, shift=0.2, ss=2.0)
    mc.kernel(mo)
    return float(mc.e_tot + NEVPT(mc).kernel())


def ch4_energy(level="dft"):
    a, e, spin, conv, s2 = relax_lowest(ch4_atoms(), tag="ch4")
    if level == "dft":
        return e
    no, ne, mo = avas.avas(uks(build(a, spin)), ["C 2p", "H 1s"], threshold=0.5)
    mc = mcscf.CASSCF(uks(build(a, spin)), no, ne)
    mc.kernel(mo)
    return float(mc.e_tot + NEVPT(mc).kernel())


def rohf(mol):
    mf = scf.ROHF(mol); mf.conv_tol = 1e-8
    mfn = scf.newton(mf); mfn.kernel()
    if not mfn.converged:
        m2 = scf.ROHF(mol); m2.level_shift = 0.4; m2.max_cycle = 400; m2.kernel()
        mfn = scf.newton(m2); mfn.kernel(m2.make_rdm1())
    return mfn


def cas_nevpt2(atoms, spin, avas_labels):
    """ROHF→AVAS→CASSCF→SC-NEVPT2 на фиксированной геометрии."""
    mol = build(atoms, spin); mol.max_memory = 12000
    mf = rohf(mol)
    ncas, nelec, mo = avas.avas(mf, avas_labels, threshold=0.2)
    ss = spin / 2 * (spin / 2 + 1)
    mc0 = mcscf.CASCI(mf, ncas, nelec)
    fci.addons.fix_spin_(mc0.fcisolver, shift=0.2, ss=ss)
    mc0.kernel(mo)
    mc = mcscf.CASSCF(mf, ncas, nelec)
    fci.addons.fix_spin_(mc.fcisolver, shift=0.2, ss=ss)
    mc.max_cycle_macro = 150; mc.conv_tol = 1e-7
    mc.kernel(mc0.cas_natorb()[0])
    e_pt = NEVPT(mc).kernel()
    return {"cas": [int(np.sum(nelec)), int(ncas)],
            "e_casscf_h": float(mc.e_tot), "casscf_converged": bool(mc.converged),
            "e_nevpt2_h": float(mc.e_tot + e_pt)}


def stage_dft(PAIR):
    M, Sol, dMM, dMC = PANEL[PAIR]
    path = os.path.join(DIR, f"pyrolysis{SUF}_{PAIR}.json")
    out = json.load(open(path)) if os.path.exists(path) else {}
    out.update({"pair": PAIR, "active": M, "solvent": Sol,
                "model": "hetero-bimetallic dimer stand-in for molten alloy"})
    t0 = time.time()
    dm_a, dm_e, dm_s, dm_c, dm_s2 = relax_lowest(dimer(M, Sol, dMM), tag=f"{PAIR}_dimer")
    dc_a, dc_e, dc_s, dc_c, dc_s2 = relax_lowest(dimer_C(M, Sol, dMM, dMC),
                                                 tag=f"{PAIR}_dimerC")
    ch_a, ch_e, ch_s, ch_c, ch_s2 = relax_lowest(dimer_CH3_H(M, Sol, dMM, dMC),
                                                  tag=f"{PAIR}_CH3H")
    ec = c_atom_energy("dft")
    ech4 = ch4_energy("dft")
    ebind = (dc_e - dm_e - ec) * HARTREE_KCAL          # <0 = связь C (коксование)
    ediss = (ch_e - dm_e - ech4) * HARTREE_KCAL         # <0 = активация CH4 выгодна
    out.update({
        "dimer_atoms": [[s, [round(x, 5) for x in c]] for s, c in dm_a],
        "dimerC_atoms": [[s, [round(x, 5) for x in c]] for s, c in dc_a],
        "ch3h_atoms": [[s, [round(x, 5) for x in c]] for s, c in ch_a],
        "spin_dimer": dm_s, "spin_dimerC": dc_s, "spin_ch3h": ch_s,
        "e_dimer_dft_h": dm_e, "e_dimerC_dft_h": dc_e, "e_ch3h_dft_h": ch_e,
        "e_C_dft_h": ec, "e_ch4_dft_h": ech4,
        "conv": {"dimer": dm_c, "dimerC": dc_c, "ch3h": ch_c},
        "s2": {"dimer": round(dm_s2, 2), "dimerC": round(dc_s2, 2)},
        "Ebind_C_dft_kcal": round(ebind, 2),
        "dEdiss_CH4_dft_kcal": round(ediss, 2),
        "wall_s": round(time.time() - t0, 1)})
    json.dump(out, open(path, "w"), indent=1)
    print(f"[{PAIR}] DFT  E_bind(C) = {ebind:.1f} | dE_diss(CH4) = {ediss:.1f} "
          f"kcal/mol (spins d/dC/ch3h = {dm_s}/{dc_s}/{ch_s})", flush=True)


def stage_corr(PAIR):
    M, Sol, dMM, dMC = PANEL[PAIR]
    path = os.path.join(DIR, f"pyrolysis{SUF}_{PAIR}.json")
    out = json.load(open(path))
    cfg_a = ACTIVE[M]; cfg_s = SOLVENT[Sol]
    t0 = time.time()
    # коррелированная E_bind(C): димер (M nd + Sol вал) и димер+C (M nd + Sol вал + C 2p)
    dm = cas_nevpt2([(s, tuple(c)) for s, c in out["dimer_atoms"]],
                    out["spin_dimer"], [cfg_a["orb"], cfg_s["orb"]])
    dc = cas_nevpt2([(s, tuple(c)) for s, c in out["dimerC_atoms"]],
                    out["spin_dimerC"], [cfg_a["orb"], cfg_s["orb"], "C 2p"])
    ec = c_atom_energy("corr")
    for lvl in ("casscf", "nevpt2"):
        eb = (dc[f"e_{lvl}_h"] - dm[f"e_{lvl}_h"] - ec) * HARTREE_KCAL
        out[f"Ebind_C_{lvl}_kcal"] = round(eb, 2)
    out["dimer_corr"] = dm; out["dimerC_corr"] = dc; out["e_C_corr_h"] = ec
    out["corr_wall_s"] = round(time.time() - t0, 1)
    json.dump(out, open(path, "w"), indent=1)
    print(f"[{PAIR}] E_bind(C) kcal/mol: DFT {out.get('Ebind_C_dft_kcal')} | "
          f"CASSCF {out['Ebind_C_casscf_kcal']} | NEVPT2 {out['Ebind_C_nevpt2_kcal']}",
          flush=True)


def stage_merge():
    rows = []
    for PAIR in PANEL:
        p = os.path.join(DIR, f"pyrolysis{SUF}_{PAIR}.json")
        if os.path.exists(p):
            d = json.load(open(p))
            rows.append({"pair": PAIR,
                         "Ebind_C_dft": d.get("Ebind_C_dft_kcal"),
                         "Ebind_C_casscf": d.get("Ebind_C_casscf_kcal"),
                         "Ebind_C_nevpt2": d.get("Ebind_C_nevpt2_kcal"),
                         "dEdiss_CH4_dft": d.get("dEdiss_CH4_dft_kcal")})
    res = {"descriptor": "molten-metal pyrolysis: C-binding volcano E_bind(C) on "
                         "active-solvent dimer + CH4 dissociative activation dE_diss",
           "note": "minimal hetero-bimetallic dimer stand-in (NOT periodic melt); "
                   "screen = relative volcano. Strong C-binding => carbide/coking; "
                   "weak => no CH4 activation. NEVPT2 column is the multireference-"
                   "corrected C-binding (metal-C bond de-risk). Ni-Bi = working ref.",
           "pairs": rows}
    json.dump(res, open(os.path.join(DIR, f"pyrolysis_descriptor{SUF}_results.json"), "w"),
              indent=1)
    print(json.dumps(res, ensure_ascii=False, indent=1))


if __name__ == "__main__":
    stage = sys.argv[1]
    if stage == "dft":
        stage_dft(sys.argv[2])
    elif stage == "corr":
        stage_corr(sys.argv[2])
    elif stage == "merge":
        stage_merge()
    else:
        raise SystemExit("usage: dft|corr <PAIR>, or merge")
