#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
routes/chemloop_vacancy.py — дескриптор химического loop (маршрут 2 разведки):
энергия образования кислородной ВАКАНСИИ (отдачи решёточного O) на панели
редокс-металлов, через CASSCF/NEVPT2 — там, где DFT+U заведомо хрупок.

Идея (из METHANE_ROUTES_RESEARCH.md). В chemical-looping носитель отдаёт
решёточный O метану, потом регенерируется воздухом. Селективность управляется
ОДНИМ дескриптором — энергией образования кислородной вакансии Ev (и
нуклеофильностью решёточного O). Есть вулкан Сабатье: слишком лёгкая вакансия
переокисляет (combustion → CO2), слишком тяжёлая не активирует CH4. Нужен
ОПТИМУМ, а его положение чувствительно к точности коррелированной энергетики —
прямой повод для CASSCF/NEVPT2 на редокс-3d-центрах.

Минимальная модель (строго говоря — стенд-ин, не периодический оксид, как и Этап 15):
терминальный оксо-центр M(IV)=O с двумя OH-лигандами, [O=M(OH)2] (нейтр.).
Реакция отдачи решёточного O (2e редокс-пара M(IV)→M(II)):
    [O=M(OH)2]  →  [M(OH)2]  +  1/2 O2      Ev = E(red) + 1/2 E(O2) − E(oxo)
Скрин Ev по M = построение вулкана; наш Mn — репер, связывает с Этапом 15/16.

Дескриптор нуклеофильности решёточного O: заряд оксо-кислорода (Löwdin) и его
вклад 2p в активное пространство — печатается вместе с Ev.

Уровни: UKS-PBE0/def2-SVP (геометрия + DFT-Ev, все металлы, дёшево);
ROHF→AVAS(M 3d + oxo O 2p)→CASSCF→SC-NEVPT2 (коррелированный Ev — на AWS).
Fan-out: один металл на запуск → embarrassingly parallel на все ядра.

Запуск:
  python3 routes/chemloop_vacancy.py dft   <M>     # геометрии + DFT-Ev
  python3 routes/chemloop_vacancy.py corr  <M>     # + CASSCF/NEVPT2-Ev (AWS)
  python3 routes/chemloop_vacancy.py merge         # свести вулкан по всем M
"""
import json
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

# Панель редокс-металлов. Для каждого: 2S оксо-состояния M(IV)=O (d^n,
# высокоспиновое) и восстановленного M(II) (d^(n+2)). Ce — 4f, отдельная метка.
# (2S = число неспаренных e; высокоспиновый предел свободного иона M(IV)/M(II).)
METALS = {
    "Ti": {"Z": 22, "spin_ox": 0, "spin_red": 2, "orb": "Ti 3d"},
    "V":  {"Z": 23, "spin_ox": 1, "spin_red": 3, "orb": "V 3d"},
    "Cr": {"Z": 24, "spin_ox": 2, "spin_red": 4, "orb": "Cr 3d"},
    "Mn": {"Z": 25, "spin_ox": 3, "spin_red": 5, "orb": "Mn 3d"},
    "Fe": {"Z": 26, "spin_ox": 4, "spin_red": 4, "orb": "Fe 3d"},
    "Co": {"Z": 27, "spin_ox": 3, "spin_red": 3, "orb": "Co 3d"},
    "Ni": {"Z": 28, "spin_ox": 2, "spin_red": 2, "orb": "Ni 3d"},
    "Cu": {"Z": 29, "spin_ox": 1, "spin_red": 1, "orb": "Cu 3d"},
    "Ce": {"Z": 58, "spin_ox": 0, "spin_red": 2, "orb": "Ce 4f"},
    # расширение панели (Франкфурт): 4d/5d + главная подгруппа. W — компонент
    # нашего же Na2WO4/Mn катализатора OCM (связка с Этапом 16); Sn/Pb — классика
    # оксидных носителей POx; вместе питают карту I4 (O2--насос: окно напряжений).
    "Zr": {"Z": 40, "spin_ox": 0, "spin_red": 2, "orb": "Zr 4d"},
    "Nb": {"Z": 41, "spin_ox": 1, "spin_red": 3, "orb": "Nb 4d"},
    "Mo": {"Z": 42, "spin_ox": 2, "spin_red": 4, "orb": "Mo 4d"},
    "W":  {"Z": 74, "spin_ox": 2, "spin_red": 4, "orb": "W 5d"},
    "Sn": {"Z": 50, "spin_ox": 0, "spin_red": 0, "orb": "Sn 5p"},
    "Pb": {"Z": 82, "spin_ox": 0, "spin_red": 0, "orb": "Pb 6p"},
}
CAS_NO = 8          # целевое число активных орбиталей (M nd + oxo O 2p), фикс-состав


def oxo_atoms(M):
    """[O=M(OH)2]: oxo сверху по z, два OH в плоскости. Стартовые длины типовые."""
    return [(M, (0.0, 0.0, 0.0)), ("O", (0.0, 0.0, 1.60)),          # oxo
            ("O", (1.75, 0.0, -0.55)), ("H", (2.35, 0.60, -0.55)),
            ("O", (-1.75, 0.0, -0.55)), ("H", (-2.35, -0.60, -0.55))]


def red_atoms(M):
    """[M(OH)2] — тот же каркас без oxo (кислородная вакансия)."""
    return [(M, (0.0, 0.0, 0.0)),
            ("O", (1.75, 0.0, -0.55)), ("H", (2.35, 0.60, -0.55)),
            ("O", (-1.75, 0.0, -0.55)), ("H", (-2.35, -0.60, -0.55))]


def build(atoms, spin):
    return gto.M(atom=atoms, basis=BASIS, ecp=BASIS, charge=0, spin=spin,
                 verbose=0, max_memory=8000)


def uks(mol, dm0=None):
    mf = dft.UKS(mol).density_fit()
    mf.xc = XC
    mf.conv_tol = 1e-8
    mfn = scf.newton(mf)
    mfn.kernel(dm0=dm0) if dm0 is not None else mfn.kernel()
    if not mfn.converged:
        m2 = dft.UKS(mol).density_fit(); m2.xc = XC; m2.level_shift = 0.3
        m2.max_cycle = 300; m2.kernel(dm0=dm0)
        mfn = scf.newton(m2); mfn.kernel(m2.make_rdm1())
    return mfn


def relax(atoms, spin, tag):
    from pyscf.geomopt.geometric_solver import optimize
    mf = uks(build(atoms, spin))
    mol = optimize(mf, maxsteps=120, assert_convergence=False)
    na = [(mol.atom_symbol(i), tuple(c))
          for i, c in enumerate(mol.atom_coords(unit="Angstrom"))]
    mfx = uks(mol)
    return na, float(mfx.e_tot), bool(mfx.converged), float(mfx.spin_square()[0])


def o2_energy(level="dft"):
    """E(O2) триплет — общий репер для 1/2 O2."""
    mol = gto.M(atom=[("O", (0, 0, 0)), ("O", (0, 0, 1.208))],
                basis=BASIS, spin=2, verbose=0)
    mf = uks(mol)
    if level == "dft":
        return float(mf.e_tot)
    no, ne, mo = avas.avas(mf, ["O 2p"], threshold=0.2)
    mc = mcscf.CASSCF(mf, no, ne); mc.kernel(mo)
    return float(mc.e_tot + NEVPT(mc).kernel())


def rohf(mol):
    mf = scf.ROHF(mol); mf.conv_tol = 1e-8
    mfn = scf.newton(mf); mfn.kernel()
    if not mfn.converged:
        m2 = scf.ROHF(mol); m2.level_shift = 0.4; m2.max_cycle = 300; m2.kernel()
        mfn = scf.newton(m2); mfn.kernel(m2.make_rdm1())
    return mfn


def cas_nevpt2(atoms, spin, avas_labels):
    """ROHF→AVAS(M nd + oxo O 2p)→CASSCF→SC-NEVPT2 + заряд/нуклеофильность oxo."""
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
    # нуклеофильность oxo O (индекс 1 в oxo_atoms): Löwdin-заряд
    pop = mf.mulliken_pop(verbose=0)[1]
    q_oxo = float(pop[1]) if len(atoms) == 6 else None
    return {"cas": [int(np.sum(nelec)), int(ncas)],
            "e_casscf_h": float(mc.e_tot), "casscf_converged": bool(mc.converged),
            "e_nevpt2_h": float(mc.e_tot + e_pt), "q_oxo": q_oxo}


def stage_dft(M):
    cfg = METALS[M]
    out = {"metal": M, "level": "dft", "model": "[O=M(OH)2] minimal stand-in"}
    path = os.path.join(DIR, f"chemloop_{M}.json")
    if os.path.exists(path):
        out = json.load(open(path))
    t0 = time.time()
    ox_a, ox_e, ox_c, ox_s = relax(oxo_atoms(M), cfg["spin_ox"], f"{M}_ox")
    rd_a, rd_e, rd_c, rd_s = relax(red_atoms(M), cfg["spin_red"], f"{M}_red")
    eo2 = o2_energy("dft")
    ev = (rd_e + 0.5 * eo2 - ox_e) * HARTREE_KCAL
    out.update({"oxo_atoms": [[s, [round(x, 5) for x in c]] for s, c in ox_a],
                "red_atoms": [[s, [round(x, 5) for x in c]] for s, c in rd_a],
                "spin_ox": cfg["spin_ox"], "spin_red": cfg["spin_red"],
                "e_oxo_dft_h": ox_e, "e_red_dft_h": rd_e, "e_o2_dft_h": eo2,
                "ox_conv": ox_c, "red_conv": rd_c, "s2_ox": round(ox_s, 2),
                "Ev_dft_kcal": round(ev, 2), "wall_s": round(time.time() - t0, 1)})
    json.dump(out, open(path, "w"), indent=1)
    print(f"[{M}] DFT Ev = {ev:.1f} kcal/mol (oxo conv={ox_c} red conv={rd_c})",
          flush=True)


def stage_corr(M):
    path = os.path.join(DIR, f"chemloop_{M}.json")
    out = json.load(open(path))
    cfg = METALS[M]
    t0 = time.time()
    ox = cas_nevpt2([(s, tuple(c)) for s, c in out["oxo_atoms"]],
                    cfg["spin_ox"], [cfg["orb"], "1 O 2p"])   # M nd + oxo O 2p
    rd = cas_nevpt2([(s, tuple(c)) for s, c in out["red_atoms"]],
                    cfg["spin_red"], [cfg["orb"]])                # M nd (вакансия)
    eo2 = o2_energy("corr")
    for lvl in ("casscf", "nevpt2"):
        ev = (rd[f"e_{lvl}_h"] + 0.5 * eo2 - ox[f"e_{lvl}_h"]) * HARTREE_KCAL
        out[f"Ev_{lvl}_kcal"] = round(ev, 2)
    out["oxo_corr"] = ox; out["red_corr"] = rd; out["e_o2_corr_h"] = eo2
    out["corr_wall_s"] = round(time.time() - t0, 1)
    json.dump(out, open(path, "w"), indent=1)
    print(f"[{M}] Ev kcal/mol: DFT {out.get('Ev_dft_kcal')} | "
          f"CASSCF {out['Ev_casscf_kcal']} | NEVPT2 {out['Ev_nevpt2_kcal']} | "
          f"q_oxo {ox.get('q_oxo')}", flush=True)


def stage_merge():
    rows = []
    for M in METALS:
        p = os.path.join(DIR, f"chemloop_{M}.json")
        if os.path.exists(p):
            d = json.load(open(p))
            rows.append({"metal": M, "Ev_dft": d.get("Ev_dft_kcal"),
                         "Ev_casscf": d.get("Ev_casscf_kcal"),
                         "Ev_nevpt2": d.get("Ev_nevpt2_kcal"),
                         "q_oxo": d.get("oxo_corr", {}).get("q_oxo")})
    res = {"descriptor": "oxygen-vacancy formation energy Ev = E(red)+1/2E(O2)-E(oxo)",
           "note": "minimal [O=M(OH)2] stand-in clusters; screen = relative volcano, "
                   "not absolute. NEVPT2 column is the multireference-corrected Ev.",
           "carriers": rows}
    json.dump(res, open(os.path.join(DIR, "chemloop_vacancy_results.json"), "w"),
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
        raise SystemExit("usage: dft|corr <M>, or merge")
