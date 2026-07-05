#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
routes/acetic_radical.py — изобретение B / маршрут №7 карты: прямое сочетание
CH4 + CO2 → CH3COOH (уксусная кислота). Максимальная маржа: две БЕСПЛАТНЫЕ
отходные молекулы факела → жидкость ~$500/т в цистерну, один передел.
Ноль коммерциализаций — и мы утверждаем, что причина в DFT-слепоте ключевых TS.

Радикально-цепная схема (фаза 1, газофазный скелет без катализатора):
  (1) C–C сочетание:      CH3• + CO2  → CH3COO•   (карбоксилация радикала)
  (2) декарбоксилирование: CH3COO• → CH3• + CO2    (обратно; известный убийца —
      быстрый распад ацетоксила; σ/π-состояния RCO2• — классический
      мультиреференс-провал DFT, барьеры врут на порядки скорости)
  (3) пропагация цепи:     CH3COO• + CH4 → CH3COOH + CH3•  (HAT, замыкает цикл)
Если (3) конкурентоспособно с (2) — цепь живёт, и прямая уксусная возможна;
решают ТОЧНЫЕ барьеры (1)-(3). Считаем UKS-PBE0 vs CASSCF/NEVPT2 (AVAS: C 2p +
O 2p — π-система CO2/карбоксила + радикальная орбиталь).

Модель: газофазные молекулы, def2-SVP (скрин; TZVP — ручкой ACE_BASIS).
TS(1)/(2) — один и тот же сэдл C–C: пин r(C–C), релакс остального, скан-максимум.
TS(3) — коллинеарный [H3C···H···O(акцептор)] как в галоген/antibep скриптах.

Запуск:
  python3 routes/acetic_radical.py cc      # стадии (1)+(2): скан C–C + концы
  python3 routes/acetic_radical.py hat     # стадия (3): HAT CH4 -> ацетоксил
  python3 routes/acetic_radical.py merge   # сводка трёх чисел + вердикт цепи
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
BASIS = os.environ.get("ACE_BASIS", "def2-svp")
XC = "pbe0"


def M(atoms, spin, charge=0):
    return gto.M(atom=atoms, basis=BASIS, charge=charge, spin=spin,
                 verbose=0, max_memory=6000)


def uks(mol, dm0=None):
    mf = dft.UKS(mol); mf.xc = XC; mf.conv_tol = 1e-8; mf.max_cycle = 200
    mf.kernel(dm0=dm0) if dm0 is not None else mf.kernel()
    if not mf.converged:
        mf.level_shift = 0.3; mf.max_cycle = 300; mf.kernel(mf.make_rdm1())
    return mf


def rohf(mol):
    mf = scf.ROHF(mol); mf.conv_tol = 1e-8; mf.max_cycle = 200; mf.kernel()
    if not mf.converged:
        mf.level_shift = 0.4; mf.max_cycle = 300; mf.kernel(mf.make_rdm1())
    return mf


def relax(atoms, spin):
    from pyscf.geomopt.geometric_solver import optimize
    mf = uks(M(atoms, spin))
    mol = optimize(mf, maxsteps=100, assert_convergence=False)
    na = [(mol.atom_symbol(i), tuple(c))
          for i, c in enumerate(mol.atom_coords(unit="Angstrom"))]
    return na, uks(mol)


def nevpt2(atoms, spin, labels, thr=0.45):
    mol = M(atoms, spin)
    mf = rohf(mol)
    ncas, nelec, mo = avas.avas(mf, labels, threshold=thr)
    ss = spin / 2 * (spin / 2 + 1)
    mc = mcscf.CASSCF(mf, ncas, nelec)
    fci.addons.fix_spin_(mc.fcisolver, shift=0.2, ss=ss)
    mc.max_cycle_macro = 150; mc.conv_tol = 1e-7
    mc.kernel(mo)
    return {"e_casscf": float(mc.e_tot),
            "e_nevpt2": float(mc.e_tot + NEVPT(mc).kernel()),
            "cas": [int(np.sum(nelec)), int(ncas)], "conv": bool(mc.converged)}


# ---------------------------------------------------------------- геометрии
def ch3_atoms():
    d = 1.08
    return [("C", (0, 0, 0))] + [("H", (d*math.cos(math.radians(120*k)),
                                        d*math.sin(math.radians(120*k)), 0))
                                 for k in range(3)]


def ch4_atoms():
    d, th = 1.09, math.radians(109.47)
    a = [("C", (0, 0, 0)), ("H", (0, 0, d))]
    for k in range(3):
        phi = math.radians(120 * k)
        a.append(("H", (d*math.sin(th)*math.cos(phi), d*math.sin(th)*math.sin(phi),
                        d*math.cos(th))))
    return a


def co2_atoms():
    return [("C", (0, 0, 0)), ("O", (0, 0, 1.16)), ("O", (0, 0, -1.16))]


def ch3_co2_atoms(rcc):
    """CH3 (зонтиком, C в нуле) + CO2, C(CO2) на +z на расстоянии rcc, слегка изогнут."""
    d, th = 1.08, math.radians(100.0)
    a = [("C", (0, 0, 0))]
    for k in range(3):
        phi = math.radians(120 * k)
        a.append(("H", (d*math.sin(th)*math.cos(phi), d*math.sin(th)*math.sin(phi),
                        -d*math.cos(th))))
    # CO2 с изломом ~150° (готовится стать карбоксилом)
    a.append(("C", (0.0, 0.0, rcc)))
    a.append(("O", (1.05, 0.0, rcc + 0.55)))
    a.append(("O", (-1.05, 0.0, rcc + 0.55)))
    return a


def acetoxyl_atoms():
    """CH3COO• стартовая."""
    return [("C", (0, 0, 0)),
            ("H", (1.02, 0, -0.40)), ("H", (-0.51, 0.88, -0.40)),
            ("H", (-0.51, -0.88, -0.40)),
            ("C", (0, 0, 1.52)),
            ("O", (1.05, 0, 2.15)), ("O", (-1.05, 0, 2.15))]


def hat_ts_atoms(rCH, rOH):
    """[H3C···H···O–C(O)CH3]: ось z, метан-хвост вниз, ацетоксил вверх (излом от оси)."""
    a = [("C", (0, 0, 0)), ("H", (0.20, 0, rCH))]           # Cα, переносимый H
    d, th = 1.09, math.radians(105.0)
    for k in range(3):
        phi = math.radians(120 * k)
        a.append(("H", (d*math.sin(th)*math.cos(phi), d*math.sin(th)*math.sin(phi),
                        -d*math.cos(th))))
    zO = rCH + rOH
    a += [("O", (0.45, 0.0, zO)), ("C", (0.30, 0.0, zO + 1.32)),
          ("O", (1.30, 0.0, zO + 1.95)), ("C", (-1.05, 0.0, zO + 1.85)),
          ("H", (-1.02, 0.0, zO + 2.93)), ("H", (-1.60, 0.88, zO + 1.52)),
          ("H", (-1.60, -0.88, zO + 1.52))]
    return a


def cc_scan():
    """Пин r(C–C), релакс остального: сэдл сочетания CH3+CO2 (он же декарбоксилир.)."""
    from pyscf.geomopt.geometric_solver import optimize
    best = None
    prof = []
    for rcc in (2.30, 2.10, 1.95, 1.80, 1.65):
        cf = os.path.join(DIR, f"_ace_cc_{rcc:.2f}.txt")
        open(cf, "w").write(f"$set\ndistance 1 5 {rcc:.4f}\n")
        try:
            mf = uks(M(ch3_co2_atoms(rcc), 1))
            mol = optimize(mf, constraints=cf, maxsteps=70, assert_convergence=False)
            mfx = uks(mol)
            e = float(mfx.e_tot)
            na = [(mol.atom_symbol(i), tuple(c))
                  for i, c in enumerate(mol.atom_coords(unit="Angstrom"))]
            prof.append({"rcc": rcc, "e_h": e})
            if best is None or e > best[1]:
                best = (na, e, rcc)
        except Exception as ex:
            print(f"  [cc] r={rcc} failed: {ex}", flush=True)
        finally:
            if os.path.exists(cf):
                os.remove(cf)
    if best is None:
        raise RuntimeError("cc scan failed")
    return best, prof


def stage_cc():
    t0 = time.time()
    out = {"stage": "cc", "model": f"gas-phase radical, {BASIS}"}
    ch3, mf_ch3 = relax(ch3_atoms(), 1)
    co2, mf_co2 = relax(co2_atoms(), 0)
    ace, mf_ace = relax(acetoxyl_atoms(), 1)
    (ts_a, e_ts, rcc), prof = cc_scan()
    eR = mf_ch3.e_tot + mf_co2.e_tot
    out["dft"] = {
        "barrier_cc_kcal": round((e_ts - eR) * HARTREE_KCAL, 2),
        "dE_reac_kcal": round((mf_ace.e_tot - eR) * HARTREE_KCAL, 2),
        "barrier_decarb_kcal": round((e_ts - mf_ace.e_tot) * HARTREE_KCAL, 2),
        "ts_rcc": rcc, "profile": prof}
    lab = ["C 2p", "O 2p"]
    n_ts = nevpt2(ts_a, 1, lab)
    n_ch3 = nevpt2(ch3, 1, ["C 2p"])
    n_co2 = nevpt2(co2, 0, lab)
    n_ace = nevpt2(ace, 1, lab)
    for lvl in ("casscf", "nevpt2"):
        eRc = n_ch3[f"e_{lvl}"] + n_co2[f"e_{lvl}"]
        out[lvl] = {
            "barrier_cc_kcal": round((n_ts[f"e_{lvl}"] - eRc) * HARTREE_KCAL, 2),
            "dE_reac_kcal": round((n_ace[f"e_{lvl}"] - eRc) * HARTREE_KCAL, 2),
            "barrier_decarb_kcal": round((n_ts[f"e_{lvl}"] - n_ace[f"e_{lvl}"])
                                         * HARTREE_KCAL, 2)}
    out["ts_cas"] = n_ts["cas"]; out["wall_s"] = round(time.time() - t0, 1)
    json.dump(out, open(os.path.join(DIR, "acetic_cc.json"), "w"), indent=1)
    print(f"[cc] barrier CC/decarb kcal: DFT {out['dft']['barrier_cc_kcal']}/"
          f"{out['dft']['barrier_decarb_kcal']} | NEVPT2 "
          f"{out['nevpt2']['barrier_cc_kcal']}/{out['nevpt2']['barrier_decarb_kcal']}"
          f" ({out['wall_s']}s)", flush=True)


def stage_hat():
    """HAT: CH3COO• + CH4 → CH3COOH + CH3• — жёсткий скан (паттерн antibep)."""
    t0 = time.time()
    out = {"stage": "hat", "model": f"rigid collinear HAT scan, {BASIS}"}
    ch4, mf_ch4 = relax(ch4_atoms(), 0)
    ace, mf_ace = relax(acetoxyl_atoms(), 1)
    best = None
    for rCH, rOH in ((1.15, 1.35), (1.25, 1.20), (1.35, 1.10), (1.50, 1.00)):
        try:
            mfx = uks(M(hat_ts_atoms(rCH, rOH), 1))
            e = float(mfx.e_tot)
            if not np.isfinite(e):
                continue
        except Exception as ex:
            print(f"  [hat] rCH={rCH} failed: {ex}", flush=True)
            continue
        if best is None or e > best[1]:
            best = (hat_ts_atoms(rCH, rOH), e, rCH, rOH)
    if best is None:
        raise RuntimeError("hat scan failed")
    ts_a, e_ts, rCH, rOH = best
    eR = mf_ace.e_tot + mf_ch4.e_tot
    out["dft"] = {"barrier_hat_kcal": round((e_ts - eR) * HARTREE_KCAL, 2),
                  "ts_rCH": rCH, "ts_rOH": rOH}
    lab = ["C 2p", "O 2p", "H 1s"]
    n_ts = nevpt2(ts_a, 1, lab, thr=0.4)
    n_ch4 = nevpt2(ch4, 0, ["C 2p", "H 1s"], thr=0.6)
    n_ace = nevpt2(ace, 1, ["C 2p", "O 2p"])
    for lvl in ("casscf", "nevpt2"):
        eRc = n_ace[f"e_{lvl}"] + n_ch4[f"e_{lvl}"]
        out[lvl] = {"barrier_hat_kcal":
                    round((n_ts[f"e_{lvl}"] - eRc) * HARTREE_KCAL, 2)}
    out["ts_cas"] = n_ts["cas"]; out["wall_s"] = round(time.time() - t0, 1)
    json.dump(out, open(os.path.join(DIR, "acetic_hat.json"), "w"), indent=1)
    print(f"[hat] barrier kcal: DFT {out['dft']['barrier_hat_kcal']} | NEVPT2 "
          f"{out['nevpt2']['barrier_hat_kcal']} ({out['wall_s']}s)", flush=True)


def stage_merge():
    cc = json.load(open(os.path.join(DIR, "acetic_cc.json")))
    hat = json.load(open(os.path.join(DIR, "acetic_hat.json")))
    chain_alive = (hat["nevpt2"]["barrier_hat_kcal"]
                   <= cc["nevpt2"]["barrier_decarb_kcal"] + 3.0)
    res = {"route": "direct CH4 + CO2 -> CH3COOH (radical chain feasibility, phase 1)",
           "numbers_nevpt2": {
               "cc_coupling_barrier": cc["nevpt2"]["barrier_cc_kcal"],
               "decarboxylation_barrier": cc["nevpt2"]["barrier_decarb_kcal"],
               "hat_propagation_barrier": hat["nevpt2"]["barrier_hat_kcal"],
               "cc_dE_reac": cc["nevpt2"]["dE_reac_kcal"]},
           "numbers_dft": {
               "cc_coupling_barrier": cc["dft"]["barrier_cc_kcal"],
               "decarboxylation_barrier": cc["dft"]["barrier_decarb_kcal"],
               "hat_propagation_barrier": hat["dft"]["barrier_hat_kcal"]},
           "chain_verdict": ("ALIVE: HAT propagation outruns (or matches) acetoxyl "
                             "decarboxylation — radical carboxylation chain viable"
                             if chain_alive else
                             "DEAD in free gas phase: acetoxyl decarboxylates faster "
                             "than it harvests H from CH4 — need stabilization of "
                             "RCO2 (complexation/field/surface) to slow step 2"),
           "note": "gas-phase skeleton; DFT-vs-NEVPT2 gap on RCO2* states is the "
                   "de-risk (DFT notoriously wrong on acetoxyl). Screen-grade."}
    json.dump(res, open(os.path.join(DIR, "acetic_results.json"), "w"), indent=1)
    print(json.dumps(res, ensure_ascii=False, indent=1))


if __name__ == "__main__":
    {"cc": stage_cc, "hat": stage_hat, "merge": stage_merge}[sys.argv[1]]()
