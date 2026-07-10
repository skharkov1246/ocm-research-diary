#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
routes/lao_cr_trimer.py — вариант 2 гребёнки (ЛАО/ПАО-скид): дескриптор
селективной ТРИмеризации этилена на Cr. Этилен со скида не возится (−104°C);
олигомеризация на месте в 1-гексен/ЛАО→ПАО ($2500–4000/т, жидкость) снимает
логистику. Квант-ниша: цикл Cr(I)/Cr(III) — метациклы со спин-кроссовером,
редокс-пара, где DFT исторически путается (известная проблема литературы
триммеризации Phillips/Sasol).

СМОК-СТАДИЯ (эта): жизнеспособность оснастки, без TS —
  1) spins: лестницы UKS-PBE0 для хромациклопентана и хромациклогептана
     [Cl2Cr(C4H8)] / [Cl2Cr(C6H12)] (модель-скелет, как bare-этапы OCM);
  2) int: релакс интермедиатов + CASSCF/NEVPT2 одноточки (AVAS [Cr 3d]) —
     сходимость и NOON-мультиреференсность (есть ли вообще наша ниша).
Дескриптор (следующая стадия, по смоку): ΔΔE‡ = барьер(внедрение C2H4 в
C7-цикл → рост/C8+) − барьер(β-H элиминирование → 1-гексен); >0 = селективно
к C6. DFT vs NEVPT2 — переворачивает ли корреляция C6/C8-чтение.

Запуск: python3 routes/lao_cr_trimer.py spins|int|sanity
Env: LAO_BASIS (def2-svp).
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
BASIS = os.environ.get("LAO_BASIS", "def2-svp")
XC = "pbe0"


def M(atoms, spin, charge=0):
    return gto.M(atom=atoms, basis=BASIS, ecp=BASIS, charge=charge, spin=spin,
                 verbose=0, max_memory=8000)


def uks(mol, dm0=None):
    mf = dft.UKS(mol).density_fit()
    mf.xc = XC
    mf.conv_tol = 1e-8
    mf.max_cycle = 200
    mf.kernel(dm0=dm0) if dm0 is not None else mf.kernel()
    if not mf.converged:
        mf.level_shift = 0.3
        mf.max_cycle = 300
        mf.kernel(mf.make_rdm1())
    return mf


def rohf(mol):
    mf = scf.ROHF(mol)
    mf.conv_tol = 1e-8
    mf.max_cycle = 200
    mf.kernel()
    if not mf.converged:
        mf.level_shift = 0.4
        mf.max_cycle = 300
        mf.kernel(mf.make_rdm1())
    return mf


# ---------------------------------------------------------------- геометрии
def _ring_carbons(n, r0=2.05, ccb=1.53):
    """Позиции C метацикла Cr(CH2)n: Cr в нуле, кольцо в плоскости xz дугой.
    Грубый старт — релакс доводит."""
    pts = []
    span = math.radians(100 + 12 * n)          # раствор дуги под размер кольца
    for k in range(n):
        th = -span / 2 + span * k / (n - 1)
        rad = r0 + (ccb - r0 * 0.35) * min(k, n - 1 - k) * 0.55
        pts.append((rad * math.sin(th), 0.0, rad * math.cos(th)))
    return pts


def chromacycle_atoms(n_c):
    """[Cl2Cr(CnH2n)] метацикл: Cr + 2 Cl (тетраэдрично сверху) + кольцо CH2."""
    a = [("Cr", (0.0, 0.0, 0.0)),
         ("Cl", (1.65, 1.45, -0.85)), ("Cl", (-1.65, 1.45, -0.85))]
    cs = _ring_carbons(n_c)
    for x, y, z in cs:
        a.append(("C", (x, y, z)))
    d = 1.09
    for i, (x, y, z) in enumerate(cs):
        # два H на каждый C, разнесены из плоскости кольца
        a.append(("H", (x * 1.12, y + d * 0.82, z * 1.12)))
        a.append(("H", (x * 1.12, y - d * 0.82, z * 1.12)))
    return a


def relax(atoms, spin, maxsteps=150):
    from pyscf.geomopt.geometric_solver import optimize
    mf = uks(M(atoms, spin))
    mol = optimize(mf, maxsteps=maxsteps, assert_convergence=False)
    na = [(mol.atom_symbol(i), tuple(c))
          for i, c in enumerate(mol.atom_coords(unit="Angstrom"))]
    mfx = uks(mol, dm0=mf.make_rdm1())
    return na, mfx


def nevpt2_point(atoms, spin, thr=0.5):
    mol = M(atoms, spin)
    mf = rohf(mol)
    ncas, nelec, mo = avas.avas(mf, ["Cr 3d"], threshold=thr)
    ss = spin / 2 * (spin / 2 + 1)
    mc = mcscf.CASSCF(mf, ncas, nelec)
    fci.addons.fix_spin_(mc.fcisolver, shift=0.2, ss=ss)
    mc.max_cycle_macro = 150
    mc.conv_tol = 1e-7
    mc.kernel(mo)
    dm1 = mc.fcisolver.make_rdm1(mc.ci, mc.ncas, mc.nelecas)
    noon = sorted(np.linalg.eigvalsh(dm1), reverse=True)
    return {"e_casscf": float(mc.e_tot),
            "e_nevpt2": float(mc.e_tot + NEVPT(mc).kernel()),
            "cas": [int(np.sum(nelec)), int(ncas)],
            "conv": bool(mc.converged),
            "noon": [round(float(x), 3) for x in noon]}


RINGS = {"c5": 4, "c7": 6}     # хромациклопентан (4 CH2), -гептан (6 CH2)


def stage_spins():
    """Лестницы спинов обоих циклов. Нейтральный Cl2Cr(кольцо) = формально
    Cr(IV) d2 (кольцо −2, 2 Cl −1) — чётный счёт: кандидаты по чётности."""
    out = {}
    path = os.path.join(DIR, "lao_cr_spins.json")
    for tag, n in RINGS.items():
        atoms = chromacycle_atoms(n)
        parity = sum(gto.charge(a[0]) for a in atoms) % 2
        cands = [s for s in range(parity, parity + 6, 2)][:3]
        out[tag] = {"ladder": [], "candidates_2S": cands}
        for s2 in cands:
            t0 = time.time()
            try:
                mf = uks(M(atoms, s2))
                out[tag]["ladder"].append(
                    {"spin_2S": s2, "e_h": float(mf.e_tot),
                     "converged": bool(mf.converged),
                     "spin_square": round(float(mf.spin_square()[0]), 3),
                     "wall_s": round(time.time() - t0, 1)})
            except Exception as e:
                out[tag]["ladder"].append({"spin_2S": s2, "error": str(e)[:150]})
            print(f"[spins][{tag}] 2S={s2}: {out[tag]['ladder'][-1]}", flush=True)
            with open(path, "w") as f:
                json.dump(out, f, indent=1)
        ok = [x for x in out[tag]["ladder"] if x.get("converged")]
        if ok:
            best = min(ok, key=lambda x: x["e_h"])
            out[tag]["ground_2S"] = best["spin_2S"]
            out[tag]["rel_kcal"] = {
                str(x["spin_2S"]): round((x["e_h"] - best["e_h"]) * HARTREE_KCAL, 2)
                for x in ok}
            print(f"[spins][{tag}] ground 2S={best['spin_2S']} "
                  f"ladder={out[tag]['rel_kcal']}", flush=True)
        with open(path, "w") as f:
            json.dump(out, f, indent=1)


def stage_int():
    """Релакс интермедиатов в основном спине + CASSCF/NEVPT2-одноточка:
    сходимость + NOON (мультиреференсность = есть ли ниша)."""
    sp = json.load(open(os.path.join(DIR, "lao_cr_spins.json")))
    path = os.path.join(DIR, "lao_cr_int.json")
    out = {}
    if os.path.exists(path):
        out = json.load(open(path))
    for tag, n in RINGS.items():
        if tag in out and out[tag].get("nevpt2"):
            print(f"[int][{tag}] resume: on disk", flush=True)
            continue
        s2 = sp[tag]["ground_2S"]
        t0 = time.time()
        atoms, mf = relax(chromacycle_atoms(n), s2)
        rec = {"spin_2S": s2, "e_dft": float(mf.e_tot),
               "dft_converged": bool(mf.converged),
               "atoms": [[s, [round(x, 6) for x in c]] for s, c in atoms]}
        try:
            rec.update({"nevpt2": nevpt2_point(atoms, s2)})
        except Exception as e:
            rec["nevpt2_error"] = str(e)[:200]
        rec["wall_s"] = round(time.time() - t0, 1)
        out[tag] = rec
        with open(path, "w") as f:
            json.dump(out, f, indent=1)
        print(f"[int][{tag}] 2S={s2} dft={rec['e_dft']:.6f} "
              f"nevpt2={rec.get('nevpt2')}", flush=True)
    # смок-вердикт: обе точки сошлись и NOON мультиреференсны?
    ver = {}
    for tag in RINGS:
        nv = out.get(tag, {}).get("nevpt2") or {}
        noon = nv.get("noon") or []
        frac = sum(1 for x in noon if 0.1 < x < 1.9)
        ver[tag] = {"casscf_conv": nv.get("conv"), "n_frac_noon": frac}
    out["_smoke_verdict"] = ver
    with open(path, "w") as f:
        json.dump(out, f, indent=1)
    print("[int] smoke verdict:", json.dumps(ver), flush=True)


def stage_sanity():
    for tag, n in RINGS.items():
        atoms = chromacycle_atoms(n)
        try:
            mol = M(atoms, sum(__import__("pyscf").gto.charge(a[0]) for a in atoms) % 2)
            print(f"[sanity] {tag}: natoms={len(atoms)} nao={mol.nao} "
                  f"nelec={mol.nelectron} build OK", flush=True)
        except Exception as e:
            print(f"[sanity] {tag} FAIL: {e}", flush=True)


if __name__ == "__main__":
    st = sys.argv[1] if len(sys.argv) > 1 else "sanity"
    {"spins": stage_spins, "int": stage_int, "sanity": stage_sanity}[st]()
