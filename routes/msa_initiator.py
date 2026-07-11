#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
routes/msa_initiator.py — дизайн инициатора MSA-цепи (рычаг конверсии и
строки $100/т в msa_cost_model). Хороший инициатор: (а) O–O гомолиз мягкий
(рождает радикалы при ~50°C, BDE ~25-35 ккал уже норм при каталитических
загрузках), (б) его радикал рвёт C–H метана НЕ ХУЖЕ пропагирующего CH3SO3•
(легаси/наш DFT ~19 ккал) — тогда цепь стартует без потерь.

Считаем ДВА числа на кандидата, DFT(UKS-PBE0) vs ROHF->AVAS->CASSCF->NEVPT2:
  1) BDE(O–O): E(2 радикала) − E(пероксид)          [без TS — гомолиз]
  2) барьер HAT: X• + CH4 -> XH + CH3• (пин-скан с якорем и гейтами
     интерьерности/выбросов — дисциплина msa_chain v5)

Кандидаты:
  H2O2   -> 2 •OH        (дёшев; •OH сверхагрессивен — риск неселективности)
  H2S2O8 -> 2 •OSO3H     (Marshall/персульфат в олеуме — Grillo-класс)
  MSA-пероксид (CH3SO2)O-O(SO2CH3)? — упрощённый прокси: CH3SO3• уже посчитан
     в msa_chain (пропагация ~19) — используется как реперная линия.

Запуск: python3 routes/msa_initiator.py bde|hat|merge|sanity
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
BASIS = os.environ.get("MSAI_BASIS", "def2-svp")
XC = "pbe0"


def M(atoms, spin, charge=0):
    return gto.M(atom=atoms, basis=BASIS, charge=charge, spin=spin,
                 verbose=0, max_memory=16000)


def uks(mol, dm0=None):
    mf = dft.UKS(mol)
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


def relax(atoms, spin, maxsteps=120):
    from pyscf.geomopt.geometric_solver import optimize
    mf = uks(M(atoms, spin))
    mol = optimize(mf, maxsteps=maxsteps, assert_convergence=False)
    na = [(mol.atom_symbol(i), tuple(c))
          for i, c in enumerate(mol.atom_coords(unit="Angstrom"))]
    return na, uks(mol, dm0=mf.make_rdm1())


def nev(atoms, spin, labels, thr=0.5):
    mol = M(atoms, spin)
    mol.max_memory = 24000
    mf = rohf(mol)
    ncas, nelec, mo = avas.avas(mf, labels, threshold=thr)
    ss = spin / 2 * (spin / 2 + 1)
    mc = mcscf.CASSCF(mf, ncas, nelec)
    fci.addons.fix_spin_(mc.fcisolver, shift=0.2, ss=ss)
    mc.max_cycle_macro = 150
    mc.conv_tol = 1e-7
    mc.kernel(mo)
    return {"e_casscf": float(mc.e_tot),
            "e_nevpt2": float(mc.e_tot + NEVPT(mc).kernel()),
            "cas": [int(np.sum(nelec)), int(ncas)], "conv": bool(mc.converged)}


# ---------------------------------------------------------------- геометрии
def h2o2_atoms():
    return [("O", (0.0, 0.73, -0.05)), ("O", (0.0, -0.73, 0.05)),
            ("H", (0.83, 0.88, 0.44)), ("H", (-0.83, -0.88, 0.44))]


def oh_atoms():
    return [("O", (0.0, 0.0, 0.0)), ("H", (0.0, 0.0, 0.97))]


def hso4_atoms():
    """•OSO3H: S в центре, радикальный O сверху (+z), OH вниз."""
    return [("S", (0.0, 0.0, 0.0)), ("O", (0.0, 0.0, 1.58)),
            ("O", (1.32, 0.55, -0.55)), ("O", (-1.32, 0.55, -0.55)),
            ("O", (0.0, -1.42, -0.55)), ("H", (0.0, -1.75, -1.44))]


def h2s2o8_atoms():
    """HO3S–O–O–SO3H: два HSO4 через O–O мост (r~1.45)."""
    a = []
    for sgn in (1, -1):
        x0 = sgn * 2.25
        a += [("S", (x0, 0.0, 0.0)),
              ("O", (sgn * 0.72, 0.0, 0.35)),            # мостовой O
              ("O", (x0 + sgn * 0.6, 1.30, -0.55)),
              ("O", (x0 + sgn * 0.6, -1.30, -0.55)),
              ("O", (x0 + sgn * 1.35, 0.0, 1.05)),
              ("H", (x0 + sgn * 2.05, 0.6, 1.35))]
    return a


def ch4_atoms():
    d, th = 1.09, math.radians(109.47)
    a = [("C", (0, 0, 0)), ("H", (0, 0, d))]
    for k in range(3):
        phi = math.radians(120 * k)
        a.append(("H", (d * math.sin(th) * math.cos(phi),
                        d * math.sin(th) * math.sin(phi), d * math.cos(th))))
    return a


def hat_ts_oh(rCH, rOH):
    """[H3C···H···OH]: C=0, H_t=1, O=2, H(O)=3 + 3H метила."""
    a = [("C", (0.0, 0.0, 0.0)), ("H", (0.18, 0.0, rCH)),
         ("O", (0.36, 0.0, rCH + rOH)), ("H", (1.22, 0.0, rCH + rOH + 0.30))]
    d, th = 1.09, math.radians(107.0)
    for k in range(3):
        phi = math.radians(120 * k + 60)
        a.append(("H", (d * math.sin(th) * math.cos(phi),
                        d * math.sin(th) * math.sin(phi),
                        -d * abs(math.cos(th)) + 0.13)))
    return a


def hat_ts_hso4(rCH, rOH):
    """[H3C···H···O–SO3H]: C=0, H_t=1, радикальный O=2, дальше SO3H-хвост."""
    zO = rCH + rOH
    a = [("C", (0.0, 0.0, 0.0)), ("H", (0.18, 0.0, rCH)),
         ("O", (0.36, 0.0, zO)), ("S", (0.36, 0.0, zO + 1.60)),
         ("O", (1.60, 0.60, zO + 2.10)), ("O", (-0.88, 0.60, zO + 2.10)),
         ("O", (0.36, -1.42, zO + 2.15)), ("H", (0.36, -1.75, zO + 3.04))]
    d, th = 1.09, math.radians(107.0)
    for k in range(3):
        phi = math.radians(120 * k + 60)
        a.append(("H", (d * math.sin(th) * math.cos(phi),
                        d * math.sin(th) * math.sin(phi),
                        -d * abs(math.cos(th)) + 0.13)))
    return a


CAND = {
    "OH": dict(perox=("H2O2", h2o2_atoms, 0), rad=("OH", oh_atoms, 1),
               hat_ts=hat_ts_oh, rad_labels=["0 O 2p"],
               ts_labels=["0 C 2p", "1 H 1s", "2 O 2p"]),
    "HSO4": dict(perox=("H2S2O8", h2s2o8_atoms, 0), rad=("HSO4", hso4_atoms, 1),
                 hat_ts=hat_ts_hso4,
                 rad_labels=["0 S 3p", "1 O 2p", "2 O 2p", "3 O 2p"],
                 ts_labels=["0 C 2p", "1 H 1s", "2 O 2p"]),
}
SCAN = ((1.25, 1.15), (1.35, 1.05), (1.50, 0.98),
        (1.18, 1.24), (1.12, 1.36), (1.10, 1.48))   # якорь+поздние, потом ранние


def stage_bde():
    out = {}
    path = os.path.join(DIR, "msa_initiator_bde.json")
    if os.path.exists(path):
        out = json.load(open(path))
    for tag, c in CAND.items():
        if tag in out and "bde_nevpt2" in out[tag]:
            print(f"[bde][{tag}] resume", flush=True)
            continue
        t0 = time.time()
        pn, pa, ps = c["perox"]
        rn, ra, rs = c["rad"]
        perox, mfp = relax(pa(), ps)
        rad, mfr = relax(ra(), rs)
        bde_dft = (2 * mfr.e_tot - mfp.e_tot) * HARTREE_KCAL
        # NEVPT2: пероксид — синглет с O-O sigma/sigma* в CAS; радикал — SOMO+O2p
        np_ = nev(perox, ps, ["O 2p"], thr=0.55)
        nr = nev(rad, rs, c["rad_labels"], thr=0.5)
        rec = {"bde_dft": round(bde_dft, 2),
               "bde_casscf": round((2 * nr["e_casscf"] - np_["e_casscf"])
                                   * HARTREE_KCAL, 2),
               "bde_nevpt2": round((2 * nr["e_nevpt2"] - np_["e_nevpt2"])
                                   * HARTREE_KCAL, 2),
               "cas_perox": np_["cas"], "cas_rad": nr["cas"],
               "conv": bool(np_["conv"] and nr["conv"]),
               "rad_atoms": [[s, [round(x, 6) for x in xyz]] for s, xyz in rad],
               "wall_s": round(time.time() - t0, 1)}
        out[tag] = rec
        with open(path, "w") as f:
            json.dump(out, f, indent=1)
        print(f"[bde][{tag}] DFT {rec['bde_dft']} | NEVPT2 {rec['bde_nevpt2']} "
              f"({rec['wall_s']}s)", flush=True)


def stage_hat():
    from pyscf.geomopt.geometric_solver import optimize
    bde = json.load(open(os.path.join(DIR, "msa_initiator_bde.json")))
    path = os.path.join(DIR, "msa_initiator_hat.json")
    out = {}
    if os.path.exists(path):
        out = json.load(open(path))
    ch4, m_ch4 = relax(ch4_atoms(), 0)
    for tag, c in CAND.items():
        if tag in out and ("nevpt2" in out[tag] or "failed" in out[tag]):
            print(f"[hat][{tag}] resume", flush=True)
            continue
        t0 = time.time()
        rad = [(s, tuple(xyz)) for s, xyz in bde[tag]["rad_atoms"]]
        m_rad = uks(M(rad, 1))
        eR = m_rad.e_tot + m_ch4.e_tot
        pts, dm, start_override = [], None, None
        for i, (rCH, rOH) in enumerate(SCAN):
            cf = os.path.join(DIR, f"_msai_{tag}_{rCH:.2f}.txt")
            open(cf, "w").write(
                f"$set\ndistance 1 2 {rCH:.4f}\ndistance 2 3 {rOH:.4f}\n")
            try:
                start = start_override if (i >= 3 and start_override) \
                    else c["hat_ts"](rCH, rOH)
                mf0 = uks(M(start, 1), dm0=dm)
                mol = optimize(mf0, constraints=cf, maxsteps=180,
                               assert_convergence=False)
                mfx = uks(mol, dm0=mf0.make_rdm1())
                na = [(mol.atom_symbol(j), tuple(xyz)) for j, xyz in
                      enumerate(mol.atom_coords(unit="Angstrom"))]
                dm = mfx.make_rdm1()
                pts.append({"rCH": rCH, "rOH": rOH, "e_h": float(mfx.e_tot),
                            "conv": bool(mfx.converged), "atoms": na})
                if i == 0:
                    start_override = na       # якорь для ранней ветки
                print(f"  [hat][{tag}] rCH={rCH} rel="
                      f"{(mfx.e_tot-eR)*HARTREE_KCAL:.2f} conv={mfx.converged}",
                      flush=True)
            except Exception as ex:
                print(f"  [hat][{tag}] rCH={rCH} failed: {ex}", flush=True)
            finally:
                if os.path.exists(cf):
                    os.remove(cf)
        pts.sort(key=lambda p: p["rCH"])
        rel = [(p["e_h"] - eR) * HARTREE_KCAL for p in pts]
        # маска выбросов + интерьерность (дисциплина msa v5)
        keep = [i for i, p in enumerate(pts) if p["conv"] and not all(
            rel[i] - rel[j] > 25.0 for j in (i - 1, i + 1)
            if 0 <= j < len(pts))]
        rec = {"scan_rel_kcal": [round(x, 2) for x in rel],
               "kept": keep, "wall_s": round(time.time() - t0, 1)}
        kr = [rel[i] for i in keep]
        if len(kr) >= 4:
            imax = max(range(1, len(kr) - 1), key=lambda i: kr[i])
            interior = kr[imax] >= max(kr[0], kr[-1])
            if interior and 0 < kr[imax] < 60:
                ts = pts[keep[imax]]
                rec["dft_barrier"] = round(kr[imax], 2)
                n_ts = nev(ts["atoms"], 1, c["ts_labels"], thr=0.4)
                n_r = nev([(s, tuple(x)) for s, x in bde[tag]["rad_atoms"]],
                          1, c["rad_labels"], thr=0.5)
                n_m = nev(ch4, 0, ["0 C 2p", "1 H 1s"], thr=0.6)
                for lv in ("casscf", "nevpt2"):
                    rec[lv] = round((n_ts[f"e_{lv}"] - n_r[f"e_{lv}"]
                                     - n_m[f"e_{lv}"]) * HARTREE_KCAL, 2)
            else:
                rec["failed"] = f"TS not interior/unphysical ({kr})"
        else:
            rec["failed"] = "too few clean points"
        out[tag] = rec
        with open(path, "w") as f:
            json.dump(out, f, indent=1)
        print(f"[hat][{tag}] -> {rec.get('dft_barrier', rec.get('failed'))} "
              f"nevpt2={rec.get('nevpt2')}", flush=True)


def stage_merge():
    bde = json.load(open(os.path.join(DIR, "msa_initiator_bde.json")))
    hat = json.load(open(os.path.join(DIR, "msa_initiator_hat.json")))
    res = {"screen": "MSA initiator design: O-O BDE + HAT(X*+CH4), DFT vs NEVPT2",
           "propagation_reference_kcal": 19.0,
           "candidates": {}}
    for tag in CAND:
        b, h = bde.get(tag, {}), hat.get(tag, {})
        res["candidates"][tag] = {
            "bde_OO": {k: b.get(f"bde_{k}") for k in ("dft", "casscf", "nevpt2")},
            "hat_CH4": {"dft": h.get("dft_barrier"),
                        "casscf": h.get("casscf"), "nevpt2": h.get("nevpt2"),
                        "failed": h.get("failed")}}
    json.dump(res, open(os.path.join(DIR, "msa_initiator_results.json"), "w"),
              indent=1)
    print(json.dumps(res, indent=1, ensure_ascii=False))


def stage_sanity():
    for tag, c in CAND.items():
        for name, fn, sp in (c["perox"], c["rad"]):
            pass
        pn, pa, ps = c["perox"]
        rn, ra, rs = c["rad"]
        for name, atoms, spin in ((pn, pa(), ps), (rn, ra(), rs)):
            try:
                mol = M(atoms, spin)
                print(f"[sanity] {name}: nao={mol.nao} nelec={mol.nelectron} OK",
                      flush=True)
            except Exception as e:
                print(f"[sanity] {name} FAIL: {e}", flush=True)


if __name__ == "__main__":
    st = sys.argv[1] if len(sys.argv) > 1 else "sanity"
    {"bde": stage_bde, "hat": stage_hat, "merge": stage_merge,
     "sanity": stage_sanity}[st]()
