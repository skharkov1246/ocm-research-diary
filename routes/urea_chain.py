#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Вариант 3 (урея-синергия): газофазное ядро NH3 + CO2 → карбамат → урея.

Синергия скида: пиролиз CH4 → C + 2H2 (наш пиро-волкан), H2 → NH3, NH3 + CO2
→ урея (CO2 ПОТРЕБЛЯЕТСЯ). Считаемое ядро (скрин, PBE0/def2-SVP+DF):
  1) ts: согласованное присоединение NH3 + CO2 → NH2COOH — 2D-пин-скан
     (r(C–N) форм., r(O–H) перенос) по рецепту MSA v7 (поточечный чекпойнт,
     интерьерный TS, маска выбросов). Лит. якорь: газофазный барьер
     ~35–40 ккал/моль (без катализа второй молекулой).
  2) thermo: реакционные энергии NH3+CO2→NH2COOH и NH2COOH+NH3→урея+H2O.
  3) НАЗНАЧЕНИЕ ДЛЯ МЕТОДА M (отрицательный контроль): вся цепь —
     закрытые оболочки; NEVPT2 на R/TS должен СОВПАСТЬ с DFT (NOON ~2/0).
     Если совпал — граница применимости метода M задокументирована с двух
     сторон (радикалы/металлы: DFT врёт; закрытые оболочки: DFT ок).
Промышленный контекст: реальный синтез — жидкофазный карбамат аммония при
150 бар; газофаза здесь — скрин-якорь + контроль метода, не модель реактора.

Стадии: thermo | ts | merge.
Запуск: python3 -u routes/urea_chain.py <stage>
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
SCAN = os.path.join(DIR, "urea_ts_scan.json")
RES = os.path.join(DIR, "urea_results.json")


def M(atoms, spin=0):
    return gto.M(atom=[(s, tuple(c)) for s, c in atoms], basis="def2-svp",
                 charge=0, spin=spin, verbose=0, max_memory=8000)


def uks(mol, dm0=None):
    mf = dft.UKS(mol).density_fit()
    mf.xc = "pbe0"
    mf.conv_tol = 1e-8
    mfn = scf.newton(mf)
    mfn.kernel(dm0=dm0) if dm0 is not None else mfn.kernel()
    if not mfn.converged:
        mf2 = dft.UKS(mol).density_fit()
        mf2.xc = "pbe0"
        mf2.level_shift = 0.3
        mf2.max_cycle = 300
        mf2.kernel(dm0=dm0)
        mfn = scf.newton(mf2)
        mfn.kernel(mf2.make_rdm1())
    return mfn


def relax(atoms, maxsteps=150):
    from pyscf.geomopt.geometric_solver import optimize
    mf = uks(M(atoms))
    mol = optimize(mf, maxsteps=maxsteps, assert_convergence=False)
    mfx = uks(mol, dm0=mf.make_rdm1())
    na = [(mol.atom_symbol(i), tuple(c))
          for i, c in enumerate(mol.atom_coords(unit="Angstrom"))]
    return na, float(mfx.e_tot), bool(mfx.converged)


def nevpt2(atoms, labels, thr=0.5):
    mol = M(atoms)
    mol.max_memory = 12000
    mf = scf.RHF(mol)
    mf.conv_tol = 1e-8
    mfn = scf.newton(mf)
    mfn.kernel()
    ncas, nelec, mo = avas.avas(mfn, labels, threshold=thr)
    mc = mcscf.CASSCF(mfn, ncas, nelec)
    fci.addons.fix_spin_(mc.fcisolver, shift=0.2, ss=0.0)
    mc.max_cycle_macro = 150
    mc.conv_tol = 1e-7
    mc.kernel(mo)
    dm1 = mc.fcisolver.make_rdm1(mc.ci, mc.ncas, mc.nelecas)
    noon = sorted(np.linalg.eigvalsh(dm1), reverse=True)
    ne = (int(sum(nelec)) if hasattr(nelec, "__len__") else int(nelec))
    return {"cas": [ne, int(ncas)], "e_casscf": float(mc.e_tot),
            "e_nevpt2": float(mc.e_tot + NEVPT(mc).kernel()),
            "casscf_conv": bool(mc.converged),
            "noon": [round(float(x), 3) for x in noon],
            "n_frac_noon": int(sum(1 for x in noon if 0.1 < x < 1.9))}


# ------------------------------------------------------------- строители
def co2_atoms():
    return [("C", (0.0, 0.0, 0.0)), ("O", (1.16, 0.0, 0.0)),
            ("O", (-1.16, 0.0, 0.0))]


def nh3_atoms():
    return [("N", (0.0, 0.0, 0.0)), ("H", (0.94, 0.0, 0.38)),
            ("H", (-0.47, 0.81, 0.38)), ("H", (-0.47, -0.81, 0.38))]


def adduct_atoms(rCN):
    """CO2 в плоскости xz + NH3 снизу по z; один N–H нацелен на O1
    (перенос H). Порядок (1-based): C=1, O1=2, O2=3, N=4, H_t=5, H=6,7."""
    return [("C", (0.0, 0.0, 0.0)), ("O", (1.16, 0.0, 0.0)),
            ("O", (-1.16, 0.0, 0.0)), ("N", (0.0, 0.0, rCN)),
            ("H", (0.41, 0.0, rCN - 0.37 + 0.94)),
            ("H", (-0.48, 0.83, rCN + 0.35)),
            ("H", (-0.48, -0.83, rCN + 0.35))]


def carbamic_atoms():
    """NH2COOH стартовая (релакс доводит)."""
    return [("C", (0.0, 0.0, 0.0)), ("O", (1.20, 0.35, 0.0)),
            ("O", (-1.05, 0.75, 0.0)), ("N", (-0.25, -1.32, 0.0)),
            ("H", (-1.85, 0.28, 0.0)), ("H", (0.45, -1.95, 0.0)),
            ("H", (-1.20, -1.65, 0.0))]


def urea_atoms():
    return [("C", (0.0, 0.0, 0.0)), ("O", (0.0, 1.22, 0.0)),
            ("N", (1.18, -0.70, 0.0)), ("N", (-1.18, -0.70, 0.0)),
            ("H", (1.22, -1.70, 0.0)), ("H", (2.05, -0.22, 0.0)),
            ("H", (-1.22, -1.70, 0.0)), ("H", (-2.05, -0.22, 0.0))]


def h2o_atoms():
    return [("O", (0.0, 0.0, 0.0)), ("H", (0.96, 0.0, 0.0)),
            ("H", (-0.24, 0.93, 0.0))]


def _load_res():
    if os.path.exists(RES):
        try:
            return json.load(open(RES))
        except Exception:
            pass
    return {}


def _save_res(res):
    json.dump(res, open(RES, "w"), indent=1)


# ---------------------------------------------------------------- стадии
def stage_thermo():
    res = _load_res()
    out = res.get("thermo", {})
    species = {"co2": co2_atoms(), "nh3": nh3_atoms(),
               "carbamic": carbamic_atoms(), "urea": urea_atoms(),
               "h2o": h2o_atoms()}
    for tag, at in species.items():
        if tag in out and out[tag].get("converged"):
            print(f"[thermo] {tag}: resume", flush=True)
            continue
        t0 = time.time()
        na, e, conv = relax(at)
        out[tag] = {"e_h": e, "converged": conv,
                    "atoms": [[s, [round(x, 6) for x in c]] for s, c in na],
                    "wall_s": round(time.time() - t0, 1)}
        res["thermo"] = out
        _save_res(res)
        print(f"[thermo] {tag}: e={e:.6f} conv={conv}", flush=True)
    e = {k: out[k]["e_h"] for k in species}
    res["thermo_kcal"] = {
        "addition_NH3_CO2_to_carbamic": round(
            (e["carbamic"] - e["co2"] - e["nh3"]) * HARTREE_KCAL, 2),
        "condensation_carbamic_NH3_to_urea_H2O": round(
            (e["urea"] + e["h2o"] - e["carbamic"] - e["nh3"]) * HARTREE_KCAL, 2),
        "overall_2NH3_CO2_to_urea_H2O": round(
            (e["urea"] + e["h2o"] - e["co2"] - 2 * e["nh3"]) * HARTREE_KCAL, 2)}
    _save_res(res)
    print("[thermo]", json.dumps(res["thermo_kcal"]), flush=True)


# путь: подход C–N + перенос H на O (согласованный 4-членный TS)
GRID = ((2.60, 2.30), (2.30, 2.00), (2.05, 1.75), (1.85, 1.50),
        (1.70, 1.30), (1.60, 1.12), (1.55, 1.00), (1.52, 0.98))


def stage_ts():
    from pyscf.geomopt.geometric_solver import optimize
    res = _load_res()
    th = res.get("thermo", {})
    if not ("co2" in th and "nh3" in th):
        raise RuntimeError("ts: сначала stage thermo")
    eR = th["co2"]["e_h"] + th["nh3"]["e_h"]

    pts = []
    if os.path.exists(SCAN):
        try:
            pts = json.load(open(SCAN))["pts"]
            print(f"[ts] resume: {len(pts)} точек", flush=True)
        except Exception:
            pts = []
    done = {(p["rCN"], p["rOH"]) for p in pts}
    dm, start = None, None
    for rCN, rOH in GRID:
        if (rCN, rOH) in done:
            saved = next(p for p in pts if (p["rCN"], p["rOH"]) == (rCN, rOH))
            start = [(s, tuple(c)) for s, c in saved["atoms"]]
            continue
        cf = os.path.join(DIR, f"_urea_ts_{rCN:.2f}.txt")
        open(cf, "w").write(
            f"$set\ndistance 1 4 {rCN:.4f}\ndistance 2 5 {rOH:.4f}\n")
        try:
            st = start if start is not None else adduct_atoms(rCN)
            mf0 = uks(M(st), dm0=dm)
            mol = optimize(mf0, constraints=cf, maxsteps=150,
                           assert_convergence=False)
            mfx = uks(mol, dm0=mf0.make_rdm1())
            na = [(mol.atom_symbol(i), tuple(c)) for i, c in
                  enumerate(mol.atom_coords(unit="Angstrom"))]
            dm = mfx.make_rdm1()
            start = na
            pts.append({"rCN": rCN, "rOH": rOH, "e_h": float(mfx.e_tot),
                        "conv": bool(mfx.converged), "atoms": na})
            json.dump({"pts": pts}, open(SCAN, "w"), indent=1)
            print(f"  [ts] rCN={rCN} rOH={rOH} rel="
                  f"{(mfx.e_tot - eR) * HARTREE_KCAL:.2f} "
                  f"conv={mfx.converged}", flush=True)
        except Exception as ex:
            print(f"  [ts] rCN={rCN} failed: {ex}", flush=True)
        finally:
            if os.path.exists(cf):
                os.remove(cf)

    order = {g: k for k, g in enumerate(GRID)}
    pts.sort(key=lambda p: order.get((p["rCN"], p["rOH"]), 99))
    ok = [p for p in pts if p["conv"]]
    rel = [(p["e_h"] - eR) * HARTREE_KCAL for p in ok]
    keep = [i for i in range(len(ok)) if not all(
        rel[i] - rel[j] > 25.0 for j in (i - 1, i + 1) if 0 <= j < len(ok))]
    kr = [rel[i] for i in keep]
    out = {"scan_rel_kcal": [round(x, 2) for x in rel]}
    if len(kr) >= 4:
        imax = max(range(1, len(kr) - 1), key=lambda i: kr[i])
        interior = kr[imax] >= max(kr[0], kr[-1])
        best = ok[keep[imax]]
        out["dft"] = {"barrier_kcal": round(kr[imax], 2),
                      "ts_rCN": best["rCN"], "ts_rOH": best["rOH"],
                      "ts_interior": bool(interior)}
        if interior and 0 < kr[imax] < 100:
            lab = ["0 C 2p", "1 O 2p", "3 N 2p", "4 H 1s"]
            try:
                n_ts = nevpt2(best["atoms"], lab, thr=0.45)
                n_co2 = nevpt2([(s, tuple(c)) for s, c in
                                th["co2"]["atoms"]], ["C 2p", "O 2p"], thr=0.45)
                n_nh3 = nevpt2([(s, tuple(c)) for s, c in
                                th["nh3"]["atoms"]], ["N 2p", "H 1s"], thr=0.45)
                out["nevpt2_control"] = {
                    "ts": n_ts, "co2": n_co2, "nh3": n_nh3,
                    "barrier_casscf": round(
                        (n_ts["e_casscf"] - n_co2["e_casscf"]
                         - n_nh3["e_casscf"]) * HARTREE_KCAL, 2),
                    "barrier_nevpt2": round(
                        (n_ts["e_nevpt2"] - n_co2["e_nevpt2"]
                         - n_nh3["e_nevpt2"]) * HARTREE_KCAL, 2)}
            except Exception as e:
                out["nevpt2_error"] = str(e)[:200]
        else:
            out["ts_failed"] = "TS not interior — no verdict"
    else:
        out["ts_failed"] = "too few points"
    res = _load_res()
    res["ts"] = out
    _save_res(res)
    print("[ts]", json.dumps({k: v for k, v in out.items()},
                             default=str)[:500], flush=True)


def stage_merge():
    res = _load_res()
    ts = res.get("ts", {})
    ctl = ts.get("nevpt2_control") or {}
    verdict = None
    if ctl:
        dft_b = (ts.get("dft") or {}).get("barrier_kcal")
        nev_b = ctl.get("barrier_nevpt2")
        frac = (ctl.get("ts") or {}).get("n_frac_noon")
        if dft_b is not None and nev_b is not None:
            verdict = {
                "dft_vs_nevpt2_gap_kcal": round(nev_b - dft_b, 2),
                "ts_n_frac_noon": frac,
                "closed_shell_control": ("PASS: DFT≈NEVPT2, NOON ~2/0 — "
                                         "метод M не нужен на закрытых "
                                         "оболочках" if frac is not None
                                         and frac <= 1 and
                                         abs(nev_b - dft_b) < 8.0 else
                                         "ATTENTION: расхождение/NOON — "
                                         "смотреть руками")}
    res["method_M_negative_control"] = verdict
    res["honesty"] = ("Скрин PBE0/def2-SVP+DF, газофаза, без ZPE/T. "
                      "Промышленный синтез — жидкофазный карбамат при "
                      "~150 бар: газофазный барьер — якорь метода, не "
                      "модель реактора. Лит. газофазный барьер "
                      "NH3+CO2 ~35–40 ккал/моль (несопряжённый).")
    _save_res(res)
    print(json.dumps({"thermo_kcal": res.get("thermo_kcal"),
                      "dft": ts.get("dft"),
                      "control": verdict}, indent=1,
                     ensure_ascii=False), flush=True)


if __name__ == "__main__":
    st = sys.argv[1] if len(sys.argv) > 1 else "thermo"
    {"thermo": stage_thermo, "ts": stage_ts, "merge": stage_merge}[st]()
