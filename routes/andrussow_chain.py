#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Трек F (шок-кандидат S1, NaCN): газофазный радикальный скелет Андруссова.

Процесс Андруссова: CH₄ + NH₃ + O₂ → HCN на Pt при ~1100 °C; радикальная
газофаза — существенная часть селективности. Считаемое ядро (скрин):
  1) HAT-пропагация •NH₂ + CH₄ → NH₃ + CH₃• — коллинеарный пин-скан
     (rC–H, rN–H) по рецепту MSA v7 (поточечный чекпойнт, маска, интерьерный
     TS, sanity), DFT vs CASSCF vs NEVPT2 (симметричные AVAS-метки C/H/N).
     Валидационный якорь: лит. барьер ~13–15 ккал/моль.
  2) Рекомбинация CH₃• + •NH₂ → CH₃NH₂ (синглетная поверхность): релакс-скан
     r(C–N) — безбарьерность как санити радикального пула.
  3) Термохимия дегидрирования: CH₃NH₂ → CH₂=NH + H₂ → HCN + 2H₂ (шаги к
     HCN, которые на Pt идут каталитически) — реакционные энергии.
Модель: UKS-PBE0/def2-SVP+DF; NEVPT2 на ROHF→AVAS. Скрин-уровень; без
NOx-ветви (пламя O₂ — отдельная задача, помечено).

Стадии: hat | assoc | thermo | merge.
Запуск: python3 -u routes/andrussow_chain.py <stage>
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
SCAN = os.path.join(DIR, "andr_hat_scan.json")
RES = os.path.join(DIR, "andr_results.json")


def M(atoms, spin):
    return gto.M(atom=[(s, tuple(c)) for s, c in atoms], basis="def2-svp",
                 charge=0, spin=spin, verbose=0, max_memory=12000)


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


def relax(atoms, spin, maxsteps=150):
    from pyscf.geomopt.geometric_solver import optimize
    mf = uks(M(atoms, spin))
    mol = optimize(mf, maxsteps=maxsteps, assert_convergence=False)
    mfx = uks(mol, dm0=mf.make_rdm1())
    na = [(mol.atom_symbol(i), tuple(c))
          for i, c in enumerate(mol.atom_coords(unit="Angstrom"))]
    return na, float(mfx.e_tot), bool(mfx.converged)


def nevpt2(atoms, spin, labels, thr=0.5):
    mol = M(atoms, spin)
    mol.max_memory = 24000
    mf = scf.ROHF(mol)
    mf.conv_tol = 1e-8
    mfn = scf.newton(mf)
    mfn.kernel()
    ncas, nelec, mo = avas.avas(mfn, labels, threshold=thr)
    ss = spin / 2 * (spin / 2 + 1)
    mc = mcscf.CASSCF(mfn, ncas, nelec)
    fci.addons.fix_spin_(mc.fcisolver, shift=0.2, ss=ss)
    mc.max_cycle_macro = 150
    mc.conv_tol = 1e-7
    mc.kernel(mo)
    ept = NEVPT(mc).kernel()
    ne = (int(sum(nelec)) if hasattr(nelec, "__len__") else int(nelec))
    return {"cas": [ne, int(ncas)], "e_casscf": float(mc.e_tot),
            "e_nevpt2": float(mc.e_tot + ept),
            "casscf_conv": bool(mc.converged)}


# ------------------------------------------------------------- строители
def ch4_atoms():
    d = 1.09
    return [("C", (0, 0, 0)), ("H", (0, 0, d)),
            ("H", (d * 0.9428, 0, -d / 3)),
            ("H", (-d * 0.4714, d * 0.8165, -d / 3)),
            ("H", (-d * 0.4714, -d * 0.8165, -d / 3))]


def nh2_atoms():
    return [("N", (0, 0, 0)), ("H", (1.03, 0, 0)), ("H", (-0.35, 0.97, 0))]


def hat_ts(rCH, rNH):
    """Коллинеарный C–H···N старт: CH₄-остов + NH₂ на оси переноса."""
    d = 1.09
    at = [("C", (0.0, 0.0, 0.0)), ("H", (0.0, 0.0, rCH)),
          ("N", (0.0, 0.0, rCH + rNH)),
          ("H", (0.95, 0.35, rCH + rNH + 0.35)),
          ("H", (-0.95, 0.35, rCH + rNH + 0.35)),
          ("H", (d * 0.9428, 0, -d / 3)),
          ("H", (-d * 0.4714, d * 0.8165, -d / 3)),
          ("H", (-d * 0.4714, -d * 0.8165, -d / 3))]
    return at


# ---------------------------------------------------------------- стадии
def stage_hat():
    from pyscf.geomopt.geometric_solver import optimize
    t0 = time.time()
    out = {"stage": "hat", "model": "pinned-relax collinear HAT NH2+CH4"}
    ch4_geo, e_ch4, _ = relax(ch4_atoms(), 0)
    nh2_geo, e_nh2, _ = relax(nh2_atoms(), 1)
    eR = e_ch4 + e_nh2

    pts = []
    if os.path.exists(SCAN):
        try:
            pts = json.load(open(SCAN))["pts"]
            print(f"[hat] resume: {len(pts)} точек", flush=True)
        except Exception:
            pts = []
    done = {(p["rCH"], p["rNH"]) for p in pts}

    def save_scan():
        json.dump({"pts": pts}, open(SCAN, "w"), indent=1)

    def point(start, rCH, rNH, dm):
        cf = os.path.join(DIR, f"_andr_hat_{rCH:.2f}.txt")
        open(cf, "w").write(
            f"$set\ndistance 1 2 {rCH:.4f}\ndistance 2 3 {rNH:.4f}\n")
        try:
            mf0 = uks(M(start, 1), dm0=dm)
            mol = optimize(mf0, constraints=cf, maxsteps=150,
                           assert_convergence=False)
            mfx = uks(mol, dm0=mf0.make_rdm1())
            na = [(mol.atom_symbol(i), tuple(c)) for i, c in
                  enumerate(mol.atom_coords(unit="Angstrom"))]
            return float(mfx.e_tot), bool(mfx.converged), na, mfx.make_rdm1()
        except Exception as ex:
            print(f"  [hat] rCH={rCH} failed: {ex}", flush=True)
            return None, False, None, dm
        finally:
            if os.path.exists(cf):
                os.remove(cf)

    # путь от реагентов к продуктам (подход → перевал → NH3+CH3)
    GRID = ((1.12, 1.55), (1.18, 1.42), (1.24, 1.32), (1.30, 1.24),
            (1.38, 1.16), (1.50, 1.08), (1.65, 1.03))
    dm, start = None, None
    for rCH, rNH in GRID:
        if (rCH, rNH) in done:
            saved = next(p for p in pts
                         if (p["rCH"], p["rNH"]) == (rCH, rNH))
            start = [(s, tuple(c)) for s, c in saved["atoms"]]
            continue
        st = start if start is not None else hat_ts(rCH, rNH)
        e, conv, na, dm = point(st, rCH, rNH, dm)
        if e is not None and np.isfinite(e):
            pts.append({"rCH": rCH, "rNH": rNH, "e_h": e,
                        "scf_converged": conv, "atoms": na})
            save_scan()
            print(f"  [hat] rCH={rCH} rNH={rNH} rel="
                  f"{(e - eR) * HARTREE_KCAL:.2f} conv={conv}", flush=True)
            start = na
    pts.sort(key=lambda p: p["rCH"])
    ok = [p for p in pts if p["scf_converged"]]
    rel = [(p["e_h"] - eR) * HARTREE_KCAL for p in ok]
    out["scan_rel_kcal"] = [round(x, 2) for x in rel]
    keep = []
    for i, p in enumerate(ok):
        nb = [rel[j] for j in (i - 1, i + 1) if 0 <= j < len(ok)]
        if nb and all(rel[i] - x > 25.0 for x in nb):
            print(f"  [hat] outlier dropped rCH={p['rCH']}", flush=True)
            continue
        keep.append((p, rel[i]))
    if len(keep) < 4:
        raise RuntimeError("hat: too few points")
    imax = max(range(len(keep)), key=lambda i: keep[i][1])
    best, bar = keep[imax]
    interior = 0 < imax < len(keep) - 1
    out["dft"] = {"barrier_kcal": round(bar, 2), "ts_rCH": best["rCH"],
                  "ts_rNH": best["rNH"], "ts_interior": interior}
    if not interior or not (0.0 < bar < 100.0):
        out["hat_stage_failed"] = "TS not bracketed — no verdict"
    else:
        lab_ts = ["0 C 2p", "1 H 1s", "2 N 2p"]
        try:
            n_ts = nevpt2(best["atoms"], 1, lab_ts, thr=0.45)
            n_ch4 = nevpt2(ch4_geo, 0, ["0 C 2p", "1 H 1s"], thr=0.6)
            n_nh2 = nevpt2(nh2_geo, 1, ["0 N 2p"], thr=0.5)
            for lv in ("casscf", "nevpt2"):
                eRc = n_ch4[f"e_{lv}"] + n_nh2[f"e_{lv}"]
                out[lv] = {"barrier_kcal": round(
                    (n_ts[f"e_{lv}"] - eRc) * HARTREE_KCAL, 2)}
            out["ts_cas"] = n_ts["cas"]
        except Exception as e:
            out["nevpt2_error"] = str(e)[:200]
    out["wall_s"] = round(time.time() - t0, 1)
    json.dump(out, open(os.path.join(DIR, "andr_hat.json"), "w"), indent=1)
    print("[hat]", json.dumps({k: v for k, v in out.items()
                               if k not in ("scan_rel_kcal",)},
                              ensure_ascii=False)[:400], flush=True)


def stage_assoc():
    """CH3+NH2 рекомбинация: релакс-скан r(C–N) на синглете."""
    from pyscf.geomopt.geometric_solver import optimize
    ch3, e_ch3, _ = relax([("C", (0, 0, 0)), ("H", (1.08, 0, 0)),
                           ("H", (-0.54, 0.94, 0)), ("H", (-0.54, -0.94, 0))],
                          1)
    nh2, e_nh2, _ = relax(nh2_atoms(), 1)
    eR = e_ch3 + e_nh2
    # сборка: CH3 у нуля, NH2 сдвинут по z на стартовые 3.2 Å
    at = ([(s, (c[0], c[1], c[2])) for s, c in ch3]
          + [(s, (c[0], c[1], c[2] + 3.2)) for s, c in nh2])
    pts = []
    dm = None
    for r in (3.2, 2.8, 2.4, 2.1, 1.9, 1.7, 1.55, 1.47):
        cf = os.path.join(DIR, f"_andr_cn_{r:.2f}.txt")
        open(cf, "w").write(f"$set\ndistance 1 5 {r:.4f}\n")
        try:
            mf0 = uks(M(at, 0), dm0=dm)
            mol = optimize(mf0, constraints=cf, maxsteps=120,
                           assert_convergence=False)
            mfx = uks(mol, dm0=mf0.make_rdm1())
            at = [(mol.atom_symbol(i), tuple(c)) for i, c in
                  enumerate(mol.atom_coords(unit="Angstrom"))]
            dm = mfx.make_rdm1()
            pts.append({"r": r, "rel_kcal": round(float(
                (mfx.e_tot - eR) * HARTREE_KCAL), 2),
                "conv": bool(mfx.converged)})
            print(f"  [assoc] r={r} rel={pts[-1]['rel_kcal']}", flush=True)
        except Exception as ex:
            print(f"  [assoc] r={r} FAIL: {ex}", flush=True)
        finally:
            if os.path.exists(cf):
                os.remove(cf)
    rels = [float(p["rel_kcal"]) for p in pts if p["conv"]]
    barrierless = bool(rels) and bool(max(rels) <= 1.0)
    json.dump({"pts": pts, "barrierless_dft": barrierless},
              open(os.path.join(DIR, "andr_assoc.json"), "w"), indent=1)
    print(f"[assoc] barrierless={barrierless}", flush=True)


def stage_thermo():
    """Термохимия к HCN (шаги дегидрирования, на Pt — каталитические)."""
    res = {}
    species = {
        "ch3nh2": ([("C", (0, 0, 0)), ("N", (1.47, 0, 0)),
                    ("H", (-0.4, 1.0, 0)), ("H", (-0.4, -0.5, 0.9)),
                    ("H", (-0.4, -0.5, -0.9)), ("H", (1.85, 0.9, 0.2)),
                    ("H", (1.85, -0.8, 0.4))], 0),
        "ch2nh": ([("C", (0, 0, 0)), ("N", (1.27, 0, 0)),
                   ("H", (-0.55, 0.93, 0)), ("H", (-0.55, -0.93, 0)),
                   ("H", (1.85, 0.85, 0))], 0),
        "hcn": ([("H", (0, 0, -1.07)), ("C", (0, 0, 0)),
                 ("N", (0, 0, 1.15))], 0),
        "h2": ([("H", (0, 0, 0)), ("H", (0.74, 0, 0))], 0),
    }
    for name, (at, sp) in species.items():
        na, e, conv = relax(at, sp)
        res[name] = {"e_h": e, "conv": conv}
        print(f"  [thermo] {name}: {e:.6f} conv={conv}", flush=True)
    dE1 = (res["ch2nh"]["e_h"] + res["h2"]["e_h"]
           - res["ch3nh2"]["e_h"]) * HARTREE_KCAL
    dE2 = (res["hcn"]["e_h"] + res["h2"]["e_h"]
           - res["ch2nh"]["e_h"]) * HARTREE_KCAL
    res["steps_kcal"] = {"ch3nh2->ch2nh+h2": round(dE1, 2),
                         "ch2nh->hcn+h2": round(dE2, 2)}
    json.dump(res, open(os.path.join(DIR, "andr_thermo.json"), "w"), indent=1)
    print(f"[thermo] {res['steps_kcal']}", flush=True)


def _drop_outliers(pts, key="rel_kcal", thr=25.0):
    """Санити из stage_hat: убрать точки, отклоняющиеся >thr ккал от ОБОИХ
    соседей (SCF-скачки состояний на многоповерхностных системах)."""
    ok = [p for p in pts if p.get("conv")]
    if len(ok) < 3:
        return ok
    keep = []
    for i, p in enumerate(ok):
        nb = [ok[j][key] for j in (i - 1, i + 1) if 0 <= j < len(ok)]
        if nb and all(abs(p[key] - v) > thr for v in nb):
            continue
        keep.append(p)
    return keep


def nh_atoms():
    return [("N", (0, 0, 0)), ("H", (1.04, 0, 0))]


def o2_atoms():
    return [("O", (0, 0, 0)), ("O", (1.21, 0, 0))]


def stage_branch():
    """S1 NaCN: дескриптор ветвления HCN vs NOx — судьба азота на развилке.
    (A) C–N сочетание CH3•+NH(триплет) на дублете — путь к HCN (продукт);
    (B) окисление NH+O2 на суммарном триплете — путь к HNO/NO (потеря в NOx).
    Нитрен NH — MR-зона (триплет ground) => барьер B на DFT/CASSCF/NEVPT2.
    Дескриптор: delta_sel = Ea(B) − Ea(A); Ea(A)=0, если сочетание безбарьерно."""
    from pyscf.geomopt.geometric_solver import optimize
    t0 = time.time()
    out = {"stage": "branch",
           "model": "NH branching: CH3+NH coupling (doublet) vs NH+O2 "
                    "oxidation (triplet), gas-phase screen"}

    # --- ветвь A: CH3 + NH сочетание (дублетная поверхность) ---
    ch3, e_ch3, _ = relax([("C", (0, 0, 0)), ("H", (1.08, 0, 0)),
                           ("H", (-0.54, 0.94, 0)), ("H", (-0.54, -0.94, 0))],
                          1)
    nh, e_nh, _ = relax(nh_atoms(), 2)
    eRA = e_ch3 + e_nh
    at = ([(s, (c[0], c[1], c[2])) for s, c in ch3]
          + [(s, (c[0], c[1], c[2] + 3.2)) for s, c in nh])
    ptsA, dm = [], None
    for r in (3.2, 2.8, 2.4, 2.1, 1.9, 1.7, 1.55, 1.47):
        cf = os.path.join(DIR, f"_andr_cnh_{r:.2f}.txt")
        open(cf, "w").write(f"$set\ndistance 1 5 {r:.4f}\n")
        try:
            mf0 = uks(M(at, 1), dm0=dm)
            mol = optimize(mf0, constraints=cf, maxsteps=120,
                           assert_convergence=False)
            mfx = uks(mol, dm0=mf0.make_rdm1())
            at = [(mol.atom_symbol(i), tuple(c)) for i, c in
                  enumerate(mol.atom_coords(unit="Angstrom"))]
            dm = mfx.make_rdm1()
            ptsA.append({"r": r, "rel_kcal": round(float(
                (mfx.e_tot - eRA) * HARTREE_KCAL), 2),
                "conv": bool(mfx.converged)})
            print(f"  [branch:A] r={r} rel={ptsA[-1]['rel_kcal']}", flush=True)
        except Exception as ex:
            print(f"  [branch:A] r={r} FAIL: {ex}", flush=True)
        finally:
            if os.path.exists(cf):
                os.remove(cf)
    keptA = _drop_outliers(ptsA)
    relsA = [float(p["rel_kcal"]) for p in keptA]
    # Ea = интерьерный максимум ПОСЛЕ дропа выбросов; монотонный спуск
    # (сочетание радикалов) => безбарьерно, Ea=0
    ea_A_dft = None
    if len(relsA) >= 3:
        imaxA = max(range(1, len(relsA) - 1), key=lambda i: relsA[i])
        interiorA = relsA[imaxA] >= max(relsA[0], relsA[-1])
        ea_A_dft = max(0.0, relsA[imaxA]) if interiorA else 0.0
    out["cn"] = {"pts": ptsA, "n_dropped": len([p for p in ptsA
                                                if p["conv"]]) - len(keptA),
                 "barrierless_dft": (ea_A_dft is not None
                                     and ea_A_dft <= 1.0),
                 "ea_dft_kcal": (round(ea_A_dft, 2)
                                 if ea_A_dft is not None else None)}

    # --- ветвь B: NH + O2 окисление (триплетная суммарная поверхность) ---
    nh_r, e_nh2p, _ = relax(nh_atoms(), 2)
    o2, e_o2, _ = relax(o2_atoms(), 2)
    eRB = e_nh2p + e_o2
    # сборка: N в нуле (H назад под углом), O2 по оси z
    atB = [("N", (0.0, 0.0, 0.0)), ("H", (-0.74, -0.73, 0.0)),
           ("O", (0.0, 0.0, 2.6)), ("O", (0.0, 0.0, 3.81))]
    ptsB, dm = [], None
    for r in (2.6, 2.3, 2.0, 1.8, 1.65, 1.5, 1.4):
        cf = os.path.join(DIR, f"_andr_no_{r:.2f}.txt")
        open(cf, "w").write(f"$set\ndistance 1 3 {r:.4f}\n")
        try:
            mf0 = uks(M(atB, 2), dm0=dm)
            mol = optimize(mf0, constraints=cf, maxsteps=120,
                           assert_convergence=False)
            mfx = uks(mol, dm0=mf0.make_rdm1())
            atB = [(mol.atom_symbol(i), tuple(c)) for i, c in
                   enumerate(mol.atom_coords(unit="Angstrom"))]
            dm = mfx.make_rdm1()
            ptsB.append({"r": r, "rel_kcal": round(float(
                (mfx.e_tot - eRB) * HARTREE_KCAL), 2),
                "conv": bool(mfx.converged),
                "atoms": [[s, [round(x, 6) for x in c]] for s, c in atB]})
            print(f"  [branch:B] r={r} rel={ptsB[-1]['rel_kcal']}", flush=True)
        except Exception as ex:
            print(f"  [branch:B] r={r} FAIL: {ex}", flush=True)
        finally:
            if os.path.exists(cf):
                os.remove(cf)
    okB = _drop_outliers(ptsB)
    out["no2"] = {"pts": [{k: p[k] for k in ("r", "rel_kcal", "conv")}
                          for p in ptsB],
                  "n_dropped": len([p for p in ptsB if p["conv"]]) - len(okB)}
    if len(okB) >= 4:
        relsB = [float(p["rel_kcal"]) for p in okB]
        imax = max(range(1, len(relsB) - 1), key=lambda i: relsB[i]) \
            if len(relsB) >= 3 else None
        interior = (imax is not None
                    and relsB[imax] >= max(relsB[0], relsB[-1]))
        out["no2"]["dft_barrier_kcal"] = (round(max(0.0, relsB[imax]), 2)
                                          if imax is not None else None)
        out["no2"]["ts_interior"] = bool(interior)
        if interior:
            ts_at = [(s, tuple(c)) for s, c in okB[imax]["atoms"]]
            try:
                n_ts = nevpt2(ts_at, 2, ["0 N 2p", "2 O 2p", "3 O 2p"],
                              thr=0.45)
                n_nh = nevpt2(nh_r, 2, ["0 N 2p"], thr=0.5)
                n_o2 = nevpt2(o2, 2, ["0 O 2p", "1 O 2p"], thr=0.5)
                for lv, key in (("e_casscf", "casscf"),
                                ("e_nevpt2", "nevpt2")):
                    out["no2"][f"{key}_barrier_kcal"] = round(
                        (n_ts[lv] - n_nh[lv] - n_o2[lv]) * HARTREE_KCAL, 2)
                out["no2"]["ts_cas"] = n_ts["cas"]
            except Exception as ex:
                out["no2"]["nevpt2_error"] = str(ex)[:200]

    # --- термохимия обеих NOx-развязок ---
    try:
        hno, e_hno, _ = relax([("H", (-0.99, 0.17, 0)), ("N", (0, 0, 0)),
                               ("O", (1.21, 0, 0))], 0)
        no, e_no, _ = relax([("N", (0, 0, 0)), ("O", (1.15, 0, 0))], 1)
        oh, e_oh, _ = relax([("O", (0, 0, 0)), ("H", (0.97, 0, 0))], 1)
        mo_ = uks(M([("O", (0.0, 0.0, 0.0))], 2))
        e_o = float(mo_.e_tot)
        out["thermo_kcal"] = {
            "nh+o2->hno+o": round((e_hno + e_o - eRB) * HARTREE_KCAL, 2),
            "nh+o2->no+oh": round((e_no + e_oh - eRB) * HARTREE_KCAL, 2)}
    except Exception as ex:
        out["thermo_error"] = str(ex)[:200]

    # --- дескриптор селективности ---
    ea_A = out["cn"].get("ea_dft_kcal") or 0.0
    for lv in ("dft", "casscf", "nevpt2"):
        eb = out["no2"].get(f"{lv}_barrier_kcal")
        if eb is not None:
            out[f"selectivity_descriptor_{lv}_kcal"] = round(eb - ea_A, 2)
    out["descriptor_note"] = ("delta_sel = Ea(NH+O2 -> NOx) − Ea(NH+CH3 -> "
                              "C–N); >0 = азот уходит в HCN-русло (хорошо "
                              "для NaCN-скида), <0 = теряется в NOx")
    out["wall_s"] = round(time.time() - t0, 1)
    json.dump(out, open(os.path.join(DIR, "andr_branch.json"), "w"), indent=1)
    print("[branch]", json.dumps({k: v for k, v in out.items()
                                  if k not in ("cn", "no2")},
                                 ensure_ascii=False)[:400], flush=True)


def stage_merge():
    def _ld(f):
        p = os.path.join(DIR, f)
        if os.path.exists(p):
            return json.load(open(p))
        return {"missing": f}
    hat, asc, th, br = (_ld("andr_hat.json"), _ld("andr_assoc.json"),
                        _ld("andr_thermo.json"), _ld("andr_branch.json"))
    res = {"route": "Andrussow radical skeleton: NH2-HAT + CH3/NH2 "
                    "recombination + dehydro thermochemistry to HCN "
                    "+ NH branching (HCN vs NOx)",
           "hat_barrier": {lv: hat.get(lv, {}).get("barrier_kcal")
                           for lv in ("dft", "casscf", "nevpt2")},
           "hat_lit_anchor_kcal": "13-15 (NH2+CH4, лит.)",
           "hat_failed": hat.get("hat_stage_failed"),
           "recombination_barrierless_dft": asc.get("barrierless_dft"),
           "dehydro_steps_kcal": th.get("steps_kcal"),
           "branching_descriptor": (
               {lv: br.get(f"selectivity_descriptor_{lv}_kcal")
                for lv in ("dft", "casscf", "nevpt2")}
               if "missing" not in br else br),
           "branching_note": br.get("descriptor_note"),
           "nox_thermo_kcal": br.get("thermo_kcal"),
           "note": "скрин-уровень; Pt-поверхность заменена газофазным "
                   "скелетом — селективность каталитической стадии не "
                   "утверждается; ветвление NH считано стадией branch"}
    json.dump(res, open(RES, "w"), indent=1)
    print(json.dumps(res, indent=1, ensure_ascii=False), flush=True)


if __name__ == "__main__":
    {"hat": stage_hat, "assoc": stage_assoc, "thermo": stage_thermo,
     "branch": stage_branch, "merge": stage_merge}[sys.argv[1]]()
