#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Этап 21-B: летучесть допанта (риск A2) — Cr vs Mn vs W.

Вопрос: не улетает ли Cr из Na/W-катализатора при 800 °C водяным паром как
CrO₂(OH)₂(g) (известный механизм деградации Cr-оксидов). Сравнительный
дескриптор: ΔE_vol(M) = E[MO₂(OH)₂(g)] − E[каркас–M=O]; сравнение против Mn
(эталон стабильной системы Mn–Na₂WO₄): ΔΔ(M−Mn) = ΔE_vol(M) − ΔE_vol(Mn).
Замороженные каркасы Cr- и Mn-кластеров СОВПАДАЮТ (RMSD 0.0, проверено) ⇒
члены вакансии/H₂O/H₂ сокращаются в ΔΔ точно. Отрицательная ΔΔ(M−Mn) =
M улетает легче Mn. Kill-критерий A2: Cr легче Mn на >20 ккал ⇒ красный флаг
и fallback Mn+вольтаж.

Справочно считается и абсолютная (редокс) реакция
site–M=O + 3H₂O(g) → каркас-вакансия* + MO₂(OH)₂(g) + 2H₂(g)
без энергии вакансии (одинакова для всех M — не нужна для сравнения).

Модель: UKS-PBE0/def2-SVP + DF; сайты — из закоммиченных R-геометрий
(CH₄ вертикально удалён; M и оксо релаксируются, каркас заморожен $freeze);
спин-лестницы у всех структур; NEVPT2 (AVAS M d + оксо O 2p) на лучших
спинах сайтов и газов — DFT-vs-NEVPT2 колонка дескриптора.

Стадии: gas | site | nev | desc. Инкрементальный ocm_vol_results.json,
resume по ключам. Запуск: python3 -u routes/ocm_volatility.py <stage>
"""
import json
import os
import sys
import time

import numpy as np
from pyscf import dft, gto, scf

HARTREE_KCAL = 627.509474
DIR = os.path.dirname(os.path.abspath(__file__))
RES = os.path.join(DIR, "ocm_vol_results.json")
BASIS = "def2-svp"
_DSH = {"Cr": "3d", "Mn": "3d", "W": "5d"}


def _load():
    if os.path.exists(RES):
        try:
            return json.load(open(RES))
        except Exception:
            pass
    return {}


def _save(d):
    json.dump(d, open(RES, "w"), indent=1)


def M(atoms, spin):
    return gto.M(atom=[(s, tuple(c)) for s, c in atoms], basis=BASIS,
                 ecp=BASIS, charge=0, spin=spin, verbose=0, max_memory=12000)


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


def relax(atoms, spin, freeze_1based=None, maxsteps=150):
    """geomeTRIC-релакс; freeze_1based — список замороженных атомов."""
    from pyscf.geomopt.geometric_solver import optimize
    mf = uks(M(atoms, spin))
    if freeze_1based:
        cf = os.path.join(DIR, f"_vol_freeze_{os.getpid()}.txt")
        open(cf, "w").write(
            "$freeze\nxyz " + ",".join(str(i) for i in freeze_1based) + "\n")
        try:
            mol = optimize(mf, constraints=cf, maxsteps=maxsteps,
                           assert_convergence=False)
        finally:
            os.remove(cf)
    else:
        mol = optimize(mf, maxsteps=maxsteps, assert_convergence=False)
    mfx = uks(mol, dm0=mf.make_rdm1())
    na = [(mol.atom_symbol(i), tuple(c))
          for i, c in enumerate(mol.atom_coords(unit="Angstrom"))]
    return na, float(mfx.e_tot), bool(mfx.converged)


def ladder(build_atoms, spins, freeze=None, tag=""):
    """Спин-лестница с релаксом на каждом спине; возвращает лучшую."""
    best = None
    out = {}
    for s in spins:
        try:
            na, e, conv = relax(build_atoms, s, freeze)
            out[str(s)] = {"e_h": e, "conv": conv}
            print(f"  [{tag}] 2S={s}: E={e:.6f} conv={conv}", flush=True)
            if conv and (best is None or e < best[1]):
                best = (s, e, na)
        except Exception as ex:
            out[str(s)] = {"error": str(ex)[:120]}
            print(f"  [{tag}] 2S={s} FAIL: {ex}", flush=True)
    return best, out


# ------------------------------------------------------------- строители
def gas_mo2oh2(metal, d_oxo, d_oh):
    """MO₂(OH)₂ — тетраэдрический старт (хромовая кислота и гомологи)."""
    t = 109.47 / 2 * np.pi / 180
    z, x = np.cos(t), np.sin(t)
    o1 = (d_oxo * x, 0.0, d_oxo * z)
    o2 = (-d_oxo * x, 0.0, d_oxo * z)
    o3 = (0.0, d_oh * x, -d_oh * z)
    o4 = (0.0, -d_oh * x, -d_oh * z)
    h3 = (0.8, d_oh * x + 0.55, -d_oh * z - 0.3)
    h4 = (0.8, -d_oh * x - 0.55, -d_oh * z - 0.3)
    return [(metal, (0.0, 0.0, 0.0)), ("O", o1), ("O", o2),
            ("O", o3), ("H", h3), ("O", o4), ("H", h4)]


def load_site(metal):
    """Сайт каркас–M=O: R-геометрия минус CH₄ (вертикально)."""
    src = ("ocm_mnw_emb_cr_ch4_geom.json" if metal == "Cr"
           else "ocm_mnw_emb_ch4_geom.json")
    g = json.load(open(os.path.join(DIR, src)))
    atoms = g["points"][0]["atoms"]
    keep = [0, 1] + list(range(7, len(atoms)))
    site = [(s, tuple(c)) for s, c in (atoms[i] for i in keep)]
    if metal == "W":
        site[0] = ("W", site[0][1])
    return site


GAS_PARAMS = {"Cr": (1.60, 1.75), "Mn": (1.58, 1.76), "W": (1.73, 1.89)}
GAS_SPINS = {"Cr": (0, 2), "Mn": (1, 3), "W": (0, 2)}
SITE_SPINS = {"Cr": (2, 4, 6), "Mn": (3, 5, 7), "W": (1, 3, 5)}


# ---------------------------------------------------------------- стадии
def stage_gas():
    d = _load()
    g = d.setdefault("gas", {})
    for name, atoms, spins in (
            ("h2o", [("O", (0, 0, 0)), ("H", (0.96, 0, 0)),
                     ("H", (-0.24, 0.93, 0))], (0,)),
            ("h2", [("H", (0, 0, 0)), ("H", (0.74, 0, 0))], (0,))):
        if name in g and "e_h" in g[name]:
            continue
        best, lad = ladder(atoms, spins, tag=name)
        g[name] = {"ladder": lad}
        if best:
            g[name].update({"spin_2S": best[0], "e_h": best[1],
                            "atoms": [[s, list(c)] for s, c in best[2]]})
        _save(d)
    for m in ("Cr", "Mn", "W"):
        key = f"{m.lower()}o2oh2"
        if key in g and "e_h" in g[key]:
            continue
        do, doh = GAS_PARAMS[m]
        best, lad = ladder(gas_mo2oh2(m, do, doh), GAS_SPINS[m], tag=key)
        g[key] = {"ladder": lad}
        if best:
            g[key].update({"spin_2S": best[0], "e_h": best[1],
                           "atoms": [[s, list(c)] for s, c in best[2]]})
        _save(d)
    print("[gas] готово", flush=True)


def stage_site():
    d = _load()
    st = d.setdefault("site", {})
    for m in ("Cr", "Mn", "W"):
        if m in st and "e_h" in st[m]:
            continue
        site = load_site(m)
        freeze = list(range(3, len(site) + 1))   # всё кроме M(1) и оксо(2)
        best, lad = ladder(site, SITE_SPINS[m], freeze=freeze, tag=f"site-{m}")
        st[m] = {"ladder": lad, "frozen_frame": True}
        if best:
            st[m].update({"spin_2S": best[0], "e_h": best[1],
                          "atoms": [[s, list(c)] for s, c in best[2]]})
        _save(d)
    print("[site] готово", flush=True)


def stage_nev():
    """NEVPT2 (AVAS M d + O 2p) на лучших спинах сайтов и газов MO₂(OH)₂."""
    from pyscf import fci, mcscf
    from pyscf.mcscf import avas
    from pyscf.mrpt import NEVPT
    d = _load()
    nv = d.setdefault("nev", {})
    jobs = []
    for m in ("Cr", "Mn", "W"):
        if d.get("site", {}).get(m, {}).get("e_h") is not None:
            jobs.append((f"site_{m}", d["site"][m], m))
        gk = f"{m.lower()}o2oh2"
        if d.get("gas", {}).get(gk, {}).get("e_h") is not None:
            jobs.append((f"gas_{m}", d["gas"][gk], m))
    for tag, rec, m in jobs:
        if tag in nv and "e_nevpt2_h" in nv[tag]:
            continue
        t0 = time.time()
        try:
            atoms = [(s, tuple(c)) for s, c in rec["atoms"]]
            spin = rec["spin_2S"]
            mol = M(atoms, spin)
            mol.max_memory = 24000
            mf = scf.ROHF(mol)
            mf.conv_tol = 1e-8
            mfn = scf.newton(mf)
            mfn.kernel()
            ncas, nelec, mo = avas.avas(
                mfn, [f"{m} {_DSH[m]}", "1 O 2p"], threshold=0.5)
            ss = spin / 2 * (spin / 2 + 1)
            mc = mcscf.CASSCF(mfn, ncas, nelec)
            fci.addons.fix_spin_(mc.fcisolver, shift=0.2, ss=ss)
            mc.max_cycle_macro = 150
            mc.conv_tol = 1e-7
            mc.kernel(mo)
            ept = NEVPT(mc).kernel()
            nv[tag] = {"cas": [int(sum(nelec) if hasattr(nelec, "__len__")
                                   else nelec), int(ncas)],
                       "e_casscf_h": float(mc.e_tot),
                       "e_nevpt2_h": float(mc.e_tot + ept),
                       "casscf_conv": bool(mc.converged),
                       "wall_s": round(time.time() - t0, 1)}
            print(f"  [nev] {tag}: NEVPT2 {nv[tag]['e_nevpt2_h']:.6f} "
                  f"({nv[tag]['wall_s']}s)", flush=True)
        except Exception as ex:
            nv[tag] = {"error": str(ex)[:200]}
            print(f"  [nev] {tag} FAIL: {ex}", flush=True)
        _save(d)
    print("[nev] готово", flush=True)


def stage_desc():
    d = _load()
    out = {"descriptor": "dE_vol(M) = E[MO2(OH)2(g)] - E[site M=O]; "
                         "ddE(M-Mn) < 0 => M улетает легче Mn (риск A2)",
           "kill_criterion_A2": "ddE(Cr-Mn) < -20 kcal => красный флаг"}
    lv = {}

    def _e(kind, m, level):
        if kind == "site":
            rec = d.get("site", {}).get(m, {})
            tag = f"site_{m}"
        else:
            rec = d.get("gas", {}).get(f"{m.lower()}o2oh2", {})
            tag = f"gas_{m}"
        if level == "dft":
            return rec.get("e_h")
        return d.get("nev", {}).get(tag, {}).get(f"e_{level}_h")
    for level in ("dft", "casscf", "nevpt2"):
        row = {}
        for m in ("Cr", "Mn", "W"):
            eg, es = _e("gas", m, level), _e("site", m, level)
            row[m] = (round((eg - es) * HARTREE_KCAL, 2)
                      if eg is not None and es is not None else None)
        dd = {}
        for m in ("Cr", "W"):
            if row[m] is not None and row["Mn"] is not None:
                dd[f"{m}-Mn"] = round(row[m] - row["Mn"], 2)
        lv[level] = {"dE_vol_kcal": row, "ddE_vs_Mn_kcal": dd}
        if level in ("dft", "nevpt2") and dd.get("Cr-Mn") is not None:
            lv[level]["verdict_A2"] = (
                "КРАСНЫЙ ФЛАГ: Cr летит существенно легче Mn"
                if dd["Cr-Mn"] < -20 else
                "A2 закрыт на этом уровне: Cr не летучее Mn критично")
    out["levels"] = lv
    # справочный абсолют (редокс, 3H2O -> 2H2), только DFT
    g = d.get("gas", {})
    if all(g.get(k, {}).get("e_h") is not None for k in ("h2o", "h2")):
        absr = {}
        for m in ("Cr", "Mn", "W"):
            eg, es = _e("gas", m, "dft"), _e("site", m, "dft")
            if eg is not None and es is not None:
                absr[m] = round((eg + 2 * g["h2"]["e_h"] - es
                                 - 3 * g["h2o"]["e_h"]) * HARTREE_KCAL, 2)
        out["abs_rxn_3H2O_to_2H2_kcal_no_vacancy"] = absr
    d["desc"] = out
    _save(d)
    print(json.dumps(out, indent=1, ensure_ascii=False), flush=True)


if __name__ == "__main__":
    {"gas": stage_gas, "site": stage_site, "nev": stage_nev,
     "desc": stage_desc}[sys.argv[1]]()
