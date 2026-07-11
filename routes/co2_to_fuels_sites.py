#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
routes/co2_to_fuels_sites.py — ВОПРОС «ОТ ОБРАТНОГО»: Cu vs сплав CuAl.

Сравнивает чистый Cu-сайт и scaling-ломающий сплав-сайт CuAl на координате
сочленения *CO–*CO. На каждом сайте релаксируются (PBE/def2-SVP, RI, замороженный
металл-димер = жёсткая модель поверхности) два эндпоинта:
  separated — два *CO atop, разведены;  coupled — *OCCO (связь C–C ~1.5 Å).
В релакс-геометриях считаются PBE и CASCI(π-многообразие, 12e,10o, AVAS) + NOON.

Дескриптор: энергия сочленения ΔE_couple = E(coupled) − E(separated), из PBE и из
CASCI; квантовая поправка Δ_corr = ΔE_couple(CASCI) − ΔE_couple(PBE). Вопрос:
меняет ли мультиреференс-поправка ПОРЯДОК сайтов (Cu vs CuAl)?

Ограничения: это РЕАКЦИОННАЯ энергия сочленения (термодинамика релакс-эндпоинтов), НЕ
барьер ΔE‡ (для седла нужен поиск TS) и НЕ фарадеевская эффективность (кластер ≠
электрод, нет constant-potential). Барьер + constant-potential — следующий AWS-прогон.

Запуск: python3 -u routes/co2_to_fuels_sites.py
"""
import os, json, time
import numpy as np
from scipy.optimize import minimize
from pyscf import gto, dft, mcscf
from pyscf.mcscf import avas

HARTREE_EV = 27.211386245988
HERE = os.path.dirname(os.path.abspath(__file__))
RESULTS = os.path.join(HERE, "co2_to_fuels_results.json")

METALS = {
    "Cu2":  [("Cu", (-1.225, 0., 0.)), ("Cu", (1.225, 0., 0.))],
    "CuAl": [("Cu", (-1.30, 0., 0.)), ("Al", (1.30, 0., 0.))],
}


def atoms_str(metal, ads):
    s = ""
    for el, (x, y, z) in metal:
        s += f"{el} {x} {y} {z}; "
    for el, (x, y, z) in zip(["C", "O", "C", "O"], ads):
        s += f"{el} {x} {y} {z}; "
    return s


def sep_ads(metal):
    x1 = metal[0][1][0]; x2 = metal[1][1][0]
    return np.array([[x1, 0, 1.85], [x1, 0, 3.00], [x2, 0, 1.85], [x2, 0, 3.00]])


def coup_ads(metal):
    return np.array([[-0.75, 0, 1.80], [-1.35, 0, 2.55],
                     [0.75, 0, 1.80], [1.35, 0, 2.55]])


def _mol(metal, ads):
    return gto.M(atom=atoms_str(metal, ads), basis="def2-svp", spin=0, verbose=0)


def _scf(mol, dm0=None):
    mf = dft.RKS(mol).density_fit(); mf.xc = "pbe"; mf.verbose = 0
    mf.conv_tol = 1e-8; mf.max_cycle = 120
    mf.kernel(dm0=dm0)
    if not mf.converged:                         # robust fallback: level-shift → newton
        mfl = dft.RKS(mol).density_fit(); mfl.xc = "pbe"; mfl.verbose = 0
        mfl.level_shift = 0.4; mfl.damp = 0.3; mfl.max_cycle = 300; mfl.conv_tol = 1e-6
        mfl.kernel(dm0=dm0)
        mf = mfl.newton(); mf.kernel(mfl.make_rdm1())
    return mf


def pbe_eg(metal, ads, dm0):
    mf = _scf(_mol(metal, ads), dm0)
    g = mf.nuc_grad_method().kernel()
    return mf.e_tot, g[len(metal):].ravel(), mf.make_rdm1()


def relax(metal, ads0, maxiter=30):
    cache = {"dm": None, "n": 0}

    def fun(x):
        e, gads, dm = pbe_eg(metal, x.reshape(-1, 3), cache["dm"])
        cache["dm"] = dm; cache["n"] += 1
        return e, gads

    res = minimize(fun, ads0.ravel(), jac=True, method="L-BFGS-B",
                   options={"maxiter": maxiter, "maxfun": 2 * maxiter,
                            "ftol": 1e-7, "gtol": 1e-4})
    ads = res.x.reshape(-1, 3)
    cc = float(np.linalg.norm(ads[0] - ads[2]))          # final C···C
    return ads, float(res.fun), cc, cache["n"]


def casci_pi(metal, ads):
    mf = _scf(_mol(metal, ads))
    no, ne, orbs = avas.avas(mf, ["C 2p", "O 2px", "O 2py"], threshold=0.4, verbose=0)
    mc = mcscf.CASCI(mf, no, ne); mc.verbose = 0; mc.kernel(orbs)
    dm = mc.fcisolver.make_rdm1(mc.ci, mc.ncas, mc.nelecas)
    n = np.linalg.eigvalsh(dm)[::-1]
    nu = float(np.sum(np.minimum(n, 2 - n)))
    return float(mf.e_tot), float(mc.e_tot), round(nu, 4), (int(ne), int(no))


def run_site(tag):
    metal = METALS[tag]
    print(f"\n--- site {tag} ---")
    out = {}
    for name, ads0 in [("separated", sep_ads(metal)), ("coupled", coup_ads(metal))]:
        t0 = time.time()
        ads, e_pbe, cc, nev = relax(metal, ads0)
        e_pbe2, e_casci, nu, cas = casci_pi(metal, ads)
        out[name] = dict(e_pbe=round(e_pbe, 6), e_casci=round(e_casci, 6),
                         nu_pi=nu, cc_final=round(cc, 3), cas=f"({cas[0]}e,{cas[1]}o)")
        print(f"  {name:10s} E_pbe={e_pbe:.5f} E_casci={e_casci:.5f} n_u(π)={nu:.3f} "
              f"C···C={cc:.2f}Å  [{nev} evals, {time.time()-t0:.0f}s]")
    dEp = (out["coupled"]["e_pbe"] - out["separated"]["e_pbe"]) * HARTREE_EV
    dEc = (out["coupled"]["e_casci"] - out["separated"]["e_casci"]) * HARTREE_EV
    out["dEcouple_pbe_eV"] = round(dEp, 3)
    out["dEcouple_casci_eV"] = round(dEc, 3)
    out["quantum_corr_eV"] = round(dEc - dEp, 3)
    print(f"  ΔE_couple: PBE={dEp:+.3f} eV  CASCI={dEc:+.3f} eV  Δ_corr={dEc-dEp:+.3f} eV")
    return out


if __name__ == "__main__":
    res = json.load(open(RESULTS)) if os.path.exists(RESULTS) else {}
    sites = res.get("sites", {})
    for tag in ("Cu2", "CuAl"):
        if tag in sites:                       # resume: skip completed sites
            print(f"--- site {tag}: already done, skipping ---"); continue
        sites[tag] = run_site(tag)
        res["sites"] = sites                   # persist after each site
        json.dump(res, open(RESULTS, "w"), ensure_ascii=False, indent=1)
        print(f"  saved site {tag} → results")
    # ordering question (once both sites present)
    if "Cu2" in sites and "CuAl" in sites:
        pbe_order = sorted(("Cu2", "CuAl"), key=lambda t: sites[t]["dEcouple_pbe_eV"])
        cas_order = sorted(("Cu2", "CuAl"), key=lambda t: sites[t]["dEcouple_casci_eV"])
        flipped = list(pbe_order) != list(cas_order)
        print("\n=== ORDERING (lower ΔE_couple = easier coupling) ===")
        print(f"  PBE  : {list(pbe_order)}")
        print(f"  CASCI: {list(cas_order)}")
        print(f"  multireference correction FLIPS site order: {flipped}")
        sites["_summary"] = dict(pbe_order=list(pbe_order), casci_order=list(cas_order),
                                 order_flipped=bool(flipped),
                                 note=("ΔE_couple = coupling REACTION energy at relaxed endpoints "
                                       "(frozen-metal cluster, def2-SVP, no constant-potential); "
                                       "NOT a barrier and NOT Faradaic efficiency. Barrier/TS = AWS."))
        res["sites"] = sites
        json.dump(res, open(RESULTS, "w"), ensure_ascii=False, indent=1)
    print("\nsaved → results;  done.")
