#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
routes/co2_to_fuels_ts_sigma.py — дирадикалоид-проверка σ-фронтира В СЕДЛЕ.

π-многообразие (12e,10o) несёт распределённую корреляцию CO и не изолирует
ОБРАЗУЮЩУЮСЯ связь C–C. Литературное утверждение «TS дирадикалоиден (NOON≈1.92/0.09)»
относится именно к σ-фронтиру связи C–C. Здесь в релакс-геометриях TS и reactant
(из co2_to_fuels_ts.py) считается CAS(4e,4o) [AVAS C 2px] → σ-фронтир n_u/strong:
действительно ли седло — локализованный дирадикал расщепляющейся связи C–C?

Запуск (после co2_to_fuels_ts.py): python3 -u routes/co2_to_fuels_ts_sigma.py
"""
import os, sys, json
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from co2_to_fuels_ts import build_ads, _scf
from co2_to_fuels_sites import _mol, METALS, RESULTS
from pyscf import mcscf
from pyscf.mcscf import avas


def sigma_nu(metal, ads):
    mf = _scf(_mol(metal, ads))
    no, ne, orbs = avas.avas(mf, ["C 2px"], threshold=0.01, verbose=0)
    mc = mcscf.CASCI(mf, no, ne); mc.verbose = 0; mc.kernel(orbs)
    dm = mc.fcisolver.make_rdm1(mc.ci, mc.ncas, mc.nelecas)
    n = np.linalg.eigvalsh(dm)[::-1]
    nu = float(np.sum(np.minimum(n, 2 - n)))
    strong = int(np.sum((n > 0.1) & (n < 1.9)))
    return round(nu, 4), strong, (int(ne), int(no)), [round(float(x), 3) for x in n]


if __name__ == "__main__":
    d = json.load(open(RESULTS)); tb = d.get("ts_barrier", {})
    for site in ("Cu2", "CuAl"):
        s = tb.get(site, {})
        if "profile" not in s or s.get("d_ts") is None:
            continue
        metal = METALS[site]
        prof = {round(p["d"], 2): p["free"] for p in s["profile"]}
        for label, dval in [("ts", s["d_ts"]), ("reactant", s["d_reactant"])]:
            ads = build_ads(dval, np.array(prof[round(dval, 2)]))
            nu, strong, cas, noon = sigma_nu(metal, ads)
            s[f"{label}_sigma_nu"] = nu
            s[f"{label}_sigma_strong"] = strong
            print(f"{site:5s} {label:8s} d={dval:.2f}  σ-frontier CAS{cas} "
                  f"n_u={nu:.3f} strong={strong}  NOON={noon}")
        tb[site] = s
    d["ts_barrier"] = tb
    json.dump(d, open(RESULTS, "w"), ensure_ascii=False, indent=1)
    print("saved σ-frontier (TS/reactant) → results")
