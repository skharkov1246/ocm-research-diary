#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
routes/co2_to_fuels_ts.py — AWS-кусок: РЕЛАКС-СЕДЛО (барьер ΔE‡) сочленения *CO–*CO.

Расширяет co2_to_fuels_sites.py от РЕАКЦИОННОЙ энергии к БАРЬЕРУ. На каждом сайте
(чистый Cu₂ vs сплав CuAl) — констрейн-релакс-скан координаты C···C: углероды
разводятся по x на фикс. расстояние d (драйвер), остальные DOF адсорбата (высоты C,
оба O) релаксируются (RI-PBE/def2-SVP, замороженный металл, аналитич. градиенты).
Максимум профиля = переходное состояние (TS). В TS:
  • CASCI(π, 12e,10o, AVAS) + NOON → дирадикалоиден ли TS? (ключевое утверждение лит.)
  • ΔE‡(PBE) и ΔE‡(CASCI) = E(TS) − E(reactant); меняет ли квант ПОРЯДОК барьеров Cu/CuAl?

Ограничения: жёсткая модель (замороженный металл-кластер, def2-SVP, нет constant-potential,
драйвер = C–C по x). Это оценка барьера и его мультиреференс-поправки, НЕ фарадеевская
эффективность. Полный constant-potential + 44q-VQE остаются за AWS.

Запуск: python3 -u routes/co2_to_fuels_ts.py
"""
import os, sys, json, time
import numpy as np
from scipy.optimize import minimize
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from co2_to_fuels_sites import METALS, _mol, HARTREE_EV, RESULTS
from pyscf import dft, mcscf
from pyscf.scf import addons
from pyscf.mcscf import avas

DS = [1.8, 2.2, 2.6, 3.0, 3.4, 3.8, 4.2]   # driven C···C (x-separation), Å


def _scf(mol, dm0=None):
    """RI-PBE, Newton(SOSCF)-primary — robust+fast for the metal-adsorbate clusters
    (≈40× faster than DIIS+level-shift on the hard CuAl geometries). Fermi-smearing
    fallback breaks any residual oscillation, then re-tightens with Newton."""
    mf = dft.RKS(mol).density_fit(); mf.xc = "pbe"; mf.verbose = 0; mf.conv_tol = 1e-8
    mfn = mf.newton()
    try:
        mfn.kernel(dm0=dm0)
    except TypeError:
        mfn.kernel()
    if not mfn.converged:
        mf2 = dft.RKS(mol).density_fit(); mf2.xc = "pbe"; mf2.verbose = 0
        mf2 = addons.smearing_(mf2, sigma=0.01, method="fermi")
        mf2.max_cycle = 200; mf2.conv_tol = 1e-6; mf2.kernel(dm0=dm0)
        mfn = mf2.newton(); mfn.kernel(mf2.make_rdm1())
    return mfn


def casci_pi(metal, ads):
    """CASCI(π-manifold, 12e,10o via AVAS) + NOON n_u at a geometry (fast SCF)."""
    mf = _scf(_mol(metal, ads))
    no, ne, orbs = avas.avas(mf, ["C 2p", "O 2px", "O 2py"], threshold=0.4, verbose=0)
    mc = mcscf.CASCI(mf, no, ne); mc.verbose = 0; mc.kernel(orbs)
    dm = mc.fcisolver.make_rdm1(mc.ci, mc.ncas, mc.nelecas)
    n = np.linalg.eigvalsh(dm)[::-1]
    nu = float(np.sum(np.minimum(n, 2 - n)))
    return float(mf.e_tot), float(mc.e_tot), round(nu, 4), (int(ne), int(no))


def build_ads(d, free):
    """free = [zC1, zC2, O1x,O1y,O1z, O2x,O2y,O2z]; carbons driven at x=±d/2, y=0."""
    zC1, zC2, o1x, o1y, o1z, o2x, o2y, o2z = free
    return np.array([[-d / 2, 0, zC1], [o1x, o1y, o1z],
                     [d / 2, 0, zC2], [o2x, o2y, o2z]])


def _eg(metal, d, free, dm0):
    mf = _scf(_mol(metal, build_ads(d, free)), dm0)
    g = mf.nuc_grad_method().kernel()
    nm = len(metal)                                   # atoms: metal..., C1,O1,C2,O2
    gf = np.array([g[nm][2], g[nm + 2][2], *g[nm + 1], *g[nm + 3]])  # free-DOF gradient
    return mf.e_tot, gf, mf.make_rdm1()


def relax_at_d(metal, d, free0, dm0, maxiter=18):
    cache = {"dm": dm0}

    def fun(x):
        e, gf, dm = _eg(metal, d, x, cache["dm"]); cache["dm"] = dm
        return e, gf

    res = minimize(fun, free0, jac=True, method="L-BFGS-B",
                   options={"maxiter": maxiter, "maxfun": 2 * maxiter,
                            "ftol": 1e-7, "gtol": 5e-4})
    return res.x, float(res.fun), cache["dm"]


def _persist(res):
    json.dump(res, open(RESULTS, "w"), ensure_ascii=False, indent=1)


def run_site_ts(tag, res):
    metal = METALS[tag]
    print(f"\n--- TS scan site {tag} ---")
    site = res.setdefault("ts_barrier", {}).setdefault(tag, {})
    prof = site.get("profile", [])
    done = {round(p["d"], 2): p for p in prof}
    free = np.array(prof[-1]["free"]) if prof else None
    dm = None
    for d in DS:                                       # per-point resume
        if round(d, 2) in done:
            free = np.array(done[round(d, 2)]["free"]); continue
        t0 = time.time()
        f0 = (free.copy() if free is not None
              else np.array([1.85, 1.85, -d / 2, 0., 3.0, d / 2, 0., 3.0]))
        f0[2], f0[5] = -d / 2, d / 2                   # re-seed O-x near new carbon x
        free, e, dm = relax_at_d(metal, d, f0, dm)
        prof.append({"d": float(d), "e_pbe": round(float(e), 6),
                     "free": [round(float(x), 3) for x in free]})
        site["profile"] = prof; _persist(res)          # checkpoint after each point
        print(f"  d={d:.2f}Å  E_pbe={e:.5f}  ({time.time()-t0:.0f}s)")

    prof = sorted(prof, key=lambda p: p["d"]); site["profile"] = prof
    es = np.array([p["e_pbe"] for p in prof])
    i_react = int(np.argmin(es))                       # reactant = global min (separated well)
    lo, hi = (0, i_react) if i_react > 0 else (0, len(es) - 1)
    i_ts = lo + int(np.argmax(es[lo:hi + 1]))          # TS = max between coupled end & reactant
    react, ts = prof[i_react], prof[i_ts]
    barrier_exists = (es[i_ts] - es[i_react]) > 0.02 / HARTREE_EV
    print(f"  reactant @ d={react['d']:.2f};  TS @ d={ts['d']:.2f}  barrier_exists={barrier_exists}")
    site.update(d_reactant=react["d"], d_ts=ts["d"], barrier_exists=bool(barrier_exists))

    for label, p in [("reactant", react), ("ts", ts)]:
        ads = build_ads(p["d"], np.array(p["free"]))
        e_pbe, e_casci, nu, cas = casci_pi(metal, ads)
        site[label + "_e_pbe"] = round(e_pbe, 6)
        site[label + "_e_casci"] = round(e_casci, 6)
        site[label + "_nu_pi"] = nu
        _persist(res)
        print(f"  {label:8s} d={p['d']:.2f}  E_casci={e_casci:.5f}  n_u(π)={nu:.3f}")
    site["barrier_pbe_eV"] = round((site["ts_e_pbe"] - site["reactant_e_pbe"]) * HARTREE_EV, 3)
    site["barrier_casci_eV"] = round((site["ts_e_casci"] - site["reactant_e_casci"]) * HARTREE_EV, 3)
    site["quantum_corr_eV"] = round(site["barrier_casci_eV"] - site["barrier_pbe_eV"], 3)
    _persist(res)
    print(f"  ΔE‡: PBE={site['barrier_pbe_eV']:+.3f}  CASCI={site['barrier_casci_eV']:+.3f}  "
          f"Δ_corr={site['quantum_corr_eV']:+.3f} eV   TS n_u(π)={site['ts_nu_pi']:.3f}")
    return site


if __name__ == "__main__":
    res = json.load(open(RESULTS)) if os.path.exists(RESULTS) else {}
    tb = res.setdefault("ts_barrier", {})
    for tag in ("Cu2", "CuAl"):
        if tb.get(tag, {}).get("barrier_pbe_eV") is not None:
            print(f"--- TS {tag}: already complete, skipping ---"); continue
        run_site_ts(tag, res)
        print(f"  saved TS {tag} → results")
    if all(tb.get(t, {}).get("barrier_pbe_eV") is not None for t in ("Cu2", "CuAl")):
        po = sorted(("Cu2", "CuAl"), key=lambda t: tb[t]["barrier_pbe_eV"])
        co = sorted(("Cu2", "CuAl"), key=lambda t: tb[t]["barrier_casci_eV"])
        tb["_summary"] = dict(pbe_order=list(po), casci_order=list(co),
                              order_flipped=bool(list(po) != list(co)),
                              note=("ΔE‡ from constrained relaxed C···C scan (frozen-metal "
                                    "cluster, def2-SVP, no constant-potential). Estimate of the "
                                    "coupling barrier + its multireference correction; NOT Faradaic "
                                    "efficiency. Constant-potential + 44q VQE remain AWS."))
        _persist(res)
        print(f"\n=== BARRIER ORDERING ===\n  PBE:{po}  CASCI:{co}  "
              f"flipped={tb['_summary']['order_flipped']}")
    print("\nsaved → results;  done.")
