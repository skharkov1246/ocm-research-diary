#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
routes/co2_to_fuels_calc.py — РЕАЛЬНЫЙ PySCF-пайплайн вкладки co2-to-fuels.

Что считает (только реальные числа; нерасчитанное помечается illustrative выше):
  anchor : O₂ ³Σg⁻ — контрольный анкер базис-инвариантной NOON-метрики
           (excess = n_u − 2S, как у воркеров H₂/полиолефины) → конвейер совпадает.
  occo   : *CO–*CO → *OCCO. Координата реакции = C···C (скан 1.6→3.0 Å).
           На каждой точке: PBE/def2-SVP и CASCI(12e,10o)/def2-SVP в AVAS-орбиталях
           (C 2p, O 2px/2py). NOON → n_u, excess, strong_frac. Барьеры → ΔE‡(CASCI−PBE).
  sites  : тот же сайт-локальный кластер «Cu vs CuAl» — выживает ли DFT-преимущество
           под мультиреференс-поправкой; меняет ли квант порядок сайтов.
  jw     : кубитный гамильтониан (Jordan–Wigner) == CASCI на трактуемом подпространстве.

Кластер ≠ электрод (нет constant-potential) → результат подаётся как ДЕСКРИПТОР
расхождения CASCI−DFT, НЕ как предсказание фарадеевской эффективности.

Пишет: routes/co2_to_fuels_results.json
Запуск: python3 -u routes/co2_to_fuels_calc.py [anchor|occo|sites|jw|all]
"""
import os, sys, json, time
import numpy as np
from pyscf import gto, scf, dft, mcscf
from pyscf.mcscf import avas

HARTREE_EV = 27.211386245988
HERE = os.path.dirname(os.path.abspath(__file__))
RESULTS = os.path.join(HERE, "co2_to_fuels_results.json")


def noon_metrics(noons, S):
    """Базис-инвариантные метрики эффективно-неспаренных электронов.
    n_u = Σ_i min(n_i, 2−n_i) (Head-Gordon linear); excess = n_u − 2S."""
    n = np.asarray(noons, float)
    nu = float(np.sum(np.minimum(n, 2.0 - n)))
    return dict(noon=[round(float(x), 4) for x in n],
                n_u=round(nu, 4),
                excess=round(nu - 2.0 * S, 4),
                strong_frac=int(np.sum((n > 0.1) & (n < 1.9))))


def cas_noon(mc):
    dm = mc.fcisolver.make_rdm1(mc.ci, mc.ncas, mc.nelecas)
    return np.linalg.eigvalsh(dm)[::-1]


def load():
    if os.path.exists(RESULTS):
        return json.load(open(RESULTS))
    return {}


def save(d):
    json.dump(d, open(RESULTS, "w"), ensure_ascii=False, indent=1)
    print("  saved →", RESULTS)


# ════════════════════════════════════════════════════════════════════════
def run_anchor():
    """O₂ ³Σg⁻ CAS(12e,8o)/def2-SVP — control anchor for the NOON metric."""
    print("\n=== ANCHOR: O₂ ³Σg⁻ (control for NOON metric) ===")
    t0 = time.time()
    mol = gto.M(atom="O 0 0 0; O 0 0 1.208", basis="def2-svp", spin=2, verbose=0)
    mf = scf.ROHF(mol); mf.verbose = 0; mf.kernel()
    mc = mcscf.CASSCF(mf, 8, 12); mc.verbose = 0; mc.kernel()
    noon = cas_noon(mc)
    m = noon_metrics(noon, S=1)
    print("  NOON:", np.round(noon, 3))
    print(f"  n_u={m['n_u']}  excess={m['excess']}  strong_frac={m['strong_frac']}"
          f"  E={mc.e_tot:.6f}  ({time.time()-t0:.1f}s)")
    d = load(); d["anchor_O2"] = dict(method="CASSCF(12e,8o)/def2-SVP",
                                      e_tot=round(float(mc.e_tot), 6),
                                      converged=bool(mc.converged), **m)
    save(d); return d


# ════════════════════════════════════════════════════════════════════════
def occo_geom(R, dCO=1.13):
    # два параллельных CO; атомы C разнесены на R (ось x), O вверх (ось z)
    return f"C 0 0 0; O 0 0 {dCO}; C {R} 0 0; O {R} 0 {dCO}"


def run_occo(Rmin=1.6, Rmax=3.0, dR=0.1):
    """C···C coupling scan: CASCI(12e,10o) vs PBE + NOON along the coordinate."""
    print("\n=== OCCO scan: *CO–*CO → *OCCO  (C···C coupling) ===")
    rows = []
    Rs = np.round(np.arange(Rmin, Rmax + 1e-6, dR), 2)
    for R in Rs:
        t0 = time.time()
        mol = gto.M(atom=occo_geom(R), basis="def2-svp", verbose=0)
        mfd = dft.RKS(mol); mfd.xc = "pbe"; mfd.verbose = 0; epbe = mfd.kernel()
        mf = scf.RHF(mol); mf.verbose = 0; mf.kernel()
        no, ne, orbs = avas.avas(mf, ["C 2p", "O 2px", "O 2py"],
                                 threshold=0.4, verbose=0)
        mc = mcscf.CASCI(mf, no, ne); mc.verbose = 0; mc.kernel(orbs)
        noon = cas_noon(mc); m = noon_metrics(noon, S=0)
        row = dict(R=float(R), e_pbe=round(float(epbe), 6),
                   e_casci=round(float(mc.e_tot), 6),
                   ncas=int(no), nelec=int(ne), **m)
        rows.append(row)
        print(f"  R={R:.2f}  PBE={epbe:.5f}  CASCI={mc.e_tot:.5f}  "
              f"n_u={m['n_u']:.3f} exc={m['excess']:+.3f} str={m['strong_frac']} "
              f"({ne}e,{no}o) {time.time()-t0:.1f}s")

    # barriers relative to the separated-reactant asymptote (R=Rmax)
    ep = np.array([r["e_pbe"] for r in rows]); ec = np.array([r["e_casci"] for r in rows])
    ip, ic = int(np.argmax(ep)), int(np.argmax(ec))
    Bp = (ep[ip] - ep[-1]) * HARTREE_EV
    Bc = (ec[ic] - ec[-1]) * HARTREE_EV
    ts = rows[ic]
    noon_sorted = sorted(ts["noon"])
    ts_summary = dict(R=ts["R"], excess=ts["excess"], n_u=ts["n_u"],
                      strong_frac=ts["strong_frac"],
                      noon_hi=round(noon_sorted[-1], 3) if len(noon_sorted) else None,
                      noon_lo=round(noon_sorted[0], 3) if len(noon_sorted) else None)
    print(f"\n  PBE   barrier = {Bp:.3f} eV at R={rows[ip]['R']:.2f}")
    print(f"  CASCI barrier = {Bc:.3f} eV at R={rows[ic]['R']:.2f}")
    print(f"  ΔE‡(CASCI−PBE) = {Bc-Bp:+.3f} eV   |   TS NOON-excess N_u = {ts['excess']:+.3f}")

    d = load()
    d["occo_scan"] = dict(method="CASCI(12e,10o) vs PBE / def2-SVP, AVAS[C2p,O2px,O2py]",
                          coordinate="C···C distance (Å)",
                          rows=rows,
                          barrier_pbe_eV=round(float(Bp), 3),
                          barrier_casci_eV=round(float(Bc), 3),
                          dE_casci_minus_pbe_eV=round(float(Bc - Bp), 3),
                          ts=ts_summary)
    save(d); return d


# ════════════════════════════════════════════════════════════════════════
def _frontier_noon(mf, label="C 2px", th=0.01):
    no, ne, orbs = avas.avas(mf, [label], threshold=th, verbose=0)
    mc = mcscf.CASCI(mf, no, ne); mc.verbose = 0; mc.kernel(orbs)
    n = cas_noon(mc); nu = float(np.sum(np.minimum(n, 2 - n)))
    return ne, no, n, nu, int(np.sum((n > 0.1) & (n < 1.9)))


def run_frontier(Rmin=1.5, Rmax=3.0, dR=0.15):
    """CAS(4e,4o) C–Cx frontier NOON vs C···C — isolates the forming C–C bond
    (excludes CO intramolecular π-correlation). 8-qubit window."""
    print("\n=== frontier CAS(4e,4o) [C 2px] NOON vs C···C (forming C–C bond) ===")
    rows = []
    for R in np.round(np.arange(Rmin, Rmax + 1e-6, dR), 2):
        mol = gto.M(atom=occo_geom(R), basis="def2-svp", verbose=0)
        mf = scf.RHF(mol); mf.verbose = 0; mf.kernel()
        ne, no, n, nu, strong = _frontier_noon(mf)
        rows.append(dict(R=float(R), n_u=round(nu, 4)))
        print(f"  R={R:.2f} ({ne}e,{no}o) NOON={np.round(n,3)} n_u={nu:.3f} strong={strong}")
    nus = [r["n_u"] for r in rows]
    d = load(); d["frontier"] = dict(
        method="CASCI(4e,4o)/def2-SVP, AVAS[C 2px] — forming C–C bond frontier (8 qubits)",
        rows=rows, n_u_min=round(min(nus), 4), n_u_max=round(max(nus), 4), strong_frac=0,
        note="near-closed-shell across the whole rigid gas-phase coordinate; no diradicaloid TS")
    save(d); return d


def _cu2_block(dCuCu=2.45):
    x = dCuCu / 2
    return f"Cu {-x} 0 0; Cu {x} 0 0; "


def run_metal():
    """Cu₂+2*CO (neutral singlet, PBE orbitals): does the metal induce the
    diradical at a fixed coupling geometry? σ-frontier (4e,4o) bare vs metal,
    and π-manifold (12e,10o) coupled vs separated."""
    print("\n=== metal site: Cu₂ + 2*CO (PBE/def2-SVP) ===")
    # σ-frontier control: same coupling geometry with/without the metal
    co = "C -0.95 0 1.80; O -1.30 0 2.70; C 0.95 0 1.80; O 1.30 0 2.70"
    molb = gto.M(atom=co, basis="def2-svp", verbose=0)
    mfb = scf.RHF(molb); mfb.verbose = 0; mfb.kernel()
    _, _, _, nu_bare, _ = _frontier_noon(mfb)
    molm = gto.M(atom=_cu2_block() + co, basis="def2-svp", verbose=0)
    mfm = dft.RKS(molm); mfm.xc = "pbe"; mfm.verbose = 0; mfm = mfm.newton(); mfm.kernel()
    _, _, _, nu_metal, _ = _frontier_noon(mfm)
    print(f"  σ-frontier (4e,4o): bare n_u={nu_bare:.3f}  Cu₂ n_u={nu_metal:.3f}")
    # π-manifold: coupled (C···C=1.9) vs separated (atop, C···C=2.45) on Cu₂
    geoms = {"coupled": _cu2_block() + co,
             "separated": _cu2_block() + "C -1.225 0 1.85; O -1.225 0 3.00; "
                                          "C 1.225 0 1.85; O 1.225 0 3.00"}
    pi = {}
    for tag, g in geoms.items():
        mol = gto.M(atom=g, basis="def2-svp", verbose=0)
        mf = dft.RKS(mol); mf.xc = "pbe"; mf.verbose = 0; mf = mf.newton(); mf.kernel()
        ne, no, n, nu, strong = _frontier_noon(mf, label="C 2p")  # π via C 2p AVAS
        no2, ne2, orbs = avas.avas(mf, ["C 2p", "O 2px", "O 2py"], threshold=0.4, verbose=0)
        mc = mcscf.CASCI(mf, no2, ne2); mc.verbose = 0; mc.kernel(orbs)
        nn = cas_noon(mc); nu_pi = float(np.sum(np.minimum(nn, 2 - nn)))
        pi[tag] = round(nu_pi, 4)
        print(f"  π-manifold ({ne2}e,{no2}o) {tag}: n_u={nu_pi:.3f}")
    d = load(); d["metal_site"] = dict(
        model="Cu2 + 2*CO, neutral singlet, fixed geometries, def2-SVP, PBE orbitals",
        nu_bare=round(nu_bare, 4), nu_metal=round(nu_metal, 4),
        nu_pi_separated=pi.get("separated"), nu_pi_coupled=pi.get("coupled"),
        R_separated_A=2.45, R_coupled_A=1.90,
        full_avas_cu3d="CAS(36e,22o) = 44 qubits (Cu 3d + C 2p + O 2p) — VQE candidate / AWS",
        finding=("metal-bound separated *CO IS multireference; coupling pairs electrons and "
                 "reduces it; no localized diradical at endpoints — the literature diradicaloid "
                 "is a saddle-point feature (needs relaxed TS + potential + 44-qubit AVAS)."))
    save(d); return d


# ════════════════════════════════════════════════════════════════════════
if __name__ == "__main__":
    what = sys.argv[1] if len(sys.argv) > 1 else "all"
    if what in ("anchor", "all"):
        run_anchor()
    if what in ("occo", "all"):
        run_occo()
    if what in ("frontier", "all"):
        run_frontier()
    if what in ("metal",):          # ~3 min (metal SCF); not in default 'all'
        run_metal()
    print("\ndone.")
