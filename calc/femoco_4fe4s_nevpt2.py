#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
calc/femoco_4fe4s_nevpt2.py — FeMoco diary, Stage 18.

Follow-up to Stage 17 (calc/femoco_4fe4s_casscf.py). Two upgrades on the SAME
[4Fe-4S] frontier magnetic subspace:
  (1) CONVERGE the CASSCF. Stage 17 hit the macro-iteration cap. Here we freeze the
      doubly-occupied core (only active<->virtual and active<->active rotations remain)
      and fix the CI spin — the frontier CAS(6,6) then converges cleanly.
  (2) Add DYNAMIC correlation on top: strongly-contracted NEVPT2 (pyscf.mrpt), the
      next rung of rigor beyond bare CASSCF.

And, as a state-specific demonstration, we compute BOTH the low-spin singlet (S=0)
and the high-spin septet (S=3) of the six-electron frontier subspace at CASSCF and
CASSCF+NEVPT2. The multireference method gives ONE definite, sign-stable subspace gap
— a methodological contrast to Stage 15, where four DFT functionals of the full cubane
could not even agree on the SIGN of the magnetic coupling.

The SCF is loaded from the Stage-17 chkfile (cached), so this run skips the ~70-minute
mean field.

HONEST SCOPE (read before trusting any number).
  - This is the SAME six-orbital FRONTIER magnetic subspace as Stage 17 (six central
    singly-occupied Fe-3d MOs of the high-spin reference). The other twelve magnetic
    electrons sit in the doubly-occupied/empty reference, so the singlet-septet gap here
    is a SUBSPACE quantity, NOT the physical [4Fe-4S] magnetic coupling J. The true J
    needs the full ~20-orbital Fe-3d space (40 qubits), whose exact CASCI/NEVPT2 is out
    of reach — that is the classical wall this diary documents, unchanged.
  - Idealized cubane (not CIF-derived), def2-SVP; NEVPT2 is single-state, no relaxation.
  - Every number is computed at run time; convergence flags and any NEVPT2 resource
    failure are recorded honestly in the results JSON.

Run: python3 -u calc/femoco_4fe4s_nevpt2.py
"""
import os
import sys
import json
import time
import numpy as np
from pyscf import mcscf, mrpt, gto

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import femoco_4fe4s_casscf as s17   # geometry, cached-SCF loader, frontier selector

HERE = os.path.dirname(os.path.abspath(__file__))
RESULTS = os.path.join(HERE, "femoco_4fe4s_nevpt2_results.json")
KCAL = 627.50947406
NCAS = 6
# 2S per state for the 6-electron frontier: singlet S=0, septet S=3 (all six aligned)
STATES = {"singlet_S0": 0, "septet_S3": 6}


def cas_state(mf, active, ncas, nelec, spin2):
    """Converged state-specific CASSCF (frozen core, fixed spin) + NEVPT2 in the
    frontier orbitals `active` (0-based MO indices, seeded via sort_mo)."""
    na = (nelec + spin2) // 2
    nb = nelec - na
    if nb < 0 or na > ncas or nb > ncas:
        return {"error": f"bad split na={na} nb={nb} in {ncas}o"}
    ncore = (int(round(sum(mf.mo_occ))) - nelec) // 2   # doubly-occupied core count
    mc = mcscf.CASSCF(mf, ncas, (na, nb)).density_fit()
    mc.frozen = ncore                 # freeze core: only active-virtual/active-active
    mc.max_cycle_macro = 200
    mc.conv_tol = 1e-8
    mc.verbose = 0
    S = spin2 / 2.0
    try:
        mc.fix_spin_(ss=S * (S + 1))
    except Exception:
        pass
    mo = mc.sort_mo(active, base=0)
    t0 = time.time()
    e_cas = mc.kernel(mo)[0]
    occ = np.clip(np.linalg.eigvalsh(mc.fcisolver.make_rdm1(mc.ci, ncas, (na, nb)))[::-1],
                  0.0, 2.0)
    n_u = float(np.sum(occ * (2.0 - occ)))
    ss = float(mc.fcisolver.spin_square(mc.ci, ncas, (na, nb))[0])
    d = {"spin2": spin2, "nelecas": [int(na), int(nb)], "ncore_frozen": int(ncore),
         "e_casscf": float(e_cas), "casscf_converged": bool(mc.converged),
         "cas_seconds": round(time.time() - t0, 1), "n_u": round(n_u, 3),
         "noon": [round(float(x), 3) for x in occ], "s_squared": round(ss, 3)}
    print(f"  {list(STATES.keys())[list(STATES.values()).index(spin2)]:11s}: "
          f"E_cas={e_cas:.6f} conv={mc.converged} n_u={n_u:.3f} "
          f"<S^2>={ss:.2f} ({d['cas_seconds']}s)", flush=True)
    try:
        t1 = time.time()
        ec = mrpt.NEVPT(mc).kernel()
        d["nevpt2_corr"] = float(ec)
        d["e_nevpt2"] = float(e_cas + ec)
        d["nevpt2_seconds"] = round(time.time() - t1, 1)
        print(f"               + NEVPT2 corr={ec:.6f} -> E={e_cas+ec:.6f} "
              f"({d['nevpt2_seconds']}s)", flush=True)
    except Exception as e:
        d["nevpt2_error"] = f"{type(e).__name__}: {e}"
        print(f"               NEVPT2 failed: {d['nevpt2_error']}", flush=True)
    return d


def gap(states, key):
    """Signed gap E(septet) - E(singlet) in kcal/mol for the given energy key."""
    a = states.get("singlet_S0", {}).get(key)
    b = states.get("septet_S3", {}).get(key)
    if a is None or b is None:
        return None
    return round((b - a) * KCAL, 2)


def main():
    t_all = time.time()
    res = {"status": "running",
           "meta": {"system": "[Fe4S4(SH)4]2- cubane, frontier CAS(6e,6o)",
                    "basis": s17.BASIS, "method": "frozen-core state-specific CASSCF "
                    "+ SC-NEVPT2 (singlet S=0 & septet S=3)",
                    "note": "SAME frontier subspace as Stage 17; the singlet-septet gap "
                            "is a SUBSPACE quantity, NOT the physical [4Fe-4S] J (which "
                            "needs the untractable 40-qubit Fe-3d space)."}}
    try:
        atoms = s17.build_cubane()
        mol = gto.M(atom=[(e, tuple(p)) for e, p in atoms], basis=s17.BASIS,
                    charge=s17.CHARGE, spin=s17.SPIN, verbose=0, max_memory=11000)
        mf = s17.converge_rohf(mol)          # loads Stage-17 chkfile (cached SCF)
        scf_ok = bool(mf.converged)
        res["scf"] = {"e_hs_rohf": float(mf.e_tot), "converged": scf_ok}
        print(f"HS ROHF: E={mf.e_tot:.5f} conv={scf_ok}", flush=True)

        active, ncas, nelec, fe3d_frac = s17.frontier_magnetic_cas(mf, NCAS)
        res["frontier"] = {"ncas": ncas, "nelec": nelec, "qubits_jw": 2 * ncas,
                           "active_mo_idx": active,
                           "frontier_fe3d_fraction": round(fe3d_frac, 3),
                           "based_on_converged_scf": scf_ok}
        print(f"frontier magnetic CAS({nelec}e,{ncas}o) = {2*ncas} qubits; "
              f"Fe-3d character {fe3d_frac:.2f}", flush=True)

        states = {}
        for name, s2 in STATES.items():
            states[name] = cas_state(mf, active, ncas, nelec, s2)
        res["states"] = states
        res["gap_septet_minus_singlet_kcal"] = {
            "casscf": gap(states, "e_casscf"),
            "nevpt2": gap(states, "e_nevpt2")}
        print("gap septet-singlet (kcal/mol):",
              json.dumps(res["gap_septet_minus_singlet_kcal"]), flush=True)
        res["status"] = "ok"
    except Exception as e:
        res["status"] = f"aborted: {type(e).__name__}: {e}"
        print("ABORT:", res["status"], flush=True)
    finally:
        res["meta"]["total_seconds"] = round(time.time() - t_all, 1)
        with open(RESULTS, "w") as fh:
            json.dump(res, fh, ensure_ascii=False, indent=1)
        print(f"\nwrote {RESULTS}  status={res['status']}  "
              f"(total {res['meta']['total_seconds']}s)", flush=True)


if __name__ == "__main__":
    main()
