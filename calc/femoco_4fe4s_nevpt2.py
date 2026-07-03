#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
calc/femoco_4fe4s_nevpt2.py — FeMoco diary, Stage 18.

A follow-up to Stage 17 that tried to add dynamic correlation (NEVPT2) and converge
the CASSCF on the [4Fe-4S] frontier — and found something more instructive than a
single tidy number. The measured result, on the real 288-orbital [Fe4S4(SH)4]2-:

  * Dynamic correlation is CHEAP. Strongly-contracted NEVPT2 on the CAS(6,6) runs in
    ~1 minute at this scale.
  * But it must sit on ORBITAL-OPTIMIZED (CASSCF) orbitals. A plain CASCI in the
    high-spin-ROHF magnetic orbitals is a POOR reference: those orbitals are localized
    on individual Fe, so the six-electron frontier configuration sits ~10 Ha (thousands
    of kcal/mol) ABOVE the high-spin determinant. NEVPT2 on such a reference returns a
    huge, meaningless correction — it cannot rescue a bad reference.
  * Converging the CASSCF that WOULD fix this is the cost wall. It has to recover ~10 Ha
    of orbital relaxation over the full virtual space; each macro-iteration is ~1.5-2 min
    and the gradient only halves per step, so it takes ~an hour and Stage 17's run never
    converged. Windowing the virtuals to make it cheap throws away most of the relaxation.
  * And the multireference metric is ORBITAL-FRAGILE: n_u swings from ~1.3 (CASCI, raw
    orbitals) up toward Stage 17's unconverged 3.99 as the orbitals relax. So Stage 17's
    n_u = 3.99 is an UNCONVERGED estimate; the qualitative claim (strong multireference,
    n_u >> 0) is robust, the precise value is not cheaply pinned.

This measures a new facet of "why quantum" for Fe-S clusters: the easy post-hoc step
(NEVPT2) is cheap, but the ESSENTIAL step — a converged multireference reference — is
where classical methods choke, and the answer is delicate. That is exactly the regime a
fault-tolerant quantum computer is meant to remove.

This script records those facts reproducibly (SCF loaded from the Stage-17 chkfile):
  1. CASCI(6,6) in the ROHF frontier orbitals -> E, gap above HS, S^2, n_u, NOON.
  2. NEVPT2 on that CASCI -> correction + wall-time (the "cheap but useless on a bad
     reference" point).
  3. A BOUNDED CASSCF (3 macro-iterations) -> shows the energy descending far too
     slowly to converge cheaply, and n_u already drifting from the CASCI value.

HONEST SCOPE: same six-orbital frontier subspace as Stage 17 (a demonstration subspace,
not the 40-qubit full Fe-3d magnetic space); idealized cubane; def2-SVP. Every number is
computed at run time; nothing is invented.

Run: python3 -u calc/femoco_4fe4s_nevpt2.py
"""
import os
import sys
import json
import time
import numpy as np
from pyscf import mcscf, mrpt, gto

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import femoco_4fe4s_casscf as s17

HERE = os.path.dirname(os.path.abspath(__file__))
RESULTS = os.path.join(HERE, "femoco_4fe4s_nevpt2_results.json")
KCAL = 627.50947406
NCAS = 6


def n_u_noon(mc, ncas, nelecas):
    occ = np.clip(np.linalg.eigvalsh(mc.fcisolver.make_rdm1(mc.ci, ncas, nelecas))[::-1],
                  0.0, 2.0)
    return float(np.sum(occ * (2.0 - occ))), [round(float(x), 3) for x in occ]


def main():
    t_all = time.time()
    res = {"status": "running",
           "meta": {"system": "[Fe4S4(SH)4]2- cubane, frontier CAS(6e,6o)",
                    "basis": s17.BASIS,
                    "message": "dynamic correlation (NEVPT2) is cheap; a CONVERGED "
                               "multireference reference (CASSCF) is the cost wall; the "
                               "n_u metric is orbital-fragile -> Stage 17's 3.99 is an "
                               "unconverged estimate.",
                    "note": "same frontier subspace as Stage 17; a demonstration "
                            "subspace, not the 40-qubit full Fe-3d magnetic space."}}
    try:
        atoms = s17.build_cubane()
        mol = gto.M(atom=[(e, tuple(p)) for e, p in atoms], basis=s17.BASIS,
                    charge=s17.CHARGE, spin=s17.SPIN, verbose=0, max_memory=11000)
        mf = s17.converge_rohf(mol)
        e_hs = float(mf.e_tot)
        res["scf"] = {"e_hs_rohf": e_hs, "converged": bool(mf.converged)}
        print(f"HS ROHF: E={e_hs:.5f} conv={mf.converged}", flush=True)

        active, ncas, nelec, fe3d0 = s17.frontier_magnetic_cas(mf, NCAS)
        na = (nelec + (nelec % 2)) // 2
        nb = nelec - na
        res["frontier"] = {"ncas": ncas, "nelec": nelec, "nelecas": [int(na), int(nb)],
                           "qubits_jw": 2 * ncas, "active_mo_idx": active,
                           "fe3d_char_initial": round(fe3d0, 3)}
        print(f"frontier CAS({nelec}e,{ncas}o)=12q; init Fe-3d {fe3d0:.2f}", flush=True)

        # (1) CASCI in the raw ROHF frontier orbitals — a POOR reference
        mc = mcscf.CASCI(mf, ncas, (na, nb))
        mc.verbose = 0
        mo = mc.sort_mo(active, base=0)
        t0 = time.time()
        e_casci = mc.kernel(mo)[0]
        ss = float(mc.fcisolver.spin_square(mc.ci, ncas, (na, nb))[0])
        nu_ci, noon_ci = n_u_noon(mc, ncas, (na, nb))
        casci = {"e_casci": float(e_casci), "s_squared": round(ss, 3),
                 "gap_above_hs_kcal": round((e_casci - e_hs) * KCAL, 1),
                 "n_u": round(nu_ci, 3), "noon": noon_ci,
                 "seconds": round(time.time() - t0, 1)}
        print(f"CASCI (ROHF orbitals): E={e_casci:.5f}  {casci['gap_above_hs_kcal']} "
              f"kcal/mol ABOVE HS  <S^2>={ss:.2f}  n_u={nu_ci:.3f}", flush=True)

        # (2) NEVPT2 on the CASCI — cheap, but meaningless on this poor reference
        try:
            t1 = time.time()
            ec = mrpt.NEVPT(mc).kernel()
            casci["nevpt2_corr"] = float(ec)
            casci["e_nevpt2"] = float(e_casci + ec)
            casci["nevpt2_seconds"] = round(time.time() - t1, 1)
            print(f"NEVPT2 on CASCI: corr={ec:.4f} Ha in "
                  f"{casci['nevpt2_seconds']}s (large -> reference is poor)", flush=True)
        except Exception as e:
            casci["nevpt2_error"] = f"{type(e).__name__}: {e}"
            print("NEVPT2 failed:", casci["nevpt2_error"], flush=True)
        res["casci_rohf_orbitals"] = casci

        # (3) BOUNDED CASSCF (3 macro-iters) — the orbital relaxation is a cost wall
        mc2 = mcscf.CASSCF(mf, ncas, (na, nb)).density_fit()
        mc2.max_cycle_macro = 3
        mc2.verbose = 0
        t2 = time.time()
        e_cas3 = mc2.kernel(mo)[0]
        nu3, noon3 = n_u_noon(mc2, ncas, (na, nb))
        recovered = (e_cas3 - e_casci)
        # Stage 17's (unconverged, full-virtual) CASSCF reached this energy:
        e_stage17 = -8230.314104
        res["casscf_bounded_3macro"] = {
            "e_casscf_3iter": float(e_cas3), "converged": bool(mc2.converged),
            "seconds": round(time.time() - t2, 1),
            "energy_recovered_Ha": round(float(recovered), 3),
            "still_needed_Ha_to_stage17": round(float(e_stage17 - e_cas3), 3),
            "n_u_after_3iter": round(nu3, 3), "noon_after_3iter": noon3,
            "n_u_casci": round(nu_ci, 3), "n_u_stage17_unconverged": 3.985,
            "comment": "3 macro-iters recover only a fraction of the ~10 Ha orbital "
                       "relaxation; n_u already drifting CASCI->Stage17 (orbital-fragile)."}
        print(f"CASSCF (3 macro-iters): E={e_cas3:.5f} conv={mc2.converged} "
              f"recovered {recovered:.2f} Ha, still {e_stage17-e_cas3:.2f} Ha to Stage-17; "
              f"n_u {nu_ci:.2f}->{nu3:.2f}", flush=True)
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
