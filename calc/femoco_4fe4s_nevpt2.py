#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
calc/femoco_4fe4s_nevpt2.py — FeMoco diary, Stage 18.

Two upgrades on the Stage-17 [4Fe-4S] frontier magnetic subspace:
  (1) CONVERGE the CASSCF. Stage 17 reported n_u from a CASSCF that hit the
      macro-iteration cap (unconverged). Full orbital optimization on the real
      288-orbital cluster is the cost wall (~1 h and still not converged). Here we
      make it tractable by WINDOWING the orbital optimization: freeze the
      doubly-occupied core AND the deep virtuals, letting the six active Fe-3d
      orbitals mix only with a window of NVIRT frontier virtuals. That converges in
      minutes and — crucially — keeps the active space from drifting off the Fe-3d
      manifold (a risk in the unfrozen unconverged run).
  (2) Add DYNAMIC correlation: strongly-contracted NEVPT2 (pyscf.mrpt) on the
      converged CASSCF (NEVPT2 itself is cheap here, ~1 min).

HONEST CHECK. We report the Löwdin Fe-3d character of the FINAL (converged) active
orbitals, and compare n_u to Stage 17's unconverged 3.99 — an explicit
does-the-number-survive-convergence test, in the spirit of the diary.

HONEST SCOPE.
  - SAME six-orbital FRONTIER magnetic subspace as Stage 17 (six central singly-
    occupied Fe-3d MOs of the HS reference); a demonstration subspace, not the full
    magnetic space (Fe-3d = 40 qubits, the untractable wall).
  - Orbital optimization is WINDOWED (core + deep virtuals frozen) for tractability —
    a disclosed approximation. Idealized cubane, def2-SVP, single-state NEVPT2.
  - SCF loaded from the Stage-17 chkfile (cached). Every number computed at run time;
    convergence flags and any failure recorded honestly.

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
NVIRT = 60          # frontier virtuals kept active for orbital rotation (rest frozen)


def fe3d_character(mf, mo_active):
    """Mean Löwdin Fe-3d weight of the given active MO columns."""
    fe3d = mf.mol.search_ao_label("Fe 3d")
    sqrtS = s17._sqrtm_sym(mf.get_ovlp())
    Cl = sqrtS @ mo_active
    return float(np.mean((Cl ** 2)[fe3d, :].sum(axis=0)))


def main():
    t_all = time.time()
    res = {"status": "running",
           "meta": {"system": "[Fe4S4(SH)4]2- cubane, frontier CAS(6e,6o)",
                    "basis": s17.BASIS,
                    "method": f"windowed converged CASSCF (core + deep virtuals frozen, "
                              f"{NVIRT} frontier virtuals active) + SC-NEVPT2",
                    "note": "SAME frontier subspace as Stage 17; a demonstration "
                            "subspace, not the full [4Fe-4S] magnetic space (40 qubits)."}}
    try:
        atoms = s17.build_cubane()
        mol = gto.M(atom=[(e, tuple(p)) for e, p in atoms], basis=s17.BASIS,
                    charge=s17.CHARGE, spin=s17.SPIN, verbose=0, max_memory=11000)
        mf = s17.converge_rohf(mol)
        scf_ok = bool(mf.converged)
        res["scf"] = {"e_hs_rohf": float(mf.e_tot), "converged": scf_ok}
        print(f"HS ROHF: E={mf.e_tot:.5f} conv={scf_ok}", flush=True)

        active, ncas, nelec, fe3d0 = s17.frontier_magnetic_cas(mf, NCAS)
        na = (nelec + (nelec % 2)) // 2
        nb = nelec - na
        ncore = (int(round(sum(mf.mo_occ))) - nelec) // 2
        nao = mf.mo_coeff.shape[1]
        # window: freeze [0, ncore) core and [ncore+ncas+NVIRT, nao) deep virtuals
        top = min(ncore + ncas + NVIRT, nao)
        frozen = list(range(ncore)) + list(range(top, nao))
        res["frontier"] = {"ncas": ncas, "nelec": nelec, "nelecas": [int(na), int(nb)],
                           "qubits_jw": 2 * ncas, "active_mo_idx": active,
                           "fe3d_char_initial": round(fe3d0, 3),
                           "ncore_frozen": int(ncore), "nvirt_active": int(top - ncore - ncas),
                           "based_on_converged_scf": scf_ok}
        print(f"frontier CAS({nelec}e,{ncas}o)=12q; init Fe-3d {fe3d0:.2f}; "
              f"windowed opt: core {ncore} frozen, {top-ncore-ncas} virtuals active",
              flush=True)

        mc = mcscf.CASSCF(mf, ncas, (na, nb)).density_fit()
        mc.frozen = frozen
        mc.max_cycle_macro = 200
        mc.conv_tol = 1e-8
        mc.verbose = 0
        mo = mc.sort_mo(active, base=0)
        t0 = time.time()
        e_cas = mc.kernel(mo)[0]
        occ = np.clip(np.linalg.eigvalsh(mc.fcisolver.make_rdm1(mc.ci, ncas, (na, nb)))[::-1],
                      0.0, 2.0)
        n_u = float(np.sum(occ * (2.0 - occ)))
        ss = float(mc.fcisolver.spin_square(mc.ci, ncas, (na, nb))[0])
        fe3d_final = fe3d_character(mf, mc.mo_coeff[:, ncore:ncore + ncas])
        casscf = {"e_casscf": float(e_cas), "casscf_converged": bool(mc.converged),
                  "cas_seconds": round(time.time() - t0, 1), "n_u": round(n_u, 3),
                  "noon": [round(float(x), 3) for x in occ], "s_squared": round(ss, 3),
                  "fe3d_char_final": round(fe3d_final, 3),
                  "n_u_stage17_unconverged": 3.985}
        print(f"CASSCF: E={e_cas:.6f} conv={mc.converged} n_u={n_u:.3f} "
              f"<S^2>={ss:.2f} Fe-3d(final)={fe3d_final:.2f} "
              f"({casscf['cas_seconds']}s)", flush=True)

        try:
            t1 = time.time()
            ec = mrpt.NEVPT(mc).kernel()
            casscf["nevpt2_corr"] = float(ec)
            casscf["e_nevpt2"] = float(e_cas + ec)
            casscf["nevpt2_seconds"] = round(time.time() - t1, 1)
            print(f"       + NEVPT2 corr={ec:.6f} -> E={e_cas+ec:.6f} "
                  f"({casscf['nevpt2_seconds']}s)", flush=True)
        except Exception as e:
            casscf["nevpt2_error"] = f"{type(e).__name__}: {e}"
            print("       NEVPT2 failed:", casscf["nevpt2_error"], flush=True)

        res["frontier_casscf_nevpt2"] = casscf
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
