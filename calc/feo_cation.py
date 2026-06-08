#!/usr/bin/env python3
"""
Track-B validation anchor — FeO+ (gas-phase paradigm of TM-oxo C-H activation).

Known answer to reproduce: the ground electronic state of FeO+ is the high-spin
SEXTET (6-Sigma+); the quartet is a low-lying excited state whose crossing with
the sextet is the textbook "two-state reactivity" of FeO+ + CH4 -> Fe+ + CH3OH
(Schroder/Shaik/Schwarz; Yoshizawa JACS 1998 10.1021/ja971723u, 2000
10.1021/ja0017965).

Same pipeline as the rest of the diary: PySCF SCF -> AVAS(Fe 3d + O 2p) -> CASSCF,
base-invariant NOON metric for multireference character, JW qubit count = 2*ncas.
Honest: minimal CAS / def2-SVP, fixed adiabatic scan; quantitative gaps need
larger basis + dynamic correlation (CASPT2/NEVPT2). We only test the QUALITATIVE
known answer (which spin state is the ground state) + the multireference signal.

Nothing here is invented: every number printed is computed at run time.
"""
import json
import numpy as np
from pyscf import gto, scf, mcscf
from pyscf.mcscf import avas

BASIS = "def2-svp"
AO_LABELS = ["Fe 3d", "O 2p"]
RS = [1.55, 1.60, 1.65, 1.70, 1.75, 1.80]      # Angstrom, CASSCF scan
SPINS = {"sextet": 5, "quartet": 3}             # value = 2S = # unpaired e


def run_point(r, spin2):
    """One CASSCF point. Returns energies + active-space info, or error dict."""
    out = {"r": r, "spin2": spin2}
    try:
        mol = gto.M(atom=f"Fe 0 0 0; O 0 0 {r}", charge=1, spin=spin2,
                    basis=BASIS, verbose=0)
        mf = scf.ROHF(mol)
        mf.level_shift = 0.3
        mf.max_cycle = 300
        e_hf = mf.kernel()
        if not mf.converged:                    # second try: Newton (SOSCF)
            mf = scf.ROHF(mol).newton()
            e_hf = mf.kernel()
        out["scf_conv"] = bool(mf.converged)
        out["e_hf"] = float(e_hf)

        ncas, nelecas, mo = avas.avas(mf, AO_LABELS, canonicalize=False)
        mc = mcscf.CASSCF(mf, ncas, nelecas)
        mc.max_cycle_macro = 200
        e_cas = mc.kernel(mo)[0]
        out["cas_conv"] = bool(mc.converged)
        out["ncas"] = int(ncas)
        out["nelecas_total"] = int(nelecas if np.isscalar(nelecas) else sum(nelecas))
        out["qubits_jw"] = 2 * int(ncas)
        out["e_cas"] = float(e_cas)

        # base-invariant multireference metric from active natural occupations
        casdm1 = mc.fcisolver.make_rdm1(mc.ci, mc.ncas, mc.nelecas)
        occ = np.linalg.eigvalsh(casdm1)[::-1]
        occ = np.clip(occ, 0.0, 2.0)
        n_u = float(np.sum(occ * (2.0 - occ)))   # effectively-unpaired electrons
        out["noon"] = [round(float(x), 3) for x in occ]
        out["n_u"] = round(n_u, 3)
        out["excess"] = round(n_u - spin2, 3)    # beyond formal spin (2S = spin2)
    except Exception as e:
        out["error"] = f"{type(e).__name__}: {e}"
    return out


def main():
    results = {"meta": {"basis": BASIS, "ao_labels": AO_LABELS,
                        "scan_A": RS, "system": "FeO+ (charge +1)",
                        "known_answer": "ground state = sextet (6-Sigma+)"},
               "scan": {}, "minima": {}}
    for name, s2 in SPINS.items():
        pts = [run_point(r, s2) for r in RS]
        results["scan"][name] = pts
        good = [p for p in pts if "e_cas" in p]
        if good:
            m = min(good, key=lambda p: p["e_cas"])
            results["minima"][name] = m
            print(f"{name:8s}: min E(CASSCF) = {m['e_cas']:.6f} Ha at r={m['r']} A "
                  f"| CAS({m['nelecas_total']}e,{m['ncas']}o)={m['qubits_jw']}q "
                  f"| excess={m['excess']}")
        else:
            print(f"{name:8s}: no converged CASSCF points")

    if {"sextet", "quartet"} <= set(results["minima"]):
        e6 = results["minima"]["sextet"]["e_cas"]
        e4 = results["minima"]["quartet"]["e_cas"]
        gap_kcal = (e4 - e6) * 627.509          # quartet - sextet
        ground = "sextet" if e6 < e4 else "quartet"
        results["verdict"] = {
            "ground_state": ground,
            "sextet_minus_quartet_kcal": round((e6 - e4) * 627.509, 2),
            "quartet_above_sextet_kcal": round(gap_kcal, 2),
            "matches_known_sextet": bool(ground == "sextet"),
        }
        print("\nVERDICT")
        print(f"  ground state (our CASSCF): {ground}  "
              f"(known answer: sextet 6-Sigma+)")
        print(f"  adiabatic gap quartet-above-sextet = {gap_kcal:.2f} kcal/mol")
        print(f"  reproduces known sextet ground state: "
              f"{results['verdict']['matches_known_sextet']}")

    with open("calc/feo_cation_results.json", "w") as f:
        json.dump(results, f, indent=2)
    print("\nwrote calc/feo_cation_results.json")


if __name__ == "__main__":
    main()
