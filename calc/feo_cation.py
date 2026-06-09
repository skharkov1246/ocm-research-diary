#!/usr/bin/env python3
"""
Track-B validation anchor — FeO+ (gas-phase paradigm of TM-oxo C-H activation).

Known answer to reproduce: the ground electronic state of FeO+ is the high-spin
SEXTET (6-Sigma+); the quartet is a low-lying excited state whose crossing with
the sextet is the textbook "two-state reactivity" of FeO+ + CH4 -> Fe+ + CH3OH
(Schroder/Shaik/Schwarz; Yoshizawa JACS 1998 10.1021/ja971723u, 2000
10.1021/ja0017965).

METHOD FIX (v2): a first pass let AVAS choose the active space PER spin state and
it returned DIFFERENT active electron counts (sextet 11e/8o vs quartet 9e/8o),
so that energy comparison was not apples-to-apples. Here the active space is
defined ONCE per geometry from the sextet ROHF reference (AVAS -> Fe 3d + O 2p),
and BOTH spin states are solved as state-specific CASSCF in that SAME
(ncore, ncas, nelecas) space — only the alpha/beta split (the target spin)
changes. That is a consistent adiabatic comparison.

Honest scope: minimal CAS(Fe 3d + O 2p) / def2-SVP, CASSCF only (no dynamic
correlation). We test only the QUALITATIVE known answer (which spin state is the
ground state) + the base-invariant NOON multireference signal + JW qubit count.
A quantitative gap would need CASPT2/NEVPT2 + a larger basis. Every number is
computed at run time; nothing is invented.
"""
import json
import numpy as np
from pyscf import gto, scf, mcscf
from pyscf.mcscf import avas

BASIS = "def2-svp"
AO_LABELS = ["Fe 3d", "O 2p"]
RS = [1.55, 1.62, 1.69, 1.76, 1.83]      # Angstrom, CASSCF scan
SPINS = {"sextet": 5, "quartet": 3}       # value = 2S = # unpaired e (known GS = sextet)


def _scf_sextet(r):
    """Sextet ROHF reference that defines the active space at geometry r."""
    mol = gto.M(atom=f"Fe 0 0 0; O 0 0 {r}", charge=1, spin=5, basis=BASIS, verbose=0)
    mf = scf.ROHF(mol)
    mf.level_shift = 0.3
    mf.max_cycle = 300
    mf.kernel()
    if not mf.converged:
        mf = scf.ROHF(mol).newton()
        mf.kernel()
    return mf


def _cas(mf, ncas, nelec_total, spin2, mo):
    """State-specific CASSCF for the requested spin in a FIXED active space."""
    na = (nelec_total + spin2) // 2
    nb = nelec_total - na
    if nb < 0 or na > ncas or nb > ncas:
        return {"error": f"bad alpha/beta split na={na} nb={nb} in {ncas}o"}
    mc = mcscf.CASSCF(mf, ncas, (na, nb))
    try:
        S = spin2 / 2.0
        mc.fix_spin_(ss=S * (S + 1))
    except Exception:
        pass
    mc.max_cycle_macro = 200
    e = mc.kernel(mo)[0]
    casdm1 = mc.fcisolver.make_rdm1(mc.ci, ncas, (na, nb))
    occ = np.clip(np.linalg.eigvalsh(casdm1)[::-1], 0.0, 2.0)
    n_u = float(np.sum(occ * (2.0 - occ)))       # effectively-unpaired electrons
    return {"e_cas": float(e), "cas_conv": bool(mc.converged), "na": int(na), "nb": int(nb),
            "n_u": round(n_u, 3), "excess": round(n_u - spin2, 3),
            "noon": [round(float(x), 3) for x in occ]}


def run_geometry(r):
    out = {"r": r}
    try:
        mf = _scf_sextet(r)                       # ONE reference defines the space
        out["scf_conv"] = bool(mf.converged)
        ncas, nelec, mo = avas.avas(mf, AO_LABELS, canonicalize=False)
        ncas = int(ncas); nelec = int(nelec)
        out.update(ncas=ncas, nelec=nelec, qubits_jw=2 * ncas)
        for name, s2 in SPINS.items():            # SAME (ncas, nelec) for both spins
            out[name] = _cas(mf, ncas, nelec, s2, mo)
    except Exception as e:
        out["error"] = f"{type(e).__name__}: {e}"
    return out


def main():
    results = {"meta": {"basis": BASIS, "ao_labels": AO_LABELS, "scan_A": RS,
                        "system": "FeO+ (charge +1)", "method": "state-specific CASSCF, fixed CAS",
                        "known_answer": "ground state = sextet (6-Sigma+)"},
               "scan": [], "minima": {}}
    for r in RS:
        g = run_geometry(r)
        results["scan"].append(g)
        tag = f"r={r}"
        if "ncas" in g:
            tag += f" CAS({g['nelec']}e,{g['ncas']}o)={g['qubits_jw']}q"
        for name in SPINS:
            d = g.get(name, {})
            tag += f" | {name}={d.get('e_cas','ERR')}"
        print(tag)

    for name in SPINS:
        good = [(g["r"], g[name]) for g in results["scan"]
                if name in g and "e_cas" in g[name]]
        if good:
            r, d = min(good, key=lambda x: x[1]["e_cas"])
            results["minima"][name] = {"r": r, **d,
                                       "ncas": next(g["ncas"] for g in results["scan"] if g["r"] == r),
                                       "nelec": next(g["nelec"] for g in results["scan"] if g["r"] == r)}

    if {"sextet", "quartet"} <= set(results["minima"]):
        e6 = results["minima"]["sextet"]["e_cas"]
        e4 = results["minima"]["quartet"]["e_cas"]
        ground = "sextet" if e6 < e4 else "quartet"
        results["verdict"] = {
            "ground_state_our_casscf": ground,
            "known_answer": "sextet",
            "matches_known_sextet": bool(ground == "sextet"),
            "quartet_above_sextet_kcal": round((e4 - e6) * 627.509, 2),
            "active_space": f"CAS({results['minima']['sextet']['nelec']}e,"
                            f"{results['minima']['sextet']['ncas']}o) = "
                            f"{2*results['minima']['sextet']['ncas']} qubits (JW)",
        }
        print("\nVERDICT:", json.dumps(results["verdict"], ensure_ascii=False))

    with open("calc/feo_cation_results.json", "w") as f:
        json.dump(results, f, indent=2)
    print("wrote calc/feo_cation_results.json")


if __name__ == "__main__":
    main()
