#!/usr/bin/env python3
"""
Track-B, Fe-zeolite α-O — STEP 1: the bare ferryl core [FeO]2+  (Fe(IV)=O).

Known answer for the zeolite α-O site (Snyder et al., Nature 2016,
10.1038/nature19059): a high-spin S=2 Fe(IV)=O. IMPORTANT HONESTY NOTE: that
S=2 ground state is *enforced by the weak-field, constrained zeolite framework*.
This script computes only the BARE ferryl core [FeO]2+ in the gas phase — it
tests the engine on Fe(IV)=O spin-state + multireference physics, but it is NOT
the framework-embedded α-O. Whether the bare core is S=2 may differ from the
embedded site; the faithful α-O validation needs a zeolite cluster model (bigger
active space -> AWS), which is the next step. We report what the minimal model
gives, with that caveat explicit.

Same consistent-active-space method as the FeO+ anchor: PySCF ROHF (highest-spin
reference defines the space) -> AVAS(Fe 3d + O 2p) -> state-specific CASSCF in the
SAME (ncas, nelecas) for every spin, only the alpha/beta split changes. Minimal
CAS / def2-SVP, CASSCF only (no dynamic correlation): we test the QUALITATIVE
spin-state ordering + NOON multireference, not a quantitative gap.

Every number is computed at run time; nothing is invented.
"""
import json
import numpy as np
from pyscf import gto, scf, mcscf
from pyscf.mcscf import avas

BASIS = "def2-svp"
AO_LABELS = ["Fe 3d", "O 2p"]
RS = [1.55, 1.62, 1.69, 1.76, 1.83]            # Angstrom, CASSCF scan
SPINS = {"quintet": 4, "triplet": 2, "singlet": 0}   # value = 2S ; known α-O = S=2 (quintet)
CHARGE = 2                                       # [FeO]2+, Fe(IV)=O


def _scf_ref(r):
    """Highest-spin (quintet) ROHF reference that defines the active space."""
    mol = gto.M(atom=f"Fe 0 0 0; O 0 0 {r}", charge=CHARGE, spin=4, basis=BASIS, verbose=0)
    mf = scf.ROHF(mol); mf.level_shift = 0.3; mf.max_cycle = 300
    mf.kernel()
    if not mf.converged:
        mf = scf.ROHF(mol).newton(); mf.kernel()
    return mf


def _cas(mf, ncas, nelec_total, spin2, mo):
    na = (nelec_total + spin2) // 2
    nb = nelec_total - na
    if nb < 0 or na > ncas or nb > ncas:
        return {"error": f"bad split na={na} nb={nb} in {ncas}o"}
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
    n_u = float(np.sum(occ * (2.0 - occ)))
    return {"e_cas": float(e), "cas_conv": bool(mc.converged), "na": int(na), "nb": int(nb),
            "n_u": round(n_u, 3), "excess": round(n_u - spin2, 3),
            "noon": [round(float(x), 3) for x in occ]}


def run_geometry(r):
    out = {"r": r}
    try:
        mf = _scf_ref(r)
        out["scf_conv"] = bool(mf.converged)
        ncas, nelec, mo = avas.avas(mf, AO_LABELS, canonicalize=False)
        ncas = int(ncas); nelec = int(nelec)
        out.update(ncas=ncas, nelec=nelec, qubits_jw=2 * ncas)
        for name, s2 in SPINS.items():
            out[name] = _cas(mf, ncas, nelec, s2, mo)
    except Exception as e:
        out["error"] = f"{type(e).__name__}: {e}"
    return out


def main():
    results = {"meta": {"basis": BASIS, "ao_labels": AO_LABELS, "scan_A": RS,
                        "system": "[FeO]2+ bare ferryl core (Fe(IV)=O)",
                        "model_caveat": "BARE core, gas phase — NOT the framework-embedded zeolite α-O; "
                                        "α-O S=2 is enforced by the zeolite weak field (needs cluster model, AWS).",
                        "zeolite_known_answer": "α-O ground state = S=2 (high-spin Fe(IV)=O), Nature 2016"},
               "scan": [], "minima": {}}
    for r in RS:
        g = run_geometry(r); results["scan"].append(g)
        line = f"r={r}"
        if "ncas" in g:
            line += f" CAS({g['nelec']}e,{g['ncas']}o)={g['qubits_jw']}q"
        for name in SPINS:
            d = g.get(name, {})
            line += f" | {name[:4]}={d.get('e_cas','ERR')}"
        print(line)

    for name in SPINS:
        good = [(g["r"], g[name], g) for g in results["scan"] if name in g and "e_cas" in g[name]]
        if good:
            r, d, g = min(good, key=lambda x: x[1]["e_cas"])
            results["minima"][name] = {"r": r, **d, "ncas": g["ncas"], "nelec": g["nelec"]}

    if results["minima"]:
        gs = min(results["minima"], key=lambda k: results["minima"][k]["e_cas"])
        e_gs = results["minima"][gs]["e_cas"]
        order = sorted(((round((results["minima"][k]["e_cas"] - e_gs) * 627.509, 2), k)
                        for k in results["minima"]))
        results["verdict"] = {
            "ground_state_our_casscf": gs,
            "zeolite_alpha_O_known": "S=2 (quintet)",
            "bare_core_matches_S2": bool(gs == "quintet"),
            "relative_kcal_above_ground": {k: v for v, k in order},
            "note": "bare ferryl core != embedded α-O; this checks the engine on Fe(IV)=O, "
                    "not the framework-enforced S=2.",
        }
        print("\nVERDICT:", json.dumps(results["verdict"], ensure_ascii=False))

    with open("calc/fe_ferryl_results.json", "w") as f:
        json.dump(results, f, indent=2)
    print("wrote calc/fe_ferryl_results.json")


if __name__ == "__main__":
    main()
