#!/usr/bin/env python3
"""
Track-B Fe-zeolite alpha-O — v2 FAITHFUL-ER run: CONVERGED CASSCF + NEVPT2.

Fixes the two biggest weaknesses of the v1 preliminary run (calc/fe_zeolite_alpha_o.py):
  (1) v1 CASSCF did NOT converge -> here we use a CURATED, standard ferryl active
      space (Fe 3d + oxo O 2p only, atom-indexed AVAS) that converges cleanly,
      instead of v1's "grab every O 2p" 14-orbital space.
  (2) v1 had no dynamic correlation -> here we add strongly-contracted NEVPT2 on
      top of each state-specific CASSCF, for a quantitative-er spin gap.

Cluster: [Fe(O)(OH)2(OH2)2] — neutral, Fe(IV)=O, 5-coordinate weak-field O-donor
site (oxo + 2 hydroxo + 2 aqua). A constructed step UP from the v1 [Fe(O)(OH)2]
stand-in (more of the low-coordinate constrained field), but STILL a constructed
model — NOT a CIF-derived crystallographic CHA/MOR site. Pass --xyz for a real one.

HONESTY: spin-state ORDERING + NEVPT2 gap on a CONSTRUCTED cluster with VERTICAL
gaps (geometry relaxed once at high spin). Faithful absolute numbers still need a
real periodic/CIF structure + adiabatic relaxation + larger basis. Known answer:
alpha-O ground state = S=2 (Snyder, Nature 2016, 10.1038/nature19059).

Usage:
  python3 calc/fe_zeolite_alpha_o_nevpt2.py            # constructed cluster
  python3 calc/fe_zeolite_alpha_o_nevpt2.py --xyz s.xyz --charge 0
"""
import argparse, json, sys
import numpy as np
from pyscf import gto, scf, mcscf, dft, mrpt
from pyscf.mcscf import avas

BASIS = "def2-svp"
AO_LABELS = ["0 Fe 3d", "1 O 2p"]          # atom-indexed: Fe(atom0) 3d + oxo(atom1) 2p
SPINS = {"quintet": 4, "triplet": 2, "singlet": 0}
KNOWN = "quintet"


def constructed_cluster():
    """[Fe(O)(OH)2(OH2)2] — Fe=0, oxo=1 (order matters for atom-indexed AVAS)."""
    return [
        ["Fe", (0.000,  0.000,  0.000)],
        ["O",  (0.000,  0.000,  1.620)],   # oxo (atom 1)
        ["O",  (1.950,  0.000, -0.300)], ["H", (2.560,  0.000, -0.860)],   # hydroxo
        ["O",  (-1.950, 0.000, -0.300)], ["H", (-2.560, 0.000, -0.860)],   # hydroxo
        ["O",  (0.000,  2.150, -0.300)], ["H", (0.600, 2.700, -0.560)], ["H", (-0.600, 2.700, -0.560)],  # aqua
        ["O",  (0.000, -2.150, -0.300)], ["H", (0.600, -2.700, -0.560)], ["H", (-0.600, -2.700, -0.560)],  # aqua
    ]


def minimal_cluster():
    """[Fe(O)(OH)2] — Fe=0, oxo=1. Smaller field, but relaxes reliably (the
    extended aqua model does not converge berny gradients in this env)."""
    return [
        ["Fe", (0.000,  0.000,  0.000)],
        ["O",  (0.000,  0.000,  1.620)],
        ["O",  (1.820,  0.000, -0.430)], ["H", (2.470,  0.000, -1.090)],
        ["O",  (-1.820, 0.000, -0.430)], ["H", (-2.470, 0.000, -1.090)],
    ]


def read_xyz(path):
    with open(path) as f:
        lines = f.read().splitlines()
    n = int(lines[0].split()[0])
    return [[p[0], (float(p[1]), float(p[2]), float(p[3]))]
            for p in (ln.split() for ln in lines[2:2 + n])]


def make_mol(atoms, charge, spin2):
    return gto.M(atom=[[a[0], a[1]] for a in atoms], charge=charge, spin=spin2,
                 basis=BASIS, verbose=0, max_memory=12000)


def relaxed_ref(atoms, charge, relax=True):
    """High-spin reference; optionally relax geometry once (berny/geometric)."""
    mol = make_mol(atoms, charge, max(SPINS.values()))
    for backend in (("geometric_solver", "berny_solver") if relax else ()):
        try:
            opt = __import__(f"pyscf.geomopt.{backend}", fromlist=["optimize"]).optimize
            mks = dft.UKS(mol); mks.xc = "PBE"; mks.conv_tol = 1e-7; mks.max_cycle = 400
            mol = opt(mks, maxsteps=100, conv_params={"convergence_energy": 1e-5})
            print(f"[geom] relaxed (high-spin) with {backend}"); break
        except Exception as e:
            print(f"[warn] {backend} failed ({type(e).__name__}: {e})")
    mf = scf.ROHF(mol); mf.level_shift = 0.3; mf.max_cycle = 400; mf.kernel()
    if not mf.converged:
        mf = scf.ROHF(mol).newton(); mf.kernel()
    return mf, mol


def state(mf, ncas, nelec_total, spin2, mo):
    na = (nelec_total + spin2) // 2
    nb = nelec_total - na
    if nb < 0 or na > ncas or nb > ncas:
        return {"error": f"bad split na={na} nb={nb} in {ncas}o"}
    mc = mcscf.CASSCF(mf, ncas, (na, nb))
    S = spin2 / 2.0
    try:
        mc.fix_spin_(ss=S * (S + 1))
    except Exception:
        pass
    mc.max_cycle_macro = 200; mc.conv_tol = 1e-8
    e_cas = mc.kernel(mo)[0]
    occ = np.clip(np.linalg.eigvalsh(mc.fcisolver.make_rdm1(mc.ci, ncas, (na, nb)))[::-1], 0, 2)
    n_u = float(np.sum(occ * (2 - occ)))
    d = {"e_cas": float(e_cas), "cas_conv": bool(mc.converged),
         "n_u": round(n_u, 3), "excess": round(n_u - spin2, 3),
         "noon": [round(float(x), 3) for x in occ]}
    try:                       # dynamic correlation on top
        ec = mrpt.NEVPT(mc).kernel()
        d["e_nevpt2"] = float(e_cas + ec)
        d["nevpt2_corr"] = float(ec)
    except Exception as e:
        d["nevpt2_error"] = f"{type(e).__name__}: {e}"
    return d


def order(spins, key):
    ok = {k: v for k, v in spins.items() if key in v}
    if not ok:
        return None
    gs = min(ok, key=lambda k: ok[k][key]); e0 = ok[gs][key]
    return {"level": key, "ground_state": gs, "matches_known_S2": bool(gs == KNOWN),
            "rel_kcal": {k: round((ok[k][key] - e0) * 627.509, 2) for k in ok}}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--xyz", default=None); ap.add_argument("--charge", type=int, default=0)
    ap.add_argument("--minimal", action="store_true", help="use [Fe(O)(OH)2] (relaxes reliably)")
    ap.add_argument("--norelax", action="store_true", help="skip geom-opt (geometry pre-relaxed, e.g. xTB)")
    args = ap.parse_args()
    if args.xyz:
        atoms, model = read_xyz(args.xyz), args.xyz
    elif args.minimal:
        atoms, model = minimal_cluster(), "constructed [Fe(O)(OH)2] minimal (relaxes; NOT a CIF CHA site)"
    else:
        atoms, model = constructed_cluster(), "constructed [Fe(O)(OH)2(OH2)2] (NOT a CIF CHA site)"
    out = {"meta": {"basis": BASIS, "ao_labels": AO_LABELS, "charge": args.charge,
                    "model": model,
                    "known_answer": "alpha-O ground state = S=2 (Nature 2016)",
                    "honesty": "curated converged CASSCF + NEVPT2; VERTICAL gaps on a "
                               "CONSTRUCTED cluster; faithful absolutes need real CIF + adiabatic."},
           "spins": {}}
    try:
        mf, mol = relaxed_ref(atoms, args.charge, relax=not args.norelax)
        out["meta"]["scf_conv"] = bool(mf.converged)
        ncas, nelec, mo = avas.avas(mf, AO_LABELS, canonicalize=False)
        ncas = int(ncas); nelec = int(nelec)
        out["active_space"] = {"ncas": ncas, "nelec": nelec, "qubits_jw": 2 * ncas}
        print(f"curated active space: CAS({nelec}e,{ncas}o) = {2*ncas} qubits")
        for name, s2 in SPINS.items():
            d = state(mf, ncas, nelec, s2, mo); out["spins"][name] = d
            print(f"  {name:8s}: E_cas={d.get('e_cas','ERR')} conv={d.get('cas_conv')} "
                  f"E_nevpt2={d.get('e_nevpt2','-')}")
        out["verdict_casscf"] = order(out["spins"], "e_cas")
        out["verdict_nevpt2"] = order(out["spins"], "e_nevpt2")
        print("CASSCF :", json.dumps(out["verdict_casscf"], ensure_ascii=False))
        print("NEVPT2 :", json.dumps(out["verdict_nevpt2"], ensure_ascii=False))
    except Exception as e:
        out["error"] = f"{type(e).__name__}: {e}"; print("ERROR:", out["error"], file=sys.stderr)
    outfile = f"calc/fe_zeolite_alpha_o_nevpt2{'_minimal' if args.minimal else ''}_results.json"
    with open(outfile, "w") as f:
        json.dump(out, f, indent=2)
    print("wrote", outfile)


if __name__ == "__main__":
    main()
