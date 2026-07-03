#!/usr/bin/env python3
"""FAST real-CHA spin ordering (pyscf-only; no ase/zse needed — reads the
xTB-relaxed geometry from disk). Robust ROKS-PBE SCF (minutes, not the 1h ROHF)
+ AVAS(Fe 3d + oxo 2p) + state-specific CASCI for S=2/1/0 on the SAME orbitals.
Coarser than CASSCF+NEVPT2 (no orbital opt / no dynamic correlation) but delivers
the spin ORDERING on the real IZA-CHA relaxed geometry before the ephemeral
container cycles. Atom order in xyz: 0=Fe, 1=oxo."""
import json, sys, time
import numpy as np
from pyscf import gto, scf, mcscf, dft
from pyscf.mcscf import avas

BASIS = sys.argv[1] if len(sys.argv) > 1 else "def2-svp"


def read_xyz(path):
    with open(path) as f:
        lines = f.read().splitlines()
    n = int(lines[0].split()[0])
    return [[p[0], (float(p[1]), float(p[2]), float(p[3]))]
            for p in (ln.split() for ln in lines[2:2 + n])]


atoms = read_xyz("calc/alpha_o_cha_relaxed.xyz")
mol = gto.M(atom=[[s, xyz] for s, xyz in atoms], charge=0, spin=4,
            basis=BASIS, verbose=0, max_memory=12000)
print("basis:", BASIS, "| n atoms:", len(atoms), "| n bf:", mol.nao, flush=True)

mf = dft.ROKS(mol).density_fit(); mf.xc = "PBE"; mf.level_shift = 0.2
mf.conv_tol = 1e-6; mf.max_cycle = 300
t = time.time(); mf.kernel()
print("ROKS-PBE conv=%s E=%.4f (%.0fs)" % (mf.converged, mf.e_tot, time.time() - t), flush=True)

ncas, nelec, mo = avas.avas(mf, ["0 Fe 3d", "1 O 2p"], canonicalize=True)
ncas = int(ncas); nelec = int(nelec)
print("CAS(%de,%do) = %d qubits" % (nelec, ncas, 2 * ncas), flush=True)

out = {"meta": {"model": "REAL IZA-CHA 6-ring cluster (xTB-relaxed active site)",
                "method": "ROKS-PBE orbitals + state-specific CASCI (FAST; no orbital opt / no NEVPT2)",
                "basis": BASIS, "active_space": f"CAS({nelec}e,{ncas}o)={2*ncas}q",
                "known_answer": "alpha-O ground state = S=2",
                "honesty": "coarse (CASCI on DFT orbitals) — gives ORDERING, not a quantitative gap"},
       "spins": {}}
res = {}
for name, s2 in [("quintet", 4), ("triplet", 2), ("singlet", 0)]:
    na = (nelec + s2) // 2; nb = nelec - na
    mc = mcscf.CASCI(mf, ncas, (na, nb))
    try: mc.fix_spin_(ss=(s2 / 2) * (s2 / 2 + 1))
    except Exception: pass
    e = mc.kernel(mo)[0]; res[name] = float(e)
    out["spins"][name] = {"e_casci": float(e)}
    print("  %s: E=%.5f" % (name, e), flush=True)

gs = min(res, key=res.get); e0 = res[gs]
out["verdict"] = {"ground_state": gs, "matches_known_S2": bool(gs == "quintet"),
                  "rel_kcal": {k: round((res[k] - e0) * 627.509, 2) for k in res}}
print("VERDICT:", json.dumps(out["verdict"]), flush=True)
with open("calc/fe_zeolite_cha_fast_results.json", "w") as f:
    json.dump(out, f, indent=2)
print("wrote results", flush=True)
