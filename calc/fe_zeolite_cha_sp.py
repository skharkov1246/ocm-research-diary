#!/usr/bin/env python3
"""Spin-state CASSCF(+NEVPT2) single-point on the RELAXED real-CHA alpha-O cluster
(calc/alpha_o_cha_relaxed.xyz), using DENSITY FITTING to make the 38-atom job
tractable. Geometry already xTB-relaxed -> no geom-opt here. Curated active space
(Fe 3d + oxo 2p, atom-indexed AVAS). Known answer: alpha-O ground state = S=2.

This is the faithful test: real IZA CHA framework geometry, standard 2-Al 6-ring
[FeO]2+ model. Honesty: cluster (not periodic), def2-SVP, 2-Al siting = para,
xTB active-site geometry. Validates spin ORDERING on a real-framework geometry."""
import json, os, sys, time
import numpy as np
from pyscf import gto, scf, mcscf, mrpt
from pyscf.mcscf import avas
from ase.io import read

SPINS = {"quintet": 4, "triplet": 2, "singlet": 0}
a = read("calc/alpha_o_cha_relaxed.xyz")
SYM = a.get_chemical_symbols(); POS = a.get_positions()


def mol(spin2):
    return gto.M(atom=[[s, tuple(p)] for s, p in zip(SYM, POS)],
                 charge=0, spin=spin2, basis="def2-svp", verbose=0, max_memory=12000)


def run():
    out = {"meta": {"model": "REAL IZA CHA 6-ring, 2 Al (para), Fe(IV)=O, silanol caps; "
                              "xTB-relaxed active site; def2-SVP + density fitting",
                    "known_answer": "alpha-O ground state = S=2 (Nature 2016)",
                    "honesty": "cluster (not periodic), def2-SVP, para Al-siting, xTB geometry; "
                               "validates spin ORDERING on real-framework geometry."},
           "spins": {}}
    t = time.time()
    m = mol(4)
    out["meta"]["nbf"] = int(m.nao)
    mf = scf.ROHF(m).density_fit(); mf.chkfile = "calc/cha_scf.chk"
    if os.path.exists(mf.chkfile):
        mf.init_guess = "chkfile"                    # reuse a prior SCF to skip the slow part
    mf.level_shift = 0.3; mf.max_cycle = 300
    mf.kernel()
    if not mf.converged:
        mf = mf.newton(); mf.kernel()
    out["meta"]["scf_conv"] = bool(mf.converged)
    print(f"DF-SCF conv={mf.converged} E={mf.e_tot:.4f} ({time.time()-t:.0f}s, {m.nao} bf)")
    ncas, nelec, mo = avas.avas(mf, ["0 Fe 3d", "1 O 2p"], canonicalize=False)
    ncas, nelec = int(ncas), int(nelec)
    out["active_space"] = {"ncas": ncas, "nelec": nelec, "qubits_jw": 2 * ncas}
    print(f"curated CAS({nelec}e,{ncas}o) = {2*ncas} qubits")
    for name, s2 in SPINS.items():
        na = (nelec + s2) // 2; nb = nelec - na
        mc = mcscf.CASSCF(mf, ncas, (na, nb))
        try: mc.fix_spin_(ss=(s2/2)*(s2/2+1))
        except Exception: pass
        mc.max_cycle_macro = 200
        e = mc.kernel(mo)[0]
        d = {"e_cas": float(e), "cas_conv": bool(mc.converged)}
        if os.environ.get("CHA_NEVPT2") == "1":      # NEVPT2 is the AWS-scale part; off by default
            try:
                d["e_nevpt2"] = float(e + mrpt.NEVPT(mc).kernel())
            except Exception as ex:
                d["nevpt2_err"] = f"{type(ex).__name__}: {str(ex)[:60]}"
        out["spins"][name] = d
        print(f"  {name:8s} cas_conv={d['cas_conv']} E_cas={d['e_cas']:.4f} "
              f"E_nevpt2={d.get('e_nevpt2','-')} ({time.time()-t:.0f}s)")

    for lvl in ("e_cas", "e_nevpt2"):
        ok = {k: v for k, v in out["spins"].items() if lvl in v}
        if ok:
            gs = min(ok, key=lambda k: ok[k][lvl]); e0 = ok[gs][lvl]
            out[f"verdict_{lvl}"] = {"ground_state": gs, "matches_known_S2": bool(gs == "quintet"),
                                     "rel_kcal": {k: round((ok[k][lvl]-e0)*627.509, 2) for k in ok}}
            print(lvl, "->", json.dumps(out[f"verdict_{lvl}"], ensure_ascii=False))
    return out


if __name__ == "__main__":
    try:
        out = run()
    except Exception as e:
        out = {"error": f"{type(e).__name__}: {e}"}; print("ERROR:", out["error"], file=sys.stderr)
    with open("calc/fe_zeolite_cha_results.json", "w") as f:
        json.dump(out, f, indent=2)
    print("wrote calc/fe_zeolite_cha_results.json")
