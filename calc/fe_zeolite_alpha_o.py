#!/usr/bin/env python3
"""
Track-B, Fe-zeolite alpha-O — EMBEDDED CLUSTER validation (AWS-ready).

Goal: test whether the MatterForge pipeline reproduces the KNOWN answer for the
zeolite alpha-O active site — a high-spin S=2 Fe(IV)=O ground state
(Snyder et al., Nature 2016, 10.1038/nature19059; PNAS 2018,
10.1073/pnas.1721717115; Chem. Soc. Rev. 2026, 10.1039/D5CS00496A).

Why this is the AWS run (not the bare core): step 1 showed the bare [FeO]2+ core
is unbound and does NOT carry the framework's weak field, so its ground state is
S=1, not S=2 (calc/fe_ferryl*.py). The S=2 of alpha-O is ENFORCED by the
low-coordinate, weak-field zeolite site, so the faithful test needs an embedded
cluster: an Fe(IV)=O with framework O donors (+ optional Madelung point charges,
same idea as the [NiO6] cathode model). That cluster's active space (Fe 3d + oxo
O 2p + framework O 2p) and full VQE statevector are heavy -> run on a big box.

Method (identical philosophy to the FeO+ anchor that already validated 6-Sigma+):
  geometry (DFT relax per spin, optional) -> AVAS(Fe 3d + O 2p) ->
  state-specific CASSCF in a CONSISTENT (ncas, nelecas) for S=2 / S=1 / S=0 ->
  CASCI reference + base-invariant NOON multireference -> JW qubit count ->
  UCCSD-VQE on a tractable sub-CAS (engine gate; full-CAS VQE is the heavy lift).

HONESTY:
  * Default geometry is an IDEALIZED minimal model [Fe(O)(OH)2] (neutral, Fe(IV),
    low-coordinate, weak-field O donors). It is a stand-in, not a crystallographic
    alpha-O site. For a faithful run pass a real cluster with  --xyz STRUCTURE.xyz
    (e.g. an Fe-exchanged CHA/MOR/FER fragment, H-terminated, optionally with
    point charges). The verdict is only as good as the geometry.
  * CASSCF has no dynamic correlation; a quantitative spin gap needs CASPT2/NEVPT2
    + a larger basis. We validate the spin-state ORDERING (is the ground state
    S=2?), plus the multireference signal. Every number is computed at run time.

Usage:
  python3 calc/fe_zeolite_alpha_o.py                 # smoke: idealized geom, SCF+CASSCF, no opt
  python3 calc/fe_zeolite_alpha_o.py --full          # AWS: DFT geom-opt per spin + CASSCF + VQE gate
  python3 calc/fe_zeolite_alpha_o.py --xyz site.xyz --charge 0 --full
"""
import argparse, json, sys
import numpy as np
from pyscf import gto, scf, mcscf, dft
from pyscf.mcscf import avas

BASIS = "def2-svp"
AO_LABELS = ["Fe 3d", "O 2p"]
SPINS = {"quintet": 4, "triplet": 2, "singlet": 0}   # 2S ; known alpha-O = S=2 (quintet)
KNOWN = "quintet"


def idealized_cluster():
    """Minimal neutral Fe(IV)=O weak-field stand-in: [Fe(=O)(OH)2]. IDEALIZED."""
    return [
        ["Fe", (0.000,  0.000,  0.000)],
        ["O",  (0.000,  0.000,  1.620)],   # ferryl oxo (Fe=O)
        ["O",  (1.820,  0.000, -0.430)],   # framework-O stand-in (hydroxo)
        ["H",  (2.470,  0.000, -1.090)],
        ["O",  (-1.820, 0.000, -0.430)],
        ["H",  (-2.470, 0.000, -1.090)],
    ]


def read_xyz(path):
    atoms = []
    with open(path) as f:
        lines = f.read().splitlines()
    n = int(lines[0].split()[0])
    for ln in lines[2:2 + n]:
        p = ln.split()
        atoms.append([p[0], (float(p[1]), float(p[2]), float(p[3]))])
    return atoms


def make_mol(atoms, charge, spin2):
    return gto.M(atom=[[a[0], a[1]] for a in atoms], charge=charge, spin=spin2,
                 basis=BASIS, verbose=0)


def scf_ref(atoms, charge, full):
    """Highest-spin reference (defines the active space). Optional DFT geom-opt."""
    mol = make_mol(atoms, charge, max(SPINS.values()))
    if full:
        done = False
        for backend in ("geometric_solver", "berny_solver"):
            try:
                opt = __import__(f"pyscf.geomopt.{backend}", fromlist=["optimize"]).optimize
                mks = dft.ROKS(mol); mks.xc = "PBE0"; mks.level_shift = 0.2; mks.max_cycle = 300
                mol = opt(mks, maxsteps=80)
                print(f"[geom] relaxed (highest-spin) with {backend}"); done = True; break
            except Exception as e:
                print(f"[warn] {backend} failed ({type(e).__name__}: {e})")
        if not done:
            print("[warn] no optimizer available; using INPUT geometry (vertical gaps only)")
    mf = scf.ROHF(mol); mf.level_shift = 0.3; mf.max_cycle = 400
    mf.kernel()
    if not mf.converged:
        mf = scf.ROHF(mol).newton(); mf.kernel()
    return mf, mol


def cas_for_spin(mf, ncas, nelec_total, spin2, mo):
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
    mc.max_cycle_macro = 250
    e = mc.kernel(mo)[0]
    casdm1 = mc.fcisolver.make_rdm1(mc.ci, ncas, (na, nb))
    occ = np.clip(np.linalg.eigvalsh(casdm1)[::-1], 0.0, 2.0)
    n_u = float(np.sum(occ * (2.0 - occ)))
    return {"e_cas": float(e), "cas_conv": bool(mc.converged), "na": int(na), "nb": int(nb),
            "n_u": round(n_u, 3), "excess": round(n_u - spin2, 3),
            "noon": [round(float(x), 3) for x in occ]}


def vqe_gate(mol, sub_orbitals=4, sub_electrons=4):
    """Engine-correctness gate on a tractable sub-CAS: UCCSD-VQE vs exact (FCI-in-CAS).
       Heavy full-CAS VQE is the actual AWS lift; this confirms the JW mapping + ansatz."""
    try:
        import pennylane as qml
        from pennylane import numpy as pnp
        symbols = [a for a in (mol.elements if hasattr(mol, "elements") else
                               [mol.atom_symbol(i) for i in range(mol.natm)])]
        coords = mol.atom_coords()  # bohr
        H, nq = qml.qchem.molecular_hamiltonian(
            symbols, coords, charge=mol.charge, mult=mol.spin + 1,
            active_electrons=sub_electrons, active_orbitals=sub_orbitals,
            method="pyscf", unit="bohr")
        exact = float(np.linalg.eigvalsh(qml.matrix(H))[0])
        dev = qml.device("lightning.qubit", wires=nq)
        singles, doubles = qml.qchem.excitations(sub_electrons, nq)
        hf = qml.qchem.hf_state(sub_electrons, nq)
        @qml.qnode(dev)
        def cost(p):
            qml.UCCSD(p, wires=range(nq), s_wires=[], d_wires=[], init_state=hf,
                      s_excitations=singles, d_excitations=doubles)
            return qml.expval(H)
        theta = pnp.zeros(len(singles) + len(doubles), requires_grad=True)
        opt = qml.GradientDescentOptimizer(0.4)
        for _ in range(80):
            theta = opt.step(cost, theta)
        e_vqe = float(cost(theta))
        return {"qubits": int(nq), "exact_ha": exact, "vqe_ha": e_vqe,
                "gap_kcal": round((e_vqe - exact) * 627.509, 3)}
    except Exception as e:
        return {"error": f"{type(e).__name__}: {e}"}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--xyz", default=None, help="real cluster geometry (else idealized model)")
    ap.add_argument("--charge", type=int, default=0)
    ap.add_argument("--full", action="store_true", help="DFT geom-opt + VQE gate (AWS)")
    args = ap.parse_args()

    atoms = read_xyz(args.xyz) if args.xyz else idealized_cluster()
    model = args.xyz or "idealized [Fe(O)(OH)2] (IDEALIZED stand-in, not crystallographic)"
    out = {"meta": {"basis": BASIS, "ao_labels": AO_LABELS, "charge": args.charge,
                    "model": model, "full_run": bool(args.full),
                    "known_answer": "alpha-O ground state = S=2 (high-spin Fe(IV)=O), Nature 2016",
                    "honesty": "CASSCF (no dynamic correlation); validates spin-state ORDERING, "
                               "not a quantitative gap; geometry idealized unless --xyz given."},
           "spins": {}}
    try:
        mf, mol = scf_ref(atoms, args.charge, args.full)
        out["meta"]["scf_conv"] = bool(mf.converged)
        ncas, nelec, mo = avas.avas(mf, AO_LABELS, canonicalize=False)
        ncas = int(ncas); nelec = int(nelec)
        out["active_space"] = {"ncas": ncas, "nelec": nelec, "qubits_jw": 2 * ncas}
        print(f"active space: CAS({nelec}e,{ncas}o) = {2*ncas} qubits (JW)")
        for name, s2 in SPINS.items():
            d = cas_for_spin(mf, ncas, nelec, s2, mo)
            out["spins"][name] = d
            print(f"  {name:8s}: E={d.get('e_cas','ERR')}  excess={d.get('excess','?')}")
        ok = {k: v for k, v in out["spins"].items() if "e_cas" in v}
        if ok:
            gs = min(ok, key=lambda k: ok[k]["e_cas"])
            e0 = ok[gs]["e_cas"]
            out["verdict"] = {
                "ground_state_our_casscf": gs,
                "known_answer": "quintet (S=2)",
                "matches_known_S2": bool(gs == KNOWN),
                "rel_kcal_above_ground": {k: round((ok[k]["e_cas"] - e0) * 627.509, 2) for k in ok},
            }
            print("VERDICT:", json.dumps(out["verdict"], ensure_ascii=False))
        if args.full:
            out["vqe_gate"] = vqe_gate(mol)
            print("VQE gate:", json.dumps(out["vqe_gate"], ensure_ascii=False))
    except Exception as e:
        out["error"] = f"{type(e).__name__}: {e}"
        print("ERROR:", out["error"], file=sys.stderr)

    with open("calc/fe_zeolite_alpha_o_results.json", "w") as f:
        json.dump(out, f, indent=2)
    print("wrote calc/fe_zeolite_alpha_o_results.json")


if __name__ == "__main__":
    main()
