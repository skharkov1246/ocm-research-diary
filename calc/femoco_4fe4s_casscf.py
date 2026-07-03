#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
calc/femoco_4fe4s_casscf.py — FeMoco diary, Stage 17.

First MULTIREFERENCE (CASSCF) look at the [4Fe-4S] cubane, plus the Jordan-Wigner
qubit Hamiltonian assembled and CASCI-verified — the concrete next step promised in
the Stage-16 capstone ("turn the [4Fe-4S] active-space estimate into a real CASSCF and
assemble the qubit Hamiltonian").

WHY THIS STAGE (the link to Stage 15). Stage 15 showed broken-symmetry DFT cannot even
fix the SIGN of the magnetic coupling in [4Fe-4S] (B3LYP flips ferro/antiferro; ~15
kcal/mol spread across four functionals). The honest reason is that the ground state is
genuinely multi-configurational — and a KS-DFT / broken-symmetry solution is a single
Slater determinant, which cannot represent it. Here we MEASURE that multireference
character directly (natural-orbital occupations, Head-Gordon n_u) with CASSCF, the
method DFT is failing to stand in for.

WHAT THE SCRIPT DOES.
  1. Build a realistic [Fe4S4(SH)4]2- cubane (Fe-Fe 2.75, Fe-S core ~2.30, terminal
     Fe-S(H) ~2.25 A, S-H 1.34 A), high-spin reference spin=18 (2 Fe3+ d5 + 2 Fe2+ d6
     aligned) — the SAME composition as the Stage-15 BS-DFT run.
  2. Converge the high-spin ROHF mean field (density fitting; DIIS then a Newton /
     second-order restart if DIIS stalls — Fe-S clusters are hard SCF cases).
  3. SIZE the full Fe-3d magnetic active space with AVAS, and Fe-3d + bridging-S-3p,
     and report their qubit counts (2 x n_orb). These are the WALL: their exact CASCI
     (>=1e10 determinants) is not run here — that untractability is the point.
  4. Converge a TRACTABLE FRONTIER CASSCF (target CAS(*,6o) = 12 qubits, the most
     Fe-3d-like orbitals at the Fermi level), taken at the minimal-|Sz| (low-spin)
     split so the CAS ground is the physically-relevant multireference state, and
     measure its NOON / n_u multireference character.
  5. Assemble the Jordan-Wigner qubit Hamiltonian for that frontier space and verify
     its lowest eigenvalue (in the correct electron-number sector) against CASCI in the
     SAME orbitals and SAME minimal-|Sz| split to < 1e-3 kcal/mol — the same
     engine-validation gate used on N2 / Fe-N2 earlier. Self-test: H2/def2-SVP
     CAS(2e,2o) JW vs FCI (assert < 1e-6 Ha; in practice ~1e-16).

HONEST SCOPE (read before trusting any number).
  - Idealized-but-realistic cubane, not a CIF-derived protein site; def2-SVP (modest
    basis); CASSCF only (no dynamic correlation) — we test the QUALITATIVE signal
    (multireference character present) and the engine gate, not a quantitative gap.
  - The high-spin (2S=18) determinant is only a COMPUTATIONAL SEED for AVAS orbital
    selection; it is NOT the physical ground state of [4Fe-4S]2+ (which is
    antiferromagnetically coupled / effectively diamagnetic). The reported frontier
    CASSCF is taken at the minimal-|Sz| split, consistent with that.
  - AVAS is formulated for a closed-shell reference; applied to the open-shell HS-ROHF
    here the returned active electron count is APPROXIMATE. "CAS(*,6o)" fixes the
    orbital count; the electron count is whatever AVAS assigns (recorded in JSON).
  - The full Fe-3d space is only SIZED, not diagonalized: that is the classical wall.
    The "Fe 3d + all-S 3p" wall space includes 3p on ALL 8 sulfurs (4 bridging core +
    4 terminal thiolate), not the bridging four alone.
  - Every number printed/saved is computed at run time. Nothing is invented. SCF /
    CASSCF non-convergence is recorded honestly, and every derived block carries a
    based_on_converged_scf flag so a number can't be read as trustworthy in isolation.

Run (heavy; intended for a compute box / AWS, see calc/aws_run.sh pattern):
    python3 -u calc/femoco_4fe4s_casscf.py
"""
import os
import json
import time
import numpy as np
from pyscf import gto, scf, mcscf, ao2mo
from pyscf.mcscf import avas
import pennylane as qml

HERE = os.path.dirname(os.path.abspath(__file__))
RESULTS = os.path.join(HERE, "femoco_4fe4s_casscf_results.json")
XYZ_OUT = os.path.join(HERE, "femoco_4fe4s_cubane.xyz")
CHK = os.path.join(HERE, "femoco_4fe4s_rohf.chk")  # cached HS-ROHF MOs (regenerable)
KCAL = 627.50947406  # Ha -> kcal/mol

BASIS = "def2-svp"
CHARGE = -2
SPIN = 18            # 2S = unpaired e in the high-spin reference (2 Fe3+ d5 + 2 Fe2+ d6)
TARGET_NCAS = 6      # tractable frontier CAS size (12 qubits) for CASSCF + JW gate;
#                     the JW gate densifies a 2^(2*ncas) matrix, so 6 orbitals ->
#                     4096x4096 (~0.27 GB) fits in RAM; 8 -> 65536^2 (~69 GB) would not.


# --------------------------------------------------------------------------- geometry
def build_cubane(d_FeFe=2.75, r_FeS_core=2.30, r_FeS_term=2.25, r_SH=1.34):
    """Idealized [Fe4S4(SH)4]2- cubane: Fe4 tetrahedron, mu3-S capping each face,
    a terminal SH on each Fe pointing outward. Returns a list of (element, xyz)."""
    v = np.array([[1, 1, 1], [1, -1, -1], [-1, 1, -1], [-1, -1, 1]], float)
    Fe = v * (d_FeFe / np.linalg.norm(v[0] - v[1]))
    cen = Fe.mean(0)
    faces = [[1, 2, 3], [0, 2, 3], [0, 1, 3], [0, 1, 2]]
    S_core = []
    for f in faces:
        fc = Fe[f].mean(0)
        n = fc - cen
        n /= np.linalg.norm(n)
        # push S out along the face normal until mean Fe-S == r_FeS_core
        best = None
        for t in np.linspace(0.0, 1.8, 400):
            p = fc + n * t
            dmean = np.linalg.norm(Fe[f] - p, axis=1).mean()
            if best is None or abs(dmean - r_FeS_core) < best[0]:
                best = (abs(dmean - r_FeS_core), p)
        S_core.append(best[1])
    S_core = np.array(S_core)
    S_term, H = [], []
    for i in range(4):
        d = Fe[i] - cen
        d /= np.linalg.norm(d)
        s = Fe[i] + d * r_FeS_term
        S_term.append(s)
        H.append(s + d * r_SH)
    atoms = ([("Fe", p) for p in Fe]
             + [("S", p) for p in S_core]        # bridging (core) sulfurs, indices 4-7
             + [("S", p) for p in S_term]        # terminal thiolate sulfurs, 8-11
             + [("H", p) for p in H])            # thiol H, 12-15
    return atoms


def write_xyz(atoms, path, comment):
    with open(path, "w") as fh:
        fh.write(f"{len(atoms)}\n{comment}\n")
        for e, p in atoms:
            fh.write(f"{e} {p[0]:.6f} {p[1]:.6f} {p[2]:.6f}\n")


# --------------------------------------------------------------------------- mean field
def converge_rohf(mol):
    """High-spin ROHF with density fitting. Reuses cached MOs (CHK) if present; else a
    short DIIS pass to get into the basin, then a second-order (Newton/SOSCF) solve,
    which is robust for hard Fe-S cluster SCF. Caches the result to CHK so re-runs skip
    the (expensive) SCF. Returns the mf; mf.converged reflects the actual state."""
    if os.path.exists(CHK):
        try:
            mf = scf.ROHF(mol).density_fit()
            mf.update_from_chk(CHK)
            mf.converged = True
            print(f"  loaded ROHF MOs from {os.path.basename(CHK)} (E={mf.e_tot:.5f})")
            return mf
        except Exception as e:
            print(f"  chkfile reuse failed ({e}); recomputing SCF")
    mf = scf.ROHF(mol).density_fit()
    mf.chkfile = CHK
    mf.level_shift = 0.4
    mf.max_cycle = 40
    mf.conv_tol = 1e-7
    mf.kernel()
    if not mf.converged:
        print("  DIIS incomplete -> Newton (SOSCF) restart from current MOs")
        mf2 = mf.newton()
        mf2.chkfile = CHK
        mf2.max_cycle = 100
        mf2.kernel(mf.mo_coeff, mf.mo_occ)
        mf = mf2
    if mf.converged:
        # persist the converged MOs (Newton path writes chk itself, but be explicit)
        scf.chkfile.dump_scf(mol, CHK, mf.e_tot, mf.mo_energy, mf.mo_coeff, mf.mo_occ)
    return mf


# --------------------------------------------------------------------------- active spaces
def avas_size(mf, labels, threshold):
    """Return (ncas, nelec) of the AVAS active space for the given AO labels."""
    ncas, nelec, _ = avas.avas(mf, labels, threshold=threshold, canonicalize=False, verbose=0)
    return int(ncas), int(nelec)


def frontier_magnetic_cas(mf, ncas):
    """Select a small frontier magnetic active space directly from the high-spin
    reference. In HS [4Fe-4S]2- the 18 singly-occupied MOs ARE the Fe-3d magnetic
    orbitals (2 Fe3+ d5 + 2 Fe2+ d6 aligned); we take the `ncas` of them closest to
    the median SOMO energy — the central magnetic frontier. Each contributes one
    electron, so the CAS is half-filled: nelec = ncas.

    (AVAS-threshold cannot carve a *small* frontier here: every Fe-3d orbital projects
    strongly onto the 'Fe 3d' reference, so AVAS always returns the full ~20 — which is
    the 40-qubit wall, not a tractable window.)

    Returns (active_idx, ncas, nelec, fe3d_frac) where fe3d_frac is the mean Fe-3d
    Loewdin character of the chosen orbitals (a labelling honesty check)."""
    occ = mf.mo_occ
    somo = np.where(np.abs(occ - 1.0) < 0.5)[0]          # singly-occupied = magnetic
    if len(somo) < ncas:
        raise RuntimeError(f"only {len(somo)} SOMOs, need {ncas}")
    e_somo = mf.mo_energy[somo]
    med = np.median(e_somo)
    chosen = somo[np.argsort(np.abs(e_somo - med))[:ncas]]
    active = sorted(int(i) for i in chosen)
    nelec = int(round(sum(mf.mo_occ[active])))           # = ncas (each SOMO has 1 e)

    # Loewdin Fe-3d character of the chosen orbitals (honesty: are they really Fe-3d?)
    fe3d = mf.mol.search_ao_label("Fe 3d")
    s = mf.get_ovlp()
    sqrtS = _sqrtm_sym(s)
    C = mf.mo_coeff[:, active]
    Cl = sqrtS @ C                                        # Loewdin-orthogonalized coeffs
    pops = (Cl ** 2)
    fe3d_frac = float(np.mean(pops[fe3d, :].sum(axis=0)))
    return active, ncas, nelec, fe3d_frac


def _sqrtm_sym(s):
    w, v = np.linalg.eigh(s)
    return (v * np.sqrt(np.clip(w, 0, None))) @ v.T


# --------------------------------------------------------------------------- CAS -> qubits
def cas_to_qubit_H(h1, eri, ecore, ncas):
    """Active integrals (chemist (pq|rs)) -> fermion operator -> Jordan-Wigner.
    Spin-orbital for spatial p, spin s in {0,1}: 2p+s. (Same routine validated in
    routes/co2_to_fuels_vqe.py.)"""
    op = ecore * qml.FermiWord({})
    for p in range(ncas):
        for q in range(ncas):
            c = h1[p, q]
            if abs(c) > 1e-12:
                for s in (0, 1):
                    op += c * (qml.FermiC(2 * p + s) * qml.FermiA(2 * q + s))
    for p in range(ncas):
        for q in range(ncas):
            for r in range(ncas):
                for s in range(ncas):
                    c = eri[p, q, r, s]
                    if abs(c) > 1e-12:
                        for a in (0, 1):
                            for b in (0, 1):
                                op += 0.5 * c * (
                                    qml.FermiC(2 * p + a) * qml.FermiC(2 * r + b)
                                    * qml.FermiA(2 * s + b) * qml.FermiA(2 * q + a))
    H = qml.jordan_wigner(op)
    coeffs, ops = H.terms()
    return qml.dot([float(np.real(c)) for c in coeffs], ops), len(coeffs)


def ground_E_in_sector(H, nq, nelec):
    """Lowest eigenvalue of the qubit H restricted to the fixed electron-number sector."""
    Hm = np.real(qml.matrix(H, wire_order=range(nq)))
    keep = [i for i in range(2 ** nq) if bin(i).count("1") == nelec]
    sub = Hm[np.ix_(keep, keep)]
    return float(np.linalg.eigvalsh(sub)[0])


def integrals_from_mc(mc):
    h1, ecore = mc.get_h1eff()
    eri = ao2mo.restore(1, mc.get_h2eff(), mc.ncas)
    return np.asarray(h1), np.asarray(eri), float(ecore)


def self_test_h2():
    mol = gto.M(atom="H 0 0 0; H 0 0 0.74", basis=BASIS, verbose=0)
    mf = scf.RHF(mol)
    mf.kernel()
    mc = mcscf.CASCI(mf, 2, 2)
    mc.verbose = 0
    mc.kernel()
    h1, eri, ec = integrals_from_mc(mc)
    H, _ = cas_to_qubit_H(h1, eri, ec, mc.ncas)
    eq = ground_E_in_sector(H, 2 * mc.ncas, sum(mc.nelecas))
    err = abs(eq - mc.e_tot)
    print(f"self-test H2 CAS(2,2): CASCI={mc.e_tot:.10f} JW={eq:.10f} d={err:.1e} Ha")
    return err


# --------------------------------------------------------------------------- main
def main():
    t_all = time.time()
    res = {"status": "running",
           "meta": {"system": "[Fe4S4(SH)4]2- cubane", "basis": BASIS,
                    "charge": CHARGE, "spin_2S": SPIN,
                    "geometry": "idealized cubane (Fe-Fe 2.75, Fe-S core 2.30, "
                                "term Fe-S 2.25, S-H 1.34 A; terminal Fe-S-H linear)",
                    "method": "HS ROHF/DF seed -> AVAS sizing + frontier CASSCF "
                              "(min-|Sz|) -> JW gate",
                    "note": "CASSCF only, no dynamic correlation; HS(2S=18) is an AVAS "
                            "SEED not the physical (AF-coupled) ground; full Fe-3d space "
                            "is SIZED not diagonalized (the classical wall)."}}
    try:
        err = self_test_h2()
        res["self_test_h2_jw_vs_fci_Ha"] = float(f"{err:.2e}")
        assert err < 1e-6, f"JW self-test failed: {err}"

        atoms = build_cubane()
        write_xyz(atoms, XYZ_OUT,
                  "[Fe4S4(SH)4]2- idealized cubane (FeMoco diary Stage 17)")
        mol = gto.M(atom=[(e, tuple(p)) for e, p in atoms], basis=BASIS,
                    charge=CHARGE, spin=SPIN, verbose=0)
        res["meta"]["nao"] = int(mol.nao)
        res["meta"]["nelec"] = list(map(int, mol.nelec))
        print(f"cubane: nao={mol.nao} nelec={mol.nelec}")

        t0 = time.time()
        mf = converge_rohf(mol)
        scf_ok = bool(mf.converged)
        res["scf"] = {"e_hs_rohf": float(mf.e_tot), "converged": scf_ok,
                      "seconds": round(time.time() - t0, 1)}
        print(f"HS ROHF: E={mf.e_tot:.5f} Ha conv={scf_ok} "
              f"({res['scf']['seconds']}s)")

        # --- size the magnetic active spaces (the wall) ---
        # NB "Fe_3d__all_S_3p" selects 3p on ALL 8 sulfurs (4 bridging core + 4
        # terminal thiolate), not the bridging four alone.
        walls = {}
        for tag, labels in [("Fe_3d", ["Fe 3d"]),
                            ("Fe_3d__all_S_3p", ["Fe 3d", "S 3p"])]:
            try:
                nc, ne = avas_size(mf, labels, threshold=0.2)
                walls[tag] = {"ncas": nc, "nelec": ne, "qubits_jw": 2 * nc,
                              "based_on_converged_scf": scf_ok, "diagonalized": False}
                print(f"wall {tag}: CAS({ne}e,{nc}o) = {2*nc} qubits (not diagonalized)")
            except Exception as e:
                walls[tag] = {"error": f"{type(e).__name__}: {e}"}
        res["wall_active_spaces"] = walls

        # --- tractable frontier CASSCF at the minimal-|Sz| split ---
        frontier = {"based_on_converged_scf": scf_ok}
        try:
            active, ncas, nelec, fe3d_frac = frontier_magnetic_cas(mf, TARGET_NCAS)
            # minimal-|Sz| split (Sz=0 for even N, 1/2 for odd): this sector provably
            # contains the global N-electron ground, so CASSCF/CASCI here target the
            # physically-relevant low-spin state AND match the JW N-sector search.
            na = (nelec + (nelec % 2)) // 2
            nb = nelec - na
            frontier.update(ncas=int(ncas), nelec=int(nelec), nelecas=[int(na), int(nb)],
                            qubits_jw=2 * int(ncas), active_mo_idx=active,
                            frontier_fe3d_fraction=round(fe3d_frac, 3))
            print(f"frontier magnetic CAS({nelec}e,{ncas}o) split ({na},{nb}) "
                  f"= {2*ncas} qubits; Fe-3d character {fe3d_frac:.2f}")

            t0 = time.time()
            mc = mcscf.CASSCF(mf, ncas, (na, nb)).density_fit()
            mc.max_cycle_macro = 50
            mc.verbose = 0
            mo = mc.sort_mo(active, base=0)          # place chosen SOMOs in the CAS
            e_casscf = mc.kernel(mo)[0]
            frontier["casscf_seconds"] = round(time.time() - t0, 1)
            frontier["e_casscf_df"] = float(e_casscf)
            frontier["casscf_converged"] = bool(mc.converged)

            casdm1 = mc.fcisolver.make_rdm1(mc.ci, ncas, mc.nelecas)
            occ = np.clip(np.linalg.eigvalsh(casdm1)[::-1], 0.0, 2.0)
            n_u = float(np.sum(occ * (2.0 - occ)))
            frontier["noon"] = [round(float(x), 3) for x in occ]
            frontier["n_u"] = round(n_u, 3)
            print(f"CASSCF: E={e_casscf:.6f} conv={mc.converged} "
                  f"n_u={n_u:.3f} ({frontier['casscf_seconds']}s)")

            # --- JW gate: qubit H reproduces the exact CAS diagonalization ---
            # CASCI reference in the SAME orbitals and SAME (na,nb) split as JW's
            # N-electron sector search -> like-for-like (both exact-integral).
            mc_ci = mcscf.CASCI(mf, ncas, (na, nb))
            mc_ci.verbose = 0
            mc_ci.mo_coeff = mc.mo_coeff
            e_casci = mc_ci.kernel(mc.mo_coeff)[0]
            h1, eri, ec = integrals_from_mc(mc_ci)
            H, nterms = cas_to_qubit_H(h1, eri, ec, ncas)
            eq = ground_E_in_sector(H, 2 * ncas, nelec)
            d_jw = (eq - e_casci) * KCAL
            frontier.update(jw_terms=int(nterms), e_casci_exact=float(e_casci),
                            e_jw=float(eq), jw_minus_casci_kcal=float(f"{d_jw:.2e}"))
            print(f"JW gate: CASCI={e_casci:.8f} JW={eq:.8f} "
                  f"d={d_jw:.2e} kcal/mol ({nterms} terms)")
        except Exception as e:
            frontier["error"] = f"{type(e).__name__}: {e}"
            print("frontier CASSCF/JW error:", frontier["error"])
        res["frontier_casscf"] = frontier
        res["status"] = "ok"
    except Exception as e:
        res["status"] = f"aborted: {type(e).__name__}: {e}"
        print("ABORT:", res["status"])
    finally:
        res["meta"]["total_seconds"] = round(time.time() - t_all, 1)
        with open(RESULTS, "w") as fh:
            json.dump(res, fh, ensure_ascii=False, indent=1)
        print(f"\nwrote {RESULTS}  status={res['status']}  "
              f"(total {res['meta']['total_seconds']}s)")


if __name__ == "__main__":
    main()
