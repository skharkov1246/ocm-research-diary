#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
calc/fe4s4_casscf_aws.py — FeMoco diary, Stage 19 (AWS job): the REAL CASSCF.

Stage 17 measured the multireference character of the [4Fe-4S] cubane but its
CASSCF never converged; Stage 18 diagnosed WHY: the high-spin-ROHF frontier
orbitals are a terrible starting point (the six-electron frontier configuration
sits ~10 Ha above the HS determinant, and CASSCF must grind through that whole
orbital relaxation, ~2 min/macro-iteration with the gradient only halving per
step). The textbook fix — and this job — is to START FROM THE RIGHT ORBITALS:

  1. Broken-symmetry (BS) mean field. Converge the high-spin (2S=18) UKS-PBE
     ferromagnetic state (easy), then FLIP the spin density on two of the four
     Fe and reconverge at Ms=0. Flip pattern (documented, recorded in JSON):
         Fe0 up, Fe1 up | Fe2 down, Fe3 down
     i.e. the two-mixed-valence-pair (2 x Fe(2.5+)-Fe(2.5+)) Noodleman BS state
     of the [Fe4S4]2+ core — the standard BS solution for this cluster. The
     flip is done by swapping the alpha/beta blocks of the converged HS density
     matrix on the AOs of Fe2 and Fe3 (diagonal atomic blocks; the standard
     approximation) and feeding that as dm0 to a spin=0 UKS run.
  2. UNO. Natural orbitals of the total (alpha+beta) BS density. Their
     fractionally-occupied orbitals ARE the magnetic Fe-3d / bridging-S-3p
     space — this replaces AVAS (which, applied to the open-shell reference,
     Stage 17 showed cannot carve a small window: every Fe-3d orbital projects
     strongly). We still verify the AVAS-style targeting a posteriori: the
     Loewdin Fe-3d + bridging-S-3p character of every chosen orbital is
     computed and recorded (honesty check that the CAS is what we claim).
  3. Ladder CASSCF: CAS(*,10o) -> (*,12o) -> (*,14o). Each rung starts from
     the previous rung's CONVERGED orbitals (expanding the window by the
     highest-magnetic-character core + virtual orbital); the next rung runs
     only if the previous one converged AND the remaining wall-clock budget
     allows. Electron counts come from the UNO occupations (recorded), never
     assumed.
  4. Qubit Hamiltonian for the LARGEST CONVERGED rung: exact (non-DF) active
     integrals -> fermion operator -> Jordan-Wigner (openfermion; pennylane
     fallback), reporting qubit and Pauli-term counts. Engine gate: the same
     builder is self-tested on H2 CAS(2,2) vs PySCF CASCI (< 1e-6 Ha) before
     any big build. The 20-28-qubit Hamiltonian itself is NOT diagonalized
     here — the number sector alone is ~1.8e5-4e7 dimensional; that wall is
     the point of the whole exercise.

HONEST SCOPE.
  - Same idealized [Fe4S4(SH)4]2- cubane geometry as Stages 17/18 (so energies
    are comparable); def2-SVP; CASSCF only (no dynamic correlation here).
    FE4S4_LIGAND=SCH3 switches to the classic [Fe4S4(SCH3)4]2- Holm-type model
    (then energies are NOT comparable to Stage 17/18 — flagged in JSON).
  - Every number is computed at run time; non-convergence, fallbacks and
    timeouts are recorded as such. Derived blocks carry the convergence flags
    of what they were computed from.
  - signal.alarm stage timeouts are best-effort (a long C-level step returns
    to Python before the exception fires); the shell-level `timeout` in
    calc/aws_run_fe4s4.sh is the hard stop.

ENV KNOBS (all optional):
  FE4S4_STAGE_TIMEOUT  seconds per stage (SCF / each CASSCF rung / JW build),
                       default 5400
  FE4S4_TOTAL_BUDGET   total wall-clock budget in seconds, default 19800
                       (leaves headroom inside the 6 h shell cap)
  FE4S4_LADDER         comma list of active-space sizes, default "10,12,14"
                       (must grow in steps of 2)
  FE4S4_MF             "UKS" (PBE, default) or "UHF" mean field
  FE4S4_LIGAND         "SH" (default, = Stage 17/18) or "SCH3" (classic model)
  FE4S4_MACRO          CASSCF max macro-iterations per rung, default 60

OUTPUT (all incremental / atomic — nothing is lost on a kill):
  calc/fe4s4_casscf_results.json    rolling results (os.replace atomic writes)
  calc/fe4s4_casscf_integrals.npz   active-space integrals of the best rung
  calc/fe4s4_cubane_aws.xyz         the geometry actually used

Run (heavy; meant for AWS via calc/aws_run_fe4s4.sh):
    python3 -u calc/fe4s4_casscf_aws.py
"""
import os
import sys
import json
import time
import signal

try:
    import numpy as np
except ImportError:
    sys.exit("[deps] numpy is required:  python3 -m pip install numpy")
try:
    from pyscf import gto, scf, dft, mcscf, ao2mo
except ImportError as e:
    sys.exit(f"[deps] pyscf is required:  python3 -m pip install pyscf   ({e})")

HERE = os.path.dirname(os.path.abspath(__file__))
RESULTS = os.path.join(HERE, "fe4s4_casscf_results.json")
NPZ = os.path.join(HERE, "fe4s4_casscf_integrals.npz")
XYZ_OUT = os.path.join(HERE, "fe4s4_cubane_aws.xyz")
KCAL = 627.50947406

BASIS = "def2-svp"
CHARGE = -2
SPIN_HS = 18                    # HS seed: 2 Fe3+ (d5) + 2 Fe2+ (d6), all aligned
FLIP_ATOMS = (2, 3)             # BS pattern: Fe0 up, Fe1 up | Fe2 down, Fe3 down
BRIDGE_S_ATOMS = (4, 5, 6, 7)   # mu3-S core sulfurs (atom order of build_cubane)

STAGE_TIMEOUT = int(os.environ.get("FE4S4_STAGE_TIMEOUT", "5400"))
TOTAL_BUDGET = int(os.environ.get("FE4S4_TOTAL_BUDGET", "19800"))
LADDER = [int(x) for x in os.environ.get("FE4S4_LADDER", "10,12,14").split(",")]
MF_KIND = os.environ.get("FE4S4_MF", "UKS").upper()
LIGAND = os.environ.get("FE4S4_LIGAND", "SH").upper()
MACRO = int(os.environ.get("FE4S4_MACRO", "60"))

T0 = time.time()


# ------------------------------------------------------------------ utilities
def elapsed():
    return time.time() - T0


def remaining():
    return TOTAL_BUDGET - elapsed()


def save(res):
    """Atomic incremental save: write to tmp, then os.replace."""
    tmp = RESULTS + ".tmp"
    with open(tmp, "w") as fh:
        json.dump(res, fh, ensure_ascii=False, indent=1)
    os.replace(tmp, RESULTS)


class StageTimeout(Exception):
    pass


def _alarm_handler(signum, frame):
    raise StageTimeout()


def with_timeout(seconds, fn, *args, **kwargs):
    """Best-effort per-stage wall clock cap via SIGALRM (main thread only).
    A long C-level step delays the exception until control returns to Python;
    the shell-level `timeout` in aws_run_fe4s4.sh is the hard backstop."""
    seconds = int(seconds)
    if seconds <= 0:
        raise StageTimeout()
    old = signal.signal(signal.SIGALRM, _alarm_handler)
    signal.alarm(seconds)
    try:
        return fn(*args, **kwargs)
    finally:
        signal.alarm(0)
        signal.signal(signal.SIGALRM, old)


def _sqrtm_sym(s):
    w, v = np.linalg.eigh(s)
    w = np.clip(w, 1e-12, None)
    return (v * np.sqrt(w)) @ v.T, (v / np.sqrt(w)) @ v.T   # S^1/2, S^-1/2


# ------------------------------------------------------------------- geometry
def build_cubane(ligand="SH", d_FeFe=2.75, r_FeS_core=2.30, r_FeS_term=2.25,
                 r_SH=1.34, r_SC=1.83, r_CH=1.09):
    """Idealized [Fe4S4(SR)4]2- cubane, same core as Stage 17
    (calc/femoco_4fe4s_casscf.py): Fe4 tetrahedron, mu3-S capping each face,
    a terminal SR on each Fe pointing outward. R = H (Stage-17/18-compatible)
    or CH3 (classic Holm-type [Fe4S4(SCH3)4]2-). Atom order: Fe 0-3,
    bridging S 4-7, terminal S 8-11, then R atoms."""
    v = np.array([[1, 1, 1], [1, -1, -1], [-1, 1, -1], [-1, -1, 1]], float)
    Fe = v * (d_FeFe / np.linalg.norm(v[0] - v[1]))
    cen = Fe.mean(0)
    faces = [[1, 2, 3], [0, 2, 3], [0, 1, 3], [0, 1, 2]]
    S_core = []
    for f in faces:
        fc = Fe[f].mean(0)
        n = fc - cen
        n /= np.linalg.norm(n)
        best = None
        for t in np.linspace(0.0, 1.8, 400):
            p = fc + n * t
            dmean = np.linalg.norm(Fe[f] - p, axis=1).mean()
            if best is None or abs(dmean - r_FeS_core) < best[0]:
                best = (abs(dmean - r_FeS_core), p)
        S_core.append(best[1])
    atoms = [("Fe", p) for p in Fe] + [("S", p) for p in S_core]
    S_term, tails = [], []
    for i in range(4):
        d = Fe[i] - cen
        d /= np.linalg.norm(d)
        s = Fe[i] + d * r_FeS_term
        S_term.append(("S", s))
        if ligand == "SH":
            tails.append(("H", s + d * r_SH))
        elif ligand == "SCH3":
            c = s + d * r_SC
            t = np.array([1.0, 0, 0]) if abs(d[0]) < 0.9 else np.array([0, 1.0, 0])
            u = np.cross(d, t); u /= np.linalg.norm(u)
            w = np.cross(d, u)
            th = np.deg2rad(180.0 - 109.47)      # S-C-H = 109.47 deg
            for ph in (0.0, 2 * np.pi / 3, 4 * np.pi / 3):
                hdir = np.cos(th) * d + np.sin(th) * (np.cos(ph) * u + np.sin(ph) * w)
                tails.append(("H", c + r_CH * hdir))
            tails.insert(len(tails) - 3, ("C", c))
        else:
            raise ValueError(f"unknown FE4S4_LIGAND={ligand!r} (use SH or SCH3)")
    return atoms + S_term + tails


def write_xyz(atoms, path, comment):
    with open(path, "w") as fh:
        fh.write(f"{len(atoms)}\n{comment}\n")
        for e, p in atoms:
            fh.write(f"{e} {p[0]:.6f} {p[1]:.6f} {p[2]:.6f}\n")


# ----------------------------------------------------------------- mean field
def _new_mf(mol):
    if MF_KIND == "UHF":
        return scf.UHF(mol).density_fit()
    mf = dft.UKS(mol).density_fit()
    mf.xc = "pbe"
    return mf


def _converge(mf, dm0=None, level_shift=0.4, max_cycle=60, tag=""):
    """DIIS with level shift; on stall, a second-order (Newton/SOSCF) restart
    from the current MOs — the robust recipe for hard Fe-S SCF cases."""
    mf.level_shift = level_shift
    mf.max_cycle = max_cycle
    mf.conv_tol = 1e-7
    mf.kernel(dm0=dm0) if dm0 is not None else mf.kernel()
    if not mf.converged:
        print(f"  [{tag}] DIIS incomplete -> Newton (SOSCF) restart", flush=True)
        mf2 = mf.newton()
        mf2.max_cycle = 100
        mf2.kernel(mf.mo_coeff, mf.mo_occ)
        mf2.mol = mf.mol
        return mf2
    return mf


def flip_dm(dm, mol, flip_atoms):
    """Swap the alpha/beta density on the diagonal AO blocks of `flip_atoms`
    (the standard BS spin-flip construction; off-diagonal blocks untouched)."""
    dma, dmb = dm[0].copy(), dm[1].copy()
    aoslice = mol.aoslice_by_atom()
    for a in flip_atoms:
        sl = slice(int(aoslice[a][2]), int(aoslice[a][3]))
        blk_a = dma[sl, sl].copy()
        dma[sl, sl] = dmb[sl, sl]
        dmb[sl, sl] = blk_a
    return np.array([dma, dmb])


def fe_spin_pops(mf, natoms=4):
    """Mulliken spin population on each Fe (sanity check of the BS pattern)."""
    dma, dmb = mf.make_rdm1()
    m = (dma - dmb) @ mf.get_ovlp()
    aoslice = mf.mol.aoslice_by_atom()
    out = []
    for a in range(natoms):
        sl = slice(int(aoslice[a][2]), int(aoslice[a][3]))
        out.append(round(float(np.trace(m[sl, sl])), 3))
    return out


def converge_reference(mol_hs, mol_bs, info):
    """HS seed -> spin-flip -> BS mean field. Returns (mf, mode) where mode is
    'BS' or 'HS-fallback'. Every step's convergence is recorded in `info`."""
    print(f"[scf] HS seed ({MF_KIND}, 2S={SPIN_HS})", flush=True)
    t0 = time.time()
    mf_hs = with_timeout(min(STAGE_TIMEOUT, remaining()), _converge,
                         _new_mf(mol_hs), tag="HS")
    info["hs"] = {"e": float(mf_hs.e_tot), "converged": bool(mf_hs.converged),
                  "seconds": round(time.time() - t0, 1)}
    print(f"[scf] HS: E={mf_hs.e_tot:.5f} conv={mf_hs.converged} "
          f"({info['hs']['seconds']}s)", flush=True)

    print(f"[scf] BS: flip alpha/beta on Fe{FLIP_ATOMS} "
          f"(pattern Fe0+ Fe1+ | Fe2- Fe3-), Ms=0", flush=True)
    dm0 = flip_dm(mf_hs.make_rdm1(), mol_hs, FLIP_ATOMS)
    t0 = time.time()
    mf_bs = with_timeout(min(STAGE_TIMEOUT, remaining()), _converge,
                         _new_mf(mol_bs), dm0=dm0, level_shift=0.3, tag="BS")
    ss, mult = mf_bs.spin_square()
    pops = fe_spin_pops(mf_bs)
    info["bs"] = {"e": float(mf_bs.e_tot), "converged": bool(mf_bs.converged),
                  "seconds": round(time.time() - t0, 1),
                  "s_squared": round(float(ss), 3),
                  "fe_mulliken_spin": pops,
                  "flip_atoms": list(FLIP_ATOMS)}
    print(f"[scf] BS: E={mf_bs.e_tot:.5f} conv={mf_bs.converged} "
          f"<S^2>={ss:.2f} Fe spins={pops} ({info['bs']['seconds']}s)", flush=True)

    # BS quality gate: converged AND genuinely spin-polarized in the +,+,-,-
    # pattern (|m_Fe| well above noise). Otherwise fall back honestly.
    signs_ok = (len(pops) == 4 and pops[0] > 0 and pops[1] > 0
                and pops[2] < 0 and pops[3] < 0)
    polarized = all(abs(p) > 1.5 for p in pops)
    info["bs"]["pattern_ok"] = bool(signs_ok and polarized)
    if mf_bs.converged and signs_ok and polarized:
        return mf_bs, "BS"
    print("[scf][WARN] BS solution not usable (unconverged or collapsed) — "
          "falling back to the HS reference (Stage-17-like orbitals; the CAS "
          "start will be poorer, recorded as such)", flush=True)
    return mf_hs, "HS-fallback"


# --------------------------------------------------------- UNO + active space
def unos(mf):
    """Natural orbitals of the total (alpha+beta) density:
    diagonalize S^1/2 D S^1/2; return (occ desc, C columns to match)."""
    s = mf.get_ovlp()
    shalf, sminus = _sqrtm_sym(s)
    dma, dmb = mf.make_rdm1()
    n, u = np.linalg.eigh(shalf @ (dma + dmb) @ shalf)
    order = np.argsort(n)[::-1]
    return np.clip(n[order], 0.0, 2.0), sminus @ u[:, order]


def loewdin_frac(mol, C, ao_idx, shalf=None):
    """Loewdin population of the given AO set in each column of C."""
    if shalf is None:
        shalf, _ = _sqrtm_sym(mol.intor("int1e_ovlp"))
    pops = (shalf @ C) ** 2
    return pops[ao_idx, :].sum(axis=0)


def magnetic_ao_idx(mol):
    """AO indices of Fe 3d plus BRIDGING S 3p (atoms 4-7 only — the mu3 core
    sulfurs, not the terminal thiolates): the AVAS-style target set."""
    idx = []
    for i, lab in enumerate(mol.ao_labels()):
        parts = lab.split()
        atom, shell = int(parts[0]), parts[2]
        if parts[1] == "Fe" and shell.startswith("3d"):
            idx.append(i)
        elif parts[1] == "S" and atom in BRIDGE_S_ATOMS and shell.startswith("3p"):
            idx.append(i)
    return np.array(idx, dtype=int)


def select_active(mol, occ, C, ncas):
    """Choose the `ncas` UNOs closest to single occupation (|n-1| minimal) —
    the magnetic space of the BS solution. Electron count = rounded sum of
    their occupations (forced even so ncore is integral and na=nb).
    Returns (mo reordered core|active|virtual, nelec, diagnostics)."""
    rank = np.argsort(np.abs(occ - 1.0))
    sel = sorted(rank[:ncas].tolist(), key=lambda i: -occ[i])
    ssum = float(sum(occ[i] for i in sel))
    nelec = int(round(ssum))
    if (mol.nelectron - nelec) % 2:
        nelec += 1 if ssum > nelec else -1
    ncore = (mol.nelectron - nelec) // 2
    rest = [i for i in np.argsort(occ)[::-1].tolist() if i not in set(sel)]
    core, virt = rest[:ncore], rest[ncore:]
    mo = C[:, core + sel + virt]
    mag = magnetic_ao_idx(mol)
    frac = loewdin_frac(mol, C[:, sel], mag)
    diag = {"active_uno_occ": [round(float(occ[i]), 3) for i in sel],
            "occ_sum": round(ssum, 3), "nelec": nelec, "ncore": ncore,
            "fe3d_bridgeS3p_frac": [round(float(x), 3) for x in frac],
            "min_core_occ": round(float(occ[core[-1]]), 3) if core else None,
            "max_virt_occ": round(float(occ[virt[0]]), 3) if virt else None}
    return mo, nelec, diag


def expand_active(mol, mc, add=2):
    """Grow a converged CAS by `add`=2: promote the core and the virtual
    orbital with the highest Fe-3d + bridging-S-3p Loewdin character into the
    window, keeping the rung's relaxed orbitals (no restart from scratch)."""
    mo, ncore, ncas = mc.mo_coeff, mc.ncore, mc.ncas
    nmo = mo.shape[1]
    mag = magnetic_ao_idx(mol)
    frac = loewdin_frac(mol, mo, mag)
    best_core = int(np.argmax(frac[:ncore]))
    best_virt = int(ncore + ncas + np.argmax(frac[ncore + ncas:]))
    order = ([i for i in range(ncore) if i != best_core] + [best_core]
             + list(range(ncore, ncore + ncas)) + [best_virt]
             + [i for i in range(ncore + ncas, nmo) if i != best_virt])
    nelec = sum(mc.nelecas) + 2
    diag = {"promoted_core_frac": round(float(frac[best_core]), 3),
            "promoted_virt_frac": round(float(frac[best_virt]), 3),
            "nelec": int(nelec)}
    return mo[:, order], int(nelec), diag


# --------------------------------------------------------------------- CASSCF
def run_rung(mol, mo, ncas, nelec, alarm_s):
    """One CASSCF ladder rung (DF-accelerated), timeout-guarded."""
    na = nb = nelec // 2
    mf_r = scf.RHF(mol).density_fit()          # orbital container; mo is given
    mc = mcscf.CASSCF(mf_r, ncas, (na, nb)).density_fit()
    mc.max_cycle_macro = MACRO
    mc.verbose = 0
    info = {"ncas": int(ncas), "nelec": int(nelec), "nelecas": [na, nb],
            "qubits_jw": 2 * int(ncas)}
    t0 = time.time()
    try:
        with_timeout(alarm_s, mc.kernel, mo)
        info["converged"] = bool(mc.converged)
        info["e_casscf_df"] = float(mc.e_tot)
        dm1 = mc.fcisolver.make_rdm1(mc.ci, ncas, mc.nelecas)
        noon = np.clip(np.linalg.eigvalsh(dm1)[::-1], 0.0, 2.0)
        info["noon"] = [round(float(x), 3) for x in noon]
        info["n_u"] = round(float(np.sum(noon * (2.0 - noon))), 3)
        try:
            ss = mc.fcisolver.spin_square(mc.ci, ncas, mc.nelecas)[0]
            info["cas_s_squared"] = round(float(ss), 3)
        except Exception:
            pass
        mag = magnetic_ao_idx(mol)
        frac = loewdin_frac(mol, mc.mo_coeff[:, mc.ncore:mc.ncore + ncas], mag)
        info["active_fe3d_bridgeS3p_mean"] = round(float(frac.mean()), 3)
    except StageTimeout:
        info["converged"] = False
        info["timed_out_after_s"] = int(alarm_s)
        print(f"  [rung {ncas}o] TIMED OUT after {alarm_s}s", flush=True)
    info["seconds"] = round(time.time() - t0, 1)
    return mc, info


# ------------------------------------------------------------ qubit Hamiltonian
def qubit_backend():
    try:
        import openfermion  # noqa: F401
        return "openfermion"
    except ImportError:
        pass
    try:
        import pennylane  # noqa: F401
        return "pennylane"
    except ImportError:
        return None


def _fermionic_terms(h1, eri, ncas, tol=1e-12):
    """Spin-orbital second-quantized terms from active integrals (chemist
    (pq|rs)); spin-orbital index for spatial p, spin s in {0,1} is 2p+s —
    the same convention validated in calc/femoco_4fe4s_casscf.py (Stage 17)
    and routes/co2_to_fuels_vqe.py."""
    for p in range(ncas):
        for q in range(ncas):
            c = h1[p, q]
            if abs(c) > tol:
                for s in (0, 1):
                    yield float(c), ((2 * p + s, 1), (2 * q + s, 0))
    for p in range(ncas):
        for q in range(ncas):
            for r in range(ncas):
                for s_ in range(ncas):
                    c = eri[p, q, r, s_]
                    if abs(c) > tol:
                        for a in (0, 1):
                            for b in (0, 1):
                                yield 0.5 * float(c), ((2 * p + a, 1),
                                                       (2 * r + b, 1),
                                                       (2 * s_ + b, 0),
                                                       (2 * q + a, 0))


def build_qubit_H(h1, eri, ecore, ncas, backend):
    """Return (qubit Hamiltonian object, number of Pauli terms)."""
    if backend == "openfermion":
        import openfermion as of
        op = of.FermionOperator((), float(ecore))
        for c, term in _fermionic_terms(h1, eri, ncas):
            op += of.FermionOperator(term, c)
        H = of.transforms.jordan_wigner(op)
        H.compress(1e-12)
        return H, len(H.terms)
    if backend == "pennylane":
        import pennylane as qml
        op = float(ecore) * qml.FermiWord({})
        for c, term in _fermionic_terms(h1, eri, ncas):
            w = None
            for idx, dag in term:
                f = qml.FermiC(idx) if dag else qml.FermiA(idx)
                w = f if w is None else w * f
            op += c * w
        H = qml.jordan_wigner(op)
        coeffs, ops = H.terms()
        return qml.dot([float(np.real(c)) for c in coeffs], ops), len(coeffs)
    raise RuntimeError("no qubit backend (install openfermion or pennylane)")


def ground_E_in_sector(H, nq, nelec, backend):
    """Dense sector-restricted ground energy — used ONLY for the tiny H2
    self-test (nq=4); never called on the big cluster Hamiltonians."""
    if backend == "openfermion":
        import openfermion as of
        Hm = np.real(of.linalg.get_sparse_operator(H, n_qubits=nq).toarray())
    else:
        import pennylane as qml
        Hm = np.real(qml.matrix(H, wire_order=range(nq)))
    keep = [i for i in range(2 ** nq) if bin(i).count("1") == nelec]
    return float(np.linalg.eigvalsh(Hm[np.ix_(keep, keep)])[0])


def self_test_h2(backend):
    """H2 CAS(2,2): JW ground energy must reproduce PySCF CASCI (<1e-6 Ha)."""
    mol = gto.M(atom="H 0 0 0; H 0 0 0.74", basis="sto-3g", verbose=0)
    mf = scf.RHF(mol)
    mf.kernel()
    mc = mcscf.CASCI(mf, 2, 2)
    mc.verbose = 0
    mc.kernel()
    h1, ecore = mc.get_h1eff()
    eri = ao2mo.restore(1, mc.get_h2eff(), 2)
    H, _ = build_qubit_H(h1, eri, ecore, 2, backend)
    err = abs(ground_E_in_sector(H, 4, 2, backend) - mc.e_tot)
    print(f"[gate] H2 CAS(2,2) JW vs CASCI ({backend}): d={err:.2e} Ha", flush=True)
    return err


def exact_active_integrals(mol, mo, ncas, nelecas):
    """Exact (non-DF) active-space integrals in the converged CASSCF orbitals
    (the CASSCF itself used DF for speed; the qubit H gets exact integrals —
    noted in JSON)."""
    mc_ci = mcscf.CASCI(scf.RHF(mol), ncas, nelecas)
    mc_ci.verbose = 0
    h1, ecore = mc_ci.get_h1eff(mo)
    eri = ao2mo.restore(1, mc_ci.get_h2eff(mo), ncas)
    return np.asarray(h1), np.asarray(eri), float(ecore)


# ------------------------------------------------------------------------ main
def main():
    res = {"status": "running",
           "meta": {"system": f"[Fe4S4({LIGAND})4]2- cubane",
                    "basis": BASIS, "charge": CHARGE,
                    "mean_field": f"{MF_KIND}" + ("" if MF_KIND == "UHF" else "-PBE"),
                    "bs_pattern": "Fe0 up, Fe1 up | Fe2 down, Fe3 down (Ms=0; "
                                  "two mixed-valence Fe2.5+ pairs, Noodleman BS "
                                  "state of the [Fe4S4]2+ core)",
                    "ladder": LADDER, "macro_max": MACRO,
                    "stage_timeout_s": STAGE_TIMEOUT,
                    "total_budget_s": TOTAL_BUDGET,
                    "geometry_note": ("identical idealized cubane core to Stage "
                                      "17/18 (energies comparable)" if LIGAND == "SH"
                                      else "SCH3-terminated classic model — NOT "
                                           "energy-comparable to Stage 17/18"),
                    "note": "BS-UNO start = the fix for Stage 18's wall (CASSCF "
                            "orbital relaxation from a poor HS-ROHF frontier). "
                            "CASSCF only, no dynamic correlation. DF-accelerated "
                            "CASSCF; qubit H built from exact active integrals. "
                            "The big qubit H is counted, not diagonalized."}}
    save(res)
    try:
        if any(b - a != 2 for a, b in zip(LADDER, LADDER[1:])):
            raise ValueError(f"FE4S4_LADDER must grow in steps of 2, got {LADDER}")

        backend = qubit_backend()
        res["qubit_backend"] = backend or ("MISSING — install openfermion "
                                           "(pip install openfermion); CASSCF "
                                           "results below are unaffected")
        if backend:
            err = self_test_h2(backend)
            res["self_test_h2_jw_vs_casci_Ha"] = float(f"{err:.2e}")
            assert err < 1e-6, f"JW self-test failed: {err} Ha"
        save(res)

        atoms = build_cubane(LIGAND)
        write_xyz(atoms, XYZ_OUT, f"[Fe4S4({LIGAND})4]2- idealized cubane "
                                  "(FeMoco diary Stage 19, AWS)")
        maxmem = int(os.environ.get("PYSCF_MAX_MEMORY", "16000"))
        mol_kw = dict(atom=[(e, tuple(p)) for e, p in atoms], basis=BASIS,
                      charge=CHARGE, verbose=0, max_memory=maxmem)
        mol_hs = gto.M(spin=SPIN_HS, **mol_kw)
        mol_bs = gto.M(spin=0, **mol_kw)
        res["meta"]["nao"] = int(mol_bs.nao)
        res["meta"]["nelec_total"] = int(mol_bs.nelectron)
        print(f"[mol] nao={mol_bs.nao} nelec={mol_bs.nelectron} "
              f"maxmem={maxmem}MB", flush=True)
        save(res)

        # --- mean field: HS seed -> spin flip -> BS (with honest fallback) ---
        scf_info = {}
        mf, mode = converge_reference(mol_hs, mol_bs, scf_info)
        if mode == "HS-fallback" and MF_KIND == "UKS" and remaining() > 2 * STAGE_TIMEOUT:
            print("[scf] retrying the whole HS->BS chain with UHF", flush=True)
            globals()["MF_KIND"] = "UHF"        # recorded fallback, not silent
            scf_info["uhf_retry"] = {}
            mf, mode = converge_reference(mol_hs, mol_bs, scf_info["uhf_retry"])
        scf_info["reference_used"] = mode
        res["scf"] = scf_info
        save(res)
        ref_mol = mol_bs if mode == "BS" else mol_hs

        # --- UNO spectrum ---
        occ, C = unos(mf)
        frac_mask = (occ > 0.02) & (occ < 1.98)
        res["uno"] = {"n_fractional_0.02_1.98": int(frac_mask.sum()),
                      "fractional_occs": [round(float(x), 3)
                                          for x in occ[frac_mask][:40]],
                      "from_reference": mode}
        print(f"[uno] {int(frac_mask.sum())} fractionally-occupied naturals "
              f"(the magnetic space)", flush=True)
        save(res)

        # --- ladder CASSCF ---
        rungs = []
        best = None            # (mc, ncas, nelec)
        last_secs = None
        for ncas in LADDER:
            need = max(900, 1.5 * last_secs) if last_secs else 900
            if remaining() < need:
                rungs.append({"ncas": ncas, "skipped": f"time budget: "
                              f"{int(remaining())}s left < {int(need)}s needed"})
                print(f"[ladder] skip {ncas}o — budget", flush=True)
                break
            if best is None:
                mo, nelec, diag = select_active(ref_mol, occ, C, ncas)
            else:
                mo, nelec, diag = expand_active(ref_mol, best[0])
            alarm = int(min(STAGE_TIMEOUT, remaining() - 300))
            print(f"[ladder] CASSCF CAS({nelec}e,{ncas}o) = {2*ncas} qubits, "
                  f"alarm {alarm}s", flush=True)
            mc, info = run_rung(ref_mol, mo, ncas, nelec, alarm)
            info["selection"] = diag
            info["based_on_reference"] = mode
            rungs.append(info)
            res["casscf_ladder"] = rungs
            save(res)
            print(f"[ladder] {ncas}o: conv={info.get('converged')} "
                  f"E={info.get('e_casscf_df')} n_u={info.get('n_u')} "
                  f"({info['seconds']}s)", flush=True)
            if not info.get("converged"):
                break
            best = (mc, ncas, nelec)
            last_secs = info["seconds"]
        res["casscf_ladder"] = rungs
        save(res)

        # --- qubit Hamiltonian for the largest converged rung ---
        if best is None:
            res["qubit_hamiltonian"] = {"skipped": "no converged CASSCF rung"}
            res["status"] = "no_converged_casscf"
        else:
            mc, ncas, nelec = best
            qh = {"ncas": ncas, "nelec": nelec, "qubits_jw": 2 * ncas,
                  "mapping": "Jordan-Wigner",
                  "integrals": "exact (non-DF) active integrals in the "
                               "DF-CASSCF-converged orbitals",
                  "note": "counted, NOT diagonalized: the fixed-N sector is "
                          f"~C({2*ncas},{nelec})^2-determinant scale — the "
                          "classical wall this track is about"}
            try:
                h1, eri, ecore = exact_active_integrals(
                    ref_mol, mc.mo_coeff, ncas, tuple(mc.nelecas))
                np.savez_compressed(NPZ, h1=h1, eri=eri, ecore=ecore,
                                    ncas=ncas, nelecas=np.array(mc.nelecas),
                                    mo_active=mc.mo_coeff[:, mc.ncore:mc.ncore + ncas],
                                    e_casscf_df=mc.e_tot)
                qh["integrals_npz"] = os.path.basename(NPZ)
                if backend:
                    t0 = time.time()
                    alarm = int(min(STAGE_TIMEOUT, max(remaining() - 120, 120)))
                    _, nterms = with_timeout(alarm, build_qubit_H,
                                             h1, eri, ecore, ncas, backend)
                    qh["n_pauli_terms"] = int(nterms)
                    qh["backend"] = backend
                    qh["build_seconds"] = round(time.time() - t0, 1)
                    print(f"[qubit] {2*ncas} qubits, {nterms} Pauli terms "
                          f"({backend}, {qh['build_seconds']}s)", flush=True)
                else:
                    qh["error"] = "no qubit backend on this box"
            except StageTimeout:
                qh["error"] = "JW build timed out (integrals npz still saved)"
            except Exception as e:
                qh["error"] = f"{type(e).__name__}: {e}"
            res["qubit_hamiltonian"] = qh
            res["status"] = "ok"
    except Exception as e:
        res["status"] = f"aborted: {type(e).__name__}: {e}"
        print("ABORT:", res["status"], flush=True)
    finally:
        res["meta"]["total_seconds"] = round(elapsed(), 1)
        save(res)
        print(f"\nwrote {RESULTS}  status={res['status']}  "
              f"(total {res['meta']['total_seconds']}s)", flush=True)


if __name__ == "__main__":
    main()
