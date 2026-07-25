#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
calc/fe_ni_s_pentlandite.py — Track N6 OPENER (AWS job): the spin structure of a
Fe-Ni-S cubane of the pentlandite motif.

WHY THIS TRACK. Norilsk Nickel's ore is pentlandite (Fe,Ni)9S8: the value chain
(flotation -> leaching -> refining) is governed by how electrons sit on the mixed
Fe/Ni sulfide surface — the metal oxidation/spin state sets the S 3p covalency,
which sets the sulfur's affinity for xanthate/dithiophosphate collectors and the
kinetics of oxidative leaching. A collector or an oxidant is, at bottom, a redox
partner: to design one we need to know the electronic descriptor of the mineral
surface, and that descriptor is exactly what single-determinant DFT gets wrong on
Fe-S clusters (Stage 15/17 of the FeMoco diary: BS-DFT cannot even fix the SIGN of
the magnetic coupling; ~15 kcal/mol spread across functionals). So we do here for
a Fe2Ni2S4 cubane what calc/fe4s4_casscf_aws.py did for [Fe4S4]: measure the
multireference character honestly (CASSCF NOON / Head-Gordon n_u) and quantify how
badly DFT disagrees with itself (functional spread of the spin gap) — that spread
IS the descriptor of "where DFT lies", the Fe-Ni-S analog of the 885 cm-1
functional-spread argument we used before.

THE CLUSTER. [Fe2Ni2S4(SH)4]2- — the idealized [Fe4S4(SH)4]2- cubane of
calc/femoco_4fe4s_cubane.xyz with two of the four metal sites turned from Fe into
Ni (atoms 2 and 3; the two metals not carrying the reference "up" spins — a
diagonal Fe2/Ni2 arrangement of the M4 tetrahedron). The nuclear charge is the
only thing changed at each swapped site (Fe Z=26 -> Ni Z=28); the M4S4 framework
geometry is FROZEN from the Fe4S4 cubane so energies stay comparable to the
[Fe4S4] reference run. Formal picture at charge -2: 2 Fe(3+) (d5, S=5/2) + 2 Ni(2+)
(d8, S=1) + 4 mu3-S(2-) + 4 terminal SH(-) => metals sum +10, HS seed 2S=14.
Optional light UKS relaxation of the whole cluster is available behind
FENIS_RELAX=1 (default 0 = single point on the frozen Fe4S4 framework, flagged as
such in the JSON — no unlabeled relaxed number).

THE PROGRAM (each stage checkpoints to JSON; a re-launched box RESUMES).
  1. BS-SCF. Converge the ferromagnetic high-spin UKS-PBE seed (2S=14, easy), then
     build broken-symmetry Ms=0 guesses by flipping alpha/beta on the diagonal AO
     blocks of chosen metals, and try 3 BS patterns:
        A  Fe0 up, Fe1 up | Ni2 dn, Ni3 dn   (metal-segregated; the requested one)
        B  Fe0 up, Ni2 up | Fe1 dn, Ni3 dn   (Fe/Ni alternating)
        C  Fe0 up, Ni3 up | Fe1 dn, Ni2 dn   (the other alternating pairing)
     Each candidate is accepted only if it converges AND its Mulliken metal spins
     carry the pattern's expected signs with real polarization (pattern_ok logic,
     as in the [Fe4S4] run); the first accepted pattern is the CASSCF reference.
  2. UNO. Natural orbitals of the total BS density; fractional occupations and the
     count inside [0.02, 1.98] (the magnetic Fe-3d / Ni-3d / bridging-S-3p space).
  3. Ladder CASSCF CAS(10,10) -> CAS(12,12), each rung started from the previous
     rung's relaxed orbitals, with per-rung NOON and n_u.
  4. Spin gap. Two full spin sectors — low S (BS Ms=0) and high S (FM 2S=14) — on
     BOTH UKS-PBE and UKS-PBE0. Four energies => two gaps; the PBE-vs-PBE0 spread
     is the functional-spread descriptor. (BS is spin-contaminated, noted; this is
     the DFT self-disagreement metric, not a spin-projected reference number.)
  5. VERDICT. Does n_u SATURATE with window size, or keep GROWING as it did for
     [Fe4S4] (n_u 9.29 -> 11.84 over the window ladder in that run)? Growth = the
     active space is still too small to contain the correlation = a genuinely
     harder multireference problem.

HONEST SCOPE.
  - Idealized cubane (frozen Fe4S4 framework unless FENIS_RELAX=1), def2-SVP,
    CASSCF only (no dynamic correlation). Qualitative descriptor + DFT-spread, not
    a quantitative surface energy. Every number is computed at run time; SCF /
    CASSCF non-convergence, fallbacks and timeouts are recorded as such, and every
    derived block carries the convergence flags of what it came from.
  - signal.alarm stage caps are best-effort (a long C step returns to Python before
    the exception fires); the shell-level `timeout` in the AWS wrapper is the hard
    stop.

ENV KNOBS (all optional):
  FENIS_STAGE_TIMEOUT  seconds per stage (each SCF / each CASSCF rung), default 5400
  FENIS_TOTAL_BUDGET   total wall-clock budget seconds, default 19800
  FENIS_LADDER         comma list of active-space sizes, default "10,12" (step 2)
  FENIS_MF             spin-gap functionals: pbe | pbe0 | both (default both)
  FENIS_MACRO          CASSCF max macro-iterations per rung, default 60
  FENIS_RELAX          1 = light UKS-PBE geometry relax first (default 0 = frozen)
  FENIS_ONLY           optional: "casscf" (BS+UNO+ladder only) or "gap"
                       (spin-gap matrix only) — used by the split launcher

OUTPUT (all incremental / atomic — nothing lost on a kill; resume-aware):
  calc/fe_ni_s_pentlandite_results.json   rolling results (os.replace atomic writes)
  calc/fe_ni_s_pentlandite_ref.npz        BS reference UNO occ/C (for resume)
  calc/fe_ni_s_pentlandite_cas{N}.npz     each converged rung's orbitals (resume)
  calc/fe_ni_s_pentlandite.xyz            the geometry actually used

Run (heavy; meant for AWS via calc/fe_ni_s_aws.py):
    python3 -u calc/fe_ni_s_pentlandite.py
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
RESULTS = os.path.join(HERE, "fe_ni_s_pentlandite_results.json")
REF_NPZ = os.path.join(HERE, "fe_ni_s_pentlandite_ref.npz")
XYZ_OUT = os.path.join(HERE, "fe_ni_s_pentlandite.xyz")
SRC_XYZ = os.path.join(HERE, "femoco_4fe4s_cubane.xyz")
KCAL = 627.50947406            # Ha -> kcal/mol
CM = 219474.6313705            # Ha -> cm^-1

BASIS = "def2-svp"
CHARGE = -2
NI_ATOMS = (2, 3)              # metal sites turned Fe -> Ni (diagonal M2/M2 split)
FE_ATOMS = (0, 1)
METAL_ATOMS = (0, 1, 2, 3)
BRIDGE_S_ATOMS = (4, 5, 6, 7)  # mu3-S core sulfurs (atom order of the cubane xyz)
# HS seed: 2 Fe3+ (d5, 5 unpaired) + 2 Ni2+ (d8, 2 unpaired), all aligned -> 2S=14
SPIN_HS = 14

STAGE_TIMEOUT = int(os.environ.get("FENIS_STAGE_TIMEOUT", "5400"))
TOTAL_BUDGET = int(os.environ.get("FENIS_TOTAL_BUDGET", "19800"))
LADDER = [int(x) for x in os.environ.get("FENIS_LADDER", "10,12").split(",")]
MF_GAP = os.environ.get("FENIS_MF", "both").lower()
MACRO = int(os.environ.get("FENIS_MACRO", "60"))
RELAX = os.environ.get("FENIS_RELAX", "0") == "1"
ONLY = os.environ.get("FENIS_ONLY", "").lower()      # "", "casscf", "gap"

# BS spin-flip patterns: (name, flip_atoms, expected_signs on metals 0..3)
BS_PATTERNS = [
    ("A: Fe0 up, Fe1 up | Ni2 dn, Ni3 dn (metal-segregated)", (2, 3),
     (+1, +1, -1, -1)),
    ("B: Fe0 up, Ni2 up | Fe1 dn, Ni3 dn (Fe/Ni alternating)", (1, 3),
     (+1, -1, +1, -1)),
    ("C: Fe0 up, Ni3 up | Fe1 dn, Ni2 dn (Fe/Ni alternating')", (1, 2),
     (+1, -1, -1, +1)),
]

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


def load_prev():
    """Resume support: reload a prior (possibly partial) results JSON."""
    if os.path.exists(RESULTS):
        try:
            with open(RESULTS) as fh:
                return json.load(fh)
        except Exception:
            return None
    return None


class StageTimeout(Exception):
    pass


def _alarm_handler(signum, frame):
    raise StageTimeout()


def with_timeout(seconds, fn, *args, **kwargs):
    """Best-effort per-stage wall-clock cap via SIGALRM (main thread only).
    The shell-level `timeout` in the AWS wrapper is the hard backstop."""
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
def read_cubane_xyz(path):
    """Read the [Fe4S4(SH)4]2- cubane geometry, then swap the NI_ATOMS metal
    sites Fe -> Ni (charge on the nucleus is the only change; framework frozen).
    Returns [(element, np.array(xyz)), ...] in the source atom order."""
    with open(path) as fh:
        lines = fh.read().splitlines()
    n = int(lines[0].split()[0])
    atoms = []
    for i in range(2, 2 + n):
        parts = lines[i].split()
        el = parts[0]
        xyz = np.array([float(parts[1]), float(parts[2]), float(parts[3])])
        atoms.append([el, xyz])
    # sanity: first four atoms must be the metals we are about to (re)assign
    if [a[0] for a in atoms[:4]] != ["Fe", "Fe", "Fe", "Fe"]:
        raise ValueError("unexpected metal order in cubane xyz; "
                         f"got {[a[0] for a in atoms[:4]]}")
    for a in NI_ATOMS:
        atoms[a][0] = "Ni"
    return [(el, xyz) for el, xyz in atoms]


def write_xyz(atoms, path, comment):
    with open(path, "w") as fh:
        fh.write(f"{len(atoms)}\n{comment}\n")
        for e, p in atoms:
            fh.write(f"{e} {p[0]:.6f} {p[1]:.6f} {p[2]:.6f}\n")


def relax_geometry(atoms, maxmem):
    """Optional light UKS-PBE relaxation of the whole cluster (FENIS_RELAX=1).
    Uses the BS Ms=0 reference; best-effort — falls back to the frozen geometry
    (recorded) if geometric/convergence is unavailable."""
    try:
        from pyscf.geomopt.geometric_solver import optimize
    except ImportError:
        return atoms, {"relaxed": False, "reason": "geometric not installed"}
    mol = gto.M(atom=[(e, tuple(p)) for e, p in atoms], basis=BASIS,
                charge=CHARGE, spin=0, verbose=0, max_memory=maxmem)
    mf = _new_mf(mol, "pbe")
    dm0 = flip_dm(_seed_hs(atoms, maxmem).make_rdm1(), mol, BS_PATTERNS[0][1])
    mf.kernel(dm0=dm0)
    try:
        mol_eq = with_timeout(min(STAGE_TIMEOUT, remaining()),
                              optimize, mf, maxsteps=20)
        new = [(mol_eq.atom_symbol(i), mol_eq.atom_coord(i, unit="Angstrom"))
               for i in range(mol_eq.natm)]
        return new, {"relaxed": True, "maxsteps": 20}
    except Exception as e:
        return atoms, {"relaxed": False, "reason": f"{type(e).__name__}: {e}"}


def _seed_hs(atoms, maxmem):
    mol = gto.M(atom=[(e, tuple(p)) for e, p in atoms], basis=BASIS,
                charge=CHARGE, spin=SPIN_HS, verbose=0, max_memory=maxmem)
    mf = _new_mf(mol, "pbe")
    mf.kernel()
    return mf


# ----------------------------------------------------------------- mean field
def _new_mf(mol, xc):
    """UKS with the requested functional (density-fitted)."""
    mf = dft.UKS(mol).density_fit()
    mf.xc = xc
    return mf


def _converge(mf, dm0=None, level_shift=0.4, max_cycle=60, tag=""):
    """DIIS with level shift; on stall, a second-order (Newton/SOSCF) restart —
    the robust recipe for hard Fe-S/Ni-S SCF cases."""
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
    """Swap alpha/beta density on the diagonal AO blocks of `flip_atoms`
    (standard BS spin-flip construction; off-diagonal blocks untouched)."""
    dma, dmb = dm[0].copy(), dm[1].copy()
    aoslice = mol.aoslice_by_atom()
    for a in flip_atoms:
        sl = slice(int(aoslice[a][2]), int(aoslice[a][3]))
        blk_a = dma[sl, sl].copy()
        dma[sl, sl] = dmb[sl, sl]
        dmb[sl, sl] = blk_a
    return np.array([dma, dmb])


def metal_spin_pops(mf):
    """Mulliken spin population on each metal (Fe0,Fe1,Ni2,Ni3) — BS sanity."""
    dma, dmb = mf.make_rdm1()
    m = (dma - dmb) @ mf.get_ovlp()
    aoslice = mf.mol.aoslice_by_atom()
    out = []
    for a in METAL_ATOMS:
        sl = slice(int(aoslice[a][2]), int(aoslice[a][3]))
        out.append(round(float(np.trace(m[sl, sl])), 3))
    return out


def converge_reference(mol_hs, mol_bs, info):
    """HS seed -> try BS patterns in order -> first accepted is the reference.
    Returns (mf, mode, pattern) where mode is 'BS' or 'HS-fallback'."""
    print(f"[scf] HS seed (UKS-PBE, 2S={SPIN_HS})", flush=True)
    t0 = time.time()
    mf_hs = with_timeout(min(STAGE_TIMEOUT, remaining()), _converge,
                         _new_mf(mol_hs, "pbe"), tag="HS")
    info["hs"] = {"e": float(mf_hs.e_tot), "converged": bool(mf_hs.converged),
                  "seconds": round(time.time() - t0, 1)}
    print(f"[scf] HS: E={mf_hs.e_tot:.5f} conv={mf_hs.converged} "
          f"({info['hs']['seconds']}s)", flush=True)

    hs_dm = mf_hs.make_rdm1()
    info["bs_attempts"] = []
    for name, flip_atoms, signs in BS_PATTERNS:
        if remaining() < 900:
            info["bs_attempts"].append({"pattern": name, "skipped": "time budget"})
            break
        print(f"[scf] BS pattern {name}: flip {flip_atoms}, Ms=0", flush=True)
        dm0 = flip_dm(hs_dm, mol_hs, flip_atoms)
        t0 = time.time()
        mf_bs = with_timeout(min(STAGE_TIMEOUT, remaining()), _converge,
                             _new_mf(mol_bs, "pbe"), dm0=dm0, level_shift=0.3,
                             tag="BS")
        ss, _ = mf_bs.spin_square()
        pops = metal_spin_pops(mf_bs)
        signs_ok = all((p > 0) == (s > 0) for p, s in zip(pops, signs))
        polarized = all(abs(p) > 1.0 for p in pops)   # Ni2+ moment ~1.7, lower bar
        rec = {"pattern": name, "flip_atoms": list(flip_atoms),
               "expected_signs": list(signs),
               "e": float(mf_bs.e_tot), "converged": bool(mf_bs.converged),
               "s_squared": round(float(ss), 3),
               "metal_mulliken_spin": pops,
               "seconds": round(time.time() - t0, 1),
               "pattern_ok": bool(mf_bs.converged and signs_ok and polarized)}
        info["bs_attempts"].append(rec)
        print(f"[scf] BS {name[:1]}: E={mf_bs.e_tot:.5f} conv={mf_bs.converged} "
              f"<S^2>={ss:.2f} spins={pops} ok={rec['pattern_ok']} "
              f"({rec['seconds']}s)", flush=True)
        if rec["pattern_ok"]:
            info["reference_pattern"] = name
            return mf_bs, "BS", name

    print("[scf][WARN] no BS pattern accepted (unconverged/collapsed) — falling "
          "back to the HS reference (poorer CAS start, recorded as such)",
          flush=True)
    info["reference_pattern"] = "HS-fallback"
    return mf_hs, "HS-fallback", "HS-fallback"


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
    """AO indices of Fe 3d + Ni 3d + bridging S 3p (atoms 4-7): the magnetic
    target set of the mixed-metal cubane."""
    idx = []
    for i, lab in enumerate(mol.ao_labels()):
        parts = lab.split()
        atom, el, shell = int(parts[0]), parts[1], parts[2]
        if el in ("Fe", "Ni") and shell.startswith("3d"):
            idx.append(i)
        elif el == "S" and atom in BRIDGE_S_ATOMS and shell.startswith("3p"):
            idx.append(i)
    return np.array(idx, dtype=int)


def select_active(mol, occ, C, ncas):
    """Choose the `ncas` UNOs closest to single occupation (|n-1| minimal) — the
    magnetic space of the BS solution. Electron count = rounded sum of their
    occupations (forced so ncore is integral and na=nb)."""
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
            "metal3d_bridgeS3p_frac": [round(float(x), 3) for x in frac],
            "min_core_occ": round(float(occ[core[-1]]), 3) if core else None,
            "max_virt_occ": round(float(occ[virt[0]]), 3) if virt else None}
    return mo, nelec, diag


def expand_active(mol, mc):
    """Grow a converged CAS by 2: promote the core and virtual orbital with the
    highest metal-3d + bridging-S-3p Loewdin character into the window, keeping
    the rung's relaxed orbitals."""
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
        info["active_metal3d_bridgeS3p_mean"] = round(float(frac.mean()), 3)
    except StageTimeout:
        info["converged"] = False
        info["timed_out_after_s"] = int(alarm_s)
        print(f"  [rung {ncas}o] TIMED OUT after {alarm_s}s", flush=True)
    info["seconds"] = round(time.time() - t0, 1)
    return mc, info


# ------------------------------------------------------------------- spin gap
def _sector_energy(atoms, maxmem, xc, spin, dm0=None, bs_flip=None, tag=""):
    """Converge one UKS spin sector; if bs_flip is given, build the BS Ms=0 guess
    by flipping those metal AO blocks of a converged HS density."""
    mol = gto.M(atom=[(e, tuple(p)) for e, p in atoms], basis=BASIS,
                charge=CHARGE, spin=spin, verbose=0, max_memory=maxmem)
    mf = _new_mf(mol, xc)
    lvl = 0.3 if bs_flip is not None else 0.4
    mf = _converge(mf, dm0=dm0, level_shift=lvl, tag=tag)
    return mf


def spin_gap(atoms, maxmem, functionals, best_flip, info):
    """Low-S (BS Ms=0) vs high-S (FM 2S=14) energy on each functional; the
    PBE-vs-PBE0 spread of the gap is the DFT self-disagreement descriptor."""
    gaps = {}
    for xc in functionals:
        if remaining() < 1200:
            info[xc] = {"skipped": "time budget"}
            continue
        try:
            print(f"[gap] {xc.upper()} high-S (2S={SPIN_HS})", flush=True)
            t0 = time.time()
            mf_hi = with_timeout(min(STAGE_TIMEOUT, remaining()), _sector_energy,
                                 atoms, maxmem, xc, SPIN_HS, tag=f"{xc}-HS")
            print(f"[gap] {xc.upper()} low-S (BS Ms=0)", flush=True)
            dm0 = flip_dm(mf_hi.make_rdm1(),
                          mf_hi.mol, best_flip) if best_flip else None
            mf_lo = with_timeout(min(STAGE_TIMEOUT, remaining()), _sector_energy,
                                 atoms, maxmem, xc, 0, dm0=dm0, bs_flip=best_flip,
                                 tag=f"{xc}-BS")
            ss_lo, _ = mf_lo.spin_square()
            e_hi, e_lo = float(mf_hi.e_tot), float(mf_lo.e_tot)
            gap = e_hi - e_lo
            gaps[xc] = gap
            info[xc] = {
                "e_high_spin": e_hi, "high_spin_2S": SPIN_HS,
                "high_spin_converged": bool(mf_hi.converged),
                "e_low_spin_bs": e_lo, "low_spin_converged": bool(mf_lo.converged),
                "low_spin_s_squared": round(float(ss_lo), 3),
                "gap_high_minus_low_kcal": round(gap * KCAL, 2),
                "gap_high_minus_low_cm": round(gap * CM, 1),
                "seconds": round(time.time() - t0, 1),
                "note": "BS low-S is spin-contaminated; this is the DFT "
                        "self-disagreement metric, not a spin-projected gap"}
            print(f"[gap] {xc.upper()}: gap={gap*KCAL:.2f} kcal/mol "
                  f"({gap*CM:.1f} cm-1)", flush=True)
        except StageTimeout:
            info[xc] = {"timed_out": True}
            print(f"[gap] {xc.upper()} TIMED OUT", flush=True)
        except Exception as e:
            info[xc] = {"error": f"{type(e).__name__}: {e}"}
    if "pbe" in gaps and "pbe0" in gaps:
        spread = gaps["pbe0"] - gaps["pbe"]
        info["functional_spread_kcal"] = round(spread * KCAL, 2)
        info["functional_spread_cm"] = round(spread * CM, 1)
        info["descriptor"] = (
            "PBE-vs-PBE0 spin-gap spread = 'where DFT lies' for the Fe-Ni-S "
            "motif (Fe-Ni-S analog of the 885 cm-1 functional-spread argument); "
            "large spread => the flotation/leaching redox descriptor cannot be "
            "trusted from a single functional, motivating the CASSCF track")
    return info


# ------------------------------------------------------------------------ main
def make_verdict(rungs):
    """n_u trend across the ladder: saturating vs growing (as [Fe4S4] did,
    9.29 -> 11.84 over its window ladder)."""
    nus = [(r["ncas"], r["n_u"]) for r in rungs
           if r.get("converged") and "n_u" in r]
    if len(nus) < 2:
        return {"n_u_by_window": nus,
                "verdict": "inconclusive: fewer than 2 converged rungs"}
    grew = nus[-1][1] - nus[0][1]
    growing = grew > 0.3
    return {"n_u_by_window": nus,
            "delta_n_u": round(grew, 3),
            "fe4s4_reference_delta": "9.29 -> 11.84 (grew ~2.5 over its ladder)",
            "verdict": ("n_u GROWS with window (like [Fe4S4]): active space still "
                        "too small to hold the correlation — genuinely hard "
                        "multireference Fe-Ni-S problem"
                        if growing else
                        "n_u SATURATES with window: the chosen active space "
                        "captures the multireference character (unlike [Fe4S4])")}


def main():
    prev = load_prev()
    resume = prev is not None and prev.get("meta", {}).get("cluster")
    res = prev if resume else {"status": "running"}
    res.setdefault("meta", {})
    res["meta"].update({
        "cluster": "[Fe2Ni2S4(SH)4]2- cubane (pentlandite motif)",
        "purpose": "открытие трека N6: спиновая структура Fe-Ni-S мотива — "
                   "дескриптор для флотореагентов/выщелачивания сульфидных руд "
                   "(пентландит (Fe,Ni)9S8, Норникель)",
        "basis": BASIS, "charge": CHARGE,
        "metal_layout": {"Fe": list(FE_ATOMS), "Ni": list(NI_ATOMS)},
        "formal_state": "2 Fe(3+) d5 (S=5/2) + 2 Ni(2+) d8 (S=1); metals sum +10",
        "hs_seed_2S": SPIN_HS,
        "geometry_source": "femoco_4fe4s_cubane.xyz (Fe4S4 framework), Fe->Ni at "
                           f"atoms {list(NI_ATOMS)}",
        "geometry_note": ("light UKS-PBE relax (FENIS_RELAX=1)" if RELAX else
                          "SINGLE POINT on the FROZEN Fe4S4 framework geometry "
                          "(no relaxation; energies comparable to the [Fe4S4] run)"),
        "ladder": LADDER, "macro_max": MACRO,
        "spin_gap_functionals": MF_GAP, "only": ONLY or "full",
        "stage_timeout_s": STAGE_TIMEOUT, "total_budget_s": TOTAL_BUDGET,
        "note": "BS-UNO start for the CASSCF ladder (the [Fe4S4] recipe). CASSCF "
                "only, no dynamic correlation. DF-accelerated. Spin gap on PBE and "
                "PBE0; the functional spread is the DFT-lies descriptor. Resume-"
                "aware: a re-launched box continues from the JSON + npz checkpoints.",
        "resumed": bool(resume)})
    save(res)
    try:
        if any(b - a != 2 for a, b in zip(LADDER, LADDER[1:])):
            raise ValueError(f"FENIS_LADDER must grow in steps of 2, got {LADDER}")

        atoms = read_cubane_xyz(SRC_XYZ)
        maxmem = int(os.environ.get("PYSCF_MAX_MEMORY", "16000"))
        if RELAX and "relax" not in res:
            atoms, relax_info = relax_geometry(atoms, maxmem)
            res["relax"] = relax_info
            save(res)
        write_xyz(atoms, XYZ_OUT, "[Fe2Ni2S4(SH)4]2- pentlandite-motif cubane "
                                  "(Track N6 opener, AWS)")
        mol_kw = dict(atom=[(e, tuple(p)) for e, p in atoms], basis=BASIS,
                      charge=CHARGE, verbose=0, max_memory=maxmem)
        mol_hs = gto.M(spin=SPIN_HS, **mol_kw)
        mol_bs = gto.M(spin=0, **mol_kw)
        res["meta"]["nao"] = int(mol_bs.nao)
        res["meta"]["nelec_total"] = int(mol_bs.nelectron)
        print(f"[mol] nao={mol_bs.nao} nelec={mol_bs.nelectron} "
              f"maxmem={maxmem}MB relax={RELAX}", flush=True)
        save(res)

        do_casscf = ONLY in ("", "casscf")
        do_gap = ONLY in ("", "gap")

        # ---- reference for BOTH the CASSCF ladder and the low-S gap flip ------
        best_flip = None
        occ = C = None
        mode = res.get("scf", {}).get("reference_used")
        if resume and os.path.exists(REF_NPZ) and mode:
            z = np.load(REF_NPZ)
            occ, C = z["occ"], z["C"]
            best_flip = tuple(int(x) for x in z["flip"]) if z["flip"].size else None
            print(f"[resume] reference reloaded from {os.path.basename(REF_NPZ)} "
                  f"(mode={mode})", flush=True)
        else:
            scf_info = {}
            mf, mode, pattern = converge_reference(mol_hs, mol_bs, scf_info)
            scf_info["reference_used"] = mode
            res["scf"] = scf_info
            save(res)
            occ, C = unos(mf)
            best_flip = (dict((n, f) for n, f, _ in BS_PATTERNS).get(pattern)
                         if mode == "BS" else None)
            np.savez_compressed(REF_NPZ, occ=occ, C=C,
                                flip=np.array(best_flip if best_flip else [],
                                              dtype=int))
        ref_mol = mol_bs if mode == "BS" else mol_hs

        # ---- UNO spectrum -----------------------------------------------------
        if "uno" not in res:
            frac_mask = (occ > 0.02) & (occ < 1.98)
            res["uno"] = {"n_fractional_0.02_1.98": int(frac_mask.sum()),
                          "fractional_occs": [round(float(x), 3)
                                              for x in occ[frac_mask][:40]],
                          "from_reference": mode}
            print(f"[uno] {int(frac_mask.sum())} fractionally-occupied naturals",
                  flush=True)
            save(res)

        # ---- ladder CASSCF ----------------------------------------------------
        if do_casscf:
            rungs = res.get("casscf_ladder", [])
            done = {r["ncas"] for r in rungs if r.get("converged")}
            best = None
            last_secs = None
            for i, ncas in enumerate(LADDER):
                cas_npz = os.path.join(HERE, f"fe_ni_s_pentlandite_cas{ncas}.npz")
                if ncas in done and os.path.exists(cas_npz):
                    z = np.load(cas_npz)
                    mo_prev = z["mo"]
                    nelec_prev = int(z["nelec"])
                    # rebuild a light mc holder for expand_active on resume
                    mf_r = scf.RHF(ref_mol).density_fit()
                    mc = mcscf.CASSCF(mf_r, ncas, (nelec_prev // 2,) * 2)
                    mc.mo_coeff = mo_prev
                    best = (mc, ncas, nelec_prev)
                    print(f"[resume] rung {ncas}o reloaded", flush=True)
                    continue
                need = max(900, 1.5 * last_secs) if last_secs else 900
                if remaining() < need:
                    rungs.append({"ncas": ncas, "skipped": f"time budget: "
                                  f"{int(remaining())}s left < {int(need)}s"})
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
                rungs = [r for r in rungs if r.get("ncas") != ncas] + [info]
                rungs.sort(key=lambda r: r["ncas"])
                res["casscf_ladder"] = rungs
                save(res)
                print(f"[ladder] {ncas}o: conv={info.get('converged')} "
                      f"E={info.get('e_casscf_df')} n_u={info.get('n_u')} "
                      f"({info['seconds']}s)", flush=True)
                if not info.get("converged"):
                    break
                np.savez_compressed(cas_npz, mo=mc.mo_coeff, nelec=nelec, ncas=ncas)
                best = (mc, ncas, nelec)
                last_secs = info["seconds"]
            res["casscf_ladder"] = rungs
            res["verdict"] = make_verdict(rungs)
            save(res)

        # ---- spin-gap matrix (PBE and PBE0) -----------------------------------
        if do_gap and "spin_gap" not in res:
            funcs = {"pbe": ["pbe"], "pbe0": ["pbe0"],
                     "both": ["pbe", "pbe0"]}.get(MF_GAP, ["pbe", "pbe0"])
            gap_flip = best_flip if best_flip else BS_PATTERNS[0][1]
            gap_info = {"low_S_flip_atoms": list(gap_flip)}
            spin_gap(atoms, maxmem, funcs, gap_flip, gap_info)
            res["spin_gap"] = gap_info
            save(res)

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
