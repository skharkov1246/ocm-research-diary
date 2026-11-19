#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
calc/flotation_adsorption.py — Track N9 GATE (AWS job): adsorption of flotation
reagents on a PENTLANDITE vs a PYRRHOTITE motif.

WHY THIS JOB (APPLIED_WIDE_SCAN.md, section N9 — top of the second echelon).
The main beneficiation pain of Talnakh is separating pentlandite (Fe,Ni)9S8 from
pyrrhotite Fe7S8: the pyrrhotite drags sulfur and iron into the smelter charge,
and that sulfur is exactly the SO2 the Sulfur Programme then catches downstream
for 250+ bn RUB. Every percent of flotation selectivity works BEFORE the stack
and is cheaper than any capture. Selectivity of a collector/depressant is set by
how it adsorbs on metal sites whose spin/covalency differ between the two
sulfides — and that electronic gate we already measured
(calc/fe_ni_s_pentlandite.py: BS pattern A metal-segregated, Fe +2.03/+2.03,
Ni -1.27/-1.26; spin gap PBE -18.8 vs PBE0 -30.7 kcal/mol => ~12 kcal/mol of
pure functional disagreement). This job is the next question in the chain, and
the cheapest honest one: DOES A REAGENT ACTUALLY PREFER ONE MOTIF OVER THE
OTHER, AND IN WHICH DIRECTION?

THE MODEL — two substrates, one open metal site each.
  Both substrates are built from the SAME idealized cubane framework
  (calc/femoco_4fe4s_cubane.xyz, the geometry the whole Fe-S diary uses), so the
  two motifs differ ONLY in nuclear charge at two metal sites — the cleanest
  possible A/B comparison:
    pentlandite motif  [Fe2Ni2S4(SH)4]2-  (Ni at atoms 2,3 — the metal-segregated
                       BS pattern A that already converged in the N6 opener)
    pyrrhotite  motif  [Fe4S4(SH)4]2-     (all four metals Fe; the same cluster
                       used by calc/fe4s4_casscf_aws.py and calc/fes_cogef.py)
  A mineral SURFACE metal site is coordinatively unsaturated, so each substrate
  is the cubane with ONE terminal SH- removed from the target metal:
    pent_fe   [Fe2Ni2S4(SH)3]1-  vacancy on Fe0   (2S_HS=14, BS flip Ni2,Ni3)
    pent_ni   [Fe2Ni2S4(SH)3]1-  vacancy on Ni2   (2S_HS=14, BS flip Ni2,Ni3)
    pyrr_fe   [Fe4S4(SH)3]1-     vacancy on Fe0   (2S_HS=18, BS flip M2,M3)
  Removing an anionic SH- leaves the formal metal oxidation states untouched
  (pentlandite: 2 Fe3+ d5 + 2 Ni2+ d8, metals sum +10; pyrrhotite: 2 Fe3+ + 2
  Fe2+, metals sum +10) — only the cluster charge goes -2 -> -1. TWO pentlandite
  sites are computed on purpose: the selectivity may sit in the NATURE of the
  site (Fe vs Ni), not in the mineral as a whole, and that distinction is
  actionable for reagent design.

THE THREE ADSORBATES (small, model — this is a screen, not a passport).
  xanthate  EtOCS2-  ethyl xanthate, the workhorse collector. Monodentate
            through one thiolate S. EXPECTED: binds STRONGER on pentlandite.
  oh        OH-      the simplest depressant / alkaline-pulp proxy.
  en        NH2CH2CH2NH2  ethylenediamine — a small stand-in for DETA
            (diethylenetriamine), the industrial pyrrhotite depressant used in
            Sudbury. Monodentate through one N. EXPECTED: binds STRONGER on
            pyrrhotite.

WHAT IS COMPUTED.
  dE_ads = E(complex) - E(substrate) - E(adsorbate),  UKS-PBE/def2-SVP + density
  fitting, BROKEN-SYMMETRY reference on every cluster (HS seed -> spin-flip of
  the same metal AO blocks -> pattern check on Mulliken metal spins, the exact
  recipe of calc/fe_ni_s_pentlandite.py). The cluster framework is FROZEN; what
  is relaxed is the M-ligand bond:
    FLOT_MODE=scan (default) — a rigid scan of the M-donor distance
            (FLOT_SCAN_N points, +-FLOT_SCAN_DR A around a chemically sensible
            r0); dE_ads is taken at the LOWEST COMPUTED POINT (a real number),
            with a parabola-refined minimum reported separately as an
            interpolation.
    FLOT_MODE=opt  — geomeTRIC constrained optimisation with all substrate atoms
            frozen (adsorbate + M-L bond free). Falls back to `scan` and says so
            if geomeTRIC is missing or the step times out.
    FLOT_MODE=rigid— one single point at r0, flagged rigid.
  The free adsorbates are PBE-relaxed once (FLOT_ADS_RELAX=1, cheap).

THE MAIN DESCRIPTOR — two double differences, both of which must be NEGATIVE
for the separation to work:
  dd_E(collector)  = dE_ads(xanthate, pentlandite) - dE_ads(xanthate, pyrrhotite)
  dd_E(depressant) = dE_ads(en,       pyrrhotite)  - dE_ads(en,       pentlandite)
  dd_E(oh)         = dE_ads(OH-,      pyrrhotite)  - dE_ads(OH-,      pentlandite)
                     (baseline: how much selectivity alkalinity alone can buy)
WHY A DOUBLE DIFFERENCE AND NOT A BINDING ENERGY. In dd_E the free-adsorbate
energy CANCELS EXACTLY (same molecule on both sides), which removes the single
worst error of this level of theory — def2-SVP without diffuse functions
describes a free anion (xanthate-, OH-) badly. The two substrates are clusters
of identical size and basis, so basis-set superposition error largely cancels
too. What survives is the difference in metal-site chemistry, which is the thing
we are actually asking about.
THEN THE HYBRID CHECK. Lesson of the spin-transistor / N6 gate: PBE magnitudes
on Fe-S are not trustworthy (12 kcal/mol of functional spread on the spin gap
alone), so the SIGN of the best pair is re-checked at PBE0//PBE (single points
on the PBE geometry, PBE density as the start). If PBE and PBE0 disagree on the
SIGN, the descriptor is declared not resolvable at this level — that is a real
result, not a failure to hide.

HONEST SCOPE (read before quoting any number).
  - A 15-atom cluster is NOT a mineral lattice: no Madelung field, no surface
    reconstruction, no Fe/Ni site disorder, no S-vacancy chemistry of real
    pyrrhotite Fe7S8.
  - Gas phase: no water, no electric double layer, no pH, no pulp potential, no
    competing Ca2+/Cu2+ activation — all of which are first-order in real
    flotation.
  - Electronic energies only: no ZPE, no thermal or entropic terms, no solvation
    -> this is NOT an adsorption isotherm and NOT a free energy of adsorption.
  - Monodentate binding only (one vacancy per cluster). Real xanthate can chelate
    bidentate; DETA is a tridentate chelate and `en` is only its two-nitrogen
    proxy — chelation is precisely where a depressant gets much of its strength,
    so absolute numbers here UNDERSTATE the amine.
  - BS solutions are spin-contaminated; no spin projection is applied.
  - No counterpoise/BSSE correction (it largely cancels in dd_E, see above).
  THIS IS A SIGN SCREEN OF A DESCRIPTOR. It can say "this reagent class prefers
  this motif, and by which site" — it cannot say "adsorb X g/t at pH 9".

ENV KNOBS (all optional):
  FLOT_STAGE_TIMEOUT seconds per SCF / per relax step, default 5400
  FLOT_TOTAL_BUDGET  total wall-clock budget seconds, default 32400
  FLOT_ONLY          comma list restricting the work; accepts pair ids
                     ("xanthate@pent_fe"), adsorbate ids ("xanthate,oh"),
                     substrate ids ("pyrr_fe") or roles ("collector")
  FLOT_MF            pbe (default: PBE everywhere + PBE0 push on the best pair)
                     | pbe0 (only the PBE0 push stage, resuming a PBE run)
                     | both (PBE0 re-check on every computed pair)
  FLOT_PBE0_SET      which descriptor the PBE0 push covers:
                     collector | depressant | both (default both)
  FLOT_MODE          scan (default) | rigid | opt
  FLOT_SCAN_N        scan points, default 3       FLOT_SCAN_DR  half-range A, 0.2
  FLOT_ADS_RELAX     1 = PBE-relax the free adsorbates (default 1)
  PYSCF_MAX_MEMORY   per-process memory hint (MB)

OUTPUT (incremental / atomic — nothing lost on a kill; resume-aware):
  calc/flotation_adsorption_results.json
  calc/flotation_adsorption_geoms.xyz        best geometry of every complex
  calc/flotation_sub_<id>.npz                substrate BS density (resume)
  calc/flotation_dm_<pair>.npz               complex BS density (PBE0 push)

Run (heavy; meant for AWS via calc/flotation_aws.py):
    python3 -u calc/flotation_adsorption.py
Geometry self-check (no pyscf needed, prints the job table and clash distances):
    python3 calc/flotation_adsorption.py --geoms
"""
import os
import sys
import json
import time
import signal
import tempfile

try:
    import numpy as np
except ImportError:
    sys.exit("[deps] numpy is required:  python3 -m pip install numpy")

HAVE_PYSCF = True
PYSCF_ERR = ""
try:
    from pyscf import gto, dft
except ImportError as e:            # allow --geoms on a machine without pyscf
    HAVE_PYSCF = False
    PYSCF_ERR = str(e)

HERE = os.path.dirname(os.path.abspath(__file__))
RESULTS = os.path.join(HERE, "flotation_adsorption_results.json")
GEOMS = os.path.join(HERE, "flotation_adsorption_geoms.xyz")
SRC_XYZ = os.path.join(HERE, "femoco_4fe4s_cubane.xyz")
KCAL = 627.50947406                 # Ha -> kcal/mol

BASIS = "def2-svp"
METAL_ATOMS = (0, 1, 2, 3)          # the M4 tetrahedron of the cubane
# in the source cubane the terminal thiolate of metal i is S(8+i) with H(12+i)
TERM_S = {i: 8 + i for i in METAL_ATOMS}
TERM_H = {i: 12 + i for i in METAL_ATOMS}

STAGE_TIMEOUT = int(os.environ.get("FLOT_STAGE_TIMEOUT", "5400"))
TOTAL_BUDGET = int(os.environ.get("FLOT_TOTAL_BUDGET", "32400"))
ONLY = [x.strip() for x in os.environ.get("FLOT_ONLY", "").split(",") if x.strip()]
MF_MODE = os.environ.get("FLOT_MF", "pbe").lower()
PBE0_SET = os.environ.get("FLOT_PBE0_SET", "both").lower()
MODE = os.environ.get("FLOT_MODE", "scan").lower()
SCAN_N = max(1, int(os.environ.get("FLOT_SCAN_N", "3")))
SCAN_DR = float(os.environ.get("FLOT_SCAN_DR", "0.2"))
ADS_RELAX = os.environ.get("FLOT_ADS_RELAX", "1") == "1"

T0 = time.time()


# =============================================================== substrates ===
# Both motifs share the cubane framework; `ni_atoms` is the ONLY chemical change
# (nuclear charge Fe Z=26 -> Ni Z=28), `site` is the metal whose terminal SH- is
# removed to create the open surface site the reagent binds to.
SUBSTRATES = {
    "pent_fe": {
        "label": "pentlandite motif [Fe2Ni2S4(SH)3]1-, vacancy on Fe0",
        "motif": "pentlandite", "ni_atoms": (2, 3), "site": 0, "site_el": "Fe",
        "charge": -1, "spin_hs": 14, "bs_flip": (2, 3),
        "bs_signs": (+1, +1, -1, -1),
        "formal": "2 Fe(3+) d5 + 2 Ni(2+) d8; metals sum +10"},
    "pent_ni": {
        "label": "pentlandite motif [Fe2Ni2S4(SH)3]1-, vacancy on Ni2",
        "motif": "pentlandite", "ni_atoms": (2, 3), "site": 2, "site_el": "Ni",
        "charge": -1, "spin_hs": 14, "bs_flip": (2, 3),
        "bs_signs": (+1, +1, -1, -1),
        "formal": "2 Fe(3+) d5 + 2 Ni(2+) d8; metals sum +10"},
    "pyrr_fe": {
        "label": "pyrrhotite motif [Fe4S4(SH)3]1-, vacancy on Fe0",
        "motif": "pyrrhotite", "ni_atoms": (), "site": 0, "site_el": "Fe",
        "charge": -1, "spin_hs": 18, "bs_flip": (2, 3),
        "bs_signs": (+1, +1, -1, -1),
        "formal": "2 Fe(3+) d5 + 2 Fe(2+) d6; metals sum +10"},
}

# ================================================================ adsorbates ===
# Built-in starting geometries (A). `donor` = the atom that touches the metal;
# `lp` = the direction, in the molecule's own frame, in which the METAL should
# sit relative to the donor (the lone-pair axis). `r0` = starting M-donor
# distance per metal element.
XANTHATE = [
    ("S", (-1.508, -0.785, 0.000)),   # 0 donor thiolate S
    ("C", (0.000, 0.000, 0.000)),     # 1 CS2 carbon
    ("S", (1.508, -0.785, 0.000)),    # 2 the other S
    ("O", (0.000, 1.340, 0.000)),     # 3 ester O
    ("C", (1.272, 2.015, 0.000)),     # 4 CH2
    ("C", (1.049, 3.520, 0.000)),     # 5 CH3
    ("H", (1.834, 1.731, 0.890)),
    ("H", (1.834, 1.731, -0.890)),
    ("H", (0.996, 3.880, 1.028)),
    ("H", (1.876, 4.011, -0.515)),
    ("H", (0.115, 3.749, -0.515)),
]
HYDROXIDE = [
    ("O", (0.000, 0.000, 0.000)),     # 0 donor
    ("H", (0.970, 0.000, 0.000)),
]
ETHYLENEDIAMINE = [
    ("N", (1.261, 1.386, 0.000)),     # 0 donor N
    ("C", (0.770, 0.000, 0.000)),
    ("C", (-0.770, 0.000, 0.000)),
    ("N", (-1.261, -1.386, 0.000)),   # distal N (free — monodentate model)
    ("H", (1.133, -0.514, 0.891)),
    ("H", (1.133, -0.514, -0.891)),
    ("H", (-1.133, 0.514, 0.891)),
    ("H", (-1.133, 0.514, -0.891)),
    ("H", (0.920, 1.866, 0.833)),
    ("H", (0.920, 1.866, -0.833)),
    ("H", (-0.920, -1.866, 0.833)),
    ("H", (-0.920, -1.866, -0.833)),
]

ADSORBATES = {
    "xanthate": {
        "label": "ethyl xanthate EtOCS2- (collector)",
        "role": "collector", "atoms": XANTHATE, "donor": 0,
        "lp": (-0.676, 0.737, 0.000),          # ~105 deg M-S-C
        "charge": -1, "spin": 0,
        "r0": {"Fe": 2.30, "Ni": 2.25},
        "note": "monodentate through one thiolate S; real xanthate can chelate "
                "bidentate — absolute binding here is a lower bound"},
    "oh": {
        "label": "hydroxide OH- (simplest depressant / alkaline-pulp proxy)",
        "role": "baseline", "atoms": HYDROXIDE, "donor": 0,
        "lp": (-0.342, 0.940, 0.000),          # ~110 deg M-O-H
        "charge": -1, "spin": 0,
        "r0": {"Fe": 1.90, "Ni": 1.95},
        "note": "baseline for how much selectivity alkalinity alone can buy"},
    "en": {
        "label": "ethylenediamine NH2CH2CH2NH2 (proxy for DETA, the industrial "
                 "pyrrhotite depressant)",
        "role": "depressant", "atoms": ETHYLENEDIAMINE, "donor": 0,
        "lp": (1.000, 0.000, 0.000),           # N lone pair, in-plane
        "charge": 0, "spin": 0,
        "r0": {"Fe": 2.15, "Ni": 2.10},
        "note": "monodentate through one N; DETA is a TRIDENTATE chelate, so "
                "this proxy understates the real depressant"},
}

PAIRS = [(a, s) for a in ("xanthate", "oh", "en")
         for s in ("pent_fe", "pent_ni", "pyrr_fe")]


def pair_id(ads, sub):
    return f"{ads}@{sub}"


def selected_pairs():
    """FLOT_ONLY filter: pair ids, adsorbate ids, substrate ids or roles."""
    if not ONLY:
        return list(PAIRS)
    want = set(ONLY)
    out = []
    for a, s in PAIRS:
        role = ADSORBATES[a]["role"]
        motif = SUBSTRATES[s]["motif"]
        if (pair_id(a, s) in want or a in want or s in want
                or role in want or motif in want):
            out.append((a, s))
    return out


# ================================================================= utilities ===
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


def npz_path(kind, key):
    return os.path.join(HERE, f"flotation_{kind}_{key.replace('@', '_at_')}.npz")


# ================================================================== geometry ===
def read_cubane():
    """[Fe4S4(SH)4]2- idealized cubane, atom order preserved."""
    with open(SRC_XYZ) as fh:
        lines = fh.read().splitlines()
    n = int(lines[0].split()[0])
    atoms = []
    for i in range(2, 2 + n):
        p = lines[i].split()
        atoms.append([p[0], np.array([float(p[1]), float(p[2]), float(p[3])])])
    if [a[0] for a in atoms[:4]] != ["Fe"] * 4:
        raise ValueError(f"unexpected metal order in {SRC_XYZ}: "
                         f"{[a[0] for a in atoms[:4]]}")
    return atoms


def build_substrate(sid):
    """Cubane -> Fe->Ni swaps -> delete the terminal SH of the site metal.
    Returns (atoms, site_position, outward unit vector of the vacancy)."""
    spec = SUBSTRATES[sid]
    atoms = read_cubane()
    for a in spec["ni_atoms"]:
        atoms[a][0] = "Ni"
    site = spec["site"]
    m = atoms[site][1].copy()
    s_term = atoms[TERM_S[site]][1].copy()
    u = s_term - m
    u = u / np.linalg.norm(u)                 # the vacated M-S bond direction
    drop = {TERM_S[site], TERM_H[site]}
    kept = [(el, xyz) for i, (el, xyz) in enumerate(atoms) if i not in drop]
    if [e for e, _ in kept[:4]] != [atoms[i][0] for i in METAL_ATOMS]:
        raise ValueError("metal block must stay first in the substrate")
    return kept, m, u


def _rot_a_to_b(a, b):
    """Rotation matrix taking unit vector a onto unit vector b (Rodrigues)."""
    a = np.asarray(a, float) / np.linalg.norm(a)
    b = np.asarray(b, float) / np.linalg.norm(b)
    v = np.cross(a, b)
    c = float(np.dot(a, b))
    s = float(np.linalg.norm(v))
    if s < 1e-9:
        if c > 0:
            return np.eye(3)
        p = np.array([1.0, 0.0, 0.0])
        if abs(a[0]) > 0.9:
            p = np.array([0.0, 1.0, 0.0])
        w = np.cross(a, p)
        w /= np.linalg.norm(w)
        return -np.eye(3) + 2.0 * np.outer(w, w)
    vx = np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]])
    return np.eye(3) + vx + vx @ vx * ((1.0 - c) / (s * s))


def _rot_about(axis, ang):
    k = np.asarray(axis, float) / np.linalg.norm(axis)
    kx = np.array([[0, -k[2], k[1]], [k[2], 0, -k[0]], [-k[1], k[0], 0]])
    return (np.eye(3) + np.sin(ang) * kx
            + (1.0 - np.cos(ang)) * (kx @ kx))


def place_adsorbate(sub_atoms, m, u, ads_atoms, spec, r):
    """Put the adsorbate donor at M + r*u with its lone-pair axis pointing back
    at the metal, then pick the azimuth about u that maximises the closest
    non-donor contact with the cluster (cheap steric relief, no relaxation)."""
    donor = spec["donor"]
    pos = np.array([p for _e, p in ads_atoms], float)
    pos = pos - pos[donor]
    R = _rot_a_to_b(np.asarray(spec["lp"], float), -np.asarray(u, float))
    pos = pos @ R.T
    sub_xyz = np.array([p for _e, p in sub_atoms], float)
    anchor = np.asarray(m, float) + r * np.asarray(u, float)
    best, best_gap, best_phi = None, -1.0, 0.0
    for k in range(12):
        phi = 2.0 * np.pi * k / 12.0
        cand = pos @ _rot_about(u, phi).T + anchor
        others = np.delete(cand, donor, axis=0)
        gap = (float(np.min(np.linalg.norm(
            others[:, None, :] - sub_xyz[None, :, :], axis=2)))
            if len(others) else 99.0)
        if gap > best_gap:
            best, best_gap, best_phi = cand, gap, phi
    out = [(e, best[i]) for i, (e, _p) in enumerate(ads_atoms)]
    return out, {"r_M_donor": round(float(r), 3),
                 "azimuth_deg": round(float(np.degrees(best_phi)), 1),
                 "min_nondonor_contact_A": round(float(best_gap), 3)}


def write_xyz(atoms, path, comment, mode="w"):
    with open(path, mode) as fh:
        fh.write(f"{len(atoms)}\n{comment}\n")
        for e, p in atoms:
            fh.write(f"{e} {p[0]:.6f} {p[1]:.6f} {p[2]:.6f}\n")


# ================================================================ mean field ===
def _mol(atoms, charge, spin, maxmem):
    return gto.M(atom=[(e, tuple(float(x) for x in p)) for e, p in atoms],
                 basis=BASIS, charge=charge, spin=spin, verbose=0,
                 max_memory=maxmem)


def _new_mf(mol, xc, restricted=False):
    mf = (dft.RKS(mol) if restricted else dft.UKS(mol)).density_fit()
    mf.xc = xc
    return mf


def _converge(mf, dm0=None, level_shift=0.4, max_cycle=60, tag=""):
    """DIIS with a level shift; on stall, a second-order (Newton) restart —
    the robust recipe of the Fe-S/Ni-S runs in this repo."""
    mf.level_shift = level_shift
    mf.max_cycle = max_cycle
    mf.conv_tol = 1e-7
    if dm0 is not None:
        mf.kernel(dm0=dm0)
    else:
        mf.kernel()
    if not mf.converged:
        print(f"  [{tag}] DIIS incomplete -> Newton (SOSCF) restart", flush=True)
        mf2 = mf.newton()
        mf2.max_cycle = 100
        mf2.kernel(mf.mo_coeff, mf.mo_occ)
        mf2.mol = mf.mol
        return mf2
    return mf


def flip_dm(dm, mol, flip_atoms):
    """Swap alpha/beta density on the diagonal AO blocks of `flip_atoms`."""
    dma, dmb = dm[0].copy(), dm[1].copy()
    aoslice = mol.aoslice_by_atom()
    for a in flip_atoms:
        sl = slice(int(aoslice[a][2]), int(aoslice[a][3]))
        blk = dma[sl, sl].copy()
        dma[sl, sl] = dmb[sl, sl]
        dmb[sl, sl] = blk
    return np.array([dma, dmb])


def metal_spin_pops(mf):
    """Mulliken spin population on metals 0..3 (they stay first in every
    cluster we build) — the BS pattern sanity check."""
    dma, dmb = mf.make_rdm1()
    m = (dma - dmb) @ mf.get_ovlp()
    aoslice = mf.mol.aoslice_by_atom()
    return [round(float(np.trace(m[slice(int(aoslice[a][2]),
                                         int(aoslice[a][3])),
                                   slice(int(aoslice[a][2]),
                                         int(aoslice[a][3]))])), 3)
            for a in METAL_ATOMS]


def bs_cluster_scf(atoms, charge, spin_hs, flip, signs, xc, maxmem, tag,
                   dm_bs=None):
    """Converge the BS Ms=0 reference of a cluster.
    dm_bs given  -> warm start straight into the BS SCF (scan points, PBE0 push).
    dm_bs None   -> HS seed first, then flip the metal AO blocks (cold start).
    Returns (mf, info)."""
    info = {"xc": xc}
    mol_bs = _mol(atoms, charge, 0, maxmem)
    if dm_bs is None:
        t0 = time.time()
        mol_hs = _mol(atoms, charge, spin_hs, maxmem)
        mf_hs = with_timeout(min(STAGE_TIMEOUT, remaining()), _converge,
                             _new_mf(mol_hs, xc), tag=f"{tag}-HS")
        info["hs"] = {"e": float(mf_hs.e_tot), "2S": spin_hs,
                      "converged": bool(mf_hs.converged),
                      "seconds": round(time.time() - t0, 1)}
        dm0 = flip_dm(mf_hs.make_rdm1(), mol_hs, flip)
        lvl = 0.3
    else:
        dm0 = dm_bs
        lvl = 0.2
    t0 = time.time()
    mf = with_timeout(min(STAGE_TIMEOUT, remaining()), _converge,
                      _new_mf(mol_bs, xc), dm0=dm0, level_shift=lvl,
                      tag=f"{tag}-BS")
    ss, _ = mf.spin_square()
    pops = metal_spin_pops(mf)
    signs_ok = all((p > 0) == (s > 0) for p, s in zip(pops, signs))
    polarized = all(abs(p) > 1.0 for p in pops)
    info.update({"e": float(mf.e_tot), "converged": bool(mf.converged),
                 "s_squared": round(float(ss), 3),
                 "metal_mulliken_spin": pops,
                 "expected_signs": list(signs),
                 "bs_pattern_ok": bool(mf.converged and signs_ok and polarized),
                 "seconds": round(time.time() - t0, 1)})
    return mf, info


def molecule_scf(atoms, charge, xc, maxmem, tag, relax=False):
    """Closed-shell free adsorbate; optional cheap PBE relaxation."""
    info = {"xc": xc, "relaxed": False}
    mol = _mol(atoms, charge, 0, maxmem)
    mf = _new_mf(mol, xc, restricted=True)
    mf = _converge(mf, level_shift=0.0, tag=tag)
    if relax:
        try:
            from pyscf.geomopt.geometric_solver import optimize
            mol_eq = with_timeout(min(STAGE_TIMEOUT, remaining()),
                                  optimize, mf, maxsteps=60)
            atoms = [(mol_eq.atom_symbol(i),
                      np.asarray(mol_eq.atom_coord(i, unit="Angstrom")))
                     for i in range(mol_eq.natm)]
            mf = _converge(_new_mf(_mol(atoms, charge, 0, maxmem), xc,
                                   restricted=True), tag=f"{tag}-eq")
            info["relaxed"] = True
        except Exception as e:
            info["relax_note"] = f"not relaxed ({type(e).__name__}: {e})"
    info.update({"e": float(mf.e_tot), "converged": bool(mf.converged)})
    return atoms, info


def constrained_opt(atoms, n_frozen, charge, spin_hs, flip, signs, maxmem, tag,
                    dm_bs=None):
    """FLOT_MODE=opt: geomeTRIC relaxation of the adsorbate only (all substrate
    atoms frozen) on the BS-PBE surface. Best effort; caller falls back."""
    from pyscf.geomopt.geometric_solver import optimize
    mf, _info = bs_cluster_scf(atoms, charge, spin_hs, flip, signs, "pbe",
                               maxmem, tag, dm_bs=dm_bs)
    with tempfile.NamedTemporaryFile("w", suffix=".txt", delete=False) as fh:
        fh.write(f"$freeze\nxyz 1-{n_frozen}\n")
        cons = fh.name
    mol_eq = with_timeout(min(STAGE_TIMEOUT, remaining()), optimize, mf,
                          maxsteps=25, constraints=cons)
    new = [(mol_eq.atom_symbol(i),
            np.asarray(mol_eq.atom_coord(i, unit="Angstrom")))
           for i in range(mol_eq.natm)]
    return new


# ================================================================== the pairs ===
def scan_grid(r0):
    if MODE == "rigid" or SCAN_N == 1:
        return [r0]
    half = (SCAN_N - 1) / 2.0
    step = SCAN_DR / half if half else 0.0
    return [round(r0 + (k - half) * step, 3) for k in range(SCAN_N)]


def run_pair(ads_id, sub_id, ads_geom, maxmem):
    """One substrate-adsorbate-site job at PBE: build, scan the M-donor bond,
    keep the lowest COMPUTED point. Returns (record, best_atoms, best_dm)."""
    aspec, sspec = ADSORBATES[ads_id], SUBSTRATES[sub_id]
    sub_atoms, m, u = build_substrate(sub_id)
    r0 = aspec["r0"][sspec["site_el"]]
    charge = sspec["charge"] + aspec["charge"]
    rec = {"adsorbate": aspec["label"], "substrate": sspec["label"],
           "site_element": sspec["site_el"], "site_atom": sspec["site"],
           "motif": sspec["motif"], "role": aspec["role"],
           "charge": charge, "hs_seed_2S": sspec["spin_hs"],
           "bs_flip_atoms": list(sspec["bs_flip"]),
           "r0_guess_A": r0, "mode": MODE,
           "n_substrate_atoms": len(sub_atoms)}
    points, dm, best = [], None, None
    best_atoms = None
    for r in scan_grid(r0):
        if remaining() < 900:
            points.append({"r": r, "skipped": "time budget"})
            break
        placed, pinfo = place_adsorbate(sub_atoms, m, u, ads_geom, aspec, r)
        atoms = list(sub_atoms) + placed
        try:
            mf, info = bs_cluster_scf(atoms, charge, sspec["spin_hs"],
                                      sspec["bs_flip"], sspec["bs_signs"],
                                      "pbe", maxmem, f"{ads_id}@{sub_id}",
                                      dm_bs=dm)
        except StageTimeout:
            points.append({"r": r, "timed_out": True})
            print(f"  [{ads_id}@{sub_id}] r={r} TIMED OUT", flush=True)
            continue
        info.update(pinfo)
        points.append(info)
        print(f"  [{ads_id}@{sub_id}] r={r:.2f} E={info['e']:.5f} "
              f"conv={info['converged']} bs_ok={info['bs_pattern_ok']} "
              f"spins={info['metal_mulliken_spin']} ({info['seconds']}s)",
              flush=True)
        if info["converged"]:
            dm = mf.make_rdm1()
            if best is None or info["e"] < best["e"]:
                best, best_atoms = info, atoms
    rec["scan"] = points
    if best is None:
        rec["status"] = "no converged scan point"
        return rec, None, None
    rec["pbe"] = {"e_complex": best["e"], "r_best_A": best["r_M_donor"],
                  "bs_pattern_ok": best["bs_pattern_ok"],
                  "s_squared": best["s_squared"],
                  "metal_mulliken_spin": best["metal_mulliken_spin"],
                  "from": "lowest computed scan point"}
    ok = [(p["r_M_donor"], p["e"]) for p in points
          if p.get("converged") and "r_M_donor" in p]
    if len(ok) >= 3:
        rr = np.array([x for x, _ in ok])
        ee = np.array([y for _, y in ok])
        a, b, _c = np.polyfit(rr, ee, 2)
        if a > 0:
            rmin = float(-b / (2 * a))
            rec["pbe"]["parabola"] = {
                "r_min_A": round(rmin, 3),
                "inside_scan_range": bool(rr.min() <= rmin <= rr.max()),
                "note": "interpolation only; dE_ads uses the computed point"}
    # FLOT_MODE=opt: relax the adsorbate on top of the best rigid point
    if MODE == "opt" and best_atoms is not None and remaining() > 3600:
        try:
            n_sub = len(sub_atoms)
            relaxed = constrained_opt(best_atoms, n_sub, charge,
                                      sspec["spin_hs"], sspec["bs_flip"],
                                      sspec["bs_signs"], maxmem,
                                      f"{ads_id}@{sub_id}-opt", dm_bs=dm)
            mf, info = bs_cluster_scf(relaxed, charge, sspec["spin_hs"],
                                      sspec["bs_flip"], sspec["bs_signs"],
                                      "pbe", maxmem, f"{ads_id}@{sub_id}-opt")
            if info["converged"] and info["e"] < best["e"]:
                mpos = np.asarray(m, float)
                dpos = np.asarray(relaxed[len(sub_atoms) + aspec["donor"]][1])
                rec["pbe"].update({
                    "e_complex": info["e"],
                    "r_best_A": round(float(np.linalg.norm(dpos - mpos)), 3),
                    "bs_pattern_ok": info["bs_pattern_ok"],
                    "from": "constrained opt (substrate frozen, adsorbate free)"})
                best_atoms, dm = relaxed, mf.make_rdm1()
        except Exception as e:
            rec["opt_note"] = (f"constrained opt failed, kept the rigid scan "
                               f"minimum ({type(e).__name__}: {e})")
    else:
        rec["relaxation_note"] = (
            "rigid scan of the M-donor bond on a frozen cluster framework and a "
            "frozen adsorbate geometry: dE_ads is essentially an interaction "
            "energy (no substrate or adsorbate reorganisation)")
    rec["status"] = "ok"
    return rec, best_atoms, dm


def pbe0_push(res, maxmem):
    """PBE0//PBE single points on the pairs (and their substrates) needed for
    the chosen descriptors. PBE density is the SCF start, so no HS seed."""
    push = res.setdefault("pbe0_push", {"set": PBE0_SET, "pairs": {},
                                        "substrates": {}, "adsorbates": {}})
    wanted = []
    if PBE0_SET in ("collector", "both"):
        wanted.append("xanthate")
    if PBE0_SET in ("depressant", "both"):
        wanted.append("en")
    for ads_id in wanted:
        # best (most negative dE_ads) pentlandite site + the pyrrhotite site
        de = {s: _de(res, "pbe", ads_id, s) for s in ("pent_fe", "pent_ni")}
        de = {s: v for s, v in de.items() if v is not None}
        if not de or _de(res, "pbe", ads_id, "pyrr_fe") is None:
            push["pairs"][ads_id] = {"skipped": "PBE half of the pair missing"}
            continue
        best_sub = min(de, key=de.get)
        for sub_id in (best_sub, "pyrr_fe"):
            pid = pair_id(ads_id, sub_id)
            if pid in push["pairs"] and push["pairs"][pid].get("e"):
                continue
            if remaining() < 2400:
                push["pairs"][pid] = {"skipped": "time budget"}
                continue
            dmz = npz_path("dm", pid)
            gz = npz_path("geom", pid)
            if not (os.path.exists(dmz) and os.path.exists(gz)):
                push["pairs"][pid] = {"skipped": "no cached PBE density/geometry"}
                continue
            z = np.load(gz, allow_pickle=True)
            atoms = [(str(e), np.asarray(p, float))
                     for e, p in zip(z["el"], z["xyz"])]
            sspec = SUBSTRATES[sub_id]
            charge = sspec["charge"] + ADSORBATES[ads_id]["charge"]
            try:
                _mf, info = bs_cluster_scf(
                    atoms, charge, sspec["spin_hs"], sspec["bs_flip"],
                    sspec["bs_signs"], "pbe0", maxmem, f"pbe0-{pid}",
                    dm_bs=np.load(dmz)["dm"])
                push["pairs"][pid] = info
                print(f"[pbe0] {pid}: E={info['e']:.5f} "
                      f"conv={info['converged']}", flush=True)
            except StageTimeout:
                push["pairs"][pid] = {"timed_out": True}
            save(res)
            # matching substrate at PBE0
            if sub_id not in push["substrates"] and remaining() > 2400:
                sub_atoms, _m, _u = build_substrate(sub_id)
                sz = npz_path("sub", sub_id)
                dm0 = np.load(sz)["dm"] if os.path.exists(sz) else None
                try:
                    _mfs, sinfo = bs_cluster_scf(
                        sub_atoms, sspec["charge"], sspec["spin_hs"],
                        sspec["bs_flip"], sspec["bs_signs"], "pbe0", maxmem,
                        f"pbe0-{sub_id}", dm_bs=dm0)
                    push["substrates"][sub_id] = sinfo
                except StageTimeout:
                    push["substrates"][sub_id] = {"timed_out": True}
                save(res)
        if ads_id not in push["adsorbates"]:
            geom = res["adsorbates"][ads_id]["geometry"]
            atoms = [(e, np.array(p, float)) for e, p in geom]
            _a, ainfo = molecule_scf(atoms, ADSORBATES[ads_id]["charge"],
                                     "pbe0", maxmem, f"pbe0-{ads_id}")
            push["adsorbates"][ads_id] = ainfo
            save(res)
    return push


# ================================================================ descriptors ===
def _de(res, level, ads_id, sub_id):
    """dE_ads in kcal/mol at the requested level, or None if incomplete."""
    if level == "pbe":
        pr = res.get("pairs", {}).get(pair_id(ads_id, sub_id), {}).get("pbe")
        sb = res.get("substrates", {}).get(sub_id, {}).get("pbe")
        ad = res.get("adsorbates", {}).get(ads_id, {}).get("pbe")
        ec = pr.get("e_complex") if pr else None
    else:
        push = res.get("pbe0_push", {})
        pr = push.get("pairs", {}).get(pair_id(ads_id, sub_id))
        sb = push.get("substrates", {}).get(sub_id)
        ad = push.get("adsorbates", {}).get(ads_id)
        ec = pr.get("e") if pr else None
    if ec is None or not sb or not ad:
        return None
    if sb.get("e") is None or ad.get("e") is None:
        return None
    return round((ec - sb["e"] - ad["e"]) * KCAL, 2)


def descriptors(res, level):
    """dE_ads table + the two double differences that decide the case."""
    out = {"level": level, "dE_ads_kcal": {}, "dd_E_kcal": {}}
    for a, s in PAIRS:
        v = _de(res, level, a, s)
        if v is not None:
            out["dE_ads_kcal"][pair_id(a, s)] = v
    tbl = out["dE_ads_kcal"]

    def dd(ads_id, sign_motif):
        """sign_motif='pent' -> pentlandite-preferring is negative;
        'pyrr' -> pyrrhotite-preferring is negative."""
        pyr = tbl.get(pair_id(ads_id, "pyrr_fe"))
        sites = {s: tbl[pair_id(ads_id, s)] for s in ("pent_fe", "pent_ni")
                 if pair_id(ads_id, s) in tbl}
        if pyr is None or not sites:
            return None
        best_site = min(sites, key=sites.get)
        pent = sites[best_site]
        val = (pent - pyr) if sign_motif == "pent" else (pyr - pent)
        return {"value_kcal": round(val, 2),
                "pentlandite_site_used": best_site,
                "per_site": {s: round((v - pyr) if sign_motif == "pent"
                                      else (pyr - v), 2)
                             for s, v in sites.items()},
                "negative_means": ("collector prefers PENTLANDITE (selectivity "
                                   "works)" if sign_motif == "pent" else
                                   "depressant prefers PYRRHOTITE (depression "
                                   "works)")}

    for key, ads_id, motif in (("collector_xanthate", "xanthate", "pent"),
                               ("depressant_en", "en", "pyrr"),
                               ("baseline_oh", "oh", "pyrr")):
        d = dd(ads_id, motif)
        if d:
            out["dd_E_kcal"][key] = d
    # site contrast inside pentlandite: is the selectivity Fe- or Ni-borne?
    for a in ("xanthate", "oh", "en"):
        fe, ni = tbl.get(pair_id(a, "pent_fe")), tbl.get(pair_id(a, "pent_ni"))
        if fe is not None and ni is not None:
            out.setdefault("pentlandite_site_contrast_kcal", {})[a] = round(
                ni - fe, 2)
    return out


def make_verdict(res):
    """Read the signs. PBE decides magnitude-free direction; PBE0 arbitrates."""
    pbe = res.get("descriptors", {}).get("pbe", {}).get("dd_E_kcal", {})
    hyb = res.get("descriptors", {}).get("pbe0", {}).get("dd_E_kcal", {})
    v = {"rule": "both dd_E must be NEGATIVE for the separation to work: "
                 "collector binds pentlandite harder AND depressant binds "
                 "pyrrhotite harder",
         "signs": {}}
    for key in ("collector_xanthate", "depressant_en", "baseline_oh"):
        p = pbe.get(key, {}).get("value_kcal")
        h = hyb.get(key, {}).get("value_kcal")
        if p is None:
            continue
        rec = {"pbe_kcal": p}
        if h is not None:
            rec["pbe0_kcal"] = h
            rec["sign_agrees"] = bool((p < 0) == (h < 0))
            rec["read"] = ("PBE and PBE0 AGREE on the sign — the direction of "
                           "the descriptor survives the functional change "
                           "(magnitudes from PBE remain unreliable, cf. the "
                           "12 kcal/mol spin-gap spread of the N6 gate)"
                           if rec["sign_agrees"] else
                           "PBE and PBE0 DISAGREE on the sign — NOT resolvable "
                           "at this level; the pair needs a correlated method "
                           "(NEVPT2 on the BS-UNO-CASSCF reference), exactly "
                           "the N6 conclusion carried into adsorption")
        else:
            rec["read"] = "PBE only; the sign has NOT been checked on a hybrid"
        rec["direction"] = ("favourable (negative)" if p < 0
                            else "UNFAVOURABLE (positive) — this reagent class "
                                 "does not separate the two motifs in the "
                                 "expected direction on this model")
        v["signs"][key] = rec
    coll = v["signs"].get("collector_xanthate", {}).get("pbe_kcal")
    depr = v["signs"].get("depressant_en", {}).get("pbe_kcal")
    if coll is not None and depr is not None:
        if coll < 0 and depr < 0:
            v["gate"] = ("PASS at PBE: the xanthate/amine pair pulls in OPPOSITE "
                         "directions on the two motifs — the two-reagent scheme "
                         "(collector + amine depressant) has an electronic basis "
                         "on this model, and the next question is magnitude, not "
                         "direction")
        elif coll < 0 or depr < 0:
            v["gate"] = ("PARTIAL at PBE: only one of the two levers points the "
                         "right way; the separation would rest on that lever "
                         "alone on this model")
        else:
            v["gate"] = ("FAIL at PBE: neither lever points the expected way on "
                         "this cluster model — either the motif difference is "
                         "too small for a 15-atom cluster to carry, or the "
                         "selectivity is not molecular-adsorption-borne (surface "
                         "oxidation / activation ions / pulp chemistry)")
    else:
        v["gate"] = "incomplete: not enough converged pairs to call the gate"
    v["site_note"] = ("pentlandite_site_contrast = dE_ads(Ni site) - "
                      "dE_ads(Fe site): negative means the Ni site binds the "
                      "reagent harder, i.e. the selectivity is Ni-BORNE and a "
                      "reagent should be designed for Ni; positive means the "
                      "difference is carried by the Fe site environment instead")
    return v


HONESTY = [
    "Кластер [M4S4(SH)3]- НЕ решётка минерала: нет маделунговского поля, нет "
    "реконструкции поверхности, нет разупорядочения Fe/Ni и вакансий серы "
    "реального пирротина Fe7S8.",
    "Газовая фаза: нет воды, нет двойного электрического слоя, нет pH, нет "
    "потенциала пульпы, нет конкурирующей активации Cu2+/Ca2+ — всё это в "
    "реальной флотации первого порядка.",
    "Только электронные энергии: без ZPE, без температуры и энтропии, без "
    "сольватации. Это НЕ изотерма адсорбции и НЕ свободная энергия.",
    "Монодентатное связывание (одна вакансия на кластер). Ксантогенат бывает "
    "бидентатным, DETA — тридентатный хелат, а en только его двухазотный "
    "прокси: абсолютные числа ЗАНИЖАЮТ амин.",
    "BS-решения спин-контаминированы, спиновая проекция не применялась.",
    "BSSE (counterpoise) не считался; в dd_E он в основном сокращается, а "
    "энергия свободного аниона сокращается ТОЧНО — поэтому дескриптор и взят "
    "как двойная разность, а не как энергия связывания.",
    "Каркас кластера заморожен; в режиме scan заморожен и адсорбат, поэтому "
    "dE_ads по смыслу ближе к энергии взаимодействия, чем к энергии адсорбции "
    "с реорганизацией.",
    "PBE завышает магнитуды на Fe-S (урок N6-гейта: разброс PBE/PBE0 по "
    "спиновой щели ~12 ккал/моль) — поэтому смотрим ЗНАК, и знак лучшей пары "
    "перепроверяем на гибриде PBE0//PBE.",
    "ЭТО ЗНАКОВЫЙ СКРИН ДЕСКРИПТОРА. Он может сказать «этот класс реагентов "
    "предпочитает этот мотив и вот на каком сайте» — и не может сказать "
    "«подавать X г/т при pH 9».",
]


# ======================================================================= main ===
def geoms_selfcheck():
    """--geoms: build every job without pyscf and print the table."""
    print(f"{'pair':>22}  {'chg':>4} {'2S_HS':>5} {'nat':>4} {'r0':>5} "
          f"{'minGap':>7}")
    for a, s in selected_pairs():
        aspec, sspec = ADSORBATES[a], SUBSTRATES[s]
        sub_atoms, m, u = build_substrate(s)
        r0 = aspec["r0"][sspec["site_el"]]
        placed, pinfo = place_adsorbate(sub_atoms, m, u, aspec["atoms"],
                                        aspec, r0)
        n = len(sub_atoms) + len(placed)
        print(f"{pair_id(a, s):>22}  {sspec['charge'] + aspec['charge']:>4} "
              f"{sspec['spin_hs']:>5} {n:>4} {r0:>5.2f} "
              f"{pinfo['min_nondonor_contact_A']:>7.2f}")
    if not HAVE_PYSCF:
        print(f"\n[note] pyscf not importable here ({PYSCF_ERR}) — geometry "
              f"check only")


def main():
    if "--geoms" in sys.argv:
        geoms_selfcheck()
        return
    if not HAVE_PYSCF:
        sys.exit(f"[deps] pyscf is required:  python3 -m pip install pyscf "
                 f"({PYSCF_ERR})")

    prev = load_prev()
    resume = bool(prev and prev.get("meta", {}).get("job"))
    res = prev if resume else {"status": "running"}
    res.setdefault("meta", {})
    res.setdefault("adsorbates", {})
    res.setdefault("substrates", {})
    res.setdefault("pairs", {})
    pairs = selected_pairs()
    res["meta"].update({
        "job": "flotation reagent adsorption: pentlandite vs pyrrhotite motif",
        "purpose": "ГЕЙТ КЕЙСА N9 (APPLIED_WIDE_SCAN.md): селективный депрессор "
                   "пирротина. Разделение пентландита (Fe,Ni)9S8 и пирротина "
                   "Fe7S8 убирает серу ДО трубы — до Серной программы и дешевле "
                   "любого улавливания SO2. Вопрос гейта: предпочитает ли "
                   "реагент один мотив другому и в какую сторону.",
        "chain": "calc/fe_ni_s_pentlandite.py (N6 gate: BS pattern A, Fe "
                 "+2.03/+2.03, Ni -1.27/-1.26; spin gap PBE -18.8 / PBE0 -30.7 "
                 "kcal/mol) -> this job",
        "level": "UKS-PBE/def2-SVP + density fitting, broken-symmetry reference; "
                 "PBE0//PBE sign check on the best pair",
        "basis": BASIS, "mode": MODE, "scan_points": SCAN_N,
        "scan_half_range_A": SCAN_DR,
        "mf_mode": MF_MODE, "pbe0_set": PBE0_SET,
        "only": ONLY or "all pairs",
        "pairs_selected": [pair_id(a, s) for a, s in pairs],
        "geometry_source": "femoco_4fe4s_cubane.xyz; Fe->Ni swaps for the "
                           "pentlandite motif; one terminal SH- removed at the "
                           "site metal to create the open surface site",
        "descriptor": "dd_E = double difference of dE_ads between the two "
                      "motifs; the free-adsorbate energy cancels EXACTLY in it",
        "stage_timeout_s": STAGE_TIMEOUT, "total_budget_s": TOTAL_BUDGET,
        "resumed": resume})
    res["honesty"] = HONESTY
    save(res)

    try:
        maxmem = int(os.environ.get("PYSCF_MAX_MEMORY", "16000"))
        need_ads = sorted({a for a, _s in pairs})
        need_sub = sorted({s for _a, s in pairs})
        do_pbe = MF_MODE in ("pbe", "both")
        do_pbe0 = MF_MODE in ("pbe", "pbe0", "both")
        print(f"[plan] pairs={[pair_id(a, s) for a, s in pairs]} "
              f"mf={MF_MODE} mode={MODE} maxmem={maxmem}MB", flush=True)

        # ---- free adsorbates ------------------------------------------------
        ads_geom = {}
        for a in need_ads:
            spec = ADSORBATES[a]
            if a in res["adsorbates"] and res["adsorbates"][a].get("pbe"):
                g = res["adsorbates"][a]["geometry"]
                ads_geom[a] = [(e, np.array(p, float)) for e, p in g]
                continue
            atoms = [(e, np.array(p, float)) for e, p in spec["atoms"]]
            atoms, info = molecule_scf(atoms, spec["charge"], "pbe", maxmem,
                                       f"ads-{a}", relax=ADS_RELAX)
            ads_geom[a] = atoms
            res["adsorbates"][a] = {
                "label": spec["label"], "role": spec["role"],
                "charge": spec["charge"], "note": spec["note"],
                "donor_atom": spec["donor"],
                "geometry": [[e, [round(float(x), 6) for x in p]]
                             for e, p in atoms],
                "pbe": info}
            print(f"[ads] {a}: E={info['e']:.6f} relaxed={info['relaxed']}",
                  flush=True)
            save(res)

        # ---- substrates (BS reference on each open-site cluster) ------------
        for s in need_sub:
            spec = SUBSTRATES[s]
            if s in res["substrates"] and res["substrates"][s].get("pbe"):
                continue
            sub_atoms, _m, _u = build_substrate(s)
            write_xyz(sub_atoms, os.path.join(HERE, f"flotation_sub_{s}.xyz"),
                      spec["label"])
            print(f"[sub] {s}: {spec['label']}", flush=True)
            mf, info = bs_cluster_scf(sub_atoms, spec["charge"], spec["spin_hs"],
                                      spec["bs_flip"], spec["bs_signs"], "pbe",
                                      maxmem, f"sub-{s}")
            np.savez_compressed(npz_path("sub", s), dm=mf.make_rdm1())
            res["substrates"][s] = {
                "label": spec["label"], "motif": spec["motif"],
                "site_element": spec["site_el"], "site_atom": spec["site"],
                "charge": spec["charge"], "formal_state": spec["formal"],
                "hs_seed_2S": spec["spin_hs"],
                "bs_flip_atoms": list(spec["bs_flip"]),
                "n_atoms": len(sub_atoms), "pbe": info}
            print(f"[sub] {s}: E={info['e']:.5f} conv={info['converged']} "
                  f"bs_ok={info['bs_pattern_ok']} "
                  f"spins={info['metal_mulliken_spin']}", flush=True)
            save(res)

        # ---- the pairs ------------------------------------------------------
        if do_pbe:
            for a, s in pairs:
                pid = pair_id(a, s)
                if res["pairs"].get(pid, {}).get("status") == "ok":
                    print(f"[resume] {pid} already done", flush=True)
                    continue
                if remaining() < 1200:
                    res["pairs"][pid] = {"skipped": "time budget"}
                    save(res)
                    continue
                print(f"[pair] {pid}", flush=True)
                rec, atoms, dm = run_pair(a, s, ads_geom[a], maxmem)
                res["pairs"][pid] = rec
                if atoms is not None:
                    write_xyz(atoms, GEOMS, f"{pid} r={rec['pbe']['r_best_A']}A",
                              mode="a")
                    np.savez_compressed(npz_path("dm", pid), dm=dm)
                    np.savez_compressed(
                        npz_path("geom", pid),
                        el=np.array([e for e, _p in atoms]),
                        xyz=np.array([p for _e, p in atoms], float))
                save(res)
                res["descriptors"] = {"pbe": descriptors(res, "pbe")}
                save(res)

        res.setdefault("descriptors", {})["pbe"] = descriptors(res, "pbe")
        save(res)

        # ---- PBE0//PBE sign check ------------------------------------------
        if do_pbe0 and remaining() > 2400:
            print(f"[pbe0] push set={PBE0_SET}", flush=True)
            pbe0_push(res, maxmem)
            res["descriptors"]["pbe0"] = descriptors(res, "pbe0")
            save(res)
        elif do_pbe0:
            res.setdefault("pbe0_push", {})["skipped"] = "time budget"

        res["verdict"] = make_verdict(res)
        res["status"] = "ok"
    except Exception as e:
        res["status"] = f"aborted: {type(e).__name__}: {e}"
        print("ABORT:", res["status"], flush=True)
    finally:
        res.setdefault("meta", {})["total_seconds"] = round(elapsed(), 1)
        save(res)
        print(f"\nwrote {RESULTS}  status={res.get('status')}  "
              f"(total {res['meta']['total_seconds']}s)", flush=True)


if __name__ == "__main__":
    main()
