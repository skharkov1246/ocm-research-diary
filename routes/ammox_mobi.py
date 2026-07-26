#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
routes/ammox_mobi.py — Track L7 OPENER (AWS job): propylene ammoxidation to
acrylonitrile on a Bi-Mo-O (Sohio) site — our OCM selectivity descriptor
transplanted 1:1.

WHY THIS TRACK (APPLIED_WIDE_SCAN.md, L7). Acrylonitrile at Saratovorgsintez is
made by the Sohio process on Mo-Bi oxides. The whole economics of that plant is
one number: the selectivity of C3H6 -> acrylonitrile against total oxidation of
the intermediate/product. That is EXACTLY the question our OCM stage answered with
ddE# = barrier(C2H6) - barrier(CH4) on a Mn=O radical-oxygen site: does the active
oxygen spare the target path, or does it eat the product faster than the feed?
Here the same descriptor is computed on a Bi-O-Mo bridge. This is an OPERATING
plant — every +1% acrylonitrile yield is money with no new construction, which is
why L7 sits above most greenfield ideas in the second-echelon priority list.

THE DESCRIPTOR (identical convention to routes/ocm_mnw_selectivity.py):
    ddE# = barrier(over-oxidation) - barrier(target allylic H abstraction)
         = barrier(acrolein  C-H)  - barrier(propylene allylic C-H)
    ddE# > 0  the site SPARES the target path -> selectivity toward acrylonitrile
    ddE# < 0  the site attacks the aldehyde/allyl intermediate faster than the
              feed -> the BEP wall, combustion to COx/HCN
Target step: allylic H abstraction from propylene by the Mo=O oxygen (the
rate-determining, selectivity-determining step of the Sohio mechanism).
Over-oxidation step: H abstraction from the aldehydic C-H of acrolein — the
gateway to total oxidation of the partially oxidized intermediate.

THE CLUSTER (generated in this script from idealized bond lengths):
  a neutral Bi-O-Mo bridge, the minimal Sohio active-site proxy
      (HO)(O=)Mo(mu-O)2Bi(OH)2   =  MoBiO6H3, Mo(VI) d0 + Bi(III), charge 0
  built from: Mo=O 1.70 A (the reactive oxo, ON the HAT axis), Mo-OH 1.90 A,
  Mo-(mu-O) 1.95 A, Bi-(mu-O) 2.20 A, Bi-OH 2.15 A, O-H 0.97 A.
  Atom order is FIXED and identical for both substrates so that the atom-indexed
  AVAS and the geomeTRIC distance pins are the same object in every structure:
      0  Mo          (the oxo metal)
      1  O           (the reactive oxo — the abstracting oxygen)
      2  H           (the transferring hydrogen)
      3  C           (the carbon it comes from: propylene CH3 / acrolein CHO)
      4..            the rest of the substrate
      last 9         the Bi-O frame (2 mu-O, Bi, Mo-OH, 2 Bi-OH) — FROZEN during
                     the scan (embedded-cluster approximation: the lattice holds
                     that shell), so the frame is identical at every scan point
                     and the profile is internally consistent.

THE PROTOCOL (the OCM machine, unchanged):
  1. spins   — UKS-PBE0 ladder over the parity-correct multiplicities at the R
               geometry: the working spin state is MEASURED, not postulated.
  2. geom    — distinguished-coordinate scan: 2D pin r(C-H) + r(O-H) (geomeTRIC),
               reactive centre relaxed, frame frozen. This gives an UPPER BOUND on
               the barrier (a relaxed scan is not a saddle-point optimization) —
               stated in the honesty block, not hidden.
  3. polish  — warm-started chain + UKS internal stability along the profile
               (kills the state-hopping sawtooth of cold SCF starts).
  4. profile — ONE fixed active space for every structure on the path:
               ROHF -> fixed AVAS on Mo 4d / O 2p / transferring-H 1s ->
               CASSCF -> SC-NEVPT2 at EVERY point, so the correlated surface picks
               its OWN R and TS instead of inheriting the DFT indices.
  5. merge   — MIRAGE-DETECTOR GATES, all of them mandatory:
               * outlier mask: a point more than 25 kcal/mol above BOTH neighbours
                 is dropped (the Stage-18 Cr spike lesson);
               * interior TS: rel[ts] >= max(rel[0], rel[-1]) — an edge maximum
                 means the scan never crossed the pass;
               * sanity window: a barrier outside (0, 150) kcal/mol is not a
                 barrier, it is a broken structure;
               * NO-VERDICT IS A LEGAL OUTCOME. If any gate fails at the NEVPT2
                 level for either substrate, the JSON says so and NO ddE# is
                 published. We would rather ship nothing than ship a mirage.

HONEST SCOPE (also written into the results JSON):
  - a 19-20 atom GAS-PHASE cluster, not the periodic Mo-Bi-O catalyst: there is no
    lattice oxygen reservoir, no Mars-van-Krevelen refill step, no NH3 activation,
    no Bi/Mo phase (alpha-Bi2Mo3O12 etc.), no coverage, no temperature.
  - a relaxed distinguished-coordinate scan gives an UPPER estimate of the barrier
    (the true saddle is at or below the scan maximum); no frequencies, no ZPE.
  - def2-SVP + ECPs on Mo and Bi; CASSCF active space is fixed by construction,
    not by threshold tuning — the size is stable across structures by design.
  - the number that transfers is the SIGN and MAGNITUDE of ddE#, not the absolute
    barriers.

ENV KNOBS:
  AMMOX_STAGE_TIMEOUT  seconds per scan point / per correlated point (5400)
  AMMOX_LEVEL          dft | cas | nevpt2   (default: all = nevpt2)
  AMMOX_NPTS           points on the path including R (default 7)
  AMMOX_SPIN           force 2S (default: taken from the spins ladder, else 0)
  AMMOX_BASIS          def2-svp
  AMMOX_MAXMEM         pyscf max_memory MB (the AWS wrapper sets it)
  AMMOX_NDOCC/AMMOX_NVIR  fixed-AVAS composition (default 5/5 -> CAS(10e,10o))

STAGES (independent, incremental JSON — a re-launched box RESUMES):
  python3 -u routes/ammox_mobi.py spins
  python3 -u routes/ammox_mobi.py geom    c3h6|acrolein
  python3 -u routes/ammox_mobi.py polish  c3h6|acrolein
  python3 -u routes/ammox_mobi.py profile c3h6|acrolein
  python3 -u routes/ammox_mobi.py merge
  python3 -u routes/ammox_mobi.py all            # everything, in order (AWS)

OUTPUT:
  routes/ammox_mobi_{sub}_geom.json     scan geometries + DFT profile
  routes/ammox_mobi_{sub}_profile.json  CASSCF / NEVPT2 profile
  routes/ammox_mobi_spins.json          the spin ladder
  routes/ammox_mobi_results.json        the deliverable (descriptor + gates)
"""
import json
import math
import os
import signal
import sys
import time

try:
    import numpy as np
except ImportError:
    sys.exit("[deps] numpy is required:  python3 -m pip install numpy")
try:
    from pyscf import dft, fci, gto, mcscf, scf
    from pyscf.mrpt import NEVPT
except ImportError as e:                                  # pyscf optional at lint
    sys.exit(f"[deps] pyscf is required:  python3 -m pip install pyscf   ({e})")

HARTREE_KCAL = 627.509474
DIR = os.path.dirname(os.path.abspath(__file__))
PREF = "ammox_mobi"
BASIS = os.environ.get("AMMOX_BASIS", "def2-svp")
MAXMEM = int(os.environ.get("AMMOX_MAXMEM") or os.environ.get("PYSCF_MAX_MEMORY")
             or "12000")
XC = "pbe0"
LEVEL = os.environ.get("AMMOX_LEVEL", "all").lower()      # dft | cas | nevpt2|all
NPTS = int(os.environ.get("AMMOX_NPTS", "7"))
STAGE_TIMEOUT = int(os.environ.get("AMMOX_STAGE_TIMEOUT", "5400"))
N_DOCC = int(os.environ.get("AMMOX_NDOCC", "5"))
N_VIR = int(os.environ.get("AMMOX_NVIR", "5"))
SUBSTRATES = ("c3h6", "acrolein")

# mirage-detector constants (same numbers as routes/ocm_mnw_robust_readout.py)
OUTLIER_KCAL = 25.0
BARRIER_MIN, BARRIER_MAX = 0.0, 150.0

# --- Bi-Mo-O frame bond lengths (Angstrom) ---------------------------------
D_MO_OXO = 1.70     # Mo=O, the reactive oxo (on the HAT axis)
D_MO_OH = 1.90      # Mo-OH
D_MO_OBR = 1.95     # Mo-(mu-O)
D_BI_OBR = 2.20     # Bi-(mu-O)
D_BI_OH = 2.15      # Bi-OH
D_OH = 0.97

# HAT grid: (r_CH, r_OH) pins in Angstrom. The R point pins only r(O-H) (the
# pre-complex). Same shape as the OCM grid so the two tracks stay comparable.
R_POINT = (None, 2.20)
FULL_PATH = [(1.15, 1.70), (1.20, 1.50), (1.25, 1.35), (1.30, 1.25),
             (1.35, 1.17), (1.45, 1.08), (1.60, 1.00), (1.80, 0.97)]
PATH = FULL_PATH[:max(2, NPTS - 1)]
GEOM_VERSION = "ammox-mobi-frozen-frame-v1"

SPIN = int(os.environ["AMMOX_SPIN"]) if os.environ.get("AMMOX_SPIN") else None
CAS_TARGET = None            # set once SPIN is known: (nelec, norb)
AVAS_LABELS = ["Mo 4d", "O 2p", "2 H 1s"]   # atom 2 (0-based) = transferring H
# NB: "Mo 3d" in the track brief means the metal VALENCE d shell; for Mo that is
# 4d (3d is core, and with def2-SVP it sits inside the ECP-28). Recorded as such.


# ----------------------------------------------------------------- utilities
class StageTimeout(Exception):
    pass


def _alarm(signum, frame):
    raise StageTimeout()


def with_timeout(seconds, fn, *args, **kwargs):
    """Best-effort per-point wall-clock cap (SIGALRM, main thread only); the
    shell `timeout` in the AWS wrapper is the hard backstop."""
    seconds = int(seconds)
    if seconds <= 0:
        raise StageTimeout()
    old = signal.signal(signal.SIGALRM, _alarm)
    signal.alarm(seconds)
    try:
        return fn(*args, **kwargs)
    finally:
        signal.alarm(0)
        signal.signal(signal.SIGALRM, old)


def _u(v):
    v = np.asarray(v, dtype=float)
    return v / np.linalg.norm(v)


def jdump(obj, path):
    tmp = path + ".tmp"
    with open(tmp, "w") as f:
        json.dump(obj, f, ensure_ascii=False, indent=1)
    os.replace(tmp, path)


def jload(path, default=None):
    if os.path.exists(path):
        try:
            with open(path) as f:
                return json.load(f)
        except Exception:
            return default
    return default


# ----------------------------------------------------------------- geometry
def frame_atoms():
    """The Bi-O-Mo frame around the reactive Mo=O (generated, not hand-typed).

    Returns the 9 frame atoms in a fixed order; they are appended LAST so the
    geomeTRIC $freeze block is simply 'the last 9 atoms'. Mo sits at (0,0,zMo)
    and the reactive oxo at the origin (see start_atoms) — the frame is placed
    relative to Mo with normalized directions scaled to the target bond lengths.
    """
    mo = np.array([0.0, 0.0, -D_MO_OXO])
    u_oh = _u((-1.0, 0.0, -0.45))                 # Mo-OH, away from the oxo axis
    u_b1 = _u((0.85, 0.95, -0.55))                # Mo-(mu-O) 1
    u_b2 = _u((0.85, -0.95, -0.55))               # Mo-(mu-O) 2
    o_h = mo + D_MO_OH * u_oh
    h_oh = o_h + D_OH * _u(u_oh + np.array([0.0, 0.0, -0.5]))
    obr1 = mo + D_MO_OBR * u_b1
    obr2 = mo + D_MO_OBR * u_b2
    # Bi on the bisector of the two bridges, at |Bi-(mu-O)| = D_BI_OBR
    mid = 0.5 * (obr1 + obr2)
    half = float(np.linalg.norm(obr1 - mid))
    t2 = D_BI_OBR ** 2 - half ** 2
    if t2 <= 0:
        raise ValueError("Bi-(mu-O) length too short for the bridge geometry")
    bi = mid + math.sqrt(t2) * _u(mid - mo)
    w = _u(bi - mo)
    obi1 = bi + D_BI_OH * _u(w + np.array([0.0, 0.9, 0.0]))
    obi2 = bi + D_BI_OH * _u(w + np.array([0.0, -0.9, 0.0]))
    hbi1 = obi1 + D_OH * _u(obi1 - mo)
    hbi2 = obi2 + D_OH * _u(obi2 - mo)
    return [("O", o_h), ("H", h_oh), ("O", obr1), ("O", obr2), ("Bi", bi),
            ("O", obi1), ("H", hbi1), ("O", obi2), ("H", hbi2)]


def _propylene(zc, d=1.09):
    """CH2=CH-CH3 with the ABSTRACTED allylic H already removed (it is atom 2,
    placed on the HAT axis); C1 (the allylic carbon) sits at (0,0,zc) with the
    broken C-H bond pointing down -z toward the oxo."""
    th = math.radians(70.5)                       # 180 - 109.5 tetrahedral
    st, ct = math.sin(th), math.cos(th)
    c1 = np.array([0.0, 0.0, zc])
    out = [("C", c1)]
    c2 = c1 + 1.50 * np.array([st, 0.0, ct])
    out.append(("C", c2))
    u21 = _u(c1 - c2)
    # C3 in the xz-plane at C1-C2-C3 = 124 deg (sp2)
    ang = math.atan2(u21[2], u21[0]) - math.radians(124.0)
    u23 = np.array([math.cos(ang), 0.0, math.sin(ang)])
    c3 = c2 + 1.33 * u23
    out.append(("C", c3))
    for phi in (math.radians(120.0), math.radians(240.0)):   # 2 H left on C1
        out.append(("H", c1 + d * np.array([st * math.cos(phi),
                                            st * math.sin(phi), ct])))
    out.append(("H", c2 + 1.08 * _u(-(u21 + u23))))          # vinyl H on C2
    u32 = _u(c2 - c3)
    a32 = math.atan2(u32[2], u32[0])
    for da in (math.radians(120.0), math.radians(-120.0)):   # =CH2 on C3
        out.append(("H", c3 + 1.08 * np.array([math.cos(a32 + da), 0.0,
                                               math.sin(a32 + da)])))
    return out


def _acrolein(zc):
    """CH2=CH-CHO with the ABSTRACTED aldehydic H removed; C1 (the carbonyl
    carbon) at (0,0,zc), its broken C-H pointing down -z toward the oxo."""
    c1 = np.array([0.0, 0.0, zc])
    out = [("C", c1)]
    u_o = _u((math.sin(math.radians(120.0)), 0.0,
              math.cos(math.radians(120.0))))     # sp2, 120 deg from the C-H axis
    u_c2 = _u((-math.sin(math.radians(120.0)), 0.0,
               math.cos(math.radians(120.0))))
    o = c1 + 1.22 * u_o
    c2 = c1 + 1.47 * u_c2
    out += [("O", o), ("C", c2)]
    u21 = _u(c1 - c2)
    ang = math.atan2(u21[2], u21[0]) + math.radians(122.0)
    u23 = np.array([math.cos(ang), 0.0, math.sin(ang)])
    c3 = c2 + 1.33 * u23
    out.append(("C", c3))
    out.append(("H", c2 + 1.08 * _u(-(u21 + u23))))
    u32 = _u(c2 - c3)
    a32 = math.atan2(u32[2], u32[0])
    for da in (math.radians(120.0), math.radians(-120.0)):
        out.append(("H", c3 + 1.08 * np.array([math.cos(a32 + da), 0.0,
                                               math.sin(a32 + da)])))
    return out


def start_atoms(sub):
    """Full starting structure. FIXED order: Mo(0), oxo O(1), transferring H(2),
    C(3), rest of the substrate, then the 9 frozen frame atoms."""
    r_oh0, r_ch0 = R_POINT[1], 1.10
    atoms = [("Mo", np.array([0.0, 0.0, -D_MO_OXO])),
             ("O", np.array([0.0, 0.0, 0.0])),
             ("H", np.array([0.0, 0.0, r_oh0]))]
    zc = r_oh0 + r_ch0
    if sub == "c3h6":
        atoms += _propylene(zc)
    elif sub == "acrolein":
        atoms += _acrolein(zc)
    else:
        raise ValueError(f"unknown substrate {sub!r}")
    atoms += frame_atoms()
    return [(e, tuple(float(x) for x in p)) for e, p in atoms]


def geometry_report(sub):
    atoms = start_atoms(sub)
    xyz = np.array([p for _e, p in atoms])
    n = len(xyz)
    dmin, pair = np.inf, None
    for i in range(n):
        for j in range(i + 1, n):
            d = float(np.linalg.norm(xyz[i] - xyz[j]))
            if d < dmin:
                dmin, pair = d, (i, j)
    return {"substrate": sub, "n_atoms": n,
            "elements": "".join(e for e, _p in atoms),
            "n_frame_frozen": len(frame_atoms()),
            "min_pair_dist": round(dmin, 3), "min_pair": list(pair),
            "overlap_ok": bool(dmin > 0.85),
            "Mo_Bi_dist": round(float(np.linalg.norm(xyz[0] - xyz[n - 5])), 3)}


# ---------------------------------------------------------------- mean field
def _parity(atoms):
    return int(sum(gto.charge(e) for e, _p in atoms) % 2)


def _set_spin(s2):
    global SPIN, CAS_TARGET
    SPIN = int(s2)
    CAS_TARGET = (2 * N_DOCC + SPIN, N_DOCC + SPIN + N_VIR)


def resolve_spin(sub="c3h6"):
    """Fix the working 2S: env override > measured ladder > parity default."""
    if SPIN is not None and CAS_TARGET is not None:
        return SPIN
    if SPIN is not None:
        _set_spin(SPIN)
        return SPIN
    lad = jload(os.path.join(DIR, f"{PREF}_spins.json"), {}) or {}
    if lad.get("ground_2S") is not None:
        _set_spin(lad["ground_2S"])
        print(f"[spin] using measured ground 2S={SPIN}", flush=True)
    else:
        _set_spin(_parity(start_atoms(sub)))
        print(f"[spin] no ladder on disk — parity default 2S={SPIN}", flush=True)
    return SPIN


def build_mol(atoms):
    # ECP only for the elements that HAVE one in def2-SVP (Mo: ECP-28,
    # Bi: ECP-60); passing ecp=BASIS wholesale makes pyscf print an
    # "ECP not found" line per light element on every mol build.
    ecp = {e: BASIS for e, _p in atoms if e in ("Mo", "Bi")}
    return gto.M(atom=[(e, tuple(p)) for e, p in atoms], basis=BASIS, ecp=ecp,
                 charge=0, spin=SPIN, verbose=0, max_memory=MAXMEM)


def uks(mol, dm0=None):
    """UKS-PBE0 (DF), Newton-primary (stable on d0 Mo=O), DIIS fallback."""
    mf = dft.UKS(mol).density_fit()
    mf.xc = XC
    mf.conv_tol = 1e-8
    mfn = scf.newton(mf)
    mfn.kernel(dm0=dm0) if dm0 is not None else mfn.kernel()
    if not mfn.converged:
        mf2 = dft.UKS(mol).density_fit()
        mf2.xc = XC
        mf2.level_shift = 0.3
        mf2.max_cycle = 300
        mf2.kernel(dm0=dm0)
        mfn = scf.newton(mf2)
        mfn.kernel(mf2.make_rdm1())
    return mfn


def uks_stable(mol, dm0=None, max_hops=3):
    """UKS + internal stability: descend to the lower UKS solution if unstable
    (cures state hopping along the scan — the OCM `polish` lesson)."""
    mf = uks(mol, dm0=dm0)
    for _ in range(max_hops):
        mo_new = mf.stability()[0]
        if mo_new is mf.mo_coeff or (
                isinstance(mo_new, (tuple, list, np.ndarray))
                and np.allclose(np.asarray(mo_new[0]),
                                np.asarray(mf.mo_coeff[0]))):
            break
        mf = uks(mol, dm0=mf.make_rdm1(mo_new, mf.mo_occ))
    return mf


# --------------------------------------------------------------- stage: spins
def stage_spins():
    """UKS-PBE0 ladder over parity-correct 2S at the R geometry: the working spin
    state is measured before the pipeline, not postulated. Resume-aware."""
    path = os.path.join(DIR, f"{PREF}_spins.json")
    done = jload(path, {}) or {}
    if done.get("ground_2S") is not None:
        print(f"[spins] resume: ground 2S={done['ground_2S']} on disk "
              f"(ladder {done.get('rel_kcal')}) — skip", flush=True)
        _set_spin(done["ground_2S"])
        return
    atoms = start_atoms("c3h6")
    p = _parity(atoms)
    cands = [p + 2 * k for k in range(3)]
    out = {"candidates_2S": cands, "geometry": geometry_report("c3h6"),
           "ladder": []}
    for s2 in cands:
        t0 = time.time()
        _set_spin(s2)
        try:
            mf = with_timeout(STAGE_TIMEOUT, uks, build_mol(atoms))
            out["ladder"].append({"spin_2S": s2, "e_h": float(mf.e_tot),
                                  "converged": bool(mf.converged),
                                  "spin_square": round(float(
                                      mf.spin_square()[0]), 3),
                                  "wall_s": round(time.time() - t0, 1)})
        except StageTimeout:
            out["ladder"].append({"spin_2S": s2, "timed_out_after_s":
                                  STAGE_TIMEOUT})
        except Exception as e:
            out["ladder"].append({"spin_2S": s2, "error": str(e)[:200]})
        jdump(out, path)
        print(f"[spins] {out['ladder'][-1]}", flush=True)
    ok = [x for x in out["ladder"] if x.get("converged")]
    if ok:
        best = min(ok, key=lambda x: x["e_h"])
        out["ground_2S"] = best["spin_2S"]
        out["rel_kcal"] = {str(x["spin_2S"]):
                           round((x["e_h"] - best["e_h"]) * HARTREE_KCAL, 2)
                           for x in ok}
        _set_spin(best["spin_2S"])
        print(f"[spins] ground 2S={SPIN}; ladder {out['rel_kcal']}", flush=True)
    jdump(out, path)


# ---------------------------------------------------------------- stage: geom
def constrained_opt(atoms, rCH, rOH, tag, dm0=None, maxsteps=110):
    """Relaxed optimization with the r(C-H) / r(O-H) pins (geomeTRIC $set) and
    the Bi-O frame frozen ($freeze on the last len(frame) atoms). Warm-started
    from the previous point's density: faster AND keeps one electronic branch
    along the scan."""
    from pyscf.geomopt.geometric_solver import optimize
    cfile = os.path.join(DIR, f"_constr_{PREF}_{tag}.txt")
    lines = ["$set"]
    if rOH is not None:
        lines.append(f"distance 2 3 {rOH:.4f}")     # O(idx1)-H(idx2), 1-based
    if rCH is not None:
        lines.append(f"distance 3 4 {rCH:.4f}")     # H(idx2)-C(idx3), 1-based
    n, n_env = len(atoms), len(frame_atoms())
    lines += ["$freeze", f"xyz {n - n_env + 1}-{n}"]
    with open(cfile, "w") as f:
        f.write("\n".join(lines) + "\n")
    try:
        mf = uks(build_mol(atoms), dm0=dm0)
        conv = dict(convergence_energy=2e-6, convergence_grms=4.5e-4,
                    convergence_gmax=9e-4, convergence_drms=2.3e-3,
                    convergence_dmax=3.6e-3)
        mol_eq = optimize(mf, constraints=cfile, maxsteps=maxsteps,
                          assert_convergence=False, **conv)
    finally:
        if os.path.exists(cfile):
            os.remove(cfile)
    new_atoms = [(mol_eq.atom_symbol(i), tuple(c))
                 for i, c in enumerate(mol_eq.atom_coords(unit="Angstrom"))]
    mf_final = uks(mol_eq, dm0=mf.make_rdm1())
    return (new_atoms, float(mf_final.e_tot), bool(mf_final.converged),
            float(mf_final.spin_square()[0]), mf_final.make_rdm1())


def stage_geom(sub):
    resolve_spin(sub)
    path = os.path.join(DIR, f"{PREF}_{sub}_geom.json")
    out = {"substrate": sub, "geom_version": GEOM_VERSION, "spin_2S": SPIN,
           "geometry": geometry_report(sub),
           "protocol": f"UKS-{XC.upper()}/{BASIS}, 2D pin r(C-H)+r(O-H) "
                       "(geomeTRIC), Bi-O frame frozen",
           "points": []}
    atoms = start_atoms(sub)
    grid = [R_POINT] + PATH
    saved = jload(path)
    if saved and saved.get("geom_version") == GEOM_VERSION:
        out["points"] = saved.get("points", [])[:len(grid)]
        if out["points"]:
            atoms = [(s, tuple(c)) for s, c in out["points"][-1]["atoms"]]
            print(f"[{sub}] resume: {len(out['points'])}/{len(grid)} points",
                  flush=True)
    elif saved:
        print(f"[{sub}] geom_version mismatch — recomputing all points",
              flush=True)
    dm = None
    for i, (rCH, rOH) in enumerate(grid):
        if i < len(out["points"]):
            continue
        t0 = time.time()
        try:
            atoms, e, conv, ss, dm = with_timeout(
                STAGE_TIMEOUT, constrained_opt, atoms, rCH, rOH, f"{sub}_{i}",
                dm, 80 if i == 0 else 110)
        except StageTimeout:
            out["truncated_at"] = i
            out["truncation_reason"] = f"point {i} exceeded {STAGE_TIMEOUT}s"
            print(f"[{sub}] point {i} TIMED OUT — profile truncated", flush=True)
            break
        out["points"].append({
            "rCH_pin": rCH, "rOH_pin": rOH, "e_h": e, "scf_converged": conv,
            "spin_square": round(ss, 3), "wall_s": round(time.time() - t0, 1),
            "atoms": [[s, [round(float(x), 6) for x in c]] for s, c in atoms]})
        jdump(out, path)
        print(f"[{sub}] point {i}/{len(grid)-1} rCH={rCH} rOH={rOH} "
              f"E={e:.6f} conv={conv} <S2>={ss:.2f} "
              f"({out['points'][-1]['wall_s']}s)", flush=True)
    _finalize_dft(out, path, sub)


def _finalize_dft(out, path, sub):
    es = [p["e_h"] for p in out["points"]]
    if len(es) < 3:
        out["dft_readable"] = False
        jdump(out, path)
        print(f"[{sub}] fewer than 3 points — DFT profile unreadable", flush=True)
        return
    rel = [(e - es[0]) * HARTREE_KCAL for e in es]
    out["rel_kcal"] = [round(x, 2) for x in rel]
    r, ts, bar, interior, dropped = robust_read(rel)
    out.update({"r_index": r, "ts_index": ts, "ts_interior": interior,
                "dropped_points": dropped, "dft_barrier_kcal": bar,
                "dft_readable": bool(bar is not None and interior
                                     and BARRIER_MIN < bar < BARRIER_MAX)})
    jdump(out, path)
    print(f"[{sub}] DFT profile {out['rel_kcal']} -> barrier {bar} "
          f"(R={r} TS={ts} interior={interior} dropped={dropped})", flush=True)


def stage_polish(sub):
    """Warm-started chain + UKS stability over the saved geometries, forward and
    backward; keeps the lower solution. Removes the cold-start sawtooth."""
    resolve_spin(sub)
    path = os.path.join(DIR, f"{PREF}_{sub}_geom.json")
    g = jload(path)
    if not g or not g.get("points"):
        print(f"[{sub}] nothing to polish", flush=True)
        return
    for direction in ("fwd", "rev"):
        dm = None
        idxs = range(len(g["points"])) if direction == "fwd" \
            else range(len(g["points"]) - 1, -1, -1)
        for i in idxs:
            p = g["points"][i]
            t0 = time.time()
            try:
                mf = with_timeout(STAGE_TIMEOUT, uks_stable,
                                  build_mol([(s, tuple(c)) for s, c
                                             in p["atoms"]]), dm)
            except StageTimeout:
                print(f"[{sub}] polish-{direction} {i}: timed out — kept",
                      flush=True)
                continue
            dm = mf.make_rdm1()
            p.setdefault("e_h_raw", p["e_h"])
            if float(mf.e_tot) < p["e_h"] - 1e-6 or direction == "fwd":
                if float(mf.e_tot) < p["e_h"] - 1e-6:
                    p["e_h"] = float(mf.e_tot)
                    p["spin_square"] = round(float(mf.spin_square()[0]), 3)
                    p["scf_converged"] = bool(mf.converged)
                p["polished"] = True
            jdump(g, path)
            print(f"[{sub}] polish-{direction} {i}: {p['e_h']:.6f} "
                  f"(raw {p['e_h_raw']:.6f}) ({time.time()-t0:.0f}s)", flush=True)
    g["polish_done"] = True
    _finalize_dft(g, path, sub)


# ------------------------------------------------------------ correlated part
def rohf(mol):
    mf = scf.ROHF(mol)
    mf.conv_tol = 1e-8
    mfn = scf.newton(mf)
    mfn.kernel()
    if not mfn.converged:
        mf2 = scf.ROHF(mol)
        mf2.level_shift = 0.4
        mf2.max_cycle = 300
        mf2.kernel()
        mfn = scf.newton(mf2)
        mfn.kernel(mf2.make_rdm1())
    return mfn


def _proj_metric(mol):
    """AO-projection metric onto the AVAS reference labels. MINAO first; if the
    metal labels are absent there (4d/6p elements behind an ECP are not always in
    MINAO), fall back to the molecule's own basis. Recorded in the JSON."""
    for pb in ("minao", None):
        try:
            if pb is None:
                pmol = mol
            else:
                pmol = mol.copy()
                pmol.basis = pb
                pmol.build(False, False)
            baslst = pmol.search_ao_label(AVAS_LABELS)
        except Exception:
            continue
        if len(baslst) == 0:
            continue
        labels = pmol.ao_labels()
        if not any(" Mo " in labels[i] for i in baslst):
            continue
        s2 = pmol.intor_symmetric("int1e_ovlp")[np.ix_(baslst, baslst)]
        s21 = gto.intor_cross("int1e_ovlp", pmol, mol)[baslst]
        return s21.T @ np.linalg.solve(s2, s21), (pb or "own-basis"), len(baslst)
    raise RuntimeError("AVAS reference labels not found in minao or own basis")


def forced_avas(mf):
    """Fixed AVAS: an active space of FIXED COMPOSITION at any geometry — all
    SOMOs + the top N_DOCC closed pairs + the top N_VIR virtuals by projection on
    the reference AOs (Mo valence 4d, O 2p, transferring-H 1s). Guarantees the
    same (ne, no) by construction; threshold-picked spaces jump size between R
    and TS (the OCM Stage-14 failure) and make barriers incomparable."""
    mol = mf.mol
    sa, src, nref = _proj_metric(mol)
    occ = np.asarray(mf.mo_occ)
    C = np.asarray(mf.mo_coeff)
    idx_d, idx_s, idx_v = (np.where(occ == 2)[0], np.where(occ == 1)[0],
                           np.where(occ == 0)[0])

    def split(idx, n_keep):
        Cb = C[:, idx]
        w, v = np.linalg.eigh(Cb.T @ sa @ Cb)
        order = np.argsort(w)[::-1]
        Cr = Cb @ v[:, order]
        return Cr[:, :n_keep], Cr[:, n_keep:], w[order][:n_keep]

    act_d, core_d, w_d = split(idx_d, N_DOCC)
    act_v, rest_v, w_v = split(idx_v, N_VIR)
    mo = np.hstack([core_d, act_d, C[:, idx_s], act_v, rest_v])
    ne = 2 * N_DOCC + len(idx_s)
    no = N_DOCC + len(idx_s) + N_VIR
    if (ne, no) != CAS_TARGET:
        raise RuntimeError(f"AVAS produced CAS({ne}e,{no}o) != {CAS_TARGET}")
    return mo, [round(float(x), 3) for x in list(w_d) + list(w_v)], src, nref


def cas_nevpt2(atoms, tag):
    """ROHF -> fixed AVAS -> CASCI natural orbitals -> CASSCF -> SC-NEVPT2."""
    mol = build_mol(atoms)
    mf = rohf(mol)
    mo, wsel, src, nref = forced_avas(mf)
    ne, no = CAS_TARGET
    ss_target = SPIN / 2.0 * (SPIN / 2.0 + 1.0)
    mc0 = mcscf.CASCI(mf, no, ne)
    mc0.fcisolver = fci.direct_spin1.FCI(mol)
    fci.addons.fix_spin_(mc0.fcisolver, shift=0.2, ss=ss_target)
    mc0.kernel(mo)
    mo_nat = mc0.cas_natorb()[0]
    mc = mcscf.CASSCF(mf, no, ne)
    fci.addons.fix_spin_(mc.fcisolver, shift=0.2, ss=ss_target)
    mc.max_cycle_macro = 150
    mc.conv_tol = 1e-7
    mc.kernel(mo_nat)
    dm1 = mc.fcisolver.make_rdm1(mc.ci, mc.ncas, mc.nelecas)
    noon = sorted(np.linalg.eigvalsh(dm1), reverse=True)
    out = {"tag": tag, "cas": list(CAS_TARGET), "avas_sel_weights": wsel,
           "avas_reference_basis": src, "avas_n_reference_aos": int(nref),
           "e_rohf_h": float(mf.e_tot), "rohf_converged": bool(mf.converged),
           "e_casci_h": float(mc0.e_tot),
           "e_casscf_h": float(mc.e_tot),
           "casscf_converged": bool(mc.converged),
           "noon": [round(float(n), 3) for n in noon],
           "n_unpaired_beyond_spin": round(float(
               sum(min(n, 2 - n) for n in noon) - SPIN), 3)}
    if LEVEL in ("all", "nevpt2"):
        e_pt = NEVPT(mc).kernel()
        out["nevpt2_corr_h"] = float(e_pt)
        out["e_nevpt2_h"] = float(mc.e_tot + e_pt)
    return out


def stage_profile(sub):
    """CASSCF (+NEVPT2) at EVERY scan point: the correlated surface picks its own
    R and TS instead of inheriting the DFT indices."""
    if LEVEL == "dft":
        print(f"[{sub}] AMMOX_LEVEL=dft — correlated profile skipped", flush=True)
        return
    resolve_spin(sub)
    g = jload(os.path.join(DIR, f"{PREF}_{sub}_geom.json"))
    if not g or not g.get("points"):
        print(f"[{sub}] no geometries — run `geom` first", flush=True)
        return
    path = os.path.join(DIR, f"{PREF}_{sub}_profile.json")
    out = jload(path) or {"substrate": sub, "spin_2S": SPIN, "points": []}
    out["cas"] = list(CAS_TARGET)
    out["avas_labels"] = AVAS_LABELS
    for i, p in enumerate(g["points"]):
        if i < len(out["points"]):
            continue
        atoms = [(s, tuple(c)) for s, c in p["atoms"]]
        t0 = time.time()
        try:
            res = with_timeout(STAGE_TIMEOUT, cas_nevpt2, atoms, f"{sub}_p{i}")
        except StageTimeout:
            out["truncated_at"] = i
            out["truncation_reason"] = f"point {i} exceeded {STAGE_TIMEOUT}s"
            jdump(out, path)
            print(f"[{sub}] correlated point {i} TIMED OUT — truncated",
                  flush=True)
            break
        except Exception as e:
            out["truncated_at"] = i
            out["truncation_reason"] = f"{type(e).__name__}: {e}"
            jdump(out, path)
            print(f"[{sub}] correlated point {i} FAILED: {e}", flush=True)
            break
        res["wall_s"] = round(time.time() - t0, 1)
        out["points"].append(res)
        jdump(out, path)
        print(f"[{sub}] profile {i}: CASSCF {res['e_casscf_h']:.6f} "
              f"(conv={res['casscf_converged']}) "
              f"NEVPT2 {res.get('e_nevpt2_h', float('nan')):.6f} "
              f"({res['wall_s']}s)", flush=True)
    levels = ["casscf"] + (["nevpt2"] if LEVEL in ("all", "nevpt2") else [])
    for lvl in levels:
        key = f"e_{lvl}_h"
        es = [p[key] for p in out["points"] if key in p]
        if len(es) < 3:
            out[f"{lvl}_readable"] = False
            continue
        rel = [(e - es[0]) * HARTREE_KCAL for e in es]
        r, ts, bar, interior, dropped = robust_read(rel)
        out[f"{lvl}_rel_kcal"] = [round(x, 2) for x in rel]
        out[f"{lvl}_r_index"], out[f"{lvl}_ts_index"] = r, ts
        out[f"{lvl}_ts_interior"] = interior
        out[f"{lvl}_dropped_points"] = dropped
        out[f"{lvl}_barrier_kcal"] = bar
        out[f"{lvl}_readable"] = bool(bar is not None and interior
                                      and BARRIER_MIN < bar < BARRIER_MAX)
    out["all_casscf_converged"] = all(p.get("casscf_converged")
                                      for p in out["points"])
    jdump(out, path)
    print(f"[{sub}] correlated: CASSCF {out.get('casscf_barrier_kcal')} | "
          f"NEVPT2 {out.get('nevpt2_barrier_kcal')} kcal/mol", flush=True)


# ------------------------------------------------------- mirage-detector gates
def robust_read(rel):
    """The Stage-18 robust readout, verbatim in spirit:
    a point more than OUTLIER_KCAL above BOTH surviving neighbours is masked; on
    what is left, TS = the interior maximum, R = the minimum to its left.
    Returns (r_idx, ts_idx, barrier, interior_flag, dropped_indices) in ORIGINAL
    indexing; (None, None, None, False, dropped) if the profile is unreadable."""
    n = len(rel)
    dropped = []
    for i in range(n):
        nb = [rel[j] for j in (i - 1, i + 1) if 0 <= j < n and j not in dropped]
        if nb and all(rel[i] - x > OUTLIER_KCAL for x in nb):
            dropped.append(i)
    keep = [i for i in range(n) if i not in dropped]
    if len(keep) < 4:
        return None, None, None, False, dropped
    kr = [rel[i] for i in keep]
    imax_k = max(range(1, len(kr) - 1), key=lambda i: kr[i])
    interior = bool(kr[imax_k] >= max(kr[0], kr[-1]))
    imin_k = min(range(imax_k + 1), key=lambda i: kr[i])
    return (keep[imin_k], keep[imax_k], round(kr[imax_k] - kr[imin_k], 2),
            interior, dropped)


def gate_report(block, lvl):
    """Explicit pass/fail of every gate for one substrate at one level."""
    bar = block.get(f"{lvl}_barrier_kcal")
    interior = block.get(f"{lvl}_ts_interior")
    dropped = block.get(f"{lvl}_dropped_points", [])
    rel = block.get(f"{lvl}_rel_kcal") or block.get("rel_kcal") or []
    gates = {
        "profile_has_points": bool(len(rel) >= 4),
        "ts_interior": bool(interior),
        "barrier_in_sanity_window": bool(bar is not None
                                         and BARRIER_MIN < bar < BARRIER_MAX),
        "not_zigzag": bool(len(dropped) <= 1)}
    return {"barrier_kcal": bar, "rel_kcal": rel, "dropped_points": dropped,
            "r_index": block.get(f"{lvl}_r_index"),
            "ts_index": block.get(f"{lvl}_ts_index"),
            "gates": gates, "passed": bool(all(gates.values()))}


def stage_merge():
    """Assemble the deliverable: profiles, barriers at three levels, the gates,
    and ddE# WITH ITS SIGN — or an explicit NO-VERDICT."""
    resolve_spin()
    res = {
        "purpose": "L7: +1% выхода акрилонитрила на ДЕЙСТВУЮЩЕМ производстве "
                   "(Саратоворгсинтез, Sohio-процесс на Mo-Bi-оксидах) — перенос "
                   "нашего OCM-дескриптора селективности 1:1 на аммоксидирование "
                   "пропилена",
        "cluster": "(HO)(O=)Mo(mu-O)2Bi(OH)2 — нейтральный Bi-O-Mo мостик, "
                   "минимальный прокси активного центра Sohio; заряд 0, "
                   "Mo(VI) d0 + Bi(III)",
        "bond_lengths_A": {"Mo=O(oxo)": D_MO_OXO, "Mo-OH": D_MO_OH,
                           "Mo-(mu-O)": D_MO_OBR, "Bi-(mu-O)": D_BI_OBR,
                           "Bi-OH": D_BI_OH},
        "atom_order": "0 Mo, 1 O(oxo), 2 H(transferring), 3 C, ... substrate, "
                      "last 9 = frozen Bi-O frame",
        "spin_2S": SPIN, "basis": BASIS, "xc_geometry": XC,
        "cas": list(CAS_TARGET), "avas_labels": AVAS_LABELS,
        "avas_note": "«Mo 3d» в постановке трека = ВАЛЕНТНАЯ d-оболочка металла; "
                     "для Mo это 4d (3d сидит внутри ECP-28 базиса def2-SVP)",
        "descriptor": "ddE# = barrier(доокисление: C-H акролеина) - "
                      "barrier(целевой аллильный H-отрыв от пропилена); "
                      "> 0 — центр ЩАДИТ целевой путь => селективность к НАК; "
                      "< 0 — доокисляет интермедиат быстрее сырья (стена BEP)",
        "gates": {"outlier_mask_kcal": OUTLIER_KCAL,
                  "ts_must_be_interior": "rel[ts] >= max(rel[0], rel[-1])",
                  "barrier_sanity_window_kcal": [BARRIER_MIN, BARRIER_MAX],
                  "no_verdict_is_legal": True},
        "honesty": [
            "газофазный кластер из 19-20 атомов — НЕ периодический Mo-Bi-O "
            "катализатор: нет резервуара решёточного кислорода и стадии его "
            "восполнения (Марс–ван Кревелен), нет активации NH3, нет фазы "
            "alpha-Bi2Mo3O12, нет покрытия и температуры",
            "distinguished-coordinate релакс-скан даёт ВЕРХНЮЮ оценку барьера "
            "(истинное седло не выше максимума скана); частот и ZPE нет",
            "def2-SVP + ECP на Mo и Bi; активное пространство фиксировано "
            "конструкцией (не подбором порога) — размер одинаков во всех точках",
            "переносится ЗНАК и величина ddE#, а не абсолютные барьеры",
            "первая координационная оболочка Bi-O заморожена (embedded-"
            "приближение): каркас идеализированный, не кристаллографический"],
        "substrates": {}}

    levels = ["dft", "casscf"] + (["nevpt2"] if LEVEL in ("all", "nevpt2")
                                  else [])
    all_ok = True
    for sub in SUBSTRATES:
        g = jload(os.path.join(DIR, f"{PREF}_{sub}_geom.json"))
        pr = jload(os.path.join(DIR, f"{PREF}_{sub}_profile.json")) or {}
        if not g:
            res["substrates"][sub] = {"error": "no geom json"}
            all_ok = False
            continue
        sd = {"geometry": g.get("geometry"), "spin_2S": g.get("spin_2S"),
              "n_points": len(g.get("points", [])),
              "wall_s_per_point": [p.get("wall_s") for p in g.get("points", [])],
              "spin_square_profile": [p.get("spin_square")
                                      for p in g.get("points", [])]}
        dft_blk = {"dft_rel_kcal": g.get("rel_kcal"),
                   "dft_barrier_kcal": g.get("dft_barrier_kcal"),
                   "dft_ts_interior": g.get("ts_interior"),
                   "dft_dropped_points": g.get("dropped_points", []),
                   "dft_r_index": g.get("r_index"),
                   "dft_ts_index": g.get("ts_index")}
        sd["dft"] = gate_report(dft_blk, "dft")
        for lvl in levels[1:]:
            sd[lvl] = gate_report(pr, lvl)
        sd["all_casscf_converged"] = pr.get("all_casscf_converged")
        if pr.get("points"):
            ts = pr.get("nevpt2_ts_index", pr.get("casscf_ts_index"))
            if ts is not None and ts < len(pr["points"]):
                sd["noon_at_ts"] = pr["points"][ts].get("noon")
        for tk in ("truncated_at", "truncation_reason"):
            if g.get(tk) is not None:
                sd[f"geom_{tk}"] = g[tk]
            if pr.get(tk) is not None:
                sd[f"profile_{tk}"] = pr[tk]
        res["substrates"][sub] = sd

    top = levels[-1]
    ok_top = all(res["substrates"].get(s, {}).get(top, {}).get("passed")
                 for s in SUBSTRATES)
    res["ddE_kcal"] = {}
    for lvl in levels:
        b = {s: res["substrates"].get(s, {}).get(lvl, {}).get("barrier_kcal")
             for s in SUBSTRATES}
        passed = all(res["substrates"].get(s, {}).get(lvl, {}).get("passed")
                     for s in SUBSTRATES)
        res["ddE_kcal"][lvl] = (round(b["acrolein"] - b["c3h6"], 2)
                                if passed and None not in b.values() else None)
    if (res["ddE_kcal"].get(top) is not None
            and res["ddE_kcal"].get("dft") is not None):
        res["quantum_shift_kcal"] = round(res["ddE_kcal"][top]
                                          - res["ddE_kcal"]["dft"], 2)
    if ok_top and all_ok and res["ddE_kcal"].get(top) is not None:
        dde = res["ddE_kcal"][top]
        res["verdict"] = {
            "level": top, "ddE_kcal": dde, "sign": "+" if dde > 0 else "-",
            "reading": ("ddE# > 0: центр ЩАДИТ целевой аллильный путь — мотив "
                        "работает на селективность к акрилонитрилу"
                        if dde > 0 else
                        "ddE# < 0: центр доокисляет интермедиат быстрее, чем "
                        "отрывает целевой аллильный H — стена BEP на этом мотиве"),
            "caveat": "знак получен на газофазном прокси при верхней оценке "
                      "барьеров; переносится как ранжирующий, не абсолютный"}
    else:
        failed = {s: {lv: res["substrates"].get(s, {}).get(lv, {}).get("gates")
                      for lv in levels} for s in SUBSTRATES}
        res["verdict"] = {"level": top, "ddE_kcal": None,
                          "reading": "NO VERDICT: гейты мираж-детектора не "
                                     "пройдены (край-TS / зигзаг / барьер вне "
                                     "диапазона / профиль оборван) — дескриптор "
                                     "НЕ публикуется",
                          "gate_detail": failed}
    path = os.path.join(DIR, f"{PREF}_results.json")
    jdump(res, path)
    print(json.dumps({k: v for k, v in res.items()
                      if k in ("ddE_kcal", "quantum_shift_kcal", "verdict")},
                     ensure_ascii=False, indent=1))
    print("saved", path, flush=True)


# ------------------------------------------------------------------ driver
def stage_all():
    stage_spins()
    for sub in SUBSTRATES:
        stage_geom(sub)
        stage_polish(sub)
    if LEVEL != "dft":
        for sub in SUBSTRATES:
            stage_profile(sub)
    stage_merge()


if __name__ == "__main__":
    stage = sys.argv[1] if len(sys.argv) > 1 else "all"
    if stage == "all":
        stage_all()
    elif stage == "spins":
        stage_spins()
    elif stage == "geom":
        stage_geom(sys.argv[2])
    elif stage == "polish":
        stage_polish(sys.argv[2])
    elif stage == "profile":
        stage_profile(sys.argv[2])
    elif stage == "merge":
        stage_merge()
    elif stage == "geometry":                 # cheap self-check, no QM
        for s in SUBSTRATES:
            print(json.dumps(geometry_report(s), indent=1))
    else:
        raise SystemExit("usage: all | spins | geometry | "
                         "geom|polish|profile c3h6|acrolein | merge")
