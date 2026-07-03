#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
routes/polyolefins_chiral_aws.py — polyolefins track (GitHub issue #11):
a CHIRAL catalyst model so that the stereo-descriptor DDE++ = E++(si) - E++(re)
is no longer zero BY SYMMETRY.

WHY THIS JOB EXISTS (diary record 2026-06-08, tab "polyolefins"):
on the minimal ACHIRAL model [Cl2Ti-CH3]+ the re and si monomer faces gave
EXACTLY the same energy (CASCI and DFT: DE(si-re)=0.00 kcal/mol) — they are
enantiomeric, symmetry-degenerate transition states. That was the honest
negative result: stereoselectivity is not a property of the bare metal, it
requires a CHIRAL ligand pocket. This job builds the smallest such pocket.

MODEL CHOICE (documented, with its limits):
  [H2Si(eta5-3-Me-C5H3)2 Ti-CH3]+  +  propene   (41 atoms, 154 e-, charge +1,
                                                 closed-shell d0 Ti(IV))
  i.e. a rac-C2-symmetric ansa-titanocene mimic:
    - two cyclopentadienyl rings locked by an H2Si ansa bridge (the bridge is
      what prevents ring rotation from washing out the chirality);
    - one 3-methyl on EACH ring, the two methyls related by the C2 axis
      (Ti...Si) — the rac arrangement. This is the minimal stand-in for
      rac-Et(Ind)2ZrCl2-type isospecific catalysts: the benzo ring of indenyl
      is truncated to a single methyl, SiMe2 bridge to SiH2.
    - Ti (not Zr): keeps comparability with the achiral [Cl2Ti-CH3]+ diary
      stages and avoids switching to an ECP mid-track.
  LIMITS (honest): a 3-methyl is a much weaker stereo-differentiator than a
  fused benzo ring, so |DDE++| here is expected to be SMALL and may drown in
  method noise — that outcome is detected and flagged, not hidden (see the
  degeneracy gate below). The ligand geometry is idealized (parallel planar
  rings, in-plane ring H) and FROZEN — see below. One ligand enantiomer is
  fixed (ENANT=+1); the other gives exactly -DDE++ by mirror symmetry.

WHAT IS COMPUTED:
  1. Two diastereomeric insertion paths — propene coordinating with its si or
     its re face to the SAME chiral pocket — x 3 points each along the
     Cossee-Arlman 1,2-insertion coordinate cc = d(C_methyl...C2(propene)),
     the forming C-C bond (same coordinate as the achiral diary stages):
        cc = 2.8 A (pi-complex-like), 2.0 A (4-center TS approximation),
        1.6 A (product-like).
     Each point: constrained relaxation (geomeTRIC through PySCF) at
     RKS-PBE/def2-SVP with density fitting (POLY_MF=UKS switches to UKS,
     same spin=0), with
       - the forming C-C distance FIXED at cc ($set distance), and
       - the WHOLE ansa ligand FROZEN ($freeze xyz): only Ti, the growing
         chain CH3 and propene relax. Rationale: the frozen idealized pocket
         removes ligand-relaxation noise and keeps the si/re pockets strictly
         identical (differential errors partly cancel); the price is that
         ABSOLUTE barriers are biased and NOT publishable as such — recorded
         in the JSON as a caveat. cc=2.0 is a constrained-scan TS
         APPROXIMATION, not an optimized saddle point (no frequencies).
  2. Stereo-descriptor:  DDE_TS  = E_TS(si) - E_TS(re)   (Curtin-Hammett:
     with fast face exchange the TS energies alone decide), plus the barrier
     difference DDE++ = [E_TS-E_pi](si) - [E_TS-E_pi](re).
  3. CASCI(AVAS: Ti 3d + 2p of the three reacting carbons) on BOTH TS
     geometries -> NOON, N_u = sum n(2-n) (same metric as the diary's 0.97 on
     the achiral TS), and a multireference correction to the descriptor:
        mr_correction = DDE_TS(CASCI) - DDE_TS(DFT).
     The AVAS space is truncated to a tractable CAS (default 6e,6o — same
     size the diary used) by Loewdin character on the target AOs; the same
     (ne,no) is ENFORCED on both faces (a differing space voids the
     comparison and is flagged). CASCI is on a DF-RHF reference (as in the
     achiral stages); it captures static, not dynamic, correlation — the
     correction is a diagnostic, not a final number. All of this is in the
     output JSON.

HONESTY GATES (all recorded in results["stereo"]):
  - degenerate_within_noise: |DDE_TS(DFT)| < POLY_DEGEN_TOL (0.05 kcal/mol)
    -> "the chiral pocket did not resolve the faces" — the main model risk
    (methyl too weak), reported as a finding, not spun.
  - suspiciously_exact_zero: |DDE_TS| < 1e-4 kcal/mol -> smells like an
    accidental mirror symmetry in the construction (a bug), flagged apart
    from the physical near-degeneracy. The smoke test additionally verifies
    that the ligand is NOT superimposable on its own mirror images.
  - face_preserved: the propene face is re-identified AFTER each relaxation;
    if the olefin flipped faces during optimization the two "diastereomers"
    collapsed onto one path and DDE is meaningless (flagged).
  - convergence flags of every SCF / optimization / CASCI, and
    |d(cc) - cc| of the honored constraint.

ENV KNOBS (all optional):
  POLY_STAGE_TIMEOUT   seconds per stage (one constrained opt / one CASCI),
                       default 4800; best-effort SIGALRM — the shell-level
                       `timeout` in calc/aws_run_polyolefins.sh is the hard stop
  POLY_TOTAL_BUDGET    total wall-clock budget (s), default 19800
  POLY_CC_POINTS       ladder of scan points, default "2.8,2.0,1.6";
                       points are processed in the given order, si and re
                       interleaved per point, and remaining points are
                       SKIPPED (recorded) when the budget runs short —
                       so the crucial 2.8+2.0 pairs are done first
  POLY_TS_CC           which point is "the TS" for DDE/CASCI, default 2.0
  POLY_FACES           "si,re" (default) or one of them
  POLY_MF              "RKS" (default) | "UKS"      POLY_XC  default "pbe"
  POLY_OPT_MAXSTEPS    geomeTRIC max steps per point, default 60
  POLY_CAS             truncated CASCI space "ne,no", default "6,6"
  POLY_AVAS_THRESHOLD  AVAS threshold, default 0.2
  POLY_DEGEN_TOL       degeneracy gate, kcal/mol, default 0.05

OUTPUT (incremental, atomic — a kill loses nothing):
  routes/polyolefins_chiral_results.json   rolling results, restartable:
                                           finished (face, cc) points and
                                           finished CASCI blocks are skipped
                                           on rerun (geometries are stored in
                                           the JSON itself)
  routes/polyolefins_chiral_model.xyz      initial template (TS point, face 1)
  routes/polyolefins_chiral_ts_<face>.xyz  optimized TS-point geometries

RUN
  heavy (meant for AWS via calc/aws_run_polyolefins.sh):
      python3 -u routes/polyolefins_chiral_aws.py
  smoke (<= 3 min, 1 core: geometry assembly for all faces x points, xyz
  dumps, C2/chirality/face/constraint validation, ast self-parse, env-knob
  validation, geomeTRIC import gate, plus ONE tiny SCF+AVAS+CASCI on free
  propene to exercise the CASCI machinery; NO big SCF):
      OMP_NUM_THREADS=1 python3 routes/polyolefins_chiral_aws.py --smoke
      (--smoke --no-scf skips even the tiny SCF)
"""
import os
import sys
import ast
import json
import time
import signal
import argparse
import tempfile
import math

import numpy as np

try:
    from pyscf import gto, scf, dft, mcscf, fci, lib
    from pyscf.mcscf import avas
except ImportError as e:
    sys.exit(f"[deps] pyscf is required:  python3 -m pip install pyscf   ({e})")

HERE = os.path.dirname(os.path.abspath(__file__))
RESULTS = os.path.join(HERE, "polyolefins_chiral_results.json")
XYZ_MODEL = os.path.join(HERE, "polyolefins_chiral_model.xyz")
XYZ_TS_TPL = os.path.join(HERE, "polyolefins_chiral_ts_{face}.xyz")
KCAL = 627.50947406

BASIS = "def2-svp"
CHARGE = 1                # [L2Ti-CH3 + propene]+ , d0 Ti(IV), closed shell
SPIN = 0

# ---- idealized construction parameters (A / deg; documented, frozen) --------
TI_CEN = 2.08             # Ti - Cp centroid
ALPHA_DEG = 65.0          # centroid axis vs C2(z) axis -> Cp-Ti-Cp = 130 deg
R_RING = 1.42 / (2 * math.sin(math.pi / 5))   # Cp circumradius (C-C 1.42)
CH_RING = 1.08            # ring C-H, in-plane radial (idealized)
CC_ME = 1.50              # ring C - methyl C
CH_ME = 1.09              # methyl C-H
SI_C = 1.87               # bridge Si - ipso C
SI_H = 1.48
TI_CM = 2.10              # Ti - chain methyl C
TI_C1 = 2.30              # Ti - olefin C1 (site template)
SIGMA_DEG = 47.5          # sigma sites: +-SIGMA_DEG from -z, in the yz plane
CC_DB = 1.40              # propene C1=C2 (template)
ENANT = +1                # fixed ligand enantiomer; the mirror one gives -DDE

T0 = time.time()


# ------------------------------------------------------------------ utilities
def elapsed():
    return time.time() - T0


def say(msg):
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def atomic_save(obj, path):
    """Atomic incremental save: tmp file in the same dir + os.replace."""
    fd, tmp = tempfile.mkstemp(dir=os.path.dirname(path) or ".",
                               prefix=os.path.basename(path) + ".", suffix=".tmp")
    try:
        with os.fdopen(fd, "w") as f:
            json.dump(obj, f, ensure_ascii=False, indent=1)
            f.flush()
            os.fsync(f.fileno())
        os.replace(tmp, path)
    except BaseException:
        if os.path.exists(tmp):
            os.unlink(tmp)
        raise


class StageTimeout(Exception):
    pass


def _alarm_handler(signum, frame):
    raise StageTimeout()


def with_timeout(seconds, fn, *args, **kwargs):
    """Best-effort per-stage cap via SIGALRM (main thread only). A long
    C-level step returns to Python before the exception fires; the shell
    `timeout` in calc/aws_run_polyolefins.sh is the hard backstop."""
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


# ------------------------------------------------------------------ env knobs
ENV_DEFAULTS = {
    "POLY_STAGE_TIMEOUT": "4800", "POLY_TOTAL_BUDGET": "19800",
    "POLY_CC_POINTS": "2.8,2.0,1.6", "POLY_TS_CC": "2.0",
    "POLY_FACES": "si,re", "POLY_MF": "RKS", "POLY_XC": "pbe",
    "POLY_OPT_MAXSTEPS": "60", "POLY_CAS": "6,6",
    "POLY_AVAS_THRESHOLD": "0.2", "POLY_DEGEN_TOL": "0.05",
}


def parse_env(env):
    """Parse + validate all POLY_* knobs. Raises ValueError on nonsense
    (validated by --smoke before any AWS money is burned)."""
    g = lambda k: env.get(k, ENV_DEFAULTS[k])   # noqa: E731
    cfg = dict(
        stage_timeout=int(g("POLY_STAGE_TIMEOUT")),
        total_budget=int(g("POLY_TOTAL_BUDGET")),
        cc_points=[float(x) for x in g("POLY_CC_POINTS").split(",") if x.strip()],
        ts_cc=float(g("POLY_TS_CC")),
        faces=[s.strip().lower() for s in g("POLY_FACES").split(",") if s.strip()],
        mf=g("POLY_MF").upper(),
        xc=g("POLY_XC").lower(),
        opt_maxsteps=int(g("POLY_OPT_MAXSTEPS")),
        cas=tuple(int(x) for x in g("POLY_CAS").split(",")),
        avas_threshold=float(g("POLY_AVAS_THRESHOLD")),
        degen_tol_kcal=float(g("POLY_DEGEN_TOL")),
    )
    if cfg["stage_timeout"] <= 60 or cfg["total_budget"] <= 120:
        raise ValueError("POLY_STAGE_TIMEOUT/POLY_TOTAL_BUDGET too small")
    if not cfg["cc_points"] or not all(1.2 < c < 4.5 for c in cfg["cc_points"]):
        raise ValueError(f"POLY_CC_POINTS out of the sane 1.2-4.5 A window: "
                         f"{cfg['cc_points']}")
    if len(set(cfg["cc_points"])) != len(cfg["cc_points"]):
        raise ValueError("POLY_CC_POINTS has duplicates")
    if cfg["ts_cc"] not in cfg["cc_points"]:
        raise ValueError(f"POLY_TS_CC={cfg['ts_cc']} not in POLY_CC_POINTS")
    if not cfg["faces"] or not set(cfg["faces"]) <= {"si", "re"}:
        raise ValueError(f"POLY_FACES must be a subset of si,re: {cfg['faces']}")
    if cfg["mf"] not in ("RKS", "UKS"):
        raise ValueError(f"POLY_MF must be RKS or UKS, got {cfg['mf']}")
    if cfg["opt_maxsteps"] < 5:
        raise ValueError("POLY_OPT_MAXSTEPS < 5 makes no sense")
    ne, no = cfg["cas"]
    if not (ne >= 2 and ne % 2 == 0 and no >= ne // 2 and no <= 16):
        raise ValueError(f"POLY_CAS=({ne}e,{no}o) invalid (even ne>=2, "
                         f"ne/2<=no<=16)")
    if not (0.0 < cfg["avas_threshold"] < 1.0):
        raise ValueError("POLY_AVAS_THRESHOLD must be in (0,1)")
    return cfg


# ------------------------------------------------------- geometry: the ligand
def _unit(v):
    v = np.asarray(v, float)
    n = np.linalg.norm(v)
    if n < 1e-12:
        raise ValueError("zero-length vector in geometry construction")
    return v / n


ROT_C2Z = np.diag([-1.0, -1.0, 1.0])       # the C2 axis of the pocket = z


def make_methyl(c, axis, r_ch=CH_ME, phase=0.0):
    """3 H of a methyl at carbon `c`; `axis` = unit vector FROM the attached
    atom TOWARD c. X-C-H = 109.47 deg (same construction as calc/fe4s4...)."""
    a = _unit(axis)
    t = np.array([1.0, 0, 0]) if abs(a[0]) < 0.9 else np.array([0, 1.0, 0])
    u = _unit(np.cross(a, t))
    w = np.cross(a, u)
    th = math.radians(180.0 - 109.471)
    out = []
    for ph in (phase, phase + 2 * math.pi / 3, phase + 4 * math.pi / 3):
        hdir = math.cos(th) * a + math.sin(th) * (math.cos(ph) * u
                                                  + math.sin(ph) * w)
        out.append(("H", c + r_ch * hdir))
    return out


def _build_ring_a():
    """Ring A of the rac ligand: eta5-(3-methyl-C5H3) with the ipso carbon
    (k=0, bonded to Si) at the top-inner ring vertex and the methyl at k=2
    (ring position 3). ENANT flips the numbering sense = the other ligand
    enantiomer. Idealized: planar ring parallel to the plane normal to the
    Ti->centroid axis, ring H in-plane radial. Returns 12 atoms:
    C(k=0..4), H(k=1,3,4), C_methyl, 3 H_methyl."""
    alpha = math.radians(ALPHA_DEG)
    cen = TI_CEN * np.array([math.sin(alpha), 0.0, math.cos(alpha)])
    n = _unit(cen)
    e_y = np.array([0.0, 1.0, 0.0])                    # in-plane (n_y = 0)
    e_t = np.cross(n, e_y)                             # in-plane, +z-ish
    ring_c = []
    for k in range(5):
        phi = ENANT * 2 * math.pi / 5 * k
        ring_c.append(cen + R_RING * (math.cos(phi) * e_t
                                      + math.sin(phi) * e_y))
    atoms = [("C", p) for p in ring_c]
    for k in (1, 3, 4):                                # ring H (not ipso/Me)
        atoms.append(("H", ring_c[k] + CH_RING * _unit(ring_c[k] - cen)))
    rad = _unit(ring_c[2] - cen)                       # 3-position methyl
    cme = ring_c[2] + CC_ME * rad
    atoms.append(("C", cme))
    atoms.extend(make_methyl(cme, rad))
    return atoms, ring_c[0]                            # (atoms, ipso C pos)


def build_catalyst():
    """[H2Si(3-Me-C5H3)2 Ti-CH3]+ (32 atoms). Ring B is the EXACT C2(z) image
    of ring A -> the rac (C2-chiral) ligand by construction; a sigma(x) image
    would instead give the meso ligand (not built). Returns (atoms, idx)."""
    atoms = [("Ti", np.zeros(3))]
    ring_a, ipso_a = _build_ring_a()
    atoms += ring_a                                                # idx 1-12
    atoms += [(el, ROT_C2Z @ p) for el, p in ring_a]               # idx 13-24
    # ansa Si on the C2 axis: |Si - ipso| = SI_C, upper root (away from Ti)
    rho2 = float(ipso_a[0] ** 2 + ipso_a[1] ** 2)
    disc = SI_C ** 2 - rho2
    if disc <= 0:
        raise ValueError(f"ansa bridge unreachable: ipso rho^2={rho2:.3f} "
                         f"> Si-C^2={SI_C**2:.3f} — bad ALPHA/R_RING")
    si = np.array([0.0, 0.0, float(ipso_a[2]) + math.sqrt(disc)])
    atoms.append(("Si", si))                                       # idx 25
    # SiH2 in the yz plane (C2-symmetric pair), tilted away from the rings
    for sgn in (+1.0, -1.0):
        hdir = _unit(np.array([0.0, sgn * math.sin(math.radians(55.0)),
                               math.cos(math.radians(55.0))]))
        atoms.append(("H", si + SI_H * hdir))                      # idx 26,27
    # growing chain CH3 at the +y sigma site
    u_me = np.array([0.0, math.sin(math.radians(SIGMA_DEG)),
                     -math.cos(math.radians(SIGMA_DEG))])
    cm = TI_CM * u_me
    atoms.append(("C", cm))                                        # idx 28
    atoms.extend(make_methyl(cm, u_me))                            # idx 29-31
    idx = dict(ti=0, ligand=list(range(1, 28)), cm=28, si=25,
               ipso_a=1, ring_a_me=9)
    return atoms, idx


# --------------------------------------------- geometry: propene on a face
def _olefin_site():
    u_ol = np.array([0.0, -math.sin(math.radians(SIGMA_DEG)),
                     -math.cos(math.radians(SIGMA_DEG))])
    return TI_C1 * u_ol


def _solve_core(cc):
    """Positions of propene C1, C2 in the x=0 (equatorial) plane forming the
    Cossee-Arlman 4-center core Ti-Cm...C2=C1-Ti with the forming bond
    d(Cm,C2)=cc. Primary: C1 fixed at the olefin site, C2 = intersection of
    the circles r=CC_DB around C1 and r=cc around Cm (outer root, away from
    Ti). Fallback (small cc, circles disjoint): C2 by alternating projection
    (d to Cm = cc, Ti-C2 >= 2.05), C1 by a theta-scan (Ti-C1 ~ 2.25,
    biased toward the nominal site) — template only, relaxation refines."""
    cm = TI_CM * np.array([0.0, math.sin(math.radians(SIGMA_DEG)),
                           -math.cos(math.radians(SIGMA_DEG))])
    c1s = _olefin_site()
    p1, p2 = c1s[1:], cm[1:]                       # 2D work in (y,z), x=0
    d = float(np.linalg.norm(p2 - p1))
    if abs(cc - CC_DB) + 0.05 < d < cc + CC_DB - 0.05:
        a = (CC_DB ** 2 - cc ** 2 + d ** 2) / (2 * d)
        h = math.sqrt(max(CC_DB ** 2 - a ** 2, 0.0))
        u = (p2 - p1) / d
        base = p1 + a * u
        perp = np.array([-u[1], u[0]])
        cands = [base + h * perp, base - h * perp]
        c2_2d = max(cands, key=lambda p: np.linalg.norm(p))   # away from Ti
        c1 = c1s
        c2 = np.array([0.0, c2_2d[0], c2_2d[1]])
        mode = "circle-intersection"
    else:
        c2_2d = p2 + cc * (p1 - p2) / d
        for _ in range(60):                        # alternating projections
            r = float(np.linalg.norm(c2_2d))
            if r < 2.05:
                c2_2d = c2_2d / r * 2.05
            v = c2_2d - p2
            c2_2d = p2 + cc * v / np.linalg.norm(v)
        c2 = np.array([0.0, c2_2d[0], c2_2d[1]])
        best = None                                # theta-scan for C1
        for th in np.linspace(0.0, 2 * math.pi, 720, endpoint=False):
            c1t = c2 + CC_DB * np.array([0.0, math.cos(th), math.sin(th)])
            cost = ((np.linalg.norm(c1t) - 2.25) ** 2
                    + 0.2 * np.linalg.norm(c1t - c1s) ** 2)
            if best is None or cost < best[0]:
                best = (cost, c1t)
        c1 = best[1]
        mode = "fallback-projection"
    return c1, c2, cm, mode


# planar propene in its local frame: C1 at origin, C2 on +x, skeleton in z=0
_PROP_LOCAL = {
    "C1": np.array([0.0, 0.0, 0.0]),
    "C2": np.array([CC_DB, 0.0, 0.0]),
    "H11": 1.085 * np.array([math.cos(math.radians(120)),
                             math.sin(math.radians(120)), 0.0]),
    "H12": 1.085 * np.array([math.cos(math.radians(-120)),
                             math.sin(math.radians(-120)), 0.0]),
    "C3": np.array([CC_DB, 0.0, 0.0]) + 1.50 * np.array(
        [math.cos(math.radians(60)), math.sin(math.radians(60)), 0.0]),
    "H2": np.array([CC_DB, 0.0, 0.0]) + 1.085 * np.array(
        [math.cos(math.radians(-60)), math.sin(math.radians(-60)), 0.0]),
}


def place_propene(g_c1, g_c2, p_ti, face_sign):
    """Rigid propene with C1 -> g_c1, C2 -> g_c2 and the molecular plane
    oriented so that Ti sits on the face_sign side of the local +z normal.
    face_sign = +1/-1 selects which enantioface coordinates; the CIP label is
    computed afterwards by face_label() — never assumed. Returns 9 atoms:
    C1, C2, C3(methyl C), H11, H12, H2, 3 x H(methyl)."""
    e1 = _unit(g_c2 - g_c1)
    w = p_ti - 0.5 * (g_c1 + g_c2)
    n = w - np.dot(w, e1) * e1
    if np.linalg.norm(n) < 0.1:
        raise ValueError("Ti nearly collinear with C1=C2 — bad core template")
    zg = float(face_sign) * _unit(n)               # local +z -> Ti side (+1)
    yg = np.cross(zg, e1)                          # right-handed: y = z cross x
    R = np.column_stack([e1, yg, zg])
    g = {k: g_c1 + R @ v for k, v in _PROP_LOCAL.items()}
    atoms = [("C", g["C1"]), ("C", g["C2"]), ("C", g["C3"]),
             ("H", g["H11"]), ("H", g["H12"]), ("H", g["H2"])]
    atoms += make_methyl(g["C3"], _unit(g["C3"] - g["C2"]))
    return atoms


def face_label(p_ti, p_c1, p_c2, p_c3):
    """Geometric CIP face of the prochiral C2 as seen from Ti. Priorities at
    sp2 C2: C1 (=CH2, duplicated C) > C3 (CH3) > H. With s = ((C1-C2) x
    (C3-C2)) . (Ti-C2):  s>0 means the cross product points to Ti, i.e. an
    observer on the Ti side sees C1->C3 COUNTERclockwise = descending
    priorities counterclockwise = si; s<0 = re. Cross-checked in --smoke
    (the two face_signs must give the two different labels)."""
    s = float(np.dot(np.cross(p_c1 - p_c2, p_c3 - p_c2), p_ti - p_c2))
    return "si" if s > 0 else "re"


def build_system(cc, face_sign):
    """Full 41-atom [catalyst + propene(face)] template at scan point cc.
    Returns (atoms, idx) with idx['c1'], idx['c2'], idx['c3'] = propene.
    The chain CH3 hydrogens are re-oriented away from BOTH Ti and the
    incoming C2 (blended axis) — with the plain Ti->Cm axis they collide
    with propene at product-like cc (caught by the smoke min-distance gate);
    the relaxation refines the exact umbrella angle anyway."""
    atoms, idx = build_catalyst()
    g_c1, g_c2, cm, mode = _solve_core(cc)
    p_ti = atoms[idx["ti"]][1]
    axis = _unit(cm - 0.5 * (p_ti + g_c2))
    for k, (el, p) in enumerate(make_methyl(cm, axis)):
        atoms[idx["cm"] + 1 + k] = (el, p)
    prop = place_propene(g_c1, g_c2, p_ti, face_sign)
    base = len(atoms)
    atoms += prop
    idx.update(c1=base, c2=base + 1, c3=base + 2, core_mode=mode)
    return atoms, idx


# ------------------------------------------------------------ geometry checks
def _set_match_max(atoms_a, atoms_b):
    """max over atoms of A of the distance to the nearest same-element atom
    of B (0 for identical sets up to permutation)."""
    worst = 0.0
    for el, p in atoms_a:
        d = min(np.linalg.norm(p - q) for e2, q in atoms_b if e2 == el)
        worst = max(worst, d)
    return worst


def check_c2(atoms, idx, tol=1e-8):
    """Ti + ligand must be EXACTLY C2(z)-symmetric (by construction)."""
    sub = [atoms[0]] + [atoms[i] for i in idx["ligand"]]
    rot = [(el, ROT_C2Z @ p) for el, p in sub]
    return _set_match_max(rot, sub)


def check_chiral(atoms, idx, tol=0.3):
    """The ligand must NOT be superimposable on its sigma(yz)/sigma(xz)
    mirror images (heuristic anti-meso / anti-achiral gate: a match would
    mean the pocket cannot distinguish faces even in principle)."""
    sub = [atoms[i] for i in idx["ligand"]]
    worst = {}
    for name, M in (("sigma_yz", np.diag([-1.0, 1, 1])),
                    ("sigma_xz", np.diag([1.0, -1, 1]))):
        worst[name] = round(_set_match_max([(el, M @ p) for el, p in sub],
                                           sub), 3)
    return worst, all(v > tol for v in worst.values())


def min_pair_dist(atoms):
    P = np.array([p for _, p in atoms])
    D = np.linalg.norm(P[:, None, :] - P[None, :, :], axis=-1)
    np.fill_diagonal(D, np.inf)
    return float(D.min())


def forming_bond(atoms, idx):
    return float(np.linalg.norm(atoms[idx["cm"]][1] - atoms[idx["c2"]][1]))


def write_xyz(atoms, path, comment):
    with open(path, "w") as fh:
        fh.write(f"{len(atoms)}\n{comment}\n")
        for e, p in atoms:
            fh.write(f"{e} {p[0]:.6f} {p[1]:.6f} {p[2]:.6f}\n")


def atoms_to_lists(atoms):
    return [[el, round(float(p[0]), 6), round(float(p[1]), 6),
             round(float(p[2]), 6)] for el, p in atoms]


def lists_to_atoms(lists):
    return [(el, np.array([x, y, z])) for el, x, y, z in lists]


# ----------------------------------------------------------------- constraints
def write_constraints(path, idx, cc):
    """geomeTRIC constraints (1-based indices): freeze the whole ansa ligand,
    fix the forming C-C bond. Ti, chain CH3 and propene relax."""
    lig = idx["ligand"]
    assert lig == list(range(lig[0], lig[-1] + 1)), "ligand block not contiguous"
    text = ("$freeze\n"
            f"xyz {lig[0] + 1}-{lig[-1] + 1}\n"
            "$set\n"
            f"distance {idx['cm'] + 1} {idx['c2'] + 1} {cc:.4f}\n")
    with open(path, "w") as fh:
        fh.write(text)
    return text


# ------------------------------------------------------------------- QC layer
def make_mol(atoms, charge=CHARGE, spin=SPIN):
    maxmem = int(os.environ.get("PYSCF_MAX_MEMORY", "8000"))
    return gto.M(atom=[(el, tuple(p)) for el, p in atoms], basis=BASIS,
                 charge=charge, spin=spin, verbose=0, max_memory=maxmem)


def new_mf(mol, cfg):
    if cfg["mf"] == "UKS":
        mf = dft.UKS(mol).density_fit()
    else:
        mf = dft.RKS(mol).density_fit()
    mf.xc = cfg["xc"]
    mf.conv_tol = 1e-8
    return mf


def converge_scf(mf, tag=""):
    """DIIS; on stall a Newton (SOSCF) restart from the current MOs."""
    mf.max_cycle = 80
    mf.kernel()
    if not mf.converged:
        say(f"    [{tag}] DIIS incomplete -> Newton (SOSCF) restart")
        mfn = mf.newton()
        mfn.max_cycle = 60
        mfn.kernel(mf.mo_coeff, mf.mo_occ)
        mf.mo_coeff, mf.mo_occ = mfn.mo_coeff, mfn.mo_occ
        mf.mo_energy = mfn.mo_energy
        mf.e_tot, mf.converged = float(mfn.e_tot), bool(mfn.converged)
    return mf


def optimize_point(cfg, atoms0, idx, cc, tag, alarm_s, workdir):
    """Constrained relaxation (frozen ligand + fixed forming C-C) at
    (R/U)KS-PBE/def2-SVP/DF via geomeTRIC, then a final SCF on the relaxed
    geometry. Returns (info dict, relaxed atoms or None)."""
    from pyscf.geomopt import geometric_solver
    cons = os.path.join(workdir, f"constraints_{tag}.txt")
    write_constraints(cons, idx, cc)
    info = {"cc_target": cc, "constraints_file": os.path.basename(cons),
            "maxsteps": cfg["opt_maxsteps"], "mean_field": f"{cfg['mf']}-"
            f"{cfg['xc'].upper()}/{BASIS}+DF", "frozen": "ansa ligand "
            f"(atoms {idx['ligand'][0] + 1}-{idx['ligand'][-1] + 1}, 1-based)"}
    t0 = time.time()
    try:
        mol = make_mol(atoms0)
        mf = new_mf(mol, cfg)
        conv, mol_opt = with_timeout(
            alarm_s, geometric_solver.kernel, mf,
            assert_convergence=False, constraints=cons,
            maxsteps=cfg["opt_maxsteps"])
        info["opt_converged"] = bool(conv)
        els = [mol_opt.atom_symbol(i) for i in range(mol_opt.natm)]
        coords = mol_opt.atom_coords(unit="ANG")
        atoms_opt = [(el, np.array(c)) for el, c in zip(els, coords)]
        info["cc_actual"] = round(forming_bond(atoms_opt, idx), 4)
        info["constraint_honored"] = bool(abs(info["cc_actual"] - cc) < 0.01)
        # final SCF at the relaxed geometry (fresh, recorded)
        mf2 = converge_scf(new_mf(make_mol(atoms_opt), cfg), tag=tag)
        info["e_dft"] = float(mf2.e_tot)
        info["scf_converged"] = bool(mf2.converged)
        info["face_after_opt"] = face_label(
            atoms_opt[idx["ti"]][1], atoms_opt[idx["c1"]][1],
            atoms_opt[idx["c2"]][1], atoms_opt[idx["c3"]][1])
        info["min_pair_dist"] = round(min_pair_dist(atoms_opt), 3)
    except StageTimeout:
        info["opt_converged"] = False
        info["timed_out_after_s"] = int(alarm_s)
        atoms_opt = None
        say(f"    [{tag}] TIMED OUT after {alarm_s}s")
    info["seconds"] = round(time.time() - t0, 1)
    return info, atoms_opt


# ------------------------------------------------------- CASCI on the TS pair
def _sqrtm_sym(s):
    w, v = np.linalg.eigh(s)
    w = np.clip(w, 1e-12, None)
    return (v * np.sqrt(w)) @ v.T


def target_ao_idx(mol, idx):
    """AO indices of Ti 3d + 2p of the three reacting carbons (chain Cm,
    propene C1, C2) — the 4-center manifold, as in the achiral diary stage."""
    react_c = {idx["cm"], idx["c1"], idx["c2"]}
    out = []
    for i, lab in enumerate(mol.ao_labels()):
        parts = lab.split()
        a, el, shell = int(parts[0]), parts[1], parts[2]
        if el == "Ti" and shell.startswith("3d"):
            out.append(i)
        elif el == "C" and a in react_c and shell.startswith("2p"):
            out.append(i)
    return np.array(out, dtype=int)


def truncate_active(mol, mf, mo, ncas_avas, nelec_avas, ne_t, no_t, tao):
    """Truncate the AVAS window (core|active|virt in `mo`, occupied-active
    first) to a tractable CAS(ne_t, no_t): keep the ne_t/2 occupied-active
    and (no_t - ne_t/2) virtual-active orbitals with the HIGHEST Loewdin
    character on the target AOs; the dropped occupied-active go to the core
    (electron count stays consistent). Everything is recorded."""
    ncore = (mol.nelectron - nelec_avas) // 2
    nocc_act = nelec_avas // 2
    S = mol.intor("int1e_ovlp")
    shalf = _sqrtm_sym(S)
    act = mo[:, ncore:ncore + ncas_avas]
    char = ((shalf @ act) ** 2)[tao, :].sum(axis=0)
    dm = mf.make_rdm1()
    if dm.ndim == 3:
        dm = dm[0] + dm[1]
    occ_diag = np.einsum("pi,pq,qi->i", act, S @ dm @ S, act)
    need_o = ne_t // 2
    need_v = no_t - need_o
    if nocc_act < need_o or (ncas_avas - nocc_act) < need_v:
        need_o = min(need_o, nocc_act)
        need_v = min(need_v, ncas_avas - nocc_act)
    sel_o = sorted(np.argsort(char[:nocc_act])[::-1][:need_o].tolist())
    sel_v = sorted((nocc_act
                    + np.argsort(char[nocc_act:])[::-1][:need_v]).tolist())
    drop_o = [i for i in range(nocc_act) if i not in set(sel_o)]
    drop_v = [i for i in range(nocc_act, ncas_avas) if i not in set(sel_v)]
    order = (list(range(ncore)) + [ncore + i for i in drop_o]
             + [ncore + i for i in sel_o] + [ncore + i for i in sel_v]
             + [ncore + i for i in drop_v]
             + list(range(ncore + ncas_avas, mo.shape[1])))
    diag = {"avas_space": [int(nelec_avas), int(ncas_avas)],
            "kept_char": [round(float(char[i]), 3) for i in sel_o + sel_v],
            "kept_occ": [round(float(occ_diag[i]), 3) for i in sel_o + sel_v],
            "cas_used": [2 * need_o, need_o + need_v],
            "note": "AVAS window truncated by Loewdin character on Ti3d + "
                    "reacting-C 2p; dropped occupied-active moved to core"}
    return mo[:, order], 2 * need_o, need_o + need_v, diag


def casci_at_geometry(cfg, atoms, idx, tag, alarm_s,
                      charge=CHARGE, spin=SPIN, labels=None, cas=None):
    """DF-RHF -> AVAS(Ti 3d + reacting C 2p) -> character-truncated
    CASCI(ne,no) -> E, NOON, N_u = sum n(2-n). Used on both TS geometries
    (and, in --smoke, on free propene with labels=['C 2p'])."""
    info = {}
    t0 = time.time()
    try:
        mol = make_mol(atoms, charge=charge, spin=spin)
        mf = scf.RHF(mol).density_fit()
        mf.conv_tol = 1e-8
        with_timeout(alarm_s, converge_scf, mf, tag=f"{tag}/RHF")
        info["e_hf"] = float(mf.e_tot)
        info["hf_converged"] = bool(mf.converged)
        if labels is None:
            labels = ["Ti 3d", f"{idx['cm']} C 2p",
                      f"{idx['c1']} C 2p", f"{idx['c2']} C 2p"]
        info["avas_labels"] = labels
        ncas_a, nelec_a, mo = avas.avas(mf, labels,
                                        threshold=cfg["avas_threshold"])
        ncas_a, nelec_a = int(ncas_a), int(nelec_a)
        ne_t, no_t = cas if cas is not None else cfg["cas"]
        tao = (target_ao_idx(mol, idx) if labels[0].startswith("Ti")
               else np.array(mol.search_ao_label(labels), dtype=int))
        mo2, ne_f, no_f, diag = truncate_active(mol, mf, mo, ncas_a, nelec_a,
                                                ne_t, no_t, tao)
        info["active_space"] = diag
        mc = mcscf.CASCI(mf, no_f, ne_f).density_fit()
        mc.verbose = 0
        try:
            fci.addons.fix_spin_(mc.fcisolver, ss=0)
        except Exception as exc:
            info["fix_spin_note"] = repr(exc)
        with_timeout(max(int(alarm_s - (time.time() - t0)), 60),
                     mc.kernel, mo2)
        info["e_casci"] = float(mc.e_tot)
        info["casci_converged"] = bool(mc.converged)
        dm1 = mc.fcisolver.make_rdm1(mc.ci, mc.ncas, mc.nelecas)
        noon = np.clip(np.linalg.eigvalsh(dm1)[::-1], 0.0, 2.0)
        info["noon"] = [round(float(x), 3) for x in noon]
        info["n_u"] = round(float(np.sum(noon * (2.0 - noon))), 3)
        info["df_note"] = "DF used throughout, incl. CASCI active integrals"
    except StageTimeout:
        info["timed_out_after_s"] = int(alarm_s)
        say(f"    [{tag}] CASCI stage TIMED OUT after {alarm_s}s")
    info["seconds"] = round(time.time() - t0, 1)
    return info


# ----------------------------------------------------------------- summarizer
def summarize(res, cfg):
    """(Re)derive the stereo block from whatever points exist. Every gate is
    explicit; numbers appear ONLY if both faces delivered the needed point."""
    pts = res.get("points", {})
    ts_key = f"cc_{cfg['ts_cc']:g}"
    pi_key = f"cc_{max(cfg['cc_points']):g}"
    st = {"convention": "DDE = E(si) - E(re); positive = re-face path lower; "
                        "the mirror ligand enantiomer gives exactly -DDE",
          "ts_point": ts_key, "pi_point": pi_key}
    have = lambda f, k: pts.get(f, {}).get(k, {}).get("e_dft") is not None  # noqa: E731
    gates = {}
    if have("si", ts_key) and have("re", ts_key):
        e_si, e_re = pts["si"][ts_key]["e_dft"], pts["re"][ts_key]["e_dft"]
        st["ddE_ts_dft_kcal"] = round((e_si - e_re) * KCAL, 3)
        for f in ("si", "re"):
            p = pts[f][ts_key]
            gates[f"{f}_ts_ok"] = bool(p.get("opt_converged")
                                       and p.get("scf_converged")
                                       and p.get("constraint_honored"))
            gates[f"{f}_face_preserved"] = bool(p.get("face_after_opt") == f)
        dd = abs(st["ddE_ts_dft_kcal"])
        gates["degenerate_within_noise"] = bool(dd < cfg["degen_tol_kcal"])
        gates["suspiciously_exact_zero"] = bool(dd < 1e-4)
        if gates["degenerate_within_noise"]:
            st["finding"] = (f"|DDE_TS| = {dd} kcal/mol < "
                             f"{cfg['degen_tol_kcal']} — the 3-methyl pocket "
                             "did NOT resolve the faces at this level; honest "
                             "outcome = this minimal chiral model is too weak "
                             "(next step: bulkier 3-substituent, e.g. tBu, or "
                             "true indenyl)" + (" — AND it is exact to <1e-4, "
                             "which smells like a construction bug (check "
                             "suspiciously_exact_zero)" if
                             gates["suspiciously_exact_zero"] else ""))
        if have("si", pi_key) and have("re", pi_key):
            b_si = (pts["si"][ts_key]["e_dft"] - pts["si"][pi_key]["e_dft"]) * KCAL
            b_re = (pts["re"][ts_key]["e_dft"] - pts["re"][pi_key]["e_dft"]) * KCAL
            st["barrier_dft_kcal"] = {"si": round(b_si, 3), "re": round(b_re, 3)}
            st["ddE_barrier_dft_kcal"] = round(b_si - b_re, 3)
            st["ddE_pi_dft_kcal"] = round(
                (pts["si"][pi_key]["e_dft"] - pts["re"][pi_key]["e_dft"]) * KCAL, 3)
    cas = res.get("casci", {})
    if all(cas.get(f, {}).get("e_casci") is not None for f in ("si", "re")):
        st["ddE_ts_casci_kcal"] = round(
            (cas["si"]["e_casci"] - cas["re"]["e_casci"]) * KCAL, 3)
        st["ddE_ts_hf_kcal"] = round(
            (cas["si"]["e_hf"] - cas["re"]["e_hf"]) * KCAL, 3)
        gates["casci_same_space"] = bool(
            cas["si"]["active_space"]["cas_used"]
            == cas["re"]["active_space"]["cas_used"])
        gates["casci_converged_both"] = bool(
            cas["si"].get("casci_converged") and cas["re"].get("casci_converged"))
        if "ddE_ts_dft_kcal" in st and gates["casci_same_space"]:
            st["mr_correction_kcal"] = round(
                st["ddE_ts_casci_kcal"] - st["ddE_ts_dft_kcal"], 3)
            st["mr_note"] = ("CASCI(DF)/RHF vs DFT difference of DDE_TS; "
                             "CASCI lacks dynamic correlation and sits on an "
                             "HF (not KS) reference — a diagnostic of "
                             "multireference sensitivity, not a final number")
        st["n_u_ts"] = {f: cas[f].get("n_u") for f in ("si", "re")}
    st["gates"] = gates
    st["reliable"] = bool(gates) and all(
        v for k, v in gates.items()
        if k not in ("degenerate_within_noise", "suspiciously_exact_zero"))
    res["stereo"] = st
    return st


# ------------------------------------------------------------------- real run
def face_sign_map(cfg):
    """Determine which construction sign gives which CIP face label at the
    TS scan point — computed, never assumed."""
    m = {}
    for f in (+1, -1):
        atoms, idx = build_system(cfg["ts_cc"], f)
        lab = face_label(atoms[idx["ti"]][1], atoms[idx["c1"]][1],
                         atoms[idx["c2"]][1], atoms[idx["c3"]][1])
        m[lab] = f
    if set(m) != {"si", "re"}:
        raise RuntimeError(f"face construction degenerate: {m}")
    return m


def main_real(args):
    cfg = parse_env(os.environ)
    res = json.load(open(RESULTS)) if os.path.exists(RESULTS) else {}
    res.setdefault("meta", {})
    res["meta"].update({
        "issue": "polyolefins track, GitHub issue #11 (chiral model for a "
                 "nonzero stereo-DDE++)",
        "system": "[H2Si(eta5-3-Me-C5H3)2Ti-CH3]+ + C3H6, 41 atoms, q=+1, "
                  "closed-shell d0",
        "model_choice": "minimal rac-C2 ansa-titanocene mimic; indenyl "
                        "truncated to 3-methyl, SiMe2 bridge to SiH2; one "
                        "enantiomer fixed (ENANT=+1)",
        "model_limits": ["3-Me is a weak stereo-differentiator: |DDE| may be "
                         "below noise (flagged, not hidden)",
                         "idealized FROZEN ligand (parallel planar rings, "
                         "in-plane ring H): absolute barriers biased, only "
                         "the si-re DIFFERENCE is the deliverable",
                         "cc=ts point is a constrained-scan TS approximation, "
                         "no saddle optimization / frequencies",
                         "def2-SVP, DF everywhere, no solvent/counterion"],
        "method": f"{cfg['mf']}-{cfg['xc'].upper()}/{BASIS}+DF, geomeTRIC "
                  f"constrained relaxation (ligand frozen, d(Cm-C2) fixed); "
                  f"CASCI(AVAS Ti3d + reacting-C 2p -> "
                  f"CAS{tuple(cfg['cas'])})//DF-RHF at both TS points",
        "env": {k: os.environ.get(k, ENV_DEFAULTS[k]) for k in ENV_DEFAULTS},
        "started": time.strftime("%Y-%m-%d %H:%M:%S"),
    })
    res["status"] = "running"
    atomic_save(res, RESULTS)
    say(f"chiral DDE job: faces={cfg['faces']} points={cfg['cc_points']} "
        f"threads={lib.num_threads()} budget={cfg['total_budget']}s")

    # model sanity (cheap, every run): C2 pocket + chirality + template dump
    atoms0, idx0 = build_system(cfg["ts_cc"], +1)
    dev = check_c2(atoms0, idx0)
    mirror, chiral = check_chiral(atoms0, idx0)
    res["model"] = {"n_atoms": len(atoms0), "nelec": 154,
                    "c2_deviation": float(f"{dev:.2e}"),
                    "mirror_mismatch": mirror, "pocket_chiral": bool(chiral),
                    "min_pair_dist_template": round(min_pair_dist(atoms0), 3)}
    if dev > 1e-6 or not chiral:
        res["status"] = "aborted: pocket not C2-chiral — construction bug"
        atomic_save(res, RESULTS)
        sys.exit(res["status"])
    write_xyz(atoms0, XYZ_MODEL, "template: [H2Si(3-MeC5H3)2Ti-CH3]+ + "
              f"propene, cc={cfg['ts_cc']} (polyolefins issue #11)")
    fmap = face_sign_map(cfg)
    res["model"]["face_sign_map"] = fmap
    atomic_save(res, RESULTS)

    workdir = tempfile.mkdtemp(prefix="poly_chiral_")
    pts = res.setdefault("points", {})
    last_secs = None
    try:
        # --- constrained scan: points in given order, faces interleaved ---
        for cc in cfg["cc_points"]:
            key = f"cc_{cc:g}"
            for face in cfg["faces"]:
                fp = pts.setdefault(face, {})
                if fp.get(key, {}).get("e_dft") is not None:
                    say(f"[scan] {face}/{key} already done (restart) — skip")
                    continue
                need = max(1200, 1.5 * last_secs) if last_secs else 1200
                if cfg["total_budget"] - elapsed() < need:
                    fp[key] = {"skipped": f"budget: "
                               f"{int(cfg['total_budget'] - elapsed())}s left "
                               f"< {int(need)}s needed"}
                    atomic_save(res, RESULTS)
                    say(f"[scan] skip {face}/{key} — budget")
                    continue
                # warm start from the nearest already-relaxed point of the
                # same face; otherwise the analytic template
                prev = [v for v in fp.values() if isinstance(v, dict)
                        and v.get("coords")]
                if prev:
                    src = min(prev, key=lambda v: abs(v["cc_target"] - cc))
                    atoms = lists_to_atoms(src["coords"])
                    start = f"warm from cc_{src['cc_target']:g}"
                else:
                    atoms, _ = build_system(cc, fmap[face])
                    start = "analytic template"
                _, idx = build_system(cc, fmap[face])   # index map only
                alarm = int(min(cfg["stage_timeout"],
                                cfg["total_budget"] - elapsed() - 120))
                say(f"[scan] {face}/{key}: constrained opt ({start}), "
                    f"alarm {alarm}s")
                info, atoms_opt = optimize_point(cfg, atoms, idx, cc,
                                                 f"{face}_{key}", alarm, workdir)
                info["start"] = start
                if atoms_opt is not None:
                    info["coords"] = atoms_to_lists(atoms_opt)
                    if cc == cfg["ts_cc"]:
                        write_xyz(atoms_opt, XYZ_TS_TPL.format(face=face),
                                  f"relaxed TS-point (cc={cc}), face={face}, "
                                  f"E_dft={info.get('e_dft')}")
                fp[key] = info
                atomic_save(res, RESULTS)
                say(f"[scan] {face}/{key}: conv={info.get('opt_converged')} "
                    f"E={info.get('e_dft')} face_after="
                    f"{info.get('face_after_opt')} ({info['seconds']}s)")
                last_secs = info["seconds"]
            summarize(res, cfg)
            atomic_save(res, RESULTS)

        # --- CASCI on both TS geometries ---
        cas = res.setdefault("casci", {})
        ts_key = f"cc_{cfg['ts_cc']:g}"
        for face in cfg["faces"]:
            if cas.get(face, {}).get("e_casci") is not None:
                say(f"[casci] {face} already done (restart) — skip")
                continue
            src = pts.get(face, {}).get(ts_key, {})
            if not src.get("coords"):
                cas[face] = {"skipped": f"no relaxed {ts_key} geometry"}
                atomic_save(res, RESULTS)
                continue
            if cfg["total_budget"] - elapsed() < 600:
                cas[face] = {"skipped": "budget"}
                atomic_save(res, RESULTS)
                continue
            alarm = int(min(cfg["stage_timeout"],
                            cfg["total_budget"] - elapsed() - 60))
            say(f"[casci] {face} @ {ts_key}: RHF -> AVAS -> "
                f"CASCI{tuple(cfg['cas'])}, alarm {alarm}s")
            _, idx = build_system(cfg["ts_cc"], fmap[face])
            cas[face] = casci_at_geometry(cfg, lists_to_atoms(src["coords"]),
                                          idx, f"casci_{face}", alarm)
            atomic_save(res, RESULTS)
            say(f"[casci] {face}: E_casci={cas[face].get('e_casci')} "
                f"n_u={cas[face].get('n_u')} ({cas[face]['seconds']}s)")

        st = summarize(res, cfg)
        res["status"] = "ok" if "ddE_ts_dft_kcal" in st else "partial"
    except Exception as e:
        res["status"] = f"aborted: {type(e).__name__}: {e}"
        say("ABORT: " + res["status"])
    finally:
        res["meta"]["total_seconds"] = round(elapsed(), 1)
        summarize(res, cfg)
        atomic_save(res, RESULTS)
        say(f"wrote {RESULTS}  status={res['status']}  "
            f"(total {res['meta']['total_seconds']}s)")
    return 0 if res["status"] in ("ok", "partial") else 1


# ---------------------------------------------------------------------- smoke
def main_smoke(args):
    """Preflight, <=3 min on 1 core: geometry assembly for every face x point,
    xyz dumps, C2/chirality/face/constraint/env/ast validation, geomeTRIC
    import gate, save/restart round-trip, and ONE tiny SCF+AVAS+CASCI on free
    propene (skippable with --no-scf). NO SCF on the 41-atom system."""
    t0 = time.time()
    fails = []
    ok = lambda c, m: (say(f"  [{'ok' if c else 'FAIL'}] {m}"),          # noqa: E731
                       None if c else fails.append(m))

    say(f"SMOKE start (threads={lib.num_threads()})")
    smoked = tempfile.mkdtemp(prefix="poly_chiral_smoke_")
    say(f"  smoke dir: {smoked}")

    # 1. ast self-parse (config/code syntax gate)
    src = open(os.path.abspath(__file__)).read()
    try:
        ast.parse(src)
        ok(True, "ast.parse(self)")
    except SyntaxError as e:
        ok(False, f"ast.parse(self): {e}")

    # 2. env knobs: defaults + a representative override set must validate;
    #    a broken set must be rejected
    try:
        cfg = parse_env(os.environ)
        parse_env({"POLY_CC_POINTS": "3.0,2.1", "POLY_TS_CC": "2.1",
                   "POLY_FACES": "si", "POLY_MF": "uks", "POLY_CAS": "4,5"})
        ok(True, f"env knobs parse (defaults: points={cfg['cc_points']}, "
                 f"cas={cfg['cas']}, faces={cfg['faces']})")
    except ValueError as e:
        ok(False, f"env parse: {e}")
        return 1
    try:
        parse_env({"POLY_TS_CC": "9.9"})
        ok(False, "bad env accepted (POLY_TS_CC=9.9 not rejected)")
    except ValueError:
        ok(True, "bad env rejected (POLY_TS_CC outside POLY_CC_POINTS)")

    # 3. catalyst: exact C2 pocket, chiral (not superimposable on mirrors)
    cat, cidx = build_catalyst()
    dev = check_c2(cat + [("C", np.array([9, 9, 9.0]))] * 0, cidx)
    ok(dev < 1e-8, f"ligand+Ti C2(z) symmetry exact (max dev {dev:.1e} A)")
    atoms_t, idx_t = build_system(cfg["ts_cc"], +1)
    mirror, chiral = check_chiral(atoms_t, idx_t)
    ok(chiral, f"pocket chiral: mirror mismatch {mirror} (> 0.3 A)")
    nel = sum(gto.charge(el) for el, _ in atoms_t) - CHARGE
    ok(nel == 154 and nel % 2 == 0,
       f"electron count {nel} (even, closed-shell cation)")

    # 4. all faces x points: build, constraint, face labels, clashes, xyz
    labels = {}
    for f in (+1, -1):
        labels[f] = set()
        for cc in cfg["cc_points"]:
            atoms, idx = build_system(cc, f)
            d = forming_bond(atoms, idx)
            ok(abs(d - cc) < 2e-3,
               f"f={f:+d} cc={cc}: forming bond {d:.4f} A "
               f"({idx['core_mode']})")
            md = min_pair_dist(atoms)
            ok(md > 0.85, f"f={f:+d} cc={cc}: min pair distance {md:.3f} A")
            if md < 1.0:      # below any bond length = a real nonbonded squeeze
                say(f"    [note] f={f:+d} cc={cc}: tight contact {md:.3f} A "
                    "(steric differentiation region — relaxation will decide)")
            lab = face_label(atoms[idx["ti"]][1], atoms[idx["c1"]][1],
                             atoms[idx["c2"]][1], atoms[idx["c3"]][1])
            labels[f].add(lab)
            path = os.path.join(smoked, f"poly_chiral_f{f:+d}_cc{cc:g}.xyz")
            write_xyz(atoms, path, f"smoke template face_sign={f:+d} "
                                   f"label={lab} cc={cc}")
        ok(len(labels[f]) == 1,
           f"f={f:+d}: face label consistent across points ({labels[f]})")
    ok(labels[+1] != labels[-1] and labels[+1] | labels[-1] == {"si", "re"},
       f"face signs map to the two distinct labels "
       f"(+1->{labels[+1]}, -1->{labels[-1]})")
    say(f"  xyz templates written to {smoked}")

    # 5. constraints file: generate + reparse
    cons = os.path.join(smoked, "constraints.txt")
    text = write_constraints(cons, idx_t, cfg["ts_cc"])
    back = open(cons).read()
    ok(back == text and f"distance {idx_t['cm'] + 1} {idx_t['c2'] + 1}" in back
       and "xyz 2-28" in back, "constraints file round-trip "
       f"(freeze 2-28, distance {idx_t['cm'] + 1}-{idx_t['c2'] + 1})")

    # 6. geomeTRIC availability (needed for the real run)
    try:
        import geometric                                   # noqa: F401
        from pyscf.geomopt import geometric_solver         # noqa: F401
        ok(True, f"geomeTRIC importable (v{geometric.__version__})")
    except ImportError as e:
        ok(False, f"geomeTRIC missing: {e}")

    # 7. save/restart round-trip of the results machinery
    rp = os.path.join(smoked, "results_roundtrip.json")
    probe = {"points": {"si": {"cc_2": {"e_dft": -1.0}}}}
    atomic_save(probe, rp)
    ok(json.load(open(rp)) == probe, "atomic save/reload round-trip")

    # 8. ONE tiny SCF gate (free propene, ~seconds): exercises RHF + AVAS +
    #    character truncation + CASCI + NOON/N_u end to end
    if not args.no_scf:
        g_c1 = np.zeros(3)
        g_c2 = np.array([CC_DB, 0.0, 0.0])
        fake_ti = np.array([0.7, 0.0, 2.3])
        prop = place_propene(g_c1, g_c2, fake_ti, +1)
        pidx = dict(cm=0, c1=0, c2=1, c3=2)   # cm unused for labels here
        t1 = time.time()
        info = casci_at_geometry(cfg, prop, pidx, "smoke_propene",
                                 alarm_s=150, charge=0, spin=0,
                                 labels=["C 2p"], cas=(2, 2))
        dt = time.time() - t1
        ok(info.get("hf_converged") is True,
           f"propene RHF/def2-SVP converged (E={info.get('e_hf'):.6f})")
        ok(info.get("casci_converged") is True and
           info.get("e_casci", 1e9) <= info.get("e_hf", -1e9) + 1e-7,
           f"CASCI(2,2) variational vs RHF "
           f"(dE={(info.get('e_casci', 0) - info.get('e_hf', 0)):.6f} Ha)")
        ok(0.0 <= info.get("n_u", 9) <= 1.0,
           f"N_u sane for a pi(2,2) space: {info.get('n_u')} "
           f"(NOON {info.get('noon')})")
        say(f"  tiny-SCF gate took {dt:.1f}s")
    else:
        say("  tiny SCF gate SKIPPED (--no-scf)")

    dt = time.time() - t0
    verdict = "PASS" if not fails else f"FAIL ({len(fails)}: {fails})"
    say(f"SMOKE {verdict} in {dt:.0f}s")
    return 0 if not fails else 1


def _parse_args(argv=None):
    p = argparse.ArgumentParser(
        description="chiral ansa-titanocene stereo-DDE job (issue #11)")
    p.add_argument("--smoke", action="store_true",
                   help="preflight: geometry+config validation, <=3 min, 1 core")
    p.add_argument("--no-scf", action="store_true",
                   help="with --smoke: skip even the tiny propene SCF gate")
    return p.parse_args(argv)


if __name__ == "__main__":
    _args = _parse_args()
    sys.exit(main_smoke(_args) if _args.smoke else main_real(_args))
