#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
calc/comos_hds.py — Track L6 OPENER (AWS job): the active phase of hydrotreating
(HDS) catalysts — a MoS2 edge-motif trimer and its Co / Ni promoted variants.

WHY THIS TRACK (APPLIED_WIDE_SCAN.md, L6). The HDS/HDN active phase of a refinery
hydrotreater is a sulfided MoS2 nanoribbon edge decorated with Co or Ni promoters
("CoMoS" / "NiMoS" phase). That is literally the same spin-entangled sulfide
chemistry as our [Fe4S4] cubanes and the Fe-Ni-S pentlandite motif
(calc/fe_ni_s_pentlandite.py), where the measured PBE-vs-PBE0 spin-gap spread was
~12 kcal/mol — i.e. exactly the regime where a single-functional DFT screen of
promoter recipes is not trustworthy. Hydroprocessing catalyst is a
sanctions-critical import position for Russian refineries; the deliverable of this
track is a quantum passport of the active phase (edge spin state, multireference
character, promoter ranking by the HDS activity descriptor) for domestic catalyst
makers.

THE CLUSTER (generated in this script, no external xyz needed).
  [M3S9]2- — a three-metal row cut from the MoS2 edge motif, built from idealized
  bond lengths: M-M 3.16 A (the MoS2 lattice constant a), M-S 2.41 A, S-S 2.05 A
  (the disulfide of the S-edge). Atom order is FIXED so that a promoter swap and
  the sulfur-vacancy truncation touch nothing else:
      0,1,2   metals  M0, M1(the promoted / edge site), M2   (row along x)
      3..6    four mu2-S bridges (two above, two below the metal row)
      7..10   two eta2-S2 disulfides capping the row ends (S-S = 2.05 A)
      11      the TERMINAL EDGE S on the central metal — the sulfur removed to
              make the coordinatively unsaturated site (CUS). It is LAST, so the
              vacancy cluster is just the atom list truncated by one.
  Compositions: `mo` = Mo3S9 (unpromoted), `co` = Mo2CoS9 (M1 -> Co),
  `ni` = Mo2NiS9 (M1 -> Ni). As in the pentlandite run, ONLY the nuclear charge at
  the swapped site changes; the framework geometry is FROZEN and the total charge
  is held at -2 for all three, so the three numbers are comparable to each other.
  Formal picture of the unpromoted cluster: 3 Mo(IV) + 4 mu2-S(2-) + 2 S2(2-) +
  1 terminal S(2-) => -2. No oxidation state is imposed on Co/Ni: the SCF decides,
  and the assumption is recorded in the JSON.

THE PROGRAM (each stage checkpoints to JSON; a re-launched box RESUMES).
  1. SPIN LADDER + BS/HS-SCF. UKS-PBE/def2-SVP. For each composition the electron
     parity is computed and the parity-correct multiplicities 2S = p, p+2, p+4 are
     converged (HS seed = the top rung); the lowest is the ground state — it is not
     postulated. Then a broken-symmetry Ms = p guess is built by flipping the
     alpha/beta AO blocks of the two OUTER metals against the central one; it is
     accepted only if it converges AND the Mulliken metal spins really carry the
     antiparallel pattern.
  2. UNO. Natural orbitals of the reference (BS if accepted, else the ground HS)
     total density: fractional occupations, their count inside [0.02, 1.98] and the
     Head-Gordon-style n_u = sum n(2-n). Sulfide edge states are strongly
     correlated — that is the thesis of this track, and n_u is the number that
     either supports it or does not.
  3. HDS DESCRIPTOR — sulfur vacancy formation energy:
         cluster + H2  ->  cluster(-S) + H2S
         dE_vac = E(cluster-S) + E(H2S) - E(cluster) - E(H2)
     the classical HDS activity descriptor (an easier vacancy = a more active
     catalyst; the volcano's left flank). The vacancy cluster gets its OWN spin
     ladder — a vacancy changes the metal count of unpaired electrons. The variant
     with 1/2 coefficients on the gas terms, as written in the track brief, is also
     reported (`dE_vac_half_coeff_kcal`) but it is NOT the balanced reaction and
     the balanced one is the descriptor of record.
  4. WHERE DFT LIES. On the best composition of this run: (a) the high-spin vs
     BS spin gap on BOTH PBE and PBE0 — the functional spread is the Fe-Ni-S
     argument (~12 kcal/mol there) transplanted to CoMoS; (b) dE_vac recomputed at
     PBE0, giving the SIGN and size of the PBE->PBE0 correction to the descriptor
     itself. If the promoter ranking is smaller than that correction, a
     single-functional screen cannot rank promoters — that is the honest verdict.

HONEST SCOPE (also written into the JSON `honesty` block).
  - A 12-atom gas-phase cluster proxy is NOT a real MoS2 nanoribbon on Al2O3:
    no support, no periodic edge, no H2/H2S atmosphere, no coverage equilibrium,
    no ZPE/thermal corrections, gas references at fixed literature geometries.
  - Single points on the idealized frozen framework (no relaxation): the numbers
    are comparable to EACH OTHER (same geometry, same charge, same basis), not to
    experimental vacancy energies.
  - def2-SVP, DF-accelerated, no dynamic correlation beyond the functional.
  - signal.alarm stage caps are best-effort; the shell `timeout` in the AWS wrapper
    is the hard stop. Every derived block carries the convergence flags it came from.

ENV KNOBS (all optional):
  COMOS_STAGE_TIMEOUT  seconds per stage (each SCF), default 5400
  COMOS_TOTAL_BUDGET   total wall-clock budget seconds, default 19800
  COMOS_ONLY           mo | co | ni — run a single composition (split launcher)
  COMOS_MF             pbe | pbe0 | both (default both) — functionals for stage 4
  COMOS_SKIP_STAGE4    1 = skip the PBE0 / spin-gap block entirely
  PYSCF_MAX_MEMORY     MB for pyscf (the AWS wrapper sets it from the box RAM)

OUTPUT (incremental, atomic writes — nothing lost on a kill; resume-aware):
  calc/comos_hds_results.json      rolling results
  calc/comos_<comp>.xyz            the geometry actually used per composition

Run (heavy; meant for AWS via calc/comos_aws.py):
    python3 -u calc/comos_hds.py
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
    from pyscf import gto, scf, dft
except ImportError as e:                                   # pyscf optional at lint
    sys.exit(f"[deps] pyscf is required:  python3 -m pip install pyscf   ({e})")

HERE = os.path.dirname(os.path.abspath(__file__))
RESULTS = os.path.join(HERE, "comos_hds_results.json")
KCAL = 627.50947406            # Ha -> kcal/mol
CM = 219474.6313705            # Ha -> cm^-1

BASIS = "def2-svp"
CHARGE = -2
XC_BASE = "pbe"

# --- idealized MoS2 edge-motif bond lengths (Angstrom) ----------------------
D_MM = 3.16      # metal-metal along the edge row = MoS2 lattice constant a
D_MS = 2.41      # metal-sulfur
D_SS = 2.05      # S-S of the eta2-disulfide (S-edge motif)
EDGE_S = 11      # index of the terminal edge S (LAST atom) -> removed for the CUS

# gas-phase references at fixed literature geometries (no relaxation, no ZPE)
H2_R = 0.741                       # A
H2S_R, H2S_ANG = 1.336, 92.1       # A, degrees

STAGE_TIMEOUT = int(os.environ.get("COMOS_STAGE_TIMEOUT", "5400"))
TOTAL_BUDGET = int(os.environ.get("COMOS_TOTAL_BUDGET", "19800"))
ONLY = os.environ.get("COMOS_ONLY", "").lower().strip()
MF_GAP = os.environ.get("COMOS_MF", "both").lower()
SKIP4 = os.environ.get("COMOS_SKIP_STAGE4", "0") == "1"

COMPOSITIONS = [("mo", None, "Mo3S9(2-) — unpromoted MoS2 edge trimer"),
                ("co", "Co", "Mo2CoS9(2-) — Co-promoted edge (CoMoS phase proxy)"),
                ("ni", "Ni", "Mo2NiS9(2-) — Ni-promoted edge (NiMoS phase proxy)")]

T0 = time.time()


# ------------------------------------------------------------------ utilities
def elapsed():
    return time.time() - T0


def remaining():
    return TOTAL_BUDGET - elapsed()


def save(res):
    """Atomic incremental save: write tmp, then os.replace."""
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
    """Best-effort per-stage wall-clock cap via SIGALRM (main thread only)."""
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
def build_cluster(promoter=None, vacancy=False):
    """Generate the [M3S9]2- MoS2-edge-motif trimer from idealized bond lengths.

    Row of three metals along x (M-M = 3.16 A). Each M-M pair is doubly bridged by
    mu2-S (one above, one below the row, M-S = 2.41 A); each end metal carries an
    eta2-S2 disulfide (S-S = 2.05 A, M-S = 2.41 A) in the metal plane; the central
    metal carries the terminal edge S (LAST atom) whose removal is the vacancy.

    Returns [(element, np.array(xyz)), ...]; `promoter` (e.g. "Co") replaces the
    CENTRAL metal only; `vacancy=True` drops the last atom (the edge S).
    """
    h_br = float(np.sqrt(D_MS ** 2 - (D_MM / 2.0) ** 2))      # 1.820 A
    half = D_SS / 2.0                                          # 1.025 A
    d_out = float(np.sqrt(D_MS ** 2 - half ** 2))              # 2.181 A

    m = [np.array([i * D_MM, 0.0, 0.0]) for i in range(3)]
    atoms = [("Mo", m[0]),
             (promoter or "Mo", m[1]),                         # the promoted site
             ("Mo", m[2])]
    # 3..6 : mu2-S bridges, above and below each M-M pair
    for xc in (D_MM / 2.0, 1.5 * D_MM):
        atoms.append(("S", np.array([xc, 0.0, +h_br])))
        atoms.append(("S", np.array([xc, 0.0, -h_br])))
    # 7..10 : eta2-S2 disulfides on the two end metals (in the metal plane)
    for sgn, mm in ((-1.0, m[0]), (+1.0, m[2])):
        for dy in (+half, -half):
            atoms.append(("S", mm + np.array([sgn * d_out, dy, 0.0])))
    # 11 : terminal edge S on the central metal — the vacancy site (LAST)
    atoms.append(("S", m[1] + np.array([0.0, D_MS, 0.0])))
    if vacancy:
        atoms = atoms[:EDGE_S]
    return atoms


def geometry_diagnostics(atoms):
    """Min pairwise distance + the metal coordination counts (overlap check)."""
    xyz = np.array([p for _e, p in atoms])
    els = [e for e, _p in atoms]
    n = len(xyz)
    dmin, pair = np.inf, None
    for i in range(n):
        for j in range(i + 1, n):
            d = float(np.linalg.norm(xyz[i] - xyz[j]))
            if d < dmin:
                dmin, pair = d, (i, j)
    coord = []
    for i, e in enumerate(els):
        if e in ("Mo", "Co", "Ni"):
            coord.append(int(sum(1 for j in range(n) if els[j] == "S"
                                 and np.linalg.norm(xyz[i] - xyz[j]) < 2.7)))
    return {"n_atoms": n, "formula": "".join(sorted(set(els))),
            "min_pair_dist": round(dmin, 3), "min_pair": list(pair),
            "metal_S_coordination": coord,
            "overlap_ok": bool(dmin > 1.9)}     # 2.05 = the intended disulfide


def write_xyz(atoms, path, comment):
    with open(path, "w") as fh:
        fh.write(f"{len(atoms)}\n{comment}\n")
        for e, p in atoms:
            fh.write(f"{e} {p[0]:.6f} {p[1]:.6f} {p[2]:.6f}\n")


def parity_of(atoms, charge=CHARGE):
    """Electron parity without building a mol (all ECP cores are even)."""
    z = sum(gto.charge(e) for e, _p in atoms)
    return int((z - charge) % 2)


def spin_candidates(atoms, charge=CHARGE, n=3):
    """The n lowest parity-correct multiplicities, 2S = p, p+2, ..."""
    p = parity_of(atoms, charge)
    return [p + 2 * k for k in range(n)]


# ----------------------------------------------------------------- mean field
def build_mol(atoms, spin, maxmem, charge=CHARGE):
    # ECP only for the elements that HAVE one in def2-SVP (Mo: ECP-28); passing
    # ecp=BASIS wholesale makes pyscf print an "ECP not found" line per light
    # element on every mol build — thousands of lines over a ladder of SCFs.
    ecp = {e: BASIS for e, _p in atoms if e in ("Mo",)}
    return gto.M(atom=[(e, tuple(p)) for e, p in atoms], basis=BASIS, ecp=ecp,
                 charge=charge, spin=spin, verbose=0, max_memory=maxmem)


def _new_mf(mol, xc, restricted=False):
    mf = (dft.RKS(mol) if restricted else dft.UKS(mol)).density_fit()
    mf.xc = xc
    return mf


def _converge(mf, dm0=None, level_shift=0.4, max_cycle=80, tag=""):
    """DIIS with a level shift; on stall, a Newton (SOSCF) restart — the robust
    recipe for hard metal-sulfide SCF cases (same as the Fe-Ni-S run)."""
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
    (standard broken-symmetry spin-flip construction)."""
    dma, dmb = dm[0].copy(), dm[1].copy()
    aoslice = mol.aoslice_by_atom()
    for a in flip_atoms:
        sl = slice(int(aoslice[a][2]), int(aoslice[a][3]))
        blk_a = dma[sl, sl].copy()
        dma[sl, sl] = dmb[sl, sl]
        dmb[sl, sl] = blk_a
    return np.array([dma, dmb])


def metal_spin_pops(mf, metals=(0, 1, 2)):
    """Mulliken spin population on each metal — the BS sanity check."""
    dma, dmb = mf.make_rdm1()
    m = (dma - dmb) @ mf.get_ovlp()
    aoslice = mf.mol.aoslice_by_atom()
    return [round(float(np.trace(m[slice(int(aoslice[a][2]),
                                         int(aoslice[a][3])),
                                 slice(int(aoslice[a][2]),
                                       int(aoslice[a][3]))])), 3)
            for a in metals]


def spin_ladder(atoms, maxmem, xc, cands, label, info, charge=CHARGE):
    """UKS ladder over parity-correct multiplicities; returns (mf_best, 2S_best).
    The ground spin state is MEASURED, not postulated."""
    rungs, best = [], None
    for s2 in cands:
        if remaining() < 600:
            rungs.append({"spin_2S": s2, "skipped": "time budget"})
            continue
        t0 = time.time()
        try:
            mol = build_mol(atoms, s2, maxmem, charge)
            mf = with_timeout(min(STAGE_TIMEOUT, remaining()), _converge,
                              _new_mf(mol, xc), tag=f"{label}-2S{s2}")
            ss, _ = mf.spin_square()
            rec = {"spin_2S": s2, "e_h": float(mf.e_tot),
                   "converged": bool(mf.converged),
                   "spin_square": round(float(ss), 3),
                   "seconds": round(time.time() - t0, 1)}
            if mf.converged and (best is None or mf.e_tot < best[0].e_tot):
                best = (mf, s2)
        except StageTimeout:
            rec = {"spin_2S": s2, "timed_out_after_s": int(STAGE_TIMEOUT)}
        except Exception as e:
            rec = {"spin_2S": s2, "error": f"{type(e).__name__}: {e}"}
        rungs.append(rec)
        print(f"[{label}] 2S={s2}: {rec}", flush=True)
    info["ladder"] = rungs
    ok = [r for r in rungs if r.get("converged")]
    if ok:
        e0 = min(r["e_h"] for r in ok)
        info["rel_kcal"] = {str(r["spin_2S"]): round((r["e_h"] - e0) * KCAL, 2)
                            for r in ok}
        info["ground_2S"] = best[1]
        info["e_ground_h"] = float(best[0].e_tot)
    else:
        info["ground_2S"] = None
    return best if best else (None, None)


def broken_symmetry(atoms, maxmem, xc, mf_hs, info, flip=(0, 2)):
    """BS Ms = parity guess: flip the alpha/beta blocks of the OUTER metals
    against the central (promoted) one. Accepted only if converged AND the metal
    Mulliken spins really carry the antiparallel pattern."""
    p = parity_of(atoms)
    try:
        mol_bs = build_mol(atoms, p, maxmem)
        dm0 = flip_dm(mf_hs.make_rdm1(), mf_hs.mol, flip)
        t0 = time.time()
        mf_bs = with_timeout(min(STAGE_TIMEOUT, remaining()), _converge,
                             _new_mf(mol_bs, xc), dm0=dm0, level_shift=0.3,
                             tag="BS")
        ss, _ = mf_bs.spin_square()
        pops = metal_spin_pops(mf_bs)
        # expected: outer metals opposite in sign to the central one, all polarized
        signs_ok = bool(pops[1] * pops[0] < 0 and pops[1] * pops[2] < 0)
        polarized = all(abs(x) > 0.3 for x in pops)
        info.update({"flip_atoms": list(flip), "ms_2S": p,
                     "e_h": float(mf_bs.e_tot),
                     "converged": bool(mf_bs.converged),
                     "spin_square": round(float(ss), 3),
                     "metal_mulliken_spin": pops,
                     "seconds": round(time.time() - t0, 1),
                     "pattern_ok": bool(mf_bs.converged and signs_ok
                                        and polarized)})
        print(f"[bs] E={mf_bs.e_tot:.5f} conv={mf_bs.converged} "
              f"<S^2>={ss:.2f} spins={pops} ok={info['pattern_ok']}", flush=True)
        return mf_bs if info["pattern_ok"] else None
    except StageTimeout:
        info["timed_out"] = True
    except Exception as e:
        info["error"] = f"{type(e).__name__}: {e}"
    return None


# --------------------------------------------------------------- UNO analysis
def unos(mf):
    """Natural orbitals of the total (alpha+beta) density: diagonalize
    S^1/2 D S^1/2; returns occupations in descending order."""
    s = mf.get_ovlp()
    shalf, sminus = _sqrtm_sym(s)
    dma, dmb = mf.make_rdm1()
    n, u = np.linalg.eigh(shalf @ (dma + dmb) @ shalf)
    order = np.argsort(n)[::-1]
    return np.clip(n[order], 0.0, 2.0), sminus @ u[:, order]


def uno_block(mf, source):
    """Fractional-occupation spectrum and n_u of the reference density."""
    occ, _C = unos(mf)
    mask = (occ > 0.02) & (occ < 1.98)
    n_u = float(np.sum(occ * (2.0 - occ)))
    return {"from_reference": source,
            "n_fractional_0.02_1.98": int(mask.sum()),
            "fractional_occs": [round(float(x), 3) for x in occ[mask][:40]],
            "n_u_head_gordon": round(n_u, 3),
            "note": "n_u = sum n(2-n) over UNO occupations of the reference "
                    "density; sulfide edge states are the target — a large n_u "
                    "means a single-determinant promoter screen is on thin ice"}


# ------------------------------------------------------- gas-phase references
def gas_refs(maxmem, xc, cache):
    """E(H2) and E(H2S) at fixed literature geometries (RKS, closed shell)."""
    key = xc
    if key in cache:
        return cache[key]
    a = np.radians(H2S_ANG / 2.0)
    h2 = [("H", np.array([0.0, 0.0, 0.0])), ("H", np.array([0.0, 0.0, H2_R]))]
    h2s = [("S", np.array([0.0, 0.0, 0.0])),
           ("H", np.array([H2S_R * np.sin(a), 0.0, H2S_R * np.cos(a)])),
           ("H", np.array([-H2S_R * np.sin(a), 0.0, H2S_R * np.cos(a)]))]
    out = {}
    for name, at in (("h2", h2), ("h2s", h2s)):
        mol = build_mol(at, 0, maxmem, charge=0)
        mf = _converge(_new_mf(mol, xc, restricted=True), tag=f"{name}-{xc}")
        out[name] = {"e_h": float(mf.e_tot), "converged": bool(mf.converged)}
        print(f"[gas] {name} ({xc}): E={mf.e_tot:.6f} conv={mf.converged}",
              flush=True)
    out["geometry_note"] = (f"fixed literature geometries: H2 r={H2_R} A; "
                            f"H2S r={H2S_R} A, angle={H2S_ANG} deg; "
                            "no relaxation, no ZPE / thermal corrections")
    cache[key] = out
    return out


def dE_vac(e_cluster, e_vac, refs):
    """cluster + H2 -> cluster(-S) + H2S  (the balanced HDS vacancy reaction)."""
    e_h2, e_h2s = refs["h2"]["e_h"], refs["h2s"]["e_h"]
    bal = (e_vac + e_h2s) - (e_cluster + e_h2)
    half = (e_vac + 0.5 * e_h2s) - (e_cluster + 0.5 * e_h2)
    return {"dE_vac_kcal": round(bal * KCAL, 2),
            "dE_vac_half_coeff_kcal": round(half * KCAL, 2),
            "reaction": "cluster + H2 -> cluster(-S) + H2S (balanced; the "
                        "descriptor of record)",
            "half_coeff_note": "dE = E(vac) + 1/2 E(H2S) - E(cluster) - 1/2 E(H2) "
                               "as written in the track brief — reported for "
                               "traceability, but it is NOT mass balanced"}


# --------------------------------------------------------------- composition
def run_composition(key, promoter, title, maxmem, res, gascache):
    """Stages 1-3 for one composition; returns its result block."""
    blk = res["compositions"].get(key, {})
    if blk.get("done"):
        print(f"[{key}] resume: already done", flush=True)
        return blk
    atoms = build_cluster(promoter)
    atoms_vac = build_cluster(promoter, vacancy=True)
    blk.update({"title": title, "promoter": promoter or "none (Mo only)",
                "charge": CHARGE,
                "geometry": geometry_diagnostics(atoms),
                "geometry_vacancy": geometry_diagnostics(atoms_vac),
                "removed_atom_index": EDGE_S,
                "spin_candidates": spin_candidates(atoms),
                "assumption": "promoter swap changes ONLY the nuclear charge at "
                              "the central metal site; framework geometry frozen, "
                              "total charge held at -2 for all compositions"})
    write_xyz(atoms, os.path.join(HERE, f"comos_{key}.xyz"),
              f"L6 {title} (idealized MoS2 edge motif, generated)")
    write_xyz(atoms_vac, os.path.join(HERE, f"comos_{key}_vac.xyz"),
              f"L6 {title} minus the terminal edge S (CUS)")
    res["compositions"][key] = blk
    save(res)
    print(f"\n=== [{key}] {title} — {blk['geometry']} ===", flush=True)

    # -- stage 1: spin ladder + BS -------------------------------------------
    if "spins" not in blk:
        blk["spins"] = {}
        mf_hs, s2 = spin_ladder(atoms, maxmem, XC_BASE, blk["spin_candidates"],
                                key, blk["spins"])
        save(res)
    else:                       # resumed: recompute the ground SCF (cheap-ish)
        s2 = blk["spins"].get("ground_2S")
        mf_hs = None
        if s2 is not None:
            mf_hs = _converge(_new_mf(build_mol(atoms, s2, maxmem), XC_BASE),
                              tag=f"{key}-redo")
    if mf_hs is None:
        blk["error"] = "no converged spin state for the parent cluster"
        blk["done"] = True
        save(res)
        return blk
    blk["e_parent_h"] = float(mf_hs.e_tot)
    blk["parent_2S"] = int(s2)

    if "bs" not in blk:
        blk["bs"] = {}
        mf_bs = broken_symmetry(atoms, maxmem, XC_BASE, mf_hs, blk["bs"])
        save(res)
    else:
        mf_bs = None
    ref_mf = mf_bs if mf_bs is not None else mf_hs
    ref_tag = "BS" if mf_bs is not None else f"ground UKS 2S={s2}"
    # the BS solution, if accepted and lower, IS the reference energy
    e_bs = (float(mf_bs.e_tot) if mf_bs is not None
            else (blk["bs"].get("e_h") if blk["bs"].get("pattern_ok") else None))
    if e_bs is not None and e_bs < blk["e_parent_h"]:
        blk["e_parent_h"] = float(e_bs)
        blk["parent_state_used"] = "BS Ms=parity (lower than the spin ladder)"
    else:
        blk["parent_state_used"] = f"UKS 2S={s2} (BS not accepted or higher)"

    # -- stage 2: UNO ---------------------------------------------------------
    if "uno" not in blk:
        blk["uno"] = uno_block(ref_mf, ref_tag)
        print(f"[{key}] UNO: {blk['uno']['n_fractional_0.02_1.98']} fractional, "
              f"n_u={blk['uno']['n_u_head_gordon']}", flush=True)
        save(res)

    # -- stage 3: sulfur vacancy ---------------------------------------------
    if "vacancy" not in blk:
        blk["vacancy"] = {"spins": {}}
        mf_v, s2v = spin_ladder(atoms_vac, maxmem, XC_BASE,
                                spin_candidates(atoms_vac), f"{key}-vac",
                                blk["vacancy"]["spins"])
        if mf_v is not None:
            blk["vacancy"]["e_h"] = float(mf_v.e_tot)
            blk["vacancy"]["2S"] = int(s2v)
        save(res)
    if blk["vacancy"].get("e_h") is not None:
        refs = gas_refs(maxmem, XC_BASE, gascache)
        res["gas_refs"] = gascache
        blk["hds_descriptor_pbe"] = dE_vac(blk["e_parent_h"],
                                           blk["vacancy"]["e_h"], refs)
        print(f"[{key}] dE_vac(PBE) = "
              f"{blk['hds_descriptor_pbe']['dE_vac_kcal']} kcal/mol", flush=True)
    else:
        blk["hds_descriptor_pbe"] = {"error": "vacancy cluster did not converge"}
    blk["done"] = True
    save(res)
    return blk


# ----------------------------------------------- stage 4: where DFT lies
def stage_where_dft_lies(key, promoter, maxmem, funcs, res, gascache):
    """On the best composition: (a) HS-vs-BS spin gap on PBE and PBE0 (the
    functional spread), (b) dE_vac recomputed at PBE0 -> the SIGN and size of the
    PBE->PBE0 correction to the HDS descriptor itself."""
    atoms = build_cluster(promoter)
    atoms_vac = build_cluster(promoter, vacancy=True)
    cands = spin_candidates(atoms)
    out = {"composition": key, "functionals": funcs, "gaps": {}, "dE_vac": {}}
    hi = cands[-1]                                   # top rung = the HS sector
    for xc in funcs:
        if remaining() < 1500:
            out["gaps"][xc] = {"skipped": "time budget"}
            continue
        try:
            print(f"[stage4] {xc.upper()} high-spin (2S={hi})", flush=True)
            t0 = time.time()
            mf_hi = with_timeout(min(STAGE_TIMEOUT, remaining()), _converge,
                                 _new_mf(build_mol(atoms, hi, maxmem), xc),
                                 tag=f"{xc}-HS")
            print(f"[stage4] {xc.upper()} low-spin (BS Ms=parity)", flush=True)
            p = parity_of(atoms)
            dm0 = flip_dm(mf_hi.make_rdm1(), mf_hi.mol, (0, 2))
            mf_lo = with_timeout(min(STAGE_TIMEOUT, remaining()), _converge,
                                 _new_mf(build_mol(atoms, p, maxmem), xc),
                                 dm0=dm0, level_shift=0.3, tag=f"{xc}-BS")
            ss_lo, _ = mf_lo.spin_square()
            gap = float(mf_hi.e_tot) - float(mf_lo.e_tot)
            out["gaps"][xc] = {
                "e_high_spin": float(mf_hi.e_tot), "high_spin_2S": hi,
                "high_spin_converged": bool(mf_hi.converged),
                "e_low_spin_bs": float(mf_lo.e_tot), "low_spin_ms_2S": p,
                "low_spin_converged": bool(mf_lo.converged),
                "low_spin_s_squared": round(float(ss_lo), 3),
                "gap_high_minus_low_kcal": round(gap * KCAL, 2),
                "gap_high_minus_low_cm": round(gap * CM, 1),
                "seconds": round(time.time() - t0, 1),
                "note": "the BS low-spin sector is spin-contaminated; this is the "
                        "DFT self-disagreement metric, not a spin-projected gap"}
            print(f"[stage4] {xc.upper()} gap = {gap*KCAL:.2f} kcal/mol",
                  flush=True)
            # dE_vac on the same functional
            mf_v = None
            for s2v in spin_candidates(atoms_vac):
                if remaining() < 900:
                    break
                try:
                    cand = with_timeout(min(STAGE_TIMEOUT, remaining()),
                                        _converge,
                                        _new_mf(build_mol(atoms_vac, s2v,
                                                          maxmem), xc),
                                        tag=f"{xc}-vac{s2v}")
                    if cand.converged and (mf_v is None
                                           or cand.e_tot < mf_v.e_tot):
                        mf_v = cand
                except Exception:
                    continue
            if mf_v is not None:
                refs = gas_refs(maxmem, xc, gascache)
                res["gas_refs"] = gascache
                e_par = min(float(mf_hi.e_tot), float(mf_lo.e_tot))
                out["dE_vac"][xc] = dE_vac(e_par, float(mf_v.e_tot), refs)
                out["dE_vac"][xc]["parent_e_h"] = e_par
                out["dE_vac"][xc]["vacancy_e_h"] = float(mf_v.e_tot)
        except StageTimeout:
            out["gaps"][xc] = {"timed_out": True}
        except Exception as e:
            out["gaps"][xc] = {"error": f"{type(e).__name__}: {e}"}
        res["where_dft_lies"] = out
        save(res)

    g = {x: out["gaps"][x].get("gap_high_minus_low_kcal") for x in out["gaps"]}
    if g.get("pbe") is not None and g.get("pbe0") is not None:
        out["spin_gap_functional_spread_kcal"] = round(g["pbe0"] - g["pbe"], 2)
    d = {x: out["dE_vac"].get(x, {}).get("dE_vac_kcal") for x in out["dE_vac"]}
    if d.get("pbe") is not None and d.get("pbe0") is not None:
        shift = round(d["pbe0"] - d["pbe"], 2)
        out["dE_vac_pbe0_minus_pbe_kcal"] = shift
        out["pbe_to_pbe0_sign"] = ("PBE0 makes the vacancy HARDER (+)" if shift > 0
                                   else "PBE0 makes the vacancy EASIER (-)")
    out["descriptor"] = (
        "PBE-vs-PBE0 spread on BOTH the spin gap and dE_vac = 'where DFT lies' "
        "for the CoMoS edge; the Fe-Ni-S argument (~12 kcal/mol spin-gap spread on "
        "the pentlandite motif) transplanted to the hydrotreating active phase")
    res["where_dft_lies"] = out
    save(res)
    return out


# ----------------------------------------------------------------- verdict
def make_verdict(res):
    rank = []
    for key, blk in res.get("compositions", {}).items():
        d = blk.get("hds_descriptor_pbe", {}).get("dE_vac_kcal")
        if d is not None:
            rank.append((key, d, blk.get("uno", {}).get("n_u_head_gordon")))
    rank.sort(key=lambda r: r[1])
    v = {"ranking_by_dE_vac_pbe_kcal": [{"composition": k, "dE_vac_kcal": d,
                                         "n_u": nu} for k, d, nu in rank],
         "convention": "lower dE_vac = an easier sulfur vacancy = the more active "
                       "HDS catalyst (left flank of the vacancy volcano)"}
    if len(rank) >= 2:
        v["spread_between_promoters_kcal"] = round(rank[-1][1] - rank[0][1], 2)
        v["best"] = rank[0][0]
    w = res.get("where_dft_lies", {})
    shift = w.get("dE_vac_pbe0_minus_pbe_kcal")
    if shift is not None:
        v["pbe_to_pbe0_shift_kcal"] = shift
        v["pbe_to_pbe0_sign"] = w.get("pbe_to_pbe0_sign")
        spread = v.get("spread_between_promoters_kcal")
        if spread is not None:
            v["ranking_survives_functional_change"] = bool(abs(spread) >
                                                           abs(shift))
            v["verdict"] = (
                "the promoter ranking is LARGER than the PBE->PBE0 correction — "
                "a DFT screen can order promoters at this level"
                if abs(spread) > abs(shift) else
                "the PBE->PBE0 correction is COMPARABLE TO OR LARGER than the "
                "spread between promoters — a single-functional DFT screen CANNOT "
                "rank Co vs Ni vs bare Mo here; the multireference track is not "
                "optional but required")
    elif rank:
        v["verdict"] = ("ranking computed at PBE only; the PBE->PBE0 correction "
                        "was not obtained in this run — no claim about whether "
                        "the ranking survives a functional change")
    else:
        v["verdict"] = "NO VERDICT: no composition produced a dE_vac"
    return v


# ------------------------------------------------------------------------ main
def main():
    prev = load_prev()
    resume = prev is not None and prev.get("meta", {}).get("cluster")
    res = prev if resume else {"status": "running"}
    res.setdefault("meta", {})
    res.setdefault("compositions", {})
    res["meta"].update({
        "cluster": "[M3S9]2- MoS2 edge-motif trimer (M1 = Mo | Co | Ni)",
        "purpose": "L6: паспорта активной фазы катализаторов гидроочистки "
                   "(CoMoS/NiMoS) — санкционное импортозамещение катализаторов "
                   "гидропроцессов НПЗ; ранжир промоторов по дескриптору HDS "
                   "(энергия серной вакансии) и проверка, переживает ли ранжир "
                   "смену функционала",
        "basis": BASIS, "charge": CHARGE, "functional_base": XC_BASE,
        "bond_lengths_A": {"M-M": D_MM, "M-S": D_MS, "S-S(disulfide)": D_SS},
        "atom_order": "0-2 metals (1 = promoted site), 3-6 mu2-S bridges, "
                      "7-10 eta2-S2 disulfides, 11 terminal edge S (vacancy site)",
        "geometry_note": "SINGLE POINTS on the idealized generated framework — no "
                         "relaxation; numbers are comparable to each other, not to "
                         "experimental vacancy energies",
        "only": ONLY or "all three compositions",
        "spin_gap_functionals": MF_GAP,
        "stage_timeout_s": STAGE_TIMEOUT, "total_budget_s": TOTAL_BUDGET,
        "resumed": bool(resume)})
    res["honesty"] = [
        "кластер-прокси (12 атомов, газовая фаза) != реальная нанолента MoS2 на "
        "Al2O3: нет носителя, нет периодического края, нет H2/H2S-среды и "
        "равновесного покрытия серой",
        "геометрия идеализированная и замороженная (без релаксации), без ZPE и "
        "температурных поправок; газовые реперы H2/H2S — при фиксированных "
        "литературных геометриях",
        "промотор посажен подменой заряда ядра при фиксированном полном заряде -2: "
        "степень окисления Co/Ni НЕ навязана, её выбирает SCF",
        "def2-SVP, DF-ускорение, динамическая корреляция только на уровне "
        "функционала (CASSCF здесь не считается — только UNO-диагностика n_u)",
        "BS-решение спин-контаминировано; спин-щель здесь — метрика "
        "самонесогласия DFT, а не спин-проецированное число",
        "числа сравнимы МЕЖДУ СОБОЙ (одна геометрия, один заряд, один базис) — "
        "абсолютные значения дескриптора не переносятся на эксперимент"]
    save(res)

    try:
        maxmem = int(os.environ.get("PYSCF_MAX_MEMORY", "16000"))
        gascache = res.get("gas_refs", {})
        comps = [c for c in COMPOSITIONS if not ONLY or c[0] == ONLY]
        if ONLY and not comps:
            raise ValueError(f"COMOS_ONLY={ONLY!r} is not one of mo|co|ni")
        print(f"[run] compositions={[c[0] for c in comps]} maxmem={maxmem}MB "
              f"budget={TOTAL_BUDGET}s", flush=True)

        for key, promoter, title in comps:
            if remaining() < 1200:
                res["compositions"].setdefault(key, {})["skipped"] = "time budget"
                save(res)
                continue
            run_composition(key, promoter, title, maxmem, res, gascache)

        if not SKIP4 and remaining() > 2400:
            done = [(k, b) for k, b in res["compositions"].items()
                    if b.get("hds_descriptor_pbe", {}).get("dE_vac_kcal")
                    is not None]
            if done:
                best = min(done,
                           key=lambda kb: kb[1]["hds_descriptor_pbe"]
                           ["dE_vac_kcal"])[0]
                promoter = dict((k, p) for k, p, _t in COMPOSITIONS)[best]
                funcs = {"pbe": ["pbe"], "pbe0": ["pbe0"],
                         "both": ["pbe", "pbe0"]}.get(MF_GAP, ["pbe", "pbe0"])
                print(f"\n=== [stage4] where DFT lies, on '{best}' ===", flush=True)
                stage_where_dft_lies(best, promoter, maxmem, funcs, res, gascache)
        else:
            res["where_dft_lies"] = {"skipped": "COMOS_SKIP_STAGE4 or time budget"}

        res["verdict"] = make_verdict(res)
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
