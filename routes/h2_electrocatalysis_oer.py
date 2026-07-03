#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
OER dG-ladder on a minimal Ni(OH)2 cluster -- RELAXED geometries (CHE).
Track: h2-electrocatalysis (MatterForge / OCM research diary).

This is the relaxed-geometry successor of the preliminary single-point ladder
from the 2026-06-08 diary entry (UKS-PBE/def2-SVP single points on
[Ni(OH)2(X)] gave eta ~ 1.40 V with an unphysical negative last step).
The original track scripts were lost; this file reconstructs the pipeline
so that it is construction-compatible with what the diary entries describe:

  Model / method (as in the diary entries)
  ----------------------------------------
  * Cluster: [Ni(OH)2] as the bare site "*"; adsorbates X = OH, O, OOH bound
    on-top of Ni -> states *, OH*, O*, OOH*.  Gas references: H2O, H2.
  * UKS-PBE / def2-SVP with density fitting (auxiliary basis chosen
    automatically by PySCF).  SCF robustness: static level shift, with a
    second-order (SOSCF / newton) fallback -- the same level_shift + SOSCF
    recipe that made all 6 states converge in the 2026-06-08 entry.
  * NEW vs. the preliminary entry: geometries are relaxed with
    pyscf + geomeTRIC (pip package "geometric").
  * 2-3 spin multiplicities per state; the lowest converged energy is used.

  Core constraints during relaxation (documented deliberately)
  -------------------------------------------------------------
  In a real NiOOH electrode the Ni active site and its bridging O atoms are
  held in place by the lattice.  A free 5-9 atom gas-phase cluster would
  instead relax into unrelated gas-phase minima (fragment rearrangement,
  adsorbate migration to O, etc.), destroying comparability between the
  ladder states.  We therefore FREEZE the cartesian positions of the
  Ni(OH)2 "core" heavy atoms -- Ni and both hydroxide O (geomeTRIC
  "$freeze / xyz 1,2,4", 1-based) -- and relax the hydroxide H atoms and
  ALL adsorbate atoms.  This is a crude lattice mimic and a known
  limitation; it is recorded in the results JSON.

  O2 reference (documented deliberately)
  --------------------------------------
  We do NOT compute triplet O2 with DFT directly: GGAs (PBE included) have a
  notorious ~0.3-0.7 eV error in the O2 binding energy.  Standard CHE
  practice instead closes the cycle with the experimental free energy of
  water splitting, 2 H2O(l) -> O2(g) + 2 H2(g): dG = 4.92 eV.  Hence
      dG4 = 4.92 eV - (dG1 + dG2 + dG3),
  equivalent to G(O2) := 2 G(H2O) - 2 G(H2) + 4.92 eV.  The preliminary
  diary ladder used the same closure (its four steps sum to ~4.92 eV).

  ZPE / TS corrections (LITERATURE VALUES, not computed here)
  -----------------------------------------------------------
  Hessians for every state x multiplicity are too expensive for this stage.
  We therefore use the standard aggregate CHE corrections
  (dZPE - T*dS, 298.15 K) applied to the adsorbate formation energies,
  from the Norskov-school CHE literature (Norskov et al., J. Phys. Chem. B
  2004, 108, 17886; as compiled and reused in e.g. Man et al., ChemCatChem
  2011, 3, 1159 and Bajdich et al., JACS 2013, 135, 13521):
      dG(OH*)  = dE(OH*)  + 0.35 eV
      dG(O*)   = dE(O*)   + 0.05 eV
      dG(OOH*) = dE(OOH*) + 0.40 eV
  These aggregates already include the ZPE/TS of the H2O/H2 references.
  The results JSON marks them explicitly as literature values
  ("computed_here": false) -- they are NOT our own frequencies.

  CHE ladder (U = 0 V, pH 0, T = 298.15 K)
  ----------------------------------------
      dG1 = G(OH*) - G(*) - G(H2O) + 1/2 G(H2)
      dG2 = G(O*)  - G(OH*)        + 1/2 G(H2)
      dG3 = G(OOH*) - G(O*) - G(H2O) + 1/2 G(H2)
      dG4 = 4.92 eV - dG1 - dG2 - dG3
      eta = max(dG_i)/e - 1.23 V

  Operational features
  --------------------
  * Incremental, atomic save into routes/h2_electrocatalysis_oer_results.json
    after every completed (state, multiplicity) calculation.
  * Restart: already-converged (state, multiplicity) records are skipped.
  * Threads taken from OMP_NUM_THREADS; progress printed with flush=True.
  * --smoke: tiny pipeline check (H2O single point + minimal H2 relaxation,
    def2-SVP), writes to a separate *_smoke.json, ~1-3 min on 1 core.

Usage:
    OMP_NUM_THREADS=1 python3 routes/h2_electrocatalysis_oer.py --smoke
    OMP_NUM_THREADS=4 python3 routes/h2_electrocatalysis_oer.py

Honesty rule of the diary: every number in the results JSON comes from a
real calculation performed by this script, except the thermal corrections
and the 4.92 eV closure, which are explicitly flagged as literature values.
"""

import argparse
import json
import os
import sys
import tempfile
import time

HARTREE_EV = 27.211386245988
WATER_SPLITTING_EV = 4.92          # experimental dG, 2H2O(l) -> O2(g) + 2H2(g)
EQUILIBRIUM_POTENTIAL_V = 1.23     # ideal OER catalyst step height

# Literature aggregate (dZPE - T*dS) corrections, eV -- see module docstring.
LIT_CORR_EV = {"OH*": 0.35, "O*": 0.05, "OOH*": 0.40}

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
RESULTS_PATH = os.path.join(SCRIPT_DIR, "h2_electrocatalysis_oer_results.json")
SMOKE_RESULTS_PATH = os.path.join(
    SCRIPT_DIR, "h2_electrocatalysis_oer_results_smoke.json")

MAX_MEMORY_MB = float(os.environ.get("PYSCF_MAX_MEMORY", 4000))

# --------------------------------------------------------------------------
# Geometries (Angstrom).  Initial guesses only -- relaxed by geomeTRIC.
# Cluster atom order is fixed: Ni(1) O(2) H(3) O(4) H(5) [adsorbate...],
# so the frozen-core selection "xyz 1,2,4" (Ni + both hydroxide O) is the
# same for every cluster state.
# --------------------------------------------------------------------------
CORE_XYZ = """\
Ni  0.0000   0.0000   0.0000
O   1.8500   0.0000   0.0000
H   2.4000   0.8000   0.0000
O  -1.8500   0.0000   0.0000
H  -2.4000  -0.8000   0.0000
"""

GEOMS = {
    "*":    CORE_XYZ,
    "OH*":  CORE_XYZ + "O   0.0000   0.0000   1.8500\n"
                       "H   0.0000   0.8000   2.4000\n",
    "O*":   CORE_XYZ + "O   0.0000   0.0000   1.6500\n",
    "OOH*": CORE_XYZ + "O   0.0000   0.0000   1.9000\n"
                       "O   1.1000   0.0000   2.8000\n"
                       "H   1.9500   0.2500   2.5500\n",
    "H2O":  "O 0.0000  0.0000  0.1173\n"
            "H 0.0000  0.7572 -0.4692\n"
            "H 0.0000 -0.7572 -0.4692\n",
    "H2":   "H 0.0 0.0 0.0\nH 0.0 0.0 0.7414\n",
}

# Multiplicities (2S+1) to try per state; lowest converged energy wins.
# Electron counts: Ni(OH)2 46e (even), +OH 55e (odd), +O 54e (even),
# +OOH 63e (odd).  Ni(II) d8 -> triplet usually lowest for "*".
STATES = [
    # (name, multiplicities, frozen_core?, part of smoke?)
    {"name": "H2O",  "mults": [1],       "freeze": False},
    {"name": "H2",   "mults": [1],       "freeze": False},
    {"name": "*",    "mults": [3, 1, 5], "freeze": True},
    {"name": "OH*",  "mults": [2, 4],    "freeze": True},
    {"name": "O*",   "mults": [3, 1, 5], "freeze": True},
    {"name": "OOH*", "mults": [2, 4],    "freeze": True},
]

FREEZE_SPEC = "xyz 1,2,4"   # geomeTRIC $freeze selection, 1-based: Ni, O, O

_CONSTRAINTS_FILE = None


def constraints_file():
    """Write the geomeTRIC $freeze constraints file once, return its path."""
    global _CONSTRAINTS_FILE
    if _CONSTRAINTS_FILE is None:
        fd, path = tempfile.mkstemp(prefix="oer_freeze_", suffix=".txt")
        with os.fdopen(fd, "w") as fh:
            fh.write("$freeze\n%s\n" % FREEZE_SPEC)
        _CONSTRAINTS_FILE = path
    return _CONSTRAINTS_FILE


def atomic_write_json(path, data):
    """Atomic incremental save: write tmp file in same dir, then os.replace."""
    tmp = path + ".tmp"
    with open(tmp, "w") as fh:
        json.dump(data, fh, indent=1, ensure_ascii=False)
        fh.write("\n")
    os.replace(tmp, path)


def load_results(path):
    if os.path.exists(path):
        with open(path) as fh:
            return json.load(fh)
    return {}


def setup_threads():
    from pyscf import lib
    nt = os.environ.get("OMP_NUM_THREADS")
    if nt:
        lib.num_threads(int(nt))
    print("[setup] pyscf threads = %d" % lib.num_threads(), flush=True)


def build_mol(xyz, spin, charge=0):
    from pyscf import gto
    return gto.M(atom=xyz, unit="Angstrom", basis="def2-svp",
                 charge=charge, spin=spin, symmetry=False,
                 verbose=0, max_memory=MAX_MEMORY_MB)


def make_mf(mol):
    """DF-UKS PBE with static level shift (diary recipe, part 1)."""
    from pyscf import dft
    mf = dft.UKS(mol)
    mf.xc = "PBE"
    mf = mf.density_fit()
    mf.level_shift = 0.25
    mf.max_cycle = 200
    mf.conv_tol = 1e-9
    return mf


def run_scf(mol):
    """Robust SCF: level-shifted DF-UKS, SOSCF(newton) fallback (diary
    recipe, part 2).  Returns (mf, energy_Ha, converged, <S^2>)."""
    mf = make_mf(mol)
    e = mf.kernel()
    used_soscf = False
    if not mf.converged:
        nmf = mf.newton()
        nmf.max_cycle = 100
        e = nmf.kernel(mf.mo_coeff, mf.mo_occ)
        mf = nmf
        used_soscf = True
    ss = float(mf.spin_square()[0])
    return mf, float(e), bool(mf.converged), ss, used_soscf


def relax_geometry(mf, freeze, maxsteps):
    """Relax with geomeTRIC through pyscf.geomopt.geometric_solver.
    Returns (mol_eq_or_None, converged, error_string_or_None)."""
    from pyscf.geomopt import geometric_solver
    kwargs = {"maxsteps": maxsteps}
    if freeze:
        kwargs["constraints"] = constraints_file()
    try:
        # v2: optimize through a SOSCF(newton) solver. The plain level-shifted
        # DF-UKS scanner fails to converge at displaced geometries
        # ("Nuclear gradients ... not converged" -> the whole optimization
        # died and v1 silently fell back to initial-geometry single points).
        solver = mf.newton()
        solver.max_cycle = 100
        mol_eq = geometric_solver.optimize(solver, assert_convergence=True,
                                           **kwargs)
        return mol_eq, True, None
    except Exception as exc:                      # noqa: BLE001 -- recorded
        return None, False, "%s: %s" % (type(exc).__name__, exc)


def mol_to_xyz(mol):
    from pyscf import lib
    coords = mol.atom_coords() * lib.param.BOHR   # Bohr -> Angstrom
    lines = []
    for i in range(mol.natm):
        x, y, z = coords[i]
        lines.append("%-2s %12.6f %12.6f %12.6f"
                     % (mol.atom_symbol(i), x, y, z))
    return "\n".join(lines)


def compute_state(name, mult, freeze, relax, maxsteps):
    """One (state, multiplicity) job: [relax ->] final robust SCF.
    Returns a JSON-ready record.  All numbers are real calculation output."""
    t0 = time.time()
    spin = mult - 1
    rec = {"state": name, "multiplicity": mult, "spin_2S": spin,
           "charge": 0, "relaxed": bool(relax),
           "frozen_core_atoms": FREEZE_SPEC if (freeze and relax) else None}

    mol = build_mol(GEOMS[name], spin)
    rec["n_atoms"] = mol.natm
    rec["n_basis"] = mol.nao_nr()

    opt_err = None
    if relax:
        # Pre-converge SCF at the initial geometry so the optimizer starts
        # from good orbitals; if plain SCF fails, seed it with SOSCF orbitals.
        mf0 = make_mf(mol)
        mf0.kernel()
        if not mf0.converged:
            nmf = mf0.newton()
            nmf.max_cycle = 100
            nmf.kernel(mf0.mo_coeff, mf0.mo_occ)
            if nmf.converged:
                mf0.mo_coeff, mf0.mo_occ = nmf.mo_coeff, nmf.mo_occ
                mf0.mo_energy = nmf.mo_energy
        rec["init_scf_converged"] = bool(mf0.converged)

        mol_eq, opt_ok, opt_err = relax_geometry(mf0, freeze, maxsteps)
        rec["opt_converged"] = opt_ok
        if opt_err:
            rec["opt_error"] = opt_err
        if mol_eq is None:
            # Honest fallback: keep the initial geometry, flag it clearly.
            mol_eq = mol
            rec["relaxed"] = False
            rec["warning"] = ("geometry optimization failed; energy below is "
                              "a single point at the INITIAL geometry")
    else:
        mol_eq = mol
        rec["opt_converged"] = None

    _, e, conv, ss, used_soscf = run_scf(mol_eq)
    s_ideal = spin / 2.0
    rec.update({
        "e_final_ha": e,
        "scf_converged": conv,
        "used_soscf": used_soscf,
        "s2": round(ss, 4),
        "s2_ideal": round(s_ideal * (s_ideal + 1.0), 4),
        "xyz_final": mol_to_xyz(mol_eq),
        "wall_s": round(time.time() - t0, 1),
    })
    return rec


def best_energy(state_block):
    """Lowest converged energy over multiplicities; None if nothing usable."""
    best = None
    for mult, rec in state_block.get("mults", {}).items():
        if not rec.get("scf_converged"):
            continue
        if best is None or rec["e_final_ha"] < best["e_final_ha"]:
            best = {"multiplicity": int(mult), "e_final_ha": rec["e_final_ha"],
                    "s2": rec.get("s2"), "relaxed": rec.get("relaxed"),
                    "opt_converged": rec.get("opt_converged")}
    return best


def compute_ladder(results):
    """CHE ladder from the best energies.  Returns ladder dict or None."""
    need = ["*", "OH*", "O*", "OOH*", "H2O", "H2"]
    E = {}
    warnings = []
    for name in need:
        blk = results.get("states", {}).get(name)
        best = best_energy(blk) if blk else None
        if best is None:
            return None, ["state %s missing or not converged" % name]
        E[name] = best["e_final_ha"] * HARTREE_EV
        if best.get("relaxed") is False or best.get("opt_converged") is False:
            warnings.append("state %s: geometry NOT (fully) relaxed -- "
                            "treat ladder as preliminary" % name)

    # Electronic formation energies (eV), then literature dZPE-T*dS.
    dE_OH = E["OH*"] + 0.5 * E["H2"] - E["*"] - E["H2O"]
    dE_O = E["O*"] + E["H2"] - E["*"] - E["H2O"]
    dE_OOH = E["OOH*"] + 1.5 * E["H2"] - E["*"] - 2.0 * E["H2O"]
    dG_OH = dE_OH + LIT_CORR_EV["OH*"]
    dG_O = dE_O + LIT_CORR_EV["O*"]
    dG_OOH = dE_OOH + LIT_CORR_EV["OOH*"]

    dg = [dG_OH,
          dG_O - dG_OH,
          dG_OOH - dG_O,
          WATER_SPLITTING_EV - dG_OOH]
    labels = ["* + H2O -> OH* + (H+ + e-)",
              "OH* -> O* + (H+ + e-)",
              "O* + H2O -> OOH* + (H+ + e-)",
              "OOH* -> * + O2 + (H+ + e-)"]
    imax = max(range(4), key=lambda i: dg[i])
    ladder = {
        "conditions": {"U_V": 0.0, "pH": 0, "T_K": 298.15},
        "formation_dE_eV": {"OH*": round(dE_OH, 3), "O*": round(dE_O, 3),
                            "OOH*": round(dE_OOH, 3)},
        "formation_dG_eV": {"OH*": round(dG_OH, 3), "O*": round(dG_O, 3),
                            "OOH*": round(dG_OOH, 3)},
        "steps_eV": {"dG1": round(dg[0], 3), "dG2": round(dg[1], 3),
                     "dG3": round(dg[2], 3), "dG4": round(dg[3], 3)},
        "step_labels": labels,
        "limiting_step": labels[imax],
        "eta_V": round(max(dg) - EQUILIBRIUM_POTENTIAL_V, 3),
        "dG4_closure": "dG4 = 4.92 eV - (dG1+dG2+dG3); experimental water "
                       "splitting dG, avoids direct triplet-O2 DFT",
        "best_states": {n: best_energy(results["states"][n]) for n in need},
    }
    return ladder, warnings


def make_meta(smoke):
    import pyscf
    try:
        import geometric
        geometric_ver = getattr(geometric, "__version__", "unknown")
    except ImportError:
        geometric_ver = "not installed"
    return {
        "track": "h2-electrocatalysis",
        "what": ("Relaxed-geometry OER dG-ladder (CHE) on minimal "
                 "[Ni(OH)2(X)] cluster; successor of the 2026-06-08 "
                 "single-point ladder (eta~1.40 V, unrelaxed)"),
        "smoke": bool(smoke),
        "method": {
            "scf": "UKS-PBE/def2-SVP, density fitting (auto aux basis)",
            "scf_robustness": "level_shift=0.25 + SOSCF(newton) fallback "
                              "(same recipe as the 2026-06-08 diary entry)",
            "geometry": "geomeTRIC via pyscf.geomopt.geometric_solver; "
                        "cluster core FROZEN: Ni + both hydroxide O "
                        "(geomeTRIC '$freeze " + FREEZE_SPEC + "', 1-based); "
                        "hydroxide H and all adsorbate atoms relaxed; "
                        "H2O/H2 references fully relaxed",
            "multiplicities": "2-3 per state, lowest converged energy taken",
        },
        "o2_reference": {
            "how": "G(O2) := 2 G(H2O) - 2 G(H2) + 4.92 eV (experimental "
                   "water-splitting dG); equivalently dG4 = 4.92 - sum(dG1..3)",
            "why": "PBE triplet O2 has a notorious ~0.3-0.7 eV binding error; "
                   "direct O2 DFT deliberately avoided (standard CHE practice)",
        },
        "thermal_corrections": {
            "computed_here": False,
            "source": "literature",
            "citation": "Norskov et al. J. Phys. Chem. B 2004, 108, 17886 "
                        "(CHE); aggregate values as reused in Man et al. "
                        "ChemCatChem 2011, 3, 1159 and Bajdich et al. JACS "
                        "2013, 135, 13521",
            "values_eV_dZPE_minus_TdS": LIT_CORR_EV,
            "note": "Applied to adsorbate formation energies; gas-reference "
                    "ZPE/TS already included in these aggregates. Hessians "
                    "were NOT computed for this cluster -- these numbers are "
                    "not ours.",
        },
        "limitations": [
            "minimal cluster: no extended NiOOH lattice, Ni oxidation state "
            "not pinned to III/IV",
            "frozen Ni/O core is a crude lattice mimic",
            "single-determinant DFT; even-electron states are plain UKS "
            "(no explicit broken-symmetry guess)",
            "def2-SVP is a small basis; no dispersion correction (matches "
            "the preliminary diary entry for comparability)",
        ],
        "versions": {"pyscf": pyscf.__version__, "geometric": geometric_ver},
    }


def main():
    ap = argparse.ArgumentParser(
        description="Relaxed OER dG-ladder (CHE) on [Ni(OH)2(X)], "
                    "UKS-PBE/def2-SVP+DF, pyscf+geomeTRIC")
    ap.add_argument("--smoke", action="store_true",
                    help="tiny pipeline check: H2O single point + minimal H2 "
                         "relaxation; separate results file; ~1-3 min, 1 core")
    ap.add_argument("--maxsteps", type=int, default=80,
                    help="geomeTRIC max optimization steps (default 80)")
    ap.add_argument("--states", type=str, default=None,
                    help="comma list to restrict states, e.g. '*,OH*' "
                         "(default: all)")
    ap.add_argument("--results", type=str, default=None,
                    help="override results JSON path")
    args = ap.parse_args()

    setup_threads()

    if args.smoke:
        results_path = args.results or SMOKE_RESULTS_PATH
        # Smoke: H2O single point (no relax) + H2 mini-relaxation, so that
        # SCF, DF, gradients, geomeTRIC, atomic save and restart all get
        # exercised once, cheaply.
        plan = [("H2O", [1], False, False, 0),
                ("H2", [1], False, True, 10)]
    else:
        results_path = args.results or RESULTS_PATH
        wanted = None
        if args.states:
            wanted = {s.strip() for s in args.states.split(",")}
        plan = [(st["name"], st["mults"], st["freeze"], True, args.maxsteps)
                for st in STATES
                if wanted is None or st["name"] in wanted]

    results = load_results(results_path)
    results["meta"] = make_meta(args.smoke)
    results.setdefault("states", {})
    results["meta"]["updated"] = time.strftime("%Y-%m-%dT%H:%M:%SZ",
                                               time.gmtime())
    atomic_write_json(results_path, results)

    njobs = sum(len(m) for _, m, _, _, _ in plan)
    k = 0
    for name, mults, freeze, relax, maxsteps in plan:
        blk = results["states"].setdefault(name, {"mults": {}})
        for mult in mults:
            k += 1
            done = blk["mults"].get(str(mult))
            if done and done.get("scf_converged") \
                    and done.get("relaxed") == relax:
                print("[%d/%d] %-4s mult=%d  SKIP (already converged, "
                      "E=%.6f Ha)" % (k, njobs, name, mult,
                                      done["e_final_ha"]), flush=True)
                continue
            print("[%d/%d] %-4s mult=%d  relax=%s ..."
                  % (k, njobs, name, mult, relax), flush=True)
            try:
                rec = compute_state(name, mult, freeze, relax, maxsteps)
            except Exception as exc:              # noqa: BLE001 -- recorded
                rec = {"state": name, "multiplicity": mult,
                       "scf_converged": False, "relaxed": False,
                       "error": "%s: %s" % (type(exc).__name__, exc)}
            blk["mults"][str(mult)] = rec
            blk["best"] = best_energy(blk)
            results["meta"]["updated"] = time.strftime(
                "%Y-%m-%dT%H:%M:%SZ", time.gmtime())
            atomic_write_json(results_path, results)
            if rec.get("scf_converged"):
                print("    -> E = %.6f Ha  <S^2>=%.3f (ideal %.3f)  "
                      "opt_conv=%s  %ds"
                      % (rec["e_final_ha"], rec["s2"], rec["s2_ideal"],
                         rec.get("opt_converged"), rec["wall_s"]),
                      flush=True)
            else:
                print("    -> FAILED: %s"
                      % rec.get("error", rec.get("opt_error", "SCF not "
                                                 "converged")), flush=True)

    if args.smoke:
        # Loose sanity windows only (not published numbers): PBE/def2-SVP
        # totals should land near textbook ballpark for H2O and H2.
        ok = True
        for name, (lo, hi) in {"H2O": (-76.8, -75.8),
                               "H2": (-1.30, -1.00)}.items():
            best = best_energy(results["states"].get(name, {}))
            if best is None:
                print("[smoke] %s: FAILED (no converged energy)" % name,
                      flush=True)
                ok = False
                continue
            e = best["e_final_ha"]
            tag = "OK" if lo < e < hi else "OUT OF SANITY WINDOW"
            if tag != "OK":
                ok = False
            print("[smoke] %-3s E = %.6f Ha  [%s]" % (name, e, tag),
                  flush=True)
        print("[smoke] pipeline %s; results in %s"
              % ("PASSED" if ok else "FAILED", results_path), flush=True)
        sys.exit(0 if ok else 1)

    ladder, warnings = compute_ladder(results)
    if ladder is None:
        print("[ladder] incomplete: %s" % "; ".join(warnings), flush=True)
    else:
        ladder["warnings"] = warnings
        results["ladder"] = ladder
        atomic_write_json(results_path, results)
        print("\n[ladder] CHE OER steps (eV), U=0, pH=0:", flush=True)
        for i, lab in enumerate(ladder["step_labels"], 1):
            print("  dG%d = %+6.3f   %s" % (i, ladder["steps_eV"]["dG%d" % i],
                                            lab), flush=True)
        print("  limiting: %s" % ladder["limiting_step"], flush=True)
        print("  eta = %.3f V  (max step - 1.23)" % ladder["eta_V"],
              flush=True)
        for w in warnings:
            print("  WARNING: %s" % w, flush=True)
    print("[done] results: %s" % results_path, flush=True)


if __name__ == "__main__":
    main()
