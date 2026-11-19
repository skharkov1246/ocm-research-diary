#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
calc/fes_cogef.py — Track T3 GATE (AWS job): COGEF stretch of a terminal Fe–S
(thiolate) bond in a [Fe4S4(SH)4]2- cubane — the computational rotation gate for
"tailings as a mine in reverse" (mechanochemical opening of the sulfide lattice;
see APPLIED_UNEXPECTED_TWISTS.md, section T3).

WHY THIS JOB. The COGEF pipeline (routes/cc_cogef_probe.py) turned "break the
C–C of polyolefins with FORCE instead of 400 C pyrolysis" into a number:
F_max = 5.58 nN for the central C–C of n-butane, room-temperature kinetics at
4.5–5 nN. T3 carries that idea to Fe–Ni–S ores: a mill is not a grinder but a
selective-cleavage REACTOR that opens the sulfide lattice under leaching. Before
any mineral-lattice model, gate #1 is the cheapest honest question: does an
Fe–S bond break EASIER or HARDER than a C–C bond under mechanical force? This
script answers that on the minimal Fe–S motif we already trust — the idealized
[Fe4S4(SH)4]2- cubane of the FeMoco diary (Stage 17 geometry) — by running the
same COGEF machinery on ONE terminal Fe–S(thiolate) bond.

THE MODEL.
  Cluster : [Fe4S4(SH)4]2- idealized cubane (calc/femoco_4fe4s_cubane.xyz),
            charge -2, def2-SVP + density fitting, UKS-PBE.
  Reference spin: HS (2S=18: 2 Fe3+ d5 + 2 Fe2+ d6, all aligned) BY DEFAULT.
            This is a deliberate, honestly-flagged choice. The physical ground
            state of [Fe4S4]2- is a broken-symmetry Ms=0 state, but BS is fragile
            across a bond-stretch scan (the spin pattern can re-collapse as the
            terminal S is pulled off, and re-flipping per point is expensive).
            The single-determinant HS surface converges robustly at every stretch
            point, and the COGEF force F = dE/dd along a CONSISTENT surface is a
            valid first screen: the terminal Fe–S(thiolate) sigma bond stiffness
            dominates the curvature, not the intra-cluster magnetic coupling.
            FES_REF=bs switches to a BS Ms=0 reference (HS seed -> flip Fe2,Fe3,
            per-point) if you want the physical-ground-state surface; every point
            records which reference it used.
  Coordinate: stretch ONE terminal Fe0–S8(H) bond — atoms (0, 8) — the least-
            bound motif of the cluster (the "loose" thiolate, analog of prying a
            sulfur out of the lattice), from equilibrium d0 (~2.25 A) to
            d0 + FES_DMAX in FES_STEP steps. The rest of the geometry is
            constrained-relaxed (geomeTRIC $freeze distance 1 9). If a relaxed
            step costs more than FES_STAGE_TIMEOUT, the point falls back to a
            RIGID displacement (move only S8 and its H12 along the bond vector,
            everything else frozen) and is flagged rigid=true — the same honest
            fallback used in routes/spin_field_oxidation.py (relaxed/fallback
            bookkeeping). FES_RIGID=1 forces the rigid scan for all points.

PER POINT: E, gradient audit (max |nuclear force| residual, Ha/Bohr), Mulliken
  spin on each Fe, and a CHEAP multireference probe every 3rd point: UNO n_u
  (Head-Gordon index sum n_i(2-n_i) over the UKS natural occupations — NO CASSCF,
  just the natural orbitals of the reference density; fast).

ANALYSIS (COGEF, identical math to cc_cogef_probe):
  E(d) -> F(d) = dE/dd (numeric) -> F_max in nN (1 eV/A = 1.602 nN;
  1 Ha/Bohr = 82.387 nN); force-modified barrier ΔE‡(F) = max[E - F·(d-d0)] -
  min[...] for F = 0/1/2/3 nN (the tilted-line construction). VERDICT compares
  F_max against the C–C anchor (F_max 5.58 nN, room kinetics 4.5–5 nN): does the
  Fe–S bond break EASIER or HARDER than C–C?

HONEST SCOPE (this is screen-gate #1, not a passport):
  - The cubane is NOT the pyrrhotite Fe7S8 lattice; a terminal SH thiolate is NOT
    lattice sulfur (no lattice Madelung field, no cooperative neighbor tension).
  - HS reference (default) is not the BS ground state; absolute energies are on
    the HS surface (force screen only). PBE on a stretching M–S bond is
    single-reference-biased; UNO n_u only flags where multireference character
    turns on — a CASSCF(2,2)-style follow-up would be the paspecheck.
  - Static COGEF force != a mill impulse; no entropy / temperature.
  Every number is computed at run time; timeouts, rigid fallbacks and non-
  convergence are recorded as such. Resume-aware (per-point checkpoints).

ENV KNOBS (all optional):
  FES_STAGE_TIMEOUT  seconds per scan point (relax cap), default 3000
  FES_DMAX           max stretch beyond equilibrium in A, default 1.5
  FES_STEP           stretch step in A, default 0.1
  FES_RIGID          1 = force the rigid scan for all points (default 0)
  FES_REF            reference spin state: 'hs' (default) or 'bs'
  PYSCF_MAX_MEMORY   per-process memory hint (MB)

OUTPUT (incremental / atomic — nothing lost on a kill; resume from the JSON):
  calc/fes_cogef_results.json

Run (heavy; meant for AWS via calc/fes_cogef_aws.py):
    python3 -u calc/fes_cogef.py
Smoke (fast local sanity, 2 points, no AWS):
    OMP_NUM_THREADS=1 python3 calc/fes_cogef.py --smoke
"""
import os
import sys
import json
import time
import signal
import argparse
import tempfile

try:
    import numpy as np
except ImportError:
    sys.exit("[deps] numpy is required:  python3 -m pip install numpy")
try:
    from pyscf import gto, scf, dft, lib
except ImportError as e:
    sys.exit(f"[deps] pyscf is required:  python3 -m pip install pyscf   ({e})")

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "fes_cogef_results.json")
SRC_XYZ = os.path.join(HERE, "femoco_4fe4s_cubane.xyz")

EV = 27.211386245988
NN_PER_HA_BOHR = 82.38723498        # 1 Ha/Bohr = 82.387 nN
A2B = 1.8897259886                  # Angstrom -> Bohr

BASIS = "def2-svp"
CHARGE = -2
SPIN_HS = 18                        # HS seed: 2 Fe3+ (d5) + 2 Fe2+ (d6), aligned
FE_ATOMS = (0, 1, 2, 3)
FE_STRETCH = 0                      # Fe0 — the metal end of the stretched bond
S_STRETCH = 8                       # terminal thiolate S bonded to Fe0
H_STRETCH = 12                      # H on that terminal S (S8-H12 = 1.34 A)
BS_FLIP = (2, 3)                    # BS pattern Fe0+ Fe1+ | Fe2- Fe3-  (fe4s4 recipe)

STAGE_TIMEOUT = int(os.environ.get("FES_STAGE_TIMEOUT", "3000"))
DMAX = float(os.environ.get("FES_DMAX", "1.5"))
STEP = float(os.environ.get("FES_STEP", "0.1"))
FORCE_RIGID = os.environ.get("FES_RIGID", "0") == "1"
REF = os.environ.get("FES_REF", "hs").lower()   # 'hs' | 'bs'
REF_SPIN = 0 if REF == "bs" else SPIN_HS

FGRID_NN = [0.0, 1.0, 2.0, 3.0]     # forces for the ΔE‡(F) barrier-under-force
CC_ANCHOR_NN = 5.58                 # cc_cogef C–C F_max
CC_ROOM_NN = (4.5, 5.0)             # cc_cogef room-temperature kinetics window


def say(m):
    print(f"[{time.strftime('%H:%M:%S')}] {m}", flush=True)


def atomic_save(obj):
    fd, tmp = tempfile.mkstemp(dir=HERE, suffix=".tmp")
    with os.fdopen(fd, "w") as f:
        json.dump(obj, f, ensure_ascii=False, indent=1)
    os.replace(tmp, OUT)


# ------------------------------------------------------------------ timeout
class StageTimeout(Exception):
    pass


def _alarm_handler(signum, frame):
    raise StageTimeout()


def with_timeout(seconds, fn, *args, **kwargs):
    """Best-effort per-point wall-clock cap via SIGALRM (main thread only).
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


# ------------------------------------------------------------------ geometry
def read_cubane_xyz(path):
    """Read [Fe4S4(SH)4]2- cubane; return [[element, np.array(xyz)], ...]."""
    with open(path) as fh:
        lines = fh.read().splitlines()
    n = int(lines[0].split()[0])
    atoms = []
    for i in range(2, 2 + n):
        p = lines[i].split()
        atoms.append([p[0], np.array([float(p[1]), float(p[2]), float(p[3])])])
    if [a[0] for a in atoms[:4]] != ["Fe", "Fe", "Fe", "Fe"]:
        raise ValueError(f"unexpected metal order: {[a[0] for a in atoms[:4]]}")
    if atoms[S_STRETCH][0] != "S" or atoms[H_STRETCH][0] != "H":
        raise ValueError("terminal S/H indices do not match the cubane xyz")
    return atoms


def bond_len(atoms, i, j):
    return float(np.linalg.norm(atoms[j][1] - atoms[i][1]))


def set_bond_rigid(atoms, i, j, k, d):
    """Rigidly place the terminal-bond distance d(i,j)=d by translating atom j
    (the terminal S) and its H (k) along the i->j unit vector. Everything else
    stays put — the honest rigid-scan fallback."""
    a = [[el, xyz.copy()] for el, xyz in atoms]
    v = a[j][1] - a[i][1]
    v = v / (np.linalg.norm(v) + 1e-12)
    d0 = np.linalg.norm(atoms[j][1] - atoms[i][1])
    shift = (d - d0) * v
    a[j][1] = a[j][1] + shift
    a[k][1] = a[k][1] + shift               # H rides with its S
    return a


# ------------------------------------------------------------------ mean field
def _build_mf(mol):
    mf = dft.UKS(mol).density_fit()
    mf.xc = "pbe"
    mf.conv_tol = 1e-7
    mf.max_cycle = 200
    mf.level_shift = 0.3
    mf.verbose = 0
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


def converge(atoms):
    """Converge the reference UKS state at this geometry. For REF='bs': HS seed
    -> flip Fe2,Fe3 -> Ms=0 BS. Returns a converged mf. Newton restart on stall."""
    a = [(el, tuple(xyz)) for el, xyz in atoms]
    if REF == "bs":
        mol_hs = gto.M(atom=a, basis=BASIS, charge=CHARGE, spin=SPIN_HS,
                       verbose=0)
        mf_hs = _build_mf(mol_hs)
        mf_hs.level_shift = 0.4
        mf_hs.kernel()
        dm0 = flip_dm(mf_hs.make_rdm1(), mol_hs, BS_FLIP)
        mol = gto.M(atom=a, basis=BASIS, charge=CHARGE, spin=0, verbose=0)
        mf = _build_mf(mol)
        mf.kernel(dm0=dm0)
    else:
        mol = gto.M(atom=a, basis=BASIS, charge=CHARGE, spin=SPIN_HS, verbose=0)
        mf = _build_mf(mol)
        mf.level_shift = 0.4
        mf.kernel()
    if not mf.converged:
        mf2 = mf.newton()
        mf2.max_cycle = 100
        mf2.kernel(mf.mo_coeff, mf.mo_occ)
        mf2.mol = mf.mol
        mf2.converged = bool(mf2.converged)
        return mf2
    return mf


def fe_spin_pops(mf):
    """Mulliken spin population on each Fe (BS/HS sanity)."""
    dma, dmb = mf.make_rdm1()
    m = (dma - dmb) @ mf.get_ovlp()
    aoslice = mf.mol.aoslice_by_atom()
    out = []
    for a in FE_ATOMS:
        sl = slice(int(aoslice[a][2]), int(aoslice[a][3]))
        out.append(round(float(np.trace(m[sl, sl])), 3))
    return out


def grad_audit(mf):
    """Max and RMS residual nuclear force (Ha/Bohr): confirms the constrained
    relax settled (rigid points carry a large residual — that is the point)."""
    try:
        g = mf.nuc_grad_method().kernel()
        return round(float(np.abs(g).max()), 5), round(float(np.sqrt((g ** 2).mean())), 5)
    except Exception as e:
        return None, f"{type(e).__name__}: {str(e)[:60]}"


def uno_nu(mf):
    """Head-Gordon unpaired index n_u = sum n_i(2-n_i) over the UKS natural
    occupations (natural orbitals of Da+Db). Cheap MR probe — no CASSCF."""
    s = mf.get_ovlp()
    w, v = np.linalg.eigh(s)
    w = np.clip(w, 1e-12, None)
    shalf = (v * np.sqrt(w)) @ v.T
    dma, dmb = mf.make_rdm1()
    n = np.linalg.eigvalsh(shalf @ (dma + dmb) @ shalf)
    n = np.clip(n, 0.0, 2.0)
    return round(float(np.sum(n * (2.0 - n))), 3)


# ------------------------------------------------------------------ scan point
def relaxed_point(atoms0, i, j, k, d, timeout_s):
    """Constrained relax freezing d(i,j)=d, else rigid fallback. Returns
    (mf, geo_atoms, rigid_flag, relax_note)."""
    rigid_atoms = set_bond_rigid(atoms0, i, j, k, d)
    if FORCE_RIGID:
        return converge(rigid_atoms), rigid_atoms, True, "forced rigid"
    try:
        from pyscf.geomopt import geometric_solver
    except ImportError:
        return converge(rigid_atoms), rigid_atoms, True, "geometric not installed"
    mf0 = converge(rigid_atoms)
    fd, cpath = tempfile.mkstemp(dir=HERE, suffix=".txt")
    with os.fdopen(fd, "w") as fh:
        fh.write(f"$freeze\ndistance {i + 1} {j + 1}\n")   # 1-based indices
    try:
        mol_eq = with_timeout(timeout_s, geometric_solver.optimize,
                              mf0, constraints=cpath, maxsteps=40)
        geo = [[mol_eq.atom_symbol(a),
                np.array(mol_eq.atom_coords()[a]) / A2B] for a in range(mol_eq.natm)]
        return converge(geo), geo, False, "relaxed"
    except StageTimeout:
        say(f"    relax d={d:.2f} timed out ({timeout_s}s) -> rigid fallback")
        return converge(rigid_atoms), rigid_atoms, True, f"relax timeout {timeout_s}s"
    except Exception as e:
        say(f"    relax d={d:.2f} failed ({type(e).__name__}) -> rigid fallback")
        return converge(rigid_atoms), rigid_atoms, True, f"relax error: {type(e).__name__}"
    finally:
        try:
            os.unlink(cpath)
        except OSError:
            pass


# ------------------------------------------------------------------ COGEF math
def cogef_analysis(res):
    pts = sorted((float(k), v["e"]) for k, v in res["scan"].items()
                 if v.get("converged"))
    if len(pts) < 4:
        return
    D = np.array([p[0] for p in pts])
    E = np.array([p[1] for p in pts])
    E = E - E.min()
    dEdd = np.gradient(E, D)                             # Ha/A
    fmax_nn = float(dEdd.max() * NN_PER_HA_BOHR / A2B)   # -> nN
    d0 = float(D[0])
    rows = []
    for F in FGRID_NN:
        f_au = F / NN_PER_HA_BOHR * A2B                  # nN -> Ha/A
        Etil = E - f_au * (D - d0)
        i_ts = int(np.argmax(Etil))
        i_min = int(np.argmin(Etil[:i_ts + 1])) if i_ts > 0 else 0
        bar = float((Etil[i_ts] - Etil[i_min]) * EV)
        rows.append({"F_nN": F, "barrier_eV": round(bar, 3),
                     "d_ts_A": round(float(D[i_ts]), 2)})
    easier = fmax_nn < CC_ANCHOR_NN
    below_room = fmax_nn < CC_ROOM_NN[1]
    verdict = (
        f"Fe–S(тиолат) F_max={fmax_nn:.2f} нН "
        + ("ЛЕГЧЕ" if easier else "ТЯЖЕЛЕЕ")
        + f" разрыва C–C (якорь 5.58 нН). "
        + ("Ниже" if below_room else "Выше")
        + " окна комнатной кинетики C–C (4.5–5 нН): "
        + ("мельница-реактор вскрывает Fe–S при силах, при которых C–C уже "
           "рвётся — вскрытие решётки механохимически ДОСТУПНО (скрин-гейт "
           "пройдён, нужен переход к решётке Fe7S8)."
           if below_room else
           "Fe–S жёстче доступного окна — механохимическое вскрытие под "
           "вопросом, нужен решёточный расчёт / другой мотив.")
    )
    res["cogef"] = {
        "reference_state": REF.upper() + (f" (2S={SPIN_HS})" if REF != "bs" else " (Ms=0 BS)"),
        "F_max_nN": round(fmax_nn, 2),
        "cc_anchor_F_max_nN": CC_ANCHOR_NN,
        "cc_room_kinetics_nN": list(CC_ROOM_NN),
        "barrier_vs_force": rows,
        "verdict": verdict,
        "note": "ΔE‡(F)=max[E−F·(d−d0)]−min[...] по скану до d0+DMAX — нижняя "
                "оценка барьера (полный разрыв Fe–S дальше по координате); F_max "
                "= max dE/dd на сетке (1 эВ/Å = 1.602 нН)",
    }


# ------------------------------------------------------------------ main
def build_grid(d0):
    nsteps = int(round(DMAX / STEP))
    return [round(d0 + k * STEP, 3) for k in range(nsteps + 1)]


def main_real():
    say(f"FeS COGEF (Fe0–S8 thiolate, ref={REF.upper()}) — threads={lib.num_threads()}")
    atoms0 = read_cubane_xyz(SRC_XYZ)
    d0 = bond_len(atoms0, FE_STRETCH, S_STRETCH)
    grid = build_grid(d0)
    res = json.load(open(OUT)) if os.path.exists(OUT) else {}
    res.setdefault("meta", {})
    res["meta"].update({
        "what": "COGEF: relaxed constrained stretch of ONE terminal Fe0–S8(H) "
                "thiolate bond of the [Fe4S4(SH)4]2- cubane, UKS-PBE/def2-SVP+DF; "
                "constrained-relax (geomeTRIC) with a rigid-scan fallback; "
                "force-modified barrier ΔE‡(F)",
        "cluster": "[Fe4S4(SH)4]2- idealized cubane (FeMoco diary Stage 17)",
        "geometry_source": "femoco_4fe4s_cubane.xyz",
        "charge": CHARGE, "basis": BASIS,
        "reference_state": REF.upper() + (f" (2S={SPIN_HS})" if REF != "bs" else " (Ms=0 BS)"),
        "reference_choice_honesty": (
            "HS (2S=18) reference by default: robust single determinant that "
            "converges at every stretch point; the terminal Fe–S sigma stiffness "
            "dominates the COGEF force, so F=dE/dd on the consistent HS surface is "
            "a valid first screen. BS Ms=0 is the physical ground state but re-"
            "collapses across the scan (use FES_REF=bs to run it). Choice is "
            "recorded per run and per point." if REF != "bs" else
            "BS Ms=0 reference (Fe0+ Fe1+ | Fe2- Fe3-), HS-seeded per point; the "
            "physical ground-state surface, but per-point spin pattern is not "
            "guaranteed to hold under the stretch — pattern flagged via Fe spins."),
        "stretched_bond": {"Fe": FE_STRETCH, "S": S_STRETCH, "H": H_STRETCH,
                           "d0_A": round(d0, 3)},
        "grid_A": grid, "step_A": STEP, "dmax_A": DMAX,
        "forces_nN": FGRID_NN,
        "purpose": "Track T3 ГЕЙТ: «хвостохранилища как рудник наоборот» — "
                   "механохимическое вскрытие сульфидной решётки. Скрин-вопрос "
                   "номер один: рвётся ли связь Fe–S ЛЕГЧЕ или ТЯЖЕЛЕЕ C–C под "
                   "силой (COGEF-перенос C–C-методики на Fe–S мотив).",
        "honesty": "СКРИН-ГЕЙТ №1, НЕ ПАСПОРТ: кубан [Fe4S4] ≠ решётка "
                   "пирротина Fe7S8 (нет мадэлунговского поля решётки, нет "
                   "кооперативного натяжения соседей); терминальный тиолат SH ≠ "
                   "решёточная сера. HS-реф (по умолч.) ≠ BS-основное состояние — "
                   "абсолютные энергии на HS-поверхности (только силовой скрин). "
                   "PBE на растяжении M–S однодетерминантно смещён; UNO n_u лишь "
                   "маркирует включение мультиреференсности (CASSCF(2,2)-поверх "
                   "— следующий шаг). Статическая сила ≠ импульс мельницы; без "
                   "энтропии/температуры.",
        "mr_probe": "UNO n_u (Head-Gordon sum n_i(2-n_i)) из UKS каждую 3-ю "
                    "точку — без CASSCF, быстрый маркер мультиреференсности",
        "flip_motivation": "перенос ✅ cc-activation COGEF (C–C, F_max 5.58 нН) "
                           "на Fe–Ni–S решётки: мельница как реактор вскрытия",
        "stage_timeout_s": STAGE_TIMEOUT, "force_rigid": FORCE_RIGID,
    })
    res.setdefault("scan", {})
    atomic_save(res)

    prev = atoms0
    for gi, d in enumerate(grid):
        key = f"{d:.2f}"
        if key in res["scan"] and res["scan"][key].get("converged"):
            say(f"  d={key} Å — уже есть (рестарт)")
            continue
        t0 = time.time()
        try:
            mf, geo, rigid, note = relaxed_point(
                prev, FE_STRETCH, S_STRETCH, H_STRETCH, d, STAGE_TIMEOUT)
            e = float(mf.e_tot)
            ss = float(mf.spin_square()[0])
            pops = fe_spin_pops(mf)
            gmax, grms = grad_audit(mf)
            rec = {"d_A": d, "e": e, "converged": bool(mf.converged),
                   "s2": round(ss, 3), "fe_mulliken_spin": pops,
                   "grad_max_ha_bohr": gmax, "grad_rms_ha_bohr": grms,
                   "rigid": bool(rigid), "relax_note": note,
                   "xyz": [[el, [round(float(c), 5) for c in xyz]]
                           for el, xyz in geo],   # for a CASSCF(2,2) follow-up
                   "t_s": round(time.time() - t0, 1)}
            if gi % 3 == 0:
                try:
                    rec["uno_n_u"] = uno_nu(mf)
                except Exception as e2:
                    rec["uno_n_u_error"] = f"{type(e2).__name__}: {str(e2)[:60]}"
            res["scan"][key] = rec
            prev = [[el, np.array(xyz)] for el, xyz in geo]
            say(f"  d={key} Å: E={e:.6f} conv={mf.converged} <S²>={ss:.2f} "
                f"rigid={rigid} spins={pops} nu={rec.get('uno_n_u')} "
                f"({rec['t_s']}s)")
        except Exception as exc:
            res["scan"][key] = {"d_A": d, "converged": False,
                                "error": f"{type(exc).__name__}: {str(exc)[:120]}"}
            say(f"  d={key} Å: FAILED {exc}")
        atomic_save(res)

    cogef_analysis(res)
    done = sum(1 for v in res["scan"].values() if v.get("converged"))
    n_rigid = sum(1 for v in res["scan"].values() if v.get("rigid"))
    res["status"] = "ok" if done == len(grid) else f"partial({done}/{len(grid)})"
    res["meta"]["rigid_points"] = n_rigid
    atomic_save(res)
    say(f"status={res['status']} rigid={n_rigid}/{len(grid)} -> {OUT}")


def main_smoke():
    global OUT, DMAX, STEP, FORCE_RIGID
    OUT = os.path.join(tempfile.mkdtemp(prefix="fes_cogef_smoke_"), "smoke.json")
    DMAX, STEP, FORCE_RIGID = 0.2, 0.2, True    # 2 rigid points, HS ref
    atoms0 = read_cubane_xyz(SRC_XYZ)
    d0 = bond_len(atoms0, FE_STRETCH, S_STRETCH)
    grid = build_grid(d0)
    say(f"SMOKE: rigid HS scan, d0={d0:.3f} grid={grid}")
    es = {}
    for d in grid:
        mf, geo, rigid, note = relaxed_point(
            atoms0, FE_STRETCH, S_STRETCH, H_STRETCH, d, 60)
        es[d] = float(mf.e_tot)
        say(f"  d={d:.2f}: E={es[d]:.6f} conv={mf.converged} rigid={rigid} "
            f"spins={fe_spin_pops(mf)}")
    ok = (len(es) == 2 and es[grid[1]] > es[grid[0]])   # stretch raises energy
    say(f"SMOKE {'PASS' if ok else 'FAIL'} "
        f"(ΔE={round((es[grid[1]]-es[grid[0]])*EV,3)} эВ на +{STEP} Å)")
    return 0 if ok else 1


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--smoke", action="store_true")
    a = ap.parse_args()
    sys.exit(main_smoke() if a.smoke else (main_real() or 0))
