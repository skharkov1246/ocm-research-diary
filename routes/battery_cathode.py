#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
battery-cathode — issue #12: systematic active-space reduction for the cathode
redox centre.  How many qubits does [NiO6] REALLY need?

RECOVERED PIPELINE.  The original scripts/results of the 2026-06-08 diary
entries were lost (never merged); this file re-implements the pipeline from
the diary description.  Numbers from a new run are reported as a
"пересчёт восстановленным пайплайном" — they may differ from the lost run.

Model
-----
QM cluster:  [NiO6]^9-  (formal Ni3+ d7, six O2-), regular octahedron with
             d(Ni-O) = 1.93 A (default).
Embedding:   Madelung field of point charges on a SIMPLIFIED ideal layered
             rock-salt lattice ("NMC-like", i.e. LiMO2 with an average TM
             charge +3):
               * simple-cubic site grid with spacing d = d(Ni-O);
               * sites with (i+j+k) odd  -> O  (charge -2);
               * sites with (i+j+k) even -> cation; (111) cation planes
                 alternate TM(+3) / Li(+1)  (plane index m=(i+j+k)/2:
                 m even -> TM, m odd -> Li); the centre (0,0,0) is the QM Ni.
               * cube |i|,|j|,|k| <= L; Evjen fractional weights (x1/2 per
                 face) on the boundary [Evjen, Phys. Rev. 39 (1932) 675];
               * the 7 QM sites are removed from the charge array;
               * exact charge neutrality of (cluster + charges) is enforced
                 by uniformly smearing the small residual over the outermost
                 shell of charges ("edge charge compensation", standard
                 embedded-cluster practice).
Literature anchors (NOT computed here, markers only):
  - LiNiO2 layered R-3m (alpha-NaFeO2 type), a=2.878 A, c=14.19 A
    [Ohzuku et al., J. Electrochem. Soc. 140 (1993) 1862];  our lattice is an
    idealized cubic rock-salt abstraction of it with d(Ni-O)=1.93 A.
  - Ni(III)O6 in LiNiO2 is Jahn-Teller distorted (~4x1.91 + 2x2.09 A);
    we use a REGULAR octahedron at 1.93 A — JT ignored (assumption, see
    "risks" in meta).

Method (matches the H2-worker convention that reproduced the NiO anchor:
CASCI *without* fix_spin_; only the Sz sector is fixed, S^2 may relax —
that relaxation is itself the spin-ambiguity signal)
------------------------------------------------------------------------
1. UKS-PBE/def2-SVP SCF (density fitting, level shift, newton() fallback)
   for 2 spin states of Ni3+ d7:  2S=1 (low-spin t2g6 eg1) and
   2S=3 (high-spin t2g5 eg2).  Both stored in the JSON.
2. UKS natural orbitals (spin-summed density) -> NOON ranking by
   fractionality f = n(2-n) (frontier-distance tie-break) -> NESTED active
   windows: 4o(8q) -> 6o(12q) -> 8o(16q) -> 10o(20q) -> 12o(24q)
   [-> 14o(28q) with --include-28q].
3. CASCI in each window on the frozen UKS-NO core; for each point:
   E, CAS-NOON, n_u = sum n(2-n)  (also sum min(n,2-n) and excess=n_u-2S),
   <S^2>, ndet, wall time.
4. Plateau analysis: consecutive dE / d(n_u); heuristic plateau thresholds
   are recorded in the JSON (they are a *convention*, not physics).

Honesty / assumptions / risks (also written into results meta):
  * cluster charge -9 is unphysical in vacuum; ONLY the embedded system is
    meaningful, and only trends/differences within one embedding.
  * bare formal point charges next to QM oxygens (nearest at d=1.93 A) can
    over-polarize / leak electrons (no capping ECPs on the first cation
    shell) — a known limitation of this cheap embedding.
  * formal charges (+3/+1/-2) of an idealized lattice; no NMC cation
    disorder, no JT distortion, no geometry relaxation.
  * CASCI on frozen UKS-NO core is not orbital-optimized (no CASSCF);
    energies are upper bounds within each window; the *series trend* is the
    deliverable.
  * SCF may land on different d-configurations per spin; convergence flags
    and <S^2> are stored, failures are reported, not hidden.

Run
---
  OMP_NUM_THREADS=4 python3 routes/battery_cathode.py            # full
  OMP_NUM_THREADS=4 python3 routes/battery_cathode.py --include-28q
  OMP_NUM_THREADS=1 python3 routes/battery_cathode.py --smoke    # <=3 min
Results -> routes/battery_cathode_results.json (atomic, incremental,
restart-safe: finished CAS points are skipped on rerun; SCF/NO cached in
routes/battery_cathode_scf_spin*.npz keyed by a parameter hash).
"""
import argparse
import hashlib
import json
import os
import sys
import time
from itertools import product
from math import comb

import numpy as np

# --- threads from environment BEFORE heavy pyscf use -------------------------
from pyscf import lib
_OMP = os.environ.get('OMP_NUM_THREADS')
if _OMP:
    lib.num_threads(int(_OMP))

from pyscf import gto, scf, dft, mcscf, qmmm  # noqa: E402

HARTREE2EV = 27.211386245988
ROUTES = os.path.dirname(os.path.abspath(__file__))
RESULTS = os.path.join(ROUTES, 'battery_cathode_results.json')
SCF_CHK = os.path.join(ROUTES, 'battery_cathode_scf_spin{spin}.npz')
SCRIPT_VERSION = 'recovered-pipeline v1 (2026-07-02, issue #12)'
MAX_MEMORY = float(os.environ.get('PYSCF_MAX_MEMORY', 8000))


# =============================================================================
# JSON helpers: atomic incremental save + restart
# =============================================================================
def _jdef(o):
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, (np.floating,)):
        return float(o)
    if isinstance(o, np.ndarray):
        return o.tolist()
    return str(o)


def atomic_save(R, path=RESULTS):
    """Write JSON atomically (tmp file + os.replace) so a kill mid-write can
    never corrupt previously saved points."""
    tmp = path + '.tmp'
    with open(tmp, 'w') as f:
        json.dump(R, f, ensure_ascii=False, indent=1, default=_jdef)
        f.flush()
        os.fsync(f.fileno())
    os.replace(tmp, path)


def load_results(path=RESULTS):
    if os.path.exists(path):
        try:
            with open(path) as f:
                return json.load(f)
        except Exception as e:
            print(f'!! could not parse existing {path}: {e!r} -> starting fresh',
                  flush=True)
    return {}


def params_hash(d):
    return hashlib.sha256(
        json.dumps(d, sort_keys=True, default=_jdef).encode()).hexdigest()[:16]


# =============================================================================
# cluster + Madelung embedding
# =============================================================================
QM_SITES = [(0, 0, 0), (1, 0, 0), (-1, 0, 0),
            (0, 1, 0), (0, -1, 0), (0, 0, 1), (0, 0, -1)]


def site_charge(i, j, k):
    """Charge/type of a lattice site of the idealized layered rock salt."""
    s = i + j + k
    if s % 2 != 0:
        return -2.0, 'O'
    m = s // 2  # (111) cation plane index
    return (3.0, 'TM') if m % 2 == 0 else (1.0, 'Li')


def build_cluster(d_nio):
    """[NiO6]^9- regular octahedron."""
    atoms = [('Ni', (0.0, 0.0, 0.0))]
    for i, j, k in QM_SITES[1:]:
        atoms.append(('O', (d_nio * i, d_nio * j, d_nio * k)))
    return atoms


def build_embedding(L, d_nio):
    """Point-charge array: Evjen-weighted cube, QM sites removed, exact
    neutrality of (cluster -9) + (charges) enforced on the outer shell.

    Returns (coords[N,3] A, charges[N], info dict)."""
    assert L >= 2, 'L >= 2 required so the QM cluster is interior'
    qm = set(QM_SITES)
    coords, charges, types, outer_idx = [], [], [], []
    counts = {'O': 0, 'TM': 0, 'Li': 0}
    for i, j, k in product(range(-L, L + 1), repeat=3):
        if (i, j, k) in qm:
            continue
        q, typ = site_charge(i, j, k)
        w = 1.0
        for c in (i, j, k):
            if abs(c) == L:
                w *= 0.5  # Evjen face/edge/corner weight
        if max(abs(i), abs(j), abs(k)) == L:
            outer_idx.append(len(charges))
        coords.append((d_nio * i, d_nio * j, d_nio * k))
        charges.append(q * w)
        types.append(typ)
        counts[typ] += 1
    coords = np.asarray(coords)
    charges = np.asarray(charges)
    raw_sum = float(charges.sum())
    # cluster carries -9; removed sites carried +3-12=-9; for total neutrality
    # the charge array must sum to exactly +9.
    residual = 9.0 - raw_sum
    charges[outer_idx] += residual / len(outer_idx)
    # distance of the closest charge to any QM atom (electron-leak indicator)
    qm_xyz = np.array([(d_nio * i, d_nio * j, d_nio * k) for i, j, k in QM_SITES])
    dmin = float(min(np.linalg.norm(coords - p, axis=1).min() for p in qm_xyz))
    info = dict(
        lattice='idealized layered rock salt (NMC-like LiMO2 abstraction)',
        lattice_note=('literature anchor (NOT computed): LiNiO2 R-3m a=2.878 A '
                      'c=14.19 A [Ohzuku 1993]; here simple-cubic grid with '
                      'spacing d(Ni-O), (111)-alternating TM(+3)/Li(+1) planes'),
        d_nio_A=d_nio, L_cells=L, n_charges=int(len(charges)),
        counts=counts, evjen_weights=True,
        raw_sum_before_compensation=round(raw_sum, 6),
        edge_compensation_residual=round(residual, 6),
        n_outer_shell_sites=len(outer_idx),
        sum_charges=round(float(charges.sum()), 6),
        total_system_charge=round(float(charges.sum()) - 9.0, 6),
        min_charge_QM_distance_A=round(dmin, 4),
        assumptions=['formal charges +3/+1/-2, no cation disorder',
                     'regular octahedron, Jahn-Teller ignored',
                     'no capping ECPs on first cation shell '
                     '(electron leakage risk)'])
    return coords, charges, types, info


def make_mol(atoms, spin, basis):
    return gto.M(atom=[(el, xyz) for el, xyz in atoms], unit='Angstrom',
                 basis=basis, charge=-9, spin=spin, verbose=0,
                 max_memory=MAX_MEMORY)


# =============================================================================
# SCF (UKS-PBE, DF, level shift, newton fallback) + natural orbitals
# =============================================================================
def run_uks(mol, mm_coords, mm_charges, xc='pbe', grids_level=None,
            max_cycle=150, level_shift=0.5, newton_fallback=True):
    """Robust embedded UKS.  Stage 1: DIIS + level shift; stage 2: DIIS w/o
    shift from stage-1 density; stage 3: second-order (newton) fallback."""
    mf = dft.UKS(mol).density_fit()
    mf.xc = xc
    if grids_level is not None:
        mf.grids.level = grids_level
    mf = qmmm.mm_charge(mf, mm_coords, mm_charges, unit='Angstrom')
    mf.level_shift = level_shift
    mf.max_cycle = max_cycle
    mf.conv_tol = 1e-8
    t0 = time.time()
    mf.kernel()
    stages = [('diis+shift%.2f' % level_shift, bool(mf.converged))]
    if not mf.converged:
        dm = mf.make_rdm1()
        mf.level_shift = 0.0
        mf.kernel(dm)
        stages.append(('diis-noshift', bool(mf.converged)))
    final = mf
    if not mf.converged and newton_fallback:
        mfn = mf.newton()
        mfn.max_cycle = 60
        mfn.kernel(mf.make_rdm1())
        stages.append(('newton', bool(mfn.converged)))
        final = mfn
    dt = time.time() - t0
    ss, mult = final.spin_square()
    # Mulliken spin density on Ni (atom 0)
    dma, dmb = final.make_rdm1()
    S = mol.intor_symmetric('int1e_ovlp')
    spin_pop = np.einsum('ij,ji->i', dma - dmb, S)
    p0, p1 = mol.aoslice_by_atom()[0][2:4]
    ni_spin = float(spin_pop[p0:p1].sum())
    diag = dict(e_tot=float(final.e_tot), converged=bool(final.converged),
                stages=stages, s2=round(float(ss), 4),
                multiplicity=round(float(mult), 4),
                ni_mulliken_spin=round(ni_spin, 4),
                seconds=round(dt, 1), xc=xc,
                grids_level=grids_level, level_shift=level_shift)
    return final, diag


def natural_orbitals(mol, mf):
    """Spin-summed UKS density -> NOs via S^1/2 P S^1/2 (Lowdin)."""
    dma, dmb = mf.make_rdm1()
    P = dma + dmb
    S = mol.intor_symmetric('int1e_ovlp')
    e, U = np.linalg.eigh(S)
    Xh = (U * np.sqrt(e)) @ U.T          # S^{1/2}
    Xmh = (U / np.sqrt(e)) @ U.T         # S^{-1/2}
    w, V = np.linalg.eigh(Xh @ P @ Xh)
    idx = np.argsort(w)[::-1]
    noon = np.clip(w[idx], 0.0, 2.0)
    C = Xmh @ V[:, idx]                  # columns = NOs, occ-descending
    return C, noon


def select_active(noon, ncas, nelec, twoS):
    """NOON ranking: fractionality f=n(2-n) descending, tie-break = distance
    from the Fermi position.  Nested by construction (top-k of one ranking).
    Returns (column order for CASCI mo, ncore, nelecas, na, nb, active idx)."""
    n = len(noon)
    f = noon * (2.0 - noon)
    fermi = float(np.sum(noon >= 1.0))
    dist = np.abs(np.arange(n) - (fermi - 0.5))
    ranked = sorted(range(n), key=lambda p: (-f[p], dist[p], p))
    act = sorted(ranked[:ncas])
    act_set = set(act)
    core = [p for p in range(n) if p not in act_set and noon[p] >= 1.0]
    virt = [p for p in range(n) if p not in act_set and noon[p] < 1.0]
    nelecas = int(nelec - 2 * len(core))
    na = (nelecas + twoS) // 2
    nb = nelecas - na
    if not (0 <= nb <= na <= ncas and na - nb == twoS):
        raise ValueError(f'inconsistent window: ncas={ncas} nelecas={nelecas} '
                         f'na={na} nb={nb} twoS={twoS}')
    return core + act + virt, len(core), nelecas, na, nb, act


def run_casci_window(mol, mm_coords, mm_charges, C_no, noon, ncas, twoS,
                     max_dets=2.0e7):
    """CASCI in the given NO window on the frozen NO core (no orbital opt).
    Only the Sz sector is fixed (H2-worker convention: no fix_spin_)."""
    order, ncore, nelecas, na, nb, act = select_active(
        noon, ncas, mol.nelectron, twoS)
    ndet = comb(ncas, na) * comb(ncas, nb)
    rec = dict(ncas=ncas, qubits=2 * ncas, nelecas=nelecas, na=na, nb=nb,
               ncore=ncore, ndet=int(ndet),
               window_noon_uks=[round(float(noon[p]), 4) for p in act])
    if ndet > max_dets:
        rec.update(skipped=f'ndet {ndet:.3g} > max_dets {max_dets:.3g}')
        return rec
    mref = qmmm.mm_charge(scf.ROHF(mol), mm_coords, mm_charges,
                          unit='Angstrom')
    mref.verbose = 0
    mc = mcscf.CASCI(mref, ncas, (na, nb))
    mc.verbose = 0
    mc.canonicalization = False
    mc.fcisolver.max_cycle = 400
    t0 = time.time()
    mc.kernel(C_no[:, order])
    dt = time.time() - t0
    dm1 = mc.fcisolver.make_rdm1(mc.ci, ncas, (na, nb))
    occ = np.sort(np.linalg.eigh(dm1)[0])[::-1]
    occ = np.clip(occ, 0.0, 2.0)
    ss, mult = mc.fcisolver.spin_square(mc.ci, ncas, (na, nb))
    n_u = float(np.sum(occ * (2.0 - occ)))       # Head-Gordon, per prompt
    n_u_min = float(np.sum(np.minimum(occ, 2.0 - occ)))
    rec.update(e_tot=float(mc.e_tot),
               converged=bool(getattr(mc, 'converged', True)),
               noon=[round(float(x), 4) for x in occ],
               n_u=round(n_u, 4),
               n_u_min=round(n_u_min, 4),
               excess=round(n_u - twoS, 4),      # n_u - 2S (formal sector)
               s2=round(float(ss), 4), multiplicity=round(float(mult), 4),
               seconds=round(dt, 1))
    return rec


# =============================================================================
# analysis: where does the window converge (plateau)?
# =============================================================================
DE_PLATEAU_EV = 0.03   # heuristic conventions, recorded in JSON
DNU_PLATEAU = 0.05


def analyze_series(series):
    rows = sorted([r for r in series if 'e_tot' in r], key=lambda r: r['ncas'])
    deltas, plateau = [], None
    for a, b in zip(rows, rows[1:]):
        dE = (b['e_tot'] - a['e_tot']) * HARTREE2EV
        dnu = b['n_u'] - a['n_u']
        deltas.append(dict(from_o=a['ncas'], to_o=b['ncas'],
                           from_q=2 * a['ncas'], to_q=2 * b['ncas'],
                           dE_eV=round(dE, 4), d_n_u=round(dnu, 4)))
        if plateau is None and abs(dE) < DE_PLATEAU_EV and abs(dnu) < DNU_PLATEAU:
            plateau = b['ncas']
    return dict(deltas=deltas,
                plateau_ncas=plateau,
                plateau_qubits=(2 * plateau if plateau else None),
                plateau_criteria=dict(dE_eV=DE_PLATEAU_EV, d_n_u=DNU_PLATEAU,
                                      note='heuristic convention, not physics'),
                note=('plateau = first window whose increment vs the previous '
                      'window is below both thresholds; None = no plateau '
                      'inside the computed series'))


# =============================================================================
# SCF cache (restart without redoing SCF)
# =============================================================================
def scf_cache_load(spin, phash):
    path = SCF_CHK.format(spin=spin)
    if not os.path.exists(path):
        return None
    try:
        z = np.load(path, allow_pickle=True)
        if str(z['phash']) != phash:
            return None
        return dict(C_no=z['C_no'], noon=z['noon'],
                    diag=json.loads(str(z['diag'])))
    except Exception as e:
        print(f'!! scf cache {path} unreadable: {e!r}', flush=True)
        return None


def scf_cache_save(spin, phash, C_no, noon, diag):
    np.savez_compressed(SCF_CHK.format(spin=spin), phash=phash, C_no=C_no,
                        noon=noon, diag=json.dumps(diag, default=_jdef))


# =============================================================================
# main
# =============================================================================
def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    ap.add_argument('--smoke', action='store_true',
                    help='fast validation: cluster + reduced charge grid, '
                         'one SCF + tiny CASCI; <=3 min on 1 core')
    ap.add_argument('--d-nio', type=float, default=1.93,
                    help='Ni-O distance, A (1.90-1.95 physical range)')
    ap.add_argument('--lattice-L', type=int, default=4,
                    help='charge cube half-size in cells (full run)')
    ap.add_argument('--basis', default='def2-svp')
    ap.add_argument('--xc', default='pbe')
    ap.add_argument('--spins', default='1,3',
                    help='comma list of 2S values (1=low-spin d7, 3=high-spin)')
    ap.add_argument('--windows', default='4,6,8,10,12',
                    help='comma list of active-space sizes (orbitals)')
    ap.add_argument('--include-28q', action='store_true',
                    help='append the 14o (28-qubit) point (heavy: ~1e7 dets)')
    ap.add_argument('--force', action='store_true',
                    help='recompute even if cached/present in results JSON')
    args = ap.parse_args()

    spins = [int(s) for s in args.spins.split(',') if s.strip()]
    windows = [int(w) for w in args.windows.split(',') if w.strip()]
    if args.include_28q and 14 not in windows:
        windows.append(14)
    windows = sorted(set(windows))

    t_start = time.time()
    print(f'== battery_cathode.py {SCRIPT_VERSION} ==', flush=True)
    print(f'   threads={lib.num_threads()}  max_memory={MAX_MEMORY} MB',
          flush=True)

    # ---------------- smoke mode ----------------
    if args.smoke:
        run_smoke(args)
        return

    # ---------------- full run ------------------
    params = dict(mode='full', d_nio=args.d_nio, L=args.lattice_L,
                  basis=args.basis, xc=args.xc, cluster_charge=-9,
                  version=SCRIPT_VERSION)
    ph = params_hash(params)

    R = load_results()
    if R and R.get('meta', {}).get('params_hash') not in (None, ph):
        bak = RESULTS + '.bak'
        os.replace(RESULTS, bak)
        print(f'!! params changed vs existing results -> old file moved to '
              f'{bak}, starting fresh', flush=True)
        R = {}
    R.setdefault('meta', {}).update(
        params=params, params_hash=ph, script=SCRIPT_VERSION,
        method=('UKS-PBE(DF)/def2-SVP embedded SCF -> UKS natural orbitals -> '
                'NOON-ranked nested CASCI windows on frozen NO core; '
                'Sz fixed, no fix_spin_ (H2-worker convention)'),
        honesty=('recovered pipeline; numbers may differ from the lost '
                 '2026-06-08 run; embedded-cluster absolute energies are NOT '
                 'voltages; only within-embedding trends are meaningful'),
        last_update=time.strftime('%Y-%m-%d %H:%M:%S'))

    atoms = build_cluster(args.d_nio)
    mm_coords, mm_charges, mm_types, emb = build_embedding(
        args.lattice_L, args.d_nio)
    R['cluster'] = dict(atoms=[[el, list(x)] for el, x in atoms],
                        charge=-9, formal='Ni3+ d7 + 6 O2-',
                        d_nio_A=args.d_nio)
    R['embedding'] = emb  # generation is deterministic from params
    atomic_save(R)
    print(f'cluster [NiO6]9- d={args.d_nio} A; charges: {emb["n_charges"]} '
          f'(sum {emb["sum_charges"]}, total system '
          f'{emb["total_system_charge"]})', flush=True)

    R.setdefault('scf', {})
    R.setdefault('cas_series', {})
    R.setdefault('analysis', {})

    for twoS in spins:
        skey = f'spin{twoS}'
        mol = make_mol(atoms, twoS, args.basis)
        sph = params_hash(dict(params, spin=twoS))

        cache = None if args.force else scf_cache_load(twoS, sph)
        if cache:
            C_no, noon, diag = cache['C_no'], cache['noon'], cache['diag']
            print(f'[{skey}] SCF loaded from cache: E={diag["e_tot"]:.6f} '
                  f'converged={diag["converged"]}', flush=True)
        else:
            print(f'[{skey}] UKS-{args.xc.upper()} SCF (2S={twoS}) ...',
                  flush=True)
            mf, diag = run_uks(mol, mm_coords, mm_charges, xc=args.xc)
            C_no, noon = natural_orbitals(mol, mf)
            scf_cache_save(twoS, sph, C_no, noon, diag)
            print(f'[{skey}] E={diag["e_tot"]:.6f} conv={diag["converged"]} '
                  f'<S2>={diag["s2"]} Ni_spin={diag["ni_mulliken_spin"]} '
                  f'({diag["seconds"]} s, stages {diag["stages"]})', flush=True)
        diag['noon_frontier_uks'] = [round(float(x), 4) for x in
                                     noon[np.argsort(noon * (2 - noon))[::-1][:16]]]
        R['scf'][skey] = diag
        atomic_save(R)

        series = R['cas_series'].setdefault(skey, {})
        for ncas in windows:
            wkey = f'{ncas}o'
            if wkey in series and 'e_tot' in series[wkey] and not args.force:
                print(f'[{skey}] {wkey}: already in results, skip', flush=True)
                continue
            print(f'[{skey}] CASCI window {ncas}o ({2*ncas}q) ...', flush=True)
            try:
                rec = run_casci_window(mol, mm_coords, mm_charges,
                                       C_no, noon, ncas, twoS)
            except Exception as e:
                rec = dict(ncas=ncas, qubits=2 * ncas, error=repr(e))
            series[wkey] = rec
            if 'e_tot' in rec:
                print(f'[{skey}] {wkey}: E={rec["e_tot"]:.6f} '
                      f'CAS({rec["nelecas"]}e,{ncas}o) ndet={rec["ndet"]} '
                      f'n_u={rec["n_u"]} <S2>={rec["s2"]} '
                      f'({rec["seconds"]} s)', flush=True)
            else:
                print(f'[{skey}] {wkey}: {rec.get("skipped") or rec.get("error")}',
                      flush=True)
            R['analysis'][skey] = analyze_series(series.values())
            R['meta']['last_update'] = time.strftime('%Y-%m-%d %H:%M:%S')
            atomic_save(R)   # incremental: every point survives a kill

    # cross-spin gap per window (caveat: different frozen NO cores per spin)
    if len(spins) >= 2:
        gaps = {}
        s0, s1 = f'spin{spins[0]}', f'spin{spins[1]}'
        for wkey, r0 in R['cas_series'].get(s0, {}).items():
            r1 = R['cas_series'].get(s1, {}).get(wkey)
            if r1 and 'e_tot' in r0 and 'e_tot' in r1:
                gaps[wkey] = round((r1['e_tot'] - r0['e_tot']) * HARTREE2EV, 4)
        R['analysis']['spin_gap_eV'] = dict(
            gaps=gaps, sign=f'{s1} minus {s0}',
            caveat='different frozen NO cores per spin -> gap includes core '
                   'relaxation error; qualitative only')
        atomic_save(R)

    print(f'== done in {time.time()-t_start:.1f} s -> {RESULTS} ==', flush=True)


# =============================================================================
# smoke: cluster + charges into JSON + one reduced-grid SCF + tiny CASCI
# =============================================================================
def run_smoke(args):
    t0 = time.time()
    L = 2            # reduced charge grid (125 sites)
    twoS = 1         # low-spin doublet only
    R = load_results()
    atoms = build_cluster(args.d_nio)
    mm_coords, mm_charges, mm_types, emb = build_embedding(L, args.d_nio)
    smoke = dict(started=time.strftime('%Y-%m-%d %H:%M:%S'),
                 note=('SMOKE: pipeline validation only (reduced charge grid '
                       'L=2, coarse DFT grid, 4o CASCI); NOT production '
                       'numbers'),
                 cluster=dict(atoms=[[el, list(x)] for el, x in atoms],
                              charge=-9, spin=twoS, basis=args.basis),
                 embedding=emb,
                 charges_full_list=[[round(float(c[0]), 4), round(float(c[1]), 4),
                                     round(float(c[2]), 4), round(float(q), 6), t]
                                    for c, q, t in zip(mm_coords, mm_charges,
                                                       mm_types)])
    R['smoke'] = smoke
    R.setdefault('meta', {})['script'] = SCRIPT_VERSION
    atomic_save(R)
    print(f'[smoke] cluster + {emb["n_charges"]} charges written '
          f'(sum {emb["sum_charges"]}, total system '
          f'{emb["total_system_charge"]})', flush=True)

    mol = make_mol(atoms, twoS, args.basis)
    print(f'[smoke] UKS-{args.xc.upper()} 2S={twoS} on reduced grid ...',
          flush=True)
    mf, diag = run_uks(mol, mm_coords, mm_charges, xc=args.xc,
                       grids_level=1, max_cycle=80)
    smoke['scf'] = diag
    atomic_save(R)
    print(f'[smoke] SCF E={diag["e_tot"]:.6f} conv={diag["converged"]} '
          f'<S2>={diag["s2"]} Ni_spin={diag["ni_mulliken_spin"]} '
          f'({diag["seconds"]} s, stages {diag["stages"]})', flush=True)

    C_no, noon = natural_orbitals(mol, mf)
    frontier = noon[np.argsort(noon * (2 - noon))[::-1][:8]]
    smoke['noon_frontier_uks'] = [round(float(x), 4) for x in frontier]
    try:
        rec = run_casci_window(mol, mm_coords, mm_charges, C_no, noon,
                               ncas=4, twoS=twoS)
    except Exception as e:
        rec = dict(ncas=4, error=repr(e))
    smoke['casci_4o'] = rec
    smoke['seconds_total'] = round(time.time() - t0, 1)
    atomic_save(R)
    if 'e_tot' in rec:
        print(f'[smoke] CASCI 4o(8q): E={rec["e_tot"]:.6f} '
              f'CAS({rec["nelecas"]}e,4o) ndet={rec["ndet"]} n_u={rec["n_u"]} '
              f'({rec["seconds"]} s)', flush=True)
    else:
        print(f'[smoke] CASCI 4o failed: {rec.get("error")}', flush=True)
    print(f'[smoke] done in {smoke["seconds_total"]} s -> {RESULTS}',
          flush=True)


if __name__ == '__main__':
    main()
