#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
battery-cathode — issue #12: VQE-PoC on the reduced NOON window (12-16 qubits).

Bridge "embedded [NiO6]^9- -> qubits -> VQE" for the DOUBLET Ni3+ (2S=1 —
the lower embedded SCF in battery_cathode_results.json):

  1. Cluster / Madelung embedding / SCF / natural-orbital machinery is
     IMPORTED from routes/battery_cathode.py (same params hash, so the
     battery_cathode_scf_spin1.npz cache is reused when present; otherwise
     the embedded UKS-PBE SCF is recomputed and cached for both scripts).
  2. CASCI reference in the NOON-ranked 6o window (12 qubits; --windows 6,8
     adds the 8o/16q point) on the frozen UKS-NO core.  The window is
     IDENTICAL BY CONSTRUCTION to battery_cathode_results.json
     (bc.select_active is imported, not re-implemented) and the CASCI energy
     is cross-checked against that file.
  3. Active-space integrals -> Jordan-Wigner qubit Hamiltonian (openfermion;
     the recipe of routes/h2_oer_quantum.py that put 18q UCCSD within
     0.26 kcal/mol of CASCI on the H2 track).  Correctness gate: lowest
     eigenvalue of H restricted to the CASCI (N, Sz) sector == E_CASCI.
  4. UCCSD-VQE (PennyLane, lightning.qubit statevector, adjoint gradients).
     OPEN SHELL (odd electron count): the Sz sector is enforced by the
     ROHF-type init_state (interleaved spin orbitals, even = alpha; alpha
     fills na lowest window orbitals, beta fills nb) plus Sz-conserving
     excitation lists — the same solution as h2_oer_quantum.py.  S^2 is NOT
     constrained, matching the pipeline convention "Sz fixed, no fix_spin_".

Self-test (always, and in --smoke): H2/def2-SVP CAS(2e,2o) — JW sector
minimum == CASCI to < 1e-6 Ha.

Honesty
-------
  * inherits ALL embedding caveats of battery_cathode.py (bare formal point
    charges, no capping ECPs, Jahn-Teller ignored, frozen NO core, no
    CASSCF); absolute energies are NOT voltages;
  * the deliverable is E_VQE - E_CASCI on the SAME active-space Hamiltonian
    (ansatz/optimizer quality on 12-16q), not new chemistry;
  * if UCCSD does not converge within the iteration/wall-time budget, the
    partial result is recorded as-is (status != 'converged') and the
    JW<->CASCI gate remains the deliverable — no silent retries;
  * --smoke does NOT touch the battery SCF (full-grid embedded UKS took
    ~600 s in battery_cathode_results.json — it cannot fit a 3-minute
    single-core budget); instead it validates the identical code path on an
    open-shell stand-in (OH radical, doublet) with a toy 4o window, no VQE.

Run
---
  OMP_NUM_THREADS=1 python3 routes/battery_cathode_vqe.py --smoke   # <=3 min
  OMP_NUM_THREADS=4 python3 routes/battery_cathode_vqe.py                # 6o/12q
  OMP_NUM_THREADS=4 python3 routes/battery_cathode_vqe.py --windows 6,8  # +16q
Results -> routes/battery_cathode_vqe_results.json (atomic, incremental,
restart-safe: windows with a finished VQE record are skipped on rerun).
"""
import argparse
import json
import os
import sys
import time

import numpy as np

# --- threads from environment BEFORE heavy pyscf use -------------------------
from pyscf import lib
_OMP = os.environ.get('OMP_NUM_THREADS')
if _OMP:
    lib.num_threads(int(_OMP))

from pyscf import gto, scf, mcscf, qmmm, ao2mo  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import battery_cathode as bc  # noqa: E402  reused: cluster, embedding, SCF, NOs

RESULTS = os.path.join(HERE, 'battery_cathode_vqe_results.json')
KCAL = 627.5094740631
SCRIPT_VERSION = 'vqe-poc v1 (2026-07-03, issue #12)'


# =============================================================================
# results JSON (atomic incremental, same convention as battery_cathode.py)
# =============================================================================
def save(R):
    bc.atomic_save(R, path=RESULTS)


def load_or_init():
    R = bc.load_results(RESULTS)
    R.setdefault('meta', {}).update(
        script=SCRIPT_VERSION,
        source_pipeline=bc.SCRIPT_VERSION,
        method=('doublet embedded UKS-NO 6o/8o NOON window (imported from '
                'battery_cathode.py) -> CASCI -> JW qubit Hamiltonian '
                '(openfermion, h2_oer_quantum.py recipe) -> sector-exact gate '
                '-> UCCSD-VQE (PennyLane lightning.qubit, adjoint)'),
        honesty=('PoC of ansatz/optimizer quality on the SAME active-space '
                 'Hamiltonian as CASCI; inherits all embedding caveats of '
                 'battery_cathode.py; absolute energies are NOT voltages; '
                 'S^2 not constrained (Sz-sector convention of the pipeline)'),
        last_update=time.strftime('%Y-%m-%d %H:%M:%S'))
    return R


# =============================================================================
# CASCI on a NOON window — same convention as bc.run_casci_window
# (bc.select_active is imported so the window is identical BY CONSTRUCTION;
#  re-implemented only because we need the mc object back for the integrals)
# =============================================================================
def casci_window(mol, C_no, noon, ncas, twoS, mm_coords=None, mm_charges=None):
    order, ncore, nelecas, na, nb, act = bc.select_active(
        noon, ncas, mol.nelectron, twoS)
    mref = scf.ROHF(mol)
    if mm_coords is not None:
        mref = qmmm.mm_charge(mref, mm_coords, mm_charges, unit='Angstrom')
    mref.verbose = 0
    mc = mcscf.CASCI(mref, ncas, (na, nb))
    mc.verbose = 0
    mc.canonicalization = False
    mc.fcisolver.max_cycle = 400
    t0 = time.time()
    mc.kernel(C_no[:, order])
    meta = dict(ncas=ncas, qubits=2 * ncas, nelecas=nelecas, na=na, nb=nb,
                ncore=ncore, e_casci=float(mc.e_tot),
                converged=bool(getattr(mc, 'converged', True)),
                window_noon_uks=[round(float(noon[p]), 4) for p in act],
                seconds=round(time.time() - t0, 1))
    return mc, meta


# =============================================================================
# JW qubit Hamiltonian + sector-exact gate (recipe of h2_oer_quantum.py)
# =============================================================================
def qubit_hamiltonian(mc):
    """Active-space integrals -> openfermion JW operator.
    Spin-orbital convention: even index = alpha (spinorb_from_spatial)."""
    import openfermion as of
    from openfermion.chem.molecular_data import spinorb_from_spatial
    ncas = mc.ncas
    h1e, ecore = mc.get_h1eff()
    h2e = ao2mo.restore(1, mc.get_h2eff(), ncas)
    # PySCF chemist (pq|rs) -> OpenFermion <pq|rs>: transpose(0,2,3,1)
    one_s, two_s = spinorb_from_spatial(
        np.asarray(h1e), np.asarray(h2e).transpose(0, 2, 3, 1))
    ham = of.InteractionOperator(float(ecore), one_s, 0.5 * two_s)
    return of.jordan_wigner(of.get_fermion_operator(ham))


def exact_in_sector(qubit_op, nq, na, nb):
    """Lowest eigenvalue of the JW Hamiltonian restricted to the CASCI
    (N, Sz) sector.  The full Fock minimum may live in a different sector —
    comparing to it would be a false gate failure (h2_oer_quantum lesson).
    Sparse build is >12 GB at >=18q (H2-track lesson) — callers must cap nq."""
    import openfermion as of
    from openfermion.linalg import jw_sz_restrict_operator
    from scipy.sparse.linalg import eigsh
    sp = of.get_sparse_operator(qubit_op, n_qubits=nq)
    sec = jw_sz_restrict_operator(sp, sz_value=0.5 * (na - nb),
                                  n_electrons=na + nb, n_qubits=nq)
    if sec.shape[0] <= 2:
        return float(np.linalg.eigvalsh(sec.toarray())[0])
    return float(eigsh(sec, k=1, which='SA', return_eigenvectors=False)[0])


def self_test_h2():
    """Gate for the JW recipe: H2/def2-SVP CAS(2e,2o), |JW - CASCI| < 1e-6 Ha."""
    mol = gto.M(atom='H 0 0 0; H 0 0 0.74', basis='def2-svp', verbose=0)
    mf = scf.RHF(mol)
    mf.verbose = 0
    mf.kernel()
    mc = mcscf.CASCI(mf, 2, 2)
    mc.verbose = 0
    mc.kernel()
    e_jw = exact_in_sector(qubit_hamiltonian(mc), 4, 1, 1)
    d = abs(e_jw - mc.e_tot)
    rec = dict(e_casci=float(mc.e_tot), e_jw_sector=e_jw,
               abs_diff_ha=float(d), passed=bool(d < 1e-6))
    if not rec['passed']:
        raise RuntimeError(f'H2 self-test FAILED: |JW-CASCI| = {d:.3e} Ha')
    return rec


# =============================================================================
# UCCSD-VQE, open-shell aware (h2_oer_quantum.py convention)
# =============================================================================
def sz_excitations(nq, na, nb):
    """ROHF-type reference determinant + Sz-conserving UCCSD excitations.
    Interleaved spin orbitals (even = alpha, matching spinorb_from_spatial):
    alpha occupies the na lowest window orbitals, beta the nb lowest — the
    ROHF determinant for occupation-descending NO windows (docc -> socc)."""
    init_state = np.zeros(nq, dtype=int)
    init_state[[2 * i for i in range(na)]] = 1
    init_state[[2 * i + 1 for i in range(nb)]] = 1
    occ = [i for i in range(nq) if init_state[i]]
    virt = [i for i in range(nq) if not init_state[i]]
    szf = lambda i: 0.5 if i % 2 == 0 else -0.5  # noqa: E731
    singles = [[o, v] for o in occ for v in virt if szf(o) == szf(v)]
    doubles = [[o1, o2, v1, v2]
               for ii, o1 in enumerate(occ) for o2 in occ[ii + 1:]
               for jj, v1 in enumerate(virt) for v2 in virt[jj + 1:]
               if abs(szf(o1) + szf(o2) - szf(v1) - szf(v2)) < 1e-9]
    return init_state, singles, doubles


def vqe_uccsd(qubit_op, nq, na, nb, iters=400, step=0.05, conv_tol=1e-7,
              max_seconds=None):
    """UCCSD-VQE on a statevector simulator.  Returns the record whatever
    happens: status in {'converged', 'max_iters', 'time_budget_exceeded'}."""
    import pennylane as qml
    H = qml.qchem.import_operator(qubit_op, format='openfermion')
    init_state, singles, doubles = sz_excitations(nq, na, nb)
    s_w, d_w = qml.qchem.excitations_to_wires(singles, doubles)
    try:
        dev = qml.device('lightning.qubit', wires=nq)
        devname = 'lightning.qubit'
    except Exception:
        dev = qml.device('default.qubit', wires=nq)
        devname = 'default.qubit'

    @qml.qnode(dev, diff_method='adjoint')
    def cost(w):
        qml.UCCSD(w, wires=range(nq), s_wires=s_w, d_wires=d_w,
                  init_state=init_state)
        return qml.expval(H)

    w = qml.numpy.zeros(len(singles) + len(doubles), requires_grad=True)
    opt = qml.AdamOptimizer(step)
    t0 = time.time()
    e_prev, hist, status = 1e9, [], 'max_iters'
    it = -1
    for it in range(iters):
        w, e = opt.step_and_cost(cost, w)
        if it % 25 == 0:
            hist.append(round(float(e), 6))
        if it > 40 and abs(e - e_prev) < conv_tol:
            status = 'converged'
            break
        if max_seconds is not None and time.time() - t0 > max_seconds:
            status = 'time_budget_exceeded'
            break
        e_prev = e
    e_vqe = float(cost(w))
    dt = time.time() - t0
    return dict(device=devname, n_params=len(singles) + len(doubles),
                n_singles=len(singles), n_doubles=len(doubles),
                e_vqe=e_vqe, iters=it + 1, status=status,
                adam_step=step, conv_tol=conv_tol,
                energy_trace_every25=hist,
                sec_per_iter=round(dt / max(it + 1, 1), 2),
                seconds=round(dt, 1))


# =============================================================================
# smoke: H2 self-test + toy OPEN-SHELL 4o window through the production path
# =============================================================================
def run_smoke(args):
    t0 = time.time()
    R = load_or_init()

    rec = self_test_h2()
    R['self_test_h2'] = rec
    save(R)
    print(f"[smoke] H2 self-test: CASCI={rec['e_casci']:.8f}  "
          f"JW={rec['e_jw_sector']:.8f}  |d|={rec['abs_diff_ha']:.2e} Ha  "
          f"passed={rec['passed']}", flush=True)

    # toy open-shell 4o window (8q): OH radical, doublet — odd electron like
    # Ni3+ 2S=1; exercises select_active / CASCI(na,nb) / JW / sector gate,
    # i.e. exactly the code path of the production 6o/12q point.  No VQE.
    mol = gto.M(atom='O 0 0 0; H 0 0 0.97', basis='def2-svp', spin=1,
                charge=0, verbose=0)
    mf = scf.UHF(mol)
    mf.verbose = 0
    mf.kernel()
    C_no, noon = bc.natural_orbitals(mol, mf)
    mc, meta = casci_window(mol, C_no, noon, ncas=4, twoS=1)
    qop = qubit_hamiltonian(mc)
    nq = 2 * meta['ncas']
    e_exact = exact_in_sector(qop, nq, meta['na'], meta['nb'])
    gate = abs(e_exact - meta['e_casci']) * KCAL
    toy = dict(system='OH radical (2S=1)/def2-SVP — open-shell stand-in; '
                      'NOT battery numbers',
               uhf_converged=bool(mf.converged), **meta,
               jw_terms=len(qop.terms), e_jw_sector=e_exact,
               gate_jw_vs_casci_kcal=round(gate, 6),
               passed=bool(gate < 1e-3))
    R['smoke'] = dict(
        note=('SMOKE: code-path validation only — no battery SCF (full-grid '
              'embedded UKS was ~600 s in battery_cathode_results.json, does '
              'not fit 3 min/1 core), no VQE optimization'),
        toy_4o=toy)
    save(R)
    print(f"[smoke] toy 4o(8q) OH: CAS({meta['nelecas']}e,4o) na={meta['na']} "
          f"nb={meta['nb']}  CASCI={meta['e_casci']:.8f}  "
          f"JW-sector={e_exact:.8f}  |d|={gate:.2e} kcal/mol  "
          f"passed={toy['passed']}", flush=True)
    if not toy['passed']:
        raise RuntimeError('smoke toy 4o gate FAILED')

    # 3-step UCCSD timing probe at 8q — wallclock anchor only, energies are
    # meaningless after 3 Adam steps and are deliberately not recorded.
    probe = vqe_uccsd(qop, nq, meta['na'], meta['nb'], iters=3)
    R['smoke']['uccsd_probe_8q'] = dict(
        device=probe['device'], n_params=probe['n_params'],
        sec_per_iter=probe['sec_per_iter'],
        note='timing anchor only (3 Adam steps), no converged energy')
    R['smoke']['seconds_total'] = round(time.time() - t0, 1)
    save(R)
    print(f"[smoke] UCCSD probe 8q: {probe['n_params']} params on "
          f"{probe['device']}, {probe['sec_per_iter']} s/iter", flush=True)
    print(f"[smoke] done in {R['smoke']['seconds_total']} s -> {RESULTS}",
          flush=True)


# =============================================================================
# full run: doublet SCF (cache or recompute) -> per-window CASCI/JW/VQE
# =============================================================================
def run_full(args):
    windows = sorted({int(w) for w in args.windows.split(',') if w.strip()})
    twoS = 1  # doublet — the lower embedded SCF (battery_cathode_results.json)

    R = load_or_init()
    R['self_test_h2'] = self_test_h2()
    save(R)
    print(f"[gate] H2 self-test passed "
          f"(|d|={R['self_test_h2']['abs_diff_ha']:.2e} Ha)", flush=True)

    # params MUST mirror battery_cathode.main() exactly to share its SCF cache
    params = dict(mode='full', d_nio=args.d_nio, L=args.lattice_L,
                  basis=args.basis, xc=args.xc, cluster_charge=-9,
                  version=bc.SCRIPT_VERSION)
    sph = bc.params_hash(dict(params, spin=twoS))
    atoms = bc.build_cluster(args.d_nio)
    mm_coords, mm_charges, _types, emb = bc.build_embedding(
        args.lattice_L, args.d_nio)
    mol = bc.make_mol(atoms, twoS, args.basis)
    print(f"[emb ] {emb['n_charges']} charges, total system charge "
          f"{emb['total_system_charge']}", flush=True)

    cache = None if args.force_scf else bc.scf_cache_load(twoS, sph)
    if cache:
        C_no, noon, diag = cache['C_no'], cache['noon'], cache['diag']
        src = 'npz cache (shared with battery_cathode.py)'
        print(f"[scf ] loaded from cache: E={diag['e_tot']:.6f} "
              f"converged={diag['converged']}", flush=True)
    else:
        print('[scf ] no cache -> embedded UKS-PBE SCF (2S=1), ~10 min on '
              '4 cores ...', flush=True)
        mf, diag = bc.run_uks(mol, mm_coords, mm_charges, xc=args.xc)
        C_no, noon = bc.natural_orbitals(mol, mf)
        bc.scf_cache_save(twoS, sph, C_no, noon, diag)
        src = 'recomputed here (cached for battery_cathode.py too)'
        print(f"[scf ] E={diag['e_tot']:.6f} converged={diag['converged']} "
              f"({diag['seconds']} s)", flush=True)
    R['scf'] = dict(diag, source=src, spin=twoS, params_hash=sph)
    save(R)
    if not diag['converged']:
        raise RuntimeError('embedded SCF unconverged — refusing to build '
                           'windows on a bad reference')

    # cross-reference: recovered-pipeline CASCI energies for the same windows
    ref = {}
    try:
        with open(bc.RESULTS) as f:
            ref = json.load(f)['cas_series'][f'spin{twoS}']
    except Exception:
        pass

    R.setdefault('windows', {})
    for ncas in windows:
        wkey = f'{ncas}o'
        W = R['windows'].setdefault(wkey, {})
        if 'vqe' in W and not args.force:
            print(f'[{wkey}] VQE record already present, skip (--force to redo)',
                  flush=True)
            continue

        mc, meta = casci_window(mol, C_no, noon, ncas, twoS,
                                mm_coords, mm_charges)
        if wkey in ref and 'e_tot' in ref.get(wkey, {}):
            meta['ref_e_casci_battery_cathode'] = ref[wkey]['e_tot']
            meta['diff_vs_ref_kcal'] = round(
                (meta['e_casci'] - ref[wkey]['e_tot']) * KCAL, 4)
        W['casci'] = meta
        save(R)
        nq, na, nb = 2 * ncas, meta['na'], meta['nb']
        print(f"[{wkey}] CASCI({meta['nelecas']}e,{ncas}o)={nq}q  "
              f"E={meta['e_casci']:.6f}"
              + (f"  (vs recovered pipeline: {meta['diff_vs_ref_kcal']:+.4f} "
                 f"kcal/mol)" if 'diff_vs_ref_kcal' in meta else ''),
              flush=True)

        qop = qubit_hamiltonian(mc)
        W['jw'] = dict(terms=len(qop.terms))
        if nq <= args.gate_maxq:
            t0 = time.time()
            e_exact = exact_in_sector(qop, nq, na, nb)
            gate = abs(e_exact - meta['e_casci']) * KCAL
            W['jw'].update(e_exact_sector=e_exact,
                           gate_vs_casci_kcal=round(gate, 6),
                           gate_passed=bool(gate < 1e-3),
                           gate_seconds=round(time.time() - t0, 1))
            print(f"[{wkey}] JW {W['jw']['terms']} terms; sector-exact gate "
                  f"|d|={gate:.2e} kcal/mol passed={W['jw']['gate_passed']}",
                  flush=True)
            if not W['jw']['gate_passed']:
                W['jw']['note'] = 'gate FAILED — VQE skipped, see honesty'
                save(R)
                continue
        else:
            W['jw']['gate'] = (f'skipped at {nq}q > --gate-maxq '
                               f'{args.gate_maxq} (sparse build memory, '
                               f'H2-track lesson); recipe validated by the '
                               f'H2 self-test and smaller windows')
            print(f"[{wkey}] JW {W['jw']['terms']} terms; exact gate skipped "
                  f"({nq}q)", flush=True)
        save(R)

        if args.no_vqe:
            continue
        print(f"[{wkey}] UCCSD-VQE on {nq}q (max {args.iters} iters, "
              f"budget {args.vqe_max_seconds} s) ...", flush=True)
        v = vqe_uccsd(qop, nq, na, nb, iters=args.iters, step=args.adam_step,
                      max_seconds=args.vqe_max_seconds)
        v['e_casci'] = meta['e_casci']
        v['vqe_minus_casci_kcal'] = round(
            (v['e_vqe'] - meta['e_casci']) * KCAL, 3)
        if v['status'] != 'converged':
            v['note'] = ('NOT converged within budget — recorded as-is; the '
                         'JW<->CASCI gate above is the validated deliverable')
        W['vqe'] = v
        R['meta']['last_update'] = time.strftime('%Y-%m-%d %H:%M:%S')
        save(R)
        print(f"[{wkey}] VQE: E={v['e_vqe']:.6f}  "
              f"VQE-CASCI={v['vqe_minus_casci_kcal']:+.3f} kcal/mol  "
              f"{v['n_params']} params, {v['iters']} iters "
              f"({v['status']}), {v['seconds']} s on {v['device']}",
              flush=True)

    print(f'== done -> {RESULTS} ==', flush=True)


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    ap.add_argument('--smoke', action='store_true',
                    help='H2 self-test + toy open-shell 4o window (no battery '
                         'SCF, no VQE); <=3 min on 1 core')
    ap.add_argument('--windows', default='6',
                    help='comma list of active-space sizes; 6 -> 12q, '
                         '"6,8" adds 16q')
    ap.add_argument('--d-nio', type=float, default=1.93)
    ap.add_argument('--lattice-L', type=int, default=4)
    ap.add_argument('--basis', default='def2-svp')
    ap.add_argument('--xc', default='pbe')
    ap.add_argument('--iters', type=int, default=400)
    ap.add_argument('--adam-step', type=float, default=0.05)
    ap.add_argument('--vqe-max-seconds', type=float, default=5400,
                    help='wall-time budget per VQE window; on expiry the '
                         'partial result is recorded honestly')
    ap.add_argument('--gate-maxq', type=int, default=16,
                    help='max qubits for the sector-exact JW gate')
    ap.add_argument('--no-vqe', action='store_true',
                    help='stop after CASCI + JW gate')
    ap.add_argument('--force', action='store_true',
                    help='recompute windows even if present in results JSON')
    ap.add_argument('--force-scf', action='store_true',
                    help='ignore the npz SCF cache')
    args = ap.parse_args()

    print(f'== battery_cathode_vqe.py {SCRIPT_VERSION} ==', flush=True)
    print(f'   threads={lib.num_threads()}  reusing pipeline: '
          f'{bc.SCRIPT_VERSION}', flush=True)
    if args.smoke:
        run_smoke(args)
    else:
        run_full(args)


if __name__ == '__main__':
    main()
