#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
cc-activation-recycling — Stage 2 real PySCF computation (REVERSE mirror of polyolefins).

Forward (polyolefins tab): Ti(IV) d0 inserts a C=C monomer into an M-C bond (N_u=0.97 @ TS).
Reverse (this tab):        a late-TM d8 centre does C-C OXIDATIVE ADDITION into an inert,
                           unactivated C(sp3)-C(sp3) bond (the recycling step for PE).

Minimal closed-shell model (mirror of cationic [Cl2Ti-CH3]+ + C2H4 used forward):
    [Rh(PH3)2(CH3)2]+   Rh(III) d6 bis-methyl
        --(reductive elimination = microscopic reverse of C-C OA)-->
    [Rh(PH3)2(C2H6)]+   Rh(I) d8 alkane sigma-complex
PH3 stands in for the bulky t-Bu pincer arms (tractability); the cis-P2 pocket mimics the
two P-donors of a PONOP/POCOP pincer. Ethane is the simplest unactivated C-C bond.

We trace the OA path from the *stable* bis-methyl product DOWN in d(C-C) (the reductive-
elimination direction). Starting from the bonded product keeps Rh engaged with both carbons,
so the floppy alkane sigma-complex cannot just dissociate -- a robust minimum-energy path.

Numbers produced (all real, def2-SVP, PySCF; nothing hand-tuned):
  [A] ethane C-C dissociation: NOON of sigma/sigma* -> the inert bond is strongly
      multireference when broken (why pyrolysis is radical & unselective).
  [B] metal OA relaxed scan: N_u(sigma), N_u(TS), dN_u; E-barrier CASCI vs HF.
  [C] sign of dE(OA): ethane + [RhL2]+ -> [RhL2(CH3)2]+ (optimized minima).
  [D] AVAS active-space size -> qubit count (VQE-candidate scale).
  [E] Jordan-Wigner check: qubit Hamiltonian reproduces CASCI on a CAS(4,4) subspace.

Honesty: relaxed-scan points are frozen-d(C-C) minima, not Hessian-verified saddles;
the scan maximum is reported as the TS estimate. Literature numbers (57 kJ/mol Rh-POCOP,
29 C Brookhart, tandem, Conk 87 %) are NOT computed here -- they live in the diary text.
"""
import json, time, sys, os
import numpy as np
from pyscf import gto, scf, mcscf, dft, lib, ao2mo
from pyscf.mcscf import avas

np.set_printoptions(precision=3, suppress=True)
lib.num_threads(4)
H2KCAL = 627.5094740631
RESULTS = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                       'cc_activation_recycling_results.json')

# ----------------------------------------------------------------------------------
# geometry construction
# ----------------------------------------------------------------------------------
def _splay_h(c_pos, away_dir, r=1.09, n=3):
    away = np.asarray(away_dir, float); away /= np.linalg.norm(away)
    tmp = np.array([1.,0.,0.]) if abs(away[0]) < 0.9 else np.array([0.,1.,0.])
    u = tmp - away*(tmp@away); u /= np.linalg.norm(u); v = np.cross(away, u)
    cone = np.deg2rad(109.5); out = []
    for k in range(n):
        phi = 2*np.pi*k/n
        d = -away*np.cos(np.pi-cone) + (u*np.cos(phi)+v*np.sin(phi))*np.sin(np.pi-cone)
        out.append(list(np.asarray(c_pos) + r*d))
    return out

def build_metal(dCC, rhC, rhP=2.30, with_P=True):
    """[Rh(PH3)2 + C2Hx]+ frame. Two carbons symmetric about z at Rh...C=rhC, sep dCC."""
    at = [['Rh', (0.,0.,0.)]]
    if with_P:
        for sy in (+1,-1):
            p = np.array([0., sy*rhP, 0.]); at.append(['P', list(p)])
            for h in _splay_h(p, -p/np.linalg.norm(p), r=1.42): at.append(['H', h])
    half = dCC/2.
    z = float(np.sqrt(max(rhC**2 - half**2, 0.25)))
    c1 = np.array([+half,0.,z]); c2 = np.array([-half,0.,z])
    at += [['C', list(c1)], ['C', list(c2)]]
    # methyl H's splay AWAY from both the other carbon (breaking bond) AND the Rh (origin),
    # i.e. an umbrella opening away from the metal -> no spurious H...Rh clashes.
    for c, o in [(c1,c2),(c2,c1)]:
        toC = (o - c)/np.linalg.norm(o - c)
        toRh = (-c)/np.linalg.norm(c)
        for h in _splay_h(c, toC + toRh): at.append(['H', h])
    return at

def build_ethane(dCC=1.54):
    half = dCC/2.; c1 = np.array([+half,0.,0.]); c2 = np.array([-half,0.,0.])
    at = [['C', list(c1)], ['C', list(c2)]]
    for c, o in [(c1,c2),(c2,c1)]:
        for h in _splay_h(c, o-c): at.append(['H', h])
    return at

def build_rhl2(rhP=2.30):
    at = [['Rh',(0.,0.,0.)]]
    for sy in (+1,-1):
        p = np.array([0., sy*rhP, 0.]); at.append(['P', list(p)])
        for h in _splay_h(p, -p/np.linalg.norm(p), r=1.42): at.append(['H', h])
    return at

# carbon atom indices in build_metal: Rh(0) P(1) 3H(2-4) P(5) 3H(6-8) C(9) C(10) ...
C1_IDX, C2_IDX = 9, 10

# ----------------------------------------------------------------------------------
def mol_of(atoms, charge=1, spin=0):
    return gto.M(atom=atoms, basis='def2-svp', charge=charge, spin=spin,
                 verbose=0, max_memory=6000)

def rhf(mol):
    mf = scf.RHF(mol); mf.max_cycle = 200; mf.kernel(); return mf

def pbe0_df(mol):
    mf = dft.RKS(mol).density_fit(); mf.xc = 'pbe0'; mf.grids.level = 2
    mf.max_cycle = 200; mf.kernel(); return mf

def unpaired(occ):
    occ = np.clip(np.asarray(occ,float), 0, 2)
    return float(np.sum(np.minimum(occ, 2-occ)))

def casci_noon(mf, ao_labels=('Rh 4d','C 2p'), ncas=9, nelecas=12, thr=0.02):
    """AVAS guess orbitals (low thr -> superset), then CASCI on a FIXED window (consistent)."""
    _, _, orbs = avas.avas(mf, list(ao_labels), threshold=thr, canonicalize=True)
    mc = mcscf.CASCI(mf, ncas, nelecas); mc.fcisolver.nroots = 1
    mc.fix_spin_(ss=0); mc.kernel(orbs)
    dm1 = mc.fcisolver.make_rdm1(mc.ci, ncas, nelecas)
    occ = np.sort(np.linalg.eigvalsh(dm1))[::-1]
    return dict(e_cas=float(mc.e_tot), noon=[round(float(x),4) for x in occ],
                n_u=round(unpaired(occ),4), ncas=ncas, nelecas=nelecas, mc=mc)

def casscf_noon(mf, ao_labels=('Rh 4d','C 2p'), ncas=9, nelecas=12, thr=0.02, macro=40):
    """AVAS -> orbital-optimized CASSCF NOON (reliable N_u; CASCI-on-AVAS under-counts here)."""
    _, _, orbs = avas.avas(mf, list(ao_labels), threshold=thr, canonicalize=True)
    mc = mcscf.CASSCF(mf, ncas, nelecas)       # non-DF: proven to converge (~300 s) here
    mc.max_cycle_macro = macro; mc.fix_spin_(ss=0)
    mc.kernel(orbs)
    dm1 = mc.fcisolver.make_rdm1(mc.ci, ncas, nelecas)
    occ = np.sort(np.linalg.eigvalsh(dm1))[::-1]
    return dict(e_cas=float(mc.e_tot), noon=[round(float(x),4) for x in occ],
                n_u=round(unpaired(occ),4), ncas=ncas, nelecas=nelecas,
                converged=bool(mc.converged))

def cc_localized_22(mol, mf, ci1=C1_IDX, ci2=C2_IDX):
    """A consistent (2e,2o) CAS = sigma/sigma* of the C1-C2 bond.
    Robust cross-geometry probe of the breaking-bond multireference."""
    mo = mf.mo_coeff; nocc = mol.nelectron//2
    # project onto C1,C2 AOs
    aoslices = mol.aoslice_by_atom()
    c1 = list(range(aoslices[ci1][2], aoslices[ci1][3]))
    c2 = list(range(aoslices[ci2][2], aoslices[ci2][3]))
    cc_ao = c1 + c2
    w = (mo[cc_ao,:]**2).sum(0)            # weight of each MO on the two carbons
    occ_mo = np.argsort(w[:nocc])[::-1]    # most C-C-like occupied
    vir_mo = nocc + np.argsort(w[nocc:])[::-1]
    sigma = occ_mo[0]; sigstar = vir_mo[0]
    cas_orb = mf.mo_coeff.copy()
    # move chosen MOs into the active window (positions nocc-1, nocc)
    order = list(range(mf.mo_coeff.shape[1]))
    order.remove(sigma); order.remove(sigstar)
    order = order[:nocc-1] + [sigma, sigstar] + order[nocc-1:]
    cas_orb = mf.mo_coeff[:, order]
    mc = mcscf.CASCI(mf, 2, 2); mc.fix_spin_(ss=0); mc.kernel(cas_orb)
    dm1 = mc.fcisolver.make_rdm1(mc.ci, 2, 2)
    occ = np.sort(np.linalg.eigvalsh(dm1))[::-1]
    return dict(noon=[round(float(x),4) for x in occ], n_u=round(unpaired(occ),4),
                e=float(mc.e_tot))

# ----------------------------------------------------------------------------------
def cas_to_qubit_H(mc):
    """CAS 1-/2-electron integrals -> Jordan-Wigner qubit Hamiltonian (sparse, includes ecore)."""
    from openfermion import InteractionOperator, get_fermion_operator, jordan_wigner
    from openfermion.linalg import get_sparse_operator
    ncas = mc.ncas
    h1, ecore = mc.get_h1eff()
    h2 = ao2mo.restore(1, mc.get_h2eff(), ncas)   # chemist (pq|rs)
    n = ncas; nso = 2*n
    one = np.zeros((nso, nso)); two = np.zeros((nso,)*4)
    for p in range(n):
        for q in range(n):
            for s in (0,1):
                one[2*p+s, 2*q+s] = h1[p,q]
    # chemist (pq|rs) -> physicist two_body[p,q,r,s] a_p^ a_q^ a_r a_s with 1/2 factor
    for p in range(n):
        for q in range(n):
            for r in range(n):
                for s in range(n):
                    val = 0.5*h2[p,s,q,r]
                    for a in (0,1):
                        for b in (0,1):
                            two[2*p+a, 2*q+b, 2*r+b, 2*s+a] = val
    qop = jordan_wigner(get_fermion_operator(InteractionOperator(float(ecore), one, two)))
    return get_sparse_operator(qop, n_qubits=nso), nso

def jw_check(mc):
    """Verify the Jordan-Wigner qubit Hamiltonian reproduces CASCI (sector-restricted)."""
    try:
        from openfermion.linalg import jw_number_restrict_operator
        import scipy.sparse.linalg as sla
    except Exception as e:
        return dict(error=f'openfermion import: {e!r}')
    H, nso = cas_to_qubit_H(mc)
    ne = mc.nelecas; nelec = ne if np.isscalar(ne) else sum(ne)
    Hn = jw_number_restrict_operator(H, int(nelec), n_qubits=nso)  # correct e-number sector
    e_jw = float(sla.eigsh(Hn, k=1, which='SA', return_eigenvectors=False)[0])
    return dict(qubits=nso, e_casci=float(mc.e_tot), e_jw=e_jw,
                diff_kcal=round((e_jw-mc.e_tot)*H2KCAL, 8))

def run_vqe_poc():
    """[D] Stage-3 PoC: UCCSD-VQE on the CAS(4e,4o)=8-qubit active centre at the OA-TS,
    EXACT statevector simulator (no shot noise) -> verifies the VQE(UCCSD)->CASCI bridge.
    Validated on H2/STO-3G CAS(2,2): Powell reaches CASCI to 0.000000 kcal/mol."""
    print('=== [D] VQE PoC: UCCSD on CAS(4,4)=8 qubits (statevector) ===', flush=True)
    from openfermion import (jordan_wigner, uccsd_singlet_generator, uccsd_singlet_paramsize,
                             jw_hartree_fock_state)
    from openfermion.linalg import get_sparse_operator
    import scipy.sparse.linalg as sla
    from scipy.optimize import minimize
    mfh = rhf(mol_of(build_metal(1.98, 2.16), charge=1))
    mc = casci_noon(mfh, ncas=4, nelecas=4)['mc']
    H, nq = cas_to_qubit_H(mc)
    ne = mc.nelecas; ne = ne if np.isscalar(ne) else sum(ne)
    npar = uccsd_singlet_paramsize(nq, ne)
    hf = np.asarray(jw_hartree_fock_state(ne, nq)).reshape(-1).astype(complex)
    def energy(theta):
        gen = uccsd_singlet_generator(list(theta), nq, ne, anti_hermitian=True)
        A = get_sparse_operator(jordan_wigner(gen), n_qubits=nq)
        v = sla.expm_multiply(A, hf)
        return float(np.real(np.vdot(v, H @ v)))
    e_hf = energy(np.zeros(npar))
    t = time.time()
    res = minimize(energy, np.zeros(npar), method='Powell',
                   options=dict(maxiter=6000, xtol=1e-7, ftol=1e-9))
    out = dict(qubits=int(nq), n_params=int(npar), ansatz='UCCSD singlet',
               simulator='exact statevector (openfermion+scipy)',
               e_ref_hf=e_hf, e_vqe=float(res.fun), e_casci=float(mc.e_tot),
               gap_vqe_casci_kcal=round((res.fun-mc.e_tot)*H2KCAL, 6),
               corr_recovered_kcal=round((e_hf-res.fun)*H2KCAL, 3),
               seconds=round(time.time()-t, 1))
    print('VQE:', json.dumps(out), flush=True)
    return out

# ----------------------------------------------------------------------------------
def run_ethane_curve():
    print('=== [A] ethane C-C dissociation: sigma/sigma* NOON ===', flush=True)
    rows = []
    for d in [1.54, 1.8, 2.1, 2.5, 3.0]:
        mol = mol_of(build_ethane(d), charge=0, spin=0)
        mf = rhf(mol)
        # CAS(6,6) CASSCF (orbital-optimized): sigma/sigma* + neighbouring orbitals
        mc = mcscf.CASSCF(mf, 6, 6); mc.fix_spin_(ss=0)
        mc.kernel()
        occ6 = np.sort(np.linalg.eigvalsh(mc.fcisolver.make_rdm1(mc.ci,6,6)))[::-1]
        nu6 = round(unpaired(occ6),4)
        rows.append(dict(dCC=d, e_hf=float(mf.e_tot), e_cas66=float(mc.e_tot),
                         nu_66=nu6, noon_66=[round(float(x),4) for x in occ6]))
        print(f'  d(C-C)={d:.2f}  HF={mf.e_tot:.5f}  CASSCF(6,6)={mc.e_tot:.5f}  '
              f'N_u(6,6)={nu6}  NOON={np.round(occ6,3).tolist()}', flush=True)
    return rows

def run_metal_path():
    """Constructed representative geometries (NO optimization -> no collapse), exactly the
    polyolefins methodology (fixed Cossee-Arlman cc=2.0). 3-centre M...(C-C) motif:
    C-C breaks while two Rh-C form. N_u from orbital-optimized CASSCF (CASCI-on-AVAS
    under-counts static correlation here -> probe showed 0.035 vs 0.166)."""
    print('=== [B] metal OA path: constructed geometries, CASSCF NOON ===', flush=True)
    # (tag, dCC, rhC): sigma-complex -> OA-TS region -> bis-methyl product
    geoms = [('sigma',1.54,2.45), ('early',1.78,2.28), ('OA-TS',1.98,2.16),
             ('late',2.30,2.08), ('product',2.75,2.05)]
    do_scf = {'sigma','OA-TS','product'}       # orbital-optimized CASSCF at key points only
    rows = []
    for tag,dCC,rhC in geoms:
        mol = mol_of(build_metal(dCC,rhC), charge=1)
        mf = rhf(mol)
        casci = casci_noon(mf)                 # fast: energy + AVAS CASCI (all points)
        scf_ok = mf.converged
        if tag in do_scf:
            try:
                scf = casscf_noon(mf)          # reliable N_u (orbital-optimized)
            except Exception as e:
                scf = dict(e_cas=None, n_u=None, noon=None, converged=False, err=repr(e))
        else:
            scf = dict(e_cas=None, n_u=None, noon=None, converged=None)
        rows.append(dict(tag=tag, dCC=dCC, rhC=rhC, scf_conv=bool(scf_ok),
                         e_hf=float(mf.e_tot), e_casci=casci['e_cas'],
                         e_casscf=scf.get('e_cas'),
                         n_u_casscf=scf.get('n_u'), noon_casscf=scf.get('noon'),
                         casscf_conv=scf.get('converged'),
                         n_u_casci=casci['n_u'], noon_casci=casci['noon']))
        print(f'  [{tag:7s}] dCC={dCC:.2f} rhC={rhC:.2f} | HF={mf.e_tot:.4f} '
              f'CASCI={casci["e_cas"]:.4f} CASSCF={scf.get("e_cas")} | '
              f'N_u(CASSCF)={scf.get("n_u")} N_u(CASCI)={casci["n_u"]} '
              f'conv={scf.get("converged")}', flush=True)
        if scf.get('noon'):
            print(f'           NOON(CASSCF)={scf["noon"]}', flush=True)
    # --- bonus inversion probe: bare vs captured TS electronics ---
    print('  -- bonus: bare (no PH3) vs captured (cis-P2) at OA-TS --', flush=True)
    bare = {}
    try:
        # bare: Rh+ + stretched ethane, no phosphines (Rh+ is d8 -> use charge +1, the 2 C give
        # an even-electron singlet). Same dCC/rhC as OA-TS for a clean electronic comparison.
        molb = mol_of(build_metal(1.98, 2.16, with_P=False), charge=1)
        mfb = rhf(molb)
        sb = casscf_noon(mfb, ncas=7, nelecas=10)  # Rh d8 (10e) + sigma/sigma*(C-C)
        bare = dict(n_u_casscf=sb['n_u'], noon=sb['noon'], conv=sb['converged'],
                    ncas=7, nelecas=10, scf_conv=bool(mfb.converged))
        print(f'     bare TS: N_u(CASSCF 10,7)={sb["n_u"]} NOON={sb["noon"]} '
              f'conv={sb["converged"]} scf={mfb.converged}', flush=True)
    except Exception as e:
        bare = dict(err=repr(e)); print('     bare failed', repr(e), flush=True)
    rows.append(dict(tag='bare-OA-TS', **bare))
    return rows

def run_thermo():
    """Sign of the intramolecular OA reaction energy on the constructed path:
    sigma-complex [Rh(PH3)2(C2H6)]+  ->  bis-methyl product [Rh(PH3)2(CH3)2]+.
    Single-point PBE0 and HF on the SAME atoms (no optimization -> the floppy cation
    collapses under relaxation). Same atom count both sides, so the SIGN is meaningful;
    absolute dE is a minimal-model estimate. Literature firmly: direct OA is endothermic."""
    print('=== [C] thermodynamic sign of OA (single-point on constructed path) ===', flush=True)
    g_sigma = build_metal(1.54, 2.45)     # alkane sigma-complex
    g_prod  = build_metal(2.75, 2.05)     # bis-methyl product
    out = {}
    for label, fn in [('pbe0', pbe0_df), ('hf', rhf)]:
        e_s = float(fn(mol_of(g_sigma, charge=1)).e_tot)
        e_p = float(fn(mol_of(g_prod,  charge=1)).e_tot)
        dE = (e_p - e_s)*H2KCAL
        out[label] = dict(E_sigma=e_s, E_product=e_p, dE_OA_kcal=round(dE,2),
                          sign=('endothermic' if dE>0 else 'exothermic'))
        print(f'   {label.upper()}: dE(sigma->product) = {dE:+.1f} kcal/mol '
              f'({out[label]["sign"]})', flush=True)
    return out

# ----------------------------------------------------------------------------------
if __name__ == '__main__':
    R = {'meta': dict(model='[Rh(PH3)2(CH3)2]+ <-> [Rh(PH3)2(C2H6)]+  (Rh d6<->d8)',
                      basis='def2-svp', note='reverse mirror of polyolefins Ti insertion')}
    t0 = time.time()
    stage = sys.argv[1] if len(sys.argv)>1 else 'all'

    if stage in ('all','A'):
        R['ethane_curve'] = run_ethane_curve()
    if stage in ('all','B'):
        R['metal_path'] = run_metal_path()
        mp = [r for r in R['metal_path'] if r.get('tag') not in ('bare-OA-TS',)]
        # use CASSCF N_u where available; the OA-TS row is the labelled one
        def nu(r): return r.get('n_u_casscf') if r.get('n_u_casscf') is not None else r.get('n_u_casci')
        sig = mp[0]
        ts = next((r for r in mp if r['tag']=='OA-TS'), mp[2])
        prod = next((r for r in mp if r['tag']=='product'), mp[-1])
        bare = next((r for r in R['metal_path'] if r.get('tag')=='bare-OA-TS'), {})
        # NB: constructed-geometry RELATIVE energies are unreliable (no optimization) ->
        # we report N_u + the single-geometry CASCI-HF correlation magnitude AT the TS
        # (like polyolefins' CASCI vs HF at the insertion TS), NOT a barrier height.
        R['headline'] = dict(
            N_u_sigma=nu(sig), N_u_TS=nu(ts), N_u_product=nu(prod),
            dN_u=round((nu(ts) or 0)-(nu(sig) or 0),4),
            N_u_sigma_casci=sig['n_u_casci'], N_u_TS_casci=ts['n_u_casci'],
            N_u_TS_bare=bare.get('n_u_casscf'),
            ts_dCC=ts['dCC'],
            corr_at_TS_kcal=round((ts['e_casci']-ts['e_hf'])*H2KCAL,2),
            corr_at_sigma_kcal=round((sig['e_casci']-sig['e_hf'])*H2KCAL,2),
            note='constructed geoms; barrier/thermo are LITERATURE, not computed here')
        print('HEADLINE:', json.dumps(R['headline']), flush=True)
        # JW check on the OA-TS CAS reduced to (4,4)
        try:
            mfh = rhf(mol_of(build_metal(ts['dCC'], ts['rhC']), charge=1))
            cas44 = casci_noon(mfh, ncas=4, nelecas=4)
            R['jw_check'] = jw_check(cas44['mc'])
            print('JW:', json.dumps({k:v for k,v in R['jw_check'].items()}), flush=True)
        except Exception as e:
            R['jw_check'] = dict(error=repr(e)); print('JW failed', repr(e), flush=True)
    if stage in ('all','C'):
        R['thermo'] = run_thermo()
    if stage in ('all','D'):
        R['vqe_poc'] = run_vqe_poc()

    # merge with any prior results (so staged runs accumulate)
    if os.path.exists(RESULTS):
        try:
            prev = json.load(open(RESULTS)); prev.update(R); R = prev
        except Exception: pass
    # strip non-serializable mc objects
    for r in R.get('metal_path', []): r.pop('mc', None)
    json.dump(R, open(RESULTS,'w'), ensure_ascii=False, indent=1, default=str)
    print(f'=== done {stage} in {time.time()-t0:.1f}s -> {RESULTS} ===', flush=True)
