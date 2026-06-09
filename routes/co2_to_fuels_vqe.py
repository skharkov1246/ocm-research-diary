#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
routes/co2_to_fuels_vqe.py — Stage 3 VQE-PoC вкладки co2-to-fuels.

Мост «активный центр → кубиты → VQE» на трактуемом подпространстве C–C сочленения:
  1. CASCI(4e,4o) на фронтире *CO–*CO (8 кубитов).
  2. Кубитный гамильтониан (Jordan–Wigner) строится из активных интегралов и
     ТОЧНО воспроизводит CASCI (сверка в секторе нужного числа электронов).
  3. VQE(UCCSD) на симуляторе → CASCI (химическая точность).
Self-test: H₂/def2-SVP CAS(2e,2o) — JW == FCI до ~1e-9 Ha (как у воркеров H₂/полиолефины).

Запуск: python3 -u routes/co2_to_fuels_vqe.py
"""
import os, json, time
import numpy as np
from pyscf import gto, scf, mcscf, ao2mo
from pyscf.mcscf import avas
import pennylane as qml

HERE = os.path.dirname(os.path.abspath(__file__))
RESULTS = os.path.join(HERE, "co2_to_fuels_results.json")
KCAL = 627.50947406  # Ha → kcal/mol


def cas_to_qubit_H(h1, eri, ecore, ncas):
    """Активные интегралы (chemist (pq|rs)) → ферми-оператор → Jordan–Wigner.
    Спин-орбиталь для пространственной p, спина s∈{0,1}:  2p+s."""
    op = ecore * qml.FermiWord({})                       # константа (ядро+кор)
    for p in range(ncas):
        for q in range(ncas):
            c = h1[p, q]
            if abs(c) > 1e-12:
                for s in (0, 1):
                    op += c * (qml.FermiC(2 * p + s) * qml.FermiA(2 * q + s))
    for p in range(ncas):
        for q in range(ncas):
            for r in range(ncas):
                for s in range(ncas):
                    c = eri[p, q, r, s]
                    if abs(c) > 1e-12:
                        for a in (0, 1):
                            for b in (0, 1):
                                op += 0.5 * c * (
                                    qml.FermiC(2 * p + a) * qml.FermiC(2 * r + b)
                                    * qml.FermiA(2 * s + b) * qml.FermiA(2 * q + a))
    H = qml.jordan_wigner(op)
    coeffs, ops = H.terms()                    # Hermitian → коэффициенты вещественны
    return qml.dot([float(np.real(c)) for c in coeffs], ops)


def ground_E_in_sector(H, nq, nelec):
    """Точная нижняя энергия кубитного H в секторе нужного числа электронов."""
    Hm = np.real(qml.matrix(H, wire_order=range(nq)))
    keep = [i for i in range(2 ** nq) if bin(i).count("1") == nelec]
    sub = Hm[np.ix_(keep, keep)]
    return float(np.linalg.eigvalsh(sub)[0])


def integrals_from_casci(mc):
    h1, ecore = mc.get_h1eff()
    eri = ao2mo.restore(1, mc.get_h2eff(), mc.ncas)
    return np.asarray(h1), np.asarray(eri), float(ecore)


def self_test_h2():
    print("=== self-test: H₂ CAS(2e,2o)/def2-SVP — JW vs FCI ===")
    mol = gto.M(atom="H 0 0 0; H 0 0 0.74", basis="def2-svp", verbose=0)
    mf = scf.RHF(mol); mf.verbose = 0; mf.kernel()
    mc = mcscf.CASCI(mf, 2, 2); mc.verbose = 0; mc.kernel()
    h1, eri, ec = integrals_from_casci(mc)
    H = cas_to_qubit_H(h1, eri, ec, mc.ncas)
    eq = ground_E_in_sector(H, 2 * mc.ncas, sum(mc.nelecas))
    print(f"  CASCI={mc.e_tot:.10f}  JW={eq:.10f}  Δ={(eq-mc.e_tot)*KCAL:.2e} kcal/mol")
    return abs(eq - mc.e_tot)


def occo_geom(R, dCO=1.13):
    return f"C 0 0 0; O 0 0 {dCO}; C {R} 0 0; O {R} 0 {dCO}"


def run_poc(R=1.6):
    print(f"\n=== *CO–*CO frontier CAS(4e,4o) @ C···C={R} Å — JW & VQE(UCCSD) ===")
    mol = gto.M(atom=occo_geom(R), basis="def2-svp", verbose=0)
    mf = scf.RHF(mol); mf.verbose = 0; mf.kernel()
    no, ne, orbs = avas.avas(mf, ["C 2px"], threshold=0.01, verbose=0)
    mc = mcscf.CASCI(mf, no, ne); mc.verbose = 0; mc.kernel(orbs)
    ncas, nelec = mc.ncas, sum(mc.nelecas)
    nq = 2 * ncas
    print(f"  active space CAS({nelec}e,{ncas}o) → {nq} qubits")

    # --- JW Hamiltonian exactly reproduces CASCI ---
    h1, eri, ec = integrals_from_casci(mc)
    H = cas_to_qubit_H(h1, eri, ec, ncas)
    n_terms = len(H.terms()[0]) if hasattr(H, "terms") else None
    eq = ground_E_in_sector(H, nq, nelec)
    dJW = (eq - mc.e_tot) * KCAL
    print(f"  JW terms={n_terms}  CASCI={mc.e_tot:.8f}  JW={eq:.8f}  Δ={dJW:.2e} kcal/mol")

    # --- VQE(UCCSD) on simulator ---
    from pennylane import numpy as pnp
    singles, doubles = qml.qchem.excitations(nelec, nq)
    s_wires, d_wires = qml.qchem.excitations_to_wires(singles, doubles)
    hf = qml.qchem.hf_state(nelec, nq)
    dev = qml.device("default.qubit", wires=nq)

    @qml.qnode(dev)
    def cost(params):
        qml.UCCSD(params, wires=range(nq), s_wires=s_wires, d_wires=d_wires, init_state=hf)
        return qml.expval(H)

    theta = pnp.zeros(len(singles) + len(doubles), requires_grad=True)
    opt = qml.GradientDescentOptimizer(stepsize=0.5)
    t0 = time.time(); e = float(cost(theta))
    for it in range(60):
        theta, e = opt.step_and_cost(cost, theta)
    dVQE = (float(e) - mc.e_tot) * KCAL
    print(f"  VQE(UCCSD): {len(theta)} params  E={float(e):.8f}  "
          f"→ CASCI gap={dVQE:+.3f} kcal/mol  ({time.time()-t0:.0f}s)")

    out = dict(geometry=f"*CO-*CO C···C={R} Å", cas=f"({nelec}e,{ncas}o)", qubits=nq,
               jw_terms=int(n_terms) if n_terms else None,
               e_casci=round(float(mc.e_tot), 8), e_jw=round(float(eq), 8),
               jw_minus_casci_kcal=float(f"{dJW:.2e}"),
               e_vqe=round(float(e), 8), vqe_gap_kcal=round(float(dVQE), 3))
    if os.path.exists(RESULTS):
        d = json.load(open(RESULTS)); d["vqe_poc"] = out
        json.dump(d, open(RESULTS, "w"), ensure_ascii=False, indent=1)
        print("  saved → results")
    return out


if __name__ == "__main__":
    err = self_test_h2()
    assert err < 1e-6, f"self-test failed: {err}"
    run_poc()
    print("\ndone.")
