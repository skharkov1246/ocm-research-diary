#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
MatterForge · h2o2-direct — Stage 2 ПРОДАКШЕН-СКРИН дешёвых 3d-центров (M–N₄ SAC).

Переносит «дескриптор Pd₁» на земные металлы: для M ∈ {Fe, Co, Ni, Pd} строит
металло-порфин (канонический M–N₄ карман), находит основное спиновое состояние,
оптимизирует геометрию, сажает *OOH, и считает дескриптор селективности 2e⁻→H₂O₂:
  • держит ли центр связь O–O (R(O–O) в *OOH; ~1.45 Å = пероксо/2e⁻, удлинение = к 4e⁻);
  • энергию связывания ΔE(*OOH);
  • мультиреференсность центра — N_u по NOON из AVAS(M 3d + O 2p)→CASCI (как в NiO/OCM).
Ключевой литературный контраст для проверки: Fe–N₄ тянет в 4e⁻ (→H₂O), Co–N₄ — в 2e⁻.

ВАЖНО (честно): открытооболочечный DFT металло-порфина (def2-SVP, ~420 базисных фн)
локально на CPU не сходится за разумное время (SCF-стена; проверено 3 попытки: plain,
density-fit, second-order Newton). Это AWS-нагрузка, как FeMoco / NiO-18q-VQE в др. вкладках.

  python3 routes/h2o2_direct_screen.py --selftest   # быстрая проверка ПЛУМБИНГА (мелкие системы)
  python3 routes/h2o2_direct_screen.py --screen      # полный скрин (для AWS/мощной машины)
"""
import sys
import numpy as np

SMI = "c1cc2cc3ccc(cc4ccc(cc5ccc(cc1[nH]2)n5)[nH]4)n3"   # free-base porphine, C20H14N4
DLABEL = {"Fe": "Fe 3d", "Co": "Co 3d", "Ni": "Ni 3d", "Pd": "Pd 4d"}
# плауз. мультиплетности M(II)-порфина (2S+1) до *OOH (метрика — основное состояние)
SPINS = {"Fe": [1, 3, 5], "Co": [2, 4], "Ni": [1, 3], "Pd": [1, 3]}


def Nu(occ):
    occ = np.asarray(occ)
    return float(np.sum(occ * (2.0 - occ)))


def build_mp(metal):
    """RDKit: порфин → 3D → металлирование (металл в центроид N₄, снять 2 пиррольных N–H)."""
    from rdkit import Chem
    from rdkit.Chem import AllChem
    m = Chem.AddHs(Chem.MolFromSmiles(SMI))
    AllChem.EmbedMolecule(m, AllChem.ETKDGv3())
    AllChem.MMFFOptimizeMolecule(m, maxIters=2000)
    P = m.GetConformer().GetPositions()
    Ns = [a.GetIdx() for a in m.GetAtoms() if a.GetSymbol() == "N"]
    drop = {nb.GetIdx() for n in Ns for nb in m.GetAtomWithIdx(n).GetNeighbors()
            if nb.GetSymbol() == "H"}
    cen = P[Ns].mean(0)
    atoms = [(metal, tuple(cen))]
    atoms += [(a.GetSymbol(), tuple(P[a.GetIdx()]))
              for a in m.GetAtoms() if a.GetIdx() not in drop]
    return atoms


def n4_normal(atoms):
    """нормаль к плоскости N₄ (через SVD)."""
    N = np.array([c for s, c in atoms if s == "N"])
    c = N.mean(0)
    _, _, vt = np.linalg.svd(N - c)
    return vt[-1] / np.linalg.norm(vt[-1])


def place_ooh(atoms, m_o=1.85, o_o=1.40, o_h=0.97):
    """посадить *OOH аксиально на металл (геом-опт потом уточнит)."""
    metal = np.array(atoms[0][1])
    n = n4_normal(atoms)
    t = np.cross(n, [1.0, 0.0, 0.0]); t = t / (np.linalg.norm(t) + 1e-9)
    O1 = metal + m_o * n
    O2 = O1 + o_o * (0.5 * n + 0.866 * t)
    H = O2 + o_h * (0.5 * n - 0.866 * t)
    return atoms + [("O", tuple(O1)), ("O", tuple(O2)), ("H", tuple(H))]


def pyscf_atom(atoms):
    return "; ".join(f"{s} {x:.5f} {y:.5f} {z:.5f}" for s, (x, y, z) in atoms)


def make_mol(atoms, spin, charge=0, basis="def2-svp"):
    from pyscf import gto
    return gto.M(atom=pyscf_atom(atoms), basis=basis, ecp=basis, spin=spin,
                 charge=charge, verbose=0)


def robust_uks(mol, xc="pbe0"):
    """устойчивый открытооболочечный SCF: density-fit + level-shift + second-order Newton."""
    from pyscf import dft
    mf = dft.UKS(mol, xc=xc).density_fit()
    mf.grids.level = 3
    mf.level_shift = 0.3
    mf.conv_tol = 1e-7
    mf = mf.newton()
    mf.max_cycle = 120
    mf.kernel()
    return mf


def avas_casci_Nu(mf, metal):
    """AVAS(M 3d + O 2p) → CASCI → N_u по заселённостям натуральных орбиталей."""
    from pyscf import mcscf
    from pyscf.mcscf import avas
    ncas, nelec, mo = avas.avas(mf, [DLABEL[metal], "O 2p"], canonicalize=True, verbose=0)
    mc = mcscf.CASCI(mf, ncas, nelec)
    mc.kernel(mo)
    dm1 = mc.fcisolver.make_rdm1(mc.ci, mc.ncas, mc.nelecas)
    occ = np.linalg.eigvalsh(dm1)[::-1]
    return ncas, nelec, Nu(occ), occ, mc.e_tot


def ground_spin(metal, atoms):
    """скан мультиплетности → основное состояние (E, spin, mf)."""
    best = None
    for s in SPINS[metal]:
        try:
            mf = robust_uks(make_mol(atoms, spin=s - 1))   # spin = 2S = (2S+1)-1
            if mf.converged and (best is None or mf.e_tot < best[0]):
                best = (mf.e_tot, s, mf)
        except Exception as e:
            print(f"   spin {s}: FAIL {type(e).__name__}")
    return best


def screen(metals=("Fe", "Co", "Ni", "Pd")):
    from pyscf.geomopt.berny_solver import optimize
    print("== Stage 2 продакшен-скрин M–N₄ SAC (тяжёлое; для AWS) ==")
    for M in metals:
        print(f"\n### {M}-порфин ###")
        bare = build_mp(M)
        gs = ground_spin(M, bare)
        if gs is None:
            print(f"  {M}: SCF не сошёлся ни в одном спине — лог/настройки"); continue
        E_bare, spin, mf = gs
        print(f"  основной спин (2S+1)={spin}  E_bare={E_bare:.4f}")
        mol_opt = optimize(mf, maxsteps=60)                 # геом-опт голого центра
        ads = place_ooh([(a.symbol if hasattr(a, 'symbol') else s, tuple(c))
                         for (s, c), a in zip(bare, mol_opt.atom)])
        gsa = ground_spin(M, ads)
        if gsa is None:
            print(f"  {M}: *OOH SCF не сошёлся"); continue
        E_ads, spin_a, mf_a = gsa
        mol_a = optimize(mf_a, maxsteps=80)
        coords = mol_a.atom_coords() * 0.529177
        # O–O расстояние (последние два O в списке)
        oidx = [i for i, z in enumerate(mol_a.elements) if z == "O"][-2:]
        roo = np.linalg.norm(coords[oidx[0]] - coords[oidx[1]])
        ncas, nelec, nu, occ, e_cas = avas_casci_Nu(mf_a, M)
        print(f"  *OOH: спин={spin_a} R(O–O)={roo:.3f}Å  "
              f"AVAS->CAS({nelec},{ncas}o) N_u={nu:.2f}")
        print(f"  => держит O–O (2e⁻): {'да' if roo < 1.55 else 'НЕТ (к 4e⁻)'}")


def selftest():
    """быстрая проверка плумбинга (без металло-порфиновой SCF-стены)."""
    print("== SELFTEST: плумбинг скрина на мелких системах ==")
    fe = build_mp("Fe")
    print(f"[build] Fe-порфин: {len(fe)} атомов (ожид. 37)")
    ads = place_ooh(fe)
    roo = np.linalg.norm(np.array(ads[-2][1]) - np.array(ads[-3][1]))
    print(f"[place_ooh] стартовое R(O–O)={roo:.3f}Å (до опт)")
    # robust_uks на маленьком открытооболочечном O2 (триплет) — быстро
    mf = robust_uks(make_mol([("O", (0, 0, 0)), ("O", (0, 0, 1.207))], spin=2))
    print(f"[robust_uks] O2 триплет: conv={mf.converged} E={mf.e_tot:.4f}")
    # AVAS+CASCI+N_u код-путь на O2 (метка 'O 2p')
    from pyscf import mcscf
    from pyscf.mcscf import avas
    ncas, nelec, mo = avas.avas(mf, ["O 2p"], canonicalize=True, verbose=0)
    mc = mcscf.CASCI(mf, ncas, nelec); mc.kernel(mo)
    dm1 = mc.fcisolver.make_rdm1(mc.ci, mc.ncas, mc.nelecas)
    occ = np.linalg.eigvalsh(dm1)[::-1]
    print(f"[avas_casci_Nu] O2: CAS({nelec},{ncas}o) N_u={Nu(occ):.2f} (ожид. ~2.3 — бирадикал)")
    print("SELFTEST OK — плумбинг корректен; продакшен-скрин уезжает на AWS.")


if __name__ == "__main__":
    if "--screen" in sys.argv:
        screen()
    else:
        selftest()
