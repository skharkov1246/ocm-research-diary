#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
routes/co2_to_fuels_fullcas.py — AWS-джоб: ПОЛНОЕ AVAS-пространство CAS(36e,22o).

Запись 2026-06-08 зафиксировала: полное AVAS(Cu 3d + C 2p + O 2p) на Cu₂+2*CO —
CAS(36e,22o) = 44 кубита, «квантово-релевантный» масштаб (как FeMoco/NiO). VQE-
statevector на 44q невозможен нигде (2^44 амплитуд); но КЛАССИЧЕСКИЙ CASCI в этом
пространстве — ~53.5 млн детерминантов — на пределе, но считается на r7i.4xlarge.
Это даёт: (1) честную верхнюю планку классики для этого узла; (2) NOON/n_u полного
пространства (мультиреференсность за пределами π-многообразия); (3) сравнение
барьера CASCI(36e,22o) с CAS(12e,10o) — сходится ли поправка по размеру пространства.

Протокол: геометрии реактанта/TS Cu₂ из ts_barrier (тот же скан), SCF — холодный
newton (у Cu₂ эндпойнтов одно SCF-решение, доказано в casscf_checks), AVAS с
принудительным окном (36e,22o). CASCI//PBE-орбитали БЕЗ оптимизации орбиталей —
на 22o CASSCF неподъёмен; оговорка та же, что у CASCI-этапов, и главный интерес
здесь NOON/масштаб, а не абсолютный барьер.

Запуск (AWS, ~128 ГБ RAM): python3 -u routes/co2_to_fuels_fullcas.py
"""
import os, sys, json, time
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from co2_to_fuels_sites import METALS, _mol, HARTREE_EV
from co2_to_fuels_ts import build_ads, _scf
from pyscf import mcscf, fci
from pyscf.mcscf import avas

HERE = os.path.dirname(os.path.abspath(__file__))
RESULTS = os.path.join(HERE, "co2_to_fuels_fullcas_results.json")
TS_SOURCE = os.path.join(HERE, "co2_to_fuels_results.json")
NCAS, NELECAS = 22, 36          # полное AVAS: Cu 3d + C 2p + O 2p (44 кубита)
MAXMEM_MB = int(os.environ.get("FULLCAS_MEM_MB", 110000))


def _persist(res):
    tmp = RESULTS + ".tmp"
    json.dump(res, open(tmp, "w"), ensure_ascii=False, indent=1)
    os.replace(tmp, RESULTS)


def endpoint(site, which):
    d = site[f"d_{which}"]
    free = [p["free"] for p in site["profile"] if p["d"] == d][0]
    return d, build_ads(d, np.array(free))


def run_point(tag, which, site, res):
    out = res.setdefault(tag, {}).setdefault(which, {})
    if out.get("e_casci") is not None:
        print(f"--- {tag}/{which}: done, skip ---"); return
    d, ads = endpoint(site, which)
    print(f"\n--- fullcas {tag}/{which} (d={d}) ---", flush=True)
    t0 = time.time()
    mol = _mol(METALS[tag], ads)
    mol.max_memory = MAXMEM_MB
    mf = _scf(mol)
    out["e_pbe"] = round(float(mf.e_tot), 6)
    out["scf_converged"] = bool(mf.converged)
    print(f"  SCF E={mf.e_tot:.6f} conv={mf.converged} ({time.time()-t0:.0f}s)",
          flush=True)

    no, ne, orbs = avas.avas(mf, ["Cu 3d", "C 2p", "O 2px", "O 2py"],
                             threshold=0.4, verbose=0)
    out["avas_native"] = f"({int(ne)}e,{int(no)}o)"
    mc = mcscf.CASCI(mf, NCAS, NELECAS)
    mc.max_memory = MAXMEM_MB
    mc.fcisolver = fci.direct_spin1.FCI(mol)
    mc.fcisolver.max_memory = MAXMEM_MB
    mc.fcisolver.max_cycle = 100
    mc.fcisolver.lindep = 1e-12
    mc.verbose = 4
    t1 = time.time()
    mc.kernel(orbs)
    dm = mc.fcisolver.make_rdm1(mc.ci, mc.ncas, mc.nelecas)
    noon = np.linalg.eigvalsh(dm)[::-1]
    nu = float(np.sum(np.minimum(noon, 2 - noon)))
    out.update(e_casci=round(float(mc.e_tot), 6),
               converged=bool(getattr(mc, "converged", True)),
               cas=f"({NELECAS}e,{NCAS}o)",
               ndet_est="C(22,18)^2 ~ 5.4e7",
               noon=[round(float(x), 4) for x in noon],
               n_u=round(nu, 4),
               wall_s=round(time.time() - t1, 1))
    _persist(res)
    print(f"  CASCI({NELECAS}e,{NCAS}o) E={mc.e_tot:.6f}  n_u={nu:.3f}  "
          f"({time.time()-t1:.0f}s)", flush=True)


if __name__ == "__main__":
    src = json.load(open(TS_SOURCE))
    res = json.load(open(RESULTS)) if os.path.exists(RESULTS) else {}
    res.setdefault("_meta", {
        "what": "full-AVAS CASCI(36e,22o)=44q on Cu2 scan endpoints (reactant, TS)",
        "why": "classical ceiling of the quantum-relevant space; NOON beyond the "
               "pi-manifold; barrier vs CAS(12e,10o) size-convergence",
        "honesty": "CASCI on PBE orbitals, NO orbital optimization (CASSCF "
                   "infeasible at 22o); frozen-metal cluster, def2-SVP, "
                   "no constant-potential; singlet"})
    site = src["ts_barrier"]["Cu2"]
    for which in ("reactant", "ts"):
        run_point("Cu2", which, site, res)
    cu2 = res.get("Cu2", {})
    if all(cu2.get(w, {}).get("e_casci") is not None for w in ("reactant", "ts")):
        cu2["barrier_casci_fullcas_eV"] = round(
            (cu2["ts"]["e_casci"] - cu2["reactant"]["e_casci"]) * HARTREE_EV, 3)
        cu2["barrier_pbe_eV"] = round(
            (cu2["ts"]["e_pbe"] - cu2["reactant"]["e_pbe"]) * HARTREE_EV, 3)
        ref = src.get("casscf_barrier", {}).get("Cu2", {})
        cu2["ref_barrier_casscf_12e10o_eV"] = ref.get("barrier_casscf_eV")
        _persist(res)
        print(f"\n=== Cu2 fullcas barrier: CASCI(36e,22o) "
              f"{cu2['barrier_casci_fullcas_eV']} eV  (PBE {cu2['barrier_pbe_eV']}, "
              f"CASSCF(12e,10o) {cu2['ref_barrier_casscf_12e10o_eV']}) ===")
    print("\nsaved -> results;  done.")
