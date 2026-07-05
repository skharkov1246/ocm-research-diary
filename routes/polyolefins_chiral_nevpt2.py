#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""routes/polyolefins_chiral_nevpt2.py — квантовая поправка к стерео-ΔΔE‡ (issue #11,
последний пункт): CASSCF(6,6) + SC-NEVPT2 на TS-геометриях ОБЕИХ граней хирального
rac-C2 ansa-титаноценового кармана из polyolefins_chiral_results.json / *_ts_{si,re}.xyz
(релаксированные точки cc=2.0 прогона 2026-07-05, Франкфурт).

Зачем: CASCI(6,6)//RHF в том прогоне дал ΔΔE_TS −0.328 vs DFT −1.858 ккал/моль
(расхождение 1.53) — но без орбитальной оптимизации и динамической корреляции это
ДИАГНОСТИКА. Здесь: (1) CASSCF — орбитально-оптимизированный, (2) SC-NEVPT2 —
динамическая корреляция, (3) ОДИНАКОВОЕ пространство обеих граней: AVAS на si
задаёт (6e,6o), re стартует с project_init_guess CASSCF-орбиталей si (урок
co2_to_fuels_casscf.py).

ЧЕСТНО: TS-геометрии — DFT(PBE)-релаксация констрейн-скана (не переоптимизированы
на CASSCF-уровне); лиганд был заморожен; def2-SVP; газофазный кластер. ΔΔE‡ здесь —
вертикальная поправка на DFT-геометриях. Знак конвенции: ΔΔE=E(si)−E(re).

Запуск: OMP_NUM_THREADS=$(nproc) python3 -u routes/polyolefins_chiral_nevpt2.py
Выход:  routes/polyolefins_chiral_nevpt2_results.json (atomic, incremental)
"""
import os, json, time, tempfile
import numpy as np
from pyscf import gto, scf, mcscf, fci, mrpt
from pyscf.mcscf import avas as avas_mod, addons as mc_addons

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "polyolefins_chiral_nevpt2_results.json")
HARTREE_KCAL = 627.5094740631

def say(m): print(f"[{time.strftime('%H:%M:%S')}] {m}", flush=True)

def atomic_save(obj):
    fd, tmp = tempfile.mkstemp(dir=HERE, suffix=".tmp")
    with os.fdopen(fd, "w") as f:
        json.dump(obj, f, ensure_ascii=False, indent=1)
    os.replace(tmp, OUT)

def load_xyz(face):
    path = os.path.join(HERE, f"polyolefins_chiral_ts_{face}.xyz")
    lines = open(path).read().splitlines()
    n = int(lines[0]); comment = lines[1]
    atoms = [ln.split()[:4] for ln in lines[2:2 + n]]
    return [[a[0], (float(a[1]), float(a[2]), float(a[3]))] for a in atoms], comment

def build(face):
    atoms, comment = load_xyz(face)
    mol = gto.M(atom=atoms, basis="def2-svp", charge=1, spin=0, verbose=0,  # катионный металлоцен: референс model.nelec=154 (155 Z - 1)
                max_memory=float(os.environ.get("PYSCF_MAX_MEMORY", 40000)))
    return mol, comment

res = json.load(open(OUT)) if os.path.exists(OUT) else {}
res.setdefault("meta", {
    "what": "CASSCF(6,6)+SC-NEVPT2 vertical stereo-ddE at the DFT TS geometries "
            "(cc=2.0) of both faces; single active space defined on si, re via "
            "project_init_guess",
    "geoms": "polyolefins_chiral_ts_{si,re}.xyz (relaxed DFT PBE constrained scan, "
             "2026-07-05 Frankfurt run; ligand frozen there)",
    "basis": "def2-svp", "convention": "ddE = E(si) - E(re); mirror ligand => -ddE",
    "honesty": "vertical correction on DFT geometries; no CASSCF re-relaxation; "
               "gas-phase cluster; reference CASCI diagnostic gave ddE_TS -0.328 "
               "vs DFT -1.858 kcal/mol",
})
atomic_save(res)

space = None      # (ncas, nelec) fixed on si
mo_si = None
mol_si = None
for face in ["si", "re"]:
    blk = res.setdefault(face, {})
    if blk.get("e_nevpt2") is not None:
        say(f"{face}: already done (restart), E_nevpt2={blk['e_nevpt2']}")
        if face == "si":
            space = (blk["ncas"], blk["nelecas"])
        continue
    mol, comment = build(face)
    blk["geom_comment"] = comment
    say(f"{face}: RHF/DF nao={mol.nao} …")
    t0 = time.time()
    mf = scf.RHF(mol).density_fit()
    mf.conv_tol = 1e-9; mf.max_cycle = 200
    mf.kernel()
    if not mf.converged:
        mfn = mf.newton(); mfn.max_cycle = 60
        mfn.kernel(mf.mo_coeff, mf.mo_occ)
        mf.mo_coeff, mf.mo_occ, mf.mo_energy = mfn.mo_coeff, mfn.mo_occ, mfn.mo_energy
        mf.e_tot, mf.converged = float(mfn.e_tot), bool(mfn.converged)
    blk.update(e_hf=float(mf.e_tot), hf_converged=bool(mf.converged),
               t_hf_s=round(time.time() - t0, 1))
    say(f"{face}: E_HF={mf.e_tot:.6f} conv={mf.converged} ({blk['t_hf_s']}s)")
    atomic_save(res)

    if face == "si":
        ncas, nelec, orbs = avas_mod.avas(
            mf, ["Ti 3d", "28 C 2p", "32 C 2p", "33 C 2p"], threshold=0.4,
            canonicalize=True)
        # усечение до (6,6) NOON-подобным правилом невозможно тут заранее —
        # повторяем протокол исходного прогона: берём 6 орбиталей вокруг Ферми
        # из AVAS-окна через CASSCF-стартер (cas_used был [6,6])
        ncas0, nelec0 = int(ncas), int(nelec)
        say(f"si: AVAS raw ({nelec0}e,{ncas0}o); fixing used space to (6e,6o) "
            f"как в референс-прогоне")
        space = (6, 6)
        mc = mcscf.CASSCF(mf, 6, 6)
        # выбрать 6 орбиталей AVAS ближе всего к Ферми: sort_mo по окну AVAS
        mc.max_cycle_macro = 80; mc.conv_tol = 1e-7
        try: fci.addons.fix_spin_(mc.fcisolver, ss=0)
        except Exception: pass
        t0 = time.time()
        mc.kernel(orbs)
        mo_si, mol_si = mc.mo_coeff, mol
    else:
        mc = mcscf.CASSCF(mf, space[0], space[1])
        mc.max_cycle_macro = 80; mc.conv_tol = 1e-7
        try: fci.addons.fix_spin_(mc.fcisolver, ss=0)
        except Exception: pass
        guess = mc_addons.project_init_guess(mc, mo_si, prev_mol=mol_si)
        t0 = time.time()
        mc.kernel(guess)
    dm1 = mc.fcisolver.make_rdm1(mc.ci, mc.ncas, mc.nelecas)
    noon = np.linalg.eigvalsh(dm1)[::-1]
    n_u = float(np.sum(np.minimum(noon, 2 - noon)))
    ne = int(sum(mc.nelecas)) if np.ndim(mc.nelecas) else int(mc.nelecas)
    blk.update(e_casscf=float(mc.e_tot), casscf_converged=bool(mc.converged),
               ncas=int(mc.ncas), nelecas=ne,
               noon=[round(float(x), 4) for x in noon], n_u=round(n_u, 4),
               t_casscf_s=round(time.time() - t0, 1))
    say(f"{face}: E_CASSCF={mc.e_tot:.6f} conv={mc.converged} n_u={n_u:.3f} "
        f"({blk['t_casscf_s']}s)")
    atomic_save(res)

    say(f"{face}: SC-NEVPT2 …")
    t0 = time.time()
    try:
        ept2 = float(mrpt.NEVPT(mc).kernel())
        blk.update(nevpt2_corr=ept2, e_nevpt2=float(mc.e_tot) + ept2,
                   t_nevpt2_s=round(time.time() - t0, 1))
        say(f"{face}: E_NEVPT2={blk['e_nevpt2']:.6f} (corr {ept2:.6f}, "
            f"{blk['t_nevpt2_s']}s)")
    except Exception as exc:
        blk.update(nevpt2_error=f"{type(exc).__name__}: {exc}")
        say(f"{face}: NEVPT2 FAILED: {exc}")
    atomic_save(res)

if all(res.get(f, {}).get("e_nevpt2") is not None for f in ("si", "re")):
    dd = lambda k: round((res["si"][k] - res["re"][k]) * HARTREE_KCAL, 3)
    res["stereo_nevpt2"] = {
        "ddE_ts_hf_kcal": dd("e_hf"),
        "ddE_ts_casscf_kcal": dd("e_casscf"),
        "ddE_ts_nevpt2_kcal": dd("e_nevpt2"),
        "reference": {"ddE_ts_dft_kcal": -1.858, "ddE_ts_casci_kcal": -0.328,
                      "source": "polyolefins_chiral_results.json (2026-07-05)"},
        "note": "vertical at DFT TS geometries; ddE=E(si)-E(re)",
    }
    res["status"] = "ok"
else:
    res["status"] = "incomplete"
atomic_save(res)
say(f"status={res['status']} -> {OUT}")
