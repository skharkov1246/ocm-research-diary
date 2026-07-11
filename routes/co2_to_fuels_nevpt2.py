#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
routes/co2_to_fuels_nevpt2.py — устойчивость CASSCF-решения + динамическая корреляция.

Надстройка над co2_to_fuels_casscf_crossval.py (тот даёт ΔE‡(CASSCF) с фикс. CAS(12e,10o)).
Три вопроса, каждый закрывает конкретную дыру в аргументации:

  A. TS CASSCF из НЕЗАВИСИМОГО AVAS-гесса (не из проекции реактанта): то же решение?
     → проекция не загнала TS в свой локальный минимум.
  B. CASSCF реактанта с ХОЛОДНОЙ SCF-ветки (та, что ~10 эВ выше warm — артефакт №2
     из casci_note): восстанавливает ли полная оптимизация орбиталей ТО ЖЕ
     CASSCF-решение? → чувствительность результата к выбору SCF-ветки, измеренная.
  C. SC-NEVPT2 поверх сошедшихся CASSCF (reactant, TS) на обоих сайтах:
     ΔE‡(CASSCF+NEVPT2) — динамическая корреляция; сравнение с литературной
     CASPT2-щелью уровня-в-уровень точнее, чем голый CASSCF.

Реализация: SCF-ветка воспроизводится тем же warm-chain (см. co2_to_fuels_casscf_crossval),
CASSCF пересобирается и СВЕРЯЕТСЯ с сохранёнными энергиями casscf_barrier (это
бесплатный тест воспроизводимости); NEVPT2 — через CASCI на сошедшихся CASSCF-
орбиталях (точные интегралы; сами орбитали из DF-CASSCF — стандартная практика,
помечено в results). Орбитали персистятся в routes/co2_to_fuels_casscf_mo.npz.

Ограничения: по-прежнему кластер (замороженный металл-димер, def2-SVP, без
constant-potential); NEVPT2 state-specific на двух точках профиля, не пересканирован
весь путь — барьер NEVPT2 наследует геометрии PBE-скана.

Запуск (после co2_to_fuels_casscf_crossval.py): python3 -u routes/co2_to_fuels_nevpt2.py
"""
import os, sys, json, time
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from co2_to_fuels_sites import METALS, _mol, HARTREE_EV, RESULTS
from co2_to_fuels_ts import build_ads, _scf
from co2_to_fuels_casscf_crossval import (NCAS, NELECAS, chain_scf, casscf_fixed,
                                 avas_guess, _persist)
from pyscf import mcscf, mrpt, dft
from pyscf.scf import addons as scf_addons

HERE = os.path.dirname(os.path.abspath(__file__))
MO_NPZ = os.path.join(HERE, "co2_to_fuels_casscf_mo.npz")


def _save_mo(tag, label, mc):
    data = dict(np.load(MO_NPZ)) if os.path.exists(MO_NPZ) else {}
    data[f"{tag}_{label}"] = mc.mo_coeff
    np.savez_compressed(MO_NPZ, **data)


def _nevpt2(mf, mc, tag_label):
    """SC-NEVPT2 через CASCI на сошедшихся CASSCF-орбиталях (точные интегралы)."""
    t0 = time.time()
    ci = mcscf.CASCI(mf, NCAS, NELECAS)
    ci.verbose = 0
    ci.kernel(mc.mo_coeff)
    dE = float(mrpt.NEVPT(ci).kernel())
    print(f"  NEVPT2 {tag_label}: E_casci={ci.e_tot:.6f}  dE_pt2={dE:.6f}  "
          f"({time.time()-t0:.0f}s)", flush=True)
    return float(ci.e_tot), dE


def run_site(tag, res):
    site = res["ts_barrier"][tag]
    base = res["casscf_barrier"][tag]
    chk = res.setdefault("casscf_checks", {}).setdefault(tag, {})
    n2 = res.setdefault("nevpt2_barrier", {}).setdefault(tag, {})
    if n2.get("barrier_nevpt2_eV") is not None:
        print(f"--- checks+NEVPT2 {tag}: already complete, skipping ---"); return
    print(f"\n--- checks + NEVPT2, site {tag} ---", flush=True)

    mfs, drift = chain_scf(tag, site)
    mf_r, mf_t = mfs[site["d_reactant"]], mfs[site["d_ts"]]
    # дозаполнить casscf_barrier, если запись сделана старой версией скрипта
    base.setdefault("scf_branch_drift_mHa", drift)

    saved = dict(np.load(MO_NPZ)) if os.path.exists(MO_NPZ) else {}

    # reactant CASSCF — сверка с сохранённым = тест воспроизводимости; если
    # орбитали уже в npz (краш-resume) — стартуем с них, иначе AVAS-гесс
    mo0_r = saved.get(f"{tag}_reactant", None)
    mc_r, _, _ = casscf_fixed(mf_r, mo0_r if mo0_r is not None
                              else avas_guess(mf_r)[0], f"{tag} reactant(re)")
    chk["reactant_repro_mHa"] = round(
        (float(mc_r.e_tot) - base["reactant_e_casscf"]) * 1e3, 3)
    _save_mo(tag, "reactant", mc_r)

    # TS CASSCF из проекции (как в основном прогоне) — сверка с сохранённым
    mo0_t = saved.get(f"{tag}_ts", None)
    if mo0_t is None:
        mo0_t = mcscf.project_init_guess(
            mcscf.CASSCF(mf_t, NCAS, NELECAS).density_fit(), mc_r.mo_coeff,
            prev_mol=mf_r.mol)
    mc_t, _, _ = casscf_fixed(mf_t, mo0_t, f"{tag} TS(projected,re)")
    chk["ts_repro_mHa"] = round((float(mc_t.e_tot) - base["ts_e_casscf"]) * 1e3, 3)
    _save_mo(tag, "ts", mc_t)
    _persist(res)

    # CHECK A: TS из независимого AVAS-гесса — тот же минимум CASSCF?
    if chk.get("A_ts_avas_minus_projected_mHa") is None:
        mc_t2, _, _ = casscf_fixed(mf_t, avas_guess(mf_t)[0], f"{tag} TS(avas)")
        chk["A_ts_avas_minus_projected_mHa"] = round(
            (float(mc_t2.e_tot) - float(mc_t.e_tot)) * 1e3, 3)
        chk["A_ts_avas_converged"] = bool(mc_t2.converged)
        _persist(res)

    # CHECK B: ПОИСК других SCF-решений в геометрии реактанта (у кластера их
    # несколько; для CuAl «холодный» newton совпадает с цепочкой, поэтому одной
    # холодной пробы мало) → CASSCF с самой удалённой ветки: то же решение?
    ads_r = build_ads(site["d_reactant"],
                      np.array([p["free"] for p in site["profile"]
                                if p["d"] == site["d_reactant"]][0]))
    branches = {"cold_newton": _scf(_mol(METALS[tag], ads_r))}
    mol_r = _mol(METALS[tag], ads_r)                  # smearing-anneal ветка
    mf_sm = dft.RKS(mol_r).density_fit(); mf_sm.xc = "pbe"; mf_sm.verbose = 0
    mf_sm = scf_addons.smearing_(mf_sm, sigma=0.01, method="fermi")
    mf_sm.max_cycle = 200; mf_sm.conv_tol = 1e-6; mf_sm.kernel()
    # newton() поверх smeared-объекта падает (дробные mo_occ в unpack_uniq_var)
    # → пере-затяжка чистым newton-RKS, warm-start от smeared-плотности
    branches["smearing_anneal"] = _scf(_mol(METALS[tag], ads_r),
                                       mf_sm.make_rdm1())
    mol_r2 = _mol(METALS[tag], ads_r)                 # level-shift/damp ветка
    mf_ls = dft.RKS(mol_r2).density_fit(); mf_ls.xc = "pbe"; mf_ls.verbose = 0
    mf_ls.level_shift = 0.4; mf_ls.damp = 0.3
    mf_ls.max_cycle = 300; mf_ls.conv_tol = 1e-6; mf_ls.kernel()
    branches["level_shift"] = _scf(_mol(METALS[tag], ads_r),
                                   mf_ls.make_rdm1())
    gaps = {k: round((float(m.e_tot) - float(mf_r.e_tot)) * HARTREE_EV, 3)
            for k, m in branches.items()}
    chk["B_scf_branches_vs_chain_eV"] = gaps
    alt = max(gaps, key=lambda k: abs(gaps[k]))
    chk["B_alt_branch"] = alt
    if abs(gaps[alt]) < 0.005:
        chk["B_note"] = ("no distinct SCF solution found by cold/smearing/"
                         "level-shift probes at the reactant geometry")
        print(f"  checks: repro r/ts={chk['reactant_repro_mHa']}/{chk['ts_repro_mHa']} mHa;  "
              f"A(avas-proj)={chk['A_ts_avas_minus_projected_mHa']} mHa;  "
              f"B: no distinct branch ({gaps})", flush=True)
    else:
        mc_alt, _, _ = casscf_fixed(branches[alt], avas_guess(branches[alt])[0],
                                    f"{tag} reactant({alt})")
        chk["B_alt_casscf_minus_chain_mHa"] = round(
            (float(mc_alt.e_tot) - float(mc_r.e_tot)) * 1e3, 3)
        chk["B_alt_casscf_converged"] = bool(mc_alt.converged)
        print(f"  checks: repro r/ts={chk['reactant_repro_mHa']}/{chk['ts_repro_mHa']} mHa;  "
              f"A(avas-proj)={chk['A_ts_avas_minus_projected_mHa']} mHa;  "
              f"B: SCF branch '{alt}' gap={gaps[alt]} eV → "
              f"CASSCF gap={chk['B_alt_casscf_minus_chain_mHa']} mHa", flush=True)
    _persist(res)

    # CHECK C / NEVPT2: динамическая корреляция на обоих концах
    e_r_ci, dpt_r = _nevpt2(mf_r, mc_r, f"{tag} reactant")
    n2.update(reactant_e_casci=round(e_r_ci, 6), reactant_pt2=round(dpt_r, 6))
    _persist(res)
    e_t_ci, dpt_t = _nevpt2(mf_t, mc_t, f"{tag} TS")
    n2.update(ts_e_casci=round(e_t_ci, 6), ts_pt2=round(dpt_t, 6),
              barrier_casscf_eV=base["barrier_casscf_eV"],
              barrier_nevpt2_eV=round(
                  ((e_t_ci + dpt_t) - (e_r_ci + dpt_r)) * HARTREE_EV, 3))
    n2["quantum_corr_vs_pbe_eV"] = round(
        n2["barrier_nevpt2_eV"] - base["barrier_pbe_eV"], 3)
    _persist(res)
    print(f"  ΔE‡: PBE={base['barrier_pbe_eV']:+.3f}  CASSCF={base['barrier_casscf_eV']:+.3f}  "
          f"NEVPT2={n2['barrier_nevpt2_eV']:+.3f} eV  "
          f"(corr vs PBE {n2['quantum_corr_vs_pbe_eV']:+.3f})", flush=True)


if __name__ == "__main__":
    res = json.load(open(RESULTS))
    for tag in ("Cu2", "CuAl"):
        if res.get("casscf_barrier", {}).get(tag, {}).get("barrier_casscf_eV") is None:
            print(f"--- {tag}: casscf_barrier not ready, skip (run co2_to_fuels_casscf.py) ---")
            continue
        run_site(tag, res)
    nb = res.get("nevpt2_barrier", {})
    if all(nb.get(t, {}).get("barrier_nevpt2_eV") is not None for t in ("Cu2", "CuAl")):
        po = sorted(("Cu2", "CuAl"),
                    key=lambda t: res["casscf_barrier"][t]["barrier_pbe_eV"])
        no_ = sorted(("Cu2", "CuAl"), key=lambda t: nb[t]["barrier_nevpt2_eV"])
        nb["_summary"] = dict(
            pbe_order=list(po), nevpt2_order=list(no_),
            order_flipped=bool(list(po) != list(no_)),
            note=("SC-NEVPT2 on converged fixed-CAS(12e,10o) CASSCF states (CASCI with "
                  "exact integrals on DF-CASSCF orbitals), same warm-chain SCF branch and "
                  "PBE-scan geometries. casscf_checks: A = TS CASSCF from independent AVAS "
                  "guess vs projected (local-minimum test), B = reactant CASSCF from the "
                  "cold SCF branch (SCF-multiple-solution sensitivity, measured). Frozen-"
                  "metal cluster, def2-SVP, no constant-potential; state-specific NEVPT2 "
                  "at two scan points, not a rescanned path."))
        _persist(res)
        print(f"\n=== NEVPT2 ORDERING ===\n  PBE:{po}  NEVPT2:{no_}  "
              f"flipped={nb['_summary']['order_flipped']}")
    print("\nsaved → results;  done.")
