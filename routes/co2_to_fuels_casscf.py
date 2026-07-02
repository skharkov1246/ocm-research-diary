#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
routes/co2_to_fuels_casscf.py — квантовая поправка к барьеру сочленения, ЧЕСТНАЯ версия.

Чинит два артефакта CASCI//PBE из co2_to_fuels_ts.py (см. casci_note в results-JSON):
  1) AVAS давал РАЗНЫЕ активные пространства на разных геометриях (12e↔14e, 11o↔10o)
     → здесь активное пространство ФИКСИРОВАНО: CAS(12e,10o) на всех геометриях,
     а на TS начальные орбитали ПРОЕЦИРУЮТСЯ из сошедшегося CASSCF реактанта
     (mcscf.project_init_guess) — одинаковый характер пространства на обоих концах.
  2) SCF металл-кластера имеет несколько решений (холодный рестарт ≈10 эВ от ветки
     релаксации) → здесь SCF-ветка ВОСПРОИЗВОДИТСЯ тем же протоколом, что в скане:
     холодный старт на d=1.8 Å, затем warm-chain по релакс-геометриям профиля до
     d=4.2 Å; на каждой точке E(PBE) сверяется с сохранённой (drift в mHa пишется
     в результаты — это контроль, что мы на ТОЙ ЖЕ ветке).
  3) CASCI брал орбитали DFT как есть; CASSCF оптимизирует орбитали вариационно —
     энергия не зависит от произвола одноэлектронной ветки в активном пространстве.

Считается: E_CASSCF(reactant), E_CASSCF(TS) на обоих сайтах (Cu₂, CuAl),
ΔE‡(CASSCF), Δ_corr = ΔE‡(CASSCF) − ΔE‡(PBE), NOON/n_u(π) в седле, порядок сайтов.

ЧЕСТНО: по-прежнему жёсткая модель — замороженный металл-димер, def2-SVP, нет
constant-potential; CASSCF без динамической корреляции (NEVPT2 — следующий шаг).
Это поправка к БАРЬЕРУ кластерной модели, не фарадеевская эффективность.

Запуск: python3 -u routes/co2_to_fuels_casscf.py
"""
import os, sys, json, time
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from co2_to_fuels_sites import METALS, _mol, HARTREE_EV, RESULTS
from co2_to_fuels_ts import build_ads, _scf
from pyscf import mcscf
from pyscf.mcscf import avas

NCAS, NELECAS = 10, 12                      # фиксированное π-многообразие CAS(12e,10o)


HERE = os.path.dirname(os.path.abspath(__file__))
STATE_NPZ = os.path.join(HERE, "co2_to_fuels_casscf_state.npz")


def _persist(res):
    json.dump(res, open(RESULTS, "w"), ensure_ascii=False, indent=1)


def _save_npz(tag, **arrs):
    data = dict(np.load(STATE_NPZ)) if os.path.exists(STATE_NPZ) else {}
    data.update({f"{tag}_{k}": v for k, v in arrs.items()})
    np.savez_compressed(STATE_NPZ, **data)


def _noon(mc):
    dm = mc.fcisolver.make_rdm1(mc.ci, mc.ncas, mc.nelecas)
    n = np.linalg.eigvalsh(dm)[::-1]
    nu = float(np.sum(np.minimum(n, 2 - n)))
    return [round(float(x), 4) for x in n], round(nu, 4)


def chain_scf(tag, site):
    """Warm-chain SCF по релакс-профилю (d=1.8→4.2, как в оригинальном скане).
    Возвращает {d: mf} для reactant/TS + drift энергии от сохранённой ветки
    НА КАЖДОЙ точке (подписанный, мГа; <0 = цепочка НИЖЕ сохранённой)."""
    metal = METALS[tag]
    want = {site["d_reactant"], site["d_ts"]}
    keep, dm, drift = {}, None, {}
    for p in sorted(site["profile"], key=lambda q: q["d"]):
        t0 = time.time()
        mf = _scf(_mol(metal, build_ads(p["d"], np.array(p["free"]))), dm)
        dm = mf.make_rdm1()
        dr = (float(mf.e_tot) - p["e_pbe"]) * 1e3             # mHa, подписанный
        drift[f"{p['d']:.1f}"] = round(dr, 3)
        print(f"  chain d={p['d']:.2f}  E={mf.e_tot:.6f}  drift={dr:+.3f} mHa  "
              f"({time.time()-t0:.0f}s)", flush=True)
        if p["d"] in want:
            keep[p["d"]] = mf
    return keep, drift


def casscf_fixed(mf, mo0, tag_label, restarts=4):
    """DF-CASSCF с фикс. CAS; при несходимости — рестарт ядра с последних
    орбиталей (свежая AH-история), до `restarts` раз."""
    mc = mcscf.CASSCF(mf, NCAS, NELECAS).density_fit()
    mc.verbose = 0; mc.conv_tol = 1e-7; mc.max_cycle_macro = 120
    t0 = time.time()
    mc.kernel(mo0)
    tries = 1
    while not mc.converged and tries < restarts:
        tries += 1
        print(f"  CASSCF {tag_label}: not converged, restart #{tries} "
              f"from current orbitals", flush=True)
        mc.kernel(mc.mo_coeff)
    noon, nu = _noon(mc)
    print(f"  CASSCF {tag_label}: E={mc.e_tot:.6f}  conv={mc.converged}  "
          f"n_u={nu:.3f}  ({tries} run(s), {time.time()-t0:.0f}s)", flush=True)
    return mc, noon, nu


def avas_guess(mf):
    """AVAS только как НАЧАЛЬНЫЕ орбитали; размер пространства всё равно фиксирован."""
    no, ne, orbs = avas.avas(mf, ["C 2p", "O 2px", "O 2py"], threshold=0.4, verbose=0)
    return orbs, (int(ne), int(no))


def run_site(tag, res):
    site = res["ts_barrier"][tag]
    out = res.setdefault("casscf_barrier", {}).setdefault(tag, {})
    if out.get("barrier_casscf_eV") is not None:
        print(f"--- CASSCF {tag}: already complete, skipping ---"); return out
    print(f"\n--- CASSCF barrier, site {tag} ---", flush=True)

    mfs, drift = chain_scf(tag, site)
    out["scf_branch_drift_mHa"] = drift
    out["scf_branch_max_drift_mHa"] = max(abs(v) for v in drift.values())
    mf_r, mf_t = mfs[site["d_reactant"]], mfs[site["d_ts"]]
    # ДВА PBE-барьера честно: (а) ветка ЭТОЙ цепочки — та же, на которой стоит
    # CASSCF (у CuAl совпала со старыми «холодными» эндпойнтами, 1.053);
    # (б) релакс-профиль скана (вариационно ниже у CuAl: warm 0.992).
    # Чувствительность CASSCF к выбору ветки меряет Check B (co2_to_fuels_nevpt2).
    out["barrier_pbe_eV"] = round(
        (float(mf_t.e_tot) - float(mf_r.e_tot)) * HARTREE_EV, 3)
    prof_e = {p["d"]: p["e_pbe"] for p in site["profile"]}
    out["barrier_pbe_profile_eV"] = round(
        (prof_e[site["d_ts"]] - prof_e[site["d_reactant"]]) * HARTREE_EV, 3)
    # персист dm/mo эндпойнтов — тёплый рестарт для NEVPT2/checks без полной цепочки
    _save_npz(tag, reactant_dm=mf_r.make_rdm1(), ts_dm=mf_t.make_rdm1())

    # reactant: AVAS-guess (размер зафиксирован NCAS/NELECAS независимо от AVAS)
    orbs_r, avas_dims_r = avas_guess(mf_r)
    mc_r, noon_r, nu_r = casscf_fixed(mf_r, orbs_r, f"{tag} reactant")
    out.update(avas_dims_reactant=f"({avas_dims_r[0]}e,{avas_dims_r[1]}o)",
               reactant_e_casscf=round(float(mc_r.e_tot), 6),
               reactant_converged=bool(mc_r.converged),
               reactant_noon=noon_r, reactant_nu=nu_r)
    _save_npz(tag, reactant_mo=mc_r.mo_coeff)
    _persist(res)

    # TS: начальные орбитали = ПРОЕКЦИЯ сошедшегося CASSCF реактанта (то же пространство)
    mo_ts = mcscf.project_init_guess(mcscf.CASSCF(mf_t, NCAS, NELECAS).density_fit(),
                                     mc_r.mo_coeff, prev_mol=mf_r.mol)
    mc_t, noon_t, nu_t = casscf_fixed(mf_t, mo_ts, f"{tag} TS")
    out.update(ts_e_casscf=round(float(mc_t.e_tot), 6),
               ts_converged=bool(mc_t.converged),
               ts_noon=noon_t, ts_nu=nu_t)
    _save_npz(tag, ts_mo=mc_t.mo_coeff)

    out["barrier_casscf_eV"] = round(
        (out["ts_e_casscf"] - out["reactant_e_casscf"]) * HARTREE_EV, 3)
    out["quantum_corr_eV"] = round(out["barrier_casscf_eV"] - out["barrier_pbe_eV"], 3)
    _persist(res)
    print(f"  ΔE‡: PBE={out['barrier_pbe_eV']:+.3f}  CASSCF={out['barrier_casscf_eV']:+.3f}  "
          f"Δ_corr={out['quantum_corr_eV']:+.3f} eV   TS n_u={out['ts_nu']:.3f}", flush=True)
    return out


if __name__ == "__main__":
    res = json.load(open(RESULTS))
    for tag in ("Cu2", "CuAl"):
        run_site(tag, res)
    cb = res["casscf_barrier"]
    if all(cb.get(t, {}).get("barrier_casscf_eV") is not None for t in ("Cu2", "CuAl")):
        po = sorted(("Cu2", "CuAl"), key=lambda t: cb[t]["barrier_pbe_eV"])
        co = sorted(("Cu2", "CuAl"), key=lambda t: cb[t]["barrier_casscf_eV"])
        cb["_summary"] = dict(
            pbe_order=list(po), casscf_order=list(co),
            order_flipped=bool(list(po) != list(co)),
            note=("Fixed CAS(12e,10o) CASSCF (orbital-optimized, DF) on the reproduced "
                  "warm-chain SCF branch; TS orbitals projected from converged reactant "
                  "CASSCF → consistent active space, no per-geometry AVAS artifact. "
                  "Frozen-metal cluster, def2-SVP, no constant-potential, no dynamic "
                  "correlation (NEVPT2 = next); cluster-model barrier, not Faradaic "
                  "efficiency. Supersedes barrier_casci_eV (casci_barrier_reliable=false)."))
        _persist(res)
        print(f"\n=== BARRIER ORDERING ===\n  PBE:{po}  CASSCF:{co}  "
              f"flipped={cb['_summary']['order_flipped']}")
    print("\nsaved → results;  done.")
