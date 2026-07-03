#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
routes/co2_to_fuels_nevpt2_aws.py — SC-NEVPT2 поверх CASSCF-реабилитации (issue #13).

Динамическая корреляция к поправке барьера C–C сочетания. Референс — РОВНО тот же
протокол, что в routes/co2_to_fuels_casscf.py (функции ИМПОРТИРУЮТСЯ, не копируются):
  * удержание SCF-ветки: реактант cold Newton + internal-stability цикл,
    TS warm-стартом с dm реактанта (scf_stable);
  * фиксированное активное пространство: AVAS(π: C 2p, O 2px, O 2py) на РЕАКТАНТЕ,
    на TS те же (nelec,ncas) принудительно, старт = project_init_guess (run_casscf);
  * state-specific CASSCF c fix_spin ss=0.
Поверх каждой сошедшейся точки (реактант/TS × Cu₂/CuAl) — SC-NEVPT2:
  primary : mrpt.NEVPT(mc).kernel() прямо на CASSCF-объекте (pyscf 2.x,
            fix_spin-fcisolver наследует класс direct_spin1 → NEVPT работает);
  fallback: CASCI с ТОЧНЫМИ интегралами (plain RHF-контейнер, без DF и без
            fix_spin-обёртки) на сошедшихся CASSCF-орбиталях → NEVPT; S² решения
            проверяется — если CASCI ушёл из синглета, число НЕ доверяется.

Чекпойнты casscf-прогона (*.chk.npz) не закоммичены → на инстансе SCF+CASSCF
пересчитываются с нуля; воспроизведение сверяется с закоммиченными числами
co2_to_fuels_casscf_results.json (repro в mHa; |Δ|>TOL_REPRO_MHA → поправка
НЕ привязана к опубликованным числам, reliable=false с причиной).

Отказоустойчивость: каждая тяжёлая стадия — под env-таймаутом
(NEVPT2_STAGE_TIMEOUT_S, default 5400 c; 0 = выключить) и в try/except —
отказ NEVPT2 на одной точке пишется честным fail'ом и НЕ валит джобу.
Incremental atomic-save в routes/co2_to_fuels_nevpt2_results.json после каждой
стадии; собственные чекпойнты co2_to_fuels_nevpt2_<site>.chk.npz (не коммитить)
делают скрипт рестартопригодным.

Итог на сайт: barrier_pbe_eV / barrier_casscf_eV / barrier_nevpt2_eV,
quantum_corr_casscf_eV, quantum_corr_nevpt2_eV, nevpt2_minus_casscf_eV,
флаги честности. Порядок сайтов — в "_summary".

Запуск полный:  OMP_NUM_THREADS=8 python3 -u routes/co2_to_fuels_nevpt2_aws.py
Дым-тест (2×CO в вакууме, малый CAS + NEVPT2, ~2 мин, 1 ядро):
                OMP_NUM_THREADS=1 python3 routes/co2_to_fuels_nevpt2_aws.py --smoke

ЧЕСТНО: замороженный металл-кластер, def2-SVP, DF в SCF/CASSCF (интегралы NEVPT2 —
по реализации pyscf.mrpt), нет constant-potential; state-specific NEVPT2 в двух
точках PBE-скана, путь не пересканирован — геометрии наследуются от PBE.
Это поправка Δ(NEVPT2−PBE) к кластерному барьеру, НЕ фарадеевская эффективность.
"""
import os
import sys
import json
import time
import signal
import argparse
import tempfile
from contextlib import contextmanager

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

from pyscf import scf, mcscf, mrpt, lib                # noqa: E402
from pyscf.mcscf import avas                            # noqa: E402
from pyscf.mcscf import addons as mc_addons             # noqa: E402

import co2_to_fuels_casscf as base                      # noqa: E402
from co2_to_fuels_casscf import (                       # noqa: E402
    HARTREE_EV, SITES, AVAS_LABELS, AVAS_THRESHOLD, TOL_BARRIER_EV,
    METALS, _mol, say, atomic_save, scf_stable, run_casscf, load_site_geoms,
    _load_chk, _smoke_mol)

OUT_RESULTS = os.path.join(HERE, "co2_to_fuels_nevpt2_results.json")
CASSCF_RESULTS = base.OUT_RESULTS      # закоммиченные числа CASSCF-прогона (кросс-чек)
CHK_TPL = os.path.join(HERE, "co2_to_fuels_nevpt2_{site}.chk.npz")

STAGE_TIMEOUT_S = int(os.environ.get("NEVPT2_STAGE_TIMEOUT_S", "5400"))
TOL_REPRO_MHA = 1.0        # |E_CASSCF(здесь) − E_CASSCF(коммит)| — «то же решение»
TOL_SS_SINGLET = 0.01      # допуск S² для fallback-CASCI без fix_spin


# ════════════════════════════════════════════════════════════════════════
#  Env-таймаут на стадию (SIGALRM; только main thread — здесь так и есть)
# ════════════════════════════════════════════════════════════════════════
class StageTimeout(Exception):
    pass


@contextmanager
def stage_timeout(seconds, label):
    if seconds <= 0:
        yield
        return

    def _handler(signum, frame):
        raise StageTimeout(f"{label}: exceeded {seconds}s (NEVPT2_STAGE_TIMEOUT_S)")

    old = signal.signal(signal.SIGALRM, _handler)
    signal.alarm(int(seconds))
    try:
        yield
    finally:
        signal.alarm(0)
        signal.signal(signal.SIGALRM, old)


def _nev_ok(point):
    return bool(point.get("nevpt2", {}).get("ok"))


# ════════════════════════════════════════════════════════════════════════
#  SC-NEVPT2 на одной точке: primary NEVPT(mc), fallback exact-CASCI→NEVPT
# ════════════════════════════════════════════════════════════════════════
def nevpt2_point(mf, mc, tag, timeout_s):
    """Возвращает dict: ok/path/e_ref/e_pt2_corr/e_tot либо честный fail (ok=False,
    error=...). Никогда не бросает — отказ одной точки не валит джобу."""
    out = {}
    t0 = time.time()
    say(f"    [{tag}] NEVPT2 (state-specific, on converged CASSCF) starting …")
    try:
        with stage_timeout(timeout_s, f"{tag}/NEVPT2(casscf)"):
            e_corr = float(mrpt.NEVPT(mc).kernel())
        out.update(ok=True, path="mrpt.NEVPT(mc) on fix_spin CASSCF",
                   e_ref=round(float(mc.e_tot), 6),
                   e_pt2_corr=round(e_corr, 6),
                   e_tot=round(float(mc.e_tot) + e_corr, 6),
                   t_nevpt2_s=round(time.time() - t0, 1))
        say(f"    [{tag}] NEVPT2: E_ref={out['e_ref']:.6f}  "
            f"dE_pt2={out['e_pt2_corr']:.6f}  E_tot={out['e_tot']:.6f}  "
            f"({out['t_nevpt2_s']}s)")
        return out
    except StageTimeout as exc:
        out.update(ok=False, error=f"timeout: {exc}",
                   t_nevpt2_s=round(time.time() - t0, 1))
        say(f"    [{tag}] NEVPT2 TIMED OUT — honest fail recorded, job continues")
        return out                      # бюджет стадии исчерпан — без fallback
    except MemoryError as exc:
        out.update(ok=False, error=f"MemoryError: {exc}",
                   t_nevpt2_s=round(time.time() - t0, 1))
        say(f"    [{tag}] NEVPT2 OOM — honest fail recorded, job continues")
        return out                      # fallback не дешевле по памяти — не пытаемся
    except Exception as exc:
        out["primary_error"] = f"{type(exc).__name__}: {exc}"
        say(f"    [{tag}] NEVPT2 on CASSCF failed ({out['primary_error']}) "
            f"-> exact-integral CASCI fallback")

    # fallback: точные интегралы (plain RHF-контейнер), plain fcisolver, S²-контроль
    t1 = time.time()
    try:
        with stage_timeout(timeout_s, f"{tag}/NEVPT2(casci-fallback)"):
            mf_x = scf.RHF(mf.mol)               # контейнер БЕЗ DF: exact integrals
            mf_x.mo_coeff = mf.mo_coeff
            mf_x.mo_occ = mf.mo_occ
            mf_x.mo_energy = mf.mo_energy
            ci = mcscf.CASCI(mf_x, mc.ncas, mc.nelecas)
            ci.verbose = 0
            ci.kernel(mc.mo_coeff)
            ss = float(ci.fcisolver.spin_square(ci.ci, ci.ncas, ci.nelecas)[0])
            e_corr = float(mrpt.NEVPT(ci).kernel())
        out.update(ok=True,
                   path="mrpt.NEVPT(CASCI exact-integral, plain fcisolver) "
                        "on converged CASSCF orbitals",
                   e_ref=round(float(ci.e_tot), 6),
                   casci_minus_casscf_mHa=round(
                       (float(ci.e_tot) - float(mc.e_tot)) * 1e3, 3),
                   spin_square_casci=round(ss, 4),
                   e_pt2_corr=round(e_corr, 6),
                   e_tot=round(float(ci.e_tot) + e_corr, 6),
                   t_nevpt2_s=round(time.time() - t1, 1))
        if abs(ss) > TOL_SS_SINGLET:
            out["ok"] = False
            out["error"] = (f"fallback CASCI (no fix_spin) left the singlet "
                            f"(S^2={ss:.4f} > {TOL_SS_SINGLET}) — NEVPT2 value "
                            f"NOT trusted")
        say(f"    [{tag}] NEVPT2(fallback): E_ref={out['e_ref']:.6f}  "
            f"dE_pt2={out['e_pt2_corr']:.6f}  S^2={ss:.4f}  ok={out['ok']}  "
            f"({out['t_nevpt2_s']}s)")
    except (StageTimeout, Exception) as exc:
        out.update(ok=False, error=f"{type(exc).__name__}: {exc}",
                   t_nevpt2_s=round(time.time() - t1, 1))
        say(f"    [{tag}] NEVPT2 fallback failed too — honest fail recorded")
    return out


# ════════════════════════════════════════════════════════════════════════
#  Пара реактант/TS на одном сайте (общий путь real- и smoke-режимов)
# ════════════════════════════════════════════════════════════════════════
def nevpt2_pair(tag, mol_r, mol_t, labels, threshold, res, out_path, chk_path,
                profile=None, published=None, args=None):
    """Протокол co2_to_fuels_casscf.compute_pair + NEVPT2 на каждой точке.
    published = запись этого сайта из co2_to_fuels_casscf_results.json (кросс-чек
    воспроизведения) или None. Каждая стадия — atomic_save; рестарт по JSON+chk."""
    entry = res.setdefault(tag, {})
    entry.pop("site_error", None)
    if profile is not None:
        entry["geom"] = {k: profile[k] for k in
                         ("d_reactant", "d_ts", "e_pbe_profile_reactant",
                          "e_pbe_profile_ts", "barrier_pbe_profile_eV")}
    atomic_save(res, out_path)
    chk = _load_chk(chk_path)
    pub = published or {}
    max_macro = getattr(args, "max_macro", 100)
    verbose = getattr(args, "verbose", 0)
    tmo = getattr(args, "stage_timeout", STAGE_TIMEOUT_S)

    r = entry.setdefault("reactant", {})
    t = entry.setdefault("ts", {})

    # ---- реактант: SCF → AVAS(фикс. CAS) → CASSCF → NEVPT2 -------------------
    if not (_nev_ok(r) and "react_cas_mo" in chk):
        say(f"  [{tag}] reactant: SCF (newton + internal stability"
            + (", warm from checkpoint" if "react_dm" in chk else ", cold") + ")")
        with stage_timeout(tmo, f"{tag}/reactant SCF"):
            mf_r, sinfo = scf_stable(mol_r, dm0=chk.get("react_dm"),
                                     tag=f"{tag}/reactant", verbose=verbose)
        if profile is not None:
            sinfo["de_vs_profile_eV"] = round(
                (float(mf_r.e_tot) - profile["e_pbe_profile_reactant"]) * HARTREE_EV, 3)
        if pub.get("reactant", {}).get("e_scf") is not None:
            sinfo["scf_repro_vs_casscf_run_mHa"] = round(
                (float(mf_r.e_tot) - pub["reactant"]["e_scf"]) * 1e3, 3)
        r.update(sinfo)
        atomic_save(res, out_path)

        if "ncas" in chk:
            ncas, nelec = int(chk["ncas"]), int(chk["nelecas"])
            orbs = chk["react_cas_mo"]            # рестарт: с сошедшихся орбиталей
            say(f"  [{tag}] fixed active space from checkpoint: ({nelec}e,{ncas}o)")
        else:
            say(f"  [{tag}] AVAS({labels}, thr={threshold}) on REACTANT defines "
                f"the FIXED active space")
            with stage_timeout(tmo, f"{tag}/AVAS"):
                ncas, nelec, orbs = avas.avas(mf_r, labels, threshold=threshold,
                                              verbose=0)
            ncas, nelec = int(ncas), int(nelec)
            say(f"  [{tag}] fixed active space: ({nelec}e,{ncas}o)")
        entry["cas_space"] = dict(ncas=ncas, nelecas=nelec, avas_labels=list(labels),
                                  avas_threshold=threshold, defined_on="reactant")
        if pub.get("cas_space"):
            entry["cas_space_matches_casscf_run"] = bool(
                (pub["cas_space"]["ncas"], pub["cas_space"]["nelecas"])
                == (ncas, nelec))
        atomic_save(res, out_path)

        with stage_timeout(tmo, f"{tag}/reactant CASSCF"):
            mc_r, cinfo = run_casscf(mf_r, ncas, nelec, orbs, tag=f"{tag}/reactant",
                                     max_macro=max_macro, verbose=verbose)
        if pub.get("reactant", {}).get("e_casscf") is not None:
            cinfo["casscf_repro_vs_casscf_run_mHa"] = round(
                (float(mc_r.e_tot) - pub["reactant"]["e_casscf"]) * 1e3, 3)
        r.update(cinfo)
        chk.update(react_dm=np.asarray(mf_r.make_rdm1()),
                   react_cas_mo=np.asarray(mc_r.mo_coeff),
                   ncas=np.int64(ncas), nelecas=np.int64(nelec))
        np.savez_compressed(chk_path, **chk)
        say(f"  [{tag}] checkpoint -> {os.path.basename(chk_path)}")
        atomic_save(res, out_path)

        if not _nev_ok(r):
            r["nevpt2"] = nevpt2_point(mf_r, mc_r, f"{tag}/reactant", tmo)
            atomic_save(res, out_path)
    else:
        ncas, nelec = int(chk["ncas"]), int(chk["nelecas"])
        say(f"  [{tag}] reactant complete (restart): ({nelec}e,{ncas}o) "
            f"from checkpoint, skipping")

    # ---- TS: ТА ЖЕ SCF-ветка, ТО ЖЕ CAS, project_init_guess → NEVPT2 --------
    if not _nev_ok(t):
        say(f"  [{tag}] TS: SCF warm-started from reactant dm (same branch)")
        with stage_timeout(tmo, f"{tag}/TS SCF"):
            mf_t, sinfo = scf_stable(mol_t, dm0=chk["react_dm"], tag=f"{tag}/ts",
                                     verbose=verbose)
        if profile is not None:
            sinfo["de_vs_profile_eV"] = round(
                (float(mf_t.e_tot) - profile["e_pbe_profile_ts"]) * HARTREE_EV, 3)
        if pub.get("ts", {}).get("e_scf") is not None:
            sinfo["scf_repro_vs_casscf_run_mHa"] = round(
                (float(mf_t.e_tot) - pub["ts"]["e_scf"]) * 1e3, 3)
        t.update(sinfo)
        atomic_save(res, out_path)

        say(f"  [{tag}] TS CASSCF: project_init_guess(reactant CAS orbitals), "
            f"FORCED ({nelec}e,{ncas}o)")
        with stage_timeout(tmo, f"{tag}/TS CASSCF"):
            if "ts_cas_mo" in chk:                # рестарт: с сошедшихся орбиталей
                mo_guess = chk["ts_cas_mo"]
            else:
                mc_t0 = mcscf.CASSCF(mf_t, ncas, nelec)
                mo_guess = mc_addons.project_init_guess(
                    mc_t0, chk["react_cas_mo"], prev_mol=mol_r)
            mc_t, cinfo = run_casscf(mf_t, ncas, nelec, mo_guess, tag=f"{tag}/ts",
                                     max_macro=max_macro, verbose=verbose)
        if pub.get("ts", {}).get("e_casscf") is not None:
            cinfo["casscf_repro_vs_casscf_run_mHa"] = round(
                (float(mc_t.e_tot) - pub["ts"]["e_casscf"]) * 1e3, 3)
        t.update(cinfo)
        chk.update(ts_cas_mo=np.asarray(mc_t.mo_coeff))
        np.savez_compressed(chk_path, **chk)
        atomic_save(res, out_path)

        t["nevpt2"] = nevpt2_point(mf_t, mc_t, f"{tag}/ts", tmo)
        atomic_save(res, out_path)
    else:
        say(f"  [{tag}] TS complete (restart), skipping")

    # ---- барьеры на трёх уровнях, поправки, флаги честности ------------------
    reasons = []
    bp = bc = bn = None
    if r.get("e_scf") is not None and t.get("e_scf") is not None:
        bp = (t["e_scf"] - r["e_scf"]) * HARTREE_EV
        entry["barrier_pbe_eV"] = round(bp, 3)
    if r.get("e_casscf") is not None and t.get("e_casscf") is not None:
        bc = (t["e_casscf"] - r["e_casscf"]) * HARTREE_EV
        entry["barrier_casscf_eV"] = round(bc, 3)
        if bp is not None:
            entry["quantum_corr_casscf_eV"] = round(bc - bp, 3)
    if _nev_ok(r) and _nev_ok(t):
        bn = (t["nevpt2"]["e_tot"] - r["nevpt2"]["e_tot"]) * HARTREE_EV
        entry["barrier_nevpt2_eV"] = round(bn, 3)
        if bp is not None:
            entry["quantum_corr_nevpt2_eV"] = round(bn - bp, 3)
        if bc is not None:
            entry["nevpt2_minus_casscf_eV"] = round(bn - bc, 3)
    else:
        bad = [p for p, d in (("reactant", r), ("ts", t)) if not _nev_ok(d)]
        reasons.append(f"NEVPT2 unavailable at {'/'.join(bad)} "
                       f"(see .nevpt2.error) — no NEVPT2 barrier")

    if not (r.get("scf_converged") and t.get("scf_converged")):
        reasons.append("SCF not converged at reactant and/or TS")
    if not (r.get("internally_stable") and t.get("internally_stable")):
        reasons.append("internal SCF instability persisted")
    if not (r.get("casscf_converged") and t.get("casscf_converged")):
        reasons.append("CASSCF not converged at reactant and/or TS")
    for pname, point in (("reactant", r), ("ts", t)):
        dr = point.get("casscf_repro_vs_casscf_run_mHa")
        if dr is not None and abs(dr) > TOL_REPRO_MHA:
            reasons.append(
                f"{pname} CASSCF here differs from the committed casscf run by "
                f"{dr:+.3f} mHa (tol {TOL_REPRO_MHA}) — different solution; the "
                f"NEVPT2 correction is NOT attached to the published numbers")
    if entry.get("cas_space_matches_casscf_run") is False:
        reasons.append("AVAS-defined active space differs from the committed "
                       "casscf run — different SCF branch or AVAS drift")
    if profile is not None and bp is not None:
        db = bp - profile["barrier_pbe_profile_eV"]
        entry["d_barrier_vs_profile_eV"] = round(db, 3)
        entry["scf_branch_consistent"] = bool(abs(db) <= TOL_BARRIER_EV)
        if not entry["scf_branch_consistent"]:
            reasons.append(
                f"our stable-SCF PBE barrier deviates from the published warm-"
                f"branch barrier by {db:+.3f} eV (tol {TOL_BARRIER_EV}) — "
                f"different SCF solution family")
    entry["reliable"] = not reasons
    entry["unreliable_reasons"] = reasons
    atomic_save(res, out_path)
    say(f"  [{tag}] ΔE‡: PBE={bp if bp is None else round(bp, 3)}  "
        f"CASSCF={bc if bc is None else round(bc, 3)}  "
        f"NEVPT2={bn if bn is None else round(bn, 3)} eV  "
        f"reliable={entry['reliable']}"
        + (f"  reasons={reasons}" if reasons else ""))
    return entry


# ════════════════════════════════════════════════════════════════════════
#  Полный прогон: Cu2 + CuAl; per-site try/except — одна авария не валит всё
# ════════════════════════════════════════════════════════════════════════
def main_real(args):
    say(f"NEVPT2 on the fixed-active-space CASSCF rehabilitation — "
        f"threads={lib.num_threads()}, stage timeout={args.stage_timeout}s")
    res = json.load(open(OUT_RESULTS)) if os.path.exists(OUT_RESULTS) else {}
    res.setdefault("meta", dict(
        method=("SC-NEVPT2 // CASSCF(fixed AVAS space) // "
                f"RKS-{base.XC.upper()}(newton+internal-stability)/{base.BASIS}, "
                "DF in SCF/CASSCF"),
        protocol=("reference protocol imported 1:1 from co2_to_fuels_casscf.py "
                  "(scf_stable, AVAS-on-reactant fixed space, run_casscf with "
                  "fix_spin ss=0, TS via project_init_guess); NEVPT2 primary = "
                  "mrpt.NEVPT(mc) on the converged CASSCF, fallback = exact-"
                  "integral CASCI (plain fcisolver, S^2-checked) on CASSCF "
                  "orbitals; per-stage env timeout + honest per-point fails"),
        reference_run=os.path.basename(CASSCF_RESULTS),
        tol_casscf_repro_mHa=TOL_REPRO_MHA,
        stage_timeout_s=args.stage_timeout,
        pyscf_version=getattr(__import__("pyscf"), "__version__", "?"),
        caveats=("frozen-metal cluster, def2-SVP, DF in SCF/CASSCF, no constant-"
                 "potential; state-specific NEVPT2 at the two PBE-scan points, "
                 "path not rescanned; cluster-barrier correction, NOT Faradaic "
                 "efficiency"),
        started=time.strftime("%Y-%m-%d %H:%M:%S")))
    atomic_save(res, OUT_RESULTS)

    published_all = (json.load(open(CASSCF_RESULTS))
                     if os.path.exists(CASSCF_RESULTS) else {})
    if not published_all:
        say(f"NOTE: {os.path.basename(CASSCF_RESULTS)} not found — no committed "
            f"CASSCF numbers to cross-check against (repro fields skipped)")

    sites = [s.strip() for s in args.sites.split(",") if s.strip()]
    for site in sites:
        if site not in SITES:
            say(f"unknown site {site!r}, skipping")
            continue
        if res.get(site, {}).get("barrier_nevpt2_eV") is not None:
            say(f"--- site {site}: NEVPT2 barrier already present, skipping "
                f"(restart) ---")
            continue
        say(f"--- site {site} ---")
        try:
            g = load_site_geoms(site)
            mol_r = _mol(METALS[site], g["ads_reactant"])
            mol_t = _mol(METALS[site], g["ads_ts"])
            say(f"  [{site}] reactant d={g['d_reactant']}, TS d={g['d_ts']}, "
                f"profile warm barrier {g['barrier_pbe_profile_eV']} eV, "
                f"nao={mol_r.nao}")
            nevpt2_pair(site, mol_r, mol_t, AVAS_LABELS, AVAS_THRESHOLD,
                        res, OUT_RESULTS, CHK_TPL.format(site=site),
                        profile=g, published=published_all.get(site), args=args)
        except (StageTimeout, Exception) as exc:      # honest fail, job continues
            res.setdefault(site, {})["site_error"] = f"{type(exc).__name__}: {exc}"
            res[site]["reliable"] = False
            atomic_save(res, OUT_RESULTS)
            say(f"--- site {site} FAILED ({res[site]['site_error']}) — "
                f"recorded, moving on ---")

    done = [s for s in SITES if res.get(s, {}).get("barrier_nevpt2_eV") is not None]
    if len(done) == len(SITES):
        po = sorted(SITES, key=lambda s: res[s]["barrier_pbe_eV"])
        co = sorted(SITES, key=lambda s: res[s]["barrier_casscf_eV"])
        no_ = sorted(SITES, key=lambda s: res[s]["barrier_nevpt2_eV"])
        res["_summary"] = dict(
            pbe_order=list(po), casscf_order=list(co), nevpt2_order=list(no_),
            order_flipped_vs_pbe=bool(list(po) != list(no_)),
            order_flipped_vs_casscf=bool(list(co) != list(no_)),
            all_reliable=bool(all(res[s].get("reliable") for s in SITES)),
            barrier_pbe_eV={s: res[s]["barrier_pbe_eV"] for s in SITES},
            barrier_casscf_eV={s: res[s]["barrier_casscf_eV"] for s in SITES},
            barrier_nevpt2_eV={s: res[s]["barrier_nevpt2_eV"] for s in SITES},
            quantum_corr_nevpt2_eV={s: res[s].get("quantum_corr_nevpt2_eV")
                                    for s in SITES},
            note=("SC-NEVPT2 (dynamic correlation) on the SAME held SCF branch "
                  "and per-site FIXED active space as the committed CASSCF "
                  "rehabilitation; reproduction vs the committed run recorded "
                  "in mHa per point. Frozen-metal cluster, def2-SVP, DF, no "
                  "constant-potential. If all_reliable=false — see per-site "
                  "unreliable_reasons; numbers then must NOT be quoted as "
                  "physics."))
        atomic_save(res, OUT_RESULTS)
        say(f"=== ORDER: PBE {po}  CASSCF {co}  NEVPT2 {no_}  "
            f"flipped_vs_pbe={res['_summary']['order_flipped_vs_pbe']}  "
            f"all_reliable={res['_summary']['all_reliable']}")
    else:
        say(f"NEVPT2 barrier available for {done or 'no'} site(s) only — "
            f"summary withheld (honest partial result)")
    say(f"saved -> {OUT_RESULTS}; done.")
    return 0


# ════════════════════════════════════════════════════════════════════════
#  Дым-тест: 2×CO в вакууме, тот же протокол + NEVPT2, малый CAS, 1 ядро
# ════════════════════════════════════════════════════════════════════════
def main_smoke(args):
    t0 = time.time()
    tmpd = tempfile.mkdtemp(prefix="co2f_nevpt2_smoke_")
    out = os.path.join(tmpd, "smoke_results.json")
    chk = os.path.join(tmpd, "smoke.chk.npz")
    say(f"SMOKE: 2xCO in vacuum, same protocol + NEVPT2, small CAS; "
        f"threads={lib.num_threads()}; out={out}")
    res = {"meta": dict(mode="smoke",
                        system="2x CO, d(C..C)=4.2 (reactant-like) vs 2.6 "
                               "(TS-like), no metal",
                        avas_labels=["C 2px"], avas_threshold=0.2)}
    entry = nevpt2_pair("smoke_2CO", _smoke_mol(4.2), _smoke_mol(2.6),
                        ["C 2px"], 0.2, res, out, chk,
                        profile=None, published=None, args=args)
    ok = (entry["reliable"]
          and _nev_ok(entry["reactant"]) and _nev_ok(entry["ts"])
          and entry.get("barrier_nevpt2_eV") is not None
          and entry["reactant"]["nevpt2"]["e_tot"] < entry["reactant"]["e_casscf"])
    say("SMOKE: restart check (everything must be skipped)")
    nevpt2_pair("smoke_2CO", _smoke_mol(4.2), _smoke_mol(2.6),
                ["C 2px"], 0.2, res, out, chk,
                profile=None, published=None, args=args)
    say(f"SMOKE {'PASS' if ok else 'FAIL'} in {time.time() - t0:.0f}s  "
        f"CAS=({entry['reactant']['nelecas']}e,{entry['reactant']['ncas']}o)  "
        f"barrier_pbe={entry['barrier_pbe_eV']}  "
        f"barrier_casscf={entry['barrier_casscf_eV']}  "
        f"barrier_nevpt2={entry.get('barrier_nevpt2_eV')} eV  "
        f"nevpt2_path='{entry['reactant']['nevpt2'].get('path')}'")
    return 0 if ok else 1


def _parse_args(argv=None):
    p = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    p.add_argument("--smoke", action="store_true",
                   help="toy 2xCO pipeline test incl. NEVPT2 (~2 min, 1 core)")
    p.add_argument("--sites", default="Cu2,CuAl",
                   help="comma-separated subset of sites (default: Cu2,CuAl)")
    p.add_argument("--max-macro", type=int, default=100, dest="max_macro",
                   help="CASSCF max macro iterations (default 100)")
    p.add_argument("--stage-timeout", type=int, default=STAGE_TIMEOUT_S,
                   dest="stage_timeout",
                   help="per-stage timeout, s; 0 disables "
                        "(default: env NEVPT2_STAGE_TIMEOUT_S or 5400)")
    p.add_argument("--verbose", type=int, default=0,
                   help="pyscf verbose level for SCF/CASSCF (default 0)")
    return p.parse_args(argv)


if __name__ == "__main__":
    _args = _parse_args()
    sys.exit(main_smoke(_args) if _args.smoke else main_real(_args))
