#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
routes/co2_to_fuels_casscf.py — лекарство от UNRELIABLE CASCI-поправки барьера C–C.

Проблема (ts_barrier.casci_note в co2_to_fuels_results.json): per-geometry AVAS давал
разные активные пространства (Cu2 12e↔14e, CuAl 11o↔10o), а у металл-кластера
несколько SCF-решений (cold-restart отличается от релаксационной ветки ~10 эВ) →
barrier_casci / quantum_corr — артефакты, не физика.

Лекарство (этот скрипт):
 1) Удержание SCF-ветки. Реактант: RKS(PBE)/def2-SVP/DF, Newton(SOSCF) + цикл
    ВНУТРЕННЕЙ стабильности до internal-stable. TS: ТА ЖЕ ветка warm-стартом
    (dm реактанта) + снова стабильность. Критерий консистентности ветки: наш
    ΔE‡(PBE) совпадает с warm-барьером опубликованного профиля в пределах
    TOL_BARRIER_EV (иначе — объективный отказ, reliable=false с причиной).
 2) ФИКСИРОВАННОЕ активное пространство на сайт: AVAS(π-многообразие: C 2p,
    O 2px, O 2py) на РЕАКТАНТЕ определяет (nelec, ncas); на TS те же (nelec, ncas)
    принудительно, стартовые орбитали — project_init_guess CASSCF-орбиталей
    реактанта. Никакого повторного AVAS на TS.
 3) CASSCF (orbital-optimized, fix_spin ss=0) вместо CASCI → энергия не зависит
    от произвола AVAS-отбора на каждой геометрии.

Геометрии — реконструкция из ts_barrier.profile результата co2_to_fuels_ts.py
(build_ads: углероды на x=±d/2, free = [zC1,zC2, O1xyz, O2xyz]):
  Cu2:  реактант d=4.2, TS d=2.6  (профильный warm-барьер PBE 1.326 эВ)
  CuAl: реактант d=3.0, TS d=2.2  (профильный warm-барьер PBE 0.992 эВ)

Итог на сайт: barrier_pbe_eV (наша стабильная ветка), barrier_casscf_eV,
quantum_corr_eV = CASSCF−PBE, NOON/n_u обеих точек, флаги достоверности.
Порядок сайтов — в "_summary".

Файлы: co2_to_fuels_casscf_results.json — atomic write после КАЖДОЙ стадии
       (SCF/CASSCF × реактант/TS × сайт), рестартопригодно;
       co2_to_fuels_casscf_<site>.chk.npz — орбитали/dm для рестарта (не коммитить).

Запуск полный:  OMP_NUM_THREADS=4 python3 -u routes/co2_to_fuels_casscf.py
Дым-тест (плаумбинг на 2×CO в вакууме, малый CAS, ~1 мин, 1 ядро):
                OMP_NUM_THREADS=1 python3 routes/co2_to_fuels_casscf.py --smoke

Ограничения: замороженный металл-кластер, def2-SVP, RI/DF везде, нет constant-potential.
Это поправка Δ(CASSCF−PBE) к кластерному барьеру сочленения, НЕ фарадеевская
эффективность. CASSCF может увести активные орбитали из π-многообразия в металл —
это вариационно легально, но метку «π» тогда читать как «стартовавшее из π».
"""
import os
import sys
import json
import time
import argparse
import tempfile

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

from pyscf import gto, dft, mcscf, fci, lib          # noqa: E402
from pyscf.mcscf import avas                          # noqa: E402
from pyscf.mcscf import addons as mc_addons           # noqa: E402

HARTREE_EV = 27.211386245988
TS_RESULTS = os.path.join(HERE, "co2_to_fuels_results.json")
OUT_RESULTS = os.path.join(HERE, "co2_to_fuels_casscf_results.json")
CHK_TPL = os.path.join(HERE, "co2_to_fuels_casscf_{site}.chk.npz")

SITES = ("Cu2", "CuAl")
AVAS_LABELS = ["C 2p", "O 2px", "O 2py"]   # π-многообразие — как в co2_to_fuels_ts.py
AVAS_THRESHOLD = 0.4
XC = "pbe"
BASIS = "def2-svp"
MAX_STAB_LOOPS = 4          # rotate+reconverge попыток до признания провала
TOL_BARRIER_EV = 0.05       # |ΔE‡(наша ветка) − ΔE‡(профиль warm)| — критерий той же ветки
CONV_TOL_SCF = 1e-8
CONV_TOL_CAS = 1e-7

# ── геометрия: те же определения, что в co2_to_fuels_sites/_ts (импорт → 1:1) ──
try:
    from co2_to_fuels_sites import METALS, _mol                    # noqa: E402
    from co2_to_fuels_ts import build_ads                          # noqa: E402
except ImportError:                       # напр. нет scipy: локальные ТОЧНЫЕ копии
    METALS = {
        "Cu2":  [("Cu", (-1.225, 0., 0.)), ("Cu", (1.225, 0., 0.))],
        "CuAl": [("Cu", (-1.30, 0., 0.)), ("Al", (1.30, 0., 0.))],
    }

    def _atoms_str(metal, ads):
        s = ""
        for el, (x, y, z) in metal:
            s += f"{el} {x} {y} {z}; "
        for el, (x, y, z) in zip(["C", "O", "C", "O"], ads):
            s += f"{el} {x} {y} {z}; "
        return s

    def _mol(metal, ads):
        return gto.M(atom=_atoms_str(metal, ads), basis=BASIS, spin=0, verbose=0)

    def build_ads(d, free):
        """free = [zC1, zC2, O1x,O1y,O1z, O2x,O2y,O2z]; C на x=±d/2, y=0."""
        zC1, zC2, o1x, o1y, o1z, o2x, o2y, o2z = free
        return np.array([[-d / 2, 0, zC1], [o1x, o1y, o1z],
                         [d / 2, 0, zC2], [o2x, o2y, o2z]])


def say(msg):
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def atomic_save(obj, path):
    """Атомарная запись JSON: tmp в той же директории + os.replace."""
    fd, tmp = tempfile.mkstemp(dir=os.path.dirname(path) or ".",
                               prefix=os.path.basename(path) + ".", suffix=".tmp")
    try:
        with os.fdopen(fd, "w") as f:
            json.dump(obj, f, ensure_ascii=False, indent=1)
            f.flush()
            os.fsync(f.fileno())
        os.replace(tmp, path)
    except BaseException:
        if os.path.exists(tmp):
            os.unlink(tmp)
        raise


# ════════════════════════════════════════════════════════════════════════
#  SCF: Newton + внутренняя стабильность (удержание одной ветки)
# ════════════════════════════════════════════════════════════════════════
def _adopt(mf, mfn):
    mf.mo_coeff, mf.mo_energy, mf.mo_occ = mfn.mo_coeff, mfn.mo_energy, mfn.mo_occ
    mf.e_tot, mf.converged = float(mfn.e_tot), bool(mfn.converged)


def _newton_solve(mf, dm0=None, mo0=None, occ0=None):
    mfn = mf.newton()
    mfn.verbose = mf.verbose
    mfn.max_cycle = 80
    if mo0 is not None:
        mfn.kernel(mo0, occ0)
    else:
        try:
            mfn.kernel(dm0=dm0)
        except TypeError:                 # старый pyscf без dm0 в newton.kernel
            mfn.kernel()
    _adopt(mf, mfn)
    return mf


def _internal_stability(mf):
    """→ (mo_i, stable_i); терпимо к разным версиям API pyscf."""
    try:
        out = mf.stability(internal=True, external=False, return_status=True)
    except TypeError:                                     # нет return_status
        out = mf.stability(internal=True, external=False)
        mo_i = out[0] if isinstance(out, (tuple, list)) else out
        stable = mo_i is mf.mo_coeff or (
            np.asarray(mo_i).shape == np.asarray(mf.mo_coeff).shape
            and np.allclose(mo_i, mf.mo_coeff))
        return mo_i, bool(stable)
    if isinstance(out, (tuple, list)) and len(out) == 4:  # (mo_i, mo_e, st_i, st_e)
        return out[0], bool(out[2])
    if isinstance(out, (tuple, list)) and len(out) == 2:
        return out[0], bool(out[1])
    raise RuntimeError(f"unexpected stability() return: {type(out)}")


def _gap_ev(mf):
    e = np.asarray(mf.mo_energy)
    occ = np.asarray(mf.mo_occ)
    return float((e[occ == 0].min() - e[occ > 0].max()) * HARTREE_EV)


def scf_stable(mol, dm0=None, tag="", verbose=0):
    """RKS(PBE)/DF через Newton, затем цикл внутренней стабильности до
    internal-stable (или MAX_STAB_LOOPS). dm0 = warm-старт (удержание ветки)."""
    t0 = time.time()
    mf = dft.RKS(mol).density_fit()
    mf.xc = XC
    mf.verbose = verbose
    mf.conv_tol = CONV_TOL_SCF
    _newton_solve(mf, dm0=dm0)
    n_rot, stable = 0, False
    for _ in range(MAX_STAB_LOOPS + 1):
        mo_i, stable = _internal_stability(mf)
        if stable:
            break
        if n_rot >= MAX_STAB_LOOPS:
            break
        n_rot += 1
        say(f"    [{tag}] SCF internally UNSTABLE -> rotate + re-converge (#{n_rot})")
        _newton_solve(mf, mo0=mo_i, occ0=mf.mo_occ)
    info = dict(e_scf=round(float(mf.e_tot), 6),
                scf_converged=bool(mf.converged),
                internally_stable=bool(stable),
                stability_rotations=int(n_rot),
                gap_eV=round(_gap_ev(mf), 3),
                warm_start=bool(dm0 is not None),
                t_scf_s=round(time.time() - t0, 1))
    say(f"    [{tag}] E_SCF={mf.e_tot:.6f}  conv={info['scf_converged']} "
        f"stable={info['internally_stable']} rot={n_rot} gap={info['gap_eV']:.2f} eV "
        f"({info['t_scf_s']}s)")
    return mf, info


# ════════════════════════════════════════════════════════════════════════
#  CASSCF с фиксированным (nelec, ncas)
# ════════════════════════════════════════════════════════════════════════
def run_casscf(mf, ncas, nelecas, mo_guess, tag="", max_macro=100, verbose=0):
    say(f"    [{tag}] CASSCF({int(nelecas)}e,{int(ncas)}o) starting …")
    t0 = time.time()
    mc = mcscf.CASSCF(mf, int(ncas), int(nelecas))
    mc.verbose = verbose
    mc.max_cycle_macro = max_macro
    mc.conv_tol = CONV_TOL_CAS
    try:
        fci.addons.fix_spin_(mc.fcisolver, ss=0)          # держим синглет
    except Exception as exc:                              # строго говоря фиксируем, если нет
        say(f"    [{tag}] fix_spin_ unavailable: {exc!r}")
    mc.kernel(mo_guess)
    dm1 = mc.fcisolver.make_rdm1(mc.ci, mc.ncas, mc.nelecas)
    noon = np.linalg.eigvalsh(dm1)[::-1]
    n_u = float(np.sum(np.minimum(noon, 2.0 - noon)))
    try:
        ss = float(mc.fcisolver.spin_square(mc.ci, mc.ncas, mc.nelecas)[0])
    except Exception:
        ss = None
    ne_out = int(sum(mc.nelecas)) if np.ndim(mc.nelecas) else int(mc.nelecas)
    info = dict(e_casscf=round(float(mc.e_tot), 6),
                casscf_converged=bool(mc.converged),
                ncas=int(mc.ncas), nelecas=ne_out,
                noon=[round(float(x), 4) for x in noon],
                n_u=round(n_u, 4),
                spin_square=(round(ss, 4) if ss is not None else None),
                t_casscf_s=round(time.time() - t0, 1))
    say(f"    [{tag}] E_CASSCF={mc.e_tot:.6f}  conv={info['casscf_converged']} "
        f"n_u={n_u:.3f} S^2={info['spin_square']} ({info['t_casscf_s']}s)")
    return mc, info


# ════════════════════════════════════════════════════════════════════════
#  Пара реактант/TS на одном сайте (общий путь real- и smoke-режимов)
# ════════════════════════════════════════════════════════════════════════
def _load_chk(path):
    if not os.path.exists(path):
        return {}
    with np.load(path) as z:
        return {k: z[k] for k in z.files}


def compute_pair(tag, mol_r, mol_t, labels, threshold, res, out_path, chk_path,
                 profile=None, args=None):
    """Полный протокол на сайт. profile = dict с профильными PBE-числами (или None
    в smoke). Каждая стадия — atomic_save; рестарт по JSON + chk.npz."""
    entry = res.setdefault(tag, {})
    if profile is not None:
        entry["geom"] = {k: profile[k] for k in
                         ("d_reactant", "d_ts", "e_pbe_profile_reactant",
                          "e_pbe_profile_ts", "barrier_pbe_profile_eV")}
    atomic_save(res, out_path)
    chk = _load_chk(chk_path)
    max_macro = getattr(args, "max_macro", 100)
    verbose = getattr(args, "verbose", 0)

    # ---- реактант: SCF (cold, newton+stability) + AVAS→фикс. CAS + CASSCF ----
    r = entry.setdefault("reactant", {})
    need_r = (r.get("e_casscf") is None) or ("react_cas_mo" not in chk)
    if need_r:
        say(f"  [{tag}] reactant: SCF (newton + internal stability)")
        mf_r, sinfo = scf_stable(mol_r, dm0=None, tag=f"{tag}/reactant",
                                 verbose=verbose)
        if profile is not None:
            sinfo["de_vs_profile_eV"] = round(
                (float(mf_r.e_tot) - profile["e_pbe_profile_reactant"]) * HARTREE_EV, 3)
        r.update(sinfo)
        atomic_save(res, out_path)

        say(f"  [{tag}] AVAS({labels}, thr={threshold}) on REACTANT defines "
            f"the FIXED active space")
        ncas, nelec, orbs = avas.avas(mf_r, labels, threshold=threshold, verbose=0)
        ncas, nelec = int(ncas), int(nelec)
        entry["cas_space"] = dict(ncas=ncas, nelecas=nelec,
                                  avas_labels=labels, avas_threshold=threshold,
                                  defined_on="reactant")
        say(f"  [{tag}] fixed active space: ({nelec}e,{ncas}o)")
        atomic_save(res, out_path)

        mc_r, cinfo = run_casscf(mf_r, ncas, nelec, orbs, tag=f"{tag}/reactant",
                                 max_macro=max_macro, verbose=verbose)
        r.update(cinfo)
        atomic_save(res, out_path)

        chk = dict(react_mo=np.asarray(mf_r.mo_coeff),
                   react_occ=np.asarray(mf_r.mo_occ),
                   react_dm=np.asarray(mf_r.make_rdm1()),
                   react_cas_mo=np.asarray(mc_r.mo_coeff),
                   ncas=np.int64(ncas), nelecas=np.int64(nelec))
        np.savez_compressed(chk_path, **chk)
        say(f"  [{tag}] checkpoint -> {os.path.basename(chk_path)}")
    else:
        ncas = int(chk["ncas"])
        nelec = int(chk["nelecas"])
        say(f"  [{tag}] reactant already done (restart): "
            f"fixed space ({nelec}e,{ncas}o) from checkpoint")

    # ---- TS: ТА ЖЕ SCF-ветка (warm dm реактанта) + CASSCF в ТОМ ЖЕ CAS ------
    t = entry.setdefault("ts", {})
    if t.get("e_casscf") is None:
        say(f"  [{tag}] TS: SCF warm-started from reactant dm (same branch)")
        mf_t, sinfo = scf_stable(mol_t, dm0=chk["react_dm"], tag=f"{tag}/ts",
                                 verbose=verbose)
        if profile is not None:
            sinfo["de_vs_profile_eV"] = round(
                (float(mf_t.e_tot) - profile["e_pbe_profile_ts"]) * HARTREE_EV, 3)
        t.update(sinfo)
        atomic_save(res, out_path)

        if getattr(args, "cross_check_cold", False):
            say(f"  [{tag}] TS cross-check: independent COLD SCF (multi-solution probe)")
            _, coldinfo = scf_stable(mol_t, dm0=None, tag=f"{tag}/ts-cold",
                                     verbose=verbose)
            t["cold_e_scf"] = coldinfo["e_scf"]
            t["cold_minus_warm_eV"] = round(
                (coldinfo["e_scf"] - t["e_scf"]) * HARTREE_EV, 3)
            atomic_save(res, out_path)

        say(f"  [{tag}] TS CASSCF: project_init_guess(reactant CAS orbitals), "
            f"FORCED ({nelec}e,{ncas}o)")
        mc_t = mcscf.CASSCF(mf_t, ncas, nelec)
        mo_guess = mc_addons.project_init_guess(mc_t, chk["react_cas_mo"],
                                                prev_mol=mol_r)
        mc_t, cinfo = run_casscf(mf_t, ncas, nelec, mo_guess, tag=f"{tag}/ts",
                                 max_macro=max_macro, verbose=verbose)
        t.update(cinfo)
        atomic_save(res, out_path)
    else:
        say(f"  [{tag}] TS already done (restart), skipping")

    # ---- барьеры, поправка, флаги достоверности ---------------------------------
    bp = (t["e_scf"] - r["e_scf"]) * HARTREE_EV
    bc = (t["e_casscf"] - r["e_casscf"]) * HARTREE_EV
    entry["barrier_pbe_eV"] = round(bp, 3)
    entry["barrier_casscf_eV"] = round(bc, 3)
    entry["quantum_corr_eV"] = round(bc - bp, 3)

    reasons = []
    if not (r["scf_converged"] and t["scf_converged"]):
        reasons.append("SCF not converged at reactant and/or TS")
    if not (r["internally_stable"] and t["internally_stable"]):
        reasons.append(f"internal SCF instability persisted after "
                       f"{MAX_STAB_LOOPS} rotations")
    if not (r["casscf_converged"] and t["casscf_converged"]):
        reasons.append("CASSCF not converged at reactant and/or TS")
    if (r["ncas"], r["nelecas"]) != (t["ncas"], t["nelecas"]):
        reasons.append("active space differs between reactant and TS "
                       "(must never happen — enforced)")
    if profile is not None:
        db = bp - profile["barrier_pbe_profile_eV"]
        entry["d_barrier_vs_profile_eV"] = round(db, 3)
        entry["scf_branch_consistent"] = bool(abs(db) <= TOL_BARRIER_EV)
        if not entry["scf_branch_consistent"]:
            reasons.append(
                f"our stable-SCF PBE barrier deviates from the published "
                f"warm-branch barrier by {db:+.3f} eV (tol {TOL_BARRIER_EV}) — "
                f"different SCF solution family; Δ(CASSCF−PBE) is NOT "
                f"transferable to the published PBE barrier")
    entry["reliable"] = not reasons
    entry["unreliable_reasons"] = reasons
    atomic_save(res, out_path)
    say(f"  [{tag}] ΔE‡: PBE={bp:+.3f}  CASSCF={bc:+.3f}  "
        f"Δ_corr={bc - bp:+.3f} eV  reliable={entry['reliable']}"
        + (f"  reasons={reasons}" if reasons else ""))
    return entry


# ════════════════════════════════════════════════════════════════════════
#  Полный прогон: Cu2 + CuAl из ts_barrier.profile
# ════════════════════════════════════════════════════════════════════════
def load_site_geoms(site):
    with open(TS_RESULTS) as f:
        tb = json.load(f)["ts_barrier"][site]
    prof = {round(p["d"], 2): p for p in tb["profile"]}
    d_r, d_t = float(tb["d_reactant"]), float(tb["d_ts"])
    pr, pt = prof[round(d_r, 2)], prof[round(d_t, 2)]
    prof_barrier = round((pt["e_pbe"] - pr["e_pbe"]) * HARTREE_EV, 3)
    stored = tb.get("barrier_pbe_warm_eV", tb.get("barrier_pbe_eV"))
    if stored is not None and abs(prof_barrier - stored) > 0.005:
        say(f"  [{site}] NOTE: profile-derived barrier {prof_barrier} eV != stored "
            f"warm barrier {stored} eV — using profile-derived value")
    return dict(d_reactant=d_r, d_ts=d_t,
                e_pbe_profile_reactant=pr["e_pbe"],
                e_pbe_profile_ts=pt["e_pbe"],
                barrier_pbe_profile_eV=prof_barrier,
                ads_reactant=build_ads(d_r, np.array(pr["free"], dtype=float)),
                ads_ts=build_ads(d_t, np.array(pt["free"], dtype=float)))


def main_real(args):
    say(f"CASSCF fixed-active-space barrier correction — threads={lib.num_threads()}")
    res = json.load(open(OUT_RESULTS)) if os.path.exists(OUT_RESULTS) else {}
    res.setdefault("meta", dict(
        method=(f"CASSCF(fixed AVAS space)//RKS-{XC.upper()}"
                f"(newton+internal-stability)/{BASIS}, density fitting throughout"),
        protocol=("reactant: cold SCF -> internal-stability loop; TS: warm start "
                  "from reactant dm -> stability again; (nelec,ncas) fixed by AVAS "
                  "on the REACTANT; TS CASSCF init = project_init_guess(reactant "
                  "CASSCF orbitals); fcisolver fix_spin ss=0"),
        avas_labels=AVAS_LABELS, avas_threshold=AVAS_THRESHOLD,
        geometry_source=("co2_to_fuels_results.json ts_barrier.profile, "
                         "build_ads reconstruction (co2_to_fuels_ts.py)"),
        tol_barrier_vs_profile_eV=TOL_BARRIER_EV,
        caveats=("frozen-metal cluster, def2-SVP, RI/DF, no constant-potential; "
                 "cluster-barrier correction, NOT Faradaic efficiency; CASSCF may "
                 "rotate active orbitals out of the initial pi manifold"),
        started=time.strftime("%Y-%m-%d %H:%M:%S")))
    atomic_save(res, OUT_RESULTS)

    sites = [s.strip() for s in args.sites.split(",") if s.strip()]
    for site in sites:
        if site not in SITES:
            say(f"unknown site {site!r}, skipping")
            continue
        if res.get(site, {}).get("quantum_corr_eV") is not None and \
                res[site].get("reliable") is not None:
            say(f"--- site {site}: already complete, skipping (restart) ---")
            continue
        say(f"--- site {site} ---")
        g = load_site_geoms(site)
        metal = METALS[site]
        mol_r = _mol(metal, g["ads_reactant"])
        mol_t = _mol(metal, g["ads_ts"])
        say(f"  [{site}] reactant d={g['d_reactant']}, TS d={g['d_ts']}, "
            f"profile warm barrier {g['barrier_pbe_profile_eV']} eV")
        compute_pair(site, mol_r, mol_t, AVAS_LABELS, AVAS_THRESHOLD,
                     res, OUT_RESULTS, CHK_TPL.format(site=site),
                     profile=g, args=args)

    done = [s for s in SITES if res.get(s, {}).get("quantum_corr_eV") is not None]
    if len(done) == len(SITES):
        po = sorted(SITES, key=lambda s: res[s]["barrier_pbe_eV"])
        co = sorted(SITES, key=lambda s: res[s]["barrier_casscf_eV"])
        res["_summary"] = dict(
            pbe_order=list(po), casscf_order=list(co),
            order_flipped=bool(list(po) != list(co)),
            all_reliable=bool(all(res[s]["reliable"] for s in SITES)),
            barrier_pbe_eV={s: res[s]["barrier_pbe_eV"] for s in SITES},
            barrier_casscf_eV={s: res[s]["barrier_casscf_eV"] for s in SITES},
            quantum_corr_eV={s: res[s]["quantum_corr_eV"] for s in SITES},
            note=("Δ(CASSCF−PBE) to the coupling barrier on a single held SCF "
                  "branch with a per-site FIXED active space (the casci_note "
                  "remedy). Frozen-metal cluster, def2-SVP, DF, no "
                  "constant-potential. If all_reliable=false — see per-site "
                  "unreliable_reasons; numbers then must NOT be quoted as "
                  "physics."))
        atomic_save(res, OUT_RESULTS)
        say(f"=== ORDER: PBE {po}  CASSCF {co}  "
            f"flipped={res['_summary']['order_flipped']}  "
            f"all_reliable={res['_summary']['all_reliable']}")
    say(f"saved -> {OUT_RESULTS}; done.")


# ════════════════════════════════════════════════════════════════════════
#  Дым-тест: 2×CO в вакууме, тот же протокол, малый CAS (~4e,4o), ≤3 мин / 1 ядро
# ════════════════════════════════════════════════════════════════════════
def _smoke_mol(d):
    """Два CO, оси молекул вдоль z, разведение по x — как адсорбат без металла."""
    a = (f"C {-d / 2:.4f} 0.0 0.0; O {-d / 2:.4f} 0.0 1.1350; "
         f"C {d / 2:.4f} 0.0 0.0; O {d / 2:.4f} 0.0 1.1350")
    return gto.M(atom=a, basis=BASIS, spin=0, verbose=0)


def main_smoke(args):
    t0 = time.time()
    tmpd = tempfile.mkdtemp(prefix="co2f_casscf_smoke_")
    out = os.path.join(tmpd, "smoke_results.json")
    chk = os.path.join(tmpd, "smoke.chk.npz")
    say(f"SMOKE: 2xCO in vacuum, same protocol, small CAS; threads="
        f"{lib.num_threads()}; out={out}")
    res = {"meta": dict(mode="smoke", system="2x CO, d(C..C)=4.2 (reactant-like) "
                                            "vs 2.6 (TS-like), no metal",
                        avas_labels=["C 2px"], avas_threshold=0.2)}
    # d как у Cu2: 'реактант' 4.2, 'TS' 2.6; AVAS[C 2px] → сигма-фронтир ~ (4e,4o)
    entry = compute_pair("smoke_2CO", _smoke_mol(4.2), _smoke_mol(2.6),
                         ["C 2px"], 0.2, res, out, chk, profile=None, args=args)
    ok = (entry["reliable"]
          and entry["reactant"]["ncas"] == entry["ts"]["ncas"]
          and entry["reactant"]["nelecas"] == entry["ts"]["nelecas"]
          and entry["reactant"]["e_casscf"] < entry["reactant"]["e_scf"] + 1.0)
    # рестарт-проверка: повторный вызов должен всё пропустить по chk+JSON
    say("SMOKE: restart check (everything must be skipped)")
    compute_pair("smoke_2CO", _smoke_mol(4.2), _smoke_mol(2.6),
                 ["C 2px"], 0.2, res, out, chk, profile=None, args=args)
    say(f"SMOKE {'PASS' if ok else 'FAIL'} in {time.time() - t0:.0f}s  "
        f"CAS=({entry['reactant']['nelecas']}e,{entry['reactant']['ncas']}o)  "
        f"barrier_pbe={entry['barrier_pbe_eV']} eV  "
        f"barrier_casscf={entry['barrier_casscf_eV']} eV  "
        f"corr={entry['quantum_corr_eV']} eV")
    return 0 if ok else 1


def _parse_args(argv=None):
    p = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    p.add_argument("--smoke", action="store_true",
                   help="toy 2xCO pipeline test (~1 min, 1 core), no real run")
    p.add_argument("--sites", default="Cu2,CuAl",
                   help="comma-separated subset of sites (default: Cu2,CuAl)")
    p.add_argument("--cross-check-cold", action="store_true", dest="cross_check_cold",
                   help="also run an independent cold SCF at each TS and record "
                        "cold-minus-warm gap (multi-solution probe; extra cost)")
    p.add_argument("--max-macro", type=int, default=100, dest="max_macro",
                   help="CASSCF max macro iterations (default 100)")
    p.add_argument("--verbose", type=int, default=0,
                   help="pyscf verbose level for SCF/CASSCF (default 0)")
    return p.parse_args(argv)


if __name__ == "__main__":
    _args = _parse_args()
    sys.exit(main_smoke(_args) if _args.smoke else main_real(_args))
