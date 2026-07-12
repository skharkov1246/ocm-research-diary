#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Этап 21-F: полевое управление селективностью — ΔΔE‡(F) на NEVPT2.

Идея (изобретение I5, слияние флагмана Cr-допанта и направления I4):
внешнее электрическое поле вдоль оси M=O как РУЧКА анти-BEP-дескриптора.
Если dΔΔE‡/dF значим на коррелированном уровне — селективность OCM-центра
управляется приложенным напряжением (двойной слой/поляризация носителя),
а не только химией допанта. Никем не публиковалось для NEVPT2-дескрипторов.

Постановка (максимальное переиспользование, как TZVP-чек Этапа 20):
- геометрии НЕ пересчитываются — берутся 9 сохранённых точек Cr-профилей
  (routes/ocm_mnw_emb_cr_{ch4,c2h6}_geom.json, frozen-frame SVP);
- поле F = ±0.005 а.е. (±0.257 В/Å — масштаб двойного слоя) вдоль единичного
  вектора Cr(0)→O(1) КАЖДОЙ точки (локальная ось сайта);
- finite-field строго: get_hcore += F·r (глобальный патч класса — поле видят
  ROHF, CASSCF и NEVPT2), к энергии добавляется ядерный член −Σ Z_A F·R_A;
  кластер нейтрален ⇒ полная энергия не зависит от начала координат;
- F=0 контроль: R- и TS-точки обязаны воспроизвести закоммиченные
  e_nevpt2_h профиля (валидация полевого кода бесплатно);
- считаем NEVPT2 не на всех 9 точках, а в окне TS±2 (+ R-точка idx 0 как
  референс) вокруг робастного нуль-полевого TS (маска >25 ккал, интерьерность);
- вердикт per-field: интерьерный максимум ВНУТРИ окна; максимум на краю
  окна ⇒ NO-VERDICT для этого поля (мираж-детектор обязателен и тут).

Стадии: field (скан, поточечный чекпойнт ocm_field_scan.json) | readout.
Запуск: python3 -u routes/ocm_field_dde.py field|readout
"""
import json
import os
import sys
import time

import numpy as np

DIR = os.path.dirname(os.path.abspath(__file__))
# env ДО импорта пайплайна: Cr embedded, основной спин 2S=4 (ladder Этапа 17)
os.environ.setdefault("OCM_MODEL", "embedded")
os.environ.setdefault("OCM_METAL", "Cr")
os.environ.setdefault("OCM_SPIN", "4")
os.environ.setdefault("OCM_MAXMEM", "24000")
sys.path.insert(0, DIR)
import ocm_mnw_selectivity as sel                      # noqa: E402
from pyscf import scf                                  # noqa: E402

HK = sel.HARTREE_KCAL
FIELDS = [float(x) for x in
          os.environ.get("OCM_FIELDS", "-0.005,0.005").split(",")]
AU_FIELD_V_A = 51.42206747                             # 1 а.е. поля в В/Å
OUTLIER_KCAL = 25.0
SCAN = os.path.join(DIR, "ocm_field_scan.json")
RES = os.path.join(DIR, "ocm_field_dde_results.json")

# ---------------------------------------------------------- finite-field патч
FIELD = np.zeros(3)                                    # а.е., лабораторный вектор
_orig_hcore = scf.hf.SCF.get_hcore


def _hcore_field(self, mol=None):
    if mol is None:
        mol = self.mol
    h = _orig_hcore(self, mol)
    if np.any(FIELD):
        with mol.with_common_orig((0, 0, 0)):
            r = mol.intor_symmetric("int1e_r", comp=3)
        h = h + np.einsum("x,xij->ij", FIELD, r)
    return h


scf.hf.SCF.get_hcore = _hcore_field                    # видят UKS/ROHF/CASSCF/PT2


def nuc_field(atoms):
    """Ядерный член −Σ Z_A F·R_A (Hartree); atoms в Å."""
    if not np.any(FIELD):
        return 0.0
    mol = sel.build_mol(atoms)
    return float(-np.einsum("i,ix,x->", mol.atom_charges(),
                            mol.atom_coords(), FIELD))


def site_axis(atoms):
    """Единичный вектор Cr(0)→O(1) — локальная ось M=O данной точки."""
    a = np.array(atoms[0][1], float)
    b = np.array(atoms[1][1], float)
    u = b - a
    return u / np.linalg.norm(u)


# ------------------------------------------------------------------- окно TS
def robust_window(sub):
    """Нуль-полевой робастный TS из закоммиченного профиля: маска >25 ккал
    (выше ОБОИХ соседей), интерьерный максимум. Окно TS±2 + R(0)."""
    prof = json.load(open(os.path.join(DIR, f"ocm_mnw_emb_cr_{sub}_profile.json")))
    rel = prof["nevpt2_rel_kcal"]
    idx = list(range(len(rel)))
    keep = []
    for i in idx:
        nb = [rel[j] for j in (i - 1, i + 1) if 0 <= j < len(rel)]
        if nb and all(rel[i] - x > OUTLIER_KCAL for x in nb):
            continue
        keep.append(i)
    interior = keep[1:-1]
    its = max(interior, key=lambda i: rel[i])
    kpos = keep.index(its)
    win = sorted(set([keep[0]] + keep[max(0, kpos - 2):kpos + 3]))
    return win, its, keep


# --------------------------------------------------------------- стадия field
def stage_field():
    pts = []
    if os.path.exists(SCAN):
        try:
            pts = json.load(open(SCAN))["pts"]
            print(f"[field] resume: {len(pts)} точек с диска", flush=True)
        except Exception:
            pts = []
    done = {(p["sub"], round(p["field"], 4), p["idx"]) for p in pts}

    def save():
        json.dump({"pts": pts}, open(SCAN, "w"), indent=1)

    jobs = []
    for sub in ("ch4", "c2h6"):
        win, its, _ = robust_window(sub)
        print(f"[field] {sub}: окно {win} (нуль-полевой TS idx {its})", flush=True)
        for f in FIELDS:
            for i in win:
                jobs.append((sub, f, i))
        # F=0 контроль: R и TS (сверка с закоммиченным профилем)
        jobs += [(sub, 0.0, win[0]), (sub, 0.0, its)]
    for sub, f, i in jobs:
        if (sub, round(f, 4), i) in done:
            continue
        g = json.load(open(os.path.join(DIR, f"ocm_mnw_emb_cr_{sub}_geom.json")))
        atoms = [(s, tuple(c)) for s, c in g["points"][i]["atoms"]]
        global FIELD
        FIELD = f * site_axis(atoms)
        t0 = time.time()
        nf = nuc_field(atoms)
        mf = sel.uks(sel.build_mol(atoms))
        e_dft = float(mf.e_tot) + nf
        res = sel.cas_nevpt2(atoms, f"{sub}_F{f:+.3f}_{i}")
        for k in ("e_rohf_h", "e_casci_h", "e_casscf_h", "e_nevpt2_h"):
            res[k] = res[k] + nf
        pts.append({"sub": sub, "field": f, "idx": i, "e_dft_h": e_dft,
                    "nuc_field_h": nf, "wall_s": round(time.time() - t0, 1),
                    **{k: res[k] for k in
                       ("e_rohf_h", "e_casscf_h", "e_nevpt2_h",
                        "casscf_converged", "noon",
                        "n_unpaired_beyond_spin")}})
        save()
        print(f"[field] {sub} F={f:+.3f} idx={i}: DFT {e_dft:.6f} "
              f"NEVPT2 {res['e_nevpt2_h']:.6f} ({pts[-1]['wall_s']}s)",
              flush=True)
    print("[field] скан завершён", flush=True)


# ------------------------------------------------------------- стадия readout
def _barrier(series):
    """(барьер, интерьерность-в-окне) по списку (idx, rel_kcal)."""
    if len(series) < 3:
        return None, False
    vals = [v for _, v in series]
    imax = int(np.argmax(vals[1:])) + 1        # исключая R-референс
    interior = 0 < imax < len(vals) - 1
    return vals[imax], interior


def stage_readout():
    pts = json.load(open(SCAN))["pts"]
    out = {"model": "embedded Cr, frozen-frame SVP geometries (Этап 17/18)",
           "axis": "unit(Cr0->O1) per point", "fields_au": FIELDS,
           "field_V_per_A": [round(f * AU_FIELD_V_A, 3) for f in FIELDS],
           "protocol": "finite-field get_hcore+F·r + nuclear term; "
                       "NEVPT2 window TS±2; interior-in-window rule",
           "zero_field_control": {}, "per_field": {}, "ddE": {}}
    # контроль F=0 против закоммиченного профиля
    for sub in ("ch4", "c2h6"):
        prof = json.load(open(os.path.join(DIR,
                              f"ocm_mnw_emb_cr_{sub}_profile.json")))
        ctrl = {}
        for p in pts:
            if p["sub"] == sub and p["field"] == 0.0:
                ref = prof["points"][p["idx"]]["e_nevpt2_h"]
                ctrl[str(p["idx"])] = {
                    "e_nevpt2_h": p["e_nevpt2_h"], "committed": ref,
                    "diff_kcal": round((p["e_nevpt2_h"] - ref) * HK, 3)}
        out["zero_field_control"][sub] = ctrl
    # барьеры per (field): rel к R-точке окна, интерьерность в окне
    dde = {}
    for f in FIELDS:
        row = {}
        for sub in ("ch4", "c2h6"):
            sp = sorted([p for p in pts
                         if p["sub"] == sub and p["field"] == f],
                        key=lambda p: p["idx"])
            if not sp:
                continue
            e0 = {k: sp[0][k] for k in ("e_dft_h", "e_casscf_h", "e_nevpt2_h")}
            levels = {}
            for lvl, key in (("dft", "e_dft_h"), ("casscf", "e_casscf_h"),
                             ("nevpt2", "e_nevpt2_h")):
                series = [(p["idx"], (p[key] - e0[key]) * HK) for p in sp]
                bar, interior = _barrier(series)
                levels[lvl] = {"barrier_kcal":
                               round(bar, 2) if bar is not None else None,
                               "interior_in_window": interior,
                               "rel_kcal": [round(v, 2) for _, v in series],
                               "idx": [i for i, _ in series]}
            row[sub] = levels
        out["per_field"][f"{f:+.3f}"] = row
        if len(row) == 2:
            d = {}
            for lvl in ("dft", "casscf", "nevpt2"):
                a, b = row["ch4"][lvl], row["c2h6"][lvl]
                ok = (a["interior_in_window"] and b["interior_in_window"]
                      and a["barrier_kcal"] is not None
                      and b["barrier_kcal"] is not None)
                d[lvl] = (round(b["barrier_kcal"] - a["barrier_kcal"], 2)
                          if ok else "NO-VERDICT (edge/missing)")
            dde[f"{f:+.3f}"] = d
    out["ddE"] = dde
    # наклон dΔΔE‡/dF по двум полям, если оба с вердиктом
    if len(FIELDS) == 2 and all(isinstance(
            dde.get(f"{f:+.3f}", {}).get("nevpt2"), (int, float))
            for f in FIELDS):
        f1, f2 = sorted(FIELDS)
        for lvl in ("dft", "casscf", "nevpt2"):
            v1 = dde[f"{f1:+.3f}"][lvl]
            v2 = dde[f"{f2:+.3f}"][lvl]
            if isinstance(v1, (int, float)) and isinstance(v2, (int, float)):
                slope = (v2 - v1) / (f2 - f1)
                out.setdefault("slope_kcal_per_au", {})[lvl] = round(slope, 1)
                out.setdefault("slope_kcal_per_V_A", {})[lvl] = round(
                    slope / AU_FIELD_V_A, 2)
    json.dump(out, open(RES, "w"), indent=1)
    print(f"[readout] ddE(F): {json.dumps(dde, ensure_ascii=False)}", flush=True)
    print(f"[readout] slope: {out.get('slope_kcal_per_V_A')}", flush=True)


# ------------------------------------------------- v2: цепочка по полю
# Урок v1: холодный старт на каждой (точка, поле) сажает CASSCF на разные
# ветки (контроль F=0 на TS разошёлся с профилем на ~26 ккал; +F профили
# краевые; UKS под −F схлопнулся в другое состояние, rel −180). Лекарство —
# warm-chain ПО ПОЛЮ: на каждой точке якорь F=0 (гейт: NEVPT2 против
# закоммиченного профиля), затем ±F тёплым стартом от плотности ROHF/UKS и
# орбиталей CASSCF якоря. Геометрия и базис между полями идентичны — орбитали
# передаются без проекции.
SCAN2 = os.path.join(DIR, "ocm_field_scan2.json")
RES2 = os.path.join(DIR, "ocm_field_dde2_results.json")


def _rohf_warm(mol, dm0=None):
    mf = scf.ROHF(mol)
    mf.conv_tol = 1e-8
    mfn = scf.newton(mf)
    mfn.kernel(dm0=dm0) if dm0 is not None else mfn.kernel()
    if not mfn.converged:
        mf2 = scf.ROHF(mol)
        mf2.level_shift = 0.4
        mf2.max_cycle = 300
        mf2.kernel(dm0=dm0)
        mfn = scf.newton(mf2)
        mfn.kernel(mf2.make_rdm1())
    return mfn


def _cas_pt2(mf, mo_start, natorb_first):
    """CASSCF+NEVPT2 на данных орбиталях; natorb_first=True — через CASCI."""
    from pyscf import fci, mcscf
    from pyscf.mrpt import NEVPT
    ne, no = sel.CAS_TARGET
    ss = sel.SPIN / 2 * (sel.SPIN / 2 + 1)
    mo = mo_start
    if natorb_first:
        mc0 = mcscf.CASCI(mf, no, ne)
        fci.addons.fix_spin_(mc0.fcisolver, shift=0.2, ss=ss)
        mc0.kernel(mo)
        mo = mc0.cas_natorb()[0]
    mc = mcscf.CASSCF(mf, no, ne)
    fci.addons.fix_spin_(mc.fcisolver, shift=0.2, ss=ss)
    mc.max_cycle_macro = 150
    mc.conv_tol = 1e-7
    mc.kernel(mo)
    ept = NEVPT(mc).kernel()
    return mc, float(mc.e_tot), float(mc.e_tot + ept)


def stage_field2():
    global FIELD
    pts = []
    if os.path.exists(SCAN2):
        try:
            pts = json.load(open(SCAN2))["pts"]
            print(f"[field2] resume: {len(pts)} записей", flush=True)
        except Exception:
            pts = []
    done = {(p["sub"], p["idx"], round(p["field"], 4)) for p in pts}

    def save():
        json.dump({"pts": pts}, open(SCAN2, "w"), indent=1)

    for sub in ("ch4", "c2h6"):
        win, its, _ = robust_window(sub)
        prof = json.load(open(os.path.join(DIR,
                              f"ocm_mnw_emb_cr_{sub}_profile.json")))
        g = json.load(open(os.path.join(DIR,
                           f"ocm_mnw_emb_cr_{sub}_geom.json")))
        print(f"[field2] {sub}: окно {win}", flush=True)
        for i in win:
            need = [f for f in [0.0] + FIELDS
                    if (sub, i, round(f, 4)) not in done]
            if not need:
                continue
            atoms = [(s, tuple(c)) for s, c in g["points"][i]["atoms"]]
            axis = site_axis(atoms)
            mol = sel.build_mol(atoms)
            mol.max_memory = max(12000, sel.MAXMEM)
            # якорь F=0
            FIELD = np.zeros(3)
            t0 = time.time()
            mfd0 = sel.uks(sel.build_mol(atoms))
            mf0 = _rohf_warm(mol)
            mo0, _ = sel.forced_avas(mf0)
            mc, e_cas0, e_pt0 = _cas_pt2(mf0, mo0, natorb_first=True)
            anchor = {"sub": sub, "idx": i, "field": 0.0,
                      "e_dft_h": float(mfd0.e_tot),
                      "e_casscf_h": e_cas0, "e_nevpt2_h": e_pt0,
                      "anchor_diff_kcal": round(
                          (e_pt0 - prof["points"][i]["e_nevpt2_h"]) * HK, 3),
                      "wall_s": round(time.time() - t0, 1)}
            if (sub, i, 0.0) not in done:
                pts.append(anchor)
                save()
            print(f"[field2] {sub} idx={i} F=0: NEVPT2 {e_pt0:.6f} "
                  f"(Δпрофиль {anchor['anchor_diff_kcal']} ккал)", flush=True)
            dm_r, dm_d, mo_c = mf0.make_rdm1(), mfd0.make_rdm1(), mc.mo_coeff
            for f in FIELDS:
                if (sub, i, round(f, 4)) in done:
                    continue
                FIELD = f * axis
                nf = nuc_field(atoms)
                t0 = time.time()
                mfd = sel.uks(sel.build_mol(atoms), dm0=dm_d)
                mf = _rohf_warm(mol, dm0=dm_r)
                _, e_cas, e_pt = _cas_pt2(mf, mo_c, natorb_first=False)
                pts.append({"sub": sub, "idx": i, "field": f,
                            "e_dft_h": float(mfd.e_tot) + nf,
                            "e_casscf_h": e_cas + nf, "e_nevpt2_h": e_pt + nf,
                            "rohf_conv": bool(mf.converged),
                            "wall_s": round(time.time() - t0, 1)})
                save()
                print(f"[field2] {sub} idx={i} F={f:+.3f}: "
                      f"NEVPT2 {e_pt + nf:.6f}", flush=True)
    print("[field2] цепочка завершена", flush=True)


def stage_readout2():
    pts = json.load(open(SCAN2))["pts"]
    allf = [0.0] + FIELDS
    out = {"model": "embedded Cr, frozen-frame SVP geometries",
           "protocol": "v2 field-chain: anchor F=0 (gate vs committed "
                       "profile) -> ±F warm-started (dm + CASSCF mo); "
                       "interior-in-window rule",
           "fields_au": allf,
           "anchor_gate_kcal": {}, "per_field": {}, "ddE": {}}
    for sub in ("ch4", "c2h6"):
        out["anchor_gate_kcal"][sub] = {
            str(p["idx"]): p.get("anchor_diff_kcal") for p in pts
            if p["sub"] == sub and p["field"] == 0.0}
    dde = {}
    for f in allf:
        row = {}
        for sub in ("ch4", "c2h6"):
            sp = sorted([p for p in pts
                         if p["sub"] == sub and p["field"] == f],
                        key=lambda p: p["idx"])
            if len(sp) < 3:
                continue
            levels = {}
            for lvl, key in (("dft", "e_dft_h"), ("casscf", "e_casscf_h"),
                             ("nevpt2", "e_nevpt2_h")):
                e0 = sp[0][key]
                series = [(p["idx"], (p[key] - e0) * HK) for p in sp]
                bar, interior = _barrier(series)
                levels[lvl] = {"barrier_kcal":
                               round(bar, 2) if bar is not None else None,
                               "interior_in_window": interior,
                               "rel_kcal": [round(v, 2) for _, v in series]}
            row[sub] = levels
        out["per_field"][f"{f:+.3f}"] = row
        if len(row) == 2:
            d = {}
            for lvl in ("dft", "casscf", "nevpt2"):
                a, b = row["ch4"][lvl], row["c2h6"][lvl]
                ok = (a["interior_in_window"] and b["interior_in_window"]
                      and a["barrier_kcal"] is not None
                      and b["barrier_kcal"] is not None)
                d[lvl] = (round(b["barrier_kcal"] - a["barrier_kcal"], 2)
                          if ok else "NO-VERDICT (edge/missing)")
            dde[f"{f:+.3f}"] = d
    out["ddE"] = dde
    ff = sorted(FIELDS)
    for lvl in ("dft", "casscf", "nevpt2"):
        vs = [dde.get(f"{x:+.3f}", {}).get(lvl) for x in ff]
        if len(vs) == 2 and all(isinstance(v, (int, float)) for v in vs):
            slope = (vs[1] - vs[0]) / (ff[1] - ff[0])
            out.setdefault("slope_kcal_per_au", {})[lvl] = round(slope, 1)
            out.setdefault("slope_kcal_per_V_A", {})[lvl] = round(
                slope / AU_FIELD_V_A, 2)
    json.dump(out, open(RES2, "w"), indent=1)
    print(f"[readout2] anchor gate: {out['anchor_gate_kcal']}", flush=True)
    print(f"[readout2] ddE(F): {json.dumps(dde, ensure_ascii=False)}",
          flush=True)
    print(f"[readout2] slope: {out.get('slope_kcal_per_V_A')}", flush=True)


# -------------------------------------------- v3: якоря warm-chain вдоль пути
# Диагноз v2: холодные ЯКОРЯ F=0 в серединных точках садятся на CASSCF-ветки
# на 16–26 ккал выше профильных (c2h6 idx4 — ROHF-ветка +75), потому что
# закоммиченный профиль считался warm-chain ВДОЛЬ ПУТИ. v3 воспроизводит
# профильный протокол: якорь точки i стартует с dm и CASSCF-орбиталей якоря
# точки i-1; гейт против профиля; при Δ>2 ккал — фолбэк на холодный
# AVAS-путь, берётся НИЖНЕЕ решение (обе энергии в чекпойнте).
SCAN3 = os.path.join(DIR, "ocm_field_scan3.json")
RES3 = os.path.join(DIR, "ocm_field_dde3_results.json")


def stage_field3():
    global FIELD
    pts = []
    if os.path.exists(SCAN3):
        try:
            pts = json.load(open(SCAN3))["pts"]
            print(f"[field3] resume: {len(pts)} записей", flush=True)
        except Exception:
            pts = []
    done = {(p["sub"], p["idx"], round(p["field"], 4)) for p in pts}

    def save():
        json.dump({"pts": pts}, open(SCAN3, "w"), indent=1)

    for sub in ("ch4", "c2h6"):
        win, its, _ = robust_window(sub)
        prof = json.load(open(os.path.join(DIR,
                              f"ocm_mnw_emb_cr_{sub}_profile.json")))
        g = json.load(open(os.path.join(DIR,
                           f"ocm_mnw_emb_cr_{sub}_geom.json")))
        print(f"[field3] {sub}: окно {win}", flush=True)
        dm_chain, mo_chain, dmd_chain, mol_prev = None, None, None, None
        for i in win:
            atoms = [(s, tuple(c)) for s, c in g["points"][i]["atoms"]]
            axis = site_axis(atoms)
            mol = sel.build_mol(atoms)
            mol.max_memory = max(12000, sel.MAXMEM)
            # якорь F=0 — warm-chain вдоль пути. v4: перенос орбиталей между
            # ГЕОМЕТРИЯМИ только через project_init_guess (v3 передавал
            # mo_coeff напрямую — AO-базис едет за атомами, CASSCF
            # коллапсировал на тысячи ккал). Санити: gate < −10 ккал =
            # коллапс, решение отбрасывается.
            from pyscf import mcscf as _mc
            FIELD = np.zeros(3)
            t0 = time.time()
            mfd0 = sel.uks(sel.build_mol(atoms), dm0=dmd_chain)
            mf0 = _rohf_warm(mol, dm0=dm_chain)
            eref = prof["points"][i]["e_nevpt2_h"]
            cand = []
            if mo_chain is not None:
                try:
                    ne_, no_ = sel.CAS_TARGET
                    mct = _mc.CASSCF(mf0, no_, ne_)
                    mo_proj = _mc.addons.project_init_guess(
                        mct, mo_chain, prev_mol=mol_prev)
                    mc1, ec1, ep1 = _cas_pt2(mf0, mo_proj,
                                             natorb_first=False)
                    cand.append(("chained-proj", mc1, ec1, ep1,
                                 (ep1 - eref) * HK))
                except Exception as ex:
                    print(f"[field3] proj-chain fail idx={i}: {ex}",
                          flush=True)
            need_cold = (not cand or cand[0][4] > 2.0
                         or cand[0][4] < -10.0)
            if need_cold:
                moc, _ = sel.forced_avas(mf0)
                mc2, ec2, ep2 = _cas_pt2(mf0, moc, natorb_first=True)
                cand.append(("cold", mc2, ec2, ep2, (ep2 - eref) * HK))
            sane = [c for c in cand if c[4] > -10.0]
            pick = min(sane or cand, key=lambda c: c[3])
            branch, mc, e_cas0, e_pt0, gate = pick
            if not sane:
                branch += "-COLLAPSED"
            if (sub, i, 0.0) not in done:
                pts.append({"sub": sub, "idx": i, "field": 0.0,
                            "e_dft_h": float(mfd0.e_tot),
                            "e_casscf_h": e_cas0, "e_nevpt2_h": e_pt0,
                            "anchor_diff_kcal": round(gate, 3),
                            "branch": branch,
                            "wall_s": round(time.time() - t0, 1)})
                save()
            print(f"[field3] {sub} idx={i} F=0 [{branch}]: NEVPT2 "
                  f"{e_pt0:.6f} (Δпрофиль {gate:.2f} ккал)", flush=True)
            dm_chain, dmd_chain = mf0.make_rdm1(), mfd0.make_rdm1()
            mo_chain, mol_prev = mc.mo_coeff, mol
            for f in FIELDS:
                if (sub, i, round(f, 4)) in done:
                    continue
                FIELD = f * axis
                nf = nuc_field(atoms)
                t0 = time.time()
                mfd = sel.uks(sel.build_mol(atoms), dm0=dmd_chain)
                mf = _rohf_warm(mol, dm0=dm_chain)
                _, e_cas, e_pt = _cas_pt2(mf, mo_chain, natorb_first=False)
                pts.append({"sub": sub, "idx": i, "field": f,
                            "e_dft_h": float(mfd.e_tot) + nf,
                            "e_casscf_h": e_cas + nf,
                            "e_nevpt2_h": e_pt + nf,
                            "rohf_conv": bool(mf.converged),
                            "wall_s": round(time.time() - t0, 1)})
                save()
                print(f"[field3] {sub} idx={i} F={f:+.3f}: "
                      f"NEVPT2 {e_pt + nf:.6f}", flush=True)
    print("[field3] цепочка завершена", flush=True)


def stage_readout3():
    global SCAN2, RES2
    scan2_orig, res2_orig = SCAN2, RES2
    # readout2 читает SCAN2/RES2 — подменяем на v3-файлы
    globals()["SCAN2"], globals()["RES2"] = SCAN3, RES3
    try:
        stage_readout2()
    finally:
        globals()["SCAN2"], globals()["RES2"] = scan2_orig, res2_orig


if __name__ == "__main__":
    {"field": stage_field, "readout": stage_readout,
     "field2": stage_field2, "readout2": stage_readout2,
     "field3": stage_field3, "readout3": stage_readout3}[sys.argv[1]]()
