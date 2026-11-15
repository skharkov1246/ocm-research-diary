#!/usr/bin/env python3
"""Быстрая оценка напряжений диска (колеса) переменной толщины — для сотен тысяч оценок в поиске формы.

Модель: тонкий диск, плоское напряжённое состояние (σz = 0), осесимметрия; толщина h(r), плотность с размазанными
лопатками ρ_eff(r) = ρ·(1 + φ·hb/h) — как routes/gtu/platform/fea/disk.py; поле температур T(r), E(T), α(T), ν — из
materials.json через disk.load_material (одна точка правды с CalculiX). Метод — одномерные конечные элементы по радиусу
(линейные, 2 точки Гаусса), уравнение равновесия в слабой форме:
    ∫ (σr δεr + σθ δεθ) h r dr = ∫ ρ_eff h ω² r² δu dr + p·h_o·r_o·δu(r_o).
Выход: σr, σθ, Мизес (плоское состояние) в точках Гаусса, наибольшее σ/σ0,2(T), средние окружные напряжения и запас
по разрушающей частоте (как disk.py: n_разр/n = sqrt(k·∫σв(T)h dr / ∫σθ h dr)), ресурс по Ларсону — Миллеру (если в
материале есть блок lmp), проверка равновесия полудиска.

Ограничение (оценка): толстая ступица колеса (h ~ r) — не тонкий диск, у оси CAX8 даёт другое распределение по
толщине. Поэтому быстрый решатель — только для поиска; элиты проверяются CalculiX (disk.py), расхождение пишется
(test_blades.py сверяет с frame3_fea_results.json и profile_opt_comp2.json).
"""
from __future__ import annotations

import math
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
GTU = os.path.dirname(os.path.dirname(HERE))
FEA = os.path.join(GTU, "platform", "fea")
if FEA not in sys.path:
    sys.path.insert(0, FEA)
import disk as DK  # noqa: E402

_MAT = {}


def material(key):
    if key not in _MAT:
        _MAT[key] = DK.load_material(key)
    return _MAT[key]


def _tab(tab, x):
    a = np.asarray(tab, float)
    return np.interp(x, a[:, 0], a[:, 1])


def solve(r_nodes, h_nodes, rho_fac, T_nodes, mat_key, rpm, p_rim_Pa=0.0, T_ref=20.0, burst_k=1.0, k_peak=1.0):
    """r_nodes возрастают (r_nodes[0] = 0 — сплошной диск, u(0) = 0; иначе — свободная граница отверстия), h, rho_fac
    (множитель плотности от лопаток), T — в узлах. p_rim_Pa — растягивающее радиальное напряжение на ободе (лопатки
    турбины). k_peak — поправка пика Мизеса толстой ступицы (калибровка по CalculiX, calibrate()); средние напряжения,
    запас по разрушению и равновесие от неё не зависят. Сборка векторно, трёхдиагональная система (solve_banded)."""
    from scipy.linalg import solve_banded
    m = material(mat_key)
    r = np.asarray(r_nodes, float)
    h = np.asarray(h_nodes, float)
    rf = np.asarray(rho_fac, float) * np.ones_like(r)
    Tn = np.asarray(T_nodes, float) * np.ones_like(r)
    nu = float(m["nu"])
    rho = float(m["density"])
    om = 2 * math.pi * rpm / 60.0
    n = len(r)
    L = np.diff(r)
    g = 1 / math.sqrt(3)
    Ka = np.zeros(n)          # диагональ
    Kb = np.zeros(n - 1)      # наддиагональ (симметрично)
    f = np.zeros(n)
    store = []
    for xi in (-g, g):
        N1, N2 = (1 - xi) / 2, (1 + xi) / 2
        rr = N1 * r[:-1] + N2 * r[1:]
        hh = N1 * h[:-1] + N2 * h[1:]
        TT = N1 * Tn[:-1] + N2 * Tn[1:]
        ff = N1 * rf[:-1] + N2 * rf[1:]
        E = _tab(m["E_GPa"], TT) * 1e9
        al = _tab(m["alpha_mean_1e6"], TT) * 1e-6
        c = E / (1 - nu * nu)
        dV = hh * rr * L / 2
        # B = [[-1/L, 1/L], [N1/r, N2/r]], D = c [[1, nu], [nu, 1]]
        b11, b12 = -1 / L, 1 / L
        b21, b22 = N1 / rr, N2 / rr
        k11 = c * (b11 * b11 + 2 * nu * b11 * b21 + b21 * b21) * dV
        k22 = c * (b12 * b12 + 2 * nu * b12 * b22 + b22 * b22) * dV
        k12 = c * (b11 * b12 + nu * (b11 * b22 + b21 * b12) + b21 * b22) * dV
        Ka[:-1] += k11
        Ka[1:] += k22
        Kb += k12
        th = al * (TT - T_ref) * (1 + nu) * c        # D @ [eth, eth] = c(1+ν)αΔT·[1, 1]
        f[:-1] += (b11 + b21) * th * dV + N1 * rho * ff * om * om * rr * dV
        f[1:] += (b12 + b22) * th * dV + N2 * rho * ff * om * om * rr * dV
        store.append((N1, N2, rr, hh, TT, E, al, L / 2))
    f[-1] += p_rim_Pa * h[-1] * r[-1]
    if r[0] <= 1e-12:
        Ka[0], f[0] = 1.0, 0.0
        Kb[0] = 0.0
    ab = np.zeros((3, n))
    ab[0, 1:] = Kb
    ab[1] = Ka
    ab[2, :-1] = Kb
    if r[0] <= 1e-12:
        ab[2, 0] = 0.0           # строка 1 не должна видеть u0 (u0 = 0)
    u = solve_banded((1, 1), ab, f)
    rg, hg, Tg, sr, st, wts = [], [], [], [], [], []
    for (N1, N2, rr, hh, TT, E, al, half) in store:
        er = np.diff(u) / L
        et = (N1 * u[:-1] + N2 * u[1:]) / rr
        th = al * (TT - T_ref)
        c = E / (1 - nu * nu)
        sr.append(c * (er + nu * et - (1 + nu) * th))
        st.append(c * (et + nu * er - (1 + nu) * th))
        rg.append(rr); hg.append(hh); Tg.append(TT); wts.append(hh * half)
    rg, hg, Tg, sr, st, wts = (np.concatenate(x) for x in (rg, hg, Tg, sr, st, wts))
    vm = np.sqrt(sr ** 2 - sr * st + st ** 2)
    ys = _tab(m["ys_MPa"], Tg) * 1e6
    uts = _tab(m["uts_MPa"], Tg) * 1e6
    yr = k_peak * vm / ys
    k = int(np.argmax(yr))
    lhs = float(np.sum(st * wts))
    rhs = float(np.sum(rho * np.interp(rg, r, rf) * om * om * rg * rg * wts)) + p_rim_Pa * r[-1] * h[-1]
    burst = math.sqrt(burst_k * float(np.sum(uts * wts)) / lhs) if lhs > 0 else float("inf")
    res = dict(vm_max_MPa=float(k_peak * vm.max() / 1e6), yield_ratio_max=float(yr[k]), yield_at_r=float(rg[k]),
               yield_at_T=float(Tg[k]), hoop_max_MPa=float(st.max() / 1e6), hoop_mean_MPa=lhs / float(np.sum(wts)) / 1e6,
               burst_margin=burst, equilibrium_rel=abs(lhs - rhs) / max(abs(rhs), 1e-300), u_rim_mm=float(u[-1] * 1e3),
               k_peak=k_peak)
    lmp = m.get("_lmp")
    if lmp is not None:
        life, ex = lmp.life_h_v(k_peak * vm / 1e6, Tg)
        j = int(np.argmin(life))
        res.update(creep_life_min_h=float(life[j]), creep_at_r=float(rg[j]), creep_at_T=float(Tg[j]),
                   creep_extrapolated=bool(ex[j]))
    return res


def mass_inertia(r, h, rho, extra_rho_fac=None):
    """масса и полярный момент инерции тела вращения толщины h(r): ∫ρ h 2πr dr, ∫ρ h 2πr³ dr"""
    r = np.asarray(r, float)
    h = np.asarray(h, float)
    fac = 1.0 if extra_rho_fac is None else np.asarray(extra_rho_fac, float)
    w = rho * fac * h * 2 * math.pi * r
    return float(np.trapezoid(w, r)), float(np.trapezoid(w * r * r, r))


# =========================================================================== калибровка пика по CalculiX
CAL_JSON = os.path.join(HERE, "disk1d_calibration.json")


def _profile_case(prof, T_tab, phi, n=301):
    from scipy.interpolate import PchipInterpolator
    r0, t0, hb = (np.asarray(prof[k], float) for k in ("r", "t", "hb"))
    rr = np.linspace(0.0, r0[-1], n)
    h = PchipInterpolator(r0, t0)(rr)
    fac = 1 + phi * np.interp(rr, r0, hb) / h
    T = np.interp(rr, T_tab["r"], T_tab["T"])
    return rr, h, fac, T


def calibrate(write=True):
    """k_peak = (σ/σ0,2 по CalculiX, CAX8 160×24) / (σ/σ0,2 тонкого диска) на шести готовых расчётах:
    frame3_fea_results.json (колёса ступеней 1 и 2 без отверстия, 100 и 120 %) и organic/structure/profile_opt_comp2.json
    (кованый профиль и гладкая доводка задней поверхности, 100 %, сетка 160×24). Тот же материал, частота, поле T, φ."""
    import json
    f3p = os.path.join(GTU, "frame3", "frame3_fea_results.json")
    pop = os.path.join(GTU, "organic", "structure", "profile_opt_comp2.json")
    d = json.load(open(f3p, encoding="utf-8"))
    rows = []
    for part in ("comp1", "comp2"):
        P = d["parts"][part]
        s = P["variants"]["solid"]
        rr, h, fac, T = _profile_case(s["profile"], s["T_C"], P["blade_fill"])
        for pct in ("100", "120"):
            c = s["cases"].get(pct)
            if not c:
                continue
            fast = solve(rr, h, fac, T, "ti64", c["rpm"])
            rows.append(dict(case=f"frame3 {part} solid {pct} %", ccx_yield_ratio=c["yield_ratio_max"],
                             fast_yield_ratio=fast["yield_ratio_max"], ccx_burst=c["burst_margin"], fast_burst=fast["burst_margin"],
                             ratio=c["yield_ratio_max"] / fast["yield_ratio_max"]))
    po = json.load(open(pop, encoding="utf-8"))
    T2 = d["parts"]["comp2"]["variants"]["solid"]["T_C"]
    for key, lab in (("base", "кованый профиль frame3 (80×12)"), ("check_160x24", "гладкая доводка (160×24)")):
        c = po[key]["c100"]
        rr, h, fac, T = _profile_case(c["profile"], T2, d["parts"]["comp2"]["blade_fill"])
        fast = solve(rr, h, fac, T, "ti64", 14950.0)
        mass, J = mass_inertia(rr, h, material("ti64")["density"])
        rows.append(dict(case=f"organic comp2 {lab}, 100 %", ccx_yield_ratio=c["yield_ratio"], fast_yield_ratio=fast["yield_ratio_max"],
                         ccx_burst=c["burst_margin"], fast_burst=fast["burst_margin"], ratio=c["yield_ratio"] / fast["yield_ratio_max"],
                         ccx_mass_kg=c["mass_kg"], fast_mass_kg=mass, ccx_J=c["J_kgm2"], fast_J=J))
    ks = np.array([r_["ratio"] for r_ in rows])
    out = dict(title="Калибровка пика Мизеса быстрого решателя диска по CalculiX", method=(
        "тонкий диск (плоское напряжённое состояние, 1D МКЭ 301 узел) против осесимметричной CAX8 disk.py на тех же профилях,"
        " полях T, частотах и φ; k_peak — среднее отношение σ/σ0,2; средние окружные напряжения и запас по разрушению"
        " определяются равновесием и от k_peak не зависят"),
        inputs=["routes/gtu/frame3/frame3_fea_results.json → parts.comp1/comp2.variants.solid (profile, T_C, cases)",
                "routes/gtu/organic/structure/profile_opt_comp2.json → base.c100, check_160x24.c100"],
        cases=rows, k_peak=float(ks.mean()), k_min=float(ks.min()), k_max=float(ks.max()),
        k_spread_rel=float((ks.max() - ks.min()) / ks.mean()),
        note=("все шесть колёс — толстая ступица без отверстия, пик у оси на задней стороне; для профилей с пиком в полотне"
              " поправка консервативна; элиты проверяются CalculiX"))
    if write:
        with open(CAL_JSON, "w", encoding="utf-8") as fo:
            json.dump(out, fo, ensure_ascii=False, indent=1)
    return out


def k_peak_cal():
    import json
    if os.path.exists(CAL_JSON):
        return float(json.load(open(CAL_JSON, encoding="utf-8"))["k_peak"])
    return float(calibrate(write=False)["k_peak"])


if __name__ == "__main__":
    import json
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--calibrate", action="store_true")
    a = ap.parse_args()
    if a.calibrate:
        print(json.dumps(calibrate(), ensure_ascii=False, indent=1))
