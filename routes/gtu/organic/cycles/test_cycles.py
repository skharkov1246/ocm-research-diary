"""Проверки cycles.py (pytest, около 1–2 мин). Числа сверки — из реальных расчётов здесь же и в gtu_cycle."""
import json
import math
import os
import sys

import numpy as np
import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import cycles as C  # noqa: E402

BASE = dict(rec=0.0, ic=False, rh=False, ibc=False, pgc="none", tcr="none")


def test_laminar_constants_from_finite_differences():
    """плоский канал: f·Re = 24, Nu_H = 140/17, Nu_T = 7,541 (учебные значения Shah & London — посчитаны заново)"""
    L = C.plate_laminar_constants()
    assert abs(L["fRe"] - 24.0) < 0.02
    assert abs(L["NuH"] - 140 / 17) < 0.01
    assert abs(L["NuT"] - 7.5407) < 0.01


def test_simple_cycle_matches_gtu_cycle_and_frame3():
    x = C.crosscheck_gtu_cycle()["simple"]
    a, b, f = x["cycles"], x["gtu_cycle"], x["frame3_json"]
    assert abs(a["P_kW"] / b["P_kW"] - 1) < 0.0015          # 0,09 % — дробление политропы на 40 / 2×40 шагов
    assert abs(a["eff"] - b["eff"]) < 0.0005
    assert abs(a["EGT_C"] - b["EGT_C"]) < 0.5
    assert abs(b["P_kW"] - f["P_kW"]) < 0.05                 # gtu_cycle воспроизводит frame3_point_results.json


def test_recuperator_matches_gtu_cycle_recup():
    x = C.crosscheck_gtu_cycle()["recup"]
    assert abs(x["cycles"]["eff"] - x["gtu_cycle"]["eff"]) < 0.003
    assert abs(x["cycles"]["T3r_C"] - x["gtu_cycle"]["T3r_C"]) < 5.0
    assert abs(x["cycles"]["w_kJkg"] - x["gtu_cycle"]["w_kJkg"]) < 1.0


def test_burn_energy_balance():
    air = C.St(0.985, 700.0, 13e5, C.Y_AIR)
    h_f = C.hs(C.St(1, C.T0, 16e5, C.Y_CH4))[0]
    mf, p = C.burn(air, C.Y_CH4, h_f, 1273.15, 12.5e5, C.LHV_CH4)
    lhs = air.m * C.hs(air)[0] + mf * h_f - (1 - C.A["eta_b"]) * mf * C.LHV_CH4
    assert abs(lhs - p.m * C.hs(p)[0]) < 1e-6 * abs(lhs) + 1.0
    assert abs(p.T - 1273.15) < 0.01


@pytest.mark.parametrize("topo,PR", [
    (dict(BASE), 13.5), (dict(BASE, rec=0.9), 5.0), (dict(BASE, ic=True), 13.5), (dict(BASE, rh=True, rec=0.85), 8.0),
    (dict(BASE, ibc=True), 13.5), (dict(BASE, pgc="real"), 13.5), (dict(BASE, tcr="steam"), 13.5),
    (dict(BASE, tcr="egr", rec=0.85, ibc=True), 6.0), (dict(BASE, tcr="fh", rh=True), 13.5)])
def test_exergy_two_ways(topo, PR):
    """разрушение эксергии по узлам: баланс эксергии (A) = T0·S_gen (B); баланс всего цикла закрыт"""
    r = C.evaluate(topo, PR, detail=True)
    assert r["ok"], r.get("why")
    e = C.exergy_account(r)
    assert e["AB_max_diff_kW"] < 0.5
    assert abs(e["balance_residual_pct"]) < 0.01
    assert min(x["X_A_kW"] for x in e["rows"]) > -0.5      # второй закон: разрушение не отрицательно
    assert max(abs(x["energy_residual_kW"]) for x in e["rows"]) < 0.5


def test_recuperator_monotonic_eff_and_size():
    effs, vols = [], []
    for eps in (0.70, 0.80, 0.90):
        r = C.evaluate(dict(BASE, rec=eps), 6.0, detail=True)
        effs.append(r["eff"])
        vols.append(C.hx_list(r, "lpbf")[0]["V_core_m3"])
    assert effs[0] < effs[1] < effs[2]
    assert vols[0] < vols[1] < vols[2]


def test_hx_sizing_consistency():
    """ламинарность, потеря давления не больше допустимой, объём = A·(Dh/2 + t)"""
    r = C.evaluate(dict(BASE, rec=0.85), 6.0, detail=True)
    h = C.hx_list(r, "lpbf")[0]
    t = C.HX_TECH["lpbf"]
    gin = r["_streams"]["g_in"]
    assert h["laminar"]
    assert h["dP_hot_Pa"] <= C.A["dp_core_share"] * C.A["dp_rg"] * gin.P * 1.0001
    assert abs(h["V_core_m3"] - h["A_m2"] * (t["Dh"] / 2 + t["t"])) < 1e-9


def test_equilibrium_constants_vs_empirical():
    """константы равновесия из Cantera против эмпирических формул (Xu — Froment, по памяти) — в пределах 7 %"""
    for T in (700.0, 900.0, 1100.0):
        a, b = C.K_eq_bar(T), C.K_eq_xf_empirical(T)
        for x, y in zip(a, b):
            assert abs(math.log(x / y)) < 0.07


def test_kinetics_converge_to_equilibrium():
    T, P = 923.15, 17e5
    X = np.zeros(len(C.SP))
    X[C.I["CH4"]], X[C.I["H2O"]] = 1 / 3, 2 / 3
    Y = X * C.MW / np.dot(X, C.MW)
    Ye, _ = C.equil(Y, T, P)
    Xeq = 1 - (Ye[C.I["CH4"]] / C.MW[C.I["CH4"]]) / (Y[C.I["CH4"]] / C.MW[C.I["CH4"]])
    W, Xc, _ = C.xf_pfr(X, 3.0, T, P, 50.0)
    assert abs(Xc[-1] - Xeq) < 0.005
    assert np.all(np.diff(Xc) > -1e-6)


def test_carbon_activity_matches_multiphase_equilibrium():
    import cantera as ct
    gr = ct.Solution("graphite.yaml")
    for X, T, P in (("CH4:1,H2O:2", 923.15, 17e5), ("CH4:1,CO2:1", 1123.15, 17e5)):
        Y = C.Y_of(X)
        Ye, _ = C.equil(Y, T, P)
        aM, aB = C.carbon_activity(Ye, T, P)
        assert abs(aM / aB - 1) < 1e-6                       # в равновесном газе оба пути дают одно
        g = C.gas()
        g.TPY = T, P, Y
        m = ct.Mixture([(g, 1.0), (gr, 0.0)])
        m.T, m.P = T, P
        m.equilibrate("TP")
        assert (m.phase_moles()[1] > 1e-6) == (aM > 1.0)


def test_egr_reforming_does_not_raise_lhv():
    """кислород выхлопа (обеднённая смесь) сжигает часть метана: теплотворность не растёт"""
    for R in (1.0, 2.0, 4.0):
        r = C.evaluate(dict(BASE, tcr="egr"), 13.5, dict(R_egr=R))
        assert r["ok"]
        assert r["fuel_info"]["lhv_gain"] <= 0.0


def test_pgc_bounds():
    r0 = C.evaluate(dict(BASE), 13.5)
    rr = C.evaluate(dict(BASE, pgc="real"), 13.5)
    ri = C.evaluate(dict(BASE, pgc="ideal"), 13.5)
    assert ri["PG_ideal"] > 1.0
    assert ri["eff"] > rr["eff"]
    assert ri["eff"] > r0["eff"]


def test_pareto_function():
    rows = [dict(a=1, b=1), dict(a=2, b=2), dict(a=2, b=1), dict(a=0, b=0)]
    idx = C.pareto(rows, [("a", 1), ("b", -1)])
    assert set(idx) == {2, 3}


RES = os.path.join(HERE, "cycles_results.json")


@pytest.mark.skipif(not os.path.exists(RES), reason="нет cycles_results.json (запустить cycles.py)")
def test_results_file():
    d = json.load(open(RES, encoding="utf-8"))
    for k in ("method", "software", "date", "inputs", "assumptions"):
        assert d[k]
    assert d["date"] == "2026-11-14"
    assert d["n_topologies"] == 576
    fr = d["pareto_lpbf"]
    assert fr and all(r["limits_ok"] and not r["bound"] for r in fr)
    idx = C.pareto(fr, [("eff_el", 1), ("added", -1), ("V_added_m3", -1)])
    assert len(idx) == len(fr)                                 # фронт сам себя не доминирует
    base = d["base"]["opt"]
    assert base["limits"]["all"]
    assert all(abs(x["X_A_kW"] - x["X_B_kW"]) < 0.5 for x in base["exergy"]["rows"])
    for o in d["results"]:
        r = o["opt"]
        if r.get("ok"):
            assert r["parts"]["lpbf"]["TRL"] <= 9
            if o["topo"]["tcr"] == "steam":
                assert r["parts"]["lpbf"]["pnp"] == "нет"         # вода в пакете
