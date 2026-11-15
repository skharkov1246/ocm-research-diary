"""Проверки эксергетического анализа KV-3: баланс, совпадение с моделью цикла, независимые формулы,
повторяемость. Запуск: python3 -m pytest -q test_exergy.py (около 15 с)."""
import json
import math
import os
import xml.etree.ElementTree as ET

import pytest

import exergy as X

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.join(HERE, "exergy_results.json")


@pytest.fixture(scope="module")
def setup():
    f3 = X.load_frame3()
    env = X.Env()
    acc, ch = X.run_simple(X.base_kw(f3), env, f3)
    return f3, env, acc, ch


@pytest.fixture(scope="module")
def res():
    return json.load(open(RES, encoding="utf-8"))


def test_chain_equals_gtu_cycle(setup):
    """Цепочка состояний совпадает с gtu_cycle.cycle и с расчётной точкой frame3."""
    f3, env, acc, ch = setup
    mc = acc["model_check"]
    assert mc["P"] < 1e-9 and mc["Wf"] < 1e-9 and mc["EGT_K"] < 1e-6
    assert abs(acc["P_shaft_kW"] - f3["P_shaft_kW_ref"]) < 0.05
    assert abs(acc["eff_shaft"] - f3["eff_ref"]) < 1e-5


def test_balance_closed(setup):
    """Баланс эксергии (способ B, T0·S_gen) закрыт до 0,5 % эксергии топлива; A и B узлов совпадают."""
    f3, env, acc, ch = setup
    b = acc["balance"]
    assert abs(b["residual_B_pct_fuel_ex"]) <= 0.5
    assert abs(b["residual_A_pct_fuel_ex"]) <= 0.5
    for k, v in acc["nodes"].items():
        assert abs(v["A_minus_B_pct_fuel_ex"]) < 1e-6, k
    # доли: продукт + разрушение + уходящие потоки = 100 % + вход воздуха - невязка
    t = acc["totals"]
    s = t["product_pct"] + t["destruction_pct"] + t["outflow_pct"]
    assert abs(s - (100 + t["air_in_pct"] - b["residual_B_pct_fuel_ex"])) < 1e-6


def test_exergetic_efficiency_two_ways(setup):
    """Эксергетический КПД = КПД по LHV · LHV / эксергия топлива на кг (другой путь к тому же числу)."""
    f3, env, acc, ch = setup
    via_lhv = acc["eff_el"] * X.C.LHV / (acc["fuel"]["ex_MJ_kg"] * 1e6)
    assert abs(via_lhv - acc["ex_eff_el"]) < 1e-9


def test_pressure_loss_formula(setup):
    """Потеря давления камеры: идеальный газ при h = const — T0·M·R/μ·ln(p3/p4) (независимая формула)."""
    f3, env, acc, ch = setup
    st, m = ch["st"], ch["m"]
    g = X.ct.Solution("gri30.yaml")
    g.TPY = st["pa"]["T"], st["pa"]["P"], st["pa"]["Y"]
    R = X.ct.gas_constant / g.mean_molecular_weight
    D = X.T0 * m["M"] * R * math.log(st["pa"]["P"] / st["pb"]["P"]) / 1e3
    assert abs(D - acc["nodes"]["comb_dp"]["kW"]) / D < 1e-6


def test_mechanical_losses(setup):
    """Механика, редуктор, генератор — работа, ушедшая в теплоту: прямой пересчёт по КПД."""
    f3, env, acc, ch = setup
    w = acc["works_kW"]
    assert abs(acc["nodes"]["mech"]["kW"] - w["W_turb"] * (1 - 0.99)) < 1e-6
    P = acc["P_shaft_kW"]
    assert abs(acc["nodes"]["gearbox"]["kW"] - P * (1 - f3["eta_gb"])) < 1e-6
    assert abs(acc["nodes"]["generator"]["kW"] + acc["nodes"]["gearbox"]["kW"] + acc["P_el_kW"] - P) < 1e-6


def test_turbine_stage_split(setup):
    """Разрушение в турбине по трём ступеням frame3 сходится с цепочкой (две части по 40 шагов) до 0,5 %."""
    f3, env, acc, ch = setup
    ts = X.turbine_stages(ch, f3["turb_ER"], f3["turb_eta"])
    tot = acc["nodes"]["turbine_gg"]["kW"] + acc["nodes"]["turbine_pt"]["kW"]
    assert abs(ts["sum_D_kW"] - tot) / tot < 0.005


def test_compressor_split(setup):
    """Две ступени: сумма разрушения = узел компрессора; КПД ступени 2 близок к frame3."""
    f3, env, acc, ch = setup
    cs = acc["compressor_stages"]
    assert abs(cs["D1_kW"] + cs["D2_kW"] - acc["nodes"]["compressor"]["kW"]) < 1e-6
    assert abs(cs["eta_p2_matched"] - f3["eta_p2"]) < 0.005


def test_fuel_exergy(setup):
    f3, env, acc, ch = setup
    fx = X.fuel_exergy_check(env)
    assert fx["diff_rel"] < 1e-9
    assert 1.03 < fx["ratio_ex_LHV_model"] < 1.05
    assert abs(fx["LHV_cantera_MJ_kg"] - X.C.LHV / 1e6) / (X.C.LHV / 1e6) < 0.002
    assert abs(fx["szargut_env_298K_kJ_mol"] - 831.65) / 831.65 < 0.005


def test_cooling_model():
    f3 = X.load_frame3()
    x, info = X.cooling_model(f3, 410.37)
    assert x(1000.0, "cast") == (0.0, 0.0)
    xn, xr = x(1200.0, "cast")
    assert abs(xn - 0.055) < 1e-12 and abs(xr - 0.055) < 1e-12
    pn, pr = x(1150.0, "printed")
    cn, cr = x(1150.0, "cast")
    assert abs(pn / cn - 0.5 / 0.75) < 1e-12 and abs(pr / cr - 0.5 / 0.75) < 1e-12


def test_repeatability(setup, res):
    """Повторный расчёт расчётной точки совпадает с exergy_results.json (округление 6 знаков)."""
    f3, env, acc, ch = setup
    b = res["variants"]["base"]
    for k, v in acc["nodes"].items():
        assert abs(v["kW"] - b["nodes"][k]["kW"]) < 2e-6 * max(1.0, abs(v["kW"])), k
    assert abs(acc["ex_eff_el"] - b["ex_eff_el"]) < 1e-6
    acc2, _ = X.run_simple(X.base_kw(f3), env, f3)
    assert acc2["nodes"]["comb_reaction"]["kW"] == acc["nodes"]["comb_reaction"]["kW"]


def test_results_file(res):
    for k in ("method", "software", "date", "inputs", "assumptions"):
        assert res.get(k), k
    assert res["date"] == "2026-11-14"
    for name, v in res["variants"].items():
        assert v["balance"]["closed"] and abs(v["balance"]["residual_B_pct_fuel_ex"]) <= 0.5, name
        if v.get("model_check"):
            assert v["model_check"]["P"] < 1e-6 and v["model_check"]["Wf"] < 1e-6, name
    rb = res["variants"]["recup"]
    assert rb["recuperator"]["T_gas_in_C"] <= 700.0 + 1e-6
    for r in res["reserve"]:
        assert r.get("balance_pct") is None or abs(r["balance_pct"]) <= 0.5


def test_svgs(res):
    for f in res["svg"].values():
        root = ET.parse(os.path.join(HERE, f)).getroot()
        assert root.tag.endswith("svg")
        txt = open(os.path.join(HERE, f), encoding="utf-8").read()
        assert "nan" not in txt.lower().replace("nanosec", "")
