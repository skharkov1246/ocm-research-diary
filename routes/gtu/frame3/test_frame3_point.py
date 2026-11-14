"""Проверки расчётной точки газогенератора типоразмера № 3 (frame3_point.py). Запуск: python3 -m pytest -q test_frame3_point.py

Проверяется согласованность frame3_point_results.json и повторяемость расчёта: U2 = pi*D2*n/60, pi1*pi2 = pi,
U2 не выше предела материала, баланс мощности турбины, компрессора и вала (турбина — по таблицам gtu_explore,
компрессор и вал — по Cantera), одна ступень турбины не снимает нужную работу, нерасчётная модель на расчётном
режиме воспроизводит проект, повтор расчёта даёт те же числа."""
import json
import math
import os
import sys

import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import frame3_point as M  # noqa: E402  (импорт меняет текущую папку на папку модуля)

RES = json.load(open(os.path.join(HERE, "frame3_point_results.json")))
REL = 3e-5          # в JSON — 6 значащих цифр


def close(a, b, rel=REL):
    return abs(a - b) <= rel * max(abs(b), 1e-12)


def inp(field_part):
    return next(i["used"] for i in RES["inputs"] if field_part in i["field"])


@pytest.fixture(scope="module")
def dp():
    return M.design_point(RES["split"]["pi1"])


def test_u2_equals_pi_d2_n():
    n_out = inp(".мягкий")["n_out"]
    for key in ("stage1", "stage2"):
        s = RES[key]
        assert s["n_rpm"] == n_out
        assert close(s["U2_m_s"], math.pi * s["D2_m"] * s["n_rpm"] / 60)
    t = RES["turbine"]["options"][str(RES["turbine"]["selected_z"])]
    assert close(t["U_mean"], math.pi * t["Dm"] * t["n"] / 60)


def test_pressure_ratio_product():
    assert close(RES["stage1"]["PR"] * RES["stage2"]["PR"], RES["split"]["pi"])
    assert close(RES["split"]["pi1"] * RES["split"]["pi2"], RES["split"]["pi"])
    assert RES["split"]["pi"] == inp(".мягкий")["PR"]


def test_u2_within_material_limit():
    u_lim = inp("LIM['U_ti']")["U_ti"]
    for key in ("stage1", "stage2"):
        s = RES[key]
        assert s["U2_m_s"] <= u_lim
        assert s["U2_m_s"] <= RES["split"]["U2_design_limit"] * (1 + 1e-9)
        assert close(s["material"]["margin_U_pct"], (1 - s["U2_m_s"] / u_lim) * 100, 1e-4)


def test_power_balance():
    pb = RES["cycle"]["power_balance"]
    assert abs(pb["residual_pct"]) < 0.5
    assert close(pb["P_turb_tables_kW"] * pb["eta_m"], pb["P_comp_cantera_kW"] + pb["P_shaft_cantera_kW"], 5e-3)
    assert close(pb["P_comp_stages_tables_kW"], pb["P_comp_cantera_kW"], 5e-3)     # две ступени = сжатие в цикле
    tc = RES["cycle"]["turbine_crosscheck_tables"]
    assert close(tc["ER_tables"], tc["ER_cantera"], 1e-4)
    assert abs(tc["EGT_tables_C"] - tc["EGT_cantera_C"]) < 2.0
    assert close(tc["ER_stagewise"], tc["ER_target"], 5e-3)


def test_one_turbine_stage_is_not_enough():
    t = RES["turbine"]
    assert t["one_stage_feasible"] is False
    assert t["one_stage_work_ratio"] > 1.0
    assert t["options"]["1"]["psi"] > t["criteria"]["psi_max"]
    sel = t["options"][str(t["selected_z"])]
    assert sel["feasible"] and sel["psi"] <= t["criteria"]["psi_max"]
    # турбина снимает работу двух ступеней компрессора и полезную мощность (с механическим КПД)
    assert close(sel["w_t_kJkg"] * RES["cycle"]["frame3"]["W_exh"], RES["cycle"]["power_balance"]["P_turb_needed_kW"], 1e-4)


def test_inlet_solver_synthetic():
    area_SYNTHETIC, c_SYNTHETIC = 0.05, 150.0
    T, p, rho = M.stat(M.E.T0, M.E.P0, c_SYNTHETIC)
    c = M.inlet(M.E.T0, M.E.P0, rho * c_SYNTHETIC * area_SYNTHETIC, area_SYNTHETIC)[0]
    assert abs(c - c_SYNTHETIC) < 1e-4


def test_repeat_gives_same_numbers(dp):
    pi1 = RES["split"]["pi1"]
    a, b = M.compressor(pi1), M.compressor(pi1)
    for k in ("s1", "s2"):
        for f in ("U2", "D2", "b2", "phi", "w0", "eta_p_model"):
            assert a[k]["act"][f] == b[k]["act"][f]
        assert a[k]["eta_p"] == b[k]["eta_p"]
    for k, js in (("s1", "stage1"), ("s2", "stage2")):
        st = dp["comp"][k]
        assert close(st["act"]["U2"], RES[js]["U2_m_s"])
        assert close(st["act"]["D2"], RES[js]["D2_m"])
        assert close(st["act"]["b2"], RES[js]["b2_m"])
        assert close(st["eta_p"], RES[js]["eta_p"])
    assert close(dp["eta_pc"], RES["compressor"]["eta_p_equiv"])
    assert dp["z_sel"] == RES["turbine"]["selected_z"]
    assert close(dp["opts"][dp["z_sel"]]["cyc"]["P_MW"] * 1e3, RES["cycle"]["frame3"]["P_shaft_kW"])


def test_offdesign_model_reproduces_design(dp):
    for k in ("s1", "s2"):
        st = dp["comp"][k]
        o = M.stage_off(st, st["act"]["T01"], st["act"]["p01"], M.G, M.N)
        assert o["ok"]
        assert close(o["PRs"], st["act"]["PRs"], 1e-6)
        assert abs(o["eta_p"] - st["eta_p"]) < 1e-6
    sc = RES["offdesign_selfcheck"]
    assert close(sc["W1"], sc["G_design"], 1e-5) and close(sc["PR"], RES["split"]["pi"], 1e-5)


def test_members_listed():
    fr = inp("frames.")
    rows = RES["members"]
    assert len(rows) == fr["n"]
    for r in rows:
        assert fr["G_min"] - 0.05 <= r["G_soft"] <= fr["G_frame"] + 0.05
        if r["arch"] != "1В-2Ц" or r["cool"] != "U":
            assert r["verdict"] == "нужен другой вариант"
        else:
            assert "natural" in r and r["verdict"]
