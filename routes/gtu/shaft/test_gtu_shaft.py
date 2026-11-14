"""Тесты узла «граница расчёта, валы и муфта» (pytest).

Запуск из корня репозитория: python3 -m pytest -q routes/gtu/shaft/test_gtu_shaft.py
Около 1 мин: повторный прогон обоих конфигов (без оценок чувствительности) во временную папку.
"""
import csv, hashlib, json, math, os, sys, tomllib

import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import gtu_shaft as S  # noqa: E402

GS, MD = (os.path.join(HERE, c) for c in S.CONFIGS)
RES = os.path.join(HERE, "shaft_results.json")


@pytest.fixture(scope="module")
def res():
    return json.load(open(RES))


@pytest.fixture(scope="module")
def rerun(tmp_path_factory):
    """повтор: тот же расчёт таблиц заново, во временную папку"""
    out = tmp_path_factory.mktemp("rerun")
    r = S.run([GS, MD], str(out), with_sens=False)
    return r, out


def rows(res, key):
    return res["units"][key]["rows"]


def iso_row(rs):
    return next(r for r in rs if abs(r["T_amb_C"] - 15.0) < 1e-9 and r["load_pct"] >= 100.0)


# ---------------------------------------------------------------- тождества на муфте
@pytest.mark.parametrize("key", ["taurus60_gs", "taurus60_md"])
def test_torque_times_omega_is_power(res, key):
    for r in rows(res, key):
        w = 2 * math.pi * r["n_rpm"] / 60
        assert r["M_Nm"] * w == pytest.approx(r["P_mech_kW"] * 1e3, rel=1e-12)


@pytest.mark.parametrize("key", ["taurus60_gs", "taurus60_md"])
def test_torque_in_csv_consistent(res, key):
    with open(os.path.join(HERE, f"shaft_{key}.csv"), encoding="utf-8") as fh:
        for r in csv.DictReader(fh):
            w = 2 * math.pi * float(r["n_rpm"]) / 60
            assert float(r["M_Nm"]) * w == pytest.approx(float(r["P_mech_kW"]) * 1e3, rel=1e-4)


def test_electric_power_identity(res):
    for r in rows(res, "taurus60_gs"):
        assert r["P_el_kW"] == pytest.approx(r["P_mech_kW"] * r["eta_red"] * r["eta_gen"], rel=1e-12)
        assert r["P_el_kW"] < r["P_mech_kW"]


def test_gen_train_reference_point():
    u = S.Unit(GS)
    P_el, er, eg = u.train.el(u.train.P_mech_ref)
    assert P_el == pytest.approx(5670e3, rel=1e-9)
    assert er == pytest.approx(0.985, abs=1e-9) and eg == pytest.approx(0.972, abs=1e-9)


def test_md_has_no_electric_power(res):
    assert all(r["P_el_kW"] is None for r in rows(res, "taurus60_md"))


# ---------------------------------------------------------------- ISO и паспорт
def test_iso_matches_passport_gs(res):
    r = iso_row(rows(res, "taurus60_gs"))
    pp = res["units"]["taurus60_gs"]["passport"]["values"]
    assert r["P_el_kW"] == pytest.approx(pp["P_out_kW"], rel=1e-3)          # калибровка
    assert r["eta_el_pct"] == pytest.approx(pp["eff_pct"], abs=0.05)         # калибровка
    assert r["n_rpm"] == pp["N_out_rpm"] == 14950.0
    assert abs(r["T_exh_C"] - pp["exh_temp_C"]) < 20.0                       # проверка, не калибровка
    assert r["W_exh_kg_s"] == pytest.approx(pp["exh_flow_kg_s"], rel=0.01)


def test_iso_matches_passport_md(res):
    r = iso_row(rows(res, "taurus60_md"))
    pp = res["units"]["taurus60_md"]["passport"]["values"]
    assert r["P_mech_kW"] == pytest.approx(pp["P_out_kW"], rel=1e-3)
    assert r["eta_mech_pct"] == pytest.approx(pp["eff_pct"], abs=0.05)
    assert abs(r["T_exh_C"] - pp["exh_temp_C"]) < 20.0
    assert r["W_exh_kg_s"] == pytest.approx(pp["exh_flow_kg_s"], rel=0.01)


def test_md_calibration_close_to_gtu_study(res):
    """та же калибровка, что gtu_study.calib (там допуски грубее: d 2e-4, T41 0,3 К)"""
    st = json.load(open(os.path.join(S.GTU, "gtu_study_results.json")))["T60"]
    c = res["units"]["taurus60_md"]["calibration"]
    assert c["d"] == pytest.approx(st["d"], abs=1e-3)
    assert c["T41_C"] == pytest.approx(st["T41_C"], abs=2.0)


def test_cycle_identity_with_gtu_cycle():
    ic = S.identity_check()
    for k, v in ic["rel_diff"].items():
        assert abs(v) < 1e-5, k


# ---------------------------------------------------------------- поведение вне ISO
@pytest.mark.parametrize("key", ["taurus60_gs", "taurus60_md"])
def test_power_falls_with_ambient_temperature(res, key):
    full = sorted((r for r in rows(res, key) if r["load_pct"] >= 100.0), key=lambda r: r["T_amb_C"])
    assert len(full) == 6
    for a, b in zip(full, full[1:]):
        assert b["P_engine_max_kW"] < a["P_engine_max_kW"]
        assert b["P_mech_kW"] <= a["P_mech_kW"] * (1 + 1e-9)


@pytest.mark.parametrize("key", ["taurus60_gs", "taurus60_md"])
def test_part_load_fractions(res, key):
    rs = rows(res, key)
    for Tc in {r["T_amb_C"] for r in rs}:
        full = next(r for r in rs if r["T_amb_C"] == Tc and r["load_pct"] >= 100.0)
        for r in (r for r in rs if r["T_amb_C"] == Tc and r["load_pct"] < 100.0):
            assert r["P_mech_kW"] == pytest.approx(full["P_mech_kW"] * r["load_pct"] / 100, rel=1e-5)
            assert r["T41_C"] < full["T41_C"]


def test_gs_speed_constant(res):
    assert {r["n_rpm"] for r in rows(res, "taurus60_gs")} == {14950.0}


def test_md_speed_limits(res):
    pp = res["units"]["taurus60_md"]["passport"]["values"]
    for r in rows(res, "taurus60_md"):
        assert r["n_rpm"] <= pp["N_pt_max_rpm"] + 1e-6
        assert r["n_gg_rpm"] <= pp["N_gg_max_rpm"] * (1 + 1e-6)
        assert r["surge_margin"] > 0


def test_balances(res):
    for key in ("taurus60_gs", "taurus60_md"):
        c = res["units"][key]["checks"]
        assert c["max_abs_mass_resid_kg_s"] < 1e-9
        assert c["max_abs_energy_resid_pct_fuel"] < 1e-3


def test_generator_cap_limits_power():
    u = S.Unit(GS, {"train": {"gen_rated_kWe": {"value": 6500.0, "example": True}}})
    fr, gg, pt, e, lim, P_eng = u.gs_full(243.15)
    assert lim == "генератор"
    assert u.train.el(pt["P_shaft"])[0] == pytest.approx(6500e3, rel=1e-6)
    assert P_eng > pt["P_shaft"]


def test_map_scaling_reproduces_design_point():
    u = S.Unit(GS)
    Wc, PR, eta = u.map(1.0, u.map.R_d)
    assert Wc == pytest.approx(u.d["W1"], rel=1e-12)
    assert PR == pytest.approx(u.d["PR"], rel=1e-12)
    assert eta == pytest.approx(u.d["eta_pc"], rel=1e-12)


# ---------------------------------------------------------------- повтор
def test_repeat_gives_same_numbers(res, rerun):
    r2, out = rerun
    for key in ("taurus60_gs", "taurus60_md"):
        a, b = rows(res, key), r2["units"][key]["rows"]
        assert len(a) == len(b) == 24
        for x, y in zip(a, b):
            for k, v in x.items():
                if isinstance(v, float) and k not in ("mass_resid_kg_s", "energy_resid_pct_fuel"):
                    assert y[k] == pytest.approx(v, rel=1e-7, abs=1e-9), (key, x["regime"], k)
                elif k not in ("mass_resid_kg_s", "energy_resid_pct_fuel"):
                    assert y[k] == v, (key, x["regime"], k)
        assert open(out / f"shaft_{key}.csv", encoding="utf-8").read() == \
            open(os.path.join(HERE, f"shaft_{key}.csv"), encoding="utf-8").read()


def test_results_up_to_date(res):
    """results-JSON посчитан на текущих входах: sha256 совпадают"""
    for f in res["inputs"]:
        p = os.path.normpath(os.path.join(HERE, f["file"]))
        assert hashlib.sha256(open(p, "rb").read()).hexdigest() == f["sha256"], f["file"]


def test_results_metadata(res):
    for k in ("method", "software", "date", "inputs"):
        assert res[k]
    assert len(res["date"]) == 10 and res["date"][4] == "-"
    for key in ("taurus60_gs", "taurus60_md"):
        u = res["units"][key]
        for blk in ("validation", "sensitivity"):
            assert u[blk]["method"]


# ---------------------------------------------------------------- проверка вне калибровки (регресс)
def test_off_iso_validation(res):
    g = res["units"]["taurus60_gs"]["validation"]
    m = res["units"]["taurus60_md"]["validation"]
    assert g["off_iso_max_abs_err_pct"] < 5.0
    assert m["off_iso_max_abs_err_pct"] < 2.0
    # точка 15 °C графика — проверка оцифровки: совпадает с паспортом в пределах 0,5 %
    for v in (g, m):
        p15 = next(p for p in v["points"] if p["T_amb_C"] == 15.0)
        assert abs(p15["err_pct"]) < 0.5


# ---------------------------------------------------------------- тестовые значения помечены
def test_example_values_marked():
    for cfg in (GS, MD):
        t = tomllib.load(open(cfg, "rb"))
        n = 0
        for sec in ("engine", "train", "site"):
            for k, v in t.get(sec, {}).items():
                if isinstance(v, dict):
                    assert v.get("example") is True, (cfg, k)
                    for f in ("value", "unit", "owner", "impact", "replaced_by"):
                        assert v.get(f) not in (None, ""), (cfg, k, f)
                    n += 1
        assert n >= 5
        for it in t.get("data_needed", []):
            for f in ("parameter", "unit", "owner", "impact", "replaced_by"):
                assert it.get(f), (cfg, f)


def test_data_needed_csv():
    with open(os.path.join(HERE, "data_needed.csv"), encoding="utf-8") as fh:
        rs = list(csv.reader(fh))
    assert rs[0] == ["исполнение", "параметр", "единицы СИ", "владелец", "влияние", "чем заменено сейчас"]
    params = {r[1] for r in rs[1:]}
    for cfg in (GS, MD):
        t = tomllib.load(open(cfg, "rb"))
        for sec in ("engine", "train", "site"):
            for k, v in t.get(sec, {}).items():
                if isinstance(v, dict) and v.get("example"):
                    assert f"{sec}.{k}" in params
    assert all(all(c.strip() for c in r) for r in rs[1:])
