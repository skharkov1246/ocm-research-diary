"""Проверка вне калибровки (validate_fleet.py): загрузка паспортов, наборы КПД узлов, подбор T41,
воспроизводимость results-JSON и согласие CSV с JSON. Запуск: python3 -m pytest -q (около 15 с)."""
import csv
import json
import math

import pytest

import validate_fleet as V

RES = V.HERE / "validation_fleet_results.json"
CSV = V.HERE / "validation_fleet.csv"


@pytest.fixture(scope="module")
def cfg():
    return V.load_cfg()


@pytest.fixture(scope="module")
def cal(cfg):
    return V.load_calibration(cfg)


@pytest.fixture(scope="module")
def recs(cfg):
    return {r["id"]: r for r in V.load_records(cfg)}


def test_records_loaded_and_classified(recs, cfg):
    assert len(recs) == 32                          # все записи паспортов с мощностью, КПД, π, расходом и T выхлопа
    cls = [r["cls"] for r in recs.values()]
    assert cls.count("calibration") == len(cfg["records"]["calibration"]) == 2
    assert cls.count("related") == len(cfg["records"]["related"]) == 4
    assert cls.count("independent") == 26
    assert [r["id"] for r in recs.values() if r["recuperated"]] == ["mercury50:GS"]
    for r in recs.values():
        assert r["shafts"] in ("single", "two", "three"), r["id"]
        assert 4 < r["PR"] < 25 and 5 < r["W_exh"] < 120 and 300 < r["EGT_C"] < 600
        # КПД паспорта согласован с расходом тепла (3600 / HR)
        if r["HR_kJ_kWh"]:
            assert abs(3600 / r["HR_kJ_kWh"] * 100 - r["eff_pct"]) < 0.25, r["id"]


def test_calibration_sets_from_study(cal):
    st = json.load(open(V.HERE.parent / "gtu_study_results.json"))
    for k in ("T60", "SGT400"):
        for p in ("eta_pc", "eta_pt_gg", "eta_pt_pt", "x_ngv", "x_rot", "x_leak", "T41_C"):
            assert cal[k][p] == st[k][p]


def test_interpolation_endpoints_and_clamp(cal, recs):
    assert V.interp_t(cal["T60"]["W_cal"], cal) == 0.0
    assert V.interp_t(cal["SGT400"]["W_cal"], cal) == 1.0
    assert V.interp_t(5.0, cal) == 0.0 and V.interp_t(100.0, cal) == 1.0
    v = dict(kind="interp")
    rec = dict(recs["taurus60:MD"], W_exh=cal["T60"]["W_cal"])
    p = V.base_params(rec, v, cal)
    for k in ("eta_pc", "eta_pt_gg", "eta_pt_pt", "x_ngv", "x_rot"):
        assert math.isclose(p[k], cal["T60"][k])
    # охлаждение по T41: на калибровочной T41 доли не меняются, на +50 °C — +12 %
    vc = dict(kind="interp", cool_T41=True)
    kw0 = V.cycle_kwargs(p, vc, p["T41_ref_C"], rec)
    kw1 = V.cycle_kwargs(p, vc, p["T41_ref_C"] + 50, rec)
    assert math.isclose(kw0["x_ngv"], p["x_ngv"]) and math.isclose(kw1["x_ngv"], p["x_ngv"] * 1.12)


def test_fit_T41_synthetic():
    # монотонная «модель»: выхлоп = 0,5·T41 − 30, ниже 700 °C не считается
    def f(T):
        if T < 700:
            raise ValueError("ТВД не вращает компрессор")
        return None, 0.5 * T - 30

    T, _, n = V.fit_T41(f, 480.0, 600, 1500, 0.01)
    assert abs(T - 1020) < 0.02 and n < 20
    with pytest.raises(V.NoSolution):
        V.fit_T41(f, 200.0, 600, 1500, 0.01)      # ниже наименьшей достижимой
    with pytest.raises(V.NoSolution):
        V.fit_T41(f, 900.0, 600, 1500, 0.01)      # выше наибольшей


def test_generator_conversion(recs, cfg):
    rec = recs["taurus60:GS"]
    c = V.compare(rec, 300.0, 0.33, cfg)
    conv = cfg["model"]["eta_generator"] * cfg["model"]["eta_gearbox"]
    assert math.isclose(c["P_pred_kW"], 300.0 * rec["W_exh"] * conv)
    assert math.isclose(c["eff_pred_pct"], 33 * conv)
    assert math.isclose((1 + c["dP_pct"] / 100) * (1 + c["dW_pct"] / 100), 1.0)
    assert V.compare(recs["taurus60:MD"], 300.0, 0.33, cfg)["conv"] == 1.0


def test_real_fit_matches_results_file(recs, cal, cfg):
    # повтор одной точки: Taurus 60 привод, основной вариант
    v = next(x for x in cfg["variant"] if x.get("primary"))
    o = V.task_cycle((recs["taurus60:MD"], v, cal, cfg))
    assert o["status"] == "ok"
    assert abs(o["EGT_model_C"] - recs["taurus60:MD"]["EGT_C"]) < 0.2
    if RES.exists():
        ref = next(r for r in json.load(open(RES))["records"] if r["id"] == "taurus60:MD")["results"][v["key"]]
        assert abs(o["T41_C"] - ref["T41_C"]) < 0.05
        assert abs(o["dP_pct"] - ref["dP_pct"]) < 0.01 and abs(o["deta_pp"] - ref["deta_pp"]) < 0.005


@pytest.mark.skipif(not RES.exists(), reason="нет validation_fleet_results.json — сначала python3 validate_fleet.py")
def test_results_file_complete_and_csv_consistent(cfg):
    res = json.load(open(RES))
    for k in ("method", "software", "date", "inputs"):
        assert res[k]
    assert res["date"] == "2026-11-14"
    P = res["primary_variant"]
    assert P == "interp"
    keys = {v["key"] for v in cfg["variant"]}
    n_results = 0
    for rec in res["records"]:
        assert keys <= set(rec["results"]), rec["id"]
        for vk, row in rec["results"].items():
            n_results += 1
            assert row["status"] in ("ok", "нет решения")
            if row["status"] == "ok" and not vk.startswith("recup_"):
                assert abs(row["EGT_model_C"] - rec["EGT_C"]) < 0.5, (rec["id"], vk)
                assert abs((1 + row["dP_pct"] / 100) * (1 + row["dW_pct"] / 100) - 1) < 1e-9
    rows = list(csv.DictReader(open(CSV)))
    assert len(rows) == n_results
    by = {(r["id"], r["model_variant"]): r for r in rows}
    for rec in res["records"]:
        row = rec["results"][P]
        c = by[(rec["id"], P)]
        if row["status"] == "ok":
            assert abs(float(c["dP_pct"]) - row["dP_pct"]) < 0.006
            assert abs(float(c["deta_pp"]) - row["deta_pp"]) < 0.006
            assert abs(float(c["T41_fit_C"]) - row["T41_C"]) < 0.06
    # статистика пересчитывается из записей
    ind = [r["results"][P] for r in res["records"] if r["cls"] == "independent" and not r["recuperated"]
           and r["results"][P]["status"] == "ok"]
    s = res["summary"][P]["independent_no_recup"]
    assert s["n"] == len(ind)
    assert math.isclose(s["dP_pct"]["mean_abs"], sum(abs(r["dP_pct"]) for r in ind) / len(ind), rel_tol=1e-9)
    assert math.isclose(s["deta_pp"]["max_abs"], max(abs(r["deta_pp"]) for r in ind), rel_tol=1e-9)
