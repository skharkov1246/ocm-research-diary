"""Проверки модели массы агрегата (pytest): аналитика блоков, согласие с расчётами трека (CAD и прочность frame3, organic/cycles),
честность калибровки (LOO не лучше выборки), полнота источников базы, монотонность по способу изготовления."""
import json
import math
import os
import sys

import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in sys.path:
    sys.path.insert(0, HERE)
import mass_model as MM  # noqa: E402

DB = json.load(open(MM.DB_FILE))


@pytest.fixture(scope="module")
def built():
    return MM.build()


def _cal(built):
    return {k: dict(K=v["K"], sigma_log=v["loo"]["rms_log"]) for k, v in built["calibration"]["per_practice"].items()}


def test_impeller_shape_factor_matches_frame3_cad():
    """k·ρ·D2³ центробежного колеса воспроизводит массу CAD-тела frame3 (Ti-6Al-4V) - по построению k = V/D2³"""
    g = json.load(open(MM.F_FRAME3_GEOM))
    p = json.load(open(MM.F_FRAME3_POINT))
    for st, sides in (("stage1", 2), ("stage2", 1)):
        m_cad = g["stages"][st]["solid"]["mass_kg"]["ti_6al_4v"]
        m = MM.block_wheel_cf("x", p[st]["D2_m"], "Ti64", sides)["mass_kg"]
        assert m == pytest.approx(m_cad, rel=1e-6)
    assert 0.05 < MM.A["k_cf_single"] < 0.09 and 0.09 < MM.A["k_cf_double"] < 0.14


def test_disk_model_reproduces_frame3_fea_disks():
    """a·ρ·D³ + b·m_обода воспроизводит три диска турбины KV-3, рассчитанных на прочность, не хуже 6 %"""
    assert MM.DISK["rows"], "нет frame3_fea_results.json"
    for r in MM.DISK["rows"]:
        assert r["fit_kg"] == pytest.approx(r["disk_kg"], rel=0.06)
    # и полная модель ступени (лопатки + замки + диск) при U = 420 м/с близка к сумме диск + обод frame3_fea
    d = json.load(open(MM.F_FRAME3_FEA))
    for k in ("turb1", "turb2", "turb3"):
        geo = d["parts"][k]["geometry"]
        Dm = geo["r_hub_m"] + geo["r_tip_m"]
        h = geo["h_airfoil_m"]
        b, a, dk, Dl = MM._stage_parts(Dm, h, geo["chord_m"], 7916.0, "turb", 420.0)
        real = [r for r in MM.DISK["rows"] if r["part"] == k][0]
        assert (b + a + dk) == pytest.approx(real["disk_kg"] + real["rim_kg"], rel=0.30)
        assert Dl == pytest.approx(real["D_live"], rel=0.10)


def test_wall_thickness_pressure_and_minimum():
    D, mat, T = 1.0, "steel", 100.0
    # малое давление: стенка по способу изготовления
    assert MM.wall_t(D, 1.5, T, mat, "cast_sand") == pytest.approx(max(0.010, 0.015 * D))
    # большое давление, печать: тонкостенная оболочка p·D/(2σ)
    p = 200.0
    t = MM.wall_t(D, p, T, mat, "lpbf")
    assert t == pytest.approx((p - 1.013) * 1e5 * D / (2 * MM.sigma_allow(mat, T)))
    # допускаемое падает с температурой
    assert MM.sigma_allow("steel", 550) < MM.sigma_allow("steel", 300)


def test_shell_mass_analytic():
    r = MM.block_shell("x", 0.8, 1.2, 1.013, 20.0, "Ti64", "lpbf")
    t = r["detail"]["t_mm"] / 1e3
    assert r["mass_kg"] == pytest.approx(math.pi * 0.8 * 1.2 * t * MM.RHO["Ti64"] * MM.PROC["lpbf"]["k_flange"], rel=1e-3)


def test_hx_area_ntu_roundtrip():
    # Cr = 1: NTU = ε/(1-ε)
    A = MM.hx_area(eps=0.9, C_min_kW_K=10.0, U_W_m2K=100.0, Cr=1.0)
    assert A == pytest.approx(9.0 * 10e3 / 100.0)
    # Cr < 1: ε(NTU) противотока возвращает исходную ε
    Cr, U, C = 0.8, 150.0, 20.0
    A = MM.hx_area(eps=0.85, C_min_kW_K=C, U_W_m2K=U, Cr=Cr)
    ntu = A * U / (C * 1e3)
    e = (1 - math.exp(-ntu * (1 - Cr))) / (1 - Cr * math.exp(-ntu * (1 - Cr)))
    assert e == pytest.approx(0.85, rel=1e-9)


def test_hx_mass_consistent_with_organic_cycles():
    """масса сердцевины по блоку (средний k_hdr) в пределах 0,6-1,6 от массы organic/cycles при той же поверхности"""
    cy = json.load(open(MM.F_CYCLES))
    r = next(x for x in cy["results"] if x.get("key") == "R90-ic-RH-ibc-none-none")["pr135"]
    for typ, key in (("psr", "foil"), ("tpms_print", "lpbf")):
        h = r["hx"][key][0]
        m = MM.block_hx(dict(type=typ, A_m2=h["A_m2"]), 1)["mass_kg"]
        assert 0.6 < m / h["mass_kg"] < 1.6, (typ, m, h["mass_kg"])


def test_database_sources_complete():
    legend = set(DB["status_legend"])
    cal = [e for e in DB["engines"] if e["use"] == "calibration"]
    assert 5 <= len(cal) <= 10
    for e in DB["engines"]:
        assert e["mass_kg"] > 0 and e["P_kW"] > 0 and e["practice"] in MM.PRACTICE
        assert e["sources"], e["id"]
        for k, v in e["sources"].items():
            assert v.get("url") and v.get("status"), (e["id"], k)
            assert any(s in v["status"] for s in legend), (e["id"], k, v["status"])
        lo, hi = e["mass_range_kg"]
        assert lo <= e["mass_kg"] <= hi
    for x in DB["packages"] + DB["reference_only"]:
        src = x.get("source")
        assert src and src.get("url") and src.get("status"), x["id"]


def test_calibration_honest(built):
    cal = built["calibration"]
    for pr, v in cal["per_practice"].items():
        assert v["n"] >= 3
        assert 0.05 < v["K"] < 50
        # исключение по одной не может быть лучше выборки
        assert v["loo"]["rms_log"] >= v["insample"]["rms_log"] - 1e-9
    # без деления на практику хуже, чем с делением (практика нужна)
    worst = max(v["loo"]["rms_log"] for v in cal["per_practice"].values())
    assert cal["single_K_all"]["loo"]["rms_log"] > worst
    # таблица LOO на каждую машину калибровки
    assert len(cal["table"]) == sum(v["n"] for v in cal["per_practice"].values())


def test_frame3_sizing_check(built):
    c = built["checks"]["frame3_sizing"]
    for e in c["D2_err_pct"]:
        assert abs(e) < 5.0
    assert c["turb_Dm_model"] == pytest.approx(c["turb_Dm_frame3"], rel=0.10)


def test_mercury_recuperator_below_package_delta(built):
    m = built["checks"]["mercury50_recuperator"]
    for k in ("eps_0.85", "eps_0.9"):
        assert max(m[k]["mass_lo_mid_hi_kg"]) < m["package_delta_kg"]


def test_evaluate_band_and_criterion(built):
    calk = _cal(built)
    ev = MM.evaluate(MM.size_from_passport(MM.kv3_spec()), cal=calk)
    assert ev["total_lo_kg"] <= ev["total_kg"] <= ev["total_hi_kg"]
    assert ev["kpd_per_t"] == pytest.approx(ev["eta_pct"] / (ev["total_kg"] / 1000), rel=1e-3)
    assert ev["kW_per_kg"] == pytest.approx(MM.kv3_spec()["P_kW"] / ev["total_kg"], rel=1e-3)
    # строки блоков в сумме дают итог (static уже умножен на K)
    assert sum(r["mass_kg"] for r in ev["rows"]) == pytest.approx(ev["total_kg"], rel=1e-3)


def test_process_and_bearings_monotone(built):
    calk = _cal(built)
    m = {p: MM.evaluate(MM.size_from_passport(MM.kv3_spec(process=p)), cal=calk)["total_kg"] for p in ("cast_invest", "lpbf")}
    m["cast_sand"] = MM.evaluate(MM.size_from_passport(MM.kv3_spec()), cal=calk)["total_kg"]
    assert m["lpbf"] < m["cast_invest"] < m["cast_sand"]
    mag = MM.evaluate(MM.size_from_passport(MM.kv3_spec(bearings=dict(kind="magnetic", n_radial=2, n_thrust=1))), cal=calk)
    assert mag["bearings_extra_kg"] > 0 and mag["total_kg"] > m["cast_sand"] * 0.99


def test_hx_adds_mass_and_pareto(built):
    ex = built["examples"]
    names = [e["name"] for e in ex["items"]]
    assert any("рекуп" in n for n in names)
    for i in ex["pareto_industrial_idx"]:
        assert ex["items"][i]["practice"] == "industrial"
    # синтетика фронта
    pts = [dict(eta_pct=30, total_kg=3000), dict(eta_pct=35, total_kg=5000), dict(eta_pct=29, total_kg=4000), dict(eta_pct=40, total_kg=9000)]
    assert MM.pareto_front(pts) == [0, 1, 3]
    assert MM.kpd_per_tonne(30.0, 3000.0) == pytest.approx(10.0)


def test_written_calibration_is_current(built):
    """mass_calibration.json соответствует коду и базе (K и LOO совпадают со свежим расчётом)"""
    d = json.load(open(MM.OUT_FILE))
    for pr, v in built["calibration"]["per_practice"].items():
        assert d["calibration"]["per_practice"][pr]["K"] == pytest.approx(v["K"], rel=1e-9)
        assert d["calibration"]["per_practice"][pr]["loo"]["rms_log"] == pytest.approx(v["loo"]["rms_log"], rel=1e-9)
