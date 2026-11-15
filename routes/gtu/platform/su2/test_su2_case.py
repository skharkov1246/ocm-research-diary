"""Тесты su2_case.py: точка, оценка давления на выходе, сборка конфигурации, план стадий, разбор истории,
сравнение с точкой, план AWS (только файл, AWS не вызывается).

Все точки, сетки и истории здесь — придуманные (SYNTHETIC), не вывод SU2 и не числа трека.
Запуск SU2 (SU2_CFD -d) — только если SU2_CFD найден ($SU2_RUN или PATH), иначе пропуск.
Сверка ключей с Common/src/CConfig.cpp — только если задан $SU2_HOME (исходник SU2 v8.6.0).
Запуск: python3 -m pytest -q test_su2_case.py
"""
import csv
import gzip
import json
import math
import os
import sys

import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import su2_case as S  # noqa: E402
import run_case as RC  # noqa: E402

# придуманная точка двусторонней ступени (SYNTHETIC — не числа трека)
POINT_SYNTHETIC = {
    "label": "SYNTHETIC ступень", "inlet": "двусторонний", "PR": 4.0, "m_total_kg_s": 20.0, "m_per_side_kg_s": 10.0,
    "T01_C": 15.0, "p01_kPa": 101.325, "T02_C": 200.0, "n_rpm": 15000.0, "D2_m": 0.7, "U2_m_s": 549.78,
    "c2_m_s": 360.0, "alpha2_flow_from_radial_deg": 72.0, "blades_full": 9, "eta_is": 0.80, "eta_p": 0.83,
}


def write_point(tmp_path, point=None, key="stageX"):
    p = tmp_path / "point_SYNTHETIC.json"
    p.write_text(json.dumps({key: point or POINT_SYNTHETIC}, ensure_ascii=False), encoding="utf-8")
    return f"{p}:{key}"


def mesh_json_synthetic(omega=15000 * 2 * math.pi / 60, double=True):
    return {"d1": {"p01": 101325.0, "T01": 288.15, "r2": 0.35, "b2": 0.02, "double_sided": double},
            "contours": {"r_out": 0.455},
            "su2": {"markers": ["inlet", "outlet", "hub", "shroud", "blade", "periodic_1", "periodic_2"],
                    "periodic": {"master": "periodic_1", "slave": "periodic_2", "axis": [0, 0, 1],
                                 "center": [0, 0, 0], "angle_deg": 40.0},
                    "rotation_rate_rad_s": [0.0, 0.0, omega], "mass_flow_passage_kg_s": 10.0 / 9}}


MESH_SYNTHETIC = """NDIME= 3
NELEM= 1
10 0 1 2 3 0
NPOIN= 5
0.10 0.0 0.0 0
0.20 0.0 0.0 1
0.10 0.1 0.0 2
0.40 0.0 0.01 3
0.42 0.0 0.0 4
NMARK= 7
MARKER_TAG= inlet
MARKER_ELEMS= 1
5 0 1 2
MARKER_TAG= outlet
MARKER_ELEMS= 1
5 3 4 1
MARKER_TAG= hub
MARKER_ELEMS= 2
5 0 1 2
5 3 4 3
MARKER_TAG= shroud
MARKER_ELEMS= 2
5 0 2 3
5 3 4 4
MARKER_TAG= blade
MARKER_ELEMS= 1
5 1 2 3
MARKER_TAG= periodic_1
MARKER_ELEMS= 1
5 0 1 3
MARKER_TAG= periodic_2
MARKER_ELEMS= 1
5 0 2 3
"""


def write_mesh(tmp_path, gz=False, mj=None):
    d = tmp_path / "mesh_SYNTHETIC"
    d.mkdir(exist_ok=True)
    js = d / "impeller_custom_coarse.json"
    js.write_text(json.dumps(mj or mesh_json_synthetic()), encoding="utf-8")
    if gz:
        m = d / "passage_custom_coarse.su2.gz"
        with gzip.open(m, "wt", encoding="utf-8") as fh:
            fh.write(MESH_SYNTHETIC)
    else:
        m = d / "passage_custom_coarse.su2"
        m.write_text(MESH_SYNTHETIC, encoding="utf-8")
    return str(m), str(js)


# ---------------------------------------------------------------- точка и сетка
def test_load_point_double_sided(tmp_path):
    pt = S.load_point(write_point(tmp_path))
    assert pt["n_sides"] == 2 and pt["double_sided"]
    assert pt["m_side"] == 10.0 and pt["m_total"] == 20.0
    assert pt["p01"] == pytest.approx(101325.0) and pt["T01"] == pytest.approx(288.15)
    assert pt["omega"] == pytest.approx(15000 * 2 * math.pi / 60)
    assert pt["c_theta2"] == pytest.approx(360.0 * math.sin(math.radians(72.0)))
    assert pt["Z"] == 9 and pt["r2"] == pytest.approx(0.35)
    assert "PR" in pt["fields"] and pt["sha256"]


def test_load_point_rejects_missing_and_inconsistent(tmp_path):
    bad = dict(POINT_SYNTHETIC)
    del bad["c2_m_s"]
    with pytest.raises(SystemExit):
        S.load_point(write_point(tmp_path, bad))
    bad = dict(POINT_SYNTHETIC, m_per_side_kg_s=7.0)
    with pytest.raises(SystemExit):
        S.load_point(write_point(tmp_path, bad))


def test_consistency_ok_and_mismatch(tmp_path):
    pt = S.load_point(write_point(tmp_path))
    m, js = write_mesh(tmp_path)
    hdr = RC.read_mesh_header(m)
    c = S.consistency(pt, json.load(open(js)), hdr)
    assert c["ok"], c
    assert c["n_passages"] == 9 and c["m_passage_point_kg_s"] == pytest.approx(10.0 / 9)
    c = S.consistency(pt, mesh_json_synthetic(omega=1000.0), hdr)
    assert not c["ok"] and c["omega_rel_diff"] < -0.3
    c = S.consistency(pt, mesh_json_synthetic(double=False), hdr)
    assert not c["ok"]


def test_stage_mesh_gz_and_link(tmp_path):
    m, _ = write_mesh(tmp_path, gz=True)
    run = tmp_path / "run"
    run.mkdir()
    dst = S.stage_mesh(m, str(run))
    assert os.path.basename(dst) == "mesh.su2" and not os.path.islink(dst)
    hdr = RC.read_mesh_header(dst)
    assert hdr["NPOIN"] == 5 and len(hdr["markers"]) == 7
    m2, _ = write_mesh(tmp_path, gz=False)
    dst = S.stage_mesh(m2, str(run))
    assert open(dst).read() == MESH_SYNTHETIC


# ---------------------------------------------------------------- оценка давления на выходе
def test_estimate_p_out_back_substitution(tmp_path):
    pt = S.load_point(write_point(tmp_path))
    g, R = 1.4, 287.058
    cp, k = g * R / (g - 1), g / (g - 1)
    r_out, b_out = 0.455, 0.02
    p = S.estimate_p_out(pt, 1.0, r_out, b_out, g, R)
    # обратная подстановка: p02 по PR и T02 точки, Cθ по сохранению момента, Cm — из неразрывности при этом p
    tau = pt["T02"] / pt["T01"]
    eta = (pt["PR"] ** (1 / k) - 1) / (tau - 1)
    p02 = pt["p01"] * (1 + eta * (tau - 1)) ** k
    assert p02 / pt["p01"] == pytest.approx(pt["PR"], rel=1e-12)
    ct = pt["c_theta2"] * pt["r2"] / r_out
    T = pt["T02"] * (p / p02) ** (1 / k)
    cm = math.sqrt(2 * cp * (pt["T02"] - T) - ct ** 2)
    m = p / (R * T) * cm * 2 * math.pi * r_out * b_out
    assert m == pytest.approx(pt["m_side"], rel=1e-6)
    # частота ниже — давление ниже; при frac = 0,6 работа 0,36 от расчётной
    p06 = S.estimate_p_out(pt, 0.6, r_out, b_out, g, R)
    assert pt["p01"] < p06 < p
    with pytest.raises(ValueError):  # узкий выход — запирание
        S.estimate_p_out(pt, 1.0, r_out, 0.0005, g, R)


# ---------------------------------------------------------------- конфигурация
def _cfg(**kw):
    c = dict(name="p1", iter=100, muscl=True, restart="restart_s2.dat", omega=1570.8, p_out=250000.0,
             T01=288.15, p01=101325.0, per_master="periodic_1", per_slave="periodic_2", angle=40.0,
             shroud="slip", vaneless="slip", mesh="mesh_split.su2", numerics={}, extra=[])
    c.update(kw)
    return dict(S.build_config(c))


def test_build_config_values():
    kv = _cfg(ramp_rot=80, omega_ini=942.0, ramp_p=160, p_ini=150000.0)
    assert kv["SOLVER"] == "RANS" and kv["KIND_TURB_MODEL"] == "SST"
    assert kv["GRID_MOVEMENT"] == "ROTATING_FRAME" and kv["ROTATION_RATE"] == "0.0 0.0 1570.8"
    # периодичность: (первый, второй, центр, угол > 0 по z, сдвиг)
    per = [x.strip() for x in kv["MARKER_PERIODIC"].strip("()").split(",")]
    assert per[:2] == ["periodic_1", "periodic_2"] and float(per[7]) == 40.0 and float(per[5]) == 0.0
    # Giles: вход — p0, затем T0, направление (1, 0, 0) в осях решётки; выход — статическое давление
    g = [x.strip() for x in kv["MARKER_GILES"].strip("()").split(",")]
    assert len(g) == 18
    assert g[:7] == ["inlet", "TOTAL_CONDITIONS_PT", "101325", "288.15", "1.0", "0.0", "0.0"]
    assert g[9:12] == ["outlet", "STATIC_PRESSURE_1D", "250000"]
    assert kv["MARKER_TURBOMACHINERY"] == "(inlet, outlet)" and kv["TURBOMACHINERY_KIND"] == "AXIAL_CENTRIFUGAL"
    assert kv["MARKER_ANALYZE_AVERAGE"] == "MASSFLUX" and kv["MARKER_ANALYZE"] == "(inlet, outlet)"
    assert kv["RAMP_MOTION_FRAME"] == "YES" and kv["RAMP_MOTION_FRAME_COEFF"] == "(942, 5, 80)"
    assert kv["RAMP_OUTLET"] == "YES" and kv["RAMP_OUTLET_COEFF"] == "(150000, 5, 160)"
    assert kv["RESTART_SOL"] == "YES" and kv["SOLUTION_FILENAME"] == "restart_s2.dat"
    assert kv["MUSCL_FLOW"] == "YES" and kv["MGLEVEL"] == "0" and kv["ITER"] == "100"
    assert kv["MESH_FILENAME"] == "mesh_split.su2" and kv["RESTART_FILENAME"] == "restart_p1.dat"
    # начальное поле: статические p, T из p01, T01 при M = 0,4
    T = 288.15 / (1 + 0.2 * 0.16)
    assert float(kv["FREESTREAM_TEMPERATURE"]) == pytest.approx(T, rel=1e-9)
    assert float(kv["FREESTREAM_PRESSURE"]) == pytest.approx(101325.0 * (T / 288.15) ** 3.5, rel=1e-9)


def test_build_config_without_ramps_and_walls():
    kv = _cfg(restart=None, muscl=False)
    assert kv["RAMP_MOTION_FRAME"] == "NO" and kv["RAMP_OUTLET"] == "NO" and kv["RESTART_SOL"] == "NO"
    assert kv["MARKER_HEATFLUX"] == "(hub, 0.0, blade, 0.0)"
    assert kv["MARKER_EULER"] == "(shroud, hub_vaneless, shroud_vaneless)" and "MARKER_SHROUD" not in kv
    kv = _cfg(shroud="stationary", vaneless="as-is")
    assert kv["MARKER_HEATFLUX"] == "(hub, 0.0, blade, 0.0, shroud, 0.0)" and kv["MARKER_SHROUD"] == "(shroud)"
    assert "MARKER_EULER" not in kv
    kv = _cfg(shroud="rotating", vaneless="stationary")
    assert kv["MARKER_SHROUD"] == "(hub_vaneless, shroud_vaneless)" and "shroud, 0.0" in kv["MARKER_HEATFLUX"]
    kv = _cfg(numerics={"turb": "SA"})
    assert kv["KIND_TURB_MODEL"] == "SA" and kv["SA_OPTIONS"] == "NONE"


def test_build_config_extra_and_duplicates():
    kv = _cfg(extra=["SCREEN_WRT_FREQ_OUTER= 1"])
    assert kv["SCREEN_WRT_FREQ_OUTER"] == "1"
    with pytest.raises(ValueError):
        _cfg(extra=["ITER= 5"])
    txt = S.cfg_text(S.build_config(dict(name="x", iter=1, omega=1.0, p_out=1.0, T01=288.15, p01=101325.0,
                                         per_master="a", per_slave="b", angle=40.0, shroud="slip", vaneless="as-is",
                                         mesh="m.su2")), "заголовок\nвторая строка")
    assert txt.startswith("% заголовок\n% вторая строка\nSOLVER= RANS")


@pytest.mark.skipif(not os.environ.get("SU2_HOME") or not RC.find_cconfig(os.environ.get("SU2_HOME")),
                    reason="нет исходника SU2 ($SU2_HOME)")
def test_keys_registered_in_cconfig():
    txt = S.cfg_text(S.build_config(dict(name="x", iter=1, omega=1.0, p_out=1.0, T01=288.15, p01=101325.0,
                                         per_master="a", per_slave="b", angle=40.0, shroud="stationary",
                                         vaneless="slip", mesh="m.su2", ramp_rot=5, ramp_p=5, numerics={"turb": "SA"})))
    unknown, n = RC.check_keys(txt, RC.find_cconfig(os.environ["SU2_HOME"]))
    assert n > 700 and unknown == []


# ---------------------------------------------------------------- план
def test_make_plan_startup_and_points():
    plan = S.make_plan(0.9, [(0.75, 600), (0.85, 300)])
    assert [p["name"] for p in plan] == ["s1_ramp", "s2_ramp", "p1_0_750", "p2_0_850"]
    assert plan[0]["frac"] == pytest.approx(0.54) and plan[1]["frac"] == pytest.approx(0.9)
    assert plan[1]["pfrac"] == 0.75 and not plan[1]["muscl"] and plan[2]["muscl"]
    assert plan[2]["pramp"] == 0 and plan[3]["pramp"] == 100
    assert all(max(p.get("ramp", 0), p.get("pramp", 0)) < p["iter"] for p in plan)
    plan = S.make_plan(1.0, [(1.0, 50)], startup=False, point_pramp=100)
    assert len(plan) == 1 and plan[0]["pramp"] == 49
    plan = S.make_plan(1.0, [(0.8, 10)], iter_scale=2.0)
    assert plan[0]["iter"] == 340 and plan[0]["pramp"] == 320 and plan[2]["iter"] == 10


def test_parse_points():
    assert S.parse_points("0.75:600, 1.0:3000") == [(0.75, 600), (1.0, 3000)]
    for bad in ("", "0.75", "0.1:10", "0.8:0", "x:y"):
        with pytest.raises(SystemExit):
            S.parse_points(bad)


# ---------------------------------------------------------------- разбор истории
HEAD = ["Time_Iter", "Outer_Iter", "Inner_Iter", "Time(sec)", "Time(sec)", "rms[Rho]", "rms[RhoE]", "rms[k]",
        "Avg_Massflow(inlet)", "Avg_Massflow(outlet)", "Avg_TotalPress(inlet)", "Avg_TotalPress(outlet)",
        "Avg_TotalTemp(inlet)", "Avg_TotalTemp(outlet)", "PRTT", "TotTotEff"]


def write_history(path, n, rho0=-2.0, rho1=-7.0, m_in=1.0, m_out=-0.999, pr=4.0, tr=1.6, noise=0.0, nan_at=None):
    """Придуманная история (SYNTHETIC) в формате SU2: невязка линейно от rho0 до rho1, расход и p0, T0
    постоянны (с малым шумом noise)."""
    with open(path, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(HEAD)
        for i in range(n):
            rho = rho0 + (rho1 - rho0) * i / max(n - 1, 1)
            if nan_at is not None and i >= nan_at:
                rho = float("nan")
            s = 1 + noise * math.sin(i)
            w.writerow([0, 0, i, 1.5, 1.5 + 0.001 * i, rho, rho + 1, rho - 1, m_in, m_out, 101325.0,
                        101325.0 * pr * s, 288.15, 288.15 * tr, pr, 80.0])


def _analyse(tmp_path, **kw):
    p = tmp_path / "history_SYNTHETIC.csv"
    write_history(p, **kw)
    rows = S.read_rows(str(p))
    return S.analyse_stage(rows, 1.4, 9, 2, {}, dict(S.CONV))


def test_analyse_converged_synthetic(tmp_path):
    h = _analyse(tmp_path, n=600)
    assert h["rows"] == 600
    pf = h["performance"]
    assert pf["PR_tt"] == pytest.approx(4.0) and pf["TR_tt"] == pytest.approx(1.6)
    assert pf["eta_tt_isentropic"] == pytest.approx((4.0 ** (0.4 / 1.4) - 1) / 0.6)
    assert pf["mass_imbalance"] == pytest.approx(0.001) and pf["mass_sign_ok"]
    assert pf["m_wheel"] == pytest.approx(18.0) and pf["m_side"] == pytest.approx(9.0)
    assert pf["PRTT_su2"] == 4.0 and pf["TotTotEff_su2"] == 80.0
    assert h["sec_per_iter_median"] == pytest.approx(1.5)  # первый столбец Time(sec) — время итерации
    c = h["convergence"]
    assert c["rms_rho_drop_orders"] == pytest.approx(5.0) and c["residual_drop_ok"]
    assert c["mass_ok"] and c["stability_ok"] and c["converged"]


def test_analyse_not_converged_synthetic(tmp_path):
    c = _analyse(tmp_path, n=600, rho1=-4.0, m_out=-0.97, noise=0.01)["convergence"]
    assert not c["residual_drop_ok"] and not c["mass_ok"] and not c["stability_ok"] and not c["converged"]
    assert c["mass_imbalance"] == pytest.approx(0.03)
    c = _analyse(tmp_path, n=600, m_out=0.2)["convergence"]  # обратный ток на выходе
    assert not c["forward_flow"] and not c["converged"]
    h = _analyse(tmp_path, n=100)  # короче окна 500 итераций — устойчивость не оценивается
    assert h["stability"].get("PR_rel_change") is None and not h["convergence"]["converged"]
    assert "stability_note" in h["convergence"]
    h = _analyse(tmp_path, n=600, nan_at=550)
    assert not h["convergence"]["no_nan"] and not h["convergence"]["converged"]


def test_residual_drop_counts_from_run_maximum(tmp_path):
    p1, p2 = tmp_path / "h1_SYNTHETIC.csv", tmp_path / "h2_SYNTHETIC.csv"
    write_history(p1, 50, rho0=-1.0, rho1=-3.0)
    write_history(p2, 600, rho0=-3.0, rho1=-5.5)
    run_max = {}
    S.analyse_stage(S.read_rows(str(p1)), 1.4, 9, 1, run_max, dict(S.CONV))
    h = S.analyse_stage(S.read_rows(str(p2)), 1.4, 9, 1, run_max, dict(S.CONV))
    assert h["residuals"]["rms[Rho]"]["max_run"] == pytest.approx(-1.0)
    assert h["convergence"]["rms_rho_drop_orders"] == pytest.approx(4.5) and h["convergence"]["residual_drop_ok"]


def test_compare_with_point(tmp_path):
    pt = S.load_point(write_point(tmp_path))
    o = S.compare_with_point({"m_wheel": 18.0, "PR_tt": 3.6, "TR_tt": 1.6, "eta_tt_isentropic": 0.7}, pt)
    assert o["m_ratio"] == pytest.approx(0.9) and o["PR_ratio"] == pytest.approx(0.9)
    assert o["d_eta_points"] == pytest.approx(-0.1) and o["TR_point"] == pytest.approx(473.15 / 288.15)
    assert "не проверка" in o["note"]


# ---------------------------------------------------------------- полный проход без SU2 и с SU2 -d
def _run_args(tmp_path, extra):
    m, js = write_mesh(tmp_path)
    return ["run", "--mesh", m, "--mesh-json", js, "--point", write_point(tmp_path), "--points", "0.75:50,0.9:40",
            "--out", str(tmp_path / "run")] + extra


def test_main_without_su2_writes_chain(tmp_path, monkeypatch):
    monkeypatch.setattr(RC, "find_su2", lambda _b: None)
    assert S.main(_run_args(tmp_path, [])) == 2
    s = json.load(open(tmp_path / "run" / "case_summary.json"))
    assert [st["status"] for st in s["stages"]] == ["no_su2"] * 4 and "blocker" in s
    assert s["consistency"]["ok"] and s["mesh"]["vaneless_split_faces"] is not None
    cfg = {st["name"]: dict(ln.split("= ", 1) for ln in open(tmp_path / "run" / st["cfg"]).read().splitlines()
                            if not ln.startswith("%")) for st in s["stages"]}
    assert cfg["s1_ramp"]["RESTART_SOL"] == "NO"
    assert cfg["s2_ramp"]["SOLUTION_FILENAME"] == "restart_s1_ramp.dat"
    assert cfg["p1_0_750"]["SOLUTION_FILENAME"] == "restart_s2_ramp.dat" and cfg["p1_0_750"]["RAMP_OUTLET"] == "NO"
    assert cfg["p2_0_900"]["RAMP_OUTLET"] == "YES" and cfg["p2_0_900"]["MUSCL_FLOW"] == "YES"
    p_est = s["setup"]["p_out_estimate_design"]
    assert float(cfg["p2_0_900"]["MARKER_GILES"].split(",")[11]) == pytest.approx(0.9 * p_est, rel=1e-9)


def test_main_restart_from(tmp_path, monkeypatch):
    monkeypatch.setattr(RC, "find_su2", lambda _b: None)
    rst = tmp_path / "restart_SYNTHETIC.dat"
    rst.write_text("SYNTHETIC")
    S.main(_run_args(tmp_path, ["--restart-from", str(rst), "--restart-state", "1.0,200000"]))
    s = json.load(open(tmp_path / "run" / "case_summary.json"))
    assert [st["name"] for st in s["stages"]] == ["p1_0_750", "p2_0_900"]
    assert s["stages"][0]["restart_from"] == "restart_in.dat" and s["stages"][0]["p_out_ramp_from"] == 200000.0


@pytest.mark.skipif(RC.find_su2(None) is None, reason="SU2_CFD не найден ($SU2_RUN или PATH)")
def test_su2_dry_run_if_installed(tmp_path):
    """SU2_CFD -d разбирает конфигурацию первой стадии (на фиктивной геометрии, сетку не читает)."""
    S.main(_run_args(tmp_path, ["--dry"]))
    s = json.load(open(tmp_path / "run" / "case_summary.json"))
    assert s["stages"][0]["status"] == "dry_ok", open(tmp_path / "run" / "dry_s1_ramp.log").read()[-2000:]
    assert "v8.6.0" in (s["software"].get("su2") or "")


# ---------------------------------------------------------------- results-JSON и план AWS
def test_record_appends(tmp_path, monkeypatch):
    monkeypatch.setattr(RC, "find_su2", lambda _b: None)
    S.main(_run_args(tmp_path, []))
    s = json.load(open(tmp_path / "run" / "case_summary.json"))
    res = tmp_path / "res_SYNTHETIC.json"
    S.record(str(res), s, str(tmp_path / "run"))
    S.record(str(res), s, str(tmp_path / "run"))
    r = json.load(open(res))
    assert len(r["runs"]) == 2 and r["date"] == "2026-11-14"
    assert {"method", "software", "inputs", "title"} <= set(r)


def test_aws_plan_jobs_pass_runner_and_budget(tmp_path):
    trial = {"runs": [{"np": 4, "mesh": {"NPOIN": 250000},
                       "stages": [{"name": "s1_ramp", "rows": 100, "sec_per_iter_median": 2.0, "sec_per_iter_p10": 1.0},
                                  {"name": "s2_ramp", "rows": 300, "sec_per_iter_median": 3.0, "sec_per_iter_p10": 2.0}]}]}
    tp = tmp_path / "trial_SYNTHETIC.json"
    tp.write_text(json.dumps(trial))
    geom = {"stages": {"stage1": {"levels": {"medium": {"n_nodes": 1000000}}},
                       "stage2": {"levels": {"medium": {"n_nodes": 1100000}}}}}
    gp = tmp_path / "geom_SYNTHETIC.json"
    gp.write_text(json.dumps(geom))
    jp = tmp_path / "jobs_SYNTHETIC.json"
    assert S.main(["aws-plan", "--trial-results", str(tp), "--geom-results", str(gp), "--out", str(jp)]) == 0
    jobs = json.load(open(jp))
    assert len(jobs["jobs"]) == 10 and jobs["defaults"]["instance_type"] == "c7a.16xlarge"
    names = [j["name"] for j in jobs["jobs"]]
    assert len(set(names)) == 10 and "su2-f3-stage1-n105" in names
    for j in jobs["jobs"]:
        assert "--speed" in j["command"] and "--np \"$NPROC\"" in j["command"] and j["max_min"] <= 720
        assert "aws " not in j["command"]  # команда задачи не вызывает AWS CLI
    plan = json.load(open(tp))["aws_plan"]
    t_ref = (2.0 * 100 + 3.0 * 300) / 400
    t1 = t_ref * (1e6 / 250000) * (4 / 64) / 0.6
    assert plan["per_stage"][0]["t_iter_s_estimate"] == pytest.approx(t1)
    it = plan["iterations"]["per_line"]
    assert it == 740 + 3000 + 7 * 2000
    exp = sum(e["expected_hours_per_line"] * e["jobs"] for e in plan["per_stage"])
    assert plan["usd_expected"] == pytest.approx(exp * 0.71) and plan["within_ceiling"]
    aws = os.path.join(HERE, "..", "aws")
    sys.path.insert(0, aws)
    try:
        import runner  # noqa: PLC0415
    except ImportError as e:  # pragma: no cover
        pytest.skip(f"runner.py не импортируется: {e}")
    assert len(runner.load_jobs(str(jp))) == 10
    # замер только по выбранной стадии и по нижнему децилю
    assert S.main(["aws-plan", "--trial-results", str(tp), "--geom-results", str(gp), "--out", str(jp),
                   "--t-ref-kind", "p10", "--t-ref-stages", "s1_ramp"]) == 0
    plan = json.load(open(tp))["aws_plan"]
    assert plan["t_ref_s_per_iter"] == pytest.approx(1.0) and "s1_ramp" in plan["t_ref_source"]
    assert plan["t_ref_by_kind"]["p10"] == pytest.approx((100 * 1.0 + 300 * 2.0) / 400)
