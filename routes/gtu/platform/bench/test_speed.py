"""Тесты замеров скорости (speed_su2.py, speed_ccx.py, speed_merge.py) на придуманных данных (SYNTHETIC), без SU2 и ccx."""
import json
import math
import os
import sys

import pytest

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import speed_ccx as C  # noqa: E402
import speed_merge as M  # noqa: E402
import speed_su2 as S  # noqa: E402


def _history(path, times):
    # SYNTHETIC: заголовок как у SU2 v8.6.0 — два столбца "Time(sec)" (время итерации и среднее)
    with open(path, "w") as f:
        f.write('"Time_Iter","Outer_Iter","Inner_Iter",    "Time(sec)"   ,    "Time(sec)"   ,    "rms[Rho]"    \n')
        avg = 0.0
        for i, t in enumerate(times):
            avg = (avg * i + t) / (i + 1)
            f.write(f"0, 0, {i}, {t}, {avg}, -5.0\n")
        f.write("0, 0, 99, 1.0\n")  # недописанная строка — пропускается


def test_parse_configs_ok_and_errors():
    c = S.parse_configs(["medium=64x1,32x2", "fine=64x1"], ["medium", "fine"], 64)
    assert c == {"medium": [(64, 1), (32, 2)], "fine": [(64, 1)]}
    with pytest.raises(ValueError):
        S.parse_configs(["medium=64x2"], ["medium"], 64)       # 128 процессов > 64 ядер
    with pytest.raises(ValueError):
        S.parse_configs(["medium=64"], ["medium"], 64)         # не NPxK
    with pytest.raises(ValueError):
        S.parse_configs(["medium=64x1"], ["medium", "fine"], 64)  # нет fine


def test_history_first_time_column_and_median_of_second_half(tmp_path):
    p = tmp_path / "history_m2.csv"
    _history(p, [5.0, 3.0, 2.0, 2.0, 1.0, 1.0, 1.0, 1.2])
    t = S.read_history_times(str(p))
    assert t == [5.0, 3.0, 2.0, 2.0, 1.0, 1.0, 1.0, 1.2]  # первая колонка, а не среднее; хвост отброшен
    st = S.iter_stats(t)
    assert st["rows"] == 8 and st["s_iter_median"] == pytest.approx(1.0)
    assert st["sum_s"] == pytest.approx(16.2)


def test_summarize_and_aggregate(tmp_path):
    runs = []
    for i, base in enumerate((0.50, 0.52)):
        d = tmp_path / f"r{i}"
        d.mkdir()
        _history(d / "history_w1.csv", [base * 0.8] * 6)
        _history(d / "history_m2.csv", [base] * 10)
        json.dump({"inputs": {"mesh_header": {"NPOIN": 100000}}, "su2": {"version": "SYNTHETIC"},
                   "stages": [{"name": "w1", "status": "finished", "seconds": 30.0},
                              {"name": "m2", "status": "finished", "seconds": 40.0}]}, open(d / "run_summary.json", "w"))
        runs.append(S.summarize_run(str(d)))
    assert runs[0]["m2"]["overhead_s"] == pytest.approx(40.0 - 10 * 0.50)
    agg = S.aggregate("medium", 16, 2, runs, 100000, 1.5)
    mu = agg["muscl"]
    assert mu["s_iter"] == pytest.approx(0.51)
    assert mu["C_us_core_node_iter"] == pytest.approx(0.51 * 16 / 100000 * 1e6)
    assert mu["iters_per_hour_machine"] == pytest.approx(2 * 3600 / 0.51)


def test_su2_point_cost_picks_best_throughput_and_scales_with_nodes():
    def cfg(np_, k, s):
        return dict(level="medium", np=np_, k=k, nodes=500000, rss_peak_gb=1.0,
                    muscl=dict(measured=True, s_iter=s, iters_per_hour_machine=k * 3600 / s, C_us_core_node_iter=s * np_ / 5e5 * 1e6),
                    first_order=dict(measured=True, s_iter=s * 0.8, iters_per_hour_machine=k * 3600 / (s * 0.8)))
    doc = dict(inputs=dict(levels=["medium"]), meshes={"medium": dict(mesh_s=360, nodes=500000)},
               configs=[cfg(64, 1, 0.2), cfg(16, 4, 0.7)])  # 18000 и 20571 итер/ч — лучше 16x4
    r = M.su2_point_cost(doc, 0.75, 60, 6000, 520, 8, {"medium": {"x": 1000000}})["medium"]
    assert r["best"] == "16x4"
    mh = 6000 * 0.7 / 4 / 3600
    assert r["machine_h_per_point"] == pytest.approx(mh)
    assert r["usd_per_point"] == pytest.approx(mh * (0.75 + 60 * 0.08 / 730 + 0.005))
    assert r["scaled"]["x"]["usd_per_point"] == pytest.approx(2 * r["usd_per_point"])
    assert r["base_machine_h"] == pytest.approx(520 * 0.16 / 3600)  # первый порядок при наибольшем NP
    assert r["mesh_machine_h"] == pytest.approx(0.1)


def test_powerfit_recovers_exponent():
    xs = [1e4, 3e4, 1e5]
    f = M.powerfit(xs, [2e-6 * x ** 1.7 for x in xs])
    assert f["b"] == pytest.approx(1.7) and f["a"] == pytest.approx(2e-6)
    assert M.powerfit([1e4], [1.0]) is None


def test_ccx_case_cost_best_of_single_and_pack():
    m = dict(nodes=100000, elements=60000, mesh_s=10,
             static=[dict(threads=64, status="finished", wall_s=60.0, maxrss_gb=4.0, max_disp_m=1.0e-3),
                     dict(threads=1, status="finished", wall_s=300.0, maxrss_gb=2.0, max_disp_m=1.0e-3)],
             pack=dict(status="finished", threads=8, copies=8, cases_per_hour_machine=120.0), freq=None)
    r = M.ccx_case_cost(dict(meshes=[m]), 0.75, 60, [100000])
    row = r["meshes"][0]
    rate = 0.75 + 60 * 0.08 / 730 + 0.005
    assert row["single"]["usd"] == pytest.approx(60 / 3600 * rate)
    assert row["pack"]["usd"] == pytest.approx(rate / 120)
    assert row["usd_best"] == pytest.approx(min(60 / 3600, 1 / 120) * rate)
    assert row["max_disp_rel_spread"] == 0.0


def test_frd_and_dat_parsers(tmp_path):
    # SYNTHETIC: фрагменты .frd (формат i3, i10, 6e12.5) и .dat ccx
    frd = tmp_path / "s.frd"
    frd.write_text("    2C   2   1\n -4  DISP        4    1\n -5  D1          1    2    1    0\n"
                   " -1         1 3.00000E-03-4.00000E-03 0.00000E+00\n"
                   " -1         2 1.00000E-03 1.00000E-03 1.00000E-03\n -3\n")
    u, n = C.frd_max_disp(str(frd))
    assert n == 2 and u == pytest.approx(5.0e-3)
    dat = tmp_path / "f.dat"
    dat.write_text("     E I G E N V A L U E   O U T P U T\n\n MODE NO    EIGENVALUE    FREQUENCY\n\n"
                   "      1   0.5101E+08   0.7142E+04   0.1137E+04   0.0000E+00\n"
                   "      2   0.7689E+08   0.8769E+04   0.1395E+04   0.0000E+00\n")
    assert C.dat_frequencies(str(dat)) == [pytest.approx(1137.0), pytest.approx(1395.0)]


def test_predict_scales_quadratically_and_needs_finished():
    assert C.predict(dict(nodes=1000, wall_s=10.0, status="finished"), 2000) == pytest.approx(40.0)
    assert C.predict(dict(nodes=1000, wall_s=10.0, status="timeout"), 2000) is None
    assert math.isclose(C.predict(dict(nodes=1000, wall_s=1.0, status="finished"), 1000, exp=1.5), 1.0)
