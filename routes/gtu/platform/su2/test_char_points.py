"""Тесты char_points.py и make_jobs_su2_kv3.py без SU2 и без AWS. Все корпуса придуманы (SYNTHETIC): числа похожи
на сетки узла geom, но не взяты из расчётов."""
import csv
import json
import math
import os
import sys

import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, "..", "aws"))
sys.dont_write_bytecode = True

import char_points as C  # noqa: E402
import make_jobs_su2_kv3 as M  # noqa: E402

D1_SYNTHETIC = dict(T01=288.15, p01=101325.0, w=188000.0, PR=4.5, U2=544.0, psi=0.64, r2=0.347, m=25.5, m_side=12.75,
                    b2=0.0203, rpm=14950.0, double_sided=True)


def mesh_json(double_sided=True):
    d1 = dict(D1_SYNTHETIC, double_sided=double_sided)
    if not double_sided:
        d1.pop("m_side")
    return {"d1": d1, "contours": {"r_out": 1.3 * d1["r2"]}}


def test_parse_fracs():
    assert C.parse_fracs("0.5, 0.7,1.1") == [0.5, 0.7, 1.1]
    for bad in ("", "0.1", "0.5,0.5", "3"):
        with pytest.raises(ValueError):
            C.parse_fracs(bad)


def test_estimate_double_sided_uses_side_flow():
    e2 = C.p_out_estimate(mesh_json(True))
    e1 = C.p_out_estimate(mesh_json(False))
    assert e2["double_sided"] and e2["m_for_estimate"] == D1_SYNTHETIC["m_side"]
    # по полному расходу через ширину одной стороны скорость выше, давление ниже — так считает run_case.py
    assert e2["p_out_estimate"] > e2["p_out_estimate_run_case"]
    assert math.isclose(e1["p_out_estimate"], e1["p_out_estimate_run_case"])
    assert math.isclose(e1["p_out_estimate"], e2["p_out_estimate_run_case"])


def test_time_budget():
    assert C.time_budget(None, 123.0, 1e9, 300) == 123.0
    assert C.time_budget(360, None, 200, 300) == 360 * 60 - 480 - 200 - 300
    with pytest.raises(ValueError):
        C.time_budget(30, None, 1500, 300)
    with pytest.raises(ValueError):
        C.time_budget(None, None, 0, 300)


def test_read_uptime(tmp_path):
    p = tmp_path / "uptime"
    p.write_text("1234.56 999.0\n")
    assert C.read_uptime(str(p)) == 1234.56


def write_history(path, n, m_in=0.7, drift=0.0):
    head = ["Inner_Iter", "rms[Rho]", "Avg_Massflow(inlet)", "Avg_Massflow(outlet)", "Avg_TotalPress(inlet)",
            "Avg_TotalPress(outlet)", "Avg_TotalTemp(inlet)", "Avg_TotalTemp(outlet)"]
    with open(path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow([f'"{h}"' for h in head])
        for i in range(n):
            k = 1 + drift * i / n
            w.writerow([i, -5 - i / 1000, m_in, -m_in, 101325.0, 101325.0 * 3.2 * k, 288.15, 288.15 * 1.55])


def test_history_windows(tmp_path):
    p = tmp_path / "h.csv"
    write_history(p, 700)
    w = C.history_windows(str(p), 1.4)
    assert w["rows"] == 700 and w["converged_pi_eta"]
    write_history(p, 700, drift=0.05)
    assert not C.history_windows(str(p), 1.4)["converged_pi_eta"]
    write_history(p, 200)
    assert C.history_windows(str(p), 1.4)["rel_change"]["PR_tt"] is None


def summary_for(tmp_path, n=700, status="finished", m_in=0.7, imbalance=0.01):
    write_history(tmp_path / "history_c1_point.csv", n, m_in)
    perf = {"m_in_passage": m_in, "m_out_passage": -m_in * (1 - imbalance), "mass_imbalance": imbalance,
            "mass_sign_ok": True, "m_wheel": m_in * 9, "m_wheel_corrected": m_in * 9, "PR_tt": 3.2, "TR_tt": 1.55,
            "eta_tt_isentropic": 0.7}
    return {"inputs": {"n_passages": 9, "gamma": 1.4},
            "stages": [{"name": "c1_point", "status": status, "seconds": 100.0, "history_file": "history_c1_point.csv",
                        "history": {"rows": n, "residuals": {"rms[Rho]": {"first": -5, "last": -5.7, "max": -5}},
                                    "performance": perf},
                        "log": {"nan": False, "periodic_matched_points": 100, "bad_match": 0}}]}


def test_point_record_double_sided(tmp_path):
    r = C.point_record(summary_for(tmp_path), str(tmp_path), 0.9, 300000.0, True)
    assert math.isclose(r["performance"]["m_wheel_total"], 0.7 * 9 * 2)
    assert r["usable"] and r["converged"] and math.isclose(r["sec_per_iter"], 100 / 700)
    r1 = C.point_record(summary_for(tmp_path), str(tmp_path), 0.9, 300000.0, False)
    assert math.isclose(r1["performance"]["m_wheel_total"], 0.7 * 9)
    bad = C.point_record(summary_for(tmp_path, imbalance=0.2), str(tmp_path), 0.9, 3e5, False)
    assert not bad["usable"] and not bad["converged"]
    div = C.point_record(summary_for(tmp_path, status="diverged"), str(tmp_path), 0.9, 3e5, False)
    assert not div["usable"]


def char_summary(job, pfracs, points, label="kv3_stage1", nodes=1000):
    return {"label": label, "case": "custom", "level": "medium", "job": job, "pfracs_assigned": pfracs,
            "mesh": {"header": {"NPOIN": nodes, "NELEM": 5 * nodes}, "double_sided": True},
            "estimate": {"p_out_estimate": 3e5}, "d1_reference": {}, "settings": {}, "plans": {},
            "base": {"restart_ok": True, "statuses": ["finished"] * 3}, "verdict": "ok",
            "points": [{"pfrac": p, "status": "finished", "usable": True, "converged": p < 1.0} for p in points]}


def test_merge(tmp_path):
    a = tmp_path / "a" / "char_summary.json"
    b = tmp_path / "b" / "char_summary.json"
    a.parent.mkdir()
    b.parent.mkdir()
    a.write_text(json.dumps(char_summary("ja", [0.5, 0.9], [0.5, 0.9])))
    b.write_text(json.dumps(char_summary("jb", [0.7, 1.1], [0.7])))
    res = C.merge([str(a), str(b)], "run-SYNTHETIC")
    g = res["lines"]["kv3_stage1/medium"]
    assert [p["pfrac"] for p in g["points"]] == [0.5, 0.7, 0.9]
    assert g["missing"] == [1.1] and g["usable"] == 3 and g["converged"] == 3
    b.write_text(json.dumps(char_summary("jb", [0.9], [0.9])))
    with pytest.raises(ValueError, match="дважды"):
        C.merge([str(a), str(b)])
    b.write_text(json.dumps(char_summary("jb", [0.7], [0.8])))
    with pytest.raises(ValueError, match="не задавалась"):
        C.merge([str(a), str(b)])
    b.write_text(json.dumps(char_summary("jb", [0.7], [0.7], nodes=2000)))
    with pytest.raises(ValueError, match="разные сетки"):
        C.merge([str(a), str(b)])


def test_merge_cli(tmp_path, capsys):
    d = tmp_path / "runs" / "job1"
    d.mkdir(parents=True)
    (d / "char_summary.json").write_text(json.dumps(char_summary("job1", [0.5], [0.5])))
    lost = tmp_path / "runs" / "_lost" / "a1" / "job1"
    lost.mkdir(parents=True)
    (lost / "char_summary.json").write_text(json.dumps(char_summary("job1", [0.5], [0.5])))
    out = tmp_path / "res.json"
    assert C.main(["merge", "--in", str(tmp_path / "runs"), "--out", str(out)]) == 0
    assert len(json.loads(out.read_text())["jobs"]) == 1     # прерванная попытка (_lost) не считается


def test_run_case_argv():
    a = C.build_parser().parse_args(["run", "--case", "custom", "--level", "medium", "--geom-out", "/g", "--pfracs",
                                     "0.7", "--work", "/w", "--out", "/o", "--su2-src", "/s"])
    cmd = C.run_case_argv(a, out="/w/p", np_=16, p_out=253000.4, plan="/w/plan.json", time_limit=1000.4,
                          restart=("/w/base/restart_b3.dat", 250000.0))
    assert cmd[cmd.index("--p-out") + 1] == "253000"
    assert cmd[cmd.index("--restart-state") + 1] == "1.0,250000"
    assert "--mpi-args=--use-hwthread-cpus" in cmd and "--extra=MARKER_PLOTTING= (blade, hub)" in cmd
    assert cmd[cmd.index("--set") + 1] == "OUTPUT_FILES=RESTART, SURFACE_PARAVIEW"
    assert cmd[cmd.index("--shroud") + 1] == "slip" and cmd[cmd.index("--vaneless") + 1] == "slip"


def test_keep_files(tmp_path):
    src, dst = tmp_path / "s", tmp_path / "d"
    src.mkdir()
    for n in ("run_summary.json", "history_c1_point.csv", "su2_c1_point.log", "case_c1_point.cfg", "surface.vtu",
              "mesh_vaneless_split.su2", "restart_in.dat", "restart_c1_point.dat", "vol_c1_point.vtu"):
        (src / n).write_text("x")
    got = C.keep_files(str(src), str(dst), "restart_c1_point.dat")
    assert sorted(got) == sorted(["run_summary.json", "history_c1_point.csv", "su2_c1_point.log", "case_c1_point.cfg",
                                  "surface.vtu", "restart_c1_point.dat"])
    C.drop_big(str(src))
    assert not (src / "mesh_vaneless_split.su2").exists() and not (src / "restart_in.dat").exists()


def test_jobs_generated_load_in_runner(tmp_path):
    import runner
    for stage in ("medium", "fine"):
        doc = M.build(stage, [0.6, 0.8, 0.9, 1.0])
        p = tmp_path / f"{stage}.json"
        p.write_text(json.dumps(doc, ensure_ascii=False))
        jobs = runner.load_jobs(p)
        assigned = []
        for j in jobs:
            assert j.mode == "container" and j.image == "al2023" and j.vcpu == 64
            assert f"--job-max-min {j.max_min}" in j.command       # срок обёртки = срок сторожа задачи
            ud = runner.build_userdata(j, run_id="r-SYNTHETIC", attempt=1, repo=runner.DEFAULT_REPO,
                                       branch=M.BRANCH, bucket="bucket-EXAMPLE", key_prefix="PREFIX/gtu/r/j")
            assert len(ud.encode()) < runner.USERDATA_LIMIT
            pf = j.command.split("--pfracs ")[1].split()[0]
            line = j.command.split("--label ")[1].split()[0]
            assigned += [(line, x) for x in pf.split(",")]
        assert len(assigned) == len(set(assigned))
        for row in doc["estimate"]["jobs"]:
            assert row["total_min"] <= row["max_min"] <= 720 and row["mem_gb_points"] < M.MEM_GB


def test_medium_split_interleaved():
    jobs = M.plan_jobs("medium", [])
    kv3 = [j for j in jobs if j[0] == "kv3_stage1"]
    assert kv3[0][2] == [0.5, 0.7, 0.9, 1.1] and kv3[1][2] == [0.6, 0.8, 1.0, 1.2]
    assert all(len(j[2]) * j[3] <= 64 for j in jobs)
