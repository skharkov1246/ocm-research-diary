"""Тесты bench_parse.py и break_even.py на придуманных журналах и файлах (SYNTHETIC).

Все корпуса строятся здесь же, в tmp_path; ни одно число корпуса не является результатом SU2 или
ценой AWS. Запуск: python3 -m pytest -q test_bench.py
"""
from __future__ import annotations

import csv
import json
import math
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent))
import bench_parse as bp  # noqa: E402
import break_even as be  # noqa: E402

GAMMA, RGAS = 1.4, 287.058
CP = GAMMA * RGAS / (GAMMA - 1.0)
CASE = {"omega_rad_s": 1000.0, "n_passages": 9, "gamma": GAMMA, "gas_constant_J_kgK": RGAS}
CONV = {"residual_field": "rms[Rho]", "residual_ref_rows": 10, "min_drop_orders": 4.0,
        "mass_imbalance_max": 0.005, "window_iters": 500, "window_rel_change_max": 0.001}
TOL = {"massflow_rel": 1e-3, "torque_rel": 2e-3, "pr_rel": 1e-3, "eta_abs": 1e-3, "imbalance_abs": 1e-3,
       "iter_conv_ratio_max": 1.5}
P01, T01, ETA = 101325.0, 288.15, 0.8

HEADER = ["Time_Iter", "Outer_Iter", "Inner_Iter", "Time(sec)", "Time(sec)", "rms[Rho]", "rms[RhoE]",
          "Linear_Solver_Iterations", "EulerianWork", "MassFlowIn_1", "TotTotEff", "PRTT",
          "Avg_Massflow(inlet)", "Avg_Massflow(outlet)", "Avg_TotalTemp(inlet)", "Avg_TotalTemp(outlet)",
          "Avg_TotalPress(inlet)", "Avg_TotalPress(outlet)"]


def synthetic_row(i: int, t_iter: float, t_wall: float, scale_m: float = 1.0) -> list:
    """SYNTHETIC: невязка −2 − 0,01·i; дисбаланс 0,05·e^(−i/50); π = 4(1 + 0,02·e^(−i/100)); η = 0,8."""
    pr = 4.0 * (1.0 + 0.02 * math.exp(-i / 100.0))
    tau = 1.0 + (pr ** ((GAMMA - 1) / GAMMA) - 1.0) / ETA
    m_in = 1.0 * scale_m
    m_out = -(1.0 + 0.05 * math.exp(-i / 50.0)) * scale_m
    work = CP * T01 * (tau - 1.0)
    return [0, 0, i, t_iter, t_wall, -2.0 - 0.01 * i, 1.0 - 0.01 * i, 7, work, m_in, 99.0, pr,
            m_in, m_out, T01, T01 * tau, P01, P01 * pr]


def write_history(path: Path, rows: list[list], header=HEADER) -> None:
    """Как SU2: заголовки в кавычках с отступами, числа с отступами (формат CTablePrinter)."""
    with open(path, "w", encoding="utf-8") as f:
        f.write(",".join(f'   "{h}"   ' for h in header) + "\n")
        for r in rows:
            f.write(",".join(f"{v:>18.10g}" if isinstance(v, float) else f"{v:>10d}" for v in r) + "\n")


def synthetic_rows(n: int, first: float = 1.5, step: float = 0.5, start: int = 0, scale_m: float = 1.0):
    rows, total = [], 0.0
    for k in range(n):
        i = start + k
        t = first if k == 0 else step
        total += t
        rows.append(synthetic_row(i, t, total / (k + 1), scale_m))
    return rows


LOG = """-------------------------------------------------------------------------
|   / __| | | |_  )   Release 8.6.0 "Harrier"                           |
SYNTHETIC: придуманный журнал для теста, не вывод SU2
139121 grid points before partitioning.
-------------------------- Performance Summary --------------------------
Simulation totals:
   Wall-clock time (hrs):    0.125000 |            Core-hrs:     4.000000
                   Cores:          32 |          DOFs/point:            7
Preprocessing phase:
       Preproc. Time (s):       20.00 |   Preproc. Time (%):       4.4444
Compute phase:
        Compute Time (s):      410.00 |    Compute Time (%):      91.1111
        Iteration count:        1200 |          Avg. s/iter:     0.341667
Output phase:
         Output Time (s):       20.00 |     Output Time (%):       4.4444
------------------------- Exit Success (SU2_CFD) ------------------------
\tMaximum resident set size (kbytes): 1048576
\tMaximum resident set size (kbytes): 1048576
"""


def make_run(d: Path, n: int = 1200, log: str = LOG, **kw) -> Path:
    d.mkdir(parents=True, exist_ok=True)
    write_history(d / "history_SYNTHETIC.csv", synthetic_rows(n, **kw))
    (d / "su2_SYNTHETIC.log").write_text(log, encoding="utf-8")
    return d


def expected_conv_iter() -> int:
    """Независимый перебор тех же условий в замкнутом виде (SYNTHETIC)."""
    for i in range(500, 5000):
        drop = 0.01 * i >= 4.0
        imb = 0.05 * math.exp(-i / 50.0) <= 0.005
        pr = lambda k: 4.0 * (1.0 + 0.02 * math.exp(-k / 100.0))  # noqa: E731
        stable = (pr(i - 500) - pr(i)) / pr(i) <= 0.001
        if drop and imb and stable:
            return i
    raise AssertionError


# ------------------------------------------------------------------------------------ история
def test_two_time_columns_mapped(tmp_path):
    p = tmp_path / "history.csv"
    write_history(p, synthetic_rows(5))
    names, rows = bp.read_history(p)
    assert "ITER_TIME" in names and "WALL_TIME" in names
    assert rows[1]["ITER_TIME"] == pytest.approx(0.5) and rows[1]["WALL_TIME"] == pytest.approx(1.0)
    tp = bp.time_profile(rows)
    assert tp["sum_iter_s"] == pytest.approx(3.5) and tp["avg_method_s"] == pytest.approx(3.5)


def test_single_time_column_needs_kind(tmp_path):
    p = tmp_path / "history.csv"
    hdr = [h for k, h in enumerate(HEADER) if k != 4]
    write_history(p, [r[:4] + r[5:] for r in synthetic_rows(3)], header=hdr)
    with pytest.raises(bp.ParseError):
        bp.read_history(p)
    names, _ = bp.read_history(p, time_kind="iter")
    assert "ITER_TIME" in names and "WALL_TIME" not in names


def test_tecplot_rejected(tmp_path):
    p = tmp_path / "history.dat"
    p.write_text('VARIABLES = \\\n"Inner_Iter"\n', encoding="utf-8")
    with pytest.raises(bp.ParseError):
        bp.read_history(p)


def test_convergence_iteration_and_time(tmp_path):
    r = bp.parse_run(make_run(tmp_path / "run"), CASE, CONV)
    exp = expected_conv_iter()
    assert r["iter_conv"] == exp
    assert r["time_to_conv_s"] == pytest.approx(1.5 + 0.5 * exp, rel=1e-8)
    assert r["status"] == "ok" and r["converged"] is True
    assert r["time_check_rel"] < 1e-8
    assert r["n_iter"] == 1200 and r["linsol_iter_mean"] == 7


def test_not_converged_short(tmp_path):
    r = bp.parse_run(make_run(tmp_path / "run", n=450), CASE, CONV)
    assert r["status"] == "not_converged" and r.get("iter_conv") is None


def test_eta_and_torque_two_methods(tmp_path):
    r = bp.parse_run(make_run(tmp_path / "run"), CASE, CONV)
    assert r["eta_tt"] == pytest.approx(ETA, rel=1e-8)
    pr_last = 4.0 * (1.0 + 0.02 * math.exp(-1199 / 100.0))
    assert r["pr_tt"] == pytest.approx(pr_last, rel=1e-9)
    tau = 1.0 + (pr_last ** (0.4 / 1.4) - 1.0) / ETA
    tq = 1.0 * CP * T01 * (tau - 1.0) / 1000.0 * 9
    assert r["torque_wheel_Nm"] == pytest.approx(tq, rel=1e-8)
    assert r["torque_methods_rel"] < 1e-8
    assert r["mdot_wheel_kg_s"] == pytest.approx(9.0)


def test_log_and_aux_files(tmp_path):
    d = make_run(tmp_path / "run")
    (d / "gpu_mem_SYNTHETIC.csv").write_text(
        "timestamp, memory.used [MiB], utilization.gpu [%]\n"
        "2026/11/14 10:00:00.000, 1024 MiB, 50 %\n2026/11/14 10:00:00.500, 4096 MiB, 70 %\n"
        "2026/11/14 10:00:01.000, 2048 MiB, 90 %\n", encoding="utf-8")
    (d / "nsys_mem_SYNTHETIC.csv").write_text(
        "SYNTHETIC preamble line\n"
        '"Time (%)","Total Time (ns)","Count","Avg (ns)","Med (ns)","Min (ns)","Max (ns)","StdDev (ns)","Operation"\n'
        '80.0,2000000000,100,20000000,20000000,1,1,0,"[CUDA memcpy Host-to-Device]"\n'
        '20.0,500000000,100,5000000,5000000,1,1,0,"[CUDA memcpy Device-to-Host]"\n'
        '0.0,1000,1,1000,1000,1000,1000,0,"[CUDA memset]"\n', encoding="utf-8")
    (d / "profile.toml").write_text('method = "SYNTHETIC"\nt_assembly_s = 100.0\nt_linsolve_s = 307.5\n',
                                    encoding="utf-8")
    (d / "machine.toml").write_text("machine_s = 1800.0  # SYNTHETIC\n", encoding="utf-8")
    r = bp.parse_run(d, CASE, CONV)
    assert r["su2_version"] == "8.6.0" and r["n_points"] == 139121
    assert r["wall_total_s"] == pytest.approx(450.0) and r["compute_s"] == pytest.approx(410.0)
    assert r["cpu_peak_rss_sum_GiB"] == pytest.approx(2.0) and r["cpu_peak_rss_max_GiB"] == pytest.approx(1.0)
    assert r["gpu_peak_mem_GiB"] == pytest.approx(4.0) and r["gpu_util_mean_pct"] == pytest.approx(70.0)
    assert r["htod_s"] == pytest.approx(2.0) and r["dtoh_s"] == pytest.approx(0.5)
    assert r["f_linsolve"] == pytest.approx(307.5 / 410.0)
    assert r["machine_s"] == pytest.approx(1800.0)


def test_diverged_nan(tmp_path):
    d = make_run(tmp_path / "run", n=50)
    txt = (d / "history_SYNTHETIC.csv").read_text(encoding="utf-8").splitlines()
    parts = txt[-1].split(",")
    parts[5] = "nan"
    txt[-1] = ",".join(parts)
    (d / "history_SYNTHETIC.csv").write_text("\n".join(txt) + "\n", encoding="utf-8")
    assert bp.parse_run(d, CASE, CONV)["status"] == "diverged"


def test_incomplete_without_exit_success(tmp_path):
    r = bp.parse_run(make_run(tmp_path / "run", log="SYNTHETIC: журнал оборван\n"), CASE, CONV)
    assert r["status"] == "incomplete" and r["wall_total_s"] is None


def test_multistage_time(tmp_path):
    d = tmp_path / "run"
    d.mkdir()
    write_history(d / "history_1_SYNTHETIC.csv", synthetic_rows(200, first=2.0, step=1.0))
    write_history(d / "history_2_SYNTHETIC.csv", synthetic_rows(1200))
    (d / "su2_SYNTHETIC.log").write_text(LOG, encoding="utf-8")
    r = bp.parse_run(d, CASE, CONV)
    stage1 = 2.0 + 1.0 * 199
    assert r["n_stages"] == 2 and r["n_iter"] == 1400
    assert r["time_to_conv_s"] == pytest.approx(stage1 + 1.5 + 0.5 * expected_conv_iter(), rel=1e-8)


def test_history_not_every_iteration(tmp_path):
    d = tmp_path / "run"
    d.mkdir()
    rows = synthetic_rows(1200)
    write_history(d / "history_SYNTHETIC.csv", rows[::10])
    (d / "su2_SYNTHETIC.log").write_text(LOG, encoding="utf-8")
    r = bp.parse_run(d, CASE, CONV)
    assert r["time_hist_sum_s"] is None
    assert r["time_hist_avg_s"] == pytest.approx(1.5 + 0.5 * 1190, rel=1e-8)
    assert "не каждую итерацию" in r["notes"]


def test_compare_tolerances(tmp_path):
    ref = bp.parse_run(make_run(tmp_path / "a"), CASE, CONV)
    close = bp.parse_run(make_run(tmp_path / "b", scale_m=1.00001, step=0.25, first=1.0), CASE, CONV)
    far = bp.parse_run(make_run(tmp_path / "c", scale_m=1.005), CASE, CONV)
    rows = [dict(ref, config_id="ref"), dict(close, config_id="close"), dict(far, config_id="far")]
    out = {c["config_id"]: c for c in bp.compare(rows, "ref", TOL)}
    assert out["close"]["all_pass"] is True
    n = expected_conv_iter()
    assert out["close"]["speedup_time_to_conv"] == pytest.approx((1.5 + 0.5 * n) / (1.0 + 0.25 * n), rel=1e-8)
    assert out["far"]["pass_mdot"] is False and out["far"]["all_pass"] is False
    with pytest.raises(bp.ParseError):
        bp.compare(rows, "missing", TOL)


BENCH_SYNTH = """
[meta]
date = "2026-11-14"
region = "us-east-1"
fallback_region = "eu-north-1"
tariffs = "tariffs_SYNTHETIC.json"
mesh_stats = "mesh_SYNTHETIC.json"
reference = "cpu"
[case]
omega_rad_s = 1000.0
n_passages = 9
gamma = 1.4
gas_constant_J_kgK = 287.058
[convergence]
residual_field = "rms[Rho]"
residual_ref_rows = 10
min_drop_orders = 4.0
mass_imbalance_max = 0.005
window_iters = 500
window_rel_change_max = 0.001
[tolerances]
massflow_rel = 1.0e-3
torque_rel = 2.0e-3
pr_rel = 1.0e-3
eta_abs = 1.0e-3
imbalance_abs = 1.0e-3
iter_conv_ratio_max = 1.5
[decision]
speedup_margin = 1.2
repeats = 2
vram_fill_max = 0.85
ram_fill_max = 0.8
gpu_vram_bytes = { "gpu.SYNTH" = 80.0e9 }
series_time_limit_h = 1.0
strong_eff_min = 0.7
weak_eff_min = 0.8
amdahl_r_host = [1.0, 2.0]
amdahl_kernel_speedup = [10.0, 0.0]
[break_even]
gpu = ["gpu.SYNTH"]
cpu = ["cpu.SYNTH"]
markets = [["od", "od"], ["spot", "spot"]]
[memory_model]
nvar_flow = 5
nvar_turb = 2
krylov_restart = 10
extra_vectors = 4
index_bytes = 4
[gpu_specs."gpu.SYNTH"]
name = "SYNTHETIC GPU"
fp64_flops = 1.0e12
fp32_flops = 2.0e12
mem_bw_Bps = 1.0e12
[[config]]
id = "cpu"
role = "cpu"
instance = "cpu.SYNTH"
market = "spot"
[[config]]
id = "gpu"
role = "gpu"
instance = "gpu.SYNTH"
market = "od"
"""
TARIFFS_SYNTH = {"inputs": [], "missing": {"us-east-1": []}, "regions": {"us-east-1": {
    "cpu.SYNTH": {"od_usd_h": 2.0, "spot_est_usd_h": 1.0},
    "gpu.SYNTH": {"od_usd_h": 6.0, "spot_est_usd_h": 3.0}}}}


def test_cli_end_to_end_and_verdict(tmp_path):
    (tmp_path / "bench.toml").write_text(BENCH_SYNTH, encoding="utf-8")
    (tmp_path / "tariffs_SYNTHETIC.json").write_text(json.dumps(TARIFFS_SYNTH), encoding="utf-8")
    (tmp_path / "mesh_SYNTHETIC.json").write_text(json.dumps({"meshes": [
        {"file": "SYNTHETIC.su2", "n_points": 1000, "edges_per_point": 6.0}]}), encoding="utf-8")
    make_run(tmp_path / "runs" / "cpu")                             # 401 с до сходимости
    make_run(tmp_path / "runs" / "gpu", first=0.3, step=0.05)       # ≈ 40 с — ускорение ≈ 10
    rc = bp.main(["--bench", str(tmp_path / "bench.toml"), "--runs", str(tmp_path / "runs"),
                  "--out", str(tmp_path / "s.csv"), "--compare", str(tmp_path / "c.csv")])
    assert rc == 0
    rows = {r["config_id"]: r for r in csv.DictReader(open(tmp_path / "s.csv", encoding="utf-8"))}
    t_cpu = 1.5 + 0.5 * 799
    assert float(rows["cpu"]["cost_solver_usd"]) == pytest.approx(t_cpu / 3600 * 1.0, rel=1e-8)
    assert float(rows["gpu"]["tariff_usd_h"]) == 6.0
    res = be.compute(tmp_path / "bench.toml", tmp_path / "s.csv", tmp_path / "c.csv")
    v = res["measured"]["gpu"][0]
    s = t_cpu / (0.3 + 0.05 * 799)
    assert v["speedup"] == pytest.approx(s, rel=1e-3) and v["s_star"] == pytest.approx(6.0)
    assert v["pass_cost"] is True and v["pass_accuracy"] is True
    assert v["verdict"] == "не хватает измерений"                    # нет замера видеопамяти
    (tmp_path / "runs" / "gpu" / "gpu_mem_SYNTHETIC.csv").write_text(
        "timestamp, memory.used [MiB]\n2026/11/14 10:00:00.000, 10240 MiB\n", encoding="utf-8")
    bp.main(["--bench", str(tmp_path / "bench.toml"), "--runs", str(tmp_path / "runs"),
             "--out", str(tmp_path / "s.csv"), "--compare", str(tmp_path / "c.csv")])
    v = be.compute(tmp_path / "bench.toml", tmp_path / "s.csv", tmp_path / "c.csv")["measured"]["gpu"][0]
    assert v["pass_vram"] is True and v["verdict"] == "переносить"
    od = next(p for p in res["break_even"] if p["gpu_market"] == "od")
    assert od["s_star"] == pytest.approx(3.0) and od["s_required"] == pytest.approx(3.6)


# ------------------------------------------------------------------------------------ формулы порога
def test_break_even_and_amdahl_roundtrip():
    assert be.break_even_speedup(6.88, 1.49936) == pytest.approx(4.5886, abs=1e-4)
    for target in (1.5, 4.6, 12.0):
        for sk in (10.0, 50.0, 0.0):
            for r in (1.0, 2.0):
                f = be.amdahl_required_fraction(target, sk, r)
                if f is None:
                    assert be.amdahl_speedup(1.0, sk, r) < target
                elif f > 0:
                    assert be.amdahl_speedup(f, sk, r) == pytest.approx(target, rel=1e-9)
    assert be.amdahl_required_fraction(5.5, 5.0) is None          # ядро в 5 раз не даст 5,5
    assert be.amdahl_required_fraction(0.8, 10.0, 1.0) == 0.0     # окупается без переноса


def test_memory_and_intensity_by_hand():
    m = be.memory_per_point(6.0, 5, 8, 10, 4, 4)
    assert m["matrix"] == pytest.approx(13 * 25 * 8)
    assert m["index"] == pytest.approx(14 * 4)
    assert m["vectors"] == pytest.approx(26 * 5 * 8)
    assert m["total"] == pytest.approx(2 * (13 * 25 * 8 + 14 * 4) + 26 * 5 * 8)
    best, worst = be.spmv_intensity(6.0, 5, 8, 4)
    flops = 2 * 25 * 13
    base = 25 * 13 * 8 + 12 * 4 + 8 + 40
    assert best == pytest.approx(flops / (base + 40)) and worst == pytest.approx(flops / (base + 13 * 40))


def test_mesh_stats_two_tets(tmp_path):
    p = tmp_path / "two_tets_SYNTHETIC.su2"
    p.write_text("% SYNTHETIC\nNDIME= 3\nNELEM= 2\n10 0 1 2 3 0\n10 1 2 3 4 1\nNPOIN= 5\n"
                 "0 0 0 0\n1 0 0 1\n0 1 0 2\n0 0 1 3\n1 1 1 4\nNMARK= 0\n", encoding="utf-8")
    st = be.mesh_stats(p)
    assert st["n_points"] == 5 and st["n_elements"] == 2 and st["n_edges"] == 9


def test_tariffs_from_synthetic_price_list(tmp_path):
    pl = tmp_path / "index_SYNTHETIC.csv"
    cols = ["SKU", "OfferTermCode", "RateCode", "TermType", "PriceDescription", "EffectiveDate", "StartingRange",
            "EndingRange", "Unit", "PricePerUnit", "Currency", "Instance Type", "vCPU", "Memory", "GPU", "GPU Memory",
            "Operating System", "Tenancy", "Pre Installed S/W", "CapacityStatus", "License Model", "MarketOption",
            "Network Performance", "Physical Processor"]
    base = {"SKU": "SYNTH1", "OfferTermCode": "X", "RateCode": "SYNTH1.X.Y", "TermType": "OnDemand",
            "PriceDescription": "SYNTHETIC", "EffectiveDate": "2026-11-01", "StartingRange": "0", "EndingRange": "Inf",
            "Unit": "Hrs", "PricePerUnit": "2.0000000000", "Currency": "USD", "Instance Type": "c8i.8xlarge",
            "vCPU": "32", "Memory": "64 GiB", "GPU": "", "GPU Memory": "NA", "Operating System": "Linux",
            "Tenancy": "Shared", "Pre Installed S/W": "NA", "CapacityStatus": "Used",
            "License Model": "No License required", "MarketOption": "OnDemand", "Network Performance": "x",
            "Physical Processor": "x"}
    with open(pl, "w", newline="", encoding="utf-8") as f:
        f.write('"FormatVersion","v1.0"\n"Disclaimer","SYNTHETIC"\n"Publication Date","2026-11-14T00:00:00Z"\n'
                '"Version","SYNTHETIC1"\n"OfferCode","AmazonEC2"\n')
        w = csv.DictWriter(f, fieldnames=cols, quoting=csv.QUOTE_ALL)
        w.writeheader()
        w.writerow(base)
        w.writerow(dict(base, SKU="SYNTH2", **{"Operating System": "Windows", "PricePerUnit": "9.0"}))
        w.writerow(dict(base, SKU="SYNTH3", **{"CapacityStatus": "UnusedCapacityReservation", "PricePerUnit": "8.0"}))
    adv = tmp_path / "spot_SYNTHETIC.json"
    adv.write_text(json.dumps({"ranges": [{"index": 4, "label": ">20%"}],
                               "spot_advisor": {"us-east-1": {"Linux": {"c8i.8xlarge": {"s": 60, "r": 4}}}}}),
                   encoding="utf-8")

    class A:
        price_list = [f"us-east-1={pl}"]
        spot_advisor = str(adv)
        types = "c8i.8xlarge,p5.4xlarge"

    out = be.tariffs(A)
    rec = out["regions"]["us-east-1"]["c8i.8xlarge"]
    assert rec["od_usd_h"] == 2.0 and rec["spot_est_usd_h"] == pytest.approx(0.8)
    assert rec["spot_interruption"] == ">20%"
    assert out["missing"]["us-east-1"] == ["p5.4xlarge"]
