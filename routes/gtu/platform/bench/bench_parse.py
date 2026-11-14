#!/usr/bin/env python3
"""Разбор журналов SU2 по конфигурациям протокола CPU/GPU → сводка CSV (PROTOCOL.md, §3).

Папка прогона конфигурации (runs/<id>/, имя — поле dir или id в bench.toml) содержит:
  history*.csv     история SU2 (TABULAR_FORMAT= CSV); несколько файлов — стадии, по порядку имён;
  *.log            вывод SU2_CFD (версия, число узлов, сводка WRT_PERFORMANCE= YES, «Exit Success»),
                   сюда же могут попасть строки /usr/bin/time -v («Maximum resident set size»);
  time_v*.txt      вывод /usr/bin/time -v по процессам (необязательно);
  gpu_mem*.csv     nvidia-smi --query-gpu=timestamp,memory.used,utilization.gpu --format=csv -lms 500;
  nsys_mem*.csv    nsys stats --report cuda_gpu_mem_time_sum --format csv (обмен CPU–GPU);
  profile.toml     доли профилировщика: t_assembly_s, t_linsolve_s, method (perf, Tracy);
  machine.toml     machine_s — время машины от запуска до гашения (журнал runner.py), если известно.

Что считается: время по истории двумя способами (сумма ITER_TIME и WALL_TIME·(итерация+1)),
секунды на итерацию, итерация и время сходимости по критериям bench.toml, расход вход/выход и
дисбаланс, полная степень сжатия π и изоэнтропический КПД по средним по расходу, момент двумя
способами (ṁ·cp·ΔT₀/ω и ṁ·EulerianWork/ω), пиковая память CPU и GPU, обмен CPU–GPU, стоимость
= часы × тариф (tariffs_us-east-1.json), сверка точности с эталонной конфигурацией.

  python3 bench_parse.py --bench bench.toml --runs runs --out bench_summary.csv --compare bench_compare.csv
  python3 bench_parse.py --history history.csv --log su2.log          # одна папка → JSON на экран
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import re
import statistics
import sys
import tomllib
from pathlib import Path

HERE = Path(__file__).resolve().parent

SUMMARY_FIELDS = [
    "config_id", "role", "instance", "market", "build", "precision", "enable_cuda", "mpi_ranks", "omp_threads",
    "status", "converged", "su2_version", "n_points", "n_stages", "n_iter", "iter_conv",
    "residual_field", "res_drop_orders", "time_hist_sum_s", "time_hist_avg_s", "time_check_rel",
    "s_per_iter_median", "s_per_iter_mean", "time_to_conv_s", "wall_total_s", "preproc_s", "compute_s",
    "output_s", "linsol_iter_mean", "mdot_in_kg_s", "mdot_out_kg_s", "mdot_imbalance", "mdot_wheel_kg_s",
    "pr_tt", "tau_tt", "eta_tt", "prtt_su2", "eta_tt_su2", "torque_passage_Nm", "torque_wheel_Nm",
    "torque_turbo_wheel_Nm", "torque_methods_rel", "cpu_peak_rss_sum_GiB", "cpu_peak_rss_max_GiB",
    "gpu_peak_mem_GiB", "gpu_util_mean_pct", "htod_s", "dtoh_s", "t_assembly_s", "t_linsolve_s",
    "f_linsolve", "tariff_usd_h", "extras_usd_h", "cost_solver_usd", "machine_s", "cost_machine_usd", "notes",
]
COMPARE_FIELDS = ["config_id", "reference", "both_converged", "mdot_rel", "torque_rel", "pr_rel", "eta_abs",
                  "imbalance_abs", "iter_conv_ratio", "speedup_time_to_conv", "pass_mdot", "pass_torque",
                  "pass_pr", "pass_eta", "pass_imbalance", "pass_iter_ratio", "all_pass"]


class ParseError(RuntimeError):
    pass


# ------------------------------------------------------------------------------------ история
def read_history(path: Path, time_kind: str | None = None) -> tuple[list[str], list[dict]]:
    """CSV истории SU2: заголовки в кавычках с пробелами. Две колонки «Time(sec)» — ITER_TIME и
    WALL_TIME в порядке регистрации (COutput.cpp v8.6.0); одна — только с явным time_kind."""
    text = Path(path).read_text(encoding="utf-8", errors="replace").splitlines()
    if not text:
        raise ParseError(f"{path}: пустой файл")
    if text[0].startswith("VARIABLES"):
        raise ParseError(f"{path}: формат Tecplot — нужен TABULAR_FORMAT= CSV")
    raw = [h.strip().strip('"').strip() for h in next(csv.reader([text[0]]))]
    names: list[str] = []
    tcols = [i for i, h in enumerate(raw) if h == "Time(sec)"]
    for i, h in enumerate(raw):
        if h == "Time(sec)":
            if len(tcols) == 2:
                h = "ITER_TIME" if i == tcols[0] else "WALL_TIME"
            elif len(tcols) == 1:
                if time_kind not in ("iter", "avg"):
                    raise ParseError(f"{path}: одна колонка Time(sec) — укажите time_kind = iter или avg "
                                     "(в v8.6.0 группа WALL_TIME даёт две колонки)")
                h = "ITER_TIME" if time_kind == "iter" else "WALL_TIME"
            else:
                raise ParseError(f"{path}: {len(tcols)} колонок Time(sec)")
        if h in names:
            raise ParseError(f"{path}: повтор колонки {h}")
        names.append(h)
    rows = []
    for ln in text[1:]:
        if not ln.strip():
            continue
        vals = next(csv.reader([ln]))
        if len(vals) != len(names):
            raise ParseError(f"{path}: в строке {len(vals)} значений, в заголовке {len(names)}")
        row = {}
        for n, v in zip(names, vals):
            v = v.strip()
            try:
                row[n] = float(v)
            except ValueError:
                row[n] = float("nan")
        rows.append(row)
    return names, rows


def _f(row: dict, key: str):
    v = row.get(key)
    return None if v is None or (isinstance(v, float) and math.isnan(v)) else v


def time_profile(rows: list[dict]) -> dict:
    """Время стадии двумя способами и секунды на итерацию (стационарная часть — вторая половина)."""
    it = [r.get("Inner_Iter") for r in rows]
    out = {"n_rows": len(rows), "every_iter": True, "sum_iter_s": None, "avg_method_s": None}
    if it and all(v is not None for v in it):
        steps = {round(b - a) for a, b in zip(it, it[1:])}
        out["every_iter"] = steps <= {1}
        out["n_iter"] = int(it[-1] - it[0] + 1)
    else:
        out["n_iter"] = len(rows)
    if "ITER_TIME" in rows[0] and out["every_iter"]:
        out["sum_iter_s"] = sum(r["ITER_TIME"] for r in rows)
        half = [r["ITER_TIME"] for r in rows[len(rows) // 2:]]
        out["s_per_iter_median"] = statistics.median(half) if half else None
    if "WALL_TIME" in rows[0] and it and it[-1] is not None:
        out["avg_method_s"] = rows[-1]["WALL_TIME"] * (it[-1] + 1)
    return out


def elapsed_at(rows: list[dict], idx: int) -> float | None:
    """Секунды от старта стадии до строки idx: WALL_TIME·(Inner_Iter+1), иначе сумма ITER_TIME."""
    r = rows[idx]
    if _f(r, "WALL_TIME") is not None and _f(r, "Inner_Iter") is not None:
        return r["WALL_TIME"] * (r["Inner_Iter"] + 1)
    if "ITER_TIME" in r:
        return sum(x["ITER_TIME"] for x in rows[: idx + 1])
    return None


def performance(row: dict, gamma: float) -> dict:
    """π, τ, η по средним по расходу (MARKER_ANALYZE) и расход со знаком (вход > 0, выход < 0)."""
    out = {}
    mi, mo = _f(row, "Avg_Massflow(inlet)"), _f(row, "Avg_Massflow(outlet)")
    if mi is not None and mo is not None:
        out.update(mdot_in=mi, mdot_out=mo, imbalance=(mi + mo) / mi if mi else None,
                   mass_sign_ok=mi > 0 > mo)
    p1, p2 = _f(row, "Avg_TotalPress(inlet)"), _f(row, "Avg_TotalPress(outlet)")
    t1, t2 = _f(row, "Avg_TotalTemp(inlet)"), _f(row, "Avg_TotalTemp(outlet)")
    if p1 and p2 and t1 and t2:
        pr, tau = p2 / p1, t2 / t1
        out.update(pr=pr, tau=tau, t01=t1, t02=t2)
        if tau != 1.0 and pr > 0:
            out["eta"] = (pr ** ((gamma - 1.0) / gamma) - 1.0) / (tau - 1.0)
    for k_su2, k in (("PRTT", "prtt_su2"), ("TotTotEff", "eta_su2"), ("EulerianWork", "euler_work"),
                     ("MassFlowIn_1", "mdot_turbo_in")):
        v = _f(row, k_su2)
        if v is not None:
            out[k] = v
    return out


def convergence(rows: list[dict], conv: dict, gamma: float) -> dict:
    """Первая строка, где выполнены все критерии сходимости bench.toml."""
    field = conv["residual_field"]
    if not rows or field not in rows[0]:
        return {"iter_conv_idx": None, "res_drop": None, "note": f"нет колонки {field}"}
    res = [r[field] for r in rows]
    ref = max(res[: conv["residual_ref_rows"]])
    perf = [performance(r, gamma) for r in rows]
    w = conv["window_iters"]
    tol = conv["window_rel_change_max"]
    it = [r.get("Inner_Iter", i) for i, r in enumerate(rows)]
    hit = None
    j0 = 0
    for i in range(len(rows)):
        while it[i] - it[j0] > w:                      # начало окна [it_i − w, it_i]
            j0 += 1
        if ref - res[i] < conv["min_drop_orders"]:
            continue
        p = perf[i]
        if p.get("imbalance") is None or abs(p["imbalance"]) > conv["mass_imbalance_max"]:
            continue
        if it[i] - it[0] < w:
            continue                                   # окна ещё нет
        ok = True
        for key in ("pr", "eta"):
            vals = [perf[j].get(key) for j in range(j0, i + 1)]
            if any(v is None for v in vals):
                ok = False
                break
            if (max(vals) - min(vals)) / abs(vals[-1]) > tol:
                ok = False
                break
        if ok:
            hit = i
            break
    return {"iter_conv_idx": hit, "res_drop": ref - res[-1], "res_ref": ref}


# ------------------------------------------------------------------------------------ журналы
PERF_KEYS = {
    "Wall-clock time (hrs):": "wall_h", "Cores:": "cores", "Preproc. Time (s):": "preproc_s",
    "Compute Time (s):": "compute_s", "Output Time (s):": "output_s", "Iteration count:": "iter_count",
    "Avg. s/iter:": "avg_s_iter",
}


def parse_log(text: str) -> dict:
    out: dict = {"exit_success": "Exit Success (SU2_CFD)" in text, "bad_match": text.count("Bad match")}
    m = re.search(r'Release\s+(\d+\.\d+\.\d+)', text)
    if m:
        out["su2_version"] = m.group(1)
    m = re.search(r"(\d+)\s+grid points before partitioning", text)
    if m:
        out["n_points"] = int(m.group(1))
    for key, name in PERF_KEYS.items():
        m = re.search(re.escape(key) + r"\s+([0-9.eE+-]+)", text)
        if m:
            out[name] = float(m.group(1))
    out["gpu_ilu"] = "GPU ILU scheduling" in text
    rss = [int(v) for v in re.findall(r"Maximum resident set size \(kbytes\):\s*(\d+)", text)]
    out["rss_kb"] = rss
    return out


def parse_gpu_mem(path: Path) -> dict:
    """CSV nvidia-smi: колонка 'memory.used [MiB]' (или без единиц при nounits)."""
    with open(path, newline="", encoding="utf-8") as f:
        rd = csv.reader(f)
        head = [h.strip() for h in next(rd)]
        im = next((i for i, h in enumerate(head) if h.startswith("memory.used")), None)
        iu = next((i for i, h in enumerate(head) if h.startswith("utilization.gpu")), None)
        if im is None:
            raise ParseError(f"{path}: нет колонки memory.used")
        mem, util = [], []
        for row in rd:
            if len(row) <= im:
                continue
            v = re.sub(r"[^0-9.]", "", row[im])
            if v:
                mem.append(float(v))
            if iu is not None and len(row) > iu:
                u = re.sub(r"[^0-9.]", "", row[iu])
                if u:
                    util.append(float(u))
    return {"peak_MiB": max(mem) if mem else None, "util_mean": statistics.mean(util) if util else None}


def parse_nsys_mem(path: Path) -> dict:
    """CSV отчёта nsys cuda_gpu_mem_time_sum: 'Total Time (ns)' и 'Operation' (формат по документации
    Nsight Systems; на настоящем выводе не проверялся — тест на придуманном файле)."""
    lines = Path(path).read_text(encoding="utf-8", errors="replace").splitlines()
    start = next((i for i, ln in enumerate(lines) if "Total Time (ns)" in ln and "Operation" in ln), None)
    if start is None:
        raise ParseError(f"{path}: нет заголовка с 'Total Time (ns)' и 'Operation'")
    rd = csv.DictReader(lines[start:])
    htod = dtoh = 0.0
    for r in rd:
        op = r.get("Operation", "")
        t = float(r["Total Time (ns)"]) / 1e9
        if "Host-to-Device" in op:
            htod += t
        elif "Device-to-Host" in op:
            dtoh += t
    return {"htod_s": htod, "dtoh_s": dtoh}


# ------------------------------------------------------------------------------------ папка прогона
def parse_run(run_dir: Path, case: dict, conv: dict, time_kind: str | None = None) -> dict:
    run_dir = Path(run_dir)
    gamma = case.get("gamma", 1.4)
    rgas = case.get("gas_constant_J_kgK", 287.058)
    cp = gamma * rgas / (gamma - 1.0)
    notes: list[str] = []
    hist = sorted(run_dir.glob("history*.csv"))
    if not hist:
        return {"status": "error", "notes": f"нет history*.csv в {run_dir}"}
    stages = []
    try:
        for h in hist:
            _, rows = read_history(h, time_kind)
            if not rows:
                raise ParseError(f"{h}: нет строк")
            stages.append(rows)
    except ParseError as e:
        return {"status": "error", "notes": str(e)}
    log_text = "\n".join(p.read_text(encoding="utf-8", errors="replace") for p in sorted(run_dir.glob("*.log")))
    log = parse_log(log_text)
    for p in sorted(run_dir.glob("time_v*.txt")):
        log["rss_kb"] += parse_log(p.read_text(encoding="utf-8", errors="replace"))["rss_kb"]

    out: dict = {"n_stages": len(stages), "su2_version": log.get("su2_version"), "n_points": log.get("n_points")}
    t_sum = t_avg = 0.0
    sum_ok = avg_ok = True
    n_iter = 0
    for rows in stages:
        tp = time_profile(rows)
        n_iter += tp["n_iter"]
        if tp["sum_iter_s"] is None:
            sum_ok = False
        else:
            t_sum += tp["sum_iter_s"]
        if tp["avg_method_s"] is None:
            avg_ok = False
        else:
            t_avg += tp["avg_method_s"]
        if not tp["every_iter"]:
            notes.append("история пишется не каждую итерацию: время — только по WALL_TIME")
    last = stages[-1]
    tp_last = time_profile(last)
    out.update(n_iter=n_iter, time_hist_sum_s=t_sum if sum_ok else None, time_hist_avg_s=t_avg if avg_ok else None,
               s_per_iter_median=tp_last.get("s_per_iter_median"))
    if sum_ok and avg_ok and t_avg > 0:
        out["time_check_rel"] = abs(t_sum - t_avg) / t_avg
        if out["time_check_rel"] > 1e-3:
            notes.append("время по ITER_TIME и по WALL_TIME расходится больше 0,1 %")
    t_hist = out["time_hist_avg_s"] if avg_ok else out["time_hist_sum_s"]
    out["s_per_iter_mean"] = t_hist / n_iter if t_hist and n_iter else None

    # Расхождение (NaN в невязках или расходе) — статус «diverged».
    diverged = any(any(isinstance(v, float) and math.isnan(v) for k, v in r.items()
                       if k.startswith("rms[") or k.startswith("Avg_Massflow")) for rows in stages for r in rows)
    cv = convergence(last, conv, gamma)
    out["residual_field"] = conv["residual_field"]
    out["res_drop_orders"] = cv.get("res_drop")
    prev = sum(time_profile(r)["avg_method_s"] or time_profile(r)["sum_iter_s"] or 0.0 for r in stages[:-1])
    if cv["iter_conv_idx"] is not None:
        out["iter_conv"] = int(last[cv["iter_conv_idx"]].get("Inner_Iter", cv["iter_conv_idx"]))
        el = elapsed_at(last, cv["iter_conv_idx"])
        out["time_to_conv_s"] = prev + el if el is not None else None
    out["converged"] = cv["iter_conv_idx"] is not None and not diverged

    perf = performance(last[-1], gamma)
    om, npas = case.get("omega_rad_s"), case.get("n_passages")
    out.update(mdot_in_kg_s=perf.get("mdot_in"), mdot_out_kg_s=perf.get("mdot_out"),
               mdot_imbalance=perf.get("imbalance"), pr_tt=perf.get("pr"), tau_tt=perf.get("tau"),
               eta_tt=perf.get("eta"), prtt_su2=perf.get("prtt_su2"), eta_tt_su2=perf.get("eta_su2"))
    if perf.get("mass_sign_ok") is False:
        notes.append("знаки расхода не «вход > 0, выход < 0»: обратное течение на выходе")
    if npas and perf.get("mdot_in") is not None:
        out["mdot_wheel_kg_s"] = perf["mdot_in"] * npas
    if om and perf.get("mdot_in") is not None and perf.get("t01") is not None:
        tq = perf["mdot_in"] * cp * (perf["t02"] - perf["t01"]) / om        # адиабатные стенки: P = ṁ·Δh₀
        out["torque_passage_Nm"] = tq
        if npas:
            out["torque_wheel_Nm"] = tq * npas
        if perf.get("euler_work") is not None and perf.get("mdot_turbo_in") is not None and npas:
            tt = perf["mdot_turbo_in"] * perf["euler_work"] / om * npas
            out["torque_turbo_wheel_Nm"] = tt
            if out.get("torque_wheel_Nm"):
                out["torque_methods_rel"] = abs(tt - out["torque_wheel_Nm"]) / abs(out["torque_wheel_Nm"])

    out.update(wall_total_s=log["wall_h"] * 3600.0 if "wall_h" in log else None,
               preproc_s=log.get("preproc_s"), compute_s=log.get("compute_s"), output_s=log.get("output_s"))
    if out["wall_total_s"] is None:
        notes.append("нет сводки WRT_PERFORMANCE: полное время — по истории, без предобработки")
    lin = [r.get("Linear_Solver_Iterations") for rows in stages for r in rows if r.get("Linear_Solver_Iterations") is not None]
    out["linsol_iter_mean"] = statistics.mean(lin) if lin else None
    if log["rss_kb"]:
        out["cpu_peak_rss_sum_GiB"] = sum(log["rss_kb"]) * 1024 / 2**30
        out["cpu_peak_rss_max_GiB"] = max(log["rss_kb"]) * 1024 / 2**30
    gm = [parse_gpu_mem(p) for p in sorted(run_dir.glob("gpu_mem*.csv"))]
    if gm:
        peaks = [g["peak_MiB"] for g in gm if g["peak_MiB"] is not None]
        utils = [g["util_mean"] for g in gm if g["util_mean"] is not None]
        out["gpu_peak_mem_GiB"] = max(peaks) / 1024 if peaks else None
        out["gpu_util_mean_pct"] = statistics.mean(utils) if utils else None
    ns = [parse_nsys_mem(p) for p in sorted(run_dir.glob("nsys_mem*.csv"))]
    if ns:
        out["htod_s"] = sum(n["htod_s"] for n in ns)
        out["dtoh_s"] = sum(n["dtoh_s"] for n in ns)
    prof = run_dir / "profile.toml"
    if prof.exists():
        pr = tomllib.loads(prof.read_text(encoding="utf-8"))
        out["t_assembly_s"], out["t_linsolve_s"] = pr.get("t_assembly_s"), pr.get("t_linsolve_s")
        tot = out["compute_s"] or t_hist
        if out["t_linsolve_s"] is not None and tot:
            out["f_linsolve"] = out["t_linsolve_s"] / tot
    mach = run_dir / "machine.toml"
    if mach.exists():
        out["machine_s"] = tomllib.loads(mach.read_text(encoding="utf-8")).get("machine_s")
    if log.get("bad_match"):
        notes.append(f"периодичность: {log['bad_match']} строк «Bad match»")

    if diverged:
        out["status"] = "diverged"
    elif not out["converged"]:
        out["status"] = "not_converged"
    elif not log["exit_success"]:
        out["status"] = "incomplete"
        notes.append("в журнале нет «Exit Success»")
    else:
        out["status"] = "ok"
    out["notes"] = "; ".join(notes)
    return out


def add_cost(row: dict, tariffs: dict | None, region: str, disk_gb: float | None = None) -> None:
    """Стоимость = часы × (тариф машины + диск gp3 + публичный IPv4); время — до сходимости."""
    if not tariffs:
        return
    rec = tariffs.get("regions", {}).get(region, {}).get(row.get("instance"))
    if not rec:
        return
    p = rec["od_usd_h"] if row.get("market") == "od" else rec.get("spot_est_usd_h")
    if p is None:
        return
    row["tariff_usd_h"] = p
    ex = tariffs.get("extras", {}).get(region, {})
    if disk_gb is not None and "ebs_gp3_usd_gb_month" in ex and "public_ipv4_usd_h" in ex:
        row["extras_usd_h"] = disk_gb * ex["ebs_gp3_usd_gb_month"] / 730.0 + ex["public_ipv4_usd_h"]
    full = p + (row.get("extras_usd_h") or 0.0)
    t = row.get("time_to_conv_s") or row.get("wall_total_s")
    if t:
        row["cost_solver_usd"] = t / 3600.0 * full
    if row.get("machine_s"):
        row["cost_machine_usd"] = row["machine_s"] / 3600.0 * full


def compare(rows: list[dict], ref_id: str, tol: dict) -> list[dict]:
    ref = next((r for r in rows if r["config_id"] == ref_id), None)
    if ref is None:
        raise ParseError(f"эталонной конфигурации {ref_id} нет среди прогонов")
    out = []

    def rel(a, b):
        return abs(a - b) / abs(b) if a is not None and b not in (None, 0) else None

    def absd(a, b):
        return abs(a - b) if a is not None and b is not None else None

    for r in rows:
        if r["config_id"] == ref_id:
            continue
        c = {"config_id": r["config_id"], "reference": ref_id,
             "both_converged": bool(r.get("converged")) and bool(ref.get("converged")),
             "mdot_rel": rel(r.get("mdot_in_kg_s"), ref.get("mdot_in_kg_s")),
             "torque_rel": rel(r.get("torque_wheel_Nm") or r.get("torque_passage_Nm"),
                               ref.get("torque_wheel_Nm") or ref.get("torque_passage_Nm")),
             "pr_rel": rel(r.get("pr_tt"), ref.get("pr_tt")),
             "eta_abs": absd(r.get("eta_tt"), ref.get("eta_tt")),
             "imbalance_abs": absd(r.get("mdot_imbalance"), ref.get("mdot_imbalance"))}
        ic, icr = r.get("iter_conv"), ref.get("iter_conv")
        c["iter_conv_ratio"] = (max(ic, icr) / max(1, min(ic, icr))) if ic is not None and icr is not None else None
        t, tr = r.get("time_to_conv_s"), ref.get("time_to_conv_s")
        c["speedup_time_to_conv"] = tr / t if t and tr else None
        for key, tkey, name in (("mdot_rel", "massflow_rel", "pass_mdot"), ("torque_rel", "torque_rel", "pass_torque"),
                                ("pr_rel", "pr_rel", "pass_pr"), ("eta_abs", "eta_abs", "pass_eta"),
                                ("imbalance_abs", "imbalance_abs", "pass_imbalance"),
                                ("iter_conv_ratio", "iter_conv_ratio_max", "pass_iter_ratio")):
            c[name] = None if c[key] is None else c[key] <= tol[tkey]
        checks = [c[k] for k in ("pass_mdot", "pass_torque", "pass_pr", "pass_eta", "pass_imbalance", "pass_iter_ratio")]
        c["all_pass"] = c["both_converged"] and all(v is True for v in checks)
        out.append(c)
    return out


def fmt(v):
    if v is None:
        return ""
    if isinstance(v, float):
        return f"{v:.10g}"
    return str(v)


def write_csv(path: Path, fields: list[str], rows: list[dict]) -> None:
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fields, extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow({k: fmt(r.get(k)) for k in fields})


def run_bench(bench_path: Path, runs: Path, out: Path, compare_out: Path | None, only: list[str] | None) -> int:
    bench = tomllib.loads(bench_path.read_text(encoding="utf-8"))
    tar_path = bench_path.parent / bench["meta"]["tariffs"]
    tariffs = json.loads(tar_path.read_text(encoding="utf-8")) if tar_path.exists() else None
    rows = []
    for cfg in bench["config"]:
        if only and cfg["id"] not in only:
            continue
        d = runs / cfg.get("dir", cfg["id"])
        if not d.exists():
            print(f"[{cfg['id']}] нет папки {d} — пропуск", file=sys.stderr)
            continue
        r = parse_run(d, bench["case"], bench["convergence"], cfg.get("time_kind"))
        r.update({k: cfg.get(k) for k in ("role", "instance", "market", "build", "precision", "enable_cuda",
                                         "mpi_ranks", "omp_threads")}, config_id=cfg["id"])
        add_cost(r, tariffs, bench["meta"]["region"], bench.get("decision", {}).get("disk_gb"))
        rows.append(r)
        print(f"[{cfg['id']}] {r['status']}: итераций {r.get('n_iter')}, сходимость {r.get('iter_conv')}, "
              f"время до сходимости {fmt(r.get('time_to_conv_s'))} с")
    if not rows:
        print("нет ни одной папки прогона", file=sys.stderr)
        return 1
    write_csv(out, SUMMARY_FIELDS, rows)
    print(f"записано: {out}")
    if compare_out:
        ref = bench["meta"]["reference"]
        try:
            cmp_rows = compare(rows, ref, bench["tolerances"])
        except ParseError as e:
            print(f"сверка не выполнена: {e}", file=sys.stderr)
            return 1
        write_csv(compare_out, COMPARE_FIELDS, cmp_rows)
        print(f"записано: {compare_out}")
    return 0


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--bench", default=str(HERE / "bench.toml"))
    ap.add_argument("--runs", default="runs", help="корень папок прогонов")
    ap.add_argument("--out", default="bench_summary.csv")
    ap.add_argument("--compare", default=None, help="CSV сверки точности с эталоном")
    ap.add_argument("--only", nargs="*", help="только эти конфигурации")
    ap.add_argument("--run-dir", default=None, help="разобрать одну папку и вывести JSON")
    ap.add_argument("--time-kind", choices=["iter", "avg"], default=None,
                    help="смысл единственной колонки Time(sec) (для --run-dir)")
    a = ap.parse_args(argv)
    if a.run_dir:
        bench = tomllib.loads(Path(a.bench).read_text(encoding="utf-8"))
        r = parse_run(Path(a.run_dir), bench["case"], bench["convergence"], a.time_kind)
        print(json.dumps(r, ensure_ascii=False, indent=1, default=str))
        return 0 if r["status"] != "error" else 2
    return run_bench(Path(a.bench), Path(a.runs), Path(a.out), Path(a.compare) if a.compare else None, a.only)


if __name__ == "__main__":
    sys.exit(main())
