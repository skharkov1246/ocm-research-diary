#!/usr/bin/env python3
"""Смета первого теста в AWS для трека ГТУ: «инстанс-часы × тариф» по строкам + диск, IPv4, S3, передача, логи, ECR.

    python3 cost_estimate.py [--config cost_estimate.toml] [--json cost_estimate_results.json] [--csv cost_estimate_lines.csv]

Три уровня по каждому сценарию:
  spot_estimate — ожидаемые часы × цена спота по Spot Instance Advisor (средняя экономия, не цена на дату);
  on_demand     — ожидаемые часы × цена по требованию (AWS Price List);
  upper_bound   — каждая попытка до своего срока (тайм-аут Batch или max_min раннера), все попытки, цена по
                  требованию. Спотовая машина дороже цены по требованию не обходится: у Batch bidPercentage
                  пустой = 100 % цены по требованию (модель API batch, botocore 1.44.3), у раннера MaxPrice —
                  1,2 × максимум спота за 6 ч (README раннера) — верхняя граница здесь берётся по требованию,
                  для раннера это допущение (спот выше цены по требованию бывает редко; не проверено).
Проверка: итог плана (сценарии in_plan) ≤ ceiling_usd на каждом уровне.
Ничего в AWS не вызывается. Результат — JSON (method, software, date, inputs у итога и у каждой строки) и CSV.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import pathlib
import sys
import tomllib

HERE = pathlib.Path(__file__).resolve().parent
LEVELS = ("spot_estimate", "on_demand", "upper_bound")
GB = 1024 ** 3


def load_toml(path: pathlib.Path) -> dict:
    with open(path, "rb") as f:
        return tomllib.load(f)


def r6(x: float) -> float:
    return round(x, 6)


# ---------------------------------------------------------------- время
def s_per_iter(cfg: dict, itype: str, procs: int) -> float:
    t = cfg["su2_time"]
    return (t["ref_s_per_iter"] * (t["mesh_nodes"] / t["ref_nodes"]) * (t["ref_procs"] / procs)
            / t["parallel_eff"] * t["cpu_factor"].get(itype, 1.0))


def scenario_shape(sc: dict, cfg: dict, batch: dict, rates: dict) -> dict:
    """Часы, типы, объёмы одного сценария (без денег)."""
    region = sc.get("region", cfg["meta"]["region"])
    ec2 = rates["regions"][region]["ec2"]
    if sc["path"] == "runner":
        itype = sc["instance_type"]
        s = s_per_iter(cfg, itype, sc["procs"])
        run_h = sc["iterations"] * s / 3600
        container_h = run_h + (sc["setup_min"] + sc["upload_min"]) / 60      # всё внутри max_min раннера
        exp_h = container_h                                                   # раннер гасит машину сразу
        bound_attempt_h = sc["cap_h"] + sc["grace_min"] / 60                  # формула бюджета runner.py
        types = [itype]
        jobs, points, iters_job = sc["jobs"], sc["jobs"], sc["iterations"]
        root_gb = sc["root_gb"]
        out_gb = sc["outputs_gb_per_job"] * jobs
        attempts = sc["attempts_max"]
        cap_h = sc["cap_h"]
    else:
        p, j, c = batch["plan"], batch["job"], batch["compute"]
        mode = sc.get("mode", batch["run"]["mode"])
        types = list(c["instance_types"])
        procs = int(j["vcpu"])
        s = max(s_per_iter(cfg, t, procs) for t in types)          # медленнейший из типов среды
        n_lines = len(p["stages"]) * len(p["speed_frac"])
        n_pts = int(p["points_per_line"])
        points = n_lines * n_pts
        if mode == "line":
            jobs, pts_job = n_lines, n_pts
            iters_job = p["iters_first"] + (n_pts - 1) * p["iters_next"]
        else:
            jobs, pts_job = points, 1
            iters_job = p["iters_first"]
        run_h = iters_job * s / 3600
        container_h = run_h + pts_job * sc["upload_min_per_point"] / 60
        exp_h = sc["setup_min"] / 60 + container_h + sc["tail_min"] / 60
        cap_h = (j["attempt_duration_s"]) / 3600
        bound_attempt_h = sc["setup_min"] / 60 + cap_h + sc["tail_min"] / 60
        root_gb = c["root_volume_gb"]
        out_gb = sc["outputs_gb_per_point"] * points
        attempts = int(j["retry_attempts"])
    st = cfg["storage"]
    od = max(ec2[t]["on_demand_usd_h"] for t in types)
    spots = [ec2[t]["spot_estimate_usd_h"] for t in types if ec2[t]["spot_estimate_usd_h"] > 0]
    return {
        "region": region, "path": sc["path"], "instance_types": types, "jobs": jobs, "points": points,
        "iterations_per_job": iters_job, "s_per_iter": r6(s), "run_h_per_job": r6(run_h),
        "container_h_per_job": r6(container_h), "expected_h_per_job": r6(exp_h if sc["path"] == "batch" else container_h),
        "cap_h": r6(cap_h), "cap_margin": r6(cap_h / container_h), "attempts_max": attempts,
        "bound_h_per_attempt": r6(bound_attempt_h), "root_gb": root_gb,
        "rate_on_demand_usd_h": od, "rate_on_demand_type": max(types, key=lambda t: ec2[t]["on_demand_usd_h"]),
        "rate_spot_estimate_usd_h": max(spots) if spots else od,
        "spot_types_without_data": [t for t in types if ec2[t]["spot_estimate_usd_h"] <= 0],
        "outputs_gb": out_gb,
        "log_gb": (iters_job * st["log_bytes_per_iter"] + st["log_bytes_per_job"]) * jobs / GB,
        "puts": st["put_per_point"] * points, "gets": st["get_per_point"] * points + st["get_inputs_per_job"] * jobs,
    }


# ---------------------------------------------------------------- деньги
def rate(rates: dict, region: str, key: str) -> float:
    return rates["regions"][region][key]["usd"]


def src(rates: dict, region: str, service: str) -> str:
    s = rates["sources"][f"{service}_{region}"]
    return f"{s['url']} (публикация {s['publication']})"


def line(scenario, level, item, qty, unit, rate_usd, rate_unit, formula, source) -> dict:
    return {"scenario": scenario, "level": level, "item": item, "quantity": r6(qty), "unit": unit,
            "rate_usd": rate_usd, "rate_unit": rate_unit, "cost_usd": r6(qty * rate_usd),
            "formula": formula, "rate_source": source}


def compute_lines(name: str, sh: dict, rates: dict, hpm: float) -> list[dict]:
    reg = sh["region"]
    gp3 = rate(rates, reg, "ebs_gp3_usd_gb_month")
    ipv4 = rate(rates, reg, "public_ipv4_usd_h")
    out = []
    for level in LEVELS:
        if level == "upper_bound":
            hours = sh["jobs"] * sh["attempts_max"] * sh["bound_h_per_attempt"]
            f_h = f"{sh['jobs']} заданий × {sh['attempts_max']} попытки × {sh['bound_h_per_attempt']:.4f} ч"
        else:
            hours = sh["jobs"] * sh["expected_h_per_job"]
            f_h = f"{sh['jobs']} заданий × {sh['expected_h_per_job']:.4f} ч"
        if level == "spot_estimate":
            ri, ru = sh["rate_spot_estimate_usd_h"], "USD/ч (спот, оценка по Spot Instance Advisor)"
            rs = rates["sources"]["spot_advisor"]["url"] + " (Last-Modified " + rates["sources"]["spot_advisor"]["publication"] + ") × " + src(rates, reg, "AmazonEC2")
        else:
            ri, ru = sh["rate_on_demand_usd_h"], f"USD/ч (по требованию, {sh['rate_on_demand_type']})"
            rs = src(rates, reg, "AmazonEC2")
        out.append(line(name, level, "compute", hours, "инстанс-ч", ri, ru, f_h + " × тариф", rs))
        out.append(line(name, level, "ebs_root", hours * sh["root_gb"], "ГБ·ч", r6(gp3 / hpm), "USD/(ГБ·ч)",
                        f"{hours:.4f} ч × {sh['root_gb']} ГБ × {gp3}/{hpm:.0f}", src(rates, reg, "AmazonEC2")))
        out.append(line(name, level, "public_ipv4", hours, "адрес·ч", ipv4, "USD/ч",
                        f"{hours:.4f} ч × 1 адрес", src(rates, reg, "AmazonVPC")))
    return out


def shared_lines(group: str, shapes: list[dict], cfg: dict, rates: dict, bucket_region: str) -> list[dict]:
    """S3, передача, логи, ECR — для группы сценариев; на верхней границе объёмы × число попыток."""
    st, hpm = cfg["storage"], cfg["meta"]["hours_per_month"]
    out = []
    for level in LEVELS:
        k = (lambda sh: sh["attempts_max"]) if level == "upper_bound" else (lambda sh: 1)
        gb = sum(sh["outputs_gb"] * k(sh) for sh in shapes)
        puts = sum(sh["puts"] * k(sh) for sh in shapes)
        gets = sum(sh["gets"] * k(sh) for sh in shapes)
        logs = {}
        for sh in shapes:
            logs[sh["region"]] = logs.get(sh["region"], 0.0) + sh["log_gb"] * k(sh)
        b = bucket_region
        out.append(line(group, level, "s3_storage", gb * st["s3_months"], "ГБ·мес",
                        rate(rates, b, "s3_standard_usd_gb_month"), "USD/(ГБ·мес)",
                        f"{gb:.3f} ГБ × {st['s3_months']} мес", src(rates, b, "AmazonS3")))
        out.append(line(group, level, "s3_put", puts, "запрос", rate(rates, b, "s3_put_usd_request"), "USD/запрос",
                        f"{puts} PUT", src(rates, b, "AmazonS3")))
        out.append(line(group, level, "s3_get", gets, "запрос", rate(rates, b, "s3_get_usd_request"), "USD/запрос",
                        f"{gets} GET", src(rates, b, "AmazonS3")))
        dl = gb * st["download_fraction"]
        out.append(line(group, level, "transfer_out_internet", dl, "ГБ",
                        rate(rates, b, "transfer_out_internet_usd_gb"), "USD/ГБ",
                        f"{dl:.3f} ГБ скачивания (100 ГБ/мес бесплатных не учтены)", src(rates, b, "AWSDataTransfer")))
        for reg, lg in sorted(logs.items()):
            out.append(line(group, level, f"cw_logs_ingest_{reg}", lg, "ГБ", rate(rates, reg, "cw_logs_ingest_usd_gb"),
                            "USD/ГБ", f"{lg:.5f} ГБ вывода", src(rates, reg, "AmazonCloudWatch")))
            out.append(line(group, level, f"cw_logs_storage_{reg}", lg * st["log_months"], "ГБ·мес",
                            rate(rates, reg, "cw_logs_storage_usd_gb_month"), "USD/(ГБ·мес)",
                            f"{lg:.5f} ГБ × {st['log_months']} мес", src(rates, reg, "AmazonCloudWatch")))
        regions = sorted({sh["region"] for sh in shapes})
        for reg in regions:
            out.append(line(group, level, f"ecr_storage_{reg}", st["image_gb"], "ГБ·мес",
                            rate(rates, reg, "ecr_storage_usd_gb_month"), "USD/(ГБ·мес)",
                            f"образ {st['image_gb']} ГБ × 1 мес", src(rates, reg, "AmazonECR")))
            if reg != b:   # машины в другом регионе: результаты едут в бакет, образ — в ECR этого региона
                xgb = sum(sh["outputs_gb"] * k(sh) for sh in shapes if sh["region"] == reg)
                out.append(line(group, level, f"transfer_{reg}_to_{b}", xgb, "ГБ",
                                rate(rates, reg, f"transfer_inter_region_to_{b}_usd_gb"), "USD/ГБ",
                                f"{xgb:.3f} ГБ результатов {reg} → {b}", src(rates, reg, "AWSDataTransfer")))
                out.append(line(group, level, f"transfer_{b}_to_{reg}_image", st["image_gb"], "ГБ",
                                rate(rates, b, f"transfer_inter_region_to_{reg}_usd_gb"), "USD/ГБ",
                                f"копия образа {st['image_gb']} ГБ {b} → {reg}", src(rates, b, "AWSDataTransfer")))
    return out


def totals(lines: list[dict], scen: str | None = None) -> dict:
    return {lv: r6(sum(x["cost_usd"] for x in lines if x["level"] == lv and (scen is None or x["scenario"] == scen)))
            for lv in LEVELS}


def estimate(cfg: dict, batch: dict, rates: dict) -> dict:
    hpm = cfg["meta"]["hours_per_month"]
    bucket_region = batch["run"]["region"]
    shapes = {sc["name"]: scenario_shape(sc, cfg, batch, rates) for sc in cfg["scenario"]}
    plan = [sc["name"] for sc in cfg["scenario"] if sc["in_plan"]]
    variants = [sc["name"] for sc in cfg["scenario"] if not sc["in_plan"]]
    lines = []
    for n in plan:
        lines += compute_lines(n, shapes[n], rates, hpm)
    lines += shared_lines("plan_shared", [shapes[n] for n in plan], cfg, rates, bucket_region)
    plan_total = totals(lines)
    var = {}
    for n in variants:
        vl = compute_lines(n, shapes[n], rates, hpm) + shared_lines(f"{n}_shared", [shapes[n]], cfg, rates,
                                                                      bucket_region)
        var[n] = {"total": totals(vl), "lines": vl}
    ceil, ceil2 = cfg["meta"]["ceiling_usd"], cfg["meta"]["ceiling_runner_usd"]
    warn = [f"{n}: срок {shapes[n]['cap_h']:.2f} ч всего в {shapes[n]['cap_margin']:.2f} раза больше ожидаемого "
            f"времени задания — после замера bench_* пересчитать" for n in shapes if shapes[n]["cap_margin"] < 1.25]
    ec2 = rates["regions"][cfg["meta"]["region"]]["ec2"]
    return {
        "method": ("Смета «инстанс-часы × тариф» + диск gp3, публичный IPv4, S3 (хранение, запросы), выход в интернет, "
                   "CloudWatch Logs (приём, хранение), ECR; часы — по модели времени итерации SU2 из локального замера "
                   "(su2_time) и по сроку попытки (верхняя граница); проверка итога плана против потолка."),
        "software": f"cost_estimate.py (Python {sys.version.split()[0]}, tomllib); тарифы — fetch_rates.py",
        "date": cfg["meta"]["date"],
        "inputs": {
            "config": "cost_estimate.toml", "batch_config": cfg["meta"]["batch_config"], "rates": cfg["meta"]["rates"],
            "rates_generated_utc": rates["meta"]["generated_utc"],
            "rates_sources": {k: v["publication"] for k, v in rates["sources"].items()},
            "su2_time": cfg["su2_time"], "storage": cfg["storage"], "hours_per_month": hpm,
            "series_from_batch_toml": {"attempt_duration_s": batch["job"]["attempt_duration_s"],
                                       "retry_attempts": batch["job"]["retry_attempts"],
                                       "instance_types": batch["compute"]["instance_types"],
                                       "stages": batch["plan"]["stages"], "speed_frac": batch["plan"]["speed_frac"],
                                       "points_per_line": batch["plan"]["points_per_line"],
                                       "iters_first": batch["plan"]["iters_first"],
                                       "iters_next": batch["plan"]["iters_next"]},
            "example_fields": "все поля с пометкой EXAMPLE в cost_estimate.toml и batch.toml (example: true)",
            "example": True,
        },
        "rates_used": {t: {"on_demand_usd_h": ec2[t]["on_demand_usd_h"], "spot_estimate_usd_h": ec2[t]["spot_estimate_usd_h"],
                           "spot_savings_pct": ec2[t]["spot_savings_pct"], "spot_interrupt_bucket": ec2[t]["spot_interrupt_bucket"],
                           "vcpu": ec2[t]["vcpu"], "memory_gib": ec2[t]["memory_gib"]}
                       for t in sorted({t for n in shapes for t in shapes[n]["instance_types"]} & set(ec2))},
        "scenarios": {n: {**shapes[n], "title": next(sc["title"] for sc in cfg["scenario"] if sc["name"] == n),
                          "in_plan": n in plan, "total": totals(lines, n) if n in plan else var[n]["total"]}
                      for n in shapes},
        "plan": plan,
        "plan_shared_total": totals(lines, "plan_shared"),
        "plan_total": plan_total,
        "ceiling_check": {"ceiling_usd": ceil, **{f"{lv}_ok": plan_total[lv] <= ceil for lv in LEVELS},
                          "headroom_upper_bound_usd": r6(ceil - plan_total["upper_bound"]),
                          "ceiling_runner_usd": ceil2, "upper_bound_ok_runner": plan_total["upper_bound"] <= ceil2},
        "warnings": warn,
        "lines": lines,
        "variants": var,
    }


def write_csv(res: dict, path: pathlib.Path) -> None:
    cols = ["scenario", "level", "item", "quantity", "unit", "rate_usd", "rate_unit", "cost_usd", "formula", "rate_source"]
    rows = list(res["lines"]) + [x for v in res["variants"].values() for x in v["lines"]]
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=cols, lineterminator="\n")
        w.writeheader()
        w.writerows(rows)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--config", default=str(HERE / "cost_estimate.toml"))
    ap.add_argument("--json", default=str(HERE / "cost_estimate_results.json"))
    ap.add_argument("--csv", default=str(HERE / "cost_estimate_lines.csv"))
    a = ap.parse_args(argv)
    cfg = load_toml(pathlib.Path(a.config))
    base = pathlib.Path(a.config).resolve().parent
    res = estimate(cfg, load_toml(base / cfg["meta"]["batch_config"]), load_toml(base / cfg["meta"]["rates"]))
    pathlib.Path(a.json).write_text(json.dumps(res, ensure_ascii=False, indent=1) + "\n", encoding="utf-8")
    write_csv(res, pathlib.Path(a.csv))
    pt, cc = res["plan_total"], res["ceiling_check"]
    for n, s in res["scenarios"].items():
        t = s["total"]
        print(f"{n:26s} {'план' if s['in_plan'] else 'вариант'}: спот ≈ {t['spot_estimate']:8.2f}  "
              f"по требованию {t['on_demand']:8.2f}  граница {t['upper_bound']:8.2f} USD")
    print(f"ИТОГ плана: спот ≈ {pt['spot_estimate']:.2f}, по требованию {pt['on_demand']:.2f}, "
          f"верхняя граница {pt['upper_bound']:.2f} USD; потолок {cc['ceiling_usd']:.0f}: "
          f"{'в пределах' if cc['upper_bound_ok'] else 'ПРЕВЫШЕН'} (запас {cc['headroom_upper_bound_usd']:.2f})")
    for w in res["warnings"]:
        print("ВНИМАНИЕ:", w)
    return 0 if all(cc[f"{lv}_ok"] for lv in LEVELS) else 1


if __name__ == "__main__":
    sys.exit(main())
