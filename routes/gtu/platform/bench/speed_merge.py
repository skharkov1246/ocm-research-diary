#!/usr/bin/env python3
"""Сводка замеров скорости (speed_su2.py, speed_ccx.py) с машин AWS → цена точки SU2 и расчёта ccx (узел bench).

Читает <in>/<задача>/speed_su2_*.json и speed_ccx_*.json (папка прогона раннера runs/<run_id>/), пишет одну
небольшую сводку для ветки. Замер — время итерации и время расчёта; цена — пересчёт по формулам ниже с явными
допущениями (каждое — в блоке assumptions итога):

  SU2, точка уровня L на машине типа T:
    лучшая конфигурация — наибольшее число итераций в час на машину (K·3600 / с_итер, MUSCL);
    машино-часов на точку = итераций_на_точку · с_итер / K / 3600;
    база линии (разгон, первый порядок, 64 процесса) = итераций_базы · с_итер(первый порядок, наибольшее NP) / 3600,
      делится на точек_в_линии; сетка концепции = время сетки (замер на машине, Gmsh);
    долл. = машино-часы · (цена спота + диск·0,08/730 + IPv4 0,005);
    другое число узлов N: с_итер ∝ N при той же конфигурации (C = с_итер·NP/N постоянно) — пересчёт, не замер.
  CalculiX, расчёт на сетке из N узлов C3D10:
    один расчёт на всю машину — время при наибольшем числе потоков; упаковка — машина / (расчётов в час);
    берётся дешёвый из двух; другое N — степенная подгонка времени по замеренным сеткам (вне диапазона — пометка).

  python3 speed_merge.py --in routes/gtu/platform/aws/runs/<run_id> --run-id <run_id> \\
          --out routes/gtu/platform/bench/speed_aws_results.json
Код 0 — записано; 2 — в папке нет ни одного замера.
"""
from __future__ import annotations

import argparse
import datetime
import glob
import hashlib
import json
import math
import os
import sys

DEFAULT_PRICES = {"c7a.16xlarge": 0.75, "c7i.16xlarge": 0.94248}
PRICE_BASIS = {
    "c7a.16xlarge": ("0,75 долл./ч — верх ожидаемого спота c7a.16xlarge us-east-1 (задача владельца 14.11.2026: 0,71–0,75; "
                     "ledger.jsonl smoke2: us-east-1c 0,7064)"),
    "c7i.16xlarge": ("0,94248 долл./ч — оценка спота по Spot Instance Advisor (routes/gtu/platform/batch/aws_rates.toml, "
                     "экономия 67 % к 2,856); живой цены c7i в журнале нет"),
}
UPPER_PRICE = 1.4206
UPPER_BASIS = "MaxPrice раннера в ledger.jsonl живого прогона smoke 14.11.2026 (максимум по зонам × 1,2)"
DISK_RATE, IPV4_RATE = 0.08 / 730, 0.005


def extras_h(disk_gb: float) -> float:
    return disk_gb * DISK_RATE + IPV4_RATE


def _sha(path):
    try:
        return hashlib.sha256(open(path, "rb").read()).hexdigest()
    except OSError:
        return None


def best_config(configs: list[dict], level: str, stage: str = "muscl"):
    c = [x for x in configs if x.get("level") == level and (x.get(stage) or {}).get("measured")]
    return max(c, key=lambda x: x[stage]["iters_per_hour_machine"]) if c else None


def su2_point_cost(doc: dict, price: float, disk_gb: float, point_iters: int, base_iters: int,
                   points_per_line: int, target_nodes: dict) -> dict:
    """Цена точки SU2 по уровням из одного speed_su2_*.json."""
    out = {}
    rate = price + extras_h(disk_gb)
    upper_rate = UPPER_PRICE + extras_h(disk_gb)
    for level in doc.get("inputs", {}).get("levels", []):
        b = best_config(doc.get("configs", []), level)
        if b is None:
            out[level] = dict(measured=False)
            continue
        s = b["muscl"]["s_iter"]
        mh_point = point_iters * s / b["k"] / 3600
        # база: первый порядок, наибольшее NP на этом уровне
        fo = [x for x in doc["configs"] if x["level"] == level and x.get("first_order", {}).get("measured")]
        fo_big = max(fo, key=lambda x: x["np"]) if fo else None
        mh_base = base_iters * fo_big["first_order"]["s_iter"] / 3600 if fo_big else None
        mesh = doc.get("meshes", {}).get(level, {})
        mh_mesh = (mesh.get("mesh_s") or 0) / 3600 if not mesh.get("reused") else None
        nodes = b.get("nodes") or mesh.get("nodes")
        row = dict(measured=True, nodes=nodes, best=f"{b['np']}x{b['k']}", s_iter=s,
                   C_us_core_node_iter=b["muscl"].get("C_us_core_node_iter"),
                   iters_per_hour_machine=b["muscl"]["iters_per_hour_machine"],
                   machine_h_per_point=mh_point, usd_per_point=mh_point * rate, usd_per_point_upper=mh_point * upper_rate,
                   base_machine_h=mh_base, base_usd_per_point=(mh_base * rate / points_per_line) if mh_base else None,
                   mesh_machine_h=mh_mesh, mesh_usd=(mh_mesh * rate) if mh_mesh is not None else None,
                   all_configs=[dict(cfg=f"{x['np']}x{x['k']}", s_iter=x["muscl"].get("s_iter"),
                                     C_us=x["muscl"].get("C_us_core_node_iter"),
                                     iters_per_hour_machine=x["muscl"].get("iters_per_hour_machine"),
                                     rss_peak_gb=x.get("rss_peak_gb"))
                                for x in doc["configs"] if x["level"] == level])
        if nodes:
            row["usd_per_point_per_Mnode"] = row["usd_per_point"] / (nodes / 1e6)
            row["scaled"] = {name: dict(nodes=n, usd_per_point=row["usd_per_point"] * n / nodes,
                                        machine_h_per_point=mh_point * n / nodes,
                                        note="пересчёт ∝ узлам при той же конфигурации, не замер")
                             for name, n in target_nodes.get(level, {}).items()}
        out[level] = row
    return out


def _spread(vals):
    """(макс − мин) / макс по max|U| расчётов при разных потоках: сверка, что расчёт один и тот же."""
    v = [x for x in vals if isinstance(x, (int, float)) and x > 0]
    return (max(v) - min(v)) / max(v) if v else None


def powerfit(xs, ys):
    """y = a·x^b по МНК в логарифмах; None при < 2 точках."""
    pts = [(math.log(x), math.log(y)) for x, y in zip(xs, ys) if x and y and x > 0 and y > 0]
    if len(pts) < 2:
        return None
    n = len(pts)
    mx = sum(p[0] for p in pts) / n
    my = sum(p[1] for p in pts) / n
    sxx = sum((p[0] - mx) ** 2 for p in pts)
    if sxx == 0:
        return None
    b = sum((p[0] - mx) * (p[1] - my) for p in pts) / sxx
    return dict(a=math.exp(my - b * mx), b=b, n=n, x_min=math.exp(min(p[0] for p in pts)),
                x_max=math.exp(max(p[0] for p in pts)))


def ccx_case_cost(doc: dict, price: float, disk_gb: float, target_nodes: list[int]) -> dict:
    rate = price + extras_h(disk_gb)
    rows, xs_single, ys_single, xs_pack, ys_pack = [], [], [], [], []
    for m in doc.get("meshes", []):
        fin = [s for s in m.get("static", []) if s.get("status") == "finished"]
        if not fin:
            rows.append(dict(nodes=m.get("nodes"), measured=False))
            continue
        tmax = max(fin, key=lambda s: s["threads"])
        single_h = tmax["wall_s"] / 3600
        row = dict(nodes=m["nodes"], elements=m.get("elements"), mesh_s=m.get("mesh_s"),
                   single=dict(threads=tmax["threads"], wall_s=tmax["wall_s"], maxrss_gb=tmax["maxrss_gb"],
                               usd=single_h * rate),
                   static_by_threads={str(s["threads"]): dict(wall_s=s["wall_s"], maxrss_gb=s["maxrss_gb"])
                                      for s in fin},
                   max_disp_rel_spread=_spread([s.get("max_disp_m") for s in fin]))
        xs_single.append(m["nodes"])
        ys_single.append(tmax["wall_s"])
        p = m.get("pack") or {}
        if p.get("status") == "finished" and p.get("cases_per_hour_machine"):
            row["pack"] = dict(threads=p["threads"], copies=p["copies"], cases_per_hour_machine=p["cases_per_hour_machine"],
                               usd=rate / p["cases_per_hour_machine"])
            xs_pack.append(m["nodes"])
            ys_pack.append(3600 / p["cases_per_hour_machine"])  # машино-секунд на расчёт
        f = m.get("freq") or {}
        if f.get("status") == "finished":
            row["freq"] = dict(threads=f["threads"], modes=f.get("modes"), wall_s=f["wall_s"], maxrss_gb=f["maxrss_gb"],
                               usd=f["wall_s"] / 3600 * rate, freq_hz_first=(f.get("freq_hz") or [])[:3])
        cands = [row["single"]["usd"]] + ([row["pack"]["usd"]] if "pack" in row else [])
        row["usd_best"] = min(cands)
        rows.append(row)
    fit_s, fit_p = powerfit(xs_single, ys_single), powerfit(xs_pack, ys_pack)
    proj = []
    for n in target_nodes:
        d = dict(nodes=n)
        for key, fit in (("single", fit_s), ("pack", fit_p)):
            if fit:
                sec = fit["a"] * n ** fit["b"]
                d[key] = dict(machine_s=sec, usd=sec / 3600 * rate,
                              extrapolated=not (fit["x_min"] <= n <= fit["x_max"]))
        proj.append(d)
    return dict(meshes=rows, fit_single=fit_s, fit_pack=fit_p, projected=proj)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--in", dest="inp", required=True, help="папка прогона раннера runs/<run_id>")
    ap.add_argument("--run-id", default="")
    ap.add_argument("--out", required=True)
    ap.add_argument("--price", action="append", default=[], help="тип=долл./ч (по умолчанию — DEFAULT_PRICES)")
    ap.add_argument("--disk-gb", type=float, default=60)
    ap.add_argument("--point-iters", type=int, default=6000, help="итераций MUSCL на точку")
    ap.add_argument("--base-iters", type=int, default=520, help="итераций разгона линии (план smoke075)")
    ap.add_argument("--points-per-line", type=int, default=8)
    ap.add_argument("--targets-su2", default="medium:kv3_stage1=1137992,medium:core130=584357,"
                                             "fine:kv3_stage1=5009000,fine:core130=2455000",
                    help="уровень:имя=узлов,... для пересчёта цены точки")
    ap.add_argument("--targets-ccx", default="100000,300000,1000000")
    a = ap.parse_args(argv)
    prices = dict(DEFAULT_PRICES)
    for p in a.price:
        k, v = p.split("=")
        prices[k] = float(v)
    tsu2 = {}
    for item in a.targets_su2.split(","):
        if item.strip():
            lv, rest = item.split(":")
            name, n = rest.split("=")
            tsu2.setdefault(lv, {})[name] = int(n)
    tccx = [int(x) for x in a.targets_ccx.split(",") if x.strip()]
    files = sorted(glob.glob(os.path.join(a.inp, "*", "speed_su2_*.json")) +
                   glob.glob(os.path.join(a.inp, "*", "speed_ccx_*.json")))
    if not files:
        print(f"[speed_merge] в {a.inp} нет speed_su2_*.json / speed_ccx_*.json", file=sys.stderr)
        return 2
    su2, ccx = [], []
    for f in files:
        d = json.load(open(f, encoding="utf-8"))
        job = os.path.basename(os.path.dirname(f))
        itype = (d.get("machine") or {}).get("itype") or "unknown"
        price = prices.get(itype, prices["c7a.16xlarge"])
        base = dict(job=job, file=os.path.relpath(f, a.inp), sha256=_sha(f), itype=itype,
                    cpu_model=(d.get("machine") or {}).get("cpu_model"), nproc=(d.get("machine") or {}).get("nproc"),
                    threads_per_core=(d.get("machine") or {}).get("threads_per_core"),
                    commit=d.get("commit"), date=d.get("date"), price_usd_h=price, ok=d.get("ok"))
        if os.path.basename(f).startswith("speed_su2_"):
            base.update(case=d["inputs"]["case"], meshes=d.get("meshes"), skipped=d.get("skipped"),
                        levels=su2_point_cost(d, price, a.disk_gb, a.point_iters, a.base_iters, a.points_per_line, tsu2))
            su2.append(base)
        else:
            base.update(case=d["inputs"]["case"], ccx_version=d.get("ccx_version"), skipped=d.get("skipped"),
                        **ccx_case_cost(d, price, a.disk_gb, tccx))
            ccx.append(base)
    doc = dict(
        title="Цена точки SU2 и расчёта CalculiX на машинах раннера (замер скорости)",
        method=__doc__.split("\n\n")[1].strip(),
        run_id=a.run_id, created_utc=datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        prices={k: dict(usd_h=v, basis=PRICE_BASIS.get(k, "--price")) for k, v in prices.items()},
        upper=dict(usd_h=UPPER_PRICE, basis=UPPER_BASIS),
        assumptions=dict(point_iters=a.point_iters, base_iters=a.base_iters, points_per_line=a.points_per_line,
                         disk_gb=a.disk_gb, disk_ipv4="0,08 долл./ГБ·мес / 730 и 0,005 долл./ч — ставки раннера (допущения)",
                         point_iters_source=("routes/gtu/platform/su2/jobs_aws_su2_kv3.json → setup.base: рампа давления "
                                             "300 и всего 6000 итераций второго порядка на точку; сходимость за это число "
                                             "итераций не проверена"),
                         targets_su2=tsu2, targets_ccx=tccx),
        su2=su2, ccx=ccx)
    txt = json.dumps(doc, ensure_ascii=False, indent=1)
    if len(txt.encode()) > 5 * 1024 * 1024:
        print("[speed_merge] сводка больше 5 МБ — не записана", file=sys.stderr)
        return 2
    with open(a.out, "w", encoding="utf-8") as f:
        f.write(txt + "\n")
    for s in su2:
        for lv, r in s["levels"].items():
            if r.get("measured"):
                print(f"[speed_merge] SU2 {s['itype']} {s['case']}/{lv}: {r['best']} {r['s_iter']:.3f} с/итер, "
                      f"точка {r['machine_h_per_point']:.2f} маш.-ч = {r['usd_per_point']:.2f} долл. "
                      f"(верх {r['usd_per_point_upper']:.2f})")
    for c in ccx:
        for m in c["meshes"]:
            if m.get("usd_best") is not None:
                print(f"[speed_merge] ccx {c['itype']} {m['nodes']} узлов: {m['usd_best']:.4f} долл. за статику")
    print(f"[speed_merge] записано: {a.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
