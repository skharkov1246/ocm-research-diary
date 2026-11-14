#!/usr/bin/env python3
"""Порог окупаемости переноса ядра SU2 на GPU (протокол PROTOCOL.md).

Три шага, сеть и AWS не вызываются — только локальные файлы:

  tariffs  Price List AWS (CSV) и данные Spot Instance Advisor → tariffs_us-east-1.json
           (цена по требованию — из прайса; цена спота — оценка: OD × (1 − экономия/100));
  mesh     сетки .su2 → mesh_stats.json (узлы, ячейки, уникальные рёбра — шаблон матрицы SU2);
  compute  bench.toml + тарифы + сетки [+ сводка bench_parse.py] → break_even.json:
           порог ускорения S* = тариф GPU / тариф CPU, доля ядра по Амдалу, оценка видеопамяти,
           упор в полосу памяти, при наличии измерений — вердикт по каждой GPU-конфигурации.

Примеры:
  python3 break_even.py tariffs --price-list us-east-1=index_us.csv --price-list eu-north-1=index_eu.csv \
          --spot-advisor spot-advisor-data.json --out tariffs_us-east-1.json
  python3 break_even.py mesh ../geom/out/core130_coarse/passage_core130_coarse.su2 --out mesh_stats.json
  python3 break_even.py compute --bench bench.toml --out break_even.json [--summary bench_summary.csv]
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
import platform
import sys
import tomllib
from pathlib import Path

DATE = "2026-11-14"
HERE = Path(__file__).resolve().parent
PRICE_URL = "https://pricing.us-east-1.amazonaws.com/offers/v1.0/aws/AmazonEC2/{version}/{region}/index.csv"
SPOT_URL = "https://spot-bid-advisor.s3.amazonaws.com/spot-advisor-data.json"
VPC_URL = "https://pricing.us-east-1.amazonaws.com/offers/v1.0/aws/AmazonVPC/{version}/{region}/index.json"
DEFAULT_TYPES = ["c8i.2xlarge", "c8i.8xlarge", "c7a.16xlarge", "c7i.16xlarge", "g6e.xlarge",
                 "g6e.2xlarge", "g6e.4xlarge", "p5.4xlarge", "p5.48xlarge", "hpc7a.96xlarge",
                 "hpc8a.96xlarge"]

# Рёбра объёмных элементов SU2 (тип VTK → пары локальных узлов): шаблон матрицы SU2 строится
# по соседям через рёбра элементов (CPhysicalGeometry::SetPoint_Connectivity).
ELEM_EDGES = {
    10: [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)],                    # тетраэдр
    12: [(0, 1), (1, 2), (2, 3), (3, 0), (4, 5), (5, 6), (6, 7), (7, 4),
         (0, 4), (1, 5), (2, 6), (3, 7)],                                     # гексаэдр
    13: [(0, 1), (1, 2), (2, 0), (3, 4), (4, 5), (5, 3), (0, 3), (1, 4), (2, 5)],  # призма
    14: [(0, 1), (1, 2), (2, 3), (3, 0), (0, 4), (1, 4), (2, 4), (3, 4)],     # пирамида
}
ELEM_NODES = {10: 4, 12: 8, 13: 6, 14: 5}


class BenchError(RuntimeError):
    pass


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def rel(path: Path) -> str:
    """Путь относительно папки протокола (как его видит репозиторий)."""
    return os.path.relpath(Path(path).resolve(), HERE)


# ------------------------------------------------------------------------------------ тарифы
def read_price_list(path: Path, types: set[str]) -> tuple[dict, dict]:
    """Цена Linux по требованию, общая аренда, без ПО, из CSV Price List AWS (offer AmazonEC2)."""
    meta: dict[str, str] = {}
    found: dict[str, dict] = {}
    with open(path, newline="", encoding="utf-8") as f:
        for _ in range(5):
            row = next(csv.reader([f.readline()]))
            if len(row) >= 2:
                meta[row[0]] = row[1]
        for r in csv.DictReader(f):
            if (r.get("Product Family") == "Storage" and r.get("Volume API Name") == "gp3"
                    and r.get("TermType") == "OnDemand" and r.get("Unit") == "GB-Mo"):
                found["_ebs_gp3"] = {"usd_gb_month": float(r["PricePerUnit"]), "sku": r["SKU"],
                                     "price_description": r["PriceDescription"]}
                continue
            it = r.get("Instance Type", "")
            if it not in types:
                continue
            if not (r.get("TermType") == "OnDemand" and r.get("Operating System") == "Linux"
                    and r.get("Tenancy") == "Shared" and r.get("Pre Installed S/W") == "NA"
                    and r.get("CapacityStatus") == "Used" and r.get("License Model") == "No License required"
                    and r.get("Unit") == "Hrs" and r.get("MarketOption", "OnDemand") in ("OnDemand", "")):
                continue
            if it in found:
                raise BenchError(f"{path}: две строки цены для {it} — фильтр неоднозначен")
            found[it] = {
                "od_usd_h": float(r["PricePerUnit"]),
                "price_description": r["PriceDescription"],
                "sku": r["SKU"], "rate_code": r["RateCode"], "effective_date": r["EffectiveDate"],
                "vcpu": int(r["vCPU"]) if r.get("vCPU", "").isdigit() else r.get("vCPU"),
                "memory": r.get("Memory"), "gpu": r.get("GPU") or None, "gpu_memory": r.get("GPU Memory") or None,
                "network": r.get("Network Performance"), "processor": r.get("Physical Processor"),
            }
    return meta, found


def tariffs(args) -> dict:
    types = set(args.types.split(",")) if args.types else set(DEFAULT_TYPES)
    spot = json.loads(Path(args.spot_advisor).read_text(encoding="utf-8"))
    ranges = {r["index"]: r["label"] for r in spot.get("ranges", [])}
    out = {
        "method": ("Цена Linux по требованию (общая аренда, без лицензий ПО, CapacityStatus=Used) — строка "
                   "прайса Price List AWS (offer AmazonEC2, CSV региона). Цена спота — оценка: OD × (1 − s/100), "
                   "где s — экономия к цене по требованию из данных Spot Instance Advisor (тот же файл питает "
                   "страницу советника); r — категория частоты прерываний. Фактическая цена спота меняется по "
                   "зонам и часам: перед запуском её читает runner.py (describe_spot_price_history). "
                   "Диск gp3 (долл./ГБ·мес) — строка Storage того же прайса; публичный IPv4 (долл./ч) — offer "
                   "AmazonVPC (JSON региона)."),
        "software": f"python3 {platform.python_version()} (csv, json)",
        "date": DATE,
        "inputs": [],
        "regions": {},
        "extras": {},
        "missing": {},
    }
    for spec in args.price_list:
        region, _, p = spec.partition("=")
        path = Path(p)
        meta, found = read_price_list(path, types)
        version = meta.get("Version", "")
        out["inputs"].append({"file": path.name, "region": region, "url": PRICE_URL.format(version=version, region=region),
                              "publication_date": meta.get("Publication Date"), "version": version,
                              "sha256": sha256(path)})
        adv = spot.get("spot_advisor", {}).get(region, {}).get("Linux", {})
        reg = {}
        extras = out["extras"].setdefault(region, {})
        ebs = found.pop("_ebs_gp3", None)
        if ebs:
            extras["ebs_gp3_usd_gb_month"] = ebs["usd_gb_month"]
            extras["ebs_gp3_sku"] = ebs["sku"]
        for it in sorted(found):
            rec = dict(found[it])
            a = adv.get(it)
            if a:
                rec["spot_savings_pct"] = a["s"]
                rec["spot_interruption"] = ranges.get(a["r"], str(a["r"]))
                rec["spot_est_usd_h"] = round(rec["od_usd_h"] * (1 - a["s"] / 100.0), 6)
            else:
                rec["spot_savings_pct"] = None
                rec["spot_interruption"] = None
                rec["spot_est_usd_h"] = None
            reg[it] = rec
        out["regions"][region] = reg
        out["missing"][region] = sorted(types - set(found))
    out["inputs"].append({"file": Path(args.spot_advisor).name, "url": SPOT_URL, "sha256": sha256(Path(args.spot_advisor))})
    for spec in getattr(args, "vpc_price", None) or []:
        region, _, p = spec.partition("=")
        d = json.loads(Path(p).read_text(encoding="utf-8"))
        for sku, prod in d.get("products", {}).items():
            ut = prod.get("attributes", {}).get("usagetype", "")
            if "PublicIPv4:InUseAddress" not in ut:
                continue
            for term in d.get("terms", {}).get("OnDemand", {}).get(sku, {}).values():
                for dim in term.get("priceDimensions", {}).values():
                    if dim.get("unit") == "Hrs":
                        out["extras"].setdefault(region, {})["public_ipv4_usd_h"] = float(dim["pricePerUnit"]["USD"])
                        out["extras"][region]["public_ipv4_sku"] = sku
        out["inputs"].append({"file": Path(p).name, "region": region,
                              "url": VPC_URL.format(version=d.get("version", ""), region=region),
                              "publication_date": d.get("publicationDate"), "version": d.get("version"),
                              "sha256": sha256(Path(p))})
    return out


# ------------------------------------------------------------------------------------ сетки
def mesh_stats(path: Path) -> dict:
    """Узлы, объёмные элементы и уникальные рёбра сетки SU2 (ASCII, одна зона)."""
    import numpy as np

    with open(path, encoding="utf-8") as f:
        ndime = nelem = npoin = None
        for line in f:
            if line.startswith("NDIME="):
                ndime = int(line.split("=")[1])
            elif line.startswith("NELEM="):
                nelem = int(line.split("=")[1])
                break
        if ndime != 3 or nelem is None:
            raise BenchError(f"{path}: нужна трёхмерная сетка SU2 с NELEM")
        by_type: dict[int, list] = {}
        for _ in range(nelem):
            parts = f.readline().split()
            t = int(parts[0])
            if t not in ELEM_NODES:
                raise BenchError(f"{path}: объёмный элемент типа {t} не поддержан")
            by_type.setdefault(t, []).append([int(v) for v in parts[1:1 + ELEM_NODES[t]]])
        for line in f:
            if line.startswith("NPOIN="):
                npoin = int(line.split("=")[1].split()[0])
                break
    if npoin is None:
        raise BenchError(f"{path}: нет NPOIN")
    keys = []
    for t, conn in by_type.items():
        c = np.asarray(conn, dtype=np.int64)
        for a, b in ELEM_EDGES[t]:
            lo = np.minimum(c[:, a], c[:, b])
            hi = np.maximum(c[:, a], c[:, b])
            keys.append(lo * npoin + hi)
    n_edges = int(np.unique(np.concatenate(keys)).size)
    return {
        "file": rel(path), "sha256": sha256(path), "n_points": npoin, "n_elements": nelem,
        "element_types": {str(t): len(c) for t, c in sorted(by_type.items())},
        "n_edges": n_edges, "edges_per_point": round(n_edges / npoin, 6),
    }


# ------------------------------------------------------------------------------------ формулы
def break_even_speedup(p_gpu: float, p_cpu: float) -> float:
    """S* — ускорение (по времени до сходимости), при котором прогон на GPU стоит столько же."""
    return p_gpu / p_cpu


def amdahl_speedup(f: float, s_kernel: float, r_host: float = 1.0, x: float = 0.0) -> float:
    """Ускорение всего прогона: доля f переносится и ускоряется в s_kernel раз (0 — бесконечно),
    остальное идёт на CPU GPU-машины в r_host раз медленнее, x — обмен CPU–GPU в долях T_cpu."""
    k = 0.0 if s_kernel == 0 else f / s_kernel
    return 1.0 / ((1.0 - f) * r_host + k + x)


def amdahl_required_fraction(s_target: float, s_kernel: float, r_host: float = 1.0, x: float = 0.0):
    """Минимальная доля ядра f, при которой amdahl_speedup ≥ s_target; None — недостижимо при f ≤ 1."""
    inv_k = 0.0 if s_kernel == 0 else 1.0 / s_kernel
    den = r_host - inv_k
    if den <= 0:
        return 0.0 if 1.0 / s_target >= inv_k + x else None
    f = (r_host + x - 1.0 / s_target) / den
    if f > 1.0 + 1e-12:
        return None
    return max(0.0, min(1.0, f))


def memory_per_point(edges_per_point: float, nvar: int, scalar_bytes: int, krylov_restart: int,
                     extra_vectors: int, index_bytes: int) -> dict:
    """Нижняя граница памяти одной линейной системы SU2 v8.6.0 на GPU, байт на узел сетки.

    Матрица LDU (CSysMatrix: gpu.d, gpu.l, gpu.u): (1 + 2E/N)·nvar² скаляров; индексы row_ptr_l/u и
    col_ind_l/u: (2 + 2E/N)·index_bytes; копия ILU(0) на устройстве (gpu_ilu) — тот же объём;
    векторы FGMRES: V и Z по (m + 1) + extra_vectors, по nvar скаляров (CSysSolve.cpp)."""
    blocks = 1.0 + 2.0 * edges_per_point
    matrix = blocks * nvar * nvar * scalar_bytes
    index = (2.0 + 2.0 * edges_per_point) * index_bytes
    ilu = matrix + index
    vectors = (2 * (krylov_restart + 1) + extra_vectors) * nvar * scalar_bytes
    return {"matrix": matrix, "index": index, "ilu0_copy": ilu, "vectors": vectors,
            "total": matrix + index + ilu + vectors}


def spmv_intensity(edges_per_point: float, nvar: int, scalar_bytes: int, index_bytes: int) -> tuple[float, float]:
    """Арифметическая интенсивность блочного SpMV (флоп/байт): (без повторного чтения x, с чтением x
    для каждого блока строки). Флопы строки: 2·nvar²·блоков; байты: блоки + индексы + x + y."""
    bpr = 1.0 + 2.0 * edges_per_point
    flops = 2.0 * nvar * nvar * bpr
    base = nvar * nvar * bpr * scalar_bytes + (bpr - 1.0) * index_bytes + 2 * index_bytes + nvar * scalar_bytes
    best = flops / (base + nvar * scalar_bytes)
    worst = flops / (base + bpr * nvar * scalar_bytes)
    return best, worst


# ------------------------------------------------------------------------------------ расчёт
def extras_usd_h(tar: dict, region: str, disk_gb):
    """Диск gp3 и публичный IPv4 машины, долл./ч (месяц — 730 ч, как в runner.py); None — нет цен."""
    ex = tar.get("extras", {}).get(region, {})
    if disk_gb is None or "ebs_gp3_usd_gb_month" not in ex or "public_ipv4_usd_h" not in ex:
        return None
    return disk_gb * ex["ebs_gp3_usd_gb_month"] / 730.0 + ex["public_ipv4_usd_h"]


def price(tar: dict, region: str, itype: str, market: str):
    rec = tar["regions"].get(region, {}).get(itype)
    if rec is None:
        return None
    return rec["od_usd_h"] if market == "od" else rec.get("spot_est_usd_h")


def load_summary(path: Path) -> list[dict]:
    with open(path, newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f))


def fnum(v):
    try:
        x = float(v)
    except (TypeError, ValueError):
        return None
    return None if math.isnan(x) else x


def measured_verdicts(bench: dict, summary: list[dict], compare: list[dict] | None) -> dict:
    dec = bench["decision"]
    margin = dec["speedup_margin"]
    rows = {r["config_id"]: r for r in summary}
    cmp = {r["config_id"]: r for r in (compare or [])}
    cpu = [r for r in summary if r.get("role") == "cpu" and r.get("converged") == "True"
           and fnum(r.get("cost_solver_usd")) is not None]
    if not cpu:
        return {"status": "нет сошедшейся CPU-конфигурации со стоимостью — решение не принимается"}
    best = min(cpu, key=lambda r: fnum(r["cost_solver_usd"]))
    t_cpu, c_cpu = fnum(best["time_to_conv_s"]), fnum(best["cost_solver_usd"])
    p_cpu = fnum(best["tariff_usd_h"])
    out = {"best_cpu": {"config_id": best["config_id"], "time_to_conv_s": t_cpu, "cost_solver_usd": c_cpu},
           "gpu": []}
    for cid, r in rows.items():
        if r.get("role") != "gpu":
            continue
        t, c, p = fnum(r.get("time_to_conv_s")), fnum(r.get("cost_solver_usd")), fnum(r.get("tariff_usd_h"))
        item = {"config_id": cid, "converged": r.get("converged") == "True"}
        if t and t_cpu and p and p_cpu:
            s = t_cpu / t
            s_star = p / p_cpu
            item.update(speedup=round(s, 4), s_star=round(s_star, 4), s_required=round(margin * s_star, 4),
                        cost_ratio=round(c / c_cpu, 4) if c and c_cpu else None, pass_cost=s >= margin * s_star)
        vram = (dec.get("gpu_vram_bytes") or {}).get(r.get("instance"))
        peak = fnum(r.get("gpu_peak_mem_GiB"))
        if vram and peak is not None:
            item["vram_fill"] = round(peak * 2**30 / vram, 4)
            item["pass_vram"] = item["vram_fill"] <= dec["vram_fill_max"]
        cr = cmp.get(cid)
        if cr is not None:
            item["pass_accuracy"] = cr.get("all_pass") == "True"
        checks = [item.get(k) for k in ("converged", "pass_cost", "pass_vram", "pass_accuracy")]
        item["verdict"] = ("переносить" if all(v is True for v in checks)
                           else "не переносить" if any(v is False for v in checks) else "не хватает измерений")
        out["gpu"].append(item)
    return out


def compute(bench_path: Path, summary_path: Path | None, compare_path: Path | None) -> dict:
    bench = tomllib.loads(bench_path.read_text(encoding="utf-8"))
    base = bench_path.parent
    tar = json.loads((base / bench["meta"]["tariffs"]).read_text(encoding="utf-8"))
    meshes = json.loads((base / bench["meta"]["mesh_stats"]).read_text(encoding="utf-8"))
    region, fb = bench["meta"]["region"], bench["meta"]["fallback_region"]
    dec, be, mm = bench["decision"], bench["break_even"], bench["memory_model"]
    extras_h = extras_usd_h(tar, region, dec.get("disk_gb"))

    pairs = []
    for g in be["gpu"]:
        for c in be["cpu"]:
            for mg, mc in be["markets"]:
                pg, pc = price(tar, region, g, mg), price(tar, region, c, mc)
                item = {"gpu": g, "gpu_market": mg, "cpu": c, "cpu_market": mc,
                        "gpu_usd_h": pg, "cpu_usd_h": pc}
                if pg is None or pc is None:
                    item["s_star"] = None
                    item["note"] = "нет цены в тарифах"
                else:
                    s_star = break_even_speedup(pg, pc)
                    item["s_star"] = round(s_star, 4)
                    item["s_required"] = round(dec["speedup_margin"] * s_star, 4)
                    if extras_h is not None:
                        item["extras_usd_h"] = round(extras_h, 6)
                        item["s_star_with_extras"] = round(break_even_speedup(pg + extras_h, pc + extras_h), 4)
                    grid = []
                    for r in dec["amdahl_r_host"]:
                        for sk in dec["amdahl_kernel_speedup"]:
                            fr = amdahl_required_fraction(dec["speedup_margin"] * s_star, sk, r)
                            grid.append({"r_host": r, "kernel_speedup": "inf" if sk == 0 else sk,
                                         "f_required": None if fr is None else round(fr, 4)})
                    item["amdahl_f_required"] = grid
                pairs.append(item)

    # Сверка с цифрами задания (eu-north-1, по требованию: 4,6 и 1,5).
    task = []
    for g, claim in (("p5.4xlarge", 4.6), ("g6e.2xlarge", 1.5)):
        pg, pc = price(tar, fb, g, "od"), price(tar, fb, "c8i.8xlarge", "od")
        if pg and pc:
            task.append({"gpu": g, "cpu": "c8i.8xlarge", "region": fb, "gpu_usd_h": pg, "cpu_usd_h": pc,
                         "s_star": round(pg / pc, 4), "task_claim": claim,
                         "match_1_decimal": round(pg / pc, 1) == claim})

    # Видеопамять и упор в полосу памяти — по сеткам узла geom.
    mem, roof = [], []
    for m in meshes["meshes"]:
        epp = m["edges_per_point"]
        for prec, sb in (("fp64", 8), ("mixed", 4)):
            flow = memory_per_point(epp, mm["nvar_flow"], sb, mm["krylov_restart"], mm["extra_vectors"], mm["index_bytes"])
            turb = memory_per_point(epp, mm["nvar_turb"], sb, mm["krylov_restart"], mm["extra_vectors"], mm["index_bytes"])
            bpp = flow["total"] + turb["total"]
            # Обмен за нелинейную итерацию (SU2 v8.6.0): матрица каждой системы идёт на GPU один раз за
            # линейное решение (CSysMatrixVectorProduct → HtDTransfer), правая часть и решение — туда,
            # решение — обратно (CSysSolve::HandleTemporariesIn/Out). Без предобуславливателей на хосте.
            nv_sum = mm["nvar_flow"] + mm["nvar_turb"]
            htod = (flow["matrix"] + turb["matrix"] + 2 * nv_sum * sb) * m["n_points"]
            dtoh = nv_sum * sb * m["n_points"]
            item = {"mesh": m["file"], "n_points": m["n_points"], "edges_per_point": epp, "precision": prec,
                    "bytes_per_point_flow": round(flow["total"], 1), "bytes_per_point_turb": round(turb["total"], 1),
                    "bytes_per_point": round(bpp, 1), "device_GB_this_mesh": round(bpp * m["n_points"] / 1e9, 4),
                    "htod_GB_per_iter": round(htod / 1e9, 4), "dtoh_GB_per_iter": round(dtoh / 1e9, 5)}
            for itype, vram in dec["gpu_vram_bytes"].items():
                item[f"max_points_{itype}"] = int(dec["vram_fill_max"] * vram / bpp)
            mem.append(item)
        for itype, spec in bench["gpu_specs"].items():
            for prec, sb, peak in (("fp64", 8, spec["fp64_flops"]), ("mixed", 4, spec["fp32_flops"])):
                best, worst = spmv_intensity(epp, mm["nvar_flow"], sb, mm["index_bytes"])
                need = best * spec["mem_bw_Bps"]
                roof.append({"mesh": m["file"], "gpu": itype, "precision": prec,
                             "intensity_flop_per_B": [round(worst, 4), round(best, 4)],
                             "flops_needed_at_bw": round(need / 1e9, 1), "flops_peak_G": round(peak / 1e9, 1),
                             "memory_bound": need < peak})
    bw = []
    for itype, spec in bench["gpu_specs"].items():
        for mkt in ("od", "spot"):
            p = price(tar, region, itype, mkt)
            if p:
                bw.append({"gpu": itype, "market": mkt, "usd_h": p,
                           "TBps_per_usd_h": round(spec["mem_bw_Bps"] / 1e12 / p, 4)})

    res = {
        "method": ("Порог окупаемости по стоимости прогона до сходимости: GPU выгоден, если ускорение "
                   "S = T_CPU / T_GPU ≥ speedup_margin · S*, S* = тариф GPU / тариф CPU (тариф — долл./ч машины). "
                   "Доля ядра — закон Амдала с поправкой на CPU GPU-машины: S = 1 / ((1 − f)·r_host + f/s_k + x). "
                   "Видеопамять — нижняя граница по структурам линейного решателя SU2 v8.6.0 (матрица LDU, копия "
                   "ILU(0), векторы FGMRES) для систем течения и турбулентности, рёбра — подсчёт по сеткам узла geom. "
                   "Упор в полосу памяти — интенсивность блочного SpMV против пиковой производительности."),
        "software": f"python3 {platform.python_version()}, break_even.py",
        "date": DATE,
        "inputs": {
            "bench": {"file": rel(bench_path), "sha256": sha256(bench_path)},
            "tariffs": {"file": bench["meta"]["tariffs"], "sha256": sha256(base / bench["meta"]["tariffs"]),
                        "price_list": [i for i in tar["inputs"]]},
            "mesh_stats": {"file": bench["meta"]["mesh_stats"], "sha256": sha256(base / bench["meta"]["mesh_stats"])},
            "region": region, "speedup_margin": dec["speedup_margin"], "vram_fill_max": dec["vram_fill_max"],
            "disk_gb": dec.get("disk_gb"), "extras_usd_h": None if extras_h is None else round(extras_h, 6),
            "memory_model": mm, "gpu_specs": bench["gpu_specs"],
        },
        "break_even": pairs,
        "task_crosscheck": task,
        "missing_in_region": tar["missing"].get(region, []),
        "gpu_memory_estimate": mem,
        "spmv_roofline": roof,
        "bandwidth_per_dollar": bw,
    }
    if summary_path:
        cmp = load_summary(compare_path) if compare_path else None
        res["inputs"]["summary"] = {"file": rel(summary_path), "sha256": sha256(summary_path)}
        res["measured"] = measured_verdicts(bench, load_summary(summary_path), cmp)
    else:
        res["measured"] = {"status": "измерений нет: прогонов по протоколу не было (14.11.2026)"}
    return res


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    t = sub.add_parser("tariffs", help="цены из локальных файлов Price List и Spot Advisor")
    t.add_argument("--price-list", action="append", required=True, metavar="REGION=PATH")
    t.add_argument("--spot-advisor", required=True)
    t.add_argument("--vpc-price", action="append", default=[], metavar="REGION=PATH",
                   help="JSON offer AmazonVPC региона (цена публичного IPv4)")
    t.add_argument("--types", default="", help="типы через запятую (по умолчанию — набор протокола)")
    t.add_argument("--out", required=True)
    m = sub.add_parser("mesh", help="узлы и рёбра сеток .su2")
    m.add_argument("meshes", nargs="+")
    m.add_argument("--out", required=True)
    c = sub.add_parser("compute", help="break_even.json")
    c.add_argument("--bench", default=str(HERE / "bench.toml"))
    c.add_argument("--summary", default=None, help="bench_summary.csv из bench_parse.py")
    c.add_argument("--compare", default=None, help="bench_compare.csv из bench_parse.py")
    c.add_argument("--out", required=True)
    a = ap.parse_args(argv)
    try:
        if a.cmd == "tariffs":
            data = tariffs(a)
        elif a.cmd == "mesh":
            import numpy
            data = {"method": ("Подсчёт узлов, объёмных элементов и уникальных рёбер сетки SU2 (рёбра элементов: "
                               "тетраэдр 6, призма 9, пирамида 8, гексаэдр 12). Число внедиагональных блоков "
                               "матрицы SU2 = 2 × рёбра (без узлов-призраков разбиения и периодики)."),
                    "software": f"python3 {platform.python_version()}, numpy {numpy.__version__}",
                    "date": DATE, "inputs": [], "meshes": []}
            for p in a.meshes:
                st = mesh_stats(Path(p))
                data["inputs"].append(st["file"])
                data["meshes"].append(st)
        else:
            data = compute(Path(a.bench), Path(a.summary) if a.summary else None,
                           Path(a.compare) if a.compare else None)
    except BenchError as e:
        print(f"ошибка: {e}", file=sys.stderr)
        return 2
    Path(a.out).write_text(json.dumps(data, ensure_ascii=False, indent=1) + "\n", encoding="utf-8")
    print(f"записано: {a.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
