#!/usr/bin/env python3
"""Тарифы AWS для сметы прогонов ГТУ из официального AWS Price List Bulk API -> aws_rates.toml.

Ключи AWS не нужны: файлы предложений публичны (https://pricing.us-east-1.amazonaws.com/offers/v1.0/...).
Страницы цен aws.amazon.com/ec2/pricing из контейнера сессии закрыты прокси (403, 14.11.2026), а Price List —
открыт; это машинный источник тех же тарифов с датой публикации каждого файла.

Что берётся (регионы us-east-1 и eu-north-1):
  EC2      — Linux, Shared, без ПО, On-Demand, за час, для типов из INSTANCE_TYPES (+ vCPU, память, процессор);
  EBS      — gp3, долл. за ГБ-месяц (AmazonEC2, EBS:VolumeUsage.gp3);
  VPC      — публичный IPv4 в работе, долл./ч (AmazonVPC, PublicIPv4:InUseAddress);
  S3       — Standard, первые 50 ТБ, ГБ-месяц; запросы PUT/COPY/POST/LIST (Tier1) и GET (Tier2), за запрос;
  CloudWatch Logs — приём в стандартный класс (DataProcessing-Bytes), хранение (TimedStorage-ByteHrs);
  Передача — в интернет первые 10 ТБ (после 100 ГБ бесплатных в месяц на весь аккаунт), между регионами;
  ECR      — хранение образов, ГБ-месяц;
  Спот     — Spot Instance Advisor (spot-bid-advisor.s3.amazonaws.com/spot-advisor-data.json, тот же файл
             читает страница советника): средняя экономия к цене по требованию, % и частота прерываний.
             Это не цена на дату запуска: живую цену читает runner.py (describe_spot_price_history) в Coldun.

Запуск:  python3 fetch_rates.py [--cache DIR] [--out aws_rates.toml]
  --cache DIR  — сохранить/взять скачанные файлы из DIR (EC2 CSV — 300 МБ, читается потоком с отбором строк).
"""
from __future__ import annotations

import argparse
import csv
import datetime as dt
import io
import json
import os
import sys
import urllib.request

BASE = "https://pricing.us-east-1.amazonaws.com"
SPOT_ADVISOR = "https://spot-bid-advisor.s3.amazonaws.com/spot-advisor-data.json"
REGIONS = ("us-east-1", "eu-north-1")
INSTANCE_TYPES = ("c7a.8xlarge", "c7a.16xlarge", "c7i.8xlarge", "c7i.16xlarge",
                  "c8i.2xlarge", "c8i.8xlarge", "c8i.16xlarge")
REGION_CODE = {"us-east-1": "USE1", "eu-north-1": "EUN1"}


# ---------------------------------------------------------------- загрузка
def _get(url: str, timeout: int = 120):
    req = urllib.request.Request(url, headers={"User-Agent": "gtu-cost-estimate/1"})
    return urllib.request.urlopen(req, timeout=timeout)  # noqa: S310 (только https-адреса выше)


def region_index(service: str, cache: str | None) -> dict:
    path = os.path.join(cache, f"{service}_region_index.json") if cache else None
    if path and os.path.exists(path):
        return json.load(open(path))
    with _get(f"{BASE}/offers/v1.0/aws/{service}/current/region_index.json") as r:
        data = json.load(r)
    if path:
        json.dump(data, open(path, "w"))
    return data


def offer_json(service: str, region: str, cache: str | None) -> tuple[dict, str]:
    url = BASE + region_index(service, cache)["regions"][region]["currentVersionUrl"]
    path = os.path.join(cache, f"{service}_{region}.json") if cache else None
    if path and os.path.exists(path):
        return json.load(open(path)), url
    with _get(url, timeout=600) as r:
        data = json.load(r)
    if path:
        json.dump(data, open(path, "w"))
    return data, url


def ec2_rows(region: str, cache: str | None) -> tuple[list[dict], dict, str]:
    """Строки EC2 CSV с нужными типами и gp3. CSV читается потоком: первые 5 строк — метаданные, 6-я — шапка."""
    url = BASE + region_index("AmazonEC2", cache)["regions"][region]["currentVersionUrl"].replace(".json", ".csv")
    head_p = os.path.join(cache, f"ec2_{region}_head.csv") if cache else None
    body_p = os.path.join(cache, f"ec2_{region}_filtered.csv") if cache else None
    if head_p and os.path.exists(head_p) and os.path.exists(body_p):
        head = list(csv.reader(open(head_p, newline="")))
        body = list(csv.reader(open(body_p, newline="")))
    else:
        head, body = [], []
        with _get(url, timeout=1800) as r:
            for i, line in enumerate(io.TextIOWrapper(r, encoding="utf-8", newline="")):
                if i < 6:
                    head.append(next(csv.reader([line])))
                elif any(t + '"' in line for t in INSTANCE_TYPES) or "VolumeUsage.gp3" in line:
                    body.append(next(csv.reader([line])))
        if head_p:
            csv.writer(open(head_p, "w", newline="")).writerows(head)
            csv.writer(open(body_p, "w", newline="")).writerows(body)
    meta = {row[0]: row[1] for row in head[:5] if len(row) >= 2}
    cols = head[5]
    return [dict(zip(cols, row)) for row in body], meta, url


def spot_advisor(cache: str | None) -> tuple[dict, str]:
    path = os.path.join(cache, "spot-advisor-data.json") if cache else None
    if path and os.path.exists(path):
        data = json.load(open(path))
        lm = data.get("_last_modified", "")
    else:
        with _get(SPOT_ADVISOR) as r:
            data = json.load(r)
            lm = r.headers.get("Last-Modified", "")
        data["_last_modified"] = lm
        if path:
            json.dump(data, open(path, "w"))
    return data, lm


# ---------------------------------------------------------------- отбор тарифов
def ondemand_terms(offer: dict) -> dict:
    out = {}
    for sku, offers in offer["terms"].get("OnDemand", {}).items():
        for o in offers.values():
            for pd in o["priceDimensions"].values():
                out.setdefault(sku, []).append(pd)
    return out


def pick(offer: dict, pred, begin: str = "0") -> dict:
    """Единственная ценовая ступень с beginRange == begin у продукта, подходящего под pred."""
    terms = ondemand_terms(offer)
    found = []
    for sku, p in offer["products"].items():
        if pred(p.get("attributes", {})):
            for pd in terms.get(sku, []):
                if pd.get("beginRange", "0") == begin:
                    found.append((sku, pd))
    if len(found) != 1:
        raise ValueError(f"ожидалась одна строка тарифа, найдено {len(found)}")
    sku, pd = found[0]
    return {"usd": float(pd["pricePerUnit"]["USD"]), "unit": pd["unit"], "sku": sku,
            "description": pd["description"]}


def extract_ec2(rows: list[dict]) -> tuple[dict, dict]:
    inst, gp3 = {}, None
    for x in rows:
        t = x.get("Instance Type", "")
        if (t in INSTANCE_TYPES and x["TermType"] == "OnDemand" and x["Tenancy"] == "Shared"
                and x["Operating System"] == "Linux" and x["Pre Installed S/W"] == "NA"
                and x["CapacityStatus"] == "Used" and x["License Model"] == "No License required"
                and x["Unit"] == "Hrs"):
            if t in inst:
                raise ValueError(f"две строки цены для {t}")
            inst[t] = {"on_demand_usd_h": float(x["PricePerUnit"]), "vcpu": int(x["vCPU"]),
                       "memory_gib": float(x["Memory"].split()[0]), "processor": x["Physical Processor"],
                       "sku": x["SKU"], "effective": x["EffectiveDate"]}
        if x.get("usageType", "").endswith("EBS:VolumeUsage.gp3") and x["TermType"] == "OnDemand":
            gp3 = {"usd": float(x["PricePerUnit"]), "unit": x["Unit"], "sku": x["SKU"],
                   "description": x["PriceDescription"]}
    missing = set(INSTANCE_TYPES) - set(inst)
    if missing:
        raise ValueError(f"нет цены для {sorted(missing)}")
    if gp3 is None:
        raise ValueError("нет цены gp3")
    return inst, gp3


def extract_region(region: str, cache: str | None) -> tuple[dict, dict]:
    code = REGION_CODE[region]
    sources = {}
    rows, meta, url = ec2_rows(region, cache)
    inst, gp3 = extract_ec2(rows)
    sources[f"AmazonEC2_{region}"] = {"url": url, "publication": meta.get("Publication Date", ""),
                                      "version": meta.get("Version", "")}

    def offer(service):
        data, u = offer_json(service, region, cache)
        sources[f"{service}_{region}"] = {"url": u, "publication": data.get("publicationDate", ""),
                                          "version": data.get("version", "")}
        return data

    vpc = offer("AmazonVPC")
    ipv4 = pick(vpc, lambda a: a.get("usagetype", "").endswith("PublicIPv4:InUseAddress"))
    s3 = offer("AmazonS3")
    s3_std = pick(s3, lambda a: a.get("usagetype", "").endswith("TimedStorage-ByteHrs")
                  and a.get("volumeType") == "Standard")
    s3_put = pick(s3, lambda a: a.get("usagetype", "") in ("Requests-Tier1", f"{code}-Requests-Tier1"))
    s3_get = pick(s3, lambda a: a.get("usagetype", "") in ("Requests-Tier2", f"{code}-Requests-Tier2"))
    cw = offer("AmazonCloudWatch")
    cw_in = pick(cw, lambda a: a.get("usagetype", "") == f"{code}-DataProcessing-Bytes")
    cw_st = pick(cw, lambda a: a.get("usagetype", "") == f"{code}-TimedStorage-ByteHrs")
    ecr = offer("AmazonECR")
    ecr_st = pick(ecr, lambda a: a.get("usagetype", "") in ("TimedStorage-ByteHrs", f"{code}-TimedStorage-ByteHrs"))
    dt_ = offer("AWSDataTransfer")
    out_inet = pick(dt_, lambda a: a.get("transferType") == "AWS Outbound" and a.get("toLocation") == "External"
                    and not a.get("usagetype", "").startswith("Global"))
    free_tier = pick(dt_, lambda a: a.get("usagetype", "") == "Global-DataTransfer-Out-Bytes")
    other = [r for r in REGIONS if r != region][0]
    inter = pick(dt_, lambda a: a.get("transferType") == "InterRegion Outbound"
                 and a.get("fromRegionCode") == region and a.get("toRegionCode") == other)
    return {
        "ec2": inst,
        "ebs_gp3_usd_gb_month": gp3,
        "public_ipv4_usd_h": ipv4,
        "s3_standard_usd_gb_month": s3_std,
        "s3_put_usd_request": s3_put,
        "s3_get_usd_request": s3_get,
        "cw_logs_ingest_usd_gb": cw_in,
        "cw_logs_storage_usd_gb_month": cw_st,
        "ecr_storage_usd_gb_month": ecr_st,
        "transfer_out_internet_usd_gb": out_inet,
        "transfer_out_free_gb_month_account": {"usd": 0.0, "unit": "GB",
                                               "gb": float(free_tier.get("description", "").split("for ")[1].split("GB")[0])
                                               if "for " in free_tier.get("description", "") else 0.0,
                                               "sku": free_tier["sku"], "description": free_tier["description"]},
        f"transfer_inter_region_to_{other}_usd_gb": inter,
    }, sources


# ---------------------------------------------------------------- запись TOML
def _key(k: str) -> str:
    return k if all(c.isalnum() or c in "_-" for c in k) else json.dumps(k)


def _val(v) -> str:
    if isinstance(v, bool):
        return "true" if v else "false"
    if isinstance(v, (int, float)):
        return repr(v)
    if isinstance(v, str):
        return json.dumps(v, ensure_ascii=False)
    if isinstance(v, list):
        return "[" + ", ".join(_val(x) for x in v) + "]"
    raise TypeError(type(v))


def dump_toml(data: dict, prefix: tuple = ()) -> str:
    """Простой писатель TOML: скаляры и списки скаляров, вложенные таблицы — заголовками [a.b]."""
    lines, subs = [], []
    for k, v in data.items():
        if isinstance(v, dict):
            subs.append((k, v))
        else:
            lines.append(f"{_key(k)} = {_val(v)}")
    out = []
    if prefix and lines:
        out.append("[" + ".".join(_key(p) for p in prefix) + "]")
    out.extend(lines)
    text = "\n".join(out) + ("\n" if out else "")
    for k, v in subs:
        block = dump_toml(v, prefix + (k,))
        text += ("\n" if text else "") + block
    return text


def build(cache: str | None) -> dict:
    spot, spot_lm = spot_advisor(cache)
    ranges = {r["index"]: r["label"] for r in spot["ranges"]}
    regions, sources = {}, {}
    for region in REGIONS:
        rates, src = extract_region(region, cache)
        sources.update(src)
        adv = spot["spot_advisor"].get(region, {}).get("Linux", {})
        for t, rec in rates["ec2"].items():
            a = adv.get(t)
            rec["spot_savings_pct"] = int(a["s"]) if a else -1          # -1: в файле советника нет
            rec["spot_interrupt_bucket"] = ranges.get(a["r"], "") if a else ""
            rec["spot_estimate_usd_h"] = round(rec["on_demand_usd_h"] * (1 - a["s"] / 100), 6) if a else -1.0
        regions[region] = rates
    sources["spot_advisor"] = {"url": SPOT_ADVISOR, "publication": spot_lm, "version": ""}
    return {
        "meta": {
            "generated_utc": dt.datetime.now(dt.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
            "method": "AWS Price List Bulk API (offers/v1.0, публичные файлы предложений) + Spot Instance Advisor",
            "software": "fetch_rates.py, Python " + sys.version.split()[0],
            "note": ("Цены без налогов, USD. spot_estimate_usd_h = on_demand × (1 − экономия/100) по файлу советника "
                     "(средняя экономия, не цена на дату); -1 — нет данных. Страницы aws.amazon.com/*/pricing из "
                     "контейнера закрыты (403), поэтому источник — Price List."),
        },
        "sources": sources,
        "regions": regions,
    }


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--cache", help="папка для скачанных файлов (берутся оттуда, если уже есть)")
    ap.add_argument("--out", default=os.path.join(os.path.dirname(os.path.abspath(__file__)), "aws_rates.toml"))
    a = ap.parse_args(argv)
    if a.cache:
        os.makedirs(a.cache, exist_ok=True)
    data = build(a.cache)
    header = ("# Тарифы AWS для сметы прогонов ГТУ. Сгенерировано fetch_rates.py — не править вручную.\n"
              "# Источник каждого числа — таблица [sources.*] (адрес файла Price List, дата публикации).\n")
    open(a.out, "w").write(header + dump_toml(data))
    print(f"записано {a.out}: {sum(len(r['ec2']) for r in data['regions'].values())} цен EC2, "
          f"{len(data['sources'])} источников")
    return 0


if __name__ == "__main__":
    sys.exit(main())
