#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
routes/econ_shock_candidates.py — скрин-экономика «шок-кандидатов» S1-S3
(NaCN на руднике / CH3SH из кислого газа / нитрометан). Та же рамка, что
skid_twin: ±2-3x, все допущения явные, объективные флаги на TRL.
"""
import json
import os

DIR = os.path.dirname(os.path.abspath(__file__))
FIX_SKID = 380.0        # $/т фикс малого масштаба (рудник/ГПЗ — обслуживаемые)
WACC = 0.12

def econ(name, t_y, price, var, capex_musd, notes):
    cash = var + FIX_SKID * (1000.0 / t_y)          # фикс нормирован на 1 кт/г
    full = cash + capex_musd * 1e6 * WACC / t_y
    ebitda = (price - full) * t_y / 1e6
    return {"product": name, "t_per_year": t_y, "price_usd_t": price,
            "variable_usd_t": round(var, 0), "full_cost_usd_t": round(full, 0),
            "margin_usd_t": round(price - full, 0),
            "ebitda_musd_y": round(ebitda, 2),
            "payback_y": round(capex_musd / ebitda, 1) if ebitda > 0 else None,
            "notes": notes}

OUT = {"frame": "скрин ±2-3x; фикс $380/т@1кт/г нормируется по тоннажу; WACC 12%",
       "candidates": []}

# S1: NaCN на золотом руднике. HCN Андруссов (выход ~65%), NH3 и NaOH завозятся
# (аммиак возят рутинно), газ скважинный. На т NaCN: 0.50 т CH4, 0.50 т NH3,
# 0.82 т NaOH. Продаём по ДОСТАВЛЕННОЙ цене рудника (2800-3500) — в этом шок.
var_s1 = 0.50 * 30 + 0.50 * 600 + 0.82 * 400 + 180   # газ+NH3+NaOH+энергия/Pt-аморт
OUT["candidates"].append(econ(
    "S1 NaCN@рудник", 2500, 3000, var_s1, 6.0,
    "выход 65%; Pt-сетки в util; рудник-якорь 3-15 кт/г = 1-6 скидов; "
    "конкурент платит доставку+конвой опасного груза; TRL химии высокий "
    "(Андруссов 100 лет), новизна — масштаб/сайтинг"))

# S2: CH3SH из кислого газа (CH4+H2S -> CH3SH+H2). На т: 0.33 т CH4, 0.71 т H2S
# (оба ~бесплатны на кислом месторождении; кредит за утилизацию H2S возможен).
var_s2 = 0.33 * 10 + 0.71 * 0 + 260                   # энергия/катализ/рецикл
OUT["candidates"].append(econ(
    "S2 CH3SH из кислого газа", 1500, 1500, var_s2, 5.0,
    "Ограничения: химия CH4+H2S каталитически TRL 2-3 (термодинамика в гору, "
    "нужен сдвиг равновесия H2-отводом) — самый спекулятивный из трёх; "
    "апсайд: поле может ПЛАТИТЬ за приём H2S (+$50-140/т кредит не заложен); "
    "рынок-якорь метионин 1.5+ Мт/г"))

# S3: нитрометан CH4+HNO3 (парофаз). На т: 0.30 т CH4, 1.15 т HNO3 (селект. 40%
# классики -> наш рычаг; закладываем 55% как цель). Цена 2000.
var_s3 = 0.30 * 30 + 1.15 * 300 + 240
OUT["candidates"].append(econ(
    "S3 нитрометан", 1200, 2000, var_s3, 4.0,
    "классика: селективность ~40%, полвека без починки; наш NEVPT2-рычаг "
    "радикальной цепи *NO2/*CH3; рынок ~0.1 Мт/г — ниша больше MSA"))

with open(os.path.join(DIR, "econ_shock_candidates.json"), "w") as f:
    json.dump(OUT, f, indent=1, ensure_ascii=False)
print(json.dumps(OUT, indent=1, ensure_ascii=False))
