#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
routes/skid_twin.py — цифровой двойник скида: одна параметрическая модель
вместо трёх markdown-файлов. Крутилки: сырьё (факел/труба), продукт-вариант
(этилен-I4 / MSA / ПАО / урея-синергия), выход/конверсия, цены, hashprice
(альтернатива на факеле), карбон-кредит. Выход: выручка/EBITDA/окупаемость/
множитель к майнингу + чувствительности.

Все константы — из наших же посчитанных файлов (экономика, YSZ, кинетический
множитель); скрин-уровень ±2–3×, честно. Запуск:
  python3 routes/skid_twin.py                # базовая таблица сценариев
  python3 routes/skid_twin.py --sens         # + чувствительности
"""
import json
import math
import os
import sys

DIR = os.path.dirname(os.path.abspath(__file__))
R_KCAL = 1.987e-3

# ------------------------------------------------------------ базовые крутилки
GAS = {
    "flare": {"gas_usd_t": 30.0,   "note": "факел: $0-100/т, часто отрицательно"},
    "pipeline": {"gas_usd_t": 180.0, "note": "труба: ~$3.5/MMBtu"},
}
CH4_T_Y = 7300.0          # т/год на 1 MMscf/д

# продукт-варианты: (выход т продукта/т CH4 при базовой конверсии,
#  цена $/т, опекс $/т продукта, капекс скида $M, логистика)
PRODUCTS = {
    "ethylene_I4": dict(
        yield_t_per_t=0.35 * 28 / 32.08,   # Y=35% моль C2H4 на 2 CH4
        price=900, opex_per_t=180, capex_musd=3.2,
        logistics="криоцистерна/труба — ТОЛЬКО при близком потребителе",
        note="Y=35% НЕ доказан (риск A1: центр даёт 1.4x, нужен стек рычагов); "
             "при Y=25% выход x0.71"),
    "MSA": dict(
        # честная база: полная конверсия CH4->MSA ~2.5% (лимит: серный узел,
        # рецикл, per-pass) ≈ 1.1 кт/г — как в HIGHVALUE-оценке; 20% = upside
        yield_t_per_t=0.025 * 96.1 / 16.04,
        price=2800, opex_per_t=830, capex_musd=4.5,   # опекс вкл. серу ~$45/т MSA
        logistics="хим-автоцистерна (нерж.) — отлично",
        note="рынок ~50-70 кт/г: ниша 10-30 скидов (риск B1); "
             "upside x8 при конверсии 20% — и тогда рынок насыщается 6 скидами"),
    "PAO": dict(
        yield_t_per_t=0.35 * 28 / 32.08 * 0.85,  # этилен -> олигомеризация 85%
        price=3200, opex_per_t=650, capex_musd=5.0,
        logistics="бочки/цистерна при 25C — отлично",
        note="надстройка над I4; дескриптор селективности в работе"),
    "urea_synergy": dict(
        yield_t_per_t=1.10,   # пиролиз C (0.75) + H2->NH3->урея из CO2 (~0.35+)
        price=380, opex_per_t=140, capex_musd=6.0,
        logistics="навалом/мешки — идеально",
        note="двойной поезд (пиролиз+NH3) — капекс выше; C и урея считаются вместе"),
}

MINING = {"hashprice_usd_ph_day": 50.0, "net_musd_y": 0.8,
          "note": "чистая майнинга при $50/PH/д; ±5x с BTC"}
CARBON = {"usd_t_co2": 0.0, "t_co2_per_skid_y": 20000,
          "note": "0 в базовом кейсе (методология не зачтена — риск E4)"}


def run(product, gas_src, y_mult=1.0, price_mult=1.0, capex_mult=1.0):
    p = PRODUCTS[product]
    g = GAS[gas_src]
    prod_t = CH4_T_Y * p["yield_t_per_t"] * y_mult
    revenue = prod_t * p["price"] * price_mult / 1e6
    gas_cost = CH4_T_Y * g["gas_usd_t"] / 1e6
    opex = prod_t * p["opex_per_t"] / 1e6 + 0.15   # +фикс. O&M
    carbon = CARBON["usd_t_co2"] * CARBON["t_co2_per_skid_y"] / 1e6
    ebitda = revenue - gas_cost - opex + carbon
    capex = p["capex_musd"] * capex_mult
    payback = capex / ebitda if ebitda > 0 else float("inf")
    return dict(product=product, gas=gas_src, prod_t_y=round(prod_t, 0),
                revenue_musd=round(revenue, 2), ebitda_musd=round(ebitda, 2),
                capex_musd=round(capex, 2),
                payback_y=round(payback, 1) if payback != float("inf") else None,
                x_mining=round(ebitda / MINING["net_musd_y"], 1))


def main():
    rows = []
    for prod in PRODUCTS:
        for gas in ("flare", "pipeline"):
            rows.append(run(prod, gas))
    print(f"{'продукт':<14}{'газ':<10}{'т/г':>7}{'выручка':>9}{'EBITDA':>8}"
          f"{'капекс':>8}{'окуп.':>7}{'×майн.':>8}")
    for r in rows:
        pb = f"{r['payback_y']}г" if r['payback_y'] else "∞"
        print(f"{r['product']:<14}{r['gas']:<10}{r['prod_t_y']:>7.0f}"
              f"{r['revenue_musd']:>8.2f}M{r['ebitda_musd']:>7.2f}M"
              f"{r['capex_musd']:>7.2f}M{pb:>7}{r['x_mining']:>8}")
    out = {"scenarios": rows, "knobs": {"GAS": GAS, "PRODUCTS": {
        k: {kk: vv for kk, vv in v.items()} for k, v in PRODUCTS.items()},
        "MINING": MINING, "CARBON": CARBON},
        "honesty": "скрин-уровень ±2-3x; Y=35% этилена НЕ доказан (риск A1); "
                   "рынок MSA — ниша 10-30 скидов (риск B1)"}
    if "--sens" in sys.argv:
        sens = {}
        for prod in PRODUCTS:
            base = run(prod, "flare")["ebitda_musd"]
            sens[prod] = {
                "yield-30%": run(prod, "flare", y_mult=0.7)["ebitda_musd"],
                "price-30%": run(prod, "flare", price_mult=0.7)["ebitda_musd"],
                "capex+50%_payback": run(prod, "flare",
                                         capex_mult=1.5)["payback_y"],
                "base": base}
        out["sensitivity"] = sens
        print("\nчувствительности (EBITDA $M/г, факел):")
        for prod, s in sens.items():
            print(f"  {prod:<14} base {s['base']:>5} | yield-30% {s['yield-30%']:>5} "
                  f"| price-30% {s['price-30%']:>5} | окуп. при капекс+50%: "
                  f"{s['capex+50%_payback']}г")
    with open(os.path.join(DIR, "skid_twin.json"), "w") as f:
        json.dump(out, f, indent=1, ensure_ascii=False)
    print("\nwrote skid_twin.json")


if __name__ == "__main__":
    main()
