#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
routes/msa_cost_model.py — себестоимость тонны MSA: наш скид (факел/труба)
vs классика (меркаптан+Cl2) vs большой прямой завод (Grillo-тип).
Стехиометрия: CH4 + SO3 -> CH3SO3H (96.1): 0.167 т CH4 + 0.833 т SO3
(=0.333 т серы) на тонну. Скрин-уровень ±30-50%, допущения в JSON.
"""
import json
import os

DIR = os.path.dirname(os.path.abspath(__file__))

CH4_T, S_T = 0.167, 0.333          # т на т MSA
SULFUR_USD = 150.0                  # $/т сера (комодити 100-200)
INITIATOR = 100.0                   # $/т MSA (пероксид/электро; болевая точка -> наш рычаг)
UTILITIES = 130.0                   # $/т: SO3-узел+дистилляция+рецикл (часть закрыта теплом серы)
FIXED_SKID = 400.0                  # $/т: труд/обслуживание/накладные МАЛОГО масштаба (1.1 кт/г)
CAPEX_SKID = 4.5e6
T_Y = 1093.0
CAPITAL = CAPEX_SKID * 0.12 / T_Y   # $/т: 12% годовых (аморт.+WACC)

def skid(gas_usd_t):
    var = CH4_T * gas_usd_t + S_T * SULFUR_USD + INITIATOR + UTILITIES
    cash = var + FIXED_SKID
    full = cash + CAPITAL
    return dict(variable=round(var, 0), cash=round(cash, 0),
                full=round(full, 0),
                margin_at_2800=round(2800 - full, 0))

OUT = {
 "per_tonne_inputs": {"CH4_t": CH4_T, "S_t": S_T,
                      "sulfur_usd_t": SULFUR_USD, "initiator_usd": INITIATOR,
                      "utilities_usd": UTILITIES, "fixed_skid_usd": FIXED_SKID,
                      "capital_usd": round(CAPITAL, 0),
                      "capital_note": "капекс $4.5M, 12%/г, 1.1 кт/г"},
 "our_skid_flare": skid(30.0),
 "our_skid_pipeline": skid(180.0),
 "benchmarks": {
   "classic_mercaptan_Cl2": {"full_est": [1200, 1800],
     "note": "метанол->меркаптан->хлорокисление: дорогая цепочка сырья + "
             "стоки HCl/хлорорганика; мировой масштаб"},
   "grillo_direct_worldscale": {"full_est": [700, 1000],
     "note": "прямой CH4+SO3 на большом заводе у трубы: атомная экономия ~100%, "
             "эффект масштаба; оценка по открытым заявлениям, не наш расчёт"}},
 "honesty": "скрин ±30-50%. На трубе против big-Grillo мы НЕ дешевле по "
            "себестоимости (масштаб бьёт); наши преимущества: (1) почти "
            "бесплатный факельный газ, (2) локальные продажи без импортной "
            "логистики, (3) квантовый дизайн инициатора (строка INITIATOR — "
            "прямая цель снижения), (4) конверсия-апсайд x8 размазывает "
            "FIXED+CAPITAL по большему тоннажу: при 20% конверсии full "
            "падает к ~$700/т."}

# апсайд: конверсия 20% (8747 т/г) — фикс и капитал размазываются
T2 = 8747.0
cap2 = CAPEX_SKID * 1.3 * 0.12 / T2       # капекс x1.3 на больший серный узел
var_fl = CH4_T * 30 + S_T * SULFUR_USD + INITIATOR + UTILITIES
OUT["our_skid_flare_20pct_conv"] = {
    "full": round(var_fl + FIXED_SKID * T_Y / T2 + cap2, 0),
    "note": "рынок насыщается 6 скидами — цена поедет; для первых игроков ок"}

with open(os.path.join(DIR, "msa_cost_model.json"), "w") as f:
    json.dump(OUT, f, indent=1, ensure_ascii=False)
print(json.dumps(OUT, indent=1, ensure_ascii=False))
