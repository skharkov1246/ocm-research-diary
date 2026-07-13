#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
routes/urea_synergy_model.py — вариант 3 гребёнки (урея-синергия):
себестоимость тонны карбамида на скиде «пиролиз CH4 → C + H2 → NH3 →
(+CO2) → урея». Идея: факельный газ → удобрение на месте (локальные
агрорынки), CO2 ПОТРЕБЛЯЕТСЯ (0.73 т/т), твёрдый углерод — соль-продукт
пиро-волкана (наш NEVPT2-скрин сплавов). Скрин-уровень ±30-50%,
допущения в JSON.

Стехиометрия (массовая, на 1 т уреи CO(NH2)2, M=60.06):
  NH3 0.567 т; CO2 0.733 т (потребляется);
  H2 на NH3 (Габер, 3H2+N2): 0.178 т/т NH3 -> 0.101 т H2;
  CH4 на H2 (пиролиз CH4->C+2H2): 3.979 т/т H2 -> 0.402 т CH4;
  углерод со-продукт: 2.979 т C/т H2 -> 0.301 т C.
"""
import json
import os

DIR = os.path.dirname(os.path.abspath(__file__))

# --- стехиометрия ---
NH3_T = 2 * 17.031 / 60.055
CO2_T = 44.009 / 60.055
H2_PER_NH3 = 3 * 2.016 / (2 * 17.031)
H2_T = NH3_T * H2_PER_NH3
CH4_PER_H2 = 16.043 / (2 * 2.016)
CH4_T = H2_T * CH4_PER_H2
C_T = H2_T * 12.011 / (2 * 2.016) * (16.043 / 16.043)

# --- допущения ($/т, скрин) ---
CO2_USD = 30.0        # захват с дымтрубы соседа / покупка; может быть и 0
UTIL_PYRO = 90.0      # энергия пиролиза (эндотермика 74.9 кДж/моль CH4) + печь
UTIL_HABER = 110.0    # мини-Габер: компрессия/синтез-контур (эл-во на месте)
UTIL_UREA = 80.0      # карбаматный контур 150 бар + грануляция
FIXED_SKID = 350.0    # труд/обслуживание малого масштаба (влезает в 10 кт/г)
CAPEX = 18e6          # мини-NH3 (болевая точка!) + урея-узел + пиролиз
T_Y = 10000.0
CAPITAL = CAPEX * 0.12 / T_Y

UREA_USD = 350.0      # мировая цена 300-450
CARBON_USD = 400.0    # пиро-углерод НЕизвестного качества: 0 (отвал)
                      # ... 1500+ (carbon black grade) — главный рычаг


def cost(gas_usd_t, carbon_usd=CARBON_USD, co2_usd=CO2_USD):
    var = (CH4_T * gas_usd_t + CO2_T * co2_usd
           + UTIL_PYRO + UTIL_HABER + UTIL_UREA)
    credit = C_T * carbon_usd
    cash = var + FIXED_SKID - credit
    full = cash + CAPITAL
    return dict(variable=round(var, 0), carbon_credit=round(credit, 0),
                cash=round(cash, 0), full=round(full, 0),
                margin_vs_urea=round(UREA_USD - full, 0))


OUT = {
 "stoichiometry_per_t_urea": {
     "NH3_t": round(NH3_T, 3), "CO2_consumed_t": round(CO2_T, 3),
     "H2_t": round(H2_T, 3), "CH4_t": round(CH4_T, 3),
     "carbon_coproduct_t": round(C_T, 3)},
 "assumptions_usd": {"co2": CO2_USD, "util_pyro": UTIL_PYRO,
                     "util_haber": UTIL_HABER, "util_urea": UTIL_UREA,
                     "fixed_skid": FIXED_SKID, "capital": round(CAPITAL, 0),
                     "capex_note": "капекс $18M (мини-NH3 — болевая точка), "
                                   "12%/г, 10 кт/г", "urea_price": UREA_USD},
 "skid_flare_carbon400": cost(30.0),
 "skid_pipeline_carbon400": cost(180.0),
 "skid_flare_carbon0_dump": cost(30.0, carbon_usd=0.0),
 "skid_flare_carbon1500_cb": cost(30.0, carbon_usd=1500.0),
 "skid_flare_cb1500_plus_co2credit50": dict(
     cost(30.0, carbon_usd=1500.0),
     co2_credit_usd=round(2.33 * 50.0, 0),
     margin_with_credit=round(UREA_USD - cost(30.0, carbon_usd=1500.0)["full"]
                              + 2.33 * 50.0, 0),
     note="апсайд: 2.33 т CO2e/т × $50 — только с сертификацией"),
 "benchmarks": {
   "worldscale_urea_SMR": {"full_est": [180, 280],
     "note": "большой завод на паровой конверсии: дешевле нас по $/т, но "
             "+0.9 т CO2/т NH3 выбросов; наш маршрут CO2-НЕГАТИВЕН "
             "(0.73 т потребляется, пиролиз без CO2)"},
   "green_urea_electrolysis": {"full_est": [500, 900],
     "note": "электролизный H2: дороже нас при факельном газе"}},
 "co2_balance_per_t": {
   "consumed_urea": -round(CO2_T, 2),
   "smr_reference_emission": +1.6,
   "note": "против SMR-уреи скид экономит ~2.3 т CO2/т — кандидат на "
           "carbon-кредиты (не считаем в марже, апсайд)"},
 "honesty": "скрин ±30-50%. Болевые точки: (1) мини-Габер капексоёмок и "
            "любит масштаб — интегральные мини-NH3 скиды (Syzygy/Talus-тип) "
            "существуют, но дороги; (2) качество пиро-углерода из расплава "
            "НЕ подтверждено (наш пиро-трек считает связывание C, не "
            "морфологию) — вилка $0-1500/т ломает/делает экономику; "
            "(3) урея-узел 150 бар — не 'скид' в строгом смысле; (4) CO2 "
            "на факельной площадке ещё найти надо (дымтруба компрессора). "
            "Квант-часть варианта: газофазный якорь NH3+CO2 (urea_chain) — "
            "отрицательный контроль метода M (закрытые оболочки)."}

if __name__ == "__main__":
    with open(os.path.join(DIR, "urea_synergy_model.json"), "w") as f:
        json.dump(OUT, f, indent=1, ensure_ascii=False)
    print(json.dumps(OUT, indent=1, ensure_ascii=False))
