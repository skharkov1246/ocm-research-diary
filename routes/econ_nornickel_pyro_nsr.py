#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
routes/econ_nornickel_pyro_nsr.py — Норникель N3: пиролиз стренд-газа Норильска —
сажа на СМП (NSR = Northern Sea Route), H2 на месте.

Норильск — чистейший stranded gas РФ: газ есть, трубы наружу нет, СМП есть.
Фильтр логистики портфеля («твёрдое или жидкое >$1500/т») даёт продукты:
 (1) углерод 0.75 т/т CH4 — ЕДИНСТВЕННЫЙ экспортный продукт (биг-бэги ->
     Дудинка -> СМП). Вся экономика сидит в гейте качества: сырая сажа
     $200-400/т НЕ проходит даже майнинг-порог ($323/$751 паритет/3x,
     econ_vs_mining), carbon black $1200-2000/т проходит 3x. Гейт = наш
     соляной слой NaCl (чистота 83 -> <0.1 wt% металла, лит. +
     METHANE_ROUTES_RESEARCH.md);
 (2) H2 0.25 т/т CH4 — ТОЛЬКО на месте. Урок нашего же самопойманного
     МИРАЖА 6.9x (METHANE_INVENTION_ECONOMICS.md): H2 НЕ уезжает с устья без
     трубы/аммиака — здесь он замещает газ в энергетике по теплотворности
     (120 МДж/кг H2 / 50 МДж/кг CH4 = 2.4 т CH4-экв/т H2). При внутреннем
     газе $20-50/т это копейки ($48-120/т H2 = $0.05-0.12/кг): проект —
     САЖЕВЫЙ ЗАВОД, водород — возврат 60% теплотворности газа в энергосистему
     (чистый расход газа 0.4 т-экв на 1 т пиролизованного CH4).

LCOH-параметрика калибрована к той же литературной TEA, что и L3
(econ_lukoil_h2_refinery.py): $1.29-1.53/кг масштаб / $3.06 малый /
углерод $250/т режет ~25% — модель ОБЯЗАНА воспроизводить (anchor_reproduction).
Капекс — из капитальной строки TEA (перпетуитет capital*WACC), интенсивность
$/(т H2/год) помечена E1; северный коэффициент сверху.

Скрин-уровень ±2-3x, НЕ банковская модель. Все данные Норникеля (баланс
Норильскгазпрома, избыток газа, внутренняя цена газа, свободная мощность,
логистика Дудинка-СМП) — company_assumption, проверить у компании.
Запуск: python3 routes/econ_nornickel_pyro_nsr.py
     -> routes/econ_nornickel_pyro_nsr_results.json
"""
import json
import math
import os

DIR = os.path.dirname(os.path.abspath(__file__))

# --------------------------------------------- масс-баланс и физика (жёсткие якоря)
H2_PER_T_CH4 = 0.25            # т H2 / т CH4 (лит. TEA: CH4 -> C + 2H2)
C_PER_T_CH4 = 0.75             # т C  / т CH4
CH4_KG_PER_KG_H2 = 4.0         # = 1/0.25
C_KG_PER_KG_H2 = 3.0           # = 0.75/0.25
CH4_T_PER_MLN_M3 = 16.04 / 22.414 * 1000.0     # 715.7 ~ 716 т CH4 / млн м3
LHV_H2, LHV_CH4 = 120.0, 50.0                  # МДж/кг (низшая теплотворность)
CH4_EQ_PER_T_H2 = LHV_H2 / LHV_CH4             # 2.4 т CH4-экв / т H2


def A(value, unit, tier, comment):
    """Допущение: значение + единица + уровень достоверности + источник."""
    return {"value": value, "unit": unit, "tier": tier, "comment": comment}


# ------------------------------------------------------------ допущения (все явно)
ASSUMPTIONS = {
 # --- физика (твёрже любого якоря) ---
 "ch4_t_per_mln_m3": A(round(CH4_T_PER_MLN_M3, 0), "т CH4/млн м3", "anchor_repo",
    "физика: плотность CH4 16.04/22.414 = 0.716 кг/м3 (н.у.) — 1 млн м3 ~ 716 т; "
    "якорь ТЗ 36/107/286 кт для 50/150/400 млн м3 воспроизводится"),
 "lhv_h2_ch4_mj_kg": A([LHV_H2, LHV_CH4], "МДж/кг", "anchor_repo",
    "физика/справочник: 1 т H2 замещает 2.4 т CH4 по теплу — база H2-кредита"),
 "mass_balance": A({"h2_t_per_t_ch4": H2_PER_T_CH4, "c_t_per_t_ch4": C_PER_T_CH4},
    "т/т CH4", "literature", "лит. TEA пиролиза: CH4 -> 0.25 H2 + 0.75 C"),
 # --- якоря репозитория ---
 "wacc": A(0.12, "1/год", "anchor_repo",
    "econ_shock_candidates / msa_cost_model: 12% годовых на капитал"),
 "asset_life_y": A(15, "лет", "anchor_repo",
    "econ_vs_mining: жизнь пиролизного актива 15 лет (против 12 у OCM-скидов); "
    "NPV через аннуитет"),
 "mmscf_conversion": A(7300.0, "т CH4/год на 1 MMscf/д", "anchor_repo",
    "skid_twin: для пересчёта масштаба в скид-эквиваленты"),
 "mining_floor_usd_t": A([323, 751], "$/т продукта", "anchor_repo",
    "econ_vs_mining_results: паритет / 3x к майнингу; сырая сажа $200-400 НЕ "
    "проходит 3x (а $200 — даже паритет), carbon black $1200-2000 проходит"),
 "cu_bi_nevpt2": A({"Ebind_C_dft": -395.02, "Ebind_C_nevpt2": 0.77},
    "ккал/моль (относит. вулкан)", "anchor_repo",
    "pyrolysis_descriptor_results.json: NEVPT2 переворачивает DFT-ранжир — "
    "Cu-Bi слабейший связыватель C = некоксующийся кандидат; наше IP, гейт G1"),
 "nacl_layer": A("чистота углерода 83 -> <0.1 wt% металла", "-", "anchor_repo",
    "соляной слой NaCl поверх расплава (дескриптор адгезии Baumli/Kaptay: NaCl "
    "66 мДж/м2 — оптимум); лит. + METHANE_ROUTES_RESEARCH.md; ГЕЙТ КАЧЕСТВА — "
    "на нём висит вся разница raw $200-400 vs carbon black $1200-2000"),
 # --- литература (чужая TEA и справочники, не наш расчёт) ---
 "lcoh_scale_usd_kg": A([1.29, 1.53], "$/кг H2", "literature",
    "TEA пиролиза в расплаве, 20 бар / 1100C, мировой масштаб, БЕЗ продажи "
    "углерода — калибровочный якорь параметрики (как в L3)"),
 "lcoh_small_usd_kg": A(3.06, "$/кг H2", "literature",
    "та же TEA, малый масштаб 1-10 т H2/день, без продажи углерода"),
 "carbon_credit_lit": A("углерод $250/т режет LCOH на ~25%; $700/т — DOE $1/кг",
    "-", "literature", "модель обязана воспроизводить (anchor_reproduction)"),
 "elec_kwh_kg_h2": A([13, 17], "кВт*ч/кг H2", "literature",
    "<13 кВт*ч/кг = ~треть электролиза (~40-52) — якорь ТЗ; в TEA тепло от "
    "сжигания части хвостового газа, электрификация — опция; расчёт МВт — "
    "по нижней границе 13"),
 "carbon_price_raw_usd_t": A([200, 400], "$/т", "literature",
    "сырая пиролизная сажа без сертификации (наполнитель/металлургия)"),
 "carbon_price_cb_usd_t": A([1200, 2000], "$/т", "literature",
    "качество carbon black (шинная спецификация) — ТОЛЬКО после гейта качества"),
 "calib_td_small": A(5.0, "т H2/день", "literature",
    "середина лит. диапазона «малый масштаб 1-10 т/д» (как в L3)"),
 "calib_td_scale": A(150.0, "т H2/день", "literature",
    "«мировой масштаб» TEA; между точками — степенная интерполяция non-gas "
    "остатка (та же калибровка, что econ_lukoil_h2_refinery)"),
 "capital_share_nongas": A(0.5, "доля", "literature",
    "доля капитальной строки в non-gas остатке LCOH (в TEA типично 40-60%)"),
 "capex_intensity_e1": A([2400, 9800], "$/(т H2/год)", "literature",
    "E1: капексная интенсивность малых модулей — ПРОИЗВОДНАЯ капитальной "
    "строки лит. TEA через перпетуитет capital*WACC (не вендорная котировка!); "
    "низ — floor масштаба, верх — 5 т/д; северный коэффициент СВЕРХУ; "
    "проверить котировками (Monolith / Modern Hydrogen)"),
 "fired_extra_gas": A([0.10, 0.25], "доля от сырьевого CH4", "literature",
    "огневой обогрев хвостовым газом (как в TEA-якоре): доп. расход газа "
    "~0.5-1.0 кг CH4/кг H2 с рекуперацией — скрин; альтернатива электрификации"),
 # --- рынок ---
 "arctic_capex_mult": A([1.3, 2.0], "x к капексу", "market",
    "северный коэффициент EPC (Норильск: логистика материалов, вечная "
    "мерзлота, зимний монтаж) — индустриальная практика; уточнить по "
    "нормативам компании"),
 "rf_carbon_black_kt_y": A([700, 1100], "кт/год", "market",
    "производство техуглерода РФ, порядок величины (как в L3); внутренний "
    "спрос шинников МЕНЬШЕ — при наших объёмах экспорт по СМП обязателен"),
 # --- данные компании: ТОЛЬКО диапазоны, ПРОВЕРИТЬ У КОМПАНИИ ---
 "ngp_balance_bln_m3_y": A([3, 4], "млрд м3/год", "company_assumption",
    "добыча Норильскгазпрома на энергетику Норильского узла — проверить у "
    "компании (публичный порядок величины)"),
 "surplus_mln_m3_y": A({"pessimistic": 50, "base": 150, "optimistic": 400},
    "млн м3/год", "company_assumption",
    "доступный для пиролиза ИЗБЫТОК сверх нужд энергетики = 36/107/286 кт "
    "CH4/год — ключевое непроверенное допущение, проверить у компании"),
 "gas_internal_usd_t": A([20, 50], "$/т CH4", "company_assumption",
    "внутренняя себестоимость газа Норильскгазпрома (своя добыча, "
    "изолированная система, альтернативной продажи нет) — проверить у "
    "компании; и сырьё, и H2-кредит считаются по ней же"),
 "nsr_logistics_usd_t": A([60, 150], "$/т", "company_assumption",
    "биг-бэги/мешки: фасовка + Норильск-Дудинка + перевалка + фрахт СМП "
    "(свой арктический флот НН, Дудинка круглогодична на Мурманск) — "
    "проверить у компании"),
 "free_power_mw": A("15-116 МВт по сценариям при электрификации", "МВт",
    "company_assumption",
    "свободная мощность Норильской энергосистемы (ТЭЦ + Усть-Хантайская/"
    "Курейская ГЭС) — проверить у компании; без неё — огневой обогрев"),
 "h2_metallurgy_use": A("восстановитель / защитные атмосферы в металлургии",
    "-", "company_assumption",
    "металлургическое применение H2 — КАЧЕСТВЕННО, ценность может быть выше "
    "топливной; НЕ монетизировано, проверить у компании"),
 "h2_cofiring_share": A("до ~20 об.% в газовые котлы/горелки без глубокой "
    "переделки", "-", "company_assumption",
    "схема подмеса H2 в топливную сеть ТЭЦ — скрин-практика, проверить "
    "применимость к конкретным котлам у компании"),
}


def V(k):
    return ASSUMPTIONS[k]["value"]


WACC = V("wacc")
LIFE = V("asset_life_y")
AF = (1 - (1 + WACC) ** (-LIFE)) / WACC             # аннуитет 15 лет ~6.81
TD_SMALL, TD_SCALE = V("calib_td_small"), V("calib_td_scale")
GAS_CAL = 180.0                                     # калибровочный газ TEA (труба)
LCOH_SCALE_MID = sum(V("lcoh_scale_usd_kg")) / 2.0  # 1.41
NG_SMALL = V("lcoh_small_usd_kg") - CH4_KG_PER_KG_H2 * GAS_CAL / 1000    # 2.34
NG_SCALE = LCOH_SCALE_MID - CH4_KG_PER_KG_H2 * GAS_CAL / 1000            # 0.69
NG_FLOOR = V("lcoh_scale_usd_kg")[0] - CH4_KG_PER_KG_H2 * GAS_CAL / 1000  # 0.57
B_SCALE = math.log(NG_SMALL / NG_SCALE) / math.log(TD_SCALE / TD_SMALL)
CAP_SHARE = V("capital_share_nongas")
ELEC_KWH_KG = V("elec_kwh_kg_h2")[0]                # 13 — якорь ТЗ


def nongas(td):
    """Non-gas остаток LCOH $/кг H2 по масштабу (та же калибровка, что L3)."""
    return max(NG_SMALL * (td / TD_SMALL) ** (-B_SCALE), NG_FLOOR)


def lcoh_lit(td, gas_usd_t, carbon_usd_t=0.0, sold=0.0):
    """LCOH $/кг H2 — только для воспроизведения литературных якорей."""
    return (CH4_KG_PER_KG_H2 * gas_usd_t / 1000.0 + nongas(td)
            - C_KG_PER_KG_H2 * carbon_usd_t / 1000.0 * sold)


# ------------------------------------------------------- сценарий (вся экономика)
def scenario(mln_m3, cb_share, cb_usd_t, raw_usd_t, logistics_usd_t,
             gas_usd_t, arctic_mult, **_):
    """Один сценарий масштаба/качества: тоннажи, EBITDA, капекс (E1), NPV."""
    ch4_kt = mln_m3 * CH4_T_PER_MLN_M3 / 1000.0
    h2_kt = ch4_kt * H2_PER_T_CH4
    c_kt = ch4_kt * C_PER_T_CH4
    td = h2_kt * 1e6 / 365.0 / 1000.0               # т H2/день
    ng = nongas(td)
    capital_kg = ng * CAP_SHARE                     # $/кг H2, капитальная строка
    opex_kg = ng * (1 - CAP_SHARE)                  # $/кг H2, non-gas опекс
    capex_musd = capital_kg * h2_kt * 1e6 / WACC * arctic_mult / 1e6
    capex_int = capital_kg * 1000.0 / WACC          # $/(т H2/год) ДО севера (E1)
    # --- денежные потоки, $M/год ---
    c_price_avg = cb_share * cb_usd_t + (1 - cb_share) * raw_usd_t
    revenue_c = c_kt * c_price_avg / 1000.0
    logistics = c_kt * logistics_usd_t / 1000.0
    opex = h2_kt * opex_kg                          # кт * $/кг = $M
    h2_credit = h2_kt * CH4_EQ_PER_T_H2 * gas_usd_t / 1000.0
    gas_feed = ch4_kt * gas_usd_t / 1000.0
    ebitda = revenue_c - logistics - opex + h2_credit - gas_feed
    npv = -capex_musd + ebitda * AF
    payback = capex_musd / ebitda if ebitda > 0 else None
    # полная себестоимость тонны углерода (H2-кредит вычтен, логистика включена)
    full_c = ((gas_feed - h2_credit + opex + capex_musd * WACC + logistics)
              * 1000.0 / c_kt)
    return {
        "surplus_mln_m3": mln_m3, "ch4_kt_y": round(ch4_kt, 0),
        "mmscf_d_equiv": round(ch4_kt * 1000 / V("mmscf_conversion"), 1),
        "h2_kt_y": round(h2_kt, 1), "carbon_kt_y": round(c_kt, 0),
        "t_h2_per_day": round(td, 0),
        "carbon_price_avg_usd_t": round(c_price_avg, 0),
        "carbon_full_cost_delivered_usd_t": round(full_c, 0),
        "revenue_carbon_musd_y": round(revenue_c, 1),
        "nsr_logistics_musd_y": round(logistics, 1),
        "opex_nongas_musd_y": round(opex, 1),
        "h2_credit_musd_y": round(h2_credit, 2),
        "gas_feed_musd_y": round(gas_feed, 2),
        "ebitda_musd_y": round(ebitda, 1),
        "capex_musd": round(capex_musd, 0),
        "capex_intensity_usd_per_t_h2_y_E1": round(capex_int, -2),
        "capex_intensity_incl_arctic_usd_per_t_h2_y": round(
            capex_int * arctic_mult, -2),
        "payback_y": round(payback, 1) if payback else None,
        "npv12_15y_musd": round(npv, 0),
        "elec_mw_if_electrified": round(h2_kt * 1e6 * ELEC_KWH_KG / 8000.0
                                        / 1000.0, 0),
    }


# ------------------------------------------------- сценарии pessimistic/base/opt
SCEN = {
 "pessimistic": dict(mln_m3=50, cb_share=0.0, cb_usd_t=1200, raw_usd_t=200,
    logistics_usd_t=150, gas_usd_t=50, arctic_mult=2.0,
    note="избыток 50 млн м3 (36 кт CH4); гейт качества ПРОВАЛЕН — только "
         "сырая сажа $200/т; СМП по верху $150/т; газ $50/т; север x2.0"),
 "base": dict(mln_m3=150, cb_share=0.7, cb_usd_t=1200, raw_usd_t=300,
    logistics_usd_t=100, gas_usd_t=35, arctic_mult=1.6,
    note="избыток 150 млн м3 (107 кт CH4); гейт качества ПРОЙДЕН до денег "
         "(условие base!): 70% тоннажа - carbon black по НИЗУ $1200/т, "
         "остальное сырая $300; СМП $100/т; газ $35/т; север x1.6"),
 "optimistic": dict(mln_m3=400, cb_share=0.9, cb_usd_t=1600, raw_usd_t=400,
    logistics_usd_t=60, gas_usd_t=20, arctic_mult=1.3,
    note="избыток 400 млн м3 (286 кт CH4); 90% carbon black $1600/т; СМП "
         "$60/т (свой флот, обратная загрузка); газ $20/т; север x1.3 — "
         "УГОЛ, давление на цену CB посчитано в carbon_market"),
}


def main():
    # ------------------ воспроизведение якорей (обязано сходиться, иначе модель врёт)
    a_scale = round(lcoh_lit(TD_SCALE, GAS_CAL), 2)
    a_small = round(lcoh_lit(TD_SMALL, GAS_CAL), 2)
    a_c250 = round(lcoh_lit(TD_SMALL, GAS_CAL, 250, 1.0), 2)
    a_c700 = round(lcoh_lit(TD_SMALL, GAS_CAL, 700, 1.0), 2)
    cut250 = (a_small - a_c250) / a_small
    t716 = round(CH4_T_PER_MLN_M3, 0)
    kt_scen = {m: round(m * CH4_T_PER_MLN_M3 / 1000.0, 0)
               for m in (50, 150, 400)}
    anchors = {
     "ch4_t_per_mln_m3": dict(model=t716, anchor=716,
        ok=abs(t716 - 716) <= 1, source="физика 16.04/22.414 (якорь ТЗ)"),
     "surplus_kt_ch4": dict(model=kt_scen, anchor={50: 36, 150: 107, 400: 286},
        ok=all(abs(kt_scen[m] - a) <= 1 for m, a in
               {50: 36, 150: 107, 400: 286}.items()),
        source="ТЗ: 50/150/400 млн м3 = 36/107/286 кт CH4/год"),
     "h2_ch4_equiv": dict(model=CH4_EQ_PER_T_H2, anchor=2.4,
        ok=abs(CH4_EQ_PER_T_H2 - 2.4) < 0.01,
        source="LHV 120/50 МДж/кг: 1 т H2 = 2.4 т CH4 по теплу"),
     "lcoh_scale_no_carbon": dict(model=a_scale, anchor=[1.29, 1.53],
        ok=1.29 <= a_scale <= 1.53, source="лит. TEA, мировой масштаб"),
     "lcoh_small_no_carbon": dict(model=a_small, anchor=3.06,
        ok=abs(a_small - 3.06) < 0.01, source="лит. TEA, 1-10 т/день"),
     "carbon_250_cuts_lcoh_pct": dict(model=round(cut250 * 100, 1),
        anchor="~25%", ok=20 <= cut250 * 100 <= 30,
        source="лит.: углерод $250/т режет LCOH на ~25%"),
     "carbon_700_to_doe": dict(model=a_c700, anchor="~$1/кг (DOE)",
        ok=0.85 <= a_c700 <= 1.15, source="лит.: $700/т добивает до DOE-цели"),
     "mining_floor_cross_check": dict(
        model={"raw_200_400": [200, 400], "cb_1200_2000": [1200, 2000]},
        anchor={"parity": 323, "x3": 751},
        ok=(400 < 751 and 1200 > 751),
        source="econ_vs_mining: сырая сажа НЕ проходит 3x-порог (а низ $200 — "
               "даже паритет $323); carbon black проходит — независимое "
               "подтверждение, что кейс живёт только через гейт качества"),
    }

    # ------------------------------------------------------------- сценарии
    scenarios = {name: dict(scenario(**sc), note=sc["note"])
                 for name, sc in SCEN.items()}
    base = SCEN["base"]
    sb = scenarios["base"]

    # E1-проверка: интенсивности в литературном диапазоне (до севера)
    e1_lo, e1_hi = V("capex_intensity_e1")
    for s in scenarios.values():
        s["capex_intensity_in_E1_band"] = (
            e1_lo <= s["capex_intensity_usd_per_t_h2_y_E1"] <= e1_hi)

    # ------------------------------------------------------------- баланс газа
    bal_lo, bal_hi = V("ngp_balance_bln_m3_y")
    gas_balance = {
     "ngp_balance_bln_m3_y": [bal_lo, bal_hi],
     "surplus_share_of_balance_pct": {
        name: [round(sc["mln_m3"] / (bal_hi * 1000) * 100, 1),
               round(sc["mln_m3"] / (bal_lo * 1000) * 100, 1)]
        for name, sc in SCEN.items()},
     "note": "избыток = 1-2% (pess) / 4-5% (base) / 10-13% (opt) добычи "
             "Норильскгазпрома — фактический свободный объём проверить у "
             "компании (kill при <30 млн м3); огневой обогрев съедает ещё "
             "~10-25% сырьевого газа, если нет свободной эл. мощности",
    }

    # ---------------------------------------- рынок углерода (главный риск цены)
    rf_lo, rf_hi = V("rf_carbon_black_kt_y")
    rf_mid = (rf_lo + rf_hi) / 2.0
    carbon_market = {"rf_carbon_black_market_kt_y": [rf_lo, rf_hi]}
    for name, sc in SCEN.items():
        s = scenarios[name]
        cb_kt = s["carbon_kt_y"] * sc["cb_share"]
        share = cb_kt / rf_mid * 100
        row = dict(carbon_kt_y=s["carbon_kt_y"], cb_kt_y=round(cb_kt, 0),
                   raw_kt_y=round(s["carbon_kt_y"] - cb_kt, 0),
                   share_rf_market_pct=round(share, 1),
                   price_pressure_flag=share > 10.0)
        if share > 10.0:
            stress = scenario(**dict(sc, cb_usd_t=800))
            row["stress_cb_800_ebitda_musd"] = stress["ebitda_musd_y"]
            row["stress_cb_800_npv_musd"] = stress["npv12_15y_musd"]
        carbon_market[name] = row
    carbon_market["main_risk"] = (
        "пиролизный углерод != печной техуглерод по морфологии: шинная "
        "спецификация НЕ гарантирована — без гейта это наполнитель $200-400/т; "
        "в optimistic продаём 193 кт CB/год = ~21% рынка РФ — цена $1600 при "
        "таком объёме не удержится, экспорт по СМП в Азию ОБЯЗАТЕЛЕН "
        "(стресс при $800/т посчитан); base 56 кт = ~6% рынка — терпимо")

    # ---------------------------------------- чувствительность (base, 6 ручек)
    def sv(**over):
        s = scenario(**dict(base, **over))
        return {"ebitda_musd_y": s["ebitda_musd_y"],
                "npv_musd": s["npv12_15y_musd"], "payback_y": s["payback_y"]}

    sens = {
     "cb_usd_t (цена carbon black)": {
        "$%d/т" % p: sv(cb_usd_t=p) for p in (800, 1200, 1600, 2000)},
     "cb_share (доля CB в сбыте, остальное raw $300)": {
        "%.0f%%" % (s * 100): sv(cb_share=s) for s in (0.0, 0.3, 0.7, 1.0)},
     "nsr_logistics_usd_t (Дудинка-СМП)": {
        "$%d/т" % p: sv(logistics_usd_t=p) for p in (60, 100, 150, 250)},
     "surplus_mln_m3 (масштаб)": {
        "%d млн м3" % m: sv(mln_m3=m) for m in (50, 100, 150, 250, 400)},
     "gas_internal_usd_t (внутренняя цена газа)": {
        "$%d/т" % g: sv(gas_usd_t=g) for g in (20, 35, 50)},
     "arctic_capex_mult (+-30% к базе 1.6)": {
        "x%.2f" % m: sv(arctic_mult=m) for m in (1.12, 1.6, 2.08)},
     "notes": {
        "cb_usd_t": "главная ручка: $800 -> EBITDA падает вдвое; ниже ~$500 "
                    "средней цены углерода NPV<0",
        "cb_share": "0% = провал гейта качества = смерть (см. kill); ручка "
                    "эквивалентна цене",
        "gas": "почти плоско: нетто-газ = 0.4 x цена x тоннаж (сырьё минус "
               "H2-кредит) — $20 vs $50 меняет EBITDA на ~$1.3M из $53M; "
               "внутренняя цена газа НЕ решает судьбу проекта",
        "logistics": "терпит даже $250/т при CB-качестве; убивает только "
                     "в связке с сырой сажей (см. kill)",
     }}

    # ------------------------------------------------------------ kill-критерии
    k_raw_log = scenario(**dict(base, cb_share=0.0, raw_usd_t=300,
                                logistics_usd_t=175))
    k_raw_only = scenario(**dict(base, cb_share=0.0, raw_usd_t=300))
    k_30 = scenario(**dict(base, mln_m3=30))
    k_30_stress = scenario(**dict(base, mln_m3=30, cb_usd_t=800))
    kill = {
     "carbon_quality_gate_failed": {
        "rule": "качество сажи НЕ подтверждено пилотом: только сырой углерод "
                "$200-400/т И логистика СМП > $150/т",
        "computed_context": {
            "raw300_nsr175": dict(ebitda_musd=k_raw_log["ebitda_musd_y"],
                                  npv_musd=k_raw_log["npv12_15y_musd"]),
            "raw300_nsr100": dict(ebitda_musd=k_raw_only["ebitda_musd_y"],
                                  npv_musd=k_raw_only["npv12_15y_musd"])},
        "why": "сырая сажа не проходит майнинг-порог $751/т (econ_vs_mining) и "
               "не тянет арктический капекс: EBITDA -3.5M/год при СМП $175 — "
               "смерть; честнее ТЗ: сырой углерод убивает и ПРИ дешёвой "
               "логистике (NPV -142M) — критерий «И» консервативен, реальный "
               "kill жёстче: без гейта качества варианта нет вообще"},
     "no_free_power": {
        "rule": "нет свободной эл. мощности под электрификацию И баланс газа "
                "не допускает огневой обогрев (+10-25% газа)",
        "computed_context": {
            "mw_needed_if_electrified": {n: s["elec_mw_if_electrified"]
                                         for n, s in scenarios.items()},
            "fired_extra_gas_mln_m3_base": [
                round(base["mln_m3"] * V("fired_extra_gas")[0], 0),
                round(base["mln_m3"] * V("fired_extra_gas")[1], 0)]},
        "why": "13 кВт*ч/кг H2 = 44 МВт в base (треть электролиза, но всё "
               "равно заметно для изолированной системы); альтернатива — "
               "жечь ещё 15-38 млн м3/год: если и мощности нет, и газа нет — "
               "теплоснабжение реактора не сходится, стоп"},
     "surplus_below_30_mln_m3": {
        "rule": "подтверждённый избыток газа < 30 млн м3/год (~21 кт CH4)",
        "computed_context": {
            "at_30_base_prices": dict(ebitda_musd=k_30["ebitda_musd_y"],
                                      capex_musd=k_30["capex_musd"],
                                      npv_musd=k_30["npv12_15y_musd"]),
            "at_30_cb_800": dict(ebitda_musd=k_30_stress["ebitda_musd_y"],
                                 npv_musd=k_30_stress["npv12_15y_musd"])},
        "why": "ниже 30 млн м3 эффект масштаба умирает (капекс-интенсивность "
               "> $10k/(т H2/г) с севером): при базовых ценах NPV ~+3M — "
               "около нуля, а при CB $800 уже -28M — запаса на ошибку нет, "
               "для Арктики это смерть"},
     "verdict_scenarios": {},
    }
    for name, sc in SCEN.items():
        s = scenarios[name]
        if sc["cb_share"] == 0.0 and (sc["logistics_usd_t"] >= 150
                                      or s["ebitda_musd_y"] <= 0):
            v = ("МЁРТВ (kill 1: гейт качества провален — только сырая "
                 "сажа, EBITDA %.1fM/год)" % s["ebitda_musd_y"])
        elif sc["mln_m3"] < 30:
            v = "МЁРТВ (избыток < 30 млн м3: kill 3)"
        else:
            v = ("жив УСЛОВНО: гейт качества (NaCl) должен быть пройден ДО "
                 "денег; мощность/огневой обогрев — подтвердить (kill 2)")
            if carbon_market[name].get("price_pressure_flag"):
                v += "; флаг давления цены CB (см. carbon_market)"
        kill["verdict_scenarios"][name] = v

    # ------------------------------------------------------------- наш рычаг
    our_lever = {
     "what": "IP по катализатору и чистоте углерода, НЕ реакторостроение: "
             "(1) NEVPT2-переранжир расплавных сплавов — Cu-Bi некоксующийся "
             "кандидат (DFT -395.02 -> NEVPT2 +0.77 ккал/моль, DFT-ранжир "
             "перевёрнут); (2) соляной слой NaCl — чистота углерода "
             "83 -> <0.1 wt% металла",
     "anchor": "pyrolysis_descriptor_results.json (наш расчёт) + лит. по NaCl "
               "(METHANE_ROUTES_RESEARCH.md)",
     "econ_target": "ИЗ ЭТОЙ МОДЕЛИ: гейт качества переводит углерод из "
                    "$200-400 (мёртво, ниже майнинг-порога $751) в $1200-2000 "
                    "(carbon black) — разница ~$900-1600/т и есть ВСЯ маржа "
                    "N3; ни одна другая ручка (газ, логистика, капекс) "
                    "столько не весит",
     "mirage_lesson": "мираж 6.9x соблюдён: H2 НЕ экспортируется, только "
                      "замещение газа на месте (2.4 т CH4-экв/т H2 = "
                      "$48-120/т H2) + опция металлургии (не монетизирована)",
    }

    # ------------------------------------------------ механическая реализация
    implementation = {
     "block_flow": [
      "1. отбор избытка из сети Норильскгазпрома, узел учёта + сероочистка "
      "ZnO (стренд-газ малосернистый — состав проверить у компании)",
      "2. компрессия до ~20 бар (якорь TEA)",
      "3. рекуперативный подогрев сырья хвостовыми потоками",
      "4. пиролизные модули: барботажные колонны расплава Cu-Bi "
      "(Cu0.45Bi0.55), 1000-1100C; сверху — соляной слой NaCl (наше IP)",
      "5. сепарация углерода: всплытие через NaCl-слой, скимминг, отмывка и "
      "рецикл соли (чистота C <0.1 wt% металла — гейт качества)",
      "6. охлаждение + пылеочистка газа (H2 + непрореагировавший CH4)",
      "7. H2-узел: мембраны/облегчённый PSA -> подмес в топливную сеть ТЭЦ "
      "(до ~20 об.% — проверить у компании) / металлургия (опция); хвост "
      "CH4 -> рецикл и огневой обогрев",
      "8. углеродный узел: сушка, грануляция, фасовка в биг-бэги, склад "
      "северного исполнения",
      "9. логистика: автозимник/дорога Норильск-Дудинка -> порт Дудинка "
      "(круглогодичная навигация на Мурманск, свой флот НН) -> СМП: "
      "Мурманск/Азия",
      "10. интеграция с энергосистемой: H2-подмес замещает газ по "
      "теплотворности (2.4 т CH4-экв/т H2)"],
     "equipment_base_150_mln_m3": {
      "реакторные модули расплава Cu-Bi/NaCl (колонны, футеровка, "
      "рекуперация)": [55e6, 85e6],
      "компрессия сырьевого газа до 20 бар": [8e6, 12e6],
      "газоохлаждение/пылеочистка + H2-узел (мембраны, подмес в сеть ТЭЦ)":
        [10e6, 18e6],
      "углеродный узел (скимминг/отмывка NaCl, сушка, грануляция, фасовка, "
      "склад)": [18e6, 30e6],
      "энергоузел/резервный электрообогрев": [8e6, 15e6],
      "логистическая инфраструктура (биг-бэг/контейнерный парк, перевалка "
      "Дудинка)": [5e6, 12e6],
      "АСУ ТП + ПАЗ (H2, расплав)": [4e6, 7e6],
      "монтаж/фундаменты (мерзлота)/инжиниринг, северное исполнение":
        [30e6, 55e6],
      "note": "сумма $138-234M, середина ~$186M ~ капекс модели $%.0fM "
              "(base, E1 x север 1.6) — скрин +-2-3x сходится"
              % sb["capex_musd"]},
     "footprint_m2": [25000, 45000],
     "footprint_note": "площадка у ТЭЦ/газовой сети Норильска + склад "
                       "углерода; плюс перевалочный склад в Дудинке",
     "utilities": {
      "электричество_МВт": {
        "огневой_обогрев_хвостовым_газом": [4, 8],
        "полная_электрификация": [15, 116],
        "note": "13 кВт*ч/кг H2 (низ лит., ~треть электролиза): 15/44/116 МВт "
                "по сценариям; свободная мощность — company_assumption, "
                "kill-критерий; огневой вариант ест +10-25% газа"},
      "вода_оборотная_м3_ч": [50, 150],
      "азот_нм3_ч": [100, 300],
      "NaCl_подпитка_кт_г": [0.3, 1.5],
      "тепло_рекуперации": "сброс тепла 1000C -> теплосеть Норильска — "
                           "опция, НЕ монетизирована",
      "note": "вода/азот — общеплощадочные сервисы Норильского узла"},
     "staffing": "45-70 чел на base (вахта): 4 смены операторов реакторного "
                 "зала и углеродного узла + КИП/механика северного "
                 "исполнения + лаборатория качества сажи + логистика",
     "rollout": [
      {"stage": "FEED + сверка с компанией (баланс газа, избыток, мощность, "
                "логистика СМП) + расчётный гейт G1", "months": 9},
      {"stage": "лаб/стенд: расплав Cu-Bi + NaCl-слой, 72+ ч, чистота "
                "углерода, образцы сажи покупателям (G2-G3)", "months": 12},
      {"stage": "пилот 1-2 т H2/день у ТЭЦ: зимний пробег, подмес H2, "
                "пробная партия сажи Дудинка -> СМП (G4)", "months": 18},
      {"stage": "первая очередь 50 млн м3/год при пройденных G3-G4",
       "months": 30},
      {"stage": "расширение до 150-400 млн м3/год ПО МЕРЕ контрактов сбыта "
                "углерода (не раньше)", "months": 36}],
     "trl_gates": [
      "G1 (расчёт, наше IP, УЖЕ ЕСТЬ): NEVPT2-переранжир сплавов — Cu-Bi "
      "некоксующийся (DFT -395.02 -> +0.77, "
      "pyrolysis_descriptor_results.json); AIMD/стенд-подтверждение",
      "G2 (лаборатория): NaCl-слой даёт чистоту <0.1 wt% металла на нашем "
      "сплаве; конверсия >=85% при 1100C, 72+ ч",
      "G3 (рынок, ДО денег на очередь): образцы прошли спецификацию carbon "
      "black у шинников/трейдеров Азии ИЛИ контракт >=50% объёма по "
      ">=$800/т с учётом доставки",
      "G4 (пилот): 1000+ ч на норильском газе; зимняя логистика биг-бэгов "
      "Дудинка-СМП пройдена; компоненты модели сходятся в +-30%",
      "G5 (энергетика, company): подтверждён избыток >=50 млн м3/год, схема "
      "H2-подмеса в котлы ТЭЦ и свободная мощность (или газ на огневой "
      "обогрев)"],
    }

    OUT = {
     "model": "Норникель N3: пиролиз стренд-газа Норильска — сажа на СМП + "
              "H2 на месте",
     "tier": "screening_pm_2_3x",
     "not_bankable": True,
     "tier_legend": {
        "anchor_repo": "число из наших файлов или физика (помечено в comment)",
        "market": "рыночный диапазон, порядок величины",
        "company_assumption": "допущение о компании — проверить у неё",
        "literature": "литература/чужая TEA, не наш расчёт"},
     "frame": ("скрин +-2-3x; WACC 12%%, жизнь пиролиза 15 лет, аннуитет "
               "%.2f; LCOH-параметрика калибрована к той же лит. TEA, что L3 "
               "(econ_lukoil_h2_refinery); капекс = капитальная строка TEA "
               "через перпетуитет capital*WACC (интенсивность E1) x северный "
               "коэффициент; H2 ТОЛЬКО на месте (мираж 6.9x)" % AF),
     "consistency_anchors": {
        "skid_twin.json": "7300 т CH4/г = 1 MMscf/д (пересчёт масштаба)",
        "econ_shock_candidates.json": "WACC 12%",
        "econ_vs_mining_results.json": "порог $323/$751 паритет/3x — сырая "
                                       "сажа не проходит, carbon black "
                                       "проходит; жизнь пиролиза 15 лет",
        "econ_lukoil_h2_refinery.py": "та же LCOH-калибровка и carbon black "
                                      "$1200-2000 / raw $200-400; рынок "
                                      "техуглерода РФ 700-1100 кт/г",
        "pyrolysis_descriptor_results.json": "Cu-Bi: Ebind_C DFT -395.02 -> "
                                             "NEVPT2 +0.77 (наше IP, G1)",
        "литературная TEA": "LCOH $1.29-1.53 масштаб / $3.06 малый; углерод "
                            "$250/т -> -25% LCOH; <13 кВт*ч/кг H2"},
     "assumptions": ASSUMPTIONS,
     "anchor_reproduction": anchors,
     "gas_balance": gas_balance,
     "scenarios": scenarios,
     "carbon_market": carbon_market,
     "sensitivity": sens,
     "kill_criteria": kill,
     "our_lever": our_lever,
     "implementation": implementation,
     "honesty": "скрин +-2-3x, НЕ банковская модель. Всё по Норникелю "
                "(баланс 3-4 млрд м3, избыток 50-400 млн м3, газ $20-50/т, "
                "СМП $60-150/т, мощность) — company_assumption. Честные "
                "выводы: (1) это САЖЕВЫЙ завод — H2-кредит $0.05-0.12/кг "
                "копеечный, нетто-газ всего 0.4 т-экв/т CH4, цена газа "
                "почти не влияет; (2) ВСЯ экономика висит на гейте качества "
                "NaCl: raw $200-400 мёртв (ниже майнинг-порога $751), "
                "carbon black $1200+ даёт NPV +$204M в base; (3) optimistic "
                "NPV $1.7B — типичный скрин-мираж: 193 кт CB/г = ~21% рынка "
                "РФ, цена $1600 не удержится, экспорт в Азию обязателен "
                "(стресс $800 посчитан); (4) pessimistic МЁРТВ по kill 1 — "
                "фича модели, а не баг; (5) E1: капекс-интенсивность — "
                "производная чужой TEA, не вендорная котировка.",
    }

    # ------------------------------------------------------------------ сводка
    print("=" * 76)
    print("Норникель N3: пиролиз стренд-газа Норильска — сажа на СМП + H2 на "
          "месте")
    print("(скрин +-2-3x, не банк; H2 только на месте — урок миража 6.9x)")
    print("=" * 76)
    ok_all = all(a["ok"] for a in anchors.values())
    print("якоря воспроизведены:", "OK" if ok_all else
          "РАСХОЖДЕНИЕ! " + json.dumps(anchors, ensure_ascii=False))
    print("  716 т CH4/млн м3 | 50/150/400 млн м3 = %.0f/%.0f/%.0f кт | "
          "LCOH $%.2f/$%.2f | C$250 режет %.0f%%"
          % (kt_scen[50], kt_scen[150], kt_scen[400], a_scale, a_small,
             cut250 * 100))
    print("  майнинг-порог: raw $200-400 < $751 (мёртв), CB $1200-2000 > "
          "$751 (проходит 3x)")
    print("\n%-13s%8s%8s%7s%7s%9s%9s%8s%7s%9s" %
          ("сценарий", "млн м3", "кт CH4", "кт C", "кт H2", "EBITDA",
           "капекс", "окуп.", "NPV", "МВт(эл)"))
    for name, s in scenarios.items():
        pb = "%.1fг" % s["payback_y"] if s["payback_y"] else "-"
        print("%-13s%8d%8.0f%7.0f%7.1f%8.1fM%8.0fM%8s%6.0fM%9.0f" %
              (name, s["surplus_mln_m3"], s["ch4_kt_y"], s["carbon_kt_y"],
               s["h2_kt_y"], s["ebitda_musd_y"], s["capex_musd"], pb,
               s["npv12_15y_musd"], s["elec_mw_if_electrified"]))
    print("\nbase: углерод full cost $%d/т (delivered) vs средняя цена "
          "$%d/т; H2-кредит $%.2fM vs сырьевой газ $%.2fM (нетто-газ "
          "0.4 т-экв/т CH4)"
          % (sb["carbon_full_cost_delivered_usd_t"],
             sb["carbon_price_avg_usd_t"], sb["h2_credit_musd_y"],
             sb["gas_feed_musd_y"]))
    print("капекс-интенсивность E1: $%.0f/(т H2/г) до севера, $%.0f с "
          "севером x1.6 (лит. диапазон %d-%d)"
          % (sb["capex_intensity_usd_per_t_h2_y_E1"],
             sb["capex_intensity_incl_arctic_usd_per_t_h2_y"], e1_lo, e1_hi))
    print("\nчувствительность (base): [EBITDA $M/г | NPV $M | окуп. лет]")
    for knob, vals in sens.items():
        if knob == "notes":
            continue
        row = "  ".join("%s: %s|%s|%s" %
                        (k, v["ebitda_musd_y"], v["npv_musd"], v["payback_y"])
                        for k, v in vals.items())
        print("  %-44s %s" % (knob, row))
    print("\nрынок углерода: base %d кт CB/г = %.1f%% рынка РФ; optimistic "
          "%d кт = %.1f%% — давление цены, стресс CB $800: EBITDA $%.1fM"
          % (carbon_market["base"]["cb_kt_y"],
             carbon_market["base"]["share_rf_market_pct"],
             carbon_market["optimistic"]["cb_kt_y"],
             carbon_market["optimistic"]["share_rf_market_pct"],
             carbon_market["optimistic"]["stress_cb_800_ebitda_musd"]))
    print("kill: (1) гейт качества провален (raw И СМП>$150: EBITDA "
          "%.1fM) | (2) нет мощности 15-116 МВт и нельзя жечь газ | "
          "(3) избыток <30 млн м3 (при CB$800 NPV %.0fM)"
          % (k_raw_log["ebitda_musd_y"], k_30_stress["npv12_15y_musd"]))
    print("вердикт:", json.dumps(kill["verdict_scenarios"],
                                 ensure_ascii=False, indent=1))
    print("наш рычаг:", our_lever["econ_target"])

    path = os.path.join(DIR, "econ_nornickel_pyro_nsr_results.json")
    with open(path, "w") as f:
        json.dump(OUT, f, indent=1, ensure_ascii=False)
    print("\nwrote %s" % os.path.basename(path))


if __name__ == "__main__":
    main()
