#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
routes/econ_lukoil_h2_refinery.py — Лукойл L3: бирюзовый пиролиз-H2 на НПЗ
против собственного SMR. НПЗ — то место, где H2 потребляется на месте
(гидроочистка/гидрокрекинг), т.е. ловушка «водород не уезжает с устья» здесь
не работает.

Модель:
 (1) параметрический LCOH пиролиза, калиброванный к литературной TEA
     (tier=literature, НЕ наш расчёт): $1.29-1.53/кг на масштабе, $3.06/кг на
     1-10 т/день без продажи углерода; углерод $250/т режет LCOH на ~25%,
     $700/т добивает до DOE-цели $1/кг — модель обязана воспроизводить эти
     якоря (anchor_reproduction);
 (2) сравнение с SMR: себестоимость SMR-H2 company_assumption $1.0-1.8/кг,
     CO2-интенсивность 9-10 т/т H2 (literature), цена CO2 $0/$20/$50;
 (3) масс-баланс 1 т CH4 -> 0.25 т H2 + 0.75 т C: на 50 кт H2/год нужно
     200 кт CH4/год и рождается 150 кт углерода/год — сбыт углерода в шинную
     промышленность РФ = ГЛАВНЫЙ риск варианта, посчитан явно;
 (4) сценарии pessimistic/base/optimistic (завод 20/50/100 кт H2/год) +
     чувствительности по 6 ручкам + kill-критерии.

Наше IP (в implementation/trl_gates, не в LCOH-числах): NEVPT2-переранжир
расплавных сплавов (Cu-Bi: DFT -395 -> NEVPT2 +0.77 — некоксующийся кандидат,
routes/pyrolysis_descriptor_results.json) и соляной слой NaCl для чистоты
углерода <0.1 wt% металла (лит., METHANE_ROUTES_RESEARCH.md).

Скрин-уровень ±2-3x, НЕ банковская модель. Все цифры по Лукойлу (спрос H2 по
заводам, себестоимость SMR, цена газа на площадке) — company_assumption,
проверять у компании. Запуск:
  python3 routes/econ_lukoil_h2_refinery.py
    -> routes/econ_lukoil_h2_refinery_results.json
"""
import json
import math
import os

DIR = os.path.dirname(os.path.abspath(__file__))

# --------------------------------------------------- масс-баланс (жёсткий якорь)
H2_PER_T_CH4 = 0.25      # т H2 / т CH4 (лит. TEA: CH4 -> C + 2H2)
C_PER_T_CH4 = 0.75       # т C  / т CH4
CH4_KG_PER_KG_H2 = 4.0   # = 1/0.25
C_KG_PER_KG_H2 = 3.0     # = 0.75/0.25


def A(value, unit, tier, comment):
    """Допущение: значение + единица + уровень достоверности + источник."""
    return {"value": value, "unit": unit, "tier": tier, "comment": comment}


# ------------------------------------------------------- допущения (все явно)
ASSUMPTIONS = {
 # --- якоря репозитория ---
 "wacc": A(0.12, "1/год", "anchor_repo",
           "skid_twin/econ_shock_candidates: 12% годовых на капитал"),
 "asset_life_y": A(12, "лет", "anchor_repo",
                   "skid_twin: жизнь актива 12 лет; NPV через аннуитет"),
 "gas_pipeline_usd_t": A(180.0, "$/т CH4", "anchor_repo",
                         "skid_twin GAS.pipeline (~$3.5/MMBtu) — мировой "
                         "уровень; калибровочная точка LCOH-якорей"),
 "mining_floor_usd_t": A([323, 751], "$/т продукта", "anchor_repo",
                         "econ_vs_mining_results: паритет / 3x к майнингу; "
                         "H2 ~$1400/т формально проходит 3x-порог, но здесь "
                         "конкурент не майнинг, а собственный SMR завода"),
 "cu_bi_nevpt2": A({"Ebind_C_dft": -395.02, "Ebind_C_nevpt2": 0.77},
                   "ккал/моль (относит. вулкан)", "anchor_repo",
                   "pyrolysis_descriptor_results.json: NEVPT2 переворачивает "
                   "DFT-ранжир — Cu-Bi слабейший связыватель C = "
                   "некоксующийся кандидат; наше IP, гейт G1"),
 # --- литература (валидировано НЕ нами — чужая TEA и статьи) ---
 "lcoh_scale_usd_kg": A([1.29, 1.53], "$/кг H2", "literature",
                        "TEA пиролиза в расплаве, 20 бар / 1100C, NiBi+NaBr, "
                        "мировой масштаб, БЕЗ продажи углерода"),
 "lcoh_small_usd_kg": A(3.06, "$/кг H2", "literature",
                        "та же TEA, малый масштаб 1-10 т H2/день, без "
                        "продажи углерода"),
 "carbon_credit_lit": A("углерод $250/т режет LCOH на ~25%; $700/т добивает "
                        "до DOE-цели $1/кг", "-", "literature",
                        "модель ОБЯЗАНА воспроизводить (anchor_reproduction)"),
 "elec_kwh_kg_h2": A([13, 17], "кВт·ч/кг H2", "literature",
                     "пиролиз ест <1/3 электричества электролиза (~52 "
                     "кВт·ч/кг); в TEA-якоре тепло — от сжигания части "
                     "хвостового газа, электрификация — опция"),
 "co2_smr_t_per_t_h2": A(9.5, "т CO2/т H2", "literature",
                         "SMR 9-10 т CO2/т H2 (сырьё+топливо); в модели 9.5"),
 "smr_gas_kg_per_kg_h2": A(3.3, "кг CH4/кг H2", "literature",
                           "SMR суммарно (сырьё+топливо) ~3.0-3.6; "
                           "используется только в kill-контексте "
                           "(соподвижка SMR с ценой газа)"),
 "carbon_price_raw_usd_t": A([200, 400], "$/т", "literature",
                             "сырая пиролизная сажа без сертификации"),
 "carbon_price_cb_usd_t": A([1200, 2000], "$/т", "literature",
                            "качество carbon black (шинная спецификация) — "
                            "ТОЛЬКО после гейта G3"),
 "nacl_layer": A("чистота углерода 83 -> <0.1 wt% металла", "-", "literature",
                 "соляной слой поверх расплава; дескриптор адгезии "
                 "Baumli/Kaptay: NaCl 66 мДж/м2 — оптимум (дёшев, нетоксичен);"
                 " METHANE_ROUTES_RESEARCH.md, гейт G2"),
 "trl": A("3-8", "TRL", "literature",
          "скид-модули коммерческие (Modern Hydrogen, Monolith с $1B "
          "DOE-займом); наша новизна — катализатор (IP), не реактор"),
 # --- калибровка параметрики (скрин-выбор, честно помечен) ---
 "calib_td_small": A(5.0, "т H2/день", "literature",
                     "середина лит. диапазона «малый масштаб 1-10 т/д»"),
 "calib_td_scale": A(150.0, "т H2/день", "literature",
                     "«мировой масштаб» TEA ~ 55 кт/г — скрин-выбор точки; "
                     "между точками — степенная интерполяция non-gas остатка"),
 "capital_share_nongas": A(0.5, "доля", "literature",
                           "доля капитальной строки в non-gas остатке LCOH "
                           "(в TEA типично 40-60%); из неё выводится capex"),
 "carbon_disposal_usd_t": A(20.0, "$/т", "market",
                            "вывоз/размещение НЕпроданного углерода "
                            "(инертный V класс); лит. TEA это не считает — "
                            "мы добавляем"),
 # --- рынок ---
 "co2_price_usd_t_scen": A({"pessimistic": 0, "base": 20, "optimistic": 50},
                           "$/т CO2", "market",
                           "РФ сегодня фактически $0; сценарии — внутренняя "
                           "цена углерода / трансграничные механизмы"),
 "rf_carbon_black_kt_y": A([700, 1100], "кт/год", "market",
                           "производство техуглерода в РФ, порядок величины "
                           "(Омск/Ярославль/Нижнекамск); внутренний спрос "
                           "шинников МЕНЬШЕ (много экспорта) — уточнить"),
 # --- данные компании: ТОЛЬКО диапазоны, проверить у компании ---
 "refinery_h2_demand_kt_y": A({"pessimistic": 20, "base": 50,
                               "optimistic": 100}, "кт H2/год",
                              "company_assumption",
                              "спрос H2 гидрокрекинговых НПЗ (кандидаты: "
                              "Пермь / Волгоград / Н.Новгород-Кстово) — "
                              "фактические объёмы проверить у компании"),
 "smr_cost_usd_kg": A([1.0, 1.8], "$/кг H2", "company_assumption",
                      "себестоимость СОБСТВЕННОГО SMR-H2 на НПЗ; низ = "
                      "амортизированный SMR на дешёвом газе — проверить у "
                      "компании (главный компаратор!)"),
 "gas_at_refinery_usd_t": A([100, 180], "$/т CH4", "company_assumption",
                            "цена газа на площадке НПЗ РФ (внутренняя ниже "
                            "мировой $180) — проверить у компании; база $140"),
 "tire_offtake_share": A([0.3, 0.9], "доля углерода в сбыт",
                         "company_assumption",
                         "какую долю пиролизного углерода реально примут "
                         "шинники/резинотехника РФ — проверить у компании; "
                         "ГЛАВНЫЙ риск варианта"),
}


def V(k):
    return ASSUMPTIONS[k]["value"]


WACC = V("wacc")
AF = (1 - (1 + WACC) ** (-V("asset_life_y"))) / WACC     # аннуитет ~6.19
GAS_REF = V("gas_pipeline_usd_t")                        # калибровочный газ
TD_SMALL = V("calib_td_small")
TD_SCALE = V("calib_td_scale")
LCOH_SCALE_MID = sum(V("lcoh_scale_usd_kg")) / 2         # 1.41
# non-gas остаток = LCOH_лит - газовая строка при калибровочном газе $180/т
NG_SMALL = V("lcoh_small_usd_kg") - CH4_KG_PER_KG_H2 * GAS_REF / 1000   # 2.34
NG_SCALE = LCOH_SCALE_MID - CH4_KG_PER_KG_H2 * GAS_REF / 1000           # 0.69
NG_FLOOR = V("lcoh_scale_usd_kg")[0] - CH4_KG_PER_KG_H2 * GAS_REF / 1000  # .57
B_SCALE = math.log(NG_SMALL / NG_SCALE) / math.log(TD_SCALE / TD_SMALL)
CAP_SHARE = V("capital_share_nongas")
DISPOSAL = V("carbon_disposal_usd_t")
CO2_SMR = V("co2_smr_t_per_t_h2")


# ------------------------------------------------------------ параметрика LCOH
def lcoh(td, gas, carbon_usd_t, carbon_sold, capex_mult=1.0,
         with_disposal=True):
    """LCOH $/кг H2 по компонентам. td — т H2/день (масштаб)."""
    nongas = max(NG_SMALL * (td / TD_SMALL) ** (-B_SCALE), NG_FLOOR)
    capital = nongas * CAP_SHARE * capex_mult
    opex = nongas * (1 - CAP_SHARE)
    gas_usd = CH4_KG_PER_KG_H2 * gas / 1000.0
    credit = C_KG_PER_KG_H2 * carbon_usd_t / 1000.0 * carbon_sold
    disposal = (C_KG_PER_KG_H2 * (1 - carbon_sold) * DISPOSAL / 1000.0
                if with_disposal else 0.0)
    full = gas_usd + capital + opex - credit + disposal
    return {"gas": round(gas_usd, 3), "capital": round(capital, 3),
            "opex_other": round(opex, 3), "carbon_credit": round(-credit, 3),
            "carbon_disposal": round(disposal, 3),
            "full_usd_kg": round(full, 2),
            "cash_usd_kg": round(full - capital, 2)}


# ------------------------------------------------- один H2-трейн на НПЗ vs SMR
def train(h2_kt, gas, carbon_usd_t, carbon_sold, smr_usd_kg, co2_usd_t,
          capex_mult=1.0, **_):
    td = h2_kt * 1000.0 / 365.0
    lc = lcoh(td, gas, carbon_usd_t, carbon_sold, capex_mult)
    kg_y = h2_kt * 1e6
    # capex из капитальной строки (перпетуитет capital*WACC, как в linейке L1/L2)
    capex_musd = lc["capital"] * kg_y / WACC / 1e6
    smr_eff = smr_usd_kg + CO2_SMR * co2_usd_t / 1000.0
    savings_kg = smr_eff - lc["cash_usd_kg"]        # cash: капекс платим отдельно
    savings_musd = savings_kg * kg_y / 1e6
    npv = -capex_musd + savings_musd * AF
    c_kt = h2_kt * C_KG_PER_KG_H2
    return {"h2_kt_y": h2_kt, "t_h2_per_day": round(td, 0),
            "ch4_kt_y": round(h2_kt * CH4_KG_PER_KG_H2, 0),
            "carbon_kt_y": round(c_kt, 0),
            "carbon_sold_kt_y": round(c_kt * carbon_sold, 0),
            "lcoh": lc,
            "smr_eff_usd_kg": round(smr_eff, 2),
            "savings_usd_kg": round(savings_kg, 2),
            "savings_musd_y": round(savings_musd, 1),
            "capex_musd": round(capex_musd, 0),
            "payback_y": (round(capex_musd / savings_musd, 1)
                          if savings_musd > 0 else None),
            "npv_musd": round(npv, 0)}


# --------------------------------------------- сценарии pessimistic/base/opt
SCEN = {
 "pessimistic": dict(h2_kt=20, gas=180, carbon_usd_t=200, carbon_sold=0.3,
                     smr_usd_kg=1.0, co2_usd_t=0, capex_mult=1.4,
                     note="малый НПЗ 20 кт/г; газ мировой $180; углерод — "
                          "сырая сажа $200/т, продаётся 30% (70% в отвал); "
                          "конкурент — амортизированный SMR $1.0/кг; CO2 $0 "
                          "(РФ сегодня); капекс +40%"),
 "base": dict(h2_kt=50, gas=140, carbon_usd_t=250, carbon_sold=0.7,
              smr_usd_kg=1.4, co2_usd_t=20, capex_mult=1.0,
              note="гидрокрекинговый НПЗ 50 кт/г; газ РФ $140; углерод $250/т "
                   "(лит. якорь), сбыт 70%; SMR $1.4/кг; CO2 $20/т"),
 "optimistic": dict(h2_kt=100, gas=100, carbon_usd_t=400, carbon_sold=0.9,
                    smr_usd_kg=1.8, co2_usd_t=50, capex_mult=0.8,
                    note="крупный НПЗ 100 кт/г; газ $100; углерод $400/т "
                         "(верх сырой; carbon black $1200+ — НЕ заложен), "
                         "сбыт 90%; SMR дорогой $1.8; CO2 $50/т; капекс -20%"),
}


def main():
    # ------------------ воспроизведение литературных якорей (обязано сходиться)
    a_scale = lcoh(TD_SCALE, GAS_REF, 0, 0, with_disposal=False)["full_usd_kg"]
    a_small = lcoh(TD_SMALL, GAS_REF, 0, 0, with_disposal=False)["full_usd_kg"]
    a_c250 = lcoh(TD_SMALL, GAS_REF, 250, 1.0,
                  with_disposal=False)["full_usd_kg"]
    a_c700 = lcoh(TD_SMALL, GAS_REF, 700, 1.0,
                  with_disposal=False)["full_usd_kg"]
    cut250 = (a_small - a_c250) / a_small
    anchors = {
     "lcoh_scale_no_carbon": dict(model=a_scale, anchor=[1.29, 1.53],
        ok=1.29 <= a_scale <= 1.53, source="лит. TEA, мировой масштаб"),
     "lcoh_small_no_carbon": dict(model=a_small, anchor=3.06,
        ok=abs(a_small - 3.06) < 0.01, source="лит. TEA, 1-10 т/день"),
     "carbon_250_cuts_lcoh_pct": dict(model=round(cut250 * 100, 1),
        anchor="~25%", ok=20 <= cut250 * 100 <= 30,
        source="лит.: углерод $250/т режет LCOH на 25%"),
     "carbon_700_to_doe": dict(model=a_c700, anchor="~$1/кг (DOE)",
        ok=0.85 <= a_c700 <= 1.15,
        source="лит.: $700/т добивает до DOE-цели"),
     "mass_balance_50kt": dict(
        model=[50 * CH4_KG_PER_KG_H2, 50 * C_KG_PER_KG_H2],
        anchor=[200, 150], ok=(50 * CH4_KG_PER_KG_H2 == 200
                               and 50 * C_KG_PER_KG_H2 == 150),
        source="1 т CH4 -> 0.25 т H2 + 0.75 т C (лит.)"),
    }

    # ---------------------------------------------------------- масс-баланс
    mass_balance = {
     "reaction": "CH4 -> C(тв) + 2 H2 (расплав Cu-Bi ~1000-1100C, 20 бар; "
                 "без CO2, без синтез-газа)",
     "per_kg_h2": {"ch4_kg": CH4_KG_PER_KG_H2, "c_kg": C_KG_PER_KG_H2},
     "per_refinery_base_50kt": {
        "ch4_kt_y": 200, "h2_kt_y": 50, "carbon_kt_y": 150,
        "ch4_mmscf_d_equiv": round(200e3 / 7300.0, 1),
        "co2_avoided_kt_y": round(CO2_SMR * 50, 0),
        "note": "CO2-экономия против SMR 9.5 т/т H2; минус собственные "
                "~0.5-1.5 т/т при огневом обогреве хвостовым газом — скрин; "
                "27 MMscf/д CH4-экв — это ТРУБНЫЙ газ НПЗ, не факел"},
    }

    # ------------------------------------------------------------- сценарии
    scenarios = {}
    for name, sc in SCEN.items():
        scenarios[name] = dict(train(**sc), note=sc["note"])

    base = SCEN["base"]
    tb = scenarios["base"]

    # ------------------------------------- рынок сбыта углерода (главный риск)
    rf_lo, rf_hi = V("rf_carbon_black_kt_y")
    rf_mid = (rf_lo + rf_hi) / 2
    carbon_market = {"rf_carbon_black_market_kt_y": [rf_lo, rf_hi]}
    for name in SCEN:
        s = scenarios[name]
        share = s["carbon_sold_kt_y"] / rf_mid * 100
        carbon_market[name] = dict(
            carbon_kt_y=s["carbon_kt_y"], sold_kt_y=s["carbon_sold_kt_y"],
            landfill_kt_y=round(s["carbon_kt_y"] - s["carbon_sold_kt_y"], 0),
            share_rf_market_pct=round(share, 1),
            price_pressure_flag=share > 10.0)
        if share > 10.0:
            stress = train(**dict(SCEN[name], carbon_usd_t=200))
            carbon_market[name]["savings_at_carbon_200_musd_y"] = \
                stress["savings_musd_y"]
    carbon_market["main_risk"] = (
        "пиролизный углерод != печной техуглерод по морфологии/структуре: "
        "шинная спецификация НЕ гарантирована — без гейта G3 это наполнитель "
        "$200-400/т или отвал; в base продаём 105 кт/г = ~12% рынка РФ, в "
        "optimistic 270 кт/г = ~30% — цена $400 при таком объёме не "
        "удержится (стресс при $200/т посчитан); carbon black $1200-2000 в "
        "экономику СОЗНАТЕЛЬНО не заложен")

    # ------------------------------------------- чувствительности (база, 6 ручек)
    def sv(**kw):
        t = train(**dict(base, **kw))
        return {"lcoh_full": t["lcoh"]["full_usd_kg"],
                "savings_musd_y": t["savings_musd_y"],
                "payback_y": t["payback_y"]}
    sens = {
     "carbon_usd_t": {f"${p}/т": sv(carbon_usd_t=p)
                      for p in (0, 250, 700, 1200)},
     "gas_usd_t": {f"${g}/т": sv(gas=g) for g in (100, 140, 180, 250)},
     "smr_cost_usd_kg": {f"${c}/кг": sv(smr_usd_kg=c)
                         for c in (1.0, 1.4, 1.8)},
     "co2_usd_t": {f"${c}/т": sv(co2_usd_t=c) for c in (0, 20, 50)},
     "carbon_sold_share": {f"{s:.0%}": sv(carbon_sold=s)
                           for s in (0.3, 0.5, 0.7, 0.9)},
     "capex_mult": {f"x{m}": sv(capex_mult=m) for m in (0.8, 1.0, 1.4)},
     "notes": {
      "carbon": "$700+/т требует качества carbon black (гейт G3); при $1200 "
                "LCOH отрицательный — углерод субсидирует H2, главный "
                "продукт де-факто сажа",
      "gas": "SMR в таблице держится фикс ($1.4) — честнее: SMR тоже ест газ "
             "(~3.3 кг/кг), реальный разрыв сжимается, см. kill_criteria",
      "smr": "низ $1.0 = амортизированный SMR — самый опасный компаратор "
             "brownfield; проверить у компании",
      "co2": "$0 — РФ сегодня; выше $20/т CO2-строка SMR добавляет "
             "+$0.19-0.48/кг к его себестоимости"}}

    # -------------------------------------------------------- kill-критерии
    k1a = train(**dict(base, carbon_sold=0.3, co2_usd_t=0, smr_usd_kg=1.0))
    k1b = train(**dict(base, carbon_sold=0.5, co2_usd_t=0, smr_usd_kg=1.0))
    k1c = train(**dict(base, carbon_sold=0.5, co2_usd_t=0, smr_usd_kg=1.4))
    # газ >$250: соподвижка SMR по газу (3.3 кг/кг от базовых $140)
    def kill_gas(g, sold, smr0, co2):
        smr_adj = smr0 + V("smr_gas_kg_per_kg_h2") * (g - base["gas"]) / 1000.0
        return train(**dict(base, gas=g, carbon_sold=sold,
                            smr_usd_kg=smr_adj, co2_usd_t=co2))
    k2a = kill_gas(250, 0.7, 1.4, 20)
    k2b = kill_gas(250, 0.3, 1.0, 0)
    kill = {
     "carbon_unsold_and_co2_zero": {
        "rule": "сбыт углерода <50% (>50% в отвал) И цена CO2 = $0/т",
        "computed_context": {
         "sold30_co2_0_smr1.0": dict(savings_musd_y=k1a["savings_musd_y"],
                                     npv_musd=k1a["npv_musd"]),
         "sold50_co2_0_smr1.0": dict(savings_musd_y=k1b["savings_musd_y"],
                                     npv_musd=k1b["npv_musd"]),
         "sold50_co2_0_smr1.4": dict(savings_musd_y=k1c["savings_musd_y"],
                                     npv_musd=k1c["npv_musd"])},
        "why": "без сбыта углерода и без цены CO2 у пиролиза не остаётся ни "
               "одного из двух козырей: против амортизированного SMR $1.0/кг "
               "NPV уходит в минус (посчитано); жив только против дорогого "
               "SMR $1.4+ — ставка «конкурент неэффективен» не наша, "
               "вариант закрыть"},
     "gas_above_250_usd_t": {
        "rule": "цена газа на площадке > $250/т CH4",
        "computed_context": {
         "gas250_base_carbon_smr_co_moves": dict(
            savings_musd_y=k2a["savings_musd_y"], npv_musd=k2a["npv_musd"]),
         "gas250_plus_carbon_fail": dict(
            savings_musd_y=k2b["savings_musd_y"], npv_musd=k2b["npv_musd"])},
        "why": "пиролиз ест 4.0 кг CH4/кг H2 против ~3.3 у SMR — при дорогом "
               "газе разрыв удельного расхода работает ПРОТИВ нас "
               "(+$0.175/кг на каждые $250/т); при газе >$250 в связке с "
               "плохим сбытом углерода NPV<0 (посчитано), и сама премисса "
               "«дешёвый газ РФ» мертва — стоп"},
     "verdict_scenarios": {
        name: ("МЁРТВ (сбыт углерода <50% при CO2 $0)"
               if (sc["carbon_sold"] < 0.5 and sc["co2_usd_t"] == 0)
               else ("МЁРТВ (газ >$250/т)" if sc["gas"] > 250 else "жив"))
        for name, sc in SCEN.items()},
    }

    # ------------------------------------------------ механическая реализация
    implementation = {
     "block_flow": [
      "1. отбор трубного/заводского газа + сероочистка (аминовая/ZnO — "
      "сервисы НПЗ уже есть)",
      "2. компрессия до ~20 бар (якорь TEA)",
      "3. рекуперативный подогрев сырья хвостовыми потоками",
      "4. пиролизный реактор: барботажные колонны расплава Cu-Bi "
      "(Cu0.45Bi0.55), 1000-1100C; сверху — соляной слой NaCl",
      "5. сепарация углерода: всплытие через NaCl-слой, скимминг, отмывка и "
      "рецикл соли (чистота C <0.1 wt% металла — гейт G2)",
      "6. охлаждение + пылеочистка газа (H2 + непрореагировавший CH4)",
      "7. PSA: H2 99.9%+ -> заводской H2-коллектор; хвост CH4 -> рецикл / "
      "топливо обогрева",
      "8. углеродный узел: сушка, грануляция/пеллеты, склад, отгрузка "
      "(шинники / резинотехника / металлургия)",
      "9. интеграция: SMR остаётся hot-standby, пиролиз замещает базовую "
      "нагрузку H2-коллектора"],
     "equipment_train_50kt": {
      "реакторные модули расплава Cu-Bi/NaCl (параллельные колонны, "
      "футеровка, рекуперация)": [45e6, 70e6],
      "компрессия сырьевого газа до 20 бар": [10e6, 15e6],
      "PSA-узел H2": [15e6, 25e6],
      "углеродный узел (скимминг/отмывка NaCl, сушка, грануляция, склад)":
        [15e6, 25e6],
      "энергоузел/электронагрев": [12e6, 20e6],
      "интеграция с H2-коллектором НПЗ + эстакады": [8e6, 15e6],
      "АСУ ТП + ПАЗ (H2, расплав)": [5e6, 8e6],
      "монтаж/фундаменты/инжиниринг": [15e6, 25e6],
      "note": "сумма $125-203M, середина ~$164M ~ капекс модели $149M "
              "(base) — скрин ±2-3x сходится"},
     "footprint_m2": [15000, 30000],
     "footprint_note": "трейн 50 кт/г на действующем НПЗ, рядом с "
                       "H2-коллектором; плюс склад углерода",
     "utilities": {
      "электричество_МВт": {"огневой_обогрев_хвостовым_газом": [10, 15],
                            "полная_электрификация": [75, 100],
                            "note": "13-17 кВт·ч/кг H2 при электрификации = "
                                    "<1/3 электролиза (~52); в LCOH-якоре "
                                    "сидит огневой вариант"},
      "вода_оборотная_м3_ч": [100, 200],
      "азот_нм3_ч": [100, 300],
      "NaCl_подпитка_кт_г": [1, 3],
      "note": "пар/сероочистка/факел — общезаводские сервисы НПЗ"},
     "staffing": "30-45 чел на трейн 50 кт/г: 4 вахты x (операторы 4-5 + "
                 "углеродный узел 2) + КИП/механика + лаборатория качества "
                 "сажи; ремонты — сервисы НПЗ",
     "rollout": [
      {"stage": "FEED + расчётные гейты G1-G2 (NEVPT2/AIMD сплава, выбор "
                "соли) + заводская сверка SMR-себестоимости", "months": 9},
      {"stage": "лаб/стенд: расплав Cu-Bi + NaCl-слой, 72+ ч, чистота "
                "углерода, образцы сажи шинникам", "months": 12},
      {"stage": "пилот на площадке НПЗ 1-2 т H2/день (LCOH-пилот ~$3/кг — "
                "нормально, меряем компоненты модели)", "months": 18},
      {"stage": "первый трейн 20-50 кт/г (EPC) при выполненных G3-G4",
       "months": 30},
      {"stage": "тираж на 2-3 НПЗ (Пермь/Волгоград/Кстово — проверить "
                "фактический спрос H2)", "months": 24}],
     "trl_gates": [
      "G1 (расчёт, наше IP): NEVPT2-переранжир подтверждён AIMD/стендом — "
      "Cu-Bi не коксуется (DFT -395 -> NEVPT2 +0.77 ккал/моль: DFT-ранжир "
      "перевёрнут, pyrolysis_descriptor_results.json)",
      "G2 (лаборатория, наше IP + лит.): NaCl-слой даёт чистоту углерода "
      "<0.1 wt% металла на НАШЕМ сплаве; конверсия >=85% при 1100C, 72+ ч",
      "G3 (рынок, ДО денег на трейн): образцы углерода прошли спецификацию "
      "шинного завода ИЛИ контракт на сбыт >=50% углерода по >=$200/т",
      "G4 (пилот): 1000+ ч на заводском газе; компоненты LCOH сходятся с "
      "моделью в ±30%; деградация расплава/соли измерена",
      "G5 (интеграция): заводская сверка SMR-себестоимости ($1.0-1.8 — "
      "company_assumption) и схемы hot-standby с H2-коллектором"],
    }

    OUT = {
     "model": "Лукойл L3: бирюзовый пиролиз-H2 на НПЗ против SMR",
     "tier": "screening_pm_2_3x",
     "not_bankable": True,
     "frame": ("скрин ±2-3x; LCOH калиброван к литературной TEA (не наш "
               f"расчёт); WACC 12%, жизнь 12 лет, аннуитет {AF:.2f}; "
               "капитальная строка LCOH -> капекс через перпетуитет "
               "capital*WACC (как в L1/L2); экономия = (SMR_eff - "
               "LCOH_cash), капекс отдельно"),
     "consistency_anchors": {
      "skid_twin.json": "трубный газ $180/т; WACC 12%; жизнь 12 лет; "
                        "7300 т CH4/г = 1 MMscf/д (для пересчёта)",
      "econ_vs_mining_results.json": "пиролиз H2+C уже был лучшим по EBITDA "
                                     "на факеле; здесь другой кейс — "
                                     "captive-замещение SMR на НПЗ",
      "pyrolysis_descriptor_results.json": "Cu-Bi: Ebind_C DFT -395.02 -> "
                                           "NEVPT2 +0.77 (наше IP, гейт G1)",
      "литературная TEA": "LCOH $1.29-1.53/кг масштаб / $3.06 малый; углерод "
                          "$250/т -> -25% LCOH; $700/т -> ~$1/кг (DOE); "
                          "воспроизведено в anchor_reproduction"},
     "assumptions": ASSUMPTIONS,
     "anchor_reproduction": anchors,
     "mass_balance": mass_balance,
     "scenarios": scenarios,
     "carbon_market": carbon_market,
     "sensitivity": sens,
     "kill_criteria": kill,
     "implementation": implementation,
     "honesty": "скрин ±2-3x, НЕ банковская модель. LCOH — чужая TEA, наша "
                "только параметризация (калибровка non-gas остатка и его "
                "деление капитал/опекс 50/50 — скрин-выбор). Всё по Лукойлу "
                "(спрос H2 заводов 20/50/100 кт/г, себестоимость SMR "
                "$1.0-1.8, газ на площадке $100-180) — company_assumption, "
                "проверить у компании. Главный риск — сбыт 60-300 кт/г "
                "углерода: шинная спецификация не гарантирована, в base мы "
                "~12% рынка РФ, в optimistic ~30% (цена поедет). Наш вклад — "
                "IP по катализатору (NEVPT2 Cu-Bi) и чистоте углерода "
                "(NaCl-слой), НЕ реакторостроение: скид-модули уже "
                "коммерческие. Pessimistic-сценарий МЁРТВ по собственному "
                "kill-критерию — это фича модели, а не баг.",
    }

    # ------------------------------------------------------------- сводка
    print("=" * 74)
    print("Лукойл L3: бирюзовый пиролиз-H2 на НПЗ против SMR (скрин ±2-3x)")
    print("=" * 74)
    ok_all = all(a["ok"] for a in anchors.values())
    print("литературные якоря воспроизведены:", "OK" if ok_all else
          "РАСХОЖДЕНИЕ! " + json.dumps(anchors, ensure_ascii=False))
    print(f"  масштаб ${a_scale}/кг (лит. 1.29-1.53) | малый ${a_small}/кг "
          f"(лит. 3.06) | углерод $250 режет {cut250 * 100:.0f}% (лит. ~25%) "
          f"| $700 -> ${a_c700}/кг (DOE ~$1)")
    print(f"масс-баланс: 50 кт H2/г <- 200 кт CH4/г -> 150 кт C/г "
          f"(27.4 MMscf/д); CO2-экономия ~{CO2_SMR * 50:.0f} кт/г vs SMR")
    print(f"\n{'сценарий':<13}{'кт H2':>6}{'LCOH':>7}{'SMR_eff':>8}"
          f"{'эконом.':>9}{'капекс':>8}{'окуп.':>7}{'NPV':>7}  углерод")
    for name, s in scenarios.items():
        pb = f"{s['payback_y']}г" if s["payback_y"] else "-"
        print(f"{name:<13}{s['h2_kt_y']:>6}{s['lcoh']['full_usd_kg']:>7}"
              f"{s['smr_eff_usd_kg']:>8}{s['savings_musd_y']:>8.1f}M"
              f"{s['capex_musd']:>7.0f}M{pb:>7}{s['npv_musd']:>6.0f}M"
              f"  {s['carbon_sold_kt_y']:.0f}/{s['carbon_kt_y']:.0f} кт "
              f"в сбыт")
    print("\nчувствительность (base, 50 кт/г): [LCOH $/кг | экономия M$/г | "
          "окуп. лет]")
    for knob, vals in sens.items():
        if knob == "notes":
            continue
        row = "  ".join(f"{k}: {v['lcoh_full']}|{v['savings_musd_y']}|"
                        f"{v['payback_y']}" for k, v in vals.items())
        print(f"  {knob:<18} {row}")
    print("\nсбыт углерода (главный риск): base 105 кт/г = "
          f"{carbon_market['base']['share_rf_market_pct']}% рынка РФ; "
          f"optimistic {carbon_market['optimistic']['sold_kt_y']:.0f} кт/г = "
          f"{carbon_market['optimistic']['share_rf_market_pct']}% — давление "
          "на цену")
    print("kill: сбыт углерода <50% И CO2 $0 (vs SMR $1.0 NPV "
          f"{k1a['npv_musd']:.0f}M) | газ >$250/т (в связке с плохим сбытом "
          f"NPV {k2b['npv_musd']:.0f}M)")
    print("вердикт сценариев:", kill["verdict_scenarios"])

    path = os.path.join(DIR, "econ_lukoil_h2_refinery_results.json")
    with open(path, "w") as f:
        json.dump(OUT, f, indent=1, ensure_ascii=False)
    print(f"\nwrote {os.path.basename(path)}")


if __name__ == "__main__":
    main()
