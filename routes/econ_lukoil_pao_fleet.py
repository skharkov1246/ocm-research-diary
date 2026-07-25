#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
routes/econ_lukoil_pao_fleet.py — Лукойл L1: флот ПАО-скидов на ПНГ с captive-
сбытом через ЛЛК (базовые масла IV группы). Углубление ПАО-варианта skid_twin.py
до флотовой модели: масс-баланс по стадиям (OCM -> ЛАО -> ПАО), capex-WBS,
флот 5/20/50 скидов против рынка ПАО (~700 кт/г) и спроса ЛЛК (30-100 кт/г —
проверить у компании), NPV @ WACC 12% / 12 лет; отдельно «труба» для площадок
без факела.

Скрин-уровень ±2-3x, НЕ банковская модель. Якоря — из наших же файлов:
  skid_twin.json          7300 т CH4/г на скид; ПАО 0.26 т/т CH4, $3200/т,
                          opex $650/т, capex $5.0M; EBITDA факел $4.46M,
                          труба $3.37M; этилен 0.305 т/т (Y=35% мол, риск A1)
  econ_ysz_flux.json      YSZ trim-стек (10% потока O) $29-40k, 48.7 кВт
  econ_vs_mining_results  пол цены продукта: $323/т паритет, $751/т = 3x майнинг
Данные Лукойл/ЛЛК — ТОЛЬКО диапазоны с tier='company_assumption' и пометкой
«проверить у компании». Ничего банковского здесь нет и не задумано.

Запуск: python3 routes/econ_lukoil_pao_fleet.py
     -> routes/econ_lukoil_pao_fleet_results.json
"""
import json
import os

DIR = os.path.dirname(os.path.abspath(__file__))
MW_C2H4_PER_2CH4 = 28.0 / 32.08     # как в skid_twin: т C2H4 на т CH4 при Y=1

# ------------------------------------------------------ допущения (все крутилки)
ASSUMPTIONS = {
    # --- якоря из репозитория ---
    "ch4_t_y_per_skid": dict(value=7300.0, unit="т CH4-экв/год", tier="anchor_repo",
        source="skid_twin.json: скид 1 MMscf/д"),
    "wacc": dict(value=0.12, unit="доля/год", tier="anchor_repo",
        source="skid_twin / econ_shock_candidates: 12%"),
    "life_years": dict(value=12, unit="лет", tier="anchor_repo",
        source="skid_twin: жизнь актива 12 лет"),
    "gas_flare_usd_t": dict(value=30.0, unit="$/т CH4", tier="anchor_repo",
        source="skid_twin GAS.flare; реальный диапазон $0-100/т, бывает отрицательно"),
    "gas_pipeline_usd_t": dict(value=180.0, unit="$/т CH4", tier="anchor_repo",
        source="skid_twin GAS.pipeline (~$3.5/MMBtu)"),
    "ocm_y_mol_base": dict(value=0.35, unit="мол. доля CH4->C2H4 за проход",
        tier="anchor_repo",
        source="skid_twin ethylene_I4; Y=35% НЕ доказан — риск A1 (центр даёт "
               "1.4x, нужен стек рычагов); gate G1 ниже"),
    "oligo_total_yield": dict(value=0.85, unit="т ПАО/т C2H4", tier="anchor_repo",
        source="skid_twin PAO: 0.305 x 0.85 = 0.26 т ПАО/т CH4 сквозной"),
    "pao_price_usd_t": dict(value=3200.0, unit="$/т", tier="anchor_repo",
        source="skid_twin PAO.price"),
    "pao_opex_usd_t": dict(value=650.0, unit="$/т ПАО", tier="anchor_repo",
        source="skid_twin PAO.opex_per_t (олигомеризация+гидрофиниш+катализатор)"),
    "capex_skid_musd": dict(value=5.0, unit="M$", tier="anchor_repo",
        source="skid_twin PAO.capex_musd; WBS-разбивка ниже; риск E1 +-2-3x"),
    "fixed_om_musd_y": dict(value=0.15, unit="M$/год", tier="anchor_repo",
        source="skid_twin: фикс. O&M поверх удельного опекса"),
    "ysz_trim_stack_usd": dict(value=[29219, 39728], unit="$", tier="anchor_repo",
        source="econ_ysz_flux.json trim-10%: 1073K / 973K"),
    "ysz_trim_power_kw": dict(value=48.7, unit="кВт", tier="anchor_repo",
        source="econ_ysz_flux.json trim-10%"),
    "mining_floor_usd_t": dict(value=[323, 751], unit="$/т продукта",
        tier="anchor_repo",
        source="econ_vs_mining_results.json: паритет / 3x к чистой майнинга; "
               "ПАО $3200/т проходит 3x-порог с запасом"),
    # --- рынок ---
    "pao_market_global_kt_y": dict(value=700.0, unit="кт/год", tier="market",
        comment="мировой рынок ПАО, порядок величины; доля >5% = давление на цену"),
    "pao_price_scen_range": dict(value=[2400, 3600], unit="$/т", tier="market",
        comment="диапазон цен ПАО (низковязкие марки) для сценариев pess/opt"),
    # --- допущения о компании (ПРОВЕРИТЬ У КОМПАНИИ) ---
    "llk_group4_demand_kt_y": dict(value=[30, 100], unit="кт/год",
        tier="company_assumption",
        comment="спрос ЛЛК-Интернешнл на базовые масла IV группы — "
                "проверить у компании"),
    "lukoil_flare_sites": dict(value=[20, 100], unit="площадок >=1 MMscf/д",
        tier="company_assumption",
        comment="кол-во ПНГ/факельных площадок Лукойла с достаточным газом — "
                "проверить у компании"),
    "captive_discount": dict(value=[0.0, 0.10], unit="доля от рыночной цены",
        tier="company_assumption",
        comment="дисконт captive-поставки в ЛЛК; в базе НЕ применён (0) — "
                "проверить у компании"),
    "png_c2plus_share": dict(value=[0.2, 0.4], unit="масс. доля C2+ в ПНГ",
        tier="company_assumption",
        comment="ПНГ != чистый CH4: считаем в CH4-экв, C3+ отгоняем в "
                "топливо/ШФЛУ; состав по площадкам — проверить у компании"),
    # --- литература / скрин-разбивки ---
    "oligo_split_c2h4_to_lao": dict(value=0.92, unit="т ЛАО/т C2H4",
        tier="literature",
        comment="скрин-разбивка: 0.92 x 0.92 = 0.846 ~ 0.85; якорем закреплено "
                "ТОЛЬКО произведение 0.85 (skid_twin)"),
    "oligo_split_lao_to_pao": dict(value=0.92, unit="т ПАО/т ЛАО (с гидрофинишем)",
        tier="literature",
        comment="вторая половина разбивки; сам сплит — скрин-уровень"),
    "pao_group3_floor_usd_t": dict(value=2000.0, unit="$/т", tier="literature",
        comment="ниже ~$2000/т ПАО теряет премию к Group III — продукт-тезис "
                "умирает (kill-критерий)"),
}

def V(k):
    return ASSUMPTIONS[k]["value"]

AF = (1 - (1 + V("wacc")) ** (-V("life_years"))) / V("wacc")   # аннуитет 6.19

# --------------------------------------------- сценарии (pess / base / opt)
SCENARIOS = {
    "pessimistic": dict(y_mol=0.25, oligo=0.80, price=2400.0, opex_t=750.0,
        capex=8.0, gas=60.0,
        note="Y=25% (риск A1 частично сбылся), олиго 0.80, цена — низ рынка, "
             "capex к верху E1, газ — верх факельного диапазона"),
    "base": dict(y_mol=0.35, oligo=0.85, price=3200.0, opex_t=650.0,
        capex=5.0, gas=30.0,
        note="якорь skid_twin как есть; Y=35% НЕ доказан (A1)"),
    "optimistic": dict(y_mol=0.35, oligo=0.88, price=3600.0, opex_t=600.0,
        capex=4.0, gas=0.0,
        note="газ бесплатный (за факел платят), олиго 0.88, серийный capex"),
}
FLEET_SIZES = (5, 20, 50)

# ------------------------------------------------- capex WBS (риск E1 +-2-3x)
CAPEX_WBS = {
    "компрессия + очистка ПНГ (сепарация, осушка, S-guard)": [0.6, 1.0],
    "chemical-looping OCM реактор + YSZ trim-стек": [0.8, 1.3],
    "закалка + выделение C2H4 (компрессия, короткая колонна)": [0.5, 0.8],
    "олигомеризация C2H4 -> ЛАО -> ПАО": [0.6, 0.9],
    "гидрофиниш ПАО (реактор + H2-узел)": [0.4, 0.6],
    "хранение/налив (ёмкости, бочки, обвалование)": [0.2, 0.4],
    "электрика + АСУ ТП": [0.3, 0.5],
    "монтаж, рама скида, обвязка": [0.4, 0.7],
}
WBS_LOW = round(sum(v[0] for v in CAPEX_WBS.values()), 2)
WBS_HIGH = round(sum(v[1] for v in CAPEX_WBS.values()), 2)
WBS_MID = round((WBS_LOW + WBS_HIGH) / 2, 2)


def skid(y_mol, oligo, price, opex_t, capex, gas, **_):
    """Один скид: масс-баланс + экономика. Формула 1:1 как в skid_twin."""
    ch4 = V("ch4_t_y_per_skid")
    c2h4 = ch4 * y_mol * MW_C2H4_PER_2CH4
    pao = c2h4 * oligo
    rev = pao * price / 1e6
    gas_cost = ch4 * gas / 1e6
    opex = pao * opex_t / 1e6 + V("fixed_om_musd_y")
    ebitda = rev - gas_cost - opex
    npv = -capex + ebitda * AF
    return dict(
        c2h4_t_y=round(c2h4, 0), pao_t_y=round(pao, 0),
        through_yield_t_per_t_ch4=round(y_mol * MW_C2H4_PER_2CH4 * oligo, 3),
        revenue_musd=round(rev, 2), ebitda_musd=round(ebitda, 2),
        capex_musd=capex,
        payback_y=round(capex / ebitda, 1) if ebitda > 0 else None,
        npv_musd=round(npv, 1))


def fleet(scen, n):
    """Флот из n скидов (все на факеле, одновременный ввод — скрин)."""
    s = skid(**scen)
    out_kt = s["pao_t_y"] * n / 1000.0
    share = out_kt / V("pao_market_global_kt_y") * 100
    llk_lo, llk_hi = V("llk_group4_demand_kt_y")
    if out_kt <= llk_lo:
        captive = "влезает в captive ЛЛК даже при min спросе"
    elif out_kt <= llk_hi:
        captive = (f"нужен спрос ЛЛК >= {out_kt:.0f} кт/г (внутри диапазона "
                   f"{llk_lo}-{llk_hi}) — проверить у компании")
    else:
        captive = "НЕ влезает в captive даже при max спросе — нужны внешние продажи"
    d = dict(n_skids=n, pao_kt_y=round(out_kt, 1),
             ebitda_fleet_musd_y=round(s["ebitda_musd"] * n, 1),
             capex_fleet_musd=round(s["capex_musd"] * n, 0),
             npv_fleet_musd=round(s["npv_musd"] * n, 0),
             share_pao_market_pct=round(share, 1),
             captive_check=captive,
             price_pressure_flag=share > 5.0)
    if d["price_pressure_flag"]:
        stress = skid(**dict(scen, price=scen["price"] * 0.8))
        d["ebitda_fleet_at_price_minus20pct_musd_y"] = round(
            stress["ebitda_musd"] * n, 1)
    return d


def main():
    base = SCENARIOS["base"]
    sb = skid(**base)
    sb_pipe = skid(**dict(base, gas=V("gas_pipeline_usd_t")))

    # -------- сверка с якорями репо (обязана сходиться, иначе модель врёт)
    anchors = {
        "c2h4_t_y": dict(model=sb["c2h4_t_y"], anchor=2230.0,
                         source="skid_twin/econ_ysz_flux"),
        "pao_t_y": dict(model=sb["pao_t_y"], anchor=1896.0, source="skid_twin"),
        "ebitda_flare_musd": dict(model=sb["ebitda_musd"], anchor=4.46,
                                  source="skid_twin PAO flare"),
        "ebitda_pipeline_musd": dict(model=sb_pipe["ebitda_musd"], anchor=3.37,
                                     source="skid_twin PAO pipeline"),
        "through_yield": dict(model=sb["through_yield_t_per_t_ch4"], anchor=0.26,
                              source="skid_twin 0.26 т ПАО/т CH4"),
    }
    for a in anchors.values():
        a["ok"] = abs(a["model"] - a["anchor"]) <= max(0.01 * a["anchor"], 0.01)

    # -------- масс-баланс базового скида (детально, по стадиям)
    ch4 = V("ch4_t_y_per_skid")
    mass_balance = {
        "png_in_ch4_eq_t_y": ch4,
        "ch4_into_c2h4_t_y": round(ch4 * base["y_mol"], 0),
        "c2h4_t_y": sb["c2h4_t_y"],
        "lao_t_y": round(sb["c2h4_t_y"] * V("oligo_split_c2h4_to_lao"), 0),
        "pao_t_y": sb["pao_t_y"],
        "tail_gas_ch4_eq_t_y": round(ch4 * (1 - base["y_mol"]), 0),
        "note": "хвост (непрореагировавший CH4 + CO/CO2 переокисления) — дожиг "
                "в мини-ТЭС площадки, НЕ монетизирован; вода OCM — в дренаж; "
                "ЛАО-строка — скрин-разбивка (якорь только сквозные 0.85)",
    }
    yield_decomposition = {
        "ocm_c2h4_per_ch4": dict(value=round(base["y_mol"] * MW_C2H4_PER_2CH4, 3),
            formula="Y_mol x 28/32.08", risk="A1: Y=35% НЕ доказан, gate G1"),
        "c2h4_to_lao": V("oligo_split_c2h4_to_lao"),
        "lao_to_pao_hydrofinish": V("oligo_split_lao_to_pao"),
        "oligo_total_anchored": V("oligo_total_yield"),
        "through_pao_per_ch4": sb["through_yield_t_per_t_ch4"],
    }

    # -------- сценарии: факел и труба
    scen_out = {}
    for name, p in SCENARIOS.items():
        scen_out[name] = dict(
            params=dict(p), flare=skid(**p),
            pipeline=skid(**dict(p, gas=V("gas_pipeline_usd_t"))))

    # -------- флот 5/20/50 (по трём сценариям)
    fleet_out = {name: {f"n{n}": fleet(p, n) for n in FLEET_SIZES}
                 for name, p in SCENARIOS.items()}

    # -------- чувствительность (база, факел): 6 ручек
    sens = {
        "ocm_y_mol (риск A1)": {
            f"Y={y:.0%}": skid(**dict(base, y_mol=y))["ebitda_musd"]
            for y in (0.20, 0.25, 0.30, 0.35)},
        "pao_price": {
            f"{m:+.0%}": skid(**dict(base, price=base["price"] * (1 + m)))
            ["ebitda_musd"] for m in (-0.3, 0.0, 0.3)},
        "opex_per_t": {
            f"{m:+.0%}": skid(**dict(base, opex_t=base["opex_t"] * (1 + m)))
            ["ebitda_musd"] for m in (-0.3, 0.0, 0.3)},
        "gas_usd_t": {
            f"${g:.0f}/т": skid(**dict(base, gas=g))["ebitda_musd"]
            for g in (0.0, 30.0, 100.0, 180.0)},
        "oligo_yield": {
            f"{o}": skid(**dict(base, oligo=o))["ebitda_musd"]
            for o in (0.80, 0.85, 0.88)},
        "capex_musd (риск E1) -> [окупаемость лет, NPV M$]": {
            f"${c}M": [skid(**dict(base, capex=c))["payback_y"],
                       skid(**dict(base, capex=c))["npv_musd"]]
            for c in (WBS_LOW, 5.0, WBS_HIGH, 10.0, 15.0)},
        "unit": "EBITDA M$/год на скид, факел (если не сказано иное)",
    }

    # -------- kill-критерии (с посчитанным контекстом, без драматизации)
    k_y = skid(**dict(base, y_mol=0.20))
    k_cap = skid(**dict(base, capex=10.0, gas=V("gas_pipeline_usd_t")))
    k_pr = skid(**dict(base, price=2000.0))
    kill = {
        "ocm_yield_per_pass_mol": dict(kill_below=0.20,
            computed_context=f"при Y=20% EBITDA факел ${k_y['ebitda_musd']}M — "
                "ещё жив, но труба падает к ~$1.3M, а рецикл/энергетика "
                "съедают остальное; ниже 20% однопроходного вариант закрыть"),
        "capex_per_skid_musd": dict(kill_above=10.0,
            computed_context=f"при $10M окупаемость на трубе "
                f"{k_cap['payback_y']} г; выше — флот проигрывает "
                "большому заводу и дисциплина E1 требует стоп"),
        "pao_price_usd_t": dict(kill_below=2000.0,
            computed_context=f"при $2000/т EBITDA факел ${k_pr['ebitda_musd']}M "
                "ещё >0, но исчезает премия IV группы к Group III — "
                "умирает сам продукт-тезис (см. pao_group3_floor_usd_t)"),
        "verdict_scenarios": {
            name: "жив" if (p["y_mol"] >= 0.20 and p["capex"] <= 10.0
                            and p["price"] >= 2000.0) else "МЁРТВ"
            for name, p in SCENARIOS.items()},
    }

    # -------- механическая реализация
    implementation = {
        "block_flow": [
            "приём ПНГ: замер, входной сепаратор конденсата",
            "компрессия + очистка: осушка, S-guard (ZnO), отгон C3+ в топливо/ШФЛУ",
            "chemical-looping OCM: O-носитель + YSZ trim-стек (10% потока O, 49 кВт)",
            "закалка + выделение C2H4: дожим, короткая ректификация",
            "олигомеризация C2H4 -> ЛАО (деценовая фракция)",
            "олигомеризация ЛАО -> ПАО + гидрофиниш (H2)",
            "хранение/налив ПАО: бочки/автоцистерна при 25C -> ЛЛК",
            "хвостовой газ -> дожиг в мини-ТЭС: тепло + электро площадки",
        ],
        "equipment": [
            dict(item="компрессор ПНГ + блок очистки", cost_musd=[0.6, 1.0]),
            dict(item="CL-OCM реактор (футеровка, теплообмен)", cost_musd=[0.77, 1.26]),
            dict(item="YSZ trim-стек", cost_usd=V("ysz_trim_stack_usd"),
                 note="якорь econ_ysz_flux, trim-10%; малая доля блока OCM"),
            dict(item="узел выделения C2H4", cost_musd=[0.5, 0.8]),
            dict(item="реакторы олигомеризации (2 ступени)", cost_musd=[0.6, 0.9]),
            dict(item="гидрофиниш + H2-узел", cost_musd=[0.4, 0.6]),
            dict(item="ёмкости/налив", cost_musd=[0.2, 0.4]),
            dict(item="электрика + АСУ ТП", cost_musd=[0.3, 0.5]),
            dict(item="монтаж/рама/обвязка", cost_musd=[0.4, 0.7]),
        ],
        "footprint_m2": [400, 700],
        "utilities": {
            "electricity_kw": dict(value=[450, 700],
                split="компрессия 300-400, YSZ trim 49 (якорь), олигомеризация/"
                      "насосы ~100, АСУ/прочее ~50; частично от мини-ТЭС на хвосте"),
            "cooling_water_m3_h": [5, 10],
            "nitrogen_nm3_h": dict(value=[10, 20], note="продувки, инертизация"),
            "h2_t_y": dict(value=[5, 10],
                note="гидрофиниш ~0.3-0.5% масс от ПАО; кандидат-источник — "
                     "хвосты OCM, проверить расчётом"),
        },
        "staffing": dict(pilot="8-10 чел (вахтовые операторы + механик + инженер)",
            series="3-5 чел/скид при кластерном обслуживании и центральной "
                   "диспетчерской", tier="screening"),
        "rollout": [
            dict(stage="FEED + дескриптор-гейт G1", duration="6-9 мес"),
            dict(stage="пилот: 1 скид на факельной площадке", duration="12-18 мес",
                 note="включая 6 мес непрерывного пробега"),
            dict(stage="серия-1: 5 скидов", duration="12 мес"),
            dict(stage="серия-2: до 20 скидов", duration="18-24 мес",
                 note="gate G4: offtake ЛЛК подтверждён"),
            dict(stage="масштаб: до 50 скидов", duration="24+ мес",
                 note="только при подтверждённом спросе ЛЛК >= 40 кт/г "
                      "или внешнем сбыте"),
        ],
        "trl_gates": [
            "G1 (расчёт, до FEED): стек рычагов селективности даёт Y>=30% "
            "однопроходного C2 — снятие риска A1 расчётом",
            "G2 (лаборатория, до пилота): 100+ ч CL-OCM на реальном ПНГ "
            "(C2+, сера), Y>=25%, стабильность O-носителя",
            "G3 (пилот, до серии): 6 мес пробега; деградация YSZ-стека "
            "<2%/1000 ч; ПАО проходит спецификацию IV группы у ЛЛК "
            "(VI, NOACK, CCS); capex-факт <= $10M",
            "G4 (коммерция, до серии-2): контракт offtake ЛЛК >= 10 кт/г",
        ],
    }

    out = {
        "model": "Лукойл L1: флот ПАО-скидов на ПНГ с captive-сбытом ЛЛК",
        "tier": "screening_pm_2_3x",
        "not_bankable": True,
        "anchors_check": anchors,
        "assumptions": ASSUMPTIONS,
        "mass_balance_per_skid_base": mass_balance,
        "yield_decomposition": yield_decomposition,
        "capex_wbs_musd": dict(blocks=CAPEX_WBS, low=WBS_LOW, mid=WBS_MID,
            high=WBS_HIGH,
            note="mid ~ $5.0M якоря skid_twin; риск E1 +-2-3x => реальный "
                 "коридор $2.5-15M; kill при >$10M"),
        "scenarios": scen_out,
        "fleet": dict(note="все скиды на факеле, одновременный ввод, NPV = "
                           "n x NPV скида (скрин; поэтапный rollout снизит NPV); "
                           f"WACC 12%, 12 лет, аннуитет {AF:.2f}",
                      **fleet_out),
        "sensitivity": sens,
        "kill_criteria": kill,
        "implementation": implementation,
        "honesty": "скрин +-2-3x, НЕ банковская модель. Три главных незакрытых "
                   "вопроса: (A1) Y=35% OCM не доказан — вся экономика висит на "
                   "gate G1/G2; (E1) capex скида известен с точностью до 2-3x; "
                   "(рынок) спрос ЛЛК 30-100 кт/г — наше допущение, проверить у "
                   "компании; 50 скидов = 13.5% мирового рынка ПАО — цена $3200 "
                   "при таком объёме не удержится (стресс -20% посчитан)",
    }

    # ------------------------------------------------------------- сводка
    print("=" * 72)
    print("Лукойл L1: флот ПАО-скидов на ПНГ -> ЛЛК   (скрин +-2-3x, не банк)")
    print("=" * 72)
    print("сверка якорей:", "OK" if all(a["ok"] for a in anchors.values())
          else "РАСХОЖДЕНИЕ! " + str(anchors))
    print(f"скид (база): {ch4:.0f} т CH4 -> {sb['c2h4_t_y']:.0f} т C2H4 -> "
          f"{sb['pao_t_y']:.0f} т ПАО (сквозной {sb['through_yield_t_per_t_ch4']})")
    print(f"\n{'сценарий':<13}{'газ':<8}{'т ПАО':>7}{'EBITDA':>8}{'капекс':>8}"
          f"{'окуп.':>7}{'NPV':>8}")
    for name, s in scen_out.items():
        for g in ("flare", "pipeline"):
            r = s[g]
            pb = f"{r['payback_y']}г" if r["payback_y"] else "-"
            print(f"{name:<13}{g:<8}{r['pao_t_y']:>7.0f}{r['ebitda_musd']:>7.2f}M"
                  f"{r['capex_musd']:>7.1f}M{pb:>7}{r['npv_musd']:>7.1f}M")
    print(f"\ncapex WBS: ${WBS_LOW}-{WBS_HIGH}M (mid ${WBS_MID}M ~ якорь $5.0M);"
          " риск E1 +-2-3x")
    print(f"\nфлот (база, факел):  {'n':>4}{'кт/г':>7}{'EBITDA':>9}{'NPV':>9}"
          f"{'доля рынка':>12}")
    for n in FLEET_SIZES:
        f = fleet_out["base"][f"n{n}"]
        print(f"{'':21}{n:>4}{f['pao_kt_y']:>7.1f}"
              f"{f['ebitda_fleet_musd_y']:>8.1f}M{f['npv_fleet_musd']:>8.0f}M"
              f"{f['share_pao_market_pct']:>11.1f}%")
        print(f"{'':25}-> {f['captive_check']}")
    print("\nчувствительность (EBITDA M$/г, скид, факел):")
    for knob, vals in sens.items():
        if knob == "unit":
            continue
        row = "  ".join(f"{k}:{v}" for k, v in vals.items())
        print(f"  {knob:<45} {row}")
    print("\nkill-критерии: Y<20% однопроходный | capex>$10M/скид | ПАО<$2000/т")
    print("вердикт сценариев:", kill["verdict_scenarios"])

    path = os.path.join(DIR, "econ_lukoil_pao_fleet_results.json")
    with open(path, "w") as f:
        json.dump(out, f, indent=1, ensure_ascii=False)
    print(f"\nwrote {os.path.basename(path)}")


if __name__ == "__main__":
    main()
