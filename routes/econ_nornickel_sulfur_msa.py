#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
routes/econ_nornickel_sulfur_msa.py — Норникель N2: Серная программа -> MSA
(SO3 по отрицательной цене). Идея: Серная программа улавливает SO2 и гонит
его через SO3 в H2SO4, которую нейтрализуют известняком в гипс — чистый
РАСХОД $30-80/т серы. Отбираем боковой поток SO3 (доли процента мощности)
и делаем из него MSA (CH4 + SO3 -> CH3SO3H) на газе Норильскгазпрома —
сера превращается из пассива в маржу. Вывоз — танк-контейнер Дудинка -> СМП.

Себестоимость тонны MSA — ровно модель msa_cost_model.py (якоря: 0.167 т
CH4 + 0.333 т S на тонну, инициатор $100, утилиты $130, fixed $400,
капитальная $494 при капексе $4.5M на 1.1 кт/г, WACC 12%); заменены только
цена серы ($150 комодити -> сценарии -50/0/+50 $/т S: сера как избавление
от пассива) и цена газа ($30 факел -> $20-50 стренд НГП).

Скрин-уровень ±2-3x, НЕ банковская модель. Данные Норникеля (масштаб SO2,
опекс гипса, трансфертная цена серы, газ, СМП-надбавка) — ТОЛЬКО диапазоны
с tier='company_assumption' и пометкой «проверить у компании».

Запуск: python3 routes/econ_nornickel_sulfur_msa.py
     -> routes/econ_nornickel_sulfur_msa_results.json
"""
import json
import os

DIR = os.path.dirname(os.path.abspath(__file__))

# ------------------------------------------------------ допущения (все крутилки)
ASSUMPTIONS = {
    # --- якоря репозитория (msa_cost_model / skid_twin / econ_vs_mining) ---
    "ch4_t_per_t_msa": dict(value=0.167, unit="т CH4/т MSA",
        tier="anchor_repo", source="msa_cost_model: стехиометрия "
        "CH4+SO3->CH3SO3H (96.1)"),
    "s_t_per_t_msa": dict(value=0.333, unit="т S/т MSA (=0.833 т SO3)",
        tier="anchor_repo", source="msa_cost_model: 0.833 т SO3 = 0.333 т S "
        "на тонну MSA"),
    "initiator_usd_t": dict(value=100.0, unit="$/т MSA", tier="anchor_repo",
        source="msa_cost_model: пероксид/электро; болевая точка -> наш "
        "квантовый рычаг (gate G1: цель <$50)"),
    "utilities_usd_t": dict(value=130.0, unit="$/т MSA", tier="anchor_repo",
        source="msa_cost_model: SO3-узел+дистилляция+рецикл",
        comment="в Норильске SO3 приходит ГОТОВЫМ из контактного отделения "
                "Серной программы (свой SO3-узел не нужен — минус к цифре), "
                "но нет тепла сжигания серы и есть арктическое исполнение "
                "(плюс) — якорь оставлен как есть, скрин"),
    "fixed_skid_usd_t": dict(value=400.0, unit="$/т MSA", tier="anchor_repo",
        source="msa_cost_model: труд/обслуживание/накладные малого масштаба "
        "(1.1 кт/г); при 20% конверсии размазывается по тоннажу"),
    "capex_skid_musd": dict(value=4.5, unit="$M/скид", tier="anchor_repo",
        source="msa_cost_model: капитальная $494/т = $4.5M * 12% / 1093 т"),
    "skid_t_y": dict(value=1093.0, unit="т MSA/год/скид", tier="anchor_repo",
        source="msa_cost_model: конверсия CH4->MSA ~2.5% per-pass suite"),
    "wacc": dict(value=0.12, unit="доля/год", tier="anchor_repo",
        source="msa_cost_model / econ_shock_candidates: 12% (аморт.+WACC)"),
    "skid20_t_y": dict(value=8747.0, unit="т MSA/год/скид",
        tier="anchor_repo", source="msa_cost_model: апсайд конверсии 20% "
        "(x8 тоннаж на тот же скид)"),
    "capex20_mult": dict(value=1.3, unit="множитель капекса",
        tier="anchor_repo", source="msa_cost_model: x1.3 на больший серный "
        "узел при 20% конверсии -> full $415/т на факеле"),
    "msa_market_kt_y": dict(value=[50, 70], unit="кт/год",
        tier="anchor_repo", source="skid_twin (риск B1): мировой рынок MSA "
        "~50-70 кт/г — ниша; насыщается 6 скидами при 20% конверсии"),
    "msa_price_usd_t": dict(value=[2500, 2800], unit="$/т",
        tier="market", comment="якорь $2800 — skid_twin; низ $2500 — "
        "скидка новичка/удалёнка; kill при <$1500 (перенасыщение ниши)"),
    "mining_benchmark_usd_t": dict(value={"parity": 323, "x3": 751},
        unit="$/т продукта", tier="anchor_repo",
        source="econ_vs_mining_results: продукт >$751/т бьёт майнинг 3x; "
        "MSA $2500-2800 проходит с запасом"),
    # --- физика ---
    "s_frac_in_so2": dict(value=0.5005, unit="т S/т SO2", tier="anchor_repo",
        source="физика: M(S)/M(SO2) = 32.06/64.06 — не репо-число, твёрже "
        "репо-чисел"),
    # --- допущения о компании (ПРОВЕРИТЬ У КОМПАНИИ) ---
    "so2_capture_mt_y": dict(value=[0.9, 2.5], unit="Мт SO2/год",
        tier="company_assumption",
        comment="масштаб улавливания Серной программы (НМЗ, этапы выхода "
                "на мощность) — проверить у компании"),
    "gypsum_opex_usd_t_s": dict(value=[30, 80], unit="$/т серы",
        tier="company_assumption",
        comment="текущий путь: H2SO4 -> нейтрализация известняком -> гипс; "
                "опекс утилизации вкл. известняк, известь, гипсохранилище — "
                "проверить у компании"),
    "sulfur_transfer_usd_t_s": dict(value=[-50, 0, 50], unit="$/т серы",
        tier="company_assumption",
        comment="трансфертная цена серы (SO3 в пересчёте на S) для MSA-узла: "
                "отрицательная = Серная программа ПЛАТИТ за приём (избавление "
                "от пассива; кредит ограничен избегаемым опексом гипса "
                "$30-80/т) — проверить у компании"),
    "gas_ngp_usd_t": dict(value=[20, 50], unit="$/т CH4",
        tier="company_assumption",
        comment="газ Норильскгазпрома: стрендовый (изолированная система, "
                "альтернативной продажи нет) — проверить у компании"),
    "smp_surcharge_usd_t": dict(value=[80, 200], unit="$/т MSA",
        tier="company_assumption",
        comment="логистическая надбавка танк-контейнер Дудинка -> СМП -> "
                "хаб (Мурманск/Азия) сверх обычной хим-логистики; kill при "
                ">$400 — проверить у компании"),
    # --- литература ---
    "msa_freeze_c": dict(value=20, unit="°C", tier="literature",
        comment="Тпл безводной MSA ~+20°C — в Арктике нужны обогреваемые "
                "танк-контейнеры и склад налива (учтено в оборудовании)"),
}


def V(k):
    return ASSUMPTIONS[k]["value"]


def pick(rng, w):
    """Точка в диапазоне [lo, hi]: w=0 низ, 0.5 середина, 1 верх."""
    return rng[0] + (rng[1] - rng[0]) * w


def msa_cost(s_usd_t, gas_usd_t, conv20=False):
    """Себестоимость т MSA у ворот скида (структура = msa_cost_model)."""
    var = (V("ch4_t_per_t_msa") * gas_usd_t + V("s_t_per_t_msa") * s_usd_t
           + V("initiator_usd_t") + V("utilities_usd_t"))
    if conv20:
        t_y = V("skid20_t_y")
        fixed = V("fixed_skid_usd_t") * V("skid_t_y") / t_y
        capex_musd = V("capex_skid_musd") * V("capex20_mult")
    else:
        t_y = V("skid_t_y")
        fixed = V("fixed_skid_usd_t")
        capex_musd = V("capex_skid_musd")
    capital = capex_musd * 1e6 * V("wacc") / t_y
    return dict(variable=var, cash=var + fixed, full=var + fixed + capital,
                t_y=t_y, capex_musd=capex_musd)


def scenario(s_usd_t, gas_usd_t, log_usd_t, msa_price, skids,
             conv20=False, **_):
    """Полный сценарий: себестоимость, landed, маржа, EBITDA, доля рынка,
    доля серного потока, сравнение тонны серы гипс-vs-MSA."""
    c = msa_cost(s_usd_t, gas_usd_t, conv20)
    landed = c["full"] + log_usd_t
    margin = msa_price - landed
    t_y = skids * c["t_y"]
    ebitda = margin * t_y / 1e6
    capex = skids * c["capex_musd"]
    payback = capex / ebitda if ebitda > 0 else None
    # доля мирового рынка MSA (50-70 кт/г)
    mkt = V("msa_market_kt_y")
    mkt_share = [round(t_y / (mkt[1] * 1000) * 100, 1),
                 round(t_y / (mkt[0] * 1000) * 100, 1)]
    # доля серного потока (0.9-2.5 Мт SO2 -> т S)
    s_used = t_y * V("s_t_per_t_msa")
    cap = V("so2_capture_mt_y")
    s_flow = [cap[0] * 1e6 * V("s_frac_in_so2"),
              cap[1] * 1e6 * V("s_frac_in_so2")]
    s_share = [round(s_used / s_flow[1] * 100, 3),
               round(s_used / s_flow[0] * 100, 3)]
    # тонна серы: гипс (чистый расход) vs MSA (маржа через 0.333 т S/т)
    gypsum = -pick(V("gypsum_opex_usd_t_s"), 0.5)
    per_t_s = margin / V("s_t_per_t_msa")
    return dict(
        full_usd_t=round(c["full"], 0),
        landed_full_usd_t=round(landed, 0),
        margin_usd_t=round(margin, 0),
        msa_t_y=round(t_y, 0), ebitda_musd_y=round(ebitda, 2),
        capex_skids_musd=round(capex, 2),
        payback_y=round(payback, 1) if payback else None,
        world_market_share_pct=mkt_share,
        s_used_t_y=round(s_used, 0),
        s_flow_share_pct=s_share,
        sulfur_per_t=dict(gypsum_usd_t_s=round(gypsum, 0),
                          msa_usd_t_s=round(per_t_s, 0),
                          delta_usd_t_s=round(per_t_s - gypsum, 0)))


# --------------------------------------------- сценарии (pess / base / opt)
# Пессимизм = сера в плюс $50 (программа не делится кредитом), газ по верху,
# СМП по верху, цена MSA по низу, 1 скид. База — середины company-диапазонов,
# сера бесплатно, 3 скида. Оптимизм — сера -$50 (платят за приём), 6 скидов.
SCENARIOS = {
    "pessimistic": dict(s_usd_t=50.0, gas_usd_t=50.0, log_usd_t=200.0,
        msa_price=2500.0, skids=1,
        note="сера +$50, газ $50, СМП $200, цена $2500, 1 скид"),
    "base": dict(s_usd_t=0.0, gas_usd_t=35.0, log_usd_t=140.0,
        msa_price=2650.0, skids=3,
        note="сера $0, газ $35, СМП $140, цена $2650, 3 скида"),
    "optimistic": dict(s_usd_t=-50.0, gas_usd_t=20.0, log_usd_t=80.0,
        msa_price=2800.0, skids=6,
        note="сера -$50 (платят за приём), газ $20, СМП $80, цена $2800, "
             "6 скидов — 9-13% мирового рынка, для новичка агрессивно"),
}


def main():
    # -------- сверка с якорями (обязана сходиться, иначе модель врёт)
    mca = json.load(open(os.path.join(DIR, "msa_cost_model.json")))
    mining = json.load(open(os.path.join(DIR,
                                         "econ_vs_mining_results.json")))
    twin = json.load(open(os.path.join(DIR, "skid_twin.json")))
    c_fl = msa_cost(150.0, 30.0)
    anchors = {
        "msa_full_flare": dict(
            model=round(c_fl["full"], 0),
            anchor=mca["our_skid_flare"]["full"],
            source="msa_cost_model.json: сера $150, газ $30"),
        "msa_full_pipeline": dict(
            model=round(msa_cost(150.0, 180.0)["full"], 0),
            anchor=mca["our_skid_pipeline"]["full"],
            source="msa_cost_model.json: газ $180"),
        "msa_full_flare_20pct": dict(
            model=round(msa_cost(150.0, 30.0, conv20=True)["full"], 0),
            anchor=mca["our_skid_flare_20pct_conv"]["full"],
            source="msa_cost_model.json: конверсия 20%"),
        "capital_usd_t": dict(
            model=round(c_fl["full"] - c_fl["cash"], 0),
            anchor=mca["per_tonne_inputs"]["capital_usd"],
            source="msa_cost_model.json: $4.5M * 12% / 1093 т"),
        "msa_price_hi": dict(
            model=V("msa_price_usd_t")[1],
            anchor=twin["knobs"]["PRODUCTS"]["MSA"]["price"],
            source="skid_twin.json: цена MSA $2800"),
        "mining_3x_usd_t": dict(
            model=V("mining_benchmark_usd_t")["x3"],
            anchor=mining["verdict"]["breakeven_price_to_beat_mining_3x"],
            source="econ_vs_mining_results.json"),
    }
    for a in anchors.values():
        a["ok"] = abs(a["model"] - a["anchor"]) <= max(
            0.005 * abs(a["anchor"]), 1.0)

    scen_out = {name: dict(params=dict(p), **scenario(**p))
                for name, p in SCENARIOS.items()}
    sb = scen_out["base"]
    base = SCENARIOS["base"]

    # -------- апсайд: конверсия 20% (1 скид = 8.7 кт/г = 12-17% рынка)
    upside = dict(scenario(s_usd_t=-50.0, gas_usd_t=20.0, log_usd_t=80.0,
                           msa_price=2500.0, skids=1, conv20=True),
        note="1 скид при 20% конверсии = 8.7 кт/г = 12-17% мирового рынка "
             "— цена в модели уже срезана до $2500, реально поедет ниже "
             "(kill $1500 не случаен); даже при $1500 маржа ~$1073/т. "
             "Это главный приз, но он ЗА гейтом G1 (инициатор/конверсия) "
             "и за вопросом поглощения рынком")

    # -------- чувствительность (вокруг базы): 5 ручек
    def m(**over):
        return scenario(**dict(base, **over))

    sens = {
        "msa_price_usd_t (ГЛАВНАЯ ручка)": {
            "$%d" % p: m(msa_price=p)["margin_usd_t"]
            for p in (1500, 2500, 2650, 2800)},
        "smp_logistics_usd_t (вторая ручка)": {
            "$%d" % p: m(log_usd_t=p)["margin_usd_t"]
            for p in (80, 140, 200, 400)},
        "sulfur_usd_t_s (-50..+150, слабая)": {
            "%+d" % p: m(s_usd_t=p)["margin_usd_t"]
            for p in (-50, 0, 50, 150)},
        "gas_ngp_usd_t (почти не влияет)": {
            "$%d" % p: m(gas_usd_t=p)["margin_usd_t"]
            for p in (20, 35, 50)},
        "skids -> EBITDA $M/г": {
            "%d скид(а)" % n: m(skids=n)["ebitda_musd_y"]
            for n in (1, 3, 6)},
        "conversion -> full $/т у ворот": {
            "2.5%": sb["full_usd_t"],
            "20%": upside["full_usd_t"]},
        "unit": "маржа $/т MSA после логистики (если не сказано иное); "
                "честно: весь размах серы -50..+150 двигает маржу на "
                "~$67/т (0.333 т S/т), газ — на ~$5/т; экономику решают "
                "цена MSA, логистика и конверсия",
    }

    # -------- kill-критерии
    k_price = m(msa_price=1500.0)
    k_both = m(msa_price=1500.0, log_usd_t=400.0)
    kill = {
        "msa_price_below_usd_t": dict(kill_below=1500,
            computed_context="перенасыщение ниши (рынок 50-70 кт/г): при "
                "$1500 маржа базы $%.0f/т, EBITDA 3 скидов $%.2fM/г, "
                "payback %.0f лет — мёртв при конверсии 2.5%%; единственная "
                "защита — конверсия 20%% (full $%.0f/т)"
                % (k_price["margin_usd_t"], k_price["ebitda_musd_y"],
                   k_price["payback_y"], upside["full_usd_t"])),
        "smp_surcharge_above_usd_t": dict(kill_above=400,
            computed_context="сама по себе надбавка $400 маржу не убивает "
                "($%.0f/т в базе), но >$400 означает сломанную логистику "
                "(нет круглогодичного вывоза/ледового класса) — для "
                "спецхимии это потеря клиентов; в комбинации с ценой $1500 "
                "маржа $%.0f/т — мёртв"
                % (m(log_usd_t=400.0)["margin_usd_t"],
                   k_both["margin_usd_t"])),
        "so3_regulatory_ban": dict(kill_if="запрет SO3-узла",
            tier="company_assumption",
            comment="отказ Ростехнадзора/экологической экспертизы в "
                "согласовании узла отбора SO3/олеума на площадке Серной "
                "программы в Норильске — вариант умирает целиком; статус "
                "разрешительной базы проверить у компании"),
        "verdict_scenarios": {},
    }
    for name, s in scen_out.items():
        p = SCENARIOS[name]
        alive = (p["msa_price"] >= 1500 and p["log_usd_t"] <= 400)
        kill["verdict_scenarios"][name] = (
            "жив (маржа $%.0f/т, EBITDA $%.2fM/г)"
            % (s["margin_usd_t"], s["ebitda_musd_y"]) if alive else "МЁРТВ")

    # -------- два пути тонны серы: гипс vs MSA
    sulfur_two_paths = {
        "gypsum_usd_t_s": dict(value=[-80, -30],
            note="чистый расход: известняк + гипсохранилище "
                 "(company_assumption, проверить)"),
        "msa_usd_t_s": {name: s["sulfur_per_t"]["msa_usd_t_s"]
                        for name, s in scen_out.items()},
        "delta_usd_t_s": {name: s["sulfur_per_t"]["delta_usd_t_s"]
                          for name, s in scen_out.items()},
        "honest": "дельта ~$3.5-4.9k/т S красивая, но применима только к "
                  "0.03-0.5% серного потока (см. market_ceiling) — MSA "
                  "НЕ решает Серную программу, а монетизирует срез пассива; "
                  "остальные 99.5%+ серы по-прежнему уходят в гипс",
    }

    # -------- потолок рынка B1 (обязательный, честный)
    def ceiling_row(label, t_y, skids_note):
        s_used = t_y * V("s_t_per_t_msa")
        cap = V("so2_capture_mt_y")
        flow = [c * 1e6 * V("s_frac_in_so2") for c in cap]
        return dict(config=label, msa_kt_y=round(t_y / 1000.0, 1),
            world_market_share_pct=[round(t_y / 70000.0 * 100, 1),
                                    round(t_y / 50000.0 * 100, 1)],
            s_used_t_y=round(s_used, 0),
            s_flow_share_pct=[round(s_used / flow[1] * 100, 3),
                              round(s_used / flow[0] * 100, 3)],
            note=skids_note)

    t1 = V("skid_t_y")
    market_ceiling = {
        "world_msa_kt_y": V("msa_market_kt_y"),
        "rows": [
            ceiling_row("1 скид @2.5%", 1 * t1, "пилот — рынок не заметит"),
            ceiling_row("3 скида @2.5%", 3 * t1,
                        "база: ~5-7% рынка — предел комфорта новичка"),
            ceiling_row("6 скидов @2.5%", 6 * t1,
                        "9-13% рынка — цена начнёт ехать"),
            ceiling_row("1 скид @20% конв.", V("skid20_t_y"),
                        "тот же тоннаж потолка меньшим числом скидов; "
                        "12-17% рынка — потолок ниши, дальше только "
                        "рост самого рынка (электроника/электролиты)"),
        ],
        "honest": "утилизация серы — ДОЛИ ПРОЦЕНТА потока (0.03-0.65%): "
                  "заголовок 'MSA спасает Серную программу' был бы ложью. "
                  "Смысл кейса — маржа на нише + витрина технологии, не "
                  "экология масштаба Мт",
    }

    # -------- механическая реализация (скрин)
    implementation = {
        "block_flow": [
            "металлургические газы (НМЗ) -> Серная программа: очистка, "
            "конверсия SO2->SO3 в контактном отделении, абсорбция в H2SO4, "
            "нейтрализация известняком -> гипс (СУЩЕСТВУЕТ, без изменений)",
            "НОВОЕ: узел отбора SO3 (олеум/десорбция) на боковом потоке "
            "<1% мощности контактного отделения",
            "газ Норильскгазпрома: осушка + очистка CH4 (стренд, $20-50/т)",
            "реакторный скид: радикальное сульфонирование CH4+SO3 под "
            "давлением, рецикл CH4, дозирование инициатора",
            "дистилляция/концентрирование -> MSA 70% / 99.5% (безводная)",
            "обогреваемый налив в нерж. танк-контейнеры (Тпл MSA ~+20°C)",
            "Дудинка -> СМП -> хаб (Мурманск/Азия) -> клиенты "
            "(электроника, электролиты, биодизель-катализ)",
        ],
        "equipment": [
            dict(item="узел отбора SO3/олеума из контактного отделения",
                 cost_musd=[0.8, 1.5]),
            dict(item="3 реакторных скида MSA по 1.1 кт/г",
                 cost_musd=[13.5, 13.5],
                 note="якорь $4.5M/скид (msa_cost_model); в payback "
                      "сценариев учтён только этот капекс"),
            dict(item="парк обогреваемых нерж. танк-контейнеров 40-80 шт "
                      "(оборот СМП 30-60 дней)", cost_musd=[1.5, 3.5]),
            dict(item="обогреваемый склад налива/накопления партий",
                 cost_musd=[0.5, 1.0]),
            dict(item="лаборатория качества (99.5% для электроники)",
                 cost_musd=[0.3, 0.6]),
        ],
        "equipment_total_musd": [16.6, 20.1],
        "equipment_note": "инфра сверх скидов добавляет ~20-30% к капексу "
                          "сценариев — паybackи удлиняются пропорционально",
        "footprint_m2": [1200, 2500],
        "footprint_note": "3 скида + SO3-отбор + склад налива; на "
                          "существующей промплощадке Серной программы",
        "utilities": {
            "electricity_kw": [900, 1800],
            "steam": "дистилляция MSA — главный потребитель; интеграция с "
                     "теплом контактного отделения (экзотермика SO2->SO3)",
            "cooling_water_m3_h": [30, 90],
            "nitrogen_nm3_h": [30, 60],
            "ch4_t_y": 548,
            "note": "3 скида; CH4 = 3 x 1093 x 0.167 т/г от НГП",
        },
        "staffing": dict(
            pilot="4-6 чел: технолог, инженер-химик, 2-3 оператора "
                  "(вахта), логист СМП",
            series="10-14 чел на 3 скида вкл. логистику",
            tier="screening"),
        "rollout": [
            dict(stage="FEED + HAZOP SO3-отбора + трансфертная цена серы "
                       "с компанией", duration="6-9 мес"),
            dict(stage="лаборатория: G1 (инициатор) + G2 (примеси SO3)",
                 duration="6-12 мес", note="параллельно FEED"),
            dict(stage="пилотный скид 1.1 кт/г + первая навигация СМП",
                 duration="12-18 мес"),
            dict(stage="серия до 3 скидов", duration="12-24 мес"),
            dict(stage="апгрейд конверсии до 20% (после G1 и офф-тейка)",
                 duration="12 мес"),
        ],
        "trl_gates": [
            "G1 (расчёт, наш квантовый рычаг): дизайн инициатора "
            "радикальной цепи CH4+SO3 — строка INITIATOR $100/т -> цель "
            "<$50/т (msa_cost_model: болевая точка) + путь к конверсии "
            "20%; NEVPT2/DFT-стек",
            "G2 (лаборатория): совместимость SO3/олеума Серной программы "
            "по примесям (пыль, As, F, Cl, туман H2SO4) с радикальным "
            "сульфонированием — анализ реального потока",
            "G3 (пилот): 1000+ ч пробега скида на реальном SO3 и газе "
            "НГП: full cost <= $1300/т по счётчикам, качество 99.5%",
            "G4 (логистика): полный цикл Дудинка->СМП->хаб: надбавка "
            "<= $200/т подтверждена, качество сохранено (обогрев)",
            "G5 (рынок): офф-тейк >= 1 кт/г по >= $2500/т до серии; для "
            "апгрейда 20% — офф-тейк >= 5 кт/г",
        ],
    }

    out = {
        "model": "Норникель N2: Серная программа -> MSA (SO3 по "
                 "отрицательной цене), вывоз Дудинка -> СМП",
        "tier": "screening_pm_2_3x",
        "not_bankable": True,
        "tier_legend": {
            "anchor_repo": "число из наших файлов или физика (помечено)",
            "market": "рыночный диапазон, порядок величины",
            "company_assumption": "допущение о компании — проверить у неё",
            "literature": "литература/индустрия, не наш расчёт"},
        "anchors_check": anchors,
        "assumptions": ASSUMPTIONS,
        "scenarios": scen_out,
        "upside_20pct_conversion": upside,
        "sensitivity": sens,
        "kill_criteria": kill,
        "sulfur_two_paths": sulfur_two_paths,
        "market_ceiling": market_ceiling,
        "implementation": implementation,
        "honesty": "скрин +-2-3x. Честные выводы: (1) 'сера по "
                   "отрицательной цене' — красивый заголовок, но слабый "
                   "рычаг: весь размах -50..+150 $/т S двигает маржу лишь "
                   "на ~$67/т MSA (0.333 т S/т), газ — на ~$5/т; экономику "
                   "решают ЦЕНА MSA (маржа $230..1530/т на $1500..2800), "
                   "ЛОГИСТИКА СМП и КОНВЕРСИЯ; (2) утилизация серы — доли "
                   "процента потока (0.03-0.65%): MSA монетизирует срез "
                   "пассива, а не решает Серную программу; (3) база (3 "
                   "скида, $2650) EBITDA ~$4.5M/г, payback ~3 г — "
                   "сопоставимо с NaCN-кейсом ($4.34M, 1.4 г), но "
                   "логистика тяжелее и рынок уже; (4) MSA $2500-2800 >> "
                   "$751 (порог 3x майнинга) — тест на specialty пройден; "
                   "(5) реальный приз — конверсия 20% (full $347/т landed "
                   "$427), но он за гейтом G1 и упирается в потолок ниши "
                   "12-17% рынка.",
    }

    # ------------------------------------------------------------- сводка
    print("=" * 74)
    print("Норникель N2: Серная программа -> MSA (SO3 по отрицательной "
          "цене)")
    print("   скрин +-2-3x, не банк; вывоз Дудинка -> СМП")
    print("=" * 74)
    print("сверка якорей:", "OK" if all(a["ok"] for a in anchors.values())
          else "РАСХОЖДЕНИЕ! " + str(anchors))
    print("\n%-13s%7s%8s%8s%8s%9s%9s%7s" % ("сценарий", "full", "landed",
          "маржа", "т/г", "EBITDA", "capex", "окуп."))
    for name, s in scen_out.items():
        print("%-13s%7.0f%8.0f%8.0f%8.0f%8.2fM%8.2fM%6.1fг"
              % (name, s["full_usd_t"], s["landed_full_usd_t"],
                 s["margin_usd_t"], s["msa_t_y"], s["ebitda_musd_y"],
                 s["capex_skids_musd"], s["payback_y"]))
    print("апсайд 20%% конв.: full $%.0f, landed $%.0f, маржа $%.0f "
          "(цена уже $2500), 8.7 кт/г"
          % (upside["full_usd_t"], upside["landed_full_usd_t"],
             upside["margin_usd_t"]))
    print("\nтонна серы: гипс -$30..80 vs MSA +$%.0f (база) -> дельта "
          "~$%.0f/т S,\n  но лишь для %.2f-%.2f%% серного потока (0.9-2.5 "
          "Мт SO2/г)"
          % (sb["sulfur_per_t"]["msa_usd_t_s"],
             sb["sulfur_per_t"]["delta_usd_t_s"],
             sb["s_flow_share_pct"][0], sb["s_flow_share_pct"][1]))
    print("потолок B1: рынок MSA 50-70 кт/г; 3 скида = %.1f-%.1f%% рынка, "
          "6 = %.1f-%.1f%%"
          % (sb["world_market_share_pct"][0],
             sb["world_market_share_pct"][1],
             scen_out["optimistic"]["world_market_share_pct"][0],
             scen_out["optimistic"]["world_market_share_pct"][1]))
    print("\nчувствительность (маржа $/т, база):")
    for knob, vals in sens.items():
        if knob == "unit":
            continue
        row = "  ".join("%s:%s" % (k, v) for k, v in vals.items())
        print("  %-38s %s" % (knob, row))
    print("\nkill: цена MSA<$1500 | СМП-надбавка>$400 | запрет SO3-узла "
          "в Норильске")
    print("вердикт:", kill["verdict_scenarios"])

    path = os.path.join(DIR, "econ_nornickel_sulfur_msa_results.json")
    with open(path, "w") as f:
        json.dump(out, f, indent=1, ensure_ascii=False)
    print("\nwrote %s" % os.path.basename(path))


if __name__ == "__main__":
    main()
