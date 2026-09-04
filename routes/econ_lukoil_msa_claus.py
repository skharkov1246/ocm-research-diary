#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
routes/econ_lukoil_msa_claus.py — Лукойл L4: MSA-скид на НПЗ — сера Клауса +
трубный газ + обслуживаемая площадка. Идея: на НПЗ уже есть установка Клауса
(H2S -> жидкая сера, продаётся по $80-150/т на месте), заводской трубный газ
и ШТАТ — MSA-скид (CH4 + SO3 -> CH3SO3H) ставится в обслуживаемый периметр,
логистика — обычная ж/д/авто РФ (не СМП!), сбыт — импортозамещение MSA в РФ
($2500-3200/т, импортные паритеты) + экспорт.

Себестоимость тонны MSA — ровно модель msa_cost_model.py (якоря: 0.167 т CH4
+ 0.333 т S на тонну, инициатор $100, утилиты $130, газ трубный $180/т,
WACC 12%); отличия от якоря — ЯВНЫЕ и в минус себестоимости:
  fixed $300-400/т (не $400): площадка обслуживаемая — операторы/ремонт/ЦЗЛ
    делятся со штатом НПЗ, риск B2 (глухой скид без людей) снят по построению;
  capex $3.5-4.5M (не $4.5M): синергия НПЗ — SO3 из существующего
    сернокислотного/Клаус-хозяйства, пар и энергия заводские (E1 ±2-3x).
Газ трубный $180/т — якорь; в MSA газ <7%% себестоимости, нечувствительность
показана явно. Каннибализации нет (Лукойл MSA не производит), но есть
конкуренция с N2-скидами НАШЕЙ платформы (econ_nornickel_sulfur_msa) за ту же
нишу 50-70 кт/г — потолок посчитан СОВМЕСТНО (shared_niche).

Скрин-уровень ±2-3x, НЕ банковская модель. Цифры по Лукойлу (цена серы Клауса
на месте, наличие сернокислотного хозяйства, доля штата) — ТОЛЬКО диапазоны с
tier='company_assumption' и пометкой «проверить у компании».

Запуск: python3 routes/econ_lukoil_msa_claus.py
     -> routes/econ_lukoil_msa_claus_results.json
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
        comment="на НПЗ пар/энергия заводские и SO3 может приходить готовым "
                "из сернокислотного хозяйства (минус), но покупные тарифы "
                "НПЗ не нулевые (плюс) — якорь оставлен как есть, скрин"),
    "gas_pipeline_usd_t": dict(value=180.0, unit="$/т CH4",
        tier="anchor_repo", source="msa_cost_model our_skid_pipeline: газ "
        "трубный $180/т",
        comment="в MSA газа всего 0.167 т/т -> $30/т MSA = <3% full: "
                "нечувствительность к цене газа показана в sensitivity"),
    "fixed_skid_usd_t": dict(value=[300.0, 400.0], unit="$/т MSA",
        tier="company_assumption",
        comment="якорь msa_cost_model $400/т — труд/обслуживание/накладные "
                "ГЛУХОГО скида малого масштаба; на обслуживаемой площадке "
                "НПЗ операторы (вахты уже есть), ремонт и ЦЗЛ делятся со "
                "штатом завода -> $300-400/т; риск B2 снят по построению; "
                "долю разделяемого штата проверить у компании"),
    "capex_skid_musd": dict(value=[3.5, 4.5], unit="$M/скид",
        tier="company_assumption",
        comment="якорь msa_cost_model $4.5M — скид со своим SO3-узлом в "
                "чистом поле; на НПЗ SO3 из существующего сернокислотного/"
                "Клаус-хозяйства, пар/энергия/эстакады заводские -> "
                "$3.5-4.5M (E1 ±2-3x); наличие сернокислотного хозяйства "
                "на конкретном НПЗ проверить у компании"),
    "skid_t_y": dict(value=1093.0, unit="т MSA/год/скид", tier="anchor_repo",
        source="msa_cost_model: конверсия CH4->MSA ~2.5% per-pass suite"),
    "wacc": dict(value=0.12, unit="доля/год", tier="anchor_repo",
        source="skid_twin / econ_shock_candidates: 12% (аморт.+WACC)"),
    "skid20_t_y": dict(value=8747.0, unit="т MSA/год/скид",
        tier="anchor_repo", source="msa_cost_model: апсайд конверсии 20% "
        "(x8 тоннаж на тот же скид)"),
    "capex20_mult": dict(value=1.3, unit="множитель капекса",
        tier="anchor_repo", source="msa_cost_model: x1.3 на больший серный "
        "узел при 20% конверсии"),
    "msa_market_kt_y": dict(value=[50, 70], unit="кт/год",
        tier="anchor_repo", source="skid_twin (риск B1): мировой рынок MSA "
        "~50-70 кт/г — ниша; насыщается 6 скидами при 20% конверсии"),
    "msa_price_world_usd_t": dict(value=2800.0, unit="$/т",
        tier="anchor_repo", source="skid_twin: мировая цена MSA $2800"),
    "mining_benchmark_usd_t": dict(value={"parity": 323, "x3": 751},
        unit="$/т продукта", tier="anchor_repo",
        source="econ_vs_mining_results: продукт >$751/т бьёт майнинг 3x; "
        "MSA $2500-3200 проходит с запасом"),
    "n2_skids_base": dict(value=3, unit="скидов N2 (Норильск)",
        tier="anchor_repo", source="econ_nornickel_sulfur_msa: база 3 скида "
        "— ОДНА ниша с L4, потолок считаем совместно (shared_niche)"),
    # --- рынок ---
    "msa_price_rf_usd_t": dict(value=[2500, 3200], unit="$/т",
        tier="market",
        comment="цена MSA внутри РФ как импортозамещение: импортный паритет "
                "= мировая $2800 (якорь) + фрахт/пошлина/НДС-логистика "
                "импортёра; низ $2500 — скидка новичка/экспортный паритет, "
                "верх $3200 — премия локального поставщика с коротким "
                "плечом; kill при <$1500 (перенасыщение ниши)"),
    "rf_msa_import_kt_y": dict(value=[3, 10], unit="кт/год",
        tier="market",
        comment="внутренний спрос РФ на MSA (гальваника/электроника, "
                "электролиты, биодизель-катализ) — порядок величины по "
                "импортозамещению; сверить с таможенной статистикой: "
                "премию $3200 несёт только этот срез, остальное — экспорт "
                "по мировой цене"),
    "logistics_rf_usd_t": dict(value=[30, 80], unit="$/т MSA",
        tier="market",
        comment="обычная хим-логистика РФ: ж/д цистерны (нерж., обогрев) / "
                "автоцистерны с НПЗ до клиента — НЕ СМП; преимущество "
                "против импорта (у импортёра фрахт+пошлина в цене)"),
    # --- допущения о компании (ПРОВЕРИТЬ У КОМПАНИИ) ---
    "sulfur_claus_usd_t": dict(value=[80.0, 150.0], unit="$/т S",
        tier="company_assumption",
        comment="сера Клауса на месте: на НПЗ обычно ДЕШЕВЛЕ рынка "
                "($150 комодити-якорь msa_cost_model = наш верх) — нет "
                "логистики вывоза, сбыт серы бывает проблемным активом; "
                "трансфертную цену проверить у компании"),
    "so3_source_onsite": dict(value="сернокислотное/Клаус-хозяйство НПЗ",
        tier="company_assumption",
        comment="ключ синергии: SO3/олеум из существующего контактного "
                "производства H2SO4 (напр. регенерация кислоты алкилации) "
                "ИЛИ компактное сжигание серы Клауса при скиде; что именно "
                "есть на конкретном НПЗ — проверить у компании (двигает "
                "capex внутри $3.5-4.5M)"),
    # --- литература ---
    "msa_freeze_c": dict(value=20, unit="°C", tier="literature",
        comment="Тпл безводной MSA ~+20°C — обогреваемые цистерны/склад "
                "налива; в средней полосе РФ это стандартная хим-практика, "
                "не арктическое исполнение (в отличие от N2)"),
}


def V(k):
    return ASSUMPTIONS[k]["value"]


def msa_cost(s_usd_t, gas_usd_t, fixed_usd_t=400.0, capex_musd=4.5,
             conv20=False):
    """Себестоимость т MSA у ворот скида (структура = msa_cost_model);
    fixed и capex — ручки L4 (синергия НПЗ), дефолты = якорь."""
    var = (V("ch4_t_per_t_msa") * gas_usd_t + V("s_t_per_t_msa") * s_usd_t
           + V("initiator_usd_t") + V("utilities_usd_t"))
    if conv20:
        t_y = V("skid20_t_y")
        fixed = fixed_usd_t * V("skid_t_y") / t_y
        capex = capex_musd * V("capex20_mult")
    else:
        t_y = V("skid_t_y")
        fixed = fixed_usd_t
        capex = capex_musd
    capital = capex * 1e6 * V("wacc") / t_y
    return dict(variable=var, cash=var + fixed, full=var + fixed + capital,
                t_y=t_y, capex_musd=capex)


def scenario(s_usd_t, gas_usd_t, fixed_usd_t, capex_musd, log_usd_t,
             msa_price, skids, conv20=False, **_):
    """Полный сценарий: себестоимость, landed, маржа, EBITDA, доля рынка
    (совместно с N2-скидами платформы), доля газа в себестоимости."""
    c = msa_cost(s_usd_t, gas_usd_t, fixed_usd_t, capex_musd, conv20)
    landed = c["full"] + log_usd_t
    margin = msa_price - landed
    t_y = skids * c["t_y"]
    ebitda = margin * t_y / 1e6
    capex = skids * c["capex_musd"]
    payback = capex / ebitda if ebitda > 0 else None
    mkt = V("msa_market_kt_y")
    mkt_share = [round(t_y / (mkt[1] * 1000) * 100, 1),
                 round(t_y / (mkt[0] * 1000) * 100, 1)]
    # платформа целиком: L4 + N2 (база 3 скида) в одной нише
    plat_t_y = t_y + V("n2_skids_base") * V("skid_t_y")
    plat_share = [round(plat_t_y / (mkt[1] * 1000) * 100, 1),
                  round(plat_t_y / (mkt[0] * 1000) * 100, 1)]
    gas_share = round(V("ch4_t_per_t_msa") * gas_usd_t / c["full"] * 100, 1)
    return dict(
        full_usd_t=round(c["full"], 0),
        landed_full_usd_t=round(landed, 0),
        margin_usd_t=round(margin, 0),
        msa_t_y=round(t_y, 0), ebitda_musd_y=round(ebitda, 2),
        capex_skids_musd=round(capex, 2),
        payback_y=round(payback, 1) if payback else None,
        world_market_share_pct=mkt_share,
        platform_share_with_n2_pct=plat_share,
        gas_share_of_full_pct=gas_share)


# --------------------------------------------- сценарии (pess / base / opt)
# Пессимизм = сера по комодити $150 (НПЗ не делится скидкой), fixed и capex
# по якорю (синергии не дали), логистика по верху, цена $2500, 1 скид —
# ровно our_skid_pipeline из msa_cost_model + ж/д. База — середины
# company-диапазонов, 2 скида. Оптимизм — сера $80, полная синергия НПЗ,
# премия импортозамещения $3200, 4 скида.
SCENARIOS = {
    "pessimistic": dict(s_usd_t=150.0, gas_usd_t=180.0, fixed_usd_t=400.0,
        capex_musd=4.5, log_usd_t=80.0, msa_price=2500.0, skids=1,
        note="сера $150, синергий нет (fixed $400, capex $4.5M), ж/д $80, "
             "цена $2500, 1 скид — это ровно our_skid_pipeline якоря"),
    "base": dict(s_usd_t=115.0, gas_usd_t=180.0, fixed_usd_t=350.0,
        capex_musd=4.0, log_usd_t=55.0, msa_price=2850.0, skids=2,
        note="сера $115, fixed $350, capex $4.0M, ж/д $55, цена $2850 "
             "(середина импортного паритета), 2 скида"),
    "optimistic": dict(s_usd_t=80.0, gas_usd_t=180.0, fixed_usd_t=300.0,
        capex_musd=3.5, log_usd_t=30.0, msa_price=3200.0, skids=4,
        note="сера $80, полная синергия (fixed $300, capex $3.5M), авто "
             "$30, премия импортозамещения $3200, 4 скида — премию $3200 "
             "несёт только внутренний спрос РФ 3-10 кт/г, дальше экспорт "
             "по ~$2500-2800"),
}


def main():
    # -------- сверка с якорями (обязана сходиться, иначе модель врёт)
    mca = json.load(open(os.path.join(DIR, "msa_cost_model.json")))
    mining = json.load(open(os.path.join(DIR,
                                         "econ_vs_mining_results.json")))
    twin = json.load(open(os.path.join(DIR, "skid_twin.json")))
    c_pl = msa_cost(150.0, 180.0)
    anchors = {
        "msa_full_pipeline": dict(
            model=round(c_pl["full"], 0),
            anchor=mca["our_skid_pipeline"]["full"],
            source="msa_cost_model.json: сера $150, газ трубный $180 — "
                   "БЛИЖАЙШИЙ якорь L4 (= наш pessimistic до логистики)"),
        "msa_variable_pipeline": dict(
            model=round(c_pl["variable"], 0),
            anchor=mca["our_skid_pipeline"]["variable"],
            source="msa_cost_model.json: variable на трубном газе"),
        "msa_full_flare": dict(
            model=round(msa_cost(150.0, 30.0)["full"], 0),
            anchor=mca["our_skid_flare"]["full"],
            source="msa_cost_model.json: сера $150, газ $30"),
        "msa_full_flare_20pct": dict(
            model=round(msa_cost(150.0, 30.0, conv20=True)["full"], 0),
            anchor=mca["our_skid_flare_20pct_conv"]["full"],
            source="msa_cost_model.json: конверсия 20%"),
        "capital_usd_t": dict(
            model=round(c_pl["full"] - c_pl["cash"], 0),
            anchor=mca["per_tonne_inputs"]["capital_usd"],
            source="msa_cost_model.json: $4.5M * 12% / 1093 т"),
        "msa_price_world": dict(
            model=V("msa_price_world_usd_t"),
            anchor=twin["knobs"]["PRODUCTS"]["MSA"]["price"],
            source="skid_twin.json: мировая цена MSA $2800 — внутри нашего "
                   "РФ-диапазона 2500-3200 (импортный паритет)"),
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

    # -------- апсайд: конверсия 20% (тот же гейт G1, что у N2)
    upside = dict(scenario(s_usd_t=80.0, gas_usd_t=180.0, fixed_usd_t=300.0,
                           capex_musd=3.5, log_usd_t=30.0, msa_price=2500.0,
                           skids=1, conv20=True),
        note="1 скид при 20% конверсии = 8.7 кт/г = 12-17% мирового рынка; "
             "full у ворот ~$387/т — тот же порядок, что у якоря N2 "
             "($347/т на отрицательной сере и газе $20: разница = сера "
             "$80 vs -$50 и газ $180 vs $20, вместе ~$70/т); цена в "
             "модели уже срезана до $2500 — рынок насыщается 6 скидами "
             "ПЛАТФОРМЫ (L4+N2 вместе), приз за гейтом G1 (инициатор/"
             "конверсия) делится между площадками, а не удваивается")

    # -------- чувствительность (вокруг базы): 6 ручек
    def m(**over):
        return scenario(**dict(base, **over))

    sens = {
        "msa_price_usd_t (ГЛАВНАЯ ручка)": {
            "$%d" % p: m(msa_price=p)["margin_usd_t"]
            for p in (1500, 2500, 2850, 3200)},
        "sulfur_claus_usd_t_s (вторая по смыслу, слабая по $)": {
            "$%d" % p: m(s_usd_t=p)["margin_usd_t"]
            for p in (80, 115, 150, 250)},
        "logistics_rf_usd_t (ж/д/авто, вся вилка ~$50/т)": {
            "$%d" % p: m(log_usd_t=p)["margin_usd_t"]
            for p in (30, 55, 80, 150)},
        "gas_pipeline_usd_t (ПОКАЗАТЕЛЬНО нечувствительна)": {
            "$%d" % p: m(gas_usd_t=p)["margin_usd_t"]
            for p in (30, 100, 180, 300)},
        "fixed_skid_usd_t (обслуживаемая площадка)": {
            "$%d" % p: m(fixed_usd_t=p)["margin_usd_t"]
            for p in (300, 350, 400)},
        "capex_musd -> payback лет (E1 +-2-3x)": {
            "$%.1fM" % c: m(capex_musd=c)["payback_y"]
            for c in (3.5, 4.0, 4.5, 9.0)},
        "skids -> EBITDA $M/г": {
            "%d скид(а)" % n: m(skids=n)["ebitda_musd_y"]
            for n in (1, 2, 4)},
        "conversion -> full $/т у ворот": {
            "2.5%": sb["full_usd_t"],
            "20%": upside["full_usd_t"]},
        "unit": "маржа $/т MSA после логистики (если не сказано иное); "
                "честно: газ $30->$300 двигает маржу лишь на ~$45/т "
                "(0.167 т/т; в базе газ %.1f%% full — «трубный vs факел» "
                "для MSA почти неважно), сера $80->$250 — на ~$57/т; "
                "экономику решают цена MSA, конверсия и (для payback) "
                "капекс" % sb["gas_share_of_full_pct"],
    }

    # -------- kill-критерии
    k_price = m(msa_price=1500.0)
    kill = {
        "msa_price_below_usd_t": dict(kill_below=1500,
            computed_context="перенасыщение ниши (рынок 50-70 кт/г, и в ней "
                "уже сидят N2-скиды нашей же платформы): при $1500 маржа "
                "базы $%.0f/т, EBITDA 2 скидов $%.2fM/г, payback %.0f лет "
                "— мёртв при конверсии 2.5%%; единственная защита — "
                "конверсия 20%% (full $%.0f/т)"
                % (k_price["margin_usd_t"], k_price["ebitda_musd_y"],
                   k_price["payback_y"], upside["full_usd_t"])),
        "so3_regulatory_ban": dict(kill_if="запрет SO3-узла на НПЗ",
            tier="company_assumption",
            comment="отказ Ростехнадзора/промбезопасности в согласовании "
                "SO3/олеум-узла в периметре действующего НПЗ (ОПО I "
                "класса, изменение декларации ПБ) — вариант умирает "
                "целиком; статус разрешительной базы и позицию завода "
                "проверить у компании"),
        "platform_niche_saturated": dict(kill_above_skids=10,
            computed_context="ниша занята: >10 скидов ПЛАТФОРМЫ (L4+N2 "
                "суммарно, @2.5%%) = >%.1f кт/г = >%.0f-%.0f%% мирового "
                "рынка — цена едет к $1500 и ниже; учитывать скиды N2 "
                "(база 3) при планировании флота L4"
                % (10 * V("skid_t_y") / 1000,
                   10 * V("skid_t_y") / 700, 10 * V("skid_t_y") / 500)),
        "verdict_scenarios": {},
    }
    for name, s in scen_out.items():
        p = SCENARIOS[name]
        plat_skids = p["skids"] + V("n2_skids_base")
        alive = (p["msa_price"] >= 1500 and plat_skids <= 10)
        kill["verdict_scenarios"][name] = (
            "жив (маржа $%.0f/т, EBITDA $%.2fM/г, платформа %d скидов)"
            % (s["margin_usd_t"], s["ebitda_musd_y"], plat_skids)
            if alive else "МЁРТВ")

    # -------- совместная ниша с N2 (обязательный, честный)
    t1 = V("skid_t_y")

    def niche_row(label, t_y, note):
        return dict(config=label, msa_kt_y=round(t_y / 1000.0, 1),
            world_market_share_pct=[round(t_y / 70000.0 * 100, 1),
                                    round(t_y / 50000.0 * 100, 1)],
            note=note)

    shared_niche = {
        "world_msa_kt_y": V("msa_market_kt_y"),
        "who_shares": "L4 (НПЗ Лукойла, эта модель) и N2 (Норильск, "
                      "econ_nornickel_sulfur_msa) — ДВЕ площадки одной "
                      "платформы продают в ОДНУ нишу 50-70 кт/г; "
                      "каннибализации продуктов Лукойла нет (MSA он не "
                      "делает), но каннибализация ВНУТРИ платформы есть",
        "rows": [
            niche_row("L4 база (2) + N2 база (3) = 5 скидов @2.5%",
                      5 * t1, "~5.5 кт/г — комфорт, рынок почти не "
                      "заметит; L4 при этом закрывает импорт РФ 3-10 кт/г "
                      "лишь частично"),
            niche_row("L4 опт (4) + N2 опт (6) = 10 скидов @2.5%",
                      10 * t1, "~11 кт/г = 16-22% рынка — ПОТОЛОК "
                      "платформы, дальше цена едет (kill >10 скидов)"),
            niche_row("1 скид @20% (любая площадка)", V("skid20_t_y"),
                      "8.7 кт/г одним скидом = 12-17% рынка — после G1 "
                      "весь флот 2.5% устаревает"),
            niche_row("6 скидов @20% (вся платформа)",
                      6 * V("skid20_t_y"), "52 кт/г = 75-105% рынка — "
                      "полное насыщение: L4 и N2 НЕ могут одновременно "
                      "гнать флоты на 20%; реальный предел платформы "
                      "@20% — 2-3 скида на обе площадки"),
        ],
        "shared_niche_note": "рынок 50-70 кт/г выдерживает СОВОКУПНО ~10 "
            "скидов платформы @2.5% конверсии (16-22% рынка — верх "
            "комфорта новичка) или ~2-3 скида @20% (уже 25-52% рынка); "
            "L4 и N2 конкурируют за одну нишу, флоты планировать "
            "суммарно; преимущество L4 в этой конкуренции — дешёвая "
            "ж/д-логистика ($30-80 vs $80-200 СМП) и премия "
            "импортозамещения РФ, преимущество N2 — сера по нулевой/"
            "отрицательной цене и дешёвый газ (~$100/т full в его "
            "пользу); честный порядок: N2 берёт мировой рынок, L4 — "
            "внутренний РФ (3-10 кт/г = 1-2 скида L4 без давления на "
            "цену)",
        "rf_domestic_slice": dict(
            demand_kt_y=V("rf_msa_import_kt_y"),
            note="премию $3200 несёт только внутренний спрос РФ (импорто"
                 "замещение); 1 скид = 1.1 кт/г закрывает 11-36% импорта, "
                 "2 скида — 22-73%: оптимистичные 4 скида уже наполовину "
                 "экспортные по ~$2500-2800"),
    }

    # -------- механическая реализация (скрин)
    implementation = {
        "block_flow": [
            "НПЗ: аминовая очистка -> кислый газ H2S -> установка Клауса "
            "-> жидкая/гранулированная сера (СУЩЕСТВУЕТ, без изменений)",
            "НОВОЕ: SO3-узел — вариант (а): отбор олеума/SO3 из "
            "существующего сернокислотного хозяйства НПЗ (регенерация "
            "кислоты алкилации и т.п.); вариант (б): компактное сжигание "
            "серы Клауса + контактный аппарат при скиде",
            "трубный/заводской газ: осушка + очистка CH4 ($180/т, <3% "
            "себестоимости)",
            "реакторный скид: радикальное сульфонирование CH4+SO3 под "
            "давлением, рецикл CH4, дозирование инициатора",
            "дистилляция/концентрирование -> MSA 70% / 99.5% (безводная)",
            "обогреваемый налив (Тпл MSA ~+20°C) в ж/д цистерны (нерж.) / "
            "автоцистерны — обычная хим-логистика РФ, не СМП",
            "клиенты: импортозамещение РФ (гальваника/электроника, "
            "электролиты ХИТ, биодизель-катализ) + экспорт",
        ],
        "equipment": [
            dict(item="SO3-узел: отбор олеума из кислотного хозяйства ИЛИ "
                      "сжигание серы + компактный контактный аппарат",
                 cost_musd=[0.6, 1.8],
                 note="вариант (а) дешевле (б); что доступно — зависит от "
                      "конкретного НПЗ (company_assumption)"),
            dict(item="2 реакторных скида MSA по 1.1 кт/г",
                 cost_musd=[7.0, 9.0],
                 note="$3.5-4.5M/скид (синергия НПЗ против якоря $4.5M); "
                      "в payback сценариев учтён только этот капекс"),
            dict(item="обогреваемый склад налива + эстакада ж/д/авто "
                      "(цистерны — аренда парка)", cost_musd=[0.4, 0.8]),
            dict(item="лаборатория качества (99.5% для электроники)",
                 cost_musd=[0.2, 0.4],
                 note="частично на базе ЦЗЛ НПЗ — ещё одна синергия "
                      "обслуживаемой площадки"),
        ],
        "equipment_total_musd": [8.2, 12.0],
        "equipment_note": "инфра сверх скидов добавляет ~15-30% к капексу "
                          "сценариев — paybackи удлиняются пропорционально; "
                          "сравнить с N2: там же +20-30%, но плюс парк "
                          "арктических танк-контейнеров",
        "footprint_m2": [600, 1500],
        "footprint_note": "2 скида + SO3-узел + налив; на действующей "
                          "промплощадке НПЗ рядом с Клаусом/кислотным "
                          "хозяйством — земли и ограждающего периметра "
                          "не требуется",
        "utilities": {
            "electricity_kw": [600, 1200],
            "steam": "дистилляция MSA — главный потребитель; пар из "
                     "заводской сети НПЗ (синергия; при варианте (б) — "
                     "плюс экзотермика сжигания серы)",
            "cooling_water_m3_h": [20, 60],
            "nitrogen_nm3_h": [20, 40],
            "ch4_t_y": 365,
            "note": "2 скида; CH4 = 2 x 1093 x 0.167 т/г из заводской "
                    "сети — для НПЗ незаметно",
        },
        "staffing": dict(
            pilot="2-4 чел ИНКРЕМЕНТАЛЬНО: технолог + инженер-химик; "
                  "сменные операторы и ремонт — существующий штат НПЗ "
                  "(суть fixed $300-400 вместо $400; риск B2 снят)",
            series="6-10 чел на 2-4 скида вкл. сбыт/логистику",
            tier="screening"),
        "rollout": [
            dict(stage="FEED + HAZOP SO3-узла + декларация ПБ / "
                       "Ростехнадзор + трансфертная цена серы с заводом",
                 duration="6-9 мес",
                 note="разрешительная база на ОПО НПЗ — критический путь "
                      "(см. kill so3_regulatory_ban)"),
            dict(stage="лаборатория: G1 (инициатор) + G2 (SO3-источник)",
                 duration="6-12 мес", note="параллельно FEED"),
            dict(stage="пилотный скид 1.1 кт/г + первые контракты "
                       "импортозамещения", duration="12-18 мес"),
            dict(stage="серия до 2-4 скидов (с оглядкой на N2-флот: "
                       "платформа <=10)", duration="12-24 мес"),
            dict(stage="апгрейд конверсии до 20% (после G1 и офф-тейка; "
                       "координация с N2 — рынок насыщается 2-3 скидами)",
                 duration="12 мес"),
        ],
        "trl_gates": [
            "G1 (расчёт, наш квантовый рычаг): дизайн инициатора "
            "радикальной цепи CH4+SO3 — строка INITIATOR $100/т -> цель "
            "<$50/т + путь к конверсии 20%; NEVPT2/DFT-стек; гейт ОБЩИЙ "
            "с N2 — решается один раз на платформу",
            "G2 (лаборатория): пригодность SO3-источника НПЗ — примеси "
            "олеума кислотного хозяйства / серы Клауса (H2S, углеводороды, "
            "зола) для радикального сульфонирования",
            "G3 (пилот): 1000+ ч пробега скида на заводском газе и SO3: "
            "full cost <= $1200/т по счётчикам, качество 99.5%",
            "G4 (рынок РФ): контракты импортозамещения >= 1 кт/г по "
            ">= $2500/т (гальваника/электролиты) до серии",
            "G5 (платформа): решение о разделе ниши L4/N2 — суммарный "
            "флот <= 10 скидов @2.5% (или <= 2-3 @20%), иначе цена MSA "
            "ломается для обеих площадок",
        ],
    }

    out = {
        "model": "Лукойл L4: MSA-скид на НПЗ — сера Клауса + трубный газ + "
                 "обслуживаемая площадка; сбыт — импортозамещение РФ, "
                 "логистика ж/д/авто (не СМП)",
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
        "shared_niche": shared_niche,
        "implementation": implementation,
        "honesty": "скрин +-2-3x. Честные выводы: (1) «сера Клауса дешевле "
                   "рынка» — правильный заголовок, но слабый рычаг: весь "
                   "размах $80-250/т S двигает маржу лишь на ~$57/т MSA "
                   "(0.333 т S/т), а газ трубный vs факельный — на ~$25/т "
                   "(газ ~2.8% full: «на НПЗ газ дорогой» для MSA не "
                   "аргумент); экономику решают ЦЕНА MSA (маржа "
                   "$357..2199/т на $1500..3200), КОНВЕРСИЯ и капекс; "
                   "(2) реальные преимущества L4 против N2 — не сырьё, а "
                   "ОБСЛУЖИВАЕМАЯ ПЛОЩАДКА (fixed $300-400, риск B2 снят "
                   "по построению, staffing инкрементальный), синергия "
                   "SO3/пар (capex $3.5-4.5M) и ж/д-логистика $30-80/т "
                   "против СМП $80-200; у N2 зато сера ~$0 и газ $20-50 "
                   "(~$100/т full в его пользу) — площадки сопоставимы и "
                   "КОНКУРИРУЮТ за одну нишу 50-70 кт/г: совокупный "
                   "потолок ~10 скидов @2.5%, честный раздел — N2 на "
                   "мировой рынок, L4 на внутренний РФ (3-10 кт/г = 1-2 "
                   "скида); (3) премия $3200 — только на срезе "
                   "импортозамещения, оптимистичные 4 скида наполовину "
                   "экспортные по ~$2500-2800; (4) MSA $2500-3200 >> $751 "
                   "(порог 3x майнинга) — тест на specialty пройден; "
                   "(5) главный приз — конверсия 20% (full ~$387/т) — за "
                   "гейтом G1, ОБЩИМ с N2, и делится между площадками: "
                   "рынок насыщается 2-3 скидами @20% на всю платформу.",
    }

    # ------------------------------------------------------------- сводка
    print("=" * 74)
    print("Лукойл L4: MSA-скид на НПЗ — сера Клауса + трубный газ + "
          "обслуживаемая площадка")
    print("   скрин +-2-3x, не банк; логистика ж/д/авто РФ, сбыт — "
          "импортозамещение")
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
          "(цена уже $2500), 8.7 кт/г — гейт G1 общий с N2"
          % (upside["full_usd_t"], upside["landed_full_usd_t"],
             upside["margin_usd_t"]))
    print("\nгаз в себестоимости: %.1f%% full (база) — трубный $180 для "
          "MSA не проблема" % sb["gas_share_of_full_pct"])
    print("совместная ниша с N2: рынок 50-70 кт/г; L4 база 2 + N2 база 3 "
          "= 5 скидов\n  = %.1f-%.1f%% рынка; потолок ПЛАТФОРМЫ ~10 скидов "
          "@2.5%% (16-22%%)"
          % (sb["platform_share_with_n2_pct"][0],
             sb["platform_share_with_n2_pct"][1]))
    print("\nчувствительность (маржа $/т, база):")
    for knob, vals in sens.items():
        if knob == "unit":
            continue
        row = "  ".join("%s:%s" % (k, v) for k, v in vals.items())
        print("  %-44s %s" % (knob, row))
    print("\nkill: цена MSA<$1500 | отказ Ростехнадзора на SO3-узел в "
          "периметре НПЗ |\n      ниша занята (>10 скидов платформы "
          "L4+N2)")
    print("вердикт:", kill["verdict_scenarios"])

    path = os.path.join(DIR, "econ_lukoil_msa_claus_results.json")
    with open(path, "w") as f:
        json.dump(out, f, indent=1, ensure_ascii=False)
    print("\nwrote %s" % os.path.basename(path))


if __name__ == "__main__":
    main()
