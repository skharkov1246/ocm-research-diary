#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
routes/econ_nornickel_ew_oer.py — Норникель N1: OER-анод электроэкстракции
Ni/Cu/Co (Pb-Ca-Sn -> Ti + неблагородное каталитическое покрытие) — экономия
электроэнергии от снижения анодного перенапряжения.

Физика (твёрдый якорь, закон Фарадея): заряд на тонну Q = n*F/M:
  Ni (n=2, M=58.69) = 913 кАч/т; Cu (2, 63.55) = 843; Co (2, 58.93) = 910.
  Каждые -100 мВ на аноде = -91.3 кВт*ч/т Ni, -84.3 Cu, -91.0 Co (при CE=1 —
  консервативно; реальный выход по току 85-95% добавил бы +5-18%).

ВАЖНО (честное разделение): рычаг работает ТОЛЬКО на электроэкстракции
(нерастворимый анод, анодная реакция = OER). На электрорафинировании анод
растворимый (M -> M2+), OER там нет — эти тонны из модели ИСКЛЮЧЕНЫ через
долю ew_share (company_assumption 30-70%, проверить у компании).

Скрин-уровень ±2-3x, НЕ банковская модель. Якоря — из наших же файлов:
  h2_oer_ladder_results.json   eta=0.435 В — OER-лестница [Ni(OH)2] (щёлочь);
                               стартовая точка переноса в сульфатную кислоту
  econ_shock / msa_cost_model  WACC 12% (аннуализация анодного парка)
Данные Норникеля (объёмы, доля EW, цена э/э, площадки) — ТОЛЬКО диапазоны
с tier='company_assumption' и пометкой «проверить у компании».

Запуск: python3 routes/econ_nornickel_ew_oer.py
     -> routes/econ_nornickel_ew_oer_results.json
"""
import json
import math
import os

DIR = os.path.dirname(os.path.abspath(__file__))
F_AH_MOL = 96485.0 / 3600.0          # 26.80 А*ч/моль — постоянная Фарадея

# ------------------------------------------------------ допущения (все крутилки)
ASSUMPTIONS = {
    # --- физика (твёрже любого якоря; tier=anchor_repo с пометкой «физика») ---
    "faraday_ah_mol": dict(value=round(F_AH_MOL, 4), unit="А*ч/моль",
        tier="anchor_repo", source="физика: F = 96485 Кл/моль / 3600 — "
        "не репо-число, твёрже репо-чисел"),
    "metals_electrochem": dict(
        value={"Ni": dict(n=2, M=58.69), "Cu": dict(n=2, M=63.55),
               "Co": dict(n=2, M=58.93)},
        unit="n — электронов; M — г/моль", tier="anchor_repo",
        source="физика: Q=n*F/M -> Ni 913 / Cu 843 / Co 910 кАч/т; "
               "каждые -100 мВ = -91.3 / -84.3 / -91.0 кВт*ч/т"),
    # --- якоря репозитория ---
    "wacc": dict(value=0.12, unit="доля/год", tier="anchor_repo",
        source="econ_shock_candidates / msa_cost_model: 12%; здесь — "
               "CRF-аннуализация анодного парка"),
    "oer_eta_anchor_v": dict(value=0.435, unit="В", tier="anchor_repo",
        source="h2_oer_ladder_results.json: лестница [Ni(OH)2(X)], лимит "
               "dG4 *OOH->O2; щёлочь, минимальный кластер — стартовая точка, "
               "перенос в сульфатную кислоту запланирован (gate G1)"),
    "dft_error_ev": dict(value=[0.2, 0.3], unit="эВ", tier="literature",
        comment="типичная ошибка DFT на *O/*OH/*OOH — ниша нашего квантового "
                "стека: считать точнее там, где DFT врёт"),
    # --- допущения о компании (ПРОВЕРИТЬ У КОМПАНИИ) ---
    "volumes_kt_y": dict(
        value={"Ni": [150, 220], "Cu": [350, 450], "Co": [5, 9]},
        unit="кт/год", tier="company_assumption",
        comment="объёмы рафинирования по металлам группы — проверить у компании"),
    "ew_share": dict(value=[0.30, 0.70], unit="доля тоннажа через "
        "электроэкстракцию", tier="company_assumption",
        comment="остальное — электрорафинирование (растворимый анод, рычаг НЕ "
                "работает). Публично заявлялся перевод Ni Кольской ГМК на "
                "электроэкстракцию — доля по Ni может быть ~0.9-1.0; "
                "проверить у компании"),
    "site_of_metal": dict(value={"Ni": "kola", "Co": "kola", "Cu": "norilsk"},
        unit="-", tier="company_assumption",
        comment="привязка переделов к площадкам (Кольская / Норильск) — "
                "проверить у компании"),
    "el_price_norilsk_usd_mwh": dict(value=[15, 30], unit="$/МВт*ч",
        tier="company_assumption",
        comment="свой газ, изолированная энергосистема — низкая альтернативная "
                "ценность э/э; проверить у компании"),
    "el_price_kola_usd_mwh": dict(value=[30, 50], unit="$/МВт*ч",
        tier="company_assumption", comment="Кольская энергосистема — "
        "проверить у компании"),
    "anodes_per_cell": dict(value=50, unit="шт/ванна", tier="company_assumption",
        comment="диапазон 40-60; только для пересчёта площади в штуки/ванны — "
                "проверить у компании"),
    "anode_geom_m2": dict(value=2.0, unit="м2 эфф./анод",
        tier="company_assumption",
        comment="~1x1 м, обе рабочие стороны — проверить у компании"),
    # --- рынок ---
    "anode_catalytic_usd_m2": dict(value=[1000, 3000], unit="$/м2",
        tier="market",
        comment="MMO/каталитический анод на Ti-основе, рыночный диапазон; "
                "наша научная цель — сбить к breakeven (см. результаты)"),
    "anode_pb_usd_m2": dict(value=300.0, unit="$/м2", tier="market",
        comment="Pb-Ca-Sn бейзлайн, диапазон 200-400"),
    # --- литература ---
    "dv_scenarios_mv": dict(value=[100, 200, 300], unit="мВ",
        tier="literature",
        comment="сценарии снижения анодного перенапряжения; справка: переход "
                "Pb -> MMO/каталитический анод в индустрии даёт сотни мВ "
                "(порядка 200-500)"),
    "anode_life_catalytic_y": dict(value=[3, 7], unit="лет", tier="literature",
        comment="заявления вендоров MMO; kill при <2 лет"),
    "anode_life_pb_y": dict(value=7.0, unit="лет", tier="literature",
        comment="Pb-Ca-Sn бейзлайн, диапазон 5-9"),
    "current_density_a_m2": dict(value={"Ni": 200.0, "Cu": 300.0, "Co": 200.0},
        unit="А/м2", tier="literature",
        comment="типовые плотности тока EW (Ni/Co 180-250, Cu 250-350); "
                "для оценки площади анодного парка"),
    "current_efficiency": dict(value=0.90, unit="доля", tier="literature",
        comment="выход по току 85-95%; учтён в токе/площади парка; в формуле "
                "экономии НЕ учтён (CE=1 — консервативно, якорь ТЗ)"),
    "hours_per_year": dict(value=8000.0, unit="ч/год", tier="literature",
        comment="коэффициент использования цеха ~0.91"),
}


def V(k):
    return ASSUMPTIONS[k]["value"]


def crf(life_y):
    """Фактор возврата капитала: аннуализация по WACC 12%."""
    w = V("wacc")
    return w / (1.0 - (1.0 + w) ** (-life_y))


def pick(rng, w):
    """Точка в диапазоне [lo, hi]: w=0 низ, 0.5 середина, 1 верх."""
    return rng[0] + (rng[1] - rng[0]) * w


# --------------------------------------------- сценарии (pess / base / opt)
# Для МОДЕЛИ ЭКОНОМИИ пессимизм = малый сдвиг мВ, малая доля EW, дешёвая
# энергия, дорогой недолговечный анод. База ΔV=-200 мВ (середина справки
# Pb->MMO), середины company-диапазонов.
SCENARIOS = {
    "pessimistic": dict(dv_mv=100.0, ew_share=0.30, vol_w=0.0, price_w=0.0,
        anode_usd_m2=3000.0, anode_life_y=3.0,
        note="ΔV=-100 мВ, доля EW 30%, низ объёмов, э/э по низу, анод дорогой "
             "и живёт 3 года"),
    "base": dict(dv_mv=200.0, ew_share=0.50, vol_w=0.5, price_w=0.5,
        anode_usd_m2=2000.0, anode_life_y=5.0,
        note="ΔV=-200 мВ, доля EW 50%, середины company-диапазонов"),
    "optimistic": dict(dv_mv=300.0, ew_share=0.70, vol_w=1.0, price_w=1.0,
        anode_usd_m2=1000.0, anode_life_y=7.0,
        note="ΔV=-300 мВ, доля EW 70%, верх объёмов и цен э/э, анод по низу "
             "рынка, 7 лет жизни"),
}


def q_kah_t(metal):
    """Заряд на тонну металла, кАч/т (физика: Q = n*F/M)."""
    m = V("metals_electrochem")[metal]
    return m["n"] * F_AH_MOL / m["M"] * 1000.0


def metal_block(metal, dv_mv, ew_share, vol_w, price_w):
    """Один металл: тоннаж через EW, ГВт*ч и $M экономии, площадь анодов."""
    q = q_kah_t(metal)
    vol_kt = pick(V("volumes_kt_y")[metal], vol_w)
    t_ew = vol_kt * 1000.0 * ew_share
    saving_kwh_t = q * dv_mv / 1000.0           # кАч/т * В = кВт*ч/т
    gwh = t_ew * saving_kwh_t / 1e6
    site = V("site_of_metal")[metal]
    price = pick(V("el_price_%s_usd_mwh" % site), price_w)
    musd = gwh * price / 1000.0
    amps = t_ew * q * 1000.0 / (V("current_efficiency") *
                                V("hours_per_year"))
    area = amps / V("current_density_a_m2")[metal]
    return dict(metal=metal, site=site, vol_kt_y=round(vol_kt, 0),
                t_ew_y=round(t_ew, 0), q_kah_t=round(q, 0),
                saving_kwh_t=round(saving_kwh_t, 1),
                el_price_usd_mwh=round(price, 1),
                saving_gwh_y=round(gwh, 1), saving_musd_y=round(musd, 2),
                anode_area_m2=round(area, -2))


def scenario(dv_mv, ew_share, vol_w, price_w, anode_usd_m2, anode_life_y,
             **_):
    """Полный сценарий: три металла + площадки + анодный парк + payback."""
    metals = [metal_block(m, dv_mv, ew_share, vol_w, price_w)
              for m in ("Ni", "Cu", "Co")]
    gwh = sum(m["saving_gwh_y"] for m in metals)
    musd = sum(m["saving_musd_y"] for m in metals)
    area = sum(m["anode_area_m2"] for m in metals)
    sites = {}
    for s in ("kola", "norilsk"):
        ms = [m for m in metals if m["site"] == s]
        sites[s] = dict(
            metals=[m["metal"] for m in ms],
            saving_gwh_y=round(sum(m["saving_gwh_y"] for m in ms), 1),
            saving_musd_y=round(sum(m["saving_musd_y"] for m in ms), 2))
        sites[s]["kill_below_2musd"] = sites[s]["saving_musd_y"] < 2.0
    # анодный парк: полный своп по рыночной цене vs Pb-бейзлайн (CRF 12%)
    capex = area * anode_usd_m2 / 1e6
    ann_cat = capex * crf(anode_life_y)
    ann_pb = area * V("anode_pb_usd_m2") * crf(V("anode_life_pb_y")) / 1e6
    net = musd - (ann_cat - ann_pb)
    pb_gross = capex / musd if musd > 0 else None
    # breakeven: цена покрытия, при которой аннуализированная надбавка к Pb
    # равна экономии энергии — ЦЕЛЬ для нашего неблагородного катализатора
    be = ((musd * 1e6 / area + V("anode_pb_usd_m2") *
           crf(V("anode_life_pb_y"))) / crf(anode_life_y)) if area else None
    anodes = area / V("anode_geom_m2")
    return dict(
        per_metal=metals, sites=sites,
        saving_gwh_y=round(gwh, 1), saving_musd_y=round(musd, 2),
        anode_area_total_m2=round(area, -3),
        anodes_pcs=round(anodes, -2),
        cells=round(anodes / V("anodes_per_cell"), -1),
        swap_capex_musd=round(capex, 0),
        annualized_catalytic_musd_y=round(ann_cat, 1),
        annualized_pb_baseline_musd_y=round(ann_pb, 1),
        net_annual_musd_y=round(net, 1),
        payback_gross_y=round(pb_gross, 0) if pb_gross else None,
        breakeven_anode_usd_m2=round(be, -1) if be else None)


def main():
    # -------- сверка с якорями (обязана сходиться, иначе модель врёт)
    lad = json.load(open(os.path.join(DIR, "h2_oer_ladder_results.json")))
    anchors = {
        "q_ni_kah_t": dict(model=round(q_kah_t("Ni"), 0), anchor=913.0,
                           source="физика Q=nF/M (ТЗ)"),
        "q_cu_kah_t": dict(model=round(q_kah_t("Cu"), 0), anchor=843.0,
                           source="физика"),
        "q_co_kah_t": dict(model=round(q_kah_t("Co"), 0), anchor=910.0,
                           source="физика"),
        "kwh_per_100mv_ni": dict(model=round(q_kah_t("Ni") * 0.1, 1),
                                 anchor=91.3, source="физика"),
        "kwh_per_100mv_cu": dict(model=round(q_kah_t("Cu") * 0.1, 1),
                                 anchor=84.3, source="физика"),
        "kwh_per_100mv_co": dict(model=round(q_kah_t("Co") * 0.1, 1),
                                 anchor=91.0, source="физика"),
        "oer_eta_v": dict(model=lad["ladder"]["eta_V"], anchor=0.435,
                          source="h2_oer_ladder_results.json"),
    }
    for a in anchors.values():
        a["ok"] = abs(a["model"] - a["anchor"]) <= max(0.005 * a["anchor"],
                                                       0.05)

    scen_out = {name: dict(params=dict(p), **scenario(**p))
                for name, p in SCENARIOS.items()}
    base = SCENARIOS["base"]
    sb = scen_out["base"]

    # -------- чувствительность (вокруг базы): 6 ручек
    def tot(**over):
        return scenario(**dict(base, **over))

    sens = {
        "dv_mv (сценарии -100/-200/-300)": {
            "-%d мВ" % dv: tot(dv_mv=dv)["saving_musd_y"]
            for dv in (100, 200, 300)},
        "ew_share (доля электроэкстракции)": {
            "%.0f%%" % (s * 100): tot(ew_share=s)["saving_musd_y"]
            for s in (0.2, 0.3, 0.5, 0.7, 0.9)},
        "el_price (+-30% к серединам)": {
            "%+.0f%%" % (m * 100):
                tot(price_w=max(0.0, min(1.0, 0.5 + m * 1.5)))
                ["saving_musd_y"]
            for m in (-0.3, 0.0, 0.3)},
        "volumes (низ/середина/верх company-диапазона)": {
            w_name: tot(vol_w=w)["saving_musd_y"]
            for w_name, w in (("низ", 0.0), ("середина", 0.5), ("верх", 1.0))},
        "anode_usd_m2 -> [net $M/г, валовый payback лет]": {
            "$%d/м2" % p: [tot(anode_usd_m2=p)["net_annual_musd_y"],
                           tot(anode_usd_m2=p)["payback_gross_y"]]
            for p in (300, 1000, 2000, 3000)},
        "anode_life_y -> [net $M/г, breakeven $/м2]": {
            "%d лет" % ln: [tot(anode_life_y=ln)["net_annual_musd_y"],
                            tot(anode_life_y=ln)["breakeven_anode_usd_m2"]]
            for ln in (2, 3, 5, 7)},
        "unit": "экономия $M/год на весь дивизион (если не сказано иное); "
                "price_w для +-30% зажат в company-диапазон",
    }

    # -------- kill-критерии
    k_share = tot(ew_share=0.20)
    k_life = tot(anode_life_y=2.0)
    kill = {
        "ew_share_below": dict(kill_below=0.20,
            computed_context="при доле EW 20%% экономия $%.2fM/год на весь "
                "дивизион — тоннажный рычаг мал, вариант закрыть"
                % k_share["saving_musd_y"]),
        "catalyst_life_below_y": dict(kill_below=2.0,
            computed_context="при жизни покрытия 2 года breakeven-цена падает "
                "к $%.0f/м2 — ниже реалистичной себестоимости покрытия; "
                "кампании перезамены съедают эффект"
                % k_life["breakeven_anode_usd_m2"]),
        "saving_per_site_below_musd_y": dict(kill_below=2.0,
            computed_context="база: Кольская $%.2fM/г, Норильск $%.2fM/г — "
                "ОБЕ площадки ниже порога $2M; критерий проходит только "
                "оптимистичный угол (ΔV=-300, доля 70%%, э/э по верху)"
                % (sb["sites"]["kola"]["saving_musd_y"],
                   sb["sites"]["norilsk"]["saving_musd_y"])),
        "anode_above_breakeven_usd_m2": dict(
            kill_above=sb["breakeven_anode_usd_m2"],
            computed_context="производный критерий ИЗ модели: при рыночных "
                "$1-3k/м2 энергорычаг не отбивает парк (валовый payback "
                "%.0f лет в базе); покрытие обязано стоить ~как Pb"
                % sb["payback_gross_y"]),
        "verdict_scenarios": {},
    }
    for name, s in scen_out.items():
        p = SCENARIOS[name]
        alive = (p["ew_share"] >= 0.20 and p["anode_life_y"] >= 2.0 and
                 all(not st["kill_below_2musd"] for st in s["sites"].values()))
        flag = ("жив по kill-критериям, НО анод по рыночной цене не "
                "окупается — нужен таргет <= breakeven $%.0f/м2"
                % s["breakeven_anode_usd_m2"]) if alive else "МЁРТВ"
        kill["verdict_scenarios"][name] = flag

    # -------- доп. эффекты (качественно, НЕ монетизированы)
    side_benefits = [
        "нет Pb в электролите и катоде — выше качество катодного металла, "
        "меньше нагрузки на очистку растворов",
        "нет шлама PbO2/PbSO4 — меньше чисток ванн и обращения с Pb-отходами",
        "стабильная геометрия анода (Ti не ползёт, Pb коробится) — ровнее "
        "распределение тока, меньше коротких замыканий, потенциально +1-3% "
        "к выходу по току (НЕ монетизировано)",
        "меньше кислотного тумана над ваннами — охрана труда, вентиляция",
        "для Cu EW — отказ от дозирования CoSO4 в электролит (его льют для "
        "защиты Pb-анода) — прямая экономия реагента (НЕ монетизирована)",
    ]

    # -------- наш вклад
    our_lever = dict(
        what="дизайн НЕБЛАГОРОДНОГО OER-катализатора для сульфатной кислой "
             "среды — там, где DFT врёт на +-0.2-0.3 эВ по *O/*OH/*OOH",
        anchor="eta=0.435 В — лестница [Ni(OH)2] h2-трека "
               "(h2_oer_ladder_results.json), щёлочь, минимальный кластер",
        planned="перенос лестницы в сульфатную кислую среду — расчёт уже "
                "запланирован (gate G1)",
        econ_target="ИЗ ЭТОЙ МОДЕЛИ: покрытие без Ir/Ru по цене <= breakeven "
                    "$%.0f-%.0f/м2 (base-opt) при ΔV=-200..-300 мВ и жизни "
                    ">=5 лет — иначе экономика анодного парка не сходится"
                    % (sb["breakeven_anode_usd_m2"],
                       scen_out["optimistic"]["breakeven_anode_usd_m2"]))

    # -------- механическая реализация
    implementation = {
        "block_flow": [
            "существующий гидрометаллургический передел: сульфатный раствор "
            "Ni/Cu/Co после выщелачивания и очистки (БЕЗ изменений)",
            "цех электроэкстракции: ванны, циркуляция электролита, катодные "
            "основы (БЕЗ изменений)",
            "ЗАМЕНА: анодный парк Pb-Ca-Sn -> Ti-основа с неблагородным "
            "OER-покрытием, drop-in в штатный подвес",
            "выпрямители/шины: без замены — снижение напряжения ячейки на "
            "0.1-0.3 В напрямую режет кВт*ч/т",
            "мониторинг: анодные потенциалы, Pb и компоненты покрытия в "
            "катоде/электролите, тренд деградации",
            "recoat-цикл: изношенное покрытие -> перенанесение (мастерская "
            "на площадке или у поставщика), Ti-основа многоразовая",
        ],
        "equipment": [
            dict(item="анодный парк, полный своп (~%.0f тыс. м2, ~%.0f ванн)"
                 % (sb["anode_area_total_m2"] / 1000.0, sb["cells"]),
                 cost_musd=[round(sb["anode_area_total_m2"] *
                                  sb["breakeven_anode_usd_m2"] / 1e6, 0),
                            sb["swap_capex_musd"]],
                 note="низ — по breakeven-цене (наша цель), верх — рынок "
                      "$2k/м2; при рыночной цене НЕ окупается"),
            dict(item="пилотная секция: 1 ванна, ~50 анодов (~100 м2)",
                 cost_musd=[0.1, 0.3]),
            dict(item="стенд ускоренной деградации + электрохим. лаборатория",
                 cost_musd=[0.3, 0.6]),
            dict(item="recoat-мастерская на площадке (опция)",
                 cost_musd=[1.0, 3.0]),
            dict(item="мониторинг (анодные потенциалы, ICP-контроль Pb)",
                 cost_musd=[0.5, 1.0]),
        ],
        "footprint_m2": [200, 400],
        "footprint_note": "прирост площади ~0: drop-in в действующий цех; "
                          "200-400 м2 — recoat-мастерская и лаборатория",
        "utilities": {
            "electricity_kw": dict(
                value=[-15000, -1500],
                note="ОТРИЦАТЕЛЬНО — сам эффект: -13..-124 ГВт*ч/год по "
                     "сценариям = -1.5..-15 МВт средней нагрузки; "
                     "плюс мастерская +50-150 кВт"),
            "cooling_water_m3_h": [1, 3],
            "nitrogen_nm3_h": 0,
            "note": "новых утилит по сути нет — проект МИНУСУЕТ энергию",
        },
        "staffing": dict(
            pilot="3-4 чел: инженер-электрохимик, технолог, операторы цеха "
                  "по совместительству",
            series="кампании замены силами цеха + 1-2 инженера сопровождения",
            tier="screening"),
        "rollout": [
            dict(stage="FEED + кислотная OER-лестница (расчёт, gate G1)",
                 duration="3-6 мес"),
            dict(stage="лабораторный образец покрытия + 1000 ч стенд (G2)",
                 duration="6-9 мес"),
            dict(stage="пилот: 1 ванна в действующем цехе (G3)",
                 duration="9-12 мес", note="включая 6 мес пробега"),
            dict(stage="секция ~10% ванн одной площадки",
                 duration="12 мес"),
            dict(stage="полный парк в ШТАТНОМ цикле замены Pb-анодов",
                 duration="24-48 мес",
                 note="инкрементальный капекс -> 0, если менять по износу Pb"),
        ],
        "trl_gates": [
            "G1 (расчёт, до FEED): перенос OER-лестницы (якорь eta=0.435 В, "
            "щёлочь) в сульфатную кислую среду: eta <= 0.30 В на "
            "неблагородном мотиве + стабильность по Пурбе — расчёт уже "
            "запланирован",
            "G2 (лаборатория, до пилота): 1000 ч в синтетическом электролите "
            "(H2SO4 + сульфаты Ni/Cu/Co, следы Cl-): -150 мВ vs Pb "
            "удержано, дрейф <10 мВ/1000 ч; себестоимость покрытия в серии "
            "<= breakeven $%.0f/м2" % sb["breakeven_anode_usd_m2"],
            "G3 (пилот, до секции): 6 мес в действующей ванне: кВт*ч/т по "
            "счётчику, Pb в катоде снижен, электролит не отравлен "
            "компонентами покрытия",
            "G4 (серия): подтверждённая жизнь >=5 лет (экстраполяция "
            "деградации) и доля EW-тоннажа >=20% подтверждена компанией",
        ],
    }

    out = {
        "model": "Норникель N1: OER-анод электроэкстракции Ni/Cu/Co — "
                 "экономия электроэнергии",
        "tier": "screening_pm_2_3x",
        "not_bankable": True,
        "tier_legend": {
            "anchor_repo": "число из наших файлов или физика (закон Фарадея "
                           "— помечено в source)",
            "market": "рыночный диапазон, порядок величины",
            "company_assumption": "допущение о компании — проверить у неё",
            "literature": "литература/индустриальная практика, не наш расчёт"},
        "anchors_check": anchors,
        "assumptions": ASSUMPTIONS,
        "scenarios": scen_out,
        "sensitivity": sens,
        "kill_criteria": kill,
        "side_benefits_qualitative": side_benefits,
        "our_lever": our_lever,
        "implementation": implementation,
        "honesty": "скрин +-2-3x. Главный честный вывод: на энергии ОДНОЙ "
                   "кейс тонкий — база ~$1.5M/год на весь дивизион при "
                   "captive-энергии $22-40/МВт*ч, обе площадки НИЖЕ "
                   "kill-порога $2M/год; рыночный MMO-анод $1-3k/м2 не "
                   "окупается (валовый payback ~190 лет). Кейс живёт только "
                   "как (1) ДЕШЁВОЕ неблагородное покрытие <= breakeven "
                   "$280-390/м2 (base-opt) — это и есть цель катализатора — "
                   "ставящееся в штатный цикл замены Pb (инкрементальный "
                   "капекс -> 0), (2) плюс немонетизированные side benefits "
                   "(Pb в катоде, стабильность тока, туман). Экономия "
                   "консервативна: CE=1 (реальные 85-95% дали бы +5-18%). "
                   "Доля EW против рафинирования — ключевое непроверенное "
                   "допущение; для рафинирования рычаг НЕ работает.",
    }

    # ------------------------------------------------------------- сводка
    print("=" * 74)
    print("Норникель N1: OER-анод электроэкстракции Ni/Cu/Co"
          "   (скрин +-2-3x, не банк)")
    print("=" * 74)
    print("сверка якорей:", "OK" if all(a["ok"] for a in anchors.values())
          else "РАСХОЖДЕНИЕ! " + str(anchors))
    print("физика: Q(кАч/т) Ni %.0f / Cu %.0f / Co %.0f;  -100 мВ = "
          "-%.1f / -%.1f / -%.1f кВт*ч/т"
          % (q_kah_t("Ni"), q_kah_t("Cu"), q_kah_t("Co"),
             q_kah_t("Ni") * 0.1, q_kah_t("Cu") * 0.1, q_kah_t("Co") * 0.1))
    print("\n%-13s%7s%9s%10s%9s%10s%12s" % ("сценарий", "ΔV,мВ", "ГВт*ч/г",
          "$M/г", "тыс.м2", "capex $M", "breakeven"))
    for name, s in scen_out.items():
        p = SCENARIOS[name]
        print("%-13s%7.0f%9.1f%10.2f%9.0f%10.0f%9.0f$/м2"
              % (name, p["dv_mv"], s["saving_gwh_y"], s["saving_musd_y"],
                 s["anode_area_total_m2"] / 1000.0, s["swap_capex_musd"],
                 s["breakeven_anode_usd_m2"]))
    print("\nбаза по площадкам: Кольская (Ni+Co) $%.2fM/г | Норильск (Cu) "
          "$%.2fM/г — kill-порог $2M/площадку НЕ пройден"
          % (sb["sites"]["kola"]["saving_musd_y"],
             sb["sites"]["norilsk"]["saving_musd_y"]))
    print("анодный парк (база): ~%.0f тыс. м2 ~ %.0f ванн; своп по рынку "
          "$%.0fM -> валовый payback %.0f лет (МЁРТВ по капексу)"
          % (sb["anode_area_total_m2"] / 1000.0, sb["cells"],
             sb["swap_capex_musd"], sb["payback_gross_y"]))
    print("\nчувствительность (экономия $M/год, вся компания):")
    for knob, vals in sens.items():
        if knob == "unit":
            continue
        row = "  ".join("%s:%s" % (k, v) for k, v in vals.items())
        print("  %-46s %s" % (knob, row))
    print("\nkill: доля EW<20%% | жизнь покрытия<2 лет | <$2M/г на площадку "
          "| анод дороже breakeven $%.0f/м2"
          % sb["breakeven_anode_usd_m2"])
    print("вердикт:", kill["verdict_scenarios"])
    print("\nнаш рычаг:", our_lever["econ_target"])

    path = os.path.join(DIR, "econ_nornickel_ew_oer_results.json")
    with open(path, "w") as f:
        json.dump(out, f, indent=1, ensure_ascii=False)
    print("\nwrote %s" % os.path.basename(path))


if __name__ == "__main__":
    main()
