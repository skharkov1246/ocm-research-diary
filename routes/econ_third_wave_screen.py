#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
routes/econ_third_wave_screen.py — ТРЕТЬЯ ВОЛНА идей: один компактный скрин
10 идей (Лукойл L11-L15, Норникель N13-N17) в одном скрипте.

Уровень: screening_pm_3_5x — ГРУБЕЕ второй волны (±2-3x): по 5-8 параметров
на идею, unit economics одной установки/линии, три угла (pess/base/opt),
вердикт alive/borderline/dead + ГЛАВНЫЙ гейт и kill-критерий на идею.
Не банк. Все company_assumption — проверить у компании.

Идеи:
  L11 DMDS       — диметилдисульфид из CH3SH (кислый газ): сульфидирование
                   CoMo/NiMo (связка L6) + агрохимия
  L12 ODS        — окислительная десульфуризация судового топлива
                   H2O2 + Mo-POM (стык N7 и нашей Mo-химии)
  L13 CO2-EOR    — свой CO2 (SMR / будущий L3) в зрелые пласты
  L14 Пек        — мезофазный пек из гудрона (гейт — π-радикальная
                   конденсация, наша NEVPT2-химия)
  L15 HPPO       — пропиленоксид на Ставролене из C3H6 + H2O2 (стык N7)
  N13 Ni-электр. — активированные Ni-электроды для щелочных электролизёров
                   (ПРЯМОЙ якорь: наша OER-лестница eta=0.435 В)
  N14 Cu2O       — антифоулинг для своего флота СМП (честно микро)
  N15 Селитра    — мини-Haber-Bosch на Таймыре -> АС/эмульсионные ВВ для
                   своих рудников (классика, НЕ наша электрохимия — честно)
  N16 Ni-63      — бета-вольтаика с партнёром (Росатом): concept_no_economics
  N17 FeCl3      — коагулянт из пирротиновых хвостов (T3) + хлора Кольской

Кросс-чек якорей репозитория: h2_oer_ladder_results.json (eta 0.435 В),
calc/fes_cogef_results.json (F_max 1.97 нН — механовскрытие Fe-S хвостов),
msa_cost_model.json (формула капитального платежа), econ_vs_mining_results
($751/т — порог specialty против майнинга 3x).

Запуск: python3 routes/econ_third_wave_screen.py
     -> routes/econ_third_wave_screen_results.json  (детерминированно)
"""
import json
import os

DIR = os.path.dirname(os.path.abspath(__file__))

WACC = 0.12  # anchor_repo: msa_cost_model / skid_twin — 12% (аморт.+WACC)
MINING_X3_USD_T = 751.0  # anchor_repo: econ_vs_mining_results


def A(value, unit, tier, note=None):
    """Параметр с тиром достоверности (5-8 штук на идею, скрин ±3-5x)."""
    d = dict(value=value, unit=unit, tier=tier)
    if note:
        d["note"] = note
    return d


def econ(price, var, fixed, units_y, capex_musd, delivery=0.0):
    """Unit economics одной установки/линии: капплатёж = capex*WACC/выпуск
    (та же формула, что в msa_cost_model — сверено в anchors_check)."""
    capital = capex_musd * 1e6 * WACC / units_y
    full = var + fixed + capital
    margin = price - full - delivery
    ebitda = margin * units_y / 1e6
    return dict(
        full_usd_unit=round(full, 1), margin_usd_unit=round(margin, 1),
        units_y=round(units_y), ebitda_musd_y=round(ebitda, 2),
        capex_musd=capex_musd,
        payback_y=round(capex_musd / ebitda, 1) if ebitda > 0 else None)


# =========================================================== ЛУКОЙЛ L11-L15

def l11_dmds():
    a = {
        "dmds_price_usd_t": A([1500, 2500], "$/т", "market",
            "DMDS: сульфидирующий агент CoMo/NiMo (связка L6 — свой же "
            "гидроочистной катализатор) + агрохимия (фумигант)"),
        "line_kt_y": A([5, 10], "кт/год", "market",
            "мировой рынок DMDS ~30-50 кт/г — линия 5-10 кт заметна, "
            "но не ломает цену"),
        "meoh_usd_t": A([250, 350], "$/т", "market"),
        "meoh_t_per_t": A(0.75, "т MeOH/т DMDS", "literature",
            "стехиометрия 2 CH3SH на (CH3S)2 ~0.68 + потери"),
        "h2s_own_usd_t": A([0, 50], "$/т H2S", "company_assumption",
            "кислый газ свой (аминовая очистка НПЗ) — трансфертная цена "
            "~0, возможен даже кредит утилизации; проверить у компании"),
        "h2s_t_per_t": A(0.4, "т H2S/т DMDS", "literature"),
        "opex_fixed_usd_t": A([200, 400], "$/т", "company_assumption",
            "обслуживаемая площадка НПЗ (как L4)"),
        "capex_musd": A([15, 30], "$M", "market",
            "линия меркаптан + окислительная конденсация, E1 ±3-5x"),
    }

    def run(price, meoh, h2s, opx, kt, capex):
        return econ(price, 0.75 * meoh + 0.4 * h2s, opx, kt * 1e3, capex)

    sc = dict(
        pess=run(1500, 350, 50, 400, 5, 30),
        base=run(2000, 300, 25, 300, 7.5, 22),
        opt=run(2500, 250, 0, 200, 10, 15))
    return dict(
        owner="Лукойл", name="L11 DMDS: диметилдисульфид из кислого газа",
        link="связка L6 (CoMo/NiMo гидроочистка): свой же спрос на "
             "сульфидирование — гарантированный якорный offtake",
        assumptions=a, unit="$/т DMDS", scenarios=sc,
        anchors_used=dict(mining_x3="цена $1500-2500 >> $751 — specialty, "
                          "тест econ_vs_mining пройден"),
        gate="ГЛАВНЫЙ ГЕЙТ: качество сульфидирующего агента (примеси в "
             "CH3SH из кислого газа) + подтверждённый offtake L6/агрохимия",
        kill="цена DMDS < $1000/т (китайское перепредложение) — маржа "
             "базы падает к нулю; или запрет меркаптанового узла (ОПО)",
        verdict="alive",
        verdict_note="база: маржа $%d/т, EBITDA $%.1fM/г, окуп. %.1f г; "
                     "даже pess-угол ($1500, всё дорого) держит маржу "
                     "$%d/т > 0 — редкость в этой волне"
                     % (sc["base"]["margin_usd_unit"],
                        sc["base"]["ebitda_musd_y"],
                        sc["base"]["payback_y"],
                        sc["pess"]["margin_usd_unit"]))


def l12_ods():
    a = {
        "hsfo_vlsfo_spread_usd_t": A([80, 200], "$/т", "market",
            "экономика = спред: окислили HSFO -> продали как compliant; "
            "спред волатилен (IMO-2020 пик ~$300, спады к $60-80)"),
        "unit_kt_y": A([100, 500], "кт/год", "market"),
        "h2o2_kg_t": A([3, 6], "кг/т топлива", "literature",
            "стехиометрия S->сульфон ~2 моль H2O2/моль S + избыток"),
        "h2o2_usd_t": A([600, 900], "$/т", "market",
            "покупной 100%-эквивалент; стык N7 (прямой синтез) — ниже"),
        "cat_makeup_usd_t": A([2, 5], "$/т топлива", "company_assumption",
            "догрузка Mo-полиоксометаллата — НАША Mo-химия (спин/пероксо)"),
        "yield_loss_usd_t": A([5, 20], "$/т", "literature",
            "экстракция сульфонов уносит 1-3% топлива"),
        "opex_usd_t": A([15, 30], "$/т", "company_assumption"),
        "capex_musd": A([30, 80], "$M", "market", "на 100-500 кт/год"),
    }

    def run(spread, h2o2kg, h2o2, cat, loss, opx, kt, capex):
        return econ(spread, h2o2kg * h2o2 / 1000 + cat + loss, opx,
                    kt * 1e3, capex)

    sc = dict(
        pess=run(80, 6, 900, 5, 20, 30, 100, 30),
        base=run(140, 4.5, 750, 3, 10, 22, 300, 55),
        opt=run(200, 3, 600, 2, 5, 15, 500, 80))
    return dict(
        owner="Лукойл", name="L12 ODS: окислительная десульфуризация "
                             "судового топлива (H2O2 + Mo-POM)",
        link="стык N7 (дешёвый H2O2) + наша Mo-пероксо-химия; убирает S "
             "БЕЗ водорода — против гидроочистки это и есть смысл",
        assumptions=a, unit="$/т топлива (доход = спред HSFO/VLSFO)",
        scenarios=sc,
        gate="ГЛАВНЫЙ ГЕЙТ: устойчивый спред >= $100/т + ПРИЁМКА "
             "окисленного топлива (спецификация ISO 8217 / стабильность "
             "после экстракции сульфонов) — без неё продукта нет",
        kill="спред < $60/т дольше года (маржа базы уходит < 0) или "
             "отказ бункеровщиков от окисленного продукта",
        verdict="borderline",
        verdict_note="база: маржа $%.0f/т при спреде $140, EBITDA $%.1fM/г,"
                     " окуп. %.1f г — красиво, но ВЕСЬ доход = волатильный "
                     "спред: pess ($80) уже -$%.0f/т; это опцион на спред, "
                     "не бизнес с себестоимостной защитой"
                     % (sc["base"]["margin_usd_unit"],
                        sc["base"]["ebitda_musd_y"],
                        sc["base"]["payback_y"],
                        -sc["pess"]["margin_usd_unit"]))


def l13_co2_eor():
    a = {
        "bbl_per_t_co2": A([2, 4], "барр./т CO2", "literature",
            "инкрементальная нефть на тонну закачанного CO2 (miscible "
            "flood, мировая практика); на конкретном пласте — MMP-тест"),
        "urals_usd_bbl": A([55, 70], "$/барр.", "market"),
        "inject_usd_t": A([20, 40], "$/т CO2", "literature",
            "компрессия + скважины + рецикл CO2"),
        "co2_supply_usd_t": A([25, 50], "$/т CO2", "company_assumption",
            "захват+осушка+труба от своего SMR (или будущего L3); "
            "проверить у компании"),
        "lifting_usd_bbl": A([15, 25], "$/барр.", "literature",
            "подъём/подготовка инкрементальной нефти"),
        "tax_take": A([0.15, 0.45], "доля выручки", "company_assumption",
            "НДПИ/экспортная нагрузка; 0.15 = адресная льгота на "
            "ТрИЗ/EOR — ключевая ручка, проверить у компании"),
        "scale_mt_y": A([0.5, 1.0], "Мт CO2/год", "market"),
        "capex_musd": A([100, 200], "$M", "market",
            "труба CO2 + компрессия + переобвязка куста, E1 ±3-5x"),
    }

    def run(bbl, oil, inj, sup, lift, tax, mt, capex):
        return econ(bbl * oil * (1 - tax), inj + sup + bbl * lift, 0.0,
                    mt * 1e6, capex)

    sc = dict(
        pess=run(2, 55, 40, 50, 25, 0.45, 0.5, 100),
        base=run(3, 62, 30, 35, 20, 0.30, 0.75, 150),
        opt=run(4, 70, 20, 25, 15, 0.15, 1.0, 200))
    return dict(
        owner="Лукойл", name="L13 CO2-EOR: свой CO2 в зрелые пласты",
        link="утилизация CO2 от SMR (и будущего L3): вместо платы за "
             "выброс — баррели; углеродная отчётность бонусом",
        assumptions=a, unit="$/т CO2 (доход = инкрементальная нефть)",
        scenarios=sc,
        gate="ГЛАВНЫЙ ГЕЙТ: налоговый режим (льгота НДПИ на EOR/ТрИЗ) + "
             "подтверждение 2-4 барр./т CO2 на КОНКРЕТНОМ пласте (MMP)",
        kill="< 1.5 барр./т CO2 на пилотном участке или полная налоговая "
             "ставка без льготы — маржа < 0 при любой цене Urals $55-70",
        verdict="borderline",
        verdict_note="база: маржа $%.0f/т CO2 — В МИНУСЕ после налога "
                     "30%% и капплатежа; вся экономика живёт в налоговой "
                     "ручке: льгота 15%% даёт $%.0f/т ($%.0fM/г на 1 Мт, "
                     "окуп. %.1f г), полная ставка убивает"
                     % (sc["base"]["margin_usd_unit"],
                        sc["opt"]["margin_usd_unit"],
                        sc["opt"]["ebitda_musd_y"],
                        sc["opt"]["payback_y"]))


def l14_pitch():
    a = {
        "gudron_usd_t": A([150, 300], "$/т", "market",
            "остаток вакуумной перегонки — свой, трансферт ~топочный мазут"),
        "pitch_price_usd_t": A([800, 3000], "$/т", "market",
            "изотропный пек (связующее анодов/электродов) $800-1200; "
            "МЕЗОФАЗНЫЙ (прекурсор углеволокна/игольчатого кокса) "
            "$2000-3000 — премия и есть приз"),
        "yield_t_t": A([0.4, 0.6], "т пека/т гудрона", "literature",
            "термополиконденсация 400-440°C"),
        "distillate_credit_usd_t": A([100, 200], "$/т пека",
            "company_assumption", "отгон (крекинг-дистилляты) в топливный "
            "пул завода"),
        "opex_usd_t": A([150, 350], "$/т пека", "company_assumption"),
        "line_kt_y": A([20, 50], "кт пека/год", "market"),
        "capex_musd": A([40, 80], "$M", "market"),
    }

    def run(price, gudron, y, cred, opx, kt, capex):
        return econ(price, gudron / y - cred, opx, kt * 1e3, capex)

    sc = dict(
        pess=run(800, 300, 0.4, 100, 350, 20, 40),
        base=run(1500, 225, 0.5, 150, 250, 35, 60),
        opt=run(3000, 150, 0.6, 200, 150, 50, 70))
    return dict(
        owner="Лукойл", name="L14 Мезофазный пек из гудрона",
        link="гейт — контроль радикальной конденсации ПАУ: наша "
             "NEVPT2 π-радикальная химия (какие радикалы растить, какие "
             "глушить — от этого зависит мезофаза vs изотроп vs кокс)",
        assumptions=a, unit="$/т пека", scenarios=sc,
        anchors_used=dict(mining_x3="мезофазная цена $1500-3000 >> $751 "
                          "— specialty; изотропная $800-1200 — на грани"),
        gate="ГЛАВНЫЙ ГЕЙТ: воспроизводимое МЕЗОФАЗНОЕ качество "
             "(содержание мезофазы, QI, точка размягчения) = контроль "
             "π-радикальной конденсации; изотроп премии не несёт",
        kill="если достижим только изотропный пек И цена < $800/т — "
             "маржа < 0 (pess-угол): вариант умирает, остаётся топливо",
        verdict="alive",
        verdict_note="база (мезофаза $1500): маржа $%d/т, EBITDA $%.0fM/г, "
                     "окуп. %.1f г — самый большой приз волны; но вердикт "
                     "живёт РОВНО за гейтом качества: изотроп-угол "
                     "-$%d/т — бинарный исход"
                     % (sc["base"]["margin_usd_unit"],
                        sc["base"]["ebitda_musd_y"],
                        sc["base"]["payback_y"],
                        -sc["pess"]["margin_usd_unit"]))


def l15_hppo():
    a = {
        "po_price_usd_t": A([1500, 2200], "$/т", "market",
            "пропиленоксид, импортозамещение РФ"),
        "c3h6_usd_t": A([800, 1000], "$/т", "market",
            "свой пропилен (Ставролен)"),
        "c3h6_t_per_t": A(0.76, "т C3H6/т PO", "literature",
            "стехиометрия 0.72 + селективность ~95%"),
        "h2o2_t_per_t": A(0.76, "т H2O2(100%)/т PO", "literature"),
        "h2o2_usd_t": A([600, 900], "$/т", "market",
            "покупной; СТЫК N7: on-site прямой синтез -> $400-550"),
        "opex_usd_t": A([150, 250], "$/т PO", "company_assumption"),
        "line_kt_y": A([50, 100], "кт PO/год", "market",
            "мини-линия против world-scale 300+ кт — капекс/т тяжёлый"),
        "capex_musd": A([150, 250], "$M", "market", "E1 ±3-5x"),
    }

    def run(po, c3, h2o2, opx, kt, capex):
        return econ(po, 0.76 * c3 + 0.76 * h2o2, opx, kt * 1e3, capex)

    sc = dict(
        pess=run(1500, 1000, 900, 250, 50, 150),
        base=run(1850, 900, 750, 200, 75, 200),
        opt=run(2200, 800, 400, 150, 100, 250))
    n7 = run(1850, 900, 450, 200, 75, 200)
    return dict(
        owner="Лукойл", name="L15 HPPO: пропиленоксид на Ставролене",
        link="стык N7 (on-site H2O2) + TS-1 эпоксидирование — наша "
             "спин/пероксо-химия (Ti-OOH интермедиат)",
        assumptions=a, unit="$/т PO", scenarios=sc,
        n7_upside=dict(n7, note="база, но H2O2 on-site $450 (N7): маржа "
                       "$%d/т, окуп. %s г — жив только этот вариант"
                       % (n7["margin_usd_unit"], n7["payback_y"])),
        anchors_used=dict(mining_x3="PO $1500-2200 >> $751 — specialty"),
        gate="ГЛАВНЫЙ ГЕЙТ: дешёвый on-site H2O2 (N7, <= $550/т) — при "
             "покупном $750 маржа $%d/т и окуп. %.0f лет: мёртво" %
             (sc["base"]["margin_usd_unit"], sc["base"]["payback_y"]),
        kill="H2O2 > $850/т при PO < $1700/т — маржа < 0; или мини-масштаб "
             "не строится дешевле ~$2500/т мощности",
        verdict="borderline",
        verdict_note="база на ПОКУПНОМ H2O2: маржа $%d/т, окуп. %.0f лет — "
                     "мертво; с N7-перекисью $450 оживает до $%d/т "
                     "(окуп. %s г) — идея-сателлит N7, сама по себе не "
                     "стоит" % (sc["base"]["margin_usd_unit"],
                                sc["base"]["payback_y"],
                                n7["margin_usd_unit"], n7["payback_y"]))


# ======================================================== НОРНИКЕЛЬ N13-N17

def n13_ni_electrodes():
    a = {
        "price_usd_m2": A([200, 600], "$/м2", "market",
            "активированный (Raney-Ni / NiMo-покрытие) электрод для AWE; "
            "голая Ni-сетка ~$100-150 — премия за активацию и ресурс"),
        "ni_kg_m2": A([0.8, 2.5], "кг Ni/м2", "literature",
            "пена/сетка + покрытие"),
        "ni_usd_kg": A([15, 19], "$/кг", "market", "LME + передел; Ni свой"),
        "processing_usd_m2": A([60, 150], "$/м2", "company_assumption",
            "плазменное/гальваническое нанесение, выщелачивание Al, QC"),
        "line_k_m2_y": A([10, 50], "тыс. м2/год", "market",
            "рынок AWE растёт с ГВт-заводами электролизёров"),
        "capex_musd": A([5, 15], "$M", "market"),
        "oer_eta_V": A(0.435, "В", "anchor_repo",
            "h2_oer_ladder_results: НАША OER-лестница на Ni(OH)2 в щёлочи "
            "— прямая расчётная валидация продукта (редкий случай: якорь "
            "репозитория лежит ровно на идее)"),
    }

    def run(price, nikg, ni, proc, km2, capex):
        return econ(price, nikg * ni + proc, 0.0, km2 * 1e3, capex)

    sc = dict(
        pess=run(200, 2.5, 19, 150, 10, 5),
        base=run(400, 1.5, 17, 100, 30, 10),
        opt=run(600, 0.8, 15, 60, 50, 15))
    return dict(
        owner="Норникель", name="N13 Активированные Ni-электроды для "
                                "щелочных электролизёров",
        link="свой Ni + наша OER-лестница (eta=0.435 В) как расчётный "
             "фундамент дизайна покрытия — самая 'наша' идея волны",
        assumptions=a, unit="$/м2", scenarios=sc,
        anchors_used=dict(oer_eta="0.435 В из h2_oer_ladder_results — "
                          "сверено в anchors_check"),
        gate="ГЛАВНЫЙ ГЕЙТ: ресурс покрытия (10 000+ ч, циклы/реверсы "
             "поляризации) + вход в цепочку OEM-электролизёрщиков — "
             "премию платят за подтверждённые часы, не за eta на бумаге",
        kill="цена < $150/м2 (коммодитизация китайскими сетками без "
             "подтверждённого ресурса) — маржа уходит < 0",
        verdict="alive",
        verdict_note="база: маржа $%d/м2, EBITDA $%.1fM/г, окуп. %.1f г; "
                     "pess-угол -$%d/м2 (активированный по цене голой "
                     "сетки — маловероятная комбинация); малый капекс, "
                     "быстрый вход" % (sc["base"]["margin_usd_unit"],
                                       sc["base"]["ebitda_musd_y"],
                                       sc["base"]["payback_y"],
                                       -sc["pess"]["margin_usd_unit"]))


def n14_cu2o_antifouling():
    a = {
        "fleet_vessels": A([6, 12], "судов", "company_assumption",
            "свой флот СМП (Arc7 контейнеровозы + танкеры/баржи); "
            "проверить у компании"),
        "fuel_t_y_vessel": A([5000, 12000], "т топлива/судно/год",
            "company_assumption"),
        "fuel_usd_t": A([450, 650], "$/т", "market"),
        "incr_saving": A([0.03, 0.07], "доля", "literature",
            "брутто чистый корпус даёт 5-15%; ИНКРЕМЕНТ к обычной краске "
            "— честно 3-7%"),
        "cu2o_t_per_docking": A([5, 15], "т Cu2O/судно/докование",
            "literature"),
        "docking_interval_y": A(3, "лет", "literature"),
        "program_cost_musd_y": A([0.3, 0.8], "$M/год", "company_assumption",
            "краска (свой Cu -> Cu2O) + нанесение в доке"),
        "capex_musd": A(1.5, "$M", "market", "малый узел Cu->Cu2O + склад"),
    }

    def run(n, fuel_t, fuel, sav, prog):
        benefit = n * fuel_t * fuel * sav / 1e6
        return dict(benefit_musd_y=round(benefit, 2),
                    net_musd_y=round(benefit - prog, 2),
                    cu2o_t_y=round(n / 3.0 * 10, 1),
                    payback_y=round(1.5 / (benefit - prog), 1)
                    if benefit > prog else None)

    sc = dict(
        pess=run(6, 5000, 450, 0.03, 0.8),
        base=run(8, 8000, 550, 0.05, 0.5),
        opt=run(12, 12000, 650, 0.07, 0.3))
    return dict(
        owner="Норникель", name="N14 Cu2O-антифоулинг для флота СМП",
        link="своя медь + свой флот: замкнутая мини-цепочка; upside — "
             "продажа краски другим операторам СМП",
        assumptions=a, unit="$M/год на флот (не $/т — экономика услуги)",
        scenarios=sc,
        gate="ГЛАВНЫЙ ГЕЙТ: ледовая абразия — стойкость Cu2O-покрытия "
             "под льдом Arc7 (обычный антифоулинг лёд сдирает за "
             "навигацию); нужен ледостойкий связующий матрикс",
        kill="покрытие не живёт одну навигацию во льду — программа "
             "теряет смысл (перекраска чаще докования невозможна)",
        verdict="borderline",
        verdict_note="ЧЕСТНО МИКРО: Cu2O ~%s т/г (не тоннажный продукт), "
                     "чистый эффект базы $%.1fM/г, окуп. %s г — для НН "
                     "нематериально; жив как практика флота, не как "
                     "бизнес" % (sc["base"]["cu2o_t_y"],
                                 sc["base"]["net_musd_y"],
                                 sc["base"]["payback_y"]))


def n15_an_explosives():
    a = {
        "an_delivered_usd_t": A([400, 700], "$/т", "company_assumption",
            "альтернатива — северный завоз аммиачной селитры/эмульсий "
            "на Таймыр (цена delivered); проверить у компании"),
        "nh3_kt_y": A([20, 50], "кт NH3/год", "market",
            "мини-Haber-Bosch на газе Таймыра — КЛАССИКА малого масштаба, "
            "НЕ наша электрохимия (честно)"),
        "nh3_cash_usd_t": A([250, 450], "$/т NH3", "literature",
            "малый HB на своём дешёвом газе: газ ~$50-80/т NH3 + "
            "opex малого масштаба"),
        "nh3_t_per_t_an": A(0.43, "т NH3/т АС", "literature",
            "0.213 в АС + 0.21 через HNO3"),
        "conv_opex_usd_t": A([80, 160], "$/т АС", "literature",
            "узел HNO3 + нейтрализация + приллирование/эмульсия"),
        "capex_musd": A([120, 250], "$M", "market",
            "мини-NH3 + HNO3 + АС: капекс малого масштаба — ГЛАВНАЯ "
            "неопределённость (E1 ±3-5x); модульные мини-HB могут сдвинуть"),
        "tea_anchor_usd_kg_n": A(2.0, "$/кг N", "anchor_repo",
            "flagship #1 (screening): ~$2/кг-N — планка стоимости "
            "азота, доставленного в Арктику"),
    }

    def run(deliv, nh3kt, nh3, opx, capex):
        an_t = nh3kt * 1e3 / 0.43
        r = econ(deliv, 0.43 * nh3, opx, an_t, capex)
        r["usd_kg_n"] = round(r["full_usd_unit"] / 345.0, 2)  # АС 34.5% N
        return r

    sc = dict(
        pess=run(450, 20, 450, 160, 250),
        base=run(550, 30, 350, 120, 180),
        opt=run(700, 50, 250, 80, 120))
    return dict(
        owner="Норникель", name="N15 Селитра/эмульсионные ВВ on-site "
                                "(мини-Haber-Bosch, Таймыр)",
        link="газ Таймыра -> NH3 -> АС/эмульсии для СВОИХ рудников; "
             "доход = замещение северного завоза, не рынок",
        assumptions=a, unit="$/т АС (цена = завоз delivered)", scenarios=sc,
        anchors_used=dict(
            tea_flagship="full базы = $%.2f/кг-N против якоря ~$2/кг-N "
                         "(flagship #1) — сходится по порядку; завоз "
                         "$550/т = $1.59/кг-N" % sc["base"]["usd_kg_n"],
            mining_x3="АС $400-700 < $751 — НЕ specialty: это "
                      "замещение затрат, не продукт с премией"),
        gate="ГЛАВНЫЙ ГЕЙТ: капекс мини-HB (E1 ±3-5x) — при модульном "
             "<$1500/т мощности и завозе $650+ оживает; плюс безопасность "
             "снабжения рудников (стратегическая, не в цифрах)",
        kill="завоз стабильно < $450/т (сильный рубль/субсидия) или "
             "капекс > $250M — full выше завоза навсегда",
        verdict="borderline",
        verdict_note="база: full $%d/т ПРОТИВ завоза $550 — маржа $%d/т, "
                     "около нуля; opt (модульный капекс + завоз $700) "
                     "даёт $%d/т и окуп. %s г; решает капекс, не химия"
                     % (sc["base"]["full_usd_unit"],
                        sc["base"]["margin_usd_unit"],
                        sc["opt"]["margin_usd_unit"],
                        sc["opt"]["payback_y"]))


def n16_ni63_betavoltaics():
    a = {
        "source_price_usd_unit": A([1000, 5000], "$/шт", "market",
            "микроисточник мкВт-класса (по аналогам NanoTritium/City Labs)"),
        "arctic_iot_units_y": A([1000, 100000], "шт/год", "market",
            "датчики трубопроводов/геотехника вечной мерзлоты/СМП-буи — "
            "рамка ниши, не прогноз"),
        "ni62_feed": A("свой Ni (изотоп Ni-62 ~3.6%)", "-",
            "company_assumption", "обогащение и облучение — только Росатом"),
        "partner": A("Росатом (облучение в реакторе, изотопная лицензия)",
            "-", "company_assumption"),
        "half_life_y": A(101.0, "лет", "literature",
            "Ni-63: бета без гаммы — идеален для необслуживаемых датчиков"),
    }
    return dict(
        owner="Норникель", name="N16 Ni-63 бета-вольтаика (КОНЦЕПТ с "
                                "Росатомом)",
        link="свой никель как сырьё изотопа; арктический IoT — свой же "
             "спрос (мерзлота/трубы/СМП)",
        assumptions=a, unit="рамка рынка, БЕЗ unit economics",
        concept_no_economics=True,
        market_frame="ниша $%dM-$%dM/год в самом широком допущении "
                     "(1-100 тыс. шт по $1000-5000) — но себестоимость "
                     "обогащения Ni-62 и облучения знает только Росатом: "
                     "числа не наши, модель честно молчит"
                     % (1000 * 1000 // 10**6, 100000 * 5000 // 10**6),
        gate="ГЛАВНЫЙ ГЕЙТ: партнёрство с Росатомом (реакторное облучение "
             "+ изотопная лицензия) — без него идеи не существует",
        kill="Росатом не входит / лицензия на изотопное производство вне "
             "контура Росатома невозможна — концепт закрывается",
        verdict="concept",
        verdict_note="concept_no_economics: в скрине только рамка рынка; "
                     "не считаем alive/dead без себестоимости изотопа")


def n17_fecl3_coagulant():
    a = {
        "price_usd_t": A([250, 450], "$/т 40%-р-ра", "market",
            "коагулянт для водоканалов (замещение импорта/Al-соли)"),
        "line_kt_y": A([20, 60], "кт/год (40% р-р)", "market"),
        "cl2_t_per_t": A(0.262, "т Cl2/т 40%-р-ра", "literature",
            "стехиометрия FeCl3: 65.6% Cl x 0.4"),
        "cl2_own_usd_t": A([0, 150], "$/т Cl2", "company_assumption",
            "хлор Кольской — свой, местами проблемный актив (кредит "
            "утилизации возможен); проверить у компании"),
        "fe_units_usd_t": A([30, 90], "$/т р-ра", "company_assumption",
            "железо из пирротиновых хвостов (T3): механоактивация + "
            "выщелачивание/хлорирование"),
        "opex_usd_t": A([40, 80], "$/т", "company_assumption"),
        "delivery_usd_t": A([30, 120], "$/т", "market",
            "40% р-р = возим 60% воды: радиус ~1000-1500 км (СЗФО)"),
        "capex_musd": A([15, 40], "$M", "market"),
    }

    def run(price, cl2, fe, opx, deliv, kt, capex):
        return econ(price, 0.262 * cl2 + fe, opx, kt * 1e3, capex,
                    delivery=deliv)

    sc = dict(
        pess=run(250, 150, 90, 80, 120, 20, 15),
        base=run(350, 75, 60, 60, 60, 40, 25),
        opt=run(450, 0, 30, 40, 30, 60, 30))
    return dict(
        owner="Норникель", name="N17 FeCl3-коагулянт из пирротиновых "
                                "хвостов + хлора Кольской",
        link="ДВОЙНАЯ утилизация: железо хвостов (T3, механовскрытие "
             "Fe-S — наш COGEF-якорь) и хлор Кольской в один продукт",
        assumptions=a, unit="$/т 40%-раствора", scenarios=sc,
        anchors_used=dict(
            fes_cogef="F_max(Fe-S) = 1.97 нН << C-C 5.58 нН "
                      "(calc/fes_cogef_results): механохимическое "
                      "вскрытие пирротиновой решётки доступно — гейт T3 "
                      "пройден расчётом, сверено в anchors_check",
            mining_x3="FeCl3 $250-450 < $751 — НЕ specialty; живёт "
                      "только на почти бесплатном сырье из отходов"),
        gate="ГЛАВНЫЙ ГЕЙТ: логистика раствора (радиус сбыта) + "
             "СанПиН-примеси (Ni/Co/As из хвостов в питьевой воде — "
             "сертификация коагулянта)",
        kill="доставка > $120/т (возим воду дальше 1500 км) или "
             "непрохождение по примесям для питьевого применения",
        verdict="borderline",
        verdict_note="база: маржа $%.0f/т, EBITDA $%.1fM/г, окуп. %.1f г "
                     "— живо, но тонко и некрупно; красота идеи в "
                     "утилизации двух проблемных потоков, не в марже"
                     % (sc["base"]["margin_usd_unit"],
                        sc["base"]["ebitda_musd_y"],
                        sc["base"]["payback_y"]))


# ================================================================== main

def main():
    ideas = {
        "L11_dmds": l11_dmds(),
        "L12_ods": l12_ods(),
        "L13_co2_eor": l13_co2_eor(),
        "L14_mesophase_pitch": l14_pitch(),
        "L15_hppo": l15_hppo(),
        "N13_ni_electrodes": n13_ni_electrodes(),
        "N14_cu2o_antifouling": n14_cu2o_antifouling(),
        "N15_an_explosives": n15_an_explosives(),
        "N16_ni63_betavoltaics": n16_ni63_betavoltaics(),
        "N17_fecl3_coagulant": n17_fecl3_coagulant(),
    }

    # -------- сверка якорей репозитория (обязана сходиться, иначе врём)
    h2 = json.load(open(os.path.join(DIR, "h2_oer_ladder_results.json")))
    fes = json.load(open(os.path.join(DIR, "..", "calc",
                                      "fes_cogef_results.json")))
    msa = json.load(open(os.path.join(DIR, "msa_cost_model.json")))
    mining = json.load(open(os.path.join(DIR,
                                         "econ_vs_mining_results.json")))
    anchors = {
        "oer_eta_V": dict(
            model=0.435, anchor=h2["ladder"]["eta_V"],
            source="h2_oer_ladder_results.json: eta OER на Ni(OH)2 в "
                   "щёлочи — фундамент N13"),
        "fes_cogef_F_max_nN": dict(
            model=1.97, anchor=fes["cogef"]["F_max_nN"],
            source="calc/fes_cogef_results.json: разрыв Fe-S(тиолат) — "
                   "механовскрытие хвостов для N17 (T3)"),
        "mining_3x_usd_t": dict(
            model=MINING_X3_USD_T,
            anchor=mining["verdict"]["breakeven_price_to_beat_mining_3x"],
            source="econ_vs_mining_results.json: порог specialty; "
                   "проходят L11/L14/L15, НЕ проходят N15/N17 (помечено)"),
        "capital_charge_formula": dict(
            model=round(4.5e6 * WACC / 1093.0),
            anchor=msa["per_tonne_inputs"]["capital_usd"],
            source="msa_cost_model.json: $4.5M x 12% / 1093 т — та же "
                   "формула капплатежа, что в econ() этого скрипта"),
    }
    for a in anchors.values():
        a["ok"] = abs(a["model"] - a["anchor"]) <= max(
            0.005 * abs(a["anchor"]), 0.005)

    # -------- сводка вердиктов + топ-3 (детерминированный скоринг)
    verdicts = {k: v["verdict"] for k, v in ideas.items()}
    alive = [k for k, v in ideas.items() if v["verdict"] == "alive"]
    # скоринг только для alive: EBITDA/окупаемость базы, дисконт x0.5
    # если pess-угол уходит в минус (устойчивость)
    scores = {}
    for k in alive:
        b = ideas[k]["scenarios"]["base"]
        s = b["ebitda_musd_y"] / b["payback_y"]
        if ideas[k]["scenarios"]["pess"]["margin_usd_unit"] < 0:
            s *= 0.5
        scores[k] = round(s, 2)
    top3 = sorted(scores, key=lambda k: -scores[k])[:3]

    out = {
        "model": "Третья волна скрининга: 10 идей (Лукойл L11-L15, "
                 "Норникель N13-N17) — один скрипт, unit economics одной "
                 "установки/линии на идею",
        "tier": "screening_pm_3_5x",
        "tier_note": "ГРУБЕЕ второй волны (±2-3x): по 5-8 параметров на "
                     "идею, три угла (pess/base/opt), без implementation-"
                     "блоков; цель — отсев, не проектирование",
        "not_bankable": True,
        "tier_legend": {
            "anchor_repo": "число из наших файлов или физика (помечено)",
            "market": "рыночный диапазон, порядок величины",
            "company_assumption": "допущение о компании — проверить у неё",
            "literature": "литература/индустрия, не наш расчёт"},
        "anchors_check": anchors,
        "ideas": ideas,
        "summary": {
            "verdicts": verdicts,
            "counts": {v: list(verdicts.values()).count(v)
                       for v in ("alive", "borderline", "dead", "concept")},
            "top3_by_attractiveness": top3,
            "top3_scores_ebitda_per_payback": scores,
            "top3_note": "скоринг = EBITDA базы / окупаемость базы, "
                         "дисконт x0.5 если pess-угол в минусе; только "
                         "alive-идеи",
        },
        "honesty": "скрин ±3-5x, не банк. Честные выводы: (1) НИ ОДНОЙ "
                   "мёртвой при базовых допущениях — подозрительно добрый "
                   "результат для третьей волны; причина в том, что все "
                   "borderline живут на одной непроверенной ручке "
                   "(спред L12, налог L13, H2O2-цена L15, капекс N15, "
                   "логистика N17) — каждый гейт способен убить; "
                   "(2) alive-тройка неоднородна: L14 — самый большой "
                   "приз за БИНАРНЫМ гейтом качества (π-радикальная "
                   "конденсация, наша NEVPT2-химия), L11 — единственная "
                   "идея с положительным pess-углом, N13 — единственная с "
                   "прямым расчётным якорем репозитория (eta 0.435 В); "
                   "(3) N14 честно микро ($1-2M/г), N16 честно концепт "
                   "без экономики; (4) N15/N17 ниже порога specialty "
                   "$751/т — живут только на бесплатном/отходном сырье и "
                   "замещении завоза, это другой класс идей; (5) стыки "
                   "работают в обе стороны: L12 и L15 добавляют спрос на "
                   "H2O2 (усиливают кейс N7), но сами без N7 не стоят.",
    }

    # ------------------------------------------------------------ сводка
    print("=" * 74)
    print("ТРЕТЬЯ ВОЛНА: 10 идей одним скрином (±3-5x, не банк)")
    print("=" * 74)
    print("сверка якорей:", "OK" if all(a["ok"] for a in anchors.values())
          else "РАСХОЖДЕНИЕ! " + str(anchors))
    print()
    for k, v in ideas.items():
        b = v.get("scenarios", {}).get("base", {})
        num = ("маржа $%s/ед, EBITDA $%s M/г, окуп. %s г"
               % (b.get("margin_usd_unit", "-"),
                  b.get("ebitda_musd_y", "-"), b.get("payback_y", "-"))
               if b and "margin_usd_unit" in b else
               (("net $%s M/г, окуп. %s г" % (b["net_musd_y"],
                                              b["payback_y"]))
                if b else "concept_no_economics"))
        print("%-24s %-10s %s" % (k, v["verdict"].upper(), num))
        print("%24s гейт: %s" % ("", v["gate"].replace(
            "ГЛАВНЫЙ ГЕЙТ: ", "")[:95]))
    print()
    print("топ-3 привлекательности (среди alive):",
          " > ".join("%s (%.2f)" % (k, scores[k]) for k in top3))
    print("счёт вердиктов:", out["summary"]["counts"])

    path = os.path.join(DIR, "econ_third_wave_screen_results.json")
    with open(path, "w") as f:
        json.dump(out, f, indent=1, ensure_ascii=False)
    print("\nwrote %s" % os.path.basename(path))


if __name__ == "__main__":
    main()
