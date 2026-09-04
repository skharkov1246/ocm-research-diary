#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
routes/econ_nornickel_pd_h2o2.py — Норникель N7: распределённый H2O2 на
Pd-катализаторе — арктическая водоочистка + инженерия спроса на палладий.

Двойная история:
 (а) H2O2 НА МЕСТЕ (скид прямого синтеза H2+O2) вместо привозной перекиси:
     в Арктику перекись едет дорого (класс опасности 5.1, обогрев, сезонная
     логистика) — базовая цена $600-900/т (100%) + северная надбавка
     $200-500/т. Скид замещает delivered-тонну своей себестоимостью.
 (б) Каждый скид создаёт НОВЫЙ спрос на Pd (single-atom / эгг-шелл,
     0.1-1 кг/скид). Для №1 производителя Pd в мире это ВИТРИНА применения,
     а не выручка: спрос-эффект в тоннах — пыль (см. pd_demand_honest),
     ценность — нарратив «Pd чистит воду в Арктике» после заката
     автокатализаторов.

Технология: прямой синтез H2+O2 в микрореакторе, вне взрывного окна
(4-96 об.% H2 в O2) за счёт мембранного дозирования/разбавления — это
ГЛАВНАЯ инженерная стена трека, и квант её не двигает (честно, как во
вкладке h2o2-direct). H2 — электролиз на месте (энергия НТЭК) или от
пиролиза N3 (econ_nornickel_pyro_nsr) — синергия в отдельном поле.

Наш вклад (trl_gates): перенос Pd1-дескриптора 2e⁻-селективности (связать
*OOH, не разрывая O-O) на дешёвые 3d-центры; Mn — робастный спин-свитч из
спинового транзистора (switch=True в PBE И PBE0), field-gate барьера
оксидации dBar/dF = -1.06 (PBE) / -1.54 эВ/(В/Å) (PBE0) — ИСПРАВЛЕННЫЕ
v2-числа (spin_field_oxidation_matrix_results.json, ядерный член поля
доначислен). Вариант B «Pd-free» — хедж; честно: скрин Fe/Co/Ni-N4 ещё
в очереди AWS, реальных чисел по 3d-центрам нет.

Скрин-уровень ±2-3x, НЕ банковская модель. Данные Норникеля (потребление
площадок, северная надбавка, цена энергии) — ТОЛЬКО диапазоны
tier='company_assumption' с пометкой «проверить у компании».

Запуск: python3 routes/econ_nornickel_pd_h2o2.py
     -> routes/econ_nornickel_pd_h2o2_results.json
"""
import json
import math
import os

DIR = os.path.dirname(os.path.abspath(__file__))

M_H2O2, M_H2, M_O2 = 34.014, 2.016, 31.998

# ------------------------------------------------------ допущения (все крутилки)
ASSUMPTIONS = {
    # --- физика / якоря репозитория ---
    "h2_stoich_t_per_t": dict(value=round(M_H2 / M_H2O2, 4),
        unit="т H2/т H2O2 (стехиометрия)", tier="anchor_repo",
        source="физика: H2+O2->H2O2, 2.016/34.014=0.0593 — реальный расход "
               "выше на 1/селективность (несел. путь жжёт H2 в воду)"),
    "o2_stoich_t_per_t": dict(value=round(M_O2 / M_H2O2, 4),
        unit="т O2/т H2O2 (стехиометрия)", tier="anchor_repo",
        source="физика: 31.998/34.014=0.9407; реальный расход /селективность"),
    "wacc": dict(value=0.12, unit="доля/год", tier="anchor_repo",
        source="дом: skid_twin / econ_nornickel_* — 12% (аморт.+WACC)"),
    "n3_h2_value_usd_t": dict(value=[48, 120], unit="$/т H2",
        tier="anchor_repo",
        source="econ_nornickel_pyro_nsr: H2 на месте ценится замещением "
               "газа по теплотворности 2.4 т CH4-экв/т H2 x газ НГП "
               "$20-50/т = $48-120/т ($0.05-0.12/кг)"),
    "mn_field_gate_ev_per_va": dict(value={"pbe": -1.06, "pbe0": -1.54},
        unit="эВ/(В/Å)", tier="anchor_repo",
        source="spin_field_oxidation_matrix_results.json (ВЕРСИЯ 2, "
               "коррекция 25.07 — ядерный член поля доначислен): "
               "dБарьер/dF Mn = -1.06 (PBE) / -1.54 (PBE0, R2=0.75, "
               "шумнее); v1-числа были завышены ~2x"),
    "mn_spin_switch_robust": dict(value=True, unit="-", tier="anchor_repo",
        source="spin_field_o2_mn(_pbe0)_results.json: полевой флип "
               "основного спина O2-комплекса Mn — switch=True в PBE И "
               "PBE0 (dGap/dF -1.26/-0.73 эВ/(В/Å)); Cu/Co/Ni хрупче"),
    # --- рынок ---
    "h2o2_base_price_usd_t": dict(value=[600, 900],
        unit="$/т H2O2 (100%)", tier="market",
        comment="балк-цена перекиси (в пересчёте на 100%) до северной "
                "надбавки; kill-порог delivered <$700/т"),
    "o2_onsite_usd_t": dict(value=[40, 80], unit="$/т O2", tier="market",
        comment="PSA/VSA на месте; у метзавода может быть свой кислород "
                "дешевле — апсайд, company_assumption"),
    "pd_per_skid_kg": dict(value=[0.1, 1.0], unit="кг Pd/скид",
        tier="market",
        comment="single-atom/эгг-шелл Pd: загрузка нарочно мала — "
                "в этом и пуанта технологии (грамм вместо килограммов)"),
    "pd_price_usd_kg": dict(value=[30000, 36000], unit="$/кг Pd",
        tier="market", comment="~$950-1120/тр.унц."),
    "nn_pd_output_t_y": dict(value=[80, 90], unit="т Pd/год",
        tier="market", comment="выпуск Норникеля ~2.6-2.8 Мунц/год — "
                "публичные отчёты, порядок величины"),
    "global_pd_demand_t_y": dict(value=[270, 300], unit="т Pd/год",
        tier="market", comment="мировой спрос, ~80% — автокатализаторы "
                "(закат с электрификацией — мотив искать новый спрос)"),
    # --- допущения о компании (ПРОВЕРИТЬ У КОМПАНИИ) ---
    "arctic_surcharge_usd_t": dict(value=[200, 500], unit="$/т H2O2",
        tier="company_assumption",
        comment="северная надбавка delivered Норильск: класс 5.1, "
                "обогреваемые ёмкости, сезонность СМП/река, страховка — "
                "проверить у компании фактические закупочные цены"),
    "site_demand_kt_y": dict(value=[0.5, 3.0], unit="кт H2O2(100%)/год",
        tier="company_assumption",
        comment="потребление площадок НН: флотация (окислитель), детокс "
                "стоков (SO2/тиосоли, органика), водоподготовка/очистные "
                "— проверить у компании; kill при <300 т/год"),
    "power_usd_kwh": dict(value=[0.02, 0.05], unit="$/кВт*ч",
        tier="company_assumption",
        comment="энергия НТЭК (гидро Усть-Хантайская/Курейская + газ) — "
                "проверить внутренний тариф"),
    # --- литература ---
    "h2_selectivity": dict(value=[0.6, 0.8], unit="доля",
        tier="literature",
        comment="селективность прямого синтеза по H2 на Pd/PdAu(Sn); "
                "остальной H2 сгорает в воду / разлагает готовый H2O2 — "
                "ровно то, что чинит наш 2e⁻-дескриптор"),
    "skid_capacity_t_y": dict(value=[100, 500], unit="т H2O2(100%)/год",
        tier="literature",
        comment="микрореакторные модули прямого синтеза — пилоты/TEA; "
                "промышленных референсов скид-класса нет (сам факт — риск)"),
    "lit_skid_cost_usd_t": dict(value=[400, 700], unit="$/т H2O2",
        tier="literature",
        comment="TEA прямого синтеза + масштабный штраф; ВНИМАНИЕ: наша "
                "bottom-up себестоимость выше (см. lit_cost_tension) — "
                "лит-TEA предполагают дешёвый H2 и мягкий капитал"),
    "elec_kwh_kg_h2": dict(value=[50, 55], unit="кВт*ч/кг H2",
        tier="literature", comment="электролиз, полный (стек+BoP)"),
    "explosive_window_pct": dict(value=[4, 96], unit="об.% H2 в O2",
        tier="literature",
        comment="взрывное окно H2/O2 — ГЛАВНАЯ инженерная стена трека; "
                "обход: мембранное дозирование O2 / разбавление N2/CO2; "
                "квант её НЕ двигает (honesty вкладки h2o2-direct)"),
    "pd_lifetime_y": dict(value=[3, 5], unit="лет", tier="literature",
        comment="жизнь Pd-катализатора до перегрузки (лиганд-выщелачивание, "
                "спекание) — определяет ГОДОВОЙ Pd-спрос замещения"),
    "h2o2_freeze_c": dict(value=-33, unit="°C", tier="literature",
        comment="35 wt% замерзает ~-33°C, в Норильске бывает ниже — "
                "обогреваемое хранение и у нас, и у привозной (учтено "
                "в надбавке за доставку)"),
    # --- инженерные оценки E1 (наши, порядок величины) ---
    "capex_skid_musd": dict(value=[2, 5], unit="$M/скид", tier="E1",
        comment="электролизёр + микрореактор + PSA O2 + концентрирование "
                "+ взрывозащита/АСУ + северное исполнение; серийный скид "
                "— к низу, первый — к верху"),
    "fixed_skid_usd_t": dict(value=[100, 250], unit="$/т H2O2", tier="E1",
        comment="вахта/обслуживание/накладные малого масштаба при "
                "дистанционном мониторинге куста скидов"),
    "other_var_usd_t": dict(value=[80, 150], unit="$/т H2O2", tier="E1",
        comment="стабилизаторы, кислотный промотор/галогенид, мембраны, "
                "деионизованная вода; Pd-амортизация внутри — шум "
                "($1-9/т, см. pd_demand_honest)"),
    "h2_electrolysis_usd_t": dict(value=[2000, 4000], unit="$/т H2",
        tier="E1",
        comment="из 50-55 кВт*ч/кг x $0.02-0.05/кВт*ч = $1.0-2.75/кг "
                "энергия + стек/BoP ~$0.8-1.3/кг; сравни N3-водород "
                "$48-120/т — синергия на порядок+"),
}


def V(k):
    return ASSUMPTIONS[k]["value"]


def h2o2_cost(h2_usd_t, sel, o2_usd_t, other_usd_t, fixed_usd_t,
              t_y, capex_musd):
    """Себестоимость т H2O2(100%) у ворот скида, bottom-up."""
    h2_t = V("h2_stoich_t_per_t") / sel
    o2_t = V("o2_stoich_t_per_t") / sel
    var = h2_t * h2_usd_t + o2_t * o2_usd_t + other_usd_t
    cash = var + fixed_usd_t
    capital = capex_musd * 1e6 * V("wacc") / t_y
    return dict(h2_line=round(h2_t * h2_usd_t, 0),
                cash=round(cash, 0), capital=round(capital, 0),
                full=round(cash + capital, 0))


def scenario(h2_usd_t, sel, o2_usd_t, other_usd_t, fixed_usd_t, t_y,
             capex_musd, base_price, surcharge, site_kt_y, **_):
    """Сценарий: себестоимость vs delivered-альтернатива, экономия на
    тонне и на площадке, вариант с N3-водородом."""
    c = h2o2_cost(h2_usd_t, sel, o2_usd_t, other_usd_t, fixed_usd_t,
                  t_y, capex_musd)
    delivered = base_price + surcharge
    savings = delivered - c["full"]
    site_t_y = site_kt_y * 1000.0
    skids = max(1, math.ceil(site_t_y / t_y))
    savings_musd = savings * site_t_y / 1e6
    capex_total = skids * capex_musd
    payback = capex_total / savings_musd if savings_musd > 0 else None
    # вариант: H2 от пиролиза N3 (замещение газа по теплотворности)
    h2_n3 = 0.5 * (V("n3_h2_value_usd_t")[0] + V("n3_h2_value_usd_t")[1])
    c3 = h2o2_cost(h2_n3, sel, o2_usd_t, other_usd_t, fixed_usd_t,
                   t_y, capex_musd)
    sav3 = delivered - c3["full"]
    return dict(
        cost=c, delivered_usd_t=round(delivered, 0),
        savings_usd_t=round(savings, 0),
        site_t_y=round(site_t_y, 0), skids=skids,
        savings_site_musd_y=round(savings_musd, 2),
        capex_total_musd=round(capex_total, 1),
        payback_y=round(payback, 1) if payback else None,
        procurement_replaced_musd_y=round(delivered * site_t_y / 1e6, 2),
        with_n3_h2=dict(full_usd_t=c3["full"],
                        savings_usd_t=round(sav3, 0),
                        savings_site_musd_y=round(sav3 * site_t_y / 1e6, 2)))


# --------------------------------------------- сценарии (pess / base / opt)
# Пессимизм = логистика вдруг дешёвая (delivered $800), скид мелкий и
# дорогой (капитальная строка $4000/т), селективность по низу. База —
# полноразмерный 500 т/г модуль, серийный капекс $2.5M, середины рынка.
# Оптимизм — селективность 0.8 (наш дескриптор отработал), капекс $2M,
# надбавка по верху.
SCENARIOS = {
    "pessimistic": dict(h2_usd_t=4000.0, sel=0.60, o2_usd_t=70.0,
        other_usd_t=150.0, fixed_usd_t=250.0, t_y=150.0, capex_musd=5.0,
        base_price=600.0, surcharge=200.0, site_kt_y=0.5,
        note="мелкий скид 150 т/г за $5M, H2 $4/кг, сел. 0.60, "
             "delivered $800 — капитальная строка хоронит"),
    "base": dict(h2_usd_t=3000.0, sel=0.70, o2_usd_t=55.0,
        other_usd_t=115.0, fixed_usd_t=175.0, t_y=500.0, capex_musd=2.5,
        base_price=750.0, surcharge=350.0, site_kt_y=1.5,
        note="полный модуль 500 т/г за $2.5M, H2 электролиз $3/кг, "
             "сел. 0.70, delivered $1100, площадки 1.5 кт/г (3 скида)"),
    "optimistic": dict(h2_usd_t=2000.0, sel=0.80, o2_usd_t=40.0,
        other_usd_t=80.0, fixed_usd_t=100.0, t_y=500.0, capex_musd=2.0,
        base_price=900.0, surcharge=500.0, site_kt_y=3.0,
        note="сел. 0.80 (дескриптор отработал), капекс $2M, delivered "
             "$1400, площадки 3 кт/г (6 скидов)"),
}


def main():
    # -------- сверка с якорями репозитория (обязана сходиться)
    mtx = json.load(open(os.path.join(
        DIR, "spin_field_oxidation_matrix_results.json")))
    mn_pbe = json.load(open(os.path.join(
        DIR, "spin_field_o2_mn_results.json")))
    mn_pbe0 = json.load(open(os.path.join(
        DIR, "spin_field_o2_mn_pbe0_results.json")))
    pyro = json.load(open(os.path.join(
        DIR, "econ_nornickel_pyro_nsr_results.json")))
    msa = json.load(open(os.path.join(
        DIR, "econ_nornickel_sulfur_msa_results.json")))

    gas = pyro["assumptions"]["gas_internal_usd_t"]["value"]
    lhv = pyro["assumptions"]["lhv_h2_ch4_mj_kg"]["value"]
    ratio = lhv[0] / lhv[1]                      # 120/50 = 2.4
    anchors = {
        "mn_field_gate_pbe": dict(
            model=V("mn_field_gate_ev_per_va")["pbe"],
            anchor=mtx["combos"]["Mn-pbe"][
                "dBarrier_dF_eV_per_VpA_corrected"],
            source="spin_field_oxidation_matrix (v2 corrected)"),
        "mn_field_gate_pbe0": dict(
            model=V("mn_field_gate_ev_per_va")["pbe0"],
            anchor=mtx["combos"]["Mn-pbe0"][
                "dBarrier_dF_eV_per_VpA_corrected"],
            source="spin_field_oxidation_matrix (v2 corrected, R2=0.75)"),
        "mn_switch_pbe": dict(
            model=V("mn_spin_switch_robust"),
            anchor=mn_pbe["ground_switches_with_field"],
            source="spin_field_o2_mn_results.json"),
        "mn_switch_pbe0": dict(
            model=V("mn_spin_switch_robust"),
            anchor=mn_pbe0["ground_switches_with_field"],
            source="spin_field_o2_mn_pbe0_results.json"),
        "n3_h2_value_lo": dict(
            model=V("n3_h2_value_usd_t")[0],
            anchor=round(gas[0] * ratio, 0),
            source="econ_nornickel_pyro_nsr: газ $%d x 2.4" % gas[0]),
        "n3_h2_value_hi": dict(
            model=V("n3_h2_value_usd_t")[1],
            anchor=round(gas[1] * ratio, 0),
            source="econ_nornickel_pyro_nsr: газ $%d x 2.4" % gas[1]),
        "wacc": dict(
            model=V("wacc"),
            anchor=msa["assumptions"]["wacc"]["value"],
            source="econ_nornickel_sulfur_msa: WACC дома 12%"),
    }
    for a in anchors.values():
        if isinstance(a["model"], bool):
            a["ok"] = a["model"] == a["anchor"]
        else:
            a["ok"] = abs(a["model"] - a["anchor"]) <= max(
                0.01 * abs(a["anchor"]), 0.01)

    scen_out = {name: dict(params=dict(p), **scenario(**p))
                for name, p in SCENARIOS.items()}
    sb = scen_out["base"]
    base = SCENARIOS["base"]

    # -------- честное напряжение с литературной TEA
    lit_cost_tension = dict(
        lit_skid_cost_usd_t=V("lit_skid_cost_usd_t"),
        our_bottom_up_base_usd_t=sb["cost"]["full"],
        gap_explained="лит-TEA $400-700/т достижима только при дешёвом H2 "
            "(<$1/кг: труба/SMR — в Норильске нет) и мягком капитале; наш "
            "bottom-up при E1-капексе $2.5M и WACC 12%% даёт капитальную "
            "строку $%.0f/т + H2-строку $%.0f/т (электролиз) = full $%.0f/т; "
            "с N3-водородом cash падает до $%.0f/т, но капитальная остаётся "
            "— именно она, а не химия, решает экономику"
            % (sb["cost"]["capital"], sb["cost"]["h2_line"],
               sb["cost"]["full"],
               sb["with_n3_h2"]["full_usd_t"] - sb["cost"]["capital"]))

    # -------- чувствительность (вокруг базы)
    def m(**over):
        return scenario(**dict(base, **over))

    h2_n3_mid = 0.5 * (V("n3_h2_value_usd_t")[0]
                       + V("n3_h2_value_usd_t")[1])
    sens = {
        "t_y_per_skid (ГЛАВНАЯ ручка — загрузка/размер скида)": {
            "%d т/г" % t: m(t_y=float(t))["savings_usd_t"]
            for t in (150, 300, 500)},
        "capex_musd (вторая ручка)": {
            "$%.1fM" % c: m(capex_musd=c)["savings_usd_t"]
            for c in (2.0, 2.5, 3.5, 5.0)},
        "delivered_usd_t (цена альтернативы)": {
            "$%d" % d: m(base_price=float(d), surcharge=0.0)[
                "savings_usd_t"]
            for d in (700, 900, 1100, 1400)},
        "h2_usd_t (источник водорода)": {
            "N3 $%d" % h2_n3_mid: m(h2_usd_t=h2_n3_mid)["savings_usd_t"],
            "$2000": m(h2_usd_t=2000.0)["savings_usd_t"],
            "$3000": m(h2_usd_t=3000.0)["savings_usd_t"],
            "$4000": m(h2_usd_t=4000.0)["savings_usd_t"]},
        "selectivity (наш дескриптор)": {
            "%.2f" % s: m(sel=s)["savings_usd_t"]
            for s in (0.60, 0.70, 0.80)},
        "arctic_surcharge_usd_t": {
            "$%d" % s: m(surcharge=float(s))["savings_usd_t"]
            for s in (200, 350, 500)},
        "unit": "экономия $/т H2O2 vs delivered (standalone-электролиз, "
                "если не сказано иное); честно: экономику решают ЗАГРУЗКА "
                "и КАПЕКС (капитальная строка $480-4000/т), потом цена "
                "альтернативы, потом H2; селективность двигает ~$85/т на "
                "0.6->0.8 — важна, но не главная",
    }

    # -------- kill-критерии
    k_dem = m(t_y=300.0, site_kt_y=0.3)
    k_dem3 = scenario(**dict(base, t_y=300.0, site_kt_y=0.3))
    kill = {
        "microreactor_safety_not_certifiable": dict(
            kill_if="взрывобезопасность микрореактора H2/O2 не проходит "
                    "сертификацию (Ростехнадзор, окно 4-96 об.%)",
            tier="literature+company_assumption",
            comment="главная стена трека — инженерная, квант её не "
                    "двигает; мембранное дозирование/разбавление должно "
                    "пройти HAZOP и получить разрешение на опасный "
                    "производственный объект в Норильске; нет прецедента "
                    "— вариант мёртв целиком"),
        "site_demand_below_t_y": dict(kill_below=300,
            computed_context="при потребности 300 т/г скид 300 т/г имеет "
                "капитальную строку $%.0f/т, full $%.0f/т, экономия "
                "$%.0f/т (электролиз) / $%.0f/т (N3) — мёртв; малый спрос "
                "не тянет капекс, а возить 300 т проще, чем строить"
                % (k_dem["cost"]["capital"], k_dem["cost"]["full"],
                   k_dem["savings_usd_t"],
                   k_dem3["with_n3_h2"]["savings_usd_t"])),
        "delivered_below_usd_t": dict(kill_below=700,
            computed_context="если логистика подешевела и delivered "
                "<$700/т — лучший наш full (оптимизм + N3-водород) "
                "$%.0f/т: экономия исчезает во всех сценариях; кейс "
                "живёт ТОЛЬКО на северной надбавке"
                % scen_out["optimistic"]["with_n3_h2"]["full_usd_t"]),
        "verdict_scenarios": {},
    }
    for name, s in scen_out.items():
        alive = (s["delivered_usd_t"] >= 700
                 and SCENARIOS[name]["site_kt_y"] * 1000 >= 300
                 and s["with_n3_h2"]["savings_usd_t"] > 0)
        kill["verdict_scenarios"][name] = (
            "жив только с N3-водородом (эл-з %+.0f, N3 %+.0f $/т)"
            % (s["savings_usd_t"], s["with_n3_h2"]["savings_usd_t"])
            if alive and s["savings_usd_t"] <= 0 else
            ("жив (эл-з %+.0f, N3 %+.0f $/т, EBITDA-экв $%.2fM/г)"
             % (s["savings_usd_t"], s["with_n3_h2"]["savings_usd_t"],
                s["savings_site_musd_y"]) if alive else "МЁРТВ"))

    # -------- синергия с N3 (пиролиз: H2 по цене замещённого газа)
    h2_need = {name: round(s["site_t_y"] * V("h2_stoich_t_per_t")
                           / SCENARIOS[name]["sel"], 0)
               for name, s in scen_out.items()}
    synergy_n3 = {
        "idea": "H2 для прямого синтеза — от пиролиза N3 (сажевый завод), "
                "где H2 ценится замещением газа по теплотворности "
                "$48-120/т, а не электролизом $2000-4000/т: H2-строка "
                "себестоимости падает с $%.0f до $%.0f/т H2O2"
                % (sb["cost"]["h2_line"],
                   sb["with_n3_h2"]["full_usd_t"] - (sb["cost"]["full"]
                   - sb["cost"]["h2_line"])),
        "per_scenario": {name: dict(
            standalone_full=s["cost"]["full"],
            n3_full=s["with_n3_h2"]["full_usd_t"],
            standalone_savings=s["savings_usd_t"],
            n3_savings=s["with_n3_h2"]["savings_usd_t"])
            for name, s in scen_out.items()},
        "h2_need_t_y": h2_need,
        "capacity_check": "потребность в H2 мала (десятки-сотни т/г) — "
                          "N3 даже с одного пиролизного скида даёт "
                          "больше; ограничение не в объёме",
        "caveats": "N3-водород существует, только если N3 построен "
                   "(связка проектов); чистота пиролизного H2 "
                   "(CH4-примесь) для Pd-катализа прямого синтеза — "
                   "проверить в лаборатории (скрин); электролизный H2 "
                   "чище, но в 20-40x дороже",
        "verdict": "синергия флипает базу через ноль: %+.0f -> %+.0f $/т "
                   "— standalone-кейс на электролизе НЕ живёт в базе, "
                   "в связке с N3 живёт" % (
                       sb["savings_usd_t"],
                       sb["with_n3_h2"]["savings_usd_t"]),
    }

    # -------- Pd-спрос: честный счёт (витрина vs тонны)
    pd_kg = V("pd_per_skid_kg")
    life = V("pd_lifetime_y")
    nn_out = V("nn_pd_output_t_y")
    glob = V("global_pd_demand_t_y")

    def fleet(n_skids):
        inv = [n_skids * pd_kg[0], n_skids * pd_kg[1]]
        repl = [inv[0] / life[1], inv[1] / life[0]]
        return dict(skids=n_skids,
            pd_inventory_kg=[round(inv[0], 1), round(inv[1], 1)],
            pd_replacement_kg_y=[round(repl[0], 2), round(repl[1], 1)],
            share_of_nn_output_pct=[
                round(repl[0] / (nn_out[1] * 1000) * 100, 4),
                round(repl[1] / (nn_out[0] * 1000) * 100, 4)])

    pd_cost_per_t = [pd_kg[0] * V("pd_price_usd_kg")[0]
                     / (life[1] * 500.0),
                     pd_kg[1] * V("pd_price_usd_kg")[1]
                     / (life[0] * 100.0)]
    pd_demand_honest = {
        "per_skid": dict(pd_kg=pd_kg, lifetime_y=life,
            pd_amortization_usd_t=[round(pd_cost_per_t[0], 1),
                                   round(pd_cost_per_t[1], 1)],
            note="Pd-амортизация $0.6-120/т H2O2 (обычно единицы $) — "
                 "в себестоимости это шум, сидит в other_var"),
        "norilsk_fleet_3_6_skids": fleet(6),
        "world_showcase_100_skids": fleet(100),
        "world_aggressive_1000_skids": fleet(1000),
        "reference": dict(nn_output_t_y=nn_out,
                          global_demand_t_y=glob,
                          autocat_share="~80% мирового спроса — "
                          "автокатализаторы (закат с электрификацией)"),
        "honest_verdict":
            "спрос-эффект в ТОННАХ — пыль: даже агрессивный мировой флот "
            "1000 скидов = 0.1-1 т Pd в инвентаре и 20-333 кг/год "
            "замещения = 0.02-0.4% ГОДОВОГО выпуска одного Норникеля. И "
            "это правильно: single-atom/эгг-шелл катализ по определению "
            "минимизирует металл. Ценность для №1 производителя Pd — не "
            "тонны, а НАРРАТИВ: новое промышленное применение Pd после "
            "заката автокатализа, витрина «Pd чистит воду в Арктике». "
            "Продавать это как «инженерию спроса» в тоннах было бы ложью.",
        "internal_conflict":
            "наш же хедж — вариант B «Pd-free 3d-центр» (Mn/Fe/Co-N4) — "
            "прямо убивает Pd-нарратив; для Норникеля вариант A (Pd1) — "
            "витрина, вариант B — страховка от ЧУЖОГО Pd-free прорыва: "
            "держать оба и говорить об этом честно",
    }

    # -------- механическая реализация (скрин)
    implementation = {
        "block_flow": [
            "энергия НТЭК -> электролизёр (H2) ЛИБО H2 от пиролиза N3 "
            "(осушка/очистка от CH4 — проверить чистоту для Pd)",
            "O2: PSA/VSA на месте (или кислород метзавода — проверить)",
            "микрореакторный модуль прямого синтеза: мембранное "
            "дозирование O2 в H2-поток / разбавление — ВНЕ взрывного "
            "окна 4-96 об.%; катализатор Pd1/эгг-шелл (вариант B: "
            "Mn/Fe/Co-N4 Pd-free), кислотный промотор+стабилизатор",
            "стабилизация + концентрирование до 35-50 wt%",
            "обогреваемое хранение (35% замерзает ~-33°C)",
            "раздача по потребителям площадки: флотация, детокс стоков, "
            "водоподготовка/очистные",
        ],
        "equipment_per_skid": [
            dict(item="электролизёр 300-500 кВт с BoP (или узел очистки "
                      "N3-водорода — дешевле)", cost_musd=[0.6, 1.2]),
            dict(item="микрореакторный блок + мембранное дозирование O2",
                 cost_musd=[0.6, 1.5]),
            dict(item="PSA/VSA кислород", cost_musd=[0.2, 0.4]),
            dict(item="стабилизация/концентрирование", cost_musd=[0.3, 0.6]),
            dict(item="взрывозащита, АСУ, сертификация, северное "
                      "исполнение", cost_musd=[0.3, 0.8]),
            dict(item="обогреваемое хранение + раздача", cost_musd=[0.2, 0.5]),
        ],
        "equipment_total_musd": [2.2, 5.0],
        "equipment_note": "сумма позиций воспроизводит E1-капекс $2-5M — "
                          "самосогласованность, не независимая оценка",
        "footprint_m2_per_skid": [150, 400],
        "utilities_per_skid_500t": dict(
            electricity_kw=[300, 500],
            note="электролиз 42 т H2/г x 50-55 МВт*ч/т ~ 2.2 ГВт*ч/г "
                 "~ 250 кВт среднее + BoP; с N3-водородом падает в ~10x",
            di_water_m3_d=[1, 3], n2_purge="продувка/разбавление"),
        "staffing": dict(
            pilot="2-3 чел на пилотный скид (вахта) + дистанционный "
                  "мониторинг",
            series="4-6 чел на куст 3-6 скидов", tier="screening"),
        "rollout": [
            dict(stage="FEED + HAZOP + предварительная позиция "
                       "Ростехнадзора по микрореактору H2/O2 + сверка "
                       "потребления площадок с компанией",
                 duration="6-9 мес"),
            dict(stage="расчёты G1-G2 (AWS-скрин 3d-N4 + Mn field-gate "
                       "релаксированная версия)", duration="6-12 мес",
                 note="параллельно FEED"),
            dict(stage="лабораторный прототип микрореактора, 1000 ч вне "
                       "взрывного окна", duration="12-18 мес"),
            dict(stage="сертификация взрывобезопасности (ГЛАВНЫЙ гейт)",
                 duration="12-24 мес", note="параллельно пилоту"),
            dict(stage="пилотный скид 100 т/г на водоподготовке одной "
                       "площадки", duration="12-18 мес"),
            dict(stage="серия 3-6 скидов под фактическое потребление",
                 duration="12-24 мес"),
        ],
        "trl_gates": [
            "G1 (расчёт, наш рычаг): перенос Pd1-дескриптора "
            "2e⁻-селективности (связать *OOH без разрыва O-O) на дешёвые "
            "3d-центры Fe/Co/Ni-N4 — конвейер h2o2_direct_screen.py "
            "собран и проверен (--selftest), но скрин — AWS-класс "
            "нагрузки; ЧЕСТНО: реальных чисел по 3d-центрам ещё НЕТ, "
            "очередь AWS",
            "G2 (расчёт): Mn — робастный спин-свитч (switch=True в PBE и "
            "PBE0) и field-gate барьера оксидации dBar/dF = -1.06 (PBE) "
            "/ -1.54 эВ/(В/Å) (PBE0, R2=0.75) — v2-ИСПРАВЛЕННЫЕ числа "
            "(ядерный член поля доначислен, магнитуда упала ~2x от v1); "
            "идея: полевое подавление 4e⁻-канала/разложения H2O2; "
            "перенос ORR-дескриптора на ТЕРМИЧЕСКИЙ прямой синтез — "
            "не проверен, отдельная задача",
            "G3 (лаборатория): микрореактор 1000 ч вне взрывного окна, "
            "селективность >=60% по H2, концентрация >=8 wt% на выходе",
            "G4 (пилот): 100 т/г на реальной воде площадки, зимняя "
            "эксплуатация, full cost по счётчикам <= $1300/т",
            "G5 (рынок/компания): подтверждённое потребление >= 300 т/г "
            "на скид + сертификат взрывобезопасности + фактическая "
            "delivered-цена привозной перекиси >= $900/т",
        ],
    }

    headline = (
        "Норникель N7: H2O2 на месте против арктической перекиси "
        "delivered $800-1400/т. Экономика живёт только у полноразмерного "
        "загруженного скида (500 т/г) и в связке с N3-водородом "
        "($48-120/т вместо электролизных $2-4k/т): база %+.0f $/т "
        "standalone -> %+.0f $/т с N3. Деньги малы (<= ~$2M/г экономии "
        "на всём контуре) — кейс про автономность снабжения и витрину "
        "нового Pd-спроса: в тоннах пыль (<=0.4%% годового выпуска даже "
        "при 1000 скидов в мире), в нарративе — золото."
        % (sb["savings_usd_t"], sb["with_n3_h2"]["savings_usd_t"]))

    out = {
        "model": "Норникель N7: распределённый H2O2 на Pd-катализаторе — "
                 "арктическая водоочистка + инженерия спроса на палладий",
        "tier": "screening_pm_2_3x",
        "not_bankable": True,
        "tier_legend": {
            "anchor_repo": "число из наших файлов или физика (помечено)",
            "market": "рыночный диапазон, порядок величины",
            "company_assumption": "допущение о компании — проверить у неё",
            "literature": "литература/индустрия, не наш расчёт",
            "E1": "инженерная оценка порядка величины (наша)"},
        "headline": headline,
        "anchors_check": anchors,
        "assumptions": ASSUMPTIONS,
        "scenarios": scen_out,
        "lit_cost_tension": lit_cost_tension,
        "sensitivity": sens,
        "kill_criteria": kill,
        "synergy_n3": synergy_n3,
        "pd_demand_honest": pd_demand_honest,
        "implementation": implementation,
        "honesty": "скрин +-2-3x. Честные выводы: (1) экономику решает "
                   "КАПИТАЛЬНАЯ строка (капекс x загрузка: $480-4000/т), "
                   "а не химия — скид живёт только полноразмерным и "
                   "загруженным; (2) standalone-база на электролизе в "
                   "минусе (%+.0f $/т ~ нуль в скрине), N3-водород "
                   "флипает в %+.0f $/т — кейс живёт в СВЯЗКЕ с N3; (3) "
                   "весь денежный приз <= ~$2M/г — для НН это не "
                   "P&L-кейс, а автономность снабжения + витрина Pd; (4) "
                   "Pd-спрос в тоннах — пыль (см. pd_demand_honest), и "
                   "наш же хедж B (Pd-free) конфликтует с Pd-нарративом "
                   "— говорим об этом прямо; (5) главная стена — "
                   "взрывное окно 4-96%%: инженерия+сертификация, квант "
                   "её НЕ двигает; наш вклад — селективность (Pd1 -> 3d, "
                   "Mn спин-свитч, field-gate -1.06..-1.54 эВ/(В/Å) "
                   "v2-исправленные), и скрин Fe/Co/Ni-N4 ещё В ОЧЕРЕДИ "
                   "AWS — чисел по 3d-центрам пока нет."
                   % (sb["savings_usd_t"],
                      sb["with_n3_h2"]["savings_usd_t"]),
    }

    # ------------------------------------------------------------- сводка
    print("=" * 74)
    print("Норникель N7: распределённый H2O2 на Pd — арктическая "
          "водоочистка")
    print("   + инженерия спроса на Pd; скрин +-2-3x, не банк")
    print("=" * 74)
    print("сверка якорей:", "OK" if all(a["ok"] for a in anchors.values())
          else "РАСХОЖДЕНИЕ! " + str(
              {k: a for k, a in anchors.items() if not a["ok"]}))
    print("\n%-13s%6s%6s%6s%9s%8s%9s%7s" % ("сценарий", "cash", "кап.",
          "full", "deliv.", "эконом.", "с N3", "скидов"))
    for name, s in scen_out.items():
        print("%-13s%6.0f%6.0f%6.0f%9.0f%+8.0f%+9.0f%7d"
              % (name, s["cost"]["cash"], s["cost"]["capital"],
                 s["cost"]["full"], s["delivered_usd_t"],
                 s["savings_usd_t"], s["with_n3_h2"]["savings_usd_t"],
                 s["skids"]))
    print("\nвесь замещаемый закуп (delivered x потребление): "
          "$%.2fM/г (база) .. $%.2fM/г (оптимизм)"
          % (sb["procurement_replaced_musd_y"],
             scen_out["optimistic"]["procurement_replaced_musd_y"]))
    print("Pd честно: 1000 скидов в мире = %s кг/г замещения = "
          "%s%% выпуска НН"
          % (pd_demand_honest["world_aggressive_1000_skids"][
              "pd_replacement_kg_y"],
             pd_demand_honest["world_aggressive_1000_skids"][
              "share_of_nn_output_pct"]))
    print("\nчувствительность (экономия $/т, база, электролиз):")
    for knob, vals in sens.items():
        if knob == "unit":
            continue
        row = "  ".join("%s:%s" % (k, v) for k, v in vals.items())
        print("  %-44s %s" % (knob, row))
    print("\nkill: сертификация взрывобезопасности | потребность "
          "<300 т/г | delivered <$700/т")
    print("вердикт:", kill["verdict_scenarios"])

    path = os.path.join(DIR, "econ_nornickel_pd_h2o2_results.json")
    with open(path, "w") as f:
        json.dump(out, f, indent=1, ensure_ascii=False)
    print("\nwrote %s" % os.path.basename(path))


if __name__ == "__main__":
    main()
