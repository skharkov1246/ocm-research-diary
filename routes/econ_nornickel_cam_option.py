#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
routes/econ_nornickel_cam_option.py — Норникель N8: опцион батарейной
вертикали Ni/Co -> сульфат -> прекурсор (pCAM) -> катодный материал (CAM,
NMC811) с квантовым де-рискингом скрининга составов.

ЭТО НЕ МОДЕЛЬ ЗАВОДА. Это real-options СКРИН (дерево решений):
  [квалификация у клиента, 2-4 года, платим qual-программу]
     -> с вероятностью p_success СТРОИМ (капекс = цена исполнения опциона)
     -> иначе выходим, потеряв только qual-программу.
Ценность считаем ТОЛЬКО по надбавке ценностной лестницы (металл продаётся
в любом случае — базовая цена Ni в опцион НЕ входит).

ДВЕ строки итога, НЕ смешивать:
  (1) NPV-опцион вертикали по сценариям — ценность решения компании;
  (2) value-of-information (VoI) нашей квантовой диагностики — дельта
      NPV от +5-15 п.п. к p_success и -0.5-1 год к квалификации, плюс
      прямая экономия синтез-итераций. Это наш продукт, не их завод.

Наш рычаг (якорь репо): вклад корреляции в вертикальное окисление
Ni3+ -> Ni4+ катодного центра [NiO6] = +0.32 эВ > порога ~0.2 В, который
убивает скрининг напряжения -> DFT+U-скрининг кандидатов 5V/high-Ni
ненадёжен -> наша мультиреференс-диагностика режет список кандидатов до
синтеза. Честно: сегодня это ТОЛЬКО диагностика трендов — спиновая щель
меняет знак с размером CASCI-окна, энергия плато не достигает
(battery_cathode_results.json), абсолютные напряжения = полное
CAS(41e,23o) = 46 кубитов = квантовое железо 2029+.

Скрин-уровень +-2-3x, НЕ банковская модель. Данные Норникеля (тоннаж в
вертикаль, вероятность квалификации, статус батарейной стратегии) —
ТОЛЬКО диапазоны tier='company_assumption', проверить у компании.

Запуск: python3 routes/econ_nornickel_cam_option.py
     -> routes/econ_nornickel_cam_option_results.json
"""
import json
import os

DIR = os.path.dirname(os.path.abspath(__file__))
TABS_JSON = os.path.join(DIR, "..", "assets", "matterforge",
                         "matterforge_tabs_content.json")

# ------------------------------------------------------ допущения (все крутилки)
ASSUMPTIONS = {
    # --- якоря репозитория ---
    "wacc": dict(value=0.12, unit="доля/год", tier="anchor_repo",
        source="skid_twin / econ_shock_candidates / econ_nornickel_ew_oer: "
               "12% — единая ставка дисконтирования дома"),
    "corr_shift_ni3_ni4_ev": dict(value=0.32, unit="эВ", tier="anchor_repo",
        source="вкладка battery-cathode (запись 2026-06-08): вертикальное "
               "окисление Ni3+->Ni4+ в [NiO6]-вложении, CASCI-HF = +0.32 эВ; "
               "абсолют — не напряжение (заряд кластера), значим сдвиг метода"),
    "voltage_kill_threshold_v": dict(value=0.2, unit="В", tier="anchor_repo",
        source="вкладка battery-cathode: ошибка напряжения ~0.2 В убивает "
               "скрининг; 0.32 > 0.2 -> одно-детерминантный DFT+U обязан "
               "систематически врать на катодном центре"),
    "full_cas_qubits": dict(value=46, unit="кубитов", tier="anchor_repo",
        source="вкладка battery-cathode: полное AVAS Ni 3d + O 2p = "
               "CAS(41e,23o) = 46 кубитов (JW); 23 орбитали x 2"),
    "no_energy_plateau": dict(value=True, unit="-", tier="anchor_repo",
        source="battery_cathode_results.json (recovered pipeline): CASCI-окна "
               "8q->24q — каждое расширение опускает E ещё на 0.26-0.81 эВ, "
               "плато-критерий 0.03 эВ не выполнен ни разу -> обрезанные окна "
               "для целевой точности 0.2 В не годятся"),
    "spin_gap_sign_flips": dict(value=True, unit="-", tier="anchor_repo",
        source="battery_cathode_results.json: щель дублет/квартет по окнам "
               "-0.19 -> +1.05 -> -0.15 -> -0.98 -> -0.99 эВ — знак меняется "
               "с окном -> абсолюты сегодня НЕ считаемы, продукт = тренды"),
    # --- рынок: ценностная лестница, $/т Ni-экв ---
    "ni_metal_usd_t": dict(value=[16000, 19000], unit="$/т", tier="market",
        comment="базовая цена Ni-металла — СПРАВОЧНО: в опцион НЕ входит, "
                "металл продаётся в любом случае"),
    "sulfate_premium_usd_t": dict(value=[1000, 2000], unit="$/т Ni",
        tier="market", comment="премия сульфата Ni battery-grade над металлом"),
    "pcam_add_usd_t": dict(value=[2000, 4000], unit="$/т Ni-экв",
        tier="market", comment="добавленная стоимость прекурсора pCAM"),
    "cam_add_usd_t": dict(value=[4000, 8000], unit="$/т Ni-экв",
        tier="market", comment="добавленная стоимость CAM (NMC811) сверх pCAM"),
    # --- литература ---
    "cam_capex_30kt_busd": dict(value=[0.8, 1.5], unit="$B за завод 30 кт "
        "Ni-экв/год", tier="literature",
        comment="интегрированный pCAM+CAM; цена исполнения опциона"),
    "capex_scale_exp": dict(value=0.7, unit="-", tier="literature",
        comment="степенной закон масштабирования капекса (six-tenths rule, "
                "консервативно 0.7)"),
    "project_life_y": dict(value=20, unit="лет", tier="literature",
        comment="срок жизни завода после пуска"),
    "build_y": dict(value=2, unit="лет", tier="literature",
        comment="стройка после положительного решения (исполнения опциона)"),
    "qual_program_musd": dict(value=[30, 80], unit="$M", tier="literature",
        comment="пилотная линия pCAM/CAM + программа квалификации у 2-3 "
                "клиентов (образцы, аудиты, циклирование) — премия опциона"),
    "ebitda_capture_share": dict(value=[0.30, 0.60], unit="доля надбавки "
        "лестницы, остающаяся как EBITDA", tier="literature",
        comment="конверсионные opex (реагенты, энергия, персонал, брак) "
                "съедают 40-70% надбавки; верх диапазона — интегрированный "
                "игрок с captive-металлом и дешёвой энергией; проверить"),
    "iter_cost_musd": dict(value=[2, 10], unit="$M/итерация", tier="literature",
        comment="одна синтез-итерация состава CAM: синтез партии, ячейки, "
                "циклирование, аналитика; длительность 6-12 мес — временнАя "
                "часть учтена рычагом dt_y, здесь только прямые затраты "
                "(без двойного счёта)"),
    # --- допущения о компании (ПРОВЕРИТЬ У КОМПАНИИ) ---
    "vertical_kt_y": dict(value=[10, 30, 60], unit="кт Ni-экв/год в вертикаль",
        tier="company_assumption",
        comment="сценарный тоннаж (pess/base/opt) — доля никелевого объёма, "
                "направляемая в батарейную вертикаль; проверить у компании"),
    "p_success": dict(value=[0.30, 0.60], unit="вероятность",
        tier="company_assumption",
        comment="вероятность технологического успеха квалификации CAM у "
                "клиента (химия + процесс + аудит + контракт); проверить"),
    "t_qual_y": dict(value=[2, 4], unit="лет", tier="company_assumption",
        comment="время квалификации у клиента до решения строить"),
    "battery_strategy_alive": dict(value=None, unit="да/нет",
        tier="company_assumption",
        comment="СТАТУС батарейной стратегии компании — бинарный гейт всей "
                "модели, НЕ считаем, спрашиваем; если закрыта -> kill"),
    "iters_saved": dict(value=[1, 3], unit="итераций",
        tier="company_assumption",
        comment="сколько синтез-итераций режет наша диагностика (короче "
                "список кандидатов до синтеза); проверить на пилоте G2"),
    "dp_success_pp": dict(value=[5, 15], unit="п.п.",
        tier="company_assumption",
        comment="прирост p_success от де-рискинга выбора состава; верх (15) "
                "требует квантового железа (абсолюты, 46q, 2029+) — см. "
                "quantum_hw_dependency"),
    "dt_qual_y": dict(value=[0.5, 1.0], unit="лет",
        tier="company_assumption",
        comment="сокращение квалификации: меньше синтез-циклов по 6-12 мес"),
}


def V(k):
    return ASSUMPTIONS[k]["value"]


def pick(rng, w):
    """Точка в диапазоне [lo, hi]: w=0 низ, 0.5 середина, 1 верх."""
    return rng[0] + (rng[1] - rng[0]) * w


def annuity(life_y):
    """PV-фактор аннуитета по WACC."""
    w = V("wacc")
    return (1.0 - (1.0 + w) ** (-life_y)) / w


def df(t_y):
    """Дисконт-фактор к году t."""
    return (1.0 + V("wacc")) ** (-t_y)


# --------------------------------------------- сценарии (pess / base / opt)
# Для ОПЦИОНА пессимизм = малый тоннаж, низ лестницы и захвата маржи,
# дорогой завод, низкая вероятность и долгая квалификация. capex_w: 1 = верх
# диапазона капекса (пессимистично). VoI-рычаги (dp/dt/итерации) движутся
# согласованно со сценарием.
SCENARIOS = {
    "pessimistic": dict(kt_y=10.0, ladder_w=0.0, capture_w=0.0, capex_w=1.0,
        p_success=0.30, t_qual_y=4.0, qual_musd=80.0,
        dp_pp=5.0, dt_y=0.5, iters_saved=1.0, iter_cost_musd=2.0,
        note="10 кт, низ лестницы и захвата, капекс по верху, p=30%, "
             "квалификация 4 года"),
    "base": dict(kt_y=30.0, ladder_w=0.5, capture_w=0.5, capex_w=0.5,
        p_success=0.45, t_qual_y=3.0, qual_musd=55.0,
        dp_pp=10.0, dt_y=0.75, iters_saved=2.0, iter_cost_musd=6.0,
        note="30 кт, середины всех диапазонов, p=45%, квалификация 3 года"),
    "optimistic": dict(kt_y=60.0, ladder_w=1.0, capture_w=1.0, capex_w=0.0,
        p_success=0.60, t_qual_y=2.0, qual_musd=30.0,
        dp_pp=15.0, dt_y=1.0, iters_saved=3.0, iter_cost_musd=10.0,
        note="60 кт, верх лестницы и захвата, капекс по низу, p=60%, "
             "квалификация 2 года"),
}

# Диагностический (без квантового железа) угол VoI — фиксирован по низу:
# без валидированных абсолютов защитим только +5 п.п. и -0.5 года.
DIAG_ONLY = dict(dp_pp=5.0, dt_y=0.5)


def ladder_premium_usd_t(w):
    """Надбавка полной лестницы металл->CAM, $/т Ni-экв (без цены металла)."""
    return (pick(V("sulfate_premium_usd_t"), w) +
            pick(V("pcam_add_usd_t"), w) + pick(V("cam_add_usd_t"), w))


def project_block(kt_y, ladder_w, capture_w, capex_w):
    """Завод (цена исполнения + поток) НА МОМЕНТ РЕШЕНИЯ строить.
    НЕ модель завода: одна EBITDA-строка по захвату надбавки лестницы."""
    prem = ladder_premium_usd_t(ladder_w)
    capture = pick(V("ebitda_capture_share"), capture_w)
    ebitda = kt_y * 1000.0 * prem * capture / 1e6            # $M/год
    capex30 = pick(V("cam_capex_30kt_busd"), capex_w) * 1000.0
    capex = capex30 * (kt_y / 30.0) ** V("capex_scale_exp")
    pv_in = ebitda * annuity(V("project_life_y")) * df(V("build_y"))
    npv_dec = pv_in - capex
    # breakeven-захват: какая доля надбавки должна остаться EBITDA,
    # чтобы завод отбился при этом капексе (вопрос компании #1)
    denom = kt_y * 1000.0 * prem / 1e6 * annuity(
        V("project_life_y")) * df(V("build_y"))
    be_capture = capex / denom if denom > 0 else None
    return dict(ladder_premium_usd_t=round(prem, 0),
                ebitda_capture=round(capture, 2),
                ebitda_musd_y=round(ebitda, 1),
                capex_musd=round(capex, 0),
                pv_inflows_musd=round(pv_in, 0),
                npv_at_decision_musd=round(npv_dec, 0),
                breakeven_capture_share=round(be_capture, 2)
                if be_capture else None)


def option_value(npv_dec, p, t_qual, qual_musd):
    """Опцион (дерево решений, не Блэк-Шоулз): -PV(квалификация) +
    p * max(0, NPV_решения) * дисконт. max(0,.) — право НЕ строить."""
    qual_pv = qual_musd * df(t_qual / 2.0)      # траты размазаны по периоду
    exercise = max(0.0, npv_dec)
    raw = -qual_pv + p * exercise * df(t_qual)
    return dict(qual_pv_musd=round(qual_pv, 1),
                exercise_musd=round(exercise, 0),
                option_npv_musd=round(raw, 1),
                # право не входить в квалификацию вовсе:
                option_floor0_musd=round(max(0.0, raw), 1),
                expected_npv_no_abandon_musd=round(
                    -qual_pv + p * npv_dec * df(t_qual), 1))


def voi_block(proj, p, t_qual, qual_musd, dp_pp, dt_y, iters_saved,
              iter_cost_musd):
    """Value-of-information нашей диагностики, УСЛОВНО на том, что компания
    вертикаль ведёт (бинарный гейт стратегии — в kill_criteria).
    Три слагаемых, без двойного счёта: (1) прямая экономия итераций,
    (2) дельта от +dp п.п. к вероятности, (3) дельта от -dt лет."""
    npv_dec = proj["npv_at_decision_musd"]
    base_opt = option_value(npv_dec, p, t_qual, qual_musd)
    rnd = iters_saved * iter_cost_musd * df(t_qual / 2.0)
    p_term = (dp_pp / 100.0) * max(0.0, npv_dec) * df(t_qual)
    t_opt = option_value(npv_dec, p, t_qual - dt_y, qual_musd)
    t_term = t_opt["option_npv_musd"] - base_opt["option_npv_musd"]
    total = rnd + p_term + t_term
    return dict(rnd_iterations_musd=round(rnd, 1),
                p_uplift_musd=round(p_term, 1),
                time_shortening_musd=round(t_term, 1),
                voi_total_musd=round(total, 1),
                note="p-рычаг монетизируется ТОЛЬКО если опцион в деньгах "
                     "(max(0, NPV_решения)); при опционе вне денег ускорение "
                     "лишь приближает затраты (t-слагаемое может быть <0) — "
                     "не прячем")


def scenario(kt_y, ladder_w, capture_w, capex_w, p_success, t_qual_y,
             qual_musd, dp_pp, dt_y, iters_saved, iter_cost_musd, **_):
    """Полный сценарий: завод-на-решении + опцион + VoI (полный и
    диагностический без квантового железа)."""
    proj = project_block(kt_y, ladder_w, capture_w, capex_w)
    opt = option_value(proj["npv_at_decision_musd"], p_success, t_qual_y,
                       qual_musd)
    voi = voi_block(proj, p_success, t_qual_y, qual_musd, dp_pp, dt_y,
                    iters_saved, iter_cost_musd)
    voi_diag = voi_block(proj, p_success, t_qual_y, qual_musd,
                         DIAG_ONLY["dp_pp"], DIAG_ONLY["dt_y"],
                         iters_saved, iter_cost_musd)
    return dict(project_at_decision=proj, option=opt,
                voi_full=voi, voi_diagnostic_only=voi_diag)


def main():
    # -------- сверка с якорями (обязана сходиться, иначе модель врёт)
    bat = json.load(open(os.path.join(DIR, "battery_cathode_results.json")))
    tabs = json.load(open(TABS_JSON))
    tab = next(t for t in tabs["tabs"] if t["id"] == "battery-cathode")
    tab_text = " ".join(u["text"] for u in tab["updates"])
    gaps = bat["analysis"]["spin_gap_eV"]["gaps"]
    signs = set(g > 0 for g in gaps.values())
    anchors = {
        "corr_shift_in_tab": dict(
            model="+0.32 эВ" in tab_text, anchor=True,
            source="matterforge_tabs_content.json / battery-cathode: вклад "
                   "корреляции Ni3+->Ni4+ = +0.32 эВ (запись 2026-06-08)"),
        "corr_above_threshold": dict(
            model=V("corr_shift_ni3_ni4_ev") > V("voltage_kill_threshold_v"),
            anchor=True, source="0.32 эВ > 0.2 В — предпосылка ненадёжности "
                                "DFT+U-скрининга"),
        "qubits_46_in_tab": dict(
            model=("46 кубитов" in tab_text) and (2 * 23 == 46), anchor=True,
            source="CAS(41e,23o) -> 46 кубитов; строка есть в записях трека"),
        "spin_gap_sign_flips": dict(
            model=(len(signs) == 2), anchor=True,
            source="battery_cathode_results.json: analysis.spin_gap_eV.gaps "
                   "содержит и + и - -> знак щели зависит от окна"),
        "no_energy_plateau": dict(
            model=(bat["analysis"]["spin1"]["plateau_qubits"] is None and
                   bat["analysis"]["spin3"]["plateau_qubits"] is None),
            anchor=True,
            source="battery_cathode_results.json: плато не достигнуто ни в "
                   "одном спине -> обрезанные окна не годятся для 0.2 В"),
        "annuity_20y_12pct": dict(
            model=round(annuity(20), 3), anchor=7.469,
            source="контроль арифметики: PV-фактор 20 лет @ 12%"),
    }
    for a in anchors.values():
        if isinstance(a["model"], bool):
            a["ok"] = (a["model"] == a["anchor"])
        else:
            a["ok"] = abs(a["model"] - a["anchor"]) <= 0.005 * a["anchor"]

    scen_out = {name: dict(params=dict(p), **scenario(**p))
                for name, p in SCENARIOS.items()}
    base = SCENARIOS["base"]
    sb = scen_out["base"]

    # -------- чувствительность (вокруг базы): раздельно опцион и VoI
    def tot(**over):
        return scenario(**dict(base, **over))

    def pair(s):
        # тройка [NPV решения; опцион; VoI]: первая колонка показывает,
        # ПОЧЕМУ опцион/VoI не двигаются — база вне денег, max(0,.) режет
        return [s["project_at_decision"]["npv_at_decision_musd"],
                s["option"]["option_npv_musd"],
                s["voi_full"]["voi_total_musd"]]

    sens = {
        "p_success": {"%.0f%%" % (p * 100): pair(tot(p_success=p))
                      for p in (0.30, 0.45, 0.60)},
        "capture_share (низ/середина/верх 0.30-0.60)": {
            n: pair(tot(capture_w=w))
            for n, w in (("низ", 0.0), ("середина", 0.5), ("верх", 1.0))},
        "ladder (низ/середина/верх $7-14k/т)": {
            n: pair(tot(ladder_w=w))
            for n, w in (("низ", 0.0), ("середина", 0.5), ("верх", 1.0))},
        "kt_y (тоннаж в вертикаль)": {
            "%d кт" % k: pair(tot(kt_y=float(k))) for k in (10, 30, 60)},
        "capex_30kt (низ/середина/верх $0.8-1.5B)": {
            n: pair(tot(capex_w=w))
            for n, w in (("низ", 0.0), ("середина", 0.5), ("верх", 1.0))},
        "t_qual_y": {"%d года" % t: pair(tot(t_qual_y=float(t)))
                     for t in (2, 3, 4)},
        "dp_pp (только VoI)": {"+%d п.п." % d:
                               tot(dp_pp=float(d))["voi_full"]
                               ["voi_total_musd"] for d in (5, 10, 15)},
        "iters_saved x iter_cost (только VoI)": {
            "%dx$%dM" % (n, c): tot(iters_saved=float(n),
                                    iter_cost_musd=float(c))["voi_full"]
            ["voi_total_musd"]
            for n, c in ((1, 2), (2, 6), (3, 10))},
        "unit": "[NPV решения $M ; NPV-опцион $M ; VoI $M] — опцион и VoI — "
                "ДВЕ РАЗНЫЕ строки, не суммировать; одинаковые значения по "
                "ручке = опцион вне денег (max(0,.) режет рычаг) — это "
                "информация, не баг",
    }

    # -------- kill-критерии
    p20 = {n: option_value(
        scen_out[n]["project_at_decision"]["npv_at_decision_musd"],
        0.20, SCENARIOS[n]["t_qual_y"], SCENARIOS[n]["qual_musd"])
        ["option_npv_musd"] for n in scen_out}
    kill = {
        "battery_strategy_closed": dict(kill_if=True,
            computed_context="БИНАРНЫЙ гейт, не считается: если батарейная "
                "стратегия компании закрыта — весь опцион и весь VoI = 0; "
                "первый вопрос компании, до любых цифр"),
        "cam_market_contracted_p_below": dict(kill_below=0.20,
            computed_context="рынок CAM законтрактован западными/китайскими "
                "игроками -> барьер квалификации: при p=0.20 опцион $M = %s "
                "(base остаётся вне денег, optimistic мельчает, но жив) — "
                "если реальный доступ к клиентам ниже, входить не во что"
                % json.dumps({k: round(v, 0) for k, v in p20.items()},
                             ensure_ascii=False)),
        "voi_below_musd": dict(kill_below=5.0,
            computed_context="VoI по сценариям: pess $%.1fM (НИЖЕ порога — "
                "в пессимистичном мире наш продукт для этого трека не "
                "продаётся), base $%.1fM (проходит), opt $%.1fM"
                % (scen_out["pessimistic"]["voi_full"]["voi_total_musd"],
                   sb["voi_full"]["voi_total_musd"],
                   scen_out["optimistic"]["voi_full"]["voi_total_musd"])),
        "capture_above_breakeven_needed": dict(
            kill_above_share=sb["project_at_decision"]
            ["breakeven_capture_share"],
            computed_context="производный ИЗ модели: в базе завод отбивается "
                "только при захвате >= %.0f%% надбавки лестницы — ВЫШЕ верха "
                "литературного диапазона (60%%); т.е. базовый угол вне денег "
                "не из-за вероятности, а из-за конверсионной маржи"
                % (sb["project_at_decision"]["breakeven_capture_share"]
                   * 100.0)),
        "verdict_scenarios": {},
    }
    for name, s in scen_out.items():
        opt_alive = s["option"]["option_npv_musd"] > 0
        voi_alive = s["voi_full"]["voi_total_musd"] >= 5.0
        if opt_alive and voi_alive:
            v = ("опцион В ДЕНЬГАХ ($%.0fM) и VoI $%.0fM — и вертикаль, и "
                 "наш продукт живы"
                 % (s["option"]["option_npv_musd"],
                    s["voi_full"]["voi_total_musd"]))
        elif voi_alive:
            v = ("опцион ВНЕ денег ($%.0fM), но VoI $%.1fM >= $5M — наш "
                 "продукт жив как срез R&D-затрат квалификации"
                 % (s["option"]["option_npv_musd"],
                    s["voi_full"]["voi_total_musd"]))
        else:
            v = "МЁРТВ по обеим строкам (опцион вне денег, VoI < $5M)"
        kill["verdict_scenarios"][name] = v

    # -------- зависимость от квантового железа (честное поле)
    quantum_hw_dependency = dict(
        today="диагностика трендов: мультиреференс-метрика n_u насыщается к "
              "12-16 кубитам (считаем классически) -> флаг «DFT+U здесь "
              "ненадёжен» + ранжирование кандидатов; этого хватает на резку "
              "списка до синтеза (rnd-слагаемое VoI) и на нижний край "
              "dp=+5 п.п.",
        needs_hardware="абсолютные напряжения / спин-порядок: энергия плато "
              "не достигает, спиновая щель меняет знак с окном -> полное "
              "CAS(41e,23o) = 46 кубитов, горизонт 2029+ (таймлайн трека); "
              "верх dp=+15 п.п. защитим только с железом",
        consequence="без квантового железа VoI ограничен строкой "
              "voi_diagnostic_only: base $%.1fM vs полный $%.1fM; "
              "optimistic $%.0fM vs $%.0fM — дельта и есть цена железа "
              "для этого трека"
              % (sb["voi_diagnostic_only"]["voi_total_musd"],
                 sb["voi_full"]["voi_total_musd"],
                 scen_out["optimistic"]["voi_diagnostic_only"]
                 ["voi_total_musd"],
                 scen_out["optimistic"]["voi_full"]["voi_total_musd"]),
        honest_artifact="в базе диагностическая строка чуть ВЫШЕ полной — "
              "не опечатка: опцион вне денег, t-слагаемое отрицательно "
              "(ускорение приближает затраты), и больший dt полного рычага "
              "вредит сильнее; dp/dt монетизируются только в деньгах — "
              "см. optimistic")

    # -------- механическая реализация (НАШЕГО продукта — не завода)
    implementation = {
        "scope_note": "реализация ниже — про НАШУ диагностику (софт + "
                      "расчёты). Завод pCAM/CAM — цена исполнения опциона "
                      "компании, его FEED мы НЕ делаем",
        "block_flow": [
            "вход: лонг-лист составов клиента (high-Ni NMC / 5V LNMO, "
            "допанты, покрытия)",
            "кластер редокс-центра [MO6] + Madelung-вложение (конвейер "
            "battery_cathode.py, воспроизведён после потери контейнера)",
            "мультиреференс-диагностика: NOON/n_u на CASCI-окнах 12-16q "
            "(насыщается — считается классически, минуты-часы)",
            "флаг ненадёжности DFT+U по составу: вклад корреляции в "
            "Ni3+->Ni4+ vs порог 0.2 В (якорь +0.32 эВ)",
            "выход: шорт-лист в синтез (режем 1-3 итерации по $2-10M) + "
            "карта «где DFT+U можно верить, где нет»",
            "2029+: те же кластеры на квантовом железе, 46 кубитов -> "
            "абсолютные напряжения и спин-порядок (verify-этап)",
        ],
        "equipment": [
            dict(item="скрининговая кампания: CPU-инстансы (r7i/128ГБ, "
                      "AWS), 20-50 составов", cost_musd=[0.05, 0.3]),
            dict(item="валидационный сет: ретро-расчёт известных катодов "
                      "(LCO, NMC532/811) для калибровки трендов",
                 cost_musd=[0.05, 0.15]),
            dict(item="интеграция в скрининг-цикл клиента (API/отчёты)",
                 cost_musd=[0.1, 0.3]),
        ],
        "footprint_m2": [0, 0],
        "footprint_note": "софт и облако; физического следа нет — завод не "
                          "наш периметр",
        "staffing": dict(pilot="2-3 чел: квантовый химик, инженер конвейера, "
                               "связка с R&D клиента",
                         tier="screening"),
        "rollout": [
            dict(stage="G1 ретро-валидация: тренды напряжений известных "
                       "катодов ранжируются правильно", duration="3-6 мес"),
            dict(stage="G2 слепой тест: наши флаги ДО синтеза клиента, "
                       "сверка после — подтверждаем iters_saved>=1",
                 duration="6-12 мес"),
            dict(stage="G3 интеграция в квалификационный цикл (это и даёт "
                       "dp/dt в VoI)", duration="6-12 мес"),
            dict(stage="G4 (2029+, квантовое железо): 46q-абсолюты, "
                       "верх dp=+15 п.п.", duration="по готовности железа"),
        ],
        "trl_gates": [
            "G1: корреляция наших трендов с известными напряжениями; провал "
            "-> метод не продаётся, kill трека как продукта",
            "G2: >=1 подтверждённая сэкономленная итерация в слепом тесте — "
            "иначе rnd-слагаемое VoI обнуляется",
            "G3: клиент подтверждает сокращение цикла >=0.5 года",
            "G4: 46q-расчёт воспроизводит измеренное напряжение +-0.1 В",
        ],
    }

    # -------- ДВЕ строки итога (не смешивать)
    headline = {
        "line_1_option_npv_musd": {
            n: s["option"]["option_npv_musd"] for n, s in scen_out.items()},
        "line_1_note": "NPV-опцион батарейной вертикали (решение КОМПАНИИ): "
                       "премия = qual-программа, исполнение = капекс завода, "
                       "выплата = захват надбавки лестницы",
        "line_2_voi_musd": {
            n: s["voi_full"]["voi_total_musd"] for n, s in scen_out.items()},
        "line_2_diagnostic_only_musd": {
            n: s["voi_diagnostic_only"]["voi_total_musd"]
            for n, s in scen_out.items()},
        "line_2_note": "value-of-information НАШЕЙ диагностики (наш продукт): "
                       "экономия итераций + dp к вероятности + dt ко времени; "
                       "условно на живой стратегии компании",
        "do_not_mix": "строки НЕ складывать: (1) — ценность чужого решения, "
                      "(2) — ценность нашей информации для этого решения",
    }

    out = {
        "model": "Норникель N8: опцион батарейной вертикали Ni/Co -> pCAM -> "
                 "CAM (NMC811) с квантовым де-рискингом скрининга",
        "frame": "ЭТО НЕ МОДЕЛЬ ЗАВОДА — это модель ОПЦИОНА: real-options "
                 "скрин (дерево решений: квалификация -> p_success -> "
                 "строить/выйти). Никакого проектирования катода сегодня не "
                 "происходит и не заявляется; считается ценность ПРАВА войти "
                 "в вертикаль и ценность ИНФОРМАЦИИ, повышающей качество "
                 "этого решения. Волатильность не моделируем (не "
                 "Блэк-Шоулз) — консервативно недооцениваем опцион.",
        "tier": "screening_pm_2_3x",
        "not_bankable": True,
        "tier_legend": {
            "anchor_repo": "число/факт из наших файлов (battery_cathode_"
                           "results.json, вкладка battery-cathode)",
            "market": "рыночный диапазон, порядок величины",
            "company_assumption": "допущение о компании — проверить у неё",
            "literature": "литература/индустриальная практика, не наш расчёт"},
        "anchors_check": anchors,
        "assumptions": ASSUMPTIONS,
        "scenarios": scen_out,
        "headline": headline,
        "sensitivity": sens,
        "kill_criteria": kill,
        "quantum_hw_dependency": quantum_hw_dependency,
        "implementation": implementation,
        "honesty": "скрин +-2-3x. Главные честные выводы: (1) базовый угол "
                   "опциона ВНЕ денег — при середине лестницы и захвате 45%% "
                   "завод не отбивает капекс при WACC 12%% (нужен захват "
                   ">=%.0f%% — выше литературного верха); опцион в деньгах "
                   "только в оптимистичном углу ($%.0fM), т.е. вертикаль — "
                   "ставка на верх лестницы+захвата, не на середину. "
                   "(2) VoI нашей диагностики при этом ПОЛОЖИТЕЛЕН и в базе "
                   "($%.1fM > kill $5M) — потому что экономия синтез-итераций "
                   "не зависит от того, в деньгах ли опцион; но большие "
                   "цифры VoI ($%.0fM) живут только там, где опцион в "
                   "деньгах. (3) Сегодня продаётся ТОЛЬКО диагностика "
                   "трендов: спиновая щель меняет знак с CASCI-окном, плато "
                   "энергии нет — абсолюты = 46 кубитов = 2029+; строка "
                   "voi_diagnostic_only — наш защитимый минимум. (4) Все "
                   "ключевые входы p_success/тоннаж/статус стратегии — "
                   "company_assumption, проверить у компании; лестница — "
                   "рыночные диапазоны, Ni-экв (не тонны CAM). Числа "
                   "детерминированы, Монте-Карло нет."
                   % (sb["project_at_decision"]["breakeven_capture_share"]
                      * 100.0,
                      scen_out["optimistic"]["option"]["option_npv_musd"],
                      sb["voi_full"]["voi_total_musd"],
                      scen_out["optimistic"]["voi_full"]["voi_total_musd"]),
    }

    # ------------------------------------------------------------- сводка
    print("=" * 74)
    print("Норникель N8: ОПЦИОН батарейной вертикали + VoI квантового "
          "де-рискинга")
    print("(скрин +-2-3x, не банк; модель опциона, НЕ завода)")
    print("=" * 74)
    print("сверка якорей:", "OK" if all(a["ok"] for a in anchors.values())
          else "РАСХОЖДЕНИЕ! " + str(anchors))
    print("рычаг: вклад корреляции Ni3+->Ni4+ = +0.32 эВ > порога 0.2 В; "
          "абсолюты = CAS(41e,23o) = 46 кубитов (2029+)")
    print("\n%-13s%6s%10s%11s%11s%12s%10s" % ("сценарий", "кт/г",
          "EBITDA$M", "капекс$M", "NPVреш$M", "ОПЦИОН $M", "VoI $M"))
    for name, s in scen_out.items():
        pr, op = s["project_at_decision"], s["option"]
        print("%-13s%6.0f%10.1f%11.0f%11.0f%12.1f%10.1f"
              % (name, SCENARIOS[name]["kt_y"], pr["ebitda_musd_y"],
                 pr["capex_musd"], pr["npv_at_decision_musd"],
                 op["option_npv_musd"], s["voi_full"]["voi_total_musd"]))
    print("\nСТРОКА 1 (опцион вертикали, $M):",
          headline["line_1_option_npv_musd"])
    print("СТРОКА 2 (VoI диагностики, $M):  ",
          headline["line_2_voi_musd"])
    print("  из них без квантового железа:  ",
          headline["line_2_diagnostic_only_musd"])
    print("  (строки НЕ складывать)")
    print("\nбаза: breakeven-захват надбавки %.0f%% (> верха диапазона 60%%) "
          "-> базовый угол вне денег по МАРЖЕ, не по вероятности"
          % (sb["project_at_decision"]["breakeven_capture_share"] * 100.0))
    print("\nчувствительность [NPV решения $M ; опцион $M ; VoI $M]:")
    for knob, vals in sens.items():
        if knob == "unit":
            continue
        row = "  ".join("%s:%s" % (k, v) for k, v in vals.items())
        print("  %-44s %s" % (knob, row))
    print("\nkill: стратегия закрыта (бинарно) | p<=0.20 (рынок "
          "законтрактован) | VoI<$5M | захват < breakeven %.0f%%"
          % (sb["project_at_decision"]["breakeven_capture_share"] * 100.0))
    print("вердикт:", kill["verdict_scenarios"])

    path = os.path.join(DIR, "econ_nornickel_cam_option_results.json")
    with open(path, "w") as f:
        json.dump(out, f, indent=1, ensure_ascii=False)
    print("\nwrote %s" % os.path.basename(path))


if __name__ == "__main__":
    main()
