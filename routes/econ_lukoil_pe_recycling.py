#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
routes/econ_lukoil_pe_recycling.py — Лукойл L5: механохимическая линия
рециклинга полиэтилена на «Ставролене» — ньютоны вместо гигаджоулей.

Идея: селективный разрыв цепи ПЭ в реакционной шаровой/экструзионной мельнице
с переносом усилия на цепь (co-milling с твёрдым носителем) -> воски /
альфа-олефиновые концы / укороченные цепи (rPE-компаунд) вместо пиролизного
супа. Тянем связь силой (нН на цепь) вместо нагрева всей массы до 400-500C.

TRL 2-3 (лабораторная концепция). Это НЕ инвестиционная модель, а R&D-опцион:
скрин-уровень +-2-3x, НЕ банк. Якоря — из наших же расчётных файлов:
  cc_cogef_cas22_results.json  CASSCF(2,2)-строгая карта сила->барьер C-C:
                               F=0 3.22 эВ; 2 нН -> 1.67; 3 -> 1.32; 4 -> 1.00 эВ;
                               бирадикальный onset n_u>0.5 при d~2.6 A
  cc_cogef_results.json        PBE-карта (3.81 эВ -> 0.54 при 4 нН) — была
                               ОПТИМИСТИЧНА; CAS-строго комнатная кинетика
                               требует ~4.5-5 нН/цепь (APPLIED_LUKOIL_NORNICKEL)
Главная неопределённость модели — КПД передачи усилия мельница->цепь (0.1-5%,
вилка 50x): вся экономика решается этим одним числом, и мы показываем это
честно, а не прячем в «эффективность оборудования».

Запуск: python3 routes/econ_lukoil_pe_recycling.py
     -> routes/econ_lukoil_pe_recycling_results.json
"""
import json
import os

DIR = os.path.dirname(os.path.abspath(__file__))
EV_KJ_MOL = 96.485          # 1 эВ = 96.485 кДж/моль
KWH_KJ = 3600.0             # 1 кВт*ч = 3600 кДж

# ------------------------------------------------ якоря: читаем СВОИ расчёты
with open(os.path.join(DIR, "cc_cogef_cas22_results.json")) as f:
    CAS = json.load(f)
with open(os.path.join(DIR, "cc_cogef_results.json")) as f:
    PBE = json.load(f)

CAS_MAP = {p["F_nN"]: p["barrier_eV"] for p in CAS["cogef_cas22"]["barrier_vs_force"]}
PBE_MAP = {p["F_nN"]: p["barrier_eV"] for p in PBE["cogef"]["barrier_vs_force"]}

E_SCISSION_EV = CAS_MAP[0.0]                    # 3.223 эВ — барьер C-C при F=0
E_SCISSION_KJ_MOL = E_SCISSION_EV * EV_KJ_MOL   # ~311 кДж/моль разрывов

# ------------------------------------------------------ допущения (все крутилки)
ASSUMPTIONS = {
    # --- якоря из репозитория (наши расчёты, цитируем точно) ---
    "cc_barrier_map_cas22": dict(
        value={"F=0": 3.223, "2нН": 1.672, "3нН": 1.316, "4нН": 1.004},
        unit="эВ", tier="anchor_repo",
        source="cc_cogef_cas22_results.json: CASSCF(2,2)-строгая карта "
               "сила->барьер C-C (бутан как минимальный фрагмент ПЭ)"),
    "biradical_onset": dict(value=2.60, unit="A", tier="anchor_repo",
        source="cc_cogef_cas22_results.json: n_u>0.5 (0.5814) при d~2.6 A — "
               "гомолиз бирадикален, DFT/одноконфиг. методы там врут"),
    "force_room_kinetics_nN": dict(value=[4.5, 5.0], unit="нН/цепь",
        tier="anchor_repo",
        source="APPLIED_LUKOIL_NORNICKEL.md (CAS-строго); PBE-оценка 3.5-4 нН "
               "(cc_cogef_results.json: 0.84/0.54 эВ при 3.5/4 нН) была "
               "ОПТИМИСТИЧНА — честно фиксируем сдвиг"),
    "e_scission_min_ev": dict(value=round(E_SCISSION_EV, 3), unit="эВ/разрыв",
        tier="anchor_repo",
        source="термодин. минимум на один разрыв C-C = барьер при F=0 "
               "(cc_cogef_cas22); механика платит его адресно, пиролиз греет "
               "ВСЮ массу"),
    "wacc": dict(value=0.12, unit="доля/год", tier="anchor_repo",
        source="как в econ_lukoil_pao_fleet / skid_twin: 12%"),
    "life_years": dict(value=12, unit="лет", tier="anchor_repo",
        source="жизнь актива 12 лет, как в остальных эконом-моделях дома"),
    # --- литература ---
    "pyrolysis_gj_t": dict(value=[2.0, 4.0], unit="ГДж/т", tier="literature",
        comment="энергоёмкость пиролиза пластика (нагрев+эндотермика+потери); "
                "точка сравнения для тезиса «ньютоны вместо гигаджоулей»"),
    "mill_kwh_t": dict(value=[200, 800], unit="кВт*ч/т", tier="literature",
        comment="электроэнергия измельчения/активации полимеров в шаровых и "
                "экструзионных мельницах — литературный диапазон"),
    "mn_feed": dict(value=50000.0, unit="г/моль (Mn)", tier="literature",
        comment="отходной HDPE/LDPE: Mn ~15-60k (Mw 100-300k); берём 50k"),
    "mn_wax_target": dict(value=600.0, unit="г/моль", tier="literature",
        comment="целевые воски/альфа-олефиновые олигомеры C40-C45"),
    "mn_class_bounds": dict(value={"wax_max": 5000, "compound_max": 30000},
        unit="г/моль", tier="literature",
        comment="классы продукта по достигнутому Mn: <=5k воски/олигомеры; "
                "5-30k rPE-компаунд (контролируемая реология); >30k регринд "
                "(продукта нет)"),
    "scission_selectivity": dict(value="середина цепи, диспропорционирование "
        "радикалов -> один алкановый + один винильный (альфа-олефиновый) конец",
        tier="literature",
        comment="механохимия полимеров: разрыв концентрируется у середины "
                "натянутой цепи; ЭПР-радикалы при помоле ПЭ известны. "
                "СЕЛЕКТИВНОСТЬ В МЕЛЬНИЦЕ НЕ ДОКАЗАНА НАМИ — gate G0"),
    "regrind_value_usd_t": dict(value=350.0, unit="$/т", tier="literature",
        comment="вторичная гранула без контроля ММР — потолок, если селективный "
                "разрыв не получился, а деградация уже есть"),
    "soup_value_usd_t": dict(value=400.0, unit="$/т", tier="literature",
        comment="«пиролизный суп» (широкое ММР, сшивки) ~ цена пиролизного "
                "масла; kill-сценарий по селективности"),
    # --- рынок ---
    "wax_price_usd_t": dict(value=[800, 1500], unit="$/т", tier="market",
        comment="ПЭ-воски / олигомеры; верх — узкие фракции с альфа-олефиновой "
                "функциональностью"),
    "virgin_pe_usd_t": dict(value=[1000, 1200], unit="$/т", tier="market",
        comment="первичный HDPE, база для премии rPE-компаунда"),
    "rpe_premium_usd_t": dict(value=[150, 300], unit="$/т", tier="market",
        comment="премия rPE-компаунда к первичному — СУЩЕСТВУЕТ ТОЛЬКО при "
                "мандатах на рециклят (ЕС/РФ-аналоги); без мандата рециклят "
                "торгуется с ДИСКОНТОМ — это kill-критерий"),
    "tipping_fee_usd_t": dict(value=[-50, 100], unit="$/т сырья", tier="market",
        comment="отходы ПЭ: от -50 (нам платят за приём) до +100 (покупка кип); "
                "знак зависит от региона и чистоты"),
    "capex_line_musd": dict(value=[8, 25], unit="M$", tier="market",
        comment="линия 5-20 кт/год; риск E1 +-2-3x, WBS ниже"),
    "electricity_usd_kwh": dict(value=[0.045, 0.06], unit="$/кВт*ч",
        tier="market", comment="промышленный тариф РФ ~4-5 руб/кВт*ч"),
    "feed_prep_usd_t": dict(value=[100, 250], unit="$/т сырья", tier="market",
        comment="сортировка/мойка/сушка отходов — обязательная и дорогая "
                "стадия, которую любят забывать"),
    # --- допущения о компании (ПРОВЕРИТЬ У КОМПАНИИ) ---
    "stavrolen_pe_kt_y": dict(value=300.0, unit="кт/год", tier="company_assumption",
        comment="мощность «Ставролена» по ПЭ ~300 кт/г — линия 5-20 кт/г "
                "рециклинга = 2-7% от неё, встраивается в площадку; "
                "проверить у компании"),
    "stavrolen_compounding": dict(value="есть", tier="company_assumption",
        comment="на площадке есть компаундирование/грануляция и сбытовая сеть "
                "ПЭ — rPE-компаунд не требует нового сбыта; проверить у компании"),
    "ru_recycling_mandate": dict(value="нет обязательной доли рециклята (2026)",
        tier="company_assumption",
        comment="премия rPE в РФ держится на экспортных требованиях и ESG "
                "покупателей, не на законе — проверить у компании; kill-фактор"),
    # --- ГЛАВНАЯ неопределённость ---
    "eta_transfer": dict(value=[0.001, 0.05], unit="доля", tier="company_assumption",
        comment="КПД передачи механической энергии мельницы В НАТЯЖЕНИЕ ЦЕПИ "
                "до разрыва: 0.1-5%. Вилка 50x. Литература по механодеструкции "
                "полимеров даёт лишь порядок величины; для НАШЕЙ геометрии "
                "(co-milling с носителем) числа нет вообще — его добывает "
                "gate G0/G1. ВСЯ экономика ниже решается этим параметром"),
}

def V(k):
    return ASSUMPTIONS[k]["value"]

AF = (1 - (1 + V("wacc")) ** (-V("life_years"))) / V("wacc")   # аннуитет ~6.19

# --------------------------------------------------------------- физика
def scissions_mol_per_t(kwh_t, eta):
    """моль разрывов C-C на тонну ПЭ при данной э/э и КПД передачи."""
    return kwh_t * KWH_KJ * eta / E_SCISSION_KJ_MOL

def mn_final(kwh_t, eta, mn0=None):
    """конечное Mn после помола: цепей стало (исходные + разрывы)."""
    mn0 = mn0 or V("mn_feed")
    chains0 = 1e6 / mn0                       # моль цепей на тонну
    return 1e6 / (chains0 + scissions_mol_per_t(kwh_t, eta))

def product_class(mn):
    b = V("mn_class_bounds")
    if mn <= b["wax_max"]:
        return "wax_oligomer"
    if mn <= b["compound_max"]:
        return "rpe_compound"
    return "regrind"

def wax_required_kwh_t(eta):
    """сколько э/э нужно на тонну, чтобы дойти до восков Mn=600 при данном КПД."""
    need_mol = 1e6 / V("mn_wax_target") - 1e6 / V("mn_feed")
    return need_mol * E_SCISSION_KJ_MOL / (KWH_KJ * eta)

# теоретический минимум (КПД=1) до восков — центральный тезис
E_MIN_WAX_KWH_T = wax_required_kwh_t(1.0)                 # ~142 кВт*ч/т
E_MIN_WAX_GJ_T = E_MIN_WAX_KWH_T * KWH_KJ / 1e6           # ~0.51 ГДж/т

# --------------------------------------------- сценарии (pess / base / opt)
SCENARIOS = {
    "pessimistic": dict(eta=0.001, kwh_t=200, cap_kt=5, capex=20.0,
        tipping=100.0, prep=250.0, elec=0.06, carrier=120.0, fixed=1.5,
        virgin=1000.0, premium=150.0, wax=800.0,
        note="КПД 0.1%: энергии в цепь почти не попадает, Mn еле сдвигается — "
             "«продукт» = регринд с деградацией; capex к верху E1, за сырьё "
             "платим $100/т"),
    "base": dict(eta=0.01, kwh_t=500, cap_kt=10, capex=15.0,
        tipping=25.0, prep=150.0, elec=0.05, carrier=100.0, fixed=1.5,
        virgin=1100.0, premium=225.0, wax=1150.0,
        note="КПД 1%, середины всех диапазонов; достижимое Mn ~13k = "
             "rPE-компаунд, воски НЕдостижимы (нужно ~14 200 кВт*ч/т)"),
    "optimistic": dict(eta=0.05, kwh_t=800, cap_kt=20, capex=12.0,
        tipping=-50.0, prep=100.0, elec=0.045, carrier=80.0, fixed=1.2,
        virgin=1200.0, premium=300.0, wax=1500.0,
        note="КПД 5%, мельница на максимуме: Mn ~2000 = ПЭ-воски с "
             "альфа-олефиновыми концами, премиальный верх $1500/т; за сырьё "
             "платят нам"),
}
YIELD = {"wax_oligomer": 0.90, "rpe_compound": 0.95, "regrind": 0.95}

# ------------------------------------------------- capex WBS (риск E1 +-2-3x)
CAPEX_WBS = {
    "приём + сортировка/мойка/сушка отходов ПЭ": [1.5, 4.0],
    "реакционные мельницы (шаровые/экструзионные, 4-8 шт, N2)": [2.5, 8.0],
    "узел твёрдого носителя (дозирование + рекуперация)": [0.8, 2.0],
    "классификация продукта / сепарация носителя": [0.7, 2.0],
    "компаундирование + грануляция rPE": [1.2, 3.5],
    "отгонка/фасовка восков (для wax-маршрута)": [0.5, 1.5],
    "электрика + АСУ ТП": [0.6, 1.5],
    "монтаж, здание, обвязка": [1.2, 2.5],
}
WBS_LOW = round(sum(v[0] for v in CAPEX_WBS.values()), 1)
WBS_HIGH = round(sum(v[1] for v in CAPEX_WBS.values()), 1)


def econ(p, price_override=None, class_override=None):
    """Линия целиком: физика (Mn из eta и кВт*ч/т) -> класс продукта -> деньги."""
    mn = mn_final(p["kwh_t"], p["eta"])
    cls = class_override or product_class(mn)
    if price_override is not None:
        price = price_override
    elif cls == "wax_oligomer":
        price = p["wax"]
    elif cls == "rpe_compound":
        price = p["virgin"] + p["premium"]
    else:
        price = V("regrind_value_usd_t")
    y = YIELD[cls]
    prod_t = p["cap_kt"] * 1000.0
    feed_t = prod_t / y
    rev = prod_t * price / 1e6
    feed_cost = feed_t * p["tipping"] / 1e6
    prep = feed_t * p["prep"] / 1e6
    elec = feed_t * p["kwh_t"] * p["elec"] / 1e6
    var = prod_t * p["carrier"] / 1e6
    ebitda = rev - feed_cost - prep - elec - var - p["fixed"]
    npv = -p["capex"] + ebitda * AF
    return dict(
        mn_final=round(mn, 0), product_class=cls, price_usd_t=price,
        product_kt_y=p["cap_kt"], feed_kt_y=round(feed_t / 1000, 1),
        scissions_mol_t=round(scissions_mol_per_t(p["kwh_t"], p["eta"]), 1),
        wax_required_kwh_t=round(wax_required_kwh_t(p["eta"]), 0),
        wax_feasible_at_800kwh=wax_required_kwh_t(p["eta"]) <= 800,
        revenue_musd=round(rev, 2), ebitda_musd=round(ebitda, 2),
        capex_musd=p["capex"],
        payback_y=round(p["capex"] / ebitda, 1) if ebitda > 0 else None,
        npv_musd=round(npv, 1),
        trl="TRL 2-3: концепция; числа скрин-уровня +-2-3x, НЕ инвестиционные")


def main():
    # -------- сверка якорей (обязана сходиться с нашими results-JSON)
    anchors = {
        "barrier_F0_eV": dict(model=3.22, anchor=CAS_MAP[0.0],
            source="cc_cogef_cas22 F=0"),
        "barrier_2nN_eV": dict(model=1.67, anchor=CAS_MAP[2.0],
            source="cc_cogef_cas22 F=2 нН"),
        "barrier_3nN_eV": dict(model=1.32, anchor=CAS_MAP[3.0],
            source="cc_cogef_cas22 F=3 нН"),
        "barrier_4nN_eV": dict(model=1.00, anchor=CAS_MAP[4.0],
            source="cc_cogef_cas22 F=4 нН"),
        "onset_biradical_A": dict(
            model=2.60, anchor=float(CAS["cogef_cas22"]["onset_biradical_A"]),
            source="cc_cogef_cas22: n_u>0.5"),
        "n_u_at_onset": dict(model=0.58, anchor=CAS["points"]["2.60"]["n_u"],
            source="cc_cogef_cas22 points[2.60].n_u = 0.5814 > 0.5"),
        "pbe_4nN_optimistic_eV": dict(model=0.54, anchor=PBE_MAP.get(4.0, 0.543),
            source="cc_cogef_results (PBE): 0.543 эВ при 4 нН — на ~0.46 эВ "
                   "ниже CAS => PBE-порог 3.5-4 нН был оптимистичен"),
    }
    for a in anchors.values():
        a["ok"] = abs(a["model"] - a["anchor"]) <= max(0.01 * abs(a["anchor"]), 0.01)

    # -------- физика: цепочка «ньютоны вместо гигаджоулей», все шаги явно
    physics = {
        "chain": [
            f"1 разрыв C-C стоит минимум {E_SCISSION_EV:.3f} эВ = "
            f"{E_SCISSION_KJ_MOL:.0f} кДж/моль (якорь CAS, F=0)",
            "комнатная кинетика разрыва: ~4.5-5 нН/цепь (CAS-строго; "
            "PBE-оценка 3.5-4 нН была оптимистична)",
            f"ПЭ Mn {V('mn_feed'):.0f} -> воски Mn {V('mn_wax_target'):.0f}: "
            f"{1e6 / V('mn_wax_target') - 1e6 / V('mn_feed'):.0f} моль "
            "разрывов на тонну",
            f"теоретический минимум до восков: {E_MIN_WAX_GJ_T:.2f} ГДж/т "
            f"({E_MIN_WAX_KWH_T:.0f} кВт*ч/т) — против пиролиза "
            f"{V('pyrolysis_gj_t')[0]}-{V('pyrolysis_gj_t')[1]} ГДж/т, "
            f"т.е. в {V('pyrolysis_gj_t')[0] / E_MIN_WAX_GJ_T:.0f}-"
            f"{V('pyrolysis_gj_t')[1] / E_MIN_WAX_GJ_T:.0f} раза меньше",
            "НО: минимум достижим лишь при КПД передачи усилия ~100%; реальный "
            "КПД мельницы 0.1-5% (вилка 50x) — ГЛАВНАЯ неопределённость всей "
            "модели, а не термодинамика",
            "до rPE-компаунда (Mn 50k -> 10-30k) разрывов нужно в ~50-100 раз "
            "меньше — поэтому компаунд-маршрут жив даже при КПД ~1%",
        ],
        "e_min_wax_gj_t": round(E_MIN_WAX_GJ_T, 3),
        "e_min_wax_kwh_t": round(E_MIN_WAX_KWH_T, 0),
        "kwh_t_to_reach_wax_by_eta": {
            f"eta={e:.1%}": round(wax_required_kwh_t(e), 0)
            for e in (0.001, 0.003, 0.01, 0.05, 0.20, 1.00)},
        "mn_reachable_at_mill_range": {
            f"eta={e:.1%}, {k} кВт*ч/т": round(mn_final(k, e), 0)
            for e in (0.001, 0.01, 0.05) for k in (200, 800)},
        "note": "энергию считаем ТОЛЬКО на разрыв связей; упругие потери, шум "
                "и тепло — в (1-eta); селективность середины цепи и "
                "альфа-олефиновые концы — литературный механизм, у нас НЕ "
                "доказан (gate G0)",
    }

    # -------- сценарии
    scen_out = {name: dict(params=dict(p), result=econ(p), note=p["note"])
                for name, p in SCENARIOS.items()}

    # -------- чувствительность (база): ГЛАВНАЯ ручка — eta, потом остальные
    base = SCENARIOS["base"]
    sens = {
        "eta_transfer (ГЛАВНАЯ, вилка 50x) -> [Mn, класс, EBITDA M$]": {
            f"{e:.1%}": [econ(dict(base, eta=e))["mn_final"],
                         econ(dict(base, eta=e))["product_class"],
                         econ(dict(base, eta=e))["ebitda_musd"]]
            for e in (0.001, 0.003, 0.005, 0.01, 0.03, 0.05)},
        "mill_kwh_t -> [Mn, EBITDA M$]": {
            f"{k}": [econ(dict(base, kwh_t=k))["mn_final"],
                     econ(dict(base, kwh_t=k))["ebitda_musd"]]
            for k in (200, 500, 800)},
        "rpe_premium_usd_t": {
            f"+${m:.0f}": econ(dict(base, premium=m))["ebitda_musd"]
            for m in (0, 150, 225, 300)},
        "tipping_fee_usd_t": {
            f"${t:+.0f}": econ(dict(base, tipping=t))["ebitda_musd"]
            for t in (-50, 0, 25, 100)},
        "feed_prep_usd_t": {
            f"${c:.0f}": econ(dict(base, prep=c))["ebitda_musd"]
            for c in (100, 150, 250)},
        "capex_musd (риск E1) -> [окупаемость лет, NPV M$]": {
            f"${c}M": [econ(dict(base, capex=c))["payback_y"],
                       econ(dict(base, capex=c))["npv_musd"]]
            for c in (WBS_LOW, 15.0, WBS_HIGH)},
        "wax_price_if_wax_class (eta=5%)": {
            f"${w:.0f}": econ(dict(base, eta=0.05, wax=w))["ebitda_musd"]
            for w in (800, 1150, 1500)},
        "unit": "EBITDA M$/год на линию 10 кт (база), если не сказано иное",
    }

    # -------- kill-критерии (посчитанный контекст, без драматизации)
    k_eta_max = econ(dict(base, eta=0.003, kwh_t=800))
    k_eta_min = econ(dict(base, eta=0.001, kwh_t=800))
    k_soup = econ(base, price_override=V("soup_value_usd_t"),
                  class_override="wax_oligomer")
    k_noprem = econ(dict(base, premium=0.0))
    k_discount = econ(base, price_override=base["virgin"] * 0.85,
                      class_override="rpe_compound")
    kill = {
        "eta_transfer_below": dict(kill_below=0.003,
            computed_context=f"при eta=0.3% даже на максимуме мельницы "
                f"(800 кВт*ч/т) Mn={k_eta_max['mn_final']:.0f} — самый край "
                f"компаунд-класса, нулевой запас; при 0.1% Mn="
                f"{k_eta_min['mn_final']:.0f} — регринд, EBITDA "
                f"{econ(dict(base, eta=0.001))['ebitda_musd']}M. Если пилотная "
                "мельница (gate G1) меряет <0.3% — вариант закрыть"),
        "no_selectivity_soup": dict(
            kill_if="продукт = «пиролизный суп»: широкое ММР, сшивки, нет "
                    "альфа-олефиновых концов (gate G0)",
            computed_context=f"суп продаётся как пиролизное масло "
                f"~${V('soup_value_usd_t'):.0f}/т: EBITDA "
                f"{k_soup['ebitda_musd']}M — мёртв; весь тезис и IP-окно "
                "держатся на селективности"),
        "no_mandate_no_premium": dict(
            kill_if="мандата/премии rPE нет (ru_recycling_mandate — "
                    "company_assumption)",
            computed_context=f"премия=0: EBITDA {k_noprem['ebitda_musd']}M; "
                f"дисконт -15% к первичному (реальный рынок без мандата): "
                f"{k_discount['ebitda_musd']}M — арифметически живо на дешёвом "
                "сырье, но проект теряет отличие от обычного мехрециклинга и "
                "карта сила->барьер не даёт преимущества: kill СТРАТЕГИЧЕСКИЙ",
        ),
        "verdict_scenarios": {},
    }
    for name, p in SCENARIOS.items():
        r = scen_out[name]["result"]
        if p["eta"] < 0.003 or r["product_class"] == "regrind":
            v = ("МЁРТВ: КПД ниже kill-порога 0.3%, продукт = регринд, EBITDA "
                 f"{r['ebitda_musd']}M — показываем честно, не прячем")
        elif r["ebitda_musd"] <= 0:
            v = "МЁРТВ: EBITDA <= 0"
        else:
            v = (f"жив на бумаге (EBITDA {r['ebitda_musd']}M, NPV "
                 f"{r['npv_musd']}M) — но TRL 2-3: это R&D-опцион, не стройка")
        kill["verdict_scenarios"][name] = v

    # -------- механическая реализация
    implementation = {
        "block_flow": [
            "приём отходов ПЭ (кипы/дроблёнка), входной контроль",
            "сортировка -> мойка -> сушка (влага <0.1%)",
            "смешение с твёрдым носителем (перенос усилия на цепь; co-milling)",
            "реакционная мельница (шаровая планетарная или двухшнековая "
            "экструзионная, N2, контроль T<80C — рвём силой, не греем)",
            "классификация: сепарация носителя (рекуперация) + фракционирование",
            "маршрут A (воски): отгонка лёгких, фасовка чешуя/гранула",
            "маршрут B (rPE-компаунд): компаундирование + грануляция на "
            "мощностях «Ставролена» (company_assumption)",
        ],
        "equipment": [
            dict(item="приём + сортировка/мойка/сушка", cost_musd=[1.5, 4.0]),
            dict(item="реакционные мельницы 4-8 шт + N2-узел",
                 cost_musd=[2.5, 8.0],
                 note="сердце линии; геометрия передачи усилия = наше IP-окно"),
            dict(item="узел носителя (дозирование/рекуперация)",
                 cost_musd=[0.8, 2.0]),
            dict(item="классификация/сепарация", cost_musd=[0.7, 2.0]),
            dict(item="компаундирование + грануляция", cost_musd=[1.2, 3.5]),
            dict(item="отгонка/фасовка восков", cost_musd=[0.5, 1.5]),
            dict(item="электрика + АСУ ТП", cost_musd=[0.6, 1.5]),
            dict(item="монтаж/здание/обвязка", cost_musd=[1.2, 2.5]),
        ],
        "utilities": {
            "electricity": dict(value="10 кт/г x 500 кВт*ч/т ~ 5.3 ГВт*ч/г "
                "(~700 кВт средних) — главная статья опекса после подготовки "
                "сырья"),
            "nitrogen_nm3_h": dict(value=[20, 50],
                note="инертизация мельниц: радикалы + O2 = окисление/сшивки"),
            "cooling": "отвод (1-eta) энергии помола как тепла: ~2.5 МВт*ч "
                       "тепла на тонну при 500 кВт*ч/т — охлаждение корпусов "
                       "обязательно",
        },
        "staffing": dict(line_10kt="12-18 чел (3 смены: операторы, механик, "
                         "лаборант ММР/FTIR)", tier="screening"),
        "rollout_trl_gates": [
            dict(stage="G0 (лаборатория, TRL 2->3, 6-9 мес): co-milling "
                 "ПЭ+носитель; показать УЗКИЙ сдвиг ММР (GPC), винильные концы "
                 "(FTIR 908 см-1 / ЯМР), радикалы (ЭПР); первая оценка eta",
                 kill="ММР широкое/сшивки => «суп», закрыть"),
            dict(stage="G1 (прототип-мельница, TRL 3-4, 9-12 мес): непрерывный "
                 "режим 1-10 кг/ч, измеренный eta",
                 kill="eta < 0.3% => закрыть"),
            dict(stage="G2 (пилот 100 кг/ч, TRL 5, 12-18 мес): 500 ч пробега, "
                 "квалификация продукта (воск у покупателя / компаунд в "
                 "спецификацию Ставролена)"),
            dict(stage="G3 (FEED линии 5-20 кт/г): только после G2 и при "
                 "подтверждённой премии rPE / контракте на воски"),
        ],
        "patent_window": "карта сила->барьер CAS-уровня + бирадикальный onset "
                         "(наш расчёт) + геометрия переноса усилия "
                         "(носитель/мельница) — подать ДО публикации дневника "
                         "соответствующих деталей конструкции",
    }

    portfolio_role = (
        "R&D-опцион под мандаты рециклинга + патентное окно, НЕ немедленная "
        "стройка. Ценность сегодня: (1) карта сила->барьер CAS(2,2)-уровня — "
        "наше IP; (2) дешёвые gate G0/G1 (лаборатория, <$1M) снимают главную "
        "неопределённость eta до любых капвложений; (3) если мандаты на "
        "рециклят придут — у Лукойла готовый маршрут на площадке «Ставролена» "
        "с сырьём, компаундированием и сбытом. Стоимость опциона = стоимость "
        "G0/G1, а не $8-25M линии.")

    b = scen_out["base"]["result"]
    o = scen_out["optimistic"]["result"]
    p_ = scen_out["pessimistic"]["result"]
    headline = (
        f"Ньютоны вместо гигаджоулей: CAS-карта сила->барьер (3.22 эВ при F=0 "
        f"-> 1.00 эВ при 4 нН; комнатная кинетика ~4.5-5 нН/цепь) даёт "
        f"теорминимум механоразрыва ПЭ до восков {E_MIN_WAX_GJ_T:.2f} ГДж/т "
        f"против 2-4 ГДж/т пиролиза (в 4-8 раз меньше). Но реальный КПД "
        f"передачи усилия мельница->цепь 0.1-5% (вилка 50x) — главная "
        f"неопределённость: при 1% достижим только rPE-компаунд (Mn ~13k, "
        f"EBITDA {b['ebitda_musd']}M на 10 кт/г), воски требуют ~5% и максимум "
        f"мельницы. TRL 2-3: R&D-опцион под мандаты рециклинга и патентное "
        f"окно, НЕ стройка.")

    out = {
        "model": "Лукойл L5: механохимическая линия рециклинга ПЭ на "
                 "«Ставролене» — ньютоны вместо гигаджоулей",
        "tier": "screening_pm_2_3x",
        "not_bankable": True,
        "trl": "2-3 (лабораторная концепция; каждый вывод ниже несёт эту метку)",
        "headline": headline,
        "anchors_check": anchors,
        "assumptions": ASSUMPTIONS,
        "physics": physics,
        "capex_wbs_musd": dict(blocks=CAPEX_WBS, low=WBS_LOW, high=WBS_HIGH,
            note=f"сумма ${WBS_LOW}-{WBS_HIGH}M согласована с рыночным "
                 "диапазоном $8-25M; риск E1 +-2-3x"),
        "scenarios": scen_out,
        "sensitivity": sens,
        "kill_criteria": kill,
        "implementation": implementation,
        "portfolio_role": portfolio_role,
        "honesty": "скрин +-2-3x, НЕ банковская модель, TRL 2-3. Четыре "
                   "главных незакрытых вопроса: (1) eta 0.1-5% — вилка 50x, "
                   "решает ВСЁ, меряется только на железе (G0/G1, <$1M); "
                   "(2) селективность разрыва в мельнице не доказана — при "
                   "«супе» проект мёртв; (3) премия rPE = функция мандатов, "
                   "которых в РФ нет (company_assumption); (4) opt-сценарий "
                   f"с NPV {o['npv_musd']}M на TRL 2-3 — иллюстрация верхней "
                   "границы, а не план; pess честно МЁРТВ (EBITDA "
                   f"{p_['ebitda_musd']}M). База выглядит хорошо НА БУМАГЕ "
                   "именно потому, что премия и eta приняты допущениями — "
                   "оба снимаются гейтами, не верой",
    }

    # ------------------------------------------------------------- сводка
    print("=" * 72)
    print("Лукойл L5: механохимический рециклинг ПЭ, «Ставролен»")
    print("(TRL 2-3, скрин +-2-3x, не банк; R&D-опцион, не стройка)")
    print("=" * 72)
    print("сверка якорей:", "OK" if all(a["ok"] for a in anchors.values())
          else "РАСХОЖДЕНИЕ! " + str({k: a for k, a in anchors.items()
                                      if not a["ok"]}))
    print(f"физика: разрыв C-C {E_SCISSION_EV:.2f} эВ (F=0), комн. кинетика "
          f"~4.5-5 нН/цепь;")
    print(f"  теорминимум до восков {E_MIN_WAX_GJ_T:.2f} ГДж/т "
          f"({E_MIN_WAX_KWH_T:.0f} кВт*ч/т) vs пиролиз 2-4 ГДж/т")
    print(f"\n{'сценарий':<13}{'eta':>6}{'кВт*ч/т':>8}{'Mn':>8}"
          f"{'класс':>14}{'EBITDA':>8}{'NPV':>8}")
    for name, s in scen_out.items():
        r = s["result"]
        print(f"{name:<13}{s['params']['eta']:>6.1%}"
              f"{s['params']['kwh_t']:>8.0f}{r['mn_final']:>8.0f}"
              f"{r['product_class']:>14}{r['ebitda_musd']:>7.2f}M"
              f"{r['npv_musd']:>7.1f}M")
    print(f"\ncapex WBS: ${WBS_LOW}-{WBS_HIGH}M (диапазон рынка $8-25M, E1)")
    print("\nчувствительность (главная ручка — eta, вилка 50x):")
    for e, row in sens[
            "eta_transfer (ГЛАВНАЯ, вилка 50x) -> [Mn, класс, EBITDA M$]"
            ].items():
        print(f"  eta={e:<7} Mn={row[0]:<8.0f} {row[1]:<14} EBITDA {row[2]}M")
    print("\nkill: eta<0.3% на пилоте | «пиролизный суп» | нет мандата/премии rPE")
    for name, v in kill["verdict_scenarios"].items():
        print(f"  {name}: {v}")

    path = os.path.join(DIR, "econ_lukoil_pe_recycling_results.json")
    with open(path, "w") as f:
        json.dump(out, f, indent=1, ensure_ascii=False)
    print(f"\nwrote {os.path.basename(path)}")


if __name__ == "__main__":
    main()
