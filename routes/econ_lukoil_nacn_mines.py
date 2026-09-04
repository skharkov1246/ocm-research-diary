#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
routes/econ_lukoil_nacn_mines.py — Лукойл L2: NaCN-скиды у золотых ГОКов +
каннибализация Саратоворгсинтеза. Углубление S1 из econ_shock_candidates.py:
 (1) честный масс-баланс Андруссова CH4+NH3+1.5O2 -> HCN -> абсорбция в NaOH
     -> NaCN (30% р-р или брикеты), выход 65% (якорь);
 (2) экономика одного скида — якорь воспроизводится ($823/т variable,
     $1263/т full, маржа $1737/т при delivered $3000/т, EBITDA $4.34M);
 (3) сценарии pessimistic/base/optimistic + чувствительности: цена NH3,
     цена золота (через спрос ГОКа), плечо доставки классики, капекс;
 (4) флот скидов по ГОКам-якорям и ЧИСТЫЙ эффект для Лукойла:
     маржа скидов минус потерянная маржа Саратоворгсинтеза.
Скрин-уровень ±2-3x, НЕ банковская модель. Все цифры по Лукойлу /
Саратоворгсинтезу — company_assumption, проверять у компании. Запуск:
  python3 routes/econ_lukoil_nacn_mines.py
    -> routes/econ_lukoil_nacn_mines_results.json
"""
import json
import math
import os

DIR = os.path.dirname(os.path.abspath(__file__))

# ------------------------------------------------- молярные массы, г/моль
M_NACN, M_HCN, M_NH3, M_NAOH, M_CH4, M_O2 = 49.0, 27.0, 17.0, 40.0, 16.04, 32.0


def A(value, unit, tier, comment):
    """Допущение: значение + единица + уровень достоверности + источник."""
    return {"value": value, "unit": unit, "tier": tier, "comment": comment}


# ------------------------------------------------------- допущения (все явно)
ASSUMPTIONS = {
 # --- якоря репозитория (не трогать без пересчёта всей линейки) ---
 "wacc": A(0.12, "1/год", "anchor_repo",
           "econ_shock_candidates / msa_cost_model: 12% годовых на капитал"),
 "asset_life_y": A(12, "лет", "anchor_repo",
                   "econ_vs_mining: 12 лет у скид-линейки; капитальная строка "
                   "здесь — перпетуитет capex*WACC/т, как в "
                   "econ_shock_candidates"),
 "fix_skid_usd_t_at_1kt": A(380.0, "$/т при 1 кт/г", "anchor_repo",
                            "econ_shock_candidates: фикс малого масштаба, "
                            "нормируется по тоннажу как 380*(1000/т_г)"),
 "yield_andrussow": A(0.65, "доля", "anchor_repo",
                      "econ_shock_candidates S1: выход HCN-тракта 65%"),
 "ch4_per_t_nacn": A(0.50, "т CH4/т NaCN", "anchor_repo",
                     "стехиометрия 0.327 т / выход 0.65 = 0.504 -> якорь 0.50"),
 "nh3_per_t_nacn": A(0.50, "т NH3/т NaCN", "anchor_repo",
                     "якорь S1; отвечает эфф. использованию NH3 ~69% "
                     "(стехиометрия 0.347 т; без рецикла при 65% было бы "
                     "0.534 т — узел частичного рецикла NH3 обязателен)"),
 "naoh_per_t_nacn": A(0.82, "т NaOH(100%)/т NaCN", "anchor_repo",
                      "якорь S1 = стехиометрия 0.816 т + ~0.5% избыток"),
 "gas_usd_t": A(30.0, "$/т CH4", "anchor_repo",
                "skid_twin: факельный/скважинный газ $30/т (труба $180/т — "
                "здесь не нужна, газ у ГОКа скважинный)"),
 "energy_pt_usd_t": A(180.0, "$/т NaCN", "anchor_repo",
                      "econ_shock_candidates S1: энергия + амортизация "
                      "Pt/Rh-сеток (823-15-300-328=180)"),
 "capex_skid_musd": A(6.0, "$M", "anchor_repo",
                      "econ_shock_candidates S1; разбивка в implementation"),
 "skid_t_y": A(2500.0, "т NaCN/год", "anchor_repo", "econ_shock_candidates S1"),
 "delivered_price_usd_t": A(3000.0, "$/т NaCN delivered", "anchor_repo",
                            "econ_shock_candidates S1: диапазон 2800-3500 у "
                            "рудника; в модели = gate + плечо*ставка"),
 "gok_demand_kt_y": A([3, 15], "кт NaCN/год на ГОК", "anchor_repo",
                      "notes S1: рудник-якорь 3-15 кт/г = 1-6 скидов"),
 # --- рыночные оценки (без твёрдых котировок) ---
 "nh3_usd_t_scenarios": A({"optimistic": 300, "base": 400, "pessimistic": 600},
                          "$/т NH3 delivered", "market",
                          "привозной аммиак жд+авто; якорь S1 считан "
                          "консервативно по $600 (= наш pessimistic)"),
 "naoh_usd_t": A(400.0, "$/т NaOH (100% основа) delivered", "market",
                 "каустик как 50% р-р, пересчёт на сухой; из var_s1 якоря"),
 "competitor_gate_usd_t": A(2200.0, "$/т NaCN ex-works", "market",
                            "мировая цена брикетов 1800-2500 ex-works; "
                            "delivered у рудника = gate + логистика"),
 "hazmat_logistics_usd_t_km": A(0.5, "$/т-км", "market",
                                "перевозка цианида класс 6.1: спецконвой, "
                                "охрана, страховка; диапазон 0.3-0.8"),
 "delivery_distance_km": A({"near": 800, "base": 1600, "far": 2400}, "км",
                           "market", "плечо завод->ГОК; 2200+0.5*1600=3000 "
                           "воспроизводит delivered-якорь $3000/т"),
 "gold_to_demand_elasticity": A(1.0, "б/р", "market",
                                "скрин-связка: золото ±30% -> переработка "
                                "руды и спрос NaCN ±30% (эластичность 1:1)"),
 # --- данные компании: ТОЛЬКО диапазоны, проверить у компании ---
 "saratov_nacn_capacity_kt_y": A([20, 40], "кт NaCN/год", "company_assumption",
                                 "мощность Саратоворгсинтеза по NaCN — "
                                 "проверить у компании"),
 "plant_margin_usd_t": A([200, 500], "$/т NaCN", "company_assumption",
                         "маржа завода на тонне — проверить у компании"),
 "plant_gate_price_usd_t": A(2200.0, "$/т NaCN ex-works", "company_assumption",
                             "netback завода принят равным рыночному gate — "
                             "проверить у компании"),
 "cannibal_share": A([0.3, 0.7], "доля скид-тоннажа", "company_assumption",
                     "какая часть скид-продаж вытесняет СОБСТВЕННЫЕ продажи "
                     "завода (остальное — импорт/конкуренты) — проверить"),
 "goks_with_gas": A([2, 5], "шт", "company_assumption",
                    "число ГОКов-якорей со скважинным газом в досягаемости — "
                    "проверить по месторождениям"),
 # --- литература ---
 "andrussow_trl": A(9, "TRL", "literature",
                    "CH4+NH3+1.5O2 -> HCN+3H2O на Pt/Rh(90/10)-сетке ~1100C, "
                    "промышленность 100 лет; новизна — масштаб и сайтинг"),
}

WACC = ASSUMPTIONS["wacc"]["value"]
FIX_SKID = ASSUMPTIONS["fix_skid_usd_t_at_1kt"]["value"]
YIELD = ASSUMPTIONS["yield_andrussow"]["value"]
CH4_T = ASSUMPTIONS["ch4_per_t_nacn"]["value"]
NH3_T = ASSUMPTIONS["nh3_per_t_nacn"]["value"]
NAOH_T = ASSUMPTIONS["naoh_per_t_nacn"]["value"]
GAS_USD = ASSUMPTIONS["gas_usd_t"]["value"]
NAOH_USD = ASSUMPTIONS["naoh_usd_t"]["value"]
ENERGY_PT = ASSUMPTIONS["energy_pt_usd_t"]["value"]
GATE = ASSUMPTIONS["plant_gate_price_usd_t"]["value"]
HAZ_KM = ASSUMPTIONS["hazmat_logistics_usd_t_km"]["value"]


# ------------------------------------------------- масс-баланс (честно, моли)
def mass_balance():
    kmol = 1000.0 / M_NACN                      # 20.4 кмоль NaCN на тонну
    ch4_st = kmol * M_CH4 / 1000.0              # 0.327 т
    nh3_st = kmol * M_NH3 / 1000.0              # 0.347 т
    naoh_st = kmol * M_NAOH / 1000.0            # 0.816 т
    o2 = 1.5 * (kmol / YIELD) * M_O2 / 1000.0   # 1.51 т O2 (из воздуха)
    t_y = ASSUMPTIONS["skid_t_y"]["value"]
    return {
     "reaction": "CH4 + NH3 + 1.5 O2 -> HCN + 3 H2O (Pt/Rh ~1100C); "
                 "HCN + NaOH -> NaCN + H2O",
     "per_t_nacn": {
      "hcn_t": round(kmol * M_HCN / 1000.0, 3),
      "ch4_stoich_t": round(ch4_st, 3),
      "ch4_at_65pct_t": round(ch4_st / YIELD, 3),
      "ch4_model_t": CH4_T,
      "nh3_stoich_t": round(nh3_st, 3),
      "nh3_at_65pct_no_recycle_t": round(nh3_st / YIELD, 3),
      "nh3_model_t": NH3_T,
      "nh3_note": "0.50 т = эфф. использование 69%: непрореагировавший NH3 "
                  "частично улавливается/рециклится (стандартный узел)",
      "naoh_stoich_t": round(naoh_st, 3),
      "naoh_model_t": NAOH_T,
      "o2_t": round(o2, 2),
      "air_t": round(o2 / 0.232, 1),
      "product_as_30pct_solution_t": round(1.0 / 0.30, 2)},
     "per_skid_year_at_2500t": {
      "ch4_t": round(CH4_T * t_y),
      "ch4_mmscf_d_equiv": round(CH4_T * t_y / 7300.0, 2),
      "nh3_t": round(NH3_T * t_y),
      "naoh_t_dry": round(NAOH_T * t_y),
      "note": "CH4 ~0.17 MMscf/д — малая скважина ГОКа покрывает; NH3 1250 "
              "т/г = ~60 автоцистерн/год; NaOH завозится как 50% р-р"}}


# --------------------------------------------------- экономика одного скида
def skid(nh3_usd, price, t_y, capex):
    var = CH4_T * GAS_USD + NH3_T * nh3_usd + NAOH_T * NAOH_USD + ENERGY_PT
    fixed = FIX_SKID * 1000.0 / t_y
    capital = capex * WACC / t_y
    full = var + fixed + capital
    margin = price - full
    ebitda = margin * t_y / 1e6
    return {"t_per_year": round(t_y), "price_usd_t": round(price),
            "variable_usd_t": round(var), "fixed_usd_t": round(fixed),
            "capital_usd_t": round(capital), "full_cost_usd_t": round(full),
            "margin_usd_t": round(margin), "ebitda_musd_y": round(ebitda, 2),
            "payback_y": round(capex / 1e6 / ebitda, 1) if ebitda > 0 else None}


# --------------- флот по ГОКам + каннибализация Саратоворгсинтеза (чистый эффект)
def fleet(sc):
    demand_t = sc["goks"] * sc["kt_per_gok"] * 1000.0
    n = int(math.ceil(demand_t / sc["t_y"]))
    price = GATE + sc["dist_km"] * HAZ_KM
    s = skid(sc["nh3"], price, sc["t_y"], sc["capex"])
    ebitda_fleet = s["margin_usd_t"] * demand_t / 1e6
    cann_t = demand_t * sc["cannibal"]
    lost_rev = cann_t * GATE / 1e6
    lost_margin = cann_t * sc["plant_margin"] / 1e6
    net = ebitda_fleet - lost_margin
    capex_fleet = n * sc["capex"] / 1e6
    return {"gok_anchors": sc["goks"], "kt_per_gok": sc["kt_per_gok"],
            "demand_kt_y": round(demand_t / 1000.0, 1), "n_skids": n,
            "capex_fleet_musd": round(capex_fleet, 1),
            "ebitda_fleet_musd_y": round(ebitda_fleet, 1),
            "cannibalized_kt_y": round(cann_t / 1000.0, 1),
            "plant_revenue_lost_musd_y": round(lost_rev, 1),
            "plant_margin_lost_musd_y": round(lost_margin, 1),
            "net_effect_lukoil_musd_y": round(net, 1),
            "payback_fleet_y": round(capex_fleet / net, 1) if net > 0 else None}


# --------------------------------------------- сценарии pessimistic/base/opt
SCEN = {
 "pessimistic": dict(nh3=600, dist_km=800, t_y=2000.0, capex=7.0e6,
                     goks=2, kt_per_gok=3.0, cannibal=0.7, plant_margin=500,
                     note="NH3 $600 (= консервативный якорь S1); ГОКи близко к "
                          "заводу (плечо 800 км -> delivered лишь $2600); "
                          "мелкие якоря 3 кт/г; капекс +$1M северное "
                          "исполнение; теряем максимум заводской маржи"),
 "base": dict(nh3=400, dist_km=1600, t_y=2500.0, capex=6.0e6,
              goks=3, kt_per_gok=5.0, cannibal=0.5, plant_margin=350,
              note="delivered $3000 = якорь S1; NH3 $400 (типовой рынок); "
                   "3 ГОКа по 5 кт/г = 6 скидов; половина тоннажа — бывшие "
                   "клиенты завода"),
 "optimistic": dict(nh3=300, dist_km=2400, t_y=2500.0, capex=6.0e6,
                    goks=5, kt_per_gok=8.0, cannibal=0.3, plant_margin=200,
                    note="дальние ГОКи (delivered $3400); 5 якорей по 8 кт/г "
                         "= 40 кт/г — флот масштаба самого завода (20-40 кт "
                         "company_assumption): каннибализация выше 0.3 "
                         "возможна, проверить сбытовую карту у компании"),
}


def main():
    mb = mass_balance()

    # воспроизведение якоря S1 (NH3 $600, delivered $3000, 2500 т/г, $6.0M)
    anchor = skid(600, 3000.0, 2500.0, 6.0e6)
    ok = (anchor["variable_usd_t"] == 823 and anchor["full_cost_usd_t"] == 1263
          and anchor["margin_usd_t"] == 1737
          and anchor["ebitda_musd_y"] == 4.34 and anchor["payback_y"] == 1.4)

    scenarios = {}
    for name, sc in SCEN.items():
        price = GATE + sc["dist_km"] * HAZ_KM
        scenarios[name] = {"skid": skid(sc["nh3"], price, sc["t_y"],
                                        sc["capex"]),
                           "fleet_and_cannibalization": fleet(sc),
                           "note": sc["note"]}

    # ---------------- чувствительности (один скид, база: NH3 $400, 2500 т/г,
    # капекс $6.0M, delivered $3000); ±30% или сценарный диапазон
    def brief(s):
        return {"full_usd_t": s["full_cost_usd_t"],
                "margin_usd_t": s["margin_usd_t"],
                "ebitda_musd_y": s["ebitda_musd_y"], "payback_y": s["payback_y"]}
    sens = {
     "nh3_usd_t": {str(p): brief(skid(p, 3000.0, 2500.0, 6.0e6))
                   for p in (300, 400, 600, 800)},
     "delivery_distance_km_classic": {
        str(d): brief(skid(400, GATE + d * HAZ_KM, 2500.0, 6.0e6))
        for d in (800, 1600, 2400)},
     "gok_demand_via_gold_pm30pct": {
        f"{m - 1:+.0%}": brief(skid(400, 3000.0, 2500.0 * m, 6.0e6))
        for m in (0.7, 1.0, 1.3)},
     "capex_pm30pct": {
        f"{m - 1:+.0%}": brief(skid(400, 3000.0, 2500.0, 6.0e6 * m))
        for m in (0.7, 1.0, 1.3)},
     "notes": {
      "nh3": "$800 — kill-порог (см. kill_criteria)",
      "distance": "плечо доставки классики двигает delivered-цену скида: "
                  "delivered = gate $2200 + плечо * $0.5/т-км",
      "gold": "золото ±30% -> спрос ГОКа ±30% (эластичность 1:1, скрин); "
              "+30% требует запаса шильды скида — дебатлнек не заложен",
      "capex": "±30% вокруг якорных $6.0M"}}

    # -------------------------------------------------------- kill-критерии
    kill = {
     "nh3_delivered_above_800_usd_t": {
        "threshold": 800, "unit": "$/т NH3 delivered",
        "margin_at_threshold_usd_t": skid(800, 3000.0, 2500.0,
                                          6.0e6)["margin_usd_t"],
        "why": "маржа формально ещё положительна (~$1637/т), но NH3 дороже "
               "$800/т означает сломанную северную логистику: 1250 т NH3/г "
               "на скид едут по тем же зимникам — премисса 'опасные грузы не "
               "едут' инвертируется, и завозная классика сравнивается"},
     "gok_anchor_below_2kt_y": {
        "threshold": 2000, "unit": "т NaCN/год у ГОКа-якоря",
        "payback_at_1500_t_y": skid(600, 3000.0, 1500.0, 6.0e6)["payback_y"],
        "payback_at_1000_t_y": skid(600, 3000.0, 1000.0, 6.0e6)["payback_y"],
        "why": "ниже 2 кт/г фикс+капитал не размазываются (окупаемость "
               "уезжает к 3-6 г) и весь скид висит на одном клиенте"},
     "hcn_production_ban_at_mine": {
        "threshold": "бинарный (регуляторика)",
        "why": "запрет производства HCN на промплощадке рудника (лицензия "
               "ОПО, Ростехнадзор) убивает вариант целиком; снимается на "
               "FEED-стадии ДО денег на серию — см. trl_gates"},
    }

    # ------------------------------------------------ механическая реализация
    implementation = {
     "block_flow": [
      "1. приём скважинного газа + сероочистка (ZnO)",
      "2. склад привозного NH3 (изотермический, 2x100 т) + узел испарения",
      "3. воздушный компрессор + фильтрация",
      "4. смеситель CH4/NH3/воздух (контроль состава ниже НКПР, пламегасители)",
      "5. реактор Андруссова: Pt/Rh(90/10)-сетки, ~1100C, мс-контакт",
      "6. котёл-утилизатор (пар на собственные нужды — реакция экзотермична)",
      "7. улавливание/частичный рецикл непрореагировавшего NH3",
      "8. абсорбция HCN в 30% NaOH -> раствор NaCN 30%",
      "9. (опция) выпарка + грануляция -> брикеты в биг-бэгах",
      "10. дожиг хвостов (следы HCN) + свеча",
      "11. ёмкостной парк NaCN + отгрузка на ЗИФ (короткое плечо/трубопровод)"],
     "equipment": {
      "реактор Андруссова + первый комплект Pt/Rh-сеток": 1.5e6,
      "котёл-утилизатор": 0.5e6,
      "абсорберы (NH3-улавливание + HCN->NaOH)": 0.9e6,
      "склад NH3 (2x100 т) + насосы/испаритель": 0.6e6,
      "склад NaOH (50% р-р) + парк готового NaCN": 0.5e6,
      "компрессоры воздуха + газоподготовка": 0.5e6,
      "КИПиА/ПАЗ: газоанализ HCN/NH3, дожиг, СБ": 0.7e6,
      "монтаж/фундаменты/эстакады, северное исполнение": 0.8e6,
      "note": "сумма $6.0M = якорный капекс; Pt в сетках (~$0.2M) возвратный"},
     "footprint_m2": 1500,
     "footprint_note": "установка ~1500 м2 + санитарно-защитная зона; "
                       "площадка на промзоне ГОКа рядом с ЗИФ",
     "utilities": {
      "электричество_кВт": 500,
      "пар": "самообеспечение от котла-утилизатора",
      "вода_оборотная_м3_ч": 10,
      "азот_нм3_ч": 50,
      "прочее": "H2SO4 для NH3-улавливания (если рецикл неполный)"},
     "staffing": "10-12 чел: 4 вахтовые бригады x 2 оператора + КИП/механик "
                 "+ начальник установки; ремонты — силами ГОКа по договору",
     "rollout": [
      {"stage": "FEED + лицензирование производства HCN/NaCN на площадке",
       "months": 9},
      {"stage": "пилотный скид на 1-м ГОКе: пуск, выход >=60-65%, "
                "сертификация NaCN у ЗИФ", "months": 12},
      {"stage": "решение о серии; тираж 2-3 скида/год", "months": 12},
      {"stage": "целевой флот 6 скидов / 3 ГОКа (base)", "months": 15}],
     "trl_gates": [
      "выход HCN >=60% подтверждён на пилотной сетке (химия TRL 9, но малый "
      "масштаб сетки — мерить, не верить)",
      "CFD/кинетика смесителя: состав вне взрывных пределов во всех режимах "
      "— расчётом ДО FEED",
      "тест выщелачивания руды ГОКа нашим 30% р-ром NaCN — приёмка ЗИФ",
      "лицензия на производство HCN у рудника получена ДО денег на серию",
      "контракт зимней NH3-логистики <= $600/т подтверждён"],
     "barrier_and_safety": {
      "категория_перевозки_исчезает":
        "цианид не едет 1500-3000 км спецконвоем по дорогам общего "
        "пользования — он рождается в сотнях метров от ЗИФ и уходит "
        "внутризаводским плечом; издержки категории 'транспорт ЦН' (конвой, "
        "охрана, страховка, разрешения на маршрут) исчезают из цепочки; "
        "остаётся привоз NH3 (рутинная практика) и NaOH (класс ниже)",
      "барьер_входа":
        "чтобы отбить клиента, конкуренту нужен свой скид — а площадка, газ "
        "и интеграция с ЗИФ уже заняты; take-or-pay с ГОКом на 5-7 лет"}}

    OUT = {
     "model": "Лукойл L2: NaCN-скиды у золотых ГОКов + каннибализация "
              "Саратоворгсинтеза",
     "tier": "screening_pm_2_3x",
     "not_bankable": True,
     "frame": "скрин ±2-3x; фикс $380/т@1кт/г нормируется по тоннажу; WACC "
              "12%; капитальная строка — перпетуитет capex*WACC/т (как в "
              "econ_shock_candidates)",
     "consistency_anchors": {
      "econ_shock_candidates.json": "S1: 2500 т/г, $3000/т delivered, var "
                                    "$823, full $1263, маржа $1737, EBITDA "
                                    "$4.34M, окупаемость 1.4 г, выход 65%; "
                                    "WACC 12%",
      "skid_twin.json": "газ факел $30/т",
      "econ_vs_mining_results.json": "NaCN $3000/т >> порога 3x к майнингу "
                                     "($751/т) — категория specialty"},
     "assumptions": ASSUMPTIONS,
     "mass_balance": mb,
     "anchor_reproduction": {"inputs": "NH3 $600, delivered $3000, 2500 т/г, "
                                       "капекс $6.0M", "skid": anchor,
                             "matches_econ_shock_candidates_S1": ok},
     "scenarios": scenarios,
     "sensitivity": sens,
     "kill_criteria": kill,
     "implementation": implementation,
     "honesty": "скрин ±2-3x, НЕ банковская модель. Всё по Лукойлу/"
                "Саратоворгсинтезу (мощность 20-40 кт/г, маржа $200-500/т, "
                "netback, доля каннибализации, список ГОКов с газом) — "
                "company_assumption, проверить у компании. Выход 65% — "
                "литературный уровень Андруссова, на пилоте подтвердить. "
                "Цены NH3/NaOH/логистики — рыночные оценки без котировок. "
                "В optimistic флот 40 кт/г сопоставим с самим заводом — "
                "каннибализация там почти наверняка выше заложенных 0.3."}

    # ------------------------------------------------------------- сводка
    print("=== Лукойл L2: NaCN-скиды у золотых ГОКов (скрин ±2-3x) ===")
    print(f"якорь S1 воспроизведён: {ok} (var $823, full $1263/т, маржа "
          f"$1737/т, EBITDA $4.34M, окуп. 1.4 г)")
    print(f"масс-баланс на т NaCN: CH4 0.50 т (стех. 0.327), NH3 0.50 т "
          f"(стех. 0.347, эфф. 69% с рециклом), NaOH 0.82 т (стех. 0.816); "
          f"на скид 2500 т/г: NH3 1250 т, NaOH 2050 т, CH4 1250 т "
          f"(~0.17 MMscf/д)")
    print(f"\n{'сценарий':<13}{'full$/т':>8}{'маржа$/т':>9}{'скидов':>7}"
          f"{'EBITDA флота':>13}{'потеря завода':>14}{'чистый эфф.':>12}"
          f"{'окуп.':>7}")
    for name, s in scenarios.items():
        f_ = s["fleet_and_cannibalization"]
        sk = s["skid"]
        print(f"{name:<13}{sk['full_cost_usd_t']:>8}{sk['margin_usd_t']:>9}"
              f"{f_['n_skids']:>7}{f_['ebitda_fleet_musd_y']:>12.1f}M"
              f"{f_['plant_margin_lost_musd_y']:>13.1f}M"
              f"{f_['net_effect_lukoil_musd_y']:>11.1f}M"
              f"{f_['payback_fleet_y']:>6}г")
    b = sens["nh3_usd_t"]
    print(f"\nчувств. NH3 $/т (full$/т): 300->{b['300']['full_usd_t']}, "
          f"400->{b['400']['full_usd_t']}, 600->{b['600']['full_usd_t']}, "
          f"800->{b['800']['full_usd_t']} (kill-порог)")
    print("kill: NH3>$800/т; ГОК-якорь<2 кт/г; запрет производства HCN у "
          "рудника")

    path = os.path.join(DIR, "econ_lukoil_nacn_mines_results.json")
    with open(path, "w") as f:
        json.dump(OUT, f, indent=1, ensure_ascii=False)
    print(f"\nwrote {os.path.basename(path)}")


if __name__ == "__main__":
    main()
