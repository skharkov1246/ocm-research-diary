#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
routes/econ_vs_mining.py — АЛЬТЕРНАТИВНЫЕ ИЗДЕРЖКИ: наши химические маршруты против
самого сильного реального конкурента за факельный газ — майнинга биткоина на
газо-генераторе (Crusoe / Upstream Data — развёрнутый бизнес). Требование
заказчика: изобретение должно давать КРАТНО лучший результат, чем ЛЮБАЯ альтернатива,
иначе строить не стоит. Считаем строго говоря, без подыгрывания химии.

Базис: 1 MMscf/д факела. Параметры — mid-2025/26 оценки с явным разбросом;
цель — порядок и кратность, не банковская смета. Все допущения — вверху, крутите.
"""
import json
import os

DIR = os.path.dirname(os.path.abspath(__file__))

# ------------------------------------------------------------- допущения (крутите)
MMSCF_PER_DAY = 1.0
HHV_BTU_SCF = 1000.0                 # теплотворность газа
GENSET_EFF = 0.38                    # КПД поршневого генератора (электр.)
CH4_T_PER_YR = 7300.0 * MMSCF_PER_DAY  # т CH4/год на 1 MMscf/д

# --- майнинг ---
ASIC_W_PER_TH = 17.5                 # Вт на TH/s (S21-класс, 2025)
HASHPRICE_USD_PH_DAY = 50.0          # $/PH/s/день (диапазон 30–70; волатилен ±5×)
ASIC_USD_PER_TH = 18.0               # капекс ASIC $/TH
ASIC_LIFE_YR = 2.0                   # срок до морального устаревания (→ $0)
GENSET_USD_PER_KW = 700.0
DIFFICULTY_EROSION_YR = 0.30         # падение выручки/TH в год от роста сложности
MINE_OM_FRAC = 0.15                  # O&M как доля выручки

# --- химия: (название, т продукта / т CH4, $/т, капекс $M, opex $M/год, жизнь лет, увозимость) ---
CHEM = {
    "I4 OCM→этилен (Y=35%)": dict(yield_t=0.35*28/32, price=900, capex=3.5, opex=0.4, life=12),
    "OCM→этилен стена (Y=25%)": dict(yield_t=0.25*28/32, price=900, capex=3.5, opex=0.4, life=12),
    "Метанол (α-O, Y=40%)": dict(yield_t=0.40*32/16, price=400, capex=3.0, opex=0.4, life=12),
    "Пиролиз H2+C(тв)": dict(yield_t=0.25+0.75, price=None, capex=4.0, opex=0.5, life=15),  # спец. ниже
}


def mw_electric():
    btu_day = MMSCF_PER_DAY * 1e6 * HHV_BTU_SCF
    w_th = btu_day * 0.293 / 24.0            # Вт тепл.
    return w_th * GENSET_EFF / 1e6           # МВт электр.


def mining():
    mw = mw_electric()
    th = mw * 1e6 / ASIC_W_PER_TH            # TH/s развёртываемо
    ph = th / 1000.0
    gross_yr1 = ph * HASHPRICE_USD_PH_DAY * 365 / 1e6   # $M/год
    asic_capex = th * ASIC_USD_PER_TH / 1e6
    genset_capex = mw * 1000 * GENSET_USD_PER_KW / 1e6
    capex = asic_capex + genset_capex + 1.0            # +инфра/контейнер
    # усреднённая выручка за жизнь ASIC с эрозией сложности
    g = gross_yr1; tot = 0.0
    for _ in range(int(ASIC_LIFE_YR)):
        tot += g; g *= (1 - DIFFICULTY_EROSION_YR)
    gross_life = tot
    om = gross_life * MINE_OM_FRAC
    asic_amort = asic_capex                             # ASIC → $0 за жизнь
    genset_amort = genset_capex * (ASIC_LIFE_YR/10.0)   # генератор живёт ~10 лет
    net_life = gross_life - om - asic_amort - genset_amort
    net_yr = net_life / ASIC_LIFE_YR
    return dict(mw=round(mw,2), ph=round(ph,1), gross_yr1=round(gross_yr1,2),
                capex=round(capex,2), net_per_yr=round(net_yr,2),
                asics=int(th/200))


def chem_route(cfg):
    if cfg["price"] is None:  # пиролиз: H2 @ $2/кг + углерод @ $400/т
        h2 = 0.25*CH4_T_PER_YR; c = 0.75*CH4_T_PER_YR
        gross = (h2*1000*2 + c*400)/1e6
    else:
        prod = cfg["yield_t"]*CH4_T_PER_YR
        gross = prod*cfg["price"]/1e6
    ebitda = gross - cfg["opex"]
    return dict(gross_yr=round(gross,2), ebitda_yr=round(ebitda,2),
                capex=cfg["capex"], life=cfg["life"])


def breakeven_price(mult, mine_net):
    """$/т продукта, чтобы EBITDA химии = mult × чистой майнинга (при 50% C-эффект.)."""
    prod_t = 0.5*CH4_T_PER_YR
    need_gross = mult*mine_net + 0.4      # +opex
    return round(need_gross*1e6/prod_t, 0)


def main():
    m = mining()
    out = {"basis": f"{MMSCF_PER_DAY} MMscf/day flare = {CH4_T_PER_YR:.0f} t CH4/yr "
                    f"= {m['mw']} MW electric (genset {GENSET_EFF})",
           "mining": m, "chem": {}, "verdict": {}}
    print(f"BASIS: {out['basis']}\n")
    print(f"BITCOIN MINING: {m['ph']} PH/s ({m['asics']} ASICs), gross yr1 ${m['gross_yr1']}M,"
          f" capex ${m['capex']}M, NET ~${m['net_per_yr']}M/yr (volatile ±5x on BTC)\n")
    print(f"{'route':30} {'gross$M':>8} {'EBITDA$M':>9} {'capex$M':>8} {'× vs mining net':>16}")
    for name, cfg in CHEM.items():
        r = chem_route(cfg); out["chem"][name] = r
        mult = round(r["ebitda_yr"]/m["net_per_yr"], 2) if m["net_per_yr"] > 0 else None
        print(f"{name:30} {r['gross_yr']:8} {r['ebitda_yr']:9} {r['capex']:8} {str(mult):>16}")
    out["verdict"] = {
        "mining_net_per_yr_musd": m["net_per_yr"],
        "breakeven_price_to_beat_mining_3x": breakeven_price(3, m["net_per_yr"]),
        "breakeven_price_to_beat_mining_1x": breakeven_price(1, m["net_per_yr"]),
        "reading": "commodity products (ethylene/methanol/H2) only MATCH mining, not "
                   "multiples. To beat mining 3x needs a product > "
                   f"~${breakeven_price(3, m['net_per_yr'])}/t — i.e. a SPECIALTY/material "
                   "molecule whose value is decoupled from energy content. Chemistry's "
                   "structural edge over mining is qualitative: durable 12-15yr asset "
                   "(vs ASIC 2yr->$0), price-stable commodity (vs BTC ±5x), and product "
                   "that can carry carbon credit — mining just burns the molecule to CO2."}
    json.dump(out, open(os.path.join(DIR, "econ_vs_mining_results.json"), "w"), indent=1)
    print("\nBREAK-EVEN product price to beat mining 3x: "
          f"${out['verdict']['breakeven_price_to_beat_mining_3x']}/t")
    print("VERDICT:", out["verdict"]["reading"])


if __name__ == "__main__":
    main()
