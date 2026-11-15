"""Блок 1. Подготовка топлива: нет / подогрев топливного газа / паровой риформинг / сухой риформинг на рециркуляте.

Модели (оценка уровня облика):
  * подогрев - в цепочке (chain.py): газ после турбины греет топливо до min(T_газа - 40 K, 400 °C), потеря давления по
    газу 0,2 % (organic/cycles: A['T_fh_max'], A['dT_app'], A['dp_fh_g']);
  * паровой и сухой риформинг меняют состав топлива; таблицы свойств цепочки его не описывают. Приращение КПД, удельной
    работы, температуры выхлопа и расхода выхлопа берётся из перебора циклов (cycles.evaluate, Cantera, равновесие с
    подходом 15 K, пинч 20 K) - как разность «с риформингом» минус «без» при той же π, T41, рекуператоре, промохлаждении,
    подогреве - и переносится на концепцию (суперпозиция приращений - допущение, отмечается в выдаче).
Минус Plug&Play у парового риформинга - вода (водоподготовка, насос, бак): пакеты воды не имеют.
"""
from . import _track as T

OPTIONS = {
    "нет": dict(name="топливный газ пакета как есть", tcr="none", water=False,
                note="точка трека: метан, 15 °C (gtu_cycle, gtu_explore)"),
    "подогрев": dict(name="подогрев топливного газа выхлопом до 400 °C", tcr="fh", water=False,
                     note="cycles: +0,55 п. п. у простого KV-3 за 2 детали (README organic/cycles, оценка)"),
    "паровой риформинг": dict(name="паровой риформинг метана теплом выхлопа (химическая рекуперация, нужна вода)", tcr="steam", water=True,
                              note="cycles: выигрыш простого KV-3 почти целиком от впрыска пара, а не от химии; химия заметна при газе "
                                   "730 °C (рекуператор + подогрев) - оценка"),
    "сухой риформинг": dict(name="сухой риформинг на рециркуляте выхлопа (без воды)", tcr="egr", water=False,
                            note="cycles: кислорода в выхлопе в 3-4 раза больше, чем CO2 и H2O, - он сжигает метан, теплотворность падает, "
                                 "сажа термодинамически выгодна (активность углерода 7-10) - оценка"),
}


def trl(prep):
    CY = T.CY()
    tcr = OPTIONS[prep]["tcr"]
    if tcr == "none":
        return 9, "традиция"
    return CY.TRL[tcr], CY.TRL_NOTE[tcr]


def chain_spec(C):
    """параметры подогрева топлива для цепочки (None - без подогрева)"""
    if C["fuel"]["prep"] != "подогрев":
        return None
    A = T.CY().A
    return dict(T_max=float(C["fuel"]["T_max"] or A["T_fh_max"]), dT_app=float(A["dT_app"]), dp_g=float(A["dp_fh_g"]))


def cycles_topo(C, tcr):
    """топология перебора циклов, ближайшая к концепции (для приращения химической рекуперации)"""
    hx, cb = C["hx"], C["combustor"]
    return dict(rec=float(hx["recup"] or 0.0), ic=C["compressor"]["intercool"] != "нет", rh=bool(hx["reheat"]),
                ibc=bool(hx["ibc"]), pgc="real" if cb["scheme"] in ("RDC", "ВР") else "none", tcr=tcr)


def reform_delta(C):
    """приращение от риформинга по перебору циклов при π, T41 и прочих блоках концепции. Возврат: dict или None"""
    prep = C["fuel"]["prep"]
    tcr = OPTIONS[prep]["tcr"]
    if tcr not in ("steam", "egr"):
        return None
    CY = T.CY()
    PR = float(C["compressor"]["pi"])
    prm = dict(T41_C=float(C["turbine"]["t41"]), x_rh=float(C["hx"]["x_rh"]), sc=float(C["fuel"]["sc"]), R_egr=float(C["fuel"]["R_egr"]))
    if C["options"].get("T_rec_max"):
        prm["T_rec_max"] = float(C["options"]["T_rec_max"])
    on = CY.evaluate(cycles_topo(C, tcr), PR, prm)
    off = CY.evaluate(cycles_topo(C, "none"), PR, prm)
    if not (on.get("ok") and off.get("ok")):
        return dict(ok=False, why=f"перебор циклов: {on.get('why') or off.get('why')}", PR=PR, prm=prm)
    return dict(ok=True, PR=PR, prm=prm, topo=cycles_topo(C, tcr),
                d_eff=on["eff"] - off["eff"], k_w=on["w_kJkg"] / off["w_kJkg"], d_EGT=on["EGT_C"] - off["EGT_C"],
                d_mexh=on["W_exh"] / on["W_air"] - off["W_exh"] / off["W_air"],
                water_t_h_per_kgs=on["water_t_h"] / on["W_air"], fuel_info=on.get("fuel_info") or {},
                ref=dict(eff_on=on["eff"], eff_off=off["eff"], w_on=on["w_kJkg"], w_off=off["w_kJkg"], EGT_on=on["EGT_C"],
                         EGT_off=off["EGT_C"]),
                method="cycles.evaluate (organic/cycles, Cantera): разность «с риформингом» - «без» при тех же π, T41, ε, "
                       "промохлаждении, подогреве; перенесена на цепочку концепции (суперпозиция приращений - оценка)")


def report(C, ch, delta):
    prep = C["fuel"]["prep"]
    t, tb = trl(prep)
    out = dict(variant=prep, name=OPTIONS[prep]["name"], trl=t, trl_basis=tb, water=OPTIONS[prep]["water"], note=OPTIONS[prep]["note"],
               T_fuel_C=round(ch["T_fuel"] - 273.15, 1) if ch.get("ok") else None)
    if prep == "подогрев" and ch.get("ok"):
        out.update(q_fuel_kJ_per_kg_air=round(ch["q_fuel"] / 1e3, 3))
    if delta:
        out["reform"] = {k: (round(v, 5) if isinstance(v, float) else v) for k, v in delta.items() if k not in ("fuel_info",)}
        fi = delta.get("fuel_info") or {}
        out["reform"]["fuel_info"] = {k: (round(v, 4) if isinstance(v, float) else v) for k, v in fi.items()
                                      if isinstance(v, (int, float, str))}
    return out
