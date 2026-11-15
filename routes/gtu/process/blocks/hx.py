"""Блок 5. Теплообменники и теплообмен цикла: рекуператор (ε 0,70-0,92; конструкция - фольга с первичной поверхностью
или печатная сердцевина LPBF), повторный подогрев на одном валу (межтурбинная камера), обращённый цикл Брайтона (IBC);
здесь же размеры промохладителя (сам охладитель - блок компрессора) и подогревателя топлива.

Потери давления рекуператора - gtu_cycle.recup по умолчанию (2,5 % по воздуху, 3,5 % по газу); газ на входе - не выше
LIM['T_rec_max'] 700 °C (gtu_explore). Подогрев и IBC - параметры перебора циклов (organic/cycles A: dp_rh 3 %, x_rh;
eta_ibc_t 0,86, eta_ibc_c 0,82, eps_ibc 0,85, dp_ibc 2 %, сетка plow). Размеры сердцевин - cycles.size_hx (ламинарные плоские
каналы, f·Re и Nu посчитаны конечными разностями; технологии cycles.HX_TECH: фольга 0,1 мм / LPBF IN625 0,35 мм) по
нагрузке и температурам цепочки концепции. Детали узлов - cycles.parts (те же позиции, что у перебора циклов). Всё - оценка.
"""
import inspect
import math

import numpy as np

from . import _track as T

TECH = {"фольга": "foil", "печать": "lpbf"}
RECUP_EPS = (0.70, 0.92)           # диапазон задачи владельца для блока (перебор циклов - 0,70-0,94)


def recup_params(C):
    """потери рекуператора: значения по умолчанию gtu_cycle.recup (читаются из сигнатуры, файл не меняется)"""
    sig = inspect.signature(T.GC().recup).parameters
    return dict(eps=float(C["hx"]["recup"]), dpa=float(sig["dp_air"].default), dpg=float(sig["dp_gas"].default),
                model=C["options"].get("rec_model", "enthalpy"))


def T_rec_max(C):
    v = C["options"].get("T_rec_max")
    return float(v) if v else float(T.FS().T_REC_MAX)


def reheat_params(C):
    A = T.CY().A
    return dict(x_rh=float(C["hx"]["x_rh"]), t_rh=float(C["hx"]["t_rh"]), dp=float(A["dp_rh"]))


def ibc_params():
    A = T.CY().A
    return dict(eta_t=float(A["eta_ibc_t"]), eta_c=float(A["eta_ibc_c"]), eps=float(A["eps_ibc"]), dp=float(A["dp_ibc"]),
                plow=tuple(A["plow"]), Cr=float(A["Cr_ic"]), dp_fan=float(A["dp_fan"]), eta_fan=float(A["eta_fan"]),
                eta_motor=float(A["eta_motor"]))


def static_rules(C):
    hx = C["hx"]
    R = []
    if hx["recup"]:
        lo, hi = RECUP_EPS
        R.append(dict(code="recup_eps", hard=False, ok=lo - 1e-9 <= hx["recup"] <= hi + 1e-9, val=hx["recup"], lim=[lo, hi],
                      text="эффективность рекуператора в диапазоне задачи 0,70-0,92"))
    if hx["reheat"] and C["turbine"]["type"] == "Р" and C["rotors"] == 1:
        R.append(dict(code="reheat_radial", hard=True, ok=False,
                      text="подогрев на одном валу между ступенями радиальной турбины: межтурбинной камере нет места между "
                           "колёсами (оценка) - нужна осевая турбина"))
    if hx["ibc"] and hx["recup"] and hx["recup"] >= 0.85:
        R.append(dict(code="ibc_after_recup", hard=False, ok=False,
                      text="IBC за рекуператором с ε >= 0,85 получает холодный газ - выигрыш мал (перебор циклов: такие сочетания "
                           "часто «бесполезны»)"))
    return R


def trl_items(C):
    CY = T.CY()
    hx = C["hx"]
    out = []
    if hx["recup"] and hx["recup_make"] == "фольга":
        out.append(("рекуператор: фольга, первичная поверхность", CY.TRL["rec"], CY.TRL_NOTE["rec"]))
    if hx["reheat"]:
        out.append(("повторный подогрев (межтурбинная камера)", CY.TRL["rh"], CY.TRL_NOTE["rh"]))
    if hx["ibc"]:
        out.append(("обращённый цикл Брайтона", CY.TRL["ibc"], CY.TRL_NOTE["ibc"]))
    return out


# =========================================================================== размеры сердцевин (cycles.size_hx)
_YST = None


def Y_exh(b):
    """массовые доли (порядок cycles.SP) газа с долей стехиометрических продуктов метана b"""
    global _YST
    CY = T.CY()
    if _YST is None:
        I, MW = CY.I, CY.MW
        n = CY.Y_AIR / MW
        nf = n[I["O2"]] / 2.0
        n = n.copy()
        n[I["O2"]] = 0.0
        n[I["CO2"]] += nf
        n[I["H2O"]] += 2 * nf
        y = n * MW
        _YST = y / y.sum()
    return (1 - b) * CY.Y_AIR + b * _YST


def _size(Q, Th_in, Th_out, Tc_in, Tc_out, hot_m, hot_P, hot_Y, dPh, cold_m, cold_P, cold_Y, dPc, tech):
    CY = T.CY()
    ua = CY.ua_counterflow(Q, Th_in, Th_out, Tc_in, Tc_out)
    if not ua:
        return dict(ok=False, why="пересечение температур")
    s = CY.size_hx(ua, hot_m, 0.5 * (Th_in + Th_out), hot_P, hot_Y, dPh, cold_m, 0.5 * (Tc_in + Tc_out), cold_P, cold_Y, dPc, tech)
    s["ok"] = True
    s["Q_kW"] = Q / 1e3
    return {k: (round(v, 5) if isinstance(v, float) else v) for k, v in s.items()}


def sizes(C, ch, G):
    """сердцевины теплообменников концепции при расходе воздуха G, кг/с: рекуператор, промохладитель, охладитель IBC,
    подогреватель топлива. Возврат: {имя: {tech: размеры}}"""
    if not ch.get("ok"):
        return {}
    CY, E = T.CY(), T.E()
    out = {}
    techs = list(TECH.values())
    if C["hx"]["recup"]:
        rp = recup_params(C)
        Tg_in, Tg_out = ch["T_rec_in"], ch["T_stack_rec"]
        out["рекуператор"] = {tc: _size(ch["Q_rec"] * G, Tg_in, Tg_out, ch["T2"], ch["T3"], ch["m_exh"] * G, ch["p_exit"] / (1 + rp["dpg"]),
                                        Y_exh(ch["b_exh"]), rp["dpg"] * ch["p_exit"], ch["mc"] * G, ch["p2"], CY.Y_AIR, rp["dpa"] * ch["p2"], tc)
                              for tc in techs}
    if ch.get("ic"):
        ic = ch["ic"][0]
        icp = T.CY().A
        Tc_out = E.T0 + ic["Q"] / max(ic["m_cold"], 1e-9) / 1006.0
        out["промохладитель"] = {tc: _size(ic["Q"] * G, ic["T_in"], ic["T_out"], E.T0, Tc_out, 1.0 * G, ic["p"], CY.Y_AIR,
                                           icp["dp_ic"] * ic["p"], ic["m_cold"] * G, E.P0, CY.Y_AIR, icp["dp_fan"], tc) for tc in techs}
    b = ch.get("ibc")
    if b and not b.get("useless"):
        Tc_out = E.T0 + b["Q"] / max(b["m_cold"], 1e-9) / 1006.0
        P_hot = b["plow"] * E.P0
        out["охладитель IBC"] = {tc: _size(b["Q"] * G, b["T_t_out"], b["T_cool"], E.T0, Tc_out, ch["m_exh"] * G, P_hot, Y_exh(ch["b_exh"]),
                                           T.CY().A["dp_ibc"] * P_hot, b["m_cold"] * G, E.P0, CY.Y_AIR, T.CY().A["dp_fan"], tc)
                                 for tc in techs}
    if C["fuel"]["prep"] == "подогрев" and ch["q_fuel"] > 0:
        A = T.CY().A
        Tg_in = ch["T_turb_exit"]
        Tg_out = E.T_h(E.h(Tg_in, ch["b_exh"]) - ch["q_fuel"] / ch["m_exh"], ch["b_exh"])
        p_fuel = ch["p3"] * A["p_fuel_k"]
        out["подогреватель топлива"] = {tc: _size(ch["q_fuel"] * G, Tg_in, Tg_out, E.T0, ch["T_fuel"], ch["m_exh"] * G, ch["p_exit"],
                                                  Y_exh(ch["b_exh"]), A["dp_fh_g"] * ch["p_exit"], ch["fuel"] * G, p_fuel, CY.Y_CH4,
                                                  A["dp_ref_f"] * p_fuel, tc) for tc in techs}
    return out


def funnel_rec(C, sz, rec_funnel):
    """рекуператор концепции в формате воронки (funnel_space.derive -> g['rec']) по размерам cycles.size_hx выбранной
    технологии: объём с коллекторами и корпусом ×2 (правило cycles.parts), модули - кубы камеры EP-M2050 1,1 м (как derive),
    масса - сердцевина size_hx (стенка и коэффициент конструкции технологии). Нужен моделям деталей, обслуживания, Plug&Play
    и габарита воронки вместо её собственной оценки печатного рекуператора."""
    s = (sz.get("рекуператор") or {}).get(TECH[C["hx"]["recup_make"]])
    if not s or not s.get("ok"):
        return rec_funnel
    V = 2.0 * s["V_core_m3"]
    a_max = 1.1
    mods = max(1, math.ceil(V / a_max ** 3))
    return dict(A_m2=s["A_m2"], V_m3=V, modules=mods, side=(V / mods) ** (1.0 / 3.0), mass_kg=s["mass_kg"],
                source=f"cycles.size_hx ({TECH[C['hx']['recup_make']]}), объём ×2, модули по кубу 1,1 м",
                funnel_estimate=({k: round(v, 3) for k, v in rec_funnel.items()} if rec_funnel else None))


def ibc_wheel(C, ch, G, n_shaft):
    """колесо компрессора IBC на валу двигателя (как cycles.ibc_run): коэффициент расхода при частоте вала; больше 0,15
    (cycles A phi_max) - нужен свой вал и редуктор"""
    b = ch.get("ibc")
    if not b or b.get("useless"):
        return None
    E = T.E()
    rho1 = b["plow"] * E.P0 * (1 - T.CY().A["dp_ibc"]) / (E.Rg(ch["b_exh"]) * b["T_cool"])
    U2 = math.sqrt(b["w_c_per_kg"] / 0.70)
    D2 = 60 * U2 / (math.pi * n_shaft)
    m_abs = ch["m_exh"] * G
    phi = m_abs / (rho1 * U2 * D2 ** 2)
    return dict(U2=round(U2, 1), D2_shaft=round(D2, 4), phi_shaft=round(phi, 4), phi_max=T.CY().A["phi_max"],
                D2_opt=round(math.sqrt(m_abs / (rho1 * U2 * 0.085)), 4), own_shaft=phi > T.CY().A["phi_max"])


def aux_parts(C, sz, ibc_w, limits_all):
    """добавочные узлы теплообмена и подготовки топлива - cycles.parts по той же ведомости (на технологию рекуператора
    концепции); возврат: items, added (без колёс), pnp (да / частично / нет по правилу перебора циклов), V_added_m3, масса"""
    CY = T.CY()
    hx = C["hx"]
    tech = TECH[hx["recup_make"]] if hx["recup"] else "lpbf"
    hx_list = []
    for name, d in sz.items():
        s = d.get(tech) or {}
        if s.get("ok"):
            nm = {"подогреватель топлива": "подогреватель топлива"}.get(name, name)
            hx_list.append(dict(name=nm, modules=s.get("modules"), V_core_m3=s["V_core_m3"], mass_kg=s["mass_kg"]))
    tcr = {"нет": "none", "подогрев": "fh", "паровой риформинг": "steam", "сухой риформинг": "egr"}[C["fuel"]["prep"]]
    if tcr in ("steam", "egr"):
        hx_list.append(dict(name="риформер", modules=1, V_core_m3=0.0, mass_kg=0.0))
    topo = dict(rec=float(hx["recup"] or 0.0), ic=C["compressor"]["intercool"] != "нет", rh=bool(hx["reheat"]), ibc=bool(hx["ibc"]),
                pgc="real" if C["combustor"]["scheme"] in ("RDC", "ВР") else "none", tcr=tcr)
    res = dict(topo=topo, nc=2, nt=3, limits=dict(all=bool(limits_all)),
               ibc=(dict(phi_shaft=ibc_w["phi_shaft"]) if ibc_w else (dict(useless=True) if hx["ibc"] else None)))
    p = CY.parts(res, hx_list, tech)
    return dict(items=[dict(name=n, count=k) for n, k in p["items"]], added=sum(k for _, k in p["items"]), pnp=p["pnp"],
                V_added_m3=round(p["V_added_m3"], 3), mass_hx_kg=round(p["mass_hx_kg"], 1), tech=tech,
                maintenance=p["maintenance"], consumables=p["consumables"], water=p["water"],
                source="organic/cycles: parts() - позиции и правило Plug&Play (вода или вне допуска - нет; без теплообменников и "
                       "экзотики - да; объём ×2 до 6 м³ - частично)")


def report(C, ch, sz, ibc_w):
    hx = C["hx"]
    out = dict(recup=hx["recup"], recup_make=hx["recup_make"] if hx["recup"] else None, reheat=bool(hx["reheat"]), ibc=bool(hx["ibc"]))
    if hx["recup"]:
        rp = recup_params(C)
        out["recup_params"] = dict(rp, T_rec_max_C=T_rec_max(C), source="gtu_cycle.recup (dp_air, dp_gas по умолчанию); LIM T_rec_max")
        if ch.get("ok"):
            out.update(T_rec_in_C=round(ch["T_rec_in"] - 273.15, 1), T3_C=round(ch["T3"] - 273.15, 1),
                       Q_rec_kJ_per_kg_air=round(ch["Q_rec"] / 1e3, 2))
    if hx["reheat"]:
        out["reheat_params"] = dict(reheat_params(C), source="organic/cycles A dp_rh; x_rh - доля ln p до подогрева (одновальный)")
    if hx["ibc"]:
        b = ch.get("ibc") if ch.get("ok") else None
        out["ibc_params"] = ibc_params()
        if b:
            out["ibc_result"] = {k: (round(v, 4) if isinstance(v, float) else v) for k, v in b.items()}
            out["ibc_wheel"] = ibc_w
    out["sizes"] = sz
    return out
