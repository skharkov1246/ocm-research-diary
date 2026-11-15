#!/usr/bin/env python3
"""Оценка сквозной концепции «топливо -> вал»: блоки (process/blocks) + модель массы (process/mass) -> КПД, масса, кВт/кг,
«КПД на тонну», детали, Plug&Play, готовность, допустимость с причинами.

Критерий владельца 15.11.2026 «КПД на тонну»: для типоразмера с фиксированной мощностью на валу (стык пакета) - КПД
на клеммах (у генераторных пакетов; у приводных - на валу) и масса агрегата-заменителя: двигатель + камера +
теплообменники + корпус и рама модуля; генератор, редуктор, вход и выхлоп пакета не входят. «КПД на тонну» = КПД % /
масса т. КПД - цепочка блоков (сверка с frame3 -0,50 п. п. при калибровке трека ±2 п. п., blocks_results.json);
масса - модель, калиброванная по 8 реальным двигателям (промышленная практика Solar: K 2,79, LOO 22,8 %;
mass_calibration.json). Ничего не подгоняется под концепцию: обе модели только вызываются.

Порядок:
  1. decode генома (process_space) -> концепция блоков;
  2. быстрый отказ по статическим правилам совместимости (blocks.compat.check_static) - без расчёта цикла;
  3. blocks.evaluate: КПД, расход и выхлоп против допуска пакета, детали, Plug&Play, TRL, размеры, правила;
  4. только у допустимых по блокам и без «невозможно» в Plug&Play: паспорт для модели массы (расход, ступени и их π,
     частоты валов, выхлоп турбины, камера, подшипники, способ изготовления корпусов) -> mass_model.size_from_passport;
     к блокам массы двигателя добавляется то, чего паспортный путь не видит (всё - оценка, помечено в строках):
       межтурбинная камера подогрева (блок камеры при давлении и температуре газа перед ней), колёса IBC и улитки,
       сердцевины теплообменников (поверхность cycles.size_hx из блоков; риформер и парогенератор - cycles.hx_list,
       пересчёт по расходу), катализатор риформера (объёмная скорость, оценка), вентиляторы охладителей и рециркуляции
       (кг/кВт двигателя, оценка), газоходы вынесенных теплообменников (лист, оценка), плотность сплавов концепции
       вместо сплавов практики (колёса турбины, сопловые аппараты, осевые ступени компрессора);
     -> mass_model.evaluate (ротор + K·неподвижная часть + теплообменники + подшипники + рама), полоса - LOO и k_колл.
Выход (record): feasible, reject (коды) и reasons (тексты), КПД на клеммах и на валу, масса (полоса), кВт/кг, КПД на тонну,
то же с множителем готовности (0,5 + 0,5·TRL/9, решение владельца 15.11.2026), детали (эффективные), Plug&Play, TRL,
допуск пакета, частота вала, габариты, ниша.

  python3 process_eval.py --selfcheck       # опорные концепции: frame3 по четырём способам, рекуператор, ... (stdout JSON)
"""
import argparse
import copy
import json
import math
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
for _p in (HERE, os.path.join(HERE, "mass")):
    if _p not in sys.path:
        sys.path.insert(0, _p)
import process_space as PS  # noqa: E402
import blocks  # noqa: E402
from blocks import _track as T  # noqa: E402
from blocks import compat  # noqa: E402
from blocks.concept import normalize  # noqa: E402
import mass_model as MM  # noqa: E402

PRACTICE = "industrial"     # K промышленной практики (Solar): двигатели KV - литые/печатные промышленные; авиационная - нижняя граница
# «гибрид» для модели массы: холодная часть - точное литьё КНР (как cast_invest), горячая часть, камера и диффузор - печать LPBF
# (как lpbf); регистрируется в словаре модели массы при импорте (файл mass_model.py не меняется)
HYBRID_SET = dict(cold=("Al", "cast_invest"), inlet=("Al", "cast_invest"), warm=("steel", "cast_invest"),
                  comb_case=("IN718", "lpbf"), hot_case=("IN718", "lpbf"), diffuser=("IN625", "lpbf"), liner=("HastX", "lpbf"))
MM.PROCESS_SETS.setdefault("hybrid", HYBRID_SET)

HX_MASS_TYPE = {"foil": "psr", "lpbf": "plate_print"}     # технология сердцевины блоков (cycles.HX_TECH) -> тип модели массы
A_EVAL = dict(
    kg_per_kW_fan=6.0, duct_v=30.0, duct_L=2.0, ghsv_h=5000.0, cat_rho=900.0, ibc_rt_k=1.0,
    basis=dict(
        fan="вентилятор с асинхронным двигателем: 6 кг на кВт мощности привода (двигатели IE3 4-полюсные 50-100 кВт - около 5-7 кг/кВт, "
            "по памяти; корпус вентилятора включён) - оценка",
        duct="газоходы вынесенного теплообменника: по два на сторону (подвод и отвод), длина 2 м, скорость 30 м/с, лист стали со "
             "стенкой по давлению (mass_model.block_shell, способ sheet); группа «теплообменники» (без K практики) - оценка",
        reheat="межтурбинная камера: блок камеры модели массы (кольцевая, D_mean - средний диаметр турбины газогенератора) при давлении "
               "газа перед ней и температуре воздуха за компрессором (как основная камера: корпус охлаждается воздухом; горячий газ на "
               "входе - быстрая последовательная камера, площадь отсчёта по правилу основной) - оценка",
        ibc="IBC: колесо компрессора k·ρ·D2³ (D2 - блок hx.ibc_wheel: на валу двигателя или своё при φ > 0,15), диффузор и корпус 1,45 D2, "
            "колесо радиальной турбины того же диаметра и улитка; свой вал - вал по модели массы - оценка",
        reformer="риформер и парогенератор (сухой - риформер): сердцевина cycles.hx_list (технология концепции) при той же π, T41 и "
                 "блоках, поверхность × (расход концепции / расход перебора циклов); катализатор - объёмная скорость 5000 1/ч по "
                 "влажному сырью (паровой: CH4 + пар S/C; сухой: CH4 × (1 + R)) при 900 кг/м³ - оценка по памяти; водоподготовка, "
                 "бак и насос - вне массы агрегата (в ведомости деталей - есть)",
        density="колёса турбины, сопловые аппараты и осевые ступени компрессора: масса блока практики × ρ сплава концепции / ρ сплава "
                "практики (форма та же; для керамики Si3N4 / CMC форма была бы другой) - оценка; плотности - funnel/make_model.MAT",
        pgc_wave="волновой ротор (ВР): барабан с каналами не моделирован - добавлена масса, равная массе блока камеры (не подыгрывать) - оценка",
    ),
)

_CAL = None


def calibration():
    global _CAL
    if _CAL is None:
        _CAL = MM.load_calibration()
    return _CAL


def _rho(mat):
    if mat == "Ni":
        mat = "IN718"
    if mat in MM.RHO:
        return MM.RHO[mat]
    return T.MK().MAT[mat]["rho"]


# =========================================================================== паспорт для модели массы
def mass_spec(C, r):
    """паспорт концепции в форме mass_db (size_from_passport) по результату blocks.evaluate"""
    c, t, hx, b = C["compressor"], C["turbine"], C["hx"], C["build"]
    comp_rep = r["blocks"]["compressor"]["stages"]
    sz = r["sizing"]["engine"]
    groups = []
    ax = [s for s in comp_rep if s["type"] == "О"]
    if ax:
        groups.append(dict(type="A", n=len(ax), pi=[s["pi"] for s in ax]))
    first_c = True
    for s in comp_rep:
        if s["type"] == "О":
            continue
        sides = 2 if (c["inlet"] == "двух" and first_c and not ax) else 1
        groups.append(dict(type="C", sides=sides, pi=s["pi"], mat=("IN718" if c["mat"] == "Ni" else c["mat"])))
        first_c = False
    rotors = 2 if C["shaft"]["drive"] == "свободная СТ" else 1
    free = rotors == 2
    n_eng = sz["n_rpm"]
    pkg_n = r["package"]["n_pkg"]
    spec = dict(id="concept", name=r["name"][:80], practice=PRACTICE, process=C["_mass_process"],
                P_kW=r["P_shaft_kW"], eta_pct=r["eff_el"], m_air_kg_s=r["G_air"], PR=c["pi"],
                T_exh_C=r["T_turb_exit_C"], comp=groups,
                turb_gg=dict(type="A" if t["type"] == "О" else "R", n=int(r["nt"])),
                turb_pt=(dict(type="A", n=int(r["npt"])) if free else None), rpm_gg=n_eng, rpm_pt=(pkg_n if free else None),
                combustor=dict(kind="annular" if C["combustor"]["chamber"] == "кольцевая" else "silo", n_inj=int(sz["n_inj"]),
                               integral=(b["case"] == "оболочка")),
                bearings=dict(kind="magnetic" if b["bear"] == "магнит" else "oil_pkg", n_radial=2 * rotors, n_thrust=rotors),
                frame=dict(L=None), hx=[])
    return spec


def _hx_row(name, s, tech):
    return dict(name=name, type=HX_MASS_TYPE[tech], A_m2=float(s["A_m2"]))


def _duct_rows(name, m, p_bar, T_C, n=2):
    """газоходы: n труб длиной duct_L, диаметр по скорости duct_v (A_EVAL basis.duct)"""
    rho = p_bar * 1e5 / (287.0 * (T_C + 273.15))
    D = math.sqrt(4 * m / (rho * A_EVAL["duct_v"]) / math.pi)
    row = MM.block_shell(f"газоходы: {name}", D, A_EVAL["duct_L"] * n, p_bar, T_C, "steel", "sheet", group="hx")
    row["basis"] = A_EVAL["basis"]["duct"]
    return row


def extra_rows(C, r, geom, G):
    """блоки массы сверх паспортного пути (A_EVAL basis); возврат: строки rows двигателя (static/rotor) и список hx"""
    rows, hxl, notes = [], [], []
    hx_tech = r["sizing"]["hx_tech"]
    szh = r["sizing"]["hx"]
    blk = r["blocks"]
    # ---- сердцевины теплообменников блоков
    for name, s in szh.items():
        if s and s.get("ok"):
            hxl.append(_hx_row(name, s, hx_tech))
    comp = blk["compressor"]
    p2, T2, T3 = comp.get("p2_bar"), comp.get("T2_C"), r.get("T3_C")
    if C["hx"]["recup"] and szh.get("рекуператор"):
        if not C["hx"].get("recup_integral"):
            rows.append(_duct_rows("рекуператор, воздух (подвод)", G, p2, T2, n=1))
            rows.append(_duct_rows("рекуператор, воздух (отвод)", G, p2, T3, n=1))
        rows.append(_duct_rows("рекуператор, газ", r["W_exh"], 1.05, r["T_turb_exit_C"], n=2))
    ic = comp.get("ic")
    if ic:
        rows.append(_duct_rows("промохладитель, воздух", G, ic.get("p_bar", p2 / 2.0) if ic.get("p_bar") else max(p2 / 2.5, 1.5),
                               ic["T_in_C"], n=2))
        wf = ic["W_fan_kJ_per_kg_air"] * G
        rows.append(MM._row("вентилятор промохладителя с двигателем", "hx", A_EVAL["kg_per_kW_fan"] * wf, A_EVAL["basis"]["fan"],
                            W_fan_kW=round(wf, 1)))
    # ---- межтурбинная камера (подогрев)
    if C["hx"]["reheat"]:
        sec = blk["turbine"].get("sections") or []
        if sec:
            s0 = sec[0]
            Dm = geom["dims"].get("turb_gg", {}).get("Dm") or geom["dims"].get("turb_gg", {}).get("D") or 0.5
            crow, _ = MM.block_combustor(dict(kind="annular", m_air=G, p_bar=s0["p_out_bar"], T_C=T2 if T3 is None else min(T2, T3),
                                              D_mean=Dm, n_inj=int(r["sizing"]["engine"]["n_inj"])), PRACTICE)
            for x in crow:
                x["name"] = "подогрев: " + x["name"]
                x["basis"] = A_EVAL["basis"]["reheat"]
            rows += crow
    # ---- IBC
    if C["hx"]["ibc"]:
        w = blk["hx"].get("ibc_wheel")
        res = blk["hx"].get("ibc_result") or {}
        if w:
            D2 = w["D2_opt"] if w["own_shaft"] else w["D2_shaft"]
            for x in (MM.block_wheel_cf("IBC: колесо компрессора", D2, "IN718" if w else "Ti64", 1),
                      MM.block_cf_static("IBC: диффузор", D2, "steel"),
                      MM.block_shell("IBC: корпус компрессора", 1.45 * D2, 0.55 * D2, 1.05, 200.0, "steel", "cast_invest"),
                      MM.block_wheel_rt("IBC: колесо турбины", A_EVAL["ibc_rt_k"] * D2, "K418"),
                      MM.block_shell("IBC: улитка турбины", 1.6 * D2, 0.5 * D2, 1.05, 450.0, "steel", "cast_invest")):
                x["basis"] = A_EVAL["basis"]["ibc"]
                rows.append(x)
            if w["own_shaft"]:
                rows.append(MM.block_shaft("IBC: свой вал", 0.15 * 1.6 * D2, 1.2 * D2))
        if res.get("W_fan"):
            wf = res["W_fan"] * G / 1e3
            rows.append(MM._row("вентилятор охладителя IBC с двигателем", "hx", A_EVAL["kg_per_kW_fan"] * wf, A_EVAL["basis"]["fan"],
                                W_fan_kW=round(wf, 1)))
        rows.append(_duct_rows("охладитель IBC, газ", r["W_exh"], 0.9, r["T_turb_exit_C"], n=2))
    # ---- риформер
    prep = C["fuel"]["prep"]
    if prep in ("паровой риформинг", "сухой риформинг"):
        CY = T.CY()
        tcr = "steam" if prep == "паровой риформинг" else "egr"
        from blocks import fuel as FU
        Cn = normalize(PS.blocks_concept(C))
        prm = dict(T41_C=float(Cn["turbine"]["t41"]), x_rh=float(Cn["hx"]["x_rh"]), sc=float(Cn["fuel"]["sc"]), R_egr=float(Cn["fuel"]["R_egr"]))
        res = CY.evaluate(FU.cycles_topo(Cn, tcr), float(Cn["compressor"]["pi"]), prm, detail=True)
        k = G / float(res["W_air"]) if res.get("ok") else 0.0
        found = False
        if res.get("ok"):
            for s in CY.hx_list(res, hx_tech):
                if s.get("name") in ("риформер + парогенератор", "риформер"):
                    hxl.append(dict(name=s["name"], type=HX_MASS_TYPE[hx_tech], A_m2=float(s["A_m2"]) * k))
                    found = True
            if tcr == "egr" and res.get("W_rc_kW"):
                wrc = float(res["W_rc_kW"]) * k
                rows.append(MM._row("вентилятор рециркуляции выхлопа с двигателем", "hx", A_EVAL["kg_per_kW_fan"] * wrc,
                                    A_EVAL["basis"]["fan"], W_kW=round(wrc, 1)))
        if not found:
            notes.append("риформер: cycles.hx_list не дал сердцевину - масса сердцевины не учтена")
        ch4 = r["fuel_kg_s"]
        nm3h_ch4 = ch4 / 16.04 * 22.414 * 3600
        feed = nm3h_ch4 * (1 + (Cn["fuel"]["sc"] if tcr == "steam" else Cn["fuel"]["R_egr"]))
        V = feed / A_EVAL["ghsv_h"]
        rows.append(MM._row("катализатор риформера", "hx", V * A_EVAL["cat_rho"], A_EVAL["basis"]["reformer"], V_m3=round(V, 3)))
    return rows, hxl, notes


def _rescale_density(geom, C):
    """масса блоков ротора и аппаратов по плотности сплавов концепции (A_EVAL basis.density)"""
    t, c = C["turbine"], C["compressor"]
    pr = dict(MM.PRACTICE[PRACTICE])
    if C["_mass_process"] and MM.PROCESS_SETS.get(C["_mass_process"]):
        pr.update(MM.PROCESS_SETS[C["_mass_process"]])
    tw, hot, cm = _rho(t["mat"]), _rho(t["hot_mat"]), _rho(c["mat"])
    for row in geom["rows"]:
        n = row["name"]
        k = None
        if n.startswith("турбина газогенератора: диски") or n.startswith("силовая турбина: диски"):
            k = tw / _rho(pr["turb_rotor"]) if n.startswith("турбина газогенератора") else None
        elif n.startswith("колесо радиальной турбины"):
            k = tw / _rho(pr["vanes_t"])
        elif "сопловые аппараты" in n and n.startswith("турбина газогенератора") or n.startswith("сопловой аппарат радиальной"):
            k = hot / _rho(pr["vanes_t"])
        elif n.startswith("компрессор осевой") and "диски" in n:
            k = cm / _rho(pr["comp_rotor"])
        if k is not None and abs(k - 1.0) > 1e-9:
            row["mass_kg"] *= k
            row.setdefault("detail", {})["rho_scale"] = round(k, 4)
            row["basis"] = row["basis"] + "; " + A_EVAL["basis"]["density"]


def mass_of(C, r, detail=False):
    """масса агрегата-заменителя по результату блоков. Возврат: dict mass_model.evaluate + spec, notes"""
    spec = mass_spec(C, r)
    geom = MM.size_from_passport(spec)
    _rescale_density(geom, C)
    rows, hxl, notes = extra_rows(C, r, geom, r["G_air"])
    if C["combustor"]["scheme"] == "ВР":
        m_c = sum(x["mass_kg"] for x in geom["rows"] if x["name"].startswith("камера"))
        rows.append(MM._row("волновой ротор: барабан и окна (как камера)", "static", m_c, A_EVAL["basis"]["pgc_wave"]))
    geom["rows"] += rows
    geom["hx"] = hxl
    ev = MM.evaluate(geom, cal=calibration())
    # строки группы «теплообменники», добавленные сюда (газоходы, вентиляторы, катализатор, вентилятор рециркуляции), входят в итог
    # mass_model.evaluate, но не в его полосу (там - только сердцевины geom["hx"]): добавляются к обеим границам (с рамой)
    extra_hx = sum(x["mass_kg"] for x in rows if x["group"] == "hx")
    k_fr = 1.0 + (MM.A["cradle_k"] if geom.get("frame") else 0.0)
    lo_kg, hi_kg = ev["total_lo_kg"] + extra_hx * k_fr, ev["total_hi_kg"] + extra_hx * k_fr
    eta = r["eff_el"]
    ev["kpd_per_t_range"] = [round(MM.kpd_per_tonne(eta, hi_kg), 3), round(MM.kpd_per_tonne(eta, lo_kg), 3)]
    out = dict(total_kg=ev["total_kg"], total_lo_kg=round(lo_kg, 1), total_hi_kg=round(hi_kg, 1), engine_kg=ev["engine_kg"],
               hx_kg=ev["hx_kg"], bearings_extra_kg=ev["bearings_extra_kg"], frame_kg=ev["frame_kg"], kW_per_kg=ev.get("kW_per_kg"),
               kpd_per_t=ev.get("kpd_per_t"), kpd_per_t_range=ev.get("kpd_per_t_range"), K=ev["K"], groups_raw=ev["groups_raw"],
               dims=dict(L_engine=geom["dims"].get("L_engine"), D_case_max=geom["dims"].get("D_case_max"),
                         comp=[geom["dims"].get(f"comp{k}") for k in range(1, 5) if geom["dims"].get(f"comp{k}")],
                         turb_gg=geom["dims"].get("turb_gg"), turb_lp=geom["dims"].get("turb_lp")),
               notes=notes)
    if detail:
        out["rows"] = sorted(ev["rows"], key=lambda x: -x["mass_kg"])
        out["spec"] = spec
    return out


# =========================================================================== оценка
STATIC_ONLY_CODES = ()


def _static_fail(C):
    """быстрый отказ: жёсткие статические правила совместимости (без цикла)"""
    Cn = normalize(PS.blocks_concept(C))
    tgt = T.package_target(Cn["kv"])
    R = compat.check_static(Cn, tgt)
    bad = [x for x in R if x.get("hard") and not x["ok"]]
    return bad, R


def _compact_reasons(rules, limit=6):
    return [dict(code=x["code"], text=str(x.get("text", ""))[:160], val=x.get("val"), lim=x.get("lim")) for x in rules[:limit]]


def evaluate(x=None, kv="KV-3", concept=None, detail=False, uid="", mass_always=False):
    """геном x (или готовая концепция process_space.decode) -> запись оценки. Никаких исключений наружу: сбой модели -
    недопустимо с причиной. mass_always - масса и у недопустимой по блокам (для сверок; feasible не меняется)."""
    t0 = time.process_time()
    C = concept if concept is not None else PS.decode(x, kv)
    cell = PS.descriptors(C)
    rec = dict(uid=uid, kv=C["kv"], cell=list(cell), cell_label=PS.cell_label(cell), label=PS.summary(C), repairs=C["_repairs"],
               feasible=False, fitness=None)
    try:
        bad, _ = _static_fail(C)
    except (ValueError, KeyError) as e:
        rec.update(stage="схема", reject=["schema"], reasons=[dict(code="schema", text=str(e)[:160])], violation=1.0)
        rec["cpu_s"] = round(time.process_time() - t0, 5)
        return rec
    if bad:
        rec.update(stage="совместимость", reject=[b["code"] for b in bad], reasons=_compact_reasons(bad), violation=T.FS().violation(bad))
        rec["cpu_s"] = round(time.process_time() - t0, 5)
        return rec
    try:
        r = blocks.evaluate(PS.blocks_concept(C))
    except Exception as e:  # noqa: BLE001 - модель трека не должна ронять перебор: отказ с причиной
        rec.update(stage="блоки", reject=["blocks_error"], reasons=[dict(code="blocks_error", text=f"{type(e).__name__}: {e}"[:160])],
                   violation=10.0)
        rec["cpu_s"] = round(time.process_time() - t0, 5)
        return rec
    rec.update(eff_shaft=r.get("eff_shaft"), eff_el=r.get("eff_el"), G_air=r.get("G_air"), EGT_C=r.get("EGT_C"),
               violation=r.get("violation"), warnings=r.get("warnings", []), nt=r.get("nt"), npt=r.get("npt"))
    if r.get("interface"):
        it = r["interface"]
        rec["interface"] = dict(dG_pct=it["dG_pct"], dWexh_pct=it["dWexh_pct"], dEGT=it["dEGT"], P_margin_at_G_cap_pct=it["P_margin_at_G_cap_pct"])
    if r.get("parts"):
        rec.update(parts=r["parts"]["n_eff_total"], parts_engine=r["parts"]["engine"]["parts"], parts_added=r["parts"]["aux"]["added"],
                   pnp=r["plugplay"]["score"], pnp_blocked=r["plugplay"]["blocked"], pnp_cycles=r["plugplay"]["pnp_cycles"],
                   trl=r["trl"]["min"], k_trl=r["trl"]["k_trl"], service=r["service"]["score"], make_score=r["make"]["score"],
                   n_rpm=r["sizing"]["engine"]["n_rpm"], funnel_score=r["funnel_score"]["score"])
    bad = [x for x in r["rules"] if x.get("hard") and not x["ok"] and x.get("group") not in blocks.evaluator.RULE_GROUP_SECONDARY]
    if not r["feasible"]:
        rec.update(stage="блоки", reject=r["reject"], reasons=_compact_reasons(bad))
        if mass_always and r.get("sizing"):
            try:
                m = mass_of(C, r, detail=detail)
                rec.update(mass_t=round(m["total_kg"] / 1e3, 4), mass_lo_t=round(m["total_lo_kg"] / 1e3, 4),
                           mass_hi_t=round(m["total_hi_kg"] / 1e3, 4), engine_t=round(m["engine_kg"] / 1e3, 4),
                           hx_t=round(m["hx_kg"] / 1e3, 4), frame_t=round(m["frame_kg"] / 1e3, 4), kW_per_kg=m["kW_per_kg"],
                           kpd_per_t_if_feasible=m["kpd_per_t"])
                if detail:
                    rec["mass"] = m
            except Exception as e:  # noqa: BLE001
                rec["mass_error"] = f"{type(e).__name__}: {e}"[:160]
        rec["cpu_s"] = round(time.process_time() - t0, 5)
        if detail:
            rec["blocks"] = r
        return rec
    if r["plugplay"]["blocked"]:
        rec.update(stage="Plug&Play", reject=["pnp_blocked"], violation=1.0,
                   reasons=[dict(code="pnp_blocked", text="Plug&Play: «невозможно» по " + ", ".join(r["plugplay"]["blocked"]))])
        rec["cpu_s"] = round(time.process_time() - t0, 5)
        return rec
    try:
        m = mass_of(C, r, detail=detail)
    except Exception as e:  # noqa: BLE001
        rec.update(stage="масса", reject=["mass_error"], reasons=[dict(code="mass_error", text=f"{type(e).__name__}: {e}"[:160])],
                   violation=1.0)
        rec["cpu_s"] = round(time.process_time() - t0, 5)
        return rec
    if not (m["total_kg"] and m["total_kg"] > 0 and math.isfinite(m["total_kg"])):
        rec.update(stage="масса", reject=["mass_nan"], reasons=[dict(code="mass_nan", text="масса не определена")], violation=1.0)
        rec["cpu_s"] = round(time.process_time() - t0, 5)
        return rec
    kpt = m["kpd_per_t"]
    rec.update(feasible=True, stage="допустима", reject=[], reasons=[], violation=0.0, fitness=kpt,
               mass_t=round(m["total_kg"] / 1e3, 4), mass_lo_t=round(m["total_lo_kg"] / 1e3, 4), mass_hi_t=round(m["total_hi_kg"] / 1e3, 4),
               engine_t=round(m["engine_kg"] / 1e3, 4), hx_t=round(m["hx_kg"] / 1e3, 4), frame_t=round(m["frame_kg"] / 1e3, 4),
               kW_per_kg=m["kW_per_kg"], kpd_per_t=kpt, kpd_per_t_range=m["kpd_per_t_range"],
               fit_trl=round(kpt * rec["k_trl"], 4), L_engine=m["dims"]["L_engine"], D_case_max=m["dims"]["D_case_max"],
               mass_notes=m["notes"])
    rec["cpu_s"] = round(time.process_time() - t0, 5)
    if detail:
        rec["mass"] = m
        rec["blocks"] = r
        rec["concept"] = PS.blocks_concept(C)
    return rec


# =========================================================================== самопроверка
def reference_concepts():
    """опорные концепции: точка frame3 (блоки reference.frame3_concept) в четырёх семьях изготовления; та же с рекуператором и
    подогревом (как пример модели массы); быстроходный и радиальный варианты"""
    from blocks import reference as REF
    base = REF.frame3_concept()
    out = []
    for mk in PS.MAKE_KEYS:
        C = copy.deepcopy(base)
        M = PS.MAKE[mk]
        C["compressor"]["make"] = M["comp_make"]
        tmat = "SC" if M["turb_make"] == "литьё" else "Haynes282"
        tl = T.MK().t41_limit("О", "U", tmat, M["turb_make"])[0]
        C["turbine"].update(make=M["turb_make"], hot_make=M["hot_make"], mat=tmat, hot_mat="K438" if M["hot_make"] == "литьё" else "Haynes282",
                            t41=min(1000.0, tl))
        C["build"].update(comb_make=M["comb_make"], case_make=M["case_make"], case="разъём" if mk == "традиционная" else "картридж")
        C.setdefault("fuel", dict(prep="нет"))
        C.setdefault("hx", dict(recup=0.0, reheat=False, ibc=False))
        C.update(_make=mk, _mass_process=M["mass_process"], _repairs=[], name=f"frame3, {mk}")
        out.append(C)
    C = copy.deepcopy(out[0])
    C["compressor"]["pi"] = 13.5
    C["hx"] = dict(recup=0.90, recup_make="фольга", reheat=True, x_rh=0.35, ibc=False)
    C["name"] = "frame3 + рекуператор 0,90 (фольга) + подогрев, традиционная"
    out.append(C)
    C = copy.deepcopy(out[3])
    C["hx"] = dict(recup=0.90, recup_make="фольга", reheat=True, x_rh=0.35, ibc=False)
    C["name"] = "frame3 + рекуператор 0,90 (фольга) + подогрев, печать"
    out.append(C)
    return out


def _full(C):
    """концепция со всеми полями (DEFAULT блоков) и служебными полями process_space"""
    Cn = normalize({k: v for k, v in C.items() if not k.startswith("_") and k != "name"})
    Cn.update({k: v for k, v in C.items() if k.startswith("_")})
    Cn["name"] = C.get("name")
    Cn.pop("rotors", None)
    return Cn


def selfcheck():
    out = []
    for C in reference_concepts():
        Cf = _full(C)
        rec = evaluate(concept=Cf, detail=True, mass_always=True)
        row = {k: rec.get(k) for k in ("label", "feasible", "reject", "eff_el", "mass_t", "mass_lo_t", "mass_hi_t", "kW_per_kg", "kpd_per_t",
                                       "kpd_per_t_if_feasible", "engine_t", "hx_t", "frame_t", "parts", "pnp", "trl", "G_air", "cpu_s")}
        row["name"] = C.get("name")
        if rec.get("mass"):
            row["top_rows"] = [(x["name"], round(x["mass_kg"])) for x in rec["mass"]["rows"][:8]]
        out.append(row)
    return out


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--selfcheck", action="store_true")
    ap.add_argument("--index", type=int, default=None, help="оценить концепцию генератора с этим индексом")
    ap.add_argument("--kv", default="KV-3")
    a = ap.parse_args(argv)
    if a.selfcheck:
        print(json.dumps(selfcheck(), ensure_ascii=False, indent=1, default=str))
    if a.index is not None:
        rec = evaluate(PS.genome_at(a.index, a.kv), a.kv, detail=False, uid=str(a.index))
        print(json.dumps(rec, ensure_ascii=False, indent=1, default=str))
    return 0


if __name__ == "__main__":
    sys.exit(main())
