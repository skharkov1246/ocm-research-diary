"""Плотный перебор термодинамических циклов двигателя-заменителя KV-3 (трек «Простая литая ГТУ-заменитель»,
issue #59, направление «organic», задача владельца 14.11.2026).

Что считает (уровень облика: 0D-цикл + 1D-размер теплообменников + равновесие и кинетика риформинга; всё — оценка):
  * систематический перебор топологий на одном валу: простой цикл; рекуператор (степень регенерации 0,70–0,92);
    промежуточное охлаждение воздух-воздух с вентилятором (без воды); повторный подогрев (межтурбинная камера);
    обращённый цикл Брайтона на выхлопе (IBC); горение с ростом давления (идеализированная модель
    «эквивалентного располагаемого давления», доля реализации и штраф КПД турбины); химическая рекуперация —
    подогрев топлива, паровой риформинг метана теплом выхлопа (нужна вода), риформинг с рециркуляцией части
    выхлопа (без воды) — и все их сочетания (6 × 2 × 2 × 2 × 3 × 4 = 576 топологий);
  * для каждой — оптимум степени сжатия при T41 1000 °C (класс U, без охлаждения лопаток) и вариант «как есть»
    (π 13,5, компрессор KV-3 без переделки); расход под мощность Taurus 60 GS 5922 кВт на валу;
  * КПД на валу и на клеммах, удельная работа, расход и температура выхлопа против мягкого допуска пакета
    (+20 % расхода, +100 °C), размер и масса теплообменников (ламинарные каналы, постоянные Nu и f·Re
    посчитаны здесь же конечными разностями), число печатных модулей по камере принтера, добавочные детали
    и узлы, Plug&Play, готовность (TRL — оценка), ремонтопригодность;
  * эксергия: разрушение по узлам двумя способами (баланс эксергии узла и Гюи — Стодола T0·S_gen);
  * риформинг: равновесие Cantera, кинетика Xu — Froment (Ni/MgAl2O4) в реакторе идеального вытеснения
    с константами равновесия из термодинамики Cantera, активность углерода (сажеобразование), прирост
    теплотворности топлива, снижение разрушения эксергии в камере, пинч по составным кривым;
  * фронт Парето «КПД — добавочные детали — объём»; разветвлённый коллектор по закону Мюррея для рекуператора;
    компромисс «степень регенерации — потеря давления» (минимум производства энтропии по Бежану).
  * сверки: простой цикл — с gtu_cycle.cycle в расчётной точке frame3; рекуператор — с gtu_cycle.recup.

Запуск: python3 cycles.py  (около 10 мин на 4 ядрах) -> cycles_results.json, cycles.csv, PNG рядом со скриптом.
AWS не вызывается. Вход: ../../frame3/frame3_point_results.json, ../../gtu_family_results.json,
../../fleet/gtu_fleet_taurus60.json, ../../fleet/gtu_fleet_mercury50.json, ../../gtu_explore.py (LIM),
../../gtu_cycle.py (сверка)."""
import csv
import datetime
import json
import math
import os
import platform
import sys
import time
from concurrent.futures import ProcessPoolExecutor

import numpy as np
import cantera as ct
import scipy
from scipy.integrate import solve_ivp
from scipy.optimize import brentq, minimize_scalar

HERE = os.path.dirname(os.path.abspath(__file__))
GTU = os.path.abspath(os.path.join(HERE, "..", ".."))
REPO = os.path.abspath(os.path.join(GTU, "..", ".."))
SCRATCH = "/tmp/claude-0/-home-user/a62745e4-d2fd-56e0-974f-12746b1a2ae0/scratchpad/organic"
DATE = "2026-11-14"

F_FRAME3 = os.path.join(GTU, "frame3", "frame3_point_results.json")
F_FAMILY = os.path.join(GTU, "gtu_family_results.json")
F_T60 = os.path.join(GTU, "fleet", "gtu_fleet_taurus60.json")
F_M50 = os.path.join(GTU, "fleet", "gtu_fleet_mercury50.json")
F_EXPLORE = os.path.join(GTU, "gtu_explore.py")


def rel(p):
    return os.path.relpath(p, REPO)


# ======================================================================== вход
def load_inputs():
    """Расчётная точка KV-3 и допуски пакета — из файлов трека (каждое значение с файлом и полем)."""
    f3 = json.load(open(F_FRAME3, encoding="utf-8"))
    cyc = f3["cycle"]["frame3"]
    inp = {x["field"]: x["used"] for x in f3["inputs"]}
    soft = inp["slots[model=Taurus 60, variant=генератор].мягкий"]
    fam = json.load(open(F_FAMILY, encoding="utf-8"))
    slot = [s for s in fam["slots"] if s["model"] == "Taurus 60" and s["variant"] == "генератор"][0]
    t60 = [r for r in json.load(open(F_T60, encoding="utf-8"))["records"] if r["variant"] == "GS"][0]
    m50 = json.load(open(F_M50, encoding="utf-8"))["records"][0]
    # LIM['T_rec_max'] из gtu_explore.py читаем текстом, чтобы не делать chdir импортом
    lim_rec = None
    for line in open(F_EXPLORE, encoding="utf-8"):
        if "T_rec_max" in line and "=" in line:
            lim_rec = float(line.split("=")[1].split(",")[0])
            break
    d = dict(
        PR0=float(soft["PR"]), T41_C=float(soft["T41_C"]), G0=float(soft["G"]), n_out=float(soft["n_out"]),
        eta_pc0=float(cyc["eta_pc"]), eta_t=float(cyc["eta_t"]), P_kW_frame3=float(cyc["P_shaft_kW"]),
        eff_frame3=float(cyc["eff_shaft"]), EGT_frame3=float(cyc["EGT_C"]),
        P_shaft_kW=float(slot["P_shaft_kW"]), P_el_kW=float(slot["P_kW"]), G_air0=float(slot["G_air"]),
        EGT0=float(slot["EGT0"]), eta0=float(slot["eta0"]), pkg=slot["pkg"],
        W_exh0=float(t60["exh_flow_kg_s"]), t60_pkg_kg=float(t60["pkg_kg"]),
        m50=dict(P_kW=float(m50["power_kW"]), eff=float(m50["eff_pct"]), PR=float(m50["PR"]),
                 W_exh=float(m50["exh_flow_kg_s"]), EGT=float(m50["exh_temp_C"]), recup=m50.get("recuperator"),
                 pkg=[float(m50["pkg_L_m"]), float(m50["pkg_W_m"]), float(m50["pkg_H_m"])], pkg_kg=float(m50["pkg_kg"])),
        T_rec_max=lim_rec,
    )
    sources = [
        dict(file=rel(F_FRAME3), field="cycle.frame3 (eta_pc, eta_t, P_shaft_kW, eff_shaft, EGT_C)",
             used=dict(eta_pc=d["eta_pc0"], eta_t=d["eta_t"], P_kW=d["P_kW_frame3"], eff=d["eff_frame3"], EGT_C=d["EGT_frame3"])),
        dict(file=rel(F_FRAME3), field="inputs[slots[model=Taurus 60, variant=генератор].мягкий]",
             used=dict(PR=d["PR0"], T41_C=d["T41_C"], G=d["G0"], n_out=d["n_out"])),
        dict(file=rel(F_FAMILY), field="slots[Taurus 60, генератор]: P_kW, P_shaft_kW, G_air, EGT0, eta0, pkg",
             used=dict(P_kW=d["P_el_kW"], P_shaft_kW=d["P_shaft_kW"], G_air=d["G_air0"], EGT0=d["EGT0"], eta0=d["eta0"], pkg=d["pkg"])),
        dict(file=rel(F_T60), field="records[GS]: exh_flow_kg_s, pkg_kg", used=dict(W_exh=d["W_exh0"], pkg_kg=d["t60_pkg_kg"])),
        dict(file=rel(F_M50), field="records[0]: power_kW, eff_pct, PR, exh_flow_kg_s, exh_temp_C, recuperator, pkg_*",
             used=d["m50"]),
        dict(file=rel(F_EXPLORE), field="LIM['T_rec_max']", used=d["T_rec_max"]),
        dict(file=rel(os.path.join(GTU, "gtu_cycle.py")), field="AIR, T0, P0, LHV, poly_process, cycle(...) по умолчанию",
             used=dict(dp_comb=0.04, p_exh_ratio=1.02, eta_b=0.995, eta_m=0.99, poly_steps=40, LHV=50.03e6)),
        dict(file=rel(os.path.join(GTU, "gtu_family.py")), field="SCEN['мягкий'] (alpha 0.20, dT 100), x_leak 0.015",
             used=dict(alpha=0.2, dT=100, x_leak=0.015)),
    ]
    return d, sources


IN, IN_SOURCES = load_inputs()

# ======================================================================== принятые значения (оценка)
AIR = "O2:0.2095,N2:0.7808,AR:0.0093,CO2:0.0004"   # как gtu_cycle.AIR
P0, T0 = 101325.0, 288.15
LHV_CH4 = 50.03e6
RU = ct.gas_constant                                # Дж/(кмоль·К)

A = dict(
    T41_C=IN["T41_C"], eta_t=IN["eta_t"], eta_pc0=IN["eta_pc0"], PR0=IN["PR0"],
    x_leak=0.015, eta_m=0.99, eta_b=0.995, dp_comb=0.04, p_exh=1.02, poly_n=40,
    eta_gen=0.972, eta_gb=0.985,
    P_shaft_kW=IN["P_shaft_kW"],
    alpha=0.20, dT_exh=100.0,
    G_cap=1.2 * IN["G_air0"], W_exh_cap=1.2 * IN["W_exh0"], EGT_cap=IN["EGT0"] + 100.0,
    PRs_1stage=5.0,                     # одна центробежная ступень до π 5 (ступень 1 KV-3: 4,49 — frame3)
    psi_t=1.5, U_t=420.0,               # число ступеней турбины: работа ступени ≤ ψ·U² (U — LIM['U_ax'])
    # рекуператор
    dp_ra=0.025, dp_rg=0.035,           # как gtu_cycle.recup по умолчанию
    T_rec_max=IN["T_rec_max"],          # газ на входе в рекуператор, °C (LIM gtu_explore)
    # промохлаждение воздух-воздух
    eps_ic=0.80, dp_ic=0.02, Cr_ic=1.5, dp_fan=500.0, eta_fan=0.65, eta_motor=0.95,
    # повторный подогрев
    dp_rh=0.03, x_rh=(0.35, 0.5, 0.65),
    # IBC
    eta_ibc_t=0.86, eta_ibc_c=0.82, eps_ibc=0.85, dp_ibc=0.02, plow=(0.25, 0.3, 0.35, 0.4, 0.5, 0.6, 0.7),
    phi_max=0.15,                       # наибольший коэффициент расхода центробежного колеса на общем валу
    # горение с ростом давления
    pgc_f_real=0.25, pgc_deta_t=0.02,
    # химическая рекуперация
    p_fuel_k=1.25, dT_app=40.0, approach_eq=15.0, T_ref_max=850.0, T_fh_max=400.0, dp_ref_g=0.02, dp_fh_g=0.002, dp_ref_f=0.05,
    pinch=20.0, sc=(1.5, 2.0, 3.0), R_egr=(0.5, 1.0, 2.0, 4.0), eta_rc=0.78, eta_pump=0.6, T_rc_cool=80.0,
    cat_loading=0.05,                   # кг катализатора (с носителем) на м² стенки — покрытие стенок риформера, оценка
    eta_eff_cat=1.0,                    # фактор эффективности тонкого покрытия (для гранул ~0,05 — оценка)
    xf_h2_seed=1e-3,                    # затравка H2/CH4 на входе (в уравнениях Xu — Froment p_H2 в знаменателе)
    # теплообменники: ламинарные каналы
    dp_core_share=0.7,                  # доля допустимой потери давления на сердцевину; остальное — коллекторы
    RH_env=0.60,                        # среда ISO 3977-2 для эксергии (как organic/exergy)
)

HX_TECH = {
    "foil": dict(name="первичная поверхность из гофрированной фольги (как Mercury 50), сварка/пайка",
                 Dh=1.5e-3, t=0.10e-3, rho=8440.0, k_w=11.0, k_struct=1.0, printer=None, fill=None,
                 note="фольга 0,1 мм (Mercury 50: 347 → Alloy 625 — по выдержке поиска); Dh и коэффициент корпуса — оценка"),
    "lpbf": dict(name="печатная сердцевина TPMS (гироид), LPBF, Inconel 625",
                 Dh=2.0e-3, t=0.35e-3, rho=8440.0, k_w=11.0, k_struct=0.6,
                 printer=dict(name="BLT-S1500", box_m=(1.5, 1.5, 1.2)), fill=0.4,
                 note="стенка 0,35 мм и канал 2 мм под вынос порошка — оценка (ISO/ASTM 52911-1: темы стенок и выноса "
                      "порошка — по выдержке поиска, значения не открыты); камера BLT-S1500 1500×1500×1200 мм — по выдержке поиска"),
}

# ======================================================================== газ (Cantera)
SP = ["N2", "O2", "AR", "CO2", "H2O", "CH4", "H2", "CO"]
I = {s: i for i, s in enumerate(SP)}
_G = None


def gas():
    global _G
    if _G is None:
        full = ct.Solution("gri30.yaml")
        _G = ct.Solution(thermo="ideal-gas", species=[full.species(n) for n in SP], transport_model="mixture-averaged")
    return _G


MW = gas().molecular_weights.copy()


def Y_of(X):
    g = gas()
    g.TPX = T0, P0, X
    return g.Y.copy()


Y_AIR = Y_of(AIR)
Y_CH4 = Y_of("CH4:1")

_W = ct.Water()
_W.TQ = T0, 0.0
H_L0, P_SAT0 = _W.enthalpy_mass, _W.P
_W.TQ = T0, 1.0
H_FG0 = _W.enthalpy_mass - H_L0


def p_sat(T):
    _W.TQ = min(max(T, 274.0), 640.0), 0.0
    return _W.P


_TS = {}


def T_sat(P):
    k = round(P, -2)
    if k not in _TS:
        _W.PQ = k, 0.0
        _TS[k] = _W.T
    return _TS[k]


class St:
    """Поток: m — кг/с на 1 кг/с воздуха на входе; liq — жидкая вода."""
    __slots__ = ("m", "T", "P", "Y", "liq")

    def __init__(self, m, T, P, Y, liq=False):
        self.m, self.T, self.P, self.Y, self.liq = float(m), float(T), float(P), np.asarray(Y, float), liq

    def cp(self, **kw):
        s = St(self.m, self.T, self.P, self.Y, self.liq)
        for k, v in kw.items():
            setattr(s, k, v)
        return s


def hs(s):
    """удельные энтальпия и энтропия потока"""
    if s.liq:
        g = gas()
        g.TPY = T0, P_SAT0, Y_of("H2O:1")
        h = g.enthalpy_mass - H_FG0 + (s.P - P_SAT0) / 999.1
        return h, g.entropy_mass - H_FG0 / T0
    g = gas()
    g.TPY = s.T, s.P, s.Y
    return g.enthalpy_mass, g.entropy_mass


def T_of_h(h, P, Y):
    g = gas()
    g.HPY = h, P, Y
    return g.T


def mix(streams, P):
    m = sum(s.m for s in streams)
    Y = sum(s.m * s.Y for s in streams) / m
    h = sum(s.m * hs(s)[0] for s in streams) / m
    return St(m, T_of_h(h, P, Y), P, Y)


# ---------------------------------------------------------------- эксергия (среда ISO 3977-2: 15 °C, 101,325 кПа, RH 60 %)
def _env():
    x_w = A["RH_env"] * P_SAT0 / P0
    g = gas()
    g.TPX = T0, P0, AIR
    X = g.X.copy() * (1 - x_w)
    X[I["H2O"]] = x_w
    gs = g.standard_gibbs_RT * RU * T0          # Дж/кмоль при 1 атм = P0
    gi = {s: gs[I[s]] for s in SP}
    e = np.zeros(len(SP))
    for s in ("N2", "O2", "AR", "CO2", "H2O"):
        e[I[s]] = -RU * T0 * math.log(X[I[s]])
    eO2, eCO2, eH2O = e[I["O2"]], e[I["CO2"]], e[I["H2O"]]
    e[I["CH4"]] = -(gi["CO2"] + 2 * gi["H2O"] - gi["CH4"] - 2 * gi["O2"]) + eCO2 + 2 * eH2O - 2 * eO2
    e[I["H2"]] = -(gi["H2O"] - gi["H2"] - 0.5 * gi["O2"]) + eH2O - 0.5 * eO2
    e[I["CO"]] = -(gi["CO2"] - gi["CO"] - 0.5 * gi["O2"]) + eCO2 - 0.5 * eO2
    e_liq = RU * T0 * math.log(P_SAT0 / (x_w * P0)) / MW[I["H2O"]]
    return X, e, e_liq


X_ENV, E_SP, E_LIQ = _env()


def ex(s):
    """удельная эксергия потока (физическая + химическая), Дж/кг"""
    if s.liq:
        return E_LIQ + (s.P - P_SAT0) / 999.1
    g = gas()
    g.TPY = T0, P0, s.Y
    h0, s0 = g.enthalpy_mass, g.entropy_mass
    X = g.X.copy()
    Mm = g.mean_molecular_weight
    g.TPY = s.T, s.P, s.Y
    eph = (g.enthalpy_mass - h0) - T0 * (g.entropy_mass - s0)
    nz = X > 1e-30
    ech = (np.dot(X, E_SP) + RU * T0 * np.sum(X[nz] * np.log(X[nz]))) / Mm
    return eph + ech


def lhv(Y):
    """низшая теплотворность смеси (на кг смеси) при T0: полное сгорание, вода — пар (Cantera)"""
    g = gas()
    g.TP = T0, P0
    hm = g.standard_enthalpies_RT * RU * T0
    n = Y / MW
    C = n[I["CO2"]] + n[I["CH4"]] + n[I["CO"]]
    H = 2 * n[I["H2O"]] + 4 * n[I["CH4"]] + 2 * n[I["H2"]]
    O = 2 * n[I["O2"]] + 2 * n[I["CO2"]] + n[I["H2O"]] + n[I["CO"]]
    o2add = (2 * C + H / 2 - O) / 2
    p = np.zeros(len(SP))
    p[I["N2"]], p[I["AR"]] = n[I["N2"]], n[I["AR"]]
    p[I["CO2"]], p[I["H2O"]] = C, H / 2
    return float(np.dot(n, hm) + o2add * hm[I["O2"]] - np.dot(p, hm))


LHV_CH4_CT = lhv(Y_CH4)


# ---------------------------------------------------------------- узлы цикла
def poly(s, p2, eta, compress, n=None):
    """политропный процесс малыми шагами, состав заморожен (как gtu_cycle.poly_process); меняет s, работа Дж/кг"""
    n = n or A["poly_n"]
    g = gas()
    g.TPY = s.T, s.P, s.Y
    h0 = g.enthalpy_mass
    r = (p2 / s.P) ** (1.0 / n)
    for _ in range(n):
        ss, h, p = g.entropy_mass, g.enthalpy_mass, g.P
        g.SP = ss, p * r
        dh = g.enthalpy_mass - h
        g.HPY = h + (dh / eta if compress else dh * eta), p * r, s.Y
    s.T, s.P = g.T, g.P
    return abs(g.enthalpy_mass - h0)


def eta_stage(PRs):
    """КПД (политропный) центробежной ступени KV-3: формула трека comp_stage + поправка, при которой
    двухступенчатый компрессор π 13,5 даёт эквивалентный КПД frame3 (0,8413)."""
    d = 0.880 - 0.010 * (math.sqrt(A["PR0"]) - 2.5) - A["eta_pc0"]
    return 0.880 - 0.010 * (PRs - 2.5) - d


def n_stages(PR):
    return 1 if PR <= A["PRs_1stage"] + 1e-9 else (2 if PR <= A["PRs_1stage"] ** 2 else 3)


def n_turb(w):
    return max(1, math.ceil(w / (A["psi_t"] * A["U_t"] ** 2) - 1e-9))


def burn(air, fuel_Y, fuel_h, T_out, P_out, lhv_f):
    """полное сгорание: подбирает расход топлива (кг на кг/с воздуха на входе), чтобы продукты имели T_out."""
    g = gas()
    h_a = hs(air)[0]
    na = air.m * air.Y / MW
    nf = fuel_Y / MW
    eb = A["eta_b"]

    def prod(mf):
        n = na + mf * nf
        C = n[I["CO2"]] + n[I["CH4"]] + n[I["CO"]]
        H = 2 * n[I["H2O"]] + 4 * n[I["CH4"]] + 2 * n[I["H2"]]
        O = 2 * n[I["O2"]] + 2 * n[I["CO2"]] + n[I["H2O"]] + n[I["CO"]]
        p = np.zeros(len(SP))
        p[I["N2"]], p[I["AR"]] = n[I["N2"]], n[I["AR"]]
        p[I["CO2"]], p[I["H2O"]] = C, H / 2
        p[I["O2"]] = (O - 2 * C - H / 2) / 2
        Yp = p * MW
        return Yp / Yp.sum(), p[I["O2"]]

    def resid(mf):
        Yp, _ = prod(mf)
        hm = (air.m * h_a + mf * fuel_h - (1 - eb) * mf * lhv_f) / (air.m + mf)
        g.HPY = hm, P_out, Yp
        return g.T - T_out

    # наибольший расход топлива, при котором ещё остаётся кислород
    C = nf[I["CO2"]] + nf[I["CH4"]] + nf[I["CO"]]
    H = 2 * nf[I["H2O"]] + 4 * nf[I["CH4"]] + 2 * nf[I["H2"]]
    O = 2 * nf[I["O2"]] + 2 * nf[I["CO2"]] + nf[I["H2O"]] + nf[I["CO"]]
    need = (2 * C + H / 2 - O) / 2
    mf_max = 0.97 * na[I["O2"]] / need
    if resid(mf_max) < 0:
        raise ValueError("не хватает кислорода до T_out")
    mf = brentq(resid, 1e-8, mf_max, xtol=1e-12, rtol=1e-12)
    Yp, _ = prod(mf)
    resid(mf)
    return mf, St(air.m + mf, g.T, P_out, Yp)


def pgc_gain(air, fuel_Y, fuel_h, mf, lhv_f):
    """идеальный предел горения с ростом давления: эквивалентное располагаемое давление (EAP, Kaemming & Paxson —
    подход; здесь простейший вид): продукты горения при постоянном объёме имеют энтропию s_cv; стационарный поток
    с той же полной энтальпией и той же энтропией имеет давление P_eap. Возвращает P_eap/p3."""
    g = gas()
    h_a = hs(air)[0]
    m = air.m + mf
    Ym = (air.m * air.Y + mf * fuel_Y) / m
    hm = (air.m * h_a + mf * fuel_h) / m
    g.HPY = hm, air.P, Ym
    u_r, v_r = g.int_energy_mass, 1.0 / g.density
    # продукты того же расхода топлива (полное сгорание)
    n = (air.m * air.Y + mf * fuel_Y) / MW
    C = n[I["CO2"]] + n[I["CH4"]] + n[I["CO"]]
    H = 2 * n[I["H2O"]] + 4 * n[I["CH4"]] + 2 * n[I["H2"]]
    O = 2 * n[I["O2"]] + 2 * n[I["CO2"]] + n[I["H2O"]] + n[I["CO"]]
    p = np.zeros(len(SP))
    p[I["N2"]], p[I["AR"]] = n[I["N2"]], n[I["AR"]]
    p[I["CO2"]], p[I["H2O"]], p[I["O2"]] = C, H / 2, (O - 2 * C - H / 2) / 2
    Yp = p * MW
    Yp /= Yp.sum()
    loss = (1 - A["eta_b"]) * mf * lhv_f / m
    g.UVY = u_r - loss, v_r, Yp
    s_cv, p_cv, T_cv = g.entropy_mass, g.P, g.T
    g.HPY = hm - loss, air.P, Yp
    R = RU / g.mean_molecular_weight
    P_eap = air.P * math.exp((g.entropy_mass - s_cv) / R)
    return P_eap / air.P, dict(p_cv_over_p3=p_cv / air.P, T_cv_K=T_cv)


def recuperate(air, gasin, eps, dpa, dpg):
    """противоток, эффективность по наименьшей теплоёмкости (энтальпийная запись)"""
    ha, _ = hs(air)
    hg, _ = hs(gasin)
    ha_max = hs(air.cp(T=gasin.T))[0]
    hg_min = hs(gasin.cp(T=air.T))[0]
    Q = eps * min(air.m * (ha_max - ha), gasin.m * (hg - hg_min))
    ao = St(air.m, T_of_h(ha + Q / air.m, air.P * (1 - dpa), air.Y), air.P * (1 - dpa), air.Y)
    go = St(gasin.m, T_of_h(hg - Q / gasin.m, gasin.P / (1 + dpg), gasin.Y), gasin.P / (1 + dpg), gasin.Y)
    return ao, go, Q


def cool_with_air(hot, eps, dp_rel, condense=False):
    """охладитель воздух-воздух с вентилятором: горячий поток охлаждается до T_h − eps·(T_h − T0);
    холодный — атмосферный воздух с теплоёмкостью Cr_ic × тепловой нагрузки на ΔT. Возвращает потоки и работу."""
    Tout = hot.T - eps * (hot.T - T0)
    Pout = hot.P * (1 - dp_rel)
    h_in, _ = hs(hot)
    outs = []
    if condense:
        g = gas()
        g.TPY = Tout, Pout, hot.Y
        X = g.X.copy()
        xs = p_sat(Tout) / Pout
        if X[I["H2O"]] > xs:
            n = hot.m * hot.Y / MW
            ndry = n.sum() - n[I["H2O"]]
            nw_v = xs * ndry / (1 - xs)
            ncond = n[I["H2O"]] - nw_v
            n2 = n.copy()
            n2[I["H2O"]] = nw_v
            Yv = n2 * MW
            mv = Yv.sum()
            vap = St(mv, Tout, Pout, Yv / mv)
            liq = St(ncond * MW[I["H2O"]], Tout, Pout, Y_of("H2O:1"), liq=True)
            liq_h = hs(St(1, T0, P_SAT0, Y_of("H2O:1"), liq=True))[0] + 4180.0 * (Tout - T0)
            Q = hot.m * h_in - vap.m * hs(vap)[0] - liq.m * liq_h
            outs = [vap, liq]
        else:
            outs = [St(hot.m, Tout, Pout, hot.Y)]
            Q = hot.m * (h_in - hs(outs[0])[0])
    else:
        outs = [St(hot.m, Tout, Pout, hot.Y)]
        Q = hot.m * (h_in - hs(outs[0])[0])
    cp_air = 1006.0
    mc = A["Cr_ic"] * Q / max(hot.T - T0, 1.0) / cp_air
    cold_in = St(mc, T0, P0, Y_AIR)
    hc = hs(cold_in)[0]
    cold_out = St(mc, T_of_h(hc + Q / mc, P0, Y_AIR), P0, Y_AIR)
    rho_c = P0 / (RU / 28.96 * T0)
    W_fan = mc * A["dp_fan"] / rho_c / (A["eta_fan"] * A["eta_motor"])
    return outs, cold_in, cold_out, Q, W_fan


# ---------------------------------------------------------------- риформинг: равновесие, сажа, кинетика
def equil(Y, h_or_T, P, mode="TP"):
    g = gas()
    if mode == "TP":
        g.TPY = h_or_T, P, Y
        g.equilibrate("TP")
    else:
        g.HPY = h_or_T, P, Y
        g.equilibrate("HP")
    return g.Y.copy(), g.T


_GR = None


def carbon_activity(Y, T, P):
    """активность углерода (графит) по разложению метана и по Будуару в газе состава Y при T, P; > 1 — сажа выгодна"""
    global _GR
    if _GR is None:
        _GR = ct.Solution("graphite.yaml")
    g = gas()
    g.TPY = T, P, Y
    pa = g.X * P / ct.one_atm
    g.TP = T, ct.one_atm          # standard_gibbs_RT в Cantera — при текущем давлении; нужны значения при 1 атм
    gs = g.standard_gibbs_RT.copy()
    _GR.TP = T, ct.one_atm
    gc = _GR.standard_gibbs_RT[0]
    KM = math.exp(-(2 * gs[I["H2"]] + gc - gs[I["CH4"]]))
    KB = math.exp(-(gc + gs[I["CO2"]] - 2 * gs[I["CO"]]))
    aM = KM * pa[I["CH4"]] / max(pa[I["H2"]], 1e-30) ** 2
    aB = KB * pa[I["CO"]] ** 2 / max(pa[I["CO2"]], 1e-30)
    return aM, aB


# Xu J., Froment G.F., AIChE J. 35 (1989) 88–96, Ni/MgAl2O4. Источник подтверждён выдержкой поиска;
# численные значения — по памяти (часто воспроизводимая таблица), с оригиналом не сверены.
XF = dict(k1=(4.225e15, 240.1e3), k2=(1.955e6, 67.13e3), k3=(1.020e15, 243.9e3),
          K_CO=(8.23e-5, -70.65e3), K_H2=(6.12e-9, -82.90e3), K_CH4=(6.65e-4, -38.28e3), K_H2O=(1.77e5, 88.68e3))


def K_eq_bar(T):
    """константы равновесия реакций Xu — Froment из термодинамики Cantera (GRI-Mech 3.0), давления в барах:
    1) CH4+H2O=CO+3H2 [бар²]; 2) CO+H2O=CO2+H2 [-]; 3) CH4+2H2O=CO2+4H2 [бар²]"""
    g = gas()
    g.TP = T, P0
    gs = g.standard_gibbs_RT
    a2b = (ct.one_atm / 1e5)
    K1 = math.exp(-(gs[I["CO"]] + 3 * gs[I["H2"]] - gs[I["CH4"]] - gs[I["H2O"]])) * a2b ** 2
    K2 = math.exp(-(gs[I["CO2"]] + gs[I["H2"]] - gs[I["CO"]] - gs[I["H2O"]]))
    return K1, K2, K1 * K2


def K_eq_xf_empirical(T):
    """эмпирические формулы, которые обычно приводят вместе с моделью Xu — Froment (по памяти, для самопроверки)"""
    K1 = math.exp(-26830.0 / T + 30.114)
    K2 = math.exp(4400.0 / T - 4.036)
    return K1, K2, K1 * K2


def xf_rates(F, T, P_bar):
    """скорости r1, r2, r3, моль/(кг кат·с); F — мольные расходы (по SP)"""
    Ft = F.sum()
    p = F / Ft * P_bar
    pH2 = max(p[I["H2"]], 1e-8)
    K1, K2, K3 = K_eq_bar(T)
    k = {n: XF[n][0] * math.exp(-XF[n][1] / (RU / 1e3 * T)) for n in ("k1", "k2", "k3")}
    K = {n: XF[n][0] * math.exp(-XF[n][1] / (RU / 1e3 * T)) for n in ("K_CO", "K_H2", "K_CH4", "K_H2O")}
    DEN = 1 + K["K_CO"] * p[I["CO"]] + K["K_H2"] * pH2 + K["K_CH4"] * p[I["CH4"]] + K["K_H2O"] * p[I["H2O"]] / pH2
    r1 = k["k1"] / pH2 ** 2.5 * (p[I["CH4"]] * p[I["H2O"]] - pH2 ** 3 * p[I["CO"]] / K1) / DEN ** 2
    r2 = k["k2"] / pH2 * (p[I["CO"]] * p[I["H2O"]] - pH2 * p[I["CO2"]] / K2) / DEN ** 2
    r3 = k["k3"] / pH2 ** 3.5 * (p[I["CH4"]] * p[I["H2O"]] ** 2 - pH2 ** 4 * p[I["CO2"]] / K3) / DEN ** 2
    return np.array([r1, r2, r3]) * 1000.0 / 3600.0      # кмоль/(кг·ч) -> моль/(кг·с)


NU_R = np.zeros((3, len(SP)))
for _j, _d in enumerate(({"CH4": -1, "H2O": -1, "CO": 1, "H2": 3}, {"CO": -1, "H2O": -1, "CO2": 1, "H2": 1},
                         {"CH4": -1, "H2O": -2, "CO2": 1, "H2": 4})):
    for _s, _v in _d.items():
        NU_R[_j, I[_s]] = _v


def xf_pfr(X_in, F_tot, T, P, W_max, eta_eff=1.0, n_out=200):
    """изотермический реактор идеального вытеснения; X_in — мольные доли (O2 уже сожжён), F_tot — моль/с.
    Возвращает массу катализатора, кг, и степень превращения CH4 по длине."""
    F0 = np.asarray(X_in, float) * F_tot
    F0[I["H2"]] = max(F0[I["H2"]], A["xf_h2_seed"] * F0[I["CH4"]])
    P_bar = P / 1e5

    def rhs(W, F):
        return NU_R.T @ xf_rates(np.maximum(F, 1e-12), T, P_bar) * eta_eff

    sol = solve_ivp(rhs, (0.0, W_max), F0, method="LSODA", rtol=1e-7, atol=1e-10 * F_tot,
                    t_eval=np.geomspace(W_max * 1e-6, W_max, n_out))
    X = 1 - sol.y[I["CH4"]] / F0[I["CH4"]]
    return sol.t, X, sol.y


# ---------------------------------------------------------------- топливо: подогрев, паровой риформинг, рециркуляция
def feed_steam(sc, P=P0):
    """смесь CH4 + вода (жидкая при T0, после насоса при давлении P) на кг смеси: состав, энтальпия, доля CH4"""
    mC, mW = MW[I["CH4"]], sc * MW[I["H2O"]]
    yC = mC / (mC + mW)
    Y = np.zeros(len(SP))
    Y[I["CH4"]], Y[I["H2O"]] = yC, 1 - yC
    h = yC * hs(St(1, T0, P, Y_CH4))[0] + (1 - yC) * hs(St(1, T0, P, Y_of("H2O:1"), liq=True))[0]
    return Y, h, yC


_CC = {}


def cold_curve_steam(sc, P, T_ref):
    """составная кривая холодной стороны парового риформинга на кг смеси: [(T, Q от входа)]; равновесие с подходом.
    Запоминается по (S/C, P с шагом 100 Па, T с шагом 0,01 К): зависит только от них."""
    key = (sc, round(P, -2), round(T_ref, 2))
    if key not in _CC:
        _CC[key] = _cold_curve_steam(sc, key[1], T_ref)
    return _CC[key]


def _cold_curve_steam(sc, P, T_ref):
    Y, h_in, yC = feed_steam(sc, P)
    Ts = T_sat(P)
    yW = 1 - yC
    hC = lambda T: hs(St(1, T, P, Y_CH4))[0]
    hW_liq = lambda T: hs(St(1, T0, P, Y_of("H2O:1"), liq=True))[0] + 4180.0 * (T - T0) * (1 + 0.0007 * (T - T0))
    pts = []
    for T in np.linspace(T0, Ts, 8):
        pts.append((T, yC * hC(T) + yW * hW_liq(T) - h_in))
    for T in np.linspace(Ts + 0.01, T_ref, 14):
        Teq = max(T - A["approach_eq"], Ts)
        Ye, _ = equil(Y, Teq, P)
        pts.append((T, hs(St(1, T, P, Ye))[0] - h_in))
    return pts


def pinch_ok(hot, cold_pts, m_cold):
    """противоток: горячий поток hot (вход), холодный — точки (T, q на кг) × m_cold; минимальная разность температур"""
    h_in, _ = hs(hot)
    Qtot = cold_pts[-1][1] * m_cold
    dmin = 1e9
    for T, q in cold_pts:
        Qh = Qtot - q * m_cold                       # тепло, отданное горячим от его входа до этого сечения
        Th = T_of_h(h_in - Qh / hot.m, hot.P, hot.Y)
        dmin = min(dmin, Th - T)
    return dmin


def ua_curves(hot, cold_pts, m_cold):
    """UA = ∫ dQ/ΔT по составным кривым (противоток), Вт/К"""
    h_in, _ = hs(hot)
    Qtot = cold_pts[-1][1] * m_cold
    ua = 0.0
    prev = None
    for T, q in cold_pts:
        Qh = Qtot - q * m_cold
        Th = T_of_h(h_in - Qh / hot.m, hot.P, hot.Y)
        dT = max(Th - T, 1.0)
        if prev is not None:
            ua += abs(q - prev[0]) * m_cold / (0.5 * (dT + prev[1]))
        prev = (q, dT)
    return ua


def make_fuel(tcr, gin, mf_guess, prm, p_fuel, stack):
    """готовит топливо из выхлопа gin; возвращает dict: Y, h (на кг топлива), T, lhv, yC (доля CH4 в сырье по массе),
    q (тепло выхлопа на кг топлива), info, а также потоки сырья для учёта эксергии"""
    if tcr == "none":
        h = hs(St(1, T0, p_fuel, Y_CH4))[0]
        return dict(Y=Y_CH4, h=h, T=T0, lhv=LHV_CH4, yC=1.0, q=0.0, info={}, feeds=[("CH4", 1.0, St(1, T0, p_fuel, Y_CH4))])
    T_top = min(gin.T - A["dT_app"], A["T_ref_max"] + 273.15)
    if tcr == "fh":
        T = max(min(T_top, A["T_fh_max"] + 273.15), T0)
        h0 = hs(St(1, T0, p_fuel, Y_CH4))[0]
        h1 = hs(St(1, T, p_fuel * (1 - A["dp_ref_f"]), Y_CH4))[0]
        return dict(Y=Y_CH4, h=h1, T=T, lhv=LHV_CH4, yC=1.0, q=h1 - h0, info=dict(T_fuel_C=T - 273.15),
                    feeds=[("CH4", 1.0, St(1, T0, p_fuel, Y_CH4))], Pout=p_fuel * (1 - A["dp_ref_f"]))
    if tcr == "steam":
        sc = prm["sc"]
        Yf, h_in, yC = feed_steam(sc, p_fuel)
        T = T_top
        info = {}
        while True:
            pts = cold_curve_steam(sc, p_fuel, T)
            dmin = pinch_ok(gin, pts, mf_guess)
            if dmin >= A["pinch"] or T < T_sat(p_fuel) + 30:
                break
            T -= 10.0
        Teq = max(T - A["approach_eq"], T_sat(p_fuel))
        Ye, _ = equil(Yf, Teq, p_fuel)
        Pout = p_fuel * (1 - A["dp_ref_f"])
        h_out = hs(St(1, T, Pout, Ye))[0]
        nC0 = yC / MW[I["CH4"]]
        X = 1 - (Ye[I["CH4"]] / MW[I["CH4"]]) / nC0
        L = lhv(Ye)
        info = dict(T_ref_C=T - 273.15, X_CH4=X, sc=sc, pinch_K=dmin, lhv_gain=L / yC / LHV_CH4_CT - 1,
                    aC=carbon_activity(Ye, Teq, p_fuel), water_per_CH4_kg=(1 - yC) / yC, cold_pts=pts)
        feeds = [("CH4", yC, St(1, T0, p_fuel, Y_CH4)), ("H2O_liq", 1 - yC, St(1, T0, p_fuel, Y_of("H2O:1"), liq=True))]
        return dict(Y=Ye, h=h_out, T=T, lhv=L, yC=yC, q=h_out - h_in, info=info, feeds=feeds, Pout=Pout)
    if tcr == "egr":
        R = prm["R_egr"]
        # рециркулят: газ за рекуператором; охлаждение воздухом до T_rc_cool и две ступени сжатия с промохлаждением
        # до давления топлива (без охлаждения одна ступень дала бы на выходе 1500 °C — нереально)
        rc = stack.cp(m=1.0)
        Mrc = 1.0 / np.sum(rc.Y / MW)
        mC, mR = MW[I["CH4"]], R * Mrc
        yC = mC / (mC + mR)
        Tc = A["T_rc_cool"] + 273.15
        rho_c = P0 / (RU / 28.96 * T0)
        fan = lambda Q, Th: A["Cr_ic"] * Q / max(Th - T0, 1.0) / 1006.0 * A["dp_fan"] / rho_c / (A["eta_fan"] * A["eta_motor"])
        s1 = rc.cp(T=min(rc.T, Tc))
        q1 = hs(rc)[0] - hs(s1)[0]
        w_fan = fan(q1, rc.T)
        pm = math.sqrt(rc.P * p_fuel)
        w1 = poly(s1, pm, A["eta_rc"], True)
        s2 = s1.cp(T=min(s1.T, Tc))
        q2 = hs(s1)[0] - hs(s2)[0]
        w_fan += fan(q2, s1.T)
        w2 = poly(s2, p_fuel, A["eta_rc"], True)
        rc_c = s2
        w_rc = w1 + w2 + w_fan
        q_rc = q1 + q2 + w_fan
        feed = mix([St(yC, T0, p_fuel, Y_CH4), rc_c.cp(m=1 - yC)], p_fuel)
        h_in = hs(feed)[0]
        T = max(T_top, feed.T)
        Pout = p_fuel * (1 - A["dp_ref_f"])
        Teq = max(T - A["approach_eq"], feed.T)
        Ye, _ = equil(feed.Y, Teq, p_fuel)
        h_out = hs(St(1, T, Pout, Ye))[0]
        if h_out < h_in:                      # экзотермика (горит кислород рециркулята): реактор адиабатный
            Ye, T = equil(feed.Y, h_in, p_fuel, "HP")
            Teq = T
            h_out = h_in
        L = lhv(Ye)
        nC0 = yC / MW[I["CH4"]]
        X = 1 - (Ye[I["CH4"]] / MW[I["CH4"]]) / nC0
        g = gas()
        g.TPY = feed.T, feed.P, feed.Y
        Xf = g.X
        info = dict(T_ref_C=T - 273.15, X_CH4=X, R_egr=R, O2_per_CH4=Xf[I["O2"]] / Xf[I["CH4"]],
                    H2O_per_CH4=Xf[I["H2O"]] / Xf[I["CH4"]], CO2_per_CH4=Xf[I["CO2"]] / Xf[I["CH4"]],
                    lhv_gain=L / yC / LHV_CH4_CT - 1, aC=carbon_activity(Ye, Teq, p_fuel), w_rc_per_kg_rc=w_rc,
                    T_rc_C=rc_c.T - 273.15, q_rc_cool_per_kg_rc=q1 + q2)
        feeds = [("CH4", yC, St(1, T0, p_fuel, Y_CH4)), ("RC", 1 - yC, rc_c.cp(m=1.0))]
        return dict(Y=Ye, h=h_out, T=T, lhv=L, yC=yC, q=h_out - h_in, info=info, feeds=feeds, Pout=Pout,
                    w_rc=w_rc * (1 - yC), q_rc=q_rc * (1 - yC), rc_in=rc)
    raise ValueError(tcr)


# ======================================================================== один вариант цикла
DEFAULT_PRM = dict(x_rh=0.5, sc=2.0, R_egr=1.0)


def evaluate(topo, PR, prm=None, detail=False):
    """topo: rec (0 или степень регенерации), ic, rh, ibc (bool), pgc ('none','real','ideal'),
    tcr ('none','fh','steam','egr'). Возвращает словарь показателей на 1 кг/с воздуха и при мощности 5922 кВт."""
    prm = {**DEFAULT_PRM, **(prm or {})}
    rec, ic, rh, ibc, pgc, tcr = (topo[k] for k in ("rec", "ic", "rh", "ibc", "pgc", "tcr"))
    T41 = prm.get("T41_C", A["T41_C"]) + 273.15
    T_rec_max = prm.get("T_rec_max", A["T_rec_max"])
    eta_t = A["eta_t"] - (A["pgc_deta_t"] if pgc == "real" else 0.0)
    comps = []                          # (имя, входы, выходы, работа наружу Вт/(кг/с), тепловая потеря)

    # --- компрессор
    a = St(1.0, T0, P0, Y_AIR)
    if ic:
        nc = max(2, n_stages(PR))
        PR1 = math.sqrt(PR)
        a1 = a.cp()
        w1 = poly(a1, P0 * PR1, eta_stage(PR1), True)
        comps.append(("Компрессор, ступень 1", [a], [a1.cp()], -w1, 0.0))
        outs, cin, cout, Qic, Wf_ic = cool_with_air(a1, A["eps_ic"], A["dp_ic"])
        a1c = outs[0]
        comps.append(("Промохладитель воздух-воздух", [a1, cin], [a1c, cout], 0.0, 0.0))
        comps.append(("Вентилятор промохладителя", [], [], -Wf_ic, Wf_ic))
        a2 = a1c.cp()
        w2 = poly(a2, a1c.P * PR / PR1, eta_stage(PR / PR1), True)
        comps.append(("Компрессор, ступень 2", [a1c], [a2.cp()], -w2, 0.0))
        w_c = w1 + w2
        lost_cool = [cout]
    else:
        nc = n_stages(PR)
        a2 = a.cp()
        w_c = poly(a2, P0 * PR, eta_stage(PR ** (1.0 / nc)), True)
        comps.append(("Компрессор", [a], [a2.cp()], -w_c, 0.0))
        Wf_ic, Qic, lost_cool = 0.0, 0.0, []
    leak = a2.cp(m=A["x_leak"])
    air2 = a2.cp(m=1 - A["x_leak"])
    comps.append(("Утечки (уходят из цикла)", [a2.cp()], [air2, leak], 0.0, 0.0))
    p3 = air2.P * (1 - A["dp_ra"]) if rec else air2.P
    p_fuel = p3 * A["p_fuel_k"]
    dpg_t = {"none": 0.0, "fh": A["dp_fh_g"]}.get(tcr, A["dp_ref_g"])
    chain = (1 + dpg_t) * (1 + A["dp_rg"] if rec else 1.0)
    p5 = P0 * A["p_exh"] * chain

    fuel = make_fuel("none", None, None, prm, p_fuel, None)
    T3r = air2.T
    stack = None
    mf_prev, T5_prev = None, None
    conv = False
    for it in range(40):
        loop = []
        air3 = St(air2.m, T3r, p3, Y_AIR)
        if pgc != "none":
            mf_try = mf_prev if mf_prev else 0.015
            PG, pgi = pgc_gain(air3, fuel["Y"], fuel["h"], mf_try, fuel["lhv"])
            f = 1.0 if pgc == "ideal" else A["pgc_f_real"]
            PRb = (1 - A["dp_comb"]) + f * (PG - (1 - A["dp_comb"]))
        else:
            PG, pgi, PRb = None, None, 1 - A["dp_comb"]
        p4 = p3 * PRb
        Pf = fuel.get("Pout", p_fuel)
        mf1, st4 = burn(air3, fuel["Y"], fuel["h"], T41, p4, fuel["lhv"])
        loop.append(("Камера сгорания", [air3, St(mf1, fuel["T"], Pf, fuel["Y"])], [st4.cp()], 0.0,
                     (1 - A["eta_b"]) * mf1 * fuel["lhv"]))
        if rh:
            p_rh = p4 * (p5 / p4) ** prm["x_rh"]
            s1 = st4.cp()
            wt1 = poly(s1, p_rh, eta_t, False)
            loop.append(("Турбина до подогрева", [st4.cp()], [s1.cp()], s1.m * wt1, 0.0))
            mf2, s2 = burn(s1, fuel["Y"], fuel["h"], T41, p_rh * (1 - A["dp_rh"]), fuel["lhv"])
            loop.append(("Камера подогрева", [s1.cp(), St(mf2, fuel["T"], Pf, fuel["Y"])], [s2.cp()], 0.0,
                         (1 - A["eta_b"]) * mf2 * fuel["lhv"]))
            s5 = s2.cp()
            wt2 = poly(s5, p5, eta_t, False)
            loop.append(("Турбина после подогрева", [s2.cp()], [s5.cp()], s5.m * wt2, 0.0))
            W_t = s1.m * wt1 + s5.m * wt2
            nt = n_turb(wt1) + n_turb(wt2)
        else:
            mf2 = 0.0
            s5 = st4.cp()
            wt = poly(s5, p5, eta_t, False)
            loop.append(("Турбина", [st4.cp()], [s5.cp()], s5.m * wt, 0.0))
            W_t = s5.m * wt
            nt = n_turb(wt)
        mf = mf1 + mf2
        st5 = s5
        # --- выхлопной тракт: риформер / подогреватель топлива -> рекуператор
        g_in = st5
        if tcr != "none":
            fuel_new = make_fuel(tcr, st5, mf, prm, p_fuel, stack if stack is not None else st5.cp(T=500.0, P=P0 * A["p_exh"]))
            hg = hs(st5)[0]
            Pg = st5.P / (1 + dpg_t)
            g_ref = St(st5.m, T_of_h(hg - mf * fuel_new["q"] / st5.m, Pg, st5.Y), Pg, st5.Y)
            g_in = g_ref
        else:
            fuel_new = fuel
        if rec:
            if g_in.T < air2.T + 10.0:
                return dict(ok=False, why="рекуператор бесполезен: газ на входе холоднее воздуха за компрессором")
            air_r, g_rec, Qrec = recuperate(air2, g_in, rec, A["dp_ra"], A["dp_rg"])
            T3r_new = air_r.T
            g_out = g_rec
        else:
            T3r_new, g_out, Qrec = air2.T, g_in, 0.0
        stack_new = g_out
        dT = abs(T3r_new - T3r) + (abs(st5.T - T5_prev) if T5_prev else 1.0)
        dm = abs(mf - mf_prev) / mf if mf_prev else 1.0
        T3r, stack, mf_prev, T5_prev = T3r_new, stack_new, mf, st5.T
        fuel_used = fuel
        fuel = fuel_new
        if dT < 0.02 and dm < 1e-6 and it > 1 and (tcr == "none" or abs(fuel_new["T"] - fuel_used["T"]) < 0.05):
            conv = True
            break
    if not conv:
        return dict(ok=False, why="итерации не сошлись")
    fuel = fuel_used
    # --- учёт узлов выхлопного тракта (по сошедшемуся состоянию)
    comps += loop
    if tcr != "none":
        feeds = [f[2].cp(m=f[1] * mf) for f in fuel["feeds"]]
        fuel_out = St(mf, fuel["T"], fuel.get("Pout", p_fuel), fuel["Y"])
        if tcr == "egr":
            rc_in = fuel["rc_in"].cp(m=mf * (1 - fuel["yC"]))
            comps.append(("Компрессор рециркуляции (2 ступени, охлаждение воздухом)", [rc_in], [feeds[1]], -mf * fuel["w_rc"],
                          mf * fuel["q_rc"]))
        if tcr == "steam":
            w_pump = mf * (1 - fuel["yC"]) * (p_fuel - P0) / 999.1 / A["eta_pump"]
            comps.append(("Насос питательной воды", [St(mf * (1 - fuel["yC"]), T0, P0, Y_of("H2O:1"), liq=True)], [feeds[1]],
                          -w_pump, w_pump * (1 - A["eta_pump"])))
        else:
            w_pump = 0.0
        name = {"fh": "Подогреватель топлива", "steam": "Риформер паровой (с парогенератором)",
                "egr": "Риформер с рециркуляцией выхлопа"}[tcr]
        comps.append((name, [st5.cp()] + feeds, [g_ref.cp(), fuel_out], 0.0, 0.0))
        comps.append(("Разделение топлива по камерам", [fuel_out], [St(mf1, fuel["T"], fuel_out.P, fuel["Y"]),
                                                                  St(mf2, fuel["T"], fuel_out.P, fuel["Y"])] if rh else
                      [St(mf1, fuel["T"], fuel_out.P, fuel["Y"])], 0.0, 0.0))
    else:
        w_pump = 0.0
    if rec:
        air_r, g_rec, Qrec = recuperate(air2, g_in, rec, A["dp_ra"], A["dp_rg"])
        comps.append(("Рекуператор", [air2.cp(), g_in.cp()], [air_r, g_rec], 0.0, 0.0))
        comps.append(("Канал к камере", [air_r.cp()], [St(air2.m, T3r, p3, Y_AIR)], 0.0, 0.0))
    else:
        comps.append(("Канал к камере", [air2.cp()], [St(air2.m, T3r, p3, Y_AIR)], 0.0, 0.0))
    T_rec_in = g_in.T if rec else None
    # --- обращённый цикл Брайтона на выхлопе
    W_ibc_t = W_ibc_c = Wf_ibc = 0.0
    ibc_info = None
    exh = stack
    lost = [leak] + lost_cool
    # --- рециркуляция: часть газа из-за рекуператора (до IBC) уходит назад в риформер
    if tcr == "egr":
        m_rc = mf * (1 - fuel["yC"])
        out = exh.cp(m=exh.m - m_rc)
        comps.append(("Отбор рециркулята", [exh.cp()], [out, fuel["rc_in"].cp(m=m_rc)], 0.0, 0.0))
        exh = out
    if ibc:
        best = None
        for pl in A["plow"]:
            r = ibc_run(exh, pl)
            if r and (best is None or r["w_net"] > best["w_net"]):
                best = r
        if best is not None and best["w_net"] > 0:
            ibc_info = best
            W_ibc_t, W_ibc_c, Wf_ibc = best["W_t"], best["W_c"], best["W_fan"]
            comps += best["comps"]
            exh = best["out"]
            lost += best["lost"]
        else:
            ibc_info = dict(useless=True)
    W_rc = mf * fuel.get("w_rc", 0.0)
    W_mech = (1 - A["eta_m"]) * (W_t + W_ibc_t)
    comps.append(("Механические потери", [], [], -W_mech, W_mech))
    if ibc_info and not ibc_info.get("useless"):
        pass
    w_net = A["eta_m"] * (W_t + W_ibc_t) - w_c - W_ibc_c - W_rc - Wf_ic - Wf_ibc - w_pump
    m_ch4 = mf * fuel["yC"]
    if w_net <= 0:
        return dict(ok=False, why="нет полезной работы")
    eff = w_net / (m_ch4 * LHV_CH4)
    W1 = A["P_shaft_kW"] * 1e3 / w_net
    W_exh = W1 * exh.m
    EGT = exh.T - 273.15
    res = dict(ok=True, topo=dict(topo), PR=PR, prm={k: prm[k] for k in prm}, eff=eff,
               eff_el=eff * A["eta_gen"] * A["eta_gb"], w_kJkg=w_net / 1e3, W_air=W1, W_exh=W_exh,
               W_fuel_CH4=W1 * m_ch4, EGT_C=EGT, T5_C=st5.T - 273.15, T3_C=air2.T - 273.15, T3r_C=T3r - 273.15,
               T_rec_in_C=(T_rec_in - 273.15) if T_rec_in else None, p3_bar=p3 / 1e5, p4_bar=p4 / 1e5,
               nc=nc, nt=nt, PG_ideal=PG, PRb=PRb, pgc_info=pgi,
               Q_rec_kW=W1 * Qrec / 1e3, Q_ic_kW=W1 * Qic / 1e3, W_fans_kW=W1 * (Wf_ic + Wf_ibc) / 1e3,
               W_rc_kW=W1 * W_rc / 1e3, fuel_T_C=fuel["T"] - 273.15, fuel_info=fuel["info"],
               water_t_h=(W1 * mf * (1 - fuel["yC"]) * 3.6 if tcr == "steam" else 0.0), ibc=ibc_info, mf=mf, m_ch4=m_ch4)
    R_exh = RU / (1.0 / np.sum(exh.Y / MW))
    res["Vexh_ratio"] = (W_exh * R_exh * exh.T / (P0 * A["p_exh"])) / (IN["W_exh0"] * RU / 28.6 * (IN["EGT0"] + 273.15) / (P0 * A["p_exh"]))
    res["limits"] = dict(G_ok=W1 <= A["G_cap"] + 1e-9, Wexh_ok=W_exh <= A["W_exh_cap"] + 1e-9, EGT_ok=EGT <= A["EGT_cap"] + 1e-9,
                         rec_ok=(T_rec_in is None) or (T_rec_in - 273.15 <= T_rec_max + 1e-6))
    res["limits"]["all"] = all(res["limits"].values())
    if detail:
        res["_comps"] = comps
        res["_lost"] = lost + [exh]
        res["_streams"] = dict(air2=air2, st5=st5, g_in=g_in, stack=stack, exh=exh, fuel=fuel, mf1=mf1, mf2=mf2,
                               p_fuel=p_fuel, W_t=W_t, w_c=w_c, W_mech=W_mech, Wf_ic=Wf_ic, Wf_ibc=Wf_ibc, w_pump=w_pump,
                               W_rc=W_rc, W_ibc_t=W_ibc_t, W_ibc_c=W_ibc_c, Qic=Qic, Qrec=Qrec)
    return res


def ibc_run(exh, plow):
    """обращённый цикл Брайтона: расширение газа из трубы до plow·P0, охлаждение воздухом (с конденсацией воды),
    сжатие до противодавления трубы. Работа на 1 кг/с воздуха двигателя."""
    p_in = P0 * A["p_exh"]
    s = exh.cp(P=p_in)
    t_out = s.cp()
    wt = poly(t_out, plow * P0, A["eta_ibc_t"], False)
    W_t = s.m * wt
    outs, cin, cout, Q, Wfan = cool_with_air(t_out, A["eps_ibc"], A["dp_ibc"], condense=True)
    vap = outs[0]
    c_out = vap.cp()
    wc = poly(c_out, p_in, A["eta_ibc_c"], True)
    W_c = vap.m * wc
    w_net = A["eta_m"] * W_t - W_c - Wfan
    comps = [("Турбина IBC", [s.cp()], [t_out.cp()], W_t, 0.0),
             ("Охладитель IBC (воздух, с конденсацией)", [t_out.cp(), cin], outs + [cout], 0.0, 0.0),
             ("Вентилятор охладителя IBC", [], [], -Wfan, Wfan),
             ("Компрессор IBC", [vap.cp()], [c_out.cp()], -W_c, 0.0)]
    lost = [cout] + ([outs[1]] if len(outs) > 1 else [])
    # коэффициент расхода колеса компрессора IBC на общем валу 14 950 об/мин
    g = gas()
    g.TPY = vap.T, vap.P, vap.Y
    rho1 = g.density
    U2 = math.sqrt(wc / 0.70)
    D2 = 60 * U2 / (math.pi * IN["n_out"])
    m_abs = vap.m * 25.0
    phi = m_abs / (rho1 * U2 * D2 ** 2)
    D2_opt = math.sqrt(m_abs / (rho1 * U2 * 0.085))
    return dict(plow=plow, w_net=w_net, W_t=W_t, W_c=W_c, W_fan=Wfan, Q=Q, comps=comps, lost=lost, out=c_out.cp(),
                T_cool_C=vap.T - 273.15, cond_frac=(outs[1].m / s.m if len(outs) > 1 else 0.0),
                phi_shaft=phi, D2_shaft=D2, D2_opt=D2_opt, rpm_opt=60 * U2 / (math.pi * D2_opt), U2=U2,
                hot_in=t_out.cp(), hot_out=vap.cp(), cold=(cin, cout))


# ======================================================================== эксергия по узлам
def exergy_account(res):
    """разрушение эксергии по узлам: (A) баланс эксергии, (B) T0·S_gen; кВт при мощности 5922 кВт"""
    W1 = res["W_air"]
    rows = []
    for name, ins, outs, W, Qloss in res["_comps"]:
        Ein = sum(s.m * ex(s) for s in ins)
        Eout = sum(s.m * ex(s) for s in outs)
        Xa = Ein - Eout - W
        Sgen = sum(s.m * hs(s)[1] for s in outs) - sum(s.m * hs(s)[1] for s in ins) + Qloss / T0
        Hin = sum(s.m * hs(s)[0] for s in ins)
        Hout = sum(s.m * hs(s)[0] for s in outs)
        dH = Hin - Hout - W - Qloss
        Xb = T0 * Sgen
        rows.append(dict(node=name, X_A_kW=Xa * W1 / 1e3, X_B_kW=Xb * W1 / 1e3, energy_residual_kW=dH * W1 / 1e3))
    st = res["_streams"]
    E_fuel = sum(f[1] * res["mf"] * ex(f[2] if f[0] == "CH4" else f[2].cp(P=P0)) for f in st["fuel"]["feeds"] if f[0] != "RC")
    E_air = ex(St(1.0, T0, P0, Y_AIR))
    E_cool_in = 0.0
    for name, ins, outs, W, Q in res["_comps"]:
        for s in ins:
            if (not s.liq) and abs(s.T - T0) < 1e-6 and abs(s.P - P0) < 1e-3 and name.startswith(("Промохладитель", "Охладитель")):
                E_cool_in += s.m * ex(s)
    lost = [(s.m * ex(s)) for s in res["_lost"]]
    X_dest = sum(r["X_A_kW"] for r in rows) * 1e3 / W1
    P = res["w_kJkg"] * 1e3
    bal = (E_fuel + E_air + E_cool_in) - P - X_dest - sum(lost)
    comb = sum(r["X_A_kW"] for r in rows if r["node"].startswith("Камера"))
    return dict(rows=rows, E_fuel_kW=E_fuel * W1 / 1e3, E_lost_kW=sum(lost) * W1 / 1e3, X_dest_kW=X_dest * W1 / 1e3,
                X_comb_kW=comb, X_comb_share=comb * 1e3 / W1 / E_fuel, eta_ex=P / E_fuel,
                balance_residual_pct=100 * bal / E_fuel,
                AB_max_diff_kW=max(abs(r["X_A_kW"] - r["X_B_kW"]) for r in rows))


# ======================================================================== теплообменники: ламинарные каналы
def plate_laminar_constants(n=400):
    """плоский канал, развитое ламинарное течение: f·Re (Фаннинг, Dh = 2·зазор), Nu при постоянном потоке (H)
    и при постоянной температуре стенки (T) — конечные разности (Пуазейль + задача Гретца на собственное значение)."""
    y = np.linspace(-1, 1, n + 1)            # полузазор = 1
    h = y[1] - y[0]
    u = 1.5 * (1 - y ** 2)                   # средняя скорость 1
    # f·Re: dp/dx = -3 μ um / b² (b — полузазор); τw = 3 μ um / b; f = τw/(ρ um²/2); Re = ρ um 4b/μ -> fRe = 24
    tau_w = -np.gradient(u, h, edge_order=2)[-1]   # μ = 1, um = 1, b = 1
    fRe = tau_w / 0.5 * 4                    # f = τw/(0,5 um²) при ρ = 1; Re = 4b um/ν
    # Nu_H: θ'' = u·C, θ(±1) = 0 (симметрия), C из баланса
    N = n - 1
    D2 = (np.diag(-2 * np.ones(N)) + np.diag(np.ones(N - 1), 1) + np.diag(np.ones(N - 1), -1)) / h ** 2
    ui = u[1:-1]
    th = np.linalg.solve(D2, ui)             # θ'' = u (масштаб)
    thf = np.concatenate([[0], th, [0]])
    q_w = np.gradient(thf, h, edge_order=2)[-1]   # поток на стенке (масштаб)
    Tb = np.trapezoid(u * thf, y) / np.trapezoid(u, y)
    NuH = abs(q_w * 4 / (0 - Tb))            # Nu = q Dh/(k (Tw - Tb)), Dh = 4b
    # Nu_T: θ'' = -λ u θ, наименьшее собственное значение
    from scipy.linalg import eigh
    M = np.diag(ui)
    w, v = eigh(-D2, M)
    lam = w[0]
    thv = np.concatenate([[0], v[:, 0], [0]])
    qv = np.gradient(thv, h, edge_order=2)[-1]
    Tbv = np.trapezoid(u * thv, y) / np.trapezoid(u, y)
    NuT = abs(qv * 4 / Tbv)
    return dict(fRe=float(fRe), NuH=float(NuH), NuT=float(NuT), eigen=float(lam))


_LAM = None


def lam():
    global _LAM
    if _LAM is None:
        _LAM = plate_laminar_constants()
    return _LAM


def tprops(T, P, Y):
    g = gas()
    g.TPY = T, P, Y
    return g.density, g.viscosity, g.thermal_conductivity, g.cp_mass


def size_hx(UA, hot_m, hot_T, hot_P, hot_Y, dPh, cold_m, cold_T, cold_P, cold_Y, dPc, tech):
    """пластинчатая противоточная сердцевина, оба тракта — ламинарные плоские каналы одинакового Dh (наименьший,
    какой позволяет изготовление: меньше Dh — меньше объём). Длина — по более жёсткой потере давления.
    UA, Вт/К; расходы абсолютные, кг/с; dP — допустимые потери, Па (на сердцевину — доля dp_core_share)."""
    t = HX_TECH[tech]
    L0 = lam()
    Dh = t["Dh"]
    rh, mh, kh, _ = tprops(hot_T, hot_P, hot_Y)
    rc, mc, kc, _ = tprops(cold_T, cold_P, cold_Y)
    hh = L0["NuT"] * kh / Dh
    hc = L0["NuT"] * kc / Dh
    U = 1.0 / (1 / hh + 1 / hc + t["t"] / t["k_w"])
    Aw = UA / U
    s = A["dp_core_share"]
    Lh = math.sqrt(s * dPh * rh * Aw * Dh ** 3 / (8 * L0["fRe"] * mh * hot_m))
    Lc = math.sqrt(s * dPc * rc * Aw * Dh ** 3 / (8 * L0["fRe"] * mc * cold_m))
    L = min(Lh, Lc)
    V = Aw * (Dh / 2 + t["t"])
    Afree = Aw * Dh / (4 * L)
    Re_h = hot_m * Dh / (Afree * mh)
    Re_c = cold_m * Dh / (Afree * mc)
    dph = 8 * L0["fRe"] * mh * hot_m * L ** 2 / (rh * Aw * Dh ** 3)
    dpc = 8 * L0["fRe"] * mc * cold_m * L ** 2 / (rc * Aw * Dh ** 3)
    mass = t["rho"] * Aw * t["t"] * (1 + t["k_struct"])
    mod = None
    if t["printer"]:
        bx = t["printer"]["box_m"]
        mod = max(1, math.ceil(V / (bx[0] * bx[1] * bx[2] * t["fill"])))
    return dict(tech=tech, UA_kW_K=UA / 1e3, U=U, A_m2=Aw, V_core_m3=V, L_m=L, A_front_m2=V / L, Dh_mm=Dh * 1e3,
                beta_m2_m3=Aw / V, Re_hot=Re_h, Re_cold=Re_c, dP_hot_Pa=dph, dP_cold_Pa=dpc, mass_kg=mass, modules=mod,
                laminar=max(Re_h, Re_c) < 2300)


def ua_counterflow(Q, Th_in, Th_out, Tc_in, Tc_out):
    d1, d2 = Th_in - Tc_out, Th_out - Tc_in
    if d1 <= 0 or d2 <= 0:
        return None
    lmtd = d1 if abs(d1 - d2) < 1e-6 else (d1 - d2) / math.log(d1 / d2)
    return Q / lmtd


def hx_list(res, tech):
    """все теплообменники варианта; абсолютные расходы при мощности 5922 кВт"""
    if not res.get("_streams"):
        return []
    W1 = res["W_air"]
    st = res["_streams"]
    out = []
    t = res["topo"]
    if t["rec"]:
        air2, gin = st["air2"], st["g_in"]
        air_r, g_rec, Q = recuperate(air2, gin, t["rec"], A["dp_ra"], A["dp_rg"])
        ua = ua_counterflow(Q, gin.T, g_rec.T, air2.T, air_r.T)
        if not ua:
            raise ValueError("рекуператор: пересечение температур")
        s = size_hx(ua * W1, gin.m * W1, 0.5 * (gin.T + g_rec.T), gin.P, gin.Y, A["dp_rg"] * gin.P,
                    air2.m * W1, 0.5 * (air2.T + air_r.T), air2.P, Y_AIR, A["dp_ra"] * air2.P, tech)
        s.update(name="рекуператор", Q_kW=Q * W1 / 1e3)
        out.append(s)
    if t["ic"]:
        for name, ins, outs, W, Ql in res["_comps"]:
            if name.startswith("Промохладитель"):
                hi, ci = ins
                ho, co = outs
                Q = hi.m * (hs(hi)[0] - hs(ho)[0])
                ua = ua_counterflow(Q, hi.T, ho.T, ci.T, co.T)
                if not ua:
                    raise ValueError("промохладитель: пересечение температур")
                s = size_hx(ua * W1, hi.m * W1, 0.5 * (hi.T + ho.T), hi.P, hi.Y, A["dp_ic"] * hi.P,
                            ci.m * W1, 0.5 * (ci.T + co.T), P0, Y_AIR, A["dp_fan"], tech)
                s.update(name="промохладитель", Q_kW=Q * W1 / 1e3)
                out.append(s)
    if t["ibc"] and res["ibc"] and not res["ibc"].get("useless"):
        b = res["ibc"]
        hi, ho = b["hot_in"], b["hot_out"]
        ci, co = b["cold"]
        ua = ua_counterflow(b["Q"], hi.T, ho.T, ci.T, co.T)
        if ua:
            s = size_hx(ua * W1, hi.m * W1, 0.5 * (hi.T + ho.T), hi.P, hi.Y, A["dp_ibc"] * hi.P,
                        ci.m * W1, 0.5 * (ci.T + co.T), P0, Y_AIR, A["dp_fan"], tech)
            s.update(name="охладитель IBC", Q_kW=b["Q"] * W1 / 1e3)
            out.append(s)
    if t["tcr"] != "none":
        f = st["fuel"]
        st5, gref = st["st5"], st["g_in"] if not t["rec"] else None
        mfa = res["mf"]
        if t["tcr"] == "steam":
            ua = ua_curves(st5, f["info"]["cold_pts"], mfa)
        else:
            h_in_cold = sum(x[1] * hs(x[2])[0] for x in f["feeds"])
            T_cold_in = T0 if t["tcr"] == "fh" else mix([x[2].cp(m=x[1]) for x in f["feeds"]], st["p_fuel"]).T
            Q = mfa * f["q"]
            Tg_out = T_of_h(hs(st5)[0] - Q / st5.m, st5.P, st5.Y)
            ua = ua_counterflow(Q, st5.T, Tg_out, T_cold_in, f["T"]) if Q > 0 else None
        if ua:
            Ycold = f["Y"]
            s = size_hx(ua * W1, st5.m * W1, st5.T - 50, st5.P, st5.Y, A["dp_ref_g"] * st5.P,
                        mfa * W1, 0.5 * (f["T"] + T0), st["p_fuel"], Ycold, A["dp_ref_f"] * st["p_fuel"], tech)
            s.update(name={"fh": "подогреватель топлива", "steam": "риформер + парогенератор", "egr": "риформер"}[t["tcr"]],
                     Q_kW=mfa * f["q"] * W1 / 1e3)
            out.append(s)
    return out


def murray_tree(m, rho, mu, A_front, A_free_face, v_root=30.0, patch=0.10, exponent=3.0, K_branch=0.3):
    """дерево-распределитель от газохода до лицевой поверхности сердцевины (H-дерево, покрывающее фронт A_front):
    ветвь уровня k обслуживает площадку A_front/2^k, её длина — половина стороны площадки; деление — до площадки
    стороной patch (дальше поток растекается внутри решётки: у гироида каналы связаны во все стороны — оценка).
    Диаметры — по закону Мюррея D_родитель^a = 2·D_ребёнок^a: a = 3 — ламинар (Murray 1926), 7/3 — турбулентный
    (Uylings 1977) — по выдержке поиска; a = 2 — сохранение площади (постоянная скорость). Всё — оценка."""
    D0 = math.sqrt(4 * m / (math.pi * rho * v_root))
    levels = max(1, math.ceil(math.log2(A_front / patch ** 2)))
    rows, dp, vol = [], 0.0, 0.0
    for k in range(levels + 1):
        n = 2 ** k
        D = D0 * 2 ** (-k / exponent)
        L = 0.5 * math.sqrt(A_front / n)
        u = m / (n * rho * math.pi * D ** 2 / 4)
        Re = rho * u * D / mu
        f = 64 / Re if Re < 2300 else 0.316 * Re ** -0.25
        dpk = f * L / D * rho * u ** 2 / 2 + K_branch * rho * u ** 2 / 2
        dp += dpk
        vol += n * math.pi * D ** 2 / 4 * L
        rows.append(dict(level=k, branches=n, D_mm=D * 1e3, L_m=L, u=u, Re=Re, dP_Pa=dpk))
    u_face = m / (rho * A_free_face)
    return dict(exponent=exponent, D_root_m=D0, levels=levels, end_D_mm=rows[-1]["D_mm"], end_u=rows[-1]["u"],
                u_core_face=u_face, end_area_over_face=(2 ** levels) * math.pi * (rows[-1]["D_mm"] / 1e3) ** 2 / 4 / A_free_face,
                dP_Pa=dp, volume_m3=vol, rows=rows[:2] + rows[-2:])


# ======================================================================== детали, Plug&Play, готовность
TRL = dict(simple=9, rec=9, ic=6, rh=5, ibc=4, pgc=3, fh=8, steam=4, egr=3)
TRL_NOTE = {
    "rec": "Mercury 50 в эксплуатации (рекуператор с первичной поверхностью — fleet/gtu_fleet_mercury50.json)",
    "ic": "промохладитель воздух-воздух на стационарной ГТУ: близкие образцы — WR-21 (промохлаждение водой, отказы "
          "охладителя — по выдержке поиска); воздух-воздух малой ГТУ — оценка 6",
    "rh": "последовательное сгорание доказано на больших ГТУ (GT24/26 — по памяти); межтурбинная камера малой ГТУ — оценка 5",
    "ibc": "обращённый цикл Брайтона — лабораторные и расчётные работы (Bianchi, De Pascale 2011 — по выдержке поиска)",
    "pgc": "детонационные камеры — стенды; нестационарный поток снижает КПД турбины (Xisto и др. 2018 — по выдержке поиска)",
    "fh": "подогрев топливного газа — обычная практика (оценка)",
    "steam": "химически рекуперированная ГТУ — расчётные работы (Kesser, Hoffman, Baughn 1994; Pashchenko 2022 — по выдержке "
             "поиска); серийных малых ГТУ нет — оценка 4",
    "egr": "риформинг газом рециркуляции на ГТУ — нет известных образцов (оценка 3)",
}


def parts(res, hx, tech):
    """добавочные детали и узлы относительно KV-3 простого цикла (2 колеса компрессора + 3 ступени турбины)."""
    t = res["topo"]
    items = []
    mods = lambda nm: sum((h.get("modules") or 1) for h in hx if h["name"].startswith(nm))
    if t["rec"]:
        items += [("рекуператор: модули сердцевины", mods("рекуп")), ("газоходы воздуха к/от рекуператора", 2),
                  ("компенсаторы", 2), ("байпасный клапан (пуск, защита)", 1)]
    if t["ic"]:
        items += [("промохладитель: модули", mods("промо")), ("вентилятор + двигатель", 2), ("переходные каналы", 2),
                  ("компенсаторы", 2)]
    if t["rh"]:
        items += [("межтурбинная камера подогрева", 1), ("топливный коллектор и форсунки", 1),
                  ("регулирующий клапан топлива", 1), ("запальник", 1)]
    if t["ibc"] and res.get("ibc") and not res["ibc"].get("useless"):
        b = res["ibc"]
        items += [("ступень турбины IBC", 1), ("колесо компрессора IBC", 1), ("охладитель IBC: модули", mods("охладитель IBC")),
                  ("вентилятор + двигатель", 2), ("конденсатоотводчик", 1)]
        if b["phi_shaft"] > A["phi_max"]:
            items += [("отдельный вал IBC + редуктор (колесо не помещается на вал 14 950 об/мин)", 2)]
    if t["pgc"] != "none":
        items += [("инициирование детонации / клапаны", 2), ("гаситель пульсаций перед турбиной", 1)]
    if t["tcr"] == "fh":
        items += [("подогреватель топлива", mods("подогр")), ("клапан", 1)]
    if t["tcr"] == "steam":
        items += [("риформер с катализатором + парогенератор: модули", mods("риформер")), ("водоподготовка", 1),
                  ("питательный насос + двигатель", 2), ("клапан пара", 1), ("бак воды", 1)]
    if t["tcr"] == "egr":
        items += [("риформер с катализатором: модули", mods("риформер")), ("компрессор рециркуляции (2 колеса) + двигатель", 3),
                  ("охладители рециркулята", 2), ("вентилятор + двигатель", 2), ("клапан", 1), ("газоход рециркуляции", 1)]
    dwheels = (res["nc"] + res["nt"]) - 5
    added = sum(n for _, n in items) + dwheels
    trl = min([TRL["simple"]] + [TRL[k] for k, on in (("rec", t["rec"]), ("ic", t["ic"]), ("rh", t["rh"]),
                                                    ("ibc", t["ibc"]), ("pgc", t["pgc"] != "none")) if on]
              + ([TRL[t["tcr"]]] if t["tcr"] != "none" else []))
    V = sum(h["V_core_m3"] for h in hx) * 2.0       # с коллекторами и корпусами — ×2 (оценка)
    water = t["tcr"] == "steam"
    if water or not res["limits"]["all"]:
        pnp = "нет"
    elif V == 0 and t["pgc"] == "none" and not t["ibc"]:
        pnp = "да"
    elif V <= 6.0:
        pnp = "частично"
    else:
        pnp = "нет"
    consum = (["катализатор (замена)"] if t["tcr"] in ("steam", "egr") else []) + (["водоподготовка (реагенты)"] if water else [])
    rot_aux = (2 if t["ic"] else 0) + (2 if t["ibc"] else 0) + (2 if t["tcr"] in ("steam", "egr") else 0)
    maint = "простой" if not consum and rot_aux == 0 and V == 0 else ("средний" if not consum and rot_aux <= 2 else "сложный")
    return dict(items=items, dwheels=dwheels, added=added, TRL=trl, V_added_m3=V, mass_hx_kg=sum(h["mass_kg"] for h in hx),
                pnp=pnp, water=water, consumables=consum, maintenance=maint)


# ======================================================================== перебор
PR_GRID = [3.0, 3.5, 4.0, 4.5, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 13.5, 15.0, 17.0, 19.0, 21.0]


QUICK = bool(os.environ.get("CYCLES_QUICK"))
FINE = bool(os.environ.get("CYCLES_FINE"))       # мелкая сетка (задача AWS, aws_jobs_cycles.json)
OUT_DIR = os.environ.get("CYCLES_OUT") or (SCRATCH if QUICK else HERE)
WORKERS = int(os.environ.get("CYCLES_WORKERS", "4"))
EPS_LIST = (0.0, 0.70, 0.75, 0.80, 0.85, 0.88, 0.90, 0.92, 0.94) if FINE else (0.0, 0.70, 0.80, 0.85, 0.90, 0.92)
if FINE:
    PR_GRID = [round(3.0 + 0.5 * i, 2) for i in range(37)]
    A.update(x_rh=(0.3, 0.4, 0.5, 0.6, 0.7), sc=(1.0, 1.5, 2.0, 2.5, 3.0), R_egr=(0.25, 0.5, 1.0, 2.0, 4.0),
             plow=tuple(round(0.2 + 0.05 * i, 2) for i in range(13)))


def topologies():
    if QUICK:      # проверка вывода: 12 топологий, вывод — в scratchpad
        return [dict(rec=r, ic=False, rh=h, ibc=False, pgc="none", tcr=c) for r in (0.0, 0.9) for h in (False, True)
                for c in ("none", "fh", "steam")]
    L = []
    for rec in EPS_LIST:
        for ic in (False, True):
            for rh in (False, True):
                for ibc in (False, True):
                    for pgc in ("none", "real", "ideal"):
                        for tcr in ("none", "fh", "steam", "egr"):
                            L.append(dict(rec=rec, ic=ic, rh=rh, ibc=ibc, pgc=pgc, tcr=tcr))
    return L


def topo_name(t):
    s = []
    if t["rec"]:
        s.append(f"рекуператор {t['rec']:.2f}")
    if t["ic"]:
        s.append("промохлаждение")
    if t["rh"]:
        s.append("подогрев")
    if t["ibc"]:
        s.append("IBC")
    if t["pgc"] == "real":
        s.append("горение с ростом давления (25 %)")
    if t["pgc"] == "ideal":
        s.append("горение с ростом давления (предел)")
    s += {"none": [], "fh": ["подогрев топлива"], "steam": ["паровой риформинг"], "egr": ["риформинг рециркулятом"]}[t["tcr"]]
    return "простой" if not s else " + ".join(s)


def topo_key(t):
    return f"R{int(round(t['rec'] * 100)):02d}-{'IC' if t['ic'] else 'ic'}-{'RH' if t['rh'] else 'rh'}-" \
           f"{'IBC' if t['ibc'] else 'ibc'}-{t['pgc']}-{t['tcr']}"


def prm_grid(t):
    xs = A["x_rh"] if t["rh"] else (0.5,)
    if t["tcr"] == "steam":
        return [dict(x_rh=x, sc=s) for x in xs for s in A["sc"]]
    if t["tcr"] == "egr":
        return [dict(x_rh=x, R_egr=r) for x in xs for r in A["R_egr"]]
    return [dict(x_rh=x) for x in xs]


def _safe(t, PR, prm, detail=False):
    try:
        return evaluate(t, PR, prm, detail)
    except Exception as e:  # noqa: BLE001
        return dict(ok=False, why=f"{type(e).__name__}: {e}")


def score(r):
    """выбор точки: сначала выполнение мягкого допуска, затем наибольший КПД на клеммах"""
    if not r or not r.get("ok"):
        return (-1, -1.0)
    return (1 if r["limits"]["all"] else 0, r["eff_el"])


def optimize_topo(t, extra=None):
    """оптимум π (сетка + уточнение) и параметров; плюс точка «как есть» π 13,5"""
    best = None
    grid_p = [dict(p, **(extra or {})) for p in prm_grid(t)]
    for prm in grid_p:
        grid = [(PR, _safe(t, PR, prm)) for PR in PR_GRID]
        ok = [(PR, r) for PR, r in grid if r.get("ok")]
        if not ok:
            continue
        PRb, rb = max(ok, key=lambda x: score(x[1]))
        i = PR_GRID.index(PRb)
        lo = PR_GRID[max(0, i - 1)]
        hi = PR_GRID[min(len(PR_GRID) - 1, i + 1)]
        # не переходить через границу числа ступеней (одна/две ступени при π 5)
        if PRb <= A["PRs_1stage"]:
            hi = min(hi, A["PRs_1stage"])
        else:
            lo = max(lo, A["PRs_1stage"] + 0.01)
        if hi - lo > 0.05:
            def obj(x):
                sc_ = score(_safe(t, x, prm))
                return -sc_[1] - 10 * sc_[0]
            ref = minimize_scalar(obj,
                                  bounds=(lo, hi), method="bounded", options=dict(xatol=0.05, maxiter=8))
            rr = _safe(t, float(ref.x), prm)
            if score(rr) > score(rb):
                PRb, rb = float(ref.x), rr
        if best is None or score(rb) > score(best[1]):
            best = (PRb, rb, prm)
    fixed = None
    for prm in grid_p:
        r = _safe(t, A["PR0"], prm)
        if r.get("ok") and (fixed is None or score(r) > score(fixed[1])):
            fixed = (A["PR0"], r, prm)
    out = dict(topo=t, key=topo_key(t), name=topo_name(t))
    for tag, b in (("opt", best), ("pr135", fixed)):
        if b is None:
            out[tag] = dict(ok=False)
            continue
        r = _safe(t, b[0], b[2], detail=True)
        out[tag] = finish(r)
    return out


def finish(r):
    """размеры теплообменников, детали, эксергия — для итоговой точки"""
    if not r.get("ok"):
        return dict(ok=False, why=r.get("why"))
    hx = {tech: hx_list(r, tech) for tech in HX_TECH}
    pa = {tech: parts(r, hx[tech], tech) for tech in HX_TECH}
    exa = exergy_account(r)
    clean = {k: v for k, v in r.items() if not k.startswith("_")}
    fi = dict(clean.get("fuel_info") or {})
    fi.pop("cold_pts", None)
    clean["fuel_info"] = fi
    if clean.get("ibc") and not clean["ibc"].get("useless"):
        clean["ibc"] = {k: v for k, v in clean["ibc"].items() if k not in ("comps", "lost", "out", "hot_in", "hot_out", "cold")}
    clean["hx"] = hx
    clean["parts"] = pa
    clean["exergy"] = exa
    clean["T41_for_G_cap_C"] = None
    if not r["limits"]["G_ok"]:
        def fG(T):
            rr = _safe(r["topo"], r["PR"], dict(r["prm"], T41_C=T))
            return (rr["W_air"] - A["G_cap"]) if rr.get("ok") else 1.0
        try:
            if fG(1250.0) < 0:
                clean["T41_for_G_cap_C"] = brentq(fG, A["T41_C"], 1250.0, xtol=0.5)
        except Exception:  # noqa: BLE001
            pass
    return clean


# ======================================================================== исследования (таблицы)
def study_recup_tradeoff(PR=8.22, prm=None, rh=True):
    """компромисс «степень регенерации — потеря давления» (Бежан: минимум производства энтропии = максимум КПД
    при заданном объёме): КПД, объём и масса печатной сердцевины; по умолчанию — рекуператор + подогрев"""
    rows = []
    base = dict(A)
    for eps in (0.70, 0.80, 0.85, 0.90, 0.92):
        for k in (0.5, 1.0, 1.5, 2.0):
            A["dp_ra"], A["dp_rg"] = base["dp_ra"] * k, base["dp_rg"] * k
            r = _safe(dict(rec=eps, ic=False, rh=rh, ibc=False, pgc="none", tcr="none"), PR, prm, detail=True)
            if r.get("ok"):
                hx = hx_list(r, "lpbf")[0]
                rows.append(dict(eps=eps, dp_total_pct=100 * (A["dp_ra"] + A["dp_rg"]), eff_el=r["eff_el"], W_air=r["W_air"],
                                 EGT_C=r["EGT_C"], limits_ok=r["limits"]["all"], V_core_m3=hx["V_core_m3"], mass_kg=hx["mass_kg"],
                                 L_m=hx["L_m"], A_front_m2=hx["A_front_m2"]))
    A["dp_ra"], A["dp_rg"] = base["dp_ra"], base["dp_rg"]
    return rows


def study_reformer():
    """равновесие и кинетика риформинга при давлении топлива KV-3 и при низком π: превращение CH4, прирост LHV,
    активность углерода; кинетика Xu — Froment — масса катализатора до 95 % равновесного превращения"""
    rows = []
    for P_bar in (5.0, 17.0):
        for sc in (1.5, 2.0, 3.0):
            Yf, h_in, yC = feed_steam(sc)
            for T_C in (450, 550, 650, 750, 850):
                T = T_C + 273.15
                Ye, _ = equil(Yf, T, P_bar * 1e5)
                X = 1 - (Ye[I["CH4"]] / MW[I["CH4"]]) / (yC / MW[I["CH4"]])
                q = hs(St(1, T, P_bar * 1e5, Ye))[0] - h_in
                aM, aB = carbon_activity(Ye, T, P_bar * 1e5)
                g = gas()
                g.TPY = T, P_bar * 1e5, Ye
                Xd = g.X
                dry = 1 - Xd[I["H2O"]]
                rows.append(dict(kind="steam", P_bar=P_bar, sc=sc, T_C=T_C, X_CH4=X, lhv_gain=lhv(Ye) / yC / LHV_CH4_CT - 1,
                                 q_MJ_per_kgCH4=q / yC / 1e6, H2_dry=Xd[I["H2"]] / dry, a_C=max(aM, aB)))
    # сухой (CO2) — для сравнения: CH4 + CO2 1:1, сажа
    for P_bar in (5.0, 17.0):
        for T_C in (550, 700, 850):
            Yf = Y_of("CH4:1,CO2:1")
            Ye, _ = equil(Yf, T_C + 273.15, P_bar * 1e5)
            rows.append(dict(kind="dry CH4:CO2 1:1", P_bar=P_bar, T_C=T_C, X_CH4=1 - Ye[I["CH4"]] / Yf[I["CH4"]],
                             a_C=max(carbon_activity(Ye, T_C + 273.15, P_bar * 1e5))))
    # кинетика: паровой S/C 2 при 17 бар и 650 °C; масса катализатора на 1 моль/с CH4
    kin = []
    for P_bar, T_C, sc in ((17.0, 650.0, 2.0), (17.0, 550.0, 2.0), (5.0, 650.0, 2.0), (17.0, 650.0, 3.0)):
        T = T_C + 273.15
        Xf = np.zeros(len(SP))
        Xf[I["CH4"]], Xf[I["H2O"]] = 1 / (1 + sc), sc / (1 + sc)
        Yf = Xf * MW / np.dot(Xf, MW)
        Ye, _ = equil(Yf, T, P_bar * 1e5)
        Xeq = 1 - (Ye[I["CH4"]] / MW[I["CH4"]]) / (Yf[I["CH4"]] / MW[I["CH4"]])
        W, Xc, _ = xf_pfr(Xf, 1.0 + sc, T, P_bar * 1e5, 20.0)
        target = 0.95 * Xeq
        W95 = float(np.interp(target, Xc, W)) if Xc[-1] >= target else None
        kin.append(dict(P_bar=P_bar, T_C=T_C, sc=sc, X_eq=Xeq, X_end=float(Xc[-1]), W95_kg_per_molCH4_s=W95))
    return dict(equilibrium=rows, kinetics=kin)


def _opt_trec800(t):
    return optimize_topo(t, extra=dict(T_rec_max=800.0))


def study_wall(res_best):
    """чувствительность печатной сердцевины рекуператора к стенке и каналу (LPBF): объём, масса, модули"""
    rows = []
    o = res_best
    rr = _safe(o["topo"], o["opt"]["PR"], o["opt"]["prm"], detail=True)
    keep = dict(HX_TECH["lpbf"])
    for Dh, t in ((2.0e-3, 0.35e-3), (2.0e-3, 0.25e-3), (1.5e-3, 0.25e-3), (1.5e-3, 0.20e-3), (3.0e-3, 0.5e-3)):
        HX_TECH["lpbf"].update(Dh=Dh, t=t)
        h = [x for x in hx_list(rr, "lpbf") if x["name"] == "рекуператор"][0]
        rows.append(dict(Dh_mm=Dh * 1e3, t_mm=t * 1e3, V_core_m3=h["V_core_m3"], mass_kg=h["mass_kg"], modules=h["modules"],
                         L_m=h["L_m"], A_front_m2=h["A_front_m2"]))
    HX_TECH["lpbf"].clear()
    HX_TECH["lpbf"].update(keep)
    foil = [x for x in hx_list(rr, "foil") if x["name"] == "рекуператор"][0]
    return dict(key=o["key"], lpbf=rows, foil=dict(V_core_m3=foil["V_core_m3"], mass_kg=foil["mass_kg"], L_m=foil["L_m"]))


def crosscheck_gtu_cycle():
    """сверка: простой цикл и рекуператор против gtu_cycle (тот же Cantera, другой код)"""
    cwd = os.getcwd()
    sys.path.insert(0, GTU)
    import gtu_cycle as C  # noqa: E402  (импорт делает chdir)
    os.chdir(cwd)
    out = {}
    PR, T41 = A["PR0"], A["T41_C"]
    ref = C.cycle(PR, T41, eta_pc=A["eta_pc0"], eta_pt_gg=A["eta_t"], eta_pt_pt=A["eta_t"], x_ngv=0.0, x_rot=0.0,
                  x_leak=A["x_leak"], dp_duct=0.0, W1=IN["G0"])
    mine = evaluate(dict(rec=0.0, ic=False, rh=False, ibc=False, pgc="none", tcr="none"), PR)
    P_mine = mine["w_kJkg"] * IN["G0"]
    out["simple"] = dict(PR=PR, T41_C=T41, W1=IN["G0"], gtu_cycle=dict(P_kW=ref["P_MW"] * 1e3, eff=ref["eff"], EGT_C=ref["EGT_C"]),
                         cycles=dict(P_kW=P_mine, eff=mine["eff"], EGT_C=mine["EGT_C"]),
                         frame3_json=dict(P_kW=IN["P_kW_frame3"], eff=IN["eff_frame3"], EGT_C=IN["EGT_frame3"]))
    # рекуператор π 5, 0,85 — eta_pc одной ступени KV
    eta = eta_stage(5.0)
    refr = C.recup(5.0, T41, 0.85, eta_pc=eta, eta_pt_gg=A["eta_t"], eta_pt_pt=A["eta_t"], x_ngv=0.0, x_rot=0.0,
                   x_leak=A["x_leak"], dp_duct=0.0, W1=1.0)
    m2 = evaluate(dict(rec=0.85, ic=False, rh=False, ibc=False, pgc="none", tcr="none"), 5.0)
    out["recup"] = dict(PR=5.0, eps=0.85, eta_pc=eta, gtu_cycle=dict(eff=refr["eff"], w_kJkg=refr["sp_work_kJkg"],
                                                                  T3r_C=refr["T3r_C"], EGT_C=refr["EGT_C"]),
                        cycles=dict(eff=m2["eff"], w_kJkg=m2["w_kJkg"], T3r_C=m2["T3r_C"], T5_C=m2["T5_C"]),
                        note="gtu_cycle.recup: эффективность по температуре воздуха (T3r = T3 + eps·(T5 − T3)); здесь — по "
                             "энтальпии наименьшей теплоёмкости; разница теплоёмкостей даёт несколько кельвинов")
    return out


# ======================================================================== Парето
def pareto(rows, keys):
    """недоминируемые: keys — список (поле, знак): +1 максимизировать, −1 минимизировать"""
    idx = []
    for i, a in enumerate(rows):
        dom = False
        for j, b in enumerate(rows):
            if i == j:
                continue
            ge = all((b[k] * s) >= (a[k] * s) - 1e-12 for k, s in keys)
            gt = any((b[k] * s) > (a[k] * s) + 1e-12 for k, s in keys)
            if ge and gt:
                dom = True
                break
        if not dom:
            idx.append(i)
    return idx


def flat(o, tag, tech):
    r = o[tag]
    if not r.get("ok"):
        return None
    p = r["parts"][tech]
    fi = r.get("fuel_info") or {}
    return dict(key=o["key"], name=o["name"], point=tag, tech=tech, PR=round(r["PR"], 3), eff_el=r["eff_el"], eff=r["eff"],
                w_kJkg=r["w_kJkg"], W_air=r["W_air"], W_exh=r["W_exh"], EGT_C=r["EGT_C"], T5_C=r["T5_C"],
                T3r_C=r["T3r_C"], limits_ok=r["limits"]["all"], G_ok=r["limits"]["G_ok"], Wexh_ok=r["limits"]["Wexh_ok"],
                EGT_ok=r["limits"]["EGT_ok"], rec_ok=r["limits"]["rec_ok"], nc=r["nc"], nt=r["nt"], dwheels=p["dwheels"],
                added=p["added"], V_added_m3=p["V_added_m3"], mass_hx_kg=p["mass_hx_kg"], TRL=p["TRL"], pnp=p["pnp"],
                water_t_h=r["water_t_h"], maintenance=p["maintenance"], X_comb_kW=r["exergy"]["X_comb_kW"],
                eta_ex=r["exergy"]["eta_ex"], bound=(o["topo"]["pgc"] == "ideal"),
                degenerate=bool(o["topo"]["ibc"] and (r.get("ibc") or {}).get("useless", False)),
                lhv_gain=fi.get("lhv_gain"), X_CH4=fi.get("X_CH4"), T_ref_C=fi.get("T_ref_C"),
                x_rh=r["prm"].get("x_rh") if o["topo"]["rh"] else None,
                sc=r["prm"].get("sc") if o["topo"]["tcr"] == "steam" else None,
                R_egr=r["prm"].get("R_egr") if o["topo"]["tcr"] == "egr" else None,
                ibc_plow=(r["ibc"] or {}).get("plow") if o["topo"]["ibc"] else None)


# ======================================================================== графики
def plots(rows, front, ref):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(9, 5.6), dpi=110)
    ok = [r for r in rows if r["limits_ok"] and not r["bound"]]
    bad = [r for r in rows if not r["limits_ok"] and not r["bound"]]
    ax.scatter([r["added"] for r in bad], [100 * r["eff_el"] for r in bad], s=8, c="#8c8c8c", label="вне мягкого допуска")
    sc_ = ax.scatter([r["added"] for r in ok], [100 * r["eff_el"] for r in ok], s=[10 + 25 * r["V_added_m3"] for r in ok],
                     c=[r["TRL"] for r in ok], cmap="viridis", vmin=3, vmax=9, alpha=0.75, label="в допуске (размер — объём)")
    fr = [rows[i] for i in front]
    ax.scatter([r["added"] for r in fr], [100 * r["eff_el"] for r in fr], s=60, facecolors="none", edgecolors="#d62728",
               linewidths=1.4, label="фронт Парето")
    ax.axhline(100 * ref, color="#555", lw=0.8, ls="--")
    ax.annotate("KV-3 простой, π 13,5 (вне допуска по расходу)", xy=(ax.get_xlim()[1], 100 * ref), xytext=(-4, -11),
                textcoords="offset points", ha="right", fontsize=8, color="#444")
    for r in fr:
        ax.annotate(r["name"].replace("рекуператор", "рек.").replace("подогрев топлива", "подогр. топлива")[:38],
                    xy=(r["added"], 100 * r["eff_el"]), xytext=(6, 3), textcoords="offset points", fontsize=7, color="#222")
    from matplotlib.ticker import MaxNLocator
    ax.xaxis.set_major_locator(MaxNLocator(integer=True))
    plt.colorbar(sc_, ax=ax, label="TRL (оценка)")
    ax.set_xlabel("добавочные детали и узлы (относительно KV-3 простого)")
    ax.set_ylabel("КПД на клеммах, %")
    ax.set_title("KV-3: перебор циклов (печатные теплообменники LPBF), оптимум π")
    ax.legend(fontsize=8, loc="upper left")
    fig.tight_layout()
    fig.savefig(os.path.join(OUT_DIR, "cycles_pareto.png"), bbox_inches="tight")
    plt.close(fig)


EX_GROUPS = [("Компрессор и промохлаждение", ("Компрессор", "Промохладитель", "Вентилятор промо")),
             ("Камеры сгорания", ("Камера",)),
             ("Турбины", ("Турбина",)),
             ("Теплообменники и риформер", ("Рекуператор", "Риформер", "Подогреватель", "Охладитель", "Вентилятор охл", "Насос")),
             ("Прочее (мех. потери, каналы, утечки)", ("Механические", "Канал", "Утечки", "Отбор", "Разделение"))]


def plots_exergy(items):
    """разрушение и потери эксергии по группам узлов, % эксергии топлива; items: [(подпись, результат opt)]"""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    colors = ["#0072B2", "#D55E00", "#009E73", "#CC79A7", "#4D4D4D", "#E69F00", "#56B4E9"]   # Okabe — Ito + тёмно-серый
    labels = [x[0] for x in items]
    fig, ax = plt.subplots(figsize=(10, 0.55 * len(items) + 2.0), dpi=110)
    left = np.zeros(len(items))
    cats = [g[0] for g in EX_GROUPS] + ["Уходящие потоки (выхлоп, воздух охладителей)", "Мощность на валу"]
    for ci, cat in enumerate(cats):
        vals = []
        for _, r in items:
            e = r["exergy"]
            Ef = e["E_fuel_kW"]
            if ci < len(EX_GROUPS):
                v = sum(x["X_A_kW"] for x in e["rows"] if x["node"].startswith(EX_GROUPS[ci][1]))
            elif ci == len(EX_GROUPS):
                v = e["E_lost_kW"]
            else:
                v = r["w_kJkg"] * r["W_air"]
            vals.append(100 * v / Ef)
        vals = np.array(vals)
        ax.barh(range(len(items)), vals, left=left, color=colors[ci], edgecolor="white", linewidth=0.6, label=cat)
        for i, v in enumerate(vals):
            if v > 4:
                ax.text(left[i] + v / 2, i, f"{v:.0f}", ha="center", va="center", fontsize=7,
                        color="white" if ci in (0, 1, 2, 4) else "black")
        left += vals
    ax.set_yticks(range(len(items)))
    ax.set_yticklabels(labels, fontsize=8)
    ax.invert_yaxis()
    ax.set_xlim(0, 101)
    ax.set_xlabel("% эксергии топлива (метан + вода, если есть)")
    ax.set_title("KV-3: куда уходит эксергия топлива (фронт Парето)")
    ax.legend(fontsize=7, ncol=2, loc="upper center", bbox_to_anchor=(0.45, -0.14), frameon=False)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT_DIR, "cycles_exergy.png"), bbox_inches="tight")
    plt.close(fig)


# ======================================================================== главная
def main():
    t0 = time.time()
    L = topologies()
    print("топологий:", len(L), "сверка и исследования...")
    cross = crosscheck_gtu_cycle()
    print("сверка простой:", json.dumps(cross["simple"], ensure_ascii=False))
    print("сверка рекуп.:", json.dumps({k: cross["recup"][k] for k in ("gtu_cycle", "cycles")}, ensure_ascii=False))
    lamc = lam()
    print("ламинарные постоянные:", lamc)
    with ProcessPoolExecutor(max_workers=WORKERS) as exe:
        res = list(exe.map(optimize_topo, L, chunksize=2))
        print(f"перебор: {time.time() - t0:.0f} с")
        L8 = [t for t in L if t["rec"] and t["pgc"] == "none" and not t["ibc"] and t["tcr"] in ("none", "fh")]
        res800 = list(exe.map(_opt_trec800, L8, chunksize=2))
    print(f"T_rec_max 800 °C: {time.time() - t0:.0f} с")
    tech = "lpbf"
    rows = []
    for o in res:
        for tag in ("opt", "pr135"):
            for tc in HX_TECH:
                f = flat(o, tag, tc)
                if f:
                    rows.append(f)
    sel = [r for r in rows if r["point"] == "opt" and r["tech"] == tech and r["limits_ok"] and not r["bound"]
           and not r["degenerate"]]
    front_i = pareto(sel, [("eff_el", 1), ("added", -1), ("V_added_m3", -1)])
    front = sorted([sel[i] for i in front_i], key=lambda r: (r["added"], -r["eff_el"]))
    sel_p = [r for r in sel if r["pnp"] != "нет" and r["TRL"] >= 5]
    front_p = sorted([sel_p[i] for i in pareto(sel_p, [("eff_el", 1), ("added", -1), ("V_added_m3", -1)])],
                     key=lambda r: (r["added"], -r["eff_el"]))
    base = [o for o in res if o["key"] == "R00-ic-rh-ibc-none-none"][0]
    ref_eff = base["pr135"]["eff_el"]
    # исследования — на рекомендуемом семействе «рекуператор + подогрев» (степень регенерации 0,90)
    rec_best = next((r for r in sel if r["key"] == "R90-ic-RH-ibc-none-none"), None) or max(
        [r for r in sel if "-RH-" in r["key"] and not r["key"].startswith("R00") and "-ic-" in r["key"]],
        key=lambda r: r["eff_el"], default=None)
    o_best = [x for x in res if x["key"] == rec_best["key"]][0] if rec_best else None
    trade = study_recup_tradeoff(o_best["opt"]["PR"], o_best["opt"]["prm"]) if o_best else []
    refm = study_reformer()
    murray = None
    if o_best:
        rr = _safe(o_best["topo"], o_best["opt"]["PR"], o_best["opt"]["prm"], detail=True)
        st = rr["_streams"]
        hx = [h for h in hx_list(rr, tech) if h["name"] == "рекуператор"][0]
        W1 = rr["W_air"]
        gin = st["g_in"]
        rho, mu, _, _ = tprops(gin.T, gin.P, gin.Y)
        Afree_face = HX_TECH[tech]["Dh"] / 2 / (HX_TECH[tech]["Dh"] + 2 * HX_TECH[tech]["t"]) * hx["A_front_m2"]
        murray = dict(key=rec_best["key"], A_front_m2=hx["A_front_m2"], A_free_face_m2=Afree_face,
                      trees={str(a): murray_tree(gin.m * W1, rho, mu, hx["A_front_m2"], Afree_face, exponent=a)
                             for a in (3.0, 7.0 / 3.0, 2.0)},
                      dP_budget_gas_Pa=(1 - A["dp_core_share"]) * A["dp_rg"] * gin.P,
                      V_core_m3=hx["V_core_m3"])
    rows800 = [f for o in res800 for f in [flat(o, "opt", tech)] if f]
    sel800 = [r for r in rows800 if r["limits_ok"]]
    front800 = sorted([sel800[i] for i in pareto(sel800, [("eff_el", 1), ("added", -1), ("V_added_m3", -1)])],
                      key=lambda r: (r["added"], -r["eff_el"])) if sel800 else []
    wall = study_wall(o_best) if o_best else None
    # Plug&Play реперы: пакеты Taurus 60 GS и Mercury 50 (тот же класс мощности, с рекуператором)
    m50 = IN["m50"]
    v_t60 = IN["pkg"][0] * IN["pkg"][1] * IN["pkg"][2]
    v_m50 = m50["pkg"][0] * m50["pkg"][1] * m50["pkg"][2]
    pack_ref = dict(T60_pkg_m3=v_t60, M50_pkg_m3=v_m50, dV_m3=v_m50 - v_t60, T60_pkg_kg=IN["t60_pkg_kg"],
                    M50_pkg_kg=m50["pkg_kg"], dM_kg=m50["pkg_kg"] - IN["t60_pkg_kg"],
                    note="разница пакетов Mercury 50 (4,6 МВт, рекуператор) и Taurus 60 GS — реальный масштаб «места под "
                         "рекуператор» в пакете этого класса (паспорта fleet/*.json)")
    out = dict(
        title="KV-3: плотный перебор термодинамических циклов (organic/cycles)",
        method="0D-цикл на одном валу: свойства газа Cantera (GRI-Mech 3.0, 8 веществ), полное сгорание, политропные процессы "
               "40 шагами (как gtu_cycle); узлы КПД — KV-3 frame3; ламинарные плоские каналы теплообменников (f·Re и Nu посчитаны "
               "конечными разностями); равновесие риформинга — Cantera, кинетика — Xu & Froment (1989) в реакторе идеального "
               "вытеснения с константами равновесия Cantera; горение с ростом давления — эквивалентное располагаемое давление "
               "(простейший вид); эксергия — баланс узла и T0·S_gen (среда ISO 3977-2, RH 60 %). Всё — оценка уровня облика.",
        software=f"python {platform.python_version()}, cantera {ct.__version__}, numpy {np.__version__}, scipy {scipy.__version__}",
        date=DATE, status="черновик: расчёт облика, все КПД узлов, потери и счёт деталей — оценка",
        inputs=IN_SOURCES, assumptions={k: (list(v) if isinstance(v, tuple) else v) for k, v in A.items()},
        hx_tech=HX_TECH, laminar_constants=lamc, crosscheck=cross,
        lhv_ch4=dict(model=LHV_CH4, cantera_T0=LHV_CH4_CT),
        exergy_ch4=dict(e_ch_kJ_mol=E_SP[I["CH4"]] / 1e6, e_ch_MJ_kg=E_SP[I["CH4"]] / MW[I["CH4"]] / 1e6,
                        note="табличное 831,65 кДж/моль (Szargut) — по памяти, для сравнения"),
        trl=TRL, trl_note=TRL_NOTE, pack_reference=pack_ref,
        n_topologies=len(L), base=dict(key=base["key"], pr135=base["pr135"], opt=base["opt"]),
        results=res, pareto_lpbf=front, pareto_lpbf_pnp=front_p, recup_tradeoff=trade, reformer=refm, murray=murray,
        trec800=dict(note="газ на входе в печатный рекуператор из Inconel 625 до 800 °C вместо LIM 700 °C — оценка, нужна "
                          "проверка ползучести (навык calculix-rotor: *CREEP)", rows=rows800, pareto=front800),
        hx_wall_sensitivity=wall,
        runtime_s=time.time() - t0,
    )
    out["feed_to_funnels"] = feed_from(out)
    keep = {base["key"]} | {r["key"] for r in front + front_p} | ({rec_best["key"]} if rec_best else set())
    if OUT_DIR == HERE:          # полный файл (12 МБ) — в scratchpad, в репозиторий — сжатый
        os.makedirs(SCRATCH, exist_ok=True)
        json.dump(out, open(os.path.join(SCRATCH, "cycles_results_full.json"), "w", encoding="utf-8"), ensure_ascii=False,
                  default=_jdef)
    json.dump(compact(out, keep), open(os.path.join(OUT_DIR, "cycles_results.json"), "w", encoding="utf-8"),
              ensure_ascii=False, indent=1, default=_jdef)
    cols = list(rows[0].keys())
    with open(os.path.join(OUT_DIR, "cycles.csv"), "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=cols)
        w.writeheader()
        for r in rows:
            w.writerow({k: (round(v, 5) if isinstance(v, float) else v) for k, v in r.items()})
    plots(sel + [r for r in rows if r["point"] == "opt" and r["tech"] == tech and (not r["limits_ok"]) and not r["bound"]],
          front_i, ref_eff)
    items = [("KV-3 простой, π опт.", base["opt"])]
    for r in front[:8]:
        o = [x for x in res if x["key"] == r["key"]][0]
        if o["key"] != base["key"]:
            items.append((f"{o['name']} (π {r['PR']:.1f})"[:70], o["opt"]))
    plots_exergy(items)
    print("feed:", write_feed(out["feed_to_funnels"]))
    print(f"готово за {time.time() - t0:.0f} с; фронт Парето: {len(front)}")
    for r in front:
        print(f"  {r['added']:3d} дет. {100 * r['eff_el']:5.1f} % V {r['V_added_m3']:5.2f} м³ π {r['PR']:5.2f} TRL {r['TRL']} "
              f"P&P {r['pnp']:8s} | {r['name']}")


def _r5(x):
    """округление до 5 значащих цифр для сжатого файла"""
    if isinstance(x, (bool, np.bool_)):
        return bool(x)
    if isinstance(x, (float, np.floating)):
        return float(f"{float(x):.5g}")
    if isinstance(x, np.integer):
        return int(x)
    if isinstance(x, dict):
        return {k: _r5(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [_r5(v) for v in x]
    return x


def compact(out, keep_keys):
    """сжатый результат для репозитория: у вариантов вне фронтов и опорных — без построчной эксергии и перечня деталей
    исполнения foil; числа — 5 значащих цифр. Полный — scratchpad/organic/cycles_results_full.json"""
    c = dict(out)
    res = []
    for o in out["results"]:
        oo = dict(o)
        for tag in ("opt", "pr135"):
            r = dict(o[tag])
            if r.get("ok") and o["key"] not in keep_keys:
                ex_ = dict(r["exergy"])
                ex_.pop("rows", None)
                r["exergy"] = ex_
                r["parts"] = {t: {k: v for k, v in p.items() if k != "items" or t == "lpbf"} for t, p in r["parts"].items()}
            oo[tag] = r
        res.append(oo)
    c["results"] = res
    c["compact_note"] = ("сжатый файл: построчная эксергия — только у опорного варианта, фронтов Парето и рекомендуемого; "
                         "полный — scratchpad/organic/cycles_results_full.json (не в репозитории)")
    return _r5(c)


def _jdef(x):
    if isinstance(x, np.bool_):
        return bool(x)
    if isinstance(x, np.integer):
        return int(x)
    if isinstance(x, np.floating):
        return float(x)
    if isinstance(x, np.ndarray):
        return x.tolist()
    return str(x)


def _pick(res, pred, tag="opt", need_limits=True):
    c = [o for o in res if pred(o["topo"]) and o[tag].get("ok") and (o[tag]["limits"]["all"] or not need_limits)]
    return max(c, key=lambda o: o[tag]["eff_el"]) if c else None


def combustor_conditions(o, tag="opt"):
    """условия на входе и выходе каждой камеры варианта (для воронки камеры)"""
    r = _safe(o["topo"], o[tag]["PR"], o[tag]["prm"], detail=True)
    rows = []
    g = gas()
    for name, ins, outs, W, Q in r["_comps"]:
        if name.startswith("Камера"):
            ox, fu = ins[0], ins[1]
            g.TPY = ox.T, ox.P, ox.Y
            Xi = g.X.copy()
            g.TPY = outs[0].T, outs[0].P, outs[0].Y
            Xo = g.X.copy()
            g.TPY = fu.T, fu.P, fu.Y
            Xf = g.X.copy()
            rows.append(dict(node=name, T_in_C=ox.T - 273.15, p_in_bar=ox.P / 1e5, O2_in_vol=float(Xi[I["O2"]]),
                             O2_out_vol=float(Xo[I["O2"]]), T_out_C=outs[0].T - 273.15, fuel_air_mass=fu.m / ox.m,
                             fuel_T_C=fu.T - 273.15, fuel_H2_vol=float(Xf[I["H2"]]), fuel_CH4_vol=float(Xf[I["CH4"]])))
    return rows


def feed_from(out):
    """что подать в воронку двигателя (routes/gtu/funnel/funnel_space.py) и воронку камеры (routes/gtu/combustor):
    новые оси и значения, модели, ограничения — с числами из этого расчёта"""
    res = out["results"]
    simple = _pick(res, lambda t: not t["rec"] and not t["ic"] and not t["rh"] and not t["ibc"] and t["pgc"] == "none"
                   and t["tcr"] == "none")
    rec_rh = _pick(res, lambda t: t["rec"] and t["rh"] and not t["ic"] and not t["ibc"] and t["pgc"] == "none" and t["tcr"] == "none")
    rec_rh_fh = _pick(res, lambda t: t["rec"] and t["rh"] and not t["ic"] and not t["ibc"] and t["pgc"] == "none" and t["tcr"] == "fh")
    rec_only = _pick(res, lambda t: t["rec"] and not t["rh"] and not t["ic"] and not t["ibc"] and t["pgc"] == "none"
                     and t["tcr"] == "none", need_limits=False)
    icr = _pick(res, lambda t: t["rec"] and t["ic"] and not t["rh"] and not t["ibc"] and t["pgc"] == "none" and t["tcr"] == "none")
    fh = _pick(res, lambda t: not t["rec"] and not t["ic"] and not t["rh"] and not t["ibc"] and t["pgc"] == "none" and t["tcr"] == "fh")
    pg = [o for o in res if o["key"] in ("R00-ic-rh-ibc-real-none", "R00-ic-rh-ibc-ideal-none")]
    short = lambda o, tag="opt": None if o is None else dict(
        key=o["key"], name=o["name"], PR=round(o[tag]["PR"], 2), eff_el=round(o[tag]["eff_el"], 4), W_air=round(o[tag]["W_air"], 2),
        EGT_C=round(o[tag]["EGT_C"], 1), T_rec_in_C=None if o[tag]["T_rec_in_C"] is None else round(o[tag]["T_rec_in_C"], 1),
        added_lpbf=o[tag]["parts"]["lpbf"]["added"], V_added_m3_lpbf=round(o[tag]["parts"]["lpbf"]["V_added_m3"], 2),
        mass_hx_kg_lpbf=round(o[tag]["parts"]["lpbf"]["mass_hx_kg"]), mass_hx_kg_foil=round(o[tag]["parts"]["foil"]["mass_hx_kg"]),
        limits_ok=o[tag]["limits"]["all"], T41_for_G_cap_C=o[tag].get("T41_for_G_cap_C"), prm=o[tag]["prm"])
    comb = {}
    for tag_, o in (("рекуператор + подогрев", rec_rh), ("простой", simple), ("рекуператор + промохлаждение", icr)):
        if o is not None:
            comb[tag_] = combustor_conditions(o)
    w_need = A["P_shaft_kW"] / A["G_cap"]
    return dict(
        source="routes/gtu/organic/cycles/cycles_results.json", date=DATE,
        status="оценка уровня облика (0D-цикл); числа — из cycles_results.json",
        key_finding=f"при T41 {A['T41_C']:.0f} °C (класс U) мягкий допуск по расходу требует удельной работы ≥ {w_need:.1f} кДж/кг; "
                    "простой цикл KV-3 на пределе, рекуператор сам по себе её снижает (потери давления) и не проходит; проходит "
                    "рекуператор вместе с повторным подогревом (или с промохлаждением) — см. representative",
        representative=dict(simple=short(simple), recup_only_out_of_limits=short(rec_only), recup_reheat=short(rec_rh),
                            recup_reheat_fuelheat=short(rec_rh_fh), recup_intercool=short(icr), fuel_heat=short(fh),
                            pgc=[short(o) for o in pg]),
        engine_funnel=dict(
            target="routes/gtu/funnel/funnel_space.py (CAT, CONT, ARCH)",
            new_axes=[
                dict(axis="reheat", kind="CAT", values=["нет", "межтурбинная камера"],
                     note="повторный подогрев на одном валу: газ после первой ступени турбины снова до T41; поднимает удельную работу "
                          "и температуру газа перед рекуператором — делает рекуператор выполнимым при мягком допуске",
                     model="cycles.evaluate(topo с rh=True, PR, prm=dict(x_rh=…))"),
                dict(axis="x_rh", kind="CONT", range=[0.3, 0.7], note="доля ln(p4/p5) до подогрева; оптимум — см. representative.recup_reheat.prm"),
                dict(axis="fuel_preheat", kind="CAT", values=["нет", "до 400 °C выхлопом"],
                     note="дешёвый узел (+2 детали), прирост КПД — representative.fuel_heat против simple"),
                dict(axis="T_rec_max_C", kind="CAT", values=[700, 800],
                     note="700 — LIM gtu_explore (нержавеющая сталь); 800 — печатный Inconel 625 (оценка, нужна ползучесть) — trec800"),
                dict(axis="recup_core_tech", kind="CAT", values=list(HX_TECH),
                     note="масса сердцевины: фольга 0,1 мм против печатной стенки 0,35 мм — hx_wall_sensitivity; печать — для коллекторов "
                          "и корпуса, сердцевина — фольга или печать тонкой стенкой"),
            ],
            changes=[
                dict(axis="eps_rec", now="0,75–0,92", proposed="0,80–0,92 только вместе с reheat или промохлаждением; без них при "
                     "T41 1000 °C вариант вне допуска по расходу (T41_for_G_cap_C)"),
                dict(axis="pgain", now="RDC 0–0,10", proposed="идеальный предел EAP посчитан (pgc_ideal_PG); при реализации 25 % и "
                     "штрафе КПД турбины 2 п. п. выигрыш около нуля — весу архетипа RDC снизить", model="cycles.pgc_gain"),
            ],
            models=["cycles.evaluate — цикл любой топологии на одном валу (Cantera)", "cycles.size_hx — объём, масса, модули печати",
                    "cycles.exergy_account — разрушение эксергии по узлам двумя способами"],
            drop=["паровой риформинг: вода в пакете (Plug&Play — нет), выигрыш в основном от впрыска пара",
                  "риформинг рециркулятом выхлопа: кислород обеднённого выхлопа сжигает метан, теплотворность не растёт",
                  "IBC: колесо компрессора IBC не помещается на вал 14 950 об/мин (нужен отдельный вал), охладитель громоздкий",
                  "горение с ростом давления: TRL 3, реальный выигрыш около нуля"],
            pgc_ideal_PG={o["key"]: o["opt"].get("PG_ideal") for o in pg},
        ),
        combustor_funnel=dict(
            target="routes/gtu/combustor (targets.json)",
            conditions=comb,
            new_axes=[
                dict(axis="T3_combustor_inlet_C", values=sorted({round(r["T_in_C"], -1) for v in comb.values() for r in v
                                                                 if r["node"] == "Камера сгорания"}),
                     note="рекуперированный цикл — горячий воздух на входе в камеру: проскок пламени, NOx"),
                dict(axis="reheat_combustor", values=[r for r in comb.get("рекуператор + подогрев", []) if r["node"] == "Камера подогрева"],
                     note="межтурбинная камера в обеднённом кислородом газе (вход — после ступени турбины)"),
                dict(axis="fuel_T_C", values=[15, 400], note="подогрев топлива выхлопом; выше 400 °C — разложение метана на стенках"),
            ],
            models=["cycles.burn — полное сгорание любой смеси CH4/H2/CO/H2O/CO2", "cycles.pgc_gain — предел горения с ростом давления",
                    "cycles.carbon_activity — сажеобразование в подогретом или риформированном топливе"],
        ),
        aws=dict(jobs="routes/gtu/organic/cycles/aws_jobs_cycles.json",
                 note="мелкая сетка (CYCLES_FINE): 864 топологии, π через 0,5; проверено командой plan без AWS; выполнимо после "
                      "попадания файлов в ветку claude/lucid-allen-yvps5d"),
    )


def write_feed(fd):
    """своя копия + слияние под ключом "cycles" в organic/feed_to_funnels.json (чужие ключи не трогаются)"""
    own = os.path.join(OUT_DIR, "feed_cycles.json")
    json.dump(fd, open(own, "w", encoding="utf-8"), ensure_ascii=False, indent=1, default=_jdef)
    if QUICK or OUT_DIR != HERE:
        return own, None
    path = os.path.join(os.path.dirname(HERE), "feed_to_funnels.json")
    d = {}
    if os.path.exists(path):
        try:
            d = json.load(open(path, encoding="utf-8"))
        except Exception:  # noqa: BLE001  чужая запись не дописана — не затираем
            return own, "organic/feed_to_funnels.json не читается — записана только своя копия"
    d["cycles"] = fd
    tmp = path + ".tmp_cycles"
    json.dump(d, open(tmp, "w", encoding="utf-8"), ensure_ascii=False, indent=1, default=_jdef)
    os.replace(tmp, path)
    return own, path


def feed_only():
    """пересобрать feed_to_funnels и графики из готового cycles_results.json (без перебора)"""
    path = os.path.join(OUT_DIR, "cycles_results.json")
    out = json.load(open(path, encoding="utf-8"))
    out["feed_to_funnels"] = feed_from(out)
    json.dump(out, open(path, "w", encoding="utf-8"), ensure_ascii=False, indent=1, default=_jdef)
    print("feed:", write_feed(out["feed_to_funnels"]))
    tech = "lpbf"
    res = out["results"]
    rows = [f for o in res for f in [flat(o, "opt", tech)] if f]
    sel = [r for r in rows if r["limits_ok"] and not r["bound"] and not r["degenerate"]]
    front_i = pareto(sel, [("eff_el", 1), ("added", -1), ("V_added_m3", -1)])
    base = [o for o in res if o["key"] == "R00-ic-rh-ibc-none-none"][0]
    plots(sel + [r for r in rows if not r["limits_ok"] and not r["bound"]], front_i, base["pr135"]["eff_el"])
    items = [("KV-3 простой, π опт.", base["opt"])]
    for r in out["pareto_lpbf"][:8]:
        o = [x for x in res if x["key"] == r["key"]][0]
        if o["key"] != base["key"]:
            items.append((f"{o['name']} (π {r['PR']:.1f})"[:70], o["opt"]))
    plots_exergy(items)


if __name__ == "__main__":
    if "--feed-only" in sys.argv:
        feed_only()
    else:
        main()
