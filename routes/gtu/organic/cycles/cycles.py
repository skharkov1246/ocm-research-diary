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
    p_fuel_k=1.25, dT_app=40.0, approach_eq=15.0, T_ref_max=850.0, T_fh_max=400.0, dp_ref_g=0.02, dp_ref_f=0.05,
    pinch=20.0, sc=(1.5, 2.0, 3.0), R_egr=(0.5, 1.0, 2.0, 4.0), eta_rc=0.78, eta_pump=0.6,
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


def T_sat(P):
    _W.PQ = P, 0.0
    return _W.T


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
    _, prod = burn(air.cp(), fuel_Y, fuel_h, 1500.0, air.P, lhv_f) if False else (None, None)
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
    gs = g.standard_gibbs_RT
    _GR.TP = T, ct.one_atm
    gc = _GR.standard_gibbs_RT[0]
    pa = g.X * P / ct.one_atm
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
def feed_steam(sc):
    """смесь CH4 + вода (жидкая при T0) на кг смеси: состав, энтальпия, массовая доля CH4"""
    mC, mW = MW[I["CH4"]], sc * MW[I["H2O"]]
    yC = mC / (mC + mW)
    Y = np.zeros(len(SP))
    Y[I["CH4"]], Y[I["H2O"]] = yC, 1 - yC
    h = yC * hs(St(1, T0, P0, Y_CH4))[0] + (1 - yC) * hs(St(1, T0, P_SAT0, Y_of("H2O:1"), liq=True))[0]
    return Y, h, yC


def cold_curve_steam(sc, P, T_ref):
    """составная кривая холодной стороны парового риформинга на кг смеси: [(T, Q от входа)]; равновесие с подходом"""
    Y, h_in, yC = feed_steam(sc)
    Ts = T_sat(P)
    yW = 1 - yC
    hC = lambda T: hs(St(1, T, P, Y_CH4))[0]
    hW_liq = lambda T: hs(St(1, T0, P_SAT0, Y_of("H2O:1"), liq=True))[0] + 4180.0 * (T - T0) * (1 + 0.0007 * (T - T0))
    pts = []
    for T in np.linspace(T0, Ts, 8):
        pts.append((T, yC * hC(T) + yW * hW_liq(T) - h_in))
    for T in np.linspace(Ts + 0.01, T_ref, 25):
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
    T_top = min(gin.T - A["dT_app"], A["T_ref_max"])
    if tcr == "fh":
        T = max(min(T_top, A["T_fh_max"]), T0)
        h0 = hs(St(1, T0, p_fuel, Y_CH4))[0]
        h1 = hs(St(1, T, p_fuel * (1 - A["dp_ref_f"]), Y_CH4))[0]
        return dict(Y=Y_CH4, h=h1, T=T, lhv=LHV_CH4, yC=1.0, q=h1 - h0, info=dict(T_fuel_C=T - 273.15),
                    feeds=[("CH4", 1.0, St(1, T0, p_fuel, Y_CH4))], Pout=p_fuel * (1 - A["dp_ref_f"]))
    if tcr == "steam":
        sc = prm["sc"]
        Yf, h_in, yC = feed_steam(sc)
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
                    aC=carbon_activity(Ye, T, Pout), water_per_CH4_kg=(1 - yC) / yC, cold_pts=pts)
        feeds = [("CH4", yC, St(1, T0, p_fuel, Y_CH4)), ("H2O_liq", 1 - yC, St(1, T0, p_fuel, Y_of("H2O:1"), liq=True))]
        return dict(Y=Ye, h=h_out, T=T, lhv=L, yC=yC, q=h_out - h_in, info=info, feeds=feeds, Pout=Pout)
    if tcr == "egr":
        R = prm["R_egr"]
        # рециркулят: газ из трубы (stack), сжатие до давления топлива
        rc = stack.cp(m=1.0)
        Mrc = 1.0 / np.sum(rc.Y / MW)
        mC, mR = MW[I["CH4"]], R * Mrc
        yC = mC / (mC + mR)
        rc_c = rc.cp()
        w_rc = poly(rc_c, p_fuel, A["eta_rc"], True)
        feed = mix([St(yC, T0, p_fuel, Y_CH4), rc_c.cp(m=1 - yC)], p_fuel)
        h_in = hs(feed)[0]
        T = max(T_top, feed.T)
        Pout = p_fuel * (1 - A["dp_ref_f"])
        Ye, _ = equil(feed.Y, max(T - A["approach_eq"], feed.T), p_fuel)
        h_out = hs(St(1, T, Pout, Ye))[0]
        if h_out < h_in:                      # экзотермика (горит кислород рециркулята): реактор адиабатный
            Ye, T = equil(feed.Y, h_in, Pout, "HP")
            h_out = h_in
        L = lhv(Ye)
        nC0 = yC / MW[I["CH4"]]
        X = 1 - (Ye[I["CH4"]] / MW[I["CH4"]]) / nC0
        g = gas()
        g.TPY = feed.T, feed.P, feed.Y
        Xf = g.X
        info = dict(T_ref_C=T - 273.15, X_CH4=X, R_egr=R, O2_per_CH4=Xf[I["O2"]] / Xf[I["CH4"]],
                    H2O_per_CH4=Xf[I["H2O"]] / Xf[I["CH4"]], CO2_per_CH4=Xf[I["CO2"]] / Xf[I["CH4"]],
                    lhv_gain=L / yC / LHV_CH4_CT - 1, aC=carbon_activity(Ye, T, Pout), w_rc_per_kg_rc=w_rc,
                    T_rc_C=rc_c.T - 273.15)
        feeds = [("CH4", yC, St(1, T0, p_fuel, Y_CH4)), ("RC", 1 - yC, rc_c.cp(m=1.0))]
        return dict(Y=Ye, h=h_out, T=T, lhv=L, yC=yC, q=h_out - h_in, info=info, feeds=feeds, Pout=Pout,
                    w_rc=w_rc * (1 - yC), rc_in=rc)
    raise ValueError(tcr)


# ======================================================================== один вариант цикла
DEFAULT_PRM = dict(x_rh=0.5, sc=2.0, R_egr=1.0)


def evaluate(topo, PR, prm=None, detail=False):
    """topo: rec (0 или степень регенерации), ic, rh, ibc (bool), pgc ('none','real','ideal'),
    tcr ('none','fh','steam','egr'). Возвращает словарь показателей на 1 кг/с воздуха и при мощности 5922 кВт."""
    prm = {**DEFAULT_PRM, **(prm or {})}
    rec, ic, rh, ibc, pgc, tcr = (topo[k] for k in ("rec", "ic", "rh", "ibc", "pgc", "tcr"))
    T41 = A["T41_C"] + 273.15
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
        comps.append(("Вентилятор промохладителя", [], [], -Wf_ic, 0.0))
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
    chain = (1 + A["dp_ref_g"] if tcr != "none" else 1.0) * (1 + A["dp_rg"] if rec else 1.0)
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
            Pg = st5.P / (1 + A["dp_ref_g"])
            g_ref = St(st5.m, T_of_h(hg - mf * fuel_new["q"] / st5.m, Pg, st5.Y), Pg, st5.Y)
            g_in = g_ref
        else:
            fuel_new = fuel
        if rec:
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
            comps.append(("Компрессор рециркуляции", [rc_in], [feeds[1]], -mf * fuel["w_rc"], 0.0))
        if tcr == "steam":
            w_pump = mf * (1 - fuel["yC"]) * (p_fuel - P0) / 999.1 / A["eta_pump"]
            comps.append(("Насос питательной воды", [St(mf * (1 - fuel["yC"]), T0, P0, Y_of("H2O:1"), liq=True)], [feeds[1]],
                          -w_pump, 0.0))
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
    # --- рециркуляция: часть газа из трубы уходит назад в риформер
    if tcr == "egr":
        m_rc = mf * (1 - fuel["yC"])
        out = exh.cp(m=exh.m - m_rc)
        comps.append(("Отбор рециркулята", [exh.cp()], [out, exh.cp(m=m_rc)], 0.0, 0.0))
        exh = out
    W_rc = mf * fuel.get("w_rc", 0.0)
    W_mech = (1 - A["eta_m"]) * (W_t + W_ibc_t)
    comps.append(("Механические потери", [], [], -W_mech, 0.0))
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
    res["limits"] = dict(G_ok=W1 <= A["G_cap"] + 1e-9, Wexh_ok=W_exh <= A["W_exh_cap"] + 1e-9, EGT_ok=EGT <= A["EGT_cap"] + 1e-9,
                         rec_ok=(T_rec_in is None) or (T_rec_in - 273.15 <= A["T_rec_max"] + 1e-6))
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
             ("Вентилятор охладителя IBC", [], [], -Wfan, 0.0),
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
        Xb = T0 * Sgen if (ins or outs) else -W
        rows.append(dict(node=name, X_A_kW=Xa * W1 / 1e3, X_B_kW=Xb * W1 / 1e3, energy_residual_kW=dH * W1 / 1e3))
    st = res["_streams"]
    E_fuel = sum(f[1] * res["mf"] * ex(f[2]) for f in st["fuel"]["feeds"] if f[0] != "RC")
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
    tau_w = -np.gradient(u, h)[-1]           # μ = 1, um = 1, b = 1
    fRe = tau_w / 0.5 * 4                    # f = τw/(0,5 um²) при ρ = 1; Re = 4b um/ν
    # Nu_H: θ'' = u·C, θ(±1) = 0 (симметрия), C из баланса
    N = n - 1
    D2 = (np.diag(-2 * np.ones(N)) + np.diag(np.ones(N - 1), 1) + np.diag(np.ones(N - 1), -1)) / h ** 2
    ui = u[1:-1]
    th = np.linalg.solve(D2, ui)             # θ'' = u (масштаб)
    thf = np.concatenate([[0], th, [0]])
    q_w = np.gradient(thf, h)[-1]            # поток на стенке (масштаб)
    Tb = np.trapezoid(u * thf, y) / np.trapezoid(u, y)
    NuH = abs(q_w * 4 / (0 - Tb))            # Nu = q Dh/(k (Tw - Tb)), Dh = 4b
    # Nu_T: θ'' = -λ u θ, наименьшее собственное значение
    from scipy.linalg import eigh
    M = np.diag(ui)
    w, v = eigh(-D2, M)
    lam = w[0]
    thv = np.concatenate([[0], v[:, 0], [0]])
    qv = np.gradient(thv, h)[-1]
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


def murray_tree(m, rho, mu, Dh, A_free, v_root=30.0, lam_LD=6.0, exponent=3.0, K_branch=0.3):
    """разветвлённый коллектор от газохода до каналов сердцевины по закону Мюррея (D_родитель^a = 2·D_ребёнок^a;
    a = 3 — ламинар, Murray 1926; 7/3 — турбулентный, Uylings 1977 — по выдержке поиска). Длина ветви — lam_LD·D.
    Возвращает уровни, диаметры, Re, потерю давления и объём дерева (оценка)."""
    D0 = math.sqrt(4 * m / (math.pi * rho * v_root))
    N_ch = A_free / (math.pi * Dh ** 2 / 4)
    levels = max(1, math.ceil(math.log2(N_ch)))
    rows, dp, vol = [], 0.0, 0.0
    D = D0
    stop = None
    for k in range(levels + 1):
        n = 2 ** k
        u = m / (n * rho * math.pi * D ** 2 / 4)
        Re = rho * u * D / mu
        f = 64 / Re if Re < 2300 else 0.316 * Re ** -0.25
        L = lam_LD * D
        dpk = f * L / D * rho * u ** 2 / 2 + K_branch * rho * u ** 2 / 2
        dp += dpk
        vol += n * math.pi * D ** 2 / 4 * L
        rows.append(dict(level=k, branches=n, D_mm=D * 1e3, u=u, Re=Re, dP_Pa=dpk))
        if D <= 4 * Dh:
            stop = k
            break
        D = D * 2 ** (-1.0 / exponent)
    return dict(D_root_m=D0, N_channels=N_ch, levels_used=len(rows), end_D_mm=rows[-1]["D_mm"], dP_Pa=dp, volume_m3=vol,
                rows=rows[:3] + rows[-2:], exponent=exponent)


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
        items += [("риформер с катализатором: модули", mods("риформер")), ("компрессор рециркуляции + двигатель", 2),
                  ("клапан", 1), ("газоход рециркуляции", 1)]
    dwheels = (res["nc"] + res["nt"]) - 5
    added = sum(n for _, n in items) + dwheels
    trl = min([TRL["simple"]] + [TRL[k] for k, on in (("rec", t["rec"]), ("ic", t["ic"]), ("rh", t["rh"]),
                                                    ("ibc", t["ibc"]), ("pgc", t["pgc"] != "none"))] + ([TRL[t["tcr"]]] if t["tcr"] != "none" else []) if True else [])
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


def topologies():
    L = []
    for rec in (0.0, 0.70, 0.80, 0.85, 0.90, 0.92):
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


def optimize_topo(t):
    """оптимум π (сетка + уточнение) и параметров; плюс точка «как есть» π 13,5"""
    best = None
    for prm in prm_grid(t):
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
            ref = minimize_scalar(lambda x: -score(_safe(t, x, prm))[1] - 10 * score(_safe(t, x, prm))[0],
                                  bounds=(lo, hi), method="bounded", options=dict(xatol=0.05, maxiter=12))
            rr = _safe(t, float(ref.x), prm)
            if score(rr) > score(rb):
                PRb, rb = float(ref.x), rr
        if best is None or score(rb) > score(best[1]):
            best = (PRb, rb, prm)
    fixed = None
    for prm in prm_grid(t):
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
    return clean


# ======================================================================== исследования (таблицы)
def study_recup_tradeoff(PR=5.0):
    """компромисс «степень регенерации — потеря давления» (Бежан: минимум производства энтропии = максимум КПД
    при заданном объёме): КПД, объём и масса печатной сердцевины"""
    rows = []
    base = dict(A)
    for eps in (0.70, 0.80, 0.85, 0.90, 0.92):
        for k in (0.5, 1.0, 1.5, 2.0):
            A["dp_ra"], A["dp_rg"] = base["dp_ra"] * k, base["dp_rg"] * k
            r = _safe(dict(rec=eps, ic=False, rh=False, ibc=False, pgc="none", tcr="none"), PR, None, detail=True)
            if r.get("ok"):
                hx = hx_list(r, "lpbf")[0]
                rows.append(dict(eps=eps, dp_total_pct=100 * (A["dp_ra"] + A["dp_rg"]), eff_el=r["eff_el"], W_air=r["W_air"],
                                 EGT_C=r["EGT_C"], V_core_m3=hx["V_core_m3"], mass_kg=hx["mass_kg"], L_m=hx["L_m"],
                                 A_front_m2=hx["A_front_m2"]))
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
    ax.scatter([r["added"] for r in bad], [100 * r["eff_el"] for r in bad], s=8, c="#c9c9c9", label="вне мягкого допуска")
    sc_ = ax.scatter([r["added"] for r in ok], [100 * r["eff_el"] for r in ok], s=[10 + 25 * r["V_added_m3"] for r in ok],
                     c=[r["TRL"] for r in ok], cmap="viridis", vmin=3, vmax=9, alpha=0.75, label="в допуске (размер — объём)")
    fr = [rows[i] for i in front]
    ax.scatter([r["added"] for r in fr], [100 * r["eff_el"] for r in fr], s=60, facecolors="none", edgecolors="#d62728",
               linewidths=1.4, label="фронт Парето")
    ax.axhline(100 * ref, color="#555", lw=0.8, ls="--")
    ax.text(ax.get_xlim()[0] + 0.3, 100 * ref + 0.3, "KV-3 простой, π 13,5", fontsize=8, color="#555")
    plt.colorbar(sc_, ax=ax, label="TRL (оценка)")
    ax.set_xlabel("добавочные детали и узлы (относительно KV-3 простого)")
    ax.set_ylabel("КПД на клеммах, %")
    ax.set_title("KV-3: перебор циклов (печатные теплообменники LPBF), оптимум π")
    ax.legend(fontsize=8, loc="lower right")
    fig.tight_layout()
    fig.savefig(os.path.join(HERE, "cycles_pareto.png"))
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
    with ProcessPoolExecutor(max_workers=4) as exe:
        res = list(exe.map(optimize_topo, L, chunksize=2))
    print(f"перебор: {time.time() - t0:.0f} с")
    tech = "lpbf"
    rows = []
    for o in res:
        for tag in ("opt", "pr135"):
            for tc in HX_TECH:
                f = flat(o, tag, tc)
                if f:
                    rows.append(f)
    sel = [r for r in rows if r["point"] == "opt" and r["tech"] == tech and r["limits_ok"] and not r["bound"]]
    front_i = pareto(sel, [("eff_el", 1), ("added", -1), ("V_added_m3", -1)])
    front = sorted([sel[i] for i in front_i], key=lambda r: (r["added"], -r["eff_el"]))
    sel_p = [r for r in sel if r["pnp"] != "нет" and r["TRL"] >= 5]
    front_p = sorted([sel_p[i] for i in pareto(sel_p, [("eff_el", 1), ("added", -1), ("V_added_m3", -1)])],
                     key=lambda r: (r["added"], -r["eff_el"]))
    base = [o for o in res if o["key"] == "R00-ic-rh-ibc-none-none"][0]
    ref_eff = base["pr135"]["eff_el"]
    # исследования
    trade = study_recup_tradeoff(5.0)
    refm = study_reformer()
    # коллектор по Мюррею для рекуператора лучшего варианта с рекуператором из фронта
    rec_best = max([r for r in sel if r["key"].startswith("R") and not r["key"].startswith("R00") and r["pnp"] != "нет"],
                   key=lambda r: r["eff_el"] - 0.002 * r["added"], default=None)
    murray = None
    if rec_best:
        o = [x for x in res if x["key"] == rec_best["key"]][0]
        rr = _safe(o["topo"], o["opt"]["PR"], o["opt"]["prm"], detail=True)
        st = rr["_streams"]
        hx = [h for h in hx_list(rr, tech) if h["name"] == "рекуператор"][0]
        W1 = rr["W_air"]
        gin = st["g_in"]
        rho, mu, _, _ = tprops(gin.T, gin.P, gin.Y)
        Afree = hx["A_m2"] * HX_TECH[tech]["Dh"] / (4 * hx["L_m"])
        murray = dict(key=rec_best["key"], gas_lam=murray_tree(gin.m * W1, rho, mu, HX_TECH[tech]["Dh"], Afree, exponent=3.0),
                      gas_turb=murray_tree(gin.m * W1, rho, mu, HX_TECH[tech]["Dh"], Afree, exponent=7.0 / 3.0),
                      dP_budget_gas_Pa=(1 - A["dp_core_share"]) * A["dp_rg"] * gin.P)
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
        runtime_s=time.time() - t0,
    )
    out["feed_to_funnels"] = feed(out, rows, front, front_p)
    json.dump(out, open(os.path.join(HERE, "cycles_results.json"), "w", encoding="utf-8"), ensure_ascii=False, indent=1,
              default=lambda x: float(x) if isinstance(x, (np.floating,)) else (x.tolist() if isinstance(x, np.ndarray) else str(x)))
    cols = list(rows[0].keys())
    with open(os.path.join(HERE, "cycles.csv"), "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=cols)
        w.writeheader()
        for r in rows:
            w.writerow({k: (round(v, 5) if isinstance(v, float) else v) for k, v in r.items()})
    plots(sel + [r for r in rows if r["point"] == "opt" and r["tech"] == tech and (not r["limits_ok"]) and not r["bound"]],
          list(range(len(sel)))[:0] + front_i, ref_eff)
    print(f"готово за {time.time() - t0:.0f} с; фронт Парето: {len(front)}")
    for r in front:
        print(f"  {r['added']:3d} дет. {100 * r['eff_el']:5.1f} % V {r['V_added_m3']:5.2f} м³ π {r['PR']:5.2f} TRL {r['TRL']} "
              f"P&P {r['pnp']:8s} | {r['name']}")


def feed(out, rows, front, front_p):
    """что подать в воронку двигателя (routes/gtu/funnel) и воронку камеры (routes/gtu/combustor)"""
    return dict(
        source="routes/gtu/organic/cycles/cycles_results.json",
        engine_funnel=dict(
            new_axes=[
                dict(axis="cycle_topology", values=["простой", "рекуператор 0,85–0,92 при π 4–6", "рекуператор при π 13,5 (как есть)",
                                                    "подогрев (межтурбинная камера)", "рекуператор + подогрев", "подогрев топлива"],
                     model="cycles.evaluate(topo, PR, prm) — КПД, расход, выхлоп, температура газа на входе в рекуператор"),
                dict(axis="recuperator_eps", values=[0.0, 0.80, 0.85, 0.90], model="cycles.size_hx — объём, масса, модули LPBF"),
                dict(axis="recuperator_dp_total_pct", values=[3.0, 6.0, 9.0], model="cycles.study_recup_tradeoff"),
                dict(axis="PR_for_recuperated", values=[4.0, 4.5, 5.0, 6.0], model="одна центробежная ступень до π 5"),
                dict(axis="hx_technology", values=list(HX_TECH.keys()), model="cycles.HX_TECH"),
            ],
            constraints=["мягкий допуск пакета: расход воздуха ≤ 1,2·21,41 кг/с, выхлоп ≤ 1,2·21,77 кг/с, ≤ 610 °C",
                         "газ на входе в рекуператор ≤ LIM['T_rec_max'] (700 °C)",
                         "Plug&Play: без воды; добавочный объём теплообменников (×2 на коллекторы) ≤ 6 м³ — оценка"],
            drop=["паровой риформинг (вода в пакете)", "риформинг рециркулятом (кислород выхлопа сжигает метан)",
                  "IBC (колесо не помещается на вал, громоздкий охладитель)", "горение с ростом давления (TRL 3)"],
        ),
        combustor_funnel=dict(
            new_axes=[
                dict(axis="T3_combustor_inlet_C", values=sorted({round(r["T3r_C"], -1) for r in front}),
                     note="рекуперированный цикл — горячий воздух на входе в камеру; риск обратного проскока и NOx"),
                dict(axis="reheat_inlet", values=dict(O2_vol_pct="≈ 13–15 (газ после ступени турбины)", T_in_C="≈ 750–850"),
                     note="межтурбинная камера в обеднённом кислородом газе — новая задача для воронки камеры"),
                dict(axis="fuel_preheat_C", values=[15, 250, 400], note="подогрев топлива выхлопом; выше 400 °C — сажа на стенках"),
                dict(axis="fuel_composition", values=["CH4", "синтез-газ парового риформинга (H2 — см. reformer)"],
                     note="только для справки: вариант с водой в KV не рекомендуется"),
            ],
            models=["cycles.burn — полное сгорание любой смеси CH4/H2/CO/H2O/CO2", "cycles.pgc_gain — предел горения с ростом давления",
                    "cycles.carbon_activity — сажеобразование в подогретом топливе"],
        ),
    )


if __name__ == "__main__":
    main()
