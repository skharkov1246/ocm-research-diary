"""Типоразмер № 3 двигателя-заменителя: согласованная расчётная точка газогенератора.

Постановка (README трека, решения владельца 14.11.2026): мягкий допуск пакета, расчётная точка — верх
типоразмера № 3, Taurus 60 генераторный (gtu_family_results.json, сценарий «мягкий»): схема 1В-2Ц, класс U,
π 13,5, T41 1000 °C, расход 25,55 кг/с, вал двигателя = вход редуктора пакета 14 950 об/мин.
gtu_family подобрал колесо под собственную оптимальную частоту (17 640 об/мин); здесь — согласованная точка
при n = 14 950 об/мин:
  1) две центробежные ступени (первая — двусторонний вход): распределение π, размеры, треугольники
     скоростей, число лопаток, КПД ступени с поправкой за отход от оптимальной частоты;
  2) осевая неохлаждаемая турбина на том же валу: сколько ступеней нужно, чтобы снять работу обеих
     ступеней компрессора и полезную мощность;
  3) пересчёт цикла gtu_cycle.cycle с итоговыми КПД и сравнение с gtu_family;
  4) остальные исполнения типоразмера № 3: работает ли тот же газогенератор на их расходе и частоте вала.

Уровень — облик (0D-цикл + одномерная средняя линия ступеней). Поправка КПД — набор потерь Oh, Yoon, Chung
(1997) + трение в диффузоре (оценка), привязанный к принятой в треке формуле gtu_explore.comp_stage
при коэффициенте расхода 0,085. Все пределы и КПД узлов — оценка; нужна проверка SU2/CalculiX и эксперимент.
Вход: ../gtu_family_results.json, ../gtu_explore.py (LIM, comp_stage, eta_axial), ../gtu_family.py (ETA_GEN, COOL,
SCEN), ../platform/fea/materials.json и disk_results.json. Выход: frame3_point_results.json."""
import os
import sys
import json
import math
import datetime
import platform

HERE = os.path.dirname(os.path.abspath(__file__))
os.chdir(HERE)
GTU = os.path.dirname(HERE)
sys.path.insert(0, GTU)

import numpy as np                                    # noqa: E402
import scipy                                          # noqa: E402
import cantera as ct                                  # noqa: E402
from scipy.optimize import brentq, minimize_scalar    # noqa: E402
import gtu_cycle as C                                 # noqa: E402  (импорт меняет текущую папку)
import gtu_explore as E                               # noqa: E402
import gtu_family as F                                # noqa: E402

os.chdir(HERE)
OUT = os.path.join(HERE, "frame3_point_results.json")
R = E.RA

# ------------------------------------------------------------------ входные данные (файл и поле)
FAM_FILE = os.path.join(GTU, "gtu_family_results.json")
FEA_DIR = os.path.join(GTU, "platform", "fea")
SCEN = "мягкий"
DESIGN = ("Taurus 60", "генератор")


def load_inputs():
    fam = json.load(open(FAM_FILE))
    slot = next(s for s in fam["slots"] if (s["model"], s["variant"]) == DESIGN)
    pt = slot[SCEN]
    frame = next(f for f in fam["frames"][SCEN] if abs(f["G_frame"] - round(pt["G"], 1)) < 1e-9)
    members = [s for s in fam["slots"] if s[SCEN] and frame["G_min"] - 0.05 <= s[SCEN]["G"] <= frame["G_frame"] + 0.05]
    mat = json.load(open(os.path.join(FEA_DIR, "materials.json")))["ti64"]
    fea = {c["case"]["name"]: c["result"] for c in json.load(open(os.path.join(FEA_DIR, "disk_results.json")))["cases"]}
    return fam, slot, pt, frame, members, mat, fea


FAM, SLOT, PT, FRAME, MEMBERS, MAT, FEA = load_inputs()
assert len(MEMBERS) == FRAME["n"], (len(MEMBERS), FRAME)
PR = PT["PR"]                    # 13.5
T41_C = PT["T41"]                # 1000
G = PT["G"]                      # 25.55 кг/с воздуха
N = float(PT["n_out"])           # 14 950 об/мин
P_EL_KW = SLOT["P_kW"]           # 5670 кВт на клеммах
P_SH_KW = SLOT["P_shaft_kW"]     # 5922 кВт на валу двигателя
G_CAP = SLOT["G_air"] * (1 + F.SCEN[SCEN]["alpha"])
EGT_CAP = SLOT["EGT0"] + F.SCEN[SCEN]["dT"]
ETA_GEN, ETA_GB = 0.972, 0.985   # генератор, редуктор (gtu_family.ETA_GEN = их произведение)
assert abs(F.ETA_GEN - ETA_GEN * ETA_GB) < 1e-12
U_TI = E.LIM["U_ti"]             # 620 м/с, титан (README трека)
U_AX = E.LIM["U_ax"]             # 420 м/с, осевая неохлаждаемая лопатка
T_AX = E.LIM["T_ax_max"]         # 1050 °C
T41_MAX_U = max(F.COOL["U"]["T"])  # 1000 °C, класс U
X_LEAK = 0.015                   # как gtu_family._run
ETA_M = 0.99                     # gtu_cycle.cycle по умолчанию
DP_COMB, P_EXH = 0.04, 1.02      # gtu_cycle.cycle по умолчанию

# ------------------------------------------------------------------ принятые величины одномерного расчёта (оценка)
A = dict(
    nu=0.40,          # отношение диаметров втулки и периферии на входе: вал проходит через оба колеса
    k1=0.97,          # загромождение входа лопатками
    B2=0.08,          # загромождение выхода колеса (лопатки + след)
    spl=0.70,         # вес разделительной лопатки в эффективном числе лопаток (доля её длины)
    eps_wake=0.366,   # доля следа (Johnston & Dean), значение по умолчанию TurboFlow 0.1.18
    zeta_d=0.15,      # потеря на форму в лопаточном диффузоре и переводном канале, доля c2^2/2
    cf_d=0.005,       # коэффициент трения стенок диффузора
    r4r2=1.55,        # выход лопаточного диффузора / радиус колеса
    alpha4=40.0,      # угол потока на выходе диффузора от радиального, град
    c4c2=0.45,        # отношение скоростей на выходе и входе диффузора
    f_inc=0.6,        # потеря на угле атаки: доля энергии нормальной составляющей (Conrad: 0,5-0,7)
    phi_ref=0.085,    # коэффициент расхода, при котором comp_stage даёт свой КПД (gtu_explore.comp_stage)
    betas=(30.0, 35.0, 40.0, 45.0, 50.0),   # угол лопатки на выходе от радиального, град (диапазон перебора)
    kcm=(0.12, 0.45),  # cm2/U2 (диапазон перебора)
    lax=0.8,          # осевая длина колеса = lax*(r2 - r1) + b2
    m_thd=0.80,       # число Маха в горлах входа колеса (относительное) и лопаточного диффузора на расчётном режиме
    phi_t=0.50,       # турбина: ca/U
    psi_max=2.5,      # турбина: предел нагрузки ступени (gtu_explore.turb_axial)
    m2_max=1.2,       # турбина: предел числа Маха за сопловым аппаратом (трансзвуковая решётка)
    u_margin=0.10,    # проектный запас по окружной скорости колёс компрессора к пределу 620 м/с (19 % по напряжениям):
                      # горячая вторая ступень (около 400 °C) и сквозное отверстие (platform/fea/README.md)
)
U_DES = U_TI * (1 - A["u_margin"])


# ------------------------------------------------------------------ свойства воздуха (таблицы gtu_explore из Cantera)
def cp_air(T):
    return (E.h(T + 1.0) - E.h(T - 1.0)) / 2.0


def gam(T):
    cp = cp_air(T)
    return cp / (cp - R)


def mu_air(T):
    """вязкость воздуха по Сазерленду"""
    return 1.716e-5 * (T / 273.15) ** 1.5 * (273.15 + 110.4) / (T + 110.4)


def stat(T0, p0, c):
    T = E.T_h(E.h(T0) - 0.5 * c * c)
    p = p0 * math.exp((E.s0(T) - E.s0(T0)) / R)
    return T, p, p / (R * T)


def flux_sonic(T0, p0):
    g = gam(T0)
    return p0 * math.sqrt(g / (R * T0)) * (2 / (g + 1)) ** ((g + 1) / (2 * (g - 1)))


def inlet(T0, p0, m, area, alpha=0.0):
    """дозвуковой вход: скорость c при расходе m через площадь area, угол закрутки alpha"""
    ca = math.cos(alpha)
    c = 300.0
    for _ in range(40):
        T = E.T_h(E.h(T0) - 0.5 * c * c)
        c = math.sqrt(gam(T) * R * T)
    f = lambda v: stat(T0, p0, v)[2] * v * ca * area - m
    if f(c) < 0:
        return None
    v = brentq(f, 1.0, c, xtol=1e-7)
    return (v,) + stat(T0, p0, v)


def inducer(T01, p01, m_side, n):
    """диаметр входа по минимуму относительной скорости на периферии"""
    om = n * math.pi / 30
    nu = A["nu"]

    def w1t(d):
        s = inlet(T01, p01, m_side, A["k1"] * math.pi / 4 * d * d * (1 - nu * nu))
        return 1e9 if s is None else math.hypot(s[0], om * d / 2)
    d = minimize_scalar(w1t, bounds=(0.05, 1.5), method="bounded", options={"xatol": 1e-6}).x
    area = A["k1"] * math.pi / 4 * d * d * (1 - nu * nu)
    c, T, p, rho = inlet(T01, p01, m_side, area)
    r1t, r1h = d / 2, nu * d / 2
    r1 = math.sqrt((r1t ** 2 + r1h ** 2) / 2)
    return dict(D1t=d, D1h=nu * d, r1t=r1t, r1h=r1h, r1=r1, A1=area, c1=c, T1=T, p1=p, rho1=rho, U1t=om * r1t,
                W1t=math.hypot(c, om * r1t), M1t=math.hypot(c, om * r1t) / math.sqrt(gam(T) * R * T),
                b1t=math.atan(om * r1t / c), b1=math.atan(om * r1 / c), b1h=math.atan(om * r1h / c))


# ------------------------------------------------------------------ потери ступени
def losses(g, f, m_side, dbl, friction=True):
    """Oh, Yoon, Chung (1997): угол атаки (Conrad), нагрузка лопаток (Coppage), трение (Jansen), зазор (Jansen),
    смешение (Johnston & Dean), рециркуляция (Oh), трение диска (Daily & Nece), перетекание (Aungier) — формулы
    сверены с TurboFlow 0.1.18, centrifugal_compressor/loss_model_oh.py. Добавлено (оценка): потеря на форму в
    лопаточном диффузоре и переводном канале zeta_d*c2^2/2, трение стенок диффузора 2*cf*(L/b2)*c_ср^2/2,
    угол атаки лопаток диффузора. Возврат — Дж/кг."""
    r1t, r1h, r1, r2, b2, Z = g["r1t"], g["r1h"], g["r1"], g["r2"], g["b2"], g["Zeff"]
    U2, cm2, cth2, c2, W2, a2 = f["U2"], f["cm2"], f["cth2"], f["c2"], f["W2"], f["alpha2"]
    W1t, W1h, W1, c1, c1m, cth1 = f["W1t"], f["W1h"], f["W1"], f["c1"], f["c1m"], f["cth1"]
    rho1, rho2, dh_eu = f["rho1"], f["rho2"], f["dh_eu"]
    Df = 1 - W2 / W1t + 0.75 * (dh_eu / U2 ** 2) / ((W1t / W2) * (Z / math.pi * (1 - r1t / r2) + 2 * r1t / r2))
    L = dict(Df=Df)
    L["inc"] = A["f_inc"] * 0.5 * (W1 * math.sin(f["inc"])) ** 2
    L["bld"] = 0.05 * Df ** 2 * U2 ** 2
    cb1, cb2 = math.cos(g["b1b"]), math.cos(g["b2b"])
    Lb = math.pi / 8 * (2 * r2 - (r1t + r1h) - b2 + 2 * g["Lax"]) * (2 / (cb1 + cb2))
    Dh = 2 * r2 * cb2 / (Z / math.pi + 2 * r2 * cb2 / b2) + 0.5 * (r1t + r1h) / r2 * cb1 / (Z / math.pi + (r1t + r1h) / (r1t - r1h) * cb1)
    Wavg = (c1 + c2 + W1t + 2 * W1h + 3 * W2) / 8
    Cf = 0.0412 * (U2 * Dh * rho1 / f["mu1"]) ** (-0.1925)
    L["sf"] = 2 * Cf * Lb / Dh * Wavg ** 2
    L["cl"] = 0.6 * g["cl"] / b2 * cth2 * math.sqrt(max(0.0, 4 * math.pi / (b2 * Z) * (r1t ** 2 - r1h ** 2) / ((r2 - r1t) * (1 + rho2 / rho1)) * cth2 * c1m))
    L["mix"] = (A["eps_wake"] / (1 - A["eps_wake"])) ** 2 * 0.5 * cm2 ** 2
    L["dif"] = A["zeta_d"] * 0.5 * c2 ** 2
    am = 0.5 * (a2 + math.radians(A["alpha4"]))
    Ld = (A["r4r2"] - 1) * r2 / math.cos(am)
    L["dfr"] = (2 * A["cf_d"] * Ld / b2 * 0.5 * c2 ** 2 * (1 + A["c4c2"] ** 2) / 2) if friction else 0.0
    L["dinc"] = A["f_inc"] * 0.5 * (c2 * math.sin(a2 - g["a2v"])) ** 2
    L["rc"] = 8e-5 * math.sinh(3.5 * a2 ** 3) * Df ** 2 * U2 ** 2
    Re = rho2 * U2 * r2 / f["mu2"]
    fdf = 2.67 / Re ** 0.5 if Re < 3e5 else 0.0622 / Re ** 0.2
    # двусторонний открытый ротор: обе стороны диска облопачены, гладкой задней стороны нет -> трение диска 0
    L["df"] = 0.0 if dbl else fdf * 0.5 * (rho1 + rho2) * r2 ** 2 * U2 ** 3 / (4 * m_side)
    rbar, bbar = 0.5 * (r1 + r2), 0.5 * ((r1t - r1h) + b2)
    Lm = math.pi / 2 * math.sqrt(((r2 - r1) ** 2 + (g["Lax"] - 0.5 * b2) ** 2) / 2)
    dP = m_side * (r2 * cth2 - r1 * cth1) / (Z * rbar * bbar * Lm)
    ucl = 0.816 * math.sqrt(2 * max(dP, 0.0) / rho2)
    L["lk"] = rho2 * Z * g["cl"] * Lm * ucl * ucl * U2 / (2 * m_side)
    return L


INTERNAL = ("inc", "bld", "sf", "cl", "mix", "dif", "dfr", "dinc")
PARASITIC = ("rc", "df", "lk")


def finish(T01, L, dh_eu):
    w0 = dh_eu + sum(L[k] for k in PARASITIC)
    dh_is = dh_eu - sum(L[k] for k in INTERNAL)
    T02 = E.T_h(E.h(T01) + w0)
    PRm = math.exp((E.s0(E.T_h(E.h(T01) + dh_is)) - E.s0(T01)) / R)
    return dict(w0=w0, dh_is=dh_is, eta_is_model=dh_is / w0, PR_model=PRm,
                eta_p_model=R * math.log(PRm) / (E.s0(T02) - E.s0(T01)), T02=T02)


def pfleiderer(r1, r2, b1, b2):
    """число лопаток по Пфлейдереру; углы от меридионального направления"""
    return 6.5 * (r2 + r1) / (r2 - r1) * math.cos(0.5 * (b1 + b2))


# ------------------------------------------------------------------ ступень: проект
def design(T01, p01, m_side, PRs, eta_p, n, beta_deg, kcm, dbl, ind, zfix=None, friction=True):
    om = n * math.pi / 30
    b2b = math.radians(beta_deg)
    T02 = E.T_s(E.s0(T01) + R * math.log(PRs) / eta_p)
    w0 = E.h(T02) - E.h(T01)
    par, p02i = 0.02 * w0, p01 * PRs
    for _ in range(80):
        dh_eu = w0 - par
        zf = zfix or 12.0
        for _ in range(30):
            zeff = zf * (1 + A["spl"])
            sig = 1 - math.sqrt(math.cos(b2b)) / zeff ** 0.7          # Wiesner
            den = sig - kcm * math.tan(b2b)
            if den < 0.2 or dh_eu < 0.5 * w0:
                return None
            U2 = math.sqrt(dh_eu / den)
            r2 = U2 / om
            zn = zfix or max(6.0, pfleiderer(ind["r1"], r2, ind["b1"], b2b))
            if abs(zn - zf) < 1e-9:
                break
            zf = zn
        cm2 = kcm * U2
        cth2 = sig * U2 - cm2 * math.tan(b2b)
        c2, W2, a2 = math.hypot(cm2, cth2), math.hypot(cm2, cth2 - U2), math.atan2(cth2, cm2)
        T2 = E.T_h(E.h(T02) - 0.5 * c2 * c2)
        p2 = p02i * math.exp((E.s0(T2) - E.s0(T02)) / R)
        rho2 = p2 / (R * T2)
        b2 = m_side / (rho2 * cm2 * 2 * math.pi * r2 * (1 - A["B2"]))
        g = dict(r1t=ind["r1t"], r1h=ind["r1h"], r1=ind["r1"], r2=r2, b2=b2, Zfull=zf, Zeff=zeff, b1b=ind["b1"], b2b=b2b,
                 cl=max(0.3e-3, 0.0007 * 2 * r2), Lax=A["lax"] * (r2 - ind["r1"]) + b2, a2v=a2, sig=sig, A1=ind["A1"], dbl=dbl)
        f = dict(U2=U2, cm2=cm2, cth2=cth2, c2=c2, W2=W2, alpha2=a2, rho1=ind["rho1"], rho2=rho2, mu1=mu_air(ind["T1"]),
                 mu2=mu_air(T2), W1t=ind["W1t"], W1h=math.hypot(ind["c1"], om * ind["r1h"]), W1=math.hypot(ind["c1"], om * ind["r1"]),
                 c1=ind["c1"], c1m=ind["c1"], cth1=0.0, inc=0.0, dh_eu=dh_eu)
        L = losses(g, f, m_side, dbl, friction)
        par_new = sum(L[k] for k in PARASITIC)
        dh_imp = dh_eu - (L["inc"] + L["bld"] + L["sf"] + L["cl"])
        p02i_new = p01 * math.exp((E.s0(E.T_h(E.h(T01) + dh_imp)) - E.s0(T01)) / R)
        done = abs(par_new - par) < 1e-4 and abs(p02i_new - p02i) < 0.01
        par, p02i = par_new, p02i_new
        if done:
            break
    res = finish(T01, L, dh_eu)
    res.update(g=g, f=f, L=L, ind=ind, n=n, PRs=PRs, eta_p_in=eta_p, beta_deg=beta_deg, kcm=kcm, U2=U2, D2=2 * r2, b2=b2,
               Zfull=zf, phi=m_side / (p01 / (R * T01) * U2 * (2 * r2) ** 2), lam=res["w0"] / U2 ** 2,
               Mu2=U2 / math.sqrt(gam(T01) * R * T01), p02i=p02i, T2=T2, p2=p2, rho2=rho2, T01=T01, p01=p01, m_side=m_side)
    return res


def best(T01, p01, m_side, PRs, eta_p, n, dbl, ulim=U_TI):
    """перебор угла лопатки и cm2/U2 по наибольшему КПД модели; затем целое число лопаток"""
    ind = inducer(T01, p01, m_side, n)
    out = None
    for bb in A["betas"]:
        def neg(k, bb=bb):
            d = design(T01, p01, m_side, PRs, eta_p, n, bb, k, dbl, ind)
            return 1.0 if d is None else -d["eta_p_model"]
        k = minimize_scalar(neg, bounds=A["kcm"], method="bounded", options={"xatol": 1e-4}).x
        d = design(T01, p01, m_side, PRs, eta_p, n, bb, k, dbl, ind)
        if d is None or d["U2"] > ulim:
            continue
        if out is None or d["eta_p_model"] > out["eta_p_model"]:
            out = d
    if out is None:
        return None
    z = float(round(out["Zfull"]))
    bb = out["beta_deg"]
    neg = lambda k: (lambda d: 1.0 if d is None else -d["eta_p_model"])(design(T01, p01, m_side, PRs, eta_p, n, bb, k, dbl, ind, z))
    k = minimize_scalar(neg, bounds=A["kcm"], method="bounded", options={"xatol": 1e-4}).x
    return design(T01, p01, m_side, PRs, eta_p, n, bb, k, dbl, ind, z)


def stage(T01, p01, m_tot, PRs, n, dbl, ulim=U_DES):
    """КПД ступени = comp_stage (принятая формула трека, коэффициент расхода 0,085) + поправка набора потерь:
    КПД модели при n минус КПД модели той же ступени при частоте, где коэффициент расхода 0,085. Положительная
    поправка не засчитывается."""
    m_side = m_tot / 2 if dbl else m_tot
    cs = E.comp_stage(T01, p01, PRs, m_tot, "ti", dbl)
    eta = cs["eta"]
    for _ in range(12):
        act = best(T01, p01, m_side, PRs, eta, n, dbl, ulim)
        if act is None:
            return None
        n0 = n * math.sqrt(A["phi_ref"] / act["phi"])
        nref = brentq(lambda x: best(T01, p01, m_side, PRs, eta, x, dbl, math.inf)["phi"] - A["phi_ref"], 0.6 * n0, 1.6 * n0, xtol=0.1)
        ref = best(T01, p01, m_side, PRs, eta, nref, dbl, math.inf)
        d_raw = act["eta_p_model"] - ref["eta_p_model"]
        new = cs["eta"] + min(0.0, d_raw)
        if abs(new - eta) < 1e-6:
            eta = new
            break
        eta = new
    act = best(T01, p01, m_side, PRs, eta, n, dbl, ulim)
    # чувствительность: без трения стенок диффузора (только набор Oh + потеря на форму)
    nf = lambda d: design(T01, p01, m_side, PRs, eta, d["n"], d["beta_deg"], d["kcm"], dbl, d["ind"], d["Zfull"], False)["eta_p_model"]
    d_nofr = nf(act) - nf(ref)
    k = n / nref
    fam_pen = 1.5 * abs(math.log(k)) / math.log(1.14)       # правило gtu_family.best_for, п. КПД цикла
    return dict(act=act, ref=ref, cs=cs, eta_p=eta, nref=nref, d_raw=d_raw, d_used=min(0.0, d_raw), d_nofr=d_nofr,
                k=k, fam_pen=fam_pen, fam_valid=0.84 <= k <= 1.14, shift=eta - act["eta_p_model"])


def compressor(pi1, n=N, m=G, ulim=U_DES):
    s1 = stage(E.T0, E.P0, m, pi1, n, True, ulim)
    if s1 is None:
        return None
    s2 = stage(s1["act"]["T02"], E.P0 * pi1, m, PR / pi1, n, False, ulim)
    if s2 is None:
        return None
    w = s1["act"]["w0"] + s2["act"]["w0"]
    T3 = s2["act"]["T02"]
    return dict(pi1=pi1, pi2=PR / pi1, s1=s1, s2=s2, w=w, T3=T3, eta_eq=R * math.log(PR) / (E.s0(T3) - E.s0(E.T0)))


def select_split(ulim=U_DES):
    grid = [round(x, 2) for x in np.arange(3.4, 5.41, 0.2)]
    scan = []
    for p in grid:
        c = compressor(p, ulim=ulim)
        scan.append(dict(pi1=p, w_kJkg=None if c is None else c["w"] / 1e3, eta_eq=None if c is None else c["eta_eq"],
                         U2_1=None if c is None else c["s1"]["act"]["U2"], U2_2=None if c is None else c["s2"]["act"]["U2"]))
    ok = [s for s in scan if s["w_kJkg"] is not None]
    p0 = min(ok, key=lambda s: s["w_kJkg"])["pi1"]
    r = minimize_scalar(lambda p: (lambda c: 1e9 if c is None else c["w"])(compressor(p, ulim=ulim)), bounds=(p0 - 0.2, p0 + 0.2),
                        method="bounded", options={"xatol": 0.01})
    return round(float(r.x), 2), scan


# ------------------------------------------------------------------ ступень: нерасчётный режим (та же геометрия)
def stage_off(st, T01, p01, m_tot, n, alpha1=0.0):
    d = st["act"]
    g = d["g"]
    dbl = g["dbl"]
    m_side = m_tot / 2 if dbl else m_tot
    om = n * math.pi / 30
    s = inlet(T01, p01, m_side, g["A1"], alpha1)
    if s is None:
        return dict(ok=False, why="вход колеса заперт")
    c1, T1, p1, rho1 = s
    c1m, cth1 = c1 * math.cos(alpha1), c1 * math.sin(alpha1)
    U1, U1t, U1h = om * g["r1"], om * g["r1t"], om * g["r1h"]
    W1 = math.hypot(c1m, U1 - cth1)
    inc = math.atan2(U1 - cth1, c1m) - g["b1b"]
    W1t = math.hypot(c1m, U1t - cth1 * g["r1"] / g["r1t"])
    W1h = math.hypot(c1m, U1h - cth1 * g["r1"] / g["r1h"])
    T0r = E.T_h(E.h(T1) + 0.5 * W1 * W1)
    p0r = p1 * math.exp((E.s0(T0r) - E.s0(T1)) / R)
    choke_i = m_side / (st["A_th"] * flux_sonic(T0r, p0r))           # доля запирающего расхода горла входа
    if choke_i >= 1.0:
        return dict(ok=False, why="горло входа колеса заперто")
    U2, r2, b2, b2b = om * g["r2"], g["r2"], g["b2"], g["b2b"]
    rho2, p02i, T2 = d["rho2"], d["p02i"] * (p01 / d["p01"]), d["T2"]
    for _ in range(300):
        cm2 = m_side / (rho2 * 2 * math.pi * r2 * b2 * (1 - A["B2"]))
        cth2 = g["sig"] * U2 - cm2 * math.tan(b2b)
        dh_eu = U2 * cth2 - U1 * cth1
        if dh_eu <= 0:
            return dict(ok=False, why="нет работы")
        c2, W2, a2 = math.hypot(cm2, cth2), math.hypot(cm2, cth2 - U2), math.atan2(cth2, cm2)
        f = dict(U2=U2, cm2=cm2, cth2=cth2, c2=c2, W2=W2, alpha2=a2, rho1=rho1, rho2=rho2, mu1=mu_air(T1), W1t=W1t, W1h=W1h,
                 W1=W1, c1=c1, c1m=c1m, cth1=cth1, inc=inc, dh_eu=dh_eu, mu2=mu_air(T2))
        L = losses(g, f, m_side, dbl)
        w0 = dh_eu + sum(L[k] for k in PARASITIC)
        T02 = E.T_h(E.h(T01) + w0)
        T2 = E.T_h(E.h(T02) - 0.5 * c2 * c2)
        dh_imp = dh_eu - (L["inc"] + L["bld"] + L["sf"] + L["cl"])
        if dh_imp <= 0:
            return dict(ok=False, why="потери больше работы")
        p02i_new = p01 * math.exp((E.s0(E.T_h(E.h(T01) + dh_imp)) - E.s0(T01)) / R)
        rho2_new = p02i_new * math.exp((E.s0(T2) - E.s0(T02)) / R) / (R * T2)
        done = abs(rho2_new - rho2) < 1e-11 * rho2 and abs(p02i_new - p02i) < 1e-9 * p02i
        rho2 = 0.5 * (rho2 + rho2_new)
        p02i = p02i_new
        if done:
            break
    res = finish(T01, L, dh_eu)
    eta_p = res["eta_p_model"] + st["shift"]
    PRs = math.exp(eta_p * (E.s0(res["T02"]) - E.s0(T01)) / R)
    choke_d = m_side / (st["A_thd"] * flux_sonic(res["T02"], p02i))     # доля запирающего расхода горла диффузора
    return dict(ok=choke_d < 1.0, why=None if choke_d < 1.0 else "горло диффузора заперто", PRs=PRs, eta_p=eta_p,
                w0=res["w0"], T02=res["T02"], p02=p01 * PRs, inc_deg=math.degrees(inc), dinc_deg=math.degrees(a2 - g["a2v"]), choke_d=choke_d,
                choke_i=choke_i,
                eta_p_model=res["eta_p_model"], dh_eu=dh_eu)


def flux_M(T0, p0, M):
    g = gam(T0)
    return p0 * math.sqrt(g / (R * T0)) * M * (1 + 0.5 * (g - 1) * M * M) ** (-(g + 1) / (2 * (g - 1)))


def throats(st):
    """площади горла входа колеса (относительное движение, средний радиус) и горла лопаточного диффузора:
    число Маха в горле m_thd на расчётном режиме (оценка)"""
    d = st["act"]
    ind, om = d["ind"], d["n"] * math.pi / 30
    W1 = math.hypot(ind["c1"], om * ind["r1"])
    T0r = E.T_h(E.h(ind["T1"]) + 0.5 * W1 * W1)
    p0r = ind["p1"] * math.exp((E.s0(T0r) - E.s0(ind["T1"])) / R)
    st["A_th"] = d["m_side"] / flux_M(T0r, p0r, A["m_thd"])
    st["A_thd"] = d["m_side"] / flux_M(d["T02"], d["p02i"], A["m_thd"])


# ------------------------------------------------------------------ турбина
def turbine_stages(T41, p41, far, W41, w_t, n, z, U_mean):
    """осевая турбина: z одинаковых ступеней на постоянном среднем диаметре, ca = phi_t*U, выход без закрутки;
    КПД ступени (политропный) — gtu_explore.eta_axial(psi) (принято в треке, оценка)."""
    b = E.bfrac(far)
    Rb = E.Rg(b)
    Dm = 60 * U_mean / (math.pi * n)
    dh = w_t / z
    psi = dh / U_mean ** 2
    eta = E.eta_axial(psi)
    ca = A["phi_t"] * U_mean
    out, T01, p01 = [], T41, p41
    for k in range(z):
        h01 = E.h(T01, b)
        T03 = E.T_h(h01 - dh, b)
        ER = math.exp((E.s0(T01, b) - E.s0(T03, b)) / (eta * Rb))
        cth2 = psi * U_mean
        c2 = math.hypot(ca, cth2)
        dh_avail = h01 - E.h(E.T_s(E.s0(T01, b) - Rb * math.log(p01 / (C.P0 * P_EXH)), b), b)   # до давления выхлопа
        if 0.5 * c2 * c2 >= dh_avail:
            out.append(dict(stage=k + 1, T01_C=T01 - 273.15, psi=psi, reaction=1 - psi / 2, eta_p=eta, dh_kJkg=dh / 1e3, c2=c2,
                            M2=None, kinetic_over_available=0.5 * c2 * c2 / dh_avail, ER=ER,
                            note="скорость за соплом больше, чем даёт весь перепад до выхлопа: ступень невозможна"))
            T01, p01 = T03, p01 / ER
            continue
        T2 = E.T_h(h01 - 0.5 * c2 * c2, b)
        cp2 = (E.h(T2 + 1, b) - E.h(T2 - 1, b)) / 2
        a2 = math.sqrt(cp2 / (cp2 - Rb) * Rb * T2)
        w2 = math.hypot(ca, cth2 - U_mean)
        p2 = p01 * math.exp((E.s0(T2, b) - E.s0(T01, b)) / Rb)          # сопло без потерь (оценка)
        h2 = W41 / (p2 / (Rb * T2) * ca * math.pi * Dm)
        T3 = E.T_h(E.h(T03, b) - 0.5 * ca * ca, b)
        p03 = p01 / ER
        p3 = p03 * math.exp((E.s0(T3, b) - E.s0(T03, b)) / Rb)
        h3 = W41 / (p3 / (Rb * T3) * ca * math.pi * Dm)
        cp3 = (E.h(T3 + 1, b) - E.h(T3 - 1, b)) / 2
        out.append(dict(stage=k + 1, T01_C=T01 - 273.15, T03_C=T03 - 273.15, ER=ER, psi=psi, reaction=1 - psi / 2, eta_p=eta,
                        dh_kJkg=dh / 1e3, c2=c2, M2=c2 / a2, alpha2_deg=math.degrees(math.atan2(cth2, ca)), w2=w2, M_w2=w2 / a2,
                        T_rel_rotor_C=T2 + 0.5 * w2 * w2 / cp2 - 273.15, h_rotor_in_m=h2, h_exit_m=h3,
                        hub_tip_exit=(Dm - h3) / (Dm + h3), AN2_m2rpm2=math.pi * Dm * h3 * n * n,
                        M_w3=math.hypot(ca, U_mean) / math.sqrt(cp3 / (cp3 - Rb) * Rb * T3), kinetic_over_available=0.5 * c2 * c2 / dh_avail))
        T01, p01 = T03, p03
    ER_tot = p41 / p01
    checks = dict(psi_ok=psi <= A["psi_max"], reaction_ok=1 - psi / 2 >= 0.0,
                  M2_ok=all(s["M2"] is not None and s["M2"] <= A["m2_max"] for s in out),
                  U_ok=U_mean <= U_AX + 1e-9, T_ok=T41 - 273.15 <= T_AX)
    return dict(z=z, Dm=Dm, U_mean=U_mean, n=n, psi=psi, eta_p=eta, ER_total=ER_tot, T_exit_C=T01 - 273.15,
                stages=out, checks=checks, feasible=all(checks.values()))


def poly_air(T1, p1, p2, eta):
    a = C.gas_air()
    a.TPX = T1, p1, C.AIR
    w = C.poly_process(a, p2, eta, True)
    return w, a.T


def eta_equiv(w_target, PRt):
    """политропный КПД всего сжатия в gtu_cycle, дающий ту же работу (Cantera)"""
    return brentq(lambda e: poly_air(C.T0, C.P0, C.P0 * PRt, e)[0] - w_target, 0.6, 0.99, xtol=1e-7)


def design_point(pi1):
    comp = compressor(pi1)
    for st in (comp["s1"], comp["s2"]):
        throats(st)
    # работа ступеней на Cantera (как в gtu_cycle) и эквивалентный КПД
    w1, T2c = poly_air(C.T0, C.P0, C.P0 * pi1, comp["s1"]["eta_p"])
    w2, T3c = poly_air(T2c, C.P0 * pi1, C.P0 * PR, comp["s2"]["eta_p"])
    eta_pc = eta_equiv(w1 + w2, PR)
    # турбина: КПД зависит от нагрузки, нагрузка — от работы; итерации с gtu_cycle.cycle
    opts = {}
    for z in (1, 2, 3):
        eta_t = 0.88
        for _ in range(20):
            cyc = C.cycle(PR, T41_C, eta_pc=eta_pc, eta_pt_gg=eta_t, eta_pt_pt=eta_t, x_ngv=0.0, x_rot=0.0, x_leak=X_LEAK,
                          dp_duct=0.0, W1=G)
            W41 = cyc["W_exh"]
            far = cyc["Wf"] / (G * (1 - X_LEAK))
            w_t = (cyc["wc_kJkg"] * 1e3 * G / ETA_M + cyc["P_MW"] * 1e6 / ETA_M) / W41
            tb = turbine_stages(T41_C + 273.15, C.P0 * PR * (1 - DP_COMB), far, W41, w_t, N, z, U_AX)
            if abs(tb["eta_p"] - eta_t) < 1e-7:
                break
            eta_t = tb["eta_p"]
        opts[z] = dict(turb=tb, cyc=cyc, w_t=w_t, far=far)
    z_sel = min(z for z in opts if opts[z]["turb"]["feasible"])
    return dict(comp=comp, w_stage_cantera=(w1, w2), T_stage_cantera=(T2c, T3c), eta_pc=eta_pc, opts=opts, z_sel=z_sel)


# ------------------------------------------------------------------ газогенератор на нерасчётном режиме (быстрая модель)
def comp_map(dp, W1, n, alpha1=0.0):
    c = dp["comp"]
    a = stage_off(c["s1"], E.T0, E.P0, W1, n, alpha1)
    if not a.get("PRs"):
        return a
    b = stage_off(c["s2"], a["T02"], a["p02"], W1, n)
    if not b.get("PRs"):
        return b
    return dict(ok=a["ok"] and b["ok"], why=a["why"] or b["why"], PR=a["PRs"] * b["PRs"], w=a["w0"] + b["w0"], T3=b["T02"],
                s1=a, s2=b)


def hot_section(T3, PRc, T41, W1, eta_t):
    """камера и турбина по таблицам gtu_explore: доля топлива, выход, мощность на валу (как gtu_cycle: eta_b, eta_m)"""
    f = E.far_for(T3, T41) / 0.995
    b = E.bfrac(f)
    W41 = W1 * (1 - X_LEAK) * (1 + f)
    p41 = PRc * E.P0 * (1 - DP_COMB)
    ER = p41 / (E.P0 * P_EXH)
    T5 = E.T_s(E.s0(T41, b) - eta_t * E.Rg(b) * math.log(ER), b)
    w_t = E.h(T41, b) - E.h(T5, b)
    return dict(f=f, W41=W41, p41=p41, ER=ER, T5=T5, w_t=w_t)


def gg_point(dp, n, T41C, alpha1=0.0):
    """рабочая точка одновального газогенератора при частоте n и T41: расход из совместной работы компрессора и
    турбины (горло соплового аппарата, закон конуса Стодолы)"""
    K0, z, n0 = dp["K0"], dp["z_sel"], N
    T41 = T41C + 273.15
    U = U_AX * n / n0
    eta_t = dp["opts"][z]["turb"]["eta_p"]

    def resid(W1, eta_t):
        cm = comp_map(dp, W1, n, alpha1)
        if not cm.get("PR"):
            return 1.0          # вход заперт или нет работы — это всегда при избытке расхода
        hs = hot_section(cm["T3"], cm["PR"], T41, W1, eta_t)
        return hs["W41"] * math.sqrt(T41) / hs["p41"] / math.sqrt(1 - hs["ER"] ** -2) / K0 - 1.0
    for _ in range(6):
        W1 = brentq(lambda w: resid(w, eta_t), 0.55 * G, 1.25 * G, xtol=1e-7)
        cm = comp_map(dp, W1, n, alpha1)
        hs = hot_section(cm["T3"], cm["PR"], T41, W1, eta_t)
        psi = hs["w_t"] / z / U ** 2
        new = E.eta_axial(psi)
        if abs(new - eta_t) < 1e-8:
            break
        eta_t = new
    P = ETA_M * hs["W41"] * hs["w_t"] - W1 * cm["w"]
    Wf = hs["f"] * W1 * (1 - X_LEAK)
    return dict(n=n, T41_C=T41C, alpha1_deg=math.degrees(alpha1), W1=W1, PR=cm["PR"], w_c=cm["w"], T3_C=cm["T3"] - 273.15,
                EGT_C=hs["T5"] - 273.15, P_kW=P / 1e3, eff=P / (Wf * C.LHV), eta_t=eta_t, psi=psi, ok=cm["ok"], why=cm["why"],
                inc1_deg=cm["s1"]["inc_deg"], eta_p1=cm["s1"]["eta_p"], eta_p2=cm["s2"]["eta_p"], PR1=cm["s1"]["PRs"], PR2=cm["s2"]["PRs"])


def surge_line(dp, n, alpha1=0.0):
    """граница устойчивости (оценка): наибольшая степень сжатия на линии постоянной частоты"""
    ws = np.linspace(0.60 * G, 1.15 * G, 56)
    allp = [(w, comp_map(dp, w, n, alpha1)) for w in ws]
    pts = [(w, c) for w, c in allp if c.get("PR")]
    w_s, c_s = max(pts, key=lambda t: t[1]["PR"])
    choke = [(w, c.get("why")) for w, c in allp if w > w_s and (not c.get("PR") or not c["ok"])]
    return dict(W_surge=w_s, PR_surge=c_s["PR"], W_choke=choke[0][0] if choke else None, choke_where=choke[0][1] if choke else None,
                scan_W=[ws[0], ws[-1]])


def solve_power(dp, n, P_kW, alpha1=0.0):
    f = lambda T: gg_point(dp, n, T, alpha1)["P_kW"] - P_kW
    lo, hi = 650.0, 1250.0
    if f(hi) < 0:
        return None
    return gg_point(dp, n, brentq(f, lo, hi, xtol=0.01), alpha1)


def confirm(dp, op):
    """та же точка по gtu_cycle.cycle (Cantera): второй способ"""
    eta_pc = eta_equiv(op["w_c"], op["PR"])
    r = C.cycle(op["PR"], op["T41_C"], eta_pc=eta_pc, eta_pt_gg=op["eta_t"], eta_pt_pt=op["eta_t"], x_ngv=0.0, x_rot=0.0,
                x_leak=X_LEAK, dp_duct=0.0, W1=op["W1"])
    return dict(P_kW=r["P_MW"] * 1e3, eff=r["eff"], EGT_C=r["EGT_C"], dP_pct=(r["P_MW"] * 1e3 / op["P_kW"] - 1) * 100,
                dEGT_K=r["EGT_C"] - op["EGT_C"])


# ------------------------------------------------------------------ прочность колёс (оценка подобием по примерам platform/fea)
def interp(tab, x):
    xs, ys = [t[0] for t in tab], [t[1] for t in tab]
    return float(np.interp(x, xs, ys)), (xs[0] <= x <= xs[-1])


def material_check(st):
    d = st["act"]
    T_in, T_out = d["T01"] - 273.15, d["T02"] - 273.15
    ys, in_range = interp(MAT["ys_MPa"], T_out)
    out = dict(U2=d["U2"], U_lim=U_TI, margin_U_pct=(1 - d["U2"] / U_TI) * 100, margin_stress_pct=(1 - (d["U2"] / U_TI) ** 2) * 100,
               T_in_C=T_in, T_out_C=T_out, ys_at_Tout_MPa=ys, ys_in_data_range=in_range)
    for case, key in (("comp_T60_1C_ti", "bore"), ("comp_T60_1C_ti_solid", "solid")):
        r = FEA[case]
        vm = r["vm_max_node_MPa"] * (d["U2"] / r["U_tip"]) ** 2
        out[f"vm_{key}_MPa"] = vm
        out[f"vm_{key}_to_ys"] = vm / ys
    return out


# ------------------------------------------------------------------ сборка результата
def sig(x, k=6):
    if isinstance(x, (bool, np.bool_)):
        return bool(x)
    if isinstance(x, str) or x is None:
        return x
    if isinstance(x, (int, np.integer)):
        return int(x)
    x = float(x)
    return 0.0 if x == 0 else float(f"{x:.{k}g}")


def clean(o):
    if isinstance(o, dict):
        return {k: clean(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [clean(v) for v in o]
    return sig(o)


def stage_report(st, label):
    d, r = st["act"], st["ref"]
    g, f, ind = d["g"], d["f"], d["ind"]
    dbl = g["dbl"]
    return dict(
        label=label, inlet="двусторонний" if dbl else "односторонний", PR=d["PRs"],
        m_total_kg_s=d["m_side"] * (2 if dbl else 1), m_per_side_kg_s=d["m_side"],
        T01_C=d["T01"] - 273.15, p01_kPa=d["p01"] / 1e3, T02_C=d["T02"] - 273.15,
        n_rpm=d["n"], D1_hub_m=ind["D1h"], D1_tip_m=ind["D1t"], D1t_D2=ind["D1t"] / d["D2"],
        W1t_m_s=ind["W1t"], M1t_rel=ind["M1t"], beta1_tip_deg=math.degrees(ind["b1t"]), beta1_rms_deg=math.degrees(ind["b1"]),
        beta1_hub_deg=math.degrees(ind["b1h"]),
        D2_m=d["D2"], b2_m=d["b2"], b2_D2=d["b2"] / d["D2"], U2_m_s=d["U2"], Mu2=d["Mu2"],
        phi=d["phi"], phi2_cm2_U2=d["kcm"], work_coeff_lambda=d["lam"], slip_wiesner=g["sig"],
        blades_full=int(g["Zfull"]), blades_splitter=int(g["Zfull"]), blades_eff_for_slip=g["Zeff"],
        beta2_blade_from_radial_deg=d["beta_deg"], beta2_blade_from_tangent_deg=90 - d["beta_deg"],
        alpha2_flow_from_radial_deg=math.degrees(f["alpha2"]), c2_m_s=f["c2"], W2_W1t=f["W2"] / ind["W1t"],
        tip_clearance_m=g["cl"], axial_length_m=g["Lax"], w_kJkg=d["w0"] / 1e3,
        eta_p=st["eta_p"], eta_is=None,
        eta_basis=dict(comp_stage_eta_p=st["cs"]["eta"], comp_stage_D2=st["cs"]["D2"], comp_stage_rpm=st["cs"]["rpm"],
                       comp_stage_U2=st["cs"]["U2"],
                       model_eta_p_at_n=d["eta_p_model"], model_eta_p_at_nref=r["eta_p_model"], n_ref_rpm=st["nref"],
                       ref_D2_m=r["D2"], ref_U2_m_s=r["U2"], ref_M1t_rel=r["ind"]["M1t"],
                       d_eta_raw=st["d_raw"], d_eta_used=st["d_used"], d_eta_without_diffuser_friction=st["d_nofr"],
                       k_speed=st["k"], gtu_family_rule_penalty_pt_cycle=st["fam_pen"], gtu_family_rule_in_range=st["fam_valid"]),
        losses_pct_of_work={k: 100 * v / d["w0"] for k, v in d["L"].items() if k != "Df"}, diffusion_factor_Df=d["L"]["Df"],
        inducer_throat_m2=st["A_th"], diffuser_throat_m2=st["A_thd"])


def eta_is_from_p(T01, PRs, eta_p):
    T02 = E.T_s(E.s0(T01) + R * math.log(PRs) / eta_p)
    T02s = E.T_s(E.s0(T01) + R * math.log(PRs))
    return (E.h(T02s) - E.h(T01)) / (E.h(T02) - E.h(T01))


def run():
    pi1, scan = select_split(U_DES)
    pi1_620, scan_620 = select_split(U_TI)          # для сравнения: без проектного запаса, только предел материала
    c620 = compressor(pi1_620, ulim=U_TI)
    dp = design_point(pi1)
    z = dp["z_sel"]
    cyc = dp["opts"][z]["cyc"]
    # пропускная способность турбины (для нерасчётных режимов) — той же быстрой моделью, что и gg_point
    hs = hot_section(dp["comp"]["T3"], PR, T41_C + 273.15, G, dp["opts"][z]["turb"]["eta_p"])
    dp["K0"] = hs["W41"] * math.sqrt(T41_C + 273.15) / hs["p41"] / math.sqrt(1 - hs["ER"] ** -2)
    s1r, s2r = stage_report(dp["comp"]["s1"], "ступень 1"), stage_report(dp["comp"]["s2"], "ступень 2")
    for sr, st in ((s1r, dp["comp"]["s1"]), (s2r, dp["comp"]["s2"])):
        sr["eta_is"] = eta_is_from_p(st["act"]["T01"], st["act"]["PRs"], st["eta_p"])
        sr["material"] = material_check(st)
    # сравнение с gtu_family
    fam = C.cycle(PR, T41_C, eta_pc=F.ARCH["1В-2Ц"]["eta_pc"], eta_pt_gg=0.865, eta_pt_pt=0.885, x_ngv=0.0, x_rot=0.0, x_leak=X_LEAK, W1=G)
    duct = C.cycle(PR, T41_C, eta_pc=dp["eta_pc"], eta_pt_gg=dp["opts"][z]["turb"]["eta_p"], eta_pt_pt=dp["opts"][z]["turb"]["eta_p"],
                   x_ngv=0.0, x_rot=0.0, x_leak=X_LEAK, dp_duct=0.01, W1=G)
    P_sh = cyc["P_MW"] * 1e3
    w_sp = P_sh / G
    tb_eta = dp["opts"][z]["turb"]["eta_p"]
    T41_need = brentq(lambda t: C.cycle(PR, t, eta_pc=dp["eta_pc"], eta_pt_gg=dp["opts"][z]["turb"]["eta_p"],
                                        eta_pt_pt=dp["opts"][z]["turb"]["eta_p"], x_ngv=0.0, x_rot=0.0, x_leak=X_LEAK,
                                        dp_duct=0.0, W1=G)["P_MW"] * 1e3 - P_SH_KW, 900.0, 1150.0, xtol=0.05)
    # чувствительность к степени сжатия при тех же КПД узлов (оценка: геометрия ступеней при этом другая)
    pr_sens = []
    for prx in (12.0, 12.5, 13.0, 13.5, 14.0, 14.5):
        r = C.cycle(prx, T41_C, eta_pc=dp["eta_pc"], eta_pt_gg=tb_eta, eta_pt_pt=tb_eta, x_ngv=0.0, x_rot=0.0, x_leak=X_LEAK,
                    dp_duct=0.0, W1=G)
        pr_sens.append(dict(PR=prx, P_shaft_kW=r["P_MW"] * 1e3, eff_shaft=r["eff"], EGT_C=r["EGT_C"]))
    # баланс мощности: компрессор и вал — Cantera (gtu_cycle), турбина — независимо, по таблицам gtu_explore
    tb = dp["opts"][z]["turb"]
    hs_chk = hot_section(dp["T_stage_cantera"][1], PR, T41_C + 273.15, G, tb["eta_p"])
    P_c = cyc["wc_kJkg"] * G
    P_c_stages = dp["comp"]["w"] * G / 1e3
    P_t_tab = hs_chk["W41"] * hs_chk["w_t"] / 1e3
    P_t_need = (P_c + P_sh) / ETA_M
    # члены типоразмера
    design_op = gg_point(dp, N, T41_C)
    members = []
    for s in MEMBERS:
        m = s[SCEN]
        row = dict(model=s["model"], variant=s["variant"], G_soft=m["G"], share_of_design_pct=100 * m["G"] / G,
                   arch=m["arch"], cool=m["cool"], T41_family=m["T41"], n_out=m.get("n_out"), P_shaft_kW=s["P_shaft_kW"],
                   G_air0=s["G_air"], G_cap=s["G_air"] * (1 + F.SCEN[SCEN]["alpha"]), EGT0=s["EGT0"], EGT_cap=s["EGT0"] + F.SCEN[SCEN]["dT"])
        if (s["model"], s["variant"]) == DESIGN:
            row["role"] = "расчётная точка"
        if m["arch"] != "1В-2Ц" or m["cool"] != "U":
            why = {"2В-1Ц": "двухвальный (газогенератор + свободная силовая турбина)", "2В-2Ц": "двухвальный, два колеса",
                   "1В-1Ц": "одновальный с одним колесом"}[m["arch"]]
            row["verdict"] = "нужен другой вариант"
            row["reason"] = f"схема {m['arch']}/{m['cool']}: {why}" + ("; охлаждаемые сопловой аппарат и рабочая лопатка" if m["cool"] == "K" else
                                                                       "; охлаждаемый сопловой аппарат" if m["cool"] == "N" else "")
            if m["arch"] == "1В-1Ц" and m.get("n_out"):
                row["U2_stage2_at_member_rpm"] = s2r["U2_m_s"] * m["n_out"] / N
            members.append(row)
            continue
        n = float(m["n_out"])
        op = solve_power(dp, n, s["P_shaft_kW"])
        sl = surge_line(dp, n)
        row["natural"] = None if op is None else dict(op, surge=sl)
        if op is not None:
            row["natural"]["confirm_gtu_cycle"] = confirm(dp, op)
        P_maxU = gg_point(dp, n, T41_MAX_U)
        row["at_T41_1000"] = dict(P_kW=P_maxU["P_kW"], W1=P_maxU["W1"], EGT_C=P_maxU["EGT_C"], P_share_pct=100 * P_maxU["P_kW"] / s["P_shaft_kW"])
        notes = []
        igv = None
        if op is not None and (op["W1"] > row["G_cap"] or not op["ok"]) and op["T41_C"] <= T41_MAX_U:
            target = min(row["G_cap"], (sl["W_choke"] or 1e9) - 0.5) - 0.01   # в допуске пакета и не у запирания
            row["igv_target_W1"] = target

            def wres(a):
                o = solve_power(dp, n, s["P_shaft_kW"], math.radians(a))
                return (o["W1"] - target) if o else -1e3
            try:
                a = brentq(wres, 0.0, 60.0, xtol=0.01)
                igv = solve_power(dp, n, s["P_shaft_kW"], math.radians(a))
                igv["confirm_gtu_cycle"] = confirm(dp, igv)
                igv["surge"] = surge_line(dp, n, math.radians(a))
            except ValueError:
                igv = None
                notes.append("поворотный входной аппарат до 60° не снижает расход до допуска")
        row["with_igv"] = igv
        use = igv or op
        if op is None:
            row["verdict"] = "нет: мощность не набирается и при 1250 °C"
        else:
            okT = use["T41_C"] <= T41_MAX_U + 1e-6
            okG = use["W1"] <= row["G_cap"] * (1 + 1e-4)
            okE = use["EGT_C"] <= row["EGT_cap"]
            okS = use["W1"] > (use.get("surge") or sl)["W_surge"]
            okC = bool(use["ok"])
            wch = (use.get("surge") or sl)["W_choke"]
            row["choke_margin_pct"] = None if wch is None else (wch / use["W1"] - 1) * 100
            row["surge_margin_pct"] = (1 - (use.get("surge") or sl)["W_surge"] / use["W1"]) * 100
            row["checks"] = dict(T41_ok=okT, G_ok=okG, EGT_ok=okE, surge_ok=okS, choke_ok=okC)
            if okT and okG and okE and okS and okC:
                tight = row["choke_margin_pct"] is not None and row["choke_margin_pct"] < 2.0
                row["verdict"] = ("да" + (" (с поворотным входным аппаратом)" if igv else "") +
                                  (f", но у запирания: запас по расходу {row['choke_margin_pct']:.1f} % — нужна характеристика (SU2)" if tight else ""))
            elif not okT:
                row["verdict"] = (f"нет в классе U: нужна T41 {use['T41_C']:.0f} °C; "
                                  f"при 1000 °C — {row['at_T41_1000']['P_share_pct']:.0f} % мощности")
            else:
                row["verdict"] = "нет: " + ", ".join(k for k, v in (("расход выше допуска", not okG), ("выхлоп выше допуска", not okE),
                                                                    ("за границей устойчивости", not okS),
                                                                    ("за границей запирания", not okC)) if v)
        row["notes"] = notes
        members.append(row)
    turb = {}
    for zz, o in dp["opts"].items():
        t = o["turb"]
        turb[str(zz)] = dict(t, w_t_kJkg=o["w_t"] / 1e3, P_shaft_kW=o["cyc"]["P_MW"] * 1e3, eff=o["cyc"]["eff"], EGT_C=o["cyc"]["EGT_C"])
    res = dict(
        title="Типоразмер № 3: согласованная расчётная точка газогенератора при 14 950 об/мин",
        method=("0D-цикл gtu_cycle.cycle (Cantera GRI-Mech 3.0, метан, ISO; одновальная схема — та же модель с одинаковым КПД "
                "турбины и без переходного канала) + одномерный расчёт средней линии двух центробежных ступеней и осевой "
                "турбины. КПД ступени компрессора = gtu_explore.comp_stage (принятая формула трека, коэффициент расхода 0,085) "
                "+ поправка: разность КПД по набору потерь Oh, Yoon, Chung (1997) с трением стенок диффузора (оценка) между "
                "колесом при 14 950 об/мин и колесом той же ступени при частоте, где коэффициент расхода 0,085; положительная "
                "поправка не засчитывается. Скольжение — Wiesner (1967), число лопаток — Пфлейдерер. КПД ступени турбины — "
                "gtu_explore.eta_axial(psi). Нерасчётные режимы — та же средняя линия на постоянной геометрии, горло "
                "соплового аппарата по закону конуса Стодолы, граница устойчивости — максимум степени сжатия на линии частоты."),
        software=(f"python {platform.python_version()}, cantera {ct.__version__}, numpy {np.__version__}, scipy {scipy.__version__}; "
                  "формулы потерь сверены с TurboFlow 0.1.18 (PyPI, MIT), файл centrifugal_compressor/loss_model_oh.py"),
        date=datetime.date.today().isoformat(),
        status=("оценка уровня облика: КПД узлов, пределы окружных скоростей и прочность — принятые величины трека; "
                "нужны расчёт течения (SU2), прочности (CalculiX) и эксперимент"),
        inputs=[
            dict(file="routes/gtu/gtu_family_results.json", field=f"slots[model={DESIGN[0]}, variant={DESIGN[1]}].{SCEN}",
                 used=dict(arch=PT["arch"], cool=PT["cool"], PR=PR, T41_C=T41_C, G=G, n_out=N, D2=PT["D2"], n_gg=PT["n_gg"],
                           U2=PT["U2"], eta=PT["eta"], w=PT["w"], wheel=PT["wheel"], rpm_k=PT["rpm_k"])),
            dict(file="routes/gtu/gtu_family_results.json", field=f"slots[model={DESIGN[0]}, variant={DESIGN[1]}]",
                 used=dict(P_kW=P_EL_KW, P_shaft_kW=P_SH_KW, G_air=SLOT["G_air"], EGT0=SLOT["EGT0"])),
            dict(file="routes/gtu/gtu_family_results.json", field=f"frames.{SCEN}", used=FRAME),
            dict(file="routes/gtu/gtu_family.py", field="SCEN, COOL['U'], ARCH['1В-2Ц'].eta_pc, ETA_GEN, _run (eta_pt_gg, eta_pt_pt, x_leak)",
                 used=dict(alpha=F.SCEN[SCEN]["alpha"], dT=F.SCEN[SCEN]["dT"], T41_max_U=T41_MAX_U, eta_pc_family=F.ARCH["1В-2Ц"]["eta_pc"],
                           eta_gen=ETA_GEN, eta_gearbox=ETA_GB, x_leak=X_LEAK)),
            dict(file="routes/gtu/gtu_explore.py", field="LIM['U_ti'], LIM['U_ax'], LIM['T_ax_max'], comp_stage, eta_axial",
                 used=dict(U_ti=U_TI, U_ax=U_AX, T_ax_max=T_AX)),
            dict(file="routes/gtu/platform/fea/materials.json", field="ti64.ys_MPa", used=MAT["ys_MPa"]),
            dict(file="routes/gtu/platform/fea/disk_results.json", field="cases[comp_T60_1C_ti, comp_T60_1C_ti_solid].result.vm_max_node_MPa, U_tip",
                 used={k: dict(vm=FEA[k]["vm_max_node_MPa"], U_tip=FEA[k]["U_tip"]) for k in ("comp_T60_1C_ti", "comp_T60_1C_ti_solid")}),
        ],
        assumptions=dict(A, note="все — оценка; в JSON — как в коде (A в frame3_point.py)"),
        problem=dict(D2_family=PT["D2"], n_family=PT["n_gg"], U2_family=PT["U2"], U2_at_n_out_same_D2=math.pi * PT["D2"] * N / 60),
        split=dict(pi1=pi1, pi2=PR / pi1, pi=PR, scan=scan, U2_design_limit=U_DES,
                   criterion="наименьшая работа сжатия (наибольший КПД всего компрессора) при U2 <= 0,9*620 = 558 м/с у обоих колёс",
                   without_margin=dict(pi1=pi1_620, pi2=PR / pi1_620, U2_limit=U_TI, scan=scan_620, eta_p_equiv_tables=c620["eta_eq"],
                                       w_kJkg=c620["w"] / 1e3, eta_p1=c620["s1"]["eta_p"], eta_p2=c620["s2"]["eta_p"],
                                       U2_1=c620["s1"]["act"]["U2"], U2_2=c620["s2"]["act"]["U2"], D2_1=c620["s1"]["act"]["D2"],
                                       D2_2=c620["s2"]["act"]["D2"], beta2_1=c620["s1"]["act"]["beta_deg"], beta2_2=c620["s2"]["act"]["beta_deg"])),
        stage1=s1r, stage2=s2r,
        compressor=dict(eta_p_equiv=dp["eta_pc"], eta_p_equiv_tables=dp["comp"]["eta_eq"], w_cantera_kJkg=sum(dp["w_stage_cantera"]) / 1e3,
                        w_tables_kJkg=dp["comp"]["w"] / 1e3, T3_C=dp["T_stage_cantera"][1] - 273.15),
        turbine=dict(selected_z=z, options=turb,
                     one_stage_feasible=turb["1"]["feasible"],
                     work_needed_kJkg=dp["opts"][z]["w_t"] / 1e3,
                     one_stage_max_work_kJkg=A["psi_max"] * U_AX ** 2 / 1e3,
                     one_stage_work_ratio=dp["opts"][z]["w_t"] / (A["psi_max"] * U_AX ** 2),
                     criteria=dict(psi_max=A["psi_max"], reaction_min=0.0, M2_max=A["m2_max"], U_max=U_AX, T41_max=T_AX)),
        cycle=dict(
            frame3=dict(P_shaft_kW=P_sh, P_el_kW=P_sh * ETA_GEN * ETA_GB, eff_shaft=cyc["eff"], eff_el=cyc["eff"] * ETA_GEN * ETA_GB,
                        w_kJkg=w_sp, EGT_C=cyc["EGT_C"], T3_C=cyc["T3_C"], wc_kJkg=cyc["wc_kJkg"], ER=cyc["GG_ER"] * cyc["PT_ER"],
                        W_exh=cyc["W_exh"], Wf=cyc["Wf"], eta_pc=dp["eta_pc"], eta_t=tb["eta_p"], dp_duct=0.0),
            frame3_with_duct_loss=dict(P_shaft_kW=duct["P_MW"] * 1e3, eff_shaft=duct["eff"], dp_duct=0.01),
            gtu_family_repro=dict(P_shaft_kW=fam["P_MW"] * 1e3, eff_shaft=fam["eff"], w_kJkg=fam["sp_work_kJkg"], EGT_C=fam["EGT_C"],
                                  eta_pc=F.ARCH["1В-2Ц"]["eta_pc"], eta_t_gg=0.865, eta_t_pt=0.885, dp_duct=0.01),
            gtu_family_json=dict(w_kJkg=PT["w"], eta_pct_after_rpm_penalty=PT["eta"], EGT_C=PT["EGT"]),
            pr_sensitivity_fixed_eta=pr_sens,
            required=dict(P_el_kW=P_EL_KW, P_shaft_kW=P_SH_KW, G_for_required_power=P_SH_KW / w_sp, G_cap=G_CAP,
                          T41_for_required_power_at_G=T41_need, EGT_cap=EGT_CAP),
            power_balance=dict(P_comp_cantera_kW=P_c, P_comp_stages_tables_kW=P_c_stages, P_turb_tables_kW=P_t_tab,
                               P_turb_needed_kW=P_t_need, P_shaft_cantera_kW=P_sh, eta_m=ETA_M,
                               residual_pct=(P_t_tab * ETA_M - P_c - P_sh) / P_t_tab * 100,
                               note="P_turb*eta_m = P_comp + P_shaft; турбина посчитана по таблицам (gtu_explore), компрессор и вал — Cantera"),
            turbine_crosscheck_tables=dict(ER_tables=hs_chk["ER"], ER_cantera=cyc["GG_ER"] * cyc["PT_ER"], EGT_tables_C=hs_chk["T5"] - 273.15,
                                           EGT_cantera_C=cyc["EGT_C"], ER_stagewise=tb["ER_total"],
                                           ER_target=(C.P0 * PR * (1 - DP_COMB)) / (C.P0 * P_EXH)),
        ),
        offdesign_selfcheck=dict(design_op, G_design=G),
        members=members,
    )
    return clean(res), dp


if __name__ == "__main__":
    import time
    t0 = time.time()
    res, dp = run()
    json.dump(res, open(OUT, "w"), ensure_ascii=False, indent=1)
    s1, s2, cy = res["stage1"], res["stage2"], res["cycle"]
    print(f"π1 {res['split']['pi1']} π2 {res['split']['pi2']:.3f}")
    for s in (s1, s2):
        print(f"{s['label']}: D1 {s['D1_hub_m']:.3f}/{s['D1_tip_m']:.3f} D2 {s['D2_m']:.3f} b2 {s['b2_m']*1e3:.1f} мм U2 {s['U2_m_s']:.0f} "
              f"φ {s['phi']:.4f} Z {s['blades_full']}+{s['blades_splitter']} β2 {s['beta2_blade_from_radial_deg']:.0f}° "
              f"η_p {s['eta_p']:.4f} (Δ {s['eta_basis']['d_eta_raw']:+.4f}) запас U2 {s['material']['margin_U_pct']:.1f} %")
    for z, t in res["turbine"]["options"].items():
        m2 = [st["M2"] for st in t["stages"]]
        print(f"турбина z={z}: Dm {t['Dm']:.3f} ψ {t['psi']:.2f} R {t['stages'][0]['reaction']:.2f} M2 {m2} "
              f"η {t['eta_p']:.3f} ER {t['ER_total']:.2f} годна {t['feasible']} P {t['P_shaft_kW']:.0f}")
    print("цикл:", {k: cy["frame3"][k] for k in ("P_shaft_kW", "P_el_kW", "eff_shaft", "eff_el", "EGT_C")})
    print("gtu_family:", cy["gtu_family_repro"], cy["gtu_family_json"])
    print("нужно:", cy["required"])
    for m in res["members"]:
        print(f"{m['model']:22s} {m['variant']:9s} G {m['G_soft']:.2f} ({m['share_of_design_pct']:.0f} %) {m['arch']}/{m['cool']} → {m['verdict']}")
    print(f"{time.time() - t0:.0f} с")
