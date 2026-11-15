"""Рекуператор KV-3 на решётке TPMS (гироид, Schwarz-D) с переменной плотностью по длине — против пластинчатого.

Нагрузка — варианты «рекуператор ε + межтурбинный подогрев» из organic/cycles (cycles_results.json, ключи
R80/R85/R90-ic-RH-ibc-none-none, поле opt): расходы, температуры входа и выхода обоих потоков, давления;
допустимые потери — как там: воздух 2,5 % p3, газ 3,5 % p_газа, на сердцевину — доля 0,7.

Модель сердцевины (уровень облика, оценка): противоток, N участков по длине; у каждого участка свой размер ячейки
a(x) и смещение уровня c(x) (лист c ± t, тракт A > c + t — газ, тракт B < c − t — воздух), толщина стенки t_w
постоянна (печать). Геометрия (ε, β, Dh, толщина стенки) — tpms_geom (marching cubes), переносы f·Re(ℓ) и Nu_H1(ℓ)
— расчёт tpms_transport (MAC, предел Стокса, проверка LBM); поправка на инерцию — по LBM, если расчёт есть.
Локально: h = Nu·k/Dh, dp/dx = 2·(f·Re)·μ·U/Dh², UA' = A_fr/(1/(h_A·β_A) + t_w/(k_w·β_mid) + 1/(h_B·β_B)).
Температуры — линейная система по участкам (средняя разность на участке), свойства по местной температуре
(Cantera: воздух и продукты сгорания метана при доле топлива из cycles), 4 итерации.
Пластинчатая сердцевина — та же модель с плоскими каналами (f·Re = 24, Nu_H1 = 8,235 — Shah & London;
и Nu_T = 7,54 как в cycles для сверки).
Оптимизация — наименьшая масса сердцевины при Q ≥ Q_цикла и потерях ≤ бюджета; переменные: фронт, длина,
a и c в трёх точках по длине (кусочно-линейно); ограничения печати — t_w, Dh ≥ Dh_min, a ≤ a_max (оценки).
"""
import json
import math
import sys
import time
from pathlib import Path

import numpy as np
from scipy.optimize import brentq, differential_evolution, minimize

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
GTU = ROOT / "routes" / "gtu"
sys.path.insert(0, str(HERE))

P0 = 101325.0
PRINT = dict(t_w=0.35e-3, t_w_thin=0.25e-3, Dh_min=1.5e-3, a_max=14e-3, a_min=2.0e-3, rho_m=8440.0, k_w=11.0,
             k_struct=0.6, box_m=(1.5, 1.5, 1.2), fill=0.4,
             note="стенка 0,35 мм (как в cycles) и вариант 0,25 мм; канал Dh ≥ 1,5 мм под вынос порошка; Inconel 625 "
                  "8440 кг/м³, k = 11 Вт/(м·К) при ~500 °C; коэффициент коллекторов и корпуса 0,6 и камера BLT-S1500 "
                  "1,5 × 1,5 × 1,2 м с заполнением 0,4 — как в organic/cycles (HX_TECH['lpbf']); всё — оценки")


# ------------------------------------------------------------------ свойства (Cantera)
class Props:
    def __init__(self, kind, p, far=0.0):
        import cantera as ct
        g = ct.Solution("gri30.yaml")
        g.transport_model = "mixture-averaged"
        if kind == "air":
            X = "N2:0.7808, O2:0.2095, AR:0.0093"
            g.TPX = 300, p, X
        else:
            # продукты полного сгорания метана при массовой доле топлива far (на кг воздуха)
            air = ct.Solution("gri30.yaml")
            air.TPX = 300, p, "N2:0.7808, O2:0.2095, AR:0.0093"
            Ya = air.Y
            iCH4, iO2, iCO2, iH2O = (g.species_index(s) for s in ("CH4", "O2", "CO2", "H2O"))
            Y = Ya / (1 + far)
            mf = far / (1 + far)
            n = mf / 16.043
            Y[iO2] -= 2 * n * 31.998
            Y[iCO2] += n * 44.009
            Y[iH2O] += 2 * n * 18.015
            g.TPY = 300, p, Y
            X = None
        self.Y = g.Y.copy()
        Ts = np.linspace(280, 1100, 83)
        tab = []
        for T in Ts:
            g.TPY = T, p, self.Y
            tab.append((g.cp_mass, g.viscosity, g.thermal_conductivity, g.density, g.enthalpy_mass))
        self.Ts, self.tab, self.p = Ts, np.array(tab), p

    def __call__(self, T):
        return [np.interp(T, self.Ts, self.tab[:, i]) for i in range(5)]

    def h(self, T):
        return np.interp(T, self.Ts, self.tab[:, 4])

    def T_of_h(self, h):
        return np.interp(h, self.tab[:, 4], self.Ts)


# ------------------------------------------------------------------ нагрузка из cycles
def duties():
    d = json.loads((GTU / "organic" / "cycles" / "cycles_results.json").read_text())
    A = d["assumptions"]
    out = {}
    for x in d["results"]:
        if x["key"] in ("R80-ic-RH-ibc-none-none", "R85-ic-RH-ibc-none-none", "R90-ic-RH-ibc-none-none"):
            o = x["opt"]
            p_h = A["p_exh"] * P0 * (1 + A["dp_rg"])
            hx = o["hx"]["lpbf"][0]
            out[x["key"]] = dict(eps=x["topo"]["rec"], PR=o["PR"], m_c=o["W_air"], m_h=o["W_exh"],
                                 Tc_in=o["T3_C"] + 273.15, Tc_out=o["T3r_C"] + 273.15, Th_in=o["T_rec_in_C"] + 273.15,
                                 Th_out_cycle=o["EGT_C"] + 273.15, p_c=o["p3_bar"] * 1e5, p_h=p_h,
                                 far=o["W_fuel_CH4"] / o["W_air"], Q_kW=o["Q_rec_kW"],
                                 dp_c=A["dp_ra"] * o["p3_bar"] * 1e5 * A["dp_core_share"],
                                 dp_h=A["dp_rg"] * p_h * A["dp_core_share"],
                                 cycles_lpbf=dict(V_core_m3=hx["V_core_m3"], mass_kg=hx["mass_kg"], L_m=hx["L_m"],
                                                  A_front_m2=hx["A_front_m2"], modules=hx["modules"]),
                                 cycles_foil=dict(V_core_m3=o["hx"]["foil"][0]["V_core_m3"],
                                                  mass_kg=o["hx"]["foil"][0]["mass_kg"]),
                                 eff_el=o["eff_el"])
    return out, dict(file="routes/gtu/organic/cycles/cycles_results.json",
                     field="results[key=R80|R85|R90-ic-RH-ibc-none-none].opt: W_air, W_exh, T3_C, T3r_C, T_rec_in_C, "
                           "EGT_C, p3_bar, W_fuel_CH4, Q_rec_kW, hx.lpbf[0], hx.foil[0]; assumptions: p_exh, dp_ra, "
                           "dp_rg, dp_core_share")


# ------------------------------------------------------------------ геометрия и переносы
class TPMSModel:
    """таблица по уровню ℓ: ε(ℓ), β(ℓ) — геометрия; f·Re(ℓ), Nu(ℓ) — из расчёта переносов"""

    def __init__(self, topo, transport_rows, inertia=None, n_geom=128):
        from tpms_geom import side_props
        self.topo = topo
        self.lv = np.round(np.arange(-1.0, 1.01, 0.1), 3)
        g = [side_props(topo, l, n_geom) for l in self.lv]
        self.eps = np.array([x["eps"] for x in g])
        self.beta = np.array([x["beta"] for x in g])
        rows = [r for r in transport_rows if r["topo"] == topo]
        L = sorted(set(r["level"] for r in rows))
        self.tl = np.array(L)
        self.fRe = np.array([np.mean([r["fRe"] for r in rows if r["level"] == l]) for l in L])
        self.Nu0 = np.array([np.mean([r["Nu"][0]["Nu"] for r in rows if r["level"] == l]) for l in L])
        self.inertia = inertia

    def E(self, l):
        return np.interp(l, self.lv, self.eps)

    def B(self, l):
        return np.interp(l, self.lv, self.beta)

    def t_level(self, c, tw_over_a):
        """полуширина уровня t для заданной толщины стенки t_w/a при смещении c"""
        bm = self.B(c)

        def r(t):
            return (1 - self.E(c + t) - self.E(t - c)) / bm - tw_over_a
        hi = 0.95 - abs(c)
        if r(hi) < 0:
            return None
        return brentq(r, 1e-6, hi)

    def side(self, l):
        """тракт F > l: ε, β (1/a), Dh (a), f·Re, Nu0"""
        e, b = self.E(l), self.B(l)
        return e, b, 4 * e / b, float(np.interp(l, self.tl, self.fRe)), float(np.interp(l, self.tl, self.Nu0))

    def corr(self, Re, Pr):
        """поправка на инерцию: (f·Re)/(f·Re)_0 и Nu/Nu_0 по LBM; нет данных — 1"""
        if self.inertia is None:
            return 1.0, 1.0
        I = self.inertia
        cf = float(np.interp(Re, I["Re"], I["f_ratio"]))
        cn = float(np.interp(Re * Pr / 0.7, I["Re"], I["Nu_ratio"]))
        return cf, cn


def core_eval(du, PH, PC, geom, prof, A_fr, L, N=40, tw=None, model=None, return_profiles=False):
    """geom: 'tpms' (prof: функции a(x), c(x)) или 'plate' (prof: dict gap, NuH); возвращает Q, потери, массу"""
    tw = PRINT["t_w"] if tw is None else tw
    x = (np.arange(N) + 0.5) / N * L
    dx = L / N
    mh, mc = du["m_h"], du["m_c"]
    Th = np.linspace(du["Th_in"], du["Th_out_cycle"], N + 1)
    Tc = np.linspace(du["Tc_out"], du["Tc_in"], N + 1)
    # геометрия участков
    if geom == "tpms":
        a = prof["a"](x / L)
        c = prof["c"](x / L)
        segs = []
        for ai, ci in zip(a, c):
            t = model.t_level(ci, tw / ai)
            if t is None:
                return None
            eA, bA, dA, fA, nA = model.side(ci + t)          # газ
            eB, bB, dB, fB, nB = model.side(t - ci)          # воздух (зеркало)
            bm = model.B(ci)
            ew = 1 - eA - eB
            segs.append((eA, bA / ai, dA * ai, fA, nA, eB, bB / ai, dB * ai, fB, nB, bm / ai, ew))
        S = np.array(segs)
    else:
        g = prof["gap"]
        pitch = 2 * g + 2 * tw
        e = g / pitch
        b = 2 / pitch
        Dh = 2 * g
        S = np.tile([e, b, Dh, 24.0, prof["Nu"], e, b, Dh, 24.0, prof["Nu"], b, 2 * tw / pitch], (N, 1))
    eA, bA, dA, fA, nA, eB, bB, dB, fB, nB, bm, ew = S.T
    for it in range(4):
        Thm = 0.5 * (Th[:-1] + Th[1:])
        Tcm = 0.5 * (Tc[:-1] + Tc[1:])
        cph, muh, kh, rhoh, _ = PH(Thm)
        cpc, muc, kc, rhoc, _ = PC(Tcm)
        Uh = mh / (rhoh * eA * A_fr)
        Uc = mc / (rhoc * eB * A_fr)
        Reh = rhoh * Uh * dA / muh
        Rec = rhoc * Uc * dB / muc
        Prh, Prc = cph * muh / kh, cpc * muc / kc
        if geom == "tpms" and model.inertia is not None:
            cfh, cnh = np.vectorize(model.corr)(Reh, Prh)
            cfc, cnc = np.vectorize(model.corr)(Rec, Prc)
        else:
            cfh = cnh = cfc = cnc = 1.0
        hh = nA * cnh * kh / dA
        hc = nB * cnc * kc / dB
        UAp = A_fr / (1 / (hh * bA) + tw / (PRINT["k_w"] * bm) + 1 / (hc * bB))
        UA = UAp * dx
        # линейная система: Th[0] = Th_in, Tc[N] = Tc_in; участок j: Q = UA·((Th_j+Th_j+1) − (Tc_j+Tc_j+1))/2
        n = N + 1
        M = np.zeros((2 * n, 2 * n))
        r = np.zeros(2 * n)
        M[0, 0] = 1; r[0] = du["Th_in"]
        M[1, n + N] = 1; r[1] = du["Tc_in"]
        Ch = mh * cph
        Cc = mc * cpc
        for j in range(N):
            k = 2 + 2 * j
            # горячий: Ch(Th_j − Th_j+1) − Q = 0
            M[k, j] += Ch[j] - UA[j] / 2; M[k, j + 1] += -Ch[j] - UA[j] / 2
            M[k, n + j] += UA[j] / 2; M[k, n + j + 1] += UA[j] / 2
            # холодный: Cc(Tc_j − Tc_j+1) − Q = 0
            M[k + 1, n + j] += Cc[j] + UA[j] / 2; M[k + 1, n + j + 1] += -Cc[j] + UA[j] / 2
            M[k + 1, j] += -UA[j] / 2; M[k + 1, j + 1] += -UA[j] / 2
        sol = np.linalg.solve(M, r)
        Th, Tc = sol[:n], sol[n:]
    Q = float(mc * (PC.h(Tc[0]) - PC.h(Tc[-1])))
    dph = float(np.sum(2 * fA * cfh * muh * Uh / dA ** 2 * dx))
    dpc = float(np.sum(2 * fB * cfc * muc * Uc / dB ** 2 * dx))
    V = A_fr * L
    m_core = PRINT["rho_m"] * float(np.mean(ew)) * V
    out = dict(Q_kW=Q / 1e3, dp_h=dph, dp_c=dpc, V=V, mass_core=m_core, mass_total=m_core * (1 + PRINT["k_struct"]),
               L=L, A_fr=A_fr, Re_h=[float(Reh.min()), float(Reh.max())], Re_c=[float(Rec.min()), float(Rec.max())],
               eps_air=float((Tc[0] - du["Tc_in"]) / (du["Th_in"] - du["Tc_in"])),
               Th_out=float(Th[-1]), Tc_out=float(Tc[0]),
               Dh_h_mm=[float(dA.min() * 1e3), float(dA.max() * 1e3)], Dh_c_mm=[float(dB.min() * 1e3), float(dB.max() * 1e3)],
               beta_mean=float(np.mean(0.5 * (bA + bB))), wall_frac=float(np.mean(ew)))
    bx = PRINT["box_m"]
    out["modules"] = int(max(1, math.ceil(V / (bx[0] * bx[1] * bx[2] * PRINT["fill"]))))
    # производство энтропии: теплообмен и трение
    Sgen_T = float(np.sum(UA * ((0.5 * (Th[:-1] + Th[1:]) - 0.5 * (Tc[:-1] + Tc[1:]))) ** 2
                          / (0.5 * (Th[:-1] + Th[1:]) * 0.5 * (Tc[:-1] + Tc[1:]))))
    Sgen_p = float(mh * dph / (np.mean(rhoh) * np.mean(Thm)) + mc * dpc / (np.mean(rhoc) * np.mean(Tcm)))
    out.update(S_gen_heat_W_K=Sgen_T, S_gen_friction_W_K=Sgen_p)
    if return_profiles:
        out["_prof"] = dict(x=x, Th=Th, Tc=Tc, Reh=Reh, Rec=Rec, dA=dA, dB=dB,
                            sgen_heat=UA * ((0.5 * (Th[:-1] + Th[1:]) - 0.5 * (Tc[:-1] + Tc[1:]))) ** 2
                            / (0.5 * (Th[:-1] + Th[1:]) * 0.5 * (Tc[:-1] + Tc[1:])) / dx)
        if geom == "tpms":
            out["_prof"].update(a=prof["a"](x / L), c=prof["c"](x / L))
    return out


def pw(v):
    """кусочно-линейный профиль по трём точкам (0, 0,5, 1)"""
    v = np.asarray(v, float)
    return lambda s: np.interp(s, [0.0, 0.5, 1.0], v)


def size_tpms(du, PH, PC, model, graded=True, tw=None, seed=0, maxiter=60):
    """наименьшая масса сердцевины при Q ≥ Q_цикла, dp ≤ бюджета, Dh ≥ Dh_min, a ≤ a_max"""
    tw = PRINT["t_w"] if tw is None else tw
    Qt = du["Q_kW"]

    def unpack(z):
        if graded:
            Afr, L, a0, a1, a2, c0, c1, c2 = z
            return Afr, L, pw([a0, a1, a2]), pw([c0, c1, c2])
        Afr, L, a0, c0 = z
        return Afr, L, pw([a0] * 3), pw([c0] * 3)

    def obj(z):
        Afr, L, fa, fc = unpack(z)
        r = core_eval(du, PH, PC, "tpms", dict(a=fa, c=fc), Afr, L, tw=tw, model=model)
        if r is None:
            return 1e9
        pen = 0.0
        pen += max(0.0, 1 - r["Q_kW"] / Qt) * 1e3
        pen += max(0.0, r["dp_h"] / du["dp_h"] - 1) * 1e2
        pen += max(0.0, r["dp_c"] / du["dp_c"] - 1) * 1e2
        pen += max(0.0, 1 - min(r["Dh_h_mm"][0], r["Dh_c_mm"][0]) / (PRINT["Dh_min"] * 1e3)) * 1e2
        return r["mass_core"] / 1e3 + pen * 10
    amin, amax = PRINT["a_min"], PRINT["a_max"]
    if graded:
        bounds = [(2.0, 40.0), (0.05, 1.5)] + [(amin, amax)] * 3 + [(-0.35, 0.3)] * 3
    else:
        bounds = [(2.0, 40.0), (0.05, 1.5), (amin, amax), (-0.35, 0.3)]
    res = differential_evolution(obj, bounds, seed=seed, maxiter=maxiter, popsize=14, tol=1e-6, polish=False,
                                 updating="deferred", workers=1)
    z = res.x
    loc = minimize(obj, z, method="Nelder-Mead", options=dict(maxiter=600, xatol=1e-5, fatol=1e-7))
    if loc.fun < res.fun:
        z = loc.x
    Afr, L, fa, fc = unpack(z)
    r = core_eval(du, PH, PC, "tpms", dict(a=fa, c=fc), Afr, L, tw=tw, model=model, return_profiles=True)
    r["z"] = [float(v) for v in z]
    r["objective"] = float(obj(z))
    r["feasible"] = bool(r["Q_kW"] >= 0.995 * Qt and r["dp_h"] <= 1.005 * du["dp_h"] and r["dp_c"] <= 1.005 * du["dp_c"])
    return r


def size_plate(du, PH, PC, gap, Nu, tw=None):
    """пластинчатая сердцевина: зазор gap (Dh = 2·gap), подбор фронта и длины на наименьшую массу"""
    tw = PRINT["t_w"] if tw is None else tw
    Qt = du["Q_kW"]

    def obj(z):
        Afr, L = z
        r = core_eval(du, PH, PC, "plate", dict(gap=gap, Nu=Nu), Afr, L, tw=tw)
        pen = max(0.0, 1 - r["Q_kW"] / Qt) * 1e3 + max(0.0, r["dp_h"] / du["dp_h"] - 1) * 1e2 + \
            max(0.0, r["dp_c"] / du["dp_c"] - 1) * 1e2
        return r["mass_core"] / 1e3 + pen * 10
    res = differential_evolution(obj, [(1.0, 40.0), (0.03, 1.5)], seed=0, maxiter=80, popsize=20, tol=1e-8, polish=False)
    loc = minimize(obj, res.x, method="Nelder-Mead", options=dict(maxiter=400, xatol=1e-6, fatol=1e-8))
    z = loc.x if loc.fun < res.fun else res.x
    r = core_eval(du, PH, PC, "plate", dict(gap=gap, Nu=Nu), z[0], z[1], tw=tw, return_profiles=True)
    r["z"] = [float(v) for v in z]
    r["feasible"] = bool(r["Q_kW"] >= 0.995 * Qt and r["dp_h"] <= 1.005 * du["dp_h"] and r["dp_c"] <= 1.005 * du["dp_c"])
    return r
