#!/usr/bin/env python3
"""Узел «граница расчёта, валы и муфта» трека gas-turbine (routes/gtu/shaft).

Что считает. По паспорту установки (fleet/gtu_fleet_<модель>.json), конфигурации TOML
(shaft_<конфиг>.toml) и модели цикла трека — таблицу на муфте выходного вала двигателя по режимам:
P_мех, кВт; M = P/omega, Н·м; n, об/мин; для генераторного исполнения ещё P_эл на клеммах.
Граница замены — муфта на валу двигателя: всё, что за ней (редуктор и генератор или
технологический компрессор), сохраняется.

Модули (подробно — README.md узла):
  1. цикл            cycle_front / cycle_ggt / cycle_pt — тот же расчёт, что gtu_cycle.cycle
                     (политропные шаги на Cantera, GRI-Mech 3.0, метан, охлаждение и утечки), но с
                     входом T1, p1 и давлением за турбиной p5 вместо ISO и с кэшем объектов Cantera;
                     совпадение с gtu_cycle.cycle при ISO проверяется (identity_check);
  2. карты            CompressorMap — обобщённая характеристика компрессора ГТД (карта NPSS axi-5
                     из pyCycle 4.4.0, масштаб на точку ISO; ДОПУЩЕНИЕ); пропускная способность
                     турбин — закон эллипса Стодолы; ProcessCompressor — техкомпрессор (EXAMPLE);
                     GenTrain — редуктор и генератор (EXAMPLE);
  3. баланс на валу   Unit.gs_* (одновальный: частота постоянна, рабочая точка на линии частоты из
                     пропускной способности турбины) и Unit.md_* (двухвальный: частота газогенератора
                     из согласования с двумя турбинами, частота СТ — из нагрузки техкомпрессора);
                     полная мощность — при постоянной T41 (предел; у MD ещё предел частоты
                     газогенератора), частичная — T41 по мощности; проверка вне ISO — validate()
                     против оцифрованных кривых листов данных Solar (digitize_solar_curves.py);
  4. ступень          stage_link — gtu_explore.comp_stage для колеса заменителя: что ступень
                     получает от вала (частота, приведённая частота по T_нв, мощность);
  5. выгрузка         shaft_<конфиг>.csv, shaft_results.json, data_needed.csv.

Запуск (из любой папки): python3 gtu_shaft.py [--out DIR] [конфиг.toml ...]
По умолчанию — оба конфига Taurus 60, файлы пишутся рядом со скриптом. Около 2–4 мин на 4 ядрах.
"""
import argparse, csv, datetime, hashlib, json, math, os, platform, sys, tomllib

import numpy as np
import cantera as ct
import scipy
from scipy.optimize import brentq, root
from scipy.interpolate import RegularGridInterpolator

HERE = os.path.dirname(os.path.abspath(__file__))
GTU = os.path.dirname(HERE)
if GTU not in sys.path:
    sys.path.insert(0, GTU)
_cwd = os.getcwd()
import gtu_cycle as C          # noqa: E402  (импорт меняет текущую папку — возвращаем)
import gtu_explore as E        # noqa: E402
os.chdir(_cwd)

T_ISO, P_ISO = 288.15, 101325.0
LHV = C.LHV
DATE = datetime.date.today().isoformat()
CONFIGS = ("shaft_taurus60_gs.toml", "shaft_taurus60_md.toml")


def val(x):
    """значение параметра TOML: число или таблица {value = ..., example = true, ...}"""
    return x["value"] if isinstance(x, dict) and "value" in x else x


def rpm2w(n):
    return 2 * math.pi * n / 60.0


# =========================================================================== 1. цикл
_SOL = {}


def _sol(key):
    s = _SOL.get(key)
    if s is None:
        s = _SOL[key] = ct.Solution("gri30.yaml")
    return s


def cycle_front(T1, P1, W1, PR, eta_pc, T41, e):
    """Компрессор и камера сгорания с подмешиванием охлаждения соплового аппарата.
    Повторяет начало gtu_cycle.cycle (функция after_mix — та же), но вход T1, P1 задаётся.
    T41 — К. Возврат: состояние перед рабочим колесом 1-й ступени (станция 41)."""
    a = _sol("air")
    a.TPX = T1, P1, C.AIR
    wc = C.poly_process(a, P1 * PR, eta_pc, True)
    T3, p3 = a.T, a.P
    Wca = W1 * (1 - e["x_ngv"] - e["x_rot"] - e["x_leak"])
    g, g2 = _sol("comb"), _sol("cool")
    gi, MW = g.species_index, g.molecular_weights

    def after_mix(far):
        f = ct.Quantity(g, constant="HP"); f.TPX = e["T_fuel"], p3, "CH4:1"; f.mass = far * Wca
        ai = ct.Quantity(g, constant="HP"); ai.TPX = T3, p3, C.AIR; ai.mass = Wca
        m = ai + f
        hm = (ai.enthalpy_mass * ai.mass + f.enthalpy_mass * f.mass) / (ai.mass + f.mass)
        Y = m.Y.copy()
        nf = Y[gi("CH4")] / MW[gi("CH4")]
        Y[gi("CH4")] = 0
        Y[gi("O2")] -= 2 * nf * MW[gi("O2")]
        Y[gi("CO2")] += nf * MW[gi("CO2")]
        Y[gi("H2O")] += 2 * nf * MW[gi("H2O")]
        p4 = p3 * (1 - e["dp_comb"])
        g.HPY = hm - (1 - e["eta_b"]) * far / (1 + far) * LHV, p4, Y
        T4 = g.T
        q4 = ct.Quantity(g, constant="HP"); q4.mass = Wca * (1 + far)
        c = ct.Quantity(g2, constant="HP"); c.TPX = T3, p4, C.AIR; c.mass = W1 * e["x_ngv"]
        return q4 + c, T4

    far = brentq(lambda f: after_mix(f)[0].T - T41, 0.002, 0.06, xtol=1e-7)
    mix, T4 = after_mix(far)
    W41, P41 = mix.mass, mix.P
    return dict(T1=T1, P1=P1, W1=W1, PR=PR, eta_pc=eta_pc, wc=wc, T3=T3, P3=p3, T4=T4, far=far,
                Wf=far * Wca, W41=W41, T41=mix.T, P41=P41, Y41=mix.Y.copy(),
                FF41=W41 * math.sqrt(mix.T) / P41)


def eta_nu(eta_d, ratio_N, dhs, dhs_d, k):
    """КПД турбины вне расчётной частоты (допущение): парабола по отношению скоростей nu = U/C0,
    eta = eta_d*(1 - k*(nu/nu_d - 1)^2), nu/nu_d = (N/N_d)*sqrt(dh_s,d/dh_s). Форма — как
    gtu_explore.turb_radial (там оптимум U/C0 = 0,70). ratio_N = None -> eta_d."""
    if ratio_N is None or dhs_d is None or k == 0:
        return eta_d
    r = ratio_N * math.sqrt(dhs_d / dhs)
    return eta_d * (1.0 - k * (r - 1.0) ** 2)


def _dhs(T, P, Y, p_out):
    g = _sol("tis")
    g.TPY = T, P, Y
    h = g.enthalpy_mass
    g.SP = g.entropy_mass, p_out
    return h - g.enthalpy_mass


def cycle_ggt(fr, e, ratio_N=None):
    """Турбина привода компрессора (у одновального — часть турбины до той же работы), подмешивание
    охлаждения ротора, переходный канал. Повторяет gtu_cycle.cycle; ratio_N — N/N_d для КПД по U/C0."""
    g = _sol("ggt")
    T41, P41, Y41, W41 = fr["T41"], fr["P41"], fr["Y41"], fr["W41"]
    need = fr["wc"] * fr["W1"] / e["eta_m"] / W41
    eta_used = [e["eta_gg"]]

    def resid(p_out):
        eta = e["eta_gg"]
        if ratio_N is not None:
            eta = eta_nu(eta, ratio_N, _dhs(T41, P41, Y41, p_out), e.get("dhs_gg_d"), e["k_nu"])
        eta_used[0] = eta
        g.TPY = T41, P41, Y41
        return C.poly_process(g, p_out, eta, False) - need

    p45 = brentq(resid, e["P5"] * 1.001, P41 * 0.98, xtol=1.0)
    resid(p45)
    T_out = g.T
    q = ct.Quantity(g, constant="HP"); q.mass = W41
    g2 = _sol("cool2")
    c = ct.Quantity(g2, constant="HP"); c.TPX = fr["T3"], p45, C.AIR; c.mass = fr["W1"] * e["x_rot"]
    m2 = q + c
    p45d = p45 * (1 - e["dp_duct"])
    T45, W45, Y45 = m2.T, m2.mass, m2.Y.copy()
    return dict(P45=p45, P45d=p45d, T45=T45, W45=W45, Y45=Y45, T_ggt_out=T_out, eta_gg=eta_used[0],
                dhs_gg=_dhs(T41, P41, Y41, p45), w_gg=need, FF45=W45 * math.sqrt(T45) / p45d)


def cycle_pt(gg, fr, e, ratio_N=None):
    """Силовая турбина (у одновального — остаток турбины) до p5; мощность на валу."""
    dhs = _dhs(gg["T45"], gg["P45d"], gg["Y45"], e["P5"])
    eta = eta_nu(e["eta_pt"], ratio_N, dhs, e.get("dhs_pt_d"), e["k_nu"])
    g = _sol("pt")
    g.TPY = gg["T45"], gg["P45d"], gg["Y45"]
    w = C.poly_process(g, e["P5"], eta, False)
    P = w * gg["W45"] * e["eta_m"]
    return dict(P_shaft=P, T5=g.T, eta_pt=eta, dhs_pt=dhs, w_pt=w, eff=P / (fr["Wf"] * LHV))


def balances(fr, gg, pt, e):
    """Балансы массы и энергии точки (абсолютные энтальпии Cantera, включая теплоту образования).
    Энергия: W1*h1 + Wf*hf = W_утечки*h3 + W45*h5 + P_вал + мех. потери + недожог."""
    g = _sol("bal")
    g.TPX = fr["T1"], fr["P1"], C.AIR; h1 = g.enthalpy_mass
    g.TPX = fr["T3"], fr["P3"], C.AIR; h3 = g.enthalpy_mass
    g.TPX = e["T_fuel"], fr["P3"], "CH4:1"; hf = g.enthalpy_mass
    g.TPY = pt["T5"], e["P5"], gg["Y45"]; h5 = g.enthalpy_mass
    W1, Wf = fr["W1"], fr["Wf"]
    Wleak = W1 * e["x_leak"]
    dm = W1 + Wf - Wleak - gg["W45"]
    loss_mech = fr["wc"] * W1 * (1 / e["eta_m"] - 1) + pt["w_pt"] * gg["W45"] * (1 - e["eta_m"])
    loss_comb = (1 - e["eta_b"]) * Wf * LHV
    dE = W1 * h1 + Wf * hf - (Wleak * h3 + gg["W45"] * h5 + pt["P_shaft"] + loss_mech + loss_comb)
    return dict(mass_resid_kg_s=dm, energy_resid_pct_fuel=100 * dE / (Wf * LHV))


# =========================================================================== 2. карты оборудования
class CompressorMap:
    """Обобщённая характеристика компрессора ГТД (ДОПУЩЕНИЕ): карта из maps/*.json, срез alpha,
    масштаб NPSS/pyCycle на расчётную точку: PR = 1 + (PR_map - 1)*s_PR, Wc = Wc_map*s_W,
    eta = eta_d*eta_map/eta_map_d (отношение карты применено к политропному КПД модели цикла),
    Nc_map = Nc_rel*Nc_map_d. Wc — приведённый расход W*sqrt(T1/288,15)/(p1/101 325)."""

    def __init__(self, path, Nc_map_d, R_d, PR_d, Wc_d, eta_d, alpha=0):
        m = json.load(open(path))
        self.path, self.name = path, m.get("name")
        self.Nc = np.array(m["NcMap"], float)
        self.R = np.array(m["RlineMap"], float)
        mk = lambda k: RegularGridInterpolator((self.Nc, self.R), np.array(m[k], float)[alpha],
                                               method="pchip", bounds_error=False, fill_value=None)
        self.fW, self.fPR, self.fE = mk("WcMap"), mk("PRmap"), mk("effMap")
        self.Nc_d, self.R_d, self.eta_d = Nc_map_d, R_d, eta_d
        x = (Nc_map_d, R_d)
        self.Wm_d, self.PRm_d, self.Em_d = float(self.fW(x)), float(self.fPR(x)), float(self.fE(x))
        self.sW, self.sPR = Wc_d / self.Wm_d, (PR_d - 1) / (self.PRm_d - 1)
        self.R_stall = float(m.get("RlineStall") or self.R[0])
        self.R_lo, self.R_hi = float(self.R[0]) - 0.4, float(self.R[-1]) + 0.4   # допустимая экстраполяция по R

    def __call__(self, Nc_rel, R):
        x = (Nc_rel * self.Nc_d, R)
        return (float(self.fW(x)) * self.sW, 1 + (float(self.fPR(x)) - 1) * self.sPR,
                self.eta_d * float(self.fE(x)) / self.Em_d)

    def flags(self, Nc_rel, R):
        f = []
        Ncm = Nc_rel * self.Nc_d
        if not (self.Nc[0] <= Ncm <= self.Nc[-1]):
            f.append("Nc вне карты")
        if R < self.R_stall:
            f.append("за границей помпажа карты")
        if R > self.R[-1]:
            f.append("за стороной запирания карты")
        return f

    def info(self):
        return dict(file=os.path.relpath(self.path, HERE), name=self.name, Nc_map_d=self.Nc_d, Rline_d=self.R_d,
                    map_at_design=dict(Wc=self.Wm_d, PR=self.PRm_d, eta_is=self.Em_d),
                    scale=dict(s_Wc=self.sW, s_PR=self.sPR, s_eta=self.eta_d / self.Em_d))


def ellipse_resid(FF, P_in, P_out, K):
    """закон эллипса Стодолы: W*sqrt(T_in)/p_in = K*sqrt(1 - (p_out/p_in)^2) -> относительная невязка"""
    return FF / (K * math.sqrt(max(1e-12, 1.0 - (P_out / P_in) ** 2))) - 1.0


class GenTrain:
    """Редуктор и генератор (EXAMPLE): потери = опорные * (доля_пост + (1 - доля_пост)*(P/P_опорн)^2).
    Опорная точка — паспорт ISO на клеммах. P_эл = P_мех*eta_ред*eta_ген тождественно."""

    def __init__(self, t, P_el_ref_W):
        self.eg, self.er = val(t["eta_gen_ref"]), val(t["eta_red_ref"])
        self.fg, self.fr = val(t["gen_loss_fixed_share"]), val(t["red_loss_fixed_share"])
        self.P_el_ref = P_el_ref_W
        self.P_gin_ref = P_el_ref_W / self.eg
        self.P_mech_ref = self.P_gin_ref / self.er
        self.Lg_ref = self.P_gin_ref - P_el_ref_W
        self.Lr_ref = self.P_mech_ref - self.P_gin_ref
        self.gen_rated = val(t["gen_rated_kWe"]) * 1e3
        self.red_rated = val(t["red_rated_kW"]) * 1e3
        self.grid_Hz, self.poles = val(t["grid_Hz"]), int(val(t["gen_poles"]))
        self.gen_rpm = 120.0 * self.grid_Hz / self.poles

    def el(self, P_mech):
        Lr = self.Lr_ref * (self.fr + (1 - self.fr) * (P_mech / self.P_mech_ref) ** 2)
        P_gin = P_mech - Lr
        Lg = self.Lg_ref * (self.fg + (1 - self.fg) * (P_gin / self.P_gin_ref) ** 2)
        P_el = P_gin - Lg
        return P_el, P_gin / P_mech, P_el / P_gin

    def mech_for_el(self, P_el):
        return brentq(lambda P: self.el(P)[0] - P_el, 0.3 * P_el, 2.0 * P_el, xtol=1e-3)

    def cap(self):
        """наибольшая P_мех, которую принимает агрегат: номинал редуктора и номинал генератора"""
        Pg = self.mech_for_el(self.gen_rated)
        return (Pg, "генератор") if Pg < self.red_rated else (self.red_rated, "редуктор")


class ProcessCompressor:
    """Технологический центробежный компрессор (EXAMPLE): линия расчётной частоты H/H_d = h(x),
    eta/eta_d = 1 - k*(x - 1)^2, x = Q/(s*Q_d), s = n/n_d; другие частоты — законы подобия
    (Q ~ s, H ~ s^2, eta по x). Рабочая точка — пересечение с сетью H = H_d*(hs + (1 - hs)*(Q/Q_d)^2).
    Q_d подобран так, что при ISO и полной мощности компрессор принимает всю мощность СТ на n_d."""
    R_U = 8314.462618

    def __init__(self, t, P_design_W):
        self.n_d = val(t["n_pt_d_rpm"])
        self.s_lo, self.s_hi = val(t["speed_range"])
        self.H_d = val(t["H_p_d_kJkg"]) * 1e3
        self.eta_d = val(t["eta_p_d"])
        self.c = val(t["head_coeffs"])
        self.k = val(t["eta_k"])
        self.x_surge, self.x_choke = val(t["x_surge"]), val(t["x_choke"])
        self.P_ml = val(t["P_mloss_kW"]) * 1e3
        self.MW, self.Z, self.kappa = val(t["gas_MW"]), val(t["Z_s"]), val(t["kappa"])
        self.p_s, self.T_s = val(t["p_s_Pa"]), val(t["T_s_K"])
        self.hs = val(t["H_static_share"])
        self.rho = self.p_s * self.MW / (self.Z * self.R_U * self.T_s)
        self.m_d = (P_design_W - self.P_ml) * self.eta_d / self.H_d
        self.Q_d = self.m_d / self.rho
        self.P_design = P_design_W

    def h(self, x):
        return self.c[0] + self.c[1] * x + self.c[2] * x * x

    def point(self, n):
        s = n / self.n_d
        Hs = lambda Q: self.H_d * (self.hs + (1 - self.hs) * (Q / self.Q_d) ** 2)
        g = lambda Q: s * s * self.H_d * self.h(Q / (s * self.Q_d)) - Hs(Q)
        x_pk = -self.c[1] / (2 * self.c[2])            # вершина h(x): левее — неустойчивая ветвь
        lo, hi = x_pk * s * self.Q_d, 2.0 * s * self.Q_d
        if g(lo) < 0:
            return None                                # напор компрессора ниже статического напора сети
        Q = brentq(g, lo, hi, xtol=1e-9)
        x = Q / (s * self.Q_d)
        H = Hs(Q)
        eta = self.eta_d * (1 - self.k * (x - 1) ** 2)
        m = self.rho * Q
        P = m * H / eta + self.P_ml
        npn = (self.kappa - 1) / (self.kappa * eta)    # (n-1)/n политропы
        PR = (1 + H * self.MW / (self.Z * self.R_U * self.T_s) * npn) ** (1 / npn)
        fl = []
        if x < self.x_surge:
            fl.append("техкомпрессор левее границы помпажа")
        if x > self.x_choke:
            fl.append("техкомпрессор правее границы запирания")
        return dict(n=n, s=s, Q_m3s=Q, x=x, H_kJkg=H / 1e3, eta_p=eta, m_kg_s=m, P_W=P, PR=PR,
                    p_d_Pa=self.p_s * PR, surge_margin=x / self.x_surge - 1, flags=fl)

    def n_for_power(self, P, n_lo, n_hi):
        f = lambda n: (self.point(n) or {"P_W": -1e30})["P_W"] - P
        if f(n_hi) < 0:
            return None
        if f(n_lo) > 0:
            return n_lo
        return brentq(f, n_lo, n_hi, xtol=1e-4)

    def info(self):
        return dict(n_d_rpm=self.n_d, Q_d_m3s=self.Q_d, m_d_kg_s=self.m_d, rho_s_kg_m3=self.rho,
                    H_d_kJkg=self.H_d / 1e3, eta_d=self.eta_d, P_design_kW=self.P_design / 1e3,
                    sized="Q_d подобран: при ISO и полной мощности компрессор принимает мощность СТ на n_d")


# =========================================================================== 3. баланс мощности на валу
class Unit:
    def __init__(self, cfg_path, overrides=None):
        self.cfg_path = os.path.abspath(cfg_path)
        cfg = tomllib.load(open(self.cfg_path, "rb"))
        for sec, kv in (overrides or {}).items():
            cfg.setdefault(sec, {}).update(kv)
        self.cfg = cfg
        u = cfg["unit"]
        self.key, self.layout = u["key"], u["layout"]
        self.gen = self.layout == "single_shaft_gen"
        self.pp_path = os.path.normpath(os.path.join(HERE, u["passport"]))
        recs = json.load(open(self.pp_path))["records"]
        self.rec = next(r for r in recs if r["variant"] == u["variant"])
        self.pp = {}
        for k, path in cfg["passport_fields"].items():
            v = self.rec
            for part in path.split("."):
                v = v[part]
            self.pp[k] = float(v)
        s, en = cfg["site"], cfg["engine"]
        self.P_amb = s["P_amb_Pa"]
        self.P1 = self.P_amb - s["dp_inlet_Pa"]
        self.P5 = (self.P_amb + s["dp_exhaust_Pa"]) * en["p_exh_ratio"]
        self.e0 = dict(x_ngv=en["x_ngv"], x_rot=en["x_rot"], x_leak=en["x_leak"], dp_comb=en["dp_comb"],
                       dp_duct=en["dp_duct"], eta_b=en["eta_b"], eta_m=en["eta_m"], T_fuel=s["T_fuel_K"],
                       k_nu=en["k_nu"])
        self.calibrate()
        cm = cfg["compressor_map"]
        self.map = CompressorMap(os.path.join(HERE, cm["file"]), cm["Nc_map_d"], cm["Rline_d"],
                                 self.d["PR"], self.d["W1"], self.d["eta_pc"], cm.get("alpha_slice", 0))
        if self.gen:
            self.N = self.pp["N_out_rpm"]
            self.train = GenTrain(cfg["train"], self.pp["P_out_kW"] * 1e3)
        else:
            self.N_gg_d = val(en["N_gg_d_rpm"])
            self.gg_limit = bool(en.get("gg_speed_limit", False))
            self.proc = ProcessCompressor(cfg["train"], self.d["P"])
            self.n_hi = min(self.pp["N_pt_max_rpm"], self.proc.s_hi * self.proc.n_d)
            self.n_lo = self.proc.s_lo * self.proc.n_d
        self._x0 = [1.0, cm["Rline_d"]]

    # ------------------------------------------------------------ калибровка (как gtu_study.calib)
    def e_at(self, eta3, P5):
        return dict(self.e0, eta_pc=eta3[0], eta_gg=eta3[1], eta_pt=eta3[2], P5=P5)

    def calibrate(self):
        en = self.cfg["engine"]
        base = en["eta_base"]
        if self.gen:
            t = self.cfg["train"]
            eg, er = val(t["eta_gen_ref"]), val(t["eta_red_ref"])
            P_t = self.pp["P_out_kW"] * 1e3 / (eg * er)
            eff_t = self.pp["eff_pct"] / 100 / (eg * er)
        else:
            P_t, eff_t = self.pp["P_out_kW"] * 1e3, self.pp["eff_pct"] / 100
        W1 = self.pp["exh_flow_kg_s"] / en["W1_exh_div"]
        PR = self.pp["PR"]
        P5 = P_ISO * en["p_exh_ratio"]          # калибровка — на условиях паспорта ISO, без потерь
        Tlo, Thi = (x + 273.15 for x in en["T41_bracket_C"])

        def run(d, T41):
            e = self.e_at([b + d for b in base], P5)
            fr = cycle_front(T_ISO, P_ISO, W1, PR, e["eta_pc"], T41, e)
            gg = cycle_ggt(fr, e)
            pt = cycle_pt(gg, fr, e)
            return fr, gg, pt, e

        fitT = lambda d: brentq(lambda T: run(d, T)[2]["P_shaft"] - P_t, Tlo, Thi, xtol=1e-3)
        d = brentq(lambda d: run(d, fitT(d))[2]["eff"] - eff_t, -0.05, 0.05, xtol=1e-7)
        T41 = fitT(d)
        fr, gg, pt, e = run(d, T41)
        self.eta3 = [b + d for b in base]
        K41 = fr["FF41"] / math.sqrt(1 - ((P5 if self.gen else gg["P45"]) / fr["P41"]) ** 2)
        K45 = gg["FF45"] / math.sqrt(1 - (P5 / gg["P45d"]) ** 2)
        self.e0.update(dhs_gg_d=gg["dhs_gg"], dhs_pt_d=pt["dhs_pt"])
        self.K41, self.K45 = K41, K45
        self.d = dict(d=d, T41=T41, W1=W1, PR=PR, eta_pc=e["eta_pc"], eta_gg=e["eta_gg"], eta_pt=e["eta_pt"],
                      P=pt["P_shaft"], eff=pt["eff"], T5=pt["T5"], T45=gg["T45"], W45=gg["W45"], Wf=fr["Wf"],
                      T3=fr["T3"], P_target=P_t, eff_target=eff_t, K41=K41, K45=K45,
                      dhs_gg=gg["dhs_gg"], dhs_pt=pt["dhs_pt"])

    # ------------------------------------------------------------ общие куски
    def e_site(self):
        return self.e_at(self.eta3, self.P5)

    def front(self, T1, Nc, R, T41, e):
        if not (self.map.R_lo <= R <= self.map.R_hi):           # дальняя экстраполяция даёт ложные ветви решений
            raise ValueError(f"R вне карты: {R:.3f}")
        Wc, PR, eta = self.map(Nc, R)
        if not (PR > 1.5 and Wc > 0 and 0.3 < eta < 1.0):     # экстраполяция карты ушла за смысл
            raise ValueError(f"карта вне области: Nc={Nc:.3f} R={R:.3f} PR={PR:.3f}")
        W1 = Wc * (self.P1 / P_ISO) / math.sqrt(T1 / T_ISO)
        fr = cycle_front(T1, self.P1, W1, PR, eta, T41, e)
        fr.update(Nc=Nc, R=R)
        return fr

    def _row(self, T1, load, n, P_mech, fr, gg, pt, e, limit, P_engine, extra=None):
        w = rpm2w(n)
        row = dict(regime=f"T{T1 - 273.15:+.0f}_L{load:.0f}", T_amb_C=T1 - 273.15, load_pct=load, n_rpm=n,
                   P_mech_kW=P_mech / 1e3, M_Nm=P_mech / w, P_el_kW=None,
                   W_air_kg_s=fr["W1"], W_fuel_kg_s=fr["Wf"], T_exh_C=pt["T5"] - 273.15,
                   eta_mech_pct=100 * P_mech / (fr["Wf"] * LHV), eta_el_pct=None,
                   P_engine_max_kW=P_engine / 1e3, limit=limit, T41_C=fr["T41"] - 273.15,
                   T45_C=gg["T45"] - 273.15, PR=fr["PR"], Nc_rel=fr["Nc"], Rline=fr["R"],
                   W_exh_kg_s=gg["W45"], eta_pc=fr["eta_pc"], eta_gg=gg["eta_gg"], eta_pt=pt["eta_pt"])
        row.update(balances(fr, gg, pt, e))
        fl = self.map.flags(fr["Nc"], fr["R"])
        row.update(extra or {})
        row["flags"] = "; ".join(fl + row.pop("_flags", []))
        return row

    # ------------------------------------------------------------ одновальный генераторный
    def gs_engine(self, T1, T41, e=None):
        """частота постоянна (N = N_d): приведённая Nc = sqrt(288,15/T1); R — из пропускной способности турбины"""
        e = e or self.e_site()
        Nc = math.sqrt(T_ISO / T1)

        def res(R):
            try:
                fr = self.front(T1, Nc, R, T41, e)
            except ValueError:            # вне карты: левее — расход мал (помпаж), правее — велик (запирание)
                return -1.0 if R < self.map.R_d else 1.0
            return ellipse_resid(fr["FF41"], fr["P41"], e["P5"], self.K41)

        R = _bracket_root(res, self.map.R_lo, self.map.R_hi)
        fr = self.front(T1, Nc, R, T41, e)
        gg = cycle_ggt(fr, e, 1.0)
        pt = cycle_pt(gg, fr, e, 1.0)
        return fr, gg, pt, e

    def gs_T41_for(self, T1, P_target, T41_hi):
        f = lambda T: self.gs_engine(T1, T)[2]["P_shaft"] - P_target
        lo = T41_hi - 200.0
        while f(lo) > 0:
            lo -= 100.0
        return brentq(f, lo, T41_hi, xtol=1e-3)

    def gs_full(self, T1, control="T41"):
        fr, gg, pt, e = self.gs_engine(T1, self.d["T41"])
        P_eng = pt["P_shaft"]
        cap, who = self.train.cap()
        if P_eng > cap:
            T41 = self.gs_T41_for(T1, cap, self.d["T41"])
            fr, gg, pt, e = self.gs_engine(T1, T41)
            return fr, gg, pt, e, who, P_eng
        return fr, gg, pt, e, "T41", P_eng

    def gs_rows(self):
        rows = []
        for Tc in self.cfg["regimes"]["T_amb_C"]:
            T1 = Tc + 273.15
            full = self.gs_full(T1)
            P_full, P_eng, lim_full = full[2]["P_shaft"], full[5], full[4]
            for load in self.cfg["regimes"]["load_pct"]:
                if load >= 100.0:
                    fr, gg, pt, e, lim = full[0], full[1], full[2], full[3], lim_full
                else:
                    T41 = self.gs_T41_for(T1, P_full * load / 100, full[0]["T41"])
                    fr, gg, pt, e = self.gs_engine(T1, T41)
                    lim = "нагрузка"
                P = pt["P_shaft"]
                P_el, er, eg = self.train.el(P)
                rows.append(self._row(T1, load, self.N, P, fr, gg, pt, e, lim, P_eng,
                                      dict(P_el_kW=P_el / 1e3, eta_el_pct=100 * P_el / (fr["Wf"] * LHV),
                                           eta_red=er, eta_gen=eg, n_gen_rpm=self.train.gen_rpm)))
        return rows

    # ------------------------------------------------------------ двухвальный механический привод
    def md_gg(self, T1, T41, e=None, x0=None, _nested=False):
        """газогенератор: неизвестные Nc (приведённая частота) и R; уравнения — пропускная способность
        турбины газогенератора и СТ (эллипс Стодолы), работа турбины = работа компрессора / eta_m"""
        e = e or self.e_site()
        th = math.sqrt(T1 / T_ISO)

        def F(x):
            try:
                fr = self.front(T1, x[0], x[1], T41, e)
                gg = cycle_ggt(fr, e, x[0] * th)
            except ValueError:
                return [1.0 + abs(x[0]), 1.0 + abs(x[1])]
            return [ellipse_resid(fr["FF41"], fr["P41"], gg["P45"], self.K41),
                    ellipse_resid(gg["FF45"], gg["P45d"], e["P5"], self.K45)]

        ok = lambda s: s.success and max(abs(v) for v in s.fun) < 1e-7
        sol = None
        for start in (x0, self._x0, [1.0, self.map.R_d]):
            if start is None:
                continue
            sol = root(F, start, method="hybr", options=dict(xtol=1e-11, maxfev=400))
            if ok(sol):
                break
        if not ok(sol) and not _nested:       # продолжение по (T1, T41) от точки ISO, где решение известно
            x = [1.0, self.map.R_d]
            for k in range(1, 13):
                f = k / 12
                Tk, T41k = T_ISO + f * (T1 - T_ISO), self.d["T41"] + f * (T41 - self.d["T41"])
                x = list(self.md_gg(Tk, T41k, e, x, _nested=True)[0]["NcR"])
            sol = root(F, x, method="hybr", options=dict(xtol=1e-11, maxfev=400))
        if not ok(sol):
            raise RuntimeError(f"газогенератор не согласован: T1={T1:.2f} T41={T41:.2f} {sol.message}")
        self._x0 = [float(sol.x[0]), float(sol.x[1])]
        fr = self.front(T1, sol.x[0], sol.x[1], T41, e)
        fr["NcR"] = (float(sol.x[0]), float(sol.x[1]))
        gg = cycle_ggt(fr, e, sol.x[0] * th)
        return fr, gg, e

    def n_gg(self, fr, T1):
        """физическая частота газогенератора, об/мин: Nc_rel * sqrt(T1/288,15) * N_gg_d"""
        return fr["Nc"] * math.sqrt(T1 / T_ISO) * self.N_gg_d

    def md_gg_at_speed(self, T1, Nc, e=None, x0=None):
        """газогенератор при заданной приведённой частоте: неизвестные R и T41 (те же два уравнения)"""
        e = e or self.e_site()
        th = math.sqrt(T1 / T_ISO)

        def F(x):
            try:
                fr = self.front(T1, Nc, x[0], x[1] * 1000.0, e)
                gg = cycle_ggt(fr, e, Nc * th)
            except ValueError:
                return [1.0 + abs(x[0]), 1.0 + abs(x[1])]
            return [ellipse_resid(fr["FF41"], fr["P41"], gg["P45"], self.K41),
                    ellipse_resid(gg["FF45"], gg["P45d"], e["P5"], self.K45)]

        start = x0 or [self.map.R_d, self.d["T41"] / 1000.0 - 0.05]
        sol = root(F, start, method="hybr", options=dict(xtol=1e-11, maxfev=400))
        if not (sol.success and max(abs(v) for v in sol.fun) < 1e-7):
            raise RuntimeError(f"газогенератор на пределе частоты не согласован: T1={T1:.2f} {sol.message}")
        T41 = float(sol.x[1]) * 1000.0
        return self.md_gg(T1, T41, e, [Nc, float(sol.x[0])])

    def md_limited_gg(self, T1, T41, e=None):
        """газогенератор на пределах: T41 не выше заданной и (если включено) частота не выше паспортной.
        Если предел частоты действует, решается прямо: Nc = (N_max/N_gg_d)/sqrt(T1/288,15), неизвестные R и T41."""
        Nmax = self.pp["N_gg_max_rpm"]
        try:
            fr, gg, ee = self.md_gg(T1, T41, e)
            over = self.n_gg(fr, T1) > Nmax * (1 + 1e-9)
        except RuntimeError:
            if not self.gg_limit:
                raise
            over = True                         # без предела точка уходит за карту — предел частоты заведомо действует
        if self.gg_limit and over:
            Nc = Nmax / self.N_gg_d / math.sqrt(T1 / T_ISO)
            fr, gg, ee = self.md_gg_at_speed(T1, Nc, e)
            if fr["T41"] <= T41 + 1e-6:
                return fr, gg, ee, "частота газогенератора"
            fr, gg, ee = self.md_gg(T1, T41, e)
        return fr, gg, ee, None

    def md_engine_opt(self, T1, gg_limit=None):
        """полная мощность при оптимальной частоте СТ (основа графика изготовителя): КПД СТ = расчётный"""
        keep = self.gg_limit
        if gg_limit is not None:
            self.gg_limit = gg_limit
        try:
            fr, gg, e, lim = self.md_limited_gg(T1, self.d["T41"])
        finally:
            self.gg_limit = keep
        pt = cycle_pt(gg, fr, e, None)
        return dict(P_kW=pt["P_shaft"] / 1e3, T41_C=fr["T41"] - 273.15, n_gg_rpm=self.n_gg(fr, T1), Nc=fr["Nc"],
                    R=fr["R"], limit=lim or "T41", T45_C=gg["T45"] - 273.15, T_exh_C=pt["T5"] - 273.15)

    def md_pt(self, fr, gg, e, n):
        return cycle_pt(gg, fr, e, n / self.proc.n_d)

    def md_T41_for(self, T1, n, P_target, T41_hi, e=None):
        def f(T):
            fr, gg, ee = self.md_gg(T1, T, e)
            return self.md_pt(fr, gg, ee, n)["P_shaft"] - P_target
        lo = T41_hi - 150.0
        while f(lo) > 0:
            lo -= 100.0
        return brentq(f, lo, T41_hi, xtol=1e-3)

    def md_full(self, T1, control="T41", e=None):
        """полная мощность: T41 = предел (или T45 = расчётной при control='T45'); частота СТ — где мощность СТ
        равна мощности техкомпрессора на сети; выше верхней частоты — ограничение по частоте"""
        if control == "T45":
            try:
                T41 = brentq(lambda T: self.md_gg(T1, T, e)[1]["T45"] - self.d["T45"],
                             self.d["T41"] - 150, self.d["T41"] + 60, xtol=1e-3)
            except (RuntimeError, ValueError):
                # без предела частоты точка уходит за карту; если предел частоты действует раньше, чем T45
                # дойдёт до уставки, полная мощность задана пределом частоты
                if not self.gg_limit:
                    raise
                fr, gg, ee, lim = self.md_limited_gg(T1, self.d["T41"] + 60, e)
                if lim is None or gg["T45"] > self.d["T45"]:
                    raise
                T41 = fr["T41"]
        else:
            T41 = self.d["T41"]
        fr, gg, ee, lim_gg = self.md_limited_gg(T1, T41, e)
        T41 = fr["T41"]
        lim = lim_gg or ("T45" if control == "T45" else "T41")
        f = lambda n: self.md_pt(fr, gg, ee, n)["P_shaft"] - self.proc.point(n)["P_W"]
        if f(self.n_hi) > 0:                                   # СТ даёт больше, чем компрессор берёт на верхней частоте
            n = self.n_hi
            P_eng = self.md_pt(fr, gg, ee, n)["P_shaft"]
            P = self.proc.point(n)["P_W"]
            T41b = self.md_T41_for(T1, n, P, T41, e)
            fr, gg, ee = self.md_gg(T1, T41b, e)
            return fr, gg, ee, n, "частота СТ (верх диапазона)", P_eng
        if f(self.n_lo) < 0:
            raise RuntimeError(f"T1={T1}: мощности не хватает даже на нижней частоте")
        n = brentq(f, self.n_lo, self.n_hi, xtol=1e-4)
        return fr, gg, ee, n, lim, self.md_pt(fr, gg, ee, n)["P_shaft"]

    def md_rows(self):
        rows = []
        for Tc in self.cfg["regimes"]["T_amb_C"]:
            T1 = Tc + 273.15
            full = self.md_full(T1)
            fr0, gg0, e0, n0, lim0, P_eng = full
            P_full = self.md_pt(fr0, gg0, e0, n0)["P_shaft"]
            for load in self.cfg["regimes"]["load_pct"]:
                if load >= 100.0:
                    fr, gg, e, n, lim = fr0, gg0, e0, n0, lim0
                else:
                    P_t = P_full * load / 100
                    n = self.proc.n_for_power(P_t, self.n_lo, self.n_hi)
                    lim = "нагрузка"
                    if n is None:
                        raise RuntimeError("нагрузка выше возможностей компрессора")
                    if self.proc.point(n)["P_W"] > P_t * (1 + 1e-6):
                        lim = "нижняя частота техкомпрессора"
                        P_t = self.proc.point(n)["P_W"]
                    T41 = self.md_T41_for(T1, n, P_t, fr0["T41"])
                    fr, gg, e = self.md_gg(T1, T41)
                pt = self.md_pt(fr, gg, e, n)
                cp = self.proc.point(n)
                n_gg = self.n_gg(fr, T1)
                fl = list(cp["flags"])
                if n_gg > self.pp["N_gg_max_rpm"] * (1 + 1e-6):
                    fl.append("частота газогенератора выше паспортного предела")
                rows.append(self._row(T1, load, n, pt["P_shaft"], fr, gg, pt, e, lim, P_eng,
                                      dict(n_gg_rpm=n_gg, Q_m3s=cp["Q_m3s"], H_kJkg=cp["H_kJkg"],
                                           eta_p_proc=cp["eta_p"], m_gas_kg_s=cp["m_kg_s"], PR_proc=cp["PR"],
                                           P_proc_kW=cp["P_W"] / 1e3, surge_margin=cp["surge_margin"],
                                           _flags=fl)))
        return rows

    def rows(self):
        return self.gs_rows() if self.gen else self.md_rows()

    # ------------------------------------------------------------ полная мощность для оценок чувствительности
    def full_engine(self, T1, control="T41"):
        if self.gen:
            fr, gg, pt, e = self.gs_engine(T1, self.d["T41"])
            if control == "T45":
                T41 = brentq(lambda T: self.gs_engine(T1, T)[1]["T45"] - self.d["T45"],
                             self.d["T41"] - 150, self.d["T41"] + 60, xtol=1e-3)
                fr, gg, pt, e = self.gs_engine(T1, T41)
            return dict(P_engine_kW=pt["P_shaft"] / 1e3, n_rpm=self.N, T41_C=fr["T41"] - 273.15, R=fr["R"])
        fr, gg, e, n, lim, P_eng = self.md_full(T1, control)
        return dict(P_engine_kW=P_eng / 1e3, n_rpm=n, T41_C=fr["T41"] - 273.15, R=fr["R"], limit=lim)

    def info(self):
        d = dict(self.d)
        d.update(T41_C=d.pop("T41") - 273.15, T5_C=d.pop("T5") - 273.15, T45_C=d.pop("T45") - 273.15,
                 T3_C=d.pop("T3") - 273.15, P_kW=d.pop("P") / 1e3, eff_pct=100 * d.pop("eff"),
                 P_target_kW=d.pop("P_target") / 1e3, eff_target_pct=100 * d.pop("eff_target"))
        out = dict(calibration=d, map=self.map.info(),
                   passport=dict(file=os.path.relpath(self.pp_path, HERE), variant=self.rec["variant"],
                                 values=self.pp, fields=self.cfg["passport_fields"],
                                 sources={k: self.rec["sources"].get(f.split(".")[-1])
                                          for k, f in self.cfg["passport_fields"].items()}))
        if self.gen:
            tr = self.train
            out["train"] = dict(P_mech_ref_kW=tr.P_mech_ref / 1e3, P_el_ref_kW=tr.P_el_ref / 1e3,
                                eta_gen_ref=tr.eg, eta_red_ref=tr.er, gen_rpm=tr.gen_rpm, grid_Hz=tr.grid_Hz,
                                cap_kW=tr.cap()[0] / 1e3, cap_by=tr.cap()[1],
                                red_rated_torque_Nm=tr.red_rated / rpm2w(self.N), example=True)
        else:
            out["process_compressor"] = dict(self.proc.info(), example=True)
            out["speed_window_rpm"] = [self.n_lo, self.n_hi]
        return out


def _bracket_root(f, lo, hi):
    a, b = f(lo), f(hi)
    if a * b > 0:
        raise RuntimeError("рабочая точка на линии частоты не найдена в пределах карты (±0,4 по R)")
    return brentq(f, lo, hi, xtol=1e-9)


# =========================================================================== 4. ступень (начало meanline)
def stage_link(units):
    """Что вал передаёт модели ступени заменителя. Расчётная точка заменителя — gtu_family_results.json
    (сценарий «мягкий», решение владельца 14.11.2026); колесо — gtu_explore.comp_stage (1D)."""
    fam_path = os.path.join(GTU, "gtu_family_results.json")
    fam = json.load(open(fam_path))
    out = dict(method="gtu_explore.comp_stage (одномерный размер центробежного колеса) на расчётной точке "
                      "заменителя из gtu_family_results.json, сценарий «мягкий»; частоты и мощности — из таблиц вала",
               inputs=[os.path.relpath(fam_path, HERE)], units={})
    for u, rows in units:
        var = "генератор" if u.gen else "привод"
        slot = next(s for s in fam["slots"] if s["model"] == u.rec["model"] and s["variant"] == var)
        r = slot["мягкий"]
        nc = 2 if r["arch"].endswith("2Ц") else 1
        PRs1 = r["PR"] if nc == 1 else r["PR"] ** 0.5
        dbl = r.get("wheel") == "двусторонний вход"
        st = E.comp_stage(E.T0, E.P0, PRs1, r["G"], "ti", dbl)
        ns = sorted({row["n_rpm"] for row in rows})
        Ts = sorted({row["T_amb_C"] for row in rows})
        rec = dict(replacement=dict(arch=r["arch"], cool=r["cool"], PR=r["PR"], T41_C=r["T41"], G_kg_s=r["G"],
                                    wheel=r.get("wheel")),
                   wheel1=dict(PR_stage=PRs1, dbl=dbl, U2_m_s=st["U2"], D2_m=st["D2"], rpm_natural=st["rpm"],
                               eta=st["eta"], ok=st["ok"]),
                   coupling=dict(n_min_rpm=ns[0], n_max_rpm=ns[-1],
                                 P_mech_max_kW=max(row["P_mech_kW"] for row in rows),
                                 M_max_Nm=max(row["M_Nm"] for row in rows)))
        if u.gen:   # одновальный заменитель: колесо на валу муфты
            rec["wheel1"]["D2_at_coupling_rpm_m"] = 60 * st["U2"] / (math.pi * u.N)
            rec["wheel1"]["rpm_ratio_coupling_to_natural"] = u.N / st["rpm"]
            rec["wheel1"]["corrected_speed_rel_by_T_amb"] = {f"{t:+.0f}": math.sqrt(T_ISO / (t + 273.15)) for t in Ts}
        out["units"][u.key] = rec
    return out


# =========================================================================== проверка вне калибровки
DIGITIZED = os.path.join(HERE, "data", "solar_available_power_taurus60.json")


def validate(u, dig=None):
    """полная мощность модели против оцифрованной кривой «Available Power» листа данных (вне точки ISO).
    GS — P_эл на клеммах (через КПД агрегата EXAMPLE); MD — мощность на валу при оптимальной частоте СТ
    (основа графика), с пределом частоты газогенератора и без него."""
    dig = dig or json.load(open(DIGITIZED))
    d = dig[u.rec["variant"]]
    pts = []
    for p in d["points"]:
        T1 = p["T_amb_C"] + 273.15
        if u.gen:
            fr, gg, pt, e = u.gs_engine(T1, u.d["T41"])
            m = u.train.el(pt["P_shaft"])[0] / 1e3
            r = dict(T_amb_C=p["T_amb_C"], P_published_kW=p["P_kW"], P_model_kW=m, err_pct=100 * (m / p["P_kW"] - 1),
                     R=fr["R"], Nc=fr["Nc"])
        else:
            a = u.md_engine_opt(T1, gg_limit=True)
            b = u.md_engine_opt(T1, gg_limit=False)
            r = dict(T_amb_C=p["T_amb_C"], P_published_kW=p["P_kW"], P_model_kW=a["P_kW"],
                     err_pct=100 * (a["P_kW"] / p["P_kW"] - 1), limit=a["limit"], n_gg_rpm=a["n_gg_rpm"],
                     P_model_no_gg_limit_kW=b["P_kW"], err_no_gg_limit_pct=100 * (b["P_kW"] / p["P_kW"] - 1),
                     n_gg_no_limit_rpm=b["n_gg_rpm"])
        if "extrapolated_C" in p:
            r["published_extrapolated_C"] = p["extrapolated_C"]
        pts.append(r)
    off = [r for r in pts if abs(r["T_amb_C"] - 15.0) > 1e-9]
    errs = [r["err_pct"] for r in off]
    return dict(method="модель (калибровка только по точке ISO паспорта) против оцифрованной кривой листа данных; "
                       "точка 15 °C — проверка оцифровки, остальные — вне калибровки",
                inputs=[os.path.relpath(DIGITIZED, HERE), d["sheet"]["url"]], source_sheet=d["sheet"]["doc"],
                basis_difference="график изготовителя: влажность 60 %, природный газ 35 МДж/нм3; модель: сухой воздух, метан",
                points=pts, off_iso_max_abs_err_pct=max(abs(x) for x in errs),
                off_iso_rms_err_pct=math.sqrt(sum(x * x for x in errs) / len(errs)))


def map_selection(cfg_gs, cfg_md):
    """выбор обобщённой карты: две открытые карты NPSS в их собственных расчётных точках и сетка положений
    (только для прозрачности) против кривой GS; кривая MD — проверка выбранного"""
    dig = json.load(open(DIGITIZED))
    cands = [("maps/npss_axi5.json", 1.0, 2.0, "axi-5 в расчётной точке файла (defaults: Nc 1,00, R 2,0)"),
             ("maps/npss_hbtf_hpc.json", 0.976, 2.05, "HBTF HPC в расчётной точке файла (defaults: Nc 0,976, R 2,05)"),
             ("maps/npss_hbtf_hpc.json", 1.0, 2.0, "HBTF HPC, Nc 1,0, R 2,0 (PR карты = defaults.PR 10,894)")]
    out = dict(method="полная мощность GS (одновальный, частота постоянна) и MD (оптимальная частота СТ) при "
                      "постоянной T41 из калибровки ISO против оцифрованных кривых листов данных",
               inputs=[os.path.relpath(DIGITIZED, HERE)], candidates=[], grid_GS=[])
    for f, Nc, R, name in cands:
        ov = {"compressor_map": {"file": f, "Nc_map_d": Nc, "Rline_d": R}}
        vg, vm = validate(Unit(cfg_gs, ov), dig), validate(Unit(cfg_md, ov), dig)
        out["candidates"].append(dict(name=name, file=f, Nc_map_d=Nc, Rline_d=R,
                                      GS=dict(rms_pct=vg["off_iso_rms_err_pct"], max_abs_pct=vg["off_iso_max_abs_err_pct"],
                                              err_pct={f"{p['T_amb_C']:+.1f}": p["err_pct"] for p in vg["points"]}),
                                      MD=dict(rms_pct=vm["off_iso_rms_err_pct"], max_abs_pct=vm["off_iso_max_abs_err_pct"],
                                              err_pct={f"{p['T_amb_C']:+.1f}": p["err_pct"] for p in vm["points"]},
                                              err_no_gg_limit_pct={f"{p['T_amb_C']:+.1f}": p["err_no_gg_limit_pct"]
                                                                   for p in vm["points"]})))
    for f in ("maps/npss_axi5.json", "maps/npss_hbtf_hpc.json"):
        for Nc in (0.975, 1.0, 1.025, 1.05):
            for R in (1.6, 2.0, 2.4):
                try:
                    v = validate(Unit(cfg_gs, {"compressor_map": {"file": f, "Nc_map_d": Nc, "Rline_d": R}}), dig)
                    out["grid_GS"].append(dict(file=f, Nc_map_d=Nc, Rline_d=R, rms_pct=v["off_iso_rms_err_pct"],
                                               max_abs_pct=v["off_iso_max_abs_err_pct"]))
                except (RuntimeError, ValueError) as ex:
                    out["grid_GS"].append(dict(file=f, Nc_map_d=Nc, Rline_d=R, error=str(ex)[:120]))
    best = min((g for g in out["grid_GS"] if "rms_pct" in g), key=lambda g: g["rms_pct"])
    out["grid_GS_best"] = best
    out["chosen"] = "axi-5 в расчётной точке файла: выбор из двух карт по кривой GS без подгонки положения точки; " \
                    "сетка положений — только для прозрачности (лучшее по GS — grid_GS_best)"
    return out


# =========================================================================== чувствительность
def sensitivity(cfg_path, temps_C=(-30.0, 15.0, 45.0)):
    """влияние допущений и тестовых значений на полную мощность (ΔP_engine, % к базе)"""
    base_u = Unit(cfg_path)
    gen = base_u.gen
    variants = [("map", "карта NPSS HBTF HPC в расчётной точке файла вместо axi-5",
                 {"compressor_map": {"file": "maps/npss_hbtf_hpc.json", "Nc_map_d": 0.976, "Rline_d": 2.05}}),
                ("map", "точка карты R = 1,6", {"compressor_map": {"Rline_d": 1.6}}),
                ("map", "точка карты R = 2,4", {"compressor_map": {"Rline_d": 2.4}}),
                ("map", "точка карты Nc = 1,025", {"compressor_map": {"Nc_map_d": 1.025}}),
                ("k_nu", "КПД турбины без поправки по U/C0 (k = 0)", {"engine": {"k_nu": 0.0}}),
                ("dp", "потери площадки: вход 1,0 кПа, выхлоп 1,49 кПа (пределы TPS60CSMD/309)",
                 {"site": {"dp_inlet_Pa": 1000.0, "dp_exhaust_Pa": 1490.0}})]
    if not gen:
        variants.append(("gg_limit", "без предела частоты газогенератора", {"engine": {"gg_speed_limit": False}}))
    res = dict(method="пересчёт полной мощности двигателя (без ограничения агрегата) с одним изменённым "
                      "допущением; база — конфигурация TOML", temps_C=list(temps_C), base={}, variants=[])
    base = {f"{t:+.0f}": base_u.full_engine(t + 273.15) for t in temps_C}
    res["base"] = base
    runs = [(k, nm, ov, "T41") for k, nm, ov in variants]
    runs.append(("control", "закон управления: T45 постоянна вместо T41", {}, "T45"))
    for key, name, ov, ctrl in runs:
        u = Unit(cfg_path, ov) if ov else base_u
        pts = {}
        for t in temps_C:
            try:
                r = u.full_engine(t + 273.15, ctrl)
            except (RuntimeError, ValueError) as ex:
                pts[f"{t:+.0f}"] = dict(error=str(ex)[:160])
                continue
            b = base[f"{t:+.0f}"]
            if r["R"] < u.map.R_stall + 0.1:
                r["note"] = "рабочая точка у границы помпажа карты (R < 1,1): запаса нет"
            r["dP_pct"] = 100 * (r["P_engine_kW"] / b["P_engine_kW"] - 1)
            r["dn_pct"] = 100 * (r["n_rpm"] / b["n_rpm"] - 1)
            pts[f"{t:+.0f}"] = r
        res["variants"].append(dict(key=key, name=name, overrides=ov, control=ctrl, points=pts))
    # аналитические: КПД агрегата и частота СТ
    if gen:
        tr = base_u.train
        res["eta_gen"] = dict(name="КПД генератора ±0,01", dP_mech_pct=[100 * (tr.eg / (tr.eg + 0.01) - 1),
                                                                         100 * (tr.eg / (tr.eg - 0.01) - 1)],
                              note="P_мех = P_эл/(eta_ред*eta_ген): калибровка ведётся по P_эл паспорта")
        res["eta_red"] = dict(name="КПД редуктора ±0,01", dP_mech_pct=[100 * (tr.er / (tr.er + 0.01) - 1),
                                                                        100 * (tr.er / (tr.er - 0.01) - 1)])
    else:
        out = {}
        for f in (0.95, 1.05):
            u = Unit(cfg_path, {"train": {"n_pt_d_rpm": base_u.cfg["train"]["n_pt_d_rpm"] | {"value": base_u.proc.n_d * f}}})
            r = u.full_engine(T_ISO)
            out[f"{f:.2f}"] = dict(n_rpm=r["n_rpm"], P_engine_kW=r["P_engine_kW"],
                                   M_Nm=r["P_engine_kW"] * 1e3 / rpm2w(r["n_rpm"]))
        b = base["+15"]
        Mb = b["P_engine_kW"] * 1e3 / rpm2w(b["n_rpm"])
        res["n_pt_d"] = dict(name="расчётная частота СТ ±5 % (ISO, полная мощность)", base_M_Nm=Mb, points=out,
                             dM_pct=[100 * (out[k]["M_Nm"] / Mb - 1) for k in out])
    return res


# =========================================================================== 5. выгрузка
CSV_COLS = ["regime", "T_amb_C", "load_pct", "n_rpm", "P_mech_kW", "M_Nm", "P_el_kW", "W_air_kg_s", "W_fuel_kg_s",
            "T_exh_C", "eta_mech_pct", "eta_el_pct", "P_engine_max_kW", "limit", "T41_C", "T45_C", "PR", "Nc_rel",
            "Rline", "W_exh_kg_s", "eta_red", "eta_gen", "n_gen_rpm", "n_gg_rpm", "Q_m3s", "H_kJkg", "PR_proc",
            "surge_margin", "mass_resid_kg_s", "energy_resid_pct_fuel", "flags"]
FMT = dict(T_amb_C=1, load_pct=0, n_rpm=1, P_mech_kW=1, M_Nm=1, P_el_kW=1, W_air_kg_s=3, W_fuel_kg_s=4, T_exh_C=1,
           eta_mech_pct=2, eta_el_pct=2, P_engine_max_kW=1, T41_C=1, T45_C=1, PR=3, Nc_rel=4, Rline=3, W_exh_kg_s=3,
           eta_red=4, eta_gen=4, n_gen_rpm=0, n_gg_rpm=0, Q_m3s=4, H_kJkg=2, PR_proc=4, surge_margin=3)


def write_csv(path, rows):
    with open(path, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, lineterminator="\n")
        w.writerow(CSV_COLS)
        for r in rows:
            out = []
            for c in CSV_COLS:
                v = r.get(c)
                if v is None:
                    out.append("")
                elif c in ("mass_resid_kg_s", "energy_resid_pct_fuel"):
                    out.append(f"{v:.1e}")
                elif isinstance(v, float):
                    out.append(f"{v:.{FMT.get(c, 4)}f}")
                else:
                    out.append(str(v))
            w.writerow(out)


def identity_check():
    """мой расчёт цикла против gtu_cycle.cycle на калибровке Taurus 60 MD из gtu_study_results.json"""
    p = os.path.join(GTU, "gtu_study_results.json")
    s = json.load(open(p))["T60"]
    kw = dict(eta_pc=s["eta_pc"], eta_pt_gg=s["eta_pt_gg"], eta_pt_pt=s["eta_pt_pt"], x_ngv=s["x_ngv"],
              x_rot=s["x_rot"], x_leak=s["x_leak"])
    W1 = 21.63 / 1.004                      # как gtu_study.py (calib «Taurus 60, привод»)
    ref = C.cycle(s["PR"], s["T41_C"], W1=W1, **kw)
    e = dict(x_ngv=s["x_ngv"], x_rot=s["x_rot"], x_leak=s["x_leak"], dp_comb=0.04, dp_duct=0.01, eta_b=0.995,
             eta_m=0.99, T_fuel=T_ISO, k_nu=0.0, eta_gg=s["eta_pt_gg"], eta_pt=s["eta_pt_pt"], P5=P_ISO * 1.02)
    fr = cycle_front(T_ISO, P_ISO, W1, s["PR"], s["eta_pc"], s["T41_C"] + 273.15, e)
    gg = cycle_ggt(fr, e)
    pt = cycle_pt(gg, fr, e)
    mine = dict(P_MW=pt["P_shaft"] / 1e6, eff=pt["eff"], EGT_C=pt["T5"] - 273.15, W_exh=gg["W45"],
                T3_C=fr["T3"] - 273.15, T45_C=gg["T45"] - 273.15, Wf=fr["Wf"])
    return dict(method="тот же вход в gtu_cycle.cycle и в cycle_front/cycle_ggt/cycle_pt этого модуля",
                inputs=[os.path.relpath(p, HERE) + " → T60", "W1 = 21,63/1,004 кг/с (gtu_study.py)"],
                gtu_cycle={k: ref[k] for k in mine}, gtu_shaft=mine,
                rel_diff={k: mine[k] / ref[k] - 1 for k in mine},
                gtu_study_T60=dict(P_MW=s["P_MW"], eff=s["eff"], EGT_C=s["EGT_C"]))


def data_needed(cfg_paths, sens):
    rows = []
    for p in cfg_paths:
        cfg = tomllib.load(open(p, "rb"))
        key, var = cfg["unit"]["key"], cfg["unit"]["variant"]
        items = []
        for sec in ("engine", "train", "site"):
            for k, v in cfg.get(sec, {}).items():
                if isinstance(v, dict) and v.get("example"):
                    items.append(dict(parameter=f"{sec}.{k}", unit=v["unit"], owner=v["owner"], impact=v["impact"],
                                      replaced_by=v["replaced_by"], sens_key=v.get("sens_key")))
        items += cfg.get("data_needed", [])
        s = sens.get(key, {})
        for it in items:
            imp = it["impact"]
            sk = it.get("sens_key")
            if sk in ("map", "k_nu", "dp", "control", "gg_limit"):
                parts = []
                for v in s.get("variants", []):
                    if v["key"] == sk:
                        pts = ", ".join(f"{t} °C: {v['points'][t]['dP_pct']:+.1f} %" if "dP_pct" in v["points"][t]
                                        else f"{t} °C: нет решения" for t in v["points"])
                        parts.append(f"{v['name']} — {pts}")
                if parts:
                    imp += "; расчёт ΔP полной мощности (shaft_results.json → sensitivity): " + "; ".join(parts)
            elif sk in ("eta_gen", "eta_red") and sk in s:
                a, b = s[sk]["dP_mech_pct"]
                imp += f"; расчёт: ±0,01 → P_мех и M на муфте {a:+.2f} / {b:+.2f} %"
            elif sk == "n_pt_d" and "n_pt_d" in s:
                a, b = s["n_pt_d"]["dM_pct"]
                imp += f"; расчёт: n_d ±5 % → M на муфте при ISO {b:+.1f} / {a:+.1f} %"
            rows.append([var, it["parameter"], it["unit"], it["owner"], imp, it["replaced_by"]])
    return rows


def sha(path):
    return hashlib.sha256(open(path, "rb").read()).hexdigest()


def run(cfg_paths, out_dir, with_sens=True):
    t0 = datetime.datetime.now()
    units, res_units, sens = [], {}, {}
    for p in cfg_paths:
        u = Unit(p)
        rows = u.rows()
        write_csv(os.path.join(out_dir, f"shaft_{u.key}.csv"), rows)
        units.append((u, rows))
        res_units[u.key] = dict(title=u.cfg["unit"]["title"], layout=u.layout, config=os.path.basename(p),
                                csv=f"shaft_{u.key}.csv", **u.info(), rows=rows)
        if with_sens:
            sens[u.key] = sensitivity(p)
            res_units[u.key]["sensitivity"] = sens[u.key]
    for u, rows in units:
        res_units[u.key]["checks"] = checks(u, rows)
        res_units[u.key]["envelope"] = envelope(rows, u.gen)
        res_units[u.key]["validation"] = validate(u)
    msel = None
    if with_sens:
        lay = {u.layout: u.cfg_path for u, _ in units}
        if "single_shaft_gen" in lay and "two_shaft_md" in lay:
            msel = map_selection(lay["single_shaft_gen"], lay["two_shaft_md"])
    if with_sens:
        with open(os.path.join(out_dir, "data_needed.csv"), "w", newline="", encoding="utf-8") as fh:
            w = csv.writer(fh, lineterminator="\n")
            w.writerow(["исполнение", "параметр", "единицы СИ", "владелец", "влияние", "чем заменено сейчас"])
            w.writerows(data_needed(cfg_paths, sens))
    files = list(cfg_paths) + [u.pp_path for u, _ in units] + [os.path.join(HERE, "maps", f) for f in
                                                                ("npss_hbtf_hpc.json", "npss_axi5.json")]
    files += [os.path.join(GTU, f) for f in ("gtu_cycle.py", "gtu_explore.py", "gtu_family_results.json",
                                             "gtu_study_results.json")] + [DIGITIZED, os.path.abspath(__file__)]
    out = dict(
        method="баланс мощности на муфте выходного вала по режимам: модель цикла gtu_cycle (Cantera GRI-Mech 3.0, "
               "метан, политропные шаги) с входом T_нв; согласование узлов вне ISO — обобщённая карта компрессора "
               "ГТД (NPSS axi-5 из pyCycle 4.4.0, масштаб на точку ISO — допущение), закон эллипса Стодолы для "
               "турбин, КПД турбин по U/C0 (допущение); пределы — постоянная T41 из калибровки по паспорту ISO и "
               "(MD) паспортная частота газогенератора; одновальный — постоянная частота, двухвальный — частота СТ "
               "по нагрузке техкомпрессора (EXAMPLE); проверка вне ISO — оцифрованные кривые листов данных Solar",
        software=f"python {platform.python_version()}, cantera {ct.__version__}, numpy {np.__version__}, "
                 f"scipy {scipy.__version__}",
        date=DATE, inputs=[dict(file=os.path.relpath(f, HERE), sha256=sha(f)) for f in files],
        example_note="все параметры агрегата (редуктор, генератор, техкомпрессор, газ, сеть) — тестовые значения "
                     "EXAMPLE: example=true в TOML, перечень — data_needed.csv",
        assumptions=[
            "характеристика компрессора ГТД — чужая обобщённая карта NPSS axi-5 (pyCycle 4.4.0, Apache-2.0) в её "
            "расчётной точке, масштабированная на точку ISO; отношение КПД карты применено к политропному КПД; "
            "выбор из двух карт — по кривой GS листа данных (map_selection)",
            "MD: точка ISO лежит на пределе частоты газогенератора 15 000 об/мин (по излому опубликованной кривой "
            "при 15 °C); ниже 15 °C полная мощность ограничена частотой газогенератора",
            "ВНА и поворотные НА неподвижны в расчётном положении на всех режимах",
            "пропускная способность турбин — закон эллипса Стодолы, калибровка по точке ISO",
            "КПД турбин вне расчётной частоты — парабола по U/C0 (k = 1)",
            "полная мощность — при постоянной T41 (калибровочной при ISO); у двухвального ещё вариант T45 = const "
            "в оценке чувствительности",
            "воздух сухой, топливо — метан при 288,15 К, давление площадки 101 325 Па, потерь на входе и выхлопе нет "
            "(основа паспорта ISO)",
            "частичная нагрузка — доля полной мощности при той же T_нв"],
        identity_check=identity_check(), units=res_units, stage_link=stage_link(units), map_selection=msel,
        elapsed_s=(datetime.datetime.now() - t0).total_seconds())
    with open(os.path.join(out_dir, "shaft_results.json"), "w", encoding="utf-8") as fh:
        json.dump(out, fh, ensure_ascii=False, indent=1, default=_json_default)
    return out


def _json_default(o):
    if isinstance(o, (np.floating,)):
        return float(o)
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, np.ndarray):
        return o.tolist()
    raise TypeError(type(o))


def checks(u, rows):
    iso = next(r for r in rows if abs(r["T_amb_C"] - 15.0) < 1e-9 and r["load_pct"] >= 100.0)
    pp = u.pp
    c = dict(iso_vs_passport=dict(
        P_model_kW=iso["P_el_kW"] if u.gen else iso["P_mech_kW"], P_passport_kW=pp["P_out_kW"],
        eff_model_pct=iso["eta_el_pct"] if u.gen else iso["eta_mech_pct"], eff_passport_pct=pp["eff_pct"],
        T_exh_model_C=iso["T_exh_C"], T_exh_passport_C=pp["exh_temp_C"],
        W_exh_model_kg_s=iso["W_exh_kg_s"], W_exh_passport_kg_s=pp["exh_flow_kg_s"],
        note="P и КПД совпадают по построению (калибровка); T и расход выхлопа — проверка"))
    full = sorted((r for r in rows if r["load_pct"] >= 100.0), key=lambda r: r["T_amb_C"])
    c["P_engine_monotonic_decreasing"] = all(a["P_engine_max_kW"] > b["P_engine_max_kW"] for a, b in zip(full, full[1:]))
    c["max_abs_mass_resid_kg_s"] = max(abs(r["mass_resid_kg_s"]) for r in rows)
    c["max_abs_energy_resid_pct_fuel"] = max(abs(r["energy_resid_pct_fuel"]) for r in rows)
    c["max_abs_M_omega_minus_P_rel"] = max(abs(r["M_Nm"] * rpm2w(r["n_rpm"]) / (r["P_mech_kW"] * 1e3) - 1) for r in rows)
    if u.gen:
        c["max_abs_Pel_identity_rel"] = max(abs(r["P_el_kW"] / (r["P_mech_kW"] * r["eta_red"] * r["eta_gen"]) - 1)
                                            for r in rows)
    c["rows_with_flags"] = [dict(regime=r["regime"], flags=r["flags"]) for r in rows if r["flags"]]
    return c


def envelope(rows, gen):
    """требования к заменителю на муфте: огибающая по режимам"""
    rM = max(rows, key=lambda r: r["M_Nm"])
    rP = max(rows, key=lambda r: r["P_mech_kW"])
    e = dict(P_mech_max_kW=rP["P_mech_kW"], P_mech_max_at=rP["regime"], M_max_Nm=rM["M_Nm"], M_max_at=rM["regime"],
             n_min_rpm=min(r["n_rpm"] for r in rows), n_max_rpm=max(r["n_rpm"] for r in rows),
             P_mech_min_kW=min(r["P_mech_kW"] for r in rows))
    if gen:
        e["P_el_max_kW"] = max(r["P_el_kW"] for r in rows)
    return e


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("configs", nargs="*", help="TOML-конфиги (по умолчанию оба Taurus 60)")
    ap.add_argument("--out", default=HERE, help="папка выгрузки (по умолчанию рядом со скриптом)")
    ap.add_argument("--no-sens", action="store_true", help="без оценок чувствительности и data_needed.csv")
    a = ap.parse_args()
    cfgs = [os.path.abspath(p) for p in a.configs] or [os.path.join(HERE, c) for c in CONFIGS]
    os.makedirs(a.out, exist_ok=True)
    res = run(cfgs, os.path.abspath(a.out), with_sens=not a.no_sens)
    for k, u in res["units"].items():
        print(f"== {u['title']}  ({u['csv']})")
        print(f"{'режим':>10} {'n':>8} {'P_мех':>8} {'M':>8} {'P_эл':>8} {'T41':>7} {'Tвых':>6} {'КПД':>6}  ограничение")
        for r in u["rows"]:
            pe = f"{r['P_el_kW']:8.0f}" if r["P_el_kW"] is not None else "       —"
            print(f"{r['regime']:>10} {r['n_rpm']:8.0f} {r['P_mech_kW']:8.0f} {r['M_Nm']:8.1f} {pe} "
                  f"{r['T41_C']:7.1f} {r['T_exh_C']:6.1f} {r['eta_mech_pct']:6.2f}  {r['limit']} {r['flags']}")
        print("проверки:", json.dumps({k2: v for k2, v in u["checks"].items() if k2 != "rows_with_flags"},
                                      ensure_ascii=False, default=float))
    print("совпадение с gtu_cycle.cycle:", {k: f"{v:.2e}" for k, v in res["identity_check"]["rel_diff"].items()})
    print(f"время: {res['elapsed_s']:.0f} с")
