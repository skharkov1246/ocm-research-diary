"""Предварительный расчёт простого цикла двухвальной ГТУ со свободной силовой турбиной.
Свойства газа - Cantera (gri30), топливо - метан, условия ISO (15 C, 101.325 кПа), без потерь на входе/выхлопе.
Назначение: калибровка по паспортам аналогов и выбор целевых параметров. Уровень - облик (0D)."""
import cantera as ct, json, sys
from scipy.optimize import brentq
import os as _os; _os.chdir(_os.path.dirname(_os.path.abspath(__file__)))

AIR = "O2:0.2095,N2:0.7808,AR:0.0093,CO2:0.0004"
P0, T0 = 101325.0, 288.15
LHV = 50.03e6  # метан, Дж/кг


def gas_air():
    g = ct.Solution("gri30.yaml")
    g.TPX = T0, P0, AIR
    return g


def poly_process(g, p2, eta, compress=True, n=40):
    """Политропный процесс малыми ступенями: возвращает работу, Дж/кг (положит.)."""
    h_start = g.enthalpy_mass
    p1 = g.P
    r = (p2 / p1) ** (1.0 / n)
    X = g.X
    for _ in range(n):
        s, h, p = g.entropy_mass, g.enthalpy_mass, g.P
        g.SP = s, p * r
        dh_is = g.enthalpy_mass - h
        dh = dh_is / eta if compress else dh_is * eta
        g.HPX = h + dh, p * r, X      # состав заморожен
    return abs(g.enthalpy_mass - h_start)


def cycle(PR, T41_C, eta_pc, eta_pt_gg, eta_pt_pt, x_ngv, x_rot, x_leak,
          dp_comb=0.04, dp_duct=0.01, p_exh_ratio=1.02, eta_b=0.995, eta_m=0.99, W1=1.0):
    """W1 - расход на входе, кг/с. x_* - доли W1: охлаждение соплового аппарата (возврат до ротора 1),
    охлаждение ротора/дисков (возврат после ТВД, работы в ТВД не даёт), утечки/уплотнения (потеряны).
    T41 - температура газа перед рабочим колесом 1-й ступени (после подмешивания охлаждения СА)."""
    T41 = T41_C + 273.15
    a = gas_air()
    h1 = a.enthalpy_mass
    wc = poly_process(a, P0 * PR, eta_pc, True)
    T3, h3, p3 = a.T, a.enthalpy_mass, a.P
    Wcomb_air = W1 * (1 - x_ngv - x_rot - x_leak)

    def after_mix(far):
        """состояние перед ротором 1: продукты сгорания + воздух охлаждения СА"""
        g = ct.Solution("gri30.yaml")
        # сжигание: воздух T3 + метан 288 K -> равновесие без диссоциации заменяем полным сгоранием
        f = ct.Quantity(g, constant="HP"); f.TPX = T0, p3, "CH4:1"; f.mass = far * Wcomb_air
        ai = ct.Quantity(g, constant="HP"); ai.TPX = T3, p3, AIR; ai.mass = Wcomb_air
        m = ai + f
        hm = (ai.enthalpy_mass * ai.mass + f.enthalpy_mass * f.mass) / (ai.mass + f.mass)
        # полное сгорание CH4 + 2 O2 -> CO2 + 2 H2O (с учётом полноты eta_b - недожог как потеря тепла)
        Y = m.Y.copy(); gi = g.species_index
        yf = Y[gi("CH4")]
        MW = g.molecular_weights
        nf = yf / MW[gi("CH4")]
        Y[gi("CH4")] = 0
        Y[gi("O2")] -= 2 * nf * MW[gi("O2")]
        Y[gi("CO2")] += nf * MW[gi("CO2")]
        Y[gi("H2O")] += 2 * nf * MW[gi("H2O")]
        p4 = p3 * (1 - dp_comb)
        g.HPY = hm - (1 - eta_b) * far / (1 + far) * LHV, p4, Y
        T4 = g.T
        q4 = ct.Quantity(g, constant="HP"); q4.mass = Wcomb_air * (1 + far)
        g2 = ct.Solution("gri30.yaml")
        c = ct.Quantity(g2, constant="HP"); c.TPX = T3, p4, AIR; c.mass = W1 * x_ngv
        mix = q4 + c
        return mix, T4

    far = brentq(lambda f: after_mix(f)[0].T - T41, 0.005, 0.05, xtol=1e-7)
    mix, T4 = after_mix(far)
    g = ct.Solution("gri30.yaml"); g.TPY = mix.T, mix.P, mix.Y
    W41 = mix.mass
    Wf = far * Wcomb_air
    # ТВД: работа = работа компрессора / мех. КПД
    need = wc * W1 / eta_m / W41
    def gg_resid(p_out):
        g.TPY = mix.T, mix.P, mix.Y
        return poly_process(g, p_out, eta_pt_gg, False) - need
    p45 = brentq(gg_resid, P0 * 1.05, mix.P * 0.95, xtol=1.0)
    gg_resid(p45)
    # подмешивание охлаждения ротора после ТВД
    q = ct.Quantity(g, constant="HP"); q.mass = W41
    g2 = ct.Solution("gri30.yaml")
    c = ct.Quantity(g2, constant="HP"); c.TPX = T3, p45, AIR; c.mass = W1 * x_rot
    m2 = q + c
    g.TPY = m2.T, p45 * (1 - dp_duct), m2.Y
    T45 = g.T
    W45 = m2.mass
    w_pt = poly_process(g, P0 * p_exh_ratio, eta_pt_pt, False)
    T5 = g.T
    P_shaft = w_pt * W45 * eta_m
    return dict(PR=PR, T41_C=T41_C, T4_C=T4 - 273.15, T3_C=T3 - 273.15, T45_C=T45 - 273.15,
                EGT_C=T5 - 273.15, W_exh=W45, Wf=Wf, P_MW=P_shaft / 1e6,
                eff=P_shaft / (Wf * LHV), sp_work_kJkg=P_shaft / W1 / 1e3,
                GG_ER=mix.P / p45, PT_ER=p45 * (1 - dp_duct) / (P0 * p_exh_ratio),
                wc_kJkg=wc / 1e3)


def cycle_r(PR, T41_C, eta_pc, eta_pt_gg, eta_pt_pt, x_ngv, x_rot, x_leak,
          dp_comb=0.04, dp_duct=0.01, p_exh_ratio=1.02, eta_b=0.995, eta_m=0.99, W1=1.0,
            T3r=None, dp_rec_air=0.0):
    """W1 - расход на входе, кг/с. x_* - доли W1: охлаждение соплового аппарата (возврат до ротора 1),
    охлаждение ротора/дисков (возврат после ТВД, работы в ТВД не даёт), утечки/уплотнения (потеряны).
    T41 - температура газа перед рабочим колесом 1-й ступени (после подмешивания охлаждения СА)."""
    T41 = T41_C + 273.15
    a = gas_air()
    h1 = a.enthalpy_mass
    wc = poly_process(a, P0 * PR, eta_pc, True)
    T3, h3, p3 = a.T, a.enthalpy_mass, a.P
    Wcomb_air = W1 * (1 - x_ngv - x_rot - x_leak)

    def after_mix(far):
        """состояние перед ротором 1: продукты сгорания + воздух охлаждения СА"""
        g = ct.Solution("gri30.yaml")
        # сжигание: воздух T3 + метан 288 K -> равновесие без диссоциации заменяем полным сгоранием
        f = ct.Quantity(g, constant="HP"); f.TPX = T0, p3, "CH4:1"; f.mass = far * Wcomb_air
        ai = ct.Quantity(g, constant="HP"); ai.TPX = (T3r or T3), p3, AIR; ai.mass = Wcomb_air
        m = ai + f
        hm = (ai.enthalpy_mass * ai.mass + f.enthalpy_mass * f.mass) / (ai.mass + f.mass)
        # полное сгорание CH4 + 2 O2 -> CO2 + 2 H2O (с учётом полноты eta_b - недожог как потеря тепла)
        Y = m.Y.copy(); gi = g.species_index
        yf = Y[gi("CH4")]
        MW = g.molecular_weights
        nf = yf / MW[gi("CH4")]
        Y[gi("CH4")] = 0
        Y[gi("O2")] -= 2 * nf * MW[gi("O2")]
        Y[gi("CO2")] += nf * MW[gi("CO2")]
        Y[gi("H2O")] += 2 * nf * MW[gi("H2O")]
        p4 = p3 * (1 - dp_comb) * (1 - dp_rec_air)
        g.HPY = hm - (1 - eta_b) * far / (1 + far) * LHV, p4, Y
        T4 = g.T
        q4 = ct.Quantity(g, constant="HP"); q4.mass = Wcomb_air * (1 + far)
        g2 = ct.Solution("gri30.yaml")
        c = ct.Quantity(g2, constant="HP"); c.TPX = T3, p4, AIR; c.mass = W1 * x_ngv
        mix = q4 + c
        return mix, T4

    far = brentq(lambda f: after_mix(f)[0].T - T41, 0.005, 0.05, xtol=1e-7)
    mix, T4 = after_mix(far)
    g = ct.Solution("gri30.yaml"); g.TPY = mix.T, mix.P, mix.Y
    W41 = mix.mass
    Wf = far * Wcomb_air
    # ТВД: работа = работа компрессора / мех. КПД
    need = wc * W1 / eta_m / W41
    def gg_resid(p_out):
        g.TPY = mix.T, mix.P, mix.Y
        return poly_process(g, p_out, eta_pt_gg, False) - need
    p45 = brentq(gg_resid, P0 * 1.05, mix.P * 0.95, xtol=1.0)
    gg_resid(p45)
    # подмешивание охлаждения ротора после ТВД
    q = ct.Quantity(g, constant="HP"); q.mass = W41
    g2 = ct.Solution("gri30.yaml")
    c = ct.Quantity(g2, constant="HP"); c.TPX = T3, p45, AIR; c.mass = W1 * x_rot
    m2 = q + c
    g.TPY = m2.T, p45 * (1 - dp_duct), m2.Y
    T45 = g.T
    W45 = m2.mass
    w_pt = poly_process(g, P0 * p_exh_ratio, eta_pt_pt, False)
    T5 = g.T
    P_shaft = w_pt * W45 * eta_m
    return dict(PR=PR, T41_C=T41_C, T4_C=T4 - 273.15, T3_C=T3 - 273.15, T45_C=T45 - 273.15,
                EGT_C=T5 - 273.15, W_exh=W45, Wf=Wf, P_MW=P_shaft / 1e6,
                eff=P_shaft / (Wf * LHV), sp_work_kJkg=P_shaft / W1 / 1e3,
                GG_ER=mix.P / p45, PT_ER=p45 * (1 - dp_duct) / (P0 * p_exh_ratio),
                wc_kJkg=wc / 1e3)




def recup(PR, T41_C, eps, dp_air=0.025, dp_gas=0.035, **kw):
    """Цикл с рекуператором: eps - степень регенерации; итерации по температуре выхлопа."""
    pe = kw.pop("p_exh_ratio", 1.02) * (1 + dp_gas)
    r = cycle_r(PR, T41_C, p_exh_ratio=pe, dp_rec_air=dp_air, **kw)
    for _ in range(6):
        T3 = r["T3_C"] + 273.15; T5 = r["EGT_C"] + 273.15
        T3r = T3 + eps * (T5 - T3)
        r = cycle_r(PR, T41_C, p_exh_ratio=pe, dp_rec_air=dp_air, T3r=T3r, **kw)
    r["T3r_C"] = T3r - 273.15
    r["Tstack_C"] = r["EGT_C"] - (T3r - T3) * 0.93     # оценка: отношение теплоёмкостей потоков
    return r

if __name__ == "__main__":
    out = {}
    # --- 1. калибровка по аналогам: подбираем T41 под паспортную мощность, сравниваем КПД и Т выхлопа ---
    def fit(name, PR, W1, P_target, **kw):
        f = lambda t: cycle(PR, t, W1=W1, **kw)["P_MW"] - P_target
        t = brentq(f, 900, 1400, xtol=0.5)
        r = cycle(PR, t, W1=W1, **kw); r["name"] = name
        return r
    # SGT-400 (13 МВт, привод): 13.4 МВт, 36.2 %, PR 16.8, выхлоп 39.4 кг/с, 555 C
    out["SGT400"] = fit("SGT-400 13 МВт привод", 16.8, 39.4 / 1.006, 13.4,
                        eta_pc=0.895, eta_pt_gg=0.875, eta_pt_pt=0.895, x_ngv=0.055, x_rot=0.055, x_leak=0.015)
    # Taurus 60 привод: 5.74 МВт, 31.96 %, PR 12.2, выхлоп 21.63 кг/с, 510 C
    out["T60"] = fit("Taurus 60 привод", 12.2, 21.63 / 1.004, 5.74,
                     eta_pc=0.875, eta_pt_gg=0.865, eta_pt_pt=0.885, x_ngv=0.045, x_rot=0.045, x_leak=0.015)
    for k in ("SGT400", "T60"):
        r = out[k]
        print(f"{r['name']}: T41={r['T41_C']:.0f}C T4={r['T4_C']:.0f}C  P={r['P_MW']:.2f} МВт  КПД={r['eff']*100:.1f}%  Tвых={r['EGT_C']:.0f}C  Gвых={r['W_exh']:.1f}")
    out.update({"method": "0D-цикл Брайтона: свойства газа Cantera (GRI-Mech 3.0), метан, ISO; КПД узлов и ограничения - принятые, см. README", "software": "cantera " + __import__("cantera").__version__ + ", python3", "date": __import__("datetime").date.today().isoformat()}); json.dump(out, open("gtu_cycle_results.json", "w"), ensure_ascii=False, indent=1)
