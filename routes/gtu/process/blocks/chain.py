"""Сквозная термодинамическая цепочка «топливо -> вал» по блокам (уровень облика, 0D; всё - оценка).

Свойства газа - таблицы gtu_explore (Cantera GRI-Mech 3.0; воздух и стехиометрические продукты метана, смесь по доле
продуктов b) через funnel_space.FastGas - те же, что у быстрого цикла воронки. Логика основного ядра повторяет
funnel_space._core (тот же порядок: охлаждение соплового аппарата до ротора, охлаждение ротора после турбины, утечка
теряется; мех. КПД 0,99 на работе турбин; полнота сгорания 0,995 - добавка топлива), поэтому при тех же КПД узлов и без
добавочных блоков цепочка даёт то же, что воронка (проверяется тестом). Сверх воронки:
  * промежуточное охлаждение воздух-воздух с вентилятором (как organic/cycles: cool_with_air);
  * подогрев топлива выхлопом (как cycles make_fuel('fh'));
  * рекуператор с эффективностью по наименьшей теплоёмкости (энтальпийная запись, как cycles.recuperate) или по
    температуре воздуха (как funnel_space._core и gtu_cycle.recup) - rec_model;
  * повторный подогрев: на одном валу - деление расширения по ln p (x_rh, как cycles), у свободной силовой турбины -
    межтурбинная камера в переходном канале;
  * горение с ростом давления: коэффициент давления камеры PRb задаётся снаружи (блок камеры);
  * обращённый цикл Брайтона на выхлопе (как cycles.ibc_run, без конденсации воды - при охлаждении до T0 + 15 % разности
    точка росы обычно не достигается; сверка с перебором циклов - тест).
Химическая рекуперация (паровой и сухой риформинг) меняет состав топлива - таблицы её не описывают; её приращение
переносится из перебора циклов (Cantera) в блоке топлива, здесь не считается.
Все величины - на 1 кг/с воздуха на входе компрессора.
"""
import math

import cantera as ct
import numpy as np

from . import _track as T

TG_CH4 = np.arange(250.0, 1300.1, 5.0)
_HCH4 = None


def h_ch4(TK):
    """энтальпия метана, Дж/кг (Cantera GRI-Mech 3.0, идеальный газ) - для подогрева топлива"""
    global _HCH4
    if _HCH4 is None:
        g = ct.Solution("gri30.yaml")
        hs = []
        for t in TG_CH4:
            g.TPX = t, ct.one_atm, "CH4:1"
            hs.append(g.enthalpy_mass)
        _HCH4 = np.array(hs)
    return float(np.interp(TK, TG_CH4, _HCH4))


def far_gen(E, T_in, b_in, T_out, hf):
    """топливо на 1 кг газа состава b_in, чтобы нагреть его от T_in до T_out полным сгоранием; hf - энтальпия топлива.
    Баланс: df·[hf + AFST·h_a(T_out) - (1+AFST)·h_s(T_out)] = h(T_out, b_in) - h(T_in, b_in). Для чистого воздуха и
    hf = HF совпадает с FastGas.far_for."""
    den = hf + E.AFST * E.h(T_out) - (1 + E.AFST) * E.h(T_out, 1.0)
    return (E.h(T_out, b_in) - E.h(T_in, b_in)) / den


def b_after(E, b_in, df):
    """доля продуктов после сжигания df кг топлива на 1 кг газа состава b_in"""
    return (b_in + (1 + E.AFST) * df) / (1 + df)


def mix(E, streams):
    """смешение потоков [(m, T, b)] при одном давлении: (m, T, b)"""
    m = sum(s[0] for s in streams)
    if m <= 0:
        return 0.0, streams[0][1], streams[0][2]
    b = sum(s[0] * s[2] for s in streams) / m
    H = sum(s[0] * E.h(s[1], s[2]) for s in streams)
    return m, E.T_h(H / m, b), b


def expand(E, T1, b, ER, eta_p):
    """политропное расширение (как funnel_space._core): T2 и работа Дж/кг"""
    T2 = E.T_s(E.s0(T1, b) - eta_p * E.Rg(b) * math.log(ER), b)
    return T2, E.h(T1, b) - E.h(T2, b)


def compress(E, T1, pi, eta_p):
    T2 = E.T_s(E.s0(T1) + E.RA * math.log(pi) / eta_p)
    return T2, E.h(T2) - E.h(T1)


def fan_power(E, Q, T_hot_in, ic):
    """вентилятор охладителя воздух-воздух (как cycles.cool_with_air): расход холодного воздуха Cr·Q/ΔT/cp, перепад dp_fan"""
    cp_air = 1006.0
    m_c = ic["Cr"] * Q / max(T_hot_in - E.T0, 1.0) / cp_air
    rho_c = E.P0 / (ct.gas_constant / 28.96 * E.T0)
    return m_c * ic["dp_fan"] / rho_c / (ic["eta_fan"] * ic["eta_motor"]), m_c


def burn_to(E, m_g, T_in, b_in, T41, x_cool, T_cool, hf):
    """камера: топливо на поток m_g (состав b_in, T_in), чтобы после подмешивания x_cool воздуха охлаждения соплового
    аппарата (T_cool) газ имел T41 (как funnel_space._core). Возврат: df (на 1 кг газа), T4, m41, b41"""
    T4 = T41
    df, b4, m41, b41 = 0.0, b_in, m_g, b_in
    for _ in range(30):
        df = far_gen(E, T_in, b_in, T4, hf)
        b4 = b_after(E, b_in, df)
        mg4 = m_g * (1 + df)
        m41 = mg4 + x_cool
        b41 = b4 * mg4 / m41
        T4n = E.T_h((m41 * E.h(T41, b41) - x_cool * E.h(T_cool)) / mg4, b4)
        if abs(T4n - T4) < 1e-5:
            T4 = T4n
            break
        T4 = T4n
    return df, T4, m41, b41


def run(spec):
    """spec - словарь (см. blocks.evaluate.chain_spec):
      stages [(pi, eta_p)], ic_after [индексы ступеней, после которых охладитель], ic {eps, dp, Cr, dp_fan, eta_fan, eta_motor},
      xl, xn, xr, xn2, xr2 - доли воздуха (утечка; охлаждение соплового аппарата и ротора основной и повторной ступени),
      rec {eps, dpa, dpg, model: 'enthalpy'|'temperature'} или None, PRb - коэффициент полного давления камеры (1 - потери
      или больше 1 при росте давления), t41 (°C), shaft 'single'|'free', eta_gg, eta_pt (политропные), dp_duct (свободная СТ),
      reheat {x_rh, t_rh (°C), dp} или None, fuel {T_max (°C), dT_app, dp_g} или None, ibc {...} или None,
      p_exh 1,02, eta_m 0,99, eta_b 0,995.
    Возврат: словарь показателей на 1 кг/с воздуха (ok=False и why - если цепочка не замыкается)."""
    E = T.E()
    eta_m, eta_b, p_exh = spec.get("eta_m", 0.99), spec.get("eta_b", 0.995), spec.get("p_exh", 1.02)
    # ---------------------------------------------------------------- компрессор и промохлаждение
    Tc, pc, wc, W_fan, Q_ic, m_fan_air = E.T0, E.P0, 0.0, 0.0, 0.0, 0.0
    st_out, ic_out = [], []
    for k, (pi, eta) in enumerate(spec["stages"]):
        T2, w = compress(E, Tc, pi, eta)
        st_out.append(dict(T1=Tc, p1=pc, T2=T2, p2=pc * pi, w=w, pi=pi, eta=eta))
        wc += w
        Tc, pc = T2, pc * pi
        if k in spec.get("ic_after", ()) and k < len(spec["stages"]) - 1:
            ic = spec["ic"]
            Tin = Tc
            Tout = Tin - ic["eps"] * (Tin - E.T0)
            Q = E.h(Tin) - E.h(Tout)
            Wf, mf_air = fan_power(E, Q, Tin, ic)
            W_fan += Wf
            Q_ic += Q
            m_fan_air += mf_air
            ic_out.append(dict(after=k, T_in=Tin, T_out=Tout, p=pc, Q=Q, W_fan=Wf, m_cold=mf_air))
            Tc, pc = Tout, pc * (1 - ic["dp"])
    T2c, p2c = Tc, pc
    xl, xn, xr = spec["xl"], spec["xn"], spec["xr"]
    xn2, xr2 = spec.get("xn2", 0.0), spec.get("xr2", 0.0)
    mc = 1.0 - xl - xn - xr - xn2 - xr2
    if mc <= 0.3:
        return dict(ok=False, why="доля воздуха на горение меньше 0,3")
    rec = spec.get("rec")
    fuel = spec.get("fuel")
    rh = spec.get("reheat")
    dpa, dpg = (rec["dpa"], rec["dpg"]) if rec else (0.0, 0.0)
    dpf = fuel["dp_g"] if fuel else 0.0
    p3 = p2c * (1 - dpa)
    p_exit = E.P0 * p_exh * (1 + dpf) * (1 + dpg)
    T41 = spec["t41"] + 273.15
    T3, Tf = T2c, E.T0
    eta_gg, eta_pt = spec["eta_gg"], spec.get("eta_pt")
    res = None
    for it in range(60):
        hf = E.HF + (h_ch4(Tf) - h_ch4(E.T0))
        PRb = spec["PRb"](T3, p3, mc) if callable(spec["PRb"]) else spec["PRb"]
        p4 = p3 * PRb
        df1, T4, m41, b41 = burn_to(E, mc, T3, 0.0, T41, xn, T2c, hf)
        fuel_kg = df1 * mc
        sec = []
        if spec["shaft"] == "single":
            if rh:
                p_rh = p4 * (p_exit / p4) ** rh["x_rh"]
                if p_rh >= p4 * 0.999 or p_rh <= p_exit * 1.001:
                    return dict(ok=False, why="повторный подогрев: давление подогрева вне перепада турбины")
                Th, w1 = expand(E, T41, b41, p4 / p_rh, eta_gg)
                sec.append(dict(name="турбина до подогрева", T_in=T41, p_in=p4, T_out=Th, p_out=p_rh, m=m41, b=b41, w=w1))
                m_a, T_a, b_a = mix(E, [(m41, Th, b41), (xr, T2c, 0.0)])
                p_r2 = p_rh * (1 - rh["dp"])
                df2, T4r, m41r, b41r = burn_to(E, m_a, T_a, b_a, rh["t_rh"] + 273.15, xn2, T2c, hf)
                fuel_kg += df2 * m_a
                ER2 = p_r2 / p_exit
                T5, w2 = expand(E, rh["t_rh"] + 273.15, b41r, ER2, eta_gg)
                sec.append(dict(name="турбина после подогрева", T_in=rh["t_rh"] + 273.15, p_in=p_r2, T_out=T5, p_out=p_exit, m=m41r,
                                b=b41r, w=w2, T4=T4r))
                W_t = m41 * w1 + m41r * w2
                m_g, Tex, b_g = mix(E, [(m41r, T5, b41r), (xr2, T2c, 0.0)])
            else:
                ER = p4 / p_exit
                if ER <= 1.05:
                    return dict(ok=False, why="степень расширения турбины <= 1,05")
                T5, wt = expand(E, T41, b41, ER, eta_gg)
                sec.append(dict(name="турбина", T_in=T41, p_in=p4, T_out=T5, p_out=p_exit, m=m41, b=b41, w=wt))
                W_t = m41 * wt
                m_g, Tex, b_g = mix(E, [(m41, T5, b41), (xr, T2c, 0.0)])
            wnet = eta_m * W_t - wc - W_fan
        else:
            wgg = wc / (eta_m * m41)
            T45 = E.T_h(E.h(T41, b41) - wgg, b41)
            ERgg = math.exp((E.s0(T41, b41) - E.s0(T45, b41)) / (eta_gg * E.Rg(b41)))
            p45 = p4 / ERgg * (1 - spec.get("dp_duct", 0.01))
            sec.append(dict(name="турбина газогенератора", T_in=T41, p_in=p4, T_out=T45, p_out=p4 / ERgg, m=m41, b=b41, w=wgg))
            m45, T45m, b45 = mix(E, [(m41, T45, b41), (xr, T2c, 0.0)])
            if rh:
                p45r = p45 * (1 - rh["dp"])
                df2, T4r, m45, b45 = burn_to(E, m45, T45m, b45, rh["t_rh"] + 273.15, xn2, T2c, hf)
                fuel_kg += df2 * (m41 + xr)
                T45m, p45 = rh["t_rh"] + 273.15, p45r
            if p45 <= p_exit * 1.15:
                return dict(ok=False, why="силовой турбине не остаётся перепада (степень расширения < 1,15)")
            ERpt = p45 / p_exit
            T5, wpt = expand(E, T45m, b45, ERpt, eta_pt)
            sec.append(dict(name="силовая турбина", T_in=T45m, p_in=p45, T_out=T5, p_out=p_exit, m=m45, b=b45, w=wpt))
            W_t = m41 * wgg + m45 * wpt
            wnet = eta_m * m45 * wpt - W_fan
            m_g, Tex, b_g = mix(E, [(m45, T5, b45), (xr2, T2c, 0.0)])
        # ---------------------------------------------------------------- выхлоп: подогрев топлива -> рекуператор
        T_turb_exit = Tex
        Tg = Tex
        fuel_tot = fuel_kg / eta_b
        q_fuel = 0.0
        Tf_new = E.T0
        if fuel:
            Tf_new = max(min(Tex - fuel["dT_app"], fuel["T_max"] + 273.15), E.T0)
            q_fuel = fuel_tot * (h_ch4(Tf_new) - h_ch4(E.T0))
            Tg = E.T_h(E.h(Tex, b_g) - q_fuel / m_g, b_g)
        T_rec_in = Tg if rec else None
        T3_new, Q_rec, T_stack = T2c, 0.0, Tg
        if rec:
            if Tg < T2c + 10.0:
                return dict(ok=False, why="рекуператор бесполезен: газ на входе холоднее воздуха за компрессором")
            if rec.get("model") == "temperature":
                T3_new = T2c + rec["eps"] * (Tg - T2c)
                Q_rec = mc * (E.h(T3_new) - E.h(T2c))
            else:
                Q_max = min(mc * (E.h(Tg) - E.h(T2c)), m_g * (E.h(Tg, b_g) - E.h(T2c, b_g)))
                Q_rec = rec["eps"] * Q_max
                T3_new = E.T_h(E.h(T2c) + Q_rec / mc)
            T_stack = E.T_h(E.h(Tg, b_g) - Q_rec / m_g, b_g)
        done = abs(T3_new - T3) < 0.02 and abs(Tf_new - Tf) < 0.05 and it > 0
        T3, Tf = T3_new, Tf_new
        res = dict(T4=T4, p4=p4, PRb=PRb, m41=m41, b41=b41, sec=sec, W_t=W_t, wnet=wnet, fuel_tot=fuel_tot, m_g=m_g, b_g=b_g,
                   T_turb_exit=T_turb_exit, T_rec_in=T_rec_in, T_stack=T_stack, Q_rec=Q_rec, q_fuel=q_fuel, T_fuel=Tf_new)
        if done or (not rec and not fuel and not callable(spec["PRb"])):
            break
    else:
        return dict(ok=False, why="итерации рекуператора / подогрева топлива не сошлись")
    # ---------------------------------------------------------------- обращённый цикл Брайтона на выхлопе
    ibc_out = None
    wnet = res["wnet"]
    T_out_stack = res["T_stack"]
    if spec.get("ibc"):
        ibc_out = ibc_run(E, res["m_g"], res["T_stack"], res["b_g"], spec["ibc"], p_exh, eta_m)
        if ibc_out and ibc_out["w_net"] > 0:
            wnet += ibc_out["w_net"]
            T_out_stack = ibc_out["T_out"]
        else:
            ibc_out = dict(useless=True, **(ibc_out or {}))
    if wnet <= 0:
        return dict(ok=False, why="удельная работа <= 0")
    eff = wnet / (res["fuel_tot"] * E.LHV)
    return dict(ok=True, eff=eff, wnet=wnet, wc=wc, W_t=res["W_t"], W_fan=W_fan, Q_ic=Q_ic, m_fan_air=m_fan_air,
                fuel=res["fuel_tot"], T2=T2c, p2=p2c, T3=T3, p3=p3, T4=res["T4"], p4=res["p4"], PRb=res["PRb"], T41=T41,
                m41=res["m41"], b41=res["b41"], mc=mc, sections=res["sec"], m_exh=res["m_g"], b_exh=res["b_g"],
                T_turb_exit=res["T_turb_exit"], T_rec_in=res["T_rec_in"], Q_rec=res["Q_rec"], T_stack_rec=res["T_stack"],
                T_stack=T_out_stack, q_fuel=res["q_fuel"], T_fuel=res["T_fuel"], p_exit=p_exit, stages=st_out, ic=ic_out,
                ibc=ibc_out, xl=xl, xn=xn, xr=xr, xn2=xn2, xr2=xr2)


def ibc_run(E, m_g, T_in, b, p, p_exh, eta_m):
    """обращённый цикл Брайтона (как cycles.ibc_run, без конденсации): газ из-за двигателя (p_exh·P0) расширяется до
    plow·P0, охлаждается атмосферным воздухом (eps, вентилятор) и сжимается обратно до p_exh·P0. Лучший plow из сетки."""
    best = None
    p_in = E.P0 * p_exh
    for pl in p["plow"]:
        ER = p_in / (pl * E.P0)
        if ER <= 1.0:
            continue
        Tt, wt = expand(E, T_in, b, ER, p["eta_t"])
        Tc_out = Tt - p["eps"] * (Tt - E.T0)
        Q = E.h(Tt, b) - E.h(Tc_out, b)
        Wf, m_cold = fan_power(E, m_g * Q, Tt, p)
        pc_in = pl * E.P0 * (1 - p["dp"])
        pi_c = p_in / pc_in
        T2 = E.T_s(E.s0(Tc_out, b) + E.Rg(b) * math.log(pi_c) / p["eta_c"], b)
        wc = E.h(T2, b) - E.h(Tc_out, b)
        w = eta_m * m_g * wt - m_g * wc - Wf
        r = dict(plow=pl, w_net=w, W_t=m_g * wt, W_c=m_g * wc, W_fan=Wf, Q=m_g * Q, T_t_out=Tt, T_cool=Tc_out, T_out=T2,
                 pi_c=pi_c, ER_t=ER, m_cold=m_cold, w_c_per_kg=wc)
        if best is None or w > best["w_net"]:
            best = r
    return best
