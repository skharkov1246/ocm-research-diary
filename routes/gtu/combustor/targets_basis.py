"""Расчётное обоснование целей камеры сгорания двигателя-заменителя KV-3 (Taurus 60 GS, мягкий допуск).

Уровень — 0D (Cantera): сеть реакторов «реактор идеального смешения (зона пламени) + реактор вытеснения (зона
догорания)», адиабатическая температура, ламинарная скорость пламени, задержка самовоспламенения, срыв в реакторе
идеального смешения, разнотопливность (метан, метан — водород, н-додекан как заменитель дизеля), чувствительность
цикла к потерям давления в камере. Всё — оценка уровня облика; геометрии камеры здесь нет.

Вход (файл и поле — в results-JSON, раздел inputs):
  routes/gtu/frame3/frame3_point_results.json — расчётная точка KV-3: PR, T41, G, КПД узлов, T3, Wf;
  routes/gtu/frame3/frame3_point.py — X_LEAK, DP_COMB, P_EXH (как в расчёте точки);
  routes/gtu/gtu_cycle.py — cycle(), AIR, T0, P0, LHV;
  routes/gtu/gtu_explore.py — LIM['T_ax_max'] (предел газа для неохлаждаемой осевой лопатки).

Запуск:  python3 targets_basis.py            # около 5–10 мин на 4 ядрах -> targets_basis_results.json
AWS не вызывается.
"""
import json
import math
import os
import sys
import time
from multiprocessing import Pool

import cantera as ct
import numpy as np
from scipy.optimize import brentq

HERE = os.path.dirname(os.path.abspath(__file__))
GTU = os.path.dirname(HERE)
sys.path.insert(0, GTU)
import gtu_cycle as C  # noqa: E402  (импорт меняет текущую папку на routes/gtu)
import gtu_explore as E  # noqa: E402

os.chdir(HERE)
OUT = os.path.join(HERE, "targets_basis_results.json")
F3 = json.load(open(os.path.join(GTU, "frame3", "frame3_point_results.json")))

X_LEAK, DP_COMB, P_EXH = 0.015, 0.04, 1.02          # как frame3_point.py (X_LEAK, DP_COMB, P_EXH)
AIR = C.AIR
AIR_REITZ = "o2:0.2095,n2:0.7905"                    # в механизме Reitz нет Ar и CO2: Ar отнесён к N2
DZ_MECH = "nDodecane_Reitz.yaml::nDodecane_IG"       # идеальный газ (в файле есть и фаза Redlich-Kwong)
O2_REF = 15.0
O2_AIR = 20.9                                        # как в навыке openfoam-combustor, раздел 7
MW_NO2, V_M = 46.0055, 22.414                        # г/моль; л/моль при 0 °C и 101,325 кПа
TAU_PSR = 0.5e-3                                     # время пребывания зоны пламени, с (оценка, см. assumptions)
TAU_HOT = [3e-3, 5e-3, 10e-3, 20e-3, 30e-3]          # полное время в горячей зоне до разбавления, с
PHI_GRID = [0.45, 0.50, 0.55, 0.60, 0.65, 0.70]
UNMIX_S = [0.0, 0.05, 0.10, 0.15, 0.20]              # sigma(phi)/phi — неравномерность смешения
H2_VOL = [0.0, 0.3, 0.5, 0.7, 1.0]                   # доля водорода в топливе, объёмн.
T_FLAME_REF = 1800.0                                 # K — сравнение топлив при одной адиабатической температуре
DP_GRID = [0.02, 0.03, 0.04, 0.05, 0.06]
T41_PART = [1000.0, 950.0, 900.0, 850.0, 800.0, 750.0, 700.0, 650.0, 600.0]


# ---------------------------------------------------------------- состояние смеси
def _sol(mech):
    """'файл.yaml' или 'файл.yaml::фаза'."""
    f, _, ph = mech.partition("::")
    return ct.Solution(f, ph) if ph else ct.Solution(f)


def mixed(phi, fuel, p, T_air, T_fuel, mech="gri30.yaml", oxid=AIR):
    """Адиабатическое смешение воздуха (T_air) и топлива (T_fuel) при давлении p для коэффициента phi."""
    g = _sol(mech)
    g.set_equivalence_ratio(phi, fuel, oxid)
    fuel_sp = [s.split(":")[0] for s in fuel.split(",")]
    Z = sum(g.Y[g.species_index(s)] for s in fuel_sp)            # доля топливного потока по массе
    Y = g.Y.copy()
    ga = _sol(mech); ga.TPX = T_air, p, oxid
    gf = _sol(mech); gf.TPX = T_fuel, p, fuel
    h = Z * gf.enthalpy_mass + (1 - Z) * ga.enthalpy_mass
    g.HPY = h, p, Y
    return g, Z


def corr15(g):
    """Сухие ppm при 15 % O2: NOx = NO + NO2, CO; плюс O2 сухой, %."""
    X = dict(zip(g.species_names, g.X))
    xw = X.get("H2O", 0.0)
    dry = lambda x: x / (1 - xw)
    o2 = 100 * dry(X.get("O2", 0.0))
    k = (O2_AIR - O2_REF) / (O2_AIR - o2)
    nox = 1e6 * dry(X.get("NO", 0.0) + X.get("NO2", 0.0)) * k
    co = 1e6 * dry(X.get("CO", 0.0)) * k
    return nox, co, o2


def psr(gin, tau, mech="gri30.yaml"):
    """Реактор идеального смешения с временем пребывания tau, старт с равновесия (горящая ветвь)."""
    inlet = ct.Reservoir(gin, clone=False)
    g = ct.Solution(mech); g.TPY = gin.T, gin.P, gin.Y; g.equilibrate("HP")
    r = ct.IdealGasReactor(g, energy="on", clone=False)
    ex = ct.Reservoir(g, clone=False)
    mfc = ct.MassFlowController(inlet, r, mdot=lambda t: r.mass / tau)
    ct.PressureController(r, ex, primary=mfc, K=0.01)
    net = ct.ReactorNet([r])
    net.rtol, net.atol = 1e-9, 1e-15
    net.advance_to_steady_state()
    out = ct.Solution(mech); out.TPY = r.T, r.phase.P, r.phase.Y
    return out


def pfr(gin, times, mech="gri30.yaml"):
    """Реактор вытеснения как лагранжев объём при постоянном давлении; состояние в моменты times (с)."""
    g = ct.Solution(mech); g.TPY = gin.T, gin.P, gin.Y
    r = ct.IdealGasConstPressureReactor(g, energy="on", clone=False)
    net = ct.ReactorNet([r]); net.rtol, net.atol = 1e-9, 1e-15
    res = []
    for t in times:
        net.advance(t)
        s = ct.Solution(mech); s.TPY = r.T, r.phase.P, r.phase.Y
        res.append(s)
    return res


def crn(phi, fuel, p, T_air, T_fuel, tau_psr=TAU_PSR, tau_hot=TAU_HOT):
    gin, Z = mixed(phi, fuel, p, T_air, T_fuel)
    T_in = gin.T
    geq = ct.Solution("gri30.yaml"); geq.TPY = gin.T, gin.P, gin.Y; geq.equilibrate("HP")
    T_ad = geq.T
    gp = psr(gin, tau_psr)
    lit = gp.T > T_in + 300
    row = dict(phi=phi, T_mix_K=T_in, T_ad_K=T_ad, lit=bool(lit), T_psr_K=gp.T)
    if not lit:
        return row
    nox, co, o2 = corr15(gp)
    row.update(psr=dict(NOx_ppm15=nox, CO_ppm15=co))
    states = pfr(gp, [t - tau_psr for t in tau_hot])
    row["hot"] = []
    for t, s in zip(tau_hot, states):
        n, c, o = corr15(s)
        row["hot"].append(dict(tau_ms=t * 1e3, T_K=s.T, NOx_ppm15=n, CO_ppm15=c, O2_dry_pct=o))
    nq, cq, _ = corr15(geq)
    row["eq"] = dict(NOx_ppm15=nq, CO_ppm15=cq)
    return row


def psr_lbo_phi(fuel, p, T_air, T_fuel, tau):
    """Наименьший phi, при котором реактор идеального смешения с временем tau ещё горит (бисекция)."""
    def lit(phi):
        gin, _ = mixed(phi, fuel, p, T_air, T_fuel)
        return psr(gin, tau).T > gin.T + 300
    lo, hi = 0.15, 0.70
    if not lit(hi):
        return None
    for _ in range(14):
        mid = 0.5 * (lo + hi)
        if lit(mid):
            hi = mid
        else:
            lo = mid
    gin, _ = mixed(hi, fuel, p, T_air, T_fuel)
    geq = ct.Solution("gri30.yaml"); geq.TPY = gin.T, gin.P, gin.Y; geq.equilibrate("HP")
    return dict(tau_ms=tau * 1e3, phi_lbo=hi, T_ad_at_lbo_K=geq.T)


def phi_for_Tad(fuel, p, T_air, T_fuel, T_target, mech="gri30.yaml", oxid=AIR):
    def f(phi):
        g, _ = mixed(phi, fuel, p, T_air, T_fuel, mech, oxid)
        g.equilibrate("HP")
        return g.T - T_target
    return brentq(f, 0.2, 1.0, xtol=1e-5)


def flame_speed(fuel, phi, p, T_air, T_fuel):
    g, _ = mixed(phi, fuel, p, T_air, T_fuel)
    f = ct.FreeFlame(g, width=0.01)
    f.set_refine_criteria(ratio=3, slope=0.07, curve=0.14)
    f.transport_model = "mixture-averaged"
    f.solve(loglevel=0, auto=True)
    return float(f.velocity[0]), float(f.T[-1]), len(f.grid)


def ign_delay(phi, fuel, p, T_air, T_fuel, mech="gri30.yaml", oxid=AIR, t_end=1.0):
    """Задержка самовоспламенения в реакторе постоянного давления: время максимума dT/dt; None — нет за t_end."""
    g, _ = mixed(phi, fuel, p, T_air, T_fuel, mech, oxid)
    T0_ = g.T
    r = ct.IdealGasConstPressureReactor(g, energy="on", clone=False)
    net = ct.ReactorNet([r]); net.rtol, net.atol = 1e-9, 1e-15
    t_prev, T_prev, best = 0.0, r.T, (0.0, 0.0)
    while net.time < t_end:
        net.step()
        dTdt = (r.T - T_prev) / max(net.time - t_prev, 1e-15)
        if dTdt > best[0]:
            best = (dTdt, net.time)
        t_prev, T_prev = net.time, r.T
        if r.T > T0_ + 400 and dTdt < 0.01 * best[0]:
            break
    return (best[1] if r.T > T0_ + 400 else None), T0_


def far_for_T(T_target, p, T_air, T_fuel, fuel="CH4:1"):
    """Общее отношение топливо/воздух (по массе) для равновесной температуры смеси T_target."""
    def f(phi):
        g, _ = mixed(phi, fuel, p, T_air, T_fuel)
        g.equilibrate("HP")
        return g.T - T_target
    phi = brentq(f, 0.05, 0.9, xtol=1e-6)
    g, Z = mixed(phi, fuel, p, T_air, T_fuel)
    return phi, Z / (1 - Z)


# ---------------------------------------------------------------- задачи для пула
def job(args):
    kind = args[0]
    if kind == "crn":
        _, phi, fuel, p, Ta, Tf, tp, th = args
        return ("crn", (phi, fuel, tp), crn(phi, fuel, p, Ta, Tf, tp, th))
    if kind == "lbo":
        _, fuel, p, Ta, Tf, tau = args
        return ("lbo", (fuel, tau), psr_lbo_phi(fuel, p, Ta, Tf, tau))
    if kind == "sl":
        _, fuel, phi, p, Ta, Tf = args
        return ("sl", fuel, flame_speed(fuel, phi, p, Ta, Tf))
    raise ValueError(kind)


def main():
    t_start = time.time()
    cyc = F3["cycle"]["frame3"]
    PR, T41_C, G = F3["split"]["pi"], 1000.0, 25.55
    for i in F3["inputs"]:
        if i["field"].startswith("slots[model=Taurus 60") and i["field"].endswith(".мягкий"):
            PR, T41_C, G = i["used"]["PR"], float(i["used"]["T41_C"]), i["used"]["G"]
    eta_pc, eta_t = cyc["eta_pc"], cyc["eta_t"]

    # 1. Повтор расчётной точки KV-3 и вход камеры
    base = C.cycle(PR, T41_C, eta_pc=eta_pc, eta_pt_gg=eta_t, eta_pt_pt=eta_t, x_ngv=0.0, x_rot=0.0, x_leak=X_LEAK,
                   dp_duct=0.0, W1=G, dp_comb=DP_COMB, p_exh_ratio=P_EXH)
    P_kW = base["P_MW"] * 1e3
    assert abs(P_kW - cyc["P_shaft_kW"]) / cyc["P_shaft_kW"] < 1e-3, (P_kW, cyc["P_shaft_kW"])
    p3 = C.P0 * PR
    T3 = base["T3_C"] + 273.15
    Wair = G * (1 - X_LEAK)
    Wf = base["Wf"]
    far = Wf / Wair
    g_st = ct.Solution("gri30.yaml"); g_st.set_equivalence_ratio(1.0, "CH4:1", AIR)
    far_st = g_st.Y[g_st.species_index("CH4")] / (1 - g_st.Y[g_st.species_index("CH4")])
    phi_overall = far / far_st
    T4_C = base["T4_C"]
    inlet = dict(p3_kPa=p3 / 1e3, T3_K=T3, T3_C=T3 - 273.15, W_air_comb_kg_s=Wair, W_fuel_kg_s=Wf, FAR=far,
                 FAR_stoich=far_st, phi_overall=phi_overall, T4_C=T4_C, T41_C=T41_C, P_shaft_kW=P_kW,
                 eff_shaft=base["eff"], dp_comb=DP_COMB, heat_input_MW=Wf * C.LHV / 1e6,
                 note="T4 = T41: класс U, охлаждения соплового аппарата нет (x_ngv = 0)")

    # 2. Чувствительность цикла к потерям давления в камере
    dp_rows = []
    for dp in DP_GRID:
        r = C.cycle(PR, T41_C, eta_pc=eta_pc, eta_pt_gg=eta_t, eta_pt_pt=eta_t, x_ngv=0.0, x_rot=0.0,
                    x_leak=X_LEAK, dp_duct=0.0, W1=G, dp_comb=dp, p_exh_ratio=P_EXH)
        dp_rows.append(dict(dp_comb=dp, P_shaft_kW=r["P_MW"] * 1e3, eff_shaft=r["eff"], EGT_C=r["EGT_C"]))
    d4 = next(x for x in dp_rows if x["dp_comb"] == 0.04)
    for x in dp_rows:
        x["dP_kW_vs_4pct"] = x["P_shaft_kW"] - d4["P_shaft_kW"]
        x["d_eff_pt_vs_4pct"] = 100 * (x["eff_shaft"] - d4["eff_shaft"])
    slope_kW = (dp_rows[-1]["P_shaft_kW"] - dp_rows[0]["P_shaft_kW"]) / (100 * (DP_GRID[-1] - DP_GRID[0]))
    slope_eff = 100 * (dp_rows[-1]["eff_shaft"] - dp_rows[0]["eff_shaft"]) / (100 * (DP_GRID[-1] - DP_GRID[0]))

    # 3. Зона пламени: адиабатическая температура от phi
    tad = []
    for phi in [0.40, 0.45, 0.50, 0.55, 0.60, 0.65, 0.70, 0.80, 1.00]:
        g, _ = mixed(phi, "CH4:1", p3, T3, C.T0); Tm = g.T; g.equilibrate("HP")
        tad.append(dict(phi=phi, T_mix_K=Tm, T_ad_K=g.T, premix_air_share=phi_overall / phi))
    phi_1700 = phi_for_Tad("CH4:1", p3, T3, C.T0, 1700.0)
    phi_1800 = phi_for_Tad("CH4:1", p3, T3, C.T0, 1800.0)
    phi_1900 = phi_for_Tad("CH4:1", p3, T3, C.T0, 1900.0)

    # 4. Неполная нагрузка при одном клапане и неизменном разделении воздуха
    #    (оценка: постоянная частота вала, постоянный расход воздуха, T3 постоянна, p3 ~ sqrt(T41) — сопло заперто)
    part = []
    for T41p in T41_PART:
        k = math.sqrt((T41p + 273.15) / (T41_C + 273.15))
        PRp, p3p = PR * k, p3 * k
        r = C.cycle(PRp, T41p, eta_pc=eta_pc, eta_pt_gg=eta_t, eta_pt_pt=eta_t, x_ngv=0.0, x_rot=0.0,
                    x_leak=X_LEAK, dp_duct=0.0, W1=G, dp_comb=DP_COMB, p_exh_ratio=P_EXH)
        part.append(dict(T41_C=T41p, PR=PRp, p3_kPa=p3p / 1e3, T3_C=r["T3_C"], P_shaft_kW=r["P_MW"] * 1e3,
                         load_frac=r["P_MW"] * 1e3 / P_kW, phi_overall=(r["Wf"] / Wair) / far_st))

    # 5. Параллельные расчёты: сеть реакторов, срыв, скорости пламени
    T_air, T_fuel = T3, C.T0
    tasks = []
    for phi in PHI_GRID:
        tasks.append(("crn", phi, "CH4:1", p3, T_air, T_fuel, TAU_PSR, TAU_HOT))
    for tp in (0.2e-3, 1.0e-3):
        tasks.append(("crn", 0.55, "CH4:1", p3, T_air, T_fuel, tp, TAU_HOT))
    # неравномерность смешения: дискретная нормальная плотность phi (11 точек ±3 sigma) при phi_mean = 0.55
    PHI_M = 0.55
    z_nodes = np.linspace(-3, 3, 13)
    w_nodes = np.exp(-0.5 * z_nodes ** 2); w_nodes /= w_nodes.sum()
    unmix_phis = sorted({round(PHI_M * (1 + s * z), 4) for s in UNMIX_S[1:] for z in z_nodes} | {PHI_M})
    for phi in unmix_phis:
        if phi not in PHI_GRID:
            tasks.append(("crn", phi, "CH4:1", p3, T_air, T_fuel, TAU_PSR, TAU_HOT))
    for tau in (0.2e-3, 0.5e-3, 1.0e-3):
        tasks.append(("lbo", "CH4:1", p3, T_air, T_fuel, tau))
    # водород: phi при T_ad = 1800 K
    h2_cases = []
    for xh in H2_VOL:
        fuel = "CH4:1" if xh == 0 else ("H2:1" if xh == 1 else f"CH4:{1 - xh:.2f},H2:{xh:.2f}")
        phi = phi_for_Tad(fuel, p3, T_air, T_fuel, T_FLAME_REF)
        h2_cases.append((xh, fuel, phi))
        tasks.append(("crn", round(phi, 5), fuel, p3, T_air, T_fuel, TAU_PSR, TAU_HOT))
        tasks.append(("sl", fuel, phi, p3, T_air, T_fuel))
    with Pool(min(4, os.cpu_count() or 1)) as pool:
        out = pool.map(job, tasks, chunksize=1)
    crn_map = {k: v for kind, k, v in out if kind == "crn"}
    lbo = [v for kind, k, v in out if kind == "lbo"]
    sl = {k: v for kind, k, v in out if kind == "sl"}

    grid = [crn_map[(phi, "CH4:1", TAU_PSR)] for phi in PHI_GRID]
    tau_sens = [crn_map[(0.55, "CH4:1", tp)] | dict(tau_psr_ms=tp * 1e3) for tp in (0.2e-3, TAU_PSR, 1.0e-3)]

    # неравномерность: струйки независимы (без догорания между ними), каждая приведена к 15 % O2,
    # среднее с весом phi (≈ доля топлива, оценка)
    unmix = []
    for sig in UNMIX_S:
        res_t = []
        for it, tau in enumerate(TAU_HOT):
            num_n, num_c, den, unlit = 0.0, 0.0, 0.0, 0.0
            for z, w in zip(z_nodes, w_nodes):
                phi = round(PHI_M * (1 + sig * z), 4) if sig > 0 else PHI_M
                rr = crn_map[(phi, "CH4:1", TAU_PSR)]
                if not rr["lit"]:
                    unlit += w
                    continue
                s = rr["hot"][it]
                # к 15 % O2 каждая струйка приводится отдельно; среднее по массе топлива ~ phi (оценка)
                wt = w * phi
                num_n += wt * s["NOx_ppm15"]; num_c += wt * s["CO_ppm15"]; den += wt
            res_t.append(dict(tau_ms=tau * 1e3, NOx_ppm15=num_n / den, CO_ppm15=num_c / den,
                              unlit_weight=unlit))
        unmix.append(dict(S=sig, phi_mean=PHI_M, by_tau=res_t))

    h2 = []
    for xh, fuel, phi in h2_cases:
        rr = crn_map[(round(phi, 5), fuel, TAU_PSR)]
        S_L, Tb, npts = sl[fuel]
        tig, Tm = ign_delay(phi, fuel, p3, T_air, T_fuel, t_end=1.0)
        h2.append(dict(H2_vol=xh, fuel=fuel, phi=phi, T_ad_K=rr["T_ad_K"], S_L_m_s=S_L, flame_points=npts,
                       ign_delay_ms=None if tig is None else tig * 1e3, T_mix_K=Tm,
                       NOx_ppm15_10ms=next(x["NOx_ppm15"] for x in rr["hot"] if abs(x["tau_ms"] - 10) < 1e-9),
                       CO_ppm15_10ms=next(x["CO_ppm15"] for x in rr["hot"] if abs(x["tau_ms"] - 10) < 1e-9)))

    # 6. Неполная нагрузка: phi зоны пламени при неизменной доле воздуха на смешение (настройка на 1800 K при 100 %)
    share = phi_overall / phi_1800
    lbo05 = next(x for x in lbo if abs(x["tau_ms"] - 0.5) < 1e-9)
    for row in part:
        phi_pz = row["phi_overall"] / share
        g, _ = mixed(phi_pz, "CH4:1", row["p3_kPa"] * 1e3, row["T3_C"] + 273.15, C.T0); g.equilibrate("HP")
        row.update(phi_pz_fixed_split=phi_pz, T_ad_pz_K=g.T, above_psr_lbo_0p5ms=bool(phi_pz > lbo05["phi_lbo"]),
                   premix_share_needed_for_1800K=row["phi_overall"] / phi_for_Tad(
                       "CH4:1", row["p3_kPa"] * 1e3, row["T3_C"] + 273.15, C.T0, 1800.0))

    # 7. Дизель (н-додекан, механизм Reitz без азотной химии): T_ad и задержка самовоспламенения при T3, p3
    dz = []
    for phi in (0.5, 0.6, 1.0):
        g, _ = mixed(phi, "c12h26:1", p3, T3, T3, DZ_MECH, AIR_REITZ); Tm = g.T
        g.equilibrate("HP")
        tig, _ = ign_delay(phi, "c12h26:1", p3, T3, T3, DZ_MECH, AIR_REITZ, t_end=0.2)
        dz.append(dict(phi=phi, T_mix_K=Tm, T_ad_K=g.T, ign_delay_ms=None if tig is None else tig * 1e3))
    gch4, _ = mixed(1.0, "CH4:1", p3, T3, C.T0); gch4.equilibrate("HP")
    tig_ch4, _ = ign_delay(0.55, "CH4:1", p3, T3, C.T0, t_end=1.0)

    # 8. Объём горячей зоны для выбранных времён пребывания (масса через зону пламени при phi = phi_1800)
    gz, Zz = mixed(phi_1800, "CH4:1", p3, T3, C.T0); gz.equilibrate("HP")
    m_hot = Wair * share * (1 + far_st * phi_1800)
    vol = [dict(tau_ms=t * 1e3, V_hot_L=1e3 * m_hot * t / gz.density) for t in TAU_HOT]

    # 9. Допуск неравномерности на выходе: предел газа неохлаждаемой осевой лопатки
    T_ax_max = E.LIM["T_ax_max"]
    pf_max = (T_ax_max - T41_C) / (T41_C - (T3 - 273.15))

    res = dict(
        title="Камера сгорания KV-3: расчётное обоснование целей (0D, Cantera)",
        method=("0D: повтор расчётной точки KV-3 через gtu_cycle.cycle; адиабатическое равновесие (HP); сеть реакторов "
                "«реактор идеального смешения (зона пламени, tau_psr) + лагранжев реактор вытеснения при постоянном "
                "давлении (догорание до tau_hot)», GRI-Mech 3.0, NOx = NO + NO2, приведение к 15 % O2 сухого газа "
                "C15 = C·(20,9−15)/(20,9−O2); срыв — бисекция phi до погасания реактора идеального смешения; "
                "ламинарная скорость — FreeFlame (mixture-averaged); задержка самовоспламенения — максимум dT/dt в "
                "реакторе постоянного давления; неравномерность смешения — нормальная плотность phi (13 узлов ±3σ), "
                "струйки независимы, среднее с весом phi; дизель — н-додекан, механизм Reitz (без азота)"),
        software=f"python {sys.version.split()[0]}, cantera {ct.__version__}, numpy {np.__version__}, "
                 f"scipy {__import__('scipy').__version__}; механизмы gri30.yaml (53/325), nDodecane_Reitz.yaml (100/553)",
        date="2026-11-14",
        status="оценка уровня облика (0D): геометрии, течения и неравномерности по реальной камере нет; числа — "
               "ориентир для целей и для воронки, не прогноз выбросов конкретной конструкции",
        runtime_s=None,
        inputs=[
            dict(file="routes/gtu/frame3/frame3_point_results.json", field="inputs[field=slots[model=Taurus 60, variant=генератор].мягкий].used (PR, T41_C, G)",
                 used=dict(PR=PR, T41_C=T41_C, G=G)),
            dict(file="routes/gtu/frame3/frame3_point_results.json", field="cycle.frame3 (eta_pc, eta_t, P_shaft_kW, T3_C, Wf)",
                 used=dict(eta_pc=eta_pc, eta_t=eta_t, P_shaft_kW=cyc["P_shaft_kW"], T3_C=cyc["T3_C"], Wf=cyc["Wf"])),
            dict(file="routes/gtu/frame3/frame3_point.py", field="X_LEAK, DP_COMB, P_EXH",
                 used=dict(X_LEAK=X_LEAK, DP_COMB=DP_COMB, P_EXH=P_EXH)),
            dict(file="routes/gtu/gtu_cycle.py", field="cycle(), AIR, T0, P0, LHV",
                 used=dict(AIR=AIR, T0=C.T0, P0=C.P0, LHV=C.LHV)),
            dict(file="routes/gtu/gtu_explore.py", field="LIM['T_ax_max']", used=dict(T_ax_max=T_ax_max)),
        ],
        assumptions=dict(
            tau_psr_ms=TAU_PSR * 1e3,
            tau_psr_note="время зоны пламени 0,5 мс — оценка; чувствительность 0,2 и 1,0 мс в tau_psr_sensitivity",
            tau_hot_ms=[t * 1e3 for t in TAU_HOT],
            fuel_T_K=C.T0, fuel="метан (как gtu_cycle); водородные смеси — CH4/H2 по объёму",
            part_load=("постоянная частота вала, постоянный расход воздуха, T3 как в цикле при меньшем PR, "
                       "PR и p3 пропорциональны sqrt(T41) (запертое сопло), КПД узлов как в расчётной точке — оценка"),
            unmixedness=("S = sigma(phi)/phi_mean; без догорания между струйками — оценка сверху для NOx; "
                         "погасшие в своём реакторе бедные струйки исключены из среднего, их вес — unlit_weight "
                         "(в камере это недожог: CO и UHC)"),
            ign_delay=("null — самовоспламенения нет за t_end (газ — 1 с, дизель — 0,2 с); время — максимум dT/dt"),
            flame_speed="ламинарная S_L; турбулентное ускорение пламени водорода (числа Льюиса < 1) не учтено",
            diesel="н-додекан в газовой фазе при T3 (испарение и его охлаждение не учтены); азотной химии в механизме "
                   "нет — NOx дизеля здесь не считается",
            o2_ref=dict(O2_air=O2_AIR, O2_ref=O2_REF),
        ),
        inlet=inlet,
        cycle_check=dict(P_shaft_kW_repro=P_kW, P_shaft_kW_frame3=cyc["P_shaft_kW"],
                         rel_diff=(P_kW - cyc["P_shaft_kW"]) / cyc["P_shaft_kW"]),
        dp_sensitivity=dict(rows=dp_rows, dP_kW_per_pct=slope_kW, d_eff_pt_per_pct=slope_eff),
        flame_temperature=dict(rows=tad, phi_for_Tad_1700K=phi_1700, phi_for_Tad_1800K=phi_1800,
                               phi_for_Tad_1900K=phi_1900, premix_share_at_1800K=share),
        crn_grid=grid,
        tau_psr_sensitivity=tau_sens,
        psr_lbo=lbo,
        unmixedness=unmix,
        hydrogen=h2,
        part_load_fixed_split=part,
        diesel=dict(rows=dz, T_ad_stoich_CH4_K=gch4.T,
                    ign_delay_CH4_phi055_ms=None if tig_ch4 is None else tig_ch4 * 1e3,
                    T_air_K=T3, p_kPa=p3 / 1e3),
        hot_zone_volume=dict(phi=phi_1800, m_hot_kg_s=m_hot, rho_hot_kg_m3=gz.density, T_hot_K=gz.T, rows=vol),
        exit_uniformity=dict(T_ax_max_C=T_ax_max, T41_C=T41_C, T3_C=T3 - 273.15, pattern_factor_max=pf_max,
                             note="PF = (Tmax − T41)/(T41 − T3): предел местного максимума на входе в сопловой "
                                  "аппарат, чтобы газ не превысил предел неохлаждаемой лопатки (оценка)"),
        units=dict(mg_Nm3_per_ppm_NO2=MW_NO2 / V_M,
                   note="1 ppm NOx (как NO2) при 15 % O2 = 2,053 мг/нм3 при 0 °C, 101,325 кПа, сухой газ"),
    )
    res["runtime_s"] = round(time.time() - t_start, 1)

    def rnd(x):
        if isinstance(x, float):
            return float(f"{x:.6g}")
        if isinstance(x, dict):
            return {k: rnd(v) for k, v in x.items()}
        if isinstance(x, (list, tuple)):
            return [rnd(v) for v in x]
        return x
    json.dump(rnd(res), open(OUT, "w"), ensure_ascii=False, indent=1)
    print("written", OUT, res["runtime_s"], "s")


if __name__ == "__main__":
    main()
