#!/usr/bin/env python3
"""Оценка кандидата семейства «лопатки и колёса»: КПД ступени и запас помпажа (модель уровня 1 воронки + поправки
на приёмы), масса и момент инерции, напряжения (быстрая оценка; элиты — CalculiX, ccx_elites.py), число деталей,
признаки вида. Всё — оценка уровня облика.

Расчётная точка KV-3 — routes/gtu/frame3/frame3_point_results.json (stage1, stage2, turbine, cycle).
Модели — routes/gtu/funnel/meanline_l1.py (cc_stage, cc_off, ax_turbine — потери Oh/Yoon/Chung и Kacker — Okapuu),
пороги срыва — funnel/offdesign_l1.py, правила — funnel/funnel_l1.py, материалы — platform/fea/materials.json.
"""
from __future__ import annotations

import json
import math
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
GTU = os.path.dirname(os.path.dirname(HERE))
FUN = os.path.join(GTU, "funnel")
for _p in (HERE, FUN):
    if _p not in sys.path:
        sys.path.insert(0, _p)
import meanline_l1 as ML          # noqa: E402
import offdesign_l1 as OD          # noqa: E402
import funnel_l1 as FL             # noqa: E402
import blade_geom as BG            # noqa: E402
import core_local as C             # noqa: E402
import disk1d as D1                # noqa: E402

FRAME3 = os.path.join(GTU, "frame3", "frame3_point_results.json")
FRAME3_FEA = os.path.join(GTU, "frame3", "frame3_fea_results.json")
U2_LIMIT = 0.9 * 620.0            # frame3_point: U2 не выше 0,9·620 м/с (проектный запас)
K_ROUTE = 0.9                     # допускаемое печатного/литого Ti-6Al-4V к деформируемому (оценка, как organic k_AM)
BURST_MIN = 1.30                  # запас по разрушающей частоте при 100 % (оценка: на 8 % выше разгона 120 %; у frame3 —
                                  # колесо 2: 1,90, диск турбины 1: 1,40)
BLADE_BEND_MAX = 0.5              # изгиб лопатки колеса от навала при 120 % не выше 0,5·σ0,2 (оценка: запас на вибрацию)
LIFE_H = 25000.0                  # ресурс по длительной прочности, ч (как funnel/strength_l1.RUPTURE_H)
SM_CAP = 0.30                     # запас помпажа сверх 30 % в пригодность не засчитывается

# =========================================================================== приёмы: поправки и их статус
TECH = {
    "tubercles": dict(
        d_eta_per_Ac=-0.10, d_stall_deg_per_Ac=60.0, applies="колесо компрессора и ротор турбины",
        status="оценка по аналогии; перенос с крыла на решётку колеса — гипотеза",
        src=("Miklosovic, Murray, Howle, Fish (2004), Phys. Fluids 16(5) L39: модель плавника горбатого кита — угол срыва"
             " больше примерно на 40 %; Johari et al. (2007), AIAA J. 45(11): до срыва подъёмная сила чуть ниже, после — выше"
             " (по памяти, первоисточники в сессии не открывались)"),
        model="КПД −0,10·(A/c) (−0,5 п. при A/c 0,05); порог срыва по углу атаки на периферии входа +60°·(A/c) (+3° при 0,05)"),
    "sweep_le": dict(
        d_stall_deg_per_sw=5.0, status="гипотеза", src=("стреловидность вперёд расширяет запас устойчивости трансзвуковых"
                                                         " вентиляторов — Wadia et al. (1998), ASME J. Turbomach. 120 (по памяти)"),
        model="порог срыва на периферии входа +5°·sw при sw > 0 (периферия впереди), КПД 0"),
    "splitter_offset": dict(
        d_eta_max=0.003, opt=0.42, width=0.12, status="оценка",
        src="смещение разделительной лопатки к стороне разрежения основной — порядка +0,3…0,5 п. в расчётах CFD (по памяти)",
        model="КПД +0,003·(1 − ((x − 0,42)/0,12)²), не ниже −0,004; x — положение по шагу"),
    "splitter_recamber": dict(status="гипотеза, в цель не входит", src="выгиб входа разделительной под местный поток (по смыслу)",
                              model="0"),
    "lean_imp": dict(status="КПД — гипотеза (0); напряжения считаются", src="навал на выходе меняет скольжение и вторичные течения (по смыслу)",
                     model="КПД 0; изгиб σ = 3ρω²r2·b2²·sin γ / t (консоль под центробежной нагрузкой)"),
    "serrations": dict(d_eta_per=-0.0067, status="оценка", src=("серрации задней кромки снижают шум кромки (Howe 1991, J. Sound Vib.,"
                                                                " по памяти); аэродинамика колеса — потеря на смешение у зубцов, оценка"),
                       model="КПД −0,0067·(глубина/b2) (−0,2 п. при 0,3)"),
    "tandem": dict(d_eta_max=0.005, opt=0.15, width=0.15, d_stall_deg=1.0, status="гипотеза",
                   src="тандемные (с разрезом индуктора) колёса — в открытых работах около +1 п. и шире характеристика (по памяти, источник не открыт)",
                   model="КПД +0,005·(1 − ((сдвиг − 0,15)/0,15)²), не ниже −0,005; порог срыва +1°"),
    "riblets": dict(s_plus=15.0, gain=0.06, min_feature={"print": 150e-6, "cast": 300e-6}, status="расчёт применимости",
                    src="Walsh (1983), NASA: −4…8 % трения при s+ ≈ 10–20 (по памяти)",
                    model="нужный шаг s = s+·ν/u_τ сравнивается с наименьшим элементом печати/литья; не выходит — выигрыш 0"),
    "fishbone": dict(status="масса — расчёт; прочность — эквивалентная толщина (консервативно для окружных), рёбра не засчитаны",
                     src="рёбра под лопатками и карманы между ними на задней поверхности (по смыслу «рыбьей кости»); нужна 3D-проверка CalculiX",
                     model="h_eff = h − 0,5·глубина (среднее sin²)"),
    "turb_bow": dict(d_eta_max=0.004, opt=15.0, width=15.0, status="оценка",
                     src="Harrison (1990), ASME J. Turbomach. 112:184 — составной навал перераспределяет вторичные потери к середине (по памяти)",
                     model="КПД ступени +0,004·(1 − ((γ − 15°)/15°)²), не ниже −0,006; γ — угол у торцов"),
    "turb_vortex": dict(k=-0.003, status="оценка", src="отход от свободного вихря — радиальные течения (по смыслу)",
                        model="КПД −0,003·(n + 1)²; реактивность у втулки 1 − (ψ/2)(r_h/r_m)^(n−1) — расчёт"),
    "turb_hollow": dict(status="масса и напряжения — расчёт по сечениям; только печать; охлаждения нет (класс U)",
                        src="полость в профиле снимает центробежную нагрузку (по смыслу)", model="A·(1 − f), I·(1 − f²)"),
    "turb_lean_sweep": dict(status="КПД — гипотеза (0); изгиб от смещения центров — расчёт", src="по смыслу", model="0"),
}


def _clip(x, lo, hi):
    return max(lo, min(hi, x))


def _bump(x, opt, width, amp, floor):
    return max(floor, amp * (1 - ((x - opt) / width) ** 2))


# =========================================================================== расчётная точка
_PT = {}


def point():
    if not _PT:
        d = json.load(open(FRAME3, encoding="utf-8"))
        s1, s2 = d["stage1"], d["stage2"]
        cyc = d["cycle"]["frame3"]
        W = cyc["W_exh"]
        Wf = cyc["Wf"]
        _PT.update(
            stage1=dict(T01=s1["T01_C"] + 273.15, p01=s1["p01_kPa"] * 1e3, m_side=s1["m_per_side_kg_s"], PR=s1["PR"], n=s1["n_rpm"],
                        dbl=True, T02=s1["T02_C"] + 273.15, mu=s1["work_coeff_lambda"], kcm=s1["phi2_cm2_U2"]),
            stage2=dict(T01=s2["T01_C"] + 273.15, p01=s2["p01_kPa"] * 1e3, m_side=s2["m_per_side_kg_s"], PR=s2["PR"], n=s2["n_rpm"],
                        dbl=False, T02=s2["T02_C"] + 273.15, mu=s2["work_coeff_lambda"], kcm=s2["phi2_cm2_U2"]),
            turb=dict(T0in=1273.15, p0in=101325.0 * cyc["ER"], m=W, far=Wf / (W - Wf), w=d["turbine"]["work_needed_kJkg"] * 1e3,
                      nst=d["turbine"]["selected_z"], n=d["turbine"]["options"][str(d["turbine"]["selected_z"])]["n"]),
            src=dict(file="routes/gtu/frame3/frame3_point_results.json",
                     fields="stage1/stage2: T01_C, p01_kPa, m_per_side_kg_s, PR, n_rpm, T02_C, work_coeff_lambda, phi2_cm2_U2; "
                            "cycle.frame3: ER, W_exh, Wf; turbine: work_needed_kJkg, selected_z, options.n; T41 1000 °C"))
        fea = json.load(open(FRAME3_FEA, encoding="utf-8"))
        _PT["fea"] = fea["parts"]
    return _PT


# =========================================================================== компрессор


def surge_margin(st, T01, p01, m_d, n, stall_bonus=0.0, n_pts=24):
    """граница устойчивости на линии частоты одной ступени (как offdesign_l1.surge, но для ступени): идём от расчётного
    расхода вниз; граница — вершина π или первый признак срыва (угол атаки лопаток диффузора > DINC_STALL, на периферии
    входа > TIP_INC_STALL + бонус приёмов, угол потока за колесом > ALPHA2_STALL). Запас — (π_s/m_s)/(π_d/m_d) − 1"""
    d0 = ML.cc_off(st, T01, p01, m_d, n)
    if not d0.get("PRs"):
        return dict(ok=False, why="нет решения в расчётной точке")
    PRd = d0["PRs"]

    def stall(r):
        w = []
        if r["dinc_deg"] > OD.DINC_STALL:
            w.append("диффузор")
        if r["inc_tip_deg"] > OD.TIP_INC_STALL + stall_bonus:
            w.append("периферия входа")
        if r["alpha2_deg"] > OD.ALPHA2_STALL:
            w.append("α2")
        return w
    if stall(d0):
        return dict(ok=True, SM=0.0, why="срыв уже в расчётной точке", crit=",".join(stall(d0)), PRd=PRd, m_s=m_d, PR_s=PRd)
    flows = m_d * np.linspace(1.0, 0.40, n_pts + 1)[1:]
    m_prev, PR_prev = m_d, PRd
    crit = "край сетки"
    m_s, PR_s = flows[-1], None
    for m in flows:
        r = ML.cc_off(st, T01, p01, m, n)
        if not r.get("PRs"):
            m_s, PR_s, crit = m_prev, PR_prev, "нет решения"
            break
        w = stall(r)
        if w:
            lo, hi = m, m_prev                    # lo — срыв, hi — без срыва
            PR_hi = PR_prev
            for _ in range(6):
                mid = 0.5 * (lo + hi)
                rm = ML.cc_off(st, T01, p01, mid, n)
                if rm.get("PRs") and not stall(rm):
                    hi, PR_hi = mid, rm["PRs"]
                else:
                    lo = mid
            m_s, PR_s, crit = hi, PR_hi, "срыв: " + ",".join(w)
            break
        if r["PRs"] < PR_prev:
            m_s, PR_s, crit = m_prev, PR_prev, "вершина π"
            break
        m_prev, PR_prev = m, r["PRs"]
    if PR_s is None:
        PR_s = PR_prev
        m_s = m_prev
    SM = (PR_s / m_s) / (PRd / m_d) - 1
    return dict(ok=True, SM=float(SM), crit=crit, PRd=float(PRd), m_s=float(m_s), PR_s=float(PR_s), eta_off_d=d0["eta_p"])


def rules_comp(st):
    bad = []
    if st["M1t"] > FL.M1T_MAX:
        bad.append(("M1t", st["M1t"], FL.M1T_MAX))
    if st["W2_W1t"] < FL.DH_MIN:
        bad.append(("W2/W1t", st["W2_W1t"], FL.DH_MIN))
    if st["D1t_D2"] > FL.D1D2_MAX:
        bad.append(("D1t/D2", st["D1t_D2"], FL.D1D2_MAX))
    if st["alpha2_deg"] > FL.ALPHA2_MAX:
        bad.append(("alpha2", st["alpha2_deg"], FL.ALPHA2_MAX))
    if "forward" in st["flags"] or st["beta2b_deg"] > FL.BETA2_MAX:
        bad.append(("beta2b", st["beta2b_deg"], FL.BETA2_MAX))
    if st["U2"] > U2_LIMIT:
        bad.append(("U2", st["U2"], U2_LIMIT))
    return bad


def _viol(bad):
    v = 0.0
    for _, val, lim in bad:
        v += abs(val - lim) / max(abs(lim), 1e-9) + 1e-3
    return v


def riblet_check(st, route):
    """нужный шаг риблетов в колесе: u_τ = W·sqrt(Cf/2), Cf = 0,0412·Re^−0,1925 (как cc_losses), s = s+·ν/u_τ"""
    f, g = st["f"], st["g"]
    W = 0.5 * (f["W1t"] + f["W2"])
    rho = 0.5 * (f["rho1"] + f["rho2"])
    mu = 0.5 * (f["mu1"] + f["mu2"])
    Re = f["U2"] * 2 * g["r2"] * rho / mu
    Cf = 0.0412 * Re ** (-0.1925)
    ut = W * math.sqrt(Cf / 2)
    s = TECH["riblets"]["s_plus"] * (mu / rho) / ut
    feat = TECH["riblets"]["min_feature"][route]
    ok = s >= feat
    gain = TECH["riblets"]["gain"] * st["losses_pct"]["sf"] / 100 if ok else 0.0
    return dict(s_needed_um=s * 1e6, min_feature_um=feat * 1e6, feasible=bool(ok), d_eta=gain)


def comp_corrections(p, st, route):
    d = {}
    if p["tub_A"] > 0:
        d["tubercles"] = TECH["tubercles"]["d_eta_per_Ac"] * p["tub_A"]
    if p["lsp"] > 0:
        o = TECH["splitter_offset"]
        d["splitter_offset"] = _bump(p["sp_off"], o["opt"], o["width"], o["d_eta_max"], -0.004)
    if p["ser_A"] > 0:
        d["serrations"] = TECH["serrations"]["d_eta_per"] * p["ser_A"]
    if p["tand"] > 0:
        o = TECH["tandem"]
        d["tandem"] = _bump(p["tand_clk"], o["opt"], o["width"], o["d_eta_max"], -0.005)
    rib = riblet_check(st, route)
    if rib["d_eta"]:
        d["riblets"] = rib["d_eta"]
    bonus = TECH["tubercles"]["d_stall_deg_per_Ac"] * p["tub_A"] + TECH["sweep_le"]["d_stall_deg_per_sw"] * max(0.0, p["sweep"])
    if p["tand"] > 0:
        bonus += TECH["tandem"]["d_stall_deg"]
    return d, bonus, rib


def comp_temperature(stage, r, r2):
    """поле температур колеса как frame3_fea: от T01 у оси к T02 на ободе по квадрату радиуса (оценка)"""
    P = point()["stage%d" % stage]
    T01, T02 = P["T01"] - 273.15, P["T02"] - 273.15
    return T01 + (T02 - T01) * (np.asarray(r) / r2) ** 2


_REF = {}
_KP = []


_KT = {}


def _ktable(comp):
    if comp not in _KT:
        _KT[comp] = D1.k_table_cal(comp)
    return _KT[comp]


def baseline_genes(comp):
    """облик frame3 в генах семейства: μ и cm2/U2 точки, 9 + 9 лопаток, разделительная 0,70 посередине шага, ν 0,40,
    Безье 0,5/0,5, линейный закон углов, без приёмов; задняя поверхность ступени 2 — кованая (база диска)"""
    if comp == "turb":
        p = dict(phi=0.5, U=420.0, zn=40, zt=60, taper=0.6, hwall=0.0, tmax_c=0.20, vortex=-1.0, bow=0.0, lean=0.0, sweep=0.0,
                 tub_A=0.0, tub_n=4, w1=1.0, w2=1.0, w3=1.0, w4=1.0, fish=0.0, route=1)
        return BG.encode("turb", p)
    P = point()["stage1" if comp == "comp1" else "stage2"]
    p = dict(mu=P["mu"], kcm=P["kcm"], Z=9, lsp=0.70, sp_off=0.5, sp_dbeta=0.0, nu=0.40, hub_a=0.5, hub_c=0.5, sh_a=0.5, sh_c=0.5,
             law1=1 / 3, law2=2 / 3, lean=0.0, sweep=0.0, tub_A=0.0, tub_n=4, ser_A=0.0, ser_n=5, tand=0.0, tand_m=0.3, tand_clk=0.2,
             t_hub_k=1.0, t_tip_k=1.0, ffd1=0.0, ffd2=0.0, ffd3=0.0, ffd4=0.0, route=0)
    if comp == "comp2":
        p.update({f"zb{i}": 0.0 for i in range(6)})
        x = BG.encode(comp, dict(p, fish=0.0, fish_r0=0.4, fish_r1=0.8, boss=1.0))
        x = base_backface(x)
        return x
        p.update(fish=0.0, fish_r0=0.4, fish_r1=0.8)
    return BG.encode(comp, p)


def base_backface(x):
    """гены zb0..zb5, при которых задняя поверхность ступени 2 совпадает с кованой frame3 (подобием)"""
    import scipy.optimize as so  # noqa: F401
    p = BG.decode("comp2", x)
    P = point()["stage2"]
    st = ML.cc_stage(P["T01"], P["p01"], P["m_side"], P["PR"], P["n"], p["mu"], p["kcm"], p["Z"] * (1 + p["lsp"]) / (1 + ML.A["spl"]),
                     False, kD=1.9, nu=p["nu"])
    g = st["g"]
    r2 = g["r2"]
    fb = BG.FRAME3_BACK()
    r_knots = np.r_[0.0, 0.058, 0.07, 0.11, 0.15, 0.20, 0.25, 0.30, fb["r2"]] / fb["r2"] * r2
    base = np.interp(r_knots / r2 * fb["r2"], fb["r"], fb["depth"]) * r2 / fb["r2"]
    d_end = 0.02 * 2 * r2
    zb = []
    prev = base[1]
    for k_ in range(2, 8):
        q = (base[k_] - d_end) / max(prev - d_end, 1e-9)
        zb.append((min(1.0, q) - BG.Q_LO) / (1 - BG.Q_LO))
        prev = d_end + (prev - d_end) * (BG.Q_LO + (1 - BG.Q_LO) * min(1.0, max(0.0, zb[-1])))
    zb = np.asarray(zb)
    names = [gg.name for gg in BG.genes_for("comp2")]
    x = np.array(x, float)
    for i in range(6):
        x[names.index(f"zb{i}")] = float(np.clip(zb[i], 0, 1))
    return x


def reference(comp):
    if comp not in _REF:
        r = evaluate(comp, baseline_genes(comp), use_ref=False)
        _REF[comp] = dict(mass=r["meta"]["mass_kg"], J=r["meta"]["J_kgm2"], eta=r["meta"].get("eta"), raw=r)
    return _REF[comp]


def evaluate(comp, x, want_geom=False, use_ref=True):
    if comp in ("comp1", "comp2"):
        return eval_comp(comp, x, want_geom, use_ref)
    return eval_turb(x, want_geom, use_ref)


def _infeasible(why, viol, d, meta=None):
    m = dict(why=why, viol=viol)
    if meta:
        m.update(meta)
    return dict(f=-1000.0 - viol, d=d, feasible=False, meta=m)


def eval_comp(comp, x, want_geom=False, use_ref=True):
    stage = 1 if comp == "comp1" else 2
    P = point()["stage%d" % stage]
    p = BG.decode(comp, x)
    route = ("print", "cast")[int(p["route"])]
    Z = int(p["Z"])
    Zeff = Z * (1 + p["lsp"])
    d_nominal = np.array([45.0, 60.0, p["lsp"], p["tub_A"], p["lean"]])
    try:
        st = ML.cc_stage(P["T01"], P["p01"], P["m_side"], P["PR"], P["n"], p["mu"], p["kcm"], Zeff / (1 + ML.A["spl"]), P["dbl"],
                         kD=1.9, nu=p["nu"], d_pen=0.01 if P["dbl"] else 0.0)
    except (ValueError, ZeroDivisionError, OverflowError) as e:
        return _infeasible(f"средняя линия: {e}", 10.0, d_nominal)
    if not st.get("ok"):
        return _infeasible(st.get("why"), 10.0, d_nominal)
    st["g"]["Zfull_used"] = Z
    bad = rules_comp(st)
    geo = BG.impeller(st, p, stage, build_mesh=want_geom)
    d = np.array([st["beta2b_deg"], geo["wrap_shroud_deg"], p["lsp"], p["tub_A"], p["lean"]])
    corr, bonus, rib = comp_corrections(p, st, route)
    sm = surge_margin(st, P["T01"], P["p01"], P["m_side"], P["n"], bonus)
    eta = st["eta_p"] + sum(corr.values())
    # прочность диска (быстро) при 100 и 120 %
    dk = geo["disk"]
    r = dk["r"]
    hb = geo["hb"]
    fac = 1 + geo["blade_fill_phi"] * hb / np.maximum(dk["h_eff"], 1e-6)
    T = comp_temperature(stage, r, r[-1])
    kt = _ktable(comp)
    if kt is not None:
        s100 = D1.solve(r, dk["h_eff"], fac, T, "ti64", P["n"], k_table=kt)
        s120 = D1.solve(r, dk["h_eff"], fac, T, "ti64", 1.2 * P["n"], k_table=kt)
    else:                     # нет калибровки по CalculiX — единая поправка пика (консервативно)
        kp = D1.k_peak_cal()
        s100 = D1.solve(r, dk["h_eff"], fac, T, "ti64", P["n"], k_peak=kp)
        s120 = D1.solve(r, dk["h_eff"], fac, T, "ti64", 1.2 * P["n"], k_peak=kp)
    u_disk = max(1.5 * s100["yield_ratio_max"], s120["yield_ratio_max"]) / K_ROUTE
    # изгиб лопатки от навала на выходе при 120 %
    mat = D1.material("ti64")
    om = 1.2 * P["n"] * math.pi / 30
    t_tip = 0.0045 * st["D2"] * p["t_tip_k"]
    sig_b = 3 * mat["density"] * om * om * st["D2"] / 2 * st["b2"] ** 2 * abs(math.sin(math.radians(p["lean"]))) / t_tip
    ys_T = float(np.interp(P["T02"] - 273.15, *np.asarray(mat["ys_MPa"], float).T)) * 1e6
    u_blade = sig_b / (BLADE_BEND_MAX * ys_T * K_ROUTE)
    checks = dict(geo["checks"])
    cons = list(bad)
    if not checks["valid"]:
        cons.append(("геометрия", 1.0, 0.0))
    if u_disk > 1:
        cons.append(("загрузка диска", u_disk, 1.0))
    if s100["burst_margin"] < BURST_MIN:
        cons.append(("запас по разрушению", s100["burst_margin"], BURST_MIN))
    if u_blade > 1:
        cons.append(("изгиб лопатки", u_blade, 1.0))
    if not sm.get("ok"):
        cons.append(("помпаж", 0.0, 0.01))
    SM = sm.get("SM", 0.0)
    mass, J = geo["mass_kg"], geo["J_kgm2"]
    ref = reference(comp) if use_ref else dict(mass=mass, J=J)
    feats = geo["features"]
    nov = novelty_vector(comp, p, d)
    meta = dict(comp=comp, eta=eta, eta_level1=st["eta_p"], eta_uncapped=st["eta_uncapped"], capped=st["capped"], corrections=corr,
                SM=SM, surge=sm, mass_kg=mass, J_kgm2=J, mass_blades_kg=geo["mass_blades_kg"], mass_disk_kg=geo["mass_disk_kg"],
                u_disk=u_disk, yr100=s100["yield_ratio_max"], yr120=s120["yield_ratio_max"], burst=s100["burst_margin"],
                sig_bend_MPa=sig_b / 1e6, u_blade=u_blade, parts=1, route=route, D2=st["D2"], b2=st["b2"], U2=st["U2"],
                beta2b=st["beta2b_deg"], Zeff=Zeff, M1t=st["M1t"], riblets=rib, checks=_jsonable(checks), features=feats,
                novelty_vec=nov.tolist(), rules_bad=[b[0] for b in bad], blade_fill_phi=geo["blade_fill_phi"])
    if cons:
        out = _infeasible("; ".join(c[0] for c in cons), _viol(cons), d, meta)
    else:
        f = 100 * eta + 50 * min(SM, SM_CAP) - 20 * (J / ref["J"] - 1) - 10 * (mass / ref["mass"] - 1)
        out = dict(f=float(f), d=d, feasible=True, meta=meta)
    if want_geom:
        out["geom"] = geo
        out["ml"] = st
    return out


def novelty_vector(comp, p, d):
    """нормированный вектор признаков формы для новизны (0..1 по каждому)"""
    if comp == "turb":
        return np.array([(d[0] - 50) / 20, d[1] / 80, d[2] / 45, (d[3] + 15) / 30, d[4] / 0.05, 1.0 if p.get("hwall", 0) > 0 else 0.0, (p["lean"] + 15) / 30,
                         (p["taper"] - 0.25) / 0.65, p["fish"] / 0.5])
    return np.array([(d[0] - 20) / 40, _clip(d[1], 0, 150) / 150, d[2] / 0.85, d[3] / 0.06, (d[4] + 25) / 50, (p["sweep"] + 0.2) / 0.4,
                     p["ser_A"] / 0.3, 1.0 if p["tand"] > 0 else 0.0, p.get("fish", 0.0) / 0.6, abs(p["sp_off"] - 0.5) / 0.2 if p["lsp"] else 0.0,
                     math.sqrt(sum(p[f"ffd{i}"] ** 2 for i in range(1, 5))) / 2])


def oem_reference_vectors(comp):
    """типовые формы OEM (по смыслу общеизвестных колёс и лопаток, не по конкретному изделию): колесо — загнутость 40–50°,
    разделительные 0,7 посередине, без навала, бугорков, серраций, тандема и карманов; турбина — радиальная укладка,
    свободный вихрь, без навала и бугорков, сплошная лопатка"""
    out = []
    if comp == "turb":
        for b3 in (60.0, 63.4, 66.0):
            p = dict(hwall=0.0, lean=0.0, taper=0.6, fish=0.0)
            out.append(novelty_vector("turb", p, np.array([b3, 25.0, 0.0, 0.0, 0.0])))
        return np.array(out)
    for b2b in (40.0, 45.0, 50.0):
        for wrap in (50.0, 70.0):
            p = dict(sweep=0.0, ser_A=0.0, tand=0.0, fish=0.0, sp_off=0.5, lsp=0.7, ffd1=0, ffd2=0, ffd3=0, ffd4=0)
            out.append(novelty_vector(comp, p, np.array([b2b, wrap, 0.7, 0.0, 0.0])))
    return np.array(out)


def _jsonable(o):
    if isinstance(o, dict):
        return {k: _jsonable(v) for k, v in o.items()}
    if isinstance(o, (np.bool_,)):
        return bool(o)
    if isinstance(o, (np.floating,)):
        return float(o)
    if isinstance(o, (np.integer,)):
        return int(o)
    return o


# =========================================================================== турбина


def turb_disk(k, bl, zt, p, rpm, n=161):
    """блиск ступени k: профиль диска frame3 (turb{k}, равнопрочный Стодолы, ступица по осевому шагу) подобием по радиусу
    до корня лопатки; ёлочного замка и «мёртвой» зоны у цельного колеса нет — обод (полоса w_rim) доходит до корня;
    толщина полотна × PCHIP(w1..w4); «рыбья кость» снимает долю полотна (эквивалентная толщина). Поле T — frame3
    (ступень 1 — cooled_750, ступени 2, 3 — T_rel), подобием по радиусу"""
    from scipy.interpolate import PchipInterpolator
    F = point()["fea"]["turb%d" % k]
    geo = F["geometry"]
    pr = F["profile"]
    var = F["variants"]["cooled_750" if k == 1 else "T_rel"]
    r0 = np.asarray(pr["r"], float)
    t0 = np.asarray(pr["t"], float)
    r_live0 = geo["r_live_rim_m"]
    r_hub0 = geo["r_hub_m"]
    w_rim = geo["w_rim_m"]
    s = bl["rh"] / r_hub0
    r_live = r_live0 * s
    r = np.linspace(0.0, bl["rh"], n)
    h_web = np.interp(r / s, r0, t0)
    wk = PchipInterpolator([0.0, 0.35, 0.7, 1.0], [p["w1"], p["w2"], p["w3"], p["w4"]])(np.clip(r / r_live, 0, 1))
    h = np.where(r <= r_live, h_web * wk, w_rim)
    h = np.maximum(h, 0.006)
    fish = np.zeros_like(r)
    if p["fish"] > 0:
        w = BG.smoothstep((r / r_live - 0.35) / 0.08) * BG.smoothstep((0.85 - r / r_live) / 0.08)
        fish = p["fish"] * np.maximum(h - 0.008, 0) * w
    h_eff = h - 0.5 * fish
    T = np.interp(r / s, np.asarray(var["T_C"]["r"]), np.asarray(var["T_C"]["T"]))
    p_rim = zt * bl["Fc"] / (2 * math.pi * bl["rh"] * h[-1])
    res = D1.solve(r, h_eff, 1.0, T, "in713c", rpm, p_rim_Pa=p_rim)
    res120 = D1.solve(r, h_eff, 1.0, T, "in713c", 1.2 * rpm, p_rim_Pa=p_rim * 1.44)
    mass, J = D1.mass_inertia(r, h_eff, D1.material("in713c")["density"])
    uts120 = res120["burst_margin"]
    return dict(r=r, h=h, h_eff=h_eff, fish=fish, T=T, p_rim_MPa=p_rim / 1e6, mass=mass, J=J, s100=res, s120=res120,
                burst=res["burst_margin"], life=res.get("creep_life_min_h", float("inf")), w_rim=w_rim, r_live=r_live,
                burst120=uts120)


def eval_turb(x, want_geom=False, use_ref=True):
    P = point()["turb"]
    p = BG.decode("turb", x)
    route = ("print", "cast")[int(p["route"])]
    E = ML._E()
    b = E.bfrac(P["far"])
    d_nom = np.array([math.degrees(math.atan(1 / p["phi"])), 0.0, math.degrees(math.atan(4 * p["bow"])), p["sweep"], p["tub_A"]])
    try:
        tm = ML.ax_turbine(P["T0in"], P["p0in"], b, P["m"], p["U"], P["n"], P["w"], P["nst"], p["phi"], int(p["zn"]), int(p["zt"]))
    except (ValueError, ZeroDivisionError, OverflowError) as e:
        return _infeasible(f"средняя линия турбины: {e}", 10.0, d_nom)
    cons = []
    corr = {}
    gam_b = math.degrees(math.atan(4 * p["bow"]))
    o = TECH["turb_bow"]
    corr["turb_bow"] = _bump(gam_b, o["opt"], o["width"], o["d_eta_max"], -0.006)
    corr["turb_vortex"] = TECH["turb_vortex"]["k"] * (p["vortex"] + 1) ** 2
    if p["tub_A"] > 0:
        corr["tubercles"] = TECH["tubercles"]["d_eta_per_Ac"] * p["tub_A"]
    eta = tm["eta_p"] + sum(corr.values())
    stages, blades, disks = [], [], []
    rpm = P["n"]
    lmp = D1.material("in713c")["_lmp"]
    uts_tab = np.asarray(D1.material("in713c")["uts_MPa"], float)
    nv = p["vortex"]
    for s in tm["stages"]:
        k = s["stage"]
        rm = s["Dm"] / 2
        Hr = 0.5 * (s["H2"] + s["H3"])
        rh = rm - Hr / 2
        R_hub = 1 - (s["psi"] / 2) * (rh / rm) ** (nv - 1)
        ca = s["phi"] * tm["U"]
        c2h = math.hypot(ca, s["psi"] * tm["U"] * (rh / rm) ** nv)
        M2h = s["M2_hub"] * c2h / math.hypot(ca, s["psi"] * tm["U"] * rm / max(rm - s["H2"] / 2, 1e-6))
        for nm_, val, lim, up in (("M2", s["M2"], FL.M2T_MAX, True), ("Mw3", s["Mw3"], FL.MW3_MAX, True), ("R", s["reaction"], FL.R_MIN, False),
                                  ("M2h", M2h, ML.HUB_M2_MAX, True), ("Rh", R_hub, ML.R_HUB_MIN, False),
                                  ("Mw3t", s["Mw3_tip"], FL.MW3_MAX, True)):
            if (up and val > lim) or (not up and val < lim):
                cons.append((f"{nm_}{k}", val, lim if lim != 0 else 0.01))
        st = dict(s, m_dot=P["m"], zt=int(p["zt"]))
        bl = BG.turbine_rotor_blade(st, p, rpm)
        T_met = s["T_rel_in"] - 273.15
        sig = bl["sig_c_MPa"] + max(0.0, bl["sig_b_peak_MPa"])
        life, ex = lmp.life_h(max(sig, 1.0), T_met)
        uts = float(np.interp(T_met, uts_tab[:, 0], uts_tab[:, 1]))
        if life < LIFE_H:
            cons.append((f"ресурс лопатки {k}", life, LIFE_H))
        if 1.2 ** 2 * sig > 0.8 * uts:
            cons.append((f"σ лопатки {k} при 120 %", 1.44 * sig, 0.8 * uts))
        if not bl["checks"]["valid"]:
            cons.append((f"геометрия лопатки {k}", 1.0, 0.0))
        dk = turb_disk(k, bl, int(p["zt"]), p, rpm)
        if dk["life"] < LIFE_H:
            cons.append((f"ресурс диска {k}", dk["life"], LIFE_H))
        if dk["burst"] < BURST_MIN:
            cons.append((f"разрушение диска {k}", dk["burst"], BURST_MIN))
        blades.append(bl)
        disks.append(dk)
        stages.append(dict(stage=k, R_hub=R_hub, M2_hub=M2h, sig_c_MPa=bl["sig_c_MPa"], sig_b_MPa=bl["sig_b_peak_MPa"],
                           sig_b_gas_MPa=bl["sig_b_gas_MPa"], T_metal_C=T_met, life_h=life, life_extrap=ex, uts_MPa=uts,
                           disk_life_h=dk["life"], disk_burst=dk["burst"], disk_p_rim_MPa=dk["p_rim_MPa"], m_blade_kg=bl["mass_blade"],
                           m_disk_kg=dk["mass"], H=bl["H"], rh=bl["rh"], chord_root=float(bl["chord"][0]), twist_in_deg=bl["features"]["twist_in_deg"]))
    zt = int(p["zt"])
    mass = sum(dk["mass"] + zt * bl["mass_blade"] for dk, bl in zip(disks, blades))
    J = sum(dk["J"] + zt * bl["J_blade"] for dk, bl in zip(disks, blades))
    d = np.array([math.degrees(math.atan(1 / p["phi"])), abs(blades[-1]["features"]["twist_in_deg"]), gam_b, p["sweep"], p["tub_A"]])
    ref = reference("turb") if use_ref else dict(mass=mass, J=J)
    meta = dict(comp="turb", eta=eta, eta_level1=tm["eta_p"], corrections=corr, mass_kg=mass, J_kgm2=J, parts=3,
                parts_inserted_equiv=3 * (1 + zt), route=route, stages=stages, Dm=tm["Dm"], length_m=tm["length"],
                Dtip=tm["Dtip"], hollow=float(np.mean([bl_["hollow"] for bl_ in blades])), novelty_vec=novelty_vector("turb", p, d).tolist(),
                features=blades[-1]["features"], SM=None)
    if cons:
        out = _infeasible("; ".join(c[0] for c in cons), _viol(cons), d, meta)
    else:
        f = 100 * eta - 20 * (J / ref["J"] - 1) - 10 * (mass / ref["mass"] - 1)
        out = dict(f=float(f), d=d, feasible=True, meta=meta)
    if want_geom:
        out["geom"] = dict(blades=blades, disks=disks, tm=tm)
    return out


# =========================================================================== дескрипторы и пределы архива
DESC = {
    "comp1": dict(names=["угол лопатки на выходе β2b, град", "закрутка (wrap) на периферии, град", "доля разделительных",
                         "амплитуда бугорков A/c", "навал на выходе, град"],
                  lo=[20.0, 0.0, 0.0, 0.0, -25.0], hi=[60.0, 150.0, 0.85, 0.06, 25.0]),
    "turb": dict(names=["угол выхода ротора β3, град", "закрутка входа ротора ст. 3, град", "составной навал у торцов, град",
                        "стреловидность, град", "амплитуда бугорков A/c"],
                 lo=[50.0, 0.0, 0.0, -15.0, 0.0], hi=[70.0, 80.0, 45.0, 15.0, 0.05]),
}
DESC["comp2"] = DESC["comp1"]
