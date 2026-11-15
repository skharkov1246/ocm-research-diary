"""Воронка двигателя-заменителя, уровень 1: углублённый расчёт отбора уровня 0 по нишам - часть отбора (на AWS).

Для каждого варианта отбора уровня 0 (merge_l0: funnel_l0_select*.json; концепция восстанавливается funnel_space.concept по
idx, n, kv, seed - отпечаток выборки и кода уровня 0 сверяется):
  1. оптимизация непрерывных осей внутри дискретной схемы варианта: несколько стартов, старт = дифференциальная эволюция (scipy,
     ограничение по правилам допустимости Лампинена, в первой популяции - точка уровня 0, зерно = f(SEED, idx, kv, старт)) +
     доводка COBYLA; цель - взвешенная оценка уровня 1, ограничение - сумма превышений правил = 0; оси и границы - var_spec;
     прирост давления RDC/волнового ротора не оптимизируется (это свойство технологии, а не выбор конструктора);
  2. расчётная точка уровня 1: цикл funnel_space._core с КПД из одномерных моделей ступеней (meanline_l1: центробежная -
     Oh/Yoon/Chung, осевая турбина - Soderberg + зазор Kacker - Okapuu, радиальная - вид Baines), геометрия funnel_space.derive,
     рекуператор по эффективности (размер и потери давления, recup_l1), быстрая прочность (strength_l1), правила уровня 0 и 1;
  3. для лучших точек стартов - нерасчётные режимы −30…+45 °C и 50-100 % нагрузки, запас по помпажу, пуск (offdesign_l1);
     первая точка, прошедшая все правила, - итог; иначе - наименее нарушающая, с причинами;
  4. баллы: детали, Plug&Play, монтаж, изготовимость - модели уровня 0 без изменений (пересчитаны для итоговой геометрии);
     КПД и размер - по уровню 1; взвешенная оценка по весам владельца; геометрия крупных узлов для эскизов (geometry).
Запуск части (AWS: несколько частей одним пулом на все ядра; результат не зависит от числа процессов):
  python funnel_l1.py --input l0_aws_full/funnel_l0_select_KV-3_n200000.json --top 4000 --shards 40 --shard 0,6,12 \
      --starts 6 --workers "$NPROC" --out "$GTU_OUT"
Слияние, отбор 50 разных и масштабирование на KV-1/KV-6 - merge_l1.py. Всё - оценка (README узла, раздел «Уровень 1»).
"""
import argparse, datetime, gzip, hashlib, json, math, os, sys, time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
GTU = os.path.dirname(HERE)
if HERE not in sys.path:
    sys.path.insert(0, HERE)
import funnel_space as FS        # noqa: E402
import funnel_l0 as L0           # noqa: E402
import make_model as MK          # noqa: E402
import parts_model as PM         # noqa: E402
import plugplay_model as PP      # noqa: E402
import service_model as SM       # noqa: E402
import meanline_l1 as ML         # noqa: E402
import strength_l1 as ST         # noqa: E402
import offdesign_l1 as OD        # noqa: E402

DATE = FS.DATE
SEED = FS.SEED
LEVEL = "L1"
PART_PREFIX = "funnel_l1"
DEFAULT_STARTS = 4               # прогонов дифференциальной эволюции на вариант
DEFAULT_MAXFEV = 200             # вычислений цели на доводку COBYLA после каждого прогона
DE_POP = 8                       # популяция DE = DE_POP × число осей
DE_GENS = 25                     # поколений DE на старт (проба 14.11.2026: COBYLA из случайных точек находила допустимую точку у
#                                  1 из 12 вариантов - сходится к недопустимым локальным минимумам; DE с правилами Лампинена - у 3 из 3
#                                  проверенных, оценка 6,5-7,3; README, раздел «Уровень 1»)
PENALTY = 20.0                   # штраф за превышение правил в цели COBYLA (в единицах оценки 0-10) - вместе с ограничением
OD_CANDIDATES = 4                # точек на проверку нерасчётных режимов (по убыванию оценки)
LOCAL_MAX_VARIANTS, LOCAL_MAX_STARTS = 12, 2      # решение владельца 14.11.2026: локально - не больше 12 вариантов и 2 стартов

# =========================================================================== правила уровня 1 (оценка)
M1T_MAX = 1.50          # относительное число Маха на периферии входа колеса (трансзвуковые колёса π 6-8 - около 1,3-1,5, по памяти;
#                         эталон frame3 - 1,27)
DH_MIN = 0.53           # W2/W1t колеса: W1t/W2 <= 1,9 (оценка по памяти о рекомендации Japikse; frame3 - 0,73)
D1D2_MAX = 0.80         # D1t/D2 (оценка)
ALPHA2_MAX = 78.0       # угол потока на выходе колеса от радиального, град (оценка; срыв в диффузоре)
BETA2_MAX = 60.0        # загнутость лопаток назад, град (оценка: технологичность)
M2T_MAX = 1.20          # число Маха за сопловым аппаратом осевой турбины (frame3_point.A m2_max)
MW3_MAX = 1.10          # относительное число Маха на выходе ротора осевой турбины (оценка)
R_MIN = 0.0             # реактивность ступени осевой турбины на среднем радиусе (frame3_point: reaction_ok)
M4_MAX = 1.10           # число Маха за сопловым аппаратом радиальной турбины (оценка)
MW5_MAX = 0.95          # относительное число Маха на периферии выхода радиального колеса (оценка)
T_TI_MAX = 427.0        # Ti-6Al-4V: предел данных materials.json, °C

# рекуператор (оценка): печатные каналы, ламинарное течение
REC = dict(DH=1.5e-3, NU=6.0, FRE=20.0, POROSITY=0.8, T_WALL=1.5e-4, K_WALL=15.0,
           basis="оценка: печатные каналы D_h 1,5 мм, Nu 6 (между квадратным каналом 3,0 и щелью 7,5 с интенсификацией), f·Re 20 "
                 "(Фаннинг), открытая доля объёма 0,8, стенка IN625 0,15 мм; потери давления - бюджет цикла (воздух 2,5 %, газ 3,5 % - "
                 "gtu_cycle.recup); коллекторы и корпус ×3 по объёму и массе (REC_HDR уровня 0)")


# =========================================================================== непрерывные оси
def var_spec(c):
    """оси оптимизации (имя, нижняя, верхняя, вид) в схеме концепции c"""
    st = FS.stage_types(c)
    tl = MK.t41_limit(c["turb"], c["cool"], c["mat_tw"], c["m_tw"])[0]
    lo_t = min(FS.T41_MIN, tl - 50.0)
    V = [("t41", lo_t, tl, "lin"), ("pi",) + tuple(FS.pi_range(c)) + ("log",)]
    if c["comp"] in FS.SPLIT_RANGE:
        V.append(("split",) + tuple(FS.SPLIT_RANGE[c["comp"]]) + ("lin",))
    V.append(("mu1",) + tuple(FS.MU_RANGE[st[0]]) + ("lin",))
    if len(st) > 1:
        V.append(("mu2",) + tuple(FS.MU_RANGE[st[-1]]) + ("lin",))
    if c["rotors"] == 2:
        k0 = next(t for t in st if t != "О")
        V.append(("phi1",) + tuple(FS.PHI_RANGE[k0]) + ("log",))
    V.append(("lt",) + ((0.9, FS.PSI_MAX) if c["turb"] == "О" else FS.NU_RANGE) + ("lin",))
    V += [("kcm", 0.20, 0.35, "lin"), ("zc", 7, 15, "int")]
    V += [("zt", 9, 17, "int"), ("zn", 12, 24, "int")] if c["turb"] == "Р" else [("zt", 30, 80, "int"), ("zn", 20, 60, "int")]
    if c["turb"] == "О":
        V.append(("phit", 0.45, 0.80, "lin"))
    if c["rotors"] == 2:
        V += [("kUpt", 0.85, 1.0, "lin"), ("phipt", 0.45, 0.80, "lin")]
    if c["recup"] == "печатный":
        V.append(("eps_rec", 0.75, 0.92, "lin"))
    V.append(("kD_case", 1.7, 2.1, "lin"))
    # число ступеней турбины газогенератора и силовой - среди разрешённых архетипом (как design() уровня 0): с проверкой у втулки
    # (свободный вихрь) двум ступеням часто не хватает реактивности и числа Маха - выбор 2/3 ступеней - дело оптимизатора, цена -
    # детали (исправление 15.11.2026; ниша от числа ступеней не зависит)
    nts = sorted(FS._axis_vals(c, c["arch"], "nt"))
    if len(nts) > 1:
        V.append(("nt", nts[0], nts[-1], "int"))
    if c["rotors"] == 2:
        npts = sorted(FS._axis_vals(c, c["arch"], "npt"))
        if len(npts) > 1:
            V.append(("npt", npts[0], npts[-1], "int"))
    return V


def x_to_c(c0, x, spec):
    c = dict(c0)
    for (name, lo, hi, kind), xi in zip(spec, x):
        u = min(1.0, max(0.0, float(xi)))
        if kind == "log":
            v = math.exp(math.log(lo) + (math.log(hi) - math.log(lo)) * u)
        else:
            v = lo + (hi - lo) * u
        c[name] = int(round(v)) if kind == "int" else v
    return c


def c_to_x(c, spec):
    x = []
    for name, lo, hi, kind in spec:
        v = float(c.get(name) or lo)
        if kind == "log":
            u = (math.log(max(v, 1e-12)) - math.log(lo)) / (math.log(hi) - math.log(lo))
        else:
            u = (v - lo) / (hi - lo) if hi > lo else 0.5
        x.append(min(1.0, max(0.0, u)))
    return np.array(x)


# =========================================================================== расчётная точка уровня 1
def _rule(code, ok, val=None, lim=None, text=""):
    return FS._rule(code, ok, val, lim, text)


def stage_pen(typ, mat, k, inlet):
    """поправки КПД ступени как уровень 0 (funnel_space.eta_comp, оценка): диагональная −0,01, колесо Ni −0,02, двусторонний вход
    1-й ступени −0,01 (до 15.11.2026 на уровне 1 двусторонний вход не учитывался, а Ni и «Д» вычитались только для цикла)"""
    return (0.01 if typ == "Д" else 0.0) + (0.02 if mat == "Ni" else 0.0) + (0.01 if (k == 0 and inlet == "двух") else 0.0)


def solve_l1(c, tgt, iters=6):
    """цикл с КПД ступеней из одномерных моделей: КПД -> цикл -> частота и размеры -> модели ступеней -> КПД (до сходимости)"""
    E = FS._E()
    typ, pis, mats, mus = FS.stage_types(c), FS.stages_pi(c), FS._stage_mats(c), FS.stage_mus(c)
    eta_c = [FS.eta_comp(t, p, m, k, c["inlet"], mu) for k, (p, t, m, mu) in enumerate(zip(pis, typ, mats, mus))]
    eg, ep = FS._eta_turb_nom(c)
    eg -= FS.ETA_NOM["ax"] if c["turb"] == "О" else FS.ETA_NOM["rad"]
    Upt = c["kUpt"] * FS.U_PT_MAX
    ccs, tm, pm, r, n, comp = [], None, None, None, None, None
    for it in range(iters):
        comp = FS._comp(list(zip(pis, eta_c)))
        r = FS._core(c, comp, eg, ep)
        if "fail" in r:
            return dict(ok=False, fail=r["fail"])
        G = tgt["P_shaft_kW"] * 1e3 / r["wnet"]
        Gc = G / c["cores"]
        stn = []
        for k, (s, t, m, mu) in enumerate(zip(comp, typ, mats, mus)):
            Gs = Gc / 2 if (k == 0 and c["inlet"] == "двух") else Gc
            stn.append(dict(type=t, Gs=Gs, rho1=s["p1"] / (E.RA * s["T1"]), U=math.sqrt(s["w"] / mu), mu=mu, mat=m))
        if c["rotors"] == 1:
            n = tgt["n_pkg"]
        else:
            s0 = next(s for s in stn if s["type"] != "О")
            n = 60 * s0["U"] / (math.pi * math.sqrt(s0["Gs"] / (s0["rho1"] * s0["U"] * c["phi1"])))
        new_c, ccs = [], []
        for k, (s, sn) in enumerate(zip(comp, stn)):
            if sn["type"] == "О":
                new_c.append(FS.eta_comp("О", s["pi"], sn["mat"], k, c["inlet"], sn["mu"]))
                ccs.append(None)
                continue
            cc = ML.cc_stage(s["T1"], s["p1"], sn["Gs"], s["pi"], n, sn["mu"], c["kcm"], float(c["zc"]),
                             k == 0 and c["inlet"] == "двух", c["kD_case"], eta0=eta_c[k], d_pen=stage_pen(sn["type"], sn["mat"], k, c["inlet"]),
                             typ=sn["type"])
            if not cc["ok"]:
                return dict(ok=False, fail=f"ступень компрессора {k + 1}: {cc['why']}")
            new_c.append(cc["eta_p"])                     # поправки ступени - внутри cc_stage (тот же КПД у геометрии и нерасчётной модели)
            ccs.append(cc)
        t = r["t"]
        Ut, Ulim, _ = FS.turb_speed(c, t)
        mg = t["m"] * Gc
        if c["turb"] == "О":
            tm = ML.ax_turbine(t["T_in"], t["p_in"], t["b"], mg, Ut, n, t["w"], c["nt"], c["phit"], c["zn"], c["zt"])
        else:
            tm = ML.rad_turbine(t["T_in"], t["p_in"], t["b"], mg, Ut, t["ER"], n, c["nt"], c["zt"], c["zn"])
        eg_new = tm["eta_p"]
        ep_new = ep
        if r["pt"]:
            p = r["pt"]
            pm = ML.ax_turbine(p["T_in"], p["p_in"], p["b"], p["m"] * G, Upt, tgt["n_pkg"], p["w"], c["npt"], c["phipt"], c["zn"], c["zt"])
            ep_new = pm["eta_p"]
        done = max([abs(a - b) for a, b in zip(new_c, eta_c)] + [abs(eg_new - eg), abs(ep_new - ep)]) < 2e-4
        eta_c, eg, ep = new_c, max(0.4, eg_new), max(0.4, ep_new)
        if done:
            break
    comp = FS._comp(list(zip(pis, eta_c)))
    r = FS._core(c, comp, eg, ep)
    if "fail" in r:
        return dict(ok=False, fail=r["fail"])
    G = tgt["P_shaft_kW"] * 1e3 / r["wnet"]
    t = r["t"]
    Ut, Ulim, _ = FS.turb_speed(c, t)
    if c["turb"] == "О":
        load = dict(psi=t["w"] / (c["nt"] * Ut ** 2))
    else:
        load = dict(nu=tm["nu"], dhs=tm["dhs"], eta_is=tm["eta_is"])
    if r["pt"]:
        load["psi_pt"] = r["pt"]["w"] / (c["npt"] * Upt ** 2)
    perf = dict(ok=True, comp=comp, r=r, eta_gg=eg, eta_pt=ep if r["pt"] else None, load=load, G=G, Ut=Ut, Ulim_t=Ulim, Upt=Upt,
                eta_th=r["eta"], w=r["wnet"] / 1e3, EGT=r["Tstack"] - 273.15, T5=r["T5"] - 273.15, T2=r["T2"] - 273.15)
    return dict(ok=True, perf=perf, ccs=ccs, tm=tm, pm=pm, n=n, eta_c=eta_c)


def recup_l1(c, perf):
    """рекуператор по эффективности: NTU противотока -> UA -> площадь (Nu ламинарных каналов) -> проходные сечения из бюджета
    потерь давления (f·Re) -> длина и объём (оценка; REC)"""
    E = FS._E()
    r = perf["r"]
    G = perf["G"]
    mc = 1.0 - r["xn"] - r["xr"] - r["xl"]
    ma = G * mc
    mg = G * (r["t"]["m"] + r["xr"]) if c["rotors"] == 1 else G * r["pt"]["m"]
    b = r["t"]["b"]
    Ta = 0.5 * (r["T2"] + r["T3"]); Tg = r["T5"] - 0.5 * (r["T3"] - r["T2"])
    Ca, Cg = ma * ML.cp(Ta), mg * ML.cp(Tg, b)
    Cmin, Cmax = min(Ca, Cg), max(Ca, Cg)
    Cr = Cmin / Cmax
    eps = c["eps_rec"]
    ntu = (math.log((1 - eps * Cr) / (1 - eps)) / (1 - Cr)) if Cr < 0.999 else eps / (1 - eps)
    UA = ntu * Cmin
    ha = REC["NU"] * ML.k_gas(Ta) / REC["DH"]; hg = REC["NU"] * ML.k_gas(Tg) / REC["DH"]
    U = 1 / (1 / ha + 1 / hg + REC["T_WALL"] / REC["K_WALL"])
    A = UA / U
    pa = r["p2"]; pg = E.P0 * 1.02 * 1.035
    dpa, dpg = 0.025 * pa, 0.035 * E.P0 * 1.02
    rho_a, rho_g = pa / (E.RA * Ta), pg / (E.Rg(b) * Tg)
    Aca = math.sqrt(REC["FRE"] * ML.mu_gas(Ta) * ma * A / (2 * rho_a * dpa * REC["DH"]))
    Acg = math.sqrt(REC["FRE"] * ML.mu_gas(Tg) * mg * A / (2 * rho_g * dpg * REC["DH"]))
    La, Lg = A * REC["DH"] / (4 * Aca), A * REC["DH"] / (4 * Acg)
    Vc = (Aca * La + Acg * Lg) / REC["POROSITY"]
    V = FS.REC_HDR * Vc
    a_max = 1.1
    mods = max(1, math.ceil(V / a_max ** 3))
    return dict(A_m2=A, UA_W_K=UA, NTU=ntu, Cr=Cr, U_W_m2K=U, V_core_m3=Vc, V_m3=V, modules=mods, side=(V / mods) ** (1 / 3),
                mass_kg=FS.REC_HDR * A * 0.00015 * 7900, L_air_m=La, L_gas_m=Lg, Ac_air_m2=Aca, Ac_gas_m2=Acg, v_air=ma / (rho_a * Aca),
                v_gas=mg / (rho_g * Acg), dp_air_pct=2.5, dp_gas_pct=3.5, eps=eps,
                L0_V_m3=FS.REC_HDR * (eps / (1 - eps)) * G * 1100.0 / FS.REC_K / FS.REC_BETA)


def _env_D(c, g):
    """D_env как funnel_space.derive (до рекуператора), затем рекуператор в корпусе"""
    D = g["D_case"] + 0.5 if c["cores"] == 1 else g["D_case"] * math.sqrt(c["cores"]) * 1.1 + 0.5
    if g["pt"]:
        D = max(D, 1.3 * g["pt"]["Dtip"] + 0.5)
    if g["rec"] and c["recup_integral"]:
        D = math.sqrt(D ** 2 + 4 * g["rec"]["V_m3"] / (math.pi * max(D, 0.5)))
    return D


def strength(c, tgt, perf, g, sol):
    """быстрая прочность (strength_l1) и правила"""
    out, R = dict(comp=[], turb=None, pt=None), []
    n = sol["n"]
    solid = c["rotor_join"] in ("заодно", "сварка")
    for k, (s, cc) in enumerate(zip(g["stages"], sol["ccs"])):
        if cc is None:
            continue
        T1, T2 = cc["T01"] - 273.15, cc["T02"] - 273.15
        rho = {"Ti64": 4430.0, "Al": 2800.0, "Ni": 8000.0}[s["mat"]]
        w = ST.wheel("comp", s["mat"], rho, cc["U2"], 0.5 * (T1 + T2), solid=solid)
        w.update(stage=k + 1, mat=s["mat"], T_out_C=T2, T_in_C=T1)
        out["comp"].append(w)
        if w.get("burst_margin_120") is not None:
            R.append(_rule(f"burst_c{k + 1}", w["burst_margin_120"] >= 1.0, w["burst_margin_120"], 1.0,
                           "запас по разрушающей частоте колеса компрессора при 120 % (подобие disk.py, k = 1)"))
        if s["mat"] == "Ti64":
            R.append(_rule(f"T_ti{k + 1}", T2 <= T_TI_MAX, T2, T_TI_MAX, "Ti-6Al-4V за колесом: предел данных materials.json"))
    E = FS._E()
    tg = g["tg"]
    mat = c["mat_tw"]
    rho_t = MK.MAT.get(mat, {}).get("rho", 8000.0)
    cooled = MK.rotor_cooled(c["cool"], c["m_tw"])
    dT, basis = ST.creep_material(mat)
    tm = sol["tm"]
    if c["turb"] == "Р":
        cpg = ML.cp(perf["r"]["T41"], perf["r"]["t"]["b"])
        T_rim = c["t41"] - tg["U"] ** 2 / (2 * cpg) - (ST.DT_COOL if cooled else 0.0)
        if mat == "Si3N4":
            w = ST.wheel("rturb", "Si3N4", rho_t, tg["U"], T_rim - 200.0, solid=True)
            wb = ST.weibull_pf(w["vm_MPa"], tg["Dtip"])
            w.update(weibull=wb, T_rim_C=T_rim)
            R.append(_rule("weibull_t", wb["pf"] <= ST.PF_MAX, wb["pf"], ST.PF_MAX,
                           "вероятность разрушения керамического колеса (Вейбулл; m - оценка)"))
        else:
            w = ST.wheel("rturb", mat, rho_t, tg["U"], T_rim - 200.0, solid=True)
            ex = ST._data()["cases"]["turb130_in713c"]["result"]
            s_rim = ex["creep"]["local_at"]["vm_MPa"] * w["scale_stress"]
            life, extr = ST.rupture_hours(s_rim, T_rim, dT)
            s_al, _ = ST.sigma_allow(T_rim, dT_cap=dT)
            w.update(T_rim_C=T_rim, sigma_rim_MPa=s_rim, sigma_allow_MPa=s_al, life_h=life, life_extrapolated=extr, creep_basis=basis)
            R.append(_rule("burst_t", w["burst_margin_120"] >= 1.0, w["burst_margin_120"], 1.0,
                           "запас по разрушающей частоте колеса турбины при 120 % (подобие disk.py, k = 1)"))
            R.append(_rule("creep_t", s_rim <= s_al, s_rim, s_al,
                           f"напряжение на ободе радиального колеса против длительной прочности за {ST.RUPTURE_H:.0f} ч, МПа "
                           f"(Ларсон - Миллер; {basis})"))
        out["turb"] = w
    else:
        s1 = tm["stages"][0]
        Tb = s1["T_rel_in"] - 273.15 - (ST.DT_COOL if cooled else 0.0)
        AN2 = max(s["AN2"] for s in tm["stages"])
        sig = ST.blade_root(rho_t, AN2)
        uts, _ = ST._uts("in713c", Tb - dT)
        life, extr = ST.rupture_hours(sig, Tb, dT)
        s_al, _ = ST.sigma_allow(Tb, dT_cap=dT)
        m_bl = 0.4 * g["mass"]["tw"] / c["nt"]
        dk = ST.disc_axial(rho_t if g["blisk"][0] else 8190.0, n, tg["D"], tg["h"], m_bl, Tb - 150.0,
                           "IN718" if not g["blisk"][0] else mat, dT if g["blisk"][0] else MK.MAT["IN718"]["dT_cap"])
        out["turb"] = dict(T_blade_C=Tb, AN2=AN2, sigma_root_MPa=sig, sigma_allow_MPa=s_al, overspeed_ratio=sig * ST.OVERSPEED ** 2 / uts,
                           life_h=life, life_extrapolated=extr, creep_basis=basis, disc=dk)
        R.append(_rule("blade_t", sig * ST.OVERSPEED ** 2 <= uts, sig * ST.OVERSPEED ** 2, uts,
                       "корень лопатки турбины при 120 % против σв (кривая IN-713C со сдвигом сплава)"))
        R.append(_rule("creep_bl", sig <= s_al, sig, s_al,
                       f"корень лопатки турбины против длительной прочности за {ST.RUPTURE_H:.0f} ч, МПа ({basis})"))
        R.append(_rule("disc_t", dk["ratio"] <= ST.DISC_RATIO_MAX, dk["ratio"], ST.DISC_RATIO_MAX,
                       "равнопрочный диск турбины: толщина ступицы / обода (оценка)"))
    pm = sol["pm"]
    if pm:
        s1 = pm["stages"][0]
        Tb = s1["T_rel_in"] - 273.15
        AN2 = max(s["AN2"] for s in pm["stages"])
        sig = ST.blade_root(8000.0, AN2)
        uts, _ = ST._uts("in713c", Tb)
        life, extr = ST.rupture_hours(sig, Tb, 0.0)
        s_al, _ = ST.sigma_allow(Tb)
        dk = ST.disc_axial(8190.0, tgt["n_pkg"], g["pt"]["D"], g["pt"]["h"], 0.2 * g["mass"]["pt"] / c["npt"], Tb - 150.0, "IN718",
                           MK.MAT["IN718"]["dT_cap"])
        out["pt"] = dict(T_blade_C=Tb, AN2=AN2, sigma_root_MPa=sig, overspeed_ratio=sig * ST.OVERSPEED ** 2 / uts, life_h=life,
                         life_extrapolated=extr, disc=dk)
        R.append(_rule("blade_pt", sig * ST.OVERSPEED ** 2 <= uts, sig * ST.OVERSPEED ** 2, uts, "корень лопатки СТ при 120 % (K418)"))
        R.append(_rule("creep_pt", sig <= s_al, sig, s_al, f"корень лопатки СТ против длительной прочности за {ST.RUPTURE_H:.0f} ч, МПа"))
        R.append(_rule("disc_pt", dk["ratio"] <= ST.DISC_RATIO_MAX, dk["ratio"], ST.DISC_RATIO_MAX, "равнопрочный диск СТ (оценка)"))
    return out, R


def ml_rules(c, sol):
    R = []
    for k, cc in enumerate(sol["ccs"]):
        if cc is None:
            continue
        R.append(_rule(f"M1t{k + 1}", cc["M1t"] <= M1T_MAX, cc["M1t"], M1T_MAX, "число Маха на периферии входа колеса (оценка)"))
        R.append(_rule(f"dH{k + 1}", cc["W2_W1t"] >= DH_MIN, cc["W2_W1t"], DH_MIN, "торможение в колесе W2/W1t (оценка)"))
        R.append(_rule(f"D1D2_{k + 1}", cc["D1t_D2"] <= D1D2_MAX, cc["D1t_D2"], D1D2_MAX, "D1t/D2 колеса (оценка)"))
        R.append(_rule(f"a2_{k + 1}", cc["alpha2_deg"] <= ALPHA2_MAX, cc["alpha2_deg"], ALPHA2_MAX, "угол потока за колесом (оценка)"))
        fl = cc["flags"]
        R.append(_rule(f"b2b{k + 1}", "forward" not in fl and cc["beta2b_deg"] <= BETA2_MAX,
                       -1.0 if "forward" in fl else cc["beta2b_deg"], [0.0, BETA2_MAX],
                       "угол лопатки на выходе от радиального (по скольжению Wiesner): загнутость вперёд или больше предела"))
    tm = sol["tm"]
    if c["turb"] == "О":
        for s in tm["stages"]:
            k = s["stage"]
            R.append(_rule(f"tM2_{k}", s["M2"] <= M2T_MAX, s["M2"], M2T_MAX, "число Маха за сопловым аппаратом турбины"))
            R.append(_rule(f"tMw3_{k}", s["Mw3"] <= MW3_MAX, s["Mw3"], MW3_MAX, "относительное число Маха за ротором турбины (оценка)"))
            R.append(_rule(f"tR_{k}", s["reaction"] >= R_MIN, s["reaction"], R_MIN, "реактивность ступени >= 0 (выход без закрутки)"))
            R.append(_rule(f"tM2h_{k}", s["M2_hub"] <= ML.HUB_M2_MAX, s["M2_hub"], ML.HUB_M2_MAX,
                           "число Маха за сопловым аппаратом у втулки (свободный вихрь; оценка)"))
            R.append(_rule(f"tRh_{k}", s["reaction_hub"] >= ML.R_HUB_MIN, s["reaction_hub"], ML.R_HUB_MIN,
                           "реактивность у втулки (свободный вихрь; оценка)"))
            R.append(_rule(f"tMw3t_{k}", s["Mw3_tip"] <= MW3_MAX, s["Mw3_tip"], MW3_MAX,
                           "относительное число Маха за ротором на периферии (свободный вихрь; оценка)"))
    else:
        for s in tm["stages"]:
            R.append(_rule("rM4", s["M4"] <= M4_MAX, s["M4"], M4_MAX, "число Маха за сопловым аппаратом радиальной турбины (оценка)"))
            R.append(_rule("rMw5", s["Mw5t"] <= MW5_MAX, s["Mw5t"], MW5_MAX, "число Маха на периферии выхода радиального колеса (оценка)"))
    if sol["pm"]:
        for s in sol["pm"]["stages"]:
            k = s["stage"]
            R.append(_rule(f"ptM2_{k}", s["M2"] <= M2T_MAX, s["M2"], M2T_MAX, "число Маха за сопловым аппаратом СТ"))
            R.append(_rule(f"ptR_{k}", s["reaction"] >= R_MIN, s["reaction"], R_MIN, "реактивность ступени СТ >= 0"))
            R.append(_rule(f"ptM2h_{k}", s["M2_hub"] <= ML.HUB_M2_MAX, s["M2_hub"], ML.HUB_M2_MAX, "число Маха за соплом СТ у втулки"))
            R.append(_rule(f"ptRh_{k}", s["reaction_hub"] >= ML.R_HUB_MIN, s["reaction_hub"], ML.R_HUB_MIN, "реактивность СТ у втулки"))
    return R


def viol_raw(bad):
    """сумма относительных превышений (как funnel_space.violation, без округления) + 1e-3 на каждое непройденное правило -
    ограничение оптимизатора: точка на самой границе правила считается нарушающей"""
    v = 0.0
    for r in bad:
        val, lim = r.get("val"), r.get("lim")
        if isinstance(lim, list) and isinstance(val, (int, float)):
            sc = max(abs(lim[0]), abs(lim[1]), 1e-9)                 # предел 0 (угол лопатки) - по наибольшему из пределов
            v += max(0.0, lim[0] - val) / sc + max(0.0, val - lim[1]) / sc
        elif isinstance(val, (int, float)) and isinstance(lim, (int, float)):
            v += abs(val - lim) / (abs(lim) if abs(lim) > 1e-9 else 1.0)
        else:
            v += 1.0
    return v + 1e-3 * len(bad)


def evaluate_point(c, tgt):
    """расчётная точка уровня 1 (без нерасчётных режимов): оценка, превышение, правила, всё для записи"""
    try:
        sol = solve_l1(c, tgt)
        if not sol["ok"]:
            bad = [_rule("perf", False, text=sol["fail"])]
            return dict(ok=False, score=0.0, viol=5.0 + 1e-3, bad=bad, sol=None)
        perf = sol["perf"]
        g = FS.derive(c, tgt, perf)
        rec = None
        if c["recup"] == "печатный":
            rec = recup_l1(c, perf)
            g["rec"] = dict(A_m2=rec["A_m2"], V_m3=rec["V_m3"], modules=rec["modules"], side=rec["side"], mass_kg=rec["mass_kg"])
        g["D_env"] = _env_D(c, g)
        g["parts_mk"], g["blisk"], g["pt_blisk"], g["ngv_ring"], g["case_pieces"] = FS._parts_for_make(c, g)
        rules = FS.check(c, tgt, perf, g) + ml_rules(c, sol)
        stg, R2 = strength(c, tgt, perf, g, sol)
        rules += R2
    except (ValueError, ZeroDivisionError, OverflowError) as e:
        return dict(ok=False, score=0.0, viol=5.0 + 1e-3, bad=[_rule("perf", False, text=f"не строится: {type(e).__name__}: {e}")], sol=None)
    bad = [r for r in rules if not r["ok"]]
    viol = viol_raw(bad)
    lay = layout(c, tgt, perf, g, sol, rec)
    bom, sv, pp, mk = FS.model_inputs(c, tgt, perf, g)
    pr, sr, ppr, mr = PM.count(bom), SM.score(sv), PP.score(pp), MK.score(mk)
    eta_el = perf["eta_th"] * 100 * 0.985 * 0.972
    s_eta = L0.s_eta_of(eta_el, tgt["eta0"])
    V_bay = min(tgt["pkg_W"], tgt["pkg_H"]) ** 2 * tgt["pkg_L"] * L0.BAY_L_FRAC
    s_size = L0.s_size_of(lay["envelope"]["V"], V_bay)
    s_em, s_es = L0.effmake(s_eta, s_size, mr["score"])
    s = dict(parts=pr["score"], plugplay=ppr["score"], service=sr["score"], effmake=s_em)
    score = L0.weighted(s)
    return dict(ok=True, score=score, viol=viol, bad=bad, rules=rules, sol=sol, perf=perf, g=g, rec=rec, strength=stg, layout=lay,
                scores=dict(s, make=mr["score"], eta=s_eta, size=s_size, etasize=s_es), eta_el=eta_el,
                models=dict(pr=pr, sr=sr, ppr=ppr, mr=mr, sv=sv, mk=mk))


# =========================================================================== геометрия для эскизов (оценка)
def layout(c, tgt, perf, g, sol, rec):
    """размеры и положение крупных узлов, м: ось x - от входного фланца к выхлопу; ось двигателя - на высоте 0,5·H над рамой"""
    E = FS._E()
    G = perf["G"]; Gc = G / c["cores"]
    P = tgt["P_shaft_kW"] * 1e3
    n_out = sol["n"] if c["rotors"] == 1 else tgt["n_pkg"]
    T = P / (2 * math.pi * n_out / 60)
    d_sh = (16 * T / (math.pi * 60e6)) ** (1 / 3)                # выходной вал: τ 60 МПа (оценка)
    items, x = [], 0.0
    ccs = sol["ccs"]
    first = next((cc for cc in ccs if cc), None)
    D1 = first["D1t"] if first else g["Dc"] * 0.6
    if c["inlet"] == "двух":
        inl = dict(type="двусторонний (радиальная камера с двух сторон колеса)", x0=x, L=1.0 * g["Dc"], W=1.15 * g["D_case"],
                   H=1.15 * g["D_case"])
    else:
        inl = dict(type="односторонний (осевой вход с коллектором)", x0=x, L=0.8 * D1, D=1.6 * D1)
    rho0 = E.P0 / (E.RA * E.T0)
    A_in = G / (rho0 * 25.0)
    inl.update(flange_area_m2=A_in, flange_D_equiv=math.sqrt(4 * A_in / math.pi), v_duct_m_s=25.0)
    x += inl["L"]
    comp = []
    for k, (s, cc) in enumerate(zip(g["stages"], ccs)):
        if cc is None:
            L = 0.12 * s["D"]
            comp.append(dict(stage=k + 1, type="осевая", x0=x, L=L, D=s["D"], hub_tip=s.get("htr")))
        else:
            L = cc["Lax"] + 0.30 * cc["D2"]
            comp.append(dict(stage=k + 1, type={"Ц": "центробежная", "Д": "диагональная"}[s["type"]], x0=x, L=L, D2=cc["D2"],
                             b2=cc["b2"], D1t=cc["D1t"], D1h=cc["D1h"], D_diffuser_out=cc["D4"], D_case=c["kD_case"] * cc["D2"],
                             diffuser="лопаточный" + (", переводной канал" if k < len(ccs) - 1 else ", сборная улитка / коллектор"),
                             blades=int(cc["Zfull"]), splitters=int(cc["Zfull"])))
        x += L
    comb = dict(type=c["comb"], D=g["D_comb"], L=g["L_comb"], injectors=g["n_inj"])
    if c["comb"] == "выносная":
        comb.update(placement="над корпусом, ось вертикальная (оценка)", x0=x, z_center=0.5 * g["D_case"] + 0.6 * g["D_comb"])
        x += 0.25 * g["Dc"]
    else:
        comb.update(placement="кольцевая вокруг вала" if c["comb"] != "ВР" else "барабан волнового ротора на валу", x0=x)
        x += g["L_comb"] + (0.6 * g["Dc"] if c["comb"] == "ВР" else 0.0)
    tm = sol["tm"]
    turb = dict(type="радиальная" if c["turb"] == "Р" else "осевая", stages=c["nt"], x0=x,
                L=tm["length"], D=(tm["D4"] if c["turb"] == "Р" else tm["Dm"]), Dtip=(tm["D4"] if c["turb"] == "Р" else tm["Dtip"]),
                D_hot_case=g["D_hot"], blades=c["zt"], vanes=c["zn"])
    x += turb["L"]
    pt = None
    if sol["pm"]:
        pm = sol["pm"]
        x += 0.5 * g["pt"]["Dtip"] * 0.6                         # переходный канал
        pt = dict(stages=c["npt"], x0=x, L=pm["length"], Dm=pm["Dm"], Dtip=pm["Dtip"], n_rpm=tgt["n_pkg"])
        x += pt["L"]
    rho_ex = E.P0 / (E.Rg(perf["r"]["t"]["b"]) * (perf["EGT"] + 273.15))
    A_ex = G * 1.02 / (rho_ex * 40.0)
    Dlast = (pt or turb)["Dtip"]
    exh = dict(type="радиальный коллектор, выход вверх (направление пакета - нет данных, оценка)", x0=x, L=0.9 * Dlast,
               D=1.6 * Dlast, flange_area_m2=A_ex, flange_D_equiv=math.sqrt(4 * A_ex / math.pi), v_m_s=40.0,
               T_C=perf["EGT"], from_recuperator=bool(rec))
    x += exh["L"]
    L_core = x
    D_core = g["D_env"] if c["cores"] == 1 else (g["D_case"] + 0.5)
    shaft = dict(output_end="холодный (со стороны входа)" if c["rotors"] == 1 else "горячий (за силовой турбиной)",
                 n_rpm=n_out, torque_kNm=T / 1e3, d_shaft=d_sh, coupling_flange_D=2.4 * d_sh, coupling_L=0.35,
                 axis_height_above_frame=None, axis_height_note="высота оси пакета - нет данных (запрос заказчику, навык gtu-dropin); "
                                                                "в эскизе - 0,5·H двигателя")
    supports = [dict(name="передняя опора (лапы)", x=0.15 * L_core), dict(name="задняя опора (лапы)", x=0.85 * L_core)]
    if L_core > 3.0:
        supports.append(dict(name="промежуточная опора", x=0.5 * L_core))
    mb = None
    if c["bear"] == "магнит":
        d_r = 0.2 * g["Dc"]
        mb = [dict(name="радиальный магнитный подшипник (вход)", x=0.05 * L_core, D=3.0 * d_r, L=0.25),
              dict(name="радиальный магнитный подшипник (турбина)", x=turb["x0"] - 0.1, D=3.0 * d_r, L=0.25),
              dict(name="осевой магнитный подшипник", x=inl["L"] * 0.5, D=4.0 * d_r, L=0.15)]
    W = H = D_core
    recg = None
    if rec:
        if c["recup_integral"]:
            recg = dict(shape="кольцевой в корпусе двигателя", D_outer=g["D_env"], L=0.8 * L_core, modules=rec["modules"],
                        V_m3=rec["V_m3"])
        else:
            nrow = max(1, int(L_core // max(rec["side"], 1e-6)))
            ncol = math.ceil(rec["modules"] / nrow)
            recg = dict(shape="блок кубических модулей сбоку от двигателя", module_side=rec["side"], modules=rec["modules"],
                        rows_along=min(nrow, rec["modules"]), cols=ncol, W=ncol * rec["side"] + 0.3, L=min(nrow, rec["modules"]) * rec["side"],
                        H=rec["side"], V_m3=rec["V_m3"])
            W += recg["W"]
    if c["comb"] == "выносная":
        H += 1.2 * g["D_comb"]
    cores = None
    if c["cores"] > 1:
        k = math.ceil(math.sqrt(c["cores"]))
        cores = dict(n=c["cores"], arrangement=f"{k}×{math.ceil(c['cores'] / k)} параллельно перед общей СТ", D_core=g["D_case"] + 0.5,
                     L_core=turb["x0"] + turb["L"])
        W = max(W, k * (g["D_case"] + 0.5)); H = max(H, math.ceil(c["cores"] / k) * (g["D_case"] + 0.5))
    L = L_core + shaft["coupling_L"]
    m = g["mass"]
    core_m = m["rotor"] + m["case"] + m["hot"] + m["comb"]
    mass = L0.MASS_K * (c["cores"] * core_m + m["pt"] + (rec["mass_kg"] if rec else 0.0))
    oil = c["bear"] == "масло"
    conn = [dict(kind="топливо", where="коллектор камеры" if c["comb"] != "выносная" else "горелка выносной камеры", x=comb["x0"],
                 n=1 if c["m_comb"] == "печать" else g["n_inj"])]
    if oil:
        conn += [dict(kind="масло: подвод и слив", where=s["name"], x=s["x"], n=2) for s in supports[:2]]
    conn += [dict(kind="электрика: жгут датчиков (один разъём)" if c["conn"] == "быстроразъём" else "электрика: жгут датчиков",
                  where="блок датчиков", x=0.3 * L_core, n=1)]
    if not oil:
        conn.append(dict(kind="электрика: силовые кабели магнитных подшипников", where="контроллер", x=0.1 * L_core, n=2))
    conn.append(dict(kind="воздух: пусковой перепуск / отбор", where="переводной канал или за компрессором", x=comp[-1]["x0"] if comp else 0.0, n=1))
    return dict(inlet=inl, compressor=comp, combustor=comb, turbine=turb, power_turbine=pt, exhaust=exh, recuperator=recg,
                output_shaft=shaft, supports=supports, magnetic_bearings=mb, cores=cores, connections=conn,
                envelope=dict(L=L, W=W, H=H, V=L * W * H), mass_kg_estimate=mass,
                note="уровень 1: размеры - одномерные модели и формулы derive; положения и обвязка - оценка для эскиза")


# =========================================================================== оптимизация
def _rng_seed(idx, kv):
    return [SEED, int(idx), int(kv.split("-")[1])]


def optimize(c0, tgt, starts=DEFAULT_STARTS, maxfev=DEFAULT_MAXFEV, idx=0, de_gens=None):
    """старт = прогон дифференциальной эволюции (scipy, ограничение «превышение правил <= 0» по правилам допустимости Лампинена:
    допустимая точка лучше недопустимой, из недопустимых - с меньшим превышением; популяция DE_POP·d, de_gens поколений,
    зерно = f(SEED, idx, kv, номер старта); в первую популяцию входит точка уровня 0) + доводка COBYLA с тем же ограничением от
    лучшей точки (maxfev вычислений). Кэш вычислений по x. Возврат: кандидаты (x, оценка, превышение), сводка стартов"""
    from scipy.optimize import minimize, differential_evolution, NonlinearConstraint
    de_gens = DE_GENS if de_gens is None else de_gens
    spec = var_spec(c0)
    d = len(spec)
    cache = {}
    best_feas = []

    def ev(x):
        x = np.clip(np.asarray(x, float), 0.0, 1.0)
        key = tuple(np.round(x, 9))
        if key not in cache:
            c = x_to_c(c0, x, spec)
            r = evaluate_point(c, tgt)
            cache[key] = (r["score"], r["viol"])
            if r["ok"] and r["viol"] == 0:
                best_feas.append((r["score"], key))
        return cache[key]
    x0 = c_to_x(c0, spec)
    runs = []
    for k in range(starts):
        n_before = len(cache)
        rng = np.random.default_rng(_rng_seed(idx, c0["kv"]) + [k])
        pop = rng.uniform(0.0, 1.0, (DE_POP * d, d))
        if k == 0:
            pop[0] = x0
        xd = x0
        try:
            rd = differential_evolution(lambda x: -ev(x)[0], [(0.0, 1.0)] * d,
                                        constraints=NonlinearConstraint(lambda x: ev(x)[1], -np.inf, 0.0), init=pop,
                                        maxiter=de_gens, seed=int(rng.integers(2 ** 31)), polish=False, tol=1e-6,
                                        mutation=(0.5, 1.0), recombination=0.7, updating="deferred", workers=1)
            xd = np.clip(rd.x, 0, 1)
        except (ValueError, FloatingPointError):
            pass
        n_de = len(cache) - n_before
        xf = xd
        try:
            res = minimize(lambda x: -ev(x)[0] + PENALTY * ev(x)[1], xd, method="COBYLA",
                           constraints=[dict(type="ineq", fun=lambda x: -ev(x)[1])],
                           bounds=[(0.0, 1.0)] * d, options=dict(rhobeg=0.05, maxiter=maxfev, tol=1e-5))
            xc = np.clip(res.x, 0, 1)
            if (ev(xc)[1], -ev(xc)[0]) <= (ev(xd)[1], -ev(xd)[0]):
                xf = xc
        except (ValueError, FloatingPointError):
            pass
        s, v = ev(xf)
        sd, vd = ev(xd)
        runs.append(dict(start=k, de_score=round(sd, 5), de_viol=round(vd, 5), score=round(s, 5), viol=round(v, 5),
                         x=[round(float(a), 5) for a in xf], nfev=len(cache) - n_before, nfev_de=n_de, with_l0_point=(k == 0)))
    cands = {}
    for key in [tuple(np.round(np.clip(x0, 0, 1), 9))] + [tuple(np.round(np.array(r["x"]), 9)) for r in runs] + \
               [k for _, k in sorted(best_feas, reverse=True)[:OD_CANDIDATES]]:
        cands[key] = ev(np.array(key))
    order = sorted(cands.items(), key=lambda kv: (kv[1][1] > 0, kv[1][1] if kv[1][1] > 0 else -kv[1][0]))
    return spec, order, runs, len(cache)


# =========================================================================== вариант целиком
def design_data(c, tgt, ev):
    """данные расчётной точки для нерасчётных режимов (offdesign_l1)"""
    E = FS._E()
    sol, perf, g = ev["sol"], ev["perf"], ev["g"]
    r = perf["r"]
    n = sol["n"]
    Gc = perf["G"] / c["cores"]
    stages = []
    for k, (s, cc, so) in enumerate(zip(g["stages"], sol["ccs"], perf["comp"])):
        if cc is None:
            U = s["U"]
            stages.append(dict(type="О", W_d=s["Gs"], T1_d=so["T1"], p1_d=so["p1"], n_d=n, U_d=U, psi_d=so["w"] / U ** 2, eta_d=so["eta"]))
        else:
            stages.append(dict(type=s["type"], ml=cc, dbl=cc["dbl"]))
    t = r["t"]
    K_gg = OD.capacity(t["m"] * Gc, t["T_in"], t["p_in"], t["ER"])
    nu_gg = OD._nu(perf["Ut"], t["T_in"], t["b"], t["ER"] ** (1.0 / c["nt"]))
    D = dict(c=c, cores=c["cores"], rotors=c["rotors"], stages=stages, n_d=n, n_pkg=tgt["n_pkg"], G=perf["G"], K_gg=K_gg,
             eta_gg_d=perf["eta_gg"], U_gg_d=perf["Ut"], nu_gg_d=nu_gg, nt=c["nt"], npt=c.get("npt") or 1, eta_pt_d=perf["eta_pt"],
             t41_d=c["t41"], P_rated_kW=tgt["P_shaft_kW"], rec=bool(ev["rec"]), UA=ev["rec"]["UA_W_K"] if ev["rec"] else 0.0,
             NTU_d=ev["rec"]["NTU"] if ev["rec"] else 0.0, Cr_d=ev["rec"]["Cr"] if ev["rec"] else 1.0,
             W_comp_kW=r["wc"] * perf["G"] / 1e3)
    if r["pt"]:
        p = r["pt"]
        D.update(K_pt=OD.capacity(p["m"] * perf["G"], p["T_in"], p["p_in"], p["ER"]), U_pt_d=perf["Upt"],
                 nu_pt_d=OD._nu(perf["Upt"], p["T_in"], p["b"], p["ER"] ** (1.0 / c["npt"])))
    J = 0.0
    for s in g["stages"]:
        mw = (0.0825 if s["type"] != "О" else 0.03) * {"Ti64": 4430.0, "Al": 2800.0, "Ni": 8000.0}[s["mat"]] * s["D"] ** 3
        J += 0.30 * mw * (s["D"] / 2) ** 2
    J += 0.30 * g["mass"]["tw"] * (g["tg"]["Dtip"] / 2) ** 2 + 0.5 * g["mass"]["shaft"] * (0.1 * g["Dc"]) ** 2
    D["J_kgm2"] = J * c["cores"]
    return D


OD_ISO_DT41 = 3.0        # ISO 100 % нерасчётной модели: T41 не дальше 3 K от расчётной (исправление скептика уровня 1)
OD_ISO_DP = 0.005        # и мощность не дальше 0,5 %
VIGV_Z = 15              # лопаток ВНА на один вход (оценка: малые центробежные ступени - 13-19)
VIGV_ACT_KG = 15.0       # привод ВНА (электромеханический), кг (оценка)


def full_check(c, tgt, ev):
    """нерасчётные режимы, помпаж, пуск -> правила: ISO 100 % воспроизводится (T41 и мощность), все режимы −30…+45 °C × 50-100 %
    решены (с ВНА и ветвью запирания), запас по помпажу на всех, самоподдержание при пуске"""
    D = design_data(c, tgt, ev)
    mx = OD.matrix(D)
    stt = OD.start(D)
    iso = mx["iso_check"]
    iso_ok = (iso is not None and abs(iso["T41_C"] - iso["T41_design_C"]) <= OD_ISO_DT41
              and abs(iso["P_kW"] / iso["P_design_kW"] - 1) <= OD_ISO_DP)
    R = [_rule("od_iso", iso_ok, None if iso is None else round(iso["T41_C"] - iso["T41_design_C"], 3), OD_ISO_DT41,
               "нерасчётная модель воспроизводит режим ISO 100 % (T41 ± 3 K, мощность ± 0,5 %)"),
         _rule("od_cover", mx["n_unsolved"] == 0, mx["n_unsolved"], 0,
               "все режимы −30…+45 °C × 50-100 % решены (ВНА A1_STEPS, у одновального - ветвь запирания; иначе нужен перепуск - Т-07)"),
         _rule("sm_min", mx["SM_min"] is not None and mx["SM_min"] >= OD.SM_MIN, mx["SM_min"] if mx["SM_min"] is not None else -1.0,
               OD.SM_MIN, "наименьший запас по помпажу на установившихся режимах −30…+45 °C, 50-100 % (граница - вершина π или "
                          "признак срыва; оценка)"),
         _rule("start_ss", stt["n_self_sustain"] is not None, stt["n_self_sustain"] or 0.0, 0.9,
               "самоподдержание при пуске (частота <= 0,9 расчётной, T41 пуска)")]
    return mx, stt, R


def add_vigv(c, tgt, ev, vigv=True, heat=False):
    """узлы, которых потребовали нерасчётные режимы. ВНА: лопатки, рычаги, кольцо синхронизации, привод (оценка VIGV_Z на вход,
    × 2 при двустороннем входе, × ядра), модуль замены - привод ВНА. Подогрев входа: клапан перепуска, коллектор на входе и
    2 трубы на ядро (оценка). Баллы деталей и монтажа и взвешенная оценка - заново. Модели деталей и монтажа уровня 0 не
    меняются (их отпечаток сверяется со входом) - узлы добавляются здесь"""
    m = ev["models"]
    bom, sv, pp, mk = FS.model_inputs(c, tgt, ev["perf"], ev["g"])
    pr = PM.count(bom)
    keyed = c["conn"] == "быстроразъём"
    units = []
    if vigv:
        nz = VIGV_Z * (2 if c["inlet"] == "двух" else 1) * c["cores"]
        n_act = (2 if c["inlet"] == "двух" else 1) * c["cores"]
        units.append(("ВНА", dict(parts=2 * nz + 2 * n_act, moving=nz + 2 * n_act, fasteners=nz, assemblies=n_act, cooled=0, welds=0),
                      f"лопатки ({nz}) + рычаги ({nz}) + кольцо синхронизации и привод на вход ({n_act}); штифты рычагов - крепёж (оценка)"))
        sv["modules"].append(dict(name="привод ВНА", mass_kg=VIGV_ACT_KG,
                                  joints=[dict(type="хомут") if keyed else dict(type="болты", n_bolts=4)],
                                  pipes=0, cables=1, cables_keyed=keyed, frequent=False, access=[]))
    if heat:
        nc = c["cores"]
        units.append(("подогрев входа", dict(parts=4 * nc, moving=nc, fasteners=8 * nc, assemblies=nc, cooled=0, welds=0),
                      "клапан перепуска + коллектор на входе + 2 трубы на ядро (оценка)"))
    for name, add, rule in units:
        for k, v in add.items():
            pr[k] = pr.get(k, 0) + v
        pr["by_unit"][name] = dict(add, rule=rule)
    n_eff = pr["parts"] + PM.K_ASSEMBLY * pr["assemblies"] + PM.K_FASTENER * pr["fasteners"] + PM.K_MOVING * pr["moving"] + PM.K_COOLED * pr["cooled"]
    pr.update(n_eff=round(n_eff, 2), score=round(PM.score_from_neff(n_eff), 3))
    sr = SM.score(sv)
    s = dict(ev["scores"])
    s.update(parts=pr["score"], service=sr["score"])
    lay = json.loads(json.dumps(ev["layout"], default=float))
    if vigv:
        lay["inlet"]["vigv"] = dict(vanes_per_inlet=VIGV_Z, inlets=(2 if c["inlet"] == "двух" else 1), actuator_kg=VIGV_ACT_KG,
                                    a1_range_deg=[0.0, max(OD.A1_STEPS)], note="поворотный ВНА перед 1-й ступенью (оценка)")
        lay["connections"].append(dict(kind="электрика: привод ВНА", where="вход", x=0.5 * lay["inlet"]["L"], n=1))
    if heat:
        lay["inlet"]["bleed_heat"] = dict(note="коллектор подогрева входа перепуском воздуха за компрессором (оценка)")
        lay["connections"].append(dict(kind="воздух: перепуск на подогрев входа", where="за компрессором -> вход",
                                       x=lay["compressor"][-1]["x0"] if lay["compressor"] else 0.0, n=1))
    ev = dict(ev, score=L0.weighted({k: s[k] for k in ("parts", "plugplay", "service", "effmake")}), scores=s, layout=lay,
              models=dict(m, pr=pr, sr=sr, sv=sv), vigv=bool(vigv), inlet_heat=bool(heat))
    return ev


def _r(x, k=4):
    if x is None:
        return None
    if isinstance(x, (bool, np.bool_)):
        return bool(x)
    if isinstance(x, (int, np.integer)):
        return int(x)
    if isinstance(x, (float, np.floating)):
        return round(float(x), k) if math.isfinite(x) else None
    return x


def _j(vals, k=4):
    return "/".join("" if v is None else f"{round(float(v), k):g}" for v in vals)


CONT_AXES = ["t41", "pi", "split", "mu1", "mu2", "phi1", "lt", "kcm", "zc", "zt", "zn", "kUpt", "phit", "phipt", "eps_rec", "kD_case", "pgain"]
COLUMNS = (["j", "idx", "kv", "arch", "niche", "l0_rank", "l0_score", "l0_feasible", "l0_eta_th"]
           + [a for a in FS.CAT if a != "arch"] + CONT_AXES
           + ["feasible", "reject", "reject_detail", "violation", "score", "d_score_vs_l0", "start_score", "start_feasible", "opt_gain",
              "eta_th_shaft", "eta_el", "s_parts", "s_plugplay", "s_service", "s_make", "s_eta", "s_size", "s_effmake",
              "G", "dG_pct", "EGT", "dEGT", "n_rpm", "pi_stages", "eta_c", "eta_c_model", "M1t", "beta2b", "phi_c", "D2", "U2",
              "eta_t", "eta_t_model", "M2_t", "nu_t", "ns_t", "psi_t", "eta_pt", "burst_min_120", "life_min_h", "pf_ceramic",
              "disc_ratio_max", "rec_V", "rec_V_L0", "SM_min", "SM_min_at", "od_unsolved", "vigv", "a1_max", "n_choke", "choke_x_min",
              "inlet_heat", "x_bleed_max",
              "P_avail_hot_kW", "eta_th_50_iso", "n_self_sustain",
              "P_starter_kW", "L_env", "W_env", "H_env", "V_env", "mass_kg", "n_parts", "n_assemblies", "n_modules", "steps",
              "n_tools", "trl_min", "pp_blocked", "evals", "starts", "t_opt_s", "t_s", "t_cpu_s"])


def evaluate_variant(j, entry, kv, n, seed, starts, maxfev, detail=True, c0=None, de_gens=None):
    """вариант отбора: концепция уровня 0 -> оптимизация -> проверка кандидатов с нерасчётными режимами -> запись.
    c0 - готовая концепция (масштабирование на другую точку KV: дискретные оси те же, непрерывные - окна design())"""
    t0, tc0 = time.time(), time.process_time()
    idx = int(entry["idx"])
    scaled = c0 is not None
    c0 = c0 if scaled else FS.concept(idx, n, kv, seed)
    tgt = FS.target(kv)
    spec, order, runs, nev = optimize(c0, tgt, starts, maxfev, idx, de_gens)
    t_opt = time.time() - t0
    ev0 = evaluate_point(x_to_c(c0, c_to_x(c0, spec), spec), tgt)
    chosen, tried = None, []
    for key, (s, v) in order[:OD_CANDIDATES]:
        if chosen is not None and not chosen[4] and (v > 0 or s <= chosen[1]["score"] + 1e-9):
            break            # допустимый уже найден; узел ВНА только снижает оценку - смотреть дальше незачем, если не выше
        c = x_to_c(c0, np.array(key), spec)
        ev = evaluate_point(c, tgt)
        if not ev["ok"]:
            tried.append(dict(score=round(s, 4), viol=round(v, 4), result="расчётная точка не строится"))
            continue
        mx, stt, R = full_check(c, tgt, ev) if ev["viol"] == 0 else (None, None, [])
        if mx and (mx["vigv_used"] or mx["heat_used"]):
            ev = add_vigv(c, tgt, ev, vigv=mx["vigv_used"], heat=mx["heat_used"])
        bad = ev["bad"] + [r for r in R if not r["ok"]]
        tried.append(dict(score=round(ev["score"], 4), viol=round(viol_raw(bad), 5), reject="|".join(r["code"] for r in bad)))
        cand = (c, ev, mx, stt, bad)
        if not bad:
            if chosen is None or chosen[4] or ev["score"] > chosen[1]["score"]:
                chosen = cand
            continue
        if chosen is None or (chosen[4] and viol_raw(bad) < viol_raw(chosen[4])):
            chosen = cand
    if chosen is None:
        c = x_to_c(c0, np.array(order[0][0]), spec)
        chosen = (c, evaluate_point(c, tgt), None, None, None)
    c, ev, mx, stt, bad = chosen
    if bad is None:
        bad = ev["bad"]
    feasible = int(not bad) and ev["ok"]
    rec = {k: None for k in COLUMNS}
    rec.update(j=j, idx=idx, kv=kv, arch=c["arch"], niche=L0.niche_of(c), l0_rank=entry.get("rank"), l0_score=entry.get("score"),
               l0_feasible=int(bool(entry.get("feasible"))), l0_eta_th=None)
    for a in FS.CAT:
        if a != "arch":
            rec[a] = c[a]
    for a in CONT_AXES:
        rec[a] = _r(c.get(a), 6)
    rec.update(feasible=int(feasible), reject="|".join(r["code"] for r in bad), reject_detail=L0._detail_rules(bad, 4),
               violation=round(viol_raw(bad), 4) if bad else 0.0, start_score=_r(ev0["score"]), start_feasible=int(ev0["ok"] and ev0["viol"] == 0),
               evals=nev, starts=starts)
    l0row = (L0.evaluate_l0(idx, n, kv, seed) if not scaled else
             {k: None for k in ("eta_th_shaft", "s_parts", "s_plugplay", "s_service", "s_make", "s_eta", "s_size", "V_env")})
    rec["l0_eta_th"] = l0row["eta_th_shaft"]
    det = None
    if ev["ok"]:
        perf, sol, g, lay, m = ev["perf"], ev["sol"], ev["g"], ev["layout"], ev["models"]
        ccs = [cc for cc in sol["ccs"] if cc]
        tm, pm = sol["tm"], sol["pm"]
        stg = ev["strength"]
        # запас при 120 %: колёса (подобие disk.py) и лопатки осевых турбин (σв / σ корня при 120 %); равнопрочный диск осевой
        # ступени здесь не входит - у него запас 1,0 по построению, его проверка - отношение толщин (disc_ratio_max)
        bursts = [w["burst_margin_120"] for w in stg["comp"] if w.get("burst_margin_120") is not None]
        lives = []
        for x in (stg["turb"], stg["pt"]):
            if not x:
                continue
            if x.get("burst_margin_120") is not None:
                bursts.append(x["burst_margin_120"])
            if x.get("overspeed_ratio"):
                bursts.append(1.0 / x["overspeed_ratio"])
            if x.get("life_h") is not None:
                lives.append(x["life_h"])
        rec.update(score=_r(ev["score"]), d_score_vs_l0=_r(ev["score"] - entry["score"]) if entry.get("score") is not None else None,
                   opt_gain=_r(ev["score"] - ev0["score"]) if ev0["ok"] else None,
                   eta_th_shaft=_r(perf["eta_th"] * 100, 3), eta_el=_r(ev["eta_el"], 3),
                   s_parts=ev["scores"]["parts"], s_plugplay=ev["scores"]["plugplay"], s_service=ev["scores"]["service"],
                   s_make=ev["scores"]["make"], s_eta=_r(ev["scores"]["eta"], 3), s_size=_r(ev["scores"]["size"], 3),
                   s_effmake=_r(ev["scores"]["effmake"], 3), G=_r(perf["G"], 3), dG_pct=_r((perf["G"] / tgt["G0"] - 1) * 100, 2),
                   EGT=_r(perf["EGT"], 1), dEGT=_r(perf["EGT"] - tgt["EGT0"], 1), n_rpm=_r(float(sol["n"]), 0),
                   pi_stages=_j([s["pi"] for s in perf["comp"]], 3), eta_c=_j(sol["eta_c"]),
                   eta_c_model=_j([cc["eta_p_model"] for cc in ccs]), M1t=_j([cc["M1t"] for cc in ccs], 3),
                   beta2b=_j([cc["beta2b_deg"] for cc in ccs], 1), phi_c=_j([cc["phi"] for cc in ccs]),
                   D2=_j([s["D"] for s in g["stages"]]), U2=_j([s["U"] for s in g["stages"]], 1),
                   eta_t=_r(perf["eta_gg"]), eta_t_model=_r(tm.get("eta_p_model") or tm.get("eta_is_model")),
                   M2_t=_r(max(s["M2"] for s in tm["stages"]) if c["turb"] == "О" else max(s["M4"] for s in tm["stages"]), 3),
                   nu_t=_r(tm.get("nu"), 3), ns_t=_r(tm.get("ns"), 3), psi_t=_r(tm.get("psi"), 3), eta_pt=_r(perf["eta_pt"]),
                   burst_min_120=_r(min(bursts) if bursts else None, 3), life_min_h=_r(min(lives) if lives else None, 0),
                   pf_ceramic=_r(stg["turb"].get("weibull", {}).get("pf") if stg["turb"] else None, 8),
                   disc_ratio_max=_r(max([x["disc"]["ratio"] for x in (stg["turb"], stg["pt"]) if x and x.get("disc")] or [0.0]) or None, 3),
                   rec_V=_r(ev["rec"]["V_m3"], 3) if ev["rec"] else None, rec_V_L0=_r(ev["rec"]["L0_V_m3"], 3) if ev["rec"] else None,
                   L_env=_r(lay["envelope"]["L"], 3), W_env=_r(lay["envelope"]["W"], 3), H_env=_r(lay["envelope"]["H"], 3),
                   V_env=_r(lay["envelope"]["V"], 3), mass_kg=_r(lay["mass_kg_estimate"], 0),
                   n_parts=m["pr"]["parts"], n_assemblies=m["pr"]["assemblies"], n_modules=len(m["sv"]["modules"]),
                   steps=m["sr"]["steps"], n_tools=len(m["sr"]["tools"]), trl_min=m["mr"]["trl_min"], pp_blocked="|".join(m["ppr"]["blocked"]))
        if mx:
            sa = mx["SM_min_at"]
            rec.update(od_unsolved=mx["n_unsolved"], vigv=int(bool(mx["vigv_used"])), a1_max=_r(mx["a1_max"], 1), n_choke=mx["n_choke"],
                       choke_x_min=_r(mx["choke_x_min"], 3), inlet_heat=int(bool(mx["heat_used"])), x_bleed_max=_r(mx["x_bleed_max"], 4))
            rec.update(SM_min=_r(mx["SM_min"]), SM_min_at=(f"{sa['Ta_C']:g} °C, {sa['load'] * 100:.0f} %, {sa['law']}" if sa else None),
                       P_avail_hot_kW=_r(mx["P_avail_hot_kW"], 0), eta_th_50_iso=_r(mx["eta_th_50_iso"] * 100 if mx["eta_th_50_iso"] else None, 2))
        if stt:
            rec.update(n_self_sustain=_r(stt["n_self_sustain"], 3), P_starter_kW=_r(stt["P_starter_kW"], 1))
        if detail:
            det = dict(j=j, idx=idx, kv=kv, concept={k: _r(v, 6) for k, v in c.items()}, geometry=_clean(lay),
                       meanline=dict(compressor=[_ml_cc(cc) for cc in ccs], turbine=_ml_t(tm, c["turb"]), power_turbine=_ml_t(pm, "О") if pm else None,
                                     calibration="meanline_l1.calibration_summary() - в результатах слияния"),
                       strength=_clean(stg), recuperator=_clean(ev["rec"]) if ev["rec"] else None,
                       offdesign=_clean({k: v for k, v in mx.items()}) if mx else None, start=_clean(stt) if stt else None,
                       modules=[dict(name=x["name"], mass_kg=_r(x["mass_kg"], 1)) for x in m["sv"]["modules"]],
                       service=dict(steps=m["sr"]["steps"], tools=m["sr"]["tools"], wrongable=m["sr"]["wrongable"],
                                    heaviest=m["sr"]["heaviest"], frequent=m["sr"]["frequent"]),
                       parts=dict(parts=m["pr"]["parts"], assemblies=m["pr"]["assemblies"], moving=m["pr"]["moving"], cooled=m["pr"]["cooled"],
                                  fasteners=m["pr"]["fasteners"], unique=m["pr"]["unique"],
                                  by_unit={u: {k: v for k, v in d.items() if k != "rule"} for u, d in m["pr"]["by_unit"].items()}),
                       plugplay={k: dict(state=v["state"], note=v.get("note")) for k, v in m["ppr"]["items"].items()},
                       make=dict(trl_by=m["mr"]["trl_by"], techs=m["mk"]["techs"], weeks=m["mr"]["weeks"]),
                       rules_failed=bad, optimization=dict(spec=[list(s) for s in spec], runs=runs, candidates_checked=tried, evals=nev),
                       l0=dict(score=entry.get("score"), rank=entry.get("rank"), eta_th_shaft=l0row["eta_th_shaft"], s_parts=l0row["s_parts"],
                               s_plugplay=l0row["s_plugplay"], s_service=l0row["s_service"], s_make=l0row["s_make"], s_eta=l0row["s_eta"],
                               s_size=l0row["s_size"], V_env=l0row["V_env"], reasons=entry.get("reasons")))
    rec["t_s"] = round(time.time() - t0, 2)
    rec["t_cpu_s"] = round(time.process_time() - tc0, 2)          # процессорное время - меньше зависит от чужой нагрузки (основа сметы)
    rec["t_opt_s"] = round(t_opt, 2)
    if det is not None:
        det["t_s"] = rec["t_s"]
    return {k: rec[k] for k in COLUMNS}, det


SCALE_FREE = ("inlet", "nt", "npt")     # при масштабировании design() подбирает под расход новой точки (с записью в repair)


def scale_concept(c_opt, kv, n0, seed):
    """та же дискретная схема на другой точке KV: непрерывные оси - в окнах допустимости design() по строке Соболя варианта;
    оси ниши и все прочие дискретные - без изменений, кроме SCALE_FREE (вход, число ступеней турбины и СТ - их design()
    подбирает под расход и частоту новой точки и пишет в repair)"""
    u = dict(zip(FS.AXES, FS.design_matrix(n0, seed)[int(c_opt["idx"])]))
    c = {k: c_opt[k] for k in ["arch"] + [a for a in FS.CAT if a != "arch"]}
    c.update(idx=int(c_opt["idx"]), n=int(n0), seed=int(seed), kv=kv, repair="масштабирование: дискретные оси без изменений")
    c = {k: c[k] for k in ["idx", "n", "seed", "kv", "arch"] + [a for a in FS.CAT if a != "arch"]} | {"repair": c["repair"]}
    FS.design(c, u, FS.target(kv))
    return c


def scale_variant(args):
    """масштабирование одного варианта: (j, концепция KV-3, kv, n0, seed, starts, maxfev) -> строка и разбор"""
    j, c_opt, kv, n0, seed, starts, maxfev = args
    _prep(n0, seed, kv)
    try:
        c = scale_concept(c_opt, kv, n0, seed)
    except (ValueError, ZeroDivisionError, OverflowError, KeyError) as e:
        return dict(j=j, kv=kv, feasible=0, reject="design", reject_detail=f"окна непрерывных осей не строятся: {e}"), None
    return evaluate_variant(j, dict(idx=c_opt["idx"]), kv, n0, seed, starts, maxfev, True, c0=c)


def _ml_cc(cc):
    keep = ("PRs", "eta_p", "eta_p_model", "shift", "U2", "D2", "b2", "beta2b_deg", "alpha2_deg", "M1t", "M2", "Mu2", "W2_W1t", "Df", "phi",
            "mu", "kcm", "Zfull", "D1t", "D1h", "D1t_D2", "Lax", "r4r2", "D4", "cl", "losses_pct", "flags", "dbl")
    return {k: _clean(cc[k]) for k in keep}


def _ml_t(tm, kind):
    if tm is None:
        return None
    if kind == "О":
        return dict(type="осевая", eta_p=_r(tm["eta_p"]), eta_p_model=_r(tm["eta_p_model"]), shift=_r(tm["shift"]), psi=_r(tm["psi"]),
                    phi=_r(tm["phi"]), Dm=_r(tm["Dm"]), Dtip=_r(tm["Dtip"]), length=_r(tm["length"]), stages=_clean(tm["stages"]))
    return dict(type="радиальная", eta_p=_r(tm["eta_p"]), eta_is=_r(tm["eta_is"]), eta_is_model=_r(tm["eta_is_model"]), nu=_r(tm["nu"]),
                ns=_r(tm["ns"]), D4=_r(tm["D4"]), length=_r(tm["length"]), stages=_clean(tm["stages"]))


def _clean(o):
    if isinstance(o, dict):
        return {str(k): _clean(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [_clean(v) for v in o]
    return _r(o, 5)


# =========================================================================== вход, части, пул
def load_select(path, strict=True):
    d = json.load(open(path, encoding="utf-8"))
    n, seed, kv = int(d["n"]), int(d["seed"]), d["kv"]
    fp = FS.fingerprint(n, seed)
    sha = L0.code_sha()
    problems = []
    if fp != d.get("fingerprint"):
        problems.append(f"отпечаток выборки {d.get('fingerprint')} != локальный {fp} (версии scipy/numpy)")
    if sha != d.get("code_sha"):
        problems.append(f"код уровня 0 {d.get('code_sha')} != текущий {sha} (отбор другой версии кода)")
    if problems and strict:
        raise SystemExit("вход уровня 1 отклонён: " + "; ".join(problems))
    return d, problems


def order_entries(sel, top=None):
    """порядок уровня 1 - двойной круг: по архетипам по кругу, внутри архетипа - по нишам (сначала место 1 каждой ниши по убыванию
    оценки, затем места 2, ...); затем почти допустимые. Так первые M вариантов любого M - тоже отбор с разнообразием (основа
    сходимости M/4, M/2, M и коротких проб). j - номер в этом порядке (единственное деление на части: j % shards)"""
    feas = [e for e in sel["selected"] if e.get("feasible")]
    near = sorted([e for e in sel["selected"] if not e.get("feasible")], key=lambda e: (e["violation"], e["idx"]))
    by = {}
    for e in feas:
        by.setdefault(e["niche"], []).append(e)
    for v in by.values():
        v.sort(key=lambda e: (-(e["score"] or 0), e["idx"]))
    best = {k: (-(v[0]["score"] or 0), v[0]["idx"]) for k, v in by.items()}
    queues = {}
    for k, v in by.items():
        for r, e in enumerate(v):
            queues.setdefault(e["arch"], []).append(((r,) + best[k], e))
    for a in queues:
        queues[a].sort(key=lambda x: x[0])
    archs = sorted(queues, key=lambda a: queues[a][0][0])
    out, pos = [], {a: 0 for a in archs}
    while any(pos[a] < len(queues[a]) for a in archs):
        for a in archs:
            if pos[a] < len(queues[a]):
                out.append(queues[a][pos[a]][1])
                pos[a] += 1
    out += near
    return out[:top] if top else out


def code_sha():
    h = hashlib.sha256()
    for f in ("funnel_l1.py", "meanline_l1.py", "strength_l1.py", "offdesign_l1.py"):
        h.update(open(os.path.join(HERE, f), "rb").read())
    h.update(L0.code_sha().encode())
    return h.hexdigest()[:16]


def _prep(n, seed, kv):
    FS._E(); L0.calib_ns(); FS.target(kv); FS.design_matrix(n, seed); ML.calibration_summary(); ST._data()


def _chunk(args):
    items, kv, n, seed, starts, maxfev, detail, de_gens = args
    _prep(n, seed, kv)
    return [evaluate_variant(j, e, kv, n, seed, starts, maxfev, detail, de_gens=de_gens) for j, e in items]


def evaluate_many(items, kv, n, seed, starts, maxfev, workers=1, detail=True, de_gens=None):
    """(j, запись отбора) -> (строка, разбор); порядок - по j; не зависит от числа процессов"""
    items = list(items)
    _prep(n, seed, kv)
    if workers <= 1 or len(items) < 2:
        got = _chunk((items, kv, n, seed, starts, maxfev, detail, de_gens))
    else:
        import multiprocessing as mp
        from concurrent.futures import ProcessPoolExecutor
        ctx = mp.get_context("fork") if "fork" in mp.get_all_start_methods() else None
        chunks = [[it] for it in items]                    # вариант - десятки секунд: по одному на задачу пула
        with ProcessPoolExecutor(workers, mp_context=ctx) as ex:
            got = [r for part in ex.map(_chunk, [(ch, kv, n, seed, starts, maxfev, detail, de_gens) for ch in chunks]) for r in part]
    return sorted(got, key=lambda rd: rd[0]["j"])


def default_workers():
    return L0.default_workers()


def method_text():
    return ("воронка, уровень 1: варианты отбора уровня 0 по нишам; для каждого - оптимизация непрерывных осей внутри дискретной схемы "
            "(несколько стартов: дифференциальная эволюция scipy с ограничением по правилам Лампинена + доводка COBYLA; зерно от SEED, "
            "idx, kv и старта) по взвешенной оценке уровня 1 с ограничением «превышение правил = 0»; "
            "одномерные модели ступеней (meanline_l1: Oh/Yoon/Chung, Soderberg + зазор Kacker - Okapuu, вид Baines; привязка уровня к "
            "формулам трека), цикл funnel_space._core (таблицы Cantera), рекуператор по эффективности, быстрая прочность (подобие "
            "disk.py, Ларсон - Миллер IN-713C, Вейбулл Si3N4), нерасчётные режимы −30…+45 °C и 50-100 %, запас по помпажу, пуск "
            "(offdesign_l1); баллы деталей, Plug&Play, монтажа, изготовимости - модели уровня 0; КПД и размер - уровня 1. Всё - оценка.")


def part_path(out, kv, M, shard, shards, tag="", gz=False):
    return os.path.join(out, f"{PART_PREFIX}{tag}_{kv}_m{M}_s{shard:03d}of{shards:03d}.json" + (".gz" if gz else ""))


def run_parts(input_path, top, shard_list, shards, out, starts, maxfev, workers=1, gz=False, tag="", strict=True, detail=True, de_gens=None):
    t0 = time.time()
    sel, problems = load_select(input_path, strict)
    kv, n, seed = sel["kv"], int(sel["n"]), int(sel["seed"])
    entries = order_entries(sel, top)
    M = len(entries)
    idx_by = {k: FS.shard_indices(M, k, shards) for k in shard_list}
    alljs = sorted(j for v in idx_by.values() for j in v)
    de_gens = DE_GENS if de_gens is None else de_gens
    got = evaluate_many([(j, entries[j]) for j in alljs], kv, n, seed, starts, maxfev, workers, detail, de_gens)
    by_j = {r["j"]: (r, d) for r, d in got}
    el = time.time() - t0
    os.makedirs(out, exist_ok=True)
    sha_in = hashlib.sha256(open(input_path, "rb").read()).hexdigest()
    ord_sha = hashlib.sha256(json.dumps([e["idx"] for e in entries]).encode()).hexdigest()[:16]
    meta = dict(method=method_text(), software=FS.software(), date=DATE, level=LEVEL, kv=kv, target=FS.target(kv), n0=n, seed=seed,
                fingerprint=sel.get("fingerprint"), l0_code_sha=sel.get("code_sha"), input=os.path.basename(input_path), input_sha256=sha_in,
                order_sha=ord_sha, M=M, top=top, shards=shards, starts=starts, maxfev=maxfev, de_gens=de_gens, de_pop=DE_POP, tag=tag,
                code_sha=code_sha(),
                input_problems=problems, weights=L0.WEIGHTS, workers=workers,
                call=dict(shards=shard_list, variants=len(alljs), elapsed_s=round(el, 2),
                          s_per_variant_core=round(el * workers / max(1, len(alljs)), 3),
                          s_per_variant_mean=round(float(np.mean([by_j[j][0]["t_s"] for j in alljs])), 3) if alljs else None,
                          host_cpus=os.cpu_count(), started=datetime.datetime.fromtimestamp(t0, datetime.timezone.utc).isoformat()))
    paths = []
    for k in shard_list:
        doc = dict(meta, shard=k, columns=COLUMNS, n_rows=len(idx_by[k]), rows=[[by_j[j][0][c] for c in COLUMNS] for j in idx_by[k]],
                   details=[by_j[j][1] for j in idx_by[k]])
        p = part_path(out, kv, M, k, shards, tag, gz)
        txt = json.dumps(doc, ensure_ascii=False, separators=(",", ":"))
        if gz:
            with gzip.open(p, "wt", encoding="utf-8") as fh:
                fh.write(txt)
        else:
            with open(p, "w", encoding="utf-8") as fh:
                fh.write(txt)
        paths.append(p)
    return paths, meta


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--input", required=True, help="отбор уровня 0 (merge_l0: funnel_l0_select*.json)")
    ap.add_argument("--top", type=int, default=0, help="первых вариантов в порядке ниш по кругу (0 - все)")
    ap.add_argument("--shards", type=int, required=True)
    ap.add_argument("--shard", default="all")
    ap.add_argument("--starts", type=int, default=DEFAULT_STARTS, help="стартов оптимизатора на вариант")
    ap.add_argument("--maxfev", type=int, default=DEFAULT_MAXFEV, help="вычислений цели на доводку COBYLA после каждого старта")
    ap.add_argument("--de-gens", type=int, default=DE_GENS, help="поколений дифференциальной эволюции на старт")
    ap.add_argument("--out", default=os.environ.get("GTU_OUT", "."))
    ap.add_argument("--workers", type=int, default=0)
    ap.add_argument("--gzip", action="store_true")
    ap.add_argument("--tag", default="", help="добавка к имени файлов части (например _starts2x)")
    ap.add_argument("--local", action="store_true", help="локальная проба: не больше 12 вариантов и 2 стартов (решение владельца)")
    ap.add_argument("--no-strict", action="store_true", help="не отклонять вход с другим отпечатком (только отладка)")
    a = ap.parse_args(argv)
    if a.local and ((a.top or 10 ** 9) > LOCAL_MAX_VARIANTS or a.starts > LOCAL_MAX_STARTS):
        ap.error(f"локально - не больше {LOCAL_MAX_VARIANTS} вариантов и {LOCAL_MAX_STARTS} стартов")
    sl = L0.parse_shards(a.shard, a.shards)
    w = a.workers or default_workers()
    paths, meta = run_parts(a.input, a.top or None, sl, a.shards, a.out, a.starts, a.maxfev, w, a.gzip, a.tag, not a.no_strict,
                            de_gens=a.de_gens)
    c = meta["call"]
    print(f"{meta['kv']} уровень 1: M={meta['M']}, части {sl[0]}..{sl[-1]} ({len(sl)} из {a.shards}), {c['variants']} вариантов, "
          f"{c['elapsed_s']} с, {w} процессов, {c['s_per_variant_core']} с на вариант на ядро; файлов {len(paths)} -> {a.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
