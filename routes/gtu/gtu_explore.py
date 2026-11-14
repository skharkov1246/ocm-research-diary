"""Перебор конструктивных схем простого литого ГТД класса 6 МВт.
Уровень: отбор облика (0D-цикл + 1D-размеры + ограничения материалов и литья).
Свойства газа - таблицы из Cantera (воздух и стехиометрические продукты метана), смесь по доле продуктов.
Все ограничения по окружным скоростям, температурам и размерам отливок - ПРИНЯТЫЕ (оценка), см. LIM."""
import math, json, itertools, sys
import numpy as np
import cantera as ct
import os as _os; _os.chdir(_os.path.dirname(_os.path.abspath(__file__)))

# ---------------------------------------------------------------- свойства газа
AIRX = "O2:0.2095,N2:0.7808,AR:0.0093,CO2:0.0004"
TG = np.arange(200.0, 2401.0, 2.0)
_g = ct.Solution("gri30.yaml")
def _tab(setter):
    h, s = [], []
    for t in TG:
        setter(t); h.append(_g.enthalpy_mass); s.append(_g.entropy_mass)
    return np.array(h), np.array(s), ct.gas_constant / _g.mean_molecular_weight
HA, SA, RA = _tab(lambda t: setattr(_g, "TPX", (t, ct.one_atm, AIRX)))
_g.TPX = 300, ct.one_atm, AIRX
Ya = _g.Y.copy(); gi = _g.species_index; MW = _g.molecular_weights
AFST = 2 * MW[gi("O2")] / MW[gi("CH4")] / Ya[gi("O2")]        # кг воздуха на кг метана
Ys = Ya * AFST; Ys[gi("O2")] = 0.0
Ys[gi("CO2")] += MW[gi("CO2")] / MW[gi("CH4")]; Ys[gi("H2O")] += 2 * MW[gi("H2O")] / MW[gi("CH4")]
Ys /= Ys.sum()
HS, SS, RS = _tab(lambda t: setattr(_g, "TPY", (t, ct.one_atm, Ys)))
_g.TPX = 288.15, ct.one_atm, "CH4:1"; HF = _g.enthalpy_mass
LHV = 50.03e6
T0, P0 = 288.15, 101325.0

def h(T, b=0.0): return float(np.interp(T, TG, (1 - b) * HA + b * HS))
def s0(T, b=0.0): return float(np.interp(T, TG, (1 - b) * SA + b * SS))
def Rg(b=0.0): return (1 - b) * RA + b * RS
def T_h(hv, b=0.0): return float(np.interp(hv, (1 - b) * HA + b * HS, TG))
def T_s(sv, b=0.0): return float(np.interp(sv, (1 - b) * SA + b * SS, TG))
def far_for(T3, T4):
    """отношение топливо/воздух для нагрева воздуха от T3 до T4 (полное сгорание)"""
    return (h(T4) - h(T3)) / (HF + AFST * h(T4) - (1 + AFST) * h(T4, 1.0))
def bfrac(f): return f * (1 + AFST) / (1 + f)

# ---------------------------------------------------------------- принятые ограничения (оценка)
LIM = dict(
    U_ti=620.0,        # фрезерованное титановое колесо, м/с
    U_al=520.0,        # фрезерованное алюминиевое колесо, м/с (T за колесом <= 200 C)
    U_mono=480.0,      # цельнолитой монотор из жаропрочного сплава, сторона компрессора, м/с
    T_al=473.0,
    U_rad900=500.0,    # литое радиальное колесо турбины при T газа <= 900 C; выше -0.7 м/с на градус
    T_rad_max=950.0,   # C, литое равноосное радиальное колесо без охлаждения
    U_ax=420.0, T_ax_max=1050.0,     # осевые неохлаждаемые литые лопатки (монокристалл)
    T_rec_max=700.0,   # C, газ на входе в рекуператор
    D_mono_ok=0.45, D_mono_max=0.60,     # м, цельнолитой монотор
    D_rad_ok=0.50, D_rad_max=0.70,       # м, цельнолитое радиальное колесо
    D_axcast=0.50,                       # м, цельнолитое осевое колесо; больше - диск + лопатки
    D_ti_max=0.90,
)

# ---------------------------------------------------------------- узлы
def comp_stage(T1, p1, PRs, m, kind, dbl=False):
    eta = 0.880 - 0.010 * (PRs - 2.5) - {"ti": 0.0, "al": 0.0, "mono": 0.02}[kind]
    for _ in range(2):
        T2 = T_s(s0(T1) + RA * math.log(PRs) / eta); w = h(T2) - h(T1)
        U2 = math.sqrt(w / 0.70)
        rho1 = p1 / (RA * T1); mm = m / 2 if dbl else m
        D2 = math.sqrt(mm / (rho1 * U2 * 0.085))
        rho2 = rho1 * PRs ** 0.62 / ((T1 + 0.6 * (T2 - T1)) / T1)
        b2 = mm / (math.pi * D2 * rho2 * 0.28 * U2)
        cl = max(0.3e-3, 0.0007 * D2)
        eta_adj = 0.880 - 0.010 * (PRs - 2.5) - {"ti": 0.0, "al": 0.0, "mono": 0.02}[kind] - 0.3 * cl / b2
        if abs(eta_adj - eta) < 1e-4: break
        eta = eta_adj
    ulim = {"ti": LIM["U_ti"], "al": LIM["U_al"], "mono": LIM["U_mono"]}[kind]
    ok = U2 <= ulim and not (kind == "al" and T2 > LIM["T_al"])
    return dict(T2=T2, w=w, U2=U2, D2=D2, rpm=60 * U2 / (math.pi * D2), eta=eta, ok=ok, b2=b2)

def u_rad_lim(TIT_C): return LIM["U_rad900"] - 0.7 * max(0.0, TIT_C - 900.0)

def turb_radial(T1, p1, b, w_req, U_fixed=None, TIT_C=900.0, rpm=None, eta_max=0.86):
    """радиальная центростремительная турбина. w_req - требуемая работа, Дж/кг газа; None -> задаётся ER через p_out.
    возвращает ER, T2, U, D, eta. КПД по U/C0 (оптимум 0.70)."""
    R = Rg(b); h1 = h(T1, b); s1 = s0(T1, b); ulim = u_rad_lim(TIT_C)
    def run(ER):
        T2s = T_s(s1 - R * math.log(ER), b); dhs = h1 - h(T2s, b); C0 = math.sqrt(2 * dhs)
        U = U_fixed if U_fixed else min(ulim, 0.70 * C0)
        nu = U / C0
        eta = eta_max * max(0.0, 1 - ((nu - 0.70) / 0.70) ** 2)
        return eta * dhs, U, eta, nu
    lo, hi = 1.05, 12.0
    if run(hi)[0] < w_req: return None
    for _ in range(50):
        mid = 0.5 * (lo + hi)
        if run(mid)[0] < w_req: lo = mid
        else: hi = mid
    ER = hi; w, U, eta, nu = run(ER)
    T2 = T_h(h1 - w, b)
    return dict(ER=ER, T2=T2, U=U, eta=eta, nu=nu, w=w, ok=(U <= ulim + 1e-6))

def turb_radial_to_p(T1, b, ER, U_fixed=None, TIT_C=900.0, eta_max=0.86):
    R = Rg(b); h1 = h(T1, b); T2s = T_s(s0(T1, b) - R * math.log(ER), b)
    dhs = h1 - h(T2s, b); C0 = math.sqrt(2 * dhs); ulim = u_rad_lim(TIT_C)
    U = U_fixed if U_fixed else min(ulim, 0.70 * C0); nu = U / C0
    eta = eta_max * max(0.0, 1 - ((nu - 0.70) / 0.70) ** 2); w = eta * dhs
    return dict(ER=ER, T2=T_h(h1 - w, b), U=U, eta=eta, nu=nu, w=w, ok=(U <= ulim + 1e-6))

def eta_axial(psi): return 0.885 - 0.035 * max(0.0, psi - 1.2)

def turb_axial(T1, b, n, U, w_req=None, ER=None):
    """осевая неохлаждаемая турбина, n ступеней, окружная скорость на среднем диаметре U."""
    R = Rg(b); h1 = h(T1, b); s1 = s0(T1, b)
    if w_req is not None:
        psi = w_req / n / U ** 2
        if psi > 2.5: return None
        eta = eta_axial(psi)
        T2 = T_h(h1 - w_req, b)
        ERx = math.exp((s1 - s0(T2, b)) / (eta * R))          # eta как политропный
        return dict(ER=ERx, T2=T2, eta=eta, psi=psi, w=w_req, U=U)
    T2 = T1
    for _ in range(4):                                        # eta зависит от psi -> итерации
        psi = (h1 - h(T2, b)) / n / U ** 2 if T2 < T1 else 1.5
        eta = eta_axial(psi)
        T2 = T_s(s1 - eta * R * math.log(ER), b)
    w = h1 - h(T2, b); psi = w / n / U ** 2
    if psi > 2.5: return None
    return dict(ER=ER, T2=T2, eta=eta, psi=psi, w=w, U=U)

# ---------------------------------------------------------------- цикл одной схемы
def evaluate(c, PR, TIT_C, P=6.0e6):
    """c - схема (словарь), PR - общая степень сжатия, TIT_C - температура газа перед турбиной, C."""
    N, eps = c["N"], c["eps"]; TIT = TIT_C + 273.15
    nsp, ncs = c["spools"], c["cstages"]
    kind = c["comp"]; dbl = c.get("dbl", False)
    cluster = N > 1
    dp_comb = 0.05 if cluster else 0.045
    dp_duct = 0.025 if cluster else 0.015
    x_leak = 0.04
    m_core = c.get("_m", 6.0)
    out = None
    for it in range(3):
        # --- компрессоры
        PRs = PR ** (1.0 / (nsp * ncs)); T, p = T0, P0; sp = []; ok = True
        for i in range(nsp):
            wsp = 0.0; st = []
            for j in range(ncs):
                k = kind if not (kind == "al" and (i > 0 or j > 0)) else "ti"     # вторая ступень горячее - титан
                r = comp_stage(T, p, PRs, m_core, k, dbl); ok &= r["ok"]
                wsp += r["w"]; T = r["T2"]; p *= PRs; st.append(r)
            sp.append(dict(w=wsp, st=st))
            if c.get("ic") and i < nsp - 1:
                T = 308.15; p *= 0.98
        if not ok: return None
        T2, p2 = T, p
        # --- рекуператор (воздух), камера
        dpa, dpg = (0.025, 0.035) if eps > 0 else (0.0, 0.0)
        p4 = p2 * (1 - dpa) * (1 - dp_comb)
        T3 = c.get("_T3", T2)
        f = far_for(T3, TIT) * (1 - x_leak) / 0.995                # кг топлива на кг воздуха на входе
        mg = (1 - x_leak) + f; b = bfrac(f / (1 - x_leak))          # газ через турбину на 1 кг воздуха
        # --- турбины газогенератора (от ВД к НД)
        Tg, pg = TIT, p4; tg = []; hot_first = True
        for i in reversed(range(nsp)):
            w_req = sp[i]["w"] / 0.99 / mg
            st0 = sp[i]["st"][0]
            if c["shaft"] == "single" and i == 0:
                break
            if c["ggt"] == "radial":
                Ufix = st0["U2"] * c.get("kD", 1.0) if c["rotor"] == "mono" else None
                r = turb_radial(Tg, pg, b, w_req, U_fixed=Ufix, TIT_C=Tg - 273.15)
                if r is None or not r["ok"]: return None
                r["D"] = 60 * r["U"] / (math.pi * st0["rpm"])
            else:
                n = 1 if c["ggt"] == "ax1" else 2
                if Tg - 273.15 > LIM["T_ax_max"] + 1: return None
                r = turb_axial(Tg, b, n, LIM["U_ax"], w_req=w_req)
                if r is None: return None
                r["D"] = 60 * r["U"] / (math.pi * st0["rpm"]); r["n"] = n
            if c["ggt"] == "radial" and Tg - 273.15 > LIM["T_rad_max"] + 1: return None
            tg.append(r); Tg = r["T2"]; pg /= r["ER"]
        # --- подмешивание утечки, переходный канал
        p_exit = P0 * 1.03 * (1 + dpg)
        if c["shaft"] == "single":
            ER = pg / p_exit
            st0 = sp[0]["st"][0]
            if c["ggt"] == "radial":
                Ufix = st0["U2"] * c.get("kD", 1.0) if c["rotor"] == "mono" else None
                if Tg - 273.15 > LIM["T_rad_max"] + 1: return None
                r = turb_radial_to_p(Tg, b, ER, U_fixed=Ufix, TIT_C=Tg - 273.15)
                if not r["ok"]: return None
            else:
                n = 2
                if Tg - 273.15 > LIM["T_ax_max"] + 1: return None
                r = turb_axial(Tg, b, n, LIM["U_ax"], ER=ER)
                if r is None: return None
                r["n"] = n
            r["D"] = 60 * r["U"] / (math.pi * st0["rpm"])
            tg.append(r)
            w_net = (r["w"] * mg - sp[0]["w"] / 0.99) * 0.98           # на 1 кг воздуха; 0.98 - мех. потери и собственные нужды
            T5 = r["T2"]; pt = None
        else:
            hm = (mg * h(Tg, b) + x_leak * h(T2)) / (mg + x_leak)
            b2 = b * mg / (mg + x_leak); Tg = T_h(hm, b2); pg *= (1 - dp_duct)
            ER = pg / p_exit
            if ER < 1.15: return None
            if c["pt"] == "radial":
                pt = turb_radial_to_p(Tg, b2, ER, TIT_C=Tg - 273.15); pt["n"] = 1
            else:
                n = 1 if c["pt"] == "ax1" else 2
                # окружная скорость СТ: под коэффициент нагрузки 1.6 на ступень, не выше 380 м/с
                R = Rg(b2); T2s = T_s(s0(Tg, b2) - R * math.log(ER), b2); dhs = h(Tg, b2) - h(T2s, b2)
                Upt = min(380.0, math.sqrt(0.9 * dhs / n / 1.6))
                pt = turb_axial(Tg, b2, n, Upt, ER=ER)
                if pt is None: return None
                pt["n"] = n
            w_net = pt["w"] * (mg + x_leak) * 0.98
            T5 = pt["T2"]
        if w_net <= 0: return None
        # --- рекуператор
        if eps > 0:
            if T5 - 273.15 > LIM["T_rec_max"]: return None
            T3n = T2 + eps * (T5 - T2)
            if abs(T3n - c.get("_T3", T2)) > 0.5:
                c["_T3"] = T3n; c["_m"] = m_core
                continue_flag = True
            else:
                continue_flag = False
        else:
            continue_flag = False
        eff = w_net / (f * LHV)
        m_tot = P / w_net; m_new = m_tot / N
        out = dict(eff=eff, w_net=w_net, m_tot=m_tot, m_core=m_new, T2=T2, T3=c.get("_T3", T2), T5=T5, f=f,
                   sp=sp, tg=tg, pt=pt, PR=PR, TIT_C=TIT_C, pg_pt=pg, Tg_pt=Tg, b=b, x_leak=x_leak, mg=mg)
        if abs(m_new - m_core) / m_core < 0.01 and not continue_flag: break
        m_core = m_new
    c.pop("_T3", None); c.pop("_m", None)
    return out

# ---------------------------------------------------------------- размеры, детали, литьё
def finalize(c, r, P=6.0e6):
    N = c["N"]; nsp, ncs = c["spools"], c["cstages"]
    m_core = r["m_core"]; flags = []; cast = "да"
    wheels_total = 0; uniq = 0; bladed_disc = 0
    Dmax_mono = 0; rotor_mass = 0.0
    for i in range(nsp):
        st = r["sp"][i]["st"]; t = r["tg"][nsp - 1 - i] if len(r["tg"]) > nsp - 1 - i else None
        Dc = st[0]["D2"]
        if c["rotor"] == "mono":
            Dt = t["D"] if t else Dc; D = max(Dc, Dt); Dmax_mono = max(Dmax_mono, D)
            wheels_total += N; uniq += 1; rotor_mass = max(rotor_mass, 1257 * D ** 3)
            if D > LIM["D_mono_max"]: return None
            if D > LIM["D_mono_ok"]: cast = "на пределе"
        else:
            wheels_total += N * ncs; uniq += ncs
            if Dc > LIM["D_ti_max"]: return None
            mass = sum(0.0825 * (4430 if c["comp"] == "ti" else 2800) * s["D2"] ** 3 for s in st)
            if t:
                if c["ggt"] == "radial":
                    wheels_total += N; uniq += 1
                    if t["D"] > LIM["D_rad_max"]: return None
                    if t["D"] > LIM["D_rad_ok"]: cast = "на пределе"
                    mass += 0.0825 * 8000 * t["D"] ** 3 * 1.2
                else:
                    n = t["n"]; wheels_total += N * n; uniq += n
                    if t["D"] * 1.15 > LIM["D_axcast"]: bladed_disc += n
                    mass += n * 0.03 * 8000 * (t["D"] * 1.15) ** 3 * 1.5
            rotor_mass = max(rotor_mass, mass)
    rotors = N * nsp
    pt = r["pt"]; pt_rpm = None; pt_D = None; gear = False
    if pt is not None:
        rotors += 1; n = pt["n"]; wheels_total += n; uniq += n
        rho = r["pg_pt"] / pt["ER"] ** 0.5 / (Rg(r["b"]) * 0.5 * (r["Tg_pt"] + pt["T2"]))
        mgas = r["m_tot"] * (r["mg"] + r["x_leak"])
        if c["drive"] == "direct":
            pt_rpm = 3000.0; pt_D = 60 * pt["U"] / (math.pi * pt_rpm)
            hb = mgas / (rho * 0.5 * pt["U"] * math.pi * pt_D)
            if hb / pt_D < 0.02 or pt_D > 2.4: return None
            bladed_disc += n
        else:
            A = mgas / (rho * 0.5 * pt["U"])
            if c["pt"] == "radial":
                pt_D = math.sqrt(mgas / (rho * pt["U"] * 0.25)); pt_rpm = 60 * pt["U"] / (math.pi * pt_D)
                if pt_D > LIM["D_rad_max"]: return None
                if pt_D > LIM["D_rad_ok"]: cast = "на пределе"
            else:
                pt_D = math.sqrt(A / (0.15 * math.pi)); pt_rpm = 60 * pt["U"] / (math.pi * pt_D)
                if pt_D * 1.15 > LIM["D_axcast"]: bladed_disc += n
            gear = c["use"] == "gen"
    else:
        gear = (N == 1)                      # одновальный одиночный - через редуктор; кластер - быстроходные генераторы
    hx = (1 if c["eps"] > 0 else 0) + (1 if c.get("ic") else 0)
    comb = 1
    # индекс сложности (чем меньше, тем проще) - формула в отчёте
    K = uniq + 0.5 * rotors + 2 * gear + 2 * hx + 2 * bladed_disc + (2 if (pt is None and N > 1) else 0)
    if pt is None and N > 1:
        eff_el = r["eff"] * 0.955
    elif gear:
        eff_el = r["eff"] * 0.985 * 0.972
    else:
        eff_el = r["eff"] * 0.972
    st0 = r["sp"][0]["st"][0]
    # масса рекуператора: NTU = eps/(1-eps), k = 120 Вт/м2К, лист 0.15 мм, коэффициент конструкции 3
    rec_t = 0.0
    if c["eps"] > 0:
        ntu = c["eps"] / (1 - c["eps"]); A = ntu * r["m_tot"] * 1100 / 120.0
        rec_t = A * 0.00015 * 7900 * 3 / 1000.0
    return dict(name=c["name"], fam=c["fam"], N=N, eps=c["eps"], PR=round(r["PR"], 2), TIT=r["TIT_C"],
                eff=round(r["eff"] * 100, 1), eff_el=round(eff_el * 100, 1), m=round(r["m_tot"], 1),
                T5=round(r["T5"] - 273.15), Tstack=round((r["T5"] - (r["T3"] - r["T2"]) * 0.93) - 273.15),
                rpm=round(st0["rpm"]), Dc=round(st0["D2"], 3), U2=round(st0["U2"]),
                Dt=round(r["tg"][0]["D"], 3) if r["tg"] else None,
                pt_rpm=round(pt_rpm) if pt_rpm else None, pt_D=round(pt_D, 2) if pt_D else None,
                wheels=wheels_total, uniq=uniq, rotors=rotors, gear=gear, hx=hx, bladed_disc=bladed_disc,
                cast=cast, K=round(K, 1), rotor_kg=round(rotor_mass), rec_t=round(rec_t, 1),
                use=c["use"], drive=c.get("drive"), eta_c=round(st0["eta"], 3),
                eta_t=round(r["tg"][0]["eta"], 3) if r["tg"] else None,
                eta_pt=round(pt["eta"], 3) if pt else None, np1=(N > 1))

def optimize(c):
    best = None
    prs = np.arange(2.5, 16.01, 0.25)
    tmax = LIM["T_rad_max"] if c["ggt"] == "radial" else LIM["T_ax_max"]
    for TIT in np.arange(800.0, tmax + 1, 25.0):
        for PR in prs:
            for kD in ((1.0, 1.1, 1.2) if c["rotor"] == "mono" else (1.0,)):
                c["kD"] = kD
                try:
                    r = evaluate(c, float(PR), float(TIT))
                except Exception:
                    r = None
                if r is None: continue
                fz = finalize(c, r)
                if fz is None: continue
                fz["kD"] = kD
                if best is None or fz["eff"] > best["eff"]: best = fz
    return best

# ---------------------------------------------------------------- генератор схем
def concepts():
    L = []
    def add(**k): L.append(k)
    # A. одновальные (только генератор): одиночные и кластер быстроходных турбогенераторов
    for N in (1, 4, 8):
        for comp, cst in (("ti", 1), ("ti", 2), ("al", 1)):
            for ggt in ("radial", "ax2"):
                for eps in (0.0, 0.85, 0.90):
                    add(fam="A", name=f"одновальный, {cst} центроб. ({comp}), турбина {'радиальная' if ggt=='radial' else 'осевая 2 ст.'}",
                        N=N, spools=1, cstages=cst, comp=comp, rotor="two", ggt=ggt, shaft="single", pt=None, eps=eps, use="gen")
        for eps in (0.0, 0.85, 0.90):
            add(fam="A", name="одновальный, цельнолитой монотор", N=N, spools=1, cstages=1, comp="mono", rotor="mono",
                ggt="radial", shaft="single", pt=None, eps=eps, use="gen")
    # B. двухвальные одиночные: газогенератор + свободная силовая турбина
    for comp, cst, dbl in (("ti", 1, False), ("ti", 1, True), ("ti", 2, False)):
        for ggt in ("radial", "ax1", "ax2"):
            for pt in ("ax1", "ax2"):
                for eps in (0.0, 0.85, 0.90):
                    add(fam="B", name=f"двухвальный, {cst} центроб.{' двусторонний' if dbl else ''}, ТГГ {'радиальная' if ggt=='radial' else ('осевая 1 ст.' if ggt=='ax1' else 'осевая 2 ст.')}, СТ {pt[-1]} ст.",
                        N=1, spools=1, cstages=cst, comp=comp, dbl=dbl, rotor="two", ggt=ggt, shaft="free", pt=pt, eps=eps,
                        use="both", drive="fast")
    # C. кластер турбо-газогенераторов + общая силовая турбина
    for N in (4, 6, 8, 12):
        for rotor, comp in (("mono", "mono"), ("two", "al"), ("two", "ti")):
            for eps in (0.0, 0.85, 0.90):
                for drive, pt in (("direct", "ax1"), ("direct", "ax2"), ("fast", "ax1")):
                    nm = {"mono": "цельнолитой монотор", "al": "алюм. колесо + литая турбина", "ti": "титан. колесо + литая турбина"}[comp]
                    add(fam="C", name=f"кластер {N} × ({nm}) + СТ {pt[-1]} ст. {'3000 об/мин' if drive=='direct' else 'быстроходная'}",
                        N=N, spools=1, cstages=1, comp=comp, rotor=rotor, ggt="radial", shaft="free", pt=pt, eps=eps,
                        use="gen" if drive == "direct" else "both", drive=drive)
    # E. два каскада последовательно (турбокомпаунд), с промежуточным охлаждением и без
    for N in (1, 4):
        for rotor, comp in (("mono", "mono"), ("two", "ti")):
            for ic in (False, True):
                for eps in (0.0, 0.85, 0.90):
                    add(fam="E", name=f"два каскада {'×'+str(N)+' ' if N>1 else ''}({'монотор' if rotor=='mono' else 'титан + литая турбина'}){', промохлаждение' if ic else ''} + СТ",
                        N=N, spools=2, cstages=1, comp=comp, rotor=rotor, ggt="radial", shaft="free", pt="ax1", eps=eps, ic=ic,
                        use="both", drive="fast")
    return L

if __name__ == "__main__":
    # --- сверка с прежней моделью (cycle.py): простой двухвальный, PR 6.5, 1000 C
    chk = dict(fam="chk", name="chk", N=1, spools=1, cstages=1, comp="ti", rotor="two", ggt="ax1", shaft="free", pt="ax1", eps=0.0, use="both", drive="fast")
    r = evaluate(chk, 6.5, 1000.0); print("сверка: PR 6.5, 1000 C ->", round(r["eff"] * 100, 1), "% ; уд. работа", round(r["w_net"] / 1e3), "кДж/кг; T5", round(r["T5"] - 273))
    chk["eps"] = 0.85
    r = evaluate(chk, 6.5, 1000.0); print("сверка с рекуператором ->", round(r["eff"] * 100, 1), "%")
    L = concepts(); print("схем:", len(L))
    res = []
    for i, c in enumerate(L):
        fz = optimize(c)
        if fz is None:
            res.append(dict(id=i + 1, name=c["name"], fam=c["fam"], N=c["N"], eps=c["eps"], feasible=False, use=c["use"]))
        else:
            fz["id"] = i + 1; fz["feasible"] = True; res.append(fz)
    json.dump(dict(**{"method": "0D-цикл Брайтона: свойства газа Cantera (GRI-Mech 3.0), метан, ISO; КПД узлов и ограничения - принятые, см. README", "software": "cantera " + __import__("cantera").__version__ + ", python3", "date": __import__("datetime").date.today().isoformat()}, lim=LIM, res=res), open(sys.argv[1] if len(sys.argv) > 1 else "gtu_explore_results.json", "w"), ensure_ascii=False, indent=1)
    ok = [x for x in res if x["feasible"]]
    print("выполнимых:", len(ok), "из", len(res))
    for x in sorted(ok, key=lambda x: -x["eff"])[:25]:
        print(f"{x['id']:3d} {x['fam']} N={x['N']:2d} eps={x['eps']:.2f} PR={x['PR']:5.2f} T={x['TIT']:.0f} КПД={x['eff']:.1f} K={x['K']:.1f} колёс={x['wheels']}/{x['uniq']} лит={x['cast']} | {x['name']}")
