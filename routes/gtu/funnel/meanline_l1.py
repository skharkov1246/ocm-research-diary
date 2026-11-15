"""Воронка двигателя-заменителя, уровень 1: одномерные модели ступеней по среднему радиусу (meanline). Всё - оценка.

Свойства газа - таблицы gtu_explore (Cantera GRI-Mech 3.0, метан) через funnel_space.FastGas.

  * Центробежная (и диагональная) ступень компрессора - набор потерь Oh, Yoon, Chung (1997, Proc. IMechE A 211:331-338):
    угол атаки (Conrad), нагрузка лопаток (Coppage), трение (Jansen), зазор (Jansen), смешение (Johnston & Dean),
    рециркуляция (Oh), трение диска (Daily & Nece), перетекание (Aungier). Формулы - те же, что ../frame3/frame3_point.losses
    (там сверены с открытым кодом TurboFlow 0.1.18, centrifugal_compressor/loss_model_oh.py; тест test_cc_losses_match_frame3
    сверяет совпадение); первоисточники в этой сессии не открывались. Добавлено (оценка): радиус выхода лопаточного диффузора
    по диаметру корпуса (r4/r2 = kD_case - 0,35), скорость на выходе диффузора c4/c2 = 0,45·1,55/(r4/r2) и невосстановленная
    половина её энергии (DUMP). Скольжение - Wiesner (1967), число лопаток эталона - Пфлейдерер (как frame3).
    Уровень КПД привязан к формуле трека gtu_explore.comp_stage (по ней калиброван цикл: SGT-400, Taurus 60, ±2 п. п.):
    сдвиг a + b·(π - 4,5), a и b - МНК по сетке расчётных точек comp_stage (μ 0,70, φ 0,085, cm2/U2 0,28). Причина линейного
    сдвига: набор Oh при φ = 0,085 почти не зависит от π (0,838-0,846 при π 2,5-7,5), формула трека теряет 1 п. на единицу π
    (обобщение трека; потерь в скачках у входа в наборе Oh нет, M1t доходит до 1,5) - наклон формулы переносится как есть,
    модель даёт зависимость от геометрии (φ, μ, cm2/U2, Z, b2, диффузор). Остаток сверки - заявленная точность TOL_CC.
    Вне сетки калибровки (μ 0,70, φ 0,085, cm2/U2 0,28) модель не подтверждена: КПД ступени ограничен сверху обобщённой
    формулой уровня 0 той же ступени (funnel_space.eta_comp + поправки derive: μ, φ, зазор/b2) плюс TOL_CC (исправление
    скептика уровня 1, 15.11.2026: оптимизатор уводил μ и cm2/U2 к нижним границам, где модель давала +1,4…+2,7 п. сверх
    comp_stage - экстраполяция, около половины её - добавки диффузора ∝ c2²); вниз модель штрафует плохую геометрию как есть.
    Поправки ступени (Ni −0,02, диагональная −0,01, двусторонний вход 1-й ступени −0,01 - как уровень 0) входят в сдвиг,
    поэтому геометрия (U2, D2) и нерасчётная модель видят тот же КПД, что цикл (исправление скептика: у колёс Ni режим ISO
    нерасчётной модели был при T41 на 26-27 K ниже расчётной).
    Нерасчётный режим той же геометрии - как frame3_point.stage_off (горла входа и диффузора: число Маха 0,8 в расчётной точке);
    закрутка на входе от ВНА (alpha1) - как там же; угол атаки на периферии входа - для признака срыва (offdesign_l1.surge).
  * Осевая турбина: набор потерь Kacker & Okapuu (1982) в записи открытого кода TurboFlow 0.1.18 (MIT;
    axial_turbine/loss_model_kacker_okapuu.py - открыт в этой сессии, колесо PyPI): профильные потери по подборам Aungier
    кривых Ainley - Mathieson (реактивная и активная решётки, смешение по θ_in/β_out), поправки на сжимаемость K_p, удар у втулки
    на входе, сверхзвуковой множитель 1 + 60(M - 1)², Рейнольдс; вторичные (K_s, удлинение, нагрузка Z), кромка (t_te/o).
    Кромка 0,8 мм, t_max/c 0,2 - оценка (литьё/печать). Радиальный зазор - Kacker & Okapuu: Δη/η = 0,93·(k/h)·(r_t/r_m)/cos β3
    (по памяти). Шаг - Цвейфель 0,8. До 15.11.2026 вместо профильных и вторичных потерь стоял Soderberg без членов на число Маха
    (soderberg() оставлен для сравнения; AX_MODEL): оптимизатор выводил сопла на M2 1,05-1,20, где Soderberg не видит роста
    потерь (замечание скептика уровня 1). Уровень привязан сдвигом AX_SHIFT к eta_axial(ψ = 1,2) = 0,885 трека при φ = 0,5 и
    удлинении H/b 1,3 (та же опорная ступень); потеря с выходной скоростью сверх φ = 0,5 - как уровень 0 (LEAVE_K).
    Свободный вихрь по высоте: у втулки - число Маха за соплом и реактивность, на периферии - относительное число Маха за
    ротором (правила funnel_l1.ml_rules; раньше проверялся только средний радиус).
  * Радиальная турбина: набор потерь вида Baines (2003; Moustapha et al., «Axial and Radial Turbines», 2003 - по памяти):
    сопловой аппарат (Soderberg), угол атаки ротора (оптимальный по скольжению Stanitz cθ4,opt = U4·(1 - 0,63π/Z), по Dixon & Hall
    - по памяти; до 15.11.2026 стояло 0,63π/(2Z) - потеря на угол атаки занижалась), канал ротора (K_p = 0,11), зазор
    (K_a 0,4, K_r 0,75, K_ar -0,3), трение диска (Daily & Nece), выходная скорость. Выходной радиус по быстроходности
    (r5t/r4 = 0,45 + 0,45·ns в пределах 0,55-0,78, оценка по смыслу кривых Rohlik, NASA TN D-4384 - отчёт по выдержке поиска).
    Уровень привязан сдвигом RAD_SHIFT к η_is = 0,86 при U/C0 = 0,70 и ns = 0,60 (gtu_explore.turb_radial). Многоступенчатая:
    полное давление на входе следующей ступени - статическое на выходе предыдущей плюс восстановленная доля (1 - LEAVE_K)
    энергии выходной скорости (раньше статическое шло как полное - энергия выхода терялась вся); степени расширения ступеней
    подбираются так, чтобы статическое давление за последней было заданным.
Сдвиги считаются при первом обращении (детерминированно) и пишутся в результаты; сырые значения моделей - тоже.
"""
import math, os, sys

HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in sys.path:
    sys.path.insert(0, HERE)
import funnel_space as FS        # noqa: E402

# =========================================================================== принятые величины (оценка; как frame3_point.A)
A = dict(nu=0.40, k1=0.97, B2=0.08, spl=0.70, eps_wake=0.366, zeta_d=0.15, cf_d=0.005, alpha4=40.0, f_inc=0.6, lax=0.8,
         m_thd=0.80, r4r2_ref=1.55, c4c2_ref=0.45)
K_DUMP = 0.5                 # невосстановленная доля энергии скорости на выходе диффузора (оценка)
INTERNAL = ("inc", "bld", "sf", "cl", "mix", "dif", "dfr", "dinc", "dump")
PARASITIC = ("rc", "df", "lk")
TOL_CC = 0.006               # заявленная точность сверки с comp_stage в его расчётных точках (доля КПД; остаток МНК - README)
PR_REF = 4.5                 # точка отсчёта линейного сдвига по π
ZW = 0.8                     # коэффициент Цвейфеля (оптимальный шаг)
LEAVE_K = FS.LEAVE_K
_SH = {}


def _E():
    return FS._E()


def cp(T, b=0.0):
    E = _E()
    return (E.h(T + 1.0, b) - E.h(T - 1.0, b)) / 2.0


def gam(T, b=0.0):
    c = cp(T, b)
    R = _E().Rg(b)
    return c / (c - R)


def mu_gas(T):
    """вязкость по Сазерленду (воздух; для продуктов сгорания - то же, оценка)"""
    return 1.716e-5 * (T / 273.15) ** 1.5 * (273.15 + 110.4) / (T + 110.4)


def k_gas(T):
    """теплопроводность по Сазерленду (воздух; оценка и для газа)"""
    return 0.0241 * (T / 273.15) ** 1.5 * (273.15 + 194.0) / (T + 194.0)


def stat(T0, p0, c, b=0.0):
    E = _E()
    R = E.Rg(b)
    T = E.T_h(E.h(T0, b) - 0.5 * c * c, b)
    p = p0 * math.exp((E.s0(T, b) - E.s0(T0, b)) / R)
    return T, p, p / (R * T)


def flux_M(T0, p0, M, b=0.0):
    g = gam(T0, b)
    R = _E().Rg(b)
    return p0 * math.sqrt(g / (R * T0)) * M * (1 + 0.5 * (g - 1) * M * M) ** (-(g + 1) / (2 * (g - 1)))


def inlet(T0, p0, m, area, alpha=0.0, b=0.0):
    """дозвуковое течение: скорость c при расходе m через площадь area (Ньютон по ρc - монотонно слева от корня).
    None - расход больше запирающего"""
    R = _E().Rg(b)
    G = m / (area * math.cos(alpha))
    if G >= 0.995 * flux_M(T0, p0, 1.0, b):
        return None
    c = G / (p0 / (R * T0))
    for _ in range(40):
        T, p, rho = stat(T0, p0, c, b)
        M2 = c * c / (gam(T, b) * R * T)
        if M2 >= 0.97:
            return None
        f = rho * c - G
        if abs(f) < 1e-10 * G:
            break
        c -= f / (rho * (1 - M2))
    T, p, rho = stat(T0, p0, c, b)
    return c, T, p, rho


# =========================================================================== центробежная ступень: потери (как frame3_point.losses)
def cc_losses(g, f, m_side, dbl, r4r2=None, dump=True):
    """Oh, Yoon, Chung (1997) + потеря на форму диффузора, трение стенок диффузора и угол атаки его лопаток (как frame3);
    r4r2 - выход диффузора / r2 (по умолчанию 1,55 - как frame3), dump - невосстановленная энергия выхода диффузора. Дж/кг"""
    r4r2 = A["r4r2_ref"] if r4r2 is None else r4r2
    c4c2 = A["c4c2_ref"] * A["r4r2_ref"] / r4r2
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
    Ld = (r4r2 - 1) * r2 / math.cos(am)
    L["dfr"] = 2 * A["cf_d"] * Ld / b2 * 0.5 * c2 ** 2 * (1 + c4c2 ** 2) / 2
    L["dinc"] = A["f_inc"] * 0.5 * (c2 * math.sin(a2 - g["a2v"])) ** 2
    L["dump"] = (K_DUMP * 0.5 * (c4c2 * c2) ** 2) if dump else 0.0
    L["rc"] = 8e-5 * math.sinh(3.5 * a2 ** 3) * Df ** 2 * U2 ** 2
    Re = rho2 * U2 * r2 / f["mu2"]
    fdf = 2.67 / Re ** 0.5 if Re < 3e5 else 0.0622 / Re ** 0.2
    L["df"] = 0.0 if dbl else fdf * 0.5 * (rho1 + rho2) * r2 ** 2 * U2 ** 3 / (4 * m_side)
    rbar, bbar = 0.5 * (r1 + r2), 0.5 * ((r1t - r1h) + b2)
    Lm = math.pi / 2 * math.sqrt(((r2 - r1) ** 2 + (g["Lax"] - 0.5 * b2) ** 2) / 2)
    dP = m_side * (r2 * cth2 - r1 * cth1) / (Z * rbar * bbar * Lm)
    ucl = 0.816 * math.sqrt(2 * max(dP, 0.0) / rho2)
    L["lk"] = rho2 * Z * g["cl"] * Lm * ucl * ucl * U2 / (2 * m_side)
    return L


def _finish(T01, L, dh_eu):
    E = _E()
    R = E.RA
    w0 = dh_eu + sum(L[k] for k in PARASITIC)
    dh_is = dh_eu - sum(L[k] for k in INTERNAL)
    T02 = E.T_h(E.h(T01) + w0)
    if dh_is <= 0:
        return dict(w0=w0, dh_is=dh_is, eta_p_model=0.3, PR_model=1.0, T02=T02, eta_is_model=0.0)
    PRm = math.exp((E.s0(E.T_h(E.h(T01) + dh_is)) - E.s0(T01)) / R)
    return dict(w0=w0, dh_is=dh_is, eta_is_model=dh_is / w0, PR_model=PRm, eta_p_model=R * math.log(PRm) / (E.s0(T02) - E.s0(T01)),
                T02=T02)


def pfleiderer(r1, r2, b1, b2):
    """число лопаток по Пфлейдереру; углы от меридионального направления (как frame3)"""
    return 6.5 * (r2 + r1) / (r2 - r1) * math.cos(0.5 * (b1 + b2))


def inducer(T01, p01, m_side, n, nu=None):
    """вход колеса по минимуму относительной скорости на периферии (как frame3): начальное приближение - несжимаемая
    формула d = (8a²/ω²)^(1/6), затем золотое сечение в [0,8; 1,35]·d0"""
    E = _E()
    nu = A["nu"] if nu is None else nu
    om = n * math.pi / 30
    k = A["k1"] * math.pi / 4 * (1 - nu * nu)
    a = m_side / (p01 / (E.RA * T01)) / k
    d0 = (8 * a * a / om ** 2) ** (1 / 6)

    def w1t(d):
        s = inlet(T01, p01, m_side, k * d * d)
        return 1e9 if s is None else math.hypot(s[0], om * d / 2)
    lo, hi = 0.8 * d0, 1.35 * d0
    gr = (math.sqrt(5) - 1) / 2
    x1, x2 = hi - gr * (hi - lo), lo + gr * (hi - lo)
    f1, f2 = w1t(x1), w1t(x2)
    for _ in range(16):
        if f1 < f2:
            hi, x2, f2 = x2, x1, f1
            x1 = hi - gr * (hi - lo); f1 = w1t(x1)
        else:
            lo, x1, f1 = x1, x2, f2
            x2 = lo + gr * (hi - lo); f2 = w1t(x2)
    d = 0.5 * (lo + hi)
    area = k * d * d
    s = inlet(T01, p01, m_side, area)
    if s is None:
        return None
    c, T, p, rho = s
    r1t, r1h = d / 2, nu * d / 2
    r1 = math.sqrt((r1t ** 2 + r1h ** 2) / 2)
    aa = math.sqrt(gam(T) * E.RA * T)
    return dict(D1t=d, D1h=nu * d, r1t=r1t, r1h=r1h, r1=r1, A1=area, c1=c, T1=T, p1=p, rho1=rho, U1t=om * r1t, a1=aa,
                W1t=math.hypot(c, om * r1t), M1t=math.hypot(c, om * r1t) / aa, M1=c / aa,
                b1t=math.atan(om * r1t / c), b1=math.atan(om * r1 / c), b1h=math.atan(om * r1h / c))


def _beta2(muE, kcm, Zeff):
    """угол лопатки на выходе от радиального по коэффициенту работы Эйлера: μE = σ(β) - kcm·tg β, σ - Wiesner.
    Возврат: (β, флаг): флаг "forward" - нужна загнутость вперёд (μE > σ(0)), "max" - загнутость больше 70°"""
    f = lambda b: 1 - math.sqrt(math.cos(b)) / Zeff ** 0.7 - kcm * math.tan(b) - muE
    lo, hi = 0.0, math.radians(70.0)
    if f(lo) < 0:
        return lo, "forward"
    if f(hi) > 0:
        return hi, "max"
    for _ in range(40):
        mid = 0.5 * (lo + hi)
        if f(mid) > 0:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi), ""


def eta_cap_l0(PRs, mu, phi, b2, D2, typ="Ц"):
    """верхняя граница КПД ступени без поправок материала, типа и входа: обобщённая формула уровня 0 той же ступени
    (funnel_space.eta_comp: comp_stage трека 0,880 - 0,010·(π - 2,5) и MU_ETA·(μ - 0,70); funnel_space.derive: зазор
    0,3·cl/b2 и PHI_PENALTY·ln²(φ/φ_opt)) плюс TOL_CC - то, что подтверждено калибровкой трека; сверх неё модель не засчитывается"""
    cl = max(0.3e-3, 0.0007 * D2)
    e = 0.880 - 0.010 * (PRs - 2.5) - FS.MU_ETA * (mu - 0.70) - 0.3 * cl / b2
    e -= FS.PHI_PENALTY * math.log(max(phi, 1e-6) / FS.PHI_OPT.get(typ, FS.PHI_OPT["Ц"])) ** 2
    return e + TOL_CC


def cc_stage(T01, p01, m_side, PRs, n, mu, kcm, Z, dbl, kD=1.9, eta0=0.84, shift=None, nu=None, ref=False, d_pen=0.0, cap=True,
             typ="Ц"):
    """центробежная ступень: заданы π, частота n, коэффициент работы μ = w/U2² (как уровень 0), cm2/U2, число основных
    лопаток Z (+ разделительные: Zeff = Z·(1+spl)), диаметр корпуса / D2 (диффузор). Итерации: КПД -> работа -> U2 -> угол
    лопатки по скольжению -> треугольники -> потери -> КПД модели (+ сдвиг CC_SHIFT). ref=True - число лопаток по Пфлейдереру.
    cap - КПД не выше eta_cap_l0 (вне сетки калибровки модель не подтверждена); d_pen - поправки ступени (Ni, диагональная,
    двусторонний вход), вычитаются из КПД внутри итераций: U2, D2 и нерасчётная модель (сдвиг shift в ответе) видят тот же КПД"""
    E = _E()
    R = E.RA
    shift = cc_shift(PRs) if shift is None else shift
    om = n * math.pi / 30
    ind = inducer(T01, p01, m_side, n, nu)
    if ind is None:
        return dict(ok=False, why="вход колеса заперт")
    r4r2 = min(1.80, max(1.30, kD - 0.35))
    eta, par, p02i = eta0, None, p01 * PRs
    flags = set()
    zf = Z
    L = f = g = None
    for it in range(40):
        T02 = E.T_s(E.s0(T01) + R * math.log(PRs) / eta)
        w0 = E.h(T02) - E.h(T01)
        U2 = math.sqrt(w0 / mu)
        r2 = U2 / om
        if r2 <= ind["r1t"] * 1.08:
            return dict(ok=False, why=f"колесо меньше входа: D1t/D2 {ind['r1t'] / r2:.2f}", ind=ind)
        if par is None:
            par = 0.02 * w0
        dh_eu = w0 - par
        muE = dh_eu / U2 ** 2
        zeff = zf * (1 + A["spl"])
        b2b, fl = _beta2(muE, kcm, zeff)
        if ref:
            zn = max(6.0, pfleiderer(ind["r1"], r2, ind["b1"], b2b))
            if abs(zn - zf) > 1e-6:
                zf = zn
                zeff = zf * (1 + A["spl"])
                b2b, fl = _beta2(muE, kcm, zeff)
        sig = 1 - math.sqrt(math.cos(b2b)) / zeff ** 0.7
        cm2 = kcm * U2
        cth2 = muE * U2
        c2, W2, a2 = math.hypot(cm2, cth2), math.hypot(cm2, cth2 - U2), math.atan2(cth2, cm2)
        T2 = E.T_h(E.h(T02) - 0.5 * c2 * c2)
        p2 = p02i * math.exp((E.s0(T2) - E.s0(T02)) / R)
        rho2 = p2 / (R * T2)
        b2 = m_side / (rho2 * cm2 * 2 * math.pi * r2 * (1 - A["B2"]))
        g = dict(r1t=ind["r1t"], r1h=ind["r1h"], r1=ind["r1"], r2=r2, b2=b2, Zfull=zf, Zeff=zeff, b1b=ind["b1"], b2b=b2b,
                 cl=max(0.3e-3, 0.0007 * 2 * r2), Lax=A["lax"] * (r2 - ind["r1"]) + b2, a2v=a2, sig=sig, A1=ind["A1"], dbl=dbl,
                 r4r2=r4r2)
        f = dict(U2=U2, cm2=cm2, cth2=cth2, c2=c2, W2=W2, alpha2=a2, rho1=ind["rho1"], rho2=rho2, mu1=mu_gas(ind["T1"]),
                 mu2=mu_gas(T2), W1t=ind["W1t"], W1h=math.hypot(ind["c1"], om * ind["r1h"]), W1=math.hypot(ind["c1"], om * ind["r1"]),
                 c1=ind["c1"], c1m=ind["c1"], cth1=0.0, inc=0.0, dh_eu=dh_eu, T2=T2)
        L = cc_losses(g, f, m_side, dbl, r4r2)
        par_new = sum(L[k] for k in PARASITIC)
        dh_imp = dh_eu - (L["inc"] + L["bld"] + L["sf"] + L["cl"])
        if dh_imp <= 0:
            return dict(ok=False, why="потери колеса больше работы", ind=ind)
        p02i_new = p01 * math.exp((E.s0(E.T_h(E.h(T01) + dh_imp)) - E.s0(T01)) / R)
        res = _finish(T01, L, dh_eu)
        e_aero = res["eta_p_model"] + shift
        capped = False
        if cap:
            ecap = eta_cap_l0(PRs, mu, m_side / (p01 / (R * T01) * U2 * (2 * r2) ** 2), b2, 2 * r2, typ)
            if e_aero > ecap:
                e_aero, capped = ecap, True
        eta_new = min(0.95, max(0.55, e_aero - d_pen))
        done = abs(eta_new - eta) < 2e-5 and abs(par_new - par) < 0.05 and abs(p02i_new - p02i) < 1.0
        eta = 0.5 * eta + 0.5 * eta_new if it > 15 else eta_new
        par, p02i = par_new, p02i_new
        if fl:
            flags.add(fl)
        if done:
            break
    res = _finish(T01, L, f["dh_eu"])
    T02 = E.T_s(E.s0(T01) + R * math.log(PRs) / eta)
    a2s = math.sqrt(gam(f["T2"]) * R * f["T2"])
    # горла входа колеса (относительное движение, средний радиус) и диффузора - число Маха A["m_thd"] в расчётной точке
    W1 = f["W1"]
    T0r = E.T_h(E.h(ind["T1"]) + 0.5 * W1 * W1)
    p0r = ind["p1"] * math.exp((E.s0(T0r) - E.s0(ind["T1"])) / R)
    A_th = m_side / flux_M(T0r, p0r, A["m_thd"])
    A_thd = m_side / flux_M(T02, p02i, A["m_thd"])
    rho01 = p01 / (R * T01)
    D2 = 2 * g["r2"]
    phi_s = m_side / (rho01 * f["U2"] * D2 ** 2)
    ecap = eta_cap_l0(PRs, mu, phi_s, g["b2"], D2, typ) if cap else None
    return dict(ok=True, why=None, flags=sorted(flags), eta_p=eta, eta_p_model=res["eta_p_model"], eta_is_model=res["eta_is_model"],
                shift=eta - res["eta_p_model"], shift_cal=shift, d_pen=d_pen, capped=bool(ecap is not None and res["eta_p_model"] + shift > ecap),
                eta_cap=ecap, eta_uncapped=res["eta_p_model"] + shift - d_pen,
                PRs=PRs, w0=E.h(T02) - E.h(T01), T01=T01, p01=p01, T02=T02, m_side=m_side, dbl=dbl, n=n,
                U2=f["U2"], D2=D2, b2=g["b2"], beta2b_deg=math.degrees(g["b2b"]), alpha2_deg=math.degrees(f["alpha2"]),
                c2=f["c2"], M2=f["c2"] / a2s, Mu2=f["U2"] / math.sqrt(gam(T01) * R * T01), M1t=ind["M1t"], M1=ind["M1"],
                W2_W1t=f["W2"] / ind["W1t"], Df=L["Df"], phi=m_side / (rho01 * f["U2"] * D2 ** 2), mu=mu, kcm=kcm,
                Zfull=g["Zfull"], Zeff=g["Zeff"], D1t=ind["D1t"], D1h=ind["D1h"], D1t_D2=ind["D1t"] / D2, Lax=g["Lax"], r4r2=r4r2,
                D4=r4r2 * D2, cl=g["cl"], losses_pct={k: 100 * v / max(res["w0"], 1.0) for k, v in L.items() if k != "Df"},
                g=g, f=f, ind=ind, A_th=A_th, A_thd=A_thd, rho2=f["rho2"], T2=f["T2"], p02i=p02i)


def cc_off(st, T01, p01, m_side, n, alpha1=0.0, iters=60):
    """нерасчётный режим той же геометрии (как frame3_point.stage_off): угол атаки на входе колеса и лопаток диффузора,
    запирание горла входа и диффузора; КПД = модель + тот же сдвиг"""
    E = _E()
    R = E.RA
    g = st["g"]
    dbl = g["dbl"]
    om = n * math.pi / 30
    s = inlet(T01, p01, m_side, g["A1"], alpha1)
    if s is None:
        return dict(ok=False, why="вход колеса заперт")
    c1, T1, p1, rho1 = s
    c1m, cth1 = c1 * math.cos(alpha1), c1 * math.sin(alpha1)
    U1, U1t, U1h = om * g["r1"], om * g["r1t"], om * g["r1h"]
    W1 = math.hypot(c1m, U1 - cth1)
    inc = math.atan2(U1 - cth1, c1m) - g["b1b"]
    b1t_d = st["ind"]["b1t"] if st.get("ind") else math.atan(U1t / max(c1m, 1e-6))
    inc_t = math.atan2(U1t - cth1 * g["r1"] / g["r1t"], c1m) - b1t_d          # периферия: закрутка ВНА - свободный вихрь (r·cθ)
    W1t = math.hypot(c1m, U1t - cth1 * g["r1"] / g["r1t"])
    W1h = math.hypot(c1m, U1h - cth1 * g["r1"] / g["r1h"])
    T0r = E.T_h(E.h(T1) + 0.5 * W1 * W1)
    p0r = p1 * math.exp((E.s0(T0r) - E.s0(T1)) / R)
    choke_i = m_side / (st["A_th"] * flux_M(T0r, p0r, 1.0))
    if choke_i >= 1.0:
        return dict(ok=False, why="горло входа колеса заперто", choke_i=choke_i)
    U2, r2, b2, b2b = om * g["r2"], g["r2"], g["b2"], g["b2b"]
    rho2, p02i, T2 = st["rho2"] * (p01 / st["p01"]), st["p02i"] * (p01 / st["p01"]), st["T2"]
    L = None
    for _ in range(iters):
        cm2 = m_side / (rho2 * 2 * math.pi * r2 * b2 * (1 - A["B2"]))
        cth2 = g["sig"] * U2 - cm2 * math.tan(b2b)
        dh_eu = U2 * cth2 - U1 * cth1
        if dh_eu <= 0:
            return dict(ok=False, why="нет работы")
        c2, W2, a2 = math.hypot(cm2, cth2), math.hypot(cm2, cth2 - U2), math.atan2(cth2, cm2)
        f = dict(U2=U2, cm2=cm2, cth2=cth2, c2=c2, W2=W2, alpha2=a2, rho1=rho1, rho2=rho2, mu1=mu_gas(T1), W1t=W1t, W1h=W1h,
                 W1=W1, c1=c1, c1m=c1m, cth1=cth1, inc=inc, dh_eu=dh_eu, mu2=mu_gas(T2))
        L = cc_losses(g, f, m_side, dbl, g["r4r2"])
        w0 = dh_eu + sum(L[k] for k in PARASITIC)
        T02 = E.T_h(E.h(T01) + w0)
        T2 = E.T_h(E.h(T02) - 0.5 * c2 * c2)
        dh_imp = dh_eu - (L["inc"] + L["bld"] + L["sf"] + L["cl"])
        if dh_imp <= 0:
            return dict(ok=False, why="потери больше работы")
        p02i_new = p01 * math.exp((E.s0(E.T_h(E.h(T01) + dh_imp)) - E.s0(T01)) / R)
        rho2_new = p02i_new * math.exp((E.s0(T2) - E.s0(T02)) / R) / (R * T2)
        done = abs(rho2_new - rho2) < 1e-7 * rho2 and abs(p02i_new - p02i) < 1e-7 * p02i
        rho2 = 0.5 * (rho2 + rho2_new)
        p02i = p02i_new
        if done:
            break
    res = _finish(T01, L, dh_eu)
    eta_p = max(0.3, res["eta_p_model"] + st["shift"])
    PRs = math.exp(eta_p * (E.s0(res["T02"]) - E.s0(T01)) / R)
    choke_d = m_side / (st["A_thd"] * flux_M(res["T02"], p02i, 1.0))
    return dict(ok=choke_d < 1.0, why=None if choke_d < 1.0 else "горло диффузора заперто", PRs=PRs, eta_p=eta_p, w0=res["w0"],
                T02=res["T02"], p02=p01 * PRs, inc_deg=math.degrees(inc), inc_tip_deg=math.degrees(inc_t), dinc_deg=math.degrees(a2 - g["a2v"]),
                choke_d=choke_d, choke_i=choke_i, alpha2_deg=math.degrees(a2))


def comp_stage_ref(T01, p01, PRs, m, dbl=False):
    """формула трека gtu_explore.comp_stage (титан) на быстрых таблицах: КПД, U2, D2, частота"""
    E = _E()
    R = E.RA
    eta = 0.880 - 0.010 * (PRs - 2.5)
    for _ in range(3):
        T2 = E.T_s(E.s0(T01) + R * math.log(PRs) / eta); w = E.h(T2) - E.h(T01)
        U2 = math.sqrt(w / 0.70)
        rho1 = p01 / (R * T01); mm = m / 2 if dbl else m
        D2 = math.sqrt(mm / (rho1 * U2 * 0.085))
        rho2 = rho1 * PRs ** 0.62 / ((T01 + 0.6 * (T2 - T01)) / T01)
        b2 = mm / (math.pi * D2 * rho2 * 0.28 * U2)
        cl = max(0.3e-3, 0.0007 * D2)
        eta_adj = 0.880 - 0.010 * (PRs - 2.5) - 0.3 * cl / b2
        if abs(eta_adj - eta) < 1e-4:
            break
        eta = eta_adj
    return dict(eta=eta, U2=U2, D2=D2, rpm=60 * U2 / (math.pi * D2), b2=b2)


CAL_GRID = [(pr, m, dbl) for pr in (2.5, 3.5, 4.5, 6.0, 7.5) for m, dbl in ((6.0, False), (20.0, False), (40.0, True))]


def cc_calibration():
    """сдвиг a + b·(π - PR_REF): МНК (comp_stage - модель) по сетке CAL_GRID (μ 0,70, φ 0,085, cm2/U2 0,28, лопатки по
    Пфлейдереру, корпус 1,9·D2); остаток после сдвига - сверка (заявлено TOL_CC)"""
    E = _E()
    rows = []
    for pr, m, dbl in CAL_GRID:
        cs = comp_stage_ref(E.T0, E.P0, pr, m, dbl)
        ms = m / 2 if dbl else m
        st = cc_stage(E.T0, E.P0, ms, pr, cs["rpm"], 0.70, 0.28, 10.0, dbl, kD=1.9, eta0=cs["eta"], shift=0.0, ref=True, cap=False)
        rows.append(dict(PR=pr, m=m, dbl=dbl, eta_cs=cs["eta"], eta_model=st["eta_p"] if st["ok"] else None,
                         phi_model=st.get("phi"), M1t=st.get("M1t"), Z_pfl=st.get("Zfull"), D2_cs=cs["D2"], D2_model=st.get("D2")))
    ok = [r for r in rows if r["eta_model"] is not None]
    xs = [r["PR"] - PR_REF for r in ok]; ys = [r["eta_cs"] - r["eta_model"] for r in ok]
    n = len(ok); mx = sum(xs) / n; my = sum(ys) / n
    b = sum((x - mx) * (y - my) for x, y in zip(xs, ys)) / sum((x - mx) ** 2 for x in xs)
    a = my - b * mx
    for r in rows:
        if r["eta_model"] is not None:
            r["d_after_shift"] = r["eta_model"] + a + b * (r["PR"] - PR_REF) - r["eta_cs"]
    return dict(a=a, b=b, rows=rows, max_abs_after=max(abs(r["d_after_shift"]) for r in ok), tol=TOL_CC)


def cc_shift(PRs=PR_REF):
    if "cc" not in _SH:
        _SH["cc_cal"] = cc_calibration()
        _SH["cc"] = (_SH["cc_cal"]["a"], _SH["cc_cal"]["b"])
    a, b = _SH["cc"]
    return a + b * (PRs - PR_REF)


# =========================================================================== осевая турбина (Kacker - Okapuu по TurboFlow 0.1.18)
AX_MODEL = "ko"      # "ko" - Kacker & Okapuu (1982) в записи TurboFlow; "soderberg" - прежняя модель (для сравнения)
T_TE = 0.8e-3        # толщина выходной кромки, м (оценка: литьё / печать)
TMAX_C = 0.20        # наибольшая толщина / хорда (оценка; значение по умолчанию формулы)
HUB_M2_MAX = 1.20    # число Маха за соплом у втулки (оценка: тот же предел, что на среднем радиусе, frame3_point.A m2_max)
R_HUB_MIN = 0.0      # реактивность у втулки (оценка: отрицательная реактивность у втулки - отрыв на входе ротора; по памяти, Dixon)


def soderberg(eps_deg, H, b_ax, Re, rotor):
    """коэффициент потерь энергии ряда: ζ* (Soderberg) с поправками на удлинение H/b (b - осевая хорда) и Рейнольдса"""
    z = 0.04 + 0.06 * (eps_deg / 100.0) ** 2
    ar = (0.975 + 0.075 * b_ax / H) if rotor else (0.993 + 0.021 * b_ax / H)
    z1 = (1 + z) * ar - 1
    return z1 * (1e5 / min(1e6, max(5e4, Re))) ** 0.25


def aungier_reaction(sc, ang):
    """Y_p реактивной (сопловой) решётки: подбор Aungier кривых Ainley - Mathieson (TurboFlow nozzle_blades); ang - выход от оси, град"""
    phi = 90.0 - ang
    scm = (0.46 + phi / 77) if phi < 30 else (0.614 + phi / 130)
    X = sc - scm
    A = (0.025 + (27 - phi) / 530) if phi < 27 else (0.025 + (27 - phi) / 3085)
    B = 0.1583 - phi / 1640
    C = 0.08 * ((phi / 30) ** 2 - 1)
    n = 1 + phi / 30
    return (A + B * X ** 2 + C * X ** 3) if phi < 30 else (A + B * abs(X) ** n)


def aungier_impulse(sc, ang):
    """Y_p активной решётки (TurboFlow impulse_blades)"""
    phi = 90.0 - ang
    scm = 0.224 + 1.575 * (phi / 90) - (phi / 90) ** 2
    X = sc - scm
    A = 0.242 - phi / 151 + (phi / 127) ** 2
    B = (0.3 + (30 - phi) / 50) if phi < 30 else (0.3 + (30 - phi) / 275)
    C = 0.88 - phi / 42.4 + (phi / 72.8) ** 2
    return A + B * X ** 2 - C * X ** 3


def _interp(x, xs, ys):
    x = min(xs[-1], max(xs[0], x))
    for k in range(len(xs) - 1):
        if x <= xs[k + 1]:
            return ys[k] + (ys[k + 1] - ys[k]) * (x - xs[k]) / (xs[k + 1] - xs[k])
    return ys[-1]


def ko_Y(th_in, b_out, M_in, M_out, s, c, bx, H, htr_in, Re_c, rotor, dyn_ratio, t_te=T_TE, tmax_c=TMAX_C):
    """коэффициент потерь полного давления ряда Y = Δp0/(p0 - p) на выходе: Kacker & Okapuu (1982) в записи TurboFlow 0.1.18
    (loss_model_kacker_okapuu: profile, secondary, trailing edge; зазор - отдельно). Углы - от оси, град; у ротора b_out < 0
    (поворот = th_in - b_out); dyn_ratio - отношение динамических напоров вход/выход (удар у втулки на входе)"""
    K1 = 1.0 if M_out < 0.2 else (1.0 - 1.25 * (M_out - 0.2) if M_out < 1.0 else 0.0)
    K2 = (M_in / max(M_out, 1e-6)) ** 2
    Kp = max(0.1, 1.0 - K2 * (1.0 - K1))
    ang = max(abs(b_out), 40.0)
    Yr, Yi = aungier_reaction(s / c, ang), aungier_impulse(s / c, ang)
    q = th_in / b_out
    Yp = Yr - abs(q) * q * (Yi - Yr)
    Yp = max(Yp, 0.8 * Yr)
    Yp *= (tmax_c / 0.2) ** max(0.0, -q)
    f_hub = _interp(htr_in, [0.5, 0.6, 0.7, 0.8, 0.9, 1.0],
                    [2.15, 1.7, 1.35, 1.12, 1.0, 1.0] if rotor else [1.4, 1.18, 1.05, 1.0, 1.0, 1.0])
    a = max(0.0, f_hub * M_in - 0.4)
    Ysh = max(0.0, 0.75 * a ** 1.75 * max(htr_in, 0.0) * dyn_ratio)
    Yp = 0.914 * (2.0 / 3.0 * Yp * Kp + Ysh)
    fRe = (Re_c / 2e5) ** (-0.4) if Re_c < 2e5 else (1.0 if Re_c <= 1e6 else (Re_c / 1e6) ** (-0.2))
    fMa = 1.0 + 60.0 * (M_out - 1.0) ** 2 if M_out > 1.0 else 1.0
    Yp *= fRe * fMa
    t_i, t_o = math.tan(math.radians(th_in)), math.tan(math.radians(b_out))
    am = math.atan(0.5 * (t_i + t_o))
    Zl = 4.0 * (t_i - t_o) ** 2 * math.cos(math.radians(b_out)) ** 2 / math.cos(am)
    hc = H / c
    far = (1.0 - 0.25 * math.sqrt(abs(2.0 - hc))) / hc if hc < 2.0 else 1.0 / hc
    Ks = max(0.1, 1.0 - min(1.0, (bx / H) ** 2) * (1.0 - Kp))   # K3 = (b/H)² не больше 1: при H/b < 1 поправка выходит из области
    #                                                             применимости и гасит вторичные потери до предела 0,1 (оценка)
    Ys = 1.2 * Ks * 0.0334 * far * Zl * math.cos(math.radians(b_out)) / math.cos(math.radians(th_in))
    o = s * math.cos(math.radians(b_out))
    rto = min(0.4, t_te / max(o, 1e-6))
    dr, di = _interp(rto, [0.0, 0.2, 0.4], [0.0, 0.045, 0.15]), _interp(rto, [0.0, 0.2, 0.4], [0.0, 0.025, 0.075])
    dphi = max(dr - abs(q) * q * (di - dr), di / 2.0)
    Yte = 1.0 / (1.0 - dphi) - 1.0
    return dict(Y=Yp + Ys + Yte, Yp=Yp, Ys=Ys, Yte=Yte, fMa=fMa, Kp=Kp, Ysh=Ysh)


def y_to_zeta(Y, M, g):
    """Y = (p01 - p02)/(p02 - p2) -> ζ = T·Δs / (c²/2) = ln(p01/p02) / (γM²/2)"""
    pr = (1.0 + 0.5 * (g - 1.0) * M * M) ** (g / (g - 1.0))
    return math.log(1.0 + Y * (1.0 - 1.0 / pr)) / (0.5 * g * M * M)


def _row(s, a_in, a_out, H, rho, vel, T, b=0.0, rotor=False, M_in=None, M_out=None, htr_in=0.8, dyn_ratio=0.0):
    """ряд: шаг s, углы входа/выхода (рад, знак: поворот = a_in + a_out), высота H, плотность и скорость на выходе.
    Осевая хорда - Цвейфель ZW; ζ - Kacker - Okapuu (AX_MODEL "ko", нужны M_in, M_out) или Soderberg; удлинение H/b"""
    eps = math.degrees(a_in + a_out)
    bx = 2 * s * math.cos(a_out) ** 2 * (math.tan(a_in) + math.tan(a_out)) / ZW
    bx = max(bx, 0.15 * s)
    stag = math.atan(0.5 * abs(math.tan(a_out) - math.tan(a_in)))
    chord = bx / max(0.2, math.cos(stag))
    out = dict(eps_deg=eps, b_ax=bx, chord=chord, AR=H / bx, pitch=s)
    if AX_MODEL == "ko" and M_out is not None:
        Re = rho * vel * chord / mu_gas(T)
        y = ko_Y(math.degrees(a_in), -math.degrees(a_out) if rotor else math.degrees(a_out), M_in, M_out, s, chord, bx, H, htr_in, Re,
                 rotor, dyn_ratio)
        out.update(zeta=y_to_zeta(y["Y"], M_out, gam(T, b)), Re=Re, Y=y["Y"], Yp=y["Yp"], Ys=y["Ys"], Yte=y["Yte"], fMa=y["fMa"])
    else:
        Dh = 2 * s * H * math.cos(a_out) / (s * math.cos(a_out) + H)
        Re = rho * vel * Dh / mu_gas(T)
        out.update(zeta=soderberg(eps, H, bx, Re, rotor), Re=Re)
    return out


def ax_turbine(T0in, p0in, b, m, U, n, w_tot, nst, phi, zn, zt, k_cl=None, shift=None):
    """осевая турбина: nst одинаковых ступеней на среднем диаметре Dm = 60U/(πn), ca = φ·U, выход без закрутки
    (реактивность 1 - ψ/2). Возврат: политропный КПД (с привязкой AX_SHIFT и потерей с выходной скоростью сверх φ = 0,5 как
    уровень 0), ступени (углы, числа Маха на среднем радиусе, у втулки и на периферии по свободному вихрю, высоты, хорды,
    удлинения, потери), сырой КПД модели"""
    E = _E()
    Rg = E.Rg(b)
    shift = ax_shift() if shift is None else shift
    Dm = 60 * U / (math.pi * n)
    rm = Dm / 2
    dh = w_tot / nst
    psi = dh / U ** 2
    ca = phi * U
    a2 = math.atan(psi / phi)                   # выход соплового аппарата
    b2 = math.atan((psi - 1) / phi)             # вход ротора (относительный), знак: >0 - со стороны набегания
    b3 = math.atan(1 / phi)                     # выход ротора (модуль)
    c2, w2, w3 = math.hypot(ca, psi * U), math.hypot(ca, (psi - 1) * U), math.hypot(ca, U)
    T01, p01 = T0in, p0in
    out = []
    for k in range(nst):
        h01 = E.h(T01, b)
        T1 = E.T_h(h01 - 0.5 * ca * ca, b)
        p1 = p01 * math.exp((E.s0(T1, b) - E.s0(T01, b)) / Rg)
        rho1 = p1 / (Rg * T1)
        T2 = E.T_h(h01 - 0.5 * c2 * c2, b)
        p2 = p01 * math.exp((E.s0(T2, b) - E.s0(T01, b)) / Rg)
        rho2 = p2 / (Rg * T2)
        H2 = m / (rho2 * ca * math.pi * Dm)
        T03 = E.T_h(h01 - dh, b)
        T3 = E.T_h(h01 - dh - 0.5 * ca * ca, b)
        aa1, aa2, aa3 = (math.sqrt(gam(T, b) * Rg * T) for T in (T1, T2, T3))
        M1, M2, Mw2, Mw3 = ca / aa1, c2 / aa2, w2 / aa2, w3 / aa3
        H1 = m / (rho1 * ca * math.pi * Dm)
        rowN = _row(math.pi * Dm / zn, 0.0, a2, H2, rho2, c2, T2, b, False, M_in=M1, M_out=M2,
                    htr_in=(Dm - H1) / (Dm + H1), dyn_ratio=(rho1 * ca * ca) / (rho2 * c2 * c2))
        eta = 0.9
        for _ in range(3):
            dhs = dh / eta
            T03s = E.T_h(h01 - dhs, b)
            p03 = p01 * math.exp((E.s0(T03s, b) - E.s0(T01, b)) / Rg)
            p3 = p03 * math.exp((E.s0(T3, b) - E.s0(T03, b)) / Rg)
            rho3 = p3 / (Rg * T3)
            H3 = m / (rho3 * ca * math.pi * Dm)
            Hr = 0.5 * (H2 + H3)
            rowR = _row(math.pi * Dm / zt, b2, b3, Hr, rho3, w3, T3, b, True, M_in=Mw2, M_out=Mw3,
                        htr_in=(Dm - H2) / (Dm + H2), dyn_ratio=(rho2 * w2 * w2) / (rho3 * w3 * w3))
            loss = rowN["zeta"] * 0.5 * c2 * c2 + rowR["zeta"] * 0.5 * w3 * w3
            eta0 = dh / (dh + loss)
            kc = k_cl if k_cl is not None else max(0.4e-3, 0.001 * (Dm + H3))
            rt_rm = (Dm + H3) / Dm
            d_tc = 0.93 * (kc / H3) * rt_rm / math.cos(b3)
            eta = eta0 * (1 - d_tc)
        cpg = cp(T2, b)
        # свободный вихрь (r·cθ = const, ca постоянна по высоте): втулка за соплом, втулка ротора, периферия выхода ротора
        rh2 = max(0.05 * rm, rm - H2 / 2)
        c2h = math.hypot(ca, psi * U * rm / rh2)
        T2h = E.T_h(h01 - 0.5 * c2h * c2h, b)
        M2h = c2h / math.sqrt(gam(T2h, b) * Rg * T2h)
        rhr = max(0.05 * rm, rm - Hr / 2)
        R_hub = 1 - (psi / 2) * (rm / rhr) ** 2
        rt3 = rm + H3 / 2
        Mw3t = math.hypot(ca, U * rt3 / rm) / aa3
        out.append(dict(stage=k + 1, psi=psi, phi=phi, reaction=1 - psi / 2, alpha2_deg=math.degrees(a2), beta2_deg=math.degrees(b2),
                        beta3_deg=math.degrees(b3), M1=M1, M2=M2, Mw3=Mw3, Mw2=Mw2, H2=H2, H3=H3, Dm=Dm, Dtip=Dm + H3,
                        htr=(Dm - H3) / (Dm + H3), AN2=math.pi * Dm * H3 * n * n, zeta_N=rowN["zeta"], zeta_R=rowR["zeta"],
                        Y_N=rowN.get("Y"), Y_R=rowR.get("Y"), fMa_N=rowN.get("fMa"), fMa_R=rowR.get("fMa"),
                        AR_N=rowN["AR"], AR_R=rowR["AR"], bx_N=rowN["b_ax"], bx_R=rowR["b_ax"], chord_R=rowR["chord"],
                        Re_R=rowR["Re"], k_cl=kc, d_eta_tc=d_tc * eta0, eta_tt=eta, M2_hub=M2h, reaction_hub=R_hub, Mw3_tip=Mw3t,
                        T_rel_in=T2 + 0.5 * w2 * w2 / cpg, T01=T01, T03=T03))
        T01, p01 = T03, p03
    eta_p_model = (E.s0(T0in, b) - E.s0(T01, b)) / (Rg * math.log(p0in / p01))
    eta_p = eta_p_model + shift - FS.leaving_loss(phi, U, w_tot / max(0.5, eta_p_model + shift))
    L = sum(s["bx_N"] + s["bx_R"] for s in out) * 1.35       # осевая длина: хорды + зазоры 35 % (оценка)
    return dict(eta_p=eta_p, eta_p_model=eta_p_model, shift=shift, model=AX_MODEL, psi=psi, phi=phi, Dm=Dm, U=U, stages=out, ER=p0in / p01,
                T_out=T01, length=L, Dtip=out[-1]["Dtip"], H_last=out[-1]["H3"])


AX_REF_AR = 1.3     # удлинение H/b опорной ступени: типичное для малых промышленных турбин (оценка; frame3 - eta_axial трека)


def ax_reference():
    """опорная ступень привязки: ψ 1,2, φ 0,5, U 400 м/с, Dm 0,6 м, T0 1273 K, p0 1,2 МПа, H/b = AX_REF_AR обоих рядов"""
    E = _E()
    b = E.bfrac(0.02)
    U, Dm, psi, phi = 400.0, 0.6, 1.2, 0.5
    n = 60 * U / (math.pi * Dm)
    m = 30.0
    a2 = math.atan(psi / phi); b2 = math.atan((psi - 1) / phi); b3 = math.atan(1 / phi)
    # число лопаток из H/b = AX_REF_AR -> шаг из Цвейфеля
    r = ax_turbine(1273.15, 1.2e6, b, m, U, n, psi * U * U, 1, phi, 40, 60, shift=0.0)
    H = r["stages"][0]["H2"]
    zs = []
    for a_in, a_out, Hh in ((0.0, a2, H), (b2, b3, r["stages"][0]["H3"])):
        bx = Hh / AX_REF_AR
        s = bx * ZW / (2 * math.cos(a_out) ** 2 * (math.tan(a_in) + math.tan(a_out)))
        zs.append(max(8, round(math.pi * Dm / s)))
    r = ax_turbine(1273.15, 1.2e6, b, m, U, n, psi * U * U, 1, phi, zs[0], zs[1], shift=0.0)
    return dict(eta_p_model=r["eta_p_model"], zn=zs[0], zt=zs[1], target=FS._E().eta_axial(1.2), AR_N=r["stages"][0]["AR_N"],
                AR_R=r["stages"][0]["AR_R"], M2=r["stages"][0]["M2"], model=AX_MODEL)


def ax_shift():
    if "ax" not in _SH:
        ref = ax_reference()
        _SH["ax_ref"] = ref
        _SH["ax"] = ref["target"] - ref["eta_p_model"]
    return _SH["ax"]


# =========================================================================== радиальная турбина (вид Baines; оценка)
NU5 = 0.35          # втулка / периферия на выходе ротора (оценка)
STANITZ_K = 0.63    # скольжение Stanitz на входе ротора: cθ4,opt/U4 = 1 - 0,63π/Z (Dixon & Hall, по памяти)
PHI4 = 0.25         # cm4/U4 на входе ротора (оценка; Rohlik: 0,2-0,3)
LZ_D4 = 0.35        # осевая длина ротора / D4 (оценка)


def r5_ratio(ns):
    return min(0.78, max(0.55, 0.45 + 0.45 * ns))


def rad_stage(T0in, p0in, b, m, U4, ER, n, zt, zn, k_cl=None, shift=None, eta0=0.85):
    """ступень радиальной турбины: заданы вход, расход, окружная скорость U4, степень расширения ER (полное давление входа /
    давление выхода), частота. Возврат: η_is = работа на валу / Δh_s (с привязкой RAD_SHIFT), потери, геометрия, числа Маха"""
    E = _E()
    Rg = E.Rg(b)
    shift = rad_shift() if shift is None else shift
    om = n * math.pi / 30
    r4 = U4 / om
    h0 = E.h(T0in, b)
    T5s = E.T_s(E.s0(T0in, b) - Rg * math.log(ER), b)
    dhs = h0 - E.h(T5s, b)
    C0 = math.sqrt(2 * max(dhs, 1.0))
    p5 = p0in / ER
    eta = eta0
    rho5 = p5 / (Rg * T5s)
    Q5 = m / rho5
    ns = om * math.sqrt(Q5) / max(dhs, 1.0) ** 0.75
    rr = r5_ratio(ns)
    r5t, r5h = rr * r4, NU5 * rr * r4
    r5 = math.sqrt((r5t ** 2 + r5h ** 2) / 2)
    A5 = math.pi * (r5t ** 2 - r5h ** 2) * 0.95
    Lz = LZ_D4 * 2 * r4
    kc = k_cl if k_cl is not None else max(0.4e-3, 0.001 * 2 * r4)
    L = {}
    for _ in range(30):
        w = eta * dhs
        cth4 = (w + L.get("w", 0.0)) / U4               # работа Эйлера = на валу + трение диска
        cm4 = PHI4 * U4
        c4 = math.hypot(cm4, cth4)
        a4 = math.atan2(cth4, cm4)
        T4 = E.T_h(h0 - 0.5 * c4 * c4, b)
        # сопловой аппарат: Soderberg, поворот от входа (закрутка улитки 0) до a4; шаг на радиусе 1,15 r4, b4 - высота
        p4i = p0in * math.exp((E.s0(T4, b) - E.s0(T0in, b)) / Rg)
        rho4 = p4i / (Rg * T4)
        b4 = m / (rho4 * cm4 * 2 * math.pi * r4 * 0.95)
        rowN = _row(2 * math.pi * 1.15 * r4 / zn, 0.0, a4, b4, rho4, c4, T4, b, False)
        L["N"] = rowN["zeta"] * 0.5 * c4 * c4
        # угол атаки: оптимальная закрутка по скольжению Stanitz cθ4,opt = U4·(1 - 0,63π/Z) ≈ U4·(1 - 2/Z) (Dixon & Hall, по памяти)
        wth4 = cth4 - U4
        W4 = math.hypot(cm4, wth4)
        b4f = math.atan2(wth4, cm4)
        b4o = math.atan2(U4 * (1 - STANITZ_K * math.pi / zt) - U4, cm4)
        L["inc"] = 0.5 * (W4 * math.sin(b4f - b4o)) ** 2
        # выход ротора: без закрутки
        T5 = E.T_h(h0 - w - L.get("w", 0.0), b)
        cm5 = m / (rho5 * A5)
        c5 = cm5
        T5st = E.T_h(E.h(T5, b) - 0.5 * c5 * c5, b)
        rho5 = p5 / (Rg * T5st)
        cm5 = m / (rho5 * A5)
        W5 = math.hypot(cm5, om * r5)
        W5t = math.hypot(cm5, om * r5t)
        b5b = math.atan(om * r5 / cm5)
        b5 = r5t - r5h
        chord = Lz / max(0.2, math.cos(b5b))
        Lh = math.pi / 4 * ((Lz - b4 / 2) + (r4 - r5t - b5 / 2))
        Dh = 0.5 * (4 * math.pi * r4 * b4 / (2 * math.pi * r4 + zt * b4) + 2 * math.pi * (r5t ** 2 - r5h ** 2) / (math.pi * (r5t - r5h) + zt * b5))
        L["pass"] = 0.11 * (Lh / Dh + 0.68 * (1 - (r5 / r4) ** 2) * math.cos(b5b) / (b5 / chord)) * 0.5 * (W4 ** 2 + W5 ** 2)
        Ca = (1 - r5t / r4) / (cm4 * b4)
        Cr = (r5t / r4) * max(Lz - b4, 0.0) / (cm5 * r5 * b5)
        L["cl"] = U4 ** 3 * zt / (8 * math.pi) * (0.4 * kc * Ca + 0.75 * kc * Cr - 0.3 * kc * math.sqrt(Ca * Cr))
        L["cl"] = max(0.0, L["cl"])
        Re = rho4 * U4 * r4 / mu_gas(T4)
        fdf = 0.0622 / Re ** 0.2 if Re >= 3e5 else 2.67 / Re ** 0.5
        L["w"] = fdf * 0.5 * (rho4 + rho5) * r4 ** 2 * U4 ** 3 / (4 * m)
        L["exit"] = 0.5 * c5 * c5
        eta_m = (dhs - L["N"] - L["inc"] - L["pass"] - L["cl"] - L["exit"] - L["w"]) / dhs
        eta_new = min(0.95, max(0.2, eta_m + shift))
        if abs(eta_new - eta) < 1e-6:
            eta = eta_new
            break
        eta = 0.5 * (eta + eta_new)
    Rgg = Rg
    aa4 = math.sqrt(gam(T4, b) * Rgg * T4); aa5 = math.sqrt(gam(T5st, b) * Rgg * T5st)
    cpg = cp(T4, b)
    T05 = E.T_h(h0 - eta * dhs, b)
    p05 = p5 * math.exp((E.s0(T05, b) - E.s0(E.T_h(E.h(T05, b) - 0.5 * c5 * c5, b), b)) / Rg)     # полное за ротором
    return dict(eta_is=eta, eta_is_model=eta_m, shift=shift, dhs=dhs, nu=U4 / C0, ns=ns, U4=U4, D4=2 * r4, b4=b4, D5t=2 * r5t, p5=p5,
                p05=p05, T05=T05,
                D5h=2 * r5h, Lz=Lz, alpha4_deg=math.degrees(a4), beta4_deg=math.degrees(b4f), beta4opt_deg=math.degrees(b4o),
                M4=c4 / aa4, Mw5t=W5t / aa5, cm5_U4=cm5 / U4, k_cl=kc, losses_pct={k: 100 * v / dhs for k, v in L.items()},
                T_rel_in=T4 + 0.5 * W4 * W4 / cpg, T_out=T5, ER=ER, AR_N=rowN["AR"])


def rad_turbine(T0in, p0in, b, m, U4, ER_tot, n, nst, zt, zn, shift=None):
    """nst ступеней с равной степенью расширения (полное на входе / статическое на выходе ступени, как уровень 0);
    полное давление на входе следующей ступени - статическое за предыдущей плюс доля (1 - LEAVE_K) энергии выходной скорости
    (межступенный канал, оценка); степень ступени подбирается так, чтобы статическое давление за последней было p0in/ER_tot.
    Политропный эквивалент для цикла (как _radial_stage)"""
    E = _E()
    Rg = E.Rg(b)
    ERs = ER_tot ** (1.0 / nst)
    p_end = p0in / ER_tot
    for _ in range(8 if nst > 1 else 1):
        T, p = T0in, p0in
        st = []
        for k in range(nst):
            s = rad_stage(T, p, b, m, U4, ERs, n, zt, zn, shift=shift)
            st.append(s)
            T = E.T_h(E.h(T, b) - s["eta_is"] * s["dhs"], b)
            p = s["p5"] * (s["p05"] / s["p5"]) ** (1.0 - LEAVE_K) if k < nst - 1 else s["p5"]
        if nst == 1 or abs(p / p_end - 1) < 1e-5:
            break
        ERs *= (p / p_end) ** (1.0 / nst)
    eta_p = (E.s0(T0in, b) - E.s0(T, b)) / (Rg * math.log(ER_tot))
    return dict(eta_p=eta_p, stages=st, nu=st[0]["nu"], ns=st[0]["ns"], dhs=st[0]["dhs"], eta_is=st[0]["eta_is"],
                eta_is_model=st[0]["eta_is_model"], D4=st[0]["D4"], length=sum(s["Lz"] for s in st) * 1.3, T_out=T)


def rad_reference():
    """опорная ступень привязки: U/C0 = 0,70, ns = 0,60, T0 1223 K, p0 0,45 МПа, ER 3, расход 5 кг/с, Z 13, сопел 18"""
    E = _E()
    b = E.bfrac(0.02)
    T0, p0, ER, m = 1223.15, 4.5e5, 3.0, 5.0
    Rg = E.Rg(b)
    T5s = E.T_s(E.s0(T0, b) - Rg * math.log(ER), b)
    dhs = E.h(T0, b) - E.h(T5s, b)
    U4 = 0.70 * math.sqrt(2 * dhs)
    T5 = E.T_h(E.h(T0, b) - 0.86 * dhs, b)
    Q5 = m / (p0 / ER / (Rg * T5))
    om = 0.60 * dhs ** 0.75 / math.sqrt(Q5)
    n = om * 30 / math.pi
    s = rad_stage(T0, p0, b, m, U4, ER, n, 13, 18, shift=0.0)
    return dict(eta_is_model=s["eta_is_model"], nu=s["nu"], ns=s["ns"], target=0.86, n=n)


def rad_shift():
    if "rad" not in _SH:
        ref = rad_reference()
        _SH["rad_ref"] = ref
        _SH["rad"] = ref["target"] - ref["eta_is_model"]
    return _SH["rad"]


def calibration_summary():
    """сдвиги привязки и сверки (для results-JSON)"""
    cc_shift(); ax_shift(); rad_shift()
    cal = _SH["cc_cal"]
    return dict(
        cc=dict(a=round(_SH["cc"][0], 5), b_per_PR=round(_SH["cc"][1], 5), PR_ref=PR_REF, tol=TOL_CC,
                max_abs_after=round(cal["max_abs_after"], 5),
                rows=[{k: (round(v, 5) if isinstance(v, float) else v) for k, v in r.items()} for r in cal["rows"]],
                basis="сдвиг a + b·(π - 4,5) к gtu_explore.comp_stage, МНК по сетке CAL_GRID (μ 0,70, φ 0,085, cm2/U2 0,28); вне сетки - "
                      "не выше eta_cap_l0 (формула уровня 0 той же ступени + TOL_CC)"),
        ax=dict(shift=round(_SH["ax"], 5), ref={k: (round(v, 5) if isinstance(v, float) else v) for k, v in _SH["ax_ref"].items()},
                basis="eta_axial(1,2) = 0,885 при φ 0,5, H/b 1,3 (Kacker - Okapuu по TurboFlow 0.1.18 + зазор Kacker - Okapuu)" if AX_MODEL == "ko"
                else "eta_axial(1,2) = 0,885 при φ 0,5, H/b 1,3 (Soderberg + зазор Kacker - Okapuu)"),
        rad=dict(shift=round(_SH["rad"], 5), ref={k: (round(v, 5) if isinstance(v, float) else v) for k, v in _SH["rad_ref"].items()},
                 basis="η_is 0,86 при U/C0 0,70 и ns 0,60 (gtu_explore.turb_radial)"))
