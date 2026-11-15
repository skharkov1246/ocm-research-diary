"""Воронка двигателя-заменителя, уровень 1: согласование узлов на нерасчётных режимах, запас по помпажу, пуск. Всё - оценка.

  * Компрессор: характеристика каждой центробежной (диагональной) ступени - одномерная модель той же геометрии на нерасчётном
    режиме (meanline_l1.cc_off: набор потерь Oh, Yoon, Chung 1997, предложенный авторами для расчёта характеристик; как
    ../frame3/frame3_point.stage_off); ступени - последовательно на одном валу. Осевая ступень (ОЦ) - обобщённая: ψ = 1 - (1 - ψd)·φ/φd
    (выходные углы постоянны), η = ηd·(1 - 2(φ/φd - 1)²), срыв при φ/φd < 0,75 (оценка).
  * Граница помпажа на линии постоянной частоты - наибольший расход, при котором наступает первое из: вершина π (как
    frame3_point.surge_line) или признак срыва по одномерной модели - угол атаки лопаток диффузора > DINC_STALL, угол атаки на
    периферии входа колеса > TIP_INC_STALL, угол потока за колесом > ALPHA2_STALL (пороги - по памяти, оценка; исправление
    скептика уровня 1, 15.11.2026: у набора потерь Oh нет модели срыва, вершина π лежала при 62-65 % расхода, а признаки срыва -
    при 72-85 %); запас ΔSM = (π_гр/π)·(G/G_гр) - 1 при той же частоте.
  * ВНА первой ступени (поворотный входной направляющий аппарат): закрутка α1 по ходу вращения, потери ряда - приращение
    Soderberg 0,06·(α1/100)² от c1²/2 (оценка); закон - открыт, пока режим решается с запасом SM_MIN без запирания, иначе
    прикрывается ступенями A1_STEPS. Если ВНА понадобился хоть на одном режиме - узел ВНА входит в детали и модули
    (funnel_l1.add_vigv).
  * Запирание компрессора: при расходе больше запирающего ступень работает на вертикальной ветви характеристики - расход
    постоянен, работа та же, π меньше (скачки), КПД падает (оценка; так ведёт себя характеристика за запиранием). Режим на
    этой ветви помечается (choke), глубина - отношение π к π на границе запирания.
  * Подогрев входа (холодный день): часть воздуха за компрессором xb перепускается на вход, смесь на входе компрессора -
    T1H_STEPS; xb - из баланса энтальпий h(T1h) - h(Ta) = xb·(h(T2) - h(Ta)); турбина и камера получают (1 - xb) расхода, работа
    сжатия перепуска теряется (оценка). Включается, только если без него при этой температуре не решён ни один режим.
  * Режим не решён, если нет совместной работы ни при одной T41 сетки, ни при одном положении ВНА и подогреве; такие режимы
    считаются поштучно (n_unsolved) и входят в правило полноты funnel_l1 (od_cover).
  * Турбины: пропускная способность по закону конуса Стодолы G·√T/p = K·√(1 - ER⁻²) (как frame3_point.gg_point); КПД по
    отношению скоростей η/ηd = 1 - K_NU·(ν/νd - 1)² (K_NU = 1 - как ../shaft/shaft_taurus60_gs.toml, допущение трека).
  * Камера, охлаждение, утечки - как цикл уровня 0 (funnel_space._core); рекуператор - UA постоянно (ламинарные каналы),
    ε по NTU = UA/(G·cp).
  * Режимы: температура воздуха T_AMB (−30, +15, +45 °C), нагрузка LOADS (100, 75, 50 % от доступной мощности при этой
    температуре; доступная - при расчётной T41, но не больше паспортной). Одновальный - частота вала постоянна (генератор);
    двухвальный - частота газогенератора свободна, силовая турбина: генератор - постоянная частота, «по нагрузке» -
    n ∝ нагрузка^(1/3) (закон кубической параболы привода компрессора, оценка; считается при +15 °C).
  * Пуск: частоты 0,3-0,9 расчётной при T41 пуска; самоподдержание - наименьшая частота, где мощность на валу > 0 при
    устойчивом компрессоре; мощность стартера (оценка) = 0,5·мощность компрессора·(n_сп)³ + разгон ротора двигателя за 60 с.
"""
import math, os, sys

HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in sys.path:
    sys.path.insert(0, HERE)
import funnel_space as FS        # noqa: E402
import meanline_l1 as ML         # noqa: E402

T_AMB = (-30.0, 15.0, 45.0)
LOADS = (1.00, 0.75, 0.50)
K_NU = 1.0
SM_MIN = 0.10                     # требуемый запас по помпажу на установившихся режимах (оценка: принято 10 %)
T41_START_C = 900.0               # наибольшая T41 при пуске (оценка; не выше расчётной)
N_START = (0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)
T_ACCEL_S = 60.0                  # время разгона ротора до самоподдержания (оценка)
T41_STEP = 40.0
T41_STEPS = 12
RES_TOL = 2e-3                    # невязка пропускной способности в найденном корне (иначе корень на разрыве - отказ)
DINC_STALL = 5.0                  # угол атаки лопаток диффузора, град: признак срыва (по памяти, оценка)
TIP_INC_STALL = 8.0               # угол атаки на периферии входа колеса, град: признак срыва (по памяти, оценка)
ALPHA2_STALL = 78.0               # угол потока за колесом от радиального, град (как правило расчётной точки ALPHA2_MAX)
A1_STEPS = (0.0, 15.0, 30.0, 45.0)   # положения ВНА, град закрутки (оценка: ход поворотных ВНА - до 50-60°)
CHOKE_X_MIN = 0.55                # наименьшее отношение π на ветви запирания (глубже - режима нет; оценка)
T1H_STEPS = (-10.0, 5.0)          # подогрев входа перепуском воздуха за компрессором, °C на входе компрессора (оценка) - в холодный
#                                   день, если без него ни один режим не решён: ниже приведённая частота; цена - перепуск xb


def _brent(f, a, b, fa=None, fb=None, xtol=1e-4, maxit=40):
    """корень на [a, b] при разных знаках (метод Брента, упрощённый: секущая + деление пополам)"""
    fa = f(a) if fa is None else fa
    fb = f(b) if fb is None else fb
    if fa == 0:
        return a
    if fb == 0:
        return b
    if (fa > 0) == (fb > 0):
        return None
    for _ in range(maxit):
        m = 0.5 * (a + b)
        s = b - fb * (b - a) / (fb - fa) if fb != fa else m
        if not (min(a, b) < s < max(a, b)) or abs(s - m) > 0.75 * abs(b - a):
            s = m
        fs = f(s)
        if fs == 0 or abs(b - a) < xtol * max(1.0, abs(s)):
            return s
        if (fs > 0) == (fa > 0):
            a, fa = s, fs
        else:
            b, fb = s, fs
        if abs(b - a) < xtol * max(1.0, abs(s)):
            return 0.5 * (a + b)
    return 0.5 * (a + b)


def _scan_root(f, lo, hi, guess, k=9):
    """корень f на [lo, hi]: сетка из k точек (None - нет решения в точке, например запирание); пара соседних точек с разными
    знаками, ближайшая к guess; если корень у границы области решений (точка с решением рядом с точкой без него) - граница
    уточняется делением пополам и проверяется знак на ней; затем _brent. None - корня нет"""
    xs = [lo + (hi - lo) * i / (k - 1) for i in range(k)]
    fs = [f(x) for x in xs]
    pairs = []
    for i in range(k - 1):
        a, b, fa, fb = xs[i], xs[i + 1], fs[i], fs[i + 1]
        if fa is not None and fb is not None:
            if (fa > 0) != (fb > 0):
                pairs.append((a, b, fa, fb))
            continue
        if (fa is None) == (fb is None):
            continue
        v, iv, fv = (a, b, fa) if fa is not None else (b, a, fb)        # с решением / без решения
        for _ in range(7):
            m = 0.5 * (v + iv)
            fm = f(m)
            if fm is None:
                iv = m
            else:
                if (fm > 0) != (fv > 0):
                    pairs.append((v, m, fv, fm) if v < m else (m, v, fm, fv))
                    break
                v, fv = m, fm
    if not pairs:
        return None
    a, b, fa, fb = min(pairs, key=lambda p: abs(0.5 * (p[0] + p[1]) - guess))
    return _brent(lambda x: (lambda v: (1.0 if x > 0.5 * (a + b) else -1.0) * (1 if fb > 0 else -1) if v is None else v)(f(x)),
                  a, b, fa, fb, xtol=2e-5)


# =========================================================================== компрессор
def _stall(r):
    """признаки срыва ступени по одномерной модели (центробежная / диагональная): список кодов"""
    out = []
    if r.get("dinc_deg") is not None and r["dinc_deg"] > DINC_STALL:
        out.append("диффузор")
    if r.get("inc_tip_deg") is not None and r["inc_tip_deg"] > TIP_INC_STALL:
        out.append("периферия входа")
    if r.get("alpha2_deg") is not None and r["alpha2_deg"] > ALPHA2_STALL:
        out.append("α2")
    return out


def comp_chain(D, W1, n, Ta_K, a1=0.0, choke=None, upto=None):
    """ступени последовательно; a1 - закрутка ВНА первой ступени, град; choke = (k, x) - ступень k на ветви запирания с π·x;
    upto - только ступени 0..upto. Возврат: dict(out - для funnel_space._core или None, ok, why, k_fail - номер ступени отказа,
    stall - признаки срыва)"""
    E = FS._E()
    R = E.RA
    T, p = Ta_K, E.P0
    Wc = W1 / D["cores"]
    ch = dict([choke]) if isinstance(choke, tuple) else dict(choke or {})
    out, ok, why, stall = [], True, None, []
    for k, st in enumerate(D["stages"]):
        if upto is not None and k > upto:
            break
        a = math.radians(a1) if k == 0 else 0.0
        if st["type"] == "О":
            phr = (Wc * T / p) / (st["W_d"] * st["T1_d"] / st["p1_d"]) * (st["n_d"] / n)
            phi_d = 0.5
            psi = 1 - (1 - st["psi_d"]) * phr - phr * phi_d * math.tan(a)          # закрутка ВНА снижает работу (выходные углы постоянны)
            U = st["U_d"] * n / st["n_d"]
            w = max(psi, 0.02) * U * U
            eta = max(0.4, st["eta_d"] * (1 - 2 * (phr - 1) ** 2))
            T2 = E.T_h(E.h(T) + w)
            PR = math.exp(eta * (E.s0(T2) - E.s0(T)) / R)
            if phr < 0.75:
                ok, why = False, "осевая ступень: срыв (φ/φd < 0,75)"
                stall.append(f"ст. {k + 1}: осевая")
        else:
            ms = Wc / 2 if st["dbl"] else Wc
            pk = p
            if a > 0:
                s_in = ML.inlet(T, p, ms, st["ml"]["g"]["A1"], a)
                if s_in is None:
                    return dict(out=None, ok=False, why="вход за ВНА заперт", k_fail=k, stall=stall)
                pk = p * math.exp(-0.06 * (a1 / 100.0) ** 2 * 0.5 * s_in[0] ** 2 / (R * T))      # потери ВНА сверх открытого
            r = ML.cc_off(st["ml"], T, pk, ms, n, alpha1=a)
            if not r.get("PRs"):
                return dict(out=None, ok=False, why=r.get("why"), k_fail=k, stall=stall)
            if not r["ok"] and k not in ch:
                return dict(out=None, ok=False, why=r["why"], k_fail=k, stall=stall)
            PR, eta, T2, w = r["PRs"] * (pk / p), r["eta_p"], r["T02"], r["w0"]
            if k in ch:
                PR *= ch[k]
            if pk != p or k in ch:
                eta = R * math.log(max(PR, 1.0001)) / (E.s0(T2) - E.s0(T))
            stall += [f"ст. {k + 1}: {x}" for x in _stall(r)]
        out.append(dict(T1=T, p1=p, T2=T2, p2=p * PR, w=w, pi=PR, eta=eta))
        T, p = T2, p * PR
    return dict(out=out, ok=ok, why=why, k_fail=None, stall=stall)


def comp_map(D, W1, n, Ta_K, a1=0.0, choke=None):
    """(comp_out для funnel_space._core, ok, причина) - обёртка comp_chain"""
    r = comp_chain(D, W1, n, Ta_K, a1, choke)
    return r["out"], r["ok"], r["why"]


def surge(D, n, Ta_K, W_op, a1=0.0):
    """граница помпажа на линии частоты n: наибольший расход, при котором наступает первое из - вершина π (12 точек + золотое
    сечение) или признак срыва (деление пополам между соседними точками сетки)"""
    Wd = D["G"] * math.sqrt(288.15 / Ta_K) * (n / D["n_d"])
    ws = [Wd * (0.30 + 0.08 * k) for k in range(13)]
    pts = []
    for w in ws:
        r = comp_chain(D, w, n, Ta_K, a1)
        if r["out"]:
            pts.append((w, math.prod(s["pi"] for s in r["out"]), r["stall"]))
    if not pts:
        return None
    k = max(range(len(pts)), key=lambda i: pts[i][1])
    w_s, pr_s = pts[k][0], pts[k][1]
    left_edge = k == 0
    if 0 < k < len(pts) - 1:
        a, b = pts[k - 1][0], pts[k + 1][0]
        gr = (math.sqrt(5) - 1) / 2
        f = lambda w: (lambda co: -math.prod(s["pi"] for s in co[0]) if co[0] else 0.0)(comp_map(D, w, n, Ta_K, a1))
        x1, x2 = b - gr * (b - a), a + gr * (b - a)
        f1, f2 = f(x1), f(x2)
        for _ in range(8):
            if f1 < f2:
                b, x2, f2 = x2, x1, f1
                x1 = b - gr * (b - a); f1 = f(x1)
            else:
                a, x1, f1 = x1, x2, f2
                x2 = a + gr * (b - a); f2 = f(x2)
        w_s = 0.5 * (a + b)
        co, _, _ = comp_map(D, w_s, n, Ta_K, a1)
        if co:
            pr_s = math.prod(s["pi"] for s in co)
    crit, what = "вершина π", None
    # признак срыва: первый сверху по расходу
    ks = [i for i in range(len(pts)) if pts[i][2]]
    if ks:
        i = max(ks)
        if i == len(pts) - 1 or pts[i][0] > w_s:
            lo, hi = pts[i][0], (pts[i + 1][0] if i + 1 < len(pts) else pts[i][0])
            what = pts[i][2]
            for _ in range(8 if hi > lo else 0):
                m = 0.5 * (lo + hi)
                r = comp_chain(D, m, n, Ta_K, a1)
                if r["out"] and r["stall"]:
                    lo, what = m, r["stall"]
                else:
                    hi = m
            if lo > w_s:
                r = comp_chain(D, lo, n, Ta_K, a1)
                if r["out"]:
                    w_s, pr_s, crit = lo, math.prod(s["pi"] for s in r["out"]), "признак срыва: " + ", ".join(what)
                    left_edge = False
    return dict(W=w_s, PR=pr_s, left_edge=left_edge, criterion=crit)


# =========================================================================== турбины и ядро
def _eta_nu(eta_d, nu, nu_d):
    return max(0.3, eta_d * (1 - K_NU * (nu / nu_d - 1) ** 2))


def _nu(U, T_in, b, ER):
    E = FS._E()
    dhs = E.h(T_in, b) - E.h(E.T_s(E.s0(T_in, b) - E.Rg(b) * math.log(max(ER, 1.0001)), b), b)
    return U / math.sqrt(2 * max(dhs, 1.0))


def _eps(D, W1):
    """ε противотока при UA = const (ламинарные каналы): NTU = NTU_d·G_d/G, отношение теплоёмкостей - расчётное"""
    if not D["rec"]:
        return 0.0
    ntu = D["NTU_d"] * D["G"] / W1
    cr = D["Cr_d"]
    if cr > 0.999:
        return ntu / (1 + ntu)
    e = math.exp(-ntu * (1 - cr))
    return (1 - e) / (1 - cr * e)


def core(D, co, W1, T41C, n_gg, n_pt=None):
    """ядро при заданных компрессоре, расходе и T41: КПД турбин по отношению скоростей (2 прохода)"""
    c = dict(D["c"], t41=T41C, eps_rec=_eps(D, W1))
    eg, ep = D["eta_gg_d"], D["eta_pt_d"]
    r = None
    for _ in range(3):
        r = FS._core(c, co, eg, ep if ep else 0.85)
        if "fail" in r:
            return r
        t = r["t"]
        nu = _nu(D["U_gg_d"] * n_gg / D["n_d"], t["T_in"], t["b"], t["ER"] ** (1.0 / D["nt"]))
        eg_n = _eta_nu(D["eta_gg_d"], nu, D["nu_gg_d"])
        ep_n = ep
        if r["pt"]:
            npt = n_pt or D["n_pkg"]
            p = r["pt"]
            nup = _nu(D["U_pt_d"] * npt / D["n_pkg"], p["T_in"], p["b"], p["ER"] ** (1.0 / D["npt"]))
            ep_n = _eta_nu(D["eta_pt_d"], nup, D["nu_pt_d"])
        if abs(eg_n - eg) < 1e-5 and (not r["pt"] or abs(ep_n - ep) < 1e-5):
            break
        eg, ep = eg_n, ep_n
    r["eg"], r["ep"] = eg, ep
    return r


def capacity(W, T, p, ER):
    return W * math.sqrt(T) / (p * math.sqrt(max(1e-9, 1 - ER ** -2)))


# =========================================================================== режимы
def _tc(Ta, t1h):
    """температура на входе компрессора, K: окружающая или после подогрева перепуском"""
    return Ta if t1h is None else max(Ta, t1h + 273.15)


def _bleed(co, Ta, Tc):
    """доля перепуска на подогрев входа и работа сжатия на единицу расхода через камеру: (xb, ступени с работой / (1 - xb))"""
    if Tc <= Ta + 1e-9:
        return 0.0, co
    E = FS._E()
    xb = (E.h(Tc) - E.h(Ta)) / (E.h(co[-1]["T2"]) - E.h(Ta))
    return xb, [dict(x, w=x["w"] / (1.0 - xb)) for x in co]


def point_single(D, Ta_C, T41C, n=None, W_br=None, a1=0.0, t1h=None):
    """одновальный: частота n (по умолчанию расчётная), T41 -> расход из совместной работы компрессора и турбины"""
    Ta = Ta_C + 273.15
    Tc = _tc(Ta, t1h)
    n = n or D["n_d"]
    cache = {}

    def resid(W1):
        co, ok, why = comp_map(D, W1, n, Tc, a1)
        if co is None:
            return None
        xb, cs = _bleed(co, Ta, Tc)
        r = core(D, cs, W1 * (1 - xb), T41C, n)
        if "fail" in r:
            return None
        t = r["t"]
        cache[W1] = (co, ok, why, r, xb)
        return capacity(W1 * (1 - xb) / D["cores"] * t["m"], t["T_in"], t["p_in"], t["ER"]) / D["K_gg"] - 1.0
    Wd = D["G"] * math.sqrt(288.15 / Tc) * (n / D["n_d"])
    a, b = W_br or (0.35 * Wd, 1.45 * Wd)
    W1 = _scan_root(resid, a, b, 0.5 * (a + b), k=7)
    if W1 is None:
        return None
    v = resid(W1)
    if v is None or abs(v) > RES_TOL or W1 not in cache:      # корень на разрыве (запирание входа), а не настоящий
        return None
    co, ok, why, r, xb = cache[W1]
    op = _pack(D, W1, co, ok, why, r, T41C, n, None, Ta, xb, Tc)
    op.update(a1=a1, choke=False, choke_x=None, choke_stage=None)
    return op


def choke_flow(D, n, Ta_K, a1=0.0, k=None):
    """запирающий расход линии частоты: наибольший W1, при котором ступени 0..k (k None - все) решаются; и номер ступени,
    которая запирается первой (None - не найден)"""
    Wd = D["G"] * math.sqrt(288.15 / Ta_K) * (n / D["n_d"])
    ws = [Wd * (0.4 + 0.1 * j) for j in range(12)]
    lo = hi = None
    for w in ws:
        r = comp_chain(D, w, n, Ta_K, a1, upto=k)
        if r["out"]:
            lo = w
        elif lo is not None:
            hi = w
            break
    if lo is None or hi is None:
        return None, None
    for _ in range(14):
        m = 0.5 * (lo + hi)
        if comp_chain(D, m, n, Ta_K, a1, upto=k)["out"]:
            lo = m
        else:
            hi = m
    return lo, comp_chain(D, hi, n, Ta_K, a1, upto=k)["k_fail"]


def point_choke(D, Ta_C, T41C, n=None, a1=0.0, n_pt=None, t1h=None):
    """режим на ветви запирания: расход - запирающий, ступень, которая запирается, отдаёт π·x при той же работе (x - из
    пропускной способности турбины); одновальный (частота задана)"""
    Ta = Ta_C + 273.15
    Tc = _tc(Ta, t1h)
    n = n or D["n_d"]
    for kk in range(len(D["stages"])):                 # запирается ступень kk; ступени до неё - на рабочей ветви
        if D["stages"][kk]["type"] == "О":
            continue
        op = _choke_k(D, Ta, Tc, T41C, n, a1, n_pt, kk)
        if op:
            return op
    return None


def _choke_k(D, Ta, Tc, T41C, n, a1, n_pt, kk):
    """ветвь запирания со ступени kk: расход - её запирающий; π ступени kk · x (та же работа); если следующая ступень запирается
    раньше, чем найден корень, - x ступени kk фиксируется на этой границе и запирается следующая (обе на вертикальной ветви)"""
    Wch, k = choke_flow(D, n, Tc, a1, kk)
    if Wch is None or k != kk:
        return None
    cache = {}
    cks = [j for j in range(kk, len(D["stages"])) if D["stages"][j]["type"] != "О"]
    fixed = {}

    def resid_for(j):
        def resid(x):
            chs = {**fixed, j: x}
            rr = comp_chain(D, Wch, n, Tc, a1, choke=chs)
            if rr["out"] is None:
                return None
            xb, cs = _bleed(rr["out"], Ta, Tc)
            r = core(D, cs, Wch * (1 - xb), T41C, n, n_pt)
            if "fail" in r:
                return None
            t = r["t"]
            cache[tuple(sorted(chs.items()))] = (rr, r, xb, chs)
            return capacity(Wch * (1 - xb) / D["cores"] * t["m"], t["T_in"], t["p_in"], t["ER"]) / D["K_gg"] - 1.0
        return resid
    for j in cks:
        f = resid_for(j)
        x = _scan_root(f, CHOKE_X_MIN, 1.0, 1.0, k=7)
        if x is not None:
            v = f(x)
            key = tuple(sorted({**fixed, j: x}.items()))
            if v is not None and abs(v) <= RES_TOL and key in cache:
                rr, r, xb, chs = cache[key]
                op = _pack(D, Wch, rr["out"], rr["ok"], rr["why"], r, T41C, n, n_pt, Ta, xb, Tc)
                op.update(a1=a1, choke=True, choke_x=min(chs.values()), choke_stage="/".join(str(q + 1) for q in sorted(chs)))
                return op
        # корня нет: наименьший x ступени j, при котором ступени дальше решаются, и остаток там ещё < 0 - запирается следующая
        lo, hi = CHOKE_X_MIN, 1.0
        if f(hi) is None:
            return None
        if f(lo) is not None:
            return None                                  # вся ветвь решается, корня нет - глубже CHOKE_X_MIN
        for _ in range(12):
            m = 0.5 * (lo + hi)
            if f(m) is None:
                lo = m
            else:
                hi = m
        if f(hi) is None or f(hi) > 0:
            return None
        fixed[j] = hi
    return None


def point_two(D, Ta_C, T41C, n_pt=None, N_br=None, a1=0.0, t1h=None):
    """двухвальный: частота газогенератора и расход - из пропускной способности турбины газогенератора и силовой турбины"""
    Ta = Ta_C + 273.15
    Tc = _tc(Ta, t1h)
    store = {}

    def gg(N):
        Wd = D["G"] * math.sqrt(288.15 / Tc) * (N / D["n_d"])
        cache = {}

        def resid(W1):
            co, ok, why = comp_map(D, W1, N, Tc, a1)
            if co is None:
                return None
            xb, cs = _bleed(co, Ta, Tc)
            r = core(D, cs, W1 * (1 - xb), T41C, N, n_pt)
            if "fail" in r:
                return None
            t = r["t"]
            cache[W1] = (co, ok, why, r, xb)
            return capacity(W1 * (1 - xb) / D["cores"] * t["m"], t["T_in"], t["p_in"], t["ER"]) / D["K_gg"] - 1.0
        W1 = _scan_root(resid, 0.5 * Wd, 1.3 * Wd, Wd, k=9)
        if W1 is None:
            return None
        v = resid(W1)
        if v is None or abs(v) > RES_TOL or W1 not in cache:
            return None
        co, ok, why, r, xb = cache[W1]
        p = r["pt"]
        val = capacity(W1 * (1 - xb) * p["m"], p["T_in"], p["p_in"], p["ER"]) / D["K_pt"] - 1.0
        store[N] = (W1, co, ok, why, r, xb)
        return val

    a, b = N_br or (0.45 * D["n_d"], 1.20 * D["n_d"])
    v0 = gg(0.5 * (a + b))
    if v0 is not None and abs(v0) <= RES_TOL:
        N = 0.5 * (a + b)
    else:
        N = _scan_root(gg, a, b, 0.5 * (a + b), k=5)
        if N is None:
            return None
        v = gg(N)
        if v is None or abs(v) > RES_TOL:
            return None
    W1, co, ok, why, r, xb = store[N]
    op = _pack(D, W1, co, ok, why, r, T41C, N, n_pt, Ta, xb, Tc)
    op.update(a1=a1, choke=False, choke_x=None, choke_stage=None)
    return op


def _pack(D, W1, co, ok, why, r, T41C, n, n_pt, Ta, xb=0.0, Tc=None):
    E = FS._E()
    Tc = Tc or Ta
    PR = math.prod(s["pi"] for s in co)
    P = W1 * (1 - xb) * r["wnet"]
    return dict(Ta_C=Ta - 273.15, T41_C=T41C, n=n, n_pt=n_pt, W1=W1, PR=PR, P_kW=P / 1e3, eta_th=r["eta"], EGT_C=r["Tstack"] - 273.15,
                eta_gg=r["eg"], eta_pt=r["ep"] if r["pt"] else None, comp_ok=ok, comp_why=why, co=co, r=r, x_bleed=xb,
                T1_C=Tc - 273.15, t1h=(Tc - 273.15) if Tc > Ta + 1e-9 else None,
                Wc_corr=W1 * math.sqrt(Tc / 288.15) / (E.P0 / 101325.0))


def solve(D, Ta_C, T41C, n_pt=None, prev=None, a1=0.0, t1h=None):
    """режим при T41, положении ВНА и подогреве входа: корень, ближайший к предыдущей точке сетки (или к расчётной при той же
    приведённой частоте) - расширяющиеся интервалы поиска (у двухвального уравнения пропускной способности бывает несколько
    корней); у одновального без корня на рабочей ветви - ветвь запирания (point_choke)"""
    Ta = Ta_C + 273.15
    Tc = _tc(Ta, t1h)
    if D["rotors"] == 1:
        g = prev["W1"] if prev else D["G"] * math.sqrt(288.15 / Tc)
        for w in (0.08, 0.2, 0.4, 0.6):
            r = point_single(D, Ta_C, T41C, W_br=(g * (1 - w), g * (1 + w)), a1=a1, t1h=t1h)
            if r:
                return r
        return point_choke(D, Ta_C, T41C, a1=a1, t1h=t1h)
    g = prev["n"] if prev else D["n_d"] * math.sqrt(Tc / 288.15)
    for w in (0.05, 0.12, 0.25, 0.45):
        r = point_two(D, Ta_C, T41C, n_pt, N_br=(g * (1 - w), min(g * (1 + w), 1.25 * D["n_d"])), a1=a1, t1h=t1h)
        if r:
            return r
    return None


def with_sm(D, op):
    s = surge(D, op["n"], op.get("T1_C", op["Ta_C"]) + 273.15, op["W1"], op.get("a1", 0.0))
    if s is None:
        op.update(SM=None, surge_W=None, surge_PR=None, surge_criterion=None)
        return op
    op.update(SM=(s["PR"] / op["PR"]) * (op["W1"] / s["W"]) - 1, surge_W=s["W"], surge_PR=s["PR"], surge_left_edge=s["left_edge"],
              surge_criterion=s["criterion"])
    return op


def _clean(op):
    return {k: (round(v, 5) if isinstance(v, float) else v) for k, v in op.items() if k not in ("co", "r")}


def _grid(D, Ta, a1, t1h=None):
    """сетка T41 от расчётной вниз при положении ВНА a1 и подогреве t1h: точки совместной работы (до 42 % мощности или отказа
    после точек)"""
    T41d = D["t41_d"]
    grid, prev, fails = [], None, []
    for k in range(T41_STEPS):
        T41 = T41d - k * T41_STEP
        op = solve(D, Ta, T41, prev=prev, a1=a1, t1h=t1h)
        if op is None:
            fails.append(dict(Ta_C=Ta, T41_C=T41, a1=a1, t1h=t1h, why="нет совместной работы узлов (у запирания компрессора или турбины)"))
            if grid:
                break
            continue                                  # холодный день: при расчётной T41 может не сходиться - ниже по T41
        grid.append(op); prev = op
        if op["P_kW"] < 0.42 * min(D["P_rated_kW"], grid[0]["P_kW"]):
            break
    return grid, fails


def _load_point(D, Ta, grid, tgt, P_av, npt, a1, t1h=None):
    T41d = D["t41_d"]
    T41 = _t41_for(grid, tgt)
    if T41 is None:
        return None, "мощность вне сетки T41"
    op = None
    for _ in range(3):                        # уточнение T41 по мощности (секущая по сетке)
        op = solve(D, Ta, T41, npt, prev=_near(grid, T41), a1=a1, t1h=t1h)
        if op is None:
            break
        if abs(op["P_kW"] - tgt) < 0.003 * P_av:
            break
        dP = _slope(grid, T41)
        T41 = min(T41d, T41 + (tgt - op["P_kW"]) / dP) if dP else T41
    if op is None:
        return None, "нет совместной работы узлов"
    if abs(op["P_kW"] - tgt) > 0.02 * P_av:
        return None, f"мощность {op['P_kW']:.0f} кВт вместо {tgt:.0f} (нет точки на сетке T41)"
    return with_sm(D, op), None


def matrix(D):
    """режимы: T_AMB × LOADS (+ закон «по нагрузке» силовой турбины при +15 °C); ВНА - открыт, пока режим решается с запасом
    SM_MIN без запирания, иначе прикрывается (A1_STEPS); подогрев входа (T1H_STEPS) - только если без него при этой температуре
    не решён ни один режим; запас по помпажу в каждой точке; каждый нерешённый режим - в fails"""
    P_rated = D["P_rated_kW"]
    T41d = D["t41_d"]
    out, fails = [], []
    for Ta in T_AMB:
        grids, gfails = {}, []

        def G(a1, t1h):
            if (a1, t1h) not in grids:
                grids[(a1, t1h)], f = _grid(D, Ta, a1, t1h)
                gfails.extend(f)
            return grids[(a1, t1h)]
        laws = [("генератор", None)]
        if D["rotors"] == 2 and Ta == 15.0:
            laws.append(("по нагрузке", None))
        pts, fl = None, None
        for t1h in (None,) + tuple(h for h in T1H_STEPS if h > Ta):
            P_av = None
            for a1 in A1_STEPS:
                g_ = G(a1, t1h)
                if g_:
                    P_av = min(P_rated, max(x["P_kW"] for x in g_))
                    break
            if P_av is None:
                continue
            pts, fl = [], []
            for ld in LOADS:
                tgt = ld * P_av
                for law, _ in laws:
                    npt = D["n_pkg"] * ld ** (1 / 3) if law == "по нагрузке" else None
                    best, why_last = None, None
                    for a1 in A1_STEPS:                 # ВНА прикрывается, пока нет запаса SM_MIN без запирания
                        g_ = G(a1, t1h)
                        if not g_:
                            continue
                        op, why = _load_point(D, Ta, g_, tgt, P_av, npt, a1, t1h)
                        if op is None:
                            why_last = why
                            continue
                        good = op.get("SM") is not None and op["SM"] >= SM_MIN
                        rank = (good and not op["choke"], good, 0.0 if good else (op.get("SM") if op.get("SM") is not None else -9.0))
                        if best is None or rank > best[0]:          # при равенстве - меньшее прикрытие (выше КПД)
                            best = (rank, op)
                        if good and not op["choke"]:
                            break
                    if best is None:
                        fl.append(dict(Ta_C=Ta, load=ld, law=law, t1h=t1h, why=why_last or "нет совместной работы узлов"))
                        continue
                    op = best[1]
                    op.update(load=ld, law=law, P_target_kW=tgt, P_avail_kW=P_av, derated=P_av < P_rated - 1e-6)
                    pts.append(_clean(op))
            if pts:
                break                                       # подогрев - только если без него не решён ни один режим
        fails.extend(gfails)
        if not pts:
            for ld in LOADS:
                for law, _ in laws:
                    fails.append(dict(Ta_C=Ta, load=ld, law=law, why="нет совместной работы ни при одной T41 сетки, ни при одном "
                                                                    "положении ВНА и подогреве входа"))
            continue
        out += pts
        fails += fl
    sms = [o["SM"] for o in out if o.get("SM") is not None]
    iso = next((o for o in out if o["Ta_C"] == 15.0 and o["load"] == 1.0 and o["law"] == "генератор"), None)
    n_req = len(T_AMB) * len(LOADS) + (len(LOADS) if D["rotors"] == 2 else 0)
    chk = [o for o in out if o.get("choke")]
    return dict(points=out, fails=fails, n_points=len(out), n_required=n_req,
                n_unsolved=sum(1 for f in fails if "load" in f), SM_min=min(sms) if sms else None,
                SM_req=SM_MIN, vigv_used=any((o.get("a1") or 0.0) > 0 for o in out),
                a1_max=max([o.get("a1") or 0.0 for o in out], default=0.0), n_choke=len(chk),
                choke_x_min=min([o["choke_x"] for o in chk], default=None),
                heat_used=any(o.get("t1h") is not None for o in out),
                x_bleed_max=max([o.get("x_bleed") or 0.0 for o in out], default=0.0),
                unsolved_note=("режимы без совместной работы в модели (ни при одной T41 сетки, ни при одном положении ВНА A1_STEPS и "
                               "подогреве входа T1H_STEPS, с ветвью запирания у одновального) - считаются поштучно и входят в правило "
                               "полноты od_cover"),
                SM_min_at=(min((o for o in out if o.get("SM") is not None), key=lambda o: o["SM"]) if sms else None),
                P_avail_hot_kW=next((o["P_avail_kW"] for o in out if o["Ta_C"] == max(T_AMB)), None),
                eta_th_50_iso=next((o["eta_th"] for o in out if o["Ta_C"] == 15.0 and o["load"] == 0.5 and o["law"] == "генератор"), None),
                iso_check=(dict(P_kW=iso["P_kW"], P_design_kW=P_rated, T41_C=iso["T41_C"], T41_design_C=T41d, a1=iso.get("a1"),
                                choke=iso.get("choke")) if iso else None))


def _t41_for(grid, P):
    pts = sorted(((g["P_kW"], g["T41_C"]) for g in grid))
    if not pts or P > pts[-1][0] + 1e-6 or P < pts[0][0] - 1e-6:
        if pts and P > pts[-1][0]:
            return pts[-1][1]
        return None
    for (p0, t0), (p1, t1) in zip(pts, pts[1:]):
        if p0 <= P <= p1:
            return t0 + (t1 - t0) * (P - p0) / (p1 - p0) if p1 > p0 else t0
    return pts[-1][1]


def _near(grid, T41):
    return min(grid, key=lambda g: abs(g["T41_C"] - T41))


def _slope(grid, T41):
    g = sorted(grid, key=lambda x: x["T41_C"])
    for a, b in zip(g, g[1:]):
        if a["T41_C"] <= T41 <= b["T41_C"] and b["T41_C"] > a["T41_C"]:
            return (b["P_kW"] - a["P_kW"]) / (b["T41_C"] - a["T41_C"])
    if len(g) > 1:
        return (g[-1]["P_kW"] - g[0]["P_kW"]) / (g[-1]["T41_C"] - g[0]["T41_C"])
    return None


def start(D):
    """пуск при +15 °C: самоподдержание и мощность стартера (оценка)"""
    T41s = min(T41_START_C, D["t41_d"])
    rows, n_ss = [], None
    for nf in N_START:
        if D["rotors"] == 1:
            op = point_single(D, 15.0, T41s, n=nf * D["n_d"])
        else:
            op = None                 # двухвальный: точка газогенератора на линии работы - по сетке T41 (ниже)
        if op is None:
            if D["rotors"] == 1:
                rows.append(dict(n_frac=nf, ok=False, why="нет совместной работы"))
            continue
        op = with_sm(D, op)
        good = op["P_kW"] > 0 and (op["SM"] is None or op["SM"] > 0)
        rows.append(dict(n_frac=nf, ok=True, P_kW=round(op["P_kW"], 1), SM=None if op["SM"] is None else round(op["SM"], 4),
                         T41_C=T41s, self_sustain=good))
        if good and n_ss is None:
            n_ss = nf
    if D["rotors"] == 2:
        prev = None
        for k in range(16):
            T41 = T41s - 40.0 * k
            op = solve(D, 15.0, T41, prev=prev)
            if op is None:
                break
            prev = op
            op = with_sm(D, op)
            nf = op["n"] / D["n_d"]
            good = op["SM"] is None or op["SM"] > 0
            rows.append(dict(n_frac=round(nf, 4), ok=True, P_kW=round(op["P_kW"], 1), SM=None if op["SM"] is None else round(op["SM"], 4),
                             T41_C=T41, self_sustain=good))
            if good:
                n_ss = nf if n_ss is None else min(n_ss, nf)
    P_comp = D["W_comp_kW"]
    J = D["J_kgm2"]
    om = (n_ss or 1.0) * D["n_d"] * math.pi / 30
    P_st = (0.5 * P_comp * (n_ss or 1.0) ** 3 + J * om * om / 2 / T_ACCEL_S / 1e3) if n_ss else None
    return dict(T41_start_C=T41s, rows=rows, n_self_sustain=n_ss, J_engine_kgm2=J, P_starter_kW=P_st,
                note="J - ротор двигателя (колёса и вал, оценка); у одновального генераторного пакета стартер разгоняет ещё "
                     "редуктор и генератор - их момента инерции нет в данных (нет данных)")
