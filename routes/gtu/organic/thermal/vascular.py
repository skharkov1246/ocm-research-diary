"""«Сосудистая» сеть охлаждения стенки: сопловая лопатка и жаровая труба KV-3 (уровень облика, оценка).

Модель — одномерная сеть каналов + двумерная теплопроводность стенки (как разрешено постановкой):
  стенка — пластина Lx × H (развёртка: x — от входной кромки к выходной / вдоль оси камеры, y — по высоте
  лопатки / по окружности), уравнение k·t·∇²T + h_g(T_g − T) − Σ h_c·πd·ds·(T − T_c)/A = 0, края теплоизолированы;
  охладитель в каждом элементе канала: ṁ·cp·dT_c = h_c·πd·ds·(T − T_c) (против потока — по направлению расхода);
  гидравлика — узловой баланс массы, Дарси — Вейсбах, коэффициент трения Churchill (1977), теплоотдача —
  Gnielinski (1976) в турбулентной области, Nu = 4,36 в ламинарной, плавный переход 2300…4000;
  стенка и охладитель решаются одной разреженной линейной системой; свойства воздуха — Cantera (по T),
  2–3 внешние итерации по свойствам и расходам.

Сети одинакового объёма каналов (правило Бежана: сравнивать при равном объёме):
  straight — N параллельных прямых каналов от подвода (x = 0) к выпуску (x = Lx);
  murray3  — двоичное дерево «дельта» из одной точки подвода, 2^K окончаний на выпуске, d_child = d·2^(−1/3)
             (закон Мюррея: Q ∝ d³ при ламинарном течении — Murray 1926, PNAS 12:207);
  uylings  — то же дерево, d_child = d·2^(−3/7) (показатель 7/3 для турбулентного течения — Uylings 1977,
             Bull. Math. Biol. 39:509; ссылки по выдержке поиска, см. organic/cycles/README.md);
  cco      — растущее дерево (constrained constructive optimization, Schreiner & Buxbaum 1993 — по памяти,
             не открыт): окончания (выпуски в газ, как эффузионные отверстия) разбрасываются с плотностью,
             пропорциональной тепловому потоку, каждое новое подключается к сегменту, дающему наименьший
             «объём Мюррея» Σ L·Q^(2/3); радиусы — по Мюррею от числа окончаний ниже по течению.
Производство энтропии (Бежан): S_ΔT = Σ q·(1/T_c − 1/T_w)·dA (стенка → охладитель) и S_Δp = Σ ṁ·Δp/(ρ·T).

Все числа на выходе — расчёт по этой модели; входные величины — из frame3 (точка KV-3) и оценки (assumptions).
"""
import json
import math
import sys
import time
from pathlib import Path

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
from scipy.optimize import brentq

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
GTU = ROOT / "routes" / "gtu"
SIGMA = 5.670374e-8

# ------------------------------------------------------------------ свойства воздуха (Cantera)
_AIR = {}


def air_props(T, p):
    """cp, μ, k, Pr, ρ сухого воздуха (Cantera gri30: N2/O2/Ar), табличная интерполяция по T при давлении p"""
    key = round(p, -3)
    if key not in _AIR:
        import cantera as ct
        g = ct.Solution("gri30.yaml")
        g.transport_model = "mixture-averaged"
        Ts = np.linspace(250, 2400, 88)
        tab = []
        for t in Ts:
            g.TPX = t, key, "N2:0.7808, O2:0.2095, AR:0.0093"
            tab.append((g.cp_mass, g.viscosity, g.thermal_conductivity, g.density))
        _AIR[key] = (Ts, np.array(tab))
    Ts, tab = _AIR[key]
    cp, mu, k, rho = (np.interp(T, Ts, tab[:, i]) for i in range(4))
    rho = rho * p / key
    return dict(cp=cp, mu=mu, k=k, rho=rho, Pr=cp * mu / k)


# ------------------------------------------------------------------ корреляции труб
def f_churchill(Re, rr=0.0):
    """коэффициент трения Дарси, все режимы (Churchill 1977, Chem. Eng. 84(24):91 — справочная формула, не открыт)"""
    Re = max(Re, 1e-6)
    A = (2.457 * math.log(1.0 / ((7.0 / Re) ** 0.9 + 0.27 * rr))) ** 16
    B = (37530.0 / Re) ** 16
    return 8.0 * ((8.0 / Re) ** 12 + (A + B) ** -1.5) ** (1.0 / 12.0)


def f_churchill_v(Re, rr=0.0):
    """векторный вариант f_churchill"""
    Re = np.maximum(np.asarray(Re, float), 1e-6)
    A = (2.457 * np.log(1.0 / ((7.0 / Re) ** 0.9 + 0.27 * rr))) ** 16
    B = (37530.0 / Re) ** 16
    return 8.0 * ((8.0 / Re) ** 12 + (A + B) ** -1.5) ** (1.0 / 12.0)


def seg_into(net):
    """узел → входящий сегмент (в дереве у каждого узла, кроме подвода, ровно один)"""
    if getattr(net, "_into", None) is None or len(net._into) != len(net.seg):
        d = {}
        for i, (a, b, _) in enumerate(net.seg):
            d[b] = i
        net._into = d
        net._into_len = len(net.seg)
    return net._into


def nu_pipe(Re, Pr, f):
    """Nu: ламинар 4,36 (q = const); турбулентный Gnielinski (1976; Int. Chem. Eng. 16:359 — справочная, не открыт);
    переход 2300…4000 — линейная интерполяция"""
    lam = 4.36
    if Re <= 2300:
        return lam
    def gn(R):
        ff = f_churchill(R)
        return (ff / 8) * (R - 1000) * Pr / (1 + 12.7 * math.sqrt(ff / 8) * (Pr ** (2 / 3) - 1))
    if Re >= 4000:
        return gn(Re)
    w = (Re - 2300) / 1700
    return (1 - w) * lam + w * gn(4000)


# ------------------------------------------------------------------ сети
class Net:
    """узлы (x, y), сегменты (a, b, d); узел 0 — подвод; выпуски — список узлов"""

    def __init__(self, name):
        self.name = name
        self.nodes = []
        self.seg = []
        self.outlets = []

    def node(self, x, y):
        self.nodes.append((x, y))
        return len(self.nodes) - 1

    def add(self, a, b, d):
        self.seg.append([a, b, d])

    def _ch(self):
        man = set(getattr(self, "manifold_segs", []))
        return [s for i, s in enumerate(self.seg) if i not in man]

    def volume(self):
        P = np.array(self.nodes)
        return float(sum(math.pi * d * d / 4 * np.hypot(*(P[b] - P[a])) for a, b, d in self._ch()))

    def scale_d(self, k):
        for s in self.seg:
            s[2] *= k

    def length(self):
        P = np.array(self.nodes)
        return float(sum(np.hypot(*(P[b] - P[a])) for a, b, d in self._ch()))


def net_straight(Lx, H, N, d):
    n = Net("straight")
    n.node(0.0, H / 2)                       # общий подвод (коллектор у входной кромки, без потерь)
    for i in range(N):
        y = H * (i + 0.5) / N
        a = n.node(0.0, y)
        b = n.node(Lx, y)
        n.add(0, a, 50 * d)                  # коллектор подвода — практически без сопротивления
        n.add(a, b, d)
        n.outlets.append(b)
    n.manifold_segs = list(range(0, 2 * N, 2))
    return n


def net_delta(Lx, H, K, d0, exponent, name, n_root=1):
    """n_root корней вдоль подвода (x = 0), у каждого — двоичное дерево «дельта» из K уровней:
    уровень k в x = Lx·k/K, 2^k узлов по y в полосе корня; d_child = d_parent·2^(−1/exponent).
    n_root = N, K = 0 — прямые каналы; n_root = 1 — одно дерево (развёртка по сложности, Бежан)."""
    n = Net(name)
    n.node(0.0, H / 2)                                 # общий подвод (коллектор у входной кромки)
    n.manifold_segs = []
    ratio = 2 ** (-1.0 / exponent)
    outs = []
    for r in range(n_root):
        span = H / n_root
        yc = span * (r + 0.5)
        root = n.node(0.0, yc)
        n.add(0, root, 50 * d0)
        n.manifold_segs.append(len(n.seg) - 1)
        level = [(root, yc, span)]
        d = d0
        if K == 0:
            e = n.node(Lx, yc)
            n.add(root, e, d)
            outs.append(e)
            continue
        for k in range(1, K + 1):
            x = Lx * k / K
            new = []
            for (nid, y0, sp_) in level:
                for sgn in (-1, 1):
                    y2 = y0 + sgn * sp_ / 4
                    c = n.node(x, y2)
                    n.add(nid, c, d)
                    new.append((c, y2, sp_ / 2))
            level = new
            d *= ratio
        outs += [nid for nid, _, _ in level]
    n.outlets = outs
    return n


def net_cco(Lx, H, qfun, Nterm, d_term, exponent=3.0, seed=0, n_root=4):
    """растущее дерево (CCO, упрощённо — Schreiner & Buxbaum 1993, по памяти, не открыт): n_root точек подвода у
    входной кромки; окончания — на сетке со случайным сдвигом, плотность ∝ тепловому потоку (отбор);
    каждое новое окончание подключается к сегменту леса (или к корню) с наименьшим приростом «объёма Мюррея»
    C = Σ L·N^(2/exponent) (N — число окончаний ниже по течению, Q ∝ N); прирост считается вдоль пути к корню,
    точка ветвления — взвешенная задача Ферма — Торричелли (итерации Вейсфельда), векторно по всем сегментам."""
    rng = np.random.default_rng(seed)
    gam = 2.0 / exponent
    # кандидаты окончаний: сетка с дрожанием, отбор ∝ q
    nc = int(Nterm * 2.5)
    ax_ = math.sqrt(Lx * H / nc)
    xs = np.arange(ax_ / 2, Lx, ax_)
    ys = np.arange(ax_ / 2, H, ax_)
    C = np.array([[x, y] for x in xs for y in ys]) + rng.uniform(-0.3, 0.3, (len(xs) * len(ys), 2)) * ax_
    C[:, 0] = np.clip(C[:, 0], 0.02 * Lx, Lx); C[:, 1] = np.clip(C[:, 1], 0.02 * H, 0.98 * H)
    q = qfun(C)
    q = np.maximum(q, 0)
    keep = rng.choice(len(C), size=min(Nterm, int((q > 0).sum())), replace=False, p=q / q.sum())
    pts = C[keep]
    roots = [np.array([0.0, H * (r + 0.5) / n_root]) for r in range(n_root)]
    nodes = [np.array([0.0, H / 2])] + roots
    A_, B_, N_, par = [], [], [], []      # сегменты: начало, конец (индексы узлов), число окончаний, родительский сегмент
    # порядок роста: от подвода
    pts = pts[np.argsort(pts[:, 0] + 0.2 * rng.uniform(0, Lx, len(pts)))]
    for p in pts:
        P = np.array(nodes)
        best = (np.inf, None)
        for r in range(1, n_root + 1):
            c = np.hypot(*(p - P[r]))
            if c < best[0]:
                best = (c, ("root", r))
        if A_:
            a = np.array(A_); b = np.array(B_); n = np.array(N_, float)
            Pa, Pb = P[a], P[b]
            # сумма по предкам: Σ L·((N+1)^γ − N^γ)
            L = np.hypot(*(Pb - Pa).T)
            inc = L * ((n + 1) ** gam - n ** gam)
            anc = np.zeros(len(a))
            for i in range(len(a)):                       # по родительской цепочке (дерево невелико)
                j = par[i]
                acc = 0.0
                while j is not None:
                    acc += inc[j]
                    j = par[j]
                anc[i] = acc
            w1, w2, w3 = (n + 1) ** gam, n ** gam, np.ones_like(n)
            X = (Pa + Pb + p) / 3.0
            for _ in range(30):
                d1 = np.maximum(np.hypot(*(X - Pa).T), 1e-9)
                d2 = np.maximum(np.hypot(*(X - Pb).T), 1e-9)
                d3 = np.maximum(np.hypot(*(X - p).T), 1e-9)
                W = w1 / d1 + w2 / d2 + w3 / d3
                X = (Pa * (w1 / d1)[:, None] + Pb * (w2 / d2)[:, None] + p * (w3 / d3)[:, None]) / W[:, None]
            cost = (w1 * np.hypot(*(X - Pa).T) + w2 * np.hypot(*(X - Pb).T) + np.hypot(*(X - p).T)
                    - L * n ** gam + anc)
            i = int(np.argmin(cost))
            if cost[i] < best[0]:
                best = (cost[i], ("seg", i, X[i]))
        kind = best[1]
        if kind[0] == "root":
            nodes.append(p)
            A_.append(kind[1]); B_.append(len(nodes) - 1); N_.append(1); par.append(None)
        else:
            i, X = kind[1], kind[2]
            a, b, n = A_[i], B_[i], N_[i]
            nodes += [X, p]
            ib, ip = len(nodes) - 2, len(nodes) - 1
            # предкам +1
            j = par[i]
            while j is not None:
                N_[j] += 1
                j = par[j]
            # сегмент i: a→ib (N+1); новый: ib→b (N), потомки b теперь висят на новом сегменте
            new_lower = len(A_)
            A_.append(ib); B_.append(b); N_.append(n); par.append(i)
            for k in range(len(A_) - 1):
                if par[k] == i and A_[k] == b:
                    par[k] = new_lower
            B_[i] = ib; N_[i] = n + 1
            A_.append(ib); B_.append(ip); N_.append(1); par.append(i)
    net = Net("cco")
    for qq in nodes:
        net.node(float(qq[0]), float(qq[1]))
    net.manifold_segs = []
    for r in range(1, n_root + 1):
        net.add(0, r, 50 * d_term * Nterm)
        net.manifold_segs.append(len(net.seg) - 1)
    for a, b, n in zip(A_, B_, N_):
        net.add(a, b, d_term * n ** (1.0 / exponent))
    has_child = set(A_)
    net.outlets = [b for b in B_ if b not in has_child]
    return net


# ------------------------------------------------------------------ решатель
class Panel:
    def __init__(self, Lx, H, nx, ny, k_w, t_w, hg_fun, Tg_fun, t_skin=0.0, rad=None):
        self.Lx, self.H, self.nx, self.ny = Lx, H, nx, ny
        self.dx, self.dy = Lx / nx, H / ny
        self.k_w, self.t_w, self.t_skin = k_w, t_w, t_skin
        xc = (np.arange(nx) + 0.5) * self.dx
        yc = (np.arange(ny) + 0.5) * self.dy
        self.X, self.Y = np.meshgrid(xc, yc, indexing="ij")
        self.hg = hg_fun(self.X, self.Y)
        self.Tg = Tg_fun(self.X, self.Y)
        self.rad = rad                     # (ε_g·ε_w эфф., T_g для излучения) или None

    def cell(self, x, y):
        i = min(self.nx - 1, max(0, int(x / self.dx)))
        j = min(self.ny - 1, max(0, int(y / self.dy)))
        return i * self.ny + j


def _m_of_dp(dp, rho, mu, D, A, L, Kloc=None):
    """расход по перепаду (Дарси — Вейсбах + Churchill; Kloc — местная потеря K·ρV²/2 вместо трения),
    векторно; возвращает m (со знаком) и dm/d(Δp)"""
    adp = np.maximum(np.abs(dp), 1e-9)
    if Kloc is not None:
        m = A * np.sqrt(2 * rho * adp / Kloc)
        return np.sign(dp) * m, m / (2 * adp)
    m = A * D ** 2 / (32 * mu * L) * rho * adp          # ламинарная оценка
    for _ in range(40):
        Re = m * D / (A * mu)
        f = f_churchill_v(Re)
        g = f * L / D * m ** 2 / (2 * rho * A ** 2) - adp
        # производная по m при «замороженном» наклоне f(Re): d/dm (f m²) ≈ f m (2 + dlnf/dlnm)
        eps = 1e-4
        f2 = f_churchill_v(Re * (1 + eps))
        slope = np.log(f2 / f) / np.log(1 + eps)
        dg = f * L / D * m / (2 * rho * A ** 2) * (2 + slope)
        mn = np.clip(m - g / dg, 0.2 * m, 5 * m)
        if np.max(np.abs(mn / m - 1)) < 1e-10:
            m = mn
            break
        m = mn
    Re = m * D / (A * mu)
    f = f_churchill_v(Re)
    f2 = f_churchill_v(Re * 1.0001)
    slope = np.log(f2 / f) / np.log(1.0001)
    dmdp = m / ((2 + slope) * adp)
    return np.sign(dp) * m, dmdp


def hydraulics(net, p_in, p_out, Tc_seg, p_mean, K_out=1.0, max_it=80):
    """узловые давления (метод Ньютона: якобиан — взвешенный лапласиан) и расходы; выпуск — потеря K_out·ρV²/2
    к давлению выпуска (p_out — число, функция x или net.p_drain); несжимаемо при местной плотности"""
    P = np.array(net.nodes)
    nn = len(P)
    segs = net.seg
    a_ = np.array([s[0] for s in segs]); b_ = np.array([s[1] for s in segs])
    D = np.array([d for a, b, d in segs])
    # длина не меньше диаметра: короткие сегменты у точек ветвления — это сам узел (оценка)
    L = np.array([max(np.hypot(*(P[b] - P[a])), d) for a, b, d in segs])
    A = math.pi * D ** 2 / 4
    pr = air_props(Tc_seg, p_mean)
    rho, mu = pr["rho"], pr["mu"]
    into = seg_into(net)
    outs = list(net.outlets)
    so = np.array([into[o] for o in outs])
    pdr = getattr(net, "p_drain", None)
    po = np.array([pdr if pdr is not None else (p_out(P[o][0]) if callable(p_out) else p_out) for o in outs])
    Ko = np.array([getattr(net, "K_out", {}).get(o, K_out) for o in outs])
    free = np.arange(1, nn)

    def resid(p):
        dp = p[a_] - p[b_]
        m, c = _m_of_dp(dp, rho, mu, D, A, L)
        mo, co = _m_of_dp(p[outs] - po, rho[so], mu[so], D[so], A[so], L[so], Kloc=Ko)
        R = np.zeros(nn)
        np.add.at(R, a_, -m)
        np.add.at(R, b_, m)
        np.add.at(R, outs, -mo)
        return R, m, c, mo, co

    def jac(c, co):
        rows = np.concatenate([a_, b_, a_, b_, outs])
        cols = np.concatenate([a_, b_, b_, a_, outs])
        vals = np.concatenate([c, c, -c, -c, co])
        J = sp.coo_matrix((vals, (rows, cols)), shape=(nn, nn)).tocsr()
        return J[free][:, free].tocsc()
    # начальное приближение: линейная задача (ламинарные проводимости, выпуск — по скоростному напору при 10 м/с)
    c0 = A * D ** 2 / (32 * mu * L) * rho
    co0 = A[so] * rho[so] * 10.0 / (Ko * 0.5 * rho[so] * 100.0)
    J0 = jac(c0, co0)
    b0 = np.zeros(nn)
    np.add.at(b0, b_[a_ == 0], c0[a_ == 0] * p_in)
    np.add.at(b0, a_[b_ == 0], c0[b_ == 0] * p_in)
    np.add.at(b0, outs, co0 * po)
    p = np.zeros(nn)
    p[0] = p_in
    p[free] = spla.spsolve(J0, b0[free])
    R, m, c, mo, co = resid(p)
    nr = np.linalg.norm(R[free])
    for it in range(max_it):
        dpn = spla.spsolve(jac(c, co), R[free])
        lam = 1.0
        for _ in range(30):
            pt_ = p.copy()
            pt_[free] += lam * dpn
            Rt, mt, ct, mot, cot = resid(pt_)
            nrt = np.linalg.norm(Rt[free])
            if nrt < nr or lam < 1e-6:
                break
            lam *= 0.5
        p, R, m, c, mo, co, nr = pt_, Rt, mt, ct, mot, cot, nrt
        res = nr / max(np.max(np.abs(m)), 1e-15)
        if res < 1e-7:
            break
    dp = p[a_] - p[b_]
    m, c = _m_of_dp(dp, rho, mu, D, A, L)
    Re = np.abs(m) * D / (A * mu)
    f = f_churchill_v(Re)
    m_in = float(np.sum(m[a_ == 0]))
    return dict(p=p, m=m, Re=Re, f=f, rho=rho, mu=mu, D=D, L=L, A=A, m_in=m_in, iters=it + 1,
                residual=float(res))


def thermal(panel, net, hyd, Tc_in, p_mean):
    """связанная линейная задача: температура стенки (nx·ny) + температура охладителя в узлах элементов"""
    nx, ny = panel.nx, panel.ny
    nc = nx * ny
    P = np.array(net.nodes)
    kt = panel.k_w * panel.t_w
    rows, cols, vals = [], [], []
    rhs = np.zeros(nc)
    diag = np.zeros(nc)
    Ac = panel.dx * panel.dy
    # теплопроводность (одинаковый шаг не обязателен)
    ii, jj = np.meshgrid(np.arange(nx), np.arange(ny), indexing="ij")
    idc = ii * ny + jj
    for (di, dj, cond) in ((1, 0, kt * panel.dy / panel.dx), (0, 1, kt * panel.dx / panel.dy)):
        m = (ii + di < nx) & (jj + dj < ny)
        a = idc[m]; b = (ii[m] + di) * ny + (jj[m] + dj)
        rows += [a, b, a, b]; cols += [b, a, a, b]
        vals += [-cond * np.ones(a.size), -cond * np.ones(a.size), cond * np.ones(a.size), cond * np.ones(a.size)]
    # горячая сторона (с учётом теплопроводности наружного слоя до средней поверхности)
    hg_eff = 1.0 / (1.0 / panel.hg + panel.t_skin / panel.k_w)
    rows.append(idc.ravel()); cols.append(idc.ravel()); vals.append((hg_eff * Ac).ravel())
    rhs += (hg_eff * panel.Tg * Ac).ravel()
    if panel.rad is not None:
        rows.append(idc.ravel()); cols.append(idc.ravel()); vals.append((panel.rad_lin * Ac).ravel())
        rhs += (panel.rad_rhs * Ac).ravel()
    # элементы каналов
    elems = []           # (cell, UA, seg, k_in_node, k_out_node)
    cn = 0               # число узлов охладителя
    seg_nodes = {}
    node_cool = {}       # узел сети -> индекс узла охладителя (на конце сегмента вниз по течению)
    # порядок: от подвода по дереву (сегменты, у которых вход известен)
    order = []
    mref = max(np.max(np.abs(hyd["m"])), 1e-30)
    out_of = {}
    for s_, (a, b, d) in enumerate(net.seg):
        up, dn = (a, b) if hyd["m"][s_] >= -1e-9 * mref else (b, a)
        out_of.setdefault(up, []).append((s_, up, dn))
    queue = [0]
    seen = {0}
    while queue:
        u = queue.pop(0)
        for item in out_of.get(u, []):
            order.append(item)
            if item[2] not in seen:
                seen.add(item[2])
                queue.append(item[2])
    if len(order) != len(net.seg):
        raise RuntimeError("сеть не дерево от подвода")
    node_cool[0] = -1                      # -1 означает «вход, T = Tc_in»
    cool_rows, cool_cols, cool_vals, cool_rhs = [], [], [], []
    wall_link = []                         # (cell, cool_mid_a, cool_mid_b, UA) для стенки
    seg_info = []
    for s, up, dn in order:
        a, b, d = net.seg[s]
        mdot = abs(hyd["m"][s])
        if mdot < 1e-12:
            node_cool[dn] = node_cool[up]
            continue
        pr = air_props(np.array([panel.Tc_seg[s]]), p_mean)
        cp = float(pr["cp"][0]); k = float(pr["k"][0])
        Re = hyd["Re"][s]
        Nu = nu_pipe(Re, float(pr["Pr"][0]), hyd["f"][s])
        hc = Nu * k / d
        Lseg = hyd["L"][s]
        is_manifold = s in getattr(net, "manifold_segs", [])
        nsub = 1 if is_manifold else max(2, int(math.ceil(Lseg / (0.5 * min(panel.dx, panel.dy)))))
        A0, B0 = P[up], P[dn]
        prev = node_cool[up]
        seg_nodes[s] = []
        for q in range(nsub):
            t0, t1 = q / nsub, (q + 1) / nsub
            xm = A0 + (t0 + t1) / 2 * (B0 - A0)
            c = panel.cell(*xm)
            UA = 0.0 if is_manifold else hc * math.pi * d * Lseg / nsub
            cur = cn; cn += 1
            # ṁcp(T_cur − T_prev) = UA(T_w − (T_cur+T_prev)/2)
            cool_rows.append(cur); cool_cols.append(("c", cur)); cool_vals.append(mdot * cp + UA / 2)
            if prev >= 0:
                cool_rows.append(cur); cool_cols.append(("c", prev)); cool_vals.append(-mdot * cp + UA / 2)
                cool_rhs.append((cur, 0.0))
            else:
                cool_rhs.append((cur, (mdot * cp - UA / 2) * Tc_in))
            cool_rows.append(cur); cool_cols.append(("w", c)); cool_vals.append(-UA)
            wall_link.append((c, prev, cur, UA))
            seg_nodes[s].append(cur)
            prev = cur
        node_cool[dn] = prev
        seg_info.append((s, hc, Re, Nu))
    N = nc + cn
    R = list(np.concatenate(rows)); C = list(np.concatenate(cols)); V = list(np.concatenate(vals))
    for c, pa, pb, UA in wall_link:
        R.append(c); C.append(c); V.append(UA)
        if pa >= 0:
            R.append(c); C.append(nc + pa); V.append(-UA / 2)
        else:
            rhs[c] += UA / 2 * Tc_in
        R.append(c); C.append(nc + pb); V.append(-UA / 2)
    for r_, c_, v_ in zip(cool_rows, cool_cols, cool_vals):
        R.append(nc + r_)
        C.append(nc + c_[1] if c_[0] == "c" else c_[1])
        V.append(v_)
    b = np.concatenate([rhs, np.zeros(cn)])
    for r_, v_ in cool_rhs:
        b[nc + r_] += v_
    M = sp.coo_matrix((V, (R, C)), shape=(N, N)).tocsc()
    x = spla.spsolve(M, b)
    Tw = x[:nc].reshape(nx, ny)
    Tc = x[nc:]
    # выход охладителя: средняя по выпускам, взвешенная по расходу
    Tout = []
    into = seg_into(net)
    for o in net.outlets:
        s_last = into[o]
        Tout.append((abs(hyd["m"][s_last]), Tc[node_cool[o]] if node_cool[o] >= 0 else Tc_in))
    mo = sum(w for w, _ in Tout)
    Tc_out = sum(w * t for w, t in Tout) / mo if mo > 0 else Tc_in
    # средняя температура охладителя по сегментам — для свойств
    seg_T = np.full(len(net.seg), Tc_in)
    for s_, lst in seg_nodes.items():
        if lst:
            seg_T[s_] = float(np.mean(Tc[lst]))
    return dict(Tw=Tw, Tc_nodes=Tc, Tc_out=float(Tc_out), wall_link=wall_link, nc=nc, seg_info=seg_info,
                node_cool=node_cool, seg_T=seg_T)


def meter_outlets(panel, net, p_in, p_out, Tc_in, qcell, iters=25):
    """дозирование выпусков: K_i такие, что расход окончания ∝ теплу его площадки (ячейка Вороного × q);
    K ≥ 1 (не меньше скоростного напора выпуска). Возвращает долю K, упёршихся в 1."""
    P = np.array(net.nodes)
    outs = net.outlets
    pts = P[outs]
    Xc = np.stack([panel.X.ravel(), panel.Y.ravel()], 1)
    near = np.argmin(((Xc[:, None, :] - pts[None, :, :]) ** 2).sum(-1), axis=1)
    dem = np.bincount(near, weights=qcell.ravel() * panel.dx * panel.dy, minlength=len(outs))
    w = dem / dem.sum()
    net.K_out = {o: 1.0 for o in outs}
    last = {o: None for o in outs}
    for it in range(iters):
        p_mean = 0.5 * (p_in + (p_out(panel.Lx) if callable(p_out) else p_out))
        hyd = hydraulics(net, p_in, p_out, np.full(len(net.seg), Tc_in + 100), p_mean)
        into = seg_into(net)
        mo = np.array([abs(hyd["m"][into[o]]) for o in outs])
        tgt = w * mo.sum()
        ratio = mo / np.maximum(tgt, 1e-15)
        for k, o in enumerate(outs):
            net.K_out[o] = max(1.0, net.K_out[o] * ratio[k] ** 2)
        if np.max(np.abs(ratio - 1)[ratio > 1]) < 0.02 if np.any(ratio > 1) else True:
            break
    floor = float(np.mean([net.K_out[o] <= 1.0 + 1e-9 for o in outs]))
    dev = float(np.max(np.abs(mo / np.maximum(tgt, 1e-15) - 1)))
    return dict(iterations=it + 1, share_unmetered=floor, max_dev=dev, K_max=float(max(net.K_out.values())))


def surface_T(panel, Tw):
    """температура горячей поверхности: средняя по толщине + перепад в наружном слое"""
    hg = panel.hg
    a = panel.t_skin / panel.k_w
    return (Tw + a * hg * panel.Tg) / (1 + a * hg)


def evaluate(panel, net, p_in, p_out, Tc_in, outer=3):
    p_mean = 0.5 * (p_in + (p_out(panel.Lx) if callable(p_out) else p_out))
    panel.Tc_seg = np.full(len(net.seg), Tc_in + 50.0)
    if panel.rad is not None:
        Tw = np.full((panel.nx, panel.ny), 1100.0)
    for it in range(outer):
        if panel.rad is not None:
            eps, Tr = panel.rad
            # линеаризация излучения: q = σε(T_r⁴ − T⁴) ≈ a·(T_r − T) с a по текущей T
            Tw_ = Tw if it else np.full((panel.nx, panel.ny), 1100.0)
            a = SIGMA * eps * (Tr ** 2 + Tw_ ** 2) * (Tr + Tw_)
            panel.rad_lin = a
            panel.rad_rhs = a * Tr
        hyd = hydraulics(net, p_in, p_out, panel.Tc_seg, p_mean)
        th = thermal(panel, net, hyd, Tc_in, p_mean)
        Tw = th["Tw"]
        Tc_n = th["Tc_nodes"]
        panel.Tc_seg = th["seg_T"]          # средняя температура охладителя в сегменте — для свойств
    Ts = surface_T(panel, Tw)
    # энтропия
    Ac = panel.dx * panel.dy
    S_dT = 0.0
    Q_c = 0.0
    for c, pa, pb, UA in th["wall_link"]:
        Tcm = 0.5 * ((Tc_n[pa] if pa >= 0 else Tc_in) + Tc_n[pb])
        q = UA * (Tw.ravel()[c] - Tcm)
        Q_c += q
        S_dT += q * (1.0 / Tcm - 1.0 / Tw.ravel()[c])
    hg_eff = 1.0 / (1.0 / panel.hg + panel.t_skin / panel.k_w)
    qg = hg_eff * (panel.Tg - Tw)
    if panel.rad is not None:
        qg = qg + panel.rad_lin * (panel.rad[1] - Tw)
    Q_g = float(qg.sum() * Ac)
    S_gas = float((qg * (1.0 / Tw - 1.0 / panel.Tg)).sum() * Ac)
    S_dp = float(sum(abs(hyd["m"][s]) * abs(hyd["p"][a] - hyd["p"][b]) / (hyd["rho"][s] * panel.Tc_seg[s])
                     for s, (a, b, d) in enumerate(net.seg)))
    vmax = float(np.max(np.abs(hyd["m"]) / (hyd["rho"] * hyd["A"])))
    pr = air_props(np.array([Tc_in]), p_mean)
    a_snd = math.sqrt(1.4 * 287.05 * Tc_in)
    return dict(net=net.name, m_c=hyd["m_in"], dp=p_in - (p_out(panel.Lx) if callable(p_out) else p_out),
                Tc_out=th["Tc_out"],
                T_max=float(Ts.max()), T_min=float(Ts.min()), T_mean=float(Ts.mean()),
                T_spread=float(Ts.max() - Ts.min()), T_std=float(Ts.std()),
                Tw_mid_max=float(Tw.max()), Q_gas_W=Q_g, Q_cool_W=float(Q_c),
                S_gen_wall_to_coolant_W_K=float(S_dT), S_gen_gas_to_wall_W_K=S_gas, S_gen_friction_W_K=S_dp,
                Re_max=float(hyd["Re"].max()), Re_min=float(hyd["Re"][hyd["Re"] > 0].min()), Mach_max=vmax / a_snd,
                volume_mm3=net.volume() * 1e9, length_mm=net.length() * 1e3, n_outlets=len(net.outlets),
                _Ts=Ts, _hyd=hyd)


def p_floor(net, p_out):
    """наименьшее давление подвода без обратного течения: наибольшее давление среди выпусков"""
    if getattr(net, "p_drain", None) is not None:
        return net.p_drain
    if not callable(p_out):
        return p_out
    P = np.array(net.nodes)
    return max(p_out(P[o][0]) for o in net.outlets)


def mdot_for(panel, net, p_out, Tc_in, m_target, p_hi):
    """давление подвода (дросселирование), дающее заданный расход охладителя"""
    def r(pin):
        return evaluate(panel, net, pin, p_out, Tc_in, outer=3)["m_c"] - m_target
    lo = p_floor(net, p_out) * (1 + 1e-4)
    if r(p_hi) < 0:
        return None
    pin = brentq(r, lo, p_hi, rtol=1e-5)
    return evaluate(panel, net, pin, p_out, Tc_in, outer=3)


def min_flow_for_limit(panel, net, p_out, Tc_in, T_lim, p_hi):
    """наименьший расход охладителя, при котором T_max горячей поверхности ≤ T_lim (дросселем на подводе)"""
    full = evaluate(panel, net, p_hi, p_out, Tc_in, outer=2)
    if full["T_max"] > T_lim:
        return dict(feasible=False, at_full_dp=full)

    def r(pin):
        return evaluate(panel, net, pin, p_out, Tc_in, outer=2)["T_max"] - T_lim
    lo = p_floor(net, p_out) * (1 + 1e-5)
    if r(lo) < 0:          # предел держится даже у порога обратного течения — расход нельзя уменьшить
        res = evaluate(panel, net, lo, p_out, Tc_in, outer=3)
        res.update(feasible=True, p_in=lo, at_backflow_floor=True)
        return res
    pin = brentq(r, lo, p_hi, rtol=1e-6)
    res = evaluate(panel, net, pin, p_out, Tc_in, outer=3)
    res["feasible"] = True
    res["p_in"] = pin
    return res


# ================================================================== случаи KV-3
def kv3_point():
    """расчётная точка KV-3 — frame3_point_results.json"""
    d = json.loads((GTU / "frame3" / "frame3_point_results.json").read_text())
    cyc = d["cycle"]["frame3"]
    st1 = d["turbine"]["options"]["3"]["stages"][0]
    PR = d["split"]["pi"]
    return dict(PR=PR, p3=PR * 101325.0, T3=cyc["T3_C"] + 273.15, T41=1000.0 + 273.15, G=25.55,
                c2=st1["c2"], M2=st1["M2"], h_blade=st1["h_rotor_in_m"], Dm=d["turbine"]["options"]["3"]["Dm"],
                dp_comb=0.04, src="routes/gtu/frame3/frame3_point_results.json: split.pi, cycle.frame3.T3_C, "
                                  "turbine.options.3.Dm, turbine.options.3.stages[0].c2/M2/h_rotor_in_m; G — inputs[0].used.G")


A = dict(   # допущения уровня облика (все — оценка)
    ngv=dict(chord=0.040, s_over_c=0.80, side_len_over_c=1.2, D_LE=0.004, V1=120.0, PF=0.20, k_w=22.0, t_w=1.0e-3,
             t_skin=0.6e-3, T_lim_C=950.0, feed_loss=0.02, N_straight=16, d_straight=0.6e-3, K_tree=4, N_cco=110,
             note="хорда 40 мм, шаг 0,8 хорды, длина развёртки стороны 1,2 хорды, кромка Ø4 мм, скорость на входе 120 м/с, "
                  "коэффициент неравномерности поля температуры 0,20 (горячая полоса посередине высоты), "
                  "печатный никелевый сплав k = 22 Вт/(м·К), предел поверхности 950 °C, канал 0,6 мм — оценки"),
    liner=dict(Lx=0.10, H=0.06, Tg_pz=1900.0, Tg_exit_over_T41=150.0, h_conv=400.0, eps_rad=0.20, k_w=24.0, t_w=1.5e-3,
               t_skin=0.5e-3, T_lim_C=900.0, dp_frac=0.03, N_straight=16, d_straight=1.0e-3, K_tree=4, N_cco=160,
               note="участок жаровой трубы 100 × 60 мм у первичной зоны; газ 1900 K в первичной зоне, к выходу — T41 + 150 K; "
                    "конвекция 400 Вт/(м²·К), излучение ε = 0,20 (несветящееся пламя газа); печатный Hastelloy X k = 24; "
                    "предел 900 °C; перепад на стенке 3 % p3 — оценки"),
)


def gas_props_T(T, p):
    return air_props(np.asarray(T, float), p)   # продукты сгорания ≈ воздух по переносу (оценка)


def ngv_case(pt, T41):
    a = A["ngv"]
    Lx = a["side_len_over_c"] * a["chord"]
    H = pt["h_blade"]
    p4 = pt["p3"] * (1 - pt["dp_comb"])
    g = 1.33
    cp_g = 1180.0
    Rg = 287.0
    T3, T4 = pt["T3"], T41

    def Tg(X, Y):
        eta = (Y - H / 2) / (0.22 * H)
        sh = np.exp(-eta ** 2)
        ym = np.linspace(0, H, 400)
        sbar = np.exp(-((ym - H / 2) / (0.22 * H)) ** 2).mean()
        return T4 + a["PF"] * (T4 - T3) * (sh - sbar)

    def hg(X, Y):
        x = np.maximum(X, 1e-5)
        V = a["V1"] + (pt["c2"] - a["V1"]) * np.clip(X / Lx, 0, 1)
        Ts = T4 - V ** 2 / (2 * cp_g)
        p = p4 * (Ts / T4) ** (g / (g - 1))
        rho = p / (Rg * Ts)
        pr = gas_props_T(Ts, p4)
        mu, k, Pr = pr["mu"], pr["k"], pr["Pr"]
        Re_x = rho * V * x / mu
        h_t = 0.0296 * Re_x ** 0.8 * Pr ** (1 / 3) * k / x
        # передняя кромка: Nu_D = 1,14·Re_D^0,5·Pr^0,4 (цилиндр, точка торможения — справочная формула, не открыт)
        pr1 = gas_props_T(np.array([T4]), p4)
        rho1 = p4 / (Rg * T4)
        ReD = rho1 * a["V1"] * a["D_LE"] / pr1["mu"][0]
        hLE = 1.14 * ReD ** 0.5 * pr1["Pr"][0] ** 0.4 * pr1["k"][0] / a["D_LE"]
        xs = math.pi * a["D_LE"] / 4
        return np.where(X < xs, np.maximum(hLE, h_t), np.maximum(h_t, hLE * np.exp(-((X - xs) / xs) ** 2)))

    panel = Panel(Lx, H, 120, 70, a["k_w"], a["t_w"], hg, Tg, t_skin=a["t_skin"])
    p_in = pt["p3"] * (1 - a["feed_loss"])
    def p_out(x):
        """статическое давление на поверхности в точке x (та же изоэнтропа, что для h_g): выпуск охладителя"""
        V = a["V1"] + (pt["c2"] - a["V1"]) * min(max(x / Lx, 0.0), 1.0)
        Ts = T4 - V ** 2 / (2 * cp_g)
        return p4 * (Ts / T4) ** (g / (g - 1))
    N_v = int(round(math.pi * pt["Dm"] / (a["s_over_c"] * a["chord"])))
    return panel, p_in, p_out, N_v


def liner_case(pt):
    a = A["liner"]
    T41 = pt["T41"]

    def Tg(X, Y):
        x = X / a["Lx"]
        Te = T41 + a["Tg_exit_over_T41"]
        return np.where(x < 0.3, a["Tg_pz"], a["Tg_pz"] + (Te - a["Tg_pz"]) * (x - 0.3) / 0.7)

    def hg(X, Y):
        return np.full_like(X, a["h_conv"])
    panel = Panel(a["Lx"], a["H"], 100, 60, a["k_w"], a["t_w"], hg, Tg, t_skin=a["t_skin"])
    panel.rad = (a["eps_rad"], panel.Tg)
    p_in = pt["p3"] * 0.99
    p_out = pt["p3"] * (0.99 - a["dp_frac"])
    return panel, p_in, p_out


NET_TITLE = {"straight": "прямые каналы (K = 0)", "tree_r8": "8 деревьев, K = 1 (Мюррей)",
             "tree_r4": "4 дерева, K = 2 (Мюррей)", "tree_r2": "2 дерева, K = 3 (Мюррей)",
             "tree_r1": "1 дерево, K = 4 (Мюррей)", "tree_r4_u": "4 дерева, K = 2 (Uylings 7/3)",
             "cco": "растущее дерево (CCO)"}
NET_K = {"straight": 0, "tree_r8": 1, "tree_r4": 2, "tree_r2": 3, "tree_r1": 4}


def build_nets(panel, Ns, ds, K, Ncco, qmap, meter=None):
    """сети одного объёма каналов V = Ns·(π/4)·ds²·Lx и одного числа окончаний Ns (кроме cco — Ncco выпусков)"""
    Lx, H = panel.Lx, panel.H
    nets = {}
    s = net_straight(Lx, H, Ns, ds)
    V = Ns * math.pi * ds * ds / 4 * Lx
    nets["straight"] = s
    for name, nr, ex in (("tree_r8", 8, 3.0), ("tree_r4", 4, 3.0), ("tree_r2", 2, 3.0), ("tree_r1", 1, 3.0),
                         ("tree_r4_u", 4, 7.0 / 3.0)):
        Kk = int(round(math.log2(Ns / nr)))
        t = net_delta(Lx, H, Kk, 1e-3, ex, name, n_root=nr)
        t.scale_d(math.sqrt(V / t.volume()))
        nets[name] = t

    def qfun(P):
        return np.array([qmap(x, y) for x, y in P])
    c = net_cco(Lx, H, qfun, Ncco, 1e-4, exponent=3.0, seed=1, n_root=4)
    c.scale_d(math.sqrt(V / c.volume()))
    nets["cco"] = c
    if meter is not None:
        p_in, p_out, Tc_in, qcell = meter
        if callable(p_out):
            # лопатка: окончания сливаются в коллектор выходной кромки (давление p_out(Lx)) — «артерия + слив»;
            # вариант с выпусками в газ (эффузия) снят: перепад статического давления вдоль профиля
            # (1300 → 850 кПа) уводит охладитель к выходной кромке, дозирование при K ≥ 1 не сходится,
            # у кромки возможен заброс газа (см. README, раздел «Что не получилось»)
            c.p_drain = p_out(Lx)
        else:
            c.metering = meter_outlets(panel, c, p_in, p_out, Tc_in, qcell)
    for n in nets.values():
        n.V_target = V
    return nets


def qmap_from(panel, Tc):
    hg = panel.hg
    q = hg * (panel.Tg - Tc)
    if panel.rad is not None:
        q = q + SIGMA * panel.rad[0] * (panel.rad[1] ** 4 - (Tc + 300) ** 4)
    from scipy.interpolate import RegularGridInterpolator
    f = RegularGridInterpolator((panel.X[:, 0], panel.Y[0]), q, bounds_error=False, fill_value=None)
    return lambda x, y: float(f([[x, y]])[0])


def strip(r):
    return {k: v for k, v in r.items() if not k.startswith("_")}


def run(out_dir=HERE, quick=False):
    t0 = time.time()
    pt = kv3_point()
    res = dict(title="«Сосудистая» сеть охлаждения: сопловая лопатка и жаровая труба KV-3",
               method=__doc__.strip().split("\n\n")[0] + " Подробно — docstring vascular.py.",
               software=f"python {sys.version.split()[0]}, numpy {np.__version__}, scipy, cantera (свойства воздуха)",
               date="2026-11-14", status="оценка уровня облика: 1D-сеть + 2D-стенка; нужны сопряжённый расчёт (CFD + "
                                         "теплопроводность, задача AWS) и испытание",
               inputs=[dict(file="routes/gtu/frame3/frame3_point_results.json", used=strip(pt))],
               assumptions=A, cases={})
    figs = {}
    # ---------------- сопловая лопатка
    a = A["ngv"]
    T41s = [1000.0, 1050.0, 1100.0, 1150.0] if not quick else [1000.0, 1100.0]
    ngv = dict(T41_C=T41s, rows=[])
    panel0, p_in, p_out, N_v = ngv_case(pt, pt["T41"])
    qm = qmap_from(panel0, pt["T3"] + 150)
    p11 = ngv_case(pt, 1100.0 + 273.15)
    qcell = p11[0].hg * (p11[0].Tg - (A["ngv"]["T_lim_C"] + 273.15))
    nets = build_nets(panel0, a["N_straight"], a["d_straight"], a["K_tree"], a["N_cco"], qm,
                      meter=(p11[1], p11[2], pt["T3"], np.maximum(qcell, 0)))
    ngv.update(N_vanes_estimate=N_v, p_in_kPa=p_in / 1e3, p_out_TE_kPa=p_out(panel0.Lx) / 1e3,
               p_out_LE_kPa=p_out(0.0) / 1e3, panel_mm=[panel0.Lx * 1e3, panel0.H * 1e3],
               channel_volume_mm3=nets["straight"].V_target * 1e9,
               note="сторона лопатки — одна панель; на лопатку 2 стороны; доля охладителя x_ngv = 2·N·ṁ/G")
    for T41 in T41s:
        panel, p_in, p_out, _ = ngv_case(pt, T41 + 273.15)
        for name, net in nets.items():
            r = min_flow_for_limit(panel, net, p_out, pt["T3"], a["T_lim_C"] + 273.15, p_in)
            row = dict(T41_C=T41, net=name, feasible=r["feasible"])
            if r["feasible"] and r["Mach_max"] > 0.5:
                row["warning"] = "число Маха в канале > 0,5: несжимаемая модель вне области, расход — оценка сверху"
            if r["feasible"]:
                row.update(strip(r))
                row["x_ngv_pct"] = 100 * 2 * N_v * r["m_c"] / pt["G"]
                row["p_in_kPa"] = r["p_in"] / 1e3
                if T41 == 1050.0:
                    figs[("ngv", T41, name)] = (panel, net, r["_Ts"], True)
            else:
                row["at_full_dp"] = strip(r["at_full_dp"])
                if T41 == 1050.0:
                    figs[("ngv", T41, name)] = (panel, net, r["at_full_dp"]["_Ts"], False)
            ngv["rows"].append(row)
            print("NGV", T41, name, row.get("x_ngv_pct"), row.get("T_spread"), round(time.time() - t0), flush=True)
    # одинаковый расход (как у прямых при 1100 °C): сравнение неравномерности и энтропии
    panel, p_in, p_out, _ = ngv_case(pt, 1100.0 + 273.15)
    base = [r for r in ngv["rows"] if r["T41_C"] == 1100.0 and r["net"] == "straight" and r["feasible"]]
    eq = []
    if base:
        m_eq = base[0]["m_c"]
        for name, net in nets.items():
            r = mdot_for(panel, net, p_out, pt["T3"], m_eq, p_in)
            if r is not None:
                eq.append(strip(r))
    ngv["equal_flow_T41_1100"] = eq
    res["cases"]["ngv"] = ngv
    # ---------------- жаровая труба
    a = A["liner"]
    panel, p_in, p_out = liner_case(pt)
    qm = qmap_from(panel, pt["T3"] + 300)
    qcell = panel.hg * (panel.Tg - (a["T_lim_C"] + 273.15)) + SIGMA * a["eps_rad"] * (panel.Tg ** 4 - (a["T_lim_C"] + 273.15) ** 4)
    nets_l = build_nets(panel, a["N_straight"], a["d_straight"], a["K_tree"], a["N_cco"], qm,
                        meter=(p_in, p_out, pt["T3"], np.maximum(qcell, 0)))
    lin = dict(p_in_kPa=p_in / 1e3, p_out_kPa=p_out / 1e3, panel_mm=[a["Lx"] * 1e3, a["H"] * 1e3],
               channel_volume_mm3=nets_l["straight"].V_target * 1e9, rows=[])
    for name, net in nets_l.items():
        full = evaluate(panel, net, p_in, p_out, pt["T3"], outer=3)
        r = min_flow_for_limit(panel, net, p_out, pt["T3"], a["T_lim_C"] + 273.15, p_in)
        row = dict(net=name, at_full_dp=strip(full), feasible=r["feasible"])
        if r["feasible"]:
            row["at_limit"] = strip(r)
            row["coolant_per_m2_kg_s"] = r["m_c"] / (a["Lx"] * a["H"])
        figs[("liner", 0, name)] = (panel, net, full["_Ts"], True)
        lin["rows"].append(row)
        print("liner", name, full["T_max"] - 273.15, row.get("coolant_per_m2_kg_s"), round(time.time() - t0), flush=True)
    res["cases"]["liner"] = lin
    res["runtime_s"] = round(time.time() - t0, 1)
    return res, figs, nets, nets_l


def tube_mesh(P0, P1, r, nseg=10):
    P0, P1 = np.asarray(P0, float), np.asarray(P1, float)
    ax = P1 - P0
    L = np.linalg.norm(ax)
    if L < 1e-12:
        return np.zeros((0, 3)), np.zeros((0, 3), int)
    ax /= L
    tmp = np.array([0, 0, 1.0]) if abs(ax[2]) < 0.9 else np.array([1.0, 0, 0])
    u = np.cross(ax, tmp); u /= np.linalg.norm(u)
    v = np.cross(ax, u)
    th = np.linspace(0, 2 * np.pi, nseg, endpoint=False)
    ring = np.cos(th)[:, None] * u + np.sin(th)[:, None] * v
    V = np.vstack([P0 + r * ring, P1 + r * ring, P0, P1])
    F = []
    for i in range(nseg):
        j = (i + 1) % nseg
        F += [[i, j, nseg + i], [j, nseg + j, nseg + i], [2 * nseg, j, i], [2 * nseg + 1, nseg + i, nseg + j]]
    return V, np.array(F)


def net_stl(net, path, z=0.0):
    """каналы сети — трубки (негатив для печати; поверхность для CFD), мм"""
    from tpms_geom import write_stl
    P = np.array(net.nodes) * 1e3
    Vs, Fs, off = [], [], 0
    for a, b, d in net.seg:
        if d > 0.01:       # коллектор прямых каналов не выводим
            continue
        V, F = tube_mesh([*P[a], z], [*P[b], z], d * 1e3 / 2)
        Vs.append(V); Fs.append(F + off); off += len(V)
    return write_stl(path, np.vstack(Vs), np.vstack(Fs), solid=net.name)


def plots(res, figs, out_dir):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    names = ["straight", "tree_r4", "tree_r1", "cco"]
    title = NET_TITLE
    for art, T41 in (("ngv", 1050.0), ("liner", 0)):
        fig, axs = plt.subplots(1, 4, figsize=(17, 4.6) if art == "ngv" else (17, 3.9))
        vmin = min(figs[(art, T41, n)][2].min() for n in names if (art, T41, n) in figs) - 273.15
        vmax = min(max(figs[(art, T41, n)][2].max() for n in names if (art, T41, n) in figs) - 273.15, 1150.0)
        for ax, n in zip(axs, names):
            if (art, T41, n) not in figs:
                ax.axis("off"); continue
            panel, net, Ts, ok = figs[(art, T41, n)]
            ext = [0, panel.Lx * 1e3, 0, panel.H * 1e3]
            im = ax.imshow((Ts - 273.15).T, origin="lower", extent=ext, cmap="magma", vmin=vmin, vmax=vmax, aspect="equal")
            P = np.array(net.nodes) * 1e3
            dmax = max(d for a_, b_, d in net.seg if d < 0.01)
            for a_, b_, d in net.seg:
                if d > 0.01:
                    continue
                ax.plot([P[a_, 0], P[b_, 0]], [P[a_, 1], P[b_, 1]], color="#7fd3ff", lw=0.3 + 1.7 * d / dmax,
                        solid_capstyle="round", alpha=0.85)
            rows = res["cases"][art]["rows"]
            if art == "ngv":
                rr = [r for r in rows if r["net"] == n and r["T41_C"] == T41]
                if rr and rr[0]["feasible"]:
                    lab = f"x_ngv = {rr[0]['x_ngv_pct']:.2f} %, разброс {rr[0]['T_spread']:.0f} K"
                else:
                    lab = f"не держит 950 °C (весь перепад:\nT_max {Ts.max() - 273.15:.0f} °C)"
            else:
                rr = [r for r in rows if r["net"] == n][0]
                lab = f"T_max {rr['at_full_dp']['T_max'] - 273.15:.0f} °C, разброс {rr['at_full_dp']['T_spread']:.0f} K"
            ax.set_title(f"{title[n]}\n{lab}", fontsize=8.5)
            ax.set_xlabel("x, мм (от входной кромки)" if art == "ngv" else "x, мм (по оси камеры)", fontsize=8)
            ax.tick_params(labelsize=7)
        cb = fig.colorbar(im, ax=axs, shrink=0.85, pad=0.01)
        cb.set_label("температура горячей поверхности, °C", fontsize=8)
        fig.suptitle(("Сопловая лопатка KV-3, сторона (развёртка), T41 = 1050 °C, расход — наименьший для предела 950 °C "
                      "(у не держащих предел — весь перепад)"
                      if art == "ngv" else
                      "Жаровая труба KV-3, участок у первичной зоны, весь перепад 3 % p3"), fontsize=10)
        fig.savefig(out_dir / f"vascular_{art}.png", dpi=130, bbox_inches="tight")
        plt.close(fig)
    # x_ngv(T41)
    fig, ax = plt.subplots(figsize=(6.4, 4.2))
    mk = {"straight": "o", "tree_r8": "v", "tree_r4": "s", "tree_r2": "<", "tree_r1": "^", "tree_r4_u": "x", "cco": "D"}
    for n in NET_TITLE:
        if n not in mk:
            continue
        rr = [r for r in res["cases"]["ngv"]["rows"] if r["net"] == n and r["feasible"]]
        if not rr:
            continue
        ax.plot([r["T41_C"] for r in rr], [r["x_ngv_pct"] for r in rr], marker=mk[n], label=title[n])
    ax.axhline(3.0, color="0.5", ls="--", lw=0.8)
    ax.text(1148, 2.9, "класс N в gtu_family: x_ngv = 3 %", fontsize=8, color="0.3", ha="right", va="top")
    ax.text(1148, 0.75, "K = 3, 4 и растущее дерево\nпредел не держат", fontsize=8, color="0.3", ha="right")
    ax.set_xlabel("T41, °C"); ax.set_ylabel("x_ngv, % расхода компрессора")
    ax.set_title("Печатная лопатка с сетью каналов: расход на охлаждение (оценка)", fontsize=10)
    ax.grid(alpha=0.3); ax.legend(fontsize=8, loc="upper left")
    fig.tight_layout()
    fig.savefig(out_dir / "vascular_xngv.png", dpi=130)
    plt.close(fig)
    # развёртка по сложности (Бежан): число уровней ветвления K при равном объёме и 16 окончаниях
    fig, ax = plt.subplots(figsize=(6.4, 4.0))
    for T41, c in ((1000.0, "#7f3fbf"), (1050.0, "#1f77b4"), (1100.0, "#d62728"), (1150.0, "#2ca02c")):
        xs, ys = [], []
        for n, K in NET_K.items():
            rr = [r for r in res["cases"]["ngv"]["rows"] if r["net"] == n and r["T41_C"] == T41 and r["feasible"]]
            if rr:
                xs.append(K); ys.append(rr[0]["x_ngv_pct"])
        if xs:
            ax.plot(xs, ys, marker="o", color=c, label=f"T41 = {T41:.0f} °C")
    ax.text(3.0, ax.get_ylim()[0] + 0.85 * (ax.get_ylim()[1] - ax.get_ylim()[0]),
            "K = 3, 4: предел\nне держится\nни при каком T41", ha="center", fontsize=8, color="0.3")
    ax.set_xticks(list(NET_K.values()))
    ax.set_xticklabels(["0\n(прямые)", "1", "2", "3", "4\n(одно дерево)"])
    ax.set_xlabel("уровней ветвления K (16 окончаний, равный объём каналов)")
    ax.set_ylabel("x_ngv, % расхода компрессора")
    ax.set_title("Оптимальная сложность сети (конструктальная развёртка, оценка)", fontsize=10)
    ax.grid(alpha=0.3); ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(out_dir / "vascular_complexity.png", dpi=130)
    plt.close(fig)


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--quick", action="store_true")
    ap.add_argument("--out", default=str(HERE))
    ap.add_argument("--figs", default=None, help="куда сохранить поля для перерисовки (pickle; по умолчанию — в --out)")
    a_ = ap.parse_args()
    out = Path(a_.out)
    res, figs, nets, nets_l = run(out, quick=a_.quick)
    (out / "vascular_results.json").write_text(json.dumps(res, ensure_ascii=False, indent=1, default=float))
    try:
        import pickle
        pk = Path(a_.figs) if a_.figs else out / "vascular_figs.pkl"
        with open(pk, "wb") as fh:                # поля для перерисовки картинок без пересчёта
            pickle.dump(dict(res=res, figs=figs), fh)
    except Exception:
        pass
    try:
        plots(res, figs, out)
    except Exception as e:                      # картинка не должна терять расчёт
        res["plot_error"] = repr(e)
    stl = {}
    for name, net in nets.items():
        stl[name] = net_stl(net, out / f"vascular_ngv_{name}.stl")
    for st in stl.values():
        pth = Path(st["path"])
        st["path"] = str(pth.relative_to(ROOT)) if str(pth).startswith(str(ROOT)) else pth.name
    res["geometry"] = dict(stl=stl, note="каналы — трубки в плоскости развёртки, мм; для печати — вычитаются из стенки")
    res["networks"] = {k: dict(nodes=len(n.nodes), segments=len(n.seg), outlets=len(n.outlets),
                               d_min_mm=1e3 * min(d for _, _, d in n.seg if d < 0.01),
                               d_max_mm=1e3 * max(d for _, _, d in n.seg if d < 0.01),
                               length_mm=1e3 * n.length()) for k, n in nets.items()}
    (out / "vascular_results.json").write_text(json.dumps(res, ensure_ascii=False, indent=1, default=float))
    print("готово", res["runtime_s"])
