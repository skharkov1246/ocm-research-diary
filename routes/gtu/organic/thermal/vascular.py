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

    def volume(self):
        P = np.array(self.nodes)
        return float(sum(math.pi * d * d / 4 * np.hypot(*(P[b] - P[a])) for a, b, d in self.seg))

    def scale_d(self, k):
        for s in self.seg:
            s[2] *= k

    def length(self):
        P = np.array(self.nodes)
        return float(sum(np.hypot(*(P[b] - P[a])) for a, b, d in self.seg))


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


def net_delta(Lx, H, K, d0, exponent, name):
    """двоичное дерево: уровень k в x = Lx·k/K, 2^k узлов по y; d_child = d_parent·2^(−1/exponent)"""
    n = Net(name)
    root = n.node(0.0, H / 2)
    level = [(root, H / 2, H)]
    ratio = 2 ** (-1.0 / exponent)
    d = d0
    for k in range(1, K + 1):
        x = Lx * k / K
        new = []
        for (nid, yc, span) in level:
            for sgn in (-1, 1):
                yc2 = yc + sgn * span / 4
                c = n.node(x, yc2)
                n.add(nid, c, d)
                new.append((c, yc2, span / 2))
        level = new
        d *= ratio
    n.outlets = [nid for nid, _, _ in level]
    n.manifold_segs = []
    return n


def net_cco(Lx, H, qfun, Nterm, d_term, exponent=3.0, seed=0, root=(0.0, None)):
    """растущее дерево: окончания по плотности теплового потока; подключение — к сегменту с наименьшим Σ L·Q^(2/exp... )"""
    rng = np.random.default_rng(seed)
    # выборка окончаний с плотностью ∝ q (отбор с отклонением), не ближе rmin друг к другу
    pts = []
    qs = qfun(np.array([[x, y] for x in np.linspace(0, Lx, 40) for y in np.linspace(0, H, 30)]))
    qmax = qs.max()
    rmin = 0.6 * math.sqrt(Lx * H / Nterm)
    tries = 0
    while len(pts) < Nterm and tries < 200000:
        tries += 1
        p = np.array([rng.uniform(0.15 * Lx / 10, Lx), rng.uniform(0.02 * H, 0.98 * H)])
        if rng.uniform(0, qmax) > qfun(p[None])[0]:
            continue
        if pts and min(np.hypot(*(p - q)) for q in pts) < rmin:
            continue
        pts.append(p)
    # порядок добавления — сначала ближние к подводу (как рост от корня)
    r0 = np.array([0.0, H / 2 if root[1] is None else root[1]])
    pts.sort(key=lambda p: np.hypot(*(p - r0)))
    # дерево: узлы, сегменты (родитель→потомок); число окончаний ниже сегмента
    nodes = [r0, pts[0]]
    segs = [[0, 1]]
    term_of = {1: True}

    def downstream_counts(nodes, segs):
        child = {}
        for i, (a, b) in enumerate(segs):
            child.setdefault(a, []).append(i)
        cnt = [0] * len(segs)

        def rec(i):
            a, b = segs[i]
            c = 1 if term_of.get(b) else 0
            for j in child.get(b, []):
                c += rec(j)
            cnt[i] = c
            return c
        for i, (a, b) in enumerate(segs):
            if a == 0:
                rec(i)
        return cnt

    def cost(nodes, segs):
        cnt = downstream_counts(nodes, segs)
        return sum(np.hypot(*(nodes[b] - nodes[a])) * cnt[i] ** (2.0 / exponent) for i, (a, b) in enumerate(segs))

    for p in pts[1:]:
        best = None
        for si in range(len(segs)):
            a, b = segs[si]
            A, B = nodes[a], nodes[b]
            # точка ветвления: перебор по сегменту и к окончанию (грубая оптимизация положения)
            for t in (0.3, 0.5, 0.7):
                base = A + t * (B - A)
                for w in (0.0, 0.25, 0.5):
                    bp = base + w * (p - base)
                    nn = nodes + [bp, p]
                    ib, ip = len(nodes), len(nodes) + 1
                    ss = [s[:] for s in segs]
                    ss[si] = [a, ib]
                    ss += [[ib, b], [ib, ip]]
                    term_of[ip] = True
                    c = cost(nn, ss)
                    del term_of[ip]
                    if best is None or c < best[0]:
                        best = (c, si, bp)
        c, si, bp = best
        a, b = segs[si]
        nodes += [bp, p]
        ib, ip = len(nodes) - 2, len(nodes) - 1
        segs[si] = [a, ib]
        segs += [[ib, b], [ib, ip]]
        term_of[ip] = True
    cnt = downstream_counts(nodes, segs)
    n = Net("cco")
    for q in nodes:
        n.node(float(q[0]), float(q[1]))
    for i, (a, b) in enumerate(segs):
        n.add(a, b, d_term * cnt[i] ** (1.0 / exponent))
    n.outlets = [i for i in term_of if term_of[i]]
    n.manifold_segs = []
    return n


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


def hydraulics(net, p_in, p_out, Tc_seg, p_mean, K_out=1.0, max_it=60):
    """узловые давления и расходы сегментов (несжимаемо при местной плотности; потери на выпуске K_out·ρV²/2)"""
    P = np.array(net.nodes)
    nn = len(P)
    segs = net.seg
    L = np.array([max(np.hypot(*(P[b] - P[a])), 1e-6) for a, b, d in segs])
    D = np.array([d for a, b, d in segs])
    A = math.pi * D ** 2 / 4
    pr = air_props(Tc_seg, p_mean)
    rho, mu = pr["rho"], pr["mu"]
    # начальная оценка: ламинарная проводимость
    G = A * D ** 2 / (32 * mu * L) * rho
    m = None
    outset = set(net.outlets)
    free = [i for i in range(nn) if i != 0]
    idx = {n_: k for k, n_ in enumerate(free)}
    for it in range(max_it):
        # выпускная «труба»: от узла-выпуска к давлению p_out с потерей K_out·ρV²/2 (V — в последнем сегменте)
        rows, cols, vals = [], [], []
        b = np.zeros(len(free))
        for s, (a, bb, d) in enumerate(segs):
            for (u, v) in ((a, bb), (bb, a)):
                if u == 0:
                    continue
                iu = idx[u]
                rows.append(iu); cols.append(iu); vals.append(G[s])
                if v == 0:
                    b[iu] += G[s] * p_in
                else:
                    rows.append(iu); cols.append(idx[v]); vals.append(-G[s])
        # выпуски
        Gout = np.zeros(nn)
        for o in net.outlets:
            s_last = [s for s, (a, bb, d) in enumerate(segs) if bb == o][0]
            if m is None:
                go = G[s_last] * 4
            else:
                mo = abs(m[s_last]) + 1e-12
                dpo = K_out * mo ** 2 / (2 * rho[s_last] * A[s_last] ** 2)
                go = mo / max(dpo, 1e-9)
            Gout[o] = go
            io = idx[o]
            rows.append(io); cols.append(io); vals.append(go)
            b[io] += go * p_out
        M = sp.coo_matrix((vals, (rows, cols)), shape=(len(free), len(free))).tocsc()
        pf = spla.spsolve(M, b)
        p = np.zeros(nn)
        p[0] = p_in
        for n_, k in idx.items():
            p[n_] = pf[k]
        dp = np.array([p[a] - p[bb] for a, bb, d in segs])
        m_new = G * dp
        # обновить проводимость по Churchill
        Re = np.abs(m_new) * D / (A * mu)
        f = np.array([f_churchill(r) for r in Re])
        Gn = np.where(np.abs(dp) > 1e-12,
                      np.sqrt(2 * rho * A ** 2 * D / (f * L) / np.maximum(np.abs(dp), 1e-12)), G)
        Gn = np.minimum(Gn, 1e3 * G + 1e-12)
        ch = np.max(np.abs(m_new - (m if m is not None else 0))) / max(np.max(np.abs(m_new)), 1e-15)
        m = m_new
        G = 0.5 * G + 0.5 * Gn
        if it > 3 and ch < 1e-7:
            break
    m_in = float(sum(m[s] for s, (a, bb, d) in enumerate(segs) if a == 0))
    return dict(p=p, m=m, Re=Re, f=f, rho=rho, mu=mu, D=D, L=L, A=A, m_in=m_in, iters=it + 1)


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
    known = {0}
    pending = list(range(len(net.seg)))
    while pending:
        nxt = []
        for s in pending:
            a, b, d = net.seg[s]
            up, dn = (a, b) if hyd["m"][s] >= 0 else (b, a)
            if up in known:
                order.append((s, up, dn))
                known.add(dn)
            else:
                nxt.append(s)
        if len(nxt) == len(pending):
            raise RuntimeError("сеть не дерево от подвода")
        pending = nxt
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
    for o in net.outlets:
        s_last = [s for s, (a, bb, d) in enumerate(net.seg) if bb == o][0]
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


def surface_T(panel, Tw):
    """температура горячей поверхности: средняя по толщине + перепад в наружном слое"""
    hg = panel.hg
    a = panel.t_skin / panel.k_w
    return (Tw + a * hg * panel.Tg) / (1 + a * hg)


def evaluate(panel, net, p_in, p_out, Tc_in, outer=3):
    p_mean = 0.5 * (p_in + p_out)
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
    return dict(net=net.name, m_c=hyd["m_in"], dp=p_in - p_out, Tc_out=th["Tc_out"],
                T_max=float(Ts.max()), T_min=float(Ts.min()), T_mean=float(Ts.mean()),
                T_spread=float(Ts.max() - Ts.min()), T_std=float(Ts.std()),
                Tw_mid_max=float(Tw.max()), Q_gas_W=Q_g, Q_cool_W=float(Q_c),
                S_gen_wall_to_coolant_W_K=float(S_dT), S_gen_gas_to_wall_W_K=S_gas, S_gen_friction_W_K=S_dp,
                Re_max=float(hyd["Re"].max()), Re_min=float(hyd["Re"][hyd["Re"] > 0].min()), Mach_max=vmax / a_snd,
                volume_mm3=net.volume() * 1e9, length_mm=net.length() * 1e3, n_outlets=len(net.outlets),
                _Ts=Ts, _hyd=hyd)


def mdot_for(panel, net, p_out, Tc_in, m_target, p_hi):
    """давление подвода (дросселирование), дающее заданный расход охладителя"""
    def r(pin):
        return evaluate(panel, net, pin, p_out, Tc_in, outer=1)["m_c"] - m_target
    lo = p_out * (1 + 1e-4)
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
    lo = p_out * (1 + 1e-5)
    pin = brentq(r, lo, p_hi, rtol=1e-6)
    res = evaluate(panel, net, pin, p_out, Tc_in, outer=3)
    res["feasible"] = True
    res["p_in"] = pin
    return res
