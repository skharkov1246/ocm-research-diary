"""Гладкие представления формы.

  * B-сплайны (открытые с закреплёнными концами и замкнутые периодические), кривизна, самопересечения ломаной;
  * неявные поверхности (SDF): примитивы, тела вращения по 2D-контуру, гладкие объединение/вычитание/пересечение
    (полиномиальный smooth-min — скругление радиуса k, «метод уровней со скруглениями»), оболочка;
  * решётки: листовые TPMS (гироид, Шварц P, алмаз) с переменной толщиной листа, относительная плотность по
    толщине (численно по ячейке); решётка Вороного (открытая — рёбра, закрытая — стенки);
  * деформация свободной формы B-сплайновой решёткой (FFD, Sederberg & Parry 1986), обратное отображение,
    якобиан (нет складок — нет самопересечений);
  * ветвящиеся каналы: дерево по закону Мюррея d_parent^k = Σ d_child^k (k = 3 — ламинарный, 2,33…2,41 —
    турбулентный: своё значение 2,41 — organic/structure/manifold_results.json) с гладкими (Безье) ветвями; L-система
    с черепахой на плоскости;
  * поле расстояний до ломаных на плоской сетке (периодичность по одной оси) и C²-гладкий профиль ребра со скруглением.
Источники — sources.py.
"""
import math

import numpy as np


# =========================================================================== B-сплайны
def knots_clamped(n_ctrl, p=3):
    if n_ctrl <= p:
        raise ValueError("контрольных точек должно быть больше степени")
    inner = n_ctrl - p - 1
    return np.r_[np.zeros(p + 1), np.linspace(0.0, 1.0, inner + 2)[1:-1], np.ones(p + 1)]


def _cox_de_boor(u, t, p, n_basis):
    u = np.asarray(u, float)
    m = len(t) - 1
    N = np.zeros((len(u), m))
    for i in range(m):
        N[:, i] = (t[i] <= u) & (u < t[i + 1])
    # правый конец: u = t_max относим к последнему ненулевому промежутку
    last = np.max(np.nonzero(t[:-1] < t[1:])[0])
    N[u >= t[-1], :] = 0.0
    N[u >= t[-1], last] = 1.0
    for k in range(1, p + 1):
        Nn = np.zeros((len(u), m - k))
        for i in range(m - k):
            d1, d2 = t[i + k] - t[i], t[i + k + 1] - t[i + 1]
            a = (u - t[i]) / d1 * N[:, i] if d1 > 0 else 0.0
            b = (t[i + k + 1] - u) / d2 * N[:, i + 1] if d2 > 0 else 0.0
            Nn[:, i] = a + b
        N = Nn
    return N[:, :n_basis]


def bspline_basis(u, n_ctrl, p=3, closed=False):
    """матрица базиса (len(u), n_ctrl); u ∈ [0, 1]. closed — равномерный периодический сплайн (C^{p−1} везде)"""
    u = np.clip(np.asarray(u, float), 0.0, 1.0)
    if not closed:
        return _cox_de_boor(u, knots_clamped(n_ctrl, p), p, n_ctrl)
    n = n_ctrl
    t = np.arange(n + 2 * p + 1, dtype=float)
    uu = p + np.mod(u, 1.0) * n
    Ne = _cox_de_boor(uu, t, p, n + p)
    B = Ne[:, :n].copy()
    B[:, :p] += Ne[:, n:n + p]
    return B


def bspline_eval(ctrl, u, p=3, closed=False):
    ctrl = np.asarray(ctrl, float)
    return bspline_basis(u, len(ctrl), p, closed) @ ctrl


def curvature_2d(P):
    """кривизна ломаной (N, 2) по конечным разностям относительно длины дуги; со знаком"""
    P = np.asarray(P, float)
    d = np.gradient(P, axis=0)
    s = np.sqrt((d ** 2).sum(axis=1)) + 1e-300
    T = d / s[:, None]
    dT = np.gradient(T, axis=0)
    return (T[:, 0] * dT[:, 1] - T[:, 1] * dT[:, 0]) / s


def turning_total(P):
    """полный поворот касательной ломаной, рад (сумма модулей углов между соседними отрезками)"""
    P = np.asarray(P, float)
    v = np.diff(P, axis=0)
    a = np.arctan2(v[:, 1], v[:, 0])
    da = np.diff(a)
    da = (da + np.pi) % (2 * np.pi) - np.pi
    return float(np.sum(np.abs(da)))


def polyline_self_intersections(P, closed=False, tol=1e-12):
    """число пересечений несмежных отрезков ломаной (N, 2)"""
    P = np.asarray(P, float)
    A = P if not closed else np.vstack([P, P[:1]])
    a, b = A[:-1], A[1:]
    n = len(a)
    if n < 3:
        return 0

    def orient(p, q, r):
        return (q[..., 0] - p[..., 0]) * (r[..., 1] - p[..., 1]) - (q[..., 1] - p[..., 1]) * (r[..., 0] - p[..., 0])
    i, j = np.triu_indices(n, k=2)
    if closed:
        keep = ~((i == 0) & (j == n - 1))
        i, j = i[keep], j[keep]
    p1, p2, q1, q2 = a[i], b[i], a[j], b[j]
    o1, o2 = orient(p1, p2, q1), orient(p1, p2, q2)
    o3, o4 = orient(q1, q2, p1), orient(q1, q2, p2)
    hit = (o1 * o2 < -tol) & (o3 * o4 < -tol)
    return int(np.sum(hit))


# =========================================================================== SDF: примитивы и гладкие операции
def smoothstep5(x):
    """C²-гладкая ступенька 0 → 1 на [0, 1]"""
    x = np.clip(x, 0.0, 1.0)
    return x * x * x * (10.0 - 15.0 * x + 6.0 * x * x)


def smin(a, b, k):
    """гладкий минимум (объединение со скруглением радиуса ~k), полиномиальный вариант"""
    if k <= 0:
        return np.minimum(a, b)
    h = np.maximum(k - np.abs(a - b), 0.0) / k
    return np.minimum(a, b) - h * h * k * 0.25


def smax(a, b, k):
    return -smin(-a, -b, k)


def union(a, b, k=0.0):
    return smin(a, b, k)


def subtract(a, b, k=0.0):
    """a без b (со скруглением k)"""
    return smax(a, -b, k)


def intersect(a, b, k=0.0):
    return smax(a, b, k)


def shell(d, t):
    """оболочка толщины t вокруг поверхности d = 0"""
    return np.abs(d) - 0.5 * t


def sd_sphere(P, c, r):
    return np.linalg.norm(np.asarray(P) - np.asarray(c), axis=-1) - r


def sd_box(P, c, half):
    q = np.abs(np.asarray(P) - np.asarray(c)) - np.asarray(half)
    return np.linalg.norm(np.maximum(q, 0.0), axis=-1) + np.minimum(np.max(q, axis=-1), 0.0)


def sd_capsule(P, a, b, r):
    P, a, b = np.asarray(P, float), np.asarray(a, float), np.asarray(b, float)
    ba = b - a
    h = np.clip(((P - a) @ ba) / max(ba @ ba, 1e-300), 0.0, 1.0)
    return np.linalg.norm(P - a - h[..., None] * ba, axis=-1) - r


def sd_cylinder_z(P, r, z0, z1):
    P = np.asarray(P, float)
    d_r = np.hypot(P[..., 0], P[..., 1]) - r
    d_z = np.maximum(z0 - P[..., 2], P[..., 2] - z1)
    q = np.stack([d_r, d_z], axis=-1)
    return np.linalg.norm(np.maximum(q, 0.0), axis=-1) + np.minimum(np.maximum(d_r, d_z), 0.0)


def sd_torus_z(P, R, r, z0=0.0):
    P = np.asarray(P, float)
    q = np.stack([np.hypot(P[..., 0], P[..., 1]) - R, P[..., 2] - z0], axis=-1)
    return np.linalg.norm(q, axis=-1) - r


def sd_polygon2d(Q, poly, chunk=200000):
    """знаковое расстояние от точек (N, 2) до замкнутого многоугольника (M, 2): внутри < 0"""
    Q = np.asarray(Q, float)
    poly = np.asarray(poly, float)
    a, b = poly, np.roll(poly, -1, axis=0)
    ba = b - a
    L2 = np.maximum((ba ** 2).sum(axis=1), 1e-300)
    out = np.empty(len(Q))
    for s in range(0, len(Q), chunk):
        q = Q[s:s + chunk]
        pa = q[:, None, :] - a[None, :, :]
        t = np.clip((pa * ba[None]).sum(axis=2) / L2[None], 0.0, 1.0)
        d2 = ((pa - t[..., None] * ba[None]) ** 2).sum(axis=2)
        dmin = np.sqrt(d2.min(axis=1))
        cond = (a[None, :, 1] > q[:, None, 1]) != (b[None, :, 1] > q[:, None, 1])
        xint = a[None, :, 0] + (q[:, None, 1] - a[None, :, 1]) * ba[None, :, 0] / np.where(ba[None, :, 1] == 0, 1e-300, ba[None, :, 1])
        inside = (np.sum(cond & (q[:, None, 0] < xint), axis=1) % 2) == 1
        out[s:s + chunk] = np.where(inside, -dmin, dmin)
    return out


def sd_revolve(P, sd2d):
    """тело вращения вокруг оси z: sd2d(Q (N, 2) в координатах (r, z))"""
    P = np.asarray(P, float)
    shp = P.shape[:-1]
    Q = np.stack([np.hypot(P[..., 0], P[..., 1]).ravel(), P[..., 2].ravel()], axis=1)
    return sd2d(Q).reshape(shp)


# =========================================================================== TPMS и Вороной
def tpms(P, kind, L):
    """значение и градиент неявной функции TPMS с периодом L (м)"""
    k = 2.0 * math.pi / L
    P = np.asarray(P, float)
    x, y, z = k * P[..., 0], k * P[..., 1], k * P[..., 2]
    sx, cx, sy, cy, sz, cz = np.sin(x), np.cos(x), np.sin(y), np.cos(y), np.sin(z), np.cos(z)
    if kind == "gyroid":
        f = sx * cy + sy * cz + sz * cx
        g = np.stack([cx * cy - sz * sx, -sx * sy + cy * cz, -sy * sz + cz * cx], axis=-1)
    elif kind == "schwarz_p":
        f = cx + cy + cz
        g = np.stack([-sx, -sy, -sz], axis=-1)
    elif kind == "diamond":
        f = sx * sy * sz + sx * cy * cz + cx * sy * cz + cx * cy * sz
        g = np.stack([cx * sy * sz + cx * cy * cz - sx * sy * cz - sx * cy * sz,
                      sx * cy * sz - sx * sy * cz + cx * cy * cz - cx * sy * sz,
                      sx * sy * cz - sx * cy * sz - cx * sy * sz + cx * cy * cz], axis=-1)
    else:
        raise ValueError(kind)
    return f, k * g


def tpms_sheet_sdf(P, kind, L, t):
    """приближённое расстояние до листа TPMS толщины t (м; может быть полем той же формы, что P[..., 0]):
    |f|/|∇f| − t/2 (первый порядок)"""
    f, g = tpms(P, kind, L)
    gn = np.maximum(np.linalg.norm(g, axis=-1), 1e-9 * 2 * math.pi / L)
    return np.abs(f) / gn - 0.5 * np.asarray(t)


_RHO_CACHE = {}


def tpms_density_table(kind, n=40):
    """относительная плотность листового TPMS от t/L (численно: доля точек ячейки внутри листа, сетка n³)"""
    key = (kind, n)
    if key in _RHO_CACHE:
        return _RHO_CACHE[key]
    g = (np.arange(n) + 0.5) / n
    X, Y, Z = np.meshgrid(g, g, g, indexing="ij")
    P = np.stack([X, Y, Z], axis=-1)
    f, gr = tpms(P, kind, 1.0)
    d = np.abs(f) / np.maximum(np.linalg.norm(gr, axis=-1), 1e-9)
    tl = np.linspace(0.0, 0.5, 51)
    rho = np.array([np.mean(d <= 0.5 * t) for t in tl])
    _RHO_CACHE[key] = (tl, rho)
    return tl, rho


def tpms_t_for_density(kind, rho):
    """толщина листа в долях периода для заданной относительной плотности (обратная таблица)"""
    tl, r = tpms_density_table(kind)
    r_m = np.maximum.accumulate(r + np.arange(len(r)) * 1e-12)
    return np.interp(rho, r_m, tl)


def voronoi_sdf(P, seeds, t, kind="open"):
    """решётка Вороного по семенам (S, 3): open — рёбра (стыки трёх ячеек: ≈ (F3 − F1)/2 − t/2), closed — стенки
    ((F2 − F1)/2 − t/2). Расстояния приближённые (для формы и проверок порошка, не для точной стенки)."""
    from scipy.spatial import cKDTree
    P = np.asarray(P, float)
    tree = cKDTree(np.asarray(seeds, float))
    shp = P.shape[:-1]
    d, _ = tree.query(P.reshape(-1, 3), k=3)
    if kind == "open":
        v = 0.5 * (d[:, 2] - d[:, 0])
    elif kind == "closed":
        v = 0.5 * (d[:, 1] - d[:, 0])
    else:
        raise ValueError(kind)
    return v.reshape(shp) - 0.5 * np.asarray(t)


def jittered_seeds(lo, hi, cell, rng, jitter=0.35):
    """семена Вороного: решётка с шагом cell и случайным сдвигом (детерминированно по rng)"""
    lo, hi = np.asarray(lo, float), np.asarray(hi, float)
    n = np.maximum(1, np.ceil((hi - lo) / cell).astype(int))
    g = [lo[i] + (np.arange(n[i]) + 0.5) * (hi[i] - lo[i]) / n[i] for i in range(3)]
    G = np.stack(np.meshgrid(*g, indexing="ij"), axis=-1).reshape(-1, 3)
    return G + (rng.random(G.shape) - 0.5) * 2 * jitter * ((hi - lo) / n)


# =========================================================================== FFD
class BSplineFFD:
    """деформация свободной формы: смещение D(x) = Σ B_i(u) B_j(v) B_k(w) d_ijk в коробке [lo, hi] (вне — с края)"""

    def __init__(self, lo, hi, n_ctrl=(4, 4, 4), p=3, disp=None):
        self.lo, self.hi = np.asarray(lo, float), np.asarray(hi, float)
        self.n = tuple(int(v) for v in n_ctrl)
        self.p = p
        self.d = np.zeros(self.n + (3,)) if disp is None else np.asarray(disp, float).reshape(self.n + (3,))

    def displacement(self, P, chunk=200000):
        P = np.asarray(P, float)
        shp = P.shape
        Q = P.reshape(-1, 3)
        out = np.empty_like(Q)
        for s in range(0, len(Q), chunk):
            u = (Q[s:s + chunk] - self.lo) / (self.hi - self.lo)
            Bx = bspline_basis(u[:, 0], self.n[0], self.p)
            By = bspline_basis(u[:, 1], self.n[1], self.p)
            Bz = bspline_basis(u[:, 2], self.n[2], self.p)
            out[s:s + chunk] = np.einsum("ni,nj,nk,ijkc->nc", Bx, By, Bz, self.d, optimize=True)
        return out.reshape(shp)

    def deform(self, P):
        return np.asarray(P, float) + self.displacement(P)

    def inverse(self, P, iters=12):
        """x: x + D(x) = P (простая итерация; сходится при ‖∂D/∂x‖ < 1)"""
        P = np.asarray(P, float)
        X = P.copy()
        for _ in range(iters):
            X = P - self.displacement(X)
        return X

    def jacobian_min(self, n=12, h=None):
        """наименьший определитель якобиана x → x + D(x) на сетке n³ коробки (≤ 0 — складка, самопересечение)"""
        g = [np.linspace(self.lo[i], self.hi[i], n) for i in range(3)]
        P = np.stack(np.meshgrid(*g, indexing="ij"), axis=-1).reshape(-1, 3)
        h = h or 1e-4 * float(np.max(self.hi - self.lo))
        J = np.zeros((len(P), 3, 3))
        for c in range(3):
            e = np.zeros(3); e[c] = h
            J[:, :, c] = (self.deform(P + e) - self.deform(P - e)) / (2 * h)
        return float(np.min(np.linalg.det(J)))


# =========================================================================== ветвящиеся каналы
def murray_parent(d_children, k=3.0):
    d = np.asarray(d_children, float)
    return float(np.sum(d ** k) ** (1.0 / k))


def bezier(p0, p1, p2, p3, n=12):
    t = np.linspace(0.0, 1.0, n)[:, None]
    return (1 - t) ** 3 * p0 + 3 * (1 - t) ** 2 * t * p1 + 3 * (1 - t) * t * t * p2 + t ** 3 * p3


def murray_tree(root, targets, d_leaf, k=3.0, alpha=0.5, order=None, n_pts=12):
    """дерево от корня к целям (N, dim): делим цели пополам по порядку order (например, по углу), узел — на доле alpha
    пути от родителя к центру подмножества. Диаметры — по закону Мюррея с показателем k от листьев (d_leaf).
    Ветви — кубические Безье с касательной, продолжающей ветвь родителя (C¹ в узлах).
    Возврат: список {pts (n_pts, dim), d, depth, n_leaves}."""
    root = np.asarray(root, float)
    T = np.asarray(targets, float)
    idx = list(range(len(T))) if order is None else list(order)
    segs = []
    d_leaf = np.broadcast_to(np.asarray(d_leaf, float), (len(T),))

    def build(sub, parent, tangent, depth):
        if len(sub) == 1:
            end = T[sub[0]]
            d = float(d_leaf[sub[0]])
        else:
            c = T[sub].mean(axis=0)
            end = parent + alpha * (c - parent)
            h = len(sub) // 2
            d1 = build(sub[:h], end, end - parent, depth + 1)
            d2 = build(sub[h:], end, end - parent, depth + 1)
            d = murray_parent([d1, d2], k)
        L = np.linalg.norm(end - parent)
        tn = tangent / max(np.linalg.norm(tangent), 1e-12)
        p1 = parent + tn * L / 3.0
        p2 = end - (end - parent) / 3.0
        segs.append(dict(pts=bezier(parent, p1, p2, end, n_pts), d=d, depth=depth, n_leaves=len(sub)))
        return d

    c0 = T[idx].mean(axis=0)
    build(idx, root, c0 - root, 0)
    return segs


def lsystem(axiom, rules, n):
    s = axiom
    for _ in range(n):
        s = "".join(rules.get(ch, ch) for ch in s)
    return s


def turtle2d(s, start, heading_deg, step, angle_deg, shrink=0.75, rng=None, jitter_deg=0.0, max_segments=5000):
    """черепаха на плоскости: F — шаг с рисованием, f — шаг без, + / − — поворот, [ ] — ветвление.
    Возврат: список отрезков (p0, p1, depth)."""
    x, y = float(start[0]), float(start[1])
    a = math.radians(heading_deg)
    st, depth, L = [], 0, step
    segs = []
    for ch in s:
        if ch in "Ff":
            x1, y1 = x + L * math.cos(a), y + L * math.sin(a)
            if ch == "F":
                segs.append(((x, y), (x1, y1), depth))
                if len(segs) >= max_segments:
                    break
            x, y = x1, y1
        elif ch in "+-":
            j = (rng.standard_normal() * jitter_deg) if (rng is not None and jitter_deg) else 0.0
            a += math.radians((angle_deg if ch == "+" else -angle_deg) + j)
        elif ch == "[":
            st.append((x, y, a, depth, L)); depth += 1; L *= shrink
        elif ch == "]":
            x, y, a, depth, L = st.pop()
    return segs


# =========================================================================== поле расстояний до ломаных (плоская сетка)
def segments_distance_grid(X, Y, segs, reach, period_y=None):
    """расстояние от узлов сетки (X, Y — 1D оси) до набора отрезков [(x0, y0, x1, y1, w)], только в полосе reach
    (дальше — inf). Возврат: dist (nx, ny), вес/ширина ближайшего отрезка w (nx, ny). period_y — периодичность по y."""
    X, Y = np.asarray(X, float), np.asarray(Y, float)
    D = np.full((len(X), len(Y)), np.inf)
    W = np.zeros((len(X), len(Y)))
    shifts = [0.0] if not period_y else [-period_y, 0.0, period_y]
    dx = X[1] - X[0] if len(X) > 1 else 1.0
    dy = Y[1] - Y[0] if len(Y) > 1 else 1.0
    for (x0, y0, x1, y1, w) in segs:
        for sh in shifts:
            ya, yb = y0 + sh, y1 + sh
            lo_x, hi_x = min(x0, x1) - reach, max(x0, x1) + reach
            lo_y, hi_y = min(ya, yb) - reach, max(ya, yb) + reach
            i0 = max(0, int(math.floor((lo_x - X[0]) / dx))); i1 = min(len(X), int(math.ceil((hi_x - X[0]) / dx)) + 1)
            j0 = max(0, int(math.floor((lo_y - Y[0]) / dy))); j1 = min(len(Y), int(math.ceil((hi_y - Y[0]) / dy)) + 1)
            if i0 >= i1 or j0 >= j1:
                continue
            gx, gy = np.meshgrid(X[i0:i1], Y[j0:j1], indexing="ij")
            bx, by = x1 - x0, yb - ya
            L2 = bx * bx + by * by
            t = np.clip(((gx - x0) * bx + (gy - ya) * by) / L2, 0.0, 1.0) if L2 > 0 else np.zeros_like(gx)
            d = np.hypot(gx - x0 - t * bx, gy - ya - t * by)
            blk = D[i0:i1, j0:j1]
            m = d < blk
            blk[m] = d[m]
            W[i0:i1, j0:j1][m] = w
    return D, W


def rib_profile(d, half_width, fillet):
    """C²-гладкий профиль ребра по расстоянию d до оси: 1 на полке (d ≤ hw), 0 за скруглением (d ≥ hw + fillet)"""
    return 1.0 - smoothstep5((np.asarray(d) - half_width) / max(fillet, 1e-9))
