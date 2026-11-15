#!/usr/bin/env python3
"""Общий модуль представлений формы и поиска — локальная копия семейства «лопатки и колёса».

Почему локальная: общий модуль routes/gtu/evo/core/ в момент создания семейства строило другое семейство
(core/.owner.json: «shell-family, building»), поэтому по правилу задачи здесь — свой под-модуль, без правок чужих
файлов. Интерфейс самодостаточный: когда общий core сойдётся, функции можно свести к нему по именам.

Что внутри (только гладкие представления — урок organic/structure: воксельные «живые» топологии проверку в CalculiX
на гладкой геометрии не прошли):
  * детерминированность: seed_of / rng_of / stable_hash — зерно из ключа прогона, семейства и номера острова;
  * B-сплайны (открытые равномерные узлы, производные), кривизна, энергия плавности, проверка C2 на узлах;
  * FFD — деформация свободной формы (Sederberg & Parry 1986, базис Бернштейна) и проверка отсутствия складок
    по знаку якобиана;
  * неявные поверхности SDF с гладкими объединением и вычитанием (полиномиальный smooth-min), гироид (TPMS) с плавным
    переходом к сплошному металлу, минимальная стенка по расстоянию до оси тела (scipy.ndimage);
  * ветвление по закону Мюррея и простая L-система (для каналов — общий инструмент, здесь не обязателен);
  * проверки сечений: самопересечение ломаной, площадь, периметр, сложность силуэта;
  * MAP-Elites: архивы сеткой и CVT (Vassiliades et al. 2018), операторы iso+line (Vassiliades & Mouret 2018),
    гауссова и полиномиальная мутация, SBX (Deb & Agrawal 1995);
  * NSGA-II (Deb et al. 2002) с доминированием по ограничениям;
  * CMA-ES (Hansen, «The CMA Evolution Strategy: A Tutorial», 2016) с ограничениями [0, 1] через штраф;
  * острова: деление на части по номеру острова, независимые острова (итог не зависит от разбиения на части и числа
    процессов), слияние архивов с контролем полноты;
  * новизна (k ближайших), рендер трёх видов в монохроме (pyvista, при отсутствии — matplotlib), лист-галерея.

Формулы операторов и тестовых функций — по памяти (первоисточники в этой сессии не открывались); проверка — тесты
test_blades.py (ZDT1, сфера, Розенброк, заполнение архива на тестовой функции, повторяемость).
"""
from __future__ import annotations

import hashlib
import json
import math
import os

import numpy as np

# =========================================================================== детерминированность


def seed_of(*parts) -> int:
    """целое зерно из строкового ключа (не зависит от PYTHONHASHSEED, платформы и числа процессов)"""
    s = ":".join(str(p) for p in parts).encode("utf-8")
    return int.from_bytes(hashlib.sha256(s).digest()[:8], "little") & ((1 << 63) - 1)


def rng_of(*parts) -> np.random.Generator:
    return np.random.Generator(np.random.PCG64(seed_of(*parts)))


def _canon(o):
    if isinstance(o, dict):
        return {str(k): _canon(v) for k, v in sorted(o.items(), key=lambda kv: str(kv[0]))}
    if isinstance(o, (list, tuple)):
        return [_canon(v) for v in o]
    if isinstance(o, np.ndarray):
        return _canon(o.tolist())
    if isinstance(o, (np.floating, float)):
        v = float(o)
        if not math.isfinite(v):
            return str(v)
        return float(f"{v:.12g}")
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, (np.bool_,)):
        return bool(o)
    return o


def stable_hash(obj) -> str:
    """sha256 канонического JSON (числа округлены до 12 значащих цифр — устойчиво к шуму последнего бита при печати)"""
    return hashlib.sha256(json.dumps(_canon(obj), sort_keys=True, ensure_ascii=False).encode("utf-8")).hexdigest()


def files_hash(paths) -> str:
    h = hashlib.sha256()
    for p in sorted(paths):
        with open(p, "rb") as f:
            h.update(os.path.basename(p).encode() + b"\0" + f.read())
    return h.hexdigest()


# =========================================================================== B-сплайны


def clamped_knots(n_ctrl: int, p: int) -> np.ndarray:
    """открытый равномерный узловой вектор: p+1 нулей, внутренние равномерно, p+1 единиц"""
    if n_ctrl <= p:
        raise ValueError("число управляющих точек должно быть больше степени")
    inner = np.linspace(0.0, 1.0, n_ctrl - p + 1)[1:-1]
    return np.concatenate([np.zeros(p + 1), inner, np.ones(p + 1)])


def bspline_basis(n_ctrl: int, p: int, t, knots=None) -> np.ndarray:
    """матрица базисных функций N[i, j] = N_j,p(t_i) (Кокс — де Бур, векторно); t в [0, 1]"""
    t = np.atleast_1d(np.asarray(t, float))
    U = clamped_knots(n_ctrl, p) if knots is None else np.asarray(knots, float)
    m = len(U) - 1
    tt = np.clip(t, U[0], U[-1])
    # степень 0: полуинтервалы; правый конец относится к последнему ненулевому интервалу
    N = np.zeros((len(tt), m))
    for j in range(m):
        N[:, j] = (tt >= U[j]) & (tt < U[j + 1])
    last = np.max(np.nonzero(U[1:] > U[:-1])[0])
    N[tt >= U[-1], :] = 0.0
    N[tt >= U[-1], last] = 1.0
    for k in range(1, p + 1):
        Nn = np.zeros((len(tt), m - k))
        for j in range(m - k):
            d1 = U[j + k] - U[j]
            d2 = U[j + k + 1] - U[j + 1]
            a = (tt - U[j]) / d1 * N[:, j] if d1 > 0 else 0.0
            b = (U[j + k + 1] - tt) / d2 * N[:, j + 1] if d2 > 0 else 0.0
            Nn[:, j] = a + b
        N = Nn
    return N[:, :n_ctrl]


def bspline_eval(ctrl, t, p: int = 3, der: int = 0) -> np.ndarray:
    """точки (и производные порядка der) открытого равномерного B-сплайна; ctrl (n, d) или (n,)"""
    P = np.asarray(ctrl, float)
    one_d = P.ndim == 1
    if one_d:
        P = P[:, None]
    U = clamped_knots(len(P), p)
    q = p
    for _ in range(der):
        n = len(P)
        Q = np.zeros((n - 1, P.shape[1]))
        for i in range(n - 1):
            d = U[i + q + 1] - U[i + 1]
            Q[i] = q * (P[i + 1] - P[i]) / d if d > 0 else 0.0
        P, U, q = Q, U[1:-1], q - 1
    if q < 0:
        out = np.zeros((len(np.atleast_1d(t)), P.shape[1]))
    else:
        out = bspline_basis(len(P), q, t, U) @ P
    return out[:, 0] if one_d else out


def curvature(d1: np.ndarray, d2: np.ndarray) -> np.ndarray:
    """кривизна по первой и второй производным (2D со знаком, 3D — модуль)"""
    if d1.shape[1] == 2:
        num = d1[:, 0] * d2[:, 1] - d1[:, 1] * d2[:, 0]
        return num / np.maximum(np.linalg.norm(d1, axis=1), 1e-300) ** 3
    cr = np.cross(d1, d2)
    return np.linalg.norm(cr, axis=1) / np.maximum(np.linalg.norm(d1, axis=1), 1e-300) ** 3


def spline_fairness(ctrl, p: int = 3, n: int = 400) -> dict:
    """мера плавности кривой: длина L, наибольшая |κ|·L, энергия плавности E = L³·∫(dκ/ds)² ds (безразмерная),
    скачок кривизны на внутренних узлах (слева/справа) относительно max|κ| — для p ≥ 3 равен нулю (C2)"""
    t = np.linspace(0.0, 1.0, n)
    d1 = bspline_eval(ctrl, t, p, 1)
    d2 = bspline_eval(ctrl, t, p, 2)
    sp = np.linalg.norm(d1, axis=1)
    k = curvature(d1, d2)
    ds = 0.5 * (sp[1:] + sp[:-1]) * np.diff(t)
    L = float(ds.sum())
    dk = np.diff(k) / np.maximum(ds, 1e-300)
    E = float(L ** 3 * np.sum(dk ** 2 * ds)) if L > 0 else 0.0
    kmax = float(np.max(np.abs(k))) if len(k) else 0.0
    U = clamped_knots(len(np.atleast_1d(np.asarray(ctrl)[:, 0] if np.asarray(ctrl).ndim > 1 else ctrl)), p)
    inner = np.unique(U[(U > 0) & (U < 1)])
    jump = 0.0
    if p >= 2 and len(inner):
        eps = 1e-7
        for u in inner:
            a1, a2 = bspline_eval(ctrl, [u - eps], p, 1), bspline_eval(ctrl, [u - eps], p, 2)
            b1, b2 = bspline_eval(ctrl, [u + eps], p, 1), bspline_eval(ctrl, [u + eps], p, 2)
            if a1.ndim == 1:
                a1, a2, b1, b2 = (np.c_[np.ones(1), x] for x in (a1, a2, b1, b2))
            jump = max(jump, abs(float(curvature(a1, a2)[0] - curvature(b1, b2)[0])))
    return dict(L=L, kmaxL=kmax * L, fairness_E=E, c2_jump_rel=jump / max(kmax, 1e-300))


# =========================================================================== FFD


def bernstein(n: int, t) -> np.ndarray:
    t = np.atleast_1d(np.asarray(t, float))
    i = np.arange(n + 1)
    c = np.array([math.comb(n, k) for k in i], float)
    return c[None, :] * t[:, None] ** i[None, :] * (1 - t[:, None]) ** (n - i)[None, :]


class FFD:
    """деформация свободной формы (Sederberg & Parry 1986): решётка (l+1)(m+1)(n+1) управляющих точек в коробке
    [lo, hi], смещения dP (l+1, m+1, n+1, 3) в метрах. Отображение гладкое (полиномиальное); складок нет, если
    якобиан положителен всюду (check_jacobian)."""

    def __init__(self, lo, hi, shape=(2, 2, 2)):
        self.lo = np.asarray(lo, float)
        self.hi = np.asarray(hi, float)
        self.shape = tuple(int(s) for s in shape)
        self.dP = np.zeros(tuple(s + 1 for s in self.shape) + (3,))

    def local(self, X):
        return (np.asarray(X, float) - self.lo) / (self.hi - self.lo)

    def deform(self, X) -> np.ndarray:
        X = np.asarray(X, float)
        s = np.clip(self.local(X), 0.0, 1.0)
        l, m, n = self.shape
        Bu, Bv, Bw = bernstein(l, s[:, 0]), bernstein(m, s[:, 1]), bernstein(n, s[:, 2])
        disp = np.einsum("pi,pj,pk,ijkd->pd", Bu, Bv, Bw, self.dP)
        inside = np.all((self.local(X) >= -1e-12) & (self.local(X) <= 1 + 1e-12), axis=1)
        return X + disp * inside[:, None]

    def check_jacobian(self, n: int = 5, h: float = 1e-4) -> float:
        """наименьший определитель якобиана отображения на решётке n³ точек внутри коробки (относительно исходного);
        > 0 — нет складок и самопересечений поверхностей, деформированных этим отображением"""
        g = np.linspace(0.02, 0.98, n)
        S = np.array(np.meshgrid(g, g, g, indexing="ij")).reshape(3, -1).T
        X = self.lo + S * (self.hi - self.lo)
        dX = (self.hi - self.lo) * h
        J = np.zeros((len(X), 3, 3))
        for a in range(3):
            e = np.zeros(3)
            e[a] = dX[a]
            J[:, :, a] = (self.deform(X + e) - self.deform(X - e)) / (2 * dX[a])
        return float(np.min(np.linalg.det(J)))


# =========================================================================== неявные поверхности (SDF) и решётки


def smin(a, b, k):
    """гладкий минимум (полиномиальный, Quilez): скругление радиуса ~k на стыке"""
    h = np.clip(0.5 + 0.5 * (b - a) / k, 0.0, 1.0)
    return b * (1 - h) + a * h - k * h * (1 - h)


def smooth_union(a, b, k):
    return smin(a, b, k)


def smooth_subtract(a, b, k):
    """a минус b со скруглением k"""
    return -smin(-a, b, k)


def smooth_intersect(a, b, k):
    return -smin(-a, -b, k)


def sd_sphere(X, c, r):
    return np.linalg.norm(np.asarray(X) - np.asarray(c), axis=-1) - r


def sd_capsule(X, a, b, r):
    X, a, b = np.asarray(X, float), np.asarray(a, float), np.asarray(b, float)
    ab = b - a
    t = np.clip(np.einsum("...d,d->...", X - a, ab) / max(np.dot(ab, ab), 1e-300), 0.0, 1.0)
    return np.linalg.norm(X - (a + t[..., None] * ab), axis=-1) - r


def sd_torus(X, R, r):
    """тор вокруг оси z"""
    X = np.asarray(X, float)
    q = np.stack([np.hypot(X[..., 0], X[..., 1]) - R, X[..., 2]], axis=-1)
    return np.linalg.norm(q, axis=-1) - r


def sd_box(X, c, half):
    q = np.abs(np.asarray(X, float) - np.asarray(c)) - np.asarray(half)
    return np.linalg.norm(np.maximum(q, 0.0), axis=-1) + np.minimum(np.max(q, axis=-1), 0.0)


def gyroid(X, cell):
    """поле гироида (TPMS) с периодом cell; лист — |g| < t"""
    k = 2 * np.pi / cell
    x, y, z = (np.asarray(X)[..., i] * k for i in range(3))
    return np.sin(x) * np.cos(y) + np.sin(y) * np.cos(z) + np.sin(z) * np.cos(x)


def tpms_graded(sdf_body, g, t_sheet, blend):
    """тело с гироидом внутри и сплошной оболочкой: решётка (|g| - t) плавно переходит в металл у поверхности на ширине
    blend (гладкое пересечение, без ступенек вокселей). sdf_body — расстояние до поверхности тела (< 0 внутри)"""
    lat = np.abs(g) - t_sheet
    w = np.clip(-sdf_body / max(blend, 1e-12), 0.0, 1.0)      # 0 у поверхности, 1 вглубь
    inner = (1 - w) * (-1.0) + w * lat                          # у поверхности — сплошной металл
    return np.maximum(sdf_body, inner)


def grid_points(lo, hi, n):
    lo, hi = np.asarray(lo, float), np.asarray(hi, float)
    ax = [np.linspace(lo[i], hi[i], n[i] if np.ndim(n) else n) for i in range(3)]
    G = np.stack(np.meshgrid(*ax, indexing="ij"), axis=-1)
    return G, [a[1] - a[0] for a in ax]


def min_wall_from_field(phi, h):
    """минимальная толщина стенки тела φ < 0 на равномерной сетке шага h: 2·(наименьшее расстояние до границы в точках
    оси тела — локальных максимумах евклидова расстояния внутри). Возврат (t_min, число точек оси)"""
    from scipy import ndimage
    inside = phi < 0
    if not inside.any():
        return 0.0, 0
    d = ndimage.distance_transform_edt(inside, sampling=h)
    mx = ndimage.maximum_filter(d, size=3)
    ridge = inside & (d >= mx - 1e-12) & (d > 0)
    # точки оси у самой границы (d ≤ h) — артефакт дискретизации острых кромок: оставляем, это и есть тонкие места
    if not ridge.any():
        return 0.0, 0
    return float(2 * d[ridge].min()), int(ridge.sum())


# =========================================================================== ветвление: Мюррей, L-система


def murray_children(r_parent, n_children=2, k=3.0, split=None):
    """радиусы дочерних ветвей по закону Мюррея r^k = Σ r_i^k (k = 3 — ламинарный минимум работы, Murray 1926;
    для турбулентного газа organic/structure нашёл k ≈ 2,41). split — доли расхода (по умолчанию поровну)"""
    split = np.full(n_children, 1.0 / n_children) if split is None else np.asarray(split, float) / np.sum(split)
    return list(r_parent * split ** (1.0 / k))


def lsystem(axiom: str, rules: dict, n: int) -> str:
    s = axiom
    for _ in range(n):
        s = "".join(rules.get(c, c) for c in s)
    return s


# =========================================================================== сечения и силуэты


def polyline_self_intersects(P, closed=True) -> bool:
    """есть ли пересечение несмежных отрезков ломаной (векторно, O(n²))"""
    P = np.asarray(P, float)
    if closed:
        P = np.vstack([P, P[:1]])
    A, B = P[:-1], P[1:]
    n = len(A)
    if n < 4:
        return False
    d = B - A

    def cross(u, v):
        return u[..., 0] * v[..., 1] - u[..., 1] * v[..., 0]
    i, j = np.triu_indices(n, k=2)
    if closed:
        keep = ~((i == 0) & (j == n - 1))
        i, j = i[keep], j[keep]
    r, s = d[i], d[j]
    qp = A[j] - A[i]
    den = cross(r, s)
    ok = np.abs(den) > 1e-300
    t = np.where(ok, cross(qp, s) / np.where(ok, den, 1), -1)
    u = np.where(ok, cross(qp, r) / np.where(ok, den, 1), -1)
    hit = ok & (t > 1e-12) & (t < 1 - 1e-12) & (u > 1e-12) & (u < 1 - 1e-12)
    return bool(hit.any())


def polygon_area(P) -> float:
    P = np.asarray(P, float)
    x, y = P[:, 0], P[:, 1]
    return 0.5 * float(np.dot(x, np.roll(y, -1)) - np.dot(y, np.roll(x, -1)))


def polygon_perimeter(P, closed=True) -> float:
    P = np.asarray(P, float)
    Q = np.vstack([P, P[:1]]) if closed else P
    return float(np.sum(np.linalg.norm(np.diff(Q, axis=0), axis=1)))


def silhouette_complexity(P) -> dict:
    """сложность силуэта (замкнутая ломаная): изопериметрическое отношение P/(2√(πA)) (1 — круг) и число перегибов
    (смен знака кривизны), нормированное на 100 точек контура"""
    P = np.asarray(P, float)
    A = abs(polygon_area(P))
    per = polygon_perimeter(P)
    iso = per / (2 * math.sqrt(math.pi * A)) if A > 0 else float("inf")
    d1 = np.roll(P, -1, axis=0) - P
    cr = d1[:, 0] * np.roll(d1, -1, axis=0)[:, 1] - d1[:, 1] * np.roll(d1, -1, axis=0)[:, 0]
    sg = np.sign(cr[np.abs(cr) > 1e-12 * max(np.max(np.abs(cr)), 1e-300)])
    infl = int(np.sum(sg != np.roll(sg, 1))) if len(sg) else 0
    return dict(isoperimetric=iso, inflections=infl, inflections_per100=100.0 * infl / max(len(P), 1))


# =========================================================================== операторы изменчивости (гены в [0, 1])


def clip01(X):
    return np.clip(X, 0.0, 1.0)


def mut_gauss(rng, X, sigma):
    return clip01(X + sigma * rng.standard_normal(X.shape))


def iso_line(rng, X, Y, s_iso=0.02, s_line=0.2):
    """iso+line (Vassiliades & Mouret 2018): x' = x + σ_iso·N(0, I) + σ_line·N(0, 1)·(y − x)"""
    return clip01(X + s_iso * rng.standard_normal(X.shape) + s_line * rng.standard_normal((X.shape[0], 1)) * (Y - X))


def sbx(rng, X, Y, eta=15.0, p_gene=0.5):
    """имитация двоичного скрещивания (SBX, Deb & Agrawal 1995) в [0, 1]; возвращает двух потомков"""
    u = rng.random(X.shape)
    beta = np.where(u <= 0.5, (2 * u) ** (1 / (eta + 1)), (1 / (2 * (1 - u))) ** (1 / (eta + 1)))
    swap = rng.random(X.shape) < p_gene
    c1 = np.where(swap, 0.5 * ((1 + beta) * X + (1 - beta) * Y), X)
    c2 = np.where(swap, 0.5 * ((1 - beta) * X + (1 + beta) * Y), Y)
    return clip01(c1), clip01(c2)


def mut_poly(rng, X, eta=20.0, p_gene=None):
    """полиномиальная мутация Деба в [0, 1]"""
    n = X.shape[-1]
    p_gene = 1.0 / n if p_gene is None else p_gene
    u = rng.random(X.shape)
    d1, d2 = X, 1 - X
    mpow = 1 / (eta + 1)
    dl = np.where(u < 0.5,
                  (2 * u + (1 - 2 * u) * (1 - d1) ** (eta + 1)) ** mpow - 1,
                  1 - (2 * (1 - u) + 2 * (u - 0.5) * (1 - d2) ** (eta + 1)) ** mpow)
    m = rng.random(X.shape) < p_gene
    return clip01(np.where(m, X + dl, X))


# =========================================================================== архивы MAP-Elites


class Archive:
    """общая часть архивов: на каждую ячейку — лучший геном (наибольшая пригодность)"""

    def __init__(self, n_cells, n_genes, n_desc):
        self.n_cells, self.n_genes, self.n_desc = n_cells, n_genes, n_desc
        self.G = np.full((n_cells, n_genes), np.nan)
        self.F = np.full(n_cells, -np.inf)
        self.D = np.full((n_cells, n_desc), np.nan)
        self.meta = [None] * n_cells
        self.n_add = self.n_improve = self.n_reject = 0

    def index(self, d) -> int:
        raise NotImplementedError

    def try_add(self, g, f, d, meta=None) -> str:
        if not np.isfinite(f):
            self.n_reject += 1
            return "rejected"
        c = self.index(np.asarray(d, float))
        if not np.isfinite(self.F[c]):
            self.n_add += 1
            st = "new"
        elif f > self.F[c] + 1e-15 or (abs(f - self.F[c]) <= 1e-15 and tuple(np.round(g, 12)) < tuple(np.round(self.G[c], 12))):
            # равенство пригодности разрешается порядком генома — слияние не зависит от порядка островов
            self.n_improve += 1
            st = "improved"
        else:
            self.n_reject += 1
            return "rejected"
        self.G[c], self.F[c], self.D[c], self.meta[c] = g, f, d, meta
        return st

    @property
    def filled(self):
        return np.isfinite(self.F)

    def coverage(self) -> float:
        return float(self.filled.mean())

    def qd_score(self, offset) -> float:
        f = self.F[self.filled]
        return float(np.sum(f - offset))

    def elites(self):
        return np.nonzero(self.filled)[0]

    def to_dict(self, keep_meta=True) -> dict:
        idx = self.elites()
        return dict(kind=self.kind, n_cells=self.n_cells, n_genes=self.n_genes, n_desc=self.n_desc, spec=self.spec(),
                    cells=[dict(c=int(c), g=self.G[c].tolist(), f=float(self.F[c]), d=self.D[c].tolist(),
                                meta=self.meta[c] if keep_meta else None) for c in idx])

    def load_cells(self, cells):
        for e in cells:
            self.try_add(np.asarray(e["g"], float), float(e["f"]), np.asarray(e["d"], float), e.get("meta"))
        return self

    def merge(self, other: "Archive") -> "Archive":
        for c in other.elites():
            self.try_add(other.G[c], other.F[c], other.D[c], other.meta[c])
        return self

    def digest(self) -> str:
        idx = self.elites()
        return stable_hash(dict(c=idx.tolist(), g=self.G[idx], f=self.F[idx]))


class GridArchive(Archive):
    kind = "grid"

    def __init__(self, bins, lo, hi, n_genes):
        self.bins = np.asarray(bins, int)
        self.lo, self.hi = np.asarray(lo, float), np.asarray(hi, float)
        super().__init__(int(np.prod(self.bins)), n_genes, len(self.bins))

    def spec(self):
        return dict(bins=self.bins.tolist(), lo=self.lo.tolist(), hi=self.hi.tolist())

    def index(self, d) -> int:
        k = np.floor((d - self.lo) / (self.hi - self.lo) * self.bins).astype(int)
        k = np.clip(k, 0, self.bins - 1)
        return int(np.ravel_multi_index(tuple(k), tuple(self.bins)))


def cvt_centroids(k, dim, n_samples=20000, iters=40, seed_key="cvt"):
    """центроиды CVT в [0, 1]^dim (Ллойд по детерминированной выборке; Vassiliades et al. 2018)"""
    r = rng_of(seed_key, k, dim)
    X = r.random((n_samples, dim))
    C = X[r.choice(n_samples, k, replace=False)].copy()
    for _ in range(iters):
        lab = _nearest(X, C)
        cnt = np.bincount(lab, minlength=k).astype(float)
        S = np.zeros_like(C)
        np.add.at(S, lab, X)
        nz = cnt > 0
        C[nz] = S[nz] / cnt[nz, None]
    return C


def _nearest(X, C, chunk=4096):
    out = np.empty(len(X), int)
    cc = np.sum(C * C, axis=1)
    for s in range(0, len(X), chunk):
        x = X[s:s + chunk]
        d = cc[None, :] - 2 * x @ C.T
        out[s:s + chunk] = np.argmin(d, axis=1)
    return out


class CVTArchive(Archive):
    kind = "cvt"

    def __init__(self, k, lo, hi, n_genes, seed_key="cvt", centroids=None):
        self.lo, self.hi = np.asarray(lo, float), np.asarray(hi, float)
        dim = len(self.lo)
        self.C = cvt_centroids(k, dim, seed_key=seed_key) if centroids is None else np.asarray(centroids, float)
        self.seed_key = seed_key
        super().__init__(len(self.C), n_genes, dim)

    def spec(self):
        return dict(k=self.n_cells, lo=self.lo.tolist(), hi=self.hi.tolist(), seed_key=self.seed_key,
                    centroids_hash=stable_hash(self.C))

    def index(self, d) -> int:
        x = np.clip((d - self.lo) / (self.hi - self.lo), 0.0, 1.0)
        return int(_nearest(x[None, :], self.C)[0])


def archive_from_dict(dd, centroids=None) -> Archive:
    sp = dd["spec"]
    if dd["kind"] == "grid":
        a = GridArchive(sp["bins"], sp["lo"], sp["hi"], dd["n_genes"])
    else:
        a = CVTArchive(sp["k"], sp["lo"], sp["hi"], dd["n_genes"], sp["seed_key"], centroids)
        if stable_hash(a.C) != sp["centroids_hash"]:
            raise ValueError("центроиды CVT не совпадают с записанными в архиве")
    return a.load_cells(dd["cells"])


# =========================================================================== MAP-Elites


def map_elites(evaluate, archive: Archive, n_genes, rng, n_evals, batch=32, n_init=None, s_iso=0.03, s_line=0.25,
               p_gauss=0.3, s_gauss=0.08, init=None, parents=None, log_every=None):
    """MAP-Elites (Mouret & Clune 2015) с операторами iso+line и гауссовой мутацией.
    evaluate(G (b, n)) -> list of dict(f, d, feasible, meta); недопустимые не попадают в архив.
    init — геномы начальной выборки (например, выбор владельца из галереи), parents — то же, но и как родители первого
    поколения. Возврат: журнал по поколениям."""
    n_init = max(batch, n_evals // 10) if n_init is None else n_init
    done = 0
    log = []
    seeds = []
    if parents is not None and len(parents):
        seeds.append(np.asarray(parents, float))
    if init is not None and len(init):
        seeds.append(np.asarray(init, float))
    queue = np.vstack(seeds) if seeds else np.zeros((0, n_genes))
    n_feas = 0
    while done < n_evals:
        b = min(batch, n_evals - done)
        if len(queue):
            G, queue = queue[:b], queue[b:]
            if len(G) < b:
                G = np.vstack([G, rng.random((b - len(G), n_genes))])
        elif done < n_init or archive.filled.sum() < 2:
            G = rng.random((b, n_genes))
        else:
            el = archive.elites()
            ia = el[rng.integers(0, len(el), b)]
            ib = el[rng.integers(0, len(el), b)]
            X, Y = archive.G[ia], archive.G[ib]
            G = iso_line(rng, X, Y, s_iso, s_line)
            m = rng.random(b) < p_gauss
            if m.any():
                G[m] = mut_gauss(rng, X[m], s_gauss)
        res = evaluate(G)
        for g, r in zip(G, res):
            if r.get("feasible", True):
                n_feas += 1
                archive.try_add(g, r["f"], r["d"], r.get("meta"))
        done += b
        if log_every is None or (done // batch) % max(1, log_every) == 0 or done >= n_evals:
            f = archive.F[archive.filled]
            log.append(dict(evals=done, feasible=n_feas, coverage=archive.coverage(), filled=int(archive.filled.sum()),
                            best=float(f.max()) if len(f) else None, mean=float(f.mean()) if len(f) else None))
    return log


# =========================================================================== NSGA-II


def nondominated_sort(F, CV=None):
    """быстрая сортировка по фронтам (Deb et al. 2002) с доминированием по ограничениям: допустимое лучше
    недопустимого, из недопустимых лучше меньшее нарушение. F (n, m) — минимизация, CV (n,) ≥ 0"""
    n = len(F)
    CV = np.zeros(n) if CV is None else np.asarray(CV, float)
    fe = CV <= 0
    le = np.all(F[:, None, :] <= F[None, :, :], axis=2)
    lt = np.any(F[:, None, :] < F[None, :, :], axis=2)
    dom_obj = le & lt
    dom = np.where(fe[:, None] & fe[None, :], dom_obj,
                   np.where(fe[:, None] & ~fe[None, :], True,
                            np.where(~fe[:, None] & fe[None, :], False, CV[:, None] < CV[None, :])))
    np.fill_diagonal(dom, False)
    n_dom = dom.sum(axis=0)
    rank = np.full(n, -1)
    front = np.nonzero(n_dom == 0)[0]
    r = 0
    fronts = []
    while len(front):
        rank[front] = r
        fronts.append(front)
        nd = n_dom.copy()
        for i in front:
            nd[dom[i]] -= 1
        n_dom = nd
        n_dom[rank >= 0] = -1
        front = np.nonzero(n_dom == 0)[0]
        r += 1
    return rank, fronts


def crowding(F):
    n, m = F.shape
    d = np.zeros(n)
    if n <= 2:
        return np.full(n, np.inf)
    for j in range(m):
        o = np.argsort(F[:, j], kind="stable")
        rng_ = F[o[-1], j] - F[o[0], j]
        d[o[0]] = d[o[-1]] = np.inf
        if rng_ > 0:
            d[o[1:-1]] += (F[o[2:], j] - F[o[:-2], j]) / rng_
    return d


def nsga2(evaluate, n_genes, rng, pop=100, gens=100, init=None, eta_c=2.0, eta_m=5.0, p_c=0.9):
    """NSGA-II. evaluate(G) -> (F (n, m) — минимизация, CV (n,) ≥ 0, meta list). Возврат: последняя популяция,
    её F, CV, ранги, meta и журнал. Индексы распределения SBX 2 и мутации 5 — по тесту ZDT1 (IGD 0,006 за 25 000 оценок;
    при 15/20 в этой реализации — 0,20: SBX без границ с обрезкой исследует слабее)"""
    P = rng.random((pop, n_genes))
    if init is not None and len(init):
        k = min(pop, len(init))
        P[:k] = np.asarray(init, float)[:k]
    F, CV, M = evaluate(P)
    log = []
    for g in range(gens):
        rank, fronts = nondominated_sort(F, CV)
        cd = np.zeros(len(P))
        for fr in fronts:
            cd[fr] = crowding(F[fr])
        # турнир по (ранг, скученность)
        a, b = rng.integers(0, len(P), (2, pop))
        win = np.where((rank[a] < rank[b]) | ((rank[a] == rank[b]) & (cd[a] > cd[b])), a, b)
        par = P[win]
        X, Y = par[0::2], par[1::2]
        cmask = rng.random(len(X)) < p_c
        c1, c2 = sbx(rng, X, Y, eta_c)
        c1 = np.where(cmask[:, None], c1, X)
        c2 = np.where(cmask[:, None], c2, Y)
        Q = mut_poly(rng, np.vstack([c1, c2])[:pop], eta_m)
        FQ, CVQ, MQ = evaluate(Q)
        R, FR, CVR, MR = np.vstack([P, Q]), np.vstack([F, FQ]), np.concatenate([CV, CVQ]), list(M) + list(MQ)
        rank, fronts = nondominated_sort(FR, CVR)
        keep = []
        for fr in fronts:
            if len(keep) + len(fr) <= pop:
                keep.extend(fr.tolist())
            else:
                cdf = crowding(FR[fr])
                o = np.lexsort((fr, -cdf))      # убывание скученности, при равенстве — индекс (детерминированно)
                keep.extend(fr[o[:pop - len(keep)]].tolist())
                break
        keep = np.array(keep)
        P, F, CV, M = R[keep], FR[keep], CVR[keep], [MR[i] for i in keep]
        log.append(dict(gen=g + 1, feasible=int(np.sum(CV <= 0)), front0=int(np.sum(nondominated_sort(F, CV)[0] == 0))))
    rank, _ = nondominated_sort(F, CV)
    return dict(P=P, F=F, CV=CV, rank=rank, meta=M, log=log)


def igd(F, ref):
    """обратное поколенческое расстояние: среднее расстояние точек эталонного фронта до найденного"""
    d = np.sqrt(((ref[:, None, :] - F[None, :, :]) ** 2).sum(-1))
    return float(d.min(axis=1).mean())


def zdt1(X):
    f1 = X[:, 0]
    g = 1 + 9 * X[:, 1:].mean(axis=1)
    return np.c_[f1, g * (1 - np.sqrt(f1 / g))]


# =========================================================================== CMA-ES


class CMAES:
    """(μ/μ_w, λ)-CMA-ES по учебнику Хансена (2016): кумуляция пути шага и ковариации, ранг-1 и ранг-μ обновления.
    Область [0, 1]^n: оценивается обрезанная точка + квадратичный штраф за выход (penalty·|x − clip(x)|²)."""

    def __init__(self, x0, sigma0, rng, lam=None, penalty=1e3):
        self.n = n = len(x0)
        self.m = np.asarray(x0, float).copy()
        self.sigma = float(sigma0)
        self.rng = rng
        self.lam = lam or 4 + int(3 * math.log(n))
        self.mu = self.lam // 2
        w = math.log(self.mu + 0.5) - np.log(np.arange(1, self.mu + 1))
        self.w = w / w.sum()
        self.mueff = 1.0 / np.sum(self.w ** 2)
        self.cc = (4 + self.mueff / n) / (n + 4 + 2 * self.mueff / n)
        self.cs = (self.mueff + 2) / (n + self.mueff + 5)
        self.c1 = 2 / ((n + 1.3) ** 2 + self.mueff)
        self.cmu = min(1 - self.c1, 2 * (self.mueff - 2 + 1 / self.mueff) / ((n + 2) ** 2 + self.mueff))
        self.damps = 1 + 2 * max(0.0, math.sqrt((self.mueff - 1) / (n + 1)) - 1) + self.cs
        self.pc = np.zeros(n)
        self.ps = np.zeros(n)
        self.C = np.eye(n)
        self.B = np.eye(n)
        self.Dg = np.ones(n)
        self.chiN = math.sqrt(n) * (1 - 1 / (4 * n) + 1 / (21 * n * n))
        self.penalty = penalty
        self.gen = 0
        self.best_x, self.best_f = self.m.copy(), np.inf

    def ask(self):
        Z = self.rng.standard_normal((self.lam, self.n))
        self._Y = Z @ (self.B * self.Dg).T
        return self.m + self.sigma * self._Y

    def penalized(self, X, f_raw):
        Xc = clip01(X)
        return np.asarray(f_raw, float) + self.penalty * np.sum((X - Xc) ** 2, axis=1)

    def tell(self, X, f):
        """f — минимизируемая величина (уже со штрафом penalized)"""
        f = np.asarray(f, float)
        o = np.argsort(f, kind="stable")
        if f[o[0]] < self.best_f:
            self.best_f, self.best_x = float(f[o[0]]), clip01(X[o[0]]).copy()
        Y = self._Y[o[:self.mu]]
        yw = self.w @ Y
        self.m = self.m + self.sigma * yw
        Cmhalf = (self.B / self.Dg) @ self.B.T
        self.ps = (1 - self.cs) * self.ps + math.sqrt(self.cs * (2 - self.cs) * self.mueff) * (Cmhalf @ yw)
        self.gen += 1
        hsig = np.linalg.norm(self.ps) / math.sqrt(1 - (1 - self.cs) ** (2 * self.gen)) / self.chiN < 1.4 + 2 / (self.n + 1)
        self.pc = (1 - self.cc) * self.pc + hsig * math.sqrt(self.cc * (2 - self.cc) * self.mueff) * yw
        rank_mu = (Y.T * self.w) @ Y
        self.C = ((1 - self.c1 - self.cmu) * self.C + self.c1 * (np.outer(self.pc, self.pc) + (1 - hsig) * self.cc * (2 - self.cc) * self.C)
                  + self.cmu * rank_mu)
        self.sigma *= math.exp((self.cs / self.damps) * (np.linalg.norm(self.ps) / self.chiN - 1))
        self.sigma = min(self.sigma, 1.0)
        self.C = np.triu(self.C) + np.triu(self.C, 1).T
        ev, self.B = np.linalg.eigh(self.C)
        self.Dg = np.sqrt(np.maximum(ev, 1e-30))


def cmaes_minimize(fun, x0, sigma0, rng, max_evals=2000, lam=None, ftarget=-np.inf):
    """fun(X (λ, n) в [0, 1]) -> f (λ,) минимизация; возврат (x_best, f_best, n_evals, журнал)"""
    es = CMAES(x0, sigma0, rng, lam)
    n_ev = 0
    log = []
    while n_ev < max_evals:
        X = es.ask()
        f = es.penalized(X, fun(clip01(X)))
        es.tell(X, f)
        n_ev += len(X)
        log.append(dict(evals=n_ev, best=es.best_f, sigma=es.sigma))
        if es.best_f <= ftarget or es.sigma < 1e-12:
            break
    return es.best_x, es.best_f, n_ev, log


# =========================================================================== острова и слияние


def islands_of_part(part: int, n_parts: int, n_islands: int):
    """острова части part из n_parts: подряд идущие номера (часть p — [p·N/P, (p+1)·N/P))"""
    if not (0 <= part < n_parts) or n_parts > n_islands:
        raise ValueError("нужно 0 ≤ part < n_parts ≤ n_islands")
    return list(range(part * n_islands // n_parts, (part + 1) * n_islands // n_parts))


def merge_island_archives(island_dicts, expected, config_hash, code_hash=None, allow_partial=False, centroids=None):
    """слияние архивов островов с контролем полноты: все ожидаемые острова есть, нет повторов, у всех один и тот же
    хэш конфигурации (и кода, если задан). Возврат (архив, отчёт). Неполный набор — ValueError, если не allow_partial"""
    seen, dup, bad_cfg, bad_code = {}, [], [], []
    for d in island_dicts:
        i = int(d["island"])
        if i in seen:
            dup.append(i)
            continue
        if d.get("config_hash") != config_hash:
            bad_cfg.append(i)
        if code_hash is not None and d.get("code_hash") != code_hash:
            bad_code.append(i)
        seen[i] = d
    exp = sorted(int(i) for i in expected)
    missing = [i for i in exp if i not in seen]
    extra = sorted(i for i in seen if i not in set(exp))
    rep = dict(expected=len(exp), received=len(seen), missing=missing, extra=extra, duplicates=sorted(set(dup)),
               config_mismatch=sorted(bad_cfg), code_mismatch=sorted(bad_code))
    rep["complete"] = not (missing or extra or dup or bad_cfg or bad_code)
    if not rep["complete"] and not allow_partial:
        raise ValueError(f"набор архивов неполон или несогласован: {rep}")
    arch = None
    n_ev = 0
    for i in sorted(seen):
        if i in bad_cfg or i in bad_code or i in extra:
            continue
        a = archive_from_dict(seen[i]["archive"], centroids)
        n_ev += int(seen[i].get("n_evals", 0))
        arch = a if arch is None else arch.merge(a)
    rep["n_evals_total"] = n_ev
    rep["filled"] = int(arch.filled.sum()) if arch is not None else 0
    rep["coverage"] = arch.coverage() if arch is not None else 0.0
    rep["digest"] = arch.digest() if arch is not None else None
    return arch, rep


# =========================================================================== новизна


def knn_novelty(X, k=5, ref=None):
    """новизна: среднее расстояние до k ближайших (в нормированном пространстве признаков). Без ref — внутри набора
    (себя не считаем); с ref — до набора ref (например, типовые формы OEM)"""
    X = np.asarray(X, float)
    R = X if ref is None else np.asarray(ref, float)
    d = np.sqrt(((X[:, None, :] - R[None, :, :]) ** 2).sum(-1))
    if ref is None:
        np.fill_diagonal(d, np.inf)
    kk = min(k, d.shape[1] - (1 if ref is None else 0))
    if kk <= 0:
        return np.zeros(len(X))
    return np.sort(d, axis=1)[:, :kk].mean(axis=1)


# =========================================================================== рендер (три вида, монохром)

VIEWS = (("сбоку", (0.0, -1.0, 0.0), (0.0, 0.0, 1.0)),          # ось вала — x, горизонтально; смотрим с −y
         ("с торца", (1.0, 0.0, 0.0), (0.0, 0.0, 1.0)),         # вдоль оси
         ("аксонометрия 3/4", None, (0.0, 0.0, 1.0)))           # азимут 55°, возвышение 30° (как эскизы воронки)
AXO_AZ, AXO_EL = 55.0, 30.0
FONT = "Liberation Serif"


def _axo_dir():
    az, el = math.radians(AXO_AZ), math.radians(AXO_EL)
    return (math.cos(el) * math.cos(az), -math.cos(el) * math.sin(az), math.sin(el))


def render_views(meshes, path, size=(560, 560), backend="auto"):
    """три вида в монохроме, каждый — отдельный PNG path.format(i); meshes — список (V (n,3), T (m,3), серый 0..1).
    pyvista (OSMesa, вне экрана) — если есть; иначе matplotlib. Возврат: список путей и имя движка"""
    if backend in ("auto", "pyvista"):
        try:
            import pyvista as pv  # noqa: F401
            return _render_pv(meshes, path, size), "pyvista"
        except ImportError:
            if backend == "pyvista":
                raise
    return _render_mpl(meshes, path, size), "matplotlib"


def _bounds(meshes):
    V = np.vstack([m[0] for m in meshes])
    return V.min(axis=0), V.max(axis=0)


def _render_pv(meshes, path, size):
    import pyvista as pv
    pv.OFF_SCREEN = True
    lo, hi = _bounds(meshes)
    c = 0.5 * (lo + hi)
    rad = 0.5 * float(np.linalg.norm(hi - lo))
    out = []
    for i, (_, d, up) in enumerate(VIEWS):
        d = _axo_dir() if d is None else d
        pl = pv.Plotter(off_screen=True, window_size=list(size))
        pl.set_background("white")
        for V, T, gray in meshes:
            faces = np.hstack([np.full((len(T), 1), 3), T]).ravel()
            mesh = pv.PolyData(np.asarray(V, float), faces)
            pl.add_mesh(mesh, color=(gray, gray, gray), smooth_shading=True, specular=0.15, ambient=0.25, diffuse=0.8)
            fe = mesh.extract_feature_edges(boundary_edges=True, feature_edges=True, manifold_edges=False,
                                            non_manifold_edges=False, feature_angle=50)
            if fe.n_points:
                pl.add_mesh(fe, color="black", line_width=1.0)
        pl.enable_parallel_projection()
        pl.camera.position = tuple(c + 4 * rad * np.asarray(d))
        pl.camera.focal_point = tuple(c)
        pl.camera.up = up
        pl.camera.parallel_scale = 1.05 * rad
        p = path.format(i)
        pl.screenshot(p)
        pl.close()
        out.append(p)
    return out


def _render_mpl(meshes, path, size):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d.art3d import Poly3DCollection
    lo, hi = _bounds(meshes)
    c = 0.5 * (lo + hi)
    rad = 0.5 * float(np.linalg.norm(hi - lo))
    out = []
    light = np.array([0.4, -0.6, 0.7])
    light /= np.linalg.norm(light)
    for i, (_, d, up) in enumerate(VIEWS):
        fig = plt.figure(figsize=(size[0] / 100, size[1] / 100), dpi=100)
        ax = fig.add_subplot(111, projection="3d", proj_type="ortho")
        for V, T, gray in meshes:
            tri = np.asarray(V)[np.asarray(T)]
            nrm = np.cross(tri[:, 1] - tri[:, 0], tri[:, 2] - tri[:, 0])
            nrm /= np.maximum(np.linalg.norm(nrm, axis=1, keepdims=True), 1e-300)
            sh = np.clip(0.35 + 0.65 * np.abs(nrm @ light), 0, 1) * gray
            pc = Poly3DCollection(tri, facecolors=np.c_[sh, sh, sh], edgecolors="none", linewidths=0)
            ax.add_collection3d(pc)
        ax.set_xlim(c[0] - rad, c[0] + rad)
        ax.set_ylim(c[1] - rad, c[1] + rad)
        ax.set_zlim(c[2] - rad, c[2] + rad)
        if d is None:
            ax.view_init(elev=AXO_EL, azim=-AXO_AZ)
        elif d[0] == 0:
            ax.view_init(elev=0, azim=-90)
        else:
            ax.view_init(elev=0, azim=0)
        ax.set_axis_off()
        p = path.format(i)
        fig.savefig(p, facecolor="white")
        plt.close(fig)
        out.append(p)
    return out


def gallery_sheet(tiles, out_path, title, footer, ncol=4, tile_px=(900, 330), dpi=100):
    """лист-галерея: tiles = [(подпись, [пути трёх видов], строки описания)]. Монохром, шрифт с засечками"""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib import image as mpimg
    n = len(tiles)
    nrow = max(1, math.ceil(n / ncol))
    W = ncol * tile_px[0]
    H = nrow * (tile_px[1] + 120) + 160
    fig = plt.figure(figsize=(W / dpi, H / dpi), dpi=dpi)
    fig.patch.set_facecolor("white")
    fig.text(0.01, 1 - 30 / H, title, family=FONT, fontsize=18, va="top")
    fig.text(0.01, 50 / H, footer, family=FONT, fontsize=11, va="top", wrap=True)
    for k, (cap, imgs, lines) in enumerate(tiles):
        r, c = divmod(k, ncol)
        x0 = c * tile_px[0] / W
        y_top = 1 - (90 + r * (tile_px[1] + 120)) / H
        w3 = tile_px[0] / 3 / W
        for j, p in enumerate(imgs[:3]):
            ax = fig.add_axes([x0 + j * w3, y_top - tile_px[1] / H, w3 * 0.98, tile_px[1] / H])
            ax.imshow(mpimg.imread(p), cmap="gray" if mpimg.imread(p).ndim == 2 else None)
            ax.set_axis_off()
        fig.text(x0 + 0.005, y_top + 4 / H, cap, family=FONT, fontsize=12.5, va="bottom")
        fig.text(x0 + 0.005, y_top - tile_px[1] / H - 6 / H, "\n".join(lines), family=FONT, fontsize=10, va="top")
    fig.savefig(out_path, facecolor="white", dpi=dpi)
    plt.close(fig)
    return out_path
