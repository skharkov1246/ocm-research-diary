#!/usr/bin/env python3
"""Геометрия семейства «лопатки и колёса»: центробежное колесо KV-3 (ступень 1 — двусторонний вход, ступень 2) и
осевые ступени турбины — цельные колёса (блиски) под печать или литьё, без центрального отверстия (frame3).

Только гладкие представления (урок organic/structure):
  * меридиональный канал — кубические кривые Безье втулки и покрывного контура (как platform/geom/impeller.py);
  * закон углов лопатки β(m) — кубический B-сплайн; закрутка (wrap) θ(m, s) — интеграл r·dθ/dm = tg β;
  * FFD (core_local.FFD, Бернштейн) — гладкая «живая» деформация средней поверхности по θ; складок нет, если якобиан > 0;
  * навал на выходе (lean/rake), стреловидность входной кромки (sweep), бугорки передней кромки (tubercles, синусоида
    по размаху — горбатый кит), серрации задней кромки (гладкие зубцы), разделительные лопатки со смещением по шагу и
    собственным выгибом у входа, тандем (разрез лопатки на два ряда со сдвигом по шагу);
  * задняя поверхность диска колеса — PCHIP по узлам (как organic/structure/profile_opt.py, проверенный выигрыш −14 %),
    «рыбья кость» — гладкие карманы между рёбрами под лопатками (sin² по окружности, окно по радиусу);
  * лопатка турбины — сечения: B-сплайн средней линии по углам входа/выхода, толщина NACA 4-го ряда со скруглённой
    кромкой, укладка центров тяжести с навалом (lean), составным навалом (bow), стреловидностью (sweep), закон закрутки
    cθ ∝ r^n (n = −1 — свободный вихрь), сужение площади, полость под печать (без охлаждения).

Все размеры — облик (оценка), не конструкторская документация. Расчётная точка — routes/gtu/frame3.
"""
from __future__ import annotations

import math

import numpy as np

import core_local as C

T_MIN = {"print": 1.0e-3, "cast": 1.5e-3}      # наименьшая стенка: LPBF / литьё по выплавляемым моделям (оценка)
T_EDGE = {"print": 0.5e-3, "cast": 0.8e-3}      # наименьшая толщина у самой кромки (скругление; оценка)
R_MIN_CURV = 0.5e-3                             # наименьший радиус кривизны кромок под печать/литьё (оценка)
T_TE_TURB = 0.8e-3                              # выходная кромка турбины (как funnel/meanline_l1.T_TE)
Q_LO = 0.45                                     # задняя поверхность ступени 2: доля остатка глубины на узле от 0,45 до 1
RMIN_H = 0.1                                    # наименьший радиус кривизны задней поверхности / местная толщина
SLOPE_MAX = 20.0                                # область применимости быстрой оценки диска: |dh/dr| ≤ 20 на 0,1…0,95·r2
                                                # (калибровка CalculiX, disk1d_calibration.json → domain)


# =========================================================================== гены


class Gene:
    __slots__ = ("name", "lo", "hi", "kind", "thr", "note")

    def __init__(self, name, lo, hi, kind="real", thr=0.0, note=""):
        self.name, self.lo, self.hi, self.kind, self.thr, self.note = name, lo, hi, kind, thr, note

    def value(self, x):
        x = float(min(1.0, max(0.0, x)))
        if self.kind == "switch":                 # x < thr — приём выключен (0), иначе lo..hi
            return 0.0 if x < self.thr else self.lo + (x - self.thr) / (1 - self.thr) * (self.hi - self.lo)
        v = self.lo + x * (self.hi - self.lo)
        if self.kind == "int":
            return int(round(v))
        if self.kind == "choice":
            return int(min(self.hi, math.floor(self.lo + x * (self.hi - self.lo + 1))))
        return v

    def encode(self, v):
        if self.kind == "switch":
            if v <= 0:
                return 0.5 * self.thr
            return self.thr + (v - self.lo) / (self.hi - self.lo) * (1 - self.thr)
        if self.kind == "choice":
            return (v - self.lo + 0.5) / (self.hi - self.lo + 1)
        return (v - self.lo) / (self.hi - self.lo)


GENES_COMP = [
    Gene("mu", 0.60, 0.85, note="коэффициент работы w/U2² (уровень 1 воронки)"),
    Gene("kcm", 0.16, 0.34, note="cm2/U2"),
    Gene("Z", 7, 15, "int", note="основные лопатки"),
    Gene("lsp", 0.35, 0.85, "switch", 0.15, "длина разделительной / основной (0 — нет)"),
    Gene("sp_off", 0.30, 0.70, note="положение разделительной по шагу (0,5 — посередине)"),
    Gene("sp_dbeta", -8.0, 8.0, note="собственный выгиб разделительной у входа, град"),
    Gene("nu", 0.30, 0.45, note="втулка/периферия на входе (колесо без отверстия)"),
    Gene("hub_a", 0.35, 0.65), Gene("hub_c", 0.35, 0.65), Gene("sh_a", 0.35, 0.65), Gene("sh_c", 0.35, 0.65),
    Gene("law1", 0.0, 0.7, note="закон углов β(m): 2-я управляющая точка"),
    Gene("law2", 0.3, 1.0, note="закон углов β(m): 3-я управляющая точка"),
    Gene("lean", -15.0, 15.0, note="навал на выходе (rake), град"),
    Gene("sweep", -0.20, 0.20, note="стреловидность входной кромки, доля хорды (+ — периферия впереди)"),
    Gene("tub_A", 0.008, 0.06, "switch", 0.25, "бугорки: амплитуда / хорда"),
    Gene("tub_n", 2, 7, "int", note="бугорки: число волн по размаху"),
    Gene("ser_A", 0.05, 0.30, "switch", 0.25, "серрации: глубина / b2"),
    Gene("ser_n", 3, 9, "int", note="серрации: число зубцов"),
    Gene("tand", 0.0, 1.0, "switch", 0.6, "тандем (вкл/выкл; значение не используется)"),
    Gene("tand_m", 0.22, 0.45, note="тандем: место разреза по хорде"),
    Gene("tand_clk", 0.05, 0.45, note="тандем: сдвиг второго ряда по шагу"),
    Gene("t_hub_k", 0.65, 1.35, note="толщина лопатки у втулки / принятая impeller.py"),
    Gene("t_tip_k", 0.65, 1.35, note="толщина у периферии / принятая"),
    Gene("ffd1", -1.0, 1.0), Gene("ffd2", -1.0, 1.0), Gene("ffd3", -1.0, 1.0), Gene("ffd4", -1.0, 1.0),
    Gene("route", 0, 1, "choice", note="0 — печать металлом, 1 — литьё"),
]
GENES_COMP2_EXTRA = [Gene(f"zb{i}", 0.0, 1.0, note="задняя поверхность: доля остатка глубины на узле PCHIP (монотонно к ободу)") for i in range(6)] + [
    Gene("boss", 0.8, 2.0, note="глубина ступицы у оси / кованой frame3"),
    Gene("fish", 0.10, 0.60, "switch", 0.3, "«рыбья кость»: глубина карманов / толщина"),
    Gene("fish_r0", 0.25, 0.55, note="карманы: от r/r2"), Gene("fish_r1", 0.60, 0.92, note="карманы: до r/r2")]

GENES_TURB = [
    Gene("phi", 0.40, 0.80, note="ca/U"), Gene("U", 380.0, 420.0, note="окружная скорость на среднем диаметре, м/с"),
    Gene("zn", 24, 72, "int"), Gene("zt", 30, 100, "int"),
    Gene("taper", 0.25, 0.90, note="площадь сечения у периферии / у корня"),
    Gene("hwall", 1.0, 3.0, "switch", 0.3, "полость с равной стенкой, мм (0 — сплошная; только печать)"),
    Gene("tmax_c", 0.12, 0.26), Gene("vortex", -1.0, 0.3, note="закон закрутки cθ ∝ r^n"),
    Gene("bow", 0.0, 0.25, note="составной навал: смещение середины / высоту"),
    Gene("lean", -15.0, 15.0, note="навал, град"), Gene("sweep", -15.0, 15.0, note="стреловидность, град"),
    Gene("tub_A", 0.008, 0.05, "switch", 0.3, "бугорки ротора: амплитуда / хорда"), Gene("tub_n", 3, 8, "int"),
    Gene("w1", 0.6, 1.3), Gene("w2", 0.6, 1.3), Gene("w3", 0.6, 1.3), Gene("w4", 0.6, 1.3),
    Gene("fish", 0.10, 0.50, "switch", 0.3, "«рыбья кость» полотна диска"),
    Gene("route", 0, 1, "choice", note="0 — печать, 1 — литьё"),
]


def genes_for(comp):
    if comp == "comp1":
        return GENES_COMP
    if comp == "comp2":
        return GENES_COMP + GENES_COMP2_EXTRA
    if comp == "turb":
        return GENES_TURB
    raise KeyError(comp)


def decode(comp, x):
    gs = genes_for(comp)
    if len(x) != len(gs):
        raise ValueError(f"{comp}: геном длины {len(x)}, нужно {len(gs)}")
    return {g.name: g.value(xi) for g, xi in zip(gs, x)}


def encode(comp, p):
    return np.array([min(1.0, max(0.0, g.encode(p[g.name]))) for g in genes_for(comp)])


# =========================================================================== мелочи


def bezier3(P, t):
    P = np.asarray(P, float)
    t = np.asarray(t, float)[:, None]
    return (1 - t) ** 3 * P[0] + 3 * (1 - t) ** 2 * t * P[1] + 3 * (1 - t) * t ** 2 * P[2] + t ** 3 * P[3]


def smoothstep(x):
    x = np.clip(x, 0.0, 1.0)
    return x * x * (3 - 2 * x)


def grid_tris(nu, nv, offset=0, flip=False):
    """треугольники структурной сетки nu × nv (узлы построчно)"""
    i, j = np.meshgrid(np.arange(nu - 1), np.arange(nv - 1), indexing="ij")
    a = (i * nv + j).ravel() + offset
    b = a + nv
    c = a + 1
    d = b + 1
    T = np.r_[np.c_[a, b, c], np.c_[c, b, d]]
    return T[:, ::-1] if flip else T


def thick_surface(X, t):
    """замкнутая оболочка тонкой лопатки: средняя поверхность X (nu, ns, 3), толщина t (nu, ns) → V, T"""
    nu, ns, _ = X.shape
    du = np.gradient(X, axis=0)
    ds = np.gradient(X, axis=1)
    n = np.cross(du, ds)
    n /= np.maximum(np.linalg.norm(n, axis=2, keepdims=True), 1e-300)
    A = X + 0.5 * t[..., None] * n
    B = X - 0.5 * t[..., None] * n
    V = np.vstack([A.reshape(-1, 3), B.reshape(-1, 3)])
    off = nu * ns
    T = [grid_tris(nu, ns), grid_tris(nu, ns, off, flip=True)]
    # кромки: u = 0, u = nu-1 (по s), s = 0, s = ns-1 (по u)
    for idx in (lambda k: k, lambda k: (nu - 1) * ns + k):
        e = np.array([idx(k) for k in range(ns)])
        T.append(np.c_[e[:-1], e[1:], e[:-1] + off])
        T.append(np.c_[e[1:], e[1:] + off, e[:-1] + off])
    for idx in (lambda k: k * ns, lambda k: k * ns + ns - 1):
        e = np.array([idx(k) for k in range(nu)])
        T.append(np.c_[e[:-1], e[:-1] + off, e[1:]])
        T.append(np.c_[e[1:], e[:-1] + off, e[1:] + off])
    return V, np.vstack(T)


def surface_mass(X, t, rho):
    """масса, момент инерции относительно оси x (r² = y² + z²) и центр тяжести тонкой лопатки"""
    a = X[1:, 1:] - X[:-1, :-1]
    b = X[:-1, 1:] - X[1:, :-1]
    dA = 0.5 * np.linalg.norm(np.cross(a, b), axis=2)
    tc = 0.25 * (t[1:, 1:] + t[:-1, :-1] + t[1:, :-1] + t[:-1, 1:])
    Xc = 0.25 * (X[1:, 1:] + X[:-1, :-1] + X[1:, :-1] + X[:-1, 1:])
    dm = rho * tc * dA
    r2 = Xc[..., 1] ** 2 + Xc[..., 2] ** 2
    M = float(dm.sum())
    return M, float((dm * r2).sum()), (Xc * dm[..., None]).sum((0, 1)) / max(M, 1e-300)


def revolve(xr, n_th=96, theta_mod=None):
    """поверхность вращения контура xr (n, 2: x, r) вокруг оси x; theta_mod(x, r, θ) → поправка x (карманы)"""
    th = np.linspace(0, 2 * np.pi, n_th, endpoint=False)
    x = np.repeat(xr[:, 0][:, None], n_th, 1)
    r = np.repeat(xr[:, 1][:, None], n_th, 1)
    if theta_mod is not None:
        x = x + theta_mod(x, r, th[None, :])
    V = np.stack([x, r * np.cos(th)[None, :], r * np.sin(th)[None, :]], axis=-1).reshape(-1, 3)
    n = len(xr)
    T = []
    for i in range(n - 1):
        for j in range(n_th):
            a, b = i * n_th + j, i * n_th + (j + 1) % n_th
            c, d = a + n_th, b + n_th
            T.append((a, c, b))
            T.append((b, c, d))
    return V, np.array(T)


def merge_meshes(parts):
    Vs, Ts, off = [], [], 0
    for V, T in parts:
        Vs.append(V)
        Ts.append(T + off)
        off += len(V)
    return np.vstack(Vs), np.vstack(Ts)


def rotate_x(V, a):
    c, s = math.cos(a), math.sin(a)
    R = np.array([[1, 0, 0], [0, c, -s], [0, s, c]])
    return V @ R.T


# =========================================================================== колесо компрессора


def meridional(g, p, nm=41):
    """втулка и покрывной контур (x, r): x = 0 — плоскость входа, втулка приходит в (Lz, r2), покрывной — в (Lz − b2, r2)"""
    r1h, r1t, r2, b2, Lz = g["r1h"], g["r1t"], g["r2"], g["b2"], g["Lax"]
    t = np.linspace(0, 1, nm)
    hub = bezier3([(0, r1h), (p["hub_a"] * Lz, r1h), (Lz, r2 - p["hub_c"] * (r2 - r1h)), (Lz, r2)], t)
    Ls = Lz - b2
    sh = bezier3([(0, r1t), (p["sh_a"] * Ls, r1t), (Ls, r2 - p["sh_c"] * (r2 - r1t)), (Ls, r2)], t)
    return hub, sh


def _cum_len(P):
    return np.r_[0.0, np.cumsum(np.linalg.norm(np.diff(P, axis=0), axis=1))]


def _ffd_theta(p, amp):
    """FFD в параметрическом пространстве (m, s, θ): решётка 3×1×1 по Бернштейну, смещаются только внутренние
    узлы по m (i = 1, 2) — вход и выход лопатки не сдвигаются"""
    F = C.FFD([0, 0, -1.0], [1, 1, 1.0], (3, 1, 1))
    for (i, j), key in zip(((1, 0), (2, 0), (1, 1), (2, 1)), ("ffd1", "ffd2", "ffd3", "ffd4")):
        F.dP[i, j, :, 2] = amp * p[key]
    return F


def impeller(ml, p, stage, nu_=29, ns=9, n_th=72, build_mesh=True, route=None):
    """геометрия колеса по итогу средней линии ml (meanline_l1.cc_stage) и генам p.
    Возврат: dict(mesh=[(V, T, серый)], mass..., profile для диска, checks, features)"""
    g = ml["g"]
    ind = ml["ind"]
    route = ("print", "cast")[int(p["route"])] if route is None else route
    tmin = T_MIN[route]
    r1h, r1t, r2, b2, Lz = g["r1h"], g["r1t"], g["r2"], g["b2"], g["Lax"]
    D2 = 2 * r2
    Z = int(round(g["Zfull_used"]))
    om = ml["n"] * math.pi / 30
    c1 = ind["c1"]
    hub, sh = meridional(g, p)
    nm = len(hub)
    s_ = np.linspace(0, 1, ns)
    MER = (1 - s_)[None, :, None] * hub[:, None, :] + s_[None, :, None] * sh[:, None, :]     # (nm, ns, 2)
    mt = np.linspace(0, 1, nm)
    # длина меридиональной линии на каждом размахе
    Lm = np.array([_cum_len(MER[:, j])[-1] for j in range(ns)])
    cumL = np.stack([_cum_len(MER[:, j]) for j in range(ns)], axis=1)                      # (nm, ns)
    rr = MER[..., 1]
    # закон углов
    beta1 = np.arctan(om * rr[0] / c1)                       # угол потока на входе по размаху (атака 0)
    beta2 = math.radians(ml["beta2b_deg"])
    law = C.bspline_eval([0.0, p["law1"], p["law2"], 1.0], mt, 3)
    beta = beta1[None, :] + (beta2 - beta1[None, :]) * law[:, None]
    dth = np.tan(beta) / np.maximum(rr, 1e-6)
    theta = np.zeros_like(rr)
    theta[1:] = np.cumsum(0.5 * (dth[1:] + dth[:-1]) * np.diff(cumL, axis=0), axis=0)
    # навал на выходе: θ периферии сдвинут относительно втулки на b2·tg(lean)/r2
    dth_te = math.tan(math.radians(p["lean"])) * b2 / r2
    theta = theta + (s_[None, :] - 0.0) * dth_te * smoothstep((mt[:, None] - 0.55) / 0.45)
    pitch = 2 * math.pi / Z
    # FFD (амплитуда — до 0,12 шага)
    F = _ffd_theta(p, 0.12 * pitch)
    pts = np.stack([np.repeat(mt[:, None], ns, 1), np.repeat(s_[None, :], nm, 0), np.zeros_like(rr)], -1).reshape(-1, 3)
    dth_ffd = (F.deform(pts)[:, 2]).reshape(nm, ns)
    theta = theta + dth_ffd
    jac = F.check_jacobian()
    # входная кромка: стреловидность + бугорки; выходная: серрации (в долях хорды m)
    sw = p["sweep"]
    m_le = np.where(sw >= 0, sw * (1 - s_), -sw * s_)
    tubA, tubn = p["tub_A"], int(p["tub_n"])
    m_le = m_le + tubA * (0.5 - 0.5 * np.cos(2 * np.pi * tubn * s_))
    serA, sern = p["ser_A"], int(p["ser_n"])
    tooth = 0.5 - 0.5 * np.cos(2 * np.pi * sern * s_)                # гладкий зубец 0..1
    m_te = 1.0 - serA * b2 / Lm * tooth
    # толщины (impeller.py: t_hub 0,011·D2, t_tip 0,0045·D2; кромки 0,4 и 0,7 от толщины)
    t_h = 0.011 * D2 * p["t_hub_k"]
    t_t = 0.0045 * D2 * p["t_tip_k"]

    def sample(m0, m1, dtheta_extra=None, theta_shift=0.0):
        """средняя поверхность участка лопатки [m0(s), m1(s)] → X (nu, ns, 3), t (nu, ns), θ, (x, r)"""
        u = np.linspace(0, 1, nu_)
        mm = m0[None, :] + (m1 - m0)[None, :] * u[:, None]
        X = np.zeros((nu_, ns, 3))
        TH = np.zeros((nu_, ns))
        XR = np.zeros((nu_, ns, 2))
        for j in range(ns):
            x = np.interp(mm[:, j], mt, MER[:, j, 0])
            r = np.interp(mm[:, j], mt, MER[:, j, 1])
            th = np.interp(mm[:, j], mt, theta[:, j]) + theta_shift
            if dtheta_extra is not None:
                th = th + dtheta_extra(mm[:, j], j)
            X[:, j] = np.c_[x, r * np.cos(th), r * np.sin(th)]
            TH[:, j] = th
            XR[:, j] = np.c_[x, r]
        tt = (t_h * (1 - s_) + t_t * s_)[None, :] * np.ones((nu_, 1))
        fle = np.minimum(1.0, 0.4 + 0.6 * u / 0.08)
        fte = np.minimum(1.0, 0.7 + 0.3 * (1 - u) / 0.05)
        tt = tt * (fle * fte)[:, None]
        return X, tt, TH, XR, mm

    pieces = []      # (имя, X, t, θ, xr, m)
    if p["tand"] > 0:
        mt_cut = p["tand_m"]
        gap = 0.03
        mA1 = np.minimum(np.full(ns, mt_cut - gap / 2), m_te)
        pieces.append(("основная, ряд 1",) + sample(m_le, np.maximum(mA1, m_le + 0.05)))
        pieces.append(("основная, ряд 2",) + sample(np.maximum(np.full(ns, mt_cut + gap / 2), m_le + 0.06), m_te,
                                                    theta_shift=p["tand_clk"] * pitch))
    else:
        pieces.append(("основная",) + sample(m_le, m_te))
    lsp = p["lsp"]
    if lsp > 0:
        ms0 = 1.0 - lsp
        db = math.radians(p["sp_dbeta"])
        rbar = 0.5 * (r1t + r2)

        def sp_extra(mm, j):
            xi = np.clip((mm - ms0) / lsp, 0, 1)
            return math.tan(db) * Lm[j] / rbar * lsp * xi * (1 - xi) ** 2
        m_sle = ms0 + tubA * (0.5 - 0.5 * np.cos(2 * np.pi * tubn * s_)) * lsp
        pieces.append(("разделительная",) + sample(np.full(ns, 0.0) + m_sle, m_te, sp_extra, p["sp_off"] * pitch))
    # ---------------- проверки лопаток
    checks = {}
    k0, k1 = max(1, int(0.05 * nu_)), nu_ - max(1, int(0.05 * nu_))
    t_all_min = min(float(pc[2][k0:k1].min()) for pc in pieces)          # стенка вне скруглений кромок
    t_edge_min = min(float(min(pc[2][0].min(), pc[2][-1].min())) for pc in pieces)
    checks["t_min_m"] = t_all_min
    checks["t_edge_min_m"] = t_edge_min
    checks["t_min_ok"] = bool(t_all_min >= tmin - 1e-12 and t_edge_min >= T_EDGE[route] - 1e-12)
    checks["ffd_jac_min"] = jac
    checks["ffd_ok"] = jac > 0.2
    # канал: ширина по шагу на каждом (m, s) у основной лопатки — шаг·r·cos β − t (с разделительной — половина шага)
    X0, t0, TH0, XR0, M0 = pieces[0][1:]
    beta_s = np.zeros_like(t0)
    for j in range(ns):
        beta_s[:, j] = np.interp(M0[:, j], mt, beta[:, j])
    zloc = np.where(M0 > (1 - lsp) if lsp > 0 else np.zeros_like(M0, bool), 2 * Z, Z)
    width = 2 * np.pi * XR0[..., 1] / zloc * np.cos(beta_s) - t0
    checks["passage_min_m"] = float(width.min())
    checks["passage_ok"] = bool(width.min() > 2 * tmin)
    # самопересечение: развёртка сечений (m·L, r·θ) — средняя линия ± t/2 на нескольких размахах
    si = False
    for name, X, tt, TH, XR, MM in pieces:
        for j in (0, ns // 2, ns - 1):
            L_ = _cum_len(XR[:, j])
            y = XR[:, j, 1] * TH[:, j]
            P = np.c_[L_, y]
            d = np.gradient(P, axis=0)
            nrm = np.c_[-d[:, 1], d[:, 0]] / np.maximum(np.linalg.norm(d, axis=1), 1e-300)[:, None]
            poly = np.vstack([P + 0.5 * tt[:, j:j + 1] * nrm, (P - 0.5 * tt[:, j:j + 1] * nrm)[::-1]])
            if C.polyline_self_intersects(poly):
                si = True
    # пересечение рядов: основная и разделительная/тандем не ближе толщины по окружности на общем m
    gap_min = np.inf
    if len(pieces) > 1:
        for a in range(len(pieces)):
            for b in range(a + 1, len(pieces)):
                Ma, THa, ta = pieces[a][5], pieces[a][3], pieces[a][2]
                Mb, THb, tb = pieces[b][5], pieces[b][3], pieces[b][2]
                for j in range(ns):
                    lo, hi = max(Ma[0, j], Mb[0, j]), min(Ma[-1, j], Mb[-1, j])
                    if hi <= lo:
                        continue
                    mq = np.linspace(lo, hi, 15)
                    da = np.interp(mq, Ma[:, j], THa[:, j])
                    db_ = np.interp(mq, Mb[:, j], THb[:, j])
                    rq = np.interp(mq, mt, rr[:, j])
                    dth = np.mod(db_ - da, pitch)
                    dist = rq * np.minimum(dth, pitch - dth) - 0.5 * (np.interp(mq, Ma[:, j], ta[:, j]) + np.interp(mq, Mb[:, j], tb[:, j]))
                    gap_min = min(gap_min, float(dist.min()))
    checks["row_gap_min_m"] = None if not np.isfinite(gap_min) else gap_min
    checks["rows_ok"] = (not np.isfinite(gap_min)) or gap_min > tmin
    checks["self_intersection"] = si
    # меридиональный канал: площадь 2πr·b по m не возрастает (без местной диффузии в меридиональном сечении)
    bw = np.linalg.norm(hub - sh, axis=1)
    area = 2 * np.pi * 0.5 * (hub[:, 1] + sh[:, 1]) * bw
    checks["area_ratio_out_in"] = float(area[-1] / area[0])
    checks["area_max_over_in"] = float(area.max() / area[0])
    # меридиональный канал: площадь нигде не больше входной на 30 % и выход/вход 0,3…1,1 (оценка; у облика impeller.py
    # с Безье 0,5/0,5 — 1,13 у ступени 1 и 1,24 у ступени 2)
    checks["area_monotone"] = bool(area.max() <= 1.30 * area[0] and 0.3 <= area[-1] / area[0] <= 1.1)
    # ---------------- диск колеса
    rho = 4430.0
    disk = impeller_disk(ml, p, stage, hub, Z, theta, mt, rr)
    # ---------------- масса и инерция лопаток (Z одинаковых наборов; ступень 1 — две стороны)
    nside = 2 if stage == 1 else 1
    Mb = Jb = 0.0
    for name, X, tt, *_ in pieces:
        m_, J_, _ = surface_mass(X, tt, rho)
        Mb += m_ * Z * nside
        Jb += J_ * Z * nside
    # распределение лопаток по радиусу для диска: hb(r) — осевая протяжённость лопаточной области
    rgrid = disk["r"]
    hb = np.zeros_like(rgrid)
    xh = np.interp(rgrid, hub[:, 1], hub[:, 0], left=0.0, right=Lz)
    for k, rv in enumerate(rgrid):
        if rv < r1h:
            continue
        if rv <= r1t:
            hb[k] = xh[k] - 0.0
        else:
            xs = np.interp(rv, sh[:, 1], sh[:, 0])
            hb[k] = max(0.0, xh[k] - xs)
    hb *= nside
    vol_col = float(np.trapezoid(hb * 2 * np.pi * rgrid, rgrid))
    phi = Mb / (rho * vol_col) if vol_col > 0 else 0.0
    out = dict(route=route, t_min_route_m=tmin, Z=Z, n_pieces=len(pieces), pieces=[pc[0] for pc in pieces],
               mass_blades_kg=Mb, J_blades=Jb, mass_disk_kg=disk["mass"], J_disk=disk["J"],
               mass_kg=Mb + disk["mass"], J_kgm2=Jb + disk["J"], blade_fill_phi=phi, hb=hb, disk=disk, checks=checks,
               wrap_shroud_deg=float(math.degrees(theta[-1, -1] - theta[0, -1])),
               wrap_hub_deg=float(math.degrees(theta[-1, 0] - theta[0, 0])), lean_deg=p["lean"], beta=beta, theta=theta)
    out["checks"]["disk_slope_max"] = disk["slope_max"]
    out["checks"]["disk_ok"] = disk["ok"]
    out["checks"]["valid"] = bool(checks["t_min_ok"] and checks["ffd_ok"] and checks["passage_ok"] and not si
                                  and checks["rows_ok"] and checks["area_monotone"] and disk["ok"])
    # ---------------- признаки формы (эффектность)
    out["features"] = impeller_features(p, ml, pieces, hub, sh, theta, mt, beta)
    if build_mesh:
        out["mesh"] = impeller_mesh(pieces, Z, disk, stage, n_th, Lz)
    return out


def impeller_disk(ml, p, stage, hub, Z, theta, mt, rr, n=241):
    """диск колеса: толщина h(r) тела вращения. Ступень 2: тракт (втулка) — линия hub, задняя поверхность — PCHIP по
    узлам (как organic/structure/profile_opt.py) между тонким пределом (стенка 10 мм под втулкой) и кованым профилем
    frame3 (+15 % глубины ступицы); «рыбья кость» — карманы между рёбрами. Ступень 1 — общее полотно между двумя
    трактами: h = 2·(Lz − x_hub) + t_rim (как frame3_fea, t_rim = 0,02·D2)"""
    from scipy.interpolate import PchipInterpolator
    g = ml["g"]
    r2, Lz = g["r2"], g["Lax"]
    t_rim = 0.02 * 2 * r2
    r = np.linspace(0.0, r2, n)
    xh = np.interp(r, hub[:, 1], hub[:, 0], left=0.0, right=Lz)      # для r < r1h — плоский носок x = 0
    ok = True
    info = {}
    if stage == 1:
        h = 2 * (Lz - xh) + t_rim
        xb = None
        fish = np.zeros_like(r)
    else:
        # база — задняя поверхность кованого профиля frame3 ступени 2 (frame3_fea_results.json → parts.comp2.variants.solid),
        # подобием по r2; узлы PCHIP — те же радиусы, что organic/structure/profile_opt.py (RC_KNOTS), в долях r2
        fb = FRAME3_BACK()
        r_knots = np.r_[0.0, 0.058, 0.07, 0.11, 0.15, 0.20, 0.25, 0.30, fb["r2"]] / fb["r2"] * r2
        base = np.interp(r_knots / r2 * fb["r2"], fb["r"], fb["depth"]) * r2 / fb["r2"]     # отступ за плоскость Lz
        xh_k = np.interp(r_knots, hub[:, 1], hub[:, 0], left=0.0, right=Lz)
        xb_min = xh_k + 0.010                     # стенка не тоньше 10 мм под трактом (как profile_opt)
        # глубина задней поверхности монотонно убывает от ступицы к ободу по построению: на каждом узле остаётся доля
        # q = Q_LO + (1 − Q_LO)·zb от глубины предыдущего узла над глубиной обода (перетяжки и обратный ход исключены:
        # именно они давали концентрацию, которую быстрая оценка не видит — disk1d_calibration.json → domain)
        dep = np.zeros_like(r_knots)
        dep[0] = dep[1] = p.get("boss", 1.0) * base[0]
        d_end = t_rim
        zb = np.array([p[f"zb{i}"] for i in range(6)])
        for k_ in range(2, 8):
            q = Q_LO + (1 - Q_LO) * zb[k_ - 2]
            dep[k_] = d_end + (dep[k_ - 1] - d_end) * q
        dep[-1] = d_end
        xb_k = Lz + dep
        xb_base = Lz + base
        # кубический B-сплайн по управляющим точкам (r_k, x_k): гладкий (C2) и монотонный при монотонных точках
        # (свойство уменьшения вариации) — без изломов в узлах, которые давал PCHIP
        cv = C.bspline_eval(np.c_[r_knots, xb_k], np.linspace(0, 1, 801), 3)
        xb = np.interp(r, cv[:, 0], cv[:, 1])
        xb = np.maximum(xb, xh + 0.010)
        h = xb - xh
        info["xb_knots"] = xb_k.tolist()
        info["xb_base_knots"] = xb_base.tolist()
        # «рыбья кость»: средняя по окружности глубина карманов = половина амплитуды sin²
        fish = np.zeros_like(r)
        if p.get("fish", 0) > 0:
            w = smoothstep((r / r2 - p["fish_r0"]) / 0.08) * smoothstep((p["fish_r1"] - r / r2) / 0.08)
            web = np.maximum(h - 0.010, 0.0)               # под трактом остаётся 10 мм
            fish = p["fish"] * web * w
    h_eff = h - 0.5 * fish
    mass, J = _mi(r, h_eff, 4430.0)
    m_ = (r > 0.1 * r2) & (r < 0.95 * r2)
    slope = float(np.abs(np.gradient(h_eff, r))[m_].max())
    info["slope_max"] = slope
    ok = bool(np.all(h_eff > 0.004) and slope <= SLOPE_MAX)
    if xb is not None:
        mb = (r > 0.17 * r2) & (r < 0.97 * r2)
        dxb = np.gradient(xb, r)
        kap = np.abs(np.gradient(dxb, r)) / (1 + dxb ** 2) ** 1.5
        info["back_dxdr_max"] = float(dxb[mb].max())
        info["back_rmin_over_h"] = float(np.min(1 / np.maximum(kap[mb], 1e-9) / np.maximum(h[mb], 1e-4)))
        ok = ok and info["back_dxdr_max"] <= 1e-6 and info["back_rmin_over_h"] >= RMIN_H
    return dict(r=r, h=h, h_eff=h_eff, xh=xh, xb=xb, fish=fish, mass=mass, J=J, ok=ok, t_rim=t_rim, info=info, Lz=Lz, slope_max=slope,
                stage=stage, Z=Z, theta_hub=np.interp(r, rr[:, 0], theta[:, 0]))


_FB = {}


def FRAME3_BACK():
    """задняя поверхность кованого профиля frame3 ступени 2: отступ за плоскость обода (z_top обода = 0) — depth(r)"""
    if not _FB:
        import json
        import os
        p = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))), "frame3", "frame3_fea_results.json")
        pr = json.load(open(p, encoding="utf-8"))["parts"]["comp2"]["variants"]["solid"]["profile"]
        r = np.asarray(pr["r"], float)
        zb = np.asarray(pr["zc"], float) - 0.5 * np.asarray(pr["t"], float)
        ztop_rim = pr["zc"][-1] + 0.5 * pr["t"][-1]
        _FB.update(r=r, depth=ztop_rim - zb, r2=float(r[-1]))
    return _FB


def _mi(r, h, rho):
    w = rho * h * 2 * np.pi * r
    return float(np.trapezoid(w, r)), float(np.trapezoid(w * r * r, r))


def impeller_features(p, ml, pieces, hub, sh, theta, mt, beta):
    """измеримые признаки вида: плавность (энергия плавности втулки, покрывного, закона углов), наименьший радиус
    кривизны входной кромки, сложность силуэта лопатки сбоку"""
    f = {}
    f["fair_hub"] = C.spline_fairness(np.array([hub[0], hub[len(hub) // 3], hub[2 * len(hub) // 3], hub[-1]]))["fairness_E"]
    law = C.spline_fairness(np.c_[np.linspace(0, 1, 4), [0.0, p["law1"], p["law2"], 1.0]])
    f["fair_law"] = law["fairness_E"]
    f["c2_jump_law"] = law["c2_jump_rel"]
    X, tt, TH, XR, MM = pieces[0][1:]
    le = XR[0]          # входная кромка (ns, 2)
    te = XR[-1]
    out = np.vstack([le, XR[1:-1, -1], te[::-1], XR[-2:0:-1, 0]])
    sc = C.silhouette_complexity(out)
    f["silhouette_iso"] = sc["isoperimetric"]
    f["silhouette_inflections"] = sc["inflections"]
    # радиус кривизны входной кромки (бугорки): по ломаной le
    d1 = np.gradient(le, axis=0)
    d2 = np.gradient(d1, axis=0)
    k = np.abs(C.curvature(d1, d2))
    f["le_rmin_m"] = float(1.0 / max(k.max(), 1e-12))
    f["wrap_shroud_deg"] = float(math.degrees(theta[-1, -1] - theta[0, -1]))
    return f


def impeller_mesh(pieces, Z, disk, stage, n_th, Lz):
    meshes = []
    blades = []
    for name, X, tt, *_ in pieces:
        V, T = thick_surface(X, tt)
        blades.append((V, T))
    Vb, Tb = merge_meshes(blades)
    allb = []
    for k in range(Z):
        allb.append((rotate_x(Vb, 2 * np.pi * k / Z), Tb))
    Vb, Tb = merge_meshes(allb)
    r = disk["r"]
    xh = disk["xh"]
    if stage == 1:
        xc = Lz + disk["t_rim"] / 2
        Vm = Vb.copy()
        Vm[:, 0] = 2 * xc - Vm[:, 0]
        Vb, Tb = merge_meshes([(Vb, Tb), (Vm, Tb[:, ::-1])])
        prof = np.vstack([np.c_[xh, r], np.c_[(2 * xc - xh)[::-1], r[::-1]]])
        Vd, Td = revolve(prof, n_th)
    else:
        xb = disk["xb"]
        fish = disk["fish"]
        thub = disk["theta_hub"]
        Zr = disk["Z"]

        def mod(x, rr, th):
            # карманы на задней поверхности: только у точек задней поверхности (вторая половина контура)
            return 0.0 * x
        prof_front = np.c_[xh, r]
        prof_back = np.c_[xb[::-1], r[::-1]]
        Vf, Tf = revolve(prof_front, n_th)
        # задняя поверхность с карманами: x = xb − fish(r)·sin²(Z(θ − θ_hub(r))/2)
        th = np.linspace(0, 2 * np.pi, n_th, endpoint=False)
        rb = r[::-1]
        xbk = xb[::-1][:, None] - (fish[::-1][:, None]) * np.sin(0.5 * Zr * (th[None, :] - thub[::-1][:, None] - np.pi / Zr)) ** 2
        Vbk = np.stack([xbk, rb[:, None] * np.cos(th)[None, :], rb[:, None] * np.sin(th)[None, :]], -1).reshape(-1, 3)
        _, Tbk = revolve(prof_back, n_th)
        rim = np.array([[xh[-1], r[-1]], [xb[-1], r[-1]]])
        Vr, Tr = revolve(rim, n_th)
        Vd, Td = merge_meshes([(Vf, Tf), (Vbk, Tbk), (Vr, Tr)])
    meshes.append((Vd, Td, 0.80))
    meshes.append((Vb, Tb, 0.62))
    return meshes


# =========================================================================== турбина: сечение и лопатка


NACA = np.array([0.2969, -0.1260, -0.3516, 0.2843, -0.1036])


def section(chord, beta_in, beta_out, tmax_c, n=41, le_ext=0.0):
    """сечение лопатки ротора в координатах (x — ось, y — окружность): средняя линия — квадратичный B-сплайн
    (касательные под углами входа β_in и выхода −β_out от оси), толщина NACA 4-го ряда со скруглённой входной и
    кромкой T_TE_TURB на выходе. Возврат: контур (2n−2, 2), средняя линия, толщины"""
    stag = 0.5 * (beta_in - beta_out)
    bx = chord * math.cos(stag)
    # квадратичный Безье: P0 (0,0), P1 на пересечении касательных, P2 (bx, y2)
    y2 = -bx * math.tan(stag) * 0 + bx * (math.tan(beta_in) - math.tan(beta_out)) * 0.5 * 0 + 0.0
    t_in, t_out = math.tan(beta_in), math.tan(beta_out)
    # P1 = (a, a·t_in), P2 = (bx, y2), касательная в P2: (bx − a, y2 − a·t_in) ∥ (1, −t_out)
    y2 = bx * 0.5 * (t_in - t_out)
    a = (y2 + bx * t_out) / (t_in + t_out) if abs(t_in + t_out) > 1e-9 else bx / 2
    P = np.array([[0, 0], [a, a * t_in], [bx, y2]])
    u = (1 - np.cos(np.linspace(0, np.pi, n))) / 2
    cam = (1 - u)[:, None] ** 2 * P[0] + 2 * ((1 - u) * u)[:, None] * P[1] + (u ** 2)[:, None] * P[2]
    d = 2 * (1 - u)[:, None] * (P[1] - P[0]) + 2 * u[:, None] * (P[2] - P[1])
    nrm = np.c_[-d[:, 1], d[:, 0]] / np.maximum(np.linalg.norm(d, axis=1), 1e-300)[:, None]
    th = 5 * tmax_c * chord * (NACA[0] * np.sqrt(u) + NACA[1] * u + NACA[2] * u ** 2 + NACA[3] * u ** 3 + NACA[4] * u ** 4)
    th = np.maximum(th, T_TE_TURB * u)            # конечная толщина выходной кромки
    th[-1] = max(th[-1], T_TE_TURB)
    if le_ext:
        cam = cam - np.array([le_ext * math.cos(math.atan(t_in)), le_ext * math.sin(math.atan(t_in))]) * (1 - u)[:, None] ** 3
    up = cam + 0.5 * th[:, None] * nrm
    lo = cam - 0.5 * th[:, None] * nrm
    poly = np.vstack([up, lo[-2:0:-1]])
    return poly, cam, th, nrm, u


def cavity(cam, th, nrm, u, wall, core_min=1.0e-3, u0=0.12, u1=0.88):
    """полость с равной стенкой: внутренний контур — средняя линия ± (t − 2·wall)/2 там, где t − 2·wall ≥ core_min
    (иначе стенки сходятся — полости нет); возврат контура полости или None"""
    ti = th - 2 * wall
    m = (u >= u0) & (u <= u1) & (ti >= core_min)
    if m.sum() < 3:
        return None
    idx = np.nonzero(m)[0]
    # только одна непрерывная полость — самый длинный участок
    br = np.split(idx, np.nonzero(np.diff(idx) > 1)[0] + 1)
    idx = max(br, key=len)
    if len(idx) < 3:
        return None
    c, t_, n_ = cam[idx], ti[idx], nrm[idx]
    up = c + 0.5 * t_[:, None] * n_
    lo = c - 0.5 * t_[:, None] * n_
    return np.vstack([up, lo[::-1]])


def raw_moments(poly):
    """A, Sx = ∫y dA, Sy = ∫x dA, ∫x² dA, ∫y² dA, ∫xy dA (со знаком обхода)"""
    x, y = poly[:, 0], poly[:, 1]
    x1, y1 = np.roll(x, -1), np.roll(y, -1)
    cr = x * y1 - x1 * y
    A = 0.5 * cr.sum()
    Sy = ((x + x1) * cr).sum() / 6
    Sx = ((y + y1) * cr).sum() / 6
    Ixx0 = ((y * y + y * y1 + y1 * y1) * cr).sum() / 12
    Iyy0 = ((x * x + x * x1 + x1 * x1) * cr).sum() / 12
    Ixy0 = ((x * y1 + 2 * x * y + 2 * x1 * y1 + x1 * y) * cr).sum() / 24
    s = 1.0 if A >= 0 else -1.0
    return np.array([A, Sx, Sy, Ixx0, Iyy0, Ixy0]) * s


def props_from_raw(R):
    A, Sx, Sy, Ixx0, Iyy0, Ixy0 = R
    cx, cy = Sy / A, Sx / A
    return dict(A=A, cx=cx, cy=cy, Ixx=Ixx0 - A * cy * cy, Iyy=Iyy0 - A * cx * cx, Ixy=Ixy0 - A * cx * cy)


def section_props(poly):
    """площадь, центр тяжести, моменты инерции Ixx, Iyy, Ixy относительно центра (формулы Грина)"""
    x, y = poly[:, 0], poly[:, 1]
    x1, y1 = np.roll(x, -1), np.roll(y, -1)
    cr = x * y1 - x1 * y
    A = 0.5 * cr.sum()
    cx = ((x + x1) * cr).sum() / (6 * A)
    cy = ((y + y1) * cr).sum() / (6 * A)
    Ixx = ((y * y + y * y1 + y1 * y1) * cr).sum() / 12 - A * cy * cy          # относительно оси x (изгиб по y)
    Iyy = ((x * x + x * x1 + x1 * x1) * cr).sum() / 12 - A * cx * cx
    Ixy = ((x * y1 + 2 * x * y + 2 * x1 * y1 + x1 * y) * cr).sum() / 24 - A * cx * cy
    s = 1.0 if A > 0 else -1.0
    return dict(A=abs(A), cx=cx, cy=cy, Ixx=s * Ixx, Iyy=s * Iyy, Ixy=s * Ixy)


def turbine_rotor_blade(st, p, rpm, rho=7916.0, ns=11, route=None):
    """лопатка ротора ступени st (словарь ступени meanline_l1.ax_turbine) и гены p → геометрия, масса, напряжения
    в корне (центробежные + изгиб газом + изгиб от смещения центров тяжести), признаки"""
    route = ("print", "cast")[int(p["route"])] if route is None else route
    tmin = T_MIN[route]
    rm = st["Dm"] / 2
    Hr = 0.5 * (st["H2"] + st["H3"])
    rh, rt = rm - Hr / 2, rm + Hr / 2
    s_ = np.linspace(0, 1, ns)
    r = rh + s_ * Hr
    phi, psi, nv = st["phi"], st["psi"], p["vortex"]
    tan_b2 = (psi * (r / rm) ** nv - r / rm) / phi              # вход ротора (знак — как meanline: + со стороны набегания)
    tan_b3 = (r / rm) / phi
    b_in, b_out = np.arctan(tan_b2), np.arctan(tan_b3)
    c_root = st["chord_R"] * 1.0
    k_c = math.sqrt(max(p["taper"], 1e-3))
    chord = c_root * (1 - (1 - k_c) * s_)                      # A ∝ хорда² → A_tip/A_root = taper
    tub = p["tub_A"] * chord * np.sin(2 * np.pi * int(p["tub_n"]) * s_) if p["tub_A"] > 0 else np.zeros(ns)
    H = Hr
    e_th = H * (math.tan(math.radians(p["lean"])) * s_ + 4 * p["bow"] * s_ * (1 - s_))
    e_x = H * math.tan(math.radians(p["sweep"])) * s_
    hwall = p["hwall"] * 1e-3 if route == "print" else 0.0
    secs, props, cavs = [], [], []
    for j in range(ns):
        poly, cam, th, nrm, uu = section(chord[j], b_in[j], b_out[j], p["tmax_c"], le_ext=max(tub[j], 0.0))
        Ro = raw_moments(poly)
        pr = props_from_raw(Ro)
        cav = None
        if hwall > 0 and 0.12 <= s_[j] <= 0.92:            # корень и торец сплошные (оценка)
            cav = cavity(cam, th, nrm, uu, hwall)
        if cav is not None:
            pe = props_from_raw(Ro - raw_moments(cav))
        else:
            pe = pr
        pr["A_eff"], pr["Ixx_eff"], pr["Iyy_eff"], pr["Ixy_eff"] = pe["A"], pe["Ixx"], pe["Iyy"], pe["Ixy"]
        pr["cavity_frac"] = 1 - pe["A"] / pr["A"]
        pr["t_max"] = float(th.max())
        pr["wall_min"] = hwall if cav is not None else float(th.max())
        poly_c = poly - np.array([pr["cx"], pr["cy"]]) + np.array([e_x[j], e_th[j]])
        secs.append(poly_c)
        cavs.append(None if cav is None else cav - np.array([pr["cx"], pr["cy"]]) + np.array([e_x[j], e_th[j]]))
        props.append(pr)
    hollow = float(np.mean([q["cavity_frac"] for q in props]))
    A = np.array([q["A_eff"] for q in props])
    om = rpm * math.pi / 30
    # центробежная сила в корне и среднее напряжение
    Fc = rho * om * om * float(np.trapezoid(A * r, r))
    sig_c = Fc / props[0]["A"]
    # моменты в корне: газ (окружная сила на единицу высоты) и смещения центров (M = −F e r_root/r для окружного,
    # M = −F e для осевого смещения; вывод в README)
    m_dot = st.get("m_dot", 25.53)
    Ft = m_dot * psi * st["Dm"] * 0 + m_dot * (psi * (om * rm)) / st["zt"]      # Н на лопатку (Δcθ = ψU)
    ft = Ft / H
    Mx_gas = ft * H * H / 2
    dF = rho * om * om * A * r
    Mx_cb = -float(np.trapezoid(dF * e_th * (rh / r), r))
    My_cb = -float(np.trapezoid(dF * e_x, r))
    Mx, My = Mx_gas + Mx_cb, My_cb
    q = props[0]
    den = q["Ixx_eff"] * q["Iyy_eff"] - q["Ixy_eff"] ** 2
    P0 = secs[0] - np.array([e_x[0], e_th[0]])
    xv, yv = P0[:, 0], P0[:, 1]
    sig_b = ((My * q["Ixx_eff"] + Mx * q["Ixy_eff"]) * xv - (Mx * q["Iyy_eff"] + My * q["Ixy_eff"]) * yv) / den
    sig_b_gas = ((Mx_gas * q["Ixy_eff"]) * xv - (Mx_gas * q["Iyy_eff"]) * yv) / den
    mass = rho * float(np.trapezoid(A, r))
    J = rho * float(np.trapezoid(A * r * r, r))
    # проверки
    wall = min(q_["wall_min"] for q_ in props)
    tmin_geom = min(float(T_TE_TURB), wall)
    si = any(C.polyline_self_intersects(sp) for sp in secs[::2])
    pitch_h = 2 * np.pi * rh / st["zt"]
    throat = pitch_h * np.cos(b_out[0]) - props[0]["t_max"]
    checks = dict(wall_min_m=wall, wall_ok=wall >= tmin - 1e-12, self_intersection=si, throat_hub_m=float(throat),
                  throat_ok=bool(throat > 2 * tmin), te_m=T_TE_TURB)
    checks["valid"] = bool(checks["wall_ok"] and not si and checks["throat_ok"])
    twist = float(math.degrees(b_in[-1] - b_in[0]))
    feats = dict(twist_in_deg=twist, twist_out_deg=float(math.degrees(b_out[-1] - b_out[0])),
                 bow_endwall_deg=float(math.degrees(math.atan(4 * p["bow"]))), sweep_deg=p["sweep"], lean_deg=p["lean"])
    sc = C.silhouette_complexity(np.vstack([np.c_[e_x - 0.5 * chord, r], np.c_[(e_x + 0.5 * chord)[::-1], r[::-1]]]))
    feats["silhouette_iso"] = sc["isoperimetric"]
    return dict(r=r, rh=rh, rt=rt, H=H, chord=chord, secs=secs, cavs=cavs, props=props, A=A, mass_blade=mass, J_blade=J, Fc=Fc,
                sig_c_MPa=sig_c / 1e6, sig_b_peak_MPa=float(sig_b.max() / 1e6), sig_b_gas_MPa=float(np.abs(sig_b_gas).max() / 1e6),
                sig_peak_MPa=float((sig_c + sig_b.max()) / 1e6), Mx_gas=Mx_gas, Mx_cb=Mx_cb, My_cb=My_cb, checks=checks,
                features=feats, hollow=hollow, route=route, b_in=b_in, b_out=b_out, e_th=e_th, e_x=e_x, ft=ft, tmin_route=tmin,
                tmin_geom=tmin_geom)


def blade_mesh_from_sections(bl, zt, n_th_rep=None):
    """поверхность лопатки по сечениям (лофт ломаных одинаковой длины) и её Z копий; ось вала — x"""
    secs = bl["secs"]
    r = bl["r"]
    n = len(secs[0])
    V = []
    for sp, rv in zip(secs, r):
        x = sp[:, 0]
        th = sp[:, 1] / rv
        V.append(np.c_[x, rv * np.cos(th), rv * np.sin(th)])
    V = np.array(V)                                   # (ns, n, 3)
    ns = len(V)
    Vf = V.reshape(-1, 3)
    T = []
    for i in range(ns - 1):
        for j in range(n):
            a, b = i * n + j, i * n + (j + 1) % n
            c, d = a + n, b + n
            T.append((a, b, c))
            T.append((b, d, c))
    # торцы (веером от центра)
    cen0 = len(Vf)
    cen1 = cen0 + 1
    Vf = np.vstack([Vf, V[0].mean(0), V[-1].mean(0)])
    for j in range(n):
        T.append((cen0, (j + 1) % n, j))
        T.append((cen1, (ns - 1) * n + j, (ns - 1) * n + (j + 1) % n))
    T = np.array(T)
    reps = [(rotate_x(Vf, 2 * np.pi * k / zt), T) for k in range(zt)]
    return merge_meshes(reps)
