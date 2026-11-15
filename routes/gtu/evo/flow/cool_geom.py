#!/usr/bin/env python3
"""Геометрия «сосудистых» сетей охлаждения стенки: геном [0, 1]^16 → гладкая сеть каналов на развёртке стенки.

Развёртка — как в organic/thermal/vascular.py: x — от входной кромки к выходной (лопатка соплового аппарата, СА) или
вдоль оси камеры (жаровая труба), y — по высоте лопатки / по окружности; края по y теплоизолированы (у обоих случаев,
чтобы прямой расчёт vascular.py шёл в той же постановке). Подвод — узел 0 (полость лопатки или кольцевой зазор камеры),
корни сети соединяются с ним «коллекторными» сегментами без теплообмена; окончания — выпуски в газ при местном
статическом давлении (у СА у входной кромки давление газа выше — проверяется отсутствие заброса газа).

Архитектуры (ген arch) — пять «видов», каждый со своим источником:
  parallel     — параллельные каналы, по желанию змейкой и с наклоном (база сравнения; vascular: почти оптимальны у СА);
  dichotomous  — вилочное дерево «дельта» с асимметрией деления (Horsfield) и законом Мюррея d^k (Murray 1926);
  fishbone     — перистое жилкование листа: ствол и боковые жилки под углом, изогнутые к краю (жилкование двудольных);
  htree        — Н-дерево (фрактал, заполняющий плоскость; размерность → 2 с глубиной) от точечных подводов;
  reticulate   — дерево с поперечными перемычками (анастомозы): петли сетей листа устойчивы к повреждению и колебаниям
                 нагрузки (Katifori, Szöllősi, Magnasco 2010; Corson 2010 — по памяти, первоисточники не открывались).
Гладкость: каждая ветвь — кубическая кривая Эрмита (C¹ с ветвью-родителем), змейка — синус с огибающей sin², у которой
нулевые значение и производная на концах (C¹ в узлах сохраняется); «живая» нерегулярность — гладкое поле смещений
(низкие гармоники) с детерминированным зерном из генома. Только гладкие представления — урок organic/structure.

Проверки изготовимости: диаметр ≥ d_min способа, стенка между каналами ≥ t_min (выборка точек осей и дерево cKDTree;
у общих узлов — вне зоны стыка), каналы внутри развёртки, радиус кривизны оси ≥ R_MIN_D·d, способ «литьё» — без петель и
не глубже одного ветвления (керамические стержни), см. ROUTES. Числа способов — правила print_review (по выдержке
поиска) или оценка.
"""
from __future__ import annotations

import hashlib
import math
import os
import sys

import numpy as np
from scipy.spatial import cKDTree

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
from core import shapes as S  # noqa: E402

# --------------------------------------------------------------------------- случаи (размеры — как vascular.py)
CASES = {
    "ngv": dict(Lx=1.2 * 0.040, H=0.0281987, title="сторона лопатки СА 1-й ступени (развёртка 48 × 28,2 мм)",
                d_max=1.6e-3, src="organic/thermal/vascular.py: A['ngv'] (хорда 40 мм × 1,2 — оценка), "
                                  "frame3_point_results.json → turbine.options.3.stages[0].h_rotor_in_m"),
    "liner": dict(Lx=0.10, H=0.06, title="участок жаровой трубы у первичной зоны (развёртка 100 × 60 мм)",
                  d_max=2.5e-3, src="organic/thermal/vascular.py: A['liner'] (оценка)"),
}

ARCH = ("parallel", "dichotomous", "fishbone", "htree", "reticulate")
ARCH_RU = {"parallel": "параллельные/змейка", "dichotomous": "вилочное дерево", "fishbone": "перистое (жилкование листа)",
           "htree": "Н-дерево (фрактал)", "reticulate": "сеть с петлями (анастомозы)"}

# способы изготовления стенки с каналами; d_min, t_min, шероховатость — правила print_review или оценка
ROUTES = {
    "print": dict(label="печать LPBF, каналы как напечатано", d_min=0.5e-3, t_min=0.4e-3, eps=30e-6, parts=1, ops=1,
                  src="d_min 0,5 мм: поставщик очищал 0,4 мм, надёжно от 1 мм (print_review RULES depowder_ideal_min_mm, "
                      "по выдержке поиска); t_min 0,4 мм ≈ наименьший элемент IN718 0,381 мм (min_feature_in718_mm); "
                      "ε 30 мкм — print_review roughness_ngv (оценка)"),
    "print_afm": dict(label="печать LPBF + абразивная проливка каналов", d_min=0.5e-3, t_min=0.4e-3, eps=10e-6, parts=1, ops=2,
                      src="ε ≤ 15 мкм после проливки — vascular README п. 4 (оценка)"),
    "cast": dict(label="литьё по выплавляемым моделям с керамическими стержнями", d_min=1.0e-3, t_min=0.8e-3, eps=5e-6, parts=1,
                 ops=1, K_max=1, loops=False, src="стержни тоньше 1 мм и разветвлённые стержни хрупки — оценка"),
    "hybrid": dict(label="литой корпус + печатная вставка с каналами (пайка)", d_min=0.5e-3, t_min=0.4e-3, eps=30e-6, parts=2,
                   ops=2, src="оценка: вставка печатается, корпус льётся (КНР), одна паяная поверхность"),
}
ROUTE_KEYS = tuple(ROUTES)
R_MIN_D = 1.0          # радиус кривизны оси канала не меньше диаметра (оценка; потери в гибе и вынос порошка)

# --------------------------------------------------------------------------- геном
GENES = [
    ("arch", "архитектура: " + ", ".join(ARCH)),
    ("n_root", "число корней (подводов) 1…16, у параллельных — каналов 4…28"),
    ("depth", "глубина ветвления 1…4 (у Н-дерева 1…5)"),
    ("gamma", "где ветвиться по длине: x_k = x0 + (Lx − x0)·(k/K)^γ, γ 0,5…2"),
    ("asym", "асимметрия деления окончаний q 0,5…0,8 (Horsfield)"),
    ("k_murray", "показатель закона Мюррея 2,0…3,0 (3 — ламинарный, 7/3 — турбулентный)"),
    ("v_eq", "объём каналов на площадь стенки (эквивалентная толщина) 0,06…0,40 мм — плотность сети"),
    ("tilt", "наклон каналов к оси (закрутка) −40…+40°"),
    ("wig_amp", "амплитуда змейки, доля шага 0…0,35"),
    ("wig_lam", "длина волны змейки, доля Lx 0,15…1"),
    ("lat_ang", "угол боковых жилок / разлёт ветвей 25…75°"),
    ("loops", "доля перемычек (петли) 0…1"),
    ("throttle", "дроссель подвода: p_in между порогом обратного тока и наибольшим"),
    ("route", "способ: " + ", ".join(ROUTE_KEYS)),
    ("organic", "«живая» нерегулярность — амплитуда гладкого поля смещений 0…0,3 шага"),
    ("x_feed", "положение линии подвода 0…0,45 Lx (у Н-дерева — центр подводов)"),
]
DIM = len(GENES)
GI = {n: i for i, (n, _) in enumerate(GENES)}


def _lin(u, a, b):
    return a + (b - a) * float(np.clip(u, 0.0, 1.0))


def _pick(u, n):
    return int(min(n - 1, math.floor(float(np.clip(u, 0.0, 1.0)) * n)))


def decode(g):
    g = np.clip(np.asarray(g, float), 0.0, 1.0)
    arch = ARCH[_pick(g[GI["arch"]], len(ARCH))]
    p = dict(arch=arch,
             n_root=int(round(math.exp(_lin(g[GI["n_root"]], 0.0, math.log(16.0))))),
             n_par=int(round(_lin(g[GI["n_root"]], 4, 28))),
             depth=1 + _pick(g[GI["depth"]], 5 if arch == "htree" else 4),
             gamma=math.exp(_lin(g[GI["gamma"]], math.log(0.5), math.log(2.0))),
             asym=_lin(g[GI["asym"]], 0.5, 0.8),
             k_murray=_lin(g[GI["k_murray"]], 2.0, 3.0),
             v_eq=_lin(g[GI["v_eq"]], 0.06e-3, 0.40e-3),
             tilt_deg=_lin(g[GI["tilt"]], -40.0, 40.0),
             wig_amp=_lin(g[GI["wig_amp"]], 0.0, 0.35),
             wig_lam=_lin(g[GI["wig_lam"]], 0.15, 1.0),
             lat_ang=_lin(g[GI["lat_ang"]], 25.0, 75.0),
             loops=_lin(g[GI["loops"]], 0.0, 1.0),
             throttle=float(g[GI["throttle"]]),
             route=ROUTE_KEYS[_pick(g[GI["route"]], len(ROUTE_KEYS))],
             organic=_lin(g[GI["organic"]], 0.0, 0.3),
             x_feed=_lin(g[GI["x_feed"]], 0.0, 0.45))
    if arch != "reticulate":
        p["loops"] = 0.0
    p["seed"] = int.from_bytes(hashlib.sha256(np.round(g, 6).tobytes()).digest()[:4], "little")
    return p


# --------------------------------------------------------------------------- сеть
class Network:
    """узлы (x, y); рёбра {a, b, d, pts (n, 2), kind: feed|ch|link}; узел 0 — подвод; outlets — окончания (выпуски)"""

    def __init__(self, case, params):
        self.case, self.p = case, params
        c = CASES[case]
        self.Lx, self.H = c["Lx"], c["H"]
        self.nodes = [(0.0, 0.5 * self.H)]
        self.edges = []
        self.outlets = []
        self.roots = []
        self.n_term = {}

    def node(self, x, y):
        self.nodes.append((float(x), float(y)))
        return len(self.nodes) - 1

    def edge(self, a, b, pts, kind="ch", d=None):
        self.edges.append(dict(a=a, b=b, d=d, pts=np.asarray(pts, float), kind=kind))
        return len(self.edges) - 1

    @property
    def P(self):
        return np.asarray(self.nodes, float)

    def ch_edges(self):
        return [e for e in self.edges if e["kind"] != "feed"]

    def length(self):
        return float(sum(_plen(e["pts"]) for e in self.ch_edges()))

    def volume(self):
        return float(sum(math.pi * e["d"] ** 2 / 4 * _plen(e["pts"]) for e in self.ch_edges()))


def _plen(P):
    return float(np.sum(np.linalg.norm(np.diff(P, axis=0), axis=1))) if len(P) > 1 else 0.0


def hermite(p0, t0, p1, t1, n):
    """кубическая кривая Эрмита p0 → p1 с касательными t0, t1 (направления; модуль — длина хорды)"""
    p0, p1 = np.asarray(p0, float), np.asarray(p1, float)
    L = float(np.linalg.norm(p1 - p0))
    t0 = np.asarray(t0, float); t0 = t0 / (np.linalg.norm(t0) + 1e-300) * L
    t1 = np.asarray(t1, float); t1 = t1 / (np.linalg.norm(t1) + 1e-300) * L
    s = np.linspace(0.0, 1.0, max(3, int(n)))[:, None]
    h00, h10, h01, h11 = 2 * s ** 3 - 3 * s ** 2 + 1, s ** 3 - 2 * s ** 2 + s, -2 * s ** 3 + 3 * s ** 2, s ** 3 - s ** 2
    return h00 * p0 + h10 * t0 + h01 * p1 + h11 * t1


def wiggle(P, amp, lam, phase=0.0):
    """змейка: смещение по нормали amp·sin(2π s/λ + φ) с огибающей sin²(π s/L) (C¹ на концах сохраняется)"""
    if amp <= 0 or len(P) < 3:
        return P
    s = np.r_[0.0, np.cumsum(np.linalg.norm(np.diff(P, axis=0), axis=1))]
    L = s[-1]
    if L <= 0:
        return P
    T = np.gradient(P, axis=0)
    T /= np.linalg.norm(T, axis=1, keepdims=True) + 1e-300
    N = np.stack([-T[:, 1], T[:, 0]], 1)
    off = amp * np.sin(2 * math.pi * s / lam + phase) * np.sin(math.pi * s / L) ** 2
    return P + off[:, None] * N


def _npts(L, step):
    return max(5, int(math.ceil(L / step)) + 1)


# --------------------------------------------------------------------------- построители архитектур
def _shear(net, x, y, x0):
    """наклон (закрутка): y += (x − x0)·tg β"""
    return y + (x - x0) * math.tan(math.radians(net.p["tilt_deg"]))


def _y_band(net, x0):
    """полоса корней, при которой наклонённые каналы остаются внутри развёртки"""
    dy = (net.Lx - x0) * math.tan(math.radians(net.p["tilt_deg"]))
    lo, hi = max(0.0, -dy), net.H - max(0.0, dy)
    if hi - lo < 0.35 * net.H:                      # слишком крутой наклон для этой высоты — ограничиваем
        c = 0.5 * (lo + hi)
        lo, hi = c - 0.175 * net.H, c + 0.175 * net.H
    return lo, hi


def build_parallel(net, step):
    p = net.p
    N = p["n_par"]
    x0 = net.Lx * p["x_feed"]
    lo, hi = _y_band(net, x0)
    pitch = (hi - lo) / N
    tb = math.tan(math.radians(p["tilt_deg"]))
    for i in range(N):
        y0 = lo + pitch * (i + 0.5)
        r = net.node(x0, y0)
        net.roots.append(r)
        e = net.node(net.Lx, y0 + (net.Lx - x0) * tb)
        P = np.stack([np.linspace(x0, net.Lx, _npts(net.Lx - x0, step)), np.linspace(y0, y0 + (net.Lx - x0) * tb, _npts(net.Lx - x0, step))], 1)
        P = wiggle(P, p["wig_amp"] * pitch, p["wig_lam"] * net.Lx, phase=0.0)
        net.edge(r, e, P)
        net.outlets.append(e)
        net.n_term[len(net.edges) - 1] = 1


def _tree(net, root_xy, span, K, x_levels, step, tangent0=(1.0, 0.0), parent=None):
    """вилочное дерево от корня: уровень k на x_levels[k]; окончания делятся в доле q (асимметрия)"""
    p = net.p
    q = p["asym"]
    leaves_total = 2 ** K
    # (узел, y-центр полосы, ширина полосы, касательная, число окончаний в поддереве)
    level = [(parent if parent is not None else net.node(*root_xy), root_xy[1], span, np.asarray(tangent0, float), leaves_total)]
    if parent is None:
        net.roots.append(level[0][0])
    tb = math.tan(math.radians(p["tilt_deg"]))
    spread = math.radians(p["lat_ang"]) / math.radians(50.0)
    for k in range(1, K + 1):
        x = x_levels[k]
        new = []
        for (nid, yc, sp_, tan_, nl) in level:
            x_par, y_par = net.nodes[nid]
            n1 = max(1, int(round(nl * q))) if nl > 1 else 1
            n2 = max(nl - n1, 1)
            f1 = n1 / (n1 + n2)
            for sgn, frac, nn in ((-1, f1, n1), (1, 1.0 - f1, n2)):
                sub = sp_ * frac
                yc2 = (yc - sp_ / 2 + sub / 2) if sgn < 0 else (yc + sp_ / 2 - sub / 2)
                y_t = yc + (yc2 - yc) * min(1.0, spread) + (x - x_par) * tb
                cid = net.node(x, y_t)
                d_end = np.array([1.0, tb])
                L = math.hypot(x - x_par, y_t - y_par)
                P = hermite((x_par, y_par), tan_, (x, y_t), d_end, _npts(L, step))
                P = wiggle(P, p["wig_amp"] * 0.5 * sub, p["wig_lam"] * net.Lx, phase=0.7 * cid)
                net.edge(nid, cid, P)
                net.n_term[len(net.edges) - 1] = nn
                new.append((cid, yc2 + (x - x_par) * tb, sub, d_end, nn))
        level = new
    for (nid, _, _, _, _) in level:
        net.outlets.append(nid)
    return [nid for nid, *_ in level]


def build_dichotomous(net, step):
    p = net.p
    K = p["depth"]
    nr = max(1, min(p["n_root"], 16))
    x0 = net.Lx * p["x_feed"]
    lo, hi = _y_band(net, x0)
    span = (hi - lo) / nr
    xl = [x0 + (net.Lx - x0) * (k / K) ** p["gamma"] for k in range(K + 1)]
    for r in range(nr):
        _tree(net, (x0, lo + span * (r + 0.5)), span, K, xl, step)


def build_fishbone(net, step):
    """стволы вдоль x (с наклоном), боковые жилки поочерёдно под углом φ, изогнутые к выходной кромке; глубина 2 —
    вторичные жилки на боковых"""
    p = net.p
    nr = max(1, min(p["n_root"], 8))
    x0 = net.Lx * p["x_feed"]
    lo, hi = _y_band(net, x0)
    span = (hi - lo) / nr
    tb = math.tan(math.radians(p["tilt_deg"]))
    phi = math.radians(p["lat_ang"])
    n_lat = max(2, int(round(2 + 6 * (p["gamma"] - 0.5) / 1.5)))       # число боковых жилок на сторону 2…8
    for r in range(nr):
        yc = lo + span * (r + 0.5)
        root = net.node(x0, yc)
        net.roots.append(root)
        xs = x0 + (net.Lx - x0) * (np.arange(1, 2 * n_lat + 1) / (2 * n_lat + 1)) ** (1.0 / p["gamma"] ** 0.5)
        prev, prev_t = root, np.array([1.0, tb])
        trunk_edges = []
        for j, xj in enumerate(xs):
            yj = yc + (xj - x0) * tb
            nid = net.node(xj, yj)
            L = math.hypot(xj - net.nodes[prev][0], yj - net.nodes[prev][1])
            P = hermite(net.nodes[prev], prev_t, (xj, yj), (1.0, tb), _npts(L, step))
            net.edge(prev, nid, P)
            trunk_edges.append(len(net.edges) - 1)
            sgn = 1 if j % 2 == 0 else -1
            reach = 0.45 * span
            Ll = reach / max(math.sin(phi), 0.3)
            ang = math.atan(tb) + sgn * phi
            xe = min(net.Lx, xj + Ll * math.cos(ang) * 0.8)
            ye = yj + sgn * reach
            leaf = net.node(xe, ye)
            t_end = np.array([1.0, 0.25 * sgn])               # кончик жилки поворачивает к выходной кромке
            P = hermite((xj, yj), (math.cos(ang), math.sin(ang)), (xe, ye), t_end, _npts(math.hypot(xe - xj, ye - yj), step))
            P = wiggle(P, p["wig_amp"] * 0.25 * span, p["wig_lam"] * net.Lx, phase=0.3 * leaf)
            net.edge(nid, leaf, P)
            n_lat_leaves = 1
            if p["depth"] >= 2:                                # вторичные жилки (порядок 3)
                mid = P[len(P) // 2]
                tmid = P[len(P) // 2 + 1] - P[len(P) // 2 - 1]
                # делим боковую на две части в середине
                net.edges.pop()
                mnode = net.node(*mid)
                P1, P2 = P[: len(P) // 2 + 1], P[len(P) // 2:]
                net.edge(nid, mnode, P1); e1 = len(net.edges) - 1
                net.edge(mnode, leaf, P2); e2 = len(net.edges) - 1
                ang2 = math.atan2(tmid[1], tmid[0]) - sgn * phi * 0.8
                L2 = 0.4 * Ll
                xe2 = min(net.Lx, mid[0] + L2 * math.cos(ang2)); ye2 = mid[1] + L2 * math.sin(ang2)
                leaf2 = net.node(xe2, ye2)
                P3 = hermite(mid, (math.cos(ang2), math.sin(ang2)), (xe2, ye2), (1.0, 0.0), _npts(L2, step))
                net.edge(mnode, leaf2, P3); e3 = len(net.edges) - 1
                net.n_term[e1] = 2; net.n_term[e2] = 1; net.n_term[e3] = 1
                net.outlets += [leaf, leaf2]
                n_lat_leaves = 2
            else:
                net.n_term[len(net.edges) - 1] = 1
                net.outlets.append(leaf)
            net.nodes_lat = getattr(net, "nodes_lat", {})
            net.nodes_lat[nid] = n_lat_leaves
            prev, prev_t = nid, np.array([1.0, tb])
        # хвост ствола до выходной кромки — тоже выпуск
        xe, ye = net.Lx, yc + (net.Lx - x0) * tb
        tail = net.node(xe, ye)
        L = math.hypot(xe - net.nodes[prev][0], ye - net.nodes[prev][1])
        net.edge(prev, tail, hermite(net.nodes[prev], prev_t, (xe, ye), (1.0, tb), _npts(L, step)))
        net.n_term[len(net.edges) - 1] = 1
        net.outlets.append(tail)
        # число окончаний ниже по стволу
        below = 1
        for j in range(len(trunk_edges) - 1, -1, -1):
            nid = net.edges[trunk_edges[j]]["b"]
            below += net.nodes_lat.get(nid, 1)
            net.n_term[trunk_edges[j]] = below


def build_htree(net, step):
    """Н-дерево: от центра — отрезок, на концах — перпендикулярные отрезки длиной ×1/√2; углы скруглены (кривые Эрмита
    с касательной родителя). Подводы — на решётке m × n, каждый в своей клетке развёртки."""
    p = net.p
    K = p["depth"]
    nr = max(1, min(p["n_root"], 12))
    asp = net.Lx / net.H
    m = max(1, int(round(math.sqrt(nr * asp))))
    n = max(1, int(round(nr / m)))
    cw, ch = net.Lx / m, net.H / n
    tb = math.tan(math.radians(p["tilt_deg"]))
    for i in range(m):
        for j in range(n):
            cx, cy = (i + 0.5) * cw, (j + 0.5) * ch
            root = net.node(cx, cy)
            net.roots.append(root)
            # первый отрезок — вдоль большей стороны клетки
            horiz = cw >= ch
            L0 = 0.5 * (cw if horiz else ch) * 0.92
            stack = [(root, np.array([1.0, tb]) if horiz else np.array([0.0, 1.0]), L0, 0)]
            tips = []
            while stack:
                nid, direc, L, lev = stack.pop()
                x, y = net.nodes[nid]
                if lev >= K:
                    tips.append(nid)
                    continue
                dvec = direc / np.linalg.norm(direc)
                for sgn in (-1.0, 1.0):
                    e = np.array([x, y]) + sgn * L * dvec
                    e[0] = min(max(e[0], 0.02 * net.Lx), 0.98 * net.Lx)
                    e[1] = min(max(e[1], 0.02 * net.H), 0.98 * net.H)
                    cid = net.node(*e)
                    perp = np.array([-dvec[1], dvec[0]])
                    # ветвь: из узла вдоль ±d, конец — касательная поворачивает к перпендикуляру (скругление)
                    P = hermite((x, y), sgn * dvec, e, sgn * dvec + 0.35 * perp, _npts(L, step))
                    P = wiggle(P, p["wig_amp"] * 0.15 * L, p["wig_lam"] * net.Lx, phase=0.5 * cid)
                    net.edge(nid, cid, P)
                    net.n_term[len(net.edges) - 1] = 2 ** (K - lev - 1)
                    stack.append((cid, perp, L / math.sqrt(2.0), lev + 1))
            net.outlets += tips


def add_loops(net, frac, step):
    """перемычки между соседними окончаниями дерева (у середин концевых рёбер): петли, как в жилковании листа"""
    if frac <= 0:
        return 0
    ends = []
    for ei, e in enumerate(net.edges):
        if e["kind"] == "ch" and e["b"] in net.outlets:
            mid = e["pts"][len(e["pts"]) // 2]
            ends.append((ei, mid))
    if len(ends) < 2:
        return 0
    M = np.array([m for _, m in ends])
    o = np.lexsort((M[:, 0], M[:, 1]))                  # по высоте — соседи по y
    pairs = [(o[i], o[i + 1]) for i in range(len(o) - 1)]
    k = int(round(frac * len(pairs)))
    rng = np.random.default_rng(net.p["seed"])
    chosen = sorted(rng.choice(len(pairs), size=k, replace=False).tolist()) if k else []
    n_add = 0
    split_at = {}
    for c in chosen:
        i, j = pairs[c]
        for idx in (i, j):
            ei = ends[idx][0]
            if ei not in split_at:
                split_at[ei] = None
    # делим выбранные концевые рёбра в середине (новый узел), затем соединяем пары
    for ei in sorted(split_at):
        e = net.edges[ei]
        P = e["pts"]
        h = len(P) // 2
        mnode = net.node(*P[h])
        nt = net.n_term.get(ei, 1)
        net.edges[ei] = dict(a=e["a"], b=mnode, d=None, pts=P[: h + 1], kind="ch")
        net.edge(mnode, e["b"], P[h:])
        net.n_term[ei] = nt
        net.n_term[len(net.edges) - 1] = nt
        split_at[ei] = mnode
    for c in chosen:
        i, j = pairs[c]
        a, b = split_at[ends[i][0]], split_at[ends[j][0]]
        pa, pb = np.array(net.nodes[a]), np.array(net.nodes[b])
        L = float(np.linalg.norm(pb - pa))
        if L < 1e-6:
            continue
        dirv = (pb - pa) / L
        P = hermite(pa, dirv + np.array([0.25, 0.0]), pb, dirv - np.array([0.25, 0.0]), _npts(L, step))
        net.edge(a, b, P, kind="link")
        n_add += 1
    return n_add


def organic_field(net):
    """гладкое поле смещений (гармоники 1–2 по обеим осям), амплитуда — доля характерного шага"""
    a = net.p["organic"]
    if a <= 0:
        return
    rng = np.random.default_rng(net.p["seed"] + 17)
    C = rng.standard_normal((2, 3, 3))
    pitch = net.H / max(4, len(net.outlets) ** 0.5 * 2)
    amp = a * pitch

    def disp(P):
        x = P[:, 0] / net.Lx; y = P[:, 1] / net.H
        out = np.zeros_like(P)
        for c in range(2):
            for i in range(3):
                for j in range(3):
                    out[:, c] += C[c, i, j] * np.sin(math.pi * (i + 1) * x + 0.3 * j) * np.sin(math.pi * (j + 1) * y + 0.2 * i)
        return amp * out / 3.0
    keep = set([0] + list(net.roots))
    Pn = net.P
    Dn = disp(Pn)
    for k in range(len(Pn)):
        if k in keep:
            Dn[k] = 0.0
    for e in net.edges:
        P = e["pts"]
        s = np.linspace(0.0, 1.0, len(P))[:, None]
        Dp = disp(P)
        # на концах — смещения узлов (C⁰), внутри — поле; огибающая убирает разрыв
        env = np.sin(math.pi * s) ** 2
        e["pts"] = P + env * Dp + (1 - s) * Dn[e["a"]] + s * Dn[e["b"]]
    net.nodes = [tuple(Pn[k] + Dn[k]) for k in range(len(Pn))]


def assign_diameters(net):
    """закон Мюррея от окончаний: d ∝ n_term^(1/k); перемычки — 0,7 d окончания; масштаб — под объём v_eq·Lx·H"""
    k = net.p["k_murray"]
    for ei, e in enumerate(net.edges):
        if e["kind"] == "ch":
            e["d"] = float(net.n_term.get(ei, 1)) ** (1.0 / k)
        elif e["kind"] == "link":
            e["d"] = 0.7
    V_target = net.p["v_eq"] * net.Lx * net.H
    V = sum(math.pi * e["d"] ** 2 / 4 * _plen(e["pts"]) for e in net.edges if e["kind"] != "feed")
    s = math.sqrt(V_target / max(V, 1e-30))
    for e in net.edges:
        if e["kind"] != "feed":
            e["d"] *= s
    # коллекторные сегменты от подвода к корням — без сопротивления и теплообмена
    for r in net.roots:
        P = np.array([net.nodes[0], net.nodes[r]])
        net.edges.insert(0, dict(a=0, b=r, d=10e-3, pts=P, kind="feed"))
    # n_term после вставки сдвинулся — пересобираем по рёбрам
    return s


def build(case, g_or_params, step=None):
    p = decode(g_or_params) if not isinstance(g_or_params, dict) else g_or_params
    net = Network(case, p)
    step = step or 0.02 * min(net.Lx, net.H) * 2.5
    {"parallel": build_parallel, "dichotomous": build_dichotomous, "fishbone": build_fishbone,
     "htree": build_htree, "reticulate": build_dichotomous}[p["arch"]](net, step)
    if p["arch"] == "reticulate":
        add_loops(net, p["loops"], step)
    organic_field(net)
    nt = dict(net.n_term)
    assign_diameters(net)
    net.n_term = {i + len(net.roots): v for i, v in nt.items()}
    return net


# --------------------------------------------------------------------------- проверки
def checks(net, route=None):
    route = route or net.p["route"]
    R = ROUTES[route]
    out = {}
    chs = [(i, e) for i, e in enumerate(net.edges) if e["kind"] != "feed"]
    ds = np.array([e["d"] for _, e in chs])
    out["d_min_mm"] = float(ds.min() * 1e3)
    out["d_max_mm"] = float(ds.max() * 1e3)
    out["thin_lt_1mm_frac"] = float(sum(_plen(e["pts"]) for _, e in chs if e["d"] < 1e-3) / max(net.length(), 1e-12))
    # выборка точек осей
    pts, eid, dd, sarc = [], [], [], []
    for i, e in chs:
        P = e["pts"]
        L = _plen(P)
        n = max(3, int(math.ceil(L / (0.4 * max(e["d"], 0.3e-3)))))
        s = np.r_[0.0, np.cumsum(np.linalg.norm(np.diff(P, axis=0), axis=1))]
        t = np.linspace(0.0, s[-1], n)
        Q = np.stack([np.interp(t, s, P[:, 0]), np.interp(t, s, P[:, 1])], 1)
        pts.append(Q); eid.append(np.full(n, i)); dd.append(np.full(n, e["d"])); sarc.append(t)
    Q = np.vstack(pts); E = np.concatenate(eid); D = np.concatenate(dd); A = np.concatenate(sarc)
    tree = cKDTree(Q)
    rmax = float(D.max() + R["t_min"]) * 1.05
    pr = tree.query_pairs(rmax, output_type="ndarray")
    margin = np.inf
    if len(pr):
        i, j = pr[:, 0], pr[:, 1]
        dist = np.linalg.norm(Q[i] - Q[j], axis=1)
        gap = dist - 0.5 * (D[i] + D[j])
        same = E[i] == E[j]
        # одно ребро: соседние по дуге точки — не стенка
        far_same = same & (np.abs(A[i] - A[j]) > math.pi * np.maximum(D[i], D[j]) + rmax)
        # смежные рёбра (общий узел): исключаем зону стыка
        ea = np.array([net.edges[k]["a"] for k in range(len(net.edges))])
        eb = np.array([net.edges[k]["b"] for k in range(len(net.edges))])
        Pn = net.P
        diff = ~same
        share = np.zeros(len(i), bool)
        near_j = np.zeros(len(i), bool)
        for na, nb in ((ea, ea), (ea, eb), (eb, ea), (eb, eb)):
            s_ = diff & (na[E[i]] == nb[E[j]])
            if s_.any():
                node = na[E[i]]
                share |= s_
                dn_i = np.linalg.norm(Q[i] - Pn[node], axis=1)
                dn_j = np.linalg.norm(Q[j] - Pn[node], axis=1)
                zone = 1.5 * (D[i] + D[j]) + 2 * R["t_min"]
                near_j |= s_ & ((dn_i < zone) | (dn_j < zone))
        # у перемычек и деления концевых рёбер — двухуровневое соседство (через один узел) тоже исключаем у стыка
        consider = far_same | (diff & ~near_j)
        if consider.any():
            margin = float(np.min(gap[consider]))
    out["wall_min_mm"] = float(margin * 1e3) if np.isfinite(margin) else 99.0
    # внутри развёртки
    yo = np.maximum(0.0, np.maximum(-(Q[:, 1] - 0.5 * D), (Q[:, 1] + 0.5 * D) - net.H)).max()
    xo = np.maximum(0.0, np.maximum(-Q[:, 0], Q[:, 0] - net.Lx)).max()
    out["outside_mm"] = float(max(xo, yo) * 1e3)
    # кривизна осей
    rmin = np.inf
    for _, e in chs:
        P = e["pts"]
        if len(P) >= 5 and _plen(P) > 2 * e["d"]:
            k = np.abs(S.curvature_2d(P))[2:-2]
            if len(k):
                rmin = min(rmin, 1.0 / max(float(k.max()), 1e-12) / e["d"])
    out["R_min_over_d"] = float(rmin) if np.isfinite(rmin) else 99.0
    out["self_intersections"] = int(sum(S.polyline_self_intersections(e["pts"]) for _, e in chs))
    # нарушения (доли) — для правил допустимости
    viol = {}
    viol["d_min"] = max(0.0, (R["d_min"] - ds.min()) / R["d_min"])
    viol["d_max"] = max(0.0, (ds.max() - CASES[net.case]["d_max"]) / CASES[net.case]["d_max"])
    viol["wall"] = max(0.0, (R["t_min"] - margin) / R["t_min"]) if np.isfinite(margin) else 0.0
    viol["outside"] = max(0.0, max(xo, yo) / R["t_min"])
    viol["curv"] = max(0.0, (R_MIN_D - rmin) / R_MIN_D) if np.isfinite(rmin) else 0.0
    if R.get("loops") is False and any(e["kind"] == "link" for e in net.edges):
        viol["cast_loops"] = 1.0
    if R.get("K_max") is not None and strahler(net) - 1 > R["K_max"]:
        viol["cast_depth"] = float(strahler(net) - 1 - R["K_max"])
    out["viol"] = {k: float(v) for k, v in viol.items() if v > 0}
    out["cv"] = float(sum(viol.values()))
    return out


# --------------------------------------------------------------------------- дескрипторы
def strahler(net):
    """порядок Хортона — Стралера по дереву от подвода (перемычки не учитываются)"""
    kids = {}
    for e in net.edges:
        if e["kind"] == "link":
            continue
        kids.setdefault(e["a"], []).append(e["b"])
    memo = {}

    def order(n, depth=0):
        if n in memo:
            return memo[n]
        ch = kids.get(n, [])
        if not ch or depth > 200:
            memo[n] = 1
            return 1
        os_ = [order(c, depth + 1) for c in ch]
        m = max(os_)
        r = m + 1 if os_.count(m) >= 2 else m
        memo[n] = r
        return r
    o = max((order(r) for r in net.roots), default=1)
    return int(o)


def raster(net, px=192):
    ny = max(8, int(round(px * net.H / net.Lx)))
    M = np.zeros((px, ny), bool)
    for e in net.ch_edges():
        P = e["pts"]
        L = _plen(P)
        n = max(2, int(math.ceil(L / (0.4 * net.Lx / px))))
        s = np.r_[0.0, np.cumsum(np.linalg.norm(np.diff(P, axis=0), axis=1))]
        t = np.linspace(0.0, s[-1], n)
        x = np.interp(t, s, P[:, 0]); y = np.interp(t, s, P[:, 1])
        i = np.clip((x / net.Lx * px).astype(int), 0, px - 1)
        j = np.clip((y / net.H * ny).astype(int), 0, ny - 1)
        M[i, j] = True
    return M


def box_dimension_lines(M, sizes=(1, 2, 4, 8, 16, 32)):
    """коробочная размерность множества пикселей осей (наклон log N(s) от log 1/s)"""
    ns, ss = [], []
    for s in sizes:
        if s > min(M.shape) // 2:
            break
        h, w = (M.shape[0] // s) * s, (M.shape[1] // s) * s
        b = M[:h, :w].reshape(h // s, s, w // s, s).any(axis=(1, 3))
        if b.sum() > 0:
            ns.append(b.sum()); ss.append(s)
    if len(ns) < 3:
        return 1.0
    p = np.polyfit(np.log(1.0 / np.asarray(ss, float)), np.log(np.asarray(ns, float)), 1)
    return float(p[0])


def descriptors(net):
    chs = net.ch_edges()
    L = np.array([_plen(e["pts"]) for e in chs])
    d = np.array([e["d"] for e in chs])
    dens = float(np.sum(d * L) / (net.Lx * net.H))
    ang, w = [], []
    for e in chs:
        v = np.diff(e["pts"], axis=0)
        ln = np.linalg.norm(v, axis=1)
        ang.append(np.arctan2(v[:, 1], np.abs(v[:, 0]) + 1e-300)); w.append(ln)
    ang = np.concatenate(ang); w = np.concatenate(w)
    tilt = float(abs(np.degrees(np.sum(ang * w) / max(w.sum(), 1e-300))))
    return dict(order=strahler(net), fractal=box_dimension_lines(raster(net)), density=dens, swirl_deg=tilt,
                n_outlets=len(net.outlets), n_edges=len(chs), length_mm=float(L.sum() * 1e3),
                volume_mm3=float(np.sum(math.pi * d ** 2 / 4 * L) * 1e9),
                n_links=int(sum(1 for e in net.edges if e["kind"] == "link")))


# --------------------------------------------------------------------------- 3D для рендера
def mesh_3d(net, kind=None, n_rep=None):
    """сетки для эскиза: стенка (лист) и каналы (трубки) на холодной стороне.
    ngv — развёртка на дуге (кривизна спинки лопатки, оценка радиуса 2,5 хорды); liner — сектор кольца жаровой трубы
    (развёртка повторена по окружности). Ось двигателя — z, вверх — y."""
    from core import mesh as M
    kind = kind or net.case
    meshes = []
    if kind == "liner":
        R = 0.30                                      # радиус жаровой трубы для эскиза (оценка)
        n_rep = n_rep or 6
        dth = net.H / R

        def to3(x, y, r=R, rep=0):
            th = (y / net.H + rep) * dth - 0.5 * n_rep * dth + math.pi / 2
            return np.stack([r * np.cos(th), r * np.sin(th), x], -1)
        z = np.linspace(0, net.Lx, 24)
        th = np.linspace(0, n_rep * dth, 24 * n_rep)
        Z, TH = np.meshgrid(z, th, indexing="ij")
        for r_, gray in ((R - 0.6e-3, 0.86),):
            th3 = TH - 0.5 * n_rep * dth + math.pi / 2
            V = np.stack([r_ * np.cos(th3), r_ * np.sin(th3), Z], -1).reshape(-1, 3)
            F = _grid_faces(len(z), len(th))
            meshes.append((V, F, gray))
        tubes = []
        for rep in range(n_rep):
            for e in net.ch_edges():
                P3 = to3(e["pts"][:, 0], e["pts"][:, 1], R + 0.5 * e["d"], rep)
                tubes.append(M.tube(P3, 0.5 * e["d"] * 1.0, n_sides=8))
        V, F = M.merge(tubes)
        meshes.append((V, F, 0.42))
        return meshes
    # ngv: дуга радиуса Rc в плоскости (x — по хорде, y — высота → z)
    Rc = 2.5 * 0.040

    def to3n(x, y, off=0.0):
        a = (x - 0.5 * net.Lx) / Rc
        return np.stack([(Rc + off) * np.sin(a), (Rc + off) * np.cos(a) - Rc, y], -1)
    x = np.linspace(0, net.Lx, 30); y = np.linspace(0, net.H, 20)
    X, Y = np.meshgrid(x, y, indexing="ij")
    V = to3n(X, Y, -0.5e-3).reshape(-1, 3)
    meshes.append((V, _grid_faces(len(x), len(y)), 0.86))
    tubes = [M.tube(to3n(e["pts"][:, 0], e["pts"][:, 1], 0.5 * e["d"]), 0.5 * e["d"], n_sides=8) for e in net.ch_edges()]
    V, F = M.merge(tubes)
    meshes.append((V, F, 0.42))
    # повернуть: ось двигателя z → по x развёртки (вид сбоку покажет развёртку)
    out = []
    for V, F, g in meshes:
        V2 = np.stack([V[:, 1], V[:, 2], V[:, 0]], 1)
        out.append((V2, F, g))
    return out


def _grid_faces(n0, n1):
    I, J = np.meshgrid(np.arange(n0 - 1), np.arange(n1 - 1), indexing="ij")
    a = (I * n1 + J).ravel(); b = (I * n1 + J + 1).ravel(); c = ((I + 1) * n1 + J + 1).ravel(); d = ((I + 1) * n1 + J).ravel()
    return np.vstack([np.stack([a, b, c], 1), np.stack([a, c, d], 1)])


if __name__ == "__main__":
    rng = np.random.default_rng(0)
    for arch_u in (0.1, 0.3, 0.5, 0.7, 0.9):
        g = rng.random(DIM); g[GI["arch"]] = arch_u
        net = build("ngv", g)
        print(net.p["arch"], descriptors(net), checks(net))
