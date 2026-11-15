"""Сплайновая лопатка (стойка, лопатка соплового аппарата, лопатка колеса) с законами по высоте и биомиметикой.

Сечение: средняя линия — B-сплайн прогибов по хорде, толщина — распределение NACA 4 (замкнутая кромка), откладывается
по нормали к средней линии; точки по хорде — косинусное сгущение к кромкам.
Законы по высоте (B-сплайны от относительной высоты s ∈ [0, 1]): хорда c(s), закрутка β(s), относительная толщина
t(s), линия укладки — стреловидность (sweep, вдоль хорды), навал (lean, поперёк хорды) и изгиб (bow — квадратичная
составляющая навала).
Биомиметика (статус — sources.py): бугорки передней кромки (tubercles, плавник горбатого кита: синусоидальная
модуляция положения передней кромки, амплитуда A·c, длина волны λ·c), зубчатая задняя кромка (serrations, треугольная
волна со скруглением), риблеты (продольные гребни высоты h и шага s по поверхности).
Ось высоты — z, хорда — x, толщина — y. Возврат — сетка точек (n_span, n_loop, 3) и треугольники.
"""
import math

import numpy as np

from . import shapes as S


def naca4_half(x, t):
    """полутолщина NACA 4 с замкнутой задней кромкой; x ∈ [0, 1], t — относительная толщина"""
    x = np.clip(x, 0.0, 1.0)
    return 5.0 * t * (0.2969 * np.sqrt(x) - 0.1260 * x - 0.3516 * x ** 2 + 0.2843 * x ** 3 - 0.1036 * x ** 4)


def section(chord, camber_ctrl, t_rel, n=40):
    """замкнутое сечение (2n − 2, 2) против часовой: от задней кромки по спинке к передней и по корытцу обратно"""
    beta = np.linspace(0.0, math.pi, n)
    x = 0.5 * (1.0 - np.cos(beta))
    cc = np.r_[0.0, np.asarray(camber_ctrl, float), 0.0]
    yc = S.bspline_eval(cc[:, None], x, p=min(3, len(cc) - 1))[:, 0]
    dyc = np.gradient(yc, x)
    th = np.arctan(dyc)
    yt = naca4_half(x, t_rel)
    xu, yu = x - yt * np.sin(th), yc + yt * np.cos(th)
    xl, yl = x + yt * np.sin(th), yc - yt * np.cos(th)
    P = np.vstack([np.stack([xu[::-1], yu[::-1]], 1), np.stack([xl[1:-1], yl[1:-1]], 1)])
    return P * chord


def _law(ctrl, s):
    c = np.asarray(ctrl, float)
    if c.ndim == 0 or len(c) == 1:
        return np.full_like(s, float(np.atleast_1d(c)[0]))
    return S.bspline_eval(c[:, None], s, p=min(3, len(c) - 1))[:, 0]


def blade(r_hub, r_tip, chord, twist_deg, t_rel, camber=(0.04, 0.05, 0.03), sweep=(0.0,), lean=(0.0,), bow=0.0,
          tubercle=None, serration=None, riblet=None, n_span=41, n_chord=40, x_stack=0.3):
    """chord, twist_deg, t_rel, sweep, lean — контрольные точки законов по высоте (м, градусы, доли, м, м);
    bow — прогиб линии навала в середине высоты, м; tubercle = (A_rel, lam_rel); serration = (A_rel, lam_rel);
    riblet = (h_rel, s_rel) — от хорды. Возврат: (P (n_span, n_loop, 3), F, laws)"""
    s = np.linspace(0.0, 1.0, n_span)
    z = r_hub + s * (r_tip - r_hub)
    c = _law(chord, s); tw = np.radians(_law(twist_deg, s)); tr = _law(t_rel, s)
    xs = _law(sweep, s); ys = _law(lean, s) + bow * 4.0 * s * (1.0 - s)
    loops = []
    for k in range(n_span):
        P = section(c[k], camber, tr[k], n_chord)
        xr = P[:, 0] / c[k]
        if tubercle:
            A, lam = tubercle
            dx = A * c[k] * math.sin(2 * math.pi * (z[k] - r_hub) / (lam * c[k] + 1e-12))
            P[:, 0] -= dx * (1.0 - S.smoothstep5(xr / 0.35))           # двигаем только носок, плавно
        if serration:
            A, lam = serration
            ph = ((z[k] - r_hub) / (lam * c[k] + 1e-12)) % 1.0
            tri = 1.0 - 4.0 * abs(ph - 0.5)                              # треугольная волна −1…1
            P[:, 0] += A * c[k] * tri * S.smoothstep5((xr - 0.75) / 0.25)
        if riblet:
            h, sp = riblet
            bump = np.maximum(0.0, np.cos(2 * math.pi * (z[k] - r_hub) / (sp * c[k] + 1e-12))) ** 4
            nrm = np.gradient(P, axis=0)[:, ::-1] * np.array([1.0, -1.0])
            nrm /= np.linalg.norm(nrm, axis=1, keepdims=True) + 1e-300
            P = P + h * c[k] * bump * nrm
        P = P - np.array([x_stack * c[k], 0.0])
        ca, sa = math.cos(tw[k]), math.sin(tw[k])
        Pr = np.stack([ca * P[:, 0] - sa * P[:, 1], sa * P[:, 0] + ca * P[:, 1]], 1) + np.array([xs[k], ys[k]])
        loops.append(np.column_stack([Pr, np.full(len(Pr), z[k])]))
    G = np.stack(loops)
    ns, nl = G.shape[:2]
    F = []
    for i in range(ns - 1):
        for j in range(nl):
            a, b = i * nl + j, i * nl + (j + 1) % nl
            F += [(a, b, b + nl), (a, b + nl, a + nl)]
    return G, np.asarray(F, np.int64), dict(s=s, z=z, chord=c, twist_deg=np.degrees(tw), t_rel=tr, sweep=xs, lean=ys)


def check_blade(G, laws):
    """сечения без самопересечений; измеренная закрутка (угол хорды от передней до задней кромки) против закона"""
    n_x = 0
    tw = []
    for k in range(G.shape[0]):
        P = G[k, :, :2]
        n_x += S.polyline_self_intersections(P, closed=True)
        te = P[0]
        le = P[np.argmax(np.linalg.norm(P - te, axis=1))]
        v = te - le
        tw.append(math.degrees(math.atan2(v[1], v[0])))
    tw = np.unwrap(np.radians(tw))
    return dict(self_intersections=int(n_x), twist_measured_deg=np.degrees(tw).tolist())
