#!/usr/bin/env python3
"""Параметрическое центробежное колесо компрессора: геометрия (Gmsh / Open CASCADE) и сетка
одного межлопаточного канала в формате SU2.

Узел платформы трека gas-turbine (routes/gtu/platform/geom). По данным одномерного расчёта
(расход, степень сжатия, частота вращения, D2; из gtu_explore.comp_stage или заданные) строит:
  * меридиональный контур втулки и покрышки - кубические кривые Безье;
  * закон углов лопатки beta(m), толщину, основные лопатки и полулопатки, загиб назад
    на выходе, радиальный зазор;
  * твёрдотельную модель колеса (STEP);
  * расчётную область одного межлопаточного канала (основная лопатка + полулопатка) с
    периодическими гранями, тетраэдрическую сетку и файл SU2 с маркерами
    inlet, outlet, hub, shroud, blade, periodic_1, periodic_2.

Система координат: ось z - ось вращения, поток входит вдоль +z; колесо вращается в сторону +theta
(вектор угловой скорости +z); periodic_2 = поворот periodic_1 на +360/Z вокруг z.
Запуск: python impeller.py --case core130 --level coarse
"""
import argparse
import json
import math
import os
import sys
import time

import numpy as np
from scipy.interpolate import CubicSpline
from scipy.optimize import brentq, minimize_scalar

HERE = os.path.dirname(os.path.abspath(__file__))
GTU_DIR = os.path.abspath(os.path.join(HERE, "..", ".."))
MARKERS = ("inlet", "outlet", "hub", "shroud", "blade", "periodic_1", "periodic_2")

# уровни сетки: размеры в долях ширины канала на выходе b2; n_gap - ячеек по высоте зазора
LEVELS = {
    "coarse": dict(h0=1 / 3, hb=1 / 7, n_gap=1),
    "medium": dict(h0=1 / 5, hb=1 / 12, n_gap=2),
    "fine": dict(h0=1 / 8, hb=1 / 20, n_gap=3),
}

# параметры по умолчанию (принятые, см. README)
DEFAULTS = dict(
    Z=9,                  # число основных лопаток
    splitters=True,       # полулопатки между основными
    m_split=0.30,         # начало полулопатки, доля меридиональной длины
    nu=0.35,              # D1h/D1s на входе
    incidence=0.0,        # угол атаки (лопатка минус поток) на входе, град
    beta2b=None,          # загиб назад на выходе, град от радиали; None - по формуле Визнера под psi
    beta_law=(0.0, 1.0),  # внутренние управляющие значения кубического закона F(m) (0..1)
    hub_bz=(0.5, 0.5),    # управляющие точки Безье втулки: доля осевой длины, доля радиального подъёма
    sh_bz=(0.5, 0.5),     # то же для покрышки
    t_hub=0.011,          # толщина лопатки у втулки / D2
    t_tip=0.0045,         # толщина лопатки у периферии / D2
    le_ratio=0.4,         # толщина на входной кромке / максимальной
    le_len=0.08,          # длина утонения у входной кромки, доля меридиональной длины
    te_ratio=0.7,         # толщина на выходной кромке / максимальной
    te_len=0.05,          # длина утонения у выходной кромки
    gap=None,             # радиальный зазор, м; None - как в comp_stage: max(0,3 мм; 0,0007 D2)
    rake=0.0,             # разность углов охвата втулки и покрышки на выходе, град (0 - кромка осевая)
    L_in=1.0,             # длина входного участка расчётной области / r1s
    r_out=1.30,           # радиус выхода расчётной области / r2
    t_back=0.03,          # толщина заднего диска / D2 (твёрдотельная модель)
    r_bore=0.35,          # радиус центрального отверстия / r1h (твёрдотельная модель)
    n_st=28,              # сечений лофта лопатки
    n_sheet=44,           # сечений лофта периодической поверхности
)


# ----------------------------------------------------------------------------- одномерные данные
_E = None


def load_explore():
    """gtu_explore (свойства газа Cantera и comp_stage). Модуль при импорте меняет cwd - возвращаем."""
    global _E
    if _E is None:
        cwd, dwb = os.getcwd(), sys.dont_write_bytecode
        if GTU_DIR not in sys.path:
            sys.path.insert(0, GTU_DIR)
        sys.dont_write_bytecode = True  # не писать __pycache__ в чужую папку
        try:
            import gtu_explore as E
        finally:
            os.chdir(cwd)
            sys.dont_write_bytecode = dwb
        _E = E
    return _E


def case_inputs(name):
    """Исходные данные именованного случая - из results-JSON трека, не вписаны руками."""
    if name == "core130":
        path = os.path.join(GTU_DIR, "gtu_explore_results.json")
        d = json.load(open(path, encoding="utf-8"))
        r = next(x for x in d["res"] if x.get("id") == 130)
        return dict(case=name, m=r["m"] / r["N"], PR=r["PR"], rpm=float(r["rpm"]), D2=r["Dc"],
                    material="al",
                    source=f"gtu_explore_results.json, схема 130 ({r['name']}): расход {r['m']} кг/с на "
                           f"{r['N']} ядра, PR {r['PR']}, {r['rpm']} об/мин, Dc {r['Dc']} м, U2 {r['U2']} м/с")
    if name == "taurus60":
        path = os.path.join(GTU_DIR, "gtu_dropin_results.json")
        d = json.load(open(path, encoding="utf-8"))
        r = next(x for x in d["res"] if x["arch"] == "1Ц" and x["cool"] == "N")
        t = r["T60"]
        return dict(case=name, m=t["G"], PR=r["PR"], rpm=None, D2=None, material="ti",
                    source=f"gtu_dropin_results.json, пакет Taurus 60, одно колесо (1Ц), охлаждаемый СА (N): "
                           f"расход {t['G']} кг/с, PR {r['PR']}; D2 и частота - gtu_explore.comp_stage "
                           f"(в JSON: D2 {t['D2']} м, {t['rpm']} об/мин)")
    raise ValueError(f"неизвестный случай {name!r}; есть: core130, taurus60")


def _bez(P, t):
    t = np.asarray(t, float)[..., None]
    return (1 - t) ** 3 * P[0] + 3 * t * (1 - t) ** 2 * P[1] + 3 * t ** 2 * (1 - t) * P[2] + t ** 3 * P[3]


def _bez1(c1, c2, x):
    """кубический закон 0..1 с управляющими значениями (0, c1, c2, 1) по равномерной абсциссе"""
    x = np.asarray(x, float)
    return 3 * x * (1 - x) ** 2 * c1 + 3 * x ** 2 * (1 - x) * c2 + x ** 3


def _pol(r, th, z):
    r, th, z = np.broadcast_arrays(np.asarray(r, float), np.asarray(th, float), np.asarray(z, float))
    return np.stack([r * np.cos(th), r * np.sin(th), z], axis=-1)


def _rotz(P, a):
    c, s = math.cos(a), math.sin(a)
    P = np.asarray(P, float)
    return np.stack([c * P[..., 0] - s * P[..., 1], s * P[..., 0] + c * P[..., 1], P[..., 2]], axis=-1)


class _Contour:
    """меридиональная кривая Безье (z, r) с параметризацией по нормированной длине дуги"""

    def __init__(self, P, n=4001):
        self.P = np.asarray(P, float)
        self.t = np.linspace(0.0, 1.0, n)
        self.zr = _bez(self.P, self.t)
        seg = np.hypot(*np.diff(self.zr, axis=0).T)
        self.s = np.concatenate([[0.0], np.cumsum(seg)])
        self.L = float(self.s[-1])

    def at(self, m):
        return _bez(self.P, np.interp(np.asarray(m, float) * self.L, self.s, self.t))


def _seg_dist(pts, poly):
    """расстояние от точек (N,2) до ломаной (M,2)"""
    A, B = poly[:-1], poly[1:]
    AB = B - A
    L2 = np.maximum((AB ** 2).sum(1), 1e-300)
    out = np.empty(len(pts))
    for i, p in enumerate(pts):
        t = np.clip(((p - A) * AB).sum(1) / L2, 0.0, 1.0)
        out[i] = np.sqrt((((A + t[:, None] * AB) - p) ** 2).sum(1)).min()
    return out


class Impeller:
    """Параметры, одномерная увязка и описание поверхностей колеса (без Gmsh)."""

    def __init__(self, m, PR, rpm=None, D2=None, case="custom", material="ti", source="задано вручную", **over):
        unknown = set(over) - set(DEFAULTS)
        if unknown:
            raise ValueError(f"неизвестные параметры: {sorted(unknown)}")
        self.p = dict(DEFAULTS)
        self.p.update({k: v for k, v in over.items() if v is not None or k in ("beta2b", "gap")})
        self.case, self.material, self.source = case, material, source
        self.warnings = []
        self._size_1d(m, PR, rpm, D2)
        self._contours()
        self._camber()

    # ------------------------------------------------------------------ 1D
    def _size_1d(self, m, PR, rpm, D2):
        E = load_explore()
        p = self.p
        T01, p01 = E.T0, E.P0
        st = E.comp_stage(T01, p01, PR, m, "ti")
        src = "comp_stage"
        if D2 is None and rpm is None:
            D2, rpm = st["D2"], st["rpm"]
        elif D2 is None:
            D2 = 60 * st["U2"] / (math.pi * rpm); src = "rpm задана, D2 из U2 comp_stage"
        elif rpm is None:
            rpm = 60 * st["U2"] / (math.pi * D2); src = "D2 задан, rpm из U2 comp_stage"
        else:
            src = "D2 и rpm заданы"
        r2 = D2 / 2
        Om = 2 * math.pi * rpm / 60
        U2 = Om * r2
        w, T2 = st["w"], st["T2"]
        psi = w / U2 ** 2
        RA = E.RA
        rho01 = p01 / (RA * T01)
        # плотность и Cm2/U2 на выходе - те же допущения, что в comp_stage
        rho2 = rho01 * PR ** 0.62 / ((T01 + 0.6 * (T2 - T01)) / T01)
        phi2 = 0.28
        b2 = m / (math.pi * D2 * rho2 * phi2 * U2)
        gap = p["gap"] if p["gap"] is not None else max(0.3e-3, 0.0007 * D2)

        # вход: D1s - минимум относительной скорости на периферии при заданном nu (без закрутки)
        nu = p["nu"]

        def stat(C):
            T1 = E.T_h(E.h(T01) - 0.5 * C * C)
            return T1, p01 * math.exp((E.s0(T1) - E.s0(T01)) / RA)

        def cp(T):
            return (E.h(T + 0.5) - E.h(T - 0.5))

        def flow(r1s):
            A = math.pi * r1s ** 2 * (1 - nu ** 2)
            Cs = 300.0
            for _ in range(30):  # скорость звука при M=1
                T1, _p = stat(Cs)
                g = cp(T1) / (cp(T1) - RA)
                Cs = math.sqrt(g * RA * T1)
            T1, p1 = stat(Cs)
            if p1 / (RA * T1) * Cs * A < m:
                return None
            f = lambda C: stat(C)[1] / (RA * stat(C)[0]) * C * A - m
            C = brentq(f, 1e-3, Cs)
            T1, p1 = stat(C)
            return C, T1, p1

        def W1s(r1s):
            fl = flow(r1s)
            return 1e9 if fl is None else math.hypot(fl[0], Om * r1s)

        grid = np.linspace(0.15 * r2, 0.95 * r2, 161)
        vals = [W1s(x) for x in grid]
        k = int(np.argmin(vals))
        res = minimize_scalar(W1s, bounds=(grid[max(k - 1, 0)], grid[min(k + 1, len(grid) - 1)]), method="bounded",
                              options=dict(xatol=1e-6))
        r1s = float(res.x)
        C1, T1, p1 = flow(r1s)
        r1h = nu * r1s
        g1 = cp(T1) / (cp(T1) - RA)
        a1 = math.sqrt(g1 * RA * T1)
        W1 = math.hypot(C1, Om * r1s)
        if r1s / r2 > 0.75:
            self.warnings.append(f"D1s/D2 = {r1s / r2:.3f} > 0,75: колесо с очень большим входом")
        if W1 / a1 > 1.3:
            self.warnings.append(f"относительное число Маха на периферии входа {W1 / a1:.2f} > 1,3: "
                                 "трансзвуковой вход, потери и запас по помпажу - проверять расчётом SU2")
        b1h = math.degrees(math.atan(Om * r1h / C1)) + p["incidence"]
        b1s = math.degrees(math.atan(Om * r1s / C1)) + p["incidence"]

        # выход: загиб назад по формуле Визнера sigma = 1 - sqrt(cos b)/Z^0.7 под psi = sigma (1 - phi2 tg b)
        Zex = p["Z"] * (2 if p["splitters"] else 1)
        sig = lambda b: 1 - math.sqrt(math.cos(math.radians(b))) / Zex ** 0.7
        psi_e = lambda b: sig(b) * (1 - phi2 * math.tan(math.radians(b)))
        if p["beta2b"] is None:
            if psi > psi_e(0.0):
                self.warnings.append(f"psi = {psi:.3f} больше, чем даёт радиальная лопатка ({psi_e(0):.3f}); загиб 0")
                b2b = 0.0
            elif psi < psi_e(60.0):
                self.warnings.append(f"psi = {psi:.3f} требует загиба > 60 град; принят 60")
                b2b = 60.0
            else:
                b2b = brentq(lambda b: psi_e(b) - psi, 0.0, 60.0)
        else:
            b2b = float(p["beta2b"])

        # осевая длина: корреляция Аунгира (по вторичным источникам)
        phi = m / (rho01 * math.pi * r2 ** 2 * U2)
        Lz = D2 * (0.014 + 0.023 * D2 / (2 * r1h) + 1.58 * phi)
        Ls = Lz - b2
        if Ls < 0.5 * (r2 - r1s):
            self.warnings.append("осевая длина покрышки мала, увеличена до 0,5 (r2 - r1s)")
            Ls = 0.5 * (r2 - r1s); Lz = Ls + b2

        self.d1 = dict(m=m, PR=PR, rpm=rpm, D2=D2, r2=r2, U2=U2, omega=Om, sizing=src,
                       T01=T01, p01=p01, w=w, T2=T2, psi=psi, eta_comp_stage=st["eta"],
                       comp_stage=dict(D2=st["D2"], rpm=st["rpm"], U2=st["U2"], b2=st["b2"]),
                       phi=phi, phi2=phi2, rho01=rho01, rho2=rho2, b2=b2, gap=gap,
                       r1s=r1s, r1h=r1h, D1s_D2=r1s / r2, C1=C1, T1=T1, p1=p1, W1s=W1, M1s_rel=W1 / a1,
                       beta1h=b1h, beta1s=b1s, beta2b=b2b, Z_exit=Zex, slip_wiesner=sig(b2b),
                       psi_euler=psi_e(b2b), Lz=Lz, Ls=Ls)
        if U2 > 620.0:
            self.warnings.append(f"U2 = {U2:.0f} м/с выше принятого предела титанового колеса 620 м/с")

    # ------------------------------------------------------------------ контуры
    def _contours(self):
        d, p = self.d1, self.p
        r1h, r1s, r2, Lz, Ls, b2 = d["r1h"], d["r1s"], d["r2"], d["Lz"], d["Ls"], d["b2"]
        kh1, kh2 = p["hub_bz"]
        ks1, ks2 = p["sh_bz"]
        self.Ph = np.array([[0, r1h], [kh1 * Lz, r1h], [Lz, r2 - kh2 * (r2 - r1h)], [Lz, r2]], float)
        self.Ps = np.array([[0, r1s], [ks1 * Ls, r1s], [Ls, r2 - ks2 * (r2 - r1s)], [Ls, r2]], float)
        self.hub, self.sh = _Contour(self.Ph), _Contour(self.Ps)
        self.z_in = -p["L_in"] * r1s
        self.r_out = p["r_out"] * r2
        far = 10 * r2
        self.hub_poly = np.vstack([[self.z_in - far, r1h], self.hub.zr, [Lz, self.r_out + far]])
        self.sh_poly = np.vstack([[self.z_in - far, r1s], self.sh.zr, [Ls, self.r_out + far]])
        # площадь сечения по квазиортогоналям (усечённые конусы между точками равной доли длины)
        mm = np.linspace(0, 1, 101)
        H, S = self.hub.at(mm), self.sh.at(mm)
        A = math.pi * (H[:, 1] + S[:, 1]) * np.hypot(*(S - H).T)
        self.area = dict(A1=float(A[0]), A2=float(A[-1]), A2_A1=float(A[-1] / A[0]),
                         min_over_A1=float(A.min() / A[0]), m_at_min=float(mm[A.argmin()]))
        if np.any(np.hypot(*(S - H).T) < 0.5 * b2):
            self.warnings.append("канал местами уже 0,5 b2 - проверить контуры")

    # ------------------------------------------------------------------ лопатка
    def _theta(self, cont, b_le, b_te, bump, mg):
        c1, c2 = self.p["beta_law"]
        beta = b_le + (b_te - b_le) * _bez1(c1, c2, mg) + bump * np.sin(np.pi * mg) ** 2
        zr = cont.at(mg)
        dth = -np.tan(np.radians(beta)) / zr[:, 1] * cont.L
        th = np.concatenate([[0.0], np.cumsum(0.5 * (dth[1:] + dth[:-1]) * np.diff(mg))])
        return beta, th

    def _camber(self):
        d, p = self.d1, self.p
        mg = np.linspace(0, 1, 2001)
        bs, ths = self._theta(self.sh, d["beta1s"], d["beta2b"], 0.0, mg)
        target = ths[-1] + math.radians(p["rake"])
        f = lambda bump: self._theta(self.hub, d["beta1h"], d["beta2b"], bump, mg)[1][-1] - target
        grid = np.linspace(-75.0, 45.0, 25)  # tg(beta) монотонен, пока beta < 90 град
        fv = [f(x) for x in grid]
        k = next((i for i in range(len(grid) - 1) if fv[i] * fv[i + 1] <= 0), None)
        if k is None:
            raise RuntimeError("не удаётся согласовать углы охвата втулки и покрышки на выходе (rake)")
        bump = brentq(f, grid[k], grid[k + 1])
        bh, thh = self._theta(self.hub, d["beta1h"], d["beta2b"], bump, mg)
        self.th_h, self.th_s = CubicSpline(mg, thh), CubicSpline(mg, ths)
        self.beta_h, self.beta_s, self.mg = bh, bs, mg
        self.bump = bump
        # производная угла по z на входной кромке (продолжение во входной участок)
        self.dthdz_h = -math.tan(math.radians(bh[0])) / d["r1h"]
        self.dthdz_s = -math.tan(math.radians(bs[0])) / d["r1s"]
        self.L_decay = abs(self.z_in) / 3
        self.geo = dict(wrap_hub_deg=float(-math.degrees(thh[-1])), wrap_shroud_deg=float(-math.degrees(ths[-1])),
                        hub_beta_bump_deg=float(bump), beta_hub_min_deg=float(bh.min()), beta_hub_max_deg=float(bh.max()),
                        beta_shroud_min_deg=float(bs.min()), beta_shroud_max_deg=float(bs.max()),
                        L_hub=self.hub.L, L_shroud=self.sh.L)
        if bh.min() < -5.0:
            self.warnings.append(f"угол лопатки у втулки уходит до {bh.min():.1f} град")

    def merid(self, which, m):
        """(z, r, theta) средней линии на втулке/покрышке для доли длины m (m<0 - вход, m>1 - за колесом)"""
        d = self.d1
        m = np.atleast_1d(np.asarray(m, float))
        cont, th, dthdz, zend = ((self.hub, self.th_h, self.dthdz_h, d["Lz"]) if which == "hub"
                                 else (self.sh, self.th_s, self.dthdz_s, d["Ls"]))
        z = np.empty_like(m); r = np.empty_like(m); t = np.empty_like(m)
        a, b, c = m < 0, (m >= 0) & (m <= 1), m > 1
        if a.any():
            z[a] = m[a] * cont.L
            r[a] = cont.P[0, 1]
            t[a] = dthdz * self.L_decay * (np.exp(z[a] / self.L_decay) - 1)
        if b.any():
            zr = cont.at(m[b]); z[b], r[b] = zr[:, 0], zr[:, 1]
            t[b] = th(m[b])
        if c.any():
            r[c] = d["r2"] * m[c]
            z[c] = zend
            t[c] = th(1.0) - math.tan(math.radians(d["beta2b"])) * np.log(m[c])
        return z, r, t

    def hs(self, m):
        """точки средней линии на втулке и покрышке, (N,3) каждая"""
        zh, rh, th = self.merid("hub", m)
        zs, rs, ts = self.merid("shroud", m)
        return _pol(rh, th, zh), _pol(rs, ts, zs)

    def normal(self, m, dm=1e-5):
        """нормаль к поверхности средней линии в середине размаха, со стороны +theta (давление)"""
        H, S = self.hs([m - dm, m, m + dm])
        C = 0.5 * (H + S)
        dCm = C[2] - C[0]
        dCs = S[1] - H[1]
        n = np.cross(dCm, dCs)
        n /= np.linalg.norm(n)
        c = C[1]
        et = np.array([-c[1], c[0], 0.0]) / math.hypot(c[0], c[1])
        return n if n @ et > 0 else -n

    def thickness(self, m, s, m0):
        p, D2 = self.p, self.d1["D2"]
        t = (p["t_hub"] + (p["t_tip"] - p["t_hub"]) * s) * D2
        x = (m - m0) / p["le_len"]
        g = 1.0 if x >= 1 else p["le_ratio"] + (1 - p["le_ratio"]) * (1 - (1 - max(x, 0.0)) ** 2)
        y = (1 - m) / p["te_len"]
        g2 = p["te_ratio"] if y <= 0 else (1.0 if y >= 1 else p["te_ratio"] + (1 - p["te_ratio"]) * (1 - (1 - y) ** 2))
        return t * g * g2

    def blade_sections(self, m0=0.0, rot=0.0):
        """плоские четырёхугольные сечения лопатки (корень ниже втулки, вершина выше покрышки)"""
        d, p = self.d1, self.p
        e_root = 0.25 * d["b2"]
        e_far = 0.25 * d["b2"]
        m_end = 1 + 1.5 * p["t_hub"] * d["D2"] / d["r2"]
        n = p["n_st"] if m0 == 0 else max(12, int(round(p["n_st"] * (1 - m0))))
        u = np.linspace(0, 1, n)
        u = 0.4 * u + 0.6 * 0.5 * (1 - np.cos(np.pi * u))
        ms = m0 + (1 - m0) * u
        ms = np.concatenate([ms, [1 + 0.5 * (m_end - 1), m_end]])
        secs = []
        for m in ms:
            H, S = self.hs(m)
            H, S = H[0], S[0]
            dv = S - H
            L = np.linalg.norm(dv)
            nv = self.normal(m)
            sr, sf = -e_root / L, 1 + e_far / L
            tr, tf = self.thickness(min(m, 1.0), sr, m0), self.thickness(min(m, 1.0), sf, m0)
            if tf <= 0 or tr <= 0:
                raise RuntimeError("отрицательная толщина на продолжении лопатки")
            Cr, Cf = H + sr * dv, H + sf * dv
            q = np.array([Cr + 0.5 * tr * nv, Cf + 0.5 * tf * nv, Cf - 0.5 * tf * nv, Cr - 0.5 * tr * nv])
            secs.append(_rotz(q, rot))
        return secs

    def sheet_lines(self, phi):
        """отрезки периодической поверхности: средняя линия, повёрнутая на phi, с запасом за стенки"""
        d, p = self.d1, self.p
        e = 0.25 * d["b2"]
        mz = 0.05 * abs(self.z_in)
        L_in, L_mid, L_out = abs(self.z_in) + mz, 0.5 * (self.hub.L + self.sh.L), (self.r_out - d["r2"]) * 1.1
        tot = L_in + L_mid + L_out
        n = p["n_sheet"]
        n_in, n_out = max(4, round(n * L_in / tot)), max(4, round(n * L_out / tot))
        n_mid = max(10, n - n_in - n_out)
        stations = []
        for z in np.linspace(self.z_in - mz, 0, n_in, endpoint=False):
            H = _pol(d["r1h"], self.dthdz_h * self.L_decay * (math.exp(z / self.L_decay) - 1), z)
            S = _pol(d["r1s"], self.dthdz_s * self.L_decay * (math.exp(z / self.L_decay) - 1), z)
            stations.append((H, S))
        for m in np.linspace(0, 1, n_mid):
            H, S = self.hs(m)
            stations.append((H[0], S[0]))
        for r in np.linspace(d["r2"], d["r2"] + L_out, n_out + 1)[1:]:
            H, S = self.hs(r / d["r2"])
            stations.append((H[0], S[0]))
        out = []
        for H, S in stations:
            u = (S - H) / np.linalg.norm(S - H)
            out.append(_rotz(np.array([H - e * u, S + e * u]), phi))
        return out

    def offset_shroud(self, n=81):
        """кривая покрышки, смещённая внутрь на зазор (z, r)"""
        g = self.d1["gap"]
        t = np.linspace(0, 1, n)
        P = self.Ps
        zr = _bez(P, t)
        dz = 3 * (1 - t) ** 2 * (P[1] - P[0])[:, None] + 6 * (1 - t) * t * (P[2] - P[1])[:, None] + 3 * t ** 2 * (P[3] - P[2])[:, None]
        T = dz.T / np.linalg.norm(dz.T, axis=1)[:, None]
        N = np.stack([T[:, 1], -T[:, 0]], axis=1)  # к втулке: (+z, -r) на входе
        return zr + g * N

    def summary(self):
        return dict(case=self.case, source=self.source, material=self.material,
                    params={k: (list(v) if isinstance(v, tuple) else v) for k, v in self.p.items()},
                    d1=self.d1, area=self.area, geometry=self.geo,
                    contours=dict(hub_bezier_zr=self.Ph.tolist(), shroud_bezier_zr=self.Ps.tolist(),
                                  z_in=self.z_in, r_out=self.r_out),
                    warnings=list(self.warnings))


# ----------------------------------------------------------------------------- Gmsh / OCC
def _gmsh():
    import gmsh
    return gmsh


def _pt(occ, P):
    return occ.addPoint(float(P[0]), float(P[1]), float(P[2]))


def _merid_pt(occ, z, r):
    return occ.addPoint(float(r), 0.0, float(z))


def _loft_solid(occ, secs):
    wires = []
    for q in secs:
        pts = [_pt(occ, x) for x in q]
        ls = [occ.addLine(pts[i], pts[(i + 1) % 4]) for i in range(4)]
        wires.append(occ.addWire(ls))
    out = occ.addThruSections(wires, makeSolid=True, makeRuled=False)
    return [t for dd, t in out if dd == 3][0]


def _loft_sheet(occ, lines):
    wires = []
    for a, b in lines:
        wires.append(occ.addWire([occ.addLine(_pt(occ, a), _pt(occ, b))]))
    out = occ.addThruSections(wires, makeSolid=False, makeRuled=False)
    return [t for dd, t in out if dd == 2][0]


def _revolve_face(occ, curves, th0, dth):
    loop = occ.addCurveLoop(curves)
    f = occ.addPlaneSurface([loop])
    if th0:
        occ.rotate([(2, f)], 0, 0, 0, 0, 0, 1, th0)
    out = occ.revolve([(2, f)], 0, 0, 0, 0, 0, 1, dth)
    return [t for dd, t in out if dd == 3][0]


def _bezier_curve(occ, P, reverse=False):
    Q = P[::-1] if reverse else P
    return occ.addBezier([_merid_pt(occ, z, r) for z, r in Q])


def _envelope(occ, imp, th0, dth):
    """тело вращения, которым обрезается лопатка: внутри покрышки минус зазор и внутри r2"""
    d = imp.d1
    zr = imp.offset_shroud()
    z_a = imp.z_in
    z_b = d["Lz"] + 0.25 * d["b2"] + 0.05 * d["D2"]
    r_in = 0.3 * d["r1h"]
    p0 = _merid_pt(occ, z_a, r_in)
    p1 = _merid_pt(occ, z_a, zr[0, 1])
    sp = [_merid_pt(occ, z, r) for z, r in zr]
    p4 = _merid_pt(occ, z_b, d["r2"])
    p5 = _merid_pt(occ, z_b, r_in)
    cs = [occ.addLine(p0, p1), occ.addLine(p1, sp[0]), occ.addSpline(sp), occ.addLine(sp[-1], p4),
          occ.addLine(p4, p5), occ.addLine(p5, p0)]
    return _revolve_face(occ, cs, th0, dth)


def _blade_solids(occ, imp, rot=0.0):
    """основная лопатка и полулопатка (обрезанные по зазору и r2), повёрнутые на rot"""
    p = imp.p
    secs_main = imp.blade_sections(0.0, rot)
    blades = [_loft_solid(occ, secs_main)]
    if p["splitters"]:
        blades.append(_loft_solid(occ, imp.blade_sections(p["m_split"], rot + math.pi / p["Z"])))
    allp = np.vstack([np.vstack(imp.blade_sections(0.0, rot))])
    th = np.unwrap(np.arctan2(allp[:, 1], allp[:, 0]))
    lo, hi = th.min() - 0.3, th.max() + 0.3 + (math.pi / p["Z"] if p["splitters"] else 0.0)
    env = _envelope(occ, imp, lo, hi - lo)
    out, _ = occ.intersect([(3, b) for b in blades], [(3, env)])
    vols = [t for dd, t in out if dd == 3]
    if len(vols) != len(blades):
        raise RuntimeError(f"обрезка лопаток дала {len(vols)} тел вместо {len(blades)}")
    return vols


def _keep_only(gmsh, vols):
    occ = gmsh.model.occ
    occ.synchronize()
    keep = {(3, v) for v in vols}
    f = {(2, abs(t)) for _, t in gmsh.model.getBoundary([(3, v) for v in vols], combined=False, oriented=False)}
    c = {(1, abs(t)) for _, t in gmsh.model.getBoundary(list(f), combined=False, oriented=False)}
    pt = {(0, abs(t)) for _, t in gmsh.model.getBoundary(list(c), combined=False, oriented=False)}
    keep |= f | c | pt
    for dim in (3, 2, 1, 0):
        rm = [e for e in gmsh.model.getEntities(dim) if e not in keep]
        if rm:
            occ.remove(rm)
    occ.synchronize()


def _init(gmsh, threads, verbose):
    gmsh.initialize(interruptible=False)
    gmsh.option.setNumber("General.Terminal", 1 if verbose else 0)
    gmsh.option.setNumber("General.Verbosity", 5 if verbose else 2)
    gmsh.option.setNumber("General.NumThreads", threads)
    gmsh.option.setNumber("Geometry.OCCParallel", 1)
    gmsh.logger.start()


def build_solid(imp, path, threads=4, verbose=False):
    """твёрдотельная модель колеса: втулка с задним диском + Z основных лопаток + Z полулопаток -> STEP"""
    gmsh = _gmsh()
    _init(gmsh, threads, verbose)
    try:
        occ = gmsh.model.occ
        gmsh.model.add("impeller")
        d, p = imp.d1, imp.p
        t_back = max(p["t_back"] * d["D2"], 1.5 * 0.25 * d["b2"])
        z_n, z_bk, r_b = -0.3 * d["r1h"], d["Lz"] + t_back, p["r_bore"] * d["r1h"]
        P0, P1 = _merid_pt(occ, z_n, r_b), _merid_pt(occ, z_n, d["r1h"])
        hb = _bezier_curve(occ, imp.Ph)
        occ.synchronize()
        bb = gmsh.model.getBoundary([(1, hb)], oriented=False)
        pa, pb = bb[0][1], bb[1][1]
        P4, P5 = _merid_pt(occ, z_bk, d["r2"]), _merid_pt(occ, z_bk, r_b)
        cs = [occ.addLine(P0, P1), occ.addLine(P1, pa), hb, occ.addLine(pb, P4), occ.addLine(P4, P5), occ.addLine(P5, P0)]
        hub = _revolve_face(occ, cs, 0.0, 2 * math.pi)
        base = _blade_solids(occ, imp)
        # объём по частям: втулка + Z (лопатки - их часть внутри втулки). getMass сплавленного тела
        # ненадёжен (грань втулки с 2Z вырезами интегрируется с ошибкой ~0,5 %), по частям - сверено
        # с формулой Паппа-Гульдина и с объёмом тетраэдрической сетки повторно прочитанного STEP
        occ.synchronize()
        V_hub = occ.getMass(3, hub)
        V_bl = sum(occ.getMass(3, b) for b in base)
        ov, _ = occ.intersect(occ.copy([(3, b) for b in base]), occ.copy([(3, hub)]))
        occ.synchronize()
        V_ov = sum(occ.getMass(3, t) for dd, t in ov if dd == 3)
        occ.remove(ov, recursive=True)
        V = V_hub + p["Z"] * (V_bl - V_ov)
        tools = [(3, b) for b in base]
        for k in range(1, p["Z"]):
            cp = occ.copy([(3, b) for b in base])
            occ.rotate(cp, 0, 0, 0, 0, 0, 1, 2 * math.pi * k / p["Z"])
            tools += cp
        out, _ = occ.fuse([(3, hub)], tools)
        vols = [t for dd, t in out if dd == 3]
        _keep_only(gmsh, vols)
        if len(vols) != 1:
            raise RuntimeError(f"твёрдое тело колеса: {len(vols)} тел вместо 1")
        bbox = gmsh.model.getBoundingBox(3, vols[0])
        gmsh.write(path)
        return dict(step=path, volume_m3=V, mass_kg=dict(ti_6al_4v=V * 4430.0, al_2618=V * 2760.0),
                    volume_parts_m3=dict(hub=V_hub, blade_pair=V_bl, blade_pair_inside_hub=V_ov),
                    bbox_occ=list(bbox), t_back=t_back, r_bore=r_b, z_nose=z_n,
                    n_faces=len(gmsh.model.getEntities(2)))
    finally:
        gmsh.logger.stop()
        gmsh.finalize()


def _edge_points(gmsh, face, n=5):
    pts = []
    for _, c in gmsh.model.getBoundary([(2, face)], combined=False, oriented=False):
        c = abs(c)
        lo, hi = gmsh.model.getParametrizationBounds(1, c)
        u = np.linspace(lo[0], hi[0], n + 2)[1:-1]
        pts.append(np.asarray(gmsh.model.getValue(1, c, list(u))).reshape(-1, 3))
    return np.vstack(pts) if pts else np.zeros((0, 3))


def _classify(gmsh, imp, vol, refs):
    """маркеры границ по геометрии: точки на рёбрах каждой грани проверяются на принадлежность поверхностям"""
    d = imp.d1
    tol = max(1e-6, 2e-4 * d["b2"])
    gap = d["gap"]
    groups = {k: [] for k in MARKERS}
    tips, les, tes = [], [], []
    diag = {}
    for _, f in gmsh.model.getBoundary([(3, vol)], oriented=False):
        f = abs(f)
        P = _edge_points(gmsh, f)
        z = P[:, 2]
        r = np.hypot(P[:, 0], P[:, 1])
        zr = np.stack([z, r], 1)
        dh = _seg_dist(zr, imp.hub_poly)
        ds = _seg_dist(zr, imp.sh_poly)
        if np.all(np.abs(z - imp.z_in) < tol):
            g = "inlet"
        elif np.all(np.abs(r - imp.r_out) < tol):
            g = "outlet"
        elif np.all(dh < tol):
            g = "hub"
        elif np.all(ds < tol):
            g = "shroud"
        else:
            dist = []
            for rt in refs:
                dist.append(max(np.linalg.norm(np.asarray(gmsh.model.getClosestPoint(2, rt, list(x))[0]) - x) for x in P))
            if dist[0] < tol:
                g = "periodic_1"
            elif dist[1] < tol:
                g = "periodic_2"
            else:
                g = "blade"
                if np.all(np.abs(ds - gap) < max(tol, 0.02 * gap)):
                    tips.append(f)
                elif np.all(np.abs(r - d["r2"]) < tol):
                    tes.append(f)
                elif gmsh.model.getType(2, f) == "Plane":
                    les.append(f)
        groups[g].append(f)
    return groups, dict(tip=tips, le=les, te=tes), tol


def build_passage(imp, level, path_su2, path_step=None, path_msh=None, threads=4, verbose=False):
    """область одного межлопаточного канала, сетка и запись SU2"""
    gmsh = _gmsh()
    L = LEVELS[level]
    d, p = imp.d1, imp.p
    Z = p["Z"]
    pitch = 2 * math.pi / Z
    # середина между основной лопаткой (0) и полулопаткой (pitch/2) соседнего канала
    phi1 = -0.25 * pitch if p["splitters"] else -0.5 * pitch
    phi2 = phi1 + pitch
    t0 = time.time()
    _init(gmsh, threads, verbose)
    try:
        occ = gmsh.model.occ
        gmsh.model.add("passage")
        # 1) тело вращения проточной части на угловом секторе, перекрывающем канал
        lines1 = imp.sheet_lines(phi1)
        allp = np.vstack(lines1 + imp.sheet_lines(phi2))
        th = np.unwrap(np.arctan2(allp[:, 1], allp[:, 0]))
        lo, hi = th.min() - 0.3, th.max() + 0.3
        if hi - lo >= 2 * math.pi - 0.05:
            raise RuntimeError("канал охватывает больше полного оборота; уменьшите L_in")
        A, B = _merid_pt(occ, imp.z_in, d["r1h"]), _merid_pt(occ, 0.0, d["r1h"])
        hb = _bezier_curve(occ, imp.Ph)
        sb = _bezier_curve(occ, imp.Ps, reverse=True)
        occ.synchronize()
        hb0, hb1 = [t for _, t in gmsh.model.getBoundary([(1, hb)], oriented=False)]
        sb0, sb1 = [t for _, t in gmsh.model.getBoundary([(1, sb)], oriented=False)]
        D, Ee = _merid_pt(occ, d["Lz"], imp.r_out), _merid_pt(occ, d["Ls"], imp.r_out)
        Hh = _merid_pt(occ, imp.z_in, d["r1s"])
        cs = [occ.addLine(A, hb0), hb, occ.addLine(hb1, D), occ.addLine(D, Ee), occ.addLine(Ee, sb0), sb,
              occ.addLine(sb1, Hh), occ.addLine(Hh, A)]
        body = _revolve_face(occ, cs, lo, hi - lo)
        # 2) периодические поверхности: средняя линия, повёрнутая на phi1, и её копия на +360/Z
        s1 = _loft_sheet(occ, lines1)
        s2 = occ.copy([(2, s1)])[0][1]
        occ.rotate([(2, s2)], 0, 0, 0, 0, 0, 1, pitch)
        refs = [t for _, t in occ.copy([(2, s1), (2, s2)])]
        out, omap = occ.fragment([(3, body)], [(2, s1), (2, s2)])
        occ.synchronize()
        vols = [t for dd, t in out if dd == 3]
        m1 = {t for dd, t in omap[1]}
        m2 = {t for dd, t in omap[2]}
        mid = None
        for v in vols:
            bf = {abs(t) for _, t in gmsh.model.getBoundary([(3, v)], oriented=False)}
            if bf & m1 and bf & m2:
                mid = v
        if mid is None:
            raise RuntimeError(f"не найден объём между периодическими поверхностями ({len(vols)} объёмов)")
        # 3) вычесть лопатки
        blades = _blade_solids(occ, imp)
        out, _ = occ.cut([(3, mid)], [(3, b) for b in blades])
        pv = [t for dd, t in out if dd == 3]
        if len(pv) != 1:
            raise RuntimeError(f"после вычитания лопаток {len(pv)} объёмов вместо 1")
        vol = pv[0]
        _keep_only(gmsh, [vol] + [])
        # опорные поверхности для классификации: снова строим копии (удалены _keep_only)
        rs1 = _loft_sheet(occ, lines1)
        rs2 = occ.copy([(2, rs1)])[0][1]
        occ.rotate([(2, rs2)], 0, 0, 0, 0, 0, 1, pitch)
        occ.synchronize()
        groups, special, tol = _classify(gmsh, imp, vol, [rs1, rs2])
        occ.remove([(2, rs1), (2, rs2)], recursive=True)
        _keep_only(gmsh, [vol])
        empty = [k for k in MARKERS if not groups[k]]
        if empty:
            raise RuntimeError(f"нет граней для маркеров {empty}")
        if len(groups["periodic_1"]) != len(groups["periodic_2"]):
            raise RuntimeError("число периодических граней не совпадает")
        t_geo = time.time() - t0
        if path_step:
            gmsh.write(path_step)
        # 4) периодичность и размеры элементов
        c, s = math.cos(pitch), math.sin(pitch)
        T = [c, -s, 0, 0, s, c, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1]
        gmsh.model.mesh.setPeriodic(2, groups["periodic_2"], groups["periodic_1"], T)
        b2, gap = d["b2"], d["gap"]
        h0, hb = L["h0"] * b2, L["hb"] * b2
        h_gap = gap / L["n_gap"]
        t_tip = p["t_tip"] * d["D2"]
        t_le = p["le_ratio"] * t_tip
        h_le = max(t_le, h_gap)
        fld = gmsh.model.mesh.field
        others = [f for f in groups["blade"] if f not in special["tip"]]
        f1 = fld.add("Distance"); fld.setNumbers(f1, "SurfacesList", others); fld.setNumber(f1, "Sampling", 80)
        t1 = fld.add("Threshold"); fld.setNumber(t1, "InField", f1)
        fld.setNumber(t1, "SizeMin", hb); fld.setNumber(t1, "SizeMax", h0)
        fld.setNumber(t1, "DistMin", 0.5 * hb); fld.setNumber(t1, "DistMax", 0.5 * hb + 3 * h0)
        fields = [t1]
        # зазор: расстояние до рёбер вершины лопатки
        tip_curves = sorted({abs(cc) for _, cc in gmsh.model.getBoundary([(2, f) for f in special["tip"]],
                                                                          combined=False, oriented=False)})
        if tip_curves:
            f2 = fld.add("Distance"); fld.setNumbers(f2, "CurvesList", tip_curves); fld.setNumber(f2, "Sampling", 1500)
            t2 = fld.add("Threshold"); fld.setNumber(t2, "InField", f2)
            dmin = math.hypot(0.5 * t_tip, gap) + 0.5 * h_gap
            fld.setNumber(t2, "SizeMin", h_gap); fld.setNumber(t2, "SizeMax", hb)
            fld.setNumber(t2, "DistMin", dmin); fld.setNumber(t2, "DistMax", dmin + 3 * hb)
            fields.append(t2)
        edge_curves = sorted({abs(cc) for _, cc in gmsh.model.getBoundary([(2, f) for f in special["le"] + special["te"]],
                                                                           combined=False, oriented=False)})
        if edge_curves:
            f3 = fld.add("Distance"); fld.setNumbers(f3, "CurvesList", edge_curves); fld.setNumber(f3, "Sampling", 400)
            t3 = fld.add("Threshold"); fld.setNumber(t3, "InField", f3)
            fld.setNumber(t3, "SizeMin", h_le); fld.setNumber(t3, "SizeMax", hb)
            fld.setNumber(t3, "DistMin", 0.5 * h_le); fld.setNumber(t3, "DistMax", 0.5 * h_le + 3 * hb)
            fields.append(t3)
        fm = fld.add("Min"); fld.setNumbers(fm, "FieldsList", fields)
        fld.setAsBackgroundMesh(fm)
        for k, v in (("Mesh.MeshSizeExtendFromBoundary", 0), ("Mesh.MeshSizeFromPoints", 0),
                     ("Mesh.MeshSizeFromCurvature", 0), ("Mesh.MeshSizeMin", 0.5 * h_gap),
                     ("Mesh.MeshSizeMax", h0), ("Mesh.Algorithm", 6), ("Mesh.Algorithm3D", 10),
                     ("Mesh.Optimize", 1), ("Mesh.MaxNumThreads3D", threads)):
            gmsh.option.setNumber(k, v)
        t1m = time.time()
        gmsh.model.mesh.generate(3)
        if _quality(gmsh)["minSICN"] < 0.05:
            gmsh.model.mesh.optimize("Netgen")
        snap = _snap_periodic(gmsh, groups["periodic_2"])
        q = _quality(gmsh)
        t_mesh = time.time() - t1m
        if q["minSICN"] <= 0:
            raise RuntimeError(f"вырожденные ячейки: minSICN = {q['minSICN']}")
        # 5) физические группы и запись
        gmsh.model.addPhysicalGroup(3, [vol], name="fluid")
        for k in MARKERS:
            gmsh.model.addPhysicalGroup(2, groups[k], name=k)
        gmsh.write(path_su2)
        if path_msh:
            gmsh.write(path_msh)
        stats = dict(level=level, h0=h0, h_blade=hb, h_gap=h_gap, h_edge=h_le, tol_classify=tol,
                     faces={k: groups[k] for k in MARKERS}, tip_faces=special["tip"], le_faces=special["le"],
                     te_faces=special["te"], quality=q, periodic_snap=snap,
                     t_geometry_s=round(t_geo, 1), t_mesh_s=round(t_mesh, 1),
                     phi1_rad=phi1, phi2_rad=phi2, sector_lo_hi_rad=[lo, hi],
                     gmsh_warnings=[w for w in gmsh.logger.get() if "Warning" in w][:50])
        return stats
    finally:
        gmsh.logger.stop()
        gmsh.finalize()


def _snap_periodic(gmsh, slaves):
    """узлы ведомой периодической грани (и её рёбер, вершин) ставятся точно в повёрнутые узлы ведущей.
    Gmsh переносит сетку поворотом, но затем проецирует узлы на ведомую грань, обрезанную отдельно;
    расхождение порядка допуска OCC (1e-8 м). Возвращает число узлов и наибольший сдвиг."""
    ents = {(2, s) for s in slaves}
    cur = {(1, abs(c)) for _, c in gmsh.model.getBoundary([(2, s) for s in slaves], combined=False, oriented=False)}
    ents |= cur
    ents |= {(0, abs(v)) for _, v in gmsh.model.getBoundary(list(cur), combined=False, oriented=False)}
    moved, dmax = {}, 0.0
    for dim, tag in sorted(ents, reverse=True):
        try:
            tm, nodes, nodes_m, tfo = gmsh.model.mesh.getPeriodicNodes(dim, tag)
        except Exception:  # noqa - сущность не периодическая
            continue
        if tm == tag or not len(nodes):
            continue
        A = np.asarray(tfo, float).reshape(4, 4)
        for n, nm in zip(nodes, nodes_m):
            if n in moved:
                continue
            xm = np.asarray(gmsh.model.mesh.getNode(nm)[0])
            x, uv, _, _ = gmsh.model.mesh.getNode(n)
            xn = A[:3, :3] @ xm + A[:3, 3]
            dmax = max(dmax, float(np.linalg.norm(xn - np.asarray(x))))
            gmsh.model.mesh.setNode(n, list(xn), list(uv))
            moved[n] = True
    return dict(n_nodes=len(moved), max_shift_m=dmax)


def _quality(gmsh):
    typ, tags, _ = gmsh.model.mesh.getElements(3)
    allt = np.concatenate([np.asarray(t) for t in tags]) if tags else np.zeros(0)
    sicn = np.asarray(gmsh.model.mesh.getElementQualities(list(allt), "minSICN"))
    gam = np.asarray(gmsh.model.mesh.getElementQualities(list(allt), "gamma"))
    nodes = gmsh.model.mesh.getNodes()[0]
    return dict(n_cells=int(len(allt)), cell_types=[int(t) for t in typ], n_nodes=int(len(nodes)),
                minSICN=float(sicn.min()), p01SICN=float(np.percentile(sicn, 1)), meanSICN=float(sicn.mean()),
                min_gamma=float(gam.min()), mean_gamma=float(gam.mean()),
                n_SICN_below_0_1=int((sicn < 0.1).sum()))


# ----------------------------------------------------------------------------- верхний уровень
def build(case="core130", level="coarse", out_dir=None, m=None, PR=None, rpm=None, D2=None,
          step=True, mesh=True, msh=False, threads=4, verbose=False, **over):
    if case == "custom":
        if m is None or PR is None:
            raise ValueError("для --case custom нужны --m и --pr")
        inp = dict(case="custom", m=m, PR=PR, rpm=rpm, D2=D2, material=over.pop("material", "ti"), source="задано вручную")
    else:
        inp = case_inputs(case)
        for k, v in (("m", m), ("PR", PR), ("rpm", rpm), ("D2", D2)):
            if v is not None:
                inp[k] = v; inp["source"] += f"; {k} переопределён: {v}"
    imp = Impeller(**inp, **over)
    out_dir = out_dir or os.path.join(HERE, "out", f"{case}_{level}")
    os.makedirs(out_dir, exist_ok=True)
    res = imp.summary()
    res["software"] = dict(python=sys.version.split()[0], numpy=np.__version__)
    try:
        g = _gmsh(); res["software"]["gmsh"] = g.__version__
    except Exception as e:  # noqa
        res["software"]["gmsh"] = f"не импортируется: {e}"
    res["files"] = {}
    t0 = time.time()
    if step:
        f = os.path.join(out_dir, f"impeller_{case}.step")
        res["solid"] = build_solid(imp, f, threads, verbose)
        res["files"]["impeller_step"] = f
    if mesh:
        su2 = os.path.join(out_dir, f"passage_{case}_{level}.su2")
        stp = os.path.join(out_dir, f"passage_{case}.step")
        mshf = os.path.join(out_dir, f"passage_{case}_{level}.msh") if msh else None
        res["mesh"] = build_passage(imp, level, su2, stp, mshf, threads, verbose)
        res["files"].update(passage_su2=su2, passage_step=stp)
        if mshf:
            res["files"]["passage_msh"] = mshf
        res["su2"] = dict(markers=list(MARKERS), volume="fluid",
                          periodic=dict(master="periodic_1", slave="periodic_2", axis=[0, 0, 1], center=[0, 0, 0],
                                        angle_deg=360.0 / imp.p["Z"],
                                        note="periodic_2 = поворот periodic_1 на +360/Z вокруг +z"),
                          rotation_rate_rad_s=[0.0, 0.0, imp.d1["omega"]])
    res["t_total_s"] = round(time.time() - t0, 1)
    jf = os.path.join(out_dir, f"impeller_{case}_{level}.json")
    with open(jf, "w", encoding="utf-8") as fh:
        json.dump(res, fh, ensure_ascii=False, indent=1, default=float)
    res["files"]["json"] = jf
    return res


def record(res, path=os.path.join(HERE, "impeller_results.json")):
    """краткая запись прогона в impeller_results.json папки (без абсолютных путей и больших файлов)"""
    db = {}
    if os.path.exists(path):
        with open(path, encoding="utf-8") as fh:
            db = json.load(fh)
    db.setdefault("method", "Параметрическое колесо по 1D-данным gtu_explore.comp_stage; Gmsh/Open CASCADE; "
                            "тетраэдрическая сетка одного межлопаточного канала; см. README.md")
    db.setdefault("runs", {})
    key = f"{res['case']}_{res['mesh']['level'] if 'mesh' in res else 'solid'}"
    keep = {k: v for k, v in res.items() if k not in ("files",)}
    if "mesh" in keep:
        keep["mesh"] = {k: v for k, v in keep["mesh"].items() if k != "gmsh_warnings"}
    if "solid" in keep:
        keep["solid"] = {k: v for k, v in keep["solid"].items() if k != "step"}
    keep["date"] = time.strftime("%Y-%m-%d")
    keep["files"] = {k: os.path.basename(v) for k, v in res["files"].items()}
    db["runs"][key] = keep
    with open(path, "w", encoding="utf-8") as fh:
        json.dump(db, fh, ensure_ascii=False, indent=1, default=float)
    return path


def main(argv=None):
    ap = argparse.ArgumentParser(description="Параметрическое центробежное колесо: STEP и сетка канала SU2")
    ap.add_argument("--case", default="core130", help="core130 | taurus60 | custom")
    ap.add_argument("--level", default="coarse", choices=list(LEVELS))
    ap.add_argument("--out", default=None, help="папка результатов (по умолчанию out/<case>_<level>)")
    ap.add_argument("--m", type=float, help="расход, кг/с")
    ap.add_argument("--pr", type=float, help="степень сжатия")
    ap.add_argument("--rpm", type=float, help="частота вращения, об/мин")
    ap.add_argument("--d2", type=float, help="диаметр колеса D2, м")
    ap.add_argument("--z", type=int, help="число основных лопаток")
    ap.add_argument("--no-splitters", action="store_true")
    ap.add_argument("--m-split", type=float)
    ap.add_argument("--beta2b", type=float, help="загиб назад, град от радиали")
    ap.add_argument("--gap", type=float, help="радиальный зазор, м")
    ap.add_argument("--incidence", type=float)
    ap.add_argument("--no-step", action="store_true")
    ap.add_argument("--no-mesh", action="store_true")
    ap.add_argument("--msh", action="store_true", help="также записать .msh")
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--verbose", action="store_true")
    ap.add_argument("--record", action="store_true", help="дописать краткий итог в impeller_results.json")
    a = ap.parse_args(argv)
    over = dict(Z=a.z, m_split=a.m_split, incidence=a.incidence)
    over = {k: v for k, v in over.items() if v is not None}
    if a.beta2b is not None:
        over["beta2b"] = a.beta2b
    if a.gap is not None:
        over["gap"] = a.gap
    if a.no_splitters:
        over["splitters"] = False
    res = build(a.case, a.level, a.out, a.m, a.pr, a.rpm, a.d2, step=not a.no_step, mesh=not a.no_mesh,
                msh=a.msh, threads=a.threads, verbose=a.verbose, **over)
    d = res["d1"]
    print(f"[{res['case']}] m={d['m']:.2f} кг/с PR={d['PR']:.2f} n={d['rpm']:.0f} об/мин D2={d['D2']:.4f} м "
          f"U2={d['U2']:.0f} м/с b2={d['b2'] * 1e3:.1f} мм D1s/D2={d['D1s_D2']:.3f} M1s,отн={d['M1s_rel']:.2f} "
          f"beta2b={d['beta2b']:.1f} Lz={d['Lz'] * 1e3:.0f} мм зазор={d['gap'] * 1e3:.2f} мм")
    g = res["geometry"]
    print(f"  охват: втулка {g['wrap_hub_deg']:.1f} град, покрышка {g['wrap_shroud_deg']:.1f} град; "
          f"beta втулки {g['beta_hub_min_deg']:.1f}..{g['beta_hub_max_deg']:.1f}")
    if "solid" in res:
        s = res["solid"]
        print(f"  STEP: {s['step']}  объём {s['volume_m3'] * 1e6:.0f} см3, масса Ti {s['mass_kg']['ti_6al_4v']:.1f} кг")
    if "mesh" in res:
        q = res["mesh"]["quality"]
        print(f"  SU2: {res['files']['passage_su2']}  ячеек {q['n_cells']}, узлов {q['n_nodes']}, "
              f"minSICN {q['minSICN']:.3f}, min gamma {q['min_gamma']:.3f}, геометрия {res['mesh']['t_geometry_s']} с, "
              f"сетка {res['mesh']['t_mesh_s']} с")
    for w in res["warnings"]:
        print("  ВНИМАНИЕ:", w)
    print("  JSON:", res["files"]["json"])
    if a.record:
        print("  итог дописан в", record(res))
    return res


if __name__ == "__main__":
    main()
