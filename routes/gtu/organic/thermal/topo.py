"""Оптимизация топологии коллектора на минимум потерь давления при равномерной раздаче (2D).

Метод — плотностный, Borrvall & Petersson (2003, Int. J. Numer. Meth. Fluids 41:77–107): течение Стокса —
Бринкмана, μ∇²u − α(ρ)u − ∇p = 0, ∇·u = 0, ρ = 1 — жидкость, ρ = 0 — стенка,
α(ρ) = α_max·q(1 − ρ)/(ρ + q); цель — рассеиваемая мощность Φ = ∫(μ|∇u|² + α|u|²)dA при доле жидкости ≤ V.
Равномерная раздача задаётся так же, как у Borrvall & Petersson: расходы на всех выходах равны и заданы
(профили скорости — граничное условие Дирихле), тогда задача самосопряжённая: dΦ/dρ = uᵀ(∂K/∂ρ)u.
Затем честная проверка: та же форма с выходами на одинаковое давление (как у реальных форсунок) —
расходы считаются, неравномерность выводится. Сравнение — с обычным коллектором «гребёнка» той же доли жидкости.

Дискретизация — разнесённая сетка MAC, прямой решатель SuperLU; фильтр плотности (Bourdin 2001) радиуса r,
обновление — критерий оптимальности (OC) с ограничением шага; продолжение по q (0,01 → 0,1 → 1).
Навье — Стокс для проверки формы: итерации Пикара, конвекция против потока (оценка влияния Re).
Единицы безразмерные: ширина области 1, μ = 1, средняя скорость на входе 1.
"""
import json
import math
import time
import sys
from pathlib import Path

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]


class Brinkman2D:
    """разнесённая сетка nx × ny, шаг h; отверстия на границе: список dict(side, j0, j1, kind, value)"""

    def __init__(self, nx, ny, h, openings, mu=1.0):
        self.nx, self.ny, self.h, self.mu = nx, ny, h, mu
        self.op = openings
        self.Nu = (nx + 1) * ny
        self.Nv = nx * (ny + 1)
        self.Np = nx * ny
        self._bc()
        self._base()

    def iu(self, i, j):
        return i * self.ny + j

    def iv(self, i, j):
        return self.Nu + i * (self.ny + 1) + j

    def ip(self, i, j):
        return self.Nu + self.Nv + i * self.ny + j

    def _bc(self):
        """нормальные граничные грани: Дирихле (стенка 0, вход/выход по профилю) или свободные (выход по давлению)"""
        nx, ny = self.nx, self.ny
        self.dir_val = {}
        self.free_out = {}          # грань -> (ячейка внутри, знак нормали, p_out)
        bnd = []
        for j in range(ny):
            bnd += [("W", self.iu(0, j), j), ("E", self.iu(nx, j), j)]
        for i in range(nx):
            bnd += [("S", self.iv(i, 0), i), ("N", self.iv(i, ny), i)]
        for side, k, s in bnd:
            self.dir_val[k] = 0.0
        for o in self.op:
            side, a, b = o["side"], o["j0"], o["j1"]
            n = b - a
            for s in range(a, b):
                xi = (s - a + 0.5) / n
                prof = 6 * xi * (1 - xi)                     # параболический профиль, средняя 1
                if side == "W":
                    k = self.iu(0, s); cell = (0, s); sgn = -1
                elif side == "E":
                    k = self.iu(self.nx, s); cell = (self.nx - 1, s); sgn = +1
                elif side == "S":
                    k = self.iv(s, 0); cell = (s, 0); sgn = -1
                else:
                    k = self.iv(s, self.ny); cell = (s, self.ny - 1); sgn = +1
                if o["kind"] == "velocity":
                    # value — средняя скорость внутрь области (вход > 0, выход < 0); знак по нормали
                    self.dir_val[k] = -sgn * o["value"] * prof
                else:
                    self.dir_val.pop(k)
                    self.free_out[k] = (cell, sgn, o["value"], side)

    def _base(self):
        """вязкая часть и градиент — не зависят от ρ"""
        nx, ny, h, mu = self.nx, self.ny, self.h, self.mu
        R, C, V = [], [], []
        diag = np.zeros(self.Nu + self.Nv)
        # u-грани
        for i in range(nx + 1):
            for j in range(ny):
                k = self.iu(i, j)
                if k in self.dir_val:
                    continue
                for di in (-1, 1):
                    ii = i + di
                    if 0 <= ii <= nx:
                        R.append(k); C.append(self.iu(ii, j)); V.append(-mu); diag[k] += mu
                    # за границей по нормали у свободного выхода — Нейман (нет связи)
                for dj in (-1, 1):
                    jj = j + dj
                    if 0 <= jj < ny:
                        R.append(k); C.append(self.iu(i, jj)); V.append(-mu); diag[k] += mu
                    else:
                        diag[k] += 2 * mu                     # стенка: зеркальный узел
        for i in range(nx):
            for j in range(ny + 1):
                k = self.iv(i, j)
                if k in self.dir_val:
                    continue
                for dj in (-1, 1):
                    jj = j + dj
                    if 0 <= jj <= ny:
                        R.append(k); C.append(self.iv(i, jj)); V.append(-mu); diag[k] += mu
                for di in (-1, 1):
                    ii = i + di
                    if 0 <= ii < nx:
                        R.append(k); C.append(self.iv(ii, j)); V.append(-mu); diag[k] += mu
                    else:
                        diag[k] += 2 * mu
        nvel = self.Nu + self.Nv
        R += list(range(nvel)); C += list(range(nvel)); V += list(diag)
        self.Kvis = sp.coo_matrix((V, (R, C)), shape=(nvel, nvel)).tocsr()
        # градиент: строки — грани, столбцы — ячейки; G[f, R] = +h, G[f, L] = −h
        R, C, V = [], [], []
        for i in range(nx + 1):
            for j in range(ny):
                k = self.iu(i, j)
                if i - 1 >= 0:
                    R.append(k); C.append((i - 1) * ny + j); V.append(-h)
                if i < nx:
                    R.append(k); C.append(i * ny + j); V.append(h)
        for i in range(nx):
            for j in range(ny + 1):
                k = self.iv(i, j)
                if j - 1 >= 0:
                    R.append(k); C.append(i * ny + j - 1); V.append(-h)
                if j < ny:
                    R.append(k); C.append(i * ny + j); V.append(h)
        self.G = sp.coo_matrix((V, (R, C)), shape=(nvel, self.Np)).tocsr()
        # усреднение α с ячеек на грани (внутренние — среднее двух; граничные свободные — одна ячейка)
        R, C, V = [], [], []
        for i in range(nx + 1):
            for j in range(ny):
                k = self.iu(i, j)
                cells = [c for c in ((i - 1) * ny + j if i >= 1 else None, i * ny + j if i < nx else None) if c is not None]
                for c in cells:
                    R.append(k); C.append(c); V.append(1.0 / len(cells))
        for i in range(nx):
            for j in range(ny + 1):
                k = self.iv(i, j)
                cells = [c for c in (i * ny + j - 1 if j >= 1 else None, i * ny + j if j < ny else None) if c is not None]
                for c in cells:
                    R.append(k); C.append(c); V.append(1.0 / len(cells))
        self.Pf = sp.coo_matrix((V, (R, C)), shape=(nvel, self.Np)).tocsr()
        nvel = self.Nu + self.Nv
        self.dir_idx = np.array(sorted(self.dir_val), np.int64)
        self.dir_vec = np.array([self.dir_val[k] for k in self.dir_idx])
        mask = np.ones(nvel, bool)
        mask[self.dir_idx] = False
        self.free_idx = np.where(mask)[0]
        self.pout_rhs = np.zeros(nvel)
        for k, (cell, sgn, pout, side) in self.free_out.items():
            # G p = h(p_R − p_L); у выхода справа/сверху p_R = p_out, слева/снизу p_L = p_out
            self.pout_rhs[k] = -sgn * h * pout

    def solve(self, alpha_cells, conv=None):
        """alpha_cells — α по ячейкам; conv — матрица конвекции (Пикар) или None. Возвращает u (все грани), p."""
        h = self.h
        alpha_f = self.Pf @ alpha_cells
        K = self.Kvis + sp.diags(alpha_f * h * h)
        if conv is not None:
            K = K + conv
        ufull = np.zeros(self.Nu + self.Nv)
        ufull[self.dir_idx] = self.dir_vec
        fi = self.free_idx
        Kff = K[fi][:, fi]
        rhs_u = -(K[fi][:, self.dir_idx] @ self.dir_vec) + self.pout_rhs[fi]
        Gf = self.G[fi]
        rhs_p = -(self.G[self.dir_idx].T @ self.dir_vec)       # Gᵀu = 0 → Gfᵀ u_f = −Gdᵀ u_d
        has_free_out = len(self.free_out) > 0
        if has_free_out:
            M = sp.bmat([[Kff, Gf], [Gf.T, None]], format="csc")
            rhs = np.concatenate([rhs_u, rhs_p])
        else:
            # давление до константы: закрепить ячейку 0
            GT = Gf.T.tolil(); GT[0, :] = 0
            Gc = Gf.tolil(); Gc[:, 0] = 0
            Z = sp.lil_matrix((self.Np, self.Np)); Z[0, 0] = 1.0
            M = sp.bmat([[Kff, Gc.tocsr()], [GT.tocsr(), Z.tocsr()]], format="csc")
            rhs = np.concatenate([rhs_u, rhs_p]); rhs[len(fi)] = 0.0
        x = spla.spsolve(M, rhs)
        ufull[fi] = x[:len(fi)]
        p = x[len(fi):]
        self.K_last = K
        return ufull, p

    def dissipation(self, u, alpha_cells):
        """Φ = uᵀKu (вязкая + Бринкман) на единицу глубины"""
        K = self.Kvis + sp.diags((self.Pf @ alpha_cells) * self.h ** 2)
        return float(u @ (K @ u))

    def dissipation_parts(self, u, alpha_cells):
        Kv = float(u @ (self.Kvis @ u))
        Ka = float(u @ (((self.Pf @ alpha_cells) * self.h ** 2) * u))
        return Kv, Ka

    def opening_flows(self, u):
        out = []
        for o in self.op:
            side, a, b = o["side"], o["j0"], o["j1"]
            q = 0.0
            for s in range(a, b):
                if side == "W":
                    q += u[self.iu(0, s)]
                elif side == "E":
                    q -= u[self.iu(self.nx, s)]
                elif side == "S":
                    q += u[self.iv(s, 0)]
                else:
                    q -= u[self.iv(s, self.ny)]
            out.append(q * self.h)       # > 0 — внутрь области
        return out

    def convection(self, u, rho_f=1.0):
        """матрица конвекции для итераций Пикара: противопоточная, консервативная по граням MAC (оценка)"""
        nx, ny, h = self.nx, self.ny, self.h
        nvel = self.Nu + self.Nv
        R, C, V = [], [], []
        U = u[:self.Nu].reshape(nx + 1, ny)
        Vv = u[self.Nu:].reshape(nx, ny + 1)

        def add(k, kn, F):
            # поток F наружу через грань контрольного объёма: противопоток
            if F >= 0:
                R.append(k); C.append(k); V.append(rho_f * F * h)
            elif kn is not None:
                R.append(k); C.append(kn); V.append(rho_f * F * h)
        for i in range(nx + 1):
            for j in range(ny):
                k = self.iu(i, j)
                if k in self.dir_val:
                    continue
                # восток/запад: скорость в центре ячейки
                if i < nx:
                    add(k, self.iu(i + 1, j), 0.5 * (U[i, j] + U[i + 1, j]))
                if i > 0:
                    add(k, self.iu(i - 1, j), -0.5 * (U[i, j] + U[i - 1, j]))
                ivn = (lambda ii, jj: Vv[ii, jj] if (0 <= ii < nx) else 0.0)
                Fn = 0.5 * (ivn(i - 1, j + 1) + ivn(i, j + 1))
                Fs = 0.5 * (ivn(i - 1, j) + ivn(i, j))
                add(k, self.iu(i, j + 1) if j + 1 < ny else None, Fn)
                add(k, self.iu(i, j - 1) if j - 1 >= 0 else None, -Fs)
        for i in range(nx):
            for j in range(ny + 1):
                k = self.iv(i, j)
                if k in self.dir_val:
                    continue
                if j < ny:
                    add(k, self.iv(i, j + 1), 0.5 * (Vv[i, j] + Vv[i, j + 1]))
                if j > 0:
                    add(k, self.iv(i, j - 1), -0.5 * (Vv[i, j] + Vv[i, j - 1]))
                iun = (lambda ii, jj: U[ii, jj] if (0 <= jj < ny) else 0.0)
                Fe = 0.5 * (iun(i + 1, j - 1) + iun(i + 1, j))
                Fw = 0.5 * (iun(i, j - 1) + iun(i, j))
                add(k, self.iv(i + 1, j) if i + 1 < nx else None, Fe)
                add(k, self.iv(i - 1, j) if i - 1 >= 0 else None, -Fw)
        return sp.coo_matrix((V, (R, C)), shape=(nvel, nvel)).tocsr()

    def solve_ns(self, alpha_cells, Re, iters=60, tol=1e-6, relax=0.7):
        """Навье — Стокс (ρ = Re при μ = 1, ширина 1, скорость 1), Пикар с релаксацией"""
        u, p = self.solve(alpha_cells)
        hist = []
        for it in range(iters):
            Cm = self.convection(u, rho_f=Re)
            un, pn = self.solve(alpha_cells, conv=Cm)
            d = np.linalg.norm(un - u) / max(np.linalg.norm(un), 1e-30)
            u = relax * un + (1 - relax) * u
            p = pn
            hist.append(d)
            if d < tol:
                break
        return u, p, dict(iterations=len(hist), last_change=float(hist[-1]))


def alpha_of(rho, amax, q):
    return amax * q * (1 - rho) / (rho + q)


def dalpha_of(rho, amax, q):
    return -amax * q * (1 + q) / (rho + q) ** 2


def density_filter(nx, ny, r):
    R, C, V = [], [], []
    ri = int(math.ceil(r))
    for i in range(nx):
        for j in range(ny):
            for di in range(-ri, ri + 1):
                for dj in range(-ri, ri + 1):
                    ii, jj = i + di, j + dj
                    if 0 <= ii < nx and 0 <= jj < ny:
                        w = r - math.hypot(di, dj)
                        if w > 0:
                            R.append(i * ny + j); C.append(ii * ny + jj); V.append(w)
    H = sp.coo_matrix((V, (R, C)), shape=(nx * ny, nx * ny)).tocsr()
    Hs = np.asarray(H.sum(1)).ravel()
    return H, Hs


def setup_case(case, nx):
    """case: 'comb6' — вход слева посередине, 6 выходов справа; 'fan8' — вход снизу, 8 выходов сверху и по бокам"""
    if case == "comb6":
        Lx, Ly = 1.0, 1.5
        ny = int(round(nx * Ly / Lx))
        h = Lx / nx
        w_in = 0.2
        n_out, w_out = 6, 0.08
        jin0 = int(round((Ly / 2 - w_in / 2) / h)); jin1 = int(round((Ly / 2 + w_in / 2) / h))
        ops = [dict(side="W", j0=jin0, j1=jin1, kind="velocity", value=1.0, name="вход")]
        Qin = (jin1 - jin0) * h * 1.0
        nw = int(round(w_out / h))
        for k in range(n_out):
            yc = Ly * (k + 0.5) / n_out
            j0 = int(round(yc / h - nw / 2))
            ops.append(dict(side="E", j0=j0, j1=j0 + nw, kind="velocity", value=-(Qin / n_out) / (nw * h),
                            name=f"выход {k + 1}"))
        return dict(nx=nx, ny=ny, h=h, Lx=Lx, Ly=Ly, ops=ops, Qin=Qin, n_out=n_out)
    if case == "fan8":
        Lx, Ly = 1.5, 1.0
        h = Lx / nx
        ny = int(round(Ly / h))
        w_in = 0.2
        iin0 = int(round((Lx / 2 - w_in / 2) / h)); iin1 = int(round((Lx / 2 + w_in / 2) / h))
        ops = [dict(side="S", j0=iin0, j1=iin1, kind="velocity", value=1.0, name="вход")]
        Qin = (iin1 - iin0) * h
        nw = int(round(0.07 / h))
        outs = []
        for k in range(4):                                  # сверху 4
            xc = Lx * (k + 0.5) / 4
            outs.append(("N", int(round(xc / h - nw / 2))))
        for side in ("W", "E"):                             # по бокам по 2
            for yc in (0.45, 0.8):
                outs.append((side, int(round(yc / h - nw / 2))))
        for k, (side, j0) in enumerate(outs):
            ops.append(dict(side=side, j0=j0, j1=j0 + nw, kind="velocity", value=-(Qin / len(outs)) / (nw * h),
                            name=f"выход {k + 1}"))
        return dict(nx=nx, ny=ny, h=h, Lx=Lx, Ly=Ly, ops=ops, Qin=Qin, n_out=len(outs))
    raise ValueError(case)


def with_pressure_outlets(cs):
    ops = []
    for o in cs["ops"]:
        if o["name"].startswith("выход"):
            ops.append(dict(o, kind="pressure", value=0.0))
        else:
            ops.append(dict(o))
    return ops


def optimize(case="comb6", nx=80, vol=0.30, Da=1e-5, rmin=1.5, qs=(0.01, 0.1, 1.0), iters=(60, 40, 40),
             move=0.2, verbose=True):
    cs = setup_case(case, nx)
    nx, ny, h = cs["nx"], cs["ny"], cs["h"]
    S = Brinkman2D(nx, ny, h, cs["ops"])
    amax = 1.0 / Da                                         # μ/(Da·L²), L = 1
    H, Hs = density_filter(nx, ny, rmin)
    x = np.full(nx * ny, vol)
    hist = []
    t0 = time.time()
    for q, nit in zip(qs, iters):
        for it in range(nit):
            rho = (H @ x) / Hs
            al = alpha_of(rho, amax, q)
            u, p = S.solve(al)
            Phi = S.dissipation(u, al)
            # чувствительность: dΦ/dρ_c = Σ_f u_f² h² dα_f/dρ_c (самосопряжённая задача)
            dal = dalpha_of(rho, amax, q)
            dPhi_drho = (S.Pf.T @ (u * u * h * h)) * dal
            dPhi_dx = H.T @ (dPhi_drho / Hs)
            dV = H.T @ (np.ones_like(x) / Hs)
            # OC: x_new = x·sqrt(−dΦ/(λ dV)); бисекция по λ
            l1, l2 = 1e-12, 1e12
            g = np.maximum(-dPhi_dx, 1e-30)
            while l2 / l1 > 1 + 1e-6:
                lm = math.sqrt(l1 * l2)
                xn = np.clip(x * np.sqrt(g / (lm * dV)), np.maximum(0, x - move), np.minimum(1, x + move))
                if ((H @ xn) / Hs).mean() > vol:
                    l1 = lm
                else:
                    l2 = lm
            change = float(np.abs(xn - x).max())
            x = xn
            mnd = float(4 * (rho * (1 - rho)).mean())
            hist.append(dict(q=q, it=it, Phi=Phi, change=change, mnd=mnd))
            if verbose and it % 10 == 0:
                print(f"q={q} it={it} Phi={Phi:.5g} change={change:.3f} mnd={mnd:.3f} t={time.time() - t0:.0f}s", flush=True)
    rho = (H @ x) / Hs
    return dict(case=cs, S=S, rho=rho.reshape(nx, ny), x=x, hist=hist, amax=amax, q=qs[-1], time_s=time.time() - t0)


def mma_update(x, df, dg, g, low, upp, xold1, xold2, it, move=0.2, asyinit=0.2, asyinc=1.2, asydec=0.7):
    """MMA (Svanberg 1987) для одного ограничения g(x) ≤ 0, двойственная задача — бисекция по λ ≥ 0.
    df, dg — градиенты цели и ограничения; xmin = 0, xmax = 1."""
    if it < 3:
        low = x - asyinit
        upp = x + asyinit
    else:
        sgn = (x - xold1) * (xold1 - xold2)
        fac = np.where(sgn > 0, asyinc, np.where(sgn < 0, asydec, 1.0))
        low = x - fac * (xold1 - low)
        upp = x + fac * (upp - xold1)
        low = np.clip(low, x - 10, x - 0.01)
        upp = np.clip(upp, x + 0.01, x + 10)
    alpha = np.maximum.reduce([np.zeros_like(x), low + 0.1 * (x - low), x - move])
    beta = np.minimum.reduce([np.ones_like(x), upp - 0.1 * (upp - x), x + move])
    ux2, xl2 = (upp - x) ** 2, (x - low) ** 2
    p0 = ux2 * (np.maximum(df, 0) + 1e-3 * np.abs(df) + 1e-9)
    q0 = xl2 * (np.maximum(-df, 0) + 1e-3 * np.abs(df) + 1e-9)
    p1 = ux2 * np.maximum(dg, 0)
    q1 = xl2 * np.maximum(-dg, 0)
    r1 = g - np.sum(p1 / (upp - x) + q1 / (x - low))

    def xof(lam):
        P = p0 + lam * p1
        Q = q0 + lam * q1
        xs = (np.sqrt(P) * low + np.sqrt(Q) * upp) / (np.sqrt(P) + np.sqrt(Q))
        return np.clip(xs, alpha, beta)

    def gof(xx):
        return r1 + np.sum(p1 / (upp - xx) + q1 / (xx - low))
    lo, hi = 0.0, 1.0
    while gof(xof(hi)) > 0 and hi < 1e12:
        hi *= 10
    for _ in range(80):
        mid = 0.5 * (lo + hi)
        if gof(xof(mid)) > 0:
            lo = mid
        else:
            hi = mid
    return xof(hi), low, upp


def project(rt, beta, eta=0.5):
    """проекция Хевисайда (Wang, Lazarov, Sigmund 2011) и её производная"""
    if beta <= 0:
        return rt, np.ones_like(rt)
    den = math.tanh(beta * eta) + math.tanh(beta * (1 - eta))
    rb = (math.tanh(beta * eta) + np.tanh(beta * (rt - eta))) / den
    d = beta * (1 - np.tanh(beta * (rt - eta)) ** 2) / den
    return rb, d


STAGES = ((0.01, 0, 60), (0.1, 0, 40), (0.1, 2, 30), (0.1, 4, 30), (0.1, 8, 30), (0.1, 16, 30))


STAGES_WARM = ((0.1, 0, 40), (0.1, 2, 30), (0.1, 4, 30), (0.1, 8, 30), (0.1, 16, 40))


def optimize_uniform(case="comb6", nx=80, vol=0.30, Da=1e-5, rmin=1.5, stages=STAGES,
                     w=1.0, move=0.1, verbose=True, x_init=None):
    """выходы на одинаковое давление; цель J = Φ/Φ0 + w·N·Σ(Q_k/Q̄ − 1)² при доле жидкости ≤ vol;
    Φ — самосопряжённая часть, разброс расходов — сопряжённая задача с той же (симметричной) матрицей"""
    cs = setup_case(case, nx)
    nx, ny, h = cs["nx"], cs["ny"], cs["h"]
    S = Brinkman2D(nx, ny, h, with_pressure_outlets(cs))
    amax = 1.0 / Da
    H, Hs = density_filter(nx, ny, rmin)
    x = np.full(nx * ny, vol) if x_init is None else np.clip(np.asarray(x_init, float).copy(), 0, 1)
    nout = cs["n_out"]
    Qbar = cs["Qin"] / nout
    # функционалы расхода по выходам: Q_k = c_kᵀ u (u — все грани)
    cvec = []
    for o in S.op[1:]:
        c = np.zeros(S.Nu + S.Nv)
        for s in range(o["j0"], o["j1"]):
            if o["side"] == "E":
                c[S.iu(nx, s)] = h
            elif o["side"] == "W":
                c[S.iu(0, s)] = -h
            elif o["side"] == "N":
                c[S.iv(s, ny)] = h
            else:
                c[S.iv(s, 0)] = -h
        cvec.append(c)
    cvec = np.array(cvec)
    fi = S.free_idx
    low = upp = None
    xold1 = xold2 = x.copy()
    hist = []
    Phi0 = None
    t0 = time.time()
    k = 0
    for q, beta, nit in stages:
        low = upp = None
        for it in range(nit):
            rt = (H @ x) / Hs
            rho, drp = project(rt, beta)
            al = alpha_of(rho, amax, q)
            dal = dalpha_of(rho, amax, q) * drp
            alpha_f = S.Pf @ al
            K = S.Kvis + sp.diags(alpha_f * h * h)
            Kff = K[fi][:, fi]
            Gf = S.G[fi]
            M = sp.bmat([[Kff, Gf], [Gf.T, None]], format="csc")
            rhs_u = -(K[fi][:, S.dir_idx] @ S.dir_vec) + S.pout_rhs[fi]
            rhs_p = -(S.G[S.dir_idx].T @ S.dir_vec)
            lu = spla.splu(M, permc_spec="COLAMD")
            xs = lu.solve(np.concatenate([rhs_u, rhs_p]))
            u = np.zeros(S.Nu + S.Nv)
            u[S.dir_idx] = S.dir_vec
            u[fi] = xs[:len(fi)]
            Phi = float(u @ (K @ u))
            if Phi0 is None:
                Phi0 = Phi
            Q = cvec @ u
            dev = Q / Qbar - 1
            J2 = w * nout * float((dev ** 2).sum())
            J = Phi / Phi0 + J2
            # чувствительности
            dPhi = (S.Pf.T @ (u * u * h * h)) * dal / Phi0
            gJ2 = np.zeros(len(xs))
            gu = (2 * w * nout * dev / Qbar) @ cvec          # ∂J2/∂u (все грани)
            gJ2[:len(fi)] = gu[fi]
            lam = lu.solve(gJ2, trans="T")
            lam_u = np.zeros(S.Nu + S.Nv)
            lam_u[fi] = lam[:len(fi)]
            dJ2 = -(S.Pf.T @ (lam_u * u * h * h)) * dal
            dJ_drho = dPhi + dJ2
            dJ_dx = H.T @ (dJ_drho / Hs)
            gV = float(rho.mean() - vol)
            dV = (H.T @ (drp / Hs)) / x.size
            xn, low, upp = mma_update(x, dJ_dx, dV, gV, low if low is not None else x - 0.2,
                                      upp if upp is not None else x + 0.2, xold1, xold2,
                                      it if low is not None else 0, move=move)
            xold2, xold1 = xold1, x.copy()
            change = float(np.abs(xn - x).max())
            x = xn
            mnd = float(4 * (rho * (1 - rho)).mean())
            hist.append(dict(q=q, beta=beta, it=it, J=J, Phi=Phi, J2=J2, max_dev=float(np.abs(dev).max()),
                             change=change, mnd=mnd))
            if verbose and it % 10 == 0:
                print(f"[U] q={q} b={beta} it={it} J={J:.4g} Phi={Phi:.4g} maxdev={np.abs(dev).max():.3f} mnd={mnd:.3f} "
                      f"t={time.time() - t0:.0f}s", flush=True)
            k += 1
    q, beta, _ = stages[-1]
    rho, _ = project((H @ x) / Hs, beta)
    return dict(case=cs, S=S, rho=rho.reshape(nx, ny), x=x, hist=hist, amax=amax, q=q, beta=beta,
                time_s=time.time() - t0, w=w, Phi0=Phi0)


def comb_baseline(cs, vol):
    """обычный коллектор той же доли жидкости: подводящая труба от входа, прямой вертикальный коллектор
    у правой стороны и короткие отводы к выходам; ширины подобраны под долю vol (оценка «как делают»)"""
    nx, ny, h = cs["nx"], cs["ny"], cs["h"]
    op_in = cs["ops"][0]
    outs = cs["ops"][1:]
    target = vol * nx * ny

    def make(wh):
        m = np.zeros((nx, ny), bool)
        j0, j1 = op_in["j0"], op_in["j1"]
        xh = int(round(0.70 * nx))
        wcol = max(2, int(round(wh / h)))
        m[:xh + wcol, j0:j1] = True                            # подводящий канал
        js = [o["j0"] for o in outs] + [o["j1"] for o in outs]
        m[xh:xh + wcol, min(js):max(js)] = True                # вертикальный коллектор
        for o in outs:
            m[xh:, o["j0"]:o["j1"]] = True                     # отводы
        return m
    best = None
    for wh in np.linspace(0.02, 0.4, 200):
        m = make(wh)
        d = abs(m.sum() - target)
        if best is None or d < best[0]:
            best = (d, wh, m)
    return best[2], best[1]


def _opening_center(cs, o):
    h, nx, ny = cs["h"], cs["nx"], cs["ny"]
    m = 0.5 * (o["j0"] + o["j1"]) * h
    return {"W": (0.0, m), "E": (nx * h, m), "S": (m, 0.0), "N": (m, ny * h)}[o["side"]]


def spoke_baseline(cs, vol):
    """обычный коллектор «ступица и спицы»: прямой подвод от входа к центру и прямые трубы равной ширины к выходам"""
    nx, ny, h = cs["nx"], cs["ny"], cs["h"]
    X, Y = np.meshgrid((np.arange(nx) + 0.5) * h, (np.arange(ny) + 0.5) * h, indexing="ij")
    hub = (0.5 * nx * h, 0.45 * ny * h)
    cin = _opening_center(cs, cs["ops"][0])
    win = (cs["ops"][0]["j1"] - cs["ops"][0]["j0"]) * h

    def seg_mask(P, Q, w):
        P, Q = np.array(P), np.array(Q)
        d = Q - P
        t = np.clip(((X - P[0]) * d[0] + (Y - P[1]) * d[1]) / (d @ d), 0, 1)
        return np.hypot(X - P[0] - t * d[0], Y - P[1] - t * d[1]) <= w / 2
    target = vol * nx * ny
    best = None
    for w in np.linspace(0.02, 0.3, 120):
        m = seg_mask(cin, hub, win)
        for o in cs["ops"][1:]:
            m |= seg_mask(hub, _opening_center(cs, o), w)
        d = abs(m.sum() - target)
        if best is None or d < best[0]:
            best = (d, w, m)
    return best[2], best[1]


def evaluate(cs, rho_cells, amax, q, label, Re_list=(0.0,)):
    """форма → (1) потери при заданных равных расходах; (2) выходы на одинаковое давление: расходы и потери"""
    nx, ny, h = cs["nx"], cs["ny"], cs["h"]
    al = alpha_of(np.clip(rho_cells.ravel(), 0, 1), amax, q)
    S1 = Brinkman2D(nx, ny, h, cs["ops"])
    u1, p1 = S1.solve(al)
    Phi1 = S1.dissipation(u1, al)
    # давление у каждого выхода при заданных равных расходах (среднее по ячейкам у отверстия)
    P1 = p1.reshape(nx, ny)
    pk = []
    for o in cs["ops"][1:]:
        sl = range(o["j0"], o["j1"])
        if o["side"] == "E":
            pk.append(float(np.mean([P1[nx - 1, j] for j in sl])))
        elif o["side"] == "W":
            pk.append(float(np.mean([P1[0, j] for j in sl])))
        elif o["side"] == "N":
            pk.append(float(np.mean([P1[i, ny - 1] for i in sl])))
        else:
            pk.append(float(np.mean([P1[i, 0] for i in sl])))
    pk = np.array(pk)
    Qk = cs["Qin"] / cs["n_out"]
    throttle = float(Qk * np.sum(pk - pk.min()))            # идеальные дроссели до общего давления выходов
    S2 = Brinkman2D(nx, ny, h, with_pressure_outlets(cs))
    out = dict(label=label, fluid_fraction=float(rho_cells.mean()), Phi_equal_split=Phi1,
               outlet_pressures_equal_split=[float(v) for v in pk], Phi_throttles=throttle,
               Phi_balanced=Phi1 + throttle,
               note_balanced="Φ_balanced = Φ при равных расходах + мощность идеальных дросселей, выравнивающих "
                             "давления выходов: цена равной раздачи на одинаковое противодавление")
    rows = []
    for Re in Re_list:
        if Re == 0:
            u2, p2 = S2.solve(al)
            info = None
        else:
            u2, p2, info = S2.solve_ns(al, Re)
        Q = S2.opening_flows(u2)
        Qout = -np.array(Q[1:])
        Qin = Q[0]
        # потери: работа давления на входе = Φ (выходы при p = 0); давление на входе — среднее по ячейкам входа
        o = cs["ops"][0]
        pin = float(np.mean([p2[i * ny + j] for i, j in ((0, s) for s in range(o["j0"], o["j1"]))]))
        Phi2 = S2.dissipation(u2, al)
        rows.append(dict(Re=Re, Q_in=float(Qin), Q_out=[float(v) for v in Qout],
                         nonuniformity_max_min=float(Qout.max() / Qout.min() - 1) if Qout.min() > 0 else None,
                         cv=float(Qout.std() / Qout.mean()), p_in=pin, Phi=Phi2, ns=info))
    out["pressure_outlets"] = rows
    out["_u1"], out["_u2"] = u1, u2
    return out


def contour_stl(path, rho, h, thickness=0.1, level=0.5):
    """2D-форма (изолиния ρ = 0,5) → тело коллектора: стенка вокруг каналов не строится; выдавливание
    области жидкости на толщину (вид «негатив каналов»); маршевые квадраты через skimage"""
    from skimage import measure
    from tpms_geom import write_stl
    nx, ny = rho.shape
    pad = np.pad(rho, 1, constant_values=0.0)
    # объём: выдавливаем пиксели жидкости — простое и честное тело для печати каналов как «негатива»
    vol = np.repeat(pad[:, :, None], 6, axis=2)
    vol[:, :, 0] = 0; vol[:, :, -1] = 0
    v, f, _, _ = measure.marching_cubes(vol, level=level, spacing=(h, h, thickness / 4))
    return write_stl(path, v, f, solid="manifold")


def plot(cs, evals, path, title):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    nx, ny, h = cs["nx"], cs["ny"], cs["h"]
    fig, axs = plt.subplots(1, len(evals), figsize=(4.6 * len(evals), 4.6 * cs["Ly"] / cs["Lx"] + 0.9))
    for ax, ev in zip(np.atleast_1d(axs), evals):
        rho = ev["_rho"]
        u = ev["_u2"]
        S = ev["_S"]
        U = u[:S.Nu].reshape(nx + 1, ny)
        V = u[S.Nu:].reshape(nx, ny + 1)
        uc = 0.5 * (U[1:] + U[:-1]); vc = 0.5 * (V[:, 1:] + V[:, :-1])
        spd = np.hypot(uc, vc)
        ext = [0, cs["Lx"], 0, cs["Ly"]]
        ax.imshow(rho.T, origin="lower", extent=ext, cmap="Greys_r", vmin=0, vmax=1)
        sm = np.ma.masked_where(rho < 0.5, spd)
        ax.imshow(sm.T, origin="lower", extent=ext, cmap="inferno", alpha=0.9)
        x = (np.arange(nx) + 0.5) * h; y = (np.arange(ny) + 0.5) * h
        um = np.where(rho > 0.5, uc, np.nan); vm = np.where(rho > 0.5, vc, np.nan)
        ax.streamplot(x, y, um.T, vm.T, color="#9fe0ff", linewidth=0.5, density=1.6, arrowsize=0.6)
        r = ev["pressure_outlets"][0]
        nu = r["nonuniformity_max_min"]
        ax.set_title(f"{ev['label']}\nразброс расходов при p_вых = const: {100 * nu:.1f} %\n"
                     f"Φ с дросселями до равных расходов: {ev['Phi_balanced']:.0f}", fontsize=8.5)
        ax.set_xticks([]); ax.set_yticks([])
    fig.suptitle(title, fontsize=10)
    fig.tight_layout()
    fig.savefig(path, dpi=130)
    plt.close(fig)


def murray2d_check(rho, h, level=0.5):
    """проверка «закона Мюррея» плоских каналов: для течения Пуазейля в щели Q ∝ w³, оптимум Мюррея даёт
    w₀² = w₁² + w₂² (показатель 2); ширины ветвей — по скелету и расстоянию до стенки (skimage)"""
    from skimage.morphology import skeletonize
    from scipy import ndimage
    m = rho > level
    sk = skeletonize(m)
    dist = ndimage.distance_transform_edt(m) * h
    # узлы ветвления: точки скелета с ≥ 3 соседями
    k = ndimage.convolve(sk.astype(int), np.ones((3, 3), int), mode="constant") - sk
    br = sk & (k >= 3)
    lab, n = ndimage.label(sk & ~ndimage.binary_dilation(br, iterations=2))
    widths = []
    for i in range(1, n + 1):
        sel = lab == i
        if sel.sum() < 6:
            continue
        widths.append(float(2 * np.median(dist[sel])))
    return dict(n_branch_nodes=int(ndimage.label(br)[1]), n_segments=len(widths),
                widths=sorted(widths, reverse=True)[:20])


def run_all(out_dir, nx=90, vol=0.30, cases=("comb6", "fan8")):
    from tpms_geom import write_stl
    out_dir = Path(out_dir)
    res = dict(title="Оптимизация топологии коллектора: минимум потерь при равномерной раздаче (2D, Стокс — Бринкман)",
               method=__doc__.strip(),
               software=f"python {sys.version.split()[0]}, numpy {np.__version__}, scipy (SuperLU), scikit-image (изолинии)",
               date="2026-11-14",
               status="демонстрация метода в 2D при малом Re (безразмерно); реальный коллектор KV-3 — турбулентный "
                      "(Re ~ 10⁵–10⁶): нужна трёхмерная оптимизация с CFD (задача AWS)",
               inputs=dict(note="безразмерная постановка; связь с KV-3 — только по числу выходов (6 и 8) и по смыслу узла"),
               parameters=dict(nx=nx, vol=vol, Da=1e-5, rmin_cells=1.5, stages=STAGES, w_uniform=1.0), cases={})
    for case in cases:
        t0 = time.time()
        rB = optimize(case, nx=nx, vol=vol, qs=(0.01, 0.1), iters=(60, 60), verbose=False)   # классика: расходы заданы
        rU = optimize_uniform(case, nx=nx, vol=vol, verbose=False)         # выходы на давлении + равномерность
        rW = optimize_uniform(case, nx=nx, vol=vol, verbose=False, x_init=rB["x"], stages=STAGES_WARM)  # тёплый старт
        cs = rU["case"]
        amax = rU["amax"]
        # оценка всех форм при одинаковой α (q = 0,1, как в конце rU) — честное сравнение
        q = 0.1
        evs = []
        for label, rho in (("Борвалл — Петерсон: расходы заданы", rB["rho"]),
                           ("равномерность в цели (выходы на давлении)", rU["rho"]),
                           ("Б — П, затем равномерность (тёплый старт)", rW["rho"])):
            ev = evaluate(cs, rho, amax, q, label, Re_list=(0.0, 20.0, 50.0))
            ev["_rho"] = rho
            ev["_S"] = Brinkman2D(cs["nx"], cs["ny"], cs["h"], with_pressure_outlets(cs))
            evs.append(ev)
        if case == "comb6":
            mb, wh = comb_baseline(cs, float(rU["rho"].mean()))
            blabel = "обычная гребёнка (та же доля жидкости)"
        else:
            mb, wh = spoke_baseline(cs, float(rU["rho"].mean()))
            blabel = "ступица и спицы (та же доля жидкости)"
        evb = evaluate(cs, mb.astype(float), amax, q, blabel, Re_list=(0.0, 20.0, 50.0))
        evb["_rho"] = mb.astype(float)
        evb["_S"] = Brinkman2D(cs["nx"], cs["ny"], cs["h"], with_pressure_outlets(cs))
        evs.append(evb)
        plot(cs, [evs[3], evs[0], evs[1], evs[2]], out_dir / f"topo_{case}.png",
             f"Коллектор: 1 вход → {cs['n_out']} выходов, доля жидкости {vol:.2f} (2D, Стокс — Бринкман, безразмерно)")
        # STL — форма тёплого старта: равномерна без дросселей (у формы Б — П нужны шайбы на выходах)
        bi = 2
        best_rho = rW["rho"]
        np.savez_compressed(out_dir / f"topo_{case}_rho.npz", bp=rB["rho"], uniform=rU["rho"], warm=rW["rho"],
                            conventional=mb.astype(float), h=cs["h"])
        stl = contour_stl(out_dir / f"topo_{case}.stl", best_rho, cs["h"] * 100.0, thickness=4.0)
        stl["design"] = evs[bi]["label"]
        stl["path"] = str(Path(stl["path"]).relative_to(ROOT)) if str(stl["path"]).startswith(str(ROOT)) else stl["path"]
        stl["units"] = "условные мм: ширина области 100, толщина 4 (выдавливание формы каналов)"
        rows = []
        for ev in evs:
            rows.append({k: v for k, v in ev.items() if not k.startswith("_")})
        base_phi = evb["pressure_outlets"][0]["Phi"]
        evw = evs[2]
        res["cases"][case] = dict(
            setup=dict(nx=cs["nx"], ny=cs["ny"], Lx=cs["Lx"], Ly=cs["Ly"], n_out=cs["n_out"], Q_in=cs["Qin"]),
            designs=rows,
            comb_header_width=float(wh),
            history_uniform=[dict(step=i, **{k: (float(v) if isinstance(v, (float, np.floating)) else v) for k, v in hh.items()})
                             for i, hh in enumerate(rU["hist"]) if i % 10 == 0 or i == len(rU["hist"]) - 1],
            dissipation_ratio_uniform_vs_comb=float(evs[1]["pressure_outlets"][0]["Phi"] / base_phi),
            nonuniformity=dict(comb=evb["pressure_outlets"][0]["nonuniformity_max_min"],
                               bp=evs[0]["pressure_outlets"][0]["nonuniformity_max_min"],
                               uniform=evs[1]["pressure_outlets"][0]["nonuniformity_max_min"]),
            murray2d=murray2d_check(best_rho, cs["h"]),
            stl=stl, time_s=round(time.time() - t0, 1))
        res["cases"][case]["balanced_ratio_uniform_vs_conventional"] = float(evs[1]["Phi_balanced"] / evb["Phi_balanced"])
        res["cases"][case]["balanced_ratio_bp_vs_conventional"] = float(evs[0]["Phi_balanced"] / evb["Phi_balanced"])
        res["cases"][case]["balanced_ratio_warm_vs_conventional"] = float(evw["Phi_balanced"] / evb["Phi_balanced"])
        res["cases"][case]["nonuniformity"]["warm"] = evw["pressure_outlets"][0]["nonuniformity_max_min"]
        print(case, "готово", round(time.time() - t0), res["cases"][case]["nonuniformity"],
              res["cases"][case]["balanced_ratio_uniform_vs_conventional"], flush=True)
        Path(out_dir, "topo_results.partial.json").write_text(json.dumps(res, ensure_ascii=False, indent=1, default=float))
    return res


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--nx", type=int, default=90)
    ap.add_argument("--vol", type=float, default=0.30)
    ap.add_argument("--cases", nargs="+", default=["comb6", "fan8"])
    ap.add_argument("--out", default=str(HERE))
    a = ap.parse_args()
    res = run_all(a.out, a.nx, a.vol, a.cases)
    Path(a.out, "topo_results.json").write_text(json.dumps(res, ensure_ascii=False, indent=1, default=float))
    print("записано")
