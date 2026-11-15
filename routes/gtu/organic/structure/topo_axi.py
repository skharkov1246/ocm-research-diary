#!/usr/bin/env python3
"""Оптимизация топологии осесимметричного тела вращения (диск, колесо) — свой КЭ-код на numpy/scipy.

Постановка (README узла, раздел «Метод»):
  * сетка — прямоугольная решётка Q4-элементов в меридиональной плоскости (r, z), активны только элементы
    внутри огибающей; осесимметричные деформации [εr, εz, εθ, γrz], интегрирование 2×2 по Гауссу, dV = 2πr dA;
  * нагрузки — вращение (объёмная сила ρω²r, зависит от плотности), тепловая (σ = D(ε − α(T)(T − Tref))),
    постоянные силы: давление лопаток на ободе, масса лопаток колеса на поверхности ступицы — растут как ω²;
  * свойства от температуры: E(T), α(T), предел текучести, прочности и длительной прочности — из
    routes/gtu/platform/fea/materials.json через disk.load_material (одна точка правды с CalculiX-узлом);
  * плотностный метод: x — переменные, x~ = H·x (линейный фильтр плотности радиуса R), x̄ = проекция Хевисайда
    (β, η = 0,5); интерполяции жёсткости, массы и «ослабления» напряжений задаются группой элементов:
      - «solid» (сплошной металл или пустота, вариант A): E = Emin + (1 − Emin)·x̄/(1 + q(1 − x̄)) (RAMP —
        конечное отношение массы к жёсткости при x̄ → 0, нужно при массовой нагрузке: Bruyneel & Duysinx 2005),
        масса ∝ x̄, напряжение ослабляется множителем x̄^0,5 (qp-подход: Le et al. 2010, Bruggi 2008);
      - «lattice» (решётка TPMS переменной плотности, вариант B): ρ* = ρmin + (1 − ρmin)·x̄, E*/E и σ*/σ —
        таблица гомогенизации гироида (lattice.py), напряжение решётки сравнивается с её прочностью;
      - «passive» — сплошной металл, не меняется (обшивка, обод, посадка на вал);
  * ограничения — по напряжению: p-норма загрузки u = σ_экв/σ_доп(T) по каждому расчётному случаю
    (P = 10) с адаптивной нормировкой (Le et al. 2010); производные — сопряжённым методом (одна факторизация
    матрицы жёсткости на все случаи и сопряжённые задачи);
  * цель — масса или полярный момент инерции (или их взвешенная сумма);
  * оптимизатор — MMA (Svanberg 1987): выпуклое разделимое приближение, подзадача решается через
    двойственную (m ≤ 4 множителей, L-BFGS-B), переменные y — штраф невыполнимости.

Проверки — test_structure.py: тонкий вращающийся диск и нагрев по аналитике, производные — конечными разностями.
"""
import math
import os
import sys

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
from scipy.optimize import minimize
from scipy.spatial import cKDTree

HERE = os.path.dirname(os.path.abspath(__file__))
FEA = os.path.normpath(os.path.join(HERE, "..", "..", "platform", "fea"))
if FEA not in sys.path:
    sys.path.insert(0, FEA)
import disk as DK  # noqa: E402  (routes/gtu/platform/fea/disk.py: материалы, Ларсон — Миллер, разбор ccx)

GP = (-1.0 / math.sqrt(3.0), 1.0 / math.sqrt(3.0))
MVEC = np.array([1.0, 1.0, 1.0, 0.0])


# =========================================================================== материал
def mat_tables(mat):
    """Функции свойств от температуры (векторно), МПа и Па. За краем таблицы — крайнее значение (как ccx)."""
    def tab(key, scale=1.0):
        a = np.asarray(mat[key], float)
        return lambda T: np.interp(np.asarray(T, float), a[:, 0], a[:, 1]) * scale
    out = dict(E=tab("E_GPa", 1e9), alpha=tab("alpha_mean_1e6", 1e-6), uts=tab("uts_MPa"), ys=tab("ys_MPa"),
               rho=float(mat["density"]), nu=float(mat["nu"]))
    lmp = mat.get("_lmp")
    if lmp is not None:
        def slm(T, t_h):
            T = np.asarray(T, float)
            P = (T + 273.15) * (lmp.C + math.log10(t_h))
            return 10.0 ** lmp.log_sigma_v(P)
        out["slm"] = slm
    return out


def D_unit(nu):
    c = 1.0 / ((1 + nu) * (1 - 2 * nu))
    return c * np.array([[1 - nu, nu, nu, 0], [nu, 1 - nu, nu, 0], [nu, nu, 1 - nu, 0],
                         [0, 0, 0, 0.5 * (1 - 2 * nu)]])


def von_mises4(s):
    """s: (..., 4) = [σr, σz, σθ, τrz]"""
    sr, sz, st, t = s[..., 0], s[..., 1], s[..., 2], s[..., 3]
    return np.sqrt(np.maximum(0.5 * ((sr - sz) ** 2 + (sz - st) ** 2 + (st - sr) ** 2) + 3 * t * t, 0.0))


def dvm_ds(s, vm):
    sr, sz, st, t = s[..., 0], s[..., 1], s[..., 2], s[..., 3]
    v = np.maximum(vm, 1e-12)[..., None]
    return np.stack([2 * sr - sz - st, 2 * sz - sr - st, 2 * st - sr - sz, 6 * t], axis=-1) / (2 * v)


# =========================================================================== сетка и КЭ
class AxiModel:
    """Осесимметричная модель на прямоугольной решётке nr×nz в [r0, r1]×[z0, z1].

    active: bool (nr, nz) — элементы внутри огибающей (остальные в расчёт не входят).
    group:  int (nr, nz) — 0 пассивный сплошной, 1 «solid» (A), 2 «lattice» (B).
    T_of_r: функция температуры металла от радиуса, °C (по оси не меняется — как в disk.py).
    """

    def __init__(self, r0, r1, z0, z1, nr, nz, active, group, mat, T_of_r, Tref=20.0):
        self.r0, self.r1, self.z0, self.z1, self.nr, self.nz = r0, r1, z0, z1, nr, nz
        self.dr, self.dz = (r1 - r0) / nr, (z1 - z0) / nz
        self.mat, self.mt = mat, mat_tables(mat)
        self.nu = self.mt["nu"]; self.D1 = D_unit(self.nu)
        self.Tref = Tref
        self.T_of_r = T_of_r
        I, J = np.nonzero(active)
        self.ei, self.ej = I, J                         # индексы активных элементов
        self.ne = len(I)
        self.group = group[I, J].astype(int)
        self.rc = r0 + (I + 0.5) * self.dr
        self.zc = z0 + (J + 0.5) * self.dz
        # узлы решётки (i, j) -> сквозной номер только для используемых
        nzp = nz + 1
        gn = np.stack([I * nzp + J, (I + 1) * nzp + J, (I + 1) * nzp + J + 1, I * nzp + J + 1], axis=1)
        used, inv = np.unique(gn, return_inverse=True)
        self.node_grid = used                          # номер узла решётки
        self.nn = len(used)
        self.enodes = inv.reshape(-1, 4)
        self.node_r = r0 + (used // nzp) * self.dr
        self.node_z = z0 + (used % nzp) * self.dz
        self.edof = np.empty((self.ne, 8), int)
        self.edof[:, 0::2] = 2 * self.enodes
        self.edof[:, 1::2] = 2 * self.enodes + 1
        self.ndof = 2 * self.nn
        # температура и свойства по элементу (в центре), εth — в точках Гаусса (по столбцу)
        self.Te = np.asarray(T_of_r(self.rc), float)
        self.Es = self.mt["E"](self.Te)
        self.col_mats()
        self.Ve = 2 * math.pi * self.rc * self.dr * self.dz
        self.rho = self.mt["rho"]
        ii, jj = np.meshgrid(np.arange(8), np.arange(8), indexing="ij")
        self.Ki = self.edof[:, ii].reshape(-1)
        self.Kj = self.edof[:, jj].reshape(-1)
        self.fixed = np.zeros(self.ndof, bool)
        self.f_fix0 = np.zeros(self.ndof)              # постоянные силы при расчётной частоте (растут как ω²)
        self.rpm0 = None

    def col_mats(self):
        """Матрицы элемента для E = 1 по столбцу решётки (зависят только от радиуса)."""
        nr, dr, dz, D1 = self.nr, self.dr, self.dz, self.D1
        Ke = np.zeros((nr, 8, 8)); fc = np.zeros((nr, 8)); fth = np.zeros((nr, 8)); Bc = np.zeros((nr, 4, 8))
        epsc = np.zeros(nr)
        detJ = dr * dz / 4
        for i in range(nr):
            ra = self.r0 + i * dr
            for xi in GP + (0.0,):
                for eta in (GP if xi != 0.0 else (0.0,)):
                    N = 0.25 * np.array([(1 - xi) * (1 - eta), (1 + xi) * (1 - eta), (1 + xi) * (1 + eta), (1 - xi) * (1 + eta)])
                    dNx = 0.25 * np.array([-(1 - eta), (1 - eta), (1 + eta), -(1 + eta)]) * 2 / dr
                    dNz = 0.25 * np.array([-(1 - xi), -(1 + xi), (1 + xi), (1 - xi)]) * 2 / dz
                    r = ra + (1 + xi) * dr / 2
                    B = np.zeros((4, 8))
                    B[0, 0::2] = dNx; B[1, 1::2] = dNz; B[2, 0::2] = N / r
                    B[3, 0::2] = dNz; B[3, 1::2] = dNx
                    T = float(self.T_of_r(np.array([r]))[0])
                    eth = float(self.mt["alpha"](T)) * (T - self.Tref)
                    if xi == 0.0:
                        Bc[i] = B; epsc[i] = eth
                        continue
                    w = 2 * math.pi * r * detJ
                    Ke[i] += B.T @ D1 @ B * w
                    fc[i, 0::2] += N * r * w
                    fth[i] += B.T @ D1 @ MVEC * eth * w
        self.Ke0, self.fc0, self.fth0, self.Bc0, self.epsc0 = Ke, fc, fth, Bc, epsc

    # ------------------------------------------------------------------ закрепления и нагрузки
    def fix(self, nodes, comp):
        self.fixed[2 * np.asarray(nodes, int) + comp] = True

    def nodes_where(self, cond):
        return np.nonzero(cond(self.node_r, self.node_z))[0]

    def add_force(self, nodes, fr=None, fz=None):
        """Сосредоточенные силы при расчётной частоте rpm0, Н (по всей окружности — осесимметрично)."""
        if fr is not None:
            np.add.at(self.f_fix0, 2 * np.asarray(nodes), fr)
        if fz is not None:
            np.add.at(self.f_fix0, 2 * np.asarray(nodes) + 1, fz)

    def add_edge_traction_r(self, r_face, zlo, zhi, p):
        """Растягивающее радиальное давление p (Па) на наружной грани r = r_face, zlo ≤ z ≤ zhi."""
        nds = self.nodes_where(lambda r, z: (np.abs(r - r_face) < 1e-9 * max(1, r_face) + 1e-12) &
                               (z >= zlo - 1e-12) & (z <= zhi + 1e-12))
        nds = nds[np.argsort(self.node_z[nds])]
        zs = self.node_z[nds]
        f = np.zeros(len(nds))
        for k in range(len(nds) - 1):
            L = zs[k + 1] - zs[k]
            F = p * 2 * math.pi * r_face * L
            f[k] += F / 2; f[k + 1] += F / 2
        self.add_force(nds, fr=f)
        return float(f.sum())

    def free_map(self):
        free = np.nonzero(~self.fixed)[0]
        m = -np.ones(self.ndof, int); m[free] = np.arange(len(free))
        ri, ci = m[self.Ki], m[self.Kj]
        keep = (ri >= 0) & (ci >= 0)
        self._free, self._ri, self._ci, self._keep = free, ri[keep], ci[keep], keep

    # ------------------------------------------------------------------ интерполяции по группам
    def interp(self, xb, prm):
        """η_E, η_m, η_σ и производные по x̄ (по элементу)."""
        g = self.group
        eE = np.ones(self.ne); em = np.ones(self.ne); es = np.ones(self.ne)
        dE = np.zeros(self.ne); dm = np.zeros(self.ne); ds = np.zeros(self.ne)
        s = g == 1
        if s.any():
            q, Emin = prm["q"], prm["Emin"]
            x = xb[s]
            den = 1 + q * (1 - x)
            eE[s] = Emin + (1 - Emin) * x / den
            dE[s] = (1 - Emin) * (1 + q) / den ** 2
            em[s] = x; dm[s] = 1.0
            xs = np.maximum(x, 1e-9)
            es[s] = xs ** prm["qs"]; ds[s] = prm["qs"] * xs ** (prm["qs"] - 1)
        l = g == 2
        if l.any():
            lat = prm["lattice"]
            rmin = lat["rho_min"]
            x = xb[l]
            rs = rmin + (1 - rmin) * x
            E, dEr = lat["E"](rs), lat["dE"](rs)
            S, dSr = lat["S"](rs), lat["dS"](rs)
            eE[l] = E; dE[l] = dEr * (1 - rmin)
            em[l] = rs; dm[l] = 1 - rmin
            es[l] = E / S; ds[l] = (dEr * S - E * dSr) / S ** 2 * (1 - rmin)
        return eE, em, es, dE, dm, ds

    # ------------------------------------------------------------------ расчёт
    def solve(self, xb, prm, cases, need_sens=True):
        """cases: [dict(rpm, thermal, allow (ne,) Па, mask (ne,) bool)]. Возврат: словарь с полями по случаям."""
        if not hasattr(self, "_free"):
            self.free_map()
        eE, em, es, dE, dm, ds = self.interp(xb, prm)
        Ee = eE * self.Es
        col = self.ei
        data = (Ee[:, None, None] * self.Ke0[col]).reshape(-1)[self._keep]
        nf = len(self._free)
        K = sp.csc_matrix((data, (self._ri, self._ci)), shape=(nf, nf))
        # MMD(Aᵀ+A) в симметричном режиме — в 2,4 раза быстрее COLAMD по умолчанию на этой сетке (замер 2026-11-14)
        lu = spla.splu(K, permc_spec="MMD_AT_PLUS_A", options=dict(SymmetricMode=True))
        F = np.zeros((self.ndof, len(cases)))
        for k, c in enumerate(cases):
            w2 = (2 * math.pi * c["rpm"] / 60.0) ** 2
            fe = (self.rho * w2 * em)[:, None] * self.fc0[col]
            if c.get("thermal", True):
                fe = fe + Ee[:, None] * self.fth0[col]
            np.add.at(F[:, k], self.edof, fe)
            F[:, k] += self.f_fix0 * (c["rpm"] / self.rpm0) ** 2 if self.rpm0 else 0.0
        U = np.zeros_like(F)
        U[self._free] = lu.solve(F[self._free])
        out = dict(U=U, cases=[])
        rhs_adj = np.zeros((self.ndof, len(cases)))
        per = []
        for k, c in enumerate(cases):
            ue = U[self.edof, k]
            eps = np.einsum("eij,ej->ei", self.Bc0[col], ue)
            if c.get("thermal", True):
                eps = eps - self.epsc0[col][:, None] * MVEC[None, :]
            sig = self.Es[:, None] * (eps @ self.D1.T)       # напряжение в сплошном материале
            vm = von_mises4(sig)
            sa = c["allow"]
            msk = c.get("mask", np.ones(self.ne, bool))
            u = es * vm / sa
            uu = np.where(msk, u, 0.0)
            P = prm["P"]
            umax = float(uu.max())
            PN = float(np.sum(np.maximum(uu, 1e-30) ** P) ** (1.0 / P))
            res = dict(sig=sig, vm=vm, u=u, umax=umax, PN=PN, w2=(2 * math.pi * c["rpm"] / 60.0) ** 2)
            if need_sens:
                dPN_du = np.where(msk, PN ** (1 - P) * np.maximum(uu, 1e-30) ** (P - 1), 0.0)
                res["dPN_expl"] = dPN_du * ds * vm / sa
                dPN_dsig = (dPN_du * es / sa)[:, None] * dvm_ds(sig, vm)
                # dσ/du_e = Es·D1·Bc
                ge = np.einsum("ei,ij,ejk->ek", dPN_dsig * self.Es[:, None], self.D1, self.Bc0[col])
                np.add.at(rhs_adj[:, k], self.edof, ge)
            per.append(res)
        if need_sens:
            L = np.zeros_like(rhs_adj)
            L[self._free] = lu.solve(rhs_adj[self._free])
            for k, c in enumerate(cases):
                ue = U[self.edof, k]; le = L[self.edof, k]
                dfe = (self.rho * per[k]["w2"] * dm)[:, None] * self.fc0[col]
                if c.get("thermal", True):
                    dfe = dfe + (dE * self.Es)[:, None] * self.fth0[col]
                dKu = (dE * self.Es)[:, None] * np.einsum("eij,ej->ei", self.Ke0[col], ue)
                per[k]["dPN"] = per[k]["dPN_expl"] + np.sum(le * (dfe - dKu), axis=1)
        out["cases"] = per
        out["em"] = em
        return out

    # ------------------------------------------------------------------ величины тела
    def mass(self, em):
        return float(np.sum(self.rho * self.Ve * em))

    def inertia(self, em):
        return float(np.sum(self.rho * self.Ve * em * self.rc ** 2))

    def hoop_mean(self, sig, em_mask):
        """средние окружные по меридиональному сечению сплошной части: ∫σθ dA / A"""
        A = self.dr * self.dz
        return float(np.sum(sig[em_mask, 2]) * A / (em_mask.sum() * A))

    def burst_margin(self, sig, solid):
        """запас по разрушающей частоте по средним окружным (как disk.py, k = 1)"""
        uts = self.mt["uts"](self.Te[solid]) * 1e6
        return float(math.sqrt(np.sum(uts) / max(np.sum(sig[solid, 2]), 1e-9)))


# =========================================================================== фильтр и проекция
def density_filter(model, R):
    """H (ne×ne), нормированный по строкам: x~ = H x. Вес — конус R − d (площадь элементов одинакова)."""
    pts = np.c_[model.rc, model.zc]
    tree = cKDTree(pts)
    D = tree.sparse_distance_matrix(tree, R, output_type="coo_matrix")
    w = np.maximum(R - D.data, 0.0)
    H = sp.coo_matrix((w, (D.row, D.col)), shape=(model.ne, model.ne)).tocsr()
    H = H + sp.identity(model.ne, format="csr") * R      # диагональ (d = 0 в sparse_distance_matrix пропускается)
    s = np.asarray(H.sum(axis=1)).ravel()
    return sp.diags(1.0 / s) @ H


def project(xt, beta, eta=0.5):
    if beta <= 0:
        return xt.copy(), np.ones_like(xt)
    a, b = math.tanh(beta * eta), math.tanh(beta * (1 - eta))
    t = np.tanh(beta * (xt - eta))
    return (a + t) / (a + b), beta * (1 - t ** 2) / (a + b)


# =========================================================================== MMA (двойственная подзадача)
class MMA:
    """MMA Svanberg (1987) с asymptote-правилами mmasub (2007); подзадача — через двойственную задачу
    (m множителей, L-BFGS-B). a0 = 1, a = 0, c = 1e3, d = 1: y — штраф невыполнимости."""

    def __init__(self, n, m, xmin=0.0, xmax=1.0, move=0.1, c=1e3, d=1.0):
        self.n, self.m = n, m
        self.xmin = np.full(n, xmin, float); self.xmax = np.full(n, xmax, float)
        self.move, self.c, self.d = move, np.full(m, c), np.full(m, d)
        self.xold1 = self.xold2 = None
        self.low = self.upp = None
        self.it = 0

    def step(self, x, f0, df0, g, dg):
        self.it += 1
        n, m = self.n, self.m
        xmin, xmax = self.xmin, self.xmax
        rng = xmax - xmin
        if self.it <= 2 or self.low is None:
            low = x - 0.5 * rng; upp = x + 0.5 * rng
        else:
            zzz = (x - self.xold1) * (self.xold1 - self.xold2)
            fac = np.ones(n); fac[zzz > 0] = 1.2; fac[zzz < 0] = 0.7
            low = x - fac * (self.xold1 - self.low); upp = x + fac * (self.upp - self.xold1)
            low = np.clip(low, x - 10 * rng, x - 0.01 * rng); upp = np.clip(upp, x + 0.01 * rng, x + 10 * rng)
        alfa = np.maximum.reduce([low + 0.1 * (x - low), x - self.move * rng, xmin])
        beta = np.minimum.reduce([upp - 0.1 * (upp - x), x + self.move * rng, xmax])
        ux1 = upp - x; xl1 = x - low; ux2 = ux1 ** 2; xl2 = xl1 ** 2
        xmami = np.maximum(rng, 1e-5)
        p0 = np.maximum(df0, 0); q0 = np.maximum(-df0, 0)
        pq0 = 0.001 * (p0 + q0) + 1e-5 / xmami
        p0 = (p0 + pq0) * ux2; q0 = (q0 + pq0) * xl2
        P = np.maximum(dg, 0); Q = np.maximum(-dg, 0)
        PQ = 0.001 * (P + Q) + 1e-5 / xmami
        P = (P + PQ) * ux2; Q = (Q + PQ) * xl2
        b = P @ (1 / ux1) + Q @ (1 / xl1) - g

        def xof(lam):
            p = p0 + lam @ P; q = q0 + lam @ Q
            sp_, sq = np.sqrt(p), np.sqrt(q)
            xx = (sp_ * low + sq * upp) / (sp_ + sq)
            return np.clip(xx, alfa, beta), p, q

        def negW(lam):
            xx, p, q = xof(lam)
            y = np.maximum(0.0, (lam - self.c) / self.d)
            W = np.sum(p / (upp - xx) + q / (xx - low)) - lam @ b + np.sum(self.c * y + 0.5 * self.d * y * y - lam * y)
            gr = P @ (1 / (upp - xx)) + Q @ (1 / (xx - low)) - b - y
            return -W, -gr

        r = minimize(negW, np.ones(m), jac=True, method="L-BFGS-B", bounds=[(0, None)] * m,
                     options=dict(maxiter=500, gtol=1e-10, ftol=1e-14))
        xnew, _, _ = xof(r.x)
        self.xold2 = self.xold1 if self.xold1 is not None else x.copy()
        self.xold1 = x.copy()
        self.low, self.upp = low, upp
        return xnew, r.x


# =========================================================================== оптимизация
def optimize(model, x0, design, cases, prm, objective="mass", w_J=0.0, R_filter=None, stages=None,
             log=None, u_target=None):
    """design: bool (ne,) — переменные проектирования (остальные элементы x̄ = 1).
    objective: 'mass' | 'inertia' | 'mix' (масса + w_J·момент инерции, оба нормированы на начальные).
    u_target: список предельных загрузок по случаям (по умолчанию 1)."""
    H = density_filter(model, R_filter)
    nd = int(design.sum())
    idx = np.nonzero(design)[0]
    m = len(cases)
    u_target = u_target or [1.0] * m
    mma = MMA(nd, m, move=prm.get("move", 0.1))
    x = x0[idx].copy()
    cnorm = np.ones(m)
    hist = []
    stages = stages or [dict(beta=1, it=40), dict(beta=2, it=30), dict(beta=4, it=30), dict(beta=8, it=30),
                        dict(beta=16, it=40)]
    full = np.ones(model.ne)
    # опорные величины для нормировки цели
    em_ref = np.ones(model.ne)
    M0, J0 = model.mass(em_ref), model.inertia(em_ref)
    it_total = 0
    for st in stages:
        beta = st["beta"]
        mma.move = st.get("move", prm.get("move", 0.1))
        for it in range(st["it"]):
            full[:] = 1.0
            full[idx] = x
            xt = H @ full
            xb, dxb = project(xt, beta)
            xb[~design] = 1.0; dxb[~design] = 0.0
            sol = model.solve(xb, prm, cases)
            eE, em, es, dE, dm, ds = model.interp(xb, prm)
            # цель
            Ms = model.rho * model.Ve
            if objective == "mass":
                f0 = np.sum(Ms * em) / M0; df_xb = Ms * dm / M0
            elif objective == "inertia":
                f0 = np.sum(Ms * em * model.rc ** 2) / J0; df_xb = Ms * dm * model.rc ** 2 / J0
            else:
                f0 = np.sum(Ms * em) / M0 + w_J * np.sum(Ms * em * model.rc ** 2) / J0
                df_xb = Ms * dm / M0 + w_J * Ms * dm * model.rc ** 2 / J0
            g = np.zeros(m); dg = np.zeros((m, nd))
            for k, c in enumerate(sol["cases"]):
                # адаптивная нормировка p-нормы (Le et al. 2010): c·PN ≈ max u
                ck = c["umax"] / max(c["PN"], 1e-12)
                cnorm[k] = ck if it_total == 0 else 0.5 * ck + 0.5 * cnorm[k]
                g[k] = cnorm[k] * c["PN"] / u_target[k] - 1.0
                dgx = cnorm[k] / u_target[k] * c["dPN"] * dxb
                dg[k] = (H.T @ dgx)[idx]
            df0 = (H.T @ (df_xb * dxb))[idx]
            xnew, lam = mma.step(x, f0, df0, g, dg)
            ch = float(np.max(np.abs(xnew - x)))
            x = xnew
            rec = dict(it=it_total, beta=beta, f0=float(f0), mass=model.mass(em), J=model.inertia(em),
                       umax=[c["umax"] for c in sol["cases"]], PN=[c["PN"] for c in sol["cases"]],
                       g=g.tolist(), change=ch, grey=float(np.mean(4 * xb[design] * (1 - xb[design]))))
            hist.append(rec)
            if log:
                log(rec)
            it_total += 1
            if it > 10 and ch < st.get("tol", 0.005) and max(g) < 0.01:
                break
    full[:] = 1.0; full[idx] = x
    xt = H @ full
    xb, _ = project(xt, stages[-1]["beta"])
    xb[~design] = 1.0
    return dict(x=full, xt=xt, xb=xb, hist=hist, H=H)


def evaluate(model, xb, prm, cases, threshold=None):
    """Проверочный расчёт: при threshold — двоичный проект (x̄ ≥ порога -> 1, иначе пустота Emin)."""
    xe = xb.copy()
    if threshold is not None:
        xe = np.where(xb >= threshold, 1.0, 0.0)
        xe[model.group == 0] = 1.0
    p2 = dict(prm)
    sol = model.solve(xe, p2, cases, need_sens=False)
    eE, em, es, *_ = model.interp(xe, p2)
    return xe, sol, em
