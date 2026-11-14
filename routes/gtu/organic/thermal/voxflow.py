"""Воксельный решатель течения Стокса и теплообмена в периодической ячейке (2D/3D), схема MAC.

Назначение: коэффициенты переноса решёток TPMS (гироид, Schwarz-D) в пределе малых Re —
f·Re (Фаннинг) и Nu_H1 (постоянный осевой тепловой поток, температура стенки однородна в сечении),
аналоги констант плоского канала 24 и 8,235, которыми пользуется organic/cycles.

Метод: разнесённая сетка (давление в центрах ячеек, компоненты скорости на гранях), прилипание на
гранях между жидкой и твёрдой ячейкой (нормальная компонента = 0 на грани, касательная — зеркальное
отражение: стенка посередине между узлами), периодичность по всем осям, течение задано объёмной силой
по оси 0. Седловая система [[A, G], [G^T, 0]] решается прямым методом (SuperLU, scipy).
Температура: T = γ·x + φ, φ периодична, φ = 0 на стенке (H1), уравнение u·∇φ − κ∇²φ = −γ·u_x,
конвекция — центральная или противопоточная. Единицы решётки: h = 1, μ = 1, κ = 1, γ = 1.

Проверка (test_thermal.py): плоский канал — f·Re = 24, Nu_H1 = 8,235 (точно на сетке, вдоль стенок);
квадратный канал — f·Re = 14,227, Nu_H1 = 3,608 (Shah, London 1978 — справочные значения).
"""
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
from scipy import ndimage


def _roll(a, shift, axis):
    return np.roll(a, shift, axis=axis)


def percolating(mask):
    """оставить наибольшую связную (с учётом периодичности) компоненту жидкости"""
    lab, n = ndimage.label(mask)
    if n <= 1:
        return mask.copy()
    parent = list(range(n + 1))

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    nd = mask.ndim
    for ax in range(nd):
        a = np.take(lab, 0, axis=ax)
        b = np.take(lab, -1, axis=ax)
        m = (a > 0) & (b > 0)
        for x, y in zip(a[m].ravel(), b[m].ravel()):
            rx, ry = find(int(x)), find(int(y))
            if rx != ry:
                parent[rx] = ry
    roots = np.array([find(i) for i in range(n + 1)])
    rl = roots[lab]
    rl[~mask] = 0
    sizes = np.bincount(rl.ravel())
    sizes[0] = 0
    return rl == np.argmax(sizes)


class PeriodicVoxelFlow:
    def __init__(self, fluid):
        self.fluid = np.asarray(fluid, bool)
        self.N = self.fluid.shape
        self.nd = self.fluid.ndim
        self.ncell = self.fluid.size
        self._build()

    def _build(self):
        f, nd, N = self.fluid, self.nd, self.N
        # активные грани компоненты d: грань с индексом i лежит между ячейками i-e_d и i
        self.act = [f & _roll(f, 1, d) for d in range(nd)]
        self.fid = []
        off = 0
        for d in range(nd):
            idx = -np.ones(N, np.int64)
            k = int(self.act[d].sum())
            idx[self.act[d]] = np.arange(off, off + k)
            self.fid.append(idx)
            off += k
        self.nu = off
        self.pid = -np.ones(N, np.int64)
        kc = int(f.sum())
        self.pid[f] = np.arange(off, off + kc)
        self.np_ = kc
        rows, cols, vals = [], [], []
        diag = np.zeros(self.nu)
        for d in range(nd):
            act = self.act[d]
            me = self.fid[d][act]
            for e in range(nd):
                for s in (+1, -1):
                    nb_act = _roll(self.act[d], -s, e)[act]   # грань f + s·e_e
                    nb_id = _roll(self.fid[d], -s, e)[act]
                    # связь с активной соседней гранью
                    rows.append(me[nb_act]); cols.append(nb_id[nb_act]); vals.append(-np.ones(nb_act.sum()))
                    if e == d:
                        diag[me] += 1.0                       # неактивная — стенка, u = 0 на грани
                    else:
                        diag[me] += np.where(nb_act, 1.0, 2.0)  # неактивная — зеркало, стенка посередине
            # градиент давления: + p[i] - p[i-e_d]
            pr = self.pid[act]
            pl = _roll(self.pid, 1, d)[act]
            rows += [me, me]; cols += [pr, pl]; vals += [np.ones(me.size), -np.ones(me.size)]
        rows.append(np.arange(self.nu)); cols.append(np.arange(self.nu)); vals.append(diag)
        r = np.concatenate(rows); c = np.concatenate(cols); v = np.concatenate(vals)
        n = self.nu + self.np_
        M = sp.coo_matrix((v, (r, c)), shape=(n, n)).tocsr()
        # симметричное замыкание: строки неразрывности = транспонированный градиент
        G = M[:self.nu, self.nu:]
        A = M[:self.nu, :self.nu]
        # закрепить давление в одной ячейке (один связный кластер)
        pin = 0
        GT = G.T.tolil()
        GT[pin, :] = 0
        GT = GT.tocsr()
        Zp = sp.csr_matrix((self.np_, self.np_))
        Zp = Zp.tolil(); Zp[pin, pin] = 1.0; Zp = Zp.tocsr()
        Gc = G.tolil(); Gc[:, pin] = 0; Gc = Gc.tocsr()
        self.K = sp.bmat([[A, Gc], [GT, Zp]], format="csc")
        self.A = A

    def stokes(self, axis=0, g=1.0):
        b = np.zeros(self.nu + self.np_)
        sl = slice(sum(int(self.act[d].sum()) for d in range(axis)),
                   sum(int(self.act[d].sum()) for d in range(axis + 1)))
        b[sl] = g
        lu = spla.splu(self.K, permc_spec="COLAMD")
        x = lu.solve(b)
        res = np.linalg.norm(self.K @ x - b) / np.linalg.norm(b)
        self.x = x
        U = []
        for d in range(self.nd):
            u = np.zeros(self.N)
            u[self.act[d]] = x[self.fid[d][self.act[d]]]
            U.append(u)
        p = np.zeros(self.N)
        p[self.fluid] = x[self.pid[self.fluid]]
        self.U, self.p, self.axis, self.g = U, p, axis, g
        return dict(residual=float(res))

    def mean_velocity(self):
        """средняя скорость в жидкости вдоль оси течения (межпоровая): ∫u dV / V_f"""
        u = self.U[self.axis]
        return float(u.sum() / self.fluid.sum())

    def divergence_max(self):
        div = np.zeros(self.N)
        for d in range(self.nd):
            div += _roll(self.U[d], -1, d) - self.U[d]
        return float(np.abs(div[self.fluid]).max())

    def scalar_h1(self, Pe_scale, scheme="central"):
        """φ: u·∇φ − ∇²φ = −u_x (γ = 1, κ = 1); скорость = Pe_scale · u_stokes. Возвращает φ и величины."""
        f, nd, N = self.fluid, self.nd, self.N
        cid = -np.ones(N, np.int64)
        nc = int(f.sum())
        cid[f] = np.arange(nc)
        rows, cols, vals = [], [], []
        diag = np.zeros(nc)
        U = [Pe_scale * u for u in self.U]
        me_all = cid[f]
        for d in range(nd):
            for s in (+1, -1):
                nb_f = _roll(f, -s, d)[f]
                nb = _roll(cid, -s, d)[f]
                # диффузия
                rows.append(me_all[nb_f]); cols.append(nb[nb_f]); vals.append(-np.ones(nb_f.sum()))
                diag += np.where(nb_f, 1.0, 2.0)           # стенка на грани: φ = 0 на расстоянии h/2
                # конвекция через грань (s=+1: грань i+1 — «выход» вверх по оси; s=-1: грань i — «вход»)
                uf = (_roll(U[d], -1, d) if s == +1 else U[d])[f]
                F = s * uf                                   # поток наружу из ячейки
                if scheme == "central":
                    w_me = 0.5 * F; w_nb = 0.5 * F
                else:
                    w_me = np.maximum(F, 0.0); w_nb = np.minimum(F, 0.0)
                diag += np.where(nb_f, w_me, 0.0)
                rows.append(me_all[nb_f]); cols.append(nb[nb_f]); vals.append(w_nb[nb_f])
        rows.append(np.arange(nc)); cols.append(np.arange(nc)); vals.append(diag)
        M = sp.coo_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))), shape=(nc, nc)).tocsc()
        ax = self.axis
        ux_c = 0.5 * (U[ax] + _roll(U[ax], -1, ax))[f]
        rhs = -ux_c
        phi_v = spla.spsolve(M, rhs)
        phi = np.zeros(N)
        phi[f] = phi_v
        # тепловой поток в стенку (жидкость получает тепло): сумма связей со стенкой 2·(0 − φ)
        qwall = 0.0
        for d in range(nd):
            for s in (+1, -1):
                nb_f = _roll(f, -s, d)[f]
                qwall += float((2.0 * (0.0 - phi_v))[~nb_f].sum())
        phib = float((ux_c * phi_v).sum() / ux_c.sum())
        return dict(phi=phi, phib=phib, Q_balance=float(ux_c.sum()), Q_wall=qwall,
                    Umean=float(ux_c.sum() / nc))
