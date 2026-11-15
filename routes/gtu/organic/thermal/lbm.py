"""Решёточный метод Больцмана в периодической ячейке: течение D3Q19 (TRT) и пассивный скаляр D3Q7 (TRT).

Второй, независимый от voxflow.py метод для коэффициентов переноса решёток TPMS, и единственный здесь —
с инерцией (конечные Re). Течение задаётся объёмной силой по оси 0 (Guo et al. 2002, антисимметричная часть),
прилипание — отскок на середине связи (half-way bounce-back), Λ = 3/16 (положение стенки точно посередине
для течения Пуазейля — Ginzburg, d'Humières 2003). Скаляр: T = γ·x + φ, φ = 0 на стенке — антиотскок
(anti-bounce-back), источник −γ·u_x; это та же постановка H1, что в voxflow.scalar_h1.

Единицы решётки: h = 1, Δt = 1, ρ ≈ 1. Плотность жидкости меняется на ~Ma², поэтому пределы: u < 0,1.
Сходимость — по относительному изменению средней скорости и φ_b за окно шагов.
"""
import numpy as np

C19 = np.array([[0, 0, 0],
                [1, 0, 0], [-1, 0, 0], [0, 1, 0], [0, -1, 0], [0, 0, 1], [0, 0, -1],
                [1, 1, 0], [-1, -1, 0], [1, -1, 0], [-1, 1, 0],
                [1, 0, 1], [-1, 0, -1], [1, 0, -1], [-1, 0, 1],
                [0, 1, 1], [0, -1, -1], [0, 1, -1], [0, -1, 1]])
W19 = np.array([1 / 3] + [1 / 18] * 6 + [1 / 36] * 12)
C7 = C19[:7]
W7 = np.array([1 / 4] + [1 / 8] * 6)


def _opp(C):
    return np.array([int(np.where((C == -c).all(1))[0][0]) for c in C])


OPP19, OPP7 = _opp(C19), _opp(C7)


class LBMCell:
    def __init__(self, fluid, nu, force, kappa=None, Lambda=3 / 16, Lambda_g=1 / 4):
        self.fluid = np.asarray(fluid, bool)
        self.solid = ~self.fluid
        self.N = self.fluid.shape
        self.nu, self.F = float(nu), float(force)
        tp = 3 * nu + 0.5
        tm = 0.5 + Lambda / (tp - 0.5)
        self.wp, self.wm = 1 / tp, 1 / tm
        self.f = np.empty((19,) + self.N)
        for i in range(19):
            self.f[i] = W19[i]
        # связи, ведущие из твёрдого узла в жидкий: после переноса заменяются отскоком
        self.bb = [np.roll(self.solid, tuple(C19[i]), axis=(0, 1, 2)) & self.fluid for i in range(19)]
        self.kappa = kappa
        if kappa is not None:
            tgm = 4 * kappa + 0.5                      # диффузия — нечётные моменты; cs² = 1/4 у D3Q7
            tgp = 0.5 + Lambda_g / (tgm - 0.5)
            self.gwp, self.gwm = 1 / tgp, 1 / tgm
            self.g = np.zeros((7,) + self.N)
            self.bb7 = [np.roll(self.solid, tuple(C7[i]), axis=(0, 1, 2)) & self.fluid for i in range(7)]
        self.it = 0

    def macros(self):
        f = self.f
        rho = f.sum(0)
        jx = f[1] - f[2] + f[7] - f[8] + f[9] - f[10] + f[11] - f[12] + f[13] - f[14]
        jy = f[3] - f[4] + f[7] - f[8] - f[9] + f[10] + f[15] - f[16] + f[17] - f[18]
        jz = f[5] - f[6] + f[11] - f[12] - f[13] + f[14] + f[15] - f[16] - f[17] + f[18]
        ux = (jx + 0.5 * self.F) / rho
        uy, uz = jy / rho, jz / rho
        for a in (ux, uy, uz):
            a[self.solid] = 0.0
        return rho, ux, uy, uz

    def step(self):
        f = self.f
        rho, ux, uy, uz = self.macros()
        usq = ux * ux + uy * uy + uz * uz
        fpost = np.empty_like(f)
        for i in range(19):
            j = OPP19[i]
            if j < i:
                continue
            cu = C19[i, 0] * ux + C19[i, 1] * uy + C19[i, 2] * uz
            feq_s = W19[i] * rho * (1 + 4.5 * cu * cu - 1.5 * usq)   # симметричная часть равновесия
            feq_a = W19[i] * rho * 3 * cu                              # антисимметричная
            fs = 0.5 * (f[i] + f[j])
            fa = 0.5 * (f[i] - f[j])
            Si = W19[i] * 3 * C19[i, 0] * self.F * (1 - 0.5 * self.wm)  # сила Guo, антисимметричная часть
            if i == j:
                fpost[i] = f[i] - self.wp * (fs - feq_s)
            else:
                fpost[i] = f[i] - self.wp * (fs - feq_s) - self.wm * (fa - feq_a) + Si
                fpost[j] = f[j] - self.wp * (fs - feq_s) + self.wm * (fa - feq_a) - Si
        for i in range(19):
            fi = np.roll(fpost[i], tuple(C19[i]), axis=(0, 1, 2))
            m = self.bb[i]
            fi[m] = fpost[OPP19[i]][m]
            f[i] = fi
        if self.kappa is not None:
            self._step_g(ux, uy, uz)
        self.it += 1

    def _step_g(self, ux, uy, uz):
        g = self.g
        phi = g.sum(0)
        src = -ux                                       # γ = 1
        phi = phi + 0.5 * src
        phi[self.solid] = 0.0
        gpost = np.empty_like(g)
        for i in range(7):
            j = OPP7[i]
            if j < i:
                continue
            cu = C7[i, 0] * ux + C7[i, 1] * uy + C7[i, 2] * uz
            geq_s = W7[i] * phi
            geq_a = W7[i] * phi * 4 * cu
            gs = 0.5 * (g[i] + g[j])
            ga = 0.5 * (g[i] - g[j])
            Ss = W7[i] * src * (1 - 0.5 * self.gwp)
            if i == j:
                gpost[i] = g[i] - self.gwp * (gs - geq_s) + Ss
            else:
                gpost[i] = g[i] - self.gwp * (gs - geq_s) - self.gwm * (ga - geq_a) + Ss
                gpost[j] = g[j] - self.gwp * (gs - geq_s) + self.gwm * (ga - geq_a) + Ss
        for i in range(7):
            gi = np.roll(gpost[i], tuple(C7[i]), axis=(0, 1, 2))
            m = self.bb7[i]
            gi[m] = -gpost[OPP7[i]][m]                  # антиотскок: φ = 0 на середине связи
            g[i] = gi

    def observables(self):
        rho, ux, uy, uz = self.macros()
        out = dict(U=float(ux[self.fluid].mean()), umax=float(np.sqrt(ux ** 2 + uy ** 2 + uz ** 2).max()),
                   rho_min=float(rho[self.fluid].min()), rho_max=float(rho[self.fluid].max()))
        if self.kappa is not None:
            phi = self.g.sum(0) - 0.5 * ux
            w = ux[self.fluid]
            out["phib"] = float((w * phi[self.fluid]).sum() / w.sum())
            out["Q"] = float(w.sum())                   # γ·Σu_x — тепло в жидкость на ячейку (κ-единицы)
        return out

    def run(self, max_steps=200000, check=500, tol=1e-7, min_steps=2000, verbose=False):
        last = None
        hist = []
        for s in range(max_steps):
            self.step()
            if (s + 1) % check == 0:
                ob = self.observables()
                key = np.array([ob["U"]] + ([ob["phib"]] if "phib" in ob else []))
                hist.append((self.it, ob))
                if verbose:
                    print(self.it, ob, flush=True)
                if last is not None and s + 1 >= min_steps:
                    ch = np.max(np.abs(key - last) / np.maximum(np.abs(key), 1e-30))
                    if ch < tol:
                        return dict(converged=True, steps=self.it, obs=ob)
                last = key
        return dict(converged=False, steps=self.it, obs=self.observables())
