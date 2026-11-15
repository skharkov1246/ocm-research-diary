"""CMA-ES (μ/μ_w, λ) по Hansen 2016, «The CMA Evolution Strategy: A Tutorial» (arXiv:1604.00772): весовая
рекомбинация, пути эволюции p_σ и p_c, ранг-1 и ранг-μ обновление ковариации, кумулятивное управление шагом.

Границы [0, 1]: точка выборки x оценивается в проекции clip(x) со штрафом w·‖x − clip(x)‖² (простой способ из того же
руководства, раздел «Boundaries and constraints»); w подстраивается к разбросу значений f. Минимизация.
"""
import math

import numpy as np


class CMAES:
    def __init__(self, x0, sigma0, rng, popsize=None, lb=0.0, ub=1.0):
        self.rng = rng
        self.n = n = len(x0)
        self.mean = np.asarray(x0, float).copy()
        self.sigma = float(sigma0)
        self.lb, self.ub = lb, ub
        self.lam = popsize or 4 + int(3 * math.log(n))
        self.mu = self.lam // 2
        w = math.log(self.mu + 0.5) - np.log(np.arange(1, self.mu + 1))
        self.w = w / w.sum()
        self.mueff = 1.0 / np.sum(self.w ** 2)
        self.cc = (4 + self.mueff / n) / (n + 4 + 2 * self.mueff / n)
        self.cs = (self.mueff + 2) / (n + self.mueff + 5)
        self.c1 = 2 / ((n + 1.3) ** 2 + self.mueff)
        self.cmu = min(1 - self.c1, 2 * (self.mueff - 2 + 1 / self.mueff) / ((n + 2) ** 2 + self.mueff))
        self.damps = 1 + 2 * max(0.0, math.sqrt((self.mueff - 1) / (n + 1)) - 1) + self.cs
        self.chin = math.sqrt(n) * (1 - 1 / (4 * n) + 1 / (21 * n * n))
        self.pc = np.zeros(n); self.ps = np.zeros(n)
        self.C = np.eye(n); self.B = np.eye(n); self.D = np.ones(n)
        self.eigeneval = 0; self.counteval = 0; self.gen = 0
        self.pen_w = 1e2

    def ask(self):
        Z = self.rng.standard_normal((self.lam, self.n))
        Y = Z * self.D @ self.B.T
        X = self.mean + self.sigma * Y
        return X

    def repair(self, X):
        Xc = np.clip(X, self.lb, self.ub)
        return Xc, np.sum((X - Xc) ** 2, axis=1)

    def tell(self, X, f_raw, pen):
        f_raw = np.asarray(f_raw, float)
        finite = np.isfinite(f_raw)
        worst = float(f_raw[finite].max()) if finite.any() else 0.0
        span = float(np.ptp(f_raw[finite])) if finite.sum() > 1 else 1.0
        self.pen_w = 1e2 * (1.0 + span)
        f = np.where(finite, f_raw, worst + 1e9) + self.pen_w * pen
        o = np.argsort(f, kind="stable")
        n = self.n
        self.counteval += len(X)
        old = self.mean.copy()
        Xs = X[o[: self.mu]]
        self.mean = self.w @ Xs
        y = (self.mean - old) / self.sigma
        invsqrtC = self.B @ np.diag(1 / self.D) @ self.B.T
        self.ps = (1 - self.cs) * self.ps + math.sqrt(self.cs * (2 - self.cs) * self.mueff) * invsqrtC @ y
        hsig = np.linalg.norm(self.ps) / math.sqrt(1 - (1 - self.cs) ** (2 * self.counteval / self.lam)) / self.chin < 1.4 + 2 / (n + 1)
        self.pc = (1 - self.cc) * self.pc + hsig * math.sqrt(self.cc * (2 - self.cc) * self.mueff) * y
        art = (Xs - old) / self.sigma
        self.C = ((1 - self.c1 - self.cmu) * self.C + self.c1 * (np.outer(self.pc, self.pc) + (1 - hsig) * self.cc * (2 - self.cc) * self.C)
                  + self.cmu * (art.T * self.w) @ art)
        self.sigma *= math.exp((self.cs / self.damps) * (np.linalg.norm(self.ps) / self.chin - 1))
        self.sigma = min(self.sigma, 1.0)
        if self.counteval - self.eigeneval > self.lam / (self.c1 + self.cmu) / n / 10:
            self.eigeneval = self.counteval
            self.C = np.triu(self.C) + np.triu(self.C, 1).T
            d2, self.B = np.linalg.eigh(self.C)
            self.D = np.sqrt(np.maximum(d2, 1e-20))
        self.gen += 1
        return f[o[0]], X[o[0]]


def fmin(f, x0, sigma0, rng, max_evals, ftarget=-np.inf, popsize=None, lb=0.0, ub=1.0, batch_f=None, tol_sigma=1e-12):
    """минимизация f (или batch_f(X) → значения) в [lb, ub]^n. Возврат dict(x, f, evals, history)."""
    es = CMAES(x0, sigma0, rng, popsize, lb, ub)
    best_x, best_f = np.clip(np.asarray(x0, float), lb, ub), np.inf
    hist = []
    while es.counteval < max_evals:
        X = es.ask()
        Xc, pen = es.repair(X)
        fv = np.asarray(batch_f(Xc) if batch_f else [f(x) for x in Xc], float)
        es.tell(X, fv, pen)
        i = int(np.argmin(np.where(np.isfinite(fv), fv, np.inf) + 1e30 * (pen > 0)))
        if pen[i] == 0 and fv[i] < best_f:
            best_f, best_x = float(fv[i]), Xc[i].copy()
        hist.append(dict(evals=es.counteval, best=best_f, sigma=es.sigma))
        if best_f <= ftarget or es.sigma * np.max(es.D) < tol_sigma:
            break
    return dict(x=best_x, f=best_f, evals=es.counteval, history=hist)
