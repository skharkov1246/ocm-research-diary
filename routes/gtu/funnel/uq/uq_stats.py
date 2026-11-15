"""Устойчивый отбор: статистика выборки - доля, доверительный интервал Уилсона, квантили, индексы Соболя первого порядка по
готовой выборке (оценка по разбиению, Plischke - Borgonovo - Smith, 2013, по памяти). Всё - оценка."""
import math

import numpy as np


def wilson(k, n, z=1.96):
    """95 % интервал Уилсона для доли k/n"""
    if n <= 0:
        return None
    p = k / n
    den = 1 + z * z / n
    c = (p + z * z / (2 * n)) / den
    h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / den
    return [round(max(0.0, c - h), 4), round(min(1.0, c + h), 4)]


def q(a, ps=(10, 50, 90)):
    """квантили массива без NaN; None, если пусто"""
    a = np.asarray([x for x in a if x is not None and np.isfinite(x)], dtype=float)
    if a.size == 0:
        return [None] * len(ps)
    return [round(float(v), 5) for v in np.percentile(a, ps)]


def bins_default(n):
    return int(max(4, min(32, round(math.sqrt(n) / 1.5))))


def sobol1(X, y, bins=None):
    """индексы первого порядка S_i = Var(E[y|x_i]) / Var(y) по одной выборке: x_i делится на M равных по числу точек слоёв (или по
    значениям, если их не больше M), корреляционное отношение η² исправляется на смещение (G - 1)/(N - 1) - ожидаемое η² при
    независимости. Возврат - список S_i (>= 0) или None, если y постоянен"""
    X = np.asarray(X, dtype=float)
    y = np.asarray(y, dtype=float)
    N = y.size
    if N < 8:
        return None
    v = float(np.var(y))
    if v <= 1e-14 * max(1.0, float(np.mean(y * y))):
        return None
    M = bins or bins_default(N)
    mu = float(np.mean(y))
    out = []
    for i in range(X.shape[1]):
        x = X[:, i]
        uq = np.unique(x)
        if uq.size <= M:
            groups = [y[x == u] for u in uq]
        else:
            order = np.argsort(x, kind="stable")
            groups = [y[g] for g in np.array_split(order, M)]
        G = len(groups)
        eta2 = sum(g.size * (float(np.mean(g)) - mu) ** 2 for g in groups if g.size) / (N * v)
        bias = (G - 1) / (N - 1)
        out.append(max(0.0, (eta2 - bias) / (1 - bias)))
    return out


def top(names, S, k=5):
    if S is None:
        return None
    order = sorted(range(len(S)), key=lambda i: (-S[i], i))[:k]
    return [[names[i], round(float(S[i]), 4)] for i in order]
