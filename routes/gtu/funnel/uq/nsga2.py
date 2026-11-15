"""NSGA-II (Deb, Pratap, Agarwal, Meyarivan, 2002) - своя реализация без зависимостей, детерминированная (numpy PCG64).

Переменные - в [0, 1]; цели минимизируются; ограничения - одно число «превышение» V >= 0 (0 - допустим) с правилом Деба:
допустимый лучше недопустимого, из недопустимых - меньшее превышение, из допустимых - доминирование по Парето. Отбор родителей -
двоичный турнир (ранг, затем расстояние скученности), скрещивание SBX (η_c = 15, вероятность 0,9, каждая переменная с вероятностью
0,5), мутация полиномиальная (η_m = 20, вероятность 1/n). Дискретные оси задачи кодируются числом в [0, 1] (слой значения).
Проверка - test_uq.py: фронт ZDT1 (известный: f2 = 1 − √f1) и ZDT1 с ограничением.
"""
import numpy as np


def lhs(n, d, rng):
    u = np.empty((n, d))
    for k in range(d):
        u[:, k] = (rng.permutation(n) + rng.random(n)) / n
    return u


def dominance(F):
    """D[i, j] = i доминирует j (все цели не хуже и хотя бы одна лучше)"""
    le = np.all(F[:, None, :] <= F[None, :, :], axis=2)
    lt = np.any(F[:, None, :] < F[None, :, :], axis=2)
    return le & lt


def nd_sort(F, V):
    """фронты с ограничениями (Деб): сначала фронты допустимых по Парето, затем недопустимые по возрастанию превышения"""
    n = len(F)
    feas = np.where(V <= 0)[0]
    fronts = []
    if feas.size:
        D = dominance(F[feas])
        cnt = D.sum(axis=0)
        cur = list(np.where(cnt == 0)[0])
        while cur:
            fronts.append([int(feas[i]) for i in cur])
            nxt = []
            for i in cur:
                for j in np.where(D[i])[0]:
                    cnt[j] -= 1
                    if cnt[j] == 0:
                        nxt.append(j)
            cur = nxt
    inf = np.where(V > 0)[0]
    if inf.size:
        order = inf[np.argsort(V[inf], kind="stable")]
        vals = V[order]
        k = 0
        while k < len(order):
            m = k
            while m + 1 < len(order) and vals[m + 1] == vals[k]:
                m += 1
            fronts.append([int(x) for x in order[k:m + 1]])
            k = m + 1
    rank = np.empty(n, dtype=int)
    for r, f in enumerate(fronts):
        rank[f] = r
    return fronts, rank


def crowding(F, front):
    f = np.asarray(front)
    d = np.zeros(len(f))
    if len(f) <= 2:
        d[:] = np.inf
        return d
    for m in range(F.shape[1]):
        v = F[f, m]
        o = np.argsort(v, kind="stable")
        span = v[o[-1]] - v[o[0]]
        d[o[0]] = d[o[-1]] = np.inf
        if span > 0:
            d[o[1:-1]] += (v[o[2:]] - v[o[:-2]]) / span
    return d


def _tournament(rng, rank, cd, V, n):
    a, b = rng.integers(0, len(rank), size=(2, n))
    pick = np.where(
        (V[a] <= 0) & (V[b] > 0), a, np.where(
            (V[b] <= 0) & (V[a] > 0), b, np.where(
                (V[a] > 0) & (V[b] > 0), np.where(V[a] <= V[b], a, b), np.where(
                    rank[a] < rank[b], a, np.where(rank[b] < rank[a], b, np.where(cd[a] >= cd[b], a, b))))))
    return pick


def sbx(rng, p1, p2, eta=15.0, pc=0.9):
    c1, c2 = p1.copy(), p2.copy()
    n, d = p1.shape
    do = rng.random(n) < pc
    for i in np.where(do)[0]:
        for k in range(d):
            if rng.random() > 0.5:
                continue
            x1, x2 = p1[i, k], p2[i, k]
            if abs(x1 - x2) < 1e-14:
                continue
            y1, y2 = min(x1, x2), max(x1, x2)
            r = rng.random()
            beta = 1.0 + 2.0 * y1 / (y2 - y1)
            alpha = 2.0 - beta ** -(eta + 1)
            bq = (r * alpha) ** (1 / (eta + 1)) if r <= 1 / alpha else (1 / (2 - r * alpha)) ** (1 / (eta + 1))
            a1 = 0.5 * ((y1 + y2) - bq * (y2 - y1))
            beta = 1.0 + 2.0 * (1.0 - y2) / (y2 - y1)
            alpha = 2.0 - beta ** -(eta + 1)
            bq = (r * alpha) ** (1 / (eta + 1)) if r <= 1 / alpha else (1 / (2 - r * alpha)) ** (1 / (eta + 1))
            a2 = 0.5 * ((y1 + y2) + bq * (y2 - y1))
            a1, a2 = min(max(a1, 0.0), 1.0), min(max(a2, 0.0), 1.0)
            if rng.random() <= 0.5:
                a1, a2 = a2, a1
            c1[i, k], c2[i, k] = a1, a2
    return c1, c2


def pm(rng, X, eta=20.0, p=None):
    Y = X.copy()
    n, d = X.shape
    p = 1.0 / d if p is None else p
    mask = rng.random((n, d)) < p
    r = rng.random((n, d))
    x = Y[mask]
    rr = r[mask]
    d1, d2 = x, 1.0 - x
    mp = 1.0 / (eta + 1)
    lo = rr < 0.5
    xy = np.where(lo, 1.0 - d1, 1.0 - d2)
    val = np.where(lo, 2 * rr + (1 - 2 * rr) * xy ** (eta + 1), 2 * (1 - rr) + 2 * (rr - 0.5) * xy ** (eta + 1))
    dq = np.where(lo, val ** mp - 1.0, 1.0 - val ** mp)
    Y[mask] = np.clip(x + dq, 0.0, 1.0)
    return Y


def survive(F, V, n):
    fronts, rank = nd_sort(F, V)
    keep, cd_all = [], np.zeros(len(F))
    for f in fronts:
        cd = crowding(F, f)
        cd_all[f] = cd
        if len(keep) + len(f) <= n:
            keep += f
        else:
            o = np.argsort(-cd, kind="stable")[: n - len(keep)]
            keep += [f[i] for i in o]
            break
    keep = np.array(keep)
    return keep, rank[keep], cd_all[keep]


def run(evaluate, n_var, n_obj, pop=100, gens=50, seed=0, init=None, eta_c=15.0, eta_m=20.0, pc=0.9, pmut=None, callback=None):
    """evaluate(X) -> (F n×n_obj, V n) - цели (минимум) и превышение. init - начальные точки (k×n_var, k <= pop).
    Возврат: X, F, V последней популяции, rank, история (число допустимых, размер фронта по поколениям)"""
    rng = np.random.Generator(np.random.PCG64(int(seed)))
    X = lhs(pop, n_var, rng)
    if init is not None and len(init):
        k = min(len(init), pop)
        X[:k] = np.clip(np.asarray(init, dtype=float)[:k], 0.0, 1.0)
    F, V = evaluate(X)
    keep, rank, cd = survive(F, V, pop)
    X, F, V = X[keep], F[keep], V[keep]
    hist = []
    for g in range(gens):
        a = _tournament(rng, rank, cd, V, pop)
        b = _tournament(rng, rank, cd, V, pop)
        c1, c2 = sbx(rng, X[a], X[b], eta_c, pc)
        Xo = pm(rng, np.vstack([c1, c2])[:pop], eta_m, pmut)
        Fo, Vo = evaluate(Xo)
        XX, FF, VV = np.vstack([X, Xo]), np.vstack([F, Fo]), np.concatenate([V, Vo])
        keep, rank, cd = survive(FF, VV, pop)
        X, F, V = XX[keep], FF[keep], VV[keep]
        nf = int((V <= 0).sum())
        hist.append(dict(gen=g + 1, feasible=nf, front=int(((rank == 0) & (V <= 0)).sum())))
        if callback:
            callback(g + 1, X, F, V, rank)
    return dict(X=X, F=F, V=V, rank=rank, history=hist)


def front_mask(F, V):
    """недоминируемые допустимые"""
    feas = np.where(V <= 0)[0]
    if not feas.size:
        return np.zeros(len(F), dtype=bool)
    D = dominance(F[feas])
    nd = feas[D.sum(axis=0) == 0]
    m = np.zeros(len(F), dtype=bool)
    m[nd] = True
    return m


def hv_mc(G, n=20000, seed=7):
    """гиперобъём фронта G (цели нормированы: 0 - лучшее, 1 - точка отсчёта) методом Монте-Карло по фиксированным точкам"""
    G = np.asarray(G, dtype=float)
    if G.size == 0:
        return 0.0
    G = G[np.all(G < 1.0, axis=1)]
    if not len(G):
        return 0.0
    rng = np.random.Generator(np.random.PCG64(seed))
    P = rng.random((n, G.shape[1]))
    hit = np.zeros(n, dtype=bool)
    for k in range(0, len(G), 64):
        g = G[k:k + 64]
        hit |= np.any(np.all(g[None, :, :] <= P[:, None, :], axis=2), axis=1)
    return float(hit.mean())


def hv2(F, ref):
    """точный гиперобъём для двух целей"""
    P = np.asarray(sorted([tuple(f) for f in F if f[0] < ref[0] and f[1] < ref[1]]))
    if not len(P):
        return 0.0
    hv, y = 0.0, ref[1]
    for x, f2 in P:
        if f2 < y:
            hv += (ref[0] - x) * (y - f2)
            y = f2
    return hv


# =========================================================================== тестовые задачи
def zdt1(X):
    X = np.asarray(X)
    f1 = X[:, 0]
    g = 1 + 9 * X[:, 1:].sum(axis=1) / (X.shape[1] - 1)
    f2 = g * (1 - np.sqrt(f1 / g))
    return np.stack([f1, f2], axis=1), np.zeros(len(X))


def zdt1_constrained(X, x1_min=0.3):
    F, _ = zdt1(X)
    return F, np.maximum(0.0, x1_min - np.asarray(X)[:, 0])


def igd(F, ref):
    F, ref = np.asarray(F), np.asarray(ref)
    d = np.sqrt(((ref[:, None, :] - F[None, :, :]) ** 2).sum(axis=2))
    return float(d.min(axis=1).mean())


def zdt1_front(n=500, x1_min=0.0):
    f1 = np.linspace(x1_min, 1.0, n)
    return np.stack([f1, 1 - np.sqrt(f1)], axis=1)
