"""NSGA-II (Deb, Pratap, Agarwal, Meyarivan 2002, IEEE TEC 6(2)) с ограничениями по правилам допустимости Деба:
допустимое решение лучше недопустимого; из двух недопустимых лучше то, у которого меньше суммарное нарушение CV;
из двух допустимых — по доминированию Парето. Все цели минимизируются. Геном — [0, 1]^d.

Тест — ZDT1 (Zitzler, Deb, Thiele 2000): фронт f2 = 1 − √f1; гиперобъём фронта с опорной точкой (1,1; 1,1) = 0,8767.
"""
import numpy as np

from . import archive as A


def dominates(fa, fb, cva=0.0, cvb=0.0):
    if cva <= 0.0 and cvb > 0.0:
        return True
    if cva > 0.0 and cvb > 0.0:
        return cva < cvb
    if cva > 0.0 and cvb <= 0.0:
        return False
    return bool(np.all(fa <= fb) and np.any(fa < fb))


def nondominated_sort(F, CV=None):
    """быстрая сортировка по фронтам (векторно, O(M·N²) по памяти N²). Возврат: ranks (N,), fronts — список массивов индексов"""
    F = np.asarray(F, float)
    N = len(F)
    CV = np.zeros(N) if CV is None else np.maximum(np.asarray(CV, float), 0.0)
    le = np.all(F[:, None, :] <= F[None, :, :], axis=2)
    lt = np.any(F[:, None, :] < F[None, :, :], axis=2)
    dom = le & lt
    feas = CV <= 0.0
    dom = np.where(feas[:, None] & feas[None, :], dom, False)
    dom |= feas[:, None] & ~feas[None, :]
    both_inf = ~feas[:, None] & ~feas[None, :]
    dom |= both_inf & (CV[:, None] < CV[None, :])
    n_dom = dom.sum(axis=0)              # сколькими доминируется j
    ranks = np.full(N, -1, int)
    fronts = []
    cur = np.nonzero(n_dom == 0)[0]
    r = 0
    n_dom = n_dom.copy()
    while len(cur):
        ranks[cur] = r
        fronts.append(cur)
        n_dom = n_dom - dom[cur].sum(axis=0)
        n_dom[ranks >= 0] = -1
        cur = np.nonzero(n_dom == 0)[0]
        r += 1
    return ranks, fronts


def crowding(F):
    F = np.asarray(F, float)
    n, m = F.shape
    d = np.zeros(n)
    if n <= 2:
        return np.full(n, np.inf)
    for j in range(m):
        o = np.argsort(F[:, j], kind="stable")
        fj = F[o, j]
        span = fj[-1] - fj[0]
        d[o[0]] = d[o[-1]] = np.inf
        if span > 0:
            d[o[1:-1]] += (fj[2:] - fj[:-2]) / span
    return d


def _select(F, CV, n):
    ranks, fronts = nondominated_sort(F, CV)
    keep, cd_all = [], np.zeros(len(F))
    for fr in fronts:
        cd = crowding(F[fr])
        cd_all[fr] = cd
        if len(keep) + len(fr) <= n:
            keep.extend(fr.tolist())
        else:
            o = np.lexsort((fr, -cd))       # по убыванию скученности, при равенстве — по индексу (детерминизм)
            keep.extend(fr[o[: n - len(keep)]].tolist())
            break
    keep = np.array(keep, int)
    return keep, ranks, cd_all


def _tournament(rng, ranks, cd, k):
    n = len(ranks)
    a, b = rng.integers(0, n, k), rng.integers(0, n, k)
    better_a = (ranks[a] < ranks[b]) | ((ranks[a] == ranks[b]) & (cd[a] > cd[b]))
    return np.where(better_a, a, b)


def nsga2(evaluate, dim, pop_size, n_gen, rng, init=None, eta_c=15.0, eta_m=20.0, pc=0.9, pm=None, callback=None):
    """evaluate(X (N, dim)) → (F (N, M), CV (N,)). init — начальные геномы (дополняются случайными).
    Возврат: dict(X, F, CV, rank, history, evals)."""
    X = rng.random((pop_size, dim))
    if init is not None and len(init):
        k = min(len(init), pop_size)
        X[:k] = np.clip(np.asarray(init, float)[:k], 0.0, 1.0)
    F, CV = evaluate(X)
    F, CV = np.asarray(F, float), np.asarray(CV, float)
    evals = len(X)
    hist = []
    keep, ranks, cd = _select(F, CV, pop_size)
    X, F, CV = X[keep], F[keep], CV[keep]
    ranks, _ = nondominated_sort(F, CV)
    cd = np.zeros(len(F))
    for r in range(ranks.max() + 1):
        idx = np.nonzero(ranks == r)[0]
        cd[idx] = crowding(F[idx])
    for g in range(n_gen):
        par = _tournament(rng, ranks, cd, pop_size)
        kids = []
        for i in range(0, pop_size, 2):
            p1, p2 = X[par[i]], X[par[(i + 1) % pop_size]]
            c1, c2 = A.sbx(rng, p1, p2, eta_c, pc)
            kids.append(A.poly_mutation(rng, c1, eta_m, pm))
            kids.append(A.poly_mutation(rng, c2, eta_m, pm))
        Q = np.stack(kids[:pop_size])
        FQ, CVQ = evaluate(Q)
        evals += len(Q)
        XX = np.vstack([X, Q]); FF = np.vstack([F, np.asarray(FQ, float)]); CC = np.concatenate([CV, np.asarray(CVQ, float)])
        keep, _, _ = _select(FF, CC, pop_size)
        X, F, CV = XX[keep], FF[keep], CC[keep]
        ranks, _ = nondominated_sort(F, CV)
        cd = np.zeros(len(F))
        for r in range(ranks.max() + 1):
            idx = np.nonzero(ranks == r)[0]
            cd[idx] = crowding(F[idx])
        rec = dict(gen=g + 1, evals=evals, n_front0=int(np.sum(ranks == 0)), n_feasible=int(np.sum(CV <= 0)))
        hist.append(rec)
        if callback:
            callback(rec, X, F, CV)
    return dict(X=X, F=F, CV=CV, rank=ranks, history=hist, evals=evals)


def hv2d(F, ref):
    """гиперобъём двух целей (минимизация) относительно опорной точки ref"""
    F = np.asarray(F, float)
    F = F[np.all(F < np.asarray(ref), axis=1)]
    if not len(F):
        return 0.0
    o = np.argsort(F[:, 0], kind="stable")
    F = F[o]
    hv, best2 = 0.0, ref[1]
    for i in range(len(F)):
        if F[i, 1] < best2:
            x_next = ref[0]
            hv += (x_next - F[i, 0]) * (best2 - F[i, 1])
            best2 = F[i, 1]
    return float(hv)


def hv_mc(F, ref, ideal, n, rng):
    """гиперобъём многих целей методом Монте-Карло в коробке [ideal, ref] (оценка, для отчёта)"""
    F = np.asarray(F, float)
    ref, ideal = np.asarray(ref, float), np.asarray(ideal, float)
    P = ideal + rng.random((n, len(ref))) * (ref - ideal)
    dom = np.zeros(n, bool)
    for f in F:
        dom |= np.all(P >= f, axis=1)
    return float(dom.mean() * np.prod(ref - ideal))


def zdt1(X):
    X = np.asarray(X, float)
    f1 = X[:, 0]
    g = 1.0 + 9.0 * X[:, 1:].mean(axis=1)
    f2 = g * (1.0 - np.sqrt(f1 / g))
    return np.stack([f1, f2], axis=1), np.zeros(len(X))


def igd(F, front):
    F, front = np.asarray(F, float), np.asarray(front, float)
    d = np.sqrt(((front[:, None, :] - F[None, :, :]) ** 2).sum(axis=2)).min(axis=1)
    return float(d.mean())
