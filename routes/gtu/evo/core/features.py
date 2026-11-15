"""Признаки «эффектности» формы. Красоту одним числом не измерить: здесь — измеримые признаки, а решает владелец по
галерее рендеров (его выбор — родители следующего поколения).

  * новизна — среднее расстояние до k ближайших в архиве (Lehman & Stanley 2011, поиск новизны) и наименьшее
    расстояние до типовых форм изготовителей (OEM) в стандартизованном пространстве признаков;
  * гладкость — скачки кривизны, энергия изгиба ∫κ'², число перегибов;
  * сложность силуэта — коробочная размерность контура, дефицит выпуклости, число дыр;
  * симметрия — действующий порядок вращательной симметрии по угловому распределению (Фурье) и зеркальность.
"""
import numpy as np

from . import shapes as S

_trapz = getattr(np, "trapezoid", None) or getattr(np, "trapz")   # numpy 2.x и 1.x (контейнер ubuntu:24.04)


def standardize(F, mu=None, sd=None):
    F = np.asarray(F, float)
    mu = F.mean(axis=0) if mu is None else mu
    sd = F.std(axis=0) if sd is None else sd
    return (F - mu) / np.where(sd > 0, sd, 1.0), mu, sd


def knn_novelty(Fz, pool_z, k=15, exclude_self=True):
    """средняя дистанция до k ближайших (Fz (n, m) против pool_z (p, m))"""
    D = np.sqrt(((Fz[:, None, :] - pool_z[None, :, :]) ** 2).sum(axis=2))
    if exclude_self:
        D = np.where(D < 1e-12, np.inf, D)
    k = max(1, min(k, D.shape[1] - (1 if exclude_self else 0)))
    Ds = np.sort(D, axis=1)[:, :k]
    Ds = np.where(np.isfinite(Ds), Ds, np.nan)
    return np.nanmean(Ds, axis=1)


def min_distance(Fz, ref_z):
    D = np.sqrt(((Fz[:, None, :] - ref_z[None, :, :]) ** 2).sum(axis=2))
    return D.min(axis=1)


def curve_smoothness(P):
    k = S.curvature_2d(P)
    s = np.r_[0.0, np.cumsum(np.linalg.norm(np.diff(P, axis=0), axis=1))]
    L = max(s[-1], 1e-12)
    dk = np.gradient(k, s + 1e-15 * np.arange(len(s)))
    sign = np.sign(k[np.abs(k) > 1e-6 / L])
    return dict(kappa_jump=float(np.max(np.abs(np.diff(k))) * L), fairness=float(_trapz(dk ** 2, s) * L ** 3),
                inflections=int(np.sum(sign[1:] != sign[:-1])) if len(sign) > 1 else 0, turning=S.turning_total(P))


def box_dimension(mask, sizes=(1, 2, 4, 8, 16)):
    """коробочная размерность границы бинарной картинки (контур — пиксели маски с соседом вне маски)"""
    m = np.asarray(mask, bool)
    edge = m & ~(np.roll(m, 1, 0) & np.roll(m, -1, 0) & np.roll(m, 1, 1) & np.roll(m, -1, 1))
    ns, ss = [], []
    for s in sizes:
        if s >= min(m.shape):
            break
        h, w = (m.shape[0] // s) * s, (m.shape[1] // s) * s
        b = edge[:h, :w].reshape(h // s, s, w // s, s).any(axis=(1, 3))
        if b.sum() > 0:
            ns.append(b.sum()); ss.append(s)
    if len(ns) < 2:
        return 1.0
    p = np.polyfit(np.log(1.0 / np.asarray(ss)), np.log(ns), 1)
    return float(p[0])


def convexity_deficit(points2d):
    from scipy.spatial import ConvexHull
    P = np.asarray(points2d, float)
    hull = ConvexHull(P)
    x, y = P[:, 0], P[:, 1]
    area = 0.5 * abs(np.dot(x, np.roll(y, 1)) - np.dot(y, np.roll(x, 1)))
    return float(1.0 - area / hull.volume) if hull.volume > 0 else 0.0


def angular_symmetry(rho, tol=0.08, n_max=24, rel_amp_min=0.05):
    """rho — угловое распределение на равномерной сетке θ (n,) или (n_z, n). Возврат: N_eff — наибольший порядок N ≤ n_max,
    при котором сдвиг на 2π/N совпадает с исходным с невязкой ≤ tol (относительно размаха); 0 — осесимметрично
    (нет углового рисунка или его размах меньше rel_amp_min от среднего); mirror — лучшая зеркальная невязка;
    asym — 1 − 1/N_eff (0 для осесимметричного)."""
    R = np.atleast_2d(np.asarray(rho, float))
    n = R.shape[1]
    amp = np.sqrt(np.mean((R - R.mean(axis=1, keepdims=True)) ** 2))
    if amp < 1e-9 or amp < rel_amp_min * abs(float(R.mean())):      # рябь меньше 5 % среднего — узора нет
        return dict(N_eff=0, mirror=0.0, asym=0.0, axisym=True)
    if amp < 1e-9:
        return dict(N_eff=0, mirror=0.0, asym=0.0, axisym=True)
    best = 1
    for N in range(n_max, 1, -1):
        if n % N:
            sh = n / N
            Rs = np.array([np.interp((np.arange(n) - sh) % n, np.arange(n + 1), np.r_[r, r[0]]) for r in R])
        else:
            Rs = np.roll(R, n // N, axis=1)
        err = np.sqrt(np.mean((Rs - R) ** 2)) / amp
        if err <= tol:
            best = N
            break
    mir = np.inf
    for k in range(n):
        Rm = np.roll(R[:, ::-1], k, axis=1)
        mir = min(mir, np.sqrt(np.mean((Rm - R) ** 2)) / amp)
    return dict(N_eff=int(best), mirror=float(mir), asym=float(1.0 - 1.0 / best), axisym=False)
