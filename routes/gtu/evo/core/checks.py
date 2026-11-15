"""Проверки изготовимости и гладкости формы по полю SDF на сетке.

  * min_wall — наименьшая толщина стенки: евклидово преобразование расстояния (scipy.ndimage) внутри тела,
    гребень (локальные максимумы) — срединная поверхность; толщина = 2·EDT − h/2 (поправка на дискретизацию:
    у стенки из чётного числа вокселей точно 2·EDT, из нечётного 2·EDT − h; погрешность ±h/2);
  * overhang — доля площади поверхности, обращённой вниз круче допустимого угла к плите (правило печати LPBF:
    самонесущее нависание не круче 45° к плите — organic/print_review RULES overhang_deg, по выдержке поиска);
    площадь — через сглаженную дельта-функцию δ_ε(φ)|∇φ|;
  * powder_escape — замкнутые полости (компоненты пустоты, не касающиеся края сетки) и доля пустоты, куда шар
    диаметром d не проходит снаружи (морфологическое открытие, связанное с краем): правило выноса порошка «канал
    ≥ 3 мм, при длине > 50 мм — ≥ 5 мм» (print_review RULES depowder_long_mm, по выдержке поиска);
  * smooth_profile — непрерывность кривизны и отсутствие самопересечений контура.
"""
import numpy as np
from scipy import ndimage as ndi

from . import shapes as S


def voxelize(sdf, lo, hi, h, chunk_z=None):
    """φ на узлах сетки с шагом h в коробке [lo, hi] (точки — центры вокселей). sdf(P (..., 3)) → (...)."""
    lo, hi = np.asarray(lo, float), np.asarray(hi, float)
    n = np.maximum(2, np.round((hi - lo) / h).astype(int))
    ax = [lo[i] + (np.arange(n[i]) + 0.5) * h for i in range(3)]
    phi = np.empty(tuple(n), np.float32)
    cz = chunk_z or max(1, int(4e6 // (n[0] * n[1])))
    X, Y = np.meshgrid(ax[0], ax[1], indexing="ij")
    for k0 in range(0, n[2], cz):
        zs = ax[2][k0:k0 + cz]
        P = np.empty((n[0], n[1], len(zs), 3))
        P[..., 0] = X[..., None]; P[..., 1] = Y[..., None]; P[..., 2] = zs[None, None, :]
        phi[:, :, k0:k0 + cz] = sdf(P)
    return phi, ax


def min_wall(phi, h, t_min=None):
    solid = phi < 0
    if not solid.any():
        return dict(t_min_m=None, note="пустое тело")
    edt = ndi.distance_transform_edt(solid, sampling=h)
    mx = ndi.maximum_filter(edt, size=3, mode="constant")
    ridge = solid & (edt >= mx - 1e-12) & (edt >= 0.75 * h)
    # край сетки — не стенка
    ridge[[0, -1], :, :] = False; ridge[:, [0, -1], :] = False; ridge[:, :, [0, -1]] = False
    t = 2.0 * edt[ridge] - 0.5 * h
    if not len(t):
        return dict(t_min_m=float(h), n_ridge=0, note="тоньше одного вокселя")
    out = dict(t_min_m=float(t.min()), t_p05_m=float(np.percentile(t, 5)), t_median_m=float(np.median(t)), n_ridge=int(len(t)))
    if t_min is not None:
        out["frac_below"] = float(np.mean(t < t_min))
    return out


def _grad(phi, h):
    gx, gy, gz = np.gradient(phi.astype(float), h)
    return np.stack([gx, gy, gz], axis=-1)


def overhang(phi, h, build=(0.0, 0.0, 1.0), max_angle_deg=45.0, plate_tol=None):
    """доля площади, обращённой вниз круче max_angle к плите (нужны поддержки). Нижняя грань на плите не считается."""
    b = np.asarray(build, float); b /= np.linalg.norm(b)
    g = _grad(phi, h)
    gn = np.linalg.norm(g, axis=-1) + 1e-30
    eps = 1.5 * h
    delta = np.where(np.abs(phi) < eps, 0.5 / eps * (1 + np.cos(np.pi * phi / eps)), 0.0)
    w = delta * gn * h ** 3                           # элемент площади
    n = g / gn[..., None]
    down = -(n @ b)                                   # косинус угла нормали с «вниз»
    need = down > np.cos(np.radians(max_angle_deg))
    # высота вдоль направления построения: плита — наименьшая высота тела
    idx = np.indices(phi.shape).transpose(1, 2, 3, 0) * h
    height = idx @ b
    solid = phi < 0
    hmin = height[solid].min() if solid.any() else 0.0
    plate_tol = 2.0 * h if plate_tol is None else plate_tol
    on_plate = height <= hmin + plate_tol
    A = float(w.sum())
    A_need = float(w[need & ~on_plate].sum())
    return dict(area_m2=A, need_support_m2=A_need, frac=A_need / A if A > 0 else 0.0)


def powder_escape(phi, h, d_list=(3e-3, 5e-3)):
    """замкнутые полости и недоступная шару диаметра d доля пустоты (пустота — φ > 0, наружа — край сетки)"""
    void = phi > 0
    lab, n = ndi.label(void)
    border = np.unique(np.concatenate([lab[[0, -1], :, :].ravel(), lab[:, [0, -1], :].ravel(), lab[:, :, [0, -1]].ravel()]))
    border = border[border > 0]
    closed = np.setdiff1d(np.arange(1, n + 1), border)
    vol_void = float(void.sum()) * h ** 3
    vol_closed = float(np.isin(lab, closed).sum()) * h ** 3 if len(closed) else 0.0
    edt = ndi.distance_transform_edt(void, sampling=h)
    res = dict(n_closed=int(len(closed)), closed_volume_m3=vol_closed, void_volume_m3=vol_void, blocked={}, blocked_m3={})
    for d in d_list:
        core = edt >= 0.5 * d
        lc, _ = ndi.label(core)
        bc = np.unique(np.concatenate([lc[[0, -1], :, :].ravel(), lc[:, [0, -1], :].ravel(), lc[:, :, [0, -1]].ravel()]))
        bc = bc[bc > 0]
        reach_core = np.isin(lc, bc) if len(bc) else np.zeros_like(core)
        # область, заметаемая шаром: точки пустоты в пределах d/2 от доступных центров
        if reach_core.any():
            dist_to_core = ndi.distance_transform_edt(~reach_core, sampling=h)
            swept = void & (dist_to_core <= 0.5 * d + 0.5 * h)
        else:
            swept = np.zeros_like(void)
        res["blocked"][f"{d * 1e3:g}mm"] = float(1.0 - swept.sum() / max(void.sum(), 1))       # доля всей пустоты коробки
        res["blocked_m3"][f"{d * 1e3:g}mm"] = float((void.sum() - swept.sum()) * h ** 3)      # объём, куда шар не проходит
    return res


def smooth_profile(P, closed=False):
    """контур (N, 2): самопересечения, наибольший скачок кривизны между соседними точками (× длина контура),
    энергия изгиба ∫κ'² ds (меньше — «чище» линия)"""
    P = np.asarray(P, float)
    k = S.curvature_2d(P)
    s = np.r_[0.0, np.cumsum(np.linalg.norm(np.diff(P, axis=0), axis=1))]
    L = s[-1]
    dk = np.gradient(k, s + 1e-300 * np.arange(len(s)))
    return dict(self_intersections=S.polyline_self_intersections(P, closed), kappa_jump_max=float(np.max(np.abs(np.diff(k))) * L),
                fairness=float(np.trapezoid(dk ** 2, s) * L ** 3), turning=S.turning_total(P), length=float(L))
