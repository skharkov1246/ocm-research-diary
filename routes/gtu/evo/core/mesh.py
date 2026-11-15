"""Сетки поверхностей для рендера и проверок.

  * surface_nets — поверхность уровня φ = 0 по сетке (двойственный метод «surface nets», Gibson 1998): вершина в
    каждой ячейке со сменой знака — среднее точек пересечения рёбер; грань — на каждом ребре со сменой знака;
    ориентация — наружу (в сторону роста φ). Без внешних библиотек (в среде нет scikit-image);
  * revolve_heightfield — замкнутая оболочка вращения: наружная поверхность r_out(z, θ) (рельеф рёбер), внутренняя
    r_in(z, θ), торцевые кольца; ось — z;
  * tube — трубка вдоль ломаной (стойки, каналы, пружины);
  * объём и площадь по сетке, запись двоичного STL, слияние сеток.
"""
import struct

import numpy as np


def surface_nets(phi, origin, h):
    phi = np.asarray(phi, float)
    nx, ny, nz = phi.shape
    s = phi < 0
    org = np.asarray(origin, float)
    # пересечения на рёбрах вдоль осей
    with np.errstate(divide="ignore", invalid="ignore"):
        ex = s[:-1, :, :] != s[1:, :, :]
        tx = np.where(ex, phi[:-1] / (phi[:-1] - phi[1:]), 0.0)
        ey = s[:, :-1, :] != s[:, 1:, :]
        ty = np.where(ey, phi[:, :-1] / (phi[:, :-1] - phi[:, 1:]), 0.0)
        ez = s[:, :, :-1] != s[:, :, 1:]
        tz = np.where(ez, phi[:, :, :-1] / (phi[:, :, :-1] - phi[:, :, 1:]), 0.0)
    cx, cy, cz = nx - 1, ny - 1, nz - 1
    acc = np.zeros((cx, cy, cz, 3))
    cnt = np.zeros((cx, cy, cz))
    I, J, K = np.meshgrid(np.arange(cx), np.arange(cy), np.arange(cz), indexing="ij")
    for dj in (0, 1):
        for dk in (0, 1):
            m = ex[:, dj:dj + cy, dk:dk + cz]
            t = tx[:, dj:dj + cy, dk:dk + cz]
            acc[..., 0] += np.where(m, I + t, 0); acc[..., 1] += np.where(m, J + dj, 0); acc[..., 2] += np.where(m, K + dk, 0)
            cnt += m
    for di in (0, 1):
        for dk in (0, 1):
            m = ey[di:di + cx, :, dk:dk + cz]
            t = ty[di:di + cx, :, dk:dk + cz]
            acc[..., 0] += np.where(m, I + di, 0); acc[..., 1] += np.where(m, J + t, 0); acc[..., 2] += np.where(m, K + dk, 0)
            cnt += m
    for di in (0, 1):
        for dj in (0, 1):
            m = ez[di:di + cx, dj:dj + cy, :]
            t = tz[di:di + cx, dj:dj + cy, :]
            acc[..., 0] += np.where(m, I + di, 0); acc[..., 1] += np.where(m, J + dj, 0); acc[..., 2] += np.where(m, K + t, 0)
            cnt += m
    active = cnt > 0
    vid = -np.ones((cx, cy, cz), np.int64)
    vid[active] = np.arange(int(active.sum()))
    V = org + h * (acc[active] / cnt[active][:, None])
    quads = []
    # рёбра x: соседние ячейки (i, j-1..j, k-1..k)
    ii, jj, kk = np.nonzero(ex[:, 1:-1, 1:-1]); jj += 1; kk += 1
    q = np.stack([vid[ii, jj - 1, kk - 1], vid[ii, jj, kk - 1], vid[ii, jj, kk], vid[ii, jj - 1, kk]], axis=1)
    flip = ~s[ii, jj, kk]
    q[flip] = q[flip][:, ::-1]
    quads.append(q)
    ii, jj, kk = np.nonzero(ey[1:-1, :, 1:-1]); ii += 1; kk += 1
    q = np.stack([vid[ii - 1, jj, kk - 1], vid[ii - 1, jj, kk], vid[ii, jj, kk], vid[ii, jj, kk - 1]], axis=1)
    flip = ~s[ii, jj, kk]
    q[flip] = q[flip][:, ::-1]
    quads.append(q)
    ii, jj, kk = np.nonzero(ez[1:-1, 1:-1, :]); ii += 1; jj += 1
    q = np.stack([vid[ii - 1, jj - 1, kk], vid[ii, jj - 1, kk], vid[ii, jj, kk], vid[ii - 1, jj, kk]], axis=1)
    flip = ~s[ii, jj, kk]
    q[flip] = q[flip][:, ::-1]
    quads.append(q)
    Q = np.vstack(quads)
    Q = Q[np.all(Q >= 0, axis=1)]
    F = np.vstack([Q[:, [0, 1, 2]], Q[:, [0, 2, 3]]])
    return V, F


def revolve_heightfield(z, theta, r_out, r_in):
    """z (nz,), θ (nt,) — периодическая сетка без повтора 2π; r_out, r_in — (nz, nt) или r_in (nz,).
    Возврат (V, F): замкнутая поверхность, нормали наружу из тела."""
    z, theta = np.asarray(z, float), np.asarray(theta, float)
    nz, nt = len(z), len(theta)
    r_out = np.broadcast_to(np.asarray(r_out, float), (nz, nt))
    r_in = np.asarray(r_in, float)
    r_in = np.broadcast_to(r_in[:, None] if r_in.ndim == 1 else r_in, (nz, nt))
    c, s_ = np.cos(theta), np.sin(theta)

    def grid(r):
        return np.stack([r * c[None, :], r * s_[None, :], np.broadcast_to(z[:, None], (nz, nt))], axis=-1).reshape(-1, 3)
    Vo, Vi = grid(r_out), grid(r_in)
    V = np.vstack([Vo, Vi])
    off = nz * nt
    iz, it = np.meshgrid(np.arange(nz - 1), np.arange(nt), indexing="ij")
    a = (iz * nt + it).ravel(); b = (iz * nt + (it + 1) % nt).ravel()
    c_ = ((iz + 1) * nt + (it + 1) % nt).ravel(); d = ((iz + 1) * nt + it).ravel()
    # наружная: нормаль от оси ((b − a) × (c − b) = θ̂ × ẑ = r̂); внутренняя — к оси
    Fo = np.vstack([np.stack([a, b, c_], 1), np.stack([a, c_, d], 1)])
    Fi = np.vstack([np.stack([a, d, c_], 1), np.stack([a, c_, b], 1)]) + off
    it1 = np.arange(nt); it2 = (it1 + 1) % nt
    # торец z0 (нормаль −z) и z1 (нормаль +z)
    o0, o1 = it1, it2
    i0, i1 = it1 + off, it2 + off
    F0 = np.vstack([np.stack([o0, i1, o1], 1), np.stack([o0, i0, i1], 1)])
    base = (nz - 1) * nt
    O0, O1 = base + it1, base + it2
    I0, I1 = O0 + off, O1 + off
    F1 = np.vstack([np.stack([O0, O1, I1], 1), np.stack([O0, I1, I0], 1)])
    return V, np.vstack([Fo, Fi, F0, F1])


def tube(P, radius, n_sides=10, caps=True):
    """трубка вдоль ломаной P (n, 3), радиус — число или массив (n,)"""
    P = np.asarray(P, float)
    n = len(P)
    R = np.broadcast_to(np.asarray(radius, float), (n,))
    T = np.gradient(P, axis=0)
    T /= np.linalg.norm(T, axis=1, keepdims=True) + 1e-300
    a = np.array([1.0, 0.0, 0.0]) if abs(T[0, 0]) < 0.9 else np.array([0.0, 1.0, 0.0])
    N = np.cross(T[0], a); N /= np.linalg.norm(N)
    frames = []
    for i in range(n):
        if i > 0:
            N = N - (N @ T[i]) * T[i]
            N /= np.linalg.norm(N) + 1e-300
        B = np.cross(T[i], N)
        frames.append((N.copy(), B))
    ang = np.linspace(0, 2 * np.pi, n_sides, endpoint=False)
    V = np.array([P[i] + R[i] * (np.cos(t) * frames[i][0] + np.sin(t) * frames[i][1]) for i in range(n) for t in ang])
    F = []
    for i in range(n - 1):
        for k in range(n_sides):
            a0, a1 = i * n_sides + k, i * n_sides + (k + 1) % n_sides
            b0, b1 = a0 + n_sides, a1 + n_sides
            F += [(a0, a1, b1), (a0, b1, b0)]
    if caps:
        c0 = len(V); c1 = c0 + 1
        V = np.vstack([V, P[0], P[-1]])
        for k in range(n_sides):
            F.append((c0, (k + 1) % n_sides, k))
            F.append((c1, (n - 1) * n_sides + k, (n - 1) * n_sides + (k + 1) % n_sides))
    return V, np.asarray(F, np.int64)


def merge(meshes):
    Vs, Fs, off = [], [], 0
    for V, F in meshes:
        Vs.append(np.asarray(V, float)); Fs.append(np.asarray(F, np.int64) + off); off += len(V)
    return np.vstack(Vs), np.vstack(Fs)


def volume(V, F):
    a, b, c = V[F[:, 0]], V[F[:, 1]], V[F[:, 2]]
    return float(np.einsum("ij,ij->i", a, np.cross(b, c)).sum() / 6.0)


def area(V, F):
    a, b, c = V[F[:, 0]], V[F[:, 1]], V[F[:, 2]]
    return float(0.5 * np.linalg.norm(np.cross(b - a, c - a), axis=1).sum())


def write_stl(path, V, F):
    V = np.asarray(V, np.float32); F = np.asarray(F)
    a, b, c = V[F[:, 0]], V[F[:, 1]], V[F[:, 2]]
    n = np.cross(b - a, c - a)
    n /= np.linalg.norm(n, axis=1, keepdims=True) + 1e-30
    with open(path, "wb") as f:
        f.write(b"evo core stl".ljust(80, b" "))
        f.write(struct.pack("<I", len(F)))
        rec = np.zeros(len(F), dtype=[("n", "<f4", 3), ("a", "<f4", 3), ("b", "<f4", 3), ("c", "<f4", 3), ("x", "<u2")])
        rec["n"], rec["a"], rec["b"], rec["c"] = n, a, b, c
        f.write(rec.tobytes())


def decimate_grid(V, F, keep_every=1):
    return V, F[::max(1, keep_every)]
