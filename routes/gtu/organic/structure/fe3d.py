#!/usr/bin/env python3
"""Трёхмерный воксельный КЭ (гексаэдр 8 узлов, интегрирование 2×2×2) — общий для гомогенизации решётки
(lattice.py) и трёхмерной оптимизации топологии (topo3d.py). Решатель — PARDISO (pypardiso), если установлен,
иначе SuperLU (scipy) или сопряжённые градиенты с алгебраической многосеточностью (pyamg) для больших задач."""
import math

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

GP = (-1.0 / math.sqrt(3.0), 1.0 / math.sqrt(3.0))
# локальные узлы куба: (0,0,0),(1,0,0),(1,1,0),(0,1,0),(0,0,1),(1,0,1),(1,1,1),(0,1,1)
NODE_XYZ = np.array([[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0], [0, 0, 1], [1, 0, 1], [1, 1, 1], [0, 1, 1]], float)


def D_iso(E, nu):
    c = E / ((1 + nu) * (1 - 2 * nu))
    D = np.zeros((6, 6))
    D[:3, :3] = c * nu
    np.fill_diagonal(D[:3, :3], c * (1 - nu))
    D[3:, 3:] = np.eye(3) * c * (1 - 2 * nu) / 2
    return D


def B_hex(xi, eta, zeta, h):
    s = 2 * NODE_XYZ - 1
    N = 0.125 * (1 + s[:, 0] * xi) * (1 + s[:, 1] * eta) * (1 + s[:, 2] * zeta)
    dN = np.stack([0.125 * s[:, 0] * (1 + s[:, 1] * eta) * (1 + s[:, 2] * zeta),
                   0.125 * s[:, 1] * (1 + s[:, 0] * xi) * (1 + s[:, 2] * zeta),
                   0.125 * s[:, 2] * (1 + s[:, 0] * xi) * (1 + s[:, 1] * eta)], axis=0) * 2 / h
    B = np.zeros((6, 24))
    B[0, 0::3] = dN[0]; B[1, 1::3] = dN[1]; B[2, 2::3] = dN[2]
    B[3, 1::3] = dN[2]; B[3, 2::3] = dN[1]
    B[4, 0::3] = dN[2]; B[4, 2::3] = dN[0]
    B[5, 0::3] = dN[1]; B[5, 1::3] = dN[0]
    return N, B


def hex_matrices(h, nu, E=1.0):
    """Ke (24×24), Bc в центре (6×24), fth (24) для единичной тепловой деформации, fb (24×3) — интеграл N (объёмная сила)."""
    D = D_iso(E, nu)
    Ke = np.zeros((24, 24)); fth = np.zeros(24); fb = np.zeros((24, 3))
    w = (h / 2) ** 3
    for a in GP:
        for b in GP:
            for c in GP:
                N, B = B_hex(a, b, c, h)
                Ke += B.T @ D @ B * w
                fth += B.T @ D @ np.array([1, 1, 1, 0, 0, 0.0]) * w
                for d in range(3):
                    fb[d::3, d] += N * w
    _, Bc = B_hex(0, 0, 0, h)
    return Ke, Bc, fth, fb, D


def solve_spd(K, F, tol=1e-8, B=None, method="auto"):
    """K — симметричная положительно определённая (CSR/CSC); F — (n,) или (n, k).
    method: auto — PARDISO (если есть) до 100 тыс. неизвестных, больше — AMG-CG (память: PARDISO на 3D-ячейке
    200 тыс. неизвестных занимал 7,8 ГБ и был убит OOM 2026-11-14); direct; amg."""
    n = K.shape[0]
    if method == "auto":
        method = "direct" if n < 100_000 else "amg"
    if method == "direct":
        try:
            import pypardiso
            return pypardiso.spsolve(K.tocsr(), F)
        except ImportError:
            pass
    if method == "direct" and n < 250_000:
        lu = spla.splu(K.tocsc())
        return lu.solve(F)
    import pyamg
    if B is None and n % 3 == 0:
        B = np.kron(np.ones((n // 3, 1)), np.eye(3))       # поступательные моды упругости (3 на узел)
    ml = pyamg.smoothed_aggregation_solver(K.tocsr(), B=B, max_coarse=2000)
    F2 = F.reshape(n, -1)
    X = np.zeros_like(F2)
    for k in range(F2.shape[1]):
        res = []
        X[:, k] = ml.solve(F2[:, k], tol=tol, accel="cg", maxiter=1000, residuals=res)
        rel = res[-1] / max(res[0], 1e-300)
        if rel > 1e3 * tol:
            raise RuntimeError(f"AMG-CG не сошёлся: невязка {rel:.2e}")
    return X.reshape(F.shape)


class Factor:
    """Факторизация для многих правых частей (PARDISO или SuperLU)."""

    def __init__(self, K):
        self.K = K.tocsr()
        try:
            import pypardiso
            self.ps = pypardiso.PyPardisoSolver(mtype=11)
            self.ps.factorize(self.K)
            self.kind = "pardiso"
        except ImportError:
            self.lu = spla.splu(K.tocsc())
            self.kind = "superlu"

    def solve(self, F):
        if self.kind == "pardiso":
            return self.ps.solve(self.K, np.ascontiguousarray(F))
        return self.lu.solve(F)

    def free(self):
        if self.kind == "pardiso":
            self.ps.free_memory(everything=True)
