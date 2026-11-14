#!/usr/bin/env python3
"""Решётка TPMS (гироид, листовой вариант) для полых печатных деталей: геометрия и численная гомогенизация.

Гироид (Schoen 1970): g(x,y,z) = sin kx·cos ky + sin ky·cos kz + sin kz·cos kx, k = 2π/L;
листовой гироид — металл там, где |g| ≤ t. Доля объёма ρ*(t) — счётом на мелкой решётке.
Schwarz-D: g = cos kx cos ky cos kz − sin kx sin ky sin kz (для сравнения, только доля объёма).

Гомогенизация — по Andreassen & Andreasen (2014): ячейка N³ вокселей, периодические условия, шесть единичных
макродеформаций, C_H = (1/V)·Σ (u0 − χ)ᵀ k_e (u0 − χ); пустые воксели исключены (лист гироида связен).
Прочность решётки (начало текучести, оценка снизу): при макронапряжении единичной эквивалентной величины
ищется наибольшее местное σ_экв по вокселям; для устойчивости к ступенькам вокселей берётся 99,5-й процентиль,
по пяти видам нагружения (одноосное [100] и [111], двухосное, сдвиг, всестороннее растяжение) — худший.
Выход: E*/E(ρ*), S(ρ*) = σ*/σ(ρ*), коэффициент анизотропии Зенера, толщина листа при размере ячейки L.
"""
import json
import math
import os
import sys
import time

import numpy as np
import scipy.sparse as sp

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import fe3d  # noqa: E402

RESULTS = os.path.join(HERE, "lattice_results.json")


def gyroid(x, y, z, L=1.0):
    k = 2 * math.pi / L
    return np.sin(k * x) * np.cos(k * y) + np.sin(k * y) * np.cos(k * z) + np.sin(k * z) * np.cos(k * x)


def schwarz_d(x, y, z, L=1.0):
    k = 2 * math.pi / L
    return np.cos(k * x) * np.cos(k * y) * np.cos(k * z) - np.sin(k * x) * np.sin(k * y) * np.sin(k * z)


_TAB = None


def _fraction_table(n=128, kind="gyroid"):
    s = (np.arange(n) + 0.5) / n
    X, Y, Z = np.meshgrid(s, s, s, indexing="ij")
    g = np.abs(gyroid(X, Y, Z) if kind == "gyroid" else schwarz_d(X, Y, Z)).ravel()
    g.sort()
    ts = np.linspace(0.0, g[-1], 400)
    fr = np.searchsorted(g, ts, side="right") / g.size
    return ts, fr


def sheet_threshold(rho, kind="gyroid"):
    """t такое, что доля |g| ≤ t равна rho (таблица на решётке 128³)."""
    global _TAB
    if rho >= 0.9999:
        return 10.0                                # сплошной металл (|g| ≤ 1,5 везде)
    if _TAB is None or _TAB[0] != kind:
        _TAB = (kind,) + _fraction_table(kind=kind)
    _, ts, fr = _TAB
    return float(np.interp(rho, fr, ts))


def voxel_fraction(t, N, kind="gyroid"):
    s = (np.arange(N) + 0.5) / N
    X, Y, Z = np.meshgrid(s, s, s, indexing="ij")
    g = gyroid(X, Y, Z) if kind == "gyroid" else schwarz_d(X, Y, Z)
    return float(np.mean(np.abs(g) <= t))


def sheet_thickness(rho, L, N=160):
    """средняя толщина листа: объём металла / (2 · площадь срединной поверхности g = 0); площадь — по маршу кубов"""
    from skimage import measure
    s = (np.arange(N + 1)) / N
    X, Y, Z = np.meshgrid(s, s, s, indexing="ij")
    g = gyroid(X, Y, Z)
    v, f, _, _ = measure.marching_cubes(g, 0.0, spacing=(1.0 / N,) * 3)
    A = measure.mesh_surface_area(v, f)          # на единичную ячейку
    return rho * L / A, A


# =========================================================================== гомогенизация
def homogenize(rho, N=48, nu=0.3, kind="gyroid"):
    t0 = time.time()
    t = sheet_threshold(rho, kind)
    s = (np.arange(N) + 0.5) / N
    X, Y, Z = np.meshgrid(s, s, s, indexing="ij")
    g = gyroid(X, Y, Z) if kind == "gyroid" else schwarz_d(X, Y, Z)
    solid = (np.abs(g) <= t)
    h = 1.0 / N
    Ke, Bc, _, _, D = fe3d.hex_matrices(h, nu, 1.0)
    I, J, K_ = np.nonzero(solid)
    ne = len(I)
    loc = NODE = fe3d.NODE_XYZ.astype(int)
    ni = (I[:, None] + loc[None, :, 0]) % N
    nj = (J[:, None] + loc[None, :, 1]) % N
    nk = (K_[:, None] + loc[None, :, 2]) % N
    gnode = (ni * N + nj) * N + nk
    used, inv = np.unique(gnode, return_inverse=True)
    enode = inv.reshape(ne, 8)
    nn = len(used)
    edof = np.empty((ne, 24), int)
    for d in range(3):
        edof[:, d::3] = 3 * enode + d
    ii, jj = np.meshgrid(np.arange(24), np.arange(24), indexing="ij")
    rows = edof[:, ii].ravel(); cols = edof[:, jj].ravel()
    data = np.tile(Ke.ravel(), ne)
    Kg = sp.csr_matrix((data, (rows, cols)), shape=(3 * nn, 3 * nn))
    # u0 для единичных деформаций в узлах элемента (локальные координаты узла, h)
    xyz = NODE * h
    U0 = np.zeros((24, 6))
    eps_unit = np.eye(6)
    for c in range(6):
        e = eps_unit[c]
        Em = np.array([[e[0], e[5] / 2, e[4] / 2], [e[5] / 2, e[1], e[3] / 2], [e[4] / 2, e[3] / 2, e[2]]])
        U0[:, c] = (xyz @ Em.T).ravel()
    fe = Ke @ U0                                   # (24, 6) одинаково для всех элементов
    F = np.zeros((3 * nn, 6))
    for c in range(6):
        np.add.at(F[:, c], edof, np.broadcast_to(fe[:, c], (ne, 24)))
    # закрепить узел 0 (сдвиг как целого)
    free = np.ones(3 * nn, bool); free[:3] = False
    fidx = np.nonzero(free)[0]
    chi = np.zeros((3 * nn, 6))
    chi[fidx] = fe3d.solve_spd(Kg[fidx][:, fidx], F[fidx])
    CH = np.zeros((6, 6))
    Ue = chi[edof]                                 # (ne, 24, 6)
    W = U0[None, :, :] - Ue
    KW = np.einsum("ij,ejc->eic", Ke, W)
    CH = np.einsum("eic,eid->cd", W, KW)           # объём ячейки = 1
    # местные деформации в центре вокселя для каждой единичной макродеформации
    epsl = np.einsum("ij,ejc->eic", Bc, W)         # (ne, 6, 6)
    # макронапряжения, приведённые к σ_экв = 1
    s3 = 1 / math.sqrt(3)
    loads = {"uniax100": np.array([1, 0, 0, 0, 0, 0.0]), "uniax111": None,
             "biax": np.array([1, 1, 0, 0, 0, 0.0]), "shear": np.array([0, 0, 0, 0, 0, 1 / math.sqrt(3)]),
             "hydro": np.array([1, 1, 1, 0, 0, 0.0])}
    # [111]: σ = n⊗n, n = (1,1,1)/√3 -> σij = 1/3; в Фойгте [1/3,1/3,1/3,1/3,1/3,1/3]
    loads["uniax111"] = np.array([1, 1, 1, 1, 1, 1.0]) / 3
    Sinv = np.linalg.inv(CH)
    ratio = {}
    for name, sv in loads.items():
        vm_m = _vm6(sv)
        if vm_m < 1e-12:      # всестороннее растяжение: σ_экв = 0, нормируем по среднему напряжению
            norm = abs(sv[:3].mean())
        else:
            norm = vm_m
        sv = sv / norm
        eb = Sinv @ sv
        el = np.einsum("eic,c->ei", epsl, eb)
        sl = el @ D.T
        vm = _vm6(sl)
        ratio[name] = dict(max=float(1 / vm.max()), p995=float(1 / np.percentile(vm, 99.5)), p99=float(1 / np.percentile(vm, 99)))
    # изотропные приближения
    C11, C12, C44 = CH[0, 0], CH[0, 1], CH[3, 3]
    E100 = (C11 - C12) * (C11 + 2 * C12) / (C11 + C12)
    zener = 2 * C44 / (C11 - C12)
    # E в направлении [111] для кубической симметрии: 1/E111 = S11 − 2(S11 − S12 − S44/2)/3
    S11, S12, S44 = Sinv[0, 0], Sinv[0, 1], Sinv[3, 3]
    E111 = 1 / (S11 - 2 * (S11 - S12 - S44 / 2) / 3)
    Es = 1.0
    worst = min(v["p995"] for k, v in ratio.items() if k != "hydro")
    return dict(rho=rho, N=N, kind=kind, t=t, vox_fraction=float(solid.mean()), n_solid=int(ne), dofs=int(3 * nn),
                C11=float(C11), C12=float(C12), C44=float(C44), E100=float(E100 / Es), E111=float(E111 / Es),
                Emin_dir=float(min(E100, E111) / Es), zener=float(zener), strength_ratio=ratio,
                S_worst_p995=float(worst), seconds=round(time.time() - t0, 1))


def _vm6(s):
    s = np.asarray(s)
    return np.sqrt(np.maximum(0.5 * ((s[..., 0] - s[..., 1]) ** 2 + (s[..., 1] - s[..., 2]) ** 2 + (s[..., 2] - s[..., 0]) ** 2)
                              + 3 * (s[..., 3] ** 2 + s[..., 4] ** 2 + s[..., 5] ** 2), 0))


# =========================================================================== закон для оптимизации
DEFAULT_RHO_MIN = 0.15


def law_from_points(rhos, E, S, rho_min=DEFAULT_RHO_MIN):
    """Монотонные сплайны PCHIP по log ρ — log E*/E и log S, с точкой (1, 1)."""
    from scipy.interpolate import PchipInterpolator
    r = np.r_[np.asarray(rhos, float), 1.0]
    lE = np.log(np.r_[np.asarray(E, float), 1.0]); lS = np.log(np.r_[np.asarray(S, float), 1.0])
    lr = np.log(r)
    fE = PchipInterpolator(lr, lE); fS = PchipInterpolator(lr, lS)
    dfE, dfS = fE.derivative(), fS.derivative()

    def E_(x):
        x = np.clip(np.asarray(x, float), r[0], 1.0)
        return np.exp(fE(np.log(x)))

    def dE_(x):
        x = np.clip(np.asarray(x, float), r[0], 1.0)
        return E_(x) * dfE(np.log(x)) / x

    def S_(x):
        x = np.clip(np.asarray(x, float), r[0], 1.0)
        return np.exp(fS(np.log(x)))

    def dS_(x):
        x = np.clip(np.asarray(x, float), r[0], 1.0)
        return S_(x) * dfS(np.log(x)) / x

    return dict(E=E_, dE=dE_, S=S_, dS=dS_, rho_min=rho_min, points=dict(rho=list(map(float, rhos)), E=list(map(float, E)),
                                                                          S=list(map(float, S))))


def default_law(rho_min=DEFAULT_RHO_MIN):
    """Закон из lattice_results.json (если посчитан), иначе — степенной закон Гибсона — Эшби для проверок
    (E*/E = ρ*², S = ρ*^1,5 — только для тестов кода, в результаты не идёт)."""
    if os.path.exists(RESULTS):
        with open(RESULTS, encoding="utf-8") as f:
            d = json.load(f)
        L = d["law"]
        return law_from_points(L["rho"], L["E"], L["S"], rho_min)
    rr = np.array([0.1, 0.2, 0.4, 0.7])
    return law_from_points(rr, rr ** 2, rr ** 1.5, rho_min)


def main():
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--N", type=int, default=48)
    ap.add_argument("--rhos", default="0.10,0.15,0.20,0.30,0.45,0.60,0.80")
    ap.add_argument("--conv", action="store_true", help="сходимость по сетке при ρ* = 0,2 (N = 32, 48, 64)")
    a = ap.parse_args()
    rhos = [float(v) for v in a.rhos.split(",")]
    out = dict(title="Листовой гироид: гомогенизация упругости и оценка прочности (свой 3D-КЭ)",
               method="Andreassen & Andreasen (2014), воксели hex8 N³, периодические условия; пустые воксели исключены; "
                      "прочность — начало текучести по 99,5-му процентилю местного σ_экв, худшее из пяти видов макронапряжения",
               software=dict(python=sys.version.split()[0], numpy=np.__version__, scipy=__import__("scipy").__version__),
               date="2026-11-14", inputs=dict(nu=0.3, N=a.N, rhos=rhos), points=[])
    for r in rhos:
        res = homogenize(r, a.N)
        print(json.dumps({k: res[k] for k in ("rho", "vox_fraction", "E100", "E111", "zener", "S_worst_p995", "seconds")}), flush=True)
        out["points"].append(res)
    if a.conv:
        out["convergence"] = [homogenize(0.2, n) for n in (32, 48, 64)]
    pts = out["points"]
    out["law"] = dict(rho=[p["rho"] for p in pts], E=[p["Emin_dir"] for p in pts], S=[p["S_worst_p995"] for p in pts],
                      note="E — наименьший модуль по направлениям [100]/[111] (консервативно); S — худший 99,5-й процентиль")
    # степенные подгонки на 0,1–0,45
    sel = [p for p in pts if p["rho"] <= 0.45]
    lr = np.log([p["rho"] for p in sel])
    out["fit_power"] = dict(
        E=dict(zip(("n", "lnC"), np.polyfit(lr, np.log([p["Emin_dir"] for p in sel]), 1).tolist())),
        S=dict(zip(("n", "lnC"), np.polyfit(lr, np.log([p["S_worst_p995"] for p in sel]), 1).tolist())),
        range="ρ* 0,10–0,45")
    th = {}
    for L in (4e-3, 6e-3, 8e-3):
        th[f"L{int(L * 1e3)}mm"] = {f"{r:.2f}": round(sheet_thickness(r, L)[0] * 1e3, 3) for r in (0.10, 0.15, 0.20, 0.30)}
    out["sheet_thickness_mm"] = th
    out["surface_area_per_cell_L2"] = round(sheet_thickness(0.2, 1.0)[1], 4)
    out["schwarz_d_threshold"] = {f"{r:.2f}": sheet_threshold(r, "schwarz_d") for r in (0.15, 0.3)}
    with open(RESULTS, "w", encoding="utf-8") as f:
        json.dump(out, f, ensure_ascii=False, indent=1)
    print("->", RESULTS)


if __name__ == "__main__":
    main()
