#!/usr/bin/env python3
"""Проверки узла organic/structure (pytest или python3 test_structure.py).

1. Свой осесимметричный КЭ против аналитики тонкого вращающегося диска и параболического нагрева.
2. Сопряжённые производные p-нормы загрузки — против конечных разностей (группы «solid» и «lattice»).
3. MMA на выпуклой задаче с известным решением.
4. Гироид: доля объёма неявной поверхности против счёта вокселей; закон гомогенизации монотонен.
5. CLOAD для осесимметричных элементов ccx: соглашение (сила на 2°-сектор или на радиан) — сверка с давлением.
"""
import math
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import topo_axi as TA  # noqa: E402

STEEL = dict(name="сталь (условно, для проверки)", density=7850.0, nu=0.3, E_GPa=[[20, 210.0], [1000, 210.0]],
             alpha_mean_1e6=[[20, 12.0], [1000, 12.0]], uts_MPa=[[20, 600.0]], ys_MPa=[[20, 400.0]], _lmp=None)


def disk_model(R=0.2, t=0.004, nr=200, nz=4, T=lambda r: 20.0 + 0 * r):
    act = np.ones((nr, nz), bool)
    grp = np.zeros((nr, nz), int)
    m = TA.AxiModel(0.0, R, -t / 2, t / 2, nr, nz, act, grp, STEEL, T)
    m.fix(m.nodes_where(lambda r, z: r < 1e-12), 0)
    m.fix(m.nodes_where(lambda r, z: (r < 1e-12) & (np.abs(z) < 1e-12)), 1)
    m.rpm0 = 10000.0
    return m


def test_spinning_disk_analytic():
    R, nu, rho = 0.2, 0.3, 7850.0
    m = disk_model(R)
    w = 2 * math.pi * 10000 / 60
    allow = np.full(m.ne, 1e9)
    sol = m.solve(np.ones(m.ne), dict(P=10), [dict(rpm=10000.0, thermal=False, allow=allow)], need_sens=False)
    sig = sol["cases"][0]["sig"]
    c = (3 + nu) / 8 * rho * w * w * R * R
    # центр: первый столбец (r = dr/2): σθ ≈ σr ≈ c(1 − (1+3ν)/(3+ν)·(r/R)²)
    r1 = m.dr / 2
    an_t = c * (1 - (1 + 3 * nu) / (3 + nu) * (r1 / R) ** 2)
    fe_t = float(np.mean(sig[m.ei == 0, 2]))
    an_rim = (1 - nu) / 4 * rho * w * w * R * R
    rl = R - m.dr / 2
    an_rim_t = c * (1 - (1 + 3 * nu) / (3 + nu) * (rl / R) ** 2)
    fe_rim = float(np.mean(sig[m.ei == m.nr - 1, 2]))
    assert abs(fe_t / an_t - 1) < 2e-3, (fe_t, an_t)
    assert abs(fe_rim / an_rim_t - 1) < 5e-3, (fe_rim, an_rim_t, an_rim)
    return dict(center_fe=fe_t, center_an=an_t, rim_fe=fe_rim, rim_an=an_rim_t)


def test_thermal_parabolic():
    R, dT = 0.2, 200.0
    T = lambda r: 20.0 + dT * (np.asarray(r) / R) ** 2  # noqa: E731
    m = disk_model(R, T=T)
    allow = np.full(m.ne, 1e9)
    sol = m.solve(np.ones(m.ne), dict(P=10), [dict(rpm=0.0, thermal=True, allow=allow)], need_sens=False)
    sig = sol["cases"][0]["sig"]
    E, a = 210e9, 12e-6
    # плоское напряжённое: σθ(r) = αEΔT(1/4 − 3/4 (r/R)²)
    r1 = m.dr / 2; rl = R - m.dr / 2
    an0 = a * E * dT * (0.25 - 0.75 * (r1 / R) ** 2)
    anR = a * E * dT * (0.25 - 0.75 * (rl / R) ** 2)
    fe0 = float(np.mean(sig[m.ei == 0, 2])); feR = float(np.mean(sig[m.ei == m.nr - 1, 2]))
    assert abs(fe0 / an0 - 1) < 3e-3, (fe0, an0)
    assert abs(feR / anR - 1) < 3e-3, (feR, anR)
    return dict(center_fe=fe0, center_an=an0, rim_fe=feR, rim_an=anR)


def _small_mixed():
    nr, nz = 16, 8
    act = np.ones((nr, nz), bool); act[12:, 6:] = False
    grp = np.ones((nr, nz), int); grp[:, :2] = 2; grp[14:, :6] = 0
    T = lambda r: 300.0 + 400.0 * (np.asarray(r) / 0.16) ** 2  # noqa: E731
    m = TA.AxiModel(0.0, 0.16, -0.03, 0.03, nr, nz, act, grp, STEEL, T)
    m.fix(m.nodes_where(lambda r, z: r < 1e-12), 0)
    m.fix(m.nodes_where(lambda r, z: (r < 1e-12) & (np.abs(z) < 1e-12)), 1)
    m.rpm0 = 12000.0
    m.add_edge_traction_r(0.16, -0.03, 0.0, 50e6)
    import lattice as LT
    lat = LT.default_law()
    prm = dict(q=5.0, Emin=1e-6, qs=0.5, P=8, lattice=lat)
    return m, prm


def test_adjoint_vs_fd():
    m, prm = _small_mixed()
    rng = np.random.default_rng(1)
    xb = rng.uniform(0.3, 0.9, m.ne); xb[m.group == 0] = 1.0
    allow = np.full(m.ne, 300e6)
    cases = [dict(rpm=12000.0, thermal=True, allow=allow), dict(rpm=14400.0, thermal=False, allow=allow * 1.3)]
    sol = m.solve(xb, prm, cases)
    errs = []
    for e in rng.choice(np.nonzero(m.group != 0)[0], 8, replace=False):
        h = 1e-6
        xp = xb.copy(); xp[e] += h
        xm = xb.copy(); xm[e] -= h
        sp = m.solve(xp, prm, cases, need_sens=False); sm = m.solve(xm, prm, cases, need_sens=False)
        for k in range(2):
            fd = (sp["cases"][k]["PN"] - sm["cases"][k]["PN"]) / (2 * h)
            ad = sol["cases"][k]["dPN"][e]
            errs.append(abs(ad - fd) / max(abs(fd), 1e-6 * abs(sol["cases"][k]["PN"])))
    assert max(errs) < 1e-4, errs
    return dict(max_rel_err=max(errs), n=len(errs))


def test_mma_convex():
    # min (x1−2)² + (x2−1)²  при  x1 + x2 ≤ 2, 0 ≤ x ≤ 3 -> (1.5, 0.5)
    mma = TA.MMA(2, 1, xmin=0.0, xmax=3.0, move=0.5)
    x = np.array([0.2, 0.2])
    for _ in range(80):
        f0 = (x[0] - 2) ** 2 + (x[1] - 1) ** 2
        df0 = np.array([2 * (x[0] - 2), 2 * (x[1] - 1)])
        g = np.array([x[0] + x[1] - 2]); dg = np.array([[1.0, 1.0]])
        x, _ = mma.step(x, f0, df0, g, dg)
    assert np.allclose(x, [1.5, 0.5], atol=2e-3), x
    return dict(x=x.tolist())


def test_gyroid_fraction():
    import lattice as LT
    for rs in (0.15, 0.3):
        t = LT.sheet_threshold(rs)
        frac = LT.voxel_fraction(t, 96)
        assert abs(frac - rs) < 0.01, (rs, frac)
    law = LT.default_law()
    xs = np.linspace(law["rho_min"], 1.0, 50)
    assert np.all(np.diff(law["E"](xs)) > 0) and np.all(np.diff(law["S"](xs)) > 0)
    return dict(ok=True)


def test_ccx_cload_convention():
    import ccx_check as CC
    r = CC.cload_convention()
    assert r["ok"], r
    return r


if __name__ == "__main__":
    import json
    out = {}
    for name, fn in list(globals().items()):
        if name.startswith("test_") and callable(fn):
            try:
                out[name] = fn()
                print("PASS", name, json.dumps(out[name], ensure_ascii=False, default=float)[:300])
            except Exception as ex:  # noqa: BLE001
                print("FAIL", name, repr(ex)[:500])
                out[name] = dict(fail=repr(ex)[:500])
