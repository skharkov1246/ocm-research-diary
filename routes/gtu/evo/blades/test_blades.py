#!/usr/bin/env python3
"""Проверки семейства «лопатки и колёса» и локального ядра поиска (pytest).

  python3 -m pytest -q test_blades.py            # около 1 мин на 4 ядрах
Рендер pyvista проверяется, если pyvista есть в окружении (иначе — matplotlib)."""
import json
import math
import os
import sys
import tempfile

import numpy as np
import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import core_local as C   # noqa: E402


# =========================================================================== ядро: представления формы


def test_bspline_partition_derivative_c2():
    N = C.bspline_basis(8, 3, np.linspace(0, 1, 101))
    assert np.abs(N.sum(1) - 1).max() < 1e-12
    ctrl = np.array([[0, 0], [1, 2], [2, -1], [3, 1], [4, 0], [5, 1.5]], float)
    t = np.linspace(0.05, 0.95, 37)
    d1 = C.bspline_eval(ctrl, t, 3, 1)
    h = 1e-6
    fd = (C.bspline_eval(ctrl, t + h, 3) - C.bspline_eval(ctrl, t - h, 3)) / (2 * h)
    assert np.abs(d1 - fd).max() < 1e-5
    f = C.spline_fairness(ctrl)
    assert f["c2_jump_rel"] < 1e-4          # кубический сплайн: кривизна непрерывна на узлах


def test_ffd_identity_and_fold():
    F = C.FFD([0, 0, 0], [1, 1, 1], (2, 2, 2))
    X = np.random.default_rng(0).random((50, 3))
    assert np.allclose(F.deform(X), X)
    F.dP[1, 1, 1] = [0.1, 0.05, 0.0]
    assert F.check_jacobian() > 0.5          # малая деформация — без складок
    F.dP[1, 1, 1] = [-4.0, 0.0, 0.0]
    assert F.check_jacobian() < 0            # большая — складка (якобиан меняет знак)


def test_sdf_smooth_union_and_min_wall():
    G, h = C.grid_points([-1, -1, -1], [1, 1, 1], 61)
    a = C.sd_sphere(G, [-0.3, 0, 0], 0.5)
    b = C.sd_sphere(G, [0.3, 0, 0], 0.5)
    u = C.smooth_union(a, b, 0.1)
    assert np.all(u <= np.minimum(a, b) + 1e-12)
    slab = np.maximum(np.abs(G[..., 2]) - 0.1, C.sd_sphere(G, [0, 0, 0], 0.9))
    t, n = C.min_wall_from_field(slab, h)
    assert abs(t - 0.2) <= 2 * h[0] and n > 0
    g = C.gyroid(G, 0.5)
    lat = C.tpms_graded(C.sd_sphere(G, [0, 0, 0], 0.8), g, 0.3, 0.1)
    assert (lat < 0).mean() < (C.sd_sphere(G, [0, 0, 0], 0.8) < 0).mean()   # решётка легче сплошного тела


def test_polyline_and_murray():
    sq = np.array([[0, 0], [1, 0], [1, 1], [0, 1]], float)
    bow = np.array([[0, 0], [1, 1], [1, 0], [0, 1]], float)
    assert not C.polyline_self_intersects(sq) and C.polyline_self_intersects(bow)
    rc = C.murray_children(1.0, 2, 3.0)
    assert abs(sum(r ** 3 for r in rc) - 1.0) < 1e-12
    assert C.lsystem("F", {"F": "F[+F]F"}, 2).count("F") == 9


# =========================================================================== ядро: поиск


def _rastrigin_eval(G):
    Y = G * 10.24 - 5.12
    f = -(10 * G.shape[1] + np.sum(Y ** 2 - 10 * np.cos(2 * np.pi * Y), 1))
    return [dict(f=fi, d=g[:2], feasible=True) for fi, g in zip(f, G)]


def test_map_elites_fills_archive():
    A = C.GridArchive([10, 10], [0, 0], [1, 1], 6)
    log = C.map_elites(_rastrigin_eval, A, 6, C.rng_of("t", "me"), 3000, batch=50)
    assert A.coverage() >= 0.95
    assert log[-1]["best"] > -40.0           # лучшее заметно лучше случайного (среднее случайное ≈ −110)
    B = C.CVTArchive(64, [0, 0], [1, 1], 6, seed_key="t-cvt")
    C.map_elites(_rastrigin_eval, B, 6, C.rng_of("t", "me2"), 2000, batch=50)
    assert B.coverage() >= 0.9


def test_nsga2_zdt1():
    res = C.nsga2(lambda X: (C.zdt1(X), np.zeros(len(X)), [None] * len(X)), 30, C.rng_of("t", "zdt1"), pop=100, gens=250)
    f1 = np.linspace(0, 1, 300)
    ref = np.c_[f1, 1 - np.sqrt(f1)]
    F = res["F"][res["rank"] == 0]
    assert C.igd(F, ref) < 0.02
    assert F[:, 0].min() < 0.05 and F[:, 0].max() > 0.95     # фронт покрыт целиком


def test_cmaes_sphere_rosenbrock():
    x, f, n, _ = C.cmaes_minimize(lambda X: np.sum((X - 0.3) ** 2, 1), np.full(10, 0.8), 0.3, C.rng_of("t", "cma1"),
                                  max_evals=6000, ftarget=1e-10)
    assert f < 1e-8

    def rosen(X):
        Y = 4 * X - 2
        return np.sum(100 * (Y[:, 1:] - Y[:, :-1] ** 2) ** 2 + (1 - Y[:, :-1]) ** 2, 1)
    x, f, n, _ = C.cmaes_minimize(rosen, np.full(5, 0.2), 0.2, C.rng_of("t", "cma2"), max_evals=30000, ftarget=1e-10)
    assert f < 1e-6


def test_islands_partition_and_merge_completeness():
    """острова на тестовой функции: деление на части не меняет итог слияния; неполный набор — отказ"""
    def island(i, N):
        A = C.CVTArchive(32, [0, 0], [1, 1], 6, seed_key="t-isl")
        C.map_elites(_rastrigin_eval, A, 6, C.rng_of("t", "isl", i), 400, batch=20)
        return dict(island=i, config_hash="cfg", code_hash="code", n_evals=400, archive=A.to_dict(keep_meta=False))
    N = 6
    all_ = [island(i, N) for i in range(N)]
    parts = [[island(i, N) for i in C.islands_of_part(p, 3, N)] for p in range(3)]
    a1, r1 = C.merge_island_archives(all_, range(N), "cfg", "code")
    a2, r2 = C.merge_island_archives([d for p in parts for d in p][::-1], range(N), "cfg", "code")
    assert r1["complete"] and r2["complete"] and a1.digest() == a2.digest()
    assert sorted(i for p in range(3) for i in C.islands_of_part(p, 3, N)) == list(range(N))
    with pytest.raises(ValueError):
        C.merge_island_archives(all_[:-1], range(N), "cfg", "code")
    with pytest.raises(ValueError):
        C.merge_island_archives(all_ + [all_[0]], range(N), "cfg", "code")
    bad = dict(all_[1], config_hash="other")
    with pytest.raises(ValueError):
        C.merge_island_archives([all_[0], bad] + all_[2:], range(N), "cfg", "code")
    _, rp = C.merge_island_archives(all_[:-1], range(N), "cfg", "code", allow_partial=True)
    assert not rp["complete"] and rp["missing"] == [N - 1]


# =========================================================================== семейство


@pytest.fixture(scope="module")
def BE():
    import blade_eval
    return blade_eval


def test_disk1d_analytic_uniform():
    """сплошной диск постоянной толщины без нагрева: σθ(0) = σr(0) = (3 + ν)/8·ρω²R²"""
    import disk1d as D1
    m = D1.material("ti64")
    R = 0.35
    r = np.linspace(0, R, 401)
    res = D1.solve(r, np.full_like(r, 0.02), 1.0, 20.0, "ti64", 14950.0)
    om = 14950 * math.pi / 30
    s0 = (3 + m["nu"]) / 8 * m["density"] * om ** 2 * R ** 2
    assert abs(res["hoop_max_MPa"] * 1e6 / s0 - 1) < 2e-3
    assert res["equilibrium_rel"] < 1e-5


def test_disk1d_calibration_recorded():
    """калибровка по CalculiX есть, таблицы k(h/r) по ступеням не убывают, на калибровочных профилях в области
    применимости быстрая оценка с таблицей не дальше чем в 1,5 раза от CalculiX"""
    d = json.load(open(os.path.join(HERE, "disk1d_calibration.json"), encoding="utf-8"))
    for comp in ("comp1", "comp2"):
        k = np.array(d["table"][comp]["k"])
        assert np.all(np.diff(k) >= -1e-12) and k.min() >= 1.0
    tc = d["table_check"]
    assert 0.6 < tc["ratio_min"] and tc["ratio_max"] < 1.5
    assert abs(d["k_peak"] - 1.53) < 0.05


def test_section_props_rectangle():
    import blade_geom as BG
    b, h = 0.04, 0.01
    P = np.array([[0, 0], [b, 0], [b, h], [0, h]], float)
    p = BG.props_from_raw(BG.raw_moments(P))
    assert abs(p["A"] - b * h) < 1e-15
    assert abs(p["Ixx"] - b * h ** 3 / 12) / (b * h ** 3 / 12) < 1e-10
    assert abs(p["Iyy"] - h * b ** 3 / 12) / (h * b ** 3 / 12) < 1e-10


def test_turbine_root_stress_analytic(BE):
    """лопатка без сужения, полости и смещений: σ_c = ρω²·(r_t² − r_h²)/2 (среднее в корне)"""
    import blade_geom as BG
    st = dict(Dm=0.5365, H2=0.05, H3=0.05, phi=0.5, psi=1.2, chord_R=0.04, zt=60, m_dot=25.5)
    p = dict(vortex=-1.0, taper=1.0, tub_A=0.0, tub_n=4, lean=0.0, bow=0.0, sweep=0.0, hwall=0.0, tmax_c=0.2, route=1)
    bl = BG.turbine_rotor_blade(st, p, 14950.0)
    om = 14950 * math.pi / 30
    # площадь сечения меняется по размаху (углы закрутки) — сравниваем с интегралом по фактическим площадям
    A = np.array([q["A"] for q in bl["props"]])
    s_ref = 7916 * om ** 2 * np.trapezoid(A * bl["r"], bl["r"]) / A[0]
    assert abs(bl["sig_c_MPa"] * 1e6 / s_ref - 1) < 1e-9
    s_simple = 7916 * om ** 2 * (bl["rt"] ** 2 - bl["rh"] ** 2) / 2
    assert abs(bl["sig_c_MPa"] * 1e6 / s_simple - 1) < 0.15     # при почти постоянной площади — близко к формуле


def test_baseline_reproduces_frame3(BE):
    """облик frame3 в генах семейства: ступень 1 допустима, КПД уровня 1 в пределах 1 п. от frame3_point (0,844)"""
    r = BE.evaluate("comp1", BE.baseline_genes("comp1"), use_ref=False)
    assert r["feasible"]
    assert abs(r["meta"]["eta_level1"] - 0.844) < 0.01
    assert abs(r["meta"]["D2"] - 0.6946) / 0.6946 < 0.01


def test_random_shapes_smooth_and_checked(BE):
    """случайные геномы: у допустимых форм нет самопересечений, стенка не тоньше предела маршрута, FFD без складок,
    задняя поверхность монотонна; заведомо тонкая лопатка отвергается"""
    import blade_geom as BG
    rng = C.rng_of("t", "shapes")
    n_valid = 0
    for comp in ("comp1", "comp2"):
        for _ in range(25):
            x = rng.random(len(BG.genes_for(comp)))
            r = BE.evaluate(comp, x, want_geom=True)
            ch = r["meta"].get("checks")
            if not ch or not ch.get("valid"):
                continue
            n_valid += 1
            assert not ch["self_intersection"]
            assert ch["t_min_m"] >= BG.T_MIN[r["meta"]["route"]] - 1e-12
            assert ch["ffd_jac_min"] > 0
            if comp == "comp2":
                assert r["geom"]["disk"]["info"]["back_dxdr_max"] <= 1e-6
    assert n_valid >= 5
    p = BG.decode("comp1", BE.baseline_genes("comp1"))
    p["t_tip_k"] = 0.65
    p["route"] = 1                     # литьё: 1,5 мм — тонкая периферия не проходит
    r = BE.evaluate("comp1", BG.encode("comp1", p), want_geom=True)
    assert not r["meta"]["checks"]["t_min_ok"]


def test_island_repeatable_real(BE):
    """настоящий остров (ступень 1, 48 оценок) дважды — один и тот же архив"""
    import evo_blades as EB
    with tempfile.TemporaryDirectory() as td:
        a = EB.run_island(("comp1", 3, 48, "mini", td, None))
        b = EB.run_island(("comp1", 3, 48, "mini", None, None))
    da = C.archive_from_dict(a["archive"], EB.make_archive("comp1", "mini").C).digest()
    db = C.archive_from_dict(b["archive"], EB.make_archive("comp1", "mini").C).digest()
    assert da == db


def test_render_monochrome():
    import blade_geom as BG
    prof = np.c_[np.linspace(0, 0.1, 12), 0.05 + 0.03 * np.sin(np.linspace(0, np.pi, 12))]
    V, T = BG.revolve(prof, 32)
    with tempfile.TemporaryDirectory() as td:
        paths, eng = C.render_views([(V, T, 0.7)], os.path.join(td, "v{}.png"), (200, 200), backend="matplotlib")
        from matplotlib import image as mpimg
        for p in paths:
            img = mpimg.imread(p)
            assert img.shape[0] >= 150
            rgb = img[..., :3]
            assert np.abs(rgb[..., 0] - rgb[..., 1]).max() < 1e-6 and np.abs(rgb[..., 1] - rgb[..., 2]).max() < 1e-6
            assert rgb.min() < 0.9                      # не пустой лист
        try:
            import pyvista  # noqa: F401
        except ImportError:
            return
        paths, eng = C.render_views([(V, T, 0.7)], os.path.join(td, "p{}.png"), (200, 200), backend="pyvista")
        assert eng == "pyvista" and all(os.path.getsize(p) > 1000 for p in paths)
