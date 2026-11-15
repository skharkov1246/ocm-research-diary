#!/usr/bin/env python3
"""Проверки общего модуля evo/core (pytest или python3 test_core.py).

 1. Повторяемость: один сид — один архив; деление островов на части и число процессов итог не меняют.
 2. MAP-Elites на тестовой функции заполняет архив (≥ 95 % ячеек) и даёт больший QD-score, чем случайный поиск.
 3. Слияние архивов: неполное — ошибка, отчёт называет недостающие острова; чужая конфигурация и разные дубли ловятся.
 4. NSGA-II на ZDT1: гиперобъём ≥ 0,85 при оптимуме 0,8767, IGD ≤ 0,02.
 5. CMA-ES: сфера и функция Розенброка в границах [0, 1].
 6. Гладкость форм: B-сплайн-контуры без самопересечений, кривизна непрерывна (скачок много меньше, чем у ломаной
    контрольных точек), периодический сплайн гладок на шве; FFD — без складок при малых смещениях, складка
    ловится при больших; обратное отображение точно.
 7. Наименьшая стенка по SDF: пластина 12 мм на вокселях 1 мм — 12 ± 1 мм; нависание шара — 0,146 ± 0,03
    (аналитика: доля сферического сегмента (1 − cos 45°)/2); порошок: замкнутая полость ловится, отверстие 8 мм
    пропускает шар 3 мм и не пропускает 12 мм.
 8. Решётки: плотность листового гироида монотонна по толщине и совпадает со счётом вокселей; дерево Мюррея
    выполняет закон в каждом узле и доходит до всех целей; лопатка — без самопересечений, закрутка по закону.
 9. Сетка: объём surface nets шара и оболочки вращения — против аналитики (1 %).
10. Рендер: лист из трёх видов и миниатюра создаются, повторный рендер даёт те же пиксели.
"""
import json
import math
import os
import sys
import tempfile

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
from core import archive as A  # noqa: E402
from core import blade as BL  # noqa: E402
from core import checks as C  # noqa: E402
from core import cmaes as CM  # noqa: E402
from core import features as FT  # noqa: E402
from core import islands as IS  # noqa: E402
from core import mesh as M  # noqa: E402
from core import nsga2 as N2  # noqa: E402
from core import render as R  # noqa: E402
from core import shapes as S  # noqa: E402


# ----------------------------------------------------------------------------- тестовая функция QD
def rastrigin_eval(G):
    X = 10.24 * (np.asarray(G) - 0.5)
    f = -(10 * X.shape[1] + np.sum(X ** 2 - 10 * np.cos(2 * np.pi * X), axis=1))
    return [dict(fitness=float(fi), desc=g[:2].copy()) for fi, g in zip(f, G)]


GRID = A.GridSpec.uniform(["x0", "x1"], [10, 10], [0, 0], [1, 1])


def run_island(island, base_seed=7, n_evals=3000, dim=6, cfg_hash="t"):
    rng = IS.rng_for(base_seed, island)
    ar = A.Archive(GRID)
    A.map_elites(rastrigin_eval, dim, ar, rng, n_evals, batch=50, uid_prefix=f"i{island:03d}-")
    return dict(island=island, config_hash=cfg_hash, evals=n_evals, archive=ar.to_dict())


def test_repeatability_and_parts():
    a = run_island(3); b = run_island(3)
    assert json.dumps(a, sort_keys=True) == json.dumps(b, sort_keys=True)
    one = IS.run_part(run_island, 4, 1, 0, workers=1, n_evals=1200)
    p0 = IS.run_part(run_island, 4, 2, 0, workers=2, n_evals=1200)
    p1 = IS.run_part(run_island, 4, 2, 1, workers=1, n_evals=1200)
    assert [d["island"] for d in p0] == [0, 2] and [d["island"] for d in p1] == [1, 3]
    m1, r1 = IS.merge_archives(one, range(4), "t")
    m2, r2 = IS.merge_archives(p1 + p0, range(4), "t")
    assert json.dumps(m1.to_dict()["elites"], sort_keys=True) == json.dumps(m2.to_dict()["elites"], sort_keys=True)
    assert r1["complete"] and r2["complete"]
    seeds = {IS.island_seed(7, i) for i in range(1000)}
    assert len(seeds) == 1000
    return dict(filled=len(m1), coverage=m1.coverage(), evals=r1["evals"])


def test_map_elites_fills():
    rng = IS.rng_for(11, 0)
    ar = A.Archive(GRID)
    A.map_elites(rastrigin_eval, 6, ar, rng, 6000, batch=50)
    rnd = A.Archive(GRID)
    rr = np.random.Generator(np.random.PCG64(11))
    for k in range(0, 6000, 500):
        G = rr.random((500, 6))
        for g, r in zip(G, rastrigin_eval(G)):
            rnd.try_add(g, r["fitness"], r["desc"], uid=str(k))
    off = -200.0
    assert ar.coverage() >= 0.95, ar.coverage()
    assert ar.qd_score(off) > rnd.qd_score(off), (ar.qd_score(off), rnd.qd_score(off))
    return dict(coverage=ar.coverage(), qd=ar.qd_score(off), qd_random=rnd.qd_score(off), best=ar.best()["fitness"])


def test_merge_completeness():
    parts = [run_island(i, n_evals=300) for i in range(3)]
    try:
        IS.merge_archives(parts[:2], range(3), "t")
        raise AssertionError("неполное слияние прошло")
    except ValueError as ex:
        assert "2" in str(ex)
    m, rep = IS.merge_archives(parts[:2], range(3), "t", allow_partial=True)
    assert rep["missing"] == [2] and not rep["complete"] and rep["status"].startswith("НЕПОЛНОЕ")
    bad = dict(parts[2], config_hash="другая")
    _, rep2 = IS.merge_archives(parts[:2] + [bad], range(3), "t", allow_partial=True)
    assert rep2["config_mismatch"] == [2]
    _, rep3 = IS.merge_archives(parts + [parts[0]], range(3), "t")
    assert rep3["duplicates"] == [0]
    other = run_island(0, base_seed=99, n_evals=300)
    try:
        IS.merge_archives(parts + [other], range(3), "t")
        raise AssertionError("разные дубли прошли")
    except ValueError:
        pass
    return dict(missing=rep["missing"], mismatch=rep2["config_mismatch"])


def test_nsga2_zdt1():
    rng = np.random.Generator(np.random.PCG64(5))
    res = N2.nsga2(N2.zdt1, 30, 100, 200, rng)
    F = res["F"][res["rank"] == 0]
    hv = N2.hv2d(F, (1.1, 1.1))
    f1 = np.linspace(0, 1, 500)
    igd = N2.igd(F, np.stack([f1, 1 - np.sqrt(f1)], 1))
    assert hv >= 0.85, hv
    assert igd <= 0.02, igd
    return dict(hv=hv, hv_opt=0.8767, igd=igd, evals=res["evals"])


def test_cmaes():
    rng = np.random.Generator(np.random.PCG64(3))
    c = np.full(8, 0.3)
    r1 = CM.fmin(lambda x: float(np.sum((x - c) ** 2)), np.full(8, 0.7), 0.2, rng, 6000, ftarget=1e-10)
    assert r1["f"] < 1e-8, r1["f"]

    def rosen(u):
        x = 4 * u - 2
        return float(np.sum(100 * (x[1:] - x[:-1] ** 2) ** 2 + (1 - x[:-1]) ** 2))
    r2 = CM.fmin(rosen, np.full(4, 0.2), 0.15, rng, 20000, ftarget=1e-10)
    assert r2["f"] < 1e-6, r2["f"]
    # граница: минимум за пределами куба — решение на грани
    r3 = CM.fmin(lambda x: float(np.sum((x - 1.3) ** 2)), np.full(3, 0.5), 0.2, rng, 3000)
    assert np.all(r3["x"] > 0.999), r3["x"]
    return dict(sphere=r1["f"], sphere_evals=r1["evals"], rosen=r2["f"], rosen_evals=r2["evals"])


def test_shapes_smooth_and_ffd():
    rng = np.random.Generator(np.random.PCG64(1))
    worst = 0
    for _ in range(50):
        z = np.linspace(0, 1, 7)
        r = 0.5 + 0.15 * rng.random(7)
        ctrl = np.stack([z, r], 1)
        P = S.bspline_eval(ctrl, np.linspace(0, 1, 400))
        worst = max(worst, S.polyline_self_intersections(P))
        sm = C.smooth_profile(P)
        lin = C.smooth_profile(np.vstack([np.interp(np.linspace(0, 6, 400), np.arange(7), z),
                                          np.interp(np.linspace(0, 6, 400), np.arange(7), r)]).T)
        assert sm["kappa_jump_max"] < 0.2 * max(lin["kappa_jump_max"], 1e-9) or lin["kappa_jump_max"] < 1e-6
    assert worst == 0
    th = np.linspace(0, 1, 721)
    ctrl = np.stack([np.cos(2 * np.pi * np.arange(8) / 8) * (1 + 0.2 * rng.random(8)),
                     np.sin(2 * np.pi * np.arange(8) / 8) * (1 + 0.2 * rng.random(8))], 1)
    P = S.bspline_eval(ctrl, th, closed=True)
    assert np.allclose(P[0], P[-1], atol=1e-9)
    k = S.curvature_2d(np.vstack([P[-6:-1], P, P[1:6]]))[5:-5]
    assert abs(k[0] - k[-1]) < 0.05 * np.max(np.abs(k))
    assert S.polyline_self_intersections(P[:-1], closed=True) == 0
    ffd = S.BSplineFFD((0, 0, 0), (1, 1, 1), (4, 4, 4), disp=0.03 * rng.standard_normal((4, 4, 4, 3)))
    jmin = ffd.jacobian_min()
    assert jmin > 0
    Pq = rng.random((200, 3)) * 0.8 + 0.1
    err = np.abs(ffd.deform(ffd.inverse(Pq)) - Pq).max()
    assert err < 1e-6, err
    big = S.BSplineFFD((0, 0, 0), (1, 1, 1), (4, 4, 4), disp=0.6 * rng.standard_normal((4, 4, 4, 3)))
    assert big.jacobian_min() <= 0
    return dict(jac_min_small=jmin, jac_min_big=big.jacobian_min(), inverse_err=err)


def test_min_wall_overhang_powder():
    h = 1e-3
    t = 12e-3
    phi, ax = C.voxelize(lambda P: np.abs(P[..., 2]) - t / 2, (-0.02, -0.02, -0.02), (0.02, 0.02, 0.02), h)
    mw = C.min_wall(phi, h)
    assert abs(mw["t_median_m"] - t) <= h, mw
    h2 = 0.01
    phi, ax = C.voxelize(lambda P: S.sd_sphere(P, (0, 0, 0), 0.4), (-0.5,) * 3, (0.5,) * 3, h2)
    ov = C.overhang(phi, h2, plate_tol=-1.0)
    exact = (1 - math.cos(math.radians(45))) / 2
    assert abs(ov["frac"] - exact) < 0.03, (ov, exact)
    hs = 1e-3
    hollow = lambda P: S.shell(S.sd_sphere(P, (0, 0, 0), 0.02), 3e-3)  # noqa: E731
    phi, _ = C.voxelize(hollow, (-0.03,) * 3, (0.03,) * 3, hs)
    pe = C.powder_escape(phi, hs, d_list=(3e-3,))
    v_in = 4 / 3 * math.pi * (0.02 - 1.5e-3) ** 3
    assert pe["n_closed"] == 1 and abs(pe["blocked_m3"]["3mm"] / v_in - 1) < 0.05, (pe, v_in)
    holed = lambda P: S.subtract(hollow(P), S.sd_capsule(P, (0, 0, 0.0), (0, 0, 0.04), 4e-3))  # noqa: E731
    phi, _ = C.voxelize(holed, (-0.03,) * 3, (0.03,) * 3, hs)
    pe2 = C.powder_escape(phi, hs, d_list=(3e-3, 12e-3))
    assert pe2["n_closed"] == 0
    assert pe2["blocked_m3"]["3mm"] < 0.05 * v_in
    assert abs(pe2["blocked_m3"]["12mm"] / v_in - 1) < 0.1, (pe2, v_in)
    return dict(slab_t_mm=mw["t_median_m"] * 1e3, overhang=ov["frac"], overhang_exact=exact, closed=pe["n_closed"],
                blocked_3mm_rel=pe2["blocked_m3"]["3mm"] / v_in, blocked_12mm_rel=pe2["blocked_m3"]["12mm"] / v_in)


def test_lattices_murray_blade():
    tl, rho = S.tpms_density_table("gyroid")
    assert np.all(np.diff(rho) >= -1e-12) and rho[0] == 0
    L, tt = 0.01, 0.15
    phi, _ = C.voxelize(lambda P: S.tpms_sheet_sdf(P, "gyroid", L, tt * L), (0, 0, 0), (2 * L,) * 3, L / 40)
    vf = float((phi < 0).mean())
    tab = float(np.interp(tt, tl, rho))
    assert abs(vf - tab) < 0.02, (vf, tab)
    t_back = float(S.tpms_t_for_density("gyroid", tab))
    assert abs(t_back - tt) < 0.011
    rng = np.random.Generator(np.random.PCG64(2))
    ang = np.linspace(0, 2 * np.pi, 14, endpoint=False)
    T = np.stack([np.cos(ang), np.sin(ang)], 1)
    segs = S.murray_tree((0.0, -2.0), T, 0.01, k=2.41, order=np.argsort(ang))
    leaves = [s for s in segs if s["n_leaves"] == 1]
    assert len(leaves) == 14
    ends = np.array([s["pts"][-1] for s in leaves])
    assert np.allclose(np.sort(ends[:, 0]), np.sort(T[:, 0]))
    root = max(segs, key=lambda s: s["n_leaves"])
    assert abs(root["d"] - 0.01 * 14 ** (1 / 2.41)) < 1e-9
    G, F, laws = BL.blade(0.1, 0.2, (0.05, 0.04), (10.0, 35.0), (0.12, 0.08), bow=0.01, tubercle=(0.05, 0.4), serration=(0.03, 0.3))
    chk = BL.check_blade(G, laws)
    assert chk["self_intersections"] == 0
    dtw = np.array(chk["twist_measured_deg"]) - laws["twist_deg"]
    assert np.max(np.abs(dtw - np.median(dtw))) < 1.5, np.max(np.abs(dtw - np.median(dtw)))
    sym = FT.angular_symmetry(np.cos(5 * np.linspace(0, 2 * np.pi, 360, endpoint=False)))
    assert sym["N_eff"] == 5
    return dict(gyroid_vf_voxels=vf, gyroid_vf_table=tab, root_d=root["d"], twist_dev_max=float(np.max(np.abs(dtw - np.median(dtw)))),
                sym_order=sym["N_eff"])


def test_mesh_volumes():
    h = 0.01
    phi, ax = C.voxelize(lambda P: S.sd_sphere(P, (0, 0, 0), 0.3), (-0.4,) * 3, (0.4,) * 3, h)
    V, F = M.surface_nets(phi, (ax[0][0], ax[1][0], ax[2][0]), h)
    ve = 4 / 3 * math.pi * 0.3 ** 3
    assert abs(M.volume(V, F) / ve - 1) < 0.01
    z = np.linspace(0, 1, 60); th = np.linspace(0, 2 * np.pi, 180, endpoint=False)
    V2, F2 = M.revolve_heightfield(z, th, np.full((60, 180), 0.5), np.full(60, 0.45))
    ve2 = math.pi * (0.25 - 0.2025)
    assert abs(M.volume(V2, F2) / ve2 - 1) < 0.01
    return dict(sphere=M.volume(V, F) / ve, ring=M.volume(V2, F2) / ve2)


def test_render():
    z = np.linspace(0, 1.0, 40); th = np.linspace(0, 2 * np.pi, 120, endpoint=False)
    ro = 0.6 + 0.05 * np.cos(6 * th)[None, :] * np.sin(np.pi * z)[:, None]
    V, F = M.revolve_heightfield(z, th, ro, np.full(40, 0.5))
    d = tempfile.mkdtemp(prefix="evo_render_test_")
    out = os.path.join(d, "sheet.png")
    res = R.render_sheet([(V, F, 0.8)], out, "Проверка рендера", ["строка 1"], backend="matplotlib", view_px=(360, 280))
    import matplotlib.pyplot as plt
    img = plt.imread(out)
    dark = float((img[..., :3].mean(axis=2) < 0.9).mean())
    assert img.shape[1] >= 1000 and 0.02 < dark < 0.8, (img.shape, dark)
    out2 = os.path.join(d, "sheet2.png")
    R.render_sheet([(V, F, 0.8)], out2, "Проверка рендера", ["строка 1"], backend="matplotlib", view_px=(360, 280))
    assert np.array_equal(plt.imread(out2), img)
    th_png = os.path.join(d, "axo.png")
    R.render_axo([(V, F, 0.8)], th_png, size=(200, 160), backend="matplotlib", label="A1")
    cs = R.contact_sheet([th_png, th_png], ["A1", "A2"], os.path.join(d, "cs.png"), ncols=2, title="галерея")
    pv = None
    if R.have_pyvista():
        pv = R.render_sheet([(V, F, 0.8)], os.path.join(d, "pv.png"), "pyvista", backend="pyvista", view_px=(360, 280))["backend"]
    return dict(dark_frac=dark, size=list(img.shape), contact=os.path.exists(cs), pyvista=pv, dir=d)


if __name__ == "__main__":
    out, ok = {}, True
    for name, fn in list(globals().items()):
        if name.startswith("test_") and callable(fn):
            try:
                out[name] = fn()
                print("PASS", name, json.dumps(out[name], ensure_ascii=False, default=float)[:300])
            except Exception as ex:  # noqa: BLE001
                ok = False
                import traceback
                traceback.print_exc()
                print("FAIL", name, repr(ex)[:500])
    sys.exit(0 if ok else 1)
