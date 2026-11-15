#!/usr/bin/env python3
"""Геометрия проектов: контур x̄ = 0,5 -> многоугольники сечения (огибающая точная, «живые» границы — по полю),
STL телом вращения, гироид переменной плотности в секторе (марш кубов), рендеры pyvista (без экрана).

Большие файлы (полные STL) — в <scratchpad>/organic/, в репозиторий — уменьшенные (не больше ~5 МБ).
"""
import json
import math
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
OUTDIR = os.environ.get("ORGANIC_OUT", "/tmp/claude-0/-home-user/a62745e4-d2fd-56e0-974f-12746b1a2ae0/scratchpad/organic")


def load(name):
    d = np.load(os.path.join(OUTDIR, f"disk_{name}.npz"))
    nr, nz = int(d["nr"]), int(d["nz"])
    F = np.full((nr, nz), np.nan)
    F[d["ei"], d["ej"]] = d["xb"]
    G = np.full((nr, nz), -1)
    G[d["ei"], d["ej"]] = d["group"]
    return dict(F=F, G=G, nr=nr, nz=nz, r1=float(d["r1"]), z0=float(d["z0"]), z1=float(d["z1"]), env=d["env"],
                base=d["base"], hist=json.loads(str(d["hist"])))


def section_polygons(D, level=0.5, smooth_tol=None):
    """многоугольники сплошного металла: контур поля (вне огибающей — 1, у оси — зеркало) ∩ огибающая"""
    from skimage import measure
    import shapely
    from shapely.geometry import Polygon
    from shapely.ops import unary_union
    F = np.nan_to_num(D["F"].copy(), nan=1.0)
    nr, nz = F.shape
    dr, dz = D["r1"] / nr, (D["z1"] - D["z0"]) / nz
    P = np.ones((nr + 4, nz + 4))
    P[2:-2, 2:-2] = F
    P[0:2, 2:-2] = F[1::-1]          # зеркало относительно оси
    cs = measure.find_contours(P, level)
    env = Polygon(D["env"]).buffer(0)
    polys = []
    for c in cs:
        r = (c[:, 0] - 2 + 0.5) * dr; z = D["z0"] + (c[:, 1] - 2 + 0.5) * dz
        if len(c) < 4:
            continue
        pg = Polygon(np.c_[r, z]).buffer(0)
        if pg.area > 0:
            polys.append(pg)
    # контуры от find_contours: внешние и дыры — разделяем чётностью вложенности
    polys.sort(key=lambda p: -p.area)
    solid = None
    for p in polys:
        depth = sum(1 for q in polys if q is not p and q.area > p.area and q.contains(p.representative_point()))
        if solid is None:
            solid = p if depth % 2 == 0 else None
            continue
        solid = solid.union(p) if depth % 2 == 0 else solid.difference(p)
    if solid is None:
        solid = env
    res = solid.intersection(env)
    res = res.buffer(0)
    if smooth_tol:
        res = res.simplify(smooth_tol, preserve_topology=True)
    geoms = list(res.geoms) if res.geom_type == "MultiPolygon" else [res]
    geoms = [g for g in geoms if g.area > 4 * dr * dz]
    return geoms


def poly_rings(g, min_r=0.0):
    out = [np.asarray(g.exterior.coords)[:-1]]
    for i in g.interiors:
        out.append(np.asarray(i.coords)[:-1])
    for k in range(len(out)):
        out[k][:, 0] = np.maximum(out[k][:, 0], min_r)
    return out


def section_stats(geoms):
    """площадь, объём тела вращения на 1 кг/м³ (Паппа — Гульдина), замкнутые полости, нависания для печати"""
    A = sum(g.area for g in geoms)
    V = sum(2 * math.pi * g.centroid.x * g.area for g in geoms)
    holes = sum(len(g.interiors) for g in geoms)
    # нависание: рёбра границы, у которых наружная нормаль смотрит вниз круче 45° от горизонтали (печать осью вверх, +z)
    over = {"+z": 0.0, "-z": 0.0}
    per = 0.0
    for g in geoms:
        for ring, is_hole in [(np.asarray(g.exterior.coords), False)] + [(np.asarray(i.coords), True) for i in g.interiors]:
            # shapely: внешний контур против часовой (после buffer), дыры — по часовой; наружная нормаль = (dy, −dx)
            ccw = _signed_area(ring) > 0
            for a, b in zip(ring[:-1], ring[1:]):
                d = b - a
                L = math.hypot(*d)
                if L < 1e-12 or (a[0] < 1e-9 and b[0] < 1e-9):
                    continue
                nx, nz = (d[1], -d[0]) if ccw else (-d[1], d[0])
                nz /= L
                per += L
                if nz < -math.sqrt(0.5):
                    over["+z"] += L
                if nz > math.sqrt(0.5):
                    over["-z"] += L
    return dict(area_m2=A, volume_m3=V, closed_cavities=holes, perimeter_m=per,
                overhang_lt45_m={k: round(v, 4) for k, v in over.items()},
                overhang_share={k: round(v / max(per, 1e-12), 4) for k, v in over.items()})


def _signed_area(ring):
    x, y = ring[:, 0], ring[:, 1]
    return 0.5 * np.sum(x[:-1] * y[1:] - x[1:] * y[:-1])


def revolve_stl(geoms, path, sections=180, angle=2 * math.pi):
    """тело вращения: каждое кольцо многоугольника вращается вокруг оси z; дыры — с обратной ориентацией"""
    import trimesh
    parts = []
    for g in geoms:
        for k, ring in enumerate(poly_rings(g)):
            ring = np.vstack([ring, ring[:1]])
            m = trimesh.creation.revolve(ring, angle=angle, sections=sections)
            parts.append(m)
    mesh = trimesh.util.concatenate(parts)
    mesh.merge_vertices()
    mesh.export(path)
    return mesh


# =========================================================================== решётка
def lattice_sector(D, law, L_cell, sector_deg, step, path=None, t_skin=2e-3, max_tris=None):
    """гироид в секторе тела вращения: металл = обшивка (расстояние до огибающей ≤ t_skin) ∪ лист |g| ≤ t(ρ*(r,z))"""
    from skimage import measure
    import shapely
    from shapely.geometry import Polygon
    import lattice as LT
    env = Polygon(D["env"]).buffer(0)
    F = np.nan_to_num(D["F"], nan=0.0)
    nr, nz = F.shape
    dr, dz = D["r1"] / nr, (D["z1"] - D["z0"]) / nz
    rmin = law["rho_min"]
    R = D["r1"]
    th = math.radians(sector_deg)
    xs = np.arange(0, R + step, step)
    ys = np.arange(-R * math.sin(th / 2) - step, R * math.sin(th / 2) + step, step)
    zs = np.arange(D["z0"] - step, D["z1"] + step, step)
    X, Y, Z = np.meshgrid(xs, ys, zs, indexing="ij")
    r = np.hypot(X, Y)
    ang = np.abs(np.arctan2(Y, X))
    # поле ρ*(r, z) и расстояние до границы огибающей (2D, затем выборка)
    rr = np.arange(0, R + step, step / 2); zz = np.arange(D["z0"], D["z1"] + step, step / 2)
    RR, ZZ = np.meshgrid(rr, zz, indexing="ij")
    pts = shapely.points(RR.ravel(), ZZ.ravel())
    ins = shapely.contains(env, pts).reshape(RR.shape)
    bnd = shapely.LineString(np.asarray(D["env"])[1:-1])   # без оси
    dist = shapely.distance(bnd, pts).reshape(RR.shape)
    sd = np.where(ins, -dist, dist)                         # знаковое расстояние: внутри отрицательно
    ii = np.clip((r / dr).astype(int), 0, nr - 1); jj = np.clip(((Z - D["z0"]) / dz).astype(int), 0, nz - 1)
    rho = rmin + (1 - rmin) * F[ii, jj]
    a = np.clip(np.round(r / (step / 2)).astype(int), 0, len(rr) - 1); b = np.clip(np.round((Z - D["z0"]) / (step / 2)).astype(int), 0, len(zz) - 1)
    s = sd[a, b]
    # порог листа от ρ*: таблица
    tt = np.vectorize(LT.sheet_threshold)(np.linspace(rmin, 0.999, 60))
    tr = np.interp(rho, np.linspace(rmin, 0.999, 60), tt)
    G = np.abs(LT.gyroid(X, Y, Z, L_cell))
    g_lat = (G - tr) * (L_cell / (2 * math.pi)) / 1.2       # приблизительно в метрах (|∇g| ~ 1,2·2π/L)
    f_skin = -(s + t_skin)                                    # < 0 внутри обшивки… знак: металл где f < 0
    f_skin = s + t_skin                                       # металл обшивки: −t_skin ≤ s ≤ 0  ->  max(s, −s − t_skin) ≤ 0
    f_skin = np.maximum(s, -s - t_skin)
    f_in = np.maximum(s, g_lat)                               # решётка внутри
    f = np.minimum(f_skin, f_in)
    f = np.maximum(f, ang - th / 2)                           # сектор
    v, fc, _, _ = measure.marching_cubes(f.astype(np.float32), 0.0, spacing=(step, step, step))
    v += np.array([xs[0], ys[0], zs[0]])
    import trimesh
    mesh = trimesh.Trimesh(v, fc, process=True)
    if max_tris and len(mesh.faces) > max_tris:
        mesh = decimate(mesh, max_tris)
    if path:
        mesh.export(path)
    return mesh


def decimate(mesh, target):
    try:
        import fast_simplification  # noqa: F401
        return mesh.simplify_quadric_decimation(face_count=int(target))
    except Exception:  # noqa: BLE001
        import pyvista as pv
        pm = pv.wrap(mesh)
        red = max(0.0, min(0.99, 1 - target / max(len(mesh.faces), 1)))
        pm = pm.triangulate().decimate_pro(red, preserve_topology=True)
        import trimesh
        f = pm.faces.reshape(-1, 4)[:, 1:]
        return trimesh.Trimesh(pm.points, f, process=True)


# =========================================================================== рендеры
def render_mesh(mesh, path, title=None, color="#b8b2a7", cam=None, size=(1400, 1000), clip=None, edges=False):
    import pyvista as pv
    pv.OFF_SCREEN = True
    pm = pv.wrap(mesh)
    if clip is not None:
        pm = pm.clip(normal=clip[0], origin=clip[1], invert=False)
    p = pv.Plotter(off_screen=True, window_size=size)
    p.set_background("white")
    p.add_mesh(pm, color=color, smooth_shading=True, specular=0.4, specular_power=20, show_edges=edges)
    p.enable_eye_dome_lighting() if False else None
    if cam:
        p.camera_position = cam
    else:
        p.view_isometric()
    if title:
        p.add_text(title, font_size=12, color="black")
    p.screenshot(path)
    p.close()
    return path
