"""Геометрия решёток TPMS: поля уровня, пористость, удельная поверхность, гидравлический диаметр,
изоповерхности (marching cubes) и запись STL. Единица длины — размер ячейки a = 1.

Поверхности (стандартные приближения уровня, Schoen 1970; Schnering, Nesper 1991 — формы общеизвестны):
  гироид    G = sin x cos y + sin y cos z + sin z cos x
  Schwarz-D D = sin x sin y sin z + sin x cos y cos z + cos x sin y cos z + cos x cos y sin z
  x = 2π·X/a. Лист (стенка) двухтрактного теплообменника: c − t < F < c + t; тракт A: F > c + t, тракт B: F < c − t.
Обе поверхности делят ячейку на два связных лабиринта; у гироида F(−r) = −F(r), у D F(x + a/2) = −F(x):
тракт B с уровнем ℓ — зеркальная копия тракта A с уровнем −ℓ (переносы одинаковы).
"""
import numpy as np

TOPOS = {"gyroid": "гироид", "diamond": "Schwarz-D (алмаз)"}


def field(topo, X, Y, Z):
    if topo == "gyroid":
        return np.sin(X) * np.cos(Y) + np.sin(Y) * np.cos(Z) + np.sin(Z) * np.cos(X)
    if topo == "diamond":
        sx, sy, sz, cx, cy, cz = np.sin(X), np.sin(Y), np.sin(Z), np.cos(X), np.cos(Y), np.cos(Z)
        return sx * sy * sz + sx * cy * cz + cx * sy * cz + cx * cy * sz
    raise ValueError(topo)


def grid(topo, n, centered=True, periods=(1, 1, 1)):
    """поле на n точек на ячейку; centered — центры вокселей (для решателя), иначе узлы с замыканием (для MC)"""
    axes = []
    for p in periods:
        if centered:
            x = (np.arange(n * p) + 0.5) / n
        else:
            x = np.arange(n * p + 1) / n
        axes.append(2 * np.pi * x)
    X, Y, Z = np.meshgrid(*axes, indexing="ij")
    return field(topo, X, Y, Z)


def voxel_mask(topo, level, n, shift=(0.0, 0.0, 0.0), eps_target=None):
    """воксели тракта F > level на сетке n³, сдвинутой на shift (доли вокселя); при eps_target порог
    сдвигается так, чтобы доля жидких вокселей равнялась истинной пористости (снимает ошибку объёма
    при наложении сетки на период — главную ошибку первого порядка)"""
    axes = [2 * np.pi * (np.arange(n) + 0.5 + s) / n for s in shift]
    X, Y, Z = np.meshgrid(*axes, indexing="ij")
    F = field(topo, X, Y, Z)
    thr = level
    if eps_target is not None:
        thr = float(np.quantile(F, 1.0 - eps_target))
    return F > thr, thr


def fraction_above(topo, level, n=160):
    return float((grid(topo, n) > level).mean())


def iso_mesh(F, level, spacing):
    from skimage import measure
    v, f, _, _ = measure.marching_cubes(F, level=level, spacing=spacing)
    return v, f


def mesh_area(v, f):
    a = v[f[:, 1]] - v[f[:, 0]]
    b = v[f[:, 2]] - v[f[:, 0]]
    return float(0.5 * np.linalg.norm(np.cross(a, b), axis=1).sum())


def surface_density(topo, level, n=160):
    """площадь поверхности F = level на единицу объёма ячейки, 1/a (узловая сетка с замыканием периода)"""
    F = grid(topo, n, centered=False)
    v, f = iso_mesh(F, level, (1.0 / n,) * 3)
    return mesh_area(v, f)


def side_props(topo, level, n=160):
    """тракт F > level: пористость ε (доля объёма ячейки), β (1/a), Dh = 4ε/β (в долях a)"""
    eps = fraction_above(topo, level, n)
    beta = surface_density(topo, level, n)
    return dict(eps=eps, beta=beta, Dh=4 * eps / beta)


def sheet_props(topo, c, t, n=160):
    """лист c ± t: доли трактов и стенки, удельная поверхность каждой стороны, средняя толщина стенки (в долях a)"""
    A = side_props(topo, c + t, n)
    Bm = side_props(topo, -(c - t), n)          # тракт B = зеркало тракта A с уровнем −(c − t)
    mid = surface_density(topo, c, n)
    e_wall = 1.0 - A["eps"] - Bm["eps"]
    return dict(A=A, B=Bm, eps_wall=e_wall, beta_mid=mid, t_wall=e_wall / mid)


def write_stl(path, v, f, solid="tpms"):
    """двоичный STL"""
    v = np.asarray(v, np.float32)
    tri = v[f]
    nrm = np.cross(tri[:, 1] - tri[:, 0], tri[:, 2] - tri[:, 0])
    ln = np.linalg.norm(nrm, axis=1, keepdims=True)
    ln[ln == 0] = 1
    nrm = (nrm / ln).astype(np.float32)
    rec = np.zeros(len(f), dtype=[("n", "<f4", 3), ("v", "<f4", (3, 3)), ("a", "<u2")])
    rec["n"] = nrm
    rec["v"] = tri
    with open(path, "wb") as fh:
        hdr = (solid + " binary STL").encode()[:80].ljust(80, b" ")
        fh.write(hdr)
        fh.write(np.uint32(len(f)).tobytes())
        fh.write(rec.tobytes())
    return dict(path=str(path), triangles=int(len(f)), bytes=84 + 50 * len(f))


def decimate_stl(src, dst, max_bytes=4.8e6):
    """уменьшенная копия STL для репозитория (≤ max_bytes): квадратичное упрощение (trimesh + fast-simplification);
    полный файл остаётся в scratchpad"""
    import trimesh
    m = trimesh.load(src)
    target = int(max_bytes / 50)
    if len(m.faces) > target:
        m = m.simplify_quadric_decimation(face_count=target)
    m.export(dst)
    return dict(path=str(dst), triangles=int(len(m.faces)), bytes=int(84 + 50 * len(m.faces)),
                source_triangles=int(len(trimesh.load(src).faces)))
