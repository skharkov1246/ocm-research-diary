"""Рендер элит: лист из трёх ракурсов в одном масштабе, монохром на белом, как эскизы воронки (funnel/sketch):
1 — вид сбоку (ось z горизонтально, вход слева, смотрим вдоль +x), 2 — вид с торца со стороны входа (смотрим вдоль +z),
3 — аксонометрия 3/4 (азимут 55°, возвышение 30°), параллельная проекция. Шрифт — Liberation Serif (метрический аналог
Times New Roman), без жирного и цветных плашек; подпись «эскиз формы (оценка), не является конструкторской документацией».

Два способа: pyvista (OSMesa, вне экрана; быстрее и с силуэтом) и matplotlib (запасной: есть везде, в том числе на AWS).
Сетки — список (V, F, серый 0…1). Ось модели: z — ось двигателя, y — вверх, x — вбок.
"""
import math
import os

import numpy as np

VIEWS = (("side", "1 — вид сбоку"), ("end", "2 — вид с торца со стороны входа"), ("axo", "3 — аксонометрия 3/4 (азимут 55°, возвышение 30°)"))
FOOT = "Эскиз формы по расчётной модели (оценка). Не является конструкторской документацией."


def _font():
    import matplotlib
    from matplotlib import font_manager
    for f in font_manager.findSystemFonts():
        if "LiberationSerif-Regular" in f:
            font_manager.fontManager.addfont(f)
            matplotlib.rcParams["font.family"] = "Liberation Serif"
            return "Liberation Serif"
    matplotlib.rcParams["font.family"] = "serif"
    return "serif"


def _cam(view):
    """направление от центра к камере и «вверх»"""
    if view == "side":
        return np.array([-1.0, 0.0, 0.0]), np.array([0.0, 1.0, 0.0])
    if view == "end":
        return np.array([0.0, 0.0, -1.0]), np.array([0.0, 1.0, 0.0])
    az, el = math.radians(55.0), math.radians(30.0)
    d = np.array([-math.cos(el) * math.sin(az), math.sin(el), -math.cos(el) * math.cos(az)])
    return d / np.linalg.norm(d), np.array([0.0, 1.0, 0.0])


def have_pyvista():
    try:
        import pyvista  # noqa: F401
        return True
    except Exception:  # noqa: BLE001
        return False


def _view_pyvista(meshes, view, size, center, half_h):
    import pyvista as pv
    pl = pv.Plotter(off_screen=True, window_size=list(size))
    pl.set_background("white")
    for V, F, gray in meshes:
        faces = np.hstack([np.full((len(F), 1), 3, np.int64), np.asarray(F, np.int64)]).ravel()
        m = pv.PolyData(np.asarray(V, float), faces)
        pl.add_mesh(m, color=(gray, gray, gray), smooth_shading=True, ambient=0.28, diffuse=0.72, specular=0.04)
        try:
            pl.add_silhouette(m, color="black", line_width=1.6)
        except Exception:  # noqa: BLE001
            pass
    d, up = _cam(view)
    pl.enable_parallel_projection()
    pl.camera.focal_point = tuple(center)
    pl.camera.position = tuple(center + d * 10.0)
    pl.camera.up = tuple(up)
    pl.camera.parallel_scale = half_h
    pl.camera.clipping_range = (0.01, 100.0)
    img = pl.screenshot(return_img=True)
    pl.close()
    return img


def _view_mpl(meshes, view, size, center, half_h):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d.art3d import Poly3DCollection
    d, up = _cam(view)
    right = np.cross(up, d); right /= np.linalg.norm(right)
    upv = np.cross(d, right)
    fig = plt.figure(figsize=(size[0] / 100, size[1] / 100), dpi=100)
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(-half_h * size[0] / size[1], half_h * size[0] / size[1]); ax.set_ylim(-half_h, half_h)
    ax.set_aspect("equal"); ax.axis("off")
    light = d * 0.6 + upv * 0.5 + right * 0.3; light /= np.linalg.norm(light)
    polys, cols, depth = [], [], []
    for V, F, gray in meshes:
        V = np.asarray(V, float) - center
        P2 = np.stack([V @ right, V @ upv], 1)
        dep = V @ d
        tri = np.asarray(F)
        n = np.cross(V[tri[:, 1]] - V[tri[:, 0]], V[tri[:, 2]] - V[tri[:, 0]])
        n /= np.linalg.norm(n, axis=1, keepdims=True) + 1e-30
        vis = (n @ d) > 0
        sh = np.clip(0.3 + 0.7 * np.abs(n @ light), 0, 1) * gray
        polys.append(P2[tri[vis]]); cols.append(sh[vis]); depth.append(dep[tri[vis]].mean(axis=1))
    P = np.concatenate(polys); C = np.concatenate(cols); D = np.concatenate(depth)
    o = np.argsort(D)
    from matplotlib.collections import PolyCollection
    pc = PolyCollection(P[o], facecolors=np.repeat(C[o][:, None], 3, 1), edgecolors="none", antialiased=False)
    ax.add_collection(pc)
    del Poly3DCollection
    fig.canvas.draw()
    img = np.asarray(fig.canvas.buffer_rgba())[..., :3].copy()
    plt.close(fig)
    return img


def render_sheet(meshes, out_png, title, notes=(), backend="auto", view_px=(900, 700), px_per_m=None, scale_bar_m=0.5,
                 center=None):
    """лист из трёх видов. px_per_m — общий масштаб набора (иначе — по габариту сетки). Возврат: dict(backend, px_per_m, path)"""
    _font()
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    allV = np.vstack([np.asarray(V, float) for V, _, _ in meshes])
    lo, hi = allV.min(axis=0), allV.max(axis=0)
    c = (lo + hi) / 2 if center is None else np.asarray(center, float)
    ext = float(np.linalg.norm(hi - lo))
    if px_per_m is None:
        px_per_m = 0.92 * min(view_px) / max(ext, 1e-9)
    half_h = view_px[1] / 2 / px_per_m
    be = backend
    if be == "auto":
        be = "pyvista" if have_pyvista() else "matplotlib"
    fn = _view_pyvista if be == "pyvista" else _view_mpl
    imgs = [fn(meshes, v, view_px, c, half_h) for v, _ in VIEWS]
    W, H = view_px
    fig = plt.figure(figsize=(3 * W / 100, (H + 230) / 100), dpi=100)
    fig.text(0.01, 1 - 30 / (H + 230), title, fontsize=15, va="top")
    for i, ((v, cap), img) in enumerate(zip(VIEWS, imgs)):
        ax = fig.add_axes([i / 3, 120 / (H + 230), 1 / 3, H / (H + 230)])
        ax.imshow(img); ax.axis("off")
        ax.text(10, 20, cap, fontsize=11, va="top")
        if v != "axo":
            L = scale_bar_m * px_per_m
            x0, y0 = 30, H - 30
            ax.plot([x0, x0 + L], [y0, y0], color="black", lw=2)
            ax.plot([x0, x0], [y0 - 6, y0 + 6], color="black", lw=1); ax.plot([x0 + L, x0 + L], [y0 - 6, y0 + 6], color="black", lw=1)
            ax.text(x0 + L / 2, y0 - 10, f"{scale_bar_m:g} м", ha="center", fontsize=10)
        else:
            ax.text(30, H - 20, "масштаб набора действует в плоскости проекции", fontsize=9)
    lines = list(notes) + [FOOT]
    for k, ln in enumerate(lines[:6]):
        fig.text(0.01, (100 - 17 * k) / (H + 230), ln, fontsize=10, va="top")
    os.makedirs(os.path.dirname(os.path.abspath(out_png)), exist_ok=True)
    fig.savefig(out_png, dpi=100, facecolor="white")
    plt.close(fig)
    return dict(backend=be, px_per_m=float(px_per_m), path=out_png, views=[v for v, _ in VIEWS])


def render_axo(meshes, out_png, size=(520, 420), px_per_m=None, backend="auto", label=None):
    """одна аксонометрия (миниатюра для галереи)"""
    _font()
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    allV = np.vstack([np.asarray(V, float) for V, _, _ in meshes])
    lo, hi = allV.min(axis=0), allV.max(axis=0)
    c = (lo + hi) / 2
    if px_per_m is None:
        px_per_m = 0.9 * min(size) / float(np.linalg.norm(hi - lo))
    be = "pyvista" if (backend == "auto" and have_pyvista()) else ("matplotlib" if backend == "auto" else backend)
    fn = _view_pyvista if be == "pyvista" else _view_mpl
    img = fn(meshes, "axo", size, c, size[1] / 2 / px_per_m)
    fig = plt.figure(figsize=(size[0] / 100, size[1] / 100), dpi=100)
    ax = fig.add_axes([0, 0, 1, 1]); ax.imshow(img); ax.axis("off")
    if label:
        ax.text(8, 16, label, fontsize=9, va="top")
    fig.savefig(out_png, dpi=100, facecolor="white")
    plt.close(fig)
    return dict(backend=be, px_per_m=float(px_per_m), path=out_png)


def contact_sheet(pngs, labels, out_png, ncols=6, title=None, foot=FOOT):
    """галерея миниатюр с подписями (номер элиты, ячейка, ключевые числа)"""
    _font()
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    imgs = [plt.imread(p) for p in pngs]
    if not imgs:
        raise ValueError("нет картинок")
    h, w = imgs[0].shape[:2]
    nrows = math.ceil(len(imgs) / ncols)
    fig = plt.figure(figsize=(ncols * w / 100, nrows * (h + 46) / 100 + 0.8), dpi=100)
    top = 0.5 / (nrows * (h + 46) / 100 + 0.8)
    if title:
        fig.text(0.01, 1 - 0.1 * top, title, fontsize=14, va="top")
    for k, (img, lab) in enumerate(zip(imgs, labels)):
        r, cc = divmod(k, ncols)
        ax = fig.add_axes([cc / ncols, 1 - top - (r + 1) / nrows * (1 - 2 * top), 1 / ncols, (1 - 2 * top) / nrows])
        ax.imshow(img); ax.axis("off")
        ax.set_title(lab, fontsize=8, pad=2)
    fig.text(0.01, 0.3 * top, foot, fontsize=10, va="bottom")
    fig.savefig(out_png, dpi=100, facecolor="white")
    plt.close(fig)
    return out_png
