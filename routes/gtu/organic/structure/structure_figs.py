#!/usr/bin/env python3
"""Рисунки и рендеры узла organic/structure (по results-JSON и полям в <scratchpad>/organic/).

structure_disks.png   — сечения: исходный профиль и проекты A/B диска турбины и колеса, масса и момент инерции
structure_lattice.png — гомогенизация гироида: E*/E и прочность S от ρ*
structure_shell.png   — рама оболочки: проект против 8 прямых стоек, жёсткость опоры по направлениям
render_*.png          — тела вращения и решётка (pyvista без экрана)
Палитра — категориальные слоты 1–4 (синий, оранжевый, бирюзовый, жёлтый), проверены validate_palette.js (light: PASS).
"""
import json
import math
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import organic_geom as OG  # noqa: E402

C1, C2, C3, C4 = "#2a78d6", "#eb6834", "#1baf7a", "#eda100"
INK, INK2, MUTED, SURF = "#0b0b0b", "#52514e", "#898781", "#fcfcfb"


def _style():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.size": 9, "axes.edgecolor": MUTED, "axes.labelcolor": INK2, "xtick.color": MUTED,
                         "ytick.color": MUTED, "axes.titlecolor": INK, "figure.facecolor": "white", "axes.facecolor": SURF,
                         "axes.grid": False})
    return plt


def field(name):
    D = OG.load(name)
    return D


def disks_figure(path, sel):
    """sel: [(part, [(name|None для исходного, подпись, вид 'x'|'rho')...])]"""
    plt = _style()
    from matplotlib.colors import LinearSegmentedColormap
    seq = LinearSegmentedColormap.from_list("seq", ["#ffffff", "#cde2fb", "#86b6ef", "#3987e5", "#1c5cab", "#0d366b"])
    gray = LinearSegmentedColormap.from_list("g", ["#ffffff", "#3a3936"])
    ncol = max(len(v) for _, v in sel)
    fig, axs = plt.subplots(len(sel), ncol, figsize=(3.3 * ncol, 3.6 * len(sel)), squeeze=False)
    for row, (part, items) in enumerate(sel):
        for col, (name, label, kind) in enumerate(items):
            ax = axs[row, col]
            if name is None:
                D = field(items[1][0])
                base = D["base"]
                ax.fill(base[:, 0] * 1e3, base[:, 1] * 1e3, color="#3a3936", lw=0)
            else:
                D = field(name)
                F = D["F"].copy()
                ext = [0, D["r1"] * 1e3, D["z0"] * 1e3, D["z1"] * 1e3]
                if kind == "rho":
                    G = D["G"]
                    import lattice as LT
                    rmin = LT.default_law()["rho_min"]
                    R = np.where(G == 2, rmin + (1 - rmin) * np.nan_to_num(F), np.where(G == 0, 1.0, np.nan))
                    im = ax.imshow(R.T, origin="lower", extent=ext, cmap=seq, vmin=0, vmax=1, interpolation="nearest")
                else:
                    im = ax.imshow(F.T, origin="lower", extent=ext, cmap=gray, vmin=0, vmax=1, interpolation="nearest")
                base = D["base"]
                ax.plot(base[:, 0] * 1e3, base[:, 1] * 1e3, color=C2, lw=1.0)
            ax.set_aspect("equal")
            ax.set_title(label, fontsize=8.5, loc="left")
            ax.set_xlabel("r, мм")
            if col == 0:
                ax.set_ylabel("z, мм")
            for s in ("top", "right"):
                ax.spines[s].set_visible(False)
        for col in range(len(items), ncol):
            axs[row, col].axis("off")
    fig.tight_layout()
    fig.savefig(path, dpi=130)
    plt.close(fig)


def lattice_figure(path):
    plt = _style()
    with open(os.path.join(HERE, "lattice_results.json"), encoding="utf-8") as f:
        L = json.load(f)
    P = L["points"]
    r = np.array([p["rho"] for p in P])
    E = np.array([p["Emin_dir"] for p in P]); S = np.array([p["S_worst_p995"] for p in P])
    fig, ax = plt.subplots(1, 2, figsize=(9.5, 3.8))
    fe, fs = L["fit_power"]["E"], L["fit_power"]["S"]
    xx = np.linspace(0.08, 1.0, 50)
    for a, y, fit, lab, c in ((ax[0], E, fe, "E*/E (наименьший по [100], [111])", C1), (ax[1], S, fs, "S = σ*/σ (начало текучести, 99,5 %)", C2)):
        a.loglog(r, y, "o", color=c, ms=6, label="гомогенизация, 3D-КЭ")
        a.loglog(xx, np.exp(fit["lnC"]) * xx ** fit["n"], "-", color=c, lw=1.5,
                 label=f"степенная подгонка 0,10–0,45: n = {fit['n']:.2f}")
        a.loglog(xx, xx ** 2, ":", color=MUTED, lw=1.2, label="Гибсон — Эшби, изгиб: ρ*²")
        a.loglog(xx, xx, "--", color=MUTED, lw=1.0, label="правило смеси: ρ*")
        a.set_xlabel("относительная плотность ρ*"); a.set_title(lab, fontsize=9, loc="left")
        a.legend(fontsize=7, frameon=False)
        for s in ("top", "right"):
            a.spines[s].set_visible(False)
    fig.tight_layout(); fig.savefig(path, dpi=130); plt.close(fig)


def shell_figure(path):
    plt = _style()
    from matplotlib.colors import LinearSegmentedColormap
    gray = LinearSegmentedColormap.from_list("g", ["#ffffff", "#3a3936"])
    d = np.load(os.path.join(OG.OUTDIR, "shell_frame.npz"))
    with open(os.path.join(HERE, "shell_results.json"), encoding="utf-8") as f:
        R = json.load(f)
    n = int(d["n"]); Rr = float(d["R"]) * 1e3
    fig = plt.figure(figsize=(12, 4.2))
    for k, (key, lab) in enumerate((("xe", f"оптимизация топологии: {R['optimized']['mass_kg']:.0f} кг"),
                                    ("xr", f"8 прямых стоек: {R['reference_8_spokes']['mass_kg']:.0f} кг"))):
        ax = fig.add_subplot(1, 3, k + 1)
        F = np.full((n, n), np.nan); F[d["ei"], d["ej"]] = d[key]
        ax.imshow(F.T, origin="lower", extent=[-Rr, Rr, -Rr, Rr], cmap=gray, vmin=0, vmax=1, interpolation="nearest")
        ax.set_aspect("equal"); ax.set_title(lab, fontsize=9, loc="left"); ax.set_xlabel("x, мм")
        if k == 0:
            ax.set_ylabel("y, мм")
        for s in ("top", "right"):
            ax.spines[s].set_visible(False)
    ax = fig.add_subplot(1, 3, 3, projection="polar")
    th = np.radians(np.arange(0, 360, 45)); th = np.r_[th, th[0]]
    ko = np.array(R["optimized"]["k_bearing_N_per_m"]) / 1e6; kr = np.array(R["reference_8_spokes"]["k_bearing_N_per_m"]) / 1e6
    ax.plot(th, np.r_[ko, ko[0]], "-o", color=C1, lw=2, ms=5, label="оптимизация топологии")
    ax.plot(th, np.r_[kr, kr[0]], "-s", color=C2, lw=2, ms=5, label="8 прямых стоек")
    ax.set_title("жёсткость опоры подшипника, Н/мкм", fontsize=9, loc="left")
    ax.legend(fontsize=7, frameon=False, loc="lower left", bbox_to_anchor=(-0.15, -0.18))
    fig.tight_layout(); fig.savefig(path, dpi=130); plt.close(fig)


def renders(out_dir, scratch, items):
    """items: [(stl_path, png_path, title, clip)] — тела вращения с вырезом, решётка"""
    import trimesh
    done = []
    for stl, png, title, clip in items:
        if not os.path.exists(stl):
            continue
        m = trimesh.load(stl)
        OG.render_mesh(m, png, title=None, clip=clip)
        done.append(png)
    return done
