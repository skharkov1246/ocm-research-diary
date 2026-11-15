"""Сводный расчёт TPMS: геометрия → переносы (MAC + LBM) → сверка методов → рекуператор KV-3 → STL и картинки.

Запуск (после tpms_transport.py; файлы переносов — аргументы):
  python tpms.py --mac <tpms_transport_mac.json> [--lbm <tpms_transport_lbm.json>] --out <папка>
Выход: tpms_results.json, tpms_*.png, tpms_graded_demo.stl (уменьшенный образец; полный — в scratchpad).
"""
import argparse
import json
import math
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
sys.path.insert(0, str(HERE))
from tpms_geom import TOPOS, side_props, sheet_props, field, write_stl, surface_density  # noqa: E402
import tpms_recup as R  # noqa: E402

PLATE = dict(fRe=24.0, NuH1=8.2355, NuT=7.5406, src="Shah & London 1978 (справочные значения плоского канала); "
                                                     "те же константы — organic/cycles/cycles_results.json: laminar_constants")


def geometry_table():
    out = {}
    for topo in TOPOS:
        rows = []
        for l in np.round(np.arange(-0.8, 0.81, 0.2), 3):
            p = side_props(topo, float(l), 160)
            rows.append(dict(level=float(l), eps=p["eps"], beta_per_a=p["beta"], Dh_over_a=p["Dh"]))
        sh = sheet_props(topo, 0.0, 0.15, 160)
        out[topo] = dict(sides=rows, sheet_c0_t015=dict(eps_gas=sh["A"]["eps"], eps_air=sh["B"]["eps"],
                                                       eps_wall=sh["eps_wall"], t_wall_over_a=sh["t_wall"],
                                                       beta_mid_per_a=sh["beta_mid"]),
                         minimal_surface_area_per_a=surface_density(topo, 0.0, 240))
    return out


def transport_summary(mac_rows, lbm_rows):
    """свести по уровню: f·Re (среднее по сдвигам сетки и разброс), Nu(Pe)"""
    out = {}
    for topo in TOPOS:
        rows = [r for r in mac_rows if r["topo"] == topo]
        lv = sorted(set(r["level"] for r in rows))
        tab = []
        for l in lv:
            rr = [r for r in rows if r["level"] == l]
            f = [r["fRe"] for r in rr]
            nu0 = [r["Nu"][0]["Nu"] for r in rr]
            pe = max(rr, key=lambda r: len(r["Nu"]))["Nu"]
            tab.append(dict(level=l, eps=rr[0]["eps"], Dh_over_a=rr[0]["Dh"], beta_per_a=rr[0]["beta"],
                            n=sorted(set(r["n"] for r in rr)), fRe=float(np.mean(f)),
                            fRe_spread_pct=float(100 * (max(f) - min(f)) / np.mean(f)),
                            k_over_a2=float(np.mean([r["k_a2"] for r in rr])), Nu_H1_stokes=float(np.mean(nu0)),
                            Nu_vs_Pe=[dict(Pe=x["Pe"], Nu=x["Nu"]) for x in pe],
                            ratio_to_plate=dict(fRe=float(np.mean(f)) / PLATE["fRe"], Nu=float(np.mean(nu0)) / PLATE["NuH1"])))
        out[topo] = tab
    lbm = {}
    for topo in TOPOS:
        rr = sorted([r for r in lbm_rows if r["topo"] == topo], key=lambda r: (r["level"], r["Re"]))
        if rr:
            lbm[topo] = [dict(level=r["level"], n=r["n"], Re=r["Re"], fRe=r["fRe"], Nu=r["Nu"], Ma=r["Ma"],
                              converged=r["hist"][-1]["converged"] if r.get("hist") else None, steps=r["steps"])
                         for r in rr]
    return out, lbm


def inertia_from_lbm(lbm, mac_tab, topo="gyroid", level=0.0):
    """отношения f·Re/(f·Re)_Стокс и Nu/Nu_Стокс по LBM — для поправки в рекуператоре (оценка: одна ячейка, одна
    сетка; переносится на все уровни)"""
    if topo not in lbm:
        return None
    rr = [r for r in lbm[topo] if abs(r["level"] - level) < 1e-9]
    if not rr:
        return None
    base = [t for t in mac_tab[topo] if abs(t["level"] - level) < 1e-9][0]
    Re = [0.0] + [r["Re"] for r in rr]
    fr = [1.0] + [r["fRe"] / base["fRe"] for r in rr]
    nr = [1.0] + [r["Nu"] / base["Nu_H1_stokes"] for r in rr]
    return dict(Re=Re, f_ratio=fr, Nu_ratio=nr, topo=topo, level=level,
                note="поправка по LBM n = 32 для одного уровня; за пределом расчёта по Re — постоянна (не экстраполируется)")


def crosscheck(level=0.0, n=32):
    """предел Стокса двумя методами на одной воксельной геометрии: MAC (прямо/MINRES) и LBM (TRT)"""
    from tpms_geom import voxel_mask
    from voxflow import PeriodicVoxelFlow, percolating
    from lbm import LBMCell
    g = side_props("gyroid", level, 200)
    m, _ = voxel_mask("gyroid", level, n, (0, 0, 0), eps_target=g["eps"])
    m = percolating(m)
    s = PeriodicVoxelFlow(m)
    s.stokes(solver="minres")
    us_mac = float(s.U[0].sum() / m.size)
    Dh, A = g["Dh"] * n, g["beta"] * n ** 2
    h = s.scalar_h1(1e-3 / (us_mac / g["eps"] * Dh))
    Nu_mac = h["Q_balance"] * Dh / (A * (-h["phib"]))
    nu, F = 1 / 6, 1e-6
    L = LBMCell(m, nu, F, kappa=nu / 0.7)
    r = L.run(max_steps=30000, check=500, tol=1e-9, min_steps=3000)
    ob = r["obs"]
    us_lbm = ob["U"] * m.mean() * nu / F
    Nu_lbm = ob["Q"] * Dh / (A * (nu / 0.7) * (-ob["phib"]))
    return dict(topo="gyroid", level=level, n=n, k_mac_over_a2=us_mac / n ** 2, k_lbm_over_a2=us_lbm / n ** 2,
                k_diff_pct=100 * (us_lbm / us_mac - 1), Nu_mac=float(Nu_mac), Nu_lbm=float(Nu_lbm),
                Nu_diff_pct=100 * (Nu_lbm / Nu_mac - 1), lbm_steps=r["steps"], lbm_converged=r["converged"],
                Pe_lbm=float(ob["U"] * Dh / (nu / 0.7)))


def recuperators(mac_rows, inertia, quick=False):
    du_all, src = R.duties()
    out = dict(source=src, print_rules=R.PRINT, plate_constants=PLATE, cases={})
    models = {t: R.TPMSModel(t, mac_rows, inertia=None) for t in TOPOS}
    models_in = {t: R.TPMSModel(t, mac_rows, inertia=inertia) for t in TOPOS} if inertia else {}
    for m in list(models.values()) + list(models_in.values()):
        m.build_t_table()
    keys = list(du_all) if not quick else ["R90-ic-RH-ibc-none-none"]
    profiles = {}
    for key in keys:
        du = du_all[key]
        PH = R.Props("gas", du["p_h"], du["far"])
        PC = R.Props("air", du["p_c"])
        case = dict(duty={k: v for k, v in du.items()}, designs={})
        t0 = time.time()
        pl = R.size_plate(du, PH, PC, 1.0e-3, PLATE["NuT"])
        case["designs"]["plate_NuT_as_cycles"] = strip(pl)
        pl2 = R.size_plate(du, PH, PC, 1.0e-3, PLATE["NuH1"])
        case["designs"]["plate_NuH1"] = strip(pl2)
        pl3 = R.size_plate(du, PH, PC, 0.75e-3, PLATE["NuH1"])          # тот же Dh 1,5 мм, что у TPMS
        case["designs"]["plate_NuH1_Dh15"] = strip(pl3)
        for nm, gap, Nu in (("plate_NuT_as_cycles", 1.0e-3, PLATE["NuT"]), ("plate_NuH1", 1.0e-3, PLATE["NuH1"]),
                            ("plate_NuH1_Dh15", 0.75e-3, PLATE["NuH1"])):
            mf = R.min_frontal(du, PH, PC, "plate", dict(gap=gap, Nu=Nu))
            if mf is not None:
                case["designs"][nm]["min_frontal"] = strip(mf)
        runs = [("gyroid", True, None, models), ("diamond", True, None, models)]
        if key == "R90-ic-RH-ibc-none-none":
            runs += [("gyroid", False, None, models), ("diamond", False, None, models),
                     ("gyroid", True, R.PRINT["t_w_thin"], models)]
            if models_in:
                runs += [("gyroid", True, None, models_in)]
        for topo, graded, tw, mm in runs:
            name = f"{topo}_{'graded' if graded else 'uniform'}" + ("_tw025" if tw else "") + \
                   ("_inertia" if mm is models_in else "")
            r = R.size_tpms(du, PH, PC, mm[topo], graded=graded, tw=tw, maxiter=25 if quick else 50)
            z = r["z"]
            prof = dict(a=R.pw(z[2:5]) if len(z) == 8 else R.pw([z[2]] * 3),
                        c=R.pw(z[5:8]) if len(z) == 8 else R.pw([z[3]] * 3))
            mf = R.min_frontal(du, PH, PC, "tpms", prof, tw=tw, model=mm[topo])
            case["designs"][name] = strip(r)
            if mf is not None:
                case["designs"][name]["min_frontal"] = strip(mf)
            if key == "R90-ic-RH-ibc-none-none":
                profiles[name] = r.get("_prof")
            print(key, name, round(r["mass_core"]), round(r["V"], 3), r["feasible"], round(time.time() - t0), flush=True)
        profiles[f"plate_{key}"] = pl2.get("_prof")
        out["cases"][key] = case
    return out, profiles


def strip(r):
    o = {}
    for k, v in r.items():
        if k.startswith("_"):
            continue
        if isinstance(v, (np.floating,)):
            v = float(v)
        o[k] = v
    return o


_MAC = None


def plots(res, profiles, out_dir):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    # 1) переносы против плоского канала
    fig, axs = plt.subplots(1, 2, figsize=(10, 3.8))
    for topo, mk in (("gyroid", "o"), ("diamond", "s")):
        tab = res["transport"]["mac_summary"][topo]
        axs[0].plot([t["eps"] for t in tab], [t["fRe"] for t in tab], marker=mk, label=TOPOS[topo])
        axs[1].plot([t["eps"] for t in tab], [t["Nu_H1_stokes"] for t in tab], marker=mk, label=TOPOS[topo])
    axs[0].axhline(PLATE["fRe"], color="0.4", ls="--"); axs[0].text(0.3, PLATE["fRe"] + 1, "плоский канал 24", fontsize=8)
    axs[1].axhline(PLATE["NuH1"], color="0.4", ls="--"); axs[1].text(0.3, PLATE["NuH1"] + 0.1, "плоский канал 8,235", fontsize=8)
    axs[0].set_xlabel("пористость тракта ε"); axs[0].set_ylabel("f·Re (Фаннинг, по Dh = 4ε/β)")
    axs[1].set_xlabel("пористость тракта ε"); axs[1].set_ylabel("Nu_H1 (предел Стокса)")
    for ax in axs:
        ax.grid(alpha=0.3); ax.legend(fontsize=8)
    fig.suptitle("Тракт TPMS в пределе малых Re: расчёт MAC (воксели 32³, два сдвига сетки)", fontsize=10)
    fig.tight_layout(); fig.savefig(out_dir / "tpms_transport.png", dpi=130); plt.close(fig)
    # 2) рекуператор: масса и объём
    rc = res["recuperator"]["cases"]
    fig, axs = plt.subplots(1, 3, figsize=(13, 4))
    keys = list(rc)
    labels = {"plate_NuH1": "пластины, Dh 2 мм", "plate_NuH1_Dh15": "пластины, Dh 1,5 мм",
              "gyroid_uniform": "гироид, ячейка пост.", "gyroid_graded": "гироид, переменная",
              "diamond_uniform": "Schwarz-D, пост.", "diamond_graded": "Schwarz-D, переменная",
              "gyroid_graded_tw025": "гироид, стенка 0,25", "gyroid_graded_inertia": "гироид + инерция (LBM)"}
    for ax, field_, yl in ((axs[0], "mass_core", "масса сердцевины, т"), (axs[1], "V", "объём сердцевины, м³"),
                           (axs[2], "A_fr", "наименьший фронт при тех же потерях, м²")):
        names = [n for n in labels if n in rc.get("R90-ic-RH-ibc-none-none", {}).get("designs", {})]
        D = rc["R90-ic-RH-ibc-none-none"]["designs"]
        vals = [(D[n].get("min_frontal", {}).get("A_fr", D[n]["A_fr"]) if field_ == "A_fr" else D[n][field_])
                / (1e3 if field_ == "mass_core" else 1) for n in names]
        cols = ["0.55" if n.startswith("plate") else "#3b7dd8" for n in names]
        ax.barh(range(len(names)), vals, color=cols)
        ax.set_yticks(range(len(names))); ax.set_yticklabels([labels[n] for n in names], fontsize=8)
        ax.invert_yaxis(); ax.set_xlabel(yl); ax.grid(axis="x", alpha=0.3)
        for i, v in enumerate(vals):
            ax.text(v, i, f" {v:.2f}", va="center", fontsize=7.5)
    fig.suptitle("Рекуператор KV-3, ε = 0,90 + подогрев (нагрузка organic/cycles): стенка 0,35 мм, те же потери", fontsize=10)
    fig.tight_layout(); fig.savefig(out_dir / "tpms_recuperator.png", dpi=130); plt.close(fig)
    # 3) профиль переменной ячейки и производство энтропии
    pr = profiles.get("gyroid_graded")
    pu = profiles.get("gyroid_uniform")
    if pr is not None:
        fig, axs = plt.subplots(1, 3, figsize=(13, 3.6))
        s = pr["x"] / pr["x"][-1]
        axs[0].plot(s, pr["a"] * 1e3, label="ячейка a, мм")
        axs[0].plot(s, pr["dA"] * 1e3, label="Dh газа, мм"); axs[0].plot(s, pr["dB"] * 1e3, label="Dh воздуха, мм")
        axs[0].set_xlabel("x/L (газ входит при 0)"); axs[0].legend(fontsize=8); axs[0].grid(alpha=0.3)
        axs[1].plot(s, pr["c"], color="C3"); axs[1].set_xlabel("x/L"); axs[1].set_ylabel("смещение уровня c (< 0 — газу больше)")
        axs[1].grid(alpha=0.3)
        axs[2].plot(s, pr["sgen_heat"], label="переменная ячейка")
        if pu is not None:
            axs[2].plot(pu["x"] / pu["x"][-1], pu["sgen_heat"], label="постоянная ячейка", ls="--")
        axs[2].set_xlabel("x/L"); axs[2].set_ylabel("S_gen теплообмена, Вт/(К·м)"); axs[2].legend(fontsize=8)
        axs[2].grid(alpha=0.3)
        fig.suptitle("Гироид с переменной плотностью: что выбрал оптимизатор (ε = 0,90)", fontsize=10)
        fig.tight_layout(); fig.savefig(out_dir / "tpms_graded_profile.png", dpi=130); plt.close(fig)


def render_png(stl_path, png_path, title):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import trimesh
    m = trimesh.load(stl_path)
    fig = plt.figure(figsize=(8, 4.2))
    ax = fig.add_subplot(111, projection="3d")
    v, f = m.vertices, m.faces
    sel = f[:: max(1, len(f) // 60000)]
    ax.plot_trisurf(v[:, 0], v[:, 1], v[:, 2], triangles=sel, cmap="copper", linewidth=0, antialiased=False, shade=True)
    ax.set_box_aspect(np.ptp(v, axis=0))
    ax.set_axis_off()
    ax.view_init(elev=24, azim=-58)
    ax.set_title(title, fontsize=9)
    fig.tight_layout(); fig.savefig(png_path, dpi=140); plt.close(fig)


def main():
    global _MAC
    ap = argparse.ArgumentParser()
    ap.add_argument("--mac", required=True)
    ap.add_argument("--lbm", default=None)
    ap.add_argument("--out", default=str(HERE))
    ap.add_argument("--big", default=None, help="папка для полноразмерного STL (scratchpad)")
    ap.add_argument("--quick", action="store_true")
    ap.add_argument("--no-crosscheck", action="store_true")
    ap.add_argument("--geometry-only", default=None, help="готовый tpms_results.json: только STL и картинки")
    a = ap.parse_args()
    _MAC = a.mac
    if a.geometry_only:
        print(geometry_outputs(a.geometry_only, Path(a.out), a.big or a.out))
        return
    out = Path(a.out)
    t0 = time.time()
    mac = json.loads(Path(a.mac).read_text())
    lbm = json.loads(Path(a.lbm).read_text()) if a.lbm and Path(a.lbm).exists() else dict(rows=[])
    res = dict(title="Рекуператор на решётках TPMS (гироид, Schwarz-D) с переменной плотностью — KV-3",
               method=(R.__doc__.strip() + "\n\nПереносы: " + mac.get("method", "mac") + " — voxflow.py; инерция — lbm.py."),
               software=f"python {sys.version.split()[0]}, numpy {np.__version__}, scipy, scikit-image (marching cubes), "
                        f"cantera (свойства), pyamg (MINRES-предобуславливатель), trimesh (просмотр STL)",
               date="2026-11-14",
               status="оценка уровня облика: переносы — предел Стокса (сетка 32³) с поправкой LBM на инерцию для одного "
                      "уровня; нужен расчёт ячейки с инерцией и теплообменом на мелкой сетке (задача AWS) и испытание образца",
               inputs=[dict(file="routes/gtu/organic/cycles/cycles_results.json", field="results[R80|R85|R90 + RH].opt"),
                       dict(file=str(Path(a.mac).name), note="таблица переносов MAC (tpms_transport.py)"),
                       dict(file=str(Path(a.lbm).name) if a.lbm else None, note="LBM с инерцией (tpms_transport.py --method lbm)")])
    res["geometry"] = geometry_table()
    summ, lbm_s = transport_summary(mac["rows"], lbm["rows"])
    res["transport"] = dict(mac_summary=summ, lbm=lbm_s, plate=PLATE,
                            mac_run=dict(args=mac.get("args"), runtime_s=mac.get("runtime_s")))
    if not a.no_crosscheck:
        res["transport"]["crosscheck_mac_vs_lbm"] = crosscheck(0.0, 32)
        print("сверка", res["transport"]["crosscheck_mac_vs_lbm"], flush=True)
    inertia = inertia_from_lbm(lbm_s, summ)
    res["transport"]["inertia_used"] = inertia
    rec, profiles = recuperators(mac["rows"], inertia, quick=a.quick)
    res["recuperator"] = rec
    plots(res, profiles, out)
    res["runtime_s"] = round(time.time() - t0, 1)
    (out / "tpms_results.json").write_text(json.dumps(res, ensure_ascii=False, indent=1, default=float))
    geometry_outputs(out / "tpms_results.json", out, a.big or out)
    print("записано", res["runtime_s"])


def sheet_block(topo, a_mm, c, tw_mm, cells, voxel_mm, a_prof=None):
    """лист TPMS (металл) в блоке cells = (nx, ny, nz) ячеек; a_prof(s) — переменная ячейка по оси x (иллюстрация)"""
    from skimage import measure
    m = R.TPMSModel(topo, json.loads(Path(_MAC).read_text())["rows"])
    Lx, Ly, Lz = (k * a_mm for k in cells)
    nx, ny, nz = (int(round(v / voxel_mm)) + 1 for v in (Lx, Ly, Lz))
    x = np.linspace(0, Lx, nx)
    a = np.full(nx, a_mm) if a_prof is None else a_prof(x / Lx)
    phi = 2 * np.pi * np.concatenate([[0], np.cumsum(0.5 * (1 / a[1:] + 1 / a[:-1]) * np.diff(x))])
    y = np.linspace(0, Ly, ny); z = np.linspace(0, Lz, nz)
    t = m.t_level_vec(np.full(nx, c), tw_mm / a)
    F = field(topo, phi[:, None, None], 2 * np.pi * y[None, :, None] / a[:, None, None],
              2 * np.pi * z[None, None, :] / a[:, None, None])
    D = np.abs(F - c) - t[:, None, None]
    D = np.pad(D, 1, constant_values=1.0)
    v, f, _, _ = measure.marching_cubes(D, level=0.0, spacing=(voxel_mm,) * 3)
    return v, f


def render_cells(path, title, items):
    """светотеневой вид листов TPMS (2 × 2 × 2 ячейки)"""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import LightSource
    fig = plt.figure(figsize=(5.2 * len(items), 4.6))
    ls = LightSource(azdeg=315, altdeg=35)
    for k, (lab, v, f) in enumerate(items):
        ax = fig.add_subplot(1, len(items), k + 1, projection="3d")
        tri = v[f]
        n = np.cross(tri[:, 1] - tri[:, 0], tri[:, 2] - tri[:, 0])
        n /= np.maximum(np.linalg.norm(n, axis=1, keepdims=True), 1e-12)
        shade = ls.shade_normals(n, fraction=1.0)
        from mpl_toolkits.mplot3d.art3d import Poly3DCollection
        col = plt.cm.copper(0.25 + 0.7 * shade)
        pc = Poly3DCollection(tri, facecolors=col, edgecolors="none", linewidths=0)
        ax.add_collection3d(pc)
        lo, hi = v.min(0), v.max(0)
        ax.set_xlim(lo[0], hi[0]); ax.set_ylim(lo[1], hi[1]); ax.set_zlim(lo[2], hi[2])
        ax.set_box_aspect(hi - lo)
        ax.set_axis_off(); ax.view_init(elev=22, azim=-55)
        ax.set_title(lab, fontsize=9)
    fig.suptitle(title, fontsize=10)
    fig.tight_layout()
    fig.savefig(path, dpi=140)
    plt.close(fig)


def geometry_outputs(res_path, out, big):
    """STL оптимума (ε = 0,90) и картинки: вызывается после main или отдельно (--geometry-only)"""
    from tpms_geom import decimate_stl
    res = json.loads(Path(res_path).read_text())
    des = res["recuperator"]["cases"]["R90-ic-RH-ibc-none-none"]["designs"]
    stl = {}
    items = []
    for topo in ("gyroid", "diamond"):
        d = des.get(f"{topo}_graded") or des.get(f"{topo}_uniform")
        a_mm = d["z"][2] * 1e3
        c = d["z"][5] if len(d["z"]) == 8 else d["z"][3]
        tw = R.PRINT["t_w"] * 1e3
        v, f = sheet_block(topo, a_mm, c, tw, (2, 2, 2), a_mm / 44)
        items.append((f"{TOPOS[topo]}: ячейка {a_mm:.1f} мм, стенка {tw:.2f} мм, c = {c:+.3f}", v, f))
        if topo == "gyroid":
            full = Path(big) / "tpms_gyroid_optimum_block_full.stl"
            v2, f2 = sheet_block(topo, a_mm, c, tw, (4, 4, 8), 0.06)
            st0 = write_stl(full, v2, f2, solid="gyroid_opt")
            st = decimate_stl(full, Path(out) / "tpms_gyroid_optimum_block.stl")
            st.update(full=dict(path=str(full), triangles=st0["triangles"], bytes=st0["bytes"]),
                      note=f"лист гироида оптимума ε = 0,90 (стенка 0,35 мм): блок 4 × 4 × 8 ячеек = "
                           f"{4 * a_mm:.1f} × {4 * a_mm:.1f} × {8 * a_mm:.1f} мм; воксель 0,06 мм; в репозитории — "
                           f"упрощённая копия ≤ 5 МБ", path="tpms_gyroid_optimum_block.stl")
            stl["optimum_block"] = st
            # иллюстрация переменной ячейки (метод, не оптимум): ячейка 1:2 по длине
            gp = Path(big) / "tpms_gyroid_graded_illustration_full.stl"
            v3, f3 = sheet_block(topo, a_mm, c, tw, (12, 3, 3), 0.08, a_prof=R.pw([a_mm, 1.5 * a_mm, 2 * a_mm]))
            st3 = write_stl(gp, v3, f3, solid="gyroid_graded")
            st4 = decimate_stl(gp, Path(out) / "tpms_gyroid_graded_illustration.stl")
            st4.update(full=dict(path=str(gp), triangles=st3["triangles"], bytes=st3["bytes"]),
                       note="ИЛЛЮСТРАЦИЯ метода переменной ячейки (a от 1 до 2 оптимальных по длине), не оптимум: "
                            "в ламинарном пределе оптимизатор выбрал постоянную ячейку", path="tpms_gyroid_graded_illustration.stl")
            stl["graded_illustration"] = st4
    render_cells(Path(out) / "tpms_cells.png", "Лист TPMS оптимума рекуператора KV-3 (ε = 0,90): металл между двумя "
                 "лабиринтами — газ и воздух", items)
    res["stl"] = stl
    Path(res_path).write_text(json.dumps(res, ensure_ascii=False, indent=1, default=float))
    return stl


if __name__ == "__main__":
    main()
