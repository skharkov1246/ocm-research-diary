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
            case["designs"][name] = strip(r)
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


def graded_block_stl(path, topo, a_prof, c_prof, t_w, size_mm, voxel_mm):
    """образец листа TPMS с переменной ячейкой по длине: фаза φ(x) = 2π∫dx/a(x), поперечные координаты
    масштабированы местной ячейкой; лист |F − c(x)| < t(x) (t — из таблицы для t_w/a)"""
    from skimage import measure
    Lx, Ly, Lz = size_mm
    nx, ny, nz = (int(round(v / voxel_mm)) + 1 for v in size_mm)
    x = np.linspace(0, Lx, nx)
    a = a_prof(x / Lx)
    phi = 2 * np.pi * np.concatenate([[0], np.cumsum(0.5 * (1 / a[1:] + 1 / a[:-1]) * np.diff(x))])
    y = np.linspace(0, Ly, ny)
    z = np.linspace(0, Lz, nz)
    m = R.TPMSModel(topo, json.loads(Path(_MAC).read_text())["rows"])
    c = c_prof(x / Lx)
    t = m.t_level_vec(c, t_w / a)
    PHI = phi[:, None, None]
    Y = 2 * np.pi * y[None, :, None] / a[:, None, None]
    Z = 2 * np.pi * z[None, None, :] / a[:, None, None]
    F = field(topo, PHI, Y, Z)
    D = np.abs(F - c[:, None, None]) - t[:, None, None]          # < 0 — металл
    D = np.pad(D, 1, constant_values=1.0)                          # замкнуть тело по краям
    v, f, _, _ = measure.marching_cubes(D, level=0.0, spacing=(voxel_mm,) * 3)
    return write_stl(path, v, f, solid=f"{topo}_graded")


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
    labels = {"plate_NuH1": "пластины", "gyroid_uniform": "гироид, ячейка пост.", "gyroid_graded": "гироид, переменная",
              "diamond_uniform": "Schwarz-D, пост.", "diamond_graded": "Schwarz-D, переменная",
              "gyroid_graded_tw025": "гироид, стенка 0,25", "gyroid_graded_inertia": "гироид + инерция (LBM)"}
    for ax, field_, yl in ((axs[0], "mass_core", "масса сердцевины, т"), (axs[1], "V", "объём сердцевины, м³"),
                           (axs[2], "A_fr", "фронт, м²")):
        names = [n for n in labels if n in rc.get("R90-ic-RH-ibc-none-none", {}).get("designs", {})]
        vals = [rc["R90-ic-RH-ibc-none-none"]["designs"][n][field_] / (1e3 if field_ == "mass_core" else 1) for n in names]
        cols = ["0.55"] + ["#3b7dd8"] * (len(names) - 1)
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
    a = ap.parse_args()
    _MAC = a.mac
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
    # образец переменной решётки по найденному профилю (гироид, ε = 0,90)
    d = rec["cases"]["R90-ic-RH-ibc-none-none"]["designs"].get("gyroid_graded")
    stl = {}
    if d:
        z = d["z"]
        fa = R.pw(np.array(z[2:5]) * 1e3)
        fc = R.pw(z[5:8])
        from tpms_geom import decimate_stl
        big_dir = Path(a.big) if a.big else out
        full = big_dir / "tpms_graded_demo_full.stl"
        st0 = graded_block_stl(full, "gyroid", fa, fc, R.PRINT["t_w"] * 1e3, (60.0, 20.0, 20.0), 0.12)
        st = decimate_stl(full, out / "tpms_graded_demo.stl")
        st["full"] = dict(path=str(full), triangles=st0["triangles"], bytes=st0["bytes"])
        st["note"] = ("образец 60 × 20 × 20 мм: профиль a(x), c(x) оптимума ε = 0,90 сжат по длине в 60 мм; мм; "
                      "воксель 0,12 мм (стенка 0,35 мм), в репозитории — упрощённая копия ≤ 5 МБ")
        st["path"] = str(Path(st["path"]).relative_to(ROOT)) if str(st["path"]).startswith(str(ROOT)) else st["path"]
        stl["demo"] = st
        render_png(out / "tpms_graded_demo.stl", out / "tpms_graded_demo.png",
                   "Гироид с переменной ячейкой (лист — металл, два лабиринта — газ и воздух), образец 60 мм")
        if a.big:
            big = Path(a.big) / "tpms_graded_module_section.stl"
            L = d["L"] * 1e3
            stb = graded_block_stl(big, "gyroid", fa, fc, R.PRINT["t_w"] * 1e3, (L, 20.0, 20.0), 0.2)
            stb["note"] = f"сечение 20 × 20 мм на всю длину сердцевины {L:.0f} мм — в scratchpad (большой файл)"
            stl["full_length_section"] = stb
    res["stl"] = stl
    res["runtime_s"] = round(time.time() - t0, 1)
    (out / "tpms_results.json").write_text(json.dumps(res, ensure_ascii=False, indent=1, default=float))
    print("записано", res["runtime_s"])


if __name__ == "__main__":
    main()
