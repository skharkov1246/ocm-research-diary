#!/usr/bin/env python3
"""Постобработка узла organic/structure: проверка проектов в CalculiX на гладком контуре, геометрия (STL), рендеры,
рисунки и сводный results-JSON.

Проверка в CalculiX (ccx_check.py): сечение — многоугольники контура x̄ = 0,5 с точной огибающей (тракт, обод),
сетка Gmsh CAX6 шага h_ccx, нагрузки и допускаемые — те же, что в оптимизации (disk_topo.setup); исходный профиль
frame3 считается той же цепочкой (сравнение «один к одному»). Решётка — по группам ρ* (шаг 0,01).

Запуск: python3 structure_post.py [--ccx] [--geom] [--figs] [--summary]   (без ключей — всё)
"""
import argparse
import json
import math
import os
import sys
import time
from multiprocessing import Pool

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import disk_topo as DT  # noqa: E402
import organic_geom as OG  # noqa: E402
import lattice as LT  # noqa: E402

OUTDIR = OG.OUTDIR
CCX_RES = os.path.join(HERE, "ccx_check_results.json")
H_CCX = {"turb1": 1.0e-3, "comp2": 1.5e-3}


def part_loads(part):
    """нагрузки и допускаемые для ccx_check — из той же постановки, что оптимизация"""
    d = DT.frame3()
    P = d["parts"][part]
    if part == "turb1":
        g = P["geometry"]
        Tv = P["variants"]["cooled_750"]["T_C"]
        T_of_r = lambda r: np.interp(r, Tv["r"], Tv["T"])  # noqa: E731
        p = P["variants"]["cooled_750"]["cases"]["100"]["rim_p_MPa"] * 1e6
        rim = dict(r=g["r_live_rim_m"], zlo=-g["w_rim_m"] / 2, zhi=g["w_rim_m"] / 2, p_Pa=p)
        return dict(mat="in713c", T_of_r=T_of_r, rim=rim, blade=None, fixz=(0.0, 0.0))
    s = P["variants"]["solid"]
    pr = s["profile"]
    r = np.asarray(pr["r"]); t = np.asarray(pr["t"]); zc = np.asarray(pr["zc"]); hb = np.asarray(pr["hb"])
    ztop = zc + t / 2
    Tv = s["T_C"]
    T_of_r = lambda rr: np.interp(rr, Tv["r"], Tv["T"])  # noqa: E731
    phi = P["blade_fill"]

    def f_r(rr, rpm):
        w2 = (2 * math.pi * rpm / 60) ** 2
        return phi * 4430.0 * w2 * np.interp(rr, r, hb) * 2 * math.pi * rr ** 2

    def on_hub(a, c):
        za, zc_ = np.interp(a[0], r, ztop), np.interp(c[0], r, ztop)
        return abs(a[1] - za) < 2e-5 and abs(c[1] - zc_) < 2e-5 and max(a[0], c[0]) > r[0] + 1e-9
    return dict(mat="ti64", T_of_r=T_of_r, rim=None, f_r=f_r, on_hub=on_hub, fixz=(0.0, float((zc - t / 2).min())))


def allow_fn(part, k_am, lc):
    import topo_axi as TA
    mat = TA.DK.load_material("in713c" if part == "turb1" else "ti64")
    mt = TA.mat_tables(mat)
    if part == "turb1":
        if lc == 0:
            return lambda T, rho: np.minimum(mt["slm"](T, DT.CREEP_H), mt["uts"](T)) * 1e6
        return lambda T, rho: mt["uts"](T) * 1e6
    if lc == 0:
        return lambda T, rho: k_am * mt["ys"](T) * 1e6 / 1.5
    return lambda T, rho: k_am * mt["ys"](T) * 1e6


def ccx_job(job):
    import ccx_check as CC
    name, part, polys, k_am, rho_field_npz = job["name"], job["part"], job["polys"], job["k_am"], job.get("rho_npz")
    L = part_loads(part)
    lat = None; rho_field = None
    if rho_field_npz:
        D = OG.load(rho_field_npz)
        lat = LT.default_law()
        F = np.nan_to_num(D["F"], nan=1.0)
        nr, nz = F.shape
        dr, dz = D["r1"] / nr, (D["z1"] - D["z0"]) / nz
        G = D["G"]
        rmin = lat["rho_min"]

        def rho_field(rr, zz):
            i = np.clip((rr / dr).astype(int), 0, nr - 1); j = np.clip(((zz - D["z0"]) / dz).astype(int), 0, nz - 1)
            x = F[i, j]
            return np.where(G[i, j] == 2, rmin + (1 - rmin) * x, 1.0)
    out = dict(name=name, part=part, cases=[])
    for lc, rpm in enumerate((DT.RPM0, 1.2 * DT.RPM0)):
        blade = None
        if L.get("f_r"):
            blade = dict(f=lambda rr, rpm=rpm: L["f_r"](rr, rpm), select=L["on_hub"])
        rim = None
        if L.get("rim"):
            rim = dict(L["rim"]); rim["p_Pa"] = rim["p_Pa"] * (rpm / DT.RPM0) ** 2
        res = CC.check(f"{name}_lc{lc + 1}", polys, L["mat"], L["T_of_r"], rpm, rim=rim, blade=blade, lattice=lat,
                       rho_field=rho_field, h=H_CCX[part], fix_z_point=L["fixz"], allow=allow_fn(part, k_am, lc),
                       workdir=os.path.join(OUTDIR, "ccx", f"{name}_lc{lc + 1}"))
        out["cases"].append(res)
    return out


def polys_of(name):
    D = OG.load(name)
    gs = OG.section_polygons(D, smooth_tol=0.3e-3)
    return [r.tolist() for g in gs for r in OG.poly_rings(g)] if len(gs) == 1 else None, gs, D


def run_ccx(names, workers=4):
    with open(DT.RESULTS, encoding="utf-8") as f:
        R = json.load(f)
    runs = {r["name"]: r for r in R["runs"]}
    jobs = []
    for part in ("turb1", "comp2"):
        model, design, cases, prm, meta, inb = DT.setup(part, "A", 1e-3 if part == "turb1" else 1.5e-3)
        base = np.asarray(meta["base_poly"])
        jobs.append(dict(name=f"{part}_baseline", part=part, polys=[base.tolist()], k_am=1.0))
    for nm in names:
        r = runs[nm]
        if r["variant"] == "A":
            pl, gs, D = polys_of(nm)
            if pl is None:
                print("несвязное сечение:", nm, len(gs)); continue
            jobs.append(dict(name=nm, part=r["part"], polys=pl, k_am=r["meta"]["k_am"]))
        else:
            D = OG.load(nm)
            jobs.append(dict(name=nm, part=r["part"], polys=[np.asarray(D["env"]).tolist()], k_am=r["meta"]["k_am"],
                             rho_npz=nm))
    with Pool(workers) as p:
        res = p.map(ccx_job, jobs, chunksize=1)
    out = dict(title="Проверка проектов в CalculiX на гладком контуре (CAX6, Gmsh)", date="2026-11-14",
               method="ccx_check.py: контур x̄ = 0,5 ∩ огибающая -> Gmsh CAX6 -> ccx 2.21; CENTRIF, *TEMPERATURE T(r), *CLOAD "
                      "(полная сила по окружности — соглашение проверено cload_convention: Ламе, 0,015 %); допускаемые — как в "
                      "оптимизации; решётка — группы ρ* с E·E*/E и ρ·ρ*, загрузка решётки против её прочности S(ρ*)",
               software=dict(ccx=res[0]["cases"][0]["ccx"] if res else None, python=sys.version.split()[0]),
               h_ccx_m=H_CCX, results=res)
    with open(CCX_RES, "w", encoding="utf-8") as f:
        json.dump(out, f, ensure_ascii=False, indent=1, default=float)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ccx", action="store_true")
    ap.add_argument("--names", default="")
    ap.add_argument("--workers", type=int, default=4)
    a = ap.parse_args()
    if a.ccx:
        names = [n for n in a.names.split(",") if n]
        run_ccx(names, a.workers)


if __name__ == "__main__":
    main()
