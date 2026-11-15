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


# =========================================================================== доводка порога
LEVELS = (0.5, 0.45, 0.4, 0.35, 0.3, 0.25, 0.2, 0.15, 0.1, 0.05)


def pick_threshold(name, part, variant, h, u_target):
    """наибольший уровень x̄, при котором двоичный проект (x̄ ≥ уровня -> металл) выполняет ограничения (воксельный КЭ)"""
    model, design, cases, prm, meta, inb = DT.setup(part, variant, h)
    d = np.load(os.path.join(OUTDIR, f"disk_{name}.npz"))
    assert np.array_equal(d["ei"], model.ei) and np.array_equal(d["ej"], model.ej), "сетка прогона не совпадает"
    xb = d["xb"]
    scan = []
    best = None
    for lv in LEVELS:
        xe = np.where(xb >= lv, 1.0, 0.0); xe[model.group == 0] = 1.0
        sol = model.solve(xe, prm, cases, need_sens=False)
        eE, em, *_ = model.interp(xe, prm)
        u = [float(np.max(np.where(xe > 0.5, c["u"], 0.0))) for c in sol["cases"]]
        ok = all(u[k] <= u_target[k] * 1.01 for k in range(len(u)))
        scan.append(dict(level=lv, mass_kg=model.mass(em), J_kgm2=model.inertia(em), u=u, ok=ok))
        if ok and best is None:
            best = scan[-1]
    return dict(scan=scan, chosen=best)


# =========================================================================== диск турбины: профиль -> disk.py (CalculiX)
def turb_profile(name, n_pts=48):
    """профиль t(r), zc(r) из двоичного проекта: толщина металла по столбцу и его средняя линия; сглаживание по 3 столбцам"""
    D = OG.load(name)
    F = np.nan_to_num(D["F"], nan=0.0)
    nr, nz = F.shape
    dr, dz = D["r1"] / nr, (D["z1"] - D["z0"]) / nz
    zc = D["z0"] + (np.arange(nz) + 0.5) * dz
    S = F >= 0.5
    t = S.sum(axis=1) * dz
    zm = np.array([zc[S[i]].mean() if S[i].any() else 0.0 for i in range(nr)])
    k = np.ones(3) / 3
    t_s = np.convolve(np.r_[t[0], t, t[-1]], k, mode="valid"); zm_s = np.convolve(np.r_[zm[0], zm, zm[-1]], k, mode="valid")
    rc = (np.arange(nr) + 0.5) * dr
    r = np.r_[0.0, np.linspace(rc[0], rc[-1], n_pts - 2), D["r1"]]
    tt = np.interp(r, rc, t_s); zz = np.interp(r, rc, zm_s)
    return dict(r=r.tolist(), t=tt.tolist(), zc=zz.tolist())


def disk_case(name, prof, rpm, k_thick=1.0, mesh=(120, 16), t_cap=None, r_fix=None):
    """случай disk.py: профиль с утолщением k (кроме слоя под пазами r ≥ r_fix), поле cooled_750, давление лопаток ∝ ω²"""
    sys.path.insert(0, DT.TA.FEA)
    import disk as DK
    d = DT.frame3(); P = d["parts"]["turb1"]
    Tv = P["variants"]["cooled_750"]["T_C"]
    p_rim = P["variants"]["cooled_750"]["cases"]["100"]["rim_p_MPa"]
    r = np.asarray(prof["r"]); t = np.asarray(prof["t"]).copy()
    if k_thick != 1.0:
        m = r < (r_fix if r_fix is not None else r[-1])
        t[m] = t[m] * k_thick
        if t_cap:
            t = np.minimum(t, t_cap)
    case = dict(name=name, material="in713c", rpm=rpm, profile=dict(r=r.tolist(), t=t.tolist(), zc=prof.get("zc", [0.0] * len(r))),
                profile_interp="pchip", T_C=Tv, rim=dict(p_MPa=p_rim * (rpm / DT.RPM0) ** 2), mesh=dict(nr=mesh[0], nz=mesh[1]))
    wd = os.path.join(OUTDIR, "diskpy", f"{name}_{int(rpm)}")
    res = DK.run_case(case, workdir=wd, keep=True)
    # местный σв(T)/σ_экв по узлам .frd этого прогона (как frame3_fea: uts_local)
    frd = DK.parse_frd(os.path.join(wd, "disk.frd"))
    mat = DK.load_material("in713c"); ua = np.asarray(mat["uts_MPa"], float)
    S, TN = frd["STRESS"], frd.get("NDTEMP", {})
    rat = []
    for n, sv in S.items():
        vm = DK.von_mises(*sv[:6]) / 1e6
        if vm > 1.0:
            rat.append(float(np.interp(TN[n][0], ua[:, 0], ua[:, 1])) / vm)
    res["uts_local"] = dict(ratio=min(rat))
    rho = 7916.0
    rr = np.linspace(r[0], r[-1], 4001); tt = np.interp(rr, r, t)
    mass = float(np.trapezoid(rho * tt * 2 * math.pi * rr, rr)); J = float(np.trapezoid(rho * tt * 2 * math.pi * rr ** 3, rr))
    return dict(rpm=rpm, k=k_thick, mass_kg=mass, J_kgm2=J, vm_max_MPa=res["vm_max_node_MPa"], vm_max_at=res["vm_max_at"],
                uts_ratio=(res.get("uts_local") or {}).get("ratio"), creep_local_min_h=res["creep"]["local_min_h"],
                creep_local_s1_min_h=res["creep"]["local_s1_min_h"], creep_at=res["creep"]["local_at"],
                burst_margin=res["burst"]["margin"], gate=res["quality"]["passed"], flags=res["flags"],
                profile=dict(r=r.tolist(), t=t.tolist()))


def size_turb(name, prof, crit, mesh=(120, 16)):
    """наименьшее k ≥ 1 (бисекция), при котором: местный ресурс LC1 ≥ crit['life_h'] и σв/σ при 120 % ≥ crit['uts_ratio']"""
    d = DT.frame3(); g = d["parts"]["turb1"]["geometry"]
    r_fix = g["r_live_rim_m"] - 4e-3; t_cap = d["parts"]["turb1"]["profile_info"]["t_max_m"]

    def ok(k):
        a = disk_case(name, prof, DT.RPM0, k, mesh, t_cap, r_fix); b = disk_case(name, prof, 1.2 * DT.RPM0, k, mesh, t_cap, r_fix)
        good = a["creep_local_min_h"] >= crit["life_h"] and (b["uts_ratio"] or 0) >= crit["uts_ratio"]
        return good, a, b
    log = []
    g1, a, b = ok(1.0); log.append(dict(k=1.0, ok=g1, life=a["creep_local_min_h"], uts120=b["uts_ratio"], mass=a["mass_kg"]))
    if g1:
        return dict(k=1.0, c100=a, c120=b, log=log)
    lo, hi = 1.0, 1.6
    g2, a2, b2 = ok(hi); log.append(dict(k=hi, ok=g2, life=a2["creep_local_min_h"], uts120=b2["uts_ratio"], mass=a2["mass_kg"]))
    if not g2:
        return dict(k=None, c100=a2, c120=b2, log=log, note="при k = 1,6 критерии не выполнены")
    best = (hi, a2, b2)
    for _ in range(7):
        mid = 0.5 * (lo + hi)
        gm, am, bm = ok(mid)
        log.append(dict(k=mid, ok=gm, life=am["creep_local_min_h"], uts120=bm["uts_ratio"], mass=am["mass_kg"]))
        if gm:
            hi, best = mid, (mid, am, bm)
        else:
            lo = mid
    return dict(k=best[0], c100=best[1], c120=best[2], log=log)
