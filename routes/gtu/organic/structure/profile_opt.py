#!/usr/bin/env python3
"""Диск турбины ступени 1 KV-3: оптимизация формы гладкого профиля прямо в CalculiX (disk.py) — доводка «живого»
профиля после оптимизации топологии.

Почему отдельно: воксельная оптимизация топологии (disk_topo.py) дала профиль «ступица на всю высоту — тонкое полотно»
на 14 % легче литого, но на гладком профиле disk.py показал концентрацию на уступе ступицы и у перехода к ободу,
которую воксельная сетка не разрешает. Здесь профиль гладкий с самого начала: t(r) — PCHIP по узлам R_KNOTS,
ступица у оси — t_max (осевой габарит литого), слой под пазами и обод — w (как у литого), 6 переменных.
Цель — масса; ограничения (каждый шаг — два расчёта disk.py: 100 % и 120 %, поле cooled_750):
  eq:  местный ресурс по Мизесу ≥ ресурса литого профиля и местный σв/σ при 120 % ≥ литого;
  lim: местный ресурс ≥ 25 000 ч и σв/σ при 120 % ≥ 1.
Метод — COBYLA (scipy), старт — литой профиль frame3; итог перепроверяется на сетке 160×24.
"""
import argparse
import json
import math
import os
import sys
import time

import numpy as np
from scipy.optimize import minimize

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import disk_topo as DT  # noqa: E402
import structure_post as SP  # noqa: E402

R_KNOTS = (0.06, 0.09, 0.12, 0.15, 0.175, 0.19)


def build(x, g):
    r_b, r_live, t_max, w = g["r_b"], g["r_live"], g["t_max"], g["w"]
    r = np.r_[0.0, R_KNOTS, r_b, r_live]
    t = np.r_[t_max, np.clip(x, 0.012, t_max), w, w]
    return dict(r=r.tolist(), t=t.tolist(), zc=[0.0] * len(r))


def geom():
    d = DT.frame3(); P = d["parts"]["turb1"]
    return dict(r_b=P["profile_info"]["r_b_m"], r_live=P["geometry"]["r_live_rim_m"], t_max=P["profile_info"]["t_max_m"],
                w=P["geometry"]["w_rim_m"], base=P["profile"])


def evaluate(x, g, tag, mesh=(80, 12)):
    p = build(x, g)
    a = SP.disk_case(tag, p, DT.RPM0, mesh=mesh); b = SP.disk_case(tag, p, 1.2 * DT.RPM0, mesh=mesh)
    return p, a, b


def run(mode, maxiter=90):
    t0 = time.time()
    g = geom()
    base = dict(r=g["base"]["r"], t=g["base"]["t"], zc=[0.0] * len(g["base"]["r"]))
    b100 = SP.disk_case(f"po_{mode}_base", base, DT.RPM0, mesh=(80, 12)); b120 = SP.disk_case(f"po_{mode}_base", base, 1.2 * DT.RPM0, mesh=(80, 12))
    if mode == "eq":
        L_req, U_req = b100["creep_local_min_h"], b120["uts_ratio"]
    else:
        L_req, U_req = DT.CREEP_H, 1.0
    x0 = np.interp(R_KNOTS, g["base"]["r"], g["base"]["t"])
    hist = []
    cache = {}

    def ev(x):
        key = tuple(np.round(x, 6))
        if key not in cache:
            p, a, b = evaluate(x, g, f"po_{mode}")
            cache[key] = (p, a, b)
            hist.append(dict(x=list(map(float, x)), mass=a["mass_kg"], life=a["creep_local_min_h"], uts120=b["uts_ratio"],
                             J=a["J_kgm2"]))
            print(mode, len(hist), round(a["mass_kg"], 3), round(math.log10(max(a["creep_local_min_h"], 1e-9)), 3),
                  round(b["uts_ratio"], 4), flush=True)
        return cache[key]

    def f(x):
        return ev(x)[1]["mass_kg"] / b100["mass_kg"]

    cons = [dict(type="ineq", fun=lambda x: math.log10(max(ev(x)[1]["creep_local_min_h"], 1e-9)) - math.log10(L_req)),
            dict(type="ineq", fun=lambda x: ev(x)[2]["uts_ratio"] / U_req - 1.0)]
    for i in range(len(R_KNOTS)):
        cons.append(dict(type="ineq", fun=lambda x, i=i: x[i] - 0.012))
        cons.append(dict(type="ineq", fun=lambda x, i=i: g["t_max"] - x[i]))
    r = minimize(f, x0, method="COBYLA", constraints=cons, options=dict(maxiter=maxiter, rhobeg=0.006, catol=1e-4))
    # лучший допустимый из истории (COBYLA может закончить на слегка недопустимой точке)
    feas = [h for h in hist if h["life"] >= L_req * 0.999 and h["uts120"] >= U_req * 0.999]
    best = min(feas, key=lambda h: h["mass"]) if feas else None
    out = dict(mode=mode, L_req_h=L_req, U_req=U_req, base_80x12=dict(c100=b100, c120=b120), n_eval=len(hist),
               cobyla=dict(success=bool(r.success), message=str(r.message)), best=best, history=hist,
               seconds=round(time.time() - t0, 1))
    if best:
        p, a, b = evaluate(np.asarray(best["x"]), g, f"po_{mode}_fine", mesh=(160, 24))
        bb1 = SP.disk_case(f"po_{mode}_base_fine", base, DT.RPM0, mesh=(160, 24)); bb2 = SP.disk_case(f"po_{mode}_base_fine", base, 1.2 * DT.RPM0, mesh=(160, 24))
        out["check_160x24"] = dict(profile=p, c100=a, c120=b, base_c100=bb1, base_c120=bb2)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mode", choices=["eq", "lim", "comp2"], required=True)
    ap.add_argument("--maxiter", type=int, default=90)
    a = ap.parse_args()
    if a.mode == "comp2":
        run_comp2(a.maxiter)
        return
    out = run(a.mode, a.maxiter)
    out.update(title="Диск турбины ступени 1: оптимизация формы гладкого профиля в CalculiX (disk.py)", date="2026-11-14",
               method=__doc__.strip().splitlines()[0], software=dict(ccx=DT.TA.DK.ccx_version(), python=sys.version.split()[0]),
               inputs=[dict(file="routes/gtu/frame3/frame3_fea_results.json", field="parts.turb1 (profile, geometry, cooled_750)"),
                       dict(file="routes/gtu/platform/fea/materials.json", field="in713c")],
               knots_m=list(R_KNOTS))
    with open(os.path.join(HERE, f"profile_opt_{a.mode}.json"), "w", encoding="utf-8") as f:
        json.dump(out, f, ensure_ascii=False, indent=1, default=float)



# =========================================================================== колесо ступени 2: задняя поверхность
RC_KNOTS = (0.07, 0.11, 0.15, 0.20, 0.25, 0.30)


def comp2_case(tag, zb_knots, rpm, mesh=(80, 12)):
    """профиль колеса: тракт (верх) — как у кованого frame3, задняя поверхность — PCHIP по узлам RC_KNOTS"""
    sys.path.insert(0, DT.TA.FEA)
    import disk as DK
    from scipy.interpolate import PchipInterpolator
    d = DT.frame3(); P = d["parts"]["comp2"]; s = P["variants"]["solid"]; pr = s["profile"]
    r0 = np.asarray(pr["r"]); t0 = np.asarray(pr["t"]); zc0 = np.asarray(pr["zc"]); hb0 = np.asarray(pr["hb"])
    ztop = zc0 + t0 / 2; zbot0 = zc0 - t0 / 2
    kr = np.r_[0.0, 0.058, RC_KNOTS, r0[-1]]
    kz = np.r_[zbot0[0], zbot0[0], zb_knots, zbot0[-1]]
    zb = PchipInterpolator(kr, kz)(r0)
    zb = np.minimum(zb, ztop - 0.010)
    t = ztop - zb; zc = 0.5 * (ztop + zb)
    case = dict(name=tag, material="ti64", rpm=rpm, profile=dict(r=r0.tolist(), t=t.tolist(), zc=zc.tolist(), hb=hb0.tolist()),
                profile_interp="pchip", blade_fill=P["blade_fill"], T_C=s["T_C"], mesh=dict(nr=mesh[0], nz=mesh[1]))
    wd = os.path.join(SP.OUTDIR, "diskpy", f"{tag}_{int(rpm)}")
    res = DK.run_case(case, workdir=wd, keep=False)
    rr = np.linspace(0, r0[-1], 4001); tt = np.interp(rr, r0, t); hh = np.interp(rr, r0, hb0)
    mass = float(np.trapezoid(4430 * tt * 2 * math.pi * rr, rr)); J = float(np.trapezoid(4430 * tt * 2 * math.pi * rr ** 3, rr))
    return dict(mass_kg=mass, J_kgm2=J, yield_ratio=res["yield_ratio_max"], yield_at=res["yield_ratio_at"],
                burst_margin=res["burst"]["margin"], vm_max_MPa=res["vm_max_node_MPa"], gate=res["quality"]["passed"],
                profile=case["profile"])


def run_comp2(maxiter=70):
    t0 = time.time()
    d = DT.frame3(); pr = d["parts"]["comp2"]["variants"]["solid"]["profile"]
    r0 = np.asarray(pr["r"]); zb0 = np.asarray(pr["zc"]) - np.asarray(pr["t"]) / 2
    x0 = np.interp(RC_KNOTS, r0, zb0)
    b1 = comp2_case("pc2_base", x0, DT.RPM0); b2 = comp2_case("pc2_base", x0, 1.2 * DT.RPM0)
    Y1, Y2 = b1["yield_ratio"], b2["yield_ratio"]
    hist, cache = [], {}

    def ev(x):
        key = tuple(np.round(x, 6))
        if key not in cache:
            a = comp2_case("pc2", x, DT.RPM0); b = comp2_case("pc2", x, 1.2 * DT.RPM0)
            cache[key] = (a, b)
            hist.append(dict(x=list(map(float, x)), mass=a["mass_kg"], J=a["J_kgm2"], y100=a["yield_ratio"], y120=b["yield_ratio"]))
            print("c2", len(hist), round(a["mass_kg"], 2), round(a["yield_ratio"], 4), round(b["yield_ratio"], 4), flush=True)
        return cache[key]
    zmin = float(zb0.min())
    cons = [dict(type="ineq", fun=lambda x: 1 - ev(x)[0]["yield_ratio"] / Y1),
            dict(type="ineq", fun=lambda x: 1 - ev(x)[1]["yield_ratio"] / Y2)]
    for i in range(len(RC_KNOTS)):
        cons.append(dict(type="ineq", fun=lambda x, i=i: x[i] - zmin))
    r = minimize(lambda x: ev(x)[0]["mass_kg"] / b1["mass_kg"], x0, method="COBYLA", constraints=cons,
                 options=dict(maxiter=maxiter, rhobeg=0.008, catol=1e-4))
    feas = [h for h in hist if h["y100"] <= Y1 * 1.001 and h["y120"] <= Y2 * 1.001]
    best = min(feas, key=lambda h: h["mass"]) if feas else None
    out = dict(title="Колесо ступени 2: оптимизация гладкой задней поверхности в CalculiX (disk.py)", date="2026-11-14",
               method="задняя поверхность — PCHIP по узлам RC_KNOTS, тракт и обод — как у кованого профиля frame3; масса лопаток — "
                      "в плотности столбца (как frame3); цель — масса; ограничения — наибольшее σ/σ0,2 при 100 и 120 % не выше "
                      "кованого (k = 1); COBYLA; проверка лучшего на сетке 160×24",
               software=dict(ccx=DT.TA.DK.ccx_version(), python=sys.version.split()[0]), knots_m=list(RC_KNOTS),
               inputs=[dict(file="routes/gtu/frame3/frame3_fea_results.json", field="parts.comp2.variants.solid (profile, T_C), blade_fill"),
                       dict(file="routes/gtu/platform/fea/materials.json", field="ti64")],
               base=dict(c100=b1, c120=b2), n_eval=len(hist),
               cobyla=dict(success=bool(r.success), message=str(r.message)), best=best, history=hist)
    if best:
        a = comp2_case("pc2_fine", np.asarray(best["x"]), DT.RPM0, mesh=(160, 24))
        b = comp2_case("pc2_fine", np.asarray(best["x"]), 1.2 * DT.RPM0, mesh=(160, 24))
        bb1 = comp2_case("pc2_base_fine", x0, DT.RPM0, mesh=(160, 24)); bb2 = comp2_case("pc2_base_fine", x0, 1.2 * DT.RPM0, mesh=(160, 24))
        out["check_160x24"] = dict(c100=a, c120=b, base_c100=bb1, base_c120=bb2)
    out["seconds"] = round(time.time() - t0, 1)
    with open(os.path.join(HERE, "profile_opt_comp2.json"), "w", encoding="utf-8") as f:
        json.dump(out, f, ensure_ascii=False, indent=1, default=float)
    return out


if __name__ == "__main__":
    main()
