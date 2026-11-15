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
    ap.add_argument("--mode", choices=["eq", "lim"], required=True)
    ap.add_argument("--maxiter", type=int, default=90)
    a = ap.parse_args()
    out = run(a.mode, a.maxiter)
    out.update(title="Диск турбины ступени 1: оптимизация формы гладкого профиля в CalculiX (disk.py)", date="2026-11-14",
               method=__doc__.strip().splitlines()[0], software=dict(ccx=DT.TA.DK.ccx_version(), python=sys.version.split()[0]),
               inputs=[dict(file="routes/gtu/frame3/frame3_fea_results.json", field="parts.turb1 (profile, geometry, cooled_750)"),
                       dict(file="routes/gtu/platform/fea/materials.json", field="in713c")],
               knots_m=list(R_KNOTS))
    with open(os.path.join(HERE, f"profile_opt_{a.mode}.json"), "w", encoding="utf-8") as f:
        json.dump(out, f, ensure_ascii=False, indent=1, default=float)


if __name__ == "__main__":
    main()
