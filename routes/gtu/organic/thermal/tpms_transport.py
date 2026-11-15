"""Коэффициенты переноса трактов TPMS (гироид, Schwarz-D) по уровню поверхности ℓ (тракт F > ℓ).

Два независимых метода (правило «каждая правка данных проверяется минимум дважды разными способами»):
  mac — voxflow.PeriodicVoxelFlow: предел Стокса, f·Re и Nu_H1(Pe) при стоксовой кинематике;
  lbm — lbm.LBMCell: D3Q19/D3Q7 TRT, с инерцией (конечные Re) — f·Re(Re) и Nu_H1(Re, Pr = 0,7).
Воксели — с подбором порога по истинной пористости (tpms_geom.voxel_mask) и сдвигом сетки.
Определения (как у плоского канала в organic/cycles): Dh = 4ε/β по истинной (marching cubes) поверхности,
U — средняя межпоровая скорость u_s/ε, f — Фаннинг: f·Re = (dp/dx)·Dh²/(2μU), Nu = q_w·Dh/(k·(T_w − T_b)).

Запуск:  python tpms_transport.py --method mac --n 40 56 --workers 4 --out tpms_transport_mac.json
         python tpms_transport.py --method lbm --n 48 --re 10 50 150 400 800 --workers $NPROC --out ...
"""
import argparse
import json
import os
import sys
import time
from multiprocessing import Pool
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from tpms_geom import side_props, voxel_mask, TOPOS  # noqa: E402
from voxflow import PeriodicVoxelFlow, percolating  # noqa: E402

LEVELS = (-0.4, -0.2, 0.0, 0.2, 0.4, 0.6, 0.8)
_GEOM = {}


def geom(topo, level):
    k = (topo, round(level, 6))
    if k not in _GEOM:
        _GEOM[k] = side_props(topo, level, 200)
    return _GEOM[k]


def mac_point(args):
    topo, level, n, shift, pes = args
    g = geom(topo, level)
    t0 = time.time()
    m, thr = voxel_mask(topo, level, n, shift, eps_target=g["eps"])
    m = percolating(m)
    s = PeriodicVoxelFlow(m)
    info = s.stokes(solver="auto" if n <= 24 else "minres")
    us = float(s.U[0].sum() / m.size)              # поверхностная скорость, g = μ = 1
    Dh = g["Dh"] * n
    A = g["beta"] * n ** 2
    U = us / g["eps"]
    fRe = Dh ** 2 / (2 * U)
    nus = []
    for Pe in pes:
        sc = Pe / (U * Dh)
        r = s.scalar_h1(sc, scheme="central")
        Nu = r["Q_balance"] * Dh / (A * (-r["phib"]))
        nus.append(dict(Pe=Pe, Nu=float(Nu), wall_balance=float(r["Q_wall"] / r["Q_balance"])))
    return dict(method="mac", topo=topo, level=level, n=n, shift=list(shift), thr=thr, eps_vox=float(m.mean()),
                eps=g["eps"], beta=g["beta"], Dh=g["Dh"], k_a2=us / n ** 2, fRe=float(fRe), Nu=nus,
                stokes=info, div_max=s.divergence_max(), time_s=round(time.time() - t0, 1))


def lbm_point(args):
    """одна точка с инерцией: заданное Re (по Dh и межпоровой скорости), Pr = 0,7"""
    topo, level, n, shift, Re, Pr, U_lu, max_steps = args
    from lbm import LBMCell
    g = geom(topo, level)
    t0 = time.time()
    m, thr = voxel_mask(topo, level, n, shift, eps_target=g["eps"])
    m = percolating(m)
    Dh = g["Dh"] * n
    A = g["beta"] * n ** 2
    nu = U_lu * Dh / Re
    nu = max(nu, 0.0034)                          # τ+ ≥ 0,51 (устойчивость TRT, оценка)
    # сила по стоксовой оценке с поправкой на инерцию подбирается в ходе счёта
    F = U_lu * nu * 2 * 30.0 / Dh ** 2            # f·Re ≈ 30 — начальная оценка
    L = LBMCell(m, nu, F, kappa=nu / Pr)
    hist = []
    for rep in range(6):
        r = L.run(max_steps=max_steps, check=1000, tol=2e-7, min_steps=4000)
        ob = r["obs"]
        Uc = ob["U"] * m.mean() / g["eps"]
        hist.append(dict(steps=L.it, U=Uc, converged=r["converged"]))
        if abs(Uc / U_lu - 1) < 0.15 or not r["converged"]:
            break
        L.F *= U_lu / Uc                            # подогнать силу под нужную скорость (Re)
    ob = L.observables()
    U = ob["U"] * m.mean() / g["eps"]
    Re_act = U * Dh / nu
    fRe = L.F * Dh ** 2 / (2 * nu * U)
    kappa = nu / Pr
    Nu = ob["Q"] * Dh / (A * kappa * (-ob["phib"]))
    return dict(method="lbm", topo=topo, level=level, n=n, shift=list(shift), eps=g["eps"], beta=g["beta"],
                Dh=g["Dh"], Re_target=Re, Re=float(Re_act), Pr=Pr, nu_lu=nu, U_lu=float(U), Ma=float(U * np.sqrt(3)),
                umax=ob["umax"], fRe=float(fRe), Nu=float(Nu), steps=L.it, hist=hist, rho_range=[ob["rho_min"], ob["rho_max"]],
                time_s=round(time.time() - t0, 1))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--method", choices=("mac", "lbm"), default="mac")
    ap.add_argument("--topos", nargs="+", default=list(TOPOS))
    ap.add_argument("--levels", nargs="+", type=float, default=list(LEVELS))
    ap.add_argument("--n", nargs="+", type=int, default=[40, 56])
    ap.add_argument("--shifts", nargs="+", type=float, default=[0.0, 0.5])
    ap.add_argument("--pe", nargs="+", type=float, default=[0.01, 10.0, 20.0, 40.0])
    ap.add_argument("--pe-only-finest", action="store_true", default=True)
    ap.add_argument("--re", nargs="+", type=float, default=[10.0, 50.0, 150.0])
    ap.add_argument("--pr", type=float, default=0.7)
    ap.add_argument("--u-lu", type=float, default=0.04)
    ap.add_argument("--max-steps", type=int, default=60000)
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--out", default=str(HERE / "tpms_transport_mac.json"))
    a = ap.parse_args()
    tasks = []
    if a.method == "mac":
        nmax = max(a.n)
        for topo in a.topos:
            for lv in a.levels:
                for n in a.n:
                    for s in a.shifts:
                        pes = a.pe if (n == nmax and s == a.shifts[0]) or not a.pe_only_finest else [a.pe[0]]
                        tasks.append((topo, lv, n, (s, s, s), tuple(pes)))
        fn = mac_point
    else:
        for topo in a.topos:
            for lv in a.levels:
                for n in a.n:
                    for Re in a.re:
                        tasks.append((topo, lv, n, (0.0, 0.0, 0.0), Re, a.pr, a.u_lu, a.max_steps))
        fn = lbm_point
    # геометрию считаем заранее (один раз на точку)
    for t in tasks:
        geom(t[0], t[1])
    t0 = time.time()
    with Pool(a.workers) as pool:
        rows = []
        for r in pool.imap_unordered(fn, tasks):
            rows.append(r)
            print(json.dumps({k: v for k, v in r.items() if k not in ("hist",)}, ensure_ascii=False), flush=True)
    rows.sort(key=lambda r: (r["topo"], r["level"], r["n"], str(r["shift"]), r.get("Re_target", 0)))
    out = dict(method=a.method, software=f"python {sys.version.split()[0]}, numpy {np.__version__}",
               date="2026-11-14", args=vars(a), runtime_s=round(time.time() - t0, 1), rows=rows)
    Path(a.out).write_text(json.dumps(out, ensure_ascii=False, indent=1))
    print("записано", a.out)


if __name__ == "__main__":
    main()
