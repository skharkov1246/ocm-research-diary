"""Проверки узла organic/thermal на аналитике и балансах (быстро, без внешних данных, кроме frame3 и cycles).
python3 -m pytest -q test_thermal.py
"""
import json
import math
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
sys.path.insert(0, str(HERE))


def test_voxflow_plates_second_order():
    from voxflow import PeriodicVoxelFlow
    errs = []
    for n in (8, 16):
        m = np.ones((4, n + 1), bool); m[:, 0] = False
        s = PeriodicVoxelFlow(m); s.stokes(solver="direct")
        fRe = (2 * n) ** 2 / (2 * s.mean_velocity())
        h = s.scalar_h1(1.0)
        Nu = h["Q_wall"] * 2 * n / (8 * (-h["phib"]))
        errs.append((abs(fRe / 24 - 1), abs(Nu / 8.2355 - 1)))
        assert abs(h["Q_wall"] / h["Q_balance"] - 1) < 1e-10        # баланс тепла точный
    assert errs[1][0] < 0.01 and errs[1][1] < 0.005
    assert errs[0][0] / errs[1][0] > 3.5                            # второй порядок


def test_voxflow_square_duct_richardson():
    from voxflow import PeriodicVoxelFlow
    out = {}
    for n in (16, 24):
        m = np.ones((3, n + 1, n + 1), bool); m[:, 0, :] = False; m[:, :, 0] = False
        s = PeriodicVoxelFlow(m); s.stokes(solver="direct")
        h = s.scalar_h1(1.0)
        out[n] = (n ** 2 / (2 * s.mean_velocity()), h["Q_wall"] * n / (4 * n * 3 * (-h["phib"])))
    f = (24 ** 2 * out[24][0] - 16 ** 2 * out[16][0]) / (24 ** 2 - 16 ** 2)
    nu = (24 ** 2 * out[24][1] - 16 ** 2 * out[16][1]) / (24 ** 2 - 16 ** 2)
    assert abs(f / 14.227 - 1) < 0.002 and abs(nu / 3.608 - 1) < 0.003


def test_minres_equals_direct():
    pytest.importorskip("pyamg")
    from voxflow import PeriodicVoxelFlow, percolating
    from tpms_geom import grid
    m = percolating(grid("gyroid", 16) > 0.3)
    s1 = PeriodicVoxelFlow(m); s1.stokes(solver="direct")
    s2 = PeriodicVoxelFlow(m); s2.stokes(solver="minres", tol=1e-12)
    assert abs(s2.mean_velocity() / s1.mean_velocity() - 1) < 1e-6


def test_lbm_plates():
    from lbm import LBMCell
    n = 12
    m = np.ones((2, n + 1, 1), bool); m[:, 0, :] = False
    nu, F = 1 / 6, 1e-6
    L = LBMCell(m, nu, F, kappa=nu / 0.7)
    r = L.run(max_steps=40000, check=200, tol=1e-10)
    ob = r["obs"]
    fRe = F * (2 * n) ** 2 / (2 * nu * ob["U"])
    Nu = ob["Q"] * 2 * n / (4 * (nu / 0.7) * (-ob["phib"]))
    assert abs(fRe / 24 - 1) < 0.01 and abs(Nu / 8.2355 - 1) < 0.01


def test_tpms_areas():
    from tpms_geom import surface_density
    assert abs(surface_density("gyroid", 0.0, 160) / 3.0919 - 1) < 0.003
    assert abs(surface_density("diamond", 0.0, 160) / 3.8379 - 1) < 0.003


def test_murray_delta_tree():
    from vascular import net_delta
    for ex in (3.0, 7.0 / 3.0):
        t = net_delta(0.05, 0.03, 3, 1e-3, ex, "t", n_root=2)
        # у каждого узла: d_parent^ex = Σ d_child^ex
        P = {}
        for a, b, d in t.seg:
            P.setdefault(a, []).append(d)
        man = set(t.manifold_segs)
        parent_d = {b: d for i, (a, b, d) in enumerate(t.seg) if i not in man}
        for node, ch in P.items():
            if node in parent_d and len(ch) == 2:
                assert abs(sum(c ** ex for c in ch) / parent_d[node] ** ex - 1) < 1e-9
        assert len(t.outlets) == 16


def test_vascular_energy_balance_and_monotone():
    from vascular import kv3_point, ngv_case, net_straight, evaluate
    pt = kv3_point()
    panel, p_in, p_out, N_v = ngv_case(pt, 1100 + 273.15)
    s = net_straight(panel.Lx, panel.H, 8, 0.6e-3)
    r1 = evaluate(panel, s, p_in, p_out, pt["T3"], outer=2)
    assert abs(r1["Q_gas_W"] / r1["Q_cool_W"] - 1) < 1e-6
    r2 = evaluate(panel, s, 0.5 * (p_in + p_out(panel.Lx)), p_out, pt["T3"], outer=2)
    assert r2["m_c"] < r1["m_c"] and r2["T_max"] > r1["T_max"]          # меньше расход — горячее
    assert r1["S_gen_wall_to_coolant_W_K"] > 0 and r1["S_gen_friction_W_K"] > 0


def test_brinkman2d_poiseuille():
    """пустой канал с параболическим входом и выходом: Δp = 12μUL/H² (плоский Пуазейль)"""
    from topo import Brinkman2D
    nx, ny = 60, 20
    h = 1.0 / ny
    ops = [dict(side="W", j0=0, j1=ny, kind="velocity", value=1.0, name="вход"),
           dict(side="E", j0=0, j1=ny, kind="velocity", value=-1.0, name="выход")]
    S = Brinkman2D(nx, ny, h, ops)
    u, p = S.solve(np.zeros(nx * ny))
    P = p.reshape(nx, ny)
    dpdx = (P[nx // 4].mean() - P[3 * nx // 4].mean()) / (0.5 * nx * h)
    assert abs(dpdx / 12.0 - 1) < 0.01
    Phi = S.dissipation(u, np.zeros(nx * ny))
    assert abs(Phi / (12.0 * nx * h) - 1) < 0.03                       # Φ = Δp·Q на единицу глубины


def test_counterflow_eps_ntu():
    """дискретный противоток с постоянными свойствами против ε-NTU"""
    import tpms_recup as R

    class Const:
        def __init__(self, cp):
            self.cp = cp
        def __call__(self, T):
            T = np.asarray(T, float)
            return [np.full_like(T, self.cp), np.full_like(T, 3e-5), np.full_like(T, 0.05), np.full_like(T, 1.0),
                    self.cp * T]
        def h(self, T):
            return self.cp * np.asarray(T, float)
    du = dict(m_h=1.0, m_c=1.0, Th_in=900.0, Th_out_cycle=500.0, Tc_in=400.0, Tc_out=800.0)
    PH, PC = Const(1100.0), Const(1000.0)
    r = R.core_eval(du, PH, PC, "plate", dict(gap=1e-3, Nu=8.2355), A_fr=0.05, L=0.2, N=200)
    pitch = 2e-3 + 2 * R.PRINT["t_w"]
    beta = 2 / pitch
    hh = 8.2355 * 0.05 / 2e-3
    U = 1 / (2 / hh + R.PRINT["t_w"] / (R.PRINT["k_w"]))
    UA = U * beta * 0.05 * 0.2
    Cmin, Cmax = 1000.0, 1100.0
    NTU, Cr = UA / Cmin, Cmin / Cmax
    eps = (1 - math.exp(-NTU * (1 - Cr))) / (1 - Cr * math.exp(-NTU * (1 - Cr)))
    assert abs(r["Q_kW"] * 1e3 / (eps * Cmin * 500.0) - 1) < 0.003


def test_jobs_plan():
    jf = HERE / "jobs_thermal.json"
    d = json.loads(jf.read_text())
    for j in d["jobs"]:
        for tok in j["command"].split():
            if tok.startswith("routes/gtu/organic/thermal/") and tok.endswith(".py"):
                assert (ROOT / tok).exists(), tok
    r = subprocess.run([sys.executable, str(ROOT / "routes/gtu/platform/aws/runner.py"), "plan", str(jf),
                        "--assume-price", "0.955"], capture_output=True, text=True, cwd=ROOT)
    assert r.returncode == 0, r.stderr[-2000:]
