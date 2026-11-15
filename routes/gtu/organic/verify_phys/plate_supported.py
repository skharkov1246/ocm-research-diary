"""Честность сравнения «пластины против гироида» (organic/thermal): у плоских каналов tpms_recup нет опор, а между
воздухом (p3 ≈ 8 бар) и газом (≈ 1,07 бар) — перепад ~7 бар. Здесь: каналы прямоугольные с рёбрами-опорами шагом w
из изгиба полосы стенки (защемлённая балка: σ = Δp·w²/(2·t²) ≤ σ_доп, σ_доп — оценка), f·Re и Nu_H1 прямоугольного
канала — своим конечно-разностным расчётом, масса и фронт — той же моделью сердцевины tpms_recup (core_eval)."""
import inspect, json, math, sys
import numpy as np
import scipy.sparse as sp, scipy.sparse.linalg as spla
from scipy.optimize import differential_evolution, minimize
sys.path.insert(0, sys.argv[1])
import tpms_recup as R

def duct(alpha, n=120):
    """прямоугольный канал b×a (alpha = b/a ≤ 1): f·Re (Фаннинг, Dh) и Nu_H1, сетка n по длинной стороне"""
    a = 1.0; b = alpha
    ny = max(8, int(round(n * alpha))); nx = n
    hx, hy = a / nx, b / ny
    N = nx * ny
    idx = lambda i, j: i * ny + j
    rows, cols, vals = [], [], []
    for i in range(nx):
        for j in range(ny):
            k = idx(i, j); d = 0.0
            for (ii, jj, h) in ((i - 1, j, hx), (i + 1, j, hx), (i, j - 1, hy), (i, j + 1, hy)):
                hh = h * h
                if 0 <= ii < nx and 0 <= jj < ny:
                    rows.append(k); cols.append(idx(ii, jj)); vals.append(1 / hh); d -= 1 / hh
                else:
                    d -= 2 / hh          # стенка на грани ячейки: значение 0 через зеркальную ячейку
            rows.append(k); cols.append(k); vals.append(d)
    L = sp.csr_matrix((vals, (rows, cols)), shape=(N, N))
    lu = spla.splu(L.tocsc())
    u = lu.solve(-np.ones(N))            # ∇²u = −1 (μ = 1, dp/dx = −1)
    um = u.mean()
    Dh = 4 * a * b / (2 * (a + b))
    fRe = Dh ** 2 / (2 * um)             # Фаннинг: f·Re = (−dp/dx)·Dh²/(2μU)
    th = lu.solve(u / um)                # ∇²θ = u/ū, θ = 0 на стенке (H1)
    Tb = np.sum(u * th) / np.sum(u)
    q_per_P = (a * b) / (2 * (a + b))    # ∫ u/ū dA / периметр
    Nu = q_per_P * Dh / (-Tb) if Tb < 0 else q_per_P * Dh / Tb
    return fRe, abs(Nu)

chk = {a: duct(a) for a in (1.0, 0.5, 0.25, 0.125)}
print("проверка решателя (Shah & London — справочные: 1 → 14,227/3,608; 0,5 → 15,548/4,123; 0,25 → 18,233/5,331; 0,125 → 20,585/6,490):")
for a, v in chk.items(): print(f"  α {a}: f·Re {v[0]:.3f}, Nu_H1 {v[1]:.3f}")

src = inspect.getsource(R.core_eval)
old = '        S = np.tile([e, b, Dh, 24.0, prof["Nu"], e, b, Dh, 24.0, prof["Nu"], b, 2 * tw / pitch], (N, 1))'
assert src.count(old) == 1
new = old + '''
        if "w" in prof:
            w = prof["w"]; py = 2 * g + 2 * tw; px = w + tw
            e = g * w / (py * px); b = 2 * (g + w) / (py * px); bmid = 2 * w / (py * px); Dh = 2 * g * w / (g + w)
            S = np.tile([e, b, Dh, prof["fRe"], prof["Nu"], e, b, Dh, prof["fRe"], prof["Nu"], bmid, 1 - 2 * e], (N, 1))'''
exec(compile(src.replace(old, new), R.__file__, "exec"), R.__dict__)

def size_rect(du, PH, PC, prof, tw):
    Qt = du["Q_kW"]
    def obj(z):
        r = R.core_eval(du, PH, PC, "plate", prof, z[0], z[1], tw=tw)
        pen = max(0.0, 1 - r["Q_kW"] / Qt) * 1e3 + max(0.0, r["dp_h"] / du["dp_h"] - 1) * 1e2 + max(0.0, r["dp_c"] / du["dp_c"] - 1) * 1e2
        return r["mass_core"] / 1e3 + 0.002 * z[0] + pen * 10
    res = differential_evolution(obj, [(1.0, 40.0), (0.03, 1.5)], seed=0, maxiter=80, popsize=20, tol=1e-8, polish=False)
    loc = minimize(obj, res.x, method="Nelder-Mead", options=dict(maxiter=400, xatol=1e-6, fatol=1e-8))
    z = loc.x if loc.fun < res.fun else res.x
    r = R.core_eval(du, PH, PC, "plate", prof, z[0], z[1], tw=tw)
    r["feasible"] = bool(r["Q_kW"] >= 0.995 * Qt and r["dp_h"] <= 1.005 * du["dp_h"] and r["dp_c"] <= 1.005 * du["dp_c"])
    return r

du_all, _ = R.duties()
du = du_all["R90-ic-RH-ibc-none-none"]
PH = R.Props("gas", du["p_h"], du["far"]); PC = R.Props("air", du["p_c"])
tw = R.PRINT["t_w"]; dP = du["p_c"] - du["p_h"]
out = dict(duty_key="R90-ic-RH-ibc-none-none", dP_bar=dP / 1e5, t_w_mm=tw * 1e3, solver_check={str(k): v for k, v in chk.items()}, rows=[])
ref = size_rect(du, PH, PC, dict(gap=0.75e-3, Nu=8.235), tw)
out["plate_no_support_Dh15"] = dict(mass_t=ref["mass_core"] / 1e3, V=ref["V"], A_fr=ref["A_fr"], feasible=ref["feasible"])
print("пластины без опор, Dh 1,5:", out["plate_no_support_Dh15"])
for sig in (50e6, 100e6):
    w = tw * math.sqrt(2 * sig / dP)
    Dh = 1.5e-3
    g = Dh * w / (2 * w - Dh)
    fRe, Nu = duct(g / w)
    r = size_rect(du, PH, PC, dict(gap=g, w=w, fRe=fRe, Nu=Nu), tw)
    row = dict(sigma_allow_MPa=sig / 1e6, w_mm=w * 1e3, gap_mm=g * 1e3, alpha=g / w, Dh_mm=Dh * 1e3, fRe=fRe, Nu_H1=Nu,
               mass_t=r["mass_core"] / 1e3, V_m3=r["V"], A_fr_m2=r["A_fr"], L_m=r["L"], wall_frac=r["wall_frac"],
               Re_h=r["Re_h"], feasible=r["feasible"])
    out["rows"].append(row); print(json.dumps(row, default=float))
json.dump(out, open(sys.argv[2], "w"), indent=1, default=float)
