#!/usr/bin/env python3
"""Трёхмерная оптимизация топологии на вокселях (гексаэдр 8 узлов) — для задач AWS и малых демонстраций.

Задачи:
  disk  — диск турбины ступени 1 KV-3: 1/8 модели (симметрия x = 0, y = 0, z = 0; 68 лопаток — 17 на четверть,
          узор симметричен относительно обеих плоскостей), центробежная и тепловая нагрузка, сила каждой лопатки —
          радиальная на грани обода в пределах половины шага паза; ограничения по напряжению p-нормой (как topo_axi:
          LC1 — σ_LM(T, 25 000 ч), LC2 — σв(T) при 120 %), цель — масса; цель по загрузке — из осесимметричного
          расчёта исходного профиля (disk_topo_results.json → baseline.turb1), иначе 1.
  frame — опорная рама оболочки (shell_topo.py), объёмная: рама длиной Lz, подшипник по средней части длины, лапы —
          по средней половине; min-max податливости по 8 направлениям силы подшипника и давлению, доли металла по зонам.
Решатель — PARDISO (многопоточный, pypardiso), запасной — SuperLU/AMG (fe3d.solve_spd).

Запуск: python3 topo3d.py --problem disk --h 1.5e-3 --iters 160 --out DIR
        python3 topo3d.py --problem frame --n 110 --iters 150 --out DIR
        python3 topo3d.py --demo   (обе задачи на грубой сетке, проверка кода; ~5 мин на 4 ядрах)
"""
import argparse
import json
import math
import os
import sys
import time

import numpy as np
import scipy.sparse as sp
from scipy.spatial import cKDTree

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import fe3d  # noqa: E402
import topo_axi as TA  # noqa: E402

M6 = np.array([1, 1, 1, 0, 0, 0.0])


def vm6(s):
    return np.sqrt(np.maximum(0.5 * ((s[..., 0] - s[..., 1]) ** 2 + (s[..., 1] - s[..., 2]) ** 2 + (s[..., 2] - s[..., 0]) ** 2)
                              + 3 * (s[..., 3] ** 2 + s[..., 4] ** 2 + s[..., 5] ** 2), 0))


def dvm6(s, vm):
    v = np.maximum(vm, 1e-12)[..., None]
    sx, sy, sz = s[..., 0], s[..., 1], s[..., 2]
    return np.stack([2 * sx - sy - sz, 2 * sy - sx - sz, 2 * sz - sx - sy, 6 * s[..., 3], 6 * s[..., 4], 6 * s[..., 5]], -1) / (2 * v)


class Vox:
    def __init__(self, x0, y0, z0, h, active, group, nu, E_e, eth_e, rho):
        nx, ny, nz = active.shape
        self.h = h
        I, J, K = np.nonzero(active)
        self.ne = len(I)
        self.group = group[I, J, K]
        self.xc = x0 + (I + 0.5) * h; self.yc = y0 + (J + 0.5) * h; self.zc = z0 + (K + 0.5) * h
        loc = fe3d.NODE_XYZ.astype(int)
        gn = ((I[:, None] + loc[:, 0]) * (ny + 1) + (J[:, None] + loc[:, 1])) * (nz + 1) + (K[:, None] + loc[:, 2])
        used, inv = np.unique(gn, return_inverse=True)
        self.enodes = inv.reshape(-1, 8)
        self.nn = len(used)
        k_ = used % (nz + 1); j_ = (used // (nz + 1)) % (ny + 1); i_ = used // ((nz + 1) * (ny + 1))
        self.nx_ = x0 + i_ * h; self.ny_ = y0 + j_ * h; self.nz_ = z0 + k_ * h
        self.edof = np.empty((self.ne, 24), int)
        for d in range(3):
            self.edof[:, d::3] = 3 * self.enodes + d
        self.ndof = 3 * self.nn
        self.Ke, self.Bc, self.fth, self.fb, self.D1 = fe3d.hex_matrices(h, nu, 1.0)
        self.Es = E_e[I, J, K] if np.ndim(E_e) == 3 else np.full(self.ne, float(E_e))
        self.eth = eth_e[I, J, K] if np.ndim(eth_e) == 3 else np.full(self.ne, float(eth_e))
        self.rho = rho
        self.V = h ** 3
        ii, jj = np.meshgrid(np.arange(24), np.arange(24), indexing="ij")
        self.Ki = self.edof[:, ii].ravel(); self.Kj = self.edof[:, jj].ravel()
        self.fixed = np.zeros(self.ndof, bool)
        self.f_fix0 = np.zeros(self.ndof)

    def finish(self):
        free = np.nonzero(~self.fixed)[0]
        m = -np.ones(self.ndof, int); m[free] = np.arange(len(free))
        ri, ci = m[self.Ki], m[self.Kj]
        self.keep = (ri >= 0) & (ci >= 0)
        self.ri, self.ci, self.free = ri[self.keep], ci[self.keep], free
        del self.Ki, self.Kj

    def K(self, eE):
        data = (eE[:, None, None] * self.Ke[None]).reshape(-1)[self.keep]
        nf = len(self.free)
        K = sp.csr_matrix((data, (self.ri, self.ci)), shape=(nf, nf))
        return K


def filt3(m, R):
    pts = np.c_[m.xc, m.yc, m.zc]
    tr = cKDTree(pts)
    D = tr.sparse_distance_matrix(tr, R, output_type="coo_matrix")
    H = sp.coo_matrix((np.maximum(R - D.data, 0), (D.row, D.col)), shape=(m.ne, m.ne)).tocsr()
    H = H + sp.identity(m.ne, format="csr") * R
    s = np.asarray(H.sum(axis=1)).ravel()
    return sp.diags(1 / s) @ H


# =========================================================================== задача: диск турбины
def problem_disk(h):
    import disk_topo as DT
    d = DT.frame3()
    P = d["parts"]["turb1"]
    g = P["geometry"]; pi = P["profile_info"]
    r_live, w, r_b, t_max = g["r_live_rim_m"], g["w_rim_m"], pi["r_b_m"], pi["t_max_m"]
    Nb = g["N_blades"]
    Tv = P["variants"]["cooled_750"]["T_C"]
    p_rim = P["variants"]["cooled_750"]["cases"]["100"]["rim_p_MPa"] * 1e6
    mat = TA.DK.load_material("in713c"); mt = TA.mat_tables(mat)
    n = int(math.ceil(r_live / h)); nz = int(math.ceil(t_max / 2 / h))
    I, J, K = np.meshgrid(np.arange(n), np.arange(n), np.arange(nz), indexing="ij")
    xc, yc, zc = (I + 0.5) * h, (J + 0.5) * h, (K + 0.5) * h
    r = np.hypot(xc, yc)
    active = (r <= r_live) & ((r < r_b) | (zc <= w / 2))
    group = np.where(r >= r_b, 0, 1)
    T = np.interp(r, Tv["r"], Tv["T"])
    E = mt["E"](T); eth = mt["alpha"](T) * (T - 20.0)
    m = Vox(0, 0, 0, h, active, group, mt["nu"], E, eth, mt["rho"])
    m.Te = np.interp(np.hypot(m.xc, m.yc), Tv["r"], Tv["T"])
    # симметрия
    m.fixed[3 * np.nonzero(m.nx_ < 1e-12)[0]] = True
    m.fixed[3 * np.nonzero(m.ny_ < 1e-12)[0] + 1] = True
    m.fixed[3 * np.nonzero(m.nz_ < 1e-12)[0] + 2] = True
    # лопатки: 17 на четверть, углы (k + 1/2)·360/68; сила лопатки = p·2π r w / 68, половина — на половину толщины
    F_blade = p_rim * 2 * math.pi * r_live * w / Nb / 2
    nr_ = np.hypot(m.nx_, m.ny_); nth = np.degrees(np.arctan2(m.ny_, m.nx_))
    dth = 360.0 / Nb
    rim_nodes = (nr_ >= r_live - 0.75 * h) & (m.nz_ <= w / 2 + 1e-9)
    nb = 0
    for k in range(Nb):
        a = (k + 0.5) * dth
        if a > 90.0:
            break
        sel = np.nonzero(rim_nodes & (np.abs(nth - a) <= 0.25 * dth))[0]
        if len(sel) == 0:
            continue
        fr = F_blade / len(sel)
        m.f_fix0[3 * sel] += fr * m.nx_[sel] / nr_[sel]
        m.f_fix0[3 * sel + 1] += fr * m.ny_[sel] / nr_[sel]
        nb += 1
    m.finish()
    uts = mt["uts"](m.Te) * 1e6
    a1 = np.minimum(mt["slm"](m.Te, 25000.0) * 1e6, uts)
    cases = [dict(rpm=14950.0, allow=a1), dict(rpm=17940.0, allow=uts)]
    meta = dict(problem="disk", part="turb1", h=h, n=n, nz=nz, elements=int(m.ne), dofs=int(m.ndof), blades_quarter=nb,
                F_blade_half_N=F_blade, r_live=r_live, r_b=r_b, w=w, t_max=t_max)
    return m, cases, meta


def solve_disk_lc(m, xb, cases, prm, need_sens=True):
    q, Emin, qs, P = prm["q"], prm["Emin"], prm["qs"], prm["P"]
    x = xb
    den = 1 + q * (1 - x)
    eE = np.where(m.group == 0, 1.0, Emin + (1 - Emin) * x / den)
    dE = np.where(m.group == 0, 0.0, (1 - Emin) * (1 + q) / den ** 2)
    em = np.where(m.group == 0, 1.0, x); dm = np.where(m.group == 0, 0.0, 1.0)
    xs = np.maximum(x, 1e-9)
    es = np.where(m.group == 0, 1.0, xs ** qs); ds = np.where(m.group == 0, 0.0, qs * xs ** (qs - 1))
    Ee = eE * m.Es
    Kf = m.K(Ee)
    fac = fe3d.Factor(Kf)
    F = np.zeros((m.ndof, len(cases)))
    fb_x = m.fb[:, 0]; fb_y = m.fb[:, 1]
    for k, c in enumerate(cases):
        w2 = (2 * math.pi * c["rpm"] / 60) ** 2
        fe = (m.rho * w2 * em)[:, None] * (m.xc[:, None] * fb_x[None] + m.yc[:, None] * fb_y[None])
        fe += (Ee * m.eth)[:, None] * m.fth[None]
        np.add.at(F[:, k], m.edof, fe)
        F[:, k] += m.f_fix0 * (c["rpm"] / 14950.0) ** 2
    U = np.zeros_like(F)
    U[m.free] = fac.solve(F[m.free]).reshape(len(m.free), -1)
    res = []
    rhs = np.zeros_like(F)
    for k, c in enumerate(cases):
        ue = U[m.edof, k]
        eps = ue @ m.Bc.T - m.eth[:, None] * M6[None]
        sig = m.Es[:, None] * (eps @ m.D1.T)
        vm = vm6(sig)
        u = es * vm / c["allow"]
        umax = float(u.max())
        PN = float(np.sum(np.maximum(u, 1e-30) ** P) ** (1 / P))
        r = dict(u=u, umax=umax, PN=PN, vm=vm, w2=(2 * math.pi * c["rpm"] / 60) ** 2)
        if need_sens:
            dPdu = PN ** (1 - P) * np.maximum(u, 1e-30) ** (P - 1)
            r["dexpl"] = dPdu * ds * vm / c["allow"]
            dsig = (dPdu * es / c["allow"])[:, None] * dvm6(sig, vm)
            ge = (dsig * m.Es[:, None]) @ m.D1 @ m.Bc
            np.add.at(rhs[:, k], m.edof, ge)
        res.append(r)
    if need_sens:
        L = np.zeros_like(rhs)
        L[m.free] = fac.solve(rhs[m.free]).reshape(len(m.free), -1)
        for k, c in enumerate(cases):
            ue = U[m.edof, k]; le = L[m.edof, k]
            dfe = (m.rho * res[k]["w2"] * dm)[:, None] * (m.xc[:, None] * fb_x[None] + m.yc[:, None] * fb_y[None])
            dfe += (dE * m.Es * m.eth)[:, None] * m.fth[None]
            dKu = (dE * m.Es)[:, None] * (ue @ m.Ke.T)
            res[k]["dPN"] = res[k]["dexpl"] + np.sum(le * (dfe - dKu), axis=1)
    fac.free()
    return res, em


def run_disk(h, iters, out, u_target=None, log=print):
    t0 = time.time()
    m, cases, meta = problem_disk(h)
    log(json.dumps(meta))
    if u_target is None:
        p = os.path.join(HERE, "disk_topo_results.json")
        u_target = [1.0, 1.0]
        if os.path.exists(p):
            with open(p, encoding="utf-8") as f:
                b = json.load(f)["baseline"]["turb1"]
            u_target = [c["u_max"] for c in b["cases"]]
    prm = dict(q=5.0, Emin=1e-6, qs=0.5, P=10.0)
    design = m.group == 1
    idx = np.nonzero(design)[0]
    H = filt3(m, 2.5 * h)
    x = np.ones(m.ne)
    mma = TA.MMA(len(idx), 2, move=0.1)
    M0 = float(np.sum(m.rho * m.V * np.ones(m.ne)))
    cn = np.ones(2)
    sched = [(1, 0.3), (2, 0.2), (4, 0.15), (8, 0.15), (16, 0.2)]
    hist = []
    it = 0
    tg = [t * 0.97 for t in u_target]
    for beta, share in sched:
        for _ in range(max(2, int(share * iters))):
            xt = H @ x
            xb, dxb = TA.project(xt, beta); xb[~design] = 1.0; dxb[~design] = 0.0
            res, em = solve_disk_lc(m, xb, cases, prm)
            f0 = float(np.sum(m.rho * m.V * em) / M0)
            df = (H.T @ (np.where(design, m.rho * m.V / M0, 0.0) * dxb))[idx]
            g = np.zeros(2); dg = np.zeros((2, len(idx)))
            for k, r in enumerate(res):
                ck = r["umax"] / r["PN"]
                cn[k] = ck if it == 0 else 0.5 * ck + 0.5 * cn[k]
                g[k] = cn[k] * r["PN"] / tg[k] - 1
                dg[k] = (H.T @ (cn[k] / tg[k] * r["dPN"] * dxb))[idx]
            mma.move = 0.1 if beta < 8 else 0.04
            xn, _ = mma.step(x[idx], f0, df, g, dg)
            ch = float(np.abs(xn - x[idx]).max()); x[idx] = xn
            rec = dict(it=it, beta=beta, mass_rel=f0, umax=[r["umax"] for r in res], g=g.tolist(), change=ch,
                       t=round(time.time() - t0, 1))
            hist.append(rec); log(json.dumps(rec))
            it += 1
    xt = H @ x; xb, _ = TA.project(xt, sched[-1][0]); xb[~design] = 1
    xe = np.where(xb >= 0.5, 1.0, 0.0); xe[~design] = 1
    res, em = solve_disk_lc(m, xe, cases, prm, need_sens=False)
    mass_q = float(np.sum(m.rho * m.V * em))
    out_d = dict(meta, u_target=u_target, mass_eighth_kg=mass_q, mass_full_kg=8 * mass_q,
                 u_binary=[r["umax"] for r in res], iterations=it, seconds=round(time.time() - t0, 1))
    os.makedirs(out, exist_ok=True)
    np.savez_compressed(os.path.join(out, f"topo3d_disk_h{h * 1e3:.1f}.npz"), xb=xb, xe=xe, xc=m.xc, yc=m.yc, zc=m.zc, h=h,
                        hist=json.dumps(hist))
    with open(os.path.join(out, f"topo3d_disk_h{h * 1e3:.1f}.json"), "w", encoding="utf-8") as f:
        json.dump(out_d, f, ensure_ascii=False, indent=1)
    _stl(xb, m, h, os.path.join(out, f"topo3d_disk_h{h * 1e3:.1f}.stl"), mirror="disk")
    return out_d


# =========================================================================== задача: рама оболочки
def run_frame(n, iters, out, Lz=0.20, log=print):
    import shell_topo as ST
    t0 = time.time()
    fp = ST.frame3_point()
    D2 = max(fp["stage1"]["D2_m"], fp["stage2"]["D2_m"])
    R = 1.9 * D2 / 2
    h = 2 * R / n
    nz = max(4, int(round(Lz / h)))
    g = ST.GEO
    ax = -R + (np.arange(n) + 0.5) * h
    X, Y, Z = np.meshgrid(ax, ax, (np.arange(nz) + 0.5) * h, indexing="ij")
    r = np.hypot(X, Y)
    r_si = R - g["t_band"]
    active = (r <= R) & (r >= g["r_bear_i"])
    zone = np.full(X.shape, -1)
    zone[(r < g["r_bear_o"])] = 0
    zone[(r >= g["r_bear_o"]) & (r < g["r_flow_hub"])] = 1
    zone[(r >= g["r_flow_hub"]) & (r < r_si)] = 2
    zone[(r >= r_si)] = 3
    passive = (zone == 0) | (zone == 3)
    void = np.zeros(X.shape, bool)
    for ang, d, _ in ST.RADIAL_CH:
        t = math.radians(ang); ux, uy = math.cos(t), math.sin(t)
        along = X * ux + Y * uy; perp = np.hypot(-X * uy + Y * ux, Z - nz * h / 2)
        seg = (along >= g["r_bear_o"] - h) & (along <= R - g["skin"])
        void |= seg & (perp <= d / 2)
        passive |= seg & (perp > d / 2) & (perp <= d / 2 + max(g["skin"], h))
    for ang, d, _ in ST.AXIAL_CH:
        t = math.radians(ang); rc = r_si + g["t_band"] / 2
        void |= np.hypot(X - rc * math.cos(t), Y - rc * math.sin(t)) <= max(d / 2, 0.6 * h)
    active &= ~void
    group = np.where(passive, 0, 1)
    m = Vox(-R, -R, 0, h, active, group, 0.342, 93.7e9, 0.0, 4430.0)
    I, J, K = np.nonzero(active)
    m.zone = zone[I, J, K]
    nr_ = np.hypot(m.nx_, m.ny_); nth = np.degrees(np.arctan2(m.ny_, m.nx_)) % 360
    for ang, half in ST.MOUNTS:
        dth = np.abs((nth - ang + 180) % 360 - 180)
        sel = np.nonzero((nr_ >= R - 1.5 * h) & (dth <= half) & (m.nz_ >= 0.25 * nz * h) & (m.nz_ <= 0.75 * nz * h))[0]
        for d in range(3):
            m.fixed[3 * sel + d] = True
    m.finish()
    ring = np.nonzero((nr_ <= g["r_bear_i"] + 0.75 * h) & (m.nz_ >= 0.3 * nz * h) & (m.nz_ <= 0.7 * nz * h))[0]
    F = np.zeros((m.ndof, ST.N_DIR + 1))
    for k in range(ST.N_DIR):
        a = 2 * math.pi * k / ST.N_DIR
        F[3 * ring, k] = ST.F_BEAR * math.cos(a) / len(ring); F[3 * ring + 1, k] = ST.F_BEAR * math.sin(a) / len(ring)
    pg = 13.5 * 101325.0 - 101325.0
    sh = np.nonzero(np.abs(nr_ - r_si) <= 0.6 * h)[0]
    A = 2 * math.pi * r_si * nz * h / max(len(sh), 1)
    F[3 * sh, ST.N_DIR] = pg * A * m.nx_[sh] / nr_[sh]; F[3 * sh + 1, ST.N_DIR] = pg * A * m.ny_[sh] / nr_[sh]
    meta = dict(problem="frame", n=n, nz=nz, h=h, Lz=nz * h, elements=int(m.ne), dofs=int(m.ndof))
    log(json.dumps(meta))
    design = m.group == 1
    idx = np.nonzero(design)[0]
    H = filt3(m, 2.0 * h)
    zones = {1: g["inner_frac"], 2: g["blockage"]}
    nzn = {z: np.sum((m.zone == z) & design) for z in zones}
    x = np.zeros(m.ne); x[idx] = 0.15; x[~design] = 1
    nk = F.shape[1]
    mma = TA.MMA(len(idx) + 1, nk + 2, move=0.1); mma.xmax[-1] = 10
    C0 = None; zv = 1.0; hist = []; it = 0
    sched = [(1, 0.35), (2, 0.2), (4, 0.15), (8, 0.15), (16, 0.15)]
    for beta, share in sched:
        for _ in range(max(2, int(share * iters))):
            xt = H @ x; xb, dxb = TA.project(xt, beta); xb[~design] = 1; dxb[~design] = 0
            eE = 1e-6 + (1 - 1e-6) * xb ** 3
            fac = fe3d.Factor(m.K(eE))
            U = np.zeros_like(F); U[m.free] = fac.solve(F[m.free]); fac.free()
            C = np.sum(F * U, axis=0)
            if C0 is None:
                C0 = C.copy(); zv = float(np.max(C / C0)) * 1.05
            ue = U[m.edof]
            ee = np.einsum("eik,ij,ejk->ek", ue, m.Ke, ue, optimize=True)
            dC = -(3 * xb ** 2)[:, None] * ee
            gg = np.zeros(nk + 2); dg = np.zeros((nk + 2, len(idx) + 1))
            for k in range(nk):
                gg[k] = C[k] / C0[k] - zv
                dg[k, :-1] = (H.T @ (dC[:, k] / C0[k] * dxb))[idx]; dg[k, -1] = -1
            for j, (z, vmax) in enumerate(zones.items()):
                inz = (m.zone == z) & design
                gg[nk + j] = np.sum(xb[inz]) / (vmax * nzn[z]) - 1
                dg[nk + j, :-1] = (H.T @ (np.where(inz, 1 / (vmax * nzn[z]), 0) * dxb))[idx]
            mma.move = 0.1 if beta < 8 else 0.05
            xn, _ = mma.step(np.r_[x[idx], zv], zv, np.r_[np.zeros(len(idx)), 1.0], gg, dg)
            ch = float(np.abs(xn[:-1] - x[idx]).max()); x[idx] = xn[:-1]; zv = float(xn[-1])
            rec = dict(it=it, beta=beta, z=zv, Cmax_rel=float(np.max(C / C0)), change=ch, t=round(time.time() - t0, 1))
            hist.append(rec); log(json.dumps(rec)); it += 1
    xt = H @ x; xb, _ = TA.project(xt, 16); xb[~design] = 1
    out_d = dict(meta, iterations=it, Cmax_rel_final=hist[-1]["Cmax_rel"], seconds=round(time.time() - t0, 1),
                 mass_kg=float(np.sum(np.where(xb >= 0.5, 1, 0)) * h ** 3 * 4430.0))
    os.makedirs(out, exist_ok=True)
    np.savez_compressed(os.path.join(out, f"topo3d_frame_n{n}.npz"), xb=xb, xc=m.xc, yc=m.yc, zc=m.zc, h=h, hist=json.dumps(hist))
    with open(os.path.join(out, f"topo3d_frame_n{n}.json"), "w", encoding="utf-8") as f:
        json.dump(out_d, f, ensure_ascii=False, indent=1)
    _stl(xb, m, h, os.path.join(out, f"topo3d_frame_n{n}.stl"), mirror=None)
    return out_d


def _stl(xb, m, h, path, mirror=None):
    """изоповерхность x̄ = 0,5 маршем кубов (scikit-image) -> STL (trimesh); у диска — достройка зеркалами до 360°"""
    try:
        from skimage import measure
        import trimesh
    except ImportError:
        return None
    i = np.round((m.xc - m.xc.min()) / h).astype(int); j = np.round((m.yc - m.yc.min()) / h).astype(int)
    k = np.round((m.zc - m.zc.min()) / h).astype(int)
    G = np.zeros((i.max() + 3, j.max() + 3, k.max() + 3))
    G[i + 1, j + 1, k + 1] = xb
    if mirror == "disk":
        G[0] = G[1]; G[:, 0] = G[:, 1]; G[:, :, 0] = G[:, :, 1]   # плоскости симметрии — без ложной стенки
    v, f, _, _ = measure.marching_cubes(G, 0.5, spacing=(h, h, h))
    v += np.array([m.xc.min() - 1.5 * h, m.yc.min() - 1.5 * h, m.zc.min() - 1.5 * h])
    mesh = trimesh.Trimesh(v, f, process=True)
    if mirror == "disk":
        parts = []
        for sx in (1, -1):
            for sy in (1, -1):
                for sz in (1, -1):
                    mm = mesh.copy(); mm.apply_scale([sx, sy, sz]) if False else mm.apply_transform(np.diag([sx, sy, sz, 1.0]))
                    if sx * sy * sz < 0:
                        mm.invert()
                    parts.append(mm)
        mesh = trimesh.util.concatenate(parts)
        mesh.merge_vertices()
    mesh.export(path)
    return path


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--problem", choices=["disk", "frame"])
    ap.add_argument("--h", type=float, default=1.5e-3)
    ap.add_argument("--n", type=int, default=110)
    ap.add_argument("--iters", type=int, default=150)
    ap.add_argument("--out", default=os.environ.get("GTU_OUT", os.path.join(HERE, "runs", "topo3d")))
    ap.add_argument("--demo", action="store_true")
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    logf = open(os.path.join(a.out, "topo3d_log.jsonl"), "a")

    def log(s):
        print(s, flush=True); logf.write(s + "\n"); logf.flush()
    if a.demo:
        r1 = run_disk(6e-3, 40, a.out, log=log)
        r2 = run_frame(40, 40, a.out, log=log)
        res = dict(title="topo3d: демонстрация на грубой сетке (проверка кода перед AWS)", date="2026-11-14",
                   method="воксельный SIMP/RAMP, MMA, PARDISO; см. docstring", disk=r1, frame=r2)
        with open(os.path.join(HERE, "topo3d_demo_results.json"), "w", encoding="utf-8") as f:
            json.dump(res, f, ensure_ascii=False, indent=1)
        return
    if a.problem == "disk":
        run_disk(a.h, a.iters, a.out, log=log)
    else:
        run_frame(a.n, a.iters, a.out, log=log)


if __name__ == "__main__":
    main()
