#!/usr/bin/env python3
"""Корпус-оболочка одной деталью: (1) оптимизация топологии сечения опорной рамы подшипника (2D, плоское
напряжённое состояние), (2) подсчёт исключённых деталей и соединений, (3) посадка оболочки KV-1…KV-3 в камеры
принтеров. Демонстрационная сетка; полная трёхмерная задача — topo3d.py, задача AWS (jobs_structure.json).

Сечение (README узла, раздел «Корпус-оболочка»): кольцо подшипника — стойки через полость и проточную часть —
стенка оболочки с осевыми каналами (топливный газ, отбор воздуха). Каналы масла (подвод сверху, слив снизу —
самотёком), вентиляции масляной полости и наддува уплотнений — радиальные, в стенках заданных полых стоек.
Нагрузки: радиальная сила подшипника по восьми направлениям (вращение дисбаланса, удары) и давление за
компрессором на стенку; опоры — две боковые лапы (как у пакетов Solar/Siemens — принято, оценка).
Цель — наименьшая податливость опоры подшипника по худшему направлению (min-max через переменную-границу),
ограничения — загромождение проточной части и доля металла во внутренней полости.
Сравнение — классическая рама из 8 прямых стоек с той же массой в тех же зонах, в той же модели.
"""
import json
import math
import os
import sys
import time

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import topo_axi as TA  # noqa: E402

ROOT = os.path.normpath(os.path.join(HERE, "..", "..", "..", ".."))
sys.path.insert(0, os.path.join(ROOT, "routes", "gtu", "funnel"))
import make_model as MK  # noqa: E402  (камеры принтеров с источниками — одна точка правды воронки)
import parts_model as PM  # noqa: E402  (счёт деталей воронки)

OUTDIR = os.environ.get("ORGANIC_OUT", "/tmp/claude-0/-home-user/a62745e4-d2fd-56e0-974f-12746b1a2ae0/scratchpad/organic")
RESULTS = os.path.join(HERE, "shell_results.json")
F3P = os.path.join(ROOT, "routes", "gtu", "frame3", "frame3_point_results.json")
FAM = os.path.join(ROOT, "routes", "gtu", "gtu_family_results.json")
LINEUP = os.path.join(ROOT, "routes", "gtu", "gtu_lineup_results.json")

# --------------------------------------------------------------------------- геометрия сечения KV-3 (оценка)
GEO = dict(
    R_o=None,            # наружный радиус оболочки = kD_case·D2max/2 (gtu_family: 1,9·D2)
    t_band=0.040,        # стенка с осевыми каналами, м (оценка: каналы до Ø24 мм + две обшивки)
    r_flow_hub=0.50,     # внутренний радиус проточной части в сечении рамы, м (оценка)
    r_bear_o=0.13, r_bear_i=0.09,   # кольцо подшипника, м (оценка: вал ~Ø0,12 м)
    b=0.10,              # осевая длина рамы (толщина плоской задачи), м (оценка)
    blockage=0.12,       # доля металла в проточной части (загромождение), оценка
    inner_frac=0.30,     # доля металла во внутренней полости, оценка
    skin=0.003,          # обшивки стенки и стенки каналов, м (оценка)
)
RADIAL_CH = [  # (угол, °, диаметр канала, м, назначение)
    (90.0, 0.012, "подвод масла к подшипнику (сверху)"),
    (270.0, 0.022, "слив масла самотёком (снизу)"),
    (210.0, 0.014, "вентиляция масляной полости"),
    (330.0, 0.014, "наддув уплотнений воздухом за компрессором"),
]
AXIAL_CH = [  # (угол, диаметр, назначение) — каналы в стенке оболочки
    (20.0, 0.024, "топливный газ к коллектору форсунок"), (160.0, 0.024, "топливный газ (второй подвод)"),
    (45.0, 0.020, "отбор воздуха на охлаждение диска турбины"), (135.0, 0.020, "отбор воздуха на охлаждение"),
    (240.0, 0.016, "дренаж камеры"), (300.0, 0.016, "проводка датчиков (термопары, давление)"),
]
MOUNTS = [(0.0, 6.0), (180.0, 6.0)]   # лапы: угол и полуширина, °
F_BEAR = 25e3        # Н, расчётная радиальная сила подшипника (оценка: ~250 кг ротора на опору × 10 g — демонстрация)
N_DIR = 8


def frame3_point():
    with open(F3P, encoding="utf-8") as f:
        return json.load(f)


class Plane:
    """Плоское напряжённое состояние, Q4, квадратная сетка шага h, толщина b."""

    def __init__(self, R, h, active, group, E, nu, b):
        self.h, self.b = h, b
        n = active.shape[0]
        self.n = n
        I, J = np.nonzero(active)
        self.ei, self.ej = I, J
        self.ne = len(I)
        self.group = group[I, J]
        x0 = -R
        self.rc = x0 + (I + 0.5) * h          # имена rc/zc — для density_filter (координаты центров)
        self.zc = x0 + (J + 0.5) * h
        npl = n + 1
        gn = np.stack([I * npl + J, (I + 1) * npl + J, (I + 1) * npl + J + 1, I * npl + J + 1], axis=1)
        used, inv = np.unique(gn, return_inverse=True)
        self.enodes = inv.reshape(-1, 4)
        self.nn = len(used)
        self.nx = x0 + (used // npl) * h
        self.ny = x0 + (used % npl) * h
        self.edof = np.empty((self.ne, 8), int)
        self.edof[:, 0::2] = 2 * self.enodes; self.edof[:, 1::2] = 2 * self.enodes + 1
        self.ndof = 2 * self.nn
        c = E / (1 - nu * nu)
        D = c * np.array([[1, nu, 0], [nu, 1, 0], [0, 0, (1 - nu) / 2]])
        self.D = D
        Ke = np.zeros((8, 8)); self.Bc = None
        for xi in TA.GP + (0.0,):
            for eta in (TA.GP if xi != 0.0 else (0.0,)):
                dNx = 0.25 * np.array([-(1 - eta), (1 - eta), (1 + eta), -(1 + eta)]) * 2 / h
                dNy = 0.25 * np.array([-(1 - xi), -(1 + xi), (1 + xi), (1 - xi)]) * 2 / h
                B = np.zeros((3, 8)); B[0, 0::2] = dNx; B[1, 1::2] = dNy; B[2, 0::2] = dNy; B[2, 1::2] = dNx
                if xi == 0.0:
                    self.Bc = B; continue
                Ke += B.T @ D @ B * b * (h * h / 4)
        self.Ke = Ke
        ii, jj = np.meshgrid(np.arange(8), np.arange(8), indexing="ij")
        self.Ki = self.edof[:, ii].ravel(); self.Kj = self.edof[:, jj].ravel()
        self.fixed = np.zeros(self.ndof, bool)

    def finish(self):
        free = np.nonzero(~self.fixed)[0]
        m = -np.ones(self.ndof, int); m[free] = np.arange(len(free))
        ri, ci = m[self.Ki], m[self.Kj]
        self.keep = (ri >= 0) & (ci >= 0); self.ri, self.ci, self.free = ri[self.keep], ci[self.keep], free

    def solve(self, eE, F):
        data = (eE[:, None, None] * self.Ke[None]).reshape(-1)[self.keep]
        nf = len(self.free)
        K = sp.csc_matrix((data, (self.ri, self.ci)), shape=(nf, nf))
        lu = spla.splu(K, permc_spec="MMD_AT_PLUS_A", options=dict(SymmetricMode=True))
        U = np.zeros_like(F)
        U[self.free] = lu.solve(F[self.free])
        return U

    def stress_vm(self, eE_solid, U):
        ue = U[self.edof]                                # (ne, 8, k)
        eps = np.einsum("ij,ejk->eik", self.Bc, ue)
        s = np.einsum("ij,ejk->eik", self.D, eps)
        sx, sy, txy = s[:, 0], s[:, 1], s[:, 2]
        return np.sqrt(sx * sx - sx * sy + sy * sy + 3 * txy * txy)


def build(R_o, h):
    g = dict(GEO, R_o=R_o)
    n = int(round(2 * R_o / h))
    h = 2 * R_o / n
    ax = -R_o + (np.arange(n) + 0.5) * h
    X, Y = np.meshgrid(ax, ax, indexing="ij")
    r = np.hypot(X, Y); th = np.degrees(np.arctan2(Y, X)) % 360
    active = (r <= R_o) & (r >= g["r_bear_i"])
    r_si = R_o - g["t_band"]
    zone = np.full((n, n), -1)
    zone[(r >= g["r_bear_i"]) & (r < g["r_bear_o"])] = 0                         # кольцо подшипника (пассивное)
    zone[(r >= g["r_bear_o"]) & (r < g["r_flow_hub"])] = 1                       # внутренняя полость (проект)
    zone[(r >= g["r_flow_hub"]) & (r < r_si)] = 2                                # проточная часть (проект)
    zone[(r >= r_si) & (r <= R_o)] = 3                                           # стенка оболочки
    passive = (zone == 0) | (zone == 3)
    void = np.zeros((n, n), bool)
    # радиальные каналы: пустота шириной d, стенки skin по бокам — от кольца подшипника до стенки
    for ang, d, _ in RADIAL_CH:
        t = np.radians(ang)
        ux, uy = math.cos(t), math.sin(t)
        along = X * ux + Y * uy; perp = np.abs(-X * uy + Y * ux)
        seg = (along >= g["r_bear_o"] - h) & (along <= R_o - g["skin"])
        void |= seg & (perp <= d / 2)
        passive |= seg & (perp > d / 2) & (perp <= d / 2 + g["skin"]) & (along <= r_si + h)
    for ang, d, _ in AXIAL_CH:
        t = np.radians(ang); rc = r_si + g["t_band"] / 2
        dd = np.hypot(X - rc * math.cos(t), Y - rc * math.sin(t))
        void |= dd <= d / 2
    active &= ~void
    group = np.where(passive, 0, 1)
    E = 93.7e9   # Ti-6Al-4V при 400 °C (materials.json, экстраполяция таблицы — как там помечено)
    m = Plane(R_o, h, active, group, E, 0.342, g["b"])
    m.zone = zone[m.ei, m.ej]
    m.theta = th[m.ei, m.ej]
    m.r = r[m.ei, m.ej]
    # лапы: узлы наружной окружности
    nr_ = np.hypot(m.nx, m.ny); nth = np.degrees(np.arctan2(m.ny, m.nx)) % 360
    for ang, half in MOUNTS:
        dth = np.abs((nth - ang + 180) % 360 - 180)
        sel = (nr_ >= R_o - 1.5 * h) & (dth <= half)
        m.fixed[2 * np.nonzero(sel)[0]] = True; m.fixed[2 * np.nonzero(sel)[0] + 1] = True
    m.finish()
    # нагрузки: подшипник — на узлы внутренней поверхности кольца (равномерно), давление — на внутреннюю обшивку стенки
    ring = np.nonzero(nr_ <= g["r_bear_i"] + 0.75 * h)[0]
    F = np.zeros((m.ndof, N_DIR + 1))
    for k in range(N_DIR):
        a = 2 * math.pi * k / N_DIR
        F[2 * ring, k] = F_BEAR * math.cos(a) / len(ring); F[2 * ring + 1, k] = F_BEAR * math.sin(a) / len(ring)
    p3 = None
    try:
        fp = frame3_point()
        p3 = float(fp["cycle"].get("p3_kPa") or fp["cycle"].get("p2_kPa") or 0) * 1e3
    except Exception:  # noqa: BLE001
        p3 = 0.0
    if not p3:
        p3 = 13.5 * 101325.0
    pg = p3 - 101325.0
    shell_in = np.nonzero(np.abs(nr_ - r_si) <= 0.6 * h)[0]
    L = 2 * math.pi * r_si / max(len(shell_in), 1)
    F[2 * shell_in, N_DIR] = pg * g["b"] * L * m.nx[shell_in] / nr_[shell_in]
    F[2 * shell_in + 1, N_DIR] = pg * g["b"] * L * m.ny[shell_in] / nr_[shell_in]
    return m, F, dict(g, h=h, n=n, p3_Pa=p3, p_gauge_Pa=pg, ring_nodes=int(len(ring)), shell_nodes=int(len(shell_in)))


def reference_spokes(m, vol_targets, n_sp=8, offset=22.5):
    """классическая рама: n_sp прямых стоек постоянной ширины в зонах 1 и 2; ширина — под ту же площадь, что у проекта"""
    xb = np.zeros(m.ne); xb[m.group == 0] = 1.0
    for z in (1, 2):
        inz = (m.zone == z) & (m.group == 1)
        if not inz.any():
            continue
        target = vol_targets[z]
        lo, hi = 0.0, 0.5
        for _ in range(40):
            w = 0.5 * (lo + hi)
            sel = np.zeros(m.ne, bool)
            for k in range(n_sp):
                a = math.radians(offset + 360.0 * k / n_sp)
                perp = np.abs(-m.rc * math.sin(a) + m.zc * math.cos(a)); along = m.rc * math.cos(a) + m.zc * math.sin(a)
                sel |= (perp <= w / 2) & (along > 0)
            area = np.sum(sel & inz)
            if area < target:
                lo = w
            else:
                hi = w
        xb[sel & inz] = 1.0
    return xb


def compliance_all(m, xb, F, q=3.0, Emin=1e-6):
    eE = Emin + (1 - Emin) * xb ** q
    U = m.solve(eE, F)
    C = np.sum(F * U, axis=0)
    return C, U, eE


def optimize_frame(m, F, iters=160, R_f=None, log=None):
    """SIMP (q = 3) + фильтр + проекция; min z при C_k/C_k0 ≤ z (k — направления подшипника и давление), V_z ≤ Vmax_z"""
    design = m.group == 1
    idx = np.nonzero(design)[0]
    nd = len(idx)
    H = TA.density_filter(m, R_f or 2.5 * m.h)
    zones = {1: GEO["inner_frac"], 2: GEO["blockage"]}
    nz_ = {z: np.sum((m.zone == z) & design) for z in zones}
    x = np.zeros(m.ne); x[idx] = 0.15; x[~design] = 1.0
    zvar = 1.0
    nk = F.shape[1]
    mma = TA.MMA(nd + 1, nk + len(zones), xmin=0.0, xmax=1.0, move=0.1)
    mma.xmax[-1] = 10.0; mma.xmin[-1] = 0.0
    C0 = None
    hist = []
    betas = [(1, 50), (2, 30), (4, 30), (8, 30), (16, 25)]
    it = 0
    for beta, nit in betas:
        for _ in range(nit):
            xt = H @ x
            xb, dxb = TA.project(xt, beta)
            xb[~design] = 1.0; dxb[~design] = 0.0
            q = 3.0
            eE = 1e-6 + (1 - 1e-6) * xb ** q
            U = m.solve(eE, F)
            C = np.sum(F * U, axis=0)
            if C0 is None:
                C0 = C.copy()
                zvar = float(np.max(C / C0)) * 1.05
            # dC/dxb_e = −q xb^(q−1) uᵀ k u
            ue = U[m.edof]
            ke_u = np.einsum("ij,ejk->eik", m.Ke, ue)
            ee = np.sum(ue * ke_u, axis=1)                        # (ne, nk)
            dC = -(q * (1 - 1e-6) * xb ** (q - 1))[:, None] * ee
            g = np.zeros(nk + len(zones)); dg = np.zeros((nk + len(zones), nd + 1))
            for k in range(nk):
                g[k] = C[k] / C0[k] - zvar
                dg[k, :nd] = (H.T @ (dC[:, k] / C0[k] * dxb))[idx]
                dg[k, nd] = -1.0
            for j, (z, vmax) in enumerate(zones.items()):
                inz = (m.zone == z) & design
                g[nk + j] = np.sum(xb[inz]) / (vmax * nz_[z]) - 1.0
                dv = np.where(inz, 1.0 / (vmax * nz_[z]), 0.0) * dxb
                dg[nk + j, :nd] = (H.T @ dv)[idx]
            xv = np.r_[x[idx], zvar]
            f0 = zvar; df0 = np.r_[np.zeros(nd), 1.0]
            mma.move = 0.1 if beta < 8 else 0.05
            xn, _ = mma.step(xv, f0, df0, g, dg)
            ch = float(np.max(np.abs(xn[:nd] - x[idx])))
            x[idx] = xn[:nd]; zvar = float(xn[nd])
            rec = dict(it=it, beta=beta, z=zvar, Cmax_rel=float(np.max(C / C0)), change=ch,
                       vol=[float(np.sum(xb[(m.zone == z) & design]) / nz_[z]) for z in zones])
            hist.append(rec)
            if log:
                log(rec)
            it += 1
    xt = H @ x
    xb, _ = TA.project(xt, betas[-1][0]); xb[~design] = 1.0
    return xb, hist, C0


def parts_count():
    """Ведомость воронки (funnel/parts_model.py): традиционный KV-3 (REF_BOM) против оболочки одной деталью.
    Дополнительно к правилам воронки — явные правила встроенных каналов (масло, топливо, воздух, датчики)."""
    ref = PM.count(PM.REF_BOM)
    shell_bom = dict(PM.REF_BOM)
    shell_bom.update(case=dict(method="печать", D=PM.REF_BOM["case"]["D"], pieces=1), hot_integral=True,
                     comb=dict(kind="кольцевая", method="печать", n_inj=PM.REF_BOM["comb"]["n_inj"]), conn="быстроразъём")
    shl = PM.count(shell_bom)
    two_bom = dict(shell_bom)
    two = PM.count(two_bom)
    # встроенные каналы: трубки масла (4 на ротор), коллектор и суфлёр — в стенке; подвод газа — канал в стенке
    # (гибкий рукав пакета остаётся); соединения трубок — по 2 на трубку (оценка)
    oil_ref = ref["by_unit"]["маслосистема"]["parts"]
    fuel_ref = ref["by_unit"]["топливо"]["parts"]
    emb = dict(
        oil_pipes_removed=oil_ref - 1,
        oil_rule="в оболочке: подвод, слив, вентиляция и наддув — каналы в стенках стоек (RADIAL_CH); суфлёр (покупной) "
                 "остаётся — 1 деталь; было: " + ref["by_unit"]["маслосистема"]["rule"],
        fuel_removed=fuel_ref - 1,
        fuel_rule="подвод газа — канал в стенке оболочки к встроенному коллектору; гибкий рукав пакета — 1 деталь",
        pipe_joints_removed=2 * ((oil_ref - 1) + (fuel_ref - 1)),
        pipe_joint_rule="по 2 резьбовых или фланцевых соединения на трубку (оценка)",
    )
    out = dict(
        reference=dict(parts=ref["parts"], assemblies=ref["assemblies"], fasteners=ref["fasteners"], n_eff=ref["n_eff"],
                       by_unit={k: dict(parts=v["parts"], fasteners=v["fasteners"], rule=v["rule"]) for k, v in ref["by_unit"].items()}),
        shell_one_piece=dict(parts=shl["parts"] - emb["oil_pipes_removed"] - emb["fuel_removed"], assemblies=shl["assemblies"],
                             fasteners=shl["fasteners"],
                             by_unit={k: dict(parts=v["parts"], fasteners=v["fasteners"], rule=v["rule"]) for k, v in shl["by_unit"].items()}),
        embedded_channels=emb,
        note="счёт по правилам funnel/parts_model.py (REF_BOM — традиционный KV-3: литые корпуса с разъёмом, болтовые фланцы, "
             "кольцевая камера из листа); оболочка — печать одной деталью (корпус, диффузоры, корпуса подшипников, жаровая труба и "
             "сопловые аппараты заодно: hot_integral), стык с выхлопным диффузором — быстроразъёмный хомут; плюс встроенные каналы")
    so = out["shell_one_piece"]
    out["removed"] = dict(parts=out["reference"]["parts"] - so["parts"], fasteners=out["reference"]["fasteners"] - so["fasteners"],
                          pipe_joints=emb["pipe_joints_removed"])
    return out


def printer_fit():
    """габарит оболочки KV-1…KV-3: D = 1,9·D2 (gtu_family), L — по формуле длины ротора воронки (funnel_space.py) для KV-3,
    для KV-1 и KV-2 — в пропорции L/D2 KV-3 (оценка); камеры — funnel/make_model.PRINTERS"""
    with open(FAM, encoding="utf-8") as f:
        fam = json.load(f)
    with open(LINEUP, encoding="utf-8") as f:
        lin = json.load(f)
    fp = frame3_point()
    D2s = [fp["stage1"]["D2_m"], fp["stage2"]["D2_m"]]
    Dm_t = fp["turbine"]["options"]["3"]["Dm"]
    L_rot_kv3 = sum(0.35 * d for d in D2s) + 0.15 * Dm_t * 3 + 0.5 * max(D2s)
    k_L = L_rot_kv3 / max(D2s)
    rows = []
    for fr in lin["lineup"][:3]:
        top = max(fr["members"], key=lambda m: m["P_shaft_kW"])
        slot = next(s for s in fam["slots"] if s["model"] == top["model"] and s["variant"] == top["variant"])
        D2 = slot["мягкий"]["D2"]
        if fr["name"] == "KV-3":
            D2 = max(D2s); basis = "frame3_point_results.json: D2 ступеней 1 и 2 (14 950 об/мин)"
            L = L_rot_kv3
        else:
            basis = f"gtu_family_results.json → {top['model']} ({top['variant']}), мягкий: D2"
            L = k_L * D2
        out = {}
        for kD in (1.7, 1.9, 2.1):
            D = kD * D2
            ps = MK.printers_for(D, L, False)
            out[f"kD{kD}"] = dict(D_m=round(D, 3), L_m=round(L, 3), printers=[p["name"] for p in ps])
        rows.append(dict(frame=fr["name"], model=top["model"], variant=top["variant"], D2_m=round(D2, 3), D2_basis=basis,
                         L_rule=("L_rot = Σ0,35·D2 + 0,15·Dm_t·3 + 0,5·D2max (funnel_space.py), " +
                                 ("KV-3 — по frame3" if fr["name"] == "KV-3" else f"в пропорции KV-3: L = {k_L:.3f}·D2 (оценка)")),
                         fit=out))
    prn = [dict(name=p["name"], proc=p["proc"], box_mm=p.get("box"), round_d_mm=p.get("round_d"), z_mm=p.get("z"),
                src=p["src"], status=p["status"], note=p.get("note")) for p in MK.PRINTERS]
    return dict(rows=rows, printers=prn, k_L=k_L, L_rot_kv3=L_rot_kv3,
                note="kD_case — диаметр оболочки к D2: 1,7–2,1 (funnel_space.py), 1,9 — gtu_family; оболочка без лап и бобышек")


def main():
    t0 = time.time()
    fp = frame3_point()
    D2 = max(fp["stage1"]["D2_m"], fp["stage2"]["D2_m"])
    R_o = 1.9 * D2 / 2
    h = R_o / 150
    m, F, geo = build(R_o, h)
    logf = open(os.path.join(OUTDIR, "log_shell.txt"), "w")
    xb, hist, C0 = optimize_frame(m, F, log=lambda r: (logf.write(json.dumps(r) + "\n"), logf.flush()))
    logf.close()
    xe = np.where(xb >= 0.5, 1.0, 0.0); xe[m.group == 0] = 1.0
    C_opt, U_opt, eE = compliance_all(m, xe, F)
    design = m.group == 1
    vols = {z: float(np.sum(xe[(m.zone == z) & design])) for z in (1, 2)}
    xr = reference_spokes(m, vols)
    C_ref, U_ref, eEr = compliance_all(m, xr, F)
    vm_o = m.stress_vm(None, U_opt) * np.where(xe > 0.5, 1.0, 0.0)[:, None]
    vm_r = m.stress_vm(None, U_ref) * np.where(xr > 0.5, 1.0, 0.0)[:, None]
    A = m.h * m.h
    mass = lambda xx: float(np.sum(xx) * A * GEO["b"] * 4430.0)  # noqa: E731
    # перемещение подшипника (среднее по кольцу) — жёсткость опоры
    nr_ = np.hypot(m.nx, m.ny)
    ring = np.nonzero(nr_ <= GEO["r_bear_i"] + 0.75 * m.h)[0]

    def stiff(U):
        ks = []
        for k in range(N_DIR):
            a = 2 * math.pi * k / N_DIR
            d = np.mean(U[2 * ring, k] * math.cos(a) + U[2 * ring + 1, k] * math.sin(a))
            ks.append(F_BEAR / d)
        return ks
    k_o, k_r = stiff(U_opt), stiff(U_ref)
    ys400 = 0.9 * float(np.interp(400.0, [20, 300, 427], [827, 538, 496])) * 1e6
    res = dict(
        title="Корпус-оболочка: оптимизация топологии сечения опорной рамы подшипника KV-3 (2D, демонстрационная сетка)",
        method="Плоское напряжённое состояние, Q4 (свой код), SIMP q = 3, фильтр плотности 2,5h, проекция Хевисайда β 1→16, MMA; "
               "min z при C_k/C_k0 ≤ z по 8 направлениям силы подшипника и давлению, ограничения доли металла по зонам; "
               "сравнение — 8 прямых стоек той же площади в тех же зонах (та же модель)",
        software=dict(python=sys.version.split()[0], numpy=np.__version__, scipy=__import__("scipy").__version__),
        date="2026-11-14",
        inputs=[dict(file="routes/gtu/frame3/frame3_point_results.json", field="stage1/stage2.D2_m (диаметр оболочки 1,9·D2), cycle (p3)"),
                dict(file="routes/gtu/platform/fea/materials.json", field="ti64: E при 400 °C (93,7 ГПа — экстраполяция таблицы), σ0,2(T)"),
                dict(file="routes/gtu/funnel/make_model.py", field="PRINTERS, printers_for"),
                dict(file="routes/gtu/funnel/parts_model.py", field="count, REF_BOM")],
        geometry=dict({k: v for k, v in geo.items()}, radial_channels=RADIAL_CH, axial_channels=AXIAL_CH, mounts=MOUNTS,
                      F_bear_N=F_BEAR, status="оценка: размеры рамы, сила подшипника и опоры приняты для демонстрации"),
        mesh=dict(n=geo["n"], h_m=geo["h"], elements=int(m.ne), dofs=int(m.ndof)),
        optimized=dict(mass_kg=mass(xe), C_rel=(C_opt / C0).tolist(), k_bearing_N_per_m=k_o,
                       vm_max_MPa=float(vm_o.max() / 1e6), vm_allow_MPa=ys400 / 1e6,
                       vol_frac={"inner": vols[1] / np.sum((m.zone == 1) & design), "flow": vols[2] / np.sum((m.zone == 2) & design)}),
        reference_8_spokes=dict(mass_kg=mass(xr), C_rel=(C_ref / C0).tolist(), k_bearing_N_per_m=k_r,
                                vm_max_MPa=float(vm_r.max() / 1e6)),
        stiffness_ratio_worst=float(min(k_o) / min(k_r)), compliance_ratio_worst_bearing=float(max(C_opt[:N_DIR]) / max(C_ref[:N_DIR])),
        iterations=len(hist), parts=parts_count(), printer_fit=printer_fit(), seconds=round(time.time() - t0, 1))
    np.savez_compressed(os.path.join(OUTDIR, "shell_frame.npz"), xb=xb, xe=xe, xr=xr, ei=m.ei, ej=m.ej, n=geo["n"], R=R_o,
                        h=geo["h"], vm_o=vm_o.max(axis=1), vm_r=vm_r.max(axis=1), hist=json.dumps(hist))
    with open(RESULTS, "w", encoding="utf-8") as f:
        json.dump(res, f, ensure_ascii=False, indent=1, default=float)
    print(json.dumps({k: res[k] for k in ("optimized", "reference_8_spokes", "stiffness_ratio_worst")}, default=float)[:1500])
    print("->", RESULTS)


if __name__ == "__main__":
    main()
