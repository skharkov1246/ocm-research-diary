#!/usr/bin/env python3
"""Геометрия элит семейства «корпус-оболочка и рама»: сетки для рендера, поле SDF и проверки изготовимости по нему.

  * meshes(p, variant) — оболочка вращения с рельефом рёбер (C²-профиль со скруглением) и площадками опор, рамы
    подшипников (стойки-профили core.blade: навал, изгиб, конусность, бугорки), горячая часть и пальцы клетки,
    цапфы опор, наружные трубки для невстроенных систем (чистая форма = всё встроено);
  * heightfield(p, variant) — r_out(z, θ) на сетке (рельеф рёбер и опор);
  * shell_sdf(p, variant) — поле расстояний (стенка — по профилям, рёбра — рельеф, стойки и пальцы — капсулы);
  * check_global — наименьшая стенка, нависания (ось печати ±z), замкнутые полости по SDF на грубой сетке;
  * check_patch — участок стенки «обшивка + решётка + обшивка» с отверстиями выноса порошка на мелкой сетке
    (0,4 мм): толщина листа решётки, вынос порошка шаром 3 и 5 мм — сверка с дешёвым правилом (поры решётки).
"""
import math
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.dirname(HERE))
from core import blade as BL  # noqa: E402
from core import checks as C  # noqa: E402
from core import mesh as M  # noqa: E402
from core import shapes as S  # noqa: E402
import shell_model as SM  # noqa: E402

RES = dict(draft=(110, 220), fine=(260, 640))
GRAY = dict(shell=0.86, frame=0.70, hot=0.55, cage=0.62, mount=0.75, pipe=0.45)


def heightfield(p, variant, res="draft", geo=None):
    nz, nt = RES[res] if isinstance(res, str) else res
    geo = geo or SM.evaluate_params(p, variant, with_geom=True)
    prof = geo["prof"]
    z = np.linspace(0.0, prof["z"][-1], nz)
    th = np.linspace(0.0, 2 * math.pi, nt, endpoint=False)
    r_o = np.interp(z, prof["z"], prof["r_o"])
    r_i = np.interp(z, prof["z"], prof["r_i"])
    R_ref = geo["R_ref"]
    s = th * R_ref
    h = np.zeros((nz, nt))
    if geo["segs"]:
        reach = 0.5 * max(sg[4] for sg in geo["segs"]) + p["rib_fillet"] + 0.01
        D, W = S.segments_distance_grid(z, s, geo["segs"], reach, period_y=2 * math.pi * R_ref)
        prof_r = S.rib_profile(D, 0.5 * np.maximum(W, 1e-4), p["rib_fillet"])
        h = np.where(np.isfinite(D), p["rib_h"] * prof_r, 0.0)
    # площадки опор (бобышки цапф и тяги): R 0,07 м, высота 0,04 м, скругление 0,03 м (оценка)
    pads = []
    for (zm, sm, _, _) in geo["mounts"]:
        pads.append((zm, sm, zm, sm, 0.14))
    Dp, _ = S.segments_distance_grid(z, s, pads, 0.12, period_y=2 * math.pi * R_ref)
    hp = np.where(np.isfinite(Dp), 0.04 * S.rib_profile(Dp, 0.07, 0.03), 0.0)
    h = np.maximum(h, hp)
    # кольца жёсткости рам (передняя, задняя, клетка): ring_h × ring_w, то же скругление, что у рёбер
    for zr in (SM.INP["z_fb"], SM.INP["z_rb"], SM.INP["zones"][4]["z0"] + p["hf_span"]):
        hr = p["ring_h"] * S.rib_profile(np.abs(z - zr), 0.5 * p["ring_w"], p["rib_fillet"])
        h = np.maximum(h, hr[:, None])
    # концы: рельеф гаснет у фланцев (2 см)
    fade = S.smoothstep5(np.minimum(z, z[-1] - z) / 0.02)[:, None]
    h = h * fade
    return dict(z=z, th=th, r_o=r_o, r_i=r_i, h=h, R_ref=R_ref, geo=geo)


def _frame_struts(p, geo, front=True):
    pre = "f_" if front else "r_"
    N = p[pre + "n"]
    zc = SM.INP["z_fb"] if front else SM.INP["z_rb"]
    prof = geo["prof"]
    r_case = float(np.interp(zc, prof["z"], prof["r_i"]))
    if geo.get("variant") == "картридж":
        r_case = SM.INP["zones"][0 if front else 4]["r_req"]
    L0 = r_case - SM.R_HUB
    chord, tc = p[pre + "chord"], p[pre + "tc"]
    lean = math.tan(math.radians(p[pre + "lean"])) * L0
    cone = math.tan(math.radians(p[pre + "cone"])) * L0 * (1 if front else -1)
    bow = p[pre + "bow"] * L0
    fr = geo["frames"]["front" if front else "rear"]
    out = []
    for j in range(N):
        th = 2 * math.pi * j / N - math.pi / 2
        tcj = (fr["t_service"] / chord) if j == 0 else tc
        tub = (p["f_tub"], 0.4) if (front and p["f_tub"] > 0.01) else None
        Gp, F, _ = BL.blade(SM.R_HUB, r_case + 0.004, chord, 0.0, tcj, camber=(0.0, 0.0, 0.0), sweep=(0.0, cone), lean=(0.0, lean),
                            bow=bow, tubercle=tub, n_span=17, n_chord=14, x_stack=0.5)
        P = Gp.reshape(-1, 3)
        er = np.array([math.cos(th), math.sin(th), 0.0]); et = np.array([-math.sin(th), math.cos(th), 0.0]); ez = np.array([0, 0, 1.0])
        V = (P[:, 0:1] * ez + P[:, 1:2] * et + P[:, 2:3] * er) + np.array([0, 0, zc])
        out.append((V, F))
    # ступица (корпус подшипника)
    zz = np.linspace(zc - 0.6 * chord, zc + 0.6 * chord, 6)
    Vh, Fh = M.revolve_heightfield(zz, np.linspace(0, 2 * math.pi, 48, endpoint=False), np.full((6, 48), SM.R_HUB), np.full(6, SM.R_HUB - 0.03))
    out.append((Vh, Fh))
    return M.merge(out)


def _cage(p, geo):
    """горячая часть (тело вращения, условно) и пальцы клетки (винтовые при форме > 1)"""
    prof = geo["prof"]
    z0 = SM.INP["zones"][4]["z0"]
    span = p["hf_span"]
    r_hot = SM.INP["r_hot"]
    zz = np.linspace(z0 - 0.05, SM.INP["zs"][-1] - 0.01, 12)
    Vh, Fh = M.revolve_heightfield(zz, np.linspace(0, 2 * math.pi, 64, endpoint=False), np.full((12, 64), r_hot), np.full(12, r_hot - 0.02))
    parts = [(Vh, Fh)]
    N = p["hf_n"]
    r_c = float(np.interp(z0 + span, prof["z"], prof["r_i"]))
    sweep = math.sqrt(max(p["hf_shape"] ** 2 - 1.0, 0.0)) * span / max(r_hot, 1e-3)
    tt = np.linspace(0, 1, 12)
    for j in range(N):
        th0 = 2 * math.pi * j / N
        th = th0 + sweep * tt
        r = r_hot + (r_c - r_hot) * S.smoothstep5(tt)
        Pz = z0 + span * tt
        P = np.stack([r * np.cos(th), r * np.sin(th), Pz], 1)
        parts.append(M.tube(P, max(0.5 * p["hf_t"], 0.003), n_sides=5, caps=False))
    return M.merge(parts)


def _mount_pins(p, geo, hf):
    prof = geo["prof"]
    out = []
    for (zm, sm, _, _), th in zip(geo["mounts"], geo["minfo"]["th"]):
        r0 = float(np.interp(zm, prof["z"], prof["r_o"])) + 0.03
        e = np.array([math.cos(th), math.sin(th), 0.0])
        P = np.array([e * r0 + [0, 0, zm], e * (r0 + 0.12) + [0, 0, zm]])
        out.append(M.tube(P, 0.045, n_sides=16))
    return M.merge(out)


def _pipes(p, geo):
    """невстроенные системы — трубки снаружи корпуса (условно по трассе)"""
    prof = geo["prof"]
    out = []
    emb = geo["emb"]
    if "топливо" not in emb:
        zr = geo["fuel"]["info"]["z_inj"]
        r = float(np.interp(zr, prof["z"], prof["r_o"])) + 0.05
        th = np.linspace(0, 2 * math.pi, 73)
        out.append(M.tube(np.stack([r * np.cos(th), r * np.sin(th), np.full(73, zr)], 1), 0.010, n_sides=8, caps=False))
        for k in range(SM.INP["fuel"]["n_inj"]):
            a = 2 * math.pi * (k + 0.5) / SM.INP["fuel"]["n_inj"]
            out.append(M.tube(np.array([[r * math.cos(a), r * math.sin(a), zr], [(r - 0.05) * math.cos(a), (r - 0.05) * math.sin(a), zr]]),
                              0.006, n_sides=6))
    k = 0
    for nm, d, Lc, _ in SM.CHANNELS[1:]:
        if nm in emb:
            continue
        for q in range(2):
            a = -math.pi / 2 + (0.35 + 0.18 * k + 0.08 * q)
            zz = np.linspace(0.05 + 0.2 * k, min(prof["z"][-1] - 0.05, 0.45 + 0.25 * k), 20)
            r = np.interp(zz, prof["z"], prof["r_o"]) + 0.045 + 0.02 * q
            out.append(M.tube(np.stack([r * math.cos(a), r * math.sin(a), zz], 1), 0.5 * d + 0.002, n_sides=8))
        k += 1
    return M.merge(out) if out else None


def meshes(p, variant, res="draft"):
    geo = SM.evaluate_params(p, variant, with_geom=True)
    geo["variant"] = variant
    hf = heightfield(p, variant, res, geo)
    Vs, Fs = M.revolve_heightfield(hf["z"], hf["th"], hf["r_o"][:, None] + hf["h"], hf["r_i"])
    out = [(Vs, Fs, GRAY["shell"])]
    V, F = _frame_struts(p, geo, True); out.append((V, F, GRAY["frame"]))
    V, F = _frame_struts(p, geo, False); out.append((V, F, GRAY["frame"]))
    V, F = _cage(p, geo); out.append((V, F, GRAY["cage"]))
    V, F = _mount_pins(p, geo, hf); out.append((V, F, GRAY["mount"]))
    pp = _pipes(p, geo)
    if pp is not None:
        out.append((pp[0], pp[1], GRAY["pipe"]))
    # ось двигателя z → горизонталь на виде сбоку; модель: y — вверх (цапфы на горизонтали ±x)
    return out


# =========================================================================== SDF и проверки
def shell_sdf(p, variant, hf=None):
    hf = hf or heightfield(p, variant, (220, 480))
    geo = hf["geo"]
    z_ax, th_ax, R_ref = hf["z"], hf["th"], hf["R_ref"]
    nt = len(th_ax)
    caps = []
    for front in (True, False):
        pre = "f_" if front else "r_"
        zc = SM.INP["z_fb"] if front else SM.INP["z_rb"]
        r_case = float(np.interp(zc, geo["prof"]["z"], geo["prof"]["r_i"]))
        L0 = r_case - SM.R_HUB
        lean = math.tan(math.radians(p[pre + "lean"])) * L0
        cone = math.tan(math.radians(p[pre + "cone"])) * L0 * (1 if front else -1)
        t = p[pre + "chord"] * p[pre + "tc"]
        for j in range(p[pre + "n"]):
            th = 2 * math.pi * j / p[pre + "n"] - math.pi / 2
            er = np.array([math.cos(th), math.sin(th), 0.0]); et = np.array([-math.sin(th), math.cos(th), 0.0])
            for q in range(4):
                s0, s1 = q / 4, (q + 1) / 4
                bw0 = p[pre + "bow"] * L0 * 4 * s0 * (1 - s0); bw1 = p[pre + "bow"] * L0 * 4 * s1 * (1 - s1)
                a = er * (SM.R_HUB + s0 * L0) + et * (lean * s0 + bw0) + np.array([0, 0, zc + cone * s0])
                b = er * (SM.R_HUB + s1 * L0) + et * (lean * s1 + bw1) + np.array([0, 0, zc + cone * s1])
                caps.append((a, b, 0.5 * t))

    def sdf(P):
        x, y, z = P[..., 0], P[..., 1], P[..., 2]
        r = np.hypot(x, y)
        th = np.mod(np.arctan2(y, x), 2 * math.pi)
        ri = np.interp(z, geo["prof"]["z"], geo["prof"]["r_i"])
        ro = np.interp(z, geo["prof"]["z"], geo["prof"]["r_o"])
        iz = np.clip(np.interp(z, z_ax, np.arange(len(z_ax))), 0, len(z_ax) - 1)
        it = th / (2 * math.pi) * nt
        i0 = np.floor(iz).astype(int); i1 = np.minimum(i0 + 1, len(z_ax) - 1); fz = iz - i0
        j0 = np.floor(it).astype(int) % nt; j1 = (j0 + 1) % nt; ft = it - np.floor(it)
        H = hf["h"]
        h = (1 - fz) * ((1 - ft) * H[i0, j0] + ft * H[i0, j1]) + fz * ((1 - ft) * H[i1, j0] + ft * H[i1, j1])
        d_wall = np.maximum(ri - r, r - (ro + h))
        d_end = np.maximum(-z, z - z_ax[-1])
        d = np.maximum(d_wall, d_end)
        for (a, b, rad) in caps:
            d = np.minimum(d, S.sd_capsule(P, a, b, rad))
        return d
    return sdf, hf


def check_global(p, variant, h=0.006):
    sdf, hf = shell_sdf(p, variant)
    R = float(np.max(hf["r_o"]) + p["rib_h"] + 0.08)
    L = float(hf["z"][-1])
    phi, ax = C.voxelize(sdf, (-R, -R, -0.04), (R, R, L + 0.04), h)
    mw = C.min_wall(phi, h, t_min=2 * h)
    out = dict(h_m=h, grid=list(phi.shape), min_wall=mw)
    for bd, vec in (("+z", (0, 0, 1)), ("-z", (0, 0, -1))):
        out[f"overhang_{bd}"] = C.overhang(phi, h, build=vec)
    pe = C.powder_escape(phi, h, d_list=(0.012,))
    out["closed_cavities"] = dict(n=pe["n_closed"], volume_m3=pe["closed_volume_m3"])
    return out


def check_patch(p, variant, h=4e-4, size=0.05):
    """участок стенки элиты (место наибольшей толщины решётки) — см. patch_escape"""
    prof = SM.profiles(p, variant)
    k = int(np.argmax(prof["t_core"]))
    t = float(prof["t"][k]); ts = min(p["skin"], 0.5 * t); rho = float(prof["rho"][k])
    tc = t - 2 * ts
    if tc <= 2 * h or p["route"] != "печать":
        return dict(note="решётки нет")
    r = patch_escape(p["lat_type"], p["cell"], rho, tc, ts, size, h)
    return r


def pore_rule(kind, cell, rho, t_core):
    """дешёвое правило модели (shell_model.lattice_pore)"""
    return SM.lattice_pore(kind, cell, rho, t_core)


def patch_escape(kind, L, rho, tc, ts, size=0.05, h=4e-4):
    """участок стенки: обшивки ts, между ними решётка толщины tc (тип, ячейка L, плотность ρ), 2 отверстия выноса порошка
    Ø max(5 мм, L/2) во внутренней обшивке; боковые грани закрыты, снаружи от внутренней обшивки — пустота, связанная с краем.
    Возврат: наименьший лист решётки, доля пустоты решётки, куда шар 3 и 5 мм не проходит от отверстий, правило модели."""
    t = tc + 2 * ts
    d_hole = max(0.005, 0.5 * L)
    if kind == "voronoi":
        rng = np.random.Generator(np.random.PCG64(7))
        seeds = S.jittered_seeds((-L, -L, -L), (size + L, size + L, tc + L), L, rng)
        tv = 0.12 * L * (rho / 0.2) ** 0.5
    tl = float(S.tpms_t_for_density("gyroid" if kind == "voronoi" else kind, rho))

    def sdf(P):
        x, y, zz = P[..., 0], P[..., 1], P[..., 2]
        lat = S.voronoi_sdf(P, seeds, tv, "open") if kind == "voronoi" else S.tpms_sheet_sdf(P, kind, L, tl * L)
        core_box = np.maximum(-zz + ts, zz - (ts + tc))
        d_lat = np.maximum(lat, core_box)
        skin_in = np.maximum(-zz, zz - ts)
        skin_out = np.maximum(-(zz - ts - tc), zz - t)
        d_side = np.minimum(np.minimum(x, size - x), np.minimum(y, size - y)) - 0.002
        frame = np.maximum(d_side, np.maximum(-zz, zz - t))
        d = np.minimum(np.minimum(np.minimum(d_lat, skin_in), skin_out), frame)
        for (cx, cy) in ((0.25 * size, 0.25 * size), (0.75 * size, 0.75 * size)):
            d = np.maximum(d, -S.sd_cylinder_z(P - np.array([cx, cy, 0.0]), 0.5 * d_hole, -0.01, ts + 1e-4))
        return d
    phi, ax = C.voxelize(sdf, (0.0, 0.0, -0.004), (size, size, t + 0.0004), h)
    phi[:, :, -1] = -1.0
    phi[[0, -1], :, :] = -1.0; phi[:, [0, -1], :] = -1.0
    zc = ax[2]
    band = ((zc > ts) & (zc < ts + tc))[None, None, :]
    in_core = (phi > 0) & band
    mw = C.min_wall(np.where(((zc > ts + h) & (zc < ts + tc - h))[None, None, :], phi, 1.0), h, t_min=SM.PRINT_SKIN_MIN)
    pe = C.powder_escape(phi, h, d_list=(SM.PORE_MIN, SM.PORE_LONG))
    vc = float(in_core.sum()) * h ** 3
    # доля пустоты решётки, куда шар не проходит: пересчёт blocked_m3 только по слою решётки
    rel = {k2: (min(1.0, v / vc) if vc > 0 else None) for k2, v in pe["blocked_m3"].items()}
    pr = pore_rule(kind, L, rho, tc)
    return dict(kind=kind, cell_m=L, rho=rho, skin_m=ts, t_core_m=tc, h_m=h, d_hole_m=d_hole, lattice_wall_min_m=mw.get("t_min_m"),
                lattice_wall_p05_m=mw.get("t_p05_m"), core_void_m3=vc, blocked_share=rel, n_closed=pe["n_closed"],
                pore_rule_m=pr, rule_says_ok_3mm=bool(pr >= SM.PORE_MIN), rule_says_ok_5mm=bool(pr >= SM.PORE_LONG))
