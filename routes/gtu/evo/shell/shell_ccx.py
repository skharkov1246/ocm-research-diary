#!/usr/bin/env python3
"""Проверка элит семейства «корпус-оболочка и рама» прямым расчётом CalculiX (ccx 2.21) на гладкой геометрии.

Сетка — отображённая (структурированная) по (θ, z, r): стенка — 3 слоя C3D8I по толщине (обшивка / решётка или металл /
обшивка), рельеф рёбер — четвёртый слой высотой h(z, θ) из того же C²-профиля, что у рендера (усреднён по ячейке:
масса сохраняется). Ступенек вокселей нет: поверхность стенки — B-сплайн профиля, рельеф — гладкая функция.
Решётка — однородный материал с E* = E·E_rel(ρ*) (закон гироида organic/structure lattice_results.json → law).

Стойки рам и пальцы клетки — балки B31 (прямоугольник той же площади и наименьшего момента инерции, что у профиля);
корпуса подшипников и горячая часть — «пауки» из жёстких балок; ротор — жёсткая балка с двумя массами *MASS на
расстоянии ±L/√12 от центра масс (тот же момент инерции, что у стержня в дешёвой модели); горячая часть — масса.
Опоры: цапфа A — x, y, z; цапфа B — y, z; задняя тяга — y (статически определимо). Нагрузки шага 1: давление на
внутреннюю поверхность по зонам, осевые силы диафрагм и задней крышки (p·πr²), температура стенки T(z) (линейный рост
по ступеням — как в дешёвой модели), стойки и пальцы — свои температуры, рассогласование горячей части 3,39 мм
(organic/print_review) — нагревом её «паука». Шаг 2 — собственные частоты (12 мод), доли кинетической энергии ротора,
горячей части и корпуса.

Выход: наибольшее напряжение по Мизесу в стенке (без полос у опор и мест крепления балок, они — особые точки балочной
модели), зональные окружные напряжения, частоты и их отнесение, сверка с дешёвой моделью.

  python3 shell_ccx.py --merged merged_V.json --uid <uid> --variant V [--mesh coarse|fine] --out DIR
  python3 shell_ccx.py --reference "литой с фланцами" --variant оболочка --out DIR
  python3 shell_ccx.py --canonical --out DIR          (гладкий цилиндр: частота против формулы Доннелла)
"""
import argparse
import json
import math
import os
import shutil
import subprocess
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.dirname(HERE))
import shell_geom as SG  # noqa: E402
import shell_model as SM  # noqa: E402

MESH = dict(coarse=dict(nz=56, nt=96, sub=3), fine=dict(nz=200, nt=360, sub=2))
STIFF = dict(E=2.0e13, nu=0.3, rho=1.0, alpha=13e-6)


def _mat_cards(name, key, Ts, rel_E=1.0, rel_rho=1.0):
    lines = [f"*MATERIAL,NAME={name}", "*ELASTIC"]
    for T in Ts:
        mp = SM.mat_props(key, T)
        lines.append(f"{mp['E'] * rel_E:.6e},{mp['nu']:.4f},{T:.1f}")
    lines.append(f"*DENSITY\n{SM.mat_props(key, 20.0)['rho'] * rel_rho:.4f}")
    lines.append("*EXPANSION,ZERO=20.")
    for T in Ts:
        lines.append(f"{SM.mat_props(key, T)['alpha']:.6e},{T:.1f}")
    return lines


def build_model(p, variant, mesh="coarse", canonical=None):
    """canonical = dict(R, L, t, E, nu, rho): гладкий цилиндр, торцы — шарнир (u_r = u_θ = 0), без балок и нагрузок"""
    cfg = MESH[mesh]
    nz, nt, sub = cfg["nz"], cfg["nt"], cfg["sub"]
    if canonical:
        z = np.linspace(0.0, canonical["L"], nz + 1)
        r_i = np.full(nz + 1, canonical["R"] - canonical["t"] / 2)
        r_o = np.full(nz + 1, canonical["R"] + canonical["t"] / 2)
        core = np.zeros(nz + 1); ts = np.full(nz + 1, canonical["t"] / 3)
        H = np.zeros((nz + 1, nt)); geo = None
    else:
        geo = SM.evaluate_params(p, variant, with_geom=True)
        geo["variant"] = variant
        prof = geo["prof"]
        z = np.linspace(0.0, prof["z"][-1], nz + 1)
        r_i = np.interp(z, prof["z"], prof["r_i"]); r_o = np.interp(z, prof["z"], prof["r_o"])
        core = np.interp(z, prof["z"], prof["core"])
        t = r_o - r_i
        ts = np.where(core > 0.5, np.minimum(p["skin"], t / 2.5), t / 3.0)
        hf = SG.heightfield(p, variant, (sub * nz + 1, sub * nt), geo)
        Hf = hf["h"]
        # усреднение рельефа по окну ячейки (масса рёбер сохраняется на грубой сетке)
        H = np.zeros((nz + 1, nt))
        for i in range(nz + 1):
            i0, i1 = max(0, sub * i - sub // 2), min(Hf.shape[0], sub * i + sub // 2 + 1)
            blk = Hf[i0:i1]
            cols = [(sub * j + d) % Hf.shape[1] for j in range(nt) for d in range(-(sub // 2), sub // 2 + 1)]
            H[i] = blk[:, cols].reshape(blk.shape[0], nt, -1).mean(axis=(0, 2))
    th = np.linspace(0.0, 2 * math.pi, nt, endpoint=False)
    nl = 4
    # радиусы слоёв узлов: (nz+1, nt, nl+1)
    R = np.zeros((nz + 1, nt, nl + 1))
    R[..., 0] = r_i[:, None]; R[..., 1] = (r_i + ts)[:, None]; R[..., 2] = (r_o - ts)[:, None]; R[..., 3] = r_o[:, None]
    R[..., 4] = r_o[:, None] + np.maximum(H, 5e-4)
    nid = np.arange((nz + 1) * nt * (nl + 1)).reshape(nz + 1, nt, nl + 1) + 1
    X = R * np.cos(th)[None, :, None]; Y = R * np.sin(th)[None, :, None]; Z = np.broadcast_to(z[:, None, None], R.shape)
    elems, emat = [], []
    zc = 0.5 * (z[:-1] + z[1:])
    corec = np.interp(zc, z, core)
    for i in range(nz):
        for j in range(nt):
            jj = (j + 1) % nt
            for k in range(nl):
                if k == 3 and 0.25 * (H[i, j] + H[i, jj] + H[i + 1, j] + H[i + 1, jj]) < 1e-3:
                    continue
                n = [nid[i, j, k], nid[i, jj, k], nid[i + 1, jj, k], nid[i + 1, j, k],
                     nid[i, j, k + 1], nid[i, jj, k + 1], nid[i + 1, jj, k + 1], nid[i + 1, j, k + 1]]
                elems.append(n)
                emat.append("CORE" if (k == 1 and corec[i] > 0.5) else "WALL")
    E = np.array(elems)
    used = np.unique(E)
    model = dict(z=z, th=th, nid=nid, X=X, Y=Y, Z=Z, E=E, emat=emat, used=used, nz=nz, nt=nt, nl=nl, r_i=r_i, r_o=r_o, geo=geo,
                 H=H, core=core, canonical=canonical, p=p, variant=variant)
    return model


def write_inp(model, path):
    p, variant, geo = model["p"], model["variant"], model["geo"]
    z, th, nid, nz, nt = model["z"], model["th"], model["nid"], model["nz"], model["nt"]
    X, Y, Z = model["X"], model["Y"], model["Z"]
    L = []
    L.append("*HEADING\nevo/shell: проверка элиты (CalculiX)")
    L.append("*NODE, NSET=NALL")
    used = set(model["used"].tolist())
    flat = nid.ravel()
    for n, x, y, zz in zip(flat, X.ravel(), Y.ravel(), Z.ravel()):
        if n in used:
            L.append(f"{n},{x:.6e},{y:.6e},{zz:.6e}")
    nmax = int(flat.max())
    extra = []                      # (id, x, y, z, T)
    nxt = [nmax + 1]

    def new_node(x, y, zz, T):
        k = nxt[0]; nxt[0] += 1
        extra.append((k, x, y, zz, T)); return k
    patch_beams = []                # (a, b) жёсткие балки от узла крепления к узлам пятна на стенке

    def attach(iz, it, dz_n, dt_n, T):
        """узел крепления балки к стенке: пятно (2·dz_n+1)×(2·dt_n+1) узлов внутренней поверхности и столбец по толщине —
        как корень стойки со скруглением; точечное крепление к узлу объёмной сетки давало бы ложную местную податливость"""
        P = np.array([X[iz, it, 0], Y[iz, it, 0], Z[iz, it, 0]])
        P = P - 0.006 * np.array([P[0], P[1], 0.0]) / max(math.hypot(P[0], P[1]), 1e-9)   # на 6 мм внутрь: балки пятна — ненулевой длины
        tipn = new_node(*P, T)
        for di in range(-dz_n, dz_n + 1):
            for dj in range(-dt_n, dt_n + 1):
                i2, j2 = min(max(iz + di, 0), nz), (it + dj) % nt
                patch_beams.append((tipn, int(nid[i2, j2, 0])))
        for kk in (1, 2, 3):
            patch_beams.append((tipn, int(nid[iz, it, kk])))
        return tipn, P
    E = model["E"]
    ew, ec = [], []
    for k, (n, m) in enumerate(zip(E, model["emat"])):
        (ec if m == "CORE" else ew).append((k + 1, n))
    ne = len(E)
    beams = dict(STRUT_F=[], STRUT_R=[], FINGER=[], SPIDER=[], ROTOR=[])
    masses = []
    sections = {}
    temps = {}
    can = model["canonical"]
    if not can:
        shell_mat = geo["mat"]
        # ------------------------------------------------ рамы подшипников
        rotor_nodes = {}
        for front in (True, False):
            pre = "f_" if front else "r_"
            zc = SM.INP["z_fb"] if front else SM.INP["z_rb"]
            iz = int(np.argmin(np.abs(z - zc)))
            T_hub = SM.T_HUB[0] if front else SM.T_HUB[1]
            T_st = SM.INP["zones"][0]["T"] if front else SM.T_STRUT_REAR
            hub = new_node(0.0, 0.0, z[iz], T_hub)
            rotor_nodes[front] = hub
            N = p[pre + "n"]
            fr = geo["frames"]["front" if front else "rear"]
            chord, tc = p[pre + "chord"], p[pre + "tc"]
            for jj in range(N):
                tcj = (fr["t_service"] / chord) if jj == 0 else tc
                t = tcj * chord
                A = 0.685 * chord * t; I = 0.0393 * chord * t ** 3
                h_eq = math.sqrt(12 * I / A); b_eq = A / h_eq
                thj = 2 * math.pi * jj / N - math.pi / 2
                it = int(np.round((thj % (2 * math.pi)) / (2 * math.pi) * nt)) % nt
                dzc = float(z[1] - z[0])
                tip, tipxyz = attach(iz, it, max(1, int(round(0.5 * chord / dzc))), 1, T_st)
                er = np.array([math.cos(thj), math.sin(thj), 0.0]); et = np.array([-math.sin(thj), math.cos(thj), 0.0])
                L0 = float(np.hypot(tipxyz[0], tipxyz[1])) - SM.R_HUB
                lean = math.tan(math.radians(p[pre + "lean"])) * L0
                cone = math.tan(math.radians(p[pre + "cone"])) * L0 * (1 if front else -1)
                root = er * SM.R_HUB - et * lean + np.array([0.0, 0.0, z[iz] - cone])
                rn = new_node(*root, T_hub)
                beams["SPIDER"].append((hub, rn))
                prev = rn
                for q in (1, 2, 3):
                    s = q / 4
                    bw = p[pre + "bow"] * L0 * 4 * s * (1 - s)
                    P = root + s * (tipxyz - root) + et * bw
                    m_ = new_node(*P, T_st)
                    beams["STRUT_F" if front else "STRUT_R"].append((prev, m_, jj))
                    prev = m_
                beams["STRUT_F" if front else "STRUT_R"].append((prev, tip, jj))
                sections[("F" if front else "R", jj)] = (h_eq, b_eq, et)
        # ------------------------------------------------ ротор: жёсткая балка, массы ±L/√12 от центра масс
        zcg, d = SM.INP["rotor_cg"], SM.INP["rotor_L"] / math.sqrt(12.0)
        nA = new_node(0.0, 0.0, zcg - d, 20.0); nB = new_node(0.0, 0.0, zcg + d, 20.0)
        beams["ROTOR"] += [(rotor_nodes[True], nA), (nA, nB), (nB, rotor_nodes[False])]
        masses += [(nA, 0.5 * SM.INP["m_rotor"]), (nB, 0.5 * SM.INP["m_rotor"])]
        # ------------------------------------------------ горячая часть на клетке
        z0 = SM.INP["zones"][4]["z0"]
        span = p["hf_span"]
        r_hot = SM.INP["r_hot"]
        iz1 = int(np.argmin(np.abs(z - (z0 + span))))
        T_case = float(SM.T_wall(z[iz1]))
        mp_c = SM.mat_props(shell_mat, T_case)
        # «паук» горячей части нагрет так, чтобы его кольцо ушло от корпуса на dr_hot (рассогласование print_review)
        T_hs = 20.0 + (SM.INP["dr_hot"] + mp_c["alpha"] * (T_case - 20.0) * r_hot) / (STIFF["alpha"] * r_hot)
        nh = new_node(0.0, 0.0, 0.5 * (z0 + z[-1]), T_hs)
        masses.append((nh, SM.INP["m_hot"]))
        N = p["hf_n"]
        sweep = math.sqrt(max(p["hf_shape"] ** 2 - 1.0, 0.0)) * span / r_hot
        for jj in range(N):
            th0 = 2 * math.pi * jj / N
            th1 = th0 + sweep
            it = int(np.round((th1 % (2 * math.pi)) / (2 * math.pi) * nt)) % nt
            rn = new_node(r_hot * math.cos(th0), r_hot * math.sin(th0), z0, T_hs)
            beams["SPIDER"].append((nh, rn))
            prev = rn
            tip, tipxyz = attach(iz1, it, 1, 0, T_case)
            rt = float(np.hypot(tipxyz[0], tipxyz[1]))
            for q in (1, 2, 3):
                s = q / 4
                a = th0 + sweep * s
                rr = r_hot + (rt - r_hot) * s
                m_ = new_node(rr * math.cos(a), rr * math.sin(a), z0 + span * s, 600.0)
                beams["FINGER"].append((prev, m_, jj)); prev = m_
            beams["FINGER"].append((prev, tip, jj))
        # ------------------------------------------------ температура стенки
        Tz = SM.T_wall(z)
        for i in range(nz + 1):
            for n in nid[i].ravel():
                temps[int(n)] = float(Tz[i])
    for (k, x, y, zz, T) in extra:
        L.append(f"{k},{x:.6e},{y:.6e},{zz:.6e}")
        temps[k] = T
    # ------------------------------------------------ элементы
    if ew:
        L.append("*ELEMENT, TYPE=C3D8I, ELSET=EWALL")
        L += [f"{k}," + ",".join(str(int(v)) for v in n) for k, n in ew]
    if ec:
        L.append("*ELEMENT, TYPE=C3D8I, ELSET=ECORE")
        L += [f"{k}," + ",".join(str(int(v)) for v in n) for k, n in ec]
    eid = ne + 1
    beam_sets = {}
    for name, lst in beams.items():
        for item in lst:
            a, b = item[0], item[1]
            key = name if name in ("SPIDER", "ROTOR") else f"{name}_{item[2]}"
            beam_sets.setdefault(key, []).append((eid, a, b))
            eid += 1
    for key, lst in beam_sets.items():
        L.append(f"*ELEMENT, TYPE=B31, ELSET=B_{key}")
        L += [f"{e},{a},{b}" for e, a, b in lst]
    patch_ids = []
    xyz = {int(n): (x, y, zz) for n, x, y, zz in zip(flat, X.ravel(), Y.ravel(), Z.ravel())}
    xyz.update({k: (x, y, zz) for (k, x, y, zz, T) in extra})
    for (a, b) in patch_beams:
        L.append(f"*ELEMENT, TYPE=B31, ELSET=P{eid}\n{eid},{a},{b}")
        patch_ids.append((eid, a, b)); eid += 1
    if masses:
        for k, (n, m) in enumerate(masses):
            L.append(f"*ELEMENT, TYPE=MASS, ELSET=M{k}\n{eid},{n}")
            eid += 1
    # ------------------------------------------------ материалы и сечения
    Ts = [20.0, 200.0, 400.0, 600.0]
    if can:
        L.append(f"*MATERIAL,NAME=WALL\n*ELASTIC\n{can['E']:.6e},{can['nu']}\n*DENSITY\n{can['rho']}\n*EXPANSION,ZERO=20.\n1e-5")
    else:
        L += _mat_cards("WALL", shell_mat, Ts)
        Er, _ = SM.lattice_rel(float(np.mean(geo["prof"]["rho"])))
        L += _mat_cards("CORE", shell_mat, Ts, Er, float(np.mean(geo["prof"]["rho"])))
        L += _mat_cards("STRUTF", "ti64" if p["route"] != "литьё" else shell_mat, Ts)
        L += _mat_cards("STRUTR", "in718p", Ts + [800.0])
        L += _mat_cards("FINGER", "in718p", Ts + [800.0])
        L.append(f"*MATERIAL,NAME=STIFF\n*ELASTIC\n{STIFF['E']:.4e},{STIFF['nu']}\n*DENSITY\n{STIFF['rho']}\n*EXPANSION,ZERO=20.\n{STIFF['alpha']:.4e}")
    if ew:
        L.append("*SOLID SECTION, ELSET=EWALL, MATERIAL=WALL")
    if ec:
        L.append("*SOLID SECTION, ELSET=ECORE, MATERIAL=CORE")
    for key, lst in beam_sets.items():
        if key == "SPIDER":
            L.append(f"*BEAM SECTION, ELSET=B_{key}, MATERIAL=STIFF, SECTION=RECT\n0.06,0.06\n0.,0.,1.")
        elif key == "ROTOR":
            L.append(f"*BEAM SECTION, ELSET=B_{key}, MATERIAL=STIFF, SECTION=RECT\n0.20,0.20\n1.,0.,0.")
        elif key.startswith("STRUT"):
            fr = "F" if key.startswith("STRUT_F") else "R"
            jj = int(key.split("_")[-1])
            h_eq, b_eq, et = sections[(fr, jj)]
            L.append(f"*BEAM SECTION, ELSET=B_{key}, MATERIAL={'STRUTF' if fr == 'F' else 'STRUTR'}, SECTION=RECT\n{h_eq:.5e},{b_eq:.5e}\n"
                     f"{et[0]:.5f},{et[1]:.5f},{et[2]:.5f}")
        elif key.startswith("FINGER"):
            jj = int(key.split("_")[-1])
            a = 2 * math.pi * jj / p["hf_n"]
            L.append(f"*BEAM SECTION, ELSET=B_{key}, MATERIAL=FINGER, SECTION=RECT\n{p['hf_t']:.5e},{p['hf_b']:.5e}\n"
                     f"{math.cos(a):.5f},{math.sin(a):.5f},0.")
    for (e_, a, b) in patch_ids:
        d = np.array(xyz[b]) - np.array(xyz[a]); d /= (np.linalg.norm(d) + 1e-30)
        ref = np.array([0.0, 0.0, 1.0]) if abs(d[2]) < 0.9 else np.array([1.0, 0.0, 0.0])
        n1 = np.cross(d, ref); n1 /= np.linalg.norm(n1)
        L.append(f"*BEAM SECTION, ELSET=P{e_}, MATERIAL=STIFF, SECTION=RECT\n0.012,0.012\n{n1[0]:.5f},{n1[1]:.5f},{n1[2]:.5f}")
    for k, (n, m) in enumerate(masses):
        L.append(f"*MASS, ELSET=M{k}\n{m:.4f}")
    # ------------------------------------------------ опоры
    if can:
        ends = [int(n) for n in np.concatenate([nid[0].ravel(), nid[-1].ravel()]) if n in used]
        L.append("*NSET, NSET=NENDS"); L += _chunks(ends)
        L.append("*TRANSFORM, NSET=NENDS, TYPE=C\n0.,0.,0.,0.,0.,1.")
        L.append("*BOUNDARY\nNENDS,1,2")
        mid = [int(v) for v in nid[nz // 2, :, 0]]        # кольцо посередине: осевое смещение (жёсткое смещение вдоль оси)
        L.append("*NSET, NSET=NMID"); L += _chunks(mid)
        L.append("*BOUNDARY\nNMID,3,3")
        pads = {}
    else:
        pads = {}
        for name, (zm, sm, _, _), thm in zip(("MA", "MB", "MR"), geo["mounts"], geo["minfo"]["th"]):
            sel = []
            for i in range(nz + 1):
                if abs(z[i] - zm) > 0.06:
                    continue
                for j in range(nt):
                    ds = ((th[j] - thm + math.pi) % (2 * math.pi) - math.pi) * model["r_o"][i]
                    if math.hypot(z[i] - zm, ds) <= 0.06:
                        sel.append(int(nid[i, j, 3]))
            pads[name] = sel
            L.append(f"*NSET, NSET=N{name}"); L += _chunks(sel)
        L.append("*BOUNDARY\nNMA,1,3\nNMB,2,3\nNMR,2,2")
    if not can:
        L.append("*INITIAL CONDITIONS, TYPE=TEMPERATURE\nNALL,20.")
    if extra:
        L.append("*INITIAL CONDITIONS, TYPE=TEMPERATURE")
        L += [f"{k},20." for (k, *_r) in extra]
    # ------------------------------------------------ шаг 1: давление, осевые силы, температура
    loads_rings = []
    if not can:
        L.append("*STEP\n*STATIC")
        L.append("*TEMPERATURE")
        L += [f"{n},{T:.2f}" for n, T in sorted(temps.items()) if (n in used or n > nmax)]
        L.append("*DLOAD")
        zones = SM.INP["zones"]
        for k, n in enumerate(model["E"]):
            # внутренняя грань — S1 у элементов слоя 0
            if (n[0] - 1) % (model["nl"] + 1) != 0:
                continue
            zc_ = 0.25 * (Z.ravel()[n[0] - 1] + Z.ravel()[n[1] - 1] + Z.ravel()[n[2] - 1] + Z.ravel()[n[3] - 1])
            zi = int(SM._zone_of(zc_))
            pg = zones[zi]["p_g"]
            if pg > 0:
                L.append(f"{k + 1},P1,{pg:.4e}")
        L.append("*CLOAD")
        bnds = [(zones[1]["z0"], zones[0]["p_g"], zones[1]["p_g"]), (zones[2]["z0"], zones[1]["p_g"], zones[2]["p_g"]),
                (zones[3]["z0"], zones[2]["p_g"], zones[3]["p_g"]), (zones[4]["z0"], zones[3]["p_g"], zones[4]["p_g"]),
                (z[-1], zones[4]["p_g"], 0.0)]
        for zb, pl, pr in bnds:
            F = (pl - pr) * math.pi * float(np.interp(zb, z, model["r_i"])) ** 2
            if abs(F) < 1.0:
                continue
            i = int(np.argmin(np.abs(z - zb)))
            ring = [int(nid[i, j, 0]) for j in range(nt)]
            loads_rings.append(i)
            for n in ring:
                L.append(f"{n},3,{F / nt:.5e}")
        L.append("*NODE FILE\nU\n*EL FILE\nS\n*END STEP")
    # ------------------------------------------------ шаг 2: частоты
    L.append("*STEP\n*FREQUENCY\n12\n*NODE FILE\nU\n*END STEP")
    with open(path, "w") as f:
        f.write("\n".join(L) + "\n")
    return dict(nmax=nmax, extra=len(extra), elements=int(eid - 1), solids=ne, beams=sum(len(v) for v in beam_sets.values()) + len(patch_ids),
                masses=masses, pads=pads, load_rings=loads_rings, mass_nodes=[n for n, _ in masses])


def _chunks(ids, k=12):
    return [",".join(str(i) for i in ids[s:s + k]) for s in range(0, len(ids), k)]


def parse_frd_blocks(path):
    """все блоки .frd по порядку: [(имя, {узел: [значения]})]"""
    out, cur, name = [], None, None
    with open(path, encoding="utf-8", errors="replace") as f:
        for line in f:
            if line.startswith(" -4"):
                name = line.split()[1]; cur = {}
            elif line.startswith(" -3"):
                if cur is not None:
                    out.append((name, cur))
                cur = None
            elif cur is not None and line.startswith(" -1"):
                n = int(line[3:13])
                s = line.rstrip("\n")
                vals, pos = [], 13
                while pos + 12 <= len(s):
                    vals.append(float(s[pos:pos + 12])); pos += 12
                cur[n] = vals
    return out


def parse_freq(dat):
    fr, on = [], False
    with open(dat, errors="replace") as f:
        for line in f:
            if "E I G E N V A L U E   O U T P U T" in line:
                on = True; continue
            if on:
                s = line.split()
                if len(s) == 5 and s[0].isdigit():
                    fr.append(float(s[3]))
                elif fr and not s:
                    pass
                elif fr and s and not s[0].isdigit():
                    on = False
    return fr


def run(p, variant, out, mesh="coarse", canonical=None, tag="case", threads=2):
    os.makedirs(out, exist_ok=True)
    t0 = time.time()
    model = build_model(p, variant, mesh, canonical)
    inp = os.path.join(out, tag + ".inp")
    info = write_inp(model, inp)
    env = dict(os.environ, OMP_NUM_THREADS=str(threads), CCX_NPROC_EQUATION_SOLVER=str(threads))
    ccx = shutil.which("ccx")
    r = subprocess.run([ccx, "-i", tag], cwd=out, capture_output=True, text=True, env=env, timeout=7200)
    with open(os.path.join(out, tag + ".log"), "w") as f:
        f.write(r.stdout[-20000:] + r.stderr[-5000:])
    res = dict(tag=tag, mesh=mesh, rc=r.returncode, seconds=round(time.time() - t0, 1), nodes=int(len(model["used"]) + info["extra"]),
               solids=info["solids"], beams=info["beams"])
    if r.returncode != 0 or not os.path.exists(os.path.join(out, tag + ".frd")):
        res["error"] = (r.stdout[-1500:] + r.stderr[-500:])
        return res
    blocks = parse_frd_blocks(os.path.join(out, tag + ".frd"))
    freqs = parse_freq(os.path.join(out, tag + ".dat"))
    res["freq_Hz"] = freqs
    nid, z, nt = model["nid"], model["z"], model["nt"]
    used = set(model["used"].tolist())
    if canonical:
        return res
    stress = [b for nm, b in blocks if nm == "STRESS"]
    S_ = stress[0] if stress else {}
    shell_nodes = np.array(sorted(n for n in S_ if n in used))
    vm = np.array([_vm(S_[n]) for n in shell_nodes])
    # исключения: опоры (±0,08 м), кольца осевых сил (±2 ряда), места крепления балок (±1 ряд и ±1 столбец)
    excl = set()
    geo = model["geo"]
    for i in range(model["nz"] + 1):
        for (zm, sm, _, _), thm in zip(geo["mounts"], geo["minfo"]["th"]):
            if abs(z[i] - zm) <= 0.10:
                for j in range(nt):
                    ds = ((model["th"][j] - thm + math.pi) % (2 * math.pi) - math.pi) * model["r_o"][i]
                    if math.hypot(z[i] - zm, ds) <= 0.10:
                        excl.update(int(v) for v in nid[i, j])
    for i in info["load_rings"]:
        for di in range(-2, 3):
            if 0 <= i + di <= model["nz"]:
                excl.update(int(v) for v in nid[i + di].ravel())
    for front in (True, False):
        zc = SM.INP["z_fb"] if front else SM.INP["z_rb"]
        i = int(np.argmin(np.abs(z - zc)))
        for di in range(-2, 3):
            if 0 <= i + di <= model["nz"]:
                excl.update(int(v) for v in nid[i + di].ravel())
    z0 = SM.INP["zones"][4]["z0"]
    i = int(np.argmin(np.abs(z - (z0 + model["p"]["hf_span"]))))
    for di in range(-2, 3):
        if 0 <= i + di <= model["nz"]:
            excl.update(int(v) for v in nid[i + di].ravel())
    keep = np.array([n not in excl for n in shell_nodes])
    zn = {}
    for k_ in range(5):
        zz = SM.INP["zones"][k_]
        i0 = int(np.searchsorted(z, zz["z0"] + 0.25 * (zz["z1"] - zz["z0"])))
        i1 = int(np.searchsorted(z, zz["z1"] - 0.25 * (zz["z1"] - zz["z0"])))
        ns = [int(v) for v in nid[i0:i1 + 1, :, :4].ravel() if int(v) in S_ and int(v) not in excl]
        if ns:
            vv = np.array([_vm(S_[n]) for n in ns])
            zn[zz["name"]] = dict(vm_max_MPa=float(vv.max() / 1e6), vm_p99_MPa=float(np.percentile(vv, 99) / 1e6),
                                  vm_median_MPa=float(np.median(vv) / 1e6))
    res["shell_vm_max_MPa"] = float(vm[keep].max() / 1e6) if keep.any() else None
    res["shell_vm_p999_MPa"] = float(np.percentile(vm[keep], 99.9) / 1e6) if keep.any() else None
    res["shell_vm_max_incl_singular_MPa"] = float(vm.max() / 1e6)
    res["zones"] = zn
    # балки: узлы разложения (номера больше исходных) — наибольшее напряжение в стойках и пальцах
    bn = [n for n in S_ if n > info["nmax"] + info["extra"]]
    res["beam_expanded_vm_max_MPa"] = float(max(_vm(S_[n]) for n in bn) / 1e6) if bn else None
    # отнесение мод: доля кинетической энергии в массах ротора и горячей части (масса корпуса — по объёмам элементов)
    disp = [b for nm, b in blocks if nm == "DISP"]
    modes = disp[1:] if len(disp) > len(freqs) else disp[-len(freqs):] if freqs else []
    coords = _frd_coords(os.path.join(out, tag + ".frd"))
    z_e0 = SM.INP["zones"][4]["z0"]
    rows = []
    for f_, U in zip(freqs, modes):
        ns = np.array(list(U.keys())); u2 = np.array([float(np.sum(np.square(U[n][:3]))) for n in ns])
        im = int(np.argmax(u2)); xm = coords.get(int(ns[im]), (0, 0, 0))
        rm, zm = math.hypot(xm[0], xm[1]), xm[2]
        shell_share = float(u2[[int(n) in used for n in ns]].sum() / max(u2.sum(), 1e-300))
        if rm < 0.25:
            kind = "ротор и рама подшипника"
        elif zm > z_e0 - 0.02 and rm < float(np.interp(zm, z, model["r_i"])) - 0.003:
            kind = "горячая часть и клетка"
        else:
            kind = "корпус"
        rows.append(dict(f_Hz=f_, kind=kind, max_at_r_z=[round(rm, 3), round(zm, 3)], shell_share_u2=shell_share))
    res["modes"] = rows
    return res


def _frd_coords(path):
    out, on = {}, False
    with open(path, errors="replace") as f:
        for line in f:
            if line.startswith("    2C"):
                on = True; continue
            if on:
                if line.startswith(" -3"):
                    break
                if line.startswith(" -1"):
                    out[int(line[3:13])] = (float(line[13:25]), float(line[25:37]), float(line[37:49]))
    return out


def _vm(s):
    sxx, syy, szz, sxy, syz, szx = s[:6]
    return math.sqrt(0.5 * ((sxx - syy) ** 2 + (syy - szz) ** 2 + (szz - sxx) ** 2) + 3 * (sxy ** 2 + syz ** 2 + szx ** 2))


def _nodal_mass(model):
    """масса корпуса по узлам: объём шестигранника (8 тетраэдров от центра — приближённо: |J|·8) × плотность / 8"""
    X, Y, Z = model["X"].ravel(), model["Y"].ravel(), model["Z"].ravel()
    geo = model["geo"]
    rho = SM.mat_props(geo["mat"], 20.0)["rho"]
    out = {}
    for n, m in zip(model["E"], model["emat"]):
        P = np.stack([X[n - 1], Y[n - 1], Z[n - 1]], 1)
        a = P[1] - P[0]; b = P[3] - P[0]; c = P[4] - P[0]
        vol = abs(np.dot(np.cross(a, b), c))
        mm = vol * rho * (np.mean(geo["prof"]["rho"]) if m == "CORE" else 1.0) / 8.0
        for k in n:
            out[int(k)] = out.get(int(k), 0.0) + mm
    return out


def compare(res, geo):
    o = geo["objs"]
    zr = {r["zone"]: r for r in geo["zone_rows"]}
    rows = []
    for name, v in (res.get("zones") or {}).items():
        c = zr[name]
        rows.append(dict(zone=name, ccx_vm_p99_MPa=v["vm_p99_MPa"], ccx_vm_median_MPa=v["vm_median_MPa"],
                         cheap_membrane_MPa=c["sig_vm"] / 1e6, cheap_thermal_MPa=c["sig_th"] / 1e6,
                         cheap_peak_MPa=(c["Kt"] * c["sig_vm"] + c["sig_th"]) / 1e6))
    return dict(zones=rows, cheap_f_beam_Hz=o["f_beam_Hz"], cheap_f_shell_Hz=o["f_shell_Hz"], ccx_modes=res.get("modes"),
                cheap_u_wall=o["u_wall"], cheap_u_cage=o["u_cage"], cheap_cage_sigma_MPa=geo["cage"]["sigma"] / 1e6)


def _one(job):
    """один случай пакета: {kind: elite|reference|canonical, uid, variant, genome?, reference?}"""
    out = os.path.join(job["out"], job["tag"])
    try:
        if job["kind"] == "canonical":
            mp = SM.mat_props("ti64", 20.0)
            can = dict(R=0.6, L=0.5, t=0.012, E=mp["E"], nu=mp["nu"], rho=mp["rho"])
            res = run(None, None, out, job["mesh"], can, tag="canonical", threads=job["threads"])
            Am = can["E"] * can["t"] / (1 - can["nu"] ** 2); Dm = can["E"] * can["t"] ** 3 / (12 * (1 - can["nu"] ** 2))
            f, m_, n_ = SM.shell_freq_donnell(can["R"], can["L"], Am, Am, can["nu"] * Am, can["E"] * can["t"] / (2 * (1 + can["nu"])),
                                              Dm, Dm, can["nu"] * Dm, Dm * (1 - can["nu"]) / 2, can["rho"] * can["t"])
            r = dict(uid="canonical", ccx=res, donnell_Hz=f, ratio=(res["freq_Hz"][0] / f if res.get("freq_Hz") else None))
        else:
            p = SM.reference_params(job["variant"], job["reference"]) if job["kind"] == "reference" else SM.decode(job["genome"])
            geo = SM.evaluate_params(p, job["variant"], with_geom=True)
            res = run(p, job["variant"], out, job["mesh"], None, tag="case", threads=job["threads"])
            r = dict(uid=job["uid"], variant=job["variant"], ccx=res, compare=compare(res, geo))
    except Exception as ex:  # noqa: BLE001
        r = dict(uid=job.get("uid"), error=repr(ex)[:500])
    os.makedirs(out, exist_ok=True)
    with open(os.path.join(out, "ccx_result.json"), "w", encoding="utf-8") as f:
        json.dump(r, f, ensure_ascii=False, indent=1, default=lambda o: float(o) if isinstance(o, np.floating) else o.tolist())
    return r


def run_batch(path, part, n_parts, workers, mesh, threads, out):
    """пакет случаев (verify_elites.json): случай k относится к части k mod n_parts; случаи части — параллельно в workers процессов"""
    import multiprocessing as mp
    with open(path, encoding="utf-8") as f:
        cases = json.load(f)["cases"]
    mine = [dict(c, out=out, mesh=mesh, threads=threads, tag=c.get("tag") or f"c{k:03d}") for k, c in enumerate(cases) if k % n_parts == part]
    print(f"[ccx] часть {part}/{n_parts}: случаев {len(mine)}, процессов {workers}, сетка {mesh}, потоков ccx {threads}", flush=True)
    with mp.get_context("fork").Pool(max(1, workers)) as pool:
        res = pool.map(_one, mine, chunksize=1)
    with open(os.path.join(out, f"ccx_part_{part:03d}.json"), "w", encoding="utf-8") as f:
        json.dump(dict(part=part, n_parts=n_parts, mesh=mesh, results=res), f, ensure_ascii=False,
                  default=lambda o: float(o) if isinstance(o, np.floating) else o.tolist())
    return res


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--batch"); ap.add_argument("--part", type=int, default=0); ap.add_argument("--n-parts", type=int, default=1)
    ap.add_argument("--workers", type=int, default=1)
    ap.add_argument("--merged"); ap.add_argument("--uid"); ap.add_argument("--variant", default="оболочка")
    ap.add_argument("--reference"); ap.add_argument("--canonical", action="store_true")
    ap.add_argument("--mesh", default="coarse", choices=list(MESH)); ap.add_argument("--threads", type=int, default=2)
    ap.add_argument("--out", required=True); ap.add_argument("--json")
    a = ap.parse_args()
    if a.batch:
        run_batch(a.batch, a.part, a.n_parts, a.workers, a.mesh, a.threads, a.out)
        return
    if a.canonical:
        mp = SM.mat_props("ti64", 20.0)
        can = dict(R=0.6, L=0.5, t=0.012, E=mp["E"], nu=mp["nu"], rho=mp["rho"])
        res = run(None, None, a.out, a.mesh, can, tag="canonical", threads=a.threads)
        Am = can["E"] * can["t"] / (1 - can["nu"] ** 2); Dm = can["E"] * can["t"] ** 3 / (12 * (1 - can["nu"] ** 2))
        f, m_, n_ = SM.shell_freq_donnell(can["R"], can["L"], Am, Am, can["nu"] * Am, can["E"] * can["t"] / (2 * (1 + can["nu"])),
                                          Dm, Dm, can["nu"] * Dm, Dm * (1 - can["nu"]) / 2, can["rho"] * can["t"])
        res.update(case=can, donnell_Hz=f, donnell_mode=[m_, n_], ratio=(res["freq_Hz"][0] / f if res.get("freq_Hz") else None))
        out = res
    else:
        if a.reference:
            p = SM.reference_params(a.variant, a.reference); tag = "ref"; uid = "ref:" + a.reference
        else:
            m = json.load(open(a.merged, encoding="utf-8"))
            e = next(x for x in m["archive"]["elites"] if x["uid"] == a.uid)
            p = SM.decode(e["genome"]); tag = a.uid.replace("-", "_"); uid = a.uid
        geo = SM.evaluate_params(p, a.variant, with_geom=True)
        res = run(p, a.variant, a.out, a.mesh, None, tag=tag, threads=a.threads)
        out = dict(uid=uid, variant=a.variant, ccx=res, compare=compare(res, geo))
    js = a.json or os.path.join(a.out, "ccx_result.json")
    with open(js, "w", encoding="utf-8") as f:
        json.dump(out, f, ensure_ascii=False, indent=1, default=lambda o: float(o) if isinstance(o, np.floating) else o.tolist())
    print(json.dumps(out, ensure_ascii=False, default=str)[:2500])


if __name__ == "__main__":
    main()
