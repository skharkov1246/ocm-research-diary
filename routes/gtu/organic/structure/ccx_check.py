#!/usr/bin/env python3
"""Проверка проектов в CalculiX: меридиональное сечение (многоугольники с отверстиями) -> сетка Gmsh CAX6 -> ccx.

Материал, таблицы E(T), α(T), длительная прочность и разбор .frd/.dat — из routes/gtu/platform/fea/disk.py
(load_material, parse_frd, parse_dat_stress, von_mises): та же точка правды, что у осесимметричного узла.
Нагрузки: *DLOAD CENTRIF (ω²), *TEMPERATURE по узлам (T(r)), силы на ободе и от массы лопаток — *CLOAD,
с соглашением ccx для осесимметричных элементов, проверенным cload_convention() (сила на сектор 2° или на
радиан — определяется сравнением с давлением *DLOAD на тех же гранях). Свойства решётки — по группам
относительной плотности (E·E*/E, ρ·ρ*), напряжение решётки сравнивается с её прочностью S(ρ*).
"""
import math
import os
import subprocess
import time

import numpy as np

import topo_axi as TA

DK = TA.DK
HERE = os.path.dirname(os.path.abspath(__file__))
_CONV = None


def mesh_polygons(polys, h, h_min=None):
    """polys: [outer (n,2), hole1, ...] в (r, z), м. Возврат: nodes (n,2), tri6 (m,6) — узлы CAX6 против часовой."""
    import gmsh
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    gmsh.model.add("sec")

    def loop(poly):
        pts = [gmsh.model.geo.addPoint(float(x), float(y), 0.0, h) for x, y in poly]
        n = len(pts)
        ls = [gmsh.model.geo.addLine(pts[i], pts[(i + 1) % n]) for i in range(n)]
        return gmsh.model.geo.addCurveLoop(ls)

    loops = [loop(p) for p in polys]
    gmsh.model.geo.addPlaneSurface(loops)
    gmsh.model.geo.synchronize()
    gmsh.option.setNumber("Mesh.CharacteristicLengthMax", h)
    if h_min:
        gmsh.option.setNumber("Mesh.CharacteristicLengthMin", h_min)
    gmsh.option.setNumber("Mesh.Algorithm", 6)
    gmsh.model.mesh.generate(2)
    gmsh.model.mesh.setOrder(2)
    tags, xyz, _ = gmsh.model.mesh.getNodes()
    xyz = xyz.reshape(-1, 3)
    et, etags, en = gmsh.model.mesh.getElements(2)
    k = list(et).index(9)
    conn = np.asarray(en[k], int).reshape(-1, 6)
    gmsh.finalize()
    order = np.argsort(tags)
    tags = np.asarray(tags, int)[order]; xyz = xyz[order]
    idx = {t: i for i, t in enumerate(tags)}
    conn = np.vectorize(idx.get)(conn)
    nodes = xyz[:, :2].copy()
    # ориентация против часовой
    a, b, c = nodes[conn[:, 0]], nodes[conn[:, 1]], nodes[conn[:, 2]]
    area = 0.5 * ((b[:, 0] - a[:, 0]) * (c[:, 1] - a[:, 1]) - (c[:, 0] - a[:, 0]) * (b[:, 1] - a[:, 1]))
    flip = area < 0
    conn[flip] = conn[flip][:, [0, 2, 1, 5, 4, 3]]
    nodes[np.abs(nodes[:, 0]) < 1e-12, 0] = 0.0
    return nodes, conn


def boundary_edges(conn):
    """рёбра границы (угол, середина, угол) — встречаются в одном элементе"""
    from collections import Counter
    E = []
    for t in conn:
        E += [(t[0], t[3], t[1]), (t[1], t[4], t[2]), (t[2], t[5], t[0])]
    cnt = Counter(tuple(sorted((a, c))) for a, _, c in E)
    return [e for e in E if cnt[tuple(sorted((e[0], e[2])))] == 1]


def edge_forces_r(nodes, edges, sel, fr_per_len):
    """Радиальные узловые силы (по всей окружности, Н) от погонной нагрузки q(r, z) [Н/м длины ребра по окружности 2πr
    уже включено в q] на выбранных рёбрах: квадратичное ребро — 1/6, 4/6, 1/6."""
    F = np.zeros(len(nodes))
    for a, m, c in edges:
        if not sel(nodes[a], nodes[c]):
            continue
        L = float(np.hypot(*(nodes[c] - nodes[a])))
        qa, qm, qc = fr_per_len(nodes[a]), fr_per_len(nodes[m]), fr_per_len(nodes[c])
        # интеграл N_i q по ребру (Симпсон для квадратичной q на квадратичных функциях формы — достаточно)
        F[a] += L * (qa / 6.0); F[m] += L * (4 * qm / 6.0); F[c] += L * (qc / 6.0)
    return F


def write_inp(path, nodes, conn, mat, T_node, Tref, omega, groups, Fr, fixed_r, fixed_z, conv):
    """groups: список (elset_ids, E_factor, rho_factor). Fr — радиальные силы (полная окружность), делятся на conv."""
    L = ["*HEADING\norganic ccx_check\n*NODE"]
    for i, (r, z) in enumerate(nodes):
        L.append(f"{i + 1}, {r:.9e}, {z:.9e}")
    for gi, (ids, fE, fr) in enumerate(groups):
        if len(ids) == 0:
            continue
        L.append(f"*ELEMENT, TYPE=CAX6, ELSET=G{gi}")
        for e in ids:
            L.append(f"{e + 1}, " + ", ".join(str(n + 1) for n in conn[e]))
    L.append("*ELSET, ELSET=EALL")
    L.extend(f"G{gi}," for gi, g in enumerate(groups) if len(g[0]))
    L.append("*NSET, NSET=NALL")
    ids = np.arange(len(nodes)) + 1
    for k in range(0, len(ids), 16):
        L.append(", ".join(map(str, ids[k:k + 16])) + ",")
    E_tab = [(e * 1e9, T) for T, e in mat["E_GPa"]]
    a_tab = [(a * 1e-6, T) for T, a in mat["alpha_mean_1e6"]]
    for gi, (ids_, fE, fr) in enumerate(groups):
        if len(ids_) == 0:
            continue
        L.append(f"*MATERIAL, NAME=M{gi}\n*ELASTIC")
        L.extend(f"{E * fE:.6e}, {mat['nu']}, {T}" for E, T in E_tab)
        L.append(f"*EXPANSION, ZERO={Tref}")
        L.extend(f"{a:.6e}, {T}" for a, T in a_tab)
        L.append(f"*DENSITY\n{float(mat['density']) * fr:.6f}")
        L.append(f"*SOLID SECTION, ELSET=G{gi}, MATERIAL=M{gi}")
    L.append("*BOUNDARY")
    for n in fixed_r:
        L.append(f"{n + 1}, 1, 1, 0.")
    for n in fixed_z:
        L.append(f"{n + 1}, 2, 2, 0.")
    L.append(f"*INITIAL CONDITIONS, TYPE=TEMPERATURE\nNALL, {Tref}")
    L.append("*STEP\n*STATIC\n*DLOAD")
    L.append(f"EALL, CENTRIF, {omega ** 2:.9e}, 0., 0., 0., 0., 1., 0.")
    nzF = np.nonzero(np.abs(Fr) > 0)[0]
    if len(nzF):
        L.append("*CLOAD")
        for n in nzF:
            L.append(f"{n + 1}, 1, {Fr[n] / conv:.9e}")
    L.append("*TEMPERATURE")
    for i in range(len(nodes)):
        L.append(f"{i + 1}, {T_node[i]:.4f}")
    L.append("*EL PRINT, ELSET=EALL\nS")
    L.append("*NODE FILE\nU, NT\n*EL FILE\nS\n*END STEP")
    with open(path, "w", encoding="utf-8") as f:
        f.write("\n".join(L) + "\n")


def run_ccx(workdir, job="sec"):
    env = dict(os.environ, OMP_NUM_THREADS="1", CCX_NPROC_EQUATION_SOLVER="1", CCX_NPROC_RESULTS="1",
               CCX_NPROC_STIFFNESS="1")
    for ext in (".frd", ".dat"):
        p = os.path.join(workdir, job + ext)
        if os.path.exists(p):
            os.remove(p)
    pr = subprocess.run([DK.ccx_path(), "-i", job], cwd=workdir, capture_output=True, text=True, timeout=3600, env=env)
    log = pr.stdout + pr.stderr
    with open(os.path.join(workdir, "ccx.log"), "w", encoding="utf-8") as f:
        f.write(log)
    if pr.returncode != 0 or "*ERROR" in log:
        raise RuntimeError(f"ccx: ошибка, журнал {workdir}/ccx.log")
    return DK.parse_frd(os.path.join(workdir, job + ".frd")), DK.parse_dat_stress(os.path.join(workdir, job + ".dat"))


def cload_convention(workdir=None):
    """Кольцо r 0,1…0,2 м, толщина 0,02 м, давление 10 МПа на наружной грани: *DLOAD P (через CAX6 грани) против *CLOAD
    той же полной силы. Возврат: множитель conv, на который делится полная сила для *CLOAD (180 — сектор 2°, 2π — радиан)."""
    global _CONV
    workdir = workdir or os.path.join(HERE, "runs", "cload_conv")
    os.makedirs(workdir, exist_ok=True)
    poly = np.array([[0.1, -0.01], [0.2, -0.01], [0.2, 0.01], [0.1, 0.01]])
    nodes, conn = mesh_polygons([poly], 0.004)
    edges = boundary_edges(conn)
    p = 10e6
    sel = lambda a, c: abs(a[0] - 0.2) < 1e-9 and abs(c[0] - 0.2) < 1e-9  # noqa: E731
    Fr = edge_forces_r(nodes, edges, sel, lambda x: p * 2 * math.pi * x[0])
    mat = dict(name="steel", density=7850.0, nu=0.3, E_GPa=[[20, 210.0]], alpha_mean_1e6=[[20, 12.0]])
    fz = [int(np.argmin(np.hypot(nodes[:, 0] - 0.1, nodes[:, 1])))]
    T = np.full(len(nodes), 20.0)
    res = {}
    for conv in (1.0,):
        write_inp(os.path.join(workdir, "sec.inp"), nodes, conn, mat, T, 20.0, 0.0, [(np.arange(len(conn)), 1.0, 1.0)],
                  Fr, [], fz, conv)
        frd, _ = run_ccx(workdir)
        U = frd["DISP"]
        iout = int(np.argmin(np.hypot(nodes[:, 0] - 0.2, nodes[:, 1])))
        res["u_cload_per_full"] = U[iout + 1][0]
    # аналитика Ламе (плоское напряжённое состояние, тонкое кольцо): u(b) = b/E·[(1−ν)A + (1+ν)B/b²]·…
    a_, b_, E, nu = 0.1, 0.2, 210e9, 0.3
    A = p * b_ ** 2 / (b_ ** 2 - a_ ** 2); B = A * a_ ** 2
    u_an = b_ / E * ((1 - nu) * A + (1 + nu) * B / b_ ** 2)
    ratio = u_an / res["u_cload_per_full"]
    conv = min((1.0, 2 * math.pi, 180.0, 360.0), key=lambda c: abs(math.log(ratio * c)))
    # подтверждение: повтор с найденным conv
    write_inp(os.path.join(workdir, "sec.inp"), nodes, conn, mat, T, 20.0, 0.0, [(np.arange(len(conn)), 1.0, 1.0)],
              Fr, [], fz, conv)
    frd, _ = run_ccx(workdir)
    u2 = frd["DISP"][iout + 1][0]
    _CONV = conv
    return dict(ok=bool(abs(u2 / u_an - 1) < 0.02), conv=conv, u_ccx=u2, u_lame_plane_stress=u_an,
                dev=float(u2 / u_an - 1), note="полная сила по окружности делится на conv перед *CLOAD")


def get_conv():
    global _CONV
    if _CONV is None:
        cload_convention()
    return _CONV


def check(name, polys, mat_key, T_of_r, rpm, rim=None, blade=None, lattice=None, rho_field=None, h=1.5e-3,
          fix_z_point=None, Tref=20.0, allow=None, workdir=None, h_min=None):
    """Один расчёт ccx. rim: dict(r, zlo, zhi, p_Pa) — растягивающее давление на грани r; blade: dict(f(r) Н/м радиуса,
    select(a, c) -> bool для рёбер ступицы); rho_field: функция ρ*(r, z) для решётки (None — сплошной).
    allow: функция (T_C, rho*) -> σ_доп, Па. Возврат: максимум σ_экв, загрузка, масса, момент инерции."""
    t0 = time.time()
    workdir = workdir or os.path.join(HERE, "runs", name)
    os.makedirs(workdir, exist_ok=True)
    mat = DK.load_material(mat_key)
    nodes, conn = mesh_polygons(polys, h, h_min)
    edges = boundary_edges(conn)
    omega = 2 * math.pi * rpm / 60.0
    Fr = np.zeros(len(nodes))
    if rim:
        sel = lambda a, c: (abs(a[0] - rim["r"]) < 1e-7 and abs(c[0] - rim["r"]) < 1e-7  # noqa: E731
                            and min(a[1], c[1]) >= rim["zlo"] - 1e-7 and max(a[1], c[1]) <= rim["zhi"] + 1e-7)
        Fr += edge_forces_r(nodes, edges, sel, lambda x: rim["p_Pa"] * 2 * math.pi * x[0])
    if blade:
        # погонная по радиусу сила f(r) -> на ребро ступицы по его радиальной проекции
        for a, m, c in edges:
            if not blade["select"](nodes[a], nodes[c]):
                continue
            ra, rc_ = nodes[a][0], nodes[c][0]
            if abs(rc_ - ra) < 1e-12:
                continue
            rr = np.linspace(min(ra, rc_), max(ra, rc_), 9)
            Ftot = float(np.trapezoid(blade["f"](rr), rr))
            Fr[a] += Ftot / 6; Fr[m] += 4 * Ftot / 6; Fr[c] += Ftot / 6
    cen = (nodes[conn[:, 0]] + nodes[conn[:, 1]] + nodes[conn[:, 2]]) / 3
    a, b, c = nodes[conn[:, 0]], nodes[conn[:, 1]], nodes[conn[:, 2]]
    area = 0.5 * np.abs((b[:, 0] - a[:, 0]) * (c[:, 1] - a[:, 1]) - (c[:, 0] - a[:, 0]) * (b[:, 1] - a[:, 1]))
    if rho_field is not None and lattice is not None:
        rs = np.clip(rho_field(cen[:, 0], cen[:, 1]), lattice["rho_min"], 1.0)
        bins = np.unique(np.round(rs, 2))
        rsb = np.round(rs, 2)
        groups = [(np.nonzero(rsb == bv)[0], float(lattice["E"](bv)), float(bv)) for bv in bins]
        rho_e = rsb
    else:
        groups = [(np.arange(len(conn)), 1.0, 1.0)]
        rho_e = np.ones(len(conn))
    fixed_r = np.nonzero(nodes[:, 0] < 1e-9)[0]
    zp = fix_z_point if fix_z_point is not None else (0.0, 0.0)
    fixed_z = [int(np.argmin(np.hypot(nodes[:, 0] - zp[0], nodes[:, 1] - zp[1])))]
    T_node = np.asarray(T_of_r(nodes[:, 0]), float)
    write_inp(os.path.join(workdir, "sec.inp"), nodes, conn, mat, T_node, Tref, omega, groups, Fr, fixed_r, fixed_z,
              get_conv())
    frd, ips = run_ccx(workdir)
    # напряжения в точках интегрирования
    vm_ip = {}
    for (e, ip), s in ips.items():
        vm_ip.setdefault(e - 1, []).append(DK.von_mises(*s))
    ne = len(conn)
    vmax_e = np.array([max(vm_ip.get(e, [0.0])) for e in range(ne)])
    Te = np.asarray(T_of_r(cen[:, 0]), float)
    out = dict(name=name, elements=int(ne), nodes=int(len(nodes)), rpm=rpm, conv=get_conv(),
               vm_max_ip_MPa=float(vmax_e.max() / 1e6),
               vm_max_at=dict(r_mm=float(cen[vmax_e.argmax(), 0] * 1e3), z_mm=float(cen[vmax_e.argmax(), 1] * 1e3),
                              T_C=float(Te[vmax_e.argmax()])))
    S = frd["STRESS"]
    vm_n = np.array([DK.von_mises(*S[i + 1][:6]) if (i + 1) in S else 0.0 for i in range(len(nodes))])
    out["vm_max_node_MPa"] = float(vm_n.max() / 1e6)
    rho = float(mat["density"])
    out["mass_kg"] = float(np.sum(rho * rho_e * 2 * math.pi * cen[:, 0] * area))
    out["J_kgm2"] = float(np.sum(rho * rho_e * 2 * math.pi * cen[:, 0] ** 3 * area))
    if allow is not None:
        S_l = np.ones(ne) if lattice is None or rho_field is None else lattice["S"](rho_e)
        E_l = np.ones(ne) if lattice is None or rho_field is None else lattice["E"](rho_e)
        # напряжение в решётке — макронапряжение ccx (с уже уменьшенным E) против прочности решётки
        ua = vmax_e / (allow(Te, rho_e) * S_l)
        out["u_max"] = float(ua.max())
        out["u_max_at"] = dict(r_mm=float(cen[ua.argmax(), 0] * 1e3), z_mm=float(cen[ua.argmax(), 1] * 1e3),
                               T_C=float(Te[ua.argmax()]), rho=float(rho_e[ua.argmax()]))
        out["margin_min"] = float(1 / max(out["u_max"], 1e-12))
        _ = E_l
    # средние окружные и запас по разрушению (сплошная часть; решётка — по её доле)
    sig_e = np.zeros((ne, 3))
    cnt = np.zeros(ne)
    for (e, ip), s in ips.items():
        sig_e[e - 1] += (s[0], s[1], s[2]); cnt[e - 1] += 1
    sig_e /= np.maximum(cnt, 1)[:, None]
    hoop = sig_e[:, 2]
    uts = mat_uts(mat, Te) * (np.ones(ne) if lattice is None or rho_field is None else lattice["S"](rho_e))
    out["hoop_mean_MPa"] = float(np.sum(hoop * area) / np.sum(area) / 1e6)
    out["burst_margin"] = float(math.sqrt(np.sum(uts * area) / max(np.sum(hoop * area), 1e-9)))
    lmp = mat.get("_lmp")
    if lmp is not None:
        life, ex = lmp.life_h_v(vmax_e / 1e6, Te)
        out["creep_life_min_h"] = float(np.min(life))
        out["creep_life_min_at"] = dict(r_mm=float(cen[np.argmin(life), 0] * 1e3), z_mm=float(cen[np.argmin(life), 1] * 1e3),
                                        T_C=float(Te[np.argmin(life)]), extrapolated=bool(ex[np.argmin(life)]))
    out["seconds"] = round(time.time() - t0, 1)
    out["ccx"] = DK.ccx_version()
    return out


def mat_uts(mat, T):
    a = np.asarray(mat["uts_MPa"], float)
    return np.interp(T, a[:, 0], a[:, 1]) * 1e6
