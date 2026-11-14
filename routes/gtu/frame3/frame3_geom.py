#!/usr/bin/env python3
"""Типоразмер № 3: геометрия и сетки колёс компрессора ступеней 1 и 2.

Вход - согласованная расчётная точка frame3_point_results.json (stage1, stage2). Узел геометрии -
routes/gtu/platform/geom/impeller.py (Gmsh / Open CASCADE). Для каждой ступени:
  * твёрдое тело колеса (STEP); ступень 1 - двусторонний вход: две зеркальные стороны на общем диске;
  * расчётная область одного межлопаточного канала (основная лопатка + полулопатка) и тетраэдрическая
    сетка SU2 с маркерами inlet outlet hub shroud blade periodic_1 periodic_2, как в узле; у ступени 1 -
    одна сторона на половину расхода;
  * независимая проверка записанных файлов (SU2 разбирается заново numpy, STEP читается обратно).
Большие файлы пишутся в --out (не в репозиторий); в frame3_geom_results.json - их пути и SHA-256.
AWS не вызывается.

Запуск:  python3 frame3_geom.py --out <папка> --levels coarse medium
"""
import argparse
import hashlib
import json
import math
import os
import platform
import sys
import time
from itertools import islice

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
GEOM = os.path.abspath(os.path.join(HERE, "..", "platform", "geom"))
POINT = os.path.join(HERE, "frame3_point_results.json")
RESULTS = os.path.join(HERE, "frame3_geom_results.json")
DATE = "2026-11-14"
STAGES = ("stage1", "stage2")

if GEOM not in sys.path:
    sys.path.insert(0, GEOM)
import impeller as I  # noqa: E402

# поля точки -> параметры impeller (все числа - из frame3_point_results.json, не вписаны руками)
FIELDS = dict(m="m_total_kg_s", PR="PR", rpm="n_rpm", D2="D2_m", T01="T01_C (+273,15 К)", p01="p01_kPa (·1000 Па)",
              Z="blades_full", splitters="blades_splitter > 0", nu="D1_hub_m / D1_tip_m", D1s="D1_tip_m", b2="b2_m",
              Lz="axial_length_m", beta2b="beta2_blade_from_radial_deg", beta1h="beta1_hub_deg",
              beta1s="beta1_tip_deg", gap="tip_clearance_m", double_sided="inlet == «двусторонний»")

# принятые значения этой обёртки (не из точки)
ACCEPTED = dict(
    r_bore=0.0,  # без сквозного отверстия: stage*.material - со сквозным отверстием оба колеса выше σ0,2
)


def rel(path):
    return os.path.relpath(path, REPO)


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def file_info(path):
    return dict(path=path, sha256=sha256(path), bytes=os.path.getsize(path))


def stage_kwargs(s):
    """параметры impeller.build по разделу точки; проверка согласованности самой точки"""
    if s["blades_splitter"] not in (0, s["blades_full"]):
        raise ValueError("impeller.py строит полулопатки только по одной между основными")
    dbl = s["inlet"] == "двусторонний"
    if dbl and abs(s["m_per_side_kg_s"] - s["m_total_kg_s"] / 2) > 1e-9:
        raise ValueError("расход стороны в точке не равен половине расхода колеса")
    kw = dict(m=s["m_total_kg_s"], PR=s["PR"], rpm=s["n_rpm"], D2=s["D2_m"], T01=s["T01_C"] + 273.15,
              p01=s["p01_kPa"] * 1e3, Z=int(s["blades_full"]), splitters=s["blades_splitter"] > 0,
              nu=s["D1_hub_m"] / s["D1_tip_m"], D1s=s["D1_tip_m"], b2=s["b2_m"], Lz=s["axial_length_m"],
              beta2b=s["beta2_blade_from_radial_deg"], beta1h=s["beta1_hub_deg"], beta1s=s["beta1_tip_deg"],
              gap=s["tip_clearance_m"], double_sided=dbl, material="ti")
    kw.update(ACCEPTED)
    return kw


# ----------------------------------------------------------------------------- независимая проверка файлов
def read_su2(path):
    """разбор SU2 (ASCII) построчно, блоки чисел - numpy (файл medium - сотни МБ)"""
    out = {"markers": {}}

    def block(fh, n, dtype):
        lines = list(islice(fh, n))
        ncol = len(lines[0].split())
        return np.fromstring("".join(lines), dtype=dtype, sep=" ").reshape(n, ncol)

    with open(path) as fh:
        for line in fh:
            s = line.strip()
            if s.startswith("NDIME="):
                out["ndime"] = int(s.split("=")[1])
            elif s.startswith("NELEM="):
                out["elem"] = block(fh, int(s.split("=")[1]), np.int64)
            elif s.startswith("NPOIN="):
                out["X"] = block(fh, int(s.split("=")[1].split()[0]), float)[:, :3]
            elif s.startswith("NMARK="):
                out["nmark"] = int(s.split("=")[1])
            elif s.startswith("MARKER_TAG="):
                tag = s.split("=", 1)[1].strip()
                n = int(next(fh).split("=")[1])
                out["markers"][tag] = block(fh, n, np.int64)
    return out


def check_su2(path, imp, Z):
    """размеры и качество - из самого файла SU2, без Gmsh и без JSON узла"""
    from scipy.spatial import cKDTree
    m = read_su2(path)
    E, X = m["elem"], m["X"]
    T = E[:, 1:5]
    a, b, c, d = (X[T[:, k]] for k in range(4))
    V = np.einsum("ij,ij->i", np.cross(b - a, c - a), d - a) / 6.0
    l2 = sum(((P - Q) ** 2).sum(1) for i, P in enumerate((a, b, c, d)) for Q in (a, b, c, d)[i + 1:])
    q = 12.0 * np.cbrt(3.0 * np.abs(V)) ** 2 / l2 * np.sign(V)
    del a, b, c, d, l2
    # граница замкнута и размечена полностью
    n = len(X)
    F = np.sort(np.vstack([T[:, [0, 1, 2]], T[:, [0, 1, 3]], T[:, [0, 2, 3]], T[:, [1, 2, 3]]]), axis=1)
    keys, cnt = np.unique((F[:, 0] * n + F[:, 1]) * n + F[:, 2], return_counts=True)
    del F
    bnd = keys[cnt == 1]
    tagged = np.concatenate([np.sort(e[:, 1:4], axis=1) for e in m["markers"].values()])
    tk = (tagged[:, 0] * n + tagged[:, 1]) * n + tagged[:, 2]
    closed = bool(cnt.max() <= 2 and len(np.unique(tk)) == len(tk) and np.array_equal(np.sort(tk), bnd))
    del keys, cnt, bnd, tagged, tk

    def nodes(k):
        return X[np.unique(m["markers"][k][:, 1:4])]

    def rr(P):
        return np.hypot(P[:, 0], P[:, 1])
    Pin, Pout, Pbl = nodes("inlet"), nodes("outlet"), nodes("blade")
    p1, p2 = m["markers"]["periodic_1"][:, 1:4], m["markers"]["periodic_2"][:, 1:4]
    n1, n2 = np.unique(p1), np.unique(p2)
    ang = 2 * math.pi / Z
    R = np.array([[math.cos(ang), -math.sin(ang), 0], [math.sin(ang), math.cos(ang), 0], [0, 0, 1]])
    dist, j = cKDTree(X[n2]).query(X[n1] @ R.T)
    mp = np.full(n, -1, dtype=np.int64)
    mp[n1] = n2[j]
    k1 = np.sort(mp[p1], axis=1)
    k2 = np.sort(p2, axis=1)
    tri_match = bool(set(map(tuple, k1.tolist())) == set(map(tuple, k2.tolist())))
    # радиальный зазор: узлы лопатки, ближайшие к контуру покрышки (отбор по дереву, затем точно по ломаной)
    zr = np.stack([Pbl[:, 2], rr(Pbl)], 1)
    sh = I._bez(imp.Ps, np.linspace(0, 1, 20001))
    dd, _ = cKDTree(sh).query(zr)
    cand = zr[np.argsort(dd)[:3000]]
    gap_min = float(I._seg_dist(cand, imp.sh_poly).min())
    return dict(
        n_cells=int(len(E)), n_nodes=int(n), only_tetrahedra=bool(np.all(E[:, 0] == 10)),
        markers={k: int(len(v)) for k, v in m["markers"].items()}, n_markers=int(m["nmark"]),
        cell_volume_min_m3=float(V.min()), cells_volume_sum_m3=float(V.sum()), mean_ratio_min=float(q.min()),
        boundary_closed_and_tagged=closed,
        D2_from_blade_nodes_m=float(2 * rr(Pbl).max()),
        b2_from_outlet_z_extent_m=float(Pout[:, 2].max() - Pout[:, 2].min()),
        r_outlet_m=float(rr(Pout).mean()), D1h_from_inlet_m=float(2 * rr(Pin).min()),
        D1s_from_inlet_m=float(2 * rr(Pin).max()), z_inlet_m=float(Pin[:, 2].mean()),
        gap_min_from_blade_nodes_m=gap_min,
        periodic=dict(n_nodes_1=int(len(n1)), n_nodes_2=int(len(n2)), angle_deg=360.0 / Z,
                      max_dist_after_rotation_m=float(dist.max()),
                      one_to_one=bool(len(set(mp[n1].tolist())) == len(n1)), triangles_match=tri_match))


def check_step(path, solid, double_sided, threads=4):
    """STEP колеса читается обратно: одно тело; объём тремя способами (по частям при построении, getMass OCC
    прочитанного тела, сумма объёмов его тетраэдрической сетки); наибольший радиус; зеркальность по узлам рёбер"""
    import gmsh
    from scipy.spatial import cKDTree
    gmsh.initialize()
    try:
        gmsh.option.setNumber("General.Terminal", 0)
        gmsh.option.setNumber("General.NumThreads", threads)
        gmsh.model.occ.importShapes(path)
        gmsh.model.occ.synchronize()
        vols = gmsh.model.getEntities(3)
        V = gmsh.model.occ.getMass(3, vols[0][1]) if len(vols) == 1 else float("nan")
        gmsh.option.setNumber("Mesh.MeshSizeMax", 0.004)
        gmsh.model.mesh.generate(1)
        X = gmsh.model.mesh.getNodes()[1].reshape(-1, 3)
        # объём тетраэдрической сетки (размеры - как в test_impeller.py, проверка 3)
        gmsh.model.mesh.clear()
        for k, v in (("Mesh.MeshSizeMax", 0.01), ("Mesh.MeshSizeMin", 0.0004), ("Mesh.MeshSizeFromCurvature", 16),
                     ("Mesh.Algorithm3D", 10), ("Mesh.MaxNumThreads3D", threads)):
            gmsh.option.setNumber(k, v)
        gmsh.model.mesh.generate(3)
        tags, Xn, _ = gmsh.model.mesh.getNodes()
        Xn = Xn.reshape(-1, 3)
        idx = np.zeros(int(tags.max()) + 1, dtype=np.int64)
        idx[tags.astype(np.int64)] = np.arange(len(tags))
        _, _, en = gmsh.model.mesh.getElements(3)
        Tt = idx[np.concatenate(en).astype(np.int64)].reshape(-1, 4)
        a, b, c, d = (Xn[Tt[:, k]] for k in range(4))
        Vt = float(np.abs(np.einsum("ij,ij->i", np.cross(b - a, c - a), d - a)).sum() / 6.0)
        r_max_3d = float(np.hypot(Xn[:, 0], Xn[:, 1]).max())
        n_tet = int(len(Tt))
    finally:
        gmsh.finalize()
    out = dict(n_volumes=len(vols), volume_occ_step_m3=V, volume_by_parts_m3=solid["volume_m3"],
               volume_rel_diff=V / solid["volume_m3"] - 1.0, volume_tet_mesh_m3=Vt, n_tet_check=n_tet,
               volume_tet_rel_diff=Vt / solid["volume_m3"] - 1.0,
               r_max_edge_nodes_m=float(np.hypot(X[:, 0], X[:, 1]).max()), r_max_tet_nodes_m=r_max_3d,
               z_min_m=float(X[:, 2].min()), z_max_m=float(X[:, 2].max()))
    if double_sided:
        zc = solid["z_mid_plane"]
        Y = X.copy()
        Y[:, 2] = 2 * zc - Y[:, 2]
        dist, _ = cKDTree(X).query(Y)
        out.update(z_mid_plane_m=zc, z_symmetry_error_m=float(abs(X[:, 2].min() + X[:, 2].max() - 2 * zc)),
                   mirror_edge_nodes_median_dist_m=float(np.median(dist)))
    return out


# ----------------------------------------------------------------------------- расчёт
def run_stage(key, s, levels, out, threads):
    kw = stage_kwargs(s)
    res_st = dict(label=s["label"], inlet=s["inlet"], params_from_point={k: kw[k] for k in FIELDS},
                  accepted=dict(ACCEPTED), levels={})
    imp = I.Impeller(**kw)  # те же поверхности, что у build: контур покрышки для проверки зазора по файлу
    for i, level in enumerate(levels):
        od = os.path.join(out, f"{key}_{level}")
        t0 = time.time()
        res = I.build("custom", level, out_dir=od, step=(i == 0), mesh=True, msh=(level == "coarse"), threads=threads,
                      source=f"frame3_point_results.json → {key}", **kw)
        t_build = time.time() - t0
        if i == 0:
            d1 = res["d1"]
            solid = res["solid"]
            res_st.update(
                d1={k: d1[k] for k in ("m", "m_side", "double_sided", "PR", "rpm", "D2", "r2", "U2", "omega", "T01", "p01",
                                       "b2", "gap", "r1s", "r1h", "D1s_D2", "Lz", "Ls", "beta1h", "beta1s", "beta2b",
                                       "beta1h_flow", "beta1s_flow", "C1", "M1s_rel", "Z_exit", "wiesner_eps_lim")},
                d1_comp_stage_not_used_for_geometry={k: d1[k] for k in ("w", "T2", "psi", "psi_euler", "eta_comp_stage",
                                                                        "slip_wiesner", "phi")},
                geometry=res["geometry"], area=res["area"], contours=res["contours"], warnings=res["warnings"],
                solid={k: v for k, v in solid.items() if k not in ("step", "bbox_occ")},  # рамка OCC грубая
                solid_file=file_info(res["files"]["impeller_step"]),
                su2=res["su2"])
            t1 = time.time()
            res_st["solid_check"] = check_step(res["files"]["impeller_step"], solid, kw["double_sided"], threads)
            res_st["solid_check"]["t_s"] = round(time.time() - t1, 1)
        t1 = time.time()
        chk = check_su2(res["files"]["passage_su2"], imp, kw["Z"])
        chk["t_s"] = round(time.time() - t1, 1)
        mesh = {k: v for k, v in res["mesh"].items() if k not in ("gmsh_warnings", "faces")}
        files = {k: file_info(v) for k, v in res["files"].items() if k != "json" and k != "impeller_step"}
        files["impeller_json"] = file_info(res["files"]["json"])
        res_st["levels"][level] = dict(
            n_cells=mesh["quality"]["n_cells"], n_nodes=mesh["quality"]["n_nodes"], minSICN=mesh["quality"]["minSICN"],
            min_gamma=mesh["quality"]["min_gamma"], mesh=mesh, n_gmsh_warnings=len(res["mesh"]["gmsh_warnings"]),
            t_build_s=round(t_build, 1), check_su2=chk, files=files)
        q = mesh["quality"]
        print(f"[{key} {level}] ячеек {q['n_cells']}, узлов {q['n_nodes']}, minSICN {q['minSICN']:.3f}, "
              f"{t_build:.0f} с; проверка файла {chk['t_s']} с", flush=True)
        json.dump(res_st, open(os.path.join(out, f"{key}_partial.json"), "w", encoding="utf-8"), ensure_ascii=False,
                  indent=1, default=float)
    return res_st


def verdicts(st, s):
    """сводка проверок: заданное значение против прочитанного из файлов"""
    v = {}
    D2, b2 = s["D2_m"], s["b2_m"]
    for level, L in st["levels"].items():
        c = L["check_su2"]
        v[level] = dict(
            D2_rel_err=c["D2_from_blade_nodes_m"] / D2 - 1.0,
            b2_rel_err=c["b2_from_outlet_z_extent_m"] / b2 - 1.0,
            D1s_rel_err=c["D1s_from_inlet_m"] / s["D1_tip_m"] - 1.0,
            D1h_rel_err=c["D1h_from_inlet_m"] / s["D1_hub_m"] - 1.0,
            gap_rel_err=c["gap_min_from_blade_nodes_m"] / s["tip_clearance_m"] - 1.0,
            r_outlet_over_r2=c["r_outlet_m"] / (D2 / 2),
            cells_volume_vs_occ_rel=c["cells_volume_sum_m3"] / L["mesh"]["volume_fluid_occ_m3"] - 1.0,
            minSICN_gmsh_vs_mean_ratio_numpy=[L["minSICN"], c["mean_ratio_min"]],
            n_cells_gmsh_vs_file=[L["n_cells"], c["n_cells"]],
            markers_7=sorted(c["markers"]) == sorted(I.MARKERS) and all(n > 0 for n in c["markers"].values()),
            ok=bool(c["only_tetrahedra"] and c["cell_volume_min_m3"] > 0 and c["boundary_closed_and_tagged"]
                    and c["periodic"]["one_to_one"] and c["periodic"]["triangles_match"]
                    and c["periodic"]["max_dist_after_rotation_m"] < 1e-10 * D2
                    and abs(c["D2_from_blade_nodes_m"] / D2 - 1) < 1e-5 and abs(c["b2_from_outlet_z_extent_m"] / b2 - 1) < 1e-4
                    and abs(c["gap_min_from_blade_nodes_m"] / s["tip_clearance_m"] - 1) < 0.02
                    and L["minSICN"] > 0 and c["mean_ratio_min"] > 0))
    sc = st["solid_check"]
    # допуски: getMass OCC - 1e-3 (точность интегрирования по сложным граням, README узла geom), объём сетки -
    # 2e-3 (линейные тетраэдры на кривых гранях, как test_impeller.py)
    v["solid"] = dict(one_body=sc["n_volumes"] == 1, volume_occ_rel_diff=sc["volume_rel_diff"],
                      volume_tet_rel_diff=sc["volume_tet_rel_diff"], r_max_over_r2=sc["r_max_edge_nodes_m"] / (D2 / 2),
                      ok=bool(sc["n_volumes"] == 1 and abs(sc["volume_rel_diff"]) < 1e-3
                              and abs(sc["volume_tet_rel_diff"]) < 2e-3
                              and abs(sc["r_max_edge_nodes_m"] / (D2 / 2) - 1) < 1e-5
                              and sc["r_max_tet_nodes_m"] < (D2 / 2) * (1 + 1e-6)
                              and (not st["d1"]["double_sided"] or sc["z_symmetry_error_m"] < 1e-6)))
    return v


def main(argv=None):
    ap = argparse.ArgumentParser(description="KV-3: STEP и сетки канала SU2 колёс компрессора ступеней 1 и 2")
    ap.add_argument("--out", default=os.path.join(GEOM, "out", "frame3"),
                    help="папка больших файлов (по умолчанию platform/geom/out/frame3, в .gitignore)")
    ap.add_argument("--levels", nargs="+", default=["coarse"], choices=list(I.LEVELS))
    ap.add_argument("--stages", nargs="+", default=list(STAGES), choices=list(STAGES))
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--results", default=RESULTS)
    a = ap.parse_args(argv)
    os.makedirs(a.out, exist_ok=True)
    point = json.load(open(POINT, encoding="utf-8"))
    import gmsh
    import scipy
    try:
        import cantera
        ct = cantera.__version__
    except Exception as e:  # noqa
        ct = f"не импортируется: {e}"
    res = dict(
        title="Типоразмер № 3: геометрия и сетки колёс компрессора ступеней 1 и 2",
        method=("impeller.py (Gmsh / Open CASCADE): меридиональный контур - кубические Безье, линейчатая лопатка, "
                "закон углов beta(m), полулопатки; размеры D2, D1, b2, осевая длина, углы лопатки на входе и выходе, "
                "число лопаток, зазор, частота, полные параметры на входе - из согласованной точки; ступень 1 - "
                "двусторонний вход: канал и сетка одной стороны на половину расхода, STEP - обе стороны зеркально. "
                "Тетраэдрическая сетка одного межлопаточного канала (основная лопатка + полулопатка) с периодическими "
                "гранями, формат SU2. Файлы проверены независимо: SU2 разобран numpy, STEP прочитан обратно Gmsh"),
        software=dict(python=sys.version.split()[0], numpy=np.__version__, scipy=scipy.__version__, gmsh=gmsh.__version__,
                      cantera=ct, os=platform.platform()),
        date=DATE,
        status=("геометрия уровня облика по одномерной точке (оценка): контуры, толщины, закон углов - принятые "
                "значения узла; расчёта течения (SU2) и прочности (CalculiX) по этим файлам нет"),
        inputs=[
            dict(file=rel(POINT), sha256=sha256(POINT), field="stage1, stage2: " + "; ".join(f"{k} ← {v}" for k, v in FIELDS.items()),
                 note="точка проверена (pass_after_fixes); при пересчёте точки по замечаниям проверяющего "
                      "(π1 около 4,70) геометрию надо перестроить этим скриптом"),
            dict(file=rel(os.path.join(GEOM, "impeller.py")), sha256=sha256(os.path.join(GEOM, "impeller.py")),
                 field="DEFAULTS (толщины, контуры Безье, m_split 0,30, длины участков области), LEVELS",
                 used=dict(DEFAULTS={k: (list(v) if isinstance(v, tuple) else v) for k, v in I.DEFAULTS.items()},
                           LEVELS=I.LEVELS)),
        ],
        assumptions=dict(
            r_bore="0: колёса без сквозного отверстия (frame3_point_results.json → stage*.material: со сквозным "
                   "отверстием vm_bore_to_ys 1,05 и 1,36, без него 0,58 и 0,75). В точке nu = 0,4 принят под вал, "
                   "проходящий через оба колеса (frame3_point.py, A['nu']) - противоречие точки, крепление колёс без "
                   "отверстия (фланцы, хирт) не прорабатывалось",
            double_sided="ступень 1: две стороны зеркально относительно средней плоскости общего диска, лопатки "
                         "сторон напротив друг друга (принято); толщина общего диска у r2 - max(0,03 D2; 0,75 b2) = "
                         "2·solid.t_back; канал и сетка - одна сторона на m/2; безлопаточный участок за колесом "
                         "шириной b2 одной стороны, стенка z = Lz - в маркере hub",
            m_split="0,30 (полулопатка - 0,70 длины): совпадает с весом разделительной лопатки spl = 0,7 точки",
            inlet_angles="углы лопатки на входе - углы потока точки (угол атаки 0; в точке загромождение входа "
                         "k1 = 0,97); свой треугольник impeller без загромождения - в d1.beta1*_flow",
            other="толщины лопаток (0,011 и 0,0045 D2), контуры Безье 0,5/0,5, закон углов, входной участок 1,0 r1s, "
                  "выход области 1,3 r2 - значения узла по умолчанию, не результат оптимизации и не проверены "
                  "прочностью"),
        stages={}, checks={},
    )
    t0 = time.time()
    for key in a.stages:
        st = run_stage(key, point[key], a.levels, a.out, a.threads)
        res["stages"][key] = st
        res["checks"][key] = verdicts(st, point[key])
    res["t_total_s"] = round(time.time() - t0, 1)
    res["not_done"] = [
        "SU2 эти сетки не читал (в среде нет сошедшихся расчётов колеса, см. навык su2-turbo, раздел 4); формат "
        "проверен своим разбором",
        "пристенных призматических слоёв нет (как в узле): для RANS с y+ ≈ 1 сетку надо дополнять",
        "маркер hub включает неподвижную стенку безлопаточного участка (z = Lz, r2 < r ≤ 1,3 r2): во вращающейся "
        "системе она будет вращаться. У ступени 1 вопрос проверяющего о ширине диффузора (общий 2·b2 или две "
        "половины) не решён: сетка одной стороны отвечает «двум половинам» со стенкой z = Lz; для общего диффузора "
        "с тонкой кромкой диска эту стенку нужно задать плоскостью симметрии - для этого маркер надо разделить",
        "входной участок - прямой осевой, как в узле; подвод к двустороннему колесу (кольцевые камеры с двух "
        "сторон) не моделируется",
        "уровень fine не строился (по росту coarse → medium - порядка 2·10⁷ ячеек и более, для AWS)",
        "прочность колёс без отверстия и частоты по этим STEP - узел CalculiX, не сделано",
    ]
    if any(lv not in a.levels for lv in ("coarse", "medium")):
        res["not_done"].append("уровни, не указанные в --levels: " + ", ".join(lv for lv in ("coarse", "medium")
                                                                               if lv not in a.levels))
    with open(a.results, "w", encoding="utf-8") as fh:
        json.dump(res, fh, ensure_ascii=False, indent=1, default=float)
    print("итог:", a.results, f"({res['t_total_s']} с)")
    for key, c in res["checks"].items():
        print(key, {k: v.get("ok") for k, v in c.items()})
    return res


if __name__ == "__main__":
    main()
