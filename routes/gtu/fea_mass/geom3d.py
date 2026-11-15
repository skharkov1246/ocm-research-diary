"""Массовая 3D-прочность (уровень 2 воронки): генератор трёхмерных секторов роторных деталей для CalculiX.

Что строится (Gmsh / Open CASCADE, тетраэдры второго порядка C3D10, периодическая сетка «узел в узел»):
  * wheel_sector - сектор 360/Z колеса заодно с лопатками (центробежное колесо компрессора; радиальная турбина - тем же
    генератором: радиальные лопатки на входе, эксдьюсер на месте входа колеса). Лопатки - ../platform/geom/impeller.py
    (основная + полулопатка), ступица - тело вращения из impeller.py (нос, втулка, задний диск, отверстие). Лопатки колеса
    закручены на 40-70 град при шаге 25-40 град, поэтому граница сектора по ступице - не плоскость, а закрученная
    линейчатая поверхность: отрезки постоянного угла от линии втулки к отверстию, угол - середина между основной лопаткой и
    полулопаткой соседнего канала (как периодическая поверхность канала в impeller.build_passage). Лопатки в сектор
    входят целиком: их корни лежат строго внутри сектора ступицы, периодические грани - только по ступице.
  * axial_sector - сектор 360/Z осевого колеса: диск (профиль равнопрочного диска уровня 1 - strength_l1.disc_axial),
    полка, перо (дуговая средняя линия по углам одномерной модели, свободный вихрь по высоте, толщина NACA 00xx) и
    замок: для вставных лопаток - ножка + хвостовик лопатки и выступ диска, связанные с живым ободом (контакт не
    моделируется, окружной связи между лопаткой и выступом нет - как у паза); для блиска (литое или печатное колесо) -
    перо прямо на ободе.
  * ring_sector - кольцо постоянной толщины (самопроверка на аналитике); beam - консольная балка (частоты).

Выход - объект Mesh3D: узлы, элементы C3D10 (порядок узлов ccx), наборы узлов (NSLAVE, NMASTER, NBORE, NFIX), класс
элемента (ступица / лопатка / замок) для осреднений, сведения о сетке и проверки (объём сектора = 1/Z тела вращения,
совпадение узлов периодических граней). Все размеры, кроме взятых из записи уровня 1, - оценка (словарь A).
"""
import math
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
GTU = os.path.dirname(HERE)
GEOM = os.path.join(GTU, "platform", "geom")

# ----------------------------------------------------------------------------- принятые величины (все - оценка)
A = dict(
    # сетка: размер элемента в лопатке = k_blade · средняя толщина лопатки; в ступице - k_hub · размер в лопатке,
    # рост от лопатки к ступице - grow на единицу расстояния (доля размера)
    levels=dict(coarse=dict(k_blade=1.6, k_hub=4.0, grow=0.35),
                base=dict(k_blade=1.0, k_hub=4.0, grow=0.30),
                fine=dict(k_blade=0.65, k_hub=4.0, grow=0.25)),
    core_bore_r1h=0.15,        # цельное колесо (заодно, сварка): технологическое отверстие r = 0,15 r1h (ось внутри
    #                            сектора ccx не держит; пик у этого отверстия из оценок исключается - концентрация ×2 ложная)
    tie_bore_r1h=0.35,         # колесо сборного ротора: отверстие под стяжку r = 0,35 r1h (как impeller.DEFAULTS r_bore)
    # осевое колесо
    t_max_c_hub=0.20, t_max_c_tip=0.08,    # наибольшая толщина пера / хорда у корня и на периферии
    t_te_c=0.015,              # толщина выходной кромки / хорда
    plat_t_c=0.10,             # толщина полки / хорда пера
    plat_len_bx=1.25,          # осевая длина полки / осевая хорда
    plat_gap_pitch=0.03,       # зазор между полками соседних лопаток / шаг (полки не связаны окружно)
    neck_w_pitch=0.30,         # толщина ножки / шаг на ободе (как frame3_fea turb_neck_s)
    neck_h_c=0.40,             # высота ножки / хорда
    foot_w_pitch=0.50,         # ширина хвостовика / шаг
    foot_h_bx=0.60,            # высота хвостовика / осевая хорда (как frame3_fea turb_h_ft_cax)
    lug_gap_pitch=0.04,        # зазор между хвостовиком лопатки и выступом диска / шаг
    rim_w_bx=1.15,             # ширина обода / осевая хорда (frame3_fea turb_w_rim_cax)
    disc_bore_R=0.20,          # радиус отверстия диска / радиус живого обода (как funnel_l2.R_BORE)
)

_IMP = None


def impeller_module():
    global _IMP
    if _IMP is None:
        if GEOM not in sys.path:
            sys.path.insert(0, GEOM)
        import impeller as IM            # noqa: E402
        _IMP = IM
    return _IMP


def _gmsh():
    import gmsh
    return gmsh


def _init(gmsh, threads=1):
    gmsh.initialize(interruptible=False)
    gmsh.option.setNumber("General.Terminal", 0)
    gmsh.option.setNumber("General.Verbosity", 2)
    gmsh.option.setNumber("General.NumThreads", threads)
    gmsh.logger.start()


GMSH_TO_CCX_TET10 = [0, 1, 2, 3, 4, 5, 6, 7, 9, 8]   # Gmsh: 8 - ребро 2-3, 9 - ребро 1-3; ccx/Abaqus C3D10: 8 - 1-3, 9 - 2-3


class Mesh3D:
    """сетка сектора: nodes (N,3) м, node_ids (N,), elems (E,10) номера узлов в порядке ccx, ecls (E,) класс элемента"""

    def __init__(self, nodes, node_ids, elems, ecls, sets, meta):
        self.nodes, self.node_ids, self.elems, self.ecls, self.sets, self.meta = nodes, node_ids, elems, ecls, sets, meta
        self.index = {int(t): i for i, t in enumerate(node_ids)}

    def xyz(self, ids):
        return self.nodes[[self.index[int(t)] for t in ids]]

    def write_inp(self, path):
        """узлы, элементы C3D10 по классам (ELSET E<класс>, общий EALL), наборы узлов"""
        with open(path, "w") as f:
            f.write("*NODE, NSET=NALL\n")
            for t, (x, y, z) in zip(self.node_ids, self.nodes):
                f.write(f"{int(t)}, {x:.10e}, {y:.10e}, {z:.10e}\n")
            classes = sorted(set(int(c) for c in self.ecls))
            eid = 0
            for c in classes:
                f.write(f"*ELEMENT, TYPE=C3D10, ELSET=E{CLASS_NAMES[c]}\n")
                for row in self.elems[self.ecls == c]:
                    eid += 1
                    f.write(f"{eid}, " + ", ".join(str(int(v)) for v in row) + "\n")
            f.write("*ELSET, ELSET=EALL\n" + ", ".join(f"E{CLASS_NAMES[c]}" for c in classes) + "\n")
            for name, ids in self.sets.items():
                f.write(f"*NSET, NSET={name}\n")
                ids = [int(v) for v in ids]
                for k in range(0, len(ids), 12):
                    f.write(", ".join(str(v) for v in ids[k:k + 12]) + ",\n")
        return path

    def element_order(self):
        """номер элемента ccx (1..E) -> индекс строки elems: порядок записи по классам"""
        order = []
        for c in sorted(set(int(c) for c in self.ecls)):
            order.extend(np.nonzero(self.ecls == c)[0].tolist())
        return np.asarray(order)


CLASS_HUB, CLASS_BLADE, CLASS_ATT = 0, 1, 2
CLASS_NAMES = {CLASS_HUB: "HUB", CLASS_BLADE: "BLADE", CLASS_ATT: "ATT"}


# ----------------------------------------------------------------------------- общие шаги Gmsh
def _rotT(a):
    c, s = math.cos(a), math.sin(a)
    return [c, -s, 0, 0, s, c, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1]


def _face_points(gmsh, f, n=4):
    """точки на рёбрах грани (центр масс кривой грани на ней не лежит)"""
    pts = []
    for _, c in gmsh.model.getBoundary([(2, f)], combined=False, oriented=False):
        c = abs(c)
        lo, hi = gmsh.model.getParametrizationBounds(1, c)
        u = np.linspace(lo[0], hi[0], n + 2)[1:-1]
        pts.append(np.asarray(gmsh.model.getValue(1, c, list(u))).reshape(-1, 3))
    return np.vstack(pts)


def _on_surface(gmsh, ref, P, tol):
    return max(np.linalg.norm(np.asarray(gmsh.model.getClosestPoint(2, ref, list(x))[0]) - x) for x in P) < tol


def _classify_periodic(gmsh, vol, refs_master, refs_slave, tol):
    """грани объёма на опорных поверхностях: master - на refs_master, slave - на refs_slave"""
    gmsh.model.occ.synchronize()
    master, slave = [], []
    for _, f in gmsh.model.getBoundary([(3, vol)], oriented=False):
        f = abs(f)
        P = _face_points(gmsh, f)
        if any(_on_surface(gmsh, r, P, tol) for r in refs_master):
            master.append(f)
        elif any(_on_surface(gmsh, r, P, tol) for r in refs_slave):
            slave.append(f)
    return master, slave


def _pair_faces(gmsh, master, slave, pitch, tol):
    """пары (slave, master): центр масс master, повёрнутый на pitch, совпадает с центром slave"""
    cm = {f: np.asarray(gmsh.model.occ.getCenterOfMass(2, f)) for f in master + slave}
    R = np.asarray(_rotT(pitch)).reshape(4, 4)[:3, :3]
    pairs = []
    for fs in slave:
        best = min(master, key=lambda fm: np.linalg.norm(R @ cm[fm] - cm[fs]))
        d = np.linalg.norm(R @ cm[best] - cm[fs])
        if d > tol:
            raise RuntimeError(f"периодическая грань {fs} без пары (ближайшая {best}, {d:.2e} м)")
        pairs.append((fs, best))
    if len({m for _, m in pairs}) != len(master) or len(slave) != len(master):
        raise RuntimeError(f"периодические грани не попарны: master {len(master)}, slave {len(slave)}")
    return pairs


def _snap(gmsh, slaves):
    IM = impeller_module()
    return IM._snap_periodic(gmsh, slaves)


def _quality(gmsh):
    typ, tags, _ = gmsh.model.mesh.getElements(3)
    allt = np.concatenate([np.asarray(t) for t in tags]) if tags else np.zeros(0)
    sicn = np.asarray(gmsh.model.mesh.getElementQualities(list(allt), "minSICN"))
    return dict(n_tets=int(len(allt)), minSICN=float(sicn.min()), p01SICN=float(np.percentile(sicn, 1)),
                n_SICN_below_0_1=int((sicn < 0.1).sum()))


def _mesh_volume(gmsh, vol, pairs, pitch, size_cb, h_min, h_max, threads=1, opt_netgen=True):
    """периодическая сетка C3D10 одного объёма; size_cb(x, y, z) -> размер элемента"""
    for fs, fm in pairs:
        gmsh.model.mesh.setPeriodic(2, [fs], [fm], _rotT(pitch))
    gmsh.model.mesh.setSizeCallback(lambda dim, tag, x, y, z, lc: max(h_min, min(h_max, size_cb(x, y, z))))
    for k, v in (("Mesh.MeshSizeExtendFromBoundary", 0), ("Mesh.MeshSizeFromPoints", 0), ("Mesh.MeshSizeFromCurvature", 0),
                 ("Mesh.MeshSizeMin", h_min), ("Mesh.MeshSizeMax", h_max), ("Mesh.Algorithm", 6), ("Mesh.Algorithm3D", 1),
                 ("Mesh.Optimize", 1), ("Mesh.ElementOrder", 1), ("Mesh.MaxNumThreads3D", threads)):
        gmsh.option.setNumber(k, v)
    gmsh.model.addPhysicalGroup(3, [vol], name="EVOL")
    gmsh.model.mesh.generate(3)
    q = _quality(gmsh)
    if q["minSICN"] < 0.05 and opt_netgen:
        gmsh.model.mesh.optimize("Netgen")
        q = _quality(gmsh)
    gmsh.option.setNumber("Mesh.HighOrderOptimize", 0)
    gmsh.option.setNumber("Mesh.SecondOrderLinear", 1)   # средние узлы - на серединах прямых рёбер: без вывернутых элементов
    gmsh.model.mesh.setOrder(2)
    snap = _snap(gmsh, [fs for fs, _ in pairs])
    q2 = _quality(gmsh)
    q.update(minSICN_order2=q2["minSICN"], p01SICN_order2=q2["p01SICN"])
    if q2["minSICN"] <= 0:
        raise RuntimeError(f"вырожденные элементы второго порядка: minSICN = {q2['minSICN']:.3g}")
    return q, snap


def _extract(gmsh, vol, master_faces, slave_faces, bore_faces):
    """узлы и элементы C3D10 из Gmsh -> массивы; наборы узлов граней"""
    tags, xyz, _ = gmsh.model.mesh.getNodes()
    tags = np.asarray(tags, dtype=np.int64)
    xyz = np.asarray(xyz).reshape(-1, 3)
    et, etags, enodes = gmsh.model.mesh.getElements(3, vol)
    k = [i for i, t in enumerate(et) if int(t) == 11]
    if not k:
        raise RuntimeError("в объёме нет тетраэдров второго порядка (тип 11)")
    conn = np.asarray(enodes[k[0]], dtype=np.int64).reshape(-1, 10)[:, GMSH_TO_CCX_TET10]
    used = np.unique(conn)
    # компактная нумерация 1..N только по узлам элементов
    new = {int(t): i + 1 for i, t in enumerate(used)}
    pos = {int(t): i for i, t in enumerate(tags)}
    nodes = xyz[[pos[int(t)] for t in used]]
    elems = np.vectorize(new.get)(conn)

    def nset(faces):
        s = set()
        for f in faces:
            nt = gmsh.model.mesh.getNodes(2, f, includeBoundary=True)[0]
            s |= {new[int(t)] for t in nt if int(t) in new}
        return sorted(s)
    return nodes, np.arange(1, len(used) + 1), elems, nset(master_faces), nset(slave_faces), nset(bore_faces)


def _match_periodic(nodes, master, slave, pitch, tol):
    """каждому узлу slave - узел master, повёрнутый на pitch (проверка «узел в узел»)"""
    R = np.asarray(_rotT(pitch)).reshape(4, 4)[:3, :3]
    M = nodes[np.asarray(master) - 1] @ R.T
    S = nodes[np.asarray(slave) - 1]
    from scipy.spatial import cKDTree
    d, _ = cKDTree(M).query(S)
    return dict(n_master=len(master), n_slave=len(slave), max_dist_m=float(d.max()) if len(d) else None,
                ok=bool(len(master) == len(slave) and len(d) and d.max() < tol))


def _fix_node(nodes, bore, slave, master, zc):
    """узел для закрепления статики (окружное и осевое): на отверстии, ближе всех к середине по оси, не на slave и
    не на master"""
    cand = sorted(set(bore) - set(slave) - set(master))
    if not cand:
        raise RuntimeError("нет узла отверстия вне периодических граней для закрепления статики")
    P = nodes[np.asarray(cand) - 1]
    return int(cand[int(np.argmin(np.abs(P[:, 2] - zc) + 1e-9 * np.arange(len(cand))))])


# ----------------------------------------------------------------------------- колесо заодно с лопатками
E_ROOT_B2 = 0.06              # заглубление корня лопатки в ступицу / b2 (impeller.py берёт 0,25 b2 - для канала; в твёрдом
#                               теле глубокий корень пересекал закрученную границу сектора у выхода; 0,06 b2 - только чтобы
#                               тела надёжно сплавились, оценка)


def _blade_sections(imp, m0=0.0, rot=0.0, e_root_b2=E_ROOT_B2):
    """как impeller.Impeller.blade_sections, но с заглублением корня e_root_b2·b2"""
    IM = impeller_module()
    d, p = imp.d1, imp.p
    e_root = e_root_b2 * d["b2"]
    e_far = 0.25 * d["b2"]
    m_end = 1 + 1.5 * p["t_hub"] * d["D2"] / d["r2"]
    n = p["n_st"] if m0 == 0 else max(12, int(round(p["n_st"] * (1 - m0))))
    u = np.linspace(0, 1, n)
    u = 0.4 * u + 0.6 * 0.5 * (1 - np.cos(np.pi * u))
    ms = m0 + (1 - m0) * u
    ms = np.concatenate([ms, [1 + 0.5 * (m_end - 1), m_end]])
    secs = []
    for m in ms:
        H, S = imp.hs(m)
        H, S = H[0], S[0]
        dv = S - H
        L = np.linalg.norm(dv)
        nv = imp.normal(m)
        sr, sf = -e_root / L, 1 + e_far / L
        tr, tf = imp.thickness(min(m, 1.0), sr, m0), imp.thickness(min(m, 1.0), sf, m0)
        Cr, Cf = H + sr * dv, H + sf * dv
        q = np.array([Cr + 0.5 * tr * nv, Cf + 0.5 * tf * nv, Cf - 0.5 * tf * nv, Cr - 0.5 * tr * nv])
        secs.append(IM._rotz(q, rot))
    return secs


def _blade_solids(occ, imp):
    """основная лопатка и полулопатка, обрезанные по покрышке минус зазор и по r2 (как impeller._blade_solids)"""
    IM = impeller_module()
    p = imp.p
    blades = [IM._loft_solid(occ, _blade_sections(imp, 0.0, 0.0))]
    if p["splitters"]:
        blades.append(IM._loft_solid(occ, _blade_sections(imp, p["m_split"], math.pi / p["Z"])))
    allp = np.vstack(_blade_sections(imp, 0.0, 0.0))
    th = np.unwrap(np.arctan2(allp[:, 1], allp[:, 0]))
    lo, hi = th.min() - 0.3, th.max() + 0.3 + (math.pi / p["Z"] if p["splitters"] else 0.0)
    env = IM._envelope(occ, imp, lo, hi - lo)
    out, _ = occ.intersect([(3, b) for b in blades], [(3, env)])
    vols = [t for dd, t in out if dd == 3]
    if len(vols) != len(blades):
        raise RuntimeError(f"обрезка лопаток дала {len(vols)} тел вместо {len(blades)}")
    return vols


def make_impeller(spec):
    """объект impeller.Impeller по спецификации детали (parts3d.wheel_spec)"""
    IM = impeller_module()
    kw = dict(Z=int(spec["Z"]), splitters=bool(spec["splitters"]), b2=spec["b2"], beta2b=spec["beta2b"], D1s=spec["D1s"],
              nu=spec["nu"], gap=spec.get("gap"), double_sided=bool(spec.get("double_sided")))
    for k in ("Lz", "beta1h", "beta1s", "t_hub", "t_tip"):
        if spec.get(k) is not None:
            kw[k] = spec[k]
    return IM.Impeller(m=spec["m"], PR=spec["PR"], rpm=spec["rpm"], D2=spec["D2"], case=spec.get("name", "part"), material="ti",
                       source="запись уровня 1", T01=spec.get("T01"), p01=spec.get("p01"), **kw)


def _hub_chain(imp, z_n, z_bk, r_b, e):
    """ломаная наружной границы ступицы (z, r, угол средней линии втулки) и отрезок отверстия; продолжение на e за торцы"""
    d = imp.d1
    m = np.linspace(0, 1, 121)
    zr = imp.hub.at(m)
    th = imp.th_h(m)
    nose = np.array([[z_n - e, d["r1h"]], [z_n, d["r1h"]], [0.5 * z_n, d["r1h"]]])
    back = np.array([[0.5 * (d["Lz"] + z_bk), d["r2"]], [z_bk, d["r2"]], [z_bk + e, d["r2"]]])
    outer = np.vstack([nose, zr, back])
    tho = np.concatenate([[th[0]] * len(nose), th, [th[-1]] * len(back)])
    s = np.concatenate([[0.0], np.cumsum(np.hypot(*np.diff(outer, axis=0).T))])
    s /= s[-1]
    bore_z = (z_n - e) + s * ((z_bk + e) - (z_n - e))
    return outer, tho, s, bore_z


def _sheet_lines(outer, tho, bore_z, r_b, e, phi, n=48):
    """отрезки закрученной поверхности: от точки втулки (с запасом e наружу) к отверстию (с запасом e внутрь)"""
    s_all = np.linspace(0, len(outer) - 1, n).round().astype(int)
    s_all = np.unique(s_all)
    lines = []
    for i in s_all:
        z0, r0 = outer[i]
        z1, r1 = bore_z[i], r_b
        v = np.array([z0 - z1, r0 - r1])
        L = np.hypot(*v)
        u = v / L
        a = np.array([z0, r0]) + e * u
        b = np.array([z1, r1]) - min(e, 0.5 * r_b / max(abs(u[1]), 1e-6)) * u      # внутрь отверстия, но не через ось
        th = tho[i] + phi
        P = lambda zr: np.array([zr[1] * math.cos(th), zr[1] * math.sin(th), zr[0]])
        lines.append((P(a), P(b)))
    return lines


def _hub_polygon(imp, z_n, z_bk, r_b, dbl):
    """меридиональный многоугольник ступицы (z, r) для классификации элементов; для двустороннего - обе стороны"""
    d = imp.d1
    zr = imp.hub.at(np.linspace(0, 1, 200))
    poly = [[z_n, r_b], [z_n, d["r1h"]]] + zr.tolist() + [[z_bk, d["r2"]], [z_bk, r_b]]
    poly = np.asarray(poly)
    if dbl:
        mir = poly.copy()
        mir[:, 0] = 2 * z_bk - mir[:, 0]
        poly = np.vstack([poly[:-1], mir[::-1][1:]])
    return poly


def _inside(poly, z, r):
    from matplotlib.path import Path
    return Path(poly).contains_points(np.stack([z, r], 1), radius=1e-9)


def wheel_sector(spec, level="base", threads=1, keep_step=None):
    """сектор колеса заодно с лопатками: impeller.py -> ступица-сектор (закрученные границы) + лопатки -> C3D10.
    spec: parts3d.wheel_spec (D2, b2, D1s, nu, beta2b, Z, splitters, double_sided, bore, ...)"""
    IM = impeller_module()
    gmsh = _gmsh()
    t0 = time.time()
    imp = make_impeller(spec)
    d, p = imp.d1, imp.p
    Z = int(p["Z"])
    pitch = 2 * math.pi / Z
    dbl = bool(p["double_sided"])
    if dbl:
        t_back = 0.5 * max(p["t_back"] * d["D2"], 2 * 1.5 * 0.25 * d["b2"])
    else:
        t_back = max(p["t_back"] * d["D2"], 1.5 * 0.25 * d["b2"])
    z_n, z_bk = -0.3 * d["r1h"], d["Lz"] + t_back
    solid = spec.get("bore") == "core"
    r_b = (A["core_bore_r1h"] if solid else A["tie_bore_r1h"]) * d["r1h"]
    e = 0.05 * d["D2"]
    phi1 = -0.25 * pitch if p["splitters"] else -0.5 * pitch
    outer, tho, s, bore_z = _hub_chain(imp, z_n, z_bk, r_b, e)
    lines1 = _sheet_lines(outer, tho, bore_z, r_b, e, phi1)
    P1 = np.vstack([np.vstack(x) for x in lines1])
    th = np.unwrap(np.arctan2(P1[:, 1], P1[:, 0]))
    lo, hi = th.min() - 0.3, th.max() + pitch + 0.3
    if hi - lo >= 2 * math.pi - 0.05:
        raise RuntimeError("сектор ступицы охватывает больше полного оборота")
    t_m = 0.5 * (p["t_hub"] + p["t_tip"]) * d["D2"]
    lv = A["levels"][level]
    h_b = lv["k_blade"] * t_m
    h_hub = lv["k_hub"] * h_b
    _init(gmsh, threads)
    try:
        occ = gmsh.model.occ
        gmsh.model.add("wheel")
        mp = IM._merid_pt
        P0, P1_ = mp(occ, z_n, r_b), mp(occ, z_n, d["r1h"])
        hb = IM._bezier_curve(occ, imp.Ph)
        occ.synchronize()
        bb = gmsh.model.getBoundary([(1, hb)], oriented=False)
        pa, pb = bb[0][1], bb[1][1]
        P4, P5 = mp(occ, z_bk, d["r2"]), mp(occ, z_bk, r_b)
        cs = [occ.addLine(P0, P1_), occ.addLine(P1_, pa), hb, occ.addLine(pb, P4), occ.addLine(P4, P5), occ.addLine(P5, P0)]
        body = IM._revolve_face(occ, cs, lo, hi - lo)
        occ.synchronize()
        V_body = occ.getMass(3, body)
        s1 = IM._loft_sheet(occ, lines1)
        s2 = occ.copy([(2, s1)])[0][1]
        occ.rotate([(2, s2)], 0, 0, 0, 0, 0, 1, pitch)
        out, omap = occ.fragment([(3, body)], [(2, s1), (2, s2)])
        occ.synchronize()
        vols = [t for dd, t in out if dd == 3]
        m1 = {t for dd, t in omap[1]}
        m2 = {t for dd, t in omap[2]}
        both = [v for v in vols if ({abs(t) for _, t in gmsh.model.getBoundary([(3, v)], oriented=False)} & m1)
                and ({abs(t) for _, t in gmsh.model.getBoundary([(3, v)], oriented=False)} & m2)]
        if len(vols) != 3 or len(both) != 1:
            raise RuntimeError(f"закрученные границы не делят ступицу на 3 объёма ({len(vols)}, средних {len(both)})")
        hub = both[0]
        V_hub = occ.getMass(3, hub)
        V_rel = V_hub / (V_body * pitch / (hi - lo)) - 1.0
        if abs(V_rel) > 2e-3:
            raise RuntimeError(f"объём сектора ступицы не равен 1/Z тела вращения ({V_rel:.2e})")
        rest = [(3, v) for v in vols if v != hub] + [(2, t) for dd, t in out if dd == 2]
        occ.remove(rest, recursive=False)        # без рекурсии: грани разреза общие с ступицей
        occ.synchronize()
        blades = _blade_solids(occ, imp)
        V_bl = sum(occ.getMass(3, b) for b in blades)
        parts = [(3, hub)] + [(3, b) for b in blades]
        if dbl:
            mir = occ.copy(parts)
            occ.mirror(mir, 0, 0, 1, -z_bk)
            parts += mir
        res, _ = occ.fuse([parts[0]], parts[1:])
        occ.synchronize()
        vv = [t for dd, t in res if dd == 3]
        if len(vv) != 1:
            raise RuntimeError(f"после сплавления ступицы и лопаток {len(vv)} тел вместо 1")
        vol = vv[0]
        IM._keep_only(gmsh, [vol])
        V_sector = occ.getMass(3, vol)
        # опорные поверхности для классификации граней
        r1 = IM._loft_sheet(occ, lines1)
        r2 = occ.copy([(2, r1)])[0][1]
        occ.rotate([(2, r2)], 0, 0, 0, 0, 0, 1, pitch)
        refm, refs = [r1], [r2]
        if dbl:
            rm = occ.copy([(2, r1), (2, r2)])
            occ.mirror(rm, 0, 0, 1, -z_bk)
            refm.append(rm[0][1]); refs.append(rm[1][1])
        occ.synchronize()
        tol = max(1e-6, 2e-4 * d["b2"])
        master, slave = _classify_periodic(gmsh, vol, refm, refs, tol)
        bore_f = []
        for _, f in gmsh.model.getBoundary([(3, vol)], oriented=False):
            f = abs(f)
            P = _face_points(gmsh, f)
            if np.all(np.abs(np.hypot(P[:, 0], P[:, 1]) - r_b) < tol):
                bore_f.append(f)
        occ.remove([(2, t) for t in refm + refs], recursive=True)
        occ.synchronize()
        if not master or not bore_f:
            raise RuntimeError(f"не найдены грани: периодические {len(master)}/{len(slave)}, отверстие {len(bore_f)}")
        pairs = _pair_faces(gmsh, master, slave, pitch, 10 * tol)
        if keep_step:
            gmsh.write(keep_step)
        t_geo = time.time() - t0
        poly = _hub_polygon(imp, z_n, z_bk, r_b, dbl)
        from matplotlib.path import Path
        path = Path(poly)
        hz = poly

        def size(x, y, z):
            r = math.hypot(x, y)
            if not path.contains_point((z, r), radius=-1e-9):
                return h_b
            # расстояние до наружной границы ступицы
            dmin = float(np.min(np.hypot(hz[:, 0] - z, hz[:, 1] - r)))
            return min(h_hub, h_b * (1 + lv["grow"] * dmin / h_b))
        t1 = time.time()
        q, snap = _mesh_volume(gmsh, vol, pairs, pitch, size, 0.25 * h_b, h_hub, threads)
        nodes, ids, elems, nm, ns, nb = _extract(gmsh, vol, master, slave, bore_f)
        t_mesh = time.time() - t1
    finally:
        gmsh.logger.stop()
        gmsh.finalize()
    per = _match_periodic(nodes, nm, ns, pitch, 1e-7 * d["D2"] * 100)
    if not per["ok"]:
        raise RuntimeError(f"узлы периодических граней не совпадают: {per}")
    cen = nodes[elems[:, :4] - 1].mean(axis=1)
    rc, zc = np.hypot(cen[:, 0], cen[:, 1]), cen[:, 2]
    inside = _inside(poly, zc, rc)
    ecls = np.where(inside, CLASS_HUB, CLASS_BLADE)
    z_mid = z_bk if dbl else 0.5 * (z_n + z_bk)
    nfix = _fix_node(nodes, nb, ns, nm, z_mid)
    sets = dict(NMASTER=nm, NSLAVE=ns, NBORE=[t for t in nb if t not in set(ns)], NFIX=[nfix])
    meta = dict(kind="wheel", Z=Z, pitch_rad=pitch, splitters=bool(p["splitters"]), double_sided=dbl, bore=spec.get("bore"),
                r_bore=r_b, z_nose=z_n, z_back=z_bk, t_back=t_back, D2=d["D2"], r2=d["r2"], r1h=d["r1h"], r1s=d["r1s"], b2=d["b2"],
                Lz=d["Lz"], beta1h=d["beta1h"], beta1s=d["beta1s"], beta2b=d["beta2b"], wrap_hub_deg=imp.geo["wrap_hub_deg"],
                t_hub=p["t_hub"] * d["D2"], t_tip=p["t_tip"] * d["D2"], h_blade=h_b, h_hub=h_hub, level=level,
                volume=dict(V_sector_occ=V_sector, V_hub_sector=V_hub, V_blades_trimmed=V_bl, hub_vs_body_over_Z=V_rel),
                quality=q, periodic=per, periodic_snap=snap, n_nodes=int(len(ids)), n_elems=int(len(elems)),
                n_hub_elems=int((ecls == CLASS_HUB).sum()), hub_polygon=poly.tolist(), warnings=list(imp.warnings),
                t_geometry_s=round(t_geo, 1), t_mesh_s=round(t_mesh, 1),
                note="сектор колеса: ступица - тело вращения impeller.py с закрученными границами по середине между лопатками, "
                     "лопатки impeller.py целиком; галтели у корня лопаток нет - пик Мизеса у корня зависит от сетки (оценка)")
    return Mesh3D(nodes, ids, elems, ecls, sets, meta)


# ----------------------------------------------------------------------------- осевое колесо
def _airfoil(c_ax, beta_in, beta_out, t_max_c, t_te_c, n=36):
    """сечение пера в плоскости (x - ось, y - окружное направление): средняя линия - парабола с углами beta_in на
    входе и -beta_out на выходе (от оси, град), толщина NACA 00xx по нормали; возвращает замкнутый контур (2n-2, 2)"""
    ti, to = math.tan(math.radians(beta_in)), math.tan(math.radians(beta_out))
    a, b = ti, (-to - ti) / (2 * c_ax)
    x = 0.5 * (1 - np.cos(np.linspace(0, math.pi, n))) * c_ax
    y = a * x + b * x * x
    dy = a + 2 * b * x
    L = np.concatenate([[0.0], np.cumsum(np.hypot(np.diff(x), np.diff(y)))])
    chord = L[-1]
    u = L / chord
    t = 5 * t_max_c * chord * (0.2969 * np.sqrt(u) - 0.1260 * u - 0.3516 * u ** 2 + 0.2843 * u ** 3 - 0.1015 * u ** 4)
    t = np.maximum(t, 0.5 * t_te_c * chord * u)          # выходная кромка конечной толщины
    nx, ny = -dy / np.hypot(1, dy), 1 / np.hypot(1, dy)
    up = np.stack([x + t * nx, y + t * ny], 1)
    lo = np.stack([x - t * nx, y - t * ny], 1)
    pts = np.vstack([up[::-1], lo[1:]])                 # от выходной кромки по спинке к входной и по корыту обратно
    xc = 0.5 * (x.min() + x.max())
    return pts - np.array([xc, 0.0]), chord


def _box_wedge(occ, r0, r1, L, a0, da):
    """клин: прямоугольник r0..r1 × z ±L/2 в плоскости x-z, повёрнутый на угол da от a0"""
    p = [occ.addPoint(r0, 0, -L / 2), occ.addPoint(r1, 0, -L / 2), occ.addPoint(r1, 0, L / 2), occ.addPoint(r0, 0, L / 2)]
    f = occ.addPlaneSurface([occ.addCurveLoop([occ.addLine(p[i], p[(i + 1) % 4]) for i in range(4)])])
    if a0:
        occ.rotate([(2, f)], 0, 0, 0, 0, 0, 1, a0)
    return [t for dd, t in occ.revolve([(2, f)], 0, 0, 0, 0, 0, 1, da) if dd == 3][0]


def _platform(occ, r0, t, L, w_ang, k_ch):
    """полка: тонкий слой r0..r0+t, в развёртке (z, θ) - параллелограмм ширины w_ang вдоль хорды пера (θ = z·k_ch/r)"""
    wires = []
    for r in (r0, r0 + t):
        pts = []
        for z, sgn in ((-L / 2, -1), (-L / 2, 1), (L / 2, 1), (L / 2, -1)):
            th = sgn * 0.5 * w_ang + z * k_ch / r
            pts.append(occ.addPoint(r * math.cos(th), r * math.sin(th), z))
        wires.append(occ.addWire([occ.addLine(pts[i], pts[(i + 1) % 4]) for i in range(4)]))
    return [x for dd, x in occ.addThruSections(wires, makeSolid=True, makeRuled=True) if dd == 3][0]


def axial_sector(spec, level="base", threads=1, keep_step=None):
    """сектор осевого колеса: диск (профиль r, t), обод, (замок, полка) и перо. spec - parts3d.axial_spec"""
    gmsh = _gmsh()
    t0 = time.time()
    Z = int(spec["Z"])
    pitch = 2 * math.pi / Z
    blisk = bool(spec["blisk"])
    r_h, r_t = spec["r_hub"], spec["r_tip"]
    bx, chord = spec["bx"], spec["chord"]
    R = spec["R_live"]
    w_rim = max(A["rim_w_bx"] * bx, spec.get("t_rim") or 0.0)
    prof_r, prof_t = np.asarray(spec["disc_r"], float), np.asarray(spec["disc_t"], float)
    t_mid = 0.5 * (A["t_max_c_hub"] + A["t_max_c_tip"]) * chord
    lv = A["levels"][level]
    h_b = lv["k_blade"] * 0.6 * t_mid
    h_hub = lv["k_hub"] * h_b
    t_plat = A["plat_t_c"] * chord
    secs = spec["sections"]
    k_ch = 0.5 * (math.tan(math.radians(secs[0]["beta_in"])) - math.tan(math.radians(secs[0]["beta_out"])))
    _init(gmsh, threads)
    try:
        occ = gmsh.model.occ
        gmsh.model.add("axial")
        # 1) диск и обод: меридиональный контур, симметричный относительно z = 0
        r_top = r_h if blisk else R
        up = [(0.5 * max(t, 1e-4), r) for r, t in zip(prof_r, prof_t) if r < r_top - 1e-9]
        up.append((0.5 * max(prof_t[-1], 1e-4), min(R, r_top)) if prof_r[-1] >= r_top - 1e-9 else (0.5 * prof_t[-1], prof_r[-1]))
        if r_top > up[-1][1] + 1e-9:
            up.append((0.5 * w_rim, up[-1][1]))
            up.append((0.5 * w_rim, r_top))
        contour = up + [(-z, r) for z, r in reversed(up)]
        tags = [occ.addPoint(r, 0.0, z) for z, r in contour]
        face = occ.addPlaneSurface([occ.addCurveLoop([occ.addLine(tags[i], tags[(i + 1) % len(tags)]) for i in range(len(tags))])])
        occ.rotate([(2, face)], 0, 0, 0, 0, 0, 1, -0.5 * pitch)
        disc = [t for dd, t in occ.revolve([(2, face)], 0, 0, 0, 0, 0, 1, pitch) if dd == 3][0]
        solids = [(3, disc)]
        r_root = r_h
        if not blisk:
            # 2) замок: хвостовик + ножка лопатки от живого обода; выступ диска между хвостовиками (половины у границ сектора)
            s_rim = R * pitch
            h_neck = A["neck_h_c"] * chord
            r_ft = max(R + 0.3 * (r_h - R), r_h - h_neck)
            w_foot, w_neck = A["foot_w_pitch"] * s_rim, A["neck_w_pitch"] * s_rim
            L_att = A["rim_w_bx"] * bx
            ov = 1e-4 * R
            foot = occ.addBox(R - ov, -0.5 * w_foot, -0.5 * L_att, r_ft - R + ov, w_foot, L_att)
            neck = occ.addBox(r_ft - ov, -0.5 * w_neck, -0.5 * L_att, r_h - r_ft + 2 * ov, w_neck, L_att)
            g_lug = A["lug_gap_pitch"] * s_rim
            y0 = 0.5 * w_foot + g_lug
            lug = occ.addBox(R - ov, y0, -0.5 * L_att, r_ft - R + ov, s_rim - 2 * y0, L_att)
            lug2 = occ.copy([(3, lug)])
            occ.rotate(lug2, 0, 0, 0, 0, 0, 1, -pitch)
            wv = _box_wedge(occ, R - 2 * ov, r_h + 0.1 * (r_t - r_h), 2 * L_att, -0.5 * pitch, pitch)
            lugs, _ = occ.intersect([(3, lug)] + lug2, [(3, wv)], removeTool=True)
            solids += [(3, foot), (3, neck)] + [x for x in lugs if x[0] == 3]
            # 3) полка
            plat = _platform(occ, r_h, t_plat, A["plat_len_bx"] * bx, (1 - A["plat_gap_pitch"]) * pitch, k_ch)
            solids.append((3, plat))
            r_root = r_h + t_plat
        # 4) перо: сечения (свободный вихрь), центр сечения - на оси x
        wires = []
        for i, sec in enumerate(secs):
            rr = (r_root - 0.3 * t_plat if not blisk else r_root - 0.02 * chord) if i == 0 else (r_t if i == len(secs) - 1 else sec["r"])
            pts, _ = _airfoil(sec["c_ax"], sec["beta_in"], sec["beta_out"], sec["t_max_c"], A["t_te_c"])
            pts = pts - pts.mean(axis=0)
            ptags = [occ.addPoint(rr * math.cos(yy / rr), rr * math.sin(yy / rr), xx) for xx, yy in pts]
            wires.append(occ.addWire([occ.addSpline(ptags + [ptags[0]])]))
        foil = [t for dd, t in occ.addThruSections(wires, makeSolid=True, makeRuled=False) if dd == 3]
        if len(foil) != 1:
            raise RuntimeError("перо не построено (addThruSections)")
        solids.append((3, foil[0]))
        occ.synchronize()
        V_foil = occ.getMass(3, foil[0])
        res, _ = occ.fuse([solids[0]], solids[1:])
        occ.synchronize()
        vv = [t for dd, t in res if dd == 3]
        if len(vv) != 1:
            raise RuntimeError(f"сектор осевого колеса: {len(vv)} тел после сплавления вместо 1")
        vol = vv[0]
        impeller_module()._keep_only(gmsh, [vol])
        V_sector = occ.getMass(3, vol)
        tol = max(1e-6, 1e-4 * r_t)
        master, slave, bore_f = [], [], []
        for _, f in gmsh.model.getBoundary([(3, vol)], oriented=False):
            f = abs(f)
            P = _face_points(gmsh, f)
            th = np.arctan2(P[:, 1], P[:, 0])
            r = np.hypot(P[:, 0], P[:, 1])
            if np.all(np.abs(th + 0.5 * pitch) * r < tol):
                master.append(f)
            elif np.all(np.abs(th - 0.5 * pitch) * r < tol):
                slave.append(f)
            elif np.all(np.abs(r - prof_r[0]) < tol):
                bore_f.append(f)
        if not master or not bore_f:
            raise RuntimeError(f"не найдены грани: периодические {len(master)}/{len(slave)}, отверстие {len(bore_f)}")
        pairs = _pair_faces(gmsh, master, slave, pitch, 10 * tol)
        if keep_step:
            gmsh.write(keep_step)
        t_geo = time.time() - t0

        def size(x, y, z):
            r = math.hypot(x, y)
            if r >= r_root - 0.5 * t_plat:
                return h_b
            if r >= R:
                return 1.5 * h_b
            return min(h_hub, 1.5 * h_b * (1 + lv["grow"] * (R - r) / (1.5 * h_b)))
        t1 = time.time()
        q, snap = _mesh_volume(gmsh, vol, pairs, pitch, size, 0.25 * h_b, h_hub, threads)
        nodes, ids, elems, nm, ns, nb = _extract(gmsh, vol, master, slave, bore_f)
        t_mesh = time.time() - t1
    finally:
        gmsh.logger.stop()
        gmsh.finalize()
    per = _match_periodic(nodes, nm, ns, pitch, 1e-5 * r_t)
    if not per["ok"]:
        raise RuntimeError(f"узлы периодических граней не совпадают: {per}")
    cen = nodes[elems[:, :4] - 1].mean(axis=1)
    rc = np.hypot(cen[:, 0], cen[:, 1])
    r_disc = r_h if blisk else R
    ecls = np.where(rc < r_disc, CLASS_HUB, np.where(rc < r_root, CLASS_ATT, CLASS_BLADE))
    nfix = _fix_node(nodes, nb, ns, nm, 0.0)
    sets = dict(NMASTER=nm, NSLAVE=ns, NBORE=[t for t in nb if t not in set(ns)], NFIX=[nfix])
    meta = dict(kind="axial", Z=Z, pitch_rad=pitch, blisk=blisk, r_bore=float(prof_r[0]), R_live=R, r_hub=r_h, r_root=r_root, r_tip=r_t,
                w_rim=w_rim, chord=chord, bx=bx, t_plat=t_plat, h_blade=h_b, h_hub=h_hub, level=level,
                volume=dict(V_sector_occ=V_sector, V_airfoil=V_foil), quality=q, periodic=per, periodic_snap=snap,
                n_nodes=int(len(ids)), n_elems=int(len(elems)), n_hub_elems=int((ecls == CLASS_HUB).sum()),
                n_att_elems=int((ecls == CLASS_ATT).sum()), t_geometry_s=round(t_geo, 1), t_mesh_s=round(t_mesh, 1),
                note=("блиск: перо на ободе" if blisk else
                      "вставная лопатка: хвостовик и ножка связаны с живым ободом, выступ диска отдельно (паз без окружной связи); "
                      "контакт и смятие зубьев не моделируются") + "; перо - дуговая средняя линия, NACA 00xx (оценка)")
    return Mesh3D(nodes, ids, elems, ecls, sets, meta)


# ----------------------------------------------------------------------------- самопроверка: кольцо и балка
def ring_sector(a, b, thk, Z, level_h, threads=1):
    """кольцо a..b толщиной thk, сектор 360/Z, периодическая сетка C3D10; NBORE - узлы отверстия"""
    gmsh = _gmsh()
    pitch = 2 * math.pi / Z
    _init(gmsh, threads)
    try:
        occ = gmsh.model.occ
        o = occ.addCylinder(0, 0, -thk / 2, 0, 0, thk, b, angle=pitch)
        i = occ.addCylinder(0, 0, -thk / 2, 0, 0, thk, a, angle=pitch)
        vol, _ = occ.cut([(3, o)], [(3, i)])
        occ.synchronize()
        v = vol[0][1]
        faces = [abs(f) for _, f in gmsh.model.getBoundary(vol, oriented=False)]
        cen = {f: occ.getCenterOfMass(2, f) for f in faces}
        fm = [f for f in faces if abs(cen[f][1]) < 1e-12 and abs(cen[f][2]) < 1e-9]
        fs = [f for f in faces if abs(math.atan2(cen[f][1], cen[f][0]) - pitch) < 1e-6]
        fb = [f for f in faces if abs(math.hypot(cen[f][0], cen[f][1]) - a * math.sin(pitch / 2) / (pitch / 2)) < 1e-6]
        pairs = [(fs[0], fm[0])]
        q, snap = _mesh_volume(gmsh, v, pairs, pitch, lambda x, y, z: level_h, 0.25 * level_h, level_h, threads)
        nodes, ids, elems, nm, ns, nb = _extract(gmsh, v, fm, fs, fb)
    finally:
        gmsh.logger.stop()
        gmsh.finalize()
    per = _match_periodic(nodes, nm, ns, pitch, 1e-9)
    nfix = _fix_node(nodes, nb, ns, nm, 0.0)
    ecls = np.full(len(elems), CLASS_HUB)
    return Mesh3D(nodes, ids, elems, ecls, dict(NMASTER=nm, NSLAVE=ns, NBORE=[t for t in nb if t not in set(ns)], NFIX=[nfix]),
                  dict(kind="ring", Z=Z, pitch_rad=pitch, r_bore=a, quality=q, periodic=per, n_nodes=int(len(ids)), n_elems=int(len(elems))))


def beam(L, w, h, n_h, threads=1):
    """консольная балка L × w × h вдоль x (заделка x = 0), C3D10; NROOT - узлы заделки"""
    gmsh = _gmsh()
    _init(gmsh, threads)
    try:
        occ = gmsh.model.occ
        v = occ.addBox(0, -w / 2, -h / 2, L, w, h)
        occ.synchronize()
        hs = h / n_h
        gmsh.option.setNumber("Mesh.MeshSizeMin", 0.5 * hs)
        gmsh.option.setNumber("Mesh.MeshSizeMax", hs)
        gmsh.option.setNumber("Mesh.ElementOrder", 2)
        gmsh.model.addPhysicalGroup(3, [v], name="EVOL")
        gmsh.model.mesh.generate(3)
        root = [abs(f) for _, f in gmsh.model.getBoundary([(3, v)], oriented=False)
                if abs(occ.getCenterOfMass(2, abs(f))[0]) < 1e-12]
        nodes, ids, elems, nr, _, _ = _extract(gmsh, v, root, [], [])
    finally:
        gmsh.logger.stop()
        gmsh.finalize()
    return Mesh3D(nodes, ids, elems, np.full(len(elems), CLASS_BLADE), dict(NROOT=nr), dict(kind="beam", n_nodes=int(len(ids))))
