"""Воронка камеры сгорания KV, уровни 2-3: генератор трёхмерных секторов OpenFOAM по варианту.

Уровень 2 (README, «Уровень 2») - два сектора на вариант, сетка blockMesh с дугами (ось z, угол θ, радиус r):
  * «элемент 3D» (element3d): один смесительный элемент, сектор 360/k_je с одной топливной струёй (циклические
    границы при крутке; при k_je = 1 или без крутки - половина с плоскостями симметрии). Струя - квадратное отверстие
    той же площади на стенке канала, канал, срез, зона горения - как ячейка уровня 1. RANS LTS (reactingFoam, EDC,
    GRI-Mech 3.0 + TDAC, P1). Цепочки рестартов: запас по срыву - топливо ×1,0; 0,9; 0,8 ... до погасания
    (полнота < 0,9 или T_max < T3 + 300 K): запас = phi_расч / phi_погасания; запас по проскоку - доля H2 в топливе
    0; 0,3; 0,5; 0,7; 1,0 по объёму до пламени в канале. Распыл дизеля - тот же сектор, sprayFoam (словари - spray_dicts);
  * «труба 3D с теплообменом стенки» (liner3d_cht): сектор трубы 360/n (одно отверстие разбавления на сектор), газ
    и стенка (cellZone gas / solid, толщина comb_space.WALL_T); отверстия разбавления и щели охлаждения - ячейки стенки,
    переведённые в газ (topoSet), их наружные грани - входы (createPatch); chtMultiRegionFoam: газ - реагирующий,
    стенка - heSolidThermo; снаружи стенки - теплоотдача к воздуху кожуха (h и T3 - comb_crn.wall_temperature).
Уровень 3 - тот же элемент 3D мельче, LES (WALE) нестационарно (les_dicts): пульсации и передаточная функция пламени
для comb_acoustic.py.

Только стандартная библиотека Python (как comb_cfd). Запуск: python comb_cfd3d.py --inputs l1_inputs.json --row 0
--kind element3d|liner3d_cht --out /tmp/case
"""
import argparse, json, math, os, shutil

import comb_cfd as CF

HERE = os.path.dirname(os.path.abspath(__file__))
MESH3 = {"L2": dict(nr_ch=10, nt=12, nz_ch_max=30, nz_cell=110, nr_out_max=26, nr_l=30, nz_l=90, nt_l=16, nr_w=3),
         "L3": dict(nr_ch=16, nt=20, nz_ch_max=48, nz_cell=150, nr_out_max=34, nr_l=50, nz_l=160, nt_l=28, nr_w=4)}
R_CORE = 0.02          # внутренний радиус сектора, доля наружного (ось не включается: тонкая труба симметрии)


class Blocks3D:
    """многоблочная сетка сектора по сетке отрезков z (zs), r (rs), θ (ts, градусы); mask[i][j][m]; zone[i][j][m]"""

    def __init__(self, zs, rs, ts, mask, nz, nr, nt, gz, gr, zone=None):
        self.zs, self.rs, self.ts, self.mask = zs, rs, ts, mask
        self.nz, self.nr, self.nt, self.gz, self.gr = nz, nr, nt, gz, gr
        self.zone = zone
        self.vid, self.verts = {}, []
        for k, z in enumerate(zs):
            for l, r in enumerate(rs):
                for q, t in enumerate(ts):
                    a = math.radians(t)
                    self.vid[(k, l, q)] = len(self.verts)
                    self.verts.append((r * math.cos(a), r * math.sin(a), z))

    def has(self, i, j, m):
        return (0 <= i < len(self.zs) - 1 and 0 <= j < len(self.rs) - 1 and 0 <= m < len(self.ts) - 1
                and self.mask[i][j][m])

    def n_cells(self):
        return sum(self.nz[i] * self.nr[j] * self.nt[m] for i in range(len(self.zs) - 1)
                   for j in range(len(self.rs) - 1) for m in range(len(self.ts) - 1) if self.mask[i][j][m])

    def dict_text(self, classify, ptypes, cyclic=None):
        v = self.vid
        H = ["scale 1;\n\nvertices\n("] + [f"    ({p[0]:.9g} {p[1]:.9g} {p[2]:.9g})" for p in self.verts]
        H.append(");\n\nblocks\n(")
        faces, edges = {}, set()
        for i in range(len(self.zs) - 1):
            for j in range(len(self.rs) - 1):
                for m in range(len(self.ts) - 1):
                    if not self.mask[i][j][m]:
                        continue
                    h = [v[(i, j, m)], v[(i, j + 1, m)], v[(i, j + 1, m + 1)], v[(i, j, m + 1)],
                         v[(i + 1, j, m)], v[(i + 1, j + 1, m)], v[(i + 1, j + 1, m + 1)], v[(i + 1, j, m + 1)]]
                    zn = (self.zone[i][j][m] + " ") if self.zone and self.zone[i][j][m] else ""
                    H.append(f"    hex ({' '.join(map(str, h))}) {zn}({self.nr[j]} {self.nt[m]} {self.nz[i]}) "
                             f"simpleGrading ({self.gr[j]:.6g} 1 {self.gz[i]:.6g})")
                    for k in (i, i + 1):
                        for l in (j, j + 1):
                            edges.add((k, l, m))
                    sides = dict(rin=((h[0], h[4], h[7], h[3]), (i, j - 1, m)), rout=((h[1], h[2], h[6], h[5]), (i, j + 1, m)),
                                 tlo=((h[0], h[1], h[5], h[4]), (i, j, m - 1)), thi=((h[3], h[7], h[6], h[2]), (i, j, m + 1)),
                                 zlo=((h[0], h[3], h[2], h[1]), (i - 1, j, m)), zhi=((h[4], h[5], h[6], h[7]), (i + 1, j, m)))
                    for side, (f, nb) in sides.items():
                        if self.has(*nb):
                            continue
                        name = classify(i, j, m, side)
                        faces.setdefault(name, []).append(f)
        H.append(");\n\nedges\n(")
        for k, l, m in sorted(edges):
            r, z = self.rs[l], self.zs[k]
            a = math.radians(0.5 * (self.ts[m] + self.ts[m + 1]))
            H.append(f"    arc {v[(k, l, m)]} {v[(k, l, m + 1)]} ({r * math.cos(a):.9g} {r * math.sin(a):.9g} {z:.9g})")
        H.append(");\n\nboundary\n(")
        for name, fl in faces.items():
            t = ptypes.get(name, "patch")
            extra = ""
            if cyclic and name in cyclic:
                t = "cyclic"
                extra = (f"        neighbourPatch {cyclic[name]};\n        transform rotational;\n"
                         f"        rotationAxis (0 0 1);\n        rotationCentre (0 0 0);\n")
            H.append(f"    {name}\n    {{\n        type {t};\n{extra}        faces\n        (")
            H += [f"            ({' '.join(map(str, f))})" for f in fl]
            H.append("        );\n    }")
        H.append(");\n\nmergePatchPairs\n(\n);\n")
        return "\n".join(H), sorted(faces)


def _seg_cells(L, d, nmin, nmax):
    return int(max(nmin, min(nmax, math.ceil(L / max(d, 1e-9)))))


def build_element3d_case(e, li, out, level="L2", of="2512", iters=None, fuel_scale=1.0, h2_vol=0.0, length_frac=1.0):
    """элемент 3D уровня 2 (или 3 - level L3); fuel_scale и h2_vol - шаги цепочек запаса по срыву и проскоку"""
    it = dict(CF.ITERS, **(iters or {}))
    M = MESH3[level]
    c, d, L = e["concept"], e["derived"], e["loads"][li]
    g = CF.cell_geometry(e, li)
    p3, T3, rho3 = L["p3_Pa"], L["T3_K"], L["rho3"]
    r_ch, R_c, l_ch, L_c = g["r_ch"], g["R_c"], g["l_ch"], g["L_cell"] * length_frac
    kje = max(1, int(c["k_je"]))
    swirl = c["swirl"]
    if kje >= 2 and swirl > 0.05:
        theta, mode = 360.0 / kje, "cyclic"
        t_lo, t_hi = -theta / 2, theta / 2
    else:
        theta, mode = 180.0 / kje, "symmetry"
        t_lo, t_hi = 0.0, theta
    m_fuel_jet = g["m_fuel_main"] * fuel_scale / kje
    d_j = 1e-3 * d["d_jet_mm"]
    a_side = d_j * math.sqrt(math.pi) / 2                # квадрат той же площади, что круглое отверстие
    th_j = math.degrees(a_side / r_ch)
    if mode == "cyclic":
        ts = [t_lo, -th_j / 2, th_j / 2, t_hi]
        jm = 1
    else:
        ts = [0.0, th_j / 2, t_hi]
        jm = 0
    ts = [t for t in ts if t_lo - 1e-12 <= t <= t_hi + 1e-12]
    r_core = R_CORE * R_c
    rs = [r_core, r_ch, R_c]
    z0 = min(0.15 * l_ch, 2 * r_ch)
    zs = sorted(set([0.0, z0, z0 + a_side, l_ch, l_ch + L_c]))
    dr = (r_ch - r_core) / M["nr_ch"]
    nr = [M["nr_ch"], _seg_cells(R_c - r_ch, 1.6 * dr, 4, M["nr_out_max"])]
    gr = [1.0, min(8.0, max(1.0, (R_c - r_ch) / (nr[1] * dr)))]
    nt = []
    for m in range(len(ts) - 1):
        w = math.radians(ts[m + 1] - ts[m]) * r_ch
        nt.append(_seg_cells(w, a_side / 3 if m == jm else 1.5 * dr, 2, M["nt"]))
    nz, gz = [], []
    for i in range(len(zs) - 1):
        a, b = zs[i], zs[i + 1]
        if b <= l_ch + 1e-12:
            nz.append(_seg_cells(b - a, 1.5 * dr if abs(a - z0) > 1e-12 else a_side / 3, 2, M["nz_ch_max"])); gz.append(1.0)
        else:
            nz.append(M["nz_cell"]); gz.append(min(30.0, max(1.0, (b - a) / M["nz_cell"] / (1.5 * dr))))
    mask = [[[(zs[i + 1] > l_ch + 1e-12) or j == 0 for m in range(len(ts) - 1)] for j in range(2)] for i in range(len(zs) - 1)]
    outer_wall = g["n_el_can"] <= 1.0 + 1e-9

    def classify(i, j, m, side):
        z0_, z1_ = zs[i], zs[i + 1]
        if side == "zlo" and i == 0:
            return "inlet_air"
        if side == "zhi" and i == len(zs) - 2:
            return "outlet"
        if side == "zlo":
            return "wall_plate"
        if side == "rin":
            return "core"
        if side == "tlo" and m == 0:
            return "cyc_lo" if mode == "cyclic" else "sym_lo"
        if side == "thi" and m == len(ts) - 2:
            return "cyc_hi" if mode == "cyclic" else "sym_hi"
        if side == "rout" and j == 0 and z1_ <= l_ch + 1e-12:
            return "inlet_fuel" if (abs(z0_ - z0) < 1e-12 and m == jm) else "wall_channel"
        if side == "rout" and j == 1:
            return "wall_liner" if outer_wall else "side_sym"
        return "wall_channel"

    ptypes = dict(wall_plate="wall", wall_channel="wall", wall_liner="wall", side_sym="symmetry", core="symmetry",
                  sym_lo="symmetryPlane", sym_hi="symmetryPlane")
    B = Blocks3D(zs, rs, ts, mask, nz, nr, nt, gz, gr)
    case = os.path.abspath(out)
    if os.path.exists(case):
        shutil.rmtree(case)
    os.makedirs(os.path.join(case, "system"))
    cyc = {"cyc_lo": "cyc_hi", "cyc_hi": "cyc_lo"} if mode == "cyclic" else None
    text, patches = B.dict_text(classify, ptypes, cyc)
    CF.foam_file(os.path.join(case, "system", "blockMeshDict"), "dictionary", "blockMeshDict", text)
    # входы: воздух равномерно с круткой (по окружности - вектор в декартовых осях на грани: средний угол сектора)
    u = g["m_air_main"] / (rho3 * math.pi * (r_ch ** 2 - r_core ** 2))
    tm = math.radians(0.5 * (ts[0] + ts[-1]))
    w = 2 * swirl * u * 0.5                                # крутка на половине радиуса (равномерная по сечению, оценка)
    Yf = {"CH4": 1.0 - h2_vol * 2.016 / (h2_vol * 2.016 + (1 - h2_vol) * 16.043)}
    Yf["H2"] = 1.0 - Yf["CH4"]
    MWf = 1.0 / (Yf["CH4"] / 16.043 + Yf["H2"] / 2.016)
    rho_f = p3 * MWf * 1e-3 / (CF.R_U * 288.15)
    A_hole = a_side ** 2 * (1.0 if mode == "cyclic" else 0.5)
    u_f = (m_fuel_jet * (1.0 if mode == "cyclic" else 0.5)) / (rho_f * A_hole)
    am = math.radians(ts[jm] + 0.5 * (ts[jm + 1] - ts[jm]))
    inlets = {"inlet_air": dict(U=(-w * math.sin(tm), w * math.cos(tm), u), T=T3, Y=CF.air_Y(), kind="air"),
              "inlet_fuel": dict(U=(-u_f * math.cos(am), -u_f * math.sin(am), 0.0), T=288.15, Y=Yf, kind="fuel")}
    special = {p: {"type": "cyclic"} for p in ("cyc_lo", "cyc_hi")} if mode == "cyclic" else \
        {p: {"type": "symmetryPlane"} for p in ("sym_lo", "sym_hi")}
    special["core"] = {"type": "symmetry"}
    if "side_sym" in patches:
        special["side_sym"] = {"type": "symmetry"}
    CF._write_common(case, e, L, inlets, patches, of, it, hot_box=(0.0, l_ch, R_c * 1.5, l_ch + L_c), Lscale=r_ch,
                     special=special)
    _fix_hot_box_3d(case, l_ch, l_ch + L_c, R_c)
    desc = dict(kind="element3d", level=level, load=L["load"], mode=mode, theta_deg=theta, n_cells=B.n_cells(),
                patches=patches, zs=zs, rs=rs, ts=ts, z_dump=l_ch, z_end=zs[-1], r_cell=R_c, p3=p3, T3=T3,
                fuel_scale=fuel_scale, h2_vol=h2_vol, jet=dict(side=a_side, u=u_f, rho=rho_f, m=m_fuel_jet),
                u_ch=u, length_frac=length_frac,
                geometry=dict(r_ch=r_ch, R_c=R_c, l_ch=l_ch, L_cell=L_c, n_el_can=g["n_el_can"]),
                swirl_note="крутка задана равномерной окружной скоростью 2·S·u·0,5 (оценка); при k_je = 1 крутка не "
                           "учитывается (половина элемента с плоскостями симметрии)" if swirl > 0.05 and mode != "cyclic" else "")
    json.dump(desc, open(os.path.join(case, "case.json"), "w", encoding="utf-8"), ensure_ascii=False, indent=1)
    return desc


def _fix_hot_box_3d(case, z0, z1, R):
    p = os.path.join(case, "system", "setFieldsDict")
    if os.path.exists(p):
        t = open(p).read()
        import re
        t = re.sub(r"box \([^)]*\) \([^)]*\);", f"box ({-2 * R:.9g} {-2 * R:.9g} {z0:.9g}) ({2 * R:.9g} {2 * R:.9g} {z1 + 1:.9g});", t)
        open(p, "w").write(t)


def build_liner3d_cht_case(e, cell_out, out, level="L2", of="2512", iters=None):
    """сектор трубы с теплообменом стенки (chtMultiRegionFoam): газ + стенка, отверстие разбавления и щели охлаждения"""
    it = dict(CF.ITERS, **(iters or {}))
    M = MESH3[level]
    c, d, L = e["concept"], e["derived"], e["loads"][0]
    p3, T3, W, rho3 = L["p3_Pa"], L["T3_K"], L["W_air"], L["rho3"]
    n_can = e["n_can"]
    annular = d.get("H_L") is not None
    mat = "керамика" if c["cooling"] == "керамика" else "металл"
    tw = {"металл": 0.0015, "керамика": 0.003}[mat]
    if annular:
        r_lo, r_hi = 0.5 * d["D_m"] - 0.5 * d["H_L"], 0.5 * d["D_m"] + 0.5 * d["H_L"]
    else:
        r_lo, r_hi = R_CORE * 0.5 * d["D_L"], 0.5 * d["D_L"]
    L_L = d["L_L"]
    u_j = 0.70 * math.sqrt(2 * c["dp_inj"] * p3 / rho3)
    m_dil = max(0.0, d["a_dil"] * W / n_can)
    A_dil = m_dil / (rho3 * u_j)
    n_h = int(max(6, min(48, round(math.pi ** 3 * r_hi ** 2 / (9 * max(A_dil, 1e-8))))))
    d_h = math.sqrt(4 * A_dil / (math.pi * n_h))
    theta = 360.0 / n_h
    a_side = d_h * math.sqrt(math.pi) / 2
    th_h = math.degrees(a_side / r_hi)
    z_dil = min(0.6 * L_L, CF.cell_geometry(e, 0)["L_cell"] + 2 * a_side)
    m_cool = d["a_cool"] * W / n_can
    n_rows = max(2, min(10, int(L_L / 0.04)))
    w_row = max(2e-4, m_cool / n_rows / (rho3 * u_j * 2 * math.pi * r_hi))
    z_rows = [L_L * (k + 0.5) / n_rows for k in range(n_rows)]
    z_rows = [z for z in z_rows if abs(z - z_dil) > 2 * a_side]
    zs = sorted(set([0.0, z_dil - a_side / 2, z_dil + a_side / 2, L_L] + sum(([z, z + w_row] for z in z_rows), [])))
    ts = [-theta / 2, -th_h / 2, th_h / 2, theta / 2]
    rs = [r_lo, r_hi, r_hi + tw]
    dz = (r_hi - r_lo) / M["nr_l"]
    nr = [M["nr_l"], M["nr_w"]]
    nt = [_seg_cells(math.radians(ts[1] - ts[0]) * r_hi, 1.5 * dz, 2, M["nt_l"]),
          _seg_cells(a_side, a_side / 4, 3, M["nt_l"]),
          _seg_cells(math.radians(ts[3] - ts[2]) * r_hi, 1.5 * dz, 2, M["nt_l"])]
    nz = [_seg_cells(zs[i + 1] - zs[i], 1.2 * dz, 2, M["nz_l"]) for i in range(len(zs) - 1)]
    mask = [[[True] * 3 for _ in range(2)] for _ in range(len(zs) - 1)]
    zone = [[["gas"] * 3, ["solid"] * 3] for _ in range(len(zs) - 1)]

    def classify(i, j, m, side):
        if side == "zlo":
            return "inlet_core" if j == 0 else "solid_end"
        if side == "zhi":
            return "outlet" if j == 0 else "solid_end"
        if side == "rin":
            return "wall_in" if annular else "core"
        if side == "rout":
            return "outer"
        return "cyc_lo" if side == "tlo" else "cyc_hi"

    B = Blocks3D(zs, rs, ts, mask, nz, nr, nt, [1.0] * len(nz), [1.0, 1.0], zone)
    case = os.path.abspath(out)
    if os.path.exists(case):
        shutil.rmtree(case)
    os.makedirs(os.path.join(case, "system"))
    text, patches = B.dict_text(classify, {"wall_in": "wall", "outer": "wall", "solid_end": "wall", "core": "symmetry"},
                                {"cyc_lo": "cyc_hi", "cyc_hi": "cyc_lo"})
    CF.foam_file(os.path.join(case, "system", "blockMeshDict"), "dictionary", "blockMeshDict", text)
    # topoSet: ячейки стенки в отверстии разбавления и в щелях охлаждения -> газ; затем грани «outer» газа -> входы
    big = 2 * (r_hi + tw)
    rows = "\n".join(f"        {{ name holes; type cellSet; action add; source boxToCell; box ({-big:.9g} {-big:.9g} {z:.9g}) "
                     f"({big:.9g} {big:.9g} {z + w_row:.9g}); }}" for z in z_rows)
    hole_box = (f"({r_hi - 1e-9:.9g} {-a_side / 2:.9g} {z_dil - a_side / 2:.9g}) "
                f"({r_hi + 2 * tw:.9g} {a_side / 2:.9g} {z_dil + a_side / 2:.9g})")
    CF.foam_file(os.path.join(case, "system", "topoSetDict"), "dictionary", "topoSetDict",
                 "actions\n(\n"
                 f"    {{ name holes; type cellSet; action new; source boxToCell; box {hole_box}; }}\n{rows}\n"
                 "    { name sol; type cellSet; action new; source zoneToCell; zone solid; }\n"
                 "    { name holes; type cellSet; action subset; source cellToCell; set sol; }\n"
                 "    { name solid; type cellZoneSet; action subtract; source setToCellZone; set holes; }\n"
                 "    { name gas; type cellZoneSet; action add; source setToCellZone; set holes; }\n);")
    CF.foam_file(os.path.join(case, "system", "createPatchDict.gas"), "dictionary", "createPatchDict",
                 "pointSync false;\n\npatches\n(\n    {\n        name inlet_dil;\n        patchInfo { type patch; }\n"
                 "        constructFrom set;\n        set dilFaces;\n    }\n    {\n        name inlet_cool;\n"
                 "        patchInfo { type patch; }\n        constructFrom set;\n        set coolFaces;\n    }\n);")
    CF.foam_file(os.path.join(case, "system", "topoSetDict.gas"), "dictionary", "topoSetDict",
                 "actions\n(\n    { name outerF; type faceSet; action new; source patchToFace; patch outer; }\n"
                 f"    {{ name dilFaces; type faceSet; action new; source boxToFace; box {hole_box}; }}\n"
                 "    { name dilFaces; type faceSet; action subset; source faceToFace; set outerF; }\n"
                 "    { name coolFaces; type faceSet; action new; source faceToFace; set outerF; }\n"
                 "    { name coolFaces; type faceSet; action subtract; source faceToFace; set dilFaces; }\n);")
    CF.foam_file(os.path.join(case, "system", "controlDict"), "dictionary", "controlDict",
                 f"application     chtMultiRegionFoam;\nstartFrom       latestTime;\nstartTime       0;\n"
                 f"stopAt          endTime;\nendTime         {it['cold'] + it['hot']};\ndeltaT          1;\n"
                 f"writeControl    timeStep;\nwriteInterval   {it['write']};\npurgeWrite      3;\nwriteFormat     ascii;\n"
                 f"writePrecision  9;\nwriteCompression off;\ntimeFormat      general;\ntimePrecision   8;\n"
                 f"runTimeModifiable true;")
    for nm, body in (("fvSchemes", "ddtSchemes { default steadyState; }\ngradSchemes { default Gauss linear; }\n"
                                   "divSchemes { default none; }\nlaplacianSchemes { default Gauss linear corrected; }\n"
                                   "interpolationSchemes { default linear; }\nsnGradSchemes { default corrected; }"),
                     ("fvSolution", "solvers {}\n")):
        CF.foam_file(os.path.join(case, "system", nm), "dictionary", nm, body)
    with open(os.path.join(case, "Allmesh"), "w") as f:
        f.write(ALLMESH_CHT)
    os.chmod(os.path.join(case, "Allmesh"), 0o755)
    desc = dict(kind="liner3d_cht", level=level, n_cells=B.n_cells(), patches=patches, theta_deg=theta, n_holes=n_h,
                d_hole=d_h, hole_side=a_side, z_dil=z_dil, cool_rows=z_rows, w_row=w_row, wall_t=tw, wall_mat=mat,
                annular=annular, r_lo=r_lo, r_hi=r_hi, u_jet=u_j, m_dil=m_dil, m_cool=m_cool, p3=p3, T3=T3,
                T_core_in=cell_out.get("T"), regions=["gas", "solid"],
                solver="chtMultiRegionFoam (v2512/v2606): газ - реагирующий (EDC, GRI-Mech 3.0 + TDAC, P1), стенка - "
                       "heSolidThermo; стык - compressible::turbulentTemperatureRadCoupledMixed; снаружи - "
                       "externalWallHeatFluxTemperature (h, T3 из comb_crn.wall_temperature)",
                todo="словари областей (constant/gas, constant/solid, 0/gas, 0/solid) пишет cht_region_dicts; свойства "
                     "печатного сплава (k(T), cp(T), ρ) - с источником до прогона (data/known_gaps)")
    json.dump(desc, open(os.path.join(case, "case.json"), "w", encoding="utf-8"), ensure_ascii=False, indent=1)
    return desc


ALLMESH_CHT = """#!/bin/bash
# Сетка сектора трубы уровня 2 с областями газ / стенка (comb_cfd3d.py). Окружение OpenFOAM подключено.
cd "${0%/*}" || exit 1
blockMesh > log.blockMesh 2>&1 || exit 1
topoSet > log.topoSet 2>&1 || exit 2
splitMeshRegions -cellZonesOnly -overwrite > log.splitMeshRegions 2>&1 || exit 3
cp system/topoSetDict.gas system/gas/topoSetDict 2>/dev/null || { mkdir -p system/gas; cp system/topoSetDict.gas system/gas/topoSetDict; }
cp system/createPatchDict.gas system/gas/createPatchDict
topoSet -region gas > log.topoSet.gas 2>&1 || exit 4
createPatch -region gas -overwrite > log.createPatch.gas 2>&1 || exit 5
checkMesh -region gas > log.checkMesh.gas 2>&1
checkMesh -region solid > log.checkMesh.solid 2>&1
"""


def spray_dicts(e, li=0):
    """словари распыла дизеля для sprayFoam (v2512/v2606): облако sprayCloud, жидкость C12H26, впрыск по типу подачи
    (центробежная - coneNozzleInjection + LISAAtomization; пневматическая (LDI) - coneNozzleInjection с малым
    углом, распад ReitzDiwakar; испарительная - испарение в канале, liquidEvaporationBoil). Механизм газовой фазы -
    н-додекан Reitz (Cantera nDodecane_Reitz.yaml -> yaml2ck -> chemkinToFoam с тремя правками навыка), TDAC;
    NOx - постобработка сетью реакторов Cantera по зонам поля (в механизме Reitz нет азотной химии)"""
    c, L = e["concept"], e["loads"][li]
    Wfd = L.get("W_fuel_diesel") or 1.125 * L["W_fuel"]
    kind = c["diesel"]
    inj = {"центробежная": dict(model="coneNozzleInjection", atom="LISAAtomization", breakup="ReitzDiwakar",
                                 cone=(30, 80), note="полый конус центробежной форсунки"),
           "пневматическая": dict(model="coneNozzleInjection", atom="none", breakup="ReitzDiwakar", cone=(0, 30),
                                  note="LDI: распыл воздухом смесителя, капли по SMD уровня 0"),
           "испарительная": dict(model="coneNozzleInjection", atom="none", breakup="none", cone=(0, 20),
                                 note="впрыск в канал смешения, испарение до зоны горения")}[kind]
    return dict(solver="sprayFoam", cloud="sprayCloud", liquid="C12H26", injection=inj, m_fuel_kg_s_per_element=
                Wfd / c["n_el"], SMD_m=None, mechanism="nDodecane_Reitz (Cantera), 100 веществ; NOx - постобработка "
                "Cantera GRI-Mech 3.0 по зонам", status="словари - шаблон уровня 2; прогон не проверялся")


def les_dicts():
    """уровень 3: LES сектора элемента (reactingFoam нестационарно, PIMPLE, WALE, PaSR), датчики давления для
    пульсаций и сигнал возбуждения для передаточной функции пламени (comb_acoustic.Network.F)"""
    return dict(solver="reactingFoam (нестационарно, PIMPLE)", turbulence="LES WALE (delta cubeRootVol)",
                combustion="PaSR (Cmix 0,1) или EDC v2005 - тот же, что прошёл PRECCINSTA", chemistry="GRI-Mech 3.0 + TDAC",
                time=dict(dt_rule="Co ≤ 0,5", t_flow_through=12, t_stats=8, note="длительность в временах протекания"),
                probes=["срез канала", "зона пламени (3 точки)", "выход"],
                forcing="широкополосное возбуждение расхода воздуха на входе ±5 % (сумма синусов 50-1500 Гц) -> "
                        "F(ω) = (Q'/Q̄)/(u'/ū) методом Винера - Хопфа / идентификацией КИХ-фильтра",
                status="шаблон уровня 3; прогон не проверялся")


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--inputs", required=True)
    ap.add_argument("--row", type=int, default=0)
    ap.add_argument("--load", type=int, default=0)
    ap.add_argument("--kind", default="element3d", choices=["element3d", "liner3d_cht"])
    ap.add_argument("--level", default="L2", choices=list(MESH3))
    ap.add_argument("--of", default="2512")
    ap.add_argument("--out", required=True)
    a = ap.parse_args(argv)
    e = json.load(open(a.inputs, encoding="utf-8"))["entries"][a.row]
    if a.kind == "element3d":
        dsc = build_element3d_case(e, a.load, a.out, a.level, a.of)
    else:
        h = e["loads"][0]["hot_init"]
        dsc = build_liner3d_cht_case(e, dict(T=h["T_eq"], Y=h["Y"]), a.out, a.level, a.of)
    print(json.dumps({k: dsc[k] for k in ("kind", "n_cells", "patches", "theta_deg")}, ensure_ascii=False))


if __name__ == "__main__":
    main()
