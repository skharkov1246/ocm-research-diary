"""Воронка камеры сгорания KV, уровень 1: генератор осесимметричных кейсов OpenFOAM (реагирующий RANS) по варианту.

Что строится (README, раздел «Уровень 1»). Два осесимметричных клина 5° (блочная сетка blockMesh, ось z, радиус x,
клин по y - как учебник SandiaD_LTS):
  * «ячейка» (cell) - один смесительный элемент в своей доле сечения жаровой трубы: смесительный канал (воздух с
    круткой, топливо - кольцевой щелью на стенке канала или заданным радиальным профилем смеси), срез канала (плита),
    зона пламени и догорания; внешняя граница - симметрия (соседние элементы такие же) или стенка трубы, если элемент
    один на трубу; по схеме - центральный пилот (гидравлическое деление), кольцевая полость у края ячейки (TVC),
    щель воздуха гашения (RQL). Считается на 100 / 50 / 30 % нагрузки;
  * «труба» (liner) - зона разбавления всей жаровой трубы (выносная и трубчатые - круг с осью, кольцевая и
    противоточная - кольцо вокруг оси двигателя): вход - осреднённое состояние на выходе ячейки при 100 % и слой
    воздуха охлаждения у стенок, щели разбавления и охлаждения, выход - поле T для RTDF, стенка - температура
    адиабатной стенки по газу.
Решатель и модели - навык openfoam-combustor (те же, что эталон Sandia D): reactingFoam, псевдовремя LTS,
k-epsilon, EDC v2005, полный GRI-Mech 3.0 (53/325, файлы учебника counterFlowFlame2D_GRI_TDAC openfoam.com
v2512) с TDAC (DAC + ISAT), излучение P1 (radiationProperties учебника SandiaD_LTS v2512). Сначала без химии
(развитие течения), затем горячее поле продуктов в зоне горения (setFields) и химия.

Словари - openfoam.com v2512/v2606 (образ opencfd/openfoam-default:2512); ключ --of 1912 убирает блок functions
(на apt-1912 любой function object падает - навык, раздел 9) - только для местной обкатки.
Только стандартная библиотека Python: генератор идёт и внутри контейнера OpenFOAM на машине AWS.

Запуск (проверка генератора): python comb_cfd.py --inputs l1_inputs.json --row 0 --load 0 --out /tmp/case
"""
import argparse, json, math, os, shutil

HERE = os.path.dirname(os.path.abspath(__file__))
FOAM_RES = os.path.join(HERE, "cfd", "foam")
DATE = "2026-11-14"
WEDGE_HALF_DEG = 2.5
TAN_A = math.tan(math.radians(WEDGE_HALF_DEG))
AIR_X = {"O2": 0.2095, "N2": 0.7808, "AR": 0.0093, "CO2": 0.0004}       # как gtu_cycle.AIR / comb_crn.AIR
MW = {"O2": 31.998, "N2": 28.014, "AR": 39.948, "CO2": 44.009, "CH4": 16.043, "H2": 2.016, "H2O": 18.015,
      "CO": 28.010, "OH": 17.007, "NO": 30.006}
R_U = 8.314462618
MAJOR = ["CH4", "O2", "N2", "AR", "CO2", "H2O", "CO", "H2", "OH", "NO"]  # поля с отдельными файлами; прочие - Ydefault

# Сетка: уровни (навык: не меньше двух сеток; разница итогов - в results-JSON)
MESH = {
    "coarse": dict(nr_ch=10, nz_ch_max=30, nz_cell=110, nr_out_max=30, nr_liner=36, nz_liner=90),
    "medium": dict(nr_ch=16, nz_ch_max=48, nz_cell=180, nr_out_max=48, nr_liner=56, nz_liner=140),
}
ITERS = dict(cold=1500, hot=5000, write=500)  # без химии / с химией / шаг записи (образец SandiaD_LTS: 1500 + 3500)
RING = False   # кольцевая щель топлива (k_je >= 4) выключена: щель в десятки мкм дала локальный шаг LTS ~3e-8 с, и за 3600
               # итераций поток не дошёл до выхода (обкатка 14.11.2026, comb_l1_smoke_results.json). Топливо - профиль смеси
               # с неравномерностью уровня 0; разрешение струй - уровень 2 (элемент 3D)
SLOT_MIN = 5e-4  # наименьшая ширина щелей разбавления, охлаждения, гашения, м (тот же расход - меньше скорость)


def air_Y():
    m = {k: v * MW[k] for k, v in AIR_X.items()}
    s = sum(m.values())
    return {k: v / s for k, v in m.items()}


def mix_Y(phi, far_st, fuel="CH4"):
    """свежая смесь метан - воздух с коэффициентом phi (массовые доли)"""
    a = air_Y()
    f = max(0.0, phi) * far_st
    Y = {k: v / (1 + f) for k, v in a.items()}
    Y[fuel] = f / (1 + f)
    return Y


def rho_of(p, T, Y):
    mw = 1.0 / sum(Y[k] / MW[k] for k in Y if Y[k] > 0)
    return p * mw * 1e-3 / (R_U * T)


# =========================================================================== блочная сетка 2D (r, z), клин
class Blocks:
    """многоблочная сетка по сетке отрезков: zs - границы по оси, rs - по радиусу; mask[i][j] - есть ли блок
    (отрезок z i, отрезок r j); nz[i], nr[j] - число ячеек; gz[i], gr[j] - сгущение (simpleGrading).
    Грани блоков на границе получают имя участка через classify(i, j, side), side: rin, rout, zlo, zhi."""

    def __init__(self, zs, rs, mask, nz, nr, gz, gr):
        self.zs, self.rs, self.mask, self.nz, self.nr, self.gz, self.gr = zs, rs, mask, nz, nr, gz, gr
        self.vid = {}
        self.verts = []
        for k, z in enumerate(zs):
            for l, r in enumerate(rs):
                if r <= 1e-12:
                    self.vid[(k, l, -1)] = self.vid[(k, l, 1)] = self._add((0.0, 0.0, z))
                else:
                    self.vid[(k, l, -1)] = self._add((r, -r * TAN_A, z))
                    self.vid[(k, l, 1)] = self._add((r, r * TAN_A, z))

    def _add(self, p):
        self.verts.append(p)
        return len(self.verts) - 1

    def has(self, i, j):
        return 0 <= i < len(self.zs) - 1 and 0 <= j < len(self.rs) - 1 and self.mask[i][j]

    def hexv(self, i, j):
        v = self.vid
        return [v[(i, j, -1)], v[(i, j + 1, -1)], v[(i, j + 1, 1)], v[(i, j, 1)],
                v[(i + 1, j, -1)], v[(i + 1, j + 1, -1)], v[(i + 1, j + 1, 1)], v[(i + 1, j, 1)]]

    def n_cells(self):
        return sum(self.nz[i] * self.nr[j] for i in range(len(self.zs) - 1) for j in range(len(self.rs) - 1)
                   if self.mask[i][j])

    def dict_text(self, classify, patch_types):
        H = [f"scale 1;\n\nvertices\n("]
        for n, p in enumerate(self.verts):
            H.append(f"    ({p[0]:.9g} {p[1]:.9g} {p[2]:.9g}) // {n}")
        H.append(");\n\nblocks\n(")
        faces = {}
        wedge_n, wedge_p, axis = [], [], []
        for i in range(len(self.zs) - 1):
            for j in range(len(self.rs) - 1):
                if not self.mask[i][j]:
                    continue
                h = self.hexv(i, j)
                H.append(f"    hex ({' '.join(map(str, h))}) ({self.nr[j]} 1 {self.nz[i]}) "
                         f"simpleGrading ({self.gr[j]:.6g} 1 {self.gz[i]:.6g})")
                wedge_n.append((h[0], h[1], h[5], h[4]))
                wedge_p.append((h[3], h[7], h[6], h[2]))
                sides = dict(rin=((h[0], h[4], h[7], h[3]), (i, j - 1)), rout=((h[1], h[2], h[6], h[5]), (i, j + 1)),
                             zlo=((h[0], h[3], h[2], h[1]), (i - 1, j)), zhi=((h[4], h[5], h[6], h[7]), (i + 1, j)))
                for side, (f, nb) in sides.items():
                    if self.has(*nb):
                        continue
                    if side == "rin" and self.rs[j] <= 1e-12:
                        axis.append(f)
                        continue
                    name = classify(i, j, side)
                    faces.setdefault(name, []).append(f)
        H.append(");\n\nedges\n(\n);\n\nboundary\n(")
        for name, fl in faces.items():
            H.append(f"    {name}\n    {{\n        type {patch_types.get(name, 'patch')};\n        faces\n        (")
            H += [f"            ({' '.join(map(str, f))})" for f in fl]
            H.append("        );\n    }")
        for name, fl, t in (("axis", axis, "empty"), ("wedge_neg", wedge_n, "wedge"), ("wedge_pos", wedge_p, "wedge")):
            if fl:
                H.append(f"    {name}\n    {{\n        type {t};\n        faces\n        (")
                H += [f"            ({' '.join(map(str, f))})" for f in fl]
                H.append("        );\n    }")
        H.append(");\n\nmergePatchPairs\n(\n);\n")
        return "\n".join(H), sorted(faces)

    def face_centres(self, i, j, side):
        """радиусы центров граней сетки на стороне блока zlo/zhi (для неоднородного профиля на входе)"""
        r0, r1, n, g = self.rs[j], self.rs[j + 1], self.nr[j], self.gr[j]
        return [0.5 * (a + b) for a, b in zip(*_graded_edges(r0, r1, n, g))]


def _graded_edges(a, b, n, g):
    """края ячеек отрезка [a, b] при simpleGrading g (отношение последней ячейки к первой)"""
    if n == 1 or abs(g - 1) < 1e-12:
        e = [a + (b - a) * k / n for k in range(n + 1)]
    else:
        q = g ** (1.0 / (n - 1))
        w = [q ** k for k in range(n)]
        s = sum(w)
        e = [a]
        for x in w:
            e.append(e[-1] + (b - a) * x / s)
    return e[:-1], e[1:]


def _cells_for(L, dx, nmin, nmax):
    return int(max(nmin, min(nmax, math.ceil(L / max(dx, 1e-9)))))


# =========================================================================== словари OpenFOAM
HDR = ("/*--------------------------------*- C++ -*----------------------------------*\\\n"
       "  KV combustor funnel L1 - generated by routes/gtu/combustor/comb_cfd.py\n"
       "\\*---------------------------------------------------------------------------*/\n")


def foam_file(path, cls, obj, body, location=None):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    loc = f'    location    "{location}";\n' if location else ""
    with open(path, "w", encoding="utf-8") as f:
        f.write(HDR + "FoamFile\n{\n    version     2.0;\n    format      ascii;\n"
                f"    class       {cls};\n{loc}    object      {obj};\n}}\n\n" + body.rstrip() + "\n")


def _vec_list(v):
    return "nonuniform List<vector> " + f"{len(v)}(" + " ".join(f"({a:.9g} {b:.9g} {c:.9g})" for a, b, c in v) + ")"


def _scal_list(vals):
    return "nonuniform List<scalar> " + f"{len(vals)}(" + " ".join(f"{x:.9g}" for x in vals) + ")"


def field_file(path, name, dims, internal, bfs, cls="volScalarField"):
    body = [f"dimensions      {dims};\n", f"internalField   {internal};\n", "boundaryField\n{"]
    for patch, entries in bfs.items():
        body.append(f"    {patch}\n    {{")
        for k, v in entries.items():
            body.append(f"        {k:<16}{v};")
        body.append("    }")
    body.append("}")
    foam_file(path, cls, name, "\n".join(body), location="0")


# =========================================================================== геометрия по варианту
def cell_geometry(e, li):
    """размеры ячейки элемента и расходы на элемент при нагрузке li (всё в СИ)"""
    c, d, L = e["concept"], e["derived"], e["loads"][li]
    n_can = e["n_can"]
    n_el_can = max(1, c["n_el"] / n_can)
    W, Wf, p3, T3 = L["W_air"], L["W_fuel"], L["p3_Pa"], L["T3_K"]
    rho3 = L["rho3"]
    r_ch = 0.5e-3 * d["d_air_mm"]
    A_L_can = d["A_L"]                      # площадь жаровой трубы на одну трубу (кольцевая - всё кольцо)
    R_c = max(math.sqrt(A_L_can / (math.pi * n_el_can)), 1.05 * r_ch)
    g = dict(r_ch=r_ch, R_c=R_c, n_el_can=n_el_can, crowded=bool(R_c <= 1.05 * r_ch + 1e-12))
    l_ch = c["l_mix"] * 2 * r_ch
    g["l_ch"] = l_ch
    # зоны до разбавления на элемент -> длина ячейки
    V_up = sum(z["V_L"] for z in d["zones"] if z["name"] != "разбавление") * 1e-3 / n_can
    g["L_cell"] = max(4 * R_c, V_up / n_el_can / (math.pi * R_c ** 2))
    # расходы на элемент
    s = c["scheme"]
    u_j = 0.70 * math.sqrt(2 * c["dp_inj"] * p3 / rho3)       # скорость струи воздуха через отверстие (Cd 0,7, как comb_space)
    m_zone = c["a_pz"] * W / c["n_el"]
    x = L.get("x_pilot", 0.0)
    g["pilot"] = g["cav"] = g["quench"] = None
    if c["passive"] == "гидравлика" and s in ("LPP", "MM", "FLOX"):
        m_p = c["a_aux"] * m_zone
        g["pilot"] = dict(m_air=m_p, m_fuel=x * Wf / c["n_el"], r=r_ch * math.sqrt(c["a_aux"]))
        m_main = m_zone - m_p
    else:
        m_main = m_zone
    m_fuel_main = (1 - x) * Wf / c["n_el"]
    if s == "TVC":
        m_c = c["a_aux"] * W / c["n_el"]
        A_c = m_c / (rho3 * u_j)
        h = R_c - math.sqrt(max(R_c ** 2 - A_c / math.pi, (r_ch * 1.02) ** 2))
        g["cav"] = dict(m_air=m_c, m_fuel=x * Wf / c["n_el"], h=h)
    if s == "RQL":
        m_q = c["a_aux"] * W / c["n_el"]
        V_rich = next(z["V_L"] for z in d["zones"] if z["name"] == "богатая") * 1e-3 / n_can / n_el_can
        L_rich = V_rich / (math.pi * R_c ** 2)
        w_q = max(SLOT_MIN, m_q / (rho3 * u_j * 2 * math.pi * R_c))
        g["quench"] = dict(m_air=m_q, z=l_ch + L_rich, w=w_q, u=m_q / (rho3 * w_q * 2 * math.pi * R_c))
        g["L_cell"] = max(g["L_cell"], L_rich + 3 * w_q + 2 * R_c)
    g.update(m_air_main=m_main, m_fuel_main=m_fuel_main, rho3=rho3, u_jet_air=u_j)
    r_in = g["pilot"]["r"] if g["pilot"] else 0.0
    g["r_in"] = r_in
    g["u_ch"] = m_main / (rho3 * math.pi * (r_ch ** 2 - r_in ** 2))
    # топливо: кольцевая щель (k_je >= 4) или радиальный профиль смеси (k_je < 4 - дискретные струи, L1 их не разрешает)
    g["fuel_mode"] = "ring" if RING and c["k_je"] >= 4 and s != "CAT" else "profile"
    if g["fuel_mode"] == "ring":
        Tf, pf = 288.15, p3
        rho_f = pf * MW["CH4"] * 1e-3 / (R_U * Tf)
        c_f = math.sqrt(1.31 * 518.28 * Tf)
        # скорость щели - по отношению импульсов J варианта (derive) к скорости в канале CFD, не выше 0,8 скорости звука
        u_f = g_uch = m_main / (rho3 * math.pi * (r_ch ** 2 - (g["pilot"]["r"] if g["pilot"] else 0.0) ** 2))
        u_f = min(g_uch * math.sqrt(max(d["J"], 0.05) * rho3 / rho_f), 0.8 * c_f)
        w = m_fuel_main / (rho_f * u_f * 2 * math.pi * r_ch)
        w = max(w, 2e-5)
        u_f = m_fuel_main / (rho_f * w * 2 * math.pi * r_ch)
        g["ring"] = dict(z0=min(0.15 * l_ch, 2 * r_ch), w=w, u=u_f, rho=rho_f,
                         J_ring=rho_f * u_f ** 2 / (rho3 * g["u_ch"] ** 2), J_design=d["J"])
    return g


def _profile_phi(e, li, g, radii):
    """радиальный профиль phi на входе канала (режим profile): доля сечения a_eff встречает топливо (самостратификация
    струй, comb_crn.strat_air_frac), внутри - неравномерность S (comb_space.derive, Холдеман) косинусом по радиусу;
    средняя phi по расходу воздуха совпадает с заданной"""
    L, c = e["loads"][li], e["concept"]
    phi_m = (g["m_fuel_main"] / g["m_air_main"]) / e["far_st"]
    a_eff = L.get("a_eff", 1.0)
    S = e["derived"]["S_unmix"]
    r_in, r_ch = g["r_in"], g["r_ch"]
    A = r_ch ** 2 - r_in ** 2
    r_cut = math.sqrt(r_in ** 2 + a_eff * A)
    out = []
    for r in radii:
        if r > r_cut:
            out.append(0.0)
            continue
        xi = (r ** 2 - r_in ** 2) / max(r_cut ** 2 - r_in ** 2, 1e-30)   # 0..1 по площади
        out.append(phi_m / a_eff * (1 + math.sqrt(2) * S * math.cos(math.pi * xi)))
    return [max(0.0, v) for v in out]


def build_cell_case(e, li, out, level="coarse", of="2512", iters=None):
    """кейс «ячейка» варианта e при нагрузке li -> папка out; возвращает описание (геометрия, сетка, расходы)"""
    it = dict(ITERS, **(iters or {}))
    M = MESH[level]
    g = cell_geometry(e, li)
    L = e["loads"][li]
    p3, T3 = L["p3_Pa"], L["T3_K"]
    r_in, r_ch, R_c, l_ch, L_c = g["r_in"], g["r_ch"], g["R_c"], g["l_ch"], g["L_cell"]
    # отрезки по радиусу
    rs = [0.0]
    if r_in > 0:
        rs.append(r_in)
    rs.append(r_ch)
    h_cav = g["cav"]["h"] if g["cav"] else 0.0
    if h_cav > 0 and R_c - h_cav > r_ch * 1.01:
        rs.append(R_c - h_cav)
    rs.append(R_c)
    j_ch_lo = 1 if r_in > 0 else 0                      # первый отрезок канала
    j_ch_hi = rs.index(r_ch) - 1                       # последний отрезок канала
    # отрезки по оси
    zs = [0.0]
    ring = g.get("ring")
    if ring:
        zs += [ring["z0"], ring["z0"] + ring["w"]]
    zs.append(l_ch)
    q = g["quench"]
    if q:
        zs += [q["z"], q["z"] + q["w"]]
    zs.append(l_ch + L_c)
    zs = sorted(set(round(z, 12) for z in zs))
    i_dump = zs.index(round(l_ch, 12))
    dr = (r_ch - r_in) / M["nr_ch"]
    nr = []
    gr = []
    for j in range(len(rs) - 1):
        a, b = rs[j], rs[j + 1]
        if j_ch_lo <= j <= j_ch_hi:
            nr.append(M["nr_ch"]); gr.append(1.0)
        elif j < j_ch_lo:
            nr.append(_cells_for(b - a, dr, 3, M["nr_ch"])); gr.append(1.0)
        else:
            n = _cells_for(b - a, 1.6 * dr, 4, M["nr_out_max"])
            nr.append(n); gr.append(min(8.0, max(1.0, (b - a) / (n * dr))) if b == rs[-1] or h_cav == 0 else 1.0)
    nz, gz = [], []
    for i in range(len(zs) - 1):
        a, b = zs[i], zs[i + 1]
        if b <= l_ch + 1e-12:
            nz.append(_cells_for(b - a, 1.5 * dr, 2, M["nz_ch_max"])); gz.append(1.0)
        elif q and abs(a - q["z"]) < 1e-12:
            nz.append(3); gz.append(1.0)
        else:
            n = M["nz_cell"] if b >= zs[-1] - 1e-12 else _cells_for(b - a, 2 * dr, 4, M["nz_cell"])
            gzv = min(30.0, max(1.0, (b - a) / n / (1.5 * dr)) ** 1.0) if a <= l_ch + 1e-12 else 1.0
            nz.append(n); gz.append(gzv)
    mask = [[(zs[i + 1] > l_ch + 1e-12) or (j_ch_lo <= j <= j_ch_hi) for j in range(len(rs) - 1)]
            for i in range(len(zs) - 1)]
    outer_wall = g["n_el_can"] <= 1.0 + 1e-9

    def classify(i, j, side):
        z0, z1 = zs[i], zs[i + 1]
        if side == "zlo" and i == 0:
            return "inlet_air"
        if side == "zhi" and i == len(zs) - 2:
            return "outlet"
        if side == "zlo" and i == i_dump:
            if r_in > 0 and j < j_ch_lo:
                return "inlet_pilot"
            if g["cav"] and rs[j] >= R_c - h_cav - 1e-12 and h_cav > 0:
                return "inlet_cavity"
            return "wall_plate"
        if side == "rout" and j == j_ch_hi and z1 <= l_ch + 1e-12:
            if ring and abs(z0 - ring["z0"]) < 1e-12:
                return "inlet_fuel"
            return "wall_channel"
        if side == "rin" and j == j_ch_lo and r_in > 0 and z1 <= l_ch + 1e-12:
            return "wall_channel"
        if side == "rout" and j == len(rs) - 2:
            if q and abs(z0 - q["z"]) < 1e-12:
                return "inlet_quench"
            return "wall_liner" if outer_wall else "side_sym"
        return "wall_channel"

    ptypes = dict(wall_plate="wall", wall_channel="wall", wall_liner="wall", side_sym="symmetry")
    B = Blocks(zs, rs, mask, nz, nr, gz, gr)
    case = os.path.abspath(out)
    if os.path.exists(case):
        shutil.rmtree(case)
    os.makedirs(os.path.join(case, "system"))
    text, patches = B.dict_text(classify, ptypes)
    foam_file(os.path.join(case, "system", "blockMeshDict"), "dictionary", "blockMeshDict", text)
    # профили на входе
    j_rng = range(j_ch_lo, j_ch_hi + 1)
    radii = [r for j in j_rng for r in B.face_centres(0, j, "zlo")]
    phi_prof = _profile_phi(e, li, g, radii) if g["fuel_mode"] == "profile" else [0.0] * len(radii)
    sw = e["concept"]["swirl"]
    u = g["u_ch"]
    rho_in = []
    Yprof = []
    for ph in phi_prof:
        Yp = mix_Y(ph, e["far_st"])
        Yprof.append(Yp)
        rho_in.append(rho_of(p3, T3, Yp))
    # скорость: массовый расход по профилю плотности; крутка - твёрдое вращение w = 2·S·u·r/R (S = Ω·R/(2u))
    A_faces = [math.pi * (b * b - a * a) for j in j_rng for a, b in zip(*_graded_edges(rs[j], rs[j + 1], nr[j], gr[j]))]
    m_air_unit = sum(rho * (1 - Y.get("CH4", 0.0)) * A for rho, Y, A in zip(rho_in, Yprof, A_faces))
    k_u = g["m_air_main"] / max(m_air_unit, 1e-30)      # осевая скорость: расход воздуха элемента по профилю плотности
    U_in = [(0.0, 2 * sw * u * r / r_ch, k_u) for r in radii]
    desc = dict(kind="cell", load=L["load"], level=level, of=of, zs=zs, rs=rs, nz=nz, nr=nr, n_cells=B.n_cells(),
                patches=patches, geometry={k: v for k, v in g.items() if not isinstance(v, dict) or v is None},
                pilot=g["pilot"], cavity=g["cav"], quench=g["quench"], ring=ring, fuel_mode=g["fuel_mode"],
                phi_profile=dict(r=radii, phi=phi_prof) if g["fuel_mode"] == "profile" else None,
                z_dump=l_ch, z_end=zs[-1], r_cell=R_c, p3=p3, T3=T3)
    inlets = {"inlet_air": dict(U=U_in, T=T3, Y=Yprof, kind="air")}
    if ring:
        inlets["inlet_fuel"] = dict(U=(-ring["u"], 0.0, 0.0), T=288.15, Y={"CH4": 1.0}, kind="fuel")
    if g["pilot"]:
        P = g["pilot"]
        phi_p = (P["m_fuel"] / P["m_air"]) / e["far_st"] if P["m_air"] > 0 else 0.0
        Yp = mix_Y(phi_p, e["far_st"])
        up = (P["m_air"] + P["m_fuel"]) / (rho_of(p3, T3, Yp) * math.pi * P["r"] ** 2)
        inlets["inlet_pilot"] = dict(U=(0.0, 0.0, up), T=T3, Y=Yp, kind="pilot")
        desc["phi_pilot"] = phi_p
    if g["cav"]:
        C = g["cav"]
        phi_c = (C["m_fuel"] / C["m_air"]) / e["far_st"]
        Yc = mix_Y(phi_c, e["far_st"])
        A_c = math.pi * (R_c ** 2 - (R_c - C["h"]) ** 2)
        inlets["inlet_cavity"] = dict(U=(0.0, 0.0, (C["m_air"] + C["m_fuel"]) / (rho_of(p3, T3, Yc) * A_c)), T=T3,
                                      Y=Yc, kind="cavity")
        desc["phi_cavity"] = phi_c
    if q:
        inlets["inlet_quench"] = dict(U=(-q["u"], 0.0, 0.0), T=T3, Y=air_Y(), kind="quench")
    hot = L["hot_init"]
    _write_common(case, e, L, inlets, patches, of, it, hot_box=(0.0, l_ch, R_c * 1.5, l_ch + L_c), Lscale=r_ch)
    desc["inlets"] = {k: dict(kind=v["kind"]) for k, v in inlets.items()}
    desc["mass_flows_per_element"] = dict(air_main=g["m_air_main"], fuel_main=g["m_fuel_main"],
                                          pilot=g["pilot"], cavity=g["cav"], quench=g["quench"])
    json.dump(desc, open(os.path.join(case, "case.json"), "w", encoding="utf-8"), ensure_ascii=False, indent=1)
    return desc


def build_liner_case(e, cell_out, out, level="coarse", of="2512", iters=None):
    """кейс «труба» (зона разбавления) при 100 %: вход - осреднённый выход ячейки (cell_out: T, Y, U_z) и слой
    воздуха охлаждения у стенок; щели разбавления и охлаждения на стенках"""
    it = dict(ITERS, **(iters or {}))
    M = MESH[level]
    c, d, L = e["concept"], e["derived"], e["loads"][0]
    p3, T3, W, rho3 = L["p3_Pa"], L["T3_K"], L["W_air"], L["rho3"]
    n_can = e["n_can"]
    annular = d.get("H_L") is not None
    L_L = d["L_L"]
    g = cell_geometry(e, 0)
    z0 = min(g["L_cell"], 0.6 * L_L)
    Ld = max(L_L - z0, 2 * (d["H_L"] if annular else d["D_L"]) * 0.5)
    if annular:
        r_lo, r_hi = 0.5 * d["D_m"] - 0.5 * d["H_L"], 0.5 * d["D_m"] + 0.5 * d["H_L"]
    else:
        r_lo, r_hi = 0.0, 0.5 * d["D_L"]
    H = r_hi - r_lo
    u_j = 0.70 * math.sqrt(2 * c["dp_inj"] * p3 / rho3)
    m_cool_can = d["a_cool"] * W / n_can
    m_dil_can = max(0.0, d["a_dil"] * W / n_can)
    f_up = z0 / (z0 + Ld)                                # доля охлаждения, вошедшая выше по течению (слой у стенки)
    m_film_in = m_cool_can * f_up
    m_cool_down = m_cool_can - m_film_in
    walls = ["wall_in", "wall_out"] if annular else ["wall_out"]
    # слой охлаждающего воздуха на входе у стенок: толщина по скорости ядра
    m_core = cell_out["m_dot_can"]
    rho_hot = rho_of(p3, cell_out["T"], cell_out["Y"])
    A_tot = math.pi * (r_hi ** 2 - r_lo ** 2)
    u_core_est = m_core / (rho_hot * A_tot)
    A_film = min(0.3 * A_tot, m_film_in / (rho3 * max(u_core_est, 5.0)))   # слой охлаждения со скоростью ядра
    t_film = A_film / (2 * math.pi * r_hi) / len(walls)
    # щели: разбавление - одна на стенку на 0,15·Ld; охлаждение - ряды через ~25 мм ниже по течению
    z_dil = 0.15 * Ld
    w_dil = max(SLOT_MIN, m_dil_can / len(walls) / (rho3 * u_j * 2 * math.pi * r_hi))
    n_rows = max(2, min(12, int(Ld / 0.025)))
    z_rows = [0.3 * Ld + k * (0.65 * Ld) / max(1, n_rows - 1) for k in range(n_rows)]
    w_row = max(SLOT_MIN, m_cool_down / n_rows / len(walls) / (rho3 * u_j * 2 * math.pi * r_hi))
    zs = [0.0, z_dil, z_dil + w_dil] + sum(([z, z + w_row] for z in z_rows), []) + [Ld]
    zs = sorted(set(round(z, 9) for z in zs if 0 <= z <= Ld))
    slots = {round(z_dil, 9): "inlet_dil"}
    for z in z_rows:
        slots.setdefault(round(z, 9), "inlet_cool")
    rs = [r_lo]
    if annular:
        rs.append(r_lo + t_film)
    rs += [r_hi - t_film, r_hi]
    rs = sorted(set(round(r, 9) for r in rs))
    dz0 = H / M["nr_liner"]
    nr = [max(3, int(round(M["nr_liner"] * (rs[j + 1] - rs[j]) / H))) for j in range(len(rs) - 1)]
    nz = [max(2, min(M["nz_liner"], int(math.ceil((zs[i + 1] - zs[i]) / (1.2 * dz0))))) for i in range(len(zs) - 1)]
    gz = [1.0] * len(nz)
    gr = [1.0] * len(nr)
    mask = [[True] * (len(rs) - 1) for _ in range(len(zs) - 1)]

    def classify(i, j, side):
        if side == "zlo":
            if (annular and j == 0) or j == len(rs) - 2:
                return "inlet_film"
            return "inlet_core"
        if side == "zhi":
            return "outlet"
        wall = "wall_out" if side == "rout" else "wall_in"
        z0_ = round(zs[i], 9)
        return slots.get(z0_, wall) + ("_o" if side == "rout" else "_i") if z0_ in slots else wall

    ptypes = {"wall_in": "wall", "wall_out": "wall"}
    B = Blocks(zs, rs, mask, nz, nr, gz, gr)
    case = os.path.abspath(out)
    if os.path.exists(case):
        shutil.rmtree(case)
    os.makedirs(os.path.join(case, "system"))
    text, patches = B.dict_text(classify, ptypes)
    foam_file(os.path.join(case, "system", "blockMeshDict"), "dictionary", "blockMeshDict", text)
    A_core = A_tot - A_film
    u_core = m_core / (rho_hot * A_core)
    inlets = {"inlet_core": dict(U=(0.0, 0.0, u_core), T=cell_out["T"], Y=cell_out["Y"], kind="core"),
              "inlet_film": dict(U=(0.0, 0.0, m_film_in / (rho3 * A_film)), T=T3, Y=air_Y(), kind="film")}
    for p in patches:
        if p.startswith("inlet_dil") or p.startswith("inlet_cool"):
            sign = -1.0 if p.endswith("_o") else 1.0
            r_w = r_hi if p.endswith("_o") else r_lo
            if "dil" in p:
                u_s = m_dil_can / len(walls) / (rho3 * w_dil * 2 * math.pi * r_w)
            else:
                u_s = m_cool_down / n_rows / len(walls) / (rho3 * w_row * 2 * math.pi * r_w)
            inlets[p] = dict(U=(sign * u_s, 0.0, 0.0), T=T3, Y=air_Y(), kind="dilution" if "dil" in p else "cooling")
    _write_common(case, e, L, inlets, patches, of, it, hot_box=None, Lscale=H)
    desc = dict(kind="liner", load=1.0, level=level, of=of, zs=zs, rs=rs, nz=nz, nr=nr, n_cells=B.n_cells(),
                patches=patches, annular=annular, z0_in_liner=z0, L_dil=Ld, r_lo=r_lo, r_hi=r_hi,
                m_core=m_core, m_film_in=m_film_in, m_cool_down=m_cool_down, m_dil=m_dil_can, u_jet=u_j,
                n_cool_rows=n_rows, w_dil=w_dil, w_row=w_row, t_film=t_film, p3=p3, T3=T3,
                T_core_in=cell_out["T"], inlets={k: dict(kind=v["kind"]) for k, v in inlets.items()})
    json.dump(desc, open(os.path.join(case, "case.json"), "w", encoding="utf-8"), ensure_ascii=False, indent=1)
    return desc


# =========================================================================== общие файлы кейса
def _write_common(case, e, L, inlets, patches, of, it, hot_box, Lscale, special=None, solver="reactingFoam"):
    p3, T3 = L["p3_Pa"], L["T3_K"]
    Ya = air_Y()
    # constant: механизм и модели
    cst = os.path.join(case, "constant")
    os.makedirs(cst, exist_ok=True)
    for f in ("reactionsGRI", "thermo.compressibleGasGRI", "radiationProperties"):
        shutil.copy(os.path.join(FOAM_RES, f), os.path.join(cst, f))
    foam_file(os.path.join(cst, "boundaryRadiationProperties"), "dictionary", "boundaryRadiationProperties",
              '".*"\n{\n    type            lookup;\n    emissivity      1;\n    absorptivity    0;\n}')
    foam_file(os.path.join(cst, "thermophysicalProperties"), "dictionary", "thermophysicalProperties",
              "thermoType\n{\n    type            hePsiThermo;\n    mixture         reactingMixture;\n"
              "    transport       sutherland;\n    thermo          janaf;\n    energy          sensibleEnthalpy;\n"
              "    equationOfState perfectGas;\n    specie          specie;\n}\n\ninertSpecie N2;\n\n"
              "chemistryReader foamChemistryReader;\nfoamChemistryFile \"<constant>/reactionsGRI\";\n"
              "foamChemistryThermoFile \"<constant>/thermo.compressibleGasGRI\";")
    foam_file(os.path.join(cst, "combustionProperties"), "dictionary", "combustionProperties",
              "combustionModel EDC;\n\nactive  true;\n\nEDCCoeffs\n{\n    version v2005;\n}")
    foam_file(os.path.join(cst, "turbulenceProperties"), "dictionary", "turbulenceProperties",
              "simulationType RAS;\n\nRAS\n{\n    RASModel        kEpsilon;\n    turbulence      on;\n    printCoeffs     on;\n}")
    foam_file(os.path.join(cst, "chemistryProperties"), "dictionary", "chemistryProperties",
              "chemistryType\n{\n    solver          ode;\n    method          TDAC;\n}\n\nchemistry       off;\n\n"
              "importantSpecies\n{\n    CO2;\n    H2O;\n    CH4;\n    O2;\n    NO;\n    CO;\n}\n\n"
              "initialChemicalTimeStep 1e-07;\n\nodeCoeffs\n{\n    solver          seulex;\n    absTol          1e-12;\n"
              "    relTol          1e-07;\n}\n\nreduction\n{\n    active          on;\n    log             off;\n"
              "    tolerance       0.0001;\n    method          DAC;\n    initialSet\n    {\n        CO;\n        CH4;\n"
              "        HO2;\n        NO;\n    }\n    automaticSIS    off;\n    fuelSpecies\n    {\n        CH4 1;\n    }\n}\n\n"
              "tabulation\n{\n    active          on;\n    log             off;\n    printProportion off;\n"
              "    printNumRetrieve off;\n    tolerance       0.003;\n    method          ISAT;\n    scaleFactor\n    {\n"
              "        otherSpecies    1;\n        Temperature     10000;\n        Pressure        1e+15;\n"
              "        deltaT          1;\n    }\n    maxNLeafs       5000;\n    chPMaxLifeTime  1000;\n"
              "    maxGrowth       100;\n    checkEntireTreeInterval 500;\n    maxDepthFactor  2;\n"
              "    minBalanceThreshold 30;\n    MRURetrieve     false;\n    maxMRUSize      0;\n    growPoints      true;\n"
              "    maxNumNewDim    10;\n}")
    foam_file(os.path.join(cst, "fvOptions"), "dictionary", "fvOptions",
              "radiation\n{\n    type            radiation;\n    libs            (radiationModels);\n}")
    foam_file(os.path.join(cst, "g"), "uniformDimensionedVectorField", "g",
              "dimensions      [0 1 -2 0 0 0 0];\nvalue           (0 0 0);")
    # system
    sysd = os.path.join(case, "system")
    funcs = "" if of == "1912" else (
        "\nfunctions\n{\n    moleFractions\n    {\n        type            psiReactionThermoMoleFractions;\n"
        "        libs            (fieldFunctionObjects);\n        executeControl  writeTime;\n"
        "        writeControl    writeTime;\n    }\n}\n")
    foam_file(os.path.join(sysd, "controlDict"), "dictionary", "controlDict",
              f"application     {solver};\nstartFrom       latestTime;\nstartTime       0;\nstopAt          endTime;\n"
              f"endTime         {it['cold']};\ndeltaT          1;\nwriteControl    timeStep;\nwriteInterval   {it['write']};\n"
              f"purgeWrite      3;\nwriteFormat     ascii;\nwritePrecision  9;\nwriteCompression off;\ntimeFormat      general;\n"
              f"timePrecision   8;\nrunTimeModifiable true;\n" + funcs)
    foam_file(os.path.join(sysd, "fvSchemes"), "dictionary", "fvSchemes",
              "ddtSchemes\n{\n    default         localEuler;\n}\n\ngradSchemes\n{\n    default         Gauss linear;\n}\n\n"
              "divSchemes\n{\n    default             none;\n    div(phi,U)          Gauss limitedLinearV 1;\n"
              "    div(phi,Yi)         Gauss limitedLinear01 1;\n    div(phi,h)          Gauss limitedLinear 1;\n"
              "    div(phi,K)          Gauss limitedLinear 1;\n    div(phid,p)         Gauss limitedLinear 1;\n"
              "    div(phi,epsilon)    Gauss limitedLinear 1;\n    div(phi,Yi_h)       Gauss limitedLinear01 1;\n"
              "    div(phi,k)          Gauss limitedLinear 1;\n"
              "    div(((rho*nuEff)*dev2(T(grad(U)))))     Gauss linear;\n}\n\n"
              "laplacianSchemes\n{\n    default         Gauss linear corrected;\n}\n\n"
              "interpolationSchemes\n{\n    default         linear;\n}\n\nsnGradSchemes\n{\n    default         corrected;\n}")
    foam_file(os.path.join(sysd, "fvSolution"), "dictionary", "fvSolution",
              "solvers\n{\n    \"rho.*\"\n    {\n        solver          diagonal;\n    }\n    p\n    {\n"
              "        solver          PCG;\n        preconditioner  DIC;\n        tolerance       1e-7;\n        relTol          0.01;\n"
              "    }\n    pFinal\n    {\n        $p;\n        relTol          0;\n    }\n    \"(U|h|k|epsilon)\"\n    {\n"
              "        solver          PBiCGStab;\n        preconditioner  DILU;\n        tolerance       1e-7;\n        relTol          0.05;\n"
              "    }\n    \"(U|h|k|epsilon)Final\"\n    {\n        $U;\n    }\n    Yi\n    {\n        solver          PBiCGStab;\n"
              "        preconditioner  DILU;\n        tolerance       1e-8;\n        relTol          0.05;\n    }\n    G\n    {\n"
              "        solver          PCG;\n        preconditioner  DIC;\n        tolerance       1e-5;\n        relTol          0.1;\n"
              "    }\n    GFinal\n    {\n        $G;\n        relTol          0;\n    }\n}\n\n"
              "PIMPLE\n{\n    momentumPredictor   yes;\n    nOuterCorrectors    1;\n    nCorrectors         2;\n"
              "    nNonOrthogonalCorrectors 0;\n    maxDeltaT           1e-4;\n    maxCo               0.3;\n"
              "    alphaTemp           0.05;\n    alphaY              0.05;\n    Yref\n    {\n        O2                  0.1;\n"
              "        CH4                 0.1;\n    }\n    rDeltaTSmoothingCoeff   0.1;\n    rDeltaTDampingCoeff     1;\n}\n\n"
              "relaxationFactors\n{\n    equations\n    {\n        \".*\"            1;\n    }\n}")
    if hot_box is not None:
        r0, z0, r1, z1 = hot_box
        hot = L["hot_init"]
        hv = "\n".join(f"            volScalarFieldValue {k:<6} {v:.9g}" for k, v in hot["Y"].items() if k in MAJOR)
        foam_file(os.path.join(sysd, "setFieldsDict"), "dictionary", "setFieldsDict",
                  f"defaultFieldValues\n(\n);\n\nregions\n(\n    boxToCell\n    {{\n"
                  f"        box ({-1:.9g} {-1:.9g} {z0:.9g}) ({r1:.9g} {1:.9g} {z1 + 1:.9g});\n        fieldValues\n        (\n"
                  f"            volScalarFieldValue T      {hot['T']:.9g}\n{hv}\n        );\n    }}\n);")
    # поля 0/
    z = os.path.join(case, "0")
    os.makedirs(z, exist_ok=True)
    wedge = {"type": "wedge"}
    sym = {"type": "symmetry"}
    empty = {"type": "empty"}

    def bfield(fn_inlet, fn_wall, fn_out):
        bf = {}
        for p in patches:
            if p.startswith("inlet"):
                bf[p] = fn_inlet(p)
            elif p.startswith("wall"):
                bf[p] = fn_wall(p)
            elif p == "outlet":
                bf[p] = fn_out(p)
            elif p == "side_sym":
                bf[p] = sym
        for sp_name, sp_entry in (special if special is not None else
                                  {"wedge_neg": wedge, "wedge_pos": wedge, "axis": empty}).items():
            bf[sp_name] = sp_entry
        return bf

    def U_in(p):
        v = inlets[p]["U"]
        return {"type": "fixedValue",
                "value": _vec_list(v) if isinstance(v, list) else f"uniform ({v[0]:.9g} {v[1]:.9g} {v[2]:.9g})"}

    field_file(os.path.join(z, "U"), "U", "[0 1 -1 0 0 0 0]", "uniform (0 0 0)",
               bfield(U_in, lambda p: {"type": "noSlip"},
                      lambda p: {"type": "pressureInletOutletVelocity", "value": "uniform (0 0 0)"}),
               cls="volVectorField")
    field_file(os.path.join(z, "p"), "p", "[1 -1 -2 0 0 0 0]", f"uniform {p3:.9g}",
               bfield(lambda p: {"type": "zeroGradient"}, lambda p: {"type": "zeroGradient"},
                      lambda p: {"type": "fixedValue", "value": f"uniform {p3:.9g}"}))
    field_file(os.path.join(z, "T"), "T", "[0 0 0 1 0 0 0]", f"uniform {T3:.9g}",
               bfield(lambda p: {"type": "fixedValue", "value": f"uniform {inlets[p]['T']:.9g}"},
                      lambda p: {"type": "zeroGradient"},
                      lambda p: {"type": "inletOutlet", "inletValue": f"uniform {T3:.9g}", "value": f"uniform {T3:.9g}"}))
    k0 = max(1e-3, 1.5 * (0.05 * 30.0) ** 2)
    eps0 = 0.09 ** 0.75 * k0 ** 1.5 / (0.07 * Lscale)
    field_file(os.path.join(z, "k"), "k", "[0 2 -2 0 0 0 0]", f"uniform {k0:.9g}",
               bfield(lambda p: {"type": "turbulentIntensityKineticEnergyInlet", "intensity": "0.05", "value": "uniform 1"},
                      lambda p: {"type": "kqRWallFunction", "value": f"uniform {k0:.9g}"},
                      lambda p: {"type": "inletOutlet", "inletValue": f"uniform {k0:.9g}", "value": f"uniform {k0:.9g}"}))
    field_file(os.path.join(z, "epsilon"), "epsilon", "[0 2 -3 0 0 0 0]", f"uniform {eps0:.9g}",
               bfield(lambda p: {"type": "turbulentMixingLengthDissipationRateInlet",
                                 "mixingLength": f"{0.07 * Lscale:.9g}", "value": "uniform 1"},
                      lambda p: {"type": "epsilonWallFunction", "value": f"uniform {eps0:.9g}"},
                      lambda p: {"type": "inletOutlet", "inletValue": f"uniform {eps0:.9g}", "value": f"uniform {eps0:.9g}"}))
    field_file(os.path.join(z, "nut"), "nut", "[0 2 -1 0 0 0 0]", "uniform 0",
               bfield(lambda p: {"type": "calculated", "value": "uniform 0"},
                      lambda p: {"type": "nutkWallFunction", "value": "uniform 0"},
                      lambda p: {"type": "calculated", "value": "uniform 0"}))
    field_file(os.path.join(z, "alphat"), "alphat", "[1 -1 -1 0 0 0 0]", "uniform 0",
               bfield(lambda p: {"type": "calculated", "value": "uniform 0"},
                      lambda p: {"type": "compressible::alphatWallFunction", "Prt": "0.85", "value": "uniform 0"},
                      lambda p: {"type": "calculated", "value": "uniform 0"}))
    bfG = bfield(lambda p: {"type": "MarshakRadiation", "T": "T", "value": "uniform 0"},
                 lambda p: {"type": "MarshakRadiation", "T": "T", "value": "uniform 0"},
                 lambda p: {"type": "MarshakRadiation", "T": "T", "value": "uniform 0"})
    field_file(os.path.join(z, "G"), "G", "[1 0 -3 0 0 0 0]", "uniform 0", bfG)
    for sp in MAJOR + ["Ydefault"]:
        def yin(p, sp=sp):
            Y = inlets[p]["Y"]
            if isinstance(Y, list):
                return {"type": "fixedValue", "value": _scal_list([y.get(sp, 0.0) for y in Y])}
            return {"type": "fixedValue", "value": f"uniform {(Y.get(sp, 0.0) if sp != 'Ydefault' else 0.0):.9g}"}
        y0 = Ya.get(sp, 0.0) if sp != "Ydefault" else 0.0
        field_file(os.path.join(z, sp), sp, "[0 0 0 0 0 0 0]", f"uniform {y0:.9g}",
                   bfield(yin, lambda p: {"type": "zeroGradient"},
                          lambda p, y0=y0: {"type": "inletOutlet", "inletValue": f"uniform {y0:.9g}",
                                            "value": f"uniform {y0:.9g}"}))
    # порядок запуска (тот же выполняет comb_l1.run_case)
    with open(os.path.join(case, "Allrun"), "w") as f:
        f.write(ALLRUN.format(cold=it["cold"], hot=it["cold"] + it["hot"], setf="runApplication setFields"
                              if hot_box is not None else ": # горячее поле не задаётся (вход уже горячий)"))
    os.chmod(os.path.join(case, "Allrun"), 0o755)


ALLRUN = """#!/bin/bash
# Порядок кейса уровня 1 (comb_cfd.py): сетка -> без химии -> горячее поле -> химия. Окружение OpenFOAM подключено.
cd "${{0%/*}}" || exit 1
runApplication() {{ app=$1; shift; "$app" "$@" > "log.$app" 2>&1 || {{ echo "ошибка $app"; exit 1; }}; }}
runApplication blockMesh
checkMesh > log.checkMesh 2>&1
foamDictionary constant/chemistryProperties -entry chemistry -set off > /dev/null
foamDictionary system/controlDict -entry endTime -set {cold} > /dev/null
reactingFoam > log.reactingFoam.cold 2>&1 || exit 2
{setf}
foamDictionary constant/chemistryProperties -entry chemistry -set on > /dev/null
foamDictionary system/controlDict -entry endTime -set {hot} > /dev/null
reactingFoam > log.reactingFoam.hot 2>&1 || exit 3
"""


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--inputs", required=True)
    ap.add_argument("--row", type=int, default=0)
    ap.add_argument("--load", type=int, default=0)
    ap.add_argument("--level", default="coarse", choices=list(MESH))
    ap.add_argument("--of", default="2512")
    ap.add_argument("--out", required=True)
    a = ap.parse_args(argv)
    e = json.load(open(a.inputs, encoding="utf-8"))["entries"][a.row]
    d = build_cell_case(e, a.load, a.out, a.level, a.of)
    print(json.dumps({k: d[k] for k in ("n_cells", "patches", "fuel_mode", "z_dump", "z_end", "r_cell")}, ensure_ascii=False))


if __name__ == "__main__":
    main()
