#!/usr/bin/env python3
"""Семейство «корпус-оболочка и рама» KV-3: геном → гладкая форма → дешёвые модели (мс на оценку).

Два варианта (решение владельца — считать оба): «оболочка» — цельная печатная оболочка газогенератора, ротор собирается
внутри (правило собираемости parts_model); «картридж» — наружный корпус с разъёмом, газогенератор (ротор, статорные
детали, горячая часть) вынимается целиком; канал корпуса должен пропускать картридж (монотонная огибающая).

Форма (только гладкие представления — урок organic/structure):
  * срединная поверхность стенки r_m(z) = E(z) + t(z)/2 + O(z): E — гладкая (LSE) огибающая требуемых радиусов зон,
    O ≥ 0 — B-сплайн смещений (8 точек, свойство выпуклой оболочки: O ≥ 0 ⇒ стенка не задевает внутренние узлы);
  * стенка — сплошная или «обшивка + решётка + обшивка» (гироид, алмаз, открытый Вороной) с переменной плотностью;
  * рёбра снаружи — рельеф C²-гладкого профиля со скруглением: классы «нет, осевые, кольца, ортогрид, диагрид,
    Вороной, ветвящиеся от опор (L-система, сужение по Мюррею — гипотеза)», порядок вращательной симметрии и
    асимметрия (искривление угла θ' = θ + a1·sin(θ − φ) + a2·sin(2θ − 1,3φ));
  * опоры: две боковые цапфы и задняя тяга; рамы подшипников — стойки-профили с навалом (lean), изгибом (bow) и
    конусностью (cone), бугорки на передней кромке передних стоек (гипотеза: в целях не учитываются);
  * горячая часть на «беличьей клетке» из гибких пальцев (радиально податливых, окружно жёстких) — податливая связь,
    которой не было в счёте organic/print_review (hot_integral: несовместность 3,39 мм по радиусу);
  * каналы: топливное дерево по Мюррею к 14 форсункам, масло в служебной стойке, отбор воздуха, датчики — в стенке
    или трубками снаружи.

Модели (все — оценка уровня облика, кроме взятых из файлов чисел): масса — интегрированием; частоты — оболочка
Доннелла с размазанными рёбрами по зонам и балочная модель корпуса на опорах с ротором на рамах и горячей частью на
клетке (обобщённая задача на собственные значения); прочность — мембранные напряжения от давления, температурные от
перепада по оси, концентрация у корня ребра, стойки (тепловое рассогласование с поворотом ступицы), пальцы клетки,
решётка (закон гироида — свой расчёт organic/structure); детали и крепёж — funnel/parts_model.count (одна точка правды,
с правилом собираемости) плюс правила семейства; ремонт — funnel/service_model.score; изготовимость — камеры
funnel/make_model (печать, ГИП, КТ), нависания, вынос порошка, литьё.
"""
import json
import math
import os
import sys

for _v in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")          # малые матрицы: многопоточный BLAS на занятой машине в десятки раз медленнее
import numpy as np  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
EVO = os.path.dirname(HERE)
GTU = os.path.dirname(EVO)
ROOT = os.path.dirname(os.path.dirname(GTU))
sys.path.insert(0, EVO)
sys.path.insert(0, os.path.join(GTU, "funnel"))
from core import features as FT  # noqa: E402
from core import shapes as S  # noqa: E402
import make_model as MK  # noqa: E402
import parts_model as PM  # noqa: E402
import service_model as SM  # noqa: E402

F3P = os.path.join(GTU, "frame3", "frame3_point_results.json")
MATS = os.path.join(GTU, "platform", "fea", "materials.json")
LAT = os.path.join(GTU, "organic", "structure", "lattice_results.json")
PRV = os.path.join(GTU, "organic", "print_review", "print_review_results.json")
MAN = os.path.join(GTU, "organic", "structure", "manifold_results.json")

VARIANTS = ("оболочка", "картридж")
F_MARGIN = 1.15          # частоты системы не ниже 1,15·f раб (принято для отбора, оценка; не норма API)
SF_YIELD = 1.5           # запас по пределу текучести (оценка, не нормативное значение)
K_ROUTE = {"печать": 0.9, "литьё_ti": 0.85, "литьё_сталь": 1.0}   # снижение прочности способа (оценка; печать — k_AM organic/structure)
K_MOUNT = 1.0e11         # Н/м: опоры пакета приняты жёсткими — считаются упругие моды корпуса и узлов; моды двигателя на раме
#                          пакета зависят от пакета (данных нет — запрос заказчику, навык gtu-dropin)
K_SPIGOT = 2.0e9         # Н/м, посадка картриджа в корпусе на плоскость (оценка)
C_KT = 0.15              # Kt = 1 + C_KT·√(H/r_f) у корня ребра (оценка; проверяется CalculiX)
SKIN_CHANNEL = 0.002     # м, стенка встроенного канала (оценка, как обшивка organic/structure)
T_STRUT_REAR = 430.0     # °C, стойки задней рамы в выхлопе (оценка: газ 458 °C)
T_HUB = (70.0, 150.0)    # °C, корпуса подшипников передний / задний (оценка)
R_HUB = 0.13             # м, наружный радиус корпуса подшипника (organic/structure/shell_topo GEO r_bear_o — оценка)
F_BEAR = 25e3            # Н, расчётная радиальная сила подшипника (shell_topo F_BEAR — оценка)
M_INTERNAL = 150.0       # кг, внутренние статорные детали зон компрессора (диффузоры, обратный канал; оценка)
SHOCK_G = 3.0            # перегрузка при перевозке и сейсмике для пальцев клетки (оценка)
BLOCKAGE_MAX = 0.15      # наибольшее загромождение проточной части стойками (оценка)
CAST_T_MIN = 0.008       # м, наименьшая стенка крупной стальной или титановой отливки (оценка)
PRINT_SKIN_MIN = 0.001   # м, наименьшая обшивка печатной стенки под давлением (оценка; элемент IN718 0,381 мм — print_review RULES)
PORE_MIN, PORE_LONG = 0.003, 0.005   # м, вынос порошка: ≥ 3 мм, при длине > 50 мм ≥ 5 мм (print_review RULES depowder_long_mm)


# =========================================================================== входы
def _interp_tab(tab, x):
    t = np.asarray(tab, float)
    return float(np.interp(x, t[:, 0], t[:, 1]))


def load_inputs():
    with open(F3P, encoding="utf-8") as f:
        fp = json.load(f)
    with open(MATS, encoding="utf-8") as f:
        mats = json.load(f)
    with open(LAT, encoding="utf-8") as f:
        lat = json.load(f)["law"]
    with open(PRV, encoding="utf-8") as f:
        prv = json.load(f)
    with open(MAN, encoding="utf-8") as f:
        man = json.load(f)
    s1, s2 = fp["stage1"], fp["stage2"]
    t3 = fp["turbine"]["options"][str(fp["turbine"]["selected_z"])]
    Dm = t3["Dm"]
    h3 = t3["stages"][-1]["h_exit_m"]
    D2max = max(s1["D2_m"], s2["D2_m"])
    r4r2 = fp["assumptions"]["r4r2"]
    PR = fp["split"]["pi"]
    p3 = 101.325e3 * PR
    p2 = s2["p01_kPa"] * 1e3
    cyc = fp["cycle"]["frame3"]
    # длины по правилу воронки (funnel_space: L_rot = Σ0,35·D2 + 0,15·Dm·n_t + 0,5·D2max); 0,5·D2max делится 0,3 / 0,7
    # между передней рамой и камерой (оценка)
    L_A = 0.3 * 0.5 * D2max
    L_B, L_C = 0.35 * s1["D2_m"], 0.35 * s2["D2_m"]
    L_D = 0.7 * 0.5 * D2max
    L_E = 0.15 * Dm * len(t3["stages"])
    zs = np.cumsum([0.0, L_A, L_B, L_C, L_D, L_E])
    Dtip = Dm + h3
    zones = [
        dict(name="передняя рама и вход", z0=zs[0], z1=zs[1], r_req=0.5 * 1.6 * s1["D1_tip_m"] + 0.02, p_g=0.0, T=s1["T01_C"],
             basis="r: входной коллектор D = 1,6·D1t (правило эскизов воронки) + 0,02 м; p — разрежение входа пренебрежимо"),
        dict(name="корпус ступени 1", z0=zs[1], z1=zs[2], r_req=0.5 * r4r2 * s1["D2_m"] + 0.04, p_g=p2 - 101.325e3, T=s1["T02_C"],
             basis="r: выход диффузора r4 = 1,55·r2 + 0,04 м на коллектор и обратный канал; p, T — за ступенью 1 (frame3)"),
        dict(name="корпус ступени 2", z0=zs[2], z1=zs[3], r_req=0.5 * r4r2 * s2["D2_m"] + 0.04, p_g=p3 - 101.325e3, T=s2["T02_C"],
             basis="r: r4 = 1,55·r2 + 0,04 м; p — за компрессором p3 = 101,325·πк; T — за ступенью 2 (frame3)"),
        dict(name="корпус камеры сгорания", z0=zs[3], z1=zs[4], r_req=0.5 * 1.1 * D2max + 0.05, p_g=p3 - 101.325e3, T=s2["T02_C"],
             basis="r: кольцевая камера D = 1,1·D2max (funnel_space) + 0,05 м; p, T — воздух за компрессором"),
        dict(name="корпус турбины", z0=zs[4], z1=zs[5], r_req=0.5 * 1.25 * Dtip + 0.03, p_g=p3 - 101.325e3, T=s2["T02_C"],
             basis="r: горячая часть D = 1,25·Dtip (funnel_space) + 0,03 м зазор клетки; p — воздух охлаждения за компрессором "
                   "(консервативно); T стенки — воздух за компрессором (горячая часть на клетке, оценка)"),
    ]
    # массы (формулы funnel_space.derive / gtu_explore.finalize — оценка)
    m_cw = 0.0825 * 4430.0 * (s1["D2_m"] ** 3 + s2["D2_m"] ** 3)
    m_tw = len(t3["stages"]) * 0.03 * 8000.0 * (1.15 * Dm) ** 3 * 1.5
    L_rot = float(zs[-1])
    m_shaft = 7850.0 * math.pi / 4 * (0.2 * D2max) ** 2 * L_rot * 0.5
    D_hot = 1.25 * Dtip; L_hot = 0.6 * D_hot
    m_hot = 8200.0 * math.pi * D_hot * L_hot * max(0.006, 0.012 * D_hot) + 0.02 * 8200.0 * Dm ** 3 * len(t3["stages"])
    hi = prv["hot_integral"]
    inp = dict(
        D2=[s1["D2_m"], s2["D2_m"]], D2max=D2max, Dm=Dm, Dtip=Dtip, PR=PR, p3=p3, p2=p2, n_rpm=s1["n_rpm"], f_op=s1["n_rpm"] / 60.0,
        P_shaft_kW=cyc["P_shaft_kW"], EGT=cyc["EGT_C"], m_air=s1["m_total_kg_s"], L=L_rot, zones=zones, zs=zs.tolist(),
        m_rotor=m_cw + m_tw + m_shaft, m_rotor_parts=dict(cw=m_cw, tw=m_tw, shaft=m_shaft), m_hot=m_hot,
        r_hot=0.5 * D_hot, dr_hot=hi["dr_free_mm_at_r"]["dr_mm"] * 1e-3 * (0.5 * D_hot) / hi["dr_free_mm_at_r"]["r_m"],
        dr_hot_src=hi["dr_free_mm_at_r"], z_fb=0.5 * L_A, z_rb=zs[5] - 0.25 * L_E, rotor_cg=0.5 * L_rot, rotor_L=0.95 * L_rot,
        lat=lat, mats=dict(ti64=mats["ti64"], in713c=mats["in713c"]),
        r_tur_hub=0.5 * (Dm - h3), r_tur_tip=0.5 * (Dm + h3),
        fuel=dict(n_inj=PM.REF_BOM["comb"]["n_inj"], d_leaf=man["tree"]["d_mm"][-1] * 1e-3 if "d_mm" in man.get("tree", {}) else 0.0141,
                  k=2.41),
        sources=dict(point="routes/gtu/frame3/frame3_point_results.json", materials="routes/gtu/platform/fea/materials.json",
                     lattice="routes/gtu/organic/structure/lattice_results.json → law",
                     hot="routes/gtu/organic/print_review/print_review_results.json → hot_integral",
                     parts="routes/gtu/funnel/parts_model.py", service="routes/gtu/funnel/service_model.py",
                     printers="routes/gtu/funnel/make_model.py (PRINTERS, HIP_VESSELS, SHELL_RULES)"),
    )
    return inp


INP = load_inputs()


# =========================================================================== материалы
def mat_props(key, T):
    """E (Па), ν, α (1/К), ρ, σ0,2 (Па), σв (Па) при T °C; статус — в MAT_NOTES"""
    if key == "ti64":
        m = INP["mats"]["ti64"]
        return dict(E=_interp_tab(m["E_GPa"], T) * 1e9, nu=m["nu"], alpha=_interp_tab(m["alpha_mean_1e6"], T) * 1e-6, rho=m["density"],
                    sy=_interp_tab(m["ys_MPa"], T) * 1e6, su=_interp_tab(m["uts_MPa"], T) * 1e6)
    if key == "in718p":           # печатный IN718 — свойства IN-713C из materials.json как заместитель (как print_review)
        m = INP["mats"]["in713c"]
        return dict(E=_interp_tab(m["E_GPa"], T) * 1e9, nu=m["nu"], alpha=_interp_tab(m["alpha_mean_1e6"], T) * 1e-6, rho=8190.0,
                    sy=_interp_tab(m["ys_MPa"], T) * 1e6, su=_interp_tab(m["uts_MPa"], T) * 1e6)
    if key == "zg15":             # ZG15Cr1Mo1V: σ0,2 ≥ 440 МПа при 20 °C (casting-spec, по выдержке поиска); остальное — оценка
        k = 1.0 - 0.25 * min(max((T - 20.0) / 400.0, 0.0), 1.2)
        return dict(E=(200.0 - 0.05 * (T - 20.0)) * 1e9, nu=0.3, alpha=12.5e-6 + 2e-9 * (T - 20.0), rho=7850.0, sy=440e6 * k, su=590e6 * k)
    raise KeyError(key)


MAT_NOTES = {
    "ti64": "materials.json ti64 (E при 400 °C — экстраполяция; σ при 300 и 427 °C — оценка, там же)",
    "in718p": "печатный IN718 для пальцев клетки — свойства IN-713C из materials.json как заместитель (как organic/print_review), "
              "плотность 8190 (make_model MAT IN718)",
    "zg15": "ZG15Cr1Mo1V: σ0,2 ≥ 440, σв ≥ 590 МПа при 20 °C — casting-spec (по выдержке поиска); снижение с температурой "
            "(−25 % к 420 °C), E, α, ρ — оценка",
}


def lattice_rel(rho):
    """E*/E и S*/σ (своя гомогенизация гироида, organic/structure lattice_results.json → law; для алмаза и Вороного
    принят тот же закон — оценка)"""
    L = INP["lat"]
    return float(np.interp(rho, L["rho"], L["E"])), float(np.interp(rho, L["rho"], L["S"]))


# =========================================================================== геном
G = []


def _g(name, lo, hi, kind="real", opts=None):
    G.append(dict(name=name, lo=lo, hi=hi, kind=kind, opts=opts))


for k in range(8):
    _g(f"off{k}", 0.0, 0.12)
for k in range(4):
    _g(f"t{k}", 0.004, 0.045)
for k in range(4):
    _g(f"core{k}", 0.0, 1.0)
_g("skin", 0.0008, 0.006)
_g("rho0", 0.10, 0.50); _g("rho1", 0.10, 0.50)
_g("lat_type", 0, 1, "cat", ["gyroid", "diamond", "voronoi"])
_g("cell", 0.006, 0.030)
_g("rib_type", 0, 1, "cat", ["нет", "осевые", "кольца", "ортогрид", "диагрид", "Вороной", "ветвящиеся"])
_g("rib_n", 0.0, 1.0); _g("rib_n2", 0.0, 1.0)
_g("rib_h", 0.0, 0.08); _g("rib_w", 0.004, 0.030); _g("rib_fillet", 0.002, 0.030)
_g("rib_angle", 15.0, 75.0); _g("rib_jitter", 0.0, 1.0)
_g("sym_n", 0, 1, "cat", [1, 2, 3, 4, 5, 6, 7, 8, 10, 12])
_g("asym1", 0.0, 0.6); _g("asym2", 0.0, 0.4); _g("asym_phase", 0.0, math.pi)
_g("mount_z", 0.25, 0.65); _g("mount_drop", 0.0, 45.0); _g("rear_z", 0.80, 0.97)
_g("f_n", 3, 12, "int"); _g("f_lean", 0.0, 40.0); _g("f_bow", 0.0, 0.2); _g("f_cone", 0.0, 50.0)
_g("f_chord", 0.04, 0.16); _g("f_tc", 0.08, 0.30); _g("f_tub", 0.0, 0.10)
_g("r_n", 3, 14, "int"); _g("r_lean", 0.0, 45.0); _g("r_bow", 0.0, 0.2); _g("r_cone", 0.0, 50.0)
_g("r_chord", 0.05, 0.18); _g("r_tc", 0.08, 0.30)
_g("hf_n", 6, 48, "int"); _g("hf_t", 0.0015, 0.008); _g("hf_b", 0.01, 0.06); _g("hf_shape", 1.0, 3.0); _g("hf_span", 0.05, 0.20)
_g("emb", 0, 1, "int6")
_g("fuel_port", 0.0, 1.0); _g("murray_k", 2.0, 3.0)
_g("route", 0, 1, "cat", ["печать", "литьё", "гибрид"])
_g("cast_mat", 0, 1, "cat", ["zg15", "ti64"])
_g("build_dir", 0, 1, "cat", ["+z", "-z"])
_g("extract", 0, 1, "cat", ["вперёд", "назад"])
DIM = len(G)
GENE_INDEX = {g["name"]: i for i, g in enumerate(G)}


def decode(x):
    x = np.clip(np.asarray(x, float), 0.0, 1.0)
    p = {}
    for v, g in zip(x, G):
        if g["kind"] == "real":
            p[g["name"]] = g["lo"] + v * (g["hi"] - g["lo"])
        elif g["kind"] == "int":
            p[g["name"]] = int(min(g["hi"], g["lo"] + math.floor(v * (g["hi"] - g["lo"] + 1))))
        elif g["kind"] == "int6":
            p[g["name"]] = int(min(5, math.floor(v * 6)))
        elif g["kind"] == "cat":
            o = g["opts"]
            p[g["name"]] = o[min(len(o) - 1, int(math.floor(v * len(o))))]
    return p


def encode(p):
    """обратное отображение (середина бина для дискретных генов) — для опорных форм и выбора владельца"""
    x = np.zeros(DIM)
    for i, g in enumerate(G):
        v = p[g["name"]]
        if g["kind"] == "real":
            x[i] = (v - g["lo"]) / (g["hi"] - g["lo"])
        elif g["kind"] == "int":
            x[i] = (v - g["lo"] + 0.5) / (g["hi"] - g["lo"] + 1)
        elif g["kind"] == "int6":
            x[i] = (v + 0.5) / 6
        else:
            x[i] = (g["opts"].index(v) + 0.5) / len(g["opts"])
    return np.clip(x, 0.0, 1.0)


# =========================================================================== геометрия
NZ = 161
ZG = np.linspace(0.0, INP["zs"][-1], NZ)


def _zone_of(z):
    zs = INP["zs"]
    return np.clip(np.searchsorted(zs, z, side="right") - 1, 0, len(zs) - 2)


def envelope(variant, extract="вперёд", beta=120.0, ramp=0.06):
    """гладкая огибающая требуемых радиусов: LSE от «плато» зон (≥ max по построению, превышение ≤ ln(n)/β ≈ 1,3 см).
    Картридж: плато продлено в сторону выемки — канал не сужается по пути картриджа."""
    z = ZG
    fs = []
    for k, zn in enumerate(INP["zones"]):
        a, b = zn["z0"], zn["z1"]
        if variant == "картридж":
            if extract == "вперёд":
                a = 0.0
            else:
                b = INP["zs"][-1]
        left = S.smoothstep5((z - (a - ramp)) / ramp)
        right = 1.0 - S.smoothstep5((z - b) / ramp)
        plateau = np.where(z < a, left, np.where(z > b, right, 1.0))
        fs.append(zn["r_req"] * plateau + 0.05 * (1 - plateau))
    F = np.stack(fs)
    m = F.max(axis=0)
    return m + np.log(np.exp(beta * (F - m)).sum(axis=0)) / beta


ENV = {("оболочка", "вперёд"): envelope("оболочка"), ("оболочка", "назад"): envelope("оболочка"),
       ("картридж", "вперёд"): envelope("картридж", "вперёд"), ("картридж", "назад"): envelope("картридж", "назад")}
R_REQ = np.array([INP["zones"][k]["r_req"] for k in _zone_of(ZG)])


def profiles(p, variant):
    u = ZG / ZG[-1]
    E = ENV[(variant, p["extract"] if variant == "картридж" else "вперёд")]
    O = S.bspline_eval(np.array([p[f"off{k}"] for k in range(8)])[:, None], u)[:, 0]
    t = S.bspline_eval(np.array([p[f"t{k}"] for k in range(4)])[:, None], u)[:, 0]
    c = S.bspline_eval(np.array([p[f"core{k}"] for k in range(4)])[:, None], u)[:, 0]
    rho = p["rho0"] + (p["rho1"] - p["rho0"]) * u
    r_m = E + 0.5 * t + O
    core = np.clip((c - 0.45) / 0.1, 0.0, 1.0)                 # плавная граница зон с решёткой
    if p["route"] != "печать":
        core = np.zeros_like(core)                              # решётка только в печатной стенке
    t_core = np.maximum(t - 2 * p["skin"], 0.0) * core
    return dict(z=ZG, E=E, O=O, t=t, r_m=r_m, r_i=r_m - 0.5 * t, r_o=r_m + 0.5 * t, core=core, t_core=t_core, rho=rho)


# --------------------------------------------------------------------------- рёбра (отрезки на развёртке (z, s = R_ref·θ))
def _halton(i, b):
    f, r = 1.0, 0.0
    while i > 0:
        f /= b; r += f * (i % b); i //= b
    return r


def _subdiv(z0, s0, z1, s1, w, step=0.10):
    n = max(1, int(math.ceil(math.hypot(z1 - z0, s1 - s0) / step)))
    t = np.linspace(0, 1, n + 1)
    zz, ss = z0 + t * (z1 - z0), s0 + t * (s1 - s0)
    return [(zz[i], ss[i], zz[i + 1], ss[i + 1], w) for i in range(n)]


def rib_segments(p, R_ref, mounts):
    """отрезки рёбер (z0, s0, z1, s1, w) на развёртке; s ∈ [0, 2πR_ref) по окружности, θ = s/R_ref"""
    L = INP["zs"][-1]
    C = 2 * math.pi * R_ref
    typ = p["rib_type"]
    N = int(p["sym_n"])
    sec = C / N
    w = p["rib_w"]
    segs = []
    z_lo, z_hi = 0.02, L - 0.02
    if typ in ("осевые", "ортогрид"):
        k = 1 + int(p["rib_n"] * 3)
        for j in range(N):
            for i in range(k):
                s = j * sec + (i + 0.5) * sec / k
                segs += _subdiv(z_lo, s, z_hi, s, w)
    if typ in ("кольца", "ортогрид"):
        nr = 2 + int(p["rib_n2"] * 10)
        for i in range(nr):
            z = z_lo + (i + 0.5) * (z_hi - z_lo) / nr
            segs += _subdiv(z, 0.0, z, C, w, step=0.08)
    if typ == "диагрид":
        k = 1 + int(p["rib_n"] * 2)
        tb = math.tan(math.radians(p["rib_angle"]))
        for j in range(N):
            for i in range(k):
                s0 = j * sec + i * sec / k
                for sg in (1.0, -1.0):
                    segs += _subdiv(z_lo, s0, z_hi, s0 + sg * (z_hi - z_lo) * tb, w)
    if typ == "Вороной":
        from scipy.spatial import Voronoi
        ns = 3 + int(p["rib_n2"] * 12)
        jit = p["rib_jitter"]
        pts = []
        for i in range(ns):
            hz, hs = _halton(i + 1, 2), _halton(i + 1, 3)
            jz, js = _halton(i + 7, 5) - 0.5, _halton(i + 7, 7) - 0.5
            pts.append((z_lo + (z_hi - z_lo) * min(max(hz + 0.3 * jit * jz, 0.02), 0.98), sec * ((hs + 0.6 * jit * js) % 1.0)))
        base = np.array(pts)
        allp = []
        for j in range(-1, N + 1):
            allp.append(base + np.array([0.0, j * sec]))
        P = np.vstack(allp)
        P = np.vstack([P, P * np.array([-1, 1]) + np.array([2 * z_lo, 0]), P * np.array([-1, 1]) + np.array([2 * z_hi, 0])])
        vor = Voronoi(P)
        for (a, b) in vor.ridge_vertices:
            if a < 0 or b < 0:
                continue
            va, vb = vor.vertices[a], vor.vertices[b]
            if not (z_lo - 1e-9 <= va[0] <= z_hi + 1e-9 and z_lo - 1e-9 <= vb[0] <= z_hi + 1e-9):
                continue
            mid = 0.5 * (va[1] + vb[1])
            if 0.0 <= mid < C:
                segs += _subdiv(va[0], va[1], vb[0], vb[1], w)
    if typ == "ветвящиеся":
        it = 2 + int(p["rib_n2"] * 1.99)
        rules = {"X": "F[+X]F[-X]+X"}
        s_l = S.lsystem("X", rules, it)
        step = (0.12 + 0.20 * p["rib_n"]) / (1.6 ** it)
        rng = np.random.Generator(np.random.PCG64(20261115))
        for (zm, sm, head, mirror) in mounts:
            for mir in ((1, -1) if mirror else (1,)):
                tr = S.turtle2d(s_l, (zm, sm), head, step, mir * p["rib_angle"], shrink=1.0, rng=rng, jitter_deg=25.0 * p["rib_jitter"])
                for (a, b, dpt) in tr:
                    ww = w * 2.0 ** (-dpt / 3.0)              # сужение ветвей «по Мюррею» (k = 3) — гипотеза
                    segs += _subdiv(a[0], a[1], b[0], b[1], ww)
    # симметрия/асимметрия: искривление угла
    out = []
    a1, a2, ph = p["asym1"], p["asym2"], p["asym_phase"]
    for (z0, s0, z1, s1, ww) in segs:
        if not (min(z0, z1) >= -1e-9 and max(z0, z1) <= L + 1e-9):
            z0, z1 = np.clip([z0, z1], 0.0, L)
        th0, th1 = s0 / R_ref, s1 / R_ref
        if typ != "кольца":
            th0 = th0 + a1 * math.sin(th0 - ph) + a2 * math.sin(2 * th0 - 1.3 * ph)
            th1 = th1 + a1 * math.sin(th1 - ph) + a2 * math.sin(2 * th1 - 1.3 * ph)
        out.append((float(z0), th0 * R_ref, float(z1), th1 * R_ref, ww))
    if p["rib_h"] < 0.002 or typ == "нет":
        return []
    return out


def mounts_of(p, R_ref=0.55):
    """опоры: (z, s, направление роста рёбер °, зеркалить) — две боковые цапфы (ниже горизонта на mount_drop) и задняя тяга снизу"""
    L = INP["zs"][-1]
    zm, zr = p["mount_z"] * L, p["rear_z"] * L
    th_r = -math.radians(p["mount_drop"])
    th_l = math.pi + math.radians(p["mount_drop"])
    return [(zm, (th_r % (2 * math.pi)) * R_ref, 0.0, True), (zm, (th_l % (2 * math.pi)) * R_ref, 0.0, True),
            (zr, (1.5 * math.pi) * R_ref, 180.0, True)], dict(z_side=zm, z_rear=zr, th=[th_r, th_l, 1.5 * math.pi])


def rib_stats(segs, p, R_ref, area):
    """объём, ориентационные моменты (размазанные «толщины» рёбер по направлениям), угловое распределение"""
    H, rf = p["rib_h"], p["rib_fillet"]
    if not segs:
        return dict(V=0.0, Szz=0.0, Szs=0.0, Sss=0.0, Izz=0.0, Iss=0.0, Izs=0.0, L=0.0, n=0, ang=np.zeros(72), ent=0.0,
                    overhang_len=0.0, by_zone=np.zeros((5, 7)))
    a = np.array(segs)
    dz, ds = a[:, 2] - a[:, 0], a[:, 3] - a[:, 1]
    Ls = np.hypot(dz, ds)
    A = H * (a[:, 4] + rf)                       # площадь профиля рельефа H·(w + r_f): ∫ smoothstep5 = 1/2
    c, s = dz / np.maximum(Ls, 1e-12), ds / np.maximum(Ls, 1e-12)
    wv = A * Ls / area
    Ib = (a[:, 4] + rf) * H ** 3 / 12.0
    wi = Ib * Ls / area
    ang = np.zeros(72)
    th = ((0.5 * (a[:, 1] + a[:, 3])) / R_ref) % (2 * math.pi)
    np.add.at(ang, (th / (2 * math.pi) * 72).astype(int) % 72, A * Ls)
    phi = np.mod(np.arctan2(s, c), math.pi)
    hist, _ = np.histogram(phi, bins=12, range=(0, math.pi), weights=A * Ls)
    pr = hist / max(hist.sum(), 1e-30)
    ent = float(-np.sum(pr[pr > 0] * np.log(pr[pr > 0])) / math.log(12))
    zmid = 0.5 * (a[:, 0] + a[:, 2])
    zi = _zone_of(zmid)
    # по зонам: Σ A·L·c⁴, Σ A·L·s⁴, Σ A·L, Σ A·L·c²s², Σ I·L·c⁴, Σ I·L·s⁴, Σ I·L·c²s² (c, s — косинус и синус угла ребра к оси)
    by_zone = np.zeros((5, 7))
    AL, IL = A * Ls, Ib * Ls
    for k in range(5):
        m = zi == k
        by_zone[k] = [np.sum(AL[m] * c[m] ** 4), np.sum(AL[m] * s[m] ** 4), np.sum(AL[m]), np.sum(AL[m] * c[m] ** 2 * s[m] ** 2),
                      np.sum(IL[m] * c[m] ** 4), np.sum(IL[m] * s[m] ** 4), np.sum(IL[m] * c[m] ** 2 * s[m] ** 2)]
    return dict(V=float(np.sum(A * Ls)), Szz=float(np.sum(wv * c ** 4)), Szs=float(np.sum(wv * c * c * s * s)), Sss=float(np.sum(wv * s ** 4)),
                Izz=float(np.sum(wi * c ** 4)), Iss=float(np.sum(wi * s ** 4)), Izs=float(np.sum(wi * c * c * s * s)),
                L=float(Ls.sum()), n=len(segs), ang=ang, ent=ent, c=c, s=s, Ls=Ls, A=A, by_zone=by_zone, zmid=zmid)


# =========================================================================== физика
def shell_freq_donnell(R, Lz, A11, A22, A12, A66, D11, D22, D12, D66, rho_s, m_max=3, n_max=14):
    """наименьшая частота ортотропной цилиндрической оболочки (Доннелл, шарнир по концам), Гц, и моды (m, n)"""
    best = (np.inf, 0, 0)
    det = A11 * A22 - A12 * A12
    for m in range(1, m_max + 1):
        lam = m * math.pi / Lz
        for n in range(2, n_max + 1):
            k = n / R
            Db = D11 * lam ** 4 + 2 * (D12 + 2 * D66) * lam ** 2 * k ** 2 + D22 * k ** 4
            den = A11 * k ** 4 + (det / A66 - 2 * A12) * lam ** 2 * k ** 2 + A22 * lam ** 4
            Mb = det * lam ** 4 / (R * R * den)
            w2 = (Db + Mb) / rho_s
            f = math.sqrt(max(w2, 0.0)) / (2 * math.pi)
            if f < best[0]:
                best = (f, m, n)
    return best


def strut_frame(N, L0, lean_deg, bow, cone_deg, chord, tc, E, alpha, T_case, T_hub, T_strut, r_case, r_hub, d_service=0.0, rho=4430.0):
    """рама подшипника из N стоек-профилей (NACA 00tc, площадь 0,685·c·t, I_min ≈ 0,0393·c·t³ — оценка формы профиля).
    Изгиб (bow) — прогиб e = bow·L: осевая податливость L/EA + (8/15)·e²L/EI (Кастильяно, параболическая ось).
    Навал (lean) — поворот ступицы снимает радиальное тепловое рассогласование (энергетический минимум по повороту и осевому
    сдвигу ступицы). Возврат: боковая жёсткость, напряжения стоек, загромождение, масса."""
    psi, chi = math.radians(lean_deg), math.radians(cone_deg)
    L = L0 / max(math.cos(psi) * math.cos(chi), 0.2)
    t = tc * chord
    A = 0.685 * chord * t
    I = 0.0393 * chord * t ** 3
    e = bow * L
    k_a = 1.0 / (L / (E * A) + (8.0 / 15.0) * e * e * L / (E * I))
    k_t = 12.0 * E * I / L ** 3
    a = np.array([math.cos(psi) * math.cos(chi), math.sin(psi) * math.cos(chi), math.sin(chi)])
    K1 = k_a * np.outer(a, a) + k_t * (np.eye(3) - np.outer(a, a))
    Klat = 0.0
    for j in range(N):
        th = 2 * math.pi * j / N
        # проекция глобального x на локальные оси стойки (радиус, окружность, ось)
        ex_loc = np.array([math.cos(th), -math.sin(th), 0.0])
        Klat += ex_loc @ K1 @ ex_loc
    # тепловое рассогласование: корпус уходит от ступицы на δr (> 0 — стойка растягивается)
    dr = alpha * ((T_case - 20.0) * r_case - (T_hub - 20.0) * r_hub - (T_strut - 20.0) * L0)
    # ступица может повернуться (φ) и сдвинуться по оси (w): минимум энергии одной стойки (все одинаковы)
    #   u = (−δr, r_hub·φ, w) в локальных осях ступицы относительно корпуса
    M = K1[1:, 1:]
    b = -K1[1:, 0] * (-dr)
    if lean_deg > 0.5 or cone_deg > 0.5:
        try:
            uv = np.linalg.solve(M + 1e-12 * np.eye(2), b)
        except np.linalg.LinAlgError:
            uv = np.zeros(2)
    else:
        uv = np.zeros(2)
    u = np.array([-dr, uv[0], uv[1]])
    P_ax = k_a * (u @ a)
    V_tr = k_t * np.linalg.norm(u - (u @ a) * a)
    c_half = 0.5 * t
    s_ax = abs(P_ax) / A + abs(P_ax) * e * c_half / I
    s_b = (V_tr * L / 2.0) * c_half / I
    # нагрузка подшипника: боковая сила распределяется по стойкам пропорционально жёсткости (оценка)
    s_load = F_BEAR / max(N / 2.0, 1.0) / A * (1.0 + e * c_half * A / I)
    sig = s_ax + s_b + s_load
    blockage = N * t * (r_case - r_hub) / (math.pi * (r_case ** 2 - r_hub ** 2))
    t_serv = max(t, d_service + 2 * SKIN_CHANNEL)
    mass = rho * L * (A * (N - 1) + 0.685 * chord * t_serv)
    return dict(K_lat=float(Klat), sigma=float(sig), sig_thermal=float(s_ax + s_b), sig_load=float(s_load), k_axial=float(k_a), k_trans=float(k_t),
                dr_th=float(dr), blockage=float(blockage), mass=float(mass), L=float(L), t=float(t), t_service=float(t_serv),
                service_ok=bool(t_serv <= 0.35 * chord + 1e-9))


def hot_cage(p):
    """беличья клетка горячей части: N пальцев толщины t (радиально), ширины b (по окружности), осевой пролёт L_span, форма
    s ≥ 1 — длина развёртки L = s·L_span (изогнутые/винтовые пальцы). Изгиб с заделкой и направляющей: σ = 3Etδ/L²,
    k_r = Ebt³/L³, k_θ = Etb³/L³ / s² (изогнутый палец мягче в плоскости — оценка)."""
    T = 600.0                                          # °C, температура пальцев (оценка, между корпусом и горячей частью)
    mp = mat_props("in718p", T)
    N, t, b, s, Ls = p["hf_n"], p["hf_t"], p["hf_b"], p["hf_shape"], p["hf_span"]
    L = s * Ls
    E = mp["E"]
    kr = E * b * t ** 3 / L ** 3
    kth = E * t * b ** 3 / L ** 3 / (s * s)
    K = N / 2.0 * (kr + kth)
    d = INP["dr_hot"]
    sig_th = 3.0 * E * t * d / L ** 2
    Vsh = INP["m_hot"] * 9.81 * SHOCK_G / max(N / 2.0, 1.0)
    sig_w = (Vsh * L / 2.0) * (0.5 * b) / (t * b ** 3 / 12.0)
    allow = 0.5 * mp["su"]                            # малоцикловая усталость, половина σв (оценка)
    gap_ok = t * 3 <= 0.03                             # три толщины в радиальном зазоре 0,03 м (оценка)
    mass = N * L * b * t * mp["rho"]
    return dict(K=float(K), sigma=float(sig_th + sig_w), sig_thermal=float(sig_th), sig_shock=float(sig_w), allow=float(allow), L=float(L),
                mass=float(mass), gap_ok=bool(gap_ok), E=float(E))


def beam_modes(EI, mz, z, mounts_z, K_f, K_r, K_hot, z_hot, n_modes=4):
    """изгибные частоты системы в вертикальной плоскости: корпус — балка Эйлера — Бернулли (элементы по сетке z),
    опоры — пружины K_MOUNT; ротор — жёсткое тело (смещение, поворот) на пружинах рам K_f, K_r; горячая часть —
    точечная масса на клетке K_hot. Возврат: частоты, Гц (по возрастанию)."""
    from scipy.linalg import eigh
    n = len(z)
    nd = 2 * n + 3
    K = np.zeros((nd, nd)); Mm = np.zeros((nd, nd))
    for e in range(n - 1):
        le = z[e + 1] - z[e]
        ei = 0.5 * (EI[e] + EI[e + 1]); me = 0.5 * (mz[e] + mz[e + 1])
        ke = ei / le ** 3 * np.array([[12, 6 * le, -12, 6 * le], [6 * le, 4 * le * le, -6 * le, 2 * le * le],
                                      [-12, -6 * le, 12, -6 * le], [6 * le, 2 * le * le, -6 * le, 4 * le * le]])
        mloc = me * le / 420.0 * np.array([[156, 22 * le, 54, -13 * le], [22 * le, 4 * le * le, 13 * le, -3 * le * le],
                                           [54, 13 * le, 156, -22 * le], [-13 * le, -3 * le * le, -22 * le, 4 * le * le]])
        idx = [2 * e, 2 * e + 1, 2 * e + 2, 2 * e + 3]
        K[np.ix_(idx, idx)] += ke; Mm[np.ix_(idx, idx)] += mloc

    def node(zz):
        return int(np.argmin(np.abs(z - zz)))
    for zz, k in mounts_z:
        i = node(zz); K[2 * i, 2 * i] += k
    iy, ip, ih = 2 * n, 2 * n + 1, 2 * n + 2
    mr, It, zc = INP["m_rotor"], INP["m_rotor"] * INP["rotor_L"] ** 2 / 12.0, INP["rotor_cg"]
    Mm[iy, iy] += mr; Mm[ip, ip] += It
    for zb, kb in ((INP["z_fb"], K_f), (INP["z_rb"], K_r)):
        i = node(zb)
        # относительное смещение: y_r + (zb − zc)·φ − w_i
        v = np.zeros(nd); v[iy] = 1.0; v[ip] = zb - zc; v[2 * i] = -1.0
        K += kb * np.outer(v, v)
    i = node(z_hot)
    Mm[ih, ih] += INP["m_hot"]
    v = np.zeros(nd); v[ih] = 1.0; v[2 * i] = -1.0
    K += K_hot * np.outer(v, v)
    w2 = eigh(K, Mm, eigvals_only=True, subset_by_index=[0, n_modes - 1])
    return np.sqrt(np.maximum(w2, 0.0)) / (2 * math.pi)


# =========================================================================== каналы
CHANNELS = [  # (система, d, м, длина-оценка, назначение) — порядок встраивания по гену emb
    ("топливо", None, None, "дерево Мюррея: подвод → 14 форсунок (k = 2,41 — manifold_results)"),
    ("масло, передняя опора", 0.022, 0.30, "подвод Ø12 и слив Ø22 мм в служебной стойке (shell_topo RADIAL_CH)"),
    ("масло, задняя опора", 0.022, 0.45, "то же у задней рамы"),
    ("отбор воздуха", 0.020, 0.50, "2 × Ø20 мм к охлаждению диска турбины (shell_topo AXIAL_CH)"),
    ("датчики", 0.016, 0.60, "2 × Ø16 мм проводка термопар и давления (shell_topo AXIAL_CH)"),
]


def fuel_tree(p, prof, n_pts=2):
    n = INP["fuel"]["n_inj"]
    z_inj = INP["zones"][3]["z0"] + 0.25 * (INP["zones"][3]["z1"] - INP["zones"][3]["z0"])
    R = float(np.interp(z_inj, prof["z"], prof["r_m"]))
    C = 2 * math.pi * R
    port_s = (0.75 + 0.5 * (p["fuel_port"] - 0.5)) * C       # подвод — снизу ± четверть окружности
    port_z = INP["zones"][3]["z1"] - 0.02
    ss = (np.arange(n) + 0.5) / n * C
    # развёртка с переносом через шов к ближайшей стороне подвода
    ss = np.where(ss - port_s > C / 2, ss - C, np.where(port_s - ss > C / 2, ss + C, ss))
    T = np.stack([np.full(n, z_inj), ss], 1)
    order = np.argsort(ss)
    segs = S.murray_tree((port_z, port_s), T, INP["fuel"]["d_leaf"], k=p["murray_k"] if "murray_k" in p else INP["fuel"]["k"],
                         order=order, n_pts=n_pts)
    Ltot = sum(float(np.sum(np.linalg.norm(np.diff(sg["pts"], axis=0), axis=1))) for sg in segs)
    dmax = max(sg["d"] for sg in segs)
    return segs, dict(L=Ltot, d_max=dmax, R=R, z_inj=z_inj)


# =========================================================================== оценка
def evaluate_params(p, variant, with_geom=False):
    prof = profiles(p, variant)
    z, r_m, t, r_i, r_o = prof["z"], prof["r_m"], prof["t"], prof["r_i"], prof["r_o"]
    dz = np.gradient(z)
    route = p["route"]
    shell_mat = "ti64" if route == "печать" else p["cast_mat"]
    zi = _zone_of(z)
    T_z = np.array([INP["zones"][k]["T"] for k in zi])
    T_mean = float(np.mean(T_z))
    mp = mat_props(shell_mat, T_mean)
    rho_m = mp["rho"]
    # ---------------- масса стенки
    t_s = np.where(prof["core"] > 0, np.minimum(p["skin"], 0.5 * t), 0.5 * t)
    t_c = prof["t_core"]
    t_solid_eq = t - t_c + prof["rho"] * t_c                       # металл в стенке на единицу площади / ρ
    m_wall = float(np.sum(2 * math.pi * r_m * t_solid_eq * dz) * rho_m)
    V_wall = float(np.sum(2 * math.pi * r_m * t * dz))
    V_core = float(np.sum(2 * math.pi * r_m * t_c * dz))
    lat_frac = V_core / V_wall
    area = float(np.sum(2 * math.pi * r_o * dz))
    R_ref = float(np.sum(r_o * dz) / np.sum(dz))
    mounts, minfo = mounts_of(p, R_ref)
    segs = rib_segments(p, R_ref, mounts)
    rs = rib_stats(segs, p, R_ref, area)
    m_ribs = rs["V"] * rho_m
    # ---------------- рамы
    E_case = mp["E"]
    rA = float(np.interp(INP["z_fb"], z, r_i)) if variant == "оболочка" else INP["zones"][0]["r_req"]
    rE = float(np.interp(INP["z_rb"], z, r_i)) if variant == "оболочка" else INP["zones"][4]["r_req"]
    emb = p["emb"]
    d_serv_f = 0.022 if emb >= 2 else 0.0
    d_serv_r = 0.022 if emb >= 3 else 0.0
    mf = mat_props("ti64" if route != "литьё" else shell_mat, INP["zones"][0]["T"])
    fr_f = strut_frame(p["f_n"], rA - R_HUB, p["f_lean"], p["f_bow"], p["f_cone"], p["f_chord"], p["f_tc"], mf["E"], mf["alpha"],
                       INP["zones"][0]["T"], T_HUB[0], INP["zones"][0]["T"], rA, R_HUB, d_serv_f, mf["rho"])
    mr = mat_props("in718p", T_STRUT_REAR)
    fr_r = strut_frame(p["r_n"], rE - R_HUB, p["r_lean"], p["r_bow"], p["r_cone"], p["r_chord"], p["r_tc"], mr["E"], mr["alpha"],
                       INP["zones"][4]["T"], T_HUB[1], T_STRUT_REAR, rE, R_HUB, d_serv_r, mr["rho"])
    sy_f = mf["sy"] * K_ROUTE["печать"] / SF_YIELD
    sy_r = mr["sy"] * K_ROUTE["печать"] / SF_YIELD
    cage = hot_cage(p)
    # ---------------- каналы
    fsegs, finfo = fuel_tree(p, prof)
    emb_list = [CHANNELS[k][0] for k in range(min(emb, 5))]
    t_wall_at = lambda zz: float(np.interp(zz, z, t))  # noqa: E731
    boss_m = 0.0
    ch_mass = 0.0
    pipes_ext = 0
    if "топливо" in emb_list:
        for sg in fsegs:
            Ls = float(np.sum(np.linalg.norm(np.diff(sg["pts"], axis=0), axis=1)))
            need = sg["d"] + 2 * SKIN_CHANNEL
            zz = float(sg["pts"][len(sg["pts"]) // 2, 0])
            avail = t_wall_at(zz) + (0.6 * p["rib_h"] if rs["L"] > 3.0 else 0.0)   # рёбра частой сетки несут канал (оценка)
            boss_m += max(0.0, need - avail) * need * Ls * rho_m
            ch_mass -= math.pi / 4 * sg["d"] ** 2 * Ls * rho_m * (need <= avail)
    else:
        pipes_ext += 1 + INP["fuel"]["n_inj"]                     # кольцо-коллектор + трубки к форсункам
        ch_mass += (finfo["L"] * 0.6) * math.pi * 0.012 * 0.0015 * 7850 + 0.5 * (1 + INP["fuel"]["n_inj"])
    for (nm, d, Lc, _) in CHANNELS[1:]:
        if nm in emb_list:
            ch_mass += math.pi * ((d / 2 + SKIN_CHANNEL) ** 2 - (d / 2) ** 2) * Lc * rho_m
        else:
            pipes_ext += 2
            ch_mass += 2 * (Lc * math.pi * d * 0.0015 * 7850 + 0.5)
    # ---------------- опоры (бобышки цапф и тяги: 3 × 6 кг — оценка) и итоговая масса
    m_mounts = 18.0
    mass = m_wall + m_ribs + fr_f["mass"] + fr_r["mass"] + cage["mass"] + boss_m + max(ch_mass, -0.5 * m_wall) + m_mounts
    m_cart_case = 0.0
    if variant == "картридж":
        # корпус картриджа (статор внутри) — тонкая стальная/титановая обечайка 8 мм по требуемому радиусу (оценка)
        m_cart_case = float(np.sum(2 * math.pi * R_REQ * 0.008 * dz) * 4430.0)
    # ---------------- прочность стенки по зонам
    E, nu, alpha = mp["E"], mp["nu"], mp["alpha"]
    Er, Sr = lattice_rel(float(np.mean(prof["rho"])))
    route_k = K_ROUTE["печать"] if route == "печать" else (K_ROUTE["литьё_ti"] if shell_mat == "ti64" else K_ROUTE["литьё_сталь"])
    u_wall, u_core, worst = 0.0, 0.0, None
    sig_th_max = 0.0
    zone_rows = []
    rib_bz = rs["by_zone"]
    for k, zn in enumerate(INP["zones"]):
        m = zi == k
        R = float(np.mean(r_m[m])); tt = float(np.mean(t[m])); tc = float(np.mean(t_c[m])); ts = tt - tc
        Lz = zn["z1"] - zn["z0"]
        circA = 2 * math.pi * R * Lz
        a_hoop = rib_bz[k, 1] / circA            # размазанная «толщина» рёбер по окружности (вклад в жёсткость s⁴)
        a_ax = rib_bz[k, 0] / circA
        t_hoop = ts + Er * tc + a_hoop
        t_ax = ts + Er * tc + a_ax
        sth = zn["p_g"] * R / t_hoop
        sz = zn["p_g"] * R / (2 * t_ax)
        svm = math.sqrt(sth ** 2 - sth * sz + sz ** 2)
        mpz = mat_props(shell_mat, zn["T"])
        allow = route_k * mpz["sy"] / SF_YIELD
        # перепад температур к соседней зоне на длине перехода ℓ = 0,5·min(L зон) (оценка): решение для оболочки с изломами
        # температуры (функция Грина балки на упругом основании) σ ≈ 0,40·EαΔT·√(Rt)/ℓ, не больше, чем у ступеньки (0,55·EαΔT) — README
        s_th = 0.0
        for kk in (k - 1, k + 1):
            if 0 <= kk < 5:
                dT = abs(zn["T"] - INP["zones"][kk]["T"])
                ell = 0.5 * min(Lz, INP["zones"][kk]["z1"] - INP["zones"][kk]["z0"])
                s_th = max(s_th, mpz["E"] * mpz["alpha"] * dT * min(0.55, 0.40 * math.sqrt(R * tt) / ell))
        Kt = 1.0 + (C_KT * math.sqrt(min(p["rib_h"], tt) / p["rib_fillet"]) if (rib_bz[k, 2] > 0) else 0.0)
        speak = Kt * svm + s_th
        u = speak / allow
        uc = (Er / max(Sr, 1e-9)) * (sth / mpz["E"]) * mpz["E"] / (mpz["sy"] * route_k) if tc > 1e-6 else 0.0
        zone_rows.append(dict(zone=zn["name"], R=R, t=tt, t_core=tc, sig_hoop=sth, sig_vm=svm, sig_th=s_th, Kt=Kt, u=u, u_core=uc, allow=allow))
        if u > u_wall:
            u_wall, worst = u, zn["name"]
        u_core = max(u_core, uc)
        sig_th_max = max(sig_th_max, s_th)
    # ---------------- частоты: оболочка по зонам
    f_shell = np.inf; f_shell_mode = None
    for k, zn in enumerate(INP["zones"]):
        row = zone_rows[k]
        R, tt, tc = row["R"], row["t"], row["t_core"]
        ts = tt - tc
        Lz = zn["z1"] - zn["z0"]
        Am = E / (1 - nu * nu) * (ts + Er * tc)
        if tc > 1e-6:
            tsk = 0.5 * ts
            Iw = 2 * (tsk ** 3 / 12 + tsk * (0.5 * (tt - tsk)) ** 2) + Er * tc ** 3 / 12
        else:
            Iw = tt ** 3 / 12
        Dw = E / (1 - nu * nu) * Iw
        circA = 2 * math.pi * R * Lz
        # рёбра зоны размазаны (Baruch & Singer — приём; эксцентриситет — переносом оси, оценка)
        Azz_r, Ass_r, Azs_r = rib_bz[k, 0] / circA, rib_bz[k, 1] / circA, rib_bz[k, 3] / circA
        Izz_r, Iss_r, Izs_r = rib_bz[k, 4] / circA, rib_bz[k, 5] / circA, rib_bz[k, 6] / circA
        A11 = Am + E * Azz_r; A22 = Am + E * Ass_r
        A12 = nu * Am + E * Azs_r; A66 = E / (2 * (1 + nu)) * (ts + Er * tc) + E * Azs_r
        ecc2 = (0.5 * tt + 0.5 * p["rib_h"]) ** 2
        D11 = Dw + E * (Izz_r + Azz_r * ecc2); D22 = Dw + E * (Iss_r + Ass_r * ecc2)
        D12 = nu * Dw + E * (Izs_r + Azs_r * ecc2); D66 = Dw * (1 - nu) / 2 + E * (Izs_r + Azs_r * ecc2)
        m_ns = (M_INTERNAL / 2.0) / circA if k in (1, 2) else 0.0
        rho_s = rho_m * (ts + prof["rho"].mean() * tc) + rho_m * rib_bz[k, 2] / max(circA, 1e-12) + m_ns
        f, mm, nn = shell_freq_donnell(R, Lz, A11, A22, A12, A66, D11, D22, D12, D66, rho_s)
        if f < f_shell:
            f_shell, f_shell_mode = f, dict(zone=zn["name"], m=mm, n=nn)
    # ---------------- частоты: балочная модель системы
    zb = np.linspace(0.0, z[-1], 41)
    rmb = np.interp(zb, z, r_m); tb = np.interp(zb, z, t - t_c + Er * t_c)
    axr = np.zeros_like(zb)
    for k in range(5):
        circA = 2 * math.pi * zone_rows[k]["R"] * (INP["zones"][k]["z1"] - INP["zones"][k]["z0"])
        mzk = (zb >= INP["zones"][k]["z0"]) & (zb <= INP["zones"][k]["z1"])
        axr[mzk] = rib_bz[k, 0] / max(circA, 1e-12)
    EI = E * math.pi * rmb ** 3 * (tb + axr)
    mz_lin = 2 * math.pi * rmb * rho_m * np.interp(zb, z, t_solid_eq) + rs["V"] * rho_m / z[-1]
    mz_lin = mz_lin + np.where((zb > INP["zones"][1]["z0"]) & (zb < INP["zones"][2]["z1"]), M_INTERNAL / (INP["zones"][2]["z1"] - INP["zones"][1]["z0"]), 0.0)
    if variant == "картридж":
        mz_lin = mz_lin + 2 * math.pi * np.interp(zb, z, R_REQ) * 0.008 * 4430.0
        Kf = 1.0 / (1.0 / fr_f["K_lat"] + 1.0 / K_SPIGOT)
        Kr = 1.0 / (1.0 / fr_r["K_lat"] + 1.0 / K_SPIGOT)
        Kh = 1.0 / (1.0 / cage["K"] + 1.0 / K_SPIGOT)
    else:
        Kf, Kr, Kh = fr_f["K_lat"], fr_r["K_lat"], cage["K"]
    mz_lin = mz_lin + (fr_f["mass"] + fr_r["mass"] + m_mounts) / z[-1]
    fb = beam_modes(EI, mz_lin, zb, [(minfo["z_side"], 2 * K_MOUNT), (minfo["z_rear"], K_MOUNT)], Kf, Kr, Kh,
                    INP["zones"][4]["z0"] + 0.5 * p["hf_span"])
    f_need = F_MARGIN * INP["f_op"]
    f_min = float(min(f_shell, fb[0]))
    # ---------------- детали, крепёж, ремонт
    D_front, D_rear = 2 * float(r_o[0]) + 0.04, 2 * float(r_o[-1]) + 0.04      # фланцы концов (+0,04 м на бурт — оценка)
    split_D = None
    if variant == "картридж":
        zsplit = 0.0 if p["extract"] == "вперёд" else z[-1]
        split_D = 2 * float(np.interp(zsplit, z, r_o)) + 0.04
    method = "печать" if route == "печать" else "литьё"
    D_env = 2 * float(np.max(r_o)) + 2 * p["rib_h"]
    L_env = float(z[-1])
    bom = dict(PM.REF_BOM)
    bom.update(case=dict(method=method, D=D_env, pieces=1), case_arch=variant, hot_integral=True,
               comb=dict(kind="кольцевая", method="печать", n_inj=INP["fuel"]["n_inj"]), conn="быстроразъём")
    pieces, printer = 1, None
    if route == "печать":
        lp = [q for q in MK.printers_for(D_env, L_env, False) if q["proc"] == "LPBF"]
        if lp:
            printer = lp[0]["name"]
        else:
            pieces = max(2, math.ceil(max(D_env, L_env) / MK.largest_print_dim()))
        bom["case"]["pieces"] = pieces
    pm = PM.count(bom)
    d_parts, d_fast, d_pipe_j, rules = 0, 0, 0, []

    def per_joint(D):
        return 2 if D <= PM.CLAMP_D_MAX else math.ceil(math.pi * D / PM.BOLT_PITCH)
    # правило семейства: стыки на концах — по диаметру своего фланца (сужение силуэта к концам даёт хомут), прочие — по корпусу
    f_pm = pm["by_unit"]["крепёж фланцев"]["fasteners"]
    n_j = int(round(f_pm / per_joint(D_env)))
    Dl = ([D_front, D_rear] + [D_env] * max(0, n_j - 2)) if n_j >= 2 else [min(D_front, D_rear)] * n_j
    d_fast += sum(per_joint(D) for D in Dl) - f_pm
    rules.append(f"стыков корпуса {n_j}: концы Ø{D_front:.2f} / Ø{D_rear:.2f} м — {per_joint(D_front)} / {per_joint(D_rear)} ед. крепежа, "
                 f"прочие Ø{D_env:.2f} м — {per_joint(D_env)} (правило parts_model: хомут до {PM.CLAMP_D_MAX} м)")
    if variant == "оболочка":
        not_emb_oil = sum(1 for nm in ("масло, передняя опора", "масло, задняя опора") if nm not in emb_list)
        if "топливо" not in emb_list:
            d_parts += 2; d_pipe_j += 4; rules.append("топливо не встроено: кольцо и подвод трубками (+2 детали, +4 соединения)")
        if not_emb_oil:
            d_parts += 2 * not_emb_oil; d_pipe_j += 4 * not_emb_oil
            rules.append(f"масло не встроено у {not_emb_oil} опор: подвод и слив трубками (+2 детали и +4 соединения на опору)")
    else:
        emb_oil = sum(1 for nm in ("масло, передняя опора", "масло, задняя опора") if nm in emb_list)
        if emb_oil:
            d_parts += -2 * emb_oil + emb_oil; d_pipe_j += -4 * emb_oil
            rules.append(f"картридж: масло встроено у {emb_oil} опор — минус 2 трубки, плюс уплотнительное кольцо перехода на опору")
        if "топливо" in emb_list:
            d_parts += -1; d_pipe_j += -2
            rules.append("картридж: топливо встроено — минус подводящая трубка (рукав пакета остаётся)")
        d_parts += 2
        d_fast += 2 if split_D <= PM.CLAMP_D_MAX else math.ceil(math.pi * split_D / PM.BOLT_PITCH)
        rules.append("картридж: наружный корпус и съёмная крышка (+2 детали), стык крышки: хомут или болты по D (правило parts_model)")
    for nm in ("отбор воздуха", "датчики"):
        if nm not in emb_list:
            d_parts += 2; d_pipe_j += 4
    if route == "гибрид":
        Df = 2 * rA + 0.1; Dr = 2 * rE + 0.1
        d_parts += 2
        d_fast += math.ceil(math.pi * Df / PM.BOLT_PITCH) + math.ceil(math.pi * Dr / PM.BOLT_PITCH)
        rules.append("гибрид: печатные рамы подшипников вставками в литой корпус на болтах (+2 детали)")
    if pieces > 1:
        d_parts += pieces - 1
        d_fast += (pieces - 1) * math.ceil(math.pi * D_env / PM.BOLT_PITCH)
        rules.append(f"оболочка не входит в камеру LPBF: {pieces} части на болтовых стыках")
    parts = pm["parts"] + d_parts
    fasteners = pm["fasteners"] + d_fast
    n_eff = pm["n_eff"] + d_parts + PM.K_FASTENER * d_fast + 0.5 * d_pipe_j
    # ремонт (service_model.score — одна точка правды воронки)
    m_eng = mass + INP["m_rotor"] + INP["m_hot"] + M_INTERNAL
    ext_pipes = pipes_ext + 3
    joints_ext = [dict(type="хомут") if 2 * float(r_o[0]) <= PM.CLAMP_D_MAX else dict(type="болты", n_bolts=math.ceil(math.pi * 2 * float(r_o[0]) / PM.BOLT_PITCH)),
                  dict(type="хомут") if 2 * float(r_o[-1]) <= PM.CLAMP_D_MAX else dict(type="болты", n_bolts=math.ceil(math.pi * 2 * float(r_o[-1]) / PM.BOLT_PITCH))]
    if variant == "оболочка":
        mods = [dict(name="двигатель целиком", mass_kg=m_eng, joints=joints_ext, pipes=ext_pipes, cables=3, cables_keyed=False, align=True,
                     frequent=True, repair=SM.REPAIR_FACTORY)]
    else:
        m_cart = INP["m_rotor"] + INP["m_hot"] + M_INTERNAL + m_cart_case + fr_f["mass"] + fr_r["mass"] + cage["mass"]
        jt = dict(type="хомут") if split_D <= PM.CLAMP_D_MAX else dict(type="болты", n_bolts=math.ceil(math.pi * split_D / PM.BOLT_PITCH))
        mods = [dict(name="картридж газогенератора", mass_kg=m_cart, joints=[jt], pipes=max(0, ext_pipes - 3), cables=3, cables_keyed=False,
                     align=True, frequent=True, repair=SM.REPAIR_SWAP),
                dict(name="наружный корпус", mass_kg=mass, joints=joints_ext, pipes=3, cables=0, align=False, frequent=False)]
    svc = SM.score(dict(modules=mods, amb=False, single_connector=False))
    # ---------------- изготовимость
    cons = {}
    q = 1.0
    notes = []
    tmin = float(t.min())
    if route == "печать":
        sh = MK.shell_process(D_env)
        q *= MK.MAT["Ti64"]["print"] * sh["proc_q"]
        if not sh["hip"]:
            cons["hip"] = 1.0
        notes.append(f"ГИП: камера {sh['vessel']} м" if sh["hip"] else "ГИП невозможен")
        if lat_frac > 0.01:
            cons["skin"] = max(0.0, (PRINT_SKIN_MIN - p["skin"]) / PRINT_SKIN_MIN)
            tl = float(S.tpms_t_for_density("gyroid", float(np.mean(prof["rho"]))))
            pore = p["cell"] * max(0.5 - tl, 0.05) if p["lat_type"] != "voronoi" else p["cell"] * (1 - float(np.mean(prof["rho"]))) * 0.6
            cons["pore"] = max(0.0, (PORE_MIN - pore) / PORE_MIN)
            if pore < PORE_LONG:
                q *= 0.6; notes.append(f"поры решётки {pore * 1e3:.1f} мм < 5 мм при длинном пути выноса порошка")
        else:
            pore = None
        # нависание: рёбра, наклонённые к плите круче 45° (ось печати — ось двигателя)
        if rs["L"] > 0:
            ov = float(np.sum(rs["Ls"] * (np.abs(rs["s"]) > math.cos(math.radians(45.0))))) / rs["L"]
        else:
            ov = 0.0
        flare = np.gradient(r_o, z)
        sgn = 1.0 if p["build_dir"] == "+z" else -1.0
        ov_sil = float(np.sum(dz * (sgn * flare > 1.0)) / z[-1])
        ov_tot = 0.7 * ov + 0.3 * ov_sil
        q *= (1.0 - 0.5 * ov_tot)
        sv_cone = min(p["f_cone"], p["r_cone"])
        if sv_cone < 45.0:
            q *= 0.9; notes.append("стойки рам наклонены к плите круче 45°: поддержки (конус < 45°)")
    else:
        ov_tot, pore = 0.0, None
        lim = MK.CAST["static"]
        q *= (MK.MAT["steel"]["cast"] if shell_mat == "zg15" else MK.MAT["Ti64"]["cast"])
        if D_env > lim["max"]:
            q *= 0.4; notes.append("крупнее предела литья 1,5 м: из частей")
        elif D_env > lim["ok"]:
            q *= 0.6
        cons["cast_wall"] = max(0.0, (CAST_T_MIN - tmin) / CAST_T_MIN)
        if p["rib_fillet"] < 0.003:
            q *= 0.8
        if route == "гибрид":
            q *= 0.9
    cons["wall_print"] = max(0.0, (0.004 - tmin) / 0.004)
    cons["u_wall"] = max(0.0, u_wall - 1.0)
    cons["u_core"] = max(0.0, u_core - 1.0)
    cons["u_strut_f"] = max(0.0, fr_f["sigma"] / sy_f - 1.0)
    cons["u_strut_r"] = max(0.0, fr_r["sigma"] / sy_r - 1.0)
    cons["u_cage"] = max(0.0, cage["sigma"] / cage["allow"] - 1.0)
    cons["cage_gap"] = 0.0 if cage["gap_ok"] else 0.5
    cons["freq"] = max(0.0, (f_need - f_min) / f_need)
    cons["blockage"] = max(0.0, fr_f["blockage"] - BLOCKAGE_MAX) + max(0.0, fr_r["blockage"] - BLOCKAGE_MAX)
    cons["service_strut"] = 0.0 if (fr_f["service_ok"] and fr_r["service_ok"]) else 0.3
    CV = float(sum(cons.values()))
    # ---------------- дескрипторы и признаки
    sym = FT.angular_symmetry(rs["ang"]) if rs["V"] > 0 else dict(N_eff=0, mirror=0.0, asym=0.0, axisym=True)
    if sym.get("axisym") or sym["N_eff"] >= 12 or sym["N_eff"] == 0:
        sym_cls = 0
    elif sym["N_eff"] >= 6:
        sym_cls = 1
    elif sym["N_eff"] >= 3:
        sym_cls = 2
    elif sym["N_eff"] == 2 or sym["mirror"] <= 0.08:
        sym_cls = 3
    else:
        sym_cls = 4
    rib_cls = G[GENE_INDEX["rib_type"]]["opts"].index(p["rib_type"]) if rs["V"] > 0 else 0
    Psil = np.stack([z, r_o], 1)
    sm = FT.curve_smoothness(Psil)
    E_env = ENV[(variant, p["extract"] if variant == "картридж" else "вперёд")]
    turn_env = S.turning_total(np.stack([z, E_env + 0.5 * t.mean()], 1))
    curv = max(0.0, sm["turning"] - turn_env)              # кривизна силуэта сверх неизбежной (огибающая)
    desc = np.array([lat_frac, rib_cls, sym_cls, curv])
    feats = np.array([lat_frac, rib_cls / 6.0, sym_cls / 4.0, curv, rs["L"] / area, p["rib_h"], rs["ent"], p["f_lean"] / 40.0,
                      p["f_bow"] / 0.2, p["r_lean"] / 45.0, p["hf_shape"] / 3.0, float(r_o.min() / r_o.max()), p["f_tub"] / 0.1,
                      sm["inflections"] / 10.0])
    objs = dict(mass_kg=mass, n_eff=n_eff, parts=parts, fasteners=fasteners, pipe_joints=d_pipe_j, service=svc["score"], q=q,
                f_min_Hz=f_min, f_shell_Hz=float(f_shell), f_beam_Hz=[float(v) for v in fb], u_wall=u_wall, u_core=u_core,
                u_strut_f=fr_f["sigma"] / sy_f, u_strut_r=fr_r["sigma"] / sy_r, u_cage=cage["sigma"] / cage["allow"], CV=CV,
                lat_frac=lat_frac, D_env=D_env, L_env=L_env, printer=printer, pieces=pieces, overhang=ov_tot, pore_m=pore,
                smooth_fair=sm["fairness"], turning=sm["turning"], curv_excess=curv, inflections=sm["inflections"])
    out = dict(objs=objs, desc=desc, feats=feats, cons=cons, CV=CV)
    if with_geom:
        out.update(prof=prof, segs=segs, mounts=mounts, minfo=minfo, rib=rs, frames=dict(front=fr_f, rear=fr_r), cage=cage,
                   fuel=dict(segs=fsegs, info=finfo), zone_rows=zone_rows, f_shell_mode=f_shell_mode, pm=pm, rules=rules, svc=svc,
                   notes=notes, emb=emb_list, sym=sym, mat=shell_mat, bom_case=bom["case"], conn=bom["conn"], R_ref=R_ref)
    return out


# --------------------------------------------------------------------------- опорная (типовая OEM) форма и фитнес
def reference_params(variant="оболочка", kind="литой с фланцами"):
    """типовые формы изготовителей (для нормировки и новизны): литой стальной корпус с кольцевыми фланцами и прямыми
    радиальными стойками, трубки снаружи; вариант с осевыми рёбрами; гладкая оболочка без рёбер (оценка)"""
    p = decode(np.full(DIM, 0.5))
    p.update({f"off{k}": 0.0 for k in range(8)})
    p.update({f"t{k}": 0.022 for k in range(4)})
    p.update({f"core{k}": 0.0 for k in range(4)})
    p.update(route="литьё", cast_mat="zg15", rib_type="кольца", rib_n=0.0, rib_n2=0.3, rib_h=0.04, rib_w=0.03, rib_fillet=0.006,
             sym_n=1, asym1=0.0, asym2=0.0, asym_phase=0.0, f_n=8, f_lean=0.0, f_bow=0.0, f_cone=0.0, f_chord=0.10, f_tc=0.18,
             f_tub=0.0, r_n=10, r_lean=0.0, r_bow=0.0, r_cone=0.0, r_chord=0.12, r_tc=0.18, hf_n=24, hf_t=0.004, hf_b=0.04,
             hf_shape=1.0, hf_span=0.12, emb=0, murray_k=2.41, build_dir="+z", extract="вперёд", mount_drop=0.0)
    if kind == "литой с осевыми рёбрами":
        p.update(rib_type="ортогрид", rib_n=0.7, sym_n=4, rib_h=0.03)
    elif kind == "гладкая оболочка":
        p.update(rib_type="нет", rib_h=0.0)
    return p


REF_KINDS = ("литой с фланцами", "литой с осевыми рёбрами", "гладкая оболочка")
_REF_CACHE = {}


def reference(variant):
    if variant not in _REF_CACHE:
        rs = [evaluate_params(reference_params(variant, k), variant) for k in REF_KINDS]
        _REF_CACHE[variant] = dict(mass=rs[0]["objs"]["mass_kg"], n_eff=rs[0]["objs"]["n_eff"], feats=np.stack([r["feats"] for r in rs]),
                                   rows=[dict(kind=k, **{kk: r["objs"][kk] for kk in ("mass_kg", "n_eff", "parts", "fasteners", "service", "q",
                                                                                      "f_min_Hz", "u_wall", "CV")}) for k, r in zip(REF_KINDS, rs)])
    return _REF_CACHE[variant]


W_FIT = dict(n_eff=0.35, mass=0.25, service=0.20, q=0.20)   # веса — допущение (порядок приоритетов владельца); фронт — NSGA-II


def fitness(objs, variant):
    ref = reference(variant)
    f = -(W_FIT["n_eff"] * objs["n_eff"] / ref["n_eff"] + W_FIT["mass"] * objs["mass_kg"] / ref["mass"]
          + W_FIT["service"] * (1.0 - objs["service"] / 10.0) + W_FIT["q"] * (1.0 - objs["q"]))
    if objs["CV"] > 0:
        f = -10.0 - objs["CV"]
    return float(f)


def evaluate(x, variant):
    p = decode(x)
    try:
        r = evaluate_params(p, variant)
    except (FloatingPointError, ValueError, np.linalg.LinAlgError, ZeroDivisionError) as ex:
        return dict(fitness=float("nan"), desc=None, objs=dict(error=repr(ex)[:200]))
    return dict(fitness=fitness(r["objs"], variant), desc=r["desc"], objs=r["objs"], meta=dict(feats=r["feats"].tolist()))


DESC_NAMES = ["доля решётки", "топология рёбер", "симметрия", "кривизна силуэта"]
DESC_EDGES = [[0.0, 0.02, 0.15, 0.35, 0.6, 1.0001],
              [-0.5, 0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5],
              [-0.5, 0.5, 1.5, 2.5, 3.5, 4.5],
              [0.0, 0.25, 0.6, 1.1, 1.8, 50.0]]
DESC_LABELS = [["нет", "мало", "треть", "половина", "много"],
               list(G[GENE_INDEX["rib_type"]]["opts"]),
               ["осесимм.", "N 6–11", "N 3–5", "N 2 / зеркало", "асимм."],
               ["прямой", "мягкий", "волна", "сильная", "органика"]]
