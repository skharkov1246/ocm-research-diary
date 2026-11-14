"""Воронка камеры сгорания KV, уровень 0: сети реакторов Cantera по схемам и модели модуля.

Для варианта (comb_space.concept + derive) считает на газе (метан, GRI-Mech 3.0) при 100 / 50 / 30 % нагрузки и на
дизеле (суррогат н-додекан) при 100 / 50 %:
  * сеть реакторов схемы: смешение (неравномерность - узлы Гаусса - Эрмита по phi) -> зона пламени (реактор идеального
    смешения, PSR; рециркуляция продуктов - закрутка и FLOX) -> догорание (реактор вытеснения, PFR) -> разбавление ->
    PFR; MILD/FLOX - рециркуляция по уносу струи; RQL - богатая PSR, гашение тремя PSR, бедная PFR; захваченный вихрь -
    полость PSR + основная PSR с горячим газом полости; пассивное гидравлическое деление - пилот так же; каталитическая -
    подогрев до зажигания катализатора, превращение части топлива на катализаторе до предела температуры, гомогенное
    догорание PFR; детонационная - сгорание при постоянном объёме вместо волны (оценка), расширение, PFR;
  * NOx, CO, UHC при 15 % O2 (сухой газ), полнота сгорания, запас по срыву (время пребывания PSR / наименьшее время до
    погасания, бисекция), проскок (критический градиент у стенки по S_L из таблицы), самовоспламенение в смесителе
    (задержка из таблицы / время в смесителе), потери давления, температура стенки (одномерный баланс в форме
    Lefebvre & Ballal), неравномерность на выходе (корреляция в форме Lefebvre + унаследованная неравномерность смешения);
  * модели модуля: число деталей, замена картриджа «хоть уборщиком», Plug&Play, изготовимость; соответствие целям
    targets.json (жёсткие и ранжирующие).

Дизель: суррогат н-додекан (одно вещество, как у ECN Spray A - по памяти, не сверено), механизм Wang, Ra, Jia, Reitz
(Fuel 136 (2014) 25-36, doi:10.1016/j.fuel.2014.07.028; заголовок nDodecane_Reitz.yaml) - азотной химии в нём нет:
равновесие продуктов в механизме Reitz переносится на GRI-Mech 3.0 по общим веществам, NOx - в PFR GRI (термический и
N2O-путь; быстрый NOx фронта пламени не учтён - оценка снизу). Свойства жидкости для распыла - типовой дизель (оценка).

Границы (точки нагрузки, таблицы S_L и задержек) - boundary(): считаются один раз и пишутся в comb_bc_results.json.
Все модели здесь - уровень облика: «оценка», если не указан источник. AWS не вызывается.
"""
import hashlib, inspect, json, math, os, sys, time
from multiprocessing import Pool

import cantera as ct
import numpy as np
from scipy.optimize import brentq

HERE = os.path.dirname(os.path.abspath(__file__))
GTU = os.path.dirname(HERE)
if HERE not in sys.path:
    sys.path.insert(0, HERE)
import comb_space as CS  # noqa: E402

DATE = "2026-11-14"
BC_FILE = os.path.join(HERE, "comb_bc_results.json")
TARGETS = json.load(open(os.path.join(HERE, "targets.json"), encoding="utf-8"))
TB = json.load(open(os.path.join(HERE, "targets_basis_results.json"), encoding="utf-8"))

AIR = "O2:0.2095,N2:0.7808,AR:0.0093,CO2:0.0004"     # как gtu_cycle.AIR
AIR_REITZ = "o2:0.2095,n2:0.7905"                     # в механизме Reitz нет Ar и CO2: Ar отнесён к N2 (как targets_basis)
DZ_MECH = ("nDodecane_Reitz.yaml", "nDodecane_IG")
X_LEAK, DP_COMB, P_EXH = 0.015, 0.04, 1.02            # как frame3_point.py и targets_basis.py
LOADS = [1.0, 0.5, 0.3]                               # доля мощности на валу
DIESEL_LOADS = [0, 1]                                 # индексы LOADS для дизеля (100 и 50 %)
O2_REF, O2_AIR = 15.0, 20.9
H2_GRID = [0.0, 0.3, 0.5, 1.0]
TAD_GRID = [1600.0, 1800.0, 2000.0, 2200.0]
SIGMA = 5.670374419e-8

# --------------------------------------------------------------------------- константы моделей (оценка, если без источника)
DIESEL_LIQ = dict(rho=840.0, sigma=0.028, mu=3.0e-3, T_b=550.0, L_v=2.5e5,
                  note="типовой дизель: вязкость 1,9-4,1 мм2/с при 40 °C (ASTM D975, по памяти); остальное - оценка")
DP_FUEL_PS = 3.0e6         # перепад на центробежной форсунке при 100 %, Па; при частичной ~ расход² (простая форсунка)
K_MIX_D, T_CHEM_K = 1.0, 1.0   # время смешения паров (K·l/u, один пролёт) и множитель химического времени α/S_L²
H_FUEL_DIESEL = 13.5       # массовая доля водорода дизеля, % (оценка) -> фактор светимости L = 336/H² (Lefebvre, по памяти)
K_FB = 10.0                # критический градиент турбулентного пограничного слоя / S_L²/α (порядок величины - Eichler,
                           # Sattelmayer 2011-2012, по памяти, не сверено)
K_CRZ = 0.5                # рециркуляция закрутки: K = K_CRZ·(крутка - 0,4)+
KV_MAX = 3.0               # предел степени рециркуляции FLOX (стеснение струи)
T_LO_CAT, T_CAT_MAX = 723.0, 1173.0   # зажигание катализатора и предел его температуры, K (оценка)
PHI_PB = 0.5               # phi бедной предгорелки каталитической схемы (оценка)
F_SEG = 0.3                # доля времени догорания, которую струйки неравномерности идут раздельно (оценка)
T_REC = 1.0e-3             # возраст рециркулирующих продуктов (догорание после зоны пламени), с (оценка)
ETA_PG = 0.05              # доля идеального прироста давления постоянного объёма, реализуемая RDC (оценка)
EPS_W = {"металл": 0.7, "керамика": 0.85}   # излучательная способность стенки: 0,7 Nimonic (по выдержке поиска), CMC - оценка
EPS_C = 0.4                # кожуха
K_WALL = {"металл": 25.0, "керамика": 8.0}  # теплопроводность стенки, Вт/(м·К): Hastelloy X около 1000 K / CMC - оценка
T_WALL_MAX = {"металл": 1173.0, "керамика": 1473.0}   # металл - targets wall_life; CMC SiC/SiC с покрытием - оценка
E_LAT, K_LAT, DH_LAT = 2.5, 10.0, 0.002     # печатная решётка: интенсификация, коэффициент потерь, гидр. диаметр (оценка)
ETA_HOLE = 0.6             # доля теплоёмкости воздуха эффузии, снятая в отверстиях (оценка)
TAU_D_EXIT = 2.0e-3        # время выравнивания неоднородности в зоне разбавления, с (оценка)
DYN_F = {"LPP": 1.0, "MM": 0.5, "FLOX": 0.5, "TVC": 0.7, "RQL": 0.6, "CAT": 0.4, "RDC": 10.0}   # эвристика пульсаций

# =========================================================================== границы: точки нагрузки и таблицы
_C = None


def _cycle():
    global _C
    if _C is None:
        cwd = os.getcwd()
        if GTU not in sys.path:
            sys.path.insert(0, GTU)
        import gtu_cycle as C   # импорт меняет текущую папку на routes/gtu - возвращаем
        os.chdir(cwd)
        _C = C
    return _C


def _mixed(phi, fuel, p, T_air, T_fuel, g=None, mech="gri30.yaml", oxid=AIR):
    g = g or ct.Solution(mech)
    g.set_equivalence_ratio(phi, fuel, oxid)
    Y = g.Y.copy()
    fsp = [s.split(":")[0] for s in fuel.split(",")]
    Z = sum(Y[g.species_index(s)] for s in fsp)
    g.TPX = T_air, p, oxid; ha = g.enthalpy_mass
    g.TPX = T_fuel, p, fuel; hf = g.enthalpy_mass
    g.HPY = Z * hf + (1 - Z) * ha, p, Y
    return g, Z


def _phi_for_Tad(fuel, p, T_air, T_fuel, T_target):
    g = ct.Solution("gri30.yaml")

    def f(phi):
        _mixed(phi, fuel, p, T_air, T_fuel, g)
        g.equilibrate("HP")
        return g.T - T_target
    return brentq(f, 0.15, 1.2, xtol=1e-5)


def _flame(args):
    xh, Tad, p, T3, Tf = args
    fuel = "CH4:1" if xh == 0 else ("H2:1" if xh == 1 else f"CH4:{1 - xh:.2f},H2:{xh:.2f}")
    phi = _phi_for_Tad(fuel, p, T3, Tf, Tad)
    g, _ = _mixed(phi, fuel, p, T3, Tf)
    alpha, Tu = g.thermal_conductivity / (g.density * g.cp_mass), g.T
    f = ct.FreeFlame(g, width=0.01)
    f.set_refine_criteria(ratio=3, slope=0.07, curve=0.14)
    f.transport_model = "mixture-averaged"
    f.solve(loglevel=0, auto=True)
    return dict(H2_vol=xh, T_ad_K=Tad, phi=phi, S_L=float(f.velocity[0]), alpha_u=alpha, T_u=Tu, points=len(f.grid))


def _ign(args):
    kind, phi, fuel, p, T_air, T_fuel, t_end = args
    if kind == "dz":
        g = ct.Solution(*DZ_MECH)
        _mixed(phi, fuel, p, T_air, T_fuel, g, oxid=AIR_REITZ)
    else:
        g, _ = _mixed(phi, fuel, p, T_air, T_fuel)
    T0_ = g.T
    r = ct.IdealGasConstPressureReactor(g, energy="on", clone=False)
    net = ct.ReactorNet([r]); net.rtol, net.atol = 1e-8, 1e-14
    t_prev, T_prev, best = 0.0, r.T, (0.0, 0.0)
    while net.time < t_end:
        net.step()
        dTdt = (r.T - T_prev) / max(net.time - t_prev, 1e-15)
        if dTdt > best[0]:
            best = (dTdt, net.time)
        t_prev, T_prev = net.time, r.T
        if r.T > T0_ + 400 and dTdt < 0.01 * best[0]:
            break
    return dict(kind=kind, fuel=fuel, phi=phi, T_K=T0_, p_kPa=p / 1e3, t_end_ms=t_end * 1e3,
                ign_ms=(best[1] * 1e3 if r.T > T0_ + 400 else None))


def _bc_inputs():
    F3 = json.load(open(os.path.join(GTU, "frame3", "frame3_point_results.json"), encoding="utf-8"))
    cyc = F3["cycle"]["frame3"]
    PR, T41_C, G = F3["split"]["pi"], 1000.0, 25.55
    for i in F3["inputs"]:
        if i["field"].startswith("slots[model=Taurus 60") and i["field"].endswith(".мягкий"):
            PR, T41_C, G = i["used"]["PR"], float(i["used"]["T41_C"]), i["used"]["G"]
    return dict(PR=PR, T41_C=T41_C, G=G, eta_pc=cyc["eta_pc"], eta_t=cyc["eta_t"], P_shaft_kW=cyc["P_shaft_kW"],
                D2_m=F3["stage2"]["D2_m"], X_LEAK=X_LEAK, DP_COMB=DP_COMB, P_EXH=P_EXH)


def _bc_key(inp):
    src = inspect.getsource(compute_boundary) + inspect.getsource(_flame) + inspect.getsource(_ign) + \
        inspect.getsource(_mixed) + json.dumps(inp, sort_keys=True) + ct.__version__
    return hashlib.sha256(src.encode()).hexdigest()[:16]


def compute_boundary(workers=4):
    """точки нагрузки KV-3 (постоянная частота и расход воздуха, p3 ~ sqrt(T41) - как targets_basis), таблицы."""
    t0 = time.time()
    C = _cycle()
    inp = _bc_inputs()
    kw = dict(eta_pc=inp["eta_pc"], eta_pt_gg=inp["eta_t"], eta_pt_pt=inp["eta_t"], x_ngv=0.0, x_rot=0.0, x_leak=X_LEAK,
              dp_duct=0.0, W1=inp["G"], dp_comb=DP_COMB, p_exh_ratio=P_EXH)
    PR, T41 = inp["PR"], inp["T41_C"]
    base = C.cycle(PR, T41, **kw)
    P100 = base["P_MW"] * 1e3
    assert abs(P100 - inp["P_shaft_kW"]) / inp["P_shaft_kW"] < 1e-3, (P100, inp["P_shaft_kW"])

    def at(T41p):
        k = math.sqrt((T41p + 273.15) / (T41 + 273.15))
        return C.cycle(PR * k, T41p, **{**kw}), PR * k
    loads = []
    air = ct.Solution("gri30.yaml")
    gst = ct.Solution("gri30.yaml"); gst.set_equivalence_ratio(1.0, "CH4:1", AIR)
    far_st = gst.Y[gst.species_index("CH4")] / (1 - gst.Y[gst.species_index("CH4")])
    for lf in LOADS:
        T41p = T41 if lf == 1.0 else brentq(lambda t: at(t)[0]["P_MW"] * 1e3 / P100 - lf, 640.0, T41, xtol=1e-4)
        r, PRp = at(T41p)
        p3, T3, W, Wf = C.P0 * PRp, r["T3_C"] + 273.15, inp["G"] * (1 - X_LEAK), r["Wf"]
        air.TPX = T3, p3, AIR
        g4, _ = _mixed((Wf / W) / far_st, "CH4:1", p3, T3, C.T0); g4.equilibrate("HP"); T4eq = g4.T
        g4.TP = r["T4_C"] + 273.15, p3; g4.equilibrate("TP")
        gs, _ = _mixed(1.0, "CH4:1", p3, T3, C.T0); gs.equilibrate("HP")
        loads.append(dict(load=lf, T41_C=T41p, PR=PRp, p3_Pa=p3, T3_K=T3, W_air=W, W_fuel=Wf, phi_ov=(Wf / W) / far_st,
                          P_shaft_kW=r["P_MW"] * 1e3, T4_K=r["T4_C"] + 273.15, T4_adiab_K=T4eq, rho4=g4.density,
                          rho3=air.density, mu3=air.viscosity, k3=air.thermal_conductivity, cp3=air.cp_mass,
                          T_st_K=gs.T))
    L0 = loads[0]
    p3, T3 = L0["p3_Pa"], L0["T3_K"]
    # равновесие (HP) метан - воздух при 100 %: phi -> T_ad, плотность (геометрия зон)
    eq = []
    ge = ct.Solution("gri30.yaml")
    for phi in list(np.round(np.arange(0.05, 3.01, 0.05), 3)):
        _mixed(float(phi), "CH4:1", p3, T3, C.T0, ge); ge.equilibrate("HP")
        eq.append([float(phi), ge.T, ge.density])
    # теплота сгорания (низшая, 298 K, H2O пар): метан и вещества недожога - GRI; н-додекан - Reitz
    def lhv(sol, sp, o2="O2", co2="CO2", h2o="H2O"):
        k = sol.species(sp)
        nC, nH = k.composition.get("C", 0), k.composition.get("H", 0)
        nO = k.composition.get("O", 0)
        sol.TP = 298.15, ct.one_atm
        h = lambda s: sol.species(s).thermo.h(298.15)          # Дж/кмоль
        dh = h(sp) + (nC + nH / 4 - nO / 2) * h(o2) - nC * h(co2) - nH / 2 * h(h2o)
        return dh / sol.molecular_weights[sol.species_index(sp)]
    LHV = {s: lhv(ge, s) for s in ("CH4", "CO", "H2", "C2H2", "C2H4", "C2H6", "CH2O", "C3H8", "CH3OH", "C2H5", "CH3")}
    dz = ct.Solution(*DZ_MECH)
    LHV_DZ = lhv(dz, "c12h26", "o2", "co2", "h2o")
    dz.set_equivalence_ratio(1.0, "c12h26:1", AIR_REITZ)
    far_st_dz = dz.Y[dz.species_index("c12h26")] / (1 - dz.Y[dz.species_index("c12h26")])
    for L in loads:
        L["W_fuel_diesel"] = L["W_fuel"] * LHV["CH4"] / LHV_DZ
        L["phi_ov_diesel"] = (L["W_fuel_diesel"] / L["W_air"]) / far_st_dz
    # таблицы: S_L(доля H2, T_ad) при 100 %; задержки самовоспламенения
    tasks_f = [(xh, Tad, p3, T3, C.T0) for xh in H2_GRID for Tad in TAD_GRID]
    tasks_i = [("gas", 0.5, ("CH4:1" if xh == 0 else ("H2:1" if xh == 1 else f"CH4:{1 - xh:.2f},H2:{xh:.2f}")),
                p3, T3, C.T0, 0.1) for xh in H2_GRID]
    tasks_i += [("dz", 0.6, "c12h26:1", loads[k]["p3_Pa"], loads[k]["T3_K"], loads[k]["T3_K"], 0.1) for k in DIESEL_LOADS]
    with Pool(max(1, min(workers, len(tasks_f)))) as pool:
        fl = pool.map(_flame, tasks_f, chunksize=1)
        ig = pool.map(_ign, tasks_i, chunksize=1)
    res = dict(
        title="Камера сгорания KV-3, уровень 0: границы (точки нагрузки) и таблицы для перебора",
        method=("точки нагрузки - gtu_cycle.cycle при постоянной частоте и расходе воздуха, PR и p3 ~ sqrt(T41) (запертое "
                "сопло), КПД узлов как в расчётной точке - допущения targets_basis.py (оценка вне расчётного режима); "
                "T41 для доли мощности 0,5 и 0,3 - brentq; равновесие HP Cantera; S_L - FreeFlame (mixture-averaged) при "
                "равной T_ad; задержка самовоспламенения - максимум dT/dt в реакторе постоянного давления; теплота "
                "сгорания - энтальпии образования механизмов при 298 K"),
        software=f"python {sys.version.split()[0]}, cantera {ct.__version__}, numpy {np.__version__}; "
                 "механизмы gri30.yaml (53/325), nDodecane_Reitz.yaml (100/432, фаза nDodecane_IG)",
        date=DATE, status="расчёт; точки частичной нагрузки - оценка (допущения targets_basis)",
        key=None, runtime_s=None,
        inputs=[dict(file="routes/gtu/frame3/frame3_point_results.json",
                     field="inputs[slots Taurus 60 генератор .мягкий].used (PR, T41_C, G); cycle.frame3 (eta_pc, eta_t, "
                           "P_shaft_kW); stage2.D2_m", used=inp),
                dict(file="routes/gtu/gtu_cycle.py", field="cycle(), AIR, T0, P0", used=dict(AIR=C.AIR, T0=C.T0, P0=C.P0)),
                dict(file="routes/gtu/combustor/targets_basis.py", field="X_LEAK, DP_COMB, P_EXH, допущения частичной нагрузки",
                     used=dict(X_LEAK=X_LEAK, DP_COMB=DP_COMB, P_EXH=P_EXH))],
        D2_m=inp["D2_m"], far_st_CH4=far_st, far_st_diesel=far_st_dz, LHV_J_kg=LHV, LHV_diesel_J_kg=LHV_DZ,
        loads=loads, eq_table=dict(cols=["phi", "T_ad_K", "rho_kg_m3"], rows=eq),
        flame_table=fl, ign_table=ig,
        check=dict(P_shaft_kW_100=P100, P_shaft_kW_frame3=inp["P_shaft_kW"],
                   T41_C_at_load=[L["T41_C"] for L in loads],
                   targets_basis_p3_kPa_100=TB["inlet"]["p3_kPa"], p3_kPa_100=p3 / 1e3,
                   targets_basis_T3_C_100=TB["inlet"]["T3_C"], T3_C_100=T3 - 273.15))
    res["key"] = _bc_key(inp)
    res["runtime_s"] = round(time.time() - t0, 1)
    return res


def _rnd(x):
    if isinstance(x, float):
        return float(f"{x:.6g}")
    if isinstance(x, dict):
        return {k: _rnd(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [_rnd(v) for v in x]
    if isinstance(x, (np.floating,)):
        return float(f"{float(x):.6g}")
    if isinstance(x, (np.integer,)):
        return int(x)
    if isinstance(x, np.bool_):
        return bool(x)
    return x


_BC = None


def boundary(workers=4, path=BC_FILE, write=True):
    """границы из файла (если ключ совпал) или расчёт; добавляет интерполяцию равновесия eq_T(phi) -> (T, rho)"""
    global _BC
    if _BC is not None:
        return _BC
    key = _bc_key(_bc_inputs())
    bc = None
    if os.path.exists(path):
        try:
            b = json.load(open(path, encoding="utf-8"))
            if b.get("key") == key:
                bc = b
        except (OSError, ValueError):
            bc = None
    if bc is None:
        bc = _rnd(compute_boundary(workers))
        if write:
            json.dump(bc, open(path, "w", encoding="utf-8"), ensure_ascii=False, indent=1)
    rows = np.array(bc["eq_table"]["rows"])
    bc["eq_T"] = lambda phi: (float(np.interp(phi, rows[:, 0], rows[:, 1])), float(np.interp(phi, rows[:, 0], rows[:, 2])))
    _BC = bc
    return bc


# =========================================================================== кинетика: потоки и реакторы
class Stream:
    """поток: расход m (кг/с), удельная энтальпия h (Дж/кг), массовые доли Y (порядок GRI), давление p (Па)"""
    __slots__ = ("m", "h", "Y", "p")

    def __init__(self, m, h, Y, p):
        self.m, self.h, self.Y, self.p = m, h, Y, p


class Kin:
    """решения Cantera одного процесса (создание Solution дорого - 0,1-0,2 с)"""

    def __init__(self):
        self.GIN, self.GR, self.GX, self.GP, self.GM = (ct.Solution("gri30.yaml") for _ in range(5))
        self.DZ = ct.Solution(*DZ_MECH)
        names = [s.upper() for s in self.GM.species_names]
        self.map_dz = [(i, names.index(s.upper())) for i, s in enumerate(self.DZ.species_names) if s.upper() in names]
        self.i = {s: self.GM.species_index(s) for s in ("NO", "NO2", "CO", "O2", "H2O", "CH4", "N2O")}
        self.HC = [k for k, s in enumerate(self.GM.species_names)
                   if self.GM.species(s).composition.get("C", 0) > 0 and s not in ("CO", "CO2")]
        self.nC = np.array([self.GM.species(s).composition.get("C", 0) for s in self.GM.species_names])
        self.cache = {}

    def air_fuel(self, T_air, T_fuel, p, fuel="CH4:1"):
        key = (round(T_air, 3), round(T_fuel, 3), round(p, 1), fuel)
        if key not in self.cache:
            g = self.GM
            g.TPX = T_air, p, AIR; ha, Ya = g.enthalpy_mass, g.Y.copy()
            g.TPX = T_fuel, p, fuel; hf, Yf = g.enthalpy_mass, g.Y.copy()
            self.cache[key] = (ha, Ya, hf, Yf)
        return self.cache[key]

    def fresh(self, m_air, m_fuel, T_air, T_fuel, p, fuel="CH4:1"):
        ha, Ya, hf, Yf = self.air_fuel(T_air, T_fuel, p, fuel)
        m = m_air + m_fuel
        return Stream(m, (m_air * ha + m_fuel * hf) / m, (m_air * Ya + m_fuel * Yf) / m, p)

    def air(self, m, T, p):
        ha, Ya, _, _ = self.air_fuel(T, T, p)
        return Stream(m, ha, Ya.copy(), p)

    @staticmethod
    def mix(*ss):
        ss = [s for s in ss if s is not None and s.m > 0]
        m = sum(s.m for s in ss)
        return Stream(m, sum(s.m * s.h for s in ss) / m, sum(s.m * s.Y for s in ss) / m, ss[0].p)

    def T(self, s):
        self.GM.HPY = s.h, s.p, s.Y
        return self.GM.T

    def rho(self, s):
        self.GM.HPY = s.h, s.p, s.Y
        return self.GM.density

    def T_eq(self, s):
        self.GM.HPY = s.h, s.p, s.Y
        self.GM.equilibrate("HP")
        return self.GM.T

    def psr(self, s, tau, init=None):
        """реактор идеального смешения: старт с равновесия (горящая ветвь) или с init=(T, p, Y); -> (выход, горит, состояние)"""
        self.GIN.HPY = s.h, s.p, s.Y
        T_in = self.GIN.T
        inlet = ct.Reservoir(self.GIN, clone=False)
        if init is None:
            self.GR.HPY = s.h, s.p, s.Y
            self.GR.equilibrate("HP")
        else:
            self.GR.TPY = init
        T_eq = self.T_eq(s)
        self.GX.TPY = self.GR.T, self.GR.P, self.GR.Y
        r = ct.IdealGasReactor(self.GR, energy="on", clone=False)
        ex = ct.Reservoir(self.GX, clone=False)
        mfc = ct.MassFlowController(inlet, r, mdot=lambda t: r.mass / tau)
        ct.PressureController(r, ex, primary=mfc, K=0.01)
        net = ct.ReactorNet([r]); net.rtol, net.atol = 1e-7, 1e-13
        net.advance_to_steady_state()
        Y = r.phase.Y.copy()
        y_in = s.Y[self.i["CH4"]]
        if y_in > 1e-5:                      # горит - израсходовано больше половины метана (годится и для горячего входа)
            lit = Y[self.i["CH4"]] < 0.5 * y_in and (r.T - T_in) > 0.25 * max(T_eq - T_in, 1.0)
        else:
            lit = (r.T - T_in) > 0.5 * max(T_eq - T_in, 1.0)
        return Stream(s.m, s.h, Y, s.p), bool(lit), (r.T, r.phase.P, Y)

    def pfr(self, s, t):
        """реактор вытеснения (лагранжев объём при постоянном давлении) на время t"""
        if t <= 0:
            return s
        self.GP.HPY = s.h, s.p, s.Y
        r = ct.IdealGasConstPressureReactor(self.GP, energy="on", clone=False)
        net = ct.ReactorNet([r]); net.rtol, net.atol = 1e-7, 1e-13
        net.advance(t)
        return Stream(s.m, s.h, r.phase.Y.copy(), s.p)

    def pfr_ign(self, s, t, n=24):
        """PFR с поиском времени расхода половины топлива (CH4) - для каталитической схемы"""
        self.GP.HPY = s.h, s.p, s.Y
        y0 = self.GP.Y[self.i["CH4"]]
        if y0 < 1e-7:                        # всё топливо превращено на катализаторе
            return self.pfr(s, t), 0.0
        r = ct.IdealGasConstPressureReactor(self.GP, energy="on", clone=False)
        net = ct.ReactorNet([r]); net.rtol, net.atol = 1e-7, 1e-13
        t_half = None
        for k in range(1, n + 1):
            net.advance(t * k / n)
            if t_half is None and r.phase.Y[self.i["CH4"]] < 0.5 * y0:
                t_half = t * k / n
        return Stream(s.m, s.h, r.phase.Y.copy(), s.p), t_half

    def psr_recirc(self, fresh, tau, K):
        """зона пламени с рециркуляцией продуктов догорания K (кг на кг свежей смеси): вход PSR = свежая + K·(выход
        PSR после T_REC догорания) - горячий газ у корня пламени (закрутка, FLOX); время пребывания по полному
        расходу tau/(1+K); первое приближение - равновесие без NOx; три итерации"""
        if K <= 1e-6:
            out, lit, st = self.psr(fresh, tau)
            return out, lit, st, fresh
        self.GM.HPY = fresh.h, fresh.p, fresh.Y
        self.GM.equilibrate("HP")
        Y = self.GM.Y.copy()
        for sp in ("NO", "NO2", "N2O"):
            Y[self.i[sp]] = 0.0
        rec = Stream(fresh.m * K, fresh.h, Y / Y.sum(), fresh.p)
        st, out, lit, inlet = None, None, False, None
        for _ in range(3):
            inlet = self.mix(fresh, rec)
            out, lit, st = self.psr(inlet, tau / (1 + K), st)
            b = self.pfr(out, T_REC)
            rec = Stream(fresh.m * K, b.h, b.Y, b.p)
        return Stream(fresh.m, out.h, out.Y, out.p), lit, st, inlet

    def tau_ext(self, inlet, tau_act, K=0.0, init=None):
        """наименьшее время пребывания до погасания (бисекция в логарифме, продолжение по горящему решению) -> запас"""
        tau_r = tau_act / (1 + K)
        out, lit, st = self.psr(inlet, tau_r, init)
        if not lit:
            for f in (3.0, 10.0, 30.0):
                o2, l2, s2 = self.psr(inlet, tau_r * f)
                if l2:
                    return 1.0 / f
            return 0.0
        lo, hi, s_hi = tau_r / 1000.0, tau_r, st
        o_lo, l_lo, _ = self.psr(inlet, lo, s_hi)
        if l_lo:
            return 1000.0
        for _ in range(7):
            mid = math.sqrt(lo * hi)
            o, l, s = self.psr(inlet, mid, s_hi)
            if l:
                hi, s_hi = mid, s
            else:
                lo = mid
        return tau_r / hi

    def from_dz(self, phi, p, T_air, T_fuel, m):
        """равновесие (HP) н-додекан - воздух в механизме Reitz -> поток GRI по общим веществам (T, p сохраняются)"""
        g = self.DZ
        _mixed(phi, "c12h26:1", p, T_air, T_fuel, g, oxid=AIR_REITZ)
        g.equilibrate("HP")
        Y = np.zeros(self.GM.n_species)
        for i_dz, i_gri in self.map_dz:
            Y[i_gri] += g.Y[i_dz]
        Y /= Y.sum()
        self.GM.TPY = g.T, p, Y
        return Stream(m, self.GM.enthalpy_mass, Y, p)

    def emissions(self, s, fuel_frac, LHV, LHV_fuel):
        """NOx, CO, UHC (ppmvd при 15 % O2; UHC как C1), O2 сухой, полнота сгорания"""
        self.GM.HPY = s.h, s.p, s.Y
        X = self.GM.X
        xw = X[self.i["H2O"]]
        o2 = 100 * X[self.i["O2"]] / (1 - xw)
        k = (O2_AIR - O2_REF) / max(O2_AIR - o2, 0.1)
        nox = 1e6 * (X[self.i["NO"]] + X[self.i["NO2"]]) / (1 - xw) * k
        co = 1e6 * X[self.i["CO"]] / (1 - xw) * k
        uhc = 1e6 * float(np.dot(X[self.HC], self.nC[self.HC])) / (1 - xw) * k
        loss = sum(self.GM.Y[self.GM.species_index(sp)] * v for sp, v in LHV.items())
        eta = 1.0 - loss / (fuel_frac * LHV_fuel)
        return dict(NOx=nox, CO=co, UHC=uhc, O2_dry=o2, eta_c=eta, T=self.GM.T)


_K = None


def kin():
    global _K
    if _K is None:
        _K = Kin()
    return _K


GH = {n: (lambda z, w: (z, w / w.sum()))(*np.polynomial.hermite_e.hermegauss(n)) for n in (1, 3, 5)}
N_NODES = {"LPP": 5, "MM": 5, "CAT": 5, "FLOX": 3, "TVC": 3, "RDC": 3, "RQL": 1}


# =========================================================================== пассивное деление топлива
def pilot_share(x100, m_rel):
    """доля топлива пилота при относительном расходе m_rel: пилот - отверстие (m ~ ΔP^0,5), основной - ламинарный
    канал (m ~ ΔP); общий перепад; при 100 % доля x100 (оценка: газ как несжимаемый)"""
    if m_rel >= 0.999999 or x100 <= 0:
        return x100
    cA, cB = x100, 1 - x100
    s = (-cA + math.sqrt(cA * cA + 4 * cB * m_rel)) / (2 * cB)
    return cA * s / m_rel


# =========================================================================== сети по схемам (газ)
def _node_streams(K, m_air, m_fuel, S, n_nodes, T3, Tf, p, fuel="CH4:1"):
    """свежая смесь, разбитая на узлы нормальной плотности phi (Гаусс - Эрмит): доля воздуха - вес, топливо ~ phi"""
    z, w = GH[n_nodes]
    out = []
    for zk, wk in zip(z, w):
        f = max(0.02, 1 + S * zk)
        out.append((wk, K.fresh(m_air * wk, m_fuel * wk * f, T3, Tf, p, fuel)))
    tot = sum(s.m for _, s in out)
    return out, tot


def gas_network(c, d, bc, li):
    """сеть реакторов схемы на газе при нагрузке LOADS[li] -> выбросы, температуры, запас по срыву.
    Узлы неравномерности (и пилот) раздельно проходят PSR и первую долю F_SEG догорания, затем смешиваются и догорают
    вместе (бедные струйки дожигаются в горячем газе), затем разбавление и PFR зоны разбавления."""
    K = kin()
    L = bc["loads"][li]
    p, T3, W, Wf = L["p3_Pa"], L["T3_K"], L["W_air"], L["W_fuel"]
    Tf = 288.15
    s = c["scheme"]
    m_rel = Wf / bc["loads"][0]["W_fuel"]
    tp, tb, td = 1e-3 * c["tau_pz_ms"], 1e-3 * c["tau_bo_ms"], 1e-3 * c["tau_dil_ms"]
    S = d["S_unmix"]
    a_pz, a_aux = c["a_pz"], c["a_aux"]
    res = dict(load=L["load"])
    seg, Ts, bypass_air, margins = [], [], 0.0, []
    lit_frac = 1.0
    K_crz = K_CRZ * max(0.0, c["swirl"] - 0.4) if c["injector"] == "закрутка" else 0.0
    nn = N_NODES[s]

    def nodes_chain(m_air, m_fuel, S_eff, Krec, pre=None):
        """узлы основной зоны: [горячий газ пилота/полости] + свежая -> PSR (рециркуляция продуктов) -> PFR F_SEG·tb"""
        nodes, _ = _node_streams(K, m_air, m_fuel, S_eff, nn, T3, Tf, p)
        lit_w = 0.0
        for wk, fr in nodes:
            inlet = fr if pre is None else K.mix(fr, Stream(pre.m * wk, pre.h, pre.Y, pre.p))
            o, lit, st, _ = K.psr_recirc(inlet, tp, Krec)
            o = K.pfr(o, F_SEG * tb)
            seg.append(o); Ts.append((wk, K.T(o)))
            lit_w += wk if lit else 0.0
        return lit_w

    if s in ("LPP", "MM", "FLOX", "TVC"):
        Krec, S_eff = K_crz, S
        if s == "FLOX":
            T_b = max(bc["eq_T"](L["phi_ov"] / a_pz)[0], 800.0)
            Krec = min(KV_MAX, max(0.0, 0.32 * (0.5 * d["L_L"] / (1e-3 * d["d_air_mm"])) * math.sqrt(T3 / T_b) - 1.0))
            S_eff = S / (1 + Krec)
        res["K_recirc"] = Krec
        pilot = None
        x100 = CS.fuel_split100(c, bc["loads"][0]["phi_ov"])
        if s == "TVC":
            x = pilot_share(x100, m_rel) if c["passive"] == "гидравлика" else x100
            pilot = dict(m_air=a_aux * W, m_fuel=x * Wf)
            m_air_main, m_fuel_main = a_pz * W, (1 - x) * Wf
        elif c["passive"] == "гидравлика":
            x = pilot_share(x100, m_rel)
            pilot = dict(m_air=a_aux * a_pz * W, m_fuel=x * Wf)
            m_air_main, m_fuel_main = (1 - a_aux) * a_pz * W, (1 - x) * Wf
        else:
            m_air_main, m_fuel_main = a_pz * W, Wf
        if c["passive"] == "струи":
            a_eff = min(1.0, m_rel ** c["k_pass"])
            bypass_air = (1 - a_eff) * m_air_main
            m_air_main *= a_eff
            S_eff = S_eff + 0.3 * (1 - a_eff)
            res["a_eff_frac"] = a_eff
        res["S_eff"] = S_eff
        pre = None
        if pilot is not None:
            Kp = 0.0 if s == "TVC" else K_crz
            pf = K.fresh(pilot["m_air"], pilot["m_fuel"], T3, Tf, p)
            po, plit, pst, pin = K.psr_recirc(pf, tp, Kp)
            margins.append(("полость" if s == "TVC" else "пилот", K.tau_ext(pin, tp, Kp, pst)))
            res["phi_pilot"] = (pilot["m_fuel"] / pilot["m_air"]) / bc["far_st_CH4"]
            res["fuel_pilot_share"] = pilot["m_fuel"] / Wf
            if plit:
                pre = po                       # горячий газ пилота/полости уходит в основную зону
            else:
                seg.append(K.pfr(po, F_SEG * tb))
        lit_frac = nodes_chain(m_air_main, m_fuel_main, S_eff, Krec, pre)
        fr = K.fresh(m_air_main, m_fuel_main, T3, Tf, p)
        if pre is not None:                   # основная зона с горячим газом пилота/полости (средний узел)
            fr = K.mix(fr, pre)
        _, _, st0, inl0 = K.psr_recirc(fr, tp, Krec)
        margins.append(("основная", K.tau_ext(inl0, tp, Krec, st0)))
        res["phi_main"] = (m_fuel_main / m_air_main) / bc["far_st_CH4"]
    elif s == "RQL":
        fr = K.fresh(a_pz * W, Wf, T3, Tf, p)
        ro, rlit, rst, rin = K.psr_recirc(fr, tp, K_crz)
        margins.append(("богатая", K.tau_ext(rin, tp, K_crz, rst)))
        res["phi_rich"] = (Wf / (a_pz * W)) / bc["far_st_CH4"]
        tq = 1e-3 * next(z["tau_ms"] for z in d["zones"] if z["name"] == "гашение")
        cur = ro
        for k in range(3):                      # гашение: воздух тремя долями, каждая - PSR tq/3
            cur = K.mix(cur, K.air(a_aux * W / 3, T3, p))
            cur, _, _ = K.psr(cur, tq / 3)
        seg.append(cur); Ts.append((1.0, K.T(cur)))
        lit_frac = 1.0 if rlit else 0.0
    elif s == "CAT":
        m_air = a_pz * W
        pre, f_pb, a_pb = None, 0.0, 0.0
        if T3 < T_LO_CAT:                     # бедная предгорелка phi 0,5 поднимает вход катализатора выше T_LO_CAT
            f_pb = min(0.5, L["cp3"] * (T_LO_CAT + 40.0 - T3) * m_air / (bc["LHV_J_kg"]["CH4"] * Wf))
            a_pb = f_pb * Wf / (PHI_PB * bc["far_st_CH4"])
            pre, plit, pst = K.psr(K.fresh(a_pb, f_pb * Wf, T3, Tf, p), 1e-3)
        res["preburner_fuel_share"] = f_pb
        nodes, _ = _node_streams(K, m_air - a_pb, (1 - f_pb) * Wf, S, nn, T3, Tf, p)
        conv, t_half = [], []
        G = K.GM
        iC, iO2, iCO2, iH2O = (G.species_index(z) for z in ("CH4", "O2", "CO2", "H2O"))
        MW = G.molecular_weights
        for wk, fr in nodes:
            inlet = fr if pre is None else K.mix(fr, Stream(pre.m * wk, pre.h, pre.Y, pre.p))
            T_in, T_ad = K.T(inlet), K.T_eq(inlet)
            x = min(1.0, max(0.0, (T_CAT_MAX - T_in) / max(T_ad - T_in, 1.0))) if T_in >= T_LO_CAT else 0.0
            conv.append(wk * x)
            Y = inlet.Y.copy()
            dn = x * Y[iC] / MW[iC]
            Y[iC] -= dn * MW[iC]; Y[iO2] -= 2 * dn * MW[iO2]; Y[iCO2] += dn * MW[iCO2]; Y[iH2O] += 2 * dn * MW[iH2O]
            o, th = K.pfr_ign(Stream(inlet.m, inlet.h, Y, p), F_SEG * tb + (1 - F_SEG) * tb, n=24)
            t_half.append(th)
            seg.append(o); Ts.append((wk, K.T(o)))
        res["catalyst_conversion"] = sum(conv)
        th_mid = t_half[len(t_half) // 2]
        margins.append(("гомогенное догорание", 0.0 if th_mid is None else (1000.0 if th_mid == 0 else tb / th_mid)))
        lit_frac = sum(wk for (wk, _), th in zip(nodes, t_half) if th is not None)
    elif s == "RDC":
        nodes, _ = _node_streams(K, a_pz * W, Wf, S, nn, T3, Tf, p)
        pcv = []
        for wk, fr in nodes:
            # сгорание при постоянном объёме вместо волны (цикл Хамфри, оценка): равновесие UV, NO/NO2/N2O обнулены;
            # NOx - закрытый реактор на время tau_pz; изоэнтропическое расширение (замороженный состав) до p3
            K.GP.HPY = fr.h, fr.p, fr.Y
            K.GP.equilibrate("UV")
            Y = K.GP.Y.copy()
            for sp in ("NO", "NO2", "N2O"):
                Y[K.i[sp]] = 0.0
            K.GP.TDY = K.GP.T, K.GP.density, Y / Y.sum()
            r = ct.IdealGasReactor(K.GP, energy="on", clone=False)
            net = ct.ReactorNet([r]); net.rtol, net.atol = 1e-7, 1e-13
            net.advance(tp)
            pcv.append(r.phase.P / p)
            K.GM.TPY = r.T, r.phase.P, r.phase.Y
            K.GM.SP = K.GM.entropy_mass, p
            o = K.pfr(Stream(fr.m, K.GM.enthalpy_mass, K.GM.Y.copy(), p), F_SEG * tb)
            seg.append(o); Ts.append((wk, K.T(o)))
        res["p_cv_ratio"] = float(np.mean(pcv))
        margins.append(("детонация", 1.0))
    # смешение струек, совместное догорание, разбавление, выход
    if s == "CAT":
        hot = K.mix(*seg)
    else:
        hot = K.pfr(K.mix(*seg), (1 - F_SEG) * tb)
    T_hot = K.T(hot)
    wsum = max(sum(w for w, _ in Ts), 1e-12)
    T_mean = sum(w * T for w, T in Ts) / wsum
    sigT = math.sqrt(sum(w * (T - T_mean) ** 2 for w, T in Ts) / wsum)
    dil = K.air(d["a_dil"] * W + d["a_cool"] * W + bypass_air, T3, p)
    ex = K.pfr(K.mix(hot, dil), td)
    em = K.emissions(ex, Wf / ex.m, bc["LHV_J_kg"], bc["LHV_J_kg"]["CH4"])
    res.update(NOx=em["NOx"], CO=em["CO"], UHC=em["UHC"], eta_c=em["eta_c"], O2_dry=em["O2_dry"], T_exit=em["T"],
               T_hot=T_hot, T_flame=T_mean, sigma_T_hot=sigT, m_hot=hot.m, lit_frac=lit_frac,
               margin_lbo=min(m for _, m in margins) if margins else None, margins={k: v for k, v in margins})
    if s == "RDC":
        res["p_gain"] = ETA_PG * max(0.0, res["p_cv_ratio"] - 1.0)
    return res


# =========================================================================== дизель
def _evap(K, SMD, T_g, p):
    """время испарения по закону d²: t = SMD²/λ, λ = 8·(k/cp)_г·ln(1+B)/ρ_ж, B = cp(T_г - T_кип)/L (оценка)"""
    Tf = 0.5 * (T_g + DIESEL_LIQ["T_b"])
    K.GM.TPX = Tf, p, AIR
    kg, cp = K.GM.thermal_conductivity, K.GM.cp_mass
    B = max(1e-3, cp * (T_g - DIESEL_LIQ["T_b"]) / DIESEL_LIQ["L_v"])
    lam = 8 * kg * math.log(1 + B) / (DIESEL_LIQ["rho"] * cp)
    return SMD ** 2 / lam, B


def smd(c, d, bc, li, Wfd):
    """средний заутеров диаметр: центробежная - 2,25·σ^0,25·μ^0,25·ṁ^0,25·ΔP^-0,5·ρA^-0,25; пневматическая и
    испарительная - предплёночная El-Shanawany & Lefebvre (обе формы - по памяти, не сверены с первоисточником)"""
    L = bc["loads"][li]
    lq = DIESEL_LIQ
    rhoA = L["rho3"]
    if c["diesel"] == "центробежная":
        n = max(1, c["n_el"])
        dPf = DP_FUEL_PS * (Wfd / bc["loads"][0]["W_fuel_diesel"]) ** 2
        return 2.25 * lq["sigma"] ** 0.25 * lq["mu"] ** 0.25 * (Wfd / n) ** 0.25 * dPf ** -0.5 * rhoA ** -0.25, dPf
    Dp = 1e-3 * d["d_air_mm"]
    U = math.sqrt(2 * c["dp_inj"] * L["p3_Pa"] / rhoA)
    ALR = c["a_pz"] * L["W_air"] / Wfd
    t1 = 0.33 * (lq["sigma"] / (rhoA * U ** 2 * Dp)) ** 0.6 * (lq["rho"] / rhoA) ** 0.1
    t2 = 0.068 * (lq["mu"] ** 2 / (lq["sigma"] * lq["rho"] * Dp)) ** 0.5
    return Dp * (1 + 1 / ALR) * (t1 + t2), None


def diesel_network(c, d, bc, li, gas_res):
    K = kin()
    L = bc["loads"][li]
    p, T3, W = L["p3_Pa"], L["T3_K"], L["W_air"]
    Wfd = L["W_fuel_diesel"]
    far_d = bc["far_st_diesel"]
    tp, tb, td = 1e-3 * c["tau_pz_ms"], 1e-3 * c["tau_bo_ms"], 1e-3 * c["tau_dil_ms"]
    D, dPf = smd(c, d, bc, li, Wfd)
    T_fl = max(gas_res["T_flame"], 1500.0)
    t_ev_air, B_air = _evap(K, D, T3, p)
    t_ev_fl, B_fl = _evap(K, D, T_fl, p)
    u = math.sqrt(2 * c["dp_inj"] * p / L["rho3"])
    l_mix = 1e-3 * d["d_air_mm"] if c["diesel"] != "центробежная" else d["D_h"] / math.sqrt(max(1, c["n_el"]) /
                                                                                            CS.N_CAN[c["chamber"]])
    t_mix = K_MIX_D * l_mix / u
    fl = [r for r in bc["flame_table"] if r["H2_vol"] == 0.0]
    st = min(fl, key=lambda r: abs(r["T_ad_K"] - 2200.0))
    t_chem = T_CHEM_K * st["alpha_u"] / st["S_L"] ** 2
    tau_mix_pre = 1e-3 * d["tau_mix_ms"]
    res = dict(load=L["load"], SMD_um=1e6 * D, dP_fuel_MPa=(dPf / 1e6 if dPf else None), t_ev_air_ms=1e3 * t_ev_air,
               t_ev_flame_ms=1e3 * t_ev_fl, t_mix_ms=1e3 * t_mix, t_chem_ms=1e3 * t_chem)
    x_ev = 0.0
    if c["diesel"] == "испарительная":
        x_ev = 1 - math.exp(-tau_mix_pre / t_ev_air)
        ig = next(r for r in bc["ign_table"] if r["kind"] == "dz" and abs(r["p_kPa"] - p / 1e3) < 1.0)
        res["ign_ms"] = ig["ign_ms"]
        res["ign_margin"] = (ig["ign_ms"] or ig["t_end_ms"]) / max(1e3 * tau_mix_pre, 1e-6)
    x_diff = (1 - x_ev) * t_mix / (t_mix + t_chem)
    f_unev = math.exp(-(tp + tb) / t_ev_fl)
    res.update(x_prevap=x_ev, x_diff=x_diff, f_unevaporated=f_unev)
    if c["scheme"] == "RQL":
        phi_r = (Wfd / (c["a_pz"] * W)) / far_d
        cur = K.from_dz(phi_r, p, T3, T3, c["a_pz"] * W + Wfd)
        cur = K.pfr(cur, tp)
        tq = 1e-3 * next(z["tau_ms"] for z in d["zones"] if z["name"] == "гашение")
        for k in range(3):
            cur = K.mix(cur, K.air(c["a_aux"] * W / 3, T3, p))
            cur, _, _ = K.psr(cur, tq / 3)
        hot = K.pfr(cur, tb)
        res["phi_rich"] = phi_r
    else:
        a_zone = (c["a_pz"] + (c["a_aux"] if c["scheme"] == "TVC" else 0.0)) * W
        air_d = x_diff * Wfd / far_d
        m_lean_air = max(a_zone - air_d, 0.05 * a_zone)
        phi_lean = (1 - x_diff) * Wfd / m_lean_air / far_d
        parts = []
        if phi_lean > 0.05:
            ln = K.from_dz(phi_lean, p, T3, T3, m_lean_air + (1 - x_diff) * Wfd)
            parts.append(K.pfr(ln, tp + tb))
        if x_diff > 1e-4:
            df = K.from_dz(1.0, p, T3, T3, air_d + x_diff * Wfd)
            parts.append(K.pfr(df, min(tp + tb, t_mix + t_chem)))   # пар горит у стехиометрии, пока не перемешан
        hot = K.mix(*parts)
        res["phi_lean"] = phi_lean
    T_hot = K.T(hot)
    used_air = hot.m - Wfd
    dil = K.air(max(W - used_air, 0.0), T3, p)
    ex = K.pfr(K.mix(hot, dil), td)
    em = K.emissions(ex, Wfd / ex.m, bc["LHV_J_kg"], bc["LHV_diesel_J_kg"])
    eta = em["eta_c"] - f_unev
    res.update(NOx=em["NOx"], CO=em["CO"], UHC=em["UHC"], eta_c=eta, T_hot=T_hot, T_exit=em["T"],
               margin_lbo=(gas_res["margin_lbo"] or 0.0) * tp / (tp + t_ev_fl))
    return res


# =========================================================================== проскок и самовоспламенение
def flashback(c, d, bc):
    """запас по проскоку в пограничном слое смесителя: g_стенки(u_min) / (K_FB·S_L²/α); u_min = u·w: закрутка
    w = 1/(1 + 4·крутка²) (вихревое ядро), струи и микросмесители w = 0,6 (следы струй) - оценка; S_L - таблица по T_ad
    богатого края неравномерности phi·(1 + 2S); H2_vol_max - наибольшая доля H2 таблицы с запасом >= 1"""
    L0 = bc["loads"][0]
    if c["scheme"] == "RQL":
        return dict(margin_CH4=None, H2_vol_max=1.0, note="богатая зона без предварительного смешения")
    phi_m = d["phi_pz_100"]
    phi_r = phi_m * (1 + 2 * d["S_unmix"])
    Tad = bc["eq_T"](phi_r)[0]
    w = 1 / (1 + 4 * c["swirl"] ** 2) if c["injector"] == "закрутка" else 0.6
    u = d["u_a"] * w
    nu = L0["mu3"] / L0["rho3"]
    Re = u * 1e-3 * d["d_air_mm"] / nu
    f = 0.316 * max(Re, 2300.0) ** -0.25
    g_w = f * u ** 2 / (8 * nu)
    out = dict(T_ad_rich_K=Tad, u_min=u, g_wall=g_w)
    h2max = None
    for xh in H2_GRID:
        rows = sorted([r for r in bc["flame_table"] if r["H2_vol"] == xh], key=lambda r: r["T_ad_K"])
        SL = float(np.interp(Tad, [r["T_ad_K"] for r in rows], [r["S_L"] for r in rows]))
        al = float(np.interp(Tad, [r["T_ad_K"] for r in rows], [r["alpha_u"] for r in rows]))
        m = g_w / (K_FB * SL ** 2 / al)
        ig = next(r for r in bc["ign_table"] if r["kind"] == "gas" and
                  (r["fuel"] == ("CH4:1" if xh == 0 else ("H2:1" if xh == 1 else f"CH4:{1 - xh:.2f},H2:{xh:.2f}"))))
        ig_m = (ig["ign_ms"] or ig["t_end_ms"]) / max(d["tau_mix_ms"], 1e-6)
        out[f"margin_H2_{xh}"] = m
        out[f"ign_margin_H2_{xh}"] = ig_m
        if m >= 1 and ig_m >= 3:
            h2max = xh
        else:
            break
    out["margin_CH4"] = out["margin_H2_0.0"]
    out["H2_vol_max"] = h2max if h2max is not None else 0.0
    return out


# =========================================================================== стенка, неравномерность, потери
def wall_temperature(c, d, bc, T_g, far_zone, lum=1.0):
    """одномерный баланс стенки: R1 + C1 = k/t·(Tw1 - Tw2) = R2 + C2 [+ теплоёмкость эффузии]. Форма - Lefebvre &
    Ballal, Gas Turbine Combustion, 3-е изд. (2010), гл. 9 (по памяти, не сверено; константа 290 и l_b = 3,4·V/A -
    по выдержке поиска). Плёнка - эффективность щели η = 1,10·m^0,65·(μa/μg)^0,15·(x/s)^-0,2·(t/s)^-0,2 (по памяти)."""
    K = kin()
    L0 = bc["loads"][0]
    p, T3, W = L0["p3_Pa"], L0["T3_K"], L0["W_air"]
    N = CS.N_CAN[c["chamber"]]
    mat = "керамика" if c["cooling"] == "керамика" else "металл"
    ew, kw, tw = EPS_W[mat], K_WALL[mat], CS.WALL_T[mat]
    Dh = d["D_h"]
    lb = 0.85 * Dh
    P_kPa = p / 1e3
    eg = 1 - math.exp(-290 * P_kPa * lum * (far_zone * lb) ** 0.5 * T_g ** -1.5)
    K.GM.TPX = T_g, p, "CO2:0.06,H2O:0.12,O2:0.09,N2:0.73"
    kg, mug, rhog = K.GM.thermal_conductivity, K.GM.viscosity, K.GM.density
    mg = c["a_pz"] * W / N
    A_L = d["A_L"]
    h1 = 0.020 * kg / Dh ** 0.2 * (mg / (A_L * mug)) ** 0.8 / Dh ** 0.0
    ka, mua = L0["k3"], L0["mu3"]
    m_an = W * max(0.2, 1 - c["a_pz"]) / N
    A_an = m_an / (L0["rho3"] * CS.U_ANN)
    D_an = 2 * A_an / (math.pi * (Dh if c["chamber"] in ("выносная", "трубчатые") else 2 * d["D_m"]))
    h2 = 0.020 * ka / D_an ** 0.2 * (m_an / (A_an * mua)) ** 0.8
    T_aw, q_hole_k = T_g, 0.0
    Ug = mg / (rhog * A_L)
    Uc = 0.5 * math.sqrt(2 * c["dp_inj"] * p / L0["rho3"])
    blow = L0["rho3"] * Uc / (rhog * Ug)
    if c["cooling"] in ("плёночное", "отверстия", "керамика", "решётка") and c["porosity"] > 0:
        xs = 0.5 / max(c["porosity"], 1e-4)
        eta_f = 1.10 * max(blow, 0.05) ** 0.65 * (mua / mug) ** 0.15 * xs ** -0.2 * 0.5 ** -0.2
        if c["cooling"] == "решётка":
            eta_f *= 0.5                                 # выпотевание малой пористости - половина щелевой (оценка)
        eta_f = min(0.95, max(0.0, eta_f))
        T_aw = T_g - eta_f * (T_g - T3)
    else:
        eta_f = 0.0
    if c["cooling"] == "отверстия":
        q_hole_k = ETA_HOLE * (d["a_cool"] * W / d["A_wall"]) * L0["cp3"]
    if c["cooling"] == "решётка":
        G = math.sqrt(2 * L0["rho3"] * c["dp_wall"] * p / K_LAT)
        h2 = E_LAT * 0.023 * ka / DH_LAT * (G * DH_LAT / mua) ** 0.8 * 0.7 ** 0.4
    Z = ew * EPS_C / (EPS_C + ew * (1 - EPS_C) * min(1.0, Dh / max(d["D_case"], Dh)))

    def outer(Tw2):
        return Z * SIGMA * (Tw2 ** 4 - T3 ** 4) + (h2 + q_hole_k) * (Tw2 - T3)

    def resid(Tw1):
        Tw2 = brentq(lambda t: kw / tw * (Tw1 - t) - outer(t), T3, Tw1)
        R1 = 0.5 * SIGMA * (1 + ew) * eg * T_g ** 1.5 * (T_g ** 2.5 - Tw1 ** 2.5)
        C1 = h1 * (T_aw - Tw1)
        return R1 + C1 - kw / tw * (Tw1 - Tw2), Tw2
    lo, hi = T3 + 1e-3, T_g - 1e-3
    if resid(hi)[0] > 0:
        Tw1 = hi
    else:
        Tw1 = brentq(lambda t: resid(t)[0], lo, hi, xtol=0.05)
    return dict(T_w1=Tw1, T_w2=resid(Tw1)[1], eps_g=eg, eta_film=eta_f, T_aw=T_aw, h1=h1, h2=h2, mat=mat,
                T_lim=T_WALL_MAX[mat], margin_K=T_WALL_MAX[mat] - Tw1)


def pattern_factor(c, d, bc, g100):
    """PF = sqrt(PF_разб² + PF_смеш²): PF_разб = θ·[1 - exp(-1/(A·L/D·ΔPL/qref))], A = 0,070 трубчатые и 0,050
    кольцевые (форма Lefebvre, по памяти, не сверено), θ = (T_гор - T4)/(T_стех - T4) - перенос с обычных камер на
    бедные; PF_смеш = 2·σT·(m_гор/m)·exp(-τ_разб/τ_выр)/(T4 - T3). RTDF = 0,5·PF (оценка)"""
    L0 = bc["loads"][0]
    T3, T4 = L0["T3_K"], g100["T_exit"]
    A = 0.070 if c["chamber"] in ("выносная", "трубчатые") else 0.050
    x = A * c["LD"] * d["dpL_qref"]
    pf_l = 1 - math.exp(-1 / max(x, 1e-6))
    theta = min(1.0, max(0.0, (g100["T_hot"] - T4) / (L0["T_st_K"] - T4)))
    m_tot = L0["W_air"] + L0["W_fuel"]
    pf_m = 2 * g100["sigma_T_hot"] * (g100["m_hot"] / m_tot) * math.exp(-1e-3 * c["tau_dil_ms"] / TAU_D_EXIT) / (T4 - T3)
    pf = math.sqrt((theta * pf_l) ** 2 + pf_m ** 2)
    return dict(PF=pf, RTDF=0.5 * pf, PF_dil=theta * pf_l, PF_mix=pf_m, lefebvre_x=x, theta=theta)


def pressure_loss(c, d, bc, g100):
    """полные потери: диффузор + смеситель (стенка параллельно) + решётка последовательно + горячие + переходный канал
    (- прирост давления RDC). Горячие: (γ/2)·M²·(T4/T3 - 1) по скорости входа в жаровую трубу (Рэлей, малое M)"""
    L0 = bc["loads"][0]
    N = CS.N_CAN[c["chamber"]]
    U_in = L0["W_air"] / N / (L0["rho3"] * d["A_L"])
    M = U_in / math.sqrt(1.4 * 287.0 * L0["T3_K"])
    dp_hot = 0.7 * M ** 2 * (g100["T_exit"] / L0["T3_K"] - 1)
    gain = g100.get("p_gain", 0.0)
    tot = d["dp_diff"] + c["dp_inj"] + d["dp_series_wall"] + dp_hot + d["dp_duct"] - gain
    sl = TB["dp_sensitivity"]
    return dict(dp_total=tot, dp_diff=d["dp_diff"], dp_inj=c["dp_inj"], dp_series_wall=d["dp_series_wall"], dp_hot=dp_hot,
                dp_duct=d["dp_duct"], p_gain=gain,
                dP_shaft_kW=sl["dP_kW_per_pct"] * 100 * (tot - DP_COMB),
                d_eff_pt=sl["d_eff_pt_per_pct"] * 100 * (tot - DP_COMB))


# =========================================================================== модели модуля
def _printer_segments(dims):
    """наименьшее число частей, при котором деталь (диаметр, длина) помещается в камеру одного из принтеров
    make_model.PRINTERS (воронка двигателя; камеры - по выдержке поиска)"""
    MK = _make_model()
    D, Lx = dims
    best = None
    for pr in MK.PRINTERS:
        if pr["proc"] != "LPBF":
            continue
        if "box" in pr:
            xy, z = min(pr["box"][:2]) / 1e3, pr["box"][2] / 1e3
        else:
            xy, z = pr["round_d"] / 1e3 * 0.7071, pr["z"] / 1e3
        k = 1
        if D <= xy and Lx <= z:
            k = 1
        elif D <= z and Lx <= xy:
            k = 1
        else:
            k_d = 1 if D <= xy else math.ceil(math.pi * D / (0.9 * xy))
            k_l = 1 if Lx <= z else math.ceil(Lx / z)
            k = k_d * k_l
        if best is None or k < best[0]:
            best = (k, pr["name"])
    return best


_MK = None


def _make_model():
    global _MK
    if _MK is None:
        sys.path.insert(0, os.path.join(GTU, "funnel"))
        import make_model as MK
        _MK = MK
    return _MK


def module_models(c, d, bc, fb, pl, wall_gas):
    """детали, замена картриджа, Plug&Play, изготовимость - правила подсчёта (оценка)"""
    ch, mk, inj = c["chamber"], c["make"], c["injector"]
    N = CS.N_CAN[ch]
    D_out = d["D_case"] if ch in ("выносная", "трубчатые") else d["D_env"]
    seg = None
    parts = dict()
    if mk == "печатный картридж":
        seg = _printer_segments((D_out, d["L_mod"]))
        parts["корпус-картридж (труба, смесители, коллекторы, штуцеры)"] = seg[0] * N
    elif mk == "гибрид":
        parts["литой корпус-труба"] = N
        parts["печатная головка смесителей"] = N
        seg = _printer_segments((D_out * 0.6, 0.3))
    else:
        parts["литая жаровая труба"] = N
        parts["литой купол"] = N
        parts["топливный коллектор"] = 1
        parts["форсунки (наконечники)"] = min(c["n_el"], 24) if inj == "закрутка" else max(1, min(c["n_el"], 400) // 20)
        parts["топливные штуцеры"] = 2
    parts["воспламенитель"] = 1 if ch != "трубчатые" else 2
    parts["уплотнение"] = N
    parts["датчик (пламя, пульсации, стенка)"] = 1
    if c["diesel"] == "центробежная":
        parts["центробежные форсунки дизеля"] = c["n_el"]
    elif mk == "литьё":
        parts["форсунки дизеля"] = max(1, min(c["n_el"], 400) // 20)
    if c["cooling"] == "керамика":
        parts["вкладыш CMC и крепления"] = 2 * N if ch in ("выносная", "трубчатые") else 2 + (seg[0] if seg else 4)
    if c["scheme"] == "CAT":
        parts["каталитический блок"] = N
        if bc["loads"][-1]["T3_K"] < T_LO_CAT:
            parts["подогреватель (предгорелка)"] = N
    if c["scheme"] == "RDC":
        parts["инициатор детонации"] = 1
    if c["scheme"] == "TVC" and mk == "литьё":
        parts["полость (кольцо)"] = N
    if ch == "трубчатые":
        parts["переходные патрубки"] = N
        parts["трубки переброса пламени"] = N
    n_parts = int(sum(parts.values()))
    welds = 0 if mk != "литьё" else 1 + parts.get("форсунки (наконечники)", 0)
    align = {"выносная": 0, "трубчатые": 2 * N, "кольцевая": 2, "противоточная": 2}[ch]
    # масса картриджа: стенка трубы + купол со смесителями (оценка)
    rho_m = 8220.0 if c["cooling"] != "керамика" else 2700.0
    t_eff = CS.WALL_T["металл"] * (2.2 if c["cooling"] == "решётка" else 1.0)
    A_dome = d["A_L"] * N
    mass = d["A_wall"] * t_eff * rho_m + A_dome * 0.010 * 8220.0 * 0.5
    if ch == "выносная":
        mass += 0.15 * mass                      # фланец картриджа (оценка)
    # замена картриджа (оценка)
    swap = {"выносная": 0.8, "трубчатые": 0.6 * N + 1.0, "кольцевая": 12.0, "противоточная": 6.0}[ch]
    people = 1 if mass <= 25 else (2 if mass <= 50 else 2)
    if mass > 50:
        swap += 0.5                              # таль
    tools = int(c["diesel"] == "центробежная") + int(c["cooling"] == "керамика") + int(align > 0)
    keyed = mk != "литьё"
    s_time = min(1.0, max(0.0, (8.0 - swap) / 7.0))
    s_mass = 1.0 if mass <= 50 else max(0.0, 1 - (mass - 50) / 250)
    srv = 10 * (0.4 * s_time + 0.2 * s_mass + 0.2 * (1.0 if keyed and align == 0 else 0.3) + 0.2 / (1 + tools))
    # Plug&Play (оценка)
    valves_gas = 1 + int("подогреватель (предгорелка)" in parts) + int(c["scheme"] == "RDC")
    valves_dz = 1
    pr = pl["p_fuel_ratio"]
    s_pf = 1.0 if pr <= 1.3 else max(0.0, 1 - (pr - 1.3) / 0.7)
    chg = {("закрутка", "центробежная"): 1.0, ("закрутка", "пневматическая"): 1.0, ("струи", "пневматическая"): 1.0,
           ("микросмесители", "пневматическая"): 0.8, ("струи", "центробежная"): 0.7}.get((inj, c["diesel"]), 0.6)
    if c["scheme"] in ("CAT", "RDC"):
        chg = min(chg, 0.3)
    iface = {"кольцевая": 1.0, "противоточная": 1.0, "трубчатые": 0.85, "выносная": 0.75}[ch]
    s_valves = 1.0 if (valves_gas == 1 and valves_dz == 1) else 0.3
    pp = 10 * (0.35 * s_valves + 0.25 * s_pf + 0.25 * chg + 0.15 * iface)
    # изготовимость (форма make_model воронки двигателя)
    MK = _make_model()
    trl = min(CS.TRL_SCHEME[c["scheme"]], CS.TRL_DIESEL[c["diesel"]], CS.TRL_COOL[c["cooling"]], CS.TRL_MAKE[mk])
    if mk in ("печатный картридж", "гибрид"):
        s_fit = 1.0 / (seg[0] if seg else 1)
        s_proc = MK.MAT["HastX"]["print"] * (0.8 if d["d_jet_mm"] < 0.5 else 1.0)
        vol_metal = mass / 8220.0
        rate = 900.0                                   # см3/ч: BLT-S1500 (make_model.PRINTERS, по выдержке поиска)
        weeks = 1e6 * vol_metal / rate / (24 * 7) + 4.0
    else:
        s_fit = 1.0 if D_out <= MK.CAST["static"]["ok"] else 0.5
        s_proc = MK.MAT["K438"]["cast"]
        weeks = 16.0
    if c["cooling"] == "керамика":
        s_proc = min(s_proc, MK.MAT["CMC"]["ceramic"])
    s_trl = (trl - 1) / 8
    s_time_m = 26 / (26 + weeks)
    make = 10 * (0.30 * s_fit + 0.25 * s_proc + 0.30 * s_trl + 0.15 * s_time_m)
    return dict(parts=n_parts, parts_list=parts, welds=welds, alignment_fits=align, printer=(seg[1] if seg else None),
                segments=(seg[0] if seg else None), mass_kg=mass, swap_h=swap, people=people, special_tools=tools,
                keyed=keyed, service_score=srv, valves_gas=valves_gas, valves_diesel=valves_dz, changeover=chg,
                interface=iface, plugplay_score=pp, trl_min=trl, make_weeks=weeks, make_score=make)


# =========================================================================== цели
def _tg(tid):
    return next(t for t in TARGETS["targets"] if t["id"] == tid)


def compliance(c, row):
    """соответствие целям targets.json: значение, предел, выполнено, нарушение (относительное, >= 0)"""
    g = {r["load"]: r for r in row["gas"]}
    dz = {r["load"]: r for r in row["diesel"]}
    out = {}

    def put(tid, val, lim, sense="max"):
        if val is None:
            out[tid] = dict(value=None, limit=lim, ok=False, viol=1.0)
            return
        v = (val - lim) / abs(lim) if sense == "max" else (lim - val) / abs(lim)
        out[tid] = dict(value=val, limit=lim, ok=bool(v <= 0), viol=max(0.0, v))
    nox_t, co_t = _tg("nox_gas"), _tg("co_gas")
    put("nox_gas", max(g[1.0]["NOx"], g[0.5]["NOx"]), nox_t["target"]["max"])
    put("co_gas", max(g[1.0]["CO"], g[0.5]["CO"]), co_t["target"]["max"])
    lr = max(g[0.3]["NOx"] / nox_t["also"]["max"], g[0.3]["CO"] / co_t["also"]["max"],
             1.0 / max(min(r["margin_lbo"] or 0.0 for r in row["gas"]), 1e-6))
    put("load_range", lr, 1.0)
    put("fuel_valves", max(row["module"]["valves_gas"], row["module"]["valves_diesel"]),
        _tg("fuel_valves")["target"]["control_valves_per_fuel"])
    dd = _tg("diesel_dry")["target"]
    put("diesel_dry", max(max(r["NOx"] for r in row["diesel"]) / dd["NOx_max"],
                          max(r["CO"] for r in row["diesel"]) / dd["CO_max"]), 1.0)
    put("fuel_changeover", row["module"]["changeover"], 0.7, "min")
    put("hydrogen", row["flashback"]["H2_vol_max"], _tg("hydrogen")["target"]["H2_vol_min"], "min")
    put("pressure_loss", row["dp"]["dp_total"], _tg("pressure_loss")["target"]["max"])
    put("exit_uniformity", max(row["pf"]["PF"] / _tg("exit_uniformity")["target"]["PF_max"],
                               row["pf"]["RTDF"] / _tg("exit_uniformity")["target"]["RTDF_max"]), 1.0)
    wl = max(row["wall_gas"]["T_w1"], row["wall_diesel"]["T_w1"]) / row["wall_gas"]["T_lim"]
    if c["scheme"] == "CAT":
        wl = max(wl, 30000 / 8000)           # ресурс каталитического блока - оценка 8000 ч
    put("wall_life", wl, 1.0)
    pc = _tg("parts_count")["target"]
    put("parts_count", max(row["module"]["parts"] / pc["parts_max"], 1.0 + row["module"]["welds"] +
                           row["module"]["alignment_fits"]), 1.0)
    cs = _tg("cartridge_swap")["target"]
    put("cartridge_swap", max(row["module"]["swap_h"] / cs["time_h_max"], 1.0 + row["module"]["special_tools"]), 1.0)
    put("combustion_dynamics", row["dyn_index"], 1.0)
    hard = [t for t in TARGETS["funnel_interface"]["hard"] if t in out]
    return out, hard


# =========================================================================== оценка варианта
def evaluate(i, n, kv="KV-3", seed=CS.SEED):
    """полная оценка варианта уровня 0: строка результата (без исключений наружу)"""
    t0 = time.perf_counter()
    c = CS.concept(i, n, kv, seed)
    row = dict(idx=i, concept=c, status="ok", reject=[])
    try:
        bc = boundary()
        d = CS.derive(c, bc)
        row["derived"] = {k: v for k, v in d.items() if k != "zones"}
        row["zones"] = d["zones"]
        bad = CS.check(c, d, bc)
        if bad:
            row.update(status="отсев", reject=[dict(code=b[0], value=b[1], limit=b[2], why=b[3]) for b in bad])
            row["t_s"] = time.perf_counter() - t0
            return _rnd(row)
        fb = flashback(c, d, bc)
        row["flashback"] = fb
        gas = [gas_network(c, d, bc, li) for li in range(len(LOADS))]
        row["gas"] = gas
        g100 = gas[0]
        rej = []
        if (g100["margin_lbo"] or 0.0) < 1.0:
            rej.append(dict(code="ПОГАСАНИЕ_100", value=g100["margin_lbo"], limit=1.0,
                            why="стабилизирующая зона гаснет на 100 % (время пребывания меньше времени погасания)"))
        if fb["margin_CH4"] is not None and fb["margin_CH4"] < 1.0:
            rej.append(dict(code="ПРОСКОК", value=fb["margin_CH4"], limit=1.0, why="проскок пламени в смеситель на метане"))
        dz = [diesel_network(c, d, bc, li, gas[li]) for li in DIESEL_LOADS]
        row["diesel"] = dz
        ig = [r.get("ign_margin") for r in dz if r.get("ign_margin") is not None]
        if ig and min(ig) < 1.0:
            rej.append(dict(code="САМОВОСПЛАМЕНЕНИЕ", value=min(ig), limit=1.0,
                            why="дизель самовоспламеняется в смесителе (задержка меньше времени смешения)"))
        far_pz = (bc["loads"][0]["W_fuel"] / (c["a_pz"] * bc["loads"][0]["W_air"]))
        row["wall_gas"] = wall_temperature(c, d, bc, g100["T_flame"], min(far_pz, bc["far_st_CH4"]), 1.0)
        far_dz = bc["loads"][0]["W_fuel_diesel"] / (c["a_pz"] * bc["loads"][0]["W_air"])
        row["wall_diesel"] = wall_temperature(c, d, bc, dz[0]["T_hot"], min(far_dz, bc["far_st_diesel"]),
                                              336.0 / H_FUEL_DIESEL ** 2)
        row["pf"] = pattern_factor(c, d, bc, g100)
        row["dp"] = pressure_loss(c, d, bc, g100)
        row["module"] = module_models(c, d, bc, fb, d, row["wall_gas"])
        mlbo = max(min(r["margin_lbo"] or 0.0 for r in gas), 1e-3)
        row["dyn_index"] = DYN_F[c["scheme"]] * (1.5 / min(mlbo, 3.0)) * (0.03 / c["dp_inj"]) ** 0.5 * \
            (1.5 if c["swirl"] > 0.8 else 1.0)
        row["targets"], hard = compliance(c, row)
        row["hard_ok"] = all(row["targets"][t]["ok"] for t in hard)
        row["n_hard_fail"] = sum(not row["targets"][t]["ok"] for t in hard)
        if rej:
            row.update(status="отсев", reject=rej)
        row["feed"] = dict(dp_comb=row["dp"]["dp_total"], D_comb=d["D_env"], L_comb=d["L_mod"],
                           m_comb=row["module"]["mass_kg"], n_inj=c["n_el"], parts=row["module"]["parts"],
                           NOx_100=gas[0]["NOx"], NOx_50=gas[1]["NOx"], NOx_30=gas[2]["NOx"],
                           CO_100=gas[0]["CO"], CO_50=gas[1]["CO"], CO_30=gas[2]["CO"], PF=row["pf"]["PF"],
                           RTDF=row["pf"]["RTDF"], T_metal_max_K=max(row["wall_gas"]["T_w1"], row["wall_diesel"]["T_w1"]),
                           fuel_valves=max(row["module"]["valves_gas"], row["module"]["valves_diesel"]),
                           H2_vol_max=fb["H2_vol_max"], diesel_NOx_100=dz[0]["NOx"], diesel_NOx_50=dz[1]["NOx"],
                           cartridge_swap_h=row["module"]["swap_h"], V_L=d["V_L"])
    except Exception as e:      # вариант не теряется: статус «ошибка» с причиной
        row.update(status="ошибка", error=f"{type(e).__name__}: {str(e)[:300]}")
    row["t_s"] = time.perf_counter() - t0
    return _rnd(row)


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--bc", action="store_true", help="пересчитать границы и таблицы (comb_bc_results.json)")
    ap.add_argument("--one", type=int, help="оценить один вариант (индекс) выборки --n")
    ap.add_argument("--n", type=int, default=300)
    ap.add_argument("--workers", type=int, default=4)
    a = ap.parse_args()
    if a.bc:
        if os.path.exists(BC_FILE):
            os.remove(BC_FILE)
        b = boundary(a.workers)
        print("границы:", BC_FILE, b["runtime_s"], "с")
    if a.one is not None:
        print(json.dumps(evaluate(a.one, a.n), ensure_ascii=False, indent=1)[:6000])
