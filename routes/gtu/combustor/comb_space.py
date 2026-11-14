"""Воронка камеры сгорания двигателя-заменителя KV, уровень 0: пространство поиска, генератор, геометрия, допустимость.

Задача владельца (14.11.2026): переработать горение, форсунки и каналы двигателя KV перебором на AWS - «круче и проще
всего существующего». Приоритеты: минимум деталей (главное), Plug&Play в пакет Solar/Siemens (один клапан на топливо),
монтаж «хоть уборщиком» (камера - сменный картридж), затем КПД (потери давления) и размер. Цели камеры - targets.json.

Что здесь:
  * оси (CAT - дискретные, CONT - непрерывные) и схемы (SCHEMES - допустимые значения дискретных осей для схемы);
  * генератор: скремблированная последовательность Соболя (scipy.stats.qmc) с фиксированным зерном SEED; вариант i
    выдаётся по индексу (concept); часть k из S получает индексы i % S == k (shard_indices) - без пересечений и пропусков;
    отпечаток матрицы (fingerprint) пишется в каждую часть, слияние сверяет его;
  * derive: распределение воздуха, гидравлика смесителя и топливных струй, неравномерность смешения, объёмы зон по
    временам пребывания, размеры жаровой трубы и модуля, доля воздуха на охлаждение по пористости и перепаду -
    без кинетики (только равновесие Cantera для плотностей); вход comb_crn.evaluate;
  * check: правила допустимости до кинетики (код правила, величина, предел).

Числа моделей ниже - «оценка», если рядом не указан источник. Все константы собраны в таблицах вверху файла.
Запуск: модуль библиотечный; самопроверка - python comb_space.py --selfcheck
"""
import argparse, hashlib, inspect, json, math, os, sys, warnings

import numpy as np
from scipy.stats import qmc

HERE = os.path.dirname(os.path.abspath(__file__))
SEED = 20261114
DATE = "2026-11-14"

# =========================================================================== дискретные оси
SCHEME_NAMES = {
    "LPP": "бедное предварительное смешение (малоэмиссионное премикс)",
    "MM": "микросмешение",
    "FLOX": "MILD/FLOX - струйная стабилизация с рециркуляцией продуктов",
    "RQL": "богатое - гашение - бедное (RQL)",
    "TVC": "захваченный вихрь",
    "CAT": "каталитическая (штраф готовности)",
    "RDC": "вращающаяся детонация (штраф готовности)",
}
CAT = {
    "scheme": list(SCHEME_NAMES),
    "chamber": ["выносная", "трубчатые", "кольцевая", "противоточная"],
    "injector": ["закрутка", "струи", "микросмесители"],
    "diesel": ["центробежная", "пневматическая", "испарительная"],
    "cooling": ["плёночное", "отверстия", "решётка", "керамика"],
    "make": ["печатный картридж", "литьё", "гибрид"],
    "passive": ["нет", "гидравлика", "струи"],
}
CAT_NOTE = {
    "chamber": "выносная - одна трубчатая камера вне двигателя на фланце (картридж); трубчатые - 6 жаровых труб "
               "(как SGT-300/400, STATE_OF_ART.md, разд. 4); кольцевая - встроенная; противоточная - кольцевая вокруг турбины",
    "injector": "закрутка - завихритель с топливными струями в лопатках; струи - топливо струями в сносящий воздух "
                "(каналы без закрутки); микросмесители - сотни-тысячи элементов с одной-двумя струями",
    "diesel": "центробежная - простая центробежная форсунка (давление топлива); пневматическая - распыл воздухом смесителя "
              "(предплёночная, LDI); испарительная - испарение и смешение до зоны горения (риск самовоспламенения)",
    "cooling": "плёночное - щели; отверстия - многодырчатая эффузия; решётка - печатная двойная стенка с решёткой, "
               "весь воздух камеры идёт через неё последовательно (перепад dp_wall); керамика - вкладыш CMC",
    "passive": "деление топлива по нагрузке без клапанов (targets_basis: одним клапаном зона пламени живёт только на "
               "83-100 %): гидравлика - два контура от одного клапана, пилот через отверстие (расход ~ ΔP^0,5), основной "
               "через ламинарные капилляры (~ ΔP); струи - самостратификация: проникновение струй ~ расходу топлива",
}

# Готовность (TRL) - оценка по STATE_OF_ART.md, разд. 4 (таблица «Технологии горения», TRL - «наша оценка»)
TRL_SCHEME = {"LPP": 9, "MM": 7, "FLOX": 5, "RQL": 5, "TVC": 5, "CAT": 6, "RDC": 4}
TRL_SCHEME_NOTE = ("STATE_OF_ART.md разд. 4 (оценка): LPP 9; микросмешение 7-8; FLOX 5-6; RQL 5-6 в промышленности; "
                   "TVC 5-6; каталитическая была 8 и снята с рынка - принято 6; детонационная 4-5")
TRL_DIESEL = {"центробежная": 9, "пневматическая": 8, "испарительная": 5}   # LDI 7-8 (Titan 250, STATE_OF_ART разд. 4); оценка
TRL_COOL = {"плёночное": 9, "отверстия": 9, "решётка": 5, "керамика": 6}      # make_model.TRL: печатное охлаждение 5, керамика статика 6
TRL_MAKE = {"печатный картридж": 7, "литьё": 9, "гибрид": 7}                 # make_model.TRL: печать статики 7 (оценка)

# Допустимые значения дискретных осей по схеме (как ARCH воронки двигателя): сэмпл попадает только в допустимое
_ALL = {k: list(v) for k, v in CAT.items() if k != "scheme"}
SCHEMES = {
    "LPP": dict(),
    "MM": dict(injector=["микросмесители"], diesel=["пневматическая", "испарительная"]),
    "FLOX": dict(injector=["струи", "микросмесители"]),
    "RQL": dict(injector=["закрутка", "струи"], passive=["нет"]),
    "TVC": dict(injector=["закрутка", "струи", "микросмесители"]),
    "CAT": dict(chamber=["выносная", "трубчатые", "кольцевая"], diesel=["испарительная"], passive=["нет"]),
    "RDC": dict(chamber=["выносная", "кольцевая"], injector=["струи", "микросмесители"],
                diesel=["пневматическая", "испарительная"], passive=["нет"]),
}

# =========================================================================== непрерывные оси
# смесительные элементы по типу форсунки: n_el - число элементов на всю камеру (лог-равномерно), k_je - топливных струй
# на элемент, d_jet - диаметр топливной струи, мм (лог-равномерно), swirl - число крутки, l_mix - длина смешения в
# калибрах элемента, d_air - допустимый эквивалентный диаметр воздушного канала элемента, мм (правило ЭЛЕМЕНТ). Оценка.
INJ = {
    "закрутка": dict(n_el=(1, 24), k_je=(4, 24), d_jet=(0.8, 4.0), swirl=(0.4, 1.2), l_mix=(2.0, 8.0), d_air=(15.0, 260.0)),
    "струи": dict(n_el=(8, 400), k_je=(1, 8), d_jet=(0.5, 3.0), swirl=(0.0, 0.4), l_mix=(3.0, 15.0), d_air=(3.0, 60.0)),
    "микросмесители": dict(n_el=(200, 8000), k_je=(1, 2), d_jet=(0.3, 1.2), swirl=(0.0, 0.2), l_mix=(1.0, 6.0),
                           d_air=(0.6, 8.0)),
}
# по схеме: a_pz - доля воздуха через зону горения (RQL - богатая зона), a_aux - вспомогательная доля (смысл - AUX_KIND),
# tau_pz - время зоны пламени, мс (лог), dp_inj - перепад на смесителе, доля p3. Оценка.
SCH = {
    "LPP": dict(a_pz=(0.35, 0.85), a_aux=(0.05, 0.25), tau_pz=(0.3, 3.0), dp_inj=(0.015, 0.05)),
    "MM": dict(a_pz=(0.35, 0.85), a_aux=(0.05, 0.25), tau_pz=(0.3, 3.0), dp_inj=(0.02, 0.06)),
    "FLOX": dict(a_pz=(0.40, 0.90), a_aux=(0.05, 0.25), tau_pz=(0.5, 4.0), dp_inj=(0.02, 0.06)),
    "RQL": dict(a_pz=(0.08, 0.22), a_aux=(0.25, 0.55), tau_pz=(0.5, 4.0), dp_inj=(0.015, 0.05)),
    "TVC": dict(a_pz=(0.30, 0.80), a_aux=(0.04, 0.15), tau_pz=(0.3, 3.0), dp_inj=(0.015, 0.05)),
    "CAT": dict(a_pz=(0.45, 0.90), a_aux=(0.0, 0.0), tau_pz=(0.5, 3.0), dp_inj=(0.015, 0.05)),
    "RDC": dict(a_pz=(0.25, 0.60), a_aux=(0.0, 0.0), tau_pz=(0.05, 0.3), dp_inj=(0.03, 0.15)),
}
AUX_KIND = {"LPP": "pilot", "MM": "pilot", "FLOX": "pilot", "TVC": "cavity", "RQL": "quench", "CAT": None, "RDC": None}
# пористость (доля открытой площади стенки) по типу охлаждения - оценка
POROSITY = {"плёночное": (0.010, 0.040), "отверстия": (0.005, 0.030), "решётка": (0.001, 0.010), "керамика": (0.0, 0.005)}
CONT = {
    "n_el": "число смесительных элементов на камеру (лог-равномерно в диапазоне INJ[форсунка])",
    "k_je": "топливных струй на элемент (INJ); число струй n_jet = n_el·k_je",
    "d_jet": "диаметр топливной струи, мм (лог-равномерно, INJ)",
    "swirl": "число крутки (закрутка 0,4-1,2; струи - наклон/крутка 0-0,4; микросмесители 0-0,2)",
    "a_pz": "доля воздуха через зону горения (SCH; у RQL - богатая зона)",
    "a_aux": "вспомогательная доля воздуха: пилот (доля воздуха зоны горения), полость TVC, воздух гашения RQL (SCH)",
    "k_pass": "сила пассивного деления: гидравлика - доля топлива пилота на 100 % (0,05-0,40); струи - показатель "
              "самостратификации k (0,2-1,0); TVC - доля топлива полости на 100 % (0,05-0,30)",
    "tau_pz": "время пребывания в зоне пламени, мс (лог, SCH)",
    "tau_bo": "время догорания до разбавления, мс (лог 2-20; targets.json guides - 10-20 мс)",
    "tau_dil": "время в зоне разбавления, мс (1-8)",
    "LD": "отношение длины жаровой трубы к диаметру (кольцевая - к высоте канала): 1,5-4,0; длина и диаметр - из объёма зон",
    "dp_inj": "перепад давления на смесителе (= на стенке при параллельном охлаждении), доля p3 (SCH)",
    "dp_wall": "перепад на печатной решётке стенки, включённой последовательно, доля p3: 0,3-2 % (только решётка)",
    "porosity": "пористость охлаждения - доля открытой площади стенки (POROSITY по типу охлаждения)",
    "l_mix": "длина смешения в калибрах элемента (INJ)",
}
AXES = list(CAT) + list(CONT)


def _pick(u, vals):
    return vals[min(int(u * len(vals)), len(vals) - 1)]


def _lerp(u, lo, hi):
    return lo + (hi - lo) * float(u)


def _loglerp(u, lo, hi):
    return math.exp(math.log(lo) + (math.log(hi) - math.log(lo)) * float(u))


# =========================================================================== генератор
_DM = {}


def _rng_kw(fn, seed):
    return {"rng": seed} if "rng" in inspect.signature(fn).parameters else {"seed": seed}


def design_matrix(n, seed=SEED):
    """матрица n x len(AXES): первые n точек из 2^m скремблированной последовательности Соболя"""
    key = (int(n), int(seed))
    if key not in _DM:
        m = max(1, math.ceil(math.log2(max(n, 2))))
        eng = qmc.Sobol(len(AXES), scramble=True, **_rng_kw(qmc.Sobol, seed))
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            _DM[key] = eng.random_base2(m)[:n]
    return _DM[key]


def fingerprint(n, seed=SEED):
    return hashlib.sha256(np.round(design_matrix(n, seed), 12).astype("<f8").tobytes()).hexdigest()[:16]


def shard_indices(n, shard, shards):
    """единственное разбиение: часть shard из shards получает индексы i с i % shards == shard"""
    if not (shards >= 1 and 0 <= shard < shards and n >= 1):
        raise ValueError(f"неверные n={n}, shard={shard}, shards={shards}")
    return list(range(shard, n, shards))


def code_hash():
    """отпечаток кода уровня 0 (comb_space, comb_crn): части с разным кодом слияние не примет"""
    h = hashlib.sha256()
    for f in ("comb_space.py", "comb_crn.py"):
        h.update(open(os.path.join(HERE, f), "rb").read())
    return h.hexdigest()[:16]


def concept(i, n, kv="KV-3", seed=SEED):
    """вариант с индексом i выборки объёма n: значения всех осей (детерминированно по i, n, seed)"""
    if not 0 <= i < n:
        raise IndexError(i)
    u = dict(zip(AXES, design_matrix(n, seed)[i]))
    s = _pick(u["scheme"], CAT["scheme"])
    allowed = {k: SCHEMES[s].get(k, _ALL[k]) for k in _ALL}
    c = dict(idx=int(i), n=int(n), seed=int(seed), kv=kv, scheme=s)
    c["chamber"] = _pick(u["chamber"], allowed["chamber"])
    c["injector"] = _pick(u["injector"], allowed["injector"])
    dz = allowed["diesel"]
    if c["injector"] == "микросмесители":
        dz = [x for x in dz if x != "центробежная"] or dz      # центробежную форсунку в микроэлемент не поставить
    c["diesel"] = _pick(u["diesel"], dz)
    c["cooling"] = _pick(u["cooling"], allowed["cooling"])
    mk = allowed["make"]
    if c["injector"] == "микросмесители" or c["cooling"] == "решётка":
        mk = [x for x in mk if x != "литьё"]                  # отверстия < 1 мм и решётка не отливаются
    c["make"] = _pick(u["make"], mk)
    c["passive"] = _pick(u["passive"], allowed["passive"])
    I, S = INJ[c["injector"]], SCH[s]
    c["n_el"] = int(round(_loglerp(u["n_el"], *I["n_el"])))
    c["k_je"] = int(round(_lerp(u["k_je"], *I["k_je"])))
    c["d_jet_mm"] = _loglerp(u["d_jet"], *I["d_jet"])
    c["swirl"] = _lerp(u["swirl"], *I["swirl"])
    c["a_pz"] = _lerp(u["a_pz"], *S["a_pz"])
    c["a_aux"] = _lerp(u["a_aux"], *S["a_aux"]) if AUX_KIND[s] else 0.0
    if s == "TVC":
        c["k_pass"] = _lerp(u["k_pass"], 0.05, 0.30)
    elif c["passive"] == "гидравлика":
        c["k_pass"] = _lerp(u["k_pass"], 0.05, 0.40)
    elif c["passive"] == "струи":
        c["k_pass"] = _lerp(u["k_pass"], 0.2, 1.0)
    else:
        c["k_pass"] = 0.0
    c["tau_pz_ms"] = _loglerp(u["tau_pz"], *S["tau_pz"])
    c["tau_bo_ms"] = _loglerp(u["tau_bo"], 2.0, 20.0)
    c["tau_dil_ms"] = _lerp(u["tau_dil"], 1.0, 8.0)
    c["LD"] = _lerp(u["LD"], 1.5, 4.0)
    c["dp_inj"] = _lerp(u["dp_inj"], *S["dp_inj"])
    c["dp_wall"] = _lerp(u["dp_wall"], 0.003, 0.02) if c["cooling"] == "решётка" else 0.0
    c["porosity"] = _lerp(u["porosity"], *POROSITY[c["cooling"]])
    c["l_mix"] = _lerp(u["l_mix"], *I["l_mix"])
    for k, v in list(c.items()):
        if isinstance(v, float):
            c[k] = round(v, 6)
    return c


# =========================================================================== геометрия и гидравлика (без кинетики)
# Константы - оценка, если не указан источник
CD_AIR = 0.70          # коэффициент расхода воздушных каналов смесителя и отверстий стенки
CD_FUEL = 0.80         # то же топливных отверстий
M3 = 0.20              # число Маха за компрессором перед диффузором камеры
K_DIFF = {"кольцевая": 0.40, "трубчатые": 0.55, "выносная": 0.70, "противоточная": 0.55}   # потери диффузора/подвода в q3
DP_DUCT = {"выносная": 0.008, "трубчатые": 0.004, "кольцевая": 0.0, "противоточная": 0.004}  # переходный канал / поворот, доля p3
N_CAN = {"выносная": 1, "трубчатые": 6, "кольцевая": 1, "противоточная": 1}
U_ANN = 30.0           # скорость в кольцевом зазоре (кожух - жаровая труба), м/с
DM_ANN = 0.88          # средний диаметр кольцевой жаровой трубы / D2 компрессора (= 0,8·1,1·Dc как funnel_space, оценка той воронки)
DM_REV = 1.25          # то же противоточной (вокруг турбины)
D_ENV_MAX = {"кольцевая": 1.9, "противоточная": 1.9, "трубчатые": 2.2}   # наружный диаметр / D2 (kD_case 1,9 - gtu_family)
CAN_MAX = dict(D=1.2, L=3.0)   # выносная камера: наибольший диаметр кожуха и длина модуля, м (укрытие пакета; сверить gtu-dropin)
S_FLOOR, S0, K_MIX = 0.01, 0.5, 2.5   # неравномерность смешения: S = S_FLOOR + S0·f_C·exp(-l_mix/K_MIX)·(1 - 0,3·(крутка-0,4)+)
C_OPT = 2.5            # оптимум C = (S/H)·J^0,5 для одностороннего вдува (Holdeman, PECS 19, 1993 - по памяти, не сверено)
D_JET_MIN = {"печатный картридж": 0.3, "гибрид": 0.3, "литьё": 0.8}   # наименьшее топливное отверстие, мм (LPBF / сверление после литья)
WALL_T = {"металл": 0.0015, "керамика": 0.003}   # толщина стенки жаровой трубы, м
GAS_R_CH4, GAS_G_CH4, T_FUEL = 518.28, 1.31, 288.15   # R, γ метана; температура топлива как gtu_cycle (T0)


def _orifice_p0(mdot, A, p_down, R=GAS_R_CH4, g=GAS_G_CH4, T0=T_FUEL, Cd=CD_FUEL):
    """давление перед отверстиями (изоэнтропическое истечение, докритика / запирание) для расхода mdot"""
    crit = (2 / (g + 1)) ** (g / (g - 1))

    def m_of(p0):
        r = p_down / p0
        if r <= crit:
            return Cd * A * p0 * math.sqrt(g / (R * T0)) * (2 / (g + 1)) ** ((g + 1) / (2 * (g - 1)))
        return Cd * A * p0 * math.sqrt(2 * g / ((g - 1) * R * T0) * (r ** (2 / g) - r ** ((g + 1) / g)))
    lo, hi = p_down * (1 + 1e-9), p_down * 2.0
    while m_of(hi) < mdot:
        hi *= 2.0
        if hi > 1e9:
            return float("inf"), 0.0, 0.0
    for _ in range(60):
        mid = 0.5 * (lo + hi)
        if m_of(mid) < mdot:
            lo = mid
        else:
            hi = mid
    p0 = hi
    r = max(p_down / p0, 1e-9)
    u = math.sqrt(max(0.0, 2 * g / (g - 1) * R * T0 * (1 - r ** ((g - 1) / g))))   # полное расширение до p_down
    Te = T0 * r ** ((g - 1) / g)
    return p0, u, p_down / (R * Te)


def derive(c, bc):
    """геометрия, воздух, смеситель, топливные струи, неравномерность - по точке 100 % (bc['loads'][0])."""
    L0 = bc["loads"][0]
    p3, T3, W, Wf = L0["p3_Pa"], L0["T3_K"], L0["W_air"], L0["W_fuel"]
    rho3, mu3 = L0["rho3"], L0["mu3"]
    far_st, phi_ov = bc["far_st_CH4"], L0["phi_ov"]
    s, ak = c["scheme"], AUX_KIND[c["scheme"]]
    d = {}
    # --- воздух через смесители и скорость в них
    a_mix = c["a_pz"]
    u_a = math.sqrt(2 * c["dp_inj"] * p3 / rho3)
    A_air = a_mix * W / (CD_AIR * math.sqrt(2 * rho3 * c["dp_inj"] * p3))
    n_el = c["n_el"]
    d_air = math.sqrt(4 * A_air / (math.pi * n_el))
    d.update(u_a=u_a, A_air_m2=A_air, d_air_mm=1e3 * d_air, tau_mix_ms=1e3 * c["l_mix"] * d_air / u_a)
    # --- топливные струи: давление топлива, импульс, критерий Холдемана
    n_jet = n_el * c["k_je"]
    A_f = n_jet * math.pi * (1e-3 * c["d_jet_mm"]) ** 2 / 4
    p_down = p3 * (1 - c["dp_inj"])
    p0, u_f, rho_f = _orifice_p0(Wf if s != "TVC" else Wf * (1 - c["k_pass"]), A_f, p_down)
    J = rho_f * u_f ** 2 / (rho3 * u_a ** 2) if u_f > 0 else 0.0
    SH = 2 * math.pi / c["k_je"]                              # шаг струй / полувысота канала (струи по окружности канала)
    C = SH * math.sqrt(J) if J > 0 else 0.0
    fC = 1 + 0.5 * math.log(max(C, 1e-3) / C_OPT) ** 2
    S = S_FLOOR + S0 * fC * math.exp(-c["l_mix"] / K_MIX) * (1 - 0.3 * max(0.0, c["swirl"] - 0.4))
    d.update(n_jet=n_jet, p_fuel_Pa=p0, p_fuel_ratio=p0 / p3, u_fuel=u_f, J=J, holdeman_C=C, f_C=fC, S_unmix=min(S, 0.6))
    # --- воздух охлаждения и разбавления: оценка площади стенки по объёму (равновесные плотности)
    zones = zone_masses(c, bc)
    V = sum(z["m"] * z["tau"] / z["rho"] for z in zones)
    geo = geometry(c, bc, V)
    geo["L_mod"] += c["l_mix"] * d_air                       # длина смесительных каналов
    dpL = c["dp_inj"] * p3
    m_cool = CD_AIR * c["porosity"] * geo["A_wall"] * math.sqrt(2 * rho3 * dpL)
    a_cool = m_cool / W
    a_aux_out = c["a_aux"] if ak in ("cavity", "quench") else 0.0
    a_dil = 1.0 - c["a_pz"] - a_aux_out - a_cool
    d.update(zones=[dict(name=z["name"], m=z["m"], tau_ms=1e3 * z["tau"], T_K=z["T"], V_L=1e3 * z["m"] * z["tau"] / z["rho"])
                    for z in zones], V_L=1e3 * V, a_cool=a_cool, a_aux_out=a_aux_out, a_dil=a_dil, **geo)
    # --- потери давления (без горячих - их добавляет comb_crn по температурам сети)
    q3_over_p3 = 0.5 * 1.4 * M3 ** 2 * (1 + 0.2 * M3 ** 2) ** -3.5   # q/p для числа Маха M3 (γ=1,4)
    d["dp_diff"] = K_DIFF[c["chamber"]] * q3_over_p3
    d["dp_duct"] = DP_DUCT[c["chamber"]]
    d["dp_series_wall"] = c["dp_wall"] if c["cooling"] == "решётка" else 0.0
    # --- отношение перепада на жаровой трубе к скоростному напору по эталонному сечению (корреляция PF)
    q_ref = 0.5 * rho3 * (W / N_CAN[c["chamber"]] / (rho3 * geo["A_ref"])) ** 2
    d["dpL_qref"] = dpL / q_ref
    d["phi_pz_100"] = phi_ov / c["a_pz"] if s not in ("TVC",) else phi_ov * (1 - c["k_pass"]) / c["a_pz"]
    return d


def zone_masses(c, bc):
    """зоны камеры: массовый расход, время пребывания и плотность по равновесию (точка 100 %) - для объёма"""
    L0 = bc["loads"][0]
    W, Wf, far_st = L0["W_air"], L0["W_fuel"], bc["far_st_CH4"]
    eqT = bc["eq_T"]                 # функция phi -> (T, rho) по таблице равновесия
    s = c["scheme"]
    tp, tb, td = 1e-3 * c["tau_pz_ms"], 1e-3 * c["tau_bo_ms"], 1e-3 * c["tau_dil_ms"]
    T4, rho4 = L0["T4_K"], L0["rho4"]
    Z = []
    if s == "RQL":
        phi_r = (Wf / W) / c["a_pz"] / far_st
        Tr, rr = eqT(phi_r)
        Z.append(dict(name="богатая", m=c["a_pz"] * W + Wf, tau=tp, T=Tr, rho=rr))
        m_q = (c["a_pz"] + c["a_aux"]) * W + Wf
        phi_q = (Wf / W) / (c["a_pz"] + c["a_aux"]) / far_st
        Tq, rq = eqT(phi_q)
        Z.append(dict(name="гашение", m=m_q, tau=1e-3 * quench_time_ms(c, bc), T=Tq, rho=rq))
        Z.append(dict(name="бедная", m=m_q, tau=tb, T=Tq, rho=rq))
    elif s == "TVC":
        f_cav = c["k_pass"]
        phi_c = f_cav * Wf / W / c["a_aux"] / far_st
        Tc, rc = eqT(phi_c)
        Z.append(dict(name="полость", m=c["a_aux"] * W + f_cav * Wf, tau=tp, T=Tc, rho=rc))
        phi_m = (1 - f_cav) * Wf / W / c["a_pz"] / far_st
        Tm, rm = eqT(phi_m)
        Z.append(dict(name="основная", m=c["a_pz"] * W + (1 - f_cav) * Wf, tau=tp, T=Tm, rho=rm))
        phi_b = Wf / W / (c["a_pz"] + c["a_aux"]) / far_st
        Tb, rb = eqT(phi_b)
        Z.append(dict(name="догорание", m=(c["a_pz"] + c["a_aux"]) * W + Wf, tau=tb, T=Tb, rho=rb))
    else:
        phi_p = (Wf / W) / c["a_pz"] / far_st
        Tp, rp = eqT(phi_p)
        Z.append(dict(name="зона пламени", m=c["a_pz"] * W + Wf, tau=tp, T=Tp, rho=rp))
        Z.append(dict(name="догорание", m=c["a_pz"] * W + Wf, tau=tb, T=Tp, rho=rp))
    Z.append(dict(name="разбавление", m=W + Wf, tau=td, T=T4, rho=rho4))
    return Z


def quench_time_ms(c, bc):
    """время смешения воздуха гашения RQL: 3 пролёта струи до оси трубы (оценка), D - по объёму без зоны гашения"""
    L0 = bc["loads"][0]
    u_j = math.sqrt(2 * c["dp_inj"] * L0["p3_Pa"] / L0["rho3"])
    W, Wf, far_st, eqT = L0["W_air"], L0["W_fuel"], bc["far_st_CH4"], bc["eq_T"]
    phi_r = (Wf / W) / c["a_pz"] / far_st
    Tr, rr = eqT(phi_r)
    V = (c["a_pz"] * W + Wf) * 1e-3 * c["tau_pz_ms"] / rr
    Dh = geometry(c, bc, V, quick=True)["D_h"]
    return 1e3 * 3 * (Dh / 2) / u_j


def geometry(c, bc, V, quick=False):
    """размеры жаровой трубы и модуля по объёму V (м3) и L/D; площадь стенки; эталонное сечение"""
    L0 = bc["loads"][0]
    W, rho3 = L0["W_air"], L0["rho3"]
    ch, N, LD = c["chamber"], N_CAN[c["chamber"]], c["LD"]
    D2 = bc["D2_m"]
    g = {}
    A_an_tot = W * max(0.2, 1 - c["a_pz"]) / (rho3 * U_ANN)          # кольцевой зазор: воздух мимо купола
    if ch in ("выносная", "трубчатые"):
        Vc = V / N
        D = (4 * Vc / (math.pi * LD)) ** (1 / 3)
        L = LD * D
        A_an = A_an_tot / N
        D_case = math.sqrt(D ** 2 + 4 * A_an / math.pi) + 2 * 0.01
        g.update(D_L=D, L_L=L, D_h=D, A_L=math.pi * D ** 2 / 4, A_wall=N * (math.pi * D * L + math.pi * D ** 2 / 4),
                 A_ref=math.pi * D_case ** 2 / 4, D_case=D_case, H_L=None, D_m=None)
        if ch == "трубчатые":
            D_ring = max(1.0 * D2, N * D_case * 1.15 / math.pi)
            g["D_env"] = D_ring + D_case
        else:
            g["D_env"] = D_case
    else:
        Dm = (DM_ANN if ch == "кольцевая" else DM_REV) * D2
        H = math.sqrt(V / (math.pi * Dm * LD))
        L = LD * H
        h_an = A_an_tot / (2 * math.pi * Dm)                       # два зазора - внутренний и наружный
        g.update(D_L=None, L_L=L, H_L=H, D_m=Dm, D_h=2 * H, A_L=math.pi * Dm * H,
                 A_wall=2 * math.pi * Dm * L + math.pi * Dm * H, A_ref=math.pi * Dm * (H + 2 * h_an),
                 D_case=Dm + H + 2 * h_an + 2 * 0.01, D_env=Dm + H + 2 * h_an + 2 * 0.01)
    if quick:
        return g
    g["L_mod"] = g["L_L"] + 0.10 + (0.25 if ch == "выносная" else 0.0)   # купол и фланцы (оценка); + смеситель в derive
    return g


def check(c, d, bc):
    """правила допустимости до кинетики: список (код, величина, предел, пояснение)"""
    bad = []
    I = INJ[c["injector"]]
    if d["a_dil"] < 0:
        bad.append(("ВОЗДУХ", d["a_dil"], 0.0, "воздуха на разбавление не остаётся: a_pz + a_aux + a_cool > 1"))
    if d["a_cool"] > 0.45:
        bad.append(("ОХЛАЖДЕНИЕ", d["a_cool"], 0.45, "на охлаждение уходит больше 45 % воздуха"))
    if not I["d_air"][0] <= d["d_air_mm"] <= I["d_air"][1]:
        bad.append(("ЭЛЕМЕНТ", d["d_air_mm"], I["d_air"], "канал элемента вне диапазона типа форсунки, мм"))
    if c["d_jet_mm"] < D_JET_MIN[c["make"]]:
        bad.append(("СТРУЯ", c["d_jet_mm"], D_JET_MIN[c["make"]], "топливное отверстие мельче предела способа изготовления, мм"))
    if d["p_fuel_ratio"] > 3.0:
        bad.append(("ТОПЛИВО", d["p_fuel_ratio"], 3.0, "давление топливного газа выше 3·p3"))
    if c["diesel"] == "центробежная" and c["n_el"] > 48:
        bad.append(("ДИЗЕЛЬ_ФОРСУНКИ", c["n_el"], 48, "центробежных форсунок больше 48 - малый расход на форсунку"))
    if c["chamber"] == "трубчатые" and c["n_el"] < N_CAN["трубчатые"]:
        bad.append(("ЭЛЕМЕНТОВ_НА_ТРУБУ", c["n_el"], N_CAN["трубчатые"], "меньше одного смесителя на жаровую трубу"))
    if c["chamber"] == "выносная" and (d["D_env"] > CAN_MAX["D"] or d["L_mod"] > CAN_MAX["L"]):
        bad.append(("ГАБАРИТ", [round(d["D_env"], 3), round(d["L_mod"], 3)], [CAN_MAX["D"], CAN_MAX["L"]],
                    "выносная камера не помещается (диаметр, длина), м"))
    if c["chamber"] in D_ENV_MAX and d["D_env"] > D_ENV_MAX[c["chamber"]] * bc["D2_m"]:
        bad.append(("ГАБАРИТ", round(d["D_env"], 3), round(D_ENV_MAX[c["chamber"]] * bc["D2_m"], 3),
                    "наружный диаметр камеры больше предела корпуса, м"))
    return bad


def _selfcheck():
    n = 64
    a = [concept(i, n) for i in range(n)]
    b = [concept(i, n) for i in range(n)]
    assert a == b
    allidx = sorted(sum((shard_indices(n, k, 5) for k in range(5)), []))
    assert allidx == list(range(n))
    from collections import Counter
    print(json.dumps(dict(fingerprint=fingerprint(n), schemes=Counter(x["scheme"] for x in a),
                          injectors=Counter(x["injector"] for x in a)), ensure_ascii=False))


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--selfcheck", action="store_true")
    if ap.parse_args().selfcheck:
        _selfcheck()
