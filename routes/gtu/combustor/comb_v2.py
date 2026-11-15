"""Камера сгорания KV, уровень 0, v2 - решение владельца 15.11.2026: один клапан на топливо, без приводов, частичную
нагрузку 30-100 % держит пассивный элемент (ищется перебором).

Что нового против v1 (comb_space.py + comb_crn.py - НЕ меняются: на них посчитан прогон AWS 20261115T001745Z-f1da,
code_hash fae37b9aadf85961; v2 импортирует их и добавляет слой):
  * четыре точки нагрузки 100 / 75 / 50 / 30 % (75 % добавлена тем же способом, что точки v1 - boundary_v2);
  * согласованное смешение: неравномерность S, проникновение струй и расслоение считаются на КАЖДОЙ нагрузке по
    отношению импульсов J этой нагрузки (расход топлива каждого контура, давления, скорость воздуха при неизменной
    геометрии), а не по точке 100 % (проверка скептика: v1 брал S при 100 % для всех нагрузок);
  * закон затухания неравномерности и расслоения при малом J (недопроникающие струи): длина затухания
    K_MIX·max(1, C_OPT/C)^BETA_LOWJ - ОЦЕНКА, проверять на уровне 1 (RANS); чувствительность к BETA_LOWJ - sens_v2.py;
  * гидравлика топлива одним клапаном: общий коллектор, контуры с разными законами сопротивления (сопло -
    изоэнтропическое истечение газа; пористая печатная вставка - закон Дарси - Форхгеймера для изотермического газа),
    разное давление за отверстиями (впрыск у входа канала смесителя или у его выхода) - доля топлива контуров на
    частичной нагрузке получается решением баланса расходов, а не задаётся;
  * пять пассивных элементов частичной нагрузки (PASSIVE_V2: гидравлический пилот, полость с вихрем TVC - схема,
    самостратификация струй, деление топлива каналами, горячая рециркуляция FLOX/эжектор) - сочетаются: деление топлива
    × самостратификация × рециркуляция (ALLOWED); у каждого своя часть сети реакторов Cantera на 30/50/75/100 %;
  * правило обратного тока (воздух в топливный коллектор) и множитель «жёсткости» топливной подачи в индексе пульсаций
    (малый перепад топлива на малой нагрузке связывает расход топлива с пульсациями давления) - чтобы схема со сдвигом
    давления не получала выигрыш даром;
  * цель «низкие выбросы 30-100 % на одном клапане» (lowem_30_100_1valve) - в дополнение к 13 целям targets.json;
  * каталитическая и детонационная схемы в v2 не входят: у них по два клапана на газ (v1: fuel_valves не выполнена ни
    разу из 27 083 допустимых), решению 2 они не отвечают;
  * исправления проверок скептиков 15.11.2026 (физика и код): пилот «гидравлика» - часть смесителей того же купола
    (время = tau_pz); эжектор - струйный насос в трубе не шире ячейки элемента; жёсткость подачи - по отклику расхода
    (закон контура); разброс сопротивления пористой вставки (ОБРАТНЫЙ_ТОК_РАЗБРОС); холостой ход 0 % (оценка,
    экстраполяция цикла) - обратный ток и устойчивость для load_range; времена зон при неизменном объёме; узлы
    неравномерности без прибавки топлива; запас по срыву ниже 1 - уточнённый; каналы - верхняя оценка с поджиганием
    группы B газом группы A (отбор - без неё); пометка «вне области облика»; сброс состояния Cantera перед вариантом.

Все числа модели - облик (0D + корреляции), «оценка» для отбора на уровень 1. AWS не вызывается.
Запуск: python comb_v2.py --bc (границы v2) | --targets (targets_v2.json) | --one I --n N | --selfcheck
"""
import argparse, hashlib, inspect, json, math, os, sys, time, warnings

import numpy as np
from scipy.optimize import brentq
from scipy.stats import qmc

HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in sys.path:
    sys.path.insert(0, HERE)
import comb_space as CS  # noqa: E402
import comb_crn as CR    # noqa: E402

DATE = "2026-11-15"
SEED = 20261115
LOADS = [1.0, 0.75, 0.5, 0.3]
BC_FILE = os.path.join(HERE, "comb_bc_v2_results.json")
TARGETS_FILE = os.path.join(HERE, "targets_v2.json")
NEW_TARGET = "lowem_30_100_1valve"

# =========================================================================== пространство поиска v2
SCHEME_NAMES = {k: CS.SCHEME_NAMES[k] for k in ("LPP", "MM", "FLOX", "TVC", "RQL")}
PASSIVE_V2 = ["нет", "гидравлика", "струи", "каналы", "рециркуляция"]   # элементы; полость TVC - схема TVC
PASSIVE_NOTE = {
    "нет": "один контур топлива, неизменная геометрия - опорный вариант (у TVC - полость с постоянной долей топлива: "
           "сопла полости и основного контура с одним законом и одним давлением за ними)",
    "гидравлика": "гидравлический пилот: часть смесителей (доля a_aux воздуха зоны; у TVC - полость) питается через "
                  "сопла прямо из коллектора и впрыскивает у выхода канала; основной контур - пористая печатная вставка "
                  "(Дарси - Форхгеймер, доля линейного члена x_lin) последовательно с соплами, впрыск у входа канала "
                  "(давление выше на off·dp_inj). При падении давления коллектора доля пилота растёт сама",
    "струи": "самостратификация: глубина проникновения струй ~ C = (S/H)·J^0,5; при малом расходе топливо остаётся в "
             "доле воздуха (C/C_OPT)^k_str, расслоение затухает по длине смешения медленнее при малом J (оценка)",
    "каналы": "пассивное деление топлива каналами: смесители двух групп (доля w_grp - группа A) от одного коллектора; "
              "A - сопла, впрыск у выхода канала; B - пористая вставка + сопла, впрыск у входа канала; на малой нагрузке "
              "топливо уходит в группу A (ступенчатость без клапана); группа B горит сама (отбор), верхняя оценка - "
              "с поджиганием горячим газом A (gas_cross)",
    "рециркуляция": "горячая рециркуляция (FLOX): воздушная струя элемента - эжектор, подсасывает газ догорания по "
                    "печатному каналу; кратность - меньшее из уноса Ricou - Spalding на длине трубы эжектора l_ej "
                    "калибров и струйного насоса в трубе не шире ячейки элемента (баланс импульса, потери возвратного "
                    "канала, выброс в жаровую трубу); топливо - на выходе эжектора",
}
# Готовность пассивного элемента - ОЦЕНКА (открытых оценок нет): деление по закону расхода и сдвигу давления в ГТУ не
# найдено (идея STATE_OF_ART.md, разд. 6.4 - «не считались»); самостратификация - только идея; эжекторная рециркуляция
# горячего газа в камере ГТУ - опыт FLOX DLR (5-6) без эжектора - 4.
TRL_PASSIVE = {"нет": 9, "гидравлика": 4, "струи": 3, "каналы": 4, "рециркуляция": 4}
# Элементы независимы физически и сочетаются: деление топлива (нет / гидравлика / каналы) × самостратификация (нет /
# струи) × горячая рециркуляция (нет / рециркуляция). TVC - полость и есть её элемент: деление (нет / гидравлика) ×
# струи; RQL - без пассивных элементов (опорная схема с одним клапаном).
ALLOWED = {
    "split": {"LPP": ["нет", "гидравлика", "каналы"], "MM": ["нет", "гидравлика", "каналы"],
              "FLOX": ["нет", "гидравлика", "каналы"], "TVC": ["нет", "гидравлика"], "RQL": ["нет"]},
    "strat": {"LPP": ["нет", "струи"], "MM": ["нет", "струи"], "FLOX": ["нет", "струи"], "TVC": ["нет", "струи"],
              "RQL": ["нет"]},
    "recirc": {"LPP": ["нет", "рециркуляция"], "MM": ["нет", "рециркуляция"], "FLOX": ["нет", "рециркуляция"],
               "TVC": ["нет"], "RQL": ["нет"]},
}
CAT_V2 = {"scheme": list(SCHEME_NAMES), "chamber": CS.CAT["chamber"], "injector": CS.CAT["injector"],
          "diesel": CS.CAT["diesel"], "cooling": CS.CAT["cooling"], "make": CS.CAT["make"],
          "split": ["нет", "гидравлика", "каналы"], "strat": ["нет", "струи"], "recirc": ["нет", "рециркуляция"]}


def passive_of(c):
    """строка сочетания пассивных элементов («нет», «гидравлика+струи», ...)"""
    act = [c[k] for k in ("split", "strat", "recirc") if c.get(k, "нет") != "нет"]
    return "+".join(act) if act else "нет"


def trl_passive(c):
    act = [c[k] for k in ("split", "strat", "recirc") if c.get(k, "нет") != "нет"]
    return min([TRL_PASSIVE[x] for x in act] or [9])
CONT_V2 = {
    "n_el": "число смесительных элементов (лог, CS.INJ)", "k_je": "топливных струй на элемент (CS.INJ)",
    "p_sup": "давление топливного коллектора при 100 %: положение (лог) в окне, где диаметры струй всех контуров в "
             "диапазоне типа форсунки и способа изготовления, пересечённом с 1,03-2,5·p3 (CS.P_FUEL_WIN)",
    "swirl": "крутка (CS.INJ)", "a_pz": "доля воздуха через зону горения (CS.SCH)",
    "a_aux": "пилот - доля воздуха зоны (0,05-0,40); TVC - воздух полости, доля всего воздуха (0,04-0,15)",
    "k_pass": "гидравлика и каналы - отношение phi пилота (группы A) к phi основного контура при 100 % (0,8-1,6); "
              "TVC - phi полости при 100 % (0,6-1,4)",
    "k_str": "показатель самостратификации: доля воздуха, встречающая топливо у входа канала, (C/C_OPT)^k_str (0,5-1,5)",
    "tau_pz": "время зоны пламени, мс (CS.SCH)",
    "tau_aux": "время полости TVC, мс (лог 0,5-5); у пилота «гидравлика» = tau_pz (часть смесителей того же купола)",
    "tau_bo": "время догорания, мс (лог 2-20)", "tau_dil": "время разбавления, мс (1-8)", "LD": "L/D жаровой трубы (1,5-4)",
    "dp_inj": "перепад на смесителе при 100 %, доля p3 (CS.SCH)", "dp_wall": "перепад на решётке стенки (0,3-2 %)",
    "porosity": "пористость охлаждения (CS.POROSITY)", "l_mix": "длина смешения, калибров (CS.INJ)",
    "x_lin": "доля линейного (дарсиевского) члена сопротивления пористой вставки при 100 % (0-0,95)",
    "f_por": "доля перепада основного контура на пористой вставке при 100 % (0,2-0,9)",
    "off": "сдвиг давления за соплами основного контура против пилотного: off·dp_inj (0 - оба у выхода канала, 1 - "
           "основной у входа канала)",
    "w_grp": "каналы: доля смесителей группы A (0,15-0,5)",
    "l_ej": "рециркуляция: длина трубы эжектора, калибров канала элемента (4-15)",
}
AXES_V2 = list(CAT_V2) + list(CONT_V2)

# Константы v2 - ОЦЕНКА, если не указан источник
BETA_LOWJ_NOMINAL = 1.0
BETA_LOWJ = BETA_LOWJ_NOMINAL   # затухание при недопроникновении: K = K_MIX·max(1, C_OPT/C)^BETA. Струи, не доходящие до середины
                       # канала, смешиваются турбулентностью канала, а не уносом струи (оценка; проверка - уровень 1 RANS)
DP_BACK_MIN = 0.01     # наименьший запас давления коллектора над давлением за соплами на любой нагрузке, доля p3: не
                       # меньше допустимого размаха пульсаций (targets.json combustion_dynamics 1 % p3), иначе воздух
                       # входит в топливный контур (обратный ток)
POR_TOL = 0.8          # сопротивление пористой вставки ниже расчётного на 20 % (разброс печати пористой структуры, оценка,
                       # не сверено): правило ОБРАТНЫЙ_ТОК проверяется и при нём (проверка скептика 15.11.2026)
DPF_REF = 0.03         # перепад топлива (доля p3), ниже которого расход топлива заметно откликается на пульсации
                       # давления: множитель индекса пульсаций (DPF_REF/dpf_min)^0,5 (оценка, по памяти, не сверено)
TAU_AUX = (0.5, 5.0)
SCH_AUX = {"pilot": (0.05, 0.40), "cavity": (0.04, 0.15)}   # пилот v2 шире v1 (0,05-0,25): на 30 % пилот должен
                                                               # нести почти всё топливо (оценка)
# Холостой ход без нагрузки при полной частоте (цель load_range: stable_pct [0, 100]; проверка скептика 15.11.2026):
# модель цикла не считает ниже доли топлива 0,005 (brentq в gtu_cycle), поэтому точка 0 % - ЛИНЕЙНАЯ ЭКСТРАПОЛЯЦИЯ
# p3, T3 и расхода топлива по точкам 50 и 30 % (оценка; вторым способом - квадратичная по четырём точкам, расхождение
# пишется в границы). На холостом ходу проверяются только устойчивость опорной зоны и обратный ток, не выбросы.
IDLE_LOAD = 0.0
# Область облика (оценка): за её пределами корреляции смешения экстраполированы - вариант помечается, проверка - L1 RANS
DOMAIN = dict(S_mix_max=0.6,            # S по корреляции срезается до 0,6 (mixing): выше - формула не действует
              C_rel=(0.25, 4.0),        # C/C_OPT: форма f_C v1 вне четырёхкратного отклонения - экстраполяция (оценка)
              a_eff_min=0.5)            # самостратификация: меньше половины воздуха встречает топливо (закон - оценка)
TAU_EXT_LADDER = (3.0, 10.0, 30.0, 100.0, 300.0, 1000.0)   # уточнение запаса по срыву ниже 1 (v1 - только 1/3, 1/10, 1/30)


def _pick(u, vals):
    return vals[min(int(u * len(vals)), len(vals) - 1)]


_DM = {}


def design_matrix(n, seed=SEED):
    key = (int(n), int(seed))
    if key not in _DM:
        m = max(1, math.ceil(math.log2(max(n, 2))))
        eng = qmc.Sobol(len(AXES_V2), scramble=True, **CS._rng_kw(qmc.Sobol, seed))
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            _DM[key] = eng.random_base2(m)[:n]
    return _DM[key]


def fingerprint(n, seed=SEED):
    return hashlib.sha256(np.round(design_matrix(n, seed), 12).astype("<f8").tobytes()).hexdigest()[:16]


shard_indices = CS.shard_indices


def code_hash():
    """отпечаток кода v2: comb_space, comb_crn (v1, импортируются) и comb_v2"""
    h = hashlib.sha256()
    for f in ("comb_space.py", "comb_crn.py", "comb_v2.py"):
        h.update(open(os.path.join(HERE, f), "rb").read())
    return h.hexdigest()[:16]


def concept(i, n, kv="KV-3", seed=SEED):
    """вариант i выборки n (детерминированно по i, n, seed); дискретные оси - только допустимые для схемы"""
    if not 0 <= i < n:
        raise IndexError(i)
    u = dict(zip(AXES_V2, design_matrix(n, seed)[i]))
    s = _pick(u["scheme"], CAT_V2["scheme"])
    allowed = {k: CS.SCHEMES[s].get(k, CS._ALL[k]) for k in CS._ALL if k != "passive"}
    c = dict(idx=int(i), n=int(n), seed=int(seed), kv=kv, version="v2", scheme=s)
    c["chamber"] = _pick(u["chamber"], allowed["chamber"])
    c["injector"] = _pick(u["injector"], allowed["injector"])
    dz = allowed["diesel"]
    if c["injector"] == "микросмесители":
        dz = [x for x in dz if x != "центробежная"] or dz
    c["diesel"] = _pick(u["diesel"], dz)
    c["cooling"] = _pick(u["cooling"], allowed["cooling"])
    mk = allowed["make"]
    if c["injector"] == "микросмесители" or c["cooling"] == "решётка":
        mk = [x for x in mk if x != "литьё"]
    c["make"] = _pick(u["make"], mk)
    for k in ("split", "strat", "recirc"):
        c[k] = _pick(u[k], ALLOWED[k][s])
    c["passive"] = passive_of(c)
    I, S = CS.INJ[c["injector"]], CS.SCH[s]
    spl, stra, rec = c["split"], c["strat"] == "струи", c["recirc"] == "рециркуляция"
    c["n_el"] = int(round(CS._loglerp(u["n_el"], *I["n_el"])))
    c["k_je"] = int(round(CS._lerp(u["k_je"], *I["k_je"])))
    c["p_sup_u"] = float(u["p_sup"])
    c["swirl"] = CS._lerp(u["swirl"], *I["swirl"])
    c["a_pz"] = CS._lerp(u["a_pz"], *S["a_pz"])
    if s == "TVC":
        c["a_aux"] = CS._lerp(u["a_aux"], *SCH_AUX["cavity"])
    elif s == "RQL":
        c["a_aux"] = CS._lerp(u["a_aux"], *S["a_aux"])         # воздух гашения (v1)
    elif spl == "гидравлика":
        c["a_aux"] = CS._lerp(u["a_aux"], *SCH_AUX["pilot"])
    else:
        c["a_aux"] = 0.0
    if s == "TVC":
        c["k_pass"] = CS._lerp(u["k_pass"], 0.6, 1.4)
    elif spl in ("гидравлика", "каналы"):
        c["k_pass"] = CS._lerp(u["k_pass"], 0.8, 1.6)
    else:
        c["k_pass"] = 0.0
    c["k_str"] = CS._lerp(u["k_str"], 0.5, 1.5) if stra else 0.0
    c["tau_pz_ms"] = CS._loglerp(u["tau_pz"], *S["tau_pz"])
    # пилот «гидравлика» - часть смесителей того же купола: своего времени пребывания у него нет (проверка скептика
    # 15.11.2026: tau_aux до 8·tau_pz давало пилоту запас по срыву даром); отдельное время - только у полости TVC
    c["tau_aux_ms"] = CS._loglerp(u["tau_aux"], *TAU_AUX) if s == "TVC" else (c["tau_pz_ms"] if spl == "гидравлика" else 0.0)
    c["tau_bo_ms"] = CS._loglerp(u["tau_bo"], 2.0, 20.0)
    c["tau_dil_ms"] = CS._lerp(u["tau_dil"], 1.0, 8.0)
    c["LD"] = CS._lerp(u["LD"], 1.5, 4.0)
    c["dp_inj"] = CS._lerp(u["dp_inj"], *S["dp_inj"])
    c["dp_wall"] = CS._lerp(u["dp_wall"], 0.003, 0.02) if c["cooling"] == "решётка" else 0.0
    c["porosity"] = CS._lerp(u["porosity"], *CS.POROSITY[c["cooling"]])
    c["l_mix"] = CS._lerp(u["l_mix"], *I["l_mix"])
    flow_law = spl in ("гидравлика", "каналы")
    c["x_lin"] = CS._lerp(u["x_lin"], 0.0, 0.95) if flow_law else 0.0
    c["f_por"] = CS._lerp(u["f_por"], 0.2, 0.9) if flow_law else 0.0
    c["off"] = CS._lerp(u["off"], 0.0, 1.0) if flow_law else 0.0
    c["w_grp"] = CS._lerp(u["w_grp"], 0.15, 0.5) if spl == "каналы" else 0.0
    c["l_ej"] = CS._lerp(u["l_ej"], 4.0, 15.0) if rec else 0.0
    for k, v in list(c.items()):
        if isinstance(v, float):
            c[k] = round(v, 6)
    return c


# =========================================================================== границы v2: точка 75 %
def _load_point(lf, inp, C, kw, P100, far_st, LHV_CH4, LHV_DZ, far_st_dz):
    """точка нагрузки lf тем же способом, что comb_crn.compute_boundary (постоянная частота и расход воздуха,
    PR и p3 ~ sqrt(T41), T41 - brentq по доле мощности)"""
    import cantera as ct
    PR, T41 = inp["PR"], inp["T41_C"]

    def at(T41p):
        k = math.sqrt((T41p + 273.15) / (T41 + 273.15))
        return C.cycle(PR * k, T41p, **kw), PR * k
    T41p = T41 if lf == 1.0 else brentq(lambda t: at(t)[0]["P_MW"] * 1e3 / P100 - lf, 640.0, T41, xtol=1e-4)
    r, PRp = at(T41p)
    p3, T3, W, Wf = C.P0 * PRp, r["T3_C"] + 273.15, inp["G"] * (1 - CR.X_LEAK), r["Wf"]
    air = ct.Solution("gri30.yaml")
    air.TPX = T3, p3, CR.AIR
    g4, _ = CR._mixed((Wf / W) / far_st, "CH4:1", p3, T3, C.T0); g4.equilibrate("HP"); T4eq = g4.T
    g4.TP = r["T4_C"] + 273.15, p3; g4.equilibrate("TP")
    gs, _ = CR._mixed(1.0, "CH4:1", p3, T3, C.T0); gs.equilibrate("HP")
    Wfd = Wf * LHV_CH4 / LHV_DZ
    return dict(load=lf, T41_C=T41p, PR=PRp, p3_Pa=p3, T3_K=T3, W_air=W, W_fuel=Wf, phi_ov=(Wf / W) / far_st,
                P_shaft_kW=r["P_MW"] * 1e3, T4_K=r["T4_C"] + 273.15, T4_adiab_K=T4eq, rho4=g4.density,
                rho3=air.density, mu3=air.viscosity, k3=air.thermal_conductivity, cp3=air.cp_mass, T_st_K=gs.T,
                W_fuel_diesel=Wfd, phi_ov_diesel=(Wfd / W) / far_st_dz)


def _idle_point(loads, far_st, LHV_CH4, LHV_DZ, far_st_dz, T0):
    """холостой ход без нагрузки при полной частоте (0 %): T41, PR, p3, T3 и расход топлива - линейная экстраполяция по
    точкам 50 и 30 % (ОЦЕНКА: модель цикла ниже доли топлива 0,005 не считает); второй способ - квадратичная по
    100/75/50/30 %, относительное расхождение - в check_quadratic. Свойства воздуха и равновесие на выходе - Cantera"""
    import cantera as ct
    by = {L["load"]: L for L in loads}
    a, b = by[0.5], by[0.3]
    keys = ("T41_C", "PR", "p3_Pa", "T3_K", "W_fuel")
    lin = {k: b[k] + (IDLE_LOAD - 0.3) * (a[k] - b[k]) / (0.5 - 0.3) for k in keys}
    xs = np.array(LOADS, dtype=float)
    quad = {k: float(np.polyval(np.polyfit(xs, [by[x][k] for x in LOADS], 2), IDLE_LOAD)) for k in keys}
    p3, T3, W, Wf = lin["p3_Pa"], lin["T3_K"], b["W_air"], lin["W_fuel"]
    air = ct.Solution("gri30.yaml")
    air.TPX = T3, p3, CR.AIR
    g4, _ = CR._mixed((Wf / W) / far_st, "CH4:1", p3, T3, T0); g4.equilibrate("HP")
    gs, _ = CR._mixed(1.0, "CH4:1", p3, T3, T0); gs.equilibrate("HP")
    Wfd = Wf * LHV_CH4 / LHV_DZ
    return dict(load=IDLE_LOAD, T41_C=lin["T41_C"], PR=lin["PR"], p3_Pa=p3, T3_K=T3, W_air=W, W_fuel=Wf,
                phi_ov=(Wf / W) / far_st, P_shaft_kW=0.0, T4_K=g4.T, T4_adiab_K=g4.T, rho4=g4.density,
                rho3=air.density, mu3=air.viscosity, k3=air.thermal_conductivity, cp3=air.cp_mass, T_st_K=gs.T,
                W_fuel_diesel=Wfd, phi_ov_diesel=(Wfd / W) / far_st_dz, fuel_share_of_100=Wf / by[1.0]["W_fuel"],
                method="линейная экстраполяция по точкам 50 и 30 % (оценка); T4 - адиабатическое равновесие",
                check_quadratic={k: dict(linear=lin[k], quadratic=quad[k], rel=abs(quad[k] - lin[k]) / abs(lin[k]))
                                 for k in keys})


def _bc2_key(bc1_key):
    src = bc1_key + inspect.getsource(_load_point) + inspect.getsource(compute_boundary_v2) + json.dumps(LOADS) + \
        inspect.getsource(_idle_point) + repr(IDLE_LOAD)
    return hashlib.sha256(src.encode()).hexdigest()[:16]


def compute_boundary_v2():
    """границы v1 (точки 100/50/30 %, таблицы S_L и задержек) + точка 75 %; самопроверка - точка 50 % этим же кодом"""
    t0 = time.time()
    bc1 = CR.boundary(write=False)
    C = CR._cycle()
    inp = CR._bc_inputs()
    kw = dict(eta_pc=inp["eta_pc"], eta_pt_gg=inp["eta_t"], eta_pt_pt=inp["eta_t"], x_ngv=0.0, x_rot=0.0,
              x_leak=CR.X_LEAK, dp_duct=0.0, W1=inp["G"], dp_comb=CR.DP_COMB, p_exh_ratio=CR.P_EXH)
    P100 = C.cycle(inp["PR"], inp["T41_C"], **kw)["P_MW"] * 1e3
    args = (inp, C, kw, P100, bc1["far_st_CH4"], bc1["LHV_J_kg"]["CH4"], bc1["LHV_diesel_J_kg"], bc1["far_st_diesel"])
    L75 = _load_point(0.75, *args)
    L50 = _load_point(0.5, *args)
    old = {L["load"]: L for L in bc1["loads"]}
    check50 = {k: dict(v1=old[0.5][k], v2=L50[k], rel=abs(L50[k] - old[0.5][k]) / max(abs(old[0.5][k]), 1e-30))
               for k in ("T41_C", "p3_Pa", "T3_K", "W_fuel", "rho3", "T4_K", "W_fuel_diesel")}
    loads = []
    for lf in LOADS:
        src = L75 if lf == 0.75 else old[lf]
        loads.append({k: v for k, v in src.items()})
    res = {k: v for k, v in bc1.items() if k not in ("eq_T", "loads", "key", "runtime_s", "title", "date", "check")}
    idle = _idle_point(loads, bc1["far_st_CH4"], bc1["LHV_J_kg"]["CH4"], bc1["LHV_diesel_J_kg"], bc1["far_st_diesel"],
                       C.T0)
    res.update(title="Камера сгорания KV-3, уровень 0 v2: границы (точки 100/75/50/30 % и холостой ход 0 % - оценка) и "
                     "таблицы v1",
               date=DATE, bc1_key=bc1["key"], loads=loads, idle=idle,
               check=dict(v1=bc1["check"], point_50_recomputed=check50,
                          note="точка 50 % пересчитана кодом v2 и сверена с v1 (rel - относительное расхождение)"),
               status="расчёт; точки частичной нагрузки - оценка (допущения targets_basis, как v1)")
    res["key"] = _bc2_key(bc1["key"])
    res["runtime_s"] = round(time.time() - t0, 1)
    return res


_BC = None


def boundary(path=BC_FILE, write=True):
    global _BC
    if _BC is not None:
        return _BC
    bc1 = CR.boundary(write=False)
    key = _bc2_key(bc1["key"])
    bc = None
    if os.path.exists(path):
        try:
            b = json.load(open(path, encoding="utf-8"))
            if b.get("key") == key:
                bc = b
        except (OSError, ValueError):
            bc = None
    if bc is None:
        bc = CR._rnd(compute_boundary_v2())
        if write:
            json.dump(bc, open(path, "w", encoding="utf-8"), ensure_ascii=False, indent=1)
    rows = np.array(bc["eq_table"]["rows"])
    bc["eq_T"] = lambda phi: (float(np.interp(phi, rows[:, 0], rows[:, 1])), float(np.interp(phi, rows[:, 0], rows[:, 2])))
    _BC = bc
    return bc


def li_of(load):
    return LOADS.index(load)


IDLE_LI = len(LOADS)        # индекс холостого хода в границах с холостым ходом (with_idle)


def with_idle(bc):
    """границы, у которых loads = 100/75/50/30 % + холостой ход (индекс IDLE_LI): для гидравлики и сети на 0 %"""
    return dict(bc, loads=list(bc["loads"]) + [bc["idle"]])


# =========================================================================== гидравлика топлива (газ)
R_F, G_F, T_F = CS.GAS_R_CH4, CS.GAS_G_CH4, CS.T_FUEL


def orifice_m(A, p_up, p_down):
    """расход газа через сопла площадью A (изоэнтропическое истечение с CD_FUEL, как v1)"""
    if p_up <= p_down or A <= 0:
        return 0.0
    return A * CS._orifice_flux(p_up, p_down)


def jet_state(p_up, p_down):
    """скорость и плотность струи при полном расширении до p_down (как comb_space._orifice_p0)"""
    r = max(p_down / p_up, 1e-9)
    u = math.sqrt(max(0.0, 2 * G_F / (G_F - 1) * R_F * T_F * (1 - r ** ((G_F - 1) / G_F))))
    Te = T_F * r ** ((G_F - 1) / G_F)
    return u, p_down / (R_F * Te)


def porous_p2(p1, m, a, b):
    """давление за пористой вставкой: (p1² - p2²)/(2RT) = a·m + b·m² (Дарси - Форхгеймер, изотермический газ)"""
    q = p1 * p1 - 2 * R_F * T_F * (a * m + b * m * m)
    return math.sqrt(q) if q > 0 else 0.0


def circuit_m(cir, p_s, p_down):
    """расход контура при давлении коллектора p_s: [пористая вставка a, b] -> сопла A -> p_down; -> (m, p перед соплами)"""
    A, a, b = cir["A"], cir.get("a", 0.0), cir.get("b", 0.0)
    if p_s <= p_down:
        return 0.0, p_s
    if a <= 0 and b <= 0:
        return orifice_m(A, p_s, p_down), p_s
    # наибольший расход - когда давление за вставкой падает до p_down
    F = (p_s * p_s - p_down * p_down) / (2 * R_F * T_F)
    m_hi = (F / a if b <= 0 else (-a + math.sqrt(a * a + 4 * b * F)) / (2 * b)) if F > 0 else 0.0
    g = lambda m: orifice_m(A, max(porous_p2(p_s, m, a, b), p_down), p_down) - m
    if m_hi <= 0 or g(m_hi) >= 0:
        return 0.0, p_down
    m = brentq(g, 0.0, m_hi, xtol=1e-12, rtol=1e-10)
    return m, porous_p2(p_s, m, a, b)


def solve_manifold(circuits, backs, Wf):
    """давление коллектора, при котором сумма расходов контуров равна Wf (один клапан задаёт только расход)"""
    p_lo = max(backs) * (1 + 1e-9)
    tot = lambda p: sum(circuit_m(cir, p, pb)[0] for cir, pb in zip(circuits, backs))
    if tot(p_lo) >= Wf:                    # контур с меньшим давлением за соплами один пропускает весь расход
        p_lo2 = min(backs) * (1 + 1e-9)
        p_s = brentq(lambda p: tot(p) - Wf, p_lo2, p_lo, xtol=1e-3)
    else:
        p_hi = p_lo * 1.2
        while tot(p_hi) < Wf:
            p_hi *= 1.5
            if p_hi > 1e9:
                raise ValueError("коллектор: расход недостижим")
        p_s = brentq(lambda p: tot(p) - Wf, p_lo, p_hi, xtol=1e-3)
    out = [circuit_m(cir, p_s, pb) for cir, pb in zip(circuits, backs)]
    return p_s, [m for m, _ in out], [pu for _, pu in out]


EJ_ZETA = 1.0          # потери возвратного канала эжектора, скоростных напоров подсоса (оценка)


def ejector_K(a, rs, l_ej, f, zeta=EJ_ZETA):
    """кратность эжектора как струйного насоса с выбросом в жаровую трубу (проверка скептика 15.11.2026: закон свободной
    струи Ricou - Spalding не знает стенки трубы и потерь возвратного канала): труба смешения постоянного сечения
    (1 + a)·A0, a - площадь подсоса / площадь струи (не больше ячейки элемента), баланс импульса, потери возвратного
    канала zeta, трение f·L/D, выброс в жаровую трубу; безразмерно (dp_inj·p3 = 1, плотность струи 1); оценка"""
    if a <= 0 or l_ej <= 0:
        return 0.0
    At = 1.0 + a
    LD = l_ej / math.sqrt(At)

    def resid(us):
        ds = (1 + zeta) * 0.5 * rs * us * us
        u0 = math.sqrt(2 * (1 + ds))
        ms = rs * us * a
        mm = u0 + ms
        rm = mm / (u0 + ms / rs)
        return -ds * At + 2 * (1 + ds) + rs * us * us * a - mm * mm / (rm * At) * (1 + f * LD / 2)
    if resid(1e-9) <= 0:
        return 0.0
    hi = 1.0
    while resid(hi) > 0:
        hi *= 2.0
    us = brentq(resid, 1e-9, hi, xtol=1e-12)
    return rs * us * a / math.sqrt(2 * (1 + (1 + zeta) * 0.5 * rs * us * us))


# =========================================================================== согласованное смешение
def mixing(J, k_je, l_mix, swirl, S_tol, beta=None):
    """неравномерность по отношению импульсов этой нагрузки: C = (2π/k_je)·J^0,5; f_C - форма v1; затухание по длине
    смешения с длиной K_MIX·max(1, C_OPT/C)^BETA (недопроникновение - медленнее, ОЦЕНКА)"""
    beta = BETA_LOWJ if beta is None else beta
    C = (2 * math.pi / max(1, k_je)) * math.sqrt(max(J, 0.0))
    fC = 1 + 0.5 * math.log(max(C, 1e-3) / CS.C_OPT) ** 2
    K = CS.K_MIX * max(1.0, CS.C_OPT / max(C, 1e-3)) ** beta
    S = CS.S_FLOOR + CS.S0 * fC * math.exp(-l_mix / K) * (1 - 0.3 * max(0.0, swirl - 0.4))
    return dict(C=C, f_C=fC, K_decay=K, S_mix=S, S=min(math.hypot(S, S_tol), 0.6))


def strat_frac(C, k_str, l_mix, K_decay):
    """самостратификация: доля воздуха канала, встречающая топливо: у входа (C/C_OPT)^k_str (не больше 1), затухание
    расслоения по длине смешения с той же длиной K_decay, что у неравномерности (ОЦЕНКА)"""
    a0 = min(1.0, (max(C, 1e-6) / CS.C_OPT) ** k_str)
    return 1.0 - (1.0 - a0) * math.exp(-l_mix / K_decay), a0


# =========================================================================== геометрия, контуры, нагрузки
def _d_of_A(A, n):
    return 1e3 * math.sqrt(4 * A / (math.pi * max(1, n)))


def _window(dfun, d_lo, d_hi, p_lo, p_hi):
    """интервал давления коллектора, где d_lo <= dfun(p) <= d_hi (dfun убывает с p), внутри [p_lo, p_hi]"""
    a, b = p_lo, p_hi
    if dfun(b) > d_hi or dfun(a) < d_lo:
        return None
    if dfun(a) > d_hi:
        a = brentq(lambda p: dfun(p) - d_hi, a, b, xtol=1.0)
    if dfun(b) < d_lo:
        b = brentq(lambda p: dfun(p) - d_lo, a, b, xtol=1.0)
    return (a, b) if a <= b else None


def zone_masses_v2(c, bc, cir):
    """зоны для объёма (точка 100 %, равновесные плотности v1): пилот/полость - со своим временем tau_aux"""
    L0 = bc["loads"][0]
    W, Wf, far = L0["W_air"], L0["W_fuel"], bc["far_st_CH4"]
    eqT = bc["eq_T"]
    s = c["scheme"]
    if s == "RQL":
        return CS.zone_masses(c, bc)
    tp, tb, ta, td = (1e-3 * c[k] for k in ("tau_pz_ms", "tau_bo_ms", "tau_aux_ms", "tau_dil_ms"))
    Z = []
    for k in cir:
        if k["air"] <= 0 or k["n_el"] < 1:          # контур без элементов - вариант отсеет check (ЭЛЕМЕНТЫ_КОНТУРОВ)
            continue
        phi = (k["x100"] * Wf / (k["air"] * W)) / far
        T, rho = eqT(phi)
        tau = ta if k["kind"] in ("pilot", "cavity") else tp
        Z.append(dict(name=k["kind"], m=k["air"] * W + k["x100"] * Wf, tau=tau, T=T, rho=rho))
    a_b = c["a_pz"] + (c["a_aux"] if s == "TVC" else 0.0)
    Tb, rb = eqT(Wf / W / a_b / far)
    Z.append(dict(name="догорание", m=a_b * W + Wf, tau=tb, T=Tb, rho=rb))
    Z.append(dict(name="разбавление", m=W + Wf, tau=td, T=L0["T4_K"], rho=L0["rho4"]))
    return Z


def derive(c, bc):
    """v2: воздух, контуры топлива (раздача при 100 %, сопла, пористая вставка), геометрия, потери, затем для каждой
    нагрузки - давление коллектора, расходы контуров, J, неравномерность, расслоение, запасы обратного тока"""
    L0 = bc["loads"][0]
    p3, T3, W, Wf, rho3, mu3 = (L0[k] for k in ("p3_Pa", "T3_K", "W_air", "W_fuel", "rho3", "mu3"))
    far, phi_ov = bc["far_st_CH4"], L0["phi_ov"]
    s, spl, inj = c["scheme"], c["split"], c["injector"]
    I = CS.INJ[inj]
    d = dict(passive=c["passive"])
    dp = c["dp_inj"]
    u_a = math.sqrt(2 * dp * p3 / rho3)
    A_air = c["a_pz"] * W / (CS.CD_AIR * math.sqrt(2 * rho3 * dp * p3))
    n_el = c["n_el"]
    d_air = math.sqrt(4 * A_air / (math.pi * n_el))
    d.update(u_a=u_a, A_air_m2=A_air, d_air_mm=1e3 * d_air, tau_mix_ms=1e3 * c["l_mix"] * d_air / u_a)
    # --- контуры: kind, доля воздуха (от W), доля топлива при 100 %, число струй, струй на элемент, давление за соплами
    cir = []
    if s == "RQL":
        cir.append(dict(kind="main", air=c["a_pz"], x100=1.0, n_el=n_el, k_je=c["k_je"], back="lo", porous=False))
    elif s == "TVC":
        x = min(0.6, c["k_pass"] * c["a_aux"] / phi_ov)
        hyd = spl == "гидравлика"
        cir.append(dict(kind="cavity", air=c["a_aux"], x100=x, n_el=n_el, k_je=c["k_je"], back="lo", porous=False))
        cir.append(dict(kind="main", air=c["a_pz"], x100=1 - x, n_el=n_el, k_je=c["k_je"], back="hi" if hyd else "lo",
                        porous=hyd))
    elif spl == "гидравлика":
        n_p = max(1, int(round(c["a_aux"] * n_el)))
        a = n_p / n_el
        x = a * c["k_pass"] / (1 - a + a * c["k_pass"]) if n_el - n_p >= 1 else 0.0
        cir.append(dict(kind="pilot", air=a * c["a_pz"], x100=x, n_el=n_p, k_je=c["k_je"], back="lo", porous=False))
        cir.append(dict(kind="main", air=(1 - a) * c["a_pz"], x100=1 - x, n_el=n_el - n_p, k_je=c["k_je"], back="hi",
                        porous=True))
    elif spl == "каналы":
        n_a = max(1, int(round(c["w_grp"] * n_el)))
        w = n_a / n_el
        x = w * c["k_pass"] / (1 - w + w * c["k_pass"]) if n_el - n_a >= 1 else 0.0
        cir.append(dict(kind="grpA", air=w * c["a_pz"], x100=x, n_el=n_a, k_je=c["k_je"], back="lo", porous=False))
        cir.append(dict(kind="grpB", air=(1 - w) * c["a_pz"], x100=1 - x, n_el=n_el - n_a, k_je=c["k_je"], back="hi",
                        porous=True))
    else:
        cir.append(dict(kind="main", air=c["a_pz"], x100=1.0, n_el=n_el, k_je=c["k_je"], back="lo", porous=False))
    d["n_el_bad"] = any(k["n_el"] < 1 for k in cir)
    for k in cir:
        k["n_jet"] = max(1, k["n_el"]) * k["k_je"]
        k["m100"] = k["x100"] * Wf
    d["fuel_split100"] = next((k["x100"] for k in cir if k["kind"] in ("pilot", "cavity", "grpA")), 0.0)
    # --- окно давления коллектора при 100 %
    p_lo100 = p3 * (1 - dp)
    p_hi100 = p3 * (1 - dp + c["off"] * dp)
    back100 = {"lo": p_lo100, "hi": p_hi100}
    d_min = max(I["d_jet"][0], CS.D_JET_MIN[c["make"]])
    d_max = I["d_jet"][1]

    def size(k, p_s):
        """площадь сопел контура при давлении коллектора p_s (у пористого - давление перед соплами по f_por)"""
        pb = back100[k["back"]]
        p_j = pb + (1 - c["f_por"]) * (p_s - pb) if k["porous"] else p_s
        if k["m100"] <= 0:
            return 0.0, p_j
        return k["m100"] / max(CS._orifice_flux(p_j, pb), 1e-30), p_j
    win = (CS.P_FUEL_WIN[0] * p3, CS.P_FUEL_WIN[1] * p3)
    ok_win = True
    for k in cir:
        if k["m100"] <= 0:
            continue
        w2 = _window(lambda p, k=k: _d_of_A(size(k, p)[0], k["n_jet"]), d_min, d_max, *win)
        if w2 is None:
            ok_win = False
            break
        win = (max(win[0], w2[0]), min(win[1], w2[1]))
        if win[0] > win[1]:
            ok_win = False
            break
    d["fuel_window_ok"] = ok_win
    d["p_sup_window_ratio"] = [win[0] / p3, win[1] / p3] if ok_win else None
    p_s100 = CS._loglerp(c["p_sup_u"], *win) if ok_win else CS.P_FUEL_WIN[1] * p3
    for k in cir:
        A, p_j = size(k, p_s100)
        k["A"] = A
        k["d_jet_mm"] = _d_of_A(A, k["n_jet"]) if A > 0 else 0.0
        if k["porous"] and k["m100"] > 0:
            F = (p_s100 ** 2 - p_j ** 2) / (2 * R_F * T_F)
            k["a"] = c["x_lin"] * F / k["m100"]
            k["b"] = (1 - c["x_lin"]) * F / k["m100"] ** 2
        k["S_tol"] = math.hypot(2 * CS.D_SIG_FUEL[c["make"]] / max(k["d_jet_mm"], 1e-3) / math.sqrt(max(1, k["k_je"])),
                                2 * CS.D_SIG_AIR[c["make"]] / (1e3 * d_air))
    main = next(k for k in cir if k["kind"] in ("main", "grpB"))
    d.update(circuits=[{kk: vv for kk, vv in k.items()} for k in cir], p_sup_Pa=p_s100, p_fuel_Pa=p_s100,
             p_fuel_ratio=p_s100 / p3, d_jet_mm=main["d_jet_mm"], n_jet=main["n_jet"],
             d_jet_window_mm=[d_min, d_max] if ok_win else [1.0, 0.0])
    # --- трение смесительного канала и трубы эжектора (Haaland, шероховатость по способу изготовления)
    Re_ch = CS.CD_AIR * u_a * d_air * rho3 / mu3
    eps = 1e-3 * CS.KS_CH[c["make"]] / d_air
    f_ch = (-1.8 * math.log10((eps / 3.7) ** 1.11 + 6.9 / max(Re_ch, 2300.0))) ** -2
    fr_mix = f_ch * c["l_mix"] * CS.CD_AIR ** 2 * dp
    fr_ej = f_ch * c["l_ej"] * CS.CD_AIR ** 2 * dp          # оценка: скорость в трубе эжектора ~ скорость канала
    d.update(Re_ch=Re_ch, f_ch=f_ch, dp_mix_fric=fr_mix + fr_ej, dp_ej_fric=fr_ej)
    # --- объёмы зон, геометрия, охлаждение, разбавление
    zones = zone_masses_v2(c, bc, cir)
    V = sum(z["m"] * z["tau"] / z["rho"] for z in zones)
    geo = CS.geometry(c, bc, V)
    geo["L_mod"] += (c["l_mix"] + c["l_ej"]) * d_air
    dpL = dp * p3
    m_cool = CS.CD_AIR * c["porosity"] * geo["A_wall"] * math.sqrt(2 * rho3 * dpL)
    a_cool = m_cool / W
    a_aux_out = c["a_aux"] if s in ("TVC", "RQL") else 0.0
    a_dil = 1.0 - c["a_pz"] - a_aux_out - a_cool
    d.update(zones=[dict(name=z["name"], m=z["m"], tau_ms=1e3 * z["tau"], T_K=z["T"], V_L=1e3 * z["m"] * z["tau"] / z["rho"])
                    for z in zones], V_L=1e3 * V, a_cool=a_cool, a_aux_out=a_aux_out, a_dil=a_dil, **geo)
    q3 = 0.5 * 1.4 * CS.M3 ** 2 * (1 + 0.2 * CS.M3 ** 2) ** -3.5
    d["dp_diff"] = CS.K_DIFF[c["chamber"]] * q3
    d["dp_duct"] = CS.DP_DUCT[c["chamber"]]
    d["dp_series_wall"] = c["dp_wall"] if c["cooling"] == "решётка" else 0.0
    q_ref = 0.5 * rho3 * (W / CS.N_CAN[c["chamber"]] / (rho3 * geo["A_ref"])) ** 2
    d["dpL_qref"] = dpL / q_ref
    # --- нагрузки: давление коллектора, расходы контуров, смешение
    d["loads"] = [load_state(c, d, bc, li, cir, A_air) for li in range(len(LOADS))] if ok_win and not d["n_el_bad"] \
        else None
    # холостой ход (0 %, границы - оценка): только для обратного тока и устойчивости опорной зоны (цель load_range)
    d["idle"] = load_state(c, d, with_idle(bc), IDLE_LI, cir, A_air) if d["loads"] else None
    if d["loads"]:
        l100 = d["loads"][0]
        mi = l100["circ"][main["kind"]]
        d.update(J=mi["J"], holdeman_C=mi["C"], f_C=mi["f_C"], S_mix=mi["S_mix"], S_tol=main["S_tol"], S_unmix=mi["S"],
                 u_fuel=mi["u_f"])
        prem = [(v["phi_reached"], v["S"], kk) for kk, v in l100["circ"].items() if kk != "cavity"]
        ph, Sx, kk = max(prem)
        d["phi_pz_100"], d["S_unmix_fb"], d["fb_circuit_100"] = ph, Sx, kk
    else:
        d.update(S_unmix=0.6, S_mix=0.6, S_tol=0.0, J=0.0, holdeman_C=0.0, f_C=0.0, u_fuel=0.0,
                 phi_pz_100=phi_ov / c["a_pz"], S_unmix_fb=0.6)
    return d


def load_state(c, d, bc, li, cir, A_air):
    """одна нагрузка: скорость воздуха при неизменной геометрии, давления за соплами, коллектор, J и смешение контуров"""
    L = bc["loads"][li]
    p3, rho3, W, Wf = L["p3_Pa"], L["rho3"], L["W_air"], L["W_fuel"]
    u_a = c["a_pz"] * W / (CS.CD_AIR * A_air * rho3)        # воздух делится геометрией - доли постоянны
    dp = rho3 * u_a ** 2 / (2 * p3)
    back = {"lo": p3 * (1 - dp), "hi": p3 * (1 - dp + c["off"] * dp)}
    act = [k for k in cir if k["A"] > 0]
    p_s, ms, pus = solve_manifold(act, [back[k["back"]] for k in act], Wf)
    st = dict(load=L["load"], dp_air=dp, u_a=u_a, p_sup_Pa=p_s, p_sup_ratio=p_s / p3,
              backflow_margin=(p_s - max(back[k["back"]] for k in act)) / p3, circ={})
    dpf_w, inv_eff = 0.0, 0.0
    for k, m, pu in zip(act, ms, pus):
        pb = back[k["back"]]
        u_f, rho_f = jet_state(pu, pb) if m > 0 else (0.0, 0.0)
        J = rho_f * u_f ** 2 / (rho3 * u_a ** 2)
        mx = mixing(J, k["k_je"], c["l_mix"], c["swirl"] if k["kind"] != "cavity" else 0.0, k["S_tol"])
        a_eff, a0 = 1.0, 1.0
        if c["strat"] == "струи" and k["kind"] in ("main", "grpB"):
            a_eff, a0 = strat_frac(mx["C"], c["k_str"], c["l_mix"], mx["K_decay"])
        phi = (m / (k["air"] * W)) / bc["far_st_CH4"]
        dpf = (p_s - pb) / p3
        dpf_w += dpf * m / Wf
        e = 1e-4 * max(p_s - pb, 1.0)                     # закон контура: n = dln m / dln Δp (сопло 0,5, запирание 0, Дарси до 1)
        m_a, m_b = circuit_m(k, p_s, pb - e)[0], circuit_m(k, p_s, pb + e)[0]
        n_law = ((math.log(m_a) - math.log(m_b)) / (math.log(p_s - pb + e) - math.log(p_s - pb - e))
                 if m_a > 0 and m_b > 0 else 0.5)
        inv_eff += (m / Wf) * n_law / max(dpf, 1e-6)
        st["circ"][k["kind"]] = dict(m=m, share=m / Wf, p_up=pu, p_back=pb, dpf=dpf, u_f=u_f, J=J, a_eff=a_eff, a0=a0,
                                     phi=phi, phi_reached=phi / a_eff, **mx)
    st["dpf_weighted"] = dpf_w
    st["dpf_eff"] = 0.5 / max(inv_eff, 1e-12)            # перепад сопла с тем же откликом расхода на пульсацию давления
    return st


def check(c, d, bc):
    """правила допустимости v1 (comb_space.check) + окно коллектора, число элементов контуров, обратный ток"""
    bad = list(CS.check(c, d, bc))
    if not d["fuel_window_ok"] and not any(b[0] == "ТОПЛИВО" for b in bad):
        bad.append(("ТОПЛИВО", None, list(CS.P_FUEL_WIN), "нет давления коллектора, при котором диаметры струй всех "
                    "контуров в допустимом диапазоне"))
    if d["n_el_bad"]:
        bad.append(("ЭЛЕМЕНТЫ_КОНТУРОВ", c["n_el"], 2, "на контур (пилот/основной, группы A/B) не хватает элементов"))
    if d["loads"]:
        # нагрузки 100/75/50/30 % и холостой ход (0 %, оценка): двигатель проходит холостой ход при каждом пуске, а
        # пассивное деление там слабее всего (проверка скептика 15.11.2026)
        sts = list(d["loads"]) + ([d["idle"]] if d.get("idle") else [])
        bm = min(x["backflow_margin"] for x in sts)
        bcx = with_idle(bc)
        if bm < DP_BACK_MIN:
            where = min(sts, key=lambda x: x["backflow_margin"])["load"]
            bad.append(("ОБРАТНЫЙ_ТОК", round(bm, 5), DP_BACK_MIN, "давление коллектора над давлением за соплами "
                        "основного контура меньше 1 % p3 на одной из нагрузок - воздух может войти в топливный контур "
                        f"(наименьший запас - {int(round(where * 100))} %)"))
        elif any(k.get("porous") and k.get("A", 0) > 0 for k in d["circuits"]):
            soft = [dict(k, a=k.get("a", 0.0) * POR_TOL, b=k.get("b", 0.0) * POR_TOL) if k.get("porous") else k
                    for k in d["circuits"]]
            bt = min(load_state(c, d, bcx, li, soft, d["A_air_m2"])["backflow_margin"] for li in range(len(sts)))
            if bt < DP_BACK_MIN:
                bad.append(("ОБРАТНЫЙ_ТОК_РАЗБРОС", round(bt, 5), DP_BACK_MIN, "при сопротивлении пористой вставки ниже "
                            "расчётного на 20 % (разброс печати, оценка) запас коллектора меньше 1 % p3"))
    return bad


def domain(d):
    """область облика (DOMAIN, оценка): по всем контурам с топливом на 100/75/50/30 % - S по корреляции до среза,
    C/C_OPT, доля воздуха, встречающая топливо; вне области - пометка (проверка скептика 15.11.2026: единственный
    устойчивый без RQL вариант пробы стоял вне области корреляций), отсева нет - проверка на L1 RANS"""
    if not d.get("loads"):
        return None
    S_mx, C_lo, C_hi, a_lo, S_cap = 0.0, 1e9, 0.0, 1.0, 0.0
    for st in d["loads"]:
        for kind, k in st["circ"].items():
            if k["share"] < 0.02:
                continue
            S_mx = max(S_mx, k["S_mix"])
            S_cap = max(S_cap, k["S"])
            C_lo, C_hi = min(C_lo, k["C"] / CS.C_OPT), max(C_hi, k["C"] / CS.C_OPT)
            a_lo = min(a_lo, k["a_eff"])
    flags = []
    if S_mx > DOMAIN["S_mix_max"]:
        flags.append("S_по_корреляции_выше_среза")
    if C_lo < DOMAIN["C_rel"][0] or C_hi > DOMAIN["C_rel"][1]:
        flags.append("C_вне_области")
    if a_lo < DOMAIN["a_eff_min"]:
        flags.append("расслоение_сильнее_половины")
    return dict(ok=not flags, flags=flags, S_mix_max=S_mx, S_used_max=S_cap, C_rel_min=C_lo, C_rel_max=C_hi,
                a_eff_min=a_lo, nodes_clipped=bool(S_cap > 0.98 / _GH_ZMIN), limits=DOMAIN)


# =========================================================================== сеть реакторов v2
ANCHOR = ("pilot", "cavity", "grpA")
_GH_ZMIN = float(abs(min(CR.GH[5][0])))


def _node_streams_v2(K, m_air, m_fuel, S, n_nodes, T3, Tf, p, fuel="CH4:1"):
    """узлы неравномерности как comb_crn._node_streams (Гаусс - Эрмит, phi_k ~ 1 + S·z_k, не меньше 0,02), но топливо
    сохраняется: после среза бедного узла доли топлива нормируются на Σ w·f (v1 при S > 0,343 добавлял топливо:
    +0,18 % при S = 0,4, +0,83 % при S = 0,6 - проверка скептика 15.11.2026)"""
    z, w = CR.GH[n_nodes]
    f = [max(0.02, 1 + S * zk) for zk in z]
    sw = sum(wk * fk for wk, fk in zip(w, f))
    out = [(wk, K.fresh(m_air * wk, m_fuel * wk * fk / sw, T3, Tf, p, fuel)) for wk, fk in zip(w, f)]
    return out, sum(x.m for _, x in out)


def tau_ext_v2(K, inlet, tau_act, Kr=0.0, init=None):
    """запас по срыву как v1 (comb_crn.Kin.tau_ext); ниже 1 - уточнение (проверка скептика 15.11.2026: v1 давал только
    1/3, 1/10, 1/30): лестница TAU_EXT_LADDER до первого горящего времени и 7 шагов бисекции в логарифме между
    погасшим и горящим (продолжение по горящему решению), шаг <= 1 %; 0 - не горит и при 1000·tau"""
    m = K.tau_ext(inlet, tau_act, Kr, init)
    if m >= 1.0:
        return m
    tau_r = tau_act / (1 + Kr)
    lo, hi, s_hi = 1.0, None, None
    for f in TAU_EXT_LADDER:
        _, lit, s_ = K.psr(inlet, tau_r * f)
        if lit:
            hi, s_hi = f, s_
            break
        lo = f
    if hi is None:
        return 0.0
    for _ in range(7):
        mid = math.sqrt(lo * hi)
        _, lit, s_ = K.psr(inlet, tau_r * mid, s_hi)
        if lit:
            hi, s_hi = mid, s_
        else:
            lo = mid
    return 1.0 / hi


def zone_taus(c, d, bc, li, st):
    """времена пребывания зон на нагрузке li при НЕИЗМЕННОМ объёме зон (геометрия та же, что при 100 %; проверка скептика
    15.11.2026: постоянное время пребывания завышало время богатой зоны RQL на 30 % и занижало у бедных зон):
    tau(L) = tau_100·[rho(L)/m(L)]/[rho(100)/m(100)], m - расход зоны этой нагрузки, rho - равновесие при phi зоны этой
    нагрузки по таблице 100 % × p3(L)/p3(100) (оценка: сдвиг T3 не учтён), у разбавления - rho4 точки. При 100 % -
    ровно времена осей. -> секунды: виды контуров (pilot, cavity, grpA, main, grpB) или pz/q (RQL), bo, dil"""
    eqT, far = bc["eq_T"], bc["far_st_CH4"]
    L0, L = bc["loads"][0], bc["loads"][li]
    st0 = d["loads"][0]
    dens = lambda Ld, phi: eqT(phi)[1] * Ld["p3_Pa"] / L0["p3_Pa"]
    ratio = lambda fn: fn(L, st) / fn(L0, st0)
    out = {}
    if c["scheme"] == "RQL":
        a, aq = c["a_pz"], c["a_pz"] + c["a_aux"]
        rich = lambda Ld, S_: dens(Ld, Ld["W_fuel"] / Ld["W_air"] / a / far) / (a * Ld["W_air"] + Ld["W_fuel"])
        lean = lambda Ld, S_: dens(Ld, Ld["W_fuel"] / Ld["W_air"] / aq / far) / (aq * Ld["W_air"] + Ld["W_fuel"])
        tq = 1e-3 * next(z["tau_ms"] for z in d["zones"] if z["name"] == "гашение")
        out.update(pz=1e-3 * c["tau_pz_ms"] * ratio(rich), q=tq * ratio(lean), bo=1e-3 * c["tau_bo_ms"] * ratio(lean))
    else:
        for k in d["circuits"]:
            kind = k["kind"]
            if kind not in st["circ"] or kind not in st0["circ"]:
                continue
            fn = lambda Ld, S_, k=k: dens(Ld, S_["circ"][k["kind"]]["phi"]) / (k["air"] * Ld["W_air"] +
                                                                                S_["circ"][k["kind"]]["m"])
            tau0 = c["tau_aux_ms"] if kind in ("pilot", "cavity") else c["tau_pz_ms"]
            out[kind] = 1e-3 * tau0 * ratio(fn)
        a_b = c["a_pz"] + (c["a_aux"] if c["scheme"] == "TVC" else 0.0)
        bo = lambda Ld, S_: dens(Ld, Ld["W_fuel"] / Ld["W_air"] / a_b / far) / (a_b * Ld["W_air"] + Ld["W_fuel"])
        out["bo"] = 1e-3 * c["tau_bo_ms"] * ratio(bo)
    dil = lambda Ld, S_: Ld["rho4"] / (Ld["W_air"] + Ld["W_fuel"])
    out["dil"] = 1e-3 * c["tau_dil_ms"] * ratio(dil)
    return out


def _c_at_load(c, d, T):
    """понятие и производные с временами зон этой нагрузки (для сетей v1: RQL на газе, дизель)"""
    kind = "pz" if c["scheme"] == "RQL" else ("grpB" if c["split"] == "каналы" else "main")
    c2 = dict(c, tau_pz_ms=1e3 * T.get(kind, 1e-3 * c["tau_pz_ms"]), tau_bo_ms=1e3 * T["bo"], tau_dil_ms=1e3 * T["dil"])
    d2 = dict(d)
    if "q" in T:
        d2["zones"] = [dict(z, tau_ms=1e3 * T["q"]) if z["name"] == "гашение" else z for z in d["zones"]]
    return c2, d2


def gas_network(c, d, bc, li, st=None, cross=False):
    """сеть на газе при нагрузке bc["loads"][li] (st - состояние гидравлики этой нагрузки, по умолчанию d["loads"][li]):
    опорная зона (пилот, полость, группа A) -> её горячий газ в узлы основной зоны (пилот, полость; группа A - только
    при cross=True, поджигание группы B газом группы A - верхняя оценка) -> основная зона с расслоением и
    рециркуляцией -> совместное догорание (с воздухом, не встретившим топливо) -> разбавление. Узлы неравномерности -
    5 узлов Гаусса - Эрмита у всех контуров (топливо сохраняется), S каждого контура - по J этой нагрузки; времена зон -
    при неизменном объёме (zone_taus); запасы по срыву ниже 1 - уточнённые (tau_ext_v2)."""
    s = c["scheme"]
    st = d["loads"][li] if st is None else st
    T = zone_taus(c, d, bc, li, st)
    if s == "RQL":
        c2, d2 = _c_at_load(c, d, T)
        r = CR.gas_network(c2, dict(d2, S_unmix=st["circ"]["main"]["S"]), bc, li)
        if (r["margin_lbo"] or 0.0) < 1.0:                   # уточнение запаса богатой зоны ниже 1 (v1 - 1/3, 1/10, 1/30)
            K = CR.kin()
            L = bc["loads"][li]
            K_crz = CR.K_CRZ * max(0.0, c["swirl"] - 0.4) if c["injector"] == "закрутка" else 0.0
            fr = K.fresh(c["a_pz"] * L["W_air"], L["W_fuel"], L["T3_K"], 288.15, L["p3_Pa"])
            _, _, rst, rin = K.psr_recirc(fr, T["pz"], K_crz)
            m = tau_ext_v2(K, rin, T["pz"], K_crz, rst)
            r["margins"] = dict(r["margins"], богатая=m)
            r["margin_lbo"] = min(r["margins"].values())
        r.update(margin_all=r["margin_lbo"], anchor="богатая", p_sup_ratio=st["p_sup_ratio"], dpf=st["dpf_eff"],
                 backflow_margin=st["backflow_margin"], tau_ms={k: 1e3 * v for k, v in T.items()})
        return r
    K = CR.kin()
    L = bc["loads"][li]
    p, T3, W, Wf = L["p3_Pa"], L["T3_K"], L["W_air"], L["W_fuel"]
    Tf = 288.15
    tb, td = T["bo"], T["dil"]
    nn = CR.N_NODES[s]
    res = dict(load=L["load"], p_sup_ratio=st["p_sup_ratio"], dpf=st["dpf_eff"], backflow_margin=st["backflow_margin"],
               tau_ms={k: 1e3 * v for k, v in T.items()}, cross=bool(cross))
    seg, Ts, margins = [], [], []
    K_crz = CR.K_CRZ * max(0.0, c["swirl"] - 0.4) if c["injector"] == "закрутка" else 0.0
    Krec, S_div = K_crz, 1.0
    main_kind = "grpB" if c["split"] == "каналы" else "main"
    cm = st["circ"][main_kind]
    if s == "FLOX":
        fx = CR.flox_recirc(c, d, bc, li)
        Krec, S_div = fx["K"], 1.0 + fx["K"]               # как v1: S/(1+K) только у внутренней рециркуляции FLOX
        res["K_flox_int"] = fx["K"]
    if c["recirc"] == "рециркуляция":
        T_b = max(bc["eq_T"](cm["phi_reached"])[0], 800.0)
        K_rs = 0.32 * c["l_ej"] * math.sqrt(T3 / T_b) - 1.0
        D_cell = math.sqrt(4 * d["A_L"] * CS.N_CAN[c["chamber"]] / (math.pi * max(1, c["n_el"])))
        a_max = max(0.0, (D_cell / (1e-3 * d["d_air_mm"])) ** 2 - 1.0)
        K_jp = ejector_K(a_max, T3 / T_b, c["l_ej"], d["f_ch"])
        K_ej = min(CR.KV_MAX, max(0.0, min(K_rs, K_jp)))
        res["K_ej"], res["K_ej_rs"], res["K_ej_jetpump"], res["ej_area_ratio"] = K_ej, K_rs, K_jp, a_max
        Krec = min(CR.KV_MAX, Krec + K_ej)
    res["K_recirc"] = Krec

    def chain(m_air, m_fuel, S_eff, Kr, tau, pre=None):
        nodes, _ = _node_streams_v2(K, m_air, m_fuel, S_eff, nn, T3, Tf, p)
        lit_w, outs = 0.0, []
        for wk, fr in nodes:
            inlet = fr if pre is None else K.mix(fr, CR.Stream(pre.m * wk, pre.h, pre.Y, pre.p))
            o, lit, _, _ = K.psr_recirc(inlet, tau, Kr)
            outs.append((wk, o))
            lit_w += wk if lit else 0.0
        return lit_w, outs

    pre = None
    anchor = None
    for kind in ANCHOR:
        if kind not in st["circ"]:
            continue
        k = st["circ"][kind]
        air = next(x["air"] for x in d["circuits"] if x["kind"] == kind) * W
        tau = T[kind]
        Kr = 0.0 if kind == "cavity" else (K_crz if kind == "pilot" else Krec)
        lit_w, outs = chain(air, k["m"], k["S"] / (S_div if kind == "grpA" else 1.0), Kr, tau)
        fr = K.fresh(air, k["m"], T3, Tf, p)
        _, plit, pst, pin = K.psr_recirc(fr, tau, Kr)
        margins.append((kind, tau_ext_v2(K, pin, tau, Kr, pst)))
        res[f"phi_{kind}"] = k["phi"]
        res[f"fuel_share_{kind}"] = k["share"]
        res[f"lit_{kind}"] = lit_w
        anchor = kind
        if (kind in ("pilot", "cavity") or cross) and lit_w >= 0.5:
            pre = K.mix(*[o for _, o in outs])               # горячий газ пилота/полости (группы A при cross) - в основную
        else:
            for wk, o in outs:
                o2 = K.pfr(o, CR.F_SEG * tb)
                seg.append(o2); Ts.append((wk * k["share"], K.T(o2)))
    # основная зона
    tp = T[main_kind]
    air_m = next(x["air"] for x in d["circuits"] if x["kind"] == main_kind) * W
    a_eff = cm["a_eff"]
    bypass = (1 - a_eff) * air_m
    air_r = a_eff * air_m
    res["a_eff_frac"], res["a0_frac"] = a_eff, cm["a0"]
    S_main = cm["S"] / S_div
    res["S_main"], res["J_main"], res["C_main"] = S_main, cm["J"], cm["C"]
    lit_m, outs = chain(air_r, cm["m"], S_main, Krec, tp, pre)
    for wk, o in outs:
        o2 = K.pfr(o, CR.F_SEG * tb)
        seg.append(o2); Ts.append((wk * cm["share"], K.T(o2)))
    fr = K.fresh(air_r, cm["m"], T3, Tf, p)
    if pre is not None:
        fr = K.mix(fr, pre)
    _, _, st0, inl0 = K.psr_recirc(fr, tp, Krec)
    margins.append((main_kind, tau_ext_v2(K, inl0, tp, Krec, st0)))
    res["phi_main"] = cm["phi"]
    res["phi_main_reached"] = cm["phi_reached"]
    res["fuel_share_main"] = cm["share"]
    res["lit_main"] = lit_m
    # совместное догорание: струйки + воздух, не встретивший топливо (расслоение) -> разбавление
    parts = list(seg)
    if bypass > 0:
        parts.append(K.air(bypass, T3, p))
    hot = K.pfr(K.mix(*parts), (1 - CR.F_SEG) * tb)
    T_hot = K.T(hot)
    wsum = max(sum(w for w, _ in Ts), 1e-12)
    T_mean = sum(w * T_ for w, T_ in Ts) / wsum
    sigT = math.sqrt(sum(w * (T_ - T_mean) ** 2 for w, T_ in Ts) / wsum)
    dil = K.air(d["a_dil"] * W + d["a_cool"] * W, T3, p)
    ex = K.pfr(K.mix(hot, dil), td)
    em = K.emissions(ex, Wf / ex.m, bc["LHV_J_kg"], bc["LHV_J_kg"]["CH4"])
    mg = dict(margins)
    m_anchor = mg[anchor] if anchor else mg[main_kind]
    res.update(NOx=em["NOx"], CO=em["CO"], UHC=em["UHC"], eta_c=em["eta_c"], O2_dry=em["O2_dry"], T_exit=em["T"],
               T_hot=T_hot, T_flame=T_mean, sigma_T_hot=sigT, m_hot=hot.m, lit_frac=lit_m,
               anchor=anchor or main_kind, margin_lbo=m_anchor, margin_all=min(mg.values()), margins=mg)
    return res


def flashback_v2(c, d, bc):
    """проскок: форма v1 (comb_crn.flashback) по самому богатому предварительно смешанному контуру при 100 % (доля H2)
    и по каждому контуру на каждой нагрузке (метан): phi с учётом расслоения, S этой нагрузки, скорость воздуха этой
    нагрузки. Таблица S_L - при p3, T3 точки 100 % (на малой нагрузке T3 ниже - запас занижен, оценка)"""
    if c["scheme"] == "RQL":
        return CR.flashback(c, d, bc)
    fb = CR.flashback(c, dict(d, phi_pz_100=d["phi_pz_100"], S_unmix=d["S_unmix_fb"]), bc)
    L0 = bc["loads"][0]
    nu = L0["mu3"] / L0["rho3"]
    rows = sorted([r for r in bc["flame_table"] if r["H2_vol"] == 0.0], key=lambda r: r["T_ad_K"])
    w = 1 / (1 + 4 * c["swirl"] ** 2) if c["injector"] == "закрутка" else 0.6
    worst = (1e9, None)
    for st in d["loads"]:
        for kind, k in st["circ"].items():
            if kind == "cavity" or k["m"] <= 0:
                continue
            Tad = bc["eq_T"](k["phi_reached"] * (1 + 2 * k["S"]))[0]
            u = st["u_a"] * w
            Re = u * 1e-3 * d["d_air_mm"] / nu
            g_w = 0.316 * max(Re, 2300.0) ** -0.25 * u ** 2 / (8 * nu)
            SL = float(np.interp(Tad, [r["T_ad_K"] for r in rows], [r["S_L"] for r in rows]))
            al = float(np.interp(Tad, [r["T_ad_K"] for r in rows], [r["alpha_u"] for r in rows]))
            m = g_w / (CR.K_FB * SL ** 2 / al)
            if m < worst[0]:
                worst = (m, f"{kind}@{int(st['load'] * 100)}")
    fb["margin_CH4_100"] = fb["margin_CH4"]
    fb["margin_CH4"] = min(fb["margin_CH4"], worst[0])
    fb["margin_CH4_where"] = worst[1]
    return fb


# =========================================================================== цели v2
def _tg(tid):
    return next(t for t in CR.TARGETS["targets"] if t["id"] == tid)


def lowem_value(row):
    """цель «низкие выбросы 30-100 % на одном клапане»: наибольшее из отношений к пределам на 30/50/75/100 %"""
    g = {r["load"]: r for r in row["gas"]}
    nox, co = _tg("nox_gas"), _tg("co_gas")
    tv = _tg("fuel_valves")["target"]
    terms = {}
    for ld in LOADS:
        pct = int(round(ld * 100))
        lo = ld < 0.5
        terms[f"NOx_{pct}"] = g[ld]["NOx"] / (nox["also"]["max"] if lo else nox["target"]["max"])
        terms[f"CO_{pct}"] = g[ld]["CO"] / (co["also"]["max"] if lo else co["target"]["max"])
        terms[f"lbo_{pct}"] = 1.0 / max(g[ld]["margin_lbo"] or 0.0, 1e-6)
        terms[f"eta_{pct}"] = (1.0 - g[ld]["eta_c"]) / (1.0 - CR.ETA_MIN)
    m = row["module"]
    terms["valves"] = max(m["valves_gas"], m["valves_diesel"]) / tv["control_valves_per_fuel"]
    terms["actuators"] = 1.0 + m.get("actuators", 0)          # пассивные элементы v2 без движущихся частей: 0
    stable = max(v for k, v in terms.items() if k.startswith(("lbo_", "eta_")))
    return max(terms.values()), terms, stable


def compliance(c, row):
    out, hard = CR.compliance(c, row)                         # 13 целей targets.json (load_range - по margin_lbo опорной зоны)
    # load_range требует устойчивого горения 0-100 % (stable_pct [0, 100]): v1 проверял только 30-100 %, v2 добавляет
    # холостой ход (0 %, границы - экстраполяция, оценка) - запас по срыву опорной зоны (проверка скептика 15.11.2026)
    idle = row.get("idle") or {}
    m_idle = idle.get("margin_lbo")
    lr = out["load_range"]
    if lr["value"] is not None:
        v_lr = max(lr["value"], 1.0 / max(m_idle or 0.0, 1e-6))
        out["load_range"] = dict(value=v_lr, limit=1.0, ok=bool(v_lr <= 1.0), viol=max(0.0, v_lr - 1.0),
                                 value_30_100=lr["value"], idle_margin_lbo=m_idle,
                                 note="1/запас опорной зоны на холостом ходу (0 %, оценка) добавлен к v1 (30-100 %)")
    v, terms, stable = lowem_value(row)
    out[NEW_TARGET] = dict(value=v, limit=1.0, ok=bool(v <= 1.0), viol=max(0.0, v - 1.0),
                           binding=max(terms, key=terms.get), stable_30_100=bool(stable <= 1.0),
                           stable_0_100=bool(stable <= 1.0 and (m_idle or 0.0) >= 1.0))
    strict = dict(row, gas=[dict(r, margin_lbo=r.get("margin_all", r["margin_lbo"])) for r in row["gas"]])
    vs, _, sts = lowem_value(strict)
    out[NEW_TARGET]["ok_strict_all_zones"] = bool(vs <= 1.0)
    dom = row.get("domain") or {}
    out[NEW_TARGET]["in_domain"] = bool(dom.get("ok", False))
    if row.get("gas_cross"):                                  # каналы: группа B поджигается газом группы A (верхняя оценка)
        vx, tx, sx = lowem_value(dict(row, gas=row["gas_cross"]))
        out[NEW_TARGET].update(value_cross=vx, ok_cross=bool(vx <= 1.0), stable_30_100_cross=bool(sx <= 1.0),
                               binding_cross=max(tx, key=tx.get))
    return out, hard + [NEW_TARGET]


def targets_v2():
    """targets_v2.json: 13 целей targets.json как база + новая жёсткая цель (решение владельца 15.11.2026)"""
    T = CR.TARGETS
    nox, co, fv = _tg("nox_gas"), _tg("co_gas"), _tg("fuel_valves")
    new = dict(
        id=NEW_TARGET, name="Низкие выбросы 30-100 % на одном клапане, без приводов", kind="hard",
        target=dict(loads_pct=[30, 50, 75, 100], NOx_max_ppm15={"50-100": nox["target"]["max"], "30": nox["also"]["max"]},
                    CO_max_ppm15={"50-100": co["target"]["max"], "30": co["also"]["max"]}, eta_min=CR.ETA_MIN,
                    margin_lbo_min=1.0, control_valves_per_fuel=fv["target"]["control_valves_per_fuel"],
                    actuators_in_combustor=fv["target"]["actuators_in_combustor"], bleed_valves=fv["target"]["bleed_valves"]),
        metric=("наибольшее из: NOx/9 и CO/15 на 50, 75, 100 %; NOx/15 и CO/50 на 30 %; 1/запас по срыву опорной зоны "
                "(пилот, полость, группа A; без них - основная) и (1 - полнота)/0,01 на каждой из 30/50/75/100 %; "
                "клапанов на топливо / 1; 1 + число приводов. Выполнено, если не больше 1"),
        source=("решение владельца 15.11.2026 (routes/gtu/README.md, «Решения владельца», п. 2); пределы - nox_gas, co_gas "
                "(с also для 30 %), load_range, fuel_valves этого файла; полнота 0,99 - comb_crn.ETA_MIN (оценка)"),
        status="цель проекта; уровень 0 - облик (сеть реакторов), подтверждение - уровни 1-2",
        notes=["точка 75 % добавлена к 30/50/100 % (между 50 и 100 % v1 не проверял ничего)",
               "смешение - по J каждой нагрузки (согласованная модель); закон затухания при малом J - оценка",
               "строгий вариант (все зоны с запасом >= 1, как load_range v1) пишется рядом: ok_strict_all_zones",
               "каналы: верхняя оценка с поджиганием группы B горячим газом группы A пишется рядом (value_cross, "
               "ok_cross) - отбор без неё; вне области облика (DOMAIN, оценка) - пометка in_domain = false"])
    return dict(
        schema="kv.combustor.targets/2", title="Цели камеры сгорания KV, v2: 13 целей v1 + низкие выбросы 30-100 % на "
        "одном клапане", date=DATE, base_file="routes/gtu/combustor/targets.json", base_date=T.get("date"),
        base_targets=[t["id"] for t in T["targets"]], new_targets=[new],
        base_notes=dict(load_range=("stable_pct [0, 100]: v1 проверял устойчивость только на 30-100 %; v2 добавляет "
                                    "холостой ход без нагрузки при полной частоте (0 %) - 1/запас по срыву опорной зоны; "
                                    "границы 0 % - линейная экстраполяция цикла по 50 и 30 % (оценка, модель цикла ниже "
                                    "доли топлива 0,005 не считает); правило ОБРАТНЫЙ_ТОК - тоже и на 0 % (проверка "
                                    "скептика 15.11.2026)")),
        funnel_interface=dict(hard=list(T["funnel_interface"]["hard"]) + [NEW_TARGET],
                              objectives=T["funnel_interface"]["objectives"],
                              note="жёсткие v1 + новая; ранжирование - как v1, оценка умножается на (0,5 + 0,5·TRL/9) - "
                                   "решение владельца 15.11.2026, п. 1"),
        trl_multiplier=dict(formula="0,5 + 0,5·TRL/9", trl="наименьший TRL из схемы, пассивного элемента, подачи дизеля, "
                            "охлаждения и способа изготовления (все TRL - оценка)",
                            passive=TRL_PASSIVE, scheme={k: CS.TRL_SCHEME[k] for k in SCHEME_NAMES}))


# =========================================================================== модуль: детали пассивных элементов
def passive_parts(c):
    """отдельные детали пассивного элемента: в печатном картридже и печатной головке гибрида контуры, пористая вставка
    и эжекторы - внутри детали (0); у литой камеры - отдельные коллектор и вставка (оценка)"""
    if c["make"] != "литьё":
        return {}
    N = CS.N_CAN[c["chamber"]]
    out = {}
    if c["split"] == "гидравлика":
        out.update({"второй топливный коллектор (пилот, полость)": 1, "пористая вставка основного контура": 1})
    if c["split"] == "каналы":
        out.update({"второй топливный коллектор (группа A)": 1, "пористая вставка группы B": 1})
    if c["recirc"] == "рециркуляция":
        out["каналы рециркуляции (эжекторы)"] = N
    return out


def _reset_kin():
    """одно начальное состояние решений Cantera перед каждым вариантом: equilibrate и PSR стартуют с текущего
    состояния объекта, и без сброса итог варианта зависел от того, что процесс считал до него (~1e-5 отн.)"""
    K = CR.kin()
    for g in (K.GIN, K.GR, K.GX, K.GP, K.GM):
        g.TPX = 300.0, 101325.0, CR.AIR
    K.DZ.TPX = 300.0, 101325.0, CR.AIR_REITZ
    K.cache.clear()


def evaluate(i, n, kv="KV-3", seed=SEED):
    t0 = time.perf_counter()
    c = concept(i, n, kv, seed)
    _reset_kin()
    row = dict(idx=i, concept=c, status="ok", reject=[])
    try:
        bc = boundary()
        d = derive(c, bc)
        row["derived"] = {k: v for k, v in d.items() if k not in ("zones", "loads", "circuits")}
        row["circuits"] = d["circuits"]
        row["zones"] = d["zones"]
        row["load_states"] = [dict({k: v for k, v in s.items() if k != "circ"},
                                   circ={kk: {x: vv[x] for x in ("share", "phi", "phi_reached", "J", "C", "S", "a_eff",
                                                                  "dpf", "K_decay")} for kk, vv in s["circ"].items()})
                              for s in d["loads"]] if d["loads"] else None
        if d.get("idle"):
            s0 = d["idle"]
            row["idle_state"] = dict({k: v for k, v in s0.items() if k != "circ"},
                                     circ={kk: {x: vv[x] for x in ("share", "phi", "dpf")} for kk, vv in s0["circ"].items()})
        row["domain"] = domain(d)
        bad = check(c, d, bc)
        if bad:
            row.update(status="отсев", reject=[dict(code=b[0], value=b[1], limit=b[2], why=b[3]) for b in bad])
            row["t_s"] = time.perf_counter() - t0
            return CR._rnd(row)
        fb = flashback_v2(c, d, bc)
        row["flashback"] = fb
        gas = [gas_network(c, d, bc, li) for li in range(len(LOADS))]
        row["gas"] = gas
        if c["split"] == "каналы":                          # верхняя оценка: группа B поджигается горячим газом группы A
            keep = ("load", "NOx", "CO", "eta_c", "margin_lbo", "margin_all", "anchor", "lit_main", "tau_ms")
            row["gas_cross"] = [{k: g[k] for k in keep if k in g}
                                for g in (gas_network(c, d, bc, li, cross=True) for li in range(len(LOADS)))]
        gi = gas_network(c, d, with_idle(bc), IDLE_LI, st=d["idle"])     # холостой ход (0 %, оценка): устойчивость
        row["idle"] = {k: gi[k] for k in ("load", "margin_lbo", "margin_all", "margins", "anchor", "eta_c", "CO", "NOx",
                                          "backflow_margin", "p_sup_ratio", "tau_ms") if k in gi}
        g100 = gas[0]
        rej = []
        if (g100.get("margin_all") or 0.0) < 1.0:
            rej.append(dict(code="ПОГАСАНИЕ_100", value=g100.get("margin_all"), limit=1.0,
                            why="одна из зон гаснет на 100 % (время пребывания меньше времени погасания)"))
        if fb["margin_CH4"] is not None and fb["margin_CH4"] < 1.0:
            rej.append(dict(code="ПРОСКОК", value=fb["margin_CH4"], limit=1.0,
                            why=f"проскок пламени в смеситель на метане ({fb.get('margin_CH4_where')})"))
        dz = []
        for ld in (1.0, 0.5):                               # времена зон этой нагрузки при неизменном объёме (оценка)
            li = li_of(ld)
            c2, d2 = _c_at_load(c, d, zone_taus(c, d, bc, li, d["loads"][li]))
            dz.append(CR.diesel_network(c2, d2, bc, li, gas[li]))
        row["diesel"] = dz
        ig = [r.get("ign_margin") for r in dz if r.get("ign_margin") is not None]
        if ig and min(ig) < CR.IGN_MARGIN_MIN:
            rej.append(dict(code="САМОВОСПЛАМЕНЕНИЕ", value=min(ig), limit=CR.IGN_MARGIN_MIN,
                            why="дизель может самовоспламениться в смесителе (задержка меньше 3 времён смешения)"))
        L0 = bc["loads"][0]
        far_pz = L0["W_fuel"] / (c["a_pz"] * L0["W_air"])
        row["wall_gas"] = CR.wall_temperature(c, d, bc, g100["T_flame"], min(far_pz, bc["far_st_CH4"]), 1.0)
        far_dz = L0["W_fuel_diesel"] / (c["a_pz"] * L0["W_air"])
        row["wall_diesel"] = CR.wall_temperature(c, d, bc, dz[0]["T_hot"], min(far_dz, bc["far_st_diesel"]),
                                                 336.0 / CR.H_FUEL_DIESEL ** 2)
        row["pf"] = CR.pattern_factor(c, d, bc, g100)
        row["dp"] = CR.pressure_loss(c, d, bc, g100)
        mod = CR.module_models(c, d, bc, fb, d, row["wall_gas"])
        extra = passive_parts(c)
        mod["parts_list"].update(extra)
        mod["parts"] = int(mod["parts"] + sum(extra.values()))
        mod["actuators"] = 0
        mod["trl_passive"] = trl_passive(c)
        mod["trl_min_v2"] = min(mod["trl_min"], mod["trl_passive"])
        row["module"] = mod
        mlbo = max(min(r["margin_lbo"] or 0.0 for r in gas), 1e-3)
        dpf_min = min(r["dpf"] for r in gas)
        f_stiff = max(1.0, DPF_REF / max(dpf_min, 1e-4)) ** 0.5
        row["dyn"] = dict(margin_min=mlbo, dpf_min=dpf_min, f_stiff=f_stiff)
        row["dyn_index"] = CR.DYN_F[c["scheme"]] * (1.5 / min(mlbo, 3.0)) * (0.03 / c["dp_inj"]) ** 0.5 * \
            (1.5 if c["swirl"] > 0.8 else 1.0) * f_stiff
        row["targets"], hard = compliance(c, row)
        row["hard_ok"] = all(row["targets"][t]["ok"] for t in hard)
        row["n_hard_fail"] = sum(not row["targets"][t]["ok"] for t in hard)
        if rej:
            row.update(status="отсев", reject=rej)
        row["feed"] = dict(dp_comb=row["dp"]["dp_total"], D_comb=d["D_env"], L_comb=d["L_mod"], m_comb=mod["mass_kg"],
                           n_inj=c["n_el"], parts=mod["parts"], **{f"NOx_{int(ld * 100)}": g["NOx"] for ld, g in zip(LOADS, gas)},
                           **{f"CO_{int(ld * 100)}": g["CO"] for ld, g in zip(LOADS, gas)}, PF=row["pf"]["PF"],
                           RTDF=row["pf"]["RTDF"], T_metal_max_K=max(row["wall_gas"]["T_w1"], row["wall_diesel"]["T_w1"]),
                           fuel_valves=max(mod["valves_gas"], mod["valves_diesel"]), H2_vol_max=fb["H2_vol_max"],
                           diesel_NOx_100=dz[0]["NOx"], diesel_NOx_50=dz[1]["NOx"], cartridge_swap_h=mod["swap_h"],
                           V_L=d["V_L"], passive=c["passive"], trl_min=mod["trl_min_v2"],
                           burnt=[bool(r["eta_c"] >= CR.ETA_MIN) for r in gas])
    except Exception as e:
        row.update(status="ошибка", error=f"{type(e).__name__}: {str(e)[:300]}")
    row["t_s"] = time.perf_counter() - t0
    return CR._rnd(row)


def _selfcheck():
    n = 64
    a = [concept(i, n) for i in range(n)]
    assert a == [concept(i, n) for i in range(n)]
    from collections import Counter
    print(json.dumps(dict(fingerprint=fingerprint(n), code_hash=code_hash(), v1_code_hash=CS.code_hash(),
                          schemes=Counter(x["scheme"] for x in a), passive=Counter(x["passive"] for x in a)),
                     ensure_ascii=False))


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--bc", action="store_true", help="пересчитать границы v2 (comb_bc_v2_results.json)")
    ap.add_argument("--targets", action="store_true", help="записать targets_v2.json")
    ap.add_argument("--one", type=int)
    ap.add_argument("--n", type=int, default=300)
    ap.add_argument("--selfcheck", action="store_true")
    a = ap.parse_args()
    if a.bc:
        if os.path.exists(BC_FILE):
            os.remove(BC_FILE)
        b = boundary()
        print("границы v2:", BC_FILE, b["runtime_s"], "с;", json.dumps(b["check"]["point_50_recomputed"], ensure_ascii=False))
    if a.targets:
        json.dump(targets_v2(), open(TARGETS_FILE, "w", encoding="utf-8"), ensure_ascii=False, indent=1)
        print("записано", TARGETS_FILE)
    if a.selfcheck:
        _selfcheck()
    if a.one is not None:
        print(json.dumps(evaluate(a.one, a.n), ensure_ascii=False, indent=1)[:8000])
