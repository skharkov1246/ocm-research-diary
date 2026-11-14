"""Воронка двигателя-заменителя, уровень 0: пространство поиска, генератор выборки, правила допустимости.

Задача владельца (14.11.2026): переизобрести двигатель-заменитель для пакетов Solar и Siemens до 20 МВт (линейка
KV-1 ... KV-7, ../gtu_lineup_results.json). Воронка: 20 000 концепций KV-3 и по 5 000 для KV-1 и KV-6 быстрыми моделями,
затем «лесенкой вниз». Все уровни идут через AWS: код делит выборку на части строго по индексу варианта.

Что здесь:
  * оси пространства (AXES) и архетипы (ARCH) - шаблоны допустимых значений дискретных осей;
  * генератор: скремблированная последовательность Соболя (scipy.stats.qmc) с фиксированным зерном SEED на любое N;
    вариант i выдаётся по индексу (concept), часть k из S получает индексы i % S == k (shard_indices) - без пересечений
    и пропусков на любой машине; отпечаток матрицы (fingerprint) пишется в каждую часть, слияние сверяет его;
  * быстрая оценка цикла уровня 0 (quick_perf) на таблицах свойств газа gtu_explore (Cantera, GRI-Mech 3.0, метан, ISO);
  * геометрия (derive): частота вала, диаметры и ширины колёс, коэффициенты расхода, турбина и силовая турбина, габарит,
    массы модулей, детали и способы их изготовления; вход моделей parts/service/plugplay/make;
  * правила допустимости (check): окружная скорость по материалу и способу изготовления, предел отливки, поковки,
    керамики и камеры принтера, T41 по классу охлаждения и сплаву, частота вала = стык пакета (одновальные - через
    коэффициент расхода колеса при частоте пакета), нагрузка ступеней турбины, мягкий допуск по расходу и выхлопу,
    габарит пакета. Каждое отсеивание - с кодом правила, величиной и пределом.

Запуск:
  python funnel_space.py --kv KV-3 --n 20000 --shard 3 --shards 20 --out "$GTU_OUT"   # одна часть (на AWS)
  python funnel_space.py --selfcheck                                                    # отладка на малой выборке
Результаты ниже уровня 0 - «оценка»: КПД узлов, пределы и коэффициенты - принятые величины трека и свои оценки (README).
"""
import argparse, datetime, hashlib, inspect, json, math, os, platform, sys, time, warnings

import numpy as np
import scipy
from scipy.stats import qmc

HERE = os.path.dirname(os.path.abspath(__file__))
GTU = os.path.dirname(HERE)
if HERE not in sys.path:
    sys.path.insert(0, HERE)
import make_model as MK          # noqa: E402  таблицы фактов и пределы
import parts_model as PM         # noqa: E402
import service_model as SM       # noqa: E402
import plugplay_model as PP      # noqa: E402

SEED = 20261114
DATE = "2026-11-14"
RESULTS = os.path.join(HERE, "gtu_funnel_space_results.json")

# =========================================================================== расчётные точки
KV = {"KV-3": ("Taurus 60", "генератор"), "KV-1": ("Saturn 20", "генератор"), "KV-6": ("Titan 130", "генератор")}
FLEET = {"Taurus 60": "gtu_fleet_taurus60.json", "Saturn 20": "gtu_fleet_saturn20.json", "Titan 130": "gtu_fleet_titan130.json"}
_TGT = {}


def target(kv):
    """точка из gtu_family_results.json (slots: model + variant), масса двигателя - из паспорта, если есть"""
    if kv in _TGT:
        return _TGT[kv]
    model, variant = KV[kv]
    fam = json.load(open(os.path.join(GTU, "gtu_family_results.json"), encoding="utf-8"))
    s = next(x for x in fam["slots"] if x["model"] == model and x["variant"] == variant)
    fl = json.load(open(os.path.join(GTU, "fleet", FLEET[model]), encoding="utf-8"))
    kg = [r.get("engine_kg") for r in fl["records"] if r.get("engine_kg")]
    t = dict(kv=kv, model=model, variant=variant, P_shaft_kW=float(s["P_shaft_kW"]), P_kW=float(s["P_kW"]), G0=float(s["G_air"]),
             EGT0=float(s["EGT0"]), eta0=float(s["eta0"]), n_pkg=float(s["rpm0"]), pi0=float(s["PR0"]),
             pkg_L=s["pkg"][0], pkg_W=s["pkg"][1], pkg_H=s["pkg"][2], rotors0=1 if s["shafts0"] == "single" else 2,
             engine_kg0=float(kg[0]) if kg else None,
             engine_kg0_note=("масса двигателя другого исполнения той же модели (паспорт fleet)" if kg else "нет в паспорте"))
    _TGT[kv] = t
    return t


# =========================================================================== оси и архетипы
CAT = {
    "arch": None,  # заполняется ключами ARCH
    "rotors": [1, 2],
    "comp": ["1Ц", "2Ц", "Д", "ОЦ"],
    "n_ax": [1, 2],
    "inlet": ["одно", "двух"],
    "turb": ["Р", "О"],
    "nt": [1, 2, 3],
    "npt": [1, 2],
    "cool": ["U", "N", "K", "P", "C"],
    "recup": ["нет", "печатный"],
    "comb": ["выносная", "кольцевая", "RDC", "ВР"],
    "bear": ["масло", "магнит"],
    "cores": [1, 2, 4],
    "rotor_join": ["заодно", "сварка", "сборный"],
    "m_cw": ["поковка", "литьё", "печать"],
    "m_tw": ["литьё", "печать", "керамика"],
    "m_hot": ["литьё", "печать", "керамика"],
    "m_comb": ["литьё", "печать", "обработка", "керамика"],
    "m_case": ["литьё", "печать", "обработка"],
    "mat_cw": ["Ti64", "Al", "Ni"],
    "mat_tw": ["K418", "CM247LC", "Haynes282", "IN718", "SC", "Si3N4"],
    "mat_hot": ["K438", "HastX", "IN625", "Haynes282", "CMC"],
    "hot_integral": [0, 1],
    "recup_integral": [0, 1],
    "conn": ["болты", "быстроразъём"],
    "cartridge": [0, 1],
    "split_bl": [0, 1],
}
CONT = {   # непрерывные оси: доля 0..1 переводится в величину в concept() (диапазоны могут зависеть от дискретных осей)
    "pi": "общая степень сжатия: 1Ц 3-8, Д 2,5-6, 2Ц 6-20, ОЦ (1 осевая) 4-14, ОЦ (2 осевые) 5-18 (оценка); верхняя граница "
          "не выше pi_max_U - π, при которой ступени укладываются в предел окружной скорости материала при выбранных mu1, mu2, split",
    "split": "доля ln π на первую ступень (2Ц) или на осевую часть (ОЦ): 0,35-0,65",
    "t41": "T41, °C: от 850 до предела класса охлаждения, типа турбины и сплава (make_model.t41_limit); плотность смещена к пределу: "
           "T41 = Tlim - (Tlim - 850)·u² (низкая T41 не даёт мощности при мягком допуске по расходу)",
    "mu1": "коэффициент работы первой ступени w/U²: центробежная 0,60-0,80, диагональная 0,50-0,70, осевая (по периферии) 0,20-0,35 (оценка)",
    "mu2": "то же для последней ступени (2Ц, ОЦ)",
    "phi1": "коэффициент расхода G/(ρ·U2·D2²) первой центробежной ступени свободного вала: 0,04-0,12; у вала с частотой пакета - выводится",
    "kcm": "cm2/U2 на выходе колеса - ширина b2: 0,20-0,35",
    "zc": "число основных лопаток колеса компрессора: 7-15",
    "kUt": "окружная скорость турбины в долях предела: 0,75-1,00; предел - меньшее из предела материала и (одновальная схема) "
           "π·n·Dmax/60, где Dmax - предел отливки, поковки, керамики или камеры принтера",
    "zt": "лопаток рабочего колеса турбины: радиальная 9-17, осевая 30-80 на ступень",
    "zn": "лопаток соплового аппарата: радиальная 12-24, осевая 20-60 на ступень",
    "kUpt": "окружная скорость силовой турбины в долях 380 м/с: 0,85-1,00 (ниже - лопатки последней ступени слишком длинные при частоте пакета)",
    "phit": "ca/U на выходе осевой турбины и СТ: 0,45-0,80 (frame3 - 0,5); выше 0,5 - потеря с выходной скоростью: половина прироста c²/2 (оценка)",
    "eps_rec": "эффективность рекуператора: 0,75-0,92",
    "pgain": "прирост полного давления в камере: RDC 0-0,10, волновой ротор 0-0,15 (оценка; модели обещают +3-5 % КПД у RDC - по выдержке поиска)",
    "kD_case": "диаметр корпуса / наибольший D2: 1,7-2,1 (gtu_family - 1,9)",
}
_DEF = dict(rotors=[1, 2], comp=["1Ц", "2Ц", "Д", "ОЦ"], n_ax=[1, 2], inlet=["одно", "двух"], turb=["Р", "О"], nt=[1, 2],
            npt=[1, 2], cool=["U", "N", "K"], recup=["нет", "печатный"], comb=["выносная", "кольцевая"], bear=["масло", "магнит"],
            cores=[1], rotor_join=["сварка", "сборный"], m_cw=["поковка", "литьё", "печать"], m_tw=["литьё", "печать"],
            m_hot=["литьё", "печать"], m_comb=["литьё", "печать", "обработка"], m_case=["литьё", "печать", "обработка"],
            mat_cw=["Ti64", "Al"], hot_integral=[0], recup_integral=[0], conn=["болты", "быстроразъём"], cartridge=[0, 1],
            split_bl=[0, 1])
ARCH = {
    "традиционная": dict(name="традиционная сборная (опорная точка frame3)", allowed=dict(
        rotors=[1], comp=["2Ц"], inlet=["двух"], turb=["О"], nt=[3], cool=["U"], recup=["нет"], comb=["кольцевая"], bear=["масло"],
        rotor_join=["сборный"], m_cw=["поковка"], mat_cw=["Ti64"], m_tw=["литьё"], mat_tw=["SC"], m_hot=["литьё"], m_comb=["обработка"],
        m_case=["литьё"], conn=["болты"], cartridge=[0]),
        note="три осевые ступени - только здесь; остальные архетипы - 1-2 ступени турбины"),
    "печатный монолит": dict(name="печатный монолит: ротор одной деталью, горячая часть одной деталью", allowed=dict(
        comp=["1Ц", "2Ц", "Д"], cool=["U", "P"], comb=["кольцевая"], rotor_join=["заодно"], m_tw=["печать"],
        mat_tw=["IN718", "Haynes282", "CM247LC"], m_hot=["печать"], m_comb=["печать"], m_case=["печать", "литьё"],
        hot_integral=[1], cartridge=[1])),
    "литой монолит": dict(name="литой монолит: ротор заодно или сварка трением двух отливок", allowed=dict(
        comp=["1Ц", "2Ц"], turb=["Р"], cool=["U", "N", "K"], rotor_join=["заодно", "сварка"], m_tw=["литьё"], mat_tw=["K418", "CM247LC"],
        m_hot=["литьё"], m_comb=["литьё", "обработка"], m_case=["литьё"], hot_integral=[0, 1])),
    "гибрид": dict(name="гибрид литьё + печать: литые колёса и корпуса, печатная горячая часть с печатным охлаждением", allowed=dict(
        cool=["P", "N"], rotor_join=["сварка", "сборный"], m_cw=["поковка", "литьё"], mat_cw=["Ti64"], m_tw=["литьё"],
        m_hot=["печать"], m_comb=["печать"], m_case=["литьё"], hot_integral=[0, 1])),
    "печатный рекуператор": dict(name="печатный рекуператор в корпусе", allowed=dict(
        comp=["1Ц", "2Ц", "Д"], cool=["U", "N", "P"], recup=["печатный"], recup_integral=[1], m_case=["литьё", "печать"])),
    "керамика": dict(name="керамическая горячая часть (Si3N4, CMC)", allowed=dict(
        comp=["1Ц", "2Ц"], turb=["Р"], cool=["C"], rotor_join=["сборный"], m_cw=["поковка", "литьё"], mat_cw=["Ti64"],
        m_tw=["керамика"], mat_tw=["Si3N4"], m_hot=["керамика"], mat_hot=["CMC"], m_comb=["керамика", "печать"])),
    "кластер": dict(name="кластер малых ядер + общая силовая турбина", allowed=dict(
        cores=[2, 4], rotors=[2], comp=["1Ц", "2Ц"], inlet=["одно"], turb=["Р"], cool=["U", "N", "K", "P"], recup=["нет"],
        rotor_join=["заодно", "сварка", "сборный"], m_tw=["литьё", "печать"])),
    "простой ротор": dict(name="предельно простой ротор: одно трансзвуковое колесо + одна турбина", allowed=dict(
        rotors=[1], comp=["1Ц"], nt=[1], cool=["U", "N", "K", "P"], rotor_join=["заодно", "сварка", "сборный"])),
    "безмасляный": dict(name="безмасляный ротор на магнитных подшипниках", allowed=dict(
        bear=["магнит"], cool=["U", "N", "K", "P"], rotor_join=["заодно", "сварка", "сборный"])),
    "RDC": dict(name="горение с ростом давления (вращающаяся детонация)", allowed=dict(
        comp=["1Ц", "2Ц", "Д"], cool=["N", "K", "P"], recup=["нет"], comb=["RDC"], m_comb=["печать"])),
    "волновой ротор": dict(name="волновой ротор (ударная волна сжимает газ перед турбиной)", allowed=dict(
        comp=["1Ц", "Д"], cool=["U", "N", "K", "P"], recup=["нет"], comb=["ВР"], m_comb=["печать", "литьё"])),
}
ARCH_KEYS = list(ARCH)
CAT["arch"] = ARCH_KEYS
AXES = list(CAT) + list(CONT)
MAT_TW = {"литьё": ["K418", "CM247LC", "SC"], "печать": ["IN718", "Haynes282", "CM247LC"], "керамика": ["Si3N4"]}
MAT_HOT = {"литьё": ["K438"], "печать": ["HastX", "IN625", "Haynes282"], "керамика": ["CMC"]}
PI_RANGE = {"1Ц": (3.0, 8.0), "Д": (2.5, 6.0), "2Ц": (6.0, 20.0), ("ОЦ", 1): (4.0, 14.0), ("ОЦ", 2): (5.0, 18.0)}
MU_RANGE = {"Ц": (0.60, 0.80), "Д": (0.50, 0.70), "О": (0.20, 0.35)}
T41_MIN = 850.0


def _pick(u, vals):
    return vals[min(int(u * len(vals)), len(vals) - 1)]


def _lerp(u, lo, hi):
    return lo + (hi - lo) * float(u)


# =========================================================================== генератор
_DM = {}


def _rng_kw(fn, seed):
    return {"rng": seed} if "rng" in inspect.signature(fn).parameters else {"seed": seed}


def design_matrix(n, seed=SEED):
    """матрица n x len(AXES) из скремблированной последовательности Соболя: первые n из 2^m точек"""
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


def stage_types(c):
    return {"1Ц": ["Ц"], "Д": ["Д"], "2Ц": ["Ц", "Ц"], "ОЦ": ["О"] * max(c.get("n_ax", 1), 1) + ["Ц"]}[c["comp"]]


def concept(i, n, kv="KV-3", seed=SEED):
    """концепция с индексом i выборки объёма n: значения всех осей (детерминированно по i, n, seed)"""
    if not 0 <= i < n:
        raise IndexError(i)
    u = dict(zip(AXES, design_matrix(n, seed)[i]))
    arch = _pick(u["arch"], ARCH_KEYS)
    T = ARCH[arch]["allowed"]

    def cat(name, vals=None):
        v = vals if vals is not None else T.get(name, _DEF[name])
        return _pick(u[name], v)

    c = dict(idx=int(i), n=int(n), seed=int(seed), kv=kv, arch=arch)
    c["cores"] = cat("cores")
    c["rotors"] = 2 if c["cores"] > 1 else cat("rotors")
    c["comp"] = cat("comp")
    c["n_ax"] = cat("n_ax") if c["comp"] == "ОЦ" else 0
    c["inlet"] = cat("inlet")
    c["turb"] = cat("turb")
    c["nt"] = cat("nt")
    c["npt"] = cat("npt") if c["rotors"] == 2 else 0
    c["cool"] = cat("cool")
    c["recup"] = cat("recup")
    c["comb"] = cat("comb")
    c["bear"] = cat("bear")
    c["m_tw"] = cat("m_tw")
    c["rotor_join"] = "сборный" if c["m_tw"] == "керамика" else cat("rotor_join")
    if c["rotor_join"] in ("заодно", "сварка"):
        c["m_cw"], c["mat_cw"] = c["m_tw"], "Ni"
    else:
        c["m_cw"] = cat("m_cw")
        c["mat_cw"] = cat("mat_cw", [m for m in T.get("mat_cw", _DEF["mat_cw"]) if m != "Ni"])
    tw = [m for m in MAT_TW[c["m_tw"]] if not (m == "SC" and c["turb"] != "О") and not (m == "Si3N4" and c["turb"] != "Р")]
    if "mat_tw" in T:
        tw = [m for m in tw if m in T["mat_tw"]] or tw
    c["mat_tw"] = cat("mat_tw", tw)
    c["m_hot"] = "печать" if c["cool"] == "P" else cat("m_hot")
    hv = MAT_HOT[c["m_hot"]]
    if "mat_hot" in T:
        hv = [m for m in hv if m in T["mat_hot"]] or hv
    c["mat_hot"] = cat("mat_hot", hv)
    c["m_comb"] = cat("m_comb")
    c["m_case"] = cat("m_case")
    c["hot_integral"] = cat("hot_integral")
    c["recup_integral"] = cat("recup_integral") if c["recup"] == "печатный" else 0
    c["conn"] = cat("conn")
    c["cartridge"] = cat("cartridge")
    c["split_bl"] = cat("split_bl")
    # непрерывные
    st = stage_types(c)
    c["mu1"] = _lerp(u["mu1"], *MU_RANGE[st[0]])
    c["mu2"] = _lerp(u["mu2"], *MU_RANGE[st[-1]])
    c["split"] = _lerp(u["split"], 0.35, 0.65)
    lo, hi = PI_RANGE[(c["comp"], c["n_ax"])] if c["comp"] == "ОЦ" else PI_RANGE[c["comp"]]
    c["pi_max_U"] = round(pi_max_by_speed(c), 4)
    hi_eff = min(hi, c["pi_max_U"])
    c["pi"] = _lerp(u["pi"], lo, hi_eff) if hi_eff > lo else lo
    tl = MK.t41_limit(c["turb"], c["cool"], c["mat_tw"], c["m_tw"])[0]
    c["t41"] = tl - (tl - T41_MIN) * float(u["t41"]) ** 2 if tl > T41_MIN else T41_MIN
    c["phi1"] = _lerp(u["phi1"], 0.04, 0.12)
    c["kcm"] = _lerp(u["kcm"], 0.20, 0.35)
    c["zc"] = int(round(_lerp(u["zc"], 7, 15)))
    c["kUt"] = _lerp(u["kUt"], 0.75, 1.00)
    c["zt"] = int(round(_lerp(u["zt"], *((9, 17) if c["turb"] == "Р" else (30, 80)))))
    c["zn"] = int(round(_lerp(u["zn"], *((12, 24) if c["turb"] == "Р" else (20, 60)))))
    c["kUpt"] = _lerp(u["kUpt"], 0.85, 1.00) if c["rotors"] == 2 else 0.0
    c["phit"] = _lerp(u["phit"], 0.45, 0.80)
    c["eps_rec"] = _lerp(u["eps_rec"], 0.75, 0.92) if c["recup"] == "печатный" else 0.0
    c["pgain"] = _lerp(u["pgain"], 0.0, {"RDC": 0.10, "ВР": 0.15}.get(c["comb"], 0.0))
    c["kD_case"] = _lerp(u["kD_case"], 1.7, 2.1)
    for k in CONT:
        if isinstance(c[k], float):
            c[k] = round(c[k], 6)
    return c


# =========================================================================== быстрая оценка цикла (уровень 0)
_EX = None


def _E():
    """таблицы свойств газа gtu_explore (импорт делает chdir в папку трека - возвращаемся)"""
    global _EX
    if _EX is None:
        cwd = os.getcwd()
        if GTU not in sys.path:
            sys.path.insert(0, GTU)
        import gtu_explore as ex
        os.chdir(cwd)
        _EX = ex
    return _EX


U_MARGIN, ETA_MARGIN = 0.97, 0.03     # запас при выборе верхней границы π: доля предела скорости и снижение КПД (оценка)


def d_cap_comp(c):
    """наибольший диаметр колеса компрессора по способу изготовления (make_model: CAST, FORGE, принтеры), м"""
    j, m, mat = c["rotor_join"], c["m_cw"], c["mat_cw"]
    if m == "печать":
        return MK.largest_print_dim()
    if j == "заодно":
        return MK.CAST["rotor_mono"]["max"]
    if j == "сварка":
        return MK.CAST["rotor_section"]["max"]
    kind = "wheel_ti" if mat == "Ti64" else "wheel_al"
    return (MK.FORGE if m == "поковка" else MK.CAST)[kind]["max"]


def d_cap_turb(c):
    """наибольший диаметр колеса турбины (радиальная - периферия; осевая - средний по ротору заодно), м"""
    m, j = c["m_tw"], c["rotor_join"]
    if m == "печать":
        return MK.largest_print_dim()
    if m == "керамика":
        return MK.CERAMIC["turb_radial"]["max"]
    if j == "заодно":
        return MK.CAST["rotor_mono"]["max"]
    if c["turb"] == "Р":
        return MK.CAST["turb_radial"]["max"] if j == "сборный" else MK.CAST["rotor_section"]["max"]
    return None


def u_cap_size(c, D):
    """при частоте пакета (одновальная схема) предел диаметра - это предел окружной скорости π·n·D/60"""
    if D is None or c["rotors"] != 1:
        return float("inf")
    return math.pi * target(c["kv"])["n_pkg"] * D / 60.0


def pi_max_by_speed(c):
    """наибольшая общая π, при которой каждая ступень укладывается в предел окружной скорости своего материала при
    выбранных коэффициентах работы mu1/mu2 и распределении split (с запасом U_MARGIN, ETA_MARGIN). Нужна, чтобы
    выборка не тратилась на заведомо превышенную окружную скорость; правило U_comp всё равно проверяется."""
    E = _E()
    typ = stage_types(c)
    mus = stage_mus(c)
    T = E.T0; pmax = []
    for k, (t, mu) in enumerate(zip(typ, mus)):
        mat = c["mat_cw"] if not (c["mat_cw"] == "Al" and k > 0) else "Ti64"
        Ul = min(MK.u_limit_comp(mat, c["m_cw"])[0], u_cap_size(c, d_cap_comp(c))) * U_MARGIN
        T2 = E.T_h(E.h(T) + mu * Ul ** 2)
        pis = 3.0
        for _ in range(3):
            eta = eta_comp(t, pis, mat, k, c["inlet"], mu) - ETA_MARGIN
            pis = math.exp(eta * (E.s0(T2) - E.s0(T)) / E.RA)
        if t == "О":
            pis = min(pis, AX_PI_MAX)
        if t == "Д":
            pis = min(pis, D_PI_MAX)
        pmax.append(pis); T = T2
    s = c["split"]
    if c["comp"] in ("1Ц", "Д"):
        return pmax[0]
    if c["comp"] == "2Ц":
        return min(pmax[0] ** (1 / s), pmax[1] ** (1 / (1 - s)))
    na = c["n_ax"]
    return min(min(pmax[:na]) ** (na / s), pmax[-1] ** (1 / (1 - s)))


def stage_mus(c):
    typ = stage_types(c)
    if c["comp"] == "ОЦ":
        return [c["mu1"]] * c["n_ax"] + [c["mu2"]]
    return [c["mu1"]] if len(typ) == 1 else [c["mu1"], c["mu2"]]


def stages_pi(c):
    st, pi, s = stage_types(c), c["pi"], c["split"]
    if c["comp"] in ("1Ц", "Д"):
        return [pi]
    if c["comp"] == "2Ц":
        return [pi ** s, pi ** (1 - s)]
    na = c["n_ax"]
    return [pi ** (s / na)] * na + [pi ** (1 - s)]


MU_ETA = 0.10     # ΔКПД = -0,10·(mu - 0,70): меньше загнутость назад - больше работа, но ниже КПД (оценка; comp_stage - mu 0,70)


def eta_comp(typ, pis, mat, k, inlet, mu=None):
    """политропный КПД ступени: формула comp_stage (gtu_explore) для центробежной; поправки - оценка"""
    if typ == "О":
        e = 0.885 - 0.04 * max(0.0, pis - 1.6)
    else:
        e = 0.880 - 0.010 * (pis - 2.5) - (0.02 if mat == "Ni" else 0.0) - (0.01 if typ == "Д" else 0.0)
        if mu is not None:
            e -= MU_ETA * (mu - 0.70)
    if k == 0 and inlet == "двух":
        e -= 0.01
    return max(0.70, e)


def _comp(stages):
    E = _E(); T, p = E.T0, E.P0; out = []
    for pis, eta in stages:
        T2 = E.T_s(E.s0(T) + E.RA * math.log(pis) / eta)
        out.append(dict(T1=T, p1=p, T2=T2, p2=p * pis, w=E.h(T2) - E.h(T), pi=pis, eta=eta))
        T, p = T2, p * pis
    return out


def cool_air(c):
    xn, xr = MK.COOL_AIR[c["cool"]]
    if not MK.rotor_cooled(c["cool"], c["m_tw"]):
        xr = 0.0
    return xn, xr


def _core(c, comp_out, eta_gg, eta_pt):
    """термодинамика при заданных КПД турбин (политропные). Схема и потери - как gtu_cycle.cycle: охлаждение соплового
    аппарата возвращается до ротора (T41 - после подмешивания), охлаждение ротора - после турбины (работы не даёт), утечка
    теряется; камера 4 % (выносная и кластер - 5 %), переходный канал 1 % (кластер 2,5 %), выхлоп 1,02·P0, мех. КПД 0,99,
    полнота сгорания 0,995; рекуператор - потери 2,5 % по воздуху и 3,5 % по газу (gtu_cycle.recup)."""
    E = _E()
    T2, p2 = comp_out[-1]["T2"], comp_out[-1]["p2"]
    wc = sum(s["w"] for s in comp_out)
    xn, xr = cool_air(c)
    xl = 0.015 + (0.01 if c["cores"] > 1 else 0.0)
    mc = 1.0 - xn - xr - xl
    rec = c["recup"] == "печатный"
    dpa, dpg = (0.025, 0.035) if rec else (0.0, 0.0)
    if c["comb"] in ("RDC", "ВР"):
        pfac = 1.0 + c["pgain"]
    else:
        pfac = 1.0 - (0.05 if (c["cores"] > 1 or c["comb"] == "выносная") else 0.04)
    p4 = p2 * (1 - dpa) * pfac
    p_exit = E.P0 * 1.02 * (1 + dpg)
    dp_duct = 0.025 if c["cores"] > 1 else 0.01
    T41 = c["t41"] + 273.15
    T3 = T2
    for _ in range(10 if rec else 1):
        T4 = T41
        for _ in range(5):
            f = E.far_for(T3, T4); b = E.bfrac(f)
            mg4 = mc * (1 + f); m41 = mg4 + xn; b41 = b * mg4 / m41
            T4 = E.T_h((m41 * E.h(T41, b41) - xn * E.h(T2)) / mg4, b)
        R41, s41, h41 = E.Rg(b41), E.s0(T41, b41), E.h(T41, b41)
        pt = None
        if c["rotors"] == 1:
            ER = p4 / p_exit
            if ER <= 1.05:
                return dict(fail="степень расширения турбины <= 1,05")
            T5 = E.T_s(s41 - eta_gg * R41 * math.log(ER), b41)
            wt = h41 - E.h(T5, b41)
            wnet = 0.99 * m41 * wt - wc
            m5 = m41 + xr; b5 = b41 * m41 / m5
            Tex = E.T_h((m41 * E.h(T5, b41) + xr * E.h(T2)) / m5, b5)
            t = dict(ER=ER, w=wt, T_in=T41, p_in=p4, b=b41, T_out=T5, p_out=p_exit, m=m41)
        else:
            wgg = wc / (0.99 * m41)
            T45 = E.T_h(h41 - wgg, b41)
            ERgg = math.exp((s41 - E.s0(T45, b41)) / (eta_gg * R41))
            p45 = p4 / ERgg * (1 - dp_duct)
            if p45 <= p_exit * 1.15:
                return dict(fail="силовой турбине не остаётся перепада (степень расширения < 1,15)")
            m45 = m41 + xr; b45 = b41 * m41 / m45
            T45m = E.T_h((m41 * E.h(T45, b41) + xr * E.h(T2)) / m45, b45)
            ERpt = p45 / p_exit
            T5 = E.T_s(E.s0(T45m, b45) - eta_pt * E.Rg(b45) * math.log(ERpt), b45)
            wpt = E.h(T45m, b45) - E.h(T5, b45)
            wnet = 0.99 * m45 * wpt
            Tex = T5
            t = dict(ER=ERgg, w=wgg, T_in=T41, p_in=p4, b=b41, T_out=T45, p_out=p4 / ERgg, m=m41)
            pt = dict(ER=ERpt, w=wpt, T_in=T45m, p_in=p45, b=b45, T_out=T5, p_out=p_exit, m=m45)
        if rec:
            T3n = T2 + c["eps_rec"] * (Tex - T2)
            done = abs(T3n - T3) < 0.2
            T3 = T3n
            if done:
                break
    if wnet <= 0:
        return dict(fail="удельная работа <= 0")
    f_tot = f * mc / 0.995
    return dict(wc=wc, wnet=wnet, eta=wnet / (f_tot * E.LHV), f=f_tot, T2=T2, p2=p2, T3=T3, T41=T41, p4=p4,
                T5=Tex, Tstack=(Tex - (T3 - T2) * 0.93) if rec else Tex, t=t, pt=pt, xn=xn, xr=xr, xl=xl)


def turb_speed(c):
    """окружная скорость турбины: kUt x меньшее из предела материала и предела размера при частоте пакета"""
    U, ub = MK.u_limit_turb(c["turb"], c["mat_tw"], c["t41"])
    cap = u_cap_size(c, d_cap_turb(c)) * 0.99
    return c["kUt"] * min(U, cap), U, ub


def _radial_stage(t, nt, U):
    """радиальная турбина: КПД по U/C0 (turb_radial gtu_explore: 0,86 при 0,70); степень расширения делится поровну"""
    E = _E()
    ERs = t["ER"] ** (1.0 / nt)
    dhs = E.h(t["T_in"], t["b"]) - E.h(E.T_s(E.s0(t["T_in"], t["b"]) - E.Rg(t["b"]) * math.log(ERs), t["b"]), t["b"])
    C0 = math.sqrt(2 * max(dhs, 1.0)); nu = U / C0
    return 0.86 * max(0.05, 1 - ((nu - 0.70) / 0.70) ** 2), nu, dhs


LEAVE_K = 0.5      # не восстановленная диффузором доля прироста выходной кинетической энергии (оценка)


def leaving_loss(phit, U, dhs):
    """снижение КПД осевой турбины при ca/U выше 0,5 (eta_axial отвечает 0,5): LEAVE_K·((phit·U)² - (0,5·U)²)/2 / Δh_s"""
    return LEAVE_K * max(0.0, (phit * U) ** 2 - (CX_T * U) ** 2) / 2.0 / max(dhs, 1.0)


def quick_perf(c, tgt, d_eta=None):
    """быстрый цикл концепции: КПД турбин - по нагрузке ступеней при окружной скорости концепции (итерации).
    d_eta - поправки КПД ступеней компрессора от геометрии (коэффициент расхода, зазор), из derive."""
    E = _E()
    pis = stages_pi(c); typ = stage_types(c)
    d_eta = d_eta or [0.0] * len(pis)
    mats = [c["mat_cw"] if not (c["mat_cw"] == "Al" and k > 0) else "Ti64" for k in range(len(pis))]
    mus = stage_mus(c)
    comp = _comp([(p, max(0.60, eta_comp(t, p, m, k, c["inlet"], mu) + de))
                  for k, (p, t, m, de, mu) in enumerate(zip(pis, typ, mats, d_eta, mus))])
    Ut, Ulim, _ = turb_speed(c)
    Upt = c["kUpt"] * 380.0
    eg, ep = 0.86, 0.88
    r, load = None, {}
    for _ in range(6):
        r = _core(c, comp, eg, ep)
        if "fail" in r:
            break
        if c["turb"] == "О":
            psi = r["t"]["w"] / (c["nt"] * Ut ** 2); eg_n = E.eta_axial(psi) - leaving_loss(c["phit"], Ut, r["t"]["w"] / eg)
            load = dict(psi=psi)
        else:
            eg_n, nu, dhs = _radial_stage(r["t"], c["nt"], Ut); load = dict(nu=nu, dhs=dhs)
        ep_n = ep
        if r["pt"]:
            psi_pt = r["pt"]["w"] / (c["npt"] * Upt ** 2)
            ep_n = E.eta_axial(psi_pt) - leaving_loss(c["phit"], Upt, r["pt"]["w"] / ep); load["psi_pt"] = psi_pt
        eg_n, ep_n = max(0.40, eg_n), max(0.40, ep_n)
        if abs(eg_n - eg) < 1e-4 and abs(ep_n - ep) < 1e-4:
            break
        eg, ep = eg_n, ep_n
    if r is None or "fail" in r:
        return dict(ok=False, fail=(r or {}).get("fail", "нет решения"), comp=comp)
    G = tgt["P_shaft_kW"] * 1e3 / r["wnet"]
    return dict(ok=True, comp=comp, r=r, eta_gg=eg, eta_pt=ep if r["pt"] else None, load=load, G=G, Ut=Ut, Ulim_t=Ulim, Upt=Upt,
                eta_th=r["eta"], w=r["wnet"] / 1e3, EGT=r["Tstack"] - 273.15, T5=r["T5"] - 273.15, T2=r["T2"] - 273.15)


# =========================================================================== геометрия и входы моделей
PHI_OPT = {"Ц": 0.085, "Д": 0.12}
PHI_RANGE = {"Ц": (0.02, 0.16), "Д": (0.05, 0.25)}
PHI_PENALTY = 0.024            # ΔКПД = -0,024·ln²(φ/φ_opt): подбор по двум ступеням frame3 (оценка)
REC_K, REC_BETA = 120.0, 1500.0  # Вт/(м2·К) - gtu_explore.finalize; м2/м3 печатного рекуператора - оценка
CX_T = 0.5                       # ca/U на выходе турбины - frame3 (phi_t 0,5)


def _shell(D, L, rho):
    t = max(0.006, 0.012 * D)
    return rho * math.pi * D * L * t


def derive(c, tgt, perf):
    """геометрия по быстрому циклу: частота, колёса, турбины, габарит, массы, детали; входы моделей"""
    E = _E()
    pis, typ = stages_pi(c), stage_types(c)
    G = perf["G"]; Gc = G / c["cores"]
    comp = perf["comp"]
    mats = [c["mat_cw"] if not (c["mat_cw"] == "Al" and k > 0) else "Ti64" for k in range(len(pis))]
    mus = stage_mus(c)
    st = []
    for k, (s, t, m, mu) in enumerate(zip(comp, typ, mats, mus)):
        Gs = Gc / 2 if (k == 0 and c["inlet"] == "двух") else Gc
        st.append(dict(type=t, pi=s["pi"], w=s["w"], U=math.sqrt(s["w"] / mu), mat=m, Gs=Gs, rho1=s["p1"] / (E.RA * s["T1"]),
                       rho2=s["p2"] / (E.RA * s["T2"]), T2=s["T2"] - 273.15, eta=s["eta"]))
    if c["rotors"] == 1:
        n = tgt["n_pkg"]
    else:
        k0 = next(k for k, s in enumerate(st) if s["type"] != "О")
        s0_ = st[k0]
        D2 = math.sqrt(s0_["Gs"] / (s0_["rho1"] * s0_["U"] * c["phi1"]))
        n = 60 * s0_["U"] / (math.pi * D2)
    d_eta = []
    for s in st:
        s["D"] = 60 * s["U"] / (math.pi * n)
        if s["type"] == "О":
            A = s["Gs"] / (s["rho1"] * 0.5 * s["U"])
            x = 1 - A / (math.pi * s["D"] ** 2 / 4)
            s["htr"] = math.sqrt(x) if x > 0 else 0.0
            d_eta.append(0.0)
        else:
            s["phi"] = s["Gs"] / (s["rho1"] * s["U"] * s["D"] ** 2)
            s["b2"] = s["Gs"] / (math.pi * s["D"] * s["rho2"] * c["kcm"] * s["U"])
            cl = max(0.3e-3, 0.0007 * s["D"])
            de = -0.3 * cl / s["b2"] - PHI_PENALTY * math.log(s["phi"] / PHI_OPT[s["type"]]) ** 2
            s["d_eta"] = de
            d_eta.append(de)
    Dc = max(s["D"] for s in st)
    # турбина газогенератора (или единственная)
    r = perf["r"]; t = r["t"]
    Ut = perf["Ut"]; Dt = 60 * Ut / (math.pi * n)
    rho5 = t["p_out"] / (E.Rg(t["b"]) * t["T_out"])
    Q5 = t["m"] * Gc / rho5
    tg = dict(type=c["turb"], nt=c["nt"], U=Ut, Ulim=perf["Ulim_t"], D=Dt, Q5=Q5)
    if c["turb"] == "О":
        h = Q5 / (math.pi * Dt * c["phit"] * Ut)
        tg.update(h=h, htr=(Dt - h) / (Dt + h), Dtip=Dt + h, AN2=math.pi * Dt * h * n ** 2, psi=perf["load"].get("psi"))
    else:
        dhs = perf["load"].get("dhs", 1.0)
        tg.update(h=None, Dtip=Dt, nu=perf["load"].get("nu"), ns=(2 * math.pi * n / 60) * math.sqrt(Q5) / max(dhs, 1.0) ** 0.75)
    # силовая турбина
    pt = None
    if c["rotors"] == 2:
        p = r["pt"]; Upt = perf["Upt"]; Dpt = 60 * Upt / (math.pi * tgt["n_pkg"])
        rho = p["p_out"] / (E.Rg(p["b"]) * p["T_out"]); Q = p["m"] * G / rho
        hp = Q / (math.pi * Dpt * c["phit"] * Upt)
        pt = dict(U=Upt, D=Dpt, h=hp, htr=(Dpt - hp) / (Dpt + hp), Dtip=Dpt + hp, psi=perf["load"].get("psi_pt"), npt=c["npt"])
    # габарит и длины (оценка)
    L_rot = sum(0.35 * s["D"] if s["type"] != "О" else 0.12 * s["D"] for s in st) + \
        (0.35 * Dt * c["nt"] if c["turb"] == "Р" else 0.15 * Dt * c["nt"]) + 0.5 * Dc
    D_case = c["kD_case"] * Dc
    D_env = D_case + 0.5 if c["cores"] == 1 else D_case * math.sqrt(c["cores"]) * 1.1 + 0.5
    if pt:
        D_env = max(D_env, 1.3 * pt["Dtip"] + 0.5)
    rec = None
    if c["recup"] == "печатный":
        eps = c["eps_rec"]; ntu = eps / (1 - eps)
        A = ntu * G * 1100.0 / REC_K
        V = A / REC_BETA
        a_max = 1.1                                        # наибольший куб в камере EP-M2050 (Z 1100 мм)
        mods = max(1, math.ceil(V / a_max ** 3))
        rec = dict(A_m2=A, V_m3=V, modules=mods, side=(V / mods) ** (1 / 3), mass_kg=A * 0.00015 * 7900 * 3)
        if c["recup_integral"]:
            D_env = math.sqrt(D_env ** 2 + 4 * V / (math.pi * max(D_env, 0.5)))
    # массы (формулы gtu_explore.finalize; прочее - оценка)
    rho_cw = {"Ti64": 4430.0, "Al": 2800.0, "Ni": 8000.0}
    m_cw = sum((0.0825 if s["type"] != "О" else 0.03) * rho_cw[s["mat"]] * s["D"] ** 3 for s in st)
    rho_t = MK.MAT.get(c["mat_tw"], {}).get("rho", 8000.0)
    m_tw = (0.0825 * rho_t * Dt ** 3 * 1.2 * c["nt"]) if c["turb"] == "Р" else c["nt"] * 0.03 * rho_t * (1.15 * Dt) ** 3 * 1.5
    m_shaft = 7850.0 * math.pi / 4 * (0.2 * Dc) ** 2 * L_rot * 0.5
    m_rotor = m_cw + m_tw + m_shaft
    rho_case = 4430.0 if c["m_case"] == "печать" else 7850.0
    L_case = 0.6 * L_rot
    m_case = _shell(D_case, L_case, rho_case)
    D_hot = 1.25 * tg["Dtip"]; L_hot = 0.6 * D_hot
    rho_hot = MK.MAT.get(c["mat_hot"], {}).get("rho", 8200.0)
    m_hot = _shell(D_hot, L_hot, rho_hot) + 0.02 * rho_hot * Dt ** 3 * c["nt"]
    if c["comb"] in ("кольцевая", "RDC"):
        D_comb = 1.1 * Dc; L_comb = 0.6 * D_comb
        m_comb = 2 * math.pi * 0.8 * D_comb * L_comb * 0.002 * 8220
        n_inj = max(6, round(math.pi * 0.8 * D_comb / 0.12))
    else:
        rho3 = r["p2"] / (E.RA * r["T3"])
        D_comb = 1.3 * math.sqrt(4 * (Gc / (rho3 * 25.0)) / math.pi); L_comb = 3.0 * D_comb
        m_comb = 2 * math.pi * D_comb * L_comb * 0.002 * 8220 + _shell(D_comb * 1.2, L_comb, 7850.0)
        n_inj = 1
    if c["comb"] == "ВР":
        m_comb += 0.05 * 8000 * Dc ** 3                    # барабан волнового ротора (оценка)
    m_pt = 0.0
    if pt:
        m_pt = c["npt"] * 0.03 * 8000 * (1.15 * pt["D"]) ** 3 * 1.5 + 7850 * math.pi / 4 * (0.15 * pt["D"]) ** 2 * pt["D"] + \
            _shell(1.3 * pt["Dtip"], 0.5 * pt["D"], 7850.0)
    geo = dict(n=n, stages=st, Dc=Dc, tg=tg, pt=pt, L_rot=L_rot, D_case=D_case, D_env=D_env, L_case=L_case, D_hot=D_hot,
               L_hot=L_hot, D_comb=D_comb, L_comb=L_comb, n_inj=n_inj, rec=rec,
               mass=dict(rotor=m_rotor, cw=m_cw, tw=m_tw, shaft=m_shaft, case=m_case, hot=m_hot, comb=m_comb, pt=m_pt),
               d_eta=d_eta)
    geo["parts_mk"], geo["blisk"], geo["pt_blisk"], geo["ngv_ring"], geo["case_pieces"] = _parts_for_make(c, geo)
    return geo


def _blisk(method, mat, Dtip, Dm):
    if mat == "SC":
        return False
    if method == "печать":
        return bool(MK.printers_for(Dtip, 0.15 * Dm, True))
    if method == "литьё":
        return bool(Dtip <= MK.CAST["turb_axial_blisk"]["max"])
    return False


def _parts_for_make(c, g):
    """детали со способом изготовления и размерами для make_model и правила размеров"""
    P = []
    st, tg = g["stages"], g["tg"]
    Dmax_rot = max(g["Dc"], tg["Dtip"])
    blisk = [True] * c["nt"]
    if c["turb"] == "О":
        blisk = [_blisk(c["m_tw"], c["mat_tw"], tg["Dtip"], tg["D"])] * c["nt"]
    jt = c["rotor_join"]; m = g["mass"]
    if jt == "заодно":
        kind = "rotor_mono"
        P.append(dict(name="ротор целиком", kind=kind, method=c["m_tw"], material=c["mat_tw"], D=Dmax_rot, L=g["L_rot"],
                      mass_kg=m["rotor"] - m["shaft"] * 0.5, rotating=True))
    elif jt == "сварка":
        P.append(dict(name="секция компрессора (сварка)", kind="rotor_section", method=c["m_tw"], material=c["mat_tw"], D=g["Dc"],
                      L=0.55 * g["L_rot"], mass_kg=m["cw"] + 0.5 * m["shaft"], rotating=True))
        P.append(dict(name="секция турбины (сварка)", kind="rotor_section", method=c["m_tw"], material=c["mat_tw"], D=tg["Dtip"],
                      L=0.45 * g["L_rot"], mass_kg=m["tw"] + 0.5 * m["shaft"], rotating=True))
    else:
        for k, s in enumerate(st):
            kind = "wheel_ti" if s["mat"] == "Ti64" else "wheel_al"
            P.append(dict(name=f"колесо компрессора {k + 1}", kind=kind, method=c["m_cw"], material=s["mat"], D=s["D"],
                          L=0.35 * s["D"], mass_kg=(0.0825 if s["type"] != "О" else 0.03) * {"Ti64": 4430.0, "Al": 2800.0}[s["mat"]] * s["D"] ** 3,
                          rotating=True))
        P.append(dict(name="вал", kind="shaft", method="поковка", material="steel", D=0.2 * g["Dc"], L=g["L_rot"], mass_kg=m["shaft"],
                      rotating=True))
        if c["turb"] == "Р":
            kind = "turb_radial"
            P.append(dict(name="колесо турбины", kind=kind, method=c["m_tw"], material=c["mat_tw"], D=tg["Dtip"], L=0.35 * tg["Dtip"],
                          mass_kg=m["tw"] / c["nt"], rotating=True, n=c["nt"]))
        else:
            for k in range(c["nt"]):
                if blisk[k]:
                    P.append(dict(name=f"блиск турбины {k + 1}", kind="turb_axial_blisk", method=c["m_tw"], material=c["mat_tw"],
                                  D=tg["Dtip"], L=0.15 * tg["D"], mass_kg=m["tw"] / c["nt"], rotating=True))
                else:
                    P.append(dict(name=f"диск турбины {k + 1}", kind="disc", method="поковка", material="IN718", D=tg["D"] - tg["h"],
                                  L=0.15 * tg["D"], mass_kg=0.6 * m["tw"] / c["nt"], rotating=True))
                    P.append(dict(name=f"лопатки турбины {k + 1}", kind="blade", method="литьё", material=c["mat_tw"], D=tg["h"],
                                  L=1.5 * tg["h"], mass_kg=0.4 * m["tw"] / c["nt"] / c["zt"], rotating=True, n=c["zt"]))
    # статор горячей части
    ngv_ring = []
    if c["hot_integral"]:
        L = g["L_hot"] + (g["L_comb"] if c["comb"] == "кольцевая" and c["m_comb"] == "печать" else 0.0)
        P.append(dict(name="горячая часть целиком", kind="static_hot", method=c["m_hot"], material=c["mat_hot"], D=g["D_hot"], L=L,
                      mass_kg=m["hot"] + (m["comb"] if L > g["L_hot"] else 0.0), rotating=False))
    else:
        for k in range(c["nt"]):
            if c["m_hot"] == "печать":
                ring = bool(MK.printers_for(g["D_hot"] * 0.85, 0.15 * g["D_hot"], False))
            elif c["m_hot"] == "литьё":
                ring = c["cool"] != "K" and g["D_hot"] * 0.85 <= MK.CAST["static"]["max"]
            else:
                ring = False
            ngv_ring.append(ring)
        P.append(dict(name="сопловой аппарат", kind="static_hot", method=c["m_hot"], material=c["mat_hot"], D=g["D_hot"] * 0.85,
                      L=0.15 * g["D_hot"], mass_kg=0.02 * MK.MAT.get(c["mat_hot"], {}).get("rho", 8200.0) * tg["D"] ** 3, rotating=False,
                      n=c["nt"]))
        P.append(dict(name="корпус турбины", kind="static", method=c["m_case"] if c["m_case"] != "печать" else "литьё",
                      material="steel", D=g["D_hot"], L=g["L_hot"], mass_kg=m["hot"], rotating=False))
    if not (c["hot_integral"] and c["comb"] == "кольцевая" and c["m_comb"] == "печать"):
        mc = c["m_comb"] if c["comb"] not in ("RDC",) else "печать"
        mat = {"печать": "HastX", "керамика": "CMC"}.get(mc, "steel")
        P.append(dict(name="жаровая труба" if c["comb"] != "RDC" else "кольцевой канал RDC", kind="static_hot", method=mc, material=mat,
                      D=g["D_comb"], L=g["L_comb"], mass_kg=m["comb"], rotating=False))
    # холодные корпуса
    pieces = 1
    if c["m_case"] == "печать":
        if not MK.printers_for(g["D_case"], g["L_case"], False):
            pieces = max(2, math.ceil(max(g["D_case"], g["L_case"]) / MK.largest_print_dim()))
    P.append(dict(name="корпус компрессора", kind="static", method=c["m_case"], material="Ti64" if c["m_case"] == "печать" else "steel",
                  D=g["D_case"] / (pieces ** 0.5 if pieces > 1 else 1.0), L=g["L_case"], mass_kg=m["case"], rotating=False,
                  split_ok=True))
    if g["rec"]:
        P.append(dict(name="модуль рекуператора", kind="static", method="печать", material="IN625", D=g["rec"]["side"], L=g["rec"]["side"],
                      mass_kg=g["rec"]["mass_kg"] / g["rec"]["modules"], rotating=False, n=g["rec"]["modules"]))
    pt_blisk = []
    if g["pt"]:
        mpt = c["m_tw"] if c["m_tw"] != "керамика" else "литьё"
        matpt = c["mat_tw"] if mpt == "печать" else "K418"
        pt_blisk = [_blisk(mpt, matpt, g["pt"]["Dtip"], g["pt"]["D"])] * c["npt"]
        for k in range(c["npt"]):
            if pt_blisk[k]:
                P.append(dict(name=f"колесо СТ {k + 1}", kind="turb_axial_blisk", method=mpt, material=matpt, D=g["pt"]["Dtip"],
                              L=0.15 * g["pt"]["D"], mass_kg=0.6 * m["pt"] / c["npt"], rotating=True))
            else:
                P.append(dict(name=f"диск СТ {k + 1}", kind="disc", method="поковка", material="IN718", D=g["pt"]["D"] - g["pt"]["h"],
                              L=0.15 * g["pt"]["D"], mass_kg=0.4 * m["pt"] / c["npt"], rotating=True))
                P.append(dict(name=f"лопатки СТ {k + 1}", kind="blade", method="литьё", material="K418", D=g["pt"]["h"], L=1.5 * g["pt"]["h"],
                              mass_kg=0.2 * m["pt"] / c["npt"] / c["zt"], rotating=True, n=c["zt"]))
    for p in P:
        p["n"] = p.get("n", 1) * (c["cores"] if p["name"] not in ("модуль рекуператора",) and not p["name"].startswith(("колесо СТ", "диск СТ", "лопатки СТ")) else 1)
    return P, blisk, pt_blisk, ngv_ring, pieces


def model_inputs(c, tgt, perf, g):
    """ведомость для parts_model, модули для service_model, интерфейсы для plugplay_model, детали для make_model"""
    keyed = c["conn"] == "быстроразъём"
    comb_kind = c["comb"]
    bom = dict(cores=c["cores"], rotors=c["rotors"], comp=stage_types(c), rotor_join=c["rotor_join"], bear=c["bear"], conn=c["conn"],
               inlet=c["inlet"],
               turb=dict(type=c["turb"], stages=[dict(blisk=g["blisk"][k], zt=c["zt"], zn=c["zn"],
                                                      ngv_ring=(g["ngv_ring"][k] if g["ngv_ring"] else True)) for k in range(c["nt"])]),
               pt=(dict(stages=[dict(blisk=g["pt_blisk"][k], zt=c["zt"], zn=c["zn"], ngv_ring=True) for k in range(c["npt"])])
                   if g["pt"] else None),
               case=dict(method=c["m_case"], D=g["D_case"], pieces=g["case_pieces"]),
               comb=dict(kind=comb_kind, method=c["m_comb"], n_inj=g["n_inj"]),
               hot_integral=bool(c["hot_integral"]), ngv_cooled=c["cool"] in ("N", "K", "P"),
               rotor_cooled=MK.rotor_cooled(c["cool"], c["m_tw"]),
               recup=dict(present=bool(g["rec"]), modules=g["rec"]["modules"] if g["rec"] else 0, integral=bool(c["recup_integral"])),
               new_starters=(c["cores"] if c["cores"] > 1 else (1 if c["rotors"] > tgt["rotors0"] else 0)))
    m = g["mass"]

    def J(D):
        return dict(type="хомут") if (keyed and D <= PM.CLAMP_D_MAX) else dict(type="болты", n_bolts=math.ceil(math.pi * D / PM.BOLT_PITCH))

    oil = c["bear"] == "масло"
    hot_int = bool(c["hot_integral"])
    comb_in_hot = hot_int and c["comb"] == "кольцевая" and c["m_comb"] == "печать"
    mods = []
    if c["cartridge"]:
        rotor_mod = dict(name="картридж ротора", mass_kg=1.3 * m["rotor"], joints=[J(0.5 * g["D_case"])] * 2, pipes=4 if oil else 0,
                         pipes_keyed=keyed, cables=6 if oil else 3, cables_keyed=keyed, align=False, balance=False, tools=[],
                         frequent=False, access=["горячая часть"])
    else:
        rotor_mod = dict(name="ротор", mass_kg=m["rotor"], joints=[dict(type="болты", n_bolts=math.ceil(2 * g["L_case"] / PM.BOLT_PITCH)),
                                                                   J(0.3 * g["Dc"])],
                         pipes=4 if oil else 0, pipes_keyed=keyed, cables=6 if oil else 3, cables_keyed=keyed, align=True,
                         balance=c["rotor_join"] == "сборный",
                         tools=(["гидравлический натяжитель стяжного болта"] if c["rotor_join"] == "сборный" else []),
                         frequent=False, access=["горячая часть"] + ([] if comb_in_hot or c["comb"] == "выносная" else ["камера"]))
    hot_mod = dict(name="горячая часть", mass_kg=m["hot"] + (m["comb"] if comb_in_hot else 0.0), joints=[J(g["D_hot"])] * 2, pipes=0,
                   cables=8, cables_keyed=keyed, frequent=True, access=[])
    comb_mod = None
    if not comb_in_hot:
        if c["comb"] in ("кольцевая", "RDC"):
            integ = c["m_comb"] == "печать"
            comb_mod = dict(name="камера", mass_kg=m["comb"], joints=[J(g["D_comb"])], pipes=1 if integ else g["n_inj"], pipes_keyed=keyed,
                            cables=2, cables_keyed=keyed, frequent=True, access=["горячая часть"])
        else:
            comb_mod = dict(name="камера", mass_kg=m["comb"], joints=[J(1.2 * g["D_comb"])] * 2, pipes=1, pipes_keyed=keyed,
                            cables=3 if c["comb"] == "ВР" else 2, cables_keyed=keyed, frequent=True, access=[])
    if c["cores"] > 1:
        core = dict(name="ядро", mass_kg=m["rotor"] + m["case"] + m["hot"] + m["comb"], joints=[J(g["D_case"]), J(g["D_hot"])],
                    pipes=1 + (4 if oil else 0), pipes_keyed=keyed, cables=16 if not keyed else 1, cables_keyed=keyed, align=False,
                    balance=False, tools=[], frequent=True, access=[])
        mods.append(core)
    else:
        mods += [rotor_mod, hot_mod] + ([comb_mod] if comb_mod else [])
    if g["pt"]:
        mods.append(dict(name="силовая турбина", mass_kg=m["pt"], joints=[J(1.3 * g["pt"]["Dtip"])] * 2, pipes=4 if oil else 0,
                         pipes_keyed=keyed, cables=4, cables_keyed=keyed, align=not c["cartridge"], frequent=False, access=[]))
    if g["rec"]:
        mods.append(dict(name="модуль рекуператора", mass_kg=g["rec"]["mass_kg"] / g["rec"]["modules"],
                         joints=[J(min(g["rec"]["side"], 1.5))] * (2 if c["recup_integral"] else 4), pipes=0, cables=2, cables_keyed=keyed,
                         frequent=False, access=[]))
    if not oil:
        mods.append(dict(name="контроллер магнитных подшипников", mass_kg=30.0, joints=[], pipes=0, cables=2, cables_keyed=True,
                         frequent=False, access=[], tools=["программа наладки магнитных подшипников"]))
    mods.append(dict(name="блок датчиков и жгут", mass_kg=10.0, joints=[], pipes=0, cables=1, cables_keyed=keyed, frequent=False, access=[]))
    sv = dict(modules=mods, amb=not oil, single_connector=keyed)
    pp = dict(n_out=(g["n"] if c["rotors"] == 1 else tgt["n_pkg"]), n_pkg=tgt["n_pkg"], G=perf["G"], G0=tgt["G0"], EGT=perf["EGT"],
              EGT0=tgt["EGT0"], D_env=g["D_env"], pkg_W=tgt["pkg_W"], pkg_H=tgt["pkg_H"], rotors=c["rotors"], rotors0=tgt["rotors0"],
              cores=c["cores"], inlet=c["inlet"], bear=c["bear"], recup=bool(g["rec"]), recup_integral=bool(c["recup_integral"]),
              comb=c["comb"], pi=c["pi"], pi0=tgt["pi0"], eta=perf["eta_th"] * 100 * 0.985 * 0.972, eta0=tgt["eta0"],
              module_kg_max=max(x["mass_kg"] for x in mods), engine_kg0=tgt["engine_kg0"])
    mk = dict(parts=g["parts_mk"], techs=techs(c))
    return bom, sv, pp, mk


def techs(c):
    T = {"традиция"}
    if c["rotor_join"] == "заодно":
        T.add("печать_ротор_монолит" if c["m_tw"] == "печать" else "литьё_ротор_монолит")
    elif c["rotor_join"] == "сварка":
        T.add("сварка_трением_отливок")
        if c["m_tw"] == "печать":
            T.add("печать_колесо_турбины")
    else:
        if c["m_cw"] == "печать":
            T.add("печать_колесо_компрессора")
        if c["m_tw"] == "печать":
            T.add("печать_колесо_турбины")
        if c["m_tw"] == "керамика":
            T.add("керамика_ротор")
    if c["m_hot"] == "печать" or c["m_comb"] == "печать":
        T.add("печать_статика")
    if c["hot_integral"] and c["m_hot"] == "печать":
        T.add("печать_горячая_часть_целиком")
    if c["m_hot"] == "керамика" or c["m_comb"] == "керамика":
        T.add("керамика_статика")
    if c["cool"] == "P":
        T.add("печатное_охлаждение")
    if c["recup"] == "печатный":
        T.add("печать_рекуператор")
    if c["bear"] == "магнит":
        T.add("магнитные_подшипники")
    if c["cores"] > 1:
        T.add("кластер_ядер")
    if c["comp"] == "1Ц":
        T.add("одно_колесо_до_7" if c["pi"] <= 7.0 else "одно_колесо_выше_7")
    if c["comb"] == "RDC":
        T.add("rdc")
    if c["comb"] == "ВР":
        T.add("волновой_ротор")
    if c["mat_tw"] == "SC":
        T.add("лопатки_SC")
    return sorted(T)


# =========================================================================== правила допустимости
AX_PI_MAX, D_PI_MAX = 1.8, 5.0          # степень сжатия осевой и диагональной ступени (оценка)
PSI_MAX = 2.5                            # gtu_explore.turb_axial, frame3
NU_RANGE = (0.50, 0.90)                  # U/C0 радиальной турбины (оценка; оптимум 0,70 - gtu_explore)
NS_RANGE = (0.25, 1.00)                  # быстроходность радиальной турбины ω·√Q/Δh^0,75 (оценка)
HTR_MIN_T, HTR_MIN_AX = 0.40, 0.35       # отношение втулки к периферии: последняя ступень турбины (оценка; frame3 - 0,54), осевая ступень компрессора (оценка)
T_REC_MAX = 700.0                        # LIM T_rec_max


def _rule(code, ok, val=None, lim=None, text=""):
    f = (lambda x: round(x, 4) if isinstance(x, float) else x)
    return dict(code=code, ok=bool(ok), val=f(val), lim=f(lim), text=text)


def check(c, tgt, perf, g):
    R = []
    if not perf["ok"]:
        return [_rule("perf", False, text=perf["fail"])]
    R.append(_rule("perf", True))
    for k, s in enumerate(g["stages"]):
        method = c["m_cw"]
        Ul, ub = MK.u_limit_comp(s["mat"], method)
        R.append(_rule(f"U_comp{k + 1}", s["U"] <= Ul, s["U"], Ul, f"окружная скорость ступени {k + 1} ({s['mat']}, {method}): {ub}"))
        if s["mat"] == "Al":
            R.append(_rule(f"T_al{k + 1}", s["T2"] <= MK.MAT["Al"]["T_max"], s["T2"], MK.MAT["Al"]["T_max"], "LIM T_al: алюминий за колесом <= 200 C"))
        if s["type"] == "О":
            R.append(_rule(f"pi_ax{k + 1}", s["pi"] <= AX_PI_MAX, s["pi"], AX_PI_MAX, "степень сжатия осевой ступени (оценка)"))
            R.append(_rule(f"htr_ax{k + 1}", s["htr"] >= HTR_MIN_AX, s["htr"], HTR_MIN_AX, "втулка/периферия осевой ступени компрессора (оценка)"))
        else:
            if s["type"] == "Д":
                R.append(_rule(f"pi_d{k + 1}", s["pi"] <= D_PI_MAX, s["pi"], D_PI_MAX, "степень сжатия диагональной ступени (оценка)"))
            lo, hi = PHI_RANGE[s["type"]]
            fixed = c["rotors"] == 1
            R.append(_rule(f"phi{k + 1}", lo <= s["phi"] <= hi, s["phi"], [lo, hi],
                           ("частота вала = частота пакета: коэффициент расхода колеса при этой частоте" if fixed else
                            "коэффициент расхода ступени на частоте газогенератора") + " (диапазон - оценка)"))
    R.append(_rule("rpm_pkg", True if c["rotors"] == 2 else abs(g["n"] - tgt["n_pkg"]) < 1e-6, g["n"] if c["rotors"] == 1 else tgt["n_pkg"],
                   tgt["n_pkg"], "одновальный: вал двигателя = вход редуктора пакета; двухвальный: силовая турбина на частоте пакета"))
    tl, tb = MK.t41_limit(c["turb"], c["cool"], c["mat_tw"], c["m_tw"])
    R.append(_rule("T41_class", c["t41"] <= tl + 1e-6, c["t41"], tl, f"T41 по классу охлаждения {c['cool']} и сплаву {c['mat_tw']}: {tb}"))
    tg = g["tg"]
    R.append(_rule("U_turb", tg["U"] <= tg["Ulim"] + 1e-6, tg["U"], tg["Ulim"], "окружная скорость турбины по сплаву и температуре"))
    if c["turb"] == "О":
        R.append(_rule("turb_psi", tg["psi"] <= PSI_MAX, tg["psi"], PSI_MAX, "нагрузка ступени осевой турбины (gtu_explore, frame3)"))
        R.append(_rule("turb_htr", tg["htr"] >= HTR_MIN_T, tg["htr"], HTR_MIN_T, "втулка/периферия последней ступени турбины (оценка; frame3 - 0,54)"))
    else:
        R.append(_rule("turb_nu", NU_RANGE[0] <= tg["nu"] <= NU_RANGE[1], tg["nu"], list(NU_RANGE), "U/C0 ступени радиальной турбины (оценка)"))
        R.append(_rule("turb_ns", NS_RANGE[0] <= tg["ns"] <= NS_RANGE[1], tg["ns"], list(NS_RANGE), "быстроходность радиальной турбины (оценка)"))
    if g["pt"]:
        R.append(_rule("pt_psi", g["pt"]["psi"] <= PSI_MAX, g["pt"]["psi"], PSI_MAX, "нагрузка ступени силовой турбины"))
        R.append(_rule("pt_htr", g["pt"]["htr"] >= HTR_MIN_T, g["pt"]["htr"], HTR_MIN_T, "втулка/периферия последней ступени СТ (оценка)"))
    if c["rotor_join"] in ("заодно", "сварка"):
        R.append(_rule("mono_blisk", all(g["blisk"]), sum(g["blisk"]), c["nt"], "ротор заодно или сварной требует цельных колёс турбины"))
    for p in g["parts_mk"]:
        a = MK.assess_part(p)
        if a["fit"]:
            continue
        hard = p["rotating"] or (p["name"] == "горячая часть целиком")
        R.append(_rule("size:" + p["name"], not hard, max(p["D"], p["L"]) if p["method"] == "печать" else p["D"], a["limit"],
                       f"{p['method']}: {a['basis']}" + ("" if hard else f" - делится на {a['pieces']} части")))
    if g["rec"]:
        R.append(_rule("recup_T", perf["T5"] <= T_REC_MAX, perf["T5"], T_REC_MAX, "LIM T_rec_max: газ на входе в рекуператор"))
    R.append(_rule("envelope", g["D_env"] <= min(tgt["pkg_W"], tgt["pkg_H"]), g["D_env"], min(tgt["pkg_W"], tgt["pkg_H"]),
                   "габарит двигателя (корпус + 0,5 м обвязки, gtu_family) против ширины и высоты пакета"))
    R.append(_rule("soft_G", perf["G"] <= tgt["G0"] * (1 + PP.SOFT["alpha"]), perf["G"], tgt["G0"] * (1 + PP.SOFT["alpha"]),
                   "мягкий допуск: расход до +20 % (решение владельца 14.11.2026)"))
    R.append(_rule("soft_EGT", perf["EGT"] <= tgt["EGT0"] + PP.SOFT["dT"], perf["EGT"], tgt["EGT0"] + PP.SOFT["dT"],
                   "мягкий допуск: выхлоп до +100 °C (решение владельца 14.11.2026)"))
    return R


def violation(bad):
    """сумма относительных превышений по непройденным правилам (0 - допустима); для ранжирования «почти допустимых»
    на следующих уровнях. Правило без числа (perf) - 1."""
    v = 0.0
    for r in bad:
        val, lim = r.get("val"), r.get("lim")
        if isinstance(lim, list) and isinstance(val, (int, float)):
            lo, hi = lim
            v += (lo - val) / abs(lo) if val < lo else (val - hi) / abs(hi)
        elif isinstance(val, (int, float)) and isinstance(lim, (int, float)) and lim:
            v += abs(val - lim) / abs(lim)
        else:
            v += 1.0
    return round(v, 4)


# =========================================================================== одна концепция целиком
def evaluate(i, n, kv="KV-3", seed=SEED, detail=False):
    c = concept(i, n, kv, seed)
    tgt = target(kv)
    p1 = quick_perf(c, tgt)
    g, p = None, None
    if p1["ok"]:
        try:
            g1 = derive(c, tgt, p1)
            p = quick_perf(c, tgt, g1["d_eta"])
            if p["ok"]:
                g = derive(c, tgt, p)
            else:
                p1 = p
        except (ValueError, ZeroDivisionError, OverflowError) as e:
            p1 = dict(ok=False, fail=f"геометрия не строится: {e}")
    perf = p if g else p1
    rules = check(c, tgt, perf, g) if g else [_rule("perf", False, text=perf.get("fail", "нет решения"))]
    bad = [r for r in rules if not r["ok"]]
    rec = dict(idx=i, kv=kv, arch=c["arch"], feasible=not bad, reject=[r["code"] for r in bad], violation=violation(bad), concept=c,
               rules=rules if detail else [r for r in rules if not r["ok"]],
               rules_passed=[r["code"] for r in rules if r["ok"]])
    if g:
        rec["perf"] = dict(eta_th_shaft=round(perf["eta_th"] * 100, 2), eta_el=round(perf["eta_th"] * 100 * 0.985 * 0.972, 2),
                           w_kJkg=round(perf["w"], 1), G=round(perf["G"], 3), dG_pct=round((perf["G"] / tgt["G0"] - 1) * 100, 1),
                           EGT=round(perf["EGT"], 1), dEGT=round(perf["EGT"] - tgt["EGT0"], 1), T2=round(perf["T2"], 1),
                           T5=round(perf["T5"], 1), eta_t=round(perf["eta_gg"], 4),
                           eta_pt=round(perf["eta_pt"], 4) if perf["eta_pt"] else None)
        rec["geo"] = dict(n=round(g["n"]), D2=[round(s["D"], 4) for s in g["stages"]], U2=[round(s["U"], 1) for s in g["stages"]],
                          phi=[round(s.get("phi", float("nan")), 4) if "phi" in s else None for s in g["stages"]],
                          b2_mm=[round(s["b2"] * 1e3, 1) if "b2" in s else None for s in g["stages"]],
                          eta_c=[round(s["eta"], 4) for s in g["stages"]],
                          turb=dict(type=c["turb"], nt=c["nt"], D=round(g["tg"]["D"], 4), U=round(g["tg"]["U"], 1),
                                    h_last=round(g["tg"]["h"], 4) if g["tg"]["h"] else None,
                                    psi=round(g["tg"]["psi"], 3) if g["tg"].get("psi") else None,
                                    nu=round(g["tg"]["nu"], 3) if g["tg"].get("nu") else None,
                                    ns=round(g["tg"]["ns"], 3) if g["tg"].get("ns") else None),
                          pt=(dict(D=round(g["pt"]["D"], 4), U=round(g["pt"]["U"], 1), psi=round(g["pt"]["psi"], 3), h_last=round(g["pt"]["h"], 4))
                              if g["pt"] else None),
                          D_env=round(g["D_env"], 3), L_rot=round(g["L_rot"], 3),
                          rotor_kg=round(g["mass"]["rotor"], 1),
                          rec=({k: round(v, 3) for k, v in g["rec"].items()} if g["rec"] else None))
        bom, sv, pp, mk = model_inputs(c, tgt, perf, g)
        pr, sr, ppr, mr = PM.count(bom), SM.score(sv), PP.score(pp), MK.score(mk)
        rec["scores"] = dict(parts=pr["score"], plugplay=ppr["score"], service=sr["score"], make=mr["score"])
        rec["parts"] = {k: pr[k] for k in ("parts", "assemblies", "moving", "cooled", "fasteners", "welds", "unique", "n_eff")}
        rec["service"] = {k: sr[k] for k in ("heaviest", "steps", "wrongable", "tools")}
        rec["plugplay"] = dict(blocked=ppr["blocked"], no_data=ppr["no_data"], states={k: v["state"] for k, v in ppr["items"].items()})
        rec["make"] = {k: mr[k] for k in ("trl_min", "weeks", "cost_u", "s_fit", "s_proc")}
        if detail:
            rec["detail"] = dict(parts=pr["by_unit"], service=sr, plugplay=ppr["items"], make=mr, bom=bom)
    return rec


# =========================================================================== часть выборки (AWS) и отладка
def software():
    import cantera
    return f"python {platform.python_version()}, numpy {np.__version__}, scipy {scipy.__version__}, cantera {cantera.__version__}"


def _eval_many(args):
    idx, n, kv, seed, detail = args
    return [evaluate(i, n, kv, seed, detail) for i in idx]


def run_shard(kv, n, shard, shards, out, seed=SEED, detail=False, workers=1):
    """часть shard из shards: концепции i % shards == shard; workers - процессы внутри части (на AWS - nproc)"""
    t0 = time.time()
    idx = shard_indices(n, shard, shards)
    if workers > 1 and len(idx) > 1:
        from concurrent.futures import ProcessPoolExecutor
        chunks = [idx[k::workers] for k in range(workers)]
        with ProcessPoolExecutor(workers) as ex:
            got = {r["idx"]: r for part in ex.map(_eval_many, [(ch, n, kv, seed, detail) for ch in chunks if ch]) for r in part}
        recs = [got[i] for i in idx]
    else:
        recs = _eval_many((idx, n, kv, seed, detail))
    doc = dict(method="воронка, уровень 0, часть пространства: концепции по индексу из последовательности Соболя, быстрый цикл "
                      "на таблицах gtu_explore (Cantera), геометрия и правила допустимости funnel_space, баллы parts/plugplay/"
                      "service/make; всё - оценка (README узла)",
               software=software(), date=DATE, kv=kv, target=target(kv), n=n, seed=seed, shard=shard, shards=shards,
               fingerprint=fingerprint(n, seed), indices=idx, workers=workers, elapsed_s=round(time.time() - t0, 2),
               inputs=inputs(), records=recs)
    os.makedirs(out, exist_ok=True)
    path = os.path.join(out, f"funnel_space_{kv}_n{n}_s{shard:03d}of{shards:03d}.json")
    with open(path, "w", encoding="utf-8") as fh:
        json.dump(doc, fh, ensure_ascii=False, separators=(",", ":"), default=_jd)
    return path, doc


def _jd(o):
    """numpy-типы -> обычные числа для json"""
    if isinstance(o, np.generic):
        return o.item()
    if isinstance(o, (set, tuple)):
        return list(o)
    raise TypeError(type(o))


def _sha(p):
    return hashlib.sha256(open(p, "rb").read()).hexdigest()


def inputs():
    files = ["gtu_family_results.json", "gtu_lineup_results.json", "gtu_explore.py", "gtu_cycle.py", "frame3/frame3_point_results.json",
             "fleet/gtu_fleet_taurus60.json", "fleet/gtu_fleet_saturn20.json", "fleet/gtu_fleet_titan130.json", "platform/fea/materials.json"]
    fields = {"gtu_family_results.json": "slots[model, variant=генератор]: P_shaft_kW, G_air, EGT0, eta0, rpm0, PR0, pkg, shafts0",
              "gtu_lineup_results.json": "lineup[name=KV-1/KV-3/KV-6].members - проверка принадлежности расчётных точек",
              "gtu_explore.py": "таблицы h, s0, T_h, T_s, far_for, bfrac, Rg, eta_axial; LIM - через make_model",
              "gtu_cycle.py": "сверка _core с cycle() в --selfcheck",
              "frame3/frame3_point_results.json": "stage1/stage2: φ 0,040 и 0,028, поправки КПД -> коэффициент PHI_PENALTY; turbine: htr последней ступени",
              "fleet/gtu_fleet_taurus60.json": "records[].engine_kg (MD 3290 кг) - масса исходного двигателя",
              "fleet/gtu_fleet_saturn20.json": "engine_kg - нет", "fleet/gtu_fleet_titan130.json": "engine_kg - нет",
              "platform/fea/materials.json": "ti64 - предел данных 427 C"}
    return [dict(file="routes/gtu/" + f, field=fields.get(f, ""), sha256=_sha(os.path.join(GTU, f))) for f in files]


def _stats(v):
    v = sorted(x for x in v if x is not None)
    if not v:
        return None
    q = lambda p: v[min(len(v) - 1, int(p * (len(v) - 1) + 0.5))]
    return dict(n=len(v), min=round(v[0], 3), p50=round(q(0.5), 3), max=round(v[-1], 3))


def summarize(recs):
    from collections import Counter
    out = dict(n=len(recs), feasible=sum(r["feasible"] for r in recs))
    rej = Counter()
    for r in recs:
        for code in r["reject"]:
            rej[code.split(":")[0] if code.startswith("size:") else code.rstrip("0123456789")] += 1
    out["reject_reasons"] = dict(rej.most_common())
    by = {}
    for a in ARCH_KEYS:
        rs = [r for r in recs if r["arch"] == a]
        ok = [r for r in rs if r["feasible"]]
        first = Counter(code for r in rs for code in r["reject"][:1])
        by[a] = dict(n=len(rs), feasible=len(ok), first_reject=dict(first.most_common(4)))
        if ok:
            for k in ("parts", "plugplay", "service", "make"):
                by[a][k] = _stats([r["scores"][k] for r in ok])
            by[a]["eta_el"] = _stats([r["perf"]["eta_el"] for r in ok])
            by[a]["n_parts"] = _stats([r["parts"]["parts"] for r in ok])
    out["by_arch"] = by
    return out


def cross_check():
    """_core (таблицы gtu_explore) против gtu_cycle.cycle (Cantera напрямую) при одинаковых КПД и долях воздуха"""
    _E()
    cwd = os.getcwd()
    from gtu_cycle import cycle
    os.chdir(cwd)
    rows = []
    for pi, T41, cool in ((6.75, 1000.0, "U"), (13.5, 1000.0, "U"), (11.0, 1050.0, "N"), (16.0, 1150.0, "K")):
        xn, xr = MK.COOL_AIR[cool]
        c = dict(cool=cool, m_tw="литьё", cores=1, recup="нет", comb="кольцевая", pgain=0.0, rotors=2, t41=T41, eps_rec=0.0)
        comp = _comp([(pi, 0.86)])
        r = _core(c, comp, 0.865, 0.885)
        g = cycle(pi, T41, eta_pc=0.86, eta_pt_gg=0.865, eta_pt_pt=0.885, x_ngv=xn, x_rot=xr, x_leak=0.015)
        rows.append(dict(pi=pi, T41=T41, cool=cool, eta_core=round(r["eta"] * 100, 2), eta_cycle=round(g["eff"] * 100, 2),
                         d_eta_pp=round((r["eta"] - g["eff"]) * 100, 2), w_core=round(r["wnet"] / 1e3, 1), w_cycle=round(g["sp_work_kJkg"], 1),
                         EGT_core=round(r["Tstack"] - 273.15, 1), EGT_cycle=round(g["EGT_C"], 1)))
    return rows


SELF_N = {"KV-3": 1024, "KV-1": 256, "KV-6": 512}
SELF_WORKERS = 4


def selfcheck(out_path=RESULTS):
    t0 = time.time()
    lineup = json.load(open(os.path.join(GTU, "gtu_lineup_results.json"), encoding="utf-8"))["lineup"]
    member = {kv: any(m["model"] == KV[kv][0] and m["variant"] == KV[kv][1] for f in lineup if f["name"] == kv for m in f["members"])
              for kv in KV}
    res = dict(title="Воронка KV, уровень 0, часть 1: пространство поиска и модели оценки - отладка на малой выборке",
               method="скремблированная последовательность Соболя (scipy.stats.qmc, зерно SEED) по осям AXES; архетипы ARCH - шаблоны "
                      "дискретных осей; быстрый цикл уровня 0 на таблицах свойств gtu_explore (Cantera GRI-Mech 3.0, метан, ISO) с КПД "
                      "турбин по нагрузке ступеней; геометрия по окружным скоростям и коэффициентам расхода; правила допустимости; "
                      "баллы parts_model, service_model, plugplay_model, make_model (формулы - README узла). Всё - оценка.",
               software=software(), date=DATE, seed=SEED, status="отладка части 1 на 4 ядрах; объёмы воронки (20 000 / 5 000) считаются на AWS",
               inputs=inputs(), targets={kv: target(kv) for kv in KV}, targets_in_lineup=member,
               axes=dict(categorical={k: v for k, v in CAT.items()}, continuous=CONT, order=AXES),
               archetypes={k: dict(name=v["name"], allowed=v["allowed"], note=v.get("note")) for k, v in ARCH.items()},
               defaults=_DEF,
               rules=dict(AX_PI_MAX=AX_PI_MAX, D_PI_MAX=D_PI_MAX, PSI_MAX=PSI_MAX, NU_RANGE=NU_RANGE, NS_RANGE=NS_RANGE,
                          HTR_MIN_T=HTR_MIN_T, HTR_MIN_AX=HTR_MIN_AX, T_REC_MAX=T_REC_MAX, PHI_RANGE=PHI_RANGE, PHI_OPT=PHI_OPT,
                          PHI_PENALTY=PHI_PENALTY, T41_MIN=T41_MIN, soft=PP.SOFT),
               models=dict(parts=dict(N_HALF=round(PM.N_HALF, 2), K=dict(assembly=PM.K_ASSEMBLY, fastener=PM.K_FASTENER, moving=PM.K_MOVING,
                                                                            cooled=PM.K_COOLED), ref_bom=PM.REF_BOM),
                           service=dict(weights=SM.W, step=SM.STEP, m_one_person=SM.M_ONE_PERSON, m_zero=SM.M_ZERO),
                           plugplay=dict(weights=PP.WEIGHTS, states=PP.STATE, soft=PP.SOFT),
                           make=dict(weights=MK.W)),
               facts=MK.facts())
    # сверка цикла двумя способами
    res["cross_check_cycle"] = cross_check()
    # разбиение на части: KV-3, n=256, 4 части - объединение = все индексы, пересечений нет, запись не зависит от части
    from concurrent.futures import ProcessPoolExecutor
    n, S = SELF_N["KV-3"], 8
    parts = [shard_indices(n, k, S) for k in range(S)]
    flat = sorted(i for p in parts for i in p)
    tt = time.time()
    with ProcessPoolExecutor(SELF_WORKERS) as ex:
        got = [r for part in ex.map(_eval_many, [(p, n, "KV-3", SEED, False) for p in parts]) for r in part]
    dt = time.time() - tt
    recs = sorted(got, key=lambda r: r["idx"])
    by_order = [evaluate(i, n, "KV-3") for i in reversed(parts[3][:16])]
    same = all(json.dumps(a, sort_keys=True, default=float) == json.dumps(recs[a["idx"]], sort_keys=True, default=float) for a in by_order)
    res["sharding_check"] = dict(n=n, shards=S, union_ok=flat == list(range(n)), duplicates=len(flat) - len(set(flat)),
                                 record_independent_of_shard_and_order=same, fingerprint=fingerprint(n))
    res["timing"] = dict(seconds_per_concept_per_core=round(dt * SELF_WORKERS / n, 4), workers=SELF_WORKERS,
                         note="оценка на этой машине (4 ядра); на c7a.16xlarge - 64 vCPU")
    res["debug_KV-3"] = summarize(recs)
    near = sorted([r for r in recs if not r["feasible"]], key=lambda r: r["violation"])[:10]
    res["debug_KV-3"]["nearest_infeasible"] = [dict(idx=r["idx"], arch=r["arch"], violation=r["violation"], reject=r["reject"]) for r in near]
    ex = {}
    for a in ARCH_KEYS:
        ok = [r for r in recs if r["arch"] == a and r["feasible"]]
        if ok:
            best = max(ok, key=lambda r: r["scores"]["parts"])
            ex[a] = evaluate(best["idx"], n, "KV-3", detail=True)
            ex[a]["note"] = "пример: допустимая концепция архетипа с наибольшим баллом деталей в отладочной выборке"
    res["examples_KV-3"] = ex
    for kv in ("KV-1", "KV-6"):
        m = SELF_N[kv]
        with ProcessPoolExecutor(SELF_WORKERS) as ex:
            rk = [r for part in ex.map(_eval_many, [(shard_indices(m, k, S), m, kv, SEED, False) for k in range(S)]) for r in part]
        res[f"debug_{kv}"] = summarize(rk)
        near = sorted([r for r in rk if not r["feasible"]], key=lambda r: r["violation"])[:10]
        res[f"debug_{kv}"]["nearest_infeasible"] = [dict(idx=r["idx"], arch=r["arch"], violation=r["violation"], reject=r["reject"]) for r in near]
    res["elapsed_s"] = round(time.time() - t0, 1)
    with open(out_path, "w", encoding="utf-8") as fh:
        json.dump(res, fh, ensure_ascii=False, indent=1, default=_jd)
    return res


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--kv", default="KV-3", choices=sorted(KV))
    ap.add_argument("--n", type=int, default=20000)
    ap.add_argument("--shard", type=int, default=0)
    ap.add_argument("--shards", type=int, default=1)
    ap.add_argument("--seed", type=int, default=SEED)
    ap.add_argument("--out", default=os.environ.get("GTU_OUT", "."))
    ap.add_argument("--detail", action="store_true", help="полные разборы моделей в каждой записи")
    ap.add_argument("--workers", type=int, default=int(os.environ.get("GTU_WORKERS", "1")), help="процессов внутри части (на AWS - nproc)")
    ap.add_argument("--selfcheck", action="store_true", help="отладка на малой выборке -> gtu_funnel_space_results.json")
    a = ap.parse_args(argv)
    if a.selfcheck:
        r = selfcheck()
        print(json.dumps(dict(sharding=r["sharding_check"], timing=r["timing"], kv3=dict(n=r["debug_KV-3"]["n"], feasible=r["debug_KV-3"]["feasible"]),
                              cross=r["cross_check_cycle"]), ensure_ascii=False, indent=1, default=_jd))
        return 0
    path, doc = run_shard(a.kv, a.n, a.shard, a.shards, a.out, a.seed, a.detail, max(1, a.workers))
    print(f"{path}: {len(doc['records'])} концепций, допустимых {sum(r['feasible'] for r in doc['records'])}, {doc['elapsed_s']} с, "
          f"отпечаток {doc['fingerprint']}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
