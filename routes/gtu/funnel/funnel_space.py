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
Все результаты уровня 0 - «оценка»: КПД узлов, пределы и коэффициенты - принятые величины трека и свои оценки (README).
"""
import argparse, datetime, hashlib, inspect, itertools, json, math, os, platform, sys, time, warnings

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
CONT = {   # непрерывные оси, выборка v2: доля 0..1 Соболя -> величина ВНУТРИ допустимого окна; окна строит design()
    # по номинальному циклу концепции (КПД узлов при выбранных безразмерных параметрах, запасы M_W, M_T, M_GEO), зависимые
    # величины (расход, диаметры, частота газогенератора, окружная скорость турбины) выводятся, а не разыгрываются
    "t41": "T41, °C: окно [T41_lo, Tlim]; Tlim - по классу охлаждения, типу турбины и сплаву (make_model.t41_limit); T41_lo - наименьшая "
           "T41, при которой номинальный цикл даёт w >= w_need·(1+M_W) при выхлопе <= t0+100-M_T (рекуператор: газ на входе <= 700-M_T, "
           "алюминиевое колесо: T2 <= 200-M_T); равномерно в окне",
    "pi": "общая степень сжатия: окно при выбранной T41, где выполнены те же условия, внутри диапазона типа компрессора (1Ц 3-8, Д 2,5-6, "
          "2Ц 6-20, ОЦ 4-14 / 5-18 - оценка) и не выше pi_max_U (ступени в пределе окружной скорости при лучшем делении); равномерно по ln π",
    "split": "доля ln π на первую ступень (2Ц: 0,35-0,65) или на осевую часть (ОЦ: 0,12-0,60); окно - где каждая ступень проходит предел "
             "окружной скорости, π осевой ступени <= 1,8, а у вала с частотой пакета - коэффициент расхода в хорошем диапазоне GOOD",
    "mu1": "коэффициент работы w/U² первой (ОЦ - осевых) ступени: центробежная 0,60-0,80, диагональная 0,50-0,70, осевая 0,20-0,35 (оценка); "
           "окно сужают предел окружной скорости и (вал с частотой пакета) коэффициент расхода GOOD",
    "mu2": "то же для последней ступени (2Ц, ОЦ)",
    "phi1": "вал газогенератора: частота n выводится из окна, где все колёса и турбина в пределах размеров способа изготовления, коэффициенты "
            "расхода ступеней и быстроходность радиальной турбины в GOOD, AN² и втулка осевой турбины проходят при ca/U <= 0,8; "
            "phi1 = коэффициент расхода первой центробежной ступени при этой n (равномерно по ln n); у вала с частотой пакета - выводится",
    "kcm": "cm2/U2 на выходе колеса - ширина b2: 0,20-0,35",
    "zc": "число основных лопаток колеса компрессора: 7-15",
    "lt": "нагрузка турбины газогенератора: осевая - ψ = w/(nt·U²) в окне [1,2; 2,3] (не ниже ψ при предельной окружной скорости); "
          "радиальная - U/C0 в окне [0,62; 0,78] (не выше U_max/C0); окружная скорость выводится: U = √(w/(nt·ψ)) или (U/C0)·C0",
    "zt": "лопаток рабочего колеса турбины: радиальная 9-17, осевая 30-80 на ступень",
    "zn": "лопаток соплового аппарата: радиальная 12-24, осевая 20-60 на ступень",
    "kUpt": "окружная скорость силовой турбины в долях 380 м/с: окно [max(0,85; U_min/380); 1,00], U_min - при котором AN² и втулка "
            "последней ступени СТ проходят при ca/U <= 0,8",
    "phit": "ca/U на выходе осевой турбины газогенератора: окно [max(0,45; по втулке 0,40 и AN²); 0,80]; выше 0,5 - потеря с выходной "
            "скоростью (половина прироста c²/2, оценка)",
    "phipt": "ca/U на выходе силовой турбины: окно как у phit, по D и U силовой турбины на частоте пакета",
    "eps_rec": "эффективность рекуператора: 0,75-0,92; рекуператор в корпусе - не выше ε, при которой корпус входит в габарит пакета",
    "pgain": "прирост полного давления в камере: RDC 0-0,10, волновой ротор 0-0,15 (оценка; модели обещают +3-5 % КПД у RDC - по выдержке поиска)",
    "kD_case": "диаметр корпуса / наибольший D2: окно [1,7; min(2,1; по габариту пакета)] (gtu_family - 1,9)",
}
SAMPLER = ("v2 (14.11.2026): дискретные оси - по шаблону архетипа; комбинация, номинальный цикл которой не замыкает точку ни при каких "
           "T41 и π (repair), заменяется ближайшей замыкающей по осям REPAIR_AXES с записью причины; непрерывные оси - внутри окон "
           "допустимости (CONT), зависимые величины выводятся. Первая версия (v1: все оси независимо во всём диапазоне) дала 1,2 % "
           "допустимых в KV-3 и 0 % в KV-6 (прогоны AWS 20261114T231900Z-b7cc и 20261114T232200Z-7ce3)")
GOOD = dict(phi={"Ц": (0.055, 0.13), "Д": (0.08, 0.20)}, nu=(0.62, 0.78), psi=(0.9, 2.0), ns_t=(0.35, 0.85), phit=(0.45, 0.80))
GOOD_BASIS = ("хорошие диапазоны для выборки (оценка): φ центробежной ступени 0,055-0,13 (оптимум КПД около 0,085-0,11, frame3 при 0,028-0,040 "
              "теряет 1,2-3,2 п.); U/C0 радиальной турбины 0,62-0,78 вокруг оптимума 0,70 (gtu_explore.turb_radial); ψ осевой 0,9-2,0 "
              "(eta_axial 0,885 ниже 1,2 и 0,857 при 2,0; предел 2,5); быстроходность радиальной турбины 0,35-0,85 вокруг 0,60; правила check шире")
M_W, M_T, M_GEO = 0.015, 8.0, 0.03    # запасы окон: удельная работа +1,5 %, температуры -8 K, геометрия 3 % (оценка)
ETA_NOM = dict(comp=-0.012, ax=-0.016, rad=-0.042, pt=0.86)  # номинальные КПД окон: поправки к формулам ступеней (зазор, ns, политропный
#   эквивалент радиальной, выходная скорость) и КПД СТ; подобраны по медианам «итог - номинал» короткой пробы KV-3 N = 300 (осевая -0,011,
#   радиальная -0,027 к первым значениям -0,005 и -0,015) - оценка; остаток покрывают запасы M_W, M_T
SPLIT_RANGE = {"2Ц": (0.35, 0.65), "ОЦ": (0.12, 0.60)}
N_PI_GRID, N_T41_GRID, N_SPLIT_GRID = 10, 6, 21
REPAIR_AXES = ("mat_tw", "mat_cw", "m_tw", "cool", "comp", "turb", "recup", "rotors", "rotor_join")
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
        comp=["1Ц", "2Ц"], turb=["Р", "О"], cool=["U", "N", "K"], rotor_join=["заодно", "сварка"], m_tw=["литьё"], mat_tw=["K418", "CM247LC"],
        m_hot=["литьё"], m_comb=["литьё", "обработка"], m_case=["литьё"], hot_integral=[0, 1]),
        note="v2: добавлена осевая турбина цельнолитыми колёсами (до 0,50 м, LIM D_axcast); с одной радиальной (литьё, T41 <= 1050 °C) "
             "архетип не замыкает KV-3 - closure_table в результатах"),
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


def _discrete(u, arch, ov=None):
    """дискретные оси по шаблону архетипа; ov - замены отдельных осей (repair), применяются, если значение допустимо"""
    ov = ov or {}
    T = ARCH[arch]["allowed"]

    def cat(name, vals=None):
        v = list(vals if vals is not None else T.get(name, _DEF[name]))
        if name in ov and ov[name] in v:
            return ov[name]
        return _pick(u[name], v)

    c = dict(arch=arch)
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
    c["mat_tw"] = cat("mat_tw", _mat_tw_list(c, T))
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
    return c


def _mat_tw_list(c, T):
    """сплавы колеса турбины: монокристалл - только осевая и только сборный ротор (лопатки в замках не бывают цельным колесом,
    правило mono_blisk), нитрид кремния - только радиальная"""
    tw = [m for m in MAT_TW[c["m_tw"]] if not (m == "SC" and (c["turb"] != "О" or c["rotor_join"] in ("заодно", "сварка")))
          and not (m == "Si3N4" and c["turb"] != "Р")]
    if "mat_tw" in T:
        tw = [m for m in tw if m in T["mat_tw"]] or tw
    return tw


def _axis_vals(c, arch, axis):
    T = ARCH[arch]["allowed"]
    if axis == "mat_tw":
        return _mat_tw_list(c, T)
    if axis == "mat_cw":
        return [c["mat_cw"]] if c["rotor_join"] in ("заодно", "сварка") else [m for m in T.get("mat_cw", _DEF["mat_cw"]) if m != "Ni"]
    if axis == "rotors" and c["cores"] > 1:
        return [2]
    return list(T.get(axis, _DEF[axis]))


SIG_KEYS = ("kv", "rotors", "cores", "comp", "n_ax", "turb", "cool", "m_tw", "mat_tw", "recup", "comb", "rotor_join", "m_cw", "mat_cw")


def _sig(c):
    return tuple((k, c[k]) for k in SIG_KEYS)


def w_need(tgt):
    """наименьшая удельная работа, при которой расход укладывается в мягкий допуск: P_вал / ((1 + α)·G0), Дж/кг"""
    return tgt["P_shaft_kW"] * 1e3 / (tgt["G0"] * (1 + PP.SOFT["alpha"]))


def egt_max(tgt):
    return tgt["EGT0"] + PP.SOFT["dT"]


def pi_range(c):
    return PI_RANGE[(c["comp"], c["n_ax"])] if c["comp"] == "ОЦ" else PI_RANGE[c["comp"]]


def _stage_mats(c):
    return [c["mat_cw"] if not (c["mat_cw"] == "Al" and k > 0) else "Ti64" for k in range(len(stage_types(c)))]


def _eta_turb_nom(c):
    """номинальные политропные КПД турбин для окон: формулы ступеней при выбранной нагрузке lt и поправки ETA_NOM (оценка)"""
    E = _E()
    if c["turb"] == "О":
        eg = E.eta_axial(c["lt"]) + ETA_NOM["ax"]
    else:
        eg = 0.86 * max(0.05, 1 - ((c["lt"] - 0.70) / 0.70) ** 2) + ETA_NOM["rad"]
    return eg, ETA_NOM["pt"]


def _nominal(c, tgt):
    """номинальный цикл концепции (КПД ступеней по формулам + ETA_NOM): ступени компрессора, ядро, расход"""
    pis, typ, mats, mus = stages_pi(c), stage_types(c), _stage_mats(c), stage_mus(c)
    comp = _comp([(p, max(0.60, eta_comp(t, p, m, k, c["inlet"], mu) + ETA_NOM["comp"]))
                  for k, (p, t, m, mu) in enumerate(zip(pis, typ, mats, mus))])
    eg, ep = _eta_turb_nom(c)
    r = _core(c, comp, eg, ep)
    if "fail" in r:
        return dict(ok=False, fail=r["fail"], comp=comp)
    return dict(ok=True, comp=comp, r=r, G=tgt["P_shaft_kW"] * 1e3 / r["wnet"])


def _margins(c, tgt, pi, t41):
    """запасы номинального цикла при (π, T41): min из w/w_need - 1 - M_W, (t0+100-M_T-EGT)/100, рекуператор (700-M_T-T5)/100,
    алюминиевое колесо (200-M_T-T2_1)/100; < 0 - вне окна"""
    cc = dict(c, pi=pi, t41=t41)
    nm = _nominal(cc, tgt)
    if not nm["ok"]:
        return -1.0, dict(fail=nm["fail"], w=0.0, EGT=None, T5=None)
    r = nm["r"]
    w, egt, t5 = r["wnet"], r["Tstack"] - 273.15, r["T5"] - 273.15
    m = {"w": w / w_need(tgt) - 1.0 - M_W, "EGT": (egt_max(tgt) - M_T - egt) / 100.0}
    if c["recup"] == "печатный":
        m["T5"] = (T_REC_MAX - M_T - t5) / 100.0
    if c["mat_cw"] == "Al":
        m["T_al"] = (MK.MAT["Al"]["T_max"] - M_T - (nm["comp"][0]["T2"] - 273.15)) / 100.0
    worst = min(m, key=m.get)
    return m[worst], dict(w=w / 1e3, EGT=egt, T5=t5, worst=worst, G=nm["G"])


def _run_edges(xs, ms):
    """наибольший непрерывный отрезок xs с ms >= 0; края уточняются линейной интерполяцией запаса до нуля"""
    best, cur = None, None
    for k, m in enumerate(ms):
        if m >= 0:
            cur = (cur[0], k) if cur else (k, k)
            if best is None or cur[1] - cur[0] > best[1] - best[0]:
                best = cur
        else:
            cur = None
    if best is None:
        return None
    i0, i1 = best

    def edge(a, b):                                   # a - вне окна, b - в окне
        ma, mb = ms[a], ms[b]
        t = mb / (mb - ma) if mb - ma > 1e-12 else 0.0
        return xs[b] + (xs[a] - xs[b]) * min(1.0, max(0.0, t)) * 0.999
    lo = xs[i0] if i0 == 0 else edge(i0 - 1, i0)
    hi = xs[i1] if i1 == len(xs) - 1 else edge(i1 + 1, i1)
    return lo, hi


def _pi_window(c, tgt, t41):
    """окно π при T41: сетка по ln π в [lo, min(hi, pi_max_U)]; возврат (окно (a, b) или None, лучшая точка)"""
    lo, hi = pi_range(c)
    hi = min(hi, c["pi_max_U"])
    if hi <= lo * 1.001:
        m, d = _margins(c, tgt, lo, t41)
        return ((lo, lo) if m >= 0 else None), dict(d, pi=lo, m=m)
    xs = [math.log(lo) + (math.log(hi) - math.log(lo)) * k / (N_PI_GRID - 1) for k in range(N_PI_GRID)]
    res = [_margins(c, tgt, math.exp(x), t41) for x in xs]
    ms = [m for m, _ in res]
    kb = max(range(len(ms)), key=lambda k: ms[k])
    best = dict(res[kb][1], pi=math.exp(xs[kb]), m=ms[kb])
    e = _run_edges(xs, ms)
    return ((math.exp(e[0]), math.exp(e[1])) if e else None), best


def _t41_pi_window(c, tgt, tl):
    """окно T41: строки T41 от Tlim вниз (N_T41_GRID), в строке - окно π; края уточняются делением пополам (3 шага).
    Возврат: (Ta, Tb) или None, лучшая точка (наибольший запас)"""
    if tl <= T41_MIN:
        w, best = _pi_window(c, tgt, tl)
        return ((tl, tl) if w else None), dict(best, t41=tl)
    Ts = [tl - (tl - T41_MIN) * k / (N_T41_GRID - 1) for k in range(N_T41_GRID)]
    ok, best, seen = [], None, False
    for T in Ts:
        w, b = _pi_window(c, tgt, T)
        b = dict(b, t41=T)
        if best is None or b["m"] > best["m"]:
            best = b
        ok.append(w is not None)
        if w is not None:
            seen = True
        elif seen:
            break
    if not any(ok):
        return None, best
    ks = [k for k, v in enumerate(ok) if v]
    k_top, k_bot = ks[0], ks[-1]

    def bis(a, b):                                    # a - вне окна, b - в окне
        for _ in range(3):
            mid = 0.5 * (a + b)
            if _pi_window(c, tgt, mid)[0] is not None:
                b = mid
            else:
                a = mid
        return b
    Tb = Ts[k_top] if k_top == 0 else bis(Ts[k_top - 1], Ts[k_top])
    Ta = Ts[k_bot] if (k_bot == len(Ts) - 1 or k_bot + 1 >= len(ok)) else bis(Ts[k_bot + 1], Ts[k_bot])
    return (Ta, Tb), best


_CLOSURE = {}


def best_case(c, mode="mid"):
    """непрерывные оси, лучшие для замыкания точки: μ - середина или верх диапазона, ψ 1,2 или U/C0 0,70, наибольший прирост
    давления, ε 0,92, вход односторонний, деление - лучшее по пределу окружной скорости"""
    cc = dict(c, inlet="одно")
    st = stage_types(cc)
    k = 0.5 if mode == "mid" else 1.0
    cc["mu1"] = MU_RANGE[st[0]][0] + k * (MU_RANGE[st[0]][1] - MU_RANGE[st[0]][0])
    cc["mu2"] = MU_RANGE[st[-1]][0] + k * (MU_RANGE[st[-1]][1] - MU_RANGE[st[-1]][0])
    cc["lt"] = 1.2 if cc["turb"] == "О" else 0.70
    cc["pgain"] = {"RDC": 0.10, "ВР": 0.15}.get(cc["comb"], 0.0)
    cc["eps_rec"] = 0.92 if cc["recup"] == "печатный" else 0.0
    cc["pi_max_U"], cc["split"] = pi_max_by_speed(cc, best_split=True)
    return cc


def closure(c):
    """замыкает ли дискретная комбинация точку KV: есть ли T41 <= Tlim и π в диапазоне, при которых номинальный цикл с лучшими
    в диапазонах непрерывными осями (μ - середина, ψ - 1,2 или U/C0 - 0,70, прирост давления - наибольший, ε = 0,92, деление -
    лучшее, вход односторонний) даёт w_need·(1+M_W) в мягком допуске по выхлопу (и 700 °C перед рекуператором). Строки T41 -
    от предела вниз; если в строке не хватает удельной работы, ниже её тоже не хватит (выход). Кэш по сигнатуре _sig"""
    key = _sig(c)
    if key in _CLOSURE:
        return _CLOSURE[key]
    tgt = target(c["kv"])
    tl = MK.t41_limit(c["turb"], c["cool"], c["mat_tw"], c["m_tw"])[0]
    Ts = [tl] if tl <= T41_MIN else [tl - (tl - T41_MIN) * k / (N_T41_GRID - 1) for k in range(N_T41_GRID)]
    best, ok, mu_ok = None, False, None
    st = stage_types(c)
    for mode in ("mid", "hi"):                  # μ середина диапазона, затем верх (выше π при том же пределе окружной скорости)
        cc = best_case(c, mode)
        for T in Ts:
            w, b = _pi_window(cc, tgt, T)
            b = dict(b, t41=T, mu=mode)
            if best is None or b["m"] > best["m"]:
                best = b
            if w is not None:
                ok, mu_ok = True, mode
                break
            if b.get("worst") == "w" or b.get("fail"):
                break
        if ok:
            break
    res = dict(ok=ok, mu=mu_ok, t41_lim=round(tl, 1), w_need=round(w_need(tgt) * (1 + M_W) / 1e3, 1),
               best=({k: (round(v, 3) if isinstance(v, float) else v) for k, v in best.items()} if best else None))
    if not ok:
        b = best or {}
        why = {"w": "мала удельная работа", "EGT": "горячий выхлоп", "T5": "газ перед рекуператором выше 700 °C",
               "T_al": "алюминиевое колесо выше 200 °C"}.get(b.get("worst"), b.get("fail", "цикл не решается"))
        res["why"] = b.get("worst") or "perf"
        res["reason"] = (f"{why}: лучшая точка T41 {b.get('t41', tl):.0f} °C, π {b.get('pi', 0):.2f} - w {b.get('w', 0):.0f} кДж/кг "
                         f"(нужно {res['w_need']:.0f}), выхлоп {b.get('EGT') or 0:.0f} °C (предел {egt_max(tgt) - M_T:.0f})"
                         + (f", газ перед рекуператором {b.get('T5') or 0:.0f} °C" if c["recup"] == "печатный" else ""))
    _CLOSURE[key] = res
    return res


def _repair(u, arch, c0):
    """если комбинация не замыкает точку - ближайшая замыкающая: меньше изменённых осей, затем раньше по REPAIR_AXES, затем
    ближе по индексу значения; оси вне шаблона архетипа не трогаются. Возврат: (концепция, запись)"""
    cl = closure(c0)
    if cl["ok"]:
        return c0, ""
    lists = {a: _axis_vals(c0, arch, a) for a in REPAIR_AXES}
    cands = []
    for combo in itertools.product(*[lists[a] for a in REPAIR_AXES]):
        ch = {a: v for a, v in zip(REPAIR_AXES, combo) if v != c0[a]}
        if not ch:
            continue
        dist = sum(min((lists[a].index(v) - lists[a].index(c0[a])) % len(lists[a]),
                       (lists[a].index(c0[a]) - lists[a].index(v)) % len(lists[a])) if c0[a] in lists[a] else 1 for a, v in ch.items())
        cands.append(((len(ch), max(REPAIR_AXES.index(a) for a in ch), dist, tuple(sorted(ch.items()))), ch))
    cands.sort(key=lambda x: x[0])
    for _, ov in cands:
        c = _discrete(u, arch, ov)
        c["kv"] = c0["kv"]
        if any(c[a] != v for a, v in ov.items()):
            continue
        if closure(c)["ok"]:
            diff = [f"{a}: {c0[a]} -> {c[a]}" for a in CAT if a != "arch" and c[a] != c0[a]]
            return c, "; ".join(diff) + f" (исходная не замыкает точку - {cl['reason']})"
    return c0, f"нет замыкающей комбинации в шаблоне архетипа - {cl['reason']}"


def _interval(bounds):
    lo = max([b[0] for b in bounds] or [-math.inf]); hi = min([b[1] for b in bounds] or [math.inf])
    return (lo, hi) if lo <= hi else None


def _place(u, win, log=False):
    a, b = win
    if log and a > 0:
        return math.exp(math.log(a) + (math.log(b) - math.log(a)) * float(u))
    return a + (b - a) * float(u)


def concept(i, n, kv="KV-3", seed=SEED):
    """концепция с индексом i выборки объёма n: значения всех осей (детерминированно по i, kv, seed; от n не зависит)"""
    if not 0 <= i < n:
        raise IndexError(i)
    u = dict(zip(AXES, design_matrix(n, seed)[i]))
    arch = _pick(u["arch"], ARCH_KEYS)
    c0 = _discrete(u, arch)
    c0["kv"] = kv
    c, rep = _repair(u, arch, c0)
    c.update(idx=int(i), n=int(n), seed=int(seed), kv=kv)
    c = {k: c[k] for k in ["idx", "n", "seed", "kv", "arch"] + [a for a in CAT if a != "arch"]}
    c["repair"] = rep
    design(c, u, target(kv))
    for k in list(CONT) + ["pi_max_U", "w_need", "pi_lo", "pi_hi", "t41_lo", "t41_hi", "n_gg"]:
        if isinstance(c.get(k), float):
            c[k] = round(c[k], 6)
    return c


def design(c, u, tgt):
    """непрерывные оси в окнах допустимости по номинальному циклу (выборка v2); dflags - условия, окно которых пусто
    (концепция тогда почти наверняка не пройдёт соответствующее правило check - это и есть причина отсева)"""
    E = _E()
    fl = []
    st = stage_types(c)
    fixed = c["rotors"] == 1
    c["mu1"] = _lerp(u["mu1"], *MU_RANGE[st[0]])
    c["mu2"] = _lerp(u["mu2"], *MU_RANGE[st[-1]])
    c["lt"] = _lerp(u["lt"], *(GOOD["psi"] if c["turb"] == "О" else GOOD["nu"]))
    c["eps_rec"] = _lerp(u["eps_rec"], 0.75, 0.92) if c["recup"] == "печатный" else 0.0
    c["pgain"] = _lerp(u["pgain"], 0.0, {"RDC": 0.10, "ВР": 0.15}.get(c["comb"], 0.0))
    c["kcm"] = _lerp(u["kcm"], 0.20, 0.35)
    c["zc"] = int(round(_lerp(u["zc"], 7, 15)))
    c["zt"] = int(round(_lerp(u["zt"], *((9, 17) if c["turb"] == "Р" else (30, 80)))))
    c["zn"] = int(round(_lerp(u["zn"], *((12, 24) if c["turb"] == "Р" else (20, 60)))))
    c["pi_max_U"], c["split"] = pi_max_by_speed(c, best_split=True)
    c["w_need"] = w_need(tgt) / 1e3
    # ---------------------------------------------------------------- T41 и π
    tl = MK.t41_limit(c["turb"], c["cool"], c["mat_tw"], c["m_tw"])[0]
    tw, best = _t41_pi_window(c, tgt, tl)
    if tw is None:                    # окно пусто при выпавших μ, ψ (U/C0), прирост давления, ε - сдвиг к лучшим (best_case) на α
        cl = closure(c)
        if cl["ok"]:
            bc = best_case(dict(c, inlet=c["inlet"]), cl["mu"])
            for alpha in (0.5, 1.0):
                trial = dict(c)
                for k in ("mu1", "mu2", "lt", "pgain", "eps_rec"):
                    trial[k] = c[k] + alpha * (bc[k] - c[k])
                trial["pi_max_U"], trial["split"] = pi_max_by_speed(trial, best_split=True)
                tw2, b2 = _t41_pi_window(trial, tgt, tl)
                if tw2:
                    c.update({k: trial[k] for k in ("mu1", "mu2", "lt", "pgain", "eps_rec", "pi_max_U", "split")})
                    c["repair"] = (c["repair"] + "; " if c["repair"] else "") + f"μ, нагрузка турбины, прирост давления, ε сдвинуты к лучшим на {alpha:g}"
                    tw, best = tw2, b2
                    break
    pw = None
    if tw:
        c["t41"] = _place(u["t41"], tw)
        pw, _ = _pi_window(c, tgt, c["t41"])
        if pw is None:                                    # край окна T41 уточнён с погрешностью - ближе к середине
            c["t41"] = 0.5 * (tw[0] + tw[1])
            pw, _ = _pi_window(c, tgt, c["t41"])
    if pw is None:
        fl.append("T41/π: окно пусто" + (f" ({best.get('worst', best.get('fail', ''))})" if best else ""))
        c["t41"] = best["t41"] if best else max(tl, T41_MIN)
        c["pi"] = best["pi"] if best else pi_range(c)[0]
        c["t41_lo"] = c["t41_hi"] = c["t41"]; c["pi_lo"] = c["pi_hi"] = c["pi"]
    else:
        c["t41_lo"], c["t41_hi"] = (tw[0], tw[1]) if tw else (c["t41"], c["t41"])
        c["pi_lo"], c["pi_hi"] = pw
        c["pi"] = _place(u["pi"], pw, log=True)
    # ---------------------------------------------------------------- вход, деление π, коэффициенты работы
    nm0 = _nominal(c, tgt)
    G0n = nm0["G"] if nm0["ok"] else tgt["G0"] * 1.2
    inlets = _axis_vals(c, c["arch"], "inlet")
    order = inlets[inlets.index(c["inlet"]):] + inlets[:inlets.index(c["inlet"])]
    found = {}
    for inl in order:
        found[inl] = _split_window(dict(c, inlet=inl), tgt, G0n, fixed)
    pick = next((i for i in order if found[i][1] == "хорошее"), None) or next((i for i in order if found[i][1] == "правило"), None)
    tier = "хорошее"
    if pick is None:
        fl.append("ступени компрессора: нет деления π и μ в пределах окружной скорости и коэффициента расхода")
        sw = SPLIT_RANGE.get(c["comp"], (c["split"], c["split"]))
    else:
        if pick != c["inlet"]:
            c["repair"] = (c["repair"] + "; " if c["repair"] else "") + f"inlet: {c['inlet']} -> {pick} (коэффициент расхода)"
        c["inlet"] = pick
        sw, tier = found[pick]
        if tier == "правило":
            fl.append("коэффициент расхода ступеней: только в пределах правила, не в GOOD")
    if c["comp"] in SPLIT_RANGE:
        c["split"] = _place(u["split"], sw)
    nm = _nominal(c, tgt)
    Gn = nm["G"] if nm["ok"] else G0n
    mw = _mu_windows(c, tgt, nm["comp"], Gn, fixed, tier)
    if mw[0] is None:
        fl.append("μ1: окно пусто")
    else:
        c["mu1"] = _place(u["mu1"], mw[0])
    if len(st) > 1:
        if mw[1] is None:
            fl.append("μ2: окно пусто")
        else:
            c["mu2"] = _place(u["mu2"], mw[1])
    nm = _nominal(c, tgt)
    if not nm["ok"]:
        fl.append(f"номинальный цикл: {nm['fail']}")
        _fill_defaults(c, u)
        c["dflags"] = "|".join(fl)
        return c
    r, G = nm["r"], nm["G"]
    Gc = G / c["cores"]
    stn = _stages_nom(c, nm["comp"], Gc)
    # ---------------------------------------------------------------- турбина газогенератора и частота: совместный выбор (вход, nt)
    t = r["t"]
    rho5 = t["p_out"] / (E.Rg(t["b"]) * t["T_out"]); Q5 = t["m"] * Gc / rho5
    nts = _axis_vals(c, c["arch"], "nt")
    nts = nts[nts.index(c["nt"]):] + nts[:nts.index(c["nt"])]
    inl_order = [c["inlet"]] + ([i for i in order if i != c["inlet"]] if not fixed else [])
    lt0 = c["lt"]
    trials = []
    for inl in inl_order:
        stn_i = _stages_nom(dict(c, inlet=inl), nm["comp"], Gc)
        for nt in nts:
            tr = _turb_trial(dict(c, inlet=inl, nt=nt), tgt, t, Q5, stn_i, lt0, fixed)
            trials.append((inl, nt, stn_i, tr))
            if tr["tier"] == "хорошее":
                break
        if trials[-1][3]["tier"] == "хорошее":
            break
    rank = {"хорошее": 0, "правило": 1, "нет": 2}
    inl, nt, stn, tr = min(trials, key=lambda x: rank[x[3]["tier"]])     # первый лучший по порядку перебора
    if nt != c["nt"]:
        c["repair"] = (c["repair"] + "; " if c["repair"] else "") + f"nt: {c['nt']} -> {nt} (нагрузка турбины, частота)"
    if inl != c["inlet"]:
        c["repair"] = (c["repair"] + "; " if c["repair"] else "") + f"inlet: {c['inlet']} -> {inl} (частота газогенератора)"
    c["nt"], c["inlet"], c["lt"] = nt, inl, tr["lt"]
    fl += tr["flags"]
    Ut, n = tr["Ut"], tr["n"]
    if fixed:
        c["n_gg"], c["phi1"] = n, 0.0
    else:
        n = _place(u["phi1"], tr["nwin"], log=True)
        k = next(k for k, s_ in enumerate(stn) if s_["type"] != "О")
        s_ = stn[k]
        c["phi1"] = s_["Gs"] * (math.pi * n / 60) ** 2 / (s_["rho1"] * s_["U"] ** 3)
        c["n_gg"] = n
    # ---------------------------------------------------------------- выход осевой турбины: ca/U по втулке и AN²
    if c["turb"] == "О":
        D = 60 * Ut / (math.pi * n)
        lo = max(GOOD["phit"][0], Q5 / (math.pi * D ** 2 * Ut * HTR_H * (1 - M_GEO)), n ** 2 * Q5 / (Ut * AN2_MAX * (1 - M_GEO)))
        hmax = _blisk_hmax(c, D)
        if hmax is not None:
            lo = max(lo, Q5 / (math.pi * D * Ut * hmax) if hmax > 0 else math.inf)
        if lo > GOOD["phit"][1]:
            fl.append("турбина: втулка, AN² или цельное колесо не проходят и при ca/U = 0,8")
            c["phit"] = GOOD["phit"][1]
        else:
            c["phit"] = _place(u["phit"], (lo, GOOD["phit"][1]))
    else:
        c["phit"] = _lerp(u["phit"], *GOOD["phit"])
    # ---------------------------------------------------------------- силовая турбина
    if c["rotors"] == 2:
        p = r["pt"]
        rho = p["p_out"] / (E.Rg(p["b"]) * p["T_out"]); Qp = p["m"] * G / rho
        npk = tgt["n_pkg"]
        U_an = npk ** 2 * Qp / (GOOD["phit"][1] * AN2_MAX * (1 - M_GEO))
        U_ht = (Qp * (math.pi * npk / 60) ** 2 / (math.pi * GOOD["phit"][1] * HTR_H * (1 - M_GEO))) ** (1 / 3)
        klo = max(0.85, U_an / U_PT_MAX, U_ht / U_PT_MAX)
        if klo > 1.0:
            fl.append("силовая турбина: AN² или втулка не проходят и при 380 м/с и ca/U = 0,8")
            c["kUpt"] = 1.0
        else:
            c["kUpt"] = _place(u["kUpt"], (klo, 1.0))
        Upt = c["kUpt"] * U_PT_MAX
        npts = _axis_vals(c, c["arch"], "npt")
        ok_n = [m for m in npts if p["w"] / (m * Upt ** 2) <= GOOD["psi"][1]] or [m for m in npts if p["w"] / (m * Upt ** 2) <= PSI_MAX]
        if not ok_n:
            fl.append("силовая турбина: нагрузка ψ > 2,5 при любом числе ступеней")
        elif c["npt"] not in ok_n:
            new = ok_n[0] if ok_n[0] > c["npt"] else ok_n[-1]
            c["repair"] = (c["repair"] + "; " if c["repair"] else "") + f"npt: {c['npt']} -> {new} (нагрузка СТ)"
            c["npt"] = new
        Dp = 60 * Upt / (math.pi * npk)
        lo = max(GOOD["phit"][0], Qp / (math.pi * Dp ** 2 * Upt * HTR_H * (1 - M_GEO)), npk ** 2 * Qp / (Upt * AN2_MAX * (1 - M_GEO)))
        c["phipt"] = GOOD["phit"][1] if lo > GOOD["phit"][1] else _place(u["phipt"], (lo, GOOD["phit"][1]))
    else:
        c["kUpt"], c["phipt"] = 0.0, 0.0
    # ---------------------------------------------------------------- габарит: корпус и рекуператор в корпусе
    Dc = max(60 * s["U"] / (math.pi * n) for s in stn)
    side = min(tgt["pkg_W"], tgt["pkg_H"])
    f = 1.0 if c["cores"] == 1 else math.sqrt(c["cores"]) * 1.1
    kmax = (side - 0.5) / (Dc * f)
    if kmax < 1.7:
        fl.append(f"габарит: корпус шире пакета и при D_корп/D2 = 1,7 (нужно {kmax:.2f})")
        c["kD_case"] = 1.7
    else:
        c["kD_case"] = _place(u["kD_case"], (1.7, min(2.1, kmax)))
    if c["recup"] == "печатный" and c["recup_integral"]:
        Denv = c["kD_case"] * Dc * f + 0.5
        Vmax = max(0.0, (side ** 2 - Denv ** 2) * math.pi * Denv / 4)
        A = Vmax * REC_BETA / REC_HDR
        ntu = A * REC_K / (G * 1100.0)
        emax = ntu / (1 + ntu)
        if emax < 0.75:
            fl.append(f"рекуператор в корпусе: в габарит входит ε не выше {emax:.2f}")
            c["eps_rec"] = 0.75
        else:
            c["eps_rec"] = _place(u["eps_rec"], (0.75, min(0.92, emax)))
    c["dflags"] = "|".join(fl)
    return c


def _blisk_hmax(c, D):
    """наибольшая высота последней лопатки, при которой осевое колесо - цельное: ротор заодно или сварной требует блисков
    (правило mono_blisk); литой блиск - Dtip <= 0,50 м (LIM D_axcast), печатный - камера принтера (здесь не ограничивает);
    None - цельное колесо не требуется"""
    if c["turb"] != "О" or c["rotor_join"] not in ("заодно", "сварка") or c["m_tw"] != "литьё":
        return None
    return MK.CAST["turb_axial_blisk"]["max"] * (1 - M_GEO) - D


def _turb_trial(c, tgt, t, Q5, stn, lt0, fixed):
    """турбина газогенератора при данных входе и числе ступеней: окно нагрузки lt (ψ или U/C0) - предел окружной скорости,
    (вал с частотой пакета) втулка, AN² и цельное колесо при ca/U <= 0,8, быстроходность радиальной; окружная скорость;
    частота (вал пакета) или окно частоты газогенератора. Возврат: dict(tier, lt, Ut, n, nwin, flags)"""
    fl = []
    U_lim = MK.u_limit_turb(c["turb"], c["mat_tw"], c["t41"])[0]
    U_max = min(U_lim, u_cap_size(c, d_cap_turb(c)) * (1 - M_GEO))
    n_fix = tgt["n_pkg"] if fixed else None
    tier = "хорошее"
    if c["turb"] == "О":
        U_min = 0.0
        if fixed:
            U_min = max(n_fix ** 2 * Q5 / (GOOD["phit"][1] * AN2_MAX * (1 - M_GEO)),
                        (Q5 * (math.pi * n_fix / 60) ** 2 / (math.pi * GOOD["phit"][1] * HTR_H * (1 - M_GEO))) ** (1 / 3))
        psi_lo = t["w"] / (c["nt"] * U_max ** 2) * (1 + M_GEO)
        psi_hi = t["w"] / (c["nt"] * U_min ** 2) if U_min > 0 else math.inf
        win = _interval([GOOD["psi"], (psi_lo, psi_hi)])
        if win is None:
            win = _interval([(0.0, PSI_MAX), (psi_lo, psi_hi)])
            tier = "правило" if win else "нет"
        if win is None:
            win = (min(psi_lo, PSI_MAX), min(psi_lo, PSI_MAX))
            fl.append("турбина: нагрузка ψ, втулка или AN² - окно пусто")
        lt = min(win[1], max(win[0], lt0))
        Ut = min(U_max, math.sqrt(t["w"] / (c["nt"] * lt)))
    else:
        dhs = _stage_dhs(t, c["nt"]); C0 = math.sqrt(2 * dhs)
        win = _interval([GOOD["nu"], (-math.inf, U_max / C0)])
        if win is None:
            win = _interval([NU_RANGE, (-math.inf, U_max / C0)])
            tier = "правило" if win else "нет"
        if win is None:
            win = (U_max / C0, U_max / C0)
            fl.append("турбина: U/C0 вне 0,5-0,9 при пределе окружной скорости")
        lt = min(win[1], max(win[0], lt0))
        Ut = lt * C0
        if fixed:
            ns = (2 * math.pi * n_fix / 60) * math.sqrt(Q5) / dhs ** 0.75
            if not GOOD["ns_t"][0] <= ns <= GOOD["ns_t"][1]:
                if NS_RANGE[0] <= ns <= NS_RANGE[1]:
                    tier = "правило" if tier == "хорошее" else tier
                else:
                    tier = "нет"; fl.append(f"турбина: быстроходность {ns:.2f} вне 0,25-1,00 при частоте пакета")
    if fixed:
        return dict(tier=tier, lt=lt, Ut=Ut, n=n_fix, nwin=(n_fix, n_fix), flags=fl)
    nwin, how = _n_window(c, stn, Ut, Q5, t)
    if how == "нет":
        fl.append("частота газогенератора: окно пусто (размеры, коэффициенты расхода, турбина) - взята наименее нарушающая")
        tier = "нет"
    elif how == "правило" and tier == "хорошее":
        tier = "правило"
    return dict(tier=tier, lt=lt, Ut=Ut, n=None, nwin=nwin, flags=fl)


def _fill_defaults(c, u):
    c.setdefault("phi1", 0.0); c.setdefault("n_gg", None)
    c["phit"] = _lerp(u["phit"], *GOOD["phit"])
    c["kUpt"] = _lerp(u["kUpt"], 0.85, 1.0) if c["rotors"] == 2 else 0.0
    c["phipt"] = _lerp(u["phipt"], *GOOD["phit"]) if c["rotors"] == 2 else 0.0
    c["kD_case"] = _lerp(u["kD_case"], 1.7, 2.1)


def _stage_dhs(t, nt):
    """изоэнтропический перепад одной ступени турбины при делении степени расширения поровну (от входа), Дж/кг"""
    E = _E()
    ERs = t["ER"] ** (1.0 / nt)
    return E.h(t["T_in"], t["b"]) - E.h(E.T_s(E.s0(t["T_in"], t["b"]) - E.Rg(t["b"]) * math.log(ERs), t["b"]), t["b"])


def _stages_nom(c, comp, Gc):
    E = _E()
    out = []
    for k, (s, t, m, mu) in enumerate(zip(comp, stage_types(c), _stage_mats(c), stage_mus(c))):
        Gs = Gc / 2 if (k == 0 and c["inlet"] == "двух") else Gc
        out.append(dict(type=t, mat=m, mu=mu, Gs=Gs, rho1=s["p1"] / (E.RA * s["T1"]), w=s["w"], pi=s["pi"], T2=s["T2"] - 273.15,
                        U=math.sqrt(s["w"] / mu)))
    return out


def _phi_win(typ, tier):
    if tier == "хорошее":
        return GOOD["phi"][typ]
    lo, hi = PHI_RANGE[typ]
    return lo * (1 + M_GEO), hi * (1 - M_GEO)


def _stage_mu_win(s, c, n_fixed, tier="хорошее"):
    """окно μ ступени: предел окружной скорости материала (и размера при частоте пакета); при частоте пакета - коэффициент
    расхода в GOOD (tier "хорошее") или в диапазоне правила (tier "правило") или втулка >= HTR_MIN_AX (осевая)"""
    lo, hi = MU_RANGE[s["type"]]
    Umax = MK.u_limit_comp(s["mat"], c["m_cw"])[0]
    if n_fixed:
        Umax = min(Umax, u_cap_size(c, d_cap_comp(c)))
    lo = max(lo, s["w"] / (Umax * (1 - M_GEO)) ** 2)
    if n_fixed:
        a = (math.pi * n_fixed / 60) ** 2
        if s["type"] == "О":
            Umin = (s["Gs"] / (0.5 * s["rho1"] * (1 - (HTR_MIN_AX * (1 + M_GEO)) ** 2) * math.pi / 4 / a)) ** (1 / 3)
            hi = min(hi, s["w"] / Umin ** 2)
        else:
            plo, phi_ = _phi_win(s["type"], tier)
            K = s["Gs"] * a / s["rho1"]
            lo = max(lo, s["w"] / ((K / plo) ** (1 / 3)) ** 2)
            hi = min(hi, s["w"] / ((K / phi_) ** (1 / 3)) ** 2)
    return (lo, hi) if lo <= hi else None


def _stage_ok(s, c, n_fixed, tier="хорошее"):
    if s["type"] == "О" and s["pi"] > AX_PI_MAX * (1 - 0.01):
        return None
    if s["type"] == "Д" and s["pi"] > D_PI_MAX * (1 - 0.01):
        return None
    if s["mat"] == "Al" and s["T2"] > MK.MAT["Al"]["T_max"] - 2:
        return None
    return _stage_mu_win(s, c, n_fixed, tier)


def _split_window(c, tgt, G, fixed):
    """окно деления π: каждая ступень - в пределе окружной скорости при каком-нибудь μ из диапазона, осевая π <= 1,8,
    диагональная <= 5, алюминиевое колесо <= 200 °C, (вал с частотой пакета) коэффициент расхода в GOOD, а если так нельзя -
    в диапазоне правила. Возврат: (окно или None, уровень "хорошее"|"правило"|"нет")"""
    n_fixed = tgt["n_pkg"] if fixed else None
    Gc = G / c["cores"]
    tiers = ("хорошее", "правило") if fixed else ("хорошее",)
    if c["comp"] not in SPLIT_RANGE:
        nm = _nominal(c, tgt)
        stn = _stages_nom(c, nm["comp"], Gc)
        for tier in tiers:
            if all(_stage_ok(s, c, n_fixed, tier) for s in stn):
                return (c["split"], c["split"]), tier
        return None, "нет"
    a, b = SPLIT_RANGE[c["comp"]]
    xs = [a + (b - a) * k / (N_SPLIT_GRID - 1) for k in range(N_SPLIT_GRID)]
    stns = []
    for x in xs:
        cc = dict(c, split=x)
        pis, typ, mats, mus = stages_pi(cc), stage_types(cc), _stage_mats(cc), stage_mus(cc)
        comp = _comp([(p, max(0.60, eta_comp(t, p, m, k, cc["inlet"], mu) + ETA_NOM["comp"]))
                      for k, (p, t, m, mu) in enumerate(zip(pis, typ, mats, mus))])
        stns.append((cc, _stages_nom(cc, comp, Gc)))
    for tier in tiers:
        ms = [1.0 if all(_stage_ok(s, cc, n_fixed, tier) for s in stn) else -1.0 for cc, stn in stns]
        w = _run_strict(xs, ms)
        if w:
            return w, tier
    return None, "нет"


def _run_strict(xs, ms):
    best, cur = None, None
    for k, m in enumerate(ms):
        if m > 0:
            cur = (cur[0], k) if cur else (k, k)
            if best is None or cur[1] - cur[0] > best[1] - best[0]:
                best = cur
        else:
            cur = None
    return (xs[best[0]], xs[best[1]]) if best else None


def _mu_windows(c, tgt, comp, G, fixed, tier="хорошее"):
    """окна μ1 (первая, у ОЦ - все осевые) и μ2 (последняя ступень)"""
    n_fixed = tgt["n_pkg"] if fixed else None
    stn = _stages_nom(c, comp, G / c["cores"])
    wins = [_stage_ok(s, c, n_fixed, tier) for s in stn]
    if tier == "хорошее" and not all(wins) and fixed:
        wins = [w or _stage_ok(s, c, n_fixed, "правило") for w, s in zip(wins, stn)]
    if c["comp"] == "ОЦ":
        na = c["n_ax"]
        w1 = _interval([w for w in wins[:na]]) if all(wins[:na]) else None
        return w1, wins[-1]
    return wins[0], (wins[-1] if len(wins) > 1 else None)


def _n_window(c, stn, Ut, Q5, t):
    """окно частоты газогенератора, об/мин: границы «хорошие» (GOOD) и «по правилам»; возврат (окно, "хорошее"|"правило"|"нет")"""
    good, rule = [], []

    def add(lo, hi, both=True, g=True):
        if both or g:
            good.append((lo, hi))
        if both:
            rule.append((lo, hi))
    Dcap = d_cap_comp(c) * (1 - M_GEO)
    for k, s in enumerate(stn):
        add(60 * s["U"] / (math.pi * Dcap), math.inf)
        if s["type"] == "О":
            nmax = (60 * s["U"] / math.pi) * math.sqrt(0.5 * s["rho1"] * s["U"] * (1 - (HTR_MIN_AX * (1 + M_GEO)) ** 2) * math.pi / 4 / s["Gs"])
            add(0.0, nmax)
        else:
            f = lambda phi: (60 / math.pi) * math.sqrt(phi * s["rho1"] * s["U"] ** 3 / s["Gs"])
            good.append((f(GOOD["phi"][s["type"]][0]), f(GOOD["phi"][s["type"]][1])))
            rule.append((f(PHI_RANGE[s["type"]][0]), f(PHI_RANGE[s["type"]][1])))
    Dt = d_cap_turb(c)
    if Dt is not None:
        add(60 * Ut / (math.pi * Dt * (1 - M_GEO)), math.inf)
    if c["turb"] == "Р":
        dhs = _stage_dhs(t, c["nt"])
        f = lambda ns: ns * 60 / (2 * math.pi) * dhs ** 0.75 / math.sqrt(Q5)
        good.append((f(GOOD["ns_t"][0]), f(GOOD["ns_t"][1])))
        rule.append((f(NS_RANGE[0]), f(NS_RANGE[1])))
    else:
        add(0.0, math.sqrt(AN2_MAX * (1 - M_GEO) * GOOD["phit"][1] * Ut / Q5))
        add(0.0, (60 * Ut / math.pi) / math.sqrt(Q5 / (math.pi * GOOD["phit"][1] * Ut * HTR_H * (1 - M_GEO))))
    w = _interval(good)
    if w:
        return w, "хорошее"
    w = _interval(rule)
    if w:
        return w, "правило"
    lo = min(b[0] for b in rule if b[0] > 0); hi = max(b[1] for b in rule if math.isfinite(b[1]))
    xs = [math.log(lo) + (math.log(hi) - math.log(lo)) * k / 40 for k in range(41)] if hi > lo else [math.log(lo)]

    def viol(x):
        v = 0.0
        for a, b in rule:
            v += max(0.0, math.log(max(a, 1e-9)) - x) + (max(0.0, x - math.log(b)) if math.isfinite(b) else 0.0)
        return v
    x = min(xs, key=viol)
    return (math.exp(x), math.exp(x)), "нет"


# =========================================================================== быстрая оценка цикла (уровень 0)
_EX = None


class FastGas:
    """те же таблицы свойств газа gtu_explore (Cantera GRI-Mech 3.0, равномерная сетка 200-2400 K с шагом 2 K) и та же линейная
    интерполяция, но скалярно на чистом Python: gtu_explore.h(T, b) на каждый вызов собирает массив (1-b)·HA + b·HS из 1101
    точек. Совпадение с gtu_explore - до округления (тест test_fastgas_matches_gtu_explore); прочие имена - из модуля"""

    def __init__(self, ex):
        self._ex = ex
        self.T_lo, self.dT = float(ex.TG[0]), float(ex.TG[1] - ex.TG[0])
        self.nT = len(ex.TG)
        self.HA, self.HS, self.SA, self.SS = (ex.HA.tolist(), ex.HS.tolist(), ex.SA.tolist(), ex.SS.tolist())
        self.RA, self.RS, self.HF, self.AFST, self.LHV, self.T0, self.P0 = ex.RA, ex.RS, ex.HF, ex.AFST, ex.LHV, ex.T0, ex.P0

    def __getattr__(self, k):
        return getattr(self._ex, k)

    def _f(self, A, B, T, b):
        x = (T - self.T_lo) / self.dT
        if x <= 0:
            return (1 - b) * A[0] + b * B[0]
        if x >= self.nT - 1:
            return (1 - b) * A[-1] + b * B[-1]
        i = int(x); f = x - i
        y0 = (1 - b) * A[i] + b * B[i]; y1 = (1 - b) * A[i + 1] + b * B[i + 1]
        return y0 + f * (y1 - y0)

    def _inv(self, A, B, v, b):
        lo, hi = 0, self.nT - 1
        if v <= (1 - b) * A[0] + b * B[0]:
            return self.T_lo
        if v >= (1 - b) * A[-1] + b * B[-1]:
            return self.T_lo + self.dT * (self.nT - 1)
        while hi - lo > 1:
            mid = (lo + hi) // 2
            if (1 - b) * A[mid] + b * B[mid] <= v:
                lo = mid
            else:
                hi = mid
        y0 = (1 - b) * A[lo] + b * B[lo]; y1 = (1 - b) * A[hi] + b * B[hi]
        return self.T_lo + self.dT * (lo + (v - y0) / (y1 - y0))

    def h(self, T, b=0.0):
        return self._f(self.HA, self.HS, T, b)

    def s0(self, T, b=0.0):
        return self._f(self.SA, self.SS, T, b)

    def T_h(self, hv, b=0.0):
        return self._inv(self.HA, self.HS, hv, b)

    def T_s(self, sv, b=0.0):
        return self._inv(self.SA, self.SS, sv, b)

    def Rg(self, b=0.0):
        return (1 - b) * self.RA + b * self.RS

    def far_for(self, T3, T4):
        return (self.h(T4) - self.h(T3)) / (self.HF + self.AFST * self.h(T4) - (1 + self.AFST) * self.h(T4, 1.0))

    def bfrac(self, f):
        return f * (1 + self.AFST) / (1 + f)

    @staticmethod
    def eta_axial(psi):
        return 0.885 - 0.035 * max(0.0, psi - 1.2)


def _E():
    """свойства газа: таблицы gtu_explore через FastGas (импорт gtu_explore делает chdir в папку трека - возвращаемся)"""
    global _EX
    if _EX is None:
        cwd = os.getcwd()
        if GTU not in sys.path:
            sys.path.insert(0, GTU)
        import gtu_explore as ex
        os.chdir(cwd)
        _EX = FastGas(ex)
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


DTIP_AX = 1.25        # Dtip / D осевой турбины для предела размера до расчёта высоты лопатки (оценка; frame3 - 1,2 на последней ступени)


def d_cap_turb(c):
    """наибольший диаметр колеса турбины (радиальная - периферия; осевая - средний диаметр), м; None - без предела
    (сборный ротор: диск с лопатками)"""
    m, j = c["m_tw"], c["rotor_join"]
    if m == "печать":
        return MK.largest_print_dim()
    if m == "керамика":
        return MK.CERAMIC["turb_radial"]["max"]
    ax = 1.0 / DTIP_AX if c["turb"] == "О" else 1.0          # осевая: средний диаметр при Dtip ≈ DTIP_AX·D (оценка)
    if j == "заодно":
        return MK.CAST["rotor_mono"]["max"] * ax
    if c["turb"] == "Р":
        return MK.CAST["turb_radial"]["max"] if j == "сборный" else MK.CAST["rotor_section"]["max"]
    if j == "сварка":                                          # сварной ротор требует цельнолитых осевых колёс (правило mono_blisk)
        return MK.CAST["turb_axial_blisk"]["max"] * ax
    return None


def u_cap_size(c, D):
    """при частоте пакета (одновальная схема) предел диаметра - это предел окружной скорости π·n·D/60"""
    if D is None or c["rotors"] != 1:
        return float("inf")
    return math.pi * target(c["kv"])["n_pkg"] * D / 60.0


def stage_pmax(c):
    """наибольшая π каждой ступени в пределе окружной скорости материала (и размера при частоте пакета) при её μ,
    с запасом U_MARGIN и ETA_MARGIN; осевая - не выше AX_PI_MAX, диагональная - D_PI_MAX"""
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
    return pmax


def pi_max_by_speed(c, best_split=False):
    """наибольшая общая π, при которой каждая ступень укладывается в предел окружной скорости своего материала при
    выбранных μ и делении split (best_split - при лучшем делении внутри SPLIT_RANGE; тогда возврат (π, split)).
    Нужна, чтобы выборка не тратилась на заведомо превышенную окружную скорость; правило U_comp всё равно проверяется."""
    pmax = stage_pmax(c)
    if c["comp"] in ("1Ц", "Д"):
        return (pmax[0], c.get("split", 0.5)) if best_split else pmax[0]
    if c["comp"] == "2Ц":
        la, lb = math.log(pmax[0]), math.log(pmax[1])
    else:
        na = c["n_ax"]
        la, lb = na * math.log(min(pmax[:na])), math.log(pmax[-1])
    if best_split:
        lo, hi = SPLIT_RANGE[c["comp"]]
        s = min(hi, max(lo, la / (la + lb)))
    else:
        s = c["split"]
    p = math.exp(min(la / s, lb / (1 - s)))
    return (p, s) if best_split else p


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
    Tstack = Tex
    if rec:                       # выхлоп за рекуператором - баланс энтальпий (исправление скептика физики P4 вместо T5 - 0,93·(T3 - T2))
        m_g, b_g = (m5, b5) if c["rotors"] == 1 else (m45, b45)
        Tstack = E.T_h(E.h(Tex, b_g) - mc * (E.h(T3) - E.h(T2)) / m_g, b_g)
    return dict(wc=wc, wnet=wnet, eta=wnet / (f_tot * E.LHV), f=f_tot, T2=T2, p2=p2, T3=T3, T41=T41, p4=p4,
                T5=Tex, Tstack=Tstack, t=t, pt=pt, xn=xn, xr=xr, xl=xl)


def turb_umax(c):
    """наибольшая окружная скорость турбины газогенератора: предел сплава при T41 и (вал с частотой пакета) предел размера"""
    U, ub = MK.u_limit_turb(c["turb"], c["mat_tw"], c["t41"])
    cap = u_cap_size(c, d_cap_turb(c)) * 0.99
    return min(U, cap), U, ub


def turb_speed(c, t):
    """окружная скорость турбины газогенератора из нагрузки lt (выборка v2): осевая U = √(w/(nt·ψ)), радиальная U = (U/C0)·C0
    ступени; не выше turb_umax. t - турбина из _core. Возврат: U, предел сплава, основание"""
    Umax, U_lim, ub = turb_umax(c)
    if c["turb"] == "О":
        U = math.sqrt(max(t["w"], 1.0) / (c["nt"] * c["lt"]))
    else:
        U = c["lt"] * math.sqrt(2 * max(_stage_dhs(t, c["nt"]), 1.0))
    return min(U, Umax), U_lim, ub


def _radial_stage(t, nt, U):
    """радиальная турбина: изоэнтропический КПД по U/C0 (gtu_explore.turb_radial: w = η·Δh_s, 0,86 при 0,70), степень
    расширения делится поровну; _core расширяет газ политропно, поэтому возвращается политропный эквивалент той же
    изоэнтропической величины на ступени (исправление скептика физики P1). Возврат: η_p, U/C0, Δh_s ступени, η_is"""
    E = _E()
    ERs = t["ER"] ** (1.0 / nt)
    Tin, b = t["T_in"], t["b"]
    dhs = E.h(Tin, b) - E.h(E.T_s(E.s0(Tin, b) - E.Rg(b) * math.log(ERs), b), b)
    C0 = math.sqrt(2 * max(dhs, 1.0)); nu = U / C0
    eta_is = 0.86 * max(0.05, 1 - ((nu - 0.70) / 0.70) ** 2)
    if ERs <= 1.0001:
        return eta_is, nu, dhs, eta_is
    T2 = E.T_h(E.h(Tin, b) - eta_is * dhs, b)
    eta_p = (E.s0(Tin, b) - E.s0(T2, b)) / (E.Rg(b) * math.log(ERs))
    return eta_p, nu, dhs, eta_is


LEAVE_K = 0.5      # не восстановленная диффузором доля прироста выходной кинетической энергии (оценка)


def leaving_loss(phit, U, dhs):
    """снижение КПД осевой турбины при ca/U выше 0,5 (eta_axial отвечает 0,5): LEAVE_K·((phit·U)² - (0,5·U)²)/2 / Δh_s"""
    return LEAVE_K * max(0.0, (phit * U) ** 2 - (CX_T * U) ** 2) / 2.0 / max(dhs, 1.0)


def quick_perf(c, tgt, d_eta=None, d_eta_t=0.0, iters=8):
    """быстрый цикл концепции: КПД турбин - по нагрузке ступеней при окружной скорости концепции (итерации; окружная скорость
    турбины газогенератора выводится из нагрузки lt). d_eta - поправки КПД ступеней компрессора от геометрии (из derive или
    funnel_l0.corrections), d_eta_t - поправка КПД радиальной турбины по быстроходности (funnel_l0)"""
    E = _E()
    pis = stages_pi(c); typ = stage_types(c)
    d_eta = d_eta or [0.0] * len(pis)
    mats = _stage_mats(c)
    mus = stage_mus(c)
    comp = _comp([(p, max(0.60, eta_comp(t, p, m, k, c["inlet"], mu) + de))
                  for k, (p, t, m, de, mu) in enumerate(zip(pis, typ, mats, d_eta, mus))])
    Upt = c["kUpt"] * U_PT_MAX
    eg, ep = _eta_turb_nom(c)
    eg -= ETA_NOM["ax"] if c["turb"] == "О" else ETA_NOM["rad"]
    r, load, Ut, Ulim = None, {}, None, None
    for _ in range(iters):
        r = _core(c, comp, eg, ep)
        if "fail" in r:
            break
        Ut, Ulim, _ = turb_speed(c, r["t"])
        if c["turb"] == "О":
            psi = r["t"]["w"] / (c["nt"] * Ut ** 2)
            eg_n = E.eta_axial(psi) - leaving_loss(c["phit"], Ut, r["t"]["w"] / eg)
            load = dict(psi=psi)
        else:
            eg_n, nu, dhs, eis = _radial_stage(r["t"], c["nt"], Ut)
            eg_n += d_eta_t
            load = dict(nu=nu, dhs=dhs, eta_is=eis)
        ep_n = ep
        if r["pt"]:
            psi_pt = r["pt"]["w"] / (c["npt"] * Upt ** 2)
            ep_n = E.eta_axial(psi_pt) - leaving_loss(c["phipt"], Upt, r["pt"]["w"] / ep)
            load["psi_pt"] = psi_pt
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
REC_HDR = 3.0                    # коллекторы и корпус рекуператора: один множитель на массу (×3 gtu_explore.finalize) и объём;
#                                  сердцевина 1500 м2/м3 со стенкой 0,15 мм весит 1777 кг/м3 (исправление скептика физики P3, оценка)
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
        hp = Q / (math.pi * Dpt * c["phipt"] * Upt)
        pt = dict(U=Upt, D=Dpt, h=hp, htr=(Dpt - hp) / (Dpt + hp), Dtip=Dpt + hp, psi=perf["load"].get("psi_pt"), npt=c["npt"],
                  AN2=math.pi * Dpt * hp * tgt["n_pkg"] ** 2, Q=Q)
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
        V = REC_HDR * A / REC_BETA                          # сердцевина + коллекторы и корпус (P3: тот же множитель, что у массы)
        a_max = 1.1                                        # наибольший куб в камере EP-M2050 (Z 1100 мм)
        mods = max(1, math.ceil(V / a_max ** 3))
        rec = dict(A_m2=A, V_m3=V, modules=mods, side=(V / mods) ** (1 / 3), mass_kg=REC_HDR * A * 0.00015 * 7900)
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
HTR_H = (1 - HTR_MIN_T) / (1 + HTR_MIN_T)  # h/D последней ступени при втулке HTR_MIN_T
AN2_MAX = 6.03e7                         # м2·(об/мин)2: последняя ступень frame3, проверенная прочностью (frame3_fea_results.json,
#                                          parts.turb3: AN² 6,028e7, σ в корне 469 МПа при 100 %, IN713C, T_rel 536 °C, запас к
#                                          длительной прочности 1,41, при разгоне 120 % σ/σв 1,02) - предел до уровня 2 (оценка; P2)
U_PT_MAX = 380.0                         # окружная скорость силовой турбины, м/с (gtu_explore.evaluate, gtu_family)


def _sig6(x):
    """6 значащих цифр (а не 4 знака после запятой: φ 0,019995 не должен становиться 0,02 - пределом правила)"""
    if isinstance(x, float) and x != 0 and math.isfinite(x):
        return float(f"{x:.6g}")
    return x


def _rule(code, ok, val=None, lim=None, text=""):
    f = (lambda x: [_sig6(float(v)) for v in x] if isinstance(x, (list, tuple)) else (_sig6(float(x)) if isinstance(x, float) else x))
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
        R.append(_rule("turb_ansq", tg["AN2"] <= AN2_MAX, tg["AN2"], AN2_MAX,
                       "AN² последней ступени турбины не выше проверенной прочностью ступени frame3 (оценка до уровня 2)"))
    else:
        R.append(_rule("turb_nu", NU_RANGE[0] <= tg["nu"] <= NU_RANGE[1], tg["nu"], list(NU_RANGE), "U/C0 ступени радиальной турбины (оценка)"))
        R.append(_rule("turb_ns", NS_RANGE[0] <= tg["ns"] <= NS_RANGE[1], tg["ns"], list(NS_RANGE), "быстроходность радиальной турбины (оценка)"))
    if g["pt"]:
        R.append(_rule("pt_psi", g["pt"]["psi"] <= PSI_MAX, g["pt"]["psi"], PSI_MAX, "нагрузка ступени силовой турбины"))
        R.append(_rule("pt_htr", g["pt"]["htr"] >= HTR_MIN_T, g["pt"]["htr"], HTR_MIN_T, "втулка/периферия последней ступени СТ (оценка)"))
        R.append(_rule("pt_ansq", g["pt"]["AN2"] <= AN2_MAX, g["pt"]["AN2"], AN2_MAX,
                       "AN² последней ступени СТ не выше проверенной прочностью ступени frame3 (оценка до уровня 2)"))
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
            lo, hi = lim                                    # вне диапазона - расстояние до ближней границы, не меньше 0
            v += max(0.0, lo - val) / abs(lo) + max(0.0, val - hi) / abs(hi)
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
