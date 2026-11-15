#!/usr/bin/env python3
"""Пространство сквозных концепций «от топлива до вала» (уровень облика): геном, декодирование в концепцию блоков,
дескрипторы ниш MAP-Elites, генератор концепций по индексу с делением на части.

Зачем. Фокус владельца 15.11.2026: «подбираем топовые конструкции всего процесса — от поступления топлива до вала;
глубже в расчёты не погружаемся; концентрируемся на геометрии, которая даст наибольшую эффективность и значительно
повысит КПД на тонну агрегата». Здесь — пространство поиска: дискретные выборы шести блоков (process/blocks) и способа
изготовления плюс непрерывные параметры облика (π и её распределение, T41 по классу охлаждения, эффективность
рекуператора, деление расширения при подогреве, доля риформинга, коэффициенты работы ступеней, нагрузка турбины).

Геном — вектор x в [0, 1)^D (D = len(GENES)). Дискретный ген: вариант = по накопленным весам (веса — доля начальной
выборки, оценка; MAP-Elites их не использует). Непрерывный ген: величина ВНУТРИ окна, которое зависит от уже выбранных
блоков (π - в диапазоне типа компрессора funnel_space.PI_RANGE, T41 - от 850 °C до предела класса охлаждения и сплава
make_model.t41_limit и т. д.), поэтому каждое значение гена осмысленно. Декодирование детерминировано; где выбор
бессмыслен по построению (класс P при литой горячей части, оболочка без печатного корпуса), он заменяется ближайшим
осмысленным с записью в repairs. Правила совместимости блоков (compat) здесь не применяются - их проверяет оценка.

Дескрипторы ниш (сетка 32 × 6 × 24 × 4 = 18 432): топология цикла (рекуператор, промохлаждение, повторный подогрев, IBC,
химическая рекуперация — 2^5), схема камеры (6), архитектура (компрессор 4 × турбина 2 × вал 3), изготовление (4).

Генератор по индексу: концепция i точки KV — строка i последовательности Соболя (scipy.stats.qmc.Sobol, scramble, сид
SEED + номер точки; scipy закреплён 1.18.1, как у задач воронки — скремблирование зависит от версии). Часть k из P
получает индексы [k·N//P, (k+1)·N//P) — разбиение ровно одно, не зависит от числа процессов и машин.

  python3 process_space.py --describe            # гены, варианты, окна, сетка ниш (JSON в stdout)
  python3 process_space.py --sample 5 --kv KV-3  # первые концепции генератора
Всё - оценка уровня облика; окна и веса - принятые величины трека (ссылки в GENES) и собственные оценки.
"""
import argparse
import hashlib
import json
import math
import os
import sys
import warnings

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in sys.path:
    sys.path.insert(0, HERE)
from blocks import _track as T  # noqa: E402

SEED = 20261115
KV_INDEX = {"KV-1": 1, "KV-3": 3, "KV-6": 6}

# =========================================================================== изготовление (ось ниши)
# Семья задаёт способ и сплавы колёс, горячей части, камеры и корпусов в концепции блоков и способ корпусов в модели массы
# (process/mass PROCESS_SETS; «гибрид» регистрируется process_eval: холодная часть - точное литьё КНР, горячая и камера - LPBF).
# Сплавы - списки воронки funnel_space.MAT_TW / MAT_HOT по способу.
MAKE = {
    "традиционная": dict(name="литьё в песок, кованые колёса компрессора (опорная практика Solar)", comp_make="поковка", turb_make="литьё",
                         hot_make="литьё", comb_make="обработка", case_make="литьё", mass_process="cast_sand"),
    "литьё КНР": dict(name="точное литьё КНР: колёса, корпуса, горячая часть (навык casting-spec)", comp_make="литьё", turb_make="литьё",
                      hot_make="литьё", comb_make="литьё", case_make="литьё", mass_process="cast_invest"),
    "гибрид": dict(name="гибрид: литые колёса и корпуса КНР + печатная горячая часть и камера LPBF", comp_make="литьё",
                   turb_make="литьё", hot_make="печать", comb_make="печать", case_make="литьё", mass_process="hybrid"),
    "печать": dict(name="печать LPBF: колёса, горячая часть, камера, корпуса (оболочка возможна)", comp_make="печать",
                   turb_make="печать", hot_make="печать", comb_make="печать", case_make="печать", mass_process="lpbf"),
}
MAKE_KEYS = list(MAKE)

# =========================================================================== гены
# kind: "cat" - варианты opts (веса w - доля начальной выборки, оценка); "cont" - доля окна (окно строит decode)
GENES = [
    # ---- топология цикла и топливо (блоки 1, 2, 5)
    dict(id="recup", kind="cat", opts=[0, 1], w=[0.5, 0.5], block="hx", what="рекуператор: нет / есть (ε - ген eps)"),
    dict(id="ic", kind="cat", opts=["нет", "воздух"], w=[0.6, 0.4], block="compressor", what="промохлаждение воздух-воздух"),
    dict(id="rh", kind="cat", opts=[0, 1], w=[0.6, 0.4], block="hx", what="повторный подогрев (межтурбинная камера)"),
    dict(id="ibc", kind="cat", opts=[0, 1], w=[0.85, 0.15], block="hx", what="обращённый цикл Брайтона на выхлопе"),
    dict(id="fuel", kind="cat", opts=["нет", "подогрев", "паровой риформинг", "сухой риформинг"], w=[0.4, 0.3, 0.2, 0.1],
         block="fuel", what="подготовка топлива"),
    # ---- камера (блок 3)
    dict(id="scheme", kind="cat", opts=["LPP", "FLOX", "TVC", "CAT", "RDC", "ВР"], w=None, block="combustor",
         what="схема камеры (LPP - пассивный пилот при одном клапане, решение владельца 15.11.2026)"),
    dict(id="chamber", kind="cat", opts=["кольцевая", "выносная"], w=[0.7, 0.3], block="combustor", what="тип камеры"),
    # ---- архитектура (блоки 2, 4, 6)
    dict(id="comp", kind="cat", opts=["1Ц", "2Ц", "Д", "ОЦ"], w=None, block="compressor", what="компрессор: архитектура воронки"),
    dict(id="n_ax", kind="cat", opts=[1, 2, 3], w=None, block="compressor", what="осевых ступеней (только ОЦ)"),
    dict(id="inlet", kind="cat", opts=["одно", "двух"], w=[0.6, 0.4], block="compressor", what="вход колеса одно- / двусторонний"),
    dict(id="turb", kind="cat", opts=["О", "Р"], w=None, block="turbine", what="турбина газогенератора: осевая / радиальная"),
    dict(id="nt", kind="cat", opts=["auto", 1, 2, 3], w=[0.55, 0.15, 0.15, 0.15], block="turbine",
         what="ступеней турбины («auto» - наименьшее при ψ <= 1,6 / U/C0 >= 0,6; радиальная - не больше 2)"),
    dict(id="cool", kind="cat", opts=["U", "N", "K", "P", "C"], w=[0.25, 0.2, 0.2, 0.2, 0.15], block="turbine",
         what="класс охлаждения (C - керамика Si3N4 / CMC, радиальная)"),
    dict(id="drive", kind="cat", opts=["прямой", "свободная СТ", "быстроходный"], w=None, block="shaft", what="вал и стык пакета"),
    # ---- изготовление и материалы
    dict(id="make", kind="cat", opts=MAKE_KEYS, w=None, block="build", what="семья изготовления (ось ниши)"),
    dict(id="tmat", kind="cont", block="turbine", what="сплав колеса турбины: индекс в списке способа (funnel_space.MAT_TW)"),
    dict(id="hmat", kind="cont", block="turbine", what="сплав горячей части: индекс в списке способа (funnel_space.MAT_HOT)"),
    dict(id="cmat", kind="cat", opts=["Ti64", "Al", "Ni"], w=[0.6, 0.25, 0.15], block="compressor", what="материал колёс компрессора"),
    dict(id="join", kind="cat", opts=["сборный", "сварка", "заодно"], w=None, block="build", what="соединение ротора"),
    dict(id="bear", kind="cat", opts=["масло", "магнит"], w=[0.6, 0.4], block="build", what="подшипники"),
    dict(id="case", kind="cat", opts=["картридж", "оболочка", "разъём"], w=[0.5, 0.3, 0.2], block="build",
         what="корпус (оболочка - только печатный корпус LPBF)"),
    dict(id="rmake", kind="cat", opts=["фольга", "печать"], w=[0.5, 0.5], block="hx", what="рекуператор: фольга первичной поверхности / LPBF"),
    dict(id="icpos", kind="cont", block="compressor", what="место промохладителя: после ступени k (индекс среди возможных)"),
    # ---- непрерывные параметры облика
    dict(id="pi", kind="cont", block="compressor", what="π: доля ln-диапазона типа компрессора (funnel_space.PI_RANGE)"),
    dict(id="split", kind="cont", block="compressor", what="распределение π: доля ln π первой ступени (2Ц 0,35-0,65) / осевой части (ОЦ 0,12-0,60)"),
    dict(id="t41", kind="cont", block="turbine", what="T41 по классу: от предела минус 200 K (не ниже 850 °C) до предела класса и сплава (make_model.t41_limit)"),
    dict(id="eps", kind="cont", block="hx", what="ε рекуператора 0,70-0,92 (диапазон блока hx.RECUP_EPS)"),
    dict(id="xrh", kind="cont", block="hx", what="подогрев: доля ln p до межтурбинной камеры 0,35-0,65 (organic/cycles A x_rh)"),
    dict(id="trh", kind="cont", block="hx", what="подогрев: температура после камеры от max(850, T41-150) до предела класса"),
    dict(id="reform", kind="cont", block="fuel", what="доля риформинга: S/C 1,5-3,0 (паровой) / рециркуляция R 0,5-4 по ln (сухой) - organic/cycles A"),
    dict(id="mu1", kind="cont", block="compressor", what="коэффициент работы первой (ОЦ - осевых) ступени: окно funnel_space.MU_RANGE"),
    dict(id="mu2", kind="cont", block="compressor", what="коэффициент работы последней ступени"),
    dict(id="kcm", kind="cont", block="compressor", what="cm2/U2 на выходе колеса (ширина b2) 0,20-0,35 (funnel_space CONT kcm)"),
    dict(id="lt", kind="cont", block="turbine", what="нагрузка турбины: осевая ψ 1,2-2,3 / радиальная U/C0 0,62-0,78 (CONT lt)"),
    dict(id="phit", kind="cont", block="turbine", what="ca/U на выходе турбины газогенератора 0,45-0,80 (CONT phit)"),
    dict(id="kUpt", kind="cont", block="turbine", what="окружная скорость силовой турбины в долях 380 м/с: 0,85-1,00 (CONT kUpt)"),
    dict(id="phipt", kind="cont", block="turbine", what="ca/U на выходе силовой турбины 0,45-0,80 (CONT phipt)"),
    dict(id="kD", kind="cont", block="build", what="диаметр корпуса / наибольший D2 1,7-2,1 (CONT kD_case)"),
]
GENE_IDS = [g["id"] for g in GENES]
GI = {g: i for i, g in enumerate(GENE_IDS)}
D = len(GENES)

RANGES = dict(
    eps=(0.70, 0.92), xrh=(0.35, 0.65), sc=(1.5, 3.0), R_egr=(0.5, 4.0), kcm=(0.20, 0.35), lt_ax=(1.2, 2.3), lt_rad=(0.62, 0.78),
    phit=(0.45, 0.80), kUpt=(0.85, 1.00), phipt=(0.45, 0.80), kD=(1.7, 2.1), trh_below=150.0, t41_below=200.0,
    basis="окна непрерывных осей воронки (funnel_space.CONT: kcm, lt, phit, kUpt, phipt, kD_case), блока hx (ε 0,70-0,92), "
          "перебора циклов (organic/cycles A: x_rh 0,35-0,65, sc 1,5-3, R_egr 0,5-4); T41 - от предела класса минус 200 K (не ниже 850 °C, "
          "funnel_space.T41_MIN) до предела: ниже удельной работы не хватает на допуск расхода +20 % (проба 300 концепций: 49 % отказов "
          "по расходу) - оценка; trh - на 150 K ниже T41 до предела класса (оценка)",
)

# =========================================================================== ниши
TOPO_BITS = (("recup", "Р"), ("ic", "ПО"), ("rh", "ПП"), ("ibc", "IBC"), ("tcr", "ХР"))
TOPO_LABELS = []
for _k in range(32):
    _s = [lab for b, (_, lab) in enumerate(TOPO_BITS) if _k >> b & 1]
    TOPO_LABELS.append("+".join(_s) if _s else "простой")
SCHEMES = ["LPP", "FLOX", "TVC", "CAT", "RDC", "ВР"]
COMPS, TURBS, DRIVES = ["1Ц", "2Ц", "Д", "ОЦ"], ["О", "Р"], ["прямой", "свободная СТ", "быстроходный"]
ARCH_LABELS = [f"{c}/{t}/{d}" for c in COMPS for t in TURBS for d in DRIVES]
DESC_NAMES = ["топология", "камера", "архитектура", "изготовление"]
DESC_LABELS = [TOPO_LABELS, SCHEMES, ARCH_LABELS, MAKE_KEYS]
DESC_BINS = [len(x) for x in DESC_LABELS]
N_CELLS = int(np.prod(DESC_BINS))
DESC_NOTE = ("топология: Р - рекуператор, ПО - промохлаждение, ПП - повторный подогрев, IBC - обращённый Брайтон, ХР - химическая "
             "рекуперация (паровой или сухой риформинг); подогрев топлива - вариант внутри ниши; архитектура: компрессор/турбина/вал")


def grid_dict():
    """сетка дескрипторов в формате evo/core/archive.GridSpec (целые бины: границы k-0,5)"""
    return dict(names=DESC_NAMES, edges=[[k - 0.5 for k in range(b + 1)] for b in DESC_BINS], labels=DESC_LABELS)


def descriptors(C):
    """номера бинов ниши по концепции блоков (после decode)"""
    topo = 0
    hx, cb = C["hx"], C["combustor"]
    flags = dict(recup=bool(hx.get("recup")), ic=C["compressor"].get("intercool", "нет") != "нет", rh=bool(hx.get("reheat")),
                 ibc=bool(hx.get("ibc")), tcr=C["fuel"].get("prep") in ("паровой риформинг", "сухой риформинг"))
    for b, (k, _) in enumerate(TOPO_BITS):
        if flags[k]:
            topo |= 1 << b
    arch = (COMPS.index(C["compressor"]["arch"]) * len(TURBS) + TURBS.index(C["turbine"]["type"])) * len(DRIVES) + \
        DRIVES.index(C["shaft"]["drive"])
    return (topo, SCHEMES.index(cb["scheme"]), arch, MAKE_KEYS.index(C["_make"]))


def cell_label(cell):
    return " | ".join(DESC_LABELS[j][int(c)] for j, c in enumerate(cell))


def cell_index(cell):
    i = 0
    for c, b in zip(cell, DESC_BINS):
        i = i * b + int(c)
    return i


# =========================================================================== декодирование
def _cat(g, x):
    opts = g["opts"]
    w = g.get("w")
    if not w:
        return opts[min(int(x * len(opts)), len(opts) - 1)]
    cw = np.cumsum(w) / float(sum(w))
    return opts[min(int(np.searchsorted(cw, x, side="right")), len(opts) - 1)]


def _lin(f, lo, hi):
    return lo + float(f) * (hi - lo)


def _log(f, lo, hi):
    return math.exp(math.log(lo) + float(f) * (math.log(hi) - math.log(lo)))


def _pick(lst, f):
    return lst[min(int(float(f) * len(lst)), len(lst) - 1)]


def decode(x, kv="KV-3"):
    """геном -> концепция блоков (dict для blocks.evaluate) с полями _make, _mass_process, _repairs, _genes.
    Ошибки схемы не возникают по построению: каждый ген отображается в допустимое по описанию значение."""
    FS, MK = T.FS(), T.MK()
    x = np.clip(np.asarray(x, float), 0.0, 1.0 - 1e-12)
    if x.shape != (D,):
        raise ValueError(f"геном: длина {x.shape}, нужно {D}")
    v = {g["id"]: (_cat(g, x[i]) if g["kind"] == "cat" else float(x[i])) for i, g in enumerate(GENES)}
    rep = []
    mk = MAKE[v["make"]]
    # ---- турбина: тип, класс, способ и сплавы
    ttype, cool = v["turb"], v["cool"]
    if cool == "C":
        ttype, t_make, h_make = "Р", "керамика", "керамика"
    else:
        t_make, h_make = mk["turb_make"], mk["hot_make"]
    if cool == "P" and h_make != "печать":
        cool = "N"
        rep.append("класс P при литой горячей части -> N (печатные каналы только в печатном сопловом аппарате)")
    tmats = [m for m in FS.MAT_TW[t_make] if not (m == "SC" and ttype == "Р")]
    tmat = _pick(tmats, v["tmat"])
    hmat = _pick(FS.MAT_HOT[h_make], v["hmat"])
    nt = v["nt"]
    if nt != "auto" and ttype == "Р" and nt > 2:
        nt = 2
        rep.append("радиальная турбина: не больше 2 ступеней")
    tl = MK.t41_limit(ttype, cool, tmat, t_make)[0]
    t41_lo = max(FS.T41_MIN, tl - RANGES["t41_below"])
    t41 = round(_lin(v["t41"], t41_lo, max(tl, t41_lo)), 1)
    # ---- компрессор
    arch = v["comp"]
    n_ax = int(v["n_ax"]) if arch == "ОЦ" else 0
    lo, hi = FS.PI_RANGE[(arch, n_ax)] if arch == "ОЦ" else FS.PI_RANGE[arch]
    pi = round(_log(v["pi"], lo, hi), 4)
    stypes = {"1Ц": ["Ц"], "Д": ["Д"], "2Ц": ["Ц", "Ц"], "ОЦ": ["О"] * n_ax + ["Ц"]}[arch]
    split = None
    if arch in FS.SPLIT_RANGE:
        split = round(_lin(v["split"], *FS.SPLIT_RANGE[arch]), 4)
    mu_lo1, mu_hi1 = FS.MU_RANGE[stypes[0]]
    mu_lo2, mu_hi2 = FS.MU_RANGE[stypes[-1]]
    mu1 = round(_lin(v["mu1"], mu_lo1, mu_hi1), 4)
    mu2 = round(_lin(v["mu2"], mu_lo2, mu_hi2), 4)
    mus = [mu1] * (len(stypes) - 1) + [mu2] if len(stypes) > 1 else [mu1]
    ic = v["ic"]
    ic_after = None
    if ic != "нет":
        if len(stypes) < 2:
            ic = "нет"
            rep.append("промохлаждение при одной ступени сжатия (1Ц, Д) -> нет (охладителю негде стоять)")
        else:
            ic_after = _pick(list(range(len(stypes) - 1)), v["icpos"])
    cmat = v["cmat"]
    join = v["join"]
    # ---- корпус
    case = v["case"]
    if case == "оболочка" and mk["case_make"] != "печать":
        case = "картридж"
        rep.append("оболочка без печатного корпуса -> картридж")
    # ---- теплообменники и топливо
    recup = round(_lin(v["eps"], *RANGES["eps"]), 4) if v["recup"] else 0.0
    rh = bool(v["rh"])
    trh_lo = max(FS.T41_MIN, t41 - RANGES["trh_below"])
    trh = round(_lin(v["trh"], trh_lo, max(tl, trh_lo)), 1) if rh else None
    xrh = round(_lin(v["xrh"], *RANGES["xrh"]), 4)
    prep = v["fuel"]
    sc = round(_lin(v["reform"], *RANGES["sc"]), 4)
    R_egr = round(_log(v["reform"], *RANGES["R_egr"]), 4)
    lt = round(_lin(v["lt"], *(RANGES["lt_ax"] if ttype == "О" else RANGES["lt_rad"])), 4)
    C = dict(
        kv=kv,
        fuel=dict(prep=prep, sc=sc if prep == "паровой риформинг" else 2.0, R_egr=R_egr if prep == "сухой риформинг" else 1.0),
        compressor=dict(arch=arch, n_ax=n_ax, inlet=v["inlet"], mat=cmat, make=mk["comp_make"], pi=pi, split=split, mu=mus,
                        kcm=round(_lin(v["kcm"], *RANGES["kcm"]), 4), intercool=ic, ic_after=ic_after),
        combustor=dict(scheme=v["scheme"], chamber=v["chamber"]),
        turbine=dict(type=ttype, nt=nt, cool=cool, mat=tmat, make=t_make, hot_make=h_make, hot_mat=hmat, t41=t41, lt=lt,
                     phit=round(_lin(v["phit"], *RANGES["phit"]), 4), kUpt=round(_lin(v["kUpt"], *RANGES["kUpt"]), 4),
                     phipt=round(_lin(v["phipt"], *RANGES["phipt"]), 4)),
        hx=dict(recup=recup, recup_make=v["rmake"], reheat=rh, x_rh=xrh, t_rh=trh, ibc=bool(v["ibc"])),
        shaft=dict(drive=v["drive"]),
        build=dict(rotor_join=join, bear=v["bear"], case=case, comb_make=mk["comb_make"], case_make=mk["case_make"],
                   kD_case=round(_lin(v["kD"], *RANGES["kD"]), 4)),
    )
    C["_make"] = v["make"]
    C["_mass_process"] = mk["mass_process"]
    C["_repairs"] = rep
    C["_t41_limit"] = tl
    return C


def _cat_x(gid, val):
    g = GENES[GI[gid]]
    opts = g["opts"]
    k = opts.index(val)
    w = g.get("w")
    if not w:
        return (k + 0.5) / len(opts)
    cw = np.concatenate([[0.0], np.cumsum(w) / float(sum(w))])
    return float(0.5 * (cw[k] + cw[k + 1]))


def _inv_lin(v, lo, hi):
    return 0.5 if hi <= lo else min(max((float(v) - lo) / (hi - lo), 0.0), 1.0 - 1e-9)


def _inv_log(v, lo, hi):
    return _inv_lin(math.log(v), math.log(lo), math.log(hi))


def make_family(C):
    """семья изготовления по способам концепции блоков (для encode опорных концепций): печатное колесо турбины или корпус -
    «печать»; печатная горячая часть - «гибрид»; кованые колёса компрессора - «традиционная»; иначе - «литьё КНР»"""
    t, b, c = C["turbine"], C["build"], C["compressor"]
    if t.get("make") == "печать" or b.get("case_make") == "печать":
        return "печать"
    if t.get("hot_make") == "печать":
        return "гибрид"
    if c.get("make") == "поковка":
        return "традиционная"
    return "литьё КНР"


def encode(C, make=None):
    """приближённая обратная к decode: концепция блоков (полная, после blocks.concept.normalize) -> геном. Нужна, чтобы начать
    острова с известных допустимых концепций (родители: опорные примеры блоков, выбор владельца, архив прошлого этапа).
    decode(encode(C)) воспроизводит дискретные выборы и непрерывные величины внутри окон (вне окна - ближайшая граница)."""
    FS, MK = T.FS(), T.MK()
    c, t, hx, f, b = C["compressor"], C["turbine"], C["hx"], C["fuel"], C["build"]
    fam = make or C.get("_make") or make_family(C)
    x = np.full(D, 0.5)
    x[GI["recup"]] = _cat_x("recup", 1 if hx.get("recup") else 0)
    x[GI["ic"]] = _cat_x("ic", c.get("intercool", "нет"))
    x[GI["rh"]] = _cat_x("rh", 1 if hx.get("reheat") else 0)
    x[GI["ibc"]] = _cat_x("ibc", 1 if hx.get("ibc") else 0)
    x[GI["fuel"]] = _cat_x("fuel", f.get("prep", "нет"))
    x[GI["scheme"]] = _cat_x("scheme", C["combustor"]["scheme"])
    x[GI["chamber"]] = _cat_x("chamber", C["combustor"].get("chamber", "кольцевая"))
    x[GI["comp"]] = _cat_x("comp", c["arch"])
    x[GI["n_ax"]] = _cat_x("n_ax", max(1, int(c.get("n_ax") or 1)))
    x[GI["inlet"]] = _cat_x("inlet", c.get("inlet", "одно"))
    x[GI["turb"]] = _cat_x("turb", t["type"])
    nt = t.get("nt", "auto")
    x[GI["nt"]] = _cat_x("nt", nt if nt == "auto" else int(nt))
    x[GI["cool"]] = _cat_x("cool", t["cool"])
    x[GI["drive"]] = _cat_x("drive", C["shaft"]["drive"])
    x[GI["make"]] = _cat_x("make", fam)
    x[GI["cmat"]] = _cat_x("cmat", c.get("mat", "Ti64"))
    x[GI["join"]] = _cat_x("join", b.get("rotor_join", "сборный"))
    x[GI["bear"]] = _cat_x("bear", b.get("bear", "масло"))
    x[GI["case"]] = _cat_x("case", b.get("case", "картридж"))
    x[GI["rmake"]] = _cat_x("rmake", hx.get("recup_make", "фольга"))
    mk = MAKE[fam]
    t_make = "керамика" if t["cool"] == "C" else mk["turb_make"]
    h_make = "керамика" if t["cool"] == "C" else mk["hot_make"]
    tmats = [m for m in FS.MAT_TW[t_make] if not (m == "SC" and t["type"] == "Р")]
    x[GI["tmat"]] = (tmats.index(t["mat"]) + 0.5) / len(tmats) if t["mat"] in tmats else 0.5
    hm = FS.MAT_HOT[h_make]
    x[GI["hmat"]] = (hm.index(t["hot_mat"]) + 0.5) / len(hm) if t.get("hot_mat") in hm else 0.5
    tmat = _pick(tmats, x[GI["tmat"]])
    tl = MK.t41_limit(t["type"], t["cool"], tmat, t_make)[0]
    t41_lo = max(FS.T41_MIN, tl - RANGES["t41_below"])
    x[GI["t41"]] = _inv_lin(t["t41"], t41_lo, max(tl, t41_lo))
    arch, n_ax = c["arch"], (int(c.get("n_ax") or 1) if c["arch"] == "ОЦ" else 0)
    lo, hi = FS.PI_RANGE[(arch, n_ax)] if arch == "ОЦ" else FS.PI_RANGE[arch]
    x[GI["pi"]] = _inv_log(c["pi"], lo, hi)
    if arch in FS.SPLIT_RANGE and c.get("split") is not None:
        x[GI["split"]] = _inv_lin(c["split"], *FS.SPLIT_RANGE[arch])
    stypes = {"1Ц": ["Ц"], "Д": ["Д"], "2Ц": ["Ц", "Ц"], "ОЦ": ["О"] * max(n_ax, 1) + ["Ц"]}[arch]
    mus = c.get("mu") or []
    if mus:
        x[GI["mu1"]] = _inv_lin(mus[0], *FS.MU_RANGE[stypes[0]])
        x[GI["mu2"]] = _inv_lin(mus[-1], *FS.MU_RANGE[stypes[-1]])
    if c.get("ic_after") is not None and len(stypes) > 1:
        x[GI["icpos"]] = (int(c["ic_after"]) + 0.5) / (len(stypes) - 1)
    if c.get("kcm") is not None:
        x[GI["kcm"]] = _inv_lin(c["kcm"], *RANGES["kcm"])
    if hx.get("recup"):
        x[GI["eps"]] = _inv_lin(hx["recup"], *RANGES["eps"])
    if hx.get("reheat"):
        x[GI["xrh"]] = _inv_lin(hx.get("x_rh", 0.5), *RANGES["xrh"])
        t41d = _lin(x[GI["t41"]], t41_lo, max(tl, t41_lo))
        trh_lo = max(FS.T41_MIN, t41d - RANGES["trh_below"])
        x[GI["trh"]] = _inv_lin(hx.get("t_rh") or t41d, trh_lo, max(tl, trh_lo))
    if f.get("prep") == "паровой риформинг":
        x[GI["reform"]] = _inv_lin(f.get("sc", 2.0), *RANGES["sc"])
    elif f.get("prep") == "сухой риформинг":
        x[GI["reform"]] = _inv_log(f.get("R_egr", 1.0), *RANGES["R_egr"])
    if t.get("lt") is not None:
        x[GI["lt"]] = _inv_lin(t["lt"], *(RANGES["lt_ax"] if t["type"] == "О" else RANGES["lt_rad"]))
    for gid, key in (("phit", "phit"), ("kUpt", "kUpt"), ("phipt", "phipt")):
        if t.get(key) is not None:
            x[GI[gid]] = _inv_lin(t[key], *RANGES[gid])
    if b.get("kD_case") is not None:
        x[GI["kD"]] = _inv_lin(b["kD_case"], *RANGES["kD"])
    return np.clip(x, 0.0, 1.0 - 1e-12)


SEEDS_FILE = os.path.join(HERE, "process_seeds.json")
FUNNEL_V3_SELECT = {"KV-3": ("l0_v3_probe/funnel_l0_select_KV-3_n300.json", 300), "KV-1": ("l0_v3_probe/funnel_l0_select_KV-1_n300.json", 300),
                    "KV-6": ("l0_v3_probe/funnel_l0_select_KV-6_n1000.json", 1000)}


def funnel_to_blocks(fc):
    """концепция воронки (funnel_space.concept, выборка v3) -> концепция блоков (обратное к blocks.concept.to_funnel)"""
    st = {"1Ц": ["Ц"], "Д": ["Д"], "2Ц": ["Ц", "Ц"], "ОЦ": ["О"] * max(int(fc.get("n_ax") or 1), 1) + ["Ц"]}[fc["comp"]]
    mus = [float(fc["mu1"])] * (len(st) - 1) + [float(fc["mu2"])] if len(st) > 1 else [float(fc["mu1"])]
    comb = fc["comb"]
    scheme, chamber = ({"RDC": ("RDC", "кольцевая"), "ВР": ("ВР", "кольцевая")}.get(comb) or ("LPP", comb))
    case = fc.get("case_arch") or ("картридж" if fc.get("cartridge") else "разъём")
    return dict(
        kv=fc["kv"], fuel=dict(prep="нет"),
        compressor=dict(arch=fc["comp"], n_ax=int(fc.get("n_ax") or 0), inlet=fc["inlet"], mat=fc["mat_cw"], make=fc["m_cw"], pi=float(fc["pi"]),
                        split=float(fc["split"]) if fc["comp"] in ("2Ц", "ОЦ") else None, mu=mus, kcm=float(fc["kcm"]), zc=int(fc["zc"])),
        combustor=dict(scheme=scheme, chamber=chamber),
        turbine=dict(type=fc["turb"], nt=int(fc["nt"]), cool=fc["cool"], mat=fc["mat_tw"], make=fc["m_tw"], hot_make=fc["m_hot"],
                     hot_mat=fc["mat_hot"], t41=float(fc["t41"]), lt=float(fc["lt"]), phit=float(fc["phit"]) or 0.5,
                     npt=int(fc["npt"]) if fc["rotors"] == 2 else "auto", kUpt=float(fc["kUpt"]) or 1.0, phipt=float(fc["phipt"]) or 0.5),
        hx=dict(recup=float(fc["eps_rec"]) if fc["recup"] == "печатный" else 0.0, recup_make="печать", recup_integral=int(fc["recup_integral"])),
        shaft=dict(drive="свободная СТ" if fc["rotors"] == 2 else "прямой"),
        build=dict(rotor_join=fc["rotor_join"], bear=fc["bear"], case=case, conn=fc["conn"], cores=int(fc["cores"]), comb_make=fc["m_comb"],
                   case_make=fc["m_case"], kD_case=float(fc["kD_case"])))


def make_seeds(kvs=("KV-3", "KV-1", "KV-6"), n_max=40, path=SEEDS_FILE):
    """замороженный файл опорных родителей из воронки: допустимые концепции отбора уровня 0 v3 (пробы l0_v3_probe; разные
    ниши воронки сначала, затем по оценке), переведённые в концепцию блоков и в геном. Кластер ядер - пропускается (нет в
    пространстве). Файл не меняется при правке воронки - отпечаток конфигурации прогона стабилен."""
    from blocks.concept import normalize
    FS = T.FS()
    out = {}
    for kv in kvs:
        rel, n = FUNNEL_V3_SELECT[kv]
        d = T.load_json(os.path.join(T.FUNNEL, rel))
        fe = [e for e in d["selected"] if e["feasible"]]
        fe.sort(key=lambda e: (-e["score"], e["idx"], e.get("case_arch") or ""))
        pick, niches = [], set()
        for e in fe:
            if e["niche"] not in niches:
                niches.add(e["niche"])
                pick.append(e)
        pick += [e for e in fe if e not in pick]
        rows = []
        for e in pick:
            if len(rows) >= n_max:
                break
            fc = FS.concept(int(e["idx"]), n, kv, int(d["seed"]))
            if e.get("case_arch"):
                fc = FS.apply_case(fc, e["case_arch"])
            if int(fc["cores"]) > 1:
                continue
            Cb = normalize(funnel_to_blocks(fc))
            x = encode(Cb)
            rows.append(dict(name=f"воронка v3 {kv} #{e['idx']} {e.get('case_arch') or ''}".strip(), niche=e["niche"], score=e["score"],
                             eta_el_funnel=e["eta_el"], genome=[round(float(v), 6) for v in x]))
        out[kv] = rows
    doc = dict(title="Опорные родители островов из воронки уровня 0 v3 (замороженный файл)", date="2026-11-15",
               source={kv: f"routes/gtu/funnel/{FUNNEL_V3_SELECT[kv][0]}" for kv in kvs}, n_max=n_max,
               rule="допустимые концепции отбора; разные ниши воронки, затем по оценке; кластер ядер - нет; перевод funnel_to_blocks -> encode "
                    "(семья изготовления - по способу колеса турбины, корпуса, горячей части): декодированная концепция близка к исходной, "
                    "но не обязана совпадать с ней (окна process_space), её допустимость проверяет прогон",
               space_hash=space_hash(), seeds=out)
    with open(path, "w", encoding="utf-8") as f:
        json.dump(doc, f, ensure_ascii=False, indent=1)
    return doc


def file_seeds(kv, path=SEEDS_FILE):
    try:
        d = json.load(open(path, encoding="utf-8"))
    except (OSError, ValueError):
        return []
    return [(r["name"], np.asarray(r["genome"], float)) for r in d.get("seeds", {}).get(kv, [])]


def reference_seeds(kv):
    """геномы опорных концепций блоков (reference.EXAMPLES и точка frame3) для этой точки - начальные родители островов;
    допустимость не предполагается: их оценивает прогон (бюджет оценок включает их)"""
    from blocks import reference as REF
    from blocks.concept import normalize
    out = []
    for name in REF.EXAMPLES:
        c = REF.example(name)
        if c.get("kv", "KV-3") != kv and not (kv != "KV-3" and c.get("kv") is None):
            continue
        if c.get("build", {}).get("cores", 1) > 1:
            continue
        Cn = normalize(dict(c, kv=kv))
        out.append((name, encode(Cn)))
    return out


def blocks_concept(C):
    """концепция без служебных полей (для blocks.evaluate)"""
    return {k: v for k, v in C.items() if not k.startswith("_")}


def summary(C):
    """короткая строка концепции для таблиц"""
    c, t, hx, f = C["compressor"], C["turbine"], C["hx"], C["fuel"]
    s = [f"{c['arch']}{c['n_ax'] if c['arch'] == 'ОЦ' else ''} π{c['pi']:.3g}" + (f" s{c['split']:.2f}" if c.get("split") else ""),
         f"{'осевая' if t['type'] == 'О' else 'радиальная'}×{t['nt']} {t['cool']} {t['mat']}", f"T41 {t['t41']:.0f}",
         f"{C['combustor']['scheme']} {C['combustor']['chamber']}", C["shaft"]["drive"], C["_make"]]
    if c["intercool"] != "нет":
        s.append("ПО")
    if hx["recup"]:
        s.append(f"Р {hx['recup']:.2f} {hx['recup_make']}")
    if hx["reheat"]:
        s.append(f"ПП x{hx['x_rh']:.2f} {hx['t_rh']:.0f}")
    if hx["ibc"]:
        s.append("IBC")
    if f["prep"] != "нет":
        s.append(f["prep"] + (f" S/C {f['sc']:.2f}" if f["prep"] == "паровой риформинг" else
                             (f" R {f['R_egr']:.2f}" if f["prep"] == "сухой риформинг" else "")))
    if C["build"]["bear"] == "магнит":
        s.append("магнит")
    if C["build"]["case"] == "оболочка":
        s.append("оболочка")
    return ", ".join(s)


# =========================================================================== генератор по индексу
def _sobol(kv, seed):
    from scipy.stats import qmc
    return qmc.Sobol(d=D, scramble=True, seed=int(seed) * 100 + KV_INDEX.get(kv, 0))


def genomes(start, n, kv="KV-3", seed=SEED):
    """геномы концепций [start, start+n) точки kv: строки последовательности Соболя (fast_forward)"""
    if n <= 0:
        return np.zeros((0, D))
    s = _sobol(kv, seed)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        if start:
            s.fast_forward(int(start))
        X = s.random(int(n))
    return np.clip(X, 0.0, 1.0 - 1e-12)


def genome_at(i, kv="KV-3", seed=SEED):
    return genomes(i, 1, kv, seed)[0]


def part_range(n, n_parts, part):
    """индексы части: [part·n//n_parts, (part+1)·n//n_parts) - разбиение ровно одно"""
    if not (1 <= n_parts <= max(n, 1)):
        raise ValueError(f"частей {n_parts} при n {n}")
    if not 0 <= part < n_parts:
        raise ValueError(f"часть {part} вне [0, {n_parts})")
    return part * n // n_parts, (part + 1) * n // n_parts


def space_hash():
    """отпечаток пространства (гены, окна, семьи, сетка) - входит в отпечаток конфигурации прогона"""
    s = json.dumps(dict(genes=GENES, ranges=RANGES, make=MAKE, desc=DESC_LABELS, seed=SEED), sort_keys=True, ensure_ascii=False)
    return hashlib.sha256(s.encode("utf-8")).hexdigest()[:16]


def describe():
    FS = T.FS()
    return dict(
        genes=GENES, D=D, ranges=RANGES, make=MAKE, pi_range={str(k): v for k, v in FS.PI_RANGE.items()},
        split_range=FS.SPLIT_RANGE, mu_range=FS.MU_RANGE, t41_min=FS.T41_MIN, mat_tw=FS.MAT_TW, mat_hot=FS.MAT_HOT,
        grid=dict(names=DESC_NAMES, bins=DESC_BINS, n_cells=N_CELLS, labels=DESC_LABELS, note=DESC_NOTE),
        generator=dict(kind="Sobol scramble (scipy.stats.qmc), сид SEED·100 + номер точки", seed=SEED,
                       part_rule="часть k из P: индексы [k·N//P, (k+1)·N//P)", scipy_pinned="1.18.1"),
        space_hash=space_hash(),
        kv_available=list(FS.KV),
    )


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--describe", action="store_true")
    ap.add_argument("--sample", type=int, default=0)
    ap.add_argument("--kv", default="KV-3")
    ap.add_argument("--start", type=int, default=0)
    ap.add_argument("--make-seeds", action="store_true", help="заморозить опорных родителей из воронки v3 -> process_seeds.json")
    a = ap.parse_args(argv)
    if a.make_seeds:
        d = make_seeds()
        print({kv: len(v) for kv, v in d["seeds"].items()}, "->", os.path.relpath(SEEDS_FILE, T.REPO))
    if a.describe:
        print(json.dumps(describe(), ensure_ascii=False, indent=1, default=str))
    for k, x in enumerate(genomes(a.start, a.sample, a.kv)):
        C = decode(x, a.kv)
        print(a.start + k, cell_label(descriptors(C)), "::", summary(C), C["_repairs"] or "")
    return 0


if __name__ == "__main__":
    sys.exit(main())
