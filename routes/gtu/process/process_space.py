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
    dict(id="t41", kind="cont", block="turbine", what="T41: от 850 °C (funnel_space.T41_MIN) до предела класса и сплава (make_model.t41_limit)"),
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
    phit=(0.45, 0.80), kUpt=(0.85, 1.00), phipt=(0.45, 0.80), kD=(1.7, 2.1), trh_below=150.0,
    basis="окна непрерывных осей воронки (funnel_space.CONT: kcm, lt, phit, kUpt, phipt, kD_case), блока hx (ε 0,70-0,92), "
          "перебора циклов (organic/cycles A: x_rh 0,35-0,65, sc 1,5-3, R_egr 0,5-4); trh - на 150 K ниже T41 до предела класса (оценка)",
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
    t41 = round(_lin(v["t41"], FS.T41_MIN, max(tl, FS.T41_MIN)), 1)
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
            pass      # одна ступень: правило ic_stages (оценка отклонит) - не чинится, ниша «ПО» с 1Ц/Д недостижима
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
    FS, MK = T.FS(), T.MK()
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
    a = ap.parse_args(argv)
    if a.describe:
        print(json.dumps(describe(), ensure_ascii=False, indent=1, default=str))
    for k, x in enumerate(genomes(a.start, a.sample, a.kv)):
        C = decode(x, a.kv)
        print(a.start + k, cell_label(descriptors(C)), "::", summary(C), C["_repairs"] or "")
    return 0


if __name__ == "__main__":
    sys.exit(main())
