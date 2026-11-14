"""Воронка двигателя-заменителя, уровень 0: число отдельных деталей, узлов, движущихся и охлаждаемых деталей.

Счёт - по ведомости уровня «деталь, которую отдельно изготавливают и отдельно ставят в двигатель». Покупное изделие
(подшипник скольжения, магнитный подшипник, свеча, форсунка, стартёр, датчик) считается одной деталью и одним узлом,
его внутренние детали не считаются. Сварной ротор - одна отдельная деталь (заготовок две, шов учитывается в welds).
Крепёж считается отдельно (fasteners) и в число деталей не входит. Каждое правило - в README узла, раздел
«parts_model», с основанием.

Балл:
    N_eff = parts + 2*assemblies + 0,1*fasteners + moving + cooled
    S     = 10 * N_HALF / (N_HALF + N_eff)
N_HALF - N_eff опорной традиционной схемы (REF_BOM: KV-3 по frame3 - два колеса, осевая турбина из трёх ступеней
с лопатками в замках, болтовые фланцы); у неё балл 5. Больше любых деталей - балл ниже (строго).
"""
import math

K_ASSEMBLY, K_FASTENER, K_MOVING, K_COOLED = 2.0, 0.1, 1.0, 1.0
BOLT_PITCH = 0.10          # м, шаг болтов фланца (оценка)
CLAMP_D_MAX = 1.20         # м, наибольший фланец под быстроразъёмный хомут (оценка)
SPLIT_D = 1.00             # м, корпус больше - с горизонтальным разъёмом (две половины; оценка)
TC_PER_CORE, TC_PER_ENGINE = 6, 2     # термопары выхлопа (оценка)


def _turb_rotor_parts(stages):
    """колёса турбины: цельное колесо (радиальное или блиск) - 1; диск с лопатками в замках - 1 + zt"""
    parts, blades = 0, 0
    for s in stages:
        if s["blisk"]:
            parts += 1
        else:
            parts += 1 + s["zt"]; blades += s["zt"]
    return parts, blades


def _ngv_parts(stages, cooled):
    """сопловые аппараты: цельное кольцо - 1; иначе сегменты-дуплеты zn/2"""
    n = 0
    for s in stages:
        n += 1 if s["ngv_ring"] else math.ceil(s["zn"] / 2)
    return n


def count(bom):
    """bom - ведомость концепции (README, «parts_model»). Возврат: parts, assemblies, moving, cooled, fasteners,
    welds, unique, n_eff, score, by_unit {узел: {parts, moving, cooled, fasteners, rule}}"""
    cores, rotors = bom["cores"], bom["rotors"]
    comp = bom["comp"]; turb = bom["turb"]; pt = bom.get("pt")
    join = bom["rotor_join"]; bear = bom["bear"]; conn = bom["conn"]
    case = bom["case"]; comb = bom["comb"]; rec = bom.get("recup") or dict(present=False)
    U = {}

    def add(unit, parts=0, moving=0, cooled=0, fasteners=0, assemblies=0, welds=0, rule=""):
        u = U.setdefault(unit, dict(parts=0, moving=0, cooled=0, fasteners=0, assemblies=0, welds=0, rule=[]))
        u["parts"] += parts; u["moving"] += moving; u["cooled"] += cooled; u["fasteners"] += fasteners
        u["assemblies"] += assemblies; u["welds"] += welds
        if rule:
            u["rule"].append(rule)

    # ------------------------------------------------ ротор газогенератора (на ядро)
    n_cw = len(comp)
    tparts, tblades = _turb_rotor_parts(turb["stages"])
    cooled_rotor = 0
    if bom.get("rotor_cooled"):
        cooled_rotor = sum(s["zt"] if not s["blisk"] else 1 for s in turb["stages"])
    if join == "заодно":
        rp, rw, rf = 1, 0, 0
        rule = "ротор заодно: колёса компрессора, турбины и вал - одна деталь"
    elif join == "сварка":
        rp, rw, rf = 1, 1, 0
        rule = "ротор из двух заготовок (секции компрессора и турбины), сварка трением: одна отдельная деталь, один шов"
    else:
        rp = n_cw + tparts + 1 + 1          # колёса + вал + стяжной болт
        rw, rf = 0, 1 + tblades             # гайка стяжки + замки лопаток
        rule = f"сборный ротор: {n_cw} колёс компрессора + колёса турбины ({tparts}, из них лопаток {tblades}) + вал + стяжной болт; крепёж - гайка и замки лопаток"
    if bear == "магнит":
        rp += 3; rule += "; магнитные подшипники: 2 радиальные втулки и упорный диск на роторе"
    elif join == "сборный":
        rp += 1; rule += "; упорный гребень - отдельная деталь"
    add("ротор", parts=rp * cores, moving=rp * cores, cooled=cooled_rotor * cores, fasteners=rf * cores,
        assemblies=cores, welds=rw * cores, rule=rule + (f"; на {cores} ядра" if cores > 1 else ""))
    add("ротор", parts=1, moving=1, rule="ступица выходной муфты - 1 на двигатель")

    # ------------------------------------------------ силовая турбина (двухвальная схема)
    if rotors == 2 and pt:
        pp, pb = _turb_rotor_parts(pt["stages"])
        n = pp + 1 + (1 if len(pt["stages"]) > 1 else 0) + (3 if bear == "магнит" else 1)
        add("силовая турбина", parts=n, moving=n, fasteners=pb + 1, assemblies=1,
            rule=f"колёса СТ ({pp}, лопаток {pb}) + вал + стяжка при 2 ступенях + упорный гребень или 3 детали магнитного подшипника")
        ngv = _ngv_parts(pt["stages"], False)
        add("силовая турбина", parts=ngv + 1 + len(pt["stages"]),
            rule="сопловые аппараты СТ + корпус СТ + кольцо над рабочим колесом на ступень")

    # ------------------------------------------------ корпуса (холодная часть, на ядро)
    m = case["method"]
    if m == "печать":
        cp = case.get("pieces", 1)
        crule = f"печать: вход, корпус и диффузоры одной деталью (частей по камере принтера: {cp})"
    else:
        cp = 1 + (1 if bom["inlet"] == "двух" else 0)
        for k, t in enumerate(comp):
            if t == "О":
                cp += 1                                   # кольцо статора осевой ступени
            elif m == "литьё":
                cp += 1                                   # корпус ступени с диффузором одной отливкой
            else:
                cp += 2                                   # корпус ступени + диффузор
            if t in ("Ц", "Д") and k > 0:
                cp += 1                                   # обратный канал перед второй ступенью
        crule = (f"{m}: вход{' + камера двустороннего входа' if bom['inlet'] == 'двух' else ''}; на центробежную ступень "
                 f"{'корпус с диффузором (1)' if m == 'литьё' else 'корпус + диффузор (2)'}; обратный канал перед второй ступенью; "
                 "статор осевой ступени - кольцо")
        if case["D"] > SPLIT_D:
            cp *= 2; crule += f"; корпус больше {SPLIT_D} м - с горизонтальным разъёмом (x2)"
    sumps = 0 if m == "печать" else 2
    add("корпуса", parts=(cp + sumps) * cores, assemblies=cores,
        rule=crule + f"; корпуса подшипников: {sumps} (у печатного корпуса - заодно)")
    add("корпуса", parts=1, rule="выхлопной диффузор и коллектор - 1 на двигатель")

    # ------------------------------------------------ камера сгорания (на ядро)
    kind, cm, ninj = comb["kind"], comb["method"], comb["n_inj"]
    hot_int = bom.get("hot_integral", False)
    if kind == "кольцевая":
        if cm == "печать":
            liner = 0 if hot_int else 1; inj = 1
            r = "печать: жаровая труба с фронтовым устройством и завихрителями - 1 (у цельной горячей части - в ней); форсунки с коллектором - 1"
        elif cm == "литьё":
            liner = 2; inj = ninj + 1
            r = f"литьё: жаровая труба + фронтовое устройство; {ninj} форсунок + коллектор"
        elif cm == "керамика":
            liner = 3; inj = ninj + 1
            r = f"CMC: внутренняя и наружная стенки + металлический фронт; {ninj} форсунок + коллектор; 4 податливые опоры (крепёж)"
        else:
            liner = 3 + ninj; inj = ninj + 1
            r = f"лист: две стенки + фронт + {ninj} завихрителей; {ninj} форсунок + коллектор"
        add("камера", parts=(liner + inj + 2) * cores, fasteners=(4 if cm == "керамика" else 0) * cores, assemblies=cores,
            rule=r + "; 2 свечи")
    elif kind == "выносная":
        n = (2 if cm == "печать" else 3) + 1 + 2
        add("камера", parts=n * cores, assemblies=cores,
            rule="выносная: жаровая труба с переходником (печать - 1, иначе 2) + корпус камеры + 1 форсунка + 2 свечи")
    elif kind == "RDC":
        add("камера", parts=4 * cores, assemblies=cores,
            rule="вращающаяся детонация: кольцевой канал + печатная форсуночная плита + преддетонатор + рубашка охлаждения")
    elif kind == "ВР":
        add("камера", parts=(3 + 4) * cores, moving=2 * cores, assemblies=2 * cores,
            rule="волновой ротор: барабан (движется) + 2 торцевые плиты + привод барабана (движется); камера - выносная: жаровая труба + корпус + форсунка + свеча")
    else:
        raise ValueError(kind)

    # ------------------------------------------------ турбина: статор (на ядро)
    cooled_ngv = bom.get("ngv_cooled", False)
    if hot_int:
        ngv_n = 1
        add("горячая часть", parts=cores, cooled=(1 if cooled_ngv else 0) * cores, assemblies=cores,
            rule="цельная горячая часть: сопловые аппараты + корпус турбины (+ жаровая труба кольцевой печатной камеры) - 1")
    else:
        ngv_n = _ngv_parts(turb["stages"], cooled_ngv)
        shrouds = sum(1 for _ in turb["stages"]) if turb["type"] == "О" else 0
        add("горячая часть", parts=(ngv_n + 1 + shrouds) * cores, cooled=(ngv_n if cooled_ngv else 0) * cores, assemblies=cores,
            rule="сопловой аппарат: цельное кольцо - 1, иначе сегменты-дуплеты zn/2; корпус турбины - 1; кольцо над колесом осевой ступени - 1")

    # ------------------------------------------------ подшипники и уплотнения (на ротор)
    n_rot = cores + (1 if rotors == 2 else 0)
    if bear == "масло":
        add("подшипники", parts=5 * n_rot, assemblies=3 * n_rot,
            rule="масло: 2 опорных + 1 упорный (покупные) + 2 масляных уплотнения на ротор")
    else:
        add("подшипники", parts=7 * n_rot, assemblies=3 * n_rot,
            rule="магнитные: 2 радиальных + 1 осевой + 2 страховочных подшипника + 2 кольца датчиков положения на ротор")
    seals = (n_cw if m == "обработка" else 0) + 1 + max(0, len(turb["stages"]) - 1)
    add("уплотнения", parts=seals * cores,
        rule="лабиринт покрывного диска на колесо (у литого и печатного корпуса - заодно) + разгрузочный поршень + межступенные турбины")

    # ------------------------------------------------ маслосистема двигателя, топливо, пуск
    if bear == "масло":
        add("маслосистема", parts=4 * n_rot + 2, assemblies=1,
            rule="на ротор: 2 подводящих + 2 сливных трубки; на двигатель: коллектор + суфлёр; насос и бак - пакета")
    add("топливо", parts=2 * cores + (1 if cores > 1 else 0), assemblies=1,
        rule="подвод газа + гибкий рукав на ядро; распределитель при кластере; регулирующий клапан - пакета")
    if bom.get("new_starters", 0):
        add("пуск", parts=2 * bom["new_starters"], moving=bom["new_starters"], assemblies=bom["new_starters"],
            rule="новый стартёр + переходник на каждый ротор, который стартёр пакета не проворачивает")

    # ------------------------------------------------ датчики и САУ
    sens = TC_PER_CORE * cores + TC_PER_ENGINE + 2 * n_rot + (4 * n_rot if bear == "масло" else 0) + 3 * cores + cores
    add("датчики", parts=sens + (1 if bear == "магнит" else 0), assemblies=1 + (1 if bear == "магнит" else 0),
        rule="термопары выхлопа, 2 датчика частоты на ротор, 4 датчика вибрации на ротор с маслом (у магнитных - встроены), "
             "3 датчика давления и жгут на ядро; контроллер магнитных подшипников")

    # ------------------------------------------------ рекуператор
    if rec.get("present"):
        mods = rec.get("modules", 1)
        n = mods + 1 + (0 if rec.get("integral") else 4)
        add("рекуператор", parts=n, assemblies=1,
            rule=f"печатные модули ({mods}) + клапан обвода; вынесенный - ещё 2 короба и 2 компенсатора")

    # ------------------------------------------------ крепёж фланцев (отдельно)
    static_chain = U["корпуса"]["parts"] // max(cores, 1) + (1 if not hot_int else 0) + 1
    joints = max(1, static_chain - 1) + (2 if rec.get("present") else 0)
    D = case["D"]
    if conn == "быстроразъём" and D <= CLAMP_D_MAX:
        per = 2; jr = f"быстроразъёмный хомут с замком - 2 на стык (фланец до {CLAMP_D_MAX} м)"
    else:
        per = math.ceil(math.pi * D / BOLT_PITCH); jr = f"болты фланца: π·D/{BOLT_PITCH} м = {per} на стык"
    add("крепёж фланцев", fasteners=joints * per * cores, rule=f"стыков корпусов на ядро: {joints}; {jr}")

    tot = {k: sum(u[k] for u in U.values()) for k in ("parts", "moving", "cooled", "fasteners", "assemblies", "welds")}
    n_eff = tot["parts"] + K_ASSEMBLY * tot["assemblies"] + K_FASTENER * tot["fasteners"] + K_MOVING * tot["moving"] + K_COOLED * tot["cooled"]
    unique = _unique(bom, U)
    by = {k: dict(v, rule="; ".join(v["rule"])) for k, v in U.items()}
    return dict(parts=tot["parts"], assemblies=tot["assemblies"], moving=tot["moving"], cooled=tot["cooled"],
                fasteners=tot["fasteners"], welds=tot["welds"], unique=unique, n_eff=round(n_eff, 2),
                score=round(score_from_neff(n_eff), 3), by_unit=by)


def _unique(bom, U):
    """число разных наименований: у кластера ядра одинаковые - детали ядра считаются один раз (оценка)"""
    cores = bom["cores"]
    per_core = ("ротор", "корпуса", "камера", "горячая часть", "уплотнения")
    n = 0
    for k, u in U.items():
        n += math.ceil(u["parts"] / cores) if (k in per_core and cores > 1) else u["parts"]
    return n


def score_from_neff(n_eff):
    return 10.0 * N_HALF / (N_HALF + max(n_eff, 0.0))


# опорная ведомость: KV-3, традиционная сборная схема по frame3 (2 колеса, двусторонний вход, осевая турбина
# 3 ступени, Dm 0,537 м > 0,50 м - диски с лопатками, класс U, кольцевая камера из листа, литые корпуса, болты);
# числа лопаток - середина диапазонов воронки (zt 55, zn 40; оценка), форсунок 14, D корпуса 1,9*0,703+0,5
REF_BOM = dict(cores=1, rotors=1, comp=["Ц", "Ц"], rotor_join="сборный", bear="масло", conn="болты", inlet="двух",
               turb=dict(type="О", stages=[dict(blisk=False, zt=55, zn=40, ngv_ring=False)] * 3),
               pt=None, case=dict(method="литьё", D=1.836), comb=dict(kind="кольцевая", method="обработка", n_inj=14),
               hot_integral=False, ngv_cooled=False, rotor_cooled=False, recup=dict(present=False), new_starters=0)
N_HALF = 1.0
N_HALF = count(REF_BOM)["n_eff"]
