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

v3 (решение владельца 15.11.2026, п. 3) - ось «корпус» bom["case_arch"]:
    «картридж» - корпус с разъёмом, газогенератор вынимается целиком; «разъём» - то же без картриджа (ротор вынимается отдельно,
    опорная точка frame3); «оболочка» - цельная печатная оболочка (вход, корпуса и диффузоры, корпуса подшипников, каналы масла
    и топлива в стенке) с ротором, собранным внутри, горячая часть - одной печатной деталью (разбор живой геометрии:
    routes/gtu/organic/structure/README.md, п. 4; organic/print_review: shell_assembly).
Правило собираемости (там же, решение владельца): число частей статора не меньше числа статорных элементов, запирающих колёса.
Запирающие элементы ротора из n_cw колёс компрессора и n_t ступеней турбины: передний (вход с корпусом подшипника), между колёсами
компрессора (n_cw - 1), СА 1 с внутренним корпусом камеры, СА 2 ... СА n_t, задний (выхлоп с корпусом подшипника) - всего
n_cw + n_t + 1. Недостающие части добавляются (trap): у печатной оболочки одна её часть - один элемент; внешний печатный корпус
картриджа элементом не служит (картридж выходит из него); цельная горячая часть - один элемент (СА 1). Ротор одной деталью
(заодно, сварка) - каждый элемент между колёсами - из двух половин (цельное кольцо между колёсами одной детали не поставить;
оценка): +1 деталь на элемент, если он ещё не разъёмный. Добавленные части холодного корпуса - новые стыки (как прочие части
корпуса), кольца СА - внутренние, без стыков. Проверка: REF_BOM как оболочка - 211 деталей, крепёж 284 (хомут) / 340 (болты),
-108 против 319 (organic/print_review_results.json -> shell_assembly).
Поправки скептика оценки (15.11.2026): половины цельной горячей части у ротора одной деталью стыкуются осевым разъёмом - 2 фланца
по ceil(L_hot/шаг) болтов (bom["hot"], оценка); картридж (одно ядро) - 2 стыка крепления в корпусе Ø0,5·D, те же, что у модуля
замены funnel_space.model_inputs (оценка); горизонтальный разъём корпуса (D > SPLIT_D) - половины не удваивают стыки по окружности,
разъём - 2 продольных фланца по длине корпуса case["L"]. Освобождение одноступенчатой радиальной турбины от половин отклонено: между
колёсами ротора одной детали всегда есть стенка с отверстием меньше обоих колёс (trap_rule).
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
    case_arch = bom.get("case_arch") or "разъём"
    shell = case_arch == "оболочка"
    hot_int = bom.get("hot_integral", False)
    cp_chain = None
    if m == "печать":
        cp = case.get("pieces", 1)
        crule = (f"печатная оболочка: вход, корпус, диффузоры и корпуса подшипников одной деталью (частей по камере принтера: {cp})" if shell else
                 f"печать: вход, корпус и диффузоры одной деталью (частей по камере принтера: {cp})")
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
        cp_chain = cp                                     # стыков по окружности столько же, сколько частей по длине
        if case["D"] > SPLIT_D:
            cp *= 2; crule += f"; корпус больше {SPLIT_D} м - с горизонтальным разъёмом (x2)"
    sumps = 0 if m == "печать" else 2
    trap = trap_rule(bom, cp)
    add("корпуса", parts=(cp + sumps + trap["add_cold"]) * cores, assemblies=cores,
        rule=crule + f"; корпуса подшипников: {sumps} (у печатного корпуса - заодно)" + (f"; {trap['rule_cold']}" if trap["add_cold"] else ""))
    add("корпуса", parts=1, rule="выхлопной диффузор и коллектор - 1 на двигатель")

    # ------------------------------------------------ камера сгорания (на ядро)
    kind, cm, ninj = comb["kind"], comb["method"], comb["n_inj"]
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
        halves = trap["one_piece"] and trap["add_hot"] > 0 and bom.get("hot")
        hf = 2 * math.ceil(bom["hot"]["L"] / BOLT_PITCH) if halves else 0     # осевой разъём половин корпуса турбины: 2 фланца (оценка)
        add("горячая часть", parts=(1 + trap["add_hot"]) * cores, cooled=(1 + trap["add_hot"] if cooled_ngv else 0) * cores, assemblies=cores,
            fasteners=hf * cores,
            rule="цельная горячая часть: сопловые аппараты + корпус турбины (+ жаровая труба кольцевой печатной камеры) - 1"
                 + (f"; {trap['rule_hot']}" if trap["add_hot"] else "")
                 + (f"; половины цельной горячей части - осевой разъём, 2 фланца по {hf // 2} болтов" if hf else ""))
    else:
        ngv_n = _ngv_parts(turb["stages"], cooled_ngv)
        shrouds = sum(1 for _ in turb["stages"]) if turb["type"] == "О" else 0
        add("горячая часть", parts=(ngv_n + 1 + shrouds + trap["add_hot"]) * cores,
            cooled=(ngv_n + trap["add_hot"] if cooled_ngv else 0) * cores, assemblies=cores,
            rule="сопловой аппарат: цельное кольцо - 1, иначе сегменты-дуплеты zn/2; корпус турбины - 1; кольцо над колесом осевой ступени - 1"
                 + (f"; {trap['rule_hot']}" if trap["add_hot"] else ""))

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
        if shell:
            n_pt = 1 if rotors == 2 else 0
            add("маслосистема", parts=1 + 4 * n_pt, assemblies=1,
                rule="оболочка: подвод, слив и вентиляция газогенератора - каналы в стенках стоек оболочки; суфлёр - 1"
                     + ("; силовая турбина вне оболочки: 2 подводящих + 2 сливных трубки" if n_pt else "")
                     + " (organic/structure shell_topo: встроенные каналы)")
        else:
            add("маслосистема", parts=4 * n_rot + 2, assemblies=1,
                rule="на ротор: 2 подводящих + 2 сливных трубки; на двигатель: коллектор + суфлёр; насос и бак - пакета")
    if shell:
        add("топливо", parts=cores + (1 if cores > 1 else 0), assemblies=1,
            rule="оболочка: подвод газа - канал в стенке к встроенному коллектору; гибкий рукав пакета - 1 на ядро; распределитель при кластере")
    else:
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
    # горизонтальный разъём (половины корпуса) не удваивает стыки по окружности: у стыка двух половин тот же π·D/шаг болтов; разъём -
    # два продольных фланца по длине корпуса (case["L"]; у ведомости без длины - не считаются). Прежде половины удваивали стыки
    # по окружности: у REF_BOM 978 единиц крепежа вместо 688 + продольных (поправка скептика оценки 15.11.2026)
    split_x = cp_chain is not None and case["D"] > SPLIT_D
    static_chain = (U["корпуса"]["parts"] // max(cores, 1) - (cp - cp_chain if split_x else 0)) + (1 if not hot_int else 0) + 1
    joints = max(1, static_chain - 1) + (2 if rec.get("present") else 0)
    D = case["D"]
    if conn == "быстроразъём" and D <= CLAMP_D_MAX:
        per = 2; jr = f"быстроразъёмный хомут с замком - 2 на стык (фланец до {CLAMP_D_MAX} м)"
    else:
        per = math.ceil(math.pi * D / BOLT_PITCH); jr = f"болты фланца: π·D/{BOLT_PITCH} м = {per} на стык"
    add("крепёж фланцев", fasteners=joints * per * cores, rule=f"стыков корпусов на ядро: {joints}; {jr}")
    if split_x and case.get("L"):
        nl = 2 * math.ceil(case["L"] / BOLT_PITCH)
        add("крепёж фланцев", fasteners=nl * cores, rule=f"горизонтальный разъём корпуса: 2 продольных фланца по длине {case['L']:.2f} м - {nl}")
    if case_arch == "картридж" and cores == 1:
        # картридж выходит из корпуса целиком: 2 стыка крепления газогенератора в корпусе Ø0,5·D - те же, что у модуля замены
        # funnel_space.model_inputs (joints=[J(0,5·D_case)]×2); без них картридж получал выигрыш монтажа без крепежа (оценка)
        Dj = 0.5 * D
        pj = 2 if (conn == "быстроразъём" and Dj <= CLAMP_D_MAX) else math.ceil(math.pi * Dj / BOLT_PITCH)
        add("крепёж фланцев", fasteners=2 * pj, rule=f"картридж: 2 стыка крепления в корпусе (Ø{Dj:.2f} м) по {pj}")

    tot = {k: sum(u[k] for u in U.values()) for k in ("parts", "moving", "cooled", "fasteners", "assemblies", "welds")}
    n_eff = tot["parts"] + K_ASSEMBLY * tot["assemblies"] + K_FASTENER * tot["fasteners"] + K_MOVING * tot["moving"] + K_COOLED * tot["cooled"]
    unique = _unique(bom, U)
    by = {k: dict(v, rule="; ".join(v["rule"])) for k, v in U.items()}
    return dict(parts=tot["parts"], assemblies=tot["assemblies"], moving=tot["moving"], cooled=tot["cooled"],
                fasteners=tot["fasteners"], welds=tot["welds"], unique=unique, n_eff=round(n_eff, 2),
                score=round(score_from_neff(n_eff), 3), by_unit=by, case_arch=case_arch,
                trap={k: v for k, v in trap.items() if not k.startswith("rule")})


def trap_rule(bom, cp):
    """правило собираемости (решение владельца 15.11.2026): сколько частей статора добавить на ядро, чтобы число частей было не меньше
    числа запирающих колёса элементов (модуль - docstring). cp - частей холодного корпуса по правилам способа. Возврат: n_trap,
    cover (чем закрыты), add_cold, add_hot, split (элементы из половин у ротора одной деталью), rule_cold, rule_hot"""
    comp, turb, case = bom["comp"], bom["turb"], bom["case"]
    n_cw, n_t = len(comp), len(turb["stages"])
    # оболочка; ведомость без оси (внешние вызовы прежнего вида: organic/structure/shell_topo.py) с печатным корпусом - тоже оболочка
    shell = bom.get("case_arch") == "оболочка" or (bom.get("case_arch") is None and case["method"] == "печать")
    hot_int = bool(bom.get("hot_integral"))
    one_piece = bom["rotor_join"] in ("заодно", "сварка")
    m = case["method"]
    if m == "печать":
        cold_cover = cp if shell else 0          # часть оболочки - один элемент; внешний печатный корпус картриджа не запирает
    else:
        cold_cover = n_cw                        # вход + обратные каналы и кольца статора осевых ступеней - отдельными деталями
    hot_cover = 1 if hot_int else n_t
    split = 0
    if one_piece:
        add_c = 0 if cold_cover >= 1 else 1                      # передний элемент
        if m == "печать":
            add_c += 2 * (n_cw - 1)                              # элементы между колёсами - половинами, оболочкой быть не могут
            split += n_cw - 1
        elif case["D"] <= SPLIT_D:
            add_c += n_cw - 1                                    # цельные кольца между колёсами одной детали - в половины
            split += n_cw - 1
        # одноступенчатая радиальная турбина - то же правило: сопловое кольцо снаружи колеса не запирает, но между колесом компрессора и
        # колесом турбины ротора одной детали стоит стенка с отверстием меньше обоих колёс (внутренний корпус камеры, задняя стенка
        # компрессора; вал - 0,2·Dc), и одна часть - из половин. Аналогия с корпусом турбокомпрессора не переносится: там колесо
        # компрессора снимается (гайка), ротор не одной деталью (поправка скептика оценки отклонена 15.11.2026, оценка)
        if hot_int:
            add_h = 1 + 2 * (n_t - 1)                            # цельная горячая часть (СА 1) - из двух половин; СА 2.. - половинами
            split += n_t
        else:
            add_h = sum(1 for s in turb["stages"] if s.get("ngv_ring"))   # цельные кольца СА - в половины (сегменты и так разъёмны)
            split += add_h
    else:
        add_c = max(0, n_cw - cold_cover)
        add_h = max(0, n_t - hot_cover)
    n_trap = n_cw + n_t + 1
    rc = (f"правило собираемости: запирающих элементов {n_trap} (передний, {n_cw - 1} между колёсами компрессора, СА 1-{n_t}, задний); "
          f"холодный корпус даёт {min(cold_cover, n_cw)} из {n_cw}" + (", ротор одной деталью - элементы между колёсами из половин" if one_piece else "")
          + f": +{add_c}")
    rh = (f"правило собираемости: СА {n_t} ступеней, горячая часть даёт {hot_cover}"
          + (", ротор одной деталью - цельные кольца СА из половин" if one_piece else "") + f": +{add_h} (кольца СА - внутренние, без стыков)")
    return dict(n_trap=n_trap, cover_cold=min(cold_cover, n_cw), cover_hot=hot_cover, add_cold=add_c, add_hot=add_h, split=split,
                one_piece=one_piece, rule_cold=rc, rule_hot=rh)


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
# числа лопаток - середина диапазонов воронки (zt 55, zn 40; оценка), форсунок 14, D корпуса 1,9*0,703+0,5;
# длина корпуса - формула воронки funnel_space.derive: L = 0,6·L_rot, L_rot = 0,35·(D2_1 + D2_2) + 0,15·Dm·3 + 0,5·Dc по
# frame3_point_results (D2 0,695 и 0,703 м, Dm 0,537 м) - нужна для продольных фланцев горизонтального разъёма (D > SPLIT_D)
REF_L_CASE = round(0.6 * (0.35 * (0.694575 + 0.702541) + 0.15 * 0.536549 * 3 + 0.5 * 0.702541), 3)
REF_BOM = dict(cores=1, rotors=1, comp=["Ц", "Ц"], rotor_join="сборный", bear="масло", conn="болты", inlet="двух",
               turb=dict(type="О", stages=[dict(blisk=False, zt=55, zn=40, ngv_ring=False)] * 3),
               pt=None, case=dict(method="литьё", D=1.836, L=REF_L_CASE), comb=dict(kind="кольцевая", method="обработка", n_inj=14),
               hot_integral=False, ngv_cooled=False, rotor_cooled=False, recup=dict(present=False), new_starters=0)
N_HALF = 1.0
N_HALF = count(REF_BOM)["n_eff"]
