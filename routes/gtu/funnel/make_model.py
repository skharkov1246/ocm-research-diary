"""Воронка двигателя-заменителя, уровень 0: изготовимость концепции -> балл 0-10.

Модель отвечает на вопросы: влезает ли деталь в камеру принтера или в предел отливки и поковки, печатается ли
сплав при своей температуре, нужен ли ГИП, сколько недель до первого комплекта и какова готовность технологии (TRL).
Здесь же лежат таблицы фактов, общие для воронки (принтеры, пределы литья и поковки, сплавы, TRL): у каждого факта
есть основание и статус - «открыт» (документ прочитан), «по выдержке поиска» (число из выдержки поисковой системы,
документ не открывался), «принято в треке» (README трека, LIM в gtu_explore.py) или «оценка» (своя, не проверена).

Балл (README узла, раздел «make_model»):
    S = 10 * (0,30*s_fit + 0,25*s_proc + 0,30*s_trl + 0,15*s_time)
    s_fit  - доля деталей, которые изготавливаются выбранным способом без деления на части;
    s_proc - худшее из качеств способа по деталям (печатаемость или литейность сплава, 0,2...1,0);
    s_trl  - (TRL_min - 1) / 8, TRL_min - наименьшая готовность среди технологий концепции;
    s_time = 26 / (26 + недели до первого комплекта).
Стоимость - в условных единицах, только для ранжирования, в балл не входит (цен в модели нет).
"""
import math

DATE = "2026-11-14"

# --------------------------------------------------------------------------- принтеры (камера, мм)
# скорость построения: у BLT-S1500 и FS1521M - из выдержки поиска; у остальных - оценка пропорционально
# числу лазеров по BLT-S1500 (900 см3/ч на 26 лазеров = 34,6 см3/ч на лазер).
_PER_LASER = 900.0 / 26.0
PRINTERS = [
    dict(name="EOS M 400", proc="LPBF", country="DE", box=(400, 400, 400), lasers=1, rate=_PER_LASER * 1,
         rate_basis="оценка: 34,6 см3/ч на лазер по BLT-S1500",
         src="https://www.eos.info/en-us/metal-solutions/metal-printers/eos-m-400", status="по выдержке поиска"),
    dict(name="Nikon SLM NXG 600E", proc="LPBF", country="DE", box=(600, 600, 1500), lasers=12, rate=_PER_LASER * 12,
         rate_basis="оценка: 34,6 см3/ч на лазер по BLT-S1500",
         src="https://nikon.com/company/technology/rd/pdf/nrr_vol6_2024_01_e.pdf", status="по выдержке поиска"),
    dict(name="BLT-S1500", proc="LPBF", country="CN", box=(1500, 1500, 1200), lasers=26, rate=900.0,
         rate_basis="по выдержке поиска (презентация BLT)",
         src="https://3dprintingindustry.com/news/blt-exhibits-new-metal-am-solutions-on-formnext-2023-show-floor-226468/",
         status="по выдержке поиска",
         note="TCT Magazine даёт 1500x1500x1500 мм; взята меньшая высота 1200 мм из материалов изготовителя"),
    dict(name="Farsoon FS1521M-U", proc="LPBF", country="CN", round_d=1530, z=1650, lasers=16, rate=400.0,
         rate_basis="по выдержке поиска (до 400 см3/ч у FS1521M)",
         src="https://www.tctmagazine.com/farsoon-launches-16-laser-fs1521m-metal-3d-printer-at-tct-asia-2023/",
         status="по выдержке поиска", note="круглая камера Ø1530 мм; Z 850 мм у FS1521M, 1650 мм у FS1521M-U"),
    dict(name="Eplus3D EP-M2050", proc="LPBF", country="CN", box=(2050, 2050, 1100), lasers=36, rate=_PER_LASER * 36,
         rate_basis="оценка: 34,6 см3/ч на лазер по BLT-S1500",
         src="https://www.metal-am.com/eplus3ds-new-ep-m2050-offers-64-laser-option-for-high-speed-additive-manufacturing/",
         status="по выдержке поиска"),
    dict(name="ExOne / Desktop Metal X1 160Pro", proc="BJ", country="US", box=(800, 500, 400), lasers=0, rate=10000.0,
         rate_basis="по выдержке поиска («более 10 000 см3/ч в зависимости от материала»)",
         src="https://www.exone.com/Admin/ExOne/media/Metal-Systems/211112_X1_Datasheet_X1_160Pro_EN_2021-11.pdf",
         status="по выдержке поиска", note="после спекания деталь меньше заготовки: принят линейный запас 0,85 (оценка)"),
]
BJ_SHRINK = 0.85          # оценка: заготовка binder jetting больше детали на усадку спекания
BJ_ROTATING = False       # оценка: детали binder jetting с остаточной пористостью для ротора не принимаются
DED = dict(proc="DED", max_m=3.0, rate=300.0, status="оценка",
           note="наплавка роботом без камеры: размер ограничен манипулятором (принято 3 м), только статорные детали "
                "с последующей обработкой; скорость 300 см3/ч - оценка")

# --------------------------------------------------------------------------- пределы литья и поковки, м
CAST = {
    "rotor_mono": dict(ok=0.45, max=0.60, basis="LIM D_mono_ok / D_mono_max (gtu_explore.py), принято в треке"),
    "rotor_section": dict(ok=0.50, max=0.70, basis="как литое радиальное колесо: LIM D_rad_ok / D_rad_max, принято в треке"),
    "turb_radial": dict(ok=0.50, max=0.70, basis="LIM D_rad_ok / D_rad_max, принято в треке; один завод КНР заявляет Ø500 мм и 50 кг "
                                                 "(навык casting-spec, по выдержке поиска)"),
    "turb_axial_blisk": dict(ok=0.50, max=0.50, basis="LIM D_axcast: цельнолитое осевое колесо до 0,50 м, больше - диск + лопатки"),
    "wheel_ti": dict(ok=0.60, max=0.60, basis="литые тела вращения Ti-6Al-4V до Ø600 мм (JMST 2016, навык casting-spec, по выдержке поиска)"),
    "wheel_al": dict(ok=0.60, max=0.60, basis="оценка: принято как у титана, данных по колёсам ZL114A нет (навык casting-spec)"),
    "static": dict(ok=1.20, max=1.50, basis="заявления производств КНР о конструкционных отливках 1,5-2 м (навык casting-spec, по выдержке поиска)"),
}
FORGE = {
    "wheel_ti": dict(max=0.90, basis="LIM D_ti_max и D2_MAX gtu_family.py: наибольшая поковка титанового колеса, принято в треке"),
    "wheel_al": dict(max=0.90, basis="оценка: принято как у титана"),
    "shaft": dict(max=3.0, basis="оценка: вал из поковки, ограничений на этом уровне нет"),
    "disc": dict(max=1.2, basis="оценка: диск осевой турбины из поковки"),
    "static": dict(max=3.0, basis="оценка: обработка из поковки, проката или листа"),
}
CERAMIC = {
    "turb_radial": dict(max=0.35, basis="оценка, не проверено: керамические колёса освоены в машинах класса 300 кВт (CGT302, "
                                        "по выдержке поиска); крупнее - нет данных"),
    "static_hot": dict(max=1.5, basis="оценка: как крупная статорная отливка"),
}

# --------------------------------------------------------------------------- сплавы
# U_*: предел окружной скорости колеса компрессора, м/с; dT_cap: сдвиг температурного предела неохлаждаемого
# колеса турбины относительно литого K418 (LIM T_rad_max 950 C), K; print: печатаемость 1,0 / 0,6 / 0,2.
MAT = {
    "Ti64": dict(name="Ti-6Al-4V (литьё - ZTC4 + ГИП)", rho=4430.0, T_max=427.0, print=1.0, cast=0.7,
                 basis="плотность и Tmax данных - platform/fea/materials.json и gtu_explore.finalize; литьё ZTC4 + ГИП - навык "
                       "casting-spec (тонкие сечения трудны, по выдержке поиска); печатаемость - оценка (массовый сплав LPBF)"),
    "Al": dict(name="алюминиевый сплав (литьё ZL114A, печать AlSi10Mg)", rho=2800.0, T_max=200.0, print=1.0, cast=0.8,
               basis="плотность - gtu_explore.finalize; T за колесом <= 200 C - LIM T_al; печатаемость - оценка"),
    "Ni": dict(name="никелевый сплав колеса компрессора заодно с турбиной", rho=8000.0, T_max=None, print=None, cast=None,
               basis="свойства берутся по сплаву турбины (mat_tw); плотность 8000 - gtu_explore.finalize"),
    "K418": dict(name="K418 (≈ IN713C), литьё равноосное", rho=8000.0, dT_cap=0.0, print=0.2, cast=0.8,
                 basis="предел 950 C - LIM T_rad_max; горячие трещины крупных колёс - навык casting-spec (по выдержке поиска); "
                       "печать - оценка по аналогии с CM247LC (высокая сумма Al+Ti)"),
    "CM247LC": dict(name="CM247LC / MAR-M247", rho=8500.0, dT_cap=25.0, print=0.2, cast=0.8,
                    basis="печать LPBF: сильное растрескивание, ГИП помогает частично (Birmingham 2019, Chalmers 2026 - по выдержке "
                          "поиска: research.birmingham.ac.uk/files/60415383; research.chalmers.se/publication/553916); "
                          "мелкозернистое литьё колёс MAR-M247 - навык casting-spec; сдвиг +25 K - оценка"),
    "Haynes282": dict(name="Haynes 282, печать", rho=8270.0, dT_cap=-100.0, print=0.6, cast=None,
                      basis="оценка, не проверено: печатаемость средняя, температурный предел ниже литого K418 на 100 K"),
    "IN718": dict(name="IN718, печать + ГИП", rho=8190.0, dT_cap=-225.0, print=1.0, cast=None,
                  basis="оценка, не проверено: хорошо печатается; рабочая температура металла около 650 C, для неохлаждаемой "
                        "турбины - газ до 725 C"),
    "SC": dict(name="монокристаллические лопатки (литьё)", rho=8700.0, dT_cap=None, print=None, cast=0.6,
               basis="T41 1050 C - LIM T_ax_max; поставщики, сплавы и контроль в КНР не изучены (навык casting-spec)"),
    "Si3N4": dict(name="нитрид кремния, спекание + обработка", rho=3200.0, dT_cap=None, print=None, cast=None, ceramic=0.5,
                  basis="CGT302 (Kawasaki, Kyocera, 300 кВт): КПД 37,0 % при 1282 C (1997), цель 42 % при 1350 C; рецессия Si3N4 "
                        "в продуктах сгорания (по выдержке поиска: cir.nii.ac.jp/crid/1870302167595063808, "
                        "mrs-j.org/pub/tmrsj/vol07/vol07_068.pdf); плотность и качество 0,5 - оценка"),
    "K438": dict(name="K438 (≈ IN738LC), литьё", rho=8200.0, print=0.2, cast=0.8,
                 basis="лопатки около 900 C, стоек к горячей коррозии - навык casting-spec (по выдержке поиска); печать IN738LC - "
                       "трещины, ГИП не залечивает шире 6 мкм (по выдержке поиска)"),
    "HastX": dict(name="Hastelloy X, печать", rho=8220.0, print=1.0, cast=None,
                  basis="оценка, не проверено: свариваемый жаростойкий сплав камер сгорания, хорошо печатается"),
    "IN625": dict(name="IN625, печать", rho=8440.0, print=1.0, cast=None,
                  basis="печатные теплообменники IN625 (Argonne binder jetting; Thar Energy, рекуператор sCO2) - по выдержке поиска; "
                        "стандартная термообработка не подходит для LPBF-металла (по выдержке поиска)"),
    "CMC": dict(name="керамический композит", rho=2700.0, print=None, cast=None, ceramic=0.6,
                basis="оценка, не проверено"),
    "steel": dict(name="литая или кованая сталь", rho=7850.0, print=None, cast=1.0, basis="оценка"),
}

# --------------------------------------------------------------------------- готовность технологий
TRL = {
    "традиция": dict(trl=9, basis="серийные ГТУ линейки (паспорта fleet/), литьё и поковка колёс, лопатки в замках", status="принято в треке"),
    "лопатки_SC": dict(trl=9, basis="серийные лопатки; для поставки в КНР - не изучено (casting-spec)", status="оценка"),
    "печать_статика": dict(trl=7, basis="печатные статорные детали горячего тракта", status="оценка"),
    "печать_колесо_компрессора": dict(trl=5, basis="печатное колесо компрессора", status="оценка"),
    "печать_колесо_турбины": dict(trl=4, basis="печатное колесо турбины из никелевого сплава", status="оценка"),
    "печать_ротор_монолит": dict(trl=4, basis="ротор (колёса и вал) одной печатной деталью; прецедентов в ГТУ не найдено", status="оценка"),
    "печать_горячая_часть_целиком": dict(trl=5, basis="сопловой аппарат, камера и корпус турбины одной деталью", status="оценка"),
    "печать_оболочка": dict(trl=4, basis="цельная печатная оболочка двигателя (вход, корпуса компрессора и диффузоры, корпуса подшипников, "
                                         "каналы масла и топлива в стенке) с ротором, собранным внутри: печатные статорные детали - TRL 7, "
                                         "цельная горячая часть - 5; оболочка 0,6-1,5 м крупнее и требует ГИП в камерах, которых в мире "
                                         "две (Ø1,6 и Ø2,05 м), и КТ ускорителем 9 МэВ; прецедентов в ГТУ в треке не найдено - как "
                                         "печать_ротор_монолит (organic/print_review: shell_assembly, hip_ct)", status="оценка"),
    "печатное_охлаждение": dict(trl=5, basis="охлаждающие каналы, напечатанные в сопловом аппарате или колесе", status="оценка"),
    "печать_рекуператор": dict(trl=4, basis="печатные рекуператоры IN625 испытаны на стендах sCO2 (Thar Energy, 2022), для ГТУ "
                                            "эксплуатации не найдено; sco2symposium.com/proceedings2022/175-paper.pdf",
                               status="по выдержке поиска"),
    "литьё_ротор_монолит": dict(trl=6, basis="цельнолитые колёса малых двигателей IN713LC (casting-spec); ротор заодно "
                                             "0,45-0,60 м - «на пределе» (LIM)", status="оценка"),
    "сварка_трением_отливок": dict(trl=6, basis="колесо турбины с валом сваркой трением - турбокомпрессоры; две крупные "
                                                "отливки - не проверено", status="оценка"),
    "керамика_ротор": dict(trl=5, basis="CGT302 300 кВт: 37-42 % (1997-1999), рецессия Si3N4", status="по выдержке поиска"),
    "керамика_статика": dict(trl=6, basis="статорные детали из CMC", status="оценка"),
    "магнитные_подшипники": dict(trl=6, basis="SKF: ГТУ 9 МВт на магнитных подшипниках (описание изделия); Siemens: паровая турбина "
                                              "до 10 МВт на электромагнитных подшипниках, опытная эксплуатация; turbomachinerymag.com/view/"
                                              "10mw-steam-turbine-with-magnetic-bearings-delivered-to-utility", status="по выдержке поиска"),
    "кластер_ядер": dict(trl=5, basis="несколько газогенераторов на общую силовую турбину (семейство C gtu_explore.py); "
                                      "прецеденты не искались", status="оценка"),
    "одно_колесо_до_7": dict(trl=8, basis="одно центробежное колесо до π 7: прецедент OPRA OP16 (gtu_simple.py)", status="принято в треке"),
    "одно_колесо_выше_7": dict(trl=5, basis="одно трансзвуковое колесо выше π 7", status="оценка"),
    "rdc": dict(trl=3, basis="вращающаяся детонация: «still experimental»; ракетные RDE - TRL 4 из 9; NETL 2025 - проекты "
                             "интеграции в ГТУ (netl.doe.gov, FY25 Advanced Turbines Peer Review)", status="по выдержке поиска"),
    "волновой_ротор": dict(trl=3, basis="стендовые испытания четырёхпортового волнового ротора NASA Glenn; испытание на двигателе - "
                                        "«необходимый шаг» (ntrs.nasa.gov/citations/19960038355)", status="по выдержке поиска"),
}

# --------------------------------------------------------------------------- сроки и стоимость (оценка)
WEEKS = dict(setup_print=1.0, hip=1.0, machining_rot=1.0, cast_tooling=10.0, cast_pour_ndt=4.0, cast_machining=2.0,
             forging=12.0, forge_machining=3.0, fabrication=6.0, ceramic=20.0, hip_unique=6.0, ct_linac=2.0)
COST_U = dict(print_h=1.0, hip=50.0, cast_tooling=300.0, cast_kg=0.5, forge_kg=1.5, forge_fix=100.0, fab_kg=0.3, fab_fix=50.0,
              ceramic_fix=500.0, ceramic_kg=5.0, hip_unique=200.0, ct_linac=100.0)

# --------------------------------------------------------------------------- печатная оболочка: ГИП и контроль (v3, решение владельца
# 15.11.2026, п. 3; разбор живой геометрии - routes/gtu/organic/print_review/print_review_results.json -> hip_ct)
HIP_VESSELS = [0.88, 1.60, 2.05]   # м, диаметр рабочей камеры ГИП: серийная QIH 200 URC; MTC Powder Solutions (Швеция); Metal Technology
#                                    (Япония) - по выдержке поиска (organic/print_review.py: RULES["hip_largest"]); длина камер 2,15-4,2 м
#                                    больше любой оболочки воронки
HIP_MARGIN = 0.05                  # м, зазор детали в камере ГИП (как в print_review.hip_ct)
HIP_SERIAL_D = 0.88                # камеры не больше этой - серийные; крупнее - две в мире (очередь и перевозка: WEEKS["hip_unique"], оценка)
SHELL_RULES = dict(
    hip="оболочка - корпус под давлением и опоры подшипников из печатного металла: ГИП обязателен (усталость LPBF Ti-6Al-4V: после ГИП "
        "+10-20 %, по выдержке поиска - print_review RULES fatigue_ti64_lpbf); камера ГИП - не меньше D + 0,05 м",
    hip_unique="D + 0,05 м больше 0,88 м: только камеры Ø1,6 и Ø2,05 м (две в мире, по выдержке поиска) - +6 недель очереди и перевозки "
               "(оценка), качество способа не выше 0,6 (оценка); больше 2,05 м - ГИП невозможен, качество 0,2 (оценка)",
    ct="КТ стенки по касательной - сотни миллиметров металла (KV-3: 0,46 м, print_review hip_ct): трубкой 450 кВ не просвечивается, "
       "нужен ускоритель 9 МэВ с грубым вокселем; каналы - проливкой, свидетели-купоны в той же печати (print_review, оценка): +2 недели, "
       "качество способа не выше 0,6 (оценка)")


# =========================================================================== пределы, общие для воронки
def u_limit_comp(mat, method):
    """предел окружной скорости колеса компрессора, м/с, и основание"""
    if mat == "Ti64":
        return (620.0, "LIM U_ti: фрезерованное титановое колесо") if method == "поковка" else \
               (558.0, "оценка: 0,9 x 620 (литьё + ГИП или печать + ГИП; прочность ниже деформированного - casting-spec)")
    if mat == "Al":
        return (520.0, "LIM U_al") if method == "поковка" else (468.0, "оценка: 0,9 x 520")
    if mat == "Ni":
        return (480.0, "LIM U_mono: колесо компрессора заодно с турбиной из жаропрочного сплава (для печати принято так же - оценка)")
    raise ValueError(mat)


def dT_cap(mat_tw):
    v = MAT.get(mat_tw, {}).get("dT_cap")
    return 0.0 if v is None else v


def u_limit_turb(turb, mat_tw, T41_C, pt=False):
    """предел окружной скорости турбины (радиальная - периферия, осевая - средний диаметр), м/с"""
    if pt:
        return 380.0, "силовая турбина: не выше 380 м/с (gtu_explore.evaluate, gtu_family)"
    if turb == "О":
        return 420.0, "LIM U_ax (средний диаметр)"
    if mat_tw == "Si3N4":
        return 500.0, "оценка: как литое колесо при 900 C, без снижения с температурой (не проверено)"
    U = 500.0 - 0.7 * max(0.0, T41_C - dT_cap(mat_tw) - 900.0)
    return U, "LIM U_rad900 - 0,7 м/с на градус выше 900 C, со сдвигом dT_cap сплава"


NGV_T41 = {"U": (1000.0, "неохлаждаемый сопловой аппарат: класс U gtu_family (T41 <= 1000 C)"),
           "N": (1200.0, "охлаждаемый сопловой аппарат: класс K gtu_family (1200 C)"),
           "K": (1200.0, "класс K gtu_family"),
           "P": (1250.0, "оценка: печатные каналы охлаждения (не проверено)"),
           "C": (1350.0, "оценка: цель CGT302 - 1350 C (по выдержке поиска)")}
COOL_AIR = {"U": (0.0, 0.0), "N": (0.03, 0.0), "K": (0.055, 0.055), "P": (0.025, 0.03), "C": (0.0, 0.0)}
COOL_AIR_BASIS = "U, N, K - gtu_family.COOL; P - оценка (печатные каналы: 2,5 % на сопловой аппарат, 3 % на колесо, если оно " \
                 "напечатано); C - без охлаждения"


def rotor_cooled(cool, m_tw):
    return cool == "K" or (cool == "P" and m_tw == "печать")


def t41_limit(turb, cool, mat_tw, m_tw):
    """наибольшая температура газа перед рабочим колесом T41, C: меньшее из пределов соплового аппарата и колеса"""
    ngv, nb = NGV_T41[cool]
    if mat_tw == "Si3N4":
        rot, rb = 1350.0, "оценка: керамическое колесо, цель CGT302 1350 C (по выдержке поиска)"
    elif rotor_cooled(cool, m_tw):
        if cool == "K":
            rot, rb = ((1200.0, "класс K gtu_family: охлаждаемая монокристаллическая лопатка") if mat_tw == "SC" else
                       (1100.0 if turb == "О" else 1050.0, "оценка: охлаждаемое равноосное литое колесо"))
        else:
            rot, rb = (1150.0 if turb == "О" else 1100.0), "оценка: печатное колесо с каналами охлаждения"
            rot += min(0.0, dT_cap(mat_tw))
    else:
        if turb == "О" and mat_tw == "SC":
            rot, rb = 1050.0, "LIM T_ax_max: неохлаждаемая монокристаллическая лопатка"
        else:
            rot, rb = 950.0 + dT_cap(mat_tw), "LIM T_rad_max 950 C для литого K418 со сдвигом dT_cap сплава"
    return (ngv, nb) if ngv <= rot else (rot, rb)


# =========================================================================== посадка детали в способ изготовления
def _fits_box(D, L, box_mm, k=1.0):
    x, y = sorted(box_mm[:2]); z = box_mm[2]
    x, y, z = x * k / 1000.0, y * k / 1000.0, z * k / 1000.0
    return (D <= x and L <= z) or (L <= y and D <= x and D <= z)


def _fits_round(D, L, d_mm, z_mm):
    d, z = d_mm / 1000.0, z_mm / 1000.0
    return (D <= d and L <= z) or (math.hypot(D, L) <= d and D <= z)


def printers_for(D, L, rotating):
    """принтеры, в камеру которых деталь (тело вращения D x L, м) входит целиком; по возрастанию объёма камеры"""
    out = []
    for p in PRINTERS:
        if p["proc"] == "BJ" and rotating and not BJ_ROTATING:
            continue
        if "round_d" in p:
            ok = _fits_round(D, L, p["round_d"], p["z"]); vol = math.pi / 4 * p["round_d"] ** 2 * p["z"]
        else:
            ok = _fits_box(D, L, p["box"], BJ_SHRINK if p["proc"] == "BJ" else 1.0); vol = p["box"][0] * p["box"][1] * p["box"][2]
        if ok:
            out.append((vol, p))
    return [p for _, p in sorted(out, key=lambda t: t[0])]


def largest_print_dim():
    return max(max(p.get("box", (0, 0, 0))[:2] + (p.get("round_d", 0),)) for p in PRINTERS if p["proc"] == "LPBF") / 1000.0


def assess_part(part):
    """part: dict(name, kind, method, material, D, L, mass_kg, rotating, n).
    Возврат: fit (входит без деления), pieces (на сколько частей делится, если не входит), route, limit, basis,
    proc_q (качество способа 0..1), weeks, cost_u."""
    m, kind, mat = part["method"], part["kind"], part["material"]
    D, L, rot = part["D"], part["L"], part.get("rotating", False)
    mass = max(part.get("mass_kg", 0.0), 0.0)
    rho = MAT.get(mat, MAT["steel"])["rho"]
    vol_cm3 = mass / rho * 1e6
    res = dict(name=part["name"], method=m, material=mat, D=round(D, 3), L=round(L, 3), fit=True, pieces=1, route=None,
               limit=None, basis=None, proc_q=1.0, weeks=0.0, cost_u=0.0, hip=False)
    if m == "печать":
        ps = printers_for(D, L, rot)
        q = MAT.get(mat, {}).get("print")
        res["proc_q"] = 0.6 if q is None else q
        if ps:
            p = ps[0]; res.update(route=p["name"], basis=f"камера {p.get('box') or ('Ø%d x %d' % (p['round_d'], p['z']))} мм ({p['status']})")
            rate = p["rate"]
        elif not rot and max(D, L) <= DED["max_m"]:
            res.update(route="DED", basis=DED["note"]); rate = DED["rate"]; res["proc_q"] = min(res["proc_q"], 0.6)
        else:
            big = largest_print_dim()
            res.update(fit=False, pieces=max(2, math.ceil(max(D, L) / big)), route="не входит ни в одну камеру",
                       limit=big, basis="наибольшая камера LPBF (таблица PRINTERS)")
            rate = 900.0
        hours = vol_cm3 / rate + 24.0
        res["hip"] = rot or mat in ("IN718", "Haynes282", "CM247LC", "K418", "Ni")
        res["weeks"] = WEEKS["setup_print"] + hours / 168.0 + (WEEKS["hip"] if res["hip"] else 0.0) + (WEEKS["machining_rot"] if rot else 0.0)
        res["cost_u"] = hours * COST_U["print_h"] + (COST_U["hip"] if res["hip"] else 0.0)
        if part.get("shell"):                     # печатная оболочка: ГИП обязателен, контроль - ускоритель и купоны (SHELL_RULES)
            sh = shell_process(D)
            if sh["hip"] and not res["hip"]:
                res["weeks"] += WEEKS["hip"]; res["cost_u"] += COST_U["hip"]
            res.update(hip=sh["hip"], hip_vessel=sh["vessel"], control=sh["control"])
            res["proc_q"] = min(res["proc_q"], sh["proc_q"])
            res["weeks"] += sh["weeks"]
            res["cost_u"] += sh["cost_u"]
    elif m == "литьё":
        lim = CAST.get(kind) or CAST["static"]
        res.update(limit=lim["max"], basis=lim["basis"], route="литьё по выплавляемым моделям")
        q = MAT.get(mat, {}).get("cast")
        res["proc_q"] = 0.8 if q is None else q
        if D > lim["max"]:
            res.update(fit=False, pieces=max(2, math.ceil(D / lim["max"])))
        elif D > lim["ok"]:
            res["proc_q"] = min(res["proc_q"], 0.6)
        res["hip"] = mat in ("Ti64",) or (rot and mat in ("K418", "CM247LC", "Ni"))
        res["weeks"] = WEEKS["cast_tooling"] + WEEKS["cast_pour_ndt"] + (WEEKS["hip"] if res["hip"] else 0.0) + WEEKS["cast_machining"]
        res["cost_u"] = COST_U["cast_tooling"] + mass * COST_U["cast_kg"] + (COST_U["hip"] if res["hip"] else 0.0)
    elif m in ("поковка", "обработка"):
        lim = FORGE.get(kind) or FORGE["static"]
        res.update(limit=lim["max"], basis=lim["basis"], route="поковка или прокат с обработкой")
        if D > lim["max"]:
            res.update(fit=False, pieces=max(2, math.ceil(D / lim["max"])))
        if m == "поковка":
            res["weeks"] = WEEKS["forging"] + WEEKS["forge_machining"]
            res["cost_u"] = COST_U["forge_fix"] + mass * COST_U["forge_kg"]
        else:
            res["weeks"] = WEEKS["fabrication"]
            res["cost_u"] = COST_U["fab_fix"] + mass * COST_U["fab_kg"]
    elif m == "керамика":
        lim = CERAMIC.get(kind) or CERAMIC["static_hot"]
        res.update(limit=lim["max"], basis=lim["basis"], route="керамика: спекание и обработка")
        res["proc_q"] = MAT.get(mat, {}).get("ceramic", 0.5)
        if D > lim["max"]:
            res.update(fit=False, pieces=max(2, math.ceil(D / lim["max"])))
        res["weeks"] = WEEKS["ceramic"]
        res["cost_u"] = COST_U["ceramic_fix"] + mass * COST_U["ceramic_kg"]
    else:
        raise ValueError(m)
    res["weeks"] = round(res["weeks"], 2); res["cost_u"] = round(res["cost_u"], 1)
    return res


def shell_process(D):
    """печатная оболочка диаметром D, м: ГИП (камера по HIP_VESSELS) и контроль (КТ ускорителем, проливка, купоны) - SHELL_RULES.
    Возврат: hip (возможен), vessel (диаметр камеры или None), proc_q (верхняя граница качества способа), weeks и cost_u сверх печати"""
    need = D + HIP_MARGIN
    fit = [v for v in HIP_VESSELS if v >= need]
    weeks, cost, q = WEEKS["ct_linac"], COST_U["ct_linac"], 0.6
    if not fit:
        return dict(hip=False, vessel=None, proc_q=0.2, weeks=weeks, cost_u=cost,
                    control="ГИП невозможен (нет камеры) + КТ ускорителем 9 МэВ, проливка каналов, купоны")
    v = min(fit)
    if v > HIP_SERIAL_D:
        weeks += WEEKS["hip_unique"]; cost += COST_U["hip_unique"]
    return dict(hip=True, vessel=v, proc_q=q, weeks=weeks, cost_u=cost,
                control=f"ГИП в камере Ø{v:g} м" + (" (две в мире)" if v > HIP_SERIAL_D else "") +
                        " + КТ ускорителем 9 МэВ, проливка каналов, купоны")


# =========================================================================== балл
W = dict(fit=0.30, proc=0.25, trl=0.30, time=0.15)


def score(mk):
    """mk: dict(parts=[part, ...], techs=[ключи TRL]). Возврат: dict(score, s_fit, s_proc, s_trl, s_time, trl_min,
    trl_by, weeks, cost_u, parts=[оценка по деталям])"""
    rows = [dict(assess_part(p), n=p.get("n", 1)) for p in mk["parts"]]
    n_all = sum(r["n"] for r in rows) or 1
    s_fit = sum(r["n"] for r in rows if r["fit"]) / n_all
    s_proc = min((r["proc_q"] for r in rows), default=1.0)
    techs = sorted(set(mk.get("techs") or ["традиция"]))
    trl_by = {t: TRL[t]["trl"] for t in techs}
    trl_min = min(trl_by.values())
    s_trl = (trl_min - 1) / 8.0
    weeks = max((r["weeks"] for r in rows), default=0.0)
    s_time = 26.0 / (26.0 + weeks)
    S = 10.0 * (W["fit"] * s_fit + W["proc"] * s_proc + W["trl"] * s_trl + W["time"] * s_time)
    return dict(score=round(S, 3), s_fit=round(s_fit, 3), s_proc=round(s_proc, 3), s_trl=round(s_trl, 3), s_time=round(s_time, 3),
                trl_min=trl_min, trl_by=trl_by, weeks=round(weeks, 1), cost_u=round(sum(r["cost_u"] * r["n"] for r in rows), 1),
                parts=rows)


def facts():
    """таблицы фактов для results-JSON"""
    return dict(printers=PRINTERS, bj_shrink=BJ_SHRINK, ded=DED, cast=CAST, forge=FORGE, ceramic=CERAMIC,
                materials=MAT, trl=TRL, ngv_t41=NGV_T41, cool_air=COOL_AIR, cool_air_basis=COOL_AIR_BASIS,
                weeks=WEEKS, cost_u=COST_U, weights=W,
                shell=dict(hip_vessels_m=HIP_VESSELS, hip_margin_m=HIP_MARGIN, hip_serial_d_m=HIP_SERIAL_D, rules=SHELL_RULES))
