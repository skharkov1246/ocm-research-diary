"""Сборка targets.json (цели камеры KV): значения обоснования берутся из targets_basis_results.json по ключу.
Запуск после targets_basis.py:  python3 make_targets.py ; проверка — python3 -m pytest -q test_targets_basis.py"""
import json
import math
import os
import re

ROOT = os.path.dirname(os.path.abspath(__file__)) + "/"
R = json.load(open(ROOT + "targets_basis_results.json"))


def get(path, d=R):
    """'a.b[f=v].c' — выбор элемента списка по полю."""
    cur = d
    for part in re.findall(r"[^.\[\]]+(?:\[[^\]]+\])?", path):
        m = re.match(r"([^\[]+)(?:\[([^=]+)=([^\]]+)\])?$", part)
        name, f, v = m.groups()
        cur = cur[name]
        if f:
            vv = float(v)
            cur = next(x for x in cur if abs(float(x[f]) - vv) < 1e-9)
    return cur


def calc(key, what, unit=""):
    return dict(key=key, value=get(key), unit=unit, what=what)


S_OPEN = "открыт"
S_SEARCH = "по выдаче поиска"
S_EST = "оценка"

SRC = dict(
    solar_9ppm=("Solar Turbines, страница SoLoNOx: «9ppm NOx, 15ppm CO, and 15 ppm UHC emissions warranty for natural gas fuel»",
                "https://www.solarturbines.com/en_US/solutions/carbon-reduction/soLoNOx-technology-redirect.html"),
    solar_load=("Caterpillar/Solar, оценка формальдегида SoLoNOx: газ 50–100 % (Titan 250 — 40–100 %), жидкое 65–100 %",
                "https://mienviro.michigan.gov/ncore/downloadfile/619518159227452474"),
    solar_fleet=("Solar: более 3500 установок SoLoNOx, более 350 млн ч", "https://www.solarturbines.com/en_US/about-us/sustainability.html"),
    solar_h2=("Turbomachinery: SoLoNOx до 20 % H2, цель 100 % к 2030 г.", "https://www.turbomachinerymag.com/view/vendor-spotlight-solar-turbines"),
    titan130_h2=("Solar, GRTgaz Deutschland, OGE: Titan 130 SoLoNOx — смеси до 25 % H2, более 200 ч",
                 "https://www.powermag.com/successful-test-of-hydrogen-in-gas-turbine-driving-a-natural-gas-compressor-in-germany/"),
    mars_co=("Испытание Mars T14000 SoLoNOx: CO 2209 ppmc на 35 % нагрузки", "https://collections.lib.utah.edu/details?id=8748"),
    sgt400_enh=("Siemens Energy, Enhanced DLE Upgrade: SGT-400 13 МВт — NOx ≤ 9 ppmvd; экономия топлива на частичной нагрузке до 30 % при наличии перепускных клапанов",
                "https://assets.siemens-energy.com/dam/203b52b2-3f57-4cbc-ab57-b0580088391b/Enhanced-DLE-Upgrade-v3-pdf_Original%20file.pdf"),
    sgt400_30=("Siemens (DirectIndustry): SGT-400 выполняет нормы выбросов вплоть до 30 % нагрузки, двухтопливная",
               "https://www.directindustry.com/prod/siemens-power-genereration/product-23116-2448234.html"),
    sgt600_liq=("Siemens SGT-600 DLE 3-го поколения: жидкое < 58 ppm на полной нагрузке без воды, CO < 40 ppm до 75 %",
                "https://www.sciencegate.app/document/10.1115/gt2017-63089"),
    ge_water=("GE Vernova, Water Injection for NOx Reduction: 25 ppm газ, 42 ppm жидкое",
              "https://www.gevernova.com/gas-power/services/gas-turbines/upgrades/water-injection-for-nox-reduction"),
    novalt=("Baker Hughes NovaLT16 Overview: < 15 ppm NOx при 15 % O2 на 50–100 %, 9 ppm по заказу; до 30 % H2 в DLN",
            "https://dam.bakerhughes.com/m/543e00e8456feb1f/original/NovaLT16-Overview.pdf"),
    kawa_m7a=("Kawasaki, 2011: гарантия NOx M7A-03 — 9,9 ppm (15 % O2)",
              "https://global.kawasaki.com/en/corp/newsroom/news/detail/20110927_1.html"),
    kawa_h2=("KHI / FH Aachen / B&B-AGEMA: M1A-17, 100 % H2, микросмешение, 72 ppmv при 1800 кВт (15 % O2), 1040 ч",
             "https://opus.bibliothek.fh-aachen.de/opus4/solrsearch/index/search/searchtype/authorsearch/author/Tomoyuki+Ogino"),
    xonon=("Xonon на Kawasaki M1A-13X: NOx < 2,5 ppm, CO < 6 ppm за 3000 ч", "https://eponline.com/articles/2001/12/01/knocking-out-nox.aspx"),
    capstone=("Capstone C65 (перепечатка паспорта): жидкие топлива < 35 ppmvd NOx при 15 % O2",
              "https://pdf.directindustry.com/pdf/capstone-turbine/c65-c65-ichp-liquid-fuels-product-literature/22062-58097.html"),
    opra_liq=("Klein, Austrem, Mowill, GT2002-30107: OP16, дизель — NOx < 25 ppm на 50–90 %, < 30 ppm на полной; CO, UHC < 10 ppm",
              "https://next.sciencegate.app/document/10.1115/gt2002-30107"),
    opra_h2=("TKI / OPRA FlameSheet: 100 % H2 — NOx < 25 ppm на номинальной температуре; цель < 9 ppm не достигнута",
             "https://projecten.topsectorenergie.nl/storage/app/uploads/public/61f/7a3/45a/61f7a345ac28b478860182.pdf"),
    leap=("GE: наконечник форсунки LEAP — около 20 деталей → 1, −25 % массы, ресурс ×5",
          "https://www.ge.com/news/reports/epiphany-disruption-ge-additive-chief-explains-3d-printing-will-upend-manufacturing"),
    sgt700=("Siemens Energy: печатная горелка — 13 деталей / 18 швов → 1, 26 → 3 недели, более 8000 ч на ТЭС Philippsthal",
            "https://siemens-energy.com/us/en/home/stories/better-turbines-thanks-to-3d-printing.html"),
    cmc=("Solar Centaur 50S, керамические жаровые трубы: SiC/SiC с покрытием 13 937 ч; оксидная — более 25 000 ч",
         "https://www.osti.gov/biblio/991117"),
    twi=("TWI FAQ: типовые интервалы — осмотр камеры 8000 EOH, горячий тракт 24 000, капитальный 48 000–50 000 (не модельные)",
         "https://www.twi-global.com/technical-knowledge/faqs/faq-what-are-the-factors-influencing-maintenance-intervals-for-gas-turbines"),
    sgt600_int=("Siemens (презентация 2013): интервалы осмотров SGT-600 продлены с 20 000 до 30 000 EOH",
                "https://dtsheet.com/doc/1509571/download-presentation--pdf-"),
    nasa_hx=("NASA: допустимая температура металла голой жаровой трубы HX 1260 K при толщине около 2 мм",
             "https://ntrs.nasa.gov/api/citations/19830025622/downloads/19830025622.pdf"),
    calpoly=("Cal Poly (проект при поддержке Solar, 2021): печатный наконечник форсунки Taurus 60 — цель ресурса ≥ 60 000 ч",
             "https://core.ac.uk/download/483721164.pdf"),
    lefebvre_dp=("Лекции по Lefebvre: потери давления камеры — фундаментальные около 0,2 %, смешение 3–4 %, всего 4–5 %",
                 "https://conference.isabe.org/content/scnotes/bobby/gtc9.pdf"),
    tq=("Конференционная статья по проектированию камер: показатель качества поля температуры обычно 0,05–0,30",
        "https://conference.isabe.org/content/scnotes/bobby/gtc9.pdf"),
    dyn=("Обзор динамики бедного горения (Huang, Yang, PECS 2009): амплитуды до порядка 10 % давления; > 5 % среднего давления считают неустойчивостью (другая работа)",
         "https://www.yang.gatech.edu/publications/Journal/PECS%20(Huang,%202009).pdf"),
    sgt100_fleet=("Паспорт SGT-100: «automatic changeover from primary to secondary fuel at any load»; H2 до 30 % DLE",
                  "routes/gtu/fleet/gtu_fleet_sgt100.json"),
    sgt400_fleet=("Паспорт SGT-400: NOx < 10 ppmvd (DLE, газ); 6 противоточных жаровых труб с запальником в каждой; H2 до 30 % DLE",
                  "routes/gtu/fleet/gtu_fleet_sgt400.json"),
    t60_fleet=("Паспорт Taurus 60: 12 форсунок SoLoNOx, факельный воспламенитель; перепускной клапан на корпусе камеры регулирует воздух горения",
               "routes/gtu/fleet/gtu_fleet_taurus60.json"),
    m50_fleet=("Паспорт Mercury 50: 8 форсунок, NOx 5 ppm (минимальная возможность двигателя)", "routes/gtu/fleet/gtu_fleet_mercury50.json"),
    fuel_p=("Паспорта Solar и Siemens: давление топливного газа на входе пакета не опубликовано", "routes/gtu/fleet/gtu_fleet_sgt400.json"),
    kkkk=("40 CFR 60 подраздел KKKK, таблица 1: новые генераторные ГТУ на газе ≤ 50 MMBtu/ч — 42 ppm при 15 % O2",
          "https://www.ecfr.gov/current/title-40/chapter-I/subchapter-C/part-60/subpart-KKKK"),
    gost=("ГОСТ 29328-92: NOx ≤ 50 мг/м3, CO — по заказчику (выдержка автореферата)", "https://pstu.ru/files/2/file/adm/dissertacii/nugumanov/Avtoreferat_Nugumanov.pdf"),
    titan250_ldi=("Solar, GT2021-58986: Titan 250 — двухтопливная DLE и LDI-система для жидкого топлива на базе газовой DLE",
                  "https://asme-turboexpo.secure-platform.com/a/solicitations/137/sessiongallery/6190/application/58986"),
)


UNSURE = {  # выдержка поиска не привязана к адресу однозначно: даём все адреса выдачи, из которых она могла прийти
    "solar_load": ["https://dec.alaska.gov/Applications/Air/airtoolsweb/Home/ViewAttachment/17157994/GLlbpQjeIWuJdONN_f7DAg2"],
    "sgt400_30": ["https://www.directindustry.com/prod/siemens-power-genereration/product-23116-2019620.html"],
    "xonon": ["https://webarchive.library.unt.edu/eot2008/20090119083903/http:/www.netl.doe.gov/publications/proceedings/99/99ats/2-14.pdf",
              "https://www.power-eng.com/environmental-emissions/xonon-equipped-kawasaki-cogeneration-system-installed-at-plains-exploration-production-co-site/"],
    "calpoly": ["https://www.osti.gov/servlets/purl/2483903"],
    "lefebvre_dp": ["https://www.scribd.com/document/353061582/CIT04-0316-pdf"],
    "tq": ["https://www.scribd.com/document/353061582/CIT04-0316-pdf"],
    "dyn": ["https://torroja.dmt.upm.es/congresos/asme_2011/data/pdfs/trk-4/GT2011-46115.pdf"],
    "sgt600_int": [],
    "nasa_hx": ["https://ntrs.nasa.gov/api/citations/19840024434/downloads/19840024434.pdf"],
    "cmc": ["https://impact.ornl.gov/en/publications/microstructural-and-mechanical-characterization-of-a-hybrid-oxide/"],
}


def ref(k, status):
    t, u = SRC[k]
    r = dict(fact=t, source=u, status=status)
    if k in UNSURE:
        r["source_alt"] = UNSURE[k]
        r["source_note"] = "адрес — из выдачи поиска; к какому именно документу относится выдержка, не проверено"
    return r


targets = [
    dict(
        id="nox_gas", kind="hard",
        metric="NOx, ppmvd при 15 % O2, природный газ",
        target=dict(max=9.0, load_pct=[50, 100]), stretch=dict(max=5.0, load_pct=[50, 100]),
        also=dict(max=15.0, load_pct=[30, 50], note="нижняя часть диапазона — до 15 ppm (оценка)"),
        best_in_class=[ref("solar_9ppm", S_SEARCH), ref("kawa_m7a", S_SEARCH), ref("sgt400_enh", S_SEARCH),
                       ref("novalt", S_SEARCH), ref("m50_fleet", S_OPEN), ref("xonon", S_SEARCH)],
        basis_calc=[
            calc("hydrogen[H2_vol=0].NOx_ppm15_10ms", "идеальное смешение, T_ad зоны пламени 1800 K, 10 мс горячей зоны", "ppm"),
            calc("crn_grid[phi=0.55].hot[tau_ms=10].NOx_ppm15", "phi 0,55 (T_ad 1852 K), 10 мс", "ppm"),
            calc("crn_grid[phi=0.55].hot[tau_ms=20].NOx_ppm15", "phi 0,55, 20 мс", "ppm"),
            calc("unmixedness[S=0.05].by_tau[tau_ms=10].NOx_ppm15", "неравномерность смешения S = 0,05 при phi 0,55, 10 мс", "ppm"),
            calc("unmixedness[S=0.1].by_tau[tau_ms=10].NOx_ppm15", "S = 0,10 при phi 0,55, 10 мс", "ppm"),
            calc("crn_grid[phi=0.6].hot[tau_ms=10].NOx_ppm15", "phi 0,60 (T_ad 1936 K), 10 мс", "ppm"),
        ],
        rationale=("9 ppm — лучшая серийная гарантия класса у трёх поставщиков; цель «не хуже лучших» при одном клапане. "
                   "Расчёт 0D: при идеальном смешении и T_ad 1800–1850 K NOx 3–6 ppm за 10 мс, то есть запас под "
                   "неравномерность и погрешность модели (множитель около 1,5) есть, если S ≤ 0,05 и T_ad ≤ 1850 K; "
                   "при S = 0,10 или T_ad 1936 K цель не выполняется. Вытянутая цель 5 ppm — уровень Mercury 50."),
        derived=[dict(id="unmixedness", metric="S = sigma(phi)/phi на выходе смесителя", target=dict(max=0.05),
                      status=S_EST, why="NOx при S 0,05 → 0,10 растёт с 7,0 до 12,4 ppm (10 мс), появляется 3,8 % погасших струек"),
                 dict(id="flame_zone_Tad", metric="адиабатическая температура зоны пламени на 100 %, K",
                      target=dict(min=1750, max=1850), status=S_EST,
                      why="ниже 1750 K — близко к срыву (1680 K при 0,5 мс), выше 1850 K — NOx за 9 ppm при малой неравномерности")],
        verify=dict(levels=["L0 сеть реакторов (targets_basis.py)", "L1 RANS reactingFoam + сеть реакторов по полю",
                            "L2 LES смесителя (S, проскок)", "L3 стенд при p3"],
                    method="приведение C15 = C·(20,9−15)/(20,9−O2) сухого газа, средние с весом расхода на выходе"),
        status="цель (предложение); числа обоснования — оценка 0D",
    ),
    dict(
        id="co_gas", kind="hard",
        metric="CO, ppmvd при 15 % O2, природный газ",
        target=dict(max=15.0, load_pct=[50, 100]), also=dict(max=50.0, load_pct=[30, 50]),
        best_in_class=[ref("solar_9ppm", S_SEARCH), ref("xonon", S_SEARCH), ref("mars_co", S_SEARCH)],
        basis_calc=[
            calc("crn_grid[phi=0.5].hot[tau_ms=3].CO_ppm15", "phi 0,50 (T_ad 1764 K), 3 мс", "ppm"),
            calc("crn_grid[phi=0.55].hot[tau_ms=3].CO_ppm15", "phi 0,55 (T_ad 1852 K), 3 мс — CO уже у равновесия", "ppm"),
            calc("crn_grid[phi=0.55].eq.CO_ppm15", "равновесный CO при phi 0,55", "ppm"),
            calc("unmixedness[S=0.1].by_tau[tau_ms=10].unlit_weight", "доля погасших бедных струек при S = 0,10", "доля"),
        ],
        rationale=("15 ppm — гарантия Solar на газе. По расчёту CO в горячей зоне упирается в равновесие уже через 3 мс "
                   "(4–11 ppm при 1764–1852 K); реальный риск CO — частичная нагрузка и бедные струйки (при S = 0,10 "
                   "около 4 % смеси не горит в своём объёме) и закалка стенкой и воздухом разбавления, которых в 0D нет."),
        verify=dict(levels=["L1 RANS с охлаждением стенки", "L3 стенд"], method="как nox_gas"),
        status="цель (предложение)",
    ),
    dict(
        id="load_range", kind="hard",
        metric="диапазон нагрузки (доля мощности на валу), где выполнены nox_gas и co_gas; устойчивое горение",
        target=dict(low_emission_gas_pct=[30, 100], low_emission_diesel_pct=[50, 100], stable_pct=[0, 100]),
        best_in_class=[ref("solar_load", S_SEARCH), ref("sgt400_30", S_SEARCH), ref("sgt400_enh", S_SEARCH), ref("novalt", S_SEARCH)],
        basis_calc=[
            calc("part_load_fixed_split[T41_C=950].load_frac", "доля нагрузки при T41 950 °C (оценка вне расчётного режима)"),
            calc("part_load_fixed_split[T41_C=950].T_ad_pz_K", "T_ad зоны пламени при неизменной доле воздуха, T41 950 °C", "K"),
            calc("part_load_fixed_split[T41_C=900].load_frac", "доля нагрузки при T41 900 °C"),
            calc("part_load_fixed_split[T41_C=900].T_ad_pz_K", "T_ad при неизменной доле воздуха, T41 900 °C — ниже срыва", "K"),
            calc("psr_lbo[tau_ms=0.5].T_ad_at_lbo_K", "T_ad срыва реактора идеального смешения при 0,5 мс", "K"),
            calc("part_load_fixed_split[T41_C=700].load_frac", "доля нагрузки при T41 700 °C"),
            calc("part_load_fixed_split[T41_C=700].premix_share_needed_for_1800K", "доля воздуха в зону пламени, нужная для 1800 K при T41 700 °C"),
            calc("flame_temperature.premix_share_at_1800K", "та же доля на 100 %"),
        ],
        rationale=("Лучшие серийные: 50–100 % (Solar, NovaLT16), 30 % у SGT-400 только с перепуском воздуха; на жидком "
                   "у Solar 65–100 %. Цель «круче всех» — 30–100 % на газе без перепуска и без ступенчатых клапанов. "
                   "Расчёт показывает, почему это трудно: при неизменном делении воздуха и одном клапане зона пламени "
                   "держит температуру выше срыва только до T41 ≈ 950 °C (около 88 % нагрузки); на 77 % она уже ниже "
                   "срыва. Чтобы держать 1800 K до 33 % нагрузки, доля воздуха через зону пламени должна пассивно "
                   "уменьшаться с 0,48 до 0,24 — вдвое. Это главное требование к схеме подачи (пассивное деление)."),
        verify=dict(levels=["L0 сеть реакторов по нагрузке", "L1 RANS на 30/50/100 %", "L3 стенд"],
                    method="нагрузка — доля мощности на валу по gtu_cycle при постоянной частоте (оценка вне расчётного режима)"),
        status="цель (предложение); расчёт неполной нагрузки — оценка",
    ),
    dict(
        id="fuel_valves", kind="hard",
        metric="число регулирующих клапанов на одно топливо; исполнительных механизмов в камере",
        target=dict(control_valves_per_fuel=1, actuators_in_combustor=0, bleed_valves=0),
        best_in_class=[ref("t60_fleet", S_OPEN), ref("sgt400_fleet", S_OPEN)],
        basis_calc=[calc("part_load_fixed_split[T41_C=700].premix_share_needed_for_1800K", "что должно делать пассивное деление на 33 % нагрузки")],
        rationale=("Требование владельца (Plug&Play в пакет со штатной топливной системой). У SoLoNOx — пилот и основной "
                   "контур плюс перепускной клапан воздуха; у SGT-400 — пилот и основная горелка в каждой из 6 труб. "
                   "Деление топлива или воздуха по нагрузке должно быть заложено геометрией (пассивно) — см. load_range."),
        verify=dict(levels=["L1 RANS гидравлики каналов на 30/50/100 %"], method="распределение расхода по каналам из расчёта"),
        status="требование владельца",
    ),
    dict(
        id="diesel_dry", kind="hard",
        metric="NOx и CO на дизеле без воды, ppmvd при 15 % O2",
        target=dict(NOx_max=42.0, CO_max=50.0, load_pct=[50, 100], water=0), stretch=dict(NOx_max=25.0),
        best_in_class=[ref("ge_water", S_SEARCH), ref("sgt600_liq", S_SEARCH), ref("opra_liq", S_SEARCH),
                       ref("capstone", S_SEARCH), ref("titan250_ldi", S_SEARCH), ref("kkkk", S_SEARCH)],
        basis_calc=[
            calc("diesel.rows[phi=0.5].ign_delay_ms", "задержка самовоспламенения н-додекана при T3, p3, phi 0,5 (механизм Reitz)", "мс"),
            calc("diesel.rows[phi=1].ign_delay_ms", "то же при phi 1,0 (богатые карманы)", "мс"),
            calc("diesel.rows[phi=1].T_ad_K", "стехиометрическая T_ad н-додекана при T3, p3", "K"),
            calc("diesel.T_ad_stoich_CH4_K", "стехиометрическая T_ad метана", "K"),
        ],
        rationale=("42 ppm — уровень водяного впрыска на жидком (GE) и норма KKKK для газа малых ГТУ; сухая камера SGT-600 "
                   "даёт < 58 ppm, опытная LPP-камера OP16 — 25–30 ppm на стенде. Цель — водяной уровень без воды, "
                   "вытянутая — уровень OP16. Расчёт: стехиометрическое пламя дизеля горячее метанового (2557 против "
                   "2460 K), поэтому диффузионное горение дизеля даёт много NOx; предварительное испарение и смешение "
                   "при 410 °C и 13,7 бар ограничено самовоспламенением за 9–16 мс (механизм Reitz, н-додекан; реальный "
                   "дизель и другие механизмы могут дать в разы меньше) — смеситель должен быть на порядок короче. "
                   "Предпочтение — бедный прямой впрыск (LDI), как Titan 250."),
        verify=dict(levels=["L1 sprayFoam / распыл + испарение", "L3 стенд"], method="NOx дизеля — только механизм с азотной химией (в Reitz её нет)"),
        status="цель (предложение)",
    ),
    dict(
        id="fuel_changeover", kind="hard",
        metric="смена газ ↔ дизель на ходу",
        target=dict(both_directions=True, load_pct=[30, 100], stop=False, water=False, power_dip_max_pct=5.0),
        best_in_class=[ref("sgt100_fleet", S_OPEN), ref("solar_load", S_SEARCH)],
        basis_calc=[],
        rationale=("SGT-100 — автоматическая смена топлива при любой нагрузке (паспорт); у Kawasaki жидкое — резерв. "
                   "Цель: смена в обе стороны во всём диапазоне малых выбросов газа, без воды и останова; провал "
                   "мощности ≤ 5 % — оценка, уточнить у заказчика."),
        verify=dict(levels=["L2 нестационарный расчёт перехода", "L3 стенд"], method="—"),
        status="цель (предложение); провал мощности — оценка",
    ),
    dict(
        id="hydrogen", kind="objective",
        metric="доля водорода в газе (объёмн.) при выполнении nox_gas и co_gas",
        target=dict(H2_vol_min=0.30, hardware_change=False),
        stretch=dict(H2_vol=1.0, NOx_max=25.0, hardware_change="сменный картридж с тем же стыком"),
        best_in_class=[ref("sgt100_fleet", S_OPEN), ref("sgt400_fleet", S_OPEN), ref("novalt", S_SEARCH),
                       ref("solar_h2", S_SEARCH), ref("titan130_h2", S_SEARCH), ref("kawa_h2", S_SEARCH), ref("opra_h2", S_SEARCH)],
        basis_calc=[
            calc("hydrogen[H2_vol=0].S_L_m_s", "S_L метана при T_ad 1800 K, T3, p3", "м/с"),
            calc("hydrogen[H2_vol=0.3].S_L_m_s", "S_L при 30 % H2", "м/с"),
            calc("hydrogen[H2_vol=1].S_L_m_s", "S_L водорода", "м/с"),
            calc("hydrogen[H2_vol=0.3].NOx_ppm15_10ms", "NOx при 30 % H2, T_ad 1800 K, 10 мс", "ppm"),
            calc("hydrogen[H2_vol=1].NOx_ppm15_10ms", "NOx при 100 % H2, T_ad 1800 K, 10 мс", "ppm"),
        ],
        rationale=("30 % — уровень серийных DLE (SGT-100/400, NovaLT); Solar — 20 % серийно, 25 % испытано. При той же "
                   "T_ad ламинарная скорость при 30 % H2 выше на 13 %, при 100 % — в 2,4 раза; турбулентная скорость "
                   "водорода растёт сильнее (не учтено) — главный риск проскок. Самовоспламенения смеси при T3 за 1 с "
                   "нет ни для одной доли H2. 100 % H2 — сменным картриджем микросмешения (у Kawasaki 72 ppm на "
                   "двигателе, у FlameSheet < 25 ppm на стенде)."),
        verify=dict(levels=["L0 S_L и задержка", "L2 LES проскока", "L3 стенд"], method="запас по проскоку — корреляции раздела 5 STATE_OF_ART"),
        status="цель (предложение)",
    ),
    dict(
        id="pressure_loss", kind="objective",
        metric="полные потери давления камеры (от выхода диффузора до входа в сопловой аппарат), доля p3",
        target=dict(max=0.040), stretch=dict(max=0.030),
        best_in_class=[ref("lefebvre_dp", S_SEARCH)],
        basis_calc=[
            calc("dp_sensitivity.dP_kW_per_pct", "мощность на валу KV-3 на 1 п. п. потерь", "кВт"),
            calc("dp_sensitivity.d_eff_pt_per_pct", "КПД на валу на 1 п. п. потерь", "п. п."),
            calc("dp_sensitivity.rows[dp_comb=0.03].dP_kW_vs_4pct", "выигрыш при 3 % против 4 %", "кВт"),
        ],
        rationale=("4 % — значение модели цикла трека (gtu_cycle dp_comb = 0,04) и середина обычных 4–5 %. Каждый "
                   "процент потерь стоит KV-3 около 49 кВт на валу и 0,27 п. п. КПД. Струйные схемы (FLOX, микросмешение) "
                   "требуют перепада на струях — потери надо считать, а не принимать."),
        verify=dict(levels=["L1 RANS холодный и горячий"], method="разность полных давлений, средних по расходу"),
        status="цель (предложение)",
    ),
    dict(
        id="exit_uniformity", kind="hard",
        metric="неравномерность температуры на выходе: PF = (Tmax − T41)/(T41 − T3); радиальный профиль",
        target=dict(PF_max=math.floor(get("exit_uniformity.pattern_factor_max") * 1000) / 1000, RTDF_max=0.05),
        best_in_class=[ref("tq", S_SEARCH)],
        basis_calc=[calc("exit_uniformity.pattern_factor_max", "предел PF, чтобы местный максимум не превысил 1050 °C при T41 1000 °C")],
        rationale=("KV-3 — неохлаждаемая осевая турбина, предел газа 1050 °C (gtu_explore LIM['T_ax_max']) при T41 1000 °C: "
                   "местный максимум не выше +50 K, отсюда PF ≤ 0,084 (0,0848 с округлением вниз) — в нижней части обычных 0,05–0,30. RTDF 0,05 — оценка."),
        verify=dict(levels=["L1 RANS с разбавлением", "L3 стенд (гребёнка термопар)"], method="поле T на плоскости выхода"),
        status="цель; PF — расчёт из принятых величин трека, RTDF — оценка",
    ),
    dict(
        id="wall_life", kind="hard",
        metric="ресурс картриджа до замены; температура металла стенки",
        target=dict(life_h_min=30000, T_metal_max_K=1173),
        best_in_class=[ref("twi", S_SEARCH), ref("sgt600_int", S_SEARCH), ref("cmc", S_SEARCH), ref("nasa_hx", S_SEARCH), ref("calpoly", S_SEARCH)],
        basis_calc=[],
        rationale=("Типовой осмотр камеры — 8000 EOH, интервалы осмотров SGT-600 — 30 000 EOH; керамические трубы Solar — "
                   "13 937 и более 25 000 ч. Цель 30 000 ч — замена картриджа не чаще осмотров горячего тракта (оценка). "
                   "Температура металла ≤ 1173 K (900 °C) — запас 87 K к допустимой 1260 K голой HX по NASA (оценка); "
                   "прочность и ползучесть печатного сплава после HIP — расчёт CalculiX (навык calculix-rotor) и испытания."),
        verify=dict(levels=["L2 сопряжённый теплообмен chtMultiRegionFoam", "CalculiX ползучесть и МЦУ", "L3 стенд"], method="—"),
        status="цель (предложение); числа — оценка",
    ),
    dict(
        id="parts_count", kind="objective",
        metric="число отдельно изготовленных деталей модуля камеры (картриджа); сварных швов после печати; посадок с центровкой",
        target=dict(parts_max=5, welds_after_print=0, alignment_fits=0, fuel_connections_per_fuel=1),
        best_in_class=[ref("t60_fleet", S_OPEN), ref("sgt400_fleet", S_OPEN), ref("leap", S_SEARCH), ref("sgt700", S_SEARCH)],
        basis_calc=[],
        baseline=dict(taurus60_solonox_units_min=17, sgt400_units_min=18, status=S_EST,
                      how="узлы по открытым данным без деталей внутри форсунок: 12 форсунок + жаровая труба + воспламенитель + "
                          "перепускной клапан + ≥ 2 коллектора; SGT-400 — 6 × (горелка + труба + запальник)"),
        rationale=("Состав цели: 1 печатный моноблок (жаровая труба, смесители, топливные каналы, охлаждение), 1 корпус-фланец, "
                   "1 блок топливных разъёмов (оба топлива), 1 воспламенитель, 1 датчик. Прецеденты печати: 20 → 1 (LEAP), "
                   "13 деталей и 18 швов → 1 (горелка SGT-700)."),
        verify=dict(levels=["конструкторская спецификация", "проверка печатаемости (свесы ≤ 55°, удаление порошка)"], method="подсчёт по спецификации"),
        status="цель владельца (минимум деталей), число — предложение",
    ),
    dict(
        id="cartridge_swap", kind="objective",
        metric="замена картриджа камеры на месте",
        target=dict(time_h_max=2.0, people=2, special_tools=0, alignment=False, keyed_connections=True,
                    self_diagnosis=["наличие пламени", "пульсации давления", "температура стенки"]),
        best_in_class=[],
        basis_calc=[],
        rationale=("Требование владельца «хоть уборщиком». Открытых данных о времени замены камеры у Solar и Siemens не "
                   "найдено (пробел); цель 2 ч — оценка. Соединения — только с ключом формы (собрать неправильно нельзя)."),
        verify=dict(levels=["макет и хронометраж"], method="—"),
        status="цель владельца; число — оценка",
    ),
    dict(
        id="combustion_dynamics", kind="hard",
        metric="амплитуда пульсаций давления в камере, доля p3 (размах)",
        target=dict(max=0.01),
        best_in_class=[ref("dyn", S_SEARCH)],
        basis_calc=[],
        rationale=("Бедное горение склонно к автоколебаниям: амплитуды до порядка 10 % давления, > 5 % считают "
                   "неустойчивостью. Цель 1 % — пятикратный запас к этому порогу (оценка; приёмочные пределы изготовителей "
                   "закрыты)."),
        verify=dict(levels=["L2 LES / акустическая сеть", "L3 стенд"], method="эталон PRECCINSTA"),
        status="цель (предложение); число — оценка",
    ),
    dict(
        id="fuel_supply_compatibility", kind="hard",
        metric="перепад давления топлива на форсунках при штатном давлении топливного газа пакета",
        target=dict(value=None, note="не задано: давление топливного газа пакетов Solar и Siemens не опубликовано"),
        best_in_class=[ref("fuel_p", S_OPEN)],
        basis_calc=[],
        rationale="Plug&Play требует уложиться в штатное давление топлива; число — после запроса заказчику (навык gtu-dropin).",
        verify=dict(levels=["запрос заказчику"], method="—"),
        status="пробел",
    ),
]

guides = [
    dict(id="hot_zone_residence", metric="время пребывания в горячей зоне до разбавления, мс", value=[10, 20], status=S_EST,
         basis_calc=[calc("hot_zone_volume.rows[tau_ms=10].V_hot_L", "объём горячей зоны KV-3 при 10 мс", "л"),
                     calc("hot_zone_volume.rows[tau_ms=20].V_hot_L", "то же при 20 мс", "л"),
                     calc("crn_grid[phi=0.55].hot[tau_ms=30].NOx_ppm15", "NOx при 30 мс (phi 0,55) — длинная зона вредна", "ppm")],
         why="CO доходит до равновесия за 3 мс, NOx растёт почти линейно со временем; 10–20 мс — запас на частичную нагрузку и закалку"),
    dict(id="design_point", metric="вход камеры KV-3", value=None, status="расчёт 0D (повтор frame3)",
         basis_calc=[calc("inlet.p3_kPa", "давление за компрессором", "кПа"), calc("inlet.T3_C", "температура за компрессором", "°C"),
                     calc("inlet.W_air_comb_kg_s", "воздух в камеру", "кг/с"), calc("inlet.W_fuel_kg_s", "метан", "кг/с"),
                     calc("inlet.phi_overall", "общий коэффициент phi"), calc("inlet.heat_input_MW", "тепловая мощность по LHV", "МВт")],
         why="граничные условия всех уровней воронки камеры"),
]

funnel = dict(
    note=("Поля, которые воронка камеры отдаёт воронке двигателя (routes/gtu/funnel/), и пороги отбора. "
          "hard — отсев, objective — ранжирование."),
    outputs_to_engine_funnel=["dp_comb (доля p3)", "V_hot_L и габарит модуля (длина, наружный диаметр, м)", "масса модуля, кг",
                              "parts (число деталей модуля)", "NOx_ppm15 и CO_ppm15 на 100 / 50 / 30 %", "PF и RTDF",
                              "T_metal_max_K", "fuel_valves", "H2_vol_max", "diesel_NOx_ppm15", "cartridge_swap_h (оценка)"],
    hard=[t["id"] for t in targets if t["kind"] == "hard" and t["id"] != "fuel_supply_compatibility"],
    objectives=[dict(id="parts_count", sense="min"), dict(id="pressure_loss", sense="min"), dict(id="nox_gas", sense="min"),
                dict(id="hydrogen", sense="max"), dict(id="cartridge_swap", sense="min")],
    cycle_coupling=dict(field="dp_comb", file="routes/gtu/gtu_cycle.py → cycle(dp_comb=…)",
                        sensitivity_key="dp_sensitivity.dP_kW_per_pct"),
    engine_funnel_fields=dict(
        file="routes/gtu/funnel/funnel_space.py (сверено 14.11.2026 по тексту файла)",
        dp_comb="pfac = 1 − 0,04 (кольцевая) или 1 − 0,05 (выносная, кластер) — принятая величина; заменить итогом камеры",
        D_comb_L_comb="D_comb = 1,1·Dc, L_comb = 0,6·D_comb (кольцевая) — оценка; заменить габаритом картриджа",
        m_comb="масса жаровой трубы по толщине 2 мм и плотности 8220 кг/м3 — оценка; заменить массой картриджа",
        n_inj="round(π·0,8·D_comb/0,12), не меньше 6 — оценка; заменить числом смесителей выбранной схемы",
    ),
)

doc = dict(
    schema="kv.combustor.targets/1",
    title="Измеримые цели камеры сгорания двигателя-заменителя KV (первым — KV-3)",
    date="2026-11-14",
    status=("предложение целей; утверждение — владелец. Внешние факты — с пометкой «открыт» / «по выдаче поиска», "
            "наши числа — «оценка» или ключ в targets_basis_results.json"),
    method="сопоставление с открытыми данными серийных камер (STATE_OF_ART.md) + расчёт 0D (targets_basis.py)",
    software=R["software"],
    inputs=[dict(file="routes/gtu/combustor/targets_basis_results.json", field="см. basis_calc[].key"),
            dict(file="routes/gtu/combustor/STATE_OF_ART.md", field="разделы 1–5"),
            dict(file="routes/gtu/fleet/gtu_fleet_*.json", field="combustor, nox, fuel")],
    applies_to=dict(frame="KV-3", design_point="Taurus 60 GS, мягкий допуск: π 13,5, T41 1000 °C, 25,55 кг/с",
                    file="routes/gtu/frame3/frame3_point_results.json",
                    other_frames="KV-1 … KV-7: те же нормы в ppm и долях; PF пересчитать по T41 и T3 типоразмера"),
    conventions=dict(emissions="ppmvd при 15 % O2, NOx = NO + NO2 (как NO2), C15 = C·(20,9−15)/(20,9−O2,сух)",
                     mg_per_ppm=get("units.mg_Nm3_per_ppm_NO2"),
                     load="доля мощности на валу при постоянной частоте",
                     key_syntax="путь в targets_basis_results.json: a.b[поле=значение].c"),
    targets=targets,
    guides=guides,
    funnel_interface=funnel,
    gaps=[
        "Гарантии NOx на жидком топливе для SoLoNOx, SGT-100, SGT-400 в открытых источниках не найдены",
        "Время замены камеры у Solar и Siemens не опубликовано — цель cartridge_swap без эталона",
        "Давление топливного газа на входе пакетов не опубликовано — fuel_supply_compatibility не задано",
        "Доля пилотного топлива SoLoNOx и SGT-400 не найдена",
        "Приёмочные пределы пульсаций давления изготовителей закрыты — combustion_dynamics по литературе",
        "Тексты норм (ГОСТ 29328-92, KKKK таблица 1 целиком, MCPD приложение II) из среды не открывались",
        "Все внешние числа этого файла, кроме паспортов fleet/, — по выдаче поиска: перед публикацией сверить с первоисточником",
    ],
)

json.dump(doc, open(ROOT + "targets.json", "w"), ensure_ascii=False, indent=1)
print("ok", len(targets), "targets")
