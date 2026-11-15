#!/usr/bin/env python3
"""Модель массы агрегата-заменителя «топливо -> вал» по блокам, калиброванная по реальным машинам (уровень облика).

Зачем. Критерий владельца 15.11.2026 «КПД на тонну»: для типоразмера с фиксированной мощностью на валу (стык пакета)
считаются КПД и масса агрегата-заменителя - двигатель + камера + теплообменники + корпус и рама модуля; генератор,
редуктор, вход и выхлоп пакета НЕ входят. Показатели: КПД, масса, удельная мощность кВт/кг, «КПД на тонну» = КПД % /
масса т, фронт Парето «КПД - масса». Масса не должна подыгрывать: блоки калибруются по реальным двигателям.

Что делает модуль.
  1. Блоки массы (функции block_*): колёса и роторы по размерам и материалу (центробежное колесо - коэффициент формы
     по объёму CAD-тел frame3; осевые ступени - диск + лопатки; радиальная турбина; вал), корпуса по диаметру, давлению,
     материалу и способу изготовления (литьё в песок, точное литьё КНР, поковка/сварка, печать LPBF, лист), камера
     по объёму (кольцевая, трубчатая, выносная) - жаровая труба, корпус, форсунки; теплообменники по поверхности и
     стенке (первичная поверхность, пластинчато-ребристый, печатные пластины, печатный TPMS, трубчатый, риформер);
     подшипники (на масле пакета, своя маслосистема, магнитные с электроникой, лепестковые), опоры, крепёж, рама модуля.
  2. Размеры по паспортным параметрам (size_from_passport): расход, степень сжатия, ступени, частоты, температура
     выхлопа -> диаметры, высоты лопаток, длины, давления и температуры по проточной части (простая 0D-термодинамика
     с постоянными свойствами - только для размеров).
  3. Калибровка: масса двигателя = Σ ротор (по размерам, без множителя) + K_практики × Σ неподвижная часть
     (корпуса, аппараты, камера, опоры, крепёж и всё, чего блоки не видят: обвязка, коробка приводов, трубы).
     K - один параметр на практику («промышленная» - Solar; «авиационная» - GE LM, SGT-A05), МНК в логарифмах;
     проверка - исключение по одной машине (LOO), сравнение с наивной степенной зависимостью от мощности.
  4. evaluate(geom) - масса концепции с полосой неопределённости, кВт/кг, КПД на тонну; pareto_front - фронт.

Запуск:  python3 mass_model.py            -> mass_calibration.json (калибровка, проверки, примеры KV-3)
         python3 -m pytest -q test_mass_model.py
Все коэффициенты, кроме взятых из расчётов трека (frame3, organic) и источников mass_db.json, - оценка (поле basis).
"""
import copy
import datetime
import hashlib
import json
import math
import os
import sys

import numpy as np
from scipy.optimize import minimize_scalar

HERE = os.path.dirname(os.path.abspath(__file__))
PROCESS = os.path.dirname(HERE)
GTU = os.path.dirname(PROCESS)
REPO = os.path.abspath(os.path.join(GTU, "..", ".."))
FUNNEL = os.path.join(GTU, "funnel")
if FUNNEL not in sys.path:
    sys.path.append(FUNNEL)
import make_model as MK  # noqa: E402  (только таблица плотностей MAT; файл не меняется)

DATE = "2026-11-15"
DB_FILE = os.path.join(HERE, "mass_db.json")
OUT_FILE = os.path.join(HERE, "mass_calibration.json")
F_FRAME3_POINT = os.path.join(GTU, "frame3", "frame3_point_results.json")
F_FRAME3_GEOM = os.path.join(GTU, "frame3", "frame3_geom_results.json")
F_CYCLES = os.path.join(GTU, "organic", "cycles", "cycles_results.json")
F_FLEET = os.path.join(GTU, "fleet")
G0 = 9.80665

# =========================================================================== материалы
# плотности: из funnel/make_model.MAT (одна точка правды трека) + недостающие (справочные значения, оценка)
RHO = {k: MK.MAT[k]["rho"] for k in ("Ti64", "Al", "K418", "IN718", "HastX", "IN625", "steel", "SC", "CM247LC", "Haynes282")}
RHO_BASIS = {k: f"funnel/make_model.MAT['{k}'].rho" for k in RHO}
for _k, _v, _b in (("NodIron", 7100.0, "чугун с шаровидным графитом (корпуса Solar - 'Nodular Iron', TPS); 7100 - справочное, оценка"),
                   ("SS410", 7740.0, "нерж. 410 (корпус камеры Solar - '410 Stainless Steel', TPS); 7740 - справочное, оценка"),
                   ("SS347", 7960.0, "нерж. 347 (фольга рекуператора первичной поверхности, Mercury 50 - по выдержке поиска); 7960 - справочное, оценка"),
                   ("IN800", 7940.0, "Incoloy 800 (трубы риформера); 7940 - справочное, оценка"),
                   ("SiFe", 7650.0, "электротехническая сталь (магнитные подшипники); 7650 - справочное, оценка")):
    RHO[_k] = _v
    RHO_BASIS[_k] = _b

# допускаемое окружное напряжение стенки корпуса, МПа, от температуры металла (кусочно-линейно) - ОЦЕНКА:
# примерно треть предела текучести или предел ползучести на 1e5 ч по памяти; уточняется прочностью (CalculiX) позже
ALLOW = {
    "NodIron": [(20, 100.0), (350, 100.0), (450, 40.0)],
    "steel": [(20, 140.0), (400, 140.0), (500, 80.0), (550, 40.0)],
    "SS410": [(20, 120.0), (450, 120.0), (550, 60.0)],
    "Ti64": [(20, 250.0), (300, 250.0), (400, 180.0), (450, 120.0)],
    "Al": [(20, 80.0), (150, 80.0), (200, 40.0)],
    "IN718": [(20, 300.0), (650, 300.0), (700, 150.0), (750, 60.0)],
    "HastX": [(20, 120.0), (700, 120.0), (850, 60.0), (900, 30.0)],
    "IN625": [(20, 140.0), (650, 140.0), (800, 60.0), (900, 25.0)],
    "K418": [(20, 200.0), (750, 200.0), (850, 120.0), (950, 60.0)],
}


def sigma_allow(mat, T_C):
    """допускаемое напряжение, Па (оценка; таблица ALLOW)"""
    pts = ALLOW.get(mat) or ALLOW["steel"]
    T = [p[0] for p in pts]
    s = [p[1] for p in pts]
    return 1e6 * float(np.interp(T_C, T, s))


# способ изготовления корпуса: наименьшая стенка t_min, стенка по жёсткости/литейности k_D·D, множитель фланцев и рёбер
# (оценка; литьё в песок - толстая стенка крупных отливок; точное литьё КНР; поковка/сварка - авиационные кольца;
#  печать LPBF - тонкая стенка с рёбрами, фланцы заодно; лист - кожухи и газоходы)
PROC = {
    "cast_sand": dict(t_min=0.010, k_D=0.015, k_flange=1.35, name="литьё в песок"),
    "cast_invest": dict(t_min=0.004, k_D=0.008, k_flange=1.25, name="точное литьё (КНР)"),
    "forged_fab": dict(t_min=0.003, k_D=0.006, k_flange=1.30, name="поковка, сварка (авиационная)"),
    "lpbf": dict(t_min=0.002, k_D=0.005, k_flange=1.10, name="печать LPBF"),
    "sheet": dict(t_min=0.0015, k_D=0.003, k_flange=1.15, name="лист"),
}
PROC_BASIS = "оценка: t = max(t_min, k_D·D, p·D/(2σ_доп)); масса оболочки = π·D·L·t·ρ·k_flange"

# =========================================================================== допущения (оценка, если не указано иное)
A = dict(
    cp_a=1005.0, g_a=1.4, R_a=287.05, cp_g=1150.0, g_g=1.333, R_g=287.0, T0=288.15, p0=101325.0,
    eta_pc={"industrial": 0.88, "aero": 0.90}, eta_pt=0.88, eta_cf=0.84,
    dp_comb=0.04, p_exh_k=1.03,
    G_sp={"industrial": 185.0, "aero": 200.0}, htr1=0.50, htr_hp=0.60,
    U_tip_c={"industrial": 380.0, "aero": 420.0},
    c0_c=0.012, kc_c=0.30, Lst_c=2.2,
    c0_t=0.020, kc_t=0.35, Lst_t=2.0, h_Dm_max=0.30,
    psi_t={"industrial": 1.6, "aero": 1.8}, phi_t=0.55,
    lam_cf=0.68, phi_cf=0.085, nu_rt=0.70,
    sig_c=1.3, karea_c=0.056, taper_c=0.85, datt_c=0.5,
    sig_t=1.5, karea_t=0.15, karea_t_lp=0.08, hD_karea=(0.10, 0.25), taper_t=0.80, datt_t=0.84,
    U_ref_disk=420.0, fU_min=0.2,
    cax_c=0.8, w_rim_cax=1.15, att_fill=0.52,
    a_disk=None, b_disk=None,
    k_vane=1.2, k_ngv=1.5,
    k_cf_single=None, k_cf_double=None, k_bore=0.95,
    k_rt=0.099, k_cf_dif=0.03, k_rt_noz=0.02,
    shaft_d=0.15, shaft_fill=0.5,
    U_ref={"annular": 10.0, "can": 12.0, "silo": 8.0}, LH_annular=3.0, LD_can=3.0, H_liner=0.7,
    t_liner={"sheet": 0.0015, "lpbf": 0.0012}, k_liner=1.6, k_tp=0.3,
    m_inj={"industrial": 4.0, "aero": 2.5},
    k_frame_bh=3.0, k_fastener=0.02, m_pipe_bearing=15.0,
    T_case_hot_C=450.0,
    amb_p=0.4e6, amb_margin=3.0, amb_axial=4.0, amb_ctrl_kg=30.0, amb_backup_kg=3.0,
    foil_kg=2.0, foil_rotor_max_kg=150.0,
    oil_own_base=120.0, oil_own_per_kW=0.015,
    cradle_k=0.08, cradle_per_m=60.0,
)
A_BASIS = {
    "thermo": "постоянные cp, γ; КПД политропный компрессора 0,88 / 0,90, турбины 0,88 - оценка (только для размеров проточной части)",
    "G_sp": "удельный расход на входе осевого компрессора 185 / 200 кг/(с·м²), втулочное отношение 0,50 (ВД 0,60) - оценка",
    "U_tip_c": "окружная скорость на периферии первой ступени, если частота не дана, 380 / 420 м/с - оценка",
    "blades": "лопатки: хорда c = 12 мм + 0,30 h (компрессор) / 20 мм + 0,35 h (турбина), длина ступени 2,2 / 2,0 c - оценка; "
              "венец V = π·Dm·σ·k_пл·c·h·сужение: турбина σ 1,5, k_пл 0,15, сужение 0,80 (frame3_fea: шаг по Цвейфелю 0,826 c_a, "
              "turb_k_area 0,15, сужение 0,7/0,6/0,4); компрессор σ 1,3, k_пл 0,056 (t/c 0,08 × 0,7), сужение 0,85 - оценка",
    "disks": "диск осевой ступени: m_диска = a·ρ·D_жив³ + b·m_обода (обод = лопатки + замки: кольцо π·D_обода·1,15 c_a·d_зам·0,52, "
             "d_зам 0,84 c турбина / 0,5 c компрессор; D_жив = Dm - h - 2 d_зам); a, b - МНК по трём дискам турбины KV-3, "
             "рассчитанным на прочность (frame3_fea_results.json → parts.turbN.profile, rim_mass_kg; U_m 420 м/с); для компрессора те же a, b - оценка",
    "disk_U": "диск и замки масштабируются по напряжениям от вращения: × max(0,2; (U_m/420 м/с)²) (σ ∝ ρU²; 420 м/с - средняя скорость "
              "турбины KV-3, на которой посчитаны диски frame3_fea) - оценка; лопатки турбины тоньше при длинном пере: k_пл 0,15 при h/Dm ≤ 0,10 "
              "-> 0,08 при h/Dm ≥ 0,25 (лопатки НД/силовой турбины) - оценка",
    "h_Dm_max": "высота лопатки последней ступени не больше 0,30 Dm: если по ψ выходит больше, диаметр растёт по расходу (тихоходные силовые турбины) - оценка",
    "k_cf": "центробежное колесо m = k·ρ·D2³: k по объёму CAD-тел frame3 (frame3_geom_results.json → stages.stageN.solid.volume_m3 / D2³): одностороннее - ступень 2, двустороннее - ступень 1; без отверстия; с отверстием ×0,95 (оценка)",
    "k_rt": "радиальная турбина m = 0,099·ρ·D³ (0,0825·1,2 - gtu_explore.finalize), аппарат 0,02·ρ·D³ (funnel_space.derive) - оценка",
    "turbine": "коэффициент нагрузки ψ = Δh/U² 1,6 / 1,8, расхода φ = c_a/U 0,55, U/C0 радиальной 0,70 - оценка",
    "centrifugal": "коэффициент работы λ = w/U2² 0,68 (frame3: 0,648 и 0,716 - среднее), оптимальный φ = 0,085 (funnel_space.PHI_OPT['Ц']); осевая длина 0,30 D2 (одност.) / 0,65 D2 (двуст.) - по frame3; корпус ступени 1,45 D2 - оценка",
    "combustor": "площадь отсчёта A = G/(ρ3·U_ref), U_ref 10 (кольцевая) / 12 (трубчатая) / 8 (выносная) м/с; жаровая труба высотой 0,7 H, длиной 3 H; стенка 1,5 мм (лист) / 1,2 мм (печать), ×1,6 на фронтовое устройство и кольца; форсунка 4 / 2,5 кг; переходные патрубки 0,3 массы труб - оценка",
    "frames": "опоры подшипников (кольцо, стойки, ступица): 3·ρ·D²·t·k_фланца на опору при местном диаметре корпуса, t - стенка способа "
              "изготовления; крепёж 2 % неподвижной части; трубы масла 15 кг на подшипник - оценка",
    "case_hot": "металл корпуса турбины 450 °C (охлаждаемый/теплоизолированный) - оценка",
    "amb": "магнитный подшипник: удельная несущая способность 0,4 МПа на проекцию (B²/2μ0 при 1,0-1,2 Тл = 0,40-0,57 МПа × геометрия), запас 3 по весу ротора, осевой - 4 веса ротора; статор 2,2 D; контроллер 30 кг (funnel_space.derive), страховочные подшипники 3 кг - оценка",
    "foil": "лепестковый газовый подшипник 2 кг; опыт - роторы до 150 кг (оценка, не проверено)",
    "inlet_frame_cf": "у центробежной машины входной корпус 0,9 D2, передняя опора 0,6 D2 (у глаза колеса) - оценка",
    "cradle": "рама модуля: 0,08 опираемой массы + 60 кг на метр длины, сталь - оценка (без калибровки: рамы двигателей в открытых данных не отделены от рам пакетов)",
}


def _k_cf_from_frame3():
    """коэффициенты формы центробежного колеса по объёмам CAD-тел frame3: k = V / D2³ (одностороннее - ступень 2,
    двустороннее - ступень 1); если файла нет - 0,069 / 0,118 (те же числа, записанные 15.11.2026)"""
    try:
        g = json.load(open(F_FRAME3_GEOM))
        p = json.load(open(F_FRAME3_POINT))
        out = {}
        for key, st in (("double", "stage1"), ("single", "stage2")):
            V = g["stages"][st]["solid"]["volume_m3"]
            D2 = p[st]["D2_m"]
            out[key] = dict(k=V / D2 ** 3, V=V, D2=D2, sides=g["stages"][st]["solid"]["sides"], r_bore=g["stages"][st]["solid"]["r_bore"])
        return out
    except (OSError, KeyError, ValueError):
        return {"double": dict(k=0.1176, V=None, D2=None), "single": dict(k=0.0690, V=None, D2=None)}


F_FRAME3_FEA = os.path.join(GTU, "frame3", "frame3_fea_results.json")


def _disk_from_frame3_fea():
    """a, b модели диска по трём дискам турбины KV-3 из frame3_fea (МНК без свободного члена):
    m_диска = a·ρ·D_жив³ + b·m_обода; m_диска = ρ·∫2πr·t(r)dr (профиль parts.turbN.profile), m_обода = rim_mass_kg,
    D_жив = 2·r_live_rim. Если файла нет - a 0,0719, b 1,13 (те же числа, посчитанные 15.11.2026)"""
    try:
        d = json.load(open(F_FRAME3_FEA))
        X, y, rows = [], [], []
        for k in ("turb1", "turb2", "turb3"):
            part = d["parts"][k]
            r = np.array(part["profile"]["r"])
            t = np.array(part["profile"]["t"])
            rho = 7916.0 if "713" in str(part["material"]) else 8190.0
            m_d = rho * float(np.trapezoid(2 * np.pi * r * t, r))
            m_rim = part["geometry"]["rim_mass_kg"]
            Dl = 2 * part["geometry"]["r_live_rim_m"]
            X.append([rho * Dl ** 3, m_rim])
            y.append(m_d)
            rows.append(dict(part=k, disk_kg=round(m_d, 2), rim_kg=round(m_rim, 2), D_live=round(Dl, 4), rho=rho))
        (a, b), *_ = np.linalg.lstsq(np.array(X), np.array(y), rcond=None)
        pred = np.array(X) @ np.array([a, b])
        for rw, pp in zip(rows, pred):
            rw["fit_kg"] = round(float(pp), 2)
        return dict(a=float(a), b=float(b), rows=rows, rho_basis="IN-713C 7916 (platform/fea/materials.json)")
    except (OSError, KeyError, ValueError):
        return dict(a=0.0719, b=1.13, rows=None)


DISK = _disk_from_frame3_fea()
A["a_disk"] = DISK["a"]
A["b_disk"] = DISK["b"]
K_CF = _k_cf_from_frame3()
A["k_cf_single"] = K_CF["single"]["k"]
A["k_cf_double"] = K_CF["double"]["k"]

PRACTICE = {
    "industrial": dict(
        name="промышленная (литые корпуса, тяжёлые диски; Solar)",
        cold=("NodIron", "cast_sand"), warm=("steel", "cast_sand"), comb_case=("SS410", "cast_sand"),
        hot_case=("steel", "cast_sand"), diffuser=("NodIron", "cast_sand"), inlet=("NodIron", "cast_sand"),
        comp_rotor="steel", turb_rotor="IN718", liner=("HastX", "sheet"), vanes_c="steel", vanes_t="K418",
        include_inlet=True, include_diffuser=True),
    "aero": dict(
        name="авиационная (кованые и сварные корпуса из Ti и Ni, лёгкие диски; GE LM, SGT-A05)",
        cold=("Ti64", "forged_fab"), warm=("IN718", "forged_fab"), comb_case=("IN718", "forged_fab"),
        hot_case=("IN718", "forged_fab"), diffuser=("IN718", "sheet"), inlet=("Ti64", "forged_fab"),
        comp_rotor="Ti64", turb_rotor="IN718", liner=("HastX", "sheet"), vanes_c="Ti64", vanes_t="K418",
        include_inlet=False, include_diffuser=False),
}

# способ изготовления неподвижной части концепции (замена материалов и способа во всех корпусах практики): материалы печати LPBF -
# Ti64, IN718, IN625 (печатаемые сплавы - funnel/make_model.MAT['print']); точное литьё КНР - литая сталь ZG и ZL114A (навык casting-spec)
PROCESS_SETS = {
    "cast_sand": None,
    "cast_invest": dict(cold=("Al", "cast_invest"), warm=("steel", "cast_invest"), comb_case=("steel", "cast_invest"),
                        hot_case=("steel", "cast_invest"), diffuser=("steel", "cast_invest"), inlet=("Al", "cast_invest")),
    "lpbf": dict(cold=("Ti64", "lpbf"), warm=("IN718", "lpbf"), comb_case=("IN718", "lpbf"), hot_case=("IN718", "lpbf"),
                 diffuser=("IN625", "lpbf"), inlet=("Ti64", "lpbf"), liner=("HastX", "lpbf")),
}


GROUPS = ("rotor", "static", "hx", "bearings", "frame")
GROUP_NAMES = {"rotor": "ротор (по размерам, без множителя)", "static": "неподвижная часть (× K практики)",
               "hx": "теплообменники (свои множители, полоса)", "bearings": "подшипники сверх масляных (магнитные, электроника)",
               "frame": "рама модуля (оценка)"}


# =========================================================================== блоки массы
def _row(name, group, mass, basis, **detail):
    return dict(name=name, group=group, mass_kg=float(mass), basis=basis, detail=detail)


def wall_t(D, p_bar, T_C, mat, process):
    """толщина стенки оболочки, м: max(t_min, k_D·D, p·D/(2σ)) (оценка; PROC, ALLOW)"""
    pr = PROC[process]
    dp = max(p_bar - 1.013, 0.0) * 1e5
    return max(pr["t_min"], pr["k_D"] * D, dp * D / (2.0 * sigma_allow(mat, T_C)))


def block_shell(name, D, L, p_bar=1.013, T_C=20.0, mat="steel", process="cast_sand", group="static"):
    t = wall_t(D, p_bar, T_C, mat, process)
    m = math.pi * D * L * t * RHO[mat] * PROC[process]["k_flange"]
    return _row(name, group, m, PROC_BASIS, D=round(D, 4), L=round(L, 4), t_mm=round(t * 1e3, 2), p_bar=round(p_bar, 2),
                T_C=round(T_C), mat=mat, process=process)


def block_cone(name, D1, D2, L, p_bar=1.013, T_C=20.0, mat="steel", process="cast_sand"):
    Dm = 0.5 * (D1 + D2)
    s = math.hypot(L, 0.5 * abs(D2 - D1))
    t = wall_t(max(D1, D2), p_bar, T_C, mat, process)
    m = math.pi * Dm * s * t * RHO[mat] * PROC[process]["k_flange"]
    return _row(name, "static", m, PROC_BASIS + " (конус)", D1=round(D1, 4), D2=round(D2, 4), L=round(L, 4),
                t_mm=round(t * 1e3, 2), mat=mat, process=process)


def block_wheel_cf(name, D2, mat="Ti64", sides=1, bore=False):
    k = A["k_cf_double"] if sides == 2 else A["k_cf_single"]
    m = k * RHO[mat] * D2 ** 3 * (A["k_bore"] if bore else 1.0)
    return _row(name, "rotor", m, A_BASIS["k_cf"], D2=round(D2, 4), k=round(k, 4), sides=sides, mat=mat, bore=bore)


def block_wheel_rt(name, D, mat="K418"):
    return _row(name, "rotor", A["k_rt"] * RHO[mat] * D ** 3, A_BASIS["k_rt"], D=round(D, 4), mat=mat)


def _stage_parts(Dm, h, c, rho, kind, U=None):
    """одна осевая ступень: (масса лопаток, масса замков, масса диска, D_жив) - модель A_BASIS['blades'], ['disks'], ['disk_U'];
    U - средняя окружная скорость, м/с (None - как у дисков frame3_fea)"""
    t = kind != "comp"
    if t:
        ka = float(np.interp(h / Dm, A["hD_karea"], (A["karea_t"], A["karea_t_lp"])))
        sig, tap, datt = A["sig_t"], A["taper_t"], A["datt_t"]
    else:
        sig, ka, tap, datt = A["sig_c"], A["karea_c"], A["taper_c"], A["datt_c"]
    fU = 1.0 if U is None else max(A["fU_min"], (U / A["U_ref_disk"]) ** 2)
    m_bl = rho * math.pi * Dm * sig * ka * c * h * tap
    D_rim = max(Dm - h, 0.06)
    d_att = datt * c
    m_att = rho * math.pi * D_rim * A["w_rim_cax"] * A["cax_c"] * c * d_att * A["att_fill"] * fU
    D_live = max(D_rim - 2 * d_att, 0.04)
    m_disk = (A["a_disk"] * rho * D_live ** 3 + A["b_disk"] * (m_bl + m_att)) * fU
    return m_bl, m_att, m_disk, D_live


def block_axial_rotor(name, stages, mat, kind="comp"):
    """осевые ступени: лопатки + замки + диск; stages - [dict(Dm, h, c, mat?)]"""
    mb = ma = md = 0.0
    for s in stages:
        b, a, d, _ = _stage_parts(s["Dm"], s["h"], s["c"], RHO[s.get("mat", mat)], kind, s.get("U"))
        mb, ma, md = mb + b, ma + a, md + d
    return _row(name, "rotor", mb + ma + md, A_BASIS["disks"] + "; " + A_BASIS["blades"], n=len(stages),
                blades_kg=round(mb, 1), attach_kg=round(ma, 1), disks_kg=round(md, 1), mat=mat, kind=kind)


def block_vanes(name, stages, mat, kind="comp"):
    """направляющие и сопловые аппараты: венец как у лопаток × полки (1,2 компрессор / 1,5 турбина)"""
    k = A["k_vane"] if kind == "comp" else A["k_ngv"]
    m = sum(_stage_parts(s["Dm"], s["h"], s["c"], RHO[s.get("mat", mat)], kind)[0] * k for s in stages)
    return _row(name, "static", m, A_BASIS["blades"] + f"; полки ×{k}", n=len(stages), mat=mat, kind=kind)


def block_shaft(name, D, L, mat="steel", fill=None):
    fill = A["shaft_fill"] if fill is None else fill
    return _row(name, "rotor", RHO[mat] * math.pi / 4 * D ** 2 * L * fill, "вал: сталь, заполнение 0,5 (полый/стяжка) - оценка",
                D=round(D, 4), L=round(L, 4), mat=mat)


def block_cf_static(name, D2, mat):
    return _row(name, "static", A["k_cf_dif"] * RHO[mat] * D2 ** 3, "лопаточный диффузор и обратный канал 0,03·ρ·D2³ - оценка",
                D2=round(D2, 4), mat=mat)


def block_rt_nozzle(name, D, mat="K418"):
    return _row(name, "static", A["k_rt_noz"] * RHO[mat] * D ** 3, A_BASIS["k_rt"], D=round(D, 4), mat=mat)


def block_frames(name, Ds, mat="steel", process="cast_sand"):
    """опоры подшипников: на каждую опору 3·ρ·D²·t·k_фланца (Ds - местные диаметры корпуса), t - стенка способа (A_BASIS['frames'])"""
    pr = PROC[process]
    m = sum(A["k_frame_bh"] * RHO[mat] * D ** 2 * max(pr["t_min"], pr["k_D"] * D) * pr["k_flange"] for D in Ds)
    return _row(name, "static", m, A_BASIS["frames"], n=len(Ds), D=[round(D, 3) for D in Ds], mat=mat, process=process)


def block_combustor(c, practice):
    """камера: c = dict(kind, m_air, p_bar, T_C, D_mean?, n_cans?, n_inj?, liner_process?, case=(mat, process)?, integral?)
    integral=True - корпус камеры заодно с корпусом двигателя (печатная оболочка): отдельного корпуса нет"""
    pr = PRACTICE[practice]
    kind = c.get("kind", "annular")
    rho3 = c["p_bar"] * 1e5 / (A["R_a"] * (c["T_C"] + 273.15))
    Aref = c["m_air"] / (rho3 * A["U_ref"][kind])
    lproc = c.get("liner_process", pr["liner"][1])
    t_l = A["t_liner"]["lpbf" if lproc == "lpbf" else "sheet"]
    mat_l = c.get("liner_mat", pr["liner"][0])
    case_mat, case_proc = c.get("case", pr["comb_case"])
    rows = []
    if kind == "annular":
        Dm = c["D_mean"]
        H = Aref / (math.pi * Dm)
        Hl = A["H_liner"] * H
        Ll = A["LH_annular"] * Hl
        m_l = 2 * math.pi * Dm * Ll * t_l * RHO[mat_l] * A["k_liner"]
        rows.append(_row("камера: жаровая труба", "static", m_l, A_BASIS["combustor"], kind=kind, A_ref=round(Aref, 4), H=round(H, 4),
                         D_mean=round(Dm, 4), L_liner=round(Ll, 4), t_mm=t_l * 1e3, mat=mat_l))
        L_case = 1.3 * Ll
        if not c.get("integral"):
            rows.append(block_shell("камера: наружный корпус", Dm + H + 0.06, L_case, c["p_bar"], c["T_C"], case_mat, case_proc))
            Din = Dm - H - 0.06
            if Din > 0.10:
                rows.append(block_shell("камера: внутренний корпус", Din, L_case, c["p_bar"], c["T_C"], case_mat, case_proc))
        geom = dict(L=L_case, D=Dm + H + 0.06)
    else:
        n = int(c.get("n_cans", 1)) if kind == "can" else 1
        Dc = math.sqrt(4 * Aref / n / math.pi)
        Lc = A["LD_can"] * Dc
        Dl = 0.75 * Dc
        m_l = n * math.pi * Dl * Lc * t_l * RHO[mat_l] * A["k_liner"] * (1 + A["k_tp"])
        rows.append(_row("камера: жаровые трубы и патрубки", "static", m_l, A_BASIS["combustor"], kind=kind, n=n, D_can=round(Dc, 4),
                         L=round(Lc, 4), t_mm=t_l * 1e3, mat=mat_l))
        if not c.get("integral"):
            sh = block_shell("камера: корпус", 1.1 * Dc, Lc, c["p_bar"], c["T_C"], case_mat, case_proc)
            sh["mass_kg"] *= n
            sh["detail"]["n"] = n
            rows.append(sh)
        geom = dict(L=Lc, D=1.1 * Dc)
    n_inj = int(c.get("n_inj", 12))
    rows.append(_row("камера: форсунки и коллектор", "static", n_inj * A["m_inj"][practice], A_BASIS["combustor"], n_inj=n_inj))
    return rows, geom


# --------------------------------------------------------------------------- теплообменники
HX_TYPES = {
    "psr": dict(name="первичная поверхность, гофрированная фольга (Solar Mercury, Capstone)", t=0.10e-3, mat="SS347", k_core=1.15,
                k_hdr=(1.5, 2.0, 2.5), basis="фольга 0,1 мм (Mercury 50: 'about 0.1 mm wall' - по выдержке поиска); k_core 1,15 (кромки, "
                "прокладки) и коллекторы/корпус k_hdr 1,5-2,5 - оценка"),
    "plate_fin": dict(name="пластинчато-ребристый паяный", t=0.25e-3, mat="SS347", k_core=1.25, k_hdr=(1.5, 2.0, 2.5),
                      basis="эквивалентная стенка 0,25 мм (пластина + рёбра) - оценка"),
    "plate_print": dict(name="печатные пластины LPBF с опорами", t=0.35e-3, mat="IN625", k_core=1.45, k_hdr=(1.2, 1.5, 1.8),
                        basis="стенка 0,35 мм (organic/thermal); опоры листа ×1,38-1,58 (organic/verify_phys: 4,48-5,13 т против 3,25 т "
                        "без опор) -> 1,45; печатные коллекторы заодно k_hdr 1,2-1,8 - оценка"),
    "tpms_print": dict(name="печатная решётка TPMS (гироид) LPBF", t=0.30e-3, mat="IN625", k_core=1.0, k_hdr=(1.2, 1.5, 1.8),
                       basis="стенка 0,25-0,35 мм (organic/thermal/tpms_results.json; 0,30 - середина); листовая TPMS самонесущая k_core 1,0; "
                       "печатные коллекторы заодно k_hdr 1,2-1,8 - оценка"),
    "tube": dict(name="трубчатый (кожухотрубный, промежуточный охладитель)", t=1.0e-3, mat="steel", k_core=1.3, k_hdr=(1.5, 2.0, 2.5),
                 basis="трубка 1,0 мм, решётки и перегородки ×1,3, корпус и коллекторы k_hdr 1,5-2,5 - оценка"),
    "reformer": dict(name="риформер трубчатый с катализатором", t=2.0e-3, mat="IN800", k_core=1.3, k_hdr=(1.5, 2.0, 2.5), cat_rho=900.0,
                     basis="трубы 2 мм Incoloy 800, катализатор 900 кг/м³ насыпной - оценка"),
}


def hx_area(Q_kW=None, eps=None, C_min_kW_K=None, U_W_m2K=120.0, Cr=1.0):
    """поверхность противоточного теплообменника (по одной стороне), м²: NTU по ε и Cr; U 120 Вт/(м²·К) - gtu_explore.finalize
    (REC_K funnel_space); при Cr = 1 NTU = ε/(1-ε)"""
    if eps is None or C_min_kW_K is None:
        raise ValueError("нужны eps и C_min_kW_K")
    if abs(Cr - 1.0) < 1e-9:
        ntu = eps / (1 - eps)
    else:
        ntu = math.log((1 - eps * Cr) / (1 - eps)) / (1 - Cr)
    return ntu * C_min_kW_K * 1e3 / U_W_m2K


def block_hx(h, level=1):
    """теплообменник: h = dict(name, type, A_m2, t_wall_mm?, mat?, V_cat_m3?); level 0/1/2 - нижняя/средняя/верхняя k_hdr"""
    T = HX_TYPES[h["type"]]
    t = h.get("t_wall_mm", T["t"] * 1e3) * 1e-3
    mat = h.get("mat", T["mat"])
    core = h["A_m2"] * t * RHO[mat] * T["k_core"]
    tot = core * T["k_hdr"][level]
    cat = h.get("V_cat_m3", 0.0) * T.get("cat_rho", 0.0)
    return _row(h.get("name", "теплообменник"), "hx", tot + cat, T["basis"], type=h["type"], A_m2=round(h["A_m2"], 1),
                t_mm=round(t * 1e3, 3), mat=mat, core_kg=round(core, 1), k_hdr=T["k_hdr"][level], catalyst_kg=round(cat, 1))


# --------------------------------------------------------------------------- подшипники и рама
def block_bearings(b, rotor_kg, P_kW):
    """b = dict(kind: oil_pkg | oil_own | magnetic | foil, n_radial, n_thrust). Возвращает строки (группы static и bearings)"""
    kind = b.get("kind", "oil_pkg")
    nr, nt = int(b.get("n_radial", 2)), int(b.get("n_thrust", 1))
    rows = []
    if kind in ("oil_pkg", "oil_own"):
        rows.append(_row("подшипники скольжения: трубы и уплотнения", "static", (nr + nt) * A["m_pipe_bearing"], A_BASIS["frames"], n=nr + nt))
        if kind == "oil_own":
            rows.append(_row("своя маслосистема (бак, насос, охладитель)", "bearings", A["oil_own_base"] + A["oil_own_per_kW"] * P_kW,
                             "120 кг + 0,015 кг/кВт - оценка (у пакета Solar маслосистема своя и в агрегат-заменитель не входит)"))
    elif kind == "magnetic":
        F = A["amb_margin"] * G0 * rotor_kg / max(nr, 1)
        D = math.sqrt(F / A["amb_p"])
        m_r = nr * RHO["SiFe"] * math.pi / 4 * ((2.2 * D) ** 2 - D ** 2) * D * 0.6
        m_sleeve = nr * RHO["SiFe"] * math.pi / 4 * (D ** 2 - (0.6 * D) ** 2) * D
        Fa = A["amb_axial"] * G0 * rotor_kg
        Aa = Fa / A["amb_p"]
        m_a = nt * RHO["SiFe"] * Aa * 0.04 * 3
        m_ctrl = A["amb_ctrl_kg"] + (nr + nt) * A["amb_backup_kg"]
        rows.append(_row("магнитные подшипники: статоры, втулки, осевой диск", "bearings", m_r + m_sleeve + m_a, A_BASIS["amb"],
                         D_radial=round(D, 4), F_radial_N=round(F), F_axial_N=round(Fa)))
        rows.append(_row("магнитные подшипники: контроллер и страховочные", "bearings", m_ctrl, A_BASIS["amb"]))
    elif kind == "foil":
        rows.append(_row("лепестковые газовые подшипники", "bearings", (nr + nt) * A["foil_kg"], A_BASIS["foil"],
                         beyond_experience=rotor_kg > A["foil_rotor_max_kg"]))
    else:
        raise ValueError(f"неизвестный тип подшипников {kind}")
    return rows


def block_cradle(L, m_supported):
    return _row("рама модуля", "frame", A["cradle_k"] * m_supported + A["cradle_per_m"] * L, A_BASIS["cradle"], L=round(L, 3))


# =========================================================================== термодинамика для размеров
def _comp_stage(T, p, pi, eta_p):
    T2 = T * pi ** ((A["g_a"] - 1) / (A["g_a"] * eta_p))
    return T2, p * pi, A["cp_a"] * (T2 - T)


def size_from_passport(spec):
    """размеры проточной части и блоки по паспорту. spec - запись mass_db.json (m_air_kg_s, PR, P_kW, T_exh_C, comp, turb_gg,
    turb_pt, rpm_gg, rpm_pt, combustor, practice) или концепция KV в той же форме. Возвращает geom для evaluate()."""
    pr_name = spec.get("practice", "industrial")
    P = dict(PRACTICE[pr_name])
    if spec.get("process"):
        P.update(PROCESS_SETS[spec["process"]] or {})
    eta_pc = A["eta_pc"][pr_name]
    P_shaft = spec["P_kW"] * 1e3
    if spec.get("P_is_electric"):
        P_shaft = spec["P_kW"] * 1e3 / spec.get("eta_gen", 0.985) + spec.get("loss_gear_kW", 0.0) * 1e3
    m = spec.get("m_air_kg_s")
    # группы ступеней компрессора с давлением по ступеням (равные π ступеней, если не даны)
    groups = copy.deepcopy(spec["comp"])
    n_all = sum(g.get("n", 1) for g in groups)
    pis = []
    for g in groups:
        if "pi" in g:
            pis.append([g["pi"]] * g.get("n", 1) if not isinstance(g["pi"], list) else g["pi"])
    if not pis:
        pi_s = spec["PR"] ** (1.0 / n_all)
        for g in groups:
            g["pis"] = [pi_s] * g.get("n", 1)
    else:
        for g, pl in zip(groups, pis):
            g["pis"] = pl
    # компрессор: состояния
    T, p = A["T0"], A["p0"]
    w_c = {"gg": 0.0, "lp": 0.0}
    states = []
    for g in groups:
        sp = g.get("spool", "gg")
        sp = "lp" if sp == "lp" else "gg"
        g["_states"] = []
        for pi in g["pis"]:
            eta = A["eta_cf"] if g["type"] == "C" else eta_pc
            T2, p2, w = _comp_stage(T, p, pi, eta)
            g["_states"].append((T, p, T2, p2, w))
            w_c[sp] += w
            T, p = T2, p2
        states.extend(g["_states"])
    T2c, p2c = T, p
    wc_tot = w_c["gg"] + w_c["lp"]
    p3 = p2c
    p4 = p3 * (1 - A["dp_comb"])
    p5 = A["p0"] * A["p_exh_k"]
    # расход и температура на входе в турбину
    if m is None:
        T4 = spec["T4_K"]
        er = p4 / p5
        T5 = T4 * er ** (-(A["g_g"] - 1) / A["g_g"] * A["eta_pt"])
        w_t = A["cp_g"] * (T4 - T5)
        m = P_shaft / (w_t - wc_tot)
        m_basis = "из мощности и принятой T4 (оценка)"
    else:
        m_basis = "паспорт"
        w_t = wc_tot + P_shaft / m
        T4 = spec["T_exh_C"] + 273.15 + w_t / A["cp_g"]
    # частоты
    rows, dims = [], dict(m_air=m, m_basis=m_basis, T2_C=T2c - 273.15, p3_bar=p3 / 1e5, T4_C=T4 - 273.15, wc_kJkg=wc_tot / 1e3,
                          wt_kJkg=w_t / 1e3, P_shaft_kW=P_shaft / 1e3)
    rpm = {"gg": spec.get("rpm_gg"), "lp": spec.get("rpm_pt")}
    comp_stages_geo = []
    L_comp = 0.0
    D_case_comp = []
    cf_D2 = []
    rho01 = A["p0"] / (A["R_a"] * A["T0"])
    A1 = None
    for gi, g in enumerate(groups):
        sp = "lp" if g.get("spool") == "lp" else "gg"
        T_in, p_in = g["_states"][0][0], g["_states"][0][1]
        if g["type"] == "A":
            rho_in = p_in / (A["R_a"] * T_in)
            if A1 is None:
                A1 = m / A["G_sp"][pr_name]
                Ag = A1
                htr = A["htr1"]
            else:
                Ag = A1 * rho01 / rho_in
                htr = A["htr_hp"]
            Dt = math.sqrt(4 * Ag / (math.pi * (1 - htr ** 2)))
            Dm = 0.5 * Dt * (1 + htr)
            if rpm[sp] is None:
                rpm[sp] = 60 * A["U_tip_c"][pr_name] / (math.pi * Dt)
            stg = []
            for (Ti, pi_, To, po, w) in g["_states"]:
                rho = pi_ / (A["R_a"] * Ti)
                h = Ag * rho_in / rho / (math.pi * Dm)
                c = A["c0_c"] + A["kc_c"] * h
                mat = P["comp_rotor"] if not (pr_name == "aero" and To - 273.15 > 420) else "IN718"
                stg.append(dict(Dm=Dm, h=h, c=c, T_C=To - 273.15, p_bar=po / 1e5, mat=mat, U=math.pi * Dm * rpm[sp] / 60))
                L_comp += A["Lst_c"] * c
            g["_geo"] = stg
            comp_stages_geo += stg
            rows.append(block_axial_rotor(f"компрессор осевой {gi + 1}: диски и лопатки", stg, P["comp_rotor"], "comp"))
            rows.append(block_vanes(f"компрессор осевой {gi + 1}: направляющие аппараты", stg,
                                    P["vanes_c"] if pr_name == "industrial" else P["vanes_c"], "comp"))
            # корпус: холодная и тёплая части (граница 300 °C)
            Lc_cold = sum(A["Lst_c"] * s["c"] for s in stg if s["T_C"] < 300)
            Lc_warm = sum(A["Lst_c"] * s["c"] for s in stg if s["T_C"] >= 300)
            Dcase = (Dm + stg[0]["h"]) * 1.10
            D_case_comp.append(Dcase)
            if Lc_cold > 0:
                rows.append(block_shell(f"корпус компрессора {gi + 1} (холодный)", Dcase, Lc_cold, max(s["p_bar"] for s in stg if s["T_C"] < 300),
                                        max(s["T_C"] for s in stg if s["T_C"] < 300), *P["cold"]))
            if Lc_warm > 0:
                rows.append(block_shell(f"корпус компрессора {gi + 1} (тёплый)", Dcase, Lc_warm, max(s["p_bar"] for s in stg),
                                        max(s["T_C"] for s in stg), *P["warm"]))
            dims[f"comp{gi + 1}"] = dict(type="A", Dt1=round(Dt, 4), Dm=round(Dm, 4), h1=round(stg[0]["h"], 4), h_last=round(stg[-1]["h"], 4),
                                         rpm=round(rpm[sp]), U_tip1=round(math.pi * Dt * rpm[sp] / 60, 1), n=len(stg))
        else:  # центробежная ступень (одна запись = одна ступень)
            sides = g.get("sides", 1)
            (Ti, pi_, To, po, w) = g["_states"][0]
            U2 = math.sqrt(w / A["lam_cf"])
            if rpm[sp] is None:
                rho1 = pi_ / (A["R_a"] * Ti)
                D2 = math.sqrt((m / sides) / (rho1 * U2 * A["phi_cf"]))
                rpm[sp] = 60 * U2 / (math.pi * D2)
            D2 = 60 * U2 / (math.pi * rpm[sp])
            mat = g.get("mat", "Ti64" if To - 273.15 < 420 else "IN718")
            rows.append(block_wheel_cf(f"колесо центробежное {gi + 1}", D2, mat, sides, g.get("bore", False)))
            rows.append(block_cf_static(f"диффузор и обратный канал {gi + 1}", D2, P["warm"][0] if To - 273.15 >= 300 else P["cold"][0]))
            Lst = (0.65 if sides == 2 else 0.30) * D2 + 0.25 * D2
            rows.append(block_shell(f"корпус центробежной ступени {gi + 1}", 1.45 * D2, Lst, po / 1e5, To - 273.15,
                                    *(P["warm"] if To - 273.15 >= 300 else P["cold"])))
            L_comp += Lst
            D_case_comp.append(1.45 * D2)
            cf_D2.append(D2)
            dims[f"comp{gi + 1}"] = dict(type="C", D2=round(D2, 4), U2=round(U2, 1), rpm=round(rpm[sp]), sides=sides,
                                         pi=round(po / pi_, 3), T_out_C=round(To - 273.15, 1))
    # турбины
    tg = spec.get("turb_gg") or {}
    tp = spec.get("turb_pt")
    w_pt_load = P_shaft / m
    if tp:
        w_gg, w_pt = w_c["gg"], w_c["lp"] + w_pt_load
    else:
        w_gg, w_pt = wc_tot + w_pt_load, 0.0
    turb_list = [("gg", tg, w_gg)] + ([("lp", tp, w_pt)] if tp else [])
    T_in, w_done = T4, 0.0
    L_turb = 0.0
    D_tip_t = []
    last_tip = None
    for key, t, w in turb_list:
        n_rpm = rpm[key] if rpm.get(key) else rpm["gg"]
        if key == "lp" and spec.get("rpm_pt"):
            n_rpm = spec["rpm_pt"]
        if t.get("type", "A") == "A":
            nst = int(t.get("n", 1))
            dh = w / nst
            p_in_grp = p4 * (p5 / p4) ** (w_done / w_t)
            # состояния по ступеням (не зависят от диаметра)
            sts, Tg, wd = [], T_in, w_done
            for i in range(nst):
                T_out = Tg - dh / A["cp_g"]
                wd += dh
                p_out = p4 * (p5 / p4) ** (wd / w_t)
                sts.append((T_out, p_out / (A["R_g"] * T_out)))
                Tg = T_out
            U = math.sqrt(dh / A["psi_t"][pr_name])
            Dm = 60 * U / (math.pi * n_rpm)
            # тихоходная силовая турбина: h последней ступени <= h_Dm_max·Dm -> диаметр по расходу (A_BASIS['h_Dm_max'])
            Q_last = m / sts[-1][1]
            Dm_flow = (Q_last / (math.pi * A["h_Dm_max"] * A["phi_t"] * math.pi * n_rpm / 60)) ** (1 / 3)
            dm_by = "psi"
            if Dm_flow > Dm:
                Dm, dm_by = Dm_flow, "flow"
                U = math.pi * Dm * n_rpm / 60
            stg = []
            for (T_out, rho) in sts:
                h = (m / rho) / (math.pi * Dm * A["phi_t"] * U)
                c = A["c0_t"] + A["kc_t"] * h
                stg.append(dict(Dm=Dm, h=h, c=c, T_C=T_out - 273.15, U=U))
                L_turb += A["Lst_t"] * c
            T_in, w_done = sts[-1][0], wd
            name = "турбина газогенератора" if key == "gg" else "силовая турбина"
            rows.append(block_axial_rotor(f"{name}: диски и лопатки", stg, P["turb_rotor"], "turb"))
            rows.append(block_vanes(f"{name}: сопловые аппараты", stg, P["vanes_t"], "turb"))
            Dcase = (Dm + stg[-1]["h"]) * 1.15
            Lt = sum(A["Lst_t"] * s["c"] for s in stg)
            rows.append(block_shell(f"{name}: корпус", Dcase, Lt, p_in_grp / 1e5, A["T_case_hot_C"], *P["hot_case"]))
            D_tip_t.append(Dm + stg[-1]["h"])
            last_tip = Dm + stg[-1]["h"]
            dims["turb_" + key] = dict(type="A", n=nst, Dm=round(Dm, 4), U=round(U, 1), rpm=round(n_rpm), h1=round(stg[0]["h"], 4),
                                       h_last=round(stg[-1]["h"], 4), Dm_by=dm_by, psi=round(dh / U ** 2, 3))
        else:  # радиальная турбина
            dhs = w / A["eta_pt"]
            U = A["nu_rt"] * math.sqrt(2 * dhs)
            D = 60 * U / (math.pi * n_rpm)
            T_in -= w / A["cp_g"]
            w_done += w
            rows.append(block_wheel_rt("колесо радиальной турбины", D, P["vanes_t"]))
            rows.append(block_rt_nozzle("сопловой аппарат радиальной турбины", D, P["vanes_t"]))
            rows.append(block_shell("корпус (улитка) радиальной турбины", 1.6 * D, 0.5 * D, p4 / 1e5, A["T_case_hot_C"], *P["hot_case"]))
            L_turb += 0.5 * D
            D_tip_t.append(D)
            last_tip = 0.6 * D
            dims["turb_" + key] = dict(type="R", D=round(D, 4), U=round(U, 1), rpm=round(n_rpm))
    # валы
    Dmax_gg = max([s["Dm"] + s["h"] for s in comp_stages_geo] + cf_D2 + D_tip_t[:1] + [0.1])
    rows.append(block_shaft("вал газогенератора", A["shaft_d"] * Dmax_gg, L_comp + L_turb * (0.5 if tp else 1.0) + 0.5 * Dmax_gg))
    if tp:
        rows.append(block_shaft("вал силовой турбины", A["shaft_d"] * D_tip_t[-1], L_turb * 0.5 + 0.5 * D_tip_t[-1]))
    # камера
    comb = dict(spec.get("combustor") or {"kind": "annular"})
    D_mean = comb.get("D_mean") or max(max(D_case_comp) / 1.10, D_tip_t[0] if D_tip_t else 0.3)
    comb.update(m_air=m, p_bar=p3 / 1e5, T_C=T2c - 273.15, D_mean=D_mean)
    crow, cgeo = block_combustor(comb, pr_name)
    rows += crow
    # вход и выхлопной диффузор
    first = dims.get("comp1", {})
    D_in = (0.9 * first["D2"]) if first.get("type") == "C" else max(D_case_comp[0] / 1.10, 0.3) * 1.3   # вход: у глаза колеса / у первой ступени
    L_in = 0.4 * D_in
    if P["include_inlet"]:
        rows.append(block_shell("входной корпус", D_in, L_in, 1.013, 20.0, *P["inlet"]))
    L_dif = 0.0
    if P["include_diffuser"]:
        L_dif = 1.0 * last_tip
        rows.append(block_cone("выхлопной диффузор", last_tip, 1.6 * last_tip, L_dif, 1.05, spec["T_exh_C"] if spec.get("T_exh_C") else 500.0,
                               *P["diffuser"]))
    # опоры
    n_spools = 2 if tp else 1
    D_case_max = max(D_case_comp + [1.15 * d for d in D_tip_t])
    D_front = (0.6 * first["D2"]) if first.get("type") == "C" else D_case_comp[0]   # передняя опора: у глаза колеса / у первой ступени
    Ds_fr = [D_front] + [1.15 * d for d in D_tip_t]                                # + опора за каждой турбиной
    rows.append(block_frames("опоры подшипников (корпуса, стойки)", Ds_fr, *P["warm"]))
    L_eng = (L_in if P["include_inlet"] else 0.0) + L_comp + cgeo["L"] + L_turb + L_dif
    dims.update(L_engine=round(L_eng, 3), D_case_max=round(D_case_max, 3), combustor=dict(kind=comb["kind"], L=round(cgeo["L"], 3),
                                                                                         D=round(cgeo["D"], 3)))
    geom = dict(name=spec.get("name", spec.get("id", "?")), practice=pr_name, rows=rows, dims=dims,
                bearings=spec.get("bearings", {"kind": "oil_pkg", "n_radial": 2 * n_spools, "n_thrust": n_spools}),
                hx=spec.get("hx", []), frame=spec.get("frame"), P_shaft_kW=P_shaft / 1e3, eta_pct=spec.get("eta_pct"))
    return geom


# =========================================================================== оценка массы концепции
_CAL = {}


def load_calibration(path=OUT_FILE):
    """K и σ (LOO) по практикам из mass_calibration.json"""
    if path not in _CAL:
        d = json.load(open(path))
        _CAL[path] = {k: dict(K=v["K"], sigma_log=v["loo"]["rms_log"]) for k, v in d["calibration"]["per_practice"].items()}
    return _CAL[path]


def raw_blocks(geom):
    """все строки блоков без множителя K: строки геометрии + подшипники + крепёж + теплообменники (средние) + рама"""
    rows = [dict(r) for r in geom["rows"]]
    m_rot = sum(r["mass_kg"] for r in rows if r["group"] == "rotor")
    P_kW = geom.get("P_shaft_kW") or 0.0
    rows += block_bearings(geom.get("bearings") or {"kind": "oil_pkg"}, m_rot, P_kW)
    m_st = sum(r["mass_kg"] for r in rows if r["group"] == "static")
    rows.append(_row("крепёж", "static", A["k_fastener"] * m_st, A_BASIS["frames"]))
    for h in geom.get("hx") or []:
        rows.append(block_hx(h, 1))
    return rows


def sums(rows):
    return {g: sum(r["mass_kg"] for r in rows if r["group"] == g) for g in GROUPS}


def evaluate(geom, K=None, sigma_log=None, cal=None):
    """масса агрегата: ротор + K·неподвижная + теплообменники + подшипники сверх масляных + рама; полоса: двигатель ×exp(±σ_LOO),
    теплообменники - нижний/верхний k_hdr; рама - на опираемую массу. Возвращает dict с блоками и итогами"""
    pr = geom.get("practice", "industrial")
    if K is None or sigma_log is None:
        c = (cal or load_calibration())[pr]
        K = c["K"] if K is None else K
        sigma_log = c["sigma_log"] if sigma_log is None else sigma_log
    rows = raw_blocks(geom)
    s = sums(rows)
    eng = s["rotor"] + K * s["static"]
    hx_lo = sum(block_hx(h, 0)["mass_kg"] for h in geom.get("hx") or [])
    hx_hi = sum(block_hx(h, 2)["mass_kg"] for h in geom.get("hx") or [])
    core = eng + s["hx"] + s["bearings"]
    fr = geom.get("frame")
    m_frame = 0.0
    if fr:
        L = fr.get("L") or geom.get("dims", {}).get("L_engine", 3.0)
        r = block_cradle(L, core)
        rows.append(r)
        m_frame = r["mass_kg"]
    tot = core + m_frame
    lo = eng * math.exp(-sigma_log) + hx_lo + s["bearings"]
    hi = eng * math.exp(sigma_log) + hx_hi + s["bearings"]
    if fr:
        lo += block_cradle(fr.get("L") or geom.get("dims", {}).get("L_engine", 3.0), lo)["mass_kg"]
        hi += block_cradle(fr.get("L") or geom.get("dims", {}).get("L_engine", 3.0), hi)["mass_kg"]
    out = dict(name=geom.get("name"), practice=pr, K=K, sigma_log=sigma_log, groups_raw={k: round(v, 1) for k, v in s.items()},
               engine_kg=round(eng, 1), hx_kg=round(s["hx"], 1), bearings_extra_kg=round(s["bearings"], 1), frame_kg=round(m_frame, 1),
               total_kg=round(tot, 1), total_lo_kg=round(lo, 1), total_hi_kg=round(hi, 1),
               rows=[dict(r, mass_kg=round(r["mass_kg"] * (K if r["group"] == "static" else 1.0), 2)) for r in rows])
    P = geom.get("P_shaft_kW")
    eta = geom.get("eta_pct")
    if P:
        out["kW_per_kg"] = round(P / tot, 4)
    if eta:
        out["eta_pct"] = eta
        out["kpd_per_t"] = round(kpd_per_tonne(eta, tot), 3)
        out["kpd_per_t_range"] = [round(kpd_per_tonne(eta, hi), 3), round(kpd_per_tonne(eta, lo), 3)]
    return out


def kpd_per_tonne(eta_pct, mass_kg):
    """показатель владельца: КПД, % / масса, т"""
    return eta_pct / (mass_kg / 1000.0)


def pareto_front(items, eta_key="eta_pct", mass_key="total_kg"):
    """недоминируемые по (КПД max, масса min); возвращает индексы в порядке роста массы"""
    idx = sorted(range(len(items)), key=lambda i: (items[i][mass_key], -items[i][eta_key]))
    front, best = [], -1e300
    for i in idx:
        if items[i][eta_key] > best + 1e-12:
            front.append(i)
            best = items[i][eta_key]
    return front


# =========================================================================== калибровка
def engine_model_parts(spec):
    """(m_rot, m_static_raw, geom) двигателя по паспорту, без теплообменников и рамы"""
    g = size_from_passport(spec)
    s = sums(raw_blocks(g))
    return s["rotor"], s["static"], g


def fit_K(parts, masses):
    """K практики: min Σ (ln(m_rot + K·m_st) - ln m)²; parts - [(m_rot, m_st)]"""
    def f(lk):
        K = math.exp(lk)
        return sum((math.log(r + K * s) - math.log(mm)) ** 2 for (r, s), mm in zip(parts, masses))
    res = minimize_scalar(f, bounds=(math.log(0.05), math.log(50.0)), method="bounded", options=dict(xatol=1e-10))
    return math.exp(res.x)


def _stats(errs_log):
    e = np.asarray(errs_log)
    return dict(rms_log=float(np.sqrt(np.mean(e ** 2))), mean_log=float(np.mean(e)), max_abs_pct=float(np.max(np.abs(np.expm1(e))) * 100),
                rms_pct=float(np.sqrt(np.mean(np.expm1(e) ** 2)) * 100))


def calibrate(db):
    engines = [e for e in db["engines"] if e["use"] == "calibration"]
    out = dict(per_practice={}, table=[])
    parts_all = {}
    for e in engines:
        r, s, g = engine_model_parts(e)
        parts_all[e["id"]] = (r, s, g)
    for pr in ("industrial", "aero"):
        es = [e for e in engines if e["practice"] == pr]
        parts = [parts_all[e["id"]][:2] for e in es]
        ms = [e["mass_kg"] for e in es]
        K = fit_K(parts, ms)
        ins = [math.log(r + K * s) - math.log(mm) for (r, s), mm in zip(parts, ms)]
        loo, rows = [], []
        for i, e in enumerate(es):
            Ki = fit_K(parts[:i] + parts[i + 1:], ms[:i] + ms[i + 1:])
            r, s = parts[i]
            pred = r + Ki * s
            loo.append(math.log(pred) - math.log(ms[i]))
            lo, hi = e.get("mass_range_kg", [ms[i], ms[i]])
            rows.append(dict(id=e["id"], name=e["name"], practice=pr, mass_kg=ms[i], mass_range_kg=[lo, hi], P_kW=e["P_kW"],
                             kW_per_kg=round(e["P_kW"] / ms[i], 3), rotor_kg=round(r, 1), static_raw_kg=round(s, 1),
                             rotor_share_pct=round(100 * r / (r + K * s), 1), pred_insample_kg=round(r + K * s, 1),
                             K_loo=round(Ki, 4), pred_loo_kg=round(pred, 1), err_loo_pct=round(100 * (pred / ms[i] - 1), 1),
                             within_source_range=bool(lo * 0.999 <= pred <= hi * 1.001),
                             dims=parts_all[e["id"]][2]["dims"]))
        out["per_practice"][pr] = dict(name=PRACTICE[pr]["name"], n=len(es), K=round(K, 4), insample=_stats(ins), loo=_stats(loo),
                                       machines=[e["id"] for e in es])
        out["table"] += rows
    # одна K на все машины (без деления на практику) - показать, что практика нужна
    pr_all = [parts_all[e["id"]][:2] for e in engines]
    ms_all = [e["mass_kg"] for e in engines]
    K1 = fit_K(pr_all, ms_all)
    loo1 = []
    for i in range(len(engines)):
        Ki = fit_K(pr_all[:i] + pr_all[i + 1:], ms_all[:i] + ms_all[i + 1:])
        loo1.append(math.log(pr_all[i][0] + Ki * pr_all[i][1]) - math.log(ms_all[i]))
    out["single_K_all"] = dict(K=round(K1, 4), loo=_stats(loo1), note="одна K без деления на практику - для сравнения")
    # наивные модели: m = a·P^b и m = P/q (по практике), LOO
    base = {}
    for pr in ("industrial", "aero"):
        es = [e for e in engines if e["practice"] == pr]
        P = np.array([e["P_kW"] for e in es], float)
        M = np.array([e["mass_kg"] for e in es], float)
        lp, lm = np.log(P), np.log(M)
        e_pow, e_q = [], []
        for i in range(len(es)):
            msk = np.arange(len(es)) != i
            b, a = np.polyfit(lp[msk], lm[msk], 1)
            e_pow.append(a + b * lp[i] - lm[i])
            q = np.mean(lm[msk] - lp[msk])
            e_q.append(q + lp[i] - lm[i])
        b, a = np.polyfit(lp, lm, 1)
        base[pr] = dict(power_law=dict(a=float(math.exp(a)), b=float(b), loo=_stats(e_pow)),
                        const_specific=dict(kW_per_kg=float(math.exp(-np.mean(lm - lp))), loo=_stats(e_q)))
    out["baselines"] = base
    # чувствительность: Taurus 70 = 5107 кг (фунты документа) вместо 5630 кг
    es = [e for e in engines if e["practice"] == "industrial"]
    parts = [parts_all[e["id"]][:2] for e in es]
    ms = [5107.0 if e["id"] == "taurus70_md" else e["mass_kg"] for e in es]
    K_alt = fit_K(parts, ms)
    loo_alt = []
    for i in range(len(es)):
        Ki = fit_K(parts[:i] + parts[i + 1:], ms[:i] + ms[i + 1:])
        loo_alt.append(math.log(parts[i][0] + Ki * parts[i][1]) - math.log(ms[i]))
    out["sensitivity_taurus70_5107"] = dict(K=round(K_alt, 4), loo=_stats(loo_alt),
                                            note="Taurus 70: 11 260 lb = 5107 кг вместо напечатанных 5630 кг")
    return out


# =========================================================================== проверки и примеры
def _fleet_packages():
    """удельная масса пакетов по паспортам fleet (кг/кВт): чтение, без правки"""
    rows = []
    for fn in sorted(os.listdir(F_FLEET)):
        if not fn.endswith(".json"):
            continue
        d = json.load(open(os.path.join(F_FLEET, fn)))
        for r in d.get("records", []):
            if r.get("pkg_kg") and r.get("power_kW"):
                rows.append(dict(model=r.get("model", d.get("model")), variant=r.get("variant"), P_kW=r["power_kW"], pkg_kg=r["pkg_kg"],
                                 kg_per_kW=round(r["pkg_kg"] / r["power_kW"], 2), engine_kg=r.get("engine_kg"),
                                 source=f"routes/gtu/fleet/{fn}"))
    return rows


def kv3_spec(practice="industrial", **over):
    """KV-3 простой цикл по frame3_point_results.json: 25,55 кг/с, π 4,49 × 3,007, вал 14 950 об/мин, осевая турбина 3 ступени"""
    p = json.load(open(F_FRAME3_POINT))
    cyc = p["cycle"]["frame3"]
    s1, s2 = p["stage1"], p["stage2"]
    spec = dict(id="kv3_simple", name="KV-3 простой (frame3)", practice=practice, P_kW=cyc["P_shaft_kW"], eta_pct=round(100 * cyc["eff_shaft"], 2),
                m_air_kg_s=s1["m_total_kg_s"], PR=round(s1["PR"] * s2["PR"], 3), T_exh_C=cyc["EGT_C"],
                comp=[dict(type="C", sides=2, pi=s1["PR"]), dict(type="C", sides=1, pi=s2["PR"])],
                turb_gg=dict(type="A", n=int(p["turbine"]["selected_z"])), turb_pt=None, rpm_gg=s1["n_rpm"], rpm_pt=None,
                combustor=dict(kind="annular", n_inj=12), bearings=dict(kind="oil_pkg", n_radial=2, n_thrust=1),
                frame=dict(L=None))
    spec.update(over)
    return spec


def checks(db, cal):
    out = {}
    calk = {k: dict(K=v["K"], sigma_log=v["loo"]["rms_log"]) for k, v in cal["per_practice"].items()}
    # 1. KG2-3 (радиальная промышленная, проверка)
    kg = next(e for e in db["engines"] if e["id"] == "kg2_3")
    spec = dict(kg, T4_K=1150.0, combustor=dict(kind="silo", n_inj=1))
    g = size_from_passport(spec)
    ev = evaluate(g, cal=calk)
    out["kg2_3_radial"] = dict(real_kg=kg["mass_kg"], pred_kg=ev["engine_kg"], band_kg=[ev["total_lo_kg"], ev["total_hi_kg"]],
                               err_pct=round(100 * (ev["engine_kg"] / kg["mass_kg"] - 1), 1), dims=g["dims"],
                               note="π 4,0 и T4 1150 K - оценка (в источнике нет); единственная найденная масса всё-радиальной промышленной машины")
    # 2. frame3: размеры по паспортному пути против геометрии frame3
    p = json.load(open(F_FRAME3_POINT))
    g3 = size_from_passport(kv3_spec())
    out["frame3_sizing"] = dict(D2_model=[g3["dims"]["comp1"]["D2"], g3["dims"]["comp2"]["D2"]], D2_frame3=[p["stage1"]["D2_m"], p["stage2"]["D2_m"]],
                                D2_err_pct=[round(100 * (g3["dims"]["comp1"]["D2"] / p["stage1"]["D2_m"] - 1), 1),
                                            round(100 * (g3["dims"]["comp2"]["D2"] / p["stage2"]["D2_m"] - 1), 1)],
                                turb_Dm_model=g3["dims"]["turb_gg"]["Dm"],
                                turb_Dm_frame3=p["turbine"]["options"][str(p["turbine"]["selected_z"])]["Dm"],
                                note="проверка размерного пути: λ 0,68 (среднее двух ступеней frame3) против своих λ ступеней")
    # 3. рекуператор Mercury 50: модель первичной поверхности против разности пакетов Mercury 50 - Centaur 50
    m_air = 17.69 / 1.015
    C_air = m_air * 1.07
    for eps in (0.85, 0.90):
        Aa = hx_area(eps=eps, C_min_kW_K=C_air, U_W_m2K=143.2, Cr=1 / 1.08)
        rows = [block_hx(dict(type="psr", A_m2=Aa), lv)["mass_kg"] for lv in (0, 1, 2)]
        out.setdefault("mercury50_recuperator", {})[f"eps_{eps}"] = dict(A_m2=round(Aa, 1), mass_lo_mid_hi_kg=[round(x) for x in rows])
    out["mercury50_recuperator"].update(package_delta_kg=45700 - 37785, U_basis="U 143,2 Вт/(м²·К) - фольга organic/cycles (hx foil)",
                                        eps_basis="ε рекуператора Mercury 50 не опубликован в открытых источниках - 0,85 и 0,90 оценка",
                                        verdict="модель должна давать меньше разности пакетов (в разность входят газоходы, удлинение пакета, "
                                                "разница двигателей)")
    # 4. SGT-A05: граница Siemens
    a05 = next(t for t in cal["table"] if t["id"] == "sgta05_kb5")
    out["sgta05_bound"] = dict(siemens_upper_kg=766, centrax_kg=800, pred_loo_kg=a05["pred_loo_kg"])
    # 5. удельные массы пакетов fleet (кг/кВт) и отношения «модуль / двигатель»
    pk = _fleet_packages()
    out["fleet_packages_kg_per_kW"] = dict(min=min(r["kg_per_kW"] for r in pk), max=max(r["kg_per_kW"] for r in pk), n=len(pk), rows=pk)
    out["skid_to_engine"] = [dict(id=x["id"], ratio=round(x["mass_kg"] / x["engine_kg"], 2)) for x in db["packages"] if x.get("engine_kg")]
    return out


def examples(cal):
    """KV-3: простой (frame3) в двух практиках и с разными подшипниками; KV-3 + рекуператор 0,90 + подогрев «как есть» (organic/cycles)
    с разными сердцевинами - масса, кВт/кг, КПД на тонну, фронт Парето"""
    calk = {k: dict(K=v["K"], sigma_log=v["loo"]["rms_log"]) for k, v in cal["per_practice"].items()}
    items = []
    variants = [("industrial", None, "oil_pkg", "промышл. практика, литьё в песок, масло пакета"),
                ("industrial", None, "magnetic", "промышл. практика, литьё в песок, магнитные"),
                ("industrial", "cast_invest", "oil_pkg", "промышл. K, точное литьё КНР, масло пакета"),
                ("industrial", "lpbf", "oil_pkg", "промышл. K, печать LPBF, масло пакета"),
                ("industrial", "lpbf", "magnetic", "промышл. K, печать LPBF, магнитные"),
                ("aero", None, "oil_pkg", "авиац. практика, масло пакета"),
                ("aero", None, "magnetic", "авиац. практика, магнитные")]
    for pr, proc, bk, label in variants:
        sp = kv3_spec(pr, bearings=dict(kind=bk, n_radial=2, n_thrust=1), process=proc)
        sp["name"] = f"KV-3 простой, {label}"
        ev = evaluate(size_from_passport(sp), cal=calk)
        ev["process"] = proc or ("cast_sand" if pr == "industrial" else "forged_fab")
        items.append(ev)
    # «как есть»: компрессор KV-3 π 13,5, рекуператор 0,90 + подогрев (cycles_results.json, key R90-ic-RH-ibc-none-none, pr135)
    try:
        cy = json.load(open(F_CYCLES))
        r = next(x for x in cy["results"] if x.get("key") == "R90-ic-RH-ibc-none-none")["pr135"]
        A_foil = r["hx"]["foil"][0]["A_m2"]
        A_lpbf = r["hx"]["lpbf"][0]["A_m2"]
        for typ, Aa in (("psr", A_foil), ("tpms_print", A_lpbf), ("plate_print", A_lpbf)):
            sp = kv3_spec("industrial")
            # тот же компрессор, другой расход и КПД (cycles); турбина - 3 ступени (размеры по своему перепаду)
            sp.update(name=f"KV-3 + рекуператор 0,90 + подогрев (как есть), сердцевина {typ}", m_air_kg_s=r["W_air"], eta_pct=round(100 * r["eff_el"], 2),
                      P_kW=kv3_spec()["P_kW"], T_exh_C=r["T_rec_in_C"], hx=[dict(name="рекуператор", type=typ, A_m2=Aa)],
                      eta_basis="КПД на клеммах (cycles eff_el), у простого - на валу (frame3): сравнение смещено на ~1,5 п. п. в пользу простого")
            ev = evaluate(size_from_passport(sp), cal=calk)
            ev["note"] = ("межтурбинная камера подогрева не добавлена отдельным блоком; турбина в размерах - перепад до входа в рекуператор; "
                          "мощность принята равной простому KV-3 (cycles считает на 1 кг/с - масштаб по расходу)")
            items.append(ev)
    except (OSError, StopIteration, KeyError):
        pass
    for ev in items:
        ev.pop("rows", None)
    front = pareto_front(items)
    return dict(items=items, pareto_idx=front, pareto_names=[items[i]["name"] for i in front])


def _sha(path):
    try:
        return hashlib.sha256(open(path, "rb").read()).hexdigest()[:16]
    except OSError:
        return None


def build():
    db = json.load(open(DB_FILE))
    t0 = datetime.datetime.now()
    cal = calibrate(db)
    chk = checks(db, cal)
    ex = examples(cal)
    eng_db = []
    for e in db["engines"] + db["reference_only"]:
        eng_db.append(dict(id=e["id"], name=e["name"], practice=e["practice"], use=e.get("use", "reference"), mass_kg=e["mass_kg"],
                           P_kW=e["P_kW"], kW_per_kg=round(e["P_kW"] / e["mass_kg"], 3), eta_pct=e.get("eta_pct"),
                           status="; ".join(sorted({v["status"] for v in (e.get("sources") or {"s": e.get("source", {})}).values() if v.get("status")}))))
    res = dict(
        title="Модель массы агрегата-заменителя по блокам, калиброванная по реальным двигателям (уровень облика)",
        method=("блоки массы по размерам и материалам (колёса - коэффициент формы CAD frame3, осевые ступени - диск + лопатки, корпуса - стенка "
                "по давлению, материалу и способу изготовления, камера по объёму, теплообменники по поверхности и стенке, подшипники, рама); "
                "размеры по паспорту - 0D-термодинамика с постоянными свойствами; калибровка: m = Σ ротор + K_практики·Σ неподвижная, "
                "МНК в логарифмах, проверка исключением по одной машине (LOO) и сравнение со степенной зависимостью от мощности"),
        software="Python 3, numpy, scipy (minimize_scalar)",
        date=DATE, status="оценка уровня облика: калибровка по 8 двигателям (3 промышленных, 5 авиационных); блоковое деление массы не проверено "
        "по реальным машинам (известны только итоги); теплообменники и рама не калиброваны (полоса)",
        inputs=dict(db=os.path.relpath(DB_FILE, REPO), db_sha16=_sha(DB_FILE), frame3_point=os.path.relpath(F_FRAME3_POINT, REPO),
                    frame3_geom=os.path.relpath(F_FRAME3_GEOM, REPO), cycles=os.path.relpath(F_CYCLES, REPO), fleet=os.path.relpath(F_FLEET, REPO),
                    densities="routes/gtu/funnel/make_model.py: MAT[*].rho", code_sha16=_sha(os.path.abspath(__file__))),
        assumptions=dict(A={k: v for k, v in A.items()}, basis=A_BASIS, rho=RHO, rho_basis=RHO_BASIS, allow_MPa=ALLOW, process=PROC,
                         process_basis=PROC_BASIS, practice=PRACTICE, hx_types=HX_TYPES, k_cf_frame3=K_CF),
        criterion=dict(name="КПД на тонну", formula="КПД, % / масса агрегата-заменителя, т",
                       mass_scope="двигатель + камера + теплообменники + корпус и рама модуля; генератор, редуктор, вход и выхлоп пакета не входят",
                       also=["КПД", "масса", "кВт/кг (мощность на валу / масса)", "фронт Парето КПД - масса"]),
        database=dict(engines=eng_db, packages=db["packages"], not_used=db["not_used"], legend=db["status_legend"]),
        calibration=cal, checks=chk, examples=ex,
        how_to_use=("import mass_model as MM; g = MM.size_from_passport(spec) (spec - как запись mass_db.json или MM.kv3_spec()); "
                    "ev = MM.evaluate(g) -> total_kg, total_lo_kg/total_hi_kg, kW_per_kg, kpd_per_t; свои блоки: g['rows'] - список строк "
                    "MM.block_* (группы rotor/static/hx/bearings/frame), g['hx'] - [dict(type, A_m2)], g['bearings'] - dict(kind); "
                    "MM.pareto_front(list_of_evaluations)"),
        not_done=[
            "нет открытой массы двигателя с центробежным компрессором и современной радиальной схемой (OPRA OP16, Kawasaki M1A, Capstone - только пакеты); "
            "K промышленной практики получен на осевых Solar и переносится на центробежные KV по допущению; единственная проверка - KG2-3 (1970-е, "
            "π и T4 приняты)",
            "доли блоков (ротор / корпуса / камера) не сверены с реальными машинами - известны только итоговые массы",
            "теплообменники не калиброваны по реальной массе (открыта только разность пакетов Mercury 50 - Centaur 50 как верхняя граница)",
            "рама модуля - оценка без калибровки (рамы пакетов Solar включают маслосистему, пуск, кожух - не отделить)",
            "печатные и литые КНР корпуса считаются стенкой своего способа, но множитель K - промышленной или авиационной практики: "
            "реальная масса печатной оболочки лежит между ними (показаны обе)",
            "двухвальный LM6000 сведён к средней линии на каждую группу ступеней; LM2500 и LM2500+ - силовая турбина принята 6-ступенчатой 3600 об/мин",
            "хосты geaerospace.com, wikipedia.org, vericor.com, motorsich.com, trid.trb.org, cloudfront Capstone не разрешены сетью сессии - "
            "их числа только по выдержке поиска",
        ],
        runtime_s=round((datetime.datetime.now() - t0).total_seconds(), 2),
    )
    return res


def main():
    res = build()
    with open(OUT_FILE, "w", encoding="utf-8") as f:
        json.dump(res, f, ensure_ascii=False, indent=1)
    cal = res["calibration"]
    for pr, v in cal["per_practice"].items():
        print(f"{pr:10s} n={v['n']} K={v['K']:.3f}  in-sample rms {v['insample']['rms_pct']:.1f} %  LOO rms {v['loo']['rms_pct']:.1f} %, "
              f"max {v['loo']['max_abs_pct']:.1f} %")
    for t in cal["table"]:
        print(f"  {t['id']:14s} real {t['mass_kg']:7.0f}  LOO {t['pred_loo_kg']:8.0f} ({t['err_loo_pct']:+.1f} %)  rotor {t['rotor_share_pct']:.0f} %")
    for ev in res["examples"]["items"]:
        print(f"  {ev['name'][:70]:70s} {ev['total_kg']:7.0f} кг [{ev['total_lo_kg']:.0f}-{ev['total_hi_kg']:.0f}]  "
              f"{ev.get('kW_per_kg', 0):.2f} кВт/кг  КПД/т {ev.get('kpd_per_t', 0):.2f}")
    print("записано", os.path.relpath(OUT_FILE, REPO))


if __name__ == "__main__":
    main()
