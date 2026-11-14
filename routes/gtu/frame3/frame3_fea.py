#!/usr/bin/env python3
"""Типоразмер № 3 двигателя-заменителя: прочность колёс компрессора и дисков турбины в CalculiX.

Вход — расчётная точка frame3_point_results.json (14 950 об/мин): размеры колёс двух центробежных ступеней,
температуры за ступенями, средний диаметр и высоты лопаток трёх ступеней осевой турбины. Узел прочности —
../platform/fea/disk.py (осесимметричная модель CAX8, CalculiX ccx 2.21, материалы materials.json).

Что считается (частоты 100, 115 и 120 % от 14 950 об/мин):
  1) колесо ступени 1 — Ti-6Al-4V, двусторонний вход, профиль симметричный: обе грани диска — линии втулки
     двух половин колеса (четверть эллипса от втулки на входе до выхода), сердцевина до радиуса втулки входа —
     сплошная на всю осевую длину (колесо без сквозного отверстия, валы-цапфы на торцах);
  2) колесо ступени 2 — Ti-6Al-4V, односторонний вход: передняя грань — линия втулки, задняя — диск со ступицей
     (бобышка под цапфу); температура за ступенью сравнивается с пределом длительной службы титана;
  3) диски трёх ступеней турбины — литой IN-713C, лопатки на хвостовиках (ёлочный замок), масса лопаток,
     полок, ножек и зубьев замка — нагрузкой на живом ободе; профиль — равнопрочный диск Стодолы
     t(r) = t_o·exp(ρω²(r_b² − r²)/(2σ)) с ограничением толщины ступицы осевым шагом ступени.
     Диск ступени 1 (газ 1000 °C) — в трёх полях температур, обод: 1000 °C (температура газа, верхняя граница),
     относительная температура торможения в колесе из точки (основной вариант для неохлаждаемого диска),
     750 °C (продувка полостей воздухом из-за компрессора — оценка); к ступице — минус 200 К по квадрату радиуса.
Выход по каждому расчёту: наибольшие окружные и эквивалентные (Мизес) напряжения, запас по пределу текучести
(где есть σ0,2(T)), запас по разрушающей частоте (средние окружные напряжения, k = 1), для турбины — длительная
прочность по Ларсону — Миллеру на 25 000 и 50 000 ч (блок lmp материала IN-713C в materials.json). Проверки:
ворота равновесия disk.py, две сетки, сверка с аналитическим вращающимся диском постоянной толщины (как в
test_disk.py), сверка полотна с равнопрочным диском на изотермическом расчёте, напряжение в корне лопатки — с AN².

Все профили, поля температур и масса лопаток — оценка уровня облика (принятые величины — словарь A).
AWS не вызывается. Запуск:
  python3 frame3_fea.py [--work ПАПКА] [--workers 4]   -> frame3_fea_results.json (около 2 мин на 4 ядрах)
Рабочие файлы ccx (.inp, .frd, .dat) — в ПАПКЕ (по умолчанию <tmp>/frame3_fea), не в репозитории."""
import argparse
import datetime
import json
import math
import os
import platform
import sys
import tempfile
import time
from concurrent.futures import ProcessPoolExecutor

import numpy as np
import scipy

HERE = os.path.dirname(os.path.abspath(__file__))
GTU = os.path.dirname(HERE)
FEA_DIR = os.path.join(GTU, "platform", "fea")
sys.path.insert(0, FEA_DIR)
import disk  # noqa: E402

POINT_FILE = os.path.join(HERE, "frame3_point_results.json")
OUT = os.path.join(HERE, "frame3_fea_results.json")
REL_POINT = "routes/gtu/frame3/frame3_point_results.json"
REL_MAT = "routes/gtu/platform/fea/materials.json"

SPEEDS = (1.00, 1.15, 1.20)            # рабочая частота и разгон 115 и 120 %
MESH = dict(base=(80, 12), fine=(160, 24))   # элементов по радиусу x по толщине (CAX8)

# ------------------------------------------------------------------ принятые величины (все — оценка)
A = dict(
    T_ref_C=20.0,
    # колёса компрессора
    comp_blade_fill=0.06,      # доля объёма над диском, занятая лопатками (как в примерах disk.py)
    comp_t_rim_D2=0.02,        # осевая толщина диска на выходе колеса (под лопатками), доля D2
    comp2_boss_D2=0.10,        # ступень 2: осевая толщина ступицы за плоскостью выхода, доля D2
    comp2_boss_r_r1h=1.3,      # ступень 2: радиус бобышки под цапфу / радиус втулки на входе
    comp_bore_r=0.030,         # вариант со сквозным отверстием под центральную стяжку, м (диаметр 60 мм)
    # лопатки турбины на хвостовиках и обод диска
    turb_c_ax=(0.030, 0.035, 0.045),   # осевая хорда лопатки рабочего колеса ступеней 1-3, м
    turb_zweifel=0.8,          # коэффициент нагрузки Цвайфеля -> шаг решётки
    turb_k_area=0.15,          # площадь профиля в корне / хорда² (t_max/c около 0,22)
    turb_taper=(0.7, 0.6, 0.4),        # площадь профиля на периферии / в корне
    turb_t_plat=(0.003, 0.0035, 0.004),  # толщина полки, м
    turb_plat_len_cax=1.2,     # длина полки / осевая хорда
    turb_h_shank=(0.012, 0.013, 0.015),  # высота ножки (шейки) между полкой и замком, м
    turb_neck_s=0.30,          # толщина шейки / шаг
    turb_neck_len_cax=1.1,     # осевая длина шейки / осевая хорда
    turb_h_ft_cax=0.6,         # высота ёлочного замка / осевая хорда
    turb_dead_fill=0.9,        # заполнение кольца замка металлом (зубья лопаток + выступы диска)
    turb_w_rim_cax=1.15,       # ширина обода диска (= длина замка) / осевая хорда
    turb_band_w=0.5,           # радиальная высота обода / его ширина
    turb_fillet_w=0.3,         # радиальная длина перехода полотно -> обод / ширина обода
    turb_sigma_d_uts=0.5,      # расчётное напряжение равнопрочного полотна / σв(T обода) на 100 %
    turb_t_max_w=2.1,          # наибольшая толщина ступицы / ширина обода (осевой шаг ступени около 2,4·c_ax)
    turb_dT_hub=200.0,         # перепад температуры обод -> центр, К (квадрат радиуса)
    turb_dT_sens=(100.0, 300.0),
    turb_T_rim_cooled=750.0,   # обод при продувке полостей воздухом из-за компрессора, °C
    turb_rim_mass_sens=(0.75, 1.25),
    creep_hours=(25000.0, 50000.0),
)

# пределы службы титана — выдержки поиска 14.11.2026, первоисточники не открыты (сеть среды закрыта)
TI_LIMITS = [
    dict(alloy="Ti-6Al-4V", T_C=400, what="практический верхний предел по прочности (прочность падает выше)",
         src="https://secat.net/elevated-temperature-tensile-testing-of-ti-6al-4v-ti-64/", status="выдержка поиска"),
    dict(alloy="Ti-6Al-4V", T_C=400, what="около 400 °C длительной службы (в сравнении с Ti-6242 — 540 °C)",
         src="https://cnc.bozemetal.com/blog/ti-6242-titanium-alloy-properties-specifications-applications/",
         status="выдержка поиска, страница поставщика"),
    dict(alloy="Ti-6Al-4V", T_C=[300, 350], what="предел длительной службы по ползучести для α+β-сплавов",
         src="https://www.superstainlessalloys.com/info/fatigue-strength-and-creep-properties-of-ti-103255110.html",
         status="выдержка поиска, коммерческий источник, низкая достоверность"),
    dict(alloy="Ti-6Al-2Sn-4Zr-2Mo (Ti-6242)", T_C=[500, 540],
         what="диски и лопатки компрессора выше 350 °C; около 500 °C (обзор), до 538-540 °C (поставщики)",
         src="https://titanium.com/alloys/titanium-and-titanium-alloys/ti-6al-2sn-4zr-2mo/ ; "
             "https://www.nipponsteel.com/en/tech/report/nssmc/pdf/106-05.pdf",
         status="выдержка поиска"),
]


# ------------------------------------------------------------------ вход
def load_point(path=POINT_FILE):
    with open(path, encoding="utf-8") as f:
        return json.load(f)


def inputs_from_point(pt):
    s1, s2 = pt["stage1"], pt["stage2"]
    tz = str(pt["turbine"]["selected_z"])
    t = pt["turbine"]["options"][tz]
    keys = ("D2_m", "b2_m", "D1_hub_m", "D1_tip_m", "axial_length_m", "T01_C", "T02_C", "n_rpm", "U2_m_s",
            "c2_m_s", "W1t_m_s", "W2_W1t", "blades_full", "blades_splitter", "inlet")
    st = [dict(stage=s["stage"], h_rotor_in_m=s["h_rotor_in_m"], h_exit_m=s["h_exit_m"],
               T01_C=s["T01_C"], T_rel_rotor_C=s["T_rel_rotor_C"], AN2_m2rpm2=s["AN2_m2rpm2"]) for s in t["stages"]]
    return dict(stage1={k: s1[k] for k in keys}, stage2={k: s2[k] for k in keys},
                turbine=dict(z=t["z"], Dm=t["Dm"], U_mean=t["U_mean"], n=t["n"], psi=t["psi"],
                             phi_t=pt["assumptions"]["phi_t"], stages=st))


# ------------------------------------------------------------------ колёса компрессора: профиль
def _quarter_ellipse(r, r0, r2, height):
    """подъём линии от 0 (при r = r2, касательная радиальная) до height (при r = r0, касательная осевая)"""
    c = np.clip((r2 - np.asarray(r, float)) / (r2 - r0), 0.0, 1.0)
    return height * (1.0 - np.sqrt(1.0 - c * c))


def _comp_radii(r_first, r1h, r1t, r2, n=24):
    """точки таблицы: центр / отверстие, середина сердцевины, втулка входа; к выходу — сгущение у втулки,
    где линия втулки круче всего; радиус периферии входа — точка таблицы (излом линии покрывного контура)"""
    u = np.linspace(0.0, 1.0, n)
    rs = r1h + (r2 - r1h) * (1.0 - np.cos(0.5 * math.pi * u))
    pts = set(np.round(rs, 9)) | {round(r1t, 9)}
    if r_first < r1h - 1e-6:
        pts |= {round(r_first, 9), round(0.5 * (r_first + r1h), 9)}
    return np.array(sorted(pts))


def comp_profile(s, double, r_in=0.0):
    """Профиль колеса по размерам точки. Ось z — от плоскости выхода к входу. Линия втулки — четверть эллипса
    с полуосями (r2 − r1h) и осевой длиной колеса по втулке L (axial_length_m), покрывной контур — четверть
    эллипса (r2 − r1t) x (L − b2). Двустороннее колесо: две половины зеркально, общий диск толщиной t_rim на
    выходе. Одностороннее: задняя грань — t_rim на выходе, плавный рост (3u² − 2u³) до бобышки толщиной
    comp2_boss_D2·D2 при радиусе comp2_boss_r_r1h·r1h. Возврат — словарь профиля disk.py (r, t, zc, hb)."""
    r2, b2, L = s["D2_m"] / 2, s["b2_m"], s["axial_length_m"]
    r1h, r1t = s["D1_hub_m"] / 2, s["D1_tip_m"] / 2
    t_rim = A["comp_t_rim_D2"] * s["D2_m"]
    r = _comp_radii(r_in, r1h, r1t, r2)
    z_hub = _quarter_ellipse(r, r1h, r2, L)                          # от плоскости выхода по втулке
    z_shr = b2 + _quarter_ellipse(r, r1t, r2, L - b2)
    hb_side = np.where(r >= r1h - 1e-12, np.maximum(z_shr - z_hub, 0.0), 0.0)
    if double:
        t = t_rim + 2.0 * z_hub
        zc = np.zeros_like(r)
        hb = 2.0 * hb_side
    else:
        r_boss = A["comp2_boss_r_r1h"] * r1h
        t_boss = A["comp2_boss_D2"] * s["D2_m"]
        u = np.clip((r2 - r) / (r2 - r_boss), 0.0, 1.0)
        z_back = -(t_rim + (t_boss - t_rim) * (3 * u * u - 2 * u ** 3))
        z_front = z_hub
        t = z_front - z_back
        zc = 0.5 * (z_front + z_back)
        hb = hb_side
    return dict(r=r.tolist(), t=t.tolist(), zc=zc.tolist(), hb=hb.tolist())


def comp_T_field(s, r_lo):
    """условно: металл от T01 у оси до T02 (полная температура за ступенью — консервативно) на выходе,
    по квадрату радиуса (рост относительной температуры торможения в колесе ~ U²)"""
    r2 = s["D2_m"] / 2
    rs = np.linspace(r_lo, r2, 21)
    T = s["T01_C"] + (s["T02_C"] - s["T01_C"]) * (rs / r2) ** 2
    return dict(r=rs.tolist(), T=T.tolist())


def rotor_relative_T(s):
    """относительная температура торможения газа на выходе колеса T2 + W2²/(2cp), T2 = T02 − c2²/(2cp);
    cp воздуха — Cantera air.yaml при T2 (итерация)"""
    import cantera as ct
    g = ct.Solution("air.yaml")
    T02 = s["T02_C"] + 273.15
    c2, W2 = s["c2_m_s"], s["W2_W1t"] * s["W1t_m_s"]
    T2 = T02 - 300.0
    for _ in range(20):
        g.TP = T2, ct.one_atm
        T2 = T02 - c2 * c2 / (2 * g.cp_mass)
    g.TP = T2, ct.one_atm
    return dict(T2_static_C=T2 - 273.15, T0_rel_C=T2 + W2 * W2 / (2 * g.cp_mass) - 273.15, W2_m_s=W2,
                cp_J_kgK=g.cp_mass, cantera=ct.__version__)


# ------------------------------------------------------------------ турбина: лопатки на ободе, профиль диска
def blade_root_stress(rho, omega, r_h, h, taper):
    """растягивающее напряжение в корне профиля от вращения, площадь линейно от A_r (корень) до taper·A_r:
    σ = ρω²·[r_h·h·(1+τ)/2 + h²·(1/2 − (1−τ)/3)]; при τ = 1 — ρω²(r_t² − r_h²)/2 = ρω²·A_кольца/(2π)"""
    return rho * omega ** 2 * (r_h * h * (1 + taper) / 2 + h * h * (0.5 - (1 - taper) / 3))


def turbine_stage_geometry(k, tin, rho):
    """Лопатки ступени k (0..2) и обод диска: шаг по Цвайфелю, масса профиля, полки, шейки и кольца замка,
    радиус живого обода, ширина обода. Скорость — рабочая из точки."""
    st = tin["stages"][k]
    Dm, psi, phi = tin["Dm"], tin["psi"], tin["phi_t"]
    omega = tin["n"] * math.pi / 30
    rm = Dm / 2
    h = 0.5 * (st["h_rotor_in_m"] + st["h_exit_m"])
    r_h = rm - h / 2                                        # корень профиля = верх полки
    c_ax = A["turb_c_ax"][k]
    tan_in, tan_out = (psi - 1) / phi, 1 / phi              # относительные углы на входе и выходе колеса
    s_cax = A["turb_zweifel"] / (2 * math.cos(math.atan(tan_out)) ** 2 * (tan_in + tan_out))
    N = int(round(math.pi * Dm / (s_cax * c_ax)))
    s = math.pi * Dm / N
    stagger = math.atan(abs(tan_in - tan_out) / 2)
    chord = c_ax / math.cos(stagger)
    A_root = A["turb_k_area"] * chord ** 2
    tau = A["turb_taper"][k]
    m_af = rho * A_root * h * (1 + tau) / 2
    x_cg = h * (0.5 - (1 - tau) / 3) / ((1 + tau) / 2)
    t_p = A["turb_t_plat"][k]
    m_pl = rho * t_p * s * A["turb_plat_len_cax"] * c_ax
    h_sh = A["turb_h_shank"][k]
    m_sh = rho * A["turb_neck_s"] * s * A["turb_neck_len_cax"] * c_ax * h_sh
    r_ft_top = r_h - t_p - h_sh
    h_ft = A["turb_h_ft_cax"] * c_ax
    r_lr = r_ft_top - h_ft
    w = A["turb_w_rim_cax"] * c_ax
    m_dz = rho * A["turb_dead_fill"] * math.pi * (r_ft_top ** 2 - r_lr ** 2) * w
    r_dz = 2 / 3 * (r_ft_top ** 3 - r_lr ** 3) / (r_ft_top ** 2 - r_lr ** 2)
    parts = dict(airfoils=(N * m_af, r_h + x_cg), platforms=(N * m_pl, r_h - t_p / 2),
                 necks=(N * m_sh, r_h - t_p - h_sh / 2), dead_zone=(m_dz, r_dz))
    M = sum(m for m, _ in parts.values())
    r_cg = sum(m * r for m, r in parts.values()) / M
    sig_root = blade_root_stress(rho, omega, r_h, h, tau)
    sig_untap = blade_root_stress(rho, omega, rm - st["h_exit_m"] / 2, st["h_exit_m"], 1.0)   # = ρω²·AN²-кольцо/(2π)
    return dict(stage=k + 1, N_blades=N, pitch_m=s, s_over_cax=s_cax, c_ax_m=c_ax, chord_m=chord,
                stagger_deg=math.degrees(stagger), beta_in_deg=math.degrees(math.atan(tan_in)),
                beta_out_deg=math.degrees(math.atan(tan_out)), h_airfoil_m=h, r_hub_m=r_h,
                r_tip_m=r_h + h, A_root_mm2=A_root * 1e6, taper=tau, m_airfoil_kg=m_af, m_platform_kg=m_pl,
                m_neck_kg=m_sh, r_firtree_top_m=r_ft_top, h_firtree_m=h_ft, r_live_rim_m=r_lr, w_rim_m=w,
                rim_mass_parts_kg={k2: dict(m_kg=m, r_m=r) for k2, (m, r) in parts.items()},
                rim_mass_kg=M, rim_r_cg_m=r_cg,
                sigma_root_MPa=sig_root / 1e6,
                sigma_root_untapered_exit_h_MPa=sig_untap / 1e6,
                AN2_point_m2rpm2=st["AN2_m2rpm2"],
                sigma_AN2_untapered_MPa=rho * (math.pi / 30) ** 2 * st["AN2_m2rpm2"] / (2 * math.pi) / 1e6,
                T_rel_rotor_C=st["T_rel_rotor_C"], T01_C=st["T01_C"])


def turbine_profile(g, rho, omega, sigma_d):
    """Равнопрочный диск Стодолы (σr = σθ = σ_d, тонкий диск): t(r) = t_b·exp(ρω²(r_b² − r²)/(2σ_d)) от ступицы
    до низа обода r_b = r_o − band·w; t_b — из радиальной силы на r_b: масса на ободе + кольцо обода,
    Q = (M·ω²·r_cg/(2π) + ρω²·w·(r_o³ − r_b³)/3) / r_b, t_b = Q/σ_d (опора обода по окружности не засчитана —
    полотно толще). Толщина ступицы — не больше turb_t_max_w·w. Обод — толщина w от r_b до r_o, переход —
    сплайн PCHIP на длине turb_fillet_w·w."""
    r_o, w = g["r_live_rim_m"], g["w_rim_m"]
    r_b = r_o - A["turb_band_w"] * w
    Q = (g["rim_mass_kg"] * omega ** 2 * g["rim_r_cg_m"] / (2 * math.pi)
         + rho * omega ** 2 * w * (r_o ** 3 - r_b ** 3) / 3) / r_b
    t_b = Q / sigma_d
    t_max = A["turb_t_max_w"] * w
    r_f = r_b - A["turb_fillet_w"] * w
    rs = np.linspace(0.0, r_f, 15)
    t_cs = t_b * np.exp(rho * omega ** 2 * (r_b ** 2 - rs ** 2) / (2 * sigma_d))
    t = np.minimum(t_cs, t_max)
    band_t = max(w, t_b)
    r = list(rs) + [r_b, r_o]
    tt = list(t) + [band_t, band_t]
    return dict(r=[float(x) for x in r], t=[float(x) for x in tt]), dict(
        r_b_m=r_b, t_b_m=t_b, t_max_m=t_max, sigma_d_MPa=sigma_d / 1e6, Q_N_per_m=Q,
        t_cs_centre_m=float(t_b * math.exp(rho * omega ** 2 * r_b ** 2 / (2 * sigma_d))),
        capped=bool(t_cs[0] > t_max), band_t_m=band_t)


def turb_T_field(r_o, T_rim, dT):
    rs = np.linspace(0.0, r_o, 21)
    return dict(r=rs.tolist(), T=(T_rim - dT * (1 - (rs / r_o) ** 2)).tolist())


# ------------------------------------------------------------------ расчёт и разбор
def _run(job):
    key, case, wd = job
    os.makedirs(wd, exist_ok=True)
    res = disk.run_case(case, workdir=wd, keep=True)
    return key, res


def node_fields(case, wd):
    """узловые σ и T из .frd на той же сетке, что строит disk.run_case"""
    prof = disk.Profile(case["profile"], case.get("profile_interp", "linear"))
    m = case["mesh"]
    mesh = disk.build_mesh(prof, int(m["nr"]), int(m["nz"]))
    frd = disk.parse_frd(os.path.join(wd, "disk.frd"))
    S, TN = frd["STRESS"], {k: v[0] for k, v in frd.get("NDTEMP", {}).items()}
    out = []
    for k, (r, z, i, j) in mesh["nodes"].items():
        sr, sz, st, trz = S[k][:4]
        c, d = 0.5 * (sr + sz), math.hypot(0.5 * (sr - sz), trz)
        out.append((r, z, TN[k], float(disk.von_mises(sr, sz, st, trz)), max(st, c + d), sr, st))
    return np.array(out)


def creep_margins(case, wd, res, lmp):
    """Длительная прочность (Ларсон — Миллер, типовые свойства отдельно отлитых образцов, разрыв):
    допускаемое σ_LM(T, t) в каждом узле и наименьшее отношение σ_LM/σ — по Мизесу и по наибольшему главному
    растягивающему; по σθ,ср сечения при средней и наибольшей температуре."""
    a = node_fields(case, wd)
    T, vm, s1 = a[:, 2], a[:, 3], a[:, 4]
    out = {}
    for th in A["creep_hours"]:
        P = (T + 273.15) * (lmp.C + math.log10(th))
        sall = 10.0 ** lmp.log_sigma_v(P) * 1e6
        ex = (P < lmp.Pmin) | (P > lmp.Pmax)
        rv = sall / np.maximum(vm, 1.0)
        iv = int(np.argmin(rv))
        pos = s1 > 0
        r1 = np.where(pos, sall / np.where(pos, s1, 1.0), np.inf)
        i1 = int(np.argmin(r1))
        sec = res["section"]
        Tm, Tx = sec["T_mean_C"], float(T.max())
        s_sec = sec["hoop_mean_ip_MPa"]
        out[f"{int(th)}h"] = dict(
            vm=dict(ratio=float(rv[iv]), at=dict(r_mm=a[iv, 0] * 1e3, z_mm=a[iv, 1] * 1e3, T_C=T[iv],
                                                  sigma_MPa=vm[iv] / 1e6, sigma_allow_MPa=sall[iv] / 1e6),
                    extrapolated=bool(ex[iv])),
            s1=dict(ratio=float(r1[i1]), at=dict(r_mm=a[i1, 0] * 1e3, z_mm=a[i1, 1] * 1e3, T_C=T[i1],
                                                  sigma_MPa=s1[i1] / 1e6, sigma_allow_MPa=sall[i1] / 1e6),
                    extrapolated=bool(ex[i1])),
            section_mean_T=dict(T_C=Tm, sigma_MPa=s_sec, sigma_allow_MPa=lmp.sigma(Tm, th),
                                ratio=lmp.sigma(Tm, th) / s_sec if s_sec > 0 else None),
            section_hot_T=dict(T_C=Tx, sigma_MPa=s_sec, sigma_allow_MPa=lmp.sigma(Tx, th),
                               ratio=lmp.sigma(Tx, th) / s_sec if s_sec > 0 else None))
    return out


def summary(res, n_design):
    yr = res["yield_ratio_max"]
    keep = ("rpm", "U_tip", "r_in_mm", "r_out_mm", "rim_p_MPa", "vm_max_node_MPa", "vm_max_at", "vm_max_ip_MPa",
            "vm_max_ip_at", "hoop_max_MPa", "hoop_max_at", "bore_mid", "rim_mid", "section", "yield_ratio_at",
            "creep", "quality", "flags", "mesh")
    d = {k: res[k] for k in keep}
    d["speed_pct"] = round(100 * res["rpm"] / n_design, 3)
    d["yield_ratio_max"] = yr
    d["yield_margin"] = None if not yr else 1.0 / yr
    d["burst_margin"] = res["burst"]["margin"]
    d["burst_rpm"] = res["burst"]["rpm_burst"]
    d["burst_margin_vs_design_rpm"] = None if not res["burst"]["rpm_burst"] else res["burst"]["rpm_burst"] / n_design
    d["seconds"] = res["run"]["seconds"]
    return d


def conv(fine, base):
    def rel(a, b):
        return None if (a is None or b is None or b == 0) else a / b - 1
    c = dict(vm_max_node=rel(fine["vm_max_node_MPa"], base["vm_max_node_MPa"]),
             hoop_max=rel(fine["hoop_max_MPa"], base["hoop_max_MPa"]),
             hoop_mean=rel(fine["section"]["hoop_mean_ip_MPa"], base["section"]["hoop_mean_ip_MPa"]),
             centre_hoop=rel(fine["bore_mid"]["hoop_MPa"], base["bore_mid"]["hoop_MPa"]),
             burst_margin=rel(fine["burst"]["margin"], base["burst"]["margin"]),
             yield_ratio=rel(fine["yield_ratio_max"], base["yield_ratio_max"]))
    if fine.get("creep") and base.get("creep"):
        for k in ("local_min_h", "local_s1_min_h", "local_ip_min_h"):
            c["creep_" + k] = rel(fine["creep"][k], base["creep"][k])
    return c


def analytic_uniform(rho, nu, omega, R, p=0.0):
    """тонкий сплошной диск постоянной толщины: σθ(0) = σr(0) = (3+ν)/8·ρω²R² + p; σθ(R) = (1−ν)/4·ρω²R² + p"""
    q = rho * omega ** 2 * R ** 2
    return dict(hoop_centre_MPa=((3 + nu) / 8 * q + p) / 1e6, hoop_rim_MPa=((1 - nu) / 4 * q + p) / 1e6,
                radial_rim_MPa=p / 1e6)


# ------------------------------------------------------------------ сборка расчётов
def build(pt):
    tin = inputs_from_point(pt)
    n0 = tin["stage1"]["n_rpm"]
    assert tin["stage2"]["n_rpm"] == n0 and tin["turbine"]["n"] == n0
    ti, ni = disk.load_material("ti64"), disk.load_material("in713c")
    jobs, meta = [], {}

    def add(key, case, meshes=("base", "fine")):
        for mk in meshes:
            nr, nz = MESH[mk]
            jobs.append((key + (mk,), dict(case, mesh=dict(nr=nr, nz=nz), threads=1)))

    # --- колёса компрессора
    for part, sk, double in (("comp1", "stage1", True), ("comp2", "stage2", False)):
        s = tin[sk]
        for variant, r_in in (("solid", 0.0), ("bore", A["comp_bore_r"])):
            prof = comp_profile(s, double, r_in)
            Tf = comp_T_field(s, r_in)
            meta[(part, variant)] = dict(profile=prof, T_C=Tf)
            speeds = SPEEDS if variant == "solid" else (1.00, 1.20)
            for k in speeds:
                add((part, variant, k), dict(name=f"f3_{part}_{variant}_{int(round(k * 100))}", material="ti64",
                                             rpm=n0 * k, profile=prof, profile_interp="pchip",
                                             blade_fill=A["comp_blade_fill"], T_C=Tf, T_ref_C=A["T_ref_C"],
                                             burst_k=1.0))
    # --- диски турбины
    rho = ni["density"]
    om0 = n0 * math.pi / 30
    for k in range(3):
        g = turbine_stage_geometry(k, tin["turbine"], rho)
        T_base = g["T_rel_rotor_C"]
        sig_d = A["turb_sigma_d_uts"] * disk.tab_interp(ni["uts_MPa"], T_base)[0] * 1e6
        prof, pinfo = turbine_profile(g, rho, om0, sig_d)
        part = f"turb{k + 1}"
        rim = dict(mass_kg=g["rim_mass_kg"], r_cg=g["rim_r_cg_m"])
        tvars = {"T_rel": T_base}
        if k == 0:
            tvars = {"T_gas_1000": tin["turbine"]["stages"][0]["T01_C"], "T_rel": T_base,
                     "cooled_750": A["turb_T_rim_cooled"]}
        meta[part] = dict(geometry=g, profile=prof, profile_info=pinfo, T_variants=tvars)
        for tv, Trim in tvars.items():
            Tf = turb_T_field(g["r_live_rim_m"], Trim, A["turb_dT_hub"])
            meta[(part, tv)] = dict(T_C=Tf, T_rim_C=Trim, dT_hub_K=A["turb_dT_hub"])
            for sp in SPEEDS:
                add((part, tv, sp), dict(name=f"f3_{part}_{tv}_{int(round(sp * 100))}", material="in713c",
                                         rpm=n0 * sp, profile=prof, profile_interp="pchip", rim=rim, T_C=Tf,
                                         T_ref_C=A["T_ref_C"], burst_k=1.0))
        # чувствительность: перепад температуры к центру и масса на ободе (100 %, основной вариант)
        if k == 0:
            for dT in A["turb_dT_sens"]:
                Tf = turb_T_field(g["r_live_rim_m"], T_base, dT)
                meta[(part, f"dT{int(dT)}")] = dict(T_C=Tf, T_rim_C=T_base, dT_hub_K=dT)
                add((part, f"dT{int(dT)}", 1.0), dict(name=f"f3_{part}_dT{int(dT)}", material="in713c", rpm=n0,
                                                       profile=prof, profile_interp="pchip", rim=rim, T_C=Tf,
                                                       T_ref_C=A["T_ref_C"], burst_k=1.0), meshes=("fine",))
        if k in (0, 2):
            Tf = turb_T_field(g["r_live_rim_m"], T_base, A["turb_dT_hub"])
            for f in A["turb_rim_mass_sens"]:
                add((part, f"mass{f}", 1.0), dict(name=f"f3_{part}_mass{f}", material="in713c", rpm=n0,
                                                   profile=prof, profile_interp="pchip",
                                                   rim=dict(mass_kg=g["rim_mass_kg"] * f, r_cg=g["rim_r_cg_m"]),
                                                   T_C=Tf, T_ref_C=A["T_ref_C"], burst_k=1.0), meshes=("fine",))
        # изотермический расчёт для сверки полотна с равнопрочным диском (100 %, без нагрева)
        add((part, "iso", 1.0), dict(name=f"f3_{part}_iso", material="in713c", rpm=n0, profile=prof,
                                     profile_interp="pchip", rim=rim, T_C=A["T_ref_C"], T_ref_C=A["T_ref_C"],
                                     burst_k=1.0), meshes=("fine",))
    # --- аналитика: тонкий сплошной диск постоянной толщины, тот же материал, радиус и частота (как test_disk.py)
    for part, mat, R, p in (("comp1", "ti64", tin["stage1"]["D2_m"] / 2, 0.0),
                            ("comp2", "ti64", tin["stage2"]["D2_m"] / 2, 0.0),
                            ("turb1", "in713c", meta["turb1"]["geometry"]["r_live_rim_m"], 100.0),
                            ("turb2", "in713c", meta["turb2"]["geometry"]["r_live_rim_m"], 150.0),
                            ("turb3", "in713c", meta["turb3"]["geometry"]["r_live_rim_m"], 200.0)):
        case = dict(name=f"f3_{part}_analytic", material=mat, rpm=n0, T_C=A["T_ref_C"], T_ref_C=A["T_ref_C"],
                    profile=dict(r=[0.0, R], t=[0.02 * R, 0.02 * R]), burst_k=1.0)
        if p:
            case["rim"] = dict(p_MPa=p)
        jobs.append(((part, "analytic", 1.0, "ref"), dict(case, mesh=dict(nr=40, nz=2), threads=1)))
        meta[(part, "analytic")] = dict(R=R, p_MPa=p, material=mat)
    return tin, meta, jobs


def run_all(work, workers):
    pt = load_point()
    tin, meta, jobs = build(pt)
    os.makedirs(work, exist_ok=True)
    full = []
    for key, case in jobs:
        tag = "_".join(str(x) for x in key)
        full.append((key, case, os.path.join(work, tag)))
    t0 = time.time()
    with ProcessPoolExecutor(max_workers=workers) as ex:
        out = dict(ex.map(_run, full))
    wall = time.time() - t0
    cases = {key: (case, wd) for key, case, wd in full}
    return pt, tin, meta, out, cases, wall


def assemble(pt, tin, meta, out, cases, wall):
    n0 = tin["stage1"]["n_rpm"]
    ti, ni = disk.load_material("ti64"), disk.load_material("in713c")
    lmp = ni["_lmp"]
    om0 = n0 * math.pi / 30
    gates = []

    def get(key, mk="fine"):
        r = out[key + (mk,)] if mk else out[key]
        gates.append(r["quality"]["gate_equilibrium"])
        return r

    parts = {}
    # --- компрессор
    for part, sk in (("comp1", "stage1"), ("comp2", "stage2")):
        s = tin[sk]
        e = dict(stage=sk, material="ti64 (Ti-6Al-4V, отожжённый, деформируемый — materials.json)",
                 inlet=s["inlet"], D2_m=s["D2_m"], b2_m=s["b2_m"], D1_hub_m=s["D1_hub_m"], D1_tip_m=s["D1_tip_m"],
                 axial_length_m=s["axial_length_m"], blades=f"{s['blades_full']} + {s['blades_splitter']}",
                 t_rim_m=A["comp_t_rim_D2"] * s["D2_m"], blade_fill=A["comp_blade_fill"],
                 T01_C=s["T01_C"], T02_C=s["T02_C"], rotor_exit=rotor_relative_T(s), variants={})
        if part == "comp2":
            e["boss_t_m"] = A["comp2_boss_D2"] * s["D2_m"]
            e["boss_r_m"] = A["comp2_boss_r_r1h"] * s["D1_hub_m"] / 2
        for variant in ("solid", "bore"):
            mt = meta[(part, variant)]
            speeds = SPEEDS if variant == "solid" else (1.00, 1.20)
            v = dict(profile=mt["profile"], T_C=mt["T_C"], cases={})
            for k in speeds:
                f, b = get((part, variant, k)), get((part, variant, k), "base")
                c = summary(f, n0)
                c["creep"] = "нет данных: в materials.json у ti64 блок lmp = null — ползучесть титана не оценена"
                c["convergence_fine_vs_base"] = conv(f, b)
                v["cases"][f"{int(round(k * 100))}"] = c
            e["variants"][variant] = v
        parts[part] = e
    # температура ступени 2 и титан
    s2 = parts["comp2"]
    s2["titanium_temperature"] = dict(
        T02_C=s2["T02_C"], T0_rel_rotor_exit_C=s2["rotor_exit"]["T0_rel_C"],
        ys_table_max_C=ti["ys_MPa"][-1][0], limits=TI_LIMITS,
        verdict=("выше допустимой для длительной службы Ti-6Al-4V: полная температура за ступенью "
                 f"{s2['T02_C']:.0f} °C больше практического предела около 400 °C и предела по ползучести "
                 "300–350 °C (выдержки поиска, первоисточники не открыты); металл колеса у выхода — между "
                 f"относительной температурой торможения {s2['rotor_exit']['T0_rel_C']:.0f} °C и T02 (задняя "
                 "полость). Ползучесть Ti-6Al-4V не оценена — нет данных в materials.json."),
        proposal=("Ti-6Al-2Sn-4Zr-2Mo-0,1Si (Ti-6242S, поковка): диски и колёса компрессора выше 350 °C, предел "
                  "службы около 500–540 °C (выдержки поиска), плотность около 4,54 г/см³ (против 4,43 у "
                  "Ti-6Al-4V) — напряжения от вращения выше на 2,5 %. Прочность и ползучесть Ti-6242S не "
                  "посчитаны: свойств в materials.json нет, сайты первоисточников из среды закрыты. Запасной "
                  "путь — оставить Ti-6Al-4V и снизить металл у выхода (продувка задней полости воздухом "
                  "из-за ступени 1) — нужен тепловой расчёт."))
    # --- турбина
    for k in range(3):
        part = f"turb{k + 1}"
        mt = meta[part]
        g = mt["geometry"]
        e = dict(material="in713c (IN-713C литой — materials.json)", geometry=g, profile=mt["profile"],
                 profile_info=mt["profile_info"], dT_hub_K=A["turb_dT_hub"], variants={})
        bl = dict(sigma_root_MPa=g["sigma_root_MPa"], T_C=g["T_rel_rotor_C"])
        life, ex = lmp.life_h(g["sigma_root_MPa"], g["T_rel_rotor_C"])
        bl["creep_life_h"] = disk._fin(life)
        bl["creep_life_extrapolated"] = ex
        for th in A["creep_hours"]:
            sa = lmp.sigma(g["T_rel_rotor_C"], th)
            P = (g["T_rel_rotor_C"] + 273.15) * (lmp.C + math.log10(th))
            bl[f"sigma_allow_{int(th)}h_MPa"] = sa
            bl[f"ratio_{int(th)}h"] = sa / g["sigma_root_MPa"]
            bl[f"extrapolated_{int(th)}h"] = bool(P < lmp.Pmin or P > lmp.Pmax)
        bl["sigma_root_115_MPa"] = g["sigma_root_MPa"] * 1.15 ** 2
        bl["sigma_root_120_MPa"] = g["sigma_root_MPa"] * 1.20 ** 2
        bl["uts_at_T_MPa"] = disk.tab_interp(ni["uts_MPa"], g["T_rel_rotor_C"])[0]
        bl["note"] = ("профиль лопатки — линейное убывание площади (оценка); металл профиля — при относительной "
                      "температуре торможения на входе в колесо (неохлаждаемая лопатка, без неравномерности поля "
                      "за камерой); замок и шейка не считались (нужна трёхмерная модель)")
        e["blade"] = bl
        for tv in mt["T_variants"]:
            tm = meta[(part, tv)]
            v = dict(T_rim_C=tm["T_rim_C"], T_C=tm["T_C"], cases={})
            for sp in SPEEDS:
                key = (part, tv, sp)
                f, b = get(key), get(key, "base")
                c = summary(f, n0)
                c["yield_ratio_ref_20C"] = f["vm_max_node_MPa"] / ni["ys_MPa"][0][1]
                c["convergence_fine_vs_base"] = conv(f, b)
                if abs(sp - 1.0) < 1e-12:
                    c["creep_margins"] = creep_margins(cases[key + ("fine",)][0], cases[key + ("fine",)][1], f, lmp)
                v["cases"][f"{int(round(sp * 100))}"] = c
            e["variants"][tv] = v
        sens = {}
        for kk in list(out):
            if kk[0] == part and kk[1] not in mt["T_variants"] and kk[1] not in ("analytic", "iso"):
                f = get(kk[:3])
                c = summary(f, n0)
                c["creep_margins"] = creep_margins(cases[kk][0], cases[kk][1], f, lmp)
                if kk[1].startswith("dT"):
                    c["dT_hub_K"] = meta[(part, kk[1])]["dT_hub_K"]
                else:
                    c["rim_mass_factor"] = float(kk[1][4:])
                sens[kk[1]] = c
        if sens:
            e["sensitivity_100pct"] = sens
        # сверка полотна с равнопрочным диском (изотермический расчёт)
        iso = get((part, "iso", 1.0))
        mp = iso["midplane"]
        rr = np.array(mp["r_mm"]) / 1e3
        pr = mt["profile_info"]
        web = (rr > 0.15 * pr["r_b_m"]) & (rr < pr["r_b_m"] - (A["turb_fillet_w"] + 0.15) * g["w_rim_m"])
        tprof = disk.Profile(mt["profile"], "pchip")
        web &= np.array([tprof.at(x)[0] < 0.98 * pr["t_max_m"] for x in rr])
        sh, sr_ = np.array(mp["hoop_MPa"])[web], np.array(mp["radial_MPa"])[web]
        e["constant_stress_check"] = dict(
            sigma_d_MPa=pr["sigma_d_MPa"], n_points=int(web.sum()),
            r_range_mm=[float(rr[web].min() * 1e3), float(rr[web].max() * 1e3)] if web.any() else None,
            hoop_MPa=[float(sh.min()), float(sh.max())] if web.any() else None,
            radial_MPa=[float(sr_.min()), float(sr_.max())] if web.any() else None,
            dev_max=float(max(np.abs(sh / pr["sigma_d_MPa"] - 1).max(), np.abs(sr_ / pr["sigma_d_MPa"] - 1).max()))
            if web.any() else None,
            vm_max_node_MPa=iso["vm_max_node_MPa"], vm_max_at=iso["vm_max_at"],
            centre_hoop_MPa=iso["bore_mid"]["hoop_MPa"], quality=iso["quality"],
            note=("изотермический расчёт (20 °C) на 100 %: в полотне, где профиль следует формуле Стодолы "
                  "(без ступицы, ограниченной по толщине, и без перехода к ободу), σr и σθ должны быть близки к "
                  "σ_d; отличие — от толщины (не тонкий диск), сплайна и опоры обода"))
        parts[part] = e
    # --- аналитика
    analytic = {}
    for part in ("comp1", "comp2", "turb1", "turb2", "turb3"):
        mt = meta[(part, "analytic")]
        m = ti if mt["material"] == "ti64" else ni
        r = out[(part, "analytic", 1.0, "ref")]
        gates.append(r["quality"]["gate_equilibrium"])
        an = analytic_uniform(m["density"], m["nu"], om0, mt["R"], mt["p_MPa"] * 1e6)
        fe = dict(hoop_centre_MPa=r["bore_mid"]["hoop_MPa"], hoop_rim_MPa=r["rim_mid"]["hoop_MPa"],
                  radial_rim_MPa=r["rim_mid"]["radial_MPa"])
        analytic[part] = dict(R_m=mt["R"], material=mt["material"], rpm=n0, rim_p_MPa=mt["p_MPa"],
                              fea=fe, analytic=an,
                              dev=dict(hoop_centre=fe["hoop_centre_MPa"] / an["hoop_centre_MPa"] - 1,
                                       hoop_rim=fe["hoop_rim_MPa"] / an["hoop_rim_MPa"] - 1),
                              quality=r["quality"], mesh=dict(nr=40, nz=2), thickness_over_R=0.02)
    return parts, analytic, gates


def table_rows(parts):
    rows = []
    for part in ("comp1", "comp2"):
        for variant, v in parts[part]["variants"].items():
            for sp, c in v["cases"].items():
                rows.append(dict(part=part, variant=variant, speed_pct=int(sp), rpm=c["rpm"], U_tip=c["U_tip"],
                                 vm_max_MPa=c["vm_max_node_MPa"], vm_at_r_mm=c["vm_max_at"]["r_mm"],
                                 hoop_max_MPa=c["hoop_max_MPa"], yield_margin=c["yield_margin"],
                                 burst_margin=c["burst_margin"], burst_rpm=c["burst_rpm"]))
    for k in (1, 2, 3):
        e = parts[f"turb{k}"]
        for tv, v in e["variants"].items():
            for sp, c in v["cases"].items():
                row = dict(part=f"turb{k}", variant=tv, T_rim_C=v["T_rim_C"], speed_pct=int(sp), rpm=c["rpm"],
                           U_tip=c["U_tip"], vm_max_MPa=c["vm_max_node_MPa"], vm_at_r_mm=c["vm_max_at"]["r_mm"],
                           hoop_max_MPa=c["hoop_max_MPa"], yield_margin=None, burst_margin=c["burst_margin"],
                           burst_rpm=c["burst_rpm"])
                if "creep_margins" in c:
                    cm = c["creep_margins"]
                    for th in ("25000h", "50000h"):
                        row[f"creep_s1_ratio_{th}"] = cm[th]["s1"]["ratio"]
                        row[f"creep_vm_ratio_{th}"] = cm[th]["vm"]["ratio"]
                    row["creep_life_s1_h"] = c["creep"]["local_s1_min_h"]
                    row["creep_life_vm_h"] = c["creep"]["local_min_h"]
                rows.append(row)
    return rows


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--work", default=os.path.join(tempfile.gettempdir(), "frame3_fea"),
                    help="папка рабочих файлов ccx (не в репозитории)")
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--out", default=OUT)
    a = ap.parse_args()
    pt, tin, meta, out, cases, wall = run_all(a.work, a.workers)
    parts, analytic, gates = assemble(pt, tin, meta, out, cases, wall)
    ver = disk.ccx_version()
    n0 = tin["stage1"]["n_rpm"]
    import cantera
    s1, s2 = tin["stage1"], tin["stage2"]
    doc = dict(
        title="Типоразмер № 3: прочность колёс компрессора и дисков турбины (CalculiX, осесимметрия)",
        method=("Осесимметричная КЭ-модель CAX8 (routes/gtu/platform/fea/disk.py, CalculiX): вращение (*DLOAD "
                "CENTRIF), поле температур T(r), свойства E(T), α(T) из materials.json; лопатки колёс компрессора — "
                "массой в плотности столбца под ними, лопатки турбины — массой на живом ободе (давление на ободе). "
                "Наибольшие окружные и эквивалентные напряжения — узловые (.frd) и в точках интегрирования (.dat); "
                "запас по пределу текучести — наибольшее σ_Мизес/σ0,2(T) по узлам; запас по разрушающей частоте — "
                "средние окружные напряжения, n_разр/n = sqrt(k·∫σв(T)dA/∫σθdA), k = 1 (верхняя оценка); "
                "длительная прочность IN-713C — параметр Ларсона — Миллера P = T[K]·(20 + lg t), разрыв, типовые "
                "свойства отдельно отлитых образцов. Профили колёс компрессора — по линиям втулки и покрывного "
                "контура из размеров точки (четверть эллипса); диски турбины — равнопрочный диск Стодолы, "
                "ступица ограничена осевым шагом ступени. Две сетки 80x12 и 160x24 элементов, итог — по второй."),
        software=(f"CalculiX ccx ({ver}); python {platform.python_version()}, numpy {np.__version__}, scipy "
                  f"{scipy.__version__}, cantera {cantera.__version__} (cp воздуха); disk.py узла FEA"),
        date="2026-11-14",
        status=("оценка уровня облика: профили колёс и дисков, поля температур, масса и форма лопаток турбины — "
                "принятые величины (словарь assumptions), не конструкция и не тепловой расчёт; свойства "
                "материалов — выдержки поиска (materials.json, статусы); для проектных выводов нужны чертёж, "
                "тепловой расчёт (OpenFOAM), трёхмерная модель лопаток и замков, свойства литого сечения и "
                "разгонные испытания"),
        inputs=[
            dict(file=REL_POINT, field="stage1: D2_m, b2_m, D1_hub_m, D1_tip_m, axial_length_m, T01_C, T02_C, n_rpm, "
                 "U2_m_s, c2_m_s, W1t_m_s, W2_W1t, blades_full, blades_splitter, inlet", used=s1),
            dict(file=REL_POINT, field="stage2: те же поля", used=s2),
            dict(file=REL_POINT, field="turbine.options[selected_z].Dm, U_mean, n, psi; stages[].h_rotor_in_m, "
                 "h_exit_m, T01_C, T_rel_rotor_C, AN2_m2rpm2; assumptions.phi_t", used=tin["turbine"]),
            dict(file=REL_MAT, field="ti64, in713c (density, nu, E_GPa, alpha_mean_1e6, uts_MPa, ys_MPa, lmp)",
                 used=dict(ti64_ys_MPa=disk.load_material("ti64")["ys_MPa"],
                           ti64_uts_MPa=disk.load_material("ti64")["uts_MPa"],
                           in713c_uts_MPa=disk.load_material("in713c")["uts_MPa"],
                           in713c_ys_MPa=disk.load_material("in713c")["ys_MPa"],
                           in713c_lmp_C=disk.load_material("in713c")["lmp"]["C"])),
            dict(file="routes/gtu/platform/fea/disk.py", field="run_case, LMP, tab_interp", used="без изменений"),
        ],
        assumptions=dict(A, note="все величины — оценка (облик); смысл каждой — комментарии словаря A в frame3_fea.py"),
        speeds=dict(n_design_rpm=n0, factors=list(SPEEDS), rpm=[n0 * k for k in SPEEDS]),
        mesh=dict(MESH, used_for_results="fine"),
        checks=dict(
            equilibrium_gate_max=max(gates), equilibrium_gate_limit=1e-3,
            all_gates_passed=bool(max(gates) < 1e-3),
            analytic_uniform_disk=analytic,
            prerequisites=("перед расчётом: .claude/skills/calculix-rotor/scripts/ccx_selfcheck.py — all_ok: true "
                           "(7 из 7); routes/gtu/platform/fea/test_disk.py — все 12 проверок пройдены (прогоны "
                           "14.11.2026)")),
        parts=parts,
        table=table_rows(parts),
        run=dict(cases=len(out), wall_s=round(wall, 1), workers=a.workers, ccx_threads=1,
                 workdir="вне репозитория (--work); .inp/.frd/.dat не сохраняются в git"),
    )
    with open(a.out, "w", encoding="utf-8") as f:
        json.dump(doc, f, ensure_ascii=False, indent=1, default=disk._json_default)
    print(f"расчётов {len(out)}, {wall:.0f} с; наибольшая невязка равновесия {max(gates):.1e}")
    for r in doc["table"]:
        cr = ""
        if "creep_s1_ratio_25000h" in r:
            cr = (f" LM25k s1 {r['creep_s1_ratio_25000h']:.2f} vm {r['creep_vm_ratio_25000h']:.2f};"
                  f" LM50k s1 {r['creep_s1_ratio_50000h']:.2f}")
        ym = "—" if r["yield_margin"] is None else f"{r['yield_margin']:.2f}"
        print(f"{r['part']:6s} {r['variant']:10s} {r['speed_pct']:3d}% U={r['U_tip']:.0f} vm={r['vm_max_MPa']:.0f} "
              f"(r={r['vm_at_r_mm']:.0f}) hoop={r['hoop_max_MPa']:.0f} yield×{ym} burst×{r['burst_margin']:.3f}{cr}")
    print("записано:", a.out)


if __name__ == "__main__":
    main()
