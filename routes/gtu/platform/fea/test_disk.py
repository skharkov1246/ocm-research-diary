#!/usr/bin/env python3
"""Проверки disk.py против аналитических решений (тонкий диск, плоское напряжённое
состояние). Запуск: python3 -m pytest -q test_disk.py  или  python3 test_disk.py"""
import math
import os
import sys
import tempfile

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import disk  # noqa: E402

STEEL = dict(name="сталь (условно)", density=7800.0, nu=0.3, E_GPa=[[20, 200.0]],
             alpha_mean_1e6=[[20, 12.0]], uts_MPa=[[20, 1000.0]])
RPM = 10000.0
W = 2 * math.pi * RPM / 60
TOL = 0.03          # требование задачи: расхождение не более 3 %


def _run(case):
    with tempfile.TemporaryDirectory() as d:
        res = disk.run_case(case, workdir=d, keep=False)
    assert res["quality"]["passed"], f"ворота равновесия: {res['quality']}"
    return res


def _mid_at(res, r_mm):
    mp = res["midplane"]
    k = int(np.argmin(np.abs(np.array(mp["r_mm"]) - r_mm)))
    assert abs(mp["r_mm"][k] - r_mm) < 1e-6, "в сетке нет узла на нужном радиусе"
    return {key: mp[key][k] for key in ("hoop_MPa", "radial_MPa", "vm_MPa", "T_C")}


def _rel(a, b):
    return abs(a - b) / abs(b)


def test_solid_disk_center():
    """Сплошной диск постоянной толщины: σθ(0) = σr(0) = (3+ν)/8·ρω²R²; σθ(R) = (1−ν)/4·ρω²R²."""
    R, nu, rho = 0.2, 0.3, 7800.0
    res = _run(dict(name="t_solid", material=STEEL, rpm=RPM,
                    profile=dict(r=[0.0, R], t=[0.004, 0.004]), mesh=dict(nr=40, nz=2)))
    s0 = (3 + nu) / 8 * rho * W ** 2 * R ** 2 / 1e6
    sR = (1 - nu) / 4 * rho * W ** 2 * R ** 2 / 1e6
    c = res["bore_mid"]
    print(f"  сплошной: σθ(0) = {c['hoop_MPa']:.3f} МПа, аналитика {s0:.3f} ({_rel(c['hoop_MPa'], s0):.2e}); "
          f"σr(0) = {c['radial_MPa']:.3f}; σθ(R) = {res['rim_mid']['hoop_MPa']:.3f}, аналитика {sR:.3f}")
    assert _rel(c["hoop_MPa"], s0) < TOL
    assert _rel(c["radial_MPa"], s0) < TOL
    assert _rel(res["rim_mid"]["hoop_MPa"], sR) < TOL
    # средние окружные: ∫σθ dA / A = ρω²R²/3 (точно из равновесия)
    sm = rho * W ** 2 * R ** 2 / 3 / 1e6
    assert _rel(res["section"]["hoop_mean_ip_MPa"], sm) < 0.005
    # запас по разрушению: sqrt(σв / σθ,ср)
    assert _rel(res["burst"]["margin"], math.sqrt(1000.0 / sm)) < 0.005


def test_disk_with_hole():
    """Диск с центральным отверстием: σθ(a) = (3+ν)/4·ρω²(b² + (1−ν)/(3+ν)·a²),
    σr,max = (3+ν)/8·ρω²(b−a)² при r = sqrt(ab), σr(a) = 0."""
    a, b, nu, rho = 0.05, 0.2, 0.3, 7800.0
    res = _run(dict(name="t_hole", material=STEEL, rpm=RPM,
                    profile=dict(r=[a, b], t=[0.004, 0.004]), mesh=dict(nr=60, nz=2)))
    sa = (3 + nu) / 4 * rho * W ** 2 * (b ** 2 + (1 - nu) / (3 + nu) * a ** 2) / 1e6
    srm = (3 + nu) / 8 * rho * W ** 2 * (b - a) ** 2 / 1e6
    m = _mid_at(res, math.sqrt(a * b) * 1e3)
    c = res["bore_mid"]
    print(f"  с отверстием: σθ(a) = {c['hoop_MPa']:.3f} МПа, аналитика {sa:.3f} ({_rel(c['hoop_MPa'], sa):.2e}); "
          f"σr(√ab) = {m['radial_MPa']:.3f}, аналитика {srm:.3f}; σr(a) = {c['radial_MPa']:.3f}")
    assert _rel(c["hoop_MPa"], sa) < TOL
    assert _rel(m["radial_MPa"], srm) < TOL
    assert abs(c["radial_MPa"]) < 0.01 * sa


def test_rim_load_lame():
    """Растягивающая нагрузка p на ободе, без вращения (Ламе): σθ(a) = 2p·b²/(b²−a²), σr(b) = p."""
    a, b, p = 0.05, 0.2, 100.0
    res = _run(dict(name="t_lame", material=STEEL, rpm=0.0, rim=dict(p_MPa=p),
                    profile=dict(r=[a, b], t=[0.004, 0.004]), mesh=dict(nr=60, nz=2)))
    sa = 2 * p * b ** 2 / (b ** 2 - a ** 2)
    print(f"  Ламе: σθ(a) = {res['bore_mid']['hoop_MPa']:.3f} МПа, аналитика {sa:.3f}; "
          f"σr(b) = {res['rim_mid']['radial_MPa']:.3f}, аналитика {p}")
    assert _rel(res["bore_mid"]["hoop_MPa"], sa) < TOL
    assert _rel(res["rim_mid"]["radial_MPa"], p) < TOL
    # равновесие полудиска: ∫σθ dA = p·b·t
    assert abs(res["section"]["err_ip_vs_equilibrium"]) < 0.005


def test_thermal_solid_disk():
    """Сплошной диск, T = T0 + ΔT·(r/b)², постоянные α, E:
    σθ(0) = αEΔT/4, σθ(b) = −αEΔT/2; средние окружные = 0."""
    b, dT = 0.2, 200.0
    rs = np.linspace(0.0, b, 81)
    res = _run(dict(name="t_therm", material=STEEL, rpm=0.0, T_ref_C=20.0,
                    T_C=dict(r=rs.tolist(), T=(20.0 + dT * (rs / b) ** 2).tolist()),
                    profile=dict(r=[0.0, b], t=[0.004, 0.004]), mesh=dict(nr=40, nz=2)))
    aE = 12e-6 * 200e9 / 1e6
    s0, sb = aE * dT / 4, -aE * dT / 2
    print(f"  термика: σθ(0) = {res['bore_mid']['hoop_MPa']:.3f} МПа, аналитика {s0:.3f}; "
          f"σθ(b) = {res['rim_mid']['hoop_MPa']:.3f}, аналитика {sb:.3f}; T(b) из .frd = {res['rim_mid']['T_C']:.2f}")
    assert _rel(res["bore_mid"]["hoop_MPa"], s0) < TOL
    assert _rel(res["rim_mid"]["hoop_MPa"], sb) < TOL
    assert abs(res["rim_mid"]["T_C"] - (20 + dT)) < 0.01        # температура дошла до ccx
    assert abs(res["section"]["hoop_mean_ip_MPa"]) < 0.005 * s0


def test_equilibrium_variable_profile():
    """Переменный профиль, масса лопаток, нагрузка на ободе, титан с T-зависимыми свойствами:
    ∫σθ dA (точки интегрирования) = ω²∫ρ_eff r² dA + p·r_o·t_o; узловой способ совпадает."""
    res = _run(dict(name="t_equil", material="ti64", rpm=20000.0,
                    profile=dict(r=[0.03, 0.08, 0.15], t=[0.06, 0.03, 0.01], hb=[0.0, 0.04, 0.02]),
                    blade_fill=0.08, rim=dict(mass_kg=2.0, r_cg=0.17),
                    T_C=dict(r=[0.03, 0.15], T=[50.0, 250.0]), mesh=dict(nr=60, nz=8)))
    s = res["section"]
    print(f"  равновесие: σθ ср = {s['hoop_mean_ip_MPa']:.3f} МПа (т. инт.), {s['hoop_mean_nodal_MPa']:.3f} "
          f"(узлы), {s['hoop_mean_equilibrium_MPa']:.3f} (равновесие); ошибки {s['err_ip_vs_equilibrium']:.1e}, "
          f"{s['err_nodal_vs_ip']:.1e}")
    assert abs(s["err_ip_vs_equilibrium"]) < 0.005
    assert abs(s["err_nodal_vs_ip"]) < 0.01


def test_pchip_profile_equilibrium_and_convergence():
    """Гладкий профиль (PCHIP через точки таблицы): равновесие выполняется, наибольшее
    напряжение по Мизесу на двух сетках (60×8 и 120×16) различается менее чем на 1 %."""
    base = dict(name="t_pchip", material="in713c", rpm=22000.0, profile_interp="pchip",
                profile=dict(r=[0.0, 0.04, 0.08, 0.12, 0.16, 0.2], t=[0.09, 0.085, 0.06, 0.035, 0.016, 0.008],
                             hb=[0.0, 0.07, 0.065, 0.055, 0.035, 0.024]),
                blade_fill=0.07, T_C=dict(r=[0.0, 0.1, 0.2], T=[620.0, 680.0, 860.0]))
    a = _run(dict(base, mesh=dict(nr=60, nz=8)))
    b = _run(dict(base, mesh=dict(nr=120, nz=16)))
    d = _rel(b["vm_max_node_MPa"], a["vm_max_node_MPa"])
    print(f"  PCHIP: σМизес max {a['vm_max_node_MPa']:.1f} / {b['vm_max_node_MPa']:.1f} МПа ({d:.2%}); "
          f"ворота {b['quality']['gate_equilibrium']:.1e}")
    assert d < 0.01


def test_lmp_in713c():
    """Подгонка Ларсона — Миллера воспроизводит точки данных; независимая точка 1700 °F / 30 ksi /
    76 ч предсказана в пределах 10 % по напряжению; минимум AMS 5391 лежит ниже типовой кривой."""
    m = disk.load_material("in713c")
    lmp = m["_lmp"]
    worst = 0.0
    for T, t, s, _ in m["lmp"]["data"]:
        worst = max(worst, _rel(lmp.sigma(T, t), s))
    T, t, s, _ = m["lmp"]["check"][0]
    chk = _rel(lmp.sigma(T, t), s)
    Ts, ts, ss, _ = m["lmp"]["spec_min"][0]
    life_typ, _ = lmp.life_h(ss, Ts)
    # обратимость: σ -> t -> σ
    life, _ = lmp.life_h(200.0, 900.0)
    back = lmp.sigma(900.0, life)
    print(f"  Ларсон — Миллер: худшая точка данных {worst:.1%}; проверочная 1700 °F/76 ч: "
          f"{lmp.sigma(T, t):.1f} МПа против {s} ({chk:.1%}); при 22 ksi/1800 °F типовой срок "
          f"{life_typ:.0f} ч (AMS: ≥30 ч); 200 МПа при 900 °C -> {life:.0f} ч")
    assert worst < 0.05
    assert chk < 0.10
    assert life_typ > ts
    assert _rel(back, 200.0) < 1e-6


def test_quality_gate_catches_corrupted_solution():
    """Ворота качества ловят испорченное решение: окружные напряжения в точках интегрирования
    одного столбца элементов завышены на 30 % — равновесие полудиска нарушается."""
    case = dict(name="t_gate", material=STEEL, rpm=RPM,
                profile=dict(r=[0.05, 0.2], t=[0.02, 0.004]), mesh=dict(nr=30, nz=4))
    mat = disk.load_material(STEEL)
    prof = disk.Profile(case["profile"])
    mesh = disk.build_mesh(prof, 30, 4)
    dens = disk.column_density_factor(prof, mesh["rc"], 0.0)
    with tempfile.TemporaryDirectory() as d:
        disk.run_case(case, workdir=d, keep=True)
        frd = disk.parse_frd(os.path.join(d, "disk.frd"))
        ips = disk.parse_dat_stress(os.path.join(d, "disk.dat"))
    args = (case, mat, prof, mesh, frd, None, W, 0.0, dens, 0.2, 0.004)
    good = disk.postprocess(*args[:5], ips, *args[6:])
    bad_ips = {k: ((v[0], v[1], v[2] * 1.3) + v[3:]) if k[0] <= 4 else v for k, v in ips.items()}
    bad = disk.postprocess(*args[:5], bad_ips, *args[6:])
    print(f"  ворота: исходное {good['quality']['gate_equilibrium']:.1e}, испорченное "
          f"{bad['quality']['gate_equilibrium']:.1e}")
    assert good["quality"]["passed"]
    assert not bad["quality"]["passed"]


def test_csv_profile_and_mesh_nodes_on_breakpoints():
    """Профиль из CSV читается; точки излома профиля попадают в узлы сетки."""
    with tempfile.TemporaryDirectory() as d:
        p = os.path.join(d, "p.csv")
        with open(p, "w", encoding="utf-8") as f:
            f.write("# r, t в метрах\nr,t,hb\n0.02,0.05,0\n0.07,0.03,0.02\n0.11,0.01,0.015\n")
        prof = disk.Profile(disk.read_profile_csv(p))
    rc = disk.radial_grid(prof, 20)
    for rb in prof.r:
        assert np.min(np.abs(rc - rb)) < 1e-12
    assert len(rc) - 1 == 20


def test_radial_grid_exact_count():
    """Элементов по радиусу ровно nr (раньше округление по участкам давало 119 или 121 вместо 120);
    при nr меньше числа участков — по одному на участок."""
    prof = disk.Profile(disk.example_cases()[0]["profile"], "pchip")      # 7 участков разной длины
    for nr in (8, 60, 120, 240):
        rc = disk.radial_grid(prof, nr)
        assert len(rc) - 1 == nr
        assert np.all(np.diff(rc) > 0)
        for rb in prof.r:
            assert np.min(np.abs(rc - rb)) < 1e-12
    assert len(disk.radial_grid(prof, 3)) - 1 == len(prof.r) - 1


def test_expansion_secant_zero():
    """*EXPANSION, ZERO=T_ref: КТР в таблице — средний (секущий) от T_ref (руководство CalculiX:
    ε = α(T)(T−T_ref) − α(T0)(T0−T_ref)). Свободный диск, равномерный нагрев до 500 °C, α(T) из таблицы:
    u_r(b) = α_ср(500)·480·b, напряжений нет. Тест с постоянным α этого не различает."""
    mat = dict(STEEL, alpha_mean_1e6=[[20, 10.0], [600, 16.0]])
    a, b = 0.05, 0.2
    res = _run(dict(name="t_secant", material=mat, rpm=0.0, T_C=500.0, T_ref_C=20.0,
                    profile=dict(r=[a, b], t=[0.01, 0.01]), mesh=dict(nr=20, nz=2)))
    ur = (10.0 + 6.0 * 480 / 580) * 1e-6 * 480 * b * 1e3
    print(f"  секущий КТР: u_r(b) = {res['rim_mid']['ur_mm']:.5f} мм, ожидается {ur:.5f}; "
          f"Мизес max {res['vm_max_node_MPa']:.2e} МПа")
    assert _rel(res["rim_mid"]["ur_mm"], ur) < 1e-3
    assert res["vm_max_node_MPa"] < 0.1
    assert res["burst"]["margin"] is None          # без вращения запас не определён
    assert any("нет вращения" in f for f in res["flags"])


def test_integration_point_mapping():
    """Порядок точек интегрирования в .dat: CAX8 развёрнут в C3D20, берётся средний слой (точки
    10–18, xi меняется быстрее). Проверка: σθ в точках совпадает с интерполяцией узлового поля
    в ту же точку, а при переставленном xi/eta — нет; в среднем слое σxz = σyz = 0."""
    case = dict(name="t_ipmap", material=STEEL, rpm=RPM, T_ref_C=20.0,
                T_C=dict(r=[0.05, 0.2], T=[20.0, 420.0]),
                profile=dict(r=[0.05, 0.12, 0.2], t=[0.04, 0.015, 0.006]), mesh=dict(nr=12, nz=4))
    prof = disk.Profile(case["profile"])
    mesh = disk.build_mesh(prof, 12, 4)
    with tempfile.TemporaryDirectory() as d:
        disk.run_case(case, workdir=d, keep=True)
        frd = disk.parse_frd(os.path.join(d, "disk.frd"))
        ips = disk.parse_dat_stress(os.path.join(d, "disk.dat"))
    err = {"принятый": [], "переставленный": []}
    for eid, conn, col in mesh["elems"]:
        st_n = np.array([frd["STRESS"][n][2] for n in conn])
        for jq in range(3):
            for iq in range(3):
                N, _, _ = disk.shape8(disk.G3[iq], disk.G3[jq])
                err["принятый"].append(ips[(eid, 9 + 3 * jq + iq + 1)][2] - N @ st_n)
                err["переставленный"].append(ips[(eid, 9 + 3 * iq + jq + 1)][2] - N @ st_n)
    rms = {k: float(np.sqrt(np.mean(np.square(v)))) / 1e6 for k, v in err.items()}
    shear = max(abs(v[4]) + abs(v[5]) for k, v in ips.items() if 10 <= k[1] <= 18)
    print(f"  точки интегрирования: СКО {rms['принятый']:.2f} МПа (принятый порядок), "
          f"{rms['переставленный']:.1f} МПа (переставленный); σxz+σyz в среднем слое {shear:.1e} Па")
    assert rms["принятый"] < 0.1 * rms["переставленный"]
    assert shear < 1.0


if __name__ == "__main__":
    fails = 0
    for name, fn in list(globals().items()):
        if name.startswith("test_") and callable(fn):
            try:
                fn()
                print(f"OK   {name}")
            except AssertionError as e:
                fails += 1
                print(f"FAIL {name}: {e}")
    print("итог:", "все проверки пройдены" if not fails else f"{fails} не пройдено")
    sys.exit(1 if fails else 0)
