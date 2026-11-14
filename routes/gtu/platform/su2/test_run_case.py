"""Тесты run_case.py без SU2: заполнение шаблона, оценка давления, разбор истории и журнала.

Все истории и журналы здесь — придуманные корпуса (FIXTURE), не вывод SU2.
Запуск: python -m pytest -q test_run_case.py
"""
import json
import math
import os
import sys

import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import run_case as rc  # noqa: E402

GEOM_JSON = os.path.join(HERE, "..", "geom", "out", "core130_coarse", "impeller_core130_coarse.json")

D1_FIXTURE = dict(T01=288.15, p01=101325.0, w=178870.0, PR=4.25, U2=505.6, psi=0.70, r2=0.2175, m=9.95, b2=0.0234)


def _params():
    p = dict(rc.NUMERICS)
    p.update(RESTART="NO", SOLUTION_FILENAME="solution.dat", P01="101325", T01="288.15", OMEGA_X="0.0",
             OMEGA_Y="0.0", OMEGA_Z="2324.46", PER_MASTER="periodic_1", PER_SLAVE="periodic_2", ANGLE="40",
             P_OUT="300000", ITER=10, MUSCL="NO", HISTORY_NAME="history_x", RESTART_NAME="restart_x.dat",
             VOLUME_NAME="vol_x", MESH="mesh.su2", FS_P="90000", FS_T="279", RAMP_ROT="YES", RAMP_ROT_INI="0", RAMP_ROT_FREQ=1,
             RAMP_ROT_ITER=100, RAMP_OUT="YES", RAMP_OUT_INI="101325", RAMP_OUT_FREQ=1, RAMP_OUT_ITER=100)
    return p


def test_template_fully_substituted_and_keys_unique():
    txt = rc.build_cfg(open(rc.TEMPLATE, encoding="utf-8").read(), _params())
    assert "${" not in txt
    keys = [ln.split("=")[0].strip() for ln in txt.splitlines() if ln and not ln.startswith("%") and "=" in ln]
    assert len(keys) == len(set(keys)), "ключ повторён"
    for k in ("MARKER_PERIODIC", "MARKER_GILES", "MARKER_SHROUD", "ROTATION_RATE", "MARKER_TURBOMACHINERY"):
        assert k in keys
    # маркеры сетки узла geom
    for m in ("inlet", "outlet", "hub", "shroud", "blade", "periodic_1", "periodic_2"):
        assert m in txt


def test_template_missing_param_fails():
    p = _params()
    del p["P_OUT"]
    with pytest.raises(KeyError):
        rc.build_cfg(open(rc.TEMPLATE, encoding="utf-8").read(), p)


def test_estimate_p_out_limits_and_monotonic():
    g, R = 1.4, 287.058
    r_out = 1.3 * D1_FIXTURE["r2"]
    p0 = rc.estimate_p_out(D1_FIXTURE, 0.0, r_out, D1_FIXTURE["b2"], g, R)
    assert p0 == pytest.approx(D1_FIXTURE["p01"], rel=1e-9)
    ps = [rc.estimate_p_out(D1_FIXTURE, f, r_out, D1_FIXTURE["b2"], g, R) for f in (0.25, 0.5, 0.75, 1.0)]
    assert all(b > a for a, b in zip(ps, ps[1:]))
    # статическое давление на выходе ниже полного PR*p01 и выше p01
    assert D1_FIXTURE["p01"] < ps[-1] < D1_FIXTURE["PR"] * D1_FIXTURE["p01"]
    # независимая проверка последней точки: полное давление по PR, скорость с выхода -> статика
    cp = g * R / (g - 1)
    T02 = D1_FIXTURE["T01"] + D1_FIXTURE["w"] / cp
    p02 = D1_FIXTURE["PR"] * D1_FIXTURE["p01"]
    ct = D1_FIXTURE["psi"] * D1_FIXTURE["U2"] * D1_FIXTURE["r2"] / r_out
    # подставим найденное давление обратно: неразрывность должна выполниться
    p = ps[-1]
    T = T02 * (p / p02) ** ((g - 1) / g)
    cm = math.sqrt(max(2 * cp * (T02 - T) - ct ** 2, 0))
    m = p / (R * T) * cm * 2 * math.pi * r_out * D1_FIXTURE["b2"]
    assert m == pytest.approx(D1_FIXTURE["m"], rel=1e-4)


HIST_FIXTURE = (
    '"Time_Iter","Outer_Iter","Inner_Iter","rms[Rho]","rms[RhoE]","Avg_Massflow(inlet)","Avg_Massflow(outlet)",'
    '"Avg_TotalPress(inlet)","Avg_TotalPress(outlet)","Avg_TotalTemp(inlet)","Avg_TotalTemp(outlet)","PRTT","TotTotEff"\n'
)


def _write_hist(tmp_path, rows, truncated=False):
    f = tmp_path / "history_fixture.csv"
    txt = HIST_FIXTURE + "".join(",".join(str(x) for x in r) + "\n" for r in rows)
    if truncated:
        txt += "0,0,999,-1.0"  # недописанная строка (останов по лимиту)
    f.write_text(txt, encoding="utf-8")
    return str(f)


def test_analyse_history_numbers(tmp_path):
    rows = []
    for i in range(100):
        rows.append([0, 0, i, -3.0 - 0.02 * i if i > 10 else -3.0 + 0.1 * i, 1.0 - 0.01 * i,
                     1.10, -1.089, 101325.0, 405300.0, 288.15, 461.04, 3.9, 85.0])
    h = rc.analyse_history(_write_hist(tmp_path, rows, truncated=True), 1.4, 9, window=50)
    assert h["rows"] == 100
    r = h["residuals"]["rms[Rho]"]
    assert r["max"] == pytest.approx(-2.0)
    assert r["iter_max"] == 10
    assert r["last"] == pytest.approx(-3.0 - 0.02 * 99)
    assert r["slope_per_100_last_half"] == pytest.approx(-2.0, rel=1e-6)
    pf = h["performance"]
    assert pf["mass_imbalance"] == pytest.approx(0.011 / 1.10)
    assert pf["mass_sign_ok"] is True
    assert pf["m_wheel"] == pytest.approx(9.9)
    assert pf["PR_tt"] == pytest.approx(4.0)
    tr = 461.04 / 288.15
    assert pf["eta_tt_isentropic"] == pytest.approx((4.0 ** (0.4 / 1.4) - 1) / (tr - 1))
    assert pf["PRTT_su2"] == 3.9 and pf["TotTotEff_su2"] == 85.0
    assert pf["window_rel_change"]["m_in"] == 0.0


def test_reverse_flow_flagged(tmp_path):
    rows = [[0, 0, i, -3.0, 1.0, 0.1, 0.3, 101325.0, 110000.0, 288.15, 300.0, 1.1, 50.0] for i in range(5)]
    h = rc.analyse_history(_write_hist(tmp_path, rows), 1.4, 9)
    assert h["performance"]["mass_sign_ok"] is False
    c = rc.smoke_checks([{"status": "finished", "history": h}], 0.05)
    assert c["mass_within_tol"] is False and c["passed"] is False


def test_smoke_checks_divergence_fails(tmp_path):
    rows = [[0, 0, i, -3.0 - 0.01 * i, 1.0, 1.1, -1.1, 101325.0, 405300.0, 288.15, 461.0, 4.0, 85.0] for i in range(40)]
    h = rc.analyse_history(_write_hist(tmp_path, rows), 1.4, 9)
    ok = rc.smoke_checks([{"status": "finished", "history": h}], 0.05)
    assert ok["passed"] is True
    bad = rc.smoke_checks([{"status": "diverged", "history": h}], 0.05)
    assert bad["passed"] is False and bad["no_divergence"] is False


LOG_FIXTURE = """Matching the periodic boundary points for marker pair 1.
 Matched 6270 points with a max distance of: 1.70109e-16.
Inlet area for Row 1: 56.7511 cm^2.
Oulet area for Row 1: 46.2201 cm^2.
Recomputed number of blades for Row 1: 9.05246.
Error in "void CSolver::SetResidual_RMS(const CGeometry*, const CConfig*, bool)":
SU2 has diverged (NaN detected).
"""


def test_parse_log(tmp_path):
    f = tmp_path / "su2_fixture.log"
    f.write_text(LOG_FIXTURE, encoding="utf-8")
    o = rc.parse_log(str(f))
    assert o["periodic_matched_points"] == 6270
    assert o["periodic_max_distance"] == pytest.approx(1.70109e-16)
    assert o["bad_match"] == 0
    assert o["nan"] is True
    assert o["inlet_area_cm2"] == pytest.approx(56.7511)
    assert o["outlet_area_cm2"] == pytest.approx(46.2201)
    assert o["recomputed_blades"] == pytest.approx(9.05246)
    assert "NaN" in o["error"]


def test_check_keys(tmp_path):
    src = tmp_path / "Common" / "src"
    src.mkdir(parents=True)
    (src / "CConfig.cpp").write_text('addDoubleOption("CFL_NUMBER", x);\naddEnumOption("SOLVER", y);\n',
                                     encoding="utf-8")
    unknown, n = rc.check_keys("SOLVER= RANS\nCFL_NUMBER= 1\nMARKER_NRBC= (x)\n% COMMENT= 1\n", str(src / "CConfig.cpp"))
    assert unknown == ["MARKER_NRBC"] and n == 2
    assert rc.find_cconfig(str(tmp_path)) == str(src / "CConfig.cpp")


def test_mpi_cmd_root_flag():
    cmd = rc.mpi_cmd(4, "/x/SU2_CFD", "a.cfg")
    assert cmd[1:3] == ["-np", "4"] and cmd[-2:] == ["/x/SU2_CFD", "a.cfg"]
    if os.geteuid() == 0:
        assert "--allow-run-as-root" in cmd
    assert rc.mpi_cmd(1, "/x/SU2_CFD", "a.cfg") == ["/x/SU2_CFD", "a.cfg"]


@pytest.mark.skipif(not os.path.isfile(GEOM_JSON), reason="нет сетки core130 coarse узла geom")
def test_geom_json_contract():
    """Стыковка с узлом geom: ось z, угол 360/Z, маркеры — те, что ждёт шаблон."""
    mj = json.load(open(GEOM_JSON, encoding="utf-8"))
    s2 = mj["su2"]
    assert s2["periodic"]["axis"] == [0, 0, 1]
    assert 360.0 / s2["periodic"]["angle_deg"] == pytest.approx(round(360.0 / s2["periodic"]["angle_deg"]))
    assert set(s2["markers"]) == {"inlet", "outlet", "hub", "shroud", "blade", "periodic_1", "periodic_2"}
    assert s2["rotation_rate_rad_s"][2] == pytest.approx(mj["d1"]["rpm"] * 2 * math.pi / 60, rel=1e-9)


SU2_FIXTURE = """NDIME= 3
NELEM= 1
10 0 1 2 3 0
NPOIN= 6
0.10 0.0 0.0 0
0.30 0.0 0.0 1
0.0 0.30 0.0 2
0.0 0.0 0.1 3
0.25 0.05 0.0 4
0.0 0.12 0.0 5
NMARK= 3
MARKER_TAG= hub
MARKER_ELEMS= 2
5 0 1 2
5 1 2 4
MARKER_TAG= shroud
MARKER_ELEMS= 1
5 0 3 5
MARKER_TAG= inlet
MARKER_ELEMS= 1
5 0 1 3
"""


def test_split_vaneless(tmp_path):
    """Придуманная сетка (FIXTURE): грань hub с узлами только на r >= 0,2 уходит в hub_vaneless."""
    src, dst = tmp_path / "m.su2", tmp_path / "s.su2"
    src.write_text(SU2_FIXTURE, encoding="utf-8")
    counts = rc.split_vaneless(str(src), str(dst), 0.2)
    assert counts == {"hub_vaneless": 1, "shroud_vaneless": 0}
    h = rc.read_mesh_header(str(dst))
    assert h["NMARK"] == 5
    assert h["markers"] == {"hub": 1, "hub_vaneless": 1, "shroud": 1, "shroud_vaneless": 0, "inlet": 1}
    txt = dst.read_text(encoding="utf-8")
    assert "MARKER_TAG= hub_vaneless\nMARKER_ELEMS= 1\n5 1 2 4\n" in txt
    assert txt.split("NMARK")[0] == SU2_FIXTURE.split("NMARK")[0]  # элементы и узлы не тронуты
