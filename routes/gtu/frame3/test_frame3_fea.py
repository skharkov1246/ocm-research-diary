"""Проверки прочности колёс и дисков типоразмера № 3 (frame3_fea.py). Запуск: python3 -m pytest -q test_frame3_fea.py

Проверяется frame3_fea_results.json: поля происхождения, совпадение входа с текущей расчётной точкой (иначе расчёт
устарел), U = ω·r, ворота равновесия, сверка с аналитическим диском, масштаб запаса по разрушению с частотой,
сходимость по сеткам, согласованность оценок ползучести; геометрия профилей и формула напряжения в корне лопатки —
без ccx; повтор двух расчётов ccx (колесо ступени 1 и диск турбины ступени 1, грубая сетка) — около 10 с."""
import json
import math
import os
import sys
import tempfile

import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import frame3_fea as M  # noqa: E402
import disk  # noqa: E402  (frame3_fea добавляет ../platform/fea в путь)

RES = json.load(open(os.path.join(HERE, "frame3_fea_results.json"), encoding="utf-8"))
POINT = M.load_point()
N0 = RES["speeds"]["n_design_rpm"]


def all_cases():
    for part, e in RES["parts"].items():
        for var, v in e["variants"].items():
            for sp, c in v["cases"].items():
                yield part, var, sp, c


def test_provenance_fields():
    for k in ("method", "software", "date", "inputs", "status", "assumptions", "conclusions", "limitations"):
        assert RES.get(k), k
    assert RES["date"] == "2026-11-14"
    assert "2.21" in RES["software"]
    for i in RES["inputs"]:
        assert i["file"] and i["field"]
    text = json.dumps(RES, ensure_ascii=False)
    assert "/tmp/" not in text and "/home/" not in text      # рабочие папки не попадают в репозиторий


def test_inputs_match_current_point():
    """если расчётная точка пересчитана, а прочность — нет, тест падает: FEA устарел"""
    tin = M.inputs_from_point(POINT)
    used = {i["field"].split(":")[0]: i["used"] for i in RES["inputs"]}
    assert used["stage1"] == tin["stage1"]
    assert used["stage2"] == tin["stage2"]
    t = next(i["used"] for i in RES["inputs"] if i["field"].startswith("turbine"))
    assert t == tin["turbine"]


def test_speeds_and_tip_speed():
    assert RES["speeds"]["factors"] == [1.0, 1.15, 1.2]
    for part, var, sp, c in all_cases():
        assert abs(c["rpm"] - N0 * int(sp) / 100) < 1e-6 * N0
        om = c["rpm"] * math.pi / 30
        assert abs(c["U_tip"] - om * c["r_out_mm"] / 1e3) < 1e-6 * c["U_tip"]
    s1 = RES["parts"]["comp1"]["variants"]["solid"]["cases"]["100"]
    assert abs(s1["U_tip"] - POINT["stage1"]["U2_m_s"]) < 1e-3 * POINT["stage1"]["U2_m_s"]


def test_equilibrium_gates_passed():
    assert RES["checks"]["all_gates_passed"]
    assert RES["checks"]["equilibrium_gate_max"] < 1e-3
    for part, var, sp, c in all_cases():
        assert c["quality"]["passed"], (part, var, sp)


def test_analytic_uniform_disk():
    """тонкий сплошной диск постоянной толщины того же радиуса, материала и частоты: σθ(0), σθ(R) (как test_disk.py)"""
    an = RES["checks"]["analytic_uniform_disk"]
    assert set(an) == {"comp1", "comp2", "turb1", "turb2", "turb3"}
    for part, v in an.items():
        assert abs(v["dev"]["hoop_centre"]) < 3e-3, part
        assert abs(v["dev"]["hoop_rim"]) < 3e-3, part
        if v["rim_p_MPa"]:
            assert abs(v["fea"]["radial_rim_MPa"] / v["rim_p_MPa"] - 1) < 3e-3


def test_burst_margin_scales_with_speed():
    """∫σθ dA от нагрева равен нулю (свободный диск), поэтому запас по разрушению обратно пропорционален частоте"""
    n = 0
    for part, e in RES["parts"].items():
        for var, v in e["variants"].items():
            c = v["cases"]
            if "100" in c and "120" in c:
                assert abs(c["100"]["burst_margin"] / c["120"]["burst_margin"] - 1.2) < 2e-3, (part, var)
                assert abs(c["100"]["burst_rpm"] - c["120"]["burst_rpm"]) < 2e-3 * c["100"]["burst_rpm"]
                n += 1
    assert n >= 7


def test_mesh_convergence():
    for part, var, sp, c in all_cases():
        cv = c["convergence_fine_vs_base"]
        for k in ("vm_max_node", "hoop_mean", "burst_margin", "centre_hoop"):
            assert abs(cv[k]) < 0.02, (part, var, sp, k, cv[k])


def test_compressor_bore_worse_than_solid():
    for part in ("comp1", "comp2"):
        v = RES["parts"][part]["variants"]
        for sp in ("100", "120"):
            assert v["bore"]["cases"][sp]["vm_max_node_MPa"] > v["solid"]["cases"][sp]["vm_max_node_MPa"]
            assert v["bore"]["cases"][sp]["yield_margin"] < v["solid"]["cases"][sp]["yield_margin"]


def test_creep_ratio_consistent_with_life():
    """отношение σ_LM/σ1 >= 1 тогда и только тогда, когда местный срок по наибольшему главному (disk.py) >= t"""
    n = 0
    for k in (1, 2):
        for tv, v in RES["parts"][f"turb{k}"]["variants"].items():
            c = v["cases"]["100"]
            for th in (25000, 50000):
                m = c["creep_margins"][f"{th}h"]["s1"]
                if m is None or abs(m["ratio"] - 1) < 0.02:
                    continue
                life = c["creep"]["local_s1_min_h"]
                assert (m["ratio"] >= 1) == (life is None or life >= th), (k, tv, th, m["ratio"], life)
                n += 1
    assert n >= 6


def test_titanium_stage2_verdict():
    t = RES["parts"]["comp2"]["titanium_temperature"]
    assert t["T02_C"] == POINT["stage2"]["T02_C"] and t["T02_C"] > 400
    assert t["T0_rel_rotor_exit_C"] < t["T02_C"]
    assert "Ti-6242" in t["proposal"]
    assert RES["parts"]["comp2"]["variants"]["solid"]["cases"]["100"]["creep"].startswith("нет данных")


def test_profiles_geometry():
    tin = M.inputs_from_point(POINT)
    s1, s2 = tin["stage1"], tin["stage2"]
    p = M.comp_profile(s1, double=True)
    t_rim = M.A["comp_t_rim_D2"] * s1["D2_m"]
    assert all(z == 0 for z in p["zc"])
    assert abs(p["r"][-1] - s1["D2_m"] / 2) < 1e-9
    assert abs(p["t"][-1] - t_rim) < 1e-9
    assert abs(p["hb"][-1] - 2 * s1["b2_m"]) < 1e-9
    assert abs(p["t"][0] - (t_rim + 2 * s1["axial_length_m"])) < 1e-9
    q = M.comp_profile(s2, double=False)
    assert abs(q["hb"][-1] - s2["b2_m"]) < 1e-9
    assert abs(q["t"][-1] - M.A["comp_t_rim_D2"] * s2["D2_m"]) < 1e-9
    assert abs(q["t"][0] - (s2["axial_length_m"] + M.A["comp2_boss_D2"] * s2["D2_m"])) < 1e-9
    for pr in (p, q):
        assert all(b > a for a, b in zip(pr["r"], pr["r"][1:]))
        assert min(pr["t"]) > 0
    ni = disk.load_material("in713c")
    om = N0 * math.pi / 30
    for k in range(3):
        g = M.turbine_stage_geometry(k, tin["turbine"], ni["density"])
        prof, info = M.turbine_profile(g, ni["density"], om, 300e6)
        assert max(prof["t"]) <= info["t_max_m"] + 1e-12 or info["band_t_m"] > info["t_max_m"]
        assert abs(prof["r"][-1] - g["r_live_rim_m"]) < 1e-12
        assert prof["t"][-1] == prof["t"][-2] == info["band_t_m"]
        assert g["r_live_rim_m"] < g["r_hub_m"] < g["r_tip_m"]
        assert abs(sum(v["m_kg"] for v in g["rim_mass_parts_kg"].values()) - g["rim_mass_kg"]) < 1e-9


def test_blade_root_formula_SYNTHETIC():
    """без сужения формула равна ρω²(r_t² − r_h²)/2 = ρ(π/30)²·AN²/(2π); значения — SYNTHETIC (выдуманы)"""
    SYNTHETIC = dict(rho=8000.0, rpm=12000.0, r_m=0.30, h=0.10)
    om = SYNTHETIC["rpm"] * math.pi / 30
    r_h, r_t = SYNTHETIC["r_m"] - SYNTHETIC["h"] / 2, SYNTHETIC["r_m"] + SYNTHETIC["h"] / 2
    s = M.blade_root_stress(SYNTHETIC["rho"], om, r_h, SYNTHETIC["h"], 1.0)
    assert abs(s / (SYNTHETIC["rho"] * om ** 2 * (r_t ** 2 - r_h ** 2) / 2) - 1) < 1e-12
    an2 = math.pi * 2 * SYNTHETIC["r_m"] * SYNTHETIC["h"] * SYNTHETIC["rpm"] ** 2
    assert abs(s / (SYNTHETIC["rho"] * (math.pi / 30) ** 2 * an2 / (2 * math.pi)) - 1) < 1e-12
    # сужение уменьшает напряжение; τ = 0 — треугольная площадь: ρω²·h·(r_h/2 + h/6)
    s0 = M.blade_root_stress(SYNTHETIC["rho"], om, r_h, SYNTHETIC["h"], 0.0)
    assert abs(s0 / (SYNTHETIC["rho"] * om ** 2 * SYNTHETIC["h"] * (r_h / 2 + SYNTHETIC["h"] / 6)) - 1) < 1e-12
    for k in (1, 2, 3):
        g = RES["parts"][f"turb{k}"]["geometry"]
        assert abs(g["sigma_root_untapered_exit_h_MPa"] / g["sigma_AN2_untapered_MPa"] - 1) < 1e-4
        assert abs(RES["parts"][f"turb{k}"]["blade"]["ccx_bar_check"]["dev"]) < 1e-3


@pytest.fixture(scope="module")
def rerun():
    """повтор двух расчётов на грубой сетке тем же построением, что в frame3_fea.build"""
    tin, meta, jobs = M.build(POINT)
    want = {("comp1", "solid", 1.0, "base"), ("turb1", "T_rel", 1.0, "base")}
    out = {}
    with tempfile.TemporaryDirectory() as d:
        for key, case in jobs:
            if key in want:
                out[key] = disk.run_case(case, workdir=os.path.join(d, "_".join(map(str, key))), keep=False)
    return out


def test_rerun_matches_results(rerun):
    for (part, var, sp, mk), r in rerun.items():
        c = RES["parts"][part]["variants"][var]["cases"]["100"]
        cv = c["convergence_fine_vs_base"]
        assert r["quality"]["passed"]
        base_vm = c["vm_max_node_MPa"] / (1 + cv["vm_max_node"])
        assert abs(r["vm_max_node_MPa"] / base_vm - 1) < 1e-6, (part, r["vm_max_node_MPa"], base_vm)
        base_b = c["burst_margin"] / (1 + cv["burst_margin"])
        assert abs(r["burst"]["margin"] / base_b - 1) < 1e-6
