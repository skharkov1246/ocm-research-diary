"""Балансы массы и энергии модели цикла (gtu_cycle.cycle, cycle_r, recup). Запуск: python3 -m pytest -q
из папки routes/gtu/validation (около 30 с)."""
import math

import pytest

import balances as B

G = B.G


@pytest.fixture(scope="module")
def results():
    out = {}
    for p in B.points():
        kw = B.resolve_kw(p)
        r = B.run_point(p["kind"], kw)
        out[p["id"]] = (p, kw, r, B.balance(p["kind"], kw, r))
    return out


def test_points_cover_all_functions(results):
    kinds = {p["kind"] for p, *_ in results.values()}
    assert kinds == {"cycle", "cycle_r", "recup"}
    assert len(results) >= 6


def test_mass_balance(results):
    for pid, (_, _, _, b) in results.items():
        assert abs(b["mass_residual_pct"]) < 1e-9, pid


def test_energy_balance_absolute_enthalpy(results):
    # граница — весь двигатель; невязка — доля теплоты топлива
    for pid, (_, _, _, b) in results.items():
        assert abs(b["energy_residual_abs_pct"]) < 1e-3, (pid, b["energy_residual_abs_pct"])


def test_energy_balance_lhv_form_matches_absolute(results):
    # вторая проверка: форма «теплота топлива = работа + выхлоп + потери» отличается от абсолютной
    # ровно на несогласованность LHV модели (50,03 МДж/кг) с Cantera при T0
    for pid, (_, _, _, b) in results.items():
        diff = b["energy_residual_lhv_pct"] - b["energy_residual_abs_pct"]
        assert abs(diff - b["lhv_form_expected_pct"]) < 1e-6, pid
        assert abs(b["energy_residual_lhv_pct"]) < 0.05, pid


def test_lhv_constant_consistent_with_cantera():
    assert abs(B.lhv_cantera(G.T0) / G.LHV - 1) < 5e-4


def test_terms_add_up_and_are_positive(results):
    for pid, (_, kw, r, b) in results.items():
        for k in ("Q_fuel_MW", "Q_exh_MW", "L_mech_MW", "L_comb_MW", "P_MW"):
            assert b[k] > 0, (pid, k)
        # КПД модели = P / (Wf·LHV)
        assert math.isclose(b["P_MW"] / b["Q_fuel_MW"], r["eff"], rel_tol=1e-9), pid


def test_exhaust_composition_reconstruction():
    g = B._gas()
    Y = B.exhaust_Y(20.0, 0.35, 0.015)
    assert abs(Y.sum() - 1) < 1e-12
    assert Y[g.species_index("CH4")] == 0
    assert 0.12 < Y[g.species_index("O2")] < 0.20       # бедная смесь ГТУ
    assert Y[g.species_index("CO2")] > 0 and Y[g.species_index("H2O")] > 0


def test_balance_detects_injected_error(results):
    # контроль чувствительности: лишний 1 % мощности обязан появиться в невязке
    p, kw, r, b = results["cycle:T60-калибровка"]
    r2 = dict(r, P_MW=r["P_MW"] * 1.01)
    b2 = B.balance(p["kind"], kw, r2)
    expected = -0.01 * r["P_MW"] * (1 + (1 / kw.get("eta_m", 0.99) - 1)) / b["Q_fuel_MW"] * 100
    assert math.isclose(b2["energy_residual_abs_pct"] - b["energy_residual_abs_pct"], expected, rel_tol=1e-3)
    # и лишний 1 % расхода выхлопа — в балансе массы
    r3 = dict(r, W_exh=r["W_exh"] * 1.01)
    assert abs(B.balance(p["kind"], kw, r3)["mass_residual_pct"]) > 0.9


def test_recup_stack_temperature(results):
    # точная температура за рекуператором из баланса; оценка recup (0,93) — с известной ошибкой
    for pid, (p, kw, r, b) in results.items():
        if p["kind"] != "recup":
            continue
        assert math.isclose(B.stack_temperature(r, kw), b["Tstack_exact_C"], abs_tol=1e-6)
        assert b["Tstack_exact_C"] < b["EGT_C"]
        assert b["T3r_C"] > b["T3_C"]
        # ошибка оценки фиксируется в результатах; невязка энергии с оценкой — меньше 5 %
        assert abs(b["energy_residual_recup_estimate_pct"]) < 5.0, pid
        assert abs(b["Tstack_estimate_error_K"]) > 1.0, pid   # оценка не точна — это находка проверки


def test_cycle_scales_linearly_with_airflow():
    # основа проверки вне калибровки: удельные величины и КПД не зависят от W1
    kw = dict(PR=12.2, T41_C=1050.0, eta_pc=0.874, eta_pt_gg=0.864, eta_pt_pt=0.884, x_ngv=0.045, x_rot=0.045,
              x_leak=0.015)
    r1, r20 = G.cycle(W1=1.0, **kw), G.cycle(W1=20.0, **kw)
    assert math.isclose(r20["P_MW"], 20 * r1["P_MW"], rel_tol=1e-6)
    assert math.isclose(r20["W_exh"], 20 * r1["W_exh"], rel_tol=1e-9)
    assert math.isclose(r20["eff"], r1["eff"], rel_tol=1e-6)
    assert math.isclose(r20["EGT_C"], r1["EGT_C"], abs_tol=1e-3)
