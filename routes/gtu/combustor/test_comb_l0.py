"""Тесты уровня 0 воронки камеры KV: повторяемость, части без потерь, слияние, монотонности моделей.

Запуск: python -m pytest -q routes/gtu/combustor/test_comb_l0.py   (около 1-3 мин на 4 ядрах; границы - из
comb_bc_results.json, при несовпадении ключа считаются заново, 5-8 мин)
"""
import copy, json, os, sys

import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import comb_space as CS  # noqa: E402
import comb_crn as CR    # noqa: E402
import comb_l0 as L0     # noqa: E402
import merge_l0 as M     # noqa: E402

N = 300


@pytest.fixture(scope="module")
def bc():
    return CR.boundary()


def _first_ok(scheme=None, start=0):
    for i in range(start, N):
        c = CS.concept(i, N)
        if scheme and c["scheme"] != scheme:
            continue
        r = CR.evaluate(i, N)
        if r["status"] == "ok":
            return i, r
    raise AssertionError(f"нет допустимого варианта {scheme}")


def _strip(r):
    r = copy.deepcopy(r)
    r.pop("t_s", None)
    return r


def test_concept_repeatable_and_in_range():
    a = [CS.concept(i, N) for i in range(N)]
    assert a == [CS.concept(i, N) for i in range(N)]
    for c in a:
        allowed = {k: CS.SCHEMES[c["scheme"]].get(k, CS._ALL[k]) for k in CS._ALL}
        for k in CS._ALL:
            assert c[k] in allowed[k], (c["idx"], k, c[k])
        lo, hi = CS.SCH[c["scheme"]]["a_pz"]
        assert lo - 1e-9 <= c["a_pz"] <= hi + 1e-9
        assert 0.0 <= c["d_jet_u"] < 1.0


def test_fingerprint():
    assert CS.fingerprint(N) == CS.fingerprint(N)
    assert CS.fingerprint(N) != CS.fingerprint(N + 1)
    assert CS.fingerprint(N, seed=1) != CS.fingerprint(N)


@pytest.mark.parametrize("n,S", [(300, 4), (30000, 30), (17, 5), (10, 10), (11, 50)])
def test_shard_partition_lossless(n, S):
    parts = [CS.shard_indices(n, k, S) for k in range(S)]
    flat = [i for p in parts for i in p]
    assert sorted(flat) == list(range(n))
    assert len(flat) == len(set(flat))


def test_shard_bad_args():
    with pytest.raises(ValueError):
        CS.shard_indices(10, 3, 3)
    with pytest.raises(ValueError):
        CS.shard_indices(10, 0, 0)


def _close(a, b, rel=1e-4):
    """равенство с допуском: решатели Cantera (equilibrate, steady state) стартуют с текущего состояния объекта, поэтому
    первый вызов в новом процессе и повторный расходятся в 6-й значащей цифре (замечено проверкой скептика 14.11.2026)"""
    if isinstance(a, dict):
        return a.keys() == b.keys() and all(_close(a[k], b[k], rel) for k in a)
    if isinstance(a, list):
        return len(a) == len(b) and all(_close(x, y, rel) for x, y in zip(a, b))
    if isinstance(a, float) and isinstance(b, (float, int)):
        return abs(a - b) <= rel * max(abs(a), abs(b), 1e-9)
    return a == b


def test_evaluate_repeatable(bc):
    i, r1 = _first_ok()
    r2 = CR.evaluate(i, N)
    assert _close(_strip(r1), _strip(r2))


def test_shards_merge_complete_and_detects_loss(tmp_path, bc):
    n, S = 12, 3
    out = tmp_path / "parts"
    for k in range(S):
        L0.run("KV-3", n, k, S, str(out), workers=1, quiet=True)
    complete, res = M.merge([str(out)], str(tmp_path / "m.json"), top=5)
    assert complete and res["counts"]["rows"] == n
    assert res["completeness"]["missing_indices"] == 0
    # потеря части - неполное слияние
    os.remove(out / L0.shard_name("KV-3", n, 1, S))
    complete2, res2 = M.merge([str(out)], str(tmp_path / "m2.json"))
    assert not complete2 and res2["completeness"]["missing_shards"] == [1]
    assert res2["completeness"]["missing_indices"] == len(CS.shard_indices(n, 1, S))


def test_merge_rejects_foreign_fingerprint(tmp_path, bc):
    n, S = 6, 2
    out = tmp_path / "p"
    for k in range(S):
        L0.run("KV-3", n, k, S, str(out), workers=1, quiet=True)
    f = out / L0.shard_name("KV-3", n, 0, S)
    d = json.load(open(f, encoding="utf-8"))
    d["fingerprint"] = "0" * 16
    json.dump(d, open(f, "w", encoding="utf-8"))
    complete, res = M.merge([str(out)], str(tmp_path / "m.json"))
    assert not complete and any("fingerprint" in p for p in res["completeness"]["problems"])


def test_co_falls_with_burnout_time(bc):
    """CO на выходе цепочки PSR -> PFR(t) убывает с временем догорания (phi 0,5, вход KV-3 100 %)"""
    K = CR.kin()
    L = bc["loads"][0]
    fr = K.fresh(1.0, 0.5 * bc["far_st_CH4"], L["T3_K"], 288.15, L["p3_Pa"])
    o, lit, _ = K.psr(fr, 0.5e-3)
    assert lit
    co = [K.emissions(K.pfr(o, t), fr.Y[K.i["CH4"]], bc["LHV_J_kg"], bc["LHV_J_kg"]["CH4"])["CO"]
          for t in (0.2e-3, 1e-3, 3e-3, 10e-3)]
    assert all(a > b for a, b in zip(co, co[1:])), co


def test_nox_rises_with_flame_temperature(bc):
    """NOx растёт с температурой зоны пламени (phi 0,45 -> 0,65 при том же времени)"""
    K = CR.kin()
    L = bc["loads"][0]
    nox, T = [], []
    for phi in (0.45, 0.5, 0.55, 0.6, 0.65):
        fr = K.fresh(1.0, phi * bc["far_st_CH4"], L["T3_K"], 288.15, L["p3_Pa"])
        o, lit, _ = K.psr(fr, 1e-3)
        o = K.pfr(o, 10e-3)
        e = K.emissions(o, fr.Y[K.i["CH4"]], bc["LHV_J_kg"], bc["LHV_J_kg"]["CH4"])
        nox.append(e["NOx"]); T.append(e["T"])
    assert all(a < b for a, b in zip(T, T[1:]))
    assert all(a < b for a, b in zip(nox, nox[1:])), nox


def test_network_monotonic_in_tau_and_temperature(bc):
    """на сети схемы LPP: дольше догорание - CO не больше; меньше воздуха в зоне (выше T) - NOx больше"""
    i, r = _first_ok("LPP")
    c = r["concept"]
    def run(**kw):
        cc = dict(c, **kw)
        d = CS.derive(cc, bc)
        return CR.gas_network(cc, d, bc, 0)
    g_short, g_long = run(tau_bo_ms=3.0), run(tau_bo_ms=15.0)
    assert g_long["CO"] <= g_short["CO"] * 1.01 + 0.05, (g_short["CO"], g_long["CO"])
    g_lean, g_rich = run(a_pz=min(0.85, c["a_pz"] + 0.08)), run(a_pz=max(0.35, c["a_pz"] - 0.08))
    assert g_rich["T_hot"] > g_lean["T_hot"]
    assert g_rich["NOx"] > g_lean["NOx"], (g_lean["NOx"], g_rich["NOx"])


def test_pilot_share_rises_at_low_flow():
    assert CR.pilot_share(0.2, 1.0) == pytest.approx(0.2)
    assert CR.pilot_share(0.2, 0.45) > CR.pilot_share(0.2, 0.7) > 0.2


def test_pattern_factor_falls_with_length(bc):
    i, r = _first_ok()
    c = r["concept"]
    d = CS.derive(c, bc)
    g = r["gas"][0]
    pf2 = CR.pattern_factor(dict(c, LD=2.0), d, bc, g)["PF"]
    pf4 = CR.pattern_factor(dict(c, LD=4.0), d, bc, g)["PF"]
    assert pf4 < pf2


def test_wall_hotter_with_hotter_gas(bc):
    i, r = _first_ok()
    c = r["concept"]
    d = CS.derive(c, bc)
    w1 = CR.wall_temperature(c, d, bc, 1700.0, 0.03)["T_w1"]
    w2 = CR.wall_temperature(c, d, bc, 1900.0, 0.03)["T_w1"]
    assert w2 > w1


# --------------------------------------------------------------------------- проверка скептика 14.11.2026: выгоды даром
def _first_concept(pred):
    for i in range(N):
        c = CS.concept(i, N)
        if pred(c):
            return c
    raise AssertionError("нет варианта")


def test_flox_recirculation_limited_by_confinement(bc):
    """плотная решётка струй не получает рециркуляцию; реже струи - больше уноса; не больше KV_MAX"""
    c = _first_concept(lambda c: c["scheme"] == "FLOX")
    d = CS.derive(c, bc)
    k_dense = CR.flox_recirc(dict(c, n_el=4000), d, bc, 0)
    k_sparse = CR.flox_recirc(dict(c, n_el=8), d, bc, 0)
    assert k_dense["K"] == 0.0 and k_dense["D_cell_over_d0"] < 3
    assert k_sparse["K"] >= k_dense["K"] and k_sparse["K"] <= CR.KV_MAX


def test_stratification_decays_with_mixing_length():
    c = dict(l_mix=15.0, k_pass=1.0)
    assert CR.strat_air_frac(c, 0.5) > 0.99                    # длинный канал - расслоение перемешано
    assert CR.strat_air_frac(dict(c, l_mix=3.0), 0.5) < 0.9    # короткий - остаётся
    assert CR.strat_air_frac(c, 1.0) == pytest.approx(1.0)


def test_hole_tolerance_adds_unmixedness(bc):
    """напечатанные отверстия (разброс 0,025 мм) дают больший разброс phi, чем сверлёные (0,008 мм)"""
    c = _first_concept(lambda c: c["injector"] == "струи")
    dp, dc = CS.derive(dict(c, make="печатный картридж"), bc), CS.derive(dict(c, make="литьё"), bc)
    assert dp["S_tol"] > 0 and dc["S_tol"] > 0
    assert dp["S_unmix"] >= dp["S_mix"] - 1e-12
    assert 2 * CS.D_SIG_FUEL["печатный картридж"] / dp["d_jet_mm"] > 2 * CS.D_SIG_FUEL["литьё"] / dc["d_jet_mm"]
    assert dp["dp_mix_fric"] > dc["dp_mix_fric"] > 0          # шероховатость печати - большее трение канала


def test_film_effectiveness_within_source_branch(bc):
    """η плёнки не выше ветви m > 1,3 (1,28·...) - без экстраполяции m^0,65"""
    i, r = _first_ok()
    c = r["concept"]
    if c["cooling"] == "решётка" or c["porosity"] <= 0:
        pytest.skip("нет плёнки")
    d = CS.derive(c, bc)
    w = CR.wall_temperature(c, d, bc, 1800.0, 0.03)
    xs = 0.5 / max(c["porosity"], 1e-4)
    assert w["eta_film"] <= 1.28 * 1.0 * xs ** -0.2 * 0.5 ** -0.2 + 1e-9


def test_rdc_gets_no_pressure_gain_and_friction_counted(bc):
    assert CR.ETA_PG == 0.0
    i, r = _first_ok()
    assert r["dp"]["p_gain"] == 0.0 or r["concept"]["scheme"] != "RDC"
    assert r["dp"]["dp_total"] >= r["dp"]["dp_inj"] + r["dp"]["dp_mix_fric"] - 1e-12


def test_unlit_load_does_not_pass_emission_targets(bc):
    """погасшая на 50 % зона (полнота < 0,99) - цели NOx/CO не выполнены, хотя NOx и CO около нуля"""
    for i in range(N):
        r = CR.evaluate(i, N)
        if r["status"] == "ok" and r["gas"][1]["eta_c"] < CR.ETA_MIN:
            assert not r["targets"]["nox_gas"]["ok"] and not r["targets"]["co_gas"]["ok"]
            assert r["feed"]["burnt"][1] is False
            return
    pytest.skip("нет варианта, погасшего на 50 %")
