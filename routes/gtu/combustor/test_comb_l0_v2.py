"""Тесты уровня 0 v2 камеры KV (решение владельца 15.11.2026): v1 не тронут, границы v2, гидравлика одного клапана,
согласованное смешение, пассивные элементы, новая цель, части и слияние.

Запуск: python -m pytest -q routes/gtu/combustor/test_comb_l0_v2.py   (1-3 мин на 4 ядрах)
"""
import copy, json, math, os, sys

import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import comb_space as CS  # noqa: E402
import comb_crn as CR    # noqa: E402
import comb_v2 as V      # noqa: E402
import comb_l0_v2 as L   # noqa: E402
import merge_l0_v2 as M  # noqa: E402

N = 300
V1_CODE_HASH_AWS = "fae37b9aadf85961"     # прогон AWS 20261115T001745Z-f1da (l0_aws/, completeness.ref.code_hash)


@pytest.fixture(scope="module")
def bc():
    return V.boundary()


def _find(pred, n=N, need_ok=False):
    for i in range(n):
        c = V.concept(i, n)
        if not pred(c):
            continue
        if not need_ok:
            return i, c
        r = V.evaluate(i, n)
        if r["status"] == "ok":
            return i, r
    raise AssertionError("нет варианта")


def test_v1_code_untouched():
    """v2 не меняет comb_space.py и comb_crn.py: прогон AWS остаётся воспроизводимым из ветки"""
    assert CS.code_hash() == V1_CODE_HASH_AWS
    assert V.code_hash() != CS.code_hash()


def test_boundary_v2_adds_75_and_keeps_v1(bc):
    b1 = CR.boundary(write=False)
    loads = {L_["load"]: L_ for L_ in bc["loads"]}
    assert [L_["load"] for L_ in bc["loads"]] == V.LOADS
    for L1 in b1["loads"]:
        for k in ("p3_Pa", "T3_K", "W_fuel", "T41_C"):
            assert loads[L1["load"]][k] == pytest.approx(L1[k], rel=1e-9)
    L75 = loads[0.75]
    assert loads[0.5]["T41_C"] < L75["T41_C"] < loads[1.0]["T41_C"]
    assert loads[0.5]["W_fuel"] < L75["W_fuel"] < loads[1.0]["W_fuel"]
    assert L75["W_air"] == pytest.approx(loads[1.0]["W_air"])
    # точка 50 % пересчитана кодом v2 и совпала с v1 (самопроверка способа)
    assert all(v["rel"] < 1e-4 for v in bc["check"]["point_50_recomputed"].values())


def test_concept_repeatable_and_allowed():
    a = [V.concept(i, N) for i in range(N)]
    assert a == [V.concept(i, N) for i in range(N)]
    for c in a:
        s = c["scheme"]
        assert s in V.SCHEME_NAMES and s not in ("CAT", "RDC")
        for k in ("split", "strat", "recirc"):
            assert c[k] in V.ALLOWED[k][s]
        assert c["passive"] == V.passive_of(c)
        if c["split"] not in ("гидравлика", "каналы"):
            assert c["off"] == 0 and c["f_por"] == 0 and c["x_lin"] == 0
        if c["strat"] != "струи":
            assert c["k_str"] == 0
    assert V.fingerprint(N) != V.fingerprint(N + 1)


def test_porous_and_orifice_laws():
    """пористая вставка: при малом перепаде расход ~ ΔP (дарсиевский член), сопло ~ ΔP^0,5"""
    p_d = 1.3e6
    cir_por = dict(A=1e-5, a=1e9, b=0.0)
    m1, _ = V.circuit_m(cir_por, p_d * 1.001, p_d)
    m2, _ = V.circuit_m(cir_por, p_d * 1.002, p_d)
    o1 = V.orifice_m(1e-5, p_d * 1.001, p_d)
    o2 = V.orifice_m(1e-5, p_d * 1.002, p_d)
    assert o2 / o1 == pytest.approx(math.sqrt(2), rel=0.02)
    assert m2 / m1 > o2 / o1 + 0.2                      # ближе к линейному закону
    assert V.circuit_m(cir_por, p_d * 0.99, p_d)[0] == 0.0


def test_manifold_mass_balance_and_split_rises(bc):
    """один коллектор: сумма расходов = расходу клапана; доля пилота растёт при падении нагрузки, если у основного
    контура пористая вставка и/или впрыск при большем давлении; без них доля почти постоянна"""
    cands = [V.concept(i, N) for i in range(N)]
    cands = [c for c in cands if c["split"] == "гидравлика" and c["scheme"] != "TVC"]
    shares, c = {}, None
    for c0 in cands:
        ds = {nm: V.derive(dict(c0, p_sup_u=0.0, **kw), bc) for nm, kw in
              (("plain", dict(off=0.0, x_lin=0.0, f_por=0.2)), ("law", dict(off=1.0, x_lin=0.95, f_por=0.9)))}
        if all(d["loads"] for d in ds.values()):
            c = c0
            break
    assert c is not None, "нет варианта с окном коллектора"
    for name, d in ds.items():
        for li, st in enumerate(d["loads"]):
            tot = sum(v["m"] for v in st["circ"].values())
            assert tot == pytest.approx(bc["loads"][li]["W_fuel"], rel=1e-6)
        shares[name] = [st["circ"]["pilot"]["share"] for st in d["loads"]]
    assert shares["law"][3] > shares["law"][0] * 1.3
    assert abs(shares["plain"][3] - shares["plain"][0]) < 0.25 * shares["plain"][0]


def test_mixing_consistent_with_J():
    """меньше J (частичная нагрузка) - больше неравномерность и медленнее затухание; BETA = 0 - затухание как v1"""
    J_opt = (CS.C_OPT / (2 * math.pi / 4)) ** 2                  # C = C_OPT при 4 струях на элемент
    hi = V.mixing(J_opt, 4, 5.0, 0.0, 0.0)
    lo = V.mixing(0.2 * J_opt, 4, 5.0, 0.0, 0.0)                 # расход топлива ниже в 2,2 раза
    assert lo["C"] < hi["C"] and lo["S"] > hi["S"]
    assert hi["K_decay"] == pytest.approx(CS.K_MIX) and lo["K_decay"] > CS.K_MIX
    assert V.mixing(0.2 * J_opt, 4, 5.0, 0.0, 0.0, beta=0.0)["K_decay"] == pytest.approx(CS.K_MIX)


def test_stratification_law():
    a, a0 = V.strat_frac(CS.C_OPT * 1.2, 1.0, 3.0, CS.K_MIX)
    assert a == pytest.approx(1.0) and a0 == 1.0               # струи доходят - расслоения нет
    a_s, a0_s = V.strat_frac(CS.C_OPT * 0.5, 1.0, 3.0, 2 * CS.K_MIX)
    a_l, _ = V.strat_frac(CS.C_OPT * 0.5, 1.0, 12.0, 2 * CS.K_MIX)
    assert a0_s == pytest.approx(0.5) and a0_s < a_s < a_l < 1.0  # затухает по длине смешения


def test_backflow_rule(bc):
    """большой сдвиг давления при низком давлении коллектора - на малой нагрузке основной контур останавливается:
    правило ОБРАТНЫЙ_ТОК (воздух в топливный контур)"""
    i, c = _find(lambda c: c["split"] == "каналы")
    cc = dict(c, off=1.0, x_lin=0.0, f_por=0.9, p_sup_u=0.0, dp_inj=0.06)
    d = V.derive(cc, bc)
    if not d["loads"]:
        pytest.skip("нет окна коллектора")
    bm = min(s["backflow_margin"] for s in d["loads"] + [d["idle"]])     # и на холостом ходу (0 %, оценка)
    codes = [b[0] for b in V.check(cc, d, bc)]
    assert (bm < V.DP_BACK_MIN) == ("ОБРАТНЫЙ_ТОК" in codes)
    assert d["idle"]["backflow_margin"] <= min(s["backflow_margin"] for s in d["loads"]) + 1e-9


def test_evaluate_four_loads_targets_and_repeatable(bc):
    i, r = _find(lambda c: c["scheme"] in ("LPP", "MM", "TVC"), need_ok=True)
    assert [g["load"] for g in r["gas"]] == V.LOADS
    assert V.NEW_TARGET in r["targets"] and len(r["targets"]) == 14
    t = r["targets"][V.NEW_TARGET]
    v, terms, stable = V.lowem_value(r)
    assert t["value"] == pytest.approx(v, rel=1e-4) and t["ok"] == (v <= 1.0)   # строка округлена до 6 знаков
    assert set(k.split("_")[0] for k in terms) >= {"NOx", "CO", "lbo", "eta", "valves", "actuators"}
    r2 = V.evaluate(i, N)
    a, b = copy.deepcopy(r), copy.deepcopy(r2)
    for x in (a, b):
        x.pop("t_s", None)
    assert json.dumps(a["targets"], sort_keys=True) == json.dumps(b["targets"], sort_keys=True)


def test_unburnt_load_fails_new_target(bc):
    """полнота < 0,99 на любой нагрузке - новая цель не выполнена, даже если NOx и CO около нуля"""
    for i in range(N):
        r = V.evaluate(i, N)
        if r["status"] == "ok" and min(g["eta_c"] for g in r["gas"]) < CR.ETA_MIN:
            assert not r["targets"][V.NEW_TARGET]["ok"]
            return
    pytest.skip("нет варианта с неполным сгоранием")


def test_recirculation_ejector_law(bc):
    """эжектор: кратность по Ricou - Spalding растёт с длиной трубы и не больше KV_MAX"""
    i, c = _find(lambda c: c["recirc"] == "рециркуляция")
    Ks = []
    for l_ej in (4.0, 8.0, 15.0):
        cc = dict(c, l_ej=l_ej)
        d = V.derive(cc, bc)
        if not d["loads"]:
            pytest.skip("нет окна коллектора")
        Ks.append(V.gas_network(cc, d, bc, 0)["K_ej"])
    assert Ks[0] <= Ks[1] <= Ks[2] <= CR.KV_MAX
    assert Ks[2] > 0


def test_targets_v2_file():
    t = V.targets_v2()
    base = [x for x in t["base_targets"] if x != "fuel_supply_compatibility"]   # запись без значения («не задано»)
    assert len(base) == 13
    assert t["funnel_interface"]["hard"][-1] == V.NEW_TARGET
    assert t["new_targets"][0]["target"]["control_valves_per_fuel"] == 1


def test_shards_merge_complete_and_detects_loss(tmp_path, bc):
    n, S = 12, 3
    out = tmp_path / "parts"
    for k in range(S):
        L.run("KV-3", n, k, S, str(out), workers=1, quiet=True)
    complete, res = M.merge([str(out)], str(tmp_path / "m.json"), top=5)
    assert complete and res["counts"]["rows"] == n
    assert res["convergence"][-1]["n"] == n
    os.remove(out / L.shard_name("KV-3", n, 1, S))
    complete2, res2 = M.merge([str(out)], str(tmp_path / "m2.json"))
    assert not complete2 and res2["completeness"]["missing_shards"] == [1]


def test_merge_rejects_mixed_beta(tmp_path, bc):
    n, S = 6, 2
    out = tmp_path / "p"
    L.run("KV-3", n, 0, S, str(out), workers=1, quiet=True)
    L.run("KV-3", n, 1, S, str(out), workers=1, quiet=True, beta=0.0)
    V.BETA_LOWJ = V.BETA_LOWJ_NOMINAL
    complete, res = M.merge([str(out)], str(tmp_path / "m.json"))
    assert not complete and any("beta" in p for p in res["completeness"]["problems"])


def test_hydraulic_pilot_shares_main_residence_time():
    """пилот «гидравлика» - часть смесителей того же купола: время пилотной зоны = tau_pz; своё время - только у полости TVC"""
    for i in range(N):
        c = V.concept(i, N)
        if c["split"] == "гидравлика" and c["scheme"] != "TVC":
            assert c["tau_aux_ms"] == c["tau_pz_ms"]
        if c["scheme"] == "TVC":
            assert V.TAU_AUX[0] <= c["tau_aux_ms"] <= V.TAU_AUX[1]


def test_ejector_jetpump_bounded():
    """эжектор в трубе: кратность растёт с площадью подсоса и при реальной трубе (a <= 4) меньше закона свободной струи"""
    rs, l_ej = 683.5 / 1700.0, 15.0
    ks = [V.ejector_K(a, rs, l_ej, 0.03) for a in (0.5, 1.0, 2.0, 4.0, 8.0)]
    assert all(x < y for x, y in zip(ks, ks[1:]))
    assert ks[3] < 0.32 * l_ej * math.sqrt(rs) - 1.0
    assert V.ejector_K(0.0, rs, l_ej, 0.03) == 0.0


def test_fuel_stiffness_uses_flow_law(bc):
    """мягкий пористый контур (Дарси) откликается на пульсацию сильнее сопла: перепад «по отклику» не больше среднего"""
    for i in range(N):
        c = V.concept(i, N)
        if c["split"] in ("гидравлика", "каналы") and c["x_lin"] > 0.5:
            d = V.derive(c, bc)
            if d["loads"]:
                st = d["loads"][-1]
                assert st["dpf_eff"] <= st["dpf_weighted"] * 1.0001 or max(v["dpf"] for v in st["circ"].values()) > 0.5
                return
    pytest.skip("нет варианта")


def test_evaluate_order_independent(bc):
    """итог варианта не зависит от того, что процесс считал до него (пул раздаёт варианты в разном порядке)"""
    strip = lambda r: json.dumps({k: v for k, v in r.items() if k != "t_s"}, sort_keys=True)
    a2, a59 = V.evaluate(2, 300), V.evaluate(59, 300)
    b59, b2 = V.evaluate(59, 300), V.evaluate(2, 300)
    assert strip(a2) == strip(b2) and strip(a59) == strip(b59)


def test_idle_point_extrapolated(bc):
    """холостой ход (0 %): линейная экстраполяция по 50 и 30 % (оценка), квадратичная по четырём точкам - в пределах 2 %;
    расход топлива ~ четверть полного, воздух - как при 100 % (постоянная частота)"""
    I = bc["idle"]
    assert I["load"] == V.IDLE_LOAD == 0.0
    assert 0.2 < I["fuel_share_of_100"] < 0.3
    assert I["W_air"] == pytest.approx(bc["loads"][0]["W_air"])
    assert I["p3_Pa"] < bc["loads"][-1]["p3_Pa"] and I["T3_K"] < bc["loads"][-1]["T3_K"]
    assert all(v["rel"] < 0.02 for v in I["check_quadratic"].values())
    assert len(V.with_idle(bc)["loads"]) == V.IDLE_LI + 1 and len(bc["loads"]) == len(V.LOADS)


def test_zone_taus_constant_volume(bc):
    """времена зон при неизменном объёме: при 100 % - ровно оси, на 30 % у богатой зоны RQL короче (горячее и ниже p3)"""
    i, r = _find(lambda c: c["scheme"] == "RQL", need_ok=True)
    c = r["concept"]
    d = V.derive(c, bc)
    T100 = V.zone_taus(c, d, bc, 0, d["loads"][0])
    assert T100["pz"] == pytest.approx(1e-3 * c["tau_pz_ms"], rel=1e-12)
    assert T100["bo"] == pytest.approx(1e-3 * c["tau_bo_ms"], rel=1e-12)
    T30 = V.zone_taus(c, d, bc, 3, d["loads"][3])
    assert T30["pz"] < T100["pz"]
    assert r["gas"][3]["tau_ms"]["pz"] == pytest.approx(1e3 * T30["pz"], rel=1e-4)


def test_node_streams_conserve_fuel():
    """узлы неравномерности v2 не добавляют топливо при срезе бедного узла (v1 при S = 0,6 добавлял ~0,8 %)"""
    K = CR.kin()
    for S in (0.1, 0.4, 0.6):
        nodes, _ = V._node_streams_v2(K, 1.0, 0.02, S, 5, 683.5, 288.15, 1.37e6)
        fuel = sum(st.m - 1.0 * wk for wk, st in nodes)
        assert fuel == pytest.approx(0.02, rel=1e-9)
    nodes1, _ = CR._node_streams(K, 1.0, 0.02, 0.6, 5, 683.5, 288.15, 1.37e6)
    assert sum(st.m - wk for wk, st in nodes1) > 0.02 * 1.005


def test_tau_ext_refined_below_one():
    """запас ниже 1 уточняется между ступенями лестницы: не меньше ступени v1 и меньше предыдущей ступени"""
    K = CR.kin()
    fr = K.fresh(1.0, 0.5 * 0.0581, 683.5, 288.15, 1.37e6)        # phi 0,5 при p3, T3 точки 100 %
    got = {}
    for tau in (1e-4, 1e-5):                                     # v1: 1/3 и 1/30 (замер 15.11.2026)
        m1 = K.tau_ext(fr, tau)
        assert 0.0 < m1 < 1.0
        m2 = V.tau_ext_v2(K, fr, tau)
        f1 = 1.0 / m1
        prev = max([1.0] + [f for f in V.TAU_EXT_LADDER if f < f1 - 1e-9])
        assert m1 - 1e-12 <= m2 < 1.0 / prev
        got[tau] = m2
    # время погасания одно: запас пропорционален времени пребывания (шаг бисекции <= 1 %)
    assert got[1e-4] == pytest.approx(10 * got[1e-5], rel=0.02)


def test_cross_ignition_and_domain(bc):
    """каналы: верхняя оценка с поджиганием группы B газом группы A пишется рядом (gas_cross), отбор - без неё;
    пометка области облика - у каждой строки с гидравликой"""
    for i in range(N):
        c = V.concept(i, N)
        if c["split"] != "каналы":
            continue
        r = V.evaluate(i, N)
        if r["status"] != "ok":
            continue
        assert len(r["gas_cross"]) == len(V.LOADS)
        t = r["targets"][V.NEW_TARGET]
        assert "value_cross" in t and "in_domain" in t and isinstance(r["domain"]["flags"], list)
        assert r["targets"]["load_range"]["idle_margin_lbo"] == r["idle"]["margin_lbo"]
        return
    pytest.skip("нет допустимого варианта с каналами")
