"""Проверки библиотеки блоков «от топлива до вала» (pytest).

python3 -m pytest -q routes/gtu/process/blocks/test_blocks.py
Пороги сверки: та же модель (воронка) - до округления; другая реализация свойств газа (Cantera в gtu_cycle и переборе
циклов против таблиц gtu_explore) - 0,15 п. п. КПД, 0,6 % расхода, 2 K; frame3 - калибровка трека ±2 п. п.
"""
import os
import sys

for _v in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import pytest  # noqa: E402

import blocks  # noqa: E402
from blocks import _track as T  # noqa: E402
from blocks import chain, combustor, reference  # noqa: E402
from blocks import evaluator as EV  # noqa: E402
from blocks.concept import normalize  # noqa: E402

TOL_EFF_PP, TOL_G_PCT, TOL_T_K = 0.15, 0.6, 2.0


@pytest.fixture(scope="module")
def frame3():
    return reference.check_frame3()


def test_kv3_reproduces_funnel_exactly(frame3):
    """точка KV-3: evaluate = funnel_l0.solve той же концепции воронки (та же модель цикла, КПД и геометрия)"""
    ev, fu = frame3["evaluate"], frame3["funnel"]
    assert abs(ev["eff_shaft"] - fu["eff_shaft"]) < 0.005
    assert abs(ev["G_air"] / fu["G_air"] - 1) < 1e-4
    assert abs(ev["EGT_C"] - fu["EGT_C"]) < 0.05


def test_kv3_within_track_calibration(frame3):
    """точка KV-3 против frame3: КПД на валу в пределах калибровки трека ±2 п. п., удельная работа ±3 %, выхлоп ±10 K"""
    assert abs(frame3["d_eff_vs_frame3_pp"]) <= frame3["calibration_pp"]
    assert abs(frame3["d_w_vs_frame3_pct"]) <= 3.0
    assert abs(frame3["d_EGT_vs_frame3_K"]) <= 10.0


def test_kv3_eta_frame3_vs_gtu_cycle(frame3):
    """с КПД узлов frame3 цепочка повторяет gtu_cycle.cycle (Cantera) - второй способ"""
    assert abs(frame3["d_eff_given_vs_gtu_cycle_pp"]) <= TOL_EFF_PP
    assert abs(frame3["d_w_given_vs_gtu_cycle_pct"]) <= TOL_G_PCT
    assert abs(frame3["eta_given"]["EGT_C"] - frame3["gtu_cycle"]["EGT_C"]) <= TOL_T_K


@pytest.mark.parametrize("key", reference.CYCLES_KEYS)
def test_five_topologies_vs_cycles(key):
    """5 топологий перебора циклов (простой, рекуператор, рекуператор + подогрев, рекуператор + промохлаждение, подогрев
    топлива): КПД, расход воздуха, выхлоп и газ на входе рекуператора - как cycles_results.json"""
    row = reference.check_cycles([key])[0]
    assert abs(row["d_eff_pp"]) <= TOL_EFF_PP, row
    assert abs(row["d_W_air_pct"]) <= TOL_G_PCT, row
    assert abs(row["d_EGT_K"]) <= TOL_T_K, row
    if row["cycles"]["T_rec_in_C"] is not None:
        assert abs(row["blocks"]["T_rec_in_C"] - row["cycles"]["T_rec_in_C"]) <= TOL_T_K


@pytest.mark.parametrize("key", reference.EXTRA_KEYS)
def test_extra_topologies_vs_cycles(key):
    """IBC и рост давления - в цепочке; паровой риформинг - приращение из перебора (проверка подключения)"""
    row = reference.check_cycles([key])[0]
    assert abs(row["d_eff_pp"]) <= TOL_EFF_PP, row
    assert abs(row["d_W_air_pct"]) <= TOL_G_PCT, row


def test_recup_vs_gtu_cycle_recup():
    r = reference.check_recup_gtu_cycle()
    assert abs(r["d_eff_pp"]) <= TOL_EFF_PP
    assert abs(r["d_T3r_K"]) <= TOL_T_K


def test_two_shaft_recup_equals_funnel():
    """двухвальный с рекуператором (температурная запись, как воронка): цепочка = funnel_space quick_perf"""
    C = normalize(dict(compressor=dict(arch="1Ц", pi=5.0), turbine=dict(type="Р", cool="U", mat="K418", make="литьё", t41=950),
                       shaft=dict(drive="свободная СТ"), hx=dict(recup=0.85), combustor=dict(dp_source="track"),
                       options=dict(rec_model="temperature")))
    tgt = T.package_target("KV-3")
    fc, sol, etas, nt, npt, _ = EV.funnel_pass(C, tgt)
    p = sol["perf"]
    ch = chain.run(EV.chain_spec(C, etas, combustor.PRb(C, 950.0)))
    assert abs(ch["eff"] - p["eta_th"]) * 100 < 0.02
    assert abs(ch["T_stack"] - 273.15 - p["EGT"]) < 0.5


def test_energy_balance_closes():
    """первый закон по цепочке (промохлаждение + рекуператор + подогрев + подогрев топлива): вход = выход до 0,05 %"""
    E = T.E()
    C = normalize(dict(compressor=dict(arch="2Ц", pi=8.0, intercool="воздух"), hx=dict(recup=0.85, reheat=True, x_rh=0.4),
                       fuel=dict(prep="подогрев"), combustor=dict(dp_source="track"), turbine=dict(nt=3)))
    spec = EV.chain_spec(C, dict(comp=[0.85, 0.84], gg=0.88, pt=0.88), combustor.PRb(C, 1000.0))
    ch = chain.run(spec)
    assert ch["ok"]
    fuel_burned = ch["fuel"] * spec.get("eta_b", 0.995)
    e_in = E.h(E.T0) + fuel_burned * E.HF
    e_out = (ch["wnet"] + ch["W_fan"]) + (1 - 0.99) * ch["W_t"] + ch["Q_ic"] + ch["xl"] * E.h(ch["T2"]) + \
        ch["m_exh"] * E.h(ch["T_stack"], ch["b_exh"])
    assert abs(e_in - e_out) / (fuel_burned * E.LHV) < 5e-4


def test_compat_hard_rules():
    cases = {
        "ic_stages": dict(compressor=dict(arch="1Ц", pi=6, intercool="воздух")),
        "reheat_radial": dict(compressor=dict(arch="1Ц", pi=6), turbine=dict(type="Р", mat="K418", t41=950, nt=2), hx=dict(reheat=True)),
        "cool_C": dict(turbine=dict(type="О", cool="C", mat="SC")),
        "cool_P_print": dict(turbine=dict(cool="P", mat="SC", make="литьё", hot_make="литьё", t41=1050)),
        "comb_cat_preburner": dict(combustor=dict(scheme="CAT")),
    }
    for code, c in cases.items():
        r = blocks.evaluate(c)
        assert code in r["reject"], (code, r["reject"])
        assert not r["feasible"]


def test_soft_rules_are_warnings():
    r = blocks.evaluate(dict(fuel=dict(prep="паровой риформинг")))
    assert "pnp_water" in r["warnings"] and "pnp_water" not in r["reject"]
    assert r["plugplay"]["water"] and r["plugplay"]["pnp_cycles"] == "нет"
    assert r["water_t_h"] > 0


def test_highspeed_shaft_through_package_gearbox():
    r = blocks.evaluate(dict(shaft=dict(drive="быстроходный")))
    sh = r["blocks"]["shaft"]
    assert len(sh["n_scan"]) >= 10                       # частота выбрана перебором в окне отношения частот
    from blocks.shaft import HS_RATIO
    lo, hi = HS_RATIO
    assert lo - 1e-6 <= sh["n_engine_rpm"] / sh["n_package_rpm"] <= hi + 1e-6
    assert r["plugplay"]["states"]["вал и муфта"] == "переходник"
    assert any(it["trl"] == 7 for it in r["trl"]["items"])


def test_output_contract():
    """единый вызов: КПД, мощность, расход и выхлоп против допуска, детали, готовность, размеры для модели массы"""
    r = blocks.evaluate(dict(compressor=dict(pi=8.2), hx=dict(recup=0.90, reheat=True, x_rh=0.35)))
    for k in ("eff_shaft", "eff_el", "P_shaft_kW", "G_air", "W_exh", "EGT_C", "interface", "parts", "trl", "sizing", "blocks",
              "feasible", "rules"):
        assert k in r, k
    for b in ("fuel", "compressor", "combustor", "turbine", "hx", "shaft"):
        assert b in r["blocks"]
    assert r["sizing"]["hx"]["рекуператор"]["mass_kg"] > 0
    assert r["interface"]["G_cap"] == pytest.approx(1.2 * 21.41, rel=1e-6)
    assert r["trl"]["min"] == 5          # подогрев на одном валу (cycles.TRL rh)


def test_catalog_numbers_from_sources():
    from blocks.catalog import catalog
    cat = catalog()
    l0 = T.load_json(T.F_COMB_L0)
    assert cat["blocks"]["combustor"]["options"]["LPP"]["dp_best_annular"] == l0["feed_to_engine_funnel"]["best"]["LPP|кольцевая"]["dp_comb"]
    assert cat["blocks"]["hx"]["recup"]["dp_air"] == 0.025 and cat["blocks"]["hx"]["recup"]["dp_gas"] == 0.035
    assert len(cat["compatibility"]) >= 15
