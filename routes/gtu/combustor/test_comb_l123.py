"""Тесты уровней 1-3 воронки камеры: входы, генераторы кейсов, разбор итогов, сеть акустических элементов, слияние,
пакеты задач (без AWS). Сетка OpenFOAM проверяется, если blockMesh есть в окружении (иначе тест пропускается).
Значения в рядах для разбора и слияния - выдуманные, помечены SYNTHETIC. Запуск: python -m pytest -q test_comb_l123.py
"""
import json, math, os, shutil, subprocess, sys, tempfile

import numpy as np
import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import comb_acoustic as CA  # noqa: E402
import comb_cfd as CF       # noqa: E402
import comb_cfd3d as C3     # noqa: E402
import comb_l1 as L1        # noqa: E402
import comb_post as CP      # noqa: E402
import merge_l1 as M1       # noqa: E402
import prepare_l1 as P      # noqa: E402

DEBUG = os.path.join(HERE, "comb_l0_debug_results.json")


@pytest.fixture(scope="module")
def inputs(tmp_path_factory):
    out = tmp_path_factory.mktemp("l1") / "l1_inputs.json"
    P.main([DEBUG, "--k", "12", "--out", str(out)])
    return json.load(open(out, encoding="utf-8")), str(out)


def test_prepare_deterministic_and_diverse(inputs, tmp_path):
    d, _ = inputs
    assert d["k"] == 12 and len(d["entries"]) == 12
    assert len(d["niches"]) == 12                                   # по кругу ниш: все разные
    out2 = tmp_path / "again.json"
    P.main([DEBUG, "--k", "12", "--out", str(out2)])
    assert json.load(open(out2, encoding="utf-8"))["fingerprint"] == d["fingerprint"]
    for e in d["entries"]:
        assert len(e["loads"]) == 3 and e["loads"][0]["load"] == 1.0
        assert e["loads"][0]["hot_init"]["T"] >= P.T_IGN_MIN


def test_shard_partition():
    for n, S in ((500, 25), (40, 20), (12, 6)):
        allidx = sorted(sum((L1.shard_indices(n, k, S) for k in range(S)), []))
        assert allidx == list(range(n))


def test_cell_mass_flow_matches_element(inputs, tmp_path):
    """осевая скорость на входе канала даёт расход воздуха элемента (клин -> полный элемент)"""
    d, _ = inputs
    for k, e in enumerate(d["entries"]):
        desc = CF.build_cell_case(e, 0, str(tmp_path / f"c{k}"))
        g = CF.cell_geometry(e, 0)
        assert desc["fuel_mode"] == "profile"                       # щель выключена (RING = False, README)
        U = CP.read_field(str(tmp_path / f"c{k}" / "0" / "U"))["boundary"]["inlet_air"]["value"]
        rs, nr = desc["rs"], desc["nr"]
        j0 = 1 if g["r_in"] > 0 else 0
        j1 = rs.index(g["r_ch"]) - 1
        A = [math.pi * (b * b - a * a) for j in range(j0, j1 + 1)
             for a, b in zip(*CF._graded_edges(rs[j], rs[j + 1], nr[j], 1.0))]
        m_air = 0.0
        L = e["loads"][0]
        for (ux, uy, uz), ph, a in zip(U, desc["phi_profile"]["phi"], A):
            Y = CF.mix_Y(ph, e["far_st"])
            m_air += CF.rho_of(L["p3_Pa"], L["T3_K"], Y) * (1 - Y["CH4"]) * uz * a
        assert abs(m_air / g["m_air_main"] - 1) < 1e-6


def test_profile_mean_phi(inputs):
    """профиль смеси: средняя phi по расходу воздуха = phi элемента (с самостратификацией - по доле a_eff)"""
    d, _ = inputs
    for e in d["entries"]:
        for li in range(3):
            g = CF.cell_geometry(e, li)
            r = np.linspace(g["r_in"], g["r_ch"], 4001)[1:] - 0.5 * (g["r_ch"] - g["r_in"]) / 4000
            ph = np.array(CF._profile_phi(e, li, g, r))
            w = r                                                   # площадь кольца ~ r dr
            want = (g["m_fuel_main"] / g["m_air_main"]) / e["far_st"]
            assert abs((ph * w).sum() / w.sum() - want) / want < 2e-3


def test_generators_3d(inputs, tmp_path):
    d, _ = inputs
    for k, e in enumerate(d["entries"][:4]):
        a = C3.build_element3d_case(e, 0, str(tmp_path / f"e{k}"))
        assert a["n_cells"] > 0 and "inlet_fuel" in a["patches"]
        if a["mode"] == "cyclic":
            assert {"cyc_lo", "cyc_hi"} <= set(a["patches"])
        h = e["loads"][0]["hot_init"]
        b = C3.build_liner3d_cht_case(e, dict(T=h["T_eq"], Y=h["Y"]), str(tmp_path / f"l{k}"))
        assert b["n_cells"] > 0 and os.path.exists(tmp_path / f"l{k}" / "Allmesh")


@pytest.mark.skipif(shutil.which("blockMesh") is None, reason="нет OpenFOAM")
def test_blockmesh_checkmesh(inputs, tmp_path):
    d, _ = inputs
    e = d["entries"][6]
    case = tmp_path / "c"
    CF.build_cell_case(e, 0, str(case), of="1912")
    env = dict(os.environ)
    subprocess.run(["blockMesh"], cwd=case, capture_output=True, env=env, timeout=300)
    r = subprocess.run(["checkMesh"], cwd=case, capture_output=True, text=True, env=env, timeout=300)
    if "etc entry" in (r.stdout + r.stderr):
        pytest.skip("окружение OpenFOAM не подключено")
    assert "Mesh OK." in r.stdout
    m = CP.Mesh(str(case))                                         # объём клина = аналитика
    g = CF.cell_geometry(e, 0)
    V = CF.TAN_A * (g["r_ch"] ** 2 * g["l_ch"] + g["R_c"] ** 2 * g["L_cell"])
    assert abs(m.V.sum() / V - 1) < 1e-4


def test_emissions_15pct_SYNTHETIC():
    X = {"H2O": 0.10, "O2": 0.135, "NO": 9e-6, "CO": 4.5e-6, "CO2": 0.04}     # SYNTHETIC
    em = CP.emissions(X)
    o2d = 100 * 0.135 / 0.9
    assert abs(em["O2_dry_pct"] - o2d) < 1e-12
    assert abs(em["NOx_ppm15"] - 1e6 * 9e-6 / 0.9 * (20.9 - 15) / (20.9 - o2d)) < 1e-9


def test_phi_atoms_fresh_and_burnt():
    Y = CF.mix_Y(0.6, 0.0580163)
    assert abs(CP.phi_atoms(Y) - 0.6) < 2e-3                         # GRI-веса против таблицы MW генератора
    # сгоревшая бедная смесь phi 0,6: CO2 + 2 H2O на CH4, остаток O2
    nC = Y["CH4"] / 16.043
    Yb = dict(Y)
    Yb["CH4"] = 0.0
    Yb["CO2"] = Y["CO2"] + nC * 44.009
    Yb["H2O"] = 2 * nC * 18.015
    Yb["O2"] = Y["O2"] - 2 * nC * 31.998
    assert abs(CP.phi_atoms(Yb) - 0.6) < 2e-3


def test_acoustic_selfcheck():
    r = CA.selfcheck()
    assert r["all_ok"], r


def test_ftf_identification_SYNTHETIC():
    rng = np.random.default_rng(1)
    dt, n, tau = 1e-4, 20000, 1.2e-3                                  # SYNTHETIC: n - tau, n = 0,8
    u = rng.standard_normal(n)
    k = int(round(tau / dt))
    q = np.zeros(n)
    q[k:] = 0.8 * u[:-k]
    r = CA.identify_ftf(u, q, dt, n_taps=32, freqs=[100.0, 300.0])
    for f, F in zip(r["f"], r["F"]):
        assert abs(abs(F) - 0.8) < 0.02
        assert abs(np.angle(F * np.exp(1j * 2 * np.pi * f * tau))) < 0.05


def _part(tmp, shard, shards, n, fp="f" * 16, code="c" * 16):
    rows = [dict(pos=i, idx=i, why="SYNTHETIC", cases={}, status="ошибка") for i in range(shard, n, shards)]
    p = dict(schema=L1.SCHEMA, kv="KV-3", n=n, shard=shard, shards=shards, inputs_fingerprint=fp, code_hash=code,
             level="coarse", iters=CF.ITERS, generate_only=False, rows=rows, runtime_s=0.0)
    f = os.path.join(tmp, f"comb_l1_KV-3_n{n}_s{shard:03d}of{shards:03d}.json")
    json.dump(p, open(f, "w"))
    return f


def test_merge_completeness_SYNTHETIC(tmp_path):
    for k in range(4):
        _part(str(tmp_path), k, 4, 10)
    ok, rep, rows = M1.check_complete(M1.load_parts([str(tmp_path)]))
    assert ok and len(rows) == 10
    os.remove(os.path.join(tmp_path, "comb_l1_KV-3_n10_s002of004.json"))
    ok, rep, _ = M1.check_complete(M1.load_parts([str(tmp_path)]))
    assert not ok and rep["missing_shards"] == [2]
    _part(str(tmp_path), 2, 4, 10, fp="0" * 16)
    ok, rep, _ = M1.check_complete(M1.load_parts([str(tmp_path)]))
    assert not ok and any("inputs_fingerprint" in x for x in rep["problems"])


def test_jobs_files_load_offline():
    sys.path.insert(0, os.path.join(HERE, "..", "platform", "aws"))
    import runner
    for f in ("jobs_calib.json", "jobs_l1.json", "jobs_l2.json", "jobs_l3.json"):
        p = os.path.join(HERE, f)
        if not os.path.exists(p):
            pytest.skip("пакеты ещё не собраны (make_jobs_l123.py)")
        jobs = runner.load_jobs(p)
        assert jobs and len({j.name for j in jobs}) == len(jobs)
        d = json.load(open(p, encoding="utf-8"))
        if f != "jobs_calib.json":
            assert 10 <= len(d["jobs"]) <= 50                         # правило прогонов: 10-50 частей
        assert d["plan_check"]["exit_code"] == 0
