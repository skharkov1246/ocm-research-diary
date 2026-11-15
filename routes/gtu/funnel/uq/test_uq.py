"""Тесты устойчивого отбора и фронтов Парето ниш: python3 -m pytest -q routes/gtu/funnel/uq/test_uq.py (около 2 мин на 4 ядрах).

Повторяемость, восстановление моделей после подстановки, номинальная выборка = расчётная точка уровня 1, сходимость Монте-Карло
(P10/P90 при удвоении выборки), оценка индексов Соболя на аналитической функции, NSGA-II на ZDT1 (известный фронт) и ZDT1 с
ограничением, гиперобъём, деление на части и слияние с контролем полноты, кодирование генов ниши, логика устойчивого отбора,
файлы задач (формат раннера, пути команд).
"""
import copy, glob, gzip, json, math, os, shutil, sys, tempfile

import numpy as np
import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
FUNNEL = os.path.dirname(HERE)
for p in (FUNNEL, HERE):
    if p not in sys.path:
        sys.path.insert(0, p)
import funnel_space as FS        # noqa: E402
import funnel_l1 as L1           # noqa: E402
import nsga2 as NS               # noqa: E402
import uq_factors as UF          # noqa: E402
import uq_patch as UP            # noqa: E402
import uq_stats as US            # noqa: E402
import uq_l1 as U                # noqa: E402
import merge_uq as MU            # noqa: E402

PROBE = os.path.join(FUNNEL, "funnel_l1.csv")
need_probe = pytest.mark.skipif(not os.path.exists(PROBE), reason="нет пробы уровня 1 funnel_l1.csv")


@pytest.fixture(scope="module")
def cands():
    import uq_probe as P
    return P.candidates_from_csv()


# =========================================================================== выборка и статистика
def test_lhs_stratified_and_repeatable():
    u = UF.lhs(50, 4, 123)
    for k in range(4):
        assert sorted(np.floor(u[:, k] * 50).astype(int)) == list(range(50))
    assert np.array_equal(u, UF.lhs(50, 4, 123))
    a, b = UF.sample(40), UF.sample(40)
    assert np.array_equal(a[1], b[1]) and UF.fingerprint(40) == UF.fingerprint(40)
    _, X2, _, d = UF.sample(40, dummy=True)
    assert np.array_equal(X2, a[1]) and d.shape == (40,)                 # пустышка не меняет первые столбцы


def test_factor_ranges_and_nominal():
    _, X, vals = UF.sample(400)
    for k, (name, grp, dist, par, nom, *_r) in enumerate(UF.FACTORS):
        col = X[:, k]
        if dist in ("unif",):
            assert par[0] - 1e-12 <= col.min() and col.max() <= par[1] + 1e-12
        if dist == "tri":
            assert par[0] - 1e-12 <= col.min() and col.max() <= par[2] + 1e-12
        if dist in ("norm",):
            assert par[2] - 1e-12 <= col.min() and col.max() <= par[3] + 1e-12
        if dist == "disc":
            assert set(np.unique(col)) <= set(par[0])
    assert UF.is_nominal(UF.nominal())
    assert UF.rough_factor("печать", 2.0) == 2.0 and UF.rough_factor("поковка", 2.0) == 1.0
    assert abs(UF.rough_factor("литьё", 2.0) - (1 + UF.K_CAST)) < 1e-12


def test_sobol_analytic():
    rng = np.random.Generator(np.random.PCG64(5))
    X = rng.random((4000, 3))
    y = X[:, 0] + 0.5 * X[:, 1]                    # S1 = 1 / 1,25 = 0,8; S2 = 0,2; S3 = 0
    S = US.sobol1(X, y)
    assert abs(S[0] - 0.8) < 0.04 and abs(S[1] - 0.2) < 0.04 and S[2] < 0.02
    assert US.sobol1(X, np.ones(4000)) is None
    assert US.wilson(5, 10)[0] < 0.5 < US.wilson(5, 10)[1]


# =========================================================================== подстановка в модели
@need_probe
def test_patch_restores_state(cands):
    st0 = UP.state()
    _, _, vals = UF.sample(6)
    tgt = FS.target("KV-3")
    for v in vals:
        with UP.Perturb(v):
            U.eval_one(cands[7], tgt, v)
    assert UP.state() == st0
    with pytest.raises(RuntimeError):
        with UP.Perturb(vals[0]):
            with UP.Perturb(vals[1]):
                pass
    assert UP.state() == st0


@need_probe
def test_nominal_reproduces_level1(cands):
    """номинальная выборка = расчётная точка уровня 1 (оценка и допустимость строк пробы, с узлами ВНА)"""
    tgt = FS.target("KV-3")
    with UP.Perturb(UF.nominal()):
        for cd in cands:
            r = U.eval_one(cd, tgt, UF.nominal())
            assert abs(r["score_w"] - cd["rec"]["score"]) < 1e-3, cd["idx"]
            assert r["feas"] == int(cd["rec"]["feasible"]), cd["idx"]


@need_probe
def test_factors_move_models_in_right_direction(cands):
    tgt = FS.target("KV-3")
    cd = cands[7]
    base = dict(UF.nominal())

    def run(**kw):
        v = dict(base, **kw)
        with UP.Perturb(v):
            return U.eval_one(cd, tgt, v)
    r0 = run()
    assert run(d_eta_c=-0.02)["eta_el"] < r0["eta_el"]
    assert run(d_eta_t=-0.015)["eta_el"] < r0["eta_el"]
    assert run(k_cl_c=1.5)["eta_el"] < r0["eta_el"]
    assert run(k_cl_t=1.5)["eta_el"] < r0["eta_el"]
    assert run(k_rough=2.5)["eta_el"] < r0["eta_el"]               # литое колесо: трение × (1 + 0,3·1,5)
    assert run(k_lmp=0.8)["life_min"] < r0["life_min"]
    assert run(dT_gas=30.0)["life_min"] < r0["life_min"]
    assert run(k_uts=0.9)["burst_min"] < r0["burst_min"]
    assert run(d_trl=-1)["score_f"] < r0["score_f"] < run(d_trl=1)["score_f"]


@need_probe
def test_mc_repeatable_and_converges(cands):
    """повторяемость (тот же результат при 1 и 2 процессах) и сходимость: P10/P50/P90 оценки и P(расчётная) при удвоении N"""
    tgt = FS.target("KV-3")
    sub = [cands[7], cands[8]]
    _, X, vals, d = UF.sample(24, dummy=True)
    a, na, _ = U.run_dp(sub, [0, 1], vals, tgt, 1)
    b, nb, _ = U.run_dp(sub, [0, 1], vals, tgt, 2)
    assert json.dumps(a, sort_keys=True, default=str) == json.dumps(b, sort_keys=True, default=str)
    s = {}
    for n in (100, 200):
        _, X, vals, d = UF.sample(n, dummy=True)
        r, nom, _ = U.run_dp([cands[7]], [0], vals, tgt, 2)
        s[n] = U.summarize_dp(r[0], X, nom[0], d)
    for i in (0, 1, 2):
        assert abs(s[100]["score_f"][i] - s[200]["score_f"][i]) < 0.2, (i, s[100]["score_f"], s[200]["score_f"])
    assert abs(s[100]["P_dp"] - s[200]["P_dp"]) < 0.1


# =========================================================================== части и слияние
@need_probe
def test_shards_merge_and_completeness(cands):
    tmp = tempfile.mkdtemp(prefix="uq_test_")
    try:
        sub = cands[6:9]
        meta = dict(kind="L1-csv", kv="KV-3")
        for k in (0, 1):
            U.run_parts(sub, meta, [k], 2, tmp, n=16, n_od=0, od_rule="none", workers=1, gz=True)
        found = MU.find_parts([tmp], "KV-3", 3, 2)
        assert [k for k, _ in found] == [0, 1]
        assert MU.detect_k(MU.find_parts([tmp], "KV-3", 0, 2)) == 3
        rows, ref, checks, calls = MU.merge_parts(found, "KV-3", 3, 2)
        assert [r["u"] for r in rows] == [0, 1, 2] and checks["indices_ok"]
        with pytest.raises(MU.MergeError):
            MU.merge_parts(found[:1], "KV-3", 3, 2)                      # нет части - отказ
        dup = os.path.join(tmp, "copy")
        os.makedirs(dup)
        shutil.copy(found[0][1], dup)                                    # повтор той же части - допустим
        rows2, _, ch2, _ = MU.merge_parts(MU.find_parts([tmp], "KV-3", 3, 2), "KV-3", 3, 2)
        assert ch2["duplicate_same"] and len(rows2) == 3
        d = json.load(gzip.open(found[0][1], "rt", encoding="utf-8"))
        d["rows"][0]["dp"]["P_dp"] = 0.999                                # повтор с другим содержанием - отказ
        with gzip.open(os.path.join(dup, os.path.basename(found[0][1])), "wt", encoding="utf-8") as fh:
            json.dump(d, fh)
        with pytest.raises(MU.MergeError):
            MU.merge_parts(MU.find_parts([tmp], "KV-3", 3, 2), "KV-3", 3, 2)
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


def _row(u, arch, niche, p_feas, p50, score, feasible=1, trl=6, t41=1000.0):
    return dict(u=u, src="L1", j_src=u, key=None, idx=u, arch=arch, niche=niche, P_feas=p_feas, robust=(p50 * p_feas if p_feas is not None else None),
                od_status="посчитаны", od=dict(od_fail_freq={}, P_od_dp=1.0),
                dp=dict(P_dp=p_feas or 0.0, score_f=[p50 - 0.3, p50, p50 + 0.3], viol_freq={"soft_G": 0.4}, rules=[], nominal=dict(feas=1, score_f=p50),
                        sobol=dict(y=[["d_eta_c", 0.3]])),
                rec=dict(feasible=feasible, score=score, trl_min=trl, t41=t41, pi=8.0, eps_rec=0.0, eta_th_shaft=30.0, V_env=8.0, mass_kg=1500.0,
                         arch=arch, **{a: (arch if a == "arch" else "x") for a in MU.M1.DISC_W if a != "arch"}))


def test_robust_selection_logic():
    rows = [_row(0, "A", "n0", 0.10, 6.5, 6.9), _row(1, "B", "n1", 0.90, 5.6, 6.0), _row(2, "C", "n2", 0.60, 5.4, 5.8),
            _row(3, "D", "n3", None, 6.0, 6.5)]
    A = MU.select(MU.std_records(rows))
    C = MU.select(MU.robust_records(rows), score_min=0.0)
    assert [s["j"] for s in A][0] == 0                       # обычный: лучший номинал первым
    assert [s["j"] for s in C] == [1, 2]                      # устойчивый: P >= 0,5 и P50 >= 5, порядок по P50 × P
    cmp_ = MU.compare(rows, A, MU.select(MU.std_records(rows, True)), C, MU.P_MIN, U.SCORE_MIN)
    assert cmp_["changed_vs_A"]["out"] == 2 and all(x["why"] for x in cmp_["out_vs_B"])
    assert "P(допустим) = 0.10" in [x for x in cmp_["out_vs_B"] if x["u"] == 0][0]["why"]


# =========================================================================== NSGA-II
def test_nsga2_zdt1_front():
    r = NS.run(NS.zdt1, 30, 2, pop=100, gens=200, seed=1)
    m = NS.front_mask(r["F"], r["V"])
    assert m.sum() >= 50
    assert NS.igd(r["F"][m], NS.zdt1_front()) < 0.01
    assert NS.hv2(r["F"][m], (1.1, 1.1)) > 0.98 * NS.hv2(NS.zdt1_front(3000), (1.1, 1.1))
    r2 = NS.run(NS.zdt1, 30, 2, pop=100, gens=200, seed=1)
    assert np.array_equal(r["X"], r2["X"])                    # детерминированно


def test_nsga2_constrained_zdt1():
    r = NS.run(NS.zdt1_constrained, 30, 2, pop=100, gens=200, seed=2)
    m = NS.front_mask(r["F"], r["V"])
    assert m.sum() >= 50 and r["X"][m, 0].min() >= 0.3 - 1e-9
    assert NS.igd(r["F"][m], NS.zdt1_front(500, 0.3)) < 0.01


def test_hv_mc():
    G = np.array([[0.2, 0.5, 0.3], [0.5, 0.2, 0.4]])        # объединение коробок: 0,28 + 0,24 − 0,15 = 0,37
    assert abs(NS.hv_mc(G, 40000) - 0.37) < 0.01


# =========================================================================== ниша
SEL_DEBUG = os.path.join(FUNNEL, "l0_aws_debug", "funnel_l0_select_KV-3_n2000.json")


@pytest.mark.skipif(not os.path.exists(SEL_DEBUG), reason="нет отбора уровня 0 отладки")
def test_niche_encode_build_roundtrip():
    import pareto_niche as PN
    rows, meta = PN.load_niches(SEL_DEBUG)
    L1._prep(meta["n"], meta["seed"], meta["kv"])
    for r in rows[:3]:
        N = PN.Niche(r["niche"], meta["kv"], r["ref_idx"], meta["n"], meta["seed"])
        for idx in r["l0_idx"][:3]:
            c = FS.concept(idx, meta["n"], meta["kv"], meta["seed"])
            c2 = N.build(N.encode(c))
            assert c2 is not None and PN.L0.niche_of(c2) == r["niche"]
            for a in PN.FREE_DISC:
                assert c2.get(a) == c.get(a), (idx, a)
            ev, ev2 = L1.evaluate_point(c, N.tgt), L1.evaluate_point(c2, N.tgt)
            assert abs(ev["score"] - ev2["score"]) < 1e-6 and abs(ev["viol"] - ev2["viol"]) < 1e-6


@pytest.mark.skipif(not os.path.exists(SEL_DEBUG), reason="нет отбора уровня 0 отладки")
def test_pareto_parts_merge_repeatable():
    import pareto_niche as PN, merge_pareto as MP
    tmp = tempfile.mkdtemp(prefix="par_test_")
    try:
        for k in (0, 1):
            PN.run_parts(SEL_DEBUG, [k], 2, tmp, pop=10, gens=2, check=0, workers=1, top=2)
        found = MP.find_parts([tmp], "KV-3", 2, 2)
        out = MP.run_merge(found, "KV-3", 2, 2, tmp, "_t")
        assert out["res"]["merge_checks"]["indices_ok"] and out["res"]["summary"]["niches"] == 2
        with pytest.raises(MP.MergeError):
            MP.merge_parts(found[1:], "KV-3", 2, 2)
        t2 = os.path.join(tmp, "again")
        PN.run_parts(SEL_DEBUG, [0], 2, t2, pop=10, gens=2, check=0, workers=1, top=2)
        a = json.load(open(found[0][1]))["rows"][0]
        b = json.load(open(glob.glob(os.path.join(t2, "*.json"))[0]))["rows"][0]
        for key in ("front", "hv", "evals", "feasible_evals"):
            assert a[key] == b[key]
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


# =========================================================================== файлы задач
@pytest.mark.parametrize("name", ["jobs_uq_debug.json", "jobs_uq.json", "jobs_pareto.json", "jobs_uq_pareto.json"])
def test_jobs_files(name):
    path = os.path.join(HERE, name)
    if not os.path.exists(path):
        pytest.skip("файл задач ещё не собран (make_jobs_uq.py)")
    sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(FUNNEL)), "gtu", "platform", "aws"))
    import runner
    jobs = runner.load_jobs(path)
    repo = os.path.dirname(os.path.dirname(os.path.dirname(FUNNEL)))
    d = json.load(open(path, encoding="utf-8"))
    for j in d["jobs"]:
        first = j["command"].split("\n")[0].split()
        assert os.path.exists(os.path.join(repo, first[1])), first[1]
        if "--select" in first:
            assert os.path.exists(os.path.join(repo, first[first.index("--select") + 1]))
        assert j["command"].rstrip().endswith("exit $rc")
    assert d["estimate"]["upper_bound_usd"] > d["estimate"]["usd"]["total_expected"] > 0
    assert len(jobs) == len(d["jobs"])


@need_probe
def test_offdesign_path_nominal(cands):
    """нерасчётные режимы под подстановкой: номинал допустимого варианта пробы проходит их, ВНА - как у уровня 1 (около 30 с)"""
    tgt = FS.target("KV-3")
    cd = cands[7]
    with UP.Perturb(UF.nominal()):
        r = U.eval_one(cd, tgt, UF.nominal(), od=True)
    assert r["feas"] == 1 and r["od_ok"] == 1 and r["vigv"] == int(cd["rec"]["vigv"])
    assert abs(r["score_w"] - cd["rec"]["score"]) < 1e-3
    sm = dict(P_dp=0.5)
    od = U.summarize_od(sm, {0: dict(r), 3: dict(r, od_ok=0, od_bad=["sm_min"])}, 4)
    assert od["n_run"] == 2 and od["P_od_dp"] == 0.5 and od["od_fail_freq"] == {"sm_min": 0.5}
