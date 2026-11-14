"""Проверки уровня 0, часть 2 (funnel_l0, merge_l0, make_jobs_l0): повторяемость, сумма весов, монотонность оценки по
деталям, деление на части без потерь и задвоений (объединение частей = расчёт целиком, при любом числе процессов),
отказ слияния при потере, задвоении и чужом отпечатке, отбор по нишам без потери допустимой ниши, файлы задач AWS.
Синтетические таблицы - придуманные корпуса (не данные установок). Запуск: python3 -m pytest -q test_funnel_l0.py"""
import json, math, os, shutil, sys, tempfile

import numpy as np
import pytest

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import funnel_l0 as L0          # noqa: E402
import funnel_space as FS       # noqa: E402
import merge_l0 as M            # noqa: E402
import parts_model as PM        # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))


@pytest.fixture(scope="module")
def tmp():
    d = tempfile.mkdtemp(prefix="test_l0_")
    yield d
    shutil.rmtree(d, ignore_errors=True)


# =========================================================================== повторяемость
def test_record_repeatable_and_order_free():
    a = [L0.evaluate_l0(i, 64, "KV-3") for i in range(6)]
    FS._DM.clear()
    b = [L0.evaluate_l0(i, 64, "KV-3") for i in reversed(range(6))][::-1]
    assert json.dumps(a, ensure_ascii=False) == json.dumps(b, ensure_ascii=False)
    assert list(a[0]) == L0.COLUMNS


def test_concept_does_not_depend_on_n():
    """основание сходимости: первые N/2 концепций выборки N - это выборка N/2"""
    for i in (0, 5, 63):
        c1, c2 = FS.concept(i, 64, "KV-3"), FS.concept(i, 128, "KV-3")
        c1.pop("n"); c2.pop("n")
        assert c1 == c2


def test_parse_shards():
    assert L0.parse_shards("3", 10) == [3]
    assert L0.parse_shards("0,3, 6", 10) == [0, 3, 6]
    assert L0.parse_shards("2-4,9", 10) == [2, 3, 4, 9]
    assert L0.parse_shards("all", 4) == [0, 1, 2, 3]
    for bad in ("10", "-1", ""):
        with pytest.raises(ValueError):
            L0.parse_shards(bad, 10)


# =========================================================================== веса и оценка
def test_weights_sum_and_scenarios():
    assert abs(sum(L0.WEIGHTS.values()) - 1.0) < 1e-12
    assert max(L0.WEIGHTS, key=L0.WEIGHTS.get) == "parts"               # главное - минимум деталей (0,35)
    assert L0.WEIGHTS == {"parts": 0.35, "plugplay": 0.25, "service": 0.20, "effmake": 0.20}
    assert abs(sum(L0.SPLIT_EFFMAKE.values()) - 1) < 1e-12 and abs(sum(L0.SPLIT_ETASIZE.values()) - 1) < 1e-12
    for s in M.scenarios(L0.WEIGHTS, L0.SPLIT_EFFMAKE):
        assert abs(sum(s["w"].values()) - 1.0) < 1e-12, s["name"]
        assert abs(sum(s["split"].values()) - 1.0) < 1e-12, s["name"]
    for bad in ({"parts": 0.5, "plugplay": 0.25, "service": 0.2, "effmake": 0.2},
                {"parts": 0.35, "plugplay": 0.25, "service": 0.4},
                {"parts": -0.1, "plugplay": 0.5, "service": 0.4, "effmake": 0.2}):
        with pytest.raises(ValueError):
            L0.check_weights(bad)
    with pytest.raises(ValueError):
        M.parse_weights("parts=0.4,plugplay=0.25,service=0.2,effmake=0.2")
    assert M.parse_weights("parts=0.35,plugplay=0.25,service=0.20,effmake=0.20") == L0.WEIGHTS


def test_score_monotone_in_parts():
    """больше деталей (N_eff) -> ниже балл деталей -> ниже взвешенная оценка при прочих равных"""
    base = dict(plugplay=6.0, service=5.0, effmake=4.0)
    prev = None
    for n_eff in (50, 100, 200, 400, 800, 1600):
        s = dict(base, parts=PM.score_from_neff(n_eff))
        v = L0.weighted(s)
        if prev is not None:
            assert v < prev
        prev = v
    s0 = dict(base, parts=5.0)
    s1 = dict(base, parts=5.1)
    assert math.isclose(L0.weighted(s1) - L0.weighted(s0), 0.035, rel_tol=1e-9)


def test_s_eta_s_size_scales():
    assert L0.s_eta_of(31.5, 31.5) == 5.0 and L0.s_eta_of(36.5, 31.5) == 10.0 and L0.s_eta_of(26.5, 31.5) == 0.0
    assert L0.s_eta_of(40, 31.5) == 10.0 and L0.s_eta_of(10, 31.5) == 0.0
    assert L0.s_size_of(0.0, 10.0) == 10.0 and L0.s_size_of(5.0, 10.0) == 5.0 and L0.s_size_of(20.0, 10.0) == 0.0


def test_calibration_ns_from_frame3():
    cal = L0.calib_ns()
    assert 0.5 < cal["ns_opt_c"] < 1.5 and 0.02 < cal["k_ns_c"] < 0.2
    for p in cal["points"]:                       # подбор по двум ступеням frame3 - в пределах 0,5 п. п.
        assert abs(p["fit"] - p["d_eta"]) < 0.005
    assert cal["ns_opt_d"] > cal["ns_opt_c"]


# =========================================================================== части и слияние
N_T, S_T = 40, 4


@pytest.fixture(scope="module")
def parts(tmp):
    """выборка KV-3, N=40: 4 части пулом из 2 процессов (одна часть - отдельным вызовом) и одна часть целиком в одном процессе"""
    pdir, wdir = os.path.join(tmp, "parts"), os.path.join(tmp, "whole")
    L0.run_parts("KV-3", N_T, [0, 1, 2], S_T, pdir, workers=2, gz=True)
    L0.run_parts("KV-3", N_T, [3], S_T, pdir, workers=1)
    L0.run_parts("KV-3", N_T, [0], 1, wdir, workers=1)
    return pdir, wdir


def test_parts_union_equals_whole(parts):
    pdir, wdir = parts
    fa = M.find_parts([pdir], "KV-3", N_T, S_T)
    assert [k for k, _ in fa] == list(range(S_T))
    ca, ra, _, ch, _ = M.merge_parts(fa, "KV-3", N_T, S_T)
    cb, rb, _, _, _ = M.merge_parts(M.find_parts([wdir], "KV-3", N_T, 1), "KV-3", N_T, 1)
    assert ch["indices_ok"] and ch["n_rows"] == N_T and not ch["missing"]
    assert ca == cb and json.dumps(ra, ensure_ascii=False) == json.dumps(rb, ensure_ascii=False)


def test_merge_outputs_identical_for_parts_and_whole(parts, tmp):
    pdir, wdir = parts
    kw = dict(k=2, top=5, k_near=1, near_max=0.5, target=None, detail_top=0, cycle_k=0)
    a = M.run_merge(M.find_parts([pdir], "KV-3", N_T, S_T), "KV-3", N_T, S_T, os.path.join(tmp, "ma"), **kw)
    b = M.run_merge(M.find_parts([wdir], "KV-3", N_T, 1), "KV-3", N_T, 1, os.path.join(tmp, "mb"), **kw)
    assert open(a["csv"], encoding="utf-8").read() == open(b["csv"], encoding="utf-8").read()
    assert a["sel"]["selected"] == b["sel"]["selected"]
    assert a["res"]["top"] == b["res"]["top"]


def _copy_parts(src, dst):
    shutil.rmtree(dst, ignore_errors=True)
    shutil.copytree(src, dst)
    return dst


def test_merge_rejects_lost_duplicated_and_foreign(parts, tmp):
    pdir, _ = parts
    # потеря части
    d = _copy_parts(pdir, os.path.join(tmp, "lost"))
    os.remove(L0.part_path(d, "KV-3", N_T, 2, S_T, gz=True))
    with pytest.raises(M.MergeError, match="missing"):
        M.merge_parts(M.find_parts([d], "KV-3", N_T, S_T), "KV-3", N_T, S_T)
    # совпадающая копия части (повтор задачи) - допускается; копия с другими записями - отказ
    d = _copy_parts(pdir, os.path.join(tmp, "dup"))
    os.makedirs(os.path.join(d, "attempt2"))
    src = L0.part_path(d, "KV-3", N_T, 3, S_T)
    shutil.copy(src, os.path.join(d, "attempt2", os.path.basename(src)))
    _, rows, _, ch, _ = M.merge_parts(M.find_parts([d], "KV-3", N_T, S_T), "KV-3", N_T, S_T)
    assert len(ch["duplicate_same"]) == 1 and len(rows) == N_T
    doc = json.load(open(src, encoding="utf-8"))
    doc["rows"][0][doc["columns"].index("eta_el")] = 99.0
    json.dump(doc, open(os.path.join(d, "attempt2", os.path.basename(src)), "w", encoding="utf-8"))
    with pytest.raises(M.MergeError, match="duplicate_conflict"):
        M.merge_parts(M.find_parts([d], "KV-3", N_T, S_T), "KV-3", N_T, S_T)
    # чужой отпечаток выборки (другая версия scipy на машине) и чужой индекс внутри части
    d = _copy_parts(pdir, os.path.join(tmp, "fp"))
    doc = json.load(open(L0.part_path(d, "KV-3", N_T, 3, S_T), encoding="utf-8"))
    doc["fingerprint"] = "0" * 16
    json.dump(doc, open(L0.part_path(d, "KV-3", N_T, 3, S_T), "w", encoding="utf-8"))
    with pytest.raises(M.MergeError, match="meta_mismatch"):
        M.merge_parts(M.find_parts([d], "KV-3", N_T, S_T), "KV-3", N_T, S_T)
    d = _copy_parts(pdir, os.path.join(tmp, "idx"))
    doc = json.load(open(L0.part_path(d, "KV-3", N_T, 3, S_T), encoding="utf-8"))
    doc["rows"][1][0] = 0                                              # индекс 0 - из части 0: задвоение
    json.dump(doc, open(L0.part_path(d, "KV-3", N_T, 3, S_T), "w", encoding="utf-8"))
    with pytest.raises(M.MergeError, match="index_errors"):
        M.merge_parts(M.find_parts([d], "KV-3", N_T, S_T), "KV-3", N_T, S_T)
    # обрезанный файл прерванной попытки рядом с целым - не мешает
    d = _copy_parts(pdir, os.path.join(tmp, "trunc"))
    os.makedirs(os.path.join(d, "_lost"))
    src = L0.part_path(d, "KV-3", N_T, 3, S_T)
    open(os.path.join(d, "_lost", os.path.basename(src)), "w").write(open(src).read()[:100])
    _, rows, _, ch, _ = M.merge_parts(M.find_parts([d], "KV-3", N_T, S_T), "KV-3", N_T, S_T)
    assert len(ch["unreadable"]) == 1 and len(rows) == N_T


# =========================================================================== отбор по нишам
def synth_table(seed=7, n=400, n_niche=60, p_feas=0.3):
    """придуманные концепции: ниши разного размера, оценки и превышения - случайные с фиксированным зерном"""
    rng = np.random.default_rng(seed)
    rows = []
    for i in range(n):
        r = {c: None for c in L0.COLUMNS}
        nk = int(rng.integers(0, n_niche) ** 1.0) if i >= n_niche else i            # каждая ниша хотя бы с одной концепцией
        feas = bool(rng.random() < p_feas)
        s = rng.uniform(1, 9, 4)
        r.update(idx=i, kv="KV-3", arch=f"арх{nk % 7}", niche=f"ниша{nk:03d}", feasible=int(feas), violation=0.0 if feas else
                 float(rng.uniform(0.01, 1.0)), s_parts=s[0], s_plugplay=s[1], s_service=s[2], s_effmake=s[3], s_etasize=s[3],
                 s_make=s[3], score=0.0, reject="" if feas else "soft_G")
        rows.append([r[c] for c in L0.COLUMNS])
    return M.Table(L0.COLUMNS, rows)


@pytest.mark.parametrize("k,top,target", [(1, 0, None), (3, 10, None), (1, 5, 10), (2, 50, 1000), (5, 0, 40)])
def test_selection_keeps_every_feasible_niche(k, top, target):
    T = synth_table()
    sc = M.scores(T, L0.WEIGHTS)
    fpos = np.nonzero(T.feas)[0]
    lay = np.zeros(T.n, dtype=int); lay[fpos] = M.pareto_layers(T.S[fpos])
    S = M.select(T, sc, k=k, top=top, pareto_pos=[p for p in fpos if lay[p] == 1], k_near=1, near_max=0.5, target=target)
    feas_niches = {T.niche[p] for p in fpos}
    got = {T.niche[p] for p in S["sel"] if T.feas[p]}
    assert feas_niches and feas_niches <= got
    assert set(S["feasible_niches"]) == feas_niches
    if target is not None and not S["target_exceeded"]:
        assert len(S["sel"]) <= target
    for p, why in S["sel"].items():                       # почти допустимые - только из ниш без допустимых
        if not T.feas[p]:
            assert T.niche[p] not in feas_niches and T.viol[p] <= 0.5 and any("без допустимых" in w for w in why)
    # лучший в каждой нише - выбран
    for nk in feas_niches:
        best = max((p for p in fpos if T.niche[p] == nk), key=lambda p: (sc[p], -T.idx[p]))
        assert best in S["sel"]


def test_target_picks_largest_k():
    T = synth_table(n=600, n_niche=40, p_feas=0.6)
    sc = M.scores(T, L0.WEIGHTS)
    sizes = [len(M.select(T, sc, k=K, top=0)["sel"]) for K in range(1, 8)]
    assert sizes == sorted(sizes)
    S = M.select(T, sc, top=0, target=sizes[3])
    assert S["K"] == 4 and len(S["sel"]) == sizes[3]


def test_pareto_layers_small():
    F = np.array([[3, 3], [1, 4], [4, 1], [2, 2], [1, 1], [3, 3]], dtype=float)
    lay = M.pareto_layers(F, 3)
    assert list(lay) == [1, 1, 1, 2, 3, 1]


def test_selection_on_real_sample_keeps_niches(parts, tmp):
    pdir, _ = parts
    r = M.run_merge(M.find_parts([pdir], "KV-3", N_T, S_T), "KV-3", N_T, S_T, os.path.join(tmp, "mr"), k=1, top=1, k_near=0,
                    target=None, detail_top=0, cycle_k=0)
    s = r["sel"]["summary"]
    assert s["lost_feasible_niches"] == [] and s["feasible_niches_covered"] == s["feasible_niches"]


# =========================================================================== файлы задач AWS
@pytest.mark.parametrize("fname", ["jobs_l0_debug.json", "jobs_l0.json"])
def test_jobs_files_cover_every_shard_once(fname):
    path = os.path.join(HERE, fname)
    if not os.path.exists(path):
        pytest.skip("файл задач ещё не собран (make_jobs_l0.py)")
    sys.path.insert(0, os.path.join(os.path.dirname(HERE), "platform", "aws"))
    import runner                                                     # noqa: E402
    jobs = runner.load_jobs(path)
    doc = json.load(open(path, encoding="utf-8"))
    assert len(jobs) <= 6                                             # одна волна под квоту 384 vCPU
    seen = {}
    for j in jobs:
        assert j.python == "3.12" and "scipy==1.18.1" in j.pip
        for line in j.command.splitlines():
            if "funnel_l0.py" not in line:
                continue
            tok = line.split()
            kv, n, S = tok[tok.index("--kv") + 1], int(tok[tok.index("--n") + 1]), int(tok[tok.index("--shards") + 1])
            assert '"$NPROC"' in line
            for k in L0.parse_shards(tok[tok.index("--shard") + 1], S):
                assert (kv, n, k) not in seen, (kv, n, k)
                seen[(kv, n, k)] = j.name
    for kv, n in {(kv, n) for kv, n, _ in seen}:
        S = max(k for a, b, k in seen if (a, b) == (kv, n)) + 1
        assert sorted(k for a, b, k in seen if (a, b) == (kv, n)) == list(range(S))
        assert any(f"--kv {kv} --n {n} --shards {S}" in c for c in doc["merge"]["commands"])
