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
    FS._DM.clear(); FS._CLOSURE.clear(); FS._CANDS.clear(); FS._NOCLOSE.clear()   # запись не зависит от кэшей замыкания и порядка счёта
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
    assert ch["indices_ok"] and ch["n_concepts"] == N_T and N_T <= ch["n_rows"] <= 2 * N_T and not ch["missing"]   # v3: 1-2 записи на концепцию
    assert ca == cb and json.dumps(ra, ensure_ascii=False) == json.dumps(rb, ensure_ascii=False)


def test_merge_outputs_identical_for_parts_and_whole(parts, tmp):
    pdir, wdir = parts
    kw = dict(k=2, top=5, k_near=1, near_max=0.5, target=None, detail_top=0, cycle_k=0, closure_proof=False)
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
    assert len(ch["duplicate_same"]) == 1 and len({r[0] for r in rows}) == N_T == ch["n_concepts"]
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
    assert len(ch["unreadable"]) == 1 and len({r[0] for r in rows}) == N_T
    # v3: у концепции две оболочки или оболочка без базовой записи - отказ
    d = _copy_parts(pdir, os.path.join(tmp, "case"))
    doc = json.load(open(L0.part_path(d, "KV-3", N_T, 3, S_T), encoding="utf-8"))
    jc = doc["columns"].index("case_arch")
    doc["rows"][0][jc] = "оболочка"
    json.dump(doc, open(L0.part_path(d, "KV-3", N_T, 3, S_T), "w", encoding="utf-8"))
    with pytest.raises(M.MergeError, match="index_errors"):
        M.merge_parts(M.find_parts([d], "KV-3", N_T, S_T), "KV-3", N_T, S_T)


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
                 s_make=s[3], score=0.0, reject="" if feas else "soft_G", trl_min=int(rng.integers(3, 10)), case_arch="картридж",
                 case_pick=1)
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
                    target=None, detail_top=0, cycle_k=0, closure_proof=False)
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


# =========================================================================== исправления скептиков и выборка v2
def test_fastgas_matches_gtu_explore():
    """скалярные свойства газа FastGas - те же таблицы gtu_explore и та же интерполяция (до округления)"""
    E = FS._E(); ex = E._ex
    rng = np.random.default_rng(3)
    for T, b in zip(rng.uniform(150, 2500, 300), rng.uniform(0, 1, 300)):
        assert math.isclose(E.h(T, b), ex.h(T, b), rel_tol=1e-12, abs_tol=1e-6)
        assert math.isclose(E.s0(T, b), ex.s0(T, b), rel_tol=1e-12, abs_tol=1e-9)
        hv, sv = ex.h(T, b), ex.s0(T, b)
        assert abs(E.T_h(hv, b) - ex.T_h(hv, b)) < 1e-6 and abs(E.T_s(sv, b) - ex.T_s(sv, b)) < 1e-6
    assert math.isclose(E.far_for(700.0, 1400.0), ex.far_for(700.0, 1400.0), rel_tol=1e-12)


def test_violation_never_negative_and_rounding():
    """φ 0,019995 при пределе 0,02-0,16: значение не округляется до предела, превышение > 0 (раньше было -0,875)"""
    r = FS._rule("phi2", False, 0.019995, [0.02, 0.16])
    assert r["val"] == 0.019995 and FS.violation([r]) > 0
    for val in (0.0199, 0.17, 0.5, 1e-5):
        assert FS.violation([FS._rule("phi1", False, val, [0.02, 0.16])]) > 0
    assert FS.violation([FS._rule("soft_G", False, 30.0, 25.0)]) == round(5 / 25, 4)


def test_radial_turbine_polytropic_equivalent():
    """P1: закон U/C0 - изоэнтропический (gtu_explore.turb_radial), _core расширяет политропно: η_p < η_is при ER > 1,
    и политропное расширение с η_p даёт ту же работу, что изоэнтропическое с η_is"""
    E = FS._E()
    t = dict(ER=3.5, T_in=1300.0, b=0.25)
    dhs = FS._stage_dhs(t, 1)
    U = 0.70 * math.sqrt(2 * dhs)
    eta_p, nu, dh, eta_is = FS._radial_stage(t, 1, U)
    assert abs(nu - 0.70) < 1e-9 and abs(eta_is - 0.86) < 1e-9 and eta_p < eta_is
    T2p = E.T_s(E.s0(1300.0, 0.25) - eta_p * E.Rg(0.25) * math.log(3.5), 0.25)
    assert math.isclose(E.h(1300.0, 0.25) - E.h(T2p, 0.25), eta_is * dhs, rel_tol=1e-6)


def test_recuperator_energy_balance_and_density():
    """P4: выхлоп за рекуператором - из баланса энтальпий; P3: кажущаяся плотность модуля 1777 кг/м3 (стенка 0,15 мм, 1500 м2/м3)"""
    E = FS._E()
    for rotors in (1, 2):
        c = dict(cool="N", m_tw="литьё", cores=1, recup="печатный", comb="кольцевая", pgain=0.0, rotors=rotors, t41=1050.0, eps_rec=0.85)
        r = FS._core(c, FS._comp([(4.5, 0.86)]), 0.85, 0.87)
        xn, xr = FS.cool_air(c); mc = 1.0 - xn - xr - 0.015
        m_g, b_g = (r["t"]["m"] + xr, r["t"]["b"] * r["t"]["m"] / (r["t"]["m"] + xr)) if rotors == 1 else (r["pt"]["m"], r["pt"]["b"])
        q_gas = m_g * (E.h(r["T5"], b_g) - E.h(r["Tstack"], b_g))
        q_air = mc * (E.h(r["T3"]) - E.h(r["T2"]))
        assert math.isclose(q_gas, q_air, rel_tol=1e-6), (rotors, q_gas, q_air)
    i = next(i for i in range(64) if FS.concept(i, 64)["recup"] == "печатный")
    c = FS.concept(i, 64); t = FS.target("KV-3"); p = FS.quick_perf(c, t)
    g = FS.derive(c, t, p)
    assert math.isclose(g["rec"]["mass_kg"] / g["rec"]["V_m3"], 7900 * 0.00015 * FS.REC_BETA, rel_tol=1e-9)


def test_an2_rules_present():
    """P2: AN² последней ступени турбины и СТ не выше проверенной ступени frame3 (6,03e7); код правила без цифры на конце"""
    assert FS.AN2_MAX == 6.03e7
    for code in ("turb_ansq", "pt_ansq"):
        assert code.rstrip("0123456789") == code
    n = 64
    seen = set()
    for i in range(n):
        c = FS.concept(i, n); t = FS.target("KV-3")
        sol = L0.solve(c, t)
        if sol["g"] is None:
            continue
        g = sol["g"]
        rules = {r["code"]: r for r in FS.check(c, t, sol["perf"], g)}
        if c["turb"] == "О":
            assert rules["turb_ansq"]["ok"] == (g["tg"]["AN2"] <= FS.AN2_MAX); seen.add("turb")
        if g["pt"]:
            assert rules["pt_ansq"]["ok"] == (g["pt"]["AN2"] <= FS.AN2_MAX); seen.add("pt")
    assert seen == {"turb", "pt"}


def test_sampler_v2_windows_and_repair():
    """выборка v2: точка (π, T41) - внутри окна, номинальный запас там >= 0; ремонт - только в шаблоне архетипа и замыкает
    точку; запись не зависит от кэшей и порядка счёта"""
    n = 48
    t = FS.target("KV-3")
    for i in range(n):
        c = FS.concept(i, n, "KV-3")
        T = FS.ARCH[c["arch"]]["allowed"]
        for ax in FS.REPAIR_AXES:
            if ax in T and not (ax == "rotors" and c["cores"] > 1):
                assert c[ax] in T[ax], (i, ax, c[ax])
        mg = FS.MG0 if "без запасов" in c["repair"] else FS.MG     # окно без запасов - запись в repair
        if "T41/π" not in c["dflags"]:
            m, _ = FS._margins(c, t, c["pi"], c["t41"], mg)
            assert m >= -0.02, (i, m)
        if "исходная не замыкает" in c["repair"]:
            assert FS.closure(c)["ok"]
        if c["rotor_join"] in ("заодно", "сварка"):
            assert c["mat_tw"] != "SC"
    a = [FS.concept(i, n, "KV-3") for i in range(6)]
    FS._CLOSURE.clear(); FS._CANDS.clear(); FS._NOCLOSE.clear(); FS._DM.clear()
    b = [FS.concept(i, n, "KV-3") for i in reversed(range(6))][::-1]
    assert a == b


def test_nominal_attainable_and_closure_kv6():
    """v2.1 (kv6-диагностика): номинал окна не лучше достижимого - КПД СТ не выше ETA_NOM["pt"] и падает, когда AN² на частоте
    пакета требует ca/U 0,8; нагрузка осевой турбины на валу пакета не ниже w/(nt·U_max²); замыкание - без запасов окон;
    лучший вход - только из шаблона; кластер KV-6 с 2Ц не замыкает точку (общая СТ, AN²), с ОЦ и охлаждаемыми SC - замыкает"""
    t6 = FS.target("KV-6")
    base = dict(arch="безмасляный", kv="KV-6", rotors=2, cores=1, comp="2Ц", n_ax=0, inlet="одно", turb="О", cool="K", m_tw="литьё",
                mat_tw="SC", recup="нет", comb="кольцевая", rotor_join="сборный", m_cw="поковка", mat_cw="Ti64", bear="магнит")
    cc = dict(FS.best_case(base, "mid"), pi=14.0, t41=1200.0)
    nm = FS._nominal(cc, t6)
    assert nm["ok"] and nm["att"]["ep"] <= FS.ETA_NOM["pt"] + 1e-12
    assert nm["att"]["ep"] < FS.ETA_NOM["pt"] - 0.005                 # KV-6: ca/U СТ у предела 0,8 - потеря выходной скорости
    one = dict(cc, rotors=1, arch="гибрид", nt=2)
    nm1 = FS._nominal(one, t6)
    U = min(FS.MK.u_limit_turb("О", "SC", 1200.0)[0], FS.u_cap_size(one, FS.d_cap_turb(one)) * (1 - FS.M_GEO))
    assert nm1["att"]["lt"] >= nm1["r"]["t"]["w"] / (FS._nt_max(one) * U ** 2) - 1e-9
    assert FS.closure(base)["margins"] == list(FS.CLOSURE_MG) == [0.0, 0.0, 0.0]
    assert FS.best_case(dict(base, arch="традиционная"))["inlet"] == "двух"
    cl2 = dict(base, arch="кластер", cores=2, bear="масло")
    r2 = FS.closure(cl2)
    assert not r2["ok"] and r2["why"] == "AN2"
    assert FS.closure(dict(cl2, cores=4))["best"]["m"] == r2["best"]["m"]   # потери кластера в модели не зависят от числа ядер
    assert FS.closure(dict(cl2, comp="ОЦ", n_ax=3))["ok"]


def test_closure_cache_key_turbine_alloy():
    """ключ кэша замыкания различает сплав колеса турбины: номинал v2.1 (_attain) зависит от предела окружной скорости
    радиальной турбины (dT_cap сплава), а у K418/CM247LC (литьё) и IN718/Haynes282/CM247LC (печать) с охлаждением K предел T41
    один (1050 °C). Без сплава в ключе итог замыкания зависел от порядка счёта (KV-3: RDC, печать, Р, K, 2 вала, 2Ц)"""
    import itertools
    t = FS.target("KV-3")
    base = dict(arch="RDC", kv="KV-3", rotors=2, cores=1, comp="2Ц", n_ax=0, inlet="одно", turb="Р", nt=1, npt=2, cool="K", recup="нет",
                comb="RDC", bear="масло", m_tw="печать", rotor_join="сборный", m_cw="поковка", mat_cw="Ti64")
    cs = [dict(base, mat_tw=m) for m in ("IN718", "Haynes282", "CM247LC")]
    assert len({FS._sig(c) for c in cs}) == 3
    clean = {}
    for c in cs:
        FS._CLOSURE.clear()
        clean[c["mat_tw"]] = FS.closure(c)["ok"]
    for perm in itertools.permutations(cs):
        FS._CLOSURE.clear()
        assert {c["mat_tw"]: FS.closure(c)["ok"] for c in perm} == clean, [c["mat_tw"] for c in perm]
    FS._CLOSURE.clear()


# =========================================================================== v3: решения владельца 15.11.2026
def test_trl_multiplier_one_function():
    """п. 1: итоговая оценка = взвешенная × (0,5 + 0,5·TRL/9) - одна функция final_score для словаря баллов и матрицы, запись
    уровня 0 и слияние (scores) считают ею; чувствительность к весам - тоже с множителем"""
    assert L0.trl_factor(9) == 1.0 and math.isclose(L0.trl_factor(3), 2 / 3) and math.isclose(L0.trl_factor(4.5), 0.75)
    assert np.allclose(L0.trl_factor(np.array([9, 5, 3])), [1.0, 0.5 + 0.5 * 5 / 9, 0.5 + 0.5 * 3 / 9])
    s = dict(parts=6.0, plugplay=5.0, service=4.0, effmake=7.0)
    assert math.isclose(L0.final_score(s, 9), L0.weighted(s)) and math.isclose(L0.final_score(s, 3), L0.weighted(s) * 2 / 3)
    S = np.array([[s[k] for k in L0.CRIT]] * 3)
    assert np.allclose(L0.final_score(S, np.array([9, 5, 3])), [L0.final_score(s, t) for t in (9, 5, 3)])
    r = L0.evaluate_l0(2, 64, "KV-3")
    sd = dict(parts=r["s_parts"], plugplay=r["s_plugplay"], service=r["s_service"], effmake=r["s_effmake"])
    assert abs(r["score"] - L0.final_score(sd, r["trl_min"])) < 2e-3 and abs(r["score_w"] - L0.weighted(sd)) < 2e-3
    assert math.isclose(r["k_trl"], L0.trl_factor(r["trl_min"]), abs_tol=1e-4)
    # слияние: та же функция; при равных баллах концепция с меньшим TRL ниже
    T = synth_table(n=80, n_niche=10, p_feas=1.0)
    sc = M.scores(T, L0.WEIGHTS)
    assert np.allclose(sc, L0.final_score(T.S, T.trl, L0.WEIGHTS))
    T.S[:] = 5.0
    sc = M.scores(T, L0.WEIGHTS)
    for a in range(10):
        for b in range(10):
            if T.trl[a] < T.trl[b]:
                assert sc[a] < sc[b]
    # чувствительность к весам: и сценарии, и случайные наборы - final_score
    import inspect
    src = inspect.getsource(M.sensitivity)
    assert "L0.final_score" in src and "T.S[fpos] @ M.T" not in src
    src1 = inspect.getsource(__import__("funnel_l1").evaluate_point) + inspect.getsource(__import__("funnel_l1").add_vigv)
    assert "L0.final_score" in src1 and "score = L0.weighted" not in src1


def test_case_axis_variants_kv3():
    """п. 3: у каждой концепции - базовая запись (картридж; у традиционной - разъём) и оболочка, если она входит в камеру LPBF;
    ровно один выбор концепции; оболочка - печатный корпус и цельная печатная горячая часть; ниша включает ось «корпус»; цикл тот же"""
    n, seen = 24, set()
    for i in range(n):
        rows = L0.evaluate_concept(i, n, "KV-3")
        assert rows[0]["case_arch"] == FS.case_base(FS.concept(i, n, "KV-3"))
        assert sum(r["case_pick"] for r in rows) == 1 and len(rows) in (1, 2)
        for r in rows:
            assert r["niche"].endswith("|" + r["case_arch"]) and L0.niche_v2(r["niche"]).count("|") == 6
            seen.add(r["case_arch"])
        if len(rows) == 2:
            b, sh = rows
            assert sh["case_arch"] == "оболочка" and sh["m_case"] == "печать" and sh["cartridge"] == 0
            if sh["cool"] != "C":
                assert sh["hot_integral"] == 1 and sh["m_hot"] == "печать"
            assert (b["eta_th_shaft"], b["G"], b["EGT"], b["pi"], b["t41"]) == (sh["eta_th_shaft"], sh["G"], sh["EGT"], sh["pi"], sh["t41"])
            assert "печать_оболочка" in sh["techs"] and sh["trl_min"] <= FS.MK.TRL["печать_оболочка"]["trl"]
            assert "shell_fit" not in (sh["reject"] or "")
        else:
            assert rows[0]["shell_note"]
        assert rows[0]["case_arch"] == "разъём" or rows[0]["cartridge"] == 1
    assert {"картридж", "оболочка"} <= seen


def test_shell_parts_live_geometry_and_trap_rule():
    """п. 3, правило собираемости: REF_BOM (традиционный KV-3) как оболочка - 211 деталей (−108), крепёж 340 при болтах (D 1,84 м >
    1,2 м); добавлены передняя крышка, кольца СА 2 и СА 3 (organic/print_review: shell_assembly). Внешний печатный корпус картриджа
    не запирает колёса; ротор одной деталью - элементы между колёсами из половин"""
    ref = PM.count(PM.REF_BOM)
    assert ref["parts"] == 319 and ref["trap"]["add_cold"] == ref["trap"]["add_hot"] == 0 and PM.N_HALF == ref["n_eff"]
    sh = dict(PM.REF_BOM, case=dict(method="печать", D=PM.REF_BOM["case"]["D"], pieces=1), hot_integral=True,
              comb=dict(kind="кольцевая", method="печать", n_inj=14), case_arch="оболочка", conn="быстроразъём")
    r = PM.count(sh)
    assert r["parts"] == 211 and ref["parts"] - r["parts"] == 108 and r["fasteners"] == 340
    assert (r["trap"]["n_trap"], r["trap"]["add_cold"], r["trap"]["add_hot"]) == (6, 1, 2)
    small = PM.count(dict(sh, case=dict(sh["case"], D=1.1)))                       # фланец под хомут: +2 вместо болтов
    assert small["fasteners"] == ref["by_unit"]["ротор"]["fasteners"] + 3 * 2
    cart = PM.count(dict(sh, case_arch="картридж"))                                 # печатный внешний корпус картриджа
    assert cart["trap"]["add_cold"] == 2 and cart["trap"]["cover_cold"] == 0
    mono = PM.count(dict(sh, rotor_join="заодно", turb=dict(type="О", stages=[dict(blisk=True, zt=55, zn=40, ngv_ring=True)] * 3)))
    assert mono["trap"]["add_cold"] == 2 and mono["trap"]["add_hot"] == 1 + 2 * 2 and mono["trap"]["split"] == 1 + 3


def test_shell_make_hip_and_control():
    """печать оболочки: камера принтера, ГИП (серийная камера 0,88 м; Ø1,6 и Ø2,05 м - две в мире), контроль ускорителем 9 МэВ"""
    MK = FS.MK
    base = dict(name="оболочка (холодная часть)", kind="static", method="печать", material="Ti64", mass_kg=150.0, rotating=False, shell=True)
    a = MK.assess_part(dict(base, D=0.7, L=0.5))
    assert a["hip"] and a["hip_vessel"] == 0.88 and a["proc_q"] <= 0.6 and a["fit"]
    b = MK.assess_part(dict(base, D=1.34, L=1.08))
    assert b["hip"] and b["hip_vessel"] == 1.6 and b["weeks"] > a["weeks"] + MK.WEEKS["hip_unique"] - 1e-9
    c = MK.assess_part(dict(base, D=2.02, L=1.0))
    assert not c["hip"] and c["hip_vessel"] is None and c["proc_q"] == 0.2
    plain = MK.assess_part(dict(base, D=0.7, L=0.5, shell=False))
    assert not plain["hip"] and plain["proc_q"] == 1.0


def test_service_shell_engine_swap_repair_site():
    """оболочка: на объекте меняется двигатель целиком, ремонт - на заводе; картридж - газогенератор вынимается целиком"""
    i = next(i for i in range(64) if len(L0.evaluate_concept(i, 64, "KV-3")) == 2 and FS.concept(i, 64, "KV-3")["cores"] == 1)
    c0 = FS.concept(i, 64, "KV-3"); t = FS.target("KV-3")
    sol = L0.solve(c0, t)
    cs = FS.apply_case(c0, "оболочка"); gs = FS.derive(cs, t, sol["perf"])
    _, sv, pp, _ = FS.model_inputs(cs, t, sol["perf"], gs)
    eng = next(m for m in sv["modules"] if m["name"] == "двигатель в оболочке")
    assert eng["frequent"] and eng["repair"] == FS.SM.REPAIR_FACTORY and not any(m["name"] in ("ротор", "горячая часть") for m in sv["modules"])
    sr = FS.SM.score(sv)
    assert sr["repair_sites"]["двигатель в оболочке"] == "на заводе" and pp["module_kg_max"] >= eng["mass_kg"]
    cb = FS.apply_case(c0, FS.case_base(c0))
    _, svb, _, _ = FS.model_inputs(cb, t, sol["perf"], sol["g"])
    names = {m["name"]: m for m in svb["modules"]}
    assert "картридж газогенератора" in names and names["картридж газогенератора"]["repair"] == FS.SM.REPAIR_SWAP


def test_merge_case_axis_and_trl_blocks(parts, tmp):
    """слияние v3: доли и сходимость - по концепциям; выбор концепций (оболочка / картридж) в сумме = допустимым концепциям;
    блок множителя готовности"""
    pdir, _ = parts
    r = M.run_merge(M.find_parts([pdir], "KV-3", N_T, S_T), "KV-3", N_T, S_T, os.path.join(tmp, "mc"), k=1, top=3, k_near=0,
                    target=None, detail_top=0, cycle_k=0, closure_proof=False)
    res = r["res"]
    ca = res["case_axis"]
    assert res["summary"]["n"] == N_T and ca["concepts"] == N_T and sum(ca["picked_among_feasible"].values()) == res["summary"]["feasible"]
    assert ca["shell_possible"] + ca["shell_impossible"] == N_T
    assert res["niches"]["n_feasible_v2"] <= res["niches"]["n_feasible"]
    assert res["trl_multiplier"]["top_overlap_with_no_multiplier"] <= res["trl_multiplier"]["top_n"]
    assert all(e.get("case_arch") in FS.CASE_VALUES for e in r["sel"]["selected"])


def test_b2_shift_tries_original_inlet_first(monkeypatch):
    """правка B2 с доработкой проверки: при пустом окне T41/π сдвиг к лучшим осям сначала при исходном входе (вход не меняется,
    если окно нашлось), лучший вход шаблона - только если при исходном окна нет (запись «inlet: ... (окно T41/π)»)"""
    n, kv = 64, "KV-3"
    tgt = FS.target(kv)
    M_ = FS.design_matrix(n)
    i = next(i for i in range(n) if FS._discrete(dict(zip(FS.AXES, M_[i])), FS._pick(M_[i][0], FS.ARCH_KEYS))["inlet"] == "двух"
             and "одно" in FS.ARCH[FS._pick(M_[i][0], FS.ARCH_KEYS)]["allowed"].get("inlet", FS._DEF["inlet"]))
    u = dict(zip(FS.AXES, M_[i]))
    arch = FS._pick(u["arch"], FS.ARCH_KEYS)

    def pre():
        c0 = FS._discrete(u, arch); c0["kv"] = kv
        c, rep = FS._repair(u, arch, c0)
        c.update(idx=i, n=n, seed=FS.SEED, kv=kv)
        c["repair"] = rep
        return c

    real = FS._t41_pi_window
    for good in ("двух", "одно"):
        calls = []

        def fake(c, tg, tl, mg=FS.MG, good=good, calls=calls):
            calls.append(c["inlet"])
            if len(calls) == 1 or c["inlet"] != good:
                return None, dict(m=-0.1, t41=tl, pi=5.0, worst="w")
            return (tl - 30.0, tl), dict(m=0.01, t41=tl, pi=5.0, worst="w")     # окно нашлось (дальше design - настоящий)
        c = pre()
        monkeypatch.setattr(FS, "_t41_pi_window", fake)
        monkeypatch.setattr(FS, "closure", lambda c, mg=None: dict(ok=True, mu="mid"))
        assert c["inlet"] == "двух"
        FS.design(c, u, tgt)
        monkeypatch.undo()
        if good == "двух":
            assert calls[:2] == ["двух", "двух"] and "(окно T41/π)" not in c["repair"]   # вход дальше может менять только φ и частота
        else:
            assert calls[:4] == ["двух", "двух", "двух", "одно"] and "inlet: двух -> одно (окно T41/π)" in c["repair"]
    FS._CLOSURE.clear()
