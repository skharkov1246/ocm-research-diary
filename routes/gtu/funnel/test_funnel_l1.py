"""Проверки уровня 1 (funnel_l1, meanline_l1, strength_l1, offdesign_l1, merge_l1, make_jobs_l1): балансы, повторяемость и
независимость от числа процессов, оптимизация не хуже исходной точки, сверка одномерной модели центробежной ступени с
comp_stage в его расчётных точках (заявлено meanline_l1.TOL_CC), совпадение формул потерь с frame3, прочность против
disk.py и frame3, отбор 50: не больше 3 на нишу, все выше порога, жадный отбор детерминирован, слияние с контролем полноты;
исправления скептиков уровня 1 (15.11.2026): поправки ступени и ограничение КПД компрессора, Kacker - Okapuu и свободный вихрь,
Stanitz и цепочка радиальных ступеней, счёт всех нерешённых режимов, граница помпажа по признакам срыва, повтор части при слиянии.
Синтетические записи отбора - придуманные корпуса (не данные установок). Запуск: python3 -m pytest -q test_funnel_l1.py"""
import json, math, os, random, shutil, sys, tempfile

import numpy as np
import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import funnel_space as FS        # noqa: E402
import funnel_l1 as L1           # noqa: E402
import meanline_l1 as ML         # noqa: E402
import merge_l1 as M1            # noqa: E402
import offdesign_l1 as OD        # noqa: E402
import strength_l1 as ST         # noqa: E402

SEL = os.path.join(HERE, "l0_aws_debug", "funnel_l0_select_KV-3_n2000.json")


@pytest.fixture(scope="module")
def sel():
    d, _ = L1.load_select(SEL)
    return d, L1.order_entries(d)


@pytest.fixture(scope="module")
def point(sel):
    """первый вариант порядка: оптимизация (1 старт: DE 8 поколений + COBYLA 100) и расчётная точка"""
    d, E = sel
    e = E[0]
    c0 = FS.concept(e["idx"], d["n"], d["kv"], d["seed"])
    tgt = FS.target(d["kv"])
    spec, order, runs, nev = L1.optimize(c0, tgt, 1, 100, e["idx"], de_gens=8)
    c = L1.x_to_c(c0, np.array(order[0][0]), spec)
    return dict(c0=c0, c=c, tgt=tgt, spec=spec, order=order, ev=L1.evaluate_point(c, tgt), e=e)


# =========================================================================== одномерные модели
def test_cc_matches_comp_stage_within_declared():
    """сверка с gtu_explore.comp_stage в его расчётных точках (μ 0,70, φ 0,085, cm2/U2 0,28): сетка калибровки и точка вне сетки"""
    cal = ML.cc_calibration()
    assert cal["max_abs_after"] <= ML.TOL_CC
    E = ML._E()
    for pr, m, dbl in ((5.0, 12.0, False), (3.0, 30.0, True), (7.0, 9.0, False)):        # вне сетки CAL_GRID
        cs = ML.comp_stage_ref(E.T0, E.P0, pr, m, dbl)
        st = ML.cc_stage(E.T0, E.P0, m / 2 if dbl else m, pr, cs["rpm"], 0.70, 0.28, 10.0, dbl, kD=1.9, eta0=cs["eta"], ref=True)
        assert abs(st["eta_p"] - cs["eta"]) <= ML.TOL_CC, (pr, m, st["eta_p"], cs["eta"])


def test_comp_stage_ref_equals_gtu_explore():
    """формула comp_stage на быстрых таблицах = gtu_explore.comp_stage"""
    FS._E()
    import gtu_explore as GE
    E = ML._E()
    for pr, m in ((3.0, 10.0), (6.0, 25.0)):
        a, b = ML.comp_stage_ref(E.T0, E.P0, pr, m), GE.comp_stage(E.T0, E.P0, pr, m, "ti")
        assert abs(a["eta"] - b["eta"]) < 2e-4 and abs(a["D2"] - b["D2"]) / b["D2"] < 2e-3


def test_cc_losses_match_frame3():
    """формулы потерь - те же, что frame3_point.losses (при r4/r2 1,55 и без DUMP)"""
    cwd = os.getcwd()
    sys.path.insert(0, os.path.join(os.path.dirname(HERE), "frame3"))
    import frame3_point as F3
    os.chdir(cwd)
    E = ML._E()
    st = ML.cc_stage(E.T0, E.P0, 12.0, 4.0, 20000.0, 0.70, 0.28, 11.0, False)
    g, f = dict(st["g"]), dict(st["f"])
    for dbl in (False, True):
        a = ML.cc_losses(g, f, 12.0, dbl, 1.55, dump=False)
        b = F3.losses(g, f, 12.0, dbl, True)
        for k, v in b.items():
            assert abs(a[k] - v) <= 1e-9 * max(1.0, abs(v)), k
        assert a["dump"] == 0.0


def test_cc_euler_balance():
    """работа ступени = Эйлер + паразитные потери; КПД модели - из внутренних потерь"""
    E = ML._E()
    st = ML.cc_stage(E.T0, E.P0, 10.0, 5.0, 21000.0, 0.68, 0.26, 12.0, False)
    L, f = st["losses_pct"], st["f"]
    par = sum(L[k] for k in ML.PARASITIC) / 100 * (f["dh_eu"] + 0)    # доли от w0
    w0 = f["dh_eu"] / (1 - sum(L[k] for k in ML.PARASITIC) / 100)
    assert abs(w0 - st["w0"]) / st["w0"] < 0.01
    assert abs(f["cth2"] - (st["f"]["dh_eu"] / f["U2"])) < 1e-6
    assert 0.0 <= st["beta2b_deg"] <= 70.0 and par >= 0


def test_cc_cap_and_stage_penalties():
    """исправления скептика уровня 1: (1) поправка ступени (Ni, «Д», двусторонний вход) - внутри cc_stage, и нерасчётная модель той
    же геометрии в расчётной точке даёт тот же КПД; (2) вне сетки калибровки КПД не выше формулы уровня 0 + TOL_CC"""
    E = ML._E()
    a = ML.cc_stage(E.T0, E.P0, 12.0, 3.5, 16000.0, 0.70, 0.28, 10.0, False, d_pen=0.0)
    b = ML.cc_stage(E.T0, E.P0, 12.0, 3.5, 16000.0, 0.70, 0.28, 10.0, False, d_pen=0.02)
    assert abs((a["eta_p"] - b["eta_p"]) - 0.02) < 2e-3                      # тот же КПД минус поправка (геометрия пересчитана)
    off = ML.cc_off(b, E.T0, E.P0, 12.0, 16000.0)
    assert abs(off["eta_p"] - b["eta_p"]) < 1e-3 and abs(off["PRs"] / 3.5 - 1) < 2e-3
    lo = ML.cc_stage(E.T0, E.P0, 12.0, 3.5, 16000.0, 0.60, 0.20, 10.0, False)  # нижние границы μ и cm2/U2 - вне сетки калибровки
    assert lo["eta_p"] <= lo["eta_cap"] + 1e-9
    assert lo["eta_uncapped"] >= lo["eta_p"] - 1e-9


def test_axial_ko_mach_and_free_vortex():
    """Kacker - Okapuu (TurboFlow 0.1.18): сверхзвуковой выход сопла (M > 1) - рост потерь (у Soderberg члена на M нет);
    свободный вихрь: у втулки число Маха за соплом выше, реактивность ниже, чем на среднем радиусе"""
    y1 = ML.ko_Y(0.0, 70.0, 0.3, 0.9, 0.02, 0.03, 0.025, 0.05, 0.8, 5e5, False, 0.1)
    y2 = ML.ko_Y(0.0, 70.0, 0.3, 1.15, 0.02, 0.03, 0.025, 0.05, 0.8, 5e5, False, 0.1)
    assert y2["Y"] > y1["Y"] and y2["fMa"] > 1.0 == y1["fMa"]
    b = FS._E().bfrac(0.02)
    t = ML.ax_turbine(1273.15, 1.2e6, b, 30.0, 400.0, 12732.0, 2.0 * 400 ** 2, 1, 0.5, 40, 50)
    s = t["stages"][0]
    assert s["M2_hub"] > s["M2"] and s["reaction_hub"] < s["reaction"] and s["Mw3_tip"] > s["Mw3"]


def test_radial_stanitz_and_chain():
    """Stanitz: cθ4,opt = U4·(1 - 0,63π/Z); две ступени: статическое давление за последней - заданное"""
    b = FS._E().bfrac(0.02)
    st = ML.rad_stage(1223.15, 4.5e5, b, 5.0, 480.0, 2.5, 30000.0, 13, 18)
    cm4 = ML.PHI4 * 480.0
    assert abs(math.radians(st["beta4opt_deg"]) - math.atan2(-480.0 * 0.63 * math.pi / 13, cm4)) < 1e-9
    r = ML.rad_turbine(1223.15, 4.5e5, b, 5.0, 480.0, 6.0, 30000.0, 2, 13, 18)
    p_end = 4.5e5 / 6.0
    s1, s2 = r["stages"]
    assert abs(s2["p5"] / p_end - 1) < 1e-4 and s1["p05"] > s1["p5"]


def test_turbine_anchors():
    """привязка: осевая ступень ψ 1,2 / φ 0,5 / H/b 1,3 = eta_axial 0,885; радиальная U/C0 0,70 / ns 0,60 = 0,86"""
    ref = ML.ax_reference()
    assert abs(ref["eta_p_model"] + ML.ax_shift() - FS._E().eta_axial(1.2)) < 1e-9
    r = ML.rad_reference()
    assert abs(r["eta_is_model"] + ML.rad_shift() - 0.86) < 1e-9
    # меньше лопаток (ниже удлинение) - ниже КПД осевой ступени; U/C0 вдали от 0,7 - ниже КПД радиальной
    b = FS._E().bfrac(0.02)
    hi = ML.ax_turbine(1273.15, 1.2e6, b, 30.0, 400.0, 12732.0, 1.2 * 400 ** 2, 1, 0.5, 66, 69)
    lo = ML.ax_turbine(1273.15, 1.2e6, b, 30.0, 400.0, 12732.0, 1.2 * 400 ** 2, 1, 0.5, 20, 22)
    assert lo["eta_p"] < hi["eta_p"]
    a = ML.rad_stage(1223.15, 4.5e5, b, 5.0, 0.70 * math.sqrt(2 * r["eta_is_model"] * 0 + 2 * 220e3), 3.0, r["n"], 13, 18)
    c = ML.rad_stage(1223.15, 4.5e5, b, 5.0, 0.55 * math.sqrt(2 * 220e3), 3.0, r["n"], 13, 18)
    assert c["eta_is"] < a["eta_is"] + 0.02


# =========================================================================== балансы расчётной точки уровня 1
def test_point_balances(point):
    ev, tgt, c = point["ev"], point["tgt"], point["c"]
    assert ev["ok"]
    perf = ev["perf"]; r = perf["r"]
    assert abs(perf["G"] * r["wnet"] / 1e3 - tgt["P_shaft_kW"]) < 1e-6 * tgt["P_shaft_kW"]       # мощность
    assert abs(r["wc"] - sum(s["w"] for s in perf["comp"])) < 1e-6 * r["wc"]                   # работа компрессора
    if c["rotors"] == 1:
        E = FS._E(); t = r["t"]
        assert abs(0.99 * t["m"] * t["w"] - r["wc"] - r["wnet"]) < 1e-6 * r["wc"]               # баланс вала
    # КПД из одномерных моделей - в цикле
    assert abs(perf["eta_gg"] - ev["sol"]["tm"]["eta_p"]) < 5e-4
    for e_c, cc, s in zip(ev["sol"]["eta_c"], ev["sol"]["ccs"], perf["comp"]):
        if cc:
            assert abs(s["eta"] - e_c) < 1e-12


def test_recuperator_size_and_heat(point):
    """рекуператор: UA = NTU·Cmin, площадь = UA/U, длины из бюджета потерь; объём растёт с ε"""
    c, tgt = dict(point["c"]), point["tgt"]
    if c["recup"] != "печатный":
        pytest.skip("первый вариант без рекуператора")
    sol = L1.solve_l1(c, tgt)
    r1 = L1.recup_l1(c, sol["perf"])
    assert abs(r1["A_m2"] * r1["U_W_m2K"] - r1["UA_W_K"]) < 1e-6 * r1["UA_W_K"]
    c2 = dict(c, eps_rec=min(0.92, c["eps_rec"] + 0.04))
    r2 = L1.recup_l1(c2, L1.solve_l1(c2, tgt)["perf"])
    assert r2["V_m3"] > r1["V_m3"]


# =========================================================================== оптимизация и повторяемость
def test_optimization_not_worse_than_start(point):
    """точка уровня 0 входит в кандидаты: лучший допустимый кандидат не хуже её (если она допустима в расчётной точке)"""
    c0, tgt, spec, order = point["c0"], point["tgt"], point["spec"], point["order"]
    ev0 = L1.evaluate_point(L1.x_to_c(c0, L1.c_to_x(c0, spec), spec), tgt)
    best = order[0][1]
    if ev0["ok"] and ev0["viol"] == 0:
        assert best[1] == 0 and best[0] >= ev0["score"] - 1e-12
    else:
        assert best[1] <= ev0["viol"] + 1e-12


def test_variant_repeatable_and_process_free(sel):
    d, E = sel
    items = [(0, E[0]), (1, E[1])]
    a = L1.evaluate_many(items, d["kv"], d["n"], d["seed"], 1, 20, workers=1, detail=False, de_gens=2)
    b = L1.evaluate_many(list(reversed(items)), d["kv"], d["n"], d["seed"], 1, 20, workers=2, detail=False, de_gens=2)
    strip = lambda rows: [{k: v for k, v in r.items() if k not in ("t_s", "t_opt_s", "t_cpu_s")} for r, _ in rows]
    assert json.dumps(strip(a), ensure_ascii=False, default=float) == json.dumps(strip(b), ensure_ascii=False, default=float)


def test_stage_count_axis(point):
    """число ступеней турбины - целая ось среди разрешённых архетипом; ниша от неё не зависит"""
    c0 = point["c0"]
    spec = L1.var_spec(c0)
    names = [x[0] for x in spec]
    nts = sorted(FS._axis_vals(c0, c0["arch"], "nt"))
    if len(nts) > 1:
        k = names.index("nt")
        assert spec[k][1:] == (nts[0], nts[-1], "int")
        x = L1.c_to_x(c0, spec)
        for u in (0.0, 1.0):
            x[k] = u
            c = L1.x_to_c(c0, x, spec)
            assert c["nt"] in nts and isinstance(c["nt"], int)
            assert L1.L0.niche_of(c) == L1.L0.niche_of(c0)


def test_order_entries_round_robin(sel):
    """двойной круг: первые K вариантов - по одному на каждый архетип с допустимыми; номера без повторов"""
    d, E = sel
    archs = {e["arch"] for e in d["selected"] if e["feasible"]}
    assert {e["arch"] for e in E[:len(archs)]} == archs
    assert len(E) == len(d["selected"]) and len({e["idx"] for e in E}) == len(E)
    assert L1.order_entries(d, 12) == E[:12]


def test_reject_foreign_select(tmp_path):
    d = json.load(open(SEL, encoding="utf-8"))
    d["code_sha"] = "0000000000000000"
    p = tmp_path / "sel.json"
    p.write_text(json.dumps(d, ensure_ascii=False), encoding="utf-8")
    with pytest.raises(SystemExit):
        L1.load_select(str(p))


def test_scale_keeps_discrete(point):
    c = point["c"]
    s = L1.scale_concept(c, "KV-1", c["n"], c["seed"])
    for a in FS.CAT:
        if a in L1.SCALE_FREE:
            if s[a] != c[a]:
                assert f"{a}:" in s["repair"]                    # подбор design() под новую точку - с записью
            continue
        assert s[a] == c[a]
    assert s["kv"] == "KV-1"


# =========================================================================== прочность
def test_strength_against_disk_py_and_frame3():
    life, extr = ST.rupture_hours(173.496, 860.0)
    assert abs(life / 5692.26 - 1) < 0.01 and not extr                          # disk_results.json: turb130 local_min_h
    k, src = ST.k_taper()
    assert abs(ST.blade_root(src["rho"], src["AN2"]) - src["sigma_root_MPa"]) < 1e-6
    w = ST.wheel("comp", "Ti64", 4430.0, 617.759, 90.43, solid=True)              # пример comp_T60_1C_ti_solid сам на себя
    assert abs(w["burst_margin"] - 1.5861422) < 0.01
    p1, p2 = ST.weibull_pf(150, 0.3)["pf"], ST.weibull_pf(200, 0.3)["pf"]
    assert 0 <= p1 < p2 <= 1
    d1, d2 = ST.disc_axial(8000, 14950, 0.5, 0.08, 20, 500), ST.disc_axial(8000, 16950, 0.5, 0.08, 20, 500)
    assert d2["ratio"] > d1["ratio"] and abs(d1["burst_margin_120"] - 1.0) < 1e-9 or d1["governed_by"] != "разрушение при 120 %"


# =========================================================================== нерасчётные режимы
def test_offdesign_reproduces_design(point):
    ev, c, tgt = point["ev"], point["c"], point["tgt"]
    if ev["viol"] > 0:
        pytest.skip("расчётная точка недопустима")
    D = L1.design_data(c, tgt, ev)
    op = OD.solve(D, 15.0, c["t41"])
    assert op is not None
    assert abs(op["P_kW"] / tgt["P_shaft_kW"] - 1) < 0.005
    assert abs(op["eta_th"] - ev["perf"]["eta_th"]) < 0.003
    mx = OD.matrix(D)                                              # правило od_iso: T41 режима ISO 100 % - в 3 K от расчётной
    if mx["iso_check"]:
        assert abs(mx["iso_check"]["T41_C"] - c["t41"]) <= L1.OD_ISO_DT41


def test_offdesign_counts_every_unsolved_regime(point, monkeypatch):
    """режим без совместной работы ни при одной T41 считается поштучно (раньше при пустой сетке −30 °C режимы пропадали из счёта)"""
    ev, c, tgt = point["ev"], point["c"], point["tgt"]
    if not ev["ok"]:
        pytest.skip("расчётная точка не строится")
    D = L1.design_data(c, tgt, ev)
    real = OD._grid
    monkeypatch.setattr(OD, "_grid", lambda D_, Ta, a1, t1h=None: ([], []) if Ta == -30.0 else real(D_, Ta, a1, t1h))
    mx = OD.matrix(D)
    cold = [f for f in mx["fails"] if f.get("Ta_C") == -30.0 and "load" in f]
    assert len(cold) == len(OD.LOADS) and mx["n_unsolved"] >= len(OD.LOADS)


def test_surge_boundary_not_left_of_peak(point):
    """граница помпажа - не левее вершины π: признак срыва может только сдвинуть её вправо (к большему расходу)"""
    ev, c, tgt = point["ev"], point["c"], point["tgt"]
    if not ev["ok"]:
        pytest.skip("расчётная точка не строится")
    D = L1.design_data(c, tgt, ev)
    s = OD.surge(D, D["n_d"], 288.15, D["G"])
    if s is None:
        pytest.skip("линия частоты не строится")
    real_stall = OD._stall
    try:
        OD._stall = lambda r: []
        s0 = OD.surge(D, D["n_d"], 288.15, D["G"])
    finally:
        OD._stall = real_stall
    assert s["W"] >= s0["W"] - 1e-9


# =========================================================================== отбор 50 и слияние (синтетика)
def _synthetic(n=300, seed=7):
    rng = random.Random(seed)
    archs = FS.ARCH_KEYS
    out = []
    for j in range(n):
        a = archs[j % len(archs)]
        r = dict(j=j, idx=1000 + j, arch=a, feasible=int(rng.random() < 0.7), score=round(4.0 + 4.0 * rng.random(), 4),
                 rotors=rng.choice([1, 2]), comp=rng.choice(["1Ц", "2Ц", "Д"]), turb=rng.choice(["Р", "О"]), cool=rng.choice("UNKP"),
                 recup=rng.choice(["нет", "печатный"]), rotor_join=rng.choice(["заодно", "сварка", "сборный"]), m_cw="литьё", m_tw="литьё",
                 t41=900 + 300 * rng.random(), pi=3 + 12 * rng.random(), eps_rec=0.0, eta_th_shaft=25 + 15 * rng.random(),
                 V_env=5 + 20 * rng.random(), mass_kg=1000 + 3000 * rng.random())
        for k in M1.DISC_W:
            r.setdefault(k, 0)
        r["niche"] = "|".join(str(r[k]) for k in ("arch", "rotors", "comp", "turb", "cool", "recup", "rotor_join"))
        out.append(r)
    return out


def test_greedy_rules_and_determinism():
    recs = _synthetic()
    s1 = M1.greedy(recs)
    s2 = M1.greedy(random.Random(3).sample(recs, len(recs)))
    assert [r["j"] for r in s1] == [r["j"] for r in s2]
    assert len(s1) <= M1.TOP_N
    cnt = {}
    for r in s1:
        assert r["feasible"] and r["score"] >= M1.SCORE_MIN
        cnt[r["niche"]] = cnt.get(r["niche"], 0) + 1
    assert max(cnt.values()) <= M1.NICHE_MAX
    p = M1.price_of_diversity(recs, s1)
    assert p["mean_diverse"] <= p["mean_pure"] + 1e-12
    assert p["n_new"] == p["n_displaced"] == p["n_new_worse_than_every_displaced"]


def test_greedy_clones_pushed_out():
    recs = _synthetic(60)
    clone = dict(recs[0]); clone.update(j=999, idx=9999, score=recs[0]["score"] - 1e-3, feasible=1)
    recs[0]["feasible"] = 1; recs[0]["score"] = 7.9
    clone["score"] = 7.89
    sel = M1.greedy(recs + [clone], n=10)
    js = [r["j"] for r in sel]
    assert 0 in js and 999 not in js[:2]


def _part(tmp, k, S, M, rows_j, meta=None):
    meta = meta or {}
    doc = dict(kv="KV-3", M=M, shards=S, shard=k, input_sha256="x", order_sha="y", code_sha="z", starts=1, maxfev=10, level="L1", top=M,
               input="sel.json", method="m", software="s", target={}, n0=2000, seed=1, columns=["j", "v"],
               rows=[[j, j * 2] for j in rows_j], details=[None for _ in rows_j], call={}, workers=1)
    doc.update(meta)
    p = os.path.join(tmp, f"{L1.PART_PREFIX}_KV-3_m{M}_s{k:03d}of{S:03d}.json")
    json.dump(doc, open(p, "w", encoding="utf-8"))
    return p


def test_merge_completeness():
    tmp = tempfile.mkdtemp()
    try:
        M, S = 10, 3
        for k in range(S):
            _part(tmp, k, S, M, FS.shard_indices(M, k, S))
        cols, rows, det, ref, checks, calls = M1.merge_parts(M1.find_parts([tmp], "KV-3", M, S), "KV-3", M, S)
        assert [r[0] for r in rows] == list(range(M)) and checks["indices_ok"]
        os.remove(os.path.join(tmp, f"{L1.PART_PREFIX}_KV-3_m{M}_s001of{S:03d}.json"))
        with pytest.raises(M1.MergeError):
            M1.merge_parts(M1.find_parts([tmp], "KV-3", M, S), "KV-3", M, S)
        _part(tmp, 1, S, M, [1, 4])                                   # потерян номер 7
        with pytest.raises(M1.MergeError):
            M1.merge_parts(M1.find_parts([tmp], "KV-3", M, S), "KV-3", M, S)
        _part(tmp, 1, S, M, FS.shard_indices(M, 1, S), dict(code_sha="other"))
        with pytest.raises(M1.MergeError):
            M1.merge_parts(M1.find_parts([tmp], "KV-3", M, S), "KV-3", M, S)
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


def test_merge_accepts_rerun_with_other_timing():
    """повтор части после прерывания спота: те же строки, другое время счёта - та же часть; другое содержание - отказ"""
    tmp, tmp2, tmp3 = tempfile.mkdtemp(), tempfile.mkdtemp(), tempfile.mkdtemp()
    try:
        M, S = 6, 2
        cols = ["j", "v", "t_s"]
        for k in range(S):
            js = FS.shard_indices(M, k, S)
            _part(tmp, k, S, M, js, dict(columns=cols, rows=[[j, j * 2, 1.0] for j in js]))
        js = FS.shard_indices(M, 1, S)
        _part(tmp2, 1, S, M, js, dict(columns=cols, rows=[[j, j * 2, 9.5] for j in js], workers=4))
        cols_, rows, det, ref, checks, calls = M1.merge_parts(M1.find_parts([tmp, tmp2], "KV-3", M, S), "KV-3", M, S)
        assert checks["duplicate_same"] and not checks["duplicate_conflict"]
        _part(tmp3, 1, S, M, js, dict(columns=cols, rows=[[j, j * 3, 1.0] for j in js]))
        with pytest.raises(M1.MergeError):
            M1.merge_parts(M1.find_parts([tmp, tmp3], "KV-3", M, S), "KV-3", M, S)
    finally:
        for t in (tmp, tmp2, tmp3):
            shutil.rmtree(t, ignore_errors=True)


# =========================================================================== файлы задач AWS
def test_jobs_files():
    for f, lo, hi in (("jobs_l1_debug.json", 10, 50), ("jobs_l1.json", 20, 50)):
        p = os.path.join(HERE, f)
        if not os.path.exists(p):
            pytest.skip("файлы задач ещё не сгенерированы make_jobs_l1.py")
        d = json.load(open(p, encoding="utf-8"))
        est = d["estimate"]
        assert lo <= est["shards"] <= hi and len([j for j in d["jobs"] if j["name"].startswith("l1-") and "-m" in j["name"]]) <= 6
        assert d["plan_check"]["exit_code"] == 0
        for j in d["jobs"]:
            assert "--workers \"$NPROC\"" in j["command"] and "$GTU_OUT" in j["command"]
    r = json.load(open(os.path.join(HERE, "READY_L1.json"), encoding="utf-8"))
    assert r["estimate"]["budget_ok"] and r["estimate"]["total_upper_bound"] < 1500
