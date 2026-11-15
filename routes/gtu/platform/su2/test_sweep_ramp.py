"""Тесты sweep_ramp.py и make_jobs_su2_kv3_sweep.py без SU2 и без AWS. Все истории и числа придуманы (SYNTHETIC):
похожи на историю SU2 по набору столбцов, но не взяты из расчётов."""
import csv
import json
import math
import os
import sys

import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, "..", "aws"))
sys.dont_write_bytecode = True

import run_case as RC  # noqa: E402
import sweep_ramp as SW  # noqa: E402
import make_jobs_su2_kv3_sweep as MS  # noqa: E402

COLS = ["Inner_Iter", "rms[Rho]", "rms[k]", "max[Rho]", "Avg_Massflow(inlet)", "Avg_Massflow(outlet)",
        "Avg_TotalPress(inlet)", "Avg_TotalPress(outlet)", "Avg_TotalTemp(inlet)", "Avg_TotalTemp(outlet)", "PRTT",
        "TotTotEff", "Avg CFL", "PressureOut_1", "MassFlowOut_1"]


def write_history(path, n, *, m_in=1.0, imb=0.005, rho0=-4.0, rho1=-7.5, pr=3.2, tr=1.55, backflow=False, drift=0.0,
                  nan_last=False, p_meas=210000.0):
    """SYNTHETIC история: расход m_in (с дрейфом drift на всю длину), дисбаланс imb, rms[Rho] линейно rho0 -> rho1."""
    eta = (pr ** (0.4 / 1.4) - 1) / (tr - 1)
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow([f'"{c}"' for c in COLS])
        for i in range(n):
            mi = m_in * (1 + drift * i / max(n - 1, 1))
            mo = (mi * (1 - imb)) if backflow else -mi * (1 - imb)
            rho = rho0 + (rho1 - rho0) * i / max(n - 1, 1)
            row = [i, rho, rho + 1, rho + 2, mi, mo, 101325.0, 101325.0 * pr, 288.15, 288.15 * tr, pr, 100 * eta, 12.0,
                   p_meas, -mo]
            if nan_last and i == n - 1:
                row[1] = "nan"
            w.writerow(row)


CTX = {"est": 300000.0, "fs_p": 90750.0, "pr_1d": 4.3, "omega_nom": 1500.0, "deadline": 1e12, "work": "/nonexistent",
       "m_pass_1d": 1.2}


def strat_two_steps(**kw):
    s = {"id": "t_syn", "title": "SYNTHETIC", "transforms": ["hist"], "set": {},
         "steps": [{"name": "a_spin", "frac": 1.0, "load": 0.6, "iter": 300, "ramp": 100, "pramp": 200, "gate": "intermediate",
                    "hold": {"iter": 200, "max": 2}},
                   {"name": "b_base", "frac": 1.0, "load": 0.75, "iter": 500, "pramp": 200, "muscl": True, "gate": "final",
                    "hold": {"iter": 300, "max": 3}}]}
    s.update(kw)
    return s


class FakeRunner:
    """Подставной шаг: пишет SYNTHETIC историю по сценарию script[(имя шага, порция)] -> dict."""

    def __init__(self, tmp, script=None, default=None):
        self.tmp, self.script, self.calls = tmp, script or {}, []
        self.default = default or {}

    def __call__(self, st, state, chunk):
        self.calls.append((st["name"], chunk, dict(state), st["_p_out"]))
        sc = dict(self.default, **self.script.get((st["name"], chunk), self.script.get(st["name"], {})))
        d = os.path.join(self.tmp, st["_dir"])
        os.makedirs(d, exist_ok=True)
        hp = os.path.join(d, f"history_{st['_dir']}.csv")
        n = sc.pop("n", int(st["iter"]) if chunk == 0 else int(st["hold"]["iter"]))
        status = sc.pop("status", "finished")
        write_history(hp, n, **sc)
        return {"status": status, "seconds": 1.0, "histories": [hp], "restart": os.path.join(d, "restart.dat"),
                "log": {"error": "SYNTHETIC NaN" if status == "diverged" else None}}


# ----------------------------------------------------------------------------------------------- шаблон
def test_derive_template_hist_and_massflow_on_real_template():
    tpl = open(RC.TEMPLATE, encoding="utf-8").read()
    txt, applied, skipped = SW.derive_template(tpl, ["hist", "massflow"])
    assert applied == ["hist", "massflow"] and skipped == []
    hist = [ln for ln in txt.splitlines() if ln.startswith("HISTORY_OUTPUT=")][0]
    for g in ("TURBO_PERF", "CFL_NUMBER", "LINSOL", "MAX_RES"):
        assert hist.count(g) == 1
    giles = [ln for ln in txt.splitlines() if ln.startswith("MARKER_GILES=")][0]
    assert "outlet, MASS_FLOW_OUTLET, ${P_OUT}, 0.1," in giles and "STATIC_PRESSURE_1D" not in giles
    # повторная правка hist не дублирует группы; massflow второй раз — отказ (строки давления уже нет)
    txt2, _, _ = SW.derive_template(txt, ["hist"])
    assert [ln for ln in txt2.splitlines() if ln.startswith("HISTORY_OUTPUT=")][0] == hist
    with pytest.raises(SystemExit):
        SW.derive_template(txt, ["massflow"])


def test_derive_template_soft_and_unknown():
    txt, applied, skipped = SW.derive_template("SOLVER= RANS\n", ["hist"])
    assert applied == [] and skipped == ["hist"]
    for bad in (["nope"], ["riemann"]):
        with pytest.raises(SystemExit):
            SW.derive_template("X\n", bad)


def test_massflow_template_fills_with_run_case_fields():
    """Производный шаблон заполняется теми же полями run_case.py; у Giles 9 значений на маркер (CConfig v8.6.0)."""
    tpl, _, _ = SW.derive_template(open(RC.TEMPLATE, encoding="utf-8").read(), ["hist", "massflow"])
    params = dict(RC.NUMERICS, FS_P="90750", FS_T="279", RESTART="NO", SOLUTION_FILENAME="s.dat", P01="101325", T01="288.15",
                  OMEGA_X="0", OMEGA_Y="0", OMEGA_Z="1500", PER_MASTER="periodic_1", PER_SLAVE="periodic_2", ANGLE="40",
                  P_OUT="0.72", ITER=3, MUSCL="NO", HISTORY_NAME="h", RESTART_NAME="r.dat", VOLUME_NAME="v",
                  MESH="m.su2", RAMP_ROT="NO", RAMP_ROT_INI="0", RAMP_ROT_FREQ=5, RAMP_ROT_ITER=1, RAMP_OUT="NO",
                  RAMP_OUT_INI="90750", RAMP_OUT_FREQ=5, RAMP_OUT_ITER=1)
    params.update(RC.wall_lines("slip", "slip"))
    cfg = RC.build_cfg(tpl, params)
    line = [ln for ln in cfg.splitlines() if ln.startswith("MARKER_GILES=")][0]
    vals = [v.strip() for v in line.split("(", 1)[1].rstrip(")").split(",")]
    assert len(vals) == 18 and vals[10] == "MASS_FLOW_OUTLET" and vals[11] == "0.72" and vals[12] == "0.1"


# ----------------------------------------------------------------------------------------------- давление и рампы
def test_throttle_and_fs_pressure():
    assert SW.throttle_p(1.0, 0.75, 300000.0, 90000.0) == pytest.approx(225000.0)
    assert SW.throttle_p(0.0, 0.75, 300000.0, 90000.0) == pytest.approx(90000.0)
    assert SW.throttle_p(0.5, 0.6, 300000.0, 90000.0) == pytest.approx(90000.0 + 0.25 * 90000.0)
    d1 = {"T01": 288.15, "p01": 101325.0}
    t = 288.15 / (1 + 0.2 * 0.4 ** 2)
    assert SW.fs_pressure(d1, 0.4) == pytest.approx(101325.0 * (t / 288.15) ** 3.5)


def test_ramp_value_matches_su2_update_rule():
    # рампа 8 итераций, обновление раз в 5: скачки на 5-й и 8-й (как в пробе c8f9 при масштабе 0,1)
    vals = [SW.ramp_value(i, 0.0, 0.6, 8, 5) for i in range(10)]
    assert vals[:5] == [0.0] * 5 and vals[5] == pytest.approx(0.375) and vals[7] == pytest.approx(0.375)
    assert vals[8] == vals[9] == 0.6
    assert SW.ramp_value(3, 1.0, 2.0, 0, 5) == 2.0


# ----------------------------------------------------------------------------------------------- стратегии
def test_generated_strategies_valid_and_distinct():
    ids = [s["id"] for s in MS.STRATEGIES]
    assert len(ids) == 6 and len(set(ids)) == 6
    assert len({s["axis"] for s in MS.STRATEGIES}) == 6
    for s in MS.STRATEGIES:
        SW.validate_strategy(s)
    jst = next(s for s in MS.STRATEGIES if s["id"] == "v4_jst")
    assert all(not st.get("muscl") for st in jst["steps"])
    mf = next(s for s in MS.STRATEGIES if s["id"] == "v6_massflow")
    assert [st.get("outlet", "pressure") for st in mf["steps"]].count("massflow") == 4
    assert mf["steps"][4]["load"] == "measured"
    for s in MS.STRATEGIES:   # база у всех одна — 0,75 оценки при номинале
        last = s["steps"][-1]
        assert last["load"] == MS.BASE_LOAD and last.get("frac", 1.0) == 1.0


def test_validate_strategy_rejects_bad():
    base = strat_two_steps()
    SW.validate_strategy(base)
    bad = strat_two_steps(set={"FLOW_SCHEME": "JST"})
    with pytest.raises(SystemExit):
        SW.validate_strategy(bad)                         # JST с MUSCL
    bad = strat_two_steps(transforms=["hist", "riemann"])
    with pytest.raises(SystemExit):
        SW.validate_strategy(bad)                         # Riemann: SU2 v8.6.0 падает
    bad = strat_two_steps()
    bad["steps"][0].update(outlet="massflow", m_frac=0.6, load=None)
    with pytest.raises(SystemExit):
        SW.validate_strategy(bad)                         # расход на выходе с рампой давления
    bad["steps"][0]["pramp"] = 0
    with pytest.raises(SystemExit):
        SW.validate_strategy(bad)                         # после расхода давление — только measured
    bad["steps"].insert(1, {"name": "x_to_p", "frac": 1.0, "load": "measured", "iter": 200})
    SW.validate_strategy(bad)
    bad = strat_two_steps()
    bad["steps"][0]["load"] = "measured"
    with pytest.raises(SystemExit):
        SW.validate_strategy(bad)                         # measured без расхода перед ним
    bad = strat_two_steps()
    bad["steps"][0]["iter"] = 250
    with pytest.raises(SystemExit):
        SW.validate_strategy(bad)                         # окно гейта после рампы не помещается
    bad = strat_two_steps()
    bad["steps"][1]["gate"] = "intermediate"
    with pytest.raises(SystemExit):
        SW.validate_strategy(bad)                         # последний шаг — не база
    bad = strat_two_steps()
    bad["steps"][0]["typo"] = 1
    with pytest.raises(SystemExit):
        SW.validate_strategy(bad)


def test_mini_step():
    s = SW.mini_step({"name": "a", "frac": 1.0, "load": 0.6, "iter": 900, "ramp": 400, "pramp": 800, "muscl": True,
                      "muscl_ramp": 250, "gate": "intermediate", "hold": {"iter": 300, "max": 3}}, 3)
    assert (s["iter"], s["ramp"], s["pramp"], s["muscl_ramp"], s["hold"]["max"]) == (3, 2, 2, 2, 0)
    p = SW.mini_step({"name": "p", "kind": "plan", "plan": "smoke075", "load": 0.75, "hold": {"iter": 300, "max": 40}}, 3)
    assert p["hold"] == {"iter": 3, "max": 1} and p["iter_scale"] == pytest.approx(0.015)
    RC.validate_plan(RC.PLANS["smoke075"], p["iter_scale"])   # план с таким масштабом run_case.py принимает


def test_backoff_steps():
    st = {"name": "r2_w055", "frac": 0.55, "load": 0.6, "iter": 220, "ramp": 120, "pramp": 120, "gate": "intermediate",
          "hold": {"iter": 200, "max": 3}, "cfl": 5.0}
    state = {"frac": 0.30, "p": 1e5, "load": 0.6}
    r1 = SW.backoff_steps(st, state, {})
    assert len(r1) == 1 and r1[0]["ramp"] == 240 and r1[0]["iter"] == 340 and r1[0]["_attempt"] == 1
    assert r1[0]["_set"]["CFL_ADAPT_PARAM"].split(",")[3].strip() == "10"
    r2 = SW.backoff_steps(r1[0], state, {})
    assert [x["name"] for x in r2] == ["r2_w055_m2", "r2_w055_r2"]
    assert r2[0]["frac"] == pytest.approx(0.425) and r2[0]["ramp"] == 120 and r2[1]["frac"] == 0.55
    assert SW.backoff_steps(r2[1], state, {}) == []
    same = dict(st, frac=0.30)
    assert SW.backoff_steps(dict(SW.backoff_steps(same, state, {})[0]), state, {}) == []


# ----------------------------------------------------------------------------------------------- гейт
def test_gate_metrics_and_check(tmp_path):
    p = tmp_path / "h.csv"
    write_history(p, 400)
    _, d = RC.read_history(str(p))
    m = SW.gate_metrics(d, 1.4, 200, 50)
    assert m["sign_ok"] and m["imbalance_mean"] == pytest.approx(0.005) and m["m_in_change"] == 0.0
    assert m["PR_tt"] == pytest.approx(3.2) and 0.7 < m["eta_tt"] < 0.73 and m["avg_cfl"] == 12.0
    g = SW.gate_check(m, "final", -4.0, 4.3)
    assert g["pass"] and g["success"] and g["rho_drop_orders"] == pytest.approx(3.5) and not g["converged_skill"]
    gi = SW.gate_check(SW.gate_metrics(d, 1.4, 100, 20), "intermediate", None, 4.3)
    assert gi["pass"]
    # обратный ток на выходе — не проходит
    write_history(p, 400, backflow=True)
    _, d = RC.read_history(str(p))
    g = SW.gate_check(SW.gate_metrics(d, 1.4, 100, 20), "intermediate", None, 4.3)
    assert not g["pass"] and any("обратный" in r for r in g["reasons"])
    # расход ещё меняется — не проходит
    write_history(p, 400, drift=0.1)
    _, d = RC.read_history(str(p))
    g = SW.gate_check(SW.gate_metrics(d, 1.4, 100, 20), "intermediate", None, 4.3)
    assert not g["pass"] and any("не установился" in r for r in g["reasons"])
    # нефизичный КПД — база не принимается
    write_history(p, 400, pr=0.9, tr=1.01)
    _, d = RC.read_history(str(p))
    g = SW.gate_check(SW.gate_metrics(d, 1.4, 200, 50), "final", -4.0, 4.3)
    assert not g["pass"] and not g["success"]
    # мало строк
    write_history(p, 50)
    _, d = RC.read_history(str(p))
    assert not SW.gate_check(SW.gate_metrics(d, 1.4, 100, 20), "intermediate", None, 4.3)["pass"]


def test_analyse_divergence(tmp_path):
    p = tmp_path / "h.csv"
    write_history(p, 9, backflow=True)
    _, d = RC.read_history(str(p))
    rec = {"frac": 1.0, "frac_from": 0.6, "p_out": 2.5e5, "p_from": 1.5e5, "ramp": 8, "pramp": 16, "su2_error": "x"}
    dv = SW.analyse_divergence(d, rec, 1500.0)
    assert dv["iter_of_stop"] == 9 and dv["frac_at_stop"] == 1.0 and dv["backflow_at_outlet"]
    assert dv["p_out_at_stop"] == pytest.approx(1.5e5 + 5 * 1e5 / 16)
    assert dv["fastest_growing"] is not None


# ----------------------------------------------------------------------------------------------- контроллер
def test_run_strategy_success(tmp_path):
    fr = FakeRunner(str(tmp_path))
    s = SW.run_strategy(strat_two_steps(), dict(CTX), runner=fr)
    assert s["verdict"] == "success" and [c[:2] for c in fr.calls] == [("a_spin", 0), ("b_base", 0)]
    assert fr.calls[0][3] == pytest.approx(0.6 * 300000.0) and fr.calls[1][3] == pytest.approx(0.75 * 300000.0)
    assert fr.calls[1][2]["frac"] == 1.0 and fr.calls[1][2]["restart"].endswith("restart.dat")
    assert s["criteria"]["success"] and len(s["points_full_speed"]) == 2 and s["iters_total"] == 800


def test_run_strategy_holds_then_stuck(tmp_path):
    fr = FakeRunner(str(tmp_path), default={"backflow": True})
    s = SW.run_strategy(strat_two_steps(), dict(CTX), runner=fr)
    assert s["verdict"] == "stuck" and s["verdict_step"] == "a_spin"
    assert [c[1] for c in fr.calls] == [0, 1, 2]          # две порции удержания, затем отказ
    assert fr.calls[1][2]["frac"] == 1.0                  # удержание — с рестарта своей же порции


def test_run_strategy_backoff_after_nan(tmp_path):
    strat = strat_two_steps(backoff={"max_total": 2})
    fr = FakeRunner(str(tmp_path), script={("a_spin", 0): {"status": "diverged", "n": 7, "nan_last": True}})
    s = SW.run_strategy(strat, dict(CTX), runner=fr)
    assert [c[0] for c in fr.calls] == ["a_spin", "a_spin_r1", "b_base"] and s["verdict"] == "success"
    assert s["backoffs"] == 1 and s["divergences"][0]["step"] == "a_spin" and s["divergence"] is None
    assert fr.calls[1][2]["restart"] is None              # откат — к последнему пройденному состоянию (свежий старт)


def test_run_strategy_diverged_without_backoff_and_timeout(tmp_path):
    fr = FakeRunner(str(tmp_path), script={"b_base": {"status": "diverged", "n": 5}})
    s = SW.run_strategy(strat_two_steps(), dict(CTX), runner=fr)
    assert s["verdict"] == "diverged" and s["divergence"]["step"] == "b_base"
    t = [0.0]

    def clock():
        t[0] += 1000.0
        return t[0]
    ctx = dict(CTX, deadline=1500.0)
    s = SW.run_strategy(strat_two_steps(), ctx, runner=FakeRunner(str(tmp_path / "t")), clock=clock)
    assert s["verdict"] == "timeout"


def test_run_strategy_mini_advances(tmp_path):
    fr = FakeRunner(str(tmp_path), default={"backflow": True})
    s = SW.run_strategy(strat_two_steps(), dict(CTX, mini_iter=3), runner=fr)
    assert [c[:2] for c in fr.calls] == [("a_spin", 0), ("b_base", 0), ("b_base", 1)]
    assert s["verdict"] == "reached_base"


def test_step_argv(tmp_path):
    strat = dict(MS.STRATEGIES[4])            # v5_sa
    st = dict(strat["steps"][1])              # b_order2: MUSCL с рампой
    ctx = dict(CTX, case="core130", level="coarse", geom_out="/g", np=32, mpi_args="--use-hwthread-cpus",
               shroud="slip", vaneless="slip", su2_src="/opt/SU2-src", work=str(tmp_path))
    argv = SW.step_argv(ctx, strat, st, "plan.json", 2.0e5, 600, "/o", {"restart": "/r.dat", "frac": 1.0, "p": 1.9e5})
    j = " ".join(argv)
    assert "rc --transforms hist" in j and "--turb SA" in j and "--set LIMITER=VENKATAKRISHNAN_WANG" in j
    assert "--extra=RAMP_MUSCL= YES" in j and "--extra=RAMP_MUSCL_COEFF= (0, 1, 250)" in j
    assert "--restart-from /r.dat --restart-state 1,190000" in j and "--p-out 200000" in j


# ----------------------------------------------------------------------------------------------- сводка и пакет
def test_merge(tmp_path):
    def summ(sid, lab, verdict, imb):
        d = tmp_path / sid / lab
        d.mkdir(parents=True)
        (d / "sweep_summary.json").write_text(json.dumps({
            "strategy": {"id": sid, "title": sid, "axis": "a"}, "label": lab, "verdict": verdict,
            "criteria": {"imbalance_mean": imb, "rho_drop_orders": 3.1}, "steps": [{"decision": "advance"}]}))
    summ("v1", "core130", "success", 0.004)
    summ("v1", "kv3_stage1", "diverged", None)
    summ("v2", "core130", "established", 0.008)
    summ("v2", "kv3_stage1", "success", 0.009)
    import glob
    res = SW.merge(glob.glob(str(tmp_path / "**" / "sweep_summary.json"), recursive=True))
    assert res["ranking"][0] == "v2" and res["best_per_wheel"] == {"core130": "v1", "kv3_stage1": "v2"}
    with pytest.raises(ValueError):
        SW.merge(glob.glob(str(tmp_path / "**" / "sweep_summary.json"), recursive=True) * 2)


def test_jobs_file_loads_in_runner(tmp_path):
    import runner as R
    doc = MS.build()
    assert len(doc["jobs"]) == 6
    for j in doc["jobs"]:
        assert j["max_min"] == MS.MAX_MIN and "sweep_ramp.py merge" in j["command"]
        assert j["command"].count("sweep_ramp.py run") == 2 and "--label kv3_stage1" in j["command"]
    p = tmp_path / "jobs.json"
    p.write_text(json.dumps(doc, ensure_ascii=False))
    jobs = R.load_jobs(str(p))
    assert [j.name for j in jobs] == [j["name"] for j in doc["jobs"]]
    e = doc["estimate"]
    assert e["upper_bound_usd"] < 20 and e["usd"]["total_expected"] < e["upper_bound_usd"]
    assert all(r["wheels"]["kv3_stage1"]["sizing"]["mandatory_s"] <= MS.TIME_LIMIT_S for r in e["jobs"])
    # срок стратегии + подготовка + выгрузка раннера помещаются в сторожа
    assert MS.SETUP_MIN + MS.TIME_LIMIT_S / 60 + MS.FINISH_MIN + 8 <= MS.MAX_MIN
    assert math.isfinite(e["machine_hours_expected"])


def test_run_strategy_massflow_then_measured_pressure(tmp_path):
    strat = {"id": "t_mf", "title": "SYNTHETIC", "transforms": ["hist"], "set": {},
             "steps": [{"name": "a_mf", "frac": 1.0, "load": None, "iter": 300, "ramp": 150, "outlet": "massflow",
                        "m_frac": 0.6, "hold": {"iter": 200, "max": 1}},
                       {"name": "b_to_p", "frac": 1.0, "load": "measured", "iter": 200, "hold": {"iter": 200, "max": 1}},
                       {"name": "c_base", "frac": 1.0, "load": 0.75, "iter": 500, "pramp": 200, "gate": "final",
                        "hold": {"iter": 300, "max": 1}}]}
    fr = FakeRunner(str(tmp_path), script={"a_mf": {"p_meas": 205000.0}, "b_to_p": {"p_meas": 205100.0}})
    s = SW.run_strategy(strat, dict(CTX), runner=fr)
    assert s["verdict"] == "success"
    assert fr.calls[0][3] == pytest.approx(0.6 * 1.2)            # цель — доля x расход канала x w
    assert fr.calls[1][3] == pytest.approx(205000.0)             # давление — измеренное в конце шага с расходом
    assert fr.calls[1][2]["p"] == pytest.approx(205000.0) and fr.calls[1][2]["outlet"] == "pressure"
    assert fr.calls[2][2]["p"] == pytest.approx(205000.0) and fr.calls[2][3] == pytest.approx(225000.0)
    assert s["steps"][0]["outlet"] == "massflow" and s["points_full_speed"][0]["p_out"] == pytest.approx(205000.0)


def test_step_argv_massflow_transform(tmp_path):
    strat = next(x for x in MS.STRATEGIES if x["id"] == "v6_massflow")
    ctx = dict(CTX, case="core130", level="coarse", geom_out="/g", np=32, mpi_args="--use-hwthread-cpus",
               shroud="slip", vaneless="slip", su2_src=None, work=str(tmp_path))
    j = " ".join(SW.step_argv(ctx, strat, strat["steps"][0], "p.json", 0.36, 600, "/o", {"restart": None}))
    assert "--transforms hist,massflow" in j and "--restart-from" not in j
    j = " ".join(SW.step_argv(ctx, strat, strat["steps"][4], "p.json", 2.1e5, 600, "/o", {"restart": None}))
    assert "--transforms hist " in j


def test_char_points_base_strategy(tmp_path):
    import char_points as C
    a = C.build_parser().parse_args(["run", "--case", "core130", "--level", "medium", "--geom-out", "/g", "--pfracs",
                                     "0.7,0.9", "--work", "/w", "--out", "/o", "--base-strategy", "plan_sweep_v2_decouple.json",
                                     "--su2-src", "/s"])
    argv = C.base_strategy_argv(a, nproc=64, work="/w/base", out="/o/base", time_limit=1800.4, label="core130")
    j = " ".join(argv)
    assert "sweep_ramp.py run --strategy plan_sweep_v2_decouple.json" in j and "--np 64" in j
    assert "--time-limit 1800" in j and "--su2-src /s" in j and "--work /w/base --out /o/base" in j
    # сводка стратегии -> база: годится только success / established при номинальной частоте и с рестартом
    d = tmp_path / "base"
    d.mkdir()
    (d / "final_restart.dat").write_bytes(b"SYNTHETIC")
    summ = {"verdict": "established", "final_restart": "final_restart.dat",
            "final_state": {"frac": 1.0, "p_out": 233000.0}, "strategy": {"id": "v2_decouple"}}
    (d / "sweep_summary.json").write_text(json.dumps(summ))
    base, rf, p = C.base_from_sweep(str(d / "sweep_summary.json"), 0)
    assert rf.endswith("final_restart.dat") and p == 233000.0 and base["restart_ok"]
    (d / "sweep_summary.json").write_text(json.dumps(dict(summ, verdict="reached_base")))
    base, rf, p = C.base_from_sweep(str(d / "sweep_summary.json"), 0)
    assert rf is None and p is None and not base["restart_ok"]
    base, rf, p = C.base_from_sweep(str(tmp_path / "nope.json"), 3)
    assert base["verdict"] == "no_summary" and rf is None


def test_run_strategy_deadline_on_established_base(tmp_path):
    """Срок кончился во время удержания базы, а гейт базы уже пройден (без 3 порядков) — established, не timeout."""
    t = [0.0]

    def clock():
        t[0] += 100.0
        return t[0]
    fr = FakeRunner(str(tmp_path), default={"rho1": -5.0})        # невязка упала только на порядок
    s = SW.run_strategy(strat_two_steps(), dict(CTX, deadline=1000.0), runner=fr, clock=clock)
    assert s["verdict"] == "established" and s["criteria"]["established"] and not s["criteria"]["success"]
