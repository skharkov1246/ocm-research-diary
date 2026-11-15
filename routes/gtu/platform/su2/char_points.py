#!/usr/bin/env python3
"""Линия частоты (характеристика) колеса в SU2 на одной машине и сводка таких прогонов.

Узел platform/su2 трека ГТУ. Обёртка над run_case.py, сам SU2 запускает run_case.py.

  run    — на машине (AWS, контейнер с SU2 и Python): сетка уже построена (--geom-out).
           1) база: план smoke075 run_case.py на всех ядрах — разгон частоты 0 -> 0,6 -> 1,0 и 150 итераций
              второго порядка при давлении на выходе p_base = base_pfrac x оценка (как прогон Б, README);
           2) точки: параллельно, по --np-per-point процессов MPI, каждая — рестарт от базы, рампа давления на
              выходе p_base -> p_j = pfrac_j x оценка за --point-pramp итераций и --point-iter итераций второго
              порядка всего; точку характеристики задаёт давление на выходе;
           3) char_summary.json и выбранные файлы в --out: сводки, history, журналы, cfg, итоговый рестарт точки,
              поверхность (если SU2 её писал). Сетки, копии сеток run_case.py и объёмные поля в --out не идут.
           Срок: --time-limit, или от срока сторожа задачи (--job-max-min) минус время с загрузки машины
           (/proc/uptime) и запасы — SU2 останавливается сам и успевает записать рестарт (README, «Останов»).
  merge  — в сессии: char_summary.json всех задач прогона -> один results-JSON (по случаю и уровню сетки),
           с проверкой: каждая заданная точка есть ровно один раз.

Оценка давления на выходе — run_case.estimate_p_out, но для двустороннего колеса — по расходу одной стороны
(d1.m_side): сетка канала — одна сторона, ширина b2 — одной стороны. run_case.py берёт d1.m (расход колеса) и
для KV-3 ступени 1 даёт оценку на 9 % ниже; поэтому давления точек передаются ему явно (--p-out).

Примеры:
  python char_points.py run --case core130 --level medium --geom-out W/geom/core130_medium \
         --pfracs 0.5,0.7,0.9,1.1 --job-max-min 360 --work W/runs --out "$GTU_OUT"
  python char_points.py merge --in ../aws/runs/<run_id> --out su2_char_kv3_results.json --run-id <run_id>
"""
from __future__ import annotations

import argparse
import glob
import json
import math
import os
import shutil
import signal
import subprocess
import sys
import time
from datetime import datetime, timezone

HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in sys.path:
    sys.path.insert(0, HERE)
sys.dont_write_bytecode = True
import run_case as RC  # noqa: E402

RUNNER_UPLOAD_RESERVE_S = 8 * 60   # runner.py: команда идёт до «старт + max_min − 8 мин» (upload_reserve_min)
BASE_PLAN = "smoke075"             # прошёл дымовую проверку на core130 coarse (README, прогон Б)
WINDOW_CONV = 500                  # окно сходимости навыка su2-turbo: π и η меняются < 0,1 % за 500 итераций
CONV_TOL = 1e-3
# Сетку каждая задача строит сама, а 3D-сетка Gmsh (HXT, Mesh.Algorithm3D = 10, несколько потоков) от запуска к
# запуску разная: core130 coarse — 139 121 / 139 087 / 139 010 узлов, KV-3 coarse — 258 565 / 258 411 / 258 299
# (bench, AWS 15.11.2026, миниатюры). Задачи одной линии сводятся, если узлы расходятся не больше чем на MESH_TOL;
# узлы каждой задачи — в сводке (mesh_nodes_by_job, у точки — mesh_nodes).
MESH_TOL = 0.01
KEEP_PATTERNS = ("run_summary.json", "history_*.csv", "su2_*.log", "dry_*.log", "case_*.cfg", "surface*.vtu", "surface*.csv")
SKIP_ALWAYS = ("mesh_vaneless_split.su2", "restart_in.dat")


# ----------------------------------------------------------------------------- чистые функции
def parse_fracs(text: str) -> list:
    vals = [float(x) for x in text.split(",") if x.strip()]
    if not vals or any(not (0.2 <= v <= 2.0) for v in vals):
        raise ValueError(f"доли давления на выходе — числа 0,2…2,0 через запятую: {text!r}")
    if len(set(round(v, 6) for v in vals)) != len(vals):
        raise ValueError(f"повтор доли давления: {text!r}")
    return vals


def p_out_estimate(mj: dict, gamma: float = 1.4, R: float = 287.058) -> dict:
    """Оценка статического давления на выходе при номинальной частоте (run_case.estimate_p_out); у двустороннего
    колеса — по расходу одной стороны."""
    d1 = dict(mj["d1"])
    r_out = float(mj.get("contours", {}).get("r_out", 1.3 * d1["r2"]))
    dbl = bool(d1.get("double_sided"))
    est_wheel = RC.estimate_p_out(d1, 1.0, r_out, d1["b2"], gamma, R)
    if dbl:
        d1["m"] = d1["m_side"]
    est = RC.estimate_p_out(d1, 1.0, r_out, d1["b2"], gamma, R)
    return {"p_out_estimate": est, "double_sided": dbl, "m_for_estimate": d1["m"], "r_out": r_out,
            "p_out_estimate_run_case": est_wheel,
            "note": ("двустороннее колесо: оценка по расходу стороны d1.m_side (run_case.py берёт d1.m — расход колеса)"
                     if dbl else "оценка run_case.estimate_p_out")}


def point_plan(n_iter: int, pramp: int) -> list:
    return [dict(name="c1_point", frac=1.0, pfrac=1.0, ramp=0, pramp=int(pramp), iter=int(n_iter), muscl=True)]


def time_budget(job_max_min: int | None, time_limit: float | None, uptime_s: float, reserve_s: float) -> float:
    """Секунд на базу и точки. Явный --time-limit главнее; иначе срок команды раннера (max_min − 8 мин) минус
    время с загрузки машины (оно больше времени от старта user-data — оценка с запасом) и запас на сводку."""
    if time_limit is not None:
        return float(time_limit)
    if not job_max_min:
        raise ValueError("нужен --time-limit или --job-max-min")
    left = job_max_min * 60 - RUNNER_UPLOAD_RESERVE_S - uptime_s - reserve_s
    if left < 120:
        raise ValueError(f"срок исчерпан: max_min {job_max_min}, с загрузки машины {uptime_s:.0f} с, "
                         f"запас {reserve_s:.0f} с — осталось {left:.0f} с")
    return left


def read_uptime(path: str = "/proc/uptime") -> float:
    with open(path, encoding="ascii") as f:
        return float(f.read().split()[0])


def window_change(vals: list, window: int) -> float | None:
    """(max − min)/|последнее| на последних window+1 значениях; None — мало строк или нечисла."""
    if len(vals) < window + 1:
        return None
    seg = vals[-window - 1:]
    if not all(math.isfinite(x) for x in seg) or seg[-1] == 0:
        return None
    return (max(seg) - min(seg)) / abs(seg[-1])


def history_windows(path: str, gamma: float, window: int = WINDOW_CONV) -> dict:
    """Изменение расхода, π, T0 на выходе и КПД за последние window итераций (своё определение η, как в
    run_case.analyse_history)."""
    _, d = RC.read_history(path)
    m = d.get("Avg_Massflow(inlet)", [])
    p1, p2 = d.get("Avg_TotalPress(inlet)", []), d.get("Avg_TotalPress(outlet)", [])
    t1, t2 = d.get("Avg_TotalTemp(inlet)", []), d.get("Avg_TotalTemp(outlet)", [])
    n = min(len(m), len(p1), len(p2), len(t1), len(t2))
    pr, eta = [], []
    for i in range(n):
        try:
            a, b = p2[i] / p1[i], t2[i] / t1[i]
            pr.append(a)
            eta.append((a ** ((gamma - 1) / gamma) - 1) / (b - 1) if b > 1 else float("nan"))
        except (ZeroDivisionError, ValueError, OverflowError):
            pr.append(float("nan")); eta.append(float("nan"))
    out = {"window_iters": window, "rows": n,
           "rel_change": {"m_in": window_change(m, window), "PR_tt": window_change(pr, window),
                          "tt_out": window_change(t2, window), "eta_tt": window_change(eta, window)}}
    rc = out["rel_change"]
    out["converged_pi_eta"] = bool(rc["PR_tt"] is not None and rc["eta_tt"] is not None
                                   and rc["PR_tt"] < CONV_TOL and rc["eta_tt"] < CONV_TOL)
    return out


def point_record(summ: dict, run_dir: str, pfrac: float, p_out: float, double_sided: bool) -> dict:
    """Итог точки по run_summary.json run_case.py: последняя стадия с историей, расход колеса (у двустороннего — две
    стороны), признаки сходимости."""
    stages = summ.get("stages", [])
    st = next((s for s in reversed(stages) if s.get("history", {}).get("rows")), None)
    rec = {"pfrac": pfrac, "p_out_Pa": p_out, "statuses": [s.get("status") for s in stages],
           "status": stages[-1].get("status") if stages else "no_stages", "dir": os.path.basename(run_dir)}
    if st is None:
        rec["history_rows"] = 0
        return rec
    h, perf = st["history"], dict(st["history"].get("performance", {}))
    sides = 2 if double_sided else 1
    n_pass = summ["inputs"]["n_passages"]
    if RC._finite(perf.get("m_in_passage")):
        perf["m_wheel_total"] = perf["m_in_passage"] * n_pass * sides
        if RC._finite(perf.get("m_wheel_corrected")):
            perf["m_wheel_corrected_total"] = perf["m_wheel_corrected"] * sides
    rho = h.get("residuals", {}).get("rms[Rho]", {})
    rec.update(stage=st["name"], history_rows=h.get("rows"), seconds=st.get("seconds"),
               sec_per_iter=(st["seconds"] / h["rows"]) if st.get("seconds") and h.get("rows") else None,
               rms_rho={k: rho.get(k) for k in ("first", "max", "last", "slope_per_100_last_half")},
               performance=perf, mass_sign_ok=perf.get("mass_sign_ok"),
               periodic={k: st.get("log", {}).get(k) for k in ("periodic_matched_points", "periodic_max_distance",
                                                               "bad_match")},
               nan=st.get("log", {}).get("nan"), su2_error=st.get("log", {}).get("error"))
    hp = os.path.join(run_dir, st["history_file"])
    if os.path.isfile(hp):
        rec["window"] = history_windows(hp, summ["inputs"]["gamma"])
    rec["usable"] = bool(rec["status"] in ("finished", "timeout") and not rec.get("nan")
                         and perf.get("mass_sign_ok") and RC._finite(perf.get("mass_imbalance"))
                         and perf["mass_imbalance"] <= 0.05)
    rec["converged"] = bool(rec["usable"] and rec.get("window", {}).get("converged_pi_eta"))
    return rec


# ----------------------------------------------------------------------------- запуск
_CHILDREN: list = []


def _forward(signum, frame):  # noqa: ARG001
    """SIGTERM/SIGINT/SIGHUP обёртке -> SIGTERM каждому run_case.py: он останавливает SU2 так, что рестарт пишется."""
    for p in list(_CHILDREN):
        if p.poll() is None:
            try:
                p.send_signal(signal.SIGTERM)
            except ProcessLookupError:
                pass


def run_case_argv(a, *, out: str, np_: int, p_out: float, plan: str, time_limit: float,
                  restart: tuple | None = None, iter_scale: float = 1.0) -> list:
    cmd = [sys.executable, "-u", os.path.join(HERE, "run_case.py"), "--case", a.case, "--level", a.level,
           "--geom-out", a.geom_out, "--plan", plan, "--p-out", f"{p_out:.6g}", "--np", str(np_),
           f"--mpi-args={a.mpi_args}", "--shroud", a.shroud, "--vaneless", a.vaneless,
           "--time-limit", f"{max(time_limit, 1.0):.0f}", "--out", out, "--iter-scale", f"{iter_scale:g}",
           "--set", f"OUTPUT_FILES={a.output_files}"]
    for x in a.extra:
        cmd.append(f"--extra={x}")
    if a.su2_src:
        cmd += ["--su2-src", a.su2_src]
    if a.turb:
        cmd += ["--turb", a.turb]
    if restart:
        cmd += ["--restart-from", restart[0], "--restart-state", f"1.0,{restart[1]:.6g}"]
    return cmd


def keep_files(src: str, dst: str, keep_restart: str | None) -> list:
    """Сводки, history, журналы, cfg, поверхность, итоговый рестарт (по имени) — в dst. Возвращает имена."""
    os.makedirs(dst, exist_ok=True)
    got = []
    for pat in KEEP_PATTERNS:
        for f in sorted(glob.glob(os.path.join(src, pat))):
            if os.path.basename(f) in SKIP_ALWAYS:
                continue
            shutil.copy2(f, os.path.join(dst, os.path.basename(f)))
            got.append(os.path.basename(f))
    if keep_restart and os.path.isfile(os.path.join(src, keep_restart)):
        shutil.move(os.path.join(src, keep_restart), os.path.join(dst, keep_restart))
        got.append(keep_restart)
    return got


def drop_big(run_dir: str) -> None:
    """Копия сетки (run_case --vaneless) и копия чужого рестарта — сотни МБ на точку; удаляются после точки."""
    for f in SKIP_ALWAYS:
        try:
            os.remove(os.path.join(run_dir, f))
        except FileNotFoundError:
            pass


def cmd_run(a) -> int:
    t_start = time.monotonic()
    pfracs = parse_fracs(a.pfracs)
    su2_mesh, js = RC.mesh_paths(a.case, a.level, a.geom_out)
    if not (os.path.isfile(su2_mesh) and os.path.isfile(js)):
        raise SystemExit(f"нет сетки {su2_mesh} или {js}: сетку строит команда задачи до обёртки")
    mj = json.load(open(js, encoding="utf-8"))
    est = p_out_estimate(mj)
    p_base = a.base_pfrac * est["p_out_estimate"]
    nproc = a.nproc or int(os.environ.get("NPROC") or os.cpu_count() or 1)
    np_pt = a.np_per_point or max(1, nproc // len(pfracs))
    if np_pt * len(pfracs) > nproc:
        raise SystemExit(f"{len(pfracs)} точек x {np_pt} процессов > {nproc} ядер")
    uptime = read_uptime(a.uptime_file) if a.time_limit is None else 0.0
    total = time_budget(a.job_max_min, a.time_limit, uptime, a.reserve_s)
    deadline = t_start + total
    os.makedirs(a.work, exist_ok=True)
    os.makedirs(a.out, exist_ok=True)
    for sig in (signal.SIGTERM, signal.SIGINT, signal.SIGHUP):
        signal.signal(sig, _forward)
    summ = {
        "created_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "label": a.label or f"{a.case}", "case": a.case, "level": a.level, "job": os.environ.get("GTU_JOB"),
        "run_id": os.environ.get("GTU_RUN_ID"), "itype": os.environ.get("GTU_ITYPE"),
        "az": os.environ.get("GTU_AZ_MACHINE") or os.environ.get("GTU_AZ"),
        "mesh": {"file": os.path.basename(su2_mesh), "json": os.path.basename(js),
                 "header": RC.read_mesh_header(su2_mesh), "double_sided": est["double_sided"]},
        "d1_reference": {k: mj["d1"].get(k) for k in ("m", "m_side", "PR", "rpm", "D2", "U2", "b2", "T01", "p01",
                                                      "eta_comp_stage")},
        "estimate": est, "base_pfrac": a.base_pfrac, "p_base_Pa": p_base,
        "pfracs_assigned": pfracs, "nproc": nproc, "np_per_point": np_pt,
        "budget": {"time_s": round(total, 1), "uptime_at_start_s": round(uptime, 1), "job_max_min": a.job_max_min,
                   "reserve_s": a.reserve_s, "base_max_frac": a.base_max_frac},
        "plans": {"base": a.base_plan, "point": point_plan(a.point_iter, a.point_pramp),
                  "iter_scale_base": a.base_iter_scale, "iter_scale_point": a.iter_scale},
        "settings": {"shroud": a.shroud, "vaneless": a.vaneless, "turb": a.turb or "SST",
                     "output_files": a.output_files, "extra": a.extra, "mpi_args": a.mpi_args},
        "base": None, "points": [],
    }

    def save():
        with open(os.path.join(a.out, "char_summary.json"), "w", encoding="utf-8") as f:
            json.dump(summ, f, ensure_ascii=False, indent=1)

    # 1) база на всех ядрах
    base_dir = os.path.join(a.work, "base")
    tl_base = min(a.base_max_frac * total, deadline - time.monotonic() - 60)
    print(f"[char] {summ['label']} {a.level}: оценка p_вых {est['p_out_estimate']:.0f} Па, база {p_base:.0f} Па "
          f"({a.base_pfrac:g}), точки {pfracs}, {np_pt} процессов на точку; срок {total:.0f} с, база до {tl_base:.0f} с"
          + (f"; база — стратегия {os.path.basename(a.base_strategy)}" if a.base_strategy else ""), flush=True)
    if a.base_strategy:   # база — стратегия перебора sweep_ramp.py (переходы по гейту расхода), а не план run_case.py
        summ["plans"]["base"] = os.path.basename(a.base_strategy)
        p = subprocess.Popen(base_strategy_argv(a, nproc=nproc, work=base_dir, out=os.path.join(a.out, "base"),
                                                time_limit=tl_base, label=summ["label"]))
        _CHILDREN.append(p)
        rc_base = p.wait()
        _CHILDREN.clear()
        base, restart, p_sw = base_from_sweep(os.path.join(a.out, "base", "sweep_summary.json"), rc_base)
        if p_sw is not None:
            p_base = p_sw
            summ["p_base_Pa"] = p_base
        summ["base"] = base
        save()
        if restart is None:
            print(f"[char] база-стратегия: {base.get('verdict')} — точки не считаются", flush=True)
            summ["verdict"] = "base_failed"
            save()
            return 3
        return run_points(a, summ, save, est, pfracs, np_pt, deadline, restart, p_base, t_start)
    p = subprocess.Popen(run_case_argv(a, out=base_dir, np_=nproc, p_out=p_base, plan=a.base_plan,
                                       time_limit=tl_base, iter_scale=a.base_iter_scale))
    _CHILDREN.append(p)
    rc_base = p.wait()
    _CHILDREN.clear()
    base_summ_path = os.path.join(base_dir, "run_summary.json")
    base = {"returncode": rc_base, "dir": "base"}
    restart = None
    if os.path.isfile(base_summ_path):
        bs = json.load(open(base_summ_path, encoding="utf-8"))
        base["statuses"] = [s.get("status") for s in bs["stages"]]
        base["smoke_checks"] = bs.get("smoke_checks")
        base["record"] = point_record(bs, base_dir, a.base_pfrac, p_base, est["double_sided"])
        last = bs["stages"][-1]
        rf = os.path.join(base_dir, f"restart_{last['name']}.dat")
        if all(s == "finished" for s in base["statuses"]) and os.path.isfile(rf):
            restart = rf
        base["su2_version"] = bs.get("su2", {}).get("version")
        base["key_check"] = bs.get("key_check")
    base["files"] = keep_files(base_dir, os.path.join(a.out, "base"), None)
    drop_big(base_dir)
    base["restart_ok"] = restart is not None
    summ["base"] = base
    save()
    if restart is None:
        print("[char] база не дошла до конца — точки не считаются", flush=True)
        summ["verdict"] = "base_failed"
        save()
        return 3
    return run_points(a, summ, save, est, pfracs, np_pt, deadline, restart, p_base, t_start)


def base_strategy_argv(a, *, nproc: int, work: str, out: str, time_limit: float, label: str) -> list:
    """sweep_ramp.py run для базы по стратегии перебора (plan_sweep_*.json). Ключи SU2 (--extra, --output-files,
    --turb) задаёт сама стратегия; давление базы — её итоговый шаг (load x та же оценка char_points.p_out_estimate)."""
    cmd = [sys.executable, "-u", os.path.join(HERE, "sweep_ramp.py"), "run", "--strategy", a.base_strategy,
           "--case", a.case, "--level", a.level, "--geom-out", a.geom_out, "--label", label, "--np", str(nproc),
           f"--mpi-args={a.mpi_args}", "--shroud", a.shroud, "--vaneless", a.vaneless,
           "--time-limit", f"{max(time_limit, 1.0):.0f}", "--work", work, "--out", out]
    if a.su2_src:
        cmd += ["--su2-src", a.su2_src]
    return cmd


BASE_STRATEGY_OK = ("success", "established")   # база годится для точек: расход согласован и установился, п и КПД физичны


def base_from_sweep(summary_path: str, returncode: int) -> tuple:
    """(запись базы, рестарт или None, давление базы или None) по sweep_summary.json стратегии."""
    base = {"returncode": returncode, "dir": "base", "kind": "sweep_strategy"}
    if not os.path.isfile(summary_path):
        base["verdict"] = "no_summary"
        base["restart_ok"] = False
        return base, None, None
    s = json.load(open(summary_path, encoding="utf-8"))
    fs = s.get("final_state") or {}
    base.update(verdict=s.get("verdict"), verdict_detail=s.get("verdict_detail"), verdict_step=s.get("verdict_step"),
                strategy=s.get("strategy"), criteria=s.get("criteria"), su2_version=s.get("su2_version"),
                key_check=s.get("key_check"), iters_total=s.get("iters_total"), final_state=fs)
    rf = os.path.join(os.path.dirname(summary_path), s["final_restart"]) if s.get("final_restart") else None
    ok = (s.get("verdict") in BASE_STRATEGY_OK and rf and os.path.isfile(rf) and abs(float(fs.get("frac", 0)) - 1.0) < 1e-9
          and fs.get("outlet", "pressure") == "pressure")
    base["restart_ok"] = bool(ok)
    return base, (rf if ok else None), (float(fs["p_out"]) if ok else None)


def run_points(a, summ, save, est, pfracs, np_pt, deadline, restart, p_base, t_start) -> int:
    """2) точки параллельно от рестарта базы."""
    plan_file = os.path.join(a.work, "plan_point.json")
    with open(plan_file, "w", encoding="utf-8") as f:
        json.dump(point_plan(a.point_iter, a.point_pramp), f)
    tl_pt = deadline - time.monotonic() - 30
    procs = []
    for pf in pfracs:
        tag = f"p{pf:.3f}"
        d = os.path.join(a.work, tag)
        po = pf * est["p_out_estimate"]
        pr = subprocess.Popen(run_case_argv(a, out=d, np_=np_pt, p_out=po, plan=plan_file, time_limit=tl_pt,
                                            restart=(restart, p_base), iter_scale=a.iter_scale),
                              stdout=open(os.path.join(a.work, f"{tag}.out"), "w"), stderr=subprocess.STDOUT)
        _CHILDREN.append(pr)
        procs.append((pf, po, tag, d, pr))
    print(f"[char] запущено точек: {len(procs)}, срок {tl_pt:.0f} с", flush=True)
    for pf, po, tag, d, pr in procs:
        rc = pr.wait()
        sp = os.path.join(d, "run_summary.json")
        if os.path.isfile(sp):
            rec = point_record(json.load(open(sp, encoding="utf-8")), d, pf, po, est["double_sided"])
        else:
            rec = {"pfrac": pf, "p_out_Pa": po, "status": "no_summary", "dir": tag}
        rec["returncode"] = rc
        rec["files"] = keep_files(d, os.path.join(a.out, tag), "restart_c1_point.dat" if a.keep_restart else None)
        shutil.copy2(os.path.join(a.work, f"{tag}.out"), os.path.join(a.out, tag, "run_case.out"))
        drop_big(d)
        summ["points"].append(rec)
        print(f"[char] {tag}: {rec.get('status')}, строк {rec.get('history_rows')}, "
              f"m_колеса {rec.get('performance', {}).get('m_wheel_total')}, PR {rec.get('performance', {}).get('PR_tt')}",
              flush=True)
        save()
    _CHILDREN.clear()
    summ["points"].sort(key=lambda r: r["pfrac"])
    summ["verdict"] = "ok" if any(r.get("usable") for r in summ["points"]) else "no_usable_points"
    summ["t_total_s"] = round(time.monotonic() - t_start, 1)
    save()
    return 0 if summ["verdict"] == "ok" else 4


# ----------------------------------------------------------------------------- сводка прогона
def merge(paths: list, run_id: str | None = None) -> dict:
    """char_summary.json задач -> характеристика по (метка, уровень). Отказ (ValueError), если точка задана
    дважды, задача прислала точку, которой ей не задавали, или сетки задач линии расходятся больше MESH_TOL."""
    groups, problems, jobs = {}, [], []
    for p in sorted(paths):
        s = json.load(open(p, encoding="utf-8"))
        key = f"{s['label']}/{s['level']}"
        nodes = s["mesh"]["header"].get("NPOIN")
        g = groups.setdefault(key, {"label": s["label"], "case": s["case"], "level": s["level"],
                                    "mesh_nodes": nodes,
                                    "mesh_cells": s["mesh"]["header"].get("NELEM"),
                                    "mesh_nodes_by_job": {},
                                    "double_sided": s["mesh"]["double_sided"], "estimate": s["estimate"],
                                    "d1_reference": s["d1_reference"], "settings": s["settings"],
                                    "plans": s["plans"], "assigned": [], "points": [], "bases": []})
        g["mesh_nodes_by_job"][str(s.get("job"))] = nodes
        ref = g["mesh_nodes"]
        if not (isinstance(ref, int) and isinstance(nodes, int) and ref > 0
                and abs(nodes - ref) / ref <= MESH_TOL):
            problems.append(f"{key}: у задач разные сетки ({ref} и {nodes} узлов, допуск {MESH_TOL:.0%})")
        jobs.append({"job": s.get("job"), "file": os.path.relpath(p), "label": s["label"], "level": s["level"],
                     "verdict": s.get("verdict"), "itype": s.get("itype"), "az": s.get("az"),
                     "t_total_s": s.get("t_total_s"), "budget": s.get("budget")})
        g["assigned"] += s["pfracs_assigned"]
        b = s.get("base") or {}
        g["bases"].append({"job": s.get("job"), "restart_ok": b.get("restart_ok"), "statuses": b.get("statuses"),
                           "su2_version": b.get("su2_version"), "record": b.get("record")})
        for r in s.get("points", []):
            if r["pfrac"] not in s["pfracs_assigned"]:
                problems.append(f"{key}: точка {r['pfrac']} не задавалась задаче {s.get('job')}")
            g["points"].append(dict(r, job=s.get("job"), mesh_nodes=nodes))
    for key, g in groups.items():
        seen = [round(x, 6) for x in g["assigned"]]
        dup = sorted({x for x in seen if seen.count(x) > 1})
        if dup:
            problems.append(f"{key}: точки заданы дважды: {dup}")
        got = [round(r["pfrac"], 6) for r in g["points"]]
        g["missing"] = sorted(set(seen) - set(got))
        g["points"].sort(key=lambda r: r["pfrac"])
        g["usable"] = sum(1 for r in g["points"] if r.get("usable"))
        g["converged"] = sum(1 for r in g["points"] if r.get("converged"))
    if problems:
        raise ValueError("; ".join(problems))
    return {"about": ("Линии частоты колёс в SU2 (char_points.py run на AWS, runner.py): база smoke075, точки — рестарт "
                      "от базы с рампой давления на выходе. Сетки узла geom — тетраэдры без пристенных призм; покрышка "
                      "и стенки за колесом — по settings (slip — без трения). Числа — проверка постановки и "
                      "сходимости по сетке, не паспортная характеристика; «converged» — π и η меняются < 0,1 % за 500 "
                      "итераций, «usable» — без NaN, прямое течение, дисбаланс расхода ≤ 5 %."),
            "run_id": run_id, "created_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
            "jobs": jobs, "lines": groups}


def cmd_merge(a) -> int:
    paths = sorted(glob.glob(os.path.join(a.inp, "**", "char_summary.json"), recursive=True))
    paths = [p for p in paths if "/_lost/" not in p.replace(os.sep, "/")]
    if not paths:
        print(f"нет char_summary.json в {a.inp}", file=sys.stderr)
        return 2
    try:
        res = merge(paths, a.run_id)
    except ValueError as e:
        print("ОТКАЗ:", e, file=sys.stderr)
        return 2
    txt = json.dumps(res, ensure_ascii=False, indent=1)
    if len(txt.encode()) > 5 * 1024 * 1024:
        print("ОТКАЗ: сводка больше 5 МБ", file=sys.stderr)
        return 2
    with open(a.out, "w", encoding="utf-8") as f:
        f.write(txt + "\n")
    for key, g in res["lines"].items():
        print(f"{key}: точек {len(g['points'])}, годных {g['usable']}, сошедшихся {g['converged']}, нет {g['missing']}")
        for r in g["points"]:
            pf = r.get("performance", {})
            print(f"  pfrac {r['pfrac']:.3f}: {r.get('status')}, m {pf.get('m_wheel_total')}, PR {pf.get('PR_tt')}, "
                  f"eta {pf.get('eta_tt_isentropic')}, сошлась {r.get('converged')}")
    print(f"записано: {a.out}")
    if a.expect:
        labels = {g["label"] for g in res["lines"].values()}
        miss = [x for x in a.expect.split(",") if x and x not in labels]
        if miss:
            print("нет сводок линий:", miss, file=sys.stderr)
            return 3
    if a.require_base_ok:
        bad = [f"{k}: {b.get('job')}" for k, g in res["lines"].items() for b in g["bases"] if not b.get("restart_ok")]
        if bad:
            print("база не дошла до конца:", bad, file=sys.stderr)
            return 3
    return 0


def build_parser():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    r = sub.add_parser("run", help="база + точки на одной машине")
    r.add_argument("--case", required=True, help="как в run_case.py: core130 | taurus60 | custom")
    r.add_argument("--level", required=True, choices=["coarse", "medium", "fine"])
    r.add_argument("--geom-out", required=True, help="папка сетки (passage_<case>_<level>.su2 и JSON)")
    r.add_argument("--label", default=None, help="имя линии в сводке (например kv3_stage1)")
    r.add_argument("--pfracs", required=True, help="доли оценки давления на выходе, через запятую")
    r.add_argument("--base-pfrac", type=float, default=0.75)
    r.add_argument("--base-plan", default=BASE_PLAN, help="план базы для run_case.py (имя или JSON)")
    r.add_argument("--base-iter-scale", type=float, default=1.0)
    r.add_argument("--base-strategy", default=None,
                   help="база — стратегия перебора sweep_ramp.py (plan_sweep_*.json) вместо --base-plan: переходы по "
                        "гейту расхода; точки — от её итогового рестарта при её давлении (вердикт success или established)")
    r.add_argument("--base-max-frac", type=float, default=0.3, help="доля срока на базу")
    r.add_argument("--point-iter", type=int, default=6000)
    r.add_argument("--point-pramp", type=int, default=300)
    r.add_argument("--iter-scale", type=float, default=1.0, help="множитель итераций точек (проба)")
    r.add_argument("--nproc", type=int, default=0, help="ядер всего (по умолчанию $NPROC)")
    r.add_argument("--np-per-point", type=int, default=0, help="процессов на точку (по умолчанию ядра / точки)")
    r.add_argument("--mpi-args", default="--use-hwthread-cpus",
                   help="ключи mpirun: слоты по потокам — на c7i.16xlarge (32 ядра x 2) 64 процесса без отказа")
    r.add_argument("--shroud", default="slip", choices=["stationary", "rotating", "slip"])
    r.add_argument("--vaneless", default="slip", choices=["as-is", "stationary", "slip"])
    r.add_argument("--turb", default=None, choices=["SST", "SA"])
    r.add_argument("--output-files", default="RESTART, SURFACE_PARAVIEW",
                   help="OUTPUT_FILES SU2: без объёмного PARAVIEW (поле восстанавливается из рестарта SU2_SOL)")
    r.add_argument("--extra", action="append", default=["MARKER_PLOTTING= (blade, hub)"],
                   help="строка SU2 вне шаблона (run_case.py --extra)")
    r.add_argument("--su2-src", default=None, help="исходник SU2 — сверка ключей cfg")
    r.add_argument("--time-limit", type=float, default=None, help="секунд на базу и точки (иначе от --job-max-min)")
    r.add_argument("--job-max-min", type=int, default=None, help="max_min задачи раннера (срок сторожа)")
    r.add_argument("--reserve-s", type=float, default=300.0, help="запас на останов SU2, копирование и сводку")
    r.add_argument("--uptime-file", default="/proc/uptime")
    r.add_argument("--keep-restart", action=argparse.BooleanOptionalAction, default=True,
                   help="итоговый рестарт точки — в --out (по умолчанию да)")
    r.add_argument("--work", required=True, help="рабочая папка (не выгружается)")
    r.add_argument("--out", required=True, help="папка итогов ($GTU_OUT)")
    m = sub.add_parser("merge", help="char_summary.json задач прогона -> results-JSON")
    m.add_argument("--in", dest="inp", required=True, help="папка runs/<run_id> раннера")
    m.add_argument("--out", required=True)
    m.add_argument("--run-id", default=None)
    m.add_argument("--require-base-ok", action="store_true", help="код 3, если у какой-то задачи база не дошла до конца "
                   "(проба на AWS)")
    m.add_argument("--expect", default="", help="метки линий через запятую: код 3, если какой-то нет")
    return ap


def main(argv=None) -> int:
    a = build_parser().parse_args(argv)
    if a.cmd == "run":
        return cmd_run(a)
    return cmd_merge(a)


if __name__ == "__main__":
    sys.exit(main())
