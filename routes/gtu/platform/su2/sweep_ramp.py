#!/usr/bin/env python3
"""Разгон базовой точки колеса в SU2 по стратегии с переходом между стадиями по установлению расхода и сводка
перебора стратегий (узел platform/su2 трека ГТУ). Обёртка над run_case.py: run_case.py и шаблон не меняются.

Зачем. Проба на AWS 20261115T003035Z-c8f9 (jobs_aws_su2_kv3_smoke.json): база обоих колёс разошлась. План smoke075
шёл с --base-iter-scale 0.1, а validate_plan масштабирует и длины рамп: разгон 0 -> 0,6 частоты за 8 итераций,
0,6 -> 1,0 — тоже за 8, давление на выходе — за 16; переход между стадиями — по числу итераций, а не по расходу
(разбор — jobs_aws_su2_kv3_sweep.json, блок diagnosis). Здесь каждая стадия — отдельный запуск run_case.py с рестартом
от предыдущей, а следующая стадия начинается, только когда гейт по истории пройден (расход на входе и выходе
согласован, расход установился, невязка не растёт); не пройден — ещё порция итераций при тех же условиях.

  run    — на машине: одна стратегия (plan_sweep_*.json) для одного колеса -> sweep_summary.json и файлы шагов в --out.
  merge  — сводка sweep_summary.json всех задач -> один JSON (стратегия x колесо), ранжирование стратегий.
  rc     — тонкая обёртка run_case.py: производный шаблон (история с CFL, линейным решателем и макс. невязками;
           вход и выход Riemann вместо Giles); вызывается из run, вручную не нужна.

Давление на выходе. База — как в char_points.py: p_base = доля x оценка (char_points.p_out_estimate, у двустороннего
колеса — по расходу стороны). Промежуточные состояния идут по линии дросселя p(w) = p_нач + w^2 (доля x оценка -
p_нач), где p_нач — статическое давление начального поля (run_case.py: M = INIT_MACH вдоль оси), w — доля частоты:
при малой частоте давление на выходе не уходит ниже p_нач (у доли x оценка(w) уходило бы — сильное разрежение
на выходе и запирание).

Примеры:
  python sweep_ramp.py run --strategy plan_sweep_v2_decouple.json --case core130 --level coarse \\
         --geom-out ../geom/out/core130_coarse --label core130 --np 32 --time-limit 1320 --work W --out OUT
  python sweep_ramp.py run ... --mini-iter 3 --np 2 --time-limit 1500     # миниатюра: 3 итерации на шаг, гейт не держит
  python sweep_ramp.py merge --in ../aws/runs/<run_id> --out su2_kv3_sweep.json --run-id <run_id>
"""
from __future__ import annotations

import argparse
import glob
import json
import math
import os
import re
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
import char_points as CP  # noqa: E402

RAMP_FREQ = 5          # обновление рамп частоты и давления, итераций (run_case.py --ramp-freq, как в пробе)
STEP_MARGIN_S = 25.0   # запас на сводку шага после останова SU2 по сроку
MIN_STEP_S = 90.0      # меньше — шаг не начинается (SU2 не успеет прочитать сетку и сделать итерации)

# Гейты. Пороги приняты здесь (не из документации SU2): промежуточный — дисбаланс расхода <= 2 % в среднем за
# последние 20 строк, расход на входе меняется <= 1 % за 100 итераций, невязка плотности не растёт (наклон <= +0,05
# порядка на 100 итераций); итоговый (база) — дисбаланс <= 1 % за 50 строк, расход меняется <= 0,5 % за 200 итераций,
# п и КПД физичны; «успех» — ещё и невязка плотности ниже первой итерации цепочки на 3 порядка (задача владельца:
# 3–4 порядка). Строже — навык su2-turbo: 4 порядка и п, КПД меняются < 0,1 % за 500 итераций (converged_skill).
GATES = {
    "intermediate": dict(window=100, imb_rows=20, imb_tol=0.02, m_change_tol=0.01, slope_max=0.05),
    "final": dict(window=200, imb_rows=50, imb_tol=0.01, m_change_tol=0.005, slope_max=0.05,
                  rho_drop_min=3.0, rho_drop_skill=4.0, pi_eta_window=500, pi_eta_tol=1e-3),
}
PHYS = dict(pr_max_factor=1.3, eta_min=0.5, eta_max=1.0, su2_eff_min=50.0, su2_eff_max=100.0, prtt_rel_tol=0.05)
VERDICT_RANK = {"success": 6, "established": 5, "reached_base": 4, "stuck": 3, "timeout": 2, "terminated": 2,
                "diverged": 1, "error": 0}

# ----------------------------------------------------------------------------------------------- производный шаблон
HIST_ADD = ["CFL_NUMBER", "LINSOL", "MAX_RES"]   # группы истории SU2 v8.6.0 (CFlowCompOutput.cpp, CFlowOutput.cpp)
RIEMANN_LINE = ("MARKER_RIEMANN= (inlet, TOTAL_CONDITIONS_PT, ${P01}, ${T01}, 0.0, 0.0, 1.0, "
                "outlet, STATIC_PRESSURE, ${P_OUT}, 0.0, 0.0, 0.0, 0.0)")


def derive_template(text: str, transforms: list) -> tuple:
    """Шаблон run_case.py с правками по списку. hist — мягкая (строки нет — пропуск), riemann — жёсткая (строка
    MARKER_GILES должна встретиться ровно один раз). Возвращает (текст, применённые, пропущенные)."""
    applied, skipped = [], []
    for t in transforms:
        if t == "hist":
            pat = re.compile(r"^(HISTORY_OUTPUT=\s*\()([^)\n]*)(\)\s*)$", re.M)
            m = list(pat.finditer(text))
            if len(m) != 1:
                skipped.append(t)
                continue
            have = [x.strip() for x in m[0].group(2).split(",") if x.strip()]
            new = have + [g for g in HIST_ADD if g not in have]
            text = text[:m[0].start()] + m[0].group(1) + ", ".join(new) + m[0].group(3) + text[m[0].end():]
            applied.append(t)
        elif t == "riemann":
            pat = re.compile(r"^MARKER_GILES=.*$", re.M)
            m = list(pat.finditer(text))
            if len(m) != 1:
                raise SystemExit(f"правка riemann: строк MARKER_GILES в шаблоне {len(m)}, нужна одна — шаблон изменился")
            text = text[:m[0].start()] + RIEMANN_LINE + text[m[0].end():]
            applied.append(t)
        else:
            raise SystemExit(f"неизвестная правка шаблона: {t}")
    return text, applied, skipped


def cmd_rc(argv: list) -> int:
    """sweep_ramp.py rc [--transforms a,b] --tpl-dir D -- <ключи run_case.py>"""
    if "--" not in argv:
        raise SystemExit("rc: после своих ключей нужен -- и ключи run_case.py")
    i = argv.index("--")
    ap = argparse.ArgumentParser(prog="sweep_ramp.py rc")
    ap.add_argument("--transforms", default="")
    ap.add_argument("--tpl-dir", required=True)
    o = ap.parse_args(argv[:i])
    tr = [t for t in o.transforms.split(",") if t]
    txt, applied, skipped = derive_template(open(RC.TEMPLATE, encoding="utf-8").read(), tr)
    os.makedirs(o.tpl_dir, exist_ok=True)
    path = os.path.join(o.tpl_dir, "impeller_passage_" + ("_".join(applied) or "base") + ".cfg.tpl")
    with open(path, "w", encoding="utf-8") as f:
        f.write(txt)
    print(f"[rc] шаблон {path}; правки {applied}; пропущены {skipped}", flush=True)
    RC.TEMPLATE = path
    return RC.main(argv[i + 1:])


# ----------------------------------------------------------------------------------------------- стратегия
def throttle_p(frac: float, load: float, est: float, fs_p: float) -> float:
    """Давление на выходе на линии дросселя: при frac = 1 — load x оценка, при frac -> 0 — p_нач."""
    return fs_p + frac ** 2 * (load * est - fs_p)


def fs_pressure(d1: dict, init_mach: float | None = None, gamma: float = 1.4) -> float:
    """Статическое давление начального поля — как run_case.py (p01, T01 при M = INIT_MACH)."""
    m0 = float(RC.NUMERICS["INIT_MACH"] if init_mach is None else init_mach)
    t = d1["T01"] / (1 + 0.5 * (gamma - 1) * m0 ** 2)
    return d1["p01"] * (t / d1["T01"]) ** (gamma / (gamma - 1))


STEP_KEYS = {"name", "kind", "plan", "frac", "load", "ramp", "pramp", "iter", "muscl", "muscl_ramp", "cfl",
             "flow_scheme", "limiter", "gate", "hold", "note"}


def flow_scheme_of(strat: dict, st: dict) -> str:
    return str(st.get("flow_scheme") or strat.get("set", {}).get("FLOW_SCHEME") or RC.NUMERICS["FLOW_SCHEME"])


def validate_strategy(strat: dict) -> None:
    """Проверка стратегии до запуска (без SU2): поля шагов, гейт успевает увидеть окно после рамп, JST без MUSCL
    (SU2 v8.6.0: «Centered schemes do not use MUSCL reconstruction», CConfig.cpp), Riemann без рампы давления
    (SU2: «Outlet pressure ramp only implemented for NRBC», CFluidIteration.cpp)."""
    for k in ("id", "title", "steps"):
        if k not in strat:
            raise SystemExit(f"стратегия: нет поля {k}")
    if not strat["steps"]:
        raise SystemExit("стратегия: пустой список шагов")
    names = set()
    riemann = "riemann" in strat.get("transforms", [])
    for st in strat["steps"]:
        unknown = set(st) - STEP_KEYS
        if unknown:
            raise SystemExit(f"шаг {st.get('name')}: неизвестные поля {sorted(unknown)}")
        nm = st.get("name", "")
        if not re.fullmatch(r"[a-z0-9_]{1,24}", nm) or nm in names:
            raise SystemExit(f"шаг {nm!r}: имя [a-z0-9_]{{1,24}}, без повторов")
        names.add(nm)
        kind = st.get("kind", "stage")
        if kind == "plan":
            if st.get("plan") not in RC.PLANS or "load" not in st:
                raise SystemExit(f"шаг {nm}: kind plan — нужны plan (имя плана run_case.py) и load")
            continue
        for k in ("frac", "iter"):
            if k not in st:
                raise SystemExit(f"шаг {nm}: нет поля {k}")
        it = int(st["iter"])
        longest = max(int(st.get("ramp", 0)), int(st.get("pramp", 0)), int(st.get("muscl_ramp", 0)))
        gate = st.get("gate", "intermediate")
        if gate not in ("intermediate", "final", "none"):
            raise SystemExit(f"шаг {nm}: gate — intermediate, final или none")
        if gate != "none" and longest > it - GATES[gate]["window"]:
            raise SystemExit(f"шаг {nm}: рампа {longest} итераций не оставляет окна гейта {GATES[gate]['window']} "
                             f"до конца шага ({it})")
        if st.get("muscl_ramp") and not st.get("muscl"):
            raise SystemExit(f"шаг {nm}: muscl_ramp без muscl")
        if flow_scheme_of(strat, st) == "JST" and st.get("muscl"):
            raise SystemExit(f"шаг {nm}: JST с MUSCL_FLOW= YES SU2 v8.6.0 не запускает — muscl: false")
        if riemann and int(st.get("pramp", 0)) > 0:
            raise SystemExit(f"шаг {nm}: с выходом Riemann рампы давления нет — pramp 0, давление ступенью")
    if strat["steps"][-1].get("gate", "intermediate") != "final" and strat["steps"][-1].get("kind") != "plan":
        raise SystemExit("последний шаг стратегии — с gate final (база)")


def mini_step(st: dict, n: int) -> dict:
    """Миниатюра шага: n итераций, рампы короче n, без порций удержания (кроме итогового шага — одна)."""
    s = dict(st)
    if s.get("kind") == "plan":
        s["iter_scale"] = n / 200.0
        s["hold"] = dict(s.get("hold", {}), iter=n, max=1)
        return s
    s["iter"] = n
    for k in ("ramp", "pramp", "muscl_ramp"):
        if int(s.get(k, 0)) > 0:
            s[k] = max(1, n - 1)
    s["hold"] = dict(s.get("hold", {}), iter=n, max=1 if s.get("gate") == "final" else 0)
    return s


def cfl_halved(set_fields: dict) -> dict:
    """Верхний предел адаптивного CFL вдвое (откат), не ниже нижнего."""
    out = dict(set_fields)
    vals = [float(x) for x in str(out.get("CFL_ADAPT_PARAM", RC.NUMERICS["CFL_ADAPT_PARAM"])).split(",")]
    if len(vals) >= 4:
        vals[3] = max(vals[2], vals[3] / 2.0)
    out["CFL_ADAPT_PARAM"] = ", ".join(f"{v:g}" for v in vals)
    return out


def _stretch(s: dict, k: int) -> dict:
    """Рампы шага в k раз длиннее, итераций — на столько же больше (окно гейта после рамп сохраняется)."""
    s = dict(s)
    add = 0
    for f in ("ramp", "pramp", "muscl_ramp"):
        if int(s.get(f, 0)) > 0:
            add = max(add, int(s[f]) * (k - 1))
            s[f] = int(s[f]) * k
    s["iter"] = int(s["iter"]) + add
    return s


def backoff_steps(st: dict, state: dict, strat_set: dict) -> list:
    """Откат после NaN или непройденного гейта (состояние — последнее пройденное, его рестарт цел).
    1-я попытка: тот же шаг, рампы вдвое длиннее, верхний предел адаптивного CFL вдвое ниже.
    2-я: шаг делится пополам по частоте и нагрузке (каждая половина — с исходной длиной рамп, то есть вдвое
    медленнее), предел CFL тоже вдвое ниже. Пусто — откатов больше нет."""
    attempt = int(st.get("_attempt", 0))
    orig = st.get("_orig", st)
    if st.get("kind") == "plan" or attempt >= 2:
        return []
    low_cfl = cfl_halved(dict(strat_set, **orig.get("_set", {})))
    cfl_max = float(low_cfl["CFL_ADAPT_PARAM"].split(",")[3])

    def finish(s, name, att):
        s["_set"], s["_attempt"], s["_orig"], s["name"] = low_cfl, att, orig, name[:24]
        if s.get("cfl"):
            s["cfl"] = min(float(s["cfl"]), cfl_max)
        return s
    if attempt == 0:
        return [finish(_stretch(orig, 2), f"{orig['name']}_r1", 1)]
    f0, f1 = float(state["frac"]), float(orig["frac"])
    l1, l0 = orig.get("load"), state.get("load")
    same_load = l1 is None or (l0 is not None and abs(l0 - l1) < 1e-9)
    if abs(f1 - f0) < 1e-9 and same_load:
        return []   # делить нечего (шаг без смены частоты и нагрузки — например, переход ко второму порядку)
    mid = dict(orig, frac=round((f0 + f1) / 2, 4), gate="intermediate",
               hold=dict(orig.get("hold", {}), max=max(2, int(orig.get("hold", {}).get("max", 2)))))
    if l1 is not None:
        mid["load"] = round(((l0 if l0 is not None else l1) + l1) / 2, 4)
    return [finish(mid, f"{orig['name']}_m2", 2), finish(dict(orig), f"{orig['name']}_r2", 2)]


# ----------------------------------------------------------------------------------------------- метрики и гейт
def _fin(x) -> bool:
    return isinstance(x, float) and math.isfinite(x)


def slope_per_100(vals: list) -> float | None:
    if len(vals) < 10 or not all(_fin(v) for v in vals):
        return None
    xs = range(len(vals))
    mx, my = (len(vals) - 1) / 2, sum(vals) / len(vals)
    sxx = sum((x - mx) ** 2 for x in xs)
    return 100 * sum((x - mx) * (y - my) for x, y in zip(xs, vals)) / sxx


def merge_rows(buf: dict | None, data: dict) -> dict:
    """Буфер истории состояния: порции удержания дописываются к шагу (колонки — пересечение)."""
    if not buf:
        return {k: list(v) for k, v in data.items()}
    keys = set(buf) & set(data)
    return {k: buf[k] + data[k] for k in keys}


def gate_metrics(d: dict, gamma: float, window: int, imb_rows: int, pi_eta_window: int = 500) -> dict:
    m_in, m_out = d.get("Avg_Massflow(inlet)", []), d.get("Avg_Massflow(outlet)", [])
    n = min(len(m_in), len(m_out))
    out = {"rows": n}
    if n == 0:
        return out
    a, b = m_in[-1], m_out[-1]
    out.update(m_in_last=a, m_out_last=b, sign_ok=bool(_fin(a) and _fin(b) and a > 0 > b))
    imb = [abs(x + y) / abs(x) for x, y in zip(m_in[n - imb_rows:n], m_out[n - imb_rows:n])
           if _fin(x) and _fin(y) and x != 0]
    out["imbalance_last"] = abs(a + b) / abs(a) if _fin(a) and _fin(b) and a != 0 else None
    out["imbalance_mean"] = sum(imb) / len(imb) if imb else None
    out["m_in_change"] = CP.window_change(m_in[:n], window)
    rho = d.get("rms[Rho]", [])
    if rho:
        out["rms_rho_last"] = rho[-1]
        out["rms_rho_slope_per_100"] = slope_per_100(rho[-window:]) if len(rho) >= window else None
    p1, p2 = d.get("Avg_TotalPress(inlet)", []), d.get("Avg_TotalPress(outlet)", [])
    t1, t2 = d.get("Avg_TotalTemp(inlet)", []), d.get("Avg_TotalTemp(outlet)", [])
    k = min(len(p1), len(p2), len(t1), len(t2))
    pr, eta = [], []
    for i in range(k):
        try:
            x, y = p2[i] / p1[i], t2[i] / t1[i]
            pr.append(x)
            eta.append((x ** ((gamma - 1) / gamma) - 1) / (y - 1) if y > 1 else float("nan"))
        except (ZeroDivisionError, ValueError, OverflowError):
            pr.append(float("nan"))
            eta.append(float("nan"))
    if k:
        out["PR_tt"], out["eta_tt"] = pr[-1], eta[-1]
        out["TR_tt"] = t2[-1] / t1[-1] if t1[-1] else None
        out["PR_change"] = CP.window_change(pr, pi_eta_window)
        out["eta_change"] = CP.window_change(eta, pi_eta_window)
    for col, key in (("PRTT", "PRTT_su2"), ("TotTotEff", "TotTotEff_su2"), ("Avg CFL", "avg_cfl"),
                     ("Max CFL", "max_cfl"), ("Linear_Solver_Iterations", "linsol_iter"), ("LinSolRes", "linsol_res")):
        if d.get(col):
            out[key] = d[col][-1]
    last = [v[-1] for v in d.values() if v]
    out["nan_in_last_row"] = any(isinstance(v, float) and math.isnan(v) for v in last)
    return out


def physical(m: dict, pr_1d: float) -> dict:
    pr, tr, eta = m.get("PR_tt"), m.get("TR_tt"), m.get("eta_tt")
    c = {"pr": bool(_fin(pr) and 1.0 < pr <= PHYS["pr_max_factor"] * pr_1d),
         "tr": bool(_fin(tr) and tr > 1.0),
         "eta": bool(_fin(eta) and PHYS["eta_min"] <= eta < PHYS["eta_max"])}
    su2e = m.get("TotTotEff_su2")
    c["su2_eff"] = bool(_fin(su2e) and PHYS["su2_eff_min"] <= su2e < PHYS["su2_eff_max"]) if su2e is not None else None
    prtt = m.get("PRTT_su2")
    c["prtt_consistent"] = (bool(_fin(prtt) and _fin(pr) and abs(prtt / pr - 1) <= PHYS["prtt_rel_tol"])
                            if prtt is not None else None)
    c["ok"] = all(v for v in c.values() if v is not None)
    return c


def gate_check(m: dict, kind: str, rho_first: float | None, pr_1d: float) -> dict:
    """Решение гейта по метрикам. intermediate/final: pass — можно дальше; final ещё success и converged_skill."""
    if kind == "none":
        return {"pass": True, "reasons": []}
    g = GATES[kind]
    r = []
    if m.get("rows", 0) < g["window"] + 1:
        r.append(f"строк {m.get('rows', 0)} < окна {g['window'] + 1}")
    if not m.get("sign_ok"):
        r.append("обратный ток или нет расхода (нужно вход > 0 > выход)")
    imb = m.get("imbalance_mean")
    if imb is None or imb > g["imb_tol"]:
        r.append(f"дисбаланс расхода {imb if imb is None else round(imb, 4)} > {g['imb_tol']}")
    ch = m.get("m_in_change")
    if ch is None or ch > g["m_change_tol"]:
        r.append(f"расход не установился: изменение {ch if ch is None else round(ch, 4)} > {g['m_change_tol']} "
                 f"за {g['window']} итераций")
    sl = m.get("rms_rho_slope_per_100")
    if sl is not None and sl > g["slope_max"]:
        r.append(f"невязка плотности растёт: {sl:+.3f} порядка/100")
    out = {"pass": not r, "reasons": r}
    if kind == "final":
        ph = physical(m, pr_1d)
        out["physical"] = ph
        if not ph["ok"]:
            out["pass"] = False
            out["reasons"].append("п или КПД нефизичны: " + ", ".join(k for k, v in ph.items() if v is False))
        drop = (rho_first - m["rms_rho_last"]) if (_fin(rho_first) and _fin(m.get("rms_rho_last"))) else None
        out["rho_drop_orders"] = drop
        out["established"] = out["pass"]
        out["success"] = bool(out["pass"] and drop is not None and drop >= g["rho_drop_min"])
        out["converged_skill"] = bool(out["success"] and drop >= g["rho_drop_skill"]
                                      and _fin(m.get("PR_change")) and _fin(m.get("eta_change"))
                                      and m["PR_change"] < g["pi_eta_tol"] and m["eta_change"] < g["pi_eta_tol"])
    return out


def ramp_value(i: int, ini: float, fin: float, length: int, freq: int) -> float:
    """Значение рампы SU2 на итерации i (CFluidIteration::UpdateRamps): обновление при i % freq == 0 и i < length,
    либо i == length; между обновлениями — прежнее."""
    if length <= 0 or i >= length:
        return fin
    j = (i // freq) * freq
    return ini + j * (fin - ini) / length


def analyse_divergence(d: dict, rec: dict, omega_nom: float) -> dict:
    """Что было перед NaN: последняя строка истории, рост невязок за 10 строк (какое уравнение первым), CFL,
    расход, частота и давление по рампам на итерации останова."""
    n = len(d.get("Inner_Iter", []))
    out = {"rows": n, "su2_error": rec.get("su2_error")}
    if n == 0:
        return out
    it_nan = int(d["Inner_Iter"][-1]) + 1
    out["iter_of_stop"] = it_nan
    fr0, fr1 = rec.get("frac_from", rec["frac"]), rec["frac"]
    p0, p1 = rec.get("p_from", rec["p_out"]), rec["p_out"]
    out["frac_at_stop"] = ramp_value(it_nan, fr0, fr1, int(rec.get("ramp", 0)), RAMP_FREQ)
    out["p_out_at_stop"] = ramp_value(it_nan, p0, p1, int(rec.get("pramp", 0)), RAMP_FREQ) if rec.get("pramp") else p1
    out["omega_at_stop_rad_s"] = out["frac_at_stop"] * omega_nom
    grow = {}
    for h, v in d.items():
        if (h.startswith("rms[") or h.startswith("max[")) and len(v) >= 2:
            seg = [x for x in v[-11:] if _fin(x)]
            if len(seg) >= 2:
                grow[h] = seg[-1] - min(seg)
    out["residual_growth_last10"] = dict(sorted(grow.items(), key=lambda kv: -kv[1])[:8])
    out["fastest_growing"] = max(grow, key=grow.get) if grow else None
    for col, key in (("Avg CFL", "avg_cfl_last"), ("Max CFL", "max_cfl_last"), ("Avg_Massflow(inlet)", "m_in_last"),
                     ("Avg_Massflow(outlet)", "m_out_last"), ("Linear_Solver_Iterations", "linsol_iter_last")):
        if d.get(col):
            out[key] = d[col][-1]
    if "m_in_last" in out and "m_out_last" in out:
        out["backflow_at_outlet"] = bool(out["m_out_last"] > 0)
    return out


# ----------------------------------------------------------------------------------------------- запуск шага
_CUR = {"proc": None, "signal": None}


def _forward(signum, frame):  # noqa: ARG001
    _CUR["signal"] = signum
    p = _CUR["proc"]
    if p is not None and p.poll() is None:
        try:
            p.send_signal(signal.SIGTERM)
        except ProcessLookupError:
            pass


def step_argv(ctx: dict, strat: dict, st: dict, plan_ref: str, p_out: float, time_limit: float, out_dir: str,
              state: dict, iter_scale: float = 1.0) -> list:
    tr = ",".join(strat.get("transforms", []))
    cmd = [sys.executable, "-u", os.path.abspath(__file__), "rc", "--transforms", tr, "--tpl-dir",
           os.path.join(ctx["work"], "_tpl"), "--",
           "--case", ctx["case"], "--level", ctx["level"], "--geom-out", ctx["geom_out"], "--plan", plan_ref,
           "--p-out", f"{p_out:.6g}", "--np", str(ctx["np"]), f"--mpi-args={ctx['mpi_args']}",
           "--shroud", ctx["shroud"], "--vaneless", ctx["vaneless"], "--ramp-freq", str(RAMP_FREQ),
           "--time-limit", f"{max(time_limit, 1.0):.0f}", "--out", out_dir, "--iter-scale", f"{iter_scale:g}"]
    fields = dict(strat.get("set", {}))
    fields.update(st.get("_set", {}))
    for k, v in sorted(fields.items()):
        cmd += ["--set", f"{k}={v}"]
    extra = list(strat.get("extra", []))
    if st.get("muscl_ramp"):
        extra += ["RAMP_MUSCL= YES", f"RAMP_MUSCL_COEFF= (0, 1, {int(st['muscl_ramp'])})", "KIND_MUSCL_RAMP= SMOOTH_FUNCTION"]
    for x in extra:
        cmd.append(f"--extra={x}")
    if strat.get("turb"):
        cmd += ["--turb", strat["turb"]]
    if ctx.get("su2_src"):
        cmd += ["--su2-src", ctx["su2_src"]]
    if state.get("restart"):
        cmd += ["--restart-from", state["restart"], "--restart-state", f"{state['frac']:.6g},{state['p']:.6g}"]
    return cmd


def run_step_subprocess(ctx: dict, argv: list, out_dir: str) -> dict:
    """Запуск шага (rc -> run_case.py), разбор run_summary.json. Возвращает итог для контроллера."""
    os.makedirs(out_dir, exist_ok=True)
    t0 = time.monotonic()
    with open(os.path.join(out_dir, "run_case.out"), "w", encoding="utf-8") as lf:
        p = subprocess.Popen(argv, stdout=lf, stderr=subprocess.STDOUT)
        _CUR["proc"] = p
        rc = p.wait()
        _CUR["proc"] = None
    res = {"returncode": rc, "seconds": round(time.monotonic() - t0, 1)}
    sp = os.path.join(out_dir, "run_summary.json")
    if not os.path.isfile(sp):
        res["status"] = "error"
        res["error"] = "нет run_summary.json (см. run_case.out)"
        return res
    s = json.load(open(sp, encoding="utf-8"))
    stages = s.get("stages", [])
    res["status"] = stages[-1].get("status") if stages else "error"
    res["histories"] = [os.path.join(out_dir, x["history_file"]) for x in stages
                        if os.path.isfile(os.path.join(out_dir, x["history_file"]))]
    res["stage_names"] = [x["name"] for x in stages]
    last = stages[-1] if stages else {}
    res["log"] = last.get("log", {})
    res["su2_version"] = s.get("su2", {}).get("version")
    res["key_check"] = s.get("key_check")
    res["restart"] = (os.path.join(out_dir, f"restart_{last['name']}.dat")
                      if last and os.path.isfile(os.path.join(out_dir, f"restart_{last['name']}.dat")) else None)
    return res


# ----------------------------------------------------------------------------------------------- контроллер
def run_strategy(strat: dict, ctx: dict, runner=None, clock=time.monotonic) -> dict:
    """Шаги стратегии по очереди с гейтами, удержанием и откатом. runner(st, state, chunk) -> итог шага (в тестах
    подменяется; по умолчанию — rc -> run_case.py). Возвращает сводку."""
    validate_strategy(strat)
    mini = ctx.get("mini_iter")
    enforce = not mini
    est, fs_p, pr_1d = ctx["est"], ctx["fs_p"], ctx["pr_1d"]
    state = {"frac": 0.0, "p": fs_p, "load": None, "restart": None}
    summ = {"steps": [], "verdict": None, "verdict_step": None, "verdict_detail": None, "points_full_speed": [],
            "divergences": [], "rho_first": None, "backoffs": 0, "iters_total": 0, "time_to_base_s": None,
            "iters_to_base": None}
    t_start = clock()
    backoff_cfg = strat.get("backoff") or {}
    queue = [mini_step(s, mini) if mini else dict(s) for s in strat["steps"]]

    def default_runner(st, state_, chunk):
        nm = st["_dir"]
        out_dir = os.path.join(ctx["work"], nm)
        if st.get("kind") == "plan" and chunk == 0:
            plan_ref, iscale = st["plan"], st.get("iter_scale", 1.0)
        else:
            hold = st["hold"]
            stage = {"name": nm, "frac": st["frac"], "pfrac": 1.0,
                     "ramp": int(st.get("ramp", 0)) if chunk == 0 else 0,
                     "pramp": int(st.get("pramp", 0)) if chunk == 0 else 0,
                     "iter": int(st["iter"]) if chunk == 0 else int(hold["iter"]),
                     "muscl": bool(st.get("muscl")) if chunk == 0 else bool(hold.get("muscl", st.get("muscl")))}
            cfl = st.get("cfl") if chunk == 0 else hold.get("cfl", st.get("cfl"))
            if cfl:
                stage["cfl"] = cfl
            for k in ("flow_scheme", "limiter"):
                if st.get(k):
                    stage[k] = st[k]
            os.makedirs(os.path.join(ctx["work"], "_plans"), exist_ok=True)
            plan_ref, iscale = os.path.join(ctx["work"], "_plans", f"{nm}.json"), 1.0
            with open(plan_ref, "w", encoding="utf-8") as f:
                json.dump([stage], f)
        st_run = dict(st)
        if chunk > 0:
            st_run.pop("muscl_ramp", None)   # рампа MUSCL — только в первой порции шага
        left = ctx["deadline"] - clock() - STEP_MARGIN_S
        argv = step_argv(ctx, strat, st_run, plan_ref, st["_p_out"], left, out_dir, state_, iscale)
        res = run_step_subprocess(ctx, argv, out_dir)
        res["argv"] = argv[argv.index("--") + 1:]
        return res

    run = runner or default_runner
    idx = 0
    while idx < len(queue):
        st = queue[idx]
        kind = st.get("kind", "stage")
        gate_kind = st.get("gate", "final" if kind == "plan" else "intermediate")
        if kind == "plan":
            st["frac"], st["_p_out"] = 1.0, st["load"] * est
            st.setdefault("hold", {}).setdefault("muscl", True)
        else:
            st["_p_out"] = (state["p"] if st.get("load") is None
                            else throttle_p(float(st["frac"]), float(st["load"]), est, fs_p))
        hold = dict({"iter": 300, "max": 3}, **st.get("hold", {}))
        st["hold"] = hold
        buf, chunk, decided, last_met = None, 0, None, {}
        while True:
            if _CUR["signal"]:
                decided = ("terminated", "сигнал останова")
                break
            left = ctx["deadline"] - clock()
            if left < MIN_STEP_S:
                decided = ("timeout", f"до срока {left:.0f} с — порция не начата")
                break
            st["_dir"] = f"{idx:02d}_{st['name']}" + (f"_h{chunk}" if chunk else "")
            res = run(st, state, chunk)
            data_last, first_rho, n_rows = {}, None, 0
            for hp in res.get("histories", []):
                _, dd = RC.read_history(hp)
                n_rows += len(dd.get("Inner_Iter", []))
                if first_rho is None and dd.get("rms[Rho]"):
                    first_rho = dd["rms[Rho]"][0]
                data_last = dd                      # гейт — по последней стадии шага
            if summ["rho_first"] is None and _fin(first_rho):
                summ["rho_first"] = first_rho
            summ["iters_total"] += n_rows
            if data_last:
                buf = merge_rows(buf, data_last)
            g = GATES.get(gate_kind, GATES["intermediate"])
            met = gate_metrics(buf or {}, 1.4, g["window"], g["imb_rows"], GATES["final"]["pi_eta_window"])
            last_met = met
            gate = gate_check(met, gate_kind, summ["rho_first"], pr_1d)
            log = res.get("log") or {}
            rec = {"step": st["name"], "dir": st["_dir"], "chunk": chunk, "kind": kind, "attempt": st.get("_attempt", 0),
                   "frac": st["frac"], "frac_from": state["frac"], "p_out": st["_p_out"], "p_from": state["p"],
                   "load": st.get("load"), "ramp": int(st.get("ramp", 0)) if chunk == 0 else 0,
                   "pramp": int(st.get("pramp", 0)) if chunk == 0 else 0,
                   "iter": st.get("iter") if chunk == 0 else hold["iter"],
                   "muscl": (bool(st.get("muscl")) or kind == "plan") if chunk == 0 else bool(hold.get("muscl", st.get("muscl"))),
                   "muscl_ramp": st.get("muscl_ramp") if chunk == 0 else None,
                   "cfl": st.get("cfl") if chunk == 0 else hold.get("cfl", st.get("cfl")),
                   "set": dict(strat.get("set", {}), **st.get("_set", {})),
                   "status": res.get("status"), "seconds": res.get("seconds"), "rows": n_rows,
                   "su2_error": log.get("error"), "metrics": met, "gate": gate,
                   "periodic": {k: log.get(k) for k in ("periodic_matched_points", "periodic_max_distance", "bad_match")}}
            if res.get("su2_version"):
                summ["su2_version"], summ["key_check"] = res["su2_version"], res.get("key_check")
            summ["steps"].append(rec)
            status = res.get("status")
            if status in ("diverged", "error") or met.get("nan_in_last_row"):
                rec["decision"] = "failed"
                summ["divergences"].append(dict(analyse_divergence(data_last, rec, ctx.get("omega_nom", 0.0)),
                                                step=st["name"], dir=rec["dir"], status=status))
                decided = ("error" if status == "error" else "diverged", f"{status} на шаге {st['name']}")
            elif status in ("timeout", "terminated", "no_time") or _CUR["signal"]:
                rec["decision"] = "stopped"
                if gate_kind == "final" and gate.get("established"):
                    decided = ("established", "срок на базе")
                else:
                    decided = ("terminated" if status == "terminated" or _CUR["signal"] else "timeout",
                               f"{status} на шаге {st['name']}")
            else:
                # порция дошла до конца: состояние зафиксировано, её рестарт — начало следующей
                state = {"frac": float(st["frac"]), "p": st["_p_out"], "restart": res.get("restart"),
                         "load": st.get("load") if st.get("load") is not None else state.get("load")}
                if abs(state["frac"] - 1.0) < 1e-9 and gate_kind != "none" and gate.get("pass"):
                    summ["points_full_speed"].append(
                        {"step": st["name"], "chunk": chunk, "p_out": state["p"], "load": state["load"],
                         "muscl": rec["muscl"], **{k: met.get(k) for k in ("m_in_last", "imbalance_mean", "PR_tt",
                                                                           "eta_tt", "PRTT_su2", "TotTotEff_su2",
                                                                           "rms_rho_last")}})
                if gate_kind == "final" and summ["time_to_base_s"] is None:
                    summ["time_to_base_s"], summ["iters_to_base"] = round(clock() - t_start, 1), summ["iters_total"]
                if gate["pass"] and gate_kind != "final":
                    rec["decision"] = "advance"
                elif gate_kind == "final" and gate.get("success"):
                    rec["decision"] = "success"
                    decided = ("success", "критерии базы выполнены")
                elif chunk < hold["max"]:
                    rec["decision"] = "hold"
                elif not enforce:
                    rec["decision"] = "advance_mini"
                    if gate_kind == "final":
                        decided = ("established" if gate.get("established") else "reached_base", "миниатюра")
                else:
                    rec["decision"] = "gate_failed"
                    why = "; ".join(gate["reasons"])
                    decided = (("established" if gate.get("established") else "reached_base", "удержание исчерпано: " + why)
                               if gate_kind == "final" else ("stuck", f"гейт шага {st['name']}: {why}"))
            if ctx.get("on_step"):
                ctx["on_step"](rec, res)
            if rec["decision"] == "hold":
                chunk += 1
                continue
            break
        backflow_stuck = decided and decided[0] == "reached_base" and not last_met.get("sign_ok")
        if decided and (decided[0] in ("diverged", "stuck") or backflow_stuck) and backoff_cfg and enforce \
                and summ["backoffs"] < int(backoff_cfg.get("max_total", 0)):
            new = backoff_steps(st, state, strat.get("set", {}))
            if new:
                summ["backoffs"] += 1
                summ["steps"][-1]["decision"] += "_backoff"
                queue[idx + 1:idx + 1] = new
                idx += 1
                continue
        if decided:
            summ["verdict"], summ["verdict_detail"], summ["verdict_step"] = decided[0], decided[1], st["name"]
            break
        idx += 1
    if summ["verdict"] is None:
        summ["verdict"], summ["verdict_detail"] = "reached_base", "шаги кончились"
    summ["divergence"] = summ["divergences"][-1] if summ["verdict"] in ("diverged", "error") and summ["divergences"] else None
    summ["final_state"] = {"frac": state["frac"], "p_out": state["p"], "load": state.get("load"),
                           "restart": state.get("restart")}
    summ["seconds_total"] = round(clock() - t_start, 1)
    fin = [r for r in summ["steps"] if r["gate"].get("established") is not None]
    summ["criteria"] = (({k: fin[-1]["gate"].get(k) for k in ("established", "success", "converged_skill",
                                                             "rho_drop_orders", "physical", "reasons")}
                         | {k: fin[-1]["metrics"].get(k) for k in ("imbalance_mean", "m_in_change", "sign_ok", "PR_tt",
                                                                   "eta_tt", "TR_tt", "PRTT_su2", "TotTotEff_su2",
                                                                   "rms_rho_last", "PR_change", "eta_change",
                                                                   "m_in_last", "avg_cfl")}) if fin else None)
    return summ


# ----------------------------------------------------------------------------------------------- run
def keep_step_files(work_dir: str, out_dir: str) -> list:
    got = CP.keep_files(work_dir, out_dir, None)
    if os.path.isfile(os.path.join(work_dir, "run_case.out")):
        shutil.copy2(os.path.join(work_dir, "run_case.out"), os.path.join(out_dir, "run_case.out"))
        got.append("run_case.out")
    CP.drop_big(work_dir)
    return got


def cmd_run(a) -> int:
    t0 = time.monotonic()
    strat = json.load(open(a.strategy, encoding="utf-8"))
    validate_strategy(strat)
    su2_mesh, js = RC.mesh_paths(a.case, a.level, a.geom_out)
    if not (os.path.isfile(su2_mesh) and os.path.isfile(js)):
        raise SystemExit(f"нет сетки {su2_mesh} или {js}: сетку строит команда задачи")
    mj = json.load(open(js, encoding="utf-8"))
    est = CP.p_out_estimate(mj)
    d1 = mj["d1"]
    uptime = CP.read_uptime(a.uptime_file) if a.job_max_min else 0.0
    budget = a.time_limit
    if a.job_max_min:
        budget = min(budget, CP.time_budget(a.job_max_min, None, uptime, a.reserve_s))
    os.makedirs(a.work, exist_ok=True)
    os.makedirs(a.out, exist_ok=True)
    for sig in (signal.SIGTERM, signal.SIGINT, signal.SIGHUP):
        signal.signal(sig, _forward)
    fs_p = fs_pressure(d1, float(strat.get("set", {}).get("INIT_MACH", RC.NUMERICS["INIT_MACH"])))
    ctx = {"case": a.case, "level": a.level, "geom_out": os.path.abspath(a.geom_out), "np": a.np, "mpi_args": a.mpi_args,
           "shroud": a.shroud, "vaneless": a.vaneless, "su2_src": a.su2_src, "work": os.path.abspath(a.work),
           "deadline": t0 + budget, "est": est["p_out_estimate"], "fs_p": fs_p, "pr_1d": float(d1["PR"]),
           "omega_nom": float(mj["su2"]["rotation_rate_rad_s"][2]), "mini_iter": a.mini_iter}
    head = {"created_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"), "label": a.label or a.case,
            "case": a.case, "level": a.level, "job": os.environ.get("GTU_JOB"), "run_id": os.environ.get("GTU_RUN_ID"),
            "itype": os.environ.get("GTU_ITYPE"), "az": os.environ.get("GTU_AZ_MACHINE") or os.environ.get("GTU_AZ"),
            "strategy": {k: strat.get(k) for k in ("id", "title", "axis", "mechanism", "hypothesis")}
            | {"file": os.path.basename(a.strategy)},
            "mesh": {"file": os.path.basename(su2_mesh), "header": RC.read_mesh_header(su2_mesh),
                     "double_sided": est["double_sided"]},
            "d1_reference": {k: d1.get(k) for k in ("m", "m_side", "PR", "rpm", "D2", "U2", "b2", "T01", "p01",
                                                    "eta_comp_stage")},
            "estimate": est, "fs_p": fs_p, "p_base_target": None, "np": a.np, "mini_iter": a.mini_iter,
            "budget": {"time_s": round(budget, 1), "time_limit_s": a.time_limit, "job_max_min": a.job_max_min,
                       "uptime_at_start_s": round(uptime, 1), "reserve_s": a.reserve_s},
            "gates": GATES, "physical_bounds": PHYS, "ramp_freq": RAMP_FREQ,
            "throttle_line": "p(w) = p_нач + w^2 (доля x оценка - p_нач); база — доля шага gate final"}
    finals = [s for s in strat["steps"] if s.get("gate") == "final" or s.get("kind") == "plan"]
    if finals:
        head["p_base_target"] = finals[-1]["load"] * est["p_out_estimate"]
        head["base_load"] = finals[-1]["load"]

    def save(summ):
        doc = dict(head, **summ)
        with open(os.path.join(a.out, "sweep_summary.json"), "w", encoding="utf-8") as f:
            json.dump(doc, f, ensure_ascii=False, indent=1)

    partial = {"steps": [], "verdict": "running"}

    def on_step(rec, res):
        rec["files"] = keep_step_files(os.path.join(ctx["work"], rec["dir"]), os.path.join(a.out, rec["dir"]))
        partial["steps"].append(rec)
        print(f"[sweep] {head['label']} {rec['dir']}: {rec['status']} {rec['seconds']} с, строк {rec['rows']}, "
              f"дисбаланс {rec['metrics'].get('imbalance_mean')}, гейт {rec['gate'].get('pass')} -> {rec.get('decision')}",
              flush=True)
        save(partial)

    ctx["on_step"] = on_step
    print(f"[sweep] {head['label']}: стратегия {strat['id']}, оценка p_вых {est['p_out_estimate']:.0f} Па, база "
          f"{head['p_base_target']:.0f} Па, p_нач {fs_p:.0f} Па, срок {budget:.0f} с, np {a.np}"
          + (f", миниатюра {a.mini_iter} итераций" if a.mini_iter else ""), flush=True)
    summ = run_strategy(strat, ctx)
    # рестарты: итоговый (база) и последний записанный шагом, где был NaN (поле <= WRT_FREQ итераций до NaN)
    fr = summ["final_state"].get("restart")
    if fr and os.path.isfile(fr) and summ["verdict"] in ("success", "established", "reached_base", "timeout"):
        shutil.copy2(fr, os.path.join(a.out, "final_restart.dat"))
        summ["final_restart"] = "final_restart.dat"
    if summ.get("divergence"):
        dv = summ["divergence"]
        cand = sorted(glob.glob(os.path.join(ctx["work"], dv["dir"], "restart_*.dat")))
        cand = [c for c in cand if not c.endswith("restart_in.dat")]
        if cand:
            shutil.copy2(cand[-1], os.path.join(a.out, dv["dir"], "restart_before_nan.dat"))
            dv["restart_before_nan"] = f"{dv['dir']}/restart_before_nan.dat"
    summ["t_total_s"] = round(time.monotonic() - t0, 1)
    save(summ)
    c = summ.get("criteria") or {}
    print(f"[sweep] {head['label']} {strat['id']}: {summ['verdict']} ({summ.get('verdict_detail')}); итераций "
          f"{summ['iters_total']}, {summ['seconds_total']} с; дисбаланс {c.get('imbalance_mean')}, PR {c.get('PR_tt')}, "
          f"eta {c.get('eta_tt')}, порядков {c.get('rho_drop_orders')}", flush=True)
    return 0 if summ["verdict"] not in ("error", "terminated") else 5


# ----------------------------------------------------------------------------------------------- merge
def score(w: dict) -> tuple:
    c = w.get("criteria") or {}
    imb = c.get("imbalance_mean")
    return (VERDICT_RANK.get(w.get("verdict"), 0), w.get("steps_passed", 0),
            -(imb if _fin(imb) else 9.0), c.get("rho_drop_orders") or -9.0)


def merge(paths: list, run_id: str | None = None) -> dict:
    strategies = {}
    for p in sorted(paths):
        s = json.load(open(p, encoding="utf-8"))
        sid = (s.get("strategy") or {}).get("id") or "?"
        g = strategies.setdefault(sid, {"id": sid, "title": s["strategy"].get("title"), "axis": s["strategy"].get("axis"),
                                        "mechanism": s["strategy"].get("mechanism"), "wheels": {}})
        lab = s.get("label")
        if lab in g["wheels"]:
            raise ValueError(f"{sid}/{lab}: две сводки ({p})")
        passed = sum(1 for r in s.get("steps", []) if r.get("decision") in ("advance", "success"))
        g["wheels"][lab] = {"file": os.path.relpath(p), "job": s.get("job"), "verdict": s.get("verdict"),
                            "verdict_step": s.get("verdict_step"), "verdict_detail": s.get("verdict_detail"),
                            "criteria": s.get("criteria"), "iters_total": s.get("iters_total"),
                            "seconds_total": s.get("seconds_total"), "time_to_base_s": s.get("time_to_base_s"),
                            "iters_to_base": s.get("iters_to_base"), "steps_passed": passed,
                            "steps_run": len(s.get("steps", [])), "backoffs": s.get("backoffs"),
                            "points_full_speed": s.get("points_full_speed"), "divergence": s.get("divergence"),
                            "mesh_nodes": (s.get("mesh") or {}).get("header", {}).get("NPOIN"),
                            "p_base_target": s.get("p_base_target"), "mini_iter": s.get("mini_iter"),
                            "itype": s.get("itype"), "su2_version": s.get("su2_version")}
    for g in strategies.values():
        g["score"] = [score(w) for w in g["wheels"].values()]
        g["success_wheels"] = sorted(k for k, w in g["wheels"].items() if w["verdict"] == "success")
    labels = sorted({lab for g in strategies.values() for lab in g["wheels"]})
    best = {}
    for lab in labels:
        cand = [(score(g["wheels"][lab]), sid) for sid, g in strategies.items() if lab in g["wheels"]]
        best[lab] = max(cand)[1] if cand else None
    ranking = sorted(strategies, key=lambda sid: (sum(x[0] for x in strategies[sid]["score"]),
                                                  sum(x[1] for x in strategies[sid]["score"]),
                                                  sum(x[2] for x in strategies[sid]["score"])), reverse=True)
    return {"about": ("Перебор стратегий разгона базовой точки SU2 (sweep_ramp.py на AWS): колесо KV-3 ступень 1 и "
                      "core130, сетка coarse. Вердикт колеса: success — дисбаланс расхода <= 1 %, расход установился, "
                      "невязка плотности упала >= 3 порядков от первой итерации, п и КПД физичны; established — всё, "
                      "кроме 3 порядков; reached_base — база достигнута без NaN, гейт не пройден; stuck — гейт "
                      "промежуточного шага не пройден; diverged — NaN; timeout — срок. Сетки — тетраэдры без "
                      "пристенных призм, покрышка и стенки за колесом без трения: числа — проверка разгона, не "
                      "характеристика колеса."),
            "run_id": run_id, "created_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
            "strategies": strategies, "ranking": ranking, "best_per_wheel": best}


def cmd_merge(a) -> int:
    paths = sorted(glob.glob(os.path.join(a.inp, "**", "sweep_summary.json"), recursive=True))
    paths = [p for p in paths if "/_lost/" not in p.replace(os.sep, "/")]
    if not paths:
        print(f"нет sweep_summary.json в {a.inp}", file=sys.stderr)
        return 2
    try:
        res = merge(paths, a.run_id)
    except ValueError as e:
        print("ОТКАЗ:", e, file=sys.stderr)
        return 2
    with open(a.out, "w", encoding="utf-8") as f:
        json.dump(res, f, ensure_ascii=False, indent=1)
        f.write("\n")
    for sid in res["ranking"]:
        g = res["strategies"][sid]
        line = "; ".join(f"{lab}: {w['verdict']}@{w.get('verdict_step')} дисб. "
                         f"{(w.get('criteria') or {}).get('imbalance_mean')}" for lab, w in sorted(g["wheels"].items()))
        print(f"{sid}: {line}")
    print("лучшая по колесу:", res["best_per_wheel"])
    print(f"записано: {a.out}")
    if a.expect:
        labels = {lab for g in res["strategies"].values() for lab in g["wheels"]}
        miss = [x for x in a.expect.split(",") if x and x not in labels]
        if miss:
            print("нет сводок колёс:", miss, file=sys.stderr)
            return 3
    return 0


def build_parser():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    r = sub.add_parser("run", help="одна стратегия для одного колеса")
    r.add_argument("--strategy", required=True, help="plan_sweep_*.json")
    r.add_argument("--case", required=True)
    r.add_argument("--level", required=True, choices=["coarse", "medium", "fine"])
    r.add_argument("--geom-out", required=True)
    r.add_argument("--label", default=None)
    r.add_argument("--np", type=int, required=True, help="процессов MPI на колесо")
    r.add_argument("--mpi-args", default="--use-hwthread-cpus")
    r.add_argument("--shroud", default="slip", choices=["stationary", "rotating", "slip"])
    r.add_argument("--vaneless", default="slip", choices=["as-is", "stationary", "slip"])
    r.add_argument("--su2-src", default=None)
    r.add_argument("--time-limit", type=float, required=True, help="секунд на всю стратегию")
    r.add_argument("--job-max-min", type=int, default=None, help="max_min задачи раннера: срок не позже сторожа")
    r.add_argument("--reserve-s", type=float, default=180.0)
    r.add_argument("--uptime-file", default="/proc/uptime")
    r.add_argument("--mini-iter", type=int, default=None, help="миниатюра: итераций на шаг, гейт не держит")
    r.add_argument("--work", required=True)
    r.add_argument("--out", required=True)
    m = sub.add_parser("merge", help="sweep_summary.json -> один JSON")
    m.add_argument("--in", dest="inp", required=True)
    m.add_argument("--out", required=True)
    m.add_argument("--run-id", default=None)
    m.add_argument("--expect", default="")
    return ap


def main(argv=None) -> int:
    argv = sys.argv[1:] if argv is None else argv
    if argv and argv[0] == "rc":
        return cmd_rc(argv[1:])
    a = build_parser().parse_args(argv)
    if a.cmd == "run":
        return cmd_run(a)
    return cmd_merge(a)


if __name__ == "__main__":
    sys.exit(main())
