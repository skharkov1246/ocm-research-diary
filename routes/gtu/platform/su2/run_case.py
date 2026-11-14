#!/usr/bin/env python3
"""Расчёт одного межлопаточного канала центробежного колеса в SU2 (узел platform/su2 трека ГТУ).

Что делает:
  1. берёт сетку канала из ../geom/out/<case>_<level>/ (нет — строит её: python ../geom/impeller.py);
  2. читает из JSON сетки угол периодичности, угловую скорость, p01, T01 и одномерные PR, КПД, расход;
  3. заполняет шаблон impeller_passage.cfg.tpl — одну или несколько стадий (разгон частоты и давления
     на выходе, затем расчёт на рабочей точке; каждая следующая стадия стартует с решения предыдущей);
  4. запускает SU2_CFD (mpirun -np N) с общим лимитом времени;
  5. разбирает history: невязки, расход на входе и выходе, полная степень сжатия, изоэнтропический КПД
     по полным параметрам (свой расчёт и столбцы TURBO_PERF SU2), пишет run_summary.json в папку прогона.

Примеры:
  python run_case.py --case core130 --level coarse --plan smoke           # дымовой прогон, 4 ядра, 20 мин
  python run_case.py --case core130 --plan smoke --dry                    # только cfg + SU2_CFD -d
  python run_case.py --parse runs/core130_coarse_smoke                    # разобрать готовый прогон
  python run_case.py --case taurus60 --plan design --np 64 --time-limit 7200 --out "$GTU_OUT"   # AWS

SU2_CFD ищется так: --su2-bin, затем $SU2_RUN/SU2_CFD, затем PATH. Ключи cfg сверяются с
Common/src/CConfig.cpp той же версии, если найден исходник (--su2-src или $SU2_HOME).
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import os
import re
import shutil
import signal
import string
import subprocess
import sys
import time
from datetime import datetime, timezone

HERE = os.path.dirname(os.path.abspath(__file__))
GEOM = os.path.normpath(os.path.join(HERE, "..", "geom"))
TEMPLATE = os.path.join(HERE, "impeller_passage.cfg.tpl")
RESULTS = os.path.join(HERE, "su2_results.json")

CP_AIR = None  # считается из GAMMA и R

# Численная схема по умолчанию (подобрана дымовыми прогонами core130 coarse 14.11.2026, см. README).
NUMERICS = dict(
    TURB="SST", SA_OPTIONS="", INIT_MACH=0.4, AOA=90.0, TU=0.03, MUT_RATIO=10.0,
    GAMMA=1.4, R_GAS=287.058, RELAX_IN=0.3, RELAX_OUT=0.3, NSPAN=10,
    FLOW_SCHEME="ROE", MUSCL="NO", LIMITER="VENKATAKRISHNAN", VENKAT=0.05, ENTROPY_FIX=0.05,
    CFL=5.0, CFL_ADAPT="YES", CFL_ADAPT_PARAM="0.5, 1.1, 1.0, 20.0, 0.1", LIN_PREC="ILU", LIN_ITER=10,
    RES_MIN=-10, OUTPUT_FILES="RESTART, PARAVIEW", WRT_FREQ=500, EXTRA="",
)

# Планы стадий. frac — доля номинальной частоты в конце стадии, pfrac — доля оценки статического давления
# на выходе (estimate_p_out) в конце стадии; ramp — длина линейной рампы частоты и давления в итерациях
# (от значений конца предыдущей стадии); muscl — второй порядок по пространству.
PLANS = {
    # дымовой прогон: разгон до номинала с первым порядком, затем второй порядок на рабочей точке
    "smoke": [
        dict(name="s1_ramp60", frac=0.6, pfrac=1.0, ramp=150, iter=200, muscl=False),
        dict(name="s2_ramp100", frac=1.0, pfrac=1.0, ramp=150, iter=250, muscl=False),
        dict(name="s3_point", frac=1.0, pfrac=1.0, ramp=0, iter=300, muscl=True),
    ],
    # одна стадия: разгон и расчёт подряд
    "single": [dict(name="s1", frac=1.0, pfrac=1.0, ramp=300, iter=1000, muscl=False)],
    # расчётная точка для облака: долгий второй порядок
    "design": [
        dict(name="s1_ramp", frac=1.0, pfrac=1.0, ramp=400, iter=800, muscl=False),
        dict(name="s2_point", frac=1.0, pfrac=1.0, ramp=0, iter=6000, muscl=True),
    ],
}


# ----------------------------------------------------------------------------------------------
# сетка и исходные данные
# ----------------------------------------------------------------------------------------------
def mesh_paths(case: str, level: str, geom_out: str | None = None):
    d = geom_out or os.path.join(GEOM, "out", f"{case}_{level}")
    return (os.path.join(d, f"passage_{case}_{level}.su2"),
            os.path.join(d, f"impeller_{case}_{level}.json"))


def ensure_mesh(case: str, level: str, geom_out: str | None, extra: list[str], log=print):
    su2, js = mesh_paths(case, level, geom_out)
    if os.path.isfile(su2) and os.path.isfile(js):
        return su2, js, False
    cmd = [sys.executable, os.path.join(GEOM, "impeller.py"), "--case", case, "--level", level,
           "--no-step"] + extra
    if geom_out:
        cmd += ["--out", geom_out]
    log(f"сетки нет ({su2}) — строю: {' '.join(cmd)}")
    r = subprocess.run(cmd, cwd=GEOM)
    if r.returncode != 0 or not (os.path.isfile(su2) and os.path.isfile(js)):
        raise SystemExit(f"impeller.py не построил сетку (код {r.returncode})")
    return su2, js, True


def read_mesh_header(su2: str):
    """NDIME, NELEM, NPOIN, маркеры и число их граней — без чтения всей сетки в память."""
    info, markers, tag = {}, {}, None
    with open(su2, encoding="utf-8", errors="replace") as f:
        for ln in f:
            if ln.startswith(("NDIME", "NELEM", "NPOIN", "NMARK")):
                k, v = ln.split("=")
                info[k.strip()] = int(v.split()[0])
            elif ln.startswith("MARKER_TAG"):
                tag = ln.split("=")[1].strip()
            elif ln.startswith("MARKER_ELEMS") and tag:
                markers[tag] = int(ln.split("=")[1])
    info["markers"] = markers
    return info


def estimate_p_out(d1: dict, frac: float, r_out: float, b_out: float, gamma: float, R: float) -> float:
    """Оценка статического давления на выходе области (r_out) по одномерным данным сетки.

    Полное давление — по d1.PR (ступень) с эффективным КПД, пересчитанным на долю частоты frac
    (работа ~ frac^2, расход ~ frac); скорость: Ctheta = psi*U2*r2/r_out (сохранение момента в
    безлопаточном участке), Cm — из неразрывности при ширине b_out. Это граничное условие-оценка,
    а не результат: точку характеристики задаёт само давление на выходе.
    """
    cp = gamma * R / (gamma - 1)
    k = gamma / (gamma - 1)
    T01, p01 = d1["T01"], d1["p01"]
    w_nom = d1["w"]
    tau_nom = 1 + w_nom / (cp * T01)
    eta_eff = (d1["PR"] ** (1 / k) - 1) / (tau_nom - 1)
    w = w_nom * frac ** 2
    T02 = T01 + w / cp
    p02 = p01 * (1 + eta_eff * (T02 / T01 - 1)) ** k
    U2 = d1["U2"] * frac
    ct = d1["psi"] * U2 * d1["r2"] / r_out
    m = d1["m"] * frac
    A = 2 * math.pi * r_out * b_out
    cm, p = 0.0, p02
    for _ in range(50):
        T = T02 - (ct ** 2 + cm ** 2) / (2 * cp)
        p = p02 * (T / T02) ** k
        cm_new = m / (p / (R * T) * A)
        if abs(cm_new - cm) < 1e-6:
            break
        cm = cm_new
    return p


# ----------------------------------------------------------------------------------------------
# SU2
# ----------------------------------------------------------------------------------------------
def find_su2(bin_arg: str | None):
    cands = []
    if bin_arg:
        cands.append(bin_arg)
    if os.environ.get("SU2_RUN"):
        cands.append(os.path.join(os.environ["SU2_RUN"], "SU2_CFD"))
    w = shutil.which("SU2_CFD")
    if w:
        cands.append(w)
    for c in cands:
        if os.path.isfile(c) and os.access(c, os.X_OK):
            return os.path.abspath(c)
    return None


def su2_version(exe: str) -> str:
    try:
        out = subprocess.run([exe, "--help"], capture_output=True, text=True, timeout=30).stdout
        return out.splitlines()[0].strip() if out else "?"
    except Exception as e:  # noqa: BLE001
        return f"? ({e})"


def find_cconfig(src_arg: str | None):
    for base in (src_arg, os.environ.get("SU2_HOME")):
        if base:
            p = os.path.join(base, "Common", "src", "CConfig.cpp")
            if os.path.isfile(p):
                return p
    return None


def check_keys(cfg_text: str, cconfig: str):
    """Ключи cfg, не зарегистрированные в CConfig.cpp (add...Option("KEY")."""
    src = open(cconfig, encoding="utf-8", errors="replace").read()
    known = set(re.findall(r'add[A-Za-z]+Option\("([A-Z0-9_]+)"', src))
    keys = {m.group(1) for m in re.finditer(r"^([A-Z][A-Z0-9_]*)\s*=", cfg_text, re.M)}
    return sorted(keys - known), len(known)


def mpi_cmd(np_: int, exe: str, cfg: str):
    if np_ <= 1:
        return [exe, cfg]
    mpirun = shutil.which("mpirun") or shutil.which("mpiexec")
    if not mpirun:
        raise SystemExit("mpirun не найден (apt-get install -y openmpi-bin), или --np 1")
    cmd = [mpirun, "-np", str(np_), "--bind-to", "none"]
    if hasattr(os, "geteuid") and os.geteuid() == 0:
        cmd.append("--allow-run-as-root")
    return cmd + [exe, cfg]


def su2_env(np_: int):
    """Окружение запуска. Сборка SU2 с OpenMPI 4.1 на одном процессе (без mpirun или -np 1) падает в
    MPI_Win_create («MPI_ERR_WIN: invalid window»); с компонентом osc pt2pt идёт (проверено 14.11.2026)."""
    env = os.environ.copy()
    if np_ <= 1:
        env.setdefault("OMPI_MCA_osc", "pt2pt")
    return env


def run_limited(cmd, cwd, log_path, deadline, env=None):
    """Запуск с лимитом до deadline (time.monotonic). Возвращает (код, статус, секунды)."""
    t0 = time.monotonic()
    budget = deadline - t0
    if budget <= 5:
        return None, "no_time", 0.0
    with open(log_path, "w", encoding="utf-8") as lf:
        p = subprocess.Popen(cmd, cwd=cwd, stdout=lf, stderr=subprocess.STDOUT, start_new_session=True, env=env)
        try:
            rc = p.wait(timeout=budget)
            status = "finished" if rc == 0 else "error"
        except subprocess.TimeoutExpired:
            os.killpg(p.pid, signal.SIGTERM)
            try:
                p.wait(timeout=30)
            except subprocess.TimeoutExpired:
                os.killpg(p.pid, signal.SIGKILL)
                p.wait()
            rc, status = p.returncode, "timeout"
    txt = open(log_path, encoding="utf-8", errors="replace").read()
    if "NaN detected" in txt or "has diverged" in txt:
        status = "diverged"
    return rc, status, time.monotonic() - t0


# ----------------------------------------------------------------------------------------------
# разбор
# ----------------------------------------------------------------------------------------------
def read_history(path: str):
    with open(path, encoding="utf-8") as f:
        rows = [r for r in csv.reader(f) if r]
    if not rows:
        return [], {}
    head = [h.strip().strip('"') for h in rows[0]]
    first = {}
    for j, h in enumerate(head):
        first.setdefault(h, j)  # повторённый заголовок (например, Time(sec)) — берётся первый столбец
    data = {h: [] for h in first}
    for r in rows[1:]:
        if len(r) != len(head):
            continue  # недописанная строка при остановке по лимиту
        for h, j in first.items():
            try:
                data[h].append(float(r[j]))
            except ValueError:
                data[h].append(float("nan"))
    return list(first), data


def _finite(x):
    return x is not None and isinstance(x, float) and math.isfinite(x)


def analyse_history(path: str, gamma: float, n_passages: int, window: int = 50):
    head, d = read_history(path)
    n = len(d.get("Inner_Iter", []))
    if n == 0:
        return {"rows": 0}
    out = {"rows": n, "iter_first": int(d["Inner_Iter"][0]), "iter_last": int(d["Inner_Iter"][-1])}
    res = {}
    for h in head:
        if h.startswith("rms["):
            v = d[h]
            fin = [x for x in v if math.isfinite(x)]
            if not fin:
                res[h] = {"nan": True}
                continue
            imax = max(range(len(v)), key=lambda i: v[i] if math.isfinite(v[i]) else -1e30)
            res[h] = {"first": v[0], "max": v[imax], "iter_max": int(d["Inner_Iter"][imax]),
                      "last": v[-1], "drop_from_max": v[imax] - v[-1], "drop_from_first": v[0] - v[-1]}
            # наклон на второй половине (порядков на 100 итераций)
            half = v[n // 2:]
            if len(half) >= 10 and all(math.isfinite(x) for x in half):
                xs = list(range(len(half)))
                mx, my = sum(xs) / len(xs), sum(half) / len(half)
                sxx = sum((x - mx) ** 2 for x in xs)
                res[h]["slope_per_100_last_half"] = 100 * sum((x - mx) * (y - my) for x, y in zip(xs, half)) / sxx
    out["residuals"] = res

    def last(col):
        v = d.get(col)
        return v[-1] if v else None

    def win_change(col):
        v = d.get(col)
        if not v or len(v) < window + 1:
            return None
        seg = v[-window - 1:]
        if not all(math.isfinite(x) for x in seg) or seg[-1] == 0:
            return None
        return (max(seg) - min(seg)) / abs(seg[-1])

    m_in, m_out = last("Avg_Massflow(inlet)"), last("Avg_Massflow(outlet)")
    pt1, pt2 = last("Avg_TotalPress(inlet)"), last("Avg_TotalPress(outlet)")
    tt1, tt2 = last("Avg_TotalTemp(inlet)"), last("Avg_TotalTemp(outlet)")
    perf = {"m_in_passage": m_in, "m_out_passage": m_out}
    if _finite(m_in) and _finite(m_out) and m_in != 0:
        perf["mass_imbalance"] = abs(m_in + m_out) / abs(m_in)
        perf["mass_sign_ok"] = bool(m_in > 0 > m_out)
        perf["m_wheel"] = m_in * n_passages
    if all(_finite(x) for x in (pt1, pt2, tt1, tt2)) and pt1 > 0 and tt1 > 0:
        pr, tr = pt2 / pt1, tt2 / tt1
        perf["PR_tt"] = pr
        perf["TR_tt"] = tr
        perf["eta_tt_isentropic"] = ((pr ** ((gamma - 1) / gamma) - 1) / (tr - 1)) if tr > 1 else None
        if _finite(m_in):
            perf["m_wheel_corrected"] = m_in * n_passages * math.sqrt(tt1 / 288.15) / (pt1 / 101325.0)
    for col, key in (("PRTT", "PRTT_su2"), ("TotTotEff", "TotTotEff_su2"), ("PRTS", "PRTS_su2"),
                     ("TotStaticEff", "TotStaticEff_su2"), ("MassFlowIn_1", "MassFlowIn_turbo"),
                     ("MassFlowOut_1", "MassFlowOut_turbo"), ("EulerianWork", "EulerianWork_su2")):
        if col in d:
            perf[key] = last(col)
    perf["window_iters"] = window
    perf["window_rel_change"] = {k: win_change(c) for k, c in (
        ("m_in", "Avg_Massflow(inlet)"), ("pt_out", "Avg_TotalPress(outlet)"), ("tt_out", "Avg_TotalTemp(outlet)"))}
    out["performance"] = perf
    if "Time(sec)" in d or "Cur_Time" in d:
        col = "Cur_Time" if "Cur_Time" in d else "Time(sec)"
        out["sec_per_iter_last"] = d[col][-1]
    return out


def parse_log(path: str):
    if not os.path.isfile(path):
        return {}
    txt = open(path, encoding="utf-8", errors="replace").read()
    o = {}
    m = re.findall(r"Matched (\d+) points with a max distance of: (\S+?)\.?\s*$", txt, re.M)
    if m:
        o["periodic_matched_points"] = int(m[-1][0])
        o["periodic_max_distance"] = float(m[-1][1])
    o["bad_match"] = len(re.findall(r"Bad match", txt))
    o["nan"] = ("NaN detected" in txt) or ("has diverged" in txt)
    o["exit_success"] = "Exit Success" in txt
    for key, pat in (("inlet_area_cm2", r"Inlet area for Row 1: (\S+) cm"),
                     ("outlet_area_cm2", r"Oulet area for Row 1: (\S+) cm"),
                     ("recomputed_blades", r"Recomputed number of blades for Row 1: (\S+?)\.?$")):
        mm = re.search(pat, txt, re.M)
        if mm:
            try:
                o[key] = float(mm.group(1).rstrip("."))
            except ValueError:
                pass
    errs = re.findall(r"^Error in .*$\n?^(.*)$", txt, re.M)
    if errs:
        o["error"] = errs[-1].strip()
    mrank = re.search(r"Parallel Computation with (\d+) MPI", txt) or re.search(r"(\d+) MPI ranks", txt)
    if mrank:
        o["mpi_ranks"] = int(mrank.group(1))
    return o


def smoke_checks(stages_out: list, mass_tol: float):
    """Критерии дымового прогона из задания: нет расхождения, невязки падают, расход вход/выход ±mass_tol."""
    last = next((s for s in reversed(stages_out) if s.get("history", {}).get("rows")), None)
    if not last:
        return {"passed": False, "reason": "нет истории"}
    h = last["history"]
    rho = h["residuals"].get("rms[Rho]", {})
    perf = h.get("performance", {})
    any_nan = any(s.get("status") == "diverged" for s in stages_out)
    first_rho = next((s["history"]["residuals"]["rms[Rho]"] for s in stages_out
                      if s.get("history", {}).get("residuals", {}).get("rms[Rho]")), {})
    c = {
        "no_divergence": not any_nan,
        "rms_rho_first_stage_start": first_rho.get("first"),
        "rms_rho_last": rho.get("last"),
        "rms_rho_drop_over_run": (first_rho.get("max", float("nan")) - rho.get("last", float("nan")))
        if first_rho and rho else None,
        "rms_rho_slope_last_half_per_100": rho.get("slope_per_100_last_half"),
        "mass_imbalance": perf.get("mass_imbalance"),
        "mass_sign_ok": perf.get("mass_sign_ok"),
    }
    falling = (c["rms_rho_drop_over_run"] or 0) > 0 and (c["rms_rho_slope_last_half_per_100"] or 1) < 0
    c["residuals_falling"] = bool(falling)
    c["mass_within_tol"] = bool(_finite(c["mass_imbalance"]) and c["mass_imbalance"] <= mass_tol
                                and c["mass_sign_ok"])
    c["mass_tol"] = mass_tol
    c["passed"] = bool(c["no_divergence"] and c["residuals_falling"] and c["mass_within_tol"])
    return c


# ----------------------------------------------------------------------------------------------
# основной ход
# ----------------------------------------------------------------------------------------------
def build_cfg(tpl: str, params: dict) -> str:
    txt = string.Template(tpl).substitute(params)
    # пустая строка ${SA_OPTIONS} не мешает SU2, но уберём её для чистоты
    return "\n".join(ln for ln in txt.splitlines() if ln.strip() != "") + "\n"


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--case", default="core130", help="core130 | taurus60 | custom (как в impeller.py)")
    ap.add_argument("--level", default="coarse", choices=["coarse", "medium", "fine"])
    ap.add_argument("--geom-out", default=None, help="папка сетки (по умолчанию ../geom/out/<case>_<level>)")
    ap.add_argument("--impeller-args", default="", help="доп. ключи impeller.py, если сетку надо строить")
    ap.add_argument("--plan", default="smoke", help=f"план стадий: {', '.join(PLANS)} или путь к JSON")
    ap.add_argument("--iter-scale", type=float, default=1.0, help="множитель итераций всех стадий")
    ap.add_argument("--frac", type=float, default=None, help="доля номинальной частоты (переопределяет план)")
    ap.add_argument("--p-out", type=float, default=None, help="статическое давление на выходе, Па "
                    "(по умолчанию — оценка estimate_p_out x pfrac стадии)")
    ap.add_argument("--turb", default=None, choices=["SST", "SA"])
    ap.add_argument("--cfl", type=float, default=None)
    ap.add_argument("--set", action="append", default=[], metavar="KEY=VAL",
                    help="переопределить параметр шаблона (имя поля $ИМЯ), можно несколько раз")
    ap.add_argument("--extra", action="append", default=[], metavar="'KEY= VAL'",
                    help="дописать строку в конец cfg (ключ SU2 вне шаблона), можно несколько раз")
    ap.add_argument("--np", type=int, default=4, help="число процессов MPI (1 — без mpirun)")
    ap.add_argument("--time-limit", type=float, default=1200.0, help="общий лимит, с (по умолчанию 20 мин)")
    ap.add_argument("--su2-bin", default=None)
    ap.add_argument("--su2-src", default=None, help="исходник SU2 для сверки ключей (Common/src/CConfig.cpp)")
    ap.add_argument("--out", default=None, help="папка прогона (по умолчанию runs/<case>_<level>_<plan>)")
    ap.add_argument("--dry", action="store_true", help="только cfg и SU2_CFD -d (сухой прогон)")
    ap.add_argument("--parse", default=None, metavar="RUN_DIR", help="только разобрать готовый прогон")
    ap.add_argument("--mass-tol", type=float, default=0.05, help="допуск дисбаланса расхода (0,05 = 5 %%)")
    ap.add_argument("--record", action="store_true", help="дописать краткий итог в su2_results.json")
    a = ap.parse_args(argv)

    if a.parse:
        summ_path = os.path.join(a.parse, "run_summary.json")
        summ = json.load(open(summ_path, encoding="utf-8"))
        g = summ["inputs"]["gamma"]
        nps = summ["inputs"]["n_passages"]
        for st in summ["stages"]:
            hp = os.path.join(a.parse, st["history_file"])
            if os.path.isfile(hp):
                st["history"] = analyse_history(hp, g, nps)
            st["log"] = parse_log(os.path.join(a.parse, st["log_file"]))
        summ["smoke_checks"] = smoke_checks(summ["stages"], a.mass_tol)
        json.dump(summ, open(summ_path, "w", encoding="utf-8"), ensure_ascii=False, indent=1)
        print_summary(summ)
        return 0

    plan = PLANS.get(a.plan)
    if plan is None:
        plan = json.load(open(a.plan, encoding="utf-8"))
    plan = [dict(s) for s in plan]
    for s in plan:
        s["iter"] = max(1, int(round(s["iter"] * a.iter_scale)))
        s["ramp"] = int(round(s["ramp"] * a.iter_scale))
        if a.frac is not None:
            s["frac"] = a.frac

    su2_mesh, js, built = ensure_mesh(a.case, a.level, a.geom_out, a.impeller_args.split())
    mj = json.load(open(js, encoding="utf-8"))
    d1, s2 = mj["d1"], mj["su2"]
    per = s2["periodic"]
    omega = s2["rotation_rate_rad_s"]
    if abs(omega[0]) > 0 or abs(omega[1]) > 0 or list(per["axis"]) != [0, 0, 1]:
        raise SystemExit("шаблон рассчитан на ось вращения z: рампа частоты SU2 меняет только z-компоненту")
    n_pass = int(round(360.0 / per["angle_deg"]))
    hdr = read_mesh_header(su2_mesh)
    missing = [m for m in ("inlet", "outlet", "hub", "shroud", "blade", per["master"], per["slave"])
               if m not in hdr["markers"]]
    if missing:
        raise SystemExit(f"в сетке нет маркеров {missing}: {sorted(hdr['markers'])}")

    params = dict(NUMERICS)
    if a.turb:
        params["TURB"] = a.turb
    if a.cfl:
        params["CFL"] = a.cfl
    for kv in a.set:
        k, v = kv.split("=", 1)
        params[k.strip()] = v.strip()
    if a.extra:
        params["EXTRA"] = "\n".join(a.extra)
    if params["TURB"] == "SA" and not params["SA_OPTIONS"]:
        params["SA_OPTIONS"] = "SA_OPTIONS= NONE"
    gamma, R = float(params["GAMMA"]), float(params["R_GAS"])
    # статические p, T начального поля — из полных параметров входа при INIT_MACH
    m0 = float(params["INIT_MACH"])
    fs_t = d1["T01"] / (1 + 0.5 * (gamma - 1) * m0 ** 2)
    fs_p = d1["p01"] * (fs_t / d1["T01"]) ** (gamma / (gamma - 1))
    params.update(FS_P=f"{fs_p:.6g}", FS_T=f"{fs_t:.6g}")
    r_out = 1.3 * d1["r2"]  # длина безлопаточного участка в impeller.py: до 1,3 r2
    p_out_est = estimate_p_out(d1, 1.0, r_out, d1["b2"], gamma, R)

    su2 = find_su2(a.su2_bin)
    out = os.path.abspath(a.out or os.path.join(HERE, "runs", f"{a.case}_{a.level}_{a.plan if a.plan in PLANS else 'custom'}"))
    os.makedirs(os.path.join(out, "TURBOMACHINERY"), exist_ok=True)
    tpl = open(TEMPLATE, encoding="utf-8").read()
    cconfig = find_cconfig(a.su2_src)

    summ = {
        "created_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "su2": {"exe": su2, "version": su2_version(su2) if su2 else None, "np": a.np},
        "inputs": {"case": a.case, "level": a.level, "mesh": su2_mesh, "mesh_json": js, "mesh_built_now": built,
                   "mesh_header": hdr, "n_passages": n_pass, "periodic_angle_deg": per["angle_deg"],
                   "omega_nominal_rad_s": omega[2], "rpm_nominal": d1["rpm"], "p01": d1["p01"], "T01": d1["T01"],
                   "gamma": gamma, "R": R, "r_out": r_out, "b_out": d1["b2"],
                   "p_out_estimate_nominal": p_out_est, "plan": a.plan, "time_limit_s": a.time_limit},
        "d1_reference": {"m_wheel": d1["m"], "PR_stage": d1["PR"], "eta_comp_stage": d1.get("eta_comp_stage"),
                         "U2": d1["U2"], "note": "одномерная оценка gtu_explore.comp_stage — для сравнения, не эталон"},
        "numerics": params, "stages": [],
    }

    deadline = time.monotonic() + a.time_limit
    prev_frac, prev_p = 0.0, d1["p01"]
    prev_restart = None
    key_problems = None
    for i, st in enumerate(plan):
        frac, pfrac = st["frac"], st["pfrac"]
        p_target = a.p_out if a.p_out is not None else pfrac * estimate_p_out(d1, frac, r_out, d1["b2"], gamma, R)
        ramp = st["ramp"]
        name = st["name"]
        p = dict(params)
        p.update(
            RESTART="YES" if prev_restart else "NO",
            SOLUTION_FILENAME=prev_restart or "solution.dat",
            P01=f"{d1['p01']:.6g}", T01=f"{d1['T01']:.6g}",
            OMEGA_X="0.0", OMEGA_Y="0.0", OMEGA_Z=f"{omega[2] * frac:.10g}",
            PER_MASTER=per["master"], PER_SLAVE=per["slave"], ANGLE=f"{per['angle_deg']:.10g}",
            P_OUT=f"{p_target:.6g}", ITER=st["iter"], MUSCL="YES" if st.get("muscl") else "NO",
            HISTORY_NAME=f"history_{name}", RESTART_NAME=f"restart_{name}.dat", VOLUME_NAME=f"vol_{name}",
            MESH=su2_mesh,
        )
        rot_changes = ramp > 0 and abs(frac - prev_frac) > 1e-9
        p.update(RAMP_ROT="YES" if rot_changes else "NO",
                 RAMP_ROT_INI=f"{omega[2] * prev_frac:.10g}", RAMP_ROT_FREQ=1, RAMP_ROT_ITER=max(ramp, 1))
        p_changes = ramp > 0 and abs(p_target - prev_p) > 1.0
        p.update(RAMP_OUT="YES" if p_changes else "NO",
                 RAMP_OUT_INI=f"{prev_p:.6g}", RAMP_OUT_FREQ=1, RAMP_OUT_ITER=max(ramp, 1))
        for k in ("FLOW_SCHEME", "CFL", "LIMITER", "TURB"):
            if k.lower() in st:
                p[k] = st[k.lower()]
        cfg_txt = build_cfg(tpl, p)
        cfg_name = f"case_{name}.cfg"
        with open(os.path.join(out, cfg_name), "w", encoding="utf-8") as f:
            f.write(cfg_txt)
        if cconfig and key_problems is None:
            unknown, nknown = check_keys(cfg_txt, cconfig)
            key_problems = {"cconfig": cconfig, "registered_keys": nknown, "unknown": unknown}
            if unknown:
                raise SystemExit(f"ключей нет в {cconfig}: {unknown}")
        rec = {"name": name, "cfg": cfg_name, "frac": frac, "rpm": d1["rpm"] * frac, "p_out": p_target,
               "p_out_ramp_from": prev_p if p_changes else None, "rot_ramp_from_frac": prev_frac if rot_changes else None,
               "ramp_iter": ramp, "iter": st["iter"], "muscl": bool(st.get("muscl")),
               "restart_from": prev_restart, "history_file": f"history_{name}.csv", "log_file": f"su2_{name}.log"}
        summ["stages"].append(rec)
        if su2 is None:
            rec["status"] = "no_su2"
            continue
        if a.dry and prev_restart:
            # SU2_CFD -d читает файл рестарта, а его в сухом прогоне нет: cfg только записан
            rec.update(status="cfg_only")
            prev_frac, prev_p = frac, p_target
            continue
        if a.dry:
            rc, status, dt = run_limited([su2, "-d", cfg_name], out, os.path.join(out, f"dry_{name}.log"), deadline,
                                         env=su2_env(1))
            txt = open(os.path.join(out, f"dry_{name}.log"), encoding="utf-8", errors="replace").read()
            rec.update(status="dry_ok" if "Exit Success" in txt else "dry_failed", seconds=round(dt, 1))
            prev_frac, prev_p, prev_restart = frac, p_target, f"restart_{name}.dat"  # цепочка стадий для просмотра cfg
            continue
        rc, status, dt = run_limited(mpi_cmd(a.np, su2, cfg_name), out, os.path.join(out, rec["log_file"]), deadline,
                                     env=su2_env(a.np))
        rec.update(returncode=rc, status=status, seconds=round(dt, 1))
        hp = os.path.join(out, rec["history_file"])
        if os.path.isfile(hp):
            rec["history"] = analyse_history(hp, gamma, n_pass)
        rec["log"] = parse_log(os.path.join(out, rec["log_file"]))
        print(f"[{name}] {status} за {dt:.0f} с; строк истории {rec.get('history', {}).get('rows', 0)}", flush=True)
        if status != "finished":
            break
        prev_frac, prev_p, prev_restart = frac, p_target, rec_restart(out, name)

    summ["key_check"] = key_problems or {"note": "исходник SU2 не найден — ключи не сверялись (--su2-src)"}
    if su2 is None:
        summ["blocker"] = "SU2_CFD не найден: --su2-bin, $SU2_RUN или PATH (установка — install_su2.sh)"
    if not a.dry and su2:
        summ["smoke_checks"] = smoke_checks(summ["stages"], a.mass_tol)
    with open(os.path.join(out, "run_summary.json"), "w", encoding="utf-8") as f:
        json.dump(summ, f, ensure_ascii=False, indent=1)
    print_summary(summ)
    print(f"папка прогона: {out}")
    if a.record and not a.dry and su2:
        record(summ, out)
    return 0 if (a.dry or su2) else 2


def rec_restart(out: str, name: str):
    f = f"restart_{name}.dat"
    return f if os.path.isfile(os.path.join(out, f)) else None


def print_summary(s: dict):
    print(f"SU2: {s['su2'].get('version')}  np={s['su2'].get('np')}")
    i = s["inputs"]
    print(f"случай {i['case']}/{i['level']}: узлов {i['mesh_header'].get('NPOIN')}, ячеек {i['mesh_header'].get('NELEM')}, "
          f"каналов {i['n_passages']}, {i['rpm_nominal']:.0f} об/мин ном., оценка p_вых {i['p_out_estimate_nominal']:.0f} Па")
    for st in s["stages"]:
        h = st.get("history", {})
        lg = st.get("log", {})
        line = f"  {st['name']}: {st.get('status')} {st.get('seconds', '')} с, частота {st['frac']:.2f} ном., p_вых {st['p_out']:.0f} Па"
        if h.get("rows"):
            r = h["residuals"].get("rms[Rho]", {})
            pf = h["performance"]
            line += (f"; итераций {h['rows']}, rms[Rho] {r.get('first', float('nan')):.2f} -> {r.get('last', float('nan')):.2f}"
                     f" (макс {r.get('max', float('nan')):.2f})")
            if pf.get("mass_imbalance") is not None:
                line += f"; m_вх {pf['m_in_passage']:+.4f} m_вых {pf['m_out_passage']:+.4f} кг/с, дисбаланс {pf['mass_imbalance']:.2%}"
            if pf.get("PR_tt") is not None:
                eta = pf.get("eta_tt_isentropic")
                line += f"; PR_tt {pf['PR_tt']:.3f}, eta_tt {eta:.3f}" if eta is not None else f"; PR_tt {pf['PR_tt']:.3f}"
        if lg.get("periodic_matched_points"):
            line += f"; периодичность: {lg['periodic_matched_points']} точек, {lg['periodic_max_distance']:.1e}, плохих {lg['bad_match']}"
        if lg.get("error"):
            line += f"; ошибка: {lg['error']}"
        print(line)
    if s.get("smoke_checks"):
        c = s["smoke_checks"]
        print("дымовые проверки:", {k: c[k] for k in ("no_divergence", "residuals_falling", "mass_within_tol", "passed") if k in c})


def record(summ: dict, out: str):
    """Краткий итог прогона в su2_results.json (без больших файлов)."""
    try:
        res = json.load(open(RESULTS, encoding="utf-8"))
    except (FileNotFoundError, json.JSONDecodeError):
        res = {"about": "Итоги прогонов SU2 узла platform/su2 (run_case.py --record). Числа дымовых прогонов на "
                        "грубой сетке без пристенных призм — проверка постановки, не характеристика колеса.",
               "runs": []}
    last = next((st for st in reversed(summ["stages"]) if st.get("history", {}).get("rows")), {})
    res["runs"].append({
        "created_utc": summ["created_utc"], "su2": summ["su2"]["version"], "np": summ["su2"]["np"],
        "case": summ["inputs"]["case"], "level": summ["inputs"]["level"], "plan": summ["inputs"]["plan"],
        "mesh_nodes": summ["inputs"]["mesh_header"].get("NPOIN"), "n_passages": summ["inputs"]["n_passages"],
        "p_out_estimate_nominal": summ["inputs"]["p_out_estimate_nominal"],
        "stages": [{k: st.get(k) for k in ("name", "status", "seconds", "frac", "p_out", "ramp_iter", "iter", "muscl")}
                   | {"rows": st.get("history", {}).get("rows"),
                      "rms_rho": st.get("history", {}).get("residuals", {}).get("rms[Rho]"),
                      "performance": st.get("history", {}).get("performance"),
                      "periodic": {k: st.get("log", {}).get(k) for k in ("periodic_matched_points", "periodic_max_distance", "bad_match")}}
                   for st in summ["stages"]],
        "smoke_checks": summ.get("smoke_checks"), "d1_reference": summ["d1_reference"],
        "run_dir": os.path.relpath(out, HERE),
    })
    with open(RESULTS, "w", encoding="utf-8") as f:
        json.dump(res, f, ensure_ascii=False, indent=1)
    print(f"записано в {RESULTS}")


if __name__ == "__main__":
    sys.exit(main())
