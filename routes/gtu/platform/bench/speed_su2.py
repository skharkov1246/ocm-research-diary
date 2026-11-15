#!/usr/bin/env python3
"""Замер скорости SU2 на машине: секунды на итерацию канала колеса при разных разбиениях ядер (узел bench).

Задача-замер для сметы уровней 2–3 воронки (SU2 ступеней): сколько машинных часов стоит одна рабочая точка
канала на сетках medium и fine. Сходимость здесь не нужна и не проверяется: меряется время итерации на
фиксированном числе итераций; время до сходимости даёт протокол PROTOCOL.md (заблокирован до сошедшегося
расчёта на CPU). Цена точки = (итераций на точку) × (с/итерацию) — это пересчёт, а не замер (speed_merge.py).

Что делает `run`:
  1. сетка канала каждого уровня — routes/gtu/platform/geom/impeller.py --no-step (время записывается: новую
     сетку строит каждая концепция уровня 2–3);
  2. для каждой конфигурации «NPxK» — K одновременных расчётов run_case.py по NP процессов MPI на одной сетке
     (так считаются точки характеристики: несколько точек на машине одновременно); план из двух стадий:
     первый порядок с рампой частоты и давления (--warm итераций), затем второй порядок MUSCL (--iters итераций)
     с рестарта — режим, в котором идут итерации рабочей точки;
  3. время итерации — первая колонка «Time(sec)» истории SU2 (время текущей итерации; вторая с тем же именем —
     среднее с начала), медиана второй половины стадии; накладные (запуск, чтение и разбиение сетки, запись) —
     время стадии минус сумма времён итераций; пиковая память процессов SU2_CFD — опрос /proc раз в 2 с;
  4. производные: C = с/итер · NP / узлы (мкс·ядро на узел и итерацию) и итераций в час на машину = K·3600/с.

Срок: при --job-max-min (срок сторожа задачи раннера) конец считается от времени работы машины (/proc/uptime;
в контейнере — то же ядро): (max_min − 8 − --reserve-min) мин; иначе --budget-s от старта. Конфигурация, на
которую не осталось времени, пропускается (status skipped_time), недосчитанная стадия всё равно даёт время
итерации по готовым строкам истории.

  python3 speed_su2.py run --case core130 --levels medium,fine --config medium=64x1,32x2,16x4 \\
          --config fine=64x1,32x2 --nproc "$NPROC" --job-max-min 75 --su2-src /opt/SU2-src --out "$GTU_OUT"
Выход: <out>/speed_su2_<case>.json (+ папки прогонов в --work, не выгружаются). Код 0 — каждый уровень
измерен хотя бы одной конфигурацией; 3 — нет; 2 — ошибка входа.
"""
from __future__ import annotations

import argparse
import datetime
import json
import os
import platform
import re
import shutil
import statistics
import subprocess
import sys
import threading
import time

HERE = os.path.dirname(os.path.abspath(__file__))
PLATFORM = os.path.dirname(HERE)
RUN_CASE = os.path.join(PLATFORM, "su2", "run_case.py")
IMPELLER = os.path.join(PLATFORM, "geom", "impeller.py")
CFG_RE = re.compile(r"^(\d+)x(\d+)$")


# ------------------------------------------------------------------ разбор
def parse_configs(items: list[str], levels: list[str], nproc: int) -> dict:
    """['medium=64x1,32x2', 'fine=64x1'] -> {'medium': [(64, 1), (32, 2)], 'fine': [(64, 1)]}."""
    out = {}
    for it in items:
        if "=" not in it:
            raise ValueError(f"--config {it!r}: ожидается уровень=NPxK[,NPxK]")
        lvl, spec = it.split("=", 1)
        lvl = lvl.strip()
        if lvl not in levels:
            raise ValueError(f"--config {it!r}: уровня {lvl} нет в --levels")
        cfgs = []
        for s in spec.split(","):
            m = CFG_RE.match(s.strip())
            if not m:
                raise ValueError(f"--config {it!r}: {s!r} не вида NPxK")
            np_, k = int(m.group(1)), int(m.group(2))
            if np_ < 1 or k < 1 or np_ * k > nproc:
                raise ValueError(f"--config {it!r}: {np_}x{k} = {np_ * k} процессов больше ядер ({nproc}) или ноль")
            cfgs.append((np_, k))
        out[lvl] = cfgs
    missing = [lv for lv in levels if lv not in out]
    if missing:
        raise ValueError(f"нет --config для уровней {missing}")
    return out


def read_history_times(path: str) -> list[float]:
    """Время каждой итерации: первая колонка «Time(sec)» истории SU2 (вторая с тем же именем — среднее)."""
    import csv
    if not os.path.isfile(path):
        return []  # стадия не началась (срок вышел раньше) — замера нет, это не ошибка
    with open(path, encoding="utf-8", errors="replace") as f:
        rows = [r for r in csv.reader(f) if r]
    if not rows:
        return []
    head = [h.strip().strip('"') for h in rows[0]]
    if "Time(sec)" not in head:
        return []
    j = head.index("Time(sec)")
    out = []
    for r in rows[1:]:
        if len(r) != len(head):
            continue  # недописанная строка при остановке по сроку
        try:
            v = float(r[j])
        except ValueError:
            continue
        if v == v and v > 0:
            out.append(v)
    return out


def iter_stats(times: list[float]) -> dict:
    """Медиана второй половины (установившийся режим), среднее, сумма."""
    if not times:
        return dict(rows=0)
    half = times[len(times) // 2:] if len(times) >= 4 else times
    return dict(rows=len(times), s_iter_median=statistics.median(half), s_iter_mean=statistics.fmean(times),
                s_iter_min=min(times), s_iter_max=max(times), sum_s=sum(times))


def summarize_run(run_dir: str, stages=("w1", "m2")) -> dict:
    """Итог одного run_case.py: статусы стадий, узлы сетки, время итерации и накладные по стадиям."""
    rec = dict(dir=os.path.basename(run_dir))
    try:
        s = json.load(open(os.path.join(run_dir, "run_summary.json"), encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as e:
        rec["error"] = f"run_summary.json: {type(e).__name__}"
        s = {}
    rec["nodes"] = (s.get("inputs", {}).get("mesh_header") or {}).get("NPOIN")
    rec["su2"] = (s.get("su2") or {}).get("version")
    by_name = {st.get("name"): st for st in s.get("stages", [])}
    for name in stages:
        st = by_name.get(name, {})
        t = read_history_times(os.path.join(run_dir, f"history_{name}.csv"))
        d = iter_stats(t)
        d.update(status=st.get("status"), seconds=st.get("seconds"))
        if st.get("seconds") is not None and d.get("sum_s") is not None:
            d["overhead_s"] = round(float(st["seconds"]) - d["sum_s"], 2)
        rec[name] = d
    return rec


def aggregate(level: str, np_: int, k: int, runs: list[dict], nodes, rss_peak_gb) -> dict:
    """Конфигурация: медиана по K расчётам; C = с/итер·NP/узлы (мкс), итераций в час на машину."""
    out = dict(level=level, np=np_, k=k, procs=np_ * k, nodes=nodes, runs=runs, rss_peak_gb=rss_peak_gb)
    for name, key in (("m2", "muscl"), ("w1", "first_order")):
        vals = [r[name]["s_iter_median"] for r in runs if r.get(name, {}).get("rows", 0) >= 3]
        if not vals:
            out[key] = dict(measured=False)
            continue
        s = statistics.median(vals)
        d = dict(measured=True, s_iter=s, s_iter_runs=vals, iters_per_hour_machine=k * 3600.0 / s)
        if nodes:
            d["C_us_core_node_iter"] = s * np_ / nodes * 1e6
        ov = [r[name].get("overhead_s") for r in runs if r.get(name, {}).get("overhead_s") is not None]
        if ov:
            d["overhead_s_median"] = statistics.median(ov)
        out[key] = d
    return out


# ------------------------------------------------------------------ машина
def machine_info(nproc: int) -> dict:
    info = dict(nproc=nproc, os_cpu_count=os.cpu_count(), itype=os.environ.get("GTU_ITYPE"),
                az=os.environ.get("GTU_AZ_MACHINE"), python=sys.version.split()[0], platform=platform.platform())
    try:
        txt = open("/proc/cpuinfo").read()
        m = re.search(r"^model name\s*:\s*(.+)$", txt, re.M)
        info["cpu_model"] = m.group(1).strip() if m else None
        info["threads_per_core"] = None
        sib, cores = re.search(r"^siblings\s*:\s*(\d+)", txt, re.M), re.search(r"^cpu cores\s*:\s*(\d+)", txt, re.M)
        if sib and cores:
            info["threads_per_core"] = int(sib.group(1)) // max(1, int(cores.group(1)))
    except OSError:
        pass
    try:
        m = re.search(r"^MemTotal:\s*(\d+)", open("/proc/meminfo").read(), re.M)
        info["mem_total_gb"] = round(int(m.group(1)) / 1024 ** 2, 2) if m else None
    except OSError:
        pass
    return info


def uptime_s() -> float:
    try:
        return float(open("/proc/uptime").read().split()[0])
    except (OSError, ValueError):
        return 0.0


class RssSampler(threading.Thread):
    """Сумма RSS процессов с данным именем (опрос /proc), пик за время жизни."""

    def __init__(self, comm: str = "SU2_CFD", period: float = 2.0):
        super().__init__(daemon=True)
        self.comm, self.period, self.peak_kb, self._halt = comm, period, 0, threading.Event()

    def run(self):
        while not self._halt.is_set():
            tot = 0
            for d in os.listdir("/proc"):
                if not d.isdigit():
                    continue
                try:
                    if open(f"/proc/{d}/comm").read().strip() != self.comm:
                        continue
                    for ln in open(f"/proc/{d}/status"):
                        if ln.startswith("VmRSS:"):
                            tot += int(ln.split()[1])
                            break
                except (OSError, ValueError):
                    continue
            self.peak_kb = max(self.peak_kb, tot)
            self._halt.wait(self.period)

    def stop(self) -> float:
        self._halt.set()
        self.join(5)
        return round(self.peak_kb / 1024 ** 2, 3)


# ------------------------------------------------------------------ запуск
def build_mesh(case: str, level: str, out_dir: str, threads: int, log: str, extra: list[str]) -> dict:
    t0 = time.monotonic()
    cmd = [sys.executable, IMPELLER, "--case", case, "--level", level, "--no-step", "--out", out_dir,
           "--threads", str(threads)] + extra
    with open(log, "w") as fh:
        rc = subprocess.call(cmd, stdout=fh, stderr=subprocess.STDOUT)
    dt = time.monotonic() - t0
    info = dict(rc=rc, mesh_s=round(dt, 1), dir=out_dir)
    try:
        j = json.load(open(os.path.join(out_dir, f"impeller_{case}_{level}.json"), encoding="utf-8"))
        q = (j.get("mesh") or {}).get("quality") or {}
        info.update(nodes=q.get("n_nodes"), cells=q.get("n_cells"), t_geometry_s=(j.get("mesh") or {}).get("t_geometry_s"))
    except (OSError, json.JSONDecodeError):
        pass
    return info


def run_config(a, level: str, np_: int, k: int, geom: str, plan_path: str, limit_s: int) -> tuple[list[dict], float]:
    procs = []
    sampler = RssSampler()
    sampler.start()
    for i in range(k):
        d = os.path.join(a.work, "runs", f"{a.case}_{level}_{np_}x{k}_{i}")
        shutil.rmtree(d, ignore_errors=True)
        os.makedirs(d, exist_ok=True)
        cmd = [sys.executable, RUN_CASE, "--case", a.case, "--level", level, "--geom-out", geom, "--plan", plan_path,
               "--np", str(np_), "--shroud", a.shroud, "--vaneless", a.vaneless, "--time-limit", str(limit_s),
               "--out", d]
        if a.su2_src:
            cmd += ["--su2-src", a.su2_src]
        if a.mpi_args is not None:
            cmd.append(f"--mpi-args={a.mpi_args}")  # с «=»: значение вида --use-hwthread-cpus argparse иначе примет за ключ
        fh = open(os.path.join(d, "run_case.out"), "w")
        procs.append((d, subprocess.Popen(cmd, stdout=fh, stderr=subprocess.STDOUT), fh))
    for d, p, fh in procs:
        p.wait()
        fh.close()
    rss = sampler.stop()
    return [dict(summarize_run(d), returncode=p.returncode) for d, p, _ in procs], rss


def cmd_run(a) -> int:
    t_start = time.monotonic()
    levels = [x.strip() for x in a.levels.split(",") if x.strip()]
    try:
        cfgs = parse_configs(a.config, levels, a.nproc)
    except ValueError as e:
        print(f"[speed_su2] {e}", file=sys.stderr)
        return 2
    os.makedirs(a.work, exist_ok=True)
    os.makedirs(a.out, exist_ok=True)
    if a.job_max_min:
        end_rel = (a.job_max_min - 8 - a.reserve_min) * 60 - uptime_s()
    else:
        end_rel = a.budget_s
    deadline = t_start + end_rel
    print(f"[speed_su2] срок замера {end_rel / 60:.1f} мин; уровни {levels}; конфигурации {cfgs}", flush=True)
    warm_ramp = max(1, int(a.warm * 0.75))
    plan = [dict(name="w1", frac=a.frac, pfrac=a.pfrac, ramp=warm_ramp, pramp=warm_ramp, iter=a.warm, muscl=False),
            dict(name="m2", frac=a.frac, pfrac=a.pfrac, ramp=0, iter=a.iters, muscl=True)]
    plan_path = os.path.join(a.work, "plan_speed.json")
    json.dump(plan, open(plan_path, "w"))
    doc = dict(
        method=("speed_su2.py: K одновременных run_case.py по NP процессов MPI на одной сетке канала; план — "
                f"{a.warm} итераций первого порядка с рампой и {a.iters} итераций MUSCL с рестарта; время итерации — "
                "медиана второй половины стадии по колонке Time(sec) (время итерации) истории SU2; C = с/итер·NP/узлы"),
        date=datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        machine=machine_info(a.nproc), commit=_git_head(),
        inputs=dict(case=a.case, levels=levels, configs={lv: [f"{n}x{k}" for n, k in c] for lv, c in cfgs.items()},
                    plan=plan, shroud=a.shroud, vaneless=a.vaneless, frac=a.frac, pfrac=a.pfrac,
                    job_max_min=a.job_max_min, reserve_min=a.reserve_min, config_limit_s=a.config_limit_s,
                    mesh_threads=a.mesh_threads or a.nproc),
        meshes={}, configs=[], skipped=[],
        notes=["время итерации не зависит от сходимости расчёта лишь приближённо: число линейных итераций FGMRES "
               "ограничено LINEAR_SOLVER_ITER шаблона и может меняться по ходу расчёта",
               "частота и давление на выходе — доли номинала из --frac/--pfrac: режим разгона, не рабочая точка"])
    out_json = os.path.join(a.out, f"speed_su2_{a.case}.json")

    def save():
        doc["elapsed_s"] = round(time.monotonic() - t_start, 1)
        with open(out_json, "w", encoding="utf-8") as f:
            json.dump(doc, f, ensure_ascii=False, indent=1)

    for level in levels:
        geom = os.path.join(a.geom_root, f"{a.case}_{level}") if a.geom_root else os.path.join(a.work, "geom", f"{a.case}_{level}")
        su2_mesh = os.path.join(geom, f"passage_{a.case}_{level}.su2")
        if os.path.isfile(su2_mesh):
            minfo = dict(rc=0, reused=True, dir=geom)
            try:
                q = json.load(open(os.path.join(geom, f"impeller_{a.case}_{level}.json")))["mesh"]["quality"]
                minfo.update(nodes=q.get("n_nodes"), cells=q.get("n_cells"))
            except (OSError, KeyError, TypeError, json.JSONDecodeError):
                pass
        elif time.monotonic() > deadline - 300:
            doc["skipped"].append(dict(level=level, what="mesh", reason="skipped_time"))
            save()
            continue
        else:
            print(f"[speed_su2] сетка {a.case}/{level} ...", flush=True)
            minfo = build_mesh(a.case, level, geom, a.mesh_threads or a.nproc,
                               os.path.join(a.out, f"geom_{a.case}_{level}.log"), a.impeller_args.split())
        doc["meshes"][level] = minfo
        print(f"[speed_su2] сетка {level}: {minfo}", flush=True)
        save()
        if minfo.get("rc") != 0 or not os.path.isfile(su2_mesh):
            doc["skipped"].append(dict(level=level, what="configs", reason="mesh_failed"))
            continue
        for np_, k in cfgs[level]:
            remain = deadline - time.monotonic()
            if remain < a.min_config_s:
                doc["skipped"].append(dict(level=level, np=np_, k=k, reason="skipped_time", remain_s=round(remain)))
                save()
                continue
            limit = int(min(a.config_limit_s, remain))
            print(f"[speed_su2] {level} {np_}x{k}: срок {limit} с", flush=True)
            try:
                runs, rss = run_config(a, level, np_, k, geom, plan_path, limit)
                nodes = next((r["nodes"] for r in runs if r.get("nodes")), minfo.get("nodes"))
                agg = aggregate(level, np_, k, runs, nodes, rss)
            except Exception as e:  # noqa: BLE001 — сбой одной конфигурации не теряет остальные замеры
                doc["skipped"].append(dict(level=level, np=np_, k=k, reason="error", error=f"{type(e).__name__}: {e}"))
                save()
                continue
            doc["configs"].append(agg)
            mu = agg["muscl"]
            print(f"[speed_su2] {level} {np_}x{k}: MUSCL {mu.get('s_iter')} с/итер, C {mu.get('C_us_core_node_iter')} мкс, "
                  f"{mu.get('iters_per_hour_machine')} итер/ч на машину, RSS {rss} ГБ", flush=True)
            save()
    measured = {c["level"] for c in doc["configs"] if c["muscl"].get("measured")}
    doc["levels_measured"] = sorted(measured)
    doc["ok"] = all(lv in measured for lv in levels)
    save()
    print(f"[speed_su2] итог: {out_json}; измерены уровни {sorted(measured)}", flush=True)
    return 0 if doc["ok"] else 3


def _git_head():
    try:
        return subprocess.run(["git", "-C", PLATFORM, "rev-parse", "HEAD"], capture_output=True, text=True,
                              timeout=20).stdout.strip() or None
    except (OSError, subprocess.SubprocessError):
        return None


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    r = sub.add_parser("run", help="замер на этой машине")
    r.add_argument("--case", default="core130")
    r.add_argument("--levels", default="medium,fine")
    r.add_argument("--config", action="append", default=[], help="уровень=NPxK[,NPxK...]; можно несколько раз")
    r.add_argument("--warm", type=int, default=40, help="итераций первого порядка (рампа)")
    r.add_argument("--iters", type=int, default=100, help="итераций MUSCL")
    r.add_argument("--frac", type=float, default=0.6, help="доля номинальной частоты")
    r.add_argument("--pfrac", type=float, default=0.95, help="доля оценки давления на выходе")
    r.add_argument("--shroud", default="slip")
    r.add_argument("--vaneless", default="slip")
    r.add_argument("--mpi-args", default=None, help="ключи mpirun (по умолчанию — как в run_case.py); писать с «=»: --mpi-args=--use-hwthread-cpus")
    r.add_argument("--nproc", type=int, default=int(os.environ.get("NPROC") or os.cpu_count() or 1))
    r.add_argument("--mesh-threads", type=int, default=0, help="потоков Gmsh (0 — все ядра)")
    r.add_argument("--impeller-args", default="", help="доп. ключи impeller.py")
    r.add_argument("--geom-root", default="", help="папка с готовыми сетками <case>_<level> (иначе строятся)")
    r.add_argument("--su2-src", default="", help="исходник SU2 для сверки ключей cfg")
    r.add_argument("--job-max-min", type=int, default=0, help="срок сторожа задачи раннера, мин (0 — --budget-s)")
    r.add_argument("--reserve-min", type=int, default=6, help="запас до конца срока, мин")
    r.add_argument("--budget-s", type=int, default=3600)
    r.add_argument("--config-limit-s", type=int, default=900, help="срок одной конфигурации, с")
    r.add_argument("--min-config-s", type=int, default=120, help="меньше осталось — конфигурация пропускается")
    r.add_argument("--work", default="/tmp/gtu_speed_su2")
    r.add_argument("--out", required=True)
    a = ap.parse_args(argv)
    return cmd_run(a)


if __name__ == "__main__":
    sys.exit(main())
