#!/usr/bin/env python3
"""Пакеты задач routes/gtu/platform/aws/runner.py: линии частоты колёс в SU2 — KV-3 ступень 1 (Taurus 60 GS,
двусторонний вход, frame3) и core130 — на сетке medium, затем fine.

Пишет рядом:
  jobs_aws_su2_kv3_smoke.json — этап 0: одна машина, обе линии на coarse, единицы десятков итераций (~20 мин):
                                проверка всей цепочки на AWS (сборка SU2, MPI 16 и 64 процесса, срок, выгрузка);
  jobs_aws_su2_kv3.json       — этап A: medium (4 машины c7a.16xlarge, по 4 точки на машину x 16 процессов MPI);
  jobs_aws_su2_kv3_fine.json  — этап B: fine (6 машин), запускать после разбора этапа A, точки — по его итогам.
В каждом — блоки estimate (машино-часы x цена спота + диск, источник времени итерации), merge (что забрать и куда
положить сводку в ветке) и plan_check (вывод runner.py plan, если не --no-plan).

Машина: образ al2023, контейнер ubuntu:24.04; container_setup ставит apt-пакеты (компилятор, OpenMPI 4.1.6,
библиотеки gmsh), venv /opt/venv с пакетами routes/gtu/platform/docker/requirements.txt и собирает SU2 v8.6.0 с
MPI скриптом install_su2.sh (те же шаги, что стадия su2-build Dockerfile узла docker; коммит тега сверяется).
Команда: сетка на машине (frame3_geom.py или impeller.py), затем char_points.py run.

  python make_jobs_su2_kv3.py                    # оба файла + plan_check
  python make_jobs_su2_kv3.py --mini DIR         # миниатюры команд (coarse, единицы итераций) для проверки
"""
from __future__ import annotations

import argparse
import json
import math
import os
import shutil
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", "..", "..", ".."))
RUNNER = os.path.join(REPO, "routes", "gtu", "platform", "aws", "runner.py")
OUT_MEDIUM = os.path.join(HERE, "jobs_aws_su2_kv3.json")
OUT_FINE = os.path.join(HERE, "jobs_aws_su2_kv3_fine.json")
OUT_SMOKE = os.path.join(HERE, "jobs_aws_su2_kv3_smoke.json")
SMOKE_MAX_MIN = 60
BRANCH = "claude/lucid-allen-yvps5d"

SU2_COMMIT = "7048fd35f56d4bf52e586429d6746a36a10255f8"   # тег v8.6.0 (Dockerfile узла docker, сверено при сборке)
CONTAINER = "ubuntu:24.04"
CONTAINER_SETUP = [
    "apt-get update -q",
    "apt-get install -y -q --no-install-recommends build-essential python3 python3-venv git pkg-config "
    "ca-certificates openmpi-bin libopenmpi-dev libgl1 libglu1-mesa libxrender1 libxcursor1 libxfixes3 libxft2 "
    "libfontconfig1 libxinerama1 libgomp1",
    "python3 -m venv /opt/venv",
    "/opt/venv/bin/pip install -q --no-cache-dir -r routes/gtu/platform/docker/requirements.txt",
    "/opt/venv/bin/python -c 'import cantera, numpy, scipy, gmsh; print(\"py ok\", cantera.__version__, gmsh.__version__)'",
    "JOBS=\"$(nproc)\" bash routes/gtu/platform/su2/install_su2.sh /opt/su2 /opt/SU2-src > \"$GTU_OUT/_su2_build.log\" 2>&1 "
    "|| { tail -n 40 \"$GTU_OUT/_su2_build.log\"; exit 1; }",
    f"test \"$(git -C /opt/SU2-src rev-parse HEAD)\" = {SU2_COMMIT}",
    "/opt/su2/bin/SU2_CFD --help | head -n 1 | tee /opt/su2/VERSION",
]
# Окружение команды: в ubuntu:24.04 его никто не задаёт (в образе узла docker — ENV). OMPI_MCA_osc=^ucx — без него
# SU2_CFD под mpirun в docker без особых прав падает в MPI_Win_create (проверка образа 14.11.2026, Dockerfile);
# vader single copy = none — CMA запрещён без CAP_SYS_PTRACE.
ENV_LINES = [
    "set -e",
    "export PATH=\"/opt/venv/bin:/opt/su2/bin:$PATH\" SU2_RUN=/opt/su2/bin PYTHONDONTWRITEBYTECODE=1",
    "export OMPI_MCA_btl_vader_single_copy_mechanism=none OMPI_MCA_osc=^ucx",
    "W=/tmp/gtu_work; mkdir -p \"$W\"",
    "{ cat /opt/su2/VERSION; python -V; mpirun --version | head -n 1; echo \"nproc $(nproc) NPROC $NPROC\"; free -g; df -h /tmp; } "
    "> \"$GTU_OUT/env.txt\" 2>&1 || true",
]

# Линии: как строится сетка и где она лежит.
LINES = {
    "kv3_stage1": dict(case="custom",
                       mesh=("python routes/gtu/frame3/frame3_geom.py --out \"$W/geom\" --levels {level} --stages stage1 "
                             "--gzip-levels --threads \"$NPROC\" --results \"$GTU_OUT/frame3_geom_results.json\" "
                             "> \"$GTU_OUT/geom.log\" 2>&1"),
                       geom_out="$W/geom/stage1_{level}",
                       about="KV-3 ступень 1 (Taurus 60 GS, генератор): двусторонний вход, 25,55 кг/с (12,775 на сторону), "
                             "π 4,49, 14 950 об/мин, D2 0,695 м — frame3_point_results.json → stage1; сетка — одна сторона"),
    "core130": dict(case="core130",
                    mesh=("python routes/gtu/platform/geom/impeller.py --case core130 --level {level} --no-step "
                          "--out \"$W/geom/core130_{level}\" --threads \"$NPROC\" > \"$GTU_OUT/geom.log\" 2>&1"),
                    geom_out="$W/geom/core130_{level}",
                    about="core130: схема 130 gtu_explore (ядро кластера 4 x 9,95 кг/с, π 4,25, 22 197 об/мин, D2 0,435 м)"),
}

# Узлы сеток: coarse и medium — замер, fine — оценка (рост medium/coarse).
NODES = {
    ("kv3_stage1", "coarse"): 258565, ("kv3_stage1", "medium"): 1137992,
    ("core130", "coarse"): 139121, ("core130", "medium"): 584357,
}
NODES_SRC = {
    "kv3_stage1": "frame3_geom.py --levels coarse medium, 14.11.2026 (frame3_geom_results.json → stages.stage1.levels; "
                  "1 200 627 и 5 819 064 ячеек)",
    "core130": "routes/gtu/platform/bench/mesh_stats.json (core130 coarse 139 121, medium 584 357 узлов)",
}
for ln in ("kv3_stage1", "core130"):
    NODES[(ln, "fine")] = int(round(NODES[(ln, "medium")] ** 2 / NODES[(ln, "coarse")], -3))

# Время итерации: t = C·1e-6 · узлы / процессы / EFF, C — ядро-мкс на узел и итерацию. Замеры — здесь, на общей
# 4-ядерной машине (Xeon 2,1 ГГц, загрузка 9–16 от чужих задач), поэтому грубые (TIMING). Ожидаемая смета — по
# C_EXPECTED, срок сторожа max_min — по C_SIZING (с запасом, чтобы точка не обрывалась на сроке).
TIMING = {
    "ref_cost_estimate": {"s_per_iter": 1.2255, "nodes": 139121, "procs": 4, "kind": "стена",
                          "source": "routes/gtu/platform/batch/cost_estimate.toml [su2_time]: дымовой прогон core130 coarse "
                                    "14.11.2026, среднее по 169 итерациям history_s1_ramp60.csv, 4 процесса"},
    "mini_kv3_coarse_points": {"s_per_iter": 5.98, "nodes": 258299, "procs": 1, "kind": "CPU",
                               "source": "миниатюра этой задачи 14.11.2026 (контейнер gtu-platform:local, KV-3 coarse): "
                                         "2 точки x 2 процесса, MUSCL, 10 итераций — CPU-секунды процессов SU2_CFD "
                                         "(/proc, utime+stime) на итерацию точки: 5,98 с -> 23,1 мкс/узел"},
    "mini_kv3_coarse_base": {"s_per_iter": 16.5, "nodes": 258299, "procs": 1, "kind": "CPU",
                             "source": "та же миниатюра, база b1/b2 (первый порядок, 4 процесса): 16,3–16,6 CPU-с на "
                                       "итерацию -> 63–64 мкс/узел; завышено ожиданием MPI в активном опросе при "
                                       "загрузке 10 на 4 ядрах и первыми итерациями стадий — для сметы не берётся"},
}
C_EXPECTED = 35.2   # мкс: ref_cost_estimate (стена, 4 процесса) — больше замера точек миниатюры (23,1)
C_SIZING = 50.0     # мкс: запас для max_min (между 35,2 и 63 базы миниатюры); фактический расход — по факту счёта
EFF = 0.8           # доля идеального ускорения на 16–64 процессах (допущение; замер — узел bench)
KB_PER_NODE = 13.0  # память: миниатюра — 11,7–12,2 КБ/узел (VmRSS процессов SU2, база и 2 точки), с запасом
RESTART_B_PER_NODE = 104   # restart_c1_point.dat миниатюры: 26 863 545 байт / 258 299 узлов
SMALL_MB_PER_POINT = 40    # history, журналы, cfg, run_summary, surface.vtu (8,3 МБ на coarse KV-3; ∝ узлам^(2/3))
S3_USD_GB_MONTH = 0.023    # допущение (тариф S3 Standard us-east-1, не сверен здесь)
EGRESS_USD_GB = 0.09       # допущение: collect качает всё из префикса задачи в сессию вне AWS

PRICE_EXPECTED = 0.75     # долл./ч: задача владельца 14.11.2026 (спот 0,71–0,75); ledger smoke2: us-east-1c 0,7064
PRICE_UPPER = 1.4206      # MaxPrice раннера в ledger.jsonl (14.11.2026)
DISK_USD_GB_H = 0.08 / 730
IPV4_USD_H = 0.005
VCPU = 64
MEM_GB = 128

BOOT_MIN = 3        # загрузка, dnf docker, скачивание образа, клон (smoke2: ccx-задача целиком ~2 мин)
SETUP_MIN = 12      # apt + pip + сборка SU2 на 64 ядрах (оценка: сборка на 4 ядрах 8 мин, README узла su2)
MESH_MIN = {("kv3_stage1", "medium"): 8, ("core130", "medium"): 4, ("kv3_stage1", "fine"): 30, ("core130", "fine"): 12}
BASE_ITERS = 520    # план smoke075 run_case.py: 170 + 200 + 150
POINT_ITER = 6000
POINT_PRAMP = 300
RESERVE_MIN = 5 + 8  # запас обёртки (--reserve-s 300) + запас раннера на выгрузку (8 мин)


def c_us(name: str) -> float:
    t = TIMING[name]
    return t["s_per_iter"] * t["procs"] / t["nodes"] * 1e6


def plan_jobs(stage: str, fine_pfracs: list) -> list:
    """Список задач этапа: (линия, уровень, доли давления, процессов на точку)."""
    if stage == "medium":
        pf = [0.50, 0.60, 0.70, 0.80, 0.90, 1.00, 1.10, 1.20]
        a, b = pf[0::2], pf[1::2]        # чередование: потеря одной машины не выбивает край линии
        return [("kv3_stage1", "medium", a, 16, "a"), ("kv3_stage1", "medium", b, 16, "b"),
                ("core130", "medium", a, 16, "a"), ("core130", "medium", b, 16, "b")]
    out = []
    for i, p in enumerate(fine_pfracs):          # KV-3 fine ~5 млн узлов: точка на машину, все 64 процесса
        out.append(("kv3_stage1", "fine", [p], 64, chr(ord("a") + i)))
    half = (len(fine_pfracs) + 1) // 2           # core130 fine ~2,4 млн узлов: две точки на машину по 32
    out.append(("core130", "fine", fine_pfracs[0::2][:half], 32, "a"))
    out.append(("core130", "fine", fine_pfracs[1::2], 32, "b"))
    return [j for j in out if j[2]]


def job_time(line: str, level: str, n_points: int, np_pt: int, c: float = C_EXPECTED,
             point_iter: int = POINT_ITER) -> dict:
    n = NODES[(line, level)]
    t_base = c * 1e-6 * n / VCPU / EFF
    t_pt = c * 1e-6 * n / np_pt / EFF
    base_min = t_base * BASE_ITERS / 60
    pt_min = t_pt * point_iter / 60
    total = BOOT_MIN + SETUP_MIN + MESH_MIN[(line, level)] + base_min + pt_min + RESERVE_MIN
    out_gb = n_points * (n * RESTART_B_PER_NODE / 1e9 + SMALL_MB_PER_POINT / 1e3)
    return dict(nodes=n, s_per_iter_base=round(t_base, 3), s_per_iter_point=round(t_pt, 3), base_min=round(base_min, 1),
                points_min=round(pt_min, 1), total_min=round(total, 1),
                mem_gb_points=round(n * KB_PER_NODE * 1e3 * n_points / 1e9, 1), outputs_gb=round(out_gb, 2))


def command(line: str, level: str, pfracs: list, np_pt: int, max_min: int, mini: dict | None = None) -> str:
    L = LINES[line]
    lvl = mini["level"] if mini else level
    geom_out = L["geom_out"].format(level=lvl)
    args = [f"--case {L['case']}", f"--level {lvl}", f"--geom-out \"{geom_out}\"", f"--label {line}",
            "--pfracs " + ",".join(f"{p:.2f}" for p in pfracs), "--base-pfrac 0.75",
            f"--np-per-point {mini['np_pt'] if mini else np_pt}"]
    if mini:
        args += [f"--base-plan {mini['base_plan']}", f"--point-iter {mini['point_iter']}",
                 f"--point-pramp {mini['point_pramp']}", f"--time-limit {mini['time_limit']}",
                 f"--nproc {mini['nproc']}"]
    else:
        args += [f"--point-iter {POINT_ITER}", f"--point-pramp {POINT_PRAMP}", f"--job-max-min {max_min}"]
    args += ["--su2-src /opt/SU2-src", "--work \"$W/runs\"", "--out \"$GTU_OUT\""]
    run = "python routes/gtu/platform/su2/char_points.py run " + " ".join(args)
    return "\n".join(ENV_LINES + [L["mesh"].format(level=lvl), run])


def estimate(jobs: list, stage: str) -> dict:
    rows, mh, mh_sz, upper, disk, out_gb = [], 0.0, 0.0, 0.0, 0.0, 0.0
    for j in jobs:
        t = j["_time"]
        h = t["total_min"] / 60
        mh += h
        mh_sz += j["_time_sizing"]["total_min"] / 60
        out_gb += t["outputs_gb"]
        upper += (PRICE_UPPER + j["disk_gb"] * DISK_USD_GB_H + IPV4_USD_H) * (j["max_min"] + 10) / 60
        disk += j["disk_gb"] * DISK_USD_GB_H * h
        rows.append(dict(job=j["name"], line=j["_line"], level=j["_level"], pfracs=j["_pfracs"], np_per_point=j["_np"],
                         max_min=j["max_min"], disk_gb=j["disk_gb"], **t,
                         total_min_at_C_sizing=j["_time_sizing"]["total_min"]))
    s3, egress = out_gb * S3_USD_GB_MONTH, out_gb * EGRESS_USD_GB
    return {
        "stage": stage,
        "machines": len(jobs), "vcpu_per_machine": VCPU,
        "machine_hours_expected": round(mh, 2),
        "machine_hours_at_C_sizing": round(mh_sz, 2),
        "outputs_gb": round(out_gb, 2),
        "usd": {"spot_expected": round(mh * PRICE_EXPECTED, 2), "disk": round(disk, 3),
                "ipv4": round(mh * IPV4_USD_H, 3), "s3_month": round(s3, 3), "transfer_to_session": round(egress, 3),
                "total_expected": round(mh * (PRICE_EXPECTED + IPV4_USD_H) + disk + s3 + egress, 2),
                "total_at_C_sizing": round(mh_sz * (PRICE_EXPECTED + IPV4_USD_H) + s3 + egress, 2)},
        "upper_bound_usd": round(upper, 2),
        "upper_bound_rule": "раннер: (MaxPrice 1,4206 + диск + IPv4) × (max_min + 10) / 60 на машину — если каждая "
                            "задача дойдёт до сторожа",
        "jobs": rows,
        "time_model": {
            "formula": "t_итерации = C · узлы / процессы / EFF; машина = загрузка + подготовка + сетка + база (520 итераций "
                       "на 64 процессах) + точки (point_iter итераций) + запасы",
            "C_expected_us": C_EXPECTED, "C_sizing_us": C_SIZING,
            "C_measured_us": {k: round(c_us(k), 1) for k in TIMING}, "measurements": TIMING, "EFF": EFF,
            "boot_min": BOOT_MIN, "setup_min": SETUP_MIN, "mesh_min": {f"{k[0]}/{k[1]}": v for k, v in MESH_MIN.items()},
            "base_iters": BASE_ITERS, "point_iter": POINT_ITER, "reserve_min": RESERVE_MIN,
            "nodes": {f"{k[0]}/{k[1]}": v for k, v in NODES.items()}, "nodes_source": NODES_SRC,
            "fine_nodes_rule": "fine — оценка: узлы medium² / узлы coarse (рост сетки при переходе уровня — тот же)",
            "notes": ["замеры C — на общей 4-ядерной машине (Xeon 2,1 ГГц, загрузка 5–16); для c7a.16xlarge (свои 64 "
                      "ядра) C, скорее всего, меньше; EFF — допущение; max_min — по C_sizing, ожидаемая смета — по "
                      "C_expected, верхняя граница — по max_min",
                      "выход: рестарт итоговый на точку (104 Б/узел, замер миниатюры) + ~40 МБ мелочи; runner.py "
                      "collect качает всё из префикса задачи в сессию — outputs_gb должно поместиться на её диск",
                      "boot, setup, mesh — оценки; сборка SU2 на 64 ядрах не мерилась; сетки medium и fine на AWS не строились",
                      "точка идёт до point_iter или до срока: срок обёртка считает от max_min и /proc/uptime, SU2 "
                      "останавливается сам и пишет рестарт"],
        },
        "price": {"expected": PRICE_EXPECTED, "upper": PRICE_UPPER,
                  "basis": "ожидаемая — задача владельца 14.11.2026 (спот c7a.16xlarge 0,71–0,75 долл./ч; ledger smoke2 "
                           "20261114T231610Z-9b64: us-east-1c 0,7064); верхняя — MaxPrice раннера в ledger.jsonl"},
        "disk_rate": "0,08 долл./ГБ·мес / 730 ч и IPv4 0,005 долл./ч — ставки раннера (допущения кода)",
        "memory_check": f"память точек на машину ≤ {MEM_GB} ГБ при {KB_PER_NODE} КБ/узел (миниатюра 11,7–12,2 с запасом) "
                        "— поле mem_gb_points; база на 64 процессах идёт одна, до точек",
        "s3_transfer_rates": "S3 0,023 долл./ГБ·мес, выход из AWS 0,09 долл./ГБ — допущения, не сверены",
    }


def build(stage: str, fine_pfracs: list) -> dict:
    jobs = []
    for line, level, pf, np_pt, tag in plan_jobs(stage, fine_pfracs):
        t = job_time(line, level, len(pf), np_pt)
        tz = job_time(line, level, len(pf), np_pt, c=C_SIZING)
        max_min = min(720, int(math.ceil(tz["total_min"] / 30.0) * 30))
        if t["total_min"] > 720 - 10:
            raise SystemExit(f"{line}/{level}: оценка {t['total_min']} мин не влезает в 720")
        disk = 50 if level == "medium" else 80
        jobs.append({"name": f"su2-{line.replace('_', '-')}-{level}-{tag}",
                     "command": command(line, level, pf, np_pt, max_min), "max_min": max_min, "disk_gb": disk,
                     "_time": t, "_time_sizing": tz, "_line": line, "_level": level, "_pfracs": pf, "_np": np_pt})
    est = estimate(jobs, stage)
    out_json = f"routes/gtu/platform/su2/su2_char_kv3_{stage}.json"
    doc = {
        "_comment": (f"Этап {'A' if stage == 'medium' else 'B'} ({stage}) линий частоты SU2: KV-3 ступень 1 и core130. Раннер берёт "
                     "только defaults и jobs; блоки estimate, merge, plan_check — для сессии с ключами. Код — ветка "
                     f"{BRANCH}. Сгенерировано make_jobs_su2_kv3.py; правки — в нём."
                     + ("" if stage == "medium" else " Запускать после разбора этапа A; точки — по его итогам "
                        "(make_jobs_su2_kv3.py --fine-pfracs).")),
        "defaults": {"instance_type": "c7a.16xlarge", "fallback_types": ["c7i.16xlarge"], "image": "al2023",
                     "container": CONTAINER, "container_setup": CONTAINER_SETUP},
        "jobs": [{k: v for k, v in j.items() if not k.startswith("_")} for j in jobs],
        "lines": {k: v["about"] for k, v in LINES.items()},
        "setup": {"base": "линия сетки -> база: план smoke075 run_case.py на 64 процессах при p_вых 0,75 оценки "
                          "(как прогон Б, README узла) -> точки: рестарт от базы, рампа p_вых за 300 итераций, всего "
                          f"{POINT_ITER} итераций второго порядка; SST, покрышка и стенки за колесом без трения "
                          "(--shroud slip --vaneless slip — как дымовой прогон; с неподвижной покрышкой на сетках geom "
                          "расчёт расходился)",
                  "pfracs": "доли оценки статического давления на выходе (run_case.estimate_p_out; у KV-3 — по расходу "
                            "стороны); у задачи — свои, чередованием",
                  "outputs_per_job": "char_summary.json, env.txt, geom.log, frame3_geom_results.json (KV-3), "
                                     "_su2_build.log, base/ и p<доля>/: run_summary.json, history_*.csv, su2_*.log, "
                                     "case_*.cfg, run_case.out, surface.vtu (лопатка и втулка, SURFACE_PARAVIEW), итоговый рестарт точки "
                                     "restart_c1_point.dat; сетки и объёмные поля не выгружаются (поле — SU2_SOL из рестарта)"},
        "merge": {
            "where": "сессия с ключами после run: результаты уже в routes/gtu/platform/aws/runs/<run_id>/<задача>/ "
                     "(или runner.py collect --run-id <run_id> --job <задача> по каждой)",
            "commands": [f"python3 routes/gtu/platform/su2/char_points.py merge --in routes/gtu/platform/aws/runs/<run_id> "
                         f"--out {out_json} --run-id <run_id>"],
            "to_branch": [out_json, "routes/gtu/platform/aws/ledger.jsonl"],
            "not_to_branch": "runs/<run_id>/ (рестарты точек: ~100 МБ на точку medium KV-3, ~0,5 ГБ fine — остаются в S3 "
                             "и в runs/, папка в .gitignore)",
            "checks": "merge откажет (код 2), если точка задана дважды, пришла не своей задаче или сводка > 5 МБ; "
                      "missing — заданные, но не пришедшие точки; converged — π и η меняются < 0,1 % за 500 итераций",
        },
        "estimate": est,
    }
    return doc


def smoke_command() -> str:
    """Этап 0: обе линии на coarse одной машиной; итерации — десятки (план smoke075 x 0,1 = 52, точки по 40);
    код задачи — merge --require-base-ok (база каждой линии дошла до конца, рестарт есть)."""
    lines = list(ENV_LINES)
    for line in LINES:
        L = LINES[line]
        lines.append(L["mesh"].format(level="coarse").replace("geom.log", f"geom_{line}.log"))
        args = [f"--case {L['case']}", "--level coarse", f"--geom-out \"{L['geom_out'].format(level='coarse')}\"",
                f"--label {line}", "--pfracs 0.70,0.90,1.00,1.10", "--base-pfrac 0.75", "--np-per-point 16",
                "--base-iter-scale 0.1", "--point-iter 40", "--point-pramp 10", f"--job-max-min {SMOKE_MAX_MIN}",
                "--su2-src /opt/SU2-src", f"--work \"$W/runs_{line}\"", f"--out \"$GTU_OUT/{line}\""]
        lines.append("python routes/gtu/platform/su2/char_points.py run " + " ".join(args)
                     + f" || echo \"{line} exit $?\" >> \"$GTU_OUT/smoke_exits.txt\"")
    lines.append("python routes/gtu/platform/su2/char_points.py merge --in \"$GTU_OUT\" "
                 "--out \"$GTU_OUT/smoke_merge.json\" --require-base-ok --expect " + ",".join(LINES))
    return "\n".join(lines)


def build_smoke() -> dict:
    disk = 40
    t_min = BOOT_MIN + SETUP_MIN + 2 * 5 + 2 * 3     # две сетки coarse (~2 мин здесь) и два прогона по ~3 мин
    upper = (PRICE_UPPER + disk * DISK_USD_GB_H + IPV4_USD_H) * (SMOKE_MAX_MIN + 10) / 60
    job = {"name": "su2-kv3-core130-smoke", "command": smoke_command(), "max_min": SMOKE_MAX_MIN, "disk_gb": disk}
    return {
        "_comment": ("Этап 0 линий частоты SU2: проба всей цепочки на одной машине до этапов A и B (обе линии на coarse, "
                     "итераций — десятки, числа не характеристика). Раннер берёт только defaults и jobs. Код — ветка "
                     f"{BRANCH}. Сгенерировано make_jobs_su2_kv3.py; правки — в нём."),
        "defaults": {"instance_type": "c7a.16xlarge", "fallback_types": ["c7i.16xlarge"], "image": "al2023",
                     "container": CONTAINER, "container_setup": CONTAINER_SETUP},
        "jobs": [job],
        "what_it_checks": ["container_setup на машине: apt, venv, pip по requirements.txt, сборка SU2 v8.6.0 с GitHub "
                           "(_su2_build.log, _container_setup.log, время — [gtu-time] container_setup в _job.log)",
                           "окружение команды (OMPI_MCA_osc=^ucx, vader) — SU2 под mpirun в docker: база на 64 процессах, "
                           "4 точки x 16 параллельно, --use-hwthread-cpus",
                           "сетки на машине: frame3_geom.py (KV-3, STEP и проверки) и impeller.py (core130)",
                           "срок от max_min и /proc/uptime (budget.uptime_at_start_s в char_summary.json)",
                           "выгрузка и сбор: char_summary.json, рестарты, surface.vtu; merge --require-base-ok -> код задачи",
                           "время подготовки и итерации на c7a.16xlarge — для уточнения сметы этапов A и B "
                           "(sec_per_iter в char_summary.json, [gtu-time] в _job.log)"],
        "merge": {"commands": ["python3 routes/gtu/platform/su2/char_points.py merge --in "
                               "routes/gtu/platform/aws/runs/<run_id> --out <scratchpad>/su2_char_kv3_smoke.json "
                               "--run-id <run_id> --require-base-ok"],
                  "to_branch": ["routes/gtu/platform/aws/ledger.jsonl"],
                  "note": "числа пробы в ветку не идут; время подготовки и итерации — в смету этапов A и B "
                          "(make_jobs_su2_kv3.py: SETUP_MIN, TIMING)"},
        "estimate": {"machines": 1, "machine_min_expected": t_min,
                     "usd": {"total_expected": round(t_min / 60 * (PRICE_EXPECTED + IPV4_USD_H + disk * DISK_USD_GB_H), 2)},
                     "upper_bound_usd": round(upper, 2),
                     "basis": "загрузка 3 + подготовка 12 + 2 сетки coarse по ~5 + 2 прогона по ~3 мин (оценки; сетка KV-3 "
                              "coarse здесь — 131 с, core130 coarse — 40 с); цена — как этапы A и B"},
    }


def plan_check(path: str) -> dict:
    """runner.py plan по файлу (без AWS, цены — допущение); журнал проверок plan — во временную папку."""
    args = ["plan", os.path.relpath(path, REPO), "--assume-price", "0.748",
            "--assume-zone-prices", "us-east-1a=0.955,us-east-1b=0.748,us-east-1c=0.706", "--branch", BRANCH]
    out = tempfile.mkdtemp(prefix="gtu_plan_")
    r = subprocess.run([sys.executable, RUNNER] + args + ["--out", out], capture_output=True, text=True, cwd=REPO,
                       timeout=300)
    shutil.rmtree(out, ignore_errors=True)
    keep = [ln for ln in r.stdout.splitlines() if ln.startswith(("[задача]", "[бюджет]", "[зоны]", "[план]", "[очередь]"))]
    ud = [ln.strip().rstrip(",") for ln in r.stdout.splitlines() if '"UserData"' in ln]
    return {"command": " ".join(["python3", os.path.relpath(RUNNER, REPO)] + args), "exit_code": r.returncode,
            "lines": keep, "userdata_sizes": ud, "stderr": r.stderr.strip()[-2000:]}


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--fine-pfracs", default="0.60,0.80,0.90,1.00")
    ap.add_argument("--no-plan", action="store_true")
    ap.add_argument("--mini", default=None, metavar="DIR", help="записать миниатюры команд в DIR и выйти")
    a = ap.parse_args(argv)
    if a.mini:
        os.makedirs(a.mini, exist_ok=True)
        mini = dict(level="coarse", np_pt=2, nproc=4, base_plan="/gtu/mini/plan_base_mini.json", point_iter=6,
                    point_pramp=3, time_limit=1500)
        with open(os.path.join(a.mini, "plan_base_mini.json"), "w", encoding="utf-8") as f:
            json.dump([dict(name="b1", frac=0.6, pfrac=0.95, ramp=2, pramp=3, iter=4),
                       dict(name="b2", frac=1.0, pfrac=0.75, ramp=2, pramp=3, iter=4),
                       dict(name="b3", frac=1.0, pfrac=0.75, ramp=0, iter=3, muscl=True)], f)
        with open(os.path.join(a.mini, "cmd_smoke_aws.sh"), "w", encoding="utf-8") as f:
            f.write(smoke_command() + "\n")
        for line in LINES:
            txt = command(line, "medium", [0.70, 0.90], 16, 360, mini)
            with open(os.path.join(a.mini, f"cmd_{line}.sh"), "w", encoding="utf-8") as f:
                f.write(txt + "\n")
        with open(os.path.join(a.mini, "container_setup.sh"), "w", encoding="utf-8") as f:
            f.write("\n".join(CONTAINER_SETUP) + "\n")
        print("миниатюры:", a.mini)
        return 0
    fine = [float(x) for x in a.fine_pfracs.split(",")]
    for stage, path in (("smoke", OUT_SMOKE), ("medium", OUT_MEDIUM), ("fine", OUT_FINE)):
        doc = build_smoke() if stage == "smoke" else build(stage, fine)
        with open(path, "w", encoding="utf-8") as f:
            json.dump(doc, f, ensure_ascii=False, indent=1)
        if not a.no_plan:
            doc["plan_check"] = plan_check(path)
            with open(path, "w", encoding="utf-8") as f:
                json.dump(doc, f, ensure_ascii=False, indent=1)
        e = doc["estimate"]
        print(f"{os.path.basename(path)}: задач {len(doc['jobs'])}, машино-часов {e.get('machine_hours_expected')}, "
              f"ожидается {e['usd']['total_expected']} долл., верхняя граница {e['upper_bound_usd']} долл."
              + (f"; plan: код {doc['plan_check']['exit_code']}" if "plan_check" in doc else ""))
    return 0


if __name__ == "__main__":
    sys.exit(main())
