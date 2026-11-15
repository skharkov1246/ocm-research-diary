#!/usr/bin/env python3
"""Файл задач раннера для замера скорости машин AWS (узел bench): jobs_aws_bench.json со сметой и проверкой plan.

Три задачи (по машине на каждую):
  bench-su2-c7a16 / bench-su2-c7i16 — SU2 v8.6.0 (сборка в контейнере ubuntu:24.04, как задачи узла su2):
      speed_su2.py на канале core130 — medium (64x1, 32x2, 16x4) и fine (64x1, 32x2); c7i.16xlarge — запасной тип
      раннера (64 vCPU = 32 ядра × 2 потока), замер показывает, во сколько обойдётся точка, если c7a не дадут;
  bench-ccx-c7a16 — CalculiX 2.21 (пакет Ubuntu) в контейнере: speed_ccx.py, колесо core130 C3D10 на 50, 150 и 450 тыс.
      узлов, статика при 1, 8, 64 потоках, упаковка 8 копий × 8 потоков, частоты до 250 тыс. узлов.
Слияние — speed_merge.py → speed_aws_results.json (цена точки SU2 medium/fine и расчёта ccx).

Раннер читает только defaults и jobs; блоки estimate, merge, plan_check — для сессии с ключами.
  python3 routes/gtu/platform/bench/make_jobs_bench.py [--plan-check] [--mini-dir ПАПКА]
--mini-dir пишет те же команды с малыми размерами (cmd_<задача>.sh) для локальной пробы в контейнере.
"""
from __future__ import annotations

import argparse
import json
import os
import re
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", "..", ".."))
OUT = os.path.join(HERE, "jobs_aws_bench.json")
REL = "routes/gtu/platform/bench/jobs_aws_bench.json"
BRANCH = "claude/lucid-allen-yvps5d"

# Сборка SU2 в контейнере — тот же список, что defaults.container_setup в routes/gtu/platform/su2/jobs_aws_su2_kv3.json
# (узел su2): одинаковая сборка — сравнимые замеры.
SU2_SETUP = [
    "apt-get update -q",
    "apt-get install -y -q --no-install-recommends build-essential python3 python3-venv git pkg-config ca-certificates "
    "openmpi-bin libopenmpi-dev libgl1 libglu1-mesa libxrender1 libxcursor1 libxfixes3 libxft2 libfontconfig1 "
    "libxinerama1 libgomp1",
    "python3 -m venv /opt/venv",
    "/opt/venv/bin/pip install -q --no-cache-dir -r routes/gtu/platform/docker/requirements.txt",
    "/opt/venv/bin/python -c 'import cantera, numpy, scipy, gmsh; print(\"py ok\", cantera.__version__, gmsh.__version__)'",
    "JOBS=\"$(nproc)\" bash routes/gtu/platform/su2/install_su2.sh /opt/su2 /opt/SU2-src > \"$GTU_OUT/_su2_build.log\" 2>&1 "
    "|| { tail -n 40 \"$GTU_OUT/_su2_build.log\"; exit 1; }",
    "test \"$(git -C /opt/SU2-src rev-parse HEAD)\" = 7048fd35f56d4bf52e586429d6746a36a10255f8",
    "/opt/su2/bin/SU2_CFD --help | head -n 1 | tee /opt/su2/VERSION",
]

# --use-hwthread-cpus: OpenMPI 4.1 по умолчанию считает слоты по физическим ядрам; у c7i.16xlarge их 32 при 64 vCPU,
# и «mpirun -np 64» без ключа отказывает («not enough slots»; проверено локально на 4 ядрах: -np 5 — отказ). У c7a.16xlarge
# 64 ядра без SMT — ключ ничего не меняет.

# CalculiX из пакета Ubuntu 24.04 (calculix-ccx 2.21-1) + Python-пакеты образа (gmsh для сетки колеса).
CCX_SETUP = [
    "apt-get update -q",
    "apt-get install -y -q --no-install-recommends calculix-ccx python3 python3-venv git ca-certificates libgl1 "
    "libglu1-mesa libxrender1 libxcursor1 libxfixes3 libxft2 libfontconfig1 libxinerama1 libgomp1",
    "python3 -m venv /opt/venv",
    "/opt/venv/bin/pip install -q --no-cache-dir -r routes/gtu/platform/docker/requirements.txt pytest==9.2.0",
    "/opt/venv/bin/python -c 'import numpy, scipy, gmsh; print(\"py ok\", numpy.__version__, gmsh.__version__)'",
    "{ ccx -v || true; } > /opt/ccx_version.txt 2>&1",
    "grep -q 'Version 2.21' /opt/ccx_version.txt",
]

ENV_SU2 = ("export PATH=\"/opt/venv/bin:/opt/su2/bin:$PATH\" SU2_RUN=/opt/su2/bin PYTHONDONTWRITEBYTECODE=1\n"
           "export OMPI_MCA_btl_vader_single_copy_mechanism=none OMPI_MCA_osc=^ucx")
ENV_CCX = "export PATH=\"/opt/venv/bin:$PATH\" PYTHONDONTWRITEBYTECODE=1"
ENV_TXT = ("{{ {ver}; python -V; echo \"nproc $(nproc) NPROC $NPROC itype ${{GTU_ITYPE:-}} az ${{GTU_AZ_MACHINE:-}}\"; "
           "lscpu | head -n 20; free -g; df -h /tmp; }} > \"$GTU_OUT/env.txt\" 2>&1 || true\n"
           "git rev-parse HEAD > \"$GTU_OUT/commit.txt\" 2>/dev/null || true")


def su2_cmd(levels, configs, warm, iters, deadline_args, su2_src="/opt/SU2-src"):
    cfg = " ".join(f"--config {lv}={','.join(c)}" for lv, c in configs.items())
    return "\n".join([
        "set -e", ENV_SU2, ENV_TXT.format(ver="cat /opt/su2/VERSION; mpirun --version | head -n 1"),
        f"python routes/gtu/platform/bench/speed_su2.py run --case core130 --levels {','.join(levels)} {cfg} "
        f"--warm {warm} --iters {iters} --nproc \"$NPROC\" --mpi-args=--use-hwthread-cpus {deadline_args} "
        f"--su2-src {su2_src} "
        "--work /tmp/gtu_speed_su2 --out \"$GTU_OUT\""])


def ccx_cmd(targets, threads, pack, freq_max, deadline_args):
    return "\n".join([
        "set -e", ENV_CCX, ENV_TXT.format(ver="cat /opt/ccx_version.txt"),
        f"python routes/gtu/platform/bench/speed_ccx.py run --case core130 --targets {','.join(map(str, targets))} "
        f"--threads {threads} --pack {pack} --freq-max-nodes {freq_max} --nproc \"$NPROC\" {deadline_args} "
        "--work /tmp/gtu_speed_ccx --out \"$GTU_OUT\""])


SU2_LEVELS = ["medium", "fine"]
SU2_CONFIGS = {"medium": ["64x1", "32x2", "16x4"], "fine": ["64x1", "32x2"]}
SU2_MAX_MIN = 90
CCX_TARGETS = [50000, 150000, 450000]
CCX_MAX_MIN = 100
PRICE, PRICE_UPPER = 0.75, 1.4206
DISK_RATE, IPV4 = 0.08 / 730, 0.005

# Модель времени задачи SU2 (для сметы, до замера): t_итер = C·узлы/NP/EFF — как в узле su2
# (routes/gtu/platform/su2/jobs_aws_su2_kv3.json → estimate.time_model).
C_US, EFF = 35.24, 0.8
NODES = {"medium": 584357, "fine": 2455000}
STAGE_OVERHEAD_S = {"medium": 60, "fine": 120}   # запуск, чтение и разбиение сетки, запись рестарта — оценка
MESH_MIN = {"medium": 4, "fine": 12}              # как mesh_min core130 в узле su2 — оценка
BOOT_SETUP = dict(su2=15.0, ccx=3.0)              # мин: загрузка + контейнер + пакеты (+ сборка SU2 12 мин — оценка узла su2)


def su2_estimate():
    rows, total = [], BOOT_SETUP["su2"]
    for lv in SU2_LEVELS:
        total += MESH_MIN[lv]
        for c in SU2_CONFIGS[lv]:
            np_, k = map(int, c.split("x"))
            s = C_US * 1e-6 * NODES[lv] / np_ / EFF
            t = (40 + 100) * s + 2 * STAGE_OVERHEAD_S[lv]
            rows.append(dict(level=lv, cfg=c, s_iter_model=round(s, 3), min=round(t / 60, 1)))
            total += t / 60
    total += 1.0  # выгрузка
    return rows, round(total, 1)


def ccx_estimate():
    """Время ccx по локальному замеру (контейнер сессии, 4 общих ядра под нагрузкой ~11, 15.11.2026):
    статика 1 поток — 3,52 с при 12 100 узлах и 8,6 с при 19 016 узлах → t ∝ N^1,8 (две точки, вне диапазона —
    экстраполяция; на своей машине c7a ожидается быстрее). Многопоточность SPOOLES — допущение: 64 потока в 3 раза
    быстрее одного."""
    a, b = 8.6 / 19016 ** 1.8, 1.8
    rows, total = [], BOOT_SETUP["ccx"] + 0.5
    for n in CCX_TARGETS:
        t1 = a * n ** b
        t64 = t1 / 3
        mesh = 0.5 + n / 200000  # мин, оценка (Gmsh HXT)
        parts = dict(mesh_min=mesh, t1_min=t1 / 60, t8_min=t1 / 2 / 60, t64_min=t64 / 60, pack_min=t1 / 2 / 60 * 1.5,
                     freq_min=(4 * t64 / 60) if n <= 250000 else 0)
        # как пропуск в speed_ccx.py: на большой сетке 1 поток и упаковка по прогнозу не успевают — не считаются
        if n > 200000:
            parts["t1_min"] = 0.0
            parts["pack_min"] = 0.0
        m = sum(parts.values())
        rows.append(dict(nodes=n, **{k: round(v, 1) for k, v in parts.items()}, min=round(m, 1)))
        total += m
    total += 1.0
    return rows, round(total, 1)


def job_cost(minutes, max_min, disk_gb):
    ex = disk_gb * DISK_RATE + IPV4
    return dict(machine_min_expected=minutes, usd_expected=round(minutes / 60 * (PRICE + ex), 3),
                usd_upper=round((max_min + 10) / 60 * (PRICE_UPPER + ex), 3))


def build():
    su2_rows, su2_min = su2_estimate()
    ccx_rows, ccx_min = ccx_estimate()
    jobs = [
        dict(name="bench-su2-c7a16", instance_type="c7a.16xlarge", container_setup=SU2_SETUP, max_min=SU2_MAX_MIN,
             command=su2_cmd(SU2_LEVELS, SU2_CONFIGS, 40, 100, f"--job-max-min {SU2_MAX_MIN}")),
        dict(name="bench-su2-c7i16", instance_type="c7i.16xlarge", container_setup=SU2_SETUP, max_min=SU2_MAX_MIN,
             command=su2_cmd(SU2_LEVELS, SU2_CONFIGS, 40, 100, f"--job-max-min {SU2_MAX_MIN}")),
        dict(name="bench-ccx-c7a16", instance_type="c7a.16xlarge", container_setup=CCX_SETUP, max_min=CCX_MAX_MIN,
             command=ccx_cmd(CCX_TARGETS, "1,8,64", 8, 250000, f"--job-max-min {CCX_MAX_MIN}")),
    ]
    est_jobs = {
        "bench-su2-c7a16": job_cost(su2_min, SU2_MAX_MIN, 60),
        "bench-su2-c7i16": job_cost(su2_min, SU2_MAX_MIN, 60),
        "bench-ccx-c7a16": job_cost(ccx_min, CCX_MAX_MIN, 60),
    }
    tot_e = round(sum(v["usd_expected"] for v in est_jobs.values()), 2)
    tot_u = round(sum(v["usd_upper"] for v in est_jobs.values()), 2)
    doc = {
        "_comment": ("Замер скорости машин для сметы уровней 2–3 воронки (узел bench): SU2 medium/fine на c7a.16xlarge и "
                     "c7i.16xlarge, CalculiX 3D на c7a.16xlarge. Раннер берёт только defaults и jobs; блоки estimate, merge, "
                     f"plan_check — для сессии с ключами. Код — ветка {BRANCH}. Сгенерировано make_jobs_bench.py; правки — в нём. "
                     "Запускать до длинных серий SU2 узла su2: замер уточняет их max_min и смету."),
        "defaults": {"instance_type": "c7a.16xlarge", "fallback_types": [], "image": "al2023",
                     "container": "ubuntu:24.04", "disk_gb": 60},
        "jobs": jobs,
        "what": {
            "bench-su2-c7a16": "секунды на итерацию канала core130 (medium 584 357 узлов, fine — уровень impeller.py) при "
                               "64x1, 32x2, 16x4 (medium) и 64x1, 32x2 (fine); время сетки; память SU2",
            "bench-su2-c7i16": "то же на запасном типе раннера (c7i: 64 vCPU = 32 ядра × 2 потока)",
            "bench-ccx-c7a16": "CalculiX: колесо core130, C3D10 на 50/150/450 тыс. узлов — статика CENTRIF при 1/8/64 потоках, "
                               "упаковка 8×8, частоты (10 форм) до 250 тыс. узлов; время, память, max|U| для сверки",
            "fallback_types": "пусто нарочно: замер должен идти на заданном типе (иначе сравнение типов теряет смысл)",
        },
        "merge": {
            "where": "сессия с ключами после run: результаты в routes/gtu/platform/aws/runs/<run_id>/<задача>/ "
                     "(или runner.py collect --run-id <run_id> --job <задача>)",
            "commands": ["python3 routes/gtu/platform/bench/speed_merge.py --in routes/gtu/platform/aws/runs/<run_id> "
                         "--run-id <run_id> --out routes/gtu/platform/bench/speed_aws_results.json"],
            "to_branch": ["routes/gtu/platform/bench/speed_aws_results.json", "routes/gtu/platform/aws/ledger.jsonl"],
            "not_to_branch": "runs/<run_id>/ целиком (speed_*.json задач — копия в сводке с sha256; журналы SU2 и ccx, сетки "
                             "и поля в S3 не выгружаются вовсе: остаются в /tmp машины)",
            "checks": ("speed_merge.py откажет (код 2) без замеров и не запишет сводку больше 5 МБ; в сводке — max_disp_rel_spread "
                       "(разброс max|U| статики при разных потоках: > 1e-6 — многопоточный SPOOLES дал другое решение, "
                       "см. disk.py: однажды 4 потока дали невязку равновесия 1,6 %), C и с/итер по каждой конфигурации; "
                       "в README узла — таблица «модель узла su2 против замера»"),
            "after": ("пересчитать max_min и смету задач узла su2 (make_jobs_su2_kv3.py — C_us и EFF) по C и с/итер "
                      "замера; без замера их смета держится на C = 35,24 мкс с общей 4-ядерной машины"),
        },
        "estimate": {
            "price": {"expected": PRICE, "upper": PRICE_UPPER,
                      "basis": ("ожидаемая — задача владельца 14.11.2026 (спот c7a.16xlarge 0,71–0,75 долл./ч; ledger smoke2 "
                                "20261114T231610Z-9b64 us-east-1c 0,7064); для c7i взята та же — живой цены c7i нет; "
                                "верхняя — MaxPrice раннера 1,4206 в ledger.jsonl")},
            "jobs": est_jobs,
            "total_usd_expected": tot_e,
            "total_usd_upper": tot_u,
            "upper_bound_rule": "раннер: (MaxPrice + диск + IPv4) × (max_min + 10) / 60 на машину — если задача дойдёт до сторожа",
            "machines_parallel": 3,
            "su2_time_model": dict(formula="t_итер = C·узлы/NP/EFF; задача = загрузка+сборка SU2 + сетки + Σ(140 итераций·t_итер + "
                                           "2 стадии × накладные) + выгрузка", C_us=C_US, EFF=EFF, nodes=NODES,
                                   stage_overhead_s=STAGE_OVERHEAD_S, mesh_min=MESH_MIN, boot_setup_min=BOOT_SETUP["su2"],
                                   rows=su2_rows, total_min=su2_min,
                                   source=("C и EFF — модель узла su2 (C из cost_estimate.toml: core130 coarse, 4 общих ядра, "
                                           "1,2255 с/итер); загрузка+сборка 15 мин — оценка узла su2 (сборка SU2 на 64 ядрах "
                                           "не мерилась); накладные стадии и сетки — оценки")),
            "ccx_time_model": dict(rows=ccx_rows, total_min=ccx_min, source=ccx_estimate.__doc__.strip(),
                                   boot_setup_min=BOOT_SETUP["ccx"],
                                   boot_setup_source=("ledger.jsonl: smoke-ccx (20261114T231610Z-9b64, контейнер ubuntu:24.04 + "
                                                      "apt calculix-ccx + disk.py) — 2,05 мин машины от запуска до гашения "
                                                      "(опрос раз в 60 с); сверху pip пакетов образа ~1 мин — оценка")),
            "guard": "каждый скрипт сам держит срок (--job-max-min): что не успевает — пропускается с пометкой, итог выгружается",
            "share_of_ceiling": "из 2000 долл. предохранителя раннера (уровни 0–1 — до 1500): ожидаемо ~2–3, верхняя граница ~7",
        },
    }
    return doc


def plan_check(path):
    cmd = [sys.executable, os.path.join(ROOT, "routes/gtu/platform/aws/runner.py"), "plan", path, "--assume-price", "0.748",
           "--assume-zone-prices", "us-east-1a=0.955,us-east-1b=0.748,us-east-1c=0.706",
           "--out", os.environ.get("GTU_PLAN_OUT", "/tmp/gtu_plan_bench")]
    p = subprocess.run(cmd, capture_output=True, text=True, cwd=ROOT)
    lines = [ln for ln in p.stdout.splitlines() if re.match(r"\[(задача|бюджет|зоны|план)\]", ln)]
    return dict(command=" ".join(["python3", "routes/gtu/platform/aws/runner.py", "plan", REL] + cmd[4:8]),
                exit_code=p.returncode, lines=lines,
                runner_note="раннер v2 (container, container_setup, NPROC); plan не вызывает AWS, цены — допущение")


def write_mini(d):
    """Те же команды с малыми размерами — для пробы в контейнере сессии (не для AWS)."""
    os.makedirs(d, exist_ok=True)
    mini = {
        "bench-su2": su2_cmd(["coarse"], {"coarse": ["2x1", "2x2"]}, 4, 4, "--budget-s 900"),
        "bench-ccx": ccx_cmd([8000, 20000], "1,2", 1, 25000, "--budget-s 600"),
    }
    for k, v in mini.items():
        open(os.path.join(d, f"cmd_{k}.sh"), "w").write(v + "\n")
    open(os.path.join(d, "setup_su2.sh"), "w").write("\n".join(SU2_SETUP) + "\n")
    open(os.path.join(d, "setup_ccx.sh"), "w").write("\n".join(CCX_SETUP) + "\n")
    return sorted(mini)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--plan-check", action="store_true", help="прогнать runner.py plan и записать итог в файл задач")
    ap.add_argument("--mini-dir", default="")
    ap.add_argument("--out", default=OUT)
    a = ap.parse_args(argv)
    doc = build()
    with open(a.out, "w", encoding="utf-8") as f:
        json.dump(doc, f, ensure_ascii=False, indent=1)
        f.write("\n")
    if a.plan_check:
        doc["plan_check"] = plan_check(a.out)
        with open(a.out, "w", encoding="utf-8") as f:
            json.dump(doc, f, ensure_ascii=False, indent=1)
            f.write("\n")
        print("\n".join(doc["plan_check"]["lines"]))
        print("plan exit", doc["plan_check"]["exit_code"])
    if a.mini_dir:
        print("mini:", write_mini(a.mini_dir))
    e = doc["estimate"]
    print(f"записано {a.out}: задач {len(doc['jobs'])}, ожидаемо {e['total_usd_expected']} долл., верх {e['total_usd_upper']}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
