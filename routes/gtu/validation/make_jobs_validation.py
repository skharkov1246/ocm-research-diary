#!/usr/bin/env python3
"""Файл задач раннера направления «проверка модели»: jobs_aws_validation.json со сметой и проверкой plan.

Две задачи — повтор эталонных расчётов трека на чистой машине AWS и сверка с закоммиченными результатами
(compare_results.py; правило «каждая правка данных проверяется дважды разными способами»):
  val-cycle — al2023 + Python 3.12 (uv) + cantera 3.2.0, numpy 2.5.3, scipy 1.18.1: тесты узла (pytest), balances.py,
      validate_fleet.py на всех ядрах (workers = $NPROC; на 4 общих ядрах — 8–19 мин); сверка balances_results.json и
      validation_fleet_results.json с веткой;
  val-fea — контейнер ubuntu:24.04 + calculix-ccx 2.21: test_disk.py (12 проверок против аналитики), disk.py --examples
      и frame3_fea.py (73 расчёта ccx, по расчёту на ядро); сверка disk_results.json и frame3_fea_results.json с веткой.
Эталонных кейсов SU2 (NASA Rotor 37, CC3, SRV2) здесь нет: геометрии и сеток в репозитории нет (навык su2-turbo,
раздел 9) — задача не собирается из того, чего нет.

Раннер читает только defaults и jobs; блоки estimate, merge, plan_check — для сессии с ключами.
  python3 routes/gtu/validation/make_jobs_validation.py [--plan-check] [--mini-dir ПАПКА]
"""
from __future__ import annotations

import argparse
import json
import os
import re
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
OUT = os.path.join(HERE, "jobs_aws_validation.json")
REL = "routes/gtu/validation/jobs_aws_validation.json"
BRANCH = "claude/lucid-allen-yvps5d"
sys.path.insert(0, os.path.join(ROOT, "routes", "gtu", "platform", "bench"))
from make_jobs_bench import CCX_SETUP, ENV_CCX, PRICE, PRICE_UPPER, DISK_RATE, IPV4  # noqa: E402

PIP = ["cantera==3.2.0", "numpy==2.5.3", "scipy==1.18.1", "pytest==9.2.0"]


def cycle_cmd():
    return "\n".join([
        "set -e",
        "git rev-parse HEAD > \"$GTU_OUT/commit.txt\" 2>/dev/null || true",
        "{ python -V; python -c 'import cantera, numpy, scipy; print(cantera.__version__, numpy.__version__, scipy.__version__)'; "
        "echo \"nproc $(nproc) NPROC $NPROC itype ${GTU_ITYPE:-}\"; } > \"$GTU_OUT/env.txt\" 2>&1",
        "cd routes/gtu/validation",
        "mkdir -p /tmp/gtu_ref && cp balances_results.json validation_fleet_results.json /tmp/gtu_ref/",
        "python -m pytest -q -p no:cacheprovider > \"$GTU_OUT/pytest.txt\" 2>&1 || echo \"pytest rc=$?\" >> \"$GTU_OUT/pytest.txt\"",
        "T0=$(date +%s); python balances.py > \"$GTU_OUT/balances.log\" 2>&1; "
        "echo \"balances_s $(( $(date +%s) - T0 ))\" >> \"$GTU_OUT/timing.txt\"",
        "sed -i \"s/^workers = [0-9]*$/workers = $NPROC/\" validation_fleet.toml",
        "T0=$(date +%s); python validate_fleet.py > \"$GTU_OUT/validate_fleet.log\" 2>&1; "
        "echo \"validate_fleet_s $(( $(date +%s) - T0 )) workers $NPROC\" >> \"$GTU_OUT/timing.txt\"",
        "git checkout -- validation_fleet.toml 2>/dev/null || sed -i \"s/^workers = [0-9]*$/workers = 4/\" validation_fleet.toml",
        "cp balances_results.json balances.csv validation_fleet_results.json validation_fleet.csv \"$GTU_OUT/\"",
        "python compare_results.py balances_results.json /tmp/gtu_ref/balances_results.json --label balances "
        "--out \"$GTU_OUT/repro_balances.json\" || true",
        "python compare_results.py validation_fleet_results.json /tmp/gtu_ref/validation_fleet_results.json "
        "--label validate_fleet --out \"$GTU_OUT/repro_validation_fleet.json\" || true",
    ])


def fea_cmd():
    return "\n".join([
        "set -e", ENV_CCX,
        "git rev-parse HEAD > \"$GTU_OUT/commit.txt\" 2>/dev/null || true",
        "{ cat /opt/ccx_version.txt; python -V; echo \"nproc $(nproc) NPROC $NPROC itype ${GTU_ITYPE:-}\"; } "
        "> \"$GTU_OUT/env.txt\" 2>&1 || true",
        "python -m pytest -q -p no:cacheprovider routes/gtu/platform/fea/test_disk.py > \"$GTU_OUT/pytest_disk.txt\" 2>&1 "
        "|| echo \"pytest rc=$?\" >> \"$GTU_OUT/pytest_disk.txt\"",
        "T0=$(date +%s); python routes/gtu/platform/fea/disk.py --examples --out \"$GTU_OUT/disk_results.json\" "
        "> \"$GTU_OUT/disk.log\" 2>&1; echo \"disk_examples_s $(( $(date +%s) - T0 ))\" >> \"$GTU_OUT/timing.txt\"",
        "T0=$(date +%s); python routes/gtu/frame3/frame3_fea.py --work /tmp/gtu_frame3_fea --workers \"$NPROC\" "
        "--out \"$GTU_OUT/frame3_fea_results.json\" > \"$GTU_OUT/frame3_fea.log\" 2>&1; "
        "echo \"frame3_fea_s $(( $(date +%s) - T0 )) workers $NPROC\" >> \"$GTU_OUT/timing.txt\"",
        "python routes/gtu/validation/compare_results.py \"$GTU_OUT/disk_results.json\" routes/gtu/platform/fea/disk_results.json "
        "--label disk_examples --out \"$GTU_OUT/repro_disk.json\" || true",
        "python routes/gtu/validation/compare_results.py \"$GTU_OUT/frame3_fea_results.json\" routes/gtu/frame3/frame3_fea_results.json "
        "--label frame3_fea --out \"$GTU_OUT/repro_frame3_fea.json\" || true",
    ])


def cost(minutes, max_min, disk_gb):
    ex = disk_gb * DISK_RATE + IPV4
    return dict(machine_min_expected=minutes, usd_expected=round(minutes / 60 * (PRICE + ex), 3),
                usd_upper=round((max_min + 10) / 60 * (PRICE_UPPER + ex), 3))


def build():
    cyc_min = dict(boot_setup=2.0, pytest=1.0, balances=0.5, validate_fleet=2.0, compare_upload=1.0)
    fea_min = dict(boot_setup=3.5, test_disk=0.2, disk_examples=0.3, frame3_fea=2.0, compare_upload=1.0)
    est = {"val-cycle": cost(round(sum(cyc_min.values()), 1), 30, 30),
           "val-fea": cost(round(sum(fea_min.values()), 1), 30, 30)}
    return {
        "_comment": ("Проверка модели: повтор эталонных расчётов трека на машине AWS и сверка с веткой. Раннер берёт только "
                     "defaults и jobs; блоки estimate, merge, plan_check — для сессии с ключами. Код — ветка "
                     f"{BRANCH}. Сгенерировано make_jobs_validation.py; правки — в нём."),
        "defaults": {"instance_type": "c7a.16xlarge", "fallback_types": ["c7i.16xlarge"], "image": "al2023",
                     "disk_gb": 30, "max_min": 30},
        "jobs": [
            dict(name="val-cycle", python="3.12", pip=PIP, command=cycle_cmd()),
            dict(name="val-fea", container="ubuntu:24.04", container_setup=CCX_SETUP, command=fea_cmd()),
        ],
        "what": {
            "val-cycle": "балансы массы и энергии (8 точек) и проверка вне калибровки по 32 паспортам (validate_fleet.py, "
                         "~190 расчётов цикла) на Python 3.12 — повтор чисел VALIDATION_PLAN.md разделов 1–2 на другой машине",
            "val-fea": "CalculiX: 12 проверок против аналитики (вращающийся диск), примеры колёс и прочность KV-3 (73 расчёта) — "
                       "повтор с ccx 2.21 из чистого контейнера",
            "not_here": ("эталоны SU2 (Rotor 37, CC3, SRV2): нет геометрии и сеток в репозитории; NASA Rotor 37 — архив TMR не "
                         "скачан (навык su2-turbo, раздел 9). Замер скорости SU2 и CalculiX — routes/gtu/platform/bench/jobs_aws_bench.json"),
        },
        "merge": {
            "where": "сессия с ключами после run: результаты в routes/gtu/platform/aws/runs/<run_id>/val-cycle/ и val-fea/",
            "commands": ["python3 routes/gtu/validation/compare_results.py --collect routes/gtu/platform/aws/runs/<run_id> "
                         "--run-id <run_id> --out routes/gtu/validation/aws_repro_results.json"],
            "to_branch": ["routes/gtu/validation/aws_repro_results.json", "routes/gtu/platform/aws/ledger.jsonl"],
            "not_to_branch": ("пересчитанные *_results.json и CSV с машины — только если сверка показала расхождение и после "
                              "разбора (иначе они совпадают с веткой); журналы — в runs/, папка в .gitignore"),
            "checks": ("compare_results.py: числа |new − ref| ≤ 1e-9 + 1e-6·|ref|, служебные ключи (дата, версии, время, workers) "
                       "пропускаются; код 3 — расхождение: разбор до публикации (правило «качество прежде скорости»); "
                       "pytest.txt и pytest_disk.txt — тесты узлов на Python 3.12"),
        },
        "estimate": {
            "price": {"expected": PRICE, "upper": PRICE_UPPER,
                      "basis": "как routes/gtu/platform/bench/jobs_aws_bench.json → estimate.price"},
            "jobs": est,
            "total_usd_expected": round(sum(v["usd_expected"] for v in est.values()), 2),
            "total_usd_upper": round(sum(v["usd_upper"] for v in est.values()), 2),
            "upper_bound_rule": "раннер: (MaxPrice + диск + IPv4) × (max_min + 10) / 60 на машину",
            "minutes": {"val-cycle": cyc_min, "val-fea": fea_min},
            "sources": {
                "boot_setup": ("ledger.jsonl, прогон smoke2 20261114T231610Z-9b64: smoke-py (al2023 + uv + Python 3.12 + cantera) "
                               "— 1,04 мин машины, smoke-ccx (ubuntu:24.04 + apt calculix-ccx) — 2,05 мин (от запуска до гашения, "
                               "опрос раз в 60 с); сверху pip пакетов образа для val-fea ~1 мин — оценка"),
                "validate_fleet": ("validation_fleet_results.json → runtime_s 1139,5 с при workers 4 на 4 общих ядрах = ~76 ядро-мин; "
                                   "на 64 процессах ~1,2 мин + запуск пула — оценка 2 мин"),
                "balances_pytest": "VALIDATION_PLAN.md, раздел 6: balances.py ~20 с, pytest 25–50 с (4 ядра)",
                "fea": ("README узла fea: test_disk ~4 с, --examples ~7 с; frame3_fea.py --help: ~5 мин на 4 ядрах (73 расчёта) "
                        "= ~20 ядро-мин → ~0,3 мин на 64 + последовательные стержни лопаток — оценка 2 мин"),
            },
        },
    }


def plan_check(path):
    cmd = [sys.executable, os.path.join(ROOT, "routes/gtu/platform/aws/runner.py"), "plan", path, "--assume-price", "0.748",
           "--assume-zone-prices", "us-east-1a=0.955,us-east-1b=0.748,us-east-1c=0.706",
           "--out", os.environ.get("GTU_PLAN_OUT", "/tmp/gtu_plan_validation")]
    p = subprocess.run(cmd, capture_output=True, text=True, cwd=ROOT)
    lines = [ln for ln in p.stdout.splitlines() if re.match(r"\[(задача|бюджет|зоны|план)\]", ln)]
    return dict(command=" ".join(["python3", "routes/gtu/platform/aws/runner.py", "plan", REL] + cmd[4:8]),
                exit_code=p.returncode, lines=lines,
                runner_note="раннер v2 (python/uv, container); plan не вызывает AWS, цены — допущение")


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--plan-check", action="store_true")
    ap.add_argument("--mini-dir", default="")
    ap.add_argument("--out", default=OUT)
    a = ap.parse_args(argv)
    doc = build()
    for _ in range(2 if a.plan_check else 1):
        with open(a.out, "w", encoding="utf-8") as f:
            json.dump(doc, f, ensure_ascii=False, indent=1)
            f.write("\n")
        if a.plan_check and "plan_check" not in doc:
            doc["plan_check"] = plan_check(a.out)
            print("\n".join(doc["plan_check"]["lines"]))
            print("plan exit", doc["plan_check"]["exit_code"])
    if a.mini_dir:
        os.makedirs(a.mini_dir, exist_ok=True)
        open(os.path.join(a.mini_dir, "cmd_val-cycle.sh"), "w").write(cycle_cmd() + "\n")
        open(os.path.join(a.mini_dir, "cmd_val-fea.sh"), "w").write(fea_cmd() + "\n")
    e = doc["estimate"]
    print(f"записано {a.out}: задач {len(doc['jobs'])}, ожидаемо {e['total_usd_expected']} долл., верх {e['total_usd_upper']}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
