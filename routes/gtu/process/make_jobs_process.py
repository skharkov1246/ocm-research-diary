#!/usr/bin/env python3
"""Задачи AWS перебора сквозных концепций «топливо -> вал» (MAP-Elites, process_elites.py) -> jobs_process.json:
острова по машинам, слияние с контролем полноты, смета, проверка командой раннера plan (код 0).

Профиль (решения владельца 14-15.11.2026: всё тяжёлое - на AWS, все CPU-мощности во всех регионах, потолок 2000 долл.;
KV-3 - не меньше 1-2 млн концепций, KV-1 и KV-6 - сотни тысяч):
  * машин M = 12 - одна волна под фактические квоты спота (us-east-1 384 vCPU и us-west-2 384 vCPU = 12 × c7a.16xlarge;
    us-east-2 - 5 vCPU, заявки до 1536 на рассмотрении - README раннера, проверка сессии с ключами 15.11.2026);
  * каждая машина - по одному острову каждой точки (часть j = остров j): KV-3, затем KV-1, затем KV-6; остров считает
    пакеты потомков пулом на все ядра ($NPROC); число островов и оценок острова - вход прогона (флаги ниже);
  * части - единственное деление: остров i -> часть i mod M (evo/core/islands.islands_of_part).
Слияние и сводка - блок "merge" файла (раннер его не читает) - в сессии с ключами после run: результаты раннер кладёт в
routes/gtu/platform/aws/runs/<run_id>/<задача>/.

Смета: время ЦП на концепцию - мини-прогон KV-3 и пробы KV-1, KV-6 (process_results.json -> runs.*.cpu_s_per_eval.all;
время ЦП процессов пула на общем контейнере под нагрузкой соседних задач) × запас SAFETY; КПД пула (барьер пакета) EFF_POOL;
подготовка машины SETUP_MIN; цены спота - замер сессии с ключами 15.11.2026 (README раннера). AWS не вызывается.

  python3 make_jobs_process.py                          # jobs_process.json + plan (код 0)
  python3 make_jobs_process.py --machines 12 --evals-kv3 400000 --evals-kv1 100000 --evals-kv6 100000
"""
import argparse
import json
import math
import os
import re
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
RUNNER = os.path.join(REPO, "routes", "gtu", "platform", "aws", "runner.py")
REL = "routes/gtu/process"
RESULTS = os.path.join(HERE, "process_results.json")
OUT = os.path.join(HERE, "jobs_process.json")
BRANCH = "claude/lucid-allen-yvps5d"        # ветка, куда потоки трека сводят черновики; код process/ должен быть в ней до run
PIP = ["cantera==3.2.0", "scipy==1.18.1", "numpy==2.5.3"]      # как задачи воронки: скремблирование Соболя зависит от версии scipy
VCPU = 64
SETUP_MIN = 6.0         # подготовка машины: smoke 14.11.2026 - 61 с без пакетов; + uv, Python 3.12, пакеты, клон (оценка, как у воронки 5) + импорт
UPLOAD_MIN = 1.0        # оценка
SAFETY = 1.5            # запас к времени ЦП мини-прогона (оно уже завышено: общий контейнер, нагрузка 6-10 на 4 ядрах) - оценка
EFF_POOL = 0.80         # доля полезной загрузки ядер пулом с барьером пакета (разброс времени концепции 0,2-100 мс) - оценка
BATCH = 1024            # потомков в пакете острова на машине (16 на ядро)
REGIONS = "us-east-1,us-east-2,us-west-2"
QUOTAS = "us-east-1=384,us-east-2=5,us-west-2=384"
ZONES = "us-east-1c=0.703,us-east-1b=0.747,us-east-1d=0.879,us-east-2b=0.601,us-east-2a=0.606,us-west-2a=0.866,us-west-2d=0.952"
PRICE = dict(expected=0.80, upper=0.952 * 1.2,
             basis="спот c7a.16xlarge (замер сессии с ключами 15.11.2026 10:45 UTC, README раннера): us-east-1c 0,703, 1b 0,747, 1d 0,879; "
                   "us-west-2a 0,866, 2d 0,952; us-east-2b 0,601 (квота 5 vCPU); ожидаемая - 0,80 (6 машин us-east-1 по 0,70-0,75 и 6 - "
                   "us-west-2 по 0,87-0,95); верхняя - MaxPrice раннера: максимум по зонам × 1,2")
EBS_GB, EBS_USD_GB_MONTH, IPV4_USD_H = 30, 0.08, 0.005           # ставки раннера (допущения кода)
S3_USD_GB_MONTH, XFER_USD_GB = 0.023, 0.02                        # допущения; межрегиональная передача - справочник раннера (предварительно)
PROFILE = dict(machines=12, evals={"KV-3": 400000, "KV-1": 100000, "KV-6": 100000}, init_frac=0.3)


def measured(path=RESULTS):
    """время ЦП на концепцию (с) по точкам и байт части на остров - из мини-прогона и проб"""
    d = json.load(open(path, encoding="utf-8"))
    t = {kv: r["cpu_s_per_eval"]["all"] for kv, r in d["runs"].items()}
    fs = {kv: r["feasible_share"] for kv, r in d["runs"].items()}
    b = {kv: r.get("part_bytes_per_island") for kv, r in d["runs"].items()}
    return t, fs, b, d


def build(machines, evals, init_frac):
    for kv, n in evals.items():
        if n < 1000:
            raise ValueError(f"{kv}: оценок острова {n} - мало для MAP-Elites")
    if not 10 <= machines <= 50:
        raise ValueError(f"частей (машин) {machines}: нужно 10-50 (каждая точка делится на столько частей)")
    t, fs, bpi, _ = measured()
    comp_min = sum(evals[kv] * t[kv] * SAFETY / (VCPU * EFF_POOL) / 60 for kv in evals)
    max_min = max(30, int(math.ceil(SETUP_MIN + 2 * comp_min + 10)))
    jobs = []
    for j in range(machines):
        lines = ["rc=0"]
        for kv, n in evals.items():
            lines.append(f"python {REL}/process_elites.py part --kv {kv} --n-islands {machines} --n-parts {machines} --part {j} "
                         f"--evals {n} --init-frac {init_frac} --batch {BATCH} --workers \"$NPROC\" --out \"$GTU_OUT\" || rc=$?")
        lines += ["ls -l \"$GTU_OUT\" > \"$GTU_OUT/_parts_list.txt\"", "exit $rc"]
        jobs.append(dict(name=f"process-m{j:02d}", command="\n".join(lines), max_min=max_min))
    merge = []
    for kv, n in evals.items():
        merge.append(f"python3 {REL}/process_elites.py merge --kv {kv} --n-islands {machines} --evals {n} --init-frac {init_frac} "
                     f"--batch {BATCH} --in routes/gtu/platform/aws/runs/<run_id> --out-json {REL}/process_aws_{kv}.json "
                     f"--archive-out routes/gtu/platform/aws/runs/<run_id>/merged/process_{kv}_archive.json.gz")
    return dict(jobs=jobs, merge=merge, max_min=max_min, comp_min=comp_min, t=t, fs=fs, bpi=bpi)


def estimate(machines, evals, B):
    n_tot = sum(evals.values()) * machines
    core_s = sum(evals[kv] * machines * B["t"][kv] * SAFETY for kv in evals)
    per_machine_min = B["comp_min"]
    machine_min = SETUP_MIN + per_machine_min + UPLOAD_MIN
    mh = machines * machine_min / 60
    ebs = EBS_GB * EBS_USD_GB_MONTH / 730
    part_mb = {kv: (B["bpi"].get(kv) or 5e6) / 1e6 for kv in evals}
    gb = machines * sum(part_mb.values()) / 1e3
    usd = dict(spot=round(mh * PRICE["expected"], 3), disk_ipv4=round(mh * (ebs + IPV4_USD_H), 4), s3_month=round(gb * S3_USD_GB_MONTH, 4),
               xfer=round(gb * XFER_USD_GB, 4))
    usd["total_expected"] = round(sum(usd.values()), 3)
    bound = machines * (PRICE["upper"] + ebs + IPV4_USD_H) * (B["max_min"] + 10) / 60
    return dict(concepts_total=n_tot, concepts_by_kv={kv: n * machines for kv, n in evals.items()},
                cpu_s_per_concept_measured=B["t"], feasible_share_measured=B["fs"], safety=SAFETY, pool_efficiency=EFF_POOL,
                core_hours=round(core_s / 3600, 2), machines=machines, vcpu_per_machine=VCPU, compute_min_per_machine=round(per_machine_min, 1),
                setup_min=SETUP_MIN, upload_min=UPLOAD_MIN, machine_min=round(machine_min, 1), machine_hours=round(mh, 2),
                watchdog_max_min=B["max_min"], parts_gb_gz=round(gb, 3), usd=usd, upper_bound_usd=round(bound, 2),
                upper_bound_rule="раннер: (MaxPrice + диск + IPv4) × (max_min + 10) / 60 на машину",
                price=PRICE, cap_2000_usd_ok=bound < 2000,
                note=("ожидаемо единицы долларов, а не десятки: концепция считается за десятки миллисекунд (быстрые модели уровня облика); "
                      "при той же цене N можно поднять в разы (флаги --evals-*), верхняя граница растёт линейно со сторожем"))


def plan_check(path, outdir, extra=()):
    cmd = [sys.executable, RUNNER, "plan", path, "--branch", BRANCH, "--assume-price", str(PRICE["expected"]), "--regions", REGIONS,
           "--assume-zone-prices", ZONES, "--vcpu-quota-per-region", QUOTAS, "--out", outdir, *extra]
    r = subprocess.run(cmd, capture_output=True, text=True, cwd=REPO)
    out = r.stdout + r.stderr
    keep = [ln for ln in out.splitlines() if re.match(r"^\[(задача|план|зоны|очередь|бюджет|регион|квоты|квота)\]", ln.strip())]
    shown = " ".join(x if " " not in x else f"'{x}'" for x in cmd[2:]).replace(path, os.path.relpath(path, REPO)).replace(outdir, "<папка>")
    return dict(command="python3 routes/gtu/platform/aws/runner.py " + shown, exit_code=r.returncode,
                lines=keep[:40] or out.splitlines()[-15:])


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--machines", type=int, default=PROFILE["machines"])
    ap.add_argument("--evals-kv3", type=int, default=PROFILE["evals"]["KV-3"])
    ap.add_argument("--evals-kv1", type=int, default=PROFILE["evals"]["KV-1"])
    ap.add_argument("--evals-kv6", type=int, default=PROFILE["evals"]["KV-6"])
    ap.add_argument("--init-frac", type=float, default=PROFILE["init_frac"])
    ap.add_argument("--plan-out", default=os.environ.get("GTU_SCRATCH", os.environ.get("TMPDIR", "/tmp")))
    ap.add_argument("--out", default=OUT)
    a = ap.parse_args(argv)
    evals = {"KV-3": a.evals_kv3, "KV-1": a.evals_kv1, "KV-6": a.evals_kv6}
    B = build(a.machines, evals, a.init_frac)
    est = estimate(a.machines, evals, B)
    doc = dict(
        _comment=("Перебор сквозных концепций «топливо -> вал» (routes/gtu/process, README): MAP-Elites, ниши «топология × камера × "
                  "архитектура × изготовление», цель - КПД на тонну, фронты КПД - масса в нишах. Острова по машинам: машина j считает "
                  f"остров j точек {', '.join(evals)} (часть j из {a.machines}). Раннер берёт только defaults и jobs; блоки merge, "
                  f"estimate, plan_check - для сессии с ключами. Код - ветка {BRANCH}. Сгенерировано make_jobs_process.py; правки - в нём."),
        defaults=dict(instance_type="c7a.16xlarge", fallback_types=["c7i.16xlarge"], image="al2023", disk_gb=EBS_GB, python="3.12", pip=PIP),
        jobs=B["jobs"],
        merge=dict(where="сессия с ключами после run: результаты - routes/gtu/platform/aws/runs/<run_id>/<задача>/ (или collect по каждой); "
                         "Python с теми же версиями пакетов (pip install " + " ".join(PIP) + ") - иначе отпечаток конфигурации не совпадёт",
                   commands=B["merge"],
                   checks=("слияние откажет (код 2), если нет хоть одного острова, остров пришёл дважды, у части другая конфигурация "
                           "(отпечаток: версия, гены, окна, сетка, сиды, опорные родители) или остров недосчитан"),
                   into_branch=f"в ветку - только сводки {REL}/process_aws_KV-*.json; слитые архивы (*.json.gz) - в runs/ (вне ветки) и S3",
                   next=("эпоха 2 (по желанию): слитый архив эпохи 1 - начальные родители (--parents <архив>), тот же файл задач с "
                         "новым сидом (--seed); острова внутри эпохи не обмениваются (итог не зависит от деления на части)")),
        estimate=est,
        run=dict(dryrun=f"python3 routes/gtu/platform/aws/runner.py run {REL}/jobs_process.json --branch {BRANCH} --regions {REGIONS}",
                 live=f"python3 routes/gtu/platform/aws/runner.py run {REL}/jobs_process.json --branch {BRANCH} --regions {REGIONS} --live",
                 on_demand=("+ --on-demand-fallback: если спот не даст машину ни в одном регионе - машины по требованию (us-east-1 384, "
                            "us-west-2 300 vCPU; цена справочника, дороже спота) - по решению владельца"),
                 note="раннер v3: квоты регионов - фактические (service-quotas или флаг --vcpu-quota-per-region); us-east-2 при квоте 5 vCPU "
                      "машину не получит, после одобрения заявки - получит без правки файла"),
    )
    os.makedirs(a.plan_out, exist_ok=True)
    with open(a.out, "w", encoding="utf-8") as f:
        json.dump(doc, f, ensure_ascii=False, indent=1)
    doc["plan_check"] = plan_check(a.out, a.plan_out)
    doc["plan_check_on_demand"] = plan_check(a.out, a.plan_out, ("--on-demand-fallback", "--assume-no-spot",
                                                                 "--od-vcpu-quota-per-region", "us-east-1=384,us-east-2=5,us-west-2=300"))
    with open(a.out, "w", encoding="utf-8") as f:
        json.dump(doc, f, ensure_ascii=False, indent=1)
    print(f"{os.path.relpath(a.out, REPO)}: {len(B['jobs'])} машин, концепций {est['concepts_total']:,} ".replace(",", " ") +
          f"({', '.join(f'{k} {v:,}'.replace(',', ' ') for k, v in est['concepts_by_kv'].items())}); ожидаемо {est['usd']['total_expected']} долл., "
          f"верхняя граница {est['upper_bound_usd']} долл.; сторож {B['max_min']} мин; plan: код {doc['plan_check']['exit_code']}, "
          f"по требованию: код {doc['plan_check_on_demand']['exit_code']}")
    for ln in doc["plan_check"]["lines"]:
        if ln.startswith("[план]") or ln.startswith("[бюджет]"):
            print("  ", ln)
    return 0 if doc["plan_check"]["exit_code"] == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
