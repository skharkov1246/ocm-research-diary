"""Файлы задач AWS уровня 0 воронки (формат раннера ../platform/aws/runner.py v2) со сметой и проверкой командой plan.

  python make_jobs_l0.py                      # jobs_l0_debug.json (KV-3, N=2000, 10 частей) и jobs_l0.json (полный)
  python make_jobs_l0.py --shards 40 --machines 4

Части (--shards, вход прогона; 10 у отладки, 20-50 у полного) - единственное деление выборки, по индексу варианта
(i % S == k). Части группируются в задачи-машины (--machines): машина считает свои части одним пулом процессов на все
ядра ($NPROC), по файлу на часть. Счёт части - секунды (смета ниже), поэтому машина на часть потратила бы почти всё время
на подготовку; число машин <= 6 (квота спота 384 vCPU = 6 × c7a.16xlarge) - одна волна без очереди.
Слияние и отбор по нишам - блок "merge" файла (раннер его не читает: он берёт только defaults и jobs) - выполняются в
сессии с ключами после run (раннер сам скачивает результаты в routes/gtu/platform/aws/runs/<run_id>/<задача>/).
Смета: машино-часы × тариф спота + диск + IPv4 + хранение S3 + передача; время счёта - по замеру короткой пробы
(funnel_l0_results.json, timing) с запасом SAFETY; подготовка машины - оценка SETUP_MIN. AWS не вызывается.
"""
import argparse, json, math, os, re, subprocess, sys

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))
RUNNER = os.path.join(REPO, "routes", "gtu", "platform", "aws", "runner.py")
REL = "routes/gtu/funnel"

PIP = ["cantera==3.2.0", "scipy==1.18.1", "numpy==2.5.3"]      # закреплены: скремблирование Соболя зависит от версии scipy
PROFILES = {
    "debug": dict(file="jobs_l0_debug.json", runs=[("KV-3", 2000)], shards=10, machines=2, max_min=25,
                  merge=dict(top=40, target=150, k_near=1, near_max=0.25),
                  why="отладочный прогон: путь части -> машины -> S3 -> слияние на настоящем раннере; отбор около 150"),
    "full": dict(file="jobs_l0.json", runs=[("KV-3", 200000), ("KV-1", 50000), ("KV-6", 50000)], shards=30, machines=3, max_min=30,
                 merge=dict(top=1000, target=4000, k_near=1, near_max=0.25),
                 why="полный уровень 0 (приоритет владельца 14.11.2026: KV-3 200 000, KV-1 и KV-6 по 50 000; отбор 3000-5000)"),
}
PRICE = dict(expected=0.706, alt=0.748, upper=1.4206,
             basis="ожидаемая: спот c7a.16xlarge us-east-1c 0,706 долл./ч, us-east-1b 0,748 (замер сессии с ключами 22:08 14.11.2026, "
                   "задача владельца); верхняя граница: MaxPrice раннера 1,4206 долл./ч (максимум по зонам × 1,2, ledger.jsonl "
                   "живого прогона smoke 22:14 14.11.2026)")
SETUP_MIN = 5.0     # оценка: smoke 14.11.2026 - 61 с от запуска до конца без пакетов; плюс uv, Python 3.12, пакеты, клон
UPLOAD_MIN = 1.0    # оценка
SAFETY = 3.0        # запас на скорость и загрузку ядер машины против короткой пробы (4 общих ядра)
VCPU = 64
EBS_GB, EBS_USD_GB_MONTH, IPV4_USD_H = 30, 0.08, 0.005           # ставки раннера (допущения кода, не сверены)
S3_USD_GB_MONTH, EGRESS_USD_GB = 0.023, 0.09                      # допущения (сайты AWS из контейнера не открываются)


def measured(results=os.path.join(HERE, "funnel_l0_results.json")):
    """время на концепцию (с, одно ядро) и байт на концепцию в файле части - из короткой пробы"""
    d = json.load(open(results, encoding="utf-8"))
    t = d["probe"]["timing"]["s_per_concept_serial"]
    b = d.get("probe", {}).get("part_bytes_per_concept")
    return t, b, d["probe"]["params"]


def build(profile, shards=None, machines=None, branch="claude/lucid-allen-yvps5d"):
    P = dict(PROFILES[profile])
    S = shards or P["shards"]
    M = machines or P["machines"]
    lo, hi = (10, 50) if profile == "debug" else (20, 50)
    if not lo <= S <= hi:
        raise ValueError(f"{profile}: частей {S}, нужно {lo}-{hi}")
    if not 1 <= M <= min(6, S):
        raise ValueError(f"машин {M}: от 1 до min(6, частей)")
    jobs = []
    for j in range(M):
        mine = [k for k in range(S) if k % M == j]
        lines = [f"python {REL}/funnel_l0.py --kv {kv} --n {n} --shards {S} --shard {','.join(map(str, mine))} "
                 f"--workers \"$NPROC\" --gzip --out \"$GTU_OUT\"" for kv, n in P["runs"]]
        lines.append("ls -l \"$GTU_OUT\" > \"$GTU_OUT/_parts_list.txt\"")
        jobs.append(dict(name=f"l0-{profile}-m{j}", command="\n".join(lines), max_min=P["max_min"]))
    merge = []
    for kv, n in P["runs"]:
        m = P["merge"]
        merge.append(f"python3 {REL}/merge_l0.py --kv {kv} --n {n} --shards {S} --in routes/gtu/platform/aws/runs/<run_id> "
                     f"--out-dir {REL}/l0_aws_{profile} --suffix _{kv}_n{n} --top {m['top']} --target {m['target']} "
                     f"--k-near {m['k_near']} --near-max {m['near_max']}" + (" --csv-gz" if n > 20000 else ""))
    return P, S, M, jobs, merge


def estimate(P, S, M, jobs, t_concept, bytes_per_concept):
    n_tot = sum(n for _, n in P["runs"])
    core_s = n_tot * t_concept * SAFETY
    per_machine_compute_min = core_s / M / VCPU / 60
    machine_min = SETUP_MIN + per_machine_compute_min + UPLOAD_MIN
    mh = M * machine_min / 60
    ebs = EBS_GB * EBS_USD_GB_MONTH / 730
    gb = n_tot * (bytes_per_concept or 1500) / 1e9
    gb_csv = n_tot * 700 / 1e9                                     # оценка: строка CSV слияния
    usd = dict(spot_expected=round(mh * PRICE["expected"], 3), spot_alt=round(mh * PRICE["alt"], 3),
               disk_ipv4=round(mh * (ebs + IPV4_USD_H), 4), s3_month=round((gb + gb_csv) * S3_USD_GB_MONTH, 4),
               transfer=round(gb * EGRESS_USD_GB, 4))
    usd["total_expected"] = round(usd["spot_alt"] + usd["disk_ipv4"] + usd["s3_month"] + usd["transfer"], 3)
    bound = M * (PRICE["upper"] + ebs + IPV4_USD_H) * (P["max_min"] + 10) / 60
    return dict(concepts=n_tot, s_per_concept_core_measured=t_concept, safety=SAFETY, core_hours=round(core_s / 3600, 3),
                machines=M, shards=S, vcpu_per_machine=VCPU, compute_min_per_machine=round(per_machine_compute_min, 2),
                setup_min=SETUP_MIN, upload_min=UPLOAD_MIN, machine_min=round(machine_min, 2), machine_hours=round(mh, 3),
                parts_gb_gz=round(gb, 4), usd=usd, upper_bound_usd=round(bound, 2),
                upper_bound_rule="раннер: (MaxPrice + диск + IPv4) × (max_min + 10) / 60 на машину",
                price=PRICE, assumptions=dict(setup="оценка: smoke 14.11.2026 - 61 с без пакетов; uv, Python 3.12, пакеты, клон - сверху",
                                              s3=f"{S3_USD_GB_MONTH} долл./ГБ-мес, передача {EGRESS_USD_GB} долл./ГБ - допущения",
                                              disk_ipv4="ставки раннера, допущения кода"),
                cap_300_usd_ok=usd["total_expected"] < 300 and bound < 300)


def plan_check(path, outdir):
    cmd = [sys.executable, RUNNER, "plan", path, "--assume-price", str(PRICE["alt"]),
           "--assume-zone-prices", f"us-east-1a=0.955,us-east-1b={PRICE['alt']},us-east-1c={PRICE['expected']}", "--out", outdir]
    r = subprocess.run(cmd, capture_output=True, text=True, cwd=REPO)
    out = r.stdout + r.stderr
    keep = [ln for ln in out.splitlines() if re.match(r"^\[(задача|план|зоны|очередь|бюджет)\]", ln.strip())]
    return dict(command="python3 routes/gtu/platform/aws/runner.py plan " + os.path.relpath(path, REPO) +
                        f" --assume-price {PRICE['alt']} --assume-zone-prices us-east-1a=0.955,us-east-1b={PRICE['alt']},"
                        f"us-east-1c={PRICE['expected']}",
                exit_code=r.returncode, lines=keep[:40],
                runner_note="раннер v2 в рабочей копии (поля python, pip, NPROC, самая дешёвая зона); если v2 ещё не влит в ветку - "
                            "перепроверить plan после влития")


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--shards", type=int, default=None, help="частей полного прогона (20-50)")
    ap.add_argument("--machines", type=int, default=None, help="машин полного прогона (<= 6)")
    ap.add_argument("--plan-out", default=os.environ.get("TMPDIR", "/tmp"))
    a = ap.parse_args(argv)
    t, b, probe = measured()
    for prof in ("debug", "full"):
        P, S, M, jobs, merge = build(prof, a.shards if prof == "full" else None, a.machines if prof == "full" else None)
        est = estimate(P, S, M, jobs, t, b)
        doc = dict(_comment=("Уровень 0 воронки (README узла funnel, раздел «Уровень 0»). Раннер берёт только defaults и jobs; блоки "
                             "merge, estimate, plan_check - для сессии с ключами. Код - ветка claude/lucid-allen-yvps5d (PR #62). "
                             "Сгенерировано make_jobs_l0.py; правки - в нём."),
                   defaults=dict(instance_type="c7a.16xlarge", fallback_types=["c7i.16xlarge"], image="al2023", disk_gb=EBS_GB,
                                 python="3.12", pip=PIP),
                   jobs=jobs,
                   merge=dict(where="сессия с ключами после run: результаты уже в routes/gtu/platform/aws/runs/<run_id>/<задача>/ "
                                    "(или collect --run-id <run_id> --job <задача> по каждой); Python с теми же версиями пакетов "
                                    "(pip install " + " ".join(PIP) + ") - иначе отпечаток выборки не совпадёт",
                              commands=merge, profile_why=P["why"],
                              checks="слияние откажет (код 2), если потеряна или задвоена хоть одна часть или индекс, или у частей "
                                     "разные отпечатки выборки и кода"),
                   estimate=est)
        path = os.path.join(HERE, P["file"])
        with open(path, "w", encoding="utf-8") as fh:
            json.dump(doc, fh, ensure_ascii=False, indent=1)
        doc["plan_check"] = plan_check(path, a.plan_out)
        with open(path, "w", encoding="utf-8") as fh:
            json.dump(doc, fh, ensure_ascii=False, indent=1)
        print(f"{P['file']}: {len(jobs)} машин, {S} частей, {est['concepts']} концепций, ожидаемо {est['usd']['total_expected']} долл., "
              f"верхняя граница {est['upper_bound_usd']} долл.; plan: код {doc['plan_check']['exit_code']}")
        for ln in doc["plan_check"]["lines"]:
            if ln.startswith("[план]"):
                print("  ", ln)
    return 0


if __name__ == "__main__":
    sys.exit(main())
