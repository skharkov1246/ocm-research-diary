"""Пакет задач уровня 0 v2 камеры KV для раннера AWS v2 (routes/gtu/platform/aws/runner.py): объём по пользе, смета,
холостой plan, мини-проверка команд.

Объём выборки - «до сходимости доли «держит 30-100 %»» по пробе N = 300 (comb_l0_v2_probe_results.json), правило
записано до прогона (оценка):
  К1 - если в пробе нет ни одного варианта с новой целью: N >= 3 / P_MIN (правило трёх: при нуле в выборке N верхняя
       граница доли с доверием 95 % - 3/N); P_MIN = 1e-5 - доля, ниже которой схема на уровне облика не стоит уровня 1;
  К2 - если есть h вариантов из n: N >= 100 / (h/n) (не меньше 100 вариантов - относительная ошибка доли <= 10 %);
  К3 - доля «держит 30-100 % устойчиво и полно» в ячейках «схема × сочетание пассивных элементов» без RQL: N >= 25 /
       (p_s · f_ячейки) (25 событий в ячейке - относительная ошибка <= 20 %); p_s - по пробе (сумма по ячейкам без RQL),
       f_ячейки - доля самой узкой ячейки в выборке;
  N = 2^ceil(log2(max К)) (выборка Соболя сбалансирована на степенях двойки), не больше 2^20 (слияние и выгрузка -
  проверено на 2^18 для v1). Сходимость проверяет слияние: доля по вложенным префиксам 2^k (merge_l0_v2.convergence).
Чувствительность к закону затухания при малом J (BETA_LOWJ, оценка) - две серии β = 0 и β = 2 объёмом N/4 каждая.

Раскладка: 48 частей по индексу (i % 48), 6 задач по 8 частей (квота спота 384 vCPU = 6 машин c7a.16xlarge);
серии β - по 4 задачи × 12 частей. Время варианта - процессорное время пробы (cpu_s в частях) × поправка c7a,
измеренная v1 (на AWS 0,617 с на вариант против 0,553 с процессорного времени в сессии на том же коде - jobs_aws_comb_l0.json).

Запуск: python make_jobs_l0_v2.py [--probe comb_l0_v2_probe_results.json] [--mini]   (AWS не вызывается)
"""
import argparse, json, math, os, re, shutil, subprocess, sys, tempfile, time

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
RUNNER = os.path.join(REPO, "routes", "gtu", "platform", "aws", "runner.py")
OUT = os.path.join(HERE, "jobs_aws_comb_l0_v2.json")
PROBE = os.path.join(HERE, "comb_l0_v2_probe_results.json")
V1_JOBS = os.path.join(HERE, "jobs_aws_comb_l0.json")
V1_AWS = os.path.join(HERE, "l0_aws", "comb_l0_KV-3_n262144_results.json")

P_MIN = 1e-5
N_CAP = 2 ** 20
SHARDS, JOBS = 48, 6
BETA_SHARDS, BETA_JOBS = 48, 4
VCPU = 64
SETUP_MIN, UPLOAD_MIN, SAFETY = 3.0, 2.0, 2.0
ZONES = "us-east-1a=0.9554,us-east-1b=0.7478,us-east-1c=0.7064,us-east-1d=0.8576,us-east-1f=1.1804"
PRICE = 0.75


def volume(probe):
    """объём по правилам К1-К3 (см. описание модуля)"""
    c = probe["counts"]
    n = c["n"]
    h = c["lowem_30_100_1valve"]
    by = probe["by_scheme_passive"]
    crit = {}
    crit["K1_zero_bound"] = dict(n_min=math.ceil(3 / P_MIN), applies=h == 0, rule="3 / P_MIN, P_MIN = 1e-5")
    crit["K2_holders"] = dict(n_min=math.ceil(100 * n / h) if h else None, applies=h > 0, rule="100 / доля")
    non_rql = {k: v for k, v in by.items() if not k.startswith("RQL|")}
    tot = sum(v["total"] for v in non_rql.values())
    st = sum(v["stable_30_100"] for v in non_rql.values())
    import comb_v2 as V
    # доля самой узкой ячейки «схема × сочетание»: схема 1/5, сочетаний у LPP/MM/FLOX 3·2·2 = 12
    f_cell = min((1 / len(V.SCHEME_NAMES)) / (len(V.ALLOWED["split"][s]) * len(V.ALLOWED["strat"][s]) *
                                             len(V.ALLOWED["recirc"][s])) for s in V.SCHEME_NAMES)
    p_s = st / tot if tot else 0.0
    p_s_used = p_s if st else 3.0 / max(tot, 1)              # нет событий - верхняя граница 3/n (консервативно - больше N)
    crit["K3_stable_per_cell"] = dict(n_min=math.ceil(25 / (p_s_used * f_cell)), applies=True,
                                      p_stable_non_rql=p_s, events=st, of=tot, f_cell=f_cell,
                                      rule="25 / (p_s · f_ячейки)")
    need = max(v["n_min"] for v in crit.values() if v["applies"] and v["n_min"])
    N = min(N_CAP, 2 ** math.ceil(math.log2(need)))
    return N, crit, need


def cpu_per_variant(probe):
    """процессорное время варианта v2 по пробе и поправка c7a по v1"""
    v1j = json.load(open(V1_JOBS, encoding="utf-8"))
    v1a = json.load(open(V1_AWS, encoding="utf-8"))
    s_v1_session = max(m["value_s"] for m in v1j["estimate"]["measurements"]["cpu_per_variant"][:2])
    s_v1_aws = v1a["timing"]["t_variant_mean_s"]
    k = s_v1_aws / s_v1_session
    s_v2 = probe["timing"]["cpu_variant_mean_s"]
    return dict(cpu_v2_session_s=s_v2, k_c7a=k, s_v2_aws=s_v2 * k,
                basis=dict(v1_session_cpu_s=s_v1_session, v1_aws_wall_per_variant_s=s_v1_aws,
                           note="время варианта на машине AWS v1 (стена процесса в пуле 64) / процессорное время в сессии"))


def job_cmd(n, shards, parts, beta=None):
    b = "" if beta is None else f" --beta {beta:g}"
    tag = "" if beta is None else f"_b{beta:g}"
    lst = " ".join(str(s) for s in parts)
    return (f"set -e\nfor S in {lst}; do\n"
            f"  python -u routes/gtu/combustor/comb_l0_v2.py --kv KV-3 --n {n} --shards {shards} --shard \"$S\""
            f"{b} --workers \"$NPROC\" --out \"$GTU_OUT\"\n"
            f"  gzip -6 \"$GTU_OUT\"/comb_l0v2_KV-3_n{n}{tag}_s*.json\ndone\n"
            f"ls -l \"$GTU_OUT\" > \"$GTU_OUT/_parts_list.txt\"")


def series(name, n, shards, jobs, s_var, beta=None):
    per_job = [list(range(j, shards, jobs)) for j in range(jobs)]
    var_job = n / jobs
    comp = var_job * s_var / VCPU / 60.0
    mm = int(5 * math.ceil((SETUP_MIN + SAFETY * comp + UPLOAD_MIN + 8) / 5))
    mm = max(20, mm)
    js = [dict(name=f"{name}-m{j}", command=job_cmd(n, shards, p, beta), max_min=mm) for j, p in enumerate(per_job)]
    est = dict(n=n, shards=shards, jobs=jobs, beta=beta, variants_per_job=var_job, compute_min_per_job=comp,
               max_min=mm, core_hours=n * s_var / 3600.0, machine_hours_expected=jobs * (SETUP_MIN + comp + UPLOAD_MIN) / 60)
    return js, est


def plan(path, price=PRICE, zones=ZONES):
    r = subprocess.run([sys.executable, RUNNER, "plan", path, "--assume-price", str(price), "--assume-zone-prices", zones],
                       capture_output=True, text=True, cwd=REPO, timeout=600)
    out = r.stdout + r.stderr
    m = re.search(r"прогноз пакета ([0-9.]+) долл", out)
    return dict(command=f"python3 routes/gtu/platform/aws/runner.py plan routes/gtu/combustor/{os.path.basename(path)} "
                        f"--assume-price {price} --assume-zone-prices {zones}",
                exit_code=r.returncode, forecast_usd=float(m.group(1)) if m else None,
                lines=[ln for ln in out.splitlines() if ln.startswith(("[задача]", "[план]", "[очередь]", "[бюджет]"))][:40])


def mini_check(pkg, n_mini=96):
    """команды задач с n = n_mini в копии ветки (git archive HEAD + файлы v2 рабочей папки), NPROC = 4; слияние"""
    tmp = tempfile.mkdtemp(prefix="comb_v2_mini_")
    subprocess.run(f"git archive HEAD routes/gtu | tar -x -C {tmp}", shell=True, check=True, cwd=REPO)
    for f in ("comb_v2.py", "comb_l0_v2.py", "merge_l0_v2.py", "comb_bc_v2_results.json", "targets_v2.json"):
        shutil.copy(os.path.join(HERE, f), os.path.join(tmp, "routes", "gtu", "combustor", f))
    res = dict(n=n_mini, how="копия: git archive HEAD routes/gtu + файлы v2 рабочей папки (ещё не в ветке); команды задач "
                             "пакета с n = 96, из корня копии, GTU_OUT - своя папка задачи, NPROC = 4", jobs=[])
    outs = []
    for j in pkg["jobs"]:
        if not j["name"].startswith("comb-l0v2-m"):
            continue
        cmd = re.sub(r"--n \d+", f"--n {n_mini}", j["command"]).replace("_n" + str(pkg["estimate"]["main"]["n"]),
                                                                       f"_n{n_mini}")
        out = os.path.join(tmp, "out", j["name"])
        os.makedirs(out, exist_ok=True)
        t0 = time.time()
        r = subprocess.run(["bash", "-c", cmd], cwd=tmp, env=dict(os.environ, GTU_OUT=out, NPROC="4"),
                           capture_output=True, text=True)
        files = [f for f in os.listdir(out) if f.endswith(".json.gz")]
        res["jobs"].append(dict(job=j["name"], exit_code=r.returncode, wall_s=round(time.time() - t0, 1), files=len(files)))
        outs.append(out)
    parts = os.path.join(tmp, "parts")
    os.makedirs(parts, exist_ok=True)
    k = 0
    for o in outs:
        for f in os.listdir(o):
            if f.endswith(".json.gz"):
                subprocess.run(f"gunzip -c {os.path.join(o, f)} > {os.path.join(parts, f[:-3])}", shell=True, check=True)
                k += 1
    r = subprocess.run([sys.executable, os.path.join(tmp, "routes", "gtu", "combustor", "merge_l0_v2.py"), parts,
                        "--out", os.path.join(tmp, "merged.json"), "--top", "5"], capture_output=True, text=True, cwd=tmp)
    res["merge"] = dict(parts=k, exit_code=r.returncode, tail=(r.stdout + r.stderr).strip().splitlines()[-1:])
    try:
        m = json.load(open(os.path.join(tmp, "merged.json"), encoding="utf-8"))
        res["merged"] = dict(complete=m["complete"], counts=m["counts"], code_hash=m["completeness"]["ref"]["code_hash"])
    except (OSError, ValueError, KeyError):
        res["merged"] = None
    shutil.rmtree(tmp, ignore_errors=True)
    return res


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--probe", default=PROBE)
    ap.add_argument("--mini", action="store_true", help="мини-проверка команд задач (n = 96, ~3-6 мин на 4 ядрах)")
    a = ap.parse_args()
    sys.path.insert(0, HERE)
    import cantera, numpy, scipy
    import comb_v2 as V
    probe = json.load(open(a.probe, encoding="utf-8"))
    N, crit, need = volume(probe)
    tv = cpu_per_variant(probe)
    jobs, est_main = series("comb-l0v2", N, SHARDS, JOBS, tv["s_v2_aws"])
    nb = N // 4
    sens = []
    for beta in (0.0, 2.0):
        js, e = series(f"comb-l0v2-b{beta:g}", nb, BETA_SHARDS, BETA_JOBS, tv["s_v2_aws"], beta)
        jobs += js
        sens.append(e)
    pins = [f"cantera=={cantera.__version__}", f"numpy=={numpy.__version__}", f"scipy=={scipy.__version__}"]
    pkg = dict(
        _comment=(f"Камера сгорания KV-3, уровень 0 v2 (решение владельца 15.11.2026: один клапан на топливо, без "
                  f"приводов, частичную нагрузку держит пассивный элемент): основная серия N = {N} (6 задач × 8 частей) и "
                  f"две серии чувствительности к закону затухания при малом J (β = 0 и 2, по N/4 = {nb}, 4 задачи × 12 "
                  f"частей). Раннер берёт только defaults и jobs; остальные блоки - для сессии с ключами. Код - ветка с "
                  f"comb_v2.py, comb_l0_v2.py, merge_l0_v2.py, comb_bc_v2_results.json (code_hash v2 {V.code_hash()})."),
        defaults=dict(instance_type="c7a.16xlarge", fallback_types=["c7i.16xlarge"], image="al2023", disk_gb=30,
                      python="3.12", pip=pins),
        jobs=jobs)
    est = dict(
        date=V.DATE, main=est_main, beta_sensitivity=sens,
        volume=dict(N=N, need=need, criteria=crit, cap=N_CAP,
                    probe=dict(file=os.path.relpath(a.probe, REPO), counts=probe["counts"], code_hash=probe["completeness"]
                               ["ref"]["code_hash"])),
        time_per_variant=tv,
        price=dict(planning=PRICE, zones=ZONES, basis="0,75 долл./ч - верх диапазона задачи владельца 14.11.2026; цены по "
                   "зонам - ledger.jsonl прогона камеры v1 20261115T001745Z-f1da (00:17 UTC 15.11.2026)"),
        usd_expected=round(sum(e["machine_hours_expected"] for e in [est_main] + sens) * 0.7064 * 1.05, 2),
        usd_expected_basis="машино-часы (подготовка + счёт + выгрузка) × 0,7064 долл./ч (us-east-1c) × 1,05 (диск, IPv4) - оценка",
        status="оценка; живой запуск - сессия с ключами (окружение Coldun, навык aws-gtu-runs); CPU-серия в пределах "
               "2000 долл. - без отдельного вопроса (решение владельца 15.11.2026)")
    pkg["estimate"] = est
    pkg["merge"] = dict(
        where="сессия с ключами после run (runs/<run_id>/<задача>/); тот же коммит, что в _job.log задач; Python с "
              + ", ".join(pins) + "; слияние потоковое (~1-2 КБ на строку)",
        commands=[
            "bash -euo pipefail -c 'for S in main b0 b2; do D=/tmp/comb_l0v2_<run_id>_$S; mkdir -p \"$D\"; done; "
            "find routes/gtu/platform/aws/runs/<run_id> -path \"*/_lost\" -prune -o -name \"comb_l0v2_*.json.gz\" -print | "
            "while read -r f; do b=$(basename \"$f\"); case \"$b\" in *_b0_*) S=b0;; *_b2_*) S=b2;; *) S=main;; esac; "
            "cp \"$f\" /tmp/comb_l0v2_<run_id>_$S/; done; ls /tmp/comb_l0v2_<run_id>_main | wc -l'",
            "mkdir -p routes/gtu/combustor/l0_aws_v2",
            f"python3 routes/gtu/combustor/merge_l0_v2.py /tmp/comb_l0v2_<run_id>_main --out "
            f"routes/gtu/combustor/l0_aws_v2/comb_l0v2_KV-3_n{N}_results.json --top 100",
            f"python3 routes/gtu/combustor/merge_l0_v2.py /tmp/comb_l0v2_<run_id>_b0 --out "
            f"routes/gtu/combustor/l0_aws_v2/comb_l0v2_KV-3_n{nb}_b0_results.json --top 30",
            f"python3 routes/gtu/combustor/merge_l0_v2.py /tmp/comb_l0v2_<run_id>_b2 --out "
            f"routes/gtu/combustor/l0_aws_v2/comb_l0v2_KV-3_n{nb}_b2_results.json --top 30"],
        expect=("первая команда печатает 48; каждое слияние - «полнота: True», код 0; код 2 - потеряна или задвоена часть, "
                "индекс, или разные отпечатки (fingerprint, code_hash, bc_key, beta): разбор, не публикация"),
        convergence=("доля новой цели по префиксам 2^k (convergence в итоге слияния): сошлась, если доля на N/2 и N "
                     "различается меньше двух стандартных ошибок; при нуле - верхняя граница 3/N"),
        to_branch="l0_aws_v2/*.json (рейтинг, держатели новой цели, ближайшие, сходимость, Парето) и ledger.jsonl",
        stays_in_s3="части comb_l0v2_*.json.gz - в S3 и runs/<run_id>/ (в .gitignore)")
    json.dump(pkg, open(OUT, "w", encoding="utf-8"), ensure_ascii=False, indent=1)
    pkg["plan_check"] = plan(OUT)
    if a.mini:
        pkg["mini_check"] = mini_check(pkg)
    json.dump(pkg, open(OUT, "w", encoding="utf-8"), ensure_ascii=False, indent=1)
    pc = pkg["plan_check"]
    print(json.dumps(dict(N=N, need=need, nb=nb, s_v2_aws=tv["s_v2_aws"], main=est_main, plan_exit=pc["exit_code"],
                          forecast=pc["forecast_usd"], usd_expected=est["usd_expected"],
                          mini=pkg.get("mini_check", {}).get("merge")), ensure_ascii=False))
    return pc["exit_code"]


if __name__ == "__main__":
    sys.exit(main())
