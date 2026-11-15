"""Пакет задач уровня 0 воронки камеры для раннера AWS (routes/gtu/platform/aws/runner.py) и смета.

Собирает jobs_l0.json: N вариантов KV-3 частями (по задаче на часть, имя comb-l0-pNN), образ al2023, Python через uv,
пакеты закреплены как в локальной отладке (тот же ключ границ comb_bc_results.json - части не пересчитывают таблицы).
Смета: время варианта - из локальной отладки (comb_l0_debug_results.json, timing); прогноз пакета - вывод команды
`runner.py plan` (без AWS) при допущенной цене спота. AWS не вызывается.

Запуск: python make_jobs_l0.py [--n 30000 --shards 30 --price 0.955 --price 1.5]
"""
import argparse, json, os, re, subprocess, sys

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
RUNNER = os.path.join(REPO, "routes", "gtu", "platform", "aws", "runner.py")
OUT = os.path.join(HERE, "jobs_l0.json")
DEBUG = os.path.join(HERE, "comb_l0_debug_results.json")

SETUP_MIN = 10.0     # подготовка машины: uv + Python + pip cantera/numpy/scipy + клон (оценка; smoke2 - срок 25 мин на всё)
UPLOAD_MIN = 8.0     # запас раннера на выгрузку (README раннера: команда идёт до «старт + max_min − 8 мин»)
BC_MIN = 6.0         # пересчёт границ, если ключ не совпал (локально 4,5 мин на 4 ядрах; на машине - 16 процессов)
SAFETY = 3.0         # запас на счёт частей
VCPU = 64            # c7a.16xlarge


def build(n, shards, max_min=None):
    dbg = json.load(open(DEBUG, encoding="utf-8"))
    t_mean = dbg["timing"]["t_variant_mean_s"]
    per = math_ceil(n / shards)
    compute_min = per * t_mean / VCPU / 60.0
    need = SETUP_MIN + BC_MIN + SAFETY * compute_min + UPLOAD_MIN
    mm = max_min or max(20, int(5 * math_ceil(need / 5)))
    import cantera, numpy, scipy            # версии - как у локальной отладки (ключ границ включает версию Cantera)
    pins = [f"cantera=={cantera.__version__}", f"numpy=={numpy.__version__}", f"scipy=={scipy.__version__}"]
    py = f"{sys.version_info.major}.{sys.version_info.minor}"
    jobs = [dict(name=f"comb-l0-p{k:02d}",
                 command=f"python -u routes/gtu/combustor/comb_l0.py --kv KV-3 --n {n} --shard {k} --shards {shards} "
                         f"--out \"$GTU_OUT\" --workers \"$NPROC\"")
            for k in range(shards)]
    pkg = dict(defaults=dict(instance_type="c7a.16xlarge", fallback_types=["c7i.16xlarge"], image="al2023", disk_gb=30,
                             max_min=mm, python=py, pip=pins),
               jobs=jobs)
    est = dict(
        title=f"Уровень 0 воронки камеры KV-3 на AWS: N = {n}, частей {shards}",
        date="2026-11-14",
        method=("время варианта - локальная отладка (comb_l0_debug_results.json, N = 300, 4 части, 4 ядра Xeon 2,1 ГГц); "
                "счёт части = (N/частей)·t_вар/64 процесса; срок сторожа = подготовка + границы + 3·счёт + выгрузка; "
                "деньги - прогноз команды runner.py plan (верхняя граница по формуле раннера) при допущенной цене спота"),
        inputs=[dict(file="routes/gtu/combustor/comb_l0_debug_results.json", field="timing.t_variant_mean_s",
                     used=t_mean)],
        per_shard_variants=per, workers=VCPU, compute_min_per_shard=compute_min,
        setup_min=SETUP_MIN, bc_recompute_min=BC_MIN, upload_min=UPLOAD_MIN, safety=SAFETY, max_min=mm,
        core_hours_total=n * t_mean / 3600.0,
        merge="после сбора: python routes/gtu/combustor/merge_l0.py runs/<run_id> --out "
              "routes/gtu/combustor/comb_l0_KV-3_results.json (контроль полноты; код 2 - части не хватает)",
        status="оценка; живой запуск - отдельная сессия с ключами (окружение Coldun)")
    return pkg, est


def math_ceil(x):
    return int(-(-x // 1))


def plan(path, price):
    """прогноз пакета командой runner.py plan (без AWS): число долларов из строки «[план] прогноз пакета X долл.»"""
    r = subprocess.run([sys.executable, RUNNER, "plan", path, "--assume-price", str(price)], capture_output=True,
                       text=True, cwd=REPO, timeout=600)
    m = re.search(r"прогноз пакета ([0-9.]+) долл", r.stdout + r.stderr)
    return dict(assume_price_usd_h=price, exit_code=r.returncode, forecast_usd=float(m.group(1)) if m else None,
                summary=[ln for ln in (r.stdout + r.stderr).splitlines() if ln.startswith("[") and "Token" not in ln])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=30000)
    ap.add_argument("--shards", type=int, default=30)
    ap.add_argument("--max-min", type=int, default=None)
    ap.add_argument("--price", type=float, action="append",
                    help="допущенная цена спота c7a.16xlarge, долл./ч (по умолчанию 0,955 и 1,5)")
    a = ap.parse_args()
    if not 10 <= a.shards <= 50:
        sys.exit("частей должно быть 10-50 (правило прогонов)")
    pkg, est = build(a.n, a.shards, a.max_min)
    json.dump(pkg, open(OUT, "w", encoding="utf-8"), ensure_ascii=False, indent=1)
    est["plan"] = [plan(OUT, p) for p in (a.price or [0.955, 1.5])]
    est["price_basis"] = ("0,955 - наибольшая цена спота по зонам в живом пробном прогоне раннера 14.11.2026 (навык "
                          "aws-gtu-runs: 0,706 в us-east-1c, 0,955 в us-east-1a); 1,5 - запасное допущение README раннера")
    est["ceiling_usd"] = 2000.0
    est["share_of_ceiling"] = [round(p["forecast_usd"] / 2000.0, 4) if p["forecast_usd"] else None for p in est["plan"]]
    pkg["estimate"] = est            # раннер читает только defaults и jobs; смета - рядом
    json.dump(pkg, open(OUT, "w", encoding="utf-8"), ensure_ascii=False, indent=1)
    print(json.dumps({k: est[k] for k in ("per_shard_variants", "compute_min_per_shard", "max_min", "core_hours_total")},
                     ensure_ascii=False), [(p["assume_price_usd_h"], p["forecast_usd"], p["exit_code"]) for p in est["plan"]])


if __name__ == "__main__":
    main()
