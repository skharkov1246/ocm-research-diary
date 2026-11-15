"""Воронка камеры сгорания KV, уровень 0: одна часть выборки (для AWS и локально).

Часть shard из shards получает индексы i с i % shards == shard (comb_space.shard_indices - разбиение одно); каждый
вариант считается comb_crn.evaluate (сети реакторов Cantera и модели модуля); ошибка варианта не теряет его - строка
со статусом «ошибка». Выход - JSON части с отпечатками выборки (fingerprint), кода (code_hash) и границ (bc_key):
слияние (merge_l0.py) не примет части с разными отпечатками.

Запуск:
  python comb_l0.py --kv KV-3 --n 30000 --shard 7 --shards 30 --out "$GTU_OUT"     # одна часть на машине AWS
  python comb_l0.py --n 300 --shard 0 --shards 4 --out /tmp/l0 --workers 4         # локальная отладка
Процессов по умолчанию - число ядер машины (c7a.16xlarge - 64). AWS не вызывается.
"""
import argparse, json, os, platform, statistics, sys, time
from multiprocessing import Pool

HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in sys.path:
    sys.path.insert(0, HERE)
import comb_space as CS  # noqa: E402
import comb_crn as CR    # noqa: E402

SCHEMA = "kv.combustor.l0.shard/1"


def shard_name(kv, n, shard, shards):
    return f"comb_l0_{kv}_n{n}_s{shard:03d}of{shards:03d}.json"


def _init():
    CR.boundary(write=False)
    CR.kin()


def _job(args):
    i, n, kv, seed = args
    return CR.evaluate(i, n, kv, seed)


def run(kv, n, shard, shards, out, workers=None, seed=CS.SEED, quiet=False):
    t0 = time.time()
    idx = CS.shard_indices(n, shard, shards)
    workers = max(1, workers or os.cpu_count() or 1)
    bc = CR.boundary(workers=min(workers, 16))          # из comb_bc_results.json или расчёт (1-3 мин)
    rows = []
    tasks = [(i, n, kv, seed) for i in idx]
    step = max(1, len(tasks) // 10)
    if workers == 1:
        it = map(_job, tasks)
        pool = None
    else:
        pool = Pool(workers, initializer=_init)
        it = pool.imap(_job, tasks, chunksize=1)
    try:
        for k, r in enumerate(it, 1):
            rows.append(r)
            if not quiet and (k % step == 0 or k == len(tasks)):
                print(f"[{shard}/{shards}] {k}/{len(tasks)} за {time.time() - t0:.0f} с", flush=True)
    finally:
        if pool is not None:
            pool.close(); pool.join()
    rows.sort(key=lambda r: r["idx"])
    ts = [r["t_s"] for r in rows]
    st = {}
    for r in rows:
        st[r["status"]] = st.get(r["status"], 0) + 1
    res = dict(
        schema=SCHEMA, kv=kv, n=n, seed=seed, shard=shard, shards=shards, n_rows=len(rows),
        fingerprint=CS.fingerprint(n, seed), code_hash=CS.code_hash(), bc_key=bc["key"],
        method="уровень 0 воронки камеры: comb_space (оси, Соболь, геометрия, правила) + comb_crn (сети реакторов Cantera, "
               "модели модуля, цели targets.json)",
        software=f"python {sys.version.split()[0]}, cantera {CR.ct.__version__}, numpy {CR.np.__version__}",
        date=CS.DATE, host=dict(cpu=os.cpu_count(), workers=workers, machine=platform.machine()),
        runtime_s=round(time.time() - t0, 1),
        t_variant_s=dict(mean=statistics.fmean(ts), median=statistics.median(ts), max=max(ts), sum=sum(ts)) if ts else None,
        status_counts=st, rows=rows)
    os.makedirs(out, exist_ok=True)
    path = os.path.join(out, shard_name(kv, n, shard, shards))
    tmp = path + ".tmp"
    with open(tmp, "w", encoding="utf-8") as f:
        json.dump(res, f, ensure_ascii=False, separators=(",", ":"))
    os.replace(tmp, path)
    if not quiet:
        print(f"записано {path}: {len(rows)} вариантов, {res['runtime_s']} с, {st}", flush=True)
    return path


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--kv", default="KV-3", choices=["KV-3"], help="типоразмер (границы считаются для KV-3)")
    ap.add_argument("--n", type=int, required=True, help="объём выборки")
    ap.add_argument("--shard", type=int, required=True)
    ap.add_argument("--shards", type=int, required=True, help="число частей (на AWS 10-50)")
    ap.add_argument("--out", required=True, help="папка результата ($GTU_OUT на машине AWS)")
    ap.add_argument("--workers", type=int, default=None, help="процессов (по умолчанию - число ядер)")
    ap.add_argument("--seed", type=int, default=CS.SEED)
    a = ap.parse_args(argv)
    run(a.kv, a.n, a.shard, a.shards, a.out, a.workers, a.seed)


if __name__ == "__main__":
    main()
