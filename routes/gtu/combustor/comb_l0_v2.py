"""Камера сгорания KV, уровень 0 v2 (решение владельца 15.11.2026): одна часть выборки - для AWS и локально.

Часть shard из shards получает индексы i с i % shards == shard (разбиение одно, comb_space.shard_indices); вариант -
comb_v2.evaluate (пассивные схемы частичной нагрузки, 100/75/50/30 %, согласованное смешение). Ошибка варианта не
теряет его - строка со статусом «ошибка». В части - отпечатки выборки (fingerprint), кода v2 (code_hash: comb_space,
comb_crn, comb_v2) и границ v2 (bc_key): слияние (merge_l0_v2.py) не примет части с разными отпечатками.

Запуск:
  python comb_l0_v2.py --n 262144 --shard 7 --shards 48 --out "$GTU_OUT" --workers "$NPROC"    # часть на машине AWS
  python comb_l0_v2.py --n 300 --shard 0 --shards 4 --out /tmp/l0v2 --workers 4                # локальная проба
AWS не вызывается.
"""
import argparse, json, os, platform, statistics, sys, time
from multiprocessing import Pool

HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in sys.path:
    sys.path.insert(0, HERE)
import comb_space as CS  # noqa: E402
import comb_crn as CR    # noqa: E402
import comb_v2 as V      # noqa: E402

SCHEMA = "kv.combustor.l0v2.shard/1"


def shard_name(kv, n, shard, shards, beta=None):
    b = "" if beta is None or beta == V.BETA_LOWJ_NOMINAL else f"_b{beta:g}"
    return f"comb_l0v2_{kv}_n{n}{b}_s{shard:03d}of{shards:03d}.json"


def _init(beta=None):
    if beta is not None:
        V.BETA_LOWJ = beta                  # чувствительность к закону затухания при малом J (оценка)
    V.boundary(write=False)
    CR.kin()


def _job(args):
    i, n, kv, seed = args
    t = time.process_time()
    r = V.evaluate(i, n, kv, seed)
    r["cpu_s"] = round(time.process_time() - t, 4)
    return r


def run(kv, n, shard, shards, out, workers=None, seed=V.SEED, quiet=False, beta=None):
    t0 = time.time()
    if beta is not None:
        V.BETA_LOWJ = beta
    idx = CS.shard_indices(n, shard, shards)
    workers = max(1, workers or os.cpu_count() or 1)
    bc = V.boundary()                                    # comb_bc_v2_results.json или расчёт (~10 с)
    tasks = [(i, n, kv, seed) for i in idx]
    step = max(1, len(tasks) // 10)
    rows = []
    if workers == 1:
        it, pool = map(_job, tasks), None
    else:
        pool = Pool(workers, initializer=_init, initargs=(beta,))
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
    cs = [r["cpu_s"] for r in rows]
    st = {}
    for r in rows:
        st[r["status"]] = st.get(r["status"], 0) + 1
    res = dict(
        schema=SCHEMA, kv=kv, n=n, seed=seed, shard=shard, shards=shards, n_rows=len(rows),
        fingerprint=V.fingerprint(n, seed), code_hash=V.code_hash(), bc_key=bc["key"], v1_code_hash=CS.code_hash(),
        beta=V.BETA_LOWJ,
        method="уровень 0 v2 воронки камеры: comb_v2 (пассивные схемы частичной нагрузки, 100/75/50/30 %, согласованное "
               "смешение, гидравлика топлива одним клапаном) поверх comb_space и comb_crn v1",
        software=f"python {sys.version.split()[0]}, cantera {CR.ct.__version__}, numpy {CR.np.__version__}",
        date=V.DATE, host=dict(cpu=os.cpu_count(), workers=workers, machine=platform.machine()),
        runtime_s=round(time.time() - t0, 1),
        t_variant_s=dict(mean=statistics.fmean(ts), median=statistics.median(ts), max=max(ts), sum=sum(ts)) if ts else None,
        cpu_variant_s=dict(mean=statistics.fmean(cs), median=statistics.median(cs), max=max(cs), sum=sum(cs)) if cs else None,
        status_counts=st, rows=rows)
    os.makedirs(out, exist_ok=True)
    path = os.path.join(out, shard_name(kv, n, shard, shards, beta))
    tmp = path + ".tmp"
    with open(tmp, "w", encoding="utf-8") as f:
        json.dump(res, f, ensure_ascii=False, separators=(",", ":"))
    os.replace(tmp, path)
    if not quiet:
        print(f"записано {path}: {len(rows)} вариантов, {res['runtime_s']} с, {st}", flush=True)
    return path


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--kv", default="KV-3", choices=["KV-3"])
    ap.add_argument("--n", type=int, required=True)
    ap.add_argument("--shard", type=int, required=True)
    ap.add_argument("--shards", type=int, required=True, help="число частей (на AWS 10-50)")
    ap.add_argument("--out", required=True)
    ap.add_argument("--workers", type=int, default=None)
    ap.add_argument("--seed", type=int, default=V.SEED)
    ap.add_argument("--beta", type=float, default=None,
                    help="показатель затухания при малом J (по умолчанию comb_v2.BETA_LOWJ = 1; чувствительность - 0 и 2)")
    a = ap.parse_args(argv)
    run(a.kv, a.n, a.shard, a.shards, a.out, a.workers, a.seed, beta=a.beta)


if __name__ == "__main__":
    main()
