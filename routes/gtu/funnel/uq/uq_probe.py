"""Устойчивый отбор: локальная миниатюра (4 общих ядра) на пробе уровня 1 - 12 вариантов funnel_l1.csv. Всё - оценка.

Концепция варианта восстанавливается без повтора оптимизации: funnel_space.concept(idx, n, kv, seed) входа пробы (отбор уровня 0
AWS N = 2000) + дискретные и непрерывные оси из строки funnel_l1.csv (оптимум уровня 1). Номинальная выборка сверяется с оценкой
строки (если код уровня 1 с тех пор менялся - расхождение пишется, а не прячется). Затем uq_l1 на двух частях (деление и слияние с
контролем полноты - как на AWS) и merge_uq.

  python uq_probe.py --n 128 --n-od 4 --workers 4      # -> uq_probe_results.json (+ части во временной папке)
"""
import argparse, csv, json, os, shutil, sys, tempfile, time

HERE = os.path.dirname(os.path.abspath(__file__))
FUNNEL = os.path.dirname(HERE)
for p in (FUNNEL, HERE):
    if p not in sys.path:
        sys.path.insert(0, p)
import funnel_space as FS        # noqa: E402
import funnel_l1 as L1           # noqa: E402
import uq_factors as UF          # noqa: E402
import uq_l1 as U                # noqa: E402
import merge_uq as MU            # noqa: E402

PROBE_CSV = os.path.join(FUNNEL, "funnel_l1.csv")
PROBE_N0, PROBE_KV = 2000, "KV-3"            # вход пробы: l0_aws_debug/funnel_l0_select_KV-3_n2000.json
INT_AXES = ("rotors", "n_ax", "nt", "npt", "cores", "hot_integral", "recup_integral", "cartridge", "split_bl", "zc", "zt", "zn")


def _num(v):
    if v in ("", None):
        return None
    try:
        f = float(v)
    except ValueError:
        return v
    return int(f) if f.is_integer() and "." not in v else f


def candidates_from_csv(path=PROBE_CSV, n0=PROBE_N0, kv=PROBE_KV, seed=FS.SEED):
    rows = list(csv.DictReader(open(path, encoding="utf-8")))
    L1._prep(n0, seed, kv)
    out = []
    for r in rows:
        rec = {k: _num(v) for k, v in r.items()}
        c = FS.concept(int(rec["idx"]), n0, kv, seed)
        for a in FS.CAT:
            if a != "arch" and rec.get(a) is not None:
                c[a] = int(rec[a]) if a in INT_AXES else rec[a]
        for a in L1.CONT_AXES:
            if rec.get(a) is not None:
                c[a] = int(round(rec[a])) if a in INT_AXES else float(rec[a])
        keep = {k: rec.get(k) for k in U.REC_KEYS if k in rec}
        out.append(dict(src="L1", j_src=int(rec["j"]), idx=int(rec["idx"]), kv=kv, niche=rec["niche"], arch=rec["arch"], concept=c,
                        vigv=int(rec.get("vigv") or 0), heat=int(rec.get("inlet_heat") or 0), rec=keep))
    return out


def write_candidates(path, cands=None):
    """кандидаты пробы в файл (формат uq_l1 --candidates; затравка pareto_niche --seed-cands)"""
    cands = cands if cands is not None else candidates_from_csv()
    with open(path, "w", encoding="utf-8") as fh:
        json.dump(dict(kind="L1-csv", kv=PROBE_KV, source=os.path.relpath(PROBE_CSV, FUNNEL), candidates=cands), fh, ensure_ascii=False,
                  default=U._r)
    return path


def nominal_check(cands):
    """номинальная выборка против строки пробы: оценка и допустимость расчётной точки"""
    tgt = FS.target(cands[0]["kv"])
    out = []
    with U.UP.Perturb(UF.nominal()):
        for cd in cands:
            r = U.eval_one(cd, tgt, UF.nominal())
            out.append(dict(j=cd["j_src"], idx=cd["idx"], arch=cd["arch"], csv_score=cd["rec"].get("score"), csv_feasible=cd["rec"].get("feasible"),
                            uq_nominal_score_w=U._r(r.get("score_w")), uq_nominal_feas_dp=r["feas"], uq_nominal_bad=r["bad"]))
    return out


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--n", type=int, default=128)
    ap.add_argument("--n-od", type=int, default=4)
    ap.add_argument("--od", default="auto", choices=["auto", "none", "all"])
    ap.add_argument("--od-pdp-min", type=float, default=0.05)
    ap.add_argument("--od-score-min", type=float, default=0.0)
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--shards", type=int, default=2)
    ap.add_argument("--out", default=os.path.join(HERE, "uq_probe_results.json"))
    ap.add_argument("--keep-parts", default=None)
    a = ap.parse_args(argv)
    if a.n > U.LOCAL_MAX["n"] or a.n_od > U.LOCAL_MAX["n_od"]:
        ap.error(f"локально - не больше {U.LOCAL_MAX}")
    t0 = time.time()
    cands = candidates_from_csv()
    nomc = nominal_check(cands)
    meta_in = dict(kind="L1-csv", kv=PROBE_KV, file=os.path.relpath(PROBE_CSV, FUNNEL), n0=PROBE_N0,
                   note="проба уровня 1 (12 вариантов отбора AWS N = 2000, 2 старта); концепции восстановлены по строкам CSV")
    tmp = a.keep_parts or tempfile.mkdtemp(prefix="uq_probe_")
    calls = []
    for k in range(a.shards):          # каждая часть - отдельным вызовом, как машины AWS
        paths, meta, rows = U.run_parts(cands, meta_in, [k], a.shards, tmp, a.n, a.n_od, UF.UQ_SEED, a.od, a.od_pdp_min, a.od_score_min,
                                        a.workers, gz=True)
        calls.append(meta["call"])
    found = MU.find_parts([tmp], PROBE_KV, len(cands), a.shards)
    out = MU.run_merge(found, PROBE_KV, len(cands), a.shards, tmp, "_probe")
    # контроль полноты: без одной части слияние отказывает
    miss_ok = False
    try:
        MU.merge_parts(found[1:], PROBE_KV, len(cands), a.shards)
    except MU.MergeError:
        miss_ok = True
    res = out["res"]
    res["title"] = "Устойчивый отбор: локальная миниатюра на пробе уровня 1 (12 вариантов)"
    res["status"] = (f"локальная миниатюра (4 общих ядра): N = {a.n} выборок на вариант, нерасчётные режимы - первые {a.n_od} строк у "
                     f"вариантов с P(расчётная) >= {a.od_pdp_min}; на AWS - N 500, N_od 100 (jobs_uq*.json)")
    res["probe"] = dict(params=vars(a), nominal_check=nomc, merge_rejects_missing_part=miss_ok, wall_s=round(time.time() - t0, 1),
                        rows=[dict(u=r["u"], j=r["j_src"], idx=r["idx"], arch=r["arch"], P_dp=r["dp"]["P_dp"], P_dp_ci=r["dp"]["P_dp_ci"],
                                   P_od_dp=(r.get("od") or {}).get("P_od_dp"), P_feas=r.get("P_feas"), od_status=r.get("od_status"),
                                   score_f=r["dp"]["score_f"], nominal=r["dp"]["nominal"], viol_freq=r["dp"]["viol_freq"],
                                   sobol_y=r["dp"]["sobol"]["y"], sobol_score=r["dp"]["sobol"]["score"], sobol_noise=r["dp"]["sobol"]["noise"],
                                   rules=r["dp"]["rules"], half=r["dp"]["half"])
                              for r in out["rows"]])
    with open(a.out, "w", encoding="utf-8") as fh:
        json.dump(res, fh, ensure_ascii=False, indent=1, default=U._r)
    if not a.keep_parts:
        shutil.rmtree(tmp, ignore_errors=True)
    s = res["summary"]
    print(f"миниатюра: {s['variants']} вариантов, номинально допустимы в расчётной точке {s['nominal_dp_feasible']}, нерасчётные посчитаны у "
          f"{s['od_done']}; без части слияние отказывает: {miss_ok}; {round(time.time() - t0, 1)} с -> {a.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
