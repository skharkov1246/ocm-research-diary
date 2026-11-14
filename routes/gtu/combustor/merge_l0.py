"""Воронка камеры сгорания KV, уровень 0: слияние частей, контроль полноты, рейтинг, фронт Парето, чувствительность.

Проверки: у всех частей один объём выборки, зерно, отпечаток матрицы Соболя (fingerprint), кода (code_hash) и границ
(bc_key); номера частей 0..S-1 по одному разу (повтор части после прерывания спота - одинаковые индексы, берётся одна);
объединение индексов = range(n) без повторов. Неполное слияние пишется с complete=false и кодом выхода 2.

Рейтинг - взвешенная сумма нормированных оценок по приоритетам владельца (14.11.2026): минимум деталей (главное),
Plug&Play, монтаж «хоть уборщиком», затем КПД (потери давления), размер, изготовимость, NOx; умножается на штраф
жёстких целей targets.json: C = exp(-2·Σ min(нарушение, 3)). Ранжируются только допустимые варианты (status = ok).
Фронт Парето - по шести целям (детали, Plug&Play, обслуживание, потери, объём, наибольший NOx газа) среди допустимых.
Чувствительность к весам: каждый вес x0,5 и x1,5 (с нормировкой) и 300 случайных наборов Дирихле (фиксированное зерно).

Запуск:
  python merge_l0.py runs/<run_id> --out comb_l0_KV-3_results.json        # части из папки (рекурсивно)
  python merge_l0.py /tmp/l0 --out /tmp/l0/merged.json --top 30
"""
import argparse, glob, json, math, os, sys, time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in sys.path:
    sys.path.insert(0, HERE)

DATE = "2026-11-14"
WEIGHTS = dict(parts=0.30, plugplay=0.20, service=0.15, eff=0.12, size=0.08, make=0.10, nox=0.05)
WEIGHTS_NOTE = ("приоритеты владельца 14.11.2026: детали (главное) > Plug&Play > монтаж и ремонт > КПД > размер; "
                "изготовимость и NOx - малыми весами (выбросы держат жёсткие цели); значения весов - оценка")
PENALTY_K = 2.0
PARETO_OBJ = [("parts", "min"), ("plugplay", "max"), ("service", "max"), ("dp", "min"), ("V_env", "min"), ("NOx", "min")]


def load_parts(paths):
    files = []
    for p in paths:
        if os.path.isdir(p):
            files += glob.glob(os.path.join(p, "**", "comb_l0_*.json"), recursive=True)
        else:
            files.append(p)
    files = sorted(set(f for f in files if not f.endswith(".tmp")))
    return [(f, json.load(open(f, encoding="utf-8"))) for f in files]


def check_complete(parts):
    """полнота и согласованность частей -> (complete, отчёт, строки по индексу)"""
    rep = dict(files=len(parts), problems=[])
    if not parts:
        rep["problems"].append("нет частей")
        return False, rep, {}
    keys = ("schema", "kv", "n", "seed", "shards", "fingerprint", "code_hash", "bc_key")
    ref = {k: parts[0][1].get(k) for k in keys}
    rep["ref"] = ref
    for f, p in parts:
        for k in keys:
            if p.get(k) != ref[k]:
                rep["problems"].append(f"{os.path.basename(f)}: {k}={p.get(k)!r} != {ref[k]!r}")
    n, S = ref["n"], ref["shards"]
    seen_shards = {}
    for f, p in parts:
        seen_shards.setdefault(p["shard"], []).append(f)
    rep["duplicate_shards"] = {k: v for k, v in seen_shards.items() if len(v) > 1}
    rep["missing_shards"] = sorted(set(range(S)) - set(seen_shards))
    rows = {}
    dup_idx = 0
    used = set()
    for f, p in parts:
        if p["shard"] in used:
            continue                                  # повтор части (спот): одна копия
        used.add(p["shard"])
        expect = set(range(p["shard"], n, S))
        got = [r["idx"] for r in p["rows"]]
        if set(got) != expect or len(got) != len(expect):
            rep["problems"].append(f"{os.path.basename(f)}: индексы части не совпадают с разбиением i % {S} == {p['shard']}")
        for r in p["rows"]:
            if r["idx"] in rows:
                dup_idx += 1
            rows[r["idx"]] = r
    rep["duplicate_indices"] = dup_idx
    missing = sorted(set(range(n)) - set(rows))
    rep["missing_indices"] = len(missing)
    rep["missing_indices_head"] = missing[:20]
    complete = not rep["problems"] and not rep["missing_shards"] and not missing and dup_idx == 0
    return complete, rep, rows


def _clip(x):
    return max(0.0, min(1.0, x))


def objectives(r):
    """нормированные оценки 0..1 (больше - лучше) и сырые величины для Парето"""
    m, c = r["module"], r["concept"]
    parts = m["parts"]
    V_env = math.pi / 4 * r["derived"]["D_env"] ** 2 * r["derived"]["L_mod"]
    nox = max(g["NOx"] for g in r["gas"])
    dp = r["dp"]["dp_total"]
    s = dict(parts=_clip(1 - math.log(max(parts, 3) / 3) / math.log(40 / 3)),
             plugplay=m["plugplay_score"] / 10, service=m["service_score"] / 10,
             eff=_clip((0.08 - dp) / 0.06), size=_clip(1 - math.log(max(V_env, 0.02) / 0.02) / math.log(3.0 / 0.02)),
             make=m["make_score"] / 10, nox=_clip(1 - nox / 30.0))
    raw = dict(parts=parts, plugplay=m["plugplay_score"], service=m["service_score"], dp=dp, V_env=V_env, NOx=nox)
    return s, raw


def penalty(r):
    hard = [t for t, v in r["targets"].items() if t in HARD]
    return math.exp(-PENALTY_K * sum(min(r["targets"][t]["viol"], 3.0) for t in hard))


def score(s, pen, w):
    tot = sum(w.values())
    return pen * sum(w[k] * s[k] for k in w) / tot


def pareto(F):
    """недоминируемые строки F (все цели - к минимуму)"""
    n = len(F)
    keep = np.ones(n, bool)
    for i in range(n):
        if not keep[i]:
            continue
        le = np.all(F <= F[i], axis=1) & np.any(F < F[i], axis=1)
        if le.any():
            keep[i] = False
            continue
        dom = np.all(F[i] <= F, axis=1) & np.any(F[i] < F, axis=1)
        keep[dom] = False
    return np.where(keep)[0]


def compact(r, sc=None, pen=None):
    c, g, dz, m = r["concept"], r["gas"], r["diesel"], r["module"]
    return dict(idx=r["idx"], score=sc, penalty=pen, scheme=c["scheme"], chamber=c["chamber"], injector=c["injector"],
                diesel=c["diesel"], cooling=c["cooling"], make=c["make"], passive=c["passive"],
                n_el=c["n_el"], n_jet=r["derived"]["n_jet"], d_jet_mm=r["derived"]["d_jet_mm"], d_air_mm=r["derived"]["d_air_mm"],
                a_pz=c["a_pz"], a_aux=c["a_aux"], k_pass=c["k_pass"], tau_pz_ms=c["tau_pz_ms"], tau_bo_ms=c["tau_bo_ms"],
                dp_inj=c["dp_inj"], dp_total=r["dp"]["dp_total"], NOx=[x["NOx"] for x in g], CO=[x["CO"] for x in g],
                margin_lbo=[x["margin_lbo"] for x in g], NOx_diesel=[x["NOx"] for x in dz], CO_diesel=[x["CO"] for x in dz],
                H2_vol_max=r["flashback"]["H2_vol_max"], PF=r["pf"]["PF"], T_wall_K=r["feed"]["T_metal_max_K"],
                parts=m["parts"], swap_h=m["swap_h"], mass_kg=m["mass_kg"], plugplay=m["plugplay_score"],
                service=m["service_score"], make_score=m["make_score"], D_env=r["derived"]["D_env"], L_mod=r["derived"]["L_mod"],
                V_L=r["derived"]["V_L"], hard_ok=r["hard_ok"],
                hard_fail=[t for t, v in r["targets"].items() if t in HARD and not v["ok"]])


HARD = []


def merge(paths, out, top=50):
    t0 = time.time()
    parts = load_parts(paths)
    complete, rep, rows = check_complete(parts)
    import comb_crn as CR   # noqa: E402  (цели и веса - из того же кода, что считал части)
    global HARD
    HARD = [t for t in CR.TARGETS["funnel_interface"]["hard"]]
    ok = [r for r in rows.values() if r["status"] == "ok"]
    status, reasons = {}, {}
    for r in rows.values():
        status[r["status"]] = status.get(r["status"], 0) + 1
        for b in r.get("reject", []):
            reasons[b["code"]] = reasons.get(b["code"], 0) + 1
    S, RAW, PEN = [], [], []
    for r in ok:
        s, raw = objectives(r)
        S.append(s); RAW.append(raw); PEN.append(penalty(r))
    wk = list(WEIGHTS)
    Smat = np.array([[s[k] for k in wk] for s in S]) if S else np.zeros((0, len(wk)))
    PENv = np.array(PEN)
    vscore = lambda w: PENv * (Smat @ np.array([w[k] for k in wk])) / sum(w.values()) if len(PENv) else np.array([])
    base = vscore(WEIGHTS)
    order = np.argsort(-base)
    top_rows = [compact(ok[i], float(base[i]), PEN[i]) for i in order[:top]]
    # фронт Парето
    if ok:
        F = np.array([[(raw[k] if sense == "min" else -raw[k]) for k, sense in PARETO_OBJ] for raw in RAW])
        hard_ok = np.array([r["hard_ok"] for r in ok])
        pool = np.where(hard_ok)[0] if hard_ok.any() else np.arange(len(ok))
        pf_idx = pool[pareto(F[pool])]
        pf_rows = sorted([compact(ok[i], float(base[i]), PEN[i]) for i in pf_idx], key=lambda x: -x["score"])
    else:
        pf_rows, hard_ok = [], np.array([])
    # чувствительность к весам
    sens = []
    top10 = set(int(i) for i in order[:10])
    best = int(order[0]) if len(order) else None
    for k in (WEIGHTS if ok else []):
        for f in (0.5, 1.5):
            w = dict(WEIGHTS); w[k] *= f
            sc = vscore(w)
            o = np.argsort(-sc)
            sens.append(dict(weight=k, factor=f, top10_overlap=len(top10 & set(int(i) for i in o[:10])) / 10,
                             best_idx=ok[int(o[0])]["idx"], base_best_rank=int(np.where(o == best)[0][0]) + 1))
    rng = np.random.default_rng(20261114)
    freq = {}
    for _ in range(300 if ok else 0):
        ww = rng.dirichlet([20 * WEIGHTS[k] for k in wk])
        sc = vscore(dict(zip(wk, ww)))
        for i in np.argsort(-sc)[:10]:
            freq[ok[int(i)]["idx"]] = freq.get(ok[int(i)]["idx"], 0) + 1
    robust = sorted(freq.items(), key=lambda x: -x[1])[:15]
    # по схемам
    by = {}
    for r in rows.values():
        sch = r["concept"]["scheme"]
        b = by.setdefault(sch, dict(total=0, ok=0, hard_ok=0, reject={}, best_idx=None, best_score=None))
        b["total"] += 1
        for x in r.get("reject", []):
            b["reject"][x["code"]] = b["reject"].get(x["code"], 0) + 1
    for j, r in enumerate(ok):
        b = by[r["concept"]["scheme"]]
        b["ok"] += 1
        b["hard_ok"] += int(r["hard_ok"])
        if b["best_score"] is None or base[j] > b["best_score"]:
            b["best_score"], b["best_idx"] = float(base[j]), r["idx"]
    for sch, b in by.items():
        sub = [r for r in ok if r["concept"]["scheme"] == sch]
        if sub:
            b["NOx_100_median"] = float(np.median([r["gas"][0]["NOx"] for r in sub]))
            b["NOx_30_median"] = float(np.median([r["gas"][2]["NOx"] for r in sub]))
            b["CO_30_median"] = float(np.median([r["gas"][2]["CO"] for r in sub]))
            b["diesel_NOx_100_median"] = float(np.median([r["diesel"][0]["NOx"] for r in sub]))
            b["parts_min"] = int(min(r["module"]["parts"] for r in sub))
            b["dp_min"] = float(min(r["dp"]["dp_total"] for r in sub))
            b["target_pass_share"] = {t: float(np.mean([r["targets"][t]["ok"] for r in sub])) for t in sub[0]["targets"]}
    # что отдать воронке двигателя: лучший по рейтингу в каждой паре (схема, тип камеры)
    feed = {}
    for j in order:
        r = ok[int(j)]
        key = f'{r["concept"]["scheme"]}|{r["concept"]["chamber"]}'
        if key not in feed:
            feed[key] = dict(idx=r["idx"], score=float(base[int(j)]), hard_ok=r["hard_ok"], **r["feed"])
    ts = [r["t_s"] for r in rows.values() if "t_s" in r]
    shard_rt = [p["runtime_s"] for _, p in parts]
    res = dict(
        title="Камера сгорания KV-3, уровень 0: слияние частей, рейтинг и фронт Парето",
        method=("слияние частей comb_l0.py с контролем полноты; рейтинг - взвешенная сумма нормированных оценок × штраф "
                "жёстких целей exp(-2·Σ нарушений); фронт Парето по шести целям; чувствительность - веса x0,5/x1,5 и "
                "300 наборов Дирихле"),
        software=f"python {sys.version.split()[0]}, numpy {np.__version__}", date=DATE,
        status="уровень облика (0D + корреляции): числа - оценка для отбора на L1 (RANS), не прогноз конструкции",
        complete=complete, completeness=rep,
        inputs=[dict(file=os.path.relpath(f, HERE) if f.startswith(HERE) else os.path.basename(f), shard=p["shard"],
                     rows=p["n_rows"], runtime_s=p["runtime_s"]) for f, p in parts],
        weights=WEIGHTS, weights_note=WEIGHTS_NOTE, penalty_k=PENALTY_K, hard_targets=HARD,
        counts=dict(n=rep.get("ref", {}).get("n"), rows=len(rows), status=status, reject_codes=reasons,
                    feasible=len(ok), hard_ok=int(hard_ok.sum()) if len(ok) else 0),
        timing=dict(t_variant_mean_s=float(np.mean(ts)) if ts else None, t_variant_median_s=float(np.median(ts)) if ts else None,
                    t_variant_max_s=float(np.max(ts)) if ts else None, t_variant_sum_s=float(np.sum(ts)) if ts else None,
                    shard_runtime_s=shard_rt),
        by_scheme=by, top=top_rows, pareto=pf_rows, pareto_objectives=PARETO_OBJ,
        sensitivity=dict(one_at_a_time=sens, dirichlet_top10_freq=[dict(idx=i, freq=f / 300) for i, f in robust]),
        feed_to_engine_funnel=dict(note="лучший вариант по рейтингу в каждой паре (схема, тип камеры); поля - "
                                        "targets.json funnel_interface; замена заглушек funnel_space.py (dp 4/5 %, "
                                        "D_comb = 1,1·Dc, L_comb, m_comb, n_inj)", best=feed),
        runtime_s=None)
    res["runtime_s"] = round(time.time() - t0, 1)
    import comb_crn as _cr
    json.dump(_cr._rnd(res), open(out, "w", encoding="utf-8"), ensure_ascii=False, indent=1)
    return complete, res


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("paths", nargs="+", help="папки или файлы частей")
    ap.add_argument("--out", required=True)
    ap.add_argument("--top", type=int, default=50)
    a = ap.parse_args(argv)
    complete, res = merge(a.paths, a.out, a.top)
    c = res["counts"]
    print(f"полнота: {complete}; строк {c['rows']} из {c['n']}; допустимых {c['feasible']}; все жёсткие цели - {c['hard_ok']}; "
          f"запись {a.out}")
    return 0 if complete else 2


if __name__ == "__main__":
    sys.exit(main())
