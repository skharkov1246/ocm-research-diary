"""Устойчивый отбор, часть 1: Монте-Карло по неопределённостям моделей уровня 1 для каждого варианта (на AWS). Всё - оценка.

Вариант - концепция уровня 1 с её расчётной точкой (оптимум уровня 1: части funnel_l1 с разборами, поле concept) или член фронта
Парето ниши (pareto_niche.py -> merge_pareto.py -> файл кандидатов). Геометрия варианта зафиксирована - меняются только величины,
которые модели знают с разбросом (uq_factors.FACTORS: модели потерь, поправки КПД, материалы, допуски изготовления, температура
газа, TRL, прирост давления RDC, признаки срыва). Выборка - латинский гиперкуб N строк (общая для всех вариантов, зерно UQ_SEED).

  1. Расчётная точка на всех N выборках: funnel_l1.evaluate_point (циклы, одномерные ступени, прочность, все правила расчётной
     точки) в подстановке uq_patch.Perturb; узлы ВНА и подогрева входа - как у варианта (add_vigv); оценка - взвешенная по весам
     владельца × множитель готовности 0,5 + 0,5·TRL/9 (решение 15.11.2026) с ошибкой TRL.
  2. Нерасчётные режимы (offdesign_l1 через funnel_l1.full_check: −30…+45 °C × 50-100 %, помпаж, пуск) - на первых N_od строках
     гиперкуба, где расчётная точка прошла; только у вариантов, прошедших порог (--od auto: P(расчётная) >= OD_PDP_MIN и P50
     оценки >= SCORE_MIN - OD_SCORE_MARGIN), - один расчёт около 30 с.
  3. Итог варианта: P(допустим) = P(расчётная точка) × P(нерасчётные | расчётная) с интервалами Уилсона, P10/P50/P90 оценки и КПД,
     частоты нарушенных правил, индексы Соболя первого порядка (по той же выборке), сходимость - те же числа на первых N/2 строках.
Деление на части - единственное: вариант u идёт в часть u % S. Результат не зависит от числа процессов.

  python uq_l1.py --l1-in aws/l1_full --l1-kv KV-3 --l1-m 3983 --l1-shards 40 --n 500 --n-od 100 --shards 50 --shard 0,6,12 \
      --workers "$NPROC" --gzip --out "$GTU_OUT"
  python uq_l1.py --candidates pareto_cands_KV-1.json --n 500 --n-od 64 --od-pdp-min 0.5 --shards 40 --shard all --out ...
Слияние, устойчивый отбор 50 и сравнение с обычным - merge_uq.py.
"""
import argparse, datetime, gzip, hashlib, json, math, os, sys, time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
FUNNEL = os.path.dirname(HERE)
for p in (FUNNEL, HERE):
    if p not in sys.path:
        sys.path.insert(0, p)
import funnel_space as FS        # noqa: E402
import funnel_l0 as L0           # noqa: E402
import funnel_l1 as L1           # noqa: E402
import merge_l1 as M1            # noqa: E402
import uq_factors as UF          # noqa: E402
import uq_patch as UP            # noqa: E402
import uq_stats as US            # noqa: E402

LEVEL = "UQ"
PART_PREFIX = "uq"
DATE = "2026-11-15"
DEFAULT_N = 500
DEFAULT_N_OD = 100
OD_PDP_MIN = 0.2              # порог P(расчётная точка) для нерасчётных режимов (оценка: ниже - вариант не отберётся при P_MIN 0,5)
OD_SCORE_MARGIN = 0.5         # и P50 оценки не ниже SCORE_MIN - 0,5
SCORE_MIN = M1.SCORE_MIN
LOCAL_MAX = dict(variants=12, n=128, n_od=4)      # локально (4 общих ядра): не больше 12 вариантов, 128 выборок, 4 нерасчётных
UQ_FILES = ("uq_factors.py", "uq_patch.py", "uq_stats.py", "uq_l1.py")


def trl_factor(t):
    f = getattr(L0, "trl_factor", None)
    if f is not None:
        return float(f(t))
    return 0.5 + 0.5 * float(t) / 9.0            # решение владельца 15.11.2026 (routes/gtu/README.md)


def code_sha():
    h = hashlib.sha256()
    for f in UQ_FILES:
        h.update(open(os.path.join(HERE, f), "rb").read())
    h.update(L1.code_sha().encode())
    return h.hexdigest()[:16]


# =========================================================================== кандидаты
REC_KEYS = ["j", "idx", "kv", "arch", "niche", "feasible", "score", "eta_th_shaft", "eta_el", "V_env", "mass_kg", "t41", "pi", "eps_rec",
            "trl_min", "vigv", "inlet_heat", "reject", "n_parts", "n_assemblies", "steps", "L_env", "W_env", "H_env"] + list(M1.DISC_W)


def load_l1(dirs, kv, M, shards, tag=""):
    """части уровня 1 (funnel_l1 с разборами) -> кандидаты (варианты с построенной расчётной точкой - у них есть concept);
    полнота частей - merge_l1.merge_parts"""
    found = M1.find_parts(dirs, kv, M, shards, tag)
    cols, rows, det, ref, checks, calls = M1.merge_parts(found, kv, M, shards)
    cands = []
    for r in rows:
        rec = dict(zip(cols, r))
        d = det.get(rec["j"])
        if not d or not d.get("concept"):
            continue
        cands.append(dict(src="L1", j_src=rec["j"], idx=rec["idx"], kv=kv, niche=rec["niche"], arch=rec["arch"], concept=d["concept"],
                          vigv=int(rec.get("vigv") or 0), heat=int(rec.get("inlet_heat") or 0),
                          rec={k: rec.get(k) for k in REC_KEYS if k in rec}))
    meta = dict(kind="L1", kv=kv, M=M, shards=shards, tag=tag, l1_code_sha=ref.get("code_sha"), l1_input_sha256=ref.get("input_sha256"),
                l1_order_sha=ref.get("order_sha"), l1_rows=len(rows), merge_checks={k: v for k, v in checks.items() if k != "duplicates"})
    return cands, meta


def load_candidates(path):
    d = json.load(gzip.open(path, "rt", encoding="utf-8") if path.endswith(".gz") else open(path, encoding="utf-8"))
    meta = dict(kind=d.get("kind", "pareto"), kv=d["kv"], file=os.path.basename(path), source=d.get("source"))
    return list(d["candidates"]), meta


def cand_sha(cands):
    key = [(c["src"], c.get("j_src"), c.get("key"), sorted((k, v) for k, v in c["concept"].items() if not isinstance(v, (dict, list))))
           for c in cands]
    return hashlib.sha256(json.dumps(key, ensure_ascii=False, default=str).encode()).hexdigest()[:16]


# =========================================================================== расчёт одной выборки
def _strength_min(stg):
    bursts, lives = [w["burst_margin_120"] for w in stg["comp"] if w.get("burst_margin_120") is not None], []
    for x in (stg["turb"], stg["pt"]):
        if not x:
            continue
        if x.get("burst_margin_120") is not None:
            bursts.append(x["burst_margin_120"])
        if x.get("overspeed_ratio"):
            bursts.append(1.0 / x["overspeed_ratio"])
        if x.get("life_h") is not None:
            lives.append(x["life_h"])
    return (min(bursts) if bursts else None), (min(lives) if lives else None)


def eval_one(cand, tgt, v, od=False):
    """вариант при значениях факторов v (внутри Perturb(v)): расчётная точка, при od - нерасчётные режимы; словарь чисел"""
    c = dict(cand["concept"])
    if c.get("comb") in ("RDC", "ВР") and c.get("pgain"):
        c["pgain"] = float(c["pgain"]) * float(v["k_pgain"])
    UP.set_variant(c)
    ev = L1.evaluate_point(c, tgt)
    if not ev["ok"]:
        return dict(ok=0, feas=0, viol=round(float(ev["viol"]), 6), bad=[r["code"] for r in ev["bad"]])
    out = dict(ok=1)
    vigv, heat = bool(cand.get("vigv")), bool(cand.get("heat"))
    if od and ev["viol"] == 0:
        mx, stt, R = L1.full_check(c, tgt, ev)
        odbad = [r["code"] for r in R if not r["ok"]]
        vigv, heat = bool(mx["vigv_used"]), bool(mx["heat_used"])
        out.update(od_ok=int(not odbad), od_bad=odbad, SM_min=mx["SM_min"], vigv=int(vigv), heat=int(heat),
                   n_ss=stt["n_self_sustain"])
    if vigv or heat:
        ev = L1.add_vigv(c, tgt, ev, vigv=vigv, heat=heat)
    sw = float(L0.weighted({k: ev["scores"][k] for k in L0.CRIT}))
    trl = int(ev["models"]["mr"]["trl_min"])
    trl_s = min(9, max(1, trl + int(v["d_trl"])))
    bmin, lmin = _strength_min(ev["strength"])
    out.update(feas=int(ev["viol"] == 0), viol=round(float(ev["viol"]), 6), bad=[r["code"] for r in ev["bad"]], score_w=sw,
               score_f=sw * trl_factor(trl_s), trl=trl, eta_el=float(ev["eta_el"]), V_env=float(ev["layout"]["envelope"]["V"]),
               burst_min=bmin, life_min=lmin)
    return out


# =========================================================================== пул: расчётная точка (выборки снаружи, варианты внутри)
_G = {}


def _dp_chunk(ids):
    cands, vals, tgt = _G["cands"], _G["vals"], _G["tgt"]
    res = {i: [None] * len(vals) for i in ids}
    nom = {}
    t0 = time.process_time()
    with UP.Perturb(UF.nominal()):
        for i in ids:
            nom[i] = eval_one(cands[i], tgt, UF.nominal())
    for s, v in enumerate(vals):
        with UP.Perturb(v):
            for i in ids:
                res[i][s] = eval_one(cands[i], tgt, v)
    return res, nom, time.process_time() - t0


def _od_task(args):
    i, s = args
    v = _G["vals"][s]
    t0 = time.process_time()
    with UP.Perturb(v):
        r = eval_one(_G["cands"][i], _G["tgt"], v, od=True)
    return i, s, r, time.process_time() - t0


def _pool(workers):
    import multiprocessing as mp
    from concurrent.futures import ProcessPoolExecutor
    ctx = mp.get_context("fork") if "fork" in mp.get_all_start_methods() else None
    return ProcessPoolExecutor(workers, mp_context=ctx)


def run_dp(cands, ids, vals, tgt, workers):
    _G.update(cands=cands, vals=vals, tgt=tgt)
    ids = list(ids)
    if not ids:
        return {}, {}, 0.0
    nch = max(1, min(len(ids), workers * 3))
    chunks = [ids[k::nch] for k in range(nch)]
    out, nom, cpu = {}, {}, 0.0
    if workers <= 1:
        parts = [_dp_chunk(ch) for ch in chunks]
    else:
        with _pool(workers) as ex:
            parts = list(ex.map(_dp_chunk, chunks))
    for r, n, t in parts:
        out.update(r); nom.update(n); cpu += t
    return out, nom, cpu


def run_od(tasks, cands, vals, tgt, workers):
    _G.update(cands=cands, vals=vals, tgt=tgt)
    if not tasks:
        return {}, 0.0
    if workers <= 1:
        got = [_od_task(t) for t in tasks]
    else:
        with _pool(workers) as ex:
            got = list(ex.map(_od_task, tasks, chunksize=1))
    out, cpu = {}, 0.0
    for i, s, r, t in got:
        out[(i, s)] = r
        cpu += t
    return out, cpu


# =========================================================================== итог варианта
def _freq(lists, n):
    cnt = {}
    for b in lists:
        for code in set(b):
            cnt[code] = cnt.get(code, 0) + 1
    return {k: round(v / n, 4) for k, v in sorted(cnt.items(), key=lambda kv: (-kv[1], kv[0]))[:10]}


def _r(x, k=5):
    if x is None:
        return None
    if isinstance(x, (bool, np.bool_)):
        return int(x)
    if isinstance(x, (int, np.integer)):
        return int(x)
    x = float(x)
    return round(x, k) if math.isfinite(x) else None


def summarize_dp(res, X, nom):
    """выборки расчётной точки одного варианта -> числа"""
    N = len(res)
    ok = np.array([r["ok"] for r in res], dtype=int)
    feas = np.array([r["feas"] for r in res], dtype=int)
    sf = np.array([r.get("score_f", np.nan) if r["ok"] else np.nan for r in res], dtype=float)
    sw = np.array([r.get("score_w", np.nan) if r["ok"] else np.nan for r in res], dtype=float)
    eta = np.array([r.get("eta_el", np.nan) if r["ok"] else np.nan for r in res], dtype=float)
    y = np.where(feas == 1, np.nan_to_num(sf), 0.0)
    half = N // 2
    S_y = US.sobol1(X, y)
    S_f = US.sobol1(X, feas.astype(float))
    m = ok == 1
    S_s = US.sobol1(X[m], sf[m]) if m.sum() >= 16 else None
    bm = [r.get("burst_min") for r in res if r["ok"]]
    lm = [math.log10(r["life_min"]) for r in res if r["ok"] and r.get("life_min") and r["life_min"] > 0 and math.isfinite(r["life_min"])]
    k = int(feas.sum())
    return dict(
        n=N, n_ok=int(ok.sum()), n_feas=k, P_dp=round(k / N, 4), P_dp_ci=US.wilson(k, N),
        score_f=US.q(sf), score_w=US.q(sw), y=US.q(y), y_mean=_r(float(np.mean(y))), eta_el=US.q(eta),
        V_env_p50=US.q([r.get("V_env") for r in res if r["ok"]], (50,))[0], burst_min_p10=US.q(bm, (10,))[0],
        life_min_log10_p10=US.q(lm, (10,))[0],
        viol_freq=_freq([r["bad"] for r in res], N),
        half=dict(n=half, P_dp=round(float(feas[:half].mean()), 4), score_f=US.q(sf[:half]), y=US.q(y[:half])),
        sobol=dict(y=US.top(UF.NAMES, S_y), feas=US.top(UF.NAMES, S_f), score=US.top(UF.NAMES, S_s),
                   y_all=[_r(s, 4) for s in S_y] if S_y else None, feas_all=[_r(s, 4) for s in S_f] if S_f else None,
                   bins=US.bins_default(N)),
        nominal=dict(ok=nom["ok"], feas=nom["feas"], viol=nom["viol"], bad=nom["bad"], score_w=_r(nom.get("score_w")),
                     score_f=_r(nom.get("score_f")), trl=nom.get("trl"), eta_el=_r(nom.get("eta_el")), V_env=_r(nom.get("V_env"))),
    )


def od_qualifies(sm, rule, pdp_min, score_min):
    if rule == "none":
        return False
    if rule == "all":
        return sm["n_feas"] > 0
    p50 = sm["score_f"][1]
    return sm["P_dp"] >= pdp_min and p50 is not None and p50 >= score_min


def summarize_od(sm, odres, n_od):
    """нерасчётные режимы на первых n_od строках, где расчётная точка прошла"""
    run = [(s, r) for s, r in sorted(odres.items())]
    n_run = len(run)
    n_pass = sum(1 for s, r in run if r.get("od_ok"))
    p_od = n_pass / n_run if n_run else None
    joint = sum(1 for s, r in run if r.get("od_ok") and r.get("feas"))          # совместно на первых n_od строках
    sm_vals = [r.get("SM_min") for s, r in run if r.get("SM_min") is not None]
    out = dict(rows=n_od, n_run=n_run, n_pass=n_pass, P_od_dp=_r(p_od, 4), P_od_dp_ci=US.wilson(n_pass, n_run),
               P_joint_rows=_r(joint / n_od, 4) if n_od else None, od_fail_freq=_freq([r.get("od_bad") or [] for s, r in run], max(1, n_run)),
               SM_min=US.q(sm_vals), vigv_share=_r(sum(1 for s, r in run if r.get("vigv")) / n_run, 4) if n_run else None,
               score_f_od=US.q([r.get("score_f") for s, r in run if r.get("od_ok")]))
    return out


def finalize(cand, sm, od):
    P = None
    status = "не считались (вариант не прошёл порог)"
    if od is not None:
        if od["n_run"]:
            P = round(sm["P_dp"] * od["P_od_dp"], 4)
            status = "посчитаны"
        else:
            status = "нет строк с прошедшей расчётной точкой среди первых n_od"
    p50 = sm["score_f"][1]
    rob = round(p50 * P, 4) if (P is not None and p50 is not None) else None
    rec = dict(cand.get("rec") or {})
    return dict(src=cand["src"], j_src=cand.get("j_src"), key=cand.get("key"), idx=cand.get("idx"), kv=cand.get("kv"),
                niche=cand.get("niche"), arch=cand.get("arch"), vigv=cand.get("vigv"), heat=cand.get("heat"), rec=rec,
                dp=sm, od=od, od_status=status, P_feas=P, robust=rob,
                P_feas_upper=sm["P_dp"] if P is None else P)


# =========================================================================== части
def part_path(out, kv, K, shard, shards, tag="", gz=False):
    return os.path.join(out, f"{PART_PREFIX}{tag}_{kv}_c{K}_s{shard:03d}of{shards:03d}.json" + (".gz" if gz else ""))


def run_parts(cands, meta_in, shard_list, shards, out, n=DEFAULT_N, n_od=DEFAULT_N_OD, seed=UF.UQ_SEED, od_rule="auto",
              pdp_min=OD_PDP_MIN, score_min=SCORE_MIN - OD_SCORE_MARGIN, workers=1, gz=False, tag=""):
    t0 = time.time()
    K = len(cands)
    kv = meta_in["kv"]
    tgt = FS.target(kv)
    c0 = cands[0]["concept"] if cands else None
    if c0:
        L1._prep(int(c0.get("n") or 2000), int(c0.get("seed") or FS.SEED), kv)
    u, X, vals = UF.sample(n, seed)
    n_od = min(n_od, n)
    idx_by = {k: FS.shard_indices(K, k, shards) for k in shard_list}
    ids = sorted(i for v in idx_by.values() for i in v)
    t1 = time.time()
    dp, nom, cpu_dp = run_dp(cands, ids, vals, tgt, workers)
    t_dp = time.time() - t1
    sms = {i: summarize_dp(dp[i], X, nom[i]) for i in ids}
    qual = [i for i in ids if od_qualifies(sms[i], od_rule, pdp_min, score_min)]
    tasks = [(i, s) for i in qual for s in range(n_od) if dp[i][s]["feas"]]
    t2 = time.time()
    odr, cpu_od = run_od(tasks, cands, vals, tgt, workers)
    t_od = time.time() - t2
    rows = {}
    for i in ids:
        od = None
        if i in qual:
            od = summarize_od(sms[i], {s: odr[(i, s)] for s in range(n_od) if (i, s) in odr}, n_od)
        rows[i] = dict(u=i, **finalize(cands[i], sms[i], od))
    el = time.time() - t0
    n_ev = len(ids) * (n + 1)
    meta = dict(level=LEVEL, date=DATE, software=FS.software(), method=__doc__.split("\n")[0], kv=kv, K=K, shards=shards, tag=tag,
                input=meta_in, cand_sha=cand_sha(cands), code_sha=code_sha(), l1_code_sha=L1.code_sha(), factors_sha=UF.table_sha(),
                n=n, n_od=n_od, seed=seed, sample_fingerprint=UF.fingerprint(n, seed), od_rule=od_rule, od_pdp_min=pdp_min,
                od_score_min=score_min, score_key="score_f: взвешенная × (0,5 + 0,5·TRL/9)", workers=workers,
                call=dict(shards=shard_list, variants=len(ids), od_variants=len(qual), od_runs=len(tasks), elapsed_s=round(el, 2),
                          dp_s=round(t_dp, 2), od_s=round(t_od, 2), cpu_dp_s=round(cpu_dp, 2), cpu_od_s=round(cpu_od, 2),
                          cpu_s_per_dp_eval=round(cpu_dp / max(1, n_ev), 5), cpu_s_per_od_eval=round(cpu_od / max(1, len(tasks)), 3) if tasks else None,
                          host_cpus=os.cpu_count(), started=datetime.datetime.fromtimestamp(t0, datetime.timezone.utc).isoformat()))
    os.makedirs(out, exist_ok=True)
    paths = []
    for k in shard_list:
        doc = dict(meta, shard=k, n_rows=len(idx_by[k]), rows=[rows[i] for i in idx_by[k]])
        p = part_path(out, kv, K, k, shards, tag, gz)
        txt = json.dumps(doc, ensure_ascii=False, separators=(",", ":"), default=_r)
        with (gzip.open(p, "wt", encoding="utf-8") if gz else open(p, "w", encoding="utf-8")) as fh:
            fh.write(txt)
        paths.append(p)
    return paths, meta, rows


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--l1-in", nargs="*", default=None, help="папки с частями funnel_l1 (с разборами)")
    ap.add_argument("--l1-kv", default="KV-3")
    ap.add_argument("--l1-m", type=int, default=None)
    ap.add_argument("--l1-shards", type=int, default=None)
    ap.add_argument("--l1-tag", default="")
    ap.add_argument("--candidates", default=None, help="файл кандидатов merge_pareto.py (вместо частей уровня 1)")
    ap.add_argument("--top", type=int, default=0, help="первых кандидатов (0 - все)")
    ap.add_argument("--n", type=int, default=DEFAULT_N, help="выборок латинского гиперкуба на вариант")
    ap.add_argument("--n-od", type=int, default=DEFAULT_N_OD, help="выборок с нерасчётными режимами (первые строки гиперкуба)")
    ap.add_argument("--seed", type=int, default=UF.UQ_SEED)
    ap.add_argument("--od", default="auto", choices=["auto", "none", "all"])
    ap.add_argument("--od-pdp-min", type=float, default=OD_PDP_MIN)
    ap.add_argument("--od-score-min", type=float, default=SCORE_MIN - OD_SCORE_MARGIN)
    ap.add_argument("--shards", type=int, required=True)
    ap.add_argument("--shard", default="all")
    ap.add_argument("--workers", type=int, default=0)
    ap.add_argument("--out", default=os.environ.get("GTU_OUT", "."))
    ap.add_argument("--gzip", action="store_true")
    ap.add_argument("--tag", default="")
    ap.add_argument("--local", action="store_true", help=f"локальная миниатюра: не больше {LOCAL_MAX}")
    a = ap.parse_args(argv)
    if a.candidates:
        cands, meta_in = load_candidates(a.candidates)
    elif a.l1_in:
        if not (a.l1_m and a.l1_shards):
            ap.error("--l1-in требует --l1-m и --l1-shards")
        cands, meta_in = load_l1(a.l1_in, a.l1_kv, a.l1_m, a.l1_shards, a.l1_tag)
    else:
        ap.error("нужен --l1-in или --candidates")
    if a.top:
        cands = cands[:a.top]
        meta_in["top"] = a.top
    if a.local and (len(cands) > LOCAL_MAX["variants"] or a.n > LOCAL_MAX["n"] or (a.od != "none" and a.n_od > LOCAL_MAX["n_od"])):
        ap.error(f"локально - не больше {LOCAL_MAX}")
    sl = L0.parse_shards(a.shard, a.shards)
    w = a.workers or L0.default_workers()
    paths, meta, rows = run_parts(cands, meta_in, sl, a.shards, a.out, a.n, a.n_od, a.seed, a.od, a.od_pdp_min, a.od_score_min, w, a.gzip,
                                  a.tag)
    c = meta["call"]
    print(f"{meta['kv']} UQ: K={meta['K']}, части {sl[0]}..{sl[-1]} ({len(sl)} из {a.shards}), {c['variants']} вариантов × {a.n} выборок, "
          f"нерасчётные - {c['od_variants']} вариантов, {c['od_runs']} расчётов; {c['elapsed_s']} с, {w} процессов; "
          f"процессор: {c['cpu_s_per_dp_eval']} с на расчётную точку, {c['cpu_s_per_od_eval']} с на нерасчётные; файлов {len(paths)} -> {a.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
