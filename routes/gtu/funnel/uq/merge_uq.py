"""Устойчивый отбор, часть 2: слияние частей uq_l1 с контролем полноты, «50 устойчивых разных» и сравнение с обычным отбором.
Всё - оценка.

  1. Слияние: все части 0..S-1, у каждой - ровно свои варианты (u % S), одинаковые вход (отпечаток кандидатов), код, таблица
     факторов, выборка (N, зерно, отпечаток матрицы), правило нерасчётных режимов. Повтор части с тем же содержанием (прерывание
     спота) - допустим, с другим - ошибка.
  2. Отборы (жадный, как merge_l1.greedy: штраф LAMBDA·exp(-d/D0) за близость к выбранным, не больше NICHE_MAX на нишу, расстояние
     merge_l1.distance):
       A - обычный: номинальная оценка уровня 1 (взвешенная), допустим на уровне 1, оценка >= SCORE_MIN - как merge_l1;
       B - обычный с множителем готовности (решение владельца 15.11.2026): номинальная оценка × (0,5 + 0,5·TRL/9) >= SCORE_MIN;
       C - устойчивый: значение P50(оценки) × P(допустим), допуск - нерасчётные режимы посчитаны, P(допустим) >= P_MIN и
           P50 оценки >= SCORE_MIN.
     Сравнение C с B (и A): сколько вариантов сменилось, у каждого вытесненного и нового - почему (P(допустим), чаще всего
     нарушаемое правило, главный фактор по Соболю, P50 против номинала).
  3. Сводки: распределение P(допустим), правила, которые чаще всего ломаются под неопределённостью, средние индексы Соболя по
     вариантам (какие неопределённости важнее всего снять расчётом или испытанием), сходимость (N/2 против N).

  python merge_uq.py --kv KV-3 --k 3954 --shards 50 --in routes/gtu/platform/aws/runs/<run_id> --out-dir aws/uq_l1 --suffix _KV-3_l1
"""
import argparse, csv, glob, gzip, hashlib, json, math, os, re, sys, time
from collections import Counter

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
FUNNEL = os.path.dirname(HERE)
for p in (FUNNEL, HERE):
    if p not in sys.path:
        sys.path.insert(0, p)
import merge_l1 as M1            # noqa: E402
import uq_factors as UF          # noqa: E402
import uq_l1 as U                # noqa: E402

P_MIN = 0.5                   # устойчиво допустим: не меньше половины реализаций неопределённостей (оценка; чувствительность - 0,25 и 0,75)
P_MIN_SENS = (0.25, 0.5, 0.75)
META_KEYS = ("level", "kv", "K", "shards", "tag", "cand_sha", "code_sha", "factors_sha", "n", "n_od", "seed", "sample_fingerprint",
             "od_rule", "od_pdp_min", "od_score_min")


class MergeError(Exception):
    pass


def find_parts(dirs, kv, K, shards, tag=""):
    pat = re.compile(rf"^{U.PART_PREFIX}{re.escape(tag)}_{re.escape(kv)}_c{K}_s(\d{{3}})of{shards:03d}\.json(\.gz)?$")
    out = []
    for d in dirs:
        for p in glob.glob(os.path.join(d, "**", f"{U.PART_PREFIX}{tag}_{kv}_c{K}_s*of{shards:03d}.json*"), recursive=True):
            m = pat.match(os.path.basename(p))
            if m:
                out.append((int(m.group(1)), p))
    return sorted(out)


def load_part(p):
    with (gzip.open(p, "rt", encoding="utf-8") if p.endswith(".gz") else open(p, encoding="utf-8")) as fh:
        return json.load(fh)


def content_hash(d):
    return hashlib.sha256(json.dumps(d["rows"], ensure_ascii=False, sort_keys=True).encode()).hexdigest()


def merge_parts(found, kv, K, shards):
    checks = dict(files=len(found), shards_expected=shards, missing=[], duplicate_same=[], duplicate_conflict=[], meta_mismatch=[],
                  index_errors=[], unreadable=[])
    ref, by = None, {}
    for k, p in found:
        try:
            d = load_part(p)
        except (OSError, EOFError, ValueError) as e:
            checks["unreadable"].append(dict(file=os.path.basename(p), error=type(e).__name__))
            continue
        bad = {}
        if d.get("shard") != k or d.get("kv") != kv or d.get("K") != K or d.get("shards") != shards or d.get("level") != U.LEVEL:
            bad["id"] = [d.get("shard"), d.get("kv"), d.get("K"), d.get("shards")]
        if ref is None:
            ref = {x: v for x, v in d.items() if x != "rows"}
        else:
            for x in META_KEYS:
                if d.get(x) != ref.get(x):
                    bad[x] = [d.get(x), ref.get(x)]
        if bad:
            checks["meta_mismatch"].append(dict(file=os.path.basename(p), diff=bad))
            continue
        exp = U.FS.shard_indices(K, k, shards)
        got = [r["u"] for r in d["rows"]]
        if got != exp:
            checks["index_errors"].append(dict(file=os.path.basename(p), shard=k, lost=sorted(set(exp) - set(got))[:20],
                                               foreign=sorted(set(got) - set(exp))[:20]))
            continue
        h = content_hash(d)
        if k in by:
            (checks["duplicate_same"] if by[k]["hash"] == h else checks["duplicate_conflict"]).append(dict(shard=k, file=os.path.basename(p)))
            continue
        by[k] = dict(rows=d["rows"], hash=h, call=d.get("call"))
    checks["missing"] = [k for k in range(shards) if k not in by]
    errs = {x: checks[x] for x in ("missing", "duplicate_conflict", "meta_mismatch", "index_errors") if checks[x]}
    if errs or ref is None:
        raise MergeError(json.dumps(errs or {"files": "частей не найдено"}, ensure_ascii=False))
    rows = sorted((r for v in by.values() for r in v["rows"]), key=lambda r: r["u"])
    if [r["u"] for r in rows] != list(range(K)):
        raise MergeError(f"номера после слияния не равны 0..{K - 1}")
    checks.update(indices_ok=True, n_rows=len(rows))
    calls = [dict(v["call"] or {}, shard=k) for k, v in sorted(by.items())]
    return rows, ref, checks, calls


# =========================================================================== записи для отбора
def _base(r):
    """поля записи для расстояния merge_l1.features: дискретные оси и непрерывные (t41, π, ε, КПД, объём, масса)"""
    rec = dict(r.get("rec") or {})
    rec["j"] = r["u"]
    rec["niche"] = r.get("niche") or rec.get("niche")
    rec["arch"] = r.get("arch") or rec.get("arch")
    return rec


def std_records(rows, with_trl=False):
    out = []
    for r in rows:
        rec = _base(r)
        s, t = rec.get("score"), rec.get("trl_min")
        if with_trl:
            s = s * U.trl_factor(t) if (s is not None and t is not None) else None
        out.append(dict(rec, score=s, feasible=int(bool(rec.get("feasible")))))
    return out


def robust_eligible(r, p_min=P_MIN, score_min=U.SCORE_MIN):
    p50 = r["dp"]["score_f"][1]
    return r.get("P_feas") is not None and r["P_feas"] >= p_min and p50 is not None and p50 >= score_min


def robust_records(rows, p_min=P_MIN, score_min=U.SCORE_MIN):
    """для merge_l1.greedy: score := P50 × P(допустим), feasible := допуск устойчивого отбора, порог 0"""
    out = []
    for r in rows:
        rec = _base(r)
        out.append(dict(rec, score=r.get("robust"), feasible=int(robust_eligible(r, p_min, score_min))))
    return out


def select(recs, n=M1.TOP_N, score_min=None):
    return M1.greedy(recs, n=n, score_min=M1.SCORE_MIN if score_min is None else score_min)


# =========================================================================== почему сменился
def _top_rule(r):
    vf = (r.get("od") or {}).get("od_fail_freq") or {}
    dpf = r["dp"].get("viol_freq") or {}
    items = sorted(list(dpf.items()) + [(k, v * (r["dp"]["P_dp"])) for k, v in vf.items()], key=lambda kv: -kv[1])
    return items[0] if items else None


def _top_factor(r):
    s = (r["dp"].get("sobol") or {}).get("y") or (r["dp"].get("sobol") or {}).get("feas")
    return s[0] if s else None


def why_out(r, sel_vals, p_min, score_min):
    """почему вариант обычного отбора не вошёл в устойчивый"""
    p50 = r["dp"]["score_f"][1]
    nom = r["dp"]["nominal"].get("score_f")
    tr, tf = _top_rule(r), _top_factor(r)
    tail = (f"; чаще всего нарушается {tr[0]} (в {tr[1] * 100:.0f} % выборок)" if tr else "") + \
           (f"; главный фактор {tf[0]} (S1 = {tf[1]:.2f})" if tf else "")
    if r.get("P_feas") is None:
        return "нерасчётные режимы под неопределённостью не считались: " + (r.get("od_status") or "") + \
               f"; P(расчётная точка) = {r['dp']['P_dp']:.2f}" + tail
    if r["P_feas"] < p_min:
        return f"P(допустим) = {r['P_feas']:.2f} < {p_min}" + tail
    if p50 is None or p50 < score_min:
        return f"P50 оценки {p50} < {score_min} (номинал {nom})" + tail
    return (f"устойчивое значение P50 × P = {r['robust']:.3f} (P = {r['P_feas']:.2f}) - вытеснен более устойчивыми, штрафом близости "
            f"или лимитом ниши (наименьшее значение среди отобранных {min(sel_vals):.3f})" if sel_vals else "вытеснен")


def why_in(r, std_sel_scores):
    rec = r.get("rec") or {}
    nom_f = rec.get("feasible")
    s = rec.get("score")
    base = f"P(допустим) = {r['P_feas']:.2f}, P50 оценки {r['dp']['score_f'][1]:.3f}, устойчивое значение {r['robust']:.3f}"
    if not nom_f:
        return base + f"; на уровне 1 номинально недопустим ({rec.get('reject')}), но под неопределённостью допустим чаще порога"
    if std_sel_scores and s is not None and s < min(std_sel_scores):
        return base + f"; номинальная оценка {s:.3f} ниже наименьшей в обычном отборе ({min(std_sel_scores):.3f}) - вошёл за счёт устойчивости"
    return base + "; в обычный отбор не прошёл по штрафу близости или лимиту ниши - освободилось место"


def compare(rows, selA, selB, selC, p_min, score_min):
    by = {r["u"]: r for r in rows}
    jA, jB, jC = ({s["j"] for s in x} for x in (selA, selB, selC))
    vals = [s["score"] for s in selC]
    sB = [s["score"] for s in selB]
    out_B = sorted(jB - jC)
    in_B = sorted(jC - jB)
    return dict(
        n=dict(A=len(selA), B=len(selB), C=len(selC)),
        overlap=dict(A_C=len(jA & jC), B_C=len(jB & jC), A_B=len(jA & jB)),
        changed_vs_B=dict(out=len(out_B), new=len(in_B)),
        changed_vs_A=dict(out=len(jA - jC), new=len(jC - jA)),
        out_vs_B=[dict(u=u, j_src=by[u].get("j_src"), key=by[u].get("key"), arch=by[u]["arch"], niche=by[u]["niche"], P_feas=by[u].get("P_feas"),
                       P_dp=by[u]["dp"]["P_dp"], why=why_out(by[u], vals, p_min, score_min)) for u in out_B],
        new_vs_B=[dict(u=u, j_src=by[u].get("j_src"), key=by[u].get("key"), arch=by[u]["arch"], niche=by[u]["niche"], P_feas=by[u].get("P_feas"),
                       why=why_in(by[u], sB)) for u in in_B],
        archetypes=dict(A=dict(Counter(by[j]["arch"] for j in jA)), B=dict(Counter(by[j]["arch"] for j in jB)),
                        C=dict(Counter(by[j]["arch"] for j in jC))),
        P_feas_of_B=dict(mean=_m([by[j].get("P_feas") for j in jB]), min=_mn([by[j].get("P_feas") for j in jB])),
        P_feas_of_C=dict(mean=_m([by[j].get("P_feas") for j in jC]), min=_mn([by[j].get("P_feas") for j in jC])),
        note="«сумма прячет потерю»: у каждого вытесненного - своя причина; P(допустим) обычного отбора (B) - сколько из него ломается "
             "под неопределённостью")


def _m(v):
    v = [x for x in v if x is not None]
    return round(float(np.mean(v)), 4) if v else None


def _mn(v):
    v = [x for x in v if x is not None]
    return round(float(min(v)), 4) if v else None


# =========================================================================== сводки
def sobol_summary(rows, key="y_all", subset=None):
    S = [r["dp"]["sobol"].get(key) for r in rows if (subset is None or r["u"] in subset) and r["dp"]["sobol"].get(key)]
    if not S:
        return None
    A = np.array(S, dtype=float)
    top1 = Counter(UF.NAMES[int(np.argmax(s))] for s in A)
    order = np.argsort(-A.mean(axis=0))
    return dict(n_variants=len(S), mean=[[UF.NAMES[i], round(float(A[:, i].mean()), 4), round(float(np.percentile(A[:, i], 90)), 4)]
                                         for i in order],
                top1_share={k: round(v / len(S), 4) for k, v in top1.most_common()},
                note="средний по вариантам индекс первого порядка и его P90; доля вариантов, где фактор - главный")


def rule_summary(rows):
    acc, n = Counter(), 0
    for r in rows:
        n += 1
        for k, v in (r["dp"].get("viol_freq") or {}).items():
            acc[re.sub(r"\d+$", "", k).rstrip("_")] += v
    return {k: round(v / max(1, n), 4) for k, v in acc.most_common(15)}


def p_hist(rows, key="P_feas_upper"):
    edges = [0.0, 0.01, 0.1, 0.25, 0.5, 0.75, 0.9, 1.0001]
    v = [r.get(key) for r in rows if r.get(key) is not None]
    h = np.histogram(v, bins=edges)[0].tolist() if v else []
    return dict(edges=edges[:-1] + [1.0], counts=h, n=len(v))


def convergence(rows):
    """те же числа на первых N/2 строках гиперкуба против всех N"""
    d = dict(P_dp=[], p10=[], p50=[], p90=[])
    for r in rows:
        h, f = r["dp"]["half"], r["dp"]
        d["P_dp"].append(abs(h["P_dp"] - f["P_dp"]))
        for k, i in (("p10", 0), ("p50", 1), ("p90", 2)):
            if h["score_f"][i] is not None and f["score_f"][i] is not None:
                d[k].append(abs(h["score_f"][i] - f["score_f"][i]))
    return {k: dict(mean=round(float(np.mean(v)), 4), p95=round(float(np.percentile(v, 95)), 4), max=round(float(max(v)), 4))
            if v else None for k, v in d.items()} | dict(rule="|N/2 − N| по вариантам: P(расчётная точка) и P10/P50/P90 оценки")


def by_arch(rows):
    out = {}
    for a in sorted({r["arch"] for r in rows}):
        rr = [r for r in rows if r["arch"] == a]
        P = [r.get("P_feas_upper") for r in rr]
        rob = [r["robust"] for r in rr if r.get("robust") is not None]
        out[a] = dict(n=len(rr), P_upper_p50=round(float(np.median(P)), 4), P_upper_max=round(float(max(P)), 4),
                      od_done=sum(1 for r in rr if r.get("P_feas") is not None), robust_max=round(max(rob), 4) if rob else None,
                      eligible=sum(1 for r in rr if robust_eligible(r)))
    return out


def entry(r, k, sel):
    rec = r.get("rec") or {}
    return dict(rank=k + 1, u=r["u"], src=r["src"], j_src=r.get("j_src"), key=r.get("key"), idx=r.get("idx"), arch=r["arch"], niche=r["niche"],
                robust=r["robust"], P_feas=r["P_feas"], P_dp=r["dp"]["P_dp"], P_dp_ci=r["dp"]["P_dp_ci"],
                P_od_dp=(r.get("od") or {}).get("P_od_dp"), P_od_dp_ci=(r.get("od") or {}).get("P_od_dp_ci"),
                score_f_p10_p50_p90=r["dp"]["score_f"], score_nominal=r["dp"]["nominal"].get("score_f"),
                eta_el_p10_p50_p90=r["dp"]["eta_el"], viol_freq=r["dp"]["viol_freq"], od_fail_freq=(r.get("od") or {}).get("od_fail_freq"),
                SM_min_p10_p50_p90=(r.get("od") or {}).get("SM_min"), sobol_top5=r["dp"]["sobol"].get("y"),
                l1=dict(score=rec.get("score"), feasible=rec.get("feasible"), eta_th_shaft=rec.get("eta_th_shaft"), V_env=rec.get("V_env"),
                        trl_min=rec.get("trl_min"), n_parts=rec.get("n_parts"), steps=rec.get("steps")),
                adj=round(sel.get("_adj"), 4) if sel.get("_adj") is not None else None, penalty=round(sel.get("_pen", 0.0), 4))


CSV_COLS = ["u", "src", "j_src", "key", "idx", "arch", "niche", "l1_feasible", "l1_score", "nominal_feas_dp", "nominal_score_f", "P_dp",
            "P_od_dp", "P_feas", "P_feas_upper", "od_status", "score_f_p10", "score_f_p50", "score_f_p90", "robust", "eta_el_p10", "eta_el_p50",
            "eta_el_p90", "top_rule", "top_factor", "robust_rank", "std_rank_B"]


def write_csv(path, rows, rank_c, rank_b):
    with open(path, "w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(CSV_COLS)
        for r in rows:
            rec, dp = r.get("rec") or {}, r["dp"]
            tr, tf = _top_rule(r), _top_factor(r)
            w.writerow([r["u"], r["src"], r.get("j_src"), r.get("key"), r.get("idx"), r["arch"], r["niche"], rec.get("feasible"), rec.get("score"),
                        dp["nominal"]["feas"], dp["nominal"].get("score_f"), dp["P_dp"], (r.get("od") or {}).get("P_od_dp"), r.get("P_feas"),
                        r.get("P_feas_upper"), r.get("od_status"), *dp["score_f"], r.get("robust"), *dp["eta_el"],
                        f"{tr[0]}:{tr[1]:.2f}" if tr else "", f"{tf[0]}:{tf[1]:.2f}" if tf else "", rank_c.get(r["u"]), rank_b.get(r["u"])])


def run_merge(found, kv, K, shards, out_dir, sfx="", p_min=P_MIN, score_min=U.SCORE_MIN):
    t0 = time.time()
    rows, ref, checks, calls = merge_parts(found, kv, K, shards)
    A = select(std_records(rows))
    B = select(std_records(rows, with_trl=True))
    C = select(robust_records(rows, p_min, score_min), score_min=0.0)
    cmp_ = compare(rows, A, B, C, p_min, score_min)
    sens = []
    for pm in P_MIN_SENS:
        Cp = select(robust_records(rows, pm, score_min), score_min=0.0)
        sens.append(dict(P_min=pm, eligible=sum(1 for r in rows if robust_eligible(r, pm, score_min)), selected=len(Cp),
                         overlap_with_C=len({s["j"] for s in Cp} & {s["j"] for s in C}),
                         overlap_with_B=len({s["j"] for s in Cp} & {s["j"] for s in B})))
    by = {r["u"]: r for r in rows}
    top = [entry(by[s["j"]], k, s) for k, s in enumerate(C)]
    rank_c = {s["j"]: k + 1 for k, s in enumerate(C)}
    rank_b = {s["j"]: k + 1 for k, s in enumerate(B)}
    eligible_ids = {r["u"] for r in rows if robust_eligible(r, p_min, score_min)}
    os.makedirs(out_dir, exist_ok=True)
    p_csv = os.path.join(out_dir, f"uq{sfx}.csv")
    write_csv(p_csv, rows, rank_c, rank_b)
    cpu_dp = sum((c.get("cpu_dp_s") or 0) for c in calls)
    cpu_od = sum((c.get("cpu_od_s") or 0) for c in calls)
    n_od = sum((c.get("od_runs") or 0) for c in calls)
    soft = dict(software=ref.get("software"), date=U.DATE)
    res = dict(title=f"Устойчивый отбор: Монте-Карло по неопределённостям уровня 1 - {kv}, {K} вариантов", **soft,
               method=(U.__doc__.split("\n\n")[0] + " " + __doc__.split("\n")[0]), kv=kv, K=K, n=ref["n"], n_od=ref["n_od"], seed=ref["seed"],
               sample_fingerprint=ref["sample_fingerprint"], input=ref.get("input"), cand_sha=ref["cand_sha"], code_sha=ref["code_sha"],
               l1_code_sha=ref.get("l1_code_sha"), factors=UF.table(), factors_sha=ref["factors_sha"],
               od_rule=dict(rule=ref["od_rule"], P_dp_min=ref["od_pdp_min"], score_min=ref["od_score_min"]),
               merge_checks=checks,
               timing=dict(cpu_dp_s=round(cpu_dp, 1), cpu_od_s=round(cpu_od, 1), od_runs=n_od,
                           cpu_s_per_dp_eval=round(cpu_dp / max(1, K * (ref["n"] + 1)), 5), cpu_s_per_od_eval=round(cpu_od / n_od, 3) if n_od else None,
                           calls=calls),
               summary=dict(
                   variants=K, nominal_dp_feasible=sum(1 for r in rows if r["dp"]["nominal"]["feas"]),
                   l1_feasible=sum(1 for r in rows if (r.get("rec") or {}).get("feasible")),
                   nominal_vs_l1_score_max_abs=_m([abs(r["dp"]["nominal"]["score_w"] - r["rec"]["score"]) for r in rows
                                                   if r["dp"]["nominal"].get("score_w") is not None and (r.get("rec") or {}).get("score") is not None]),
                   od_done=sum(1 for r in rows if r.get("P_feas") is not None), robust_eligible=len(eligible_ids),
                   P_hist_upper=p_hist(rows), P_hist_od=p_hist([r for r in rows if r.get("P_feas") is not None], "P_feas"),
                   rules_under_uncertainty=rule_summary(rows), by_arch=by_arch(rows)),
               sobol=dict(y=sobol_summary(rows), feas=sobol_summary(rows, "feas_all"),
                          y_eligible=sobol_summary(rows, subset=eligible_ids),
                          estimator="по готовой выборке: корреляционное отношение по слоям x_i (M = √N/1,5) с поправкой на смещение (Plischke - "
                                    "Borgonovo - Smith 2013, по памяти); y = оценка × [допустим в расчётной точке]"),
               convergence=convergence(rows),
               selection=dict(rule=dict(top_n=M1.TOP_N, niche_max=M1.NICHE_MAX, lam=M1.LAMBDA, d0=M1.D0, P_min=p_min, score_min=score_min,
                                        A="merge_l1: номинальная оценка уровня 1, допустим, >= SCORE_MIN",
                                        B="номинальная оценка × (0,5 + 0,5·TRL/9) >= SCORE_MIN (решение владельца 15.11.2026)",
                                        C="P50(оценки × множитель TRL) × P(допустим); допуск: нерасчётные посчитаны, P >= P_min, P50 >= SCORE_MIN"),
                              compare=cmp_, P_min_sensitivity=sens),
               files=dict(csv=os.path.basename(p_csv), robust50=f"uq_robust50{sfx}.json"))
    res["timing"]["merge_s"] = round(time.time() - t0, 1)
    p_res = os.path.join(out_dir, f"uq_results{sfx}.json")
    with open(p_res, "w", encoding="utf-8") as fh:
        json.dump(res, fh, ensure_ascii=False, indent=1, default=U._r)
    p_top = os.path.join(out_dir, f"uq_robust50{sfx}.json")
    with open(p_top, "w", encoding="utf-8") as fh:
        json.dump(dict(title=f"Устойчивый отбор: {len(top)} устойчивых разных - {kv}", **soft, kv=kv, rule=res["selection"]["rule"],
                       selected=top, standard_B=[dict(rank=k + 1, u=s["j"], arch=s["arch"], niche=s["niche"], score_trl=round(s["score"], 4),
                                                      P_feas=by[s["j"]].get("P_feas"), P_dp=by[s["j"]]["dp"]["P_dp"]) for k, s in enumerate(B)],
                       note="все числа - оценка уровня 1 под неопределённостями uq_factors (README узла uq)"),
                  fh, ensure_ascii=False, indent=1, default=U._r)
    return dict(results=p_res, top=p_top, csv=p_csv, res=res, rows=rows, A=A, B=B, C=C)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--kv", required=True)
    ap.add_argument("--k", type=int, required=True, help="число кандидатов K (имя частей _c<K>_)")
    ap.add_argument("--shards", type=int, required=True)
    ap.add_argument("--tag", default="")
    ap.add_argument("--in", dest="dirs", nargs="+", required=True)
    ap.add_argument("--out-dir", default=HERE)
    ap.add_argument("--suffix", default="")
    ap.add_argument("--p-min", type=float, default=P_MIN)
    a = ap.parse_args(argv)
    found = find_parts(a.dirs, a.kv, a.k, a.shards, a.tag)
    try:
        out = run_merge(found, a.kv, a.k, a.shards, a.out_dir, a.suffix, a.p_min)
    except MergeError as e:
        print("слияние отклонено:", e)
        return 2
    s, c = out["res"]["summary"], out["res"]["selection"]["compare"]
    print(f"{a.kv}: {s['variants']} вариантов, нерасчётные под неопределённостью - {s['od_done']}, устойчиво допустимых (P >= {a.p_min}) - "
          f"{s['robust_eligible']}; отбор: обычный B {c['n']['B']}, устойчивый C {c['n']['C']}, сменилось {c['changed_vs_B']['out']} "
          f"-> {out['results']}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
