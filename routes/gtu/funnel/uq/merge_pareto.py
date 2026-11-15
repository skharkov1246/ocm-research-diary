"""Фронты Парето ниш: слияние частей pareto_niche с контролем полноты, сводки и файл кандидатов для устойчивого отбора. Всё - оценка.

  1. Слияние: все части 0..S-1, у каждой - ровно свои ниши (k % S), одинаковые вход (отбор уровня 0, sha256), код, параметры NSGA-II
     (популяция, поколения, проверка, режим устойчивости); повтор части с тем же содержанием допустим, с другим - ошибка.
  2. Сводки: у скольких ниш есть допустимый фронт, размеры фронтов, гиперобъём (нормированные цели, точка отсчёта 1) и его рост со
     второй половины поколений (сходимость), сколько представителей прошли нерасчётные режимы, сравнение с вариантами уровня 1
     (доминирует ли проверенный фронт ниши их точки, прирост лучшей оценки), общий фронт по всем нишам (недоминируемые среди
     проверенных допустимых) по архетипам.
  3. Кандидаты устойчивого отбора: проверенные допустимые члены фронтов - компактно (ниша, опорная концепция, гены): uq_l1 строит
     концепцию заново (pareto_niche.Niche.build) - файл остаётся малым и для KV-1 (1864 ниши).

  python merge_pareto.py --kv KV-3 --nn 401 --shards 50 --in routes/gtu/platform/aws/runs/<run_id> --out-dir aws/pareto --suffix _KV-3
"""
import argparse, csv, glob, gzip, hashlib, json, os, re, sys, time
from collections import Counter

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
FUNNEL = os.path.dirname(HERE)
for p in (FUNNEL, HERE):
    if p not in sys.path:
        sys.path.insert(0, p)
import nsga2 as NS               # noqa: E402
import pareto_niche as PN        # noqa: E402

META_KEYS = ("level", "kv", "NN", "shards", "tag", "code_sha", "pop", "gens", "check", "robust", "p_target")


class MergeError(Exception):
    pass


def find_parts(dirs, kv, NN, shards, tag=""):
    pat = re.compile(rf"^{PN.PART_PREFIX}{re.escape(tag)}_{re.escape(kv)}_n{NN}_s(\d{{3}})of{shards:03d}\.json(\.gz)?$")
    out = []
    for d in dirs:
        for p in glob.glob(os.path.join(d, "**", f"{PN.PART_PREFIX}{tag}_{kv}_n{NN}_s*of{shards:03d}.json*"), recursive=True):
            m = pat.match(os.path.basename(p))
            if m:
                out.append((int(m.group(1)), p))
    return sorted(out)


def load_part(p):
    with (gzip.open(p, "rt", encoding="utf-8") if p.endswith(".gz") else open(p, encoding="utf-8")) as fh:
        return json.load(fh)


def content_hash(d):
    rows = [{k: v for k, v in r.items() if k not in ("t_s", "t_cpu_s")} for r in d["rows"]]
    return hashlib.sha256(json.dumps(rows, ensure_ascii=False, sort_keys=True).encode()).hexdigest()


def merge_parts(found, kv, NN, shards):
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
        if d.get("shard") != k or d.get("kv") != kv or d.get("NN") != NN or d.get("shards") != shards or d.get("level") != PN.LEVEL:
            bad["id"] = [d.get("shard"), d.get("kv"), d.get("NN"), d.get("shards")]
        if ref is None:
            ref = {x: v for x, v in d.items() if x != "rows"}
        else:
            for x in META_KEYS:
                if d.get(x) != ref.get(x):
                    bad[x] = [d.get(x), ref.get(x)]
            if (d.get("input") or {}).get("sha256") != (ref.get("input") or {}).get("sha256"):
                bad["input"] = "другой отбор уровня 0"
        if bad:
            checks["meta_mismatch"].append(dict(file=os.path.basename(p), diff=bad))
            continue
        exp = PN.FS.shard_indices(NN, k, shards)
        got = [r["k"] for r in d["rows"]]
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
    rows = sorted((r for v in by.values() for r in v["rows"]), key=lambda r: r["k"])
    if [r["k"] for r in rows] != list(range(NN)):
        raise MergeError(f"номера ниш после слияния не равны 0..{NN - 1}")
    checks.update(indices_ok=True, n_rows=len(rows))
    return rows, ref, checks, [dict(v["call"] or {}, shard=k) for k, v in sorted(by.items())]


def _st(v):
    v = [x for x in v if x is not None]
    if not v:
        return None
    return dict(n=len(v), mean=round(float(np.mean(v)), 4), p50=round(float(np.median(v)), 4), max=round(float(max(v)), 4))


def summarize(rows):
    with_front = [r for r in rows if r["front_size"]]
    conv = [r["hv"] - (r.get("half") or {}).get("hv", 0.0) for r in with_front if r.get("half")]
    cmp_ = [r["l1_compare"] for r in rows if r.get("l1_compare")]
    by = {}
    for a in sorted({r["arch"] for r in rows}):
        rr = [r for r in rows if r["arch"] == a]
        by[a] = dict(niches=len(rr), with_front=sum(1 for r in rr if r["front_size"]), checked_feasible=sum(r["checked_feasible"] for r in rr),
                     hv_max=round(max(r["hv"] for r in rr), 4), best_score_w=max([r["best_score_w"] for r in rr if r["best_score_w"] is not None],
                                                                                  default=None))
    return dict(niches=len(rows), with_feasible_front=len(with_front), front_size=_st([r["front_size"] for r in rows]), hv=_st([r["hv"] for r in rows]),
                hv_gain_second_half=_st(conv), hv_gain_gt_001=sum(1 for x in conv if x > 0.01),
                checked=sum(len(r["checked"]) for r in rows), checked_feasible=sum(r["checked_feasible"] for r in rows),
                evals=sum(r["evals"] for r in rows), robust_evals=sum(r["robust_evals"] for r in rows),
                l1_compare=dict(niches=len(cmp_), l1_feasible=sum(c["l1_feasible"] for c in cmp_),
                                dominated=sum(c["dominated_by_checked_front"] for c in cmp_),
                                gain=_st([c["gain"] for c in cmp_]), gain_positive=sum(1 for c in cmp_ if (c.get("gain") or 0) > 0)) if cmp_ else None,
                by_arch=by,
                note="гиперобъём - Монте-Карло 20 000 точек в нормированных целях (0 - лучшее, 1 - точка отсчёта); рост со второй половины "
                     "поколений - мера сходимости")


def global_front(rows):
    mem = [(r, ch) for r in rows for ch in r["checked"] if ch.get("feasible") and ch.get("G")]
    if not mem:
        return dict(n=0)
    G = np.array([ch["G"] for _, ch in mem])
    m = NS.front_mask(G, np.zeros(len(G)))
    sel = [mem[i] for i in np.where(m)[0]]
    return dict(n=len(sel), of=len(mem), by_arch=dict(Counter(r["arch"] for r, _ in sel).most_common()),
                members=[dict(niche=r["niche"], arch=r["arch"], key=f"{r['k']}:{ch['front_i']}", G=ch["G"], score_w=ch.get("score_w"),
                              score_f=ch.get("score_f"), eta_el=(ch.get("rec") or {}).get("eta_el"), V_env=(ch.get("rec") or {}).get("V_env"))
                         for r, ch in sel][:400])


def candidates(rows, ref):
    out = []
    for r in rows:
        for ch in r["checked"]:
            if not ch.get("feasible"):
                continue
            rec = dict(ch.get("rec") or {})
            out.append(dict(src="pareto" + (ref.get("tag") or ""), key=f"{r['k']}:{ch['front_i']}", idx=r["ref_idx"], kv=r["kv"], niche=r["niche"],
                            arch=r["arch"], g=ch["g"], ref_idx=r["ref_idx"], n0=ref["input"]["n"], seed=ref["input"]["seed"],
                            vigv=int(rec.get("vigv") or 0), heat=int(rec.get("inlet_heat") or 0), rec=rec))
    return out


CSV_COLS = ["k", "niche", "arch", "front_i", "feasible", "score_w", "score_f", "eta_el", "V_env", "n_parts", "steps", "trl_min", "SM_min", "reject",
            "G1_parts", "G2_plugplay", "G3_service", "G4_eta", "G5_size"]


def run_merge(found, kv, NN, shards, out_dir, sfx=""):
    t0 = time.time()
    rows, ref, checks, calls = merge_parts(found, kv, NN, shards)
    os.makedirs(out_dir, exist_ok=True)
    cands = candidates(rows, ref)
    p_c = os.path.join(out_dir, f"pareto_cands{sfx}.json.gz")
    with gzip.open(p_c, "wt", encoding="utf-8") as fh:
        json.dump(dict(kind="pareto", kv=kv, source=dict(code_sha=ref["code_sha"], l1_code_sha=ref.get("l1_code_sha"), input=ref["input"],
                                                         pop=ref["pop"], gens=ref["gens"], robust=ref["robust"], tag=ref.get("tag")),
                       genes=PN.GENES, candidates=cands), fh, ensure_ascii=False, separators=(",", ":"))
    p_csv = os.path.join(out_dir, f"pareto{sfx}.csv")
    with open(p_csv, "w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(CSV_COLS)
        for r in rows:
            for ch in r["checked"]:
                rec = ch.get("rec") or {}
                G = ch.get("G") or [None] * 5
                w.writerow([r["k"], r["niche"], r["arch"], ch["front_i"], ch.get("feasible"), ch.get("score_w"), ch.get("score_f"), rec.get("eta_el"),
                            rec.get("V_env"), rec.get("n_parts"), rec.get("steps"), rec.get("trl_min"), rec.get("SM_min"),
                            ch.get("reject") or "|".join(ch.get("od_bad") or []), *G])
    res = dict(title=f"Фронты Парето ниш (NSGA-II) - {kv}, {NN} ниш" + (f", устойчивый режим K = {ref['robust']}" if ref["robust"] else ""),
               software=ref.get("software"), date=PN.DATE, method=PN.__doc__.split("\n\n")[0], kv=kv, NN=NN, input=ref["input"], l1=ref.get("l1"),
               code_sha=ref["code_sha"], l1_code_sha=ref.get("l1_code_sha"), params=dict(pop=ref["pop"], gens=ref["gens"], check=ref["check"],
                                                                                       robust=ref["robust"], p_target=ref.get("p_target")),
               genes=ref["genes"], objectives=ref["objectives"], merge_checks=checks,
               timing=dict(cpu_s=round(sum(c.get("cpu_s") or 0 for c in calls), 1), evals=sum(c.get("evals") or 0 for c in calls),
                           robust_evals=sum(c.get("robust_evals") or 0 for c in calls), calls=calls),
               summary=summarize(rows), global_front=global_front(rows),
               niches=[{k: r.get(k) for k in ("k", "niche", "arch", "ref_idx", "n_l1_seeds", "n_l0_seeds", "evals", "feasible_evals", "robust_evals",
                                               "hv", "half", "front_size", "checked_feasible", "best_score_w", "best_score_f", "l1_compare", "t_cpu_s", "t_cpu")}
                       for r in rows],
               files=dict(candidates=os.path.basename(p_c), csv=os.path.basename(p_csv), n_candidates=len(cands)))
    res["timing"]["merge_s"] = round(time.time() - t0, 1)
    p_res = os.path.join(out_dir, f"pareto_results{sfx}.json")
    with open(p_res, "w", encoding="utf-8") as fh:
        json.dump(res, fh, ensure_ascii=False, indent=1)
    return dict(results=p_res, candidates=p_c, csv=p_csv, res=res, rows=rows)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--kv", required=True)
    ap.add_argument("--nn", type=int, required=True, help="число ниш NN (имя частей _n<NN>_)")
    ap.add_argument("--shards", type=int, required=True)
    ap.add_argument("--tag", default="")
    ap.add_argument("--in", dest="dirs", nargs="+", required=True)
    ap.add_argument("--out-dir", default=HERE)
    ap.add_argument("--suffix", default="")
    a = ap.parse_args(argv)
    try:
        out = run_merge(find_parts(a.dirs, a.kv, a.nn, a.shards, a.tag), a.kv, a.nn, a.shards, a.out_dir, a.suffix)
    except MergeError as e:
        print("слияние отклонено:", e)
        return 2
    s = out["res"]["summary"]
    print(f"{a.kv}: {s['niches']} ниш, с допустимым фронтом {s['with_feasible_front']}, проверено {s['checked']} (допустимы {s['checked_feasible']}); "
          f"кандидатов {out['res']['files']['n_candidates']} -> {out['results']}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
