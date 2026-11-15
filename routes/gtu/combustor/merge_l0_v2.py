"""Камера сгорания KV, уровень 0 v2: слияние частей (потоковое), полнота, цель «низкие выбросы 30-100 % на одном
клапане», рейтинг с множителем готовности, фронт Парето, сходимость доли по объёму выборки.

Проверки - как merge_l0.py (v1): один объём, зерно, отпечатки выборки, кода v2 и границ v2; части 0..S-1 по одному
разу (повтор части после прерывания спота - одна копия); индексы = range(n) без повторов. Неполное слияние - complete =
false и код выхода 2.

Слияние потоковое: часть читается, строки сжимаются до нужных полей, часть отпускается. Замер (строки пробы,
размноженные по индексам, 8 частей): пик 0,47 ГБ на 40 000 строк и 0,89 ГБ на 88 000 (части по 11 000 строк, как у
2^19 / 48) -> на 524 288 строк ~4-5 ГБ (линейная экстраполяция, оценка); v1 держал все строки - ~4,3 ГБ на 262 144.

Рейтинг - веса и нормировки v1 (merge_l0.WEIGHTS, objectives) × штраф жёстких целей v2 exp(-2·Σ min(нарушение, 3)) ×
множитель готовности (0,5 + 0,5·TRL/9) - решение владельца 15.11.2026, п. 1; TRL - наименьший из схемы, пассивного
элемента, подачи дизеля, охлаждения и изготовления (оценка).
Сходимость - доля вариантов, выполняющих новую цель, по вложенным префиксам выборки Соболя (первые 2^k индексов).

Запуск:
  python merge_l0_v2.py runs/<run_id> --out l0_aws_v2/comb_l0v2_KV-3_n<N>_results.json
  python merge_l0_v2.py /tmp/l0v2 --out comb_l0_v2_probe_results.json --top 30
"""
import argparse, glob, gzip, json, math, os, sys, time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in sys.path:
    sys.path.insert(0, HERE)
import merge_l0 as M1  # noqa: E402  (веса, нормировки, Парето - те же, что v1)

DATE = "2026-11-15"
KEYS = ("schema", "kv", "n", "seed", "shards", "fingerprint", "code_hash", "bc_key", "beta")


def part_files(paths):
    files = []
    for p in paths:
        if os.path.isdir(p):
            files += glob.glob(os.path.join(p, "**", "comb_l0v2_*.json"), recursive=True)
            files += glob.glob(os.path.join(p, "**", "comb_l0v2_*.json.gz"), recursive=True)
        else:
            files.append(p)
    return sorted(set(f for f in files if not f.endswith(".tmp")))


def _load(f):
    op = gzip.open if f.endswith(".gz") else open
    with op(f, "rt", encoding="utf-8") as fh:
        return json.load(fh)


def trl_mult(trl):
    return 0.5 + 0.5 * trl / 9.0


def compact(r, hard, new):
    """строка для рейтинга и отчёта (~1 КБ)"""
    c, g, m = r["concept"], r["gas"], r["module"]
    t = r["targets"]
    return dict(idx=r["idx"], scheme=c["scheme"], passive=c["passive"], split=c["split"], strat=c["strat"],
                recirc=c["recirc"], chamber=c["chamber"], injector=c["injector"],
                diesel=c["diesel"], cooling=c["cooling"], make=c["make"], n_el=c["n_el"], k_je=c["k_je"],
                a_pz=c["a_pz"], a_aux=c["a_aux"], k_pass=c["k_pass"], k_str=c["k_str"], x_lin=c["x_lin"],
                f_por=c["f_por"], off=c["off"], w_grp=c["w_grp"], l_ej=c["l_ej"], tau_pz_ms=c["tau_pz_ms"],
                tau_aux_ms=c["tau_aux_ms"], tau_bo_ms=c["tau_bo_ms"], dp_inj=c["dp_inj"], l_mix=c["l_mix"],
                d_jet_mm=r["derived"]["d_jet_mm"], d_air_mm=r["derived"]["d_air_mm"],
                p_sup_ratio=[x.get("p_sup_ratio") for x in g], fuel_share_anchor=[
                    next((x[k] for k in ("fuel_share_pilot", "fuel_share_cavity", "fuel_share_grpA") if k in x), None)
                    for x in g],
                a_eff=[x.get("a_eff_frac") for x in g], K_recirc=[x.get("K_recirc") for x in g],
                NOx=[x["NOx"] for x in g], CO=[x["CO"] for x in g], eta=[x["eta_c"] for x in g],
                margin_lbo=[x["margin_lbo"] for x in g], margin_all=[x.get("margin_all") for x in g],
                dpf=[x.get("dpf") for x in g], backflow_min=min((x.get("backflow_margin") or 9.0) for x in g),
                NOx_diesel=[x["NOx"] for x in r["diesel"]], dp_total=r["dp"]["dp_total"], PF=r["pf"]["PF"],
                T_wall_K=r["feed"]["T_metal_max_K"], H2_vol_max=r["flashback"]["H2_vol_max"], parts=m["parts"],
                swap_h=m["swap_h"], mass_kg=m["mass_kg"], plugplay=m["plugplay_score"], service=m["service_score"],
                make_score=m["make_score"], trl_min=m["trl_min_v2"], dyn_index=r["dyn_index"],
                f_stiff=r["dyn"]["f_stiff"], D_env=r["derived"]["D_env"], L_mod=r["derived"]["L_mod"],
                V_L=r["derived"]["V_L"], new_value=t[new]["value"], new_binding=t[new]["binding"],
                new_ok=t[new]["ok"], stable_30_100=t[new]["stable_30_100"], new_ok_strict=t[new]["ok_strict_all_zones"],
                load_range_ok=t["load_range"]["ok"], hard_ok=r["hard_ok"],
                hard_fail=[k for k in hard if not t[k]["ok"]], feed=r["feed"])


def merge(paths, out, top=50, keep_holders=3000):
    t0 = time.time()
    import comb_crn as CR   # noqa: E402
    import comb_v2 as V     # noqa: E402
    HARD = list(CR.TARGETS["funnel_interface"]["hard"]) + [V.NEW_TARGET]
    files = part_files(paths)
    rep = dict(files=len(files), problems=[])
    ref, seen_shards, used = None, {}, set()
    n = None
    have = None
    dup_idx = 0
    status, reasons = {}, {}
    by = {}
    feas = []                       # (idx, s-вектор, raw, penalty, trl_mult, compact)
    holders_idx = []
    stable_idx = []
    binding = {}
    ts, cpu, shard_rt = [], [], []
    inputs = []
    for f in files:
        p = _load(f)
        meta = {k: p.get(k) for k in KEYS}
        if ref is None:
            ref = meta
            n = ref["n"]
            have = bytearray(n)
            rep["ref"] = ref
        for k in KEYS:
            if meta[k] != ref[k]:
                rep["problems"].append(f"{os.path.basename(f)}: {k}={meta[k]!r} != {ref[k]!r}")
        seen_shards.setdefault(p["shard"], []).append(f)
        inputs.append(dict(file=os.path.basename(f), shard=p["shard"], rows=p["n_rows"], runtime_s=p["runtime_s"]))
        if p["shard"] in used or any(meta[k] != ref[k] for k in KEYS):
            continue
        used.add(p["shard"])
        shard_rt.append(p["runtime_s"])
        expect = set(range(p["shard"], n, ref["shards"]))
        got = [r["idx"] for r in p["rows"]]
        if set(got) != expect or len(got) != len(expect):
            rep["problems"].append(f"{os.path.basename(f)}: индексы части не совпадают с разбиением i % {ref['shards']} == "
                                   f"{p['shard']}")
        for r in p["rows"]:
            i = r["idx"]
            if have[i]:
                dup_idx += 1
                continue
            have[i] = 1
            status[r["status"]] = status.get(r["status"], 0) + 1
            for b in r.get("reject", []):
                reasons[b["code"]] = reasons.get(b["code"], 0) + 1
            c = r["concept"]
            key = f'{c["scheme"]}|{c["passive"]}'
            b = by.setdefault(key, dict(total=0, ok=0, stable_30_100=0, lowem=0, lowem_strict=0, load_range=0, hard_ok=0,
                                        reject={}, err=0))
            b["total"] += 1
            for x in r.get("reject", []):
                b["reject"][x["code"]] = b["reject"].get(x["code"], 0) + 1
            if "t_s" in r:
                ts.append(r["t_s"])
            if "cpu_s" in r:
                cpu.append(r["cpu_s"])
            if r["status"] == "ошибка":
                b["err"] += 1
            if r["status"] != "ok":
                continue
            cp = compact(r, HARD, V.NEW_TARGET)
            s, raw = M1.objectives(r)
            pen = math.exp(-M1.PENALTY_K * sum(min(r["targets"][t]["viol"], 3.0) for t in HARD))
            feas.append((i, s, raw, pen, trl_mult(cp["trl_min"]), cp))
            b["ok"] += 1
            b["stable_30_100"] += int(cp["stable_30_100"])
            b["lowem"] += int(cp["new_ok"])
            b["lowem_strict"] += int(cp["new_ok_strict"])
            b["load_range"] += int(cp["load_range_ok"])
            b["hard_ok"] += int(cp["hard_ok"])
            binding[cp["new_binding"].split("_")[0]] = binding.get(cp["new_binding"].split("_")[0], 0) + 1
            if cp["new_ok"]:
                holders_idx.append(i)
            if cp["stable_30_100"]:
                stable_idx.append(i)
        del p
    rep["duplicate_shards"] = {k: v for k, v in seen_shards.items() if len(v) > 1}
    rep["missing_shards"] = sorted(set(range(ref["shards"])) - set(seen_shards)) if ref else []
    missing = [i for i in range(n) if not have[i]] if ref else []
    rep["duplicate_indices"] = dup_idx
    rep["missing_indices"] = len(missing)
    rep["missing_indices_head"] = missing[:20]
    complete = bool(ref) and not rep["problems"] and not rep["missing_shards"] and not missing and dup_idx == 0
    # рейтинг
    W = M1.WEIGHTS
    wk = list(W)
    Smat = np.array([[x[1][k] for k in wk] for x in feas]) if feas else np.zeros((0, len(wk)))
    PEN = np.array([x[3] for x in feas])
    TM = np.array([x[4] for x in feas])
    vscore = lambda w: PEN * TM * (Smat @ np.array([w[k] for k in wk])) / sum(w.values()) if len(PEN) else np.array([])
    base = vscore(W)
    order = np.argsort(-base)
    top_rows = [dict(feas[j][5], score=float(base[j]), penalty=float(PEN[j]), trl_mult=float(TM[j]))
                for j in order[:top]]
    # фронт Парето: среди выполнивших все жёсткие, иначе среди выполнивших новую цель, иначе среди допустимых
    pf_rows, pool_name = [], None
    if feas:
        F = np.array([[(x[2][k] if sense == "min" else -x[2][k]) for k, sense in M1.PARETO_OBJ] for x in feas])
        hard_ok = np.array([x[5]["hard_ok"] for x in feas])
        new_ok = np.array([x[5]["new_ok"] for x in feas])
        if hard_ok.any():
            pool, pool_name = np.where(hard_ok)[0], "все жёсткие цели v2"
        elif new_ok.any():
            pool, pool_name = np.where(new_ok)[0], "новая цель (низкие выбросы 30-100 % на одном клапане)"
        else:
            pool, pool_name = np.arange(len(feas)), "все допустимые (ни один не выполнил новую цель)"
        pf_idx = pool[M1.pareto(F[pool])] if len(pool) else []
        pf_rows = sorted([dict(feas[j][5], score=float(base[j])) for j in pf_idx], key=lambda x: -x["score"])[:500]
    # чувствительность к весам
    sens = []
    if feas:
        top10 = set(int(j) for j in order[:10])
        best = int(order[0])
        for k in W:
            for fct in (0.5, 1.5):
                w = dict(W); w[k] *= fct
                o = np.argsort(-vscore(w))
                sens.append(dict(weight=k, factor=fct, top10_overlap=len(top10 & set(int(j) for j in o[:10])) / 10,
                                 best_idx=feas[int(o[0])][0], base_best_rank=int(np.where(o == best)[0][0]) + 1))
    rng = np.random.default_rng(20261115)
    freq = {}
    for _ in range(300 if feas else 0):
        ww = rng.dirichlet([20 * W[k] for k in wk])
        for j in np.argsort(-vscore(dict(zip(wk, ww))))[:10]:
            freq[feas[int(j)][0]] = freq.get(feas[int(j)][0], 0) + 1
    # сходимость по вложенным префиксам Соболя (первые 2^k индексов выборки)
    conv = []
    if ref:
        hset, sset = np.zeros(n, bool), np.zeros(n, bool)
        okset = np.zeros(n, bool)
        hset[holders_idx] = True
        sset[stable_idx] = True
        okset[[x[0] for x in feas]] = True
        ks = [2 ** k for k in range(6, 40) if 2 ** k < n] + [n]
        for m in ks:
            h, s_, o_ = int(hset[:m].sum()), int(sset[:m].sum()), int(okset[:m].sum())
            p_ = h / m
            se = math.sqrt(p_ * (1 - p_) / m)
            conv.append(dict(n=m, feasible=o_, stable_30_100=s_, lowem=h, share_lowem=p_, se_share=se,
                             rel_se=(se / p_) if h else None, upper95_if_zero=(3.0 / m) if h == 0 else None))
    for key, b in by.items():
        b["share_lowem_of_total"] = b["lowem"] / b["total"] if b["total"] else 0.0
    # по схеме × пассивный элемент: медианы среди допустимых
    med = {}
    for x in feas:
        cp = x[5]
        med.setdefault(f'{cp["scheme"]}|{cp["passive"]}', []).append(cp)
    for key, lst in med.items():
        b = by[key]
        for li, pct in enumerate((100, 75, 50, 30)):
            b[f"NOx_{pct}_median"] = float(np.median([r["NOx"][li] for r in lst]))
            b[f"eta_{pct}_median"] = float(np.median([r["eta"][li] for r in lst]))
        ratio = [r["fuel_share_anchor"][3] / r["fuel_share_anchor"][0] for r in lst if r["fuel_share_anchor"][0]]
        b["share_anchor_30_over_100_median"] = float(np.median(ratio)) if ratio else None
        b["share_anchor_30_over_100_max"] = float(max(ratio)) if ratio else None
        ae = [r["a_eff"][3] for r in lst if r["a_eff"][3] is not None]
        b["a_eff_30_min"] = float(min(ae)) if ae else None
        b["new_value_min"] = float(min(r["new_value"] for r in lst))
        b["best_new_idx"] = min(lst, key=lambda r: r["new_value"])["idx"]
    holders = sorted([x[5] for x in feas if x[5]["new_ok"]], key=lambda r: r["new_value"])[:keep_holders]
    nearest = sorted([x[5] for x in feas], key=lambda r: r["new_value"])[:30]
    feed = {}
    for j in order:
        cp = feas[int(j)][5]
        key = f'{cp["scheme"]}|{cp["passive"]}|{cp["chamber"]}'
        if key not in feed:
            feed[key] = dict(idx=cp["idx"], score=float(base[int(j)]), hard_ok=cp["hard_ok"], new_ok=cp["new_ok"],
                             **cp["feed"])
    res = dict(
        title="Камера сгорания KV-3, уровень 0 v2: пассивные схемы частичной нагрузки - слияние, новая цель, рейтинг",
        method=("потоковое слияние частей comb_l0_v2.py с контролем полноты; цель lowem_30_100_1valve (targets_v2.json); "
                "рейтинг - веса и нормировки v1 × штраф жёстких целей v2 exp(-2·Σ нарушений) × (0,5 + 0,5·TRL/9); фронт "
                "Парето по шести целям v1; сходимость доли - вложенные префиксы Соболя"),
        software=f"python {sys.version.split()[0]}, numpy {np.__version__}", date=DATE,
        status="уровень облика (0D + корреляции): числа - оценка для отбора на L1 (RANS), не прогноз конструкции; закон "
               "затухания расслоения и неравномерности при малом J - оценка (BETA_LOWJ), проверка - уровень 1",
        complete=complete, completeness=rep, inputs=inputs,
        weights=W, penalty_k=M1.PENALTY_K, hard_targets=HARD, trl_multiplier="0,5 + 0,5·TRL/9",
        constants=dict(BETA_LOWJ=(ref or {}).get("beta"), DP_BACK_MIN=V.DP_BACK_MIN, DPF_REF=V.DPF_REF, LOADS=V.LOADS),
        counts=dict(n=n, rows=int(sum(status.values())), status=status, reject_codes=reasons, feasible=len(feas),
                    stable_30_100=len(stable_idx), lowem_30_100_1valve=len(holders_idx),
                    lowem_strict_all_zones=int(sum(x[5]["new_ok_strict"] for x in feas)),
                    load_range_base=int(sum(x[5]["load_range_ok"] for x in feas)),
                    hard_ok=int(sum(x[5]["hard_ok"] for x in feas)), binding_term_among_feasible=binding),
        timing=dict(t_variant_mean_s=float(np.mean(ts)) if ts else None, t_variant_median_s=float(np.median(ts)) if ts else None,
                    cpu_variant_mean_s=float(np.mean(cpu)) if cpu else None,
                    cpu_variant_median_s=float(np.median(cpu)) if cpu else None,
                    cpu_variant_max_s=float(np.max(cpu)) if cpu else None, shard_runtime_s=shard_rt),
        by_scheme_passive=by, convergence=conv, holders=holders, nearest_to_new_target=nearest,
        top=top_rows, pareto=pf_rows, pareto_pool=pool_name, pareto_objectives=M1.PARETO_OBJ,
        sensitivity=dict(one_at_a_time=sens, dirichlet_top10_freq=[dict(idx=i, freq=f / 300) for i, f in
                                                                   sorted(freq.items(), key=lambda x: -x[1])[:15]]),
        feed_to_engine_funnel=dict(note="лучший по рейтингу v2 в каждой тройке (схема, пассивный элемент, тип камеры); "
                                        "поля - targets.json funnel_interface + passive, trl_min, NOx/CO на 100/75/50/30 %",
                                   best=feed),
        runtime_s=None)
    res["runtime_s"] = round(time.time() - t0, 1)
    json.dump(CR._rnd(res), open(out, "w", encoding="utf-8"), ensure_ascii=False, indent=1)
    return complete, res


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("paths", nargs="+")
    ap.add_argument("--out", required=True)
    ap.add_argument("--top", type=int, default=50)
    a = ap.parse_args(argv)
    complete, res = merge(a.paths, a.out, a.top)
    c = res["counts"]
    print(f"полнота: {complete}; строк {c['rows']} из {c['n']}; допустимых {c['feasible']}; держат 30-100 % "
          f"{c['stable_30_100']}; новая цель {c['lowem_30_100_1valve']}; все жёсткие v2 {c['hard_ok']}; запись {a.out}")
    return 0 if complete else 2


if __name__ == "__main__":
    sys.exit(main())
