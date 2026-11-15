"""Воронка двигателя-заменителя, уровень 0, часть 2: слияние частей, рейтинг, ниши и отбор на уровень 1.

  python merge_l0.py --kv KV-3 --n 200000 --shards 30 --in runs/<run_id> --out-dir l0_KV-3 --top 1000 --target 4000
  python merge_l0.py --probe          # короткая локальная проба: части funnel_l0 (KV-3, N=256, 10 частей) -> слияние

Слияние: все файлы частей funnel_l0_<KV>_n<N>_s<k>of<S>.json[.gz] из --in (рекурсивно); у всех частей одинаковые kv, n,
seed, shards, отпечаток выборки и отпечаток кода; у части k ровно индексы i % S == k; ни один индекс 0..n-1 не потерян и
не задвоен (копия части после повтора задачи допускается, только если она совпадает побайтно по записям). Любое
нарушение - отказ с перечнем (код выхода 2).

Рейтинг: критерии - баллы 0-10 на постоянных шкалах (нормировка - деление на 10; от выборки не зависит, поэтому оценка
концепции одна и та же при любом N и любом делении на части); взвешенная оценка по весам владельца; место среди
допустимых (равенство - по меньшему индексу); фронт Парето по четырём критериям; чувствительность к весам (±50 % на
каждый вес по очереди и 200 случайных наборов в тех же пределах); архив ниш (MAP-Elites); отбор на уровень 1: лучшие K
каждой допустимой ниши + общий топ + фронт Парето [+ почти допустимые из ниш без допустимых - отдельной причиной];
сходимость: топ-50 и ниши на первых N/4, N/2 и N концепциях (вариант i не зависит от N - свойство Соболя, тест).
Пишет funnel_l0_results<sfx>.json, funnel_l0<sfx>.csv (все варианты), funnel_l0_select<sfx>.json.
"""
import argparse, csv, glob, gzip, hashlib, json, math, os, re, shutil, sys, tempfile, time
from collections import Counter, defaultdict

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
GTU = os.path.dirname(HERE)
if HERE not in sys.path:
    sys.path.insert(0, HERE)
import funnel_l0 as L0           # noqa: E402
import funnel_space as FS        # noqa: E402

DATE = L0.DATE
TOP_LIST = 40          # топ в результатах (со стабильностью при разных весах)
TOP_CONV = 50          # топ для сходимости при удвоении N
N_RANDOM_W = 200       # случайных наборов весов (каждый вес × U(0,5; 1,5), затем нормировка)
EXTRA_COLS = ["rank", "niche_rank", "niche_size", "pareto", "selected", "select_reason"]


class MergeError(Exception):
    pass


# =========================================================================== чтение и проверка частей
def find_parts(dirs, kv, n, shards):
    pat = re.compile(rf"^{L0.PART_PREFIX}_{re.escape(kv)}_n{n}_s(\d{{3}})of{shards:03d}\.json(\.gz)?$")
    out = []
    for d in dirs:
        for p in glob.glob(os.path.join(d, "**", f"{L0.PART_PREFIX}_{kv}_n{n}_s*of{shards:03d}.json*"), recursive=True):
            m = pat.match(os.path.basename(p))
            if m:
                out.append((int(m.group(1)), p))
    return sorted(out, key=lambda x: (x[0], x[1]))


def load_part(path):
    op = gzip.open if path.endswith(".gz") else open
    with op(path, "rt", encoding="utf-8") as fh:
        return json.load(fh)


def merge_parts(found, kv, n, shards, seed=L0.SEED):
    """found: [(k, путь)]. Возврат: columns, rows (по индексу), meta, checks. Нарушения -> MergeError"""
    checks = dict(files=len(found), shards_expected=shards, missing=[], duplicate_same=[], duplicate_conflict=[],
                  meta_mismatch=[], index_errors=[], unreadable=[])
    ref, by = None, {}
    for k, p in found:
        try:
            d = load_part(p)
        except (OSError, EOFError, ValueError) as e:      # обрезанный файл прерванной попытки: часть возьмётся из повтора
            checks["unreadable"].append(dict(file=os.path.relpath(p), error=type(e).__name__))
            continue
        meta = {x: d.get(x) for x in ("kv", "n", "seed", "shards", "shard", "fingerprint", "code_sha", "level")}
        want = dict(kv=kv, n=n, seed=seed, shards=shards, shard=k, level=L0.LEVEL)
        bad = {x: (meta[x], v) for x, v in want.items() if meta[x] != v}
        if ref is None:
            ref = dict(meta, columns=d["columns"], target=d.get("target"), calib=d.get("calib"), weights=d.get("weights"),
                       software=d.get("software"), method=d.get("method"))
        else:
            for x in ("fingerprint", "code_sha"):
                if meta[x] != ref[x]:
                    bad[x] = (meta[x], ref[x])
            if d["columns"] != ref["columns"]:
                bad["columns"] = ("другой состав колонок", "")
        if bad:
            checks["meta_mismatch"].append(dict(file=os.path.basename(p), diff={x: list(v) for x, v in bad.items()}))
            continue
        rows = d["rows"]
        for r in rows:
            for j, v in enumerate(r):
                if isinstance(v, str):
                    r[j] = sys.intern(v)
        exp = FS.shard_indices(n, k, shards)
        got = [r[0] for r in rows]
        if got != exp:
            lost = sorted(set(exp) - set(got)); extra = sorted(set(got) - set(exp))
            checks["index_errors"].append(dict(file=os.path.basename(p), shard=k, lost=lost[:20], n_lost=len(lost),
                                               foreign=extra[:20], n_foreign=len(extra), dup=len(got) - len(set(got))))
            continue
        h = hashlib.sha256(json.dumps(rows, ensure_ascii=False, separators=(",", ":")).encode()).hexdigest()
        if k in by:
            (checks["duplicate_same"] if by[k]["hash"] == h else checks["duplicate_conflict"]).append(
                dict(shard=k, file=os.path.basename(p), kept=os.path.basename(by[k]["path"])))
            continue
        by[k] = dict(rows=rows, hash=h, path=p, call=d.get("call"), workers=d.get("workers"))
    checks["missing"] = [k for k in range(shards) if k not in by]
    errs = {x: checks[x] for x in ("missing", "duplicate_conflict", "meta_mismatch", "index_errors") if checks[x]}
    if errs or ref is None:
        raise MergeError(json.dumps(errs or {"files": "частей не найдено"}, ensure_ascii=False))
    rows = sorted((r for v in by.values() for r in v["rows"]), key=lambda r: r[0])
    idx = [r[0] for r in rows]
    if idx != list(range(n)):
        raise MergeError(f"индексы после слияния не равны 0..{n - 1}: {len(idx)} записей, уникальных {len(set(idx))}")
    checks.update(indices_ok=True, n_rows=len(rows), duplicates=0, fingerprint=ref["fingerprint"], code_sha=ref["code_sha"])
    fp_local = FS.fingerprint(n, seed)
    checks["fingerprint_local"] = dict(value=fp_local, match=fp_local == ref["fingerprint"],
                                       note="отпечаток матрицы выборки в окружении слияния; несовпадение - другие версии scipy/numpy: "
                                            "уровень 1 обязан считать с закреплёнными версиями задач (подробности - вне слияния)")
    calls = {}
    for k, v in sorted(by.items()):
        c = v["call"] or {}
        key = (c.get("started"), tuple(c.get("shards") or []))
        calls[key] = dict(shards=c.get("shards"), concepts=c.get("concepts"), elapsed_s=c.get("elapsed_s"),
                          s_per_concept_core=c.get("s_per_concept_core"), workers=v["workers"], host_cpus=c.get("host_cpus"))
    return ref["columns"], rows, ref, checks, list(calls.values())


# =========================================================================== рейтинг
def parse_weights(s):
    if not s:
        return dict(L0.WEIGHTS)
    w = {}
    for tok in s.split(","):
        k, v = tok.split("=")
        w[k.strip()] = float(v)
    return L0.check_weights(w)


class Table:
    """объединённые записи: строки-списки + числовые столбцы numpy"""

    def __init__(self, columns, rows):
        self.cols, self.rows = columns, rows
        self.C = {c: j for j, c in enumerate(columns)}
        self.n = len(rows)
        self.idx = np.array([r[0] for r in rows], dtype=np.int64)
        self.feas = np.array([bool(r[self.C["feasible"]]) for r in rows])
        self.viol = self.num("violation")
        self.has = ~np.isnan(self.num("s_parts"))
        self.S = np.column_stack([self.num(L0.CRIT_COL[k]) for k in L0.CRIT])
        self.niche = [r[self.C["niche"]] for r in rows]
        self.arch = [r[self.C["arch"]] for r in rows]

    def num(self, col):
        j = self.C[col]
        return np.array([np.nan if r[j] is None else float(r[j]) for r in self.rows])

    def get(self, pos, col):
        return self.rows[pos][self.C[col]]

    def s_effmake(self, split_em):
        """четвёртый критерий при другом внутреннем делении (для чувствительности)"""
        es, mk = self.num("s_etasize"), self.num("s_make")
        return split_em["etasize"] * es + split_em["make"] * mk


def scores(T, w, S=None):
    S = T.S if S is None else S
    return S @ np.array([w[k] for k in L0.CRIT])


def order_feasible(T, sc, mask=None):
    """позиции допустимых (и в маске) по убыванию оценки; равенство - по меньшему индексу"""
    m = T.feas.copy() if mask is None else (T.feas & mask)
    pos = np.nonzero(m)[0]
    return pos[np.lexsort((T.idx[pos], -np.round(sc[pos], 9)))]


def pareto_layers(F, max_layers=3):
    """слои недоминируемых (максимум по всем столбцам F); 0 - глубже max_layers"""
    layer = np.zeros(len(F), dtype=int)
    rem = np.arange(len(F))
    for L in range(1, max_layers + 1):
        if not len(rem):
            break
        sub = F[rem]
        order = np.argsort(-sub.sum(axis=1), kind="stable")
        front = []
        for j in order:
            if front:
                fr = sub[front]
                dom = np.all(fr >= sub[j], axis=1) & np.any(fr > sub[j], axis=1)
                if dom.any():
                    continue
            front.append(j)
        layer[rem[front]] = L
        rem = np.delete(rem, front)
    return layer


# =========================================================================== ниши и отбор
def niche_lists(T, sc):
    """ниша -> позиции допустимых по убыванию оценки; все ниши -> позиции"""
    feas_by, all_by = defaultdict(list), defaultdict(list)
    for p in order_feasible(T, sc):
        feas_by[T.niche[p]].append(int(p))
    for p in range(T.n):
        all_by[T.niche[p]].append(p)
    return feas_by, all_by


def select(T, sc, k=3, top=40, pareto_pos=(), k_near=0, near_max=0.25, target=None):
    """отбор на уровень 1. Возврат: dict(pos -> причины), K, флаги. Каждая допустимая ниша - не меньше одного варианта."""
    feas_by, all_by = niche_lists(T, sc)
    order = order_feasible(T, sc)

    def base(K):
        sel = defaultdict(list)
        for nk, lst in feas_by.items():
            for r, p in enumerate(lst[:K]):
                sel[p].append(f"ниша: место {r + 1} из {len(lst)}")
        for r, p in enumerate(order[:top]):
            sel[int(p)].append(f"общий топ: место {r + 1}")
        for p in pareto_pos:
            sel[int(p)].append("фронт Парето")
        return sel

    kmax = max([len(v) for v in feas_by.values()] or [1])
    K, exceeded = k, False
    if target is not None:
        K = 1
        if len(base(1)) > target:
            exceeded = True
        else:
            lo, hi = 1, kmax
            while lo < hi:
                mid = (lo + hi + 1) // 2
                if len(base(mid)) <= target:
                    lo = mid
                else:
                    hi = mid - 1
            K = lo
    sel = base(K)
    near = []
    if k_near > 0:
        for nk, lst in all_by.items():
            if nk in feas_by:
                continue
            cand = [p for p in lst if T.has[p] and T.viol[p] <= near_max]
            cand.sort(key=lambda p: (T.viol[p], -sc[p], T.idx[p]))
            near += cand[:k_near]
        near.sort(key=lambda p: (T.viol[p], -sc[p], T.idx[p]))
        room = (target - len(sel)) if target is not None else len(near)
        for p in near[:max(0, room)]:
            sel[p].append(f"ниша без допустимых: превышение {T.viol[p]:.3f} (уровень 1 оптимизирует непрерывные оси)")
    return dict(sel=sel, K=K, kmax=kmax, target_exceeded=exceeded, near_candidates=len(near),
                feasible_niches=sorted(feas_by), niche_sizes={nk: len(v) for nk, v in feas_by.items()})


# =========================================================================== чувствительность к весам
def scenarios(base_w, split_em):
    out = [dict(name="базовые веса", w=dict(base_w), split=dict(split_em))]
    for k in L0.CRIT:
        for f in (0.5, 1.5):
            w = dict(base_w); w[k] *= f
            s = sum(w.values())
            out.append(dict(name=f"{k} ×{f:g}", w={x: v / s for x, v in w.items()}, split=dict(split_em)))
    for f in (0.5, 1.5):
        e = split_em["etasize"] * f
        out.append(dict(name=f"внутри 4-го: КПД и размер ×{f:g}", w=dict(base_w),
                        split={"etasize": e / (e + split_em["make"]), "make": split_em["make"] / (e + split_em["make"])}))
    return out


def sensitivity(T, base_w, sel_args, top_n=TOP_LIST, n_random=N_RANDOM_W, seed=L0.SEED):
    sc0 = scores(T, base_w)
    base_order = order_feasible(T, sc0)
    base_top = [int(p) for p in base_order[:top_n]]
    base_sel = set(select(T, sc0, **sel_args)["sel"])
    oat, ranks = [], defaultdict(list)
    for s in scenarios(base_w, L0.SPLIT_EFFMAKE):
        S = T.S.copy()
        if s["split"] != L0.SPLIT_EFFMAKE:
            S[:, L0.CRIT.index("effmake")] = T.s_effmake(s["split"])
        sc = scores(T, s["w"], S)
        order = order_feasible(T, sc)
        pos_rank = {int(p): r + 1 for r, p in enumerate(order)}
        top = [int(p) for p in order[:top_n]]
        sel = set(select(T, sc, **sel_args)["sel"])
        for p in base_top:
            ranks[p].append(pos_rank.get(p))
        oat.append(dict(scenario=s["name"], weights={k: round(v, 4) for k, v in s["w"].items()},
                        split_effmake={k: round(v, 4) for k, v in s["split"].items()},
                        top1_idx=int(T.idx[order[0]]) if len(order) else None,
                        top1_same=bool(len(order) and base_top and int(order[0]) == base_top[0]),
                        top_overlap=len(set(top) & set(base_top)), top_n=min(top_n, len(order)),
                        select_jaccard=round(len(sel & base_sel) / max(1, len(sel | base_sel)), 4)))
    rng = np.random.default_rng(seed)
    bw = np.array([base_w[k] for k in L0.CRIT])
    M = rng.uniform(0.5, 1.5, size=(n_random, len(bw))) * bw
    M /= M.sum(axis=1, keepdims=True)
    fpos = np.nonzero(T.feas)[0]
    cnt = Counter()
    top1 = Counter()
    if len(fpos):
        SC = T.S[fpos] @ M.T                                     # допустимые × наборы весов
        for j in range(n_random):
            o = np.lexsort((T.idx[fpos], -np.round(SC[:, j], 9)))
            for p in fpos[o[:top_n]]:
                cnt[int(p)] += 1
            top1[int(T.idx[fpos[o[0]]])] += 1
    stab = {p: dict(random_share=round(cnt[p] / n_random, 3),
                    oat_in_top=sum(1 for r in ranks[p][1:] if r is not None and r <= top_n), oat_n=len(ranks[p]) - 1,
                    rank_min=min(r for r in ranks[p] if r is not None), rank_max=max(r for r in ranks[p] if r is not None))
            for p in base_top}
    return dict(oat=oat, random=dict(n=n_random, seed=seed, rule="каждый вес × U(0,5; 1,5), затем сумма = 1",
                                     top1_counts={str(k): v for k, v in top1.most_common(10)},
                                     base_top_always=sum(1 for p in base_top if cnt[p] == n_random),
                                     base_top_share_ge_0_8=sum(1 for p in base_top if cnt[p] >= 0.8 * n_random),
                                     outsiders_entering=sum(1 for p in cnt if p not in stab)),
                stability=stab)


# =========================================================================== сходимость при удвоении N
def convergence(T, w, sel_args, top_n=TOP_CONV):
    sc = scores(T, w)
    sizes = sorted({max(1, T.n // 4), max(1, T.n // 2), T.n})
    out, prev_top = [], None
    full_top = set(int(p) for p in order_feasible(T, sc)[:top_n])
    for m in sizes:
        mask = T.idx < m
        order = order_feasible(T, sc, mask)
        top = [int(p) for p in order[:top_n]]
        feas_n = sorted({T.niche[p] for p in np.nonzero(T.feas & mask)[0]})
        row = dict(N=int(m), feasible=int((T.feas & mask).sum()), feasible_niches=len(feas_n),
                   archetypes_feasible=len({T.arch[p] for p in np.nonzero(T.feas & mask)[0]}),
                   top_n=len(top), best_score=round(float(sc[order[0]]), 4) if len(order) else None,
                   mean_top=round(float(np.mean(sc[top])), 4) if top else None,
                   top_overlap_with_full=len(set(top) & full_top))
        if prev_top is not None:
            row["new_in_top_vs_half"] = len(set(top) - set(prev_top[0]))
            row["new_feasible_niches_vs_half"] = len(set(feas_n) - set(prev_top[1]))
        out.append(row)
        prev_top = (top, feas_n)
    return out


# =========================================================================== выборочная сверка цикла с gtu_cycle.cycle
def cycle_check(T, n, kv, seed, k_max=6):
    """двухвальные одноядерные без рекуператора, кольцевая или выносная камера: ядро уровня 0 против gtu_cycle.cycle
    (Cantera напрямую) при тех же КПД (компрессор - эквивалентный политропный), долях воздуха и потерях"""
    E = FS._E()
    cwd = os.getcwd()
    if GTU not in sys.path:
        sys.path.insert(0, GTU)
    from gtu_cycle import cycle
    os.chdir(cwd)
    C = T.C
    cand = [p for p in range(T.n) if T.has[p] and T.rows[p][C["rotors"]] == 2 and T.rows[p][C["cores"]] == 1
            and T.rows[p][C["recup"]] == "нет" and T.rows[p][C["comb"]] in ("кольцевая", "выносная")]
    sc = T.num("score")
    cand.sort(key=lambda p: (not T.feas[p], -sc[p], T.idx[p]))
    rows = []
    for p in cand[:k_max]:
        i = int(T.idx[p])
        c = FS.concept(i, n, kv, seed)
        sol = L0.solve(c, FS.target(kv))
        pf = sol["perf"]
        comp = pf["comp"]
        eta_pc = E.RA * math.log(c["pi"]) / (E.s0(comp[-1]["T2"]) - E.s0(comp[0]["T1"]))
        xn, xr = FS.cool_air(c)
        g = cycle(c["pi"], c["t41"], eta_pc, pf["eta_gg"], pf["eta_pt"], xn, xr, 0.015,
                  dp_comb=0.05 if c["comb"] == "выносная" else 0.04, dp_duct=0.01)
        rows.append(dict(idx=i, arch=c["arch"], feasible=bool(T.feas[p]), pi=round(c["pi"], 3), t41=round(c["t41"], 1),
                         eta_l0=round(pf["eta_th"] * 100, 3), eta_cycle=round(g["eff"] * 100, 3),
                         d_eta_pp=round((pf["eta_th"] - g["eff"]) * 100, 3), w_l0=round(pf["w"], 2), w_cycle=round(g["sp_work_kJkg"], 2),
                         EGT_l0=round(pf["EGT"], 1), EGT_cycle=round(g["EGT_C"], 1)))
    return dict(method="ядро уровня 0 (таблицы gtu_explore) против gtu_cycle.cycle (Cantera) при одинаковых КПД узлов",
                n=len(rows), max_abs_d_eta_pp=max([abs(r["d_eta_pp"]) for r in rows], default=None), rows=rows)


# =========================================================================== сводки
def _q(v):
    v = np.asarray([x for x in v if not np.isnan(x)])
    if not len(v):
        return None
    return dict(n=int(len(v)), min=round(float(v.min()), 3), p10=round(float(np.percentile(v, 10)), 3),
                p50=round(float(np.median(v)), 3), p90=round(float(np.percentile(v, 90)), 3), max=round(float(v.max()), 3),
                std=round(float(v.std()), 3))


def summary(T, sc):
    rej = Counter()
    for r in T.rows:
        for code in (r[T.C["reject"]] or "").split("|"):
            if code:
                rej[code.split(":")[0] if code.startswith("size:") else code.rstrip("0123456789")] += 1
    by = {}
    for a in FS.ARCH_KEYS:
        pos = [p for p in range(T.n) if T.arch[p] == a]
        fp = [p for p in pos if T.feas[p]]
        first = Counter((T.get(p, "reject") or "").split("|")[0] for p in pos if not T.feas[p])
        d = dict(n=len(pos), feasible=len(fp), first_reject=dict(first.most_common(4)))
        if "repair" in T.C:
            rp = [(T.get(p, "repair") or "") for p in pos]
            d["sampler"] = dict(
                repaired=sum(1 for x in rp if "->" in x and "исходная не замыкает" in x),
                no_closure=sum(1 for x in rp if x.startswith("нет замыкающей")),
                no_closure_reason=Counter(x.split(" - ", 1)[1].split(":")[0] for x in rp if x.startswith("нет замыкающей")).most_common(3),
                design_flags=Counter(f.split(":")[0] for p in pos for f in (T.get(p, "dflags") or "").split("|") if f).most_common(4))
        if fp:
            b = max(fp, key=lambda p: (sc[p], -T.idx[p]))
            d.update(best_idx=int(T.idx[b]), best_score=round(float(sc[b]), 4), eta_el=_q(T.num("eta_el")[fp]),
                     n_parts=_q(T.num("n_parts")[fp]))
        else:
            nz = [p for p in pos if T.has[p]]
            if nz:
                b = min(nz, key=lambda p: (T.viol[p], T.idx[p]))
                d.update(nearest_idx=int(T.idx[b]), nearest_violation=round(float(T.viol[b]), 4), nearest_reject=T.get(b, "reject"))
        by[a] = d
    out = dict(n=T.n, feasible=int(T.feas.sum()), feasible_share=round(float(T.feas.mean()), 5), reject_reasons=dict(rej.most_common()),
               archetypes_feasible=sum(1 for a in by.values() if a["feasible"]), by_arch=by)
    if "dflags" in T.C:
        out["design_flags"] = dict(Counter(f.split(":")[0] for p in range(T.n) for f in (T.get(p, "dflags") or "").split("|") if f).most_common())
        out["no_closure"] = sum(1 for p in range(T.n) if (T.get(p, "repair") or "").startswith("нет замыкающей"))
    return out


def normalization(T, w):
    stats = {k: _q(T.S[T.feas, j]) for j, k in enumerate(L0.CRIT)}
    std = {k: (stats[k]["std"] if stats[k] else 0.0) for k in L0.CRIT}
    tot = sum(w[k] * std[k] for k in L0.CRIT) or 1.0
    return dict(method="постоянные шкалы 0-10, нормировка - деление на 10 (в оценке - сумма w·s, w с суммой 1)",
                why="оценка концепции не зависит от выборки: слияние частей = расчёт целиком, место концепции сравнимо при любом N; "
                    "нормировка по размаху выборки сделала бы оценку зависимой от N и состава частей",
                scales={"parts": "parts_model: 10·N_HALF/(N_HALF + N_eff), опорная ведомость frame3 = 5",
                        "plugplay": "plugplay_model: 10·Σ w_i·s_i по 10 интерфейсам gtu-dropin",
                        "service": "service_model: модули, шаги, инструмент, защита от ошибок, самодиагностика",
                        "effmake": f"{L0.SPLIT_EFFMAKE['etasize']}·(0,5·s_eta + 0,5·s_size) + {L0.SPLIT_EFFMAKE['make']}·s_make; "
                                   f"s_eta = 5 + 5·(КПД - КПД исходной)/{L0.ETA_SPAN_PP:g} п. п.; s_size = 10·(1 - V/V_отсека), "
                                   "V_отсека = min(ш, в)²·L_пакета/3; s_make - make_model (оценка)"},
                feasible_stats=stats,
                effective_weights={k: round(w[k] * std[k] / tot, 4) for k in L0.CRIT},
                effective_note="доля разброса взвешенной оценки от каждого критерия: w·σ / Σ w·σ по допустимым - критерий с малым "
                               "разбросом влияет на выбор меньше своего веса")


def top_entries(T, sc, order, stab, layer, n, kv, seed, k=TOP_LIST, detail=True):
    out = []
    keys = ["arch", "niche", "comp", "turb", "nt", "cool", "recup", "comb", "bear", "cores", "rotor_join", "m_cw", "m_tw", "pi", "t41",
            "eta_th_shaft", "eta_el", "w_kJkg", "dG_pct", "dEGT", "n_rpm", "D2", "D_env", "L_env", "W_env", "H_env", "V_env", "mass_kg",
            "s_parts", "s_plugplay", "s_service", "s_make", "s_eta", "s_size", "s_effmake", "n_parts", "n_assemblies", "n_modules",
            "heaviest_kg", "steps", "n_tools", "trl_min", "techs"]
    for r, p in enumerate(order[:k]):
        p = int(p)
        e = dict(rank=r + 1, idx=int(T.idx[p]), score=round(float(sc[p]), 4), pareto_layer=int(layer[p]))
        e.update({x: T.get(p, x) for x in keys})
        e["stability"] = stab.get(p)
        if detail:
            _, det = L0.evaluate_l0(int(T.idx[p]), n, kv, seed, detail=True)
            e["detail"] = dict(geometry=det["geometry"], modules=det["service"]["modules"], tools=det["service"]["tools"],
                               plugplay=det["plugplay"], trl_by=det["make"]["trl_by"])
        out.append(e)
    return out


# =========================================================================== запись
def write_csv(path, T, extra):
    op = gzip.open if path.endswith(".gz") else open
    with op(path, "wt", encoding="utf-8", newline="") as fh:
        wr = csv.writer(fh)
        wr.writerow(T.cols + EXTRA_COLS)
        for p, r in enumerate(T.rows):
            wr.writerow(["" if v is None else v for v in r] + ["" if v is None else v for v in extra[p]])


def inputs(parts_meta):
    files = ["funnel/funnel_l0.py", "funnel/merge_l0.py", "funnel/funnel_space.py", "funnel/parts_model.py", "funnel/service_model.py",
             "funnel/plugplay_model.py", "funnel/make_model.py", "gtu_cycle.py", "gtu_explore.py", "gtu_family_results.json",
             "frame3/frame3_point_results.json", "fleet/gtu_fleet_taurus60.json"]
    fields = {"frame3/frame3_point_results.json": "stage1/stage2: m_per_side, T01, p01, w, eta_is, eta_basis.n_ref_rpm, k_speed, "
                                                  "d_eta_used -> ns_opt и K поправки по быстроходности",
              "gtu_family_results.json": "slots: расчётные точки KV (P, G0, EGT0, eta0, частота, пакет)",
              "gtu_cycle.py": "cycle() - выборочная сверка ядра цикла",
              "fleet/gtu_fleet_taurus60.json": "engine_kg Taurus 60 MD - сопоставление массы (только для справки)"}
    out = [dict(file="routes/gtu/" + f, field=fields.get(f, "код расчёта"),
                sha256=hashlib.sha256(open(os.path.join(GTU, f), "rb").read()).hexdigest()) for f in files]
    out.append(dict(file="части уровня 0", field=f"{parts_meta['files']} файлов, отпечаток выборки {parts_meta['fingerprint']}, "
                                                  f"код {parts_meta['code_sha']}"))
    return out


CSV_BRANCH_MAX_MB = 5.0      # CSV крупнее - в ветку не кладётся: --csv-dir вне ветки (S3, routes/gtu/platform/aws/runs/ - в .gitignore)


def run_merge(found, kv, n, shards, out_dir, sfx="", weights=None, seed=L0.SEED, k=3, top=40, use_pareto=True, k_near=1,
              near_max=0.25, target=None, csv_gz=False, detail_top=TOP_LIST, cycle_k=6, extra_results=None, csv_dir=None):
    t0 = time.time()
    cols, rows, ref, checks, calls = merge_parts(found, kv, n, shards, seed)
    w = L0.check_weights(dict(weights or L0.WEIGHTS))
    T = Table(cols, rows)
    sc = scores(T, w)
    jsc = T.C["score"]
    for p, r in enumerate(T.rows):                      # оценка при весах слияния
        r[jsc] = None if np.isnan(sc[p]) else round(float(sc[p]), 4)
    order = order_feasible(T, sc)
    rank = {int(p): r + 1 for r, p in enumerate(order)}
    fpos = np.nonzero(T.feas)[0]
    layer = np.zeros(T.n, dtype=int)
    if len(fpos):
        layer[fpos] = pareto_layers(T.S[fpos], 3)
    front = [int(p) for p in fpos if layer[p] == 1]
    sel_args = dict(k=k, top=top, pareto_pos=front if use_pareto else (), k_near=k_near, near_max=near_max, target=target)
    S = select(T, sc, **sel_args)
    sens = sensitivity(T, w, dict(sel_args, k=S["K"], target=None))
    conv = convergence(T, w, sel_args)
    feas_by, _ = niche_lists(T, sc)
    nrank = {p: (r + 1, len(v)) for v in feas_by.values() for r, p in enumerate(v)}
    extra = []
    for p in range(T.n):
        nr = nrank.get(p, (None, None))
        sr = S["sel"].get(p)
        extra.append([rank.get(p), nr[0], nr[1], int(layer[p]) if T.feas[p] else None, int(bool(sr)), "; ".join(sr) if sr else None])
    os.makedirs(out_dir, exist_ok=True)
    os.makedirs(csv_dir or out_dir, exist_ok=True)
    p_csv = os.path.join(csv_dir or out_dir, f"funnel_l0{sfx}.csv" + (".gz" if csv_gz else ""))
    p_sel = os.path.join(out_dir, f"funnel_l0_select{sfx}.json")
    p_res = os.path.join(out_dir, f"funnel_l0_results{sfx}.json")
    write_csv(p_csv, T, extra)
    soft = dict(method=ref["method"], software=ref["software"], date=DATE)
    # файл отбора
    sel_rows = []
    for p, why in S["sel"].items():
        nr = nrank.get(p, (None, None))
        sel_rows.append(dict(idx=int(T.idx[p]), niche=T.niche[p], arch=T.arch[p], niche_rank=nr[0], niche_size=nr[1], rank=rank.get(p),
                             score=None if np.isnan(sc[p]) else round(float(sc[p]), 4), feasible=bool(T.feas[p]),
                             violation=float(T.viol[p]), eta_el=T.get(p, "eta_el"), reasons=why))
    sel_rows.sort(key=lambda e: (not e["feasible"], e["rank"] or 10 ** 9, e["violation"], e["idx"]))
    covered = {e["niche"] for e in sel_rows if e["feasible"]}
    sel_summary = dict(n_selected=len(sel_rows), n_feasible=sum(e["feasible"] for e in sel_rows),
                       n_near=sum(not e["feasible"] for e in sel_rows), K_per_niche=S["K"], K_max_in_niche=S["kmax"],
                       target=target, target_exceeded=S["target_exceeded"], feasible_niches=len(S["feasible_niches"]),
                       feasible_niches_covered=len(covered & set(S["feasible_niches"])),
                       lost_feasible_niches=sorted(set(S["feasible_niches"]) - covered),
                       near_candidates=S["near_candidates"], archetypes=dict(Counter(e["arch"] for e in sel_rows).most_common()),
                       reasons=dict(Counter(x.split(":")[0] for e in sel_rows for x in e["reasons"]).most_common()))
    sel_doc = dict(title=f"Воронка, уровень 0 -> уровень 1: отбор {kv}", **soft, kv=kv, n=n, seed=seed, fingerprint=ref["fingerprint"],
                   code_sha=ref["code_sha"], weights=w,
                   params=dict(k=k, top=top, pareto=use_pareto, k_near=k_near, near_max=near_max, target=target),
                   how_to_rebuild="funnel_space.concept(idx, n, kv, seed) с закреплёнными версиями (fingerprint должен совпасть)",
                   summary=sel_summary, selected=sel_rows)
    with open(p_sel, "w", encoding="utf-8") as fh:
        json.dump(sel_doc, fh, ensure_ascii=False, indent=1)
    # результаты
    detail = detail_top if checks["fingerprint_local"]["match"] else 0
    archive = []
    for nk, pos in sorted(niche_lists(T, sc)[1].items()):
        fp = [p for p in pos if T.feas[p]]
        e = dict(niche=nk, n=len(pos), feasible=len(fp))
        if fp:
            b = max(fp, key=lambda p: (sc[p], -T.idx[p]))
            e.update(best_idx=int(T.idx[b]), best_score=round(float(sc[b]), 4), best_eta_el=T.get(b, "eta_el"))
        else:
            hz = [p for p in pos if T.has[p]]
            if hz:
                b = min(hz, key=lambda p: (T.viol[p], T.idx[p]))
                e.update(nearest_idx=int(T.idx[b]), nearest_violation=round(float(T.viol[b]), 4))
        archive.append(e)
    archive.sort(key=lambda e: (-e["feasible"], -(e.get("best_score") or 0), e.get("nearest_violation", 9e9), e["niche"]))
    res = dict(title=f"Воронка двигателя-заменителя, уровень 0: рейтинг, ниши, отбор на уровень 1 - {kv}, N = {n}",
               **soft, status=(extra_results or {}).pop("status", "слияние частей уровня 0"), kv=kv, target=ref["target"], n=n,
               seed=seed, shards=shards, inputs=inputs(checks), calib_ns=ref["calib"], merge_checks=checks,
               timing=dict(calls=calls, merge_s=None,
                           s_per_concept_core=round(float(np.mean([c["s_per_concept_core"] for c in calls if c.get("s_per_concept_core")])), 5)
                           if calls else None),
               weights=dict(owner=w, names=L0.CRIT_NAMES, split_effmake=L0.SPLIT_EFFMAKE, split_etasize=L0.SPLIT_ETASIZE,
                            sum=round(sum(w.values()), 12), basis="веса - решение владельца 14.11.2026; деление 4-го критерия - оценка"),
               normalization=normalization(T, w), summary=summary(T, sc),
               top=top_entries(T, sc, order, sens["stability"], layer, n, kv, seed, TOP_LIST, detail=bool(detail)),
               pareto=dict(criteria=L0.CRIT, n_front=len(front), layers={str(L): int((layer[fpos] == L).sum()) for L in (1, 2, 3)},
                           top_on_front=sum(1 for p in order[:TOP_LIST] if layer[p] == 1),
                           front=[dict(idx=int(T.idx[p]), arch=T.arch[p], niche=T.niche[p], rank=rank[p], score=round(float(sc[p]), 4),
                                       **{L0.CRIT_COL[c]: T.get(p, L0.CRIT_COL[c]) for c in L0.CRIT}, eta_el=T.get(p, "eta_el"))
                                  for p in sorted(front, key=lambda p: rank[p])[:300]]),
               sensitivity=dict(oat=sens["oat"], random=sens["random"],
                                top_stable=[dict(idx=int(T.idx[p]), rank=rank[p], **v) for p, v in sens["stability"].items()]),
               niches=dict(axes=L0.NICHE_AXES, n_seen=len(archive), n_feasible=sum(1 for e in archive if e["feasible"]),
                           archive=archive),
               selection=dict(file=os.path.basename(p_sel), **sel_summary),
               convergence=dict(rule=f"топ-{TOP_CONV} и допустимые ниши на первых N/4, N/2, N концепциях; концепция i одна и та же "
                                     "при любом N (префикс последовательности Соболя), поэтому первые N/2 - это выборка объёма N/2",
                                rows=conv),
               cycle_check=cycle_check(T, n, kv, seed, cycle_k) if (cycle_k and checks["fingerprint_local"]["match"]) else None,
               files=dict(csv=os.path.basename(p_csv), select=os.path.basename(p_sel)))
    csv_mb = os.path.getsize(p_csv) / 1e6
    res["files"]["csv_detail"] = dict(
        bytes=os.path.getsize(p_csv), sha256=hashlib.sha256(open(p_csv, "rb").read()).hexdigest(),
        outside_out_dir=bool(csv_dir) and os.path.abspath(csv_dir) != os.path.abspath(out_dir),
        note=(f"CSV {csv_mb:.1f} МБ - в ветку не класть (предел {CSV_BRANCH_MAX_MB:g} МБ): лежит в --csv-dir (S3 или локально у сессии с "
              "ключами); в ветку - только results и select JSON" if csv_mb > CSV_BRANCH_MAX_MB else
              f"CSV {csv_mb:.2f} МБ - можно в ветку (предел {CSV_BRANCH_MAX_MB:g} МБ)"))
    if csv_mb > CSV_BRANCH_MAX_MB and not res["files"]["csv_detail"]["outside_out_dir"]:
        print(f"ВНИМАНИЕ: {p_csv} - {csv_mb:.1f} МБ в папке результатов; в ветку не класть (нужен --csv-dir вне ветки)", file=sys.stderr)
    trad = [p for p in fpos if T.arch[p] == "традиционная"]
    tgt = ref["target"] or {}
    if tgt.get("engine_kg0") and trad:
        mk = T.num("mass_kg")[trad]
        kgkw = float(np.median(mk)) / tgt["P_shaft_kW"]
        res["mass_reference"] = dict(model_traditional_median_kg=round(float(np.median(mk))), engine_kg0=tgt["engine_kg0"],
                                     ratio=round(float(np.median(mk)) / tgt["engine_kg0"], 3),
                                     model_kg_per_kW=round(kgkw, 3), engine_kg0_per_kW=round(tgt["engine_kg0"] / tgt["P_shaft_kW"], 3),
                                     note=tgt.get("engine_kg0_note", "") +
                                     "; масса модели - оценка по оболочкам и колёсам без вспомогательных систем и занижена против "
                                     "исходного двигателя (см. ratio). Массы модулей по тем же формулам ВХОДЯТ в баллы: service_model.s_mass "
                                     "(самый тяжёлый модуль, 0,30 балла монтажа) и подъём plugplay_model (модуль против массы исходного) - "
                                     "заниженная абсолютная масса завышает эти два балла одинаково для всех концепций; сравнение концепций "
                                     "между собой - по относительной массе (исправление скептика физики: прежде было написано «в баллы не входит»)")
    if extra_results:
        res.update(extra_results)
    res["timing"]["merge_s"] = round(time.time() - t0, 2)
    with open(p_res, "w", encoding="utf-8") as fh:
        json.dump(res, fh, ensure_ascii=False, indent=1)
    return dict(results=p_res, csv=p_csv, select=p_sel, res=res, sel=sel_doc, table=T)


# =========================================================================== короткая локальная проба
PROBE = dict(kv="KV-3", n=256, shards=10, workers=4, top=40, target=150, k=3, k_near=1, near_max=0.25)


def probe(out_dir=HERE, scale=False, **kw):
    """короткая проба (владелец: локально уровень 0 - не больше 300 концепций): части funnel_l0 в пуле из 4 процессов ->
    слияние; та же выборка одной частью в одном процессе (объединение частей = расчёт целиком, независимость от числа
    процессов); замер времени на концепцию в одном процессе"""
    P = dict(PROBE, **kw)
    if P["n"] > 300:
        raise ValueError("локальная проба уровня 0 - не больше 300 концепций (решение владельца 14.11.2026)")
    tmp = tempfile.mkdtemp(prefix="funnel_l0_probe_")
    try:
        t0 = time.time()
        paths, _ = L0.run_parts(P["kv"], P["n"], list(range(P["shards"])), P["shards"], tmp, workers=P["workers"], gz=True)
        t_par = time.time() - t0
        gz_bytes = sum(os.path.getsize(x) for x in paths)
        whole = os.path.join(tmp, "whole")
        t1 = time.time()
        L0.run_parts(P["kv"], P["n"], [0], 1, whole, workers=1)
        t_ser = time.time() - t1
        a = merge_parts(find_parts([tmp], P["kv"], P["n"], P["shards"]), P["kv"], P["n"], P["shards"])[1]
        b = merge_parts(find_parts([whole], P["kv"], P["n"], 1), P["kv"], P["n"], 1)[1]
        same = json.dumps(a, ensure_ascii=False) == json.dumps(b, ensure_ascii=False)
        extra = dict(status=f"короткая локальная проба (4 общих ядра): {P['kv']}, N = {P['n']}, {P['shards']} частей; отладочный "
                            "(N = 2000) и полный прогоны - задачи AWS jobs_l0_debug.json и jobs_l0.json",
                     probe=dict(params=P, sharding_check=dict(parts_union_equals_whole=same, n=P["n"], shards=P["shards"],
                                                              workers_parts=P["workers"], workers_whole=1),
                                timing=dict(parts_wall_s=round(t_par, 2), whole_serial_wall_s=round(t_ser, 2),
                                            s_per_concept_serial=round(t_ser / P["n"], 5),
                                            parts_bytes_gz=gz_bytes, whole_bytes_json=os.path.getsize(L0.part_path(whole, P["kv"], P["n"], 0, 1)),
                                            note="whole_serial - один процесс, включая запись файла; это и есть время на вариант")),
                     )
        extra["probe"]["part_bytes_per_concept"] = round(gz_bytes / P["n"], 1)
        if scale:
            whole_doc = load_part(L0.part_path(whole, P["kv"], P["n"], 0, 1))
            meta = {x: whole_doc[x] for x in whole_doc if x not in ("rows", "columns", "n_rows")}
            extra["merge_scale_check"] = scale_check(whole_doc["rows"], whole_doc["columns"], meta)
        r = run_merge(find_parts([tmp], P["kv"], P["n"], P["shards"]), P["kv"], P["n"], P["shards"], out_dir, "",
                      k=P["k"], top=P["top"], k_near=P["k_near"], near_max=P["near_max"], target=P["target"], extra_results=extra)
        return r
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


def scale_check(src_rows, columns, meta, n=200000, shards=30):
    """проверка памяти и времени слияния на настоящем объёме без счёта: синтетические части (строки пробы по кругу,
    индексы 0..n-1 настоящие, отпечаток SYNTHETIC) -> слияние отдельным процессом; пик памяти и время"""
    import subprocess
    tmp = tempfile.mkdtemp(prefix="funnel_l0_scale_")
    try:
        for k in range(shards):
            rows = []
            for i in FS.shard_indices(n, k, shards):
                r = list(src_rows[i % len(src_rows)]); r[0] = i; rows.append(r)
            doc = dict(meta, n=n, shards=shards, shard=k, fingerprint="SYNTHETIC0000000", columns=columns, rows=rows, n_rows=len(rows))
            with gzip.open(L0.part_path(tmp, meta["kv"], n, k, shards, gz=True), "wt", encoding="utf-8") as fh:
                json.dump(doc, fh, ensure_ascii=False, separators=(",", ":"))
        gz = sum(os.path.getsize(L0.part_path(tmp, meta["kv"], n, k, shards, gz=True)) for k in range(shards))
        rep = os.path.join(tmp, "rss.json")
        r = subprocess.run([sys.executable, os.path.abspath(__file__), "--kv", meta["kv"], "--n", str(n), "--shards", str(shards),
                            "--in", tmp, "--out-dir", os.path.join(tmp, "out"), "--suffix", "_synthetic", "--top", "1000",
                            "--target", "4000", "--csv-gz", "--detail-top", "0", "--cycle-check", "0", "--rss-report", rep],
                           capture_output=True, text=True)
        d = json.load(open(rep, encoding="utf-8")) if os.path.exists(rep) else {}
        csv_gz = os.path.join(tmp, "out", "funnel_l0_synthetic.csv.gz")
        return dict(n=n, shards=shards, exit_code=r.returncode, parts_gz_MB=round(gz / 1e6, 1),
                    csv_gz_MB=round(os.path.getsize(csv_gz) / 1e6, 1) if os.path.exists(csv_gz) else None,
                    merge_wall_s=d.get("wall_s"), merge_maxrss_MB=d.get("maxrss_MB"),
                    note="синтетические части из строк пробы (SYNTHETIC): проверка механики, памяти и времени слияния, не результат")
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--probe", action="store_true", help="короткая локальная проба (KV-3, N=256, 10 частей)")
    ap.add_argument("--kv", default="KV-3", choices=sorted(FS.KV))
    ap.add_argument("--n", type=int)
    ap.add_argument("--shards", type=int)
    ap.add_argument("--seed", type=int, default=L0.SEED)
    ap.add_argument("--in", dest="inp", nargs="+", default=["."], help="папки с частями (поиск рекурсивный)")
    ap.add_argument("--out-dir", default=".")
    ap.add_argument("--suffix", default="", help="добавка к именам файлов, например _KV-3_n200000")
    ap.add_argument("--weights", default="", help="parts=0.35,plugplay=0.25,service=0.20,effmake=0.20 (сумма 1)")
    ap.add_argument("--k", type=int, default=3, help="вариантов на допустимую нишу (если нет --target)")
    ap.add_argument("--top", type=int, default=40, help="общий топ")
    ap.add_argument("--target", type=int, default=None, help="целевой размер отбора: K подбирается наибольшим, чтобы уложиться")
    ap.add_argument("--no-pareto", action="store_true")
    ap.add_argument("--k-near", type=int, default=1, help="почти допустимых на нишу без допустимых (0 - не брать)")
    ap.add_argument("--near-max", type=float, default=0.25, help="наибольшее превышение почти допустимых")
    ap.add_argument("--csv-gz", action="store_true")
    ap.add_argument("--csv-dir", default=None, help="папка CSV вне ветки (S3, routes/gtu/platform/aws/runs/<run_id>/merged): "
                                                    f"CSV крупнее {CSV_BRANCH_MAX_MB:g} МБ в ветку не кладётся")
    ap.add_argument("--detail-top", type=int, default=TOP_LIST)
    ap.add_argument("--cycle-check", type=int, default=6)
    ap.add_argument("--scale-check", action="store_true", help="с --probe: слияние синтетических частей N=200000 (память, время)")
    ap.add_argument("--rss-report", default="", help=argparse.SUPPRESS)
    a = ap.parse_args(argv)
    t0 = time.time()
    if a.probe:
        r = probe(scale=a.scale_check)
    else:
        if not (a.n and a.shards):
            ap.error("нужны --n и --shards")
        try:
            r = run_merge(find_parts(a.inp, a.kv, a.n, a.shards), a.kv, a.n, a.shards, a.out_dir, a.suffix, parse_weights(a.weights),
                          a.seed, a.k, a.top, not a.no_pareto, a.k_near, a.near_max, a.target, a.csv_gz, a.detail_top, a.cycle_check,
                          csv_dir=a.csv_dir)
        except MergeError as e:
            print(f"слияние отклонено: {e}", file=sys.stderr)
            return 2
    s = r["res"]["summary"]; sl = r["res"]["selection"]
    print(f"{r['res']['kv']} N={r['res']['n']}: допустимых {s['feasible']} ({s['feasible_share']:.2%}), ниш допустимых "
          f"{sl['feasible_niches']}, отбор {sl['n_selected']} (допустимых {sl['n_feasible']}, почти {sl['n_near']}, K={sl['K_per_niche']}); "
          f"фронт Парето {r['res']['pareto']['n_front']}; -> {r['results']}, {r['csv']}, {r['select']}")
    if a.rss_report:
        import resource
        with open(a.rss_report, "w", encoding="utf-8") as fh:
            json.dump(dict(maxrss_MB=round(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024, 1),
                           wall_s=round(time.time() - t0, 1)), fh)
    return 0


if __name__ == "__main__":
    sys.exit(main())
