#!/usr/bin/env python3
"""MAP-Elites по пространству сквозных концепций «топливо -> вал»: ниши (топология цикла × схема камеры × архитектура ×
изготовление), цель - «КПД на тонну», в каждой нише - фронт Парето «КПД - масса»; острова по машинам, слияние архивов
с контролем полноты, сводка.

Метод - MAP-Elites (Mouret & Clune 2015) с операторами evo/core/archive (гауссова мутация, iso+line, SBX, полиномиальная
мутация) плюс сброс 1-2 генов (переход в соседнюю нишу). Ограничения: в нише хранятся (а) лучшая допустимая концепция по
«КПД на тонну», (б) фронт Парето допустимых по (КПД на клеммах max, масса min) - не больше FRONT_CAP точек, (в) лучшая
недопустимая по сумме нарушений (родитель для ниш, где допустимых ещё нет). Детали и Plug&Play - вторичные цели: хранятся
у каждой точки, «невозможно» в Plug&Play - отказ (process_eval), отбор в сводке - и без ограничений, и с ними.

Остров = машина: один главный процесс ведёт архив, концепции пакета считаются пулом на все ядра ($NPROC); порядок результатов
пула фиксирован, поэтому итог острова не зависит от числа процессов. Начальная часть острова i - концепции генератора
по индексу [i·n_init, (i+1)·n_init) (process_space.genomes), дальше - потомки элит; сид острова - evo/core/islands.island_seed.
Острова части - islands_of_part (остров i -> часть i mod P). Слияние: все острова на месте, одна конфигурация (отпечаток),
число оценок острова = заданному, без повторов; иначе отказ (код 2), если не разрешено явно.

  python3 process_elites.py part  --kv KV-3 --n-islands 2 --n-parts 2 --part 0 --evals 2500 --workers 4 --out DIR
  python3 process_elites.py merge --kv KV-3 --n-islands 2 --evals 2500 --in DIR [--in DIR2 ...] --out-json F [--archive-out F.gz]
  python3 process_elites.py mini  --kv KV-3 --n 5000 --islands 2 --workers 4 --out DIR   # части по очереди + слияние (локально)
"""
import argparse
import collections
import glob
import gzip
import json
import math
import multiprocessing as mp
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
GTU = os.path.dirname(HERE)
for _p in (HERE, GTU):
    if _p not in sys.path:
        sys.path.insert(0, _p)
import process_space as PS  # noqa: E402
from evo.core import archive as EA  # noqa: E402  (операторы изменчивости; файл не меняется)
from evo.core import islands as EI  # noqa: E402  (сид острова, деление на части, отпечаток конфигурации)

FRONT_CAP = 6                 # точек фронта «КПД - масса» в нише (оценка: крайние + лучший КПД на тонну + по зазорам массы)
P_FEAS_PARENT = 0.8           # доля родителей из допустимых элит, пока есть недопустимые ниши без допустимых (оценка)
OPS = (("gauss", 0.40), ("iso_line", 0.30), ("reset", 0.20), ("sbx", 0.10))
SIGMA = 0.06                  # гауссова мутация в долях окна гена
TOURNAMENT = 3                # родитель из недопустимых ниш - лучший по нарушениям из трёх случайных (оценка)
VERSION = "process_elites v1 (15.11.2026)"
REC_KEYS = ("uid", "fit", "eff", "eff_shaft", "mass", "mass_lo", "mass_hi", "kwkg", "parts", "pnp", "trl", "ktrl", "fit_trl",
            "viol", "rej", "label", "G", "n_rpm")


# =========================================================================== записи и архив ниш
def compact(rec, genome):
    """запись архива из записи process_eval (геном хранится: концепция восстанавливается process_space.decode)"""
    return dict(uid=rec["uid"], g=[round(float(v), 6) for v in genome], fit=rec.get("fitness"), eff=rec.get("eff_el"),
                eff_shaft=rec.get("eff_shaft"), mass=rec.get("mass_t"), mass_lo=rec.get("mass_lo_t"), mass_hi=rec.get("mass_hi_t"),
                kwkg=rec.get("kW_per_kg"), parts=rec.get("parts"), pnp=rec.get("pnp"), trl=rec.get("trl"), ktrl=rec.get("k_trl"),
                fit_trl=rec.get("fit_trl"), viol=rec.get("violation"), rej=(rec.get("reject") or [])[:4], label=rec.get("label"),
                G=rec.get("G_air"), n_rpm=rec.get("n_rpm"), feas=bool(rec.get("feasible")), stage=rec.get("stage"))


def _better_best(a, b):
    """a лучше b (допустимые): больший КПД на тонну, при равенстве - меньший uid"""
    return a["fit"] > b["fit"] or (a["fit"] == b["fit"] and a["uid"] < b["uid"])


def _better_infeas(a, b):
    va = a["viol"] if a["viol"] is not None else 1e9
    vb = b["viol"] if b["viol"] is not None else 1e9
    return va < vb or (va == vb and a["uid"] < b["uid"])


def pareto(points):
    """недоминируемые по (eff max, mass min); порядок - по росту массы, при равенстве - uid"""
    pts = sorted(points, key=lambda p: (p["mass"], -p["eff"], p["uid"]))
    out, best = [], -1e300
    for p in pts:
        if p["eff"] > best + 1e-12:
            out.append(p)
            best = p["eff"]
    return out


def prune(front, cap=FRONT_CAP):
    """не больше cap точек: лёгкая, лучший КПД, лучший КПД на тонну, затем - по наибольшему зазору массы (детерминировано)"""
    if len(front) <= cap:
        return front
    keep = {front[0]["uid"], front[-1]["uid"], sorted(front, key=lambda p: (-p["fit"], p["uid"]))[0]["uid"]}
    while len(keep) < cap:
        kept = sorted([p for p in front if p["uid"] in keep], key=lambda p: p["mass"])
        gaps = []
        for a, b in zip(kept, kept[1:]):
            inner = [p for p in front if a["mass"] < p["mass"] < b["mass"] and p["uid"] not in keep]
            if inner:
                mid = 0.5 * (a["mass"] + b["mass"])
                c = min(inner, key=lambda p: (abs(p["mass"] - mid), p["uid"]))
                gaps.append((b["mass"] - a["mass"], c["uid"]))
        if not gaps:
            break
        keep.add(max(gaps)[1])
    return [p for p in front if p["uid"] in keep]


class NicheArchive:
    def __init__(self):
        self.cells = {}

    def add(self, c, cell):
        cell = tuple(int(v) for v in cell)
        e = self.cells.setdefault(cell, dict(best=None, front=[], infeas=None))
        if c["feas"] and c["fit"] is not None and math.isfinite(c["fit"]):
            changed = False
            if e["best"] is None or _better_best(c, e["best"]):
                e["best"] = c
                changed = True
            f = pareto(e["front"] + [c])
            if any(p["uid"] == c["uid"] for p in f):
                e["front"] = prune(f)
                changed = True
            return changed
        if e["infeas"] is None or _better_infeas(c, e["infeas"]):
            e["infeas"] = c
            return True
        return False

    def merge(self, cells):
        """cells: {cell: entry} другого архива; правила те же - итог не зависит от порядка (фронт - объединение и отсечение)"""
        for cell, o in cells.items():
            e = self.cells.setdefault(tuple(cell), dict(best=None, front=[], infeas=None))
            if o["best"] is not None and (e["best"] is None or _better_best(o["best"], e["best"])):
                e["best"] = o["best"]
            if o["front"]:
                seen = {p["uid"]: p for p in e["front"] + o["front"]}
                e["front"] = prune(pareto(list(seen.values())))
            if o["infeas"] is not None and (e["infeas"] is None or _better_infeas(o["infeas"], e["infeas"])):
                e["infeas"] = o["infeas"]

    def feasible_cells(self):
        return sorted(c for c, e in self.cells.items() if e["best"] is not None)

    def infeasible_only(self):
        return sorted(c for c, e in self.cells.items() if e["best"] is None and e["infeas"] is not None)

    def n_feasible(self):
        return sum(1 for e in self.cells.values() if e["best"] is not None)

    def qd(self):
        return float(sum(e["best"]["fit"] for e in self.cells.values() if e["best"] is not None))

    def best(self):
        bs = [e["best"] for e in self.cells.values() if e["best"] is not None]
        return sorted(bs, key=lambda b: (-b["fit"], b["uid"]))[0] if bs else None

    def to_list(self):
        return [dict(cell=list(c), best=e["best"], front=e["front"], infeas=e["infeas"]) for c, e in sorted(self.cells.items())]

    @classmethod
    def from_list(cls, lst):
        a = cls()
        for d in lst:
            a.cells[tuple(d["cell"])] = dict(best=d["best"], front=d["front"] or [], infeas=d["infeas"])
        return a

    def parent(self, rng, fc=None, ic=None):
        """родитель: из допустимых ниш (лучший или точка фронта) с вероятностью P_FEAS_PARENT, иначе - турнир недопустимых.
        fc, ic - списки ниш (передаются один раз на пакет: сортировка на каждого потомка дорога при тысячах ниш)"""
        fc = self.feasible_cells() if fc is None else fc
        ic = self.infeasible_only() if ic is None else ic
        if fc and (not ic or rng.random() < P_FEAS_PARENT):
            e = self.cells[fc[int(rng.integers(len(fc)))]]
            pool = [e["best"]] + [p for p in e["front"] if p["uid"] != e["best"]["uid"]]
            return np.asarray(pool[int(rng.integers(len(pool)))]["g"], float)
        cs = ic or sorted(self.cells)
        # недопустимые: турнир из TOURNAMENT ниш по сумме нарушений - давление к допустимости, пока допустимых мало
        pick = [cs[int(rng.integers(len(cs)))] for _ in range(TOURNAMENT)]
        e = min((self.cells[c] for c in pick), key=lambda e: ((e["infeas"] or e["best"])["viol"] or 0.0, (e["infeas"] or e["best"])["uid"]))
        p = e["infeas"] or e["best"]
        return np.asarray(p["g"], float)


def children(rng, arch, k):
    names = [o for o, _ in OPS]
    w = np.array([p for _, p in OPS], float)
    w /= w.sum()
    out = []
    fc, ic = arch.feasible_cells(), arch.infeasible_only()
    par = lambda: arch.parent(rng, fc, ic)  # noqa: E731
    while len(out) < k:
        op = names[int(rng.choice(len(names), p=w))]
        x = par()
        if op == "gauss":
            y = EA.gaussian(rng, x, SIGMA)
        elif op == "iso_line":
            y = EA.iso_line(rng, x, par(), s_iso=0.02, s_line=0.25)
        elif op == "reset":
            y = x.copy()
            for j in rng.choice(PS.D, size=int(rng.integers(1, 3)), replace=False):
                y[j] = rng.random()
        else:
            y, _ = EA.sbx(rng, x, par())
            y = EA.poly_mutation(rng, y)
        out.append(np.clip(y, 0.0, 1.0 - 1e-12))
    return np.stack(out)


# =========================================================================== счёт (пул процессов)
_PE = None


def _init_worker():
    global _PE
    import process_eval
    _PE = process_eval


def _eval_one(args):
    x, kv, uid = args
    global _PE
    if _PE is None:
        _init_worker()
    rec = _PE.evaluate(np.asarray(x, float), kv, uid=uid)
    return compact(rec, x), rec["cell"], rec.get("cpu_s", 0.0), rec.get("stage"), rec.get("reject") or []


class Evaluator:
    def __init__(self, workers):
        self.workers = max(1, int(workers))
        self.pool = None
        if self.workers > 1:
            ctx = mp.get_context("fork") if "fork" in mp.get_all_start_methods() else mp.get_context()
            self.pool = ctx.Pool(self.workers, initializer=_init_worker)

    def map(self, G, kv, uids):
        jobs = [(list(map(float, g)), kv, u) for g, u in zip(G, uids)]
        if self.pool is None:
            return [_eval_one(j) for j in jobs]
        cs = max(1, len(jobs) // (self.workers * 4))
        return self.pool.map(_eval_one, jobs, chunksize=cs)

    def close(self):
        if self.pool is not None:
            self.pool.close()
            self.pool.join()


def seed_genomes(kv, parents=None):
    """начальные родители острова: опорные концепции блоков (process_space.reference_seeds) + геномы из файла parents
    (слитый архив прошлого этапа или выбор владельца: JSON/JSON.gz со списком archive[].best/front[] или genomes[])"""
    out = [x for _, x in PS.reference_seeds(kv)]
    if parents:
        op = gzip.open if parents.endswith(".gz") else open
        with op(parents, "rt", encoding="utf-8") as f:
            d = json.load(f)
        for e in d.get("archive", []):
            for p in ([e.get("best")] + (e.get("front") or [])):
                if p and p.get("g"):
                    out.append(np.asarray(p["g"], float))
        out += [np.asarray(g, float) for g in d.get("genomes", [])]
    uniq, seen = [], set()
    for x in out:
        k = tuple(np.round(x, 6))
        if k not in seen:
            seen.add(k)
            uniq.append(np.clip(np.asarray(x, float), 0.0, 1.0 - 1e-12))
    return uniq


def run_config(kv, n_islands, evals, init_frac, batch, seed, parents=None):
    seeds = seed_genomes(kv, parents)
    sh = EI.config_hash([list(np.round(x, 6)) for x in seeds])
    return dict(version=VERSION, kv=kv, n_islands=int(n_islands), evals=int(evals), init_frac=float(init_frac), batch=int(batch),
                seed=int(seed), space=PS.space_hash(), front_cap=FRONT_CAP, ops=OPS, sigma=SIGMA, p_feas=P_FEAS_PARENT,
                tournament=TOURNAMENT, seeds=len(seeds), seeds_hash=sh)


def run_island(island, kv, evals, init_frac, batch, seed, ev, log=None, seeds=()):
    rng = EI.rng_for(seed, island, salt=kv)
    arch = NicheArchive()
    seeds = list(seeds)[:evals]
    n_init = int(round(init_frac * evals))
    done = 0
    stages, rejects = collections.Counter(), collections.Counter()
    cpu = collections.defaultdict(float)
    hist = []
    t0 = time.time()
    feas_n = 0
    while done < evals:
        k = min(batch, evals - done)
        if done < len(seeds):
            k = min(k, len(seeds) - done)
            G = np.stack(seeds[done:done + k])
            src = "опорные"
        elif done < len(seeds) + n_init:
            k = min(k, len(seeds) + n_init - done)
            G = PS.genomes(island * n_init + done - len(seeds), k, kv, seed)
            src = "индекс"
        elif arch.cells:
            G = children(rng, arch, k)
            src = "потомки"
        else:
            G = rng.random((k, PS.D))
            src = "случайные"
        uids = [f"{kv}:{island:03d}:{done + j:08d}" for j in range(k)]
        for g, (c, cell, cs, st, rj) in zip(G, ev.map(G, kv, uids)):
            arch.add(c, cell)
            stages[st] += 1
            cpu[st] += cs
            feas_n += c["feas"]
            for code in rj:
                rejects[code] += 1
        done += k
        b = arch.best()
        hist.append(dict(evals=done, src=src, feasible_cells=arch.n_feasible(), cells=len(arch.cells), feasible=feas_n,
                         best=(round(b["fit"], 4) if b else None), qd=round(arch.qd(), 2), wall_s=round(time.time() - t0, 1)))
        if log:
            log(f"[{kv} остров {island}] {done}/{evals} ({src}): ниш с допустимыми {arch.n_feasible()}, допустимых {feas_n}, "
                f"лучший КПД/т {b['fit'] if b else None}, {time.time() - t0:.0f} с")
    return dict(kind="process_island", island=int(island), kv=kv, evals=done, n_init=n_init, feasible=feas_n,
                stages=dict(stages), cpu_s_by_stage={k: round(v, 3) for k, v in cpu.items()}, rejects=dict(rejects.most_common(40)),
                history=hist, wall_s=round(time.time() - t0, 1), archive=arch.to_list())


# =========================================================================== часть и слияние
def part_file(out, kv, part):
    return os.path.join(out, f"process_{kv}_part{part:03d}.json.gz")


def run_part(kv, n_islands, n_parts, part, evals, init_frac=0.3, batch=512, seed=PS.SEED, workers=1, out=".", log=print, parents=None):
    cfg = run_config(kv, n_islands, evals, init_frac, batch, seed, parents)
    h = EI.config_hash(cfg)
    ids = EI.islands_of_part(n_islands, n_parts, part)
    seeds = seed_genomes(kv, parents)
    ev = Evaluator(workers)
    t0 = time.time()
    try:
        isl = [run_island(i, kv, evals, init_frac, batch, seed, ev, log=log, seeds=seeds) for i in ids]
    finally:
        ev.close()
    arch = NicheArchive()
    for d in isl:
        arch.merge(NicheArchive.from_list(d["archive"]).cells)
    doc = dict(kind="process_part", version=VERSION, kv=kv, part=part, n_parts=n_parts, islands=ids, config=cfg, config_hash=h,
               island_stats={str(d["island"]): {k: d[k] for k in ("evals", "n_init", "feasible", "stages", "cpu_s_by_stage", "rejects",
                                                                  "history", "wall_s")} for d in isl},
               workers=int(workers), wall_s=round(time.time() - t0, 1), archive=arch.to_list())
    os.makedirs(out, exist_ok=True)
    fn = part_file(out, kv, part)
    with gzip.open(fn, "wt", encoding="utf-8") as f:
        json.dump(doc, f, ensure_ascii=False, separators=(",", ":"))
    if log:
        log(f"[{kv} часть {part}/{n_parts}] острова {ids}: ниш с допустимыми {arch.n_feasible()}, {time.time() - t0:.0f} с -> {fn}")
    return fn, doc


def merge_parts(paths, kv, n_islands, evals, allow_partial=False, init_frac=0.3, batch=512, seed=PS.SEED, parents=None):
    """слияние частей с контролем полноты: все острова 0..n_islands-1, одна конфигурация, оценок острова = evals, без повторов"""
    cfg = run_config(kv, n_islands, evals, init_frac, batch, seed, parents)
    h = EI.config_hash(cfg)
    seen, dup, foreign, badh, short = {}, [], [], [], []
    arch = NicheArchive()
    parts = []
    for p in sorted(paths):
        with gzip.open(p, "rt", encoding="utf-8") as f:
            d = json.load(f)
        if d.get("kind") != "process_part" or d.get("kv") != kv:
            continue
        if d.get("config_hash") != h:
            badh.append(os.path.basename(p))
            continue
        parts.append(d["part"])
        for i in d["islands"]:
            st = d["island_stats"][str(i)]
            if i in seen:
                dup.append(i)
                continue
            if not 0 <= i < n_islands:
                foreign.append(i)
                continue
            if st["evals"] != evals:
                short.append(dict(island=i, evals=st["evals"]))
            seen[i] = st
        arch.merge(NicheArchive.from_list(d["archive"]).cells)
    missing = [i for i in range(n_islands) if i not in seen]
    rep = dict(expected_islands=n_islands, present=len(seen), missing=missing, duplicates=sorted(set(dup)), foreign=sorted(set(foreign)),
               config_mismatch_files=badh, short_islands=short, parts=sorted(parts), config_hash=h,
               evals_total=int(sum(s["evals"] for s in seen.values())), evals_expected=int(n_islands * evals))
    rep["complete"] = not (missing or dup or foreign or badh or short)
    rep["status"] = "полное" if rep["complete"] else (f"НЕПОЛНОЕ: нет островов {missing[:20]}, повторы {sorted(set(dup))[:20]}, "
                                                      f"чужая конфигурация {badh[:5]}, недосчитаны {short[:5]}")
    if not rep["complete"] and not allow_partial:
        raise SystemExit(f"слияние {kv}: {rep['status']} (разрешить - --allow-partial)")
    return arch, seen, rep, cfg


# =========================================================================== сводка
def _row(c, cell=None):
    r = {k: c.get(k) for k in REC_KEYS}
    if cell is not None:
        r["cell"] = list(cell)
        r["niche"] = PS.cell_label(cell)
    return r


def reference_parts(kv):
    """эффективные детали опорной точки frame3 (блоки) - порог «не больше деталей, чем у frame3» для отбора с ограничениями"""
    try:
        import blocks
        from blocks import reference as REF
        r = blocks.evaluate(dict(REF.frame3_concept(), kv=kv))
        return r["parts"]["n_eff_total"]
    except Exception:  # noqa: BLE001
        return None


PNP_MIN = 7.0     # Plug&Play не ниже 7 из 10 для отбора с ограничениями (оценка; 10 - всё сохраняется)


def summarize(arch, stats, rep, cfg, kv, top=25, detail_top=8, t_run=None):
    import process_eval as PE
    fc = arch.feasible_cells()
    bests = [(c, arch.cells[c]["best"]) for c in fc]
    allpts = [(c, p) for c in fc for p in arch.cells[c]["front"]]
    n_ev = int(sum(s["evals"] for s in stats.values()))
    n_feas = int(sum(s["feasible"] for s in stats.values()))
    stages, rejects, cpu = collections.Counter(), collections.Counter(), collections.Counter()
    for s in stats.values():
        stages.update(s["stages"])
        rejects.update(s["rejects"])
        cpu.update(s["cpu_s_by_stage"])
    ref_parts = reference_parts(kv)
    ordered = sorted(bests, key=lambda cb: (-cb[1]["fit"], cb[1]["uid"]))
    cons = [cb for cb in ordered if (cb[1]["pnp"] or 0) >= PNP_MIN and (ref_parts is None or (cb[1]["parts"] or 1e9) <= ref_parts)]
    by_trl = sorted(bests, key=lambda cb: (-(cb[1]["fit_trl"] or 0), cb[1]["uid"]))
    front = pareto([dict(p, _cell=c) for c, p in allpts])
    marg = {}
    for j, name in enumerate(PS.DESC_NAMES):
        d = {}
        for c, b in bests:
            lab = PS.DESC_LABELS[j][c[j]]
            if lab not in d or _better_best(b, d[lab][1]):
                d[lab] = (c, b)
        marg[name] = {lab: dict(fit=round(b["fit"], 3), eff=b["eff"], mass=b["mass"], niche=PS.cell_label(c), label=b["label"], uid=b["uid"])
                      for lab, (c, b) in sorted(d.items(), key=lambda kv_: -kv_[1][1]["fit"])}
    cov = {}
    for j, name in enumerate(PS.DESC_NAMES):
        used = sorted({c[j] for c in fc})
        cov[name] = dict(bins=PS.DESC_BINS[j], used=len(used), labels=[PS.DESC_LABELS[j][u] for u in used])
    # подробности лучших (повторная оценка с деталями)
    det = []
    seen_top = []
    for c, b in ordered[:detail_top]:
        rec = PE.evaluate(np.asarray(b["g"], float), kv, detail=True, uid=b["uid"])
        m = rec.get("mass") or {}
        bl = rec.get("blocks") or {}
        det.append(dict(uid=b["uid"], niche=PS.cell_label(c), label=rec.get("label"), eff_el=rec.get("eff_el"), eff_shaft=rec.get("eff_shaft"),
                        mass_t=rec.get("mass_t"), mass_band_t=[rec.get("mass_lo_t"), rec.get("mass_hi_t")], kW_per_kg=rec.get("kW_per_kg"),
                        kpd_per_t=rec.get("kpd_per_t"), kpd_per_t_range=rec.get("kpd_per_t_range"), engine_t=rec.get("engine_t"),
                        hx_t=rec.get("hx_t"), frame_t=rec.get("frame_t"), parts=rec.get("parts"), pnp=rec.get("pnp"), trl=rec.get("trl"),
                        interface=rec.get("interface"), n_rpm=rec.get("n_rpm"), L_engine=rec.get("L_engine"), D_case_max=rec.get("D_case_max"),
                        reproduced=abs((rec.get("kpd_per_t") or 0) - b["fit"]) < 1e-6,
                        mass_top_rows=[[x["name"], round(x["mass_kg"])] for x in (m.get("rows") or [])[:8]],
                        compressor=[{k: s.get(k) for k in ("type", "pi", "eta_p", "D2_m", "U2_ms", "mat")}
                                    for s in (bl.get("blocks", {}).get("compressor", {}).get("stages") or [])],
                        turbine=bl.get("blocks", {}).get("turbine", {}).get("geometry"), trl_items=(bl.get("trl") or {}).get("items", [])[:4],
                        warnings=rec.get("warnings"), mass_notes=rec.get("mass_notes")))
        seen_top.append(b["uid"])
    return dict(
        kv=kv, config=cfg, merge=rep, evals=n_ev, feasible=n_feas, feasible_share=round(n_feas / max(n_ev, 1), 5),
        stages=dict(stages), rejects_top=dict(rejects.most_common(20)),
        cpu_s_per_eval=dict(all=round(sum(cpu.values()) / max(n_ev, 1), 5),
                            by_stage={k: round(cpu[k] / max(stages[k], 1), 5) for k in stages}),
        niches=dict(total=PS.N_CELLS, touched=len(arch.cells), feasible=len(fc), coverage=round(len(fc) / PS.N_CELLS, 5),
                    qd_score=round(arch.qd(), 2), per_descriptor=cov, front_points=len(allpts)),
        best=_row(ordered[0][1], ordered[0][0]) if ordered else None,
        top=[_row(b, c) for c, b in ordered[:top]],
        top_constrained=dict(rule=f"Plug&Play >= {PNP_MIN:g} и эффективных деталей не больше, чем у опорной точки frame3 ({ref_parts})",
                             rows=[_row(b, c) for c, b in cons[:top]]),
        top_trl=dict(rule="КПД на тонну × (0,5 + 0,5·TRL/9) - множитель готовности владельца 15.11.2026",
                     rows=[_row(b, c) for c, b in by_trl[:top]]),
        pareto_global=dict(rule="недоминируемые по (КПД на клеммах max, масса min) среди фронтов всех ниш",
                           rows=[_row(p, p["_cell"]) for p in front][:60]),
        best_by_descriptor=marg,
        detail_top=det,
        history={str(i): s["history"][-1] for i, s in stats.items()},
        run_wall_s=t_run,
    )


def write_archive(arch, path, meta):
    with gzip.open(path, "wt", encoding="utf-8") as f:
        json.dump(dict(meta, archive=arch.to_list()), f, ensure_ascii=False, separators=(",", ":"))


# =========================================================================== интерфейс командной строки
def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)

    def common(sp):
        sp.add_argument("--kv", default="KV-3")
        sp.add_argument("--n-islands", type=int, required=False)
        sp.add_argument("--evals", type=int, required=False, help="оценок на остров")
        sp.add_argument("--init-frac", type=float, default=0.3)
        sp.add_argument("--batch", type=int, default=512)
        sp.add_argument("--seed", type=int, default=PS.SEED)
        sp.add_argument("--parents", help="начальные родители: слитый архив прошлого этапа / выбор владельца (JSON или JSON.gz)")

    sp = sub.add_parser("part", help="острова одной части (машины)")
    common(sp)
    sp.add_argument("--n-parts", type=int, required=True)
    sp.add_argument("--part", type=int, required=True)
    sp.add_argument("--workers", default="1")
    sp.add_argument("--out", required=True)
    sp = sub.add_parser("merge", help="слияние частей с контролем полноты и сводка")
    common(sp)
    sp.add_argument("--in", dest="inp", action="append", required=True, help="папка с частями (можно несколько; поиск вглубь)")
    sp.add_argument("--allow-partial", action="store_true")
    sp.add_argument("--out-json", required=True)
    sp.add_argument("--archive-out")
    sp.add_argument("--top", type=int, default=25)
    sp.add_argument("--detail-top", type=int, default=8)
    sp = sub.add_parser("mini", help="локальный прогон: части по очереди, слияние, сводка")
    common(sp)
    sp.add_argument("--n", type=int, required=True, help="всего оценок")
    sp.add_argument("--islands", type=int, default=2)
    sp.add_argument("--workers", default="1")
    sp.add_argument("--out", required=True)
    sp.add_argument("--out-json", required=True)
    a = ap.parse_args(argv)
    workers = lambda w: int(os.environ.get("NPROC", "1")) if str(w).startswith("$") else int(w)  # noqa: E731
    if a.cmd == "part":
        run_part(a.kv, a.n_islands, a.n_parts, a.part, a.evals, a.init_frac, a.batch, a.seed, workers(a.workers), a.out, parents=a.parents)
        return 0
    if a.cmd == "mini":
        ev_isl = a.n // a.islands
        t0 = time.time()
        for p in range(a.islands):
            run_part(a.kv, a.islands, a.islands, p, ev_isl, a.init_frac, a.batch, a.seed, workers(a.workers), a.out, parents=a.parents)
        a.inp, a.n_islands, a.evals, a.allow_partial, a.top, a.detail_top, a.archive_out = [a.out], a.islands, ev_isl, False, 25, 8, None
        t_run = round(time.time() - t0, 1)
    else:
        t_run = None
    paths = []
    for d in a.inp:
        paths += glob.glob(os.path.join(d, "**", f"process_{a.kv}_part*.json.gz"), recursive=True)
    arch, stats, rep, cfg = merge_parts(paths, a.kv, a.n_islands, a.evals, a.allow_partial, a.init_frac, a.batch, a.seed, a.parents)
    s = summarize(arch, stats, rep, cfg, a.kv, a.top, a.detail_top, t_run)
    with open(a.out_json, "w", encoding="utf-8") as f:
        json.dump(s, f, ensure_ascii=False, indent=1, default=str)
    if a.archive_out:
        write_archive(arch, a.archive_out, dict(kind="process_merged", kv=a.kv, merge=rep, config=cfg))
    b = s["best"]
    print(f"{a.kv}: оценок {s['evals']}, допустимых {s['feasible']} ({100 * s['feasible_share']:.2f} %), ниш с допустимыми "
          f"{s['niches']['feasible']} из {s['niches']['total']}, слияние: {rep['status']}; лучший КПД/т "
          f"{b['fit'] if b else None} ({b['niche'] if b else ''})")
    return 0 if rep["complete"] else 2


if __name__ == "__main__":
    sys.exit(main())
