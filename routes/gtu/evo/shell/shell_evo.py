#!/usr/bin/env python3
"""Поиск форм семейства «корпус-оболочка и рама»: MAP-Elites на островах → слияние с контролем полноты → NSGA-II
по фронту (масса, детали, ремонт, изготовимость, новизна против типовых форм OEM) → CMA-ES-доводка элит в своей ячейке
→ галерея → сводка. Вариант («оболочка» / «картридж») — параметр прогона: архивы вариантов не смешиваются.

Команды:
  run-part  --variant V --n-islands N --n-parts P --part k --evals E [--seed S] [--workers W] [--picks f.json] --out DIR
  merge     --variant V --n-islands N --in DIR [DIR …] --out merged_V.json [--allow-partial]
  nsga      --variant V --merged merged_V.json --pop 120 --gens 60 --out nsga_V.json
  polish    --variant V --merged merged_V.json --top 3 --evals 1500 --out polish_V.json
  gallery   --variant V --merged … [--nsga …] [--polish …] --out DIR   (рендер: pyvista, если есть, иначе matplotlib)
  report    --dir DIR --out evo_shell_results.json
  picks     --merged merged_V.json --uids u1,u2 --out picks.json      (выбор владельца → родители следующего прогона)

Острова: сид — sha256(seed | номер | вариант); часть — номер mod частей; смесь операторов — по номеру острова mod 4.
Код модели входит в хэш конфигурации: архивы, посчитанные разным кодом, не сливаются.
"""
import os

for _v in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")
import argparse  # noqa: E402
import glob  # noqa: E402
import hashlib  # noqa: E402
import json  # noqa: E402
import math  # noqa: E402
import sys  # noqa: E402
import tempfile  # noqa: E402
import time  # noqa: E402

import numpy as np  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
EVO = os.path.dirname(HERE)
sys.path.insert(0, HERE)
sys.path.insert(0, EVO)
from core import archive as A  # noqa: E402
from core import cmaes as CM  # noqa: E402
from core import features as FT  # noqa: E402
from core import islands as IS  # noqa: E402
from core import nsga2 as N2  # noqa: E402
import shell_model as SM  # noqa: E402

OUT = os.environ.get("EVO_OUT") or os.path.join(tempfile.gettempdir(), "evo_shell")
GRID = A.GridSpec(SM.DESC_NAMES, SM.DESC_EDGES, SM.DESC_LABELS)
ALIAS = {"shell": "оболочка", "cartridge": "картридж", "оболочка": "оболочка", "картридж": "картридж"}   # латиница — для команд AWS
OPS = [  # смесь операторов и шаг мутации острова (по номеру mod 4)
    dict(ops=(("iso_line", 0.7), ("gauss", 0.2), ("sbx", 0.1)), sigma=0.05),
    dict(ops=(("iso_line", 0.2), ("gauss", 0.7), ("sbx", 0.1)), sigma=0.10),
    dict(ops=(("iso_line", 0.3), ("gauss", 0.2), ("sbx", 0.5)), sigma=0.05),
    dict(ops=(("iso_line", 0.4), ("gauss", 0.4), ("sbx", 0.2)), sigma=0.03),
]


def jd(o):
    if isinstance(o, (np.floating,)):
        return float(o)
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, (np.bool_,)):
        return bool(o)
    if isinstance(o, np.ndarray):
        return o.tolist()
    raise TypeError(type(o))


def dump(obj, path):
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    tmp = path + ".tmp"
    with open(tmp, "w", encoding="utf-8") as f:
        json.dump(obj, f, ensure_ascii=False, default=jd)
    os.replace(tmp, path)


def load(path):
    with open(path, encoding="utf-8") as f:
        return json.load(f)


def code_hash():
    h = hashlib.sha256()
    for fn in (os.path.join(HERE, "shell_model.py"), os.path.join(EVO, "core", "archive.py"), os.path.join(EVO, "core", "shapes.py"),
               os.path.join(EVO, "core", "features.py")):
        with open(fn, "rb") as f:
            h.update(f.read())
    return h.hexdigest()[:16]


def config(variant, evals, batch, seed, n_islands, picks_digest=""):
    return dict(family="shell", variant=variant, evals_per_island=evals, batch=batch, seed=seed, n_islands=n_islands, dim=SM.DIM,
                genes=[g["name"] for g in SM.G], desc=SM.DESC_EDGES, w_fit=SM.W_FIT, code=code_hash(), picks=picks_digest)


def init_genomes(variant, picks):
    """начальные родители: опорные формы OEM (для сравнения в той же сетке) и выбор владельца"""
    G0 = [SM.encode(SM.reference_params(variant, k)) for k in SM.REF_KINDS]
    if picks:
        G0 += [np.asarray(g, float) for g in picks.get("genomes", [])]
    return G0


def run_island(island, variant, evals, batch, seed, cfg_hash, picks=None):
    t0 = time.time()
    rng = IS.rng_for(seed, island, salt=variant)
    ar = A.Archive(GRID)
    o = OPS[island % len(OPS)]

    def ev(Gm):
        return [SM.evaluate(x, variant) for x in Gm]
    hist = A.map_elites(ev, SM.DIM, ar, rng, evals, batch=batch, init=init_genomes(variant, picks), ops=o["ops"], sigma=o["sigma"],
                        uid_prefix=f"{'o' if variant == 'оболочка' else 'k'}{island:04d}-")
    step = max(1, len(hist) // 40)
    return dict(island=island, variant=variant, config_hash=cfg_hash, evals=evals, archive=ar.to_dict(), history=hist[::step] + hist[-1:],
                seconds=round(time.time() - t0, 1), ops=o)


def cmd_run_part(a):
    picks = load(a.picks) if a.picks else None
    pd = hashlib.sha256(json.dumps(picks, sort_keys=True).encode()).hexdigest()[:12] if picks else ""
    cfg = config(a.variant, a.evals, a.batch, a.seed, a.n_islands, pd)
    ch = IS.config_hash(cfg)
    ids = IS.islands_of_part(a.n_islands, a.n_parts, a.part)
    print(f"[part] вариант {a.variant}: часть {a.part}/{a.n_parts}, острова {ids[:8]}{'…' if len(ids) > 8 else ''} ({len(ids)}), "
          f"оценок на остров {a.evals}, конфигурация {ch}", flush=True)
    t0 = time.time()
    res = IS.run_part(run_island, a.n_islands, a.n_parts, a.part, workers=a.workers, variant=a.variant, evals=a.evals,
                      batch=a.batch, seed=a.seed, cfg_hash=ch, picks=picks)
    d = os.path.join(a.out, "islands", a.variant)
    for r in res:
        dump(r, os.path.join(d, f"island_{r['island']:04d}.json"))
    man = dict(variant=a.variant, part=a.part, n_parts=a.n_parts, n_islands=a.n_islands, islands=ids, config=cfg, config_hash=ch,
               evals=int(sum(r["evals"] for r in res)), seconds=round(time.time() - t0, 1),
               filled=[len(r["archive"]["elites"]) for r in res])
    dump(man, os.path.join(a.out, "islands", a.variant, f"part_{a.part:03d}_of_{a.n_parts:03d}.json"))
    print(f"[part] готово: {man['evals']} оценок за {man['seconds']} с; заполнено ячеек по островам {man['filled'][:8]}", flush=True)


def cmd_merge(a):
    files = []
    for d in a.inp:
        files += sorted(glob.glob(os.path.join(d, "**", "island_*.json"), recursive=True))
    isl = [load(f) for f in files]
    isl = [d for d in isl if d.get("variant") == a.variant]
    if not isl:
        raise SystemExit(f"нет островов варианта {a.variant} в {a.inp}")
    ch = isl[0]["config_hash"]
    merged, rep = IS.merge_archives(isl, range(a.n_islands), ch, allow_partial=a.allow_partial)
    F = np.array([e["meta"]["feats"] for e in merged.elites.values()])
    mu, sd = F.mean(axis=0), F.std(axis=0)
    out = dict(variant=a.variant, config_hash=ch, report=rep, archive=merged.to_dict(), feats_mu=mu.tolist(), feats_sd=sd.tolist(),
               files=len(files), n_feasible=int(sum(1 for e in merged.elites.values() if e["objs"].get("CV", 1) == 0)))
    dump(out, a.out)
    print(f"[merge] {a.variant}: {rep['status']}; ячеек {rep['merged_filled']} из {GRID.n_cells} ({100 * rep['merged_coverage']:.1f} %), "
          f"допустимых {out['n_feasible']}, оценок {rep['evals']}", flush=True)


def _merged(path):
    m = load(path)
    return m, A.Archive.from_dict(m["archive"])


def effect_ctx(m, variant):
    mu, sd = np.asarray(m["feats_mu"]), np.asarray(m["feats_sd"])
    ref = SM.reference(variant)["feats"]
    ref_z, _, _ = FT.standardize(ref, mu, sd)
    ar = A.Archive.from_dict(m["archive"])
    pool = np.array([e["meta"]["feats"] for e in ar.elites.values() if e["objs"].get("CV", 1) == 0] or
                    [e["meta"]["feats"] for e in ar.elites.values()])
    pool_z, _, _ = FT.standardize(pool, mu, sd)
    return dict(mu=mu, sd=sd, ref_z=ref_z, pool_z=pool_z)


def novelty(feats, ctx):
    Fz, _, _ = FT.standardize(np.atleast_2d(feats), ctx["mu"], ctx["sd"])
    return FT.min_distance(Fz, ctx["ref_z"]), FT.knn_novelty(Fz, ctx["pool_z"], k=15, exclude_self=True)


def cmd_nsga(a):
    m, ar = _merged(a.merged)
    ctx = effect_ctx(m, a.variant)
    ref = SM.reference(a.variant)

    def ev(X):
        F, CV = [], []
        for x in X:
            r = SM.evaluate_params(SM.decode(x), a.variant)
            o = r["objs"]
            nov, _ = novelty(r["feats"], ctx)
            F.append([o["mass_kg"] / ref["mass"], o["n_eff"] / ref["n_eff"], 1.0 - o["service"] / 10.0, 1.0 - o["q"], -float(nov[0])])
            CV.append(r["CV"])
        return np.array(F), np.array(CV)
    el = sorted(ar.elites.values(), key=lambda e: (-e["fitness"], e["uid"]))
    init = np.stack([e["genome"] for e in el[:a.pop]])
    rng = IS.rng_for(a.seed, 0, salt="nsga-" + a.variant)
    t0 = time.time()
    res = N2.nsga2(ev, SM.DIM, a.pop, a.gens, rng, init=init)
    front = [i for i in range(len(res["X"])) if res["rank"][i] == 0 and res["CV"][i] <= 0]
    rows = []
    for i in front:
        r = SM.evaluate_params(SM.decode(res["X"][i]), a.variant)
        nov, knn = novelty(r["feats"], ctx)
        rows.append(dict(genome=res["X"][i].tolist(), F=res["F"][i].tolist(), objs=r["objs"], desc=r["desc"].tolist(),
                         cell=list(GRID.cell(r["desc"])), novelty_oem=float(nov[0]), novelty_knn=float(knn[0]),
                         fitness=SM.fitness(r["objs"], a.variant)))
    Ff = np.array([r["F"] for r in rows]) if rows else np.zeros((0, 5))
    hv = None
    if len(Ff):
        ideal = Ff.min(axis=0) - 1e-9
        refp = np.array([1.2, 1.2, 1.0, 1.0, 0.0])
        hv = N2.hv_mc(Ff, refp, np.minimum(ideal, refp - 1e-6), 200000, IS.rng_for(1, 0, "hv"))
    out = dict(variant=a.variant, pop=a.pop, gens=a.gens, evals=res["evals"], seconds=round(time.time() - t0, 1), front=rows,
               objectives=["масса / типовая", "N_eff / типовая", "1 − ремонт/10", "1 − изготовимость", "− новизна против типовых форм OEM"],
               hv_mc=hv, hv_ref=[1.2, 1.2, 1.0, 1.0, 0.0], history=res["history"][::max(1, len(res["history"]) // 30)])
    dump(out, a.out)
    print(f"[nsga] {a.variant}: фронт {len(rows)} допустимых, оценок {res['evals']}, гиперобъём (МК) {hv}", flush=True)


def cmd_polish(a):
    m, ar = _merged(a.merged)
    el = sorted([e for e in ar.elites.values() if e["objs"].get("CV", 1) == 0], key=lambda e: (-e["fitness"], e["uid"]))
    picks = load(a.picks)["uids"] if a.picks else []
    chosen = el[:a.top] + [e for e in ar.elites.values() if e["uid"] in picks and e not in el[:a.top]]
    rows = []
    for k, e in enumerate(chosen):
        cell0 = tuple(e["cell"])

        def f(x):
            r = SM.evaluate(x, a.variant)
            if r["desc"] is None or not np.isfinite(r["fitness"]):
                return 1e3
            pen = 0.0 if GRID.cell(r["desc"]) == cell0 else 5.0
            return -r["fitness"] + pen
        rng = IS.rng_for(a.seed, k, salt="cma-" + a.variant)
        t0 = time.time()
        res = CM.fmin(f, e["genome"], 0.015, rng, a.evals)
        r = SM.evaluate(res["x"], a.variant)
        rows.append(dict(uid=e["uid"], cell=list(cell0), fitness_before=e["fitness"], fitness_after=r["fitness"],
                         same_cell=GRID.cell(r["desc"]) == cell0, genome=res["x"].tolist(), objs_before=e["objs"], objs_after=r["objs"],
                         evals=res["evals"], seconds=round(time.time() - t0, 1)))
        print(f"[polish] {a.variant} {e['uid']}: {e['fitness']:.4f} → {r['fitness']:.4f} (та же ячейка: {rows[-1]['same_cell']})", flush=True)
    dump(dict(variant=a.variant, rows=rows), a.out)


def select_gallery(ar, n_max=36):
    """галерея: лучшая допустимая элита в каждом сочетании (топология рёбер, симметрия); затем — лучшие по доле решётки и кривизне"""
    feas = [e for e in ar.elites.values() if e["objs"].get("CV", 1) == 0]
    best = {}
    for e in feas:
        key = (e["cell"][1], e["cell"][2])
        if key not in best or (e["fitness"], e["uid"]) > (best[key]["fitness"], best[key]["uid"]):
            best[key] = e
    sel = sorted(best.values(), key=lambda e: (e["cell"][1], e["cell"][2]))
    seen = {e["uid"] for e in sel}
    for j in (0, 3):
        for b in range(GRID.bins[j]):
            cand = [e for e in feas if e["cell"][j] == b and e["uid"] not in seen]
            if cand:
                e = max(cand, key=lambda e: (e["fitness"], e["uid"]))
                sel.append(e); seen.add(e["uid"])
    return sel[:n_max]


def cmd_gallery(a):
    import shell_geom as SG
    from core import render as R
    m, ar = _merged(a.merged)
    ctx = effect_ctx(m, a.variant)
    sel = select_gallery(ar, a.n)
    os.makedirs(a.out, exist_ok=True)
    pngs, labels, rows = [], [], []
    tag = "o" if a.variant == "оболочка" else "k"
    px = a.px_per_m
    for k, e in enumerate(sel):
        p = SM.decode(e["genome"])
        meshes = SG.meshes(p, a.variant, res=a.res)
        nov, knn = novelty(e["meta"]["feats"], ctx)
        lab = (f"{tag}{k + 1:02d} {SM.DESC_LABELS[1][e['cell'][1]]}, {SM.DESC_LABELS[2][e['cell'][2]]}; "
               f"{e['objs']['mass_kg']:.0f} кг, N_eff {e['objs']['n_eff']:.0f}")
        fn = os.path.join(a.out, f"thumb_{tag}{k + 1:02d}.png")
        R.render_axo(meshes, fn, size=(420, 340), px_per_m=px, backend=a.backend)
        pngs.append(fn); labels.append(lab)
        rows.append(dict(n=f"{tag}{k + 1:02d}", uid=e["uid"], cell=e["cell"], fitness=e["fitness"], objs=e["objs"], novelty_oem=float(nov[0]),
                         novelty_knn=float(knn[0]), thumb=os.path.basename(fn), params={kk: p[kk] for kk in
                         ("route", "rib_type", "sym_n", "rib_h", "rib_w", "rib_fillet", "f_n", "f_lean", "f_bow", "r_n", "r_lean", "r_bow",
                          "hf_n", "hf_shape", "lat_type", "cell", "emb", "cast_mat")}))
        print(f"[gallery] {lab}", flush=True)
    cs = R.contact_sheet(pngs, labels, os.path.join(a.out, f"gallery_{tag}.png"), ncols=6,
                         title=f"Корпус-оболочка KV-3, вариант «{a.variant}»: элиты MAP-Elites (лучшая в сочетании «рёбра × симметрия»), "
                               f"аксонометрия в одном масштабе")
    sheets = []
    hl = sorted(rows, key=lambda r: -r["fitness"])[:2] + sorted(rows, key=lambda r: -r["novelty_oem"])[:2]
    done = set()
    for r in hl:
        if r["uid"] in done:
            continue
        done.add(r["uid"])
        e = next(x for x in sel if x["uid"] == r["uid"])
        p = SM.decode(e["genome"])
        fn = os.path.join(a.out, f"sheet_{r['n']}.png")
        o = e["objs"]
        R.render_sheet(SG.meshes(p, a.variant, res=a.res), fn,
                       f"Корпус-оболочка KV-3, «{a.variant}», элита {r['n']} ({e['uid']})",
                       [f"ячейка: {', '.join(SM.DESC_LABELS[j][c] for j, c in enumerate(e['cell']))}; способ: {p['route']}",
                        f"масса {o['mass_kg']:.0f} кг; деталей {o['parts']}, крепёж {o['fasteners']}, N_eff {o['n_eff']:.0f}; ремонт {o['service']:.2f}/10; "
                        f"изготовимость {o['q']:.2f}",
                        f"наименьшая частота {o['f_min_Hz']:.0f} Гц (нужно ≥ {SM.F_MARGIN * SM.INP['f_op']:.0f}); загрузка стенки {o['u_wall']:.2f}, "
                        f"клетки {o['u_cage']:.2f}; новизна против типовых форм {r['novelty_oem']:.2f}",
                        "все числа — дешёвые модели уровня облика (оценка); проверка — CalculiX (shell_ccx.py)"],
                       px_per_m=1.4 * px, backend=a.backend, view_px=(760, 600))
        sheets.append(os.path.basename(fn))
    dump(dict(variant=a.variant, rows=rows, contact=os.path.basename(cs), sheets=sheets, px_per_m=px), os.path.join(a.out, f"gallery_{tag}.json"))


def cmd_picks(a):
    m, ar = _merged(a.merged)
    uids = [u.strip() for u in a.uids.split(",") if u.strip()]
    by = {e["uid"]: e for e in ar.elites.values()}
    miss = [u for u in uids if u not in by]
    if miss:
        raise SystemExit(f"нет таких элит: {miss}")
    dump(dict(variant=m["variant"], uids=uids, genomes=[by[u]["genome"].tolist() for u in uids], note="выбор владельца по галерее"), a.out)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    sp = ap.add_subparsers(dest="cmd", required=True)
    p = sp.add_parser("run-part")
    p.add_argument("--variant", required=True, type=lambda v: ALIAS[v], choices=SM.VARIANTS)
    p.add_argument("--n-islands", type=int, required=True); p.add_argument("--n-parts", type=int, required=True)
    p.add_argument("--part", type=int, required=True); p.add_argument("--evals", type=int, required=True)
    p.add_argument("--batch", type=int, default=64); p.add_argument("--seed", type=int, default=20261115)
    p.add_argument("--workers", type=int, default=1); p.add_argument("--picks"); p.add_argument("--out", default=OUT)
    p = sp.add_parser("merge")
    p.add_argument("--variant", required=True, type=lambda v: ALIAS[v], choices=SM.VARIANTS); p.add_argument("--n-islands", type=int, required=True)
    p.add_argument("--in", dest="inp", nargs="+", required=True); p.add_argument("--out", required=True)
    p.add_argument("--allow-partial", action="store_true")
    p = sp.add_parser("nsga")
    p.add_argument("--variant", required=True, type=lambda v: ALIAS[v], choices=SM.VARIANTS); p.add_argument("--merged", required=True)
    p.add_argument("--pop", type=int, default=120); p.add_argument("--gens", type=int, default=60)
    p.add_argument("--seed", type=int, default=20261115); p.add_argument("--out", required=True)
    p = sp.add_parser("polish")
    p.add_argument("--variant", required=True, type=lambda v: ALIAS[v], choices=SM.VARIANTS); p.add_argument("--merged", required=True)
    p.add_argument("--top", type=int, default=3); p.add_argument("--evals", type=int, default=1500)
    p.add_argument("--picks"); p.add_argument("--seed", type=int, default=20261115); p.add_argument("--out", required=True)
    p = sp.add_parser("gallery")
    p.add_argument("--variant", required=True, type=lambda v: ALIAS[v], choices=SM.VARIANTS); p.add_argument("--merged", required=True)
    p.add_argument("--n", type=int, default=36); p.add_argument("--res", default="draft"); p.add_argument("--backend", default="auto")
    p.add_argument("--px-per-m", type=float, default=230.0); p.add_argument("--out", required=True)
    p = sp.add_parser("picks")
    p.add_argument("--merged", required=True); p.add_argument("--uids", required=True); p.add_argument("--out", required=True)
    a = ap.parse_args()
    dict(**{"run-part": cmd_run_part, "merge": cmd_merge, "nsga": cmd_nsga, "polish": cmd_polish, "gallery": cmd_gallery,
            "picks": cmd_picks})[a.cmd](a)


if __name__ == "__main__":
    main()
