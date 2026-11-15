#!/usr/bin/env python3
"""Поиск формы семейства «лопатки и колёса»: острова MAP-Elites → слияние с контролем полноты → NSGA-II (Парето) →
CMA-ES (доводка элит в своей ячейке) → элиты для CalculiX и галереи. Интерактивная эволюция: выбор владельца из
галереи (picks.json) становится родителями следующего поколения (--parents).

Команды:
  python3 evo_blades.py islands --comp comp2 --n-islands 64 --part 0 --n-parts 8 --evals 20000 --out DIR
  python3 evo_blades.py merge   --comp comp2 --n-islands 64 --in DIR [DIR ...] --out merged_comp2.json
  python3 evo_blades.py polish  --comp comp2 --merged merged_comp2.json --nsga-pop 64 --nsga-gens 30 --cma-elites 8 --out DIR
  python3 evo_blades.py elites  --comp comp2 --merged merged_comp2.json --k 16 --out elites_comp2.json [--mesh-dir DIR]
  python3 evo_blades.py mini    --out DIR          # мини-прогон на этой машине (сотни оценок на компонент)
  python3 evo_blades.py part    --part P --n-parts N ...   # всё для одной части на машине AWS (острова → свои элиты → ccx)

Острова независимы: зерно — seed_of(run_key, «blades», comp, island), миграции нет, поэтому итог слияния не зависит
от разбиения на части и числа процессов (тест test_blades.py::test_partition_invariance). Архив — CVT (Vassiliades et
al. 2018), пригодность одна на все острова (иначе слияние «лучшего в ячейке» сравнивало бы разные величины).
"""
from __future__ import annotations

import argparse
import glob
import json
import math
import os
import sys
import time
from multiprocessing import Pool

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in sys.path:
    sys.path.insert(0, HERE)
import core_local as C   # noqa: E402

SCRATCH = os.environ.get("EVO_OUT", os.path.join(os.environ.get("TMPDIR", "/tmp"), "evo_blades"))
CODE_FILES = [os.path.join(HERE, f) for f in ("core_local.py", "blade_geom.py", "blade_eval.py", "disk1d.py", "evo_blades.py")]

# =========================================================================== конфигурация (входит в хэш)
CONFIG = dict(
    family="blades", run_key="kv3-blades-v1",
    archive=dict(kind="cvt", k={"mini": 128, "aws": 1024}),
    me=dict(batch=32, s_iso=(0.01, 0.06), s_line=(0.1, 0.4), p_gauss=0.3, s_gauss=0.08, init_frac=0.1),
    fitness=dict(comp="100·η + 50·min(SM, 0,30) − 20·(J/J_ref − 1) − 10·(m/m_ref − 1); недопустимые: −1000 − нарушение",
                 turb="100·η − 20·(J/J_ref − 1) − 10·(m/m_ref − 1); недопустимые: −1000 − нарушение"),
    point="routes/gtu/frame3/frame3_point_results.json",
)


def config_hash(scale):
    c = dict(CONFIG)
    c["scale"] = scale
    return C.stable_hash(c)


def code_hash():
    return C.files_hash(CODE_FILES)


def make_archive(comp, scale):
    import blade_eval as BE
    D = BE.DESC[comp]
    k = CONFIG["archive"]["k"][scale]
    import blade_geom as BG
    key = f"{CONFIG['run_key']}:cvt:{comp}:{k}"
    return C.CVTArchive(k, D["lo"], D["hi"], len(BG.genes_for(comp)), seed_key=key, centroids=_centroids(key, k, len(D["lo"])))


def _centroids(key, k, dim):
    """центроиды CVT детерминированы ключом; кэш на диске (k = 1024 — около 15 с на ядро)"""
    os.makedirs(SCRATCH, exist_ok=True)
    p = os.path.join(SCRATCH, f"cvt_{C.seed_of(key, k, dim):x}.npy")
    if os.path.exists(p):
        Cc = np.load(p)
        if Cc.shape == (k, dim):
            return Cc
    Cc = C.cvt_centroids(k, dim, seed_key=key)
    tmp = p + f".{os.getpid()}.tmp.npy"
    np.save(tmp, Cc)
    os.replace(tmp, p)
    return Cc


def slim_meta(m):
    """что хранить в архиве на ячейку (полный расчёт элиты повторяется детерминированно по геному)"""
    keep = ("comp", "eta", "eta_level1", "SM", "mass_kg", "J_kgm2", "u_disk", "burst", "route", "corrections", "D2", "beta2b", "Zeff",
            "parts", "novelty_vec", "sig_bend_MPa", "hollow", "Dm", "Dtip", "why", "viol")
    out = {k: m.get(k) for k in keep if k in m}
    if "stages" in m:
        out["stages"] = [dict(stage=s["stage"], life_h=s["life_h"], sig_c_MPa=s["sig_c_MPa"], sig_b_MPa=s["sig_b_MPa"],
                              disk_life_h=s["disk_life_h"], disk_burst=s["disk_burst"], R_hub=s["R_hub"]) for s in m["stages"]]
    if "features" in m:
        out["features"] = {k: v for k, v in m["features"].items() if isinstance(v, (int, float))}
    if "surge" in m and isinstance(m["surge"], dict):
        out["surge_crit"] = m["surge"].get("crit")
    return C._canon(out)


def make_evaluator(comp):
    import blade_eval as BE

    def ev(G):
        out = []
        for g in G:
            r = BE.evaluate(comp, g)
            out.append(dict(f=r["f"], d=r["d"], feasible=True, meta=slim_meta(r["meta"]), ok=r["feasible"]))
        return out
    return ev


# =========================================================================== остров


def run_island(args):
    comp, island, n_evals, scale, out_dir, parents = args
    t0 = time.time()
    rng = C.rng_of(CONFIG["run_key"], CONFIG["family"], comp, island)
    me = CONFIG["me"]
    # шаги мутации острова — детерминированно из номера (разнообразие поиска без разных пригодностей)
    q = (island * 0.618033988749895) % 1.0
    s_iso = me["s_iso"][0] + q * (me["s_iso"][1] - me["s_iso"][0])
    s_line = me["s_line"][0] + ((island * 0.414213562) % 1.0) * (me["s_line"][1] - me["s_line"][0])
    arch = make_archive(comp, scale)
    ev = make_evaluator(comp)
    import blade_eval as BE
    init = [BE.baseline_genes(comp)] if island == 0 else []
    par = np.asarray(parents, float) if parents is not None and len(parents) else None
    log = C.map_elites(ev, arch, arch.n_genes, rng, n_evals, batch=me["batch"], n_init=max(me["batch"], int(me["init_frac"] * n_evals)),
                       s_iso=s_iso, s_line=s_line, p_gauss=me["p_gauss"], s_gauss=me["s_gauss"], init=np.array(init) if init else None,
                       parents=par, log_every=max(1, n_evals // me["batch"] // 20))
    feas = int(np.sum([m is not None and m.get("why") is None for m in arch.meta]))
    res = dict(family=CONFIG["family"], comp=comp, island=island, seed=C.seed_of(CONFIG["run_key"], CONFIG["family"], comp, island),
               config_hash=config_hash(scale), code_hash=code_hash(), scale=scale, n_evals=n_evals, s_iso=s_iso, s_line=s_line,
               archive=arch.to_dict(), log=log, feasible_cells=feas, seconds=round(time.time() - t0, 2),
               parents=len(par) if par is not None else 0)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)
        with open(os.path.join(out_dir, f"island_{comp}_{island:04d}.json"), "w", encoding="utf-8") as f:
            json.dump(res, f, ensure_ascii=False)
    return res


def run_islands(comp, islands, n_evals, scale, out_dir, workers=None, parents=None):
    jobs = [(comp, i, n_evals, scale, out_dir, parents) for i in islands]
    workers = workers or os.cpu_count() or 1
    if workers <= 1 or len(jobs) == 1:
        return [run_island(j) for j in jobs]
    with Pool(min(workers, len(jobs))) as pool:
        return list(pool.imap(run_island, jobs))      # порядок — по номеру острова


def load_islands(dirs, comp):
    out = []
    for d in dirs:
        for p in sorted(glob.glob(os.path.join(d, "**", f"island_{comp}_*.json"), recursive=True)):
            out.append(json.load(open(p, encoding="utf-8")))
    return out


def merge(comp, dirs, n_islands, scale, allow_partial=False, check_code=True):
    isl = load_islands(dirs, comp)
    a, rep = C.merge_island_archives(isl, range(n_islands), config_hash(scale), code_hash() if check_code else None,
                                     allow_partial=allow_partial, centroids=make_archive(comp, scale).C)
    feas = [c for c in a.elites() if a.meta[c] and a.meta[c].get("why") is None] if a is not None else []
    rep["feasible_cells"] = len(feas)
    rep["comp"] = comp
    return a, rep


# =========================================================================== NSGA-II и CMA-ES


def objectives(comp, r, ref):
    m = r["meta"]
    if comp == "turb":
        F = [-m.get("eta", 0.0), m.get("J_kgm2", 1e3) / ref["J"], m.get("mass_kg", 1e4) / ref["mass"]]
    else:
        F = [-m.get("eta", 0.0), -min(m.get("SM", 0.0) or 0.0, 0.30), m.get("J_kgm2", 1e3) / ref["J"], m.get("mass_kg", 1e4) / ref["mass"],
             m.get("u_disk", 9.0)]
    cv = 0.0 if r["feasible"] else float(m.get("viol", 1.0))
    return F, cv


def nsga_stage(comp, init, pop, gens, seed_tag="nsga"):
    import blade_eval as BE
    ref = BE.reference(comp)

    def ev(G):
        Fs, CVs, Ms = [], [], []
        for g in G:
            r = BE.evaluate(comp, g)
            F, cv = objectives(comp, r, ref)
            Fs.append(F)
            CVs.append(cv)
            Ms.append(dict(slim_meta(r["meta"]), f=r["f"], d=np.asarray(r["d"]).tolist(), feasible=r["feasible"]))
        return np.array(Fs), np.array(CVs), Ms
    n = len(init[0]) if len(init) else None
    res = C.nsga2(ev, n, C.rng_of(CONFIG["run_key"], seed_tag, comp), pop=pop, gens=gens, init=init)
    front = [dict(g=res["P"][i].tolist(), F=res["F"][i].tolist(), meta=res["meta"][i]) for i in range(len(res["P"]))
             if res["rank"][i] == 0 and res["CV"][i] <= 0]
    return dict(front=front, log=res["log"], n_evals=pop * (gens + 1))


def cma_polish(comp, arch, cell, max_evals, seed_tag="cma"):
    """доводка элиты в её ячейке: CMA-ES по −f; выход из ячейки — штраф (разнообразие архива сохраняется)"""
    import blade_eval as BE
    x0 = arch.G[cell].copy()

    def fun(X):
        out = []
        for x in X:
            r = BE.evaluate(comp, x)
            pen = 0.0 if arch.index(np.asarray(r["d"], float)) == cell else 50.0
            out.append(-(r["f"] if r["feasible"] else -1000.0) + pen)
        return np.array(out)
    xb, fb, n, log = C.cmaes_minimize(fun, x0, 0.05, C.rng_of(CONFIG["run_key"], seed_tag, comp, cell), max_evals=max_evals, lam=8)
    r = BE.evaluate(comp, xb)
    return dict(cell=int(cell), f0=float(arch.F[cell]), f1=float(r["f"]), feasible=r["feasible"], g=xb.tolist(), n_evals=n,
                meta=slim_meta(r["meta"]), d=np.asarray(r["d"]).tolist())


# =========================================================================== элиты


def select_elites(comp, arch, k=16, top_frac=0.5):
    """элиты для CalculiX и галереи: среди допустимых ячеек верхняя доля по пригодности, из неё k самых непохожих
    (выбор дальней точки в нормированном векторе признаков формы, первая — лучшая по пригодности)"""
    cells = [c for c in arch.elites() if arch.meta[c] and arch.meta[c].get("why") is None]
    if not cells:
        return []
    cells = sorted(cells, key=lambda c: -arch.F[c])
    top = cells[:max(k, int(math.ceil(top_frac * len(cells))))]
    V = np.array([arch.meta[c]["novelty_vec"] for c in top], float)
    chosen = [0]
    dmin = np.linalg.norm(V - V[0], axis=1)
    while len(chosen) < min(k, len(top)):
        j = int(np.argmax(dmin))
        chosen.append(j)
        dmin = np.minimum(dmin, np.linalg.norm(V - V[j], axis=1))
    return [top[j] for j in chosen]


def effect_metrics(comp, arch, cells):
    """эффектность — измеримые признаки (решает владелец по галерее): новизна в архиве (k = 5 ближайших), расстояние до
    типовых форм OEM, гладкость (энергия плавности закона углов, C2 на узлах), сложность силуэта"""
    import blade_eval as BE
    allc = [c for c in arch.elites() if arch.meta[c] and arch.meta[c].get("novelty_vec")]
    Vall = np.array([arch.meta[c]["novelty_vec"] for c in allc], float)
    oem = BE.oem_reference_vectors(comp)
    out = {}
    for c in cells:
        v = np.array(arch.meta[c]["novelty_vec"], float)[None, :]
        d = np.sqrt(((Vall - v) ** 2).sum(1))
        d = np.sort(d[d > 1e-12])[:5]
        f = arch.meta[c].get("features", {})
        out[int(c)] = dict(novelty_archive=float(d.mean()) if len(d) else 0.0,
                           novelty_oem=float(np.min(np.sqrt(((oem - v) ** 2).sum(1)))),
                           fairness_law=f.get("fair_law"), c2_jump_law=f.get("c2_jump_law"),
                           silhouette_iso=f.get("silhouette_iso"), silhouette_inflections=f.get("silhouette_inflections"),
                           le_rmin_m=f.get("le_rmin_m"))
    return out


def export_elite_meshes(comp, arch, cells, mesh_dir):
    """сетки элит для рендера (NPZ: V, T, серый на каждую часть) — рендер идёт отдельным окружением с pyvista"""
    import blade_eval as BE
    import blade_geom as BG
    os.makedirs(mesh_dir, exist_ok=True)
    files = []
    for c in cells:
        r = BE.evaluate(comp, arch.G[c], want_geom=True)
        if comp == "turb":
            parts = []
            zt = int(BG.decode("turb", arch.G[c])["zt"])
            x0 = 0.0
            for k, (bl, dk) in enumerate(zip(r["geom"]["blades"], r["geom"]["disks"])):
                Vb, Tb = BG.blade_mesh_from_sections(bl, zt)
                prof = np.vstack([np.c_[-0.5 * dk["h_eff"], dk["r"]], np.c_[(0.5 * dk["h_eff"])[::-1], dk["r"][::-1]]])
                Vd, Td = BG.revolve(prof, 96)
                pitch = 0.16 + 0.03 * k
                Vb[:, 0] += x0
                Vd[:, 0] += x0
                parts += [(Vd, Td, 0.80), (Vb, Tb, 0.60)]
                x0 += pitch
            meshes = parts
        else:
            meshes = r["geom"]["mesh"]
        p = os.path.join(mesh_dir, f"{comp}_cell{int(c):04d}.npz")
        arrs = {}
        for i, (V, T, gray) in enumerate(meshes):
            arrs[f"V{i}"] = V.astype(np.float32)
            arrs[f"T{i}"] = T.astype(np.int32)
            arrs[f"g{i}"] = np.array([gray])
        np.savez_compressed(p, n=np.array([len(meshes)]), **arrs)
        files.append(p)
    return files


def elite_records(comp, arch, cells):
    import blade_eval as BE
    import blade_geom as BG
    em = effect_metrics(comp, arch, cells)
    out = []
    for rank, c in enumerate(cells):
        r = BE.evaluate(comp, arch.G[c])
        m = r["meta"]
        rec = dict(id=f"{comp}-{int(c):04d}", cell=int(c), rank=rank, f=r["f"], feasible=r["feasible"],
                   genes=BG.decode(comp, arch.G[c]), g=arch.G[c].tolist(), d=np.asarray(r["d"]).tolist(),
                   desc_names=BE.DESC[comp]["names"], effect=em[int(c)], meta=C._canon({k: v for k, v in m.items() if k not in ("novelty_vec",)}))
        out.append(rec)
    return out


# =========================================================================== CLI


def _parents(path, comp):
    if not path:
        return None
    d = json.load(open(path, encoding="utf-8"))
    ids = d.get(comp, [])
    G = []
    for e in d.get("elites", []):
        if e.get("id") in ids or e.get("cell") in ids:
            G.append(e["g"])
    return np.array(G) if G else None


def cmd_islands(a):
    isl = C.islands_of_part(a.part, a.n_parts, a.n_islands)
    res = run_islands(a.comp, isl, a.evals, a.scale, a.out, a.workers, _parents(a.parents, a.comp))
    print(json.dumps(dict(comp=a.comp, islands=isl, filled=[r["archive"]["cells"].__len__() for r in res],
                          feasible=[r["feasible_cells"] for r in res], seconds=[r["seconds"] for r in res]), ensure_ascii=False))


def cmd_merge(a):
    arch, rep = merge(a.comp, a.inp, a.n_islands, a.scale, a.allow_partial)
    out = dict(report=rep, archive=arch.to_dict())
    with open(a.out, "w", encoding="utf-8") as f:
        json.dump(out, f, ensure_ascii=False)
    print(json.dumps(rep, ensure_ascii=False))


def _load_merged(path, comp, scale):
    d = json.load(open(path, encoding="utf-8"))
    return C.archive_from_dict(d["archive"], make_archive(comp, scale).C), d["report"]


def cmd_polish(a):
    arch, rep = _load_merged(a.merged, a.comp, a.scale)
    cells = select_elites(a.comp, arch, k=max(a.nsga_pop, a.cma_elites))
    out = dict(comp=a.comp)
    if a.nsga_gens > 0 and cells:
        init = np.array([arch.G[c] for c in cells][:a.nsga_pop])
        ns = nsga_stage(a.comp, init, a.nsga_pop, a.nsga_gens)
        out["nsga"] = ns
        for e in ns["front"]:          # фронт возвращается в архив (в свою ячейку, если лучше)
            if e["meta"]["feasible"]:
                arch.try_add(np.asarray(e["g"]), e["meta"]["f"], np.asarray(e["meta"]["d"]), {k: v for k, v in e["meta"].items() if k not in ("f", "d", "feasible")})
    cm = []
    for c in select_elites(a.comp, arch, k=a.cma_elites):
        pr = cma_polish(a.comp, arch, c, a.cma_evals)
        cm.append(pr)
        if pr["feasible"] and pr["f1"] > arch.F[c]:
            arch.try_add(np.asarray(pr["g"]), pr["f1"], np.asarray(pr["d"]), pr["meta"])
    out["cma"] = cm
    out["archive"] = arch.to_dict()
    out["report"] = rep
    os.makedirs(a.out, exist_ok=True)
    with open(os.path.join(a.out, f"polished_{a.comp}.json"), "w", encoding="utf-8") as f:
        json.dump(C._canon(out), f, ensure_ascii=False)
    print(json.dumps(dict(nsga_front=len(out.get("nsga", {}).get("front", [])), cma=[(p["cell"], round(p["f0"], 3), round(p["f1"], 3)) for p in cm]),
                     ensure_ascii=False))


def cmd_elites(a):
    d = json.load(open(a.merged, encoding="utf-8"))
    arch = C.archive_from_dict(d["archive"], make_archive(a.comp, a.scale).C)
    cells = select_elites(a.comp, arch, a.k)
    recs = elite_records(a.comp, arch, cells)
    if a.mesh_dir:
        export_elite_meshes(a.comp, arch, cells, a.mesh_dir)
    with open(a.out, "w", encoding="utf-8") as f:
        json.dump(C._canon(dict(comp=a.comp, elites=recs)), f, ensure_ascii=False, indent=1)
    print(f"{len(recs)} элит → {a.out}")


def cmd_part(a):
    """одна часть на машине AWS: острова всех компонент → слияние своих островов (неполный набор — это норма части) →
    доводка → свои элиты → проверка элит в CalculiX (если ccx есть)"""
    t0 = time.time()
    summary = dict(part=a.part, n_parts=a.n_parts, comps={})
    for comp in a.comps.split(","):
        isl = C.islands_of_part(a.part, a.n_parts, a.n_islands)
        d_isl = os.path.join(a.out, "islands")
        run_islands(comp, isl, a.evals, a.scale, d_isl, a.workers)
        arch, rep = merge(comp, [d_isl], a.n_islands, a.scale, allow_partial=True)
        cells = select_elites(comp, arch, a.k)
        recs = elite_records(comp, arch, cells)
        p = os.path.join(a.out, f"elites_{comp}_part{a.part:03d}.json")
        with open(p, "w", encoding="utf-8") as f:
            json.dump(C._canon(dict(comp=comp, part=a.part, elites=recs, report=rep)), f, ensure_ascii=False)
        cc = None
        if a.ccx:
            import ccx_elites as CE
            cc = CE.check_file(p, os.path.join(a.out, "ccx"), workers=a.workers)
        summary["comps"][comp] = dict(report=rep, elites=len(recs), ccx=None if cc is None else cc.get("summary"))
    summary["seconds"] = round(time.time() - t0, 1)
    with open(os.path.join(a.out, f"part_{a.part:03d}_summary.json"), "w", encoding="utf-8") as f:
        json.dump(C._canon(summary), f, ensure_ascii=False, indent=1)
    print(json.dumps(C._canon(summary), ensure_ascii=False)[:2000])


def cmd_mini(a):
    """мини-прогон на этой машине: на каждый компонент n_islands × evals оценок MAP-Elites, слияние с проверкой полноты,
    NSGA-II и CMA-ES в малом объёме, элиты и сетки для рендера"""
    t0 = time.time()
    out = dict(title="Семейство «лопатки и колёса»: мини-прогон (локально)", date=time.strftime("%Y-%m-%d"), comps={})
    os.makedirs(a.out, exist_ok=True)
    for comp in a.comps.split(","):
        tc = time.time()
        d_isl = os.path.join(a.out, "islands")
        run_islands(comp, range(a.n_islands), a.evals, "mini", d_isl, a.workers, _parents(a.parents, comp))
        arch, rep = merge(comp, [d_isl], a.n_islands, "mini")
        merged_p = os.path.join(a.out, f"merged_{comp}.json")
        with open(merged_p, "w", encoding="utf-8") as f:
            json.dump(dict(report=rep, archive=arch.to_dict()), f, ensure_ascii=False)
        ns = None
        cells0 = select_elites(comp, arch, k=a.nsga_pop)
        if a.nsga_gens > 0 and len(cells0) >= 4:
            init = np.array([arch.G[c] for c in cells0])
            ns = nsga_stage(comp, init, a.nsga_pop, a.nsga_gens)
            for e in ns["front"]:
                arch.try_add(np.asarray(e["g"]), e["meta"]["f"], np.asarray(e["meta"]["d"]),
                             {k: v for k, v in e["meta"].items() if k not in ("f", "d", "feasible")})
        cm = []
        for c in select_elites(comp, arch, k=a.cma_elites):
            pr = cma_polish(comp, arch, c, a.cma_evals)
            cm.append(pr)
            if pr["feasible"] and pr["f1"] > arch.F[c]:
                arch.try_add(np.asarray(pr["g"]), pr["f1"], np.asarray(pr["d"]), pr["meta"])
        cells = select_elites(comp, arch, a.k)
        recs = elite_records(comp, arch, cells)
        export_elite_meshes(comp, arch, cells, os.path.join(a.out, "meshes"))
        import blade_eval as BE
        ref = BE.reference(comp)
        feas = [c for c in arch.elites() if arch.meta[c] and arch.meta[c].get("why") is None]
        fv = np.array([arch.F[c] for c in feas]) if feas else np.array([])
        out["comps"][comp] = dict(
            merge_report=rep, n_evals_me=rep["n_evals_total"], n_evals_nsga=ns["n_evals"] if ns else 0,
            n_evals_cma=sum(p["n_evals"] for p in cm), cells_total=arch.n_cells, cells_filled=int(arch.filled.sum()),
            cells_feasible=len(feas), coverage_feasible=len(feas) / arch.n_cells, f_best=float(fv.max()) if len(fv) else None,
            f_median=float(np.median(fv)) if len(fv) else None,
            reference=dict(mass_kg=ref["mass"], J_kgm2=ref["J"], eta=ref["eta"], feasible=ref["raw"]["feasible"],
                           why=ref["raw"]["meta"].get("why"), f=ref["raw"]["f"],
                           SM=ref["raw"]["meta"].get("SM")),
            nsga_front=[dict(F=e["F"], eta=e["meta"].get("eta"), SM=e["meta"].get("SM"), mass_kg=e["meta"].get("mass_kg"),
                             J_kgm2=e["meta"].get("J_kgm2")) for e in (ns["front"] if ns else [])],
            nsga_log=ns["log"] if ns else None,
            cma=[dict(cell=p["cell"], f0=p["f0"], f1=p["f1"], feasible=p["feasible"], n_evals=p["n_evals"]) for p in cm],
            elites=recs, seconds=round(time.time() - tc, 1))
        merged_p2 = os.path.join(a.out, f"merged_{comp}_polished.json")
        with open(merged_p2, "w", encoding="utf-8") as f:
            json.dump(dict(report=rep, archive=arch.to_dict()), f, ensure_ascii=False)
    out["seconds"] = round(time.time() - t0, 1)
    out["config_hash"] = config_hash("mini")
    out["code_hash"] = code_hash()
    with open(os.path.join(a.out, "mini_results.json"), "w", encoding="utf-8") as f:
        json.dump(C._canon(out), f, ensure_ascii=False, indent=1)
    print(json.dumps({k: dict(cells_feasible=v["cells_feasible"], f_best=v["f_best"], evals=v["n_evals_me"] + v["n_evals_nsga"] + v["n_evals_cma"],
                              seconds=v["seconds"]) for k, v in out["comps"].items()}, ensure_ascii=False))


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    s = sub.add_parser("islands")
    s.add_argument("--comp", required=True, choices=["comp1", "comp2", "turb"])
    s.add_argument("--n-islands", type=int, required=True)
    s.add_argument("--part", type=int, default=0)
    s.add_argument("--n-parts", type=int, default=1)
    s.add_argument("--evals", type=int, default=2000)
    s.add_argument("--scale", default="aws", choices=["mini", "aws"])
    s.add_argument("--workers", type=int, default=None)
    s.add_argument("--parents", default=None)
    s.add_argument("--out", required=True)
    s.set_defaults(fn=cmd_islands)
    s = sub.add_parser("merge")
    s.add_argument("--comp", required=True, choices=["comp1", "comp2", "turb"])
    s.add_argument("--n-islands", type=int, required=True)
    s.add_argument("--in", dest="inp", nargs="+", required=True)
    s.add_argument("--scale", default="aws", choices=["mini", "aws"])
    s.add_argument("--allow-partial", action="store_true")
    s.add_argument("--out", required=True)
    s.set_defaults(fn=cmd_merge)
    s = sub.add_parser("polish")
    s.add_argument("--comp", required=True, choices=["comp1", "comp2", "turb"])
    s.add_argument("--merged", required=True)
    s.add_argument("--scale", default="aws", choices=["mini", "aws"])
    s.add_argument("--nsga-pop", type=int, default=64)
    s.add_argument("--nsga-gens", type=int, default=30)
    s.add_argument("--cma-elites", type=int, default=8)
    s.add_argument("--cma-evals", type=int, default=240)
    s.add_argument("--out", required=True)
    s.set_defaults(fn=cmd_polish)
    s = sub.add_parser("elites")
    s.add_argument("--comp", required=True, choices=["comp1", "comp2", "turb"])
    s.add_argument("--merged", required=True)
    s.add_argument("--scale", default="aws", choices=["mini", "aws"])
    s.add_argument("--k", type=int, default=16)
    s.add_argument("--mesh-dir", default=None)
    s.add_argument("--out", required=True)
    s.set_defaults(fn=cmd_elites)
    s = sub.add_parser("part")
    s.add_argument("--part", type=int, required=True)
    s.add_argument("--n-parts", type=int, required=True)
    s.add_argument("--n-islands", type=int, required=True)
    s.add_argument("--comps", default="comp1,comp2,turb")
    s.add_argument("--evals", type=int, default=4000)
    s.add_argument("--scale", default="aws", choices=["mini", "aws"])
    s.add_argument("--k", type=int, default=8)
    s.add_argument("--workers", type=int, default=None)
    s.add_argument("--ccx", action="store_true")
    s.add_argument("--out", required=True)
    s.set_defaults(fn=cmd_part)
    s = sub.add_parser("mini")
    s.add_argument("--comps", default="comp1,comp2,turb")
    s.add_argument("--n-islands", type=int, default=4)
    s.add_argument("--evals", type=int, default=300)
    s.add_argument("--nsga-pop", type=int, default=16)
    s.add_argument("--nsga-gens", type=int, default=4)
    s.add_argument("--cma-elites", type=int, default=2)
    s.add_argument("--cma-evals", type=int, default=48)
    s.add_argument("--k", type=int, default=8)
    s.add_argument("--workers", type=int, default=None)
    s.add_argument("--parents", default=None)
    s.add_argument("--out", default=os.path.join(SCRATCH, "mini"))
    s.set_defaults(fn=cmd_mini)
    a = ap.parse_args(argv)
    a.fn(a)


if __name__ == "__main__":
    main()
