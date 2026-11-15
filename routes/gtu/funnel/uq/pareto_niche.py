"""Многокритериальный поиск внутри каждой ниши (NSGA-II) - фронт Парето ниши (на AWS). Всё - оценка.

Ниша - как у отбора уровня 0 (funnel_l0.niche_of: архетип × роторы × компрессор × турбина × охлаждение × рекуператор × изготовление
ротора [× корпус у v3]). Внутри ниши свободны остальные дискретные оси шаблона архетипа (ядра, вход, камера, подшипники, сплавы и
способы изготовления горячей части, камеры и корпуса, соединения, картридж, разъём лопаток) и все непрерывные оси оптимизации
уровня 1 (funnel_l1.var_spec, включая число ступеней турбины). Ген - число в [0, 1]: дискретная ось - слой значения (тот же выбор,
что у генератора уровня 0 funnel_space._discrete с закреплёнными осями ниши), непрерывная - доля окна var_spec. Прирост давления
RDC / волнового ротора - не ось (свойство технологии): середина диапазона, RDC 0,05, ВР 0,075 (оценка).

Цели (все - минимум, нормированы 0-1): 1 − детали/10, 1 − Plug&Play/10, 1 − монтаж/10 (баллы моделей уровня 0 для геометрии уровня 1),
(45 − КПД на клеммах, %)/30, ln(объём габарита, м³)/ln 60. Ограничение - правила расчётной точки уровня 1 (funnel_l1.evaluate_point:
сумма относительных превышений = 0). Режим --robust K: ещё P(расчётная точка допустима) по K общим выборкам неопределённостей
(uq_factors, первые K строк гиперкуба) >= P_TARGET - ограничение «устойчиво допустим».
Начальная популяция: варианты уровня 1 этой ниши (если даны части funnel_l1), концепции отбора уровня 0 этой ниши, остальное -
латинский гиперкуб. Итог ниши: фронт (недоминируемые допустимые последней популяции), до K_CHECK представителей (крайние по
каждой цели, лучшие по взвешенной и итоговой оценке, остальные - наиболее удалённые) проходят нерасчётные режимы (full_check) с
узлами ВНА; сравнение с вариантами уровня 1 ниши (доминирует ли фронт их точки); гиперобъём на половине и в конце (сходимость).
Деление на части - единственное: ниша k (порядок - по убыванию лучшей оценки уровня 0) идёт в часть k % S.

  python pareto_niche.py --select ../l0_aws_full/funnel_l0_select_KV-3_n200000.json --pop 200 --gens 100 --shards 50 --shard 0,6 \
      --workers "$NPROC" --gzip --out "$GTU_OUT" [--l1-in ../aws/l1_full --l1-m 3983 --l1-shards 40] [--robust 16 --tag _rob16]
Слияние и файл кандидатов для устойчивого отбора - merge_pareto.py.
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
import nsga2 as NS               # noqa: E402
import uq_factors as UF          # noqa: E402
import uq_patch as UP            # noqa: E402
import uq_l1 as U                # noqa: E402

LEVEL = "PARETO"
PART_PREFIX = "pareto"
DATE = U.DATE
FREE_DISC = ["cores", "n_ax", "inlet", "comb", "bear", "mat_cw", "mat_tw", "m_hot", "mat_hot", "m_comb", "m_case", "hot_integral",
             "recup_integral", "conn", "cartridge", "split_bl"]
CONT_GENES = ["t41", "pi", "split", "mu1", "mu2", "phi1", "lt", "kcm", "zc", "zt", "zn", "phit", "kUpt", "phipt", "eps_rec", "kD_case",
              "nt", "npt"]
GENES = FREE_DISC + CONT_GENES
NG = len(GENES)
OBJ = ["1 − s_parts/10", "1 − s_plugplay/10", "1 − s_service/10", "(45 − η_el)/30", "ln V_env / ln 60"]
PGAIN_MID = {"RDC": 0.05, "ВР": 0.075}
P_TARGET = 0.5                # режим --robust: P(расчётная точка) по K выборкам не ниже (оценка: как P_MIN устойчивого отбора)
K_CHECK = 24                  # представителей фронта на нерасчётные режимы (около 30 с каждый)
L0_SEEDS = 16                 # концепций отбора уровня 0 ниши в начальной популяции (concept() - около 0,5 с)
DEFAULT_POP, DEFAULT_GENS = 200, 100
LOCAL_MAX = dict(niches=2, pop=24, gens=6, check=2)
PARETO_FILES = ("nsga2.py", "pareto_niche.py")
_UGRID = sorted({(k + 0.5) / L for L in range(1, 9) for k in range(L)})


def code_sha():
    h = hashlib.sha256()
    for f in PARETO_FILES:
        h.update(open(os.path.join(HERE, f), "rb").read())
    h.update(U.code_sha().encode())
    return h.hexdigest()[:16]


def niche_axes():
    return list(getattr(L0, "NICHE_AXES", ["arch", "rotors", "comp", "turb", "cool", "recup", "rotor_make"]))


def niche_values(niche):
    """строка ниши -> значения закреплённых осей (ov генератора) и прочие оси ниши (корпус у v3)"""
    vals = dict(zip(niche_axes(), niche.split("|")))
    ov = dict(rotors=int(vals["rotors"].rstrip("р")), comp=vals["comp"], turb=vals["turb"], cool=vals["cool"], recup=vals["recup"])
    rm = vals["rotor_make"]
    join, make = rm.split(":", 1)
    ov["rotor_join"] = join
    if "/" in make:
        ov["m_cw"], ov["m_tw"] = make.split("/", 1)
    else:
        ov["m_tw"] = ov["m_cw"] = make
    extra = {k: v for k, v in vals.items() if k not in ("arch", "rotors", "comp", "turb", "cool", "recup", "rotor_make")}
    return vals["arch"], ov, extra


class Niche:
    """одна ниша: опорная концепция (лучшая концепция отбора уровня 0 ниши - производные поля), построение концепции из генов"""

    def __init__(self, niche, kv, ref_idx, n0, seed):
        self.niche, self.kv, self.ref_idx, self.n0, self.seed = niche, kv, int(ref_idx), int(n0), int(seed)
        self.arch, self.ov, self.extra = niche_values(niche)
        self.tgt = FS.target(kv)
        self.ref = FS.concept(self.ref_idx, self.n0, kv, self.seed)
        self.cache = {}

    def build(self, g):
        u = {a: 0.5 for a in FS.CAT}
        u.update({a: float(x) for a, x in zip(FREE_DISC, g[:len(FREE_DISC)])})
        d = FS._discrete(u, self.arch, ov=self.ov)
        c = dict(self.ref)
        c.update(d)
        c.update(self.extra)
        c["kv"] = self.kv
        c["pgain"] = PGAIN_MID.get(c["comb"], 0.0)
        c["repair"] = "фронт Парето ниши (pareto_niche)"
        spec = L1.var_spec(c)
        gc = dict(zip(CONT_GENES, g[len(FREE_DISC):]))
        x = np.array([gc.get(s[0], 0.5) for s in spec])
        c = L1.x_to_c(c, x, spec)
        if L0.niche_of(c) != self.niche:
            return None
        return c

    def encode(self, c):
        """концепция ниши -> гены (дискретные - слой, который даёт то же значение; непрерывные - доля окна var_spec)"""
        g = np.full(NG, 0.5)
        u = {a: 0.5 for a in FS.CAT}
        for k, a in enumerate(FREE_DISC):
            for x in _UGRID:
                u[a] = x
                if FS._discrete(u, self.arch, ov=self.ov).get(a) == c.get(a):
                    break
            g[k] = u[a]
        spec = L1.var_spec(c)
        x = L1.c_to_x(c, spec)
        for (name, *_), xi in zip(spec, x):
            g[len(FREE_DISC) + CONT_GENES.index(name)] = xi
        return g

    def evaluate(self, g):
        key = tuple(np.round(g, 9))
        if key in self.cache:
            return self.cache[key]
        c = self.build(g)
        if c is None:
            r = dict(ok=False, viol=10.0)
        else:
            ev = L1.evaluate_point(c, self.tgt)
            r = objectives(ev)
            r["c"] = c
        self.cache[key] = r
        return r


def objectives(ev):
    if not ev["ok"]:
        return dict(ok=False, viol=float(ev["viol"]))
    s = ev["scores"]
    eta = float(ev["eta_el"])
    V = float(ev["layout"]["envelope"]["V"])
    G = [1 - s["parts"] / 10, 1 - s["plugplay"] / 10, 1 - s["service"] / 10, (45.0 - eta) / 30.0, math.log(max(V, 1e-3)) / math.log(60.0)]
    sw = float(L0.weighted({k: s[k] for k in L0.CRIT}))
    trl = int(ev["models"]["mr"]["trl_min"])
    return dict(ok=True, viol=float(ev["viol"]), G=[float(x) for x in G], score_w=sw, score_f=sw * U.trl_factor(trl), trl=trl, eta_el=eta,
                V_env=V)


def evaluator(N, robust_vals=None, stats=None):
    """X -> (F, V) для NSGA-II; robust_vals - список значений факторов (K выборок): V += max(0, P_TARGET − P̂)"""
    def f(X):
        R = [N.evaluate(x) for x in X]
        F = np.array([r["G"] if r["ok"] else [2.0] * 5 for r in R])
        V = np.array([r["viol"] if r["ok"] else max(5.0, r["viol"]) for r in R], dtype=float)
        if stats is not None:
            stats["evals"] += len(X)
            stats["feasible"] += int((V <= 0).sum())
        if robust_vals:
            idx = [i for i, r in enumerate(R) if r["ok"] and r["viol"] <= 0 and "p_rob" not in r]
            if idx:
                hits = np.zeros(len(idx))
                for v in robust_vals:
                    with UP.Perturb(v):
                        for m, i in enumerate(idx):
                            c = dict(R[i]["c"])
                            if c.get("comb") in PGAIN_MID:
                                c["pgain"] = c["pgain"] * float(v["k_pgain"])
                            UP.set_variant(c)
                            ev = L1.evaluate_point(c, N.tgt)
                            hits[m] += 1.0 if (ev["ok"] and ev["viol"] == 0) else 0.0
                for m, i in enumerate(idx):
                    R[i]["p_rob"] = float(hits[m] / len(robust_vals))
                if stats is not None:
                    stats["robust_evals"] += len(idx) * len(robust_vals)
            for i, r in enumerate(R):
                if r["ok"] and r["viol"] <= 0:
                    V[i] = max(0.0, P_TARGET - r.get("p_rob", 0.0))
        return F, V
    return f


# =========================================================================== представители фронта и проверка
def representatives(G, scores_w, scores_f, k):
    n = len(G)
    if n <= k:
        return list(range(n))
    pick = []
    for m in range(G.shape[1]):
        i = int(np.argmin(G[:, m]))
        if i not in pick:
            pick.append(i)
    for s in (scores_w, scores_f):
        i = int(np.argmax(s))
        if i not in pick:
            pick.append(i)
    d = np.min(np.linalg.norm(G[:, None, :] - G[None, pick, :], axis=2), axis=1)
    while len(pick) < k:
        i = int(np.argmax(d))
        pick.append(i)
        d = np.minimum(d, np.linalg.norm(G - G[i], axis=1))
    return pick[:k]


def full(N, c):
    """нерасчётные режимы и узлы ВНА для члена фронта -> итоговые цели и допустимость (как отбор уровня 1)"""
    ev = L1.evaluate_point(c, N.tgt)
    if not ev["ok"] or ev["viol"] > 0:
        return dict(feasible=0, reject="|".join(r["code"] for r in ev["bad"]))
    mx, stt, R = L1.full_check(c, N.tgt, ev)
    bad = [r["code"] for r in R if not r["ok"]]
    if mx["vigv_used"] or mx["heat_used"]:
        ev = L1.add_vigv(c, N.tgt, ev, vigv=mx["vigv_used"], heat=mx["heat_used"])
    o = objectives(ev)
    perf, lay, m = ev["perf"], ev["layout"], ev["models"]
    rec = dict(kv=N.kv, niche=N.niche, arch=c["arch"], feasible=int(not bad), score=round(o["score_w"], 5), eta_th_shaft=round(perf["eta_th"] * 100, 3),
               eta_el=round(o["eta_el"], 3), V_env=round(o["V_env"], 4), mass_kg=round(float(lay["mass_kg_estimate"]), 0),
               t41=round(float(c["t41"]), 3), pi=round(float(c["pi"]), 4), eps_rec=round(float(c.get("eps_rec") or 0.0), 4), trl_min=o["trl"],
               vigv=int(bool(mx["vigv_used"])), inlet_heat=int(bool(mx["heat_used"])), reject="|".join(bad), n_parts=m["pr"]["parts"],
               n_assemblies=m["pr"]["assemblies"], steps=m["sr"]["steps"], s_parts=ev["scores"]["parts"], s_plugplay=ev["scores"]["plugplay"],
               s_service=ev["scores"]["service"], SM_min=mx["SM_min"])
    rec.update({a: c.get(a) for a in M1.DISC_W})
    return dict(feasible=int(not bad), G=o["G"], score_w=o["score_w"], score_f=o["score_f"], rec=rec, od_bad=bad)


def l1_G(rec):
    try:
        return [1 - rec["s_parts"] / 10, 1 - rec["s_plugplay"] / 10, 1 - rec["s_service"] / 10, (45.0 - rec["eta_el"]) / 30.0,
                math.log(max(rec["V_env"], 1e-3)) / math.log(60.0)]
    except (KeyError, TypeError):
        return None


def _r(x, k=5):
    return U._r(x, k)


def run_niche(task):
    """одна ниша: NSGA-II, фронт, проверка представителей, сравнение с уровнем 1"""
    t0, tc0 = time.time(), time.process_time()
    k, niche, kv, ref_idx, n0, seed, l0_idx, l1_seeds, P = task
    N = Niche(niche, kv, ref_idx, n0, seed)
    rseed = int(hashlib.sha256(f"{UF.UQ_SEED}|{kv}|{niche}|{P['tag']}".encode()).hexdigest()[:8], 16)
    init = []
    for cc in l1_seeds:
        init.append(N.encode(cc["concept"]))
    for idx in l0_idx[:L0_SEEDS]:
        try:
            init.append(N.encode(FS.concept(int(idx), n0, kv, seed)))
        except (ValueError, ZeroDivisionError, OverflowError, KeyError):
            pass
    rv = None
    if P["robust"]:
        _, _, vals = UF.sample(max(P["robust"], 2), UF.UQ_SEED)
        rv = vals[:P["robust"]]
    stats = dict(evals=0, feasible=0, robust_evals=0)
    ev_fn = evaluator(N, rv, stats)
    half = {}

    def cb(g, X, F, V, rank):
        if g == max(1, P["gens"] // 2):
            m = NS.front_mask(F, V)
            half.update(gen=g, hv=NS.hv_mc(F[m]) if m.any() else 0.0, front=int(m.sum()))
    res = NS.run(ev_fn, NG, 5, pop=P["pop"], gens=P["gens"], seed=rseed, init=np.array(init) if init else None, callback=cb)
    X, F, V = res["X"], res["F"], res["V"]
    m = NS.front_mask(F, V)
    fi = np.where(m)[0]
    hv = NS.hv_mc(F[m]) if m.any() else 0.0
    front = []
    for i in fi:
        r = N.evaluate(X[i])
        front.append(dict(g=[round(float(x), 5) for x in X[i]], G=[round(x, 5) for x in r["G"]], score_w=_r(r["score_w"]), score_f=_r(r["score_f"]),
                          trl=r["trl"], eta_el=_r(r["eta_el"], 3), V_env=_r(r["V_env"], 4), p_rob=_r(r.get("p_rob"), 4)))
    checked = []
    if front:
        G = np.array([f["G"] for f in front])
        reps = representatives(G, np.array([f["score_w"] for f in front]), np.array([f["score_f"] for f in front]), P["check"])
        for q in reps:
            c = N.build(np.array(front[q]["g"]))
            fr = full(N, c) if c is not None else dict(feasible=0, reject="ниша не совпала")
            checked.append(dict(front_i=int(q), g=front[q]["g"], G_nominal=front[q]["G"], **{x: (_r(v) if isinstance(v, float) else v)
                                                                                                for x, v in fr.items()}))
    cmp_ = None
    l1f = [cc for cc in l1_seeds if (cc.get("rec") or {}).get("feasible")]
    if l1f:
        ok = [ch for ch in checked if ch.get("feasible")]
        Gc = np.array([ch["G"] for ch in ok]) if ok else np.zeros((0, 5))
        dom = 0
        for cc in l1f:
            g1 = l1_G(cc["rec"])
            if g1 is not None and len(Gc) and np.any(np.all(Gc <= np.array(g1), axis=1) & np.any(Gc < np.array(g1), axis=1)):
                dom += 1
        best_l1 = max(cc["rec"]["score"] for cc in l1f)
        best_fr = max([ch["score_w"] for ch in ok], default=None)
        cmp_ = dict(l1_feasible=len(l1f), dominated_by_checked_front=dom, best_l1_score=_r(best_l1), best_front_score_checked=_r(best_fr),
                    gain=_r(best_fr - best_l1) if best_fr is not None else None)
    nfe = sum(1 for ch in checked if ch.get("feasible"))
    return dict(k=k, niche=niche, kv=kv, arch=N.arch, ref_idx=ref_idx, n_l1_seeds=len(l1_seeds), n_l0_seeds=min(len(l0_idx), L0_SEEDS),
                evals=stats["evals"], feasible_evals=stats["feasible"], robust_evals=stats["robust_evals"], cache=len(N.cache),
                history=res["history"][::max(1, len(res["history"]) // 10)], hv=round(hv, 5), half=half,
                front_size=len(front), front=front, checked=checked, checked_feasible=nfe,
                best_score_w=_r(max([f["score_w"] for f in front], default=None)), best_score_f=_r(max([f["score_f"] for f in front], default=None)),
                l1_compare=cmp_, t_s=round(time.time() - t0, 2), t_cpu_s=round(time.process_time() - tc0, 2))


# =========================================================================== ниши, части, пул
def load_niches(select_path):
    d = json.load(open(select_path, encoding="utf-8"))
    by = {}
    for e in d["selected"]:
        by.setdefault(e["niche"], []).append(e)
    rows = []
    for niche, es in by.items():
        es.sort(key=lambda e: (not e.get("feasible"), -(e.get("score") or 0.0), e["idx"]))
        rows.append(dict(niche=niche, ref_idx=es[0]["idx"], l0_idx=[e["idx"] for e in es], best=es[0].get("score") or 0.0,
                         feasible=bool(es[0].get("feasible"))))
    rows.sort(key=lambda r: (not r["feasible"], -r["best"], r["niche"]))
    meta = dict(file=os.path.basename(select_path), sha256=hashlib.sha256(open(select_path, "rb").read()).hexdigest(), kv=d["kv"], n=int(d["n"]),
                seed=int(d["seed"]), l0_code_sha=d.get("code_sha"), fingerprint=d.get("fingerprint"), niches=len(rows))
    return rows, meta


def part_path(out, kv, NN, shard, shards, tag="", gz=False):
    return os.path.join(out, f"{PART_PREFIX}{tag}_{kv}_n{NN}_s{shard:03d}of{shards:03d}.json" + (".gz" if gz else ""))


def run_parts(select_path, shard_list, shards, out, pop=DEFAULT_POP, gens=DEFAULT_GENS, check=K_CHECK, robust=0, workers=1, gz=False,
              tag="", l1=None, top=None):
    t0 = time.time()
    niches, meta_in = load_niches(select_path)
    if top:
        niches = niches[:top]
    NN = len(niches)
    kv, n0, seed = meta_in["kv"], meta_in["n"], meta_in["seed"]
    L1._prep(n0, seed, kv)
    l1_by = {}
    l1_meta = None
    if l1:
        cands, l1_meta = U.load_l1(l1["dirs"], kv, l1["M"], l1["shards"], l1.get("tag", ""))
        for cd in cands:
            l1_by.setdefault(cd["niche"], []).append(cd)
    P = dict(pop=pop, gens=gens, check=check, robust=robust, tag=tag)
    idx_by = {k: FS.shard_indices(NN, k, shards) for k in shard_list}
    ks = sorted(k for v in idx_by.values() for k in v)
    tasks = []
    for k in ks:
        r = niches[k]
        seeds = sorted(l1_by.get(r["niche"], []), key=lambda cd: -((cd.get("rec") or {}).get("score") or 0.0))[:max(0, pop // 4)]
        seeds = [dict(concept=cd["concept"], rec=dict(cd.get("rec") or {})) for cd in seeds]
        tasks.append((k, r["niche"], kv, r["ref_idx"], n0, seed, r["l0_idx"], seeds, P))
    if workers <= 1 or len(tasks) < 2:
        got = [run_niche(t) for t in tasks]
    else:
        with U._pool(workers) as ex:
            got = list(ex.map(run_niche, tasks, chunksize=1))
    by = {g["k"]: g for g in got}
    el = time.time() - t0
    meta = dict(level=LEVEL, date=DATE, software=FS.software(), method=__doc__.split("\n")[0], kv=kv, NN=NN, shards=shards, tag=tag,
                input=meta_in, l1=l1_meta, code_sha=code_sha(), l1_code_sha=L1.code_sha(), factors_sha=UF.table_sha(), pop=pop, gens=gens,
                check=check, robust=robust, p_target=P_TARGET if robust else None, genes=GENES, objectives=OBJ, workers=workers,
                call=dict(shards=shard_list, niches=len(ks), elapsed_s=round(el, 2), evals=sum(g["evals"] for g in got),
                          robust_evals=sum(g["robust_evals"] for g in got), checked=sum(len(g["checked"]) for g in got),
                          cpu_s=round(sum(g["t_cpu_s"] for g in got), 2), host_cpus=os.cpu_count(),
                          started=datetime.datetime.fromtimestamp(t0, datetime.timezone.utc).isoformat()))
    os.makedirs(out, exist_ok=True)
    paths = []
    for k in shard_list:
        doc = dict(meta, shard=k, n_rows=len(idx_by[k]), rows=[by[i] for i in idx_by[k]])
        p = part_path(out, kv, NN, k, shards, tag, gz)
        txt = json.dumps(doc, ensure_ascii=False, separators=(",", ":"), default=_r)
        with (gzip.open(p, "wt", encoding="utf-8") if gz else open(p, "w", encoding="utf-8")) as fh:
            fh.write(txt)
        paths.append(p)
    return paths, meta, got


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--select", required=True, help="отбор уровня 0 (merge_l0: funnel_l0_select*.json) - список ниш и опорных концепций")
    ap.add_argument("--pop", type=int, default=DEFAULT_POP)
    ap.add_argument("--gens", type=int, default=DEFAULT_GENS)
    ap.add_argument("--check", type=int, default=K_CHECK, help="представителей фронта на нерасчётные режимы")
    ap.add_argument("--robust", type=int, default=0, help="K выборок неопределённостей в ограничении (0 - номинальный фронт)")
    ap.add_argument("--l1-in", nargs="*", default=None)
    ap.add_argument("--l1-m", type=int, default=None)
    ap.add_argument("--l1-shards", type=int, default=None)
    ap.add_argument("--l1-tag", default="")
    ap.add_argument("--top", type=int, default=0, help="первых ниш (0 - все)")
    ap.add_argument("--shards", type=int, required=True)
    ap.add_argument("--shard", default="all")
    ap.add_argument("--workers", type=int, default=0)
    ap.add_argument("--out", default=os.environ.get("GTU_OUT", "."))
    ap.add_argument("--gzip", action="store_true")
    ap.add_argument("--tag", default="")
    ap.add_argument("--local", action="store_true", help=f"локальная миниатюра: не больше {LOCAL_MAX}")
    a = ap.parse_args(argv)
    if a.local and ((a.top or 10 ** 9) > LOCAL_MAX["niches"] or a.pop > LOCAL_MAX["pop"] or a.gens > LOCAL_MAX["gens"] or a.check > LOCAL_MAX["check"]):
        ap.error(f"локально - не больше {LOCAL_MAX}")
    l1 = dict(dirs=a.l1_in, M=a.l1_m, shards=a.l1_shards, tag=a.l1_tag) if a.l1_in else None
    sl = L0.parse_shards(a.shard, a.shards)
    w = a.workers or L0.default_workers()
    paths, meta, got = run_parts(a.select, sl, a.shards, a.out, a.pop, a.gens, a.check, a.robust, w, a.gzip, a.tag, l1, a.top or None)
    c = meta["call"]
    print(f"{meta['kv']} фронты ниш: {meta['NN']} ниш, части {sl[0]}..{sl[-1]} ({len(sl)} из {a.shards}), {c['niches']} ниш, {c['evals']} "
          f"вычислений (+{c['robust_evals']} под неопределённостью), проверено {c['checked']}; {c['elapsed_s']} с, {w} процессов -> {a.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
