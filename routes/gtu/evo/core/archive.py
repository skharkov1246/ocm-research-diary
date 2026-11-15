"""MAP-Elites: сетка дескрипторов, архив элит и операторы изменчивости.

Метод — Mouret & Clune 2015 (arXiv:1504.04909, «Illuminating search spaces by mapping elites»); оператор iso+line —
Vassiliades & Mouret 2018 (GECCO, «Discovering the elite hypervolume by leveraging interspecies correlation»);
SBX и полиномиальная мутация — Deb & Agrawal 1995, Deb 2001 (см. sources.py, статус ссылок там).

Геном — вектор в [0, 1]^d. Фитнес — чем больше, тем лучше; nan/inf в архив не попадают.
Детерминизм: выбор родителей — по отсортированному списку ячеек; при равном фитнесе остаётся элита с меньшим uid.
"""
import math

import numpy as np


class GridSpec:
    """Сетка дескрипторов. Для каждого дескриптора — границы бинов edges (длина bins + 1); значения за краями
    попадают в крайние бины (клип)."""

    def __init__(self, names, edges, labels=None):
        if len(names) != len(edges):
            raise ValueError("names и edges разной длины")
        self.names = list(names)
        self.edges = [np.asarray(e, float) for e in edges]
        for e in self.edges:
            if e.ndim != 1 or len(e) < 2 or np.any(np.diff(e) <= 0):
                raise ValueError("edges: возрастающие границы, не меньше двух")
        self.bins = [len(e) - 1 for e in self.edges]
        self.labels = labels or [[f"{e[i]:g}…{e[i + 1]:g}" for i in range(len(e) - 1)] for e in self.edges]

    @classmethod
    def uniform(cls, names, bins, lows, highs):
        return cls(names, [np.linspace(lo, hi, b + 1) for b, lo, hi in zip(bins, lows, highs)])

    @property
    def n_cells(self):
        return int(np.prod(self.bins))

    def cell(self, desc):
        d = np.asarray(desc, float)
        if d.shape != (len(self.edges),) or not np.all(np.isfinite(d)):
            return None
        out = []
        for v, e, b in zip(d, self.edges, self.bins):
            i = int(np.searchsorted(e[1:-1], v, side="right"))
            out.append(min(max(i, 0), b - 1))
        return tuple(out)

    def to_dict(self):
        return dict(names=self.names, edges=[e.tolist() for e in self.edges], labels=self.labels)

    @classmethod
    def from_dict(cls, d):
        return cls(d["names"], d["edges"], d.get("labels"))


class Archive:
    """Архив MAP-Elites: ячейка → элита {genome, fitness, desc, objs, meta, uid, cell}."""

    def __init__(self, grid):
        self.grid = grid
        self.elites = {}
        self.n_tried = 0
        self.n_added = 0

    def __len__(self):
        return len(self.elites)

    def try_add(self, genome, fitness, desc, objs=None, meta=None, uid=""):
        self.n_tried += 1
        if fitness is None or not np.isfinite(fitness):
            return False
        c = self.grid.cell(desc)
        if c is None:
            return False
        cur = self.elites.get(c)
        if cur is not None:
            if fitness < cur["fitness"] or (fitness == cur["fitness"] and str(uid) >= str(cur["uid"])):
                return False
        self.elites[c] = dict(genome=np.asarray(genome, float).copy(), fitness=float(fitness), desc=np.asarray(desc, float).copy(),
                              objs=objs or {}, meta=meta or {}, uid=str(uid), cell=c)
        self.n_added += 1
        return True

    def cells(self):
        return sorted(self.elites)

    def coverage(self):
        return len(self.elites) / self.grid.n_cells

    def qd_score(self, offset=0.0):
        """QD-score: сумма (фитнес − offset) по заполненным ячейкам (offset — нижняя граница фитнеса)."""
        return float(sum(e["fitness"] - offset for e in self.elites.values()))

    def best(self):
        if not self.elites:
            return None
        return min(self.elites.values(), key=lambda e: (-e["fitness"], e["uid"]))

    def sample_genomes(self, rng, k):
        cs = self.cells()
        idx = rng.integers(0, len(cs), size=k)
        return np.stack([self.elites[cs[i]]["genome"] for i in idx])

    def genomes(self):
        cs = self.cells()
        return np.stack([self.elites[c]["genome"] for c in cs]) if cs else np.zeros((0, 0))

    def to_dict(self):
        return dict(grid=self.grid.to_dict(), n_tried=self.n_tried, n_added=self.n_added,
                    elites=[dict(cell=list(c), genome=self.elites[c]["genome"].tolist(), fitness=self.elites[c]["fitness"],
                                 desc=self.elites[c]["desc"].tolist(), objs=self.elites[c]["objs"], meta=self.elites[c]["meta"],
                                 uid=self.elites[c]["uid"]) for c in self.cells()])

    @classmethod
    def from_dict(cls, d):
        a = cls(GridSpec.from_dict(d["grid"]))
        for e in d["elites"]:
            c = tuple(e["cell"])
            a.elites[c] = dict(genome=np.asarray(e["genome"], float), fitness=float(e["fitness"]), desc=np.asarray(e["desc"], float),
                               objs=e.get("objs", {}), meta=e.get("meta", {}), uid=str(e["uid"]), cell=c)
        a.n_tried, a.n_added = d.get("n_tried", 0), d.get("n_added", 0)
        return a


# --------------------------------------------------------------------------- операторы изменчивости (геном в [0, 1]^d)
def gaussian(rng, x, sigma):
    return np.clip(x + sigma * rng.standard_normal(x.shape), 0.0, 1.0)


def iso_line(rng, x1, x2, s_iso=0.01, s_line=0.2):
    """iso+line (Vassiliades & Mouret 2018): x1 + σ_iso·N(0, I) + σ_line·N(0, 1)·(x2 − x1)."""
    return np.clip(x1 + s_iso * rng.standard_normal(x1.shape) + s_line * rng.standard_normal() * (x2 - x1), 0.0, 1.0)


def sbx(rng, p1, p2, eta=15.0, prob=0.9):
    """SBX с границами [0, 1] (Deb & Agrawal 1995; вариант с ограниченной областью — Deb 2001). Два потомка."""
    c1, c2 = p1.copy(), p2.copy()
    if rng.random() > prob:
        return c1, c2
    for i in range(len(p1)):
        if rng.random() > 0.5 or abs(p1[i] - p2[i]) < 1e-14:
            continue
        y1, y2 = min(p1[i], p2[i]), max(p1[i], p2[i])
        u = rng.random()
        out = []
        for beta in (1.0 + 2.0 * y1 / (y2 - y1), 1.0 + 2.0 * (1.0 - y2) / (y2 - y1)):
            alpha = 2.0 - beta ** -(eta + 1.0)
            if u <= 1.0 / alpha:
                bq = (u * alpha) ** (1.0 / (eta + 1.0))
            else:
                bq = (1.0 / (2.0 - u * alpha)) ** (1.0 / (eta + 1.0))
            out.append(bq)
        a = 0.5 * ((y1 + y2) - out[0] * (y2 - y1))
        b = 0.5 * ((y1 + y2) + out[1] * (y2 - y1))
        a, b = min(max(a, 0.0), 1.0), min(max(b, 0.0), 1.0)
        if rng.random() <= 0.5:
            a, b = b, a
        c1[i], c2[i] = a, b
    return c1, c2


def poly_mutation(rng, x, eta=20.0, pm=None):
    """полиномиальная мутация с границами [0, 1] (Deb 2001)"""
    y = x.copy()
    d = len(x)
    pm = 1.0 / d if pm is None else pm
    for i in range(d):
        if rng.random() > pm:
            continue
        u = rng.random()
        d1, d2 = y[i], 1.0 - y[i]
        if u < 0.5:
            xy = 1.0 - d1
            val = 2.0 * u + (1.0 - 2.0 * u) * xy ** (eta + 1.0)
            dq = val ** (1.0 / (eta + 1.0)) - 1.0
        else:
            xy = 1.0 - d2
            val = 2.0 * (1.0 - u) + 2.0 * (u - 0.5) * xy ** (eta + 1.0)
            dq = 1.0 - val ** (1.0 / (eta + 1.0))
        y[i] = min(max(y[i] + dq, 0.0), 1.0)
    return y


DEFAULT_OPS = (("iso_line", 0.5), ("gauss", 0.3), ("sbx", 0.2))


def make_children(rng, archive, n, ops=DEFAULT_OPS, sigma=0.05, s_iso=0.01, s_line=0.2, pm=None):
    """n потомков из архива; оператор — по вероятностям ops"""
    names = [o for o, _ in ops]
    p = np.array([w for _, w in ops], float); p /= p.sum()
    out = []
    while len(out) < n:
        op = names[int(rng.choice(len(names), p=p))]
        if op == "gauss":
            x = archive.sample_genomes(rng, 1)[0]
            out.append(poly_mutation(rng, gaussian(rng, x, sigma), pm=pm) if pm else gaussian(rng, x, sigma))
        elif op == "iso_line":
            x1, x2 = archive.sample_genomes(rng, 2)
            out.append(iso_line(rng, x1, x2, s_iso, s_line))
        elif op == "sbx":
            x1, x2 = archive.sample_genomes(rng, 2)
            c1, c2 = sbx(rng, x1, x2)
            out.append(poly_mutation(rng, c1, pm=pm))
            if len(out) < n:
                out.append(poly_mutation(rng, c2, pm=pm))
        else:
            raise ValueError(op)
    return np.stack(out[:n])


def map_elites(evaluate, dim, archive, rng, n_evals, batch=64, n_init=None, init=None, ops=DEFAULT_OPS, sigma=0.05,
               s_iso=0.01, s_line=0.2, uid_prefix="", log_every=0, log=None, eval_offset=0):
    """Основной цикл MAP-Elites.

    evaluate(G (B, dim)) → список словарей {fitness, desc, objs?, meta?} той же длины.
    init — начальные геномы (родители: выбор владельца, архив прошлого этапа); затем случайные до n_init.
    Возврат — история [{evals, filled, coverage, qd, best}] по пакетам.
    """
    n_init = max(batch, 2) if n_init is None else n_init
    hist, done = [], 0
    pending = [np.clip(np.asarray(g, float), 0.0, 1.0) for g in (init if init is not None else [])]
    offset = None

    def run(G):
        nonlocal done, offset
        res = evaluate(G)
        if len(res) != len(G):
            raise RuntimeError("evaluate вернул не столько результатов, сколько геномов")
        for g, r in zip(G, res):
            uid = f"{uid_prefix}{eval_offset + done:09d}"
            archive.try_add(g, r.get("fitness"), r.get("desc"), r.get("objs"), r.get("meta"), uid)
            done += 1
        fits = [e["fitness"] for e in archive.elites.values()]
        if fits:
            offset = min(fits) if offset is None else min(offset, min(fits))
        b = archive.best()
        rec = dict(evals=done, filled=len(archive), coverage=archive.coverage(), best=(b["fitness"] if b else None))
        hist.append(rec)
        if log and log_every and (len(hist) % log_every == 0):
            log(rec)

    while pending and done < n_evals:
        k = min(batch, len(pending), n_evals - done)
        run(np.stack(pending[:k]))
        pending = pending[k:]
    while done < min(n_init, n_evals) or (len(archive) == 0 and done < n_evals):
        k = min(batch, n_evals - done)
        run(rng.random((k, dim)))
    while done < n_evals:
        k = min(batch, n_evals - done)
        run(make_children(rng, archive, k, ops, sigma, s_iso, s_line))
    for rec in hist:
        rec["qd_offset"] = offset
    return hist


def qd_offset_default(archives):
    fits = [e["fitness"] for a in archives for e in a.elites.values()]
    return min(fits) if fits else 0.0


def coverage_table(archive):
    """заполненность по каждому дескриптору (доля бинов, где есть хоть одна элита)"""
    out = {}
    for j, (name, b) in enumerate(zip(archive.grid.names, archive.grid.bins)):
        used = sorted({c[j] for c in archive.elites})
        out[name] = dict(bins=b, used=len(used), which=used)
    return out


def _isfinite(x):
    return x is not None and math.isfinite(x)
