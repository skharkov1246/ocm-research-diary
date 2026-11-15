"""Воронка камеры сгорания KV, уровни 1-3: извлечение итогов из кейса OpenFOAM без function objects.

Читает сетку (polyMesh: points, faces, owner, neighbour, boundary) и поля последнего (и предпоследнего) момента записи в
формате ascii; работает одинаково на apt-1912 (где function objects падают - навык openfoam-combustor, раздел 9) и на
v2512/v2606. Значение на участке границы - записанное (fixedValue, inletOutlet, calculated) или значение пристенной
ячейки (zeroGradient - по определению). Только numpy.

Что считается (имена полей results-JSON - в README, раздел «Уровень 1»):
  * выход: средние по расходу (вес - поток phi) T, состав, O2 сухой, NOx = NO + NO2 и CO в ppm сухого газа при 15 % O2
    (C15 = C_сух·(20,9 − 15)/(20,9 − O2,сух) - навык, раздел 7), несгоревший CH4, полнота сгорания по теплоте сгорания
    CH4, CO, H2 (значения - из входов кейса), радиальный профиль T и RTDF = (T_max профиля − T_ср)/(T_ср − T3);
  * полное давление p0 = p + ρ|U|²/2 на входах и выходе (по расходу), потери (p0_вх − p0_вых)/p0_вх;
  * баланс массы вход/выход; сходимость - изменение T, NOx, CO на выходе между двумя последними записями;
  * ячейка: неравномерность на срезе канала S = σ/среднее phi по атомам C, H, O (вес - расход через слой ячеек),
    наибольшая T в канале (проскок / самовоспламенение), доля обратного тока в зоне горения на 2 калибрах за срезом;
  * труба: температура газа в пристенных ячейках стенок (адиабатная стенка) по длине.
Запуск: python comb_post.py <кейс> [--kind cell|liner]
"""
import argparse, glob, json, math, os, re

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
FOAM_RES = os.path.join(HERE, "cfd", "foam")
R_U = 8.314462618
ATOMS = {"C": 12.011, "H": 1.008, "O": 15.999, "N": 14.007, "AR": 39.948}


def _strip_header(t):
    t = re.sub(r"/\*.*?\*/", "", t, flags=re.S)
    t = re.sub(r"//[^\n]*", "", t)
    m = re.search(r"FoamFile\s*\{.*?\}", t, flags=re.S)
    return t[m.end():] if m else t


def _list_block(t, start=0):
    """число и содержимое первого списка вида N ( ... ) начиная с позиции start"""
    m = re.compile(r"(\d+)\s*\(").search(t, start)
    if not m:
        raise ValueError("нет списка")
    n = int(m.group(1))
    depth, i = 1, m.end()
    while depth and i < len(t):
        ch = t[i]
        if ch == "(":
            depth += 1
        elif ch == ")":
            depth -= 1
        i += 1
    return n, t[m.end():i - 1], i


def read_points(path):
    n, body, _ = _list_block(_strip_header(open(path).read()))
    a = np.array(re.findall(r"\(\s*([-+0-9.eE]+)\s+([-+0-9.eE]+)\s+([-+0-9.eE]+)\s*\)", body), float)
    assert len(a) == n
    return a


def read_faces(path):
    n, body, _ = _list_block(_strip_header(open(path).read()))
    faces = [list(map(int, f.split())) for f in re.findall(r"\d+\s*\(([^)]*)\)", body)]
    assert len(faces) == n
    return faces


def read_labels(path):
    n, body, _ = _list_block(_strip_header(open(path).read()))
    a = np.array(body.split(), int)
    assert len(a) == n
    return a


def read_boundary(path):
    t = _strip_header(open(path).read())
    _, body, _ = _list_block(t)
    out = {}
    for m in re.finditer(r"(\w+)\s*\{([^}]*)\}", body):
        d = dict(re.findall(r"(\w+)\s+([^;]+);", m.group(2)))
        out[m.group(1)] = dict(type=d.get("type"), nFaces=int(d["nFaces"]), startFace=int(d["startFace"]))
    return out


def _parse_value(s, ncomp):
    s = s.strip()
    if s.startswith("uniform"):
        v = s[len("uniform"):].strip()
        if v.startswith("("):
            return np.array(v.strip("()").split(), float)
        return float(v)
    if s.startswith("nonuniform"):
        n, body, _ = _list_block(s)
        if ncomp == 1:
            a = np.array(body.split(), float)
        else:
            a = np.array(re.findall(r"\(([^)]*)\)", body), dtype=object)
            a = np.array([list(map(float, x.split())) for x in a], float)
        assert len(a) == n, (len(a), n)
        return a
    raise ValueError(s[:40])


def read_field(path):
    """поле: внутреннее (массив или число) и значения на участках (если записаны)"""
    t = _strip_header(open(path).read())
    cls = re.search(r"class\s+(\w+);", open(path).read()).group(1)
    ncomp = 3 if "Vector" in cls else 1
    m = re.search(r"internalField\s+(.*?);\s*\n\s*boundaryField", t, flags=re.S)
    internal = _parse_value(m.group(1), ncomp)
    bf = {}
    bstart = t.find("boundaryField")
    btxt = t[bstart:]
    for pm in re.finditer(r"\n\s*(\w+)\s*\n\s*\{", btxt):
        name = pm.group(1)
        j = pm.end()
        depth = 1
        k = j
        while depth and k < len(btxt):
            if btxt[k] == "{":
                depth += 1
            elif btxt[k] == "}":
                depth -= 1
            k += 1
        blk = btxt[j:k - 1]
        vm = re.search(r"\bvalue\s+((?:uniform|nonuniform).*?);\s*(?:\n|$)", blk, flags=re.S)
        tm = re.search(r"type\s+([\w:]+);", blk)
        bf[name] = dict(type=tm.group(1) if tm else None, value=_parse_value(vm.group(1), ncomp) if vm else None)
    return dict(cls=cls, internal=internal, boundary=bf, ncomp=ncomp)


class Mesh:
    def __init__(self, case):
        pm = os.path.join(case, "constant", "polyMesh")
        self.points = read_points(os.path.join(pm, "points"))
        self.faces = read_faces(os.path.join(pm, "faces"))
        self.owner = read_labels(os.path.join(pm, "owner"))
        self.neigh = read_labels(os.path.join(pm, "neighbour"))
        self.boundary = read_boundary(os.path.join(pm, "boundary"))
        self.nc = int(max(self.owner.max(), self.neigh.max() if len(self.neigh) else 0)) + 1
        nf = len(self.faces)
        Cf = np.zeros((nf, 3))
        Sf = np.zeros((nf, 3))
        for i, f in enumerate(self.faces):
            P = self.points[f]
            c0 = P.mean(axis=0)
            s = np.zeros(3)
            ca = np.zeros(3)
            a_tot = 0.0
            for k in range(len(f)):
                a, b = P[k], P[(k + 1) % len(f)]
                tri = 0.5 * np.cross(a - c0, b - c0)
                at = np.linalg.norm(tri)
                s += tri
                ca += at * (a + b + c0) / 3.0
                a_tot += at
            Cf[i] = ca / a_tot if a_tot > 0 else c0
            Sf[i] = s
        self.Cf, self.Sf = Cf, Sf
        cEst = np.zeros((self.nc, 3))
        cnt = np.zeros(self.nc)
        nin = len(self.neigh)
        np.add.at(cEst, self.owner, Cf)
        np.add.at(cnt, self.owner, 1)
        np.add.at(cEst, self.neigh, Cf[:nin])
        np.add.at(cnt, self.neigh, 1)
        cEst /= cnt[:, None]
        V = np.zeros(self.nc)
        C = np.zeros((self.nc, 3))
        pv = np.einsum("ij,ij->i", Sf, Cf - cEst[self.owner]) / 3.0
        np.add.at(V, self.owner, pv)
        np.add.at(C, self.owner, pv[:, None] * (0.75 * Cf + 0.25 * cEst[self.owner]))
        pvn = np.einsum("ij,ij->i", -Sf[:nin], Cf[:nin] - cEst[self.neigh]) / 3.0
        np.add.at(V, self.neigh, pvn)
        np.add.at(C, self.neigh, pvn[:, None] * (0.75 * Cf[:nin] + 0.25 * cEst[self.neigh]))
        self.V = V
        self.C = C / np.where(np.abs(V) > 0, V, 1.0)[:, None]

    def patch(self, name):
        b = self.boundary[name]
        s, n = b["startFace"], b["nFaces"]
        return slice(s, s + n), self.owner[s:s + n]


def latest_times(case, k=2):
    ts = []
    for d in os.listdir(case):
        try:
            t = float(d)
        except ValueError:
            continue
        if t > 0 and os.path.isdir(os.path.join(case, d)):
            ts.append((t, d))
    return [d for _, d in sorted(ts)[-k:]]


_MW = None


def mol_weights():
    """молярные массы веществ механизма (thermo.compressibleGasGRI) и состав по атомам"""
    global _MW
    if _MW is None:
        t = open(os.path.join(FOAM_RES, "thermo.compressibleGasGRI")).read()
        mw, el = {}, {}
        for m in re.finditer(r"\n([A-Za-z0-9()\-,]+)\s*\n\{(.*?)\n\}", t, flags=re.S):
            body = m.group(2)
            w = re.search(r"molWeight\s+([0-9.eE+-]+);", body)
            if w:
                mw[m.group(1)] = float(w.group(1))
            e = re.search(r"elements\s*\{([^}]*)\}", body)
            if e:
                el[m.group(1)] = {k.upper(): float(v) for k, v in re.findall(r"(\w+)\s+([0-9.]+);", e.group(1))}
        _MW = (mw, el)
    return _MW


class Case:
    def __init__(self, case, time=None):
        self.case = case
        self.mesh = Mesh(case)
        ts = latest_times(case, 2)
        self.time = time or (ts[-1] if ts else None)
        self.prev = ts[0] if len(ts) == 2 and time is None else None
        self.desc = json.load(open(os.path.join(case, "case.json"), encoding="utf-8"))
        self._cache = {}

    def field(self, name, time=None):
        t = time or self.time
        key = (name, t)
        if key not in self._cache:
            p = os.path.join(self.case, t, name)
            self._cache[key] = read_field(p) if os.path.exists(p) else None
        return self._cache[key]

    def species(self, time=None):
        mw, _ = mol_weights()
        t = time or self.time
        return [s for s in mw if os.path.exists(os.path.join(self.case, t, s))]

    def cellY(self, time=None):
        out = {}
        for s in self.species(time):
            f = self.field(s, time)
            v = f["internal"]
            out[s] = np.full(self.mesh.nc, v) if np.isscalar(v) else v
        return out

    def patch_vals(self, name, patch, time=None):
        """значения поля на гранях участка: записанные или пристенной ячейки"""
        f = self.field(name, time)
        sl, own = self.mesh.patch(patch)
        n = sl.stop - sl.start
        if f is None:
            return None
        b = f["boundary"].get(patch, {})
        v = b.get("value")
        if v is not None:
            return np.broadcast_to(v, (n, f["ncomp"])) .copy() if f["ncomp"] == 3 else np.broadcast_to(v, (n,)).copy()
        iv = f["internal"]
        if np.isscalar(iv):
            return np.full(n, iv)
        return iv[own].copy()


def x_from_Y(Y):
    mw, _ = mol_weights()
    n = {s: Y[s] / mw[s] for s in Y}
    tot = sum(n.values())
    return {s: v / tot for s, v in n.items()}, 1.0 / tot


def phi_atoms(Y):
    """phi по атомам: (2·C + H/2)/O - годится для свежей, горящей и сгоревшей смеси"""
    _, el = mol_weights()
    mw, _ = mol_weights()
    nC = nH = nO = 0.0
    for s, y in Y.items():
        e = el.get(s, {})
        k = y / mw[s]
        nC = nC + k * e.get("C", 0)
        nH = nH + k * e.get("H", 0)
        nO = nO + k * e.get("O", 0)
    return (2 * nC + 0.5 * nH) / np.maximum(nO, 1e-30)


def emissions(X, O2_ref=15.0):
    """мольные доли (влажные) -> сухие ppm при 15 % O2 (NOx = NO + NO2)"""
    xh2o = X.get("H2O", 0.0)
    dry = 1.0 - xh2o
    o2d = 100.0 * X.get("O2", 0.0) / dry
    corr = (20.9 - O2_ref) / max(20.9 - o2d, 1e-6)
    ppm = lambda s: 1e6 * X.get(s, 0.0) / dry * corr
    return dict(O2_dry_pct=o2d, NOx_ppm15=ppm("NO") + ppm("NO2"), NO_ppm15=ppm("NO"), NO2_ppm15=ppm("NO2"),
                CO_ppm15=ppm("CO"), CH4_ppm15=ppm("CH4"), corr15=corr, note="ppm сухого газа при 15 % O2")


def _patch_state(cs, patch, time=None):
    phi = cs.patch_vals("phi", patch, time)
    T = cs.patch_vals("T", patch, time)
    p = cs.patch_vals("p", patch, time)
    U = cs.patch_vals("U", patch, time)
    Y = {s: cs.patch_vals(s, patch, time) for s in cs.species(time)}
    return phi, T, p, U, Y


def mixed_mean(cs, patches, time=None):
    """средние по расходу через участки: T, состав (массовые доли), p0"""
    W, sT, sY, sp0, n = 0.0, 0.0, {}, 0.0, 0
    mw, _ = mol_weights()
    for pt in patches:
        phi, T, p, U, Y = _patch_state(cs, pt, time)
        if phi is None:
            continue
        w = np.abs(phi)
        Wt = w.sum()
        if Wt <= 0:
            continue
        inv_mw = sum(Y[s] / mw[s] for s in Y)
        rho = p / (R_U * 1e3 * inv_mw * T)
        p0 = p + 0.5 * rho * np.einsum("ij,ij->i", U, U)
        W += Wt
        sT += (w * T).sum()
        sp0 += (w * p0).sum()
        for s in Y:
            sY[s] = sY.get(s, 0.0) + (w * Y[s]).sum()
        n += len(w)
    if W <= 0:
        return None
    return dict(m=W, T=sT / W, Y={s: v / W for s, v in sY.items()}, p0=sp0 / W, faces=n)


def outlet_metrics(cs, T3, lhv, fuel_in_frac, time=None):
    phi, T, p, U, Y = _patch_state(cs, "outlet", time)
    w = np.clip(phi, 0, None)
    W = w.sum()
    if W <= 0:
        return dict(m=0.0, status="нет расхода через выход")
    Ym = {s: float((w * Y[s]).sum() / W) for s in Y}
    Tm = float((w * T).sum() / W)
    X, mwm = x_from_Y(Ym)
    em = emissions(X)
    q_left = sum(Ym.get(s, 0.0) * lhv.get(s, 0.0) for s in ("CH4", "CO", "H2"))
    eta = 1.0 - q_left / max(fuel_in_frac * lhv["CH4"], 1e-30)
    # радиальный профиль: грани выхода - радиальные полосы клина
    r = np.hypot(cs.mesh.Cf[cs.mesh.patch("outlet")[0], 0], 0.0)
    order = np.argsort(r)
    prof = dict(r=r[order].tolist(), T=T[order].tolist(), w=w[order].tolist())
    rtdf = (float(T.max()) - Tm) / max(Tm - T3, 1e-9)
    return dict(m=float(W), T=Tm, T_max_face=float(T.max()), RTDF=rtdf, Y_major={s: Ym[s] for s in
                ("CH4", "O2", "CO2", "H2O", "CO", "H2", "NO", "NO2", "OH") if s in Ym},
                eta_c=float(eta), MW=float(mwm), profile=prof, **em)


def convergence(cs, T3, lhv, fuel_in_frac):
    if not cs.prev:
        return dict(ok=None, note="одна запись - сходимость не оценена")
    a = outlet_metrics(cs, T3, lhv, fuel_in_frac, cs.prev)
    b = outlet_metrics(cs, T3, lhv, fuel_in_frac, cs.time)
    if "T" not in a or "T" not in b:
        return dict(ok=False, note="нет расхода через выход в одной из записей", times=[cs.prev, cs.time])
    rel = lambda x, y, floor: abs(x - y) / max(abs(y), floor)
    d = dict(dT=rel(a["T"], b["T"], 1.0), dNOx=rel(a["NOx_ppm15"], b["NOx_ppm15"], 0.5),
             dCO=rel(a["CO_ppm15"], b["CO_ppm15"], 1.0), times=[cs.prev, cs.time])
    d["ok"] = bool(d["dT"] < 0.002 and d["dNOx"] < 0.05 and d["dCO"] < 0.10)
    d["rule"] = "между двумя последними записями: T < 0,2 %, NOx < 5 % (не меньше 0,5 ppm), CO < 10 % (не меньше 1 ppm)"
    return d


def log_info(case):
    out = {}
    for f in sorted(glob.glob(os.path.join(case, "log.reactingFoam*"))):
        t = open(f, errors="replace").read()
        et = re.findall(r"ExecutionTime = ([0-9.]+) s", t)
        res = {}
        for m in re.finditer(r"Solving for (\w+), Initial residual = ([0-9.eE+-]+)", t[-20000:]):
            res[m.group(1)] = float(m.group(2))
        out[os.path.basename(f)] = dict(ended=t.rstrip().endswith("End"), exec_s=float(et[-1]) if et else None,
                                        fatal="FOAM FATAL" in t, last_initial_residual=res)
    return out


def cell_metrics(cs):
    d = cs.desc
    C = cs.mesh.C
    z, r = C[:, 2], np.hypot(C[:, 0], 0.0)
    T = cs.field("T")["internal"]
    U = cs.field("U")["internal"]
    Y = cs.cellY()
    zd, rch = d["z_dump"], d["geometry"]["r_ch"]
    out = {}
    _ = r
    ch = z < zd
    out["T_max_channel"] = float(T[ch].max()) if ch.any() else None
    out["flashback_or_autoignition"] = bool(out["T_max_channel"] is not None and out["T_max_channel"] > d["T3"] + 200.0)
    # слой ячеек на срезе канала (последний ряд перед срезом)
    if ch.any():
        zl = z[ch].max()
        lay = ch & (np.abs(z - zl) < 1e-9 + 1e-6 * zd)
        ph = phi_atoms({s: v[lay] for s, v in Y.items()})
        m = np.abs(U[lay, 2]) * cs.mesh.V[lay]          # расход через слой ~ u_z·V (ячейки слоя одной толщины)
        mean = float((m * ph).sum() / m.sum())
        sd = float(math.sqrt(max(0.0, (m * (ph - mean) ** 2).sum() / m.sum())))
        out.update(phi_exit_mean=mean, S_exit=sd / max(mean, 1e-12), S_exit_note="σ/среднее phi по атомам на срезе канала, вес u_z·V")
    # обратный ток на 2 калибрах за срезом
    z2 = zd + 2 * 2 * rch
    if z2 < d["z_end"]:
        zc = z[np.argmin(np.abs(z - z2))]
        lay = np.abs(z - zc) < 1e-9 + 1e-6 * zd
        mw, _ = mol_weights()
        inv_mw = sum(Y[s][lay] / mw[s] for s in Y)
        pc = cs.field("p")["internal"]
        pc = np.full(lay.sum(), pc) if np.isscalar(pc) else pc[lay]
        rho = pc / (R_U * 1e3 * inv_mw * T[lay])
        g = rho * U[lay, 2] * cs.mesh.V[lay]                 # расход через слой ~ ρ·u_z·V (слой одной толщины)
        out["recirc_ratio_2cal"] = float(-g[g < 0].sum() / max(g[g > 0].sum(), 1e-30))
        out["recirc_note"] = "обратный расход / прямой через слой ячеек на 2 калибрах канала за срезом"
    hot = T > d["T3"] + 400.0
    out["V_hot_frac"] = float(cs.mesh.V[hot].sum() / cs.mesh.V.sum())
    out["T_max"] = float(T.max())
    # масса газа в области и в горячей зоне (клин) -> время пребывания для сверки сетью реакторов Cantera
    mw, _ = mol_weights()
    inv_mw = sum(Y[s] / mw[s] for s in Y)
    pc = cs.field("p")["internal"]
    pc = np.full(cs.mesh.nc, pc) if np.isscalar(pc) else pc
    rho = pc / (R_U * 1e3 * inv_mw * T)
    out["mass_domain"] = float((rho * cs.mesh.V).sum())
    out["mass_hot"] = float((rho * cs.mesh.V)[hot].sum())
    out["mass_downstream"] = float((rho * cs.mesh.V)[~ch].sum())
    return out


def liner_metrics(cs):
    out = {}
    T = cs.field("T")["internal"]
    for wp in [p for p in cs.mesh.boundary if p.startswith("wall")]:
        sl, own = cs.mesh.patch(wp)
        zc = cs.mesh.Cf[sl, 2]
        Tw = T[own]
        o = np.argsort(zc)
        out[wp] = dict(T_aw_max=float(Tw.max()), z_at_max=float(zc[np.argmax(Tw)]), T_aw_mean=float(Tw.mean()),
                       z=zc[o].tolist(), T_aw=Tw[o].tolist())
    out["T_aw_max"] = max(v["T_aw_max"] for v in out.values() if isinstance(v, dict)) if out else None
    return out


def analyse(case, inputs):
    """итоги кейса; inputs: dict(T3, lhv{CH4,CO,H2}, fuel_in_frac) из описания варианта"""
    cs = Case(case)
    d = cs.desc
    res = dict(case=os.path.basename(case.rstrip("/")), kind=d["kind"], load=d["load"], level=d["level"],
               n_cells=d["n_cells"], time=cs.time, log=log_info(case))
    if cs.time is None:
        res["status"] = "нет записи"
        return res
    inl = [p for p in cs.mesh.boundary if p.startswith("inlet")]
    mi = mixed_mean(cs, inl)
    mo = mixed_mean(cs, ["outlet"])
    res["mass_in"], res["mass_out"] = mi["m"], mo["m"]
    res["mass_imbalance"] = (mi["m"] - mo["m"]) / mi["m"]
    res["p0_in"], res["p0_out"] = mi["p0"], mo["p0"]
    res["dp0_rel"] = (mi["p0"] - mo["p0"]) / mi["p0"]
    fin = mi["Y"].get("CH4", 0.0)
    res["fuel_in_frac"] = fin
    res["outlet"] = outlet_metrics(cs, inputs["T3"], inputs["lhv"], fin)
    res["outlet_mixed_state"] = dict(T=mo["T"], Y={s: v for s, v in mo["Y"].items() if v > 1e-7},
                                     m=mo["m"], Uz=None)
    res["convergence"] = convergence(cs, inputs["T3"], inputs["lhv"], fin)
    res.update(cell_metrics(cs) if d["kind"] == "cell" else liner_metrics(cs))
    fatal = any(v["fatal"] for v in res["log"].values())
    res["status"] = "ошибка решателя" if fatal else ("сошлось" if res["convergence"].get("ok") else "не сошлось")
    return res


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("case")
    ap.add_argument("--T3", type=float, required=True)
    ap.add_argument("--lhv", default="CH4=50.0e6,CO=10.1e6,H2=120.0e6",
                    help="теплота сгорания, Дж/кг (по умолчанию - для проверки разбора; в расчёте - из входов варианта)")
    a = ap.parse_args(argv)
    lhv = {k: float(v) for k, v in (x.split("=") for x in a.lhv.split(","))}
    r = analyse(a.case, dict(T3=a.T3, lhv=lhv))
    r["outlet"].pop("profile", None)
    print(json.dumps(r, ensure_ascii=False, indent=1, default=float)[:5000])


if __name__ == "__main__":
    main()
