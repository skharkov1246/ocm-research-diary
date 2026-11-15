#!/usr/bin/env python3
"""Проверка элит семейства «лопатки и колёса» прямым расчётом CalculiX и калибровка быстрой оценки по зонам.

  * Колесо компрессора: осесимметричная модель routes/gtu/platform/fea/disk.py (CAX8, CENTRIF, поле T(r), лопатки —
    массой в плотности столбца, как frame3 и organic) по профилю элиты: тракт (линия втулки) и задняя поверхность,
    у ступени 1 — полотно между двумя трактами; «рыбья кость» — эквивалентной толщиной (как быстрая оценка).
  * Диск турбины (блиск): то же, материал IN-713C, нагрузка лопаток — давлением на ободе, ресурс по Ларсону — Миллеру.
  * Из узлов .frd: σ_Мизес/σ0,2(T) по зонам «у оси» (r < 0,25·r2) и «полотно»; запас по разрушению и проверка
    равновесия — из disk.py.

Калибровка (--calibrate): набор профилей ступеней 1 и 2 (облик frame3 и варианты задней поверхности, ступицы, втулки)
→ k_axis, k_web = средние отношения ccx / тонкий диск по зонам; пишется в disk1d_calibration.json → zones.
Команды:
  python3 ccx_elites.py --calibrate [--n 12] [--mesh 80 12]
  python3 ccx_elites.py --elites elites_comp2.json --out DIR
"""
from __future__ import annotations

import argparse
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
import disk1d as D1        # noqa: E402
DK = D1.DK
WORK = os.environ.get("EVO_OUT", os.path.join(os.environ.get("TMPDIR", "/tmp"), "evo_blades"))


def comp_case(comp, g, rpm_factor=1.0, mesh=(80, 12), tag="c"):
    """case для disk.py по геному колеса; возврат (case, данные быстрой оценки)"""
    import blade_eval as BE
    r = BE.evaluate(comp, np.asarray(g, float), want_geom=True)
    geo = r["geom"]
    dk = geo["disk"]
    rr = dk["r"]
    stage = 1 if comp == "comp1" else 2
    if stage == 1:
        xc = dk["Lz"] + dk["t_rim"] / 2
        zc = np.full_like(rr, xc)
    else:
        zc = dk["xh"] + 0.5 * dk["h"]
    t = dk["h_eff"]
    idx = np.unique(np.r_[np.linspace(0, len(rr) - 1, 61).astype(int)])
    P = BE.point()["stage%d" % stage]
    T = BE.comp_temperature(stage, rr, rr[-1])
    case = dict(name=f"{tag}_{comp}", material="ti64", rpm=P["n"] * rpm_factor,
                profile=dict(r=rr[idx].tolist(), t=t[idx].tolist(), zc=zc[idx].tolist(), hb=geo["hb"][idx].tolist()),
                profile_interp="pchip", blade_fill=float(geo["blade_fill_phi"]), T_C=dict(r=rr.tolist(), T=T.tolist()),
                burst_k=1.0, mesh=dict(nr=mesh[0], nz=mesh[1]))
    return case, r


def columns(case, res_dir):
    """наибольшее σ_Мизес/σ0,2(T) по узлам каждого столбца сетки (одинаковый радиус): (r, yr)"""
    mat = DK.load_material(case["material"])
    prof = DK.Profile(case["profile"], case.get("profile_interp", "linear"))
    mesh = DK.build_mesh(prof, case["mesh"]["nr"], case["mesh"]["nz"])
    frd = DK.parse_frd(os.path.join(res_dir, "disk.frd"))
    S, TN = frd["STRESS"], frd.get("NDTEMP", {})
    ys = np.asarray(mat["ys_MPa"], float)
    col, colr = {}, {}
    for n, (r, z, i, j) in mesh["nodes"].items():
        if n not in S:
            continue
        s = S[n]
        vm = DK.von_mises(*s[:6])
        Tn = TN[n][0] if n in TN else DK.temperature_field(case, r)
        yr = vm / 1e6 / float(np.interp(Tn, ys[:, 0], ys[:, 1]))
        if yr > col.get(i, -1):
            col[i], colr[i] = yr, r
    ii = sorted(col)
    return np.array([colr[i] for i in ii]), np.array([col[i] for i in ii])


X_BINS = [0.0, 0.1, 0.2, 0.35, 0.5, 0.75, 1.0, 1.5, 2.0, 3.0, 5.0, 1e9]
YR_CAL_MAX = 1.0           # таблица — по профилям вблизи допустимости (σ/σ0,2 по ccx ≤ 1): там нужна точность
TABLE_Q = float(os.environ.get("EVO_TABLE_Q", "0.75"))   # квантиль отношения в ячейке h/r (1 — огибающая по максимуму)


def zone_max(case, res_dir, r_split):
    """σ_Мизес/σ0,2(T) по узлам .frd в двух зонах"""
    mat = DK.load_material(case["material"])
    prof = DK.Profile(case["profile"], case.get("profile_interp", "linear"))
    mesh = DK.build_mesh(prof, case["mesh"]["nr"], case["mesh"]["nz"])
    frd = DK.parse_frd(os.path.join(res_dir, "disk.frd"))
    S, TN = frd["STRESS"], frd.get("NDTEMP", {})
    ys = np.asarray(mat["ys_MPa"], float)
    out = {"axis": (0.0, None), "web": (0.0, None)}
    for n, (r, z, i, j) in mesh["nodes"].items():
        if n not in S:
            continue
        s = S[n]
        vm = DK.von_mises(s[0], s[1], s[2], s[3], s[4] if len(s) > 4 else 0.0, s[5] if len(s) > 5 else 0.0)
        Tn = TN[n][0] if n in TN else DK.temperature_field(case, r)
        yr = vm / 1e6 / float(np.interp(Tn, ys[:, 0], ys[:, 1]))
        zk = "axis" if r < r_split else "web"
        if yr > out[zk][0]:
            out[zk] = (yr, r)
    return dict(yr_axis=out["axis"][0], r_axis=out["axis"][1], yr_web=out["web"][0], r_web=out["web"][1])


def run_comp(job):
    comp, g, tag, mesh, rpm_factor, keep = job
    t0 = time.time()
    case, r = comp_case(comp, g, rpm_factor, mesh, tag)
    wd = os.path.join(WORK, "ccx", tag)
    res = DK.run_case(case, workdir=wd, keep=True)
    rs = D1.R_SPLIT_FRAC * case["profile"]["r"][-1]
    zm = zone_max(case, wd, rs)
    # быстрая оценка того же профиля (k = 1 — сырые зоны)
    import blade_eval as BE
    geo = r["geom"]
    dk = geo["disk"]
    fac = 1 + geo["blade_fill_phi"] * geo["hb"] / np.maximum(dk["h_eff"], 1e-6)
    stage = 1 if comp == "comp1" else 2
    T = BE.comp_temperature(stage, dk["r"], dk["r"][-1])
    fast = D1.solve(dk["r"], dk["h_eff"], fac, T, "ti64", case["rpm"], r_split=rs, return_field=True)
    rc, yc = columns(case, wd)
    F = fast.pop("field")
    o = np.argsort(F["r"])
    fy = np.interp(rc, F["r"][o], F["yr_raw"][o])
    hx = np.interp(rc, F["r"][o], F["h"][o]) / np.maximum(rc, 1e-3)
    col = dict(r=rc.tolist(), x=hx.tolist(), ccx=yc.tolist(), fast=fy.tolist())
    kt = D1.k_table_cal(comp)
    fast_cal = None
    if kt is not None:
        fast_cal = D1.solve(dk["r"], dk["h_eff"], fac, T, "ti64", case["rpm"], k_table=kt)["yield_ratio_max"]
    if not keep:
        import shutil
        shutil.rmtree(wd, ignore_errors=True)
    return dict(tag=tag, comp=comp, rpm=case["rpm"], ccx=dict(yield_ratio_max=res["yield_ratio_max"], yield_at=res["yield_ratio_at"],
                burst_margin=res["burst"]["margin"] if isinstance(res.get("burst"), dict) else res.get("burst_margin"),
                gate=res["quality"]["passed"], seconds=res["run"]["seconds"], **zm),
                fast=dict(yr_axis_raw=fast["yr_axis_raw"], yr_web_raw=fast["yr_web_raw"], burst_margin=fast["burst_margin"],
                          r_web_peak=fast["r_web_peak"], yield_ratio_cal=fast_cal), columns=col,
                mass_kg=r["meta"].get("mass_kg"), feasible_fast=r["feasible"], seconds=round(time.time() - t0, 1))


def calibration_set(n=12, n_out=4):
    """профили для калибровки: облик frame3 обеих ступеней и варианты (задняя поверхность тоньше/толще, ступица,
    втулка, «рыбья кость») — геномы детерминированы"""
    import blade_eval as BE
    import blade_geom as BG
    import core_local as C
    out = []
    n_o = 0
    for comp in ("comp1", "comp2"):
        x0 = BE.baseline_genes(comp)
        out.append((comp, x0, f"cal_{comp}_base"))
    names = [g.name for g in BG.genes_for("comp2")]
    x0 = BE.baseline_genes("comp2")
    rng = C.rng_of("blades-calibration", n)
    k = 0
    while len(out) < n:
        comp = "comp2" if k % 3 else "comp1"
        x = BE.baseline_genes(comp).copy()
        gn = [g.name for g in BG.genes_for(comp)]
        for nm in ("hub_a", "hub_c", "sh_a", "sh_c", "t_hub_k", "t_tip_k"):
            x[gn.index(nm)] = np.clip(x[gn.index(nm)] + 0.35 * rng.standard_normal(), 0, 1)
        if comp == "comp2":
            z = rng.uniform(0.1, 0.9)
            for i in range(6):
                x[names.index(f"zb{i}")] = np.clip(z + 0.15 * rng.standard_normal(), 0, 1)
            x[names.index("boss")] = rng.uniform(0, 1)
            x[names.index("fish")] = rng.uniform(0, 1)
        r = BE.evaluate(comp, x, want_geom=True)
        sl = r["geom"]["disk"]["slope_max"]
        dok = r["geom"]["disk"]["ok"]
        if dok and len([o for o in out if not o[2].startswith("out_")]) < n:
            out.append((comp, x, f"cal_{comp}_{k:02d}"))
        elif not dok and n_o < n_out:
            out.append((comp, x, f"out_{comp}_{k:02d}"))
            n_o += 1
        k += 1
        if k > 2000 or (len(out) >= n + n_out):
            break
    return out


def calibrate(n=12, mesh=(80, 12), workers=4):
    cset = calibration_set(n)
    jobs = [(c, x.tolist(), tag, mesh, 1.0, False) for c, x, tag in cset]
    with Pool(min(workers, len(jobs))) as pool:
        rows = pool.map(run_comp, jobs)
    def build_table(comp_):
        X, K = [], []
        for r_ in rows:
            if (r_["comp"] != comp_ or not r_["ccx"]["gate"] or r_["tag"].startswith("out_")
                    or r_["ccx"]["yield_ratio_max"] > YR_CAL_MAX):
                continue
            cc = r_["columns"]
            for x, a, b, rv in zip(cc["x"], cc["ccx"], cc["fast"], cc["r"]):
                if a >= 0.15 and b > 0 and rv < 0.97 * cc["r"][-1]:      # у обода напряжения малы — отношение неустойчиво
                    X.append(x)
                    K.append(a / b)
        X, K = np.array(X), np.array(K)
        xs, ks, ns = [], [], []
        for lo_, hi_ in zip(X_BINS[:-1], X_BINS[1:]):
            m = (X >= lo_) & (X < hi_)
            if m.any():
                xs.append(float(np.median(X[m])))
                ks.append(float(np.quantile(K[m], TABLE_Q)))
                ns.append(int(m.sum()))
        ks = np.maximum.accumulate(np.array(ks)).tolist()        # не убывает с толщиной (огибающая)
        return dict(x=[0.0] + xs + [1e6], k=[ks[0]] + ks + [ks[-1]], n_points=ns, bins=X_BINS[:-1],
                    n_profiles=len({r_["tag"] for r_ in rows if r_["comp"] == comp_ and not r_["tag"].startswith("out_")
                                    and r_["ccx"]["yield_ratio_max"] <= YR_CAL_MAX}),
                    quantile=TABLE_Q,
                    method="квантиль TABLE_Q отношения ccx/тонкий диск по столбцам сетки в ячейках местного h/r; не убывает с h/r")
    table = {cm: build_table(cm) for cm in ("comp1", "comp2")}
    ka = np.array([r_["ccx"]["yr_axis"] / r_["fast"]["yr_axis_raw"] for r_ in rows if r_["ccx"]["gate"] and not r_["tag"].startswith("out_")])
    kw = np.array([r_["ccx"]["yr_web"] / r_["fast"]["yr_web_raw"] for r_ in rows if r_["ccx"]["gate"] and not r_["tag"].startswith("out_")])
    zones = dict(method=("CalculiX disk.py (CAX8 %d×%d, 100 %% частоты) против тонкого диска на %d профилях: облик frame3 и варианты "
                         "задней поверхности, ступицы, втулки, «рыбьей кости»; k — наибольшее из отношений по зонам "
                         "(консервативно), среднее и разброс — для сведения") % (mesh[0], mesh[1], len(rows)),
                 r_split_frac=D1.R_SPLIT_FRAC, k_axis=float(ka.max()), k_web=float(kw.max()), k_axis_mean=float(ka.mean()),
                 k_web_mean=float(kw.mean()), k_axis_range=[float(ka.min()), float(ka.max())], k_web_range=[float(kw.min()), float(kw.max())],
                 cases=rows, date=time.strftime("%Y-%m-%d"))
    d = json.load(open(D1.CAL_JSON, encoding="utf-8")) if os.path.exists(D1.CAL_JSON) else D1.calibrate(write=False)
    for r_ in zones["cases"]:
        r_.pop("columns", None)
    d["zones"] = zones
    d["table"] = table
    with open(D1.CAL_JSON, "w", encoding="utf-8") as f:
        json.dump(d, f, ensure_ascii=False, indent=1, default=float)
    # проверка таблицы на тех же профилях: max_r k(h/r)·быстрая против max ccx
    chk = []
    for c_, x_, tag in cset:
        kt = (np.asarray(table[c_]["x"]), np.asarray(table[c_]["k"]))
        case, r = comp_case(c_, x_, 1.0, mesh, tag)
        geo = r["geom"]
        dk = geo["disk"]
        fac = 1 + geo["blade_fill_phi"] * geo["hb"] / np.maximum(dk["h_eff"], 1e-6)
        T = __import__("blade_eval").comp_temperature(1 if c_ == "comp1" else 2, dk["r"], dk["r"][-1])
        f_ = D1.solve(dk["r"], dk["h_eff"], fac, T, "ti64", case["rpm"], k_table=kt)["yield_ratio_max"]
        ccx_ = [q for q in rows if q["tag"] == tag][0]["ccx"]["yield_ratio_max"]
        sl = geo["disk"]["slope_max"]
        raw = [q for q in rows if q["tag"] == tag][0]["fast"]
        chk.append(dict(tag=tag, ccx=ccx_, fast_table=f_, fast_raw=max(raw["yr_axis_raw"], raw["yr_web_raw"]), ratio=ccx_ / f_, slope_max=sl,
                        in_domain=not tag.startswith("out_"), used_for_table=(not tag.startswith("out_")) and ccx_ <= YR_CAL_MAX,
                        verdict_ccx=bool(ccx_ <= 0.9), verdict_fast=bool(f_ <= 0.9)))
    ind = [q for q in chk if q["in_domain"]]
    d["domain"] = dict(slope_max_limit=__import__("blade_geom").SLOPE_MAX, rmin_over_h=__import__("blade_geom").RMIN_H,
                       rule=("|dh/dr| ≤ %g на 0,1…0,95·r2; задняя поверхность монотонна (глубина убывает к ободу), R_min/h ≥ %g"
                             % (__import__("blade_geom").SLOPE_MAX, __import__("blade_geom").RMIN_H)),
                       out_of_domain=[q for q in chk if not q["in_domain"]],
                       note="вне области (крутой переход ступицы в тонкое полотно) отношение ccx/тонкий диск неустойчиво — такие профили недопустимы в поиске")
    agree = [q["verdict_ccx"] == q["verdict_fast"] for q in ind]
    unsafe = [q["tag"] for q in ind if q["verdict_fast"] and not q["verdict_ccx"]]
    d["table_check"] = dict(cases=chk, ratio_min=min(q["ratio"] for q in ind), ratio_max=max(q["ratio"] for q in ind),
                            verdict_threshold="σ/σ0,2 ≤ 0,9 при 100 % (порядок предела 120 %/k: 0,9·σ0,2 при k = 0,9 и запасе 1,2² — справочно)",
                            verdict_agree_share=float(np.mean(agree)), unsafe_cases=unsafe,
                            note="ratio ≤ 1 — быстрая оценка с таблицей не ниже CalculiX (консервативна) на калибровочных профилях")
    with open(D1.CAL_JSON, "w", encoding="utf-8") as f:
        json.dump(d, f, ensure_ascii=False, indent=1, default=float)
    zones["table"] = table
    zones["table_check"] = d["table_check"]
    return zones


# =========================================================================== элиты


def turb_cases(g, mesh=(80, 12), tag="t"):
    import blade_eval as BE
    r = BE.evaluate("turb", np.asarray(g, float), want_geom=True)
    out = []
    for k, (bl, dk) in enumerate(zip(r["geom"]["blades"], r["geom"]["disks"]), start=1):
        idx = np.linspace(0, len(dk["r"]) - 1, 61).astype(int)
        case = dict(name=f"{tag}_turb{k}", material="in713c", rpm=BE.point()["turb"]["n"],
                    profile=dict(r=dk["r"][idx].tolist(), t=dk["h_eff"][idx].tolist()), profile_interp="pchip",
                    T_C=dict(r=dk["r"].tolist(), T=dk["T"].tolist()), rim=dict(p_MPa=dk["p_rim_MPa"]), burst_k=1.0,
                    mesh=dict(nr=mesh[0], nz=mesh[1]))
        out.append((case, dk))
    return out, r


def run_turb(job):
    g, tag, mesh = job
    cases, r = turb_cases(g, mesh, tag)
    rows = []
    for case, dk in cases:
        wd = os.path.join(WORK, "ccx", case["name"])
        res = DK.run_case(case, workdir=wd, keep=False)
        creep = res.get("creep") if isinstance(res.get("creep"), dict) else {}
        rows.append(dict(stage=case["name"][-1], ccx_burst=res["burst"]["margin"] if isinstance(res.get("burst"), dict) else res.get("burst_margin"),
                         ccx_life_h=creep.get("local_min_h") if creep else None, ccx_gate=res["quality"]["passed"],
                         ccx_vm_max_MPa=res.get("vm_max_node_MPa"), fast_burst=dk["burst"], fast_life_h=dk["life"]))
    return dict(tag=tag, stages=rows)


def check_file(path, out_dir, workers=4, mesh=(80, 12)):
    """проверка элит из файла evo_blades.py elites/part: колёса — disk.py при 100 и 120 %, турбина — диски блисков"""
    os.makedirs(out_dir, exist_ok=True)
    d = json.load(open(path, encoding="utf-8"))
    comp = d["comp"]
    rows = []
    if comp == "turb":
        jobs = [(e["g"], e["id"], mesh) for e in d["elites"]]
        with Pool(max(1, min(workers, len(jobs)))) as pool:
            rows = pool.map(run_turb, jobs)
    else:
        jobs = []
        for e in d["elites"]:
            for f in (1.0, 1.2):
                jobs.append((comp, e["g"], f"{e['id']}_{int(f * 100)}", mesh, f, False))
        with Pool(max(1, min(workers, len(jobs)))) as pool:
            rows = pool.map(run_comp, jobs)
        for r_ in rows:
            r_["fast_cal"] = r_["fast"]["yield_ratio_cal"]
            r_["ratio_ccx_to_fast_cal"] = r_["ccx"]["yield_ratio_max"] / r_["fast_cal"] if r_["fast_cal"] else None
            r_.pop("columns", None)
    summ = dict(comp=comp, n=len(rows))
    if comp != "turb" and rows:
        rr = np.array([r_["ratio_ccx_to_fast_cal"] for r_ in rows if r_["ratio_ccx_to_fast_cal"] is not None])
        summ.update(ratio_min=float(rr.min()), ratio_max=float(rr.max()), all_gates=all(r_["ccx"]["gate"] for r_ in rows),
                    conservative_share=float(np.mean(rr <= 1.0)))
    out = dict(source=os.path.basename(path), rows=rows, summary=summ)
    with open(os.path.join(out_dir, f"ccx_{os.path.basename(path)}"), "w", encoding="utf-8") as f:
        json.dump(out, f, ensure_ascii=False, indent=1, default=float)
    return out


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--calibrate", action="store_true")
    ap.add_argument("--n", type=int, default=12)
    ap.add_argument("--mesh", type=int, nargs=2, default=[80, 12])
    ap.add_argument("--elites", default=None)
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--out", default=os.path.join(WORK, "ccx_check"))
    a = ap.parse_args()
    if a.calibrate:
        z = calibrate(a.n, tuple(a.mesh), a.workers)
        print(json.dumps({k: v for k, v in z.items() if k != "cases"}, ensure_ascii=False, indent=1))
    if a.elites:
        o = check_file(a.elites, a.out, a.workers, tuple(a.mesh))
        print(json.dumps(o["summary"], ensure_ascii=False))
