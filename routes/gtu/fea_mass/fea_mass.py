#!/usr/bin/env python3
"""Массовая 3D-прочность вариантов уровня 1 воронки (уровень 2 по прочности): сектор каждой роторной детали в CalculiX -
статика 100/115/120 %, запас по разрушающей частоте, частоты с преднагружением и диаграмма Кэмпбелла, ползучесть по
Ларсону - Миллеру на 25 000 ч, МЦУ пуск - останов; вердикт варианта с причинами; слияние частей с контролем полноты;
возврат вердиктов к записям уровня 1 для отбора 50.

  python3 fea_mass.py prepare --l1-in <папка частей уровня 1> --kv KV-3 --m 3983 --shards 40 --out cases/fea3d_cases_KV-3_m3983.json
  python3 fea_mass.py prepare --top50 ../funnel/funnel_l1_top50.json --out cases/fea3d_cases_probe.json
  python3 fea_mass.py run --cases <файл случаев> --shards S --shard 0,7 --workers "$NPROC" --out "$GTU_OUT"
  python3 fea_mass.py merge --cases <файл случаев> --in <папки частей> --out-dir <папка>
  python3 fea_mass.py select --l1-in <части уровня 1> --kv KV-3 --m 3983 --shards 40 --verdicts <fea3d_verdicts.json> --out-dir <папка>
  python3 fea_mass.py probe            # миниатюра: 2 варианта funnel_l1_top50.json локально

Деление на части - одно, по месту варианта в файле случаев: вариант r -> часть r % S; все его детали - в той же части.
Внутри части - пул процессов по деталям (--workers; ccx в один поток): 8-16 вариантов на машине идут одновременно.
AWS этот скрипт не вызывает: задачи - jobs_fea_mass.json (make_jobs_fea_mass.py), запуск - routes/gtu/platform/aws/runner.py.
Все величины, кроме взятых из записей уровня 1 и расчёта ccx, - оценка (geom3d.A, parts3d, post3d.CRIT).
"""
import argparse
import glob
import gzip
import hashlib
import json
import math
import os
import re
import shutil
import sys
import tempfile
import time
import traceback
from concurrent.futures import ProcessPoolExecutor, as_completed

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
GTU = os.path.dirname(HERE)
FUNNEL = os.path.join(GTU, "funnel")
for p in (HERE,):
    if p not in sys.path:
        sys.path.insert(0, p)
import geom3d as G           # noqa: E402
import ccx3d as C            # noqa: E402
import materials3d as M      # noqa: E402
import parts3d as P          # noqa: E402
import post3d as Q           # noqa: E402

DATE = "2026-11-15"
PART_PREFIX = "fea3d"
LEVEL_TAG = "FEA3D"
CODE_FILES = ("geom3d.py", "ccx3d.py", "materials3d.py", "parts3d.py", "post3d.py", "fea_mass.py")
DEFAULTS = dict(level="base", static_fracs=[1.0, 1.15, 1.2], freq_fracs=[0.0, 0.5, 1.0, 1.2], nmodes=12, cold_zero=False,
                conv_every=0, conv_level="fine", timeout_s=7200)


def code_sha():
    h = hashlib.sha256()
    for f in CODE_FILES:
        h.update(open(os.path.join(HERE, f), "rb").read())
    return h.hexdigest()[:16]


def _sha(path):
    return hashlib.sha256(open(path, "rb").read()).hexdigest()


def _r(x, k=5):
    if isinstance(x, float):
        if not math.isfinite(x):
            return None if math.isnan(x) else (1e30 if x > 0 else -1e30)
        return float(f"{x:.{k}g}")
    if isinstance(x, (np.floating,)):
        return _r(float(x), k)
    if isinstance(x, (np.integer,)):
        return int(x)
    if isinstance(x, dict):
        return {str(a): _r(b, k) for a, b in x.items()}
    if isinstance(x, (list, tuple)):
        return [_r(b, k) for b in x]
    return x


# =========================================================================== поле температур по узлам
def temperature_nodes(mesh, part):
    r = np.hypot(mesh.nodes[:, 0], mesh.nodes[:, 1])
    ts = part["T_spec"]
    if ts["type"] == "quad":
        r_in = mesh.meta.get("r_bore", 0.0)
        r_out = ts["r_out"]
        x = np.clip((r - r_in) / max(r_out - r_in, 1e-9), 0.0, 1.0)
        return ts["T_in"] + (ts["T_out"] - ts["T_in"]) * x ** 2
    rb, R, rr = ts["r_bore"], ts["R"], ts["r_root"]
    T = np.empty(len(r))
    d = r < R
    T[d] = ts["T_hub"] + (ts["T_rim"] - ts["T_hub"]) * ((r[d] - rb) / max(R - rb, 1e-9)) ** 2
    a = (~d) & (r < rr)
    T[a] = ts["T_rim"] + (ts["T_blade"] - ts["T_rim"]) * (r[a] - R) / max(rr - R, 1e-9)
    T[r >= rr] = ts["T_blade"]
    return T


# =========================================================================== одна деталь
def build_mesh(part, level):
    if part["kind"] == "wheel":
        return G.wheel_sector(part["geom"], level=level)
    return G.axial_sector(part["geom"], level=level)


def part_materials(part):
    blade = M.for_part(part["alloy"])
    disc = M.for_part(part.get("disc_alloy") or part["alloy"])
    return {G.CLASS_HUB: disc, G.CLASS_BLADE: blade, G.CLASS_ATT: blade}


def _mat_cards(mats, mesh):
    present = sorted(set(int(c) for c in mesh.ecls))
    if all(mats[c].alloy == mats[present[0]].alloy and mats[c].rho == mats[present[0]].rho for c in present):
        return C.material_cards(mats[present[0]], "MAT", "EALL")
    out = []
    for c in present:
        out += C.material_cards(mats[c], f"M{G.CLASS_NAMES[c]}", f"E{G.CLASS_NAMES[c]}")
    return out


def run_part(part, workdir, opts, level=None, freq_fracs=None, static_fracs=None, keep=False):
    """деталь: сетка -> статика -> частоты -> критерии -> вердикт. Ошибка не роняет вариант: статус «ошибка» с причиной"""
    level = level or opts["level"]
    t0 = time.time()
    res = dict(id=part["id"], kind=part["kind"], role=part["role"], stage=part["stage"], rpm=part["rpm"], alloy=part["alloy"],
               disc_alloy=part.get("disc_alloy"), level=level, status="ok")
    os.makedirs(workdir, exist_ok=True)
    try:
        mesh = build_mesh(part, level)
    except Exception as e:                                   # noqa: BLE001
        res.update(status="ошибка", stage_failed="geometry", error=f"{type(e).__name__}: {e}"[:600], t_total_s=round(time.time() - t0, 1))
        return res
    meta = {k: v for k, v in mesh.meta.items() if k != "hub_polygon"}
    res["mesh"] = meta
    try:
        mats = part_materials(part)
        res["materials"] = {G.CLASS_NAMES[c]: m.describe() for c, m in mats.items() if c in set(int(x) for x in mesh.ecls)}
        T = temperature_nodes(mesh, part)
        mesh.write_inp(os.path.join(workdir, "mesh.inp"))
        cards = _mat_cards(mats, mesh)
        N = mesh.meta["Z"]
        excl = 2.0 * mesh.meta["r_bore"] if (part["kind"] == "wheel" and part["geom"].get("bore") == "core") else None
        res["T_field"] = dict(spec=part["T_spec"], T_min_C=float(T.min()), T_max_C=float(T.max()))
        # ------------------------------------------------ статика
        sf = static_fracs or opts["static_fracs"]
        st_run = C.run_ccx(workdir, "static", C.static_deck(N, cards, part["rpm"], sf, mesh.node_ids, T), timeout=opts["timeout_s"])
        res["ccx_static"] = st_run
        if not st_run["ok"]:
            raise RuntimeError(f"ccx статика: {st_run['errors'] or st_run['rc']}")
        steps = C.read_frd_steps(os.path.join(workdir, "static.frd"))
        if len(steps) < len(sf):
            raise RuntimeError(f"в static.frd {len(steps)} шагов из {len(sf)}")
        sres, ev100, S100 = [], None, None
        for k, f in enumerate(sf):
            ids, S = steps[k]["STRESS"]
            Sn = Q.align(mesh, ids, S[:, :6])
            ev_k, ev = Q.static_eval(mesh, mats, T, Sn, f, excl)
            ids_u, U = steps[k]["DISP"]
            Un = Q.align(mesh, ids_u, U[:, :3])
            ev_k["u_max_mm"] = float(np.nanmax(np.linalg.norm(Un, axis=1)) * 1e3)
            sres.append(ev_k)
            if k == 0:
                ev100, S100 = ev, Sn
            if abs(f - 1.2) < 1e-9:
                rob120 = Q.robinson(mesh, mats[G.CLASS_HUB], Sn, T)
        res["static"] = dict(fractions=sf, steps=sres)
        rob = Q.robinson(mesh, mats[G.CLASS_HUB], S100, T)
        if rob and rob.get("margin") is not None:
            rob = dict(rob, margin_120=rob["margin"] / 1.2,
                       margin_120_direct=(rob120 or {}).get("margin") if 1.2 in sf else None)
            dm = mats[G.CLASS_HUB]
            if dm.creep:
                lm, _ = dm.rupture_h(rob["hoop_mean_MPa"], rob["T_mean_C"])
                lh, _ = dm.rupture_h(rob["hoop_mean_MPa"], rob["T_hot_C"])
                rob.update(section_life_mean_h=lm, section_life_h=lh)
        res["burst"] = rob
        res["creep"] = Q.creep_eval(mesh, mats, ev100, T, mesh.meta)
        res["lcf"] = Q.lcf_eval(mats, sres[0])
        if part["alloy"] == "Si3N4":
            res["weibull"] = Q.weibull_eval(mesh, S100, N)
        if not keep:
            for f in ("static.frd",):
                try:
                    os.remove(os.path.join(workdir, f))
                except OSError:
                    pass
        # ------------------------------------------------ частоты
        ff = freq_fracs if freq_fracs is not None else opts["freq_fracs"]
        if ff:
            fr = C.run_ccx(workdir, "freq", C.freq_deck(N, cards, part["rpm"], ff, mesh.node_ids, T, opts["nmodes"], N // 2,
                                                       opts.get("cold_zero")), timeout=opts["timeout_s"])
            res["ccx_freq"] = fr
            if not fr["ok"]:
                raise RuntimeError(f"ccx частоты: {fr['errors'] or fr['rc']}")
            eig = C.read_eigen_steps(os.path.join(workdir, "freq.dat"))
            if opts.get("cold_zero"):
                res["freq_cold"] = {str(k): v for k, v in eig[0].items()}
                eig = eig[1:]
            if len(eig) != len(ff):
                raise RuntimeError(f"в freq.dat {len(eig)} таблиц частот из {len(ff)}")
            table = list(zip(ff, eig))
            res["campbell"] = Q.campbell(table, part["rpm"], N, part["eo_sources"], part["op"], part.get("free_stator"))
        v, hard, soft = Q.verdict(res)
        res.update(verdict=v, reasons_hard=hard, reasons_soft=soft)
    except Exception as e:                                   # noqa: BLE001
        res.update(status="ошибка", error=f"{type(e).__name__}: {e}"[:600], trace=traceback.format_exc()[-1500:])
    finally:
        if not keep:
            for f in glob.glob(os.path.join(workdir, "*")):
                if not f.endswith((".log", ".sta")):
                    try:
                        os.remove(f)
                    except OSError:
                        pass
    res["t_total_s"] = round(time.time() - t0, 1)
    res["cpu_ccx_s"] = round(sum((res.get(k) or {}).get("cpu_s", 0.0) for k in ("ccx_static", "ccx_freq")), 1)
    return _r(res)


def _conv(base, fine):
    """сходимость по двум сеткам: изменение ключевых величин (fine/base - 1)"""
    def g(r, path):
        x = r
        for k in path:
            if x is None:
                return None
            x = x[k] if isinstance(k, int) and isinstance(x, list) and k < len(x) else (x.get(k) if isinstance(x, dict) else None)
        return x
    out = {}
    for name, path in (("burst_margin", ("burst", "margin")), ("hoop_mean", ("burst", "hoop_mean_MPa")),
                       ("vm_hub_avg", ("static", "steps", 0, "regions", "hub", "vm_peak_avg_MPa")),
                       ("vm_blade_avg", ("static", "steps", 0, "regions", "blade", "vm_peak_avg_MPa")),
                       ("vm_blade_node", ("static", "steps", 0, "regions", "blade", "vm_peak_node_MPa"))):
        a, b = g(base, path), g(fine, path)
        if isinstance(a, (int, float)) and isinstance(b, (int, float)) and a:
            out[name] = b / a - 1
    fa = ((base.get("campbell") or {}).get("table") or [{}])
    fb = ((fine.get("campbell") or {}).get("table") or [{}])
    if fa and fb and fa[-1].get("nd") and fb[-1].get("nd"):
        devs = []
        for nd, v in fa[-1]["nd"].items():
            w = fb[-1]["nd"].get(nd) or []
            for x, y in zip(v[:3], w[:3]):
                devs.append(y / x - 1)
        if devs:
            out["freq_first3_max_abs"] = max(abs(d) for d in devs)
    return out


# =========================================================================== вариант
def _part_job(args):
    j, part, root, opts, mode = args
    wd = os.path.join(root, f"j{j}_{part['id']}_{mode}")
    try:
        if mode == "conv":
            r = run_part(part, wd, opts, level=opts["conv_level"], freq_fracs=[1.0], static_fracs=[1.0])
        else:
            r = run_part(part, wd, opts)
    finally:
        shutil.rmtree(wd, ignore_errors=True)
    return j, part["id"], mode, r


def variant_verdict(parts_res, pending):
    hard = [f"{p['id']}: {x}" for p in parts_res for x in p.get("reasons_hard") or []]
    soft = [f"{p['id']}: {x}" for p in parts_res for x in p.get("reasons_soft") or []]
    errs = [f"{p['id']}: {p.get('error')}" for p in parts_res if p.get("status") != "ok"]
    if hard:
        v = "отказ"
    elif errs:
        v = "ошибка"
    elif soft:
        v = "условно"
    elif not parts_res:
        v = "нет деталей"
    else:
        v = "годен"
    return v, hard, soft, errs


def summarize_variant(c, parts_res, pending, conv=None, t_s=None):
    v, hard, soft, errs = variant_verdict(parts_res, pending)

    def mins(fn):
        vals = [x for x in (fn(p) for p in parts_res) if isinstance(x, (int, float))]
        return min(vals) if vals else None
    b120 = mins(lambda p: (p.get("burst") or {}).get("margin_120"))
    yr = max([x for x in ((((p.get("static") or {}).get("steps") or [{}])[0]).get("yield_ratio_avg_max") for p in parts_res)
              if isinstance(x, (int, float))], default=None)
    life = mins(lambda p: min([x for x in [v_.get("life_avg_h") if k != "blade_root" else v_.get("life_h")
                                           for k, v_ in (p.get("creep") or {}).items()] if isinstance(x, (int, float))], default=None))
    lcf = mins(lambda p: min([x.get("cycles") for x in (p.get("lcf") or {}).values() if isinstance(x.get("cycles"), (int, float))],
                             default=None))
    nres = sum(1 for p in parts_res for x in ((p.get("campbell") or {}).get("crossings") or []) if x["in_range"])
    return _r(dict(j=c["j"], idx=c["idx"], kv=c.get("kv"), arch=c.get("arch"), niche=c.get("niche"), score_l1=c.get("score"),
                   feasible_l1=c.get("feasible"), verdict=v, reasons_hard=hard, reasons_soft=soft, errors=errs, pending=pending,
                   summary=dict(burst_min_120=b120, yield_ratio_max=yr, life_min_h=life, lcf_min_cycles=lcf, resonances_in_range=nres,
                                parts=len(parts_res), parts_ok=sum(1 for p in parts_res if p.get("status") == "ok")),
                   parts=parts_res, conv=conv, t_s=t_s))


def run_variants(cases, root, opts, workers=1, conv_every=0):
    """все детали всех вариантов - один пул процессов; итог - по вариантам в порядке cases"""
    jobs, meta = [], {}
    for r, c in enumerate(cases):
        parts, pending = P.variant_parts(c)
        meta[c["j"]] = dict(c=c, n=len(parts), pending=pending, order=[p["id"] for p in parts], t0=None)
        for p in parts:
            jobs.append((c["j"], p, root, opts, "main"))
            if conv_every and r % conv_every == 0:
                jobs.append((c["j"], p, root, opts, "conv"))
                meta[c["j"]]["conv"] = True
    got = {}
    t0 = time.time()
    if workers > 1 and len(jobs) > 1:
        with ProcessPoolExecutor(workers) as ex:
            futs = [ex.submit(_part_job, jb) for jb in jobs]
            for f in as_completed(futs):
                j, pid, mode, r = f.result()
                got[(j, pid, mode)] = r
    else:
        for jb in jobs:
            j, pid, mode, r = _part_job(jb)
            got[(j, pid, mode)] = r
    out = []
    for c in cases:
        m = meta[c["j"]]
        pr = [got[(c["j"], pid, "main")] for pid in m["order"]]
        conv = None
        if m.get("conv"):
            conv = {pid: _conv(got[(c["j"], pid, "main")], got[(c["j"], pid, "conv")]) for pid in m["order"]
                    if got.get((c["j"], pid, "conv"), {}).get("status") == "ok" and got[(c["j"], pid, "main")].get("status") == "ok"}
        out.append(summarize_variant(c, pr, m["pending"], conv, round(sum(p.get("t_total_s") or 0 for p in pr), 1)))
    return out, round(time.time() - t0, 1)


# =========================================================================== случаи (вход задач)
def load_cases(path):
    d = json.load(open(path, encoding="utf-8"))
    return d, d["cases"]


def prepare_from_l1(l1_dirs, kv, M_, shards, tag="", feasible_only=True, score_min=None):
    if FUNNEL not in sys.path:
        sys.path.insert(0, FUNNEL)
    import merge_l1 as ML1                       # noqa: E402
    found = ML1.find_parts(l1_dirs, kv, M_, shards, tag)
    cols, rows, det, ref, checks, calls = ML1.merge_parts(found, kv, M_, shards)
    recs = [dict(zip(cols, r)) for r in rows]
    sel = []
    for r in recs:
        d = det.get(r["j"])
        if not d or (feasible_only and not r.get("feasible")):
            continue
        if score_min is not None and (r.get("score") is None or r["score"] < score_min):
            continue
        v = dict(d, j=r["j"], idx=r["idx"], kv=kv, arch=r["arch"], niche=r.get("niche"), score=r.get("score"), feasible=r.get("feasible"),
                 n_rpm=r["n_rpm"], G=r["G"])
        sel.append(P.compact(v))
    sel.sort(key=lambda v: (-(v.get("score") or 0.0), v["j"]))
    src = dict(kind="l1_parts", dirs=l1_dirs, kv=kv, M=M_, shards=shards, input_sha256=ref.get("input_sha256"), code_sha_l1=ref.get("code_sha"),
               merge_checks={k: v for k, v in checks.items() if k in ("files", "missing", "duplicate_conflict", "indices_ok", "n_rows")},
               n_records=len(recs), n_feasible=sum(1 for r in recs if r.get("feasible")), n_with_details=len(det))
    return sel, src


def prepare_from_top50(path):
    d = json.load(open(path, encoding="utf-8"))
    sel = [P.compact(dict(v, kv=v.get("kv") or d.get("kv"), feasible=1, score=v.get("score"), G=v["G"], n_rpm=v["n_rpm"]))
           for v in d["selected"]]
    return sel, dict(kind="top50", file=os.path.relpath(path, GTU), sha256=_sha(path), kv=d.get("kv"), M=d.get("M"))


def write_cases(sel, src, out):
    os.makedirs(os.path.dirname(os.path.abspath(out)), exist_ok=True)
    doc = dict(title="массовая 3D-прочность: случаи (варианты уровня 1)", date=DATE, source=src, n=len(sel), code_sha=code_sha(),
               order="по оценке уровня 1 (убывание), затем j; часть варианта r = r % S", cases=sel)
    with open(out, "w", encoding="utf-8") as fh:
        json.dump(_r(doc, 8), fh, ensure_ascii=False, separators=(",", ":"))
    return out


# =========================================================================== части и слияние
def part_name(kv, n, k, S, tag=""):
    return f"{PART_PREFIX}{tag}_{kv}_n{n}_s{k:03d}of{S:03d}.json.gz"


def run_shards(cases_path, shards, shard_list, out, opts, workers, tag=""):
    doc, cases = load_cases(cases_path)
    kv = doc["source"].get("kv") or (cases[0]["kv"] if cases else "KV")
    os.makedirs(out, exist_ok=True)
    root = tempfile.mkdtemp(prefix="fea3d_", dir=os.environ.get("FEA3D_TMP") or None)
    written = []
    try:
        mine = [(k, [c for r, c in enumerate(cases) if r % shards == k]) for k in shard_list]
        flat = [c for _, cs in mine for c in cs]
        res, wall = run_variants(flat, root, opts, workers, opts.get("conv_every") or 0)
        by_j = {r["j"]: r for r in res}
        for k, cs in mine:
            rows = [by_j[c["j"]] for c in cs]
            meta = dict(level=LEVEL_TAG, kv=kv, n=len(cases), shards=shards, shard=k, cases_sha256=_sha(cases_path),
                        cases_file=os.path.basename(cases_path), code_sha=code_sha(), opts=opts, ccx=C.ccx_version(),
                        j=[c["j"] for c in cs], workers=workers, wall_s=wall, host_nproc=os.cpu_count())
            p = os.path.join(out, part_name(kv, len(cases), k, shards, tag))
            with gzip.open(p, "wt", encoding="utf-8") as fh:
                json.dump(dict(meta, results=rows), fh, ensure_ascii=False)
            written.append(p)
    finally:
        shutil.rmtree(root, ignore_errors=True)
    return written


class MergeError(Exception):
    pass


def _content_hash(rows):
    def strip(o):
        if isinstance(o, dict):
            return {k: strip(v) for k, v in o.items() if k not in ("t_s", "t_total_s", "t_geometry_s", "t_mesh_s", "seconds", "cpu_s",
                                                                   "cpu_ccx_s", "ccx_static", "ccx_freq", "trace")}
        if isinstance(o, list):
            return [strip(v) for v in o]
        return o
    return hashlib.sha256(json.dumps(strip(rows), ensure_ascii=False, sort_keys=True).encode()).hexdigest()


def merge(cases_path, dirs, tag=""):
    """слияние частей: у всех частей тот же файл случаев (sha256), код, параметры; у части k ровно варианты r % S == k;
    ни один вариант не потерян и не задвоен (повтор части - то же содержание без времени счёта). Нарушение - MergeError"""
    doc, cases = load_cases(cases_path)
    n = len(cases)
    sha = _sha(cases_path)
    files = sorted({p for d in dirs for p in glob.glob(os.path.join(d, "**", f"{PART_PREFIX}{tag}_*_n{n}_s*of*.json.gz"), recursive=True)})
    checks = dict(files=len(files), missing=[], duplicate_same=[], duplicate_conflict=[], meta_mismatch=[], index_errors=[], unreadable=[])
    by, ref = {}, None
    for p in files:
        try:
            with gzip.open(p, "rt", encoding="utf-8") as fh:
                d = json.load(fh)
        except (OSError, EOFError, ValueError) as e:
            checks["unreadable"].append(dict(file=os.path.basename(p), error=type(e).__name__))
            continue
        meta = {k: d.get(k) for k in ("level", "n", "shards", "cases_sha256", "code_sha", "opts")}
        if ref is None:
            ref = meta
        bad = {k: [meta[k], ref[k]] for k in meta if meta[k] != ref[k]}
        if meta["cases_sha256"] != sha:
            bad["cases_sha256"] = [meta["cases_sha256"], sha]
        if bad:
            checks["meta_mismatch"].append(dict(file=os.path.basename(p), diff=bad))
            continue
        k, S = d["shard"], d["shards"]
        exp = [c["j"] for r, c in enumerate(cases) if r % S == k]
        got = [x["j"] for x in d["results"]]
        if got != exp or d.get("j") != exp:
            checks["index_errors"].append(dict(file=os.path.basename(p), shard=k, lost=sorted(set(exp) - set(got))[:20],
                                               foreign=sorted(set(got) - set(exp))[:20]))
            continue
        h = _content_hash(d["results"])
        if k in by:
            (checks["duplicate_same"] if by[k]["hash"] == h else checks["duplicate_conflict"]).append(dict(shard=k, file=os.path.basename(p)))
            continue
        by[k] = dict(rows=d["results"], hash=h, wall_s=d.get("wall_s"), workers=d.get("workers"), ccx=d.get("ccx"))
    S = ref["shards"] if ref else None
    if S is not None:
        checks["missing"] = [k for k in range(S) if k not in by]
    errs = {x: checks[x] for x in ("missing", "duplicate_conflict", "meta_mismatch", "index_errors", "unreadable") if checks[x]}
    if errs or ref is None:
        raise MergeError(json.dumps(errs or {"files": "частей не найдено"}, ensure_ascii=False))
    rows = {x["j"]: x for v in by.values() for x in v["rows"]}
    if sorted(rows) != sorted(c["j"] for c in cases) or len(rows) != n:
        raise MergeError("после слияния набор вариантов не равен файлу случаев")
    ordered = [rows[c["j"]] for c in cases]
    for c, rr in zip(cases, ordered):
        if rr["idx"] != c["idx"]:
            raise MergeError(f"j {c['j']}: idx {rr['idx']} не совпадает со случаем {c['idx']}")
    checks.update(indices_ok=True, n=n, shards=S)
    return ordered, ref, checks, by


def summary(rows):
    from collections import Counter
    vc = Counter(r["verdict"] for r in rows)
    reasons = Counter()
    for r in rows:
        for x in r["reasons_hard"]:
            reasons[re.sub(r"[\d.,]+", "#", x.split(": ", 1)[1] if ": " in x else x)[:80]] += 1
    by_arch = {}
    for r in rows:
        a = by_arch.setdefault(r.get("arch") or "?", Counter())
        a[r["verdict"]] += 1
    parts_t = [p.get("t_total_s") or 0 for r in rows for p in r["parts"]]
    cpu = [p.get("cpu_ccx_s") or 0 for r in rows for p in r["parts"]]
    conv = [v for r in rows if r.get("conv") for v in r["conv"].values()]
    conv_s = {}
    for key in ("burst_margin", "hoop_mean", "vm_hub_avg", "vm_blade_avg", "vm_blade_node", "freq_first3_max_abs"):
        vals = [abs(c[key]) for c in conv if key in c]
        if vals:
            conv_s[key] = dict(n=len(vals), p50=float(np.median(vals)), max=float(max(vals)))
    return dict(n=len(rows), verdicts=dict(vc), hard_reasons=dict(reasons.most_common(30)), by_arch={a: dict(v) for a, v in by_arch.items()},
                part_time_s=dict(n=len(parts_t), mean=float(np.mean(parts_t)) if parts_t else None, max=float(max(parts_t)) if parts_t else None),
                ccx_cpu_s=dict(mean=float(np.mean(cpu)) if cpu else None, total=float(sum(cpu))), mesh_convergence=conv_s)


def verdicts_doc(rows, cases_path, ref):
    """вердикты для возврата к записям уровня 1: по j (и idx для сверки); блок results - в формате funnel_l1_fea_results.json
    (ключ «j<номер>:<деталь>», burst_margin_120) - его читает merge_l1.apply_fea_results"""
    res = []
    for r in rows:
        for p in r["parts"]:
            res.append(dict(key=f"j{r['j']}:{p['id']}", result=dict(burst_margin_120=(p.get("burst") or {}).get("margin_120"),
                                                                  verdict=p.get("verdict"), status=p.get("status"))))
    return dict(title="массовая 3D-прочность: вердикты вариантов уровня 1", date=DATE, cases=os.path.basename(cases_path),
                cases_sha256=_sha(cases_path), code_sha=ref.get("code_sha"), criteria=Q.CRIT, criteria_basis=Q.CRIT_BASIS,
                variants=[dict(j=r["j"], idx=r["idx"], verdict=r["verdict"], reasons_hard=r["reasons_hard"], reasons_soft=r["reasons_soft"],
                               errors=r["errors"], summary=r["summary"]) for r in rows], results=res)


def write_merge(cases_path, dirs, out_dir, sfx="", tag=""):
    rows, ref, checks, by = merge(cases_path, dirs, tag)
    os.makedirs(out_dir, exist_ok=True)
    doc, cases = load_cases(cases_path)
    s = summary(rows)
    full = dict(title="массовая 3D-прочность вариантов уровня 1 (CalculiX, секторы с циклической симметрией)", date=DATE,
                method=__doc__.split("\n")[0], software=dict(ccx=next((v["ccx"] for v in by.values() if v.get("ccx")), None),
                                                            python=sys.version.split()[0], numpy=np.__version__),
                cases=os.path.basename(cases_path), cases_sha256=_sha(cases_path), source=doc.get("source"), code_sha=ref["code_sha"],
                opts=ref["opts"], merge_checks=checks, criteria=Q.CRIT, criteria_basis=Q.CRIT_BASIS, geometry_assumptions=G.A,
                parts_rules=dict(low_eo=P.LOW_EO, op_range=P.OP_RANGE, op_basis=P.OP_BASIS, diffuser_z=P.DIFF_Z_RANGE),
                summary=s, variants=rows,
                note="всё - оценка: геометрия колёс и дисков по записям уровня 1 (одномерные модели), поле температур по радиусу, "
                     "свойства материалов - materials.json и оценки materials3d.py; запас по разрушению k = 1 - верхняя оценка")
    p1 = os.path.join(out_dir, f"fea3d_results{sfx}.json.gz")
    with gzip.open(p1, "wt", encoding="utf-8") as fh:
        json.dump(full, fh, ensure_ascii=False)
    p2 = os.path.join(out_dir, f"fea3d_verdicts{sfx}.json")
    with open(p2, "w", encoding="utf-8") as fh:
        json.dump(verdicts_doc(rows, cases_path, ref), fh, ensure_ascii=False, indent=1)
    p3 = os.path.join(out_dir, f"fea3d_summary{sfx}.json")
    with open(p3, "w", encoding="utf-8") as fh:
        json.dump(dict(cases=os.path.basename(cases_path), merge_checks=checks, summary=s, criteria=Q.CRIT), fh, ensure_ascii=False, indent=1)
    return p1, p2, p3, s


# =========================================================================== возврат к уровню 1 и отбор 50
def apply_verdicts(recs, path, reject=("отказ",)):
    """вердикты 3D до отбора 50: вариант с вердиктом из reject становится недопустимым на уровне 1 с причиной fea3d_<вид> и
    расшифровкой; «ошибка» и «условно» - не отсеиваются (пометка в reject_detail). Возврат: результаты по j и число отсеянных"""
    doc = json.load(open(path, encoding="utf-8"))
    by = {int(v["j"]): v for v in doc["variants"]}
    n_bad = 0
    for r in recs:
        v = by.get(int(r["j"]))
        if not v:
            continue
        if v.get("idx") is not None and int(v["idx"]) != int(r["idx"]):
            raise MergeError(f"j {r['j']}: idx уровня 1 {r['idx']} не совпадает с вердиктом 3D {v['idx']}")
        if v["verdict"] in reject and r.get("feasible"):
            r["feasible"] = 0
            kinds = sorted({_reason_kind(x) for x in v["reasons_hard"]})
            r["reject"] = "|".join(x for x in ((r.get("reject") or ""), *[f"fea3d_{k}" for k in kinds]) if x)
            r["reject_detail"] = "; ".join(x for x in ((r.get("reject_detail") or ""), "CalculiX 3D: " + "; ".join(v["reasons_hard"][:4])) if x)
            n_bad += 1
        elif v["verdict"] in ("условно", "ошибка"):
            r["reject_detail"] = "; ".join(x for x in ((r.get("reject_detail") or ""), f"CalculiX 3D {v['verdict']}: " +
                                                       "; ".join((v["reasons_soft"] or v["errors"])[:3])) if x)
    return {int(v["j"]): dict(verdict=v["verdict"], summary=v["summary"]) for v in doc["variants"]}, n_bad


def _reason_kind(x):
    for key, k in (("разрушающей", "burst"), ("σ/σ0,2", "yield"), ("ползуч", "creep"), ("резонанс", "resonance"), ("керамик", "weibull"),
                   ("корень пера при 120", "blade")):
        if key in x:
            return k
    return "other"


def select(l1_dirs, kv, M_, shards, verdicts, out_dir, sfx="", workers=1, no_scale=True):
    """отбор 50 разных уровня 1 с вердиктами 3D: merge_l1.run_merge, в котором отсев по CalculiX (apply_fea_results) заменён
    на apply_verdicts (код воронки не меняется - подмена функции на время вызова)"""
    if FUNNEL not in sys.path:
        sys.path.insert(0, FUNNEL)
    import merge_l1 as ML1                       # noqa: E402
    found = ML1.find_parts(l1_dirs, kv, M_, shards)
    orig = ML1.apply_fea_results
    ML1.apply_fea_results = lambda recs, path, burst_min=None: apply_verdicts(recs, path)
    try:
        r = ML1.run_merge(found, kv, M_, shards, out_dir, sfx, () if no_scale else ("KV-1", "KV-6"), None, None, workers, 0,
                          fea_results=verdicts)
    finally:
        ML1.apply_fea_results = orig
    return r


# =========================================================================== миниатюра
def probe(out_dir=None, workers=None, level="coarse", n=2):
    """локальная миниатюра: n вариантов funnel_l1_top50.json, части 2, слияние, вердикты"""
    out_dir = out_dir or os.path.join(HERE, "probe")
    os.makedirs(out_dir, exist_ok=True)
    sel, src = prepare_from_top50(os.path.join(FUNNEL, "funnel_l1_top50.json"))
    sel = sel[:n]
    cases = write_cases(sel, src, os.path.join(out_dir, "fea3d_cases_probe.json"))
    opts = dict(DEFAULTS, level=level)
    work = tempfile.mkdtemp(prefix="fea3d_probe_")
    t0 = time.time()
    try:
        run_shards(cases, 2, [0, 1], work, opts, workers or min(4, os.cpu_count() or 1))
        p1, p2, p3, s = write_merge(cases, [work], out_dir, "_probe")
    finally:
        shutil.rmtree(work, ignore_errors=True)
    return dict(files=[p1, p2, p3], summary=s, wall_s=round(time.time() - t0, 1))


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    a = sub.add_parser("prepare")
    a.add_argument("--l1-in", nargs="*")
    a.add_argument("--top50")
    a.add_argument("--kv", default="KV-3")
    a.add_argument("--m", type=int)
    a.add_argument("--shards", type=int)
    a.add_argument("--tag", default="")
    a.add_argument("--all-with-details", action="store_true", help="не только допустимые уровня 1")
    a.add_argument("--score-min", type=float, default=None)
    a.add_argument("--out", required=True)
    b = sub.add_parser("run")
    b.add_argument("--cases", required=True)
    b.add_argument("--shards", type=int, required=True)
    b.add_argument("--shard", required=True, help="номера частей через запятую или all")
    b.add_argument("--workers", type=int, default=int(os.environ.get("NPROC", "1") or 1))
    b.add_argument("--level", default=DEFAULTS["level"], choices=list(G.A["levels"]))
    b.add_argument("--conv-every", type=int, default=0, help="каждый K-й вариант - ещё сетка conv-level (сходимость)")
    b.add_argument("--freq-fracs", default=None, help="доли частоты вращения для частот, через запятую")
    b.add_argument("--timeout", type=int, default=DEFAULTS["timeout_s"])
    b.add_argument("--tag", default="")
    b.add_argument("--out", default=os.environ.get("GTU_OUT", "."))
    c = sub.add_parser("merge")
    c.add_argument("--cases", required=True)
    c.add_argument("--in", dest="inp", nargs="+", required=True)
    c.add_argument("--out-dir", required=True)
    c.add_argument("--suffix", default="")
    c.add_argument("--tag", default="")
    d = sub.add_parser("select")
    d.add_argument("--l1-in", nargs="+", required=True)
    d.add_argument("--kv", default="KV-3")
    d.add_argument("--m", type=int, required=True)
    d.add_argument("--shards", type=int, required=True)
    d.add_argument("--verdicts", required=True)
    d.add_argument("--out-dir", required=True)
    d.add_argument("--suffix", default="")
    d.add_argument("--scale", action="store_true")
    e = sub.add_parser("probe")
    e.add_argument("--out-dir", default=None)
    e.add_argument("--workers", type=int, default=None)
    e.add_argument("--level", default="coarse")
    e.add_argument("--n", type=int, default=2)
    x = ap.parse_args(argv)
    if x.cmd == "prepare":
        if x.top50:
            sel, src = prepare_from_top50(x.top50)
        else:
            sel, src = prepare_from_l1(x.l1_in, x.kv, x.m, x.shards, x.tag, not x.all_with_details, x.score_min)
        print(write_cases(sel, src, x.out), len(sel))
    elif x.cmd == "run":
        doc, cases = load_cases(x.cases)
        shard_list = list(range(x.shards)) if x.shard == "all" else [int(s) for s in x.shard.split(",") if s != ""]
        opts = dict(DEFAULTS, level=x.level, conv_every=x.conv_every, timeout_s=x.timeout)
        if x.freq_fracs:
            opts["freq_fracs"] = [float(s) for s in x.freq_fracs.split(",")]
        for p in run_shards(x.cases, x.shards, shard_list, x.out, opts, x.workers, x.tag):
            print(p)
    elif x.cmd == "merge":
        try:
            p1, p2, p3, s = write_merge(x.cases, x.inp, x.out_dir, x.suffix, x.tag)
        except MergeError as e:
            print("слияние отклонено:", e, file=sys.stderr)
            return 2
        print(p1, p2, p3)
        print(json.dumps(s["verdicts"], ensure_ascii=False))
    elif x.cmd == "select":
        r = select(x.l1_in, x.kv, x.m, x.shards, x.verdicts, x.out_dir, x.suffix, no_scale=not x.scale)
        print(json.dumps(r["res"]["selection"]["n"] if isinstance(r, dict) and "res" in r else "готово", ensure_ascii=False))
    elif x.cmd == "probe":
        r = probe(x.out_dir, x.workers, x.level, x.n)
        print(json.dumps(r, ensure_ascii=False, indent=1)[:4000])
    return 0


if __name__ == "__main__":
    sys.exit(main())
