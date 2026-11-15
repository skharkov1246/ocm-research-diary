#!/usr/bin/env python3
"""Сводка семейства «корпус-оболочка и рама» → evo_shell_results.json (числа — из файлов прогона, ничего не вписывается руками).

  python3 shell_report.py --run DIR [--ccx f1.json …] [--patch patch_scan.json] [--sdf sdf.json] [--gallery-o g.json] [--gallery-k g.json]
                          [--plan plan.txt] --out evo_shell_results.json
DIR — папка мини-прогона (merged_o/k.json, nsga_o/k.json, polish_o/k.json).
"""
import argparse
import json
import os
import platform
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.dirname(HERE))
import shell_model as SM  # noqa: E402

KEYS = ("mass_kg", "n_eff", "parts", "fasteners", "pipe_joints", "service", "q", "f_min_Hz", "u_wall", "u_core", "u_strut_f",
        "u_strut_r", "u_cage", "lat_frac", "D_env", "printer", "pieces", "overhang", "pore_m", "curv_excess")


def _o(o):
    return {k: (round(float(o[k]), 4) if isinstance(o.get(k), (int, float, np.floating)) and not isinstance(o.get(k), bool) else o.get(k))
            for k in KEYS if k in o}


def variant_block(run, tag, variant):
    m = json.load(open(os.path.join(run, f"merged_{tag}.json"), encoding="utf-8"))
    el = m["archive"]["elites"]
    fe = sorted([e for e in el if e["objs"].get("CV", 1) == 0], key=lambda e: (-e["fitness"], e["uid"]))
    ref = SM.reference(variant)
    per = {}
    for j, name in enumerate(SM.DESC_NAMES):
        used = sorted({e["cell"][j] for e in fe})
        per[name] = dict(bins=len(SM.DESC_LABELS[j]), feasible_bins=[SM.DESC_LABELS[j][b] for b in used])
    best = fe[0] if fe else None
    out = dict(variant=variant, merge=dict((k, m["report"][k]) for k in ("status", "expected", "present", "missing", "evals",
                                                                         "merged_filled", "merged_coverage")),
               cells_total=int(np.prod([len(x) for x in SM.DESC_LABELS])), feasible_elites=len(fe), descriptors_feasible=per,
               reference=ref["rows"],
               best=dict(uid=best["uid"], fitness=best["fitness"], cell=[SM.DESC_LABELS[j][c] for j, c in enumerate(best["cell"])],
                         objs=_o(best["objs"]), params={k: SM.decode(best["genome"])[k] for k in
                                                        ("route", "rib_type", "sym_n", "lat_type", "f_n", "f_lean", "f_bow", "r_n", "r_lean",
                                                         "r_bow", "hf_n", "hf_shape", "hf_t", "hf_b", "ring_h", "emb", "cast_mat")}) if best else None,
               top10=[dict(uid=e["uid"], fitness=round(e["fitness"], 4), cell=[SM.DESC_LABELS[j][c] for j, c in enumerate(e["cell"])],
                           objs=_o(e["objs"])) for e in fe[:10]])
    if best:
        b = best["objs"]
        r0 = ref["rows"][0]
        out["vs_reference"] = dict(mass_ratio=b["mass_kg"] / r0["mass_kg"], parts_delta=b["parts"] - r0["parts"],
                                   fasteners_delta=b["fasteners"] - r0["fasteners"], n_eff_ratio=b["n_eff"] / r0["n_eff"],
                                   service_delta=b["service"] - r0["service"])
    fn = os.path.join(run, f"nsga_{tag}.json")
    if os.path.exists(fn):
        n = json.load(open(fn, encoding="utf-8"))
        fr = n["front"]
        F = np.array([r["F"] for r in fr]) if fr else np.zeros((0, 5))
        out["nsga"] = dict(front=len(fr), evals=n["evals"], hv_mc=n["hv_mc"], objectives=n["objectives"],
                           ranges={o: [float(F[:, j].min()), float(F[:, j].max())] for j, o in enumerate(n["objectives"])} if len(F) else {},
                           extremes={nm: dict(objs=_o(min(fr, key=lambda r: r["F"][j])["objs"]),
                                              novelty_oem=min(fr, key=lambda r: r["F"][j])["novelty_oem"])
                                     for j, nm in ((0, "наименьшая масса"), (1, "наименьший N_eff"), (4, "наибольшая новизна"))} if fr else {})
    fn = os.path.join(run, f"polish_{tag}.json")
    if os.path.exists(fn):
        pl = json.load(open(fn, encoding="utf-8"))
        out["polish"] = [dict(uid=r["uid"], before=r["fitness_before"], after=r["fitness_after"], same_cell=r["same_cell"], evals=r["evals"],
                              mass_before=r["objs_before"].get("mass_kg"), mass_after=r["objs_after"].get("mass_kg")) for r in pl["rows"]]
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--run", required=True); ap.add_argument("--ccx", nargs="*", default=[]); ap.add_argument("--patch")
    ap.add_argument("--sdf"); ap.add_argument("--gallery-o"); ap.add_argument("--gallery-k"); ap.add_argument("--plan")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    res = dict(
        title="KV-3, семейство «корпус-оболочка и рама»: MAP-Elites + NSGA-II + CMA-ES на гладких формах, мини-прогон",
        date="2026-11-15",
        status="черновик: мини-прогон локально (2 варианта × 8 островов × 4000 оценок); все числа — дешёвые модели уровня облика "
               "(оценка), сверка с CalculiX грубой сетки — ниже; AWS не вызывался (только runner.py plan)",
        software=dict(python=platform.python_version(), numpy=np.__version__, scipy=__import__("scipy").__version__),
        inputs=SM.INP["sources"], zones=[{k: (round(v, 4) if isinstance(v, float) else v) for k, v in z.items()} for z in SM.INP["zones"]],
        masses=dict(rotor_kg=SM.INP["m_rotor"], hot_kg=SM.INP["m_hot"], dr_hot_mm=SM.INP["dr_hot"] * 1e3,
                    basis="ротор и горячая часть — формулы funnel_space.derive (оценка); рассогласование горячей части — print_review"),
        criteria=dict(f_margin=SM.F_MARGIN, f_op_Hz=SM.INP["f_op"], sf_yield=SM.SF_YIELD, k_route=SM.K_ROUTE, w_fit=SM.W_FIT,
                      note="все — принятые для отбора (оценка), не нормы"),
        descriptors=dict(names=SM.DESC_NAMES, edges=SM.DESC_EDGES, labels=SM.DESC_LABELS),
        genes=[dict(name=g["name"], lo=g["lo"], hi=g["hi"], kind=g["kind"], opts=g["opts"]) for g in SM.G],
        variants={})
    for tag, v in (("o", "оболочка"), ("k", "картридж")):
        if os.path.exists(os.path.join(a.run, f"merged_{tag}.json")):
            res["variants"][v] = variant_block(a.run, tag, v)
    cc = []
    for fn in a.ccx:
        d = json.load(open(fn, encoding="utf-8"))
        c = d.get("ccx", d)
        cc.append(dict(case=d.get("uid") or c.get("tag"), variant=d.get("variant"), mesh=c.get("mesh"), seconds=c.get("seconds"), solids=c.get("solids"),
                       beams=c.get("beams"), freq_Hz=c.get("freq_Hz"), modes=c.get("modes"), shell_vm_max_MPa=c.get("shell_vm_max_MPa"),
                       shell_vm_p999_MPa=c.get("shell_vm_p999_MPa"), zones=(d.get("compare") or {}).get("zones"),
                       cheap_f_beam_Hz=(d.get("compare") or {}).get("cheap_f_beam_Hz"), donnell_Hz=d.get("donnell_Hz"), ratio=d.get("ratio")))
    res["ccx"] = cc
    if a.patch:
        rows = json.load(open(a.patch, encoding="utf-8"))
        agree_new = agree_old = 0
        for r in rows:
            pr = SM.lattice_pore(r["kind"], r["cell_m"], r["rho"], r["t_core_m"])
            v3, v5 = r["blocked_share"]["3mm"] <= 0.6, r["blocked_share"]["5mm"] <= 0.6
            agree_new += ((pr >= SM.PORE_MIN) == v3) + ((pr >= SM.PORE_LONG) == v5)
            agree_old += ((r["pore_rule_m"] >= SM.PORE_MIN) == v3) + ((r["pore_rule_m"] >= SM.PORE_LONG) == v5)
        res["patch_scan"] = dict(cases=len(rows), decisions=2 * len(rows), agree_old_rule=agree_old, agree_new_rule=agree_new,
                                 rows=[dict(kind=r["kind"], cell_mm=r["cell_m"] * 1e3, t_core_mm=r["t_core_m"] * 1e3, rho=r["rho"],
                                            blocked_3mm=r["blocked_share"]["3mm"], blocked_5mm=r["blocked_share"]["5mm"]) for r in rows])
    if a.sdf and os.path.exists(a.sdf):
        res["sdf_checks"] = json.load(open(a.sdf, encoding="utf-8"))
    for tag, fn in (("o", a.gallery_o), ("k", a.gallery_k)):
        if fn and os.path.exists(fn):
            g = json.load(open(fn, encoding="utf-8"))
            res.setdefault("gallery", {})[tag] = dict(contact=g["contact"], sheets=g["sheets"], px_per_m=g["px_per_m"],
                                                      rows=[dict(n=r["n"], uid=r["uid"], cell=[SM.DESC_LABELS[j][c] for j, c in enumerate(r["cell"])],
                                                                 fitness=round(r["fitness"], 4), novelty_oem=round(r["novelty_oem"], 3),
                                                                 mass_kg=round(r["objs"]["mass_kg"], 1), n_eff=r["objs"]["n_eff"],
                                                                 params=r["params"]) for r in g["rows"]])
    if a.plan and os.path.exists(a.plan):
        res["aws_plan"] = open(a.plan, encoding="utf-8").read().splitlines()[-40:]
    with open(a.out, "w", encoding="utf-8") as f:
        json.dump(res, f, ensure_ascii=False, indent=1, default=lambda o: float(o) if isinstance(o, np.floating) else
                  (o.tolist() if isinstance(o, np.ndarray) else str(o)))
    print("->", a.out)


if __name__ == "__main__":
    main()
