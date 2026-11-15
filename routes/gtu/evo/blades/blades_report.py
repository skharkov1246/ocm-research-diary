#!/usr/bin/env python3
"""Сводка семейства «лопатки и колёса» → blades_mini_results.json (в репозитории): мини-прогон (evo_blades.py mini),
проверка элит CalculiX (ccx_elites.py), калибровка быстрой оценки (disk1d_calibration.json), план AWS (runner.py plan),
статусы приёмов (blade_eval.TECH), входы и версии ПО.

  python3 blades_report.py --mini DIR/mini_results.json --ccx DIR/ccx_check --render DIR/render/render_report.json \\
      --plan-usd 36.25 --out blades_mini_results.json
"""
import argparse
import glob
import hashlib
import json
import os
import platform
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)


def sha(p):
    return hashlib.sha256(open(p, "rb").read()).hexdigest()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mini", required=True)
    ap.add_argument("--ccx", required=True)
    ap.add_argument("--render", default=None)
    ap.add_argument("--plan-usd", type=float, default=None)
    ap.add_argument("--out", default=os.path.join(HERE, "blades_mini_results.json"))
    a = ap.parse_args()
    import numpy
    import scipy
    import blade_eval as BE
    import core_local as C
    import disk1d as D1
    mini = json.load(open(a.mini, encoding="utf-8"))
    cal = json.load(open(D1.CAL_JSON, encoding="utf-8"))
    ccx = {}
    for p in sorted(glob.glob(os.path.join(a.ccx, "ccx_elites_*.json"))):
        d = json.load(open(p, encoding="utf-8"))
        ccx[d["summary"]["comp"]] = d
    comps = {}
    for comp, v in mini["comps"].items():
        el = []
        cc = ccx.get(comp)
        for e in v["elites"]:
            m = e["meta"]
            rec = dict(id=e["id"], f=e["f"], feasible=e["feasible"], descriptors=dict(zip(e["desc_names"], e["d"])), genes=e["genes"],
                       effect=e["effect"],
                       numbers={k: m.get(k) for k in ("eta", "eta_level1", "eta_hyp", "SM", "SM_hyp", "mass_kg", "J_kgm2", "u_disk", "burst",
                                                     "sig_bend_MPa", "D2", "U2", "beta2b", "Zeff", "parts", "parts_inserted_equiv", "route",
                                                     "hollow", "Dm", "Dtip", "length_m") if k in m},
                       corrections=m.get("corrections"), corrections_hyp=m.get("corrections_hyp"))
            if comp == "turb":
                rec["stages"] = m.get("stages")
            if cc:
                rows = [r for r in cc["rows"] if r["tag"].startswith(e["id"])]
                if comp == "turb":
                    rec["ccx"] = rows[0]["stages"] if rows else None
                else:
                    rec["ccx"] = [dict(speed_pct=int(r["tag"].rsplit("_", 1)[1]), ccx_yield_ratio=r["ccx"]["yield_ratio_max"],
                                       fast_yield_ratio=r["fast_cal"], ratio=r["ratio_ccx_to_fast_cal"], burst_ccx=r["ccx"]["burst_margin"],
                                       gate=r["ccx"]["gate"]) for r in rows]
                    y = {x["speed_pct"]: x["ccx_yield_ratio"] for x in rec["ccx"]}
                    if 100 in y and 120 in y:
                        rec["u_disk_ccx"] = max(1.5 * y[100], y[120]) / BE.K_ROUTE
                        rec["ccx_pass"] = rec["u_disk_ccx"] <= 1.0
            el.append(rec)
        comps[comp] = dict({k: v[k] for k in ("merge_report", "n_evals_me", "n_evals_nsga", "n_evals_cma", "cells_total", "cells_filled",
                                             "cells_feasible", "coverage_feasible", "f_best", "f_median", "reference", "nsga_front", "cma",
                                             "seconds", "timing") if k in v},
                           elites=el, ccx_summary=cc["summary"] if cc else None)
    out = dict(
        title="Семейство «лопатки и колёса» KV-3: мини-прогон поиска формы (MAP-Elites → NSGA-II → CMA-ES), проверка элит CalculiX",
        date=mini["date"], status=("оценка уровня облика: КПД и запас помпажа — модель уровня 1 воронки + поправки приёмов (статусы в tech);"
                                   " прочность — быстрая оценка, калиброванная по CalculiX, элиты проверены CalculiX (осесимметрия);"
                                   " AWS не вызывался (только runner.py plan)"),
        method=("острова MAP-Elites (архив CVT, операторы iso+line и гауссова мутация), слияние островов с контролем полноты,"
                " NSGA-II по (−η, −запас помпажа, J, m, загрузка диска) для колёс и (−η, J, m) для турбины, CMA-ES — доводка элит в"
                " своей ячейке; пригодность и политика поправок — CONFIG evo_blades.py"),
        software=dict(python=platform.python_version(), numpy=numpy.__version__, scipy=scipy.__version__, ccx=D1.DK.ccx_version()),
        inputs=[dict(file="routes/gtu/frame3/frame3_point_results.json", sha256=sha(BE.FRAME3), field=BE.point()["src"]["fields"]),
                dict(file="routes/gtu/frame3/frame3_fea_results.json", sha256=sha(BE.FRAME3_FEA),
                     field="parts.comp1/comp2 (профили, поля T, ccx — калибровка), turb1..3 (профили дисков, поля T, геометрия)"),
                dict(file="routes/gtu/organic/structure/profile_opt_comp2.json", field="base, check_160x24 (калибровка пика)"),
                dict(file="routes/gtu/funnel/meanline_l1.py", field="cc_stage, cc_off, ax_turbine (модель уровня 1)"),
                dict(file="routes/gtu/funnel/offdesign_l1.py", field="DINC_STALL, TIP_INC_STALL, ALPHA2_STALL"),
                dict(file="routes/gtu/funnel/funnel_l1.py", field="M1T_MAX, DH_MIN, D1D2_MAX, ALPHA2_MAX, BETA2_MAX, M2T_MAX, MW3_MAX, R_MIN"),
                dict(file="routes/gtu/platform/fea/materials.json", field="ti64, in713c"),
                dict(file="routes/gtu/platform/fea/disk.py", field="run_case, load_material, LMP (CalculiX)")],
        config=dict(config_hash=mini.get("config_hash"), code_hash=mini.get("code_hash"),
                    fitness=__import__("evo_blades").CONFIG["fitness"]),
        tech=BE.TECH,
        calibration=dict(k_peak=cal["k_peak"], k_peak_range=[cal["k_min"], cal["k_max"]], table=cal.get("table"),
                         table_check={k: v for k, v in cal.get("table_check", {}).items() if k != "cases"},
                         domain=cal.get("domain", {}).get("rule")),
        comps=comps,
        render=json.load(open(a.render, encoding="utf-8")) if a.render and os.path.exists(a.render) else None,
        aws=dict(jobs="routes/gtu/evo/blades/jobs_evo_blades.json", plan_usd_at_1usd_h=a.plan_usd,
                 note="runner.py plan --assume-price 1.0: код 0; AWS не вызывался"),
        seconds_mini=mini.get("seconds"))
    with open(a.out, "w", encoding="utf-8") as f:
        json.dump(C._canon(out), f, ensure_ascii=False, indent=1)
    print(a.out, os.path.getsize(a.out) // 1024, "КБ")


if __name__ == "__main__":
    main()
