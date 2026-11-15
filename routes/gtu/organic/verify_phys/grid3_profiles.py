"""Независимость от сетки проверенных профилей (disk.py / CalculiX): 80x12, 160x24, 320x48.
Колесо ступени 2 (исходное и profile_eq) и диск турбины 1 (исходный, eq, lim)."""
import json, os, sys, time
import numpy as np
S = sys.argv[1]  # папка routes/gtu/organic/structure (в раскладке репозитория); argv[2] — выход .json
# запуск: ORGANIC_OUT=<рабочая папка> python3 grid3_profiles.py <structure> <выход.json>  (~4 мин на 1 ядре)
sys.path.insert(0, S)
import disk_topo as DT
import structure_post as SP
import profile_opt as PO
out = {}
c2 = json.load(open(os.path.join(S, "profile_opt_comp2.json")))
d = DT.frame3(); pr = d["parts"]["comp2"]["variants"]["solid"]["profile"]
r0 = np.asarray(pr["r"]); zb0 = np.asarray(pr["zc"]) - np.asarray(pr["t"]) / 2
x0 = np.interp(PO.RC_KNOTS, r0, zb0)
xb = np.asarray(c2["best"]["x"])
grids = [(80, 12), (160, 24), (320, 48)]
for tag, x in (("comp2_base", x0), ("comp2_eq", xb)):
    for g in grids:
        t0 = time.time()
        a = PO.comp2_case(f"g3_{tag}_{g[0]}", x, DT.RPM0, mesh=g)
        b = PO.comp2_case(f"g3_{tag}_{g[0]}", x, 1.2 * DT.RPM0, mesh=g)
        out[f"{tag}_{g[0]}x{g[1]}"] = dict(y100=a["yield_ratio"], y120=b["yield_ratio"], at100=a["yield_at"], at120=b["yield_at"],
                                          vm100=a["vm_max_MPa"], burst=a["burst_margin"], mass=a["mass_kg"], s=time.time() - t0)
        print(tag, g, json.dumps(out[f"{tag}_{g[0]}x{g[1]}"], default=float), flush=True)
g_ = PO.geom()
base = dict(r=g_["base"]["r"], t=g_["base"]["t"], zc=[0.0] * len(g_["base"]["r"]))
for f, tag in (("profile_opt_eq.json", "turb_eq"), ("profile_opt_lim.json", "turb_lim"), (None, "turb_base")):
    prof = base if f is None else PO.build(np.asarray(json.load(open(os.path.join(S, f)))["best"]["x"]), g_)
    for g in grids:
        t0 = time.time()
        a = SP.disk_case(f"g3_{tag}_{g[0]}", prof, DT.RPM0, mesh=g)
        b = SP.disk_case(f"g3_{tag}_{g[0]}", prof, 1.2 * DT.RPM0, mesh=g)
        out[f"{tag}_{g[0]}x{g[1]}"] = dict(life100=a["creep_local_min_h"], uts120=b["uts_ratio"], vm100=a["vm_max_MPa"],
                                          burst=a["burst_margin"], mass=a["mass_kg"], creep_at=a["creep_at"], s=time.time() - t0)
        print(tag, g, json.dumps(out[f"{tag}_{g[0]}x{g[1]}"], default=float), flush=True)
json.dump(out, open(sys.argv[2] if len(sys.argv) > 2 else os.path.join(os.path.dirname(os.path.abspath(__file__)), "grid3_profiles.json"), "w"), indent=1, default=float)
