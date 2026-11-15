"""Сходимость оптимизации топологии disk_topo.py по журналам итераций (log_<прогон>.txt в $ORGANIC_OUT) и сравнение
с прогоном на вдвое крупной сетке (disk_topo.py --only … --h-t 2e-3 --h-c 3e-3 --results <coarse.json>).
Аргументы: <папка журналов основного прогона> <disk_topo_results.json основной> <coarse.json> <папка журналов coarse> <выход.json>"""
import glob, json, os, sys
logs, fine_json, coarse_json, coarse_logs, out_path = sys.argv[1:6]

def stats(path, n=15):
    L = [json.loads(l) for l in open(path) if l.startswith('{"it"')]
    if not L:
        return None
    last = L[-n:]
    ms = [x["mass"] for x in last]; um = [max(x["umax"]) for x in last]
    return dict(last_it=L[-1]["it"], beta=L[-1].get("beta"), change_last=L[-1]["change"], mass_min=min(ms), mass_max=max(ms),
                dmass_per_it=(last[-1]["mass"] - last[0]["mass"]) / len(last), umax_min=min(um), umax_max=max(um),
                grey_last=L[-1].get("grey"), n_window=len(last))

F = json.load(open(fine_json)); C = json.load(open(coarse_json))
fr = {r["name"]: r for r in F["runs"]}
out = dict(runs={}, coarse_vs_fine={})
for r in F["runs"]:
    p = os.path.join(logs, f"log_{r['name']}.txt")
    out["runs"][r["name"]] = dict(iterations=r["iterations"], constraint_ok=r["constraint_ok"], mass_kg=r["mass_kg"],
                                  R_filter_m=r["R_filter_m"], log_window=stats(p) if os.path.exists(p) else None)
for r in C["runs"]:
    o = fr[r["name"]]
    ub = lambda x: x["repairs"][-1]["u_binary"] if isinstance(x["repairs"], list) and x["repairs"] else None
    p = os.path.join(coarse_logs, f"log_{r['name']}.txt")
    out["coarse_vs_fine"][r["name"]] = dict(
        coarse=dict(h_m=C["assumptions"]["h_turb_m" if r["part"] == "turb1" else "h_comp_m"], R_filter_m=r["R_filter_m"],
                    mass_kg=r["mass_kg"], J_kgm2=r["J_kgm2"], u_binary=ub(r), constraint_ok=r["constraint_ok"],
                    log_window=stats(p) if os.path.exists(p) else None),
        fine=dict(h_m=F["assumptions"]["h_turb_m" if r["part"] == "turb1" else "h_comp_m"], R_filter_m=o["R_filter_m"],
                  mass_kg=o["mass_kg"], J_kgm2=o["J_kgm2"], u_binary=ub(o), constraint_ok=o["constraint_ok"]),
        d_mass_rel=(r["mass_kg"] - o["mass_kg"]) / o["mass_kg"])
out["coarse_baseline"] = {k: dict(mass_kg=v["mass_kg"], u=[c["u_max"] for c in v["cases"]]) for k, v in C["baseline"].items()}
out["fine_baseline"] = {k: dict(mass_kg=v["mass_kg"], u=[c["u_max"] for c in v["cases"]]) for k, v in F["baseline"].items()}
json.dump(out, open(out_path, "w"), ensure_ascii=False, indent=1, default=float)
print(json.dumps(out["coarse_vs_fine"], ensure_ascii=False, default=float)[:900])
