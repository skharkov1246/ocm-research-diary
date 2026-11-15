"""Вторая сетка для 2D-оптимизации коллектора (thermal/topo.py): случай comb6 при nx 80 (как в репозитории) и nx 120.
Радиус фильтра в topo.py задан в ячейках (1,5), поэтому на мелкой сетке он меньше в метрах — это не чистая проверка сетки.
Аргументы: <topo_results nx80> <topo_results nx120> <выход.json>"""
import json, sys
a = json.load(open(sys.argv[1])); b = json.load(open(sys.argv[2]))
keys = ("nonuniformity", "balanced_ratio_uniform_vs_conventional", "balanced_ratio_bp_vs_conventional",
        "balanced_ratio_warm_vs_conventional", "dissipation_ratio_uniform_vs_comb")
out = dict(note=__doc__.splitlines()[1], parameters_nx80=a["parameters"], parameters_nx120=b["parameters"], comb6={})
for k in keys:
    out["comb6"][k] = dict(nx80=a["cases"]["comb6"].get(k), nx120=b["cases"]["comb6"].get(k))
json.dump(out, open(sys.argv[3], "w"), ensure_ascii=False, indent=1, default=float)
print(json.dumps(out["comb6"], ensure_ascii=False, default=float))
