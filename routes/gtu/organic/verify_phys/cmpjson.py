import json, sys, math
a = json.load(open(sys.argv[1])); b = json.load(open(sys.argv[2]))
skip = set(sys.argv[3].split(",")) if len(sys.argv) > 3 else set()
diffs = []
def walk(x, y, p):
    if p.split("/")[-1] in skip: return
    if isinstance(x, dict) and isinstance(y, dict):
        for k in set(x) | set(y):
            if k not in x or k not in y: diffs.append((p + "/" + k, "missing", None)); continue
            walk(x[k], y[k], p + "/" + k)
    elif isinstance(x, list) and isinstance(y, list):
        if len(x) != len(y): diffs.append((p, "len", (len(x), len(y)))); return
        for i, (u, v) in enumerate(zip(x, y)): walk(u, v, f"{p}[{i}]")
    elif isinstance(x, (int, float)) and isinstance(y, (int, float)) and not isinstance(x, bool):
        if x != y:
            rel = abs(x - y) / max(abs(x), abs(y), 1e-30)
            diffs.append((p, rel, (x, y)))
    elif x != y:
        diffs.append((p, "val", (str(x)[:60], str(y)[:60])))
walk(a, b, "")
num = [d for d in diffs if isinstance(d[1], float)]
oth = [d for d in diffs if not isinstance(d[1], float)]
print("числовых расхождений:", len(num), "макс. отн.:", max([d[1] for d in num], default=0))
for d in sorted(num, key=lambda d: -d[1])[:12]: print("  ", d)
print("прочих:", len(oth))
for d in oth[:15]: print("  ", d)
