"""пересчитать только reformer.catalyst_estimate в готовом сжатом cycles_results.json исправленной функцией"""
import json, sys, os
sys.path.insert(0, sys.argv[1])
import cycles as C
p = sys.argv[2]
raw = open(p, encoding="utf-8").read()
out = json.loads(raw)
old = out["reformer"]["catalyst_estimate"]
new = C._r5(C.catalyst_estimate(out))
print("было:", old); print("стало:", new)
out["reformer"]["catalyst_estimate"] = new
if len(sys.argv) > 3 and sys.argv[3] == "--write":
    json.dump(out, open(p, "w", encoding="utf-8"), ensure_ascii=False, separators=(",", ":"), default=C._jdef)
    print("записано", p)
