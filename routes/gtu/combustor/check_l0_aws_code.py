"""На каком коде шёл уровень 0 камеры на AWS (прогон 20261115T001745Z-f1da, n = 262 144) - проверка двумя способами.

1) отпечаток кода: code_hash частей (completeness.ref.code_hash итога слияния - слияние сверяет его у всех 48 частей)
   против sha256(comb_space.py + comb_crn.py) каждого коммита истории этих файлов (git show);
2) пересчёт: первые K строк рейтинга AWS заново - кодом рабочей копии и кодом другой копии (--other, например
   `git archive 8321dd4` - код до исправлений скептика 14.11.2026), сравнение NOx, CO, запасов по срыву, потерь, PF,
   стенки, деталей.
Числа пересчёта - уровень облика (как в итоге AWS). AWS не вызывается.

Запуск: python check_l0_aws_code.py --other <папка combustor другой копии> --out l0_aws/comb_l0_aws_code_check_results.json
"""
import argparse, hashlib, json, os, subprocess, sys

HERE = os.path.dirname(os.path.abspath(__file__))
AWS = os.path.join(HERE, "l0_aws", "comb_l0_KV-3_n262144_results.json")
REL = "routes/gtu/combustor"


def history_hashes():
    repo = subprocess.run(["git", "rev-parse", "--show-toplevel"], cwd=HERE, capture_output=True, text=True).stdout.strip()
    log = subprocess.run(["git", "log", "--format=%h|%ad|%s", "--date=iso", "--", f"{REL}/comb_space.py",
                          f"{REL}/comb_crn.py"], cwd=repo, capture_output=True, text=True).stdout.splitlines()
    out = []
    for ln in log:
        h, date, subj = ln.split("|", 2)
        blobs = []
        for f in ("comb_space.py", "comb_crn.py"):
            r = subprocess.run(["git", "show", f"{h}:{REL}/{f}"], cwd=repo, capture_output=True)
            if r.returncode != 0:
                blobs = None
                break
            blobs.append(r.stdout)
        if blobs is None:
            continue
        out.append(dict(commit=h, date=date, subject=subj, code_hash=hashlib.sha256(b"".join(blobs)).hexdigest()[:16]))
    return out


def recompute(code_dir, idxs, n):
    """пересчёт строк в отдельном процессе с кодом из code_dir (модули v1 кэшируются - один процесс на копию)"""
    src = r"""
import json, sys, os
code_dir, n, idxs = sys.argv[1], int(sys.argv[2]), json.loads(sys.argv[3])
sys.path.insert(0, code_dir); os.chdir(code_dir)
import comb_space as CS, comb_crn as CR
out = dict(code_hash=CS.code_hash(), bc_key=CR.boundary(write=False)["key"], rows=[])
for i in idxs:
    r = CR.evaluate(i, n)
    g = r.get("gas") or []
    out["rows"].append(dict(idx=i, status=r["status"], NOx=[x["NOx"] for x in g], CO=[x["CO"] for x in g],
                            margin_lbo=[x["margin_lbo"] for x in g], dp_total=(r.get("dp") or {}).get("dp_total"),
                            PF=(r.get("pf") or {}).get("PF"), parts=(r.get("module") or {}).get("parts"),
                            T_wall_K=(r.get("feed") or {}).get("T_metal_max_K")))
print(json.dumps(out))
"""
    r = subprocess.run([sys.executable, "-c", src, code_dir, str(n), json.dumps(idxs)], capture_output=True, text=True,
                       timeout=3600)
    if r.returncode != 0:
        raise RuntimeError(r.stderr[-2000:])
    return json.loads(r.stdout.strip().splitlines()[-1])


def rel_diff(a, b):
    if isinstance(a, list):
        return max((rel_diff(x, y) for x, y in zip(a, b)), default=0.0)
    if a is None or b is None:
        return None if a != b else 0.0
    return abs(a - b) / max(abs(a), abs(b), 1e-9)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--other", help="папка combustor другой копии кода (например git archive 8321dd4)")
    ap.add_argument("--k", type=int, default=5)
    ap.add_argument("--out", default=os.path.join(HERE, "l0_aws", "comb_l0_aws_code_check_results.json"))
    a = ap.parse_args()
    aws = json.load(open(AWS, encoding="utf-8"))
    ref = aws["completeness"]["ref"]
    hist = history_hashes()
    match = [h for h in hist if h["code_hash"] == ref["code_hash"]]
    top = aws["top"][:a.k]
    idxs = [r["idx"] for r in top]
    keys = dict(NOx="NOx", CO="CO", margin_lbo="margin_lbo", dp_total="dp_total", PF="PF", parts="parts", T_wall_K="T_wall_K")
    runs = {"working_copy": recompute(HERE, idxs, ref["n"])}
    if a.other:
        runs["other"] = recompute(os.path.abspath(a.other), idxs, ref["n"])
    cmp = {}
    for name, rr in runs.items():
        cmp[name] = []
        for r, t in zip(rr["rows"], top):
            cmp[name].append(dict(idx=r["idx"], max_rel_diff={k: rel_diff(r[k], t[v]) for k, v in keys.items()}, recomputed=r,
                                  aws={k: t[v] for k, v in keys.items()}))
    res = dict(
        title="Камера KV-3, уровень 0: на каком коде шёл прогон AWS 20261115T001745Z-f1da (n = 262 144)",
        method=("1) отпечаток кода частей (completeness.ref.code_hash) против sha256(comb_space.py + comb_crn.py) каждого "
                "коммита истории; 2) пересчёт первых K строк рейтинга AWS кодом рабочей копии и кодом другой копии"),
        aws_ref=ref, history=hist, aws_code_matches_commits=match,
        recompute=dict(k=a.k, idx=idxs, runs={k: dict(code_hash=v["code_hash"], bc_key=v["bc_key"]) for k, v in runs.items()},
                       compare=cmp),
        verdict=None)
    w = cmp["working_copy"]
    same = runs["working_copy"]["code_hash"] == ref["code_hash"] and \
        all(all((x is not None and x < 1e-4) for x in r["max_rel_diff"].values()) for r in w)
    res["verdict"] = ("код прогона AWS = код рабочей копии (отпечаток совпал, строки воспроизведены до 1e-4)" if same else
                      "код прогона AWS отличается от рабочей копии - разбор")
    json.dump(res, open(a.out, "w", encoding="utf-8"), ensure_ascii=False, indent=1)
    print(res["verdict"], "| совпадения отпечатка с коммитами:", [m["commit"] for m in match])
    for name in cmp:
        print(name, runs[name]["code_hash"], [round(max(v for v in r["max_rel_diff"].values() if v is not None), 4)
                                             for r in cmp[name]])


if __name__ == "__main__":
    main()
