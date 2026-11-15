"""Камера KV-3, уровень 0 v2: чувствительность пробы к закону затухания при малом J (BETA_LOWJ, оценка).

Собирает comb_l0_v2_beta_sens_results.json из трёх итогов слияния (merge_l0_v2.py) одной и той же выборки при
β = 0 (затухание как v1), 1 (номинал) и 2: счётчики, ячейки «схема × пассивный элемент», ближайшие к новой цели.
Сами итоги считаются командами:
  python comb_l0_v2.py --n 300 --shard S --shards 4 --beta B --out DIR_B --workers 4   (S = 0..3)
  python merge_l0_v2.py DIR_B --out merged_bB.json --top 30
Запуск: python beta_sens_v2.py merged_b0.json merged_b1.json merged_b2.json [--out comb_l0_v2_beta_sens_results.json]
AWS не вызывается.
"""
import argparse, json, os, sys

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "comb_l0_v2_beta_sens_results.json")
CELL_KEYS = ("total", "ok", "stable_30_100", "stable_in_domain", "stable_0_100", "lowem", "new_value_min",
             "a_eff_30_min", "share_anchor_30_over_100_max")


def build(paths):
    res = dict(title="Камера KV-3, уровень 0 v2: чувствительность пробы N = 300 к закону затухания при малом J",
               method="та же выборка (зерно 20261115, N = 300), BETA_LOWJ = 0 (затухание как v1), 1 (номинал), 2; "
                      "comb_l0_v2.py --beta, merge_l0_v2.py, сборка - beta_sens_v2.py",
               status="оценка: закон затухания при малом J не сверен, проверка - уровень 1 (RANS)",
               code_hash={}, data_hash={}, bc_key={}, complete={}, beta={}, counts={}, by_cell={}, nearest={})
    for p in paths:
        m = json.load(open(p, encoding="utf-8"))
        ref = m["completeness"]["ref"]
        b = ref["beta"]
        key = f"{b:g}"
        res["beta"][key] = b
        res["code_hash"][key] = ref["code_hash"]
        res["data_hash"][key] = ref.get("data_hash")
        res["bc_key"][key] = ref["bc_key"]
        res["complete"][key] = m["complete"]
        res["counts"][key] = m["counts"]
        for cell, v in m["by_scheme_passive"].items():
            res["by_cell"].setdefault(cell, {})[key] = {k: v.get(k) for k in CELL_KEYS}
        res["nearest"][key] = [dict(idx=r["idx"], scheme=r["scheme"], passive=r["passive"], new_value=r["new_value"],
                                    binding=r["new_binding"], in_domain=r.get("in_domain")) for r in
                               m["nearest_to_new_target"][:10]]
    hs = set(res["code_hash"].values())
    res["same_code"] = len(hs) == 1
    return res


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("merged", nargs=3)
    ap.add_argument("--out", default=OUT)
    a = ap.parse_args(argv)
    res = build(a.merged)
    json.dump(res, open(a.out, "w", encoding="utf-8"), ensure_ascii=False, indent=1)
    print(json.dumps(dict(same_code=res["same_code"], complete=res["complete"],
                          stable={k: v["stable_30_100"] for k, v in res["counts"].items()},
                          lowem={k: v["lowem_30_100_1valve"] for k, v in res["counts"].items()}), ensure_ascii=False))
    return 0 if res["same_code"] and all(res["complete"].values()) else 2


if __name__ == "__main__":
    sys.exit(main())
