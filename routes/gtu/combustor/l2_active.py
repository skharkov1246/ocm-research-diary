"""Воронка камеры сгорания KV, уровень 2: активный дорасчёт по модели-заменителю (узел platform/surrogate).

Что делает (навык surrogate-opt, шаги 1-6): по итогам уровня 1 (merge_l1 -> rows) для семейства «схема × тип камеры ×
форсунка» строит problem.json (входы - непрерывные оси comb_space в границах схемы и форсунки; выходы - показатели
уровня 1), таблицу runs.csv (ключ прогона, порция, вариант), вызывает surrogate.py fit (ансамбль MLP + гауссовский
процесс, проверка по блокам) и surrogate.py next (ожидаемое улучшение с вероятностью выполнить ограничения), и
переводит предложенные точки в варианты: дискретные оси - семейства, непрерывные - из next.csv, геометрия - comb_space.
derive. Цель - потери давления (КПД, 49,4 кВт на процент у KV-3 - targets_basis), ограничения - жёсткие цели L1:
NOx_100 ≤ 9, CO_100 ≤ 15, NOx_30 ≤ 15, RTDF ≤ 0,05, T_metal ≤ 1173 K (targets.json).
Предложенные точки сначала идут на уровень 1 (дёшево), прошедшие - на уровень 2; остановка - максимум EI меньше допуска
цели две порции подряд или исчерпан бюджет уровня (навык, шаг 6). Прогноз модели - не результат: число публикуется
только после прямого расчёта (навык, правило 1).

Запуск:
  python l2_active.py comb_l1_KV-3_results.json --inputs l1_inputs_KV-3.json --family "FLOX|выносная|струи" \\
         --k 12 --run-id <run_id> --out l2_active/FLOX_vyn
  python l2_active.py ... --check-plumbing    # проверка связки на любых данных (не результат)
"""
import argparse, csv, json, math, os, subprocess, sys

HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in sys.path:
    sys.path.insert(0, HERE)
SURR = os.path.join(HERE, "..", "platform", "surrogate", "surrogate.py")
CONSTRAINTS = ["NOx_100<=9", "CO_100<=15", "NOx_30<=15", "RTDF<=0.05", "T_metal<=1173"]
OUTPUTS = ["dp_total", "NOx_100", "CO_100", "NOx_30", "RTDF", "T_metal"]


def bounds_for(scheme, injector):
    """границы непрерывных осей семейства - как в comb_space (SCH, INJ)"""
    import comb_space as CS
    S, I = CS.SCH[scheme], CS.INJ[injector]
    b = dict(a_pz=S["a_pz"], dp_inj=S["dp_inj"], tau_pz_ms=S["tau_pz"], tau_bo_ms=(2.0, 20.0), LD=(1.5, 4.0),
             l_mix=I["l_mix"], d_jet_u=(0.0, 1.0))
    if CS.AUX_KIND[scheme]:
        b["a_aux"] = S["a_aux"]
    if I["swirl"][1] > I["swirl"][0]:
        b["swirl"] = I["swirl"]
    return {k: [float(v[0]), float(v[1])] for k, v in b.items()}


def table(rows, ents, names, run_id):
    out = []
    for s in rows:
        e = ents[s["pos"]]
        c = e["concept"]
        lin, c30 = s.get("liner_100") or {}, s.get("cell_30") or {}
        y = dict(dp_total=s.get("dp_total"), NOx_100=lin.get("NOx_ppm15"), CO_100=lin.get("CO_ppm15"),
                 NOx_30=c30.get("NOx_ppm15"), RTDF=lin.get("RTDF"), T_metal=(s.get("wall") or {}).get("T_metal_K"))
        r = {k: c[k] for k in names}
        r.update({k: ("" if v is None else v) for k, v in y.items()})
        r.update(batch=0, run_id=run_id, job=f"l1-pos{s['pos']}", idx=s["idx"])
        out.append(r)
    return out


def write_csv(path, rows):
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)


def surrogate(*args):
    r = subprocess.run([sys.executable, SURR, *args], capture_output=True, text=True, cwd=os.path.dirname(SURR))
    if r.returncode != 0:
        raise SystemExit(f"surrogate.py {args[0]}: код {r.returncode}\n{r.stderr[-2000:]}")
    return r.stdout


def points_to_entries(next_csv, family_entry, names):
    """точки next.csv -> концепции (дискретные оси семейства, непрерывные - из таблицы) с геометрией derive"""
    import comb_space as CS
    import comb_crn as CR
    bc = CR.boundary()
    out = []
    with open(next_csv, encoding="utf-8") as f:
        for k, row in enumerate(csv.DictReader(f)):
            c = dict(family_entry["concept"])
            for n in names:
                c[n] = float(row[n])
            c["idx"], c["origin"] = -1 - k, "l2_active"
            d = CS.derive(c, bc)
            bad = CS.check(c, d, bc)
            out.append(dict(concept=c, ok=not bad, reject=[b[0] for b in bad],
                            pred={kk: row[kk] for kk in row if kk.startswith(("pred_", "sd_", "ei", "pof", "mode"))}))
    return out


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("l1_results")
    ap.add_argument("--inputs", required=True)
    ap.add_argument("--family", help="схема|камера|форсунка; по умолчанию - самое многочисленное")
    ap.add_argument("--k", type=int, default=12, help="точек в порции (по числу параллельных машин × кейсов)")
    ap.add_argument("--run-id", default="local")
    ap.add_argument("--out", required=True)
    ap.add_argument("--check-plumbing", action="store_true", help="проверка связки: модель по любым строкам")
    a = ap.parse_args(argv)
    res = json.load(open(a.l1_results, encoding="utf-8"))
    ents = json.load(open(a.inputs, encoding="utf-8"))["entries"]
    fam = {}
    for s in res["rows"]:
        key = f'{s["scheme"]}|{s["chamber"]}|{s["injector"]}'
        fam.setdefault(key, []).append(s)
    key = a.family or max(fam, key=lambda k: len(fam[k]))
    rows = fam[key] if not a.check_plumbing else res["rows"]
    sch, _, inj = key.split("|")
    b = bounds_for(sch, inj)
    os.makedirs(a.out, exist_ok=True)
    prob = dict(inputs=b, outputs=OUTPUTS)
    pj = os.path.join(a.out, "problem.json")
    json.dump(prob, open(pj, "w", encoding="utf-8"), ensure_ascii=False, indent=1)
    tab = table(rows, ents, list(b), a.run_id)
    if len(tab) < 2 * len(b):
        print(f"[внимание] строк {len(tab)} меньше 2·d = {2 * len(b)}: модель ненадёжна (навык: n ≈ 10·d)")
    rc = os.path.join(a.out, "runs.csv")
    write_csv(rc, tab)
    model = os.path.join(a.out, "model.joblib")
    rep = os.path.join(a.out, "report.json")
    print(surrogate("fit", "--problem", pj, "--data", rc, "--blocks", "kmeans", "--n-blocks", str(min(20, len(tab))),
                    "--folds", str(min(5, max(2, len(tab) // 3))), "--kinds", "mlp,gp", "--jobs", "4", "--out", model,
                    "--report", rep))
    nxt = os.path.join(a.out, "next.csv")
    cons = sum((["--constraint", x] for x in CONSTRAINTS), [])
    print(surrogate("next", "--model", model, "--objective", "dp_total", "--min", *cons, "--k", str(a.k), "--kind", "gp",
                    "--out", nxt))
    ref = ents[rows[0]["pos"]]
    pts = points_to_entries(nxt, ref, list(b))
    json.dump(dict(family=key, run_id=a.run_id, n_train=len(tab), constraints=CONSTRAINTS, objective="dp_total (min)",
                   status="прогноз модели - не результат; точки идут на прямой расчёт L1, затем L2", points=pts,
                   check_plumbing=a.check_plumbing),
              open(os.path.join(a.out, "l2_points.json"), "w", encoding="utf-8"), ensure_ascii=False, indent=1, default=float)
    print(f"семейство {key}: обучение на {len(tab)} строках, предложено {len(pts)}, годных по правилам "
          f"{sum(p['ok'] for p in pts)}; {a.out}")


if __name__ == "__main__":
    main()
