"""Покрытие дискретных осей уровня 0 после ремонта выборки v2 (замечания критика полноты, 15.11.2026).

  python3 coverage_l0.py [--workers 4]      # -> funnel_l0_coverage.json (минуты на 4 ядрах; AWS не нужен)

Что считается (код уровня 0 не меняется - отпечаток code_sha прежний, входы уровня 1 принимаются):
1. По каждой точке (KV-3, KV-1, KV-6; слитые прогоны AWS l0_aws_full/, выборка v2) и архетипу: какие значения осей ниши
   (роторы, компрессор, турбина, охлаждение, рекуператор, изготовление ротора) разрешает шаблон ARCH и какие остались среди
   допустимых ниш. Для каждого пропавшего значения - перебор всех сочетаний шаблона с этим значением функцией
   funnel_space.closure (тот же перебор, что closure_table, но по всем архетипам, а не только без допустимых): сколько
   сочетаний замыкают точку и чего не хватает. Итог - причина пропажи: «ремонт уводит» (ни одно сочетание не замыкает - выборка
   его не порождает), «правила check()» (порождается, но все концепции отсеяны) или «не встретилось» (замыкает, но в выборку не
   попало).
2. KV-1: допустимые ниши выборки v1 (git f47b2cb, прогон 20261114T232200Z-7ce3), которых нет среди допустимых v2; для каждой -
   сколько сочетаний её шаблона замыкают точку по closure() v2. «closure строже check»: в v1 ниша имела допустимую концепцию по
   правилам check(), а closure() v2 не замыкает ни одного её сочетания - v2 такую нишу не порождает.
Всё - оценка уровня 0 (README узла, разделы «Выборка v2…» и «Что не закрыто»).
"""
import argparse, hashlib, itertools, json, os, subprocess, sys, time
from collections import Counter, defaultdict
from concurrent.futures import ProcessPoolExecutor

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))
sys.path.insert(0, HERE)
import funnel_space as FS          # noqa: E402
import funnel_l0 as L0             # noqa: E402

DATE = "2026-11-14"
FILES = {"KV-3": "l0_aws_full/funnel_l0_results_KV-3_n200000.json", "KV-1": "l0_aws_full/funnel_l0_results_KV-1_n50000.json",
         "KV-6": "l0_aws_full/funnel_l0_results_KV-6_n50000.json"}
V1 = dict(commit="f47b2cb", run_id="20261114T232200Z-7ce3", file="routes/gtu/funnel/l0_aws_full/funnel_l0_results_KV-1_n50000.json")
AXES = ("rotors", "comp", "turb", "cool", "recup", "rotor_make")
WHY = {"w": "мала удельная работа", "EGT": "горячий выхлоп", "T5": "газ перед рекуператором выше 700 °C", "T_al": "алюминиевое колесо выше 200 °C",
       "turb": "турбина на валу пакета не существует в пределах ψ, U/C0, AN², втулки", "AN2": "AN² последней ступени выше 6,03e7", "perf": "цикл не решается"}


def parse_niche(s):
    p = s.split("|")
    return dict(arch=p[0], rotors=int(p[1].rstrip("р")), comp=p[2], turb=p[3], cool=p[4], recup=p[5], rotor_make=p[6])


def combos(arch, kv):
    """все сочетания осей цикла шаблона архетипа - тот же перебор, что funnel_space.closure_table (детерминирован)"""
    T = FS.ARCH[arch]["allowed"]
    L = lambda a: list(T.get(a, FS._DEF[a]))
    seen = set()
    for cores, rot, comp, turb, cool, recup, comb, m_tw, rj in itertools.product(
            L("cores"), L("rotors"), L("comp"), L("turb"), L("cool"), L("recup"), L("comb"), L("m_tw"), L("rotor_join")):
        rot = 2 if cores > 1 else rot
        rj = "сборный" if m_tw == "керамика" else rj
        for n_ax in (L("n_ax") if comp == "ОЦ" else [0]):
            c = dict(arch=arch, kv=kv, cores=cores, rotors=rot, comp=comp, n_ax=n_ax, turb=turb, cool=cool, recup=recup, comb=comb,
                     m_tw=m_tw, rotor_join=rj)
            if rj in ("заодно", "сварка"):
                c["m_cw"], c["mat_cw"] = m_tw, "Ni"
                mcws = [(m_tw, "Ni")]
            else:
                # колесо компрессора - все способы шаблона (ниша различает изготовление ротора «сборный:m_cw/m_tw»)
                mcws = [(m, "Ti64" if "Ti64" in L("mat_cw") else [x for x in L("mat_cw") if x != "Ni"][0]) for m in L("m_cw")]
            for m_cw, mat_cw in mcws:
                c2 = dict(c, m_cw=m_cw, mat_cw=mat_cw)
                for mt in FS._mat_tw_list(c2, T):
                    cc = dict(c2, mat_tw=mt)
                    k = (FS._sig(cc), L0.rotor_make(cc))
                    if k in seen:
                        continue
                    seen.add(k)
                    yield cc


def _eval_arch(args):
    arch, kv = args
    rows = []
    for cc in combos(arch, kv):
        r = FS.closure(cc)
        b = r.get("best") or {}
        rows.append(dict(arch=arch, rotors=cc["rotors"], comp=cc["comp"], turb=cc["turb"], cool=cc["cool"], recup=cc["recup"],
                         rotor_make=L0.rotor_make(cc), cores=cc["cores"], comb=cc["comb"], mat_tw=cc["mat_tw"], ok=bool(r["ok"]),
                         why=None if r["ok"] else r.get("why", "perf"), reason=None if r["ok"] else r.get("reason"), m=b.get("m")))
    return arch, kv, rows


def _summ(rows):
    ok = sum(1 for r in rows if r["ok"])
    why = Counter(r["why"] for r in rows if not r["ok"])
    fails = [r for r in rows if not r["ok"]]
    near = max(fails, key=lambda r: (r["m"] if r["m"] is not None else -9)) if fails else None
    return dict(combos=len(rows), closable=ok, fail_reasons={WHY.get(k, k): v for k, v in why.most_common()},
                nearest_fail=(dict(combo={a: near[a] for a in AXES + ("cores", "comb", "mat_tw")}, reason=near["reason"]) if near else None))


def coverage(kv, res, table):
    """по архетипам: значения осей ниши в шаблоне, среди допустимых, пропавшие и почему"""
    arc = [dict(a, **parse_niche(a["niche"])) for a in res["niches"]["archive"]]
    out = {}
    for arch in FS.ARCH_KEYS:
        rows = table.get(arch, [])
        mine = [a for a in arc if a["arch"] == arch]
        feas = [a for a in mine if a["feasible"] > 0]
        ba = res["summary"]["by_arch"].get(arch, {})
        rep = (ba.get("sampler") or {}).get("repaired")
        d = dict(n=ba.get("n"), feasible=ba.get("feasible"), repaired=rep,
                 repaired_share=(round(rep / ba["n"], 3) if rep is not None and ba.get("n") else None),
                 first_reject=ba.get("first_reject"), niches_seen=len(mine), niches_feasible=len(feas),
                 niches_feasible_list=sorted(a["niche"] for a in feas), closure=_summ(rows), axes={})
        for ax in AXES:
            tmpl = sorted({r[ax] for r in rows}, key=str)
            have = sorted({a[ax] for a in feas}, key=str)
            lost = []
            for v in tmpl:
                if v in have:
                    continue
                sub = [r for r in rows if r[ax] == v]
                seen = [a for a in mine if a[ax] == v]
                s = _summ(sub)
                if s["closable"] == 0:
                    verdict = "ремонт уводит: ни одно сочетание шаблона с этим значением не замыкает точку - выборка его не порождает"
                elif seen:
                    verdict = (f"правила check(): порождается ({sum(a['n'] for a in seen)} концепций в {len(seen)} нишах), допустимых нет - "
                               "причины отсева по концепциям в CSV прогона (S3)")
                else:
                    verdict = "замыкает, но в выборку не попало (ремонт выбирает ближайшую замыкающую замену, другие оси меняются первыми)"
                lost.append(dict(value=v, combos=s["combos"], closable=s["closable"], seen_concepts=sum(a["n"] for a in seen),
                                 seen_niches=len(seen), verdict=verdict, fail_reasons=s["fail_reasons"], nearest_fail=s["nearest_fail"]))
            d["axes"][ax] = dict(template=tmpl, feasible=have, lost=lost)
        out[arch] = d
    return out


def v1_lost(res2, table, v1doc):
    a2 = {a["niche"]: a for a in res2["niches"]["archive"]}
    f2 = {k for k, a in a2.items() if a["feasible"] > 0}
    f1 = {a["niche"]: a for a in v1doc["niches"]["archive"] if a["feasible"] > 0}
    lost = []
    for k, a in sorted(f1.items(), key=lambda kv: -kv[1]["best_score"]):
        if k in f2:
            continue
        nz = parse_niche(k)
        rows = [r for r in table.get(nz["arch"], []) if all(r[x] == nz[x] for x in AXES)]
        s = _summ(rows)
        seen = a2.get(k)
        if s["combos"] == 0:
            verdict = "сочетания нет в шаблоне v2 (шаблон архетипа изменён)"
        elif s["closable"] == 0:
            verdict = "closure строже check: в v1 допустима по правилам check(), closure() v2 не замыкает ни одного сочетания - v2 нишу не порождает"
        elif seen:
            verdict = f"порождается v2 ({seen['n']} концепций), допустимых нет - правила check() при окнах v2"
        else:
            verdict = "замыкает по closure(), но в выборку v2 не попало"
        lost.append(dict(niche=k, v1_feasible=a["feasible"], v1_best_score=a["best_score"], v1_best_eta_el=a.get("best_eta_el"),
                         v2_seen=(seen or {}).get("n", 0), v2_feasible=(seen or {}).get("feasible", 0), combos=s["combos"],
                         closable=s["closable"], verdict=verdict, fail_reasons=s["fail_reasons"],
                         nearest_fail=(s["nearest_fail"] or {}).get("reason")))
    by = Counter(x["verdict"].split(":")[0] for x in lost)
    return dict(v1=dict(commit=V1["commit"], run_id=V1["run_id"], file=V1["file"], feasible_niches=len(f1)),
                v2_feasible_niches=len(f2), lost=len(lost), lost_not_seen_in_v2=sum(1 for x in lost if x["v2_seen"] == 0),
                by_verdict=dict(by), niches=lost,
                note="ниши v1 с 1-3 допустимыми концепциями из 50 000 - редкие; v2 даёт больше допустимых ниш (столбец v2_feasible_niches), "
                     "но эти теряет: закрыть - сверить closure() с check() (README, «Что не закрыто»)")


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--workers", type=int, default=min(4, os.cpu_count() or 1))
    ap.add_argument("--out", default=os.path.join(HERE, "funnel_l0_coverage.json"))
    a = ap.parse_args(argv)
    t0 = time.time()
    res = {kv: json.load(open(os.path.join(HERE, f), encoding="utf-8")) for kv, f in FILES.items()}
    for kv, d in res.items():
        if d.get("code_sha") not in (None, L0.code_sha()) and d.get("merge_checks", {}).get("code_sha") not in (None, L0.code_sha()):
            print(f"[внимание] {kv}: результаты другой версии кода уровня 0", file=sys.stderr)
    tasks = [(arch, kv) for kv in FILES for arch in FS.ARCH_KEYS]
    table = defaultdict(dict)
    with ProcessPoolExecutor(max(1, a.workers)) as ex:
        for arch, kv, rows in ex.map(_eval_arch, tasks):
            table[kv][arch] = rows
    v1doc = json.loads(subprocess.run(["git", "-C", REPO, "show", f"{V1['commit']}:{V1['file']}"], check=True, capture_output=True).stdout)
    cov = {kv: coverage(kv, res[kv], table[kv]) for kv in FILES}
    lost = v1_lost(res["KV-1"], table["KV-1"], v1doc)
    sha = lambda p: hashlib.sha256(open(p, "rb").read()).hexdigest()
    doc = dict(title="Уровень 0: покрытие осей ниш после ремонта выборки v2 и ниши v1, потерянные v2", date=DATE, status="оценка",
               method=("значения осей ниши в шаблоне ARCH против допустимых ниш слитых прогонов AWS (выборка v2); пропавшие значения - перебор "
                       "сочетаний шаблона funnel_space.closure (номинальный цикл с лучшими непрерывными осями; тот же перебор, что "
                       "closure_table); KV-1 - допустимые ниши v1 (git f47b2cb), которых нет среди допустимых v2, и closure() v2 по их "
                       "сочетаниям. Код уровня 0 не менялся (code_sha " + L0.code_sha() + ")"),
               software=f"python {sys.version.split()[0]}", code_sha_l0=L0.code_sha(),
               inputs=[dict(file=f"routes/gtu/funnel/{f}", sha256=sha(os.path.join(HERE, f))) for f in FILES.values()]
               + [dict(file=V1["file"], git=V1["commit"], field="выборка v1, прогон " + V1["run_id"]),
                  dict(file="routes/gtu/funnel/funnel_space.py", sha256=sha(os.path.join(HERE, "funnel_space.py")))],
               axes=list(AXES), coverage=cov, kv1_v1_lost_niches=lost,
               summary={kv: {arch: dict(feasible=c["feasible"], niches_feasible=c["niches_feasible"], repaired_share=c["repaired_share"],
                                         closable=f"{c['closure']['closable']}/{c['closure']['combos']}",
                                         lost_values={ax: [x["value"] for x in c["axes"][ax]["lost"]] for ax in AXES if c["axes"][ax]["lost"]})
                             for arch, c in cov[kv].items()} for kv in FILES},
               seconds=round(time.time() - t0, 1))
    with open(a.out, "w", encoding="utf-8") as fh:
        json.dump(doc, fh, ensure_ascii=False, indent=1)
    print(f"{a.out}: {doc['seconds']} с; KV-1 потеряно ниш v1: {lost['lost']} ({lost['by_verdict']})")
    return 0


if __name__ == "__main__":
    sys.exit(main())
