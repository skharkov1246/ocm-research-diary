"""Файлы задач AWS устойчивого отбора и фронтов Парето ниш (формат раннера ../../platform/aws/runner.py v2) со сметой и plan.

  python make_jobs_uq.py      # jobs_uq_debug.json, jobs_uq.json, jobs_pareto.json, jobs_uq_pareto.json, READY_UQ.json

Порядок (каждый шаг - после того, как вход лежит в ветке):
  1. jobs_uq_debug.json - Монте-Карло отладочного уровня 1 (aws/l1_debug, 150 вариантов): время на выборку при 64 процессах;
  2. jobs_uq.json - Монте-Карло всех вариантов полного уровня 1 KV-3 (aws/l1_full, 3983): хрупкость оптимумов уровня 1, запасы по
     правилам, нужные для 90 % выборок, главные неопределённости - вход для режима устойчивости фронтов;
  3. jobs_pareto.json - NSGA-II в каждой нише отбора уровня 0: KV-3 (затравка - варианты уровня 1), затем KV-1, KV-6; в каждой точке -
     номинальный фронт и устойчивый (ограничение «допустим в 75 % из 16 выборок»);
  4. jobs_uq_pareto.json - Монте-Карло объединённого набора (варианты уровня 1 + проверенные члены фронтов) -> merge_uq: 50 устойчивых
     разных и сравнение с обычным отбором.
Смета: процессорное время на вычисление - из локальных миниатюр (uq_probe_results.json, pareto_probe*_results.json) и пробы уровня 1
(нерасчётные режимы, funnel_l1_results.json) с запасом SAFETY; доли вариантов, доходящих до нерасчётных режимов, и число кандидатов
фронтов - допущения (поле assume у каждой группы, «оценка»). Сторож max_min каждой машины - 2 × расчётное время + установка + 10 мин;
верхняя граница - правило раннера. AWS не вызывается.
"""
import argparse, json, math, os, sys

HERE = os.path.dirname(os.path.abspath(__file__))
FUNNEL = os.path.dirname(HERE)
REPO = os.path.dirname(os.path.dirname(os.path.dirname(FUNNEL)))
for p in (FUNNEL, HERE):
    if p not in sys.path:
        sys.path.insert(0, p)
import make_jobs_l0 as J0        # noqa: E402
import make_jobs_l1 as J1        # noqa: E402

REL = "routes/gtu/funnel/uq"
FREL = "routes/gtu/funnel"
SAFETY = 1.5
CEILING, SPENT_NOTE = 2000.0, "израсходовано около 3 долл. (задание 15.11.2026; точное - runner.py status по ledger.jsonl)"
MAX_COMPUTE_MIN = 150         # расчёт на машину не дольше (оценка: сторож вдвое больше - запас на медленную машину)
TAIL = "\nrc=$?\nls -l \"$GTU_OUT\" > \"$GTU_OUT/_parts_list.txt\"\nexit $rc"
RUNS = "routes/gtu/platform/aws/runs/<run_id>"
SELECT = {"KV-3": f"{FREL}/l0_aws_full/funnel_l0_select_KV-3_n200000.json", "KV-1": f"{FREL}/l0_aws_full/funnel_l0_select_KV-1_n50000.json",
          "KV-6": f"{FREL}/l0_aws_full/funnel_l0_select_KV-6_n50000.json"}
L1_FULL = dict(dir=f"{FREL}/aws/l1_full", M=3983, shards=40)
L1_DEBUG = dict(dir=f"{FREL}/aws/l1_debug", M=150, shards=10)
NOM = dict(pop=200, gens=100, check=24)
ROB = dict(pop=120, gens=80, check=24, robust=16, p_target=0.75, tag="_rob16")
ROB_TOP = {"KV-1": 400}       # KV-1: устойчивый режим - первые 400 ниш по лучшей оценке уровня 0 (из 1864; по пользе, оценка)
ASSUME = dict(od_frac_l1=0.10, p_dp_od_l1=0.4, od_frac_debug=0.2, cands_nom_per_niche=8, cands_rob_per_niche=6, od_frac_pareto=0.3,
              p_dp_od_pareto=0.75,
              basis="оценка: доли вариантов, доходящих до нерасчётных режимов, и P(расчётная точка) у них - по миниатюре (у оптимумов уровня 1 "
                    "P около 0,03-0,08, порог 0,2 проходят немногие); проверенных допустимых членов фронта на нишу - по миниатюре 2 из 2 "
                    "(на AWS проверяется до 24); фактические числа даст jobs_uq_debug.json и первая волна")


def _load(name):
    p = os.path.join(HERE, name)
    return json.load(open(p, encoding="utf-8")) if os.path.exists(p) else None


def measured():
    """процессорное время, с: расчётная точка под неопределённостью, вычисление NSGA-II (номинальное и устойчивое), проверка члена
    фронта (нерасчётные режимы), нерасчётные режимы варианта уровня 1"""
    u = _load("uq_probe_results.json") or {}
    pn = _load("pareto_probe_results.json") or {}
    pr = _load("pareto_probe_rob_results.json") or {}

    def per(doc, key):
        v = [((n.get("t_cpu") or {}).get(key)) for n in doc.get("niches", [])]
        v = [x for x in v if x]
        return sum(v) / len(v) if v else None
    m1 = J1.measured(os.path.join(FUNNEL, "funnel_l1_results.json"))
    out = dict(dp_eval=(u.get("timing") or {}).get("cpu_s_per_dp_eval") or 0.025,
               nsga_eval=per(pn, "per_eval") or 0.012, nsga_eval_rob=per(pr, "per_eval") or 0.02,
               check=per(pn, "per_check") or m1["t_od_mean"], od_eval=m1["t_od_mean"],
               seed_niche=per(pn, "seed") or 8.0,
               basis=("локальные миниатюры на 4 общих ядрах (uq_probe_results.json, pareto_probe_results.json, pareto_probe_rob_results.json) и "
                      "проба уровня 1 (нерасчётные режимы - funnel_l1_results.json, " + m1["basis"] + ")"))
    return {k: (round(v, 5) if isinstance(v, float) else v) for k, v in out.items()}


def niches_of(kv):
    import pareto_niche as PN
    rows, meta = PN.load_niches(os.path.join(REPO, SELECT[kv]))
    return len(rows)


def machines_for(core_s, shards, cap=None):
    m = max(1, math.ceil(core_s / (J0.VCPU * 60 * MAX_COMPUTE_MIN)))
    return min(m, shards, cap or shards)


def group_jobs(name, cmd_of, shards, machines, core_s):
    per_min = core_s / machines / J0.VCPU / 60
    max_min = max(30, int(math.ceil(J0.SETUP_MIN + 2 * per_min + 10)))
    jobs = []
    for j in range(machines):
        mine = [k for k in range(shards) if k % machines == j]
        jobs.append(dict(name=f"{name}-m{j}", command=cmd_of(",".join(map(str, mine))) + TAIL, max_min=max_min))
    return jobs, dict(core_hours=round(core_s / 3600, 2), machines=machines, compute_min_per_machine=round(per_min, 1), max_min=max_min)


def usd(jobs, core_s, gb=0.05):
    nmach = len(jobs)
    per = core_s / max(1, nmach) / J0.VCPU / 60
    mh = nmach * (J0.SETUP_MIN + per + J0.UPLOAD_MIN) / 60
    ebs = J0.EBS_GB * J0.EBS_USD_GB_MONTH / 730
    P = J0.PRICE
    u = dict(spot_expected=round(mh * P["expected"], 3), spot_alt=round(mh * P["alt"], 3), disk_ipv4=round(mh * (ebs + J0.IPV4_USD_H), 4),
             s3_month=round(gb * J0.S3_USD_GB_MONTH, 4), transfer=round(gb * J0.EGRESS_USD_GB, 4))
    u["total_expected"] = round(u["spot_alt"] + u["disk_ipv4"] + u["s3_month"] + u["transfer"], 3)
    bound = sum((P["upper"] + ebs + J0.IPV4_USD_H) * (jb["max_min"] + 10) / 60 for jb in jobs)
    return dict(machine_hours=round(mh, 3), usd=u, upper_bound_usd=round(bound, 2))


# =========================================================================== группы
def uq_group(name, kv, M, n, n_od, od_frac, p_dp, shards, l1=None, cand_files=(), cands=None, pdp=0.2, smin=4.5, m=None):
    K = cands if cands is not None else M
    core = SAFETY * (K * (n + 1) * m["dp_eval"] + K * od_frac * n_od * p_dp * m["od_eval"])
    parts = []
    if l1:
        parts.append(f"--l1-in {l1['dir']} --l1-kv {kv} --l1-m {l1['M']} --l1-shards {l1['shards']}")
    if cand_files:
        parts.append("--candidates " + " ".join(cand_files))
    src = " ".join(parts)

    def cmd(sh):
        return (f"python {REL}/uq_l1.py {src} --n {n} --n-od {n_od} --od auto --od-pdp-min {pdp} --od-score-min {smin} --shards {shards} "
                f"--shard {sh} --workers \"$NPROC\" --gzip --out \"$GTU_OUT\"")
    jobs, info = group_jobs(name, cmd, shards, machines_for(core, shards, 6 if K < 5000 else None), core)
    tag = name.replace("uq-", "")
    merge = (f"python3 {REL}/merge_uq.py --kv {kv} --shards {shards} --in {RUNS} --out-dir {REL}/aws/{tag} --suffix _{tag}")
    info.update(kv=kv, candidates_assumed=K, n=n, n_od=n_od, od_rule=dict(pdp_min=pdp, score_min=smin),
                assume=dict(od_frac=od_frac, p_dp_at_od=p_dp), core_s_basis="SAFETY × (K·(N + 1)·t_расч + K·доля·N_od·P·t_нерасч)")
    return jobs, info, merge, core


def pareto_group(name, kv, NN, P, shards, m, l1=None, top=None):
    nn = min(NN, top) if top else NN
    evals = P["pop"] * (P["gens"] + 1)
    t_eval = evals * m["nsga_eval"] + (evals * P["robust"] * m["dp_eval"] if P.get("robust") else 0.0)   # устойчивый: каждое построенное × K
    core = SAFETY * nn * (t_eval + P["check"] * m["check"] + m["seed_niche"])
    extra = (f" --l1-in {l1['dir']} --l1-m {l1['M']} --l1-shards {l1['shards']} --l1-optional" if l1 else "") + \
            (f" --robust {P['robust']} --p-target {P['p_target']} --tag {P['tag']}" if P.get("robust") else "") + (f" --top {top}" if top else "")

    def cmd(sh):
        return (f"python {REL}/pareto_niche.py --select {SELECT[kv]} --pop {P['pop']} --gens {P['gens']} --check {P['check']}{extra} "
                f"--shards {shards} --shard {sh} --workers \"$NPROC\" --gzip --out \"$GTU_OUT\"")
    jobs, info = group_jobs(name, cmd, shards, machines_for(core, shards), core)
    tg = P.get("tag", "")
    merge = (f"python3 {REL}/merge_pareto.py --kv {kv} --nn {nn} --shards {shards}" + (f" --tag {tg}" if tg else "") +
             f" --in {RUNS} --out-dir {REL}/aws/pareto --suffix {tg}_{kv}")
    info.update(kv=kv, niches=nn, pop=P["pop"], gens=P["gens"], evals_per_niche=evals, robust=P.get("robust", 0), p_target=P.get("p_target"),
                check=P["check"], core_s_basis="SAFETY × ниши × (вычисления × t + проверка × t_нерасч + затравка)" +
                                                (" ; устойчивый: каждое построенное решение ещё × K выборок" if P.get("robust") else ""))
    return jobs, info, merge, core


def doc_for(groups, comment, extra=None):
    jobs = [j for g in groups for j in g[0]]
    core = sum(g[3] for g in groups)
    est = usd(jobs, core)
    est.update(core_hours=round(core / 3600, 2), machines=len(jobs), groups={g[0][0]["name"].rsplit("-m", 1)[0]: g[1] for g in groups},
               safety=SAFETY, price=J0.PRICE, upper_bound_rule="раннер: (MaxPrice + диск + IPv4) × (max_min + 10) / 60 на машину")
    d = dict(_comment=comment, defaults=dict(instance_type="c7a.16xlarge", fallback_types=["c7i.16xlarge"], image="al2023", disk_gb=J0.EBS_GB,
                                             python="3.12", pip=J0.PIP),
             jobs=jobs, merge={g[0][0]["name"].rsplit("-m", 1)[0]: g[2] for g in groups}, estimate=est)
    if extra:
        d.update(extra)
    return d


def build(m):
    out = {}
    g = uq_group("uq-l1_debug", "KV-3", L1_DEBUG["M"], 500, 100, ASSUME["od_frac_debug"], 0.4, 10, l1=L1_DEBUG, m=m)
    out["jobs_uq_debug.json"] = doc_for([g], "Устойчивый отбор, отладка: Монте-Карло 150 вариантов отладочного уровня 1 (aws/l1_debug - после "
                                             "прогона jobs_l1_debug.json). Сгенерировано make_jobs_uq.py.",
                                        dict(needs=[f"{L1_DEBUG['dir']}/funnel_l1_KV-3_m150_s*of010.json.gz в ветке"]))
    g = uq_group("uq-l1_full", "KV-3", L1_FULL["M"], 500, 100, ASSUME["od_frac_l1"], ASSUME["p_dp_od_l1"], 50, l1=L1_FULL, m=m)
    out["jobs_uq.json"] = doc_for([g], "Устойчивый отбор: Монте-Карло всех вариантов полного уровня 1 KV-3 (N = 500 выборок, нерасчётные - "
                                       "100 первых у прошедших порог). Сгенерировано make_jobs_uq.py.",
                                  dict(needs=[f"{L1_FULL['dir']}/funnel_l1_KV-3_m3983_s*of040.json.gz в ветке (полный уровень 1)"]))
    pg, NNs = [], {}
    for kv in ("KV-3", "KV-1", "KV-6"):
        NNs[kv] = niches_of(kv)
        sh = 50 if NNs[kv] >= 200 else 20
        l1 = L1_FULL if kv == "KV-3" else None
        pg.append(pareto_group(f"par-{kv.lower()}", kv, NNs[kv], NOM, sh, m, l1=l1))
        pg.append(pareto_group(f"par-{kv.lower()}-rob16", kv, NNs[kv], ROB, sh, m, l1=l1, top=ROB_TOP.get(kv)))
    out["jobs_pareto.json"] = doc_for(pg, "Фронты Парето ниш (NSGA-II): KV-3, затем KV-1, KV-6; номинальный и устойчивый режимы. Раннер ставит "
                                          "задачи по порядку очереди под квоту. Сгенерировано make_jobs_uq.py.",
                                      dict(needs=["отборы уровня 0 l0_aws_full/*.json - уже в ветке",
                                                  f"KV-3: затравка - {L1_FULL['dir']} (--l1-optional: без частей ниши идут без затравки)"],
                                           niches=NNs))
    ug = []
    for kv in ("KV-3", "KV-1", "KV-6"):
        nn_rob = min(NNs[kv], ROB_TOP.get(kv, NNs[kv]))
        K = int(NNs[kv] * ASSUME["cands_nom_per_niche"] + nn_rob * ASSUME["cands_rob_per_niche"]) + (L1_FULL["M"] if kv == "KV-3" else 0)
        files = [f"{REL}/aws/pareto/pareto_cands_{kv}.json.gz", f"{REL}/aws/pareto/pareto_cands_rob16_{kv}.json.gz"]
        sh = 50 if K > 3000 else 20
        ug.append(uq_group(f"uq-union_{kv}", kv, K, 500, 48, ASSUME["od_frac_pareto"], ASSUME["p_dp_od_pareto"], sh,
                           l1=L1_FULL if kv == "KV-3" else None, cand_files=files, cands=K, pdp=0.5, smin=5.0, m=m))
    out["jobs_uq_pareto.json"] = doc_for(ug, "Устойчивый отбор объединённого набора: варианты уровня 1 (KV-3) + проверенные члены фронтов "
                                             "(номинальных и устойчивых) - N = 500, нерасчётные - 48 первых у P(расчётная) >= 0,5 и P50 >= 5. "
                                             "Сгенерировано make_jobs_uq.py.",
                                         dict(needs=[f"{REL}/aws/pareto/pareto_cands*_<KV>.json.gz в ветке (merge_pareto.py после "
                                                     "jobs_pareto.json)", f"KV-3: {L1_FULL['dir']}"]))
    return out


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--plan-out", default=os.environ.get("TMPDIR", "/tmp"))
    ap.add_argument("--no-plan", action="store_true")
    a = ap.parse_args(argv)
    m = measured()
    docs = build(m)
    tot_e = tot_u = 0.0
    ready = dict(title="Устойчивый отбор и фронты Парето ниш: готовность к AWS", date="2026-11-15", measured=m, assume=ASSUME, files={},
                 order=["1. jobs_uq_debug.json (после l1_debug) -> merge_uq -> время на выборку при 64 процессах; при расхождении с миниатюрой "
                        "больше чем вдвое - перегенерировать make_jobs_uq.py",
                        "2. jobs_uq.json (после l1_full) -> merge_uq: хрупкость оптимумов уровня 1 и запасы по правилам (rules_margins)",
                        "3. jobs_pareto.json -> merge_pareto по каждой группе (merge.*) -> pareto_cands*.json.gz в ветку",
                        "4. jobs_uq_pareto.json -> merge_uq по каждой точке: uq_robust50_*.json - 50 устойчивых разных"])
    for name, d in docs.items():
        path = os.path.join(HERE, name)
        with open(path, "w", encoding="utf-8") as fh:
            json.dump(d, fh, ensure_ascii=False, indent=1)
        if not a.no_plan:
            d["plan_check"] = J0.plan_check(path, a.plan_out)
            with open(path, "w", encoding="utf-8") as fh:
                json.dump(d, fh, ensure_ascii=False, indent=1)
        e = d["estimate"]
        tot_e += e["usd"]["total_expected"]; tot_u += e["upper_bound_usd"]
        ready["files"][name] = dict(machines=e["machines"], core_hours=e["core_hours"], expected_usd=e["usd"]["total_expected"],
                                    upper_bound_usd=e["upper_bound_usd"], plan_exit=(d.get("plan_check") or {}).get("exit_code"))
        print(f"{name}: {e['machines']} машин, {e['core_hours']} ядро-ч, ожидаемо {e['usd']['total_expected']} долл., верхняя граница "
              f"{e['upper_bound_usd']} долл.; plan: {(d.get('plan_check') or {}).get('exit_code')}")
    ready.update(total_expected_usd=round(tot_e, 2), total_upper_bound_usd=round(tot_u, 2), ceiling_usd=CEILING, spent_note=SPENT_NOTE,
                 within_ceiling=tot_u < CEILING - 50, ready_for_aws=all((v["plan_exit"] in (0, None)) for v in ready["files"].values()),
                 how_to_regenerate=f"python3 {REL}/make_jobs_uq.py --plan-out <папка>")
    with open(os.path.join(HERE, "READY_UQ.json"), "w", encoding="utf-8") as fh:
        json.dump(ready, fh, ensure_ascii=False, indent=1)
    print(f"READY_UQ.json: всего ожидаемо {ready['total_expected_usd']} долл., верхняя граница {ready['total_upper_bound_usd']} долл. "
          f"(потолок {CEILING}, {SPENT_NOTE})")
    return 0


if __name__ == "__main__":
    sys.exit(main())
