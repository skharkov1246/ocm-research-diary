"""Файлы задач AWS уровня 1 воронки (формат раннера ../platform/aws/runner.py v2) со сметой и проверкой командой plan.

  python make_jobs_l1.py        # jobs_l1_debug.json, jobs_l1.json, jobs_l1_fea.json, READY_L1.json
  python make_jobs_l1.py --v3   # v3 (решения владельца 15.11.2026): jobs_l1_v3.json и jobs_l1_v3_fea.json по отбору уровня 0 v3
                                # (l0_aws_v3/funnel_l0_select_KV-3_n200000.json); v2-файлы не трогает

Вход уровня 1 - отбор слитого уровня 0 (merge_l0): отладка - l0_aws_debug/funnel_l0_select_KV-3_n2000.json (выборка v2, уже в
ветке, 307 вариантов), полный - l0_aws_full/funnel_l0_select_KV-3_n200000.json (выборка v2, прогон 20261115T010625Z-8893, в
ветке, 3983 варианта; funnel_l1.py откажет, если отпечаток выборки или кода уровня 0 не тот). Число вариантов M частей -
min(--top, записей во входе): от него зависят имена частей и команда слияния. Части - единственное деление по номеру варианта j
(порядок ниш по кругу), группируются в машины (<= 6 - квота 384 vCPU). Смета: машино-часы × тариф спота + диск + IPv4 + хранение
+ передача; время на вариант - процессорное время короткой пробы (funnel_l1_results.json, timing.cpu_s_per_variant) в одном
процессе, пересчитанное на число стартов и вычислений, с запасом SAFETY.
Объём и старты выбраны по пользе (сходимость отбора - таблица convergence и задача l1-starts2x), бюджет уровней 0 + 1 - до 1500 долл.
AWS не вызывается.
"""
import argparse, json, math, os, re, subprocess, sys

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))
RUNNER = os.path.join(REPO, "routes", "gtu", "platform", "aws", "runner.py")
REL = "routes/gtu/funnel"
sys.path.insert(0, HERE)
import make_jobs_l0 as J0        # noqa: E402

PIP = J0.PIP
PRICE = J0.PRICE
VCPU = J0.VCPU
SETUP_MIN, UPLOAD_MIN = J0.SETUP_MIN, J0.UPLOAD_MIN
SAFETY = 1.5                     # запас на процессорное время пробы (до 15.11.2026 - ×3 на настенное время общего контейнера, которое
#                                  само завышено чужой нагрузкой в 2,5 раза - замер скептика отбора; смета была завышена около 7 раз)
BUDGET_L01 = 1500.0
PROFILES = {
    "debug": dict(file="jobs_l1_debug.json", input=f"{REL}/l0_aws_debug/funnel_l0_select_KV-3_n2000.json", M=150, shards=10, machines=2,
                  starts=4, maxfev=200,
                  why="отладка уровня 1 на AWS: около 150 вариантов отбора N = 2000 (выборка v2) - время на вариант при 64 процессах, доли "
                      "допустимых, отбор 50 разных; вход уже в ветке"),
    "full": dict(file="jobs_l1.json", input=f"{REL}/l0_aws_full/funnel_l0_select_KV-3_n200000.json", M=4000, shards=40, machines=5,
                 starts=4, maxfev=200, conv=dict(M=640, starts=8, shards=10),
                 why="полный уровень 1: отбор слитого уровня 0 KV-3 N = 200 000 (прогон 20261115T010625Z-8893, merge_l0 --target 4000 -> "
                     "3983 записи); 4 старта (DE 25 поколений + доводка COBYLA 200); проверка сходимости - те же первые 640 вариантов с 8 "
                     "стартами (задача l1-starts2x)"),
}

PROFILES["v3"] = dict(file="jobs_l1_v3.json", input=f"{REL}/l0_aws_v3/funnel_l0_select_KV-3_n200000.json", M=4000, shards=40, machines=5,
                      starts=4, maxfev=200, conv=dict(M=640, starts=8, shards=10), fea_dir="l1_v3_fea", v3=True,
                      why="уровень 1 v3: отбор слитого уровня 0 v3 KV-3 N = 200 000 (jobs_l0_v3.json, merge_l0 --target 4000; ниша - с осью "
                          "«корпус», вариант корпуса - из записи отбора); итоговая оценка × (0,5 + 0,5·TRL/9) и в цели оптимизации; 4 старта "
                          "(DE 25 поколений + доводка COBYLA 200); проверка сходимости - первые 640 вариантов с 8 стартами (l1-starts2x)")
RESCORE_V2 = (f"python3 {REL}/merge_l1.py --rescore --kv KV-3 --m 3983 --shards 40 --in routes/gtu/platform/aws/runs/<run_id_l1_full_v2> "
              f"--out-dir {REL}/aws/l1_full_rescore_v3 --suffix _KV-3_m3983_rescore_v3 --fea-results {REL}/aws/l1_fea/funnel_l1_fea_results.json")


def measured(results=os.path.join(HERE, "funnel_l1_results.json")):
    """время на вариант в одном процессе (с) при стартах и вычислениях пробы и доля оптимизации: процессорное время (меньше
    зависит от чужой нагрузки на общем контейнере), если записано; иначе настенное"""
    d = json.load(open(results, encoding="utf-8"))
    t = dict(d["timing"]["s_per_variant"])
    cpu = d["timing"].get("cpu_s_per_variant") or {}
    basis = "настенное время пробы"
    if cpu.get("mean"):
        k = cpu["mean"] / t["mean"]
        t = dict(mean=cpu["mean"], max=cpu["max"])
        basis = f"процессорное время пробы (настенное × {k:.2f})"
        dd = dict(d.get("timing_detail") or {})
        if dd.get("t_opt_mean") is not None:
            dd["t_opt_mean"] = dd["t_opt_mean"] * k
            dd["t_od_mean"] = dd["t_od_mean"] * k
        d["timing_detail"] = dd
    rows = dict(d.get("timing_detail") or {})
    if not rows.get("evals_per_start"):
        import csv
        recs = list(csv.DictReader(open(os.path.join(os.path.dirname(results), d["files"]["csv"]), encoding="utf-8")))
        ev = [float(r["evals"]) / float(r["starts"]) for r in recs if r.get("evals") and r.get("starts")]
        rows["evals_per_start"] = sum(ev) / len(ev) if ev else None
    return dict(t_mean=t["mean"], t_max=t["max"], starts=d["starts"], maxfev=d["maxfev"], opt_share=rows.get("opt_share", 0.6),
                t_opt_mean=rows.get("t_opt_mean"), t_od_mean=rows.get("t_od_mean"), evals_per_start=rows.get("evals_per_start"),
                de_gens=d.get("de_gens"), basis=basis, code_sha=d.get("code_sha"))


def entries_available(input_rel):
    """число вариантов во входе (если файл отбора уже в ветке): M частей = min(--top, число записей) - иначе имена частей
    (m<M>) и команда слияния разойдутся (до 15.11.2026: --m 4000 при 3983 записях - слияние не нашло бы частей)"""
    path = os.path.join(REPO, input_rel)
    if not os.path.exists(path):
        return None
    import funnel_l1 as L1
    sel, problems = L1.load_select(path, strict=False)
    return dict(n=len(L1.order_entries(sel)), problems=problems)


def t_variant(m, starts, maxfev):
    """пересчёт на другие старты: оптимизация ∝ starts (DE с теми же поколениями; доводка COBYLA - малая доля), нерасчётные режимы -
    постоянны; доводка пересчитывается по доле maxfev в вычислениях старта (оценка)"""
    t_opt = m["t_opt_mean"] if m.get("t_opt_mean") is not None else m["t_mean"] * m["opt_share"]
    t_od = m["t_od_mean"] if m.get("t_od_mean") is not None else m["t_mean"] - t_opt
    ev_start = m.get("evals_per_start") or 3000.0
    k_fev = (ev_start - m["maxfev"] + maxfev) / ev_start
    return t_opt * starts / m["starts"] * k_fev + t_od


def build(profile, branch="claude/lucid-allen-yvps5d"):
    P = dict(PROFILES[profile])
    av = entries_available(P["input"])
    if av is not None:
        if av["problems"]:
            raise ValueError(f"{profile}: вход {P['input']} не принят уровнем 1: {av['problems']}")
        P["M_requested"], P["M"] = P["M"], min(P["M"], av["n"])
        P["input_entries"] = av["n"]
    S, Mm = P["shards"], P["machines"]
    lo, hi = (10, 50) if profile == "debug" else (20, 50)
    if not lo <= S <= hi:
        raise ValueError(f"{profile}: частей {S}, нужно {lo}-{hi}")
    jobs = []
    m = measured()
    tv = t_variant(m, P["starts"], P["maxfev"]) * SAFETY
    comp_min = P["M"] * tv / Mm / VCPU / 60
    P["max_min"] = max(30, int(math.ceil(SETUP_MIN + 2 * comp_min + 10)))
    # Раннер выполняет команду как `bash cmd.sh` без -e: код выхода задачи - код последней строки. Код funnel_l1.py сохраняется
    # и возвращается после списка частей - иначе упавший расчёт дал бы exit_code 0 и status ok в summary.json (проверка 15.11.2026).
    tail = "\nrc=$?\nls -l \"$GTU_OUT\" > \"$GTU_OUT/_parts_list.txt\"\nexit $rc"
    for j in range(Mm):
        mine = [k for k in range(S) if k % Mm == j]
        cmd = (f"python {REL}/funnel_l1.py --input {P['input']} --top {P['M']} --shards {S} --shard {','.join(map(str, mine))} "
               f"--starts {P['starts']} --maxfev {P['maxfev']} --workers \"$NPROC\" --gzip --out \"$GTU_OUT\"" + tail)
        jobs.append(dict(name=f"l1-{profile}-m{j}", command=cmd, max_min=P["max_min"]))
    conv = P.get("conv")
    if conv:
        tvc = t_variant(m, conv["starts"], P["maxfev"]) * SAFETY
        mm = max(30, int(math.ceil(SETUP_MIN + 2 * conv["M"] * tvc / VCPU / 60 + 10)))
        jobs.append(dict(name="l1-starts2x", max_min=mm,
                         command=(f"python {REL}/funnel_l1.py --input {P['input']} --top {conv['M']} --shards {conv['shards']} --shard all "
                                  f"--starts {conv['starts']} --maxfev {P['maxfev']} --workers \"$NPROC\" --gzip --tag _starts2x --out \"$GTU_OUT\""
                                  + tail)))
        P["conv_max_min"] = mm
    base = (f"python3 {REL}/merge_l1.py --kv KV-3 --m {P['M']} --shards {S} --in routes/gtu/platform/aws/runs/<run_id> "
            f"--out-dir {REL}/aws/l1_{profile} --suffix _KV-3_m{P['M']} --fea-top 0" + (" --compare-in routes/gtu/platform/aws/runs/<run_id>" if conv else ""))
    fea_dir = P.get("fea_dir", "l1_fea")
    # отладка - без масштабирования на KV-1/KV-6: оно повторяет оптимизацию отобранных (до 50 × 2 точки × starts) в процессе слияния,
    # то есть на машине сессии, а не на AWS (решение владельца: локально - не больше 12 вариантов и 2 стартов)
    merge = dict(first=base + (f" --no-scale --fea-cases-out {REL}/aws/l1_{profile}/funnel_l1_fea_cases.json" if conv else " --no-scale"),
                 final=(base + f" --fea-results {REL}/aws/{fea_dir}/funnel_l1_fea_results.json") if conv else None)
    return P, jobs, merge, tv, m


def estimate(P, jobs, tv, m):
    core_s = P["M"] * tv
    if P.get("conv"):
        core_s += P["conv"]["M"] * t_variant(m, P["conv"]["starts"], P["maxfev"]) * SAFETY
    nmach = len(jobs)
    per_machine = core_s / nmach / VCPU / 60
    machine_min = SETUP_MIN + per_machine + UPLOAD_MIN
    mh = nmach * machine_min / 60
    ebs = J0.EBS_GB * J0.EBS_USD_GB_MONTH / 730
    gb = P["M"] * 25e3 / 1e9                                  # оценка: разбор варианта в части ~25 кБ сжатого JSON
    usd = dict(spot_expected=round(mh * PRICE["expected"], 3), spot_alt=round(mh * PRICE["alt"], 3), disk_ipv4=round(mh * (ebs + J0.IPV4_USD_H), 4),
               s3_month=round(gb * J0.S3_USD_GB_MONTH, 4), transfer=round(gb * J0.EGRESS_USD_GB, 4))
    usd["total_expected"] = round(usd["spot_alt"] + usd["disk_ipv4"] + usd["s3_month"] + usd["transfer"], 3)
    bound = sum((PRICE["upper"] + ebs + J0.IPV4_USD_H) * (jb["max_min"] + 10) / 60 for jb in jobs)
    return dict(variants=P["M"], starts=P["starts"], maxfev=P["maxfev"], s_per_variant_core_with_safety=round(tv, 2),
                measured=m, safety=SAFETY, time_basis=m.get("basis"), core_hours=round(core_s / 3600, 2), machines=nmach, shards=P["shards"],
                compute_min_per_machine=round(per_machine, 1), machine_hours=round(mh, 3), usd=usd, upper_bound_usd=round(bound, 2),
                upper_bound_rule="раннер: (MaxPrice + диск + IPv4) × (max_min + 10) / 60 на машину", price=PRICE)


def fea_jobs(l1_dir="l1_full", name="l1-fea300"):
    return dict(_comment=("CalculiX (disk.py) для лучших 300 вариантов уровня 1: вход - funnel_l1_fea_cases.json, который пишет merge_l1.py "
                          "--fea-cases-out (положить в ветку до запуска). Контейнер ubuntu:24.04 с calculix-ccx и numpy; ccx - в один поток, "
                          "случаи параллельно ($NPROC)."),
                defaults=dict(instance_type="c7a.16xlarge", fallback_types=["c7i.16xlarge"], image="al2023", disk_gb=30),
                jobs=[dict(name=name, container="ubuntu:24.04",
                           container_setup=["apt-get update -q",
                                            "apt-get install -y -q --no-install-recommends calculix-ccx python3 python3-numpy"],
                           command=f"python3 {REL}/strength_l1.py --run-cases {REL}/aws/{l1_dir}/funnel_l1_fea_cases.json --workers \"$NPROC\" --out \"$GTU_OUT\"",
                           max_min=40)],
                estimate=dict(cases="до 300 вариантов × 1-3 колеса", s_per_case="5-20 с (disk.py, сетка 40×6, один поток - оценка по примерам)",
                              usd_expected=round(0.5 * PRICE["alt"], 2),
                              usd_expected_basis="оценка: машина около 30 мин (как make_jobs_l23.py, строка jobs_l1_fea.json сметы)",
                              usd_upper_bound=round((PRICE["upper"] + J0.IPV4_USD_H) * 50 / 60, 2)))


def write_v3(plan_out):
    """jobs_l1_v3.json (по отбору уровня 0 v3) и jobs_l1_v3_fea.json; блоки merge (с пересчётом прогона v2 без нового счёта) и estimate"""
    P, jobs, merge, tv, m = build("v3")
    est = estimate(P, jobs, tv, m)
    pending = P.get("input_entries") is None
    doc = dict(_comment=("Уровень 1 воронки, v3 (README узла funnel, раздел «v3: решения владельца 15.11.2026»). Раннер берёт только defaults "
                         "и jobs; блоки input, merge, estimate, plan_check - для сессии с ключами. Код - ветка claude/lucid-allen-yvps5d. "
                         "Сгенерировано make_jobs_l1.py --v3; правки - в нём."),
               defaults=dict(instance_type="c7a.16xlarge", fallback_types=["c7i.16xlarge"], image="al2023", disk_gb=J0.EBS_GB,
                             python="3.12", pip=PIP),
               jobs=jobs,
               input=dict(file=P["input"], entries=P.get("input_entries"), M=P["M"], M_requested=P.get("M_requested"), pending=pending,
                          note=(("вход ещё не посчитан: отбор уровня 0 v3 даст jobs_l0_v3.json (merge_l0 --target 4000). После слияния - в "
                                 "ветку и ПЕРЕГЕНЕРИРОВАТЬ этот файл (python3 routes/gtu/funnel/make_jobs_l1.py --v3): M частей = min(4000, "
                                 "записей отбора) - иначе имена частей и команда слияния разойдутся; funnel_l1.py откажет, если отпечаток "
                                 "выборки или кода уровня 0 не тот") if pending else
                                "в ветке; отпечатки выборки и кода уровня 0 совпадают с текущими (funnel_l1.load_select)")),
               merge=dict(command=merge["first"], final=merge["final"],
                          where="сессия с ключами после run (раннер скачивает части в routes/gtu/platform/aws/runs/<run_id>/)",
                          then=("первое слияние - без масштабирования, пишет funnel_l1_fea_cases.json -> в ветку -> jobs_l1_v3_fea.json "
                                "(CalculiX лучших 300) -> итог в aws/l1_v3_fea/ -> окончательное слияние (final): отсев по CalculiX до отбора 50 "
                                "и масштабирование 50 на KV-1/KV-6"),
                          rescore_v2_l1_full=dict(command=RESCORE_V2,
                                                  what=("пересчёт отбора 50 по уже посчитанным записям уровня 1 v2 (прогон на коде f5a4bfe, "
                                                        "aws/l1_full/) с множителем TRL - без нового счёта и без AWS: оценка - funnel_l0.final_score "
                                                        "по частным баллам записей; ось «корпус» в записях v2 нет (картридж/разъём по оси "
                                                        "cartridge выборки v2); --m и --shards - как у слияния прогона v2")),
                          profile_why=P["why"]),
               estimate=est)
    path = os.path.join(HERE, P["file"])
    with open(path, "w", encoding="utf-8") as fh:
        json.dump(doc, fh, ensure_ascii=False, indent=1)
    doc["plan_check"] = J0.plan_check(path, plan_out)
    with open(path, "w", encoding="utf-8") as fh:
        json.dump(doc, fh, ensure_ascii=False, indent=1)
    print(f"{P['file']}: {len(jobs)} машин, {P['M']} вариантов{' (вход ждёт уровня 0 v3)' if pending else ''}, {P['starts']} стартов × "
          f"{P['maxfev']}, ожидаемо {est['usd']['total_expected']} долл., верхняя граница {est['upper_bound_usd']} долл.; plan: код "
          f"{doc['plan_check']['exit_code']}")
    fj = fea_jobs("l1_v3", "l1-v3-fea300")
    pf = os.path.join(HERE, "jobs_l1_v3_fea.json")
    with open(pf, "w", encoding="utf-8") as fh:
        json.dump(fj, fh, ensure_ascii=False, indent=1)
    fj["plan_check"] = J0.plan_check(pf, plan_out)
    with open(pf, "w", encoding="utf-8") as fh:
        json.dump(fj, fh, ensure_ascii=False, indent=1)
    print(f"jobs_l1_v3_fea.json: plan код {fj['plan_check']['exit_code']}")
    return doc, fj


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--plan-out", default=os.environ.get("TMPDIR", "/tmp"))
    ap.add_argument("--v3", action="store_true", help="только jobs_l1_v3.json и jobs_l1_v3_fea.json")
    a = ap.parse_args(argv)
    if a.v3:
        doc, fj = write_v3(a.plan_out)
        return 0 if doc["plan_check"]["exit_code"] == 0 and fj["plan_check"]["exit_code"] == 0 else 1
    docs = {}
    for prof in ("debug", "full"):
        P, jobs, merge, tv, m = build(prof)
        est = estimate(P, jobs, tv, m)
        doc = dict(_comment=("Уровень 1 воронки (README узла funnel, раздел «Уровень 1»). Раннер берёт только defaults и jobs; блоки merge, "
                             "estimate, plan_check - для сессии с ключами. Код - ветка claude/lucid-allen-yvps5d (PR #62). "
                             "Сгенерировано make_jobs_l1.py; правки - в нём."),
                   defaults=dict(instance_type="c7a.16xlarge", fallback_types=["c7i.16xlarge"], image="al2023", disk_gb=J0.EBS_GB,
                                 python="3.12", pip=PIP),
                   jobs=jobs, input=dict(file=P["input"], entries=P.get("input_entries"), M=P["M"], M_requested=P.get("M_requested"),
                                         note=("уже в ветке (выборка v2); отпечатки выборки и кода уровня 0 совпадают с текущими "
                                               "(funnel_l1.load_select)")),
                   merge=dict(command=merge["first"], final=merge["final"],
                              where="сессия с ключами после run (раннер скачивает части в routes/gtu/platform/aws/runs/<run_id>/)",
                              then=("первое слияние - без масштабирования (быстро), пишет funnel_l1_fea_cases.json -> в ветку -> jobs_l1_fea.json "
                                    "(CalculiX лучших 300) -> итог в aws/l1_fea/ -> окончательное слияние (final): отсев по CalculiX до отбора 50 и "
                                    "масштабирование 50 на KV-1/KV-6 (RUNBOOK_L0_L1.md)") if prof == "full" else None,
                              profile_why=P["why"]),
                   estimate=est)
        path = os.path.join(HERE, P["file"])
        with open(path, "w", encoding="utf-8") as fh:
            json.dump(doc, fh, ensure_ascii=False, indent=1)
        doc["plan_check"] = J0.plan_check(path, a.plan_out)
        with open(path, "w", encoding="utf-8") as fh:
            json.dump(doc, fh, ensure_ascii=False, indent=1)
        print(f"{P['file']}: {len(jobs)} машин, {P['M']} вариантов, {P['starts']} стартов × {P['maxfev']}, ожидаемо {est['usd']['total_expected']} "
              f"долл., верхняя граница {est['upper_bound_usd']} долл.; plan: код {doc['plan_check']['exit_code']}")
        docs[prof] = doc
    fj = fea_jobs()
    pf = os.path.join(HERE, "jobs_l1_fea.json")
    with open(pf, "w", encoding="utf-8") as fh:
        json.dump(fj, fh, ensure_ascii=False, indent=1)
    fj["plan_check"] = J0.plan_check(pf, a.plan_out)
    with open(pf, "w", encoding="utf-8") as fh:
        json.dump(fj, fh, ensure_ascii=False, indent=1)
    print(f"jobs_l1_fea.json: plan код {fj['plan_check']['exit_code']}")
    write_ready(docs, fj)
    return 0


def write_ready(docs, fj):
    r0 = json.load(open(os.path.join(HERE, "READY_L0.json"), encoding="utf-8"))
    l0 = r0["estimate_total_usd"]
    l1_exp = round(sum(d["estimate"]["usd"]["total_expected"] for d in docs.values()) + fj["estimate"].get("usd_expected", 0.0), 2)
    l1_up = round(sum(d["estimate"]["upper_bound_usd"] for d in docs.values()) + fj["estimate"]["usd_upper_bound"], 2)
    plan_ok = all(d["plan_check"]["exit_code"] == 0 for d in docs.values()) and fj["plan_check"]["exit_code"] == 0
    doc = dict(title="Уровень 1 воронки двигателя-заменителя: готовность к AWS", date="2026-11-14", ready_for_aws=plan_ok,
               files=dict(debug=f"{REL}/jobs_l1_debug.json", full=f"{REL}/jobs_l1.json", fea=f"{REL}/jobs_l1_fea.json"),
               order=["1. уровень 0 прошёл: отладка 20261115T010201Z-1a6e (l0_aws_debug/), полный 20261115T010625Z-8893 (l0_aws_full/, "
                      f"отбор KV-3 - {docs['full']['input'].get('entries')} записей) - оба входа уровня 1 в ветке",
                      f"2. jobs_l1_debug.json ({docs['debug']['estimate']['variants']} вариантов, 2 машины) -> merge_l1 (команда merge.command) -> "
                      "проверить время на вариант при 64 процессах, доли допустимых, отбор 50; при большом расхождении с пробой - "
                      "перегенерировать make_jobs_l1.py",
                      f"3. jobs_l1.json ({docs['full']['estimate']['variants']} вариантов, 5 машин + l1-starts2x) -> merge_l1 --compare-in "
                      "(сходимость по стартам) -> funnel_l1_fea_cases.json в ветку -> jobs_l1_fea.json (CalculiX, лучшие 300) -> итог 50 разных"],
               plan_exit_codes={k: d["plan_check"]["exit_code"] for k, d in docs.items()} | {"fea": fj["plan_check"]["exit_code"]},
               plan_note=docs["full"]["plan_check"]["runner_note"],
               estimate=dict(level0=l0, level1=dict(expected=l1_exp, upper_bound=l1_up,
                                                     by_file={d_: dict(expected=v["estimate"]["usd"]["total_expected"],
                                                                       upper_bound=v["estimate"]["upper_bound_usd"],
                                                                       machine_hours=v["estimate"]["machine_hours"],
                                                                       core_hours=v["estimate"]["core_hours"]) for d_, v in docs.items()}),
                             total_expected=round(l0["expected"] + l1_exp, 2), total_upper_bound=round(l0["upper_bound"] + l1_up, 2),
                             budget_l01_usd=BUDGET_L01, budget_ok=(l0["upper_bound"] + l1_up) < BUDGET_L01,
                             note="уровни 0 + 1 по верхней границе раннера - доли процента бюджета 1500 долл.: объём ограничен пользой "
                                  "(сходимость отбора), а не деньгами; запас бюджета - уровням 2-3, камере и живой геометрии"),
               known_limitations=["время на вариант - процессорное время короткой пробы (12 вариантов, 2 старта) с запасом ×1,5; "
                                  "фактическое даст jobs_l1_debug.json; сторож max_min и предохранитель раннера ограничивают перерасход",
                                  "раннер v2 - в ветке claude/lucid-allen-yvps5d; если к запуску он изменится, повторить plan",
                                  "окончательное слияние (merge.final в jobs_l1.json) масштабирует 50 отобранных на KV-1 и KV-6 повтором "
                                  "оптимизации (100 вариантов × 4 старта) внутри merge_l1.py, то есть на машине сессии, а не на AWS; "
                                  "в смете уровня 1 его нет"],
               how_to_regenerate="python3 routes/gtu/funnel/make_jobs_l1.py")
    with open(os.path.join(HERE, "READY_L1.json"), "w", encoding="utf-8") as fh:
        json.dump(doc, fh, ensure_ascii=False, indent=1)
    print(f"READY_L1.json: ready={plan_ok}, уровни 0+1 ожидаемо {doc['estimate']['total_expected']} долл., верхняя граница "
          f"{doc['estimate']['total_upper_bound']} долл.")


if __name__ == "__main__":
    sys.exit(main())
