"""Камера KV-3, уровень 0 v2: файл готовности к AWS (READY_V2.json) после исправлений проверок скептиков 15.11.2026.

Все числа - из файлов: пакет задач (jobs_aws_comb_l0_v2.json: смета, plan, мини-проверка, команды слияния), проба
N = 300 (comb_l0_v2_probe_results.json), чувствительность к β (comb_l0_v2_beta_sens_results.json), отпечатки кода,
данных и границ - из кода ветки. Итог тестов передаётся строкой (--tests) из вывода pytest.
Запуск: python ready_v2.py --tests "test_comb_l0_v2.py: N passed; test_comb_l0.py: M passed" [--tests-ok]
AWS не вызывается.
"""
import argparse, json, os, re, sys

HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in sys.path:
    sys.path.insert(0, HERE)
OUT = os.path.join(HERE, "READY_V2.json")
JOBS = os.path.join(HERE, "jobs_aws_comb_l0_v2.json")
PROBE = os.path.join(HERE, "comb_l0_v2_probe_results.json")
BETA = os.path.join(HERE, "comb_l0_v2_beta_sens_results.json")
CEILING = 2000.0

FIXES = [
    dict(src="физика", what="пилот «гидравлика» - время пребывания = tau_pz (часть смесителей того же купола); своё "
                            "время - только у полости TVC", status="принято"),
    dict(src="физика", what="эжектор - меньшее из уноса Ricou - Spalding и струйного насоса в трубе не шире ячейки "
                            "элемента (баланс импульса, zeta = 1 - оценка, трение, выброс в жаровую трубу)",
         status="принято"),
    dict(src="физика", what="жёсткость подачи в индексе пульсаций - по отклику расхода (закон контура)", status="принято"),
    dict(src="физика", what="разброс сопротивления пористой вставки -20 % (оценка) - правило ОБРАТНЫЙ_ТОК_РАЗБРОС",
         status="принято"),
    dict(src="физика", what="холостой ход 0 % (экстраполяция цикла, оценка): ОБРАТНЫЙ_ТОК и запас опорной зоны в "
                            "load_range; новая цель - 30-100 % (решение владельца 2)", status="принято"),
    dict(src="физика", what="README: предел роста доли пилота - закон расхода, не правило обратного тока; физический "
                            "довод против каталитической схемы (2,84 раза при пределе 1,9)", status="принято"),
    dict(src="физика", what="узлы неравномерности без прибавки топлива при срезе бедного узла", status="принято"),
    dict(src="физика", what="запас по срыву ниже 1 - лестница до 1000 времён и 7 шагов бисекции (v1: 1/3, 1/10, 1/30)",
         status="принято"),
    dict(src="физика", what="времена зон на частичной нагрузке - при неизменном объёме (было постоянное время)",
         status="принято"),
    dict(src="физика", what="каналы: поджигание группы B газом группы A", status="частично",
         why="истина между пределами: отбор - без поджигания, верхняя оценка с поджиганием - рядом (gas_cross, "
             "value_cross, счётчики слияния); проверка - уровень 1"),
    dict(src="физика", what="варианты вне области корреляций (S по корреляции выше среза, C/C_OPT, расслоение)",
         status="частично", why="границы области - оценка: пометка in_domain и счётчики «в области облика» вместо "
                                "отсева; такие варианты уровень 1 (RANS) проверяет первыми"),
    dict(src="код", what="сброс состояния Cantera перед каждым вариантом (итог не зависит от порядка в пуле)",
         status="принято"),
    dict(src="код", what="data_hash (targets.json, targets_basis_results.json) в частях и в сверке слияния",
         status="принято"),
    dict(src="код", what="слияние: только файлы частей (не итоги, не _lost), повтор части с другим содержимым - "
                         "проблема полноты, код 3 при строках «ошибка»", status="принято"),
]

LIMITS = [
    "все числа уровня 0 - облик (0D + корреляции), оценка для отбора на уровень 1, не прогноз конструкции",
    "закон затухания неравномерности и расслоения при малом J (BETA_LOWJ) - оценка; серии β = 0 и 2 - чувствительность",
    "холостой ход - экстраполяция цикла по 50 и 30 % (модель цикла ниже доли топлива 0,005 не считает)",
    "времена зон при неизменном объёме: плотность - равновесие по таблице 100 % × p3(L)/p3(100), сдвиг T3 не учтён",
    "эжектор: потери возвратного канала zeta = 1 и скорость в трубе ~ скорость канала - оценка",
    "разброс пористой вставки -20 % и константа жёсткости 0,03 p3 - оценка; TRL элементов - оценка",
    "каналы: отбор без поджигания группы B газом группы A (нижний предел), верхняя оценка - рядом",
    "область облика (S до 0,6, C/C_OPT 0,25-4, доля воздуха с топливом >= 0,5) - оценка; вне неё - пометка без отсева",
    "стенка - только на 100 %; горячее пятно пилота на малой нагрузке не моделируется",
    "слиянию 2^19 строк нужно не меньше 6 ГБ свободной памяти (экстраполяция замера, оценка)",
]


def build(tests, tests_ok):
    import comb_v2 as V
    import comb_l0_v2 as L
    import comb_space as CS
    pkg = json.load(open(JOBS, encoding="utf-8"))
    probe = json.load(open(PROBE, encoding="utf-8"))
    beta = json.load(open(BETA, encoding="utf-8")) if os.path.exists(BETA) else None
    code, data, bc = V.code_hash(), L.data_hash(), V.boundary(write=False)["key"]
    m = re.search(r"code_hash v2 ([0-9a-f]{16})", pkg["_comment"])
    pkg_code = m.group(1) if m else None
    ref = probe["completeness"]["ref"]
    est = pkg["estimate"]
    pc, mc = pkg.get("plan_check", {}), pkg.get("mini_check", {})
    checks = dict(
        tests=dict(summary=tests, ok=bool(tests_ok)),
        probe=dict(file=os.path.basename(PROBE), complete=probe["complete"], code_hash=ref["code_hash"],
                   data_hash=ref.get("data_hash"), bc_key=ref["bc_key"], same_code=ref["code_hash"] == code,
                   counts={k: probe["counts"][k] for k in ("n", "feasible", "stable_30_100", "stable_30_100_in_domain",
                                                            "stable_0_100", "lowem_30_100_1valve", "hard_ok",
                                                            "load_range_base", "reject_codes", "cross_ignition_upper")},
                   cpu_variant_mean_s=probe["timing"]["cpu_variant_mean_s"]),
        beta_sensitivity=(dict(file=os.path.basename(BETA), same_code=beta.get("same_code"),
                               code_hash=beta["code_hash"], complete=beta.get("complete"),
                               stable_30_100={k: v["stable_30_100"] for k, v in beta["counts"].items()},
                               lowem={k: v["lowem_30_100_1valve"] for k, v in beta["counts"].items()})
                          if beta else None),
        plan=dict(command=pc.get("command"), exit_code=pc.get("exit_code"), forecast_usd=pc.get("forecast_usd")),
        mini=dict(n=mc.get("n"), jobs_exit=[j["exit_code"] for j in mc.get("jobs", [])], merge=mc.get("merge"),
                  merged=mc.get("merged")),
        package_code_hash=pkg_code, package_same_code=pkg_code == code)
    mini_ok = bool(mc) and all(j["exit_code"] == 0 for j in mc.get("jobs", [])) and \
        (mc.get("merge") or {}).get("exit_code") == 0 and (mc.get("merged") or {}).get("complete") is True and \
        (mc.get("merged") or {}).get("code_hash") == code
    ready = bool(tests_ok and pc.get("exit_code") == 0 and mini_ok and probe["complete"] and ref["code_hash"] == code
                 and pkg_code == code and (beta is None or (beta.get("same_code") and all(beta["complete"].values())
                                                            and code in beta["code_hash"].values())))
    main, sens = est["main"], est["beta_sensitivity"]
    order = [dict(order=1, series="основная", jobs=[j["name"] for j in pkg["jobs"] if re.fullmatch(r"comb-l0v2-m\d+", j["name"])],
                  n=main["n"], shards=main["shards"], beta=V.BETA_LOWJ_NOMINAL, max_min=main["max_min"],
                  why="доля «держит 30-100 %» и держатели новой цели до сходимости (правила К1-К3)")]
    for k, e in enumerate(sens, 2):
        order.append(dict(order=k, series=f"чувствительность β = {e['beta']:g}",
                          jobs=[j["name"] for j in pkg["jobs"] if j["name"].startswith(f"comb-l0v2-b{e['beta']:g}-")],
                          n=e["n"], shards=e["shards"], beta=e["beta"], max_min=e["max_min"],
                          why="закон затухания при малом J - оценка; серия показывает, держится ли вывод"))
    jobs_rel = "routes/gtu/combustor/" + os.path.basename(JOBS)
    res = dict(
        schema="kv.combustor.ready/1",
        title="Камера сгорания KV-3, уровень 0 v2: готовность к AWS после исправлений проверок скептиков 15.11.2026",
        date=V.DATE, ready_for_aws=ready,
        ready_note=("готово: тесты, проба, plan (код 0) и мини-проверка команд прошли на одном коде" if ready else
                    "не готово: см. checks"),
        code=dict(code_hash_v2=code, data_hash=data, bc_key=bc, v1_code_hash=CS.code_hash(),
                  v1_note="comb_space.py и comb_crn.py не менялись: прогон AWS v1 20261115T001745Z-f1da (262 144, "
                          "отпечаток fae37b9aadf85961) остаётся итогом v1; исправления этого выпуска - только v2",
                  commit_first=["routes/gtu/combustor/" + f for f in ("comb_v2.py", "comb_l0_v2.py", "merge_l0_v2.py",
                                "comb_bc_v2_results.json", "targets_v2.json", "test_comb_l0_v2.py", "beta_sens_v2.py",
                                "make_jobs_l0_v2.py", "ready_v2.py", "jobs_aws_comb_l0_v2.json",
                                "comb_l0_v2_probe_results.json", "comb_l0_v2_beta_sens_results.json", "README.md",
                                "READY_V2.json")],
                  commit_note="после коммита: python3 -c 'import sys; sys.path.insert(0, \"routes/gtu/combustor\"); "
                              "import comb_v2 as V; print(V.code_hash())' должен дать code_hash_v2; иначе - заново "
                              "make_jobs_l0_v2.py --mini и ready_v2.py"),
        checks=checks,
        series_order=order,
        commands=dict(
            plan=pc.get("command"),
            run_check=f"python3 routes/gtu/platform/aws/runner.py run {jobs_rel} --branch <ветка с коммитом>",
            run_live=f"python3 routes/gtu/platform/aws/runner.py run {jobs_rel} --branch <ветка с коммитом> --live",
            status="python3 routes/gtu/platform/aws/runner.py status",
            merge=pkg["merge"]["commands"], merge_expect=pkg["merge"]["expect"],
            beta_sens_after_merge=("python3 routes/gtu/combustor/beta_sens_v2.py <слияние b0> <слияние основной> "
                                   "<слияние b2> --out routes/gtu/combustor/l0_aws_v2/comb_l0v2_beta_sens_results.json "
                                   "(основная серия при β = 1 - своя выборка N, серии β - N/4: сравнение долей, не строк)")),
        estimate=dict(core_hours=round(main["core_hours"] + sum(e["core_hours"] for e in sens), 1),
                      core_hours_by_series=[round(main["core_hours"], 1)] + [round(e["core_hours"], 1) for e in sens],
                      machine_hours_expected=round(main["machine_hours_expected"] +
                                                   sum(e["machine_hours_expected"] for e in sens), 2),
                      usd_expected=est["usd_expected"], usd_expected_basis=est["usd_expected_basis"],
                      usd_plan_forecast_upper=pc.get("forecast_usd"), ceiling_usd=CEILING,
                      share_of_ceiling=round((pc.get("forecast_usd") or 0.0) / CEILING, 4),
                      time_per_variant=est["time_per_variant"],
                      note="CPU-серия в пределах 2000 долл. - без отдельного вопроса (решение владельца 15.11.2026); "
                           "отчёт после каждых 500 долл. - навык aws-gtu-runs"),
        fixes=FIXES, rejected=[], limitations=LIMITS,
        after_aws=["слияние трёх серий (команды merge), сходимость доли по префиксам 2^k",
                   "если держателей нет и на 2^19 - уровень 1 (RANS) для ближайших вариантов каждой ячейки, вне области "
                   "облика - первыми; верхняя оценка каналов (value_cross) - отдельный список",
                   "решение владельца о постановке - только окном, через основную сессию"])
    return res


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--tests", required=True, help="итог pytest строкой")
    ap.add_argument("--tests-ok", action="store_true", help="все тесты прошли")
    ap.add_argument("--out", default=OUT)
    a = ap.parse_args(argv)
    res = build(a.tests, a.tests_ok)
    json.dump(res, open(a.out, "w", encoding="utf-8"), ensure_ascii=False, indent=1)
    print(json.dumps(dict(ready_for_aws=res["ready_for_aws"], code=res["code"]["code_hash_v2"],
                          plan=res["checks"]["plan"]["exit_code"], usd=res["estimate"]["usd_expected"],
                          forecast=res["estimate"]["usd_plan_forecast_upper"]), ensure_ascii=False))
    return 0 if res["ready_for_aws"] else 2


if __name__ == "__main__":
    sys.exit(main())
