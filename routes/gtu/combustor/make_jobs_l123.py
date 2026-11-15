"""Воронка камеры сгорания KV, уровни 1-3, каналы и калибровка: пакеты задач раннера AWS, смета, проверка plan,
сводная смета против общего потолка 2000 долл. вместе с воронкой двигателя и другими пакетами трека.

Пишет: jobs_calib.json (эталон Sandia D до L1), jobs_l1.json (L1, 500 вариантов), jobs_l2.json (L2: элемент 3D с
запасами; ожидающие - теплообмен стенки, распыл, каналы, топология, PRECCINSTA), jobs_l3.json (L3: LES элемента с
возбуждением и без), comb_budget_results.json (сводная смета). Раннер читает только defaults и jobs; блоки estimate,
pending_jobs, merge, plan_check - для людей и сессии с ключами.

Основа сметы (каждое число - с источником):
  * время ячейки·итерации reactingFoam (без химии / с химией, GRI-Mech 3.0 + TDAC) - местная обкатка уровня 1
    (comb_l1_smoke_results.json; apt-1912, 4 общих ядра под нагрузкой 20-34 - оценка сверху для c7a);
  * размеры сеток - генератор (comb_cfd / comb_cfd3d) на 12 вариантах отбора comb_l0_debug_results.json;
  * цены: ожидаемая 0,7064 долл./ч (спот c7a.16xlarge us-east-1c, ledger.jsonl прогона 20261114T231610Z-9b64),
    планируемая 0,75 (верх диапазона задачи владельца 14.11.2026), верхняя граница - формула раннера
    (MaxPrice × (сторож + 10 мин)), её считает runner.py plan;
  * S3 0,023 долл./ГБ·мес, передача 0,09 долл./ГБ, gp3 0,08 долл./ГБ·мес, IPv4 0,005 долл./ч - Price List AWS
    (routes/gtu/platform/batch/aws_rates.toml).
AWS не вызывается. Запуск: python make_jobs_l123.py [--smoke-dir <папка кейсов обкатки>]
"""
import argparse, glob, json, math, os, re, statistics, subprocess, sys, tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
GTU = os.path.dirname(HERE)
RUNNER = os.path.join(GTU, "platform", "aws", "runner.py")
LEDGER = os.path.join(GTU, "platform", "aws", "ledger.jsonl")
SMOKE = os.path.join(HERE, "comb_l1_smoke_results.json")
BUDGET = os.path.join(HERE, "comb_budget_results.json")
DATE = "2026-11-14"
sys.path.insert(0, HERE)

IMG = "opencfd/openfoam-default:2512"      # Docker Hub: тег 2512 от 23.01.2026, 0,48 ГБ (запрос API 14.11.2026)
CSETUP = ["apt-get update -q", "apt-get install -y -q --no-install-recommends python3 python3-numpy curl ca-certificates"]
PRICE = dict(expected=0.7064, planning=0.75)
ZONES = "us-east-1a=0.9554,us-east-1b=0.7478,us-east-1c=0.7064,us-east-1d=0.8576,us-east-1f=1.1804"
RATES = dict(s3_gb_month=0.023, transfer_gb=0.09, ebs_gb_month=0.08, ipv4_h=0.005)
VCPU = 64
SETUP_MIN, UPLOAD_MIN = 15.0, 8.0          # образ 0,48 ГБ + apt + клон; выгрузка (оценка; smoke2 - 61 с без пакетов)
ROOT_ENV = "export OMPI_ALLOW_RUN_AS_ROOT=1 OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1"


def ceil5(x):
    return int(5 * math.ceil(x / 5.0))


# =========================================================================== основа: время и сетки
def smoke_timing(smoke_dir=None):
    """время ячейки·итерации из журналов местной обкатки (или из сохранённого comb_l1_smoke_results.json)"""
    if smoke_dir is None:
        return json.load(open(SMOKE, encoding="utf-8"))["timing"]
    rows = []
    for case in sorted(glob.glob(os.path.join(smoke_dir, "**", "case.json"), recursive=True)):
        d = os.path.dirname(case)
        desc = json.load(open(case, encoding="utf-8"))
        for ph in ("cold", "hot"):
            f = os.path.join(d, f"log.reactingFoam.{ph}")
            if not os.path.exists(f):
                continue
            t = open(f, errors="replace").read()
            et = [float(x) for x in re.findall(r"ExecutionTime = ([0-9.]+) s", t)]
            ti = [float(x) for x in re.findall(r"\nTime = ([0-9.]+)", t)]
            if len(et) > 20 and len(ti) > 20:
                n_it = ti[-1] - ti[9]
                rows.append(dict(case=os.path.relpath(d, smoke_dir), phase=ph, cells=desc["n_cells"], iters=n_it,
                                 cpu_s=et[-1] - et[9], s_per_cell_iter=(et[-1] - et[9]) / n_it / desc["n_cells"],
                                 ended=t.rstrip().endswith("End")))
    tc = [r["s_per_cell_iter"] for r in rows if r["phase"] == "cold"]
    th = [r["s_per_cell_iter"] for r in rows if r["phase"] == "hot"]
    return dict(rows=rows, cold_s_per_cell_iter=max(tc), hot_s_per_cell_iter=max(th),
                cold_median=statistics.median(tc), hot_median=statistics.median(th),
                rule="в смету - наибольшее по кейсам обкатки (между 10-й и последней итерацией журнала)")


def mesh_stats():
    """сетки генератора на вариантах отбора уровня 0 (comb_l0_debug_results.json, 12 вариантов)"""
    import prepare_l1 as P
    import comb_cfd as CF
    import comb_cfd3d as C3
    import comb_crn as CR
    full, compact, meta = P._rows_from([os.path.join(HERE, "comb_l0_debug_results.json")])
    sel = P.select(full, compact, 12)
    bc, K = CR.boundary(), CR.kin()
    ents = [P.entry(i, meta["n"], meta["seed"], "KV-3", bc, K, w, s, l0) for i, w, s, l0 in sel]
    out = dict(cell=[], cell_medium=[], liner=[], element3d_L2=[], element3d_L3_les=[], liner3d_cht=[], idx=[])
    with tempfile.TemporaryDirectory() as t:
        for e in ents:
            h = e["loads"][0]["hot_init"]
            fake = dict(T=h["T_eq"], Y=h["Y"], m_dot_can=(e["concept"]["a_pz"] * e["loads"][0]["W_air"]) / e["n_can"])
            out["idx"].append(e["idx"])
            out["cell"].append(CF.build_cell_case(e, 0, t + "/a")["n_cells"])
            out["cell_medium"].append(CF.build_cell_case(e, 0, t + "/b", "medium")["n_cells"])
            out["liner"].append(CF.build_liner_case(e, fake, t + "/c")["n_cells"])
            out["element3d_L2"].append(C3.build_element3d_case(e, 0, t + "/d")["n_cells"])
            out["element3d_L3_les"].append(C3.build_element3d_case(e, 0, t + "/e", "L3", length_frac=0.35)["n_cells"])
            out["liner3d_cht"].append(C3.build_liner3d_cht_case(e, fake, t + "/f")["n_cells"])
    return {k: (v if k == "idx" else dict(values=v, mean=statistics.fmean(v), max=max(v))) for k, v in out.items()}


# =========================================================================== пакеты
def money(machine_hours, gb_out, disk_gb, months=3.0):
    sp = machine_hours * PRICE["expected"]
    pl = machine_hours * PRICE["planning"]
    disk = machine_hours * disk_gb * RATES["ebs_gb_month"] / 730.0
    ip = machine_hours * RATES["ipv4_h"]
    s3 = gb_out * RATES["s3_gb_month"] * months
    tr = gb_out * RATES["transfer_gb"]
    return dict(spot_expected=sp, spot_planning=pl, disk=disk, ipv4=ip, s3_3months=s3, transfer_to_session=tr,
                total_expected=sp + disk + ip + s3 + tr, total_planning=pl + disk + ip + s3 + tr)


def pkg_l1(T, MS, n=500, shards=25, grid_every=10):
    per = math.ceil(n / shards)
    tc, th = T["cold_s_per_cell_iter"], T["hot_s_per_cell_iter"]
    it = dict(cold=600, hot=3000)
    case_s = lambda cells: cells * (it["cold"] * tc + it["hot"] * th)
    mean_cell, max_cell = MS["cell"]["mean"], MS["cell"]["max"]
    mean_lin, max_lin = MS["liner"]["mean"], MS["liner"]["max"]
    n_med = math.ceil(per / grid_every)
    cpu = per * 3 * case_s(mean_cell) + per * case_s(mean_lin) + n_med * case_s(MS["cell_medium"]["mean"])
    second = 0.3 * per * 3 * mean_cell * 3000 * th           # второй проход 30 % несошедшихся ячеек (+3000 итераций)
    crit = case_s(max_cell) + case_s(max_lin)                 # путь: ячейка 100 % -> труба
    comp_exp = max((cpu + second) / VCPU, crit) / 60.0
    comp_plan = 2.0 * max((cpu + second) / VCPU, crit + case_s(max_cell)) / 60.0
    mm = ceil5(SETUP_MIN + comp_plan + UPLOAD_MIN)
    jobs = [dict(name=f"comb-l1-p{k:02d}",
                 command=f"{ROOT_ENV}\npython3 -u routes/gtu/combustor/comb_l1.py --inputs routes/gtu/combustor/"
                         f"l1_inputs_KV-3.json --n {n} --shard {k} --shards {shards} --out \"$GTU_OUT\" --workers \"$NPROC\" "
                         f"--grid-every {grid_every} --timeout-min {int(0.4 * mm)}") for k in range(shards)]
    mh = shards * (SETUP_MIN + comp_exp + UPLOAD_MIN) / 60.0
    mh_plan = shards * (SETUP_MIN + comp_plan + UPLOAD_MIN) / 60.0
    gb = shards * (per * 4 + n_med) * 0.6e-3                 # итоги кейса ~0,6 МБ (сжатые: логи, сетка, основные поля)
    est = dict(level="L1", variants=n, shards=shards, variants_per_shard=per, cases_per_shard=per * 4 + n_med,
               s_per_cell_iter=dict(cold=tc, hot=th), iterations=it, mean_cells=dict(cell=mean_cell, liner=mean_lin),
               cpu_s_per_shard=cpu + second, critical_path_s=crit, compute_min=dict(expected=comp_exp, planning=comp_plan),
               setup_min=SETUP_MIN, upload_min=UPLOAD_MIN, max_min=mm, core_hours=shards * (cpu + second) / 3600,
               machine_hours=dict(expected=mh, planning=mh_plan), outputs_gb=gb,
               usd=money(mh, gb, 50), usd_planning_hours=money(mh_plan, gb, 50)["total_planning"],
               rule="кейс - одно ядро; часть - пул на 64 ядра; труба после ячейки 100 %; срок сторожа = подготовка + "
                    "2 × max(ЦП части / 64, путь ячейка + труба + ещё ячейка) + выгрузка")
    pk = dict(_comment="Уровень 1 воронки камеры KV-3 (README, «Уровень 1»). Раннер берёт только defaults и jobs. Перед "
                       "запуском: эталон Sandia D пройден (jobs_calib.json), l1_inputs_KV-3.json в ветке (prepare_l1.py "
                       "по частям уровня 0). Сгенерировано make_jobs_l123.py.",
              defaults=dict(instance_type="c7a.16xlarge", fallback_types=["c7i.16xlarge"], image="al2023", disk_gb=50,
                            max_min=mm, container=IMG, container_setup=CSETUP),
              jobs=jobs,
              merge=dict(where="сессия с ключами после run (или collect), Python с cantera==3.2.0",
                         commands=["python3 routes/gtu/combustor/merge_l1.py routes/gtu/platform/aws/runs/<run_id> --inputs "
                                   "routes/gtu/combustor/l1_inputs_KV-3.json --out routes/gtu/combustor/comb_l1_KV-3_results.json"],
                         expect="полнота True, строк 500 из 500; код 2 - потеря части или чужой отпечаток: разбор, не публикация"),
              estimate=est)
    return pk


def pkg_l2(T, MS, n=40, shards=20, np_=8):
    per = math.ceil(n / shards)
    tc, th = T["cold_s_per_cell_iter"], T["hot_s_per_cell_iter"]
    base = lambda c: c * (800 * tc + 3500 * th)
    rst = lambda c: c * 1500 * th
    m, mx = MS["element3d_L2"]["mean"], MS["element3d_L2"]["max"]
    eff = 0.8                                                  # параллельная эффективность 8 ядер на кейс (оценка)
    lbo = lambda c: base(c) + 5 * rst(c)
    fb = lambda c: base(c) + 4 * rst(c)
    cpu_var = 2 * base(m) + lbo(m) + fb(m)
    cpu = per * cpu_var / eff
    crit = lbo(mx) / np_ / eff
    comp_exp = max(cpu / VCPU, lbo(m) / np_ / eff) / 60.0
    comp_plan = 1.5 * max(cpu / VCPU, crit) / 60.0
    mm = min(720, ceil5(SETUP_MIN + comp_plan + UPLOAD_MIN))
    jobs = [dict(name=f"comb-l2-p{k:02d}",
                 command=f"{ROOT_ENV}\npython3 -u routes/gtu/combustor/comb_l23.py --level L2 --inputs routes/gtu/combustor/"
                         f"l2_inputs_KV-3.json --n {n} --shard {k} --shards {shards} --out \"$GTU_OUT\" --np {np_} "
                         f"--workers \"$NPROC\" --timeout-min {int(0.9 * mm)}") for k in range(shards)]
    mh = shards * (SETUP_MIN + comp_exp + UPLOAD_MIN) / 60.0
    mh_plan = shards * (SETUP_MIN + comp_plan + UPLOAD_MIN) / 60.0
    gb = shards * per * 4 * 5e-3
    est = dict(level="L2 (элемент 3D: 100/50/30 %, цепочки запаса по срыву и проскоку)", variants=n, shards=shards,
               np_per_case=np_, parallel_eff=eff, mean_cells=m, max_cells=mx, cpu_s_per_variant=cpu_var / eff,
               compute_min=dict(expected=comp_exp, planning=comp_plan), max_min=mm, core_hours=shards * cpu / 3600,
               machine_hours=dict(expected=mh, planning=mh_plan), outputs_gb=gb, usd=money(mh, gb, 80),
               usd_planning_hours=money(mh_plan, gb, 80)["total_planning"])
    # ожидающие (не готовы к запуску) - со сметой
    cht_core_h = n * MS["liner3d_cht"]["mean"] * 4000 * th * 1.2 / 3600
    spray_core_h = n * 2 * m * 3000 * th * 3.0 / 3600
    act_core_h = 0.3 * shards * cpu / 3600 + 0.07 * 0.0
    pend = [dict(name="comb-l2-cht", what="сектор трубы с теплообменом стенки, chtMultiRegionFoam, 40 вариантов",
                 ready="сетка и области - готово и проверено на apt-1912 (comb_cfd3d.build_liner3d_cht_case, Allmesh: "
                       "газ 73 356 ячеек, стенка 6 504, входы разбавления и охлаждения); НЕ готово: словари областей "
                       "(constant/gas, constant/solid, 0/gas, 0/solid, regionProperties), свойства печатного сплава k(T), cp(T)",
                 core_hours=cht_core_h, basis="40 × средняя сетка × 4000 итераций × 1,2·t_хим (оценка)"),
            dict(name="comb-l2-spray", what="распыл дизеля sprayFoam, элемент 3D, 100 и 50 %",
                 ready="НЕ готово: механизм н-додекана Reitz в формате OpenFOAM (yaml2ck -> chemkinToFoam с правками навыка), "
                       "словари облака (comb_cfd3d.spray_dicts - шаблон)",
                 core_hours=spray_core_h, basis="40 × 2 × средняя сетка × 3000 итераций × 3·t_хим (механизм 100 веществ, оценка)"),
            dict(name="comb-l2-active", what="активный дорасчёт: 3 порции по 12 точек (l2_active.py) -> L1 + L2",
                 ready="связка проверена (l2_active.py --check-plumbing); запуск - после L1",
                 core_hours=act_core_h, basis="30 % счёта L2 (оценка)")]
    pk = dict(_comment="Уровень 2 воронки камеры KV-3 (README, «Уровень 2»). Раннер берёт только defaults и jobs; "
                       "pending_jobs - не готовы к запуску (что дописать - поле ready). Перед запуском: l2_inputs_KV-3.json "
                       "(merge_l1 -> select_l2) в ветке. Сгенерировано make_jobs_l123.py.",
              defaults=dict(instance_type="c7a.16xlarge", fallback_types=["c7i.16xlarge"], image="al2023", disk_gb=80,
                            max_min=mm, container=IMG, container_setup=CSETUP),
              jobs=jobs, pending_jobs=pend, estimate=est)
    return pk


def pkg_l3(T, MS, n=5, np_=32):
    th = T["hot_s_per_cell_iter"]
    m, mx = MS["element3d_L3_les"]["mean"], MS["element3d_L3_les"]["max"]
    steps = 90000                                              # 20 времён протекания зоны пламени при Co 0,5 (оценка)
    s_step = 1.5 * th                                          # PIMPLE 3 коррекции против 2 в RANS (оценка)
    eff = 0.75
    run_s = lambda c: c * steps * s_step / eff
    comp_plan_min = min(700.0, 1.3 * run_s(m) / np_ / 60.0)
    mm = min(720, ceil5(SETUP_MIN + comp_plan_min + UPLOAD_MIN))
    shards = 2 * n                                             # задача = вариант × вид (собственный режим, возбуждение)
    jobs = [dict(name=f"comb-l3-p{k:02d}",
                 command=f"{ROOT_ENV}\npython3 -u routes/gtu/combustor/comb_l23.py --level L3 --inputs routes/gtu/combustor/"
                         f"l3_inputs_KV-3.json --n {n} --shard {k} --shards {shards} --out \"$GTU_OUT\" --np {np_} "
                         f"--workers {np_} --timeout-min {int(0.95 * mm)}") for k in range(shards)]
    core_h = n * 2 * run_s(m) / 3600
    mh = n * 2 * (run_s(m) / np_ / 3600 * (VCPU / np_) ** 0 + (SETUP_MIN + UPLOAD_MIN) / 60.0)
    gb = n * 2 * 0.2
    pre_core_h = 1.0e6 * 45000 * s_step / eff / 3600
    est = dict(level="L3 (LES элемента: собственный режим и с возбуждением)", variants=n, np_per_case=np_, steps=steps,
               s_per_cell_step=s_step, parallel_eff=eff, mean_cells=m, max_cells=mx, core_hours=core_h,
               max_min=mm, machine_hours=dict(expected=mh, planning=1.5 * mh), outputs_gb=gb, usd=money(mh, gb, 100),
               usd_planning_hours=money(1.5 * mh, gb, 100)["total_planning"],
               note="задача - один кейс (вариант × вид) на 32 ядрах машины (половина ядер простаивает - LES на 64 "
                    "ядрах при 0,3 млн ячеек упирается в обмен; оценка); сторож 720 мин - если LES не уложится, раннер "
                    "получит код 124: продолжение с последней записи вторым пакетом (контрольные точки в $GTU_OUT)",
               session_side=["comb_acoustic.py: сеть элементов с F(ω) из LES (identify_ftf) - секунды",
                             "CalculiX: прочность и ползучесть картриджа (навык calculix-rotor) - минуты, AWS не нужен"])
    pend = [dict(name="calib-preccinsta-les", what="эталон PRECCINSTA (LES), устойчивый и неустойчивый режимы",
                 ready="НЕ готово: геометрия (по Meier et al., 2007 - не открыт) и данные DLR/TNF не загружены",
                 core_hours=2 * pre_core_h, basis="1 млн ячеек × 45 000 шагов × 1,5·t_хим / 0,75 на режим (оценка)")]
    pk = dict(_comment="Уровень 3 воронки камеры KV-3 (README, «Уровень 3»). Раннер берёт только defaults и jobs. Перед "
                       "запуском: PRECCINSTA пройден (pending), l3_inputs_KV-3.json (3-5 лучших L2) в ветке.",
              defaults=dict(instance_type="c7a.16xlarge", fallback_types=["c7i.16xlarge"], image="al2023", disk_gb=100,
                            max_min=mm, container=IMG, container_setup=CSETUP),
              jobs=jobs, pending_jobs=pend, estimate=est)
    return pk


def pkg_calib(T):
    th, tc = T["hot_s_per_cell_iter"], T["cold_s_per_cell_iter"]
    cells = 5 * 20 + 5 * 20 + 5 * 70 + 1 * 70 + 5 * 70 + 60 * 70     # blockMeshDict SandiaD_LTS v2512 (подсчёт по блокам)
    cpu = cells * (1500 * tc + 3500 * th)
    mm = ceil5(SETUP_MIN + 3 * cpu / 60.0 + UPLOAD_MIN)
    jobs = [dict(name="calib-sandiaD", max_min=mm,
                 command="bash routes/gtu/combustor/cfd/calib_sandiaD.sh 1.44 &\nbash routes/gtu/combustor/cfd/calib_sandiaD.sh 1.6 &\nwait")]
    mh = (SETUP_MIN + cpu / 60.0 + UPLOAD_MIN) / 60.0
    est = dict(cells=cells, cpu_s_per_case=cpu, cases=2, max_min=mm, machine_hours=dict(expected=mh, planning=3 * mh),
               usd=money(mh, 0.05, 30))
    crit = dict(sandiaD=dict(
        before_run=True, status="критерии проекта (оценка), записаны до прогона; данные TNF - загрузить (tnfworkshop.org)",
        stations="ось и радиальные профили z/D = 7,5; 15; 30; 45; 60; 75 (sampleDict учебника)",
        pass_if=["T на оси: положение максимума в пределах ±10 % по z/D, величина максимума ±100 K",
                 "радиальные профили T и доли смешения на z/D = 15, 30, 45: СКО отклонения ≤ 10 % максимума профиля",
                 "CO: максимум профиля в пределах множителя 1,5",
                 "NO на z/D = 30-60: в пределах множителя 2 (RANS + EDC обычно завышает - по выдержке поиска)"],
        then="тот же отклонение - в results-JSON L1 рядом с числами камеры; не прошёл - L1 не запускается, разбор моделей "
             "(C1 1,6; PaSR вместо EDC)"),
        preccinsta=dict(
        status="критерии проекта (оценка); данные DLR (Meier et al., Combust. Flame 150, 2007 - по выдержке поиска)",
        pass_if=["средние скорости на высотах 1,5; 5; 10; 20; 30 мм над срезом: отклонение ≤ 15 % наибольшей скорости",
                 "длина внутренней зоны обратного тока ±15 %",
                 "неустойчивый режим (φ 0,70, 25,1 кВт): частота автоколебаний ±10 %, размах - тот же порядок"],
        then="не прошёл - LES уровня 3 и сеть элементов не публикуются, только тренды"))
    pk = dict(_comment="Калибровка до перебора (README, «Калибровка»): Sandia D (RANS) - ворота уровня 1. Та же версия "
                       "OpenFOAM и модели, что кейсы L1. Сгенерировано make_jobs_l123.py.",
              defaults=dict(instance_type="c7a.16xlarge", fallback_types=["c7i.16xlarge"], image="al2023", disk_gb=30,
                            max_min=mm, container=IMG, container_setup=CSETUP),
              jobs=jobs, criteria=crit, estimate=est)
    return pk


# =========================================================================== plan и сводная смета
def plan(path):
    r = subprocess.run([sys.executable, RUNNER, "plan", path, "--assume-price", str(PRICE["planning"]),
                        "--assume-zone-prices", ZONES], capture_output=True, text=True, cwd=REPO, timeout=900)
    txt = r.stdout + r.stderr
    m = re.search(r"прогноз пакета ([0-9.]+) долл", txt)
    sp = re.search(r"израсходовано ([0-9.]+)", txt)
    return dict(command=f"python3 routes/gtu/platform/aws/runner.py plan {os.path.relpath(path, REPO)} --assume-price "
                        f"{PRICE['planning']} --assume-zone-prices {ZONES}",
                exit_code=r.returncode, forecast_usd=float(m.group(1)) if m else None,
                spent_in_ledger=float(sp.group(1)) if sp else None,
                lines=[ln for ln in txt.splitlines() if ln.startswith("[") and ("бюджет" in ln or "план" in ln
                                                                                or "внимание" in ln)][-6:])


def _get(d, path):
    for k in path:
        d = d[k]
    return d


def other_packages():
    """сметы других пакетов трека в репозитории (тот же потолок)"""
    out = []
    spec = [
        ("воронка двигателя L0, отладка (запущен 14.11.2026)", "funnel/jobs_l0_debug.json",
         ["estimate", "usd", "total_expected"], ["estimate", "upper_bound_usd"], "расход уже в ledger.jsonl"),
        ("воронка двигателя L0, полный (запущен 14.11.2026)", "funnel/jobs_l0.json",
         ["estimate", "usd", "total_expected"], ["estimate", "upper_bound_usd"], "расход уже в ledger.jsonl"),
        ("камера L0, N = 262 144 (jobs_aws_comb_l0.json)", "combustor/jobs_aws_comb_l0.json",
         ["estimate", "usd", "total_expected_spot_0_706"], ["plan_check", "forecast_usd"], "альтернатива jobs_l0.json"),
        ("камера L0, N = 30 000 (jobs_l0.json)", "combustor/jobs_l0.json", None, ["estimate", "plan", 0, "forecast_usd"],
         "альтернатива jobs_aws_comb_l0.json; в сводную - большее из двух"),
        ("SU2 KV-3, пробный", "platform/su2/jobs_aws_su2_kv3_smoke.json", ["estimate", "usd", "total_expected"],
         ["estimate", "upper_bound_usd"], ""),
        ("SU2 KV-3, medium", "platform/su2/jobs_aws_su2_kv3.json", ["estimate", "usd", "total_expected"],
         ["estimate", "upper_bound_usd"], ""),
        ("SU2 KV-3, fine", "platform/su2/jobs_aws_su2_kv3_fine.json", ["estimate", "usd", "total_expected"],
         ["estimate", "upper_bound_usd"], "")]
    for name, f, pe, pu, note in spec:
        p = os.path.join(GTU, f)
        if not os.path.exists(p):
            continue
        d = json.load(open(p, encoding="utf-8"))
        exp = up = None
        try:
            exp = _get(d, pe) if pe else None
        except (KeyError, IndexError, TypeError):
            pass
        try:
            up = _get(d, pu)
        except (KeyError, IndexError, TypeError):
            pass
        out.append(dict(name=name, file=f"routes/gtu/{f}", expected_usd=exp, upper_usd=up, note=note))
    return out


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--smoke-dir", help="папка кейсов местной обкатки L1 (журналы reactingFoam) - обновить время")
    a = ap.parse_args(argv)
    T = smoke_timing(a.smoke_dir)
    if a.smoke_dir:
        sm = json.load(open(SMOKE, encoding="utf-8")) if os.path.exists(SMOKE) else {}
        sm["timing"] = T
        json.dump(sm, open(SMOKE, "w", encoding="utf-8"), ensure_ascii=False, indent=1)
    MS = mesh_stats()
    pk = dict(calib=pkg_calib(T), l1=pkg_l1(T, MS), l2=pkg_l2(T, MS), l3=pkg_l3(T, MS))
    files = dict(calib="jobs_calib.json", l1="jobs_l1.json", l2="jobs_l2.json", l3="jobs_l3.json")
    checks = {}
    for k, p in pk.items():
        path = os.path.join(HERE, files[k])
        json.dump(p, open(path, "w", encoding="utf-8"), ensure_ascii=False, indent=1)
        checks[k] = plan(path)
        p["plan_check"] = checks[k]
        json.dump(p, open(path, "w", encoding="utf-8"), ensure_ascii=False, indent=1)
    ch = json.load(open(os.path.join(HERE, "channels_plan_results.json"), encoding="utf-8"))
    rows = []
    for k in ("calib", "l1", "l2", "l3"):
        e = pk[k]["estimate"]
        u = e["usd"]
        rows.append(dict(item=f"камера {k.upper()}", file=f"routes/gtu/combustor/{files[k]}", expected_usd=u["total_expected"],
                         planning_usd=e.get("usd_planning_hours", u["total_planning"]), upper_runner_usd=checks[k]["forecast_usd"],
                         plan_exit_code=checks[k]["exit_code"]))
        for pj in pk[k].get("pending_jobs", []):
            mh = pj["core_hours"] / VCPU
            rows.append(dict(item=f"камера {k.upper()} - {pj['name']} (не готово к запуску)", expected_usd=mh * PRICE["expected"],
                             planning_usd=1.5 * mh * PRICE["planning"], upper_runner_usd=None, basis=pj["basis"]))
    for x in ch["channels"] + ch["topology"]:
        rows.append(dict(item=f"каналы/печать - {x['id']} (постановка, сетки не готовы)", expected_usd=x["usd_expected"],
                         planning_usd=x["usd_planning"], upper_runner_usd=None, basis=x["basis"]))
    others = other_packages()
    comb_l0_upper = max([o["upper_usd"] or 0 for o in others if o["name"].startswith("камера L0")] or [0])
    other_upper = sum((o["upper_usd"] or 0) for o in others if not o["name"].startswith("камера L0")) + comb_l0_upper
    other_exp = sum((o["expected_usd"] or 0) for o in others if not o["name"].startswith("камера L0")) + \
        max([o["expected_usd"] or 0 for o in others if o["name"].startswith("камера L0")] or [0])
    spent = checks["l1"]["spent_in_ledger"] or 0.0
    comb_exp = sum(r["expected_usd"] for r in rows)
    comb_plan = sum(r["planning_usd"] for r in rows)
    comb_upper = sum((r["upper_runner_usd"] or r["planning_usd"]) for r in rows)
    total_upper = spent + comb_upper + other_upper
    stops = [
        "0. Перед всем: runner.py plan каждого файла (сделано - plan_check), ветка запушена, ledger.jsonl закоммичен",
        "1. jobs_calib.json (Sandia D): не прошёл критерии - стоп, разбор моделей; уровень 1 не запускается",
        "2. jobs_l1.json: до запуска prepare_l1.py по слитым частям уровня 0 (l1_inputs_KV-3.json в ветке); после - "
        "merge_l1.py: полнота 500/500, флаги сверки Cantera < 10 % вариантов, сходимость > 70 % кейсов; иначе - стоп "
        "(второй проход несошедшихся, разбор)",
        "3. отчёт раннера после каждых 500 долл. общего расхода - решение владельца о продолжении",
        "4. jobs_l2.json: отбор select_l2 (40) + активный дорасчёт; проход L2 - запас по срыву ≥ 1,5 и H2 ≥ 0,3 по объёму "
        "без проскока при выполненных NOx/CO; меньше 3 вариантов прошли - стоп и пересмотр схем (уровень 0)",
        "5. PRECCINSTA (pending) пройден - только тогда jobs_l3.json (3-5 вариантов)",
        "6. L3: размах пульсаций ≤ 1 % p3 (targets.json) по LES и рост мод сети элементов ≤ 0 с F(ω) из LES; "
        "стенка ≤ 1173 K и ресурс 30 000 ч по CalculiX - кандидат в стенд"]
    res = dict(title="Сводная смета воронки камеры KV-3 (уровни 1-3, каналы, калибровка) против общего потолка AWS",
               method="машино-часы × цена спота c7a.16xlarge + диск + IPv4 + хранение S3 + передача; верхняя граница - "
                      "прогноз runner.py plan (MaxPrice × (сторож + 10 мин)); другие пакеты - их поля estimate",
               software=f"python {sys.version.split()[0]}; runner.py plan", date=DATE,
               status="смета - оценка; ни один пакет уровней 1-3 не запускался",
               inputs=[dict(file="routes/gtu/combustor/comb_l1_smoke_results.json", field="timing", used=dict(
                   cold=T["cold_s_per_cell_iter"], hot=T["hot_s_per_cell_iter"])),
                       dict(file="routes/gtu/platform/aws/ledger.jsonl", field="израсходовано (строка plan)", used=spent),
                       dict(file="routes/gtu/platform/batch/aws_rates.toml", field="s3, transfer, ebs, ipv4", used=RATES),
                       dict(file="routes/gtu/combustor/channels_plan_results.json", field="channels, topology")],
               prices=dict(PRICE, zones=ZONES, basis="0,7064 - us-east-1c (ledger 20261114T231610Z-9b64); 0,75 - верх "
                                                     "диапазона задачи владельца; зоны - ledger того же прогона"),
               mesh_stats=MS, timing=T, combustor_rows=rows, other_packages=others,
               totals=dict(spent_ledger=spent, combustor_expected=comb_exp, combustor_planning=comb_plan,
                           combustor_upper=comb_upper, other_expected=other_exp, other_upper=other_upper,
                           all_upper=total_upper, ceiling=2000.0, headroom_upper=2000.0 - total_upper,
                           share_upper=total_upper / 2000.0,
                           engine_funnel_L1plus="сметы уровней 1+ воронки двигателя в репозитории нет; должна уложиться "
                                                "в headroom_upper"),
               run_order=stops)
    json.dump(res, open(BUDGET, "w", encoding="utf-8"), ensure_ascii=False, indent=1, default=float)
    print(json.dumps({k: (round(v, 2) if isinstance(v, float) else v) for k, v in res["totals"].items()}, ensure_ascii=False))
    for k, c in checks.items():
        print(k, "plan код", c["exit_code"], "прогноз", c["forecast_usd"])


if __name__ == "__main__":
    main()
