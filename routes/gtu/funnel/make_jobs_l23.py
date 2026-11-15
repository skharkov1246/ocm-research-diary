"""Задачи AWS уровней 2 и 3 воронки, сводная смета уровней 0-3 против общего потолка 2000 долл. и индекс прогонов aws/RUNS.json.

  python3 make_jobs_l23.py [--l3-slots 5] [--plan-out DIR]

Пишет (AWS не вызывается; runner.py plan - без AWS):
- jobs_l2.json - уровень 2 для 50 разных слитого уровня 1: геометрия колёс и осесимметричная прочность всех колёс и дисков
  (l2-geomfea-*, контейнер ubuntu:24.04 + CalculiX + gmsh), грубый CFD центробежных ступеней в SU2 (l2-su2-*, контейнер со сборкой
  SU2 v8.6.0, как пакеты узла su2); в pending_jobs - то, для чего нет кода (CFD турбин, осевых и диагональных ступеней, камера 1D);
- jobs_l3.json - уровень 3 для 5 (--l3-slots) лучших уровня 2: полная линия частоты центробежных ступеней на сетке medium
  (l3-su2-*); в pending_jobs - сопряжённый теплообмен горячей части (OpenFOAM chtMultiRegionFoam), трёхмерная прочность, частоты и
  ползучесть (CalculiX), CFD турбин; камера сгорания (OpenFOAM) - пакет воронки камеры ../combustor/jobs_l3.json;
- funnel_budget.json - сводная смета уровней 0-3 (машино-часы × тариф спота c7a.16xlarge us-east-1 + диск + IPv4 + хранение +
  передача), фактический расход по ledger.jsonl, другие пакеты трека под тем же потолком, порядок запуска, точки остановки
  (предохранитель раннера 2000 долл., отчёт после каждых 500 - решение владельца 14.11.2026), критерии прохода уровней, что
  можно сделать локально;
- aws/RUNS.json - прогоны AWS воронки: что сделано (run_id, расход, где результаты), что ждёт и куда класть результаты.
Числа сметы - оценка; основа - в блоках estimate и в README узла.
"""
import argparse, json, math, os, sys
from collections import Counter, defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))
GTU = os.path.dirname(HERE)
REL = "routes/gtu/funnel"
sys.path.insert(0, HERE)
import make_jobs_l0 as J0          # noqa: E402

DATE = "2026-11-14"
BRANCH = "claude/lucid-allen-yvps5d"
PRICE, VCPU = J0.PRICE, J0.VCPU
EBS_H = J0.EBS_GB * J0.EBS_USD_GB_MONTH / 730
CEILING, REPORT_STEP, BUDGET_L01 = 2000.0, 500.0, 1500.0
GRACE_MIN = 10                                  # раннер: срок машины = max_min + grace_min
SU2_JOBS = os.path.join(GTU, "platform", "su2", "jobs_aws_su2_kv3.json")
_SU2 = json.load(open(SU2_JOBS, encoding="utf-8"))
TM = _SU2["estimate"]["time_model"]
C_US = dict(expected=TM["C_expected_us"], sizing=TM["C_sizing_us"])        # мкс на узел·итерацию на процесс (узел su2, замеры)
EFF = TM["EFF"]
NODES = dict(coarse=dict(expected=110240, sizing=TM["nodes"]["kv3_stage1/coarse"],
                         basis="ожидаемо - сетка coarse колеса 1-й ступени варианта idx 1164 пробы уровня 1 (impeller.py --case custom, "
                               "15.11.2026, 110 240 узлов, 37 с на 4 потоках); для запаса - KV-3 ступень 1 coarse (двусторонний вход) 258 565"),
             medium=dict(expected=TM["nodes"]["core130/medium"], sizing=TM["nodes"]["kv3_stage1/medium"],
                         basis="core130 medium 584 357 (mesh_stats.json) и KV-3 ступень 1 medium 1 137 992 (frame3_geom) - узел su2"))
MESH_MIN = dict(coarse=dict(expected=1.0, sizing=3.0), medium=dict(expected=TM["mesh_min"]["core130/medium"], sizing=TM["mesh_min"]["kv3_stage1/medium"]))
BASE_ITERS = 2000            # база - стратегия перебора разгона (plan_sweep_*): обязательных итераций до 2000 на колесо (jobs_aws_su2_kv3_sweep.json)
POINT_ITER = dict(coarse=3000, medium=6000)
POINTS_PER_JOB = 4           # точки идут одновременно по 16 процессов
BOOT_MIN, SETUP_SU2_MIN, RESERVE_MIN = TM["boot_min"], TM["setup_min"], TM["reserve_min"]
SETUP_GEOM_MIN = 8.0         # оценка: apt (ccx, библиотеки gmsh) + venv по docker/requirements.txt, без сборки SU2
N_L2 = 50
SHARDS_GEOMFEA, MACH_GEOMFEA = 10, 2
SHARDS_SU2_L2 = 17           # >= 10 частей; ~3 варианта (до 6 колёс) на часть; машины идут очередью под квоту 6
L2_CASES = f"{REL}/aws/l2/l2_cases_KV-3.json"
L3_CASES = f"{REL}/aws/l3/l3_cases_KV-3.json"
L1_TOP50 = f"{REL}/aws/l1_full/funnel_l1_top50_KV-3_m3983.json"
SWEEP_PLAN = "routes/gtu/platform/su2/plan_sweep_WINNER.json"   # заглушка: файла нет - задача упадёт сразу, пока не подставлен победитель
SU2_ENV = ("set -e\nexport PATH=\"/opt/venv/bin:/opt/su2/bin:$PATH\" SU2_RUN=/opt/su2/bin PYTHONDONTWRITEBYTECODE=1\n"
           "export OMPI_MCA_btl_vader_single_copy_mechanism=none OMPI_MCA_osc=^ucx\nW=/tmp/gtu_work; mkdir -p \"$W\"\n"
           "{ cat /opt/su2/VERSION; python -V; mpirun --version | head -n 1; echo \"nproc $(nproc) NPROC $NPROC\"; } > \"$GTU_OUT/env.txt\" 2>&1 || true\n")
GEOM_SETUP = ["apt-get update -q",
              "apt-get install -y -q --no-install-recommends python3 python3-venv calculix-ccx libgl1 libglu1-mesa libxrender1 libxcursor1 "
              "libxfixes3 libxft2 libfontconfig1 libxinerama1 libgomp1 ca-certificates",
              "python3 -m venv /opt/venv",
              "/opt/venv/bin/pip install -q --no-cache-dir -r routes/gtu/platform/docker/requirements.txt",
              "/opt/venv/bin/python -c 'import numpy, scipy, gmsh; print(\"py ok\", gmsh.__version__)'",
              "ccx -v || true"]


def _usd(mh, gb=0.0):
    u = dict(spot_expected=round(mh * PRICE["expected"], 2), spot_alt=round(mh * PRICE["alt"], 2), disk_ipv4=round(mh * (EBS_H + J0.IPV4_USD_H), 3),
             s3_3months=round(gb * J0.S3_USD_GB_MONTH * 3, 3), transfer=round(gb * J0.EGRESS_USD_GB, 3))
    u["total_expected"] = round(u["spot_alt"] + u["disk_ipv4"] + u["s3_3months"] + u["transfer"], 2)
    return u


def _bound(jobs, defaults):
    return round(sum((PRICE["upper"] + j.get("disk_gb", defaults.get("disk_gb", J0.EBS_GB)) * J0.EBS_USD_GB_MONTH / 730 + J0.IPV4_USD_H)
                     * (j.get("max_min", defaults.get("max_min", 30)) + GRACE_MIN) / 60 for j in jobs), 2)


def stage_minutes(level, which, procs_base=64, procs_point=16):
    """минуты одной ступени: сетка + база (BASE_ITERS на 64 процессах) + 4 точки одновременно (POINT_ITER на 16)"""
    t_it = lambda procs: C_US[which] * 1e-6 * NODES[level][which] / procs / EFF
    return MESH_MIN[level][which] + (BASE_ITERS * t_it(procs_base) + POINT_ITER[level] * t_it(procs_point)) / 60


def mean_centrifugal():
    """среднее число центробежных колёс на вариант - по допустимым записям отбора уровня 0 (вход полного уровня 1)"""
    d = json.load(open(os.path.join(HERE, "l0_aws_full", "funnel_l0_select_KV-3_n200000.json"), encoding="utf-8"))
    E = [e for e in d["selected"] if e["feasible"]]
    nc = {"1Ц": 1, "2Ц": 2, "Д": 0, "ОЦ": 1}
    comp = Counter(e["niche"].split("|")[2] for e in E)
    return dict(mean=round(sum(nc[k] * v for k, v in comp.items()) / len(E), 3), max=2, by_comp=dict(comp),
                non_centrifugal_share=round((comp["Д"] + comp["ОЦ"]) / len(E), 3), n=len(E),
                basis="допустимые записи отбора уровня 0 KV-3 N = 200 000 (l0_aws_full/funnel_l0_select_KV-3_n200000.json): 1Ц - 1 колесо, "
                      "2Ц - 2, ОЦ - 1 центробежное + осевые, Д - диагональное (impeller.py его не строит)")


# =========================================================================== уровень 2
def level2():
    mc = mean_centrifugal()
    # геометрия + прочность
    gf_jobs = []
    for m in range(MACH_GEOMFEA):
        mine = [k for k in range(SHARDS_GEOMFEA) if k % MACH_GEOMFEA == m]
        gf_jobs.append(dict(name=f"l2-geomfea-m{m}", container="ubuntu:24.04", container_setup=GEOM_SETUP, max_min=45,
                            command=(f"/opt/venv/bin/python {REL}/funnel_l2.py run --cases {L2_CASES} --what geom,fea --level coarse "
                                     f"--shards {SHARDS_GEOMFEA} --shard {','.join(map(str, mine))} --workers \"$NPROC\" --out \"$GTU_OUT\"\n"
                                     # раннер: bash cmd.sh без -e - код задачи = код последней строки; код расчёта - после списка частей
                                     "rc=$?\nls -l \"$GTU_OUT\" > \"$GTU_OUT/_parts_list.txt\"\nexit $rc")))
    gf_mh = MACH_GEOMFEA * (BOOT_MIN + SETUP_GEOM_MIN + 10.0) / 60          # оценка: сетки coarse параллельно (37 с на 4 потоках на колесо), диски по 4,5 с
    # SU2 coarse
    per_shard_var = math.ceil(N_L2 / SHARDS_SU2_L2)
    st_exp, st_siz = stage_minutes("coarse", "expected"), stage_minutes("coarse", "sizing")
    mm = int(min(720, math.ceil(BOOT_MIN + SETUP_SU2_MIN + per_shard_var * mc["max"] * st_siz + RESERVE_MIN)))
    su2_jobs = [dict(name=f"l2-su2-s{k:02d}", max_min=mm, disk_gb=50,
                     command=(SU2_ENV + f"python {REL}/funnel_l2.py run --cases {L2_CASES} --what su2 --level coarse --shards {SHARDS_SU2_L2} "
                                        f"--shard {k} --workers \"$NPROC\" --job-max-min {mm} --base-strategy {SWEEP_PLAN} --work \"$W\" --out \"$GTU_OUT\""))
                for k in range(SHARDS_SU2_L2)]
    n_st = N_L2 * mc["mean"]
    su2_mh = (n_st * st_exp + SHARDS_SU2_L2 * (BOOT_MIN + SETUP_SU2_MIN)) / 60
    su2_mh_siz = (n_st * st_siz + SHARDS_SU2_L2 * (BOOT_MIN + SETUP_SU2_MIN)) / 60
    defaults = dict(instance_type="c7a.16xlarge", fallback_types=["c7i.16xlarge"], image="al2023", disk_gb=J0.EBS_GB)
    su2_def = {k: v for k, v in _SU2["defaults"].items()}
    # в одном файле - один defaults: у задач SU2 контейнер и container_setup свои (как в пакетах узла su2)
    for j in su2_jobs:
        j["container"] = su2_def["container"]
        j["container_setup"] = su2_def["container_setup"]
    jobs = gf_jobs + su2_jobs
    pend_turb_mh = N_L2 * 3.0 * st_exp / 60          # оценка: 3 ступени турбин на вариант (газогенератор + силовая), как колесо
    pend_ax_mh = N_L2 * mc["non_centrifugal_share"] * 1.5 * st_exp / 60
    pending = [dict(name="l2-turb-cfd", what="грубый CFD ступеней турбины газогенератора и силовой турбины (SU2, mixing-plane)",
                    ready="НЕ готово: нет генератора геометрии турбинной решётки и сетки канала (узел geom - только центробежное колесо); "
                          "постановка - навык su2-turbo (Giles, mixing-plane)",
                    machine_hours=dict(expected=round(pend_turb_mh, 1), sizing=round(pend_turb_mh * st_siz / st_exp, 1)),
                    basis="50 вариантов × 3 ступени × время ступени coarse как у колеса компрессора (оценка)"),
               dict(name="l2-axdiag-cfd", what="CFD осевых ступеней (ОЦ) и диагональных колёс",
                    ready="НЕ готово: impeller.py строит только центробежное колесо",
                    machine_hours=dict(expected=round(pend_ax_mh, 1), sizing=round(pend_ax_mh * st_siz / st_exp, 1)),
                    basis=f"доля вариантов с ОЦ или Д {mc['non_centrifugal_share']} × 1,5 ступени (оценка)"),
               dict(name="l2-comb1d", what="одномерная камера сгорания",
                    ready="НЕ готово: стык с воронкой камеры (../combustor: merge_l1 -> feed_to_engine_funnel, поля targets.json "
                          "funnel_interface) не написан; счёт камеры - в её пакетах (смета там)",
                    machine_hours=dict(expected=0.0, sizing=0.0), basis="AWS не нужен: подбор по готовым итогам воронки камеры")]
    for p in pending:
        p["usd"] = dict(expected=round(p["machine_hours"]["expected"] * PRICE["alt"], 2), planning=round(p["machine_hours"]["sizing"] * PRICE["alt"], 2))
    gb = n_st * 0.25 + 0.05                               # оценка: сводки, history, поверхности, рестарты coarse ~0,25 ГБ на колесо
    est = dict(level="L2", variants=N_L2, centrifugal_stages=dict(mean_per_variant=mc["mean"], total_expected=round(n_st, 1), basis=mc["basis"]),
               su2_time_model=dict(formula="t_итерации = C · узлы / процессы / EFF (узел su2); ступень = сетка + база BASE_ITERS на 64 + "
                                           "POINT_ITER на 16 (4 точки одновременно)", C_us=C_US, EFF=EFF, nodes=NODES["coarse"],
                                   base_iters=BASE_ITERS, point_iter=POINT_ITER["coarse"],
                                   stage_min=dict(expected=round(st_exp, 1), sizing=round(st_siz, 1)), setup_min=SETUP_SU2_MIN, boot_min=BOOT_MIN),
               machine_hours=dict(geomfea=round(gf_mh, 2), su2_expected=round(su2_mh, 1), su2_sizing=round(su2_mh_siz, 1)),
               outputs_gb=round(gb, 1))
    est["usd"] = _usd(gf_mh + su2_mh, gb)
    est["usd_planning"] = _usd(gf_mh + su2_mh_siz, gb)["total_expected"]
    est["upper_bound_usd"] = dict(geomfea=_bound(gf_jobs, defaults), su2=_bound(su2_jobs, defaults))
    est["upper_bound_usd"]["total"] = round(sum(est["upper_bound_usd"].values()), 2)
    est["pending_usd"] = dict(expected=round(sum(p["usd"]["expected"] for p in pending), 2), planning=round(sum(p["usd"]["planning"] for p in pending), 2))
    doc = dict(_comment=("Уровень 2 воронки (README узла funnel, «Уровни 2-3»). Раннер берёт только defaults и jobs; блоки input, gates, merge, "
                         "pending_jobs, estimate, plan_check - для сессии с ключами. Код - ветка " + BRANCH + " (PR #62). Сгенерировано "
                         "make_jobs_l23.py; правки - в нём. Задачи l2-su2-* запускать только после ворот gates.su2."),
               defaults=defaults, jobs=jobs,
               input=dict(top50=L1_TOP50, cases=L2_CASES,
                          make_cases=f"python3 {REL}/funnel_l2.py prepare --input {L1_TOP50} --top {N_L2} --out {L2_CASES}",
                          note="файлы входа появятся после слияния полного уровня 1 (RUNBOOK_L0_L1.md, шаг 7); l2_cases - в ветку до run"),
               gates=dict(geomfea="слит полный уровень 1, l2_cases_KV-3.json в ветке", su2=(
                   "перебор стратегий разгона SU2 (platform/su2/jobs_aws_su2_kv3_sweep.json) дал вердикт success на KV-3 coarse; имя "
                   "плана-победителя подставить вместо " + SWEEP_PLAN + " и пересобрать (make_jobs_l23.py --sweep-plan); "
                   "до этого база SU2 на сетках geom расходится (README узла su2) - деньги на 17 машин без вердикта не тратить"),
                   order="сначала l2-geomfea-* (дёшево, отсев по прочности), затем l2-su2-* только для прошедших прочность (prepare по ним)"),
               merge=dict(commands=[f"python3 {REL}/funnel_l2.py merge --cases {L2_CASES} --in routes/gtu/platform/aws/runs/<run_id> "
                                    f"--what geom,fea --shards {SHARDS_GEOMFEA} --out-dir {REL}/aws/l2",
                                    f"python3 {REL}/funnel_l2.py merge --cases {L2_CASES} --in routes/gtu/platform/aws/runs/<run_id> "
                                    f"--what su2 --shards {SHARDS_SU2_L2} --out-dir {REL}/aws/l2/su2"],
                          where="сессия с ключами после run; в ветку - funnel_l2_results.json и части (без сеток и рестартов)"),
               pending_jobs=pending, estimate=est)
    return doc


# =========================================================================== уровень 3
def level3(slots):
    mc = mean_centrifugal()
    pf = [0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2]
    halves = dict(a=pf[0::2], b=pf[1::2])
    st_exp, st_siz = stage_minutes("medium", "expected"), stage_minutes("medium", "sizing")
    mm = int(min(720, math.ceil(BOOT_MIN + SETUP_SU2_MIN + st_siz + RESERVE_MIN)))
    jobs = []
    for s in range(slots):
        for stg in (1, 2):
            for h, fr in halves.items():
                jobs.append(dict(name=f"l3-su2-v{s}-st{stg}{h}", max_min=mm, disk_gb=60, container=_SU2["defaults"]["container"],
                                 container_setup=_SU2["defaults"]["container_setup"],
                                 command=(SU2_ENV + f"python {REL}/funnel_l2.py run --cases {L3_CASES} --what su2 --level medium --shards {slots} "
                                                    f"--shard {s} --stage {stg} --pfracs {','.join(map(str, fr))} --tag {h} --workers \"$NPROC\" "
                                                    f"--job-max-min {mm} --base-strategy {SWEEP_PLAN} --work \"$W\" --out \"$GTU_OUT\"")))
    n_full = slots * mc["mean"] * 2                       # задач с колесом (у варианта с одним колесом задачи st2 пусты)
    n_empty = len(jobs) - n_full
    mh = (n_full * (BOOT_MIN + SETUP_SU2_MIN + st_exp) + n_empty * (BOOT_MIN + SETUP_SU2_MIN)) / 60
    mh_siz = (n_full * (BOOT_MIN + SETUP_SU2_MIN + st_siz) + n_empty * (BOOT_MIN + SETUP_SU2_MIN)) / 60
    defaults = dict(instance_type="c7a.16xlarge", fallback_types=["c7i.16xlarge"], image="al2023", disk_gb=60)
    turb_mh = slots * 3.0 * (BOOT_MIN + SETUP_SU2_MIN + st_exp) / 60      # турбина: одна задача (4 точки) на ступень - ветвь запирания короче
    ccx_mh, cht_mh = slots * 2.0, slots * 1.7
    pending = [dict(name="l3-turb-cfd", what="характеристики ступеней турбины на сетке medium (SU2)",
                    ready="НЕ готово: нет геометрии и сетки турбинной решётки (как l2-turb-cfd)",
                    machine_hours=dict(expected=round(turb_mh, 1), sizing=round(turb_mh * st_siz / st_exp, 1)),
                    basis="3 ступени × одна задача (4 точки по 16 процессов) × время задачи medium, как у колеса (оценка)"),
               dict(name="l3-ccx3d", what="трёхмерная прочность, частоты с преднагружением, диаграмма Кэмпбелла и ползучесть колёс и дисков "
                                          "(CalculiX, сектор из STEP impeller.py через Gmsh, циклическая симметрия; навык calculix-rotor)",
                    ready="НЕ готово: кода узла 3D-прочности нет (platform/fea - осесимметричный disk.py; platform/bench/speed_ccx.py - "
                          "только замер скорости ccx на C3D10); STEP колёс есть (impeller.py), турбинных лопаток - нет",
                    machine_hours=dict(expected=round(ccx_mh, 1), sizing=round(ccx_mh * 2, 1)),
                    basis="2 машино-часа на вариант: ~4 колеса × (статика 3 частоты + частоты 10 форм на 5 частотах + ползучесть), "
                          "сектор ~450 тыс. узлов C3D10 (оценка; замера ccx на AWS нет)"),
               dict(name="l3-cht", what="сопряжённый теплообмен печатной/литой горячей части (сопловой аппарат и корпус с каналами "
                                        "охлаждения; OpenFOAM chtMultiRegionFoam, навык openfoam-combustor)",
                    ready="НЕ готово: нет геометрии каналов охлаждения горячей части и словарей областей (у воронки камеры тот же "
                          "пробел - comb-l2-cht); свойства печатных сплавов k(T), cp(T) не собраны",
                    machine_hours=dict(expected=round(cht_mh, 1), sizing=round(cht_mh * 3, 1)),
                    basis="3 млн ячеек × 5000 итераций × 2e-5 ядро-с на ячейку·итерацию / (64 × 0,75) ≈ 1,7 ч на вариант (оценка, замера нет)"),
               dict(name="l3-comb", what="камера сгорания (OpenFOAM: RANS/LES элемента, NOx и CO, температура стенки)",
                    ready="пакет воронки камеры ../combustor/jobs_l3.json (comb_l23.py --level L3) по входу l3_inputs_KV-3.json; "
                          "связь с вариантами двигателя - по типу камеры и p3/T3/T41 (стык не написан); смета - в воронке камеры",
                    machine_hours=dict(expected=0.0, sizing=0.0), basis="считается в смете воронки камеры (comb_budget_results.json)")]
    for p in pending:
        p["usd"] = dict(expected=round(p["machine_hours"]["expected"] * PRICE["alt"], 2), planning=round(p["machine_hours"]["sizing"] * PRICE["alt"], 2))
    gb = n_full * 0.5 * 1.1                               # оценка: рестарты medium (104 Б/узел) и сводки на задачу ~0,5 ГБ
    est = dict(level="L3", slots=slots, su2=dict(level="medium", pfracs=pf, halves=halves, stage_min=dict(expected=round(st_exp, 1), sizing=round(st_siz, 1)),
                                                nodes=NODES["medium"], point_iter=POINT_ITER["medium"], base_iters=BASE_ITERS, jobs=len(jobs),
                                                jobs_with_wheel_expected=round(n_full, 1)),
               machine_hours=dict(expected=round(mh, 1), sizing=round(mh_siz, 1)), outputs_gb=round(gb, 1))
    est["usd"] = _usd(mh, gb)
    est["usd_planning"] = _usd(mh_siz, gb)["total_expected"]
    est["upper_bound_usd"] = dict(total=_bound(jobs, defaults))
    est["pending_usd"] = dict(expected=round(sum(p["usd"]["expected"] for p in pending), 2), planning=round(sum(p["usd"]["planning"] for p in pending), 2))
    est["extension_10"] = dict(note="расширение до 10 вариантов - make_jobs_l23.py --l3-slots 10 (задач вдвое больше, смета ~×2)",
                               usd_expected=round(est["usd"]["total_expected"] * 10 / slots, 2),
                               upper_bound_usd=round(est["upper_bound_usd"]["total"] * 10 / slots, 2))
    doc = dict(_comment=("Уровень 3 воронки (README узла funnel, «Уровни 2-3»). Раннер берёт только defaults и jobs. Задачи - по одной ступени и "
                         "половине точек (a/b) на машину: 4 точки по 16 процессов, сетка medium; вариант с одним колесом - задачи st2 "
                         "завершаются сразу. Сгенерировано make_jobs_l23.py."),
               defaults=defaults, jobs=jobs,
               input=dict(cases=L3_CASES, make_cases=(f"python3 {REL}/funnel_l2.py prepare --input {L1_TOP50} --select {REL}/aws/l2/funnel_l2_results.json "
                                                      f"--top {slots} --out {L3_CASES}"),
                          note="l3_select - прошедшие уровень 2 (не больше 2 на архетип, до 10); меньше 5 - точка остановки"),
               gates=dict(su2="уровень 2 SU2 прошёл (линия coarse получена у >= 3 из 4 точек и расхождение с одномерной моделью в пределах); "
                              "сверка постановки на открытом эталоне (Rotor 37 или CC3, узел su2) - до того числа medium - не характеристика"),
               merge=dict(command=f"python3 {REL}/funnel_l2.py merge --cases {L3_CASES} --in routes/gtu/platform/aws/runs/<run_id> --what su2 "
                                  f"--shards {slots} --out-dir {REL}/aws/l3"),
               pending_jobs=pending, estimate=est)
    return doc


# =========================================================================== сводная смета
def ledger_runs():
    rows = [json.loads(x) for x in open(os.path.join(GTU, "platform", "aws", "ledger.jsonl"), encoding="utf-8")]
    launch = {r["iid"]: r for r in rows if r["ev"] == "launch"}
    by = defaultdict(lambda: dict(usd=0.0, machine_min=0.0, jobs=[]))
    for r in rows:
        if r["ev"] != "end":
            continue
        l = launch[r["iid"]]
        dur = max(60.0, r["t_end"] - l["t_launch"])
        b = by[r["run_id"]]
        b["usd"] += (l["spot_usd_h"] + l["ebs_usd_h"] + l["ip_usd_h"]) * dur / 3600
        b["machine_min"] += dur / 60
        b["jobs"].append(r["job"])
    return {k: dict(usd=round(v["usd"], 3), machine_min=round(v["machine_min"], 1), jobs=sorted(v["jobs"])) for k, v in by.items()}


FUNNEL_RUNS = {   # прогоны воронки двигателя в журнале (README узла, разделы уровня 0)
    "20261114T231900Z-b7cc": dict(level="L0", what="отладка, выборка v1 (заменена)", results="снято (выборка v1 заменена v2)"),
    "20261114T232200Z-7ce3": dict(level="L0", what="полный, выборка v1 (заменена)", results="S3: <префикс ALPHA_O_S3>/gtu/funnel/l0_full_20261114T232200Z-7ce3/"),
    "20261115T010201Z-1a6e": dict(level="L0", what="отладка, выборка v2, N = 2000 × 3 точки", results=f"{REL}/l0_aws_debug/"),
    "20261115T010625Z-8893": dict(level="L0", what="полный, выборка v2: KV-3 200 000, KV-1 и KV-6 по 50 000", results=(
        f"{REL}/l0_aws_full/ (JSON); CSV - S3 <префикс ALPHA_O_S3>/gtu/funnel/l0v2_full_20261115T010625Z-8893/")),
}


def pkg(path):
    d = json.load(open(os.path.join(REPO, path), encoding="utf-8"))
    e = d.get("estimate") or {}
    u = e.get("usd") or {}
    exp = u.get("total_expected") if isinstance(u, dict) else None
    up = e.get("upper_bound_usd")
    if isinstance(up, dict):
        up = up.get("total")
    up = up if up is not None else (e.get("usd_upper_bound") or (d.get("plan_check") or {}).get("forecast_usd"))
    return dict(file=path, jobs=len(d["jobs"]), expected=exp, upper=up, planning=e.get("usd_planning"),
                pending_expected=(e.get("pending_usd") or {}).get("expected", 0.0), pending_planning=(e.get("pending_usd") or {}).get("planning", 0.0))


def budget(l2, l3, plans):
    runs = ledger_runs()
    spent = round(sum(v["usd"] for v in runs.values()), 3)
    funnel_spent = round(sum(runs.get(k, {}).get("usd", 0.0) for k in FUNNEL_RUNS), 3)
    rows = []
    for lev, f, status in [("L0", f"{REL}/jobs_l0_debug.json", "сделан: 20261115T010201Z-1a6e"), ("L0", f"{REL}/jobs_l0.json", "сделан: 20261115T010625Z-8893"),
                           ("L1", f"{REL}/jobs_l1_debug.json", "ждёт"), ("L1", f"{REL}/jobs_l1.json", "ждёт"), ("L1", f"{REL}/jobs_l1_fea.json", "ждёт входа"),
                           ("L2", f"{REL}/jobs_l2.json", "ждёт уровня 1 и ворот SU2"), ("L3", f"{REL}/jobs_l3.json", "ждёт уровня 2")]:
        p = pkg(f)
        p.update(level=lev, status=status, plan_exit_code=plans.get(os.path.basename(f)))
        if lev == "L0":
            rid = status.split(": ")[1]
            p["actual_usd"] = runs.get(rid, {}).get("usd")
            p["actual_machine_min"] = runs.get(rid, {}).get("machine_min")
        rows.append(p)
    l1fea = next(r for r in rows if r["file"].endswith("jobs_l1_fea.json"))
    l1fea["expected"] = l1fea["expected"] or round(0.5 * PRICE["alt"], 2)        # оценка: машина ~30 мин (300 вариантов × до 3 колёс по 5-20 с)
    eng = dict(expected=round(sum((r["actual_usd"] if r.get("actual_usd") is not None else (r["expected"] or 0)) for r in rows), 2),
               planning=round(sum((r["actual_usd"] if r.get("actual_usd") is not None else (r.get("planning") or r["expected"] or 0)) + r["pending_planning"] for r in rows), 2),
               upper=round(sum((r["actual_usd"] if r.get("actual_usd") is not None else (r["upper"] or 0)) + r["pending_planning"] for r in rows), 2),
               pending_planning=round(sum(r["pending_planning"] for r in rows), 2))
    l01_upper = round(sum((r["actual_usd"] if r.get("actual_usd") is not None else (r["upper"] or 0)) for r in rows if r["level"] in ("L0", "L1")), 2)
    comb = json.load(open(os.path.join(GTU, "combustor", "comb_budget_results.json"), encoding="utf-8"))["totals"]
    others = [dict(name="воронка камеры: калибровка + уровни 1-3 + каналы (comb_budget_results.json)", expected=round(comb["combustor_expected"], 2),
                   planning=round(comb["combustor_planning"], 2), upper=round(comb["combustor_upper"], 2))]
    for name, f in [("камера L0, N = 262 144 (сделан 20261115T001745Z-f1da)", "routes/gtu/combustor/jobs_aws_comb_l0.json"),
                    ("SU2 KV-3: перебор стратегий разгона", "routes/gtu/platform/su2/jobs_aws_su2_kv3_sweep.json"),
                    ("SU2 KV-3: medium", "routes/gtu/platform/su2/jobs_aws_su2_kv3.json"),
                    ("SU2 KV-3: fine", "routes/gtu/platform/su2/jobs_aws_su2_kv3_fine.json"),
                    ("органические циклы (сделан 20261115T001745Z-f55e)", "routes/gtu/organic/jobs_aws_organic.json")]:
        p = pkg(f)
        done = "сделан" in name
        rid = name.split("сделан ")[1].rstrip(")") if done else None
        others.append(dict(name=name, file=f, expected=(runs.get(rid, {}).get("usd") if done else p["expected"]),
                           planning=(runs.get(rid, {}).get("usd") if done else p["expected"]),
                           upper=(runs.get(rid, {}).get("usd") if done else p["upper"])))
    oth = {k: round(sum((o[k] or 0) for o in others), 2) for k in ("expected", "planning", "upper")}
    other_spent = round(spent - funnel_spent - sum(runs.get(r, {}).get("usd", 0) for r in ("20261115T001745Z-f1da", "20261115T001745Z-f55e")), 3)
    total = {k: round(eng[k] + oth[k] + other_spent, 2) for k in ("expected", "planning", "upper")}
    # порядок запуска и пересечение порогов отчёта (накопленно, верхняя граница и ожидаемо)
    order = [("L1 отладка", "jobs_l1_debug.json"), ("L1 полный + l1-starts2x", "jobs_l1.json"), ("L1 CalculiX лучших 300", "jobs_l1_fea.json"),
             ("L2 геометрия + прочность", "l2.geomfea"), ("L2 SU2 coarse (после ворот)", "l2.su2"), ("L3 SU2 medium, 5 вариантов", "jobs_l3.json")]
    cum_e, cum_u, steps = spent, spent, []
    l2e = l2["estimate"]
    vals = {"jobs_l1_debug.json": (pkg(f"{REL}/jobs_l1_debug.json")["expected"], pkg(f"{REL}/jobs_l1_debug.json")["upper"]),
            "jobs_l1.json": (pkg(f"{REL}/jobs_l1.json")["expected"], pkg(f"{REL}/jobs_l1.json")["upper"]),
            "jobs_l1_fea.json": (l1fea["expected"], l1fea["upper"]),
            "l2.geomfea": (round(l2e["machine_hours"]["geomfea"] * PRICE["alt"], 2), l2e["upper_bound_usd"]["geomfea"]),
            "l2.su2": (round(l2e["usd"]["total_expected"] - l2e["machine_hours"]["geomfea"] * PRICE["alt"], 2), l2e["upper_bound_usd"]["su2"]),
            "jobs_l3.json": (l3["estimate"]["usd"]["total_expected"], l3["estimate"]["upper_bound_usd"]["total"])}
    for name, key in order:
        e, u = vals[key]
        cum_e += e or 0
        cum_u += u or 0
        steps.append(dict(step=name, expected=e, upper=u, cum_expected=round(cum_e, 2), cum_upper=round(cum_u, 2),
                          report_threshold_crossed=[t for t in (500, 1000, 1500, 2000) if cum_u - (u or 0) < t <= cum_u]))
    return dict(spent_ledger=spent, spent_funnel_runs=funnel_spent, runs=runs, engine_rows=rows, engine=eng, level01_upper=l01_upper,
                level01_budget=BUDGET_L01, level01_ok=l01_upper < BUDGET_L01, others=others, others_sum=oth, other_spent_outside_lists=other_spent,
                total=total, ceiling=CEILING, report_step=REPORT_STEP, headroom_upper=round(CEILING - total["upper"], 2),
                engine_launch_order=steps)


CRITERIA = {
    "L0": ["полнота: все части и индексы 0..N-1 ровно по разу, у всех частей один отпечаток выборки и кода (merge_l0 - иначе код 2)",
           "доля допустимых по точке отличается от отладки N = 2000 не больше чем вдвое (тот же код и выборка; KV-3 v2: 54,3 %, 393 ниши; "
           "KV-1 и KV-6 - по их results)",
           "сходимость: новых допустимых ниш при удвоении N - не больше 1 % (KV-3: +1 из 393 - сошлось); топ-50 по чистой оценке не обязан "
           "сходиться (27 новых при удвоении при приросте средней оценки +0,025 - плато): уровень 1 оптимизирует непрерывные оси внутри ниш",
           "отбор на уровень 1: каждая допустимая ниша хотя бы одним вариантом, 3000-5000 записей (KV-3: 3983)"],
    "L1": ["полнота: 3983 из 3983 вариантов, части 0..39 ровно по разу, один отпечаток кода и входа (merge_l1 - иначе код 2)",
           "время на вариант на AWS не больше 2 × сметного (307 с на ядро при 4 стартах) - иначе остановка и пересчёт сметы",
           "допустимых с оценкой >= 5,0 - не меньше 50 (иначе итог - меньше 50 вариантов, причины - reject_detail; это результат моделей)",
           "сходимость по стартам (l1-starts2x, 640 вариантов, 4 против 8 стартов): у >= 90 % вариантов оценка меняется меньше чем на 0,05 и "
           "допустимость не меняется; иначе - повтор с 8 стартами",
           "сходимость по объёму (таблица convergence, M/2 против M): общих ниш в отборе 50 >= 80 % и средняя оценка 50 разных меняется < 0,05; "
           "иначе - merge_l0 --target 10000 и повтор (смета ×2,5)",
           "цена разнообразия: средняя оценка 50 разных против чистого топ-50 и число вытесненных - записать (критерия нет: это мера)"],
    "L2": ["полнота: все колёса и диски 50 вариантов ровно по разу (funnel_l2 merge, иначе код 2)",
           "геометрия: impeller.py построил каждое центробежное колесо (код 0); предупреждения impeller.py просмотрены",
           f"прочность (disk.py, сетка 80×12): запас по разрушающей частоте при 120 % >= 1,0 у всех колёс и дисков; ресурс по ползучести >= 25 000 ч",
           "SU2 coarse: линия частоты получена у >= 3 из 4 точек каждой ступени; в расчётной точке |Δη| <= 2 п. п. и |Δπ|/π <= 5 % против "
           "одномерной модели уровня 1 - иначе КПД варианта пересчитывается по SU2 и оценка обновляется",
           "на уровень 3 - прошедшие по убыванию оценки, не больше 2 на архетип; меньше 5 прошедших - стоп и разбор (уровень 3 не запускается)"],
    "L3": ["SU2 medium: запас по помпажу по линии частоты >= 10 % (граница - последняя устойчивая точка), КПД в расчётной точке отличается "
           "от coarse не больше 1 п. п. (иначе сетка не сошлась - fine для одного колеса)",
           "3D-прочность (когда будет код): запас по разрушению при 120 % >= 1,0; резонансы с порядками 1-4 и числом лопаток соплового аппарата "
           "вне ±10 % рабочих частот 90-105 % на диаграмме Кэмпбелла; ползучесть 25 000 ч",
           "сопряжённый теплообмен (когда будет код): температура металла не выше предела класса охлаждения и сплава (make_model.t41_limit, "
           "dT_cap печатных сплавов)",
           "камера (воронка камеры): её критерии - рост мод <= 0, пульсации <= 1 % p3, стенка <= 1173 K, ресурс 30 000 ч",
           "итог - 1-2 варианта на стенд"],
}
LOCAL = ["тесты: python3 -m pytest -q test_models.py test_funnel_l0.py test_funnel_l1.py test_funnel_l2.py (минуты)",
         "проба уровня 1: python3 merge_l1.py --probe (12 вариантов, 2 старта, ~15 мин на 4 ядрах)",
         "покрытие осей уровня 0: python3 coverage_l0.py (closure по всем шаблонам трёх точек, ~2 мин)",
         "уровень 2 на пробе: funnel_l2.py prepare по funnel_l1_top50.json, run --what geom,fea --no-mesh (STEP и disk.py, ~2 мин на вариант), "
         "run --what su2 --dry (команды); одна сетка coarse impeller.py - 37 с на 4 потоках",
         "проверка входов: funnel_l1.load_select для отборов l0_aws_debug и l0_aws_full (отпечатки выборки и кода)",
         "plan всех файлов задач: make_jobs_l1.py и make_jobs_l23.py (без AWS)",
         "нельзя локально: SU2 (сборки нет, 4 ядра), объёмы уровней 0-1 больше N 300 и 12 вариантов (решение владельца), 3D CalculiX и OpenFOAM"]


def runs_index(b, l2, l3):
    r = b["runs"]
    done = [dict(run_id=k, level=v["level"], what=v["what"], results=v["results"], usd=r.get(k, {}).get("usd"),
                 machine_min=r.get(k, {}).get("machine_min"), jobs=r.get(k, {}).get("jobs")) for k, v in FUNNEL_RUNS.items()]
    todo = [dict(level="L1", file=f"{REL}/jobs_l1_debug.json", results_dir=f"{REL}/aws/l1_debug/", s3="<префикс ALPHA_O_S3>/gtu/<run_id>/l1-debug-m*/"),
            dict(level="L1", file=f"{REL}/jobs_l1.json", results_dir=f"{REL}/aws/l1_full/", s3="<префикс ALPHA_O_S3>/gtu/<run_id>/l1-full-m*/, l1-starts2x/"),
            dict(level="L1", file=f"{REL}/jobs_l1_fea.json", results_dir=f"{REL}/aws/l1_fea/", s3="<префикс ALPHA_O_S3>/gtu/<run_id>/l1-fea300/"),
            dict(level="L2", file=f"{REL}/jobs_l2.json", results_dir=f"{REL}/aws/l2/ (su2 - aws/l2/su2/)", s3="<префикс ALPHA_O_S3>/gtu/<run_id>/l2-*/"),
            dict(level="L3", file=f"{REL}/jobs_l3.json", results_dir=f"{REL}/aws/l3/", s3="<префикс ALPHA_O_S3>/gtu/<run_id>/l3-*/")]
    return dict(title="воронка двигателя-заменителя: прогоны AWS", date=DATE, ledger="routes/gtu/platform/aws/ledger.jsonl",
                rule=("в ветку - слитые results-JSON, top50 и части уровня 1 (JSON .gz до 5 МБ); CSV больше 5 МБ, сетки, рестарты SU2 - в S3 "
                      "(<префикс ALPHA_O_S3>/gtu/<run_id>/<задача>/) и runs/<run_id>/ (в .gitignore); после каждого боевого прогона - "
                      "закоммитить ledger.jsonl и дописать run_id сюда (make_jobs_l23.py, словарь FUNNEL_RUNS)"),
                done=done, todo=todo, how_to_regenerate=f"python3 {REL}/make_jobs_l23.py")


def main(argv=None):
    global SWEEP_PLAN
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--l3-slots", type=int, default=5)
    ap.add_argument("--sweep-plan", default=None, help="план-победитель перебора разгона SU2 (путь в репозитории)")
    ap.add_argument("--plan-out", default=os.environ.get("TMPDIR", "/tmp"))
    a = ap.parse_args(argv)
    if a.sweep_plan:
        SWEEP_PLAN = a.sweep_plan
    if not 5 <= a.l3_slots <= 10:
        raise SystemExit("уровень 3 - от 5 до 10 вариантов")
    docs = {"jobs_l2.json": level2(), "jobs_l3.json": level3(a.l3_slots)}
    plans = {}
    for f in ("jobs_l0_debug.json", "jobs_l0.json", "jobs_l1_debug.json", "jobs_l1.json", "jobs_l1_fea.json"):
        plans[f] = J0.plan_check(os.path.join(HERE, f), a.plan_out)["exit_code"]
    for f, doc in docs.items():
        path = os.path.join(HERE, f)
        with open(path, "w", encoding="utf-8") as fh:
            json.dump(doc, fh, ensure_ascii=False, indent=1)
        doc["plan_check"] = J0.plan_check(path, a.plan_out)
        plans[f] = doc["plan_check"]["exit_code"]
        with open(path, "w", encoding="utf-8") as fh:
            json.dump(doc, fh, ensure_ascii=False, indent=1)
        e = doc["estimate"]
        print(f"{f}: {len(doc['jobs'])} задач, ожидаемо {e['usd']['total_expected']} долл., верхняя граница {e['upper_bound_usd']['total']} долл., "
              f"ожидают кода {len(doc['pending_jobs'])} (план. {e['pending_usd']['planning']} долл.); plan: код {plans[f]}")
    b = budget(docs["jobs_l2.json"], docs["jobs_l3.json"], plans)
    out = dict(title="Воронка двигателя-заменителя: сводная смета уровней 0-3 против потолка 2000 долл.", date=DATE, status="оценка",
               method=("машино-часы × тариф спота c7a.16xlarge us-east-1 (ожидаемо 0,748 долл./ч - us-east-1b, запас к 0,706 us-east-1c; "
                       "верхняя граница - MaxPrice раннера 1,4206 × (max_min + 10 мин) на машину) + диск gp3 и IPv4 (ставки раннера) + "
                       "хранение S3 3 мес. и передача (0,023 долл./ГБ·мес, 0,09 долл./ГБ - допущения); сделанные прогоны - по ledger.jsonl "
                       "(цена спота фактической зоны × время машины); задачи без кода (pending) - по планируемой сумме"),
               software=f"python {sys.version.split()[0]}; runner.py plan (раннер v2 в рабочей копии)",
               inputs=[dict(file=f"{REL}/{f}") for f in ("jobs_l0_debug.json", "jobs_l0.json", "jobs_l1_debug.json", "jobs_l1.json", "jobs_l1_fea.json",
                                                          "jobs_l2.json", "jobs_l3.json")]
               + [dict(file="routes/gtu/platform/aws/ledger.jsonl"), dict(file="routes/gtu/combustor/comb_budget_results.json"),
                  dict(file="routes/gtu/platform/su2/jobs_aws_su2_kv3.json", field="estimate.time_model")],
               price=PRICE, plan_exit_codes=plans, budget=b, criteria=CRITERIA, local_before_aws=LOCAL,
               stop_points=["предохранитель раннера: итого (израсходовано + в работе + новая) <= 2000 долл. перед каждым запуском - жёсткий стоп",
                            "отчёт владельцу после каждых 500 долл. общего расхода (report_step раннера; решение владельца 14.11.2026) - "
                            "продолжение после ответа",
                            "уровни 0 + 1 - до 1500 долл. (решение владельца); по смете - доли процента, объём ограничен пользой",
                            "после отладки уровня 1: фактическое время на вариант > 2 × сметы - стоп, пересчёт make_jobs_l1.py",
                            "после двух первых задач l2-su2: нет линии частоты (база разошлась) - стоп SU2 уровня 2 (остальные 15 машин не запускать)",
                            "уровень 2 прошли меньше 5 - уровень 3 не запускается, разбор",
                            "расширение уровня 3 до 10 вариантов - только если после уровня 3 на 5 запас до 2000 по верхней границе больше её удвоения"])
    with open(os.path.join(HERE, "funnel_budget.json"), "w", encoding="utf-8") as fh:
        json.dump(out, fh, ensure_ascii=False, indent=1)
    os.makedirs(os.path.join(HERE, "aws"), exist_ok=True)
    with open(os.path.join(HERE, "aws", "RUNS.json"), "w", encoding="utf-8") as fh:
        json.dump(runs_index(b, docs["jobs_l2.json"], docs["jobs_l3.json"]), fh, ensure_ascii=False, indent=1)
    t = b["total"]
    print(f"funnel_budget.json: воронка ожидаемо {b['engine']['expected']}, план. {b['engine']['planning']}, верхняя граница {b['engine']['upper']} долл.; "
          f"всё вместе (с камерой и SU2 KV-3) по верхней границе {t['upper']} из {CEILING} (запас {b['headroom_upper']}); израсходовано {b['spent_ledger']}")
    return 0 if all(v == 0 for v in plans.values()) else 1


if __name__ == "__main__":
    sys.exit(main())
