#!/usr/bin/env python3
"""Пакет routes/gtu/platform/aws/runner.py: перебор шести стратегий разгона базовой точки SU2 — по стратегии на машину
c7a.16xlarge, на каждой оба колеса (KV-3 ступень 1 и core130), сетка coarse, ~20–30 мин на машину.

Пишет рядом:
  plan_sweep_v1_control.json … plan_sweep_v6_riemann.json — стратегии для sweep_ramp.py run (шаги, гейты, ключи SU2);
  jobs_aws_su2_kv3_sweep.json — 6 задач раннера + разбор провала пробы 20261115T003035Z-c8f9 (diagnosis), таблица
                                стратегий, критерии успеха, смета, сбор (merge) и вывод runner.py plan (plan_check).
Команда задачи: обе сетки параллельно (frame3_geom.py, impeller.py), затем sweep_ramp.py run для двух колёс
параллельно по NPROC/2 процессов MPI, затем sweep_ramp.py merge по $GTU_OUT.

  python make_jobs_su2_kv3_sweep.py                 # стратегии, пакет, plan_check
  python make_jobs_su2_kv3_sweep.py --no-plan       # без runner.py plan
  python make_jobs_su2_kv3_sweep.py --mini DIR      # миниатюры команд (3 итерации на шаг, 2 процесса на колесо)
"""
from __future__ import annotations

import argparse
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in sys.path:
    sys.path.insert(0, HERE)
sys.dont_write_bytecode = True
import make_jobs_su2_kv3 as M  # noqa: E402
import run_case as RC  # noqa: E402
import sweep_ramp as SW  # noqa: E402

OUT_JOBS = os.path.join(HERE, "jobs_aws_su2_kv3_sweep.json")
SMOKE_DIR = os.path.join(HERE, "aws_smoke_20261115T003035Z-c8f9")
REL = "routes/gtu/platform/su2"
MAX_MIN = 45            # сторож машины; команда должна кончиться к max_min − 8 мин (выгрузка раннера)
TIME_LIMIT_S = 1320     # 22 мин на стратегию у каждого колеса (оба колеса — параллельно); успех — раньше
BASE_LOAD = 0.75        # база — как char_points.py --base-pfrac 0.75 (точки линии частоты стартуют от неё)
LIGHT_LOAD = 0.60       # «лёгкий» выход: при 0,5 одномерная оценка даёт радиальный M на r = 1,3 r2 до 1,07 (без потерь)
DISK_GB = 40

# ----------------------------------------------------------------------------------------------- стратегии
OUT_SET = {"OUTPUT_FILES": "RESTART", "WRT_FREQ": "100"}   # рестарт каждые 100 итераций: поле <= 100 итераций до NaN
WANG = {"LIMITER": "VENKATAKRISHNAN_WANG", "VENKAT": "0.05"}


def _st(name, frac, load, iter_, ramp=0, pramp=0, muscl=False, muscl_ramp=0, gate="intermediate", hold=None, **kw):
    s = {"name": name, "frac": frac, "load": load, "iter": iter_, "ramp": ramp, "pramp": pramp, "muscl": muscl,
         "gate": gate, "hold": hold or {"iter": 300, "max": 3}}
    if muscl_ramp:
        s["muscl_ramp"] = muscl_ramp
    s.update(kw)
    return s


FINAL_HOLD = {"iter": 300, "max": 40}

STRATEGIES = [
    {
        "id": "v1_control", "axis": "контроль",
        "title": "Контроль: план smoke075 без сжатия итераций (путь прогона Б целиком)",
        "mechanism": ("Рампы частоты и давления идут одновременно, давление отстаёт (рампа вдвое длиннее): 0 -> 0,6 "
                      "частоты за 80 итераций, давление до 0,95 оценки при 0,6 за 160 (всего 170), 0,6 -> 1,0 за 80 и "
                      "давление до 0,75 оценки за 160 (200), затем 150 итераций MUSCL; переходы — по числу итераций, как "
                      "в пробе, но масштаб 1,0 вместо 0,1. Дальше — порции по 300 итераций при базе до критериев."),
        "hypothesis": ("Провал пробы вызван только сжатием рамп в 10 раз (--base-iter-scale 0.1): core130 по этому "
                       "пути в прогоне Б (README) дошёл до дисбаланса 1,4 %; KV-3 по нему не считался ни разу."),
        "differs_from_smoke": "масштаб итераций 1,0 вместо 0,1; удержание при базе с гейтом",
        "transforms": ["hist"], "set": dict(OUT_SET),
        "steps": [{"name": "smoke075", "kind": "plan", "plan": "smoke075", "load": BASE_LOAD, "gate": "final",
                   "hold": {"iter": 300, "max": 40, "muscl": True}}],
    },
    {
        "id": "v2_decouple", "axis": "путь нагружения",
        "title": "Сначала частота при лёгком выходе, затем второй порядок, затем нагружение — переходы по гейту",
        "mechanism": ("a: частота 0 -> 1,0 за 450 итераций при давлении на выходе, растущем вдвое дольше (900) только до "
                      "0,60 оценки (лёгкий режим, далеко от помпажа), первый порядок, CFL 5 (1–20); гейт. b: второй "
                      "порядок плавно — RAMP_MUSCL косинусом за 250 итераций, ограничитель VENKATAKRISHNAN_WANG "
                      "(K = 0,05 от размаха); гейт. c: давление 0,60 -> 0,75 оценки за 400 итераций (порция 600); итоговый гейт."),
        "hypothesis": ("Расходится, когда давление на выходе поднимают раньше, чем установился расход (в пробе к концу "
                       "s1 расход на выходе уже обратный); частоту и нагрузку надо развести."),
        "differs_from_smoke": "частота и давление разведены, лёгкий старт 0,60, рампы 450/900/400, плавный MUSCL, гейты",
        "transforms": ["hist"], "set": dict(OUT_SET, **WANG),
        "steps": [_st("a_spin_light", 1.0, LIGHT_LOAD, 1000, ramp=450, pramp=900, cfl=5.0),
                  _st("b_order2", 1.0, LIGHT_LOAD, 400, muscl=True, muscl_ramp=250),
                  _st("c_base", 1.0, BASE_LOAD, 600, pramp=400, muscl=True, gate="final", hold=FINAL_HOLD)],
    },
    {
        "id": "v3_ladder", "axis": "обратная связь",
        "title": "Лестница по линии дросселя с гейтом на каждой ступени и откатом при NaN",
        "mechanism": ("Частота ступенями 0,30 -> 0,55 -> 0,75 -> 0,90 -> 1,00, давление на выходе на каждой ступени — по "
                      "линии дросселя p = p_нач + w^2 (0,60 оценки − p_нач), рампа ступени 120 итераций, порция 220, гейт "
                      "и до трёх порций по 200 на установление. NaN или непройденный гейт -> откат к последнему "
                      "пройденному рестарту: та же ступень с рампой вдвое длиннее и пределом CFL вдвое ниже, затем — "
                      "ступень пополам (всего откатов до 4). Потом плавный второй порядок и давление до 0,75."),
        "hypothesis": ("Любое заранее заданное расписание может оказаться слишком быстрым на каком-то отрезке; "
                       "квазистационарный путь с проверкой каждого состояния и возвратом при срыве его найдёт."),
        "differs_from_smoke": "5 ступеней частоты с давлением ~ w^2, гейт на каждой, автоматический откат",
        "transforms": ["hist"], "set": dict(OUT_SET, **WANG), "backoff": {"max_total": 4},
        "steps": ([_st(f"r{k}_w{int(round(w * 100)):03d}", w, LIGHT_LOAD, 220, ramp=120, pramp=120,
                       hold={"iter": 200, "max": 3}) for k, w in enumerate((0.30, 0.55, 0.75, 0.90, 1.00), 1)]
                  + [_st("f_order2", 1.0, LIGHT_LOAD, 400, muscl=True, muscl_ramp=250, hold={"iter": 200, "max": 2}),
                     _st("g_base", 1.0, BASE_LOAD, 500, pramp=300, muscl=True, gate="final", hold=FINAL_HOLD)]),
    },
    {
        "id": "v4_jst", "axis": "численная схема",
        "title": "Центральная схема JST без переключения порядка, малый CFL с медленным ростом, LU_SGS",
        "mechanism": ("JST (второй порядок с искусственной вязкостью, JST_SENSOR_COEFF (0,5; 0,02) шаблона) с первой "
                      "итерации — переключения порядка нет (MUSCL с JST SU2 v8.6.0 не запускает); CFL 1 -> до 8 "
                      "(множитель роста 1,05, сброса 0,5), предобусловливатель LU_SGS — как в образцах турбомашин SU2 "
                      "(Aachen — JST, APU — LU_SGS); частота за 600, давление до 0,60 за 1100 (порция 1200), затем до 0,75 "
                      "за 500."),
        "hypothesis": ("Разгон Roe первого порядка с CFL до 20 и резким включением MUSCL слишком агрессивен для "
                       "тетраэдров без призм при U2 ~ 500–540 м/с; другая схема и осторожный CFL дойдут."),
        "differs_from_smoke": "схема JST вместо Roe, без MUSCL, CFL 1–8, LU_SGS, длинные рампы, гейты",
        "transforms": ["hist"],
        "set": dict(OUT_SET, FLOW_SCHEME="JST", LIN_PREC="LU_SGS", CFL_ADAPT_PARAM="0.5, 1.05, 0.5, 8.0, 0.1"),
        "steps": [_st("a_spin_light", 1.0, LIGHT_LOAD, 1200, ramp=600, pramp=1100, cfl=1.0),
                  _st("b_base", 1.0, BASE_LOAD, 700, pramp=500, cfl=2.0, gate="final",
                      hold={"iter": 300, "max": 40, "cfl": 4.0})],
    },
    {
        "id": "v5_sa", "axis": "модель турбулентности",
        "title": "Модель Спаларта — Аллмараса вместо SST на том же пути, что v2",
        "mechanism": ("Путь и численная схема v2 (частота при лёгком выходе -> плавный второй порядок -> нагружение), "
                      "но SA (одно уравнение, SA_OPTIONS= NONE) вместо SST. Вместо «сначала Эйлер или ламинарное "
                      "решение, потом RANS», которое в SU2 v8.6.0 с этим шаблоном не делается (см. diagnosis.not_used)."),
        "hypothesis": ("Срыв вызывает уравнение для w модели SST на сетке без пристенных призм (в истории пробы "
                       "rms[w] не записывался по уравнениям до NaN); SA мягче."),
        "differs_from_smoke": "SA вместо SST; остальное — как v2",
        "transforms": ["hist"], "set": dict(OUT_SET, **WANG), "turb": "SA",
        "steps": [_st("a_spin_light", 1.0, LIGHT_LOAD, 1000, ramp=450, pramp=900, cfl=5.0),
                  _st("b_order2", 1.0, LIGHT_LOAD, 400, muscl=True, muscl_ramp=250),
                  _st("c_base", 1.0, BASE_LOAD, 600, pramp=400, muscl=True, gate="final", hold=FINAL_HOLD)],
    },
    {
        "id": "v6_riemann", "axis": "граничное условие",
        "title": "Вход и выход Riemann вместо Giles, давление на выходе ступенями (лестница «частота, потом давление»)",
        "mechanism": ("MARKER_RIEMANN: вход TOTAL_CONDITIONS_PT (p0, T0, направление +z), выход STATIC_PRESSURE; "
                      "MARKER_TURBOMACHINERY остаётся (рампа частоты и TURBO_PERF). Рампы давления у Riemann в SU2 "
                      "v8.6.0 нет — давление ступенью между порциями: частота 0 -> 0,5 при p_нач, давление на линию "
                      "дросселя; 0,5 -> 0,8, давление; 0,8 -> 1,0, давление 0,60; плавный второй порядок; 0,675; 0,75. "
                      "Каждое состояние — с гейтом."),
        "hypothesis": ("Неотражающее условие Giles (осреднение 1D, релаксация 0,3) ведёт себя плохо при обратном токе "
                       "на выходе, который был в пробе к концу s1; характеристическое условие Riemann его терпит."),
        "differs_from_smoke": "Riemann вместо Giles на входе и выходе, давление ступенями, лестница, гейты",
        "transforms": ["hist", "riemann"], "set": dict(OUT_SET, **WANG),
        "steps": [_st("a_spin050", 0.5, None, 300, ramp=200, hold={"iter": 200, "max": 2}),
                  _st("b_p050", 0.5, LIGHT_LOAD, 200, hold={"iter": 200, "max": 3}),
                  _st("c_spin080", 0.8, None, 250, ramp=150, hold={"iter": 200, "max": 2}),
                  _st("d_p080", 0.8, LIGHT_LOAD, 200, hold={"iter": 200, "max": 3}),
                  _st("e_spin100", 1.0, None, 200, ramp=100, hold={"iter": 200, "max": 2}),
                  _st("f_p100", 1.0, LIGHT_LOAD, 200, hold={"iter": 200, "max": 3}),
                  _st("g_order2", 1.0, LIGHT_LOAD, 350, muscl=True, muscl_ramp=250, hold={"iter": 200, "max": 2}),
                  _st("h_p0675", 1.0, 0.675, 200, muscl=True, hold={"iter": 200, "max": 3}),
                  _st("i_base", 1.0, BASE_LOAD, 250, muscl=True, gate="final", hold=FINAL_HOLD)],
    },
]

# ----------------------------------------------------------------------------------------------- время и смета
WHEELS = {
    "kv3_stage1": dict(case="custom", geom_out="$W/geom/stage1_coarse", nodes=258411,
                       nodes_src="char_summary_kv3_stage1.json пробы c8f9 (сетка, построенная на AWS)"),
    "core130": dict(case="core130", geom_out="$W/geom/core130_coarse", nodes=139087,
                    nodes_src="char_summary_core130.json пробы c8f9"),
}
STEP_OVERHEAD_S = 20     # на запуск шага: run_case.py (копия сетки с разделёнными маркерами, сверка ключей) + чтение
                         # сетки и разбиение в SU2; оценка: в пробе 9 итераций KV-3 вместе с запуском — 7,1 с на 64 процессах
SETUP_MIN = 4.5          # загрузка + container_setup (сборка SU2 141 с) + сетки: в пробе char_points стартовал на 248-й
                         # секунде (budget.uptime_at_start_s, KV-3), здесь сетки строятся параллельно
FINISH_MIN = 1.5         # merge, выгрузка сводок, гашение (вся проба c8f9 — 305 с от запуска до конца по ledger)


def strategy_iters(strat: dict) -> tuple:
    """(итераций обязательной части, запусков SU2) — без порций удержания."""
    it, n = 0, 0
    for s in strat["steps"]:
        if s.get("kind") == "plan":
            it += sum(x["iter"] for x in RC.validate_plan(RC.PLANS[s["plan"]], 1.0))
        else:
            it += int(s["iter"])
        n += 1
    return it, n


def wheel_time(nodes: int, np_: int, iters: int, launches: int, c_us: float) -> dict:
    t_it = c_us * 1e-6 * nodes / np_ / M.EFF
    return {"s_per_iter": round(t_it, 3), "mandatory_s": round(launches * STEP_OVERHEAD_S + iters * t_it),
            "iters_in_limit": int((TIME_LIMIT_S - launches * STEP_OVERHEAD_S) / t_it)}


def estimate(jobs: list) -> dict:
    np_ = M.VCPU // 2
    rows, mach_min = [], 0.0
    for j in jobs:
        s = j["_strategy"]
        iters, launches = strategy_iters(s)
        w = {lab: {"expected": wheel_time(W["nodes"], np_, iters, launches, M.C_EXPECTED),
                   "sizing": wheel_time(W["nodes"], np_, iters, launches, M.C_SIZING)} for lab, W in WHEELS.items()}
        worst = max(v["sizing"]["mandatory_s"] for v in w.values())
        if worst > TIME_LIMIT_S:
            raise SystemExit(f"{s['id']}: обязательная часть {worst} с при C_sizing не влезает в срок {TIME_LIMIT_S} с")
        m_exp = SETUP_MIN + TIME_LIMIT_S / 60 + FINISH_MIN
        m_low = SETUP_MIN + max(v["expected"]["mandatory_s"] for v in w.values()) / 60 + FINISH_MIN
        mach_min += m_exp
        rows.append({"job": j["name"], "strategy": s["id"], "mandatory_iters_per_wheel": iters, "su2_launches_min": launches,
                     "wheels": w, "machine_min_expected": round(m_exp, 1), "machine_min_if_success_early": round(m_low, 1),
                     "max_min": j["max_min"]})
    usd_h = M.PRICE_EXPECTED + M.IPV4_USD_H + DISK_GB * M.DISK_USD_GB_H
    upper = len(jobs) * (M.PRICE_UPPER + DISK_GB * M.DISK_USD_GB_H + M.IPV4_USD_H) * (MAX_MIN + 10) / 60
    out_gb = len(jobs) * 2 * 0.06
    return {
        "machines": len(jobs), "instance_type": "c7a.16xlarge (запасной c7i.16xlarge)", "vcpu_per_machine": M.VCPU,
        "np_per_wheel": np_, "time_limit_per_wheel_s": TIME_LIMIT_S,
        "machine_hours_expected": round(mach_min / 60, 2),
        "usd": {"total_expected": round(mach_min / 60 * usd_h + out_gb * (M.S3_USD_GB_MONTH + M.EGRESS_USD_GB), 2),
                "spot_expected": round(mach_min / 60 * M.PRICE_EXPECTED, 2),
                "ipv4_disk": round(mach_min / 60 * (M.IPV4_USD_H + DISK_GB * M.DISK_USD_GB_H), 3),
                "s3_and_transfer": round(out_gb * (M.S3_USD_GB_MONTH + M.EGRESS_USD_GB), 3)},
        "upper_bound_usd": round(upper, 2),
        "upper_bound_rule": "раннер: (MaxPrice 1,4206 + диск + IPv4) x (max_min + 10)/60 на машину — если все дойдут до сторожа",
        "jobs": rows,
        "time_model": {
            "formula": ("машина = подготовка (загрузка, контейнер, сборка SU2, обе сетки параллельно) + срок стратегии "
                        f"{TIME_LIMIT_S} с (оба колеса параллельно, по {np_} процессов; при успехе — раньше) + сбор; "
                        "итерация t = C x узлы / процессы / EFF; запуск шага + "
                        f"{STEP_OVERHEAD_S} с"),
            "C_expected_us": M.C_EXPECTED, "C_sizing_us": M.C_SIZING, "EFF": M.EFF, "setup_min": SETUP_MIN,
            "finish_min": FINISH_MIN, "step_overhead_s": STEP_OVERHEAD_S,
            "nodes": {k: v["nodes"] for k, v in WHEELS.items()}, "nodes_source": {k: v["nodes_src"] for k, v in WHEELS.items()},
            "basis": ["C — замеры make_jobs_su2_kv3.py (общая 4-ядерная машина, грубо); на c7a.16xlarge, скорее всего, "
                      "быстрее — тогда в срок войдёт больше порций удержания, а не меньше",
                      "подготовка — ledger.jsonl пробы c8f9 (запуск -> конец 305 с на всю пробу) и budget.uptime_at_start_s "
                      "(248 с до старта SU2 у KV-3); сборка SU2 141 с — журнал прогона по заданию (в ветке нет)",
                      "STEP_OVERHEAD_S — оценка, не замер", "S3 и выход из AWS — ставки-допущения make_jobs_su2_kv3.py"],
        },
        "price": {"expected_usd_h": M.PRICE_EXPECTED, "upper_usd_h": M.PRICE_UPPER,
                  "basis": "ledger.jsonl: спот c7a.16xlarge us-east-1c 0,7064 (проба c8f9); MaxPrice раннера 1,4206"},
        "outputs": "на колесо: sweep_summary.json, история и журнал каждой порции, cfg, итоговый рестарт (~14–27 МБ), "
                   "рестарт перед NaN, если был; оценка выхода 6 x 2 x ~60 МБ",
    }


# ----------------------------------------------------------------------------------------------- разбор провала
def diagnosis() -> dict:
    kv = json.load(open(os.path.join(SMOKE_DIR, "char_summary_kv3_stage1.json"), encoding="utf-8"))
    co = json.load(open(os.path.join(SMOKE_DIR, "char_summary_core130.json"), encoding="utf-8"))
    scaled = RC.validate_plan(RC.PLANS["smoke075"], 0.1)
    full = RC.validate_plan(RC.PLANS["smoke075"], 1.0)
    rf = 5   # run_case.py --ramp-freq по умолчанию (char_points его не меняет)

    def sched(st, f0):
        return {i: round(SW.ramp_value(i, f0, st["frac"], st["ramp"], rf), 4) for i in range(0, st["iter"])
                if i in (0, rf, st["ramp"]) or i == st["iter"] - 1}

    def last(s):
        r = s["base"]["record"]
        p = r["performance"]
        return {"stage": r["stage"], "history_rows": r["history_rows"], "seconds_incl_su2_start": r["seconds"],
                "su2_error": r["su2_error"], "rms_rho": r["rms_rho"], "m_in_passage": p["m_in_passage"],
                "m_out_passage": p["m_out_passage"], "MassFlowOut_turbo": p["MassFlowOut_turbo"],
                "mass_imbalance": p["mass_imbalance"], "PR_tt": p["PR_tt"], "eta_tt": p["eta_tt_isentropic"],
                "PRTT_su2": p["PRTT_su2"], "TotTotEff_su2": p["TotTotEff_su2"],
                "periodic": r["periodic"], "nproc": s["nproc"], "p_base_Pa": s["p_base_Pa"],
                "p_out_estimate_Pa": s["estimate"]["p_out_estimate"], "mesh_nodes": s["mesh"]["header"]["NPOIN"]}
    return {
        "run": {"run_id": "20261115T003035Z-c8f9", "jobs_file": "jobs_aws_su2_kv3_smoke.json",
                "machine": "1 x c7a.16xlarge us-east-1c, спот 0,7064 долл./ч; запуск -> конец 305 с (ledger.jsonl)",
                "files": f"{REL}/aws_smoke_20261115T003035Z-c8f9/ (char_summary обоих колёс, smoke_merge.json, smoke_exits.txt)"},
        "what_was_run": {
            "command": "char_points.py run ... --base-iter-scale 0.1 --base-pfrac 0.75 (план базы smoke075, 64 процесса)",
            "plan_smoke075_full": full, "plan_smoke075_as_run_x0_1": scaled,
            "omega_frac_by_iteration_su2": {"s1_ramp60": sched(scaled[0], 0.0), "s2_ramp100_p075": sched(scaled[1], 0.6)},
            "omega_rule": ("CFluidIteration::UpdateRamps v8.6.0: значение = начало + i x (конец − начало)/длина, обновление при "
                           "i % 5 == 0 и i < длины или i == длине: при рампе 8 частота меняется скачками на 5-й и 8-й "
                           "итерациях (s1: 0 -> 0,375 -> 0,6 номинала, s2: 0,6 -> 0,85 -> 1,0)"),
        },
        "observed": {
            "end_of_s1_from_task": {"kv3_stage1": "17 итераций, m_вх 0,377 / m_вых 0,319 кг/с (оба > 0 — обратный ток на "
                                                  "выходе), дисбаланс 185 %, PR_tt 0,843, rms[Rho] ~ −4,3",
                                    "core130": "17 итераций, дисбаланс 165 %",
                                    "source": "журнал прогона по заданию владельца; в char_summary — только последняя стадия"},
            "last_stage_kv3_stage1": last(kv), "last_stage_core130": last(co),
            "notes": ["NaN на 9-й (KV-3) и 10-й (core130) итерации s2 — через 1–2 итерации после того, как частота "
                      "скачком дошла до 1,0 на 8-й итерации",
                      "до NaN невязка плотности в истории не росла лавиной (rms[Rho] KV-3: начало −4,29, максимум −3,93, "
                      "последняя −4,32): срыв местный и внезапный, а не постепенный",
                      "в последней строке s2 расход на выходе обратный у обоих колёс (Avg_Massflow(outlet) > 0, "
                      "MassFlowOut_1 < 0), PR_tt 0,99 и 1,11, КПД нефизичен (TotTotEff KV-3 −200 %)",
                      "для сравнения — прогон Б (README, core130 coarse, тот же путь целиком, по частям): после s1 (170 "
                      "итераций) дисбаланс 43,5 %, после s2 (200) — 23,2 %, после s3 (150, MUSCL) — 1,4 %: даже полный "
                      "s1 не устанавливает расход, но и не расходится"]},
        "cause": [
            "Главная: --base-iter-scale 0.1 сжал не только число итераций, но и длины рамп (run_case.validate_plan "
            "масштабирует ramp и pramp вместе с iter): разгон 0 -> 0,6 -> 1,0 номинала прошёл за 8 + 8 итераций скачками "
            "по 0,25–0,375 номинала (обновление рампы раз в 5 итераций) — по сути импульсный пуск колеса. В прогоне Б "
            "тот же путь занимал 80 + 80 итераций, и core130 его выдержал.",
            "Вторая: переход s1 -> s2 — по числу итераций. s2 стартовал с поля, где расход не установился и на выходе уже "
            "шёл обратный ток (дисбаланс 165–185 %), и за 16 итераций поднимал давление до 0,75 оценки.",
            "Сопутствующая: каждый рестарт начинает адаптивный CFL заново с 5 (рост x1,1 за итерацию до 20) — к 9-й "
            "итерации s2 CFL мог быть до ~12 (5 x 1,1^9; в истории CFL не писался — проверить нельзя).",
        ],
        "not_established": [
            "место NaN в поле: объёмное поле и рестарт перед NaN не записывались (OUTPUT_WRT_FREQ 500, стадия короче)",
            "CFL и линейный решатель по итерациям: в HISTORY_OUTPUT шаблона нет групп CFL_NUMBER и LINSOL",
            "какое уравнение сорвалось первым: максимальных невязок (MAX_RES) в истории нет, по rms до NaN роста нет",
            "полные истории s1 обоих колёс в сессии не лежат (в ветке — только char_summary последней стадии)",
            "KV-3 по полному пути smoke075 не считался ни разу — то, что он выдержит полный путь, не проверено",
        ],
        "fixed_in_sweep": [
            "история с группами CFL_NUMBER, LINSOL, MAX_RES (производный шаблон sweep_ramp.py rc; run_case.py и шаблон "
            "не тронуты)",
            "рестарт каждые 100 итераций (OUTPUT_WRT_FREQ): при NaN остаётся поле <= 100 итераций до срыва "
            "(restart_before_nan.dat), разбор — SU2_SOL в сессии",
            "при NaN — divergence в sweep_summary.json: итерация, частота и давление по рампам на ней, какие невязки "
            "росли за 10 строк, CFL, знак расхода",
            "переход между стадиями — по гейту расхода и невязки, а не по числу итераций",
        ],
        "checked_in_su2_source_v8_6_0": [
            "RAMP_ROTATING_FRAME и RAMP_OUTLET_PRESSURE из постановки — не ключи v8.6.0: старые имена RAMP_ROTATION_FRAME, "
            "RAMP_OUTLET_PRESSURE отклоняются («… is deprectaed», CConfig.cpp); нужны RAMP_MOTION_FRAME и RAMP_OUTLET "
            "(шаблон их и использует)",
            "RAMP_OUTLET ставит значение Giles для STATIC_PRESSURE, STATIC_PRESSURE_1D, RADIAL_EQUILIBRIUM и "
            "MASS_FLOW_OUTLET; с Riemann — ошибка «Outlet pressure ramp only implemented for NRBC» (CFluidIteration.cpp)",
            "RAMP_MUSCL, RAMP_MUSCL_COEFF (начало, частота обновления, длина), KIND_MUSCL_RAMP (ITERATION | "
            "SMOOTH_FUNCTION), RAMP_MUSCL_POWER: множитель реконструкции 0 -> 1 (CEulerSolver.cpp, CFluidIteration.cpp)",
            "JST с явно заданным MUSCL_FLOW= YES — ошибка «Centered schemes do not use MUSCL reconstruction» (CConfig.cpp); "
            "шаблон пишет MUSCL_FLOW всегда",
            "Venkatakrishnan: eps^2 = (REF_ELEM_LENGTH x K)^3, REF_ELEM_LENGTH = 1,0 (СИ) — при K = 0,05 это 1,25e-4 против "
            "квадратов перепадов давления в Па: коэффициент на ограничитель почти не влияет; VENKATAKRISHNAN_WANG: "
            "eps^2 = (K x размах переменной)^2 (CLimiterDetails.hpp)",
            "CFL_ADAPT_PARAM = (множитель сброса, множитель роста, мин., макс.[, допустимая невязка линейного "
            "решателя][, итерация начала])",
            "повтор ключа в cfg — ошибка «option appears twice»: ключи шаблона меняются только его полями (--set) или "
            "производным шаблоном, --extra — только новые ключи",
            "начальное поле уже свободный поток по оси: run_case.py — M = 0,4 вдоль +z (AOA 90), p и T из p01, T01",
        ],
        "not_used": {
            "euler_or_laminar_then_rans": (
                "CTurbSolver::LoadRestart читает турбулентные переменные из столбцов сразу после консервативных: в "
                "рестарте Эйлера или ламинарного решения с VOLUME_OUTPUT (…, PRIMITIVE) шаблона там давление и "
                "температура, без PRIMITIVE — нули, и SST получает k = w = 0 (мю_t = 0/0). Эмуляция «RANS со всеми "
                "стенками без трения» тоже не годится: без вязких стенок SU2 ставит расстояние до стенки 0 во всей "
                "области (CGeometry::ComputeWallDistance), а SST и SA на него делят. Вместо этого — ось «модель "
                "турбулентности» (v5_sa)."),
            "giles_mass_flow_outlet": (
                "MASS_FLOW_OUTLET у Giles в v8.6.0 есть и рампуется, но образца в TestCases нет, коэффициент релаксации "
                "давления (второе значение) размерный и не подобран, знак AverageMassFlowRate выхода сверен только "
                "косвенно (MassFlowOut_1 > 0 при прямом течении, прогон Б) — кандидат следующего перебора."),
            "classic_venkat_coefficient_scan": "коэффициент классического Venkatakrishnan в размерном расчёте не действует (см. выше)",
        },
    }


# ----------------------------------------------------------------------------------------------- команда задачи
def command(strat_file: str, mini: dict | None = None) -> str:
    env = list(M.ENV_LINES)
    if mini:
        env = [ln if not ln.startswith("W=") else "W=/gtu/work; mkdir -p \"$W\"" for ln in env]
    lines = env + ["H=$(( NPROC / 2 )); [ \"$H\" -ge 1 ] || H=1" if not mini else f"H={mini['np']}"]
    mesh_kv = ("python routes/gtu/frame3/frame3_geom.py --out \"$W/geom\" --levels coarse --stages stage1 --gzip-levels "
               "--threads \"$H\" --results \"$GTU_OUT/frame3_geom_results.json\" > \"$GTU_OUT/geom_kv3_stage1.log\" 2>&1")
    mesh_co = ("python routes/gtu/platform/geom/impeller.py --case core130 --level coarse --no-step "
               "--out \"$W/geom/core130_coarse\" --threads \"$H\" > \"$GTU_OUT/geom_core130.log\" 2>&1")
    if mini:   # миниатюры по частям: сетки строятся один раз и лежат в примонтированной $W
        mesh_kv = f"[ -f \"$W/geom/stage1_coarse/passage_custom_coarse.su2\" ] || {mesh_kv}"
        mesh_co = f"[ -f \"$W/geom/core130_coarse/passage_core130_coarse.su2\" ] || {mesh_co}"
    lines += [f"( {mesh_kv} ) & G1=$!", f"( {mesh_co} ) & G2=$!", "wait \"$G1\"; wait \"$G2\"",
              f"S={REL}/{strat_file}"]
    for k, (lab, W) in enumerate(WHEELS.items(), 1):
        args = ["--strategy \"$S\"", f"--case {W['case']}", "--level coarse", f"--geom-out \"{W['geom_out']}\"",
                f"--label {lab}", "--np \"$H\"", "--su2-src /opt/SU2-src"]
        if mini:
            args += [f"--mini-iter {mini['iter']}", f"--time-limit {mini['time_limit']}"]
        else:
            args += [f"--time-limit {TIME_LIMIT_S}", f"--job-max-min {MAX_MIN}"]
        args += [f"--work \"$W/sw_{lab}\"", f"--out \"$GTU_OUT/{lab}\""]
        lines.append(f"python {REL}/sweep_ramp.py run " + " ".join(args) + f" > \"$GTU_OUT/sweep_{lab}.out\" 2>&1 & P{k}=$!")
    for k, lab in enumerate(WHEELS, 1):
        lines.append(f"wait \"$P{k}\" || echo \"{lab} exit $?\" >> \"$GTU_OUT/sweep_exits.txt\"")
    lines.append(f"python {REL}/sweep_ramp.py merge --in \"$GTU_OUT\" --out \"$GTU_OUT/sweep_merge.json\" "
                 "--expect " + ",".join(WHEELS))
    return "\n".join(lines)


def strategy_doc(s: dict) -> dict:
    doc = {k: v for k, v in s.items()}
    doc["_comment"] = ("Стратегия разгона базовой точки для sweep_ramp.py run (перебор jobs_aws_su2_kv3_sweep.json). "
                       "Сгенерировано make_jobs_su2_kv3_sweep.py; правки — в нём. load — доля оценки давления на выходе "
                       "при номинальной частоте (char_points.p_out_estimate); при частоте w < 1 давление — по линии "
                       "дросселя p_нач + w^2 (load x оценка − p_нач); load null — давление прежнее. set — поля шаблона "
                       "run_case.py (--set), transforms — правки производного шаблона (sweep_ramp.py rc).")
    return doc


def criteria_doc() -> dict:
    est = {lab: json.load(open(os.path.join(SMOKE_DIR, f"char_summary_{lab}.json"), encoding="utf-8"))["estimate"]
           ["p_out_estimate"] for lab in WHEELS}
    pr = {"kv3_stage1": 4.49, "core130": 4.25}
    return {
        "base_point": {lab: {"frac": 1.0, "load": BASE_LOAD, "p_out_Pa": round(BASE_LOAD * est[lab]),
                             "p_out_estimate_Pa": round(est[lab]), "PR_1d": pr[lab]} for lab in WHEELS},
        "success": ["дисбаланс расхода вход/выход <= 1 % (среднее за 50 строк) при прямом течении (вход > 0 > выход)",
                    "расход на входе меняется <= 0,5 % за 200 итераций, невязка плотности не растёт",
                    "rms[Rho] ниже первой итерации цепочки >= 3 порядков (цель — 4)",
                    "п и КПД физичны: 1 < PR_tt <= 1,3 x PR одномерной оценки, TR_tt > 1, 0,5 <= eta_tt < 1, "
                    "TotTotEff SU2 в [50; 100) %, PRTT SU2 и свой PR_tt расходятся <= 5 %; без NaN"],
        "stricter_skill_su2_turbo": "converged_skill: >= 4 порядков и п, КПД меняются < 0,1 % за 500 итераций",
        "gate_between_stages": "дисбаланс <= 2 % (среднее за 20 строк), расход меняется <= 1 % за 100 итераций, "
                               "невязка не растёт (наклон <= +0,05 порядка/100), прямое течение",
        "verdicts": "success > established (всё, кроме 3 порядков) > reached_base > stuck > timeout > diverged > error",
        "thresholds_source": "пороги гейтов приняты здесь (sweep_ramp.GATES), не из документации SU2",
    }


def build() -> dict:
    jobs = []
    for s in STRATEGIES:
        SW.validate_strategy(s)
        fn = f"plan_sweep_{s['id']}.json"
        with open(os.path.join(HERE, fn), "w", encoding="utf-8") as f:
            json.dump(strategy_doc(s), f, ensure_ascii=False, indent=1)
            f.write("\n")
        jobs.append({"name": f"su2-kv3-sweep-{s['id'].replace('_', '-')}", "command": command(fn), "max_min": MAX_MIN,
                     "disk_gb": DISK_GB, "_strategy": s, "_file": fn})
    est = estimate(jobs)
    table = []
    for j in jobs:
        s = j["_strategy"]
        er = next(r for r in est["jobs"] if r["job"] == j["name"])
        table.append({"job": j["name"], "file": f"{REL}/{j['_file']}", "id": s["id"], "axis": s["axis"],
                      "title": s["title"], "mechanism": s["mechanism"], "hypothesis": s["hypothesis"],
                      "differs_from_smoke": s["differs_from_smoke"], "steps": [x["name"] for x in s["steps"]],
                      "turbulence": s.get("turb", "SST"), "scheme": s.get("set", {}).get("FLOW_SCHEME", "ROE"),
                      "outlet": "Riemann STATIC_PRESSURE" if "riemann" in s.get("transforms", []) else "Giles STATIC_PRESSURE_1D",
                      "backoff": bool(s.get("backoff")),
                      "expected_min_machine": er["machine_min_expected"],
                      "mandatory_min_kv3": round(er["wheels"]["kv3_stage1"]["expected"]["mandatory_s"] / 60, 1),
                      "mandatory_min_core130": round(er["wheels"]["core130"]["expected"]["mandatory_s"] / 60, 1)})
    out_json = "<scratchpad>/su2_kv3_sweep_<run_id>.json"
    return {
        "_comment": ("Перебор стратегий разгона базовой точки SU2: 6 машин c7a.16xlarge, по стратегии на машину, на каждой "
                     "оба колеса (KV-3 ступень 1 и core130, coarse). Раннер берёт только defaults и jobs; остальное — для "
                     f"сессии. Код — ветка {M.BRANCH} (запуск — после того, как sweep_ramp.py и plan_sweep_*.json "
                     "окажутся в ветке). Сгенерировано make_jobs_su2_kv3_sweep.py; правки — в нём."),
        "defaults": {"instance_type": "c7a.16xlarge", "fallback_types": ["c7i.16xlarge"], "image": "al2023",
                     "container": M.CONTAINER, "container_setup": M.CONTAINER_SETUP},
        "jobs": [{k: v for k, v in j.items() if not k.startswith("_")} for j in jobs],
        "diagnosis": diagnosis(),
        "strategies": table,
        "success_criteria": criteria_doc(),
        "outputs_per_job": ("$GTU_OUT: env.txt, geom_*.log, frame3_geom_results.json, sweep_<колесо>.out, sweep_exits.txt "
                            "(если колесо кончилось не 0), sweep_merge.json; <колесо>/: sweep_summary.json (шаги, гейты, "
                            "метрики, вердикт, divergence), <NN_шаг[_hK]>/: history_*.csv (с CFL, LINSOL, MAX_RES), su2_*.log, "
                            "case_*.cfg, run_summary.json, run_case.out; final_restart.dat; restart_before_nan.dat при NaN"),
        "merge": {"where": "сессия с ключами после run: routes/gtu/platform/aws/runs/<run_id>/<задача>/",
                  "commands": [f"python3 {REL}/sweep_ramp.py merge --in routes/gtu/platform/aws/runs/<run_id> "
                               f"--out {out_json} --run-id <run_id> --expect {','.join(WHEELS)}"],
                  "to_branch": ["routes/gtu/platform/aws/ledger.jsonl"],
                  "then": ("победитель по колесу (best_per_wheel) -> база этапа A: char_points.py run --base-strategy "
                           f"{REL}/plan_sweep_<id>.json (линии частоты на medium; jobs_aws_su2_kv3.json пересобрать)"),
                  "note": "числа перебора — проверка разгона на coarse без призм, в публикацию не идут"},
        "estimate": est,
    }


def mini_files(d: str) -> None:
    os.makedirs(d, exist_ok=True)
    mini = {"np": 2, "iter": 3, "time_limit": 3000}
    for s in STRATEGIES:
        with open(os.path.join(d, f"cmd_{s['id']}.sh"), "w", encoding="utf-8") as f:
            f.write(command(f"plan_sweep_{s['id']}.json", mini) + "\n")
    with open(os.path.join(d, "run_in_container.sh"), "w", encoding="utf-8") as f:
        f.write("""#!/bin/bash
# Эмуляция RUN_CONTAINER раннера (как миниатюры make_jobs_su2_kv3.py): root, копия репозитория (git archive HEAD +
# новые файлы) -> /gtu/repo, $OUT -> /gtu/out, рабочая папка с сетками -> /gtu/work (общая для миниатюр — сетки строятся
# один раз); окружение образа узла docker снимается — на AWS в ubuntu:24.04 его нет, команда задаёт его сама.
# Запуск: bash run_in_container.sh <id стратегии> <папка repo> <папка out> <папка work> <исходник SU2>
ID=$1; REPO=$2; OUTD=$3; WORK=$4; SRC=$5; MINI=$(cd "$(dirname "$0")" && pwd)
mkdir -p "$OUTD" "$WORK"
docker run --rm --name "gtu-sweep-$ID" --user 0:0 --entrypoint /bin/bash --cpus 4 --shm-size 2g -w /gtu/repo \\
  -v "$REPO:/gtu/repo" -v "$OUTD:/gtu/out" -v "$WORK:/gtu/work" -v "$MINI:/gtu/mini:ro" -v "$SRC:/opt/SU2-src:ro" \\
  -e GTU_OUT=/gtu/out -e GTU_RUN_ID=mini-local -e GTU_JOB="su2-kv3-sweep-$ID" -e GTU_ATTEMPT=1 -e NPROC=4 \\
  -e OMP_NUM_THREADS=4 -e GTU_ITYPE=local -e GTU_AZ_MACHINE=local -e DEBIAN_FRONTEND=noninteractive \\
  -e OMPI_ALLOW_RUN_AS_ROOT=1 -e OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1 \\
  gtu-platform:local -c "env -u OMPI_MCA_osc -u OMPI_MCA_btl_vader_single_copy_mechanism -u SU2_RUN -u GTU_HOME \\
    -u PYTHONDONTWRITEBYTECODE PATH=/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin bash /gtu/mini/cmd_$ID.sh"
echo "exit $?"
""")
    print("миниатюры:", d)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--no-plan", action="store_true")
    ap.add_argument("--mini", default=None, metavar="DIR")
    a = ap.parse_args(argv)
    if a.mini:
        mini_files(a.mini)
        return 0
    doc = build()
    with open(OUT_JOBS, "w", encoding="utf-8") as f:
        json.dump(doc, f, ensure_ascii=False, indent=1)
    if not a.no_plan:
        doc["plan_check"] = M.plan_check(OUT_JOBS)
        with open(OUT_JOBS, "w", encoding="utf-8") as f:
            json.dump(doc, f, ensure_ascii=False, indent=1)
    e = doc["estimate"]
    print(f"{os.path.basename(OUT_JOBS)}: задач {len(doc['jobs'])}, машино-часов {e['machine_hours_expected']}, ожидается "
          f"{e['usd']['total_expected']} долл., верхняя граница {e['upper_bound_usd']} долл."
          + (f"; plan: код {doc['plan_check']['exit_code']}" if "plan_check" in doc else ""))
    for r in doc["strategies"]:
        print(f"  {r['id']}: {r['axis']}; машина ~{r['expected_min_machine']} мин; обязательная часть KV-3 "
              f"{r['mandatory_min_kv3']} мин, core130 {r['mandatory_min_core130']} мин")
    return 0 if doc.get("plan_check", {}).get("exit_code", 0) == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
