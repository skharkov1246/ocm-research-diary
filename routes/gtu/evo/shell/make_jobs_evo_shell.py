#!/usr/bin/env python3
"""Генератор jobs_evo_shell.json — задачи AWS семейства «корпус-оболочка и рама» для раннера routes/gtu/platform/aws/runner.py
и списка элит verify_elites.json для проверки в CalculiX.

Пакет 1 «острова» — MAP-Elites: 2 варианта × 10 частей, 640 островов на вариант (64 на машину, по процессу на vCPU),
30 000 оценок на остров: 19,2 млн оценок на вариант. Часть — острова с номером ≡ part (mod 10); итог части —
архивы островов в $GTU_OUT/islands/<вариант>/. Слияние (shell_evo.py merge) — после collect, локально: оно проверяет,
что пришли все 640 островов с одной конфигурацией (хэш кода модели), иначе — ошибка с перечнем недостающих.
Пакет 2 «проверка» — CalculiX (контейнер ubuntu:24.04 + calculix-ccx), сетка fine, элиты из verify_elites.json
(лучшие, самые новые и края фронта NSGA-II мини-прогона + типовые формы + гладкий цилиндр): случай k — в части k mod 4.

Оценка времени — по мини-прогону (≈ 6–10 мс на оценку на занятой машине 4 ядра; c7a — не медленнее, оценка) и по
расчётам ccx грубой сетки (≈ 70–120 с в 2 потока; fine — 6,4× элементов, время ≈ ×20–25 — оценка); сторож — с запасом ×3.
Проверка без AWS: python3 routes/gtu/platform/aws/runner.py plan routes/gtu/evo/shell/jobs_evo_shell.json --assume-price 1.0
Перед запуском: закоммитить routes/gtu/evo (код, verify_elites.json) в ветку прогона.
"""
import argparse
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
OUT = os.path.join(HERE, "jobs_evo_shell.json")
VER = os.path.join(HERE, "verify_elites.json")
S = "routes/gtu/evo/shell"
PIP = ["numpy==2.5.3", "scipy==1.18.1"]
N_ISL, N_PARTS, EVALS = 640, 10, 30000
CCX_PARTS, CCX_WORKERS, CCX_THREADS = 4, 6, 10


def verify_cases(merged, nsga, top=12, novel=6):
    import numpy as np
    import shell_evo as SE
    import shell_model as SM
    cases = []
    for tag, variant in (("o", "оболочка"), ("k", "картридж")):
        if tag not in merged:
            continue
        m = json.load(open(merged[tag], encoding="utf-8"))
        ctx = SE.effect_ctx(m, variant)
        fe = [e for e in m["archive"]["elites"] if e["objs"].get("CV", 1) == 0]
        fe.sort(key=lambda e: (-e["fitness"], e["uid"]))
        chosen = {e["uid"]: e for e in fe[:top]}
        nov = sorted(fe, key=lambda e: -float(SE.novelty(np.array(e["meta"]["feats"]), ctx)[0][0]))
        for e in nov:
            if len(chosen) >= top + novel:
                break
            chosen.setdefault(e["uid"], e)
        for uid, e in chosen.items():
            cases.append(dict(kind="elite", uid=uid, variant=variant, genome=e["genome"], fitness=e["fitness"], cell=e["cell"]))
        if tag in nsga:
            fr = json.load(open(nsga[tag], encoding="utf-8"))["front"]
            for j, name in ((0, "наименьшая масса"), (1, "наименьший N_eff"), (4, "наибольшая новизна")):
                if fr:
                    r = min(fr, key=lambda r: r["F"][j])
                    cases.append(dict(kind="elite", uid=f"nsga-{tag}-{j}", variant=variant, genome=r["genome"], note=name))
        for k in SM.REF_KINDS:
            cases.append(dict(kind="reference", uid=f"ref-{tag}-{k}", variant=variant, reference=k))
    cases.append(dict(kind="canonical", uid="canonical"))
    for k, c in enumerate(cases):
        c["tag"] = f"c{k:03d}"
    return cases


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--merged-o"); ap.add_argument("--merged-k"); ap.add_argument("--nsga-o"); ap.add_argument("--nsga-k")
    a = ap.parse_args()
    if a.merged_o or a.merged_k:
        merged = {k: v for k, v in (("o", a.merged_o), ("k", a.merged_k)) if v}
        nsga = {k: v for k, v in (("o", a.nsga_o), ("k", a.nsga_k)) if v}
        cases = verify_cases(merged, nsga)
        with open(VER, "w", encoding="utf-8") as f:
            json.dump(dict(_comment="элиты мини-прогона для проверки в CalculiX (shell_ccx.py --batch); сгенерировано make_jobs_evo_shell.py",
                           cases=cases), f, ensure_ascii=False)
        print("->", VER, len(cases), "случаев")
    jobs = []
    for v in ("shell", "cartridge"):
        for part in range(N_PARTS):
            jobs.append(dict(name=f"evo-shell-{v}-p{part:02d}", max_min=45, python="3.13", pip=PIP,
                             command=f"python -u {S}/shell_evo.py run-part --variant {v} --n-islands {N_ISL} --n-parts {N_PARTS} "
                                     f"--part {part} --evals {EVALS} --workers \"$NPROC\" --out \"$GTU_OUT\""))
    for part in range(CCX_PARTS):
        jobs.append(dict(name=f"evo-shell-ccx-p{part}", max_min=150, container="ubuntu:24.04",
                         container_setup=["apt-get update -q",
                                          "apt-get install -y -q --no-install-recommends calculix-ccx python3 python3-numpy python3-scipy"],
                         command=f"python3 -u {S}/shell_ccx.py --batch {S}/verify_elites.json --part {part} --n-parts {CCX_PARTS} "
                                 f"--workers {CCX_WORKERS} --threads {CCX_THREADS} --mesh fine --out \"$GTU_OUT\""))
    d = {"_comment": "Задачи семейства routes/gtu/evo/shell (README семейства, раздел «AWS»). Раннер берёт только defaults и jobs. "
                     "Пакет 1 — острова MAP-Elites (20 задач), пакет 2 — CalculiX элит мини-прогона (4 задачи). Перед прогоном "
                     "закоммитить routes/gtu/evo вместе с verify_elites.json. Сгенерировано make_jobs_evo_shell.py; правки — в нём.",
         "defaults": {"instance_type": "c7a.16xlarge", "fallback_types": ["c7i.16xlarge"], "image": "al2023", "disk_gb": 40},
         "jobs": jobs,
         "estimate": {"islands": f"{2 * N_ISL} островов × {EVALS} оценок = {2 * N_ISL * EVALS / 1e6:.1f} млн оценок",
                      "islands_basis": "6–10 мс на оценку (мини-прогон, 4 занятых ядра); 64 острова на машину — ≈ 3–5 мин счёта + установка "
                                       "≈ 5 мин (оценка)",
                      "ccx": "случаев verify_elites.json по частям 4 × ≈ 11; сетка fine (120 × 240 × 4 слоя); ≈ 15–40 мин на случай "
                             "в 10 потоков (оценка по грубой сетке)"}}
    with open(OUT, "w", encoding="utf-8") as f:
        json.dump(d, f, ensure_ascii=False, indent=1)
    print("->", OUT, len(jobs), "задач")


if __name__ == "__main__":
    main()
