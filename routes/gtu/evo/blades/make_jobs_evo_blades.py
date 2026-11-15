#!/usr/bin/env python3
"""Задачи AWS семейства «лопатки и колёса» → jobs_evo_blades.json (раннер routes/gtu/platform/aws/runner.py).

Стадия 1 (этот файл): N_PARTS машин c7a.16xlarge; на каждой — часть островов по номеру (islands_of_part), на каждый
компонент (колёса ступеней 1 и 2, турбина) по одному острову на ядро, EVALS оценок на остров; затем слияние островов
части, свои элиты и их проверка в CalculiX (disk.py) — всё в контейнере ubuntu:24.04 (CalculiX из apt, numpy, scipy,
cantera — в venv теми же версиями, что у задач воронки). После сбора — слияние всех частей локально с контролем
полноты (evo_blades.py merge-parts), галерея, выбор владельца → стадия 2 (--parents).

  python3 make_jobs_evo_blades.py [--parts 24] [--evals 20000]
  python3 ../../platform/aws/runner.py plan jobs_evo_blades.json --assume-price 1.0     # без AWS
"""
import argparse
import json
import os

HERE = os.path.dirname(os.path.abspath(__file__))
REL = "routes/gtu/evo/blades"
NPROC = 64                     # c7a.16xlarge: 64 vCPU — по острову на ядро
PIP = ["numpy==2.5.3", "scipy==1.18.1", "cantera==3.2.0"]     # как funnel/jobs_l1.json
# время ЦП оценки на одном ядре — мини-прогон 15.11.2026 (blades_mini_results.json → comps.*.timing.ms_per_eval_cpu)
T_EVAL_S = {"comp1": 0.0214, "comp2": 0.0208, "turb": 0.0374}


def build(parts, evals, k_elites, comps=("comp1", "comp2", "turb")):
    n_isl = parts * NPROC
    compute_min = evals * sum(T_EVAL_S[c] for c in comps) / 60.0 * 1.25      # запас 25 % на разброс и слияние
    max_min = int(min(720, max(30, round(compute_min + 15 + 10))))           # + установка в контейнере, CalculiX элит
    setup = ["apt-get update -q",
             "apt-get install -y -q --no-install-recommends calculix-ccx python3 python3-venv python3-pip",
             "python3 -m venv /opt/evo",
             "/opt/evo/bin/pip install -q " + " ".join(PIP)]
    jobs = []
    for p in range(parts):
        cmd = (f"export EVO_OUT=/tmp/evo_blades\n"
               f"/opt/evo/bin/python -u {REL}/evo_blades.py part --part {p} --n-parts {parts} --n-islands {n_isl} "
               f"--comps {','.join(comps)} --evals {evals} --k {k_elites} --workers \"$NPROC\" --ccx --scale aws --out \"$GTU_OUT\"")
        jobs.append(dict(name=f"blades-p{p:03d}", container="ubuntu:24.04", container_setup=setup, command=cmd, max_min=max_min))
    total = n_isl * evals * len(comps)
    return dict(
        _comment=("Семейство «лопатки и колёса» (routes/gtu/evo/blades, README): стадия 1 — острова MAP-Elites по частям, элиты части → "
                  f"CalculiX. Островов на компонент {n_isl} ({parts} частей × {NPROC}), оценок на остров {evals}, всего оценок "
                  + f"{total:,}".replace(",", " ") + f". Срок задачи {max_min} мин (расчёт ≈ {compute_min:.0f} мин + установка + CalculiX). "
                  "Сгенерировано make_jobs_evo_blades.py; правки — в нём. После сбора: evo_blades.py merge-parts (контроль полноты)."),
        defaults=dict(instance_type="c7a.16xlarge", fallback_types=["c7i.16xlarge"], image="al2023", disk_gb=30),
        jobs=jobs)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--parts", type=int, default=24)
    ap.add_argument("--evals", type=int, default=20000)
    ap.add_argument("--k", type=int, default=8)
    ap.add_argument("--out", default=os.path.join(HERE, "jobs_evo_blades.json"))
    a = ap.parse_args()
    d = build(a.parts, a.evals, a.k)
    with open(a.out, "w", encoding="utf-8") as f:
        json.dump(d, f, ensure_ascii=False, indent=1)
    print(d["_comment"])


if __name__ == "__main__":
    main()
