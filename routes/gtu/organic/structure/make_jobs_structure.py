#!/usr/bin/env python3
"""Генератор jobs_structure.json — задачи AWS узла organic/structure для раннера routes/gtu/platform/aws/runner.py.

Здесь (4 ядра, соседние сессии) — малые демонстрации; на AWS (c7a.16xlarge: 64 vCPU, 128 ГБ) — то, что не помещается:
  struct-lattice-n80     — гомогенизация гироида на сетке 80³ (сходимость прочности: здесь N 48, плотные — 32);
  struct-axi-fine        — все осесимметричные варианты на вдвое более мелкой сетке (сходимость по сетке);
  struct-topo3d-disk-h20 — 3D диск турбины ступени 1, 1/8 модели, шаг 2,0 мм (≈0,6 млн степеней свободы);
  struct-topo3d-disk-h15 — то же, шаг 1,5 мм (≈1,3 млн; PARDISO в 64 потока, память ~30–60 ГБ — оценка);
  struct-topo3d-frame-n160 — 3D рама оболочки, 160 вокселей по диаметру (≈1,3 млн степеней свободы).
Оценка времени — по локальной демонстрации (topo3d_demo_results.json) и росту прямого решателя ~N^1,5 (оценка);
сторож max_min — с запасом ×2. Перед прогоном код и lattice_results.json, disk_topo_results.json должны быть в ветке.
Проверка без AWS: python3 routes/gtu/platform/aws/runner.py plan routes/gtu/organic/structure/jobs_structure.json --assume-price 1.0
"""
import json
import os

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "jobs_structure.json")
S = "routes/gtu/organic/structure"
PIP = ["numpy==2.5.3", "scipy==1.18.1", "scikit-image==0.26.0", "trimesh==5.2.1", "shapely==2.2.0",
       "matplotlib==3.11.2", "pypardiso==0.4.7", "pyamg==5.4.1"]


def main():
    jobs = [
        dict(name="struct-lattice-n80", max_min=90,
             command=f"python -u {S}/lattice.py --N 80 --n-dense 80 --conv --out \"$GTU_OUT/lattice_results_n80.json\""),
        dict(name="struct-axi-fine", max_min=300,
             # радиус фильтра — тот же в метрах, что у основного прогона (2,5 / 3,75 мм): иначе при 2,5h он уменьшается
             # вдвое вместе с сеткой и меняется сама постановка, а не только сетка (проверка 15.11.2026)
             command=f"ORGANIC_OUT=\"$GTU_OUT\" python -u {S}/disk_topo.py --all --h-t 0.5e-3 --h-c 0.75e-3 "
                     f"--r-filter-t 2.5e-3 --r-filter-c 3.75e-3 --workers 9 "
                     f"--results \"$GTU_OUT/disk_topo_results_fine.json\""),
        dict(name="struct-topo3d-disk-h20", max_min=360,
             command=f"python -u {S}/topo3d.py --problem disk --h 2.0e-3 --iters 160 --out \"$GTU_OUT\""),
        dict(name="struct-topo3d-disk-h15", max_min=720,
             command=f"python -u {S}/topo3d.py --problem disk --h 1.5e-3 --iters 160 --out \"$GTU_OUT\""),
        dict(name="struct-topo3d-frame-n160", max_min=480,
             command=f"python -u {S}/topo3d.py --problem frame --n 160 --iters 150 --out \"$GTU_OUT\""),
    ]
    d = {"_comment": "Задачи узла routes/gtu/organic/structure (README узла, раздел «AWS»). Раннер берёт только defaults и jobs. "
                     "Код — ветка claude/lucid-allen-yvps5d; перед прогоном закоммитить узел вместе с lattice_results.json и "
                     "disk_topo_results.json (цель по загрузке 3D-диска). Сгенерировано make_jobs_structure.py; правки — в нём.",
         "defaults": {"instance_type": "c7a.16xlarge", "fallback_types": ["c7i.16xlarge"], "image": "al2023", "disk_gb": 40,
                      "python": "3.13", "pip": PIP},
         "jobs": jobs}
    with open(OUT, "w", encoding="utf-8") as f:
        json.dump(d, f, ensure_ascii=False, indent=1)
    print("->", OUT, len(jobs), "задач")


if __name__ == "__main__":
    main()
