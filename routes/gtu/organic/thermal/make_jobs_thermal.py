"""Генератор jobs_thermal.json — задачи AWS узла organic/thermal (формат routes/gtu/platform/aws/runner.py).

Запускаемые задачи — только код этой папки, проверенный локально на малых размерах (README, «Проверки»):
  tpms-mac      — таблица переносов TPMS в пределе Стокса на сетках 48³ и 64³, три сдвига сетки (сходимость);
  tpms-lbm-48   — течение с инерцией и теплообмен в ячейке (LBM D3Q19/D3Q7), Re 10…800, сетка 48³;
  tpms-lbm-64   — то же на 64³ для двух уровней (проверка сетки);
  topo2d-hires  — оптимизация топологии коллектора на сетке 240 по ширине, доли жидкости 0,2/0,3/0,4, оба случая;
  vascular      — сети охлаждения (воспроизведение локального расчёта на машине: та же версия кода).
Трёхмерные задачи с CFD и сопряжённым теплообменом — блок planned_3d: постановка и оценка, без команды
(кода сборки кейсов OpenFOAM здесь нет — он не проверен и в пакет не включён).
Цена спота — из routes/gtu/platform/aws/README.md (прогон 14.11.2026: 0,706–0,955 долл./ч); оценки времени —
по локальным замерам ×(запас), всё — оценка.
"""
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
P = "routes/gtu/organic/thermal"
PIP = ["numpy==2.5.3", "scipy==1.18.1", "scikit-image==0.26.0", "pyamg==5.4.1", "cantera==3.2.0",
       "matplotlib==3.11.2", "trimesh==5.2.1"]
PRICE = dict(low=0.706, high=0.955, src="routes/gtu/platform/aws/README.md (цены спота c7a.16xlarge 14.11.2026)")


def jobs():
    J = []
    J.append(dict(
        name="thermal-tpms-mac",
        command=(f"python {P}/tpms_transport.py --method mac --n 48 64 --shifts 0 0.33 0.67 "
                 f"--levels -0.6 -0.4 -0.2 0.0 0.2 0.4 0.6 0.8 --pe 0.01 10 20 40 --workers 24 "
                 f"--out \"$GTU_OUT/tpms_transport_mac64.json\""),
        max_min=150, disk_gb=30))
    re_all = "10 50 150 400 800"
    J.append(dict(
        name="thermal-tpms-lbm-48",
        command=(f"python {P}/tpms_transport.py --method lbm --n 48 --levels -0.2 0.0 0.2 0.4 --re {re_all} "
                 f"--u-lu 0.04 --max-steps 120000 --workers \"$NPROC\" --out \"$GTU_OUT/tpms_transport_lbm48.json\""),
        max_min=360, disk_gb=30))
    J.append(dict(
        name="thermal-tpms-lbm-64",
        command=(f"python {P}/tpms_transport.py --method lbm --n 64 --levels 0.0 0.4 --re {re_all} "
                 f"--u-lu 0.04 --max-steps 160000 --workers \"$NPROC\" --out \"$GTU_OUT/tpms_transport_lbm64.json\""),
        max_min=600, disk_gb=30))
    runs = []
    for case in ("comb6", "fan8"):
        for vol in ("0.20", "0.30", "0.40"):
            runs.append(f"mkdir -p \"$GTU_OUT/{case}_{vol}\" && python {P}/topo.py --nx 240 --vol {vol} --cases {case} "
                        f"--out \"$GTU_OUT/{case}_{vol}\" > \"$GTU_OUT/{case}_{vol}.log\" 2>&1")
    J.append(dict(
        name="thermal-topo2d-hires",
        command=" &\n".join(runs) + " &\nwait\nls -lR \"$GTU_OUT\" > \"$GTU_OUT/_list.txt\"",
        max_min=240, disk_gb=30))
    J.append(dict(
        name="thermal-vascular",
        command=f"python {P}/vascular.py --out \"$GTU_OUT\"",
        max_min=90, disk_gb=20))
    return J


def estimate(J):
    rows = []
    hours = dict(**{"thermal-tpms-mac": 2.0, "thermal-tpms-lbm-48": 5.0, "thermal-tpms-lbm-64": 9.0,
                    "thermal-topo2d-hires": 3.0, "thermal-vascular": 0.6})
    tot_hi = 0.0
    for j in J:
        h = hours[j["name"]]
        cap = (j["max_min"] + 10) / 60
        rows.append(dict(job=j["name"], expected_h=h, cap_h=round(cap, 2), usd_expected=round(h * PRICE["low"], 2),
                         usd_cap=round(cap * PRICE["high"] * 1.2, 2)))
        tot_hi += cap * PRICE["high"] * 1.2
    return dict(rows=rows, usd_cap_total=round(tot_hi, 2),
                usd_expected_total=round(sum(r["usd_expected"] for r in rows), 2), price=PRICE,
                note="ожидаемые часы — локальные замеры при свободной машине, пересчитанные на 64 ядра с запасом ×2–3 "
                     "(оценка); верхняя граница — срок сторожа × максимум цены × 1,2, как считает предохранитель раннера")


PLANNED_3D = [
    dict(name="manifold-3d-topo-cfd",
         what="трёхмерная оптимизация топологии топливного/воздушного коллектора KV-3: турбулентное течение (RANS, "
              "Re ~ 10⁵–10⁶), цель — потери полного давления при равных расходах на выходах (разброс ≤ 2 %), объём ≤ V",
         how="плотностный метод (Бринкман) с сопряжённым решателем: OpenFOAM v2412 adjointOptimisationFoam, режим "
             "topology (наличие режима в этой версии — не проверено) или собственный 3D-код (не написан)",
         size="сетка ~2–5 млн ячеек, 200–400 итераций оптимизации, 64 ядра MPI", est_hours="8–20 (оценка)",
         blocked_by="нет проверенного кода сборки кейса; нет контейнера OpenFOAM в пакете платформы"),
    dict(name="vane-cht-vascular",
         what="сопряжённый теплообмен лопатки СА с сетью каналов (прямые / 8 деревьев / растущее дерево): газ RANS + "
              "охладитель + теплопроводность металла; проверка 1D+2D-модели vascular.py (x_ngv, T_max, разброс)",
         how="chtMultiRegionFoam (OpenFOAM, навык openfoam-combustor) или SU2 CHT; геометрия каналов — vascular_ngv_*.stl",
         size="~3–8 млн ячеек на лопатку с периодичностью по шагу", est_hours="6–12 на вариант (оценка)",
         blocked_by="нет сборщика кейса и профиля лопатки (есть только развёртка стороны)"),
    dict(name="tpms-module-cht",
         what="модуль рекуператора с распределителями: оба лабиринта гироида с переменной ячейкой + стенка; проверка "
              "сердцевины и коллекторов (неравномерность по фронту)",
         how="chtMultiRegionFoam на сетке snappyHexMesh по tpms_graded_module_section.stl (scratchpad)",
         size="сечение 40 × 40 мм на всю длину, ~20–40 млн ячеек", est_hours="10–30 (оценка)",
         blocked_by="нет сборщика кейса; большой STL хранится вне репозитория"),
]


def main():
    J = jobs()
    out = {
        "_comment": ("Узел organic/thermal (README.md). Раннер берёт только defaults и jobs; блоки estimate, planned_3d, "
                     "merge — для сессии с ключами. Код — ветка claude/lucid-allen-yvps5d. Сгенерировано "
                     "make_jobs_thermal.py; правки — в нём. AWS из этой сессии не вызывался."),
        "defaults": {"instance_type": "c7a.16xlarge", "fallback_types": ["c7i.16xlarge"], "image": "al2023",
                     "disk_gb": 30, "python": "3.12", "pip": PIP},
        "jobs": J,
        "estimate": estimate(J),
        "merge": {"where": "сессия с ключами после run: runs/<run_id>/<задача>/",
                  "commands": [
                      f"python3 {P}/tpms.py --mac routes/gtu/platform/aws/runs/<run_id>/thermal-tpms-mac/tpms_transport_mac64.json "
                      f"--lbm routes/gtu/platform/aws/runs/<run_id>/thermal-tpms-lbm-48/tpms_transport_lbm48.json --out {P}",
                  ],
                  "checks": "сверить f·Re и Nu на 32³ (этот узел) и 48³/64³ — разница — оценка ошибки сетки; LBM при "
                            "Re 10 — с пределом Стокса MAC (±5 %)"},
        "planned_3d": PLANNED_3D,
    }
    (HERE / "jobs_thermal.json").write_text(json.dumps(out, ensure_ascii=False, indent=1))
    print(json.dumps(out["estimate"], ensure_ascii=False, indent=1))


if __name__ == "__main__":
    main()
