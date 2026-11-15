"""Сборка verify_phys_results.json — проверка физики узлов routes/gtu/organic (exergy, cycles, structure, thermal),
перезапуск в чистой копии и независимые проверки (15.11.2026). Скрипты проверок лежат рядом:
  grid3_profiles.py   — проверенные профили колеса и диска в disk.py/CalculiX на сетках 80×12, 160×24, 320×48;
  topo_convergence.py — сходимость disk_topo.py по журналам итераций и вторая (вдвое крупнее) сетка;
  reformer_check.py   — равновесие SMR своим решателем на эмпирических K, тепловой баланс риформера, активность углерода
                        по составной кривой, кинетика Xu — Froment при условиях варианта;
  plate_supported.py  — пластинчатый рекуператор с рёбрами-опорами против перепада давления (f·Re и Nu_H1
                        прямоугольного канала своим конечно-разностным расчётом), модель сердцевины tpms_recup;
  exergy_h2o.py       — влияние конденсации воды выхлопа на эксергию и эксергия CH4 по табличным ΔG_f;
  topo2d_grids.py     — коллектор thermal/topo.py на второй сетке (nx 120 против 80);
  cmpjson.py          — сравнение results-JSON перезапуска с репозиторием.
Запуск: python3 collect.py <папка с выходами проверок> <выход.json>"""
import json, os, platform, sys
import numpy as np, scipy, cantera
D, OUT = sys.argv[1], sys.argv[2]
L = lambda f: json.load(open(os.path.join(D, f)))
res = dict(
    title="Проверка физики organic: эксергия, циклы и риформинг, топология и прочность, теплообмен (перезапуск в чистой копии)",
    method="Перезапуск exergy.py, cycles.py, vascular.py, tpms.py (с готовой таблицей переносов mac32), ccx-проверки "
           "воксельных проектов (structure_post.py --ccx) и тестов в чистой копии routes/gtu; сравнение results-JSON с "
           "репозиторием (cmpjson.py); независимые проверки — скрипты рядом (см. docstring collect.py). Всё — оценка уровня "
           "облика, как и проверяемые узлы.",
    software=f"python {platform.python_version()}, numpy {np.__version__}, scipy {scipy.__version__}, cantera {cantera.__version__}, "
             "ccx 2.21, gmsh 4.15.2",
    date="2026-11-15",
    inputs=[dict(file="routes/gtu/organic/exergy/exergy_results.json", field="весь файл (сравнение перезапуска)"),
            dict(file="routes/gtu/organic/cycles/cycles_results.json", field="results[*] (сравнение), reformer, results[*steam*].opt"),
            dict(file="routes/gtu/organic/structure/profile_opt_comp2.json", field="best.x"),
            dict(file="routes/gtu/organic/structure/profile_opt_eq.json", field="best.x"),
            dict(file="routes/gtu/organic/structure/profile_opt_lim.json", field="best.x"),
            dict(file="routes/gtu/organic/structure/disk_topo_results.json", field="runs, baseline"),
            dict(file="routes/gtu/organic/structure/ccx_check_results.json", field="results (сравнение)"),
            dict(file="routes/gtu/organic/thermal/tpms_results.json", field="recuperator.cases[R90-ic-RH-ibc-none-none]"),
            dict(file="routes/gtu/organic/thermal/vascular_results.json", field="весь файл (сравнение)"),
            dict(file="<scratchpad>/organic/log_<прогон>.txt", field="журналы итераций disk_topo.py (вне репозитория)")],
    rerun=dict(
        exergy=dict(numeric_diffs=0, note="совпадение до последнего знака; tests 11/12 до исправления (date бралась как дата "
                    "перезапуска), 12/12 после закрепления DATE = 2026-11-14 в exergy.py", tests_after_fix="12 passed"),
        cycles=dict(max_rel_diff=1.75e-4, where="results[522].pr135.fuel_info.pinch_K (шаг поиска пинча)", tests="22 passed",
                    runtime_s=579),
        structure_ccx=dict(note="structure_post.py --ccx для comp2_A_mass_eq и turb1_A_mass_eq и исходных профилей — загрузки "
                           "совпали до последнего знака: колесо 1,4374 / 1,3430, диск 1,2749 / 1,4142", tests="test_structure 6/6"),
        thermal=dict(vascular_numeric_diffs=0, tpms_recuperator_numeric_diffs=0, topo_nx80_numeric_diffs=0, tests="10 passed",
                     note="в vascular_results.json репозитория были поля networks.*.volume_mm3, cases.ngv.equal_flow_note и другой "
                          "текст geometry.note, которых vascular.py не писал (правка файла мимо скрипта) — скрипт дополнен; "
                          "у tpms.py иллюстративный STL в репозитории прорежен постобработкой (40 000 против 96 000 треугольников)")),
    exergy=L("exergy_h2o.json"),
    cycles_reformer=L("reformer_check.json"),
    structure_grid3=L("grid3_profiles.json"),
    structure_topology=L("topo_convergence.json"),
    thermal_plate_supported=L("plate_supported.json"),
)
if os.path.exists(os.path.join(D, "vascular_rerun2_cmp.json")):
    res["rerun"]["thermal"]["vascular_after_fix"] = L("vascular_rerun2_cmp.json")
if os.path.exists(os.path.join(D, "topo2d_grids.json")):
    res["thermal_topo2d_grids"] = L("topo2d_grids.json")
json.dump(res, open(OUT, "w", encoding="utf-8"), ensure_ascii=False, indent=1, default=float)
print("->", OUT)
