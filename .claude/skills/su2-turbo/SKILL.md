---
name: su2-turbo
description: Аэродинамика ступеней трека «Простая литая ГТУ-заменитель» в SU2 — установка (сборка v8.6.0 «Harrier»), сетка одного межлопаточного канала из routes/gtu/platform/geom/impeller.py (Gmsh / Open CASCADE, формат .su2, маркеры), конфигурация сжимаемого RANS (SST/SA, вращающаяся система, периодичность, полные параметры на входе, статическое давление на выходе, Giles, mixing-plane для ступени), сходимость и контроль расхода, полная степень сжатия и изоэнтропический КПД, характеристика и граница устойчивости, проверка на эталонах NASA Rotor 37 / Stage 35, NASA CC3, DLR SRV2; каждый ключ сверен с исходником SU2 v8.6.0. Используй, когда готовишь или разбираешь расчёт SU2 колеса компрессора, турбины или ступени, строишь характеристику или сверяешь SU2 с эталоном. Триггеры — «SU2», «SU2_CFD», «расчёт колеса», «межлопаточный канал», «mixing plane», «характеристика компрессора», «граница помпажа», «Rotor 37», «CC3», «SRV2», «КПД ступени».
---

# SU2: ступени центробежного компрессора и турбины

Узел платформы трека (issue #59): одномерные размеры (`gtu_explore.comp_stage`, навык `gtu-cycle`) → сетка
канала (`platform/geom/impeller.py`) → SU2 → степень сжатия, КПД, характеристика. Состояние на 14.11.2026:
**SU2 в контейнере не установлен, ни одного прогона SU2 не было.** Ключи и поведение сверены с исходником
SU2 v8.6.0 и исходником документации (раздел «Источники»). Пометки: «(код)» — по исходнику, «(док.)» — по
документации, «(пример)» — по учебной конфигурации SU2, «не проверено» — как есть.

## 1. Что где лежит

| что | где |
|---|---|
| геометрия и сетка колеса | `routes/gtu/platform/geom/impeller.py`; описание и ограничения — `platform/geom/README.md` |
| вход случаев `core130`, `taurus60` | `gtu_explore_results.json` (схема 130), `gtu_dropin_results.json` (1Ц, класс N) |
| числа для SU2 после сетки | `out/<case>_<level>/impeller_<case>_<level>.json` → `su2.periodic.angle_deg`, `su2.rotation_rate_rad_s`, `d1.p01`, `d1.T01`, `d1.PR`, `d1.eta_comp_stage` |
| шаблон конфигурации колеса | `.claude/skills/su2-turbo/templates/impeller_rans.cfg` (58 ключей, все есть в v8.6.0) |
| π, η и расход из `history.csv` | `.claude/skills/su2-turbo/scripts/su2_pr_eta.py` |
| облачные прогоны | `platform/aws/runner.py` (DryRun по умолчанию, живой запуск — `--live`), README там же |
| узел SU2 трека | `routes/gtu/platform/su2/` — папки ещё нет; имя из примера задачи в `platform/aws/README.md` |

Сетку турбины репозиторий не строит: `impeller.py` делает только колесо компрессора.

## 2. Установка

- **Пакета apt и pip нет:** Ubuntu 24.04 `apt-cache search --names-only su2` пусто (индексы noble на месте,
  `gmsh` в них есть); PyPI `su2` и `pysu2` — 404. Проверено 14.11.2026.
- **Сборка из исходников** (док. Build-SU2-Linux-MacOS; v8.6.0 «Harrier», тег 30.10.2026, meson ≥ 1.8.2;
  `meson.py` сам скачивает meson, ninja и внешние библиотеки с GitHub):
  ```bash
  git clone --depth 1 -b v8.6.0 https://github.com/su2code/SU2.git && cd SU2
  ./meson.py setup build -Dwith-mpi=disabled --prefix=$HOME/su2   # MPI: -Dwith-mpi=enabled; OpenMP: -Dwith-omp=true
  ./ninja -C build -j4 install
  export SU2_RUN=$HOME/su2/bin PATH=$HOME/su2/bin:$PATH PYTHONPATH=$HOME/su2/bin:$PYTHONPATH
  ```
  Здесь сборка не запускалась (отклонена системой разрешений сессии): время, пакеты apt для MPI и папка
  `bin` под `--prefix` не проверены.
- **Готовые архивы:** страница загрузки берёт из последнего выпуска GitHub файлы `*linux64-omp.zip` и
  `*linux64-mpi.zip` (код `download.html`); есть ли они у v8.6.0 — не проверено (API GitHub из сессии закрыт).
- Запуск (док. Execution, код `SU2_CFD.cpp`): `SU2_CFD case.cfg`; `mpirun -n 4 SU2_CFD case.cfg`; сборка с
  OpenMP — `SU2_CFD -t 4 case.cfg`; `SU2_CFD -d case.cfg` — сухой прогон: разбор конфигурации на фиктивной геометрии.

## 3. Сетка одного канала

```bash
cd routes/gtu/platform/geom && python impeller.py --case core130 --level coarse
grep -n "NDIME\|NMARK\|MARKER_TAG" out/core130_coarse/passage_core130_coarse.su2
```

- Маркеры `inlet outlet hub shroud blade periodic_1 periodic_2`, объём `fluid`, только тетраэдры. Gmsh пишет
  имя физической группы в `MARKER_TAG`, тетраэдр — код 10 (проверено записью куба в Gmsh 4.15.2). Сетка в метрах.
- Ось машины — `z`: радиус сечений по высоте SU2 считает от оси `z` (док., код).
- **Периодичность — узел в узел** (код `CGeometry::MatchPeriodic`): точки первого маркера поворачиваются на угол
  (в градусах), донор ищется в пределах 1·10⁻⁶ единиц сетки. Промах — «Bad match for point …», расчёт идёт
  дальше («Computation will continue, but be cautious»). В журнале нужна строка «Matched N points with a max
  distance of: …» без «Bad match». У `impeller.py` `periodic_2` = поворот `periodic_1` на +360/Z вокруг +z,
  значит порядок `(periodic_1, periodic_2, …, угол > 0)`.
- Пристенных призм нет (README узла geom): для RANS — пристенные функции или доработка сетки. В Gmsh
  `extrudeBoundaryLayer` есть только у встроенного ядра (`gmsh.model.geo`), у `gmsh.model.occ` нет (проверено
  `dir()`); как совместить с областью OCC из `impeller.py` — не проверено.
- Маркер `hub` идёт от носка входа через колесо до r = 1,3·r2 (`impeller.py`, `hub_poly`): во вращающейся системе
  вращается и стенка безлопаточного диффузора. Пока маркер не разделён (не сделано), числа по диффузору — с оговоркой.

## 4. Конфигурация: одна решётка (колесо)

Основа — `templates/impeller_rans.cfg`: `<…>` заменить числами из JSON сетки (раздел 1).

- `GRID_MOVEMENT= ROTATING_FRAME`: сетка неподвижна, `ROTATION_RATE` — вектор в рад/с, `MACH_MOTION` во всех
  примерах 0,35 (док.). Источник вращения в коде — ω×ρu от консервативных переменных, то есть скорость в решении
  абсолютная (вывод из `flow_sources.cpp`, прямой строки в документации нет): полные параметры на входе и в
  `SURFACE_TOTAL_PRESSURE` — абсолютные. Знак ω — только из JSON сетки.
- `MARKER_SHROUD`: скорость сетки там обнуляется, покрышка не вращается (док., код); в `MARKER_HEATFLUX` она тоже
  нужна, как в Aachen (пример).
- `INIT_OPTION= TD_CONDITIONS` с малым `MACH_NUMBER`: умолчание `REYNOLDS` пересчитывает давление из числа
  Рейнольдса (док.). Модели: SST по умолчанию — `V2003m`, SA — SA-noft2, есть `SA_OPTIONS= ROTATION` (док., код).
- **С турбомашинными маркерами** (нужны для Giles и `TURBO_PERF`): `MARKER_TURBOMACHINERY= (inlet, outlet)`,
  `TURBOMACHINERY_KIND= AXIAL_CENTRIFUGAL`, `TURBO_PERF_KIND= COMPRESSOR`; вместо `MARKER_INLET/OUTLET` —
  `MARKER_GILES= (inlet, TOTAL_CONDITIONS_PT, <p01>, <T01>, 1.0, 0.0, 0.0, 1.0, 0.0, outlet, STATIC_PRESSURE_1D, <p>, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0)`
  (9 значений: маркер, тип, два значения, три направления, релаксация среднего и Фурье — док.);
  `AVERAGE_PROCESS_KIND= MIXEDOUT`, `PERFORMANCE_AVERAGE_PROCESS_KIND= MIXEDOUT`, `TURBO_PERF` в `HISTORY_OUTPUT`.
  На тетраэдрах — `SPANWISE_KIND= EQUISPACED` и `NUM_SPANWISE_SECTIONS`: `AUTOMATIC` считает узлы на стыке входа и
  периодической грани (док., «to be reviewed»). До запуска — `mkdir TURBOMACHINERY`, иначе профили по высоте не
  пишутся (док.). Большие π и ω — рампы `RAMP_OUTLET`, `RAMP_MOTION_FRAME` с `*_COEFF` = (начало, шаг, длина),
  задавать всегда (док.); рампа ω работает только с `MARKER_TURBOMACHINERY` (док.).
- **`TURBOMACHINERY_KIND` по коду** (`CPhysicalGeometry.cpp`: по какой координате мерится высота): колесо с осевым
  входом и радиальным выходом — `AXIAL_CENTRIFUGAL`; радиальный диффузор — `CENTRIFUGAL`; рабочее колесо
  радиальной турбины — `CENTRIPETAL_AXIAL`; её сопловой аппарат — `CENTRIPETAL`; осевой венец — `AXIAL`.

## 5. Ступень: mixing plane

Образцы в репозитории SU2: `TestCases/turbomachinery/Aachen_turbine/aachen_3D_MP_restart.cfg` (3 венца, 3D),
`APU_turbocharger/Jones_restart.cfg` (радиальная турбина: статор + ротор).

- Главный файл: `MULTIZONE= YES`, `CONFIG_LIST= (zone_1.cfg, zone_2.cfg)` по потоку, `OUTER_ITER` вместо `ITER`;
  файл зоны — только отличия: `GRID_MOVEMENT`, `ROTATION_RATE`, рампа (док., пример).
- Стык — парой (верх по потоку, низ) в `MARKER_ZONE_INTERFACE` и `MARKER_MIXINGPLANE_INTERFACE`, на нём Giles
  `MIXING_OUT` и `MIXING_IN`; `MIXINGPLANE_INTERFACE_KIND= LINEAR_INTERPOLATION`, `TURBULENT_MIXINGPLANE= YES` (примеры).
- `MARKER_PERIODIC` и `MARKER_TURBOMACHINERY` — по зонам по порядку: сечения зоны n ищутся по n-й периодической
  паре (код: `PeriodicBoundary == iZone + 1`). `TURBOMACHINERY_KIND`, `TURBO_PERF_KIND` — по значению на зону,
  иначе «Insufficient TURBO_PERF_KIND options» (док.).
- Поля — с номером зоны: `RMS_DENSITY[1]`, `TotalTotalEfficiency[1]`; `TURBO_PERF` SU2 сам относит к последней
  зоне (код `COutput.cpp`). Тип расчёта КПД ступени — из `TURBO_PERF_KIND` зоны 0 (код `CTurboOutput.cpp`).

## 6. Сходимость и расход

- Останов — невязки (`CONV_FIELD= RMS_DENSITY`, `CONV_RESIDUAL_MINVAL`) или Коши по полям `MARKER_ANALYZE`
  (`CONV_CAUCHY_ELEMS`, `CONV_CAUCHY_EPS`); поля `TURBO_PERF` в `CONV_FIELD` не работают (док.).
- Расход — столбцы `Avg_Massflow(inlet)`, `Avg_Massflow(outlet)` в `history.csv`, знаки противоположны: нормаль
  границы в SU2 смотрит внутрь области (код `CEulerSolver.cpp`). `Avg_Massflow` без скобок — сумма со знаком,
  то есть дисбаланс; `Avg_TotalPress` без скобок — сумма давлений, не среднее (код `CFlowOutput.cpp`).
- Принято в этом навыке (не из документации): точка сошлась, если невязки упали не меньше чем на 4 порядка,
  дисбаланс расхода ≤ 0,5 %, π и η за последние 500 итераций меняются меньше чем на 0,1 %.
- Giles и осреднённый Riemann работают только на самой мелкой сетке; неустойчиво с многосеточным — `MGLEVEL= 0` (док.).

## 7. Степень сжатия и КПД

`python3 .claude/skills/su2-turbo/scripts/su2_pr_eta.py history.csv 1.4` — последняя строка файла:
π₀ = p₀₂/p₀₁ и η = (π₀^((γ−1)/γ) − 1)/(T₀₂/T₀₁ − 1) по столбцам со скобками (осреднение по расходу —
умолчание `MARKER_ANALYZE_AVERAGE`, код), модули расходов и дисбаланс. Формула — для идеального газа с постоянным
γ (как `IDEAL_GAS`). Разбор проверен на искусственном файле в формате писателя SU2 (по коду: заголовки в
кавычках, столбцы выровнены пробелами), на настоящем `history.csv` — нет.

- `TURBO_PERF` (код `CTurboOutput.cpp`): `PressureRatioTT`, `TotalTotalEfficiency` (в процентах); изоэнтропическая
  полная энтальпия там — h(s₁, p₂) + V₂²/2, а не h(s₁, p₀₂). В отчёт — оба КПД, каждый со своим определением.
- Приведённый расход характеристики — ṁ·√(T₀₁/288,15)/(p₀₁/101 325).
- С одномерной оценкой (`d1.PR`, `d1.eta_comp_stage`) — сравнение, не проверка: расхождение записать, одномерную
  модель под SU2 не подгонять.

## 8. Характеристика и граница устойчивости

- Линия частоты: давление на выходе от низкого (запирание: расход перестаёт расти) вверх; каждая точка — с решения
  предыдущей: `RESTART_SOL= YES`, `SOLUTION_FILENAME=` её `restart` без расширения (SU2 добавляет расширение, в
  мультизоне — и номер зоны; код `CConfig::GetFilename`).
- У вершины характеристики, где π почти не зависит от давления, — Giles `MASS_FLOW_OUTLET` (значение 1 — расход,
  значение 2 — релаксация; док.).
- Последняя сошедшаяся стационарная точка — граница устойчивости стационарного расчёта, а не помпаж: помпаж
  зависит от системы (объёмы, сеть) и требует нестационарного расчёта, которого в навыке нет. Так и подписывать.
- Другая частота — пересчёт `ROTATION_RATE`; каждая точка — отдельная задача `runner.py`.

## 9. Проверка на открытых эталонах

Числа таблицы — **предварительно**: из выдержек поиска, первоисточники (NTRS, AGARD, ERCOFTAC, DLR) из контейнера
не открылись. Отдельно сверена только арифметика U = π·D·n/60.

| эталон | что проверяет | числа (предварительно) | геометрия и данные |
|---|---|---|---|
| NASA Rotor 37 | трансзвуковой венец, ROTATING_FRAME, скачки, зазор | 36 лопаток; 17 188,7 об/мин; 454,14 м/с (D ≈ 0,505 м); π 2,106; 20,19 кг/с; зазор 0,356 мм (AGARD-AR-355); в TP-1659 проектная π 2,05 | страница NASA TMR даёт архив геометрии и данных на nasa.gov — не скачан |
| NASA Stage 35 | mixing plane, осевая ступень | проектная π 1,82; пиковые КПД ротора и ступени 0,872 и 0,845 при π 1,875 и 1,842; запас 21,8 % (TP-1338) | NASA TP-1338 — не открыт |
| NASA CC3 | колесо с полулопатками, ближе всего к `core130` | 4,54 кг/с; π 4; 21 789 об/мин; 15 + 15 лопаток; загиб 50°; диффузор 24 лопатки | координаты — NASA CR-204134 (1997), не открыт; NASA HECC (2026): изготовленное колесо расходится с проектом |
| DLR SRV2-O | трансзвуковой вход, ближе к `taurus60` (M₁s отн. 1,41) | 13 + 13; 50 000 об/мин; 2,55 кг/с; π проектная 6,1–6,2, измеренная ступени 5,7; 586 м/с; r2 112 мм (586,4 м/с по арифметике) | ERCOFTAC KB AC 6-08 — не открыт, доступность геометрии не проверена |

Порядок: Rotor 37 (открыты геометрия и данные) — отладка постановки; затем CC3 и SRV2-O. Допусков согласия
документация не задаёт: в PR — отклонение π, η и расхода в расчётной точке и форма линии частоты против опыта,
без подгонки констант модели под эталон.

## 10. Обязательные правила

- Ключи — только из `Common/src/CConfig.cpp` той же версии; `config_template.cfg` неполон (нет, например,
  `MULTIZONE`). Неизвестный ключ останавливает SU2: «invalid option name. Check current SU2 options in
  config_template.cfg» (код). Проверка имён до запуска (значения — `SU2_CFD -d case.cfg`):
  ```bash
  curl -sS https://raw.githubusercontent.com/su2code/SU2/v8.6.0/Common/src/CConfig.cpp -o /tmp/CConfig.cpp
  grep -oE '^[A-Z][A-Z0-9_]*' case.cfg | sort -u | while read k; do grep -q "\"$k\"" /tmp/CConfig.cpp || echo "нет ключа: $k"; done
  ```
- Публикуемое число — из results-JSON рядом со скриптом (навык `gtu-cycle`), вместе с конфигурацией, журналом и
  версией SU2; «SU2 показал» без журнала прогона не писать.
- Локально — не дольше 20 минут; дольше — `runner.py` (us-east-1, спот, потолок 300 долл., запись только в
  `s3://$ALPHA_O_S3/<префикс>/gtu/…`), живой запуск только флагом `--live`. Ключей и секретов в файлах нет.

## 11. Типовые ошибки

- Скопированы `TestCases/turbomachinery/centrifugal_blade` и `centrifugal_stage`: в них удалённые ключи
  `MARKER_NRBC`, `GRID_MOVEMENT_KIND`, `ROTATION_RATE_Z`, `MARKER_TURBO_PERFORMANCE`, `MIXING_PROCESS_TYPE`,
  `MARKER_MIXINGPLANE` (поймано проверкой выше); в регрессионных тестах их нет.
- Старые имена `RAMP_OUTLET_PRESSURE`, `RAMP_ROTATION_FRAME` — ошибка «… is deprectaed» (код); нужны `RAMP_OUTLET`,
  `RAMP_MOTION_FRAME`. Рампа давления с выходом Riemann — «Outlet pressure ramp only implemented for NRBC» (док.).
- Перепутан порядок: `MARKER_INLET` — T₀, затем p₀; Giles и Riemann `TOTAL_CONDITIONS_PT` — p₀, затем T₀.
- Пары `MARKER_PERIODIC` в обратном порядке или угол со знаком «минус»: SU2 не останавливается, печатает «Bad match».
- `TURBOMACHINERY_KIND= CENTRIFUGAL` у колеса с осевым входом: высота на входе мерится по `z` — нужно `AXIAL_CENTRIFUGAL`.
- `MARKER_GILES` без `MARKER_TURBOMACHINERY` — «Giles Boundary conditions can only be used with turbomachinery markers» (док.).
- y⁺ не проверен на сетке без призм: поле `Y_PLUS` (группа `PRIMITIVE`) — в `VOLUME_OUTPUT`.

## 12. Проверка результата

1. Журнал: «Matched … max distance» без «Bad match»; маркеры и число узлов — как в JSON сетки.
2. Невязки и дисбаланс расхода — по разделу 6; расход на входе и выходе по модулю совпадают.
3. π и η — двумя способами (раздел 7: свой расчёт и `TURBO_PERF`), расхождение объяснено определением.
4. Сетка coarse → medium (`impeller.py --level medium`): изменение π и η записано; на одном уровне вывода нет.
5. До колёс трека — Rotor 37 (раздел 9) с записанным отклонением от опыта.

## Источники (сверено 14.11.2026)

- Исходник SU2 v8.6.0 «Harrier» (клон тега): https://github.com/su2code/SU2/tree/v8.6.0 — `Common/src/CConfig.cpp`,
  `Common/include/option_structure.hpp` и `.inl`, `Common/src/geometry/CGeometry.cpp`, `CPhysicalGeometry.cpp`,
  `SU2_CFD/src/output/CTurboOutput.cpp`, `CFlowOutput.cpp`, `COutput.cpp`, `SU2_CFD/include/output/COutput.hpp`,
  `SU2_CFD/src/numerics/flow/flow_sources.cpp`, `SU2_CFD/src/solvers/CEulerSolver.cpp`, `SU2_CFD/src/SU2_CFD.cpp`,
  `config_template.cfg`, `meson.build`, `meson_options.txt`, `TestCases/turbomachinery/*`, `TestCases/serial_regression.py`.
- Документация — исходник сайта https://github.com/su2code/su2code.github.io (коммит 67d0f94, 11.11.2026); сам
  сайт из контейнера не открылся (DNS): https://su2code.github.io/docs_v7/Setup-Guide-Turbomachinery/ ,
  https://su2code.github.io/docs_v7/Markers-and-BC/ , https://su2code.github.io/docs_v7/Build-SU2-Linux-MacOS/ ,
  https://su2code.github.io/docs_v7/SU2-Linux-MacOS/ , https://su2code.github.io/docs_v7/Execution/ ,
  https://su2code.github.io/docs_v7/Turbulence-Models/ , https://su2code.github.io/tutorials/Aachen_Turbine/ ,
  https://su2code.github.io/download.html ; сетки примеров — https://github.com/su2code/TestCases (список `turbomachinery/`).
- Gmsh 4.15.2: справка API установленного пакета; запись `.su2` проверена на кубе.
- NASA TMR, Rotor 37: https://github.com/tmbwg/turbmodels/blob/main/Other_exp_Data/rotor37_exp.html
  (страница https://tmbwg.github.io/turbmodels/Other_exp_Data/rotor37_exp.html).
- Только выдержки поиска: AGARD-AR-355 https://www.sto.nato.int/publications/AGARD/AGARD-AR-355/AGARD-AR-355.pdf ;
  TP-1338 https://ntrs.nasa.gov/citations/19790001889 ; TP-1659 https://ntrs.nasa.gov/citations/19800012840 ;
  CR-204134 https://ntrs.nasa.gov/citations/19970024917 ; CC3 — https://www.jafmonline.net/article_2322_e303ff6a8d63de9b462bc46b639af732.pdf ,
  https://ntrs.nasa.gov/citations/20260005056 ; SRV2 — https://kbwiki.ercoftac.org/w/index.php/AC_6-08 ,
  https://torroja.dmt.upm.es/congresos/asme_2011/data/pdfs/trk-33/GT2011-46505.pdf ,
  https://www.researchgate.net/publication/224796803_Investigations_of_the_Flow_Through_a_High_Pressure_Ratio_Centrifugal_Impeller
