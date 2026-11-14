---
name: su2-turbo
description: Аэродинамика ступеней трека «Простая литая ГТУ-заменитель» в SU2 — установка (сборка v8.6.0 «Harrier»), сетка одного межлопаточного канала из routes/gtu/platform/geom/impeller.py (Gmsh / Open CASCADE, формат .su2, маркеры), конфигурация сжимаемого RANS (SST/SA, вращающаяся система, периодичность, полные параметры на входе, статическое давление на выходе, Giles, mixing-plane для ступени), сходимость и контроль расхода, полная степень сжатия и изоэнтропический КПД, характеристика и граница устойчивости, проверка на эталонах NASA Rotor 37 / Stage 35, NASA CC3, DLR SRV2; каждый ключ сверен с исходником SU2 v8.6.0. Используй, когда готовишь или разбираешь расчёт SU2 колеса компрессора, турбины или ступени, строишь характеристику или сверяешь SU2 с эталоном. Триггеры — «SU2», «SU2_CFD», «расчёт колеса», «межлопаточный канал», «mixing plane», «характеристика компрессора», «граница помпажа», «Rotor 37», «CC3», «SRV2», «КПД ступени».
---

# SU2: ступени центробежного компрессора и турбины

Узел платформы трека (issue #59): одномерные размеры (`gtu_explore.comp_stage`, навык `gtu-cycle`) →
сетка канала (`platform/geom/impeller.py`) → SU2 → степень сжатия, КПД, характеристика. Состояние на
14.11.2026: **SU2 в контейнере не установлен, ни одного прогона SU2 не было.** Ключи и поведение ниже
сверены с исходником SU2 v8.6.0 и исходником документации (раздел «Источники»); пометка «(код)» — по
исходнику, «(док.)» — по документации, «(пример)» — по учебной конфигурации SU2, «не проверено» — как есть.

## 1. Что где лежит

| что | где |
|---|---|
| геометрия и сетка колеса | `routes/gtu/platform/geom/impeller.py`, описание и ограничения — `platform/geom/README.md` |
| вход случаев `core130`, `taurus60` | `gtu_explore_results.json` (схема 130), `gtu_dropin_results.json` (1Ц, класс N) |
| параметры для SU2 после сетки | `out/<case>_<level>/impeller_<case>_<level>.json` → `su2.periodic.angle_deg`, `su2.rotation_rate_rad_s`, `d1.p01`, `d1.T01`, `d1.PR`, `d1.eta_comp_stage` |
| облачные прогоны | `platform/aws/runner.py` (DryRun по умолчанию, живой запуск — `--live`), README там же |
| узел SU2 (конфигурации, обработка) | `routes/gtu/platform/su2/` — папки ещё нет; имя взято из примера задачи в `platform/aws/README.md` |

Сетка турбины в репозитории не строится: `impeller.py` делает только колесо компрессора.

## 2. Установка

- **apt и pip пакета нет.** Ubuntu 24.04: `apt-cache search --names-only su2` пусто (индексы noble есть,
  `gmsh` в них находится); PyPI: `su2` и `pysu2` — 404. Проверено 14.11.2026.
- **Сборка из исходников** (док. Build-SU2-Linux-MacOS; v8.6.0 «Harrier», тег от 30.10.2026, meson ≥ 1.8.2;
  `meson.py` сам скачивает meson, ninja и внешние библиотеки с GitHub):
  ```bash
  git clone --depth 1 -b v8.6.0 https://github.com/su2code/SU2.git && cd SU2
  ./meson.py setup build -Dwith-mpi=disabled --prefix=$HOME/su2   # MPI: -Dwith-mpi=enabled; OpenMP: -Dwith-omp=true
  ./ninja -C build -j4 install
  export SU2_RUN=$HOME/su2/bin PATH=$HOME/su2/bin:$PATH PYTHONPATH=$HOME/su2/bin:$PYTHONPATH
  ```
  В этой сессии сборка не запускалась (отклонена системой разрешений сессии): время сборки, нужные
  пакеты apt и путь `bin` под `--prefix` не проверены.
- **Готовые архивы:** страница загрузки берёт из последнего выпуска GitHub файлы `*linux64-omp.zip` и
  `*linux64-mpi.zip` (код `download.html`); есть ли они у v8.6.0 — не проверено (API GitHub из сессии закрыт).
- Запуск (док. Execution; код `SU2_CFD.cpp`): `SU2_CFD case.cfg`; MPI — `mpirun -n 4 SU2_CFD case.cfg`;
  сборка с OpenMP — `SU2_CFD -t 4 case.cfg`; `SU2_CFD -d case.cfg` — сухой прогон: только разбор
  конфигурации на фиктивной геометрии.

## 3. Сетка одного канала

```bash
cd routes/gtu/platform/geom
python impeller.py --case core130 --level coarse      # out/core130_coarse/passage_core130_coarse.su2
grep -n "NDIME\|NMARK\|MARKER_TAG" out/core130_coarse/passage_core130_coarse.su2
```

- Маркеры: `inlet outlet hub shroud blade periodic_1 periodic_2`, объём `fluid`, только тетраэдры. Gmsh
  пишет имя физической группы в `MARKER_TAG`, тетраэдр — код 10 (проверено записью куба в Gmsh 4.15.2).
- Ось машины — `z`: радиус для сечений по высоте SU2 считает от оси `z` (док., код). Сетка в метрах.
- **Периодичность — узел в узел** (код `CGeometry::MatchPeriodic`): точки первого маркера поворачиваются на
  заданный угол (в градусах) и ищется узел второго в пределах 1·10⁻⁶ (единицы сетки). Промах печатается
  «Bad match for point …», расчёт идёт дальше («Computation will continue, but be cautious»). В журнале
  должна быть строка «Matched N points with a max distance of: …» без «Bad match». У `impeller.py`
  `periodic_2` = поворот `periodic_1` на +360/Z вокруг +z, поэтому порядок — `(periodic_1, periodic_2, …)`.
- Пристенных призм нет (README узла geom): для RANS — пристенные функции или своя доработка сетки. В Gmsh
  `extrudeBoundaryLayer` есть только у встроенного ядра (`gmsh.model.geo`), у `gmsh.model.occ` такой функции
  нет (проверено `dir()`); как совместить с областью OCC из `impeller.py` — не проверено.
- Маркер `hub` идёт от носка входа через колесо до r = 1,3·r2 (`impeller.py`, `hub_poly`): во вращающейся
  системе вращается и стенка безлопаточного диффузора за r2. До разделения маркера в `impeller.py` (не
  сделано) числа по диффузору — с этой оговоркой.

## 4. Конфигурация: одна решётка (колесо)

Все ключи блока есть в `Common/src/CConfig.cpp` v8.6.0 (проверено). Численная схема — из учебных
конфигураций SU2 (`axial_stage_2D`, Aachen), под колесо не подбиралась. `<…>` — из JSON сетки (раздел 1).

```
SOLVER= RANS
KIND_TURB_MODEL= SST                       % или SA
MATH_PROBLEM= DIRECT
MACH_NUMBER= 0.05
FREESTREAM_OPTION= TEMPERATURE_FS
FREESTREAM_PRESSURE= <d1.p01>
FREESTREAM_TEMPERATURE= <d1.T01>
FREESTREAM_TURBULENCEINTENSITY= 0.03
FREESTREAM_TURB2LAMVISCRATIO= 100.0
INIT_OPTION= TD_CONDITIONS
REF_DIMENSIONALIZATION= DIMENSIONAL
FLUID_MODEL= IDEAL_GAS
GAMMA_VALUE= 1.4
GAS_CONSTANT= 287.058
VISCOSITY_MODEL= SUTHERLAND
MU_REF= 1.716E-5
MU_T_REF= 273.15
SUTHERLAND_CONSTANT= 110.4
CONDUCTIVITY_MODEL= CONSTANT_PRANDTL
GRID_MOVEMENT= ROTATING_FRAME
MOTION_ORIGIN= 0.0 0.0 0.0
ROTATION_RATE= <su2.rotation_rate_rad_s: 0 0 ω>
MACH_MOTION= 0.35
MARKER_HEATFLUX= (hub, 0.0, blade, 0.0, shroud, 0.0)
MARKER_SHROUD= (shroud)
MARKER_WALL_FUNCTIONS= (hub, STANDARD_WALL_FUNCTION, blade, STANDARD_WALL_FUNCTION, shroud, STANDARD_WALL_FUNCTION)
MARKER_PERIODIC= (periodic_1, periodic_2, 0.0, 0.0, 0.0, 0.0, 0.0, <su2.periodic.angle_deg>, 0.0, 0.0, 0.0)
INLET_TYPE= TOTAL_CONDITIONS
MARKER_INLET= (inlet, <T01>, <p01>, 0.0, 0.0, 1.0)
MARKER_OUTLET= (outlet, <статическое давление на выходе, Па>)
MARKER_ANALYZE= (inlet, outlet)
MARKER_ANALYZE_AVERAGE= MASSFLUX
NUM_METHOD_GRAD= WEIGHTED_LEAST_SQUARES
CFL_NUMBER= 2.0
LINEAR_SOLVER= FGMRES
LINEAR_SOLVER_PREC= LU_SGS
LINEAR_SOLVER_ERROR= 1E-4
LINEAR_SOLVER_ITER= 10
CONV_NUM_METHOD_FLOW= ROE
MUSCL_FLOW= YES
RAMP_MUSCL= YES
RAMP_MUSCL_COEFF= (0.0, 1.0, 20.0)
SLOPE_LIMITER_FLOW= VAN_ALBADA_EDGE
ENTROPY_FIX_COEFF= 0.01
CONV_NUM_METHOD_TURB= SCALAR_UPWIND
MUSCL_TURB= NO
ITER= <число итераций>
CONV_FIELD= RMS_DENSITY
CONV_RESIDUAL_MINVAL= <порог, lg>
SCREEN_OUTPUT= (INNER_ITER, RMS_DENSITY, RMS_ENERGY, SURFACE_MASSFLOW, SURFACE_TOTAL_PRESSURE)
HISTORY_OUTPUT= (ITER, RMS_RES, FLOW_COEFF, FLOW_COEFF_SURF)
VOLUME_OUTPUT= (COORDINATES, SOLUTION, PRIMITIVE)
OUTPUT_FILES= (RESTART, PARAVIEW)
TABULAR_FORMAT= CSV
MESH_FILENAME= passage_core130_coarse.su2
MESH_FORMAT= SU2
```

- `ROTATING_FRAME`: сетка неподвижна, `ROTATION_RATE` — вектор в рад/с, `MACH_MOTION` во всех примерах
  0,35 (док.). Источник вращения в коде — ω×ρu от консервативных переменных, то есть скорость в решении
  абсолютная (вывод из `flow_sources.cpp`; прямой строки в документации нет): полные параметры на входе и в
  `SURFACE_TOTAL_PRESSURE` — абсолютные. Знак ω — только из JSON сетки, не вручную.
- `MARKER_SHROUD`: в этих узлах скорость сетки обнуляется, покрышка не вращается (док., код); в
  `MARKER_HEATFLUX` она тоже нужна, как в Aachen (пример).
- `INIT_OPTION= TD_CONDITIONS`: умолчание `REYNOLDS` пересчитывает давление из числа Рейнольдса (док.).
- Модели: SST по умолчанию — `V2003m`; SA по умолчанию — SA-noft2; у SA есть `SA_OPTIONS= ROTATION` (док., код).
- **Вариант с турбомашинными маркерами** (нужен для Giles, `TURBO_PERF`, наклонов выхода):
  `MARKER_TURBOMACHINERY= (inlet, outlet)`, `TURBOMACHINERY_KIND= AXIAL_CENTRIFUGAL`,
  `TURBO_PERF_KIND= COMPRESSOR`, вместо `MARKER_INLET/OUTLET` —
  `MARKER_GILES= (inlet, TOTAL_CONDITIONS_PT, <p01>, <T01>, 1.0, 0.0, 0.0, 1.0, 0.0, outlet, STATIC_PRESSURE_1D, <p>, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0)`
  (9 значений: маркер, тип, два значения, три направления, релаксация среднего и Фурье — док.),
  `AVERAGE_PROCESS_KIND= MIXEDOUT`, `PERFORMANCE_AVERAGE_PROCESS_KIND= MIXEDOUT`, в `HISTORY_OUTPUT` — `TURBO_PERF`.
  Для тетраэдров — `SPANWISE_KIND= EQUISPACED` и `NUM_SPANWISE_SECTIONS`: `AUTOMATIC` считает узлы на
  стыке входа и периодической грани (док., пометка «to be reviewed»). Перед запуском — `mkdir TURBOMACHINERY`
  (иначе профили по высоте не пишутся, док.). Пуск на большой π и ω — `RAMP_OUTLET`, `RAMP_MOTION_FRAME` и
  их `*_COEFF` = (начало, шаг обновления, длина), задавать всегда (док.).
- **`TURBOMACHINERY_KIND` по коду** (`CPhysicalGeometry.cpp`, по какой координате мерится высота):
  колесо с осевым входом и радиальным выходом — `AXIAL_CENTRIFUGAL`; радиальный диффузор (вход и выход
  радиальные) — `CENTRIFUGAL`; рабочее колесо радиальной турбины — `CENTRIPETAL_AXIAL`; её сопловой
  аппарат — `CENTRIPETAL`; осевой венец — `AXIAL`.

## 5. Ступень: mixing plane

Образцы: `TestCases/turbomachinery/Aachen_turbine/aachen_3D_MP_restart.cfg` (3 венца, 3D) и
`APU_turbocharger/Jones_restart.cfg` (радиальная турбина: статор + ротор) в репозитории SU2.

- Главный файл: `MULTIZONE= YES`, `CONFIG_LIST= (zone_1.cfg, zone_2.cfg)` по потоку; `OUTER_ITER` вместо `ITER`.
  Файл зоны — только отличия: `GRID_MOVEMENT`, `ROTATION_RATE`, рампа (док., пример).
- Каждый стык — парой (верх по потоку, низ): `MARKER_ZONE_INTERFACE` и `MARKER_MIXINGPLANE_INTERFACE`, на
  них Giles `MIXING_OUT` и `MIXING_IN`; `MIXINGPLANE_INTERFACE_KIND= LINEAR_INTERPOLATION`,
  `TURBULENT_MIXINGPLANE= YES` (все примеры).
- `MARKER_PERIODIC` и `MARKER_TURBOMACHINERY` — по зонам по порядку: сечения зоны n ищутся по n-й
  периодической паре (код: `PeriodicBoundary == iZone + 1`). `TURBOMACHINERY_KIND` и `TURBO_PERF_KIND` —
  по одному значению на зону, иначе «Insufficient TURBO_PERF_KIND options» (док.).
- Поля выхода — с номером зоны: `RMS_DENSITY[1]`, `TotalTotalEfficiency[1]`; `TURBO_PERF` SU2 сам относит
  к последней зоне (код `COutput.cpp`). Тип расчёта КПД ступени берётся из `TURBO_PERF_KIND` зоны 0 (код).

## 6. Сходимость и расход

- Критерий останова — невязки (`CONV_FIELD= RMS_DENSITY`, `CONV_RESIDUAL_MINVAL`) или Коши по полям
  `MARKER_ANALYZE` (`CONV_CAUCHY_ELEMS`, `CONV_CAUCHY_EPS`); поля `TURBO_PERF` в `CONV_FIELD` не работают (док.).
- Расход — столбцы `Avg_Massflow(inlet)` и `Avg_Massflow(outlet)` в `history.csv`. Знаки противоположны:
  нормаль границы в SU2 смотрит внутрь области (код `CEulerSolver.cpp`). `Avg_Massflow` без скобок — сумма
  по маркерам со знаком, то есть дисбаланс; `Avg_TotalPress` без скобок — сумма давлений, не среднее (код
  `CFlowOutput.cpp`): брать только столбцы со скобками.
- Принято в этом навыке (не из документации SU2): точка сошлась, если невязки упали не меньше чем на
  4 порядка, дисбаланс расхода ≤ 0,5 %, π и η за последние 500 итераций меняются меньше 0,1 %.
- Giles и осреднённый Riemann работают только на самой мелкой сетке; неустойчиво с многосеточным — `MGLEVEL= 0` (док.).

## 7. Степень сжатия и КПД

Свой расчёт — из столбцов `history.csv` (осреднение по расходу, `MARKER_ANALYZE_AVERAGE= MASSFLUX` — умолчание, код):

```python
import sys   # python3 su2_pr_eta.py history.csv 1.4
g = float(sys.argv[2]) if len(sys.argv) > 2 else 1.4          # GAMMA_VALUE из конфигурации
lines = [l for l in open(sys.argv[1]) if l.strip()]
head = [h.strip().strip('"') for h in lines[0].split(",")]
row = dict(zip(head, (float(v) for v in lines[-1].split(","))))
pt1, pt2 = row["Avg_TotalPress(inlet)"], row["Avg_TotalPress(outlet)"]
tt1, tt2 = row["Avg_TotalTemp(inlet)"], row["Avg_TotalTemp(outlet)"]
m1, m2 = abs(row["Avg_Massflow(inlet)"]), abs(row["Avg_Massflow(outlet)"])
pr = pt2 / pt1
eta = (pr ** ((g - 1) / g) - 1) / (tt2 / tt1 - 1)
print(f"PR_tt {pr:.4f}  eta_tt {eta:.4f}  m_in {m1:.4f}  m_out {m2:.4f}  дисбаланс {abs(m1 - m2) / m1:.3%}")
```

Разбор проверен на искусственном файле в формате писателя SU2 (по коду: заголовки в кавычках, столбцы
выровнены пробелами); на настоящем `history.csv` — не проверен. Формула η — для идеального газа с
постоянным γ, как `IDEAL_GAS` в конфигурации. Приведённый расход для характеристики — ṁ·√(T₀₁/288,15)/(p₀₁/101 325).

- `TURBO_PERF` (код `CTurboOutput.cpp`): `PressureRatioTT`, `TotalTotalEfficiency` (в процентах); изоэнтропическая
  полная энтальпия там — h(s₁, p₂) + V₂²/2, а не h(s₁, p₀₂). В отчёт — оба КПД, каждый со своим определением.
- Сравнение с одномерной оценкой (`d1.PR`, `d1.eta_comp_stage`) — сравнение, а не проверка: расхождение
  записывается, одномерную модель под SU2 не подгонять.

## 8. Характеристика и граница устойчивости

- Линия частоты: давление на выходе от низкого (запирание: расход перестаёт расти) вверх; каждая точка — с
  решения предыдущей: `RESTART_SOL= YES`, `SOLUTION_FILENAME=` путь к её `restart` без расширения (SU2
  добавляет расширение, в мультизоне — и номер зоны; код `CConfig::GetFilename`).
- У вершины характеристики, где π почти не зависит от давления, — Giles `MASS_FLOW_OUTLET` (значение 1 —
  расход, значение 2 — релаксация; док.).
- Последняя сошедшаяся стационарная точка — граница устойчивости стационарного расчёта, а не помпаж: помпаж
  зависит от системы (объёмы, сеть) и требует нестационарного расчёта, которого в навыке нет. Так и подписывать.
- Другие частоты — пересчёт `ROTATION_RATE`; каждая точка — отдельная задача `runner.py`.

## 9. Проверка на открытых эталонах

Все числа таблицы — **предварительно**: из выдержек поиска, первоисточники (NTRS, AGARD, ERCOFTAC, DLR) из
контейнера не открылись. Арифметика U = π·D·n/60 сверена отдельно.

| эталон | что проверяет | числа (предварительно) | геометрия и данные |
|---|---|---|---|
| NASA Rotor 37 | одиночный трансзвуковой венец, ROTATING_FRAME, ударные волны, зазор | 36 лопаток; 17 188,7 об/мин; 454,14 м/с на периферии (D ≈ 0,505 м по арифметике); π 2,106; 20,19 кг/с; зазор 0,356 мм (AGARD-AR-355); в NASA TP-1659 проектная π 2,05 | NASA TMR: страница Rotor 37 (исходник в github.com/tmbwg/turbmodels) даёт архив с геометрией и данными на nasa.gov — не скачан |
| NASA Stage 35 | mixing plane на осевой ступени | проектная π 1,82; пиковые КПД ротора и ступени 0,872 и 0,845 при π 1,875 и 1,842; запас устойчивости 21,8 % (TP-1338) | NASA TP-1338 — не открыт |
| NASA CC3 | центробежное колесо с полулопатками, ближе всего к `core130` | 4,54 кг/с; π 4:1; 21 789 об/мин; 15 + 15 лопаток; загиб 50°; лопаточный диффузор 24 лопатки | координаты — NASA CR-204134 (McKain, Holbrook, 1997) — не открыт; NASA HECC (2026) пишет о расхождении изготовленного колеса с проектом |
| DLR SRV2-O | трансзвуковой вход колеса, ближе всего к `taurus60` (M₁s отн. 1,41) | 13 + 13 лопаток; 50 000 об/мин; 2,55 кг/с; π проектная 6,1–6,2, измеренная ступени 5,7; 586 м/с; r2 112 мм (π·0,224·50 000/60 = 586,4 м/с) | ERCOFTAC KB AC 6-08 — не открыт; доступность геометрии не проверена |

Порядок: Rotor 37 (открыты и геометрия, и данные) — отладка постановки; затем CC3 и SRV2-O. Допуски
согласия документация не задаёт: в PR записать отклонение π, η и расхода в расчётной точке и форму линии
частоты против опыта, без подгонки констант модели под эталон.

## 10. Обязательные правила

- Ключи — только из `Common/src/CConfig.cpp` той же версии; `config_template.cfg` неполон (нет, например,
  `MULTIZONE`). Неизвестный ключ останавливает SU2: «invalid option name. Check current SU2 options in
  config_template.cfg» (код). Проверка до запуска (из папки с конфигурацией):
  ```bash
  curl -sS https://raw.githubusercontent.com/su2code/SU2/v8.6.0/Common/src/CConfig.cpp -o /tmp/CConfig.cpp
  grep -oE '^[A-Z][A-Z0-9_]*' case.cfg | sort -u | while read k; do grep -q "\"$k\"" /tmp/CConfig.cpp || echo "нет ключа: $k"; done
  ```
  Проверяет имена, не значения; значения — `SU2_CFD -d case.cfg`.
- Публикуемое число — из results-JSON рядом со скриптом (как в навыке `gtu-cycle`), с конфигурацией,
  журналом и версией SU2; «SU2 показал» без журнала прогона не писать.
- Локально — не дольше 20 минут; дольше — `runner.py` (us-east-1, спот, потолок 300 долл., запись только
  в `s3://$ALPHA_O_S3/<префикс>/gtu/…`), живой запуск только флагом `--live`. Ключей в файлах нет.

## 11. Типовые ошибки

- Скопированы `TestCases/turbomachinery/centrifugal_blade` и `centrifugal_stage`: там удалённые ключи
  `MARKER_NRBC`, `GRID_MOVEMENT_KIND`, `ROTATION_RATE_Z`, `MARKER_TURBO_PERFORMANCE`, `MIXING_PROCESS_TYPE`,
  `MARKER_MIXINGPLANE` (проверено проверкой выше); в регрессионных тестах их нет.
- Старые имена рамп `RAMP_OUTLET_PRESSURE`, `RAMP_ROTATION_FRAME` — ошибка «deprectaed» (код); нужны `RAMP_OUTLET`,
  `RAMP_MOTION_FRAME`. Рампа давления с выходом Riemann — «Outlet pressure ramp only implemented for NRBC» (док.).
- Перепутан порядок: `MARKER_INLET` — T₀, затем p₀; Giles и Riemann `TOTAL_CONDITIONS_PT` — p₀, затем T₀.
- Пары `MARKER_PERIODIC` в обратном порядке или угол со знаком «минус»: SU2 не останавливается, только
  печатает «Bad match» — читать журнал.
- `TURBOMACHINERY_KIND= CENTRIFUGAL` у колеса с осевым входом: высота на входе мерится по `z` — нужно `AXIAL_CENTRIFUGAL`.
- Нет `MARKER_TURBOMACHINERY` при `MARKER_GILES` — «Giles Boundary conditions can only be used with
  turbomachinery markers» (док.).
- y⁺ не проверен: сетка без призм. Поле `Y_PLUS` (группа `PRIMITIVE`) — в `VOLUME_OUTPUT`.

## 12. Проверка результата

1. Журнал: «Matched … max distance» без «Bad match»; маркеры и число узлов как в JSON сетки.
2. Невязки и дисбаланс расхода — по разделу 6; расход на входе и выходе по модулю совпадают.
3. π и η посчитаны двумя способами (раздел 7: свой и `TURBO_PERF`), расхождение объяснено определением.
4. Сетка: coarse → medium (`impeller.py --level medium`) — изменение π и η записано; на одном уровне вывода нет.
5. До SU2 на колесах трека — Rotor 37 (раздел 9) с записанным отклонением от опыта.

## Источники (сверено 14.11.2026)

- Исходник SU2 v8.6.0 «Harrier» (клон тега, коммит от 30.10.2026): https://github.com/su2code/SU2/tree/v8.6.0 —
  `Common/src/CConfig.cpp` (опции, «invalid option name», устаревшие рампы), `Common/include/option_structure.hpp`
  и `.inl` (значения перечислений, `MARKER_PERIODIC` в градусах), `Common/src/geometry/CGeometry.cpp`
  (`MatchPeriodic`), `Common/src/geometry/CPhysicalGeometry.cpp` (высота по `TURBOMACHINERY_KIND`),
  `SU2_CFD/src/output/CTurboOutput.cpp`, `CFlowOutput.cpp`, `COutput.cpp`, `include/output/COutput.hpp`,
  `SU2_CFD/src/numerics/flow/flow_sources.cpp`, `SU2_CFD/src/solvers/CEulerSolver.cpp`, `SU2_CFD/src/SU2_CFD.cpp`,
  `config_template.cfg`, `meson.build`, `meson_options.txt`, `TestCases/turbomachinery/*`, `TestCases/serial_regression.py`.
- Документация SU2 — исходник сайта https://github.com/su2code/su2code.github.io (коммит 67d0f94 от 11.11.2026);
  сам https://su2code.github.io из контейнера не открылся (DNS). Страницы:
  https://su2code.github.io/docs_v7/Setup-Guide-Turbomachinery/ ,
  https://su2code.github.io/docs_v7/Markers-and-BC/ , https://su2code.github.io/docs_v7/Build-SU2-Linux-MacOS/ ,
  https://su2code.github.io/docs_v7/SU2-Linux-MacOS/ , https://su2code.github.io/docs_v7/Execution/ ,
  https://su2code.github.io/docs_v7/Turbulence-Models/ , https://su2code.github.io/tutorials/Aachen_Turbine/ ,
  https://su2code.github.io/download.html
- Сетки учебных случаев: https://github.com/su2code/TestCases (список файлов `turbomachinery/`).
- Gmsh 4.15.2: справка API установленного пакета (`gmsh.model.geo.extrudeBoundaryLayer`), запись `.su2` проверена.
- NASA TMR, Rotor 37: исходник https://github.com/tmbwg/turbmodels/blob/main/Other_exp_Data/rotor37_exp.html
  (страница https://tmbwg.github.io/turbmodels/Other_exp_Data/rotor37_exp.html).
- Только выдержки поиска (не открывались): AGARD-AR-355 https://www.sto.nato.int/publications/AGARD/AGARD-AR-355/AGARD-AR-355.pdf ;
  NASA TP-1338 https://ntrs.nasa.gov/citations/19790001889 ; NASA TP-1659 https://ntrs.nasa.gov/citations/19800012840 ;
  NASA CR-204134 https://ntrs.nasa.gov/citations/19970024917 ; CC3 — https://www.jafmonline.net/article_2322_e303ff6a8d63de9b462bc46b639af732.pdf ,
  https://ntrs.nasa.gov/citations/20260005056 ; SRV2 — https://kbwiki.ercoftac.org/w/index.php/AC_6-08 ,
  https://torroja.dmt.upm.es/congresos/asme_2011/data/pdfs/trk-33/GT2011-46505.pdf ,
  https://www.researchgate.net/publication/224796803_Investigations_of_the_Flow_Through_a_High_Pressure_Ratio_Centrifugal_Impeller
