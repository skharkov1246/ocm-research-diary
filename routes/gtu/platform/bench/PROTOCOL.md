# Протокол сравнения CPU и GPU: один межлопаточный канал

Пункт 5 задания «Расчёт промышленной газовой турбины на AWS» (14.11.2026): по каким измерениям решать,
переносить ли вычислительное ядро на CUDA и когда переходить к нескольким узлам (Hpc8a + EFA), с проверкой
точности и сходимости каждой конфигурации. Трек «Простая литая ГТУ-заменитель», issue #59.

**Состояние на 14.11.2026.** Прогонов по протоколу не было; в облаке ничего не запускалось, AWS не
вызывался. Числа ниже — из `break_even.json`, `tariffs_us-east-1.json`, `mesh_stats.json` (скрипт
`break_even.py`, входные файлы с sha256) и из исходников решателей. Измеренных ускорений нет — их даст
этот протокол. Бенчмарк заблокирован, пока на CPU нет сошедшегося расчёта канала (README узла `su2`,
раздел «Дымовой прогон»): сравнивать время до сходимости, которой нет, нельзя.

| файл | что это |
|---|---|
| `PROTOCOL.md` | этот протокол |
| `bench.toml` | одна задача, критерии сходимости, допуски, правила решения, 8 конфигураций |
| `bench_parse.py` | разбор журналов SU2 по конфигурациям → `bench_summary.csv`, `bench_compare.csv` |
| `break_even.py` | тарифы (`tariffs`), сетки (`mesh`), порог окупаемости (`compute`) → `break_even.json` |
| `test_bench.py` | 17 тестов на придуманных журналах (SYNTHETIC) |

## 1. Что умеют наши решатели на GPU (проверено по первоисточникам)

### 1.1. SU2 v8.6.0 «Harrier»

Проверено по клону тега `v8.6.0` (коммит 7048fd3) и исходнику документации su2code.github.io (коммит
67d0f94 от 11.11.2026).

| что | как в v8.6.0 | где в исходнике |
|---|---|---|
| сборка | `-Denable-cuda=true`; `-Dcuda-arch=auto` определяет GPU машины через `nvidia-smi`, иначе указать (например `sm_89`); `-Dwith-omp=true` обязателен по смыслу (см. ниже) | `meson_options.txt`, `meson.build` (строки 19–60, 269–272) |
| включение | `ENABLE_CUDA= YES` (по умолчанию `NO`) | `Common/src/CConfig.cpp:1240` |
| что идёт на GPU | весь линейный решатель: векторы Крылова, произведение блочной матрицы на вектор (LDU), предобуславливатели Jacobi, ILU(0) (раскраска + итерационная факторизация, `LINEAR_SOLVER_ILU_GPU_SWEEPS`, по умолчанию 2; треугольные решения по уровням), LU-SGS и их варианты с int8-хранением вне диагонали (`Q_JACOBI`, `Q_LU_SGS`, `Q_IDENTITY`) | `Common/src/linear_algebra/CSysMatrixGPU.cu`, `CSysVectorGPU.cu`, `CSysMatrix.cpp:221–472` |
| что остаётся на CPU | сборка якобиана и невязок, градиенты, ограничители, граничные условия (Giles), турбомашинные осреднения, вывод: в `SU2_CFD/src` CUDA упоминается только в `CNewtonIntegration.cpp` | поиск `cuda`/`gpu` по `SU2_CFD/` |
| обмен CPU–GPU | матрица идёт на GPU один раз за линейное решение (конструктор `CSysMatrixVectorProduct` → `HtDTransfer`, копия `d`, `l`, `u`), правая часть и решение — туда, решение — обратно; предобуславливатели LINELET и PaStiX выполняются на хосте с пересылкой векторов при каждом применении | `CMatrixVectorProduct.hpp:91`, `CSysSolve.hpp:255–272`, `CPreconditioner.hpp:40–90` |
| ограничения | только один процесс MPI: «ENABLE_CUDA= YES is not supported with more than one MPI rank … Use a single rank with OpenMP threads» — хост-часть идёт потоками OpenMP; только прямой решатель `SU2_CFD` (не сопряжённый AD); блоки матрицы квадратные | `CConfig.cpp:3712–3733`, `CSysMatrix.cpp:298–422` |
| точность | `-Denable-mixedprec=true` — линейный решатель и предобуславливатель в FP32 (`su2mixedfloat = float`), остальное в FP64; документация: «for fluid simulations this does not reduce accuracy since the solution is iterative» | `Common/include/code_config.hpp:185–197`, `Build-SU2-Linux-MacOS.md` |
| документация | «FGMRES on GPUs with CUDA: build with `-Denable-cuda=true` and set `ENABLE_CUDA= YES` (v8.3.0)»; в коде v8.6.0 векторы переводятся на GPU для всех решателей `CSysSolve::Solve` — на GPU проверять FGMRES | `_docs_v7/Guide-to-v8.md` |
| профилирование | Tracy (`-Dtracy_enable=true`, зоны `SU2_ZONE_SCOPED` в 66 файлах, в том числе `HtDTransfer` и `MatrixVectorProductGPU`), gprof (`-Denable-gprof=true`); отдельного вывода «время сборки / время решения» нет: в истории — `ITER_TIME` и `WALL_TIME` (обе колонки называются `Time(sec)`), в сводке `WRT_PERFORMANCE= YES` — предобработка, счёт, вывод | `COutput.cpp:2196–2200`, `CDriver.cpp:488–520`, `_docs_v7/Tracy-Integration.md` |

Вывод для протокола. «Перенос ядра на CUDA» для SU2 — это не своя разработка, а сборка с CUDA и
`ENABLE_CUDA= YES`: ядро = линейный решатель. Сборка якобиана остаётся на CPU GPU-машины, причём в
одном процессе с потоками OpenMP (на p5.4xlarge — 16 vCPU, на g6e.2xlarge — 8). Поэтому решают три
измерения: доля линейного решателя во времени CPU (f), ускорение самого решателя на GPU (s_k) и замедление
хост-части на CPU GPU-машины (r_host). Своя CUDA-реализация сборки — отдельное решение (раздел 6.5).

### 1.2. OpenFOAM

- **openfoam.org (OpenFOAM-dev, коммит fed424ac от 13.11.2026):** в `src/` и `applications/` нет ни одного
  файла с `cuda`; линейные решатели ядра — `GAMG`, `PBiCG`, `PBiCGStab`, `PCG`, `smoothSolver`, `diagonalSolver`
  (только CPU). Проверено клоном.
- **openfoam.com (ESI):** GPU — только внешними модулями. `petsc4Foam` (модуль `modules/external-solver`)
  передаёт линейную систему в PETSc; через PETSc доступны его бэкенды CUDA и предобуславливатель AmgX
  (`PCAMGX` — «Interface to NVIDIA's AmgX algebraic multigrid», исходник PETSc `src/ksp/pc/impls/amgx/amgx.cxx`,
  ветка main). По материалам NVIDIA на конференции ESI 2021 (выдержка поиска): AmgX через преобразование
  матрицы LDU → CSR (`foam2csr`) при каждом вызове решателя, ускорение решателя давления до 9 раз на A100
  против GAMG-PCG на узле CPU — **предварительно**. Сборка матрицы остаётся на CPU. Сам сайт
  develop.openfoam.com из контейнера недоступен (DNS не разрешается) — README модуля и текущий статус GPU не
  сверены; чтобы сверить, хост develop.openfoam.com нужно разрешить в настройках сети сессии.
- **OGL** (обёртка Ginkgo, github.com/hpsim/OGL): README — «PROOF OF CONCEPT … to provide GPGPU capabilities
  to OpenFOAM», исполнитель `cuda`/`hip`/`omp`, параметр `updateRHS` — копирование матрицы на устройство при
  каждом вызове решателя, значки версий v2212 и 10.
- **AmgX** (github.com/NVIDIA/AMGX, README): fp32, fp64 и смешанная точность, блочные системы, MPI, AMG
  классический и агрегационный, Крылов CG/BiCGSTAB/GMRES.
- Для нашей камеры сгорания (навык `openfoam-combustor`: `reactingFoam`, `chtMultiRegionFoam`) время обычно
  уходит не только в линейный решатель, но и в химию — это оценка, а не замер: долю f мерить тем же
  протоколом (раздел 4.2) до любых затрат на GPU.

### 1.3. CalculiX

- **Официальный репозиторий** (github.com/Dhondtguido/CalculiX, коммит 57a4e71 от 14.11.2026): решатель
  PaStiX умеет GPU — переменная окружения `PASTIX_GPU=1`, планировщик PaRSEC, смешанная точность по
  умолчанию (`PASTIX_MIXED_PRECISION`, при неудаче — возврат к двойной; `src/pastix.c`). GPU работает только
  в совмещённом вызове `pastix_main` (статика `linstatic`, нелинейные шаги `nonlingeo`, `electromagnetics`,
  устойчивость `arpackbu`); в отдельной факторизации `pastix_factor_main*` (частоты `arpack`, динамика `dyna`
  и др.) GPU принудительно выключен. Сборка: CMake `-DPASTIX=ON` (PaStiX v6.4.0, «supporting PaStiX builds with
  CUDA») или `-DPASTIX4CALCULIX=ON` (`INSTALL.md`). Руководство (`doc/CalculiX.tex`): PaStiX «really fast and
  can use the GPU … an acceleration of a factor between 3 and 8 compared to PARDISO» — утверждение автора,
  не наш замер.
- **Пакет Ubuntu 24.04 `calculix-ccx` 2.21-1** (установлен в контейнере): собран только со SPOOLES —
  в бинарнике строки «the PARDISO / PASTIX / SGI / TAUCS library is not linked».
- Наши задачи прочности (осесимметричный диск, сектор колеса — навык `calculix-rotor`) малы, а частоты
  (диаграмма Кэмпбелла) идут через `arpack` без GPU: GPU для CalculiX в протокол не входит. Вернуться, если
  сектор с контактом и ползучестью даст статику, где SPOOLES упрётся в память или время.

### 1.4. FP64 на H100 и L40S

| | H100 (p5.4xlarge) | L40S (g6e.2xlarge) | источник |
|---|---:|---:|---|
| память | 80 ГБ HBM3 | 48 ГБ GDDR6 | Price List AWS 14.11.2026 (поле `GPU Memory`) |
| полоса памяти | 3,35 ТБ/с | 0,864 ТБ/с | паспорта NVIDIA (выдержки поиска; www.nvidia.com из контейнера недоступен) |
| FP64 | 34 TFLOPS (Tensor Core — 67) | ≈ 1,43 TFLOPS | H100 SXM — страница и паспорт NVIDIA; L40S — FP32 91,6 TFLOPS / 64 |
| FP32 | — | 91,6 TFLOPS | паспорт L40S |

- В паспорте L40S числа FP64 нет; 1,43 TFLOPS — расчёт по белой книге Ada: «The FP64 TFLOP rate is 1/64th the
  TFLOP rate of FP32 operations», два блока FP64 на SM «to ensure any programs with FP64 code operate
  correctly» (выдержки поиска). Ранний обзор Hopper давал 30 и 60 TFLOPS FP64 — в расчёт взяты числа паспорта.
- Что p5 несёт именно H100 SXM, — вывод из «80 GB HBM3» в прайсе (у вариантов PCIe и NVL другая память), по
  странице AWS P5 не сверено (aws.amazon.com недоступен).
- **Что это значит для SU2 v8.6.0.** На GPU идёт только линейная алгебра, а она упирается в полосу памяти:
  интенсивность блочного SpMV на наших сетках — 0,20–0,24 флоп/байт в FP64 и 0,40–0,47 в FP32
  (`break_even.json`, `spmv_roofline`). На полосе L40S FP64 нужно лишь 205 GFLOPS из 1 431 доступных, на H100 —
  796 из 34 000: слабый FP64 у L40S для SpMV не ограничение, ограничение — полоса (у H100 в 3,9 раза шире).
  FP64 L40S может сказаться в факторизации ILU и обращении диагональных блоков (кубические по nvar операции) —
  это проверяет конфигурация `l40s_mixed` против `h100_mixed`. Полоса на доллар: H100 по требованию —
  0,49 ТБ/с на долл./ч, L40S по требованию — 0,39, L40S на споте (оценка) — 0,84.

## 2. Задача бенчмарка: одна сетка, одни условия, одни критерии

1. **Сетка** — одна, канал расчётной точки KV-3 (Taurus 60 генераторный, frame3), фиксируется sha256
   (`case.mesh_sha256`). Пока узел geom выдаёт только тетраэдры без пристенного слоя, а сошедшегося расчёта
   нет; пробный кандидат — `taurus60_coarse` (165 236 узлов), для слабого масштабирования —
   `core130_coarse` / `core130_medium` (139 121 / 584 357 узлов, `mesh_stats.json`).
2. **Начальное поле** — общий рестарт после разгона на CPU (`case.restart`, sha256): разгон (рампы частоты и
   давления, `run_case.py`) в сравнение не входит, он неустойчив и к GPU отношения не имеет.
3. **Один файл конфигурации SU2** на все конфигурации. Разрешено различаться только: `ENABLE_CUDA`, число
   процессов MPI и потоков OpenMP (`SU2_CFD -t N`), сборке (`fp64` / `mixed`). Одинаковы: граничные условия
   (Giles, p₀₁, T₀₁, p вых), модель (SST), схема (Roe + MUSCL), `CFL_NUMBER`/`CFL_ADAPT`, `LINEAR_SOLVER= FGMRES`,
   `LINEAR_SOLVER_PREC= ILU`, `LINEAR_SOLVER_ERROR`, `LINEAR_SOLVER_ITER`, `ITER`, `CONV_FIELD`, `CONV_RESIDUAL_MINVAL`.
   ILU на GPU строится иначе (итерационная факторизация по раскраске), а у MPI ILU блочный по процессам, —
   поэтому число линейных итераций разное по построению и записывается, а не подгоняется.
4. **Вывод:** `HISTORY_OUTPUT= (ITER, WALL_TIME, RMS_RES, LINSOL, FLOW_COEFF, FLOW_COEFF_SURF, TURBO_PERF)`
   (в шаблоне узла su2 нет `LINSOL` — добавить), `HISTORY_WRT_FREQ_INNER= 1`, `TABULAR_FORMAT= CSV`,
   `WRT_PERFORMANCE= YES`, `MARKER_ANALYZE= (inlet, outlet)`; рестарт — каждые N итераций (`OUTPUT_WRT_FREQ`),
   чтобы прерывание спота не стоило прогона. Все ключи сверены с `CConfig.cpp` v8.6.0.
5. **Сходимость** (`[convergence]`, из навыка `su2-turbo`, раздел 6): `rms[Rho]` упала не меньше чем на 4
   порядка от максимума первых 10 строк; |ṁвх + ṁвых| / ṁвх ≤ 0,5 %; π и η за последние 500 итераций меняются
   меньше чем на 0,1 %. Итерация сходимости — первая, где выполнено всё; время до сходимости — по истории.

## 3. Конфигурации

| id | машина | рынок | сборка | MPI × OpenMP | зачем |
|---|---|---|---|---|---|
| `cpu_c8i_mpi` | c8i.8xlarge, 32 vCPU | спот | CPU FP64 | 32 × 1 | эталон задания (точность и время) |
| `cpu_c8i_omp` | c8i.8xlarge | спот | CPU FP64 | 1 × 32 | отделяет «OpenMP вместо MPI» от «GPU» |
| `cpu_c7a_mpi` | c7a.16xlarge, 64 vCPU | спот | CPU FP64 | 64 × 1 | стандартная машина `runner.py` |
| `h100_host` | p5.4xlarge, 16 vCPU | по требованию | CUDA, `ENABLE_CUDA= NO` | 1 × 16 | r_host и f на CPU GPU-машины |
| `h100_fp64` | p5.4xlarge | по требованию | CUDA FP64 | 1 × 16 | ускорение в двойной точности |
| `h100_mixed` | p5.4xlarge | по требованию | CUDA mixed | 1 × 16 | цена смешанной точности по точности |
| `l40s_host` | g6e.2xlarge, 8 vCPU | по требованию | CUDA mixed, `ENABLE_CUDA= NO` | 1 × 8 | r_host на 8 vCPU |
| `l40s_mixed` | g6e.2xlarge | по требованию | CUDA mixed | 1 × 8 | L40S — только если mixed прошёл на H100 и задача влезает в 48 ГБ |

На c8i.8xlarge дополнительно (один раз) — 16 × 1 против 32 × 1: берётся более быстрое. Первый GPU-прогон —
по требованию с заданной длительностью (задание, п. 5); CPU — спот (решение владельца). Необязательно:
`h100_q` с `LINEAR_SOLVER_PREC= Q_LU_SGS` (int8 вне диагонали уменьшает трафик памяти) — предобуславливатель
другой, поэтому только с полной проверкой точности.

## 4. Что мерить на каждой конфигурации

| величина | как | файл в папке прогона | колонка `bench_summary.csv` |
|---|---|---|---|
| время на итерацию | `ITER_TIME` истории, медиана второй половины; среднее = время / итерации | `history*.csv` | `s_per_iter_median`, `s_per_iter_mean` |
| время по истории (два способа) | Σ `ITER_TIME` и `WALL_TIME`·(итерация + 1) — должны совпасть до 0,1 % | `history*.csv` | `time_hist_sum_s`, `time_hist_avg_s`, `time_check_rel` |
| итерации и время до сходимости | первая строка, где выполнены критерии раздела 2 | `history*.csv` | `iter_conv`, `time_to_conv_s` |
| полное время, предобработка, счёт, вывод | сводка `WRT_PERFORMANCE= YES` | `*.log` | `wall_total_s`, `preproc_s`, `compute_s`, `output_s` |
| сборка системы и решение | профилировщик (4.2) | `profile.toml` | `t_assembly_s`, `t_linsolve_s`, `f_linsolve` |
| линейные итерации | `LINSOL` в истории | `history*.csv` | `linsol_iter_mean` |
| пиковая память CPU | `/usr/bin/time -v` на каждый процесс: «Maximum resident set size»; сумма по процессам — верхняя граница | `*.log`, `time_v*.txt` | `cpu_peak_rss_sum_GiB`, `cpu_peak_rss_max_GiB` |
| пиковая память GPU, загрузка | `nvidia-smi --query-gpu=timestamp,memory.used,utilization.gpu --format=csv -lms 500` в фоне | `gpu_mem*.csv` | `gpu_peak_mem_GiB`, `gpu_util_mean_pct` |
| обмен CPU–GPU | Nsight Systems: `nsys profile -t cuda,osrt -o run SU2_CFD -t N case.cfg`, затем `nsys stats -r cuda_gpu_mem_time_sum -f csv run.nsys-rep` (имя отчёта сверить с `nsys stats --help-reports` установленной версии) | `nsys_mem*.csv` | `htod_s`, `dtoh_s` |
| стоимость | часы × (тариф машины + диск gp3 + публичный IPv4); время — до сходимости и машинное (запуск → гашение, журнал `runner.py`) | `machine.toml` | `tariff_usd_h`, `extras_usd_h`, `cost_solver_usd`, `cost_machine_usd` |
| точность | расход, момент, π, η, дисбаланс, невязки — раздел 5 | `history*.csv` | `mdot_*`, `torque_*`, `pr_tt`, `eta_tt`, `res_drop_orders` |

### 4.1. Повторы и шум

Решающие конфигурации (лучшая CPU и лучшая GPU) — по 2 прогона, время — медиана. Машина — целиком своя
(не общая, как контейнер, где 4 ядра делились с чужими задачами: там 0,5–2,3 с на итерацию на coarse —
ориентир, не замер). Сборка SU2 и установка пакетов — в образе, не в счёте прогона.

### 4.2. Доля линейного решателя f (два способа)

SU2 не выводит отдельно время сборки и решения, поэтому:

1. **perf** (без пересборки): `perf record -F 499 -g -o perf.$OMPI_COMM_WORLD_RANK.data SU2_CFD case.cfg` под
   `mpirun`; затем `perf report --children --stdio -i perf.0.data` — включительная доля символов
   `CSysSolve<…>::Solve` (решение с построением предобуславливателя) против `C*Solver::*Residual`,
   `C*Solver::*Iteration` (сборка). Команды на SU2 в этой сессии не запускались.
2. **Tracy** (`-Dtracy_enable=true`, зоны уже расставлены): на машине без экрана — запись `tracy-capture`
   и выгрузка `tracy-csvexport` из `subprojects/tracy` (утилиты есть в проекте Tracy; на SU2 не проверено).

Расхождение двух способов больше 5 п. п. — разбор, не публикация. Проверка третьим способом после GPU-прогона:
T_GPU должно сходиться с (1 − f)·r_host·T_CPU + время решателя на GPU (из `nsys` — `cuda_gpu_kern_sum`).

### 4.3. Видеопамять и обмен: оценка до аренды

Нижняя граница памяти линейных систем SU2 на GPU (`break_even.json`, `gpu_memory_estimate`): матрица LDU
(1 + 2E/N блоков nvar² на узел), её индексы, копия ILU(0) на устройстве и 22 + 4 вектора FGMRES, для систем
течения (nvar = 5) и SST (nvar = 2). Рёбра посчитаны по сеткам узла geom: 6,0–6,3 на узел.

| сетка | узлов | FP64, байт/узел | mixed, байт/узел | FP64 на GPU, ГБ | на GPU за итерацию (FP64), ГБ | макс. узлов H100 / L40S (FP64) |
|---|---:|---:|---:|---:|---:|---:|
| core130 coarse | 139 121 | 7 709 | 3 966 | 1,07 | 0,435 | 8,8 млн / 5,3 млн |
| core130 medium | 584 357 | 8 031 | 4 133 | 4,69 | 1,918 | 8,5 млн / 5,1 млн |
| taurus60 coarse | 165 236 | 7 805 | 4 016 | 1,29 | 0,524 | 8,7 млн / 5,2 млн |

Это нижняя граница (без контекста CUDA, рабочих массивов и предобуславливателя на хосте): по ней канал
влезает в любую из двух GPU с запасом; правило решения — по замеру (`gpu_peak_mem_GiB`) с экстраполяцией по
числу узлов (байт/узел из двух сеток).

## 5. Проверка точности и сходимости каждой конфигурации

1. Каждая конфигурация сама проходит критерии раздела 2 (`converged`) и не расходится (NaN → `diverged`);
   в журнале — «Exit Success», «Matched … points» без «Bad match».
2. Против эталона `cpu_c8i_mpi` (`bench_compare.csv`, допуски `[tolerances]`):

| величина | допуск | как считается |
|---|---|---|
| расход ṁвх | 0,1 % отн. | `Avg_Massflow(inlet)`, канал |
| момент колеса | 0,2 % отн. | ṁ·cp·(T₀₂ − T₀₁)/ω × Z (адиабатные стенки) |
| полная степень сжатия π | 0,1 % отн. | `Avg_TotalPress(outlet)/(inlet)`, средние по расходу |
| изоэнтропический КПД η | 0,1 п. п. абс. | (π^((γ−1)/γ) − 1)/(τ − 1) |
| дисбаланс расхода | 0,1 п. п. абс. | (ṁвх + ṁвых)/ṁвх |
| итерации до сходимости | не больше чем в 1,5 раза | иначе сходимость другая — разбор, время не сравнивается |

   Основание допусков: разбиение сетки MPI на 4 и 2 процесса уже дало 6,1·10⁻⁵ по расходу (README узла su2);
   допуски в 16–33 раза шире этого шума и много уже погрешности отбора схем (±2 п. п. по КПД).
3. Внутренние проверки разбора: время двумя способами (≤ 0,1 %), момент двумя способами — по энтальпии и
   по `EulerianWork` из `TURBO_PERF` (`torque_methods_rel`); на пробном прогоне SU2 core130 (не сошедшемся)
   они разошлись на 0,08 % при дисбалансе 1,4 % и на 21 % при дисбалансе 88 % — расхождение методов само
   показывает несошедшийся расчёт.
4. `h100_mixed` против `h100_fp64` — отдельно: если смешанная точность не проходит допуски, L40S не
   испытывается (по заданию — только «если допустима смешанная точность»).

## 6. Правило решения: переносить ли ядро на CUDA

### 6.1. Окупаемость тарифа

Прогон на GPU дешевле, если ускорение по времени до сходимости S = T_CPU / T_GPU больше порога
S* = тариф GPU / тариф CPU. Принято: переносить при S ≥ 1,2·S* (запас на разброс времени и цены спота).
Сравнение — с самой дешёвой по стоимости прогона CPU-конфигурацией, а не с одной машиной.

Тарифы us-east-1 (Price List AWS, версия 20261114012926; спот — оценка по Spot Instance Advisor):

| машина | по требованию, долл./ч | спот (оценка), долл./ч | экономия спота | прерывания (советник) |
|---|---:|---:|---:|---|
| c8i.8xlarge | 1,49936 | 0,58475 | 61 % | > 20 % |
| c7a.16xlarge | 3,28448 | 1,47802 | 55 % | 10–15 % |
| c7i.16xlarge | 2,856 | 0,94248 | 67 % | > 20 % |
| g6e.2xlarge (L40S) | 2,24208 | 1,03136 | 54 % | 5–10 % |
| p5.4xlarge (H100) | 6,88 | 6,2608 | 9 % | < 5 % |

Порог S* (`break_even.json`, `break_even`):

| GPU \ CPU | c8i.8xlarge | c7a.16xlarge | c7i.16xlarge |
|---|---:|---:|---:|
| H100 по требованию / CPU по требованию | 4,59 | 2,09 | 2,41 |
| H100 спот / CPU спот | 10,71 | 4,24 | 6,64 |
| **H100 по требованию / CPU спот** | **11,77** | **4,65** | **7,30** |
| L40S по требованию / CPU по требованию | 1,50 | 0,68 | 0,79 |
| L40S спот / CPU спот | 1,76 | 0,70 | 1,09 |
| **L40S по требованию / CPU спот** | **3,83** | **1,52** | **2,38** |

- **Сверка с заданием.** Цифры задания (4,6 для H100 и 1,5 для L40S) — eu-north-1, по требованию: прайс AWS
  даёт там 7,3616 / 1,60528 = 4,59 и 2,37822 / 1,60528 = 1,48, все пять тарифов задания совпадают с прайсом
  до округления (`task_crosscheck`). Для нас они не годятся дважды: регион us-east-1, а CPU идёт на споте —
  тогда порог для H100 по требованию 11,8 против c8i.8xlarge и 4,7 против c7a.16xlarge.
- Диск 50 ГБ gp3 (0,08 долл./ГБ·мес) и публичный IPv4 (0,005 долл./ч) добавляют 0,0105 долл./ч к каждой
  машине (прайс EC2 и VPC; ставки совпали с допущениями `runner.py`): S* с ними — 11,58 вместо 11,77.
- **Спот GPU не держать за константу.** Советник даёт p5.4xlarge экономию 9 %, а сторонний трекер (выдержка
  поиска, предварительно) — спот 1,34–3,76 долл./ч, то есть в 2–5 раз дешевле. Правило применяется к ценам на
  момент решения: фактическая цена спота — из `runner.py` (describe_spot_price_history) и `ledger.jsonl`.

### 6.2. Доля ядра (закон Амдала с поправкой на CPU GPU-машины)

S = 1 / ((1 − f)·r_host + f/s_k + x), где f — доля линейного решателя во времени CPU-эталона (раздел 4.2),
s_k — ускорение решателя на GPU, r_host — во сколько раз хост-часть медленнее на CPU GPU-машины
(`h100_host`/`l40s_host` против эталона), x — обмен CPU–GPU в долях T_CPU. Нужная доля f для S = 1,2·S*
(`amdahl_f_required`; «—» — недостижимо при f ≤ 1):

| пара | S требуемое | f при s_k = ∞, r_host = 1 / 2 / 4 | f при s_k = 20, r_host = 1 / 2 / 4 |
|---|---:|---|---|
| H100 по треб. / c8i по треб. | 5,51 | 0,82 / 0,91 / 0,95 | 0,86 / 0,93 / 0,97 |
| H100 по треб. / c8i спот | 14,12 | 0,93 / 0,96 / 0,98 | 0,98 / 0,99 / 0,99 (при s_k = 10 — недостижимо) |
| H100 по треб. / c7a спот | 5,59 | 0,82 / 0,91 / 0,96 | 0,86 / 0,93 / 0,97 |
| L40S по треб. / c8i по треб. | 1,79 | 0,44 / 0,72 / 0,86 | 0,47 / 0,74 / 0,87 |
| L40S по треб. / c8i спот | 4,60 | 0,78 / 0,89 / 0,95 | 0,82 / 0,91 / 0,96 |
| L40S спот / c8i спот | 2,12 | 0,53 / 0,76 / 0,88 | 0,56 / 0,78 / 0,89 |

Чтение: H100 окупается против CPU на споте, только если линейный решатель занимает больше 93 % времени
CPU и хост-часть на 16 vCPU не медленнее эталона. Поэтому f мерится **до** аренды GPU.

### 6.3. Влезает ли задача в видеопамять

Пик `gpu_peak_mem_GiB` ≤ 0,85 объёма (80 ГБ H100, 48 ГБ L40S; номинал NVIDIA в ГБ — 10⁹ байт, это строже
фактического объёма в ГиБ). Для расчётной сетки — экстраполяция по замеру двух сеток (байт на узел); оценка
4.3 — верхняя граница числа узлов, а не разрешение.

### 6.4. Порядок ворот (каждое «нет» — остановка без затрат на GPU)

1. **CPU сходится** на выбранной сетке (раздел 2) — сейчас нет.
2. **f на CPU** (оба способа 4.2) ≥ f из таблицы 6.2 при s_k = ∞ и r_host = 1 для самой дешёвой пары. Меньше —
   GPU не окупится ни при каком ускорении: не арендовать.
3. **Память:** оценка 4.3 по расчётной сетке ≤ 0,85 VRAM.
4. **Один прогон H100 по требованию** с заданной длительностью: `h100_host`, `h100_fp64`, `h100_mixed`;
   ворота точности (раздел 5).
5. **Окупаемость** по измеренным S, тарифу на дату и цене спота из журнала → `break_even.py compute
   --summary bench_summary.csv --compare bench_compare.csv` → вердикт по каждой GPU-конфигурации («переносить»
   — только если сошлась, точность в допусках, S ≥ 1,2·S*, память в пределах; не хватает замера — «не хватает
   измерений»).
6. **L40S** — только после того, как `h100_mixed` прошёл точность, и если по 6.1–6.3 L40S может окупиться.

GPU быстрее, но дороже (1 < S < 1,2·S*), — выбор ради срока серии, а не стоимости: решает владелец.

### 6.5. Своя CUDA-реализация сборки якобиана

Имеет смысл, только если после 6.4 хост-часть (1 − f)·r_host — главная доля времени GPU-прогона, а
экономия по плану серий N·(C_CPU − C_GPU) больше стоимости разработки. Это решение владельца; в протокол
входит только замер, который его обосновывает.

## 7. Когда переходить к нескольким узлам (Hpc8a + EFA)

1. **Память:** байт/узел сетки на CPU (замер `cpu_peak_rss_sum_GiB`, две сетки) × узлов расчётной сетки
   > 0,8 × память самой большой подходящей одной машины.
2. **Срок:** время одного расчёта > `series_time_limit_h` (задаёт владелец; сейчас EXAMPLE 12 ч), и сильное
   масштабирование внутри одной машины уже исчерпано (эффективность < 0,7).
3. **Независимые точки** (линия частоты, план для суррогата) — не причина для MPI на нескольких узлах: они
   идут параллельно отдельными задачами `runner.py` на одиночных машинах (квота 384 vCPU = 6 × 64).
4. **GPU и несколько узлов:** в SU2 v8.6.0 CUDA — только один процесс MPI, многоузловой GPU-расчёт не
   поддержан; несколько узлов — только CPU MPI.
5. **Регион.** По прайсу AWS от 14.11.2026 Hpc8a и Hpc7a в us-east-1 нет (из HPC там только Graviton hpc7g);
   в eu-north-1 есть: hpc8a.96xlarge — 8,49769 долл./ч (192 vCPU, 768 ГиБ, 300 Гбит/с), hpc7a.96xlarge —
   7,7252. Переход к Hpc8a — это запасной регион eu-north-1, а значит решение владельца; машины us-east-1 с
   EFA — не проверено (страница типов AWS недоступна).

**Проверка масштабирования** (по сводке `WRT_PERFORMANCE`: «Points/core», «Core-s/iter/Mpoints»):

- **сильное:** одна сетка, 1, 2, 4, 8, 16, 32, 64 процесса на c7a.16xlarge, затем 2 и 4 узла;
  E_s(N) = T₁ / (N·T_N); граница — 0,7;
- **слабое:** число узлов сетки на процесс постоянно (core130 coarse → medium: 139 121 → 584 357 узлов,
  в 4,2 раза; процессов — в 4 раза больше); E_w = T₁ / T_N; граница — 0,8;
- время — на итерацию (медиана) и до сходимости; меняется число итераций — записать: разбиение меняет
  предобуславливатель.

## 8. Как пользоваться

```bash
# 0) тарифы и сетки (публичные файлы скачиваются вручную, без ключей AWS)
curl -sSO https://pricing.us-east-1.amazonaws.com/offers/v1.0/aws/AmazonEC2/current/us-east-1/index.csv
curl -sSO https://spot-bid-advisor.s3.amazonaws.com/spot-advisor-data.json
python3 break_even.py tariffs --price-list us-east-1=index.csv --spot-advisor spot-advisor-data.json \
        --vpc-price us-east-1=vpc_us-east-1.json --out tariffs_us-east-1.json
python3 break_even.py mesh ../geom/out/taurus60_coarse/passage_taurus60_coarse.su2 --out mesh_stats.json
python3 break_even.py compute --bench bench.toml --out break_even.json
# 1) после прогонов: runs/<id>/ с history*.csv, *.log, gpu_mem*.csv, nsys_mem*.csv, profile.toml
python3 bench_parse.py --bench bench.toml --runs runs --out bench_summary.csv --compare bench_compare.csv
python3 break_even.py compute --bench bench.toml --summary bench_summary.csv --compare bench_compare.csv \
        --out break_even.json
python3 -m pytest -q test_bench.py
```

Смета протокола — формула «часы × тариф» (`break_even.json`), без живых цен спота её не публиковать.
Пример (EXAMPLE, не план): 2 ч на p5.4xlarge по требованию = 2 × (6,88 + 0,0105) = 13,78 долл.

## 9. Что не проверено и чего не хватает

- Сошедшегося расчёта канала на CPU нет — ворота 1 закрыты.
- `runner.py` запускает только спот и образы `ubuntu2404`/`al2023` без драйвера NVIDIA: для GPU-прогона
  нужны рынок «по требованию», образ с драйвером и CUDA и квоты семейств P и G (16 vCPU для p5.4xlarge,
  8 vCPU для g6e.2xlarge; одобренные 384 vCPU — квота стандартных семейств).
- Сборка SU2 с CUDA и сам CUDA-прогон в контейнере не проверялись (GPU нет).
- Команды `perf`, Tracy (`tracy-capture`, `tracy-csvexport`) и `nsys stats` на SU2 не запускались; формат CSV
  `nsys` в `bench_parse.py` — по документации, тест — на придуманном файле.
- develop.openfoam.com, www.nvidia.com, aws.amazon.com, docs.aws.amazon.com из контейнера недоступны: статус
  GPU в petsc4Foam, паспортные числа NVIDIA и принадлежность p5 к H100 SXM — по выдержкам поиска и выводу.
- Цена спота — оценка советника; для p5.4xlarge она расходится со сторонним трекером в разы.

## Источники (проверено 14.11.2026)

- SU2 v8.6.0: https://github.com/su2code/SU2/tree/v8.6.0 (клон, коммит 7048fd3) — `meson_options.txt`, `meson.build`,
  `Common/src/CConfig.cpp`, `Common/include/code_config.hpp`, `Common/src/linear_algebra/CSysMatrixGPU.cu`,
  `CSysVectorGPU.cu`, `CSysMatrix.cpp`, `CSysSolve.cpp`, `Common/include/linear_algebra/CPreconditioner.hpp`,
  `CMatrixVectorProduct.hpp`, `CSysSolve.hpp`, `SU2_CFD/src/output/COutput.cpp`, `CFlowOutput.cpp`,
  `CFlowCompOutput.cpp`, `SU2_CFD/src/drivers/CDriver.cpp`; документация — https://github.com/su2code/su2code.github.io
  (коммит 67d0f94): `_docs_v7/Guide-to-v8.md`, `Build-SU2-Linux-MacOS.md`, `Tracy-Integration.md`.
- OpenFOAM-dev: https://github.com/OpenFOAM/OpenFOAM-dev (клон, коммит fed424ac). PETSc `PCAMGX`:
  https://gitlab.com/petsc/petsc/-/raw/main/src/ksp/pc/impls/amgx/amgx.cxx . OGL:
  https://raw.githubusercontent.com/hpsim/OGL/main/README.md . AmgX: https://raw.githubusercontent.com/NVIDIA/AMGX/main/README.md .
  Выдержки поиска: https://www.esi-group.com/sites/default/files/resource/other/2611/OpenFOAM_Conference_2021_Abstract_Martineau_NVIDIA.pdf ,
  https://www.esi-group.com/sites/default/files/resource/other/1672/8th_OpenFOAM_Conference_NVIDIA_Posey.pdf ,
  https://develop.openfoam.com/modules/external-solver/-/issues/25 .
- CalculiX: https://github.com/Dhondtguido/CalculiX (клон, коммит 57a4e71) — `src/pastix.c`, `INSTALL.md`,
  `CMakeLists.txt`, `doc/CalculiX.tex`; пакет Ubuntu `calculix-ccx` 2.21-1 (`strings /usr/bin/ccx`).
- NVIDIA (выдержки поиска): https://www.nvidia.com/en-us/data-center/h100/ ,
  https://resources.nvidia.com/en-us-gpu-resources/h100-datasheet-24306 , https://developer.nvidia.com/blog/nvidia-hopper-architecture-in-depth/ ,
  https://www.nvidia.com/en-us/data-center/l40s/ ,
  https://images.nvidia.com/aem-dam/en-zz/Solutions/technologies/NVIDIA-ADA-GPU-PROVIZ-Architecture-Whitepaper_1.1.pdf .
- Тарифы: Price List AWS — https://pricing.us-east-1.amazonaws.com/offers/v1.0/aws/AmazonEC2/20261114012926/us-east-1/index.csv ,
  то же для eu-north-1, https://pricing.us-east-1.amazonaws.com/offers/v1.0/aws/AmazonVPC/20261105181945/us-east-1/index.json ;
  Spot Instance Advisor — https://spot-bid-advisor.s3.amazonaws.com/spot-advisor-data.json (sha256 — в
  `tariffs_us-east-1.json`). Сторонний трекер спота (выдержка поиска): https://compute.doit.com/spot/us-east-1/p5.4xlarge .
- Трек: навыки `su2-turbo` (разделы 6, 7), `aws-gtu-runs`, `calculix-rotor`, `openfoam-combustor`; README узлов
  `platform/su2`, `platform/aws`, `platform/geom`; сетки `platform/geom/out/` (sha256 — в `mesh_stats.json`).
