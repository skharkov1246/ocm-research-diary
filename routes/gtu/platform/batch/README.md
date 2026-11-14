# Первый тест в AWS: образ, AWS Batch, S3, смета — `routes/gtu/platform/batch/`

Узел платформы трека «Простая литая ГТУ-заменитель» (issue #59) по пункту 4 задания владельца «Расчёт
промышленной газовой турбины на AWS» (14.11.2026): Docker-образ расчёта, шаблоны AWS Batch для серий независимых
рабочих точек, структура входов и выходов в S3, смета «инстанс-часы × тариф» и порядок локальной подготовки.

Решения владельца 14.11.2026, на которых стоит узел: регион us-east-1 (квоты 384 vCPU обычные и спот одобрены,
бакет и роль машин есть); Стокгольм eu-north-1 — только запасной вариант параметром; сначала CPU, GPU — по
протоколу; потолок расходов задания — 300 долл., спот; расчётная точка первой — Taurus 60 в генераторном
исполнении (типоразмер KV-3, «№ 3»). В README трека с вечера 14.11.2026 записан и предохранитель раннера
2000 долл. — смета проверяется по обоим порогам, основной — 300.

**Состояние на 14.11.2026.** В AWS ничего не создавалось и не вызывалось (ключей в этом контейнере нет и не
нужно). Файлы и тесты готовы; образ собирался локально (раздел 1); шаблоны Batch сверены с моделью API AWS Batch
из botocore 1.44.3 (раздел 2); тарифы взяты из официального AWS Price List (раздел 4). Документация
docs.aws.amazon.com в сети сессии не разрешена — что сверено по другим источникам, помечено у каждого поля.

## Файлы

| Файл | Что это |
|---|---|
| `../docker/Dockerfile`, `../docker/Dockerfile.dockerignore`, `../docker/requirements.txt` | образ: Python 3.12, cantera, numpy, scipy, gmsh, boto3, CalculiX 2.21, SU2 v8.6.0 с MPI, код трека |
| `entrypoint.py` | точка входа образа: задание по JSON, индекс массива → строки плана, выгрузка после каждой точки |
| `batch.toml` | конфигурация серии: регион, типы машин, квота, тайм-аут, повторы, план, решатель (TOML) |
| `templates/*.json` | шаблоны запросов: `compute_environment`, `job_queue`, `job_definition`, `submit_array_job` |
| `render.py` | шаблоны → готовые запросы, `plan.csv`, `job.json`, `manifest.json`, `commands.sh`; сверка с моделью API; AWS не вызывает |
| `fetch_rates.py` → `aws_rates.toml` | тарифы us-east-1 и eu-north-1 из AWS Price List Bulk API + Spot Instance Advisor, с источниками |
| `cost_estimate.py`, `cost_estimate.toml` | смета; результат `cost_estimate_results.json` и `cost_estimate_lines.csv` |
| `bench_runner_jobs.json` | файл задач `runner.py` для первого замера (c7a.16xlarge и c8i.8xlarge) |
| `test_batch.py` | 27 тестов без сети и AWS |
| `LOCAL_FIRST.md` | что сделать локально до заказа машин; чем первый прогон в раннере отличается от Batch |
| `runs/` | вывод `render.py` (в `.gitignore`) |

```bash
cd routes/gtu/platform/batch
python3 -m pytest -q test_batch.py          # 27 passed
python3 render.py                           # us-east-1; --region eu-north-1 — запасной; --mode point — по точкам
python3 cost_estimate.py                    # выход 1, если итог плана выше потолка
python3 fetch_rates.py --cache /tmp/prices  # обновить тарифы (EC2 CSV — 300 МБ, читается потоком)
```

## 1. Образ

Собран и проверен локально 14.11.2026 (раздел «Что проверено»). Две стадии: `su2-build` собирает SU2 v8.6.0 с MPI (`meson -Dwith-mpi=enabled -Dbuildtype=release`, как
`../su2/install_su2.sh`), `runtime` — Ubuntu 24.04 с пакетами и кодом трека. Пользователь uid 1000 (OpenMPI не
запускается от root без особых переменных), `OMPI_MCA_btl_vader_single_copy_mechanism=none` (в контейнере без
CAP_SYS_PTRACE обмен MPI идёт через `/dev/shm`, её размер задаёт `sharedMemorySize` задания).

- **SU2 собирается, а не скачивается.** У выпуска v8.6.0 нет архивов под именами прошлых выпусков:
  `SU2-v8.6.0-linux64.zip`, `-linux64-mpi.zip`, `-linux64-omp.zip` — 404, тогда как у v8.5.0 `-linux64-mpi.zip` и
  `-linux64-omp.zip` отдаются (302), у v8.0.0–v8.3.0 — `-linux64.zip` и `-linux64-mpi.zip` (запросы к
  `github.com/su2code/SU2/releases/download/…`, 14.11.2026; список файлов выпуска через API GitHub закрыт — 403).
  Сборка из тега уже проверена узлом su2 на тех же Ubuntu 24.04 / OpenMPI 4.1.6 и связывает SU2 с той же OpenMPI,
  что в образе. Коммит тега `7048fd35f56d4bf52e586429d6746a36a10255f8` (`git ls-remote`) сверяется при сборке.
- **Основа — Ubuntu 24.04**, а не `python:3.12-slim` (Debian 12): системный Python 3.12, `calculix-ccx` 2.21
  (как в навыке calculix-rotor) и OpenMPI 4.1.6 — те же, на которых проверены узлы платформы.
- **Повторяемость:** версии Python закреплены; версии apt не закрепляются (архив Ubuntu убирает старые сборки),
  фактические пишутся в `/opt/gtu/PACKAGES.txt`; хеш коммита — `/opt/gtu/GIT_SHA`, `$GTU_GIT_SHA` и метка
  `org.opencontainers.image.revision`. Собирать из чистой копии коммита (`git worktree add --detach`).
- **hadolint 2.15.1** (колесо `hadolint-bin` с PyPI): без замечаний; два осознанных исключения помечены в файле
  (DL3008 — версии apt, DL3003 — `cd` в одной команде с секретом сборки).
- **`ccx -v` выходит с кодом 201** при нормальной работе (проверено) — проверка версии в образе сделана без пайпа.

Порядок и флаги сборки за прокси — `LOCAL_FIRST.md`, шаг 1. Код трека в образе принадлежит uid 1000: скрипты пишут рядом с собой (`run_case.py` строит сетку в `../geom/out/`).

## 2. AWS Batch: шаблоны и сверка полей

`render.py` подставляет значения из `batch.toml` в `templates/*.json` и проверяет каждый запрос:
1. **ParamValidator botocore 1.44.3** по модели API `batch` 2016-08-10 (обязательные поля, типы, лишние поля) —
   это машинное описание, из которого собран справочник API AWS Batch; плюс допустимые значения всех перечислений;
2. **ограничения из текста документации** (таблица ниже) и решений владельца (регион, квота);
3. тесты ломают запрос (неверное перечисление, нет `maxvCpus`, лишнее поле, `maxvCpus` 384, тайм-аут 30 с,
   11 попыток, массив из 1) — каждая поломка ловится.

Источники в таблице: **API** — текст документации поля в модели API `batch` botocore 1.44.3
(`botocore/data/batch/2016-08-10/service-2.json.gz`); **архив** — официальный архив руководства AWS Batch на GitHub
(`awsdocs/aws-batch-user-guide`, снимок 2023 г., может отставать); **поиск** — выдержка поиска 14.11.2026, страница
не открыта. Адреса справочника для сверки владельцем — в конце раздела.

### compute_environment.json (CreateComputeEnvironment)

| поле | значение | источник и цитата |
|---|---|---|
| `computeEnvironmentName` | `gtu-ce-<run_id>` | API: «up to 128 characters … letters, numbers, hyphens (-), and underscores (_)» |
| `type`, `state` | `MANAGED`, `ENABLED` | API: перечисление CEType; «A compute environment must be created in the ENABLED state» |
| `computeResources.type` | `SPOT` | API: CRType {EC2, SPOT, FARGATE, FARGATE_SPOT, ECS_MANAGED_INSTANCES} |
| `allocationStrategy` | `SPOT_PRICE_CAPACITY_OPTIMIZED` | API: «looks at both price and capacity to select the Spot Instance pools that are the least likely to be interrupted and have the lowest possible price … only available for Spot Instance compute resources» |
| `minvCpus` | 0 | API: «minimum number of vCPUs that a compute environment should maintain (even if … DISABLED)» — 0: простой бесплатен |
| `maxvCpus` | 320 | API (обязательное): «Batch never exceeds maxvCpus by more than a single instance» → 320 + 64 = 384 = квота |
| `instanceTypes` | `c7a.16xlarge`, `c7i.16xlarge`, `c8i.16xlarge` | API: «specify instance families … or specific sizes»; выбраны размеры: одно задание — одна машина, MPI на 64 vCPU |
| `subnets`, `securityGroupIds` | заглушки | API: subnets «required for … EC2, SPOT»; «One or more security groups must be specified» |
| `instanceRole` | заглушка | API: «The Amazon ECS instance profile … required for Amazon EC2 instances types … short name or full ARN» |
| `spotIamFleetRole` | заглушка | API, **разночтение**: у `type` — «If you choose SPOT, you must also specify … spotIamFleetRole», у самого поля — «required if the allocation strategy set to BEST_FIT or … isn't specified». Задано: удовлетворяет обоим; «must use the newer AmazonEC2SpotFleetTaggingRole managed policy» — иначе машины без тегов |
| `ec2Configuration` | `[{imageType: ECS_AL2023}]` | API: «If Ec2Configuration isn't specified, the default is ECS_AL2023»; поддержка образов AL2 для ECS кончилась 30.06.2026 — задано явно |
| `computeResources.tags` | `Project`, `RunId`, `Name` | API: «tags to be applied to Amazon EC2 resources that are launched» — теги прогона на машинах |
| `bidPercentage` | не задано | API: «If you leave this field empty, the default value is 100% of the On-Demand price. For most use cases, we recommend leaving this field empty» — потолок цены спота = цена по требованию (на нём верхняя граница сметы) |
| `serviceRole` | не задано | API: «If your account already created the Batch service-linked role, that role is used by default … attempts to create the Batch service-linked role» |
| `tags` (среды) | `Project`, `RunId` | API: «These tags don't propagate to the underlying compute resources» — поэтому теги и в `computeResources.tags` |

### job_queue.json (CreateJobQueue)

| поле | значение | источник и цитата |
|---|---|---|
| `jobQueueName`, `state` | `gtu-jq-<run_id>`, `ENABLED` | API: до 128 знаков; JQState |
| `priority` | 1 | API (обязательное): «higher integer value … evaluated first» — очередь одна |
| `computeEnvironmentOrder` | `[{order: 1, computeEnvironment: <ARN среды>}]` | API: «up to three compute environments»; «must be in the VALID state before you can associate them» |
| `jobStateTimeLimitActions` | RUNNABLE, 7200 с, CANCEL | API: «minimum value for maxTimeSeconds is 600 … maximum 86,400»; state «only supported value is RUNNABLE»; для ECS «only supported value is CANCEL» — задание, ждущее ёмкости или квоты 2 ч, отменяется, а не висит |

### job_definition.json (RegisterJobDefinition)

| поле | значение | источник и цитата |
|---|---|---|
| `type`, `platformCapabilities` | `container`, `["EC2"]` | API: JobDefinitionType; «If no value is specified, it defaults to EC2» |
| `containerProperties.image` | `<аккаунт>.dkr.ecr.<регион>.amazonaws.com/gtu-platform:<sha12>` | API: «Images in Amazon ECR repositories use the full registry and repository URI (for example, 123456789012.dkr.ecr.<region-name>.amazonaws.com/<repository-name>)» |
| `command`, `parameters` | `["--job", "Ref::job_uri"]`, `{job_uri: s3://…/inputs/job.json}` | API: «Default parameter substitution placeholders»; подстановка `Ref::` — архив `job_definition_parameters.md` |
| `jobRoleArn` | заглушка | API: «IAM role that the container can assume for Amazon Web Services permissions» — S3 из контейнера без ключей |
| `resourceRequirements` | VCPU 64, MEMORY 32 768 МиБ | API: ResourceType {GPU, VCPU, MEMORY}; «memory hard limit (in MiB) … at least 4 MiB»; поля `vcpus`/`memory` устарели. Память — меньше памяти машины: агент ECS регистрирует меньше установленной (архив `memory-management.md`) |
| `environment` | `GTU_RUN_ID`, `AWS_DEFAULT_REGION` | API: «Environment variables cannot start with "AWS_BATCH"» |
| `linuxParameters` | `initProcessEnabled: true`, `sharedMemorySize: 4096` | API: «init process … that forwards signals and reaps processes»; «size (in MiB) of the /dev/shm volume» |
| `logConfiguration` | `awslogs`, `awslogs-region`, `awslogs-stream-prefix: gtu` | API: LogDriver; «awslogs Specifies the Amazon CloudWatch Logs logging driver»; группа по умолчанию `aws/batch/job` — архив `using_awslogs.md` |
| `retryStrategy.attempts` | 2 | API: «between 1 and 10 attempts». 2 — один повтор, как `--retries 1` раннера; держит верхнюю границу сметы в 300 долл. |
| `evaluateOnExit` | `Host EC2*` → RETRY; `*` → EXIT | API: «up to 5 objects …; If none of the listed conditions match, then the job is retried» — поэтому последнее правило EXIT; `onStatusReason: "Host EC2*"` для прерывания спота — поиск: страница AWS «Use Amazon EC2 Spot best practices for AWS Batch» (`…/batch/latest/userguide/bestpractice6.html`) — не открыта; что именно совпадает с «Host EC2*», сверить по `statusReason` первого прерывания |
| `timeout.attemptDurationSeconds` | 12 600 | API: «minimum value for the timeout is 60 seconds. For array jobs, the timeout applies to the child jobs»; «If a job is terminated due to a timeout, it isn't retried»; архив `job_timeouts.md`: SIGTERM, через 30 с SIGKILL. Точка входа останавливает расчёт сама за 600 с до срока и выгружает точку |
| `propagateTags` | true | API: теги задания — на задачу ECS; «over 50, the job is moved to the FAILED state» (у нас 2) |

### submit_array_job.json (SubmitJob)

| поле | значение | источник и цитата |
|---|---|---|
| `arrayProperties.size` | 10 (режим line: 2 ступени × 5 частот) или 80 (point) | API: «The array size can be between 2 and 10,000» |
| индекс дочернего задания | `AWS_BATCH_JOB_ARRAY_INDEX` | архив `array_jobs.md`: «set to the container's corresponding job array index number. The first array job index is numbered 0»; ID дочернего — `<ID>:<индекс>` |
| попытка | `AWS_BATCH_JOB_ATTEMPT` | архив `job_retries.md`: «The first attempt is numbered 1» |
| `timeout` | 12 600 | API: «For array jobs, child jobs have the same timeout configuration as the parent job» |

**Индекс → строка плана.** `plan.csv`: строка = рабочая точка (index, row_id, line, ступень, доля частоты, точка,
доля и давление на выходе, `iter_scale`, `state`). Режим `point`: дочернее задание N считает строку index = N
(80 заданий, каждая точка с холодного старта). Режим `line` (основной): задание N считает строки line = N подряд —
линию частоты из 8 точек, каждая следующая с рестарта предыдущей (`--restart-from`, `--restart-state` узла su2);
по навыку su2-turbo точки одной линии идут подряд. На первой неудачной точке линия останавливается
(граница устойчивости стационарного расчёта), остальные — `skipped`.

**Запасной регион.** `python3 render.py --region eu-north-1` меняет ARN среды, адрес ECR и регион журналов;
типы те же (все три 16xlarge есть в Price List eu-north-1). Бакет остаётся в us-east-1: результаты и копия образа
идут между регионами (0,02 долл./ГБ, строки сметы `transfer_*`). Других регионов `render.py` не принимает.

**Роли — создаёт владелец** (агентские ключи этого не умеют и не должны). По тексту API: роль Spot Fleet с
управляемой политикой `AmazonEC2SpotFleetTaggingRole`; служебная роль Batch создаётся сама. Не проверено в сессии
(документация недоступна): профиль экземпляра ECS с политикой `AmazonEC2ContainerServiceforEC2Role` (регистрация
агента ECS, скачивание из ECR, журналы) — существующий `ALPHA_O_INSTANCE_PROFILE` раннера даёт только запись в S3
и, вероятно, не подойдёт; роль задания с доверием `ecs-tasks.amazonaws.com` и правами `s3:GetObject`,
`s3:ListBucket` на `<префикс>/gtu/*` и `s3:PutObject` на `<префикс>/gtu/<run_id>/*`.

**У API Batch нет DryRun** (в модели API такого параметра нет): проверка до расходов — только `render.py` и тесты.
Команды для владельца — `commands.sh` в папке вывода (создать, отправить, следить, разобрать); `render.py` их не
выполняет.

Справочник для сверки владельцем: API_ComputeResource, API_CreateComputeEnvironment, API_CreateJobQueue,
API_RegisterJobDefinition, API_SubmitJob, API_RetryStrategy, API_JobTimeout — `https://docs.aws.amazon.com/batch/latest/APIReference/`;
руководство: `array_jobs.html`, `job_timeouts.html`, `job_retries.html`, `bestpractice6.html`, `using_awslogs.html` —
`https://docs.aws.amazon.com/batch/latest/userguide/`.

## 3. Структура S3

```
s3://<бакет>/<префикс>/gtu/<run_id>/
├── manifest.json                    паспорт прогона: run_id, регион, git sha и «грязность» рабочей копии,
│                                    образ (URI; дайджест — после push), решатель (SU2 v8.6.0, коммит, сборка,
│                                    план стадий, сходимость, турбулентность), план (файл, sha256, строк,
│                                    режим, размер массива), Batch (имена, тайм-аут, повторы, vCPU, память, типы)
├── inputs/                          только чтение для заданий; скачивается в /work/inputs ($GTU_INPUTS)
│   ├── job.json                     задание: команда, режим, np, свой лимит времени
│   ├── plan.csv                     индекс → строка(и) плана
│   └── …                            сетки и конфигурации — по желанию
├── outputs/<row_id>/attempt-<n>/    всё из $GTU_OUT точки (run_summary.json, history, рестарт, поле SU2);
│   └── _status.json                 выгружается ПОСЛЕДНИМ: метка готовности точки
└── logs/<row_id>/attempt-<n>.log    вывод команды точки (тот же текст — в CloudWatch Logs)
```

- `<бакет>/<префикс>` — из `ALPHA_O_S3` (роль пишет только под префиксом). Ключ прогона `run_id` —
  `[A-Za-z0-9._-]`, до 48 знаков, как у раннера; один `run_id` не использовать и в раннере, и в Batch (раннер
  кладёт `<префикс>/gtu/<run_id>/<job>/`).
- Кто пишет: `manifest.json` и `inputs/` — `render.py` локально, в S3 их копирует владелец (`commands.sh`, шаг 1)
  до отправки; `outputs/` и `logs/` — точка входа после каждой точки, так что тайм-аут или прерывание теряют не
  больше одной точки. Повтор пишет в `attempt-2/`, первая попытка не затирается (у ключей агента нет `DeleteObject`).
- Версия кода — `manifest.json → code.git_sha`, она же в образе (`/opt/gtu/GIT_SHA`) и в каждом `_status.json`
  (`git_sha`); настройки решателя — `manifest.json → solver` и `run_summary.json` каждой точки (все параметры
  шаблона SU2, версия, сетка, по стадиям — статус, невязки, расход, π, η).
- `_status.json`: row_id, index, line, attempt, run_id, command, git_sha, batch_job_id, batch_ce, batch_jq,
  started_utc, finished_utc, exit_code, status (`ok`, `failed`, `timeout`, `terminated`, `no_time`), seconds.

## 4. Смета

`cost_estimate.py` → `cost_estimate_results.json` (у итога и у каждой строки — метод, ПО, дата, входы, формула,
адрес тарифа) и `cost_estimate_lines.csv`. Формула строк: вычисления = инстанс-часы × тариф; диск gp3 = ч × ГБ ×
ставка / 730; публичный IPv4 = ч × ставка; S3 = ГБ·мес × ставка + запросы; передача = ГБ × ставка; CloudWatch
Logs = ГБ × приём + ГБ·мес × хранение; ECR = ГБ·мес × ставка. Три уровня: **спот (оценка)** — ожидаемые часы ×
цена спота по Spot Instance Advisor; **по требованию** — ожидаемые часы × цена Price List; **верхняя граница** —
каждая попытка до своего срока, все попытки, цена по требованию (потолок спота у Batch при пустом
`bidPercentage`), объёмы данных × число попыток.

**Тарифы us-east-1** — AWS Price List Bulk API (публикация EC2 2026-11-14T01:29:26Z; S3 2026-11-09; CloudWatch
2026-11-04; VPC 2026-11-05; передача 2026-11-11; ECR 2026-11-13); страницы цен aws.amazon.com из контейнера
закрыты (403), поэтому источник — машинные файлы тех же цен. Адрес каждого файла — `aws_rates.toml`, `[sources]`.

| тариф | us-east-1 | eu-north-1 |
|---|---|---|
| c7a.16xlarge, по требованию / спот (оценка, экономия, прерывания) | 3,28448 / 1,47802 долл./ч (55 %, 10–15 %) | 3,51616 / 0,70323 (80 %, 5–10 %) |
| c7i.16xlarge | 2,856 / 0,94248 (67 %, > 20 %) | 3,0576 / 0,55037 (82 %, 5–10 %) |
| c8i.16xlarge | 2,99872 / 1,13951 (62 %, > 20 %) | 3,21056 / нет данных советника |
| c8i.8xlarge (32 vCPU, 64 ГиБ) | 1,49936 / 0,58475 (61 %, > 20 %) | 1,60528 / нет данных |
| диск gp3 | 0,08 долл./(ГБ·мес) | 0,0836 |
| публичный IPv4 в работе | 0,005 долл./ч | 0,005 |
| S3 Standard (первые 50 ТБ); PUT; GET | 0,023 долл./(ГБ·мес); 5·10⁻⁶; 4·10⁻⁷ долл./запрос | те же |
| выход в интернет (после 100 ГБ/мес бесплатных) | 0,09 долл./ГБ | 0,09 |
| между us-east-1 и eu-north-1 | 0,02 долл./ГБ | 0,02 |
| CloudWatch Logs: приём; хранение | 0,50 долл./ГБ; 0,03 долл./(ГБ·мес) | 0,54; 0,028 |
| ECR, хранение | 0,10 долл./(ГБ·мес) | 0,10 |

Сверка с документом владельца (там eu-north-1 по сайту costfluent): c8i.8xlarge 1,61 и c8i.2xlarge 0,40 долл./ч —
в Price List eu-north-1 1,60528 и 0,40132. Цена спота — **оценка**: файл советника (Last-Modified 20.03.2026)
даёт среднюю экономию к цене по требованию, а не цену на дату; живую цену читает `runner.py`
(`describe_spot_price_history`) в окружении Coldun перед каждым запуском.

**Итог (расчёт `cost_estimate.py`, 14.11.2026; часы — по модели времени итерации, см. «Входы сметы»)**

| сценарий | часов на задание | спот (оценка) | по требованию | верхняя граница |
|---|---|---|---|---|
| замер канала SU2, c7a.16xlarge (раннер, 64 процесса) | 1,26 (срок 2,5) | 1,87 | 4,14 | 17,58 |
| замер канала SU2, c8i.8xlarge (раннер, 32 процесса) | 2,32 (срок 3,5) | 1,38 | 3,50 | 11,08 |
| серия: 2 ступени × 5 частот × 8 точек, Batch, 10 линий | 2,43 (срок 3,5) | 36,09 | 79,95 | 246,96 |
| S3, передача, журналы, ECR (на весь план) | — | 2,55 | 2,55 | 4,96 |
| **итог плана** | | **41,90** | **90,15** | **280,58** |
| вариант: серия по точкам (80 заданий с холодного старта) | 1,24 | 150,04 | 329,36 | 1 980,51 |
| вариант: серия в eu-north-1 (бакет и образ в us-east-1) | 2,43 | 20,15 | 88,45 | 269,91 |

**Проверка потолка:** итог плана в пределах 300 долл. на всех трёх уровнях; запас верхней границы — 19,42 долл.
(и в пределах 2000 предохранителя раннера). Серия по линиям в 4,1 раза дешевле серии по точкам (по требованию:
79,95 против 329,36) — поэтому режим `line` основной; режим `point` при том же сроке 3,5 ч выходит за потолок по
верхней границе. Верхняя граница держится на `retry_attempts = 2` и сроке 3,5 ч: 5 попыток (как в примере AWS)
дали бы для серии 5 × 10 × 3,75 ч × 3,28 = 616 долл. — только с потолком 2000.

### Входы сметы, которых нет (EXAMPLE)

| вход | значение | единицы | владелец данных | влияние на смету |
|---|---|---|---|---|
| время итерации эталона `ref_s_per_iter` | 1,2255 (замер, грубый) | с | узел su2: замер на ненагруженной машине (LOCAL_FIRST, шаг 6) | время и стоимость ∝ |
| процессов эталона `ref_procs` | 4 | шт. | узел su2 (в журнале SU2 не печатается) | ∝ |
| доля идеального ускорения `parallel_eff` | 0,6 | 1 | замер `bench_*` | ∝ 1/η; при 0,4 серия не уложится в срок 3,5 ч — пересчитать срок |
| скорость vCPU по типам `cpu_factor` | 1,0 | 1 | замер `bench_*` (c7a — vCPU = ядро, c7i/c8i — гиперпоток; не проверено) | ∝ |
| узлов сетки канала `mesh_nodes` | 600 000 | шт. | узел geom (сетки taurus60 medium ещё нет; core130 medium — 584 357) | ∝ |
| итераций точки с рестарта `iters_next` | 1 000 | итераций | узел su2 | линия = 6370 + 7 × iters_next |
| подготовка машины `setup_min` | 15 / 20 (раннер, со сборкой SU2), 5 (Batch) | мин | первый прогон `smoke` | + часы |
| простой машины Batch до сжатия `tail_min` | 10 | мин | AWS (в документации не указан) | + часы |
| корневой диск машины Batch | 30 | ГБ | AWS (ECS AL2023 по умолчанию — не проверено) | доли цента |
| результатов на точку | 0,25 | ГБ | узел su2 | S3 и выход: 21 ГБ → 2,4 долл. |
| размер образа | 1,451 (замер `docker image inspect`, без сжатия) | ГиБ | сборка образа | ECR: 0,15 долл./мес |
| доли частоты (план) | 0,70; 0,80; 0,90; 1,00; 1,05 | 1 | gtu-cycle / su2 | на смету не влияют; положение линий |
| давление на выходе ступени (план) | 3,0·10⁵; 1,1·10⁶ | Па | gtu-cycle / su2 (одномерная оценка d1) | на смету не влияет; положение точек |
| геометрия второй ступени | то же колесо `taurus60` | — | узел geom (второго колеса нет) | ступень 2 пока не своей геометрией |

Каждое значение помечено EXAMPLE в `cost_estimate.toml` / `batch.toml`, строки плана — `example: true`.

## 5. Локально до заказа машин

`LOCAL_FIRST.md`: сборка образа, самопроверка, сетка, грубый прогон SU2 на 4 ядрах, проверка сходимости, замер
времени итерации, точка входа без AWS, шаблоны и смета, план раннера; таблица «раннер против Batch».

## Что проверено (14.11.2026, контейнер сессии)

- `test_batch.py` — 27 passed: рендер для us-east-1 и eu-north-1 в режимах line и point без ошибок сверки;
  поломки запросов ловятся; точка входа на локальном хранилище — выгрузка каждой точки, рестарт от предыдущей,
  остановка линии, режим point, ошибки задания, сухой прогон, тайм-аут молчащей команды, SIGTERM → выгрузка и
  код 143; смета совпадает с пересчётом, арифметика строк, потолок; файл задач замера принимает `runner.load_jobs`;
  разбор Price List на образце; статические проверки Dockerfile.
- `runner.py` не менялся: его тесты — 33 passed; `runner.py plan bench_runner_jobs.json` проходит.
- hadolint 2.15.1 — без замечаний.
- **Сборка образа** (`docker build`, BuildKit, основа `mirror.gcr.io/library/ubuntu:24.04` — Docker Hub ответил 429):
  код 0; стадия SU2 — 1205 с на 4 нагруженных ядрах (835 единиц компиляции), остальное — 134 с с кэшем стадии SU2.
  Образ 1 558 026 163 байт (1,451 ГиБ без сжатия), код трека в нём — 3,8 МБ, SU2 — 70 МБ. Пробные сборки поймали
  и исправили: нет `libGL.so.1` для gmsh (добавлен `libgl1`), uid 1000 в ubuntu:24.04 занят пользователем `ubuntu`.
  Собрано из рабочей копии с чужими незакоммиченными файлами (`GIT_SHA=<HEAD>-dirty`) — для AWS собирать из
  чистой копии коммита.
- **Образ в работе:** uid 1000; `/opt/su2/VERSION` = `v8.6.0 7048fd35…`; ccx `Version 2.21`; Python 3.12.3,
  cantera 3.2.0, numpy 2.5.3, scipy 1.18.1, gmsh 4.15.2, boto3 1.44.3; OpenMPI 4.1.6; `mpirun -np 2 SU2_CFD --help` —
  «SU2 v8.6.0 "Harrier"» без ошибок OpenMPI. Точка входа в контейнере (`--init`, без сети, хранилище — папка):
  режим line, индексы 0 и 1 → 2 + 2 точки ok, рестарт от предыдущей точки и состояние переданы, MPI 2 процесса,
  `git_sha` в `_status.json`; `docker stop` (SIGTERM через init) → статус `terminated`, вывод выгружен, код 143 за 1 с.
  Повтор в ту же папку с тем же `run_id` и попыткой дописывает поверх — один `run_id` на один прогон.

## Что не проверено

- Ничего не создавалось в AWS: роли, среда, очередь, определение, отправка, прерывание спота — не пробовались.
  У Batch нет DryRun; первая живая отправка — с владельцем, после замера `bench_*`.
- Документация AWS Batch не открыта (docs.aws.amazon.com в сети сессии не разрешён): поля сверены по модели API
  botocore, архиву руководства и одной выдержке поиска (`Host EC2*`). Сеть — решение владельца.
- Время итерации SU2 на машинах AWS, доля ускорения, скорость c8i против c7a — только модель; это и меряет
  первый прогон. Цены спота — оценка по советнику, не цена на дату.
- Публичный IPv4 у машин Batch (зависит от подсети), корневой диск ECS AL2023 30 ГБ, простой до сжатия среды —
  приняты.

## Источники

- Модель API AWS Batch: botocore 1.44.3, `botocore/data/batch/2016-08-10/service-2.json.gz` — https://pypi.org/project/botocore/1.44.3/
- Архив руководства AWS Batch: https://github.com/awsdocs/aws-batch-user-guide (`doc_source/array_jobs.md`,
  `job_timeouts.md`, `job_retries.md`, `job_env_vars.md`, `using_awslogs.md`, `memory-management.md`,
  `job_definition_parameters.md`, `service_limits.md`)
- Выдержка поиска 14.11.2026: «Use Amazon EC2 Spot best practices for AWS Batch» —
  https://docs.aws.amazon.com/batch/latest/userguide/bestpractice6.html ; «Introducing retry strategies for AWS Batch» —
  https://aws.amazon.com/blogs/compute/introducing-retry-strategies-for-aws-batch
- AWS Price List Bulk API: https://pricing.us-east-1.amazonaws.com/offers/v1.0/aws/index.json (файлы и даты — `aws_rates.toml`)
- Spot Instance Advisor: https://spot-bid-advisor.s3.amazonaws.com/spot-advisor-data.json
- SU2: тег v8.6.0 (`git ls-remote https://github.com/su2code/SU2`), выпуски — `github.com/su2code/SU2/releases/download/<тег>/…`
