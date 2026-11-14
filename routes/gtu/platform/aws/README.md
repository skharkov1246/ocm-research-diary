# Прогоны платформы ГТУ на AWS — `runner.py` (версия 2)

Узел платформы «Простая литая ГТУ-заменитель»: пакет независимых задач (SU2, OpenFOAM,
CalculiX, Gmsh, счёт на Python и любые скрипты репозитория) идёт на спотовых машинах в us-east-1.
Основа — `calc/aws_orchestrate.py` и `calc/aws_setup/`; отличия: спот, много задач
параллельно, бюджетный предохранитель, журнал расходов, сборщик только своих машин.

Состояние на 14.11.2026:
- **Первый живой пробный прогон** (сессия с ключами, 22:14, задача `smoke`, версия 1) прошёл весь путь:
  спот c7a.16xlarge → user-data → клон ветки → выгрузка 6 объектов под ролью машины в
  `<префикс ALPHA_O_S3>/gtu/<run>/smoke/` → `_DONE.json` → `collect` → гашение за 61 с → `ledger.jsonl`.
  Найдено: (1) образ `ubuntu2404` недоступен — `AccessDenied` на `ssm:GetParameter` параметра Canonical
  (дополнение политики не применено, рассчитывать на него нельзя); `al2023` работает, но системный Python там 3.9;
  (2) раннер не задавал подсеть — машина встала в us-east-1a по 0,955 долл./ч, хотя в 1c было 0,706, в 1b — 0,748
  (+25…35 %); (3) бюджетная проверка по максимуму цены по зонам — консервативна и верна, оставлена.
- **Версия 2** (этот код) отвечает на это: самая дешёвая зона с подсетью по умолчанию, образ `al2023` по
  умолчанию и Python нужной версии через uv, контейнеры Docker для SU2 / OpenFOAM / CalculiX, очередь под квоту,
  `NPROC` на машине. Проверено только тестами на заглушках (73 passed) и `plan`; живого запуска v2 не было —
  первый будет `jobs.smoke2.json` (раздел «Второй пробный прогон»).

## Файлы

| Файл | Что это |
|---|---|
| `runner.py` | Исполнитель: `plan`, `run` (по умолчанию DryRun; живой запуск — `--live`), `status`, `reap`, `collect` |
| `test_runner.py` | 73 теста без сети: botocore Stubber и свои заглушки AWS; user-data выполняется локально на подменённых командах (`dnf`, `apt-get`, `git`, `uv`, `docker`, `curl` IMDS) |
| `jobs.example.json` | Пример пакета: `smoke` (проверка машины) и `explore-249` (`routes/gtu/gtu_explore.py`, Python 3.12) |
| `jobs.smoke2.json` | Второй пробный прогон: `smoke-py` (al2023 + Python 3.12 + cantera) и `smoke-ccx` (контейнер ubuntu:24.04 + CalculiX) |
| `smoke2.py` | Проверки второго пробного прогона: версии, точка `gtu_cycle.cycle` и пул на `NPROC`, `ccx -v` и диск `disk.py`, зона машины |
| `iam_policy_addendum.json` | Дополнение к `calc/aws_setup/user_policy.json`: спот, цена спота, образ Ubuntu, необязательные `DescribeSubnets`/`DescribeImages` (не применено; v2 работает без него) |
| `ledger.jsonl` | Журнал расходов (коммитить после каждого боевого прогона) |
| `runs/<run_id>/` | Скачанные результаты, `summary.json`, журнал DRYRUN (в `.gitignore`) |

## Как идёт прогон

1. Задачи читаются из JSON и проверяются (имена, пакеты, срок, образ контейнера, отсутствие ключей в команде
   и в `container_setup`).
2. Перед каждым запуском: цена спота по зонам (`describe_spot_price_history`, окно 6 ч) → проверка бюджета
   **по максимуму цены за окно по всем зонам** → `run_instances`. Одновременно не больше `--max-parallel`
   машин (6) и не больше `--max-vcpu` vCPU (384 = 6 × 64; прежнее имя ключа `--vcpu-quota` тоже работает).
3. **Самая дешёвая зона.** Для типа машины зоны упорядочиваются по возрастанию последней цены спота
   (зоны без истории цены — тип там не продаётся — не пробуются). Машина просится в первую; при
   `InsufficientInstanceCapacity`, `SpotMaxPriceTooLow`, `Unsupported` — в следующую зону, затем — запасной тип
   `c7i.16xlarge` (со своей проверкой бюджета и своим порядком зон). Зона со сбоем ёмкости или прерыванием спота
   уходит в конец очереди на 30 мин (`avoid_s`). При квотных ошибках (`MaxSpotInstanceCountExceeded`,
   `VcpuLimitExceeded`, `InstanceLimitExceeded`) задача ждёт (п. 9). `MaxPrice` — по максимуму за окно × 1,2,
   поэтому верхняя граница бюджета от выбора зоны не зависит.
4. **Подсеть.** `describe_subnets` с фильтром `default-for-az=true` (один раз за прогон) → `SubnetId` подсети
   по умолчанию выбранной зоны. Подсеть по умолчанию сама выдаёт публичный IPv4 (`MapPublicIpOnLaunch`) — как и
   раньше без подсети; если у подсети это выключено, раннер задаёт сетевой интерфейс с
   `AssociatePublicIpAddress=true` (машине нужен интернет: git, pip, uv, docker, S3). Нет прав на
   `DescribeSubnets` (дополнение политики не применено) — зона задаётся `Placement.AvailabilityZone`, подсеть по
   умолчанию этой зоны выбирает EC2 (по модели API EC2 в botocore: без `SubnetId` берётся подсеть по умолчанию;
   сочетание с `Placement` вживую не проверено). Зона без подсети по умолчанию пропускается.
5. Параметры машины: разовый спот (`SpotInstanceType=one-time`, прерывание = terminate, `MaxPrice` задан),
   `InstanceInitiatedShutdownBehavior=terminate`, IMDSv2 обязателен (`HttpTokens=required`, hop limit 1), роль
   экземпляра из `ALPHA_O_INSTANCE_PROFILE`, без ключа SSH; группа безопасности — по умолчанию; gp3 с удалением
   при гашении; теги `Project=<ALPHA_O_TAG>`, `RunId=<run_id>`, `Job=<имя>`, `Name=gtu-<имя>` на машине и диске.
   `ClientToken` — свой в каждом процессе и для каждой зоны.
6. На машине (user-data): первой командой `shutdown -P +<max_min>` (сторож); база (git, python3-pip; для
   контейнера — docker); venv с boto3 (для выгрузки); пакеты задачи; `git clone --depth 1 -b <ветка>`;
   команда задачи из корня клона с `timeout` до остатка срока минус 8 мин на выгрузку. Три режима:
   - **system** — системный `python3` (как в версии 1);
   - **python** (`"python": "3.12"`) — uv (из PyPI в отдельный venv системного python3, при сбое — архив
     последнего выпуска с GitHub), `uv python install 3.12`, venv на этой версии (`uv venv --seed`, т. е. с pip),
     пакеты `pip` и файл `requirements` через `uv pip`; `PATH` указывает на venv;
   - **container** (`"container": "<образ>"`) — на al2023 ставится и запускается docker, образ скачивается,
     команды `container_setup` и затем команда задачи идут внутри контейнера (`--cpus` = все ядра,
     `--shm-size` = половина памяти, клон и `$GTU_OUT` примонтированы, пользователь root, OpenMPI под root
     разрешён `OMPI_ALLOW_RUN_AS_ROOT=1`). Выгрузка, `_DONE.json` и `poweroff` — снаружи; в ловушке выхода
     контейнер снимается (`docker kill`), поэтому машина гасится и при падении контейнера.
   Время каждой фазы пишется в журнал машины строками `[gtu-time] <фаза> <с>s total=<с>s`; зона машины по IMDS —
   строкой `[gtu] machine az=… planned=… cheapest=…`.
7. Выгрузка — под ролью машины в `s3://<бакет>/<префикс>/gtu/<run_id>/<job>/` (повтор N — в
   `…/gtu/<run_id>/_attemptN/<job>/`): всё из `$GTU_OUT`, `_job.log`, `_stdout.log` (в контейнере ещё
   `_container_setup.log`) и последним — метка `_DONE.json` (фаза, код выхода, время старта, команды и конца,
   `nproc`, зона машины по IMDS, запрошенная зона, тип). Затем `poweroff` (= terminate). Результаты выгружаются не
   дольше «до сторожа − 2 мин» (и не дольше 30 мин), метка — не дольше 60 с.
8. Исполнитель раз в `--poll-s` (60 с) читает метку своей попытки. Метка есть — скачивает результаты в
   `runs/<run_id>/<job>/` и гасит машину. Машина погашена без итога — прерывание спота (нет метки, код причины
   `Server.SpotInstanceTermination` или метка с кодом 143 от ловушки `TERM`): частичные файлы скачиваются в
   `runs/<run_id>/_lost/a<N>/<job>/`, повтор (`--retries 1`) идёт **в префикс новой попытки и другую зону** —
   итог задачи берётся только из одной попытки, результаты не смешиваются и не задваиваются. Метка с номером
   другой попытки не принимается. Истёк срок `max_min + grace_min` (10 мин) — гасит сам.
9. **Очередь под квоту.** Если машин пакета больше, чем помещается в `--max-vcpu` (или `--max-parallel`),
   остальные задачи ждут и запускаются по мере освобождения (по метке или гашению); бюджет проверяется перед
   каждым запуском. Задача, которой не хватило бюджета, пока работают машины прогона, не снимается: бюджет
   пересчитывается, когда закроется машина (верхняя граница «в работе» сменяется фактическим расходом). Квотная
   ошибка AWS при пустом списке своих машин (квоту держат чужие машины) — пауза `poll_s` и повтор, до
   `--quota-wait-polls` (15) раз, затем статус `quota`.
10. В `finally` (конец, ошибка, Ctrl-C, SIGTERM, SIGHUP) — внешний сборщик: ищет машины с тегами
   `Project=<ALPHA_O_TAG>` **и** `RunId=<run_id>`, сверяет теги ещё раз на своей стороне и гасит только их
   и известные ему id. Машины соседних прогонов с тем же `Project` не трогает (урок
   `calc/aws_orchestrate.py`: широкий сборщик однажды убил чужой прогон). `TerminateInstances`
   отвергает весь пакет из-за одного негодного id — тогда сборщик гасит по одной. Машина, гашение которой AWS не
   принял, остаётся в журнале открытой (в бюджете — верхней границей) до сверки.

## Формат задач

```json
{
  "defaults": {"instance_type": "c7a.16xlarge", "fallback_types": ["c7i.16xlarge"], "disk_gb": 50, "max_min": 30},
  "jobs": [
    {"name": "cycle-sweep", "python": "3.12", "pip": ["cantera", "numpy", "scipy"],
     "command": "python routes/gtu/gtu_explore.py \"$GTU_OUT/out.json\""},
    {"name": "ccx-disk", "container": "ubuntu:24.04",
     "container_setup": ["apt-get update -q", "apt-get install -y -q --no-install-recommends calculix-ccx python3-numpy python3-scipy"],
     "command": "python3 routes/gtu/platform/fea/disk.py --examples --out \"$GTU_OUT/disk.json\""}
  ]
}
```

| Поле | Что это |
|---|---|
| `name` | обязательно, `[A-Za-z0-9._-]`, не начинается с `_` |
| `command` | обязательно, bash (несколько строк допустимы), выполняется из корня клона |
| `instance_type`, `fallback_types` | тип и запасные типы (по умолчанию `c7a.16xlarge`, `c7i.16xlarge`) |
| `max_min` | 20…720, срок сторожа |
| `image` | `al2023` (по умолчанию) или `ubuntu2404` (по умолчанию для старых задач с полем `apt`) |
| `disk_gb` | 8…2000 |
| `apt` / `dnf` | системные пакеты (по образу) |
| `pip` | пакеты Python (в режиме python — через `uv pip`) |
| `python` | `"3.12"` или `"3.12.7"` — строкой: venv на этой версии через uv |
| `requirements` | путь файла зависимостей в репозитории, ставится после клонирования |
| `container` | образ Docker (`ubuntu:24.04`, `opencfd/openfoam-default:<тег>`, `public.ecr.aws/docker/library/ubuntu:24.04`, …); только на `al2023`; с ним поля `apt`, `pip`, `python`, `requirements` не применяются |
| `container_setup` | команды внутри контейнера до команды задачи (`bash -e`: первая ошибка — стоп) |
| `vcpu` | число vCPU для типов вне таблицы `VCPU` в коде |

Неизвестное поле — ошибка (ловит опечатки). Старые файлы задач работают без правок: `plan` проходит на
`jobs.example.json`, `jobs.smoke2.json`, `su2/aws_job.example.json`, `batch/bench_runner_jobs.json` (тест
`test_all_repo_job_files_plan`); но задачи с `"image": "ubuntu2404"` вживую упрутся в образ (см. «Типовые сбои»).

В команде доступны: `GTU_OUT`, `GTU_RUN_ID`, `GTU_JOB`, `GTU_ATTEMPT`, **`NPROC`** (число ядер — для пула
процессов), `OMP_NUM_THREADS=$NPROC`, `GTU_ITYPE`, `GTU_AZ` (зона, в которую просили), `GTU_AZ_CHEAPEST`
(самая дешёвая по цене при запуске), `GTU_AZ_MACHINE` (зона по IMDS). В контейнере `GTU_OUT=/gtu/out`, клон —
`/gtu/repo`; IMDS из контейнера недоступен (hop limit 1) — ключей роли у команды нет, зона — из `GTU_AZ_MACHINE`.
Docker Hub ограничивает анонимные скачивания по IP; у каждой машины свой адрес, но при отказе (429) образ можно
брать с зеркала `public.ecr.aws/docker/library/…`.

Окружение образа в контейнере (например, `bashrc` OpenFOAM) команда подключает сама: точка входа — `/bin/bash`,
`container_setup` и команда идут разными процессами.

## Бюджетный предохранитель

Перед каждым запуском:

```
итого = израсходовано + в работе + пакет (только DRYRUN/plan) + новая  ≤  потолок (--ceiling, 2000; отчёт владельцу после каждых 500 — report_step)
новая        = (MaxPrice + диск + IPv4) × (max_min + grace_min) / 60
MaxPrice     = максимум цены спота за 6 ч по всем зонам × 1,2 (--price-margin), вверх до 0,0001
израсходовано = Σ по закрытым машинам (цена спота в её фактической зоне при запуске + диск + IPv4) × max(60 с, конец − запуск)
               + --spent-before (расход вне журнала, ≥ 0)
в работе     = Σ по открытым машинам (MaxPrice + диск + IPv4) × весь отпущенный срок
диск         = disk_gb × 0,08 / 730 долл./ч (--ebs-usd-gb-month; ставка gp3 — допущение, в сессии не проверялась)
IPv4         = 0,005 долл./ч (--ipv4-usd-h): публичный адрес; ставка — допущение, в сессии не проверялась
```

Выбор самой дешёвой зоны снижает «израсходовано» (в журнале — цена фактической зоны: `spot_usd_h`, `az`,
`az_rank`, `az_cheapest`, порядок `zones`), но не прогноз: «новая» и «в работе» остаются по максимуму. Спотовая
машина не платит больше своего `MaxPrice`, сторож и сборщик гасят её не позже отпущенного срока, поэтому «в
работе» и «новая» — верхние границы (при верных ставках диска и IPv4; трафик и запросы S3 не учтены).
«Израсходовано» — оценка: цена спота за время работы может подняться до `MaxPrice`. Отказ пишется в журнал
(`ev=budget`, `ok=false`), машина не запускается. Пороги предохранителя проверяются при запуске: `--spent-before`,
ставки ≥ 0, `--grace-min` ≥ 5, `--poll-s` ≥ 5, `--price-margin` ≥ 1. В начале боевого прогона и после `reap`
открытые записи сверяются с AWS: погашенные закрываются по верхней границе срока.

Журнал `ledger.jsonl` живёт в контейнере сессии, а контейнер эфемерен: **после каждого боевого прогона коммитить
`ledger.jsonl`**, иначе следующая сессия начнёт счёт с нуля (тогда — `run --spent-before`). В журнале и выводе нет
ключей, имени бакета (в нём номер аккаунта; печатается как `<BUCKET>`) и id подсетей (журнал — в открытом
репозитории; пишется только способ задания зоны `net`: `subnet`, `eni`, `placement`).

## Где действуют ключи

Ключи ограниченного IAM-пользователя заданы переменными окружения в окружении **Coldun**:
`AWS_ACCESS_KEY_ID`, `AWS_SECRET_ACCESS_KEY`, `AWS_DEFAULT_REGION`, `ALPHA_O_S3=s3://<бакет>/<префикс>`,
`ALPHA_O_INSTANCE_PROFILE`, `ALPHA_O_TAG`. Запускать `run`, `reap`, `collect` — из сессии, открытой в этом
окружении. В других контейнерах ключи не действуют: только тесты и `plan`. Ключи не печатать и не вставлять в чат,
файлы, команды задач: исполнитель отказывается, если видит ключ в команде, `container_setup`, user-data или записи
журнала. На машинах ключей нет — S3 через роль экземпляра по IMDSv2.

## Как запускать из сессии

Из корня репозитория. Ветка должна быть запушена в публичный репозиторий — машина клонирует её
(`run` проверяет это через `git ls-remote`).

```bash
pip install boto3 pytest                          # если нет
python3 -m pytest -q routes/gtu/platform/aws/test_runner.py                        # 73 passed
python3 routes/gtu/platform/aws/runner.py plan routes/gtu/platform/aws/jobs.smoke2.json \
        --assume-price 0.955 --assume-zone-prices us-east-1a=0.955,us-east-1b=0.748,us-east-1c=0.706 \
        --show-userdata                           # без AWS; цены — допущение (здесь — числа прогона 14.11.2026)
python3 routes/gtu/platform/aws/runner.py status  # израсходовано / в работе / остаток
python3 routes/gtu/platform/aws/runner.py run routes/gtu/platform/aws/jobs.smoke2.json \
        --branch <ветка>                          # по умолчанию DryRun: права, бюджет, настоящий порядок зон
python3 routes/gtu/platform/aws/runner.py run routes/gtu/platform/aws/jobs.smoke2.json \
        --branch <ветка> --live                   # живой прогон (или GTU_LIVE=1); расходы
```

Без `--live` (или `GTU_LIVE=1`) `run` машин не создаёт: `run_instances` идёт с `DryRun=True` (с подсетью самой
дешёвой зоны). `--dryrun`, `GTU_DRYRUN=1`, `ALPHA_O_DRYRUN=1` включают проверку даже при `--live`.
`plan` печатает строку `[зоны]` (порядок зон и ожидаемая цена) и параметры `run_instances`; без
`--assume-zone-prices` зону в `plan` выбирает EC2 (в `run` — по настоящим ценам).

Долгий прогон — в фоне с журналом в файл (навык `aws-gtu-runs`, шаг 5: `setsid python3 -u …`). Если сессия
оборвётся, машины погасит сторож (`max_min`), а `reap` по `run_id` погасит остатки и сверкой закроет записи.

`plan --show-userdata` печатает user-data с маской `<BUCKET>` вместо имени бакета — это только вывод;
на машину уходит текст с настоящим (закавыченным) именем, его синтаксис проверяют тесты (`bash -n`).

## Второй пробный прогон (`jobs.smoke2.json`)

Цель — за минуты и центы проверить на живой машине всё новое в v2. Перед запуском в сессии с ключами:
закоммитить и запушить ветку с `runner.py`, `smoke2.py`, `jobs.smoke2.json` (машина клонирует ветку и запускает
`smoke2.py` из неё), прогнать тесты, затем `run` без `--live` и только потом `--live`.

| Задача | Что проверяет | Где смотреть |
|---|---|---|
| `smoke-py` | al2023 + uv + Python 3.12 + `pip` cantera, numpy, scipy; точка `gtu_cycle.cycle` (PR 16,8, T41 1150 °C) против локального расчёта 14.11.2026 (12,1664 МВт, КПД 0,3623; расхождение пишется, счёт не валит); пул процессов на `NPROC` точек; зона | `runs/<run>/smoke-py/smoke_py.json`, `az_check.json`, `nproc.txt`, `commit.txt`, `_job.log` (`[gtu-time] venv …`) |
| `smoke-ccx` | docker на al2023, контейнер `ubuntu:24.04`, `apt-get install calculix-ccx` внутри, `ccx -v` (выходит с кодом 201 — норма), диск `disk.py` (колесо схемы 130, сетка 30 × 4); `NPROC` внутри; зона | `runs/<run>/smoke-ccx/smoke_ccx.json`, `az_check.json`, `_container_setup.log`, `_job.log` (`[gtu-time] container_setup …`) |

Ожидается: обе задачи `status: ok`; в `summary.json` у обеих `cheapest_zone: true`, `az` = `az_cheapest` и
`az_machine` = `az` (зона из IMDS машины); в `az_check.json` — `"ok": true`; `nproc.txt` = 64; `status` — открытых
машин 0. Код 3 у задачи — счёт прошёл, но машина не в самой дешёвой зоне (причина — в `az_check.json` и строки
`launch_error` журнала: сбой ёмкости в дешёвой зоне — законный откат). Время подготовки — строки `[gtu-time]`.
Прогноз пакета при ценах 14.11.2026 — 1,35 долл. (верхняя граница по максимуму цены); ожидаемый расход — центы.

Локально (14.11.2026, без AWS) `smoke2.py py` и `smoke2.py ccx` прошли: Python 3.13, cantera 3.2.0 — точка
совпала с эталоном до 3·10⁻⁶; ccx 2.21 — диск за 0,56 с.

## Порядок проверки

1. Тесты: `python3 -m pytest -q routes/gtu/platform/aws/test_runner.py` — 73 passed.
2. `plan` — прочитать user-data, параметры `run_instances` и строку `[зоны]`.
3. `run` без `--live` в окружении Coldun: `запуск разрешён (DryRunOperation)` по каждой задаче, строка `[подсети]`
   (или отказ `describe_subnets` — тогда зона через `Placement`) и `чтение префикса результатов: OK`.
4. Живой `jobs.smoke2.json` (раздел выше), затем `status` и `reap --run-id <id> --list-only`.
5. Рабочий пакет; `explore-249` — сравнить скачанный `gtu_explore_results.json` с `routes/gtu/gtu_explore_results.json`.
6. Закоммитить `ledger.jsonl`.

Дополнение политики (`iam_policy_addendum.json`) применяет только владелец с админ-доступом; v2 от него не зависит:
```bash
sed -e s/REGION/us-east-1/g -e s/ACCOUNT_ID/<номер аккаунта>/g \
    routes/gtu/platform/aws/iam_policy_addendum.json | grep -v '"_comment"' > /tmp/gtu_add.json
aws iam put-user-policy --user-name alpha-o-agent --policy-name gtu-spot --policy-document file:///tmp/gtu_add.json
```

Сверено с моделью API EC2 в botocore 1.44.3: `MaxPrice` — потолок почасовой цены спота; base64 для `UserData` SDK
делает сам; `ClientToken` — до 64 ASCII; `InstanceInitiatedShutdownBehavior` по умолчанию `stop`; `HttpTokens=required`
отключает IMDSv1; `DescribeSpotPriceHistory` с окном возвращает и последнюю цену до начала окна; `TerminateInstances`
при одном негодном id не гасит ни одной машины пакета; фильтр `default-for-az` у `DescribeSubnets`; без `SubnetId`
EC2 берёт подсеть по умолчанию; `AssociatePublicIpAddress` — только у нового интерфейса eth0; `Placement.AvailabilityZone`.

**Не проверено вживую** (v2): uv и Python 3.12 на al2023; docker из `dnf` на al2023 и скачивание образа с Docker Hub;
`--shm-size` в половину памяти; сочетание `Placement` без `SubnetId`; права ключа на `DescribeSubnets` и
`DescribeImages`; путь `DescribeImages` для Ubuntu; код причины `Server.SpotInstanceTermination` и метка 143 при
прерывании; квотные коды при очереди; успевает ли выгрузка при прерывании.

## Типовые сбои

| Симптом | Причина и действие |
|---|---|
| статус `image_unavailable` | образ `ubuntu2404` не разрешён ключу (SSM Canonical, затем `DescribeImages`): `image: al2023`, или `container: ubuntu:24.04` с пакетами в `container_setup` |
| фаза `base` 10/11, `venv` 12/13, `packages` 16/17, `clone` 14/15 | сбой подготовки — `_job.log` |
| фаза `venv`, код 18 / 19 | uv не поставился (ни PyPI, ни выпуск GitHub) / `uv python install` не скачал Python |
| фаза `requirements`, код 23 | файла `requirements` нет в ветке или пакеты не ставятся |
| фаза `docker` 21, `pull` 22 | docker не запустился / образ не скачан (429 Docker Hub — зеркало `public.ecr.aws`) |
| фаза `container`, код 20 | упала команда `container_setup` — `_container_setup.log`; 125–127 — docker не запустил контейнер |
| `no_time` (124), `command` 124 / 137 | не хватило `max_min` |
| `[подсети] describe_subnets: UnauthorizedOperation` | нет прав — зона через `Placement`; это не ошибка |
| `[очередь] … квота занята (своих машин нет)` | квоту держат чужие машины; после `--quota-wait-polls` — статус `quota` |
| `cheapest_zone: false` | в дешёвой зоне не было ёмкости (`launch_error` с зоной) или зона после прерывания |

## Откат

- Погасить машины прогона: `python3 routes/gtu/platform/aws/runner.py reap --run-id <id>`
  (`--list-only` — только показать).
- Вручную, если исполнителя нет:
  ```bash
  aws ec2 describe-instances --region us-east-1 \
      --filters Name=tag:Project,Values=$ALPHA_O_TAG Name=tag:RunId,Values=<id> \
                Name=instance-state-name,Values=pending,running,stopping,stopped \
      --query 'Reservations[].Instances[].InstanceId' --output text
  aws ec2 terminate-instances --region us-east-1 --instance-ids <id1> <id2>
  ```
  Гашение разрешено политикой только для машин с тегом `Project=<ALPHA_O_TAG>`.
- Результаты в S3 удаляет владелец (у ключа агента нет `DeleteObject`): префикс `<префикс>/gtu/<run_id>/`
  (с подпапками `_attemptN/`).
- Скачать результаты попытки: `runner.py collect --run-id <id> --job <имя> [--attempt N]`.
- Права спота отозвать: `aws iam delete-user-policy --user-name alpha-o-agent --policy-name gtu-spot`.
- Ключ после работ — удалить (`calc/aws_setup/README.md`).
- Локально: `runs/` можно удалять; `ledger.jsonl` не удалять и не править вручную — он держит учёт.
