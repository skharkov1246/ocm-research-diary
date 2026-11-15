---
name: aws-gtu-runs
description: Прогоны трека «Простая литая ГТУ-заменитель» на AWS через доступ MatterForge скриптом routes/gtu/platform/aws/runner.py — где действуют ключи (окружение Coldun, дочерняя сессия), us-east-1, квоты 384 vCPU и очередь под них, спот c7a.16xlarge с запасным c7i.16xlarge в самой дешёвой зоне (подсеть по умолчанию зоны), образ al2023 с Python нужной версии через uv, контейнеры Docker (SU2, OpenFOAM, CalculiX), NPROC на машине, бакет из ALPHA_O_S3 с подпапкой gtu/, IMDSv2, автогашение тремя способами, бюджетный предохранитель 2000 долл. (решение владельца 14.11.2026; отчёт после каждых 500) и журнал ledger.jsonl, пакет задач, сбор и проверка результатов, прерывание спота, запрет ключей в файлах и выводе; версия 3 — соседние мощности: регионы us-east-2 и us-west-2 флагом --regions (фактическая квота каждого региона, образ al2023 через DescribeImages вне us-east-1, один бакет, MaxPrice и бюджет по максимуму по регионам) и запасные машины по требованию --on-demand-fallback. Используй, когда ставишь расчёт SU2, OpenFOAM, CalculiX или скрипт трека на облачные машины, проверяешь права и бюджет, гасишь машины прогона, забираешь результаты или разбираешь сбой прогона. Триггеры — «прогон на AWS», «запусти на споте», «runner.py», «облачный расчёт ГТУ», «Coldun», «ALPHA_O_S3», «погаси машины», «бюджет AWS», «прерывание спота», «ledger.jsonl», «дешёвая зона», «контейнер на AWS», «Python 3.12 на машине», «очередь под квоту», «smoke2», «соседние мощности», «соседний регион», «--regions», «us-west-2», «us-east-2», «машины по требованию», «on-demand».
---

# Прогоны трека ГТУ на AWS

Дорогие расчёты платформы (SU2 — ступени, OpenFOAM — камера, CalculiX — диски, `gtu_explore.py` и любой
скрипт репозитория) идут пакетом независимых задач на спотовых машинах us-east-1. Исполнитель —
`routes/gtu/platform/aws/runner.py`; основа — `calc/aws_orchestrate.py` и `calc/aws_setup/`. Состояние на
14.11.2026: первый живой пробный прогон версии 1 (`smoke`, 22:14) прошёл весь путь до гашения и журнала; по его
находкам сделана версия 2 — самая дешёвая зона, al2023 + uv, контейнеры, очередь, `NPROC` (README узла, раздел
«Состояние»). Версия 2 проверена тестами на заглушках — 94 passed (прогнано 14.11.2026, после правок трёх
проверок: неподтверждённое гашение, сеть в сборщике, машина за ошибкой 5xx, сверка, ёмкость, метка 143,
окружение OpenFOAM) — и `plan`; живого запуска v2 ещё не было, первый — `jobs.smoke2.json`.
15.11.2026 — версия 3 «соседние мощности» (решение владельца: «можешь занять ещё соседние мощности»; потолок
2000 долл. и отчёт каждые 500 — прежние): флаги `--regions` и `--on-demand-fallback`, раздел 12. Без новых флагов
раннер работает точно как версия 2. Тесты — 112 passed (прогнано 15.11.2026); живого запуска v3 не было.

## 1. Что где лежит

| что | где |
|---|---|
| исполнитель: `plan`, `run` (DryRun по умолчанию, `--live`), `status`, `reap`, `collect` | `routes/gtu/platform/aws/runner.py` |
| тесты без сети (botocore Stubber и свои заглушки), пример пакета `smoke` + `explore-249` | `test_runner.py`, `jobs.example.json` там же |
| второй пробный прогон: `smoke-py` (al2023 + Python 3.12) и `smoke-ccx` (контейнер ubuntu:24.04 + CalculiX), проверки | `jobs.smoke2.json`, `smoke2.py` там же |
| порядок проверки, формат задач, формула бюджета, откат | `routes/gtu/platform/aws/README.md` |
| права спота (дополнение, применяет владелец) | `routes/gtu/platform/aws/iam_policy_addendum.json` |
| справочник версии 3: цены по требованию по регионам и межрегиональная передача (дата, источник, «предварительно») | `routes/gtu/platform/aws/aws_prices_ref.json` |
| журнал расходов (появится после первого боевого прогона; коммитить) | `routes/gtu/platform/aws/ledger.jsonl` |
| результаты, `summary.json`, журнал проверок `journal.dryrun.jsonl` (в `.gitignore`) | `routes/gtu/platform/aws/runs/<run_id>/` |
| разовая настройка аккаунта (владелец, админ-ключи, не агент) | `calc/aws_setup/setup.sh`, `user_policy.json`, `instance_role_policy.json` |
| предшественник: одна машина on-demand для alpha-O, не для ГТУ | `calc/aws_orchestrate.py` |
| флот eu-central-1 других треков — не трогать | `calc/fleet_status.py` |

## 2. Где действуют ключи

- **Только в окружении Coldun.** Сверено 14.11.2026 (`get_environment`): личное окружение, сеть без ограничений,
  переменные `AWS_ACCESS_KEY_ID`, `AWS_SECRET_ACCESS_KEY`, `AWS_DEFAULT_REGION`, `ALPHA_O_S3`,
  `ALPHA_O_INSTANCE_PROFILE`, `ALPHA_O_TAG` (видны только имена), скрипта настройки и репозиториев нет.
- В любом другом контейнере ключи не действуют: только тесты и `plan`, вызовов AWS не делать.
- Идущую сессию в Coldun не перевести (`switch_environment` не переводит на неограниченную сеть) — нужна
  **дочерняя сессия**: `list_environments` → id окружения с именем Coldun (id в файлы не писать) →
  `create_session` с `environment_id`, `source_url` репозитория, `source_revision` = запушенная ветка, задачей
  в `prompt`; `permission_mode` не `plan` (иначе дочерняя ждёт одобрения). Ход — `get_session`
  (`status_bucket`) и `list_events`; чисто закончившаяся сессия может не отчитаться сама.
- Согласие на живой запуск. Сообщение `send_message` приходит в дочернюю сессию как от сессии, а не от
  человека, и само по себе новых затратных действий не разрешает (описание инструмента); как помечается
  начальный `prompt` `create_session` — не проверено. Поэтому задача дочерней — шаги 0–4, а `--live`
  владелец разрешает своим сообщением в самой дочерней сессии (claude.ai/code).
- По документации облачных окружений: значения переменных читает любой, кто пользуется окружением;
  идущая сессия держит прочитанные значения до перезапуска своей ВМ. ВМ сессии засыпает через несколько
  минут без активности; не даёт ли ей уснуть фоновый процесс — не проверено. Уснул исполнитель — машины
  гасит сторож (`max_min`), результаты остаются в S3: `collect`, затем `reap --run-id` (сверка закроет журнал).
- Доступ проверять вызовом API, не печатая значений ключей (ни длины, ни начала): строка «вызов от имени» в
  `run` (номер аккаунта скрыт) или `aws sts get-caller-identity --query Arn --output text | sed -E
  's/[0-9]{12}/<аккаунт>/'` (есть ли `aws` в Coldun — не проверено). В имени бакета номер аккаунта — в выводе
  `<BUCKET>`.

## 3. Параметры, которые нельзя менять молча

| параметр | значение | откуда |
|---|---|---|
| регион | us-east-1; соседние — только флагом `--regions` (версия 3, раздел 12) | `REGION` в коде; eu-central-1 занят флотом `fleet_status.py` |
| регион ручных команд | `--region` региона машины (us-east-1, если журнал не пишет `region`) | значение `AWS_DEFAULT_REGION` в Coldun не читалось |
| квоты | us-east-1: 384 vCPU спот (заявка до 1536 на рассмотрении) и 384 по требованию — проверка 15.11.2026 | `--max-vcpu 384` (= `--vcpu-quota`) = 6 машин × 64; `--max-parallel 6`; сверх — очередь; версия 3 — квота каждого региона (раздел 12) |
| типы | `c7a.16xlarge`, запас `c7i.16xlarge` (по 64 vCPU) | `DEFAULT_TYPE`, `DEFAULT_FALLBACK`, `VCPU` |
| зона | самая дешёвая по последней цене спота; сбой ёмкости — следующая зона, затем запасной тип | `zone_order`, подсеть по умолчанию (`describe_subnets`), без прав — `Placement` |
| образ | `al2023` по умолчанию; `ubuntu2404` — только старые задачи с `apt` (SSM Canonical агенту запрещён) | `DEFAULT_IMAGE`, `IMAGES` |
| спот | `MarketType=spot`, `SpotInstanceType=one-time`, `InstanceInterruptionBehavior=terminate`, `MaxPrice` | `build_params` |
| IMDSv2 | `HttpTokens=required`, `HttpPutResponseHopLimit=1`, `HttpEndpoint=enabled` | на машине нет ключей, S3 — роль |
| запись | `s3://<бакет>/<префикс>/gtu/<run_id>/<job>/` | роль машины: только `s3:PutObject` в `<префикс>/*` |
| чтение | `ListBucket` по префиксу, `GetObject`; `DeleteObject` нет | `user_policy.json` |
| теги | `Project=<ALPHA_O_TAG>`, `RunId`, `Job`, `Name=gtu-<job>` на машине и диске | политика гасит только машины с тегом `Project` |
| потолок | 2000 долл. (`--ceiling`), отчёт после каждых 500 (`report_step`) | журнал `ledger.jsonl` |

Квоты EC2 считаются в vCPU (архив руководства EC2, `using-spot-limits.md`). Коды: `L-1216C47A` — Running
On-Demand Standard (A, C, D, H, I, M, R, T, Z) instances (сверено по примеру awscli 1.46.3
`examples/service-quotas/get-service-quota.rst`); `L-34B43A08` — All Standard (A, C, D, H, I, M, R, T, Z) Spot
Instance Requests (по выдержкам поиска, не сверено). Прав Service Quotas у ключа агента нет — смотрит владелец:
`aws service-quotas get-service-quota --region us-east-1 --service-code ec2 --quota-code L-34B43A08`.

## 4. Порядок работы

**Шаг 0.** Начало сессии по `CLAUDE.md`; прочитать README узла и навык расчёта (`su2-turbo`,
`openfoam-combustor`, `calculix-rotor`, `gtu-cycle`, `surrogate-opt`).

**Шаг 1. Пакет задач.** JSON `{"defaults": {...}, "jobs": [...]}`; поля задачи: `name` (`[A-Za-z0-9._-]`),
`command` (bash, из корня клона), `instance_type`, `fallback_types`, `max_min` (20…720 — срок сторожа),
`image` (`al2023` по умолчанию; `ubuntu2404` — старые задачи с `apt`), `disk_gb`, `apt`/`dnf`, `pip`, `vcpu`,
`python` (`"3.12"` строкой — venv через uv), `requirements` (файл в репозитории), `container` (образ Docker) и
`container_setup` (команды внутри контейнера). Неизвестное поле — отказ. Какой режим брать:
- счёт на Python (Cantera, numpy, scipy, `gtu_explore.py`, суррогат) — `"python": "3.12"` + `pip`;
- решатели из apt (CalculiX, OpenMPI, Gmsh) — `"container": "ubuntu:24.04"` + `container_setup`
  (`apt-get update -q`, `apt-get install -y -q --no-install-recommends …`); OpenFOAM — образ `opencfd/openfoam-*`
  или Ubuntu; SU2 — сборка `install_su2.sh` в контейнере или образ из `routes/gtu/platform/docker/Dockerfile`, когда
  он опубликован. Точка входа образа заменена на `/bin/bash` (в образе нужен `bash`; `alpine` не годится);
  окружение образов OpenCFD (`/openfoam/bash.rc`) раннер подключает сам до `container_setup`, окружение других
  образов (`/opt/openfoam<N>/etc/bashrc` у openfoam.org) команда подключает сама. `apt-get` в команде без
  `container` на al2023 не сработает — `plan` и `run` печатают `[внимание]`;
- пакет из десятков частей — по задаче на часть, имя `part01…`; раннер сам держит очередь под квоту.
Команда пишет всё в `$GTU_OUT`; доступны `GTU_RUN_ID`, `GTU_JOB`, `GTU_ATTEMPT`, `NPROC` (ядра — для пула
процессов), `OMP_NUM_THREADS=$NPROC`, `GTU_AZ`, `GTU_AZ_CHEAPEST`, `GTU_AZ_MACHINE`, `GTU_ITYPE`, `python` из venv
машины. В контейнере `GTU_OUT=/gtu/out`, клон — `/gtu/repo`, пользователь root, IMDS недоступен.

**Шаг 2. Офлайн-проверка (любой контейнер).**
```bash
python3 -m pytest -q -p no:cacheprovider routes/gtu/platform/aws/test_runner.py   # 112 passed
python3 routes/gtu/platform/aws/runner.py plan <jobs.json> --assume-price 1.5 --show-userdata \
        [--assume-zone-prices us-east-1a=0.955,us-east-1b=0.748,us-east-1c=0.706] \
        --out <scratchpad>/runs                                                    # цены — допущение
```
Прочитать user-data, параметры `run_instances` и строку `[зоны]` (порядок зон, ожидаемая цена). `plan` в `ledger.jsonl` не пишет (проверки бюджета — в
`<out>/<run_id>/journal.dryrun.jsonl`); прогноз при допущенной цене — не число для отчёта. Код выхода `plan`:
0 — все задачи спланированы; 1 — хоть одна нет (бюджет, нет цены; строка `НЕ спланированы`); 2 — отказ проверки
(`--assume-price` ≤ 0, user-data любой попытки, включая повтор, больше 16 КБ, ошибка в задаче).

**Шаг 3. Ветка запушена.** Машина клонирует публичный репозиторий (`git clone --depth 1 -b <ветка>`);
`run` сверяет ветку через `git ls-remote`.

**Шаг 4. Проверка прав (Coldun).** `python3 -c 'import boto3'` (нет — `pip install boto3`), затем
`python3 routes/gtu/platform/aws/runner.py run <jobs.json> --branch <ветка>` без `--live`. Ожидается
`запуск разрешён (DryRunOperation)` по каждой задаче и `чтение префикса результатов: OK`; код выхода 1, если
хоть одна задача не `dryrun_ok`. Цена спота по зонам, подсети по умолчанию (`[подсети]`; отказ `describe_subnets` —
не ошибка, зона пойдёт через `Placement`) и образ из SSM при этом читаются настоящими вызовами.

**Шаг 5. Живой прогон** — только с явного разрешения владельца (п. 2). Для версии 2 первый — `jobs.smoke2.json`
(две задачи по 25 мин сторожа, прогноз 1,35 долл. по максимуму, расход — центы; что смотреть — README узла,
раздел «Второй пробный прогон»), затем рабочий пакет:
```bash
RUN=gtu-$(date -u +%Y%m%dT%H%M%SZ)-smoke          # run_id: [A-Za-z0-9._-], до 48 знаков
setsid python3 -u routes/gtu/platform/aws/runner.py run <jobs.json> --branch <ветка> \
       --run-id "$RUN" --live > "/tmp/gtu_$RUN.log" 2>&1 < /dev/null &
echo "$RUN"
```
- `--run-id` задать заранее и записать: иначе он виден только в строке `[прогон] <run_id> режим=run …`
  (вторая, после «вызов от имени») и в имени `runs/<run_id>/`, а эта папка появляется лишь при первом сборе.
- `-u`: без него вывод в файл буферизуется и строки журнала приходят с задержкой.
- `setsid`: фоновая команда инструмента Bash (`run_in_background`) останавливается по своему `timeout`
  (по умолчанию 30 мин, не больше 2 ч), и прогон обрывается (по SIGTERM исполнитель гасит машины в `finally`).
  Отвязанный `setsid` процесс переживает конец вызова Bash (проверено в контейнере 14.11.2026).
- Ход — `tail "/tmp/gtu_$RUN.log"` и `runner.py status`.

**Шаг 6. Сбор.** По метке `_DONE.json` исполнитель сам скачивает `<префикс>/gtu/<run_id>/<job>/` в
`runs/<run_id>/<job>/` и гасит машину. Вручную: `runner.py collect --run-id <id> --job <имя>`.

**Шаг 7. Учёт.** `runner.py status` (открытых машин 0), затем `ledger.jsonl` и проверенные results-JSON — в
ветку и PR по правилам `CLAUDE.md`; контейнер эфемерен, незакоммиченный журнал пропадает. Пуш из дочерней
сессии не проверялся — проверить `git push` до живого прогона.

## 5. Автогашение — три независимых способа

1. **Выключение = гашение.** `InstanceInitiatedShutdownBehavior=terminate` (по умолчанию у EC2 — `stop`) и
   разовый спот с `terminate` при прерывании. `halt` машину не выключает: процессор встаёт в `HLT`, машина
   продолжает работать — в командах только `poweroff` или `shutdown -P`.
2. **Сторож на машине.** Первая команда user-data (после `#!/bin/bash`) — `shutdown -P +<max_min>` (`+m` —
   через m минут). Ловушка `EXIT` выгружает `$GTU_OUT` (с `_job.log`, `_stdout.log`) не дольше «до сторожа −
   2 мин» и не дольше 30 мин, затем метку `_DONE.json` (до 60 с), затем `poweroff`. Команда задачи идёт до
   «старт скрипта + `max_min` − 8 мин»; подготовка машины входит в этот срок.
3. **Внешний сборщик.** В `finally` (конец, ошибка, Ctrl-C, SIGTERM, SIGHUP) — гашение машин с тегами
   `Project` **и** `RunId` этого прогона; срок задачи у сборщика — `max_min + grace_min` (10). Сетевой сбой
   у `TerminateInstances` не обрывает сборщик: повтор по одной и строка `[гашение][ОШИБКА] … погасите вручную:
   aws ec2 terminate-instances … <id>`. Гашение не подтверждено посреди прогона — ещё одна попытка, не вышло —
   новые машины не запускаются (остаток очереди — `halted`), чтобы живых vCPU не стало больше `--max-vcpu`. Если
   сессия умерла, 1 и 2 всё равно гасят; не сработали (user-data не стартовал) — сверка в начале следующего
   прогона или после `reap` гасит машину, живую после срока, и считает её по факту. Остатки —
   `runner.py reap --run-id <id>` (`--list-only` — только показать) или вручную командами из раздела «Откат» README:
   ```bash
   aws ec2 describe-instances --region us-east-1 \
       --filters Name=tag:Project,Values=$ALPHA_O_TAG Name=tag:RunId,Values=<id> \
                 Name=instance-state-name,Values=pending,running,stopping,stopped \
       --query 'Reservations[].Instances[].InstanceId' --output text
   aws ec2 terminate-instances --region us-east-1 --instance-ids <id1> <id2>
   ```

## 6. Бюджетный предохранитель

Перед каждым запуском (код `budget_check`, `Ledger.summary`):
```
итого    = израсходовано + в работе + пакет (только DryRun/plan) + новая  ≤  2000 (--ceiling)
новая    = (MaxPrice + диск + IPv4) × (max_min + grace_min) / 60,   grace_min = 10
MaxPrice = максимум цены спота за 6 ч по всем зонам × 1,2 (--price-margin), вверх до 0,0001
диск     = disk_gb × 0,08 / 730 долл./ч;   IPv4 = 0,005 долл./ч   (обе ставки — допущения кода, не сверены)
```
Открытые машины — по верхней границе `(MaxPrice + диск + IPv4) × весь срок` (в живом прогоне после сверки —
× max(срок, сколько машина уже живёт): машина, которую не погасили ни сторож, ни сборщик, ни сверка, идёт по
факту; `status` — по сроку); закрытые — по цене спота в своей
фактической зоне на момент запуска и фактическому времени (не меньше 60 с) — это оценка, не граница. Самая
дешёвая зона снижает фактический расход (прогон 14.11.2026: 0,706 в 1c против 0,955 в 1a), но проверка бюджета
и `MaxPrice` остаются по максимуму цены по зонам — не менять. В очереди задача, которой не хватило бюджета при
своих машинах в работе, ждёт: бюджет пересчитывается, когда закрывается машина. Трафик и запросы S3
не учтены. Сверка формулы: `plan` `jobs.example.json` при допущенной цене 1,5 даёт 0,91 (`smoke`) + 1,51
(`explore-249`) = 2,41 долл. Отказ пишется в журнал (`ev=budget`, `ok=false`), машина не создаётся. Расход вне
журнала — `run --spent-before`.

Внешняя страховка — AWS Budgets, ставит владелец (у ключа агента прав нет). Budgets только оповещает и машин
не гасит; без фильтров считает весь аккаунт, включая флот других треков:
```bash
aws budgets create-budget --account-id <аккаунт> \
  --budget '{"BudgetName":"gtu","BudgetLimit":{"Amount":"2000","Unit":"USD"},"TimeUnit":"MONTHLY","BudgetType":"COST"}' \
  --notifications-with-subscribers '[{"Notification":{"NotificationType":"ACTUAL","ComparisonOperator":"GREATER_THAN","Threshold":80,"ThresholdType":"PERCENTAGE"},"Subscribers":[{"SubscriptionType":"EMAIL","Address":"<почта владельца>"}]}]'
```

## 7. Прерывание спота

- AWS предупреждает примерно за 2 минуты: на машине появляется `/latest/meta-data/spot/instance-action`
  (до этого — 404), в EventBridge — событие «EC2 Spot Instance Interruption Warning»; раньше может прийти
  рекомендация `events/recommendations/rebalance`.
- Исполнитель: машина `shutting-down`/`terminated` без итога — код причины `Server.Spot…`
  (`SpotInstanceTermination`, `SpotInstanceShutdown`), либо причины нет, а `_DONE.json` нет или в нём код 143
  (ловушка `TERM` успела выгрузить). Причина `Client.…` (сторож `shutdown -P`, `poweroff`, наш terminate) — не
  прерывание: метка 143 тогда — итог `failed` без повтора (повтор зависшей задачи удвоил бы расход). Метка 143 при
  ещё живой машине — раннер ждёт причину до трёх опросов. Прерывание → частичные файлы
  скачиваются в `runs/<run>/_lost/a<N>/<job>/` → один повтор (`--retries 1`) **в префикс новой попытки**
  `…/gtu/<run>/_attemptN/<job>/` и в другую зону (прерванная — в конце очереди на 30 мин); дальше итог `lost`.
  Итог задачи — только из одной попытки (`runs/<run>/<job>/`), метка чужой попытки не принимается: результаты не
  смешиваются и не задваиваются. Файлы прерванной попытки в S3 не удаляются (у ключа нет `DeleteObject`).
  Вручную — `runner.py collect --run-id <id> --job <имя> --attempt N`.
- Выгрузка при прерывании не гарантирована: ловушка `TERM` запускает выгрузку, успеет ли она — не проверено.
  Ранней выгрузки из команды нет: `S3_BUCKET` и `S3_KEY` в user-data не экспортируются. Длинную задачу — делить
  на части или писать контрольные точки в `$GTU_OUT`. Сторож уведомления внутри команды (в `runner.py` не
  встроен; проверен локально на подменённом `curl`, на машине — нет):
  ```bash
  ( while sleep 5; do
      T=$(curl -sf -X PUT http://169.254.169.254/latest/api/token -H "X-aws-ec2-metadata-token-ttl-seconds: 21600")
      if N=$(curl -sf -H "X-aws-ec2-metadata-token: $T" http://169.254.169.254/latest/meta-data/spot/instance-action)
      then printf '%s\n' "$N" > "$GTU_OUT/_spot_notice.json"; break; fi
    done ) &
  SPOT_W=$!
  # ... основная работа задачи; появился _spot_notice.json — сохранить контрольную точку ...
  kill "$SPOT_W" 2>/dev/null
  ```
  Запись вида `curl … > файл && break` не годится: файл создаётся пустым и при 404.
- `MaxPrice` держит верхнюю границу бюджета, но по документации EC2 повышает частоту прерываний. Частые
  прерывания — поднять `--price-margin` (бюджет учтёт), а не убирать `MaxPrice` без правки формулы.

## 8. Обязательные правила

- Ключи — только в переменных окружения Coldun: не в файлах, командах задач, user-data, журнале, чате и
  выводе; исполнитель отказывает, найдя ключ. Фиктивное значение в тесте — со словом SYNTHETIC, FIXTURE,
  EXAMPLE или PLACEHOLDER (сканер `gate.yml`).
- Живой запуск — только `--live` (или `GTU_LIVE=1`) и с разрешения владельца; программный `Config.dryrun`
  по умолчанию живой — не вызывать `Runner` из своего кода без `dryrun=True`.
- Запись — только под префиксом `ALPHA_O_S3`, подпапка `gtu/<run_id>/`; бакет один (us-east-1). Машины — в us-east-1,
  в соседних регионах — только флагом `--regions` (решение владельца 15.11.2026), по требованию — только
  `--on-demand-fallback`.
- Гасить только свой `run_id`. Прежний сборщик `aws_orchestrate.py` по одному `Project` убил соседний прогон
  fe4s4 (сейчас он только под `ALPHA_O_REAP_ALL=1`) — не повторять.
- Только владелец: политики IAM, служебная роль `AWSServiceRoleForEC2Spot`, квоты, Budgets, удаление
  объектов S3 и ключа доступа. Изменение сети и переменных окружения — сначала строка владельцу.
- Местный расчёт — не дольше 20 минут, длиннее — сюда. В дневник — только числа из скачанных results-JSON
  с `run_id`; предварительное помечать.

## 9. Типовые ошибки

| симптом | причина и действие |
|---|---|
| `UnauthorizedOperation` в DryRun | не применено дополнение политики: владелец — шаг 3 порядка проверки README |
| «ветки … нет» | ветка не запушена; машина клонирует только публичный репозиторий |
| `MaxSpotInstanceCountExceeded`, `VcpuLimitExceeded`, `InstanceLimitExceeded` | квота: задача ждёт, пока освободится машина прогона; своих машин нет (квоту держат чужие) — пауза и повтор до `--quota-wait-polls` (15), затем итог `quota`; поднять квоту — владелец |
| `InsufficientInstanceCapacity`, `Unsupported` | нет ёмкости или конфигурации в зоне: следующая зона по цене, затем `c7i.16xlarge`; нигде нет — повтор через `--poll-s` до `--capacity-wait-polls` (15), затем итог `capacity` с перечнем `тип@зона:код` в `detail` |
| итог `no_price` | у типов нет истории цены спота за окно (тип не продаётся, опечатка) — машина вслепую не запускается |
| итог `halted` | гашение машины прогона не подтверждено, очередь остановлена: погасить по строке `[гашение][ОШИБКА]` или `reap --run-id` |
| итог `failed`, `detail: RunInstances: <код>` | иная ошибка запуска (права, параметры); машина, созданная за ошибкой 5xx, находится по `ClientToken` и учитывается как запущенная |
| фаза `container`, код 127 / фаза `done`, код 127 в контейнере | в образе нет `bash` / команда не найдена — окружение образа не подключено (кроме `/openfoam/bash.rc`) |
| `SpotMaxPriceTooLow` | цена выше `MaxPrice`: тоже следующая зона, затем запасной тип (код по документации не сверен) |
| `image_unavailable` | образ `ubuntu2404` не разрешён ключу (SSM Canonical, `DescribeImages`): `image: al2023` или `container: ubuntu:24.04` с `apt-get` в `container_setup` |
| `[подсети] describe_subnets: UnauthorizedOperation` | нет прав — зона задаётся `Placement`; не ошибка |
| `[квоты] <регион> спот: 0 vCPU … регион не используется` | версия 3: квоту соседнего региона не дали ни флаг, ни service-quotas — `--vcpu-quota-per-region <регион>=<vCPU>` (не выше фактической) |
| `[регион] <регион> …: отказ в правах … исключён до конца прогона` | версия 3: ключ агента не может запускать (или читать цену) в регионе; остальные регионы работают; запись `ev=region_denied` |
| `образ al2023 в <регион> недоступен` | версия 3: `ssm:GetParameter` вне us-east-1 запрещён, а `ec2:DescribeImages` тоже не дал образ — регион пропускается |
| итог `capacity` / `quota` с перечнем регионов в `detail` | версия 3: ни один регион не дал машину; `--on-demand-fallback` добавит машины по требованию (если квота по требованию есть) |
| `cheapest_zone: false` в итоге | в дешёвой зоне не было ёмкости (`launch_error` с зоной) или зона после прерывания |
| фаза в `_DONE.json`: `base`/`venv`/`packages`/`clone` | сбой подготовки, коды 10–17; читать `_job.log` |
| фаза `venv`, код 18 / 19; `requirements`, 23 | uv не поставился / Python не скачан; файла зависимостей нет в ветке |
| фаза `docker` 21, `pull` 22, `container` 20 | docker не запустился / образ не скачан (429 Docker Hub — зеркало `public.ecr.aws/docker/library/…`) / упала команда `container_setup` (`_container_setup.log`) |
| `no_time` (код 124), `command` с кодом 124 или 137 | не хватило `max_min` (137 — `timeout` добил KILL): поднять срок или делить задачу |
| user-data больше 16 КБ | длинную команду — в скрипт репозитория, в задаче — его вызов; проверяется для всех попыток (повтор длиннее) до первого запуска |
| в журнале прогона нет строк | вывод в файл буферизован: запускать `python3 -u` |
| прогон оборвался через 30 мин или 2 ч | остановлена фоновая команда инструмента Bash по `timeout`: запускать через `setsid` (шаг 5) |
| `status` начал счёт с нуля | `ledger.jsonl` не закоммичен в прошлый раз: `run --spent-before` с расходом из консоли |
| пакет apt устарел (OpenFOAM 1912 на Ubuntu 24.04) | ставить нужную версию командой задачи, см. `openfoam-combustor` |

## 10. Проверка результата

- Итог `run`: у каждой задачи `status: ok`; в `_DONE.json` `phase: done`, `exit_code: 0`.
- `smoke`: `nproc.txt` = 64, `commit.txt` совпадает с последним коммитом запушенной ветки, в `_job.log` время
  подготовки (строки `[gtu-time] <фаза> <с>s`).
- `jobs.smoke2.json`: у обеих задач `cheapest_zone: true`, `az` = `az_cheapest` = `az_machine` в `summary.json`,
  `"ok": true` в `az_check.json`; `smoke_py.json` — Python 3.12.x, версии cantera/numpy/scipy, точка цикла и
  `cycle_ref_rel_diff`, пул на `NPROC`; `smoke_ccx.json` — `ccx_version` (код 201 у `ccx -v` — норма) и диск. Код 3 —
  счёт прошёл, машина не в самой дешёвой зоне (разбор — `az_check.json`, `launch_error` в журнале).
- После прогона: `runner.py status` — открытых машин 0; `runner.py reap --run-id <id> --list-only` — живых нет.
- Повторяемость: `explore-249` сравнить с `routes/gtu/gtu_explore_results.json`; расхождение — разбор, не публикация.

## 11. Ссылки трека

`routes/gtu/platform/aws/` (`runner.py`, `test_runner.py`, `README.md`, `jobs.example.json`,
`iam_policy_addendum.json`, `aws_prices_ref.json`); `calc/aws_orchestrate.py`; `calc/aws_setup/` (`setup.sh`, `README.md`, политики);
`routes/gtu/gtu_cycle.py`, `gtu_explore.py`, `gtu_dropin.py` и их `*_results.json`; навыки `gtu-cycle`,
`gtu-dropin`, `su2-turbo`, `openfoam-combustor`, `calculix-rotor`, `surrogate-opt`.

## 12. Соседние мощности (версия 3)

Решение владельца 15.11.2026: «можешь занять ещё соседние мощности»; потолок 2000 долл. со стоп-краном и отчётом
каждые 500 — без изменений. Регионы и запасной путь — флаги запуска; файлы задач не меняются. Без новых флагов —
точно версия 2 (те же вызовы AWS и журнал; проверено 94 прежними тестами и дифференциальным сравнением с кодом
версии 2: 22 сценария заглушки, `plan` по 37 файлам задач, `status` по живому журналу).

Проверка сессии с ключами 15.11.2026 10:45 UTC (без трат), спот c7a.16xlarge, долл./ч:

| регион | спот, vCPU | по требованию, vCPU | лучшие зоны | SSM | DryRun |
|---|---|---|---|---|---|
| us-east-1 | 384 (заявка 1536) | 384 | 1c 0,703 · 1b 0,747 · 1d 0,879 | да | да |
| us-east-2 | 5 (заявка 1536) | 5 | 2b 0,601 · 2a 0,606 · 2c 0,940 | AccessDenied | да |
| us-west-2 | 384 | 300 | 2a 0,866 · 2d 0,952 · 2c 1,422 | AccessDenied | да |
| eu-north-1 | 5 | 5 | 1a 0,601 | AccessDenied | да |

Сегодня это значит: +6 машин спота в us-west-2; us-east-2 и eu-north-1 — после одобрения заявок (с квотой 5 vCPU
машина 64 vCPU не помещается, раннер туда не просится); `--on-demand-fallback` — ещё до 10 машин по требованию
(6 в us-east-1, 4 в us-west-2), только когда спот не дал машину нигде.

Как работает (подробно — README узла, раздел «Соседние мощности»):
- **квота региона** — service-quotas (`L-34B43A08` спот, `L-1216C47A` по требованию); флаг
  `--vcpu-quota-per-region` / `--od-vcpu-quota-per-region` только ограничивает сверху; service-quotas не ответил —
  флаг, иначе запасное: us-east-1 — `--max-vcpu`, соседний — `--vcpu-quota-fallback` (0 — не использовать). Есть ли у
  ключа агента право `servicequotas:GetServiceQuota` — неизвестно: задавать квоты флагом по таблице выше;
- **регион задачи** — наименьшая цена спота среди регионов со свободной квотой; соседний — с запасом
  `--region-price-margin` (0,05); регион с отказом в правах исключается, с квотной ошибкой AWS — на паузу;
- **MaxPrice и бюджет** — по максимуму цены за окно по всем разрешённым регионам; передача результатов в бакет
  us-east-1 из соседнего региона — 0,02 долл./ГБ × `min(disk_gb, --xfer-gb 10)` в смете (`xfer_usd` в журнале);
- **образ** — SSM своего региона, при отказе — `DescribeImages` (amazon, `al2023-ami-2023*-x86_64`, самый свежий);
- **по требованию** — без `InstanceMarketOptions`, тег `Market=on-demand`, гашение — те же три способа; в бюджете
  — максимум цены справочника `aws_prices_ref.json` по регионам × `--price-margin` (MaxPrice у такой машины нет).
  Справочник — **предварительный** (выдержки поиска: c7a.16xlarge 3,2845 в трёх регионах США, c7i.16xlarge 2,856
  только us-east-1): сверить со страницей цен AWS перед первым живым прогоном по требованию;
- **журнал** пишет `region` и `market` (запуск, окончание, ошибки, DryRun), проверка бюджета — `regions`;
  `status` — разбивка по регионам и рынкам; `reap` и сверка — по всем регионам прогона из журнала; `collect` — как
  было (бакет один). Журнал после прогона с `--regions` читать только раннером версии 3.

```bash
python3 routes/gtu/platform/aws/runner.py plan <jobs.json> --assume-price 0.95 --regions us-east-1,us-east-2,us-west-2 \
        --assume-zone-prices us-east-1c=0.703,us-east-2b=0.601,us-west-2a=0.866 \
        --vcpu-quota-per-region us-east-1=384,us-east-2=5,us-west-2=384          # регион каждой задачи и волны
python3 routes/gtu/platform/aws/runner.py run <jobs.json> --branch <ветка> --regions us-east-1,us-east-2,us-west-2 \
        --vcpu-quota-per-region us-east-1=384,us-east-2=5,us-west-2=384          # DryRun: права в каждом регионе
# живой: + --on-demand-fallback --od-vcpu-quota-per-region us-east-1=384,us-west-2=300 --live (п. 2 — разрешение)
```

Не проверено вживую: запуск в us-west-2 и us-east-2, `DescribeImages` и service-quotas у ключа агента, разрешён
ли ключу запуск по требованию (политика писалась под спот), цены справочника.

## Источники (проверено 14.11.2026; версия 3 — 15.11.2026)

docs.aws.amazon.com из контейнера не открывается; параметры API сверены по модели API в botocore 1.44.3,
имена ключей CLI — по awscli 1.46.3 (разбор команд без ключей, до отказа «Unable to locate credentials»),
поведение — по официальному архиву руководства EC2 (снимок 15.06.2023, может отставать от текущей документации).

- botocore `ec2` service-2: `SpotMarketOptions` (one-time/persistent, terminate, предупреждение о `MaxPrice`),
  `InstanceMetadataOptionsRequest` (hop limit 1–64), `RunInstances` (`InstanceInitiatedShutdownBehavior`,
  `DryRun` → `DryRunOperation`/`UnauthorizedOperation`), фильтры `tag:<key>`; `budgets` `CreateBudget` (тело
  команды из п. 6 прошло `ParamValidator`); `service-quotas` `GetServiceQuota` — https://pypi.org/project/botocore/1.44.3/
- awscli 1.46.3, `awscli/examples/service-quotas/get-service-quota.rst` — код и имя квоты `L-1216C47A`.
- Архив руководства EC2 (https://github.com/awsdocs/amazon-ec2-user-guide, коммит e3707c4, `doc_source/`):
  `spot-interruptions.md` (2 минуты), `spot-instance-termination-notices.md`, `rebalance-recommendations.md`,
  `terminating-instances.md` (`halt`), `instancedata-add-user-data.md` (16 КБ до base64), `using-spot-limits.md`,
  `troubleshooting-launch.md` (`InstanceLimitExceeded`, `InsufficientInstanceCapacity`, `Unsupported`),
  `spot-requests.md` (`AWSServiceRoleForEC2Spot`), `ExamplePolicies_EC2.md` (`spot-instances-request`)
- systemd `shutdown` (`+m`, `-P`, `-c`): https://github.com/systemd/systemd/blob/main/man/shutdown.xml
- Образы Ubuntu в SSM: https://documentation.ubuntu.com/aws/aws-how-to/instances/find-ubuntu-images/
- Облачные окружения Claude Code (переменные, сон ВМ): https://code.claude.com/docs/en/cloud-environments
- OpenFOAM 1912 в Ubuntu 24.04: `apt-cache policy openfoam` (noble/universe, 1912.200626-2build3)
- Версия 3: квоты, цены спота, SSM и DryRun по регионам — проверка сессии с ключами 15.11.2026 10:45 UTC (без трат);
  цены по требованию — выдержки поиска 15.11.2026 (calculator.holori.com, aws-pricing.com), страница цен AWS и
  агрегаторы из контейнера агента не открываются — `aws_prices_ref.json`, статус «предварительно»
- 64 vCPU у c7a/c7i.16xlarge — выдержки поиска по сторонним каталогам
  (https://instances.vantage.sh/aws/ec2/c7a.16xlarge); страница типов AWS недоступна из контейнера

**Не проверено:** квоты 384 vCPU в аккаунте и код `L-34B43A08`; коды `SpotMaxPriceTooLow`,
`MaxSpotInstanceCountExceeded`, `VcpuLimitExceeded` (в архиве нет; `VcpuLimitExceeded` упомянут в
`aws_orchestrate.py` как случившийся) и код причины гашения при прерывании; ресурс `spot-instances-request/*` —
по архиву 2023 г. проверяется, только если запрос спота помечается тегом при создании (`runner.py` его не
помечает), по текущей документации не сверено; есть ли в аккаунте `AWSServiceRoleForEC2Spot` (по архиву для
запуска спота через API она должна существовать); ставки gp3 и IPv4; успевает ли выгрузка при прерывании;
установка пакетов и `venv` на Ubuntu; не даёт ли фоновый процесс уснуть ВМ сессии; пуш из дочерней сессии.
Версия 2 вживую: uv и Python 3.12 на al2023; docker из `dnf` на al2023, скачивание образа с Docker Hub,
`--shm-size` в половину памяти; `Placement` без `SubnetId` (подсеть по умолчанию зоны — по модели API botocore);
права ключа на `DescribeSubnets` и `DescribeImages`; метка 143 при прерывании; квотные коды при очереди. На
al2023 (живой прогон 14.11.2026) работали `dnf`, `python3 -m venv`, boto3 и выгрузка.
