---
name: aws-gtu-runs
description: Прогоны трека «Простая литая ГТУ-заменитель» на AWS через доступ MatterForge скриптом routes/gtu/platform/aws/runner.py — где действуют ключи (окружение Coldun, дочерняя сессия), us-east-1, квоты 384 vCPU, спот c7a.16xlarge с запасным c7i.16xlarge, бакет из ALPHA_O_S3 с подпапкой gtu/, IMDSv2, автогашение тремя способами, бюджетный предохранитель 300 долл. и журнал ledger.jsonl, пакет задач, сбор и проверка результатов, прерывание спота, запрет ключей в файлах и выводе. Используй, когда ставишь расчёт SU2, OpenFOAM, CalculiX или скрипт трека на облачные машины, проверяешь права и бюджет, гасишь машины прогона, забираешь результаты или разбираешь сбой прогона. Триггеры — «прогон на AWS», «запусти на споте», «runner.py», «облачный расчёт ГТУ», «Coldun», «ALPHA_O_S3», «погаси машины», «бюджет AWS», «прерывание спота», «ledger.jsonl».
---

# Прогоны трека ГТУ на AWS

Дорогие расчёты платформы (SU2 — ступени, OpenFOAM — камера, CalculiX — диски, `gtu_explore.py` и любой
скрипт репозитория) идут пакетом независимых задач на спотовых машинах us-east-1. Исполнитель —
`routes/gtu/platform/aws/runner.py`; основа — `calc/aws_orchestrate.py` и `calc/aws_setup/`. Состояние на
14.11.2026: боевых запусков не было (README узла), тесты на заглушках — 25 passed (прогнано 14.11.2026).

## 1. Что где лежит

| что | где |
|---|---|
| исполнитель: `plan`, `run` (DryRun по умолчанию, `--live`), `status`, `reap`, `collect` | `routes/gtu/platform/aws/runner.py` |
| тесты без сети (botocore Stubber), пример пакета `smoke` + `explore-249` | `test_runner.py`, `jobs.example.json` там же |
| порядок проверки, формат задач, формула бюджета, откат | `routes/gtu/platform/aws/README.md` |
| права спота (дополнение, применяет владелец) | `routes/gtu/platform/aws/iam_policy_addendum.json` |
| журнал расходов (коммитить) / скачанные результаты (в `.gitignore`) | `ledger.jsonl` / `runs/<run_id>/<job>/` там же |
| разовая настройка аккаунта (владелец, админ-ключи, не агент) | `calc/aws_setup/setup.sh`, `user_policy.json`, `instance_role_policy.json` |
| предшественник: одна машина on-demand для alpha-O, не для ГТУ | `calc/aws_orchestrate.py` |
| флот eu-central-1 других треков — не трогать | `calc/fleet_status.py` |

## 2. Где действуют ключи

- **Только в окружении Coldun.** Сверено 14.11.2026 (`get_environment`): личное окружение, сеть без ограничений,
  переменные `AWS_ACCESS_KEY_ID`, `AWS_SECRET_ACCESS_KEY`, `AWS_DEFAULT_REGION`, `ALPHA_O_S3`,
  `ALPHA_O_INSTANCE_PROFILE`, `ALPHA_O_TAG` (видны только имена), скрипта настройки и репозиториев нет.
- В любом другом контейнере ключи не действуют: только тесты и `plan`, вызовов AWS не делать.
- Идущую сессию в Coldun не перевести (окружение с неограниченной сетью) — нужна **дочерняя сессия**:
  `list_environments` → id окружения с именем Coldun (id в файлы не писать) → `create_session` с
  `environment_id`, `source_url` репозитория, `source_revision` = запушенная ветка, задачей в `prompt`;
  `permission_mode` не `plan` (иначе дочерняя ждёт одобрения). Ход — `get_session` (`status_bucket`) и
  `list_events`; чисто закончившаяся сессия может не отчитаться сама.
- Переменные окружения читает любой, кто им пользуется (документация облачных окружений); изменённые
  значения сессия видит после перезапуска своей машины.
- Доступ проверять вызовом API (`aws sts get-caller-identity` или строка «вызов от имени» в `run`), не
  печатая значений: ни длины, ни начала ключа. В имени бакета номер аккаунта — в выводе `<BUCKET>`.

## 3. Параметры, которые нельзя менять молча

| параметр | значение | откуда |
|---|---|---|
| регион | us-east-1 | `REGION` в коде; eu-central-1 занят флотом `fleet_status.py` |
| квоты | 384 vCPU on-demand и спот (данные владельца) | `--vcpu-quota 384` = 6 машин × 64; `--max-parallel 6` |
| типы | `c7a.16xlarge`, запас `c7i.16xlarge` (по 64 vCPU) | `DEFAULT_TYPE`, `DEFAULT_FALLBACK`, `VCPU` |
| спот | `MarketType=spot`, `SpotInstanceType=one-time`, `InstanceInterruptionBehavior=terminate`, `MaxPrice` | `build_params` |
| IMDSv2 | `HttpTokens=required`, `HttpPutResponseHopLimit=1`, `HttpEndpoint=enabled` | на машине нет ключей, S3 — роль |
| запись | `s3://<бакет>/<префикс>/gtu/<run_id>/<job>/` | роль машины: только `s3:PutObject` в `<префикс>/*` |
| чтение | `ListBucket` по префиксу, `GetObject`; `DeleteObject` нет | `user_policy.json` |
| теги | `Project=<ALPHA_O_TAG>`, `RunId`, `Job`, `Name=gtu-<job>` на машине и диске | политика гасит только машины с тегом `Project` |
| потолок | 300 долл. (`--ceiling`) | журнал `ledger.jsonl` |

Коды квот (по выдержкам поиска из документации EC2): `L-1216C47A` — Running On-Demand Standard
(A, C, D, H, I, M, R, T, Z) instances, `L-34B43A08` — All Standard Spot Instance Requests, обе в vCPU.
Прав Service Quotas у ключа агента нет — смотрит владелец:
`aws service-quotas get-service-quota --region us-east-1 --service-code ec2 --quota-code L-34B43A08`.

## 4. Порядок работы

**Шаг 0.** Начало сессии по `CLAUDE.md`; прочитать README узла и навык расчёта (`su2-turbo`,
`openfoam-combustor`, `gtu-cycle`, `surrogate-opt`).

**Шаг 1. Пакет задач.** JSON `{"defaults": {...}, "jobs": [...]}`; поля задачи: `name` (`[A-Za-z0-9._-]`),
`command` (bash, из корня клона), `instance_type`, `fallback_types`, `max_min` (20…720 — срок сторожа),
`image` (`ubuntu2404` или `al2023`), `disk_gb`, `apt`/`dnf`, `pip`, `vcpu`. Неизвестное поле — отказ. Команда
пишет всё в `$GTU_OUT`; доступны `GTU_RUN_ID`, `GTU_JOB`, `GTU_ATTEMPT`, `OMP_NUM_THREADS=$(nproc)`.

**Шаг 2. Офлайн-проверка (любой контейнер).**
```bash
python3 -m pytest -q -p no:cacheprovider routes/gtu/platform/aws/test_runner.py   # 25 passed
python3 routes/gtu/platform/aws/runner.py plan <jobs.json> --assume-price 1.5 --show-userdata \
        --out <scratchpad>/runs --ledger <scratchpad>/ledger.jsonl                 # цена — допущение
```
Прочитать user-data и параметры `run_instances`; прогноз пакета при допущенной цене — не число для отчёта.

**Шаг 3. Ветка запушена.** Машина клонирует публичный репозиторий (`git clone --depth 1 -b <ветка>`);
`run` сверяет ветку через `git ls-remote`.

**Шаг 4. Проверка прав (Coldun).** `python3 -c 'import boto3'` (нет — `pip install boto3`), затем
`python3 routes/gtu/platform/aws/runner.py run <jobs.json> --branch <ветка>` без `--live`. Ожидается
`запуск разрешён (DryRunOperation)` по каждой задаче и `чтение префикса результатов: OK`.

**Шаг 5. Живой прогон** — только с явного разрешения владельца в этой сессии. Первый — один `smoke`
(`max_min` 20), затем рабочий пакет. Долгий — в фоне с журналом в файл:
`python3 routes/gtu/platform/aws/runner.py run <jobs.json> --branch <ветка> --live > /tmp/gtu_run.log 2>&1 &`.
`run_id` — первая строка журнала и имя папки `runs/<run_id>/`; записать его сразу.

**Шаг 6. Сбор.** По метке `_DONE.json` исполнитель сам скачивает `<префикс>/gtu/<run_id>/<job>/` в
`runs/<run_id>/<job>/` и гасит машину. Вручную: `runner.py collect --run-id <id> --job <имя>`.

**Шаг 7. Учёт.** `runner.py status` (открытых машин 0), затем `ledger.jsonl` и проверенные results-JSON — в
ветку и PR по правилам `CLAUDE.md`; контейнер эфемерен, незакоммиченный журнал пропадает.

## 5. Автогашение — три независимых способа

1. **Выключение = гашение.** `InstanceInitiatedShutdownBehavior=terminate` (по умолчанию у EC2 — `stop`) и
   разовый спот с `terminate` при прерывании. `halt` машину не выключает и она продолжает работать —
   в командах только `poweroff` или `shutdown -P`.
2. **Сторож на машине.** Первая строка user-data — `shutdown -P +<max_min>` (`+m` — через m минут); ловушка
   `EXIT` выгружает `$GTU_OUT`, `_job.log`, `_stdout.log`, последней — `_DONE.json`, затем `poweroff`.
   Команде задачи отводится `max_min` минус 8 мин на выгрузку.
3. **Внешний сборщик.** В `finally` (конец, ошибка, Ctrl-C, SIGTERM, SIGHUP) — гашение машин с тегами
   `Project` **и** `RunId` этого прогона; срок сборщика — `max_min + grace_min` (10). Если сессия умерла,
   1 и 2 всё равно гасят; остатки — `runner.py reap --run-id <id>` (`--list-only` — только показать) или
   вручную командами из раздела «Откат» README (`aws ec2 describe-instances --filters
   Name=tag:Project,... Name=tag:RunId,...` → `aws ec2 terminate-instances --instance-ids ...`).

## 6. Бюджетный предохранитель

Перед каждым запуском: `израсходовано + в работе + новая ≤ 300`. Новая = `(MaxPrice + диск) × (max_min + 10) / 60`,
`MaxPrice` = максимум цены спота за 6 ч по всем зонам × 1,2 (`--price-margin`). Открытые машины считаются по
верхней границе, закрытые — по цене спота в своей зоне и фактическому времени (не меньше 60 с). Ставка gp3
0,08 долл./ГБ·мес — допущение кода. Отказ пишется в журнал (`ev=budget`, `ok=false`), машина не создаётся.
Расход вне журнала — `--spent-before`. Внешняя страховка — бюджет AWS Budgets, ставит владелец (у ключа
агента прав нет): `aws budgets create-budget --account-id <аккаунт> --budget
'{"BudgetName":"gtu","BudgetLimit":{"Amount":"300","Unit":"USD"},"TimeUnit":"MONTHLY","BudgetType":"COST"}'`.

## 7. Прерывание спота

- AWS предупреждает примерно за 2 минуты: на машине появляется `/latest/meta-data/spot/instance-action`
  (до этого — 404), в EventBridge — событие «EC2 Spot Instance Interruption Warning»; раньше может прийти
  рекомендация `events/recommendations/rebalance`.
- Исполнитель: машина `terminated` без `_DONE.json` → скачать, что успело лечь в S3, → один повтор
  (`--retries 1`) в тот же префикс, номер попытки — в `_DONE.json`; дальше итог `lost`.
- Выгрузка при прерывании не гарантирована: ловушка `TERM` запускает выгрузку, успеет ли она за 2 минуты —
  не проверено. Ранней выгрузки из команды нет: `S3_BUCKET` и `S3_KEY` в user-data не экспортируются.
  Длинную задачу — делить на части или писать контрольные точки в `$GTU_OUT`. Сторож уведомления
  внутри команды (на машине не проверялся, в `runner.py` не встроен):
  ```bash
  ( while sleep 5; do T=$(curl -s -X PUT http://169.254.169.254/latest/api/token \
        -H "X-aws-ec2-metadata-token-ttl-seconds: 21600")
      curl -sf -H "X-aws-ec2-metadata-token: $T" http://169.254.169.254/latest/meta-data/spot/instance-action \
        > "$GTU_OUT/_spot_notice.json" && break; done ) &
  ```
- `MaxPrice` держит верхнюю границу бюджета, но по документации EC2 повышает частоту прерываний. Частые
  прерывания — поднять `--price-margin` (бюджет учтёт), а не убирать `MaxPrice` без правки формулы.

## 8. Обязательные правила

- Ключи — только в переменных окружения Coldun: не в файлах, командах задач, user-data, журнале, чате и
  выводе; исполнитель отказывает, найдя ключ. Фиктивное значение в тесте — со словом SYNTHETIC или EXAMPLE.
- Живой запуск — только `--live` (или `GTU_LIVE=1`) и с разрешения владельца; программный `Config.dryrun`
  по умолчанию живой — не вызывать `Runner` из своего кода без `dryrun=True`.
- Запись — только под префиксом `ALPHA_O_S3`, подпапка `gtu/<run_id>/`; регион — только us-east-1.
- Гасить только свой `run_id`. Широкий сборщик по одному `Project` уже убил чужой прогон
  (`ALPHA_O_REAP_ALL` в `aws_orchestrate.py`) — не повторять.
- Только владелец: политики IAM, служебная роль `AWSServiceRoleForEC2Spot`, квоты, Budgets, удаление
  объектов S3 и ключа доступа. Изменение сети и переменных окружения — сначала строка владельцу.
- Местный расчёт — не дольше 20 минут, длиннее — сюда. В дневник — только числа из скачанных results-JSON
  с `run_id`; предварительное помечать.

## 9. Типовые ошибки

| симптом | причина и действие |
|---|---|
| `UnauthorizedOperation` в DryRun | не применено дополнение политики: владелец — шаг 3 порядка проверки README |
| «ветки … нет» | ветка не запушена; машина клонирует только публичный репозиторий |
| `MaxSpotInstanceCountExceeded`, `VcpuLimitExceeded`, `InstanceLimitExceeded` | квота: задача ждёт свободной машины; поднять квоту — владелец |
| `InsufficientInstanceCapacity`, `Unsupported` | нет ёмкости или конфигурации в зоне: исполнитель берёт `c7i.16xlarge` |
| `SpotMaxPriceTooLow` | цена выше `MaxPrice`: тоже переход на запасной тип (код по документации не сверен) |
| фаза в `_DONE.json`: `base`/`venv`/`packages`/`clone` | сбой подготовки, код 10–17; читать `_job.log` |
| `no_time` (код 124) или `command` с кодом 124 | не хватило `max_min`: поднять срок или делить задачу |
| user-data больше 16 КБ | длинную команду — в скрипт репозитория, в задаче — его вызов |
| `status` начал счёт с нуля | `ledger.jsonl` не закоммичен в прошлый раз: `--spent-before` с расходом из консоли |
| пакет apt устарел (OpenFOAM 1912 на Ubuntu 24.04) | ставить нужную версию командой задачи, см. `openfoam-combustor` |

## 10. Проверка результата

- Итог `run`: у каждой задачи `status: ok`; в `_DONE.json` `phase: done`, `exit_code: 0`.
- `smoke`: `nproc.txt` = 64, `commit.txt` совпадает с запушенной веткой, в `_job.log` время подготовки.
- После прогона: `runner.py status` — открытых машин 0; `runner.py reap --run-id <id> --list-only` — живых нет.
- Повторяемость: `explore-249` сравнить с `routes/gtu/gtu_explore_results.json`; расхождение — разбор, не публикация.

## 11. Ссылки трека

`routes/gtu/platform/aws/` (`runner.py`, `test_runner.py`, `README.md`, `jobs.example.json`,
`iam_policy_addendum.json`); `calc/aws_orchestrate.py`; `calc/aws_setup/` (`setup.sh`, `README.md`, политики);
`routes/gtu/gtu_cycle.py`, `gtu_explore.py`, `gtu_dropin.py` и их `*_results.json`; навыки `gtu-cycle`,
`gtu-dropin`, `su2-turbo`, `openfoam-combustor`, `surrogate-opt`.

## Источники (проверено 14.11.2026)

docs.aws.amazon.com из контейнера не открывается; параметры API сверены по модели API в botocore 1.44.3
(выпуск на PyPI 13.11.2026), имена ключей CLI — по awscli 1.46.3, поведение — по официальному архиву
руководства EC2 (снимок 15.06.2023, может отставать от текущей документации).

- botocore `ec2` service-2: `SpotMarketOptions` (one-time/persistent, terminate, предупреждение о `MaxPrice`),
  `InstanceMetadataOptionsRequest` (hop limit 1–64), `RunInstances` (`InstanceInitiatedShutdownBehavior`,
  `DryRun` → `DryRunOperation`/`UnauthorizedOperation`), фильтры `tag:<key>`; `budgets` `CreateBudget`;
  `service-quotas` `GetServiceQuota` — https://pypi.org/project/botocore/1.44.3/
- Архив руководства EC2 (https://github.com/awsdocs/amazon-ec2-user-guide, коммит e3707c4):
  `spot-instance-termination-notices.md`, `spot-best-practices.md`, `terminating-instances.md` (`halt`), `instancedata-add-user-data.md` (16 КБ),
  `using-spot-limits.md`, `troubleshooting-launch.md`, `spot-requests.md` (`AWSServiceRoleForEC2Spot`),
  `ExamplePolicies_EC2.md` (`spot-instances-request/*`), `rebalance-recommendations.md`
- Квоты (выдержки поиска): https://docs.aws.amazon.com/AWSEC2/latest/UserGuide/using-spot-limits.html,
  https://docs.aws.amazon.com/AWSEC2/latest/UserGuide/ec2-on-demand-instances.html
- systemd `shutdown` (`+m`, `-P`, `-c`): https://github.com/systemd/systemd/blob/main/man/shutdown.xml
- Образы Ubuntu в SSM: https://documentation.ubuntu.com/aws/aws-how-to/instances/find-ubuntu-images/
- Облачные окружения Claude Code: https://code.claude.com/docs/en/cloud-environments
- 64 vCPU у c7a/c7i.16xlarge — выдержки поиска по сторонним каталогам
  (https://instances.vantage.sh/aws/ec2/c7a.16xlarge); страница типов AWS недоступна из контейнера

**Не проверено:** квоты 384 vCPU в аккаунте; код `SpotMaxPriceTooLow` и код причины гашения при прерывании;
нужны ли ресурс `spot-instances-request/*` и служебная роль именно этому аккаунту; ставка gp3; успевает ли
выгрузка за 2 минуты; установка пакетов и `venv` на реальных образах; `curl` в образе Ubuntu.
