---
name: aws-gtu-runs
description: Прогоны трека «Простая литая ГТУ-заменитель» на AWS через доступ MatterForge скриптом routes/gtu/platform/aws/runner.py — где действуют ключи (окружение Coldun, дочерняя сессия), us-east-1, квоты 384 vCPU, спот c7a.16xlarge с запасным c7i.16xlarge, бакет из ALPHA_O_S3 с подпапкой gtu/, IMDSv2, автогашение тремя способами, бюджетный предохранитель 300 долл. и журнал ledger.jsonl, пакет задач, сбор и проверка результатов, прерывание спота, запрет ключей в файлах и выводе. Используй, когда ставишь расчёт SU2, OpenFOAM, CalculiX или скрипт трека на облачные машины, проверяешь права и бюджет, гасишь машины прогона, забираешь результаты или разбираешь сбой прогона. Триггеры — «прогон на AWS», «запусти на споте», «runner.py», «облачный расчёт ГТУ», «Coldun», «ALPHA_O_S3», «погаси машины», «бюджет AWS», «прерывание спота», «ledger.jsonl».
---

# Прогоны трека ГТУ на AWS

Дорогие расчёты платформы (SU2 — ступени, OpenFOAM — камера, CalculiX — диски, `gtu_explore.py` и любой
скрипт репозитория) идут пакетом независимых задач на спотовых машинах us-east-1. Исполнитель —
`routes/gtu/platform/aws/runner.py`; основа — `calc/aws_orchestrate.py` и `calc/aws_setup/`. Состояние на
14.11.2026: боевых запусков не было (README узла), тесты на заглушках — 32 passed (прогнано 14.11.2026).

## 1. Что где лежит

| что | где |
|---|---|
| исполнитель: `plan`, `run` (DryRun по умолчанию, `--live`), `status`, `reap`, `collect` | `routes/gtu/platform/aws/runner.py` |
| тесты без сети (botocore Stubber), пример пакета `smoke` + `explore-249` | `test_runner.py`, `jobs.example.json` там же |
| порядок проверки, формат задач, формула бюджета, откат | `routes/gtu/platform/aws/README.md` |
| права спота (дополнение, применяет владелец) | `routes/gtu/platform/aws/iam_policy_addendum.json` |
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
| регион | us-east-1 | `REGION` в коде; eu-central-1 занят флотом `fleet_status.py` |
| регион ручных команд | всегда `--region us-east-1` | значение `AWS_DEFAULT_REGION` в Coldun не читалось |
| квоты | 384 vCPU on-demand и спот (данные владельца) | `--vcpu-quota 384` = 6 машин × 64; `--max-parallel 6` |
| типы | `c7a.16xlarge`, запас `c7i.16xlarge` (по 64 vCPU) | `DEFAULT_TYPE`, `DEFAULT_FALLBACK`, `VCPU` |
| спот | `MarketType=spot`, `SpotInstanceType=one-time`, `InstanceInterruptionBehavior=terminate`, `MaxPrice` | `build_params` |
| IMDSv2 | `HttpTokens=required`, `HttpPutResponseHopLimit=1`, `HttpEndpoint=enabled` | на машине нет ключей, S3 — роль |
| запись | `s3://<бакет>/<префикс>/gtu/<run_id>/<job>/` | роль машины: только `s3:PutObject` в `<префикс>/*` |
| чтение | `ListBucket` по префиксу, `GetObject`; `DeleteObject` нет | `user_policy.json` |
| теги | `Project=<ALPHA_O_TAG>`, `RunId`, `Job`, `Name=gtu-<job>` на машине и диске | политика гасит только машины с тегом `Project` |
| потолок | 300 долл. (`--ceiling`) | журнал `ledger.jsonl` |

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
`image` (`ubuntu2404` или `al2023`), `disk_gb`, `apt`/`dnf`, `pip`, `vcpu`. Неизвестное поле — отказ. Команда
пишет всё в `$GTU_OUT`; доступны `GTU_RUN_ID`, `GTU_JOB`, `GTU_ATTEMPT`, `OMP_NUM_THREADS=$(nproc)`, `python` из
venv машины.

**Шаг 2. Офлайн-проверка (любой контейнер).**
```bash
python3 -m pytest -q -p no:cacheprovider routes/gtu/platform/aws/test_runner.py   # 32 passed
python3 routes/gtu/platform/aws/runner.py plan <jobs.json> --assume-price 1.5 --show-userdata \
        --out <scratchpad>/runs                                                    # цена — допущение
```
Прочитать user-data и параметры `run_instances`. `plan` в `ledger.jsonl` не пишет (проверки бюджета — в
`<out>/<run_id>/journal.dryrun.jsonl`); прогноз при допущенной цене — не число для отчёта.

**Шаг 3. Ветка запушена.** Машина клонирует публичный репозиторий (`git clone --depth 1 -b <ветка>`);
`run` сверяет ветку через `git ls-remote`.

**Шаг 4. Проверка прав (Coldun).** `python3 -c 'import boto3'` (нет — `pip install boto3`), затем
`python3 routes/gtu/platform/aws/runner.py run <jobs.json> --branch <ветка>` без `--live`. Ожидается
`запуск разрешён (DryRunOperation)` по каждой задаче и `чтение префикса результатов: OK`; код выхода 1, если
хоть одна задача не `dryrun_ok`. Цена спота и образ из SSM при этом читаются настоящими вызовами.

**Шаг 5. Живой прогон** — только с явного разрешения владельца (п. 2). Первый — один `smoke` в отдельном
файле задач (`max_min` 20), затем рабочий пакет:
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
   `Project` **и** `RunId` этого прогона; срок задачи у сборщика — `max_min + grace_min` (10). Если сессия
   умерла, 1 и 2 всё равно гасят; остатки — `runner.py reap --run-id <id>` (`--list-only` — только показать)
   или вручную командами из раздела «Откат» README:
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
итого    = израсходовано + в работе + пакет (только DryRun/plan) + новая  ≤  300 (--ceiling)
новая    = (MaxPrice + диск + IPv4) × (max_min + grace_min) / 60,   grace_min = 10
MaxPrice = максимум цены спота за 6 ч по всем зонам × 1,2 (--price-margin), вверх до 0,0001
диск     = disk_gb × 0,08 / 730 долл./ч;   IPv4 = 0,005 долл./ч   (обе ставки — допущения кода, не сверены)
```
Открытые машины — по верхней границе `(MaxPrice + диск + IPv4) × весь срок`; закрытые — по цене спота в своей
зоне на момент запуска и фактическому времени (не меньше 60 с) — это оценка, не граница. Трафик и запросы S3
не учтены. Сверка формулы: `plan` `jobs.example.json` при допущенной цене 1,5 даёт 0,91 (`smoke`) + 1,51
(`explore-249`) = 2,41 долл. Отказ пишется в журнал (`ev=budget`, `ok=false`), машина не создаётся. Расход вне
журнала — `run --spent-before`.

Внешняя страховка — AWS Budgets, ставит владелец (у ключа агента прав нет). Budgets только оповещает и машин
не гасит; без фильтров считает весь аккаунт, включая флот других треков:
```bash
aws budgets create-budget --account-id <аккаунт> \
  --budget '{"BudgetName":"gtu","BudgetLimit":{"Amount":"300","Unit":"USD"},"TimeUnit":"MONTHLY","BudgetType":"COST"}' \
  --notifications-with-subscribers '[{"Notification":{"NotificationType":"ACTUAL","ComparisonOperator":"GREATER_THAN","Threshold":80,"ThresholdType":"PERCENTAGE"},"Subscribers":[{"SubscriptionType":"EMAIL","Address":"<почта владельца>"}]}]'
```

## 7. Прерывание спота

- AWS предупреждает примерно за 2 минуты: на машине появляется `/latest/meta-data/spot/instance-action`
  (до этого — 404), в EventBridge — событие «EC2 Spot Instance Interruption Warning»; раньше может прийти
  рекомендация `events/recommendations/rebalance`.
- Исполнитель: машина `shutting-down`/`terminated` без `_DONE.json` → скачать, что успело лечь в S3, → один
  повтор (`--retries 1`) в тот же префикс, номер попытки — в `_DONE.json`; дальше итог `lost`. Файлы первой
  попытки не удаляются (у ключа нет `DeleteObject`) и лежат вперемешку со второй.
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
- Запись — только под префиксом `ALPHA_O_S3`, подпапка `gtu/<run_id>/`; регион — только us-east-1.
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
| `MaxSpotInstanceCountExceeded`, `VcpuLimitExceeded`, `InstanceLimitExceeded` | квота: задача ждёт, пока освободится машина прогона; своих машин нет — итог `quota`; поднять квоту — владелец |
| `InsufficientInstanceCapacity`, `Unsupported` | нет ёмкости или конфигурации в зоне: исполнитель берёт `c7i.16xlarge` |
| `SpotMaxPriceTooLow` | цена выше `MaxPrice`: тоже переход на запасной тип (код по документации не сверен) |
| фаза в `_DONE.json`: `base`/`venv`/`packages`/`clone` | сбой подготовки, коды 10–17; читать `_job.log` |
| `no_time` (код 124), `command` с кодом 124 или 137 | не хватило `max_min` (137 — `timeout` добил KILL): поднять срок или делить задачу |
| user-data больше 16 КБ | длинную команду — в скрипт репозитория, в задаче — его вызов |
| в журнале прогона нет строк | вывод в файл буферизован: запускать `python3 -u` |
| прогон оборвался через 30 мин или 2 ч | остановлена фоновая команда инструмента Bash по `timeout`: запускать через `setsid` (шаг 5) |
| `status` начал счёт с нуля | `ledger.jsonl` не закоммичен в прошлый раз: `run --spent-before` с расходом из консоли |
| пакет apt устарел (OpenFOAM 1912 на Ubuntu 24.04) | ставить нужную версию командой задачи, см. `openfoam-combustor` |

## 10. Проверка результата

- Итог `run`: у каждой задачи `status: ok`; в `_DONE.json` `phase: done`, `exit_code: 0`.
- `smoke`: `nproc.txt` = 64, `commit.txt` совпадает с последним коммитом запушенной ветки, в `_job.log` время
  подготовки.
- После прогона: `runner.py status` — открытых машин 0; `runner.py reap --run-id <id> --list-only` — живых нет.
- Повторяемость: `explore-249` сравнить с `routes/gtu/gtu_explore_results.json`; расхождение — разбор, не публикация.

## 11. Ссылки трека

`routes/gtu/platform/aws/` (`runner.py`, `test_runner.py`, `README.md`, `jobs.example.json`,
`iam_policy_addendum.json`); `calc/aws_orchestrate.py`; `calc/aws_setup/` (`setup.sh`, `README.md`, политики);
`routes/gtu/gtu_cycle.py`, `gtu_explore.py`, `gtu_dropin.py` и их `*_results.json`; навыки `gtu-cycle`,
`gtu-dropin`, `su2-turbo`, `openfoam-combustor`, `calculix-rotor`, `surrogate-opt`.

## Источники (проверено 14.11.2026)

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
- 64 vCPU у c7a/c7i.16xlarge — выдержки поиска по сторонним каталогам
  (https://instances.vantage.sh/aws/ec2/c7a.16xlarge); страница типов AWS недоступна из контейнера

**Не проверено:** квоты 384 vCPU в аккаунте и код `L-34B43A08`; коды `SpotMaxPriceTooLow`,
`MaxSpotInstanceCountExceeded`, `VcpuLimitExceeded` (в архиве нет; `VcpuLimitExceeded` упомянут в
`aws_orchestrate.py` как случившийся) и код причины гашения при прерывании; ресурс `spot-instances-request/*` —
по архиву 2023 г. проверяется, только если запрос спота помечается тегом при создании (`runner.py` его не
помечает), по текущей документации не сверено; есть ли в аккаунте `AWSServiceRoleForEC2Spot` (по архиву для
запуска спота через API она должна существовать); ставки gp3 и IPv4; успевает ли выгрузка при прерывании;
установка пакетов и `venv` на реальных образах; `curl` в образе Ubuntu; не даёт ли фоновый процесс уснуть ВМ
сессии; пуш из дочерней сессии.
