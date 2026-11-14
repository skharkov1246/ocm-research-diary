# Прогоны платформы ГТУ на AWS — `runner.py`

Узел платформы «Простая литая ГТУ-заменитель»: пакет независимых задач (SU2, OpenFOAM,
CalculiX, Gmsh и любые скрипты репозитория) идёт на спотовых машинах в us-east-1.
Основа — `calc/aws_orchestrate.py` и `calc/aws_setup/`; отличия: спот, много задач
параллельно, бюджетный предохранитель, журнал расходов, сборщик только своих машин.

Состояние на 14.11.2026: код и тесты на заглушках готовы. Боевых запусков не было:
в контейнере, где писался узел, ключи AWS не действуют, вызовов AWS не делалось.

## Файлы

| Файл | Что это |
|---|---|
| `runner.py` | Исполнитель: `plan`, `run` (по умолчанию DryRun; живой запуск — `--live`), `status`, `reap`, `collect` |
| `test_runner.py` | 25 тестов без сети: botocore Stubber вместо AWS; user-data выполняется локально на подменённых командах |
| `jobs.example.json` | Пример пакета: `smoke` (проверка машины) и `explore-249` (`routes/gtu/gtu_explore.py`) |
| `iam_policy_addendum.json` | Дополнение к `calc/aws_setup/user_policy.json`: спот, цена спота, образ Ubuntu 24.04 |
| `ledger.jsonl` | Журнал расходов (появится после первого боевого прогона; коммитить) |
| `runs/<run_id>/` | Скачанные результаты, `summary.json`, журнал DRYRUN (в `.gitignore`) |

## Как идёт прогон

1. Задачи читаются из JSON и проверяются (имена, пакеты, срок, отсутствие ключей в команде).
2. Перед каждым запуском: цена спота (`describe_spot_price_history`, окно 6 ч, все зоны) →
   проверка бюджета → `run_instances`. Одновременно не больше `--max-parallel` машин (6) и
   не больше `--vcpu-quota` vCPU (384 = 6 × 64).
3. Тип `c7a.16xlarge`; при `InsufficientInstanceCapacity`, `SpotMaxPriceTooLow`, `Unsupported` —
   запасной `c7i.16xlarge` (со своей проверкой бюджета). При квотных ошибках
   (`MaxSpotInstanceCountExceeded`, `VcpuLimitExceeded`) задача ждёт, пока освободится машина.
4. Параметры машины: разовый спот (`SpotInstanceType=one-time`, прерывание = terminate,
   `MaxPrice` задан), `InstanceInitiatedShutdownBehavior=terminate`, IMDSv2 обязателен
   (`HttpTokens=required`, hop limit 1), роль экземпляра из `ALPHA_O_INSTANCE_PROFILE`,
   без ключа SSH и входящих правил, gp3 с удалением при гашении, теги
   `Project=<ALPHA_O_TAG>`, `RunId=<run_id>`, `Job=<имя>`, `Name=gtu-<имя>` на машине и диске.
5. На машине (user-data): первой командой `shutdown -P +<max_min>` (сторож); затем git,
   python3-venv, boto3; пакеты задачи (`apt`/`dnf`, `pip`); `git clone --depth 1 -b <ветка>`
   публичного репозитория; команда задачи из корня клона с `timeout` до остатка срока минус
   8 мин на выгрузку. Результаты команда пишет в каталог `$GTU_OUT`.
6. Выгрузка — под ролью машины в `s3://<бакет>/<префикс>/gtu/<run_id>/<job>/`: всё из `$GTU_OUT`,
   `_job.log`, `_stdout.log` и последним — метка `_DONE.json` (фаза, код выхода, время). Затем `poweroff`
   (= terminate). Выгрузка и метка делаются и при сбое (ловушка EXIT), если успели поставиться boto3.
7. Исполнитель раз в `--poll-s` (60 с) читает метку. Метка есть — скачивает результаты в
   `runs/<run_id>/<job>/` и гасит машину. Машина погашена без метки (прерывание спота) — один повтор
   (`--retries 1`). Истёк срок `max_min + grace_min` (10 мин) — гасит сам.
8. В `finally` (конец, ошибка, Ctrl-C, SIGTERM, SIGHUP) — внешний сборщик: ищет машины с тегами
   `Project=<ALPHA_O_TAG>` **и** `RunId=<run_id>`, сверяет теги ещё раз на своей стороне и гасит только их
   и известные ему id. Машины соседних прогонов с тем же `Project` не трогает (урок
   `calc/aws_orchestrate.py`: широкий сборщик однажды убил чужой прогон).

## Формат задач

```json
{
  "defaults": {"instance_type": "c7a.16xlarge", "fallback_types": ["c7i.16xlarge"],
               "image": "ubuntu2404", "disk_gb": 50, "max_min": 30},
  "jobs": [
    {"name": "su2-comp130", "command": "python routes/gtu/platform/su2/run.py --out \"$GTU_OUT\"",
     "apt": ["gmsh"], "pip": ["numpy"], "max_min": 120}
  ]
}
```

Поля задачи: `name` (обязательно, `[A-Za-z0-9._-]`), `command` (обязательно, bash; несколько строк
допустимы), `instance_type`, `fallback_types`, `max_min` (20…720, срок сторожа), `image`
(`ubuntu2404` или `al2023`), `disk_gb`, `apt` или `dnf` (по образу), `pip`, `vcpu` (для типов вне
таблицы `VCPU` в коде). Неизвестное поле — ошибка (ловит опечатки). В команде доступны `GTU_OUT`,
`GTU_RUN_ID`, `GTU_JOB`, `GTU_ATTEMPT`, `OMP_NUM_THREADS=$(nproc)`; `python` — из venv машины.

## Бюджетный предохранитель

Перед каждым запуском:

```
итого = израсходовано + в работе + пакет (только DRYRUN/plan) + новая  ≤  потолок (--ceiling, 300)
новая        = (MaxPrice + диск) × (max_min + grace_min) / 60
MaxPrice     = максимум цены спота за 6 ч по всем зонам × 1,2 (--price-margin), вверх до 0,0001
израсходовано = Σ по закрытым машинам (цена спота в её зоне при запуске + диск) × max(60 с, конец − запуск)
               + --spent-before (расход вне журнала)
в работе     = Σ по открытым машинам (MaxPrice + диск) × весь отпущенный срок
диск         = disk_gb × 0,08 / 730 долл./ч (--ebs-usd-gb-month; ставка gp3 — допущение, в сессии не проверялась)
```

Спотовая машина не платит больше своего `MaxPrice`, сторож и сборщик гасят её не позже отпущенного
срока, поэтому «в работе» и «новая» — верхние границы. Отказ пишется в журнал (`ev=budget`, `ok=false`),
машина не запускается. В начале боевого прогона открытые записи прошлых прогонов сверяются с AWS:
погашенные закрываются по верхней границе срока.

Журнал `ledger.jsonl` живёт в контейнере, а контейнер эфемерен: **после каждого боевого прогона
коммитить `ledger.jsonl`**, иначе следующая сессия начнёт счёт с нуля. Внешняя страховка — бюджет
AWS Budgets (команда в конце `calc/aws_setup/setup.sh`). В журнале и выводе нет ключей и имени бакета
(в имени бакета номер аккаунта): бакет печатается как `<BUCKET>`.

## Где действуют ключи

Ключи ограниченного IAM-пользователя заданы переменными окружения Claude Code в окружении **Coldun**:
`AWS_ACCESS_KEY_ID`, `AWS_SECRET_ACCESS_KEY`, `AWS_DEFAULT_REGION`, `ALPHA_O_S3=s3://<бакет>/<префикс>`,
`ALPHA_O_INSTANCE_PROFILE`, `ALPHA_O_TAG`. Запускать `run`, `reap`, `collect` — из сессии, открытой в этом
окружении. В других контейнерах (как тот, где писался узел) ключи не действуют. Ключи не печатать и не
вставлять в чат, файлы, команды задач: исполнитель отказывается, если видит ключ в команде, user-data или
записи журнала. На машинах ключей нет — S3 через роль экземпляра по IMDSv2.

## Как запускать из сессии

Из корня репозитория. Ветка должна быть запушена в публичный репозиторий — машина клонирует её
(`run` проверяет это через `git ls-remote`).

```bash
pip install boto3 pytest                          # если нет
python3 -m pytest -q routes/gtu/platform/aws/test_runner.py
python3 routes/gtu/platform/aws/runner.py plan routes/gtu/platform/aws/jobs.example.json \
        --assume-price 1.5 --show-userdata        # без AWS; цена — допущение
python3 routes/gtu/platform/aws/runner.py status  # израсходовано / в работе / остаток
python3 routes/gtu/platform/aws/runner.py run routes/gtu/platform/aws/jobs.example.json \
        --branch <ветка>                          # по умолчанию DryRun: права и бюджет, без запуска
python3 routes/gtu/platform/aws/runner.py run routes/gtu/platform/aws/jobs.example.json \
        --branch <ветка> --live                   # живой прогон (или GTU_LIVE=1); расходы
```

Без `--live` (или `GTU_LIVE=1`) `run` машин не создаёт: `run_instances` идёт с `DryRun=True`.
`--dryrun`, `GTU_DRYRUN=1`, `ALPHA_O_DRYRUN=1` включают проверку даже при `--live`.

Долгий прогон — в фоне с журналом в файл: `... run ... --live > /tmp/gtu_run.log 2>&1 &`. Если сессия
оборвётся, машины погасит сторож (`max_min`), а `reap` по `run_id` (в первой строке вывода и в
`runs/<run_id>/`) погасит остатки и закроет записи.

## Порядок проверки (до первого боевого прогона)

1. Тесты: `python3 -m pytest -q routes/gtu/platform/aws/test_runner.py` — 25 passed.
2. `plan` — прочитать user-data и параметры `run_instances`.
3. Владелец с админ-доступом (не агент) применяет дополнение политики и создаёт
   служебную роль спота (если её нет):
   ```bash
   sed -e s/REGION/us-east-1/g -e s/ACCOUNT_ID/<номер аккаунта>/g \
       routes/gtu/platform/aws/iam_policy_addendum.json | grep -v '"_comment"' > /tmp/gtu_add.json
   aws iam put-user-policy --user-name alpha-o-agent --policy-name gtu-spot \
       --policy-document file:///tmp/gtu_add.json
   aws iam create-service-linked-role --aws-service-name spot.amazonaws.com   # ошибка «уже есть» — норма
   ```
   Не проверено в этой сессии (документация AWS недоступна из контейнера): что для спота через
   `RunInstances` нужен ресурс `spot-instances-request/*` и служебная роль `AWSServiceRoleForEC2Spot`.
   Путь SSM образа Ubuntu 24.04 сверен с документацией Canonical (шаблон
   `ubuntu/server/<выпуск>/stable/current/<арх>/hvm/ebs-gp3/ami-id`, gp3 — для 23.10 и новее).
4. `run` без `--live` в окружении Coldun: ожидается `запуск разрешён (DryRunOperation)` для каждой задачи
   и `чтение префикса результатов: OK`. `UnauthorizedOperation` — не хватает прав (п. 3).
5. Первый живой прогон (`--live`) — только `smoke` (вынести в отдельный файл задач, `max_min` 20): проверить
   `runs/<run_id>/smoke/` (`nproc.txt` = 64, `commit.txt`, `_job.log`, `_DONE.json`), `status`
   (открытых машин 0) и в консоли EC2 — что машина в `terminated`.
6. Затем `explore-249`: сравнить скачанный `gtu_explore_results.json` с `routes/gtu/gtu_explore_results.json`.
7. Закоммитить `ledger.jsonl`.

Не проверено на машине: установка пакетов, `python3 -m venv` и выгрузка на реальных образах;
профиль `al2023` (нужен ли отдельный пакет для venv); время подготовки машины (`smoke` даёт его в
`_job.log`).

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
- Результаты в S3 удаляет владелец (у ключа агента нет `DeleteObject`): префикс `<префикс>/gtu/<run_id>/`.
- Права спота отозвать: `aws iam delete-user-policy --user-name alpha-o-agent --policy-name gtu-spot`.
- Ключ после работ — удалить (`calc/aws_setup/README.md`).
- Локально: `runs/` можно удалять; `ledger.jsonl` не удалять и не править вручную — он держит учёт.
