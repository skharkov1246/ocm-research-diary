# Что сделать локально до заказа машин

Узел «первый тест в AWS» трека «Простая литая ГТУ-заменитель» (issue #59), пункт 4 задания владельца
«Расчёт промышленной газовой турбины на AWS» (14.11.2026). Всё ниже идёт на своей машине или в контейнере
сессии, без ключей AWS и без расходов. Правило трека: местный расчёт — не дольше 20 минут, длиннее — на AWS.

Состояние на 14.11.2026: шаги 1, 2, 7, 8, 9 выполнены в контейнере сессии (итоги — в README узла, раздел
«Что проверено»); шаги 3–6 — порядок действий, в этой сессии не выполнялись (SU2 узла su2 занят соседним
потоком, машина нагружена).

## 1. Собрать образ

Из корня чистой копии коммита (в образ не должны попасть чужие незакоммиченные файлы):

```bash
git worktree add --detach /tmp/gtu-img HEAD && cd /tmp/gtu-img
docker build -f routes/gtu/platform/docker/Dockerfile \
    --build-arg GIT_SHA="$(git rev-parse HEAD)" -t gtu-platform:"$(git rev-parse --short=12 HEAD)" .
```

- Docker Hub отвечает 429 на анонимные скачивания с общих адресов — тогда `--build-arg BASE_IMAGE=<зеркало>/ubuntu:24.04`
  (в сессии сработало `mirror.gcr.io/library/ubuntu:24.04`; ECR Public отдал манифест, но слои с CloudFront
  прокси сессии не пропустил).
- За TLS-прокси: `--network=host`, `--build-arg HTTPS_PROXY --build-arg NO_PROXY` и
  `--secret id=extra_ca,src=/etc/ssl/certs/ca-certificates.crt` — набор сертификатов доступен только командам
  сборки и в слои образа не попадает.
- Время: сборка SU2 — основная часть (835 единиц компиляции; на 4 нагруженных ядрах — около 15 мин).

## 2. Самопроверка образа

```bash
IMG=gtu-platform:<sha12>
docker run --rm $IMG                                    # справка entrypoint.py
docker run --rm --entrypoint bash $IMG -c 'cat /opt/su2/VERSION /opt/gtu/GIT_SHA; ccx -v || true;
  python -c "import cantera, numpy, scipy, gmsh, boto3; print(cantera.__version__, gmsh.__version__)";
  mpirun -np 2 SU2_CFD --help | head -3'
```

Ожидается: `v8.6.0 7048fd35…`, `Version 2.21` (ccx при этом выходит с кодом 201 — так у ccx 2.21), версии пакетов
из `docker/requirements.txt`, заставка SU2 «Release 8.6.0 "Harrier"» без ошибок OpenMPI.

## 3. Сетка одного межлопаточного канала

```bash
cd routes/gtu/platform/geom
python impeller.py --case taurus60 --level coarse    # 165 236 узлов, 38 с на 4 ядрах (README узла geom)
python impeller.py --case taurus60 --level medium    # уровень серии; core130 medium — 584 357 узлов, 73 с
```

Проверить: семь маркеров, периодичность узел в узел (README узла geom, проверки 7–8). Число узлов medium
записать в `cost_estimate.toml` → `su2_time.mesh_nodes` (сейчас EXAMPLE 600 000).

## 4. Грубый прогон SU2 на 4 ядрах

```bash
cd routes/gtu/platform/su2
python run_case.py --case taurus60 --level coarse --plan smoke --np 4 --time-limit 1200 --out runs/local_smoke
```

`smoke` — 520 итераций (разгон 0 → 0,6 → 1,0 номинала и 150 итераций второго порядка). Это проверка постановки,
не число для публикации: на этой сетке номинальная частота пока не достигнута (навык `su2-turbo`, раздел 4).

## 5. Проверка сходимости

По `run_summary.json` и `history_*.csv` прогона (критерии навыка `su2-turbo`, раздел 6):

- невязки упали не меньше чем на 4 порядка, без NaN;
- дисбаланс расхода вход/выход ≤ 0,5 %, знаки: вход > 0, выход < 0;
- π и η за последние 500 итераций меняются меньше чем на 0,1 %;
- в журнале «Matched N points» без «Bad match».

Не сошлось локально — на AWS не отправлять: машина не исправит постановку.

## 6. Замер времени итерации

Столбец `Time(sec)` (ITER_TIME) в `history_*.csv` — время одной итерации. Мерить на ненагруженной машине,
с числом процессов, равным числу ядер. Записать в `cost_estimate.toml`:

```toml
[su2_time]
ref_s_per_iter = <среднее ITER_TIME>   # сейчас 1,2255 с — дымовой прогон узла su2 на нагруженной машине
ref_nodes = <узлов сетки>              # сейчас 139 121 (core130 coarse)
ref_procs = <процессов>                # сейчас 4 (EXAMPLE)
```

и пересчитать смету: `python3 routes/gtu/platform/batch/cost_estimate.py` (выход 1 — потолок превышен).
После облачного замера (`bench_*`) туда же — `parallel_eff` и `cpu_factor` по типам.

## 7. Точка входа без AWS

Хранилище S3 подменяется папкой: `s3://bucket/key` ↔ `<папка>/bucket/key`.

```bash
cd routes/gtu/platform/batch
python3 render.py --out /tmp/gtu_render                      # план, job.json, manifest.json, запросы Batch
mkdir -p /tmp/gtu_s3/PLACEHOLDER-bucket/PLACEHOLDER-prefix/gtu/gtu-EXAMPLE-kv3-char
cp -r /tmp/gtu_render/s3/* /tmp/gtu_s3/PLACEHOLDER-bucket/PLACEHOLDER-prefix/gtu/gtu-EXAMPLE-kv3-char/
python3 entrypoint.py --job s3://PLACEHOLDER-bucket/PLACEHOLDER-prefix/gtu/gtu-EXAMPLE-kv3-char/inputs/job.json \
    --index 3 --local-root /tmp/gtu_s3 --dry-run             # строки линии 3 и команды, без запуска
```

То же в контейнере (команда своя, чтобы не запускать SU2):
`docker run --rm -v /tmp/gtu_s3:/s3 $IMG --job s3://…/inputs/job.json --index 3 --local-root /s3 --dry-run`.

## 8. Шаблоны Batch и смета

```bash
python3 render.py                       # us-east-1; сверка с моделью API batch (botocore) и ограничениями
python3 render.py --region eu-north-1   # запасной регион
python3 cost_estimate.py                # cost_estimate_results.json, cost_estimate_lines.csv; проверка потолка
python3 -m pytest -q test_batch.py
```

## 9. План раннера для первого замера

```bash
python3 routes/gtu/platform/aws/runner.py plan routes/gtu/platform/batch/bench_runner_jobs.json \
        --assume-price 1.5 --show-userdata --out /tmp/gtu_runs   # цена — допущение, не число для отчёта
```

## Чем первый прогон отличается: раннер и Batch

Первый замер (один канал SU2 на c7a.16xlarge и c8i.8xlarge) идёт через `runner.py` — EC2 спот напрямую, он уже
есть и проверен на заглушках. Batch — для серий (характеристика: 80 точек), после замера.

| | `runner.py` (первый прогон) | AWS Batch (серии) |
|---|---|---|
| что запускается | машина Ubuntu 24.04 (SSM), пакеты apt/pip и сборка SU2 на машине, `git clone` ветки | образ `gtu-platform` из ECR на машинах ECS AL2023 |
| выбор машины | точный тип задачи (`instance_type`), запасные — `fallback_types` (в замере пусто: сравнение чистое) | список типов среды, выбирает Batch (`SPOT_PRICE_CAPACITY_OPTIMIZED`) |
| рынок | только спот (`MarketType=spot` в коде) | спот (`type: SPOT`); по требованию — отдельная среда `type: EC2` |
| c8i.8xlarge | нужно поле `"vcpu": 32` (типа нет в таблице `VCPU` раннера) | в серии не используется (только 16xlarge) |
| параллельно | `--max-parallel 6`, `--vcpu-quota 384` | `maxvCpus` 320 (+ одна машина сверх = 384) |
| срок | `max_min` (сторож `shutdown -P` на машине) + `grace_min` | `attemptDurationSeconds` 12 600 с; точка входа сама останавливает расчёт за 600 с до него |
| повтор | `--retries 1` при потере машины без метки | `retryStrategy`: 2 попытки, повтор только при `Host EC2*`; тайм-аут не повторяется |
| бюджет | предохранитель `--ceiling`, журнал `ledger.jsonl`, живая цена спота перед каждым запуском | встроенного нет: смета `cost_estimate.py` до отправки, `maxvCpus`, отмена RUNNABLE через 2 ч, AWS Budgets |
| проверка без затрат | `run` без `--live` = DryRun EC2 (права и бюджет) | у API Batch нет DryRun (в модели API нет такого параметра): только `render.py` |
| журналы | `_job.log`, `_stdout.log` в S3 | CloudWatch Logs (`awslogs`) + `logs/<row_id>/attempt-<n>.log` в S3 |
| S3 | `<префикс>/gtu/<run_id>/<job>/` | `<префикс>/gtu/<run_id>/{inputs,outputs,logs,manifest.json}` |
| что создаёт владелец | дополнение политики, служебная роль спота | репозиторий ECR, роль машин ECS, роль задания, роль Spot Fleet, среда, очередь, определение (`commands.sh`) |
