#!/usr/bin/env python3
"""Прогоны расчётной платформы ГТУ на спотовых машинах AWS (us-east-1 и соседние регионы), версия 3.

Запускается ИЗ СЕССИИ, где действуют ключи ограниченного IAM-пользователя (окружение
Coldun). Машины работают под ролью экземпляра: ключей на машине и в user-data нет.
Подробности, порядок проверки и отката — README.md рядом.

  python3 runner.py plan    jobs.json --assume-price 1.5   # без AWS: user-data, параметры, бюджет
  python3 runner.py run     jobs.json --branch B            # по умолчанию: проверка прав (DryRun=True), без запуска
  python3 runner.py run     jobs.json --branch B --live     # живой запуск спотовых машин (расходы!)
  python3 runner.py status                                  # журнал: израсходовано / в работе / остаток
  python3 runner.py reap    --run-id ID                     # погасить машины прогона по тегам
  python3 runner.py collect --run-id ID --job NAME          # скачать результаты задачи

Защита от утечки машин — тройная, как в calc/aws_orchestrate.py:
  1) InstanceInitiatedShutdownBehavior=terminate и разовый спот (прерывание = terminate);
  2) сторож на машине: `shutdown -P +N` — первая команда user-data;
  3) внешний сборщик: в finally гасит машины с тегами Project=<ALPHA_O_TAG> и
     RunId=<run_id> — только этого прогона, соседние прогоны с тем же Project не трогает.
Бюджет: перед каждым запуском «израсходовано + в работе (верхняя граница) + новая задача»
против потолка. Цена — describe_spot_price_history; MaxPrice спота ставится равной оценке
по максимуму цены за окно по всем зонам, поэтому почасовая цена машины сверху ограничена
числом, вошедшим в прогноз.

Версия 2 (после первого живого прогона 14.11.2026):
  А) самая дешёвая зона: зоны типа машины — по возрастанию последней цены спота, подсеть по
     умолчанию зоны (describe_subnets, default-for-az); при InsufficientInstanceCapacity /
     SpotMaxPriceTooLow / Unsupported — следующая зона, затем запасной тип; бюджет — по максимуму;
  Б) образ по умолчанию al2023; поле "python": "3.12" — venv на этой версии через uv;
  В) поле "container" (образ Docker) и "container_setup" — команда задачи внутри контейнера,
     выгрузка, метка и гашение — снаружи;
  Г) очередь под квоту (--max-vcpu, 384): запуск по мере освобождения машин, бюджет — перед
     каждым запуском; повтор после прерывания спота — в свой префикс попытки (без задвоения);
  Д) на машине переменная NPROC (число ядер).
Правки после трёх проверок v2: машина с неподтверждённым гашением останавливает очередь (--max-vcpu);
сетевой сбой гашения не уходит из finally; машина за ошибкой 5xx ищется по ClientToken; сверка гасит машину,
живую после срока; DryRun-успех гасится; user-data всех попыток — до первого запуска; статусы capacity
(с повтором), no_price, halted; метка 143 при гашении изнутри (Client.…) — не прерывание; окружение OpenFOAM
(/openfoam/bash.rc) в контейнере.

Версия 3 — соседние мощности (решение владельца 15.11.2026). Без новых флагов — точно как версия 2: те же
вызовы AWS, тот же журнал. Новые флаги:
  Е) --regions us-east-1,us-east-2,us-west-2: на регион — свои клиенты EC2/SSM, образ из SSM своего региона,
     подсеть по умолчанию самой дешёвой зоны, группа безопасности по умолчанию и своя квота
     (--vcpu-quota-per-region, иначе service-quotas, иначе запасное значение); профиль роли (IAM глобален)
     и бакет (в --region, us-east-1) — общие; регион очередной задачи — с наименьшей ценой спота среди
     регионов со свободной квотой (соседний — с запасом цены --region-price-margin); бюджет и MaxPrice — по
     максимуму цены по всем разрешённым регионам; межрегиональная передача результатов — в смете;
  Ж) --on-demand-fallback: спот недоступен во всех регионах (нет ёмкости, цены, прав или квоты спота) —
     машина по требованию по цене из справочника aws_prices_ref.json, при своей квоте по требованию;
  З) status, reap, collect и сверка — по всем регионам, где у прогона были машины; журнал пишет регион и рынок.
"""
from __future__ import annotations

import argparse
import datetime as dt
import hashlib
import json
import math
import os
import re
import secrets as _secrets
import shlex
import signal
import subprocess
import sys
import time
from collections import deque
from dataclasses import dataclass, field
from pathlib import Path

HERE = Path(__file__).resolve().parent
REGION = "us-east-1"
DEFAULT_TYPE = "c7a.16xlarge"
DEFAULT_FALLBACK = ("c7i.16xlarge",)
DEFAULT_REPO = "https://github.com/skharkov1246/ocm-research-diary"
DEFAULT_TAG = "sandbox-agent"          # как в calc/aws_setup/setup.sh
DEFAULT_IMAGE = "al2023"               # ubuntu2404: параметр SSM Canonical агенту не разрешён (прогон 14.11.2026)
PRODUCT = "Linux/UNIX"
USERDATA_LIMIT = 16384                 # байт до base64 (ограничение EC2)
LIVE_STATES = ["pending", "running", "stopping", "stopped"]
SPOT_REASON = "Server.SpotInstanceTermination"   # код причины гашения при прерывании (по архиву EC2, не сверен)
# Прерывание спота — причина Server.Spot… (Termination или Shutdown). Причина Client.… — машину погасили
# изнутри или мы сами (сторож shutdown -P, poweroff после итога, наш terminate): это не прерывание,
# даже если ловушка TERM успела выгрузить метку с кодом 143.
SPOT_REASON_PREFIX = "Server.Spot"
OWN_REASON_PREFIX = "Client."

# vCPU известных типов (квота спота считается в vCPU). Другие типы — поле "vcpu" в задаче.
VCPU = {f"c7{f}.{s}": n for f in "ai"
        for s, n in (("4xlarge", 16), ("8xlarge", 32), ("16xlarge", 64))}

IMAGES = {
    "al2023": {"ssm": "/aws/service/ami-amazon-linux-latest/al2023-ami-kernel-default-x86_64",
               "root": "/dev/xvda", "pkg": "dnf"},
    # Параметр SSM Canonical в живом прогоне 14.11.2026 — AccessDenied (дополнение политики не применено).
    # Запасной путь — DescribeImages по владельцу Canonical (099720109477); вживую не проверен: есть ли
    # у ключа агента ec2:DescribeImages — неизвестно. Не вышло — задача получает статус image_unavailable.
    "ubuntu2404": {"ssm": "/aws/service/canonical/ubuntu/server/24.04/stable/current/"
                          "amd64/hvm/ebs-gp3/ami-id",
                   "root": "/dev/sda1", "pkg": "apt", "owner": "099720109477",
                   "name": "ubuntu/images/hvm-ssd-gp3/ubuntu-noble-24.04-amd64-server-*"},
}

CAPACITY_ERRORS = {"InsufficientInstanceCapacity", "SpotMaxPriceTooLow", "Unsupported"}
QUOTA_ERRORS = {"MaxSpotInstanceCountExceeded", "VcpuLimitExceeded", "InstanceLimitExceeded"}
# Версия 3: отказ в правах в регионе (политика ключа агента закреплена за us-east-1 условием
# aws:RequestedRegion; регион не включён в аккаунте) — регион исключается до конца прогона для этого рынка.
DENY_ERRORS = {"UnauthorizedOperation", "AccessDenied", "AccessDeniedException", "AuthFailure", "OptInRequired",
               "UnrecognizedClientException", "InvalidClientTokenId"}

# --------------------------------------------------------------------------- версия 3: регионы
REGION_RE = re.compile(r"^[a-z]{2}(?:-[a-z]+)+-\d$")
MARKETS = ("spot", "on-demand")
# Коды квот EC2 в vCPU: L-1216C47A — Running On-Demand Standard (A, C, D, H, I, M, R, T, Z) instances (сверено по
# примеру awscli 1.46.3 examples/service-quotas/get-service-quota.rst); L-34B43A08 — All Standard (A, C, D, H, I,
# M, R, T, Z) Spot Instance Requests (по выдержкам поиска, не сверено).
QUOTA_CODES = {"spot": "L-34B43A08", "on-demand": "L-1216C47A"}
QUOTA_FALLBACK = 64                    # квота соседнего региона, если её не задали и service-quotas не ответил
PRICE_REF = HERE / "aws_prices_ref.json"   # справочник: цены по требованию и межрегиональная передача, с датой

NAME_RE = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._-]{0,63}$")
BRANCH_RE = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._/-]{0,199}$")
REPO_RE = re.compile(r"^https://github\.com/[A-Za-z0-9_.-]+/[A-Za-z0-9_.-]+$")
PKG_RE = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._+\-\[\],:=<>!~]{0,127}$")
RUN_ID_RE = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._-]{0,47}$")
AZ_RE = re.compile(r"^[a-z]{2}(?:-[a-z]+)+-\d[a-z]$")
PY_RE = re.compile(r"^3\.\d{1,2}(?:\.\d{1,2})?$")
REQ_RE = re.compile(r"^[A-Za-z0-9_][A-Za-z0-9._/-]{0,199}$")
IMAGE_REF_RE = re.compile(
    r"^(?:[a-z0-9][a-z0-9.-]*(?::\d{1,5})?/)?[a-z0-9]+(?:[._-][a-z0-9]+)*(?:/[a-z0-9]+(?:[._-][a-z0-9]+)*)*"
    r"(?::[A-Za-z0-9_][A-Za-z0-9_.-]{0,127})?(?:@sha256:[0-9a-f]{64})?$")
HEREDOC_MARKS = ("GTU_CMD_EOF", "GTU_SETUP_EOF", "GTU_ENTRY_EOF", "GTU_PY_EOF")
# apt-get / apt install|update в команде задачи без контейнера на образе al2023 (dnf) — почти наверняка
# задача написана под прежний образ Ubuntu: предупреждение в plan и run (не отказ — совместимость).
APT_IN_CMD_RE = re.compile(r"(?:^|[\s;&|(])(?:sudo\s+)?apt(?:-get)?\s+(?:-\S+\s+)*(?:install|update)\b", re.M)
KEY_ID_RE = re.compile(r"(?<![A-Z0-9])(AKIA|ASIA)[A-Z0-9]{16}(?![A-Z0-9])")
ARN_ACCOUNT_RE = re.compile(r"(arn:aws[a-z-]*:[a-z0-9-]*:[a-z0-9-]*:)\d{12}(?=:)")   # номер аккаунта в ARN
SECRET_WORDS = ("aws_secret_access_key", "aws_session_token", "aws_access_key_id")
SECRET_ENV = ("AWS_ACCESS_KEY_ID", "AWS_SECRET_ACCESS_KEY", "AWS_SESSION_TOKEN")


class RunnerError(Exception):
    """Ошибка конфигурации или отказ (без вывода секретов)."""


class ImageUnavailable(RunnerError):
    """Образ машины не найден или не разрешён ключу агента."""


# --------------------------------------------------------------------------- секреты
def find_secret(text: str) -> str | None:
    """Вернуть описание найденного секрета (без самого значения) или None."""
    for name in SECRET_ENV:
        v = os.environ.get(name)
        if v and len(v) >= 8 and v in text:
            return f"значение переменной {name}"
    if KEY_ID_RE.search(text):
        return "строка вида идентификатора ключа AWS (AKIA…/ASIA…)"
    low = text.lower()
    for w in SECRET_WORDS:
        if w in low:
            return f"имя переменной ключа ({w})"
    return None


def assert_no_secret(text: str, what: str) -> None:
    hit = find_secret(text)
    if hit:
        raise RunnerError(f"{what}: найдено {hit} — отказ, ключи не должны попадать в {what}")


class Redactor:
    """Маскирует имя бакета (в нём номер аккаунта) и значения ключей в выводе и журнале."""

    def __init__(self, bucket: str | None = None):
        self.pairs = []
        for name in SECRET_ENV:
            v = os.environ.get(name)
            if v and len(v) >= 8:
                self.pairs.append((v, "<скрыто>"))
        if bucket:
            self.pairs.append((bucket, "<BUCKET>"))

    def __call__(self, text: str) -> str:
        for s, r in self.pairs:
            text = text.replace(s, r)
        text = ARN_ACCOUNT_RE.sub(r"\1<аккаунт>", text)    # в тексте ошибок AWS (AccessDenied и т. п.)
        return KEY_ID_RE.sub("<скрыто>", text)


# --------------------------------------------------------------------------- задачи
@dataclass
class Job:
    name: str
    command: str
    instance_type: str = DEFAULT_TYPE
    fallback_types: list = field(default_factory=lambda: list(DEFAULT_FALLBACK))
    max_min: int = 60
    image: str = ""                 # "" — al2023; у старых задач с полем apt — ubuntu2404 (как было)
    disk_gb: int = 50
    apt: list = field(default_factory=list)
    dnf: list = field(default_factory=list)
    pip: list = field(default_factory=list)
    vcpu: int = 0
    python: str = ""                # "3.12" — venv на этой версии через uv; "" — системный python3
    requirements: str = ""          # файл зависимостей в репозитории, ставится после клонирования
    container: str = ""             # образ Docker: команда задачи идёт внутри контейнера
    container_setup: list = field(default_factory=list)   # команды внутри контейнера до команды задачи

    @property
    def types(self) -> list:
        return [self.instance_type] + [t for t in self.fallback_types if t != self.instance_type]

    @property
    def mode(self) -> str:
        return "container" if self.container else ("python" if self.python else "system")


JOB_KEYS = {f for f in Job.__dataclass_fields__}


def load_jobs(path, max_job_min: int = 720) -> list:
    """Прочитать задачи: {"defaults": {...}, "jobs": [...]} или просто [...]."""
    try:
        data = json.loads(Path(path).read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as e:
        raise RunnerError(f"файл задач не прочитан: {type(e).__name__}") from None
    if isinstance(data, list):
        defaults, items = {}, data
    elif isinstance(data, dict):
        defaults, items = data.get("defaults", {}), data.get("jobs")
    else:
        raise RunnerError("файл задач: ожидается список задач или объект {defaults, jobs}")
    if not isinstance(defaults, dict):
        raise RunnerError("defaults — объект")
    if not isinstance(items, list) or not items:
        raise RunnerError("в файле задач нет непустого списка jobs")
    jobs, seen = [], set()
    for raw in items:
        if not isinstance(raw, dict):
            raise RunnerError(f"задача — объект, а не {type(raw).__name__}")
        merged = {**defaults, **raw}
        unknown = set(merged) - JOB_KEYS
        if unknown:
            raise RunnerError(f"неизвестные поля задачи: {sorted(unknown)}")
        if "name" not in merged or "command" not in merged:
            raise RunnerError("у задачи обязательны поля name и command")
        job = Job(**merged)
        validate_job(job, max_job_min)
        if job.name in seen:
            raise RunnerError(f"повтор имени задачи: {job.name}")
        seen.add(job.name)
        jobs.append(job)
    return jobs


def validate_job(job: Job, max_job_min: int = 720) -> None:
    if not isinstance(job.name, str) or not NAME_RE.match(job.name):
        raise RunnerError(f"имя задачи {job.name!r}: только [A-Za-z0-9._-], до 64 знаков")
    if not isinstance(job.command, str) or not job.command.strip():
        raise RunnerError(f"{job.name}: пустая команда")
    if any(line.strip() == "GTU_CMD_EOF" for line in job.command.splitlines()):
        raise RunnerError(f"{job.name}: команда содержит служебную строку GTU_CMD_EOF")
    assert_no_secret(job.command, f"команда задачи {job.name}")
    if not isinstance(job.fallback_types, list):
        raise RunnerError(f"{job.name}: fallback_types — список типов")
    for lst in (job.apt, job.dnf, job.pip, job.container_setup):
        if not isinstance(lst, list):
            raise RunnerError(f"{job.name}: apt, dnf, pip, container_setup — списки")
    for k in ("python", "requirements", "container"):
        if not isinstance(getattr(job, k), str):
            raise RunnerError(f"{job.name}: {k} — строка (версию Python — в кавычках: \"3.12\")")
    if job.image in ("", None):
        # Старые задачи с apt рассчитаны на Ubuntu (прежний образ по умолчанию) — для них образ прежний.
        job.image = "ubuntu2404" if job.apt else DEFAULT_IMAGE
    if isinstance(job.max_min, bool) or not isinstance(job.max_min, int) or not 20 <= job.max_min <= max_job_min:
        raise RunnerError(f"{job.name}: max_min должен быть целым в [20, {max_job_min}]")
    if job.image not in IMAGES:
        raise RunnerError(f"{job.name}: image — одно из {sorted(IMAGES)}")
    if not isinstance(job.disk_gb, int) or not 8 <= job.disk_gb <= 2000:
        raise RunnerError(f"{job.name}: disk_gb в [8, 2000]")
    pkg = IMAGES[job.image]["pkg"]
    if pkg == "apt" and job.dnf or pkg == "dnf" and job.apt:
        raise RunnerError(f"{job.name}: для image={job.image} пакеты задаются в поле {pkg}")
    for p in [*job.apt, *job.dnf, *job.pip]:
        if not isinstance(p, str) or not PKG_RE.match(p):
            raise RunnerError(f"{job.name}: недопустимое имя пакета {p!r}")
    if job.python and not PY_RE.match(job.python):
        raise RunnerError(f"{job.name}: python — версия вида \"3.12\" или \"3.12.7\"")
    if job.requirements and (not REQ_RE.match(job.requirements) or ".." in job.requirements.split("/")
                             or "//" in job.requirements):
        raise RunnerError(f"{job.name}: requirements — относительный путь файла в репозитории")
    for c in job.container_setup:
        if not isinstance(c, str) or not c.strip():
            raise RunnerError(f"{job.name}: container_setup — список непустых команд")
        if any(line.strip() in HEREDOC_MARKS for line in c.splitlines()):
            raise RunnerError(f"{job.name}: container_setup содержит служебную строку")
        assert_no_secret(c, f"container_setup задачи {job.name}")
    if job.container:
        if len(job.container) > 255 or not IMAGE_REF_RE.match(job.container):
            raise RunnerError(f"{job.name}: container — имя образа вида ubuntu:24.04 или org/image:tag")
        if job.image != "al2023":
            raise RunnerError(f"{job.name}: container идёт на образе машины al2023 (docker из dnf)")
        extra = [k for k in ("apt", "pip", "python", "requirements") if getattr(job, k)]
        if extra:
            raise RunnerError(f"{job.name}: с container поля {extra} не применяются — "
                              f"ставьте пакеты командами container_setup внутри контейнера")
    elif job.container_setup:
        raise RunnerError(f"{job.name}: container_setup без container")
    for t in job.types:
        if not isinstance(t, str) or not re.match(r"^[a-z0-9]+\.[a-z0-9]+$", t):
            raise RunnerError(f"{job.name}: недопустимый тип машины {t!r}")
    if isinstance(job.vcpu, bool) or not isinstance(job.vcpu, int) or not 0 <= job.vcpu <= 1024:
        raise RunnerError(f"{job.name}: vcpu — целое в [0, 1024] (0 — по таблице типов)")
    missing = [t for t in job.types if t not in VCPU]
    if missing and not job.vcpu:
        raise RunnerError(f"{job.name}: нет числа vCPU для {missing} — задайте поле vcpu")
    # Для типов из таблицы число vCPU берётся из неё: заниженное поле vcpu не обходит квоту.
    job.vcpu = max([job.vcpu] + [VCPU[t] for t in job.types if t in VCPU])


def job_warnings(job: Job) -> list:
    """Предупреждения по задаче (не отказ): то, что почти наверняка не сработает на машине."""
    out = []
    if not job.container and IMAGES[job.image]["pkg"] == "dnf" and APT_IN_CMD_RE.search(job.command):
        out.append(f"{job.name}: команда вызывает apt-get, а образ машины {job.image} (dnf, apt-get нет) — "
                   f"задайте \"container\": \"ubuntu:24.04\" и перенесите установку в container_setup")
    return out


# --------------------------------------------------------------------------- user-data
UPLOAD_PY = r'''import os, sys, boto3
bucket, prefix, src = sys.argv[1], sys.argv[2].rstrip("/"), sys.argv[3]
s3, bad = boto3.client("s3"), 0
pairs = ([(os.path.join(r, f), os.path.relpath(os.path.join(r, f), src)) for r, _, fs in os.walk(src) for f in sorted(fs)]
         if os.path.isdir(src) else [(src, os.path.basename(src))])
for path, rel in pairs:
    try:
        s3.upload_file(path, bucket, prefix + "/" + rel.replace(os.sep, "/"))
    except Exception as e:
        bad += 1; print("upload failed:", rel, type(e).__name__)
sys.exit(1 if bad else 0)
'''

# Сборка в два прохода: сначала блоки режима (@@@ИМЯ@@@, текст раннера), затем значения (@@ИМЯ@@)
# одним проходом — текст команды и container_setup повторно не сканируется.
USERDATA_TEMPLATE = r'''#!/bin/bash
shutdown -P +@@MAX_MIN@@
# gtu runner v2: run=@@RUN_ID@@ job=@@JOB@@ attempt=@@ATTEMPT@@ mode=@@MODE@@; the watchdog above is the first command
exec >@@LOG@@ 2>&1
set -u
export HOME="${HOME:-/root}" DEBIAN_FRONTEND=noninteractive AWS_DEFAULT_REGION=@@REGION@@
T0=$(date +%s); DEADLINE=$((T0 + @@BUDGET_S@@)); WD_END=$((T0 + @@MAX_MIN@@ * 60)); TP=$T0; T_CMD=0
G=@@ROOT@@; OUT=$G/out; PHASE=setup; AZ_M=""; IT_M=""
mkdir -p "$OUT"
NPROC=$(nproc)
export GTU_OUT="$OUT" GTU_RUN_ID=@@RUN_ID_Q@@ GTU_JOB=@@JOB_Q@@ GTU_ATTEMPT=@@ATTEMPT@@ NPROC
export GTU_AZ=@@AZ_Q@@ GTU_AZ_CHEAPEST=@@AZ_CHEAP_Q@@ GTU_ITYPE=@@ITYPE_Q@@
S3_BUCKET=@@BUCKET_Q@@; S3_KEY=@@KEY_Q@@
cat >"$G/upload.py" <<'GTU_PY_EOF'
@@UPLOAD_PY@@GTU_PY_EOF
phase() {  # log how long the previous phase took (machine log), then switch to $1
  local n; n=$(date +%s)
  echo "[gtu-time] $PHASE $((n - TP))s total=$((n - T0))s"
  TP=$n; PHASE=$1
}
imds() {  # IMDSv2 metadata path -> stdout, empty if unavailable
  local t
  t=$(curl -sf -m 2 -X PUT -H 'X-aws-ec2-metadata-token-ttl-seconds: 300' http://169.254.169.254/latest/api/token) || return 0
  curl -sf -m 2 -H "X-aws-ec2-metadata-token: $t" "http://169.254.169.254/latest/meta-data/$1" || true
}
up() {  # $1 - path, $2 - time limit in seconds
  if [ -x "$G/venv/bin/python" ] && "$G/venv/bin/python" -c 'import boto3' 2>/dev/null; then
    timeout "$2s" "$G/venv/bin/python" "$G/upload.py" "$S3_BUCKET" "$S3_KEY" "$1"
  elif command -v aws >/dev/null 2>&1; then
    if [ -d "$1" ]; then timeout "$2s" aws s3 cp --recursive --only-show-errors "$1" "s3://$S3_BUCKET/$S3_KEY/"
    else timeout "$2s" aws s3 cp --only-show-errors "$1" "s3://$S3_BUCKET/$S3_KEY/"; fi
  else echo "UPLOAD-UNAVAILABLE"; fi
}
finish() {
  RC=$?; set +e; trap - EXIT
@@@STOP@@@
  echo "[gtu-time] $PHASE $(( $(date +%s) - TP ))s total=$(( $(date +%s) - T0 ))s"
  echo "finish phase=$PHASE rc=$RC"
  [ -n "$AZ_M" ] || AZ_M=$(imds placement/availability-zone | tr -cd 'a-z0-9-')
  cp @@LOG@@ "$OUT/_job.log"
  # results must finish before the watchdog: keep 120 s for the marker and poweroff
  LEFT=$((WD_END - $(date +%s) - 120)); [ "$LEFT" -gt 1800 ] && LEFT=1800; [ "$LEFT" -lt 30 ] && LEFT=30
  up "$OUT" "$LEFT"
  printf '{"run_id":"%s","job":"%s","attempt":%s,"phase":"%s","exit_code":%d,"t_start":%d,"t_cmd":%d,"t_end":%d,"nproc":%d,"az":"%s","az_planned":"%s","itype":"%s"}\n' \
    "$GTU_RUN_ID" "$GTU_JOB" "$GTU_ATTEMPT" "$PHASE" "$RC" "$T0" "$T_CMD" "$(date +%s)" "$NPROC" "$AZ_M" "$GTU_AZ" "$IT_M" >"$G/_DONE.json"
  up "$G/_DONE.json" 60
  poweroff
}
trap finish EXIT
trap 'exit 143' TERM HUP INT
phase base
@@@BASE@@@
AZ_M=$(imds placement/availability-zone | tr -cd 'a-z0-9-'); IT_M=$(imds instance-type | tr -cd 'a-z0-9.')
export GTU_AZ_MACHINE="$AZ_M"
echo "[gtu] machine az=$AZ_M planned=$GTU_AZ cheapest=$GTU_AZ_CHEAPEST type=$IT_M nproc=$NPROC"
phase venv
@@@VENV@@@
phase packages
@@@PACKAGES@@@
@@@DOCKER@@@
phase clone
timeout 10m git clone -q --depth 1 -b @@BRANCH_Q@@ @@REPO_Q@@ "$G/repo" || exit 14
cd "$G/repo" || exit 15
echo "[gtu] commit $(git rev-parse HEAD 2>/dev/null)"
@@@REQUIREMENTS@@@
@@@PREP@@@
cat >@@@CMD_FILE@@@ <<'GTU_CMD_EOF'
@@COMMAND@@
GTU_CMD_EOF
export PATH="$G/venv/bin:$PATH" OMP_NUM_THREADS="$NPROC"
REMAIN=$((DEADLINE - $(date +%s)))
[ "$REMAIN" -gt 120 ] || { PHASE=no_time; exit 124; }
phase command
T_CMD=$(date +%s)
@@@RUN@@@
phase done
exit $RC
'''

VENV_SYSTEM = r'''python3 -m venv "$G/venv" || exit 12
timeout 20m "$G/venv/bin/pip" install -q --upgrade pip boto3 || exit 13'''

# uv: из PATH, иначе pip в отдельный venv системного python3, иначе архив выпуска GitHub.
VENV_UV = r'''UV=$(command -v uv 2>/dev/null || true)
if [ -z "$UV" ] && python3 -m venv "$G/uvboot" && timeout 10m "$G/uvboot/bin/pip" install -q uv; then
  UV="$G/uvboot/bin/uv"
fi
if [ ! -x "$UV" ]; then
  echo "[gtu] uv from PyPI failed - trying the GitHub release"
  mkdir -p "$G/uvrel"
  timeout 5m curl -fsSL "https://github.com/astral-sh/uv/releases/latest/download/uv-$(uname -m)-unknown-linux-gnu.tar.gz" | tar -xz -C "$G/uvrel" --strip-components=1
  UV="$G/uvrel/uv"
fi
[ -x "$UV" ] || exit 18
export UV_PYTHON_INSTALL_DIR="$G/python" UV_PYTHON_BIN_DIR="$G/python/bin" UV_CACHE_DIR="$G/uv-cache" UV_NO_PROGRESS=1
echo "[gtu] $("$UV" --version 2>&1)"
timeout 15m "$UV" python install @@PY@@ || exit 19
timeout 10m "$UV" venv --seed -q --python @@PY@@ "$G/venv" || exit 12
timeout 20m "$UV" pip install -q --python "$G/venv/bin/python" boto3 || exit 13
echo "[gtu] venv $("$G/venv/bin/python" -V 2>&1)"'''

DOCKER_BLOCK = r'''phase docker
timeout 5m systemctl start docker || exit 21
for _ in $(seq 90); do docker info >/dev/null 2>&1 && break; sleep 2; done
docker info >/dev/null 2>&1 || exit 21
phase pull
timeout 30m docker pull -q @@IMG_Q@@ || exit 22
echo "[gtu] image digests: $(docker image inspect --format '{{json .RepoDigests}}' @@IMG_Q@@ 2>/dev/null)"'''

PREP_CONTAINER = r'''mkdir -p "$G/ctr" "$G/meta"
cat >"$G/ctr/setup.sh" <<'GTU_SETUP_EOF'
@@SETUP@@
GTU_SETUP_EOF
cat >"$G/ctr/entry.sh" <<'GTU_ENTRY_EOF'
#!/bin/bash
set -u
# images opencfd/openfoam-*: their entrypoint (replaced by bash here) sources the OpenFOAM environment
for _rc in /openfoam/bash.rc; do [ -f "$_rc" ] && { set +u; . "$_rc" >/dev/null 2>&1; set -u; }; done
T=$(date +%s)
if ! bash -e @@CR@@/ctr/setup.sh >"$GTU_OUT/_container_setup.log" 2>&1; then
  echo "[gtu] container_setup failed - see _container_setup.log"; exit 20
fi
echo $(( $(date +%s) - T )) >@@CR@@/meta/setup_s
touch @@CR@@/meta/setup_ok
exec bash @@CR@@/ctr/cmd.sh
GTU_ENTRY_EOF
chmod +x "$G/ctr/entry.sh"'''

RUN_HOST = r'''timeout -k 60 "${REMAIN}s" bash "$G/cmd.sh" >"$OUT/_stdout.log" 2>&1
RC=$?'''

# Контейнер: все ядра (--cpus), /dev/shm — половина памяти (MPI), клон и $GTU_OUT примонтированы,
# пользователь root (container_setup ставит пакеты), OpenMPI под root разрешён переменными.
# IMDS из контейнера недоступен (hop limit 1): ключей роли у команды нет, зона — из GTU_AZ_MACHINE.
# Без --init: есть ли docker-init в пакете docker al2023 — не проверено; контейнер снимает docker kill в finish.
RUN_CONTAINER = r'''read -r _ MEM_KB _ </proc/meminfo
timeout -k 60 "${REMAIN}s" docker run --rm --name gtu-job --user 0:0 --entrypoint /bin/bash \
  --cpus "$NPROC" --shm-size "$((MEM_KB / 2))k" -w @@CR@@/repo \
  -v "$G/repo:@@CR@@/repo" -v "$OUT:@@CR@@/out" -v "$G/ctr:@@CR@@/ctr:ro" -v "$G/meta:@@CR@@/meta" \
  -e GTU_OUT=@@CR@@/out -e GTU_RUN_ID -e GTU_JOB -e GTU_ATTEMPT -e GTU_AZ -e GTU_AZ_CHEAPEST -e GTU_AZ_MACHINE \
  -e GTU_ITYPE -e NPROC -e OMP_NUM_THREADS="$NPROC" -e DEBIAN_FRONTEND=noninteractive \
  -e OMPI_ALLOW_RUN_AS_ROOT=1 -e OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1 \
  @@IMG_Q@@ @@CR@@/ctr/entry.sh >"$OUT/_stdout.log" 2>&1
RC=$?
[ -f "$G/meta/setup_s" ] && echo "[gtu-time] container_setup $(cat "$G/meta/setup_s")s (inside the command time)"
if [ ! -f "$G/meta/setup_ok" ]; then PHASE=container; [ "$RC" -eq 0 ] && RC=20; exit "$RC"; fi'''

STOP_CONTAINER = "  timeout 30 docker kill gtu-job >/dev/null 2>&1"


def _fill(text: str, vals: dict, pattern: str) -> str:
    return re.sub(pattern, lambda m: vals[m.group(1)], text)


def build_userdata(job: Job, *, run_id: str, attempt: int, repo: str, branch: str,
                   bucket: str, key_prefix: str, region: str = REGION,
                   upload_reserve_min: int = 8, az: str = "", az_cheapest: str = "", itype: str = "",
                   _root: str = "/opt/gtu", _log: str = "/var/log/gtu-job.log", _ctr_root: str = "/gtu") -> str:
    """Скрипт первой загрузки машины. Первая команда — сторож `shutdown -P +max_min`.

    `_root`, `_log`, `_ctr_root` — только для локального теста скрипта на подменённых командах.
    """
    q = shlex.quote
    if not REPO_RE.match(repo):
        raise RunnerError("repo: ожидается https://github.com/<владелец>/<репозиторий>")
    if not BRANCH_RE.match(branch) or ".." in branch:
        raise RunnerError(f"недопустимое имя ветки {branch!r}")
    if not RUN_ID_RE.match(run_id):
        raise RunnerError(f"недопустимый run_id {run_id!r}")
    for z in (az, az_cheapest):
        if z and not AZ_RE.match(z):
            raise RunnerError(f"недопустимая зона {z!r}")
    img = IMAGES[job.image]
    uv = bool(job.python)
    pip_cmd = '"$UV" pip install -q --python "$G/venv/bin/python"' if uv else '"$G/venv/bin/pip" install -q'
    # Сначала база (git, venv, boto3) — чтобы при сбое пакетов задачи журнал всё равно выгрузился.
    if img["pkg"] == "apt":
        apt = "timeout {m}m apt-get -o DPkg::Lock::Timeout=600 -q"
        base = (f"{apt.format(m=20)} update || exit 10\n"
                f"{apt.format(m=20)} install -y git python3-venv python3-pip || exit 11")
        sys_pk = (f"{apt.format(m=40)} install -y {' '.join(q(p) for p in job.apt)} || exit 16"
                  if job.apt else "")
    else:
        base = f"timeout 20m dnf -q install -y git python3-pip{' docker' if job.container else ''} || exit 11"
        sys_pk = (f"timeout 40m dnf -q install -y {' '.join(q(p) for p in job.dnf)} || exit 16"
                  if job.dnf else "")
    pip_pk = (f'timeout 40m {pip_cmd} {" ".join(q(p) for p in job.pip)} || exit 17' if job.pip else "")
    req = ""
    if job.requirements:
        r = q(job.requirements)
        req = (f'phase requirements\n[ -f {r} ] || {{ echo "[gtu] requirements file not found: "{r}; exit 23; }}\n'
               f"timeout 40m {pip_cmd} -r {r} || exit 23")
    blocks = {
        "BASE": base,
        "VENV": VENV_UV if uv else VENV_SYSTEM,
        "PACKAGES": "\n".join(x for x in (sys_pk, pip_pk) if x) or ":",
        "DOCKER": DOCKER_BLOCK if job.container else "",
        "REQUIREMENTS": req,
        "PREP": PREP_CONTAINER if job.container else "",
        "CMD_FILE": '"$G/ctr/cmd.sh"' if job.container else '"$G/cmd.sh"',
        "RUN": RUN_CONTAINER if job.container else RUN_HOST,
        "STOP": STOP_CONTAINER if job.container else "",
    }
    budget_s = job.max_min * 60 - upload_reserve_min * 60
    if budget_s < 600:
        raise RunnerError(f"{job.name}: max_min слишком мал для запаса на выгрузку")
    vals = {
        "MAX_MIN": str(job.max_min), "RUN_ID": run_id, "JOB": job.name, "ATTEMPT": str(int(attempt)),
        "MODE": job.mode, "LOG": q(_log), "REGION": q(region), "BUDGET_S": str(budget_s), "ROOT": q(_root),
        "RUN_ID_Q": q(run_id), "JOB_Q": q(job.name), "BUCKET_Q": q(bucket), "KEY_Q": q(key_prefix),
        "AZ_Q": q(az), "AZ_CHEAP_Q": q(az_cheapest), "ITYPE_Q": q(itype),
        "UPLOAD_PY": UPLOAD_PY, "BRANCH_Q": q(branch), "REPO_Q": q(repo),
        "PY": q(job.python), "IMG_Q": q(job.container), "CR": q(_ctr_root),
        "COMMAND": job.command.rstrip("\n"), "SETUP": "\n".join(job.container_setup) or ":",
    }
    ud = _fill(_fill(USERDATA_TEMPLATE, blocks, r"@@@([A-Z_]+)@@@"), vals, r"@@([A-Z_]+)@@")
    assert_no_secret(ud, "user-data")
    if len(ud.encode("utf-8")) > USERDATA_LIMIT:
        raise RunnerError(f"{job.name}: user-data {len(ud.encode())} байт > {USERDATA_LIMIT}")
    return ud


# --------------------------------------------------------------------------- журнал
class Ledger:
    """Журнал JSONL: запуски, окончания, проверки бюджета. В журнале нет бакета и ключей."""

    def __init__(self, path, redact: Redactor | None = None, clock=time.time):
        self.path = Path(path)
        self.redact = redact or Redactor()
        self.clock = clock

    def append(self, ev: str, **kw) -> dict:
        rec = {"ts": round(self.clock(), 1), "ev": ev, **kw}
        line = self.redact(json.dumps(rec, ensure_ascii=False, sort_keys=True))
        assert_no_secret(line, "журнал")
        self.path.parent.mkdir(parents=True, exist_ok=True)
        with self.path.open("a", encoding="utf-8") as fh:
            fh.write(line + "\n")
        return rec

    def records(self) -> list:
        if not self.path.exists():
            return []
        out = []
        for line in self.path.read_text(encoding="utf-8").splitlines():
            try:
                out.append(json.loads(line))
            except json.JSONDecodeError:
                continue
        return out

    def summary(self, now: float | None = None) -> dict:
        """Израсходовано (закрытые машины) и верхняя граница по открытым.

        Закрытая: (спот в её зоне на момент запуска + диск + IPv4) × max(60 с, конец − запуск) —
        оценка, не граница: цена спота за время работы может подняться (до MaxPrice).
        Открытая (нет записи end): (MaxPrice + диск + IPv4) × весь отпущенный срок — верхняя
        граница, потому что сторож и сборщик гасят машину не позже этого срока.
        `now` (живой прогон после сверки): открытая машина старше отпущенного срока — значит, её не
        погасили ни сторож, ни сборщик, ни сверка — считается по фактическому времени, а не по сроку.
        Без `now` (status, plan, DryRun — сверки не было) — по сроку, как раньше.
        Версия 3: машина по требованию пишет свою цену в spot_usd_h и max_usd_h (журнал читает и версия 2);
        межрегиональная передача результатов (xfer_usd, разово) прибавляется и к закрытой, и к открытой.
        """
        recs = self.records()
        launches = {r["iid"]: r for r in recs if r.get("ev") == "launch"}
        ends = {r["iid"]: r for r in recs if r.get("ev") == "end"}
        spent = committed = 0.0
        open_ = []
        for iid, L in launches.items():
            extra = L["ebs_usd_h"] + L.get("ip_usd_h", 0.0)
            if iid in ends:
                dur = max(60.0, ends[iid]["t_end"] - L["t_launch"])
                spent += (L["spot_usd_h"] + extra) * dur / 3600.0 + L.get("xfer_usd", 0.0)
            else:
                hours = L["allot_min"] / 60.0
                if now is not None:
                    hours = max(hours, (now - L["t_launch"]) / 3600.0)
                committed += (L["max_usd_h"] + extra) * hours + L.get("xfer_usd", 0.0)
                open_.append(iid)
        return {"spent": spent, "committed": committed, "open": open_,
                "closed": len(launches) - len(open_)}


# --------------------------------------------------------------------------- конфигурация
@dataclass
class Config:
    bucket: str
    prefix: str
    instance_profile: str
    run_id: str
    tag: str = DEFAULT_TAG
    repo: str = DEFAULT_REPO
    branch: str = "main"
    region: str = REGION
    ceiling: float = 2000.0          # решение владельца 14.11.2026: до 2000 долл. со стоп-краном
    report_step: float = 500.0       # отчёт владельцу при пересечении каждых 500 долл.
    max_parallel: int = 6
    vcpu_quota: int = 384            # --max-vcpu: квота спота в vCPU (6 × c7a.16xlarge)
    quota_wait_polls: int = 15       # квотная ошибка без своих машин в работе: ждать столько опросов
    capacity_wait_polls: int = 15    # нет ёмкости спота ни в одной зоне ни на одном типе: столько повторов
    price_margin: float = 1.2
    price_window_h: float = 6.0
    ebs_usd_gb_month: float = 0.08
    ipv4_usd_h: float = 0.005            # публичный IPv4 подсети по умолчанию (ставка — допущение)
    grace_min: int = 10
    upload_reserve_min: int = 8
    poll_s: float = 60.0
    retries: int = 1
    avoid_s: float = 1800.0          # зона после сбоя ёмкости или прерывания — в конец очереди на столько секунд
    spent_before: float = 0.0
    ledger_path: Path = HERE / "ledger.jsonl"
    out_dir: Path = HERE / "runs"
    dryrun: bool = False
    offline_price: float | None = None   # режим plan: цена спота задана вручную, AWS не вызывается
    offline_zones: dict | None = None    # режим plan: цены по зонам (допущение) для порядка зон
    # ---- версия 3 (без этих полей — точно версия 2)
    regions: list = field(default_factory=list)          # разрешённые регионы машин; [] — только region
    vcpu_quota_region: dict = field(default_factory=dict)      # явная квота спота по регионам, vCPU
    od_vcpu_quota_region: dict = field(default_factory=dict)   # явная квота по требованию по регионам, vCPU
    od_vcpu_quota: int = 384             # запасная квота по требованию в region (данные владельца: 384)
    vcpu_quota_fallback: int = QUOTA_FALLBACK   # запасная квота соседнего региона (спот и по требованию)
    on_demand_fallback: bool = False     # спот недоступен везде — машина по требованию
    region_price_margin: float = 0.05    # соседний регион берётся, только если дешевле с этим запасом
    xfer_gb: float = 10.0                # объём результатов задачи для сметы межрегиональной передачи, ГБ
    price_ref: Path = PRICE_REF          # справочник цен по требованию и передачи
    max_parallel_auto: bool = False      # --max-parallel не задан: в режиме регионов — по сумме квот
    offline_no_spot: bool = False        # режим plan: спот недоступен (проверка пути по требованию)

    @property
    def region_list(self) -> list:
        return list(self.regions) or [self.region]

    @property
    def multi(self) -> bool:
        """Режим версии 3: регионы сверх region или запас по требованию."""
        return self.region_list != [self.region] or self.on_demand_fallback


def parse_s3_uri(uri: str) -> tuple:
    m = re.match(r"^s3://([a-z0-9][a-z0-9.-]{1,61}[a-z0-9])/(.+)$", uri or "")
    if not m or not m.group(2).strip("/"):
        raise RunnerError("ALPHA_O_S3 должен быть вида s3://<бакет>/<префикс> (префикс обязателен)")
    return m.group(1), m.group(2).strip("/")


def parse_zone_prices(text: str) -> dict:
    """'us-east-1a=0.955,us-east-1c=0.706' -> {зона: цена}."""
    out = {}
    for part in (text or "").split(","):
        if not part.strip():
            continue
        az, _, p = part.partition("=")
        az = az.strip()
        try:
            v = float(p)
        except ValueError:
            raise RunnerError(f"--assume-zone-prices: цена {p!r} не число") from None
        if not AZ_RE.match(az) or not v > 0:
            raise RunnerError(f"--assume-zone-prices: ожидается зона=цена, получено {part.strip()!r}")
        out[az] = v
    if not out:
        raise RunnerError("--assume-zone-prices: пусто")
    return out


def parse_regions(text: str) -> list:
    """'us-east-1,us-east-2' -> ['us-east-1', 'us-east-2'] (порядок сохраняется, повторы отброшены)."""
    out = []
    for part in (text or "").split(","):
        r = part.strip()
        if not r:
            continue
        if not REGION_RE.match(r):
            raise RunnerError(f"--regions: недопустимый регион {r!r}")
        if r not in out:
            out.append(r)
    if not out:
        raise RunnerError("--regions: пусто")
    return out


def parse_region_ints(text: str, flag: str) -> dict:
    """'us-east-1=384,us-east-2=256' -> {регион: vCPU}."""
    out = {}
    for part in (text or "").split(","):
        if not part.strip():
            continue
        r, _, v = part.partition("=")
        r = r.strip()
        try:
            n = int(v)
        except ValueError:
            raise RunnerError(f"{flag}: значение {v!r} не целое") from None
        if not REGION_RE.match(r) or n < 0:
            raise RunnerError(f"{flag}: ожидается регион=vCPU, получено {part.strip()!r}")
        out[r] = n
    if not out:
        raise RunnerError(f"{flag}: пусто")
    return out


def load_price_ref(path) -> dict:
    """Справочник цен (версия 3): по требованию {тип: {регион: долл./ч}}, передача {"откуда>куда" | default: долл./ГБ}."""
    try:
        data = json.loads(Path(path).read_text(encoding="utf-8"))
        od, xf = data["on_demand_usd_h"], data["xfer_usd_gb"]
        od = {t: {r: float(p) for r, p in m.items()} for t, m in od.items()}
        xf = {k: float(v) for k, v in xf.items()}
    except (OSError, ValueError, KeyError, TypeError, AttributeError) as e:
        raise RunnerError(f"справочник цен {Path(path).name} не прочитан: {type(e).__name__}") from None
    if "default" not in xf or any(p <= 0 for m in od.values() for p in m.values()) or any(v < 0 for v in xf.values()):
        raise RunnerError(f"справочник цен {Path(path).name}: нужны цены > 0 и xfer_usd_gb.default")
    return {"on_demand_usd_h": od, "xfer_usd_gb": xf, "date": str(data.get("date", "?")),
            "status": str(data.get("status", ""))}


def new_run_id(clock=time.time) -> str:
    t = dt.datetime.fromtimestamp(clock(), tz=dt.timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    return f"{t}-{_secrets.token_hex(2)}"


@dataclass
class Active:
    job: Job
    attempt: int
    iid: str
    itype: str
    az: str
    t_launch: float
    deadline: float
    az_req: str = ""          # зона, в которую просили
    az_rank: int = -1         # место зоны по цене (0 — самая дешёвая)
    az_cheapest: str = ""     # самая дешёвая зона типа на момент запуска
    spot_usd_h: float = 0.0   # цена спота в фактической зоне на момент запуска (по требованию — его цена)
    t_mark: float = 0.0       # когда впервые увидели метку с кодом 143 при ещё живой машине
    region: str = ""          # версия 3: регион машины ("" — region конфигурации)
    market: str = "spot"      # версия 3: spot | on-demand


@dataclass
class _Try:
    """Итоги попыток запуска одной задачи (по типам, регионам, рынкам) — для статуса и пояснения."""
    refused: bool = False
    no_price: list = field(default_factory=list)       # типы без цены спота ни в одном регионе
    cap: list = field(default_factory=list)
    skipped: list = field(default_factory=list)
    quota: list = field(default_factory=list)          # версия 3: квота (ошибка, пауза региона, нет свободной)
    denied: list = field(default_factory=list)         # версия 3: отказ в правах в регионе
    failed: list = field(default_factory=list)         # версия 3: иная ошибка запуска в регионе
    image_bad: list = field(default_factory=list)      # версия 3: образ недоступен в регионе
    image_ok: bool = False
    od_no_price: list = field(default_factory=list)    # версия 3: нет цены по требованию в справочнике
    ami: dict = field(default_factory=dict)            # регион -> образ (или ImageUnavailable)

    def spot_unavailable(self) -> bool:
        """Спот не дал машину по причине, которую машина по требованию может обойти."""
        if self.refused or (self.image_bad and not self.image_ok):
            return False
        return bool(self.cap or self.quota or self.no_price or self.denied)


@dataclass
class PriceInfo:
    itype: str
    latest: dict          # зона -> последняя цена
    window_max: float     # максимум за окно по всем зонам

    def az_price(self, az: str) -> float:
        return self.latest.get(az, max(self.latest.values()) if self.latest else self.window_max)

    def by_price(self) -> list:
        """[(зона, цена)] по возрастанию цены; псевдозоны режима plan отброшены."""
        return sorted(((z, p) for z, p in self.latest.items() if AZ_RE.match(z)), key=lambda x: (x[1], x[0]))


# --------------------------------------------------------------------------- исполнитель
class Runner:
    def __init__(self, cfg: Config, ec2=None, s3=None, ssm=None, sts=None,
                 clock=time.time, sleep=time.sleep, out=print, region_clients=None):
        """`region_clients(регион) -> {"ec2", "ssm", "servicequotas"}` — клиенты регионов, кроме cfg.region
        (для него — ec2/ssm из аргументов). Вызывается лениво: только когда регион действительно нужен."""
        self.cfg, self.ec2, self.s3, self.ssm, self.sts = cfg, ec2, s3, ssm, sts
        self.clock, self.sleep, self._out = clock, sleep, out
        self.region_clients = region_clients
        self.regions = cfg.region_list
        self.v3 = cfg.multi
        self._rc: dict = {}             # регион -> клиенты
        self._quota: dict = {}          # (регион, рынок) -> (vCPU, источник)
        self._qblock: dict = {}         # (регион, рынок) -> до какого времени пауза после квотной ошибки
        self._denied: dict = {}         # (регион, рынок) -> код отказа в правах
        self._iid_region: dict = {}     # машина -> регион
        self._ledger_regions = None
        self._ref = None                # справочник цен (версия 3), читается лениво
        self.redact = Redactor(cfg.bucket)
        self.ledger = Ledger(cfg.ledger_path, self.redact, clock)          # общий учёт
        side = cfg.dryrun or cfg.offline_price is not None
        self.journal = (Ledger(cfg.out_dir / cfg.run_id / "journal.dryrun.jsonl", self.redact, clock)
                        if side else self.ledger)                          # куда пишем
        self.active: dict = {}
        self.results: dict = {}
        self.virtual = 0.0          # DRYRUN/plan: прогноз уже «запущенных» задач пакета
        self.last_params = None
        self._prices: dict = {}     # (регион, тип) -> (время, PriceInfo)
        self._ami: dict = {}        # (регион, образ) -> id | ImageUnavailable
        self._subnets: dict = {}    # регион -> dict зона -> подсеть; False — отказ в правах (нет ключа — не читались)
        self._avoid: dict = {}      # (тип, зона) -> до какого времени зона в конце очереди
        self._why: dict = {}        # задача -> пояснение к статусу
        self._quota_wait: dict = {}
        self._cap_wait: dict = {}
        self._budget_wait: dict = {}
        # Гашение не подтверждено: машина, возможно, жива и занимает квоту — новых машин прогон не запускает.
        self._unconfirmed: dict = {}
        self._reconciled = False    # сверка журнала прошла: бюджет считает открытые машины по факту, если дольше срока
        self._closed_n = 0          # счётчик закрытых машин: бюджет пересчитывается, когда он меняется
        self._seq = 0
        # ClientToken уникален в пределах процесса: повторный запуск с тем же --run-id не получит
        # чужой (старый) экземпляр или IdempotentParameterMismatch по совпавшему токену.
        self._nonce = _secrets.token_hex(4)

    # ---- вывод
    def say(self, msg: str) -> None:
        self._out(self.redact(msg))

    @property
    def offline(self) -> bool:
        return self.cfg.offline_price is not None

    def key_prefix(self, job_name: str, run_id: str | None = None, attempt: int = 1) -> str:
        """Префикс результатов попытки: первая — <run>/<job>, повтор N — <run>/_attemptN/<job>.

        Имя задачи не начинается с «_», поэтому префикс повтора не совпадёт ни с одной задачей,
        а файлы прерванной попытки не смешаются с файлами следующей.
        """
        base = f"{self.cfg.prefix}/gtu/{run_id or self.cfg.run_id}"
        return f"{base}/{job_name}" if attempt <= 1 else f"{base}/_attempt{int(attempt)}/{job_name}"

    # ---- регионы (версия 3)
    def client(self, region: str | None, kind: str):
        """Клиент EC2 / SSM / Service Quotas региона. Для cfg.region — переданные в конструктор (как v2)."""
        region = region or self.cfg.region
        if region == self.cfg.region and kind in ("ec2", "ssm") and getattr(self, kind) is not None:
            return getattr(self, kind)
        if region not in self._rc:
            if self.region_clients is None:
                raise RunnerError(f"нет клиентов AWS для региона {region}")
            self._rc[region] = self.region_clients(region)
        return self._rc[region].get(kind)

    def _rf(self, region: str | None, market: str | None = None) -> dict:
        """Поля журнала версии 3 (регион, рынок); в режиме версии 2 — пусто: журнал как прежде."""
        if not self.v3:
            return {}
        out = {"region": region or self.cfg.region}
        if market:
            out["market"] = market
        return out

    def _where(self, region: str | None) -> str:
        return f" @{region or self.cfg.region}" if self.v3 else ""

    def _region_of(self, iid: str) -> str:
        for d in (self.active, self._unconfirmed):
            a = d.get(iid)
            if a is not None and a.region:
                return a.region
        if iid in self._iid_region:
            return self._iid_region[iid]
        if self._ledger_regions is None:
            self._ledger_regions = {r["iid"]: r["region"] for r in self.ledger.records()
                                    if r.get("ev") == "launch" and r.get("region") and r.get("iid")}
        return self._ledger_regions.get(iid, self.cfg.region)

    def _by_region(self, iids) -> dict:
        """{регион: [id в исходном порядке]}; в режиме версии 2 — один регион."""
        out: dict = {}
        for iid in iids:
            out.setdefault(self._region_of(iid), []).append(iid)
        return out

    @property
    def ref(self) -> dict:
        if self._ref is None:
            self._ref = load_price_ref(self.cfg.price_ref)
        return self._ref

    def od_prices(self, itype: str) -> dict:
        """Цены по требованию типа по разрешённым регионам из справочника: {регион: долл./ч}."""
        m = self.ref["on_demand_usd_h"].get(itype, {})
        return {r: m[r] for r in self.regions if r in m}

    def xfer_usd(self, job: Job, region: str | None) -> float:
        """Смета передачи результатов из региона машины в бакет (region): объём — min(диск, --xfer-gb)."""
        region = region or self.cfg.region
        if not self.v3 or region == self.cfg.region:
            return 0.0
        xf = self.ref["xfer_usd_gb"]
        rate = xf.get(f"{region}>{self.cfg.region}", xf["default"])
        return min(float(job.disk_gb), self.cfg.xfer_gb) * rate

    def quota(self, region: str | None, market: str = "spot") -> int:
        region = region or self.cfg.region
        if not self.v3:
            return self.cfg.vcpu_quota
        if (region, market) not in self._quota:
            self._quota[(region, market)] = self._read_quota(region, market)
        return self._quota[(region, market)][0]

    def _read_quota(self, region: str, market: str) -> tuple:
        """(vCPU, источник): флаг → service-quotas → запасное значение (только чтение, без изменений)."""
        cfg = self.cfg
        explicit = cfg.vcpu_quota_region if market == "spot" else cfg.od_vcpu_quota_region
        if region in explicit:
            return int(explicit[region]), "флаг"
        if region == cfg.region:
            fb = cfg.vcpu_quota if market == "spot" else cfg.od_vcpu_quota
        else:
            fb = cfg.vcpu_quota_fallback
        if self.offline:
            return fb, "запасное (plan: service-quotas не вызывается)"
        try:
            sq = self.client(region, "servicequotas")
            if sq is None:
                return fb, "запасное (нет клиента service-quotas)"
            v = sq.get_service_quota(ServiceCode="ec2", QuotaCode=QUOTA_CODES[market])["Quota"]["Value"]
            return int(v), f"service-quotas {QUOTA_CODES[market]}"
        except Exception as e:      # нет прав (у ключа агента их нет), сеть, неизвестный ответ — запасное
            resp = getattr(e, "response", None)
            code = resp.get("Error", {}).get("Code") if isinstance(resp, dict) else None
            return fb, f"запасное (service-quotas: {code or type(e).__name__})"

    def markets(self) -> list:
        return ["spot", "on-demand"] if self.cfg.on_demand_fallback else ["spot"]

    def load_quotas(self) -> None:
        """Версия 3: квоты по регионам — в начале прогона и plan, с источником; --max-parallel по сумме квот."""
        if not self.v3:
            return
        rows = []
        for r in self.regions:
            for m in self.markets():
                self.quota(r, m)
                v, src = self._quota[(r, m)]
                rows.append([r, m, v, src])
                self.say(f"[квоты] {r} {'спот' if m == 'spot' else 'по требованию'}: {v} vCPU ({src})")
        self.journal.append("quotas", run_id=self.cfg.run_id, quotas=rows)
        if self.cfg.max_parallel_auto:
            n = max(self.cfg.max_parallel, sum(self.quota(r, "spot") for r in self.regions) // 64)
            if n != self.cfg.max_parallel:
                self.say(f"[очередь] --max-parallel не задан: {n} машин (сумма квот спота по регионам / 64 vCPU)")
            self.cfg.max_parallel = n

    def total_quota(self) -> int:
        if not self.v3:
            return self.cfg.vcpu_quota
        return sum(self.quota(r, m) for r in self.regions for m in self.markets())

    def _used(self, region: str, market: str = "spot") -> int:
        return sum(a.job.vcpu for a in self.active.values()
                   if (a.region or self.cfg.region) == region and a.market == market)

    def _fits(self, job: Job) -> bool:
        """Есть ли регион (и рынок), где задача помещается в квоту вместе с машинами прогона."""
        if not self.v3:
            return sum(a.job.vcpu for a in self.active.values()) + job.vcpu <= self.cfg.vcpu_quota
        return any(self._used(r, m) + job.vcpu <= self.quota(r, m) for r in self.regions for m in self.markets())

    def _blocked(self, region: str, market: str) -> str:
        """Почему регион сейчас не годится для рынка: отказ в правах / пауза после квотной ошибки / ''."""
        if (region, market) in self._denied:
            return f"нет прав ({self._denied[(region, market)]})"
        if self._qblock.get((region, market), 0.0) > self.clock():
            return "квота (пауза после ошибки)"
        return ""

    # ---- цены и бюджет
    def price(self, itype: str, region: str | None = None) -> PriceInfo:
        region = region or self.cfg.region
        now = self.clock()
        hit = self._prices.get((region, itype))
        if hit and now - hit[0] < 300:
            return hit[1]
        if self.offline:
            p = float(self.cfg.offline_price)
            zp = dict(self.cfg.offline_zones or {})
            if self.v3:                      # зоны своего региона (имя зоны = регион + буква)
                zp = {z: v for z, v in zp.items() if z[:-1] == region}
            info = PriceInfo(itype, zp or {"assumed": p}, max([p] + list(zp.values())))
            if self.cfg.offline_no_spot:     # plan: спот недоступен — проверка пути по требованию
                info = PriceInfo(itype, {}, 0.0)
        else:
            end = dt.datetime.fromtimestamp(now, tz=dt.timezone.utc)
            start = end - dt.timedelta(hours=self.cfg.price_window_h)
            kw = dict(InstanceTypes=[itype], ProductDescriptions=[PRODUCT], StartTime=start, EndTime=end)
            latest, ts, wmax = {}, {}, 0.0
            ec2 = self.client(region, "ec2")
            while True:
                r = ec2.describe_spot_price_history(**kw)
                for e in r.get("SpotPriceHistory", []):
                    if e.get("InstanceType") != itype:
                        continue
                    p, az, t = float(e["SpotPrice"]), e["AvailabilityZone"], e["Timestamp"]
                    wmax = max(wmax, p)
                    if az not in ts or t > ts[az]:
                        ts[az], latest[az] = t, p
                tok = r.get("NextToken")
                if not tok:
                    break
                kw["NextToken"] = tok
            info = PriceInfo(itype, latest, wmax)
        self._prices[itype] = (now, info)
        return info

    def max_price(self, info: PriceInfo) -> float:
        return math.ceil(info.window_max * self.cfg.price_margin * 10000) / 10000.0

    def ebs_usd_h(self, job: Job) -> float:
        return job.disk_gb * self.cfg.ebs_usd_gb_month / 730.0

    def allot_min(self, job: Job) -> int:
        return job.max_min + self.cfg.grace_min

    def budget_check(self, job: Job, itype: str, info: PriceInfo) -> tuple:
        """(разрешено?, MaxPrice, цифры). Пишет запись budget в журнал."""
        s = self.ledger.summary(now=self.clock() if self._reconciled else None)
        maxp = self.max_price(info)
        new = (maxp + self.ebs_usd_h(job) + self.cfg.ipv4_usd_h) * self.allot_min(job) / 60.0
        spent = s["spent"] + self.cfg.spent_before
        total = spent + s["committed"] + self.virtual + new
        ok = info.window_max > 0 and total <= self.cfg.ceiling + 1e-9
        nums = {"spent": round(spent, 4), "committed": round(s["committed"], 4),
                "virtual": round(self.virtual, 4), "new": round(new, 4), "total": round(total, 4),
                "ceiling": self.cfg.ceiling}
        self.journal.append("budget", run_id=self.cfg.run_id, job=job.name, itype=itype,
                            spot_usd_h=round(max(info.latest.values()) if info.latest else 0.0, 4),
                            max_usd_h=maxp, ok=ok, **nums)
        verdict = "можно" if ok else "ОТКАЗ: превышение потолка"
        step = self.cfg.report_step
        before = total - new
        if ok and step > 0 and int(total // step) > int(before // step):
            порог = int(total // step) * step
            self.journal.append("budget_report", run_id=self.cfg.run_id, job=job.name, threshold=порог,
                                total=round(total, 4), ceiling=self.cfg.ceiling)
            self.say(f"[бюджет] ОТЧЁТ ВЛАДЕЛЬЦУ: с этой задачей расходы пересекают {порог:.0f} долл. "
                     f"(всего {total:.2f} из {self.cfg.ceiling:.0f})")
        self.say(f"[бюджет] {job.name} на {itype}: израсходовано {spent:.2f} + в работе "
                 f"{s['committed']:.2f} + пакет {self.virtual:.2f} + новая {new:.2f} = {total:.2f} "
                 f"из {self.cfg.ceiling:.2f} долл. (MaxPrice {maxp:.4f}/ч) — {verdict}")
        return ok, maxp, nums

    # ---- зоны и подсети
    def _avoided(self, itype: str, az: str) -> bool:
        return self._avoid.get((itype, az), 0.0) > self.clock()

    def zone_order(self, itype: str, info: PriceInfo) -> list:
        """Порядок попыток: по цене, зоны со свежим сбоем ёмкости или прерыванием — в конце."""
        zs = info.by_price()
        return [z for z in zs if not self._avoided(itype, z[0])] + [z for z in zs if self._avoided(itype, z[0])]

    def _say_zones(self, job: Job, itype: str, order: list, by_price: list, maxp: float) -> None:
        if not order:
            hint = " (в plan — --assume-zone-prices)" if self.offline else ""
            self.say(f"[зоны] {job.name} {itype}: цены по зонам неизвестны{hint} — зону выберет EC2")
            return
        txt = " → ".join(f"{z} {p:.4f}" + (" (после сбоя)" if self._avoided(itype, z) else "") for z, p in order)
        tag = " (допущение plan)" if self.offline else ""
        self.say(f"[зоны] {job.name} {itype}: {txt}{tag}; ожидается {order[0][0]} по {order[0][1]:.4f} долл./ч, "
                 f"самая дешёвая {by_price[0][0]}; MaxPrice {maxp:.4f} и бюджет — по максимуму цены за окно")

    def subnets(self):
        """Подсети по умолчанию по зонам (describe_subnets, default-for-az=true); False — нет прав."""
        if self._subnets is not None:
            return self._subnets
        from botocore.exceptions import ClientError
        try:
            out, kw = {}, {"Filters": [{"Name": "default-for-az", "Values": ["true"]}]}
            while True:
                r = self.ec2.describe_subnets(**kw)
                for s in r.get("Subnets", []):
                    if s.get("State", "available") == "available" and s.get("AvailabilityZone"):
                        out[s["AvailabilityZone"]] = {"id": s["SubnetId"],
                                                      "public": bool(s.get("MapPublicIpOnLaunch", False))}
                if not r.get("NextToken"):
                    break
                kw["NextToken"] = r["NextToken"]
            self._subnets = out
            self.say(f"[подсети] по умолчанию в зонах: {', '.join(sorted(out)) or 'нет'}")
        except ClientError as e:
            code = e.response.get("Error", {}).get("Code", "?")
            self._subnets = False
            self.journal.append("subnets", run_id=self.cfg.run_id, ok=False, code=code)
            self.say(f"[подсети] describe_subnets: {code} — зона задаётся Placement.AvailabilityZone "
                     f"(подсеть по умолчанию этой зоны выберет EC2)")
        return self._subnets

    def placement(self, az: str | None) -> tuple:
        """('subnet', подсеть) | ('az', None) | ('none', None) | ('skip', None)."""
        if not az:
            return "none", None
        if self.offline:
            return "subnet", {"id": f"<подсеть по умолчанию {az}>", "public": True}
        sn = self.subnets()
        if not sn:                       # нет прав или нет подсетей по умолчанию — зона через Placement
            return "az", None
        if az in sn:
            return "subnet", sn[az]
        return "skip", None

    # ---- параметры запуска
    def ami(self, image: str) -> str:
        if self.offline:
            return f"<ssm:{IMAGES[image]['ssm']}>"
        hit = self._ami.get(image)
        if isinstance(hit, ImageUnavailable):
            raise hit
        if hit:
            return hit
        from botocore.exceptions import ClientError
        img = IMAGES[image]
        try:
            self._ami[image] = self.ssm.get_parameter(Name=img["ssm"])["Parameter"]["Value"]
            return self._ami[image]
        except ClientError as e:
            code = e.response.get("Error", {}).get("Code", "?")
        err = None
        if img.get("owner"):
            self.say(f"[образ] {image}: ssm:GetParameter {code} — ищу образ владельца {img['owner']} (DescribeImages)")
            try:
                r = self.ec2.describe_images(Owners=[img["owner"]], Filters=[
                    {"Name": "name", "Values": [img["name"]]}, {"Name": "state", "Values": ["available"]},
                    {"Name": "architecture", "Values": ["x86_64"]}])
                found = sorted(r.get("Images", []), key=lambda i: i.get("CreationDate", ""))
                if found:
                    self._ami[image] = found[-1]["ImageId"]
                    self.say(f"[образ] {image}: {found[-1].get('Name', '')}")
                    return self._ami[image]
                err = f"ssm:GetParameter {code}, DescribeImages — образов нет"
            except ClientError as e2:
                err = f"ssm:GetParameter {code}, ec2:DescribeImages {e2.response.get('Error', {}).get('Code', '?')}"
        self._ami[image] = ImageUnavailable(
            f"образ {image} недоступен ({err or 'ssm:GetParameter ' + code}) — используйте image al2023 "
            f"или container ubuntu:24.04 (пакеты apt — командами container_setup)")
        raise self._ami[image]

    def build_params(self, job: Job, itype: str, ami: str, maxp: float, attempt: int, *,
                     az: str | None = None, subnet: dict | None = None, az_cheapest: str = "") -> dict:
        cfg = self.cfg
        ud = build_userdata(job, run_id=cfg.run_id, attempt=attempt, repo=cfg.repo, branch=cfg.branch,
                            bucket=cfg.bucket, key_prefix=self.key_prefix(job.name, attempt=attempt),
                            region=cfg.region, upload_reserve_min=cfg.upload_reserve_min,
                            az=az or "", az_cheapest=az_cheapest, itype=itype)
        tags = [{"Key": "Project", "Value": cfg.tag}, {"Key": "RunId", "Value": cfg.run_id},
                {"Key": "Job", "Value": job.name}, {"Key": "Name", "Value": f"gtu-{job.name}"}]
        self._seq += 1
        token = "gtu-" + hashlib.sha256(
            f"{cfg.run_id}|{job.name}|{attempt}|{itype}|{az}|{self._seq}|{self._nonce}".encode()).hexdigest()[:40]
        params = dict(
            ImageId=ami, InstanceType=itype, MinCount=1, MaxCount=1,
            IamInstanceProfile={"Name": cfg.instance_profile},
            InstanceInitiatedShutdownBehavior="terminate",
            InstanceMarketOptions={"MarketType": "spot",
                                   "SpotOptions": {"MaxPrice": f"{maxp:.4f}",
                                                   "SpotInstanceType": "one-time",
                                                   "InstanceInterruptionBehavior": "terminate"}},
            MetadataOptions={"HttpTokens": "required", "HttpPutResponseHopLimit": 1,
                             "HttpEndpoint": "enabled"},
            BlockDeviceMappings=[{"DeviceName": IMAGES[job.image]["root"],
                                  "Ebs": {"VolumeSize": job.disk_gb, "VolumeType": "gp3",
                                          "DeleteOnTermination": True}}],
            TagSpecifications=[{"ResourceType": "instance", "Tags": tags},
                               {"ResourceType": "volume", "Tags": tags}],
            UserData=ud,
            ClientToken=token,
        )
        if subnet:
            if subnet["public"]:          # подсеть по умолчанию выдаёт публичный IPv4 сама — как раньше
                params["SubnetId"] = subnet["id"]
            else:                         # иначе публичный IPv4 явно: машине нужен интернет (git, pip, S3)
                params["NetworkInterfaces"] = [{"DeviceIndex": 0, "SubnetId": subnet["id"],
                                                "AssociatePublicIpAddress": True, "DeleteOnTermination": True}]
        elif az:
            params["Placement"] = {"AvailabilityZone": az}
        return params

    # ---- запуск одной задачи
    def launch(self, job: Job, attempt: int = 1) -> tuple:
        """('launched', Active) | ('refused'|'defer'|'capacity'|'no_price'|'failed'|'image_unavailable'|
        'dryrun_ok'|'dryrun_denied'|'planned', None).

        'capacity' — ни одна зона ни одного типа не приняла (нет ёмкости спота): очередь повторит позже;
        'no_price' — ни у одного типа нет истории цены: без цены нет ни MaxPrice, ни бюджета.
        Причина отказа — в self._why[имя задачи] (уходит в поле detail итога).
        """
        from botocore.exceptions import ClientError
        refused, ami = False, None
        no_price, cap, skipped = [], [], []
        self._why.pop(job.name, None)             # пояснение прошлой попытки к этой не относится
        for itype in job.types:
            info = self.price(itype)
            if info.window_max <= 0:
                no_price.append(itype)
                self.journal.append("launch_error", run_id=self.cfg.run_id, job=job.name, itype=itype,
                                    attempt=attempt, code="NoSpotPrice")
                self.say(f"[цена] нет истории цены спота для {itype} — тип пропущен")
                continue
            ok, maxp, nums = self.budget_check(job, itype, info)
            if not ok:
                refused = True
                continue
            if ami is None:                       # образ — после бюджета: отказ по бюджету не зовёт SSM
                try:
                    ami = self.ami(job.image)
                except ImageUnavailable as e:
                    self._why[job.name] = str(e)
                    self.journal.append("launch_error", run_id=self.cfg.run_id, job=job.name, attempt=attempt,
                                        code="ImageUnavailable", detail=str(e))
                    self.say(f"[запуск] {job.name}: {e}")
                    return "image_unavailable", None
            by_price = info.by_price()
            rank = {z: i for i, (z, _) in enumerate(by_price)}
            cheapest = by_price[0][0] if by_price else ""
            order = self.zone_order(itype, info)
            self._say_zones(job, itype, order, by_price, maxp)
            for az, _p in (order or [(None, None)]):
                how, subnet = self.placement(az)
                if how == "skip":
                    skipped.append(f"{itype}@{az}")
                    self.journal.append("zone_skip", run_id=self.cfg.run_id, job=job.name, itype=itype, az=az,
                                        attempt=attempt, reason="no_default_subnet")
                    self.say(f"[зоны] {job.name}: в {az} нет подсети по умолчанию — зона пропущена")
                    continue
                params = self.build_params(job, itype, ami, maxp, attempt, az=az, subnet=subnet,
                                           az_cheapest=cheapest)
                self.last_params = params
                if self.offline:
                    self.virtual += nums["new"]
                    return "planned", None
                if self.cfg.dryrun:
                    try:
                        rr = self.ec2.run_instances(DryRun=True, **params)
                    except ClientError as e:
                        code = e.response.get("Error", {}).get("Code", "?")
                        if code == "DryRunOperation":
                            self.virtual += nums["new"]
                            self.say(f"[dryrun] {job.name} на {itype} в {az or 'зоне EC2'}: запуск разрешён (DryRunOperation)")
                            self.journal.append("dryrun", run_id=self.cfg.run_id, job=job.name, itype=itype,
                                                az=az or "", ok=True)
                            return "dryrun_ok", None
                        self.say(f"[dryrun] {job.name} на {itype} в {az or 'зоне EC2'}: отказ {code}")
                        self.journal.append("dryrun", run_id=self.cfg.run_id, job=job.name, itype=itype,
                                            az=az or "", ok=False, code=code)
                        if code in CAPACITY_ERRORS:
                            cap.append(f"{itype}@{az or 'EC2'}:{code}")
                            continue
                        return "dryrun_denied", None
                    # DryRun проигнорирован и машина создана: в режиме DryRun сборщика в finally нет — гасим здесь
                    ids = [i.get("InstanceId") for i in (rr or {}).get("Instances", [])]
                    self.say(f"[dryrun] {job.name}: неожиданный успех без DryRunOperation — гашу созданное {ids}")
                    self.journal.append("dryrun", run_id=self.cfg.run_id, job=job.name, itype=itype,
                                        az=az or "", ok=False, code="UnexpectedSuccess")
                    self._terminate(ids, "DryRun проигнорирован")
                    return "dryrun_denied", None
                try:
                    r = self.ec2.run_instances(**params)
                except ClientError as e:
                    code = e.response.get("Error", {}).get("Code", "?")
                    self.journal.append("launch_error", run_id=self.cfg.run_id, job=job.name, itype=itype,
                                        az=az or "", attempt=attempt, code=code)
                    self.say(f"[запуск] {job.name} на {itype} в {az or 'зоне EC2'}: {code}")
                    if code in CAPACITY_ERRORS:
                        cap.append(f"{itype}@{az or 'EC2'}:{code}")
                        if az:
                            self._avoid[(itype, az)] = self.clock() + self.cfg.avoid_s
                        continue                            # следующая зона, затем следующий тип
                    if code in QUOTA_ERRORS:
                        return "defer", None
                    # Иная ошибка (InternalError / Unavailable после повторов botocore): EC2 мог машину
                    # всё же создать — ищем по ClientToken, чтобы она сразу вошла в учёт, квоту и бюджет.
                    r = self._find_by_token(params["ClientToken"])
                    if r is None:
                        self._why[job.name] = f"RunInstances: {code}"
                        return "failed", None
                    self.say(f"[запуск] {job.name}: несмотря на {code} машина создана — беру в учёт")
                inst = r["Instances"][0]
                iid = inst["InstanceId"]
                az_act = inst.get("Placement", {}).get("AvailabilityZone", "") or (az or "")
                spot = info.az_price(az_act)
                now = self.clock()
                act = Active(job, attempt, iid, itype, az_act, now, now + self.allot_min(job) * 60,
                             az_req=az or "", az_rank=rank.get(az_act, -1), az_cheapest=cheapest, spot_usd_h=spot)
                self.ledger.append("launch", run_id=self.cfg.run_id, job=job.name, attempt=attempt, iid=iid,
                                   itype=itype, az=az_act, az_req=az or "", az_rank=act.az_rank,
                                   az_cheapest=cheapest, zones=[[z, round(p, 4)] for z, p in by_price],
                                   # как задана зона (id подсети не пишем: журнал коммитится в открытый репозиторий)
                                   net=("subnet" if subnet["public"] else "eni") if subnet else ("placement" if az else "auto"),
                                   image=job.image, mode=job.mode,
                                   spot_usd_h=round(spot, 4), max_usd_h=maxp,
                                   ebs_usd_h=round(self.ebs_usd_h(job), 5), ip_usd_h=self.cfg.ipv4_usd_h,
                                   t_launch=round(now, 1), allot_min=self.allot_min(job))
                where = "самая дешёвая" if az_act and az_act == cheapest else f"место по цене {act.az_rank + 1}"
                self.say(f"[запуск] {job.name} (попытка {attempt}): {iid} {itype} {az_act} ({where}, "
                         f"{spot:.4f} долл./ч), сторож {job.max_min} мин, сборщик через {self.allot_min(job)} мин")
                return "launched", act
            if order:
                self.say(f"[зоны] {job.name}: на {itype} ни одна зона не приняла — следующий тип")
        if refused:
            return "refused", None
        if len(no_price) == len(job.types):
            self._why[job.name] = (f"нет истории цены спота за {self.cfg.price_window_h:g} ч для "
                                   f"{', '.join(no_price)} — без цены машина не запускается "
                                   f"(MaxPrice и бюджет не посчитать)")
            return "no_price", None
        why = []
        if cap:
            why.append("нет ёмкости спота: " + ", ".join(cap))
        if skipped:
            why.append("нет подсети по умолчанию: " + ", ".join(skipped))
        if no_price:
            why.append("нет цены: " + ", ".join(no_price))
        if why:
            self._why[job.name] = "; ".join(why)
        return ("capacity" if cap else "failed"), None

    def _find_by_token(self, token: str):
        """Машина, созданная запуском с этим ClientToken (не terminated), в виде ответа run_instances; иначе None."""
        from botocore.exceptions import BotoCoreError, ClientError
        try:
            rr = self.ec2.describe_instances(Filters=[{"Name": "client-token", "Values": [token]}])
        except (ClientError, BotoCoreError):
            return None
        found = [i for res in rr.get("Reservations", []) for i in res.get("Instances", [])
                 if i.get("State", {}).get("Name") != "terminated"]
        return {"Instances": found[:1]} if found else None

    # ---- наблюдение
    def describe(self, iids: list) -> dict:
        out, kw = {}, {"Filters": [{"Name": "instance-id", "Values": list(iids)}]}
        while True:
            r = self.ec2.describe_instances(**kw)
            for res in r.get("Reservations", []):
                for i in res.get("Instances", []):
                    out[i["InstanceId"]] = {"state": i["State"]["Name"],
                                            "reason": i.get("StateReason", {}).get("Code", ""),
                                            "az": i.get("Placement", {}).get("AvailabilityZone", "")}
            if not r.get("NextToken"):
                return out
            kw["NextToken"] = r["NextToken"]

    def marker(self, job_name: str, run_id: str | None = None, attempt: int = 1) -> dict | None:
        from botocore.exceptions import ClientError
        try:
            body = self.s3.get_object(Bucket=self.cfg.bucket,
                                      Key=self.key_prefix(job_name, run_id, attempt) + "/_DONE.json")["Body"].read()
        except ClientError as e:
            code = e.response.get("Error", {}).get("Code", "")
            if code not in ("NoSuchKey", "404"):
                self.say(f"[s3] {job_name}: чтение метки завершения — {code}")
            return None
        try:
            return json.loads(body)
        except json.JSONDecodeError:
            return {"phase": "unparseable", "exit_code": -1}

    def collect(self, job_name: str, run_id: str | None = None, dest: Path | None = None,
                attempt: int = 1) -> list:
        """Скачать всё из префикса попытки (key_prefix) в dest (по умолчанию runs/<run_id>/<job>)."""
        run_id = run_id or self.cfg.run_id
        pfx = self.key_prefix(job_name, run_id, attempt) + "/"
        dest = Path(dest or self.cfg.out_dir / run_id / job_name)
        got, kw = [], {"Bucket": self.cfg.bucket, "Prefix": pfx}
        while True:
            r = self.s3.list_objects_v2(**kw)
            for obj in r.get("Contents", []):
                rel = obj["Key"][len(pfx):]
                parts = rel.split("/")
                if not rel or rel.startswith("/") or any(p in ("", ".", "..") for p in parts):
                    self.say(f"[сбор] пропущен ключ с недопустимым путём: {obj['Key']}")
                    continue
                target = dest.joinpath(*parts)
                target.parent.mkdir(parents=True, exist_ok=True)
                body = self.s3.get_object(Bucket=self.cfg.bucket, Key=obj["Key"])["Body"]
                with target.open("wb") as fh:
                    for chunk in iter(lambda: body.read(1 << 20), b""):
                        fh.write(chunk)
                got.append(rel)
            if not r.get("IsTruncated"):
                break
            kw["ContinuationToken"] = r["NextContinuationToken"]
        self.say(f"[сбор] {job_name} (попытка {attempt}): {len(got)} файл(ов) -> {dest}")
        return got

    # ---- гашение
    def _terminate(self, iids: list, why: str) -> set:
        """Погасить машины; вернуть множество id, гашение которых AWS принял.

        TerminateInstances отвергает весь пакет, если хоть один id не годится (например,
        InvalidInstanceID.NotFound у только что запущенной машины). Поэтому при отказе
        пакета — повтор по одной, чтобы один плохой id не оставил живыми остальные.
        Сетевой сбой (EndpointConnectionError, ReadTimeoutError — BotoCoreError) обрабатывается так же:
        иначе исключение ушло бы из finally мимо повтора по одной и команды ручного гашения.
        """
        from botocore.exceptions import BotoCoreError, ClientError

        def _code(e):
            return e.response.get("Error", {}).get("Code", "?") if isinstance(e, ClientError) else type(e).__name__
        iids = sorted(set(i for i in iids if i))
        if not iids:
            return set()
        self.say(f"[гашение] {iids} ({why})")
        try:
            self.ec2.terminate_instances(InstanceIds=iids)
            return set(iids)
        except (ClientError, BotoCoreError) as e:
            code = _code(e)
            if len(iids) > 1:
                self.say(f"[гашение] пакет отвергнут ({code}) — гашу по одной")
            else:
                self._manual(iids, code)
                return set()
        ok = set()
        for iid in iids:
            try:
                self.ec2.terminate_instances(InstanceIds=[iid])
                ok.add(iid)
            except (ClientError, BotoCoreError) as e:
                self._manual([iid], _code(e))
        return ok

    def _manual(self, iids: list, code: str) -> None:
        self.say(f"[гашение][ОШИБКА] {code} — погасите вручную: "
                 f"aws ec2 terminate-instances --region {self.cfg.region} --instance-ids {' '.join(iids)}")

    def _close(self, a: Active, outcome: str, state: str = "", reason: str = "") -> None:
        self.active.pop(a.iid, None)
        self._closed_n += 1
        self.ledger.append("end", run_id=self.cfg.run_id, job=a.job.name, iid=a.iid,
                           t_end=round(self.clock(), 1), outcome=outcome, state=state, reason=reason)

    def _drop(self, a: Active) -> None:
        """Гашение не подтверждено: снять с наблюдения, но запись в журнале оставить открытой.

        Открытая запись входит в бюджет верхней границей; закроет её сверка (reconcile)
        в начале следующего прогона или после `reap`, когда AWS покажет terminated.
        Машина, возможно, жива и занимает квоту: она остаётся в self._unconfirmed — очередь
        новых машин не запускает (иначе живых vCPU станет больше --max-vcpu), сборщик в finally
        гасит её вместе с остальными.
        """
        self.active.pop(a.iid, None)
        self._unconfirmed[a.iid] = a
        self.say(f"[журнал] {a.iid} ({a.job.name}): гашение не подтверждено — запись остаётся открытой, "
                 f"новые машины прогона не запускаются")

    def _retry_unconfirmed(self) -> None:
        """Ещё раз погасить машины с неподтверждённым гашением; принятые — закрыть в журнале."""
        ids = sorted(self._unconfirmed)
        done = self._terminate(ids, "повтор: гашение прежде не подтверждено")
        for iid in ids:
            if iid in done:
                self._close(self._unconfirmed.pop(iid), "reaped", reason="повторное гашение")

    def find_own(self, run_id: str | None = None) -> list:
        """Машины с тегами Project=<tag> и RunId=<run_id> в живых состояниях (с проверкой тегов на месте)."""
        run_id = run_id or self.cfg.run_id
        kw = {"Filters": [{"Name": "tag:Project", "Values": [self.cfg.tag]},
                          {"Name": "tag:RunId", "Values": [run_id]},
                          {"Name": "instance-state-name", "Values": LIVE_STATES}]}
        own = []
        while True:
            r = self.ec2.describe_instances(**kw)
            for res in r.get("Reservations", []):
                for i in res.get("Instances", []):
                    tags = {t["Key"]: t["Value"] for t in i.get("Tags", [])}
                    if tags.get("Project") == self.cfg.tag and tags.get("RunId") == run_id:
                        own.append(i["InstanceId"])
                    else:
                        self.say(f"[сборщик] {i['InstanceId']} — не из прогона {run_id}, не трогаю")
            if not r.get("NextToken"):
                return own
            kw["NextToken"] = r["NextToken"]

    def reap(self, why: str, run_id: str | None = None, known: list = (), list_only: bool = False) -> list:
        """Внешний сборщик: погасить машины этого прогона (по тегам) и известные нам id."""
        try:
            own = self.find_own(run_id)
        except Exception as e:     # не мешаем погасить хотя бы известные id
            self.say(f"[сборщик] поиск по тегам не удался: {type(e).__name__}")
            own = []
        targets = sorted(set(own) | set(known))
        if not targets:
            self.say(f"[сборщик] живых машин прогона {run_id or self.cfg.run_id} нет")
            return []
        if list_only:
            self.say(f"[сборщик] найдены (не гашу): {targets}")
            return targets
        done = self._terminate(targets, why)
        for iid in targets:
            if iid in self.active:
                if iid in done:
                    self._close(self.active[iid], "reaped", reason=why)
                else:
                    self._drop(self.active[iid])
            elif iid in self._unconfirmed and iid in done:
                self._close(self._unconfirmed.pop(iid), "reaped", reason=why)
        return targets

    def reconcile(self) -> None:
        """Сверка журнала с AWS по открытым записям (после падения сессии и т. п.).

        Погашена — закрыть верхней границей срока. Жива, но срок (max_min + запас) истёк — сторож и
        сборщик не сработали (user-data не стартовал, сессия упала): журнал доказывает, что машина наша, —
        гасим и закрываем по фактическому времени. Жива в пределах срока — «в работе».
        """
        self._reconciled = True
        open_ = self.ledger.summary()["open"]
        if not open_:
            return
        launches = {r["iid"]: r for r in self.ledger.records() if r.get("ev") == "launch"}
        info = self.describe(open_)
        now = self.clock()
        for iid in open_:
            st = info.get(iid, {}).get("state", "not-found")
            L = launches[iid]
            if st in ("shutting-down", "terminated", "not-found"):
                t_end = min(now, L["t_launch"] + L["allot_min"] * 60)   # верхняя граница
                self.ledger.append("end", run_id=L["run_id"], job=L["job"], iid=iid, t_end=round(t_end, 1),
                                   outcome="reconciled", state=st, reason="закрыто сверкой")
            elif now > L["t_launch"] + L["allot_min"] * 60:
                if iid in self._terminate([iid], f"сверка: {st} после срока {L['allot_min']} мин"):
                    self.ledger.append("end", run_id=L["run_id"], job=L["job"], iid=iid, t_end=round(now, 1),
                                       outcome="reaped_overdue", state=st, reason="жива после срока")
                else:
                    self.say(f"[сверка] {iid} ({L.get('job', '?')}) жива после срока и не погашена — в бюджете "
                             f"прогона считается по фактическому времени ({(now - L['t_launch']) / 3600:.1f} ч)")
            else:
                self.say(f"[сверка] {iid} ({L.get('job', '?')}) ещё {st} — учитывается как «в работе»")

    # ---- основной цикл
    def check_userdata(self, jobs: list) -> None:
        """Собрать user-data всех задач и всех попыток до первого запуска.

        Повтор пишет в префикс _attemptN — user-data длиннее; RunnerError (предел 16 КБ и т. п.)
        посреди пакета остановил бы прогон, и сборщик в finally погасил бы работающие машины.
        Зона и тип — самой длинной записи (длина зон региона одинакова).
        """
        cfg = self.cfg
        az = f"{cfg.region}a"
        for j in jobs:
            for att in range(1, cfg.retries + 2):
                build_userdata(j, run_id=cfg.run_id, attempt=att, repo=cfg.repo, branch=cfg.branch,
                               bucket=cfg.bucket, key_prefix=self.key_prefix(j.name, attempt=att),
                               region=cfg.region, upload_reserve_min=cfg.upload_reserve_min,
                               az=az, az_cheapest=az, itype=max(j.types, key=len))

    def run(self, jobs: list) -> dict:
        cfg = self.cfg
        mode = "plan" if self.offline else ("dryrun" if cfg.dryrun else "run")
        self.check_userdata(jobs)              # до первого запуска, а не посреди пакета
        self.journal.append("run_start", run_id=cfg.run_id, mode=mode, branch=cfg.branch,
                            jobs=[j.name for j in jobs], ceiling=cfg.ceiling, max_parallel=cfg.max_parallel,
                            max_vcpu=cfg.vcpu_quota)
        need = sum(j.vcpu for j in jobs)
        self.say(f"[прогон] {cfg.run_id} режим={mode} задач={len(jobs)} ветка={cfg.branch} "
                 f"параллельно≤{cfg.max_parallel} квота={cfg.vcpu_quota} vCPU потолок={cfg.ceiling} долл."
                 + (f" — очередь: {need} vCPU пакета > квоты, запуск по мере освобождения" if need > cfg.vcpu_quota else ""))
        for j in jobs:
            for w in job_warnings(j):
                self.say(f"[внимание] {w}")
        pending = deque((j, 1) for j in jobs)
        try:
            if mode == "run":
                self.reconcile()
            while pending or self.active:
                waiting = self._launch_some(pending, sideline=(mode != "run"))
                if not self.active:
                    if waiting:                  # квотная ошибка без своих машин в работе — пауза и повтор
                        self.sleep(cfg.poll_s)
                    continue          # без ожидания _launch_some всегда снимает задачу
                self.sleep(cfg.poll_s)
                self._poll(pending)
        finally:
            if mode == "run":
                self.reap("finally: конец прогона / ошибка / прерывание",
                          known=list(self.active) + list(self._unconfirmed))
            self.journal.append("run_end", run_id=cfg.run_id, mode=mode, results=self.results)
            summ = cfg.out_dir / cfg.run_id / ("summary.json" if mode == "run" else f"summary.{mode}.json")
            summ.parent.mkdir(parents=True, exist_ok=True)
            summ.write_text(self.redact(json.dumps({"run_id": cfg.run_id, "mode": mode, "branch": cfg.branch,
                                                    "results": self.results}, ensure_ascii=False, indent=1)),
                            encoding="utf-8")
        return self.results

    def _launch_some(self, pending: deque, sideline: bool = False) -> bool:
        """Запустить, что помещается в квоту, бюджет и --max-parallel.

        Возвращает True, если головная задача ждёт квоты при пустом списке своих машин (тогда
        цикл делает паузу poll_s); иначе False. Задача, упёршаяся в квоту или бюджет при своих
        машинах в работе, остаётся в очереди: квота освободится, а бюджет пересчитается, когда
        закроется машина (верхняя граница «в работе» сменится фактическим расходом).
        Нет ёмкости спота ни в одной зоне ни на одном типе — повтор через poll_s, не больше
        capacity_wait_polls раз (ёмкость обычно возвращается). Есть машина с неподтверждённым
        гашением — ещё одна попытка погасить, не вышло — очередь снимается со статусом halted.
        """
        cfg = self.cfg
        if self._unconfirmed and pending:
            self._retry_unconfirmed()
            if self._unconfirmed:
                why = (f"гашение {sorted(self._unconfirmed)} не подтверждено — машина, возможно, жива; "
                       f"новые машины прогона не запускаются")
                for job, attempt in pending:
                    self.results[job.name] = {"status": "halted", "attempts": attempt, "detail": why}
                self.say(f"[очередь] СТОП: {why}; снято задач: {len(pending)}")
                pending.clear()
                return False
        while pending:
            if len(self.active) >= cfg.max_parallel:
                return False
            job, attempt = pending[0]
            used = sum(a.job.vcpu for a in self.active.values())
            if used + job.vcpu > cfg.vcpu_quota:
                if self.active:
                    return False
                self.results[job.name] = {"status": "quota", "attempts": attempt,
                                          "detail": f"{job.vcpu} vCPU > квоты {cfg.vcpu_quota}"}
                pending.popleft()
                continue
            if self.active and self._budget_wait.get(job.name) == self._closed_n:
                return False                             # с прошлого отказа бюджет не менялся
            status, act = self.launch(job, attempt)
            if status == "launched":
                self.active[act.iid] = act
                self._budget_wait.pop(job.name, None)
                self._quota_wait.pop(job.name, None)
                self._cap_wait.pop(job.name, None)
                pending.popleft()
                continue
            if status == "defer":
                if self.active:
                    return False                         # ждём освобождения квоты
                if not sideline:
                    n = self._quota_wait[job.name] = self._quota_wait.get(job.name, 0) + 1
                    if n <= cfg.quota_wait_polls:
                        self.say(f"[очередь] {job.name}: квота занята (своих машин нет) — повтор через "
                                 f"{cfg.poll_s:.0f} с ({n}/{cfg.quota_wait_polls})")
                        return True
            if status == "capacity" and not sideline:
                n = self._cap_wait[job.name] = self._cap_wait.get(job.name, 0) + 1
                if n <= cfg.capacity_wait_polls:
                    self.say(f"[очередь] {job.name}: нет ёмкости спота ни в одной зоне ни на одном типе — "
                             f"повтор через {cfg.poll_s:.0f} с ({n}/{cfg.capacity_wait_polls})")
                    return not self.active               # без своих машин — пауза poll_s в цикле run
            if status == "refused" and not sideline and self.active:
                self._budget_wait[job.name] = self._closed_n
                self.say(f"[очередь] {job.name}: бюджет не пропускает — жду окончания машины прогона")
                return False
            pending.popleft()
            res = {"status": "quota" if status == "defer" else status, "attempts": attempt}
            if job.name in self._why:
                res["detail"] = self._why[job.name]
            self.results[job.name] = res
        return False

    def _res(self, a: Active) -> dict:
        return {"attempts": a.attempt, "itype": a.itype, "iid": a.iid, "az": a.az, "az_rank": a.az_rank,
                "az_cheapest": a.az_cheapest, "spot_usd_h": round(a.spot_usd_h, 4),
                "cheapest_zone": bool(a.az_cheapest) and a.az == a.az_cheapest}

    def _poll(self, pending: deque) -> None:
        info = self.describe(list(self.active))
        now = self.clock()
        for iid in list(self.active):
            a = self.active[iid]
            st = info.get(iid, {})
            state, reason = st.get("state", "pending?"), st.get("reason", "")
            if st.get("az") and not a.az:
                a.az = st["az"]
            done = self.marker(a.job.name, attempt=a.attempt)
            if done is not None and done.get("attempt") not in (None, a.attempt):
                self.say(f"[итог] {a.job.name}: метка попытки {done.get('attempt')} — не этой ({a.attempt}), не принимаю")
                done = None
            gone = state in ("shutting-down", "terminated")
            # Метка 143 (TERM) при ещё живой машине: причины гашения EC2 пока не показал — ждём до трёх
            # опросов, иначе прерванную часть сочли бы окончательным сбоем и не повторили.
            if (done is not None and done.get("exit_code") == 143 and not gone and now <= a.deadline
                    and (not a.t_mark or now - a.t_mark < 3 * self.cfg.poll_s)):
                a.t_mark = a.t_mark or now
                self.say(f"[ждём] {a.job.name}: метка с кодом 143 (TERM), машина ещё {state} — жду причину гашения")
                continue
            ok = done is not None and done.get("phase") == "done" and done.get("exit_code") == 0
            # Прерывание спота: машина гаснет не по нашей команде; ловушка TERM могла успеть выгрузить
            # метку с кодом 143 — такая метка не итог, а повод повторить часть в префикс новой попытки.
            # Причина Client.… (сторож shutdown -P, poweroff, наш terminate) — не прерывание: повтор
            # зависшей задачи только удвоил бы расход.
            spot = reason.startswith(SPOT_REASON_PREFIX)
            own = reason.startswith(OWN_REASON_PREFIX)
            interrupted = gone and not ok and (spot or (not own and (done or {}).get("exit_code") == 143))
            retry = a.attempt <= self.cfg.retries and (interrupted or (gone and done is None))
            res = self._res(a)
            if done is not None and not retry:
                files = self.collect(a.job.name, attempt=a.attempt)
                if iid in self._terminate([iid], "задача завершена"):
                    self._close(a, "done", state)
                else:
                    self._drop(a)
                res.update(status="ok" if ok else "failed", phase=done.get("phase"),
                           exit_code=done.get("exit_code"), files=len(files), az_machine=done.get("az", ""))
                if interrupted:
                    res["interrupted"] = True
                self.results[a.job.name] = res
                where = "самая дешёвая" if res["cheapest_zone"] else f"место по цене {a.az_rank + 1}"
                self.say(f"[итог] {a.job.name}: {res['status']} (фаза {res['phase']}, код {res['exit_code']}), "
                         f"зона {a.az} ({where}, {a.spot_usd_h:.4f} долл./ч)")
            elif gone:
                if retry:
                    lost_dir = self.cfg.out_dir / self.cfg.run_id / "_lost" / f"a{a.attempt}" / a.job.name
                    self.collect(a.job.name, attempt=a.attempt, dest=lost_dir)
                    self._close(a, "lost", state, reason)
                    if a.az:
                        self._avoid[(a.itype, a.az)] = now + self.cfg.avoid_s
                    pending.appendleft((a.job, a.attempt + 1))
                    self.say(f"[итог] {a.job.name}: машина {state} без итога ({reason or 'метка ' + str((done or {}).get('exit_code'))}) "
                             f"— повтор в префикс попытки {a.attempt + 1}, зона {a.az} — в конце очереди")
                else:
                    files = self.collect(a.job.name, attempt=a.attempt)
                    self._close(a, "lost", state, reason)
                    res.update(status="lost", reason=reason, files=len(files))
                    self.results[a.job.name] = res
                    self.say(f"[итог] {a.job.name}: машина {state} без метки завершения ({reason})")
            elif now > a.deadline:
                gone_ = iid in self._terminate([iid], "срок сборщика истёк")
                files = self.collect(a.job.name, attempt=a.attempt)
                if gone_:
                    self._close(a, "deadline", state)
                else:
                    self._drop(a)
                res.update(status="deadline", files=len(files))
                self.results[a.job.name] = res
            else:
                self.say(f"[ждём] {a.job.name}: {state}, {int((now - a.t_launch) / 60)} мин")


# --------------------------------------------------------------------------- CLI
def _clients(region: str):
    import boto3
    from botocore.config import Config as BCfg
    c = BCfg(region_name=region, retries={"max_attempts": 8, "mode": "adaptive"})
    return (boto3.client("ec2", config=c), boto3.client("s3", config=c),
            boto3.client("ssm", config=c), boto3.client("sts", config=c))


def _check_branch(repo: str, branch: str) -> None:
    try:
        r = subprocess.run(["git", "ls-remote", "--heads", repo, branch], capture_output=True,
                           text=True, timeout=60)
    except Exception as e:
        raise RunnerError(f"проверка ветки не удалась: {type(e).__name__}") from None
    # Точное совпадение ссылки: шаблон ls-remote «x» находит и refs/heads/xy/x.
    refs = {ln.split("\t", 1)[1].strip() for ln in r.stdout.splitlines() if "\t" in ln}
    if r.returncode != 0 or f"refs/heads/{branch}" not in refs:
        raise RunnerError(f"ветки {branch!r} нет в {repo} — машина не сможет её клонировать")


def _config_from_args(a, need_env: bool) -> Config:
    if need_env:
        missing = [k for k in ("AWS_ACCESS_KEY_ID", "AWS_SECRET_ACCESS_KEY", "ALPHA_O_S3",
                               "ALPHA_O_INSTANCE_PROFILE") if not os.environ.get(k)]
        if missing:
            raise RunnerError(f"не заданы переменные окружения: {missing} (значения не печатаются)")
        bucket, prefix = parse_s3_uri(os.environ["ALPHA_O_S3"])
        profile = os.environ["ALPHA_O_INSTANCE_PROFILE"]
    else:
        uri = os.environ.get("ALPHA_O_S3")
        bucket, prefix = parse_s3_uri(uri) if uri else ("example-bucket", "PREFIX")
        profile = os.environ.get("ALPHA_O_INSTANCE_PROFILE", "INSTANCE_PROFILE")
    run_id = getattr(a, "run_id", None) or new_run_id()
    if not RUN_ID_RE.match(run_id):
        raise RunnerError(f"недопустимый run_id {run_id!r}")
    cfg = Config(bucket=bucket, prefix=prefix, instance_profile=profile, run_id=run_id,
                 tag=os.environ.get("ALPHA_O_TAG", DEFAULT_TAG), region=a.region, ceiling=a.ceiling,
                 ledger_path=Path(a.ledger))
    for name in ("repo", "branch", "max_parallel", "vcpu_quota", "quota_wait_polls", "capacity_wait_polls", "price_margin",
                 "price_window_h", "ebs_usd_gb_month", "ipv4_usd_h", "grace_min", "poll_s", "retries",
                 "spent_before"):
        if getattr(a, name, None) is not None:
            setattr(cfg, name, getattr(a, name))
    if getattr(a, "out", None):
        cfg.out_dir = Path(a.out)
    validate_config(cfg)
    return cfg


def validate_config(cfg: Config) -> None:
    """Пороги, на которых держится бюджетный предохранитель и верхняя граница «в работе»."""
    bad = []
    if cfg.max_parallel < 1 or cfg.vcpu_quota < 1:
        bad.append("max_parallel, vcpu_quota (--max-vcpu) ≥ 1")
    if not cfg.ceiling > 0:
        bad.append("ceiling > 0")
    if not cfg.price_margin >= 1.0 or not cfg.price_window_h > 0:
        bad.append("price_margin ≥ 1.0, price_window_h > 0")
    if not (cfg.spent_before >= 0 and cfg.ebs_usd_gb_month >= 0 and cfg.ipv4_usd_h >= 0):
        bad.append("spent_before, ebs_usd_gb_month, ipv4_usd_h ≥ 0 (иначе потолок обходится)")
    if cfg.grace_min < 5:
        bad.append("grace_min ≥ 5 (запас на загрузку машины до старта сторожа)")
    if not cfg.poll_s >= 5 or cfg.retries < 0 or cfg.quota_wait_polls < 0 or cfg.capacity_wait_polls < 0:
        bad.append("poll_s ≥ 5, retries ≥ 0, quota_wait_polls ≥ 0, capacity_wait_polls ≥ 0")
    if bad:
        raise RunnerError("; ".join(bad))


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description="Прогоны ГТУ-платформы на спотовых машинах AWS (us-east-1).")
    sub = p.add_subparsers(dest="cmd", required=True)

    def common(sp):
        sp.add_argument("--ledger", default=str(HERE / "ledger.jsonl"), help="журнал расходов (JSONL)")
        sp.add_argument("--ceiling", type=float, default=2000.0, help="потолок расходов, долл. (2000, решение владельца 14.11.2026)")
        sp.add_argument("--region", default=REGION)

    def knobs(sp):
        sp.add_argument("jobs", help="JSON-файл задач")
        sp.add_argument("--branch", default=os.environ.get("GTU_BRANCH", "main"))
        sp.add_argument("--repo", default=DEFAULT_REPO)
        sp.add_argument("--run-id")
        sp.add_argument("--max-parallel", type=int, default=6, help="машин одновременно (6)")
        sp.add_argument("--max-vcpu", "--vcpu-quota", dest="vcpu_quota", type=int, default=384,
                        help="квота спота, vCPU (384): сверх неё задачи ждут в очереди")
        sp.add_argument("--quota-wait-polls", type=int, default=15,
                        help="квотная ошибка AWS без своих машин в работе: опросов ожидания (15)")
        sp.add_argument("--capacity-wait-polls", type=int, default=15,
                        help="нет ёмкости спота ни в одной зоне ни на одном типе: повторов через --poll-s (15)")
        sp.add_argument("--price-margin", type=float, default=1.2, help="MaxPrice = макс. цена за окно × запас")
        sp.add_argument("--price-window-h", type=float, default=6.0)
        sp.add_argument("--ebs-usd-gb-month", type=float, default=0.08)
        sp.add_argument("--ipv4-usd-h", type=float, default=0.005, help="публичный IPv4, долл./ч (допущение)")
        sp.add_argument("--grace-min", type=int, default=10, help="запас сборщика сверх сторожа, мин")
        sp.add_argument("--poll-s", type=float, default=60.0)
        sp.add_argument("--retries", type=int, default=1, help="повторов при потере машины (прерывание спота)")
        sp.add_argument("--spent-before", type=float, default=0.0, help="расход вне журнала, долл.")
        sp.add_argument("--max-job-min", type=int, default=720)
        sp.add_argument("--out", help="куда скачивать результаты (runs/)")
        common(sp)

    sp = sub.add_parser("run", help="по умолчанию проверка прав (DryRun); --live — живой запуск")
    knobs(sp)
    sp.add_argument("--dryrun", action="store_true", help="DryRun=True: проверка прав без запуска (и так по умолчанию)")
    sp.add_argument("--live", action="store_true", help="живой запуск спотовых машин (расходы!)")
    sp.add_argument("--skip-branch-check", action="store_true")
    sp = sub.add_parser("plan", help="без AWS: user-data, параметры run_instances, порядок зон, бюджет")
    knobs(sp)
    sp.add_argument("--assume-price", type=float, required=True, help="цена спота, долл./ч (допущение)")
    sp.add_argument("--assume-zone-prices", help="цены по зонам (допущение): us-east-1a=0.955,us-east-1c=0.706")
    sp.add_argument("--show-userdata", action="store_true")
    sp = sub.add_parser("status", help="журнал: израсходовано, в работе, остаток")
    common(sp)
    sp = sub.add_parser("reap", help="погасить машины прогона по тегам")
    sp.add_argument("--run-id", required=True)
    sp.add_argument("--list-only", action="store_true")
    common(sp)
    sp = sub.add_parser("collect", help="скачать результаты задачи")
    sp.add_argument("--run-id", required=True)
    sp.add_argument("--job", required=True)
    sp.add_argument("--attempt", type=int, default=1, help="номер попытки (повтор после прерывания — 2)")
    sp.add_argument("--out")
    common(sp)
    return p


def main(argv=None) -> int:
    a = build_parser().parse_args(argv)

    def _sig(signum, _frame):
        raise KeyboardInterrupt(f"сигнал {signum}")
    for s in (signal.SIGTERM, signal.SIGHUP):
        signal.signal(s, _sig)

    try:
        if a.cmd == "status":
            led = Ledger(a.ledger)
            s = led.summary()
            total = s["spent"] + s["committed"]
            print(f"журнал: {a.ledger}\nзакрыто машин: {s['closed']}, израсходовано ≈ {s['spent']:.2f} долл.\n"
                  f"открыто: {len(s['open'])} {s['open']}, верхняя граница по ним {s['committed']:.2f} долл.\n"
                  f"итого {total:.2f} из {a.ceiling:.2f}, остаток {a.ceiling - total:.2f} долл.")
            return 0
        if a.cmd == "plan":
            if not a.assume_price > 0:
                raise RunnerError("--assume-price: цена спота должна быть > 0 (иначе ни одна задача не планируется)")
            cfg = _config_from_args(a, need_env=False)
            cfg.offline_price = a.assume_price
            if a.assume_zone_prices:
                cfg.offline_zones = parse_zone_prices(a.assume_zone_prices)
            jobs = load_jobs(a.jobs, a.max_job_min)
            r = Runner(cfg)
            r.check_userdata(jobs)            # все попытки: повтор пишет в более длинный префикс
            for job in jobs:
                extra = {"python": f"Python {job.python} (uv)", "system": "системный python3",
                         "container": f"контейнер {job.container}, команд container_setup: {len(job.container_setup)}"}
                r.say(f"[задача] {job.name}: образ {job.image}, {extra[job.mode]}, {job.vcpu} vCPU, "
                      f"сторож {job.max_min} мин")
                for w in job_warnings(job):
                    r.say(f"[внимание] {w}")
                st, _ = r.launch(job)
                r.results[job.name] = {"status": st}
                if job.name in r._why:
                    r.results[job.name]["detail"] = r._why[job.name]
                if r.last_params and st == "planned":
                    shown = dict(r.last_params)
                    ud = shown.pop("UserData")
                    shown["UserData"] = f"<{len(ud.encode())} байт>"
                    r.say(json.dumps(shown, ensure_ascii=False, indent=1))
                    if a.show_userdata:
                        r.say(ud)
                r.last_params = None
            need = sum(j.vcpu for j in jobs)
            if need > cfg.vcpu_quota:
                r.say(f"[очередь] пакет {need} vCPU > квоты {cfg.vcpu_quota}: одновременно не больше "
                      f"{cfg.vcpu_quota} vCPU и {cfg.max_parallel} машин, остальные ждут освобождения")
            r.say(f"[план] прогноз пакета {r.virtual:.2f} долл. (цена спота — допущение {a.assume_price}/ч)")
            bad = {k: v for k, v in r.results.items() if v["status"] != "planned"}
            if bad:
                r.say(f"[план] НЕ спланированы: {bad}")
                return 1
            return 0
        if a.cmd == "run":
            # по умолчанию — проверка прав без запуска (DryRun=True); живой запуск только с --live или GTU_LIVE=1
            live = bool(getattr(a, "live", False)) or os.environ.get("GTU_LIVE") == "1"
            dry = (not live) or a.dryrun or os.environ.get("GTU_DRYRUN") == "1" or os.environ.get("ALPHA_O_DRYRUN") == "1"
            if dry and not a.dryrun:
                print("[умолчание] режим проверки прав (DryRun). Для живого запуска: --live или GTU_LIVE=1")
            cfg = _config_from_args(a, need_env=True)
            cfg.dryrun = dry
            jobs = load_jobs(a.jobs, a.max_job_min)
            if not a.skip_branch_check:
                _check_branch(cfg.repo, cfg.branch)
            ec2, s3, ssm, sts = _clients(cfg.region)
            arn = sts.get_caller_identity()["Arn"]
            print("вызов от имени:", re.sub(r"::\d{12}:", "::<аккаунт>:", arn))
            r = Runner(cfg, ec2, s3, ssm, sts)
            if dry:   # чтение префикса результатов — проверка s3:ListBucket
                from botocore.exceptions import ClientError
                try:
                    s3.list_objects_v2(Bucket=cfg.bucket, Prefix=f"{cfg.prefix}/gtu/", MaxKeys=1)
                    r.say("[dryrun] чтение префикса результатов: OK")
                except ClientError as e:
                    r.say(f"[dryrun] чтение префикса результатов: {e.response.get('Error', {}).get('Code', '?')}")
            res = r.run(jobs)
            print(json.dumps(res, ensure_ascii=False, indent=1))
            bad = [k for k, v in res.items() if v.get("status") not in ("ok", "dryrun_ok")]
            return 1 if bad else 0
        cfg = _config_from_args(a, need_env=True)
        ec2, s3, ssm, sts = _clients(cfg.region)
        r = Runner(cfg, ec2, s3, ssm, sts)
        if a.cmd == "reap":
            r.reap("ручной сборщик", run_id=a.run_id, list_only=a.list_only)
            if not a.list_only:
                r.reconcile()     # закрыть в журнале погашенные машины (верхней границей срока)
            return 0
        if a.cmd == "collect":
            r.collect(a.job, run_id=a.run_id, dest=Path(a.out) if a.out else None, attempt=a.attempt)
            return 0
    except RunnerError as e:
        print("ОТКАЗ:", e, file=sys.stderr)
        return 2
    except Exception as e:   # ошибки AWS: код и текст без бакета и ключей, машины уже погашены в finally
        if type(e).__module__.startswith("botocore"):
            red = Redactor(parse_s3_uri(os.environ["ALPHA_O_S3"])[0] if os.environ.get("ALPHA_O_S3") else None)
            print("ОШИБКА AWS:", red(f"{type(e).__name__}: {e}"), file=sys.stderr)
            return 3
        raise
    return 0


if __name__ == "__main__":
    sys.exit(main())
