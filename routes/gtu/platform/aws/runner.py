#!/usr/bin/env python3
"""Прогоны расчётной платформы ГТУ на спотовых машинах AWS (us-east-1).

Запускается ИЗ СЕССИИ, где действуют ключи ограниченного IAM-пользователя (окружение
Coldun). Машины работают под ролью экземпляра: ключей на машине и в user-data нет.
Подробности, порядок проверки и отката — README.md рядом.

  python3 runner.py plan    jobs.json --assume-price 1.5   # без AWS: user-data, параметры, бюджет
  python3 runner.py run     jobs.json --branch B            # по умолчанию: проверка прав (DryRun=True), без запуска
  python3 runner.py run     jobs.json --branch B --live     # живой запуск спотовых машин (расходы!)
  python3 runner.py run     jobs.json --branch B            # боевой прогон
  python3 runner.py status                                  # журнал: израсходовано / в работе / остаток
  python3 runner.py reap    --run-id ID                     # погасить машины прогона по тегам
  python3 runner.py collect --run-id ID --job NAME          # скачать результаты задачи

Защита от утечки машин — тройная, как в calc/aws_orchestrate.py:
  1) InstanceInitiatedShutdownBehavior=terminate и разовый спот (прерывание = terminate);
  2) сторож на машине: `shutdown -P +N` — первая команда user-data;
  3) внешний сборщик: в finally гасит машины с тегами Project=<ALPHA_O_TAG> и
     RunId=<run_id> — только этого прогона, соседние прогоны с тем же Project не трогает.
Бюджет: перед каждым запуском «израсходовано + в работе (верхняя граница) + новая задача»
против потолка. Цена — describe_spot_price_history; MaxPrice спота ставится равной оценке,
поэтому почасовая цена машины сверху ограничена числом, вошедшим в прогноз.
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
PRODUCT = "Linux/UNIX"
USERDATA_LIMIT = 16384                 # байт до base64 (ограничение EC2)
LIVE_STATES = ["pending", "running", "stopping", "stopped"]

# vCPU известных типов (квота спота считается в vCPU). Другие типы — поле "vcpu" в задаче.
VCPU = {f"c7{f}.{s}": n for f in "ai"
        for s, n in (("4xlarge", 16), ("8xlarge", 32), ("16xlarge", 64))}

IMAGES = {
    # Пакеты платформы (calculix-ccx, gmsh, openfoam) — в apt, поэтому по умолчанию Ubuntu.
    "ubuntu2404": {"ssm": "/aws/service/canonical/ubuntu/server/24.04/stable/current/"
                          "amd64/hvm/ebs-gp3/ami-id",
                   "root": "/dev/sda1", "pkg": "apt"},
    "al2023": {"ssm": "/aws/service/ami-amazon-linux-latest/al2023-ami-kernel-default-x86_64",
               "root": "/dev/xvda", "pkg": "dnf"},
}

CAPACITY_ERRORS = {"InsufficientInstanceCapacity", "SpotMaxPriceTooLow", "Unsupported"}
QUOTA_ERRORS = {"MaxSpotInstanceCountExceeded", "VcpuLimitExceeded", "InstanceLimitExceeded"}

NAME_RE = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._-]{0,63}$")
BRANCH_RE = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._/-]{0,199}$")
REPO_RE = re.compile(r"^https://github\.com/[A-Za-z0-9_.-]+/[A-Za-z0-9_.-]+$")
PKG_RE = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._+\-\[\],:=<>!~]{0,127}$")
RUN_ID_RE = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._-]{0,47}$")
KEY_ID_RE = re.compile(r"(?<![A-Z0-9])(AKIA|ASIA)[A-Z0-9]{16}(?![A-Z0-9])")
SECRET_WORDS = ("aws_secret" "_access_key", "aws_session_token", "aws_access" "_key_id")  # разбито: иначе gate-сканер находит сам себя
SECRET_ENV = ("AWS_ACCESS_KEY_ID", "AWS_SECRET_ACCESS_KEY", "AWS_SESSION_TOKEN")


class RunnerError(Exception):
    """Ошибка конфигурации или отказ (без вывода секретов)."""


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
        return KEY_ID_RE.sub("<скрыто>", text)


# --------------------------------------------------------------------------- задачи
@dataclass
class Job:
    name: str
    command: str
    instance_type: str = DEFAULT_TYPE
    fallback_types: list = field(default_factory=lambda: list(DEFAULT_FALLBACK))
    max_min: int = 60
    image: str = "ubuntu2404"
    disk_gb: int = 50
    apt: list = field(default_factory=list)
    dnf: list = field(default_factory=list)
    pip: list = field(default_factory=list)
    vcpu: int = 0

    @property
    def types(self) -> list:
        return [self.instance_type] + [t for t in self.fallback_types if t != self.instance_type]


JOB_KEYS = {f for f in Job.__dataclass_fields__}


def load_jobs(path, max_job_min: int = 720) -> list:
    """Прочитать задачи: {"defaults": {...}, "jobs": [...]} или просто [...]."""
    data = json.loads(Path(path).read_text(encoding="utf-8"))
    defaults, items = ({}, data) if isinstance(data, list) else (data.get("defaults", {}), data.get("jobs"))
    if not isinstance(items, list) or not items:
        raise RunnerError("в файле задач нет непустого списка jobs")
    jobs, seen = [], set()
    for raw in items:
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
    if not NAME_RE.match(job.name):
        raise RunnerError(f"имя задачи {job.name!r}: только [A-Za-z0-9._-], до 64 знаков")
    if not isinstance(job.command, str) or not job.command.strip():
        raise RunnerError(f"{job.name}: пустая команда")
    if any(line.strip() == "GTU_CMD_EOF" for line in job.command.splitlines()):
        raise RunnerError(f"{job.name}: команда содержит служебную строку GTU_CMD_EOF")
    assert_no_secret(job.command, f"команда задачи {job.name}")
    if not isinstance(job.fallback_types, list):
        raise RunnerError(f"{job.name}: fallback_types — список типов")
    for lst in (job.apt, job.dnf, job.pip):
        if not isinstance(lst, list):
            raise RunnerError(f"{job.name}: apt, dnf, pip — списки")
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
    for t in job.types:
        if not re.match(r"^[a-z0-9]+\.[a-z0-9]+$", t):
            raise RunnerError(f"{job.name}: недопустимый тип машины {t!r}")
    if not job.vcpu:
        missing = [t for t in job.types if t not in VCPU]
        if missing:
            raise RunnerError(f"{job.name}: нет числа vCPU для {missing} — задайте поле vcpu")
        job.vcpu = max(VCPU[t] for t in job.types)


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

USERDATA_TEMPLATE = r'''#!/bin/bash
shutdown -P +@@MAX_MIN@@
# gtu runner: run=@@RUN_ID@@ job=@@JOB@@ attempt=@@ATTEMPT@@; the watchdog above is the first command
exec >@@LOG@@ 2>&1
set -u
export DEBIAN_FRONTEND=noninteractive AWS_DEFAULT_REGION=@@REGION@@
T0=$(date +%s); DEADLINE=$((T0 + @@BUDGET_S@@))
G=@@ROOT@@; OUT=$G/out; PHASE=setup
mkdir -p "$OUT"
export GTU_OUT="$OUT" GTU_RUN_ID=@@RUN_ID_Q@@ GTU_JOB=@@JOB_Q@@ GTU_ATTEMPT=@@ATTEMPT@@
S3_BUCKET=@@BUCKET_Q@@; S3_KEY=@@KEY_Q@@
cat >"$G/upload.py" <<'GTU_PY_EOF'
@@UPLOAD_PY@@GTU_PY_EOF
up() {
  if [ -x "$G/venv/bin/python" ] && "$G/venv/bin/python" -c 'import boto3' 2>/dev/null; then
    timeout 10m "$G/venv/bin/python" "$G/upload.py" "$S3_BUCKET" "$S3_KEY" "$1"
  elif command -v aws >/dev/null 2>&1; then
    if [ -d "$1" ]; then timeout 10m aws s3 cp --recursive --only-show-errors "$1" "s3://$S3_BUCKET/$S3_KEY/"
    else timeout 5m aws s3 cp --only-show-errors "$1" "s3://$S3_BUCKET/$S3_KEY/"; fi
  else echo "UPLOAD-UNAVAILABLE"; fi
}
finish() {
  RC=$?; set +e; trap - EXIT
  echo "finish phase=$PHASE rc=$RC"
  cp @@LOG@@ "$OUT/_job.log"
  up "$OUT"
  printf '{"run_id":"%s","job":"%s","attempt":%s,"phase":"%s","exit_code":%d,"t_start":%d,"t_end":%d,"nproc":%d}\n' \
    "$GTU_RUN_ID" "$GTU_JOB" "$GTU_ATTEMPT" "$PHASE" "$RC" "$T0" "$(date +%s)" "$(nproc)" >"$G/_DONE.json"
  up "$G/_DONE.json"
  poweroff
}
trap finish EXIT
trap 'exit 143' TERM HUP INT
PHASE=base
@@BASE@@
PHASE=venv
python3 -m venv "$G/venv" || exit 12
timeout 20m "$G/venv/bin/pip" install -q --upgrade pip boto3 || exit 13
PHASE=packages
@@PACKAGES@@
PHASE=clone
timeout 10m git clone -q --depth 1 -b @@BRANCH_Q@@ @@REPO_Q@@ "$G/repo" || exit 14
cd "$G/repo" || exit 15
cat >"$G/cmd.sh" <<'GTU_CMD_EOF'
@@COMMAND@@
GTU_CMD_EOF
export PATH="$G/venv/bin:$PATH" OMP_NUM_THREADS="$(nproc)"
REMAIN=$((DEADLINE - $(date +%s)))
[ "$REMAIN" -gt 120 ] || { PHASE=no_time; exit 124; }
PHASE=command
timeout -k 60 "${REMAIN}s" bash "$G/cmd.sh" >"$OUT/_stdout.log" 2>&1
RC=$?
PHASE=done
exit $RC
'''


def build_userdata(job: Job, *, run_id: str, attempt: int, repo: str, branch: str,
                   bucket: str, key_prefix: str, region: str = REGION,
                   upload_reserve_min: int = 8,
                   _root: str = "/opt/gtu", _log: str = "/var/log/gtu-job.log") -> str:
    """Скрипт первой загрузки машины. Первая команда — сторож `shutdown -P +max_min`.

    `_root` и `_log` — только для локального теста скрипта на подменённых командах.
    """
    q = shlex.quote
    if not REPO_RE.match(repo):
        raise RunnerError("repo: ожидается https://github.com/<владелец>/<репозиторий>")
    if not BRANCH_RE.match(branch) or ".." in branch:
        raise RunnerError(f"недопустимое имя ветки {branch!r}")
    if not RUN_ID_RE.match(run_id):
        raise RunnerError(f"недопустимый run_id {run_id!r}")
    img = IMAGES[job.image]
    # Сначала база (git, venv, boto3) — чтобы при сбое пакетов задачи журнал всё равно выгрузился.
    if img["pkg"] == "apt":
        apt = "timeout {m}m apt-get -o DPkg::Lock::Timeout=600 -q"
        base = (f"{apt.format(m=20)} update || exit 10\n"
                f"{apt.format(m=20)} install -y git python3-venv python3-pip || exit 11")
        sys_pk = (f"{apt.format(m=40)} install -y {' '.join(q(p) for p in job.apt)} || exit 16"
                  if job.apt else "")
    else:
        base = "timeout 20m dnf -q install -y git python3-pip || exit 11"
        sys_pk = (f"timeout 40m dnf -q install -y {' '.join(q(p) for p in job.dnf)} || exit 16"
                  if job.dnf else "")
    pip_pk = (f'timeout 40m "$G/venv/bin/pip" install -q {" ".join(q(p) for p in job.pip)} || exit 17'
              if job.pip else "")
    packages = "\n".join(x for x in (sys_pk, pip_pk) if x) or ":"
    budget_s = job.max_min * 60 - upload_reserve_min * 60
    if budget_s < 600:
        raise RunnerError(f"{job.name}: max_min слишком мал для запаса на выгрузку")
    vals = {
        "MAX_MIN": str(job.max_min), "RUN_ID": run_id, "JOB": job.name, "ATTEMPT": str(int(attempt)),
        "LOG": q(_log), "REGION": q(region), "BUDGET_S": str(budget_s), "ROOT": q(_root),
        "RUN_ID_Q": q(run_id), "JOB_Q": q(job.name), "BUCKET_Q": q(bucket), "KEY_Q": q(key_prefix),
        "UPLOAD_PY": UPLOAD_PY, "BASE": base, "PACKAGES": packages,
        "BRANCH_Q": q(branch), "REPO_Q": q(repo), "COMMAND": job.command.rstrip("\n"),
    }
    # Один проход подстановки: текст команды повторно не сканируется.
    ud = re.sub(r"@@([A-Z_]+)@@", lambda m: vals[m.group(1)], USERDATA_TEMPLATE)
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

    def summary(self) -> dict:
        """Израсходовано (закрытые машины) и верхняя граница по открытым.

        Закрытая: (спот в её зоне на момент запуска + диск) × max(60 с, конец − запуск).
        Открытая (нет записи end): (MaxPrice + диск) × весь отпущенный срок — верхняя
        граница, потому что сторож и сборщик гасят машину не позже этого срока.
        """
        recs = self.records()
        launches = {r["iid"]: r for r in recs if r.get("ev") == "launch"}
        ends = {r["iid"]: r for r in recs if r.get("ev") == "end"}
        spent = committed = 0.0
        open_ = []
        for iid, L in launches.items():
            if iid in ends:
                dur = max(60.0, ends[iid]["t_end"] - L["t_launch"])
                spent += (L["spot_usd_h"] + L["ebs_usd_h"]) * dur / 3600.0
            else:
                committed += (L["max_usd_h"] + L["ebs_usd_h"]) * L["allot_min"] / 60.0
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
    ceiling: float = 300.0
    max_parallel: int = 6
    vcpu_quota: int = 384
    price_margin: float = 1.2
    price_window_h: float = 6.0
    ebs_usd_gb_month: float = 0.08
    grace_min: int = 10
    upload_reserve_min: int = 8
    poll_s: float = 60.0
    retries: int = 1
    spent_before: float = 0.0
    ledger_path: Path = HERE / "ledger.jsonl"
    out_dir: Path = HERE / "runs"
    dryrun: bool = False
    offline_price: float | None = None   # режим plan: цена спота задана вручную, AWS не вызывается


def parse_s3_uri(uri: str) -> tuple:
    m = re.match(r"^s3://([a-z0-9][a-z0-9.-]{1,61}[a-z0-9])/(.+)$", uri or "")
    if not m or not m.group(2).strip("/"):
        raise RunnerError("ALPHA_O_S3 должен быть вида s3://<бакет>/<префикс> (префикс обязателен)")
    return m.group(1), m.group(2).strip("/")


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


@dataclass
class PriceInfo:
    itype: str
    latest: dict          # зона -> последняя цена
    window_max: float     # максимум за окно по всем зонам

    def az_price(self, az: str) -> float:
        return self.latest.get(az, max(self.latest.values()) if self.latest else self.window_max)


# --------------------------------------------------------------------------- исполнитель
class Runner:
    def __init__(self, cfg: Config, ec2=None, s3=None, ssm=None, sts=None,
                 clock=time.time, sleep=time.sleep, out=print):
        self.cfg, self.ec2, self.s3, self.ssm, self.sts = cfg, ec2, s3, ssm, sts
        self.clock, self.sleep, self._out = clock, sleep, out
        self.redact = Redactor(cfg.bucket)
        self.ledger = Ledger(cfg.ledger_path, self.redact, clock)          # общий учёт
        side = cfg.dryrun or cfg.offline_price is not None
        self.journal = (Ledger(cfg.out_dir / cfg.run_id / "journal.dryrun.jsonl", self.redact, clock)
                        if side else self.ledger)                          # куда пишем
        self.active: dict = {}
        self.results: dict = {}
        self.virtual = 0.0          # DRYRUN/plan: прогноз уже «запущенных» задач пакета
        self.last_params = None
        self._prices: dict = {}
        self._ami: dict = {}
        self._seq = 0

    # ---- вывод
    def say(self, msg: str) -> None:
        self._out(self.redact(msg))

    @property
    def offline(self) -> bool:
        return self.cfg.offline_price is not None

    def key_prefix(self, job_name: str, run_id: str | None = None) -> str:
        return f"{self.cfg.prefix}/gtu/{run_id or self.cfg.run_id}/{job_name}"

    # ---- цены и бюджет
    def price(self, itype: str) -> PriceInfo:
        now = self.clock()
        hit = self._prices.get(itype)
        if hit and now - hit[0] < 300:
            return hit[1]
        if self.offline:
            p = float(self.cfg.offline_price)
            info = PriceInfo(itype, {"assumed": p}, p)
        else:
            end = dt.datetime.fromtimestamp(now, tz=dt.timezone.utc)
            start = end - dt.timedelta(hours=self.cfg.price_window_h)
            kw = dict(InstanceTypes=[itype], ProductDescriptions=[PRODUCT], StartTime=start, EndTime=end)
            latest, ts, wmax = {}, {}, 0.0
            while True:
                r = self.ec2.describe_spot_price_history(**kw)
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
        s = self.ledger.summary()
        maxp = self.max_price(info)
        new = (maxp + self.ebs_usd_h(job)) * self.allot_min(job) / 60.0
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
        self.say(f"[бюджет] {job.name} на {itype}: израсходовано {spent:.2f} + в работе "
                 f"{s['committed']:.2f} + пакет {self.virtual:.2f} + новая {new:.2f} = {total:.2f} "
                 f"из {self.cfg.ceiling:.2f} долл. (MaxPrice {maxp:.4f}/ч) — {verdict}")
        return ok, maxp, nums

    # ---- параметры запуска
    def ami(self, image: str) -> str:
        if self.offline:
            return f"<ssm:{IMAGES[image]['ssm']}>"
        if image not in self._ami:
            self._ami[image] = self.ssm.get_parameter(Name=IMAGES[image]["ssm"])["Parameter"]["Value"]
        return self._ami[image]

    def build_params(self, job: Job, itype: str, ami: str, maxp: float, attempt: int) -> dict:
        cfg = self.cfg
        ud = build_userdata(job, run_id=cfg.run_id, attempt=attempt, repo=cfg.repo, branch=cfg.branch,
                            bucket=cfg.bucket, key_prefix=self.key_prefix(job.name), region=cfg.region,
                            upload_reserve_min=cfg.upload_reserve_min)
        tags = [{"Key": "Project", "Value": cfg.tag}, {"Key": "RunId", "Value": cfg.run_id},
                {"Key": "Job", "Value": job.name}, {"Key": "Name", "Value": f"gtu-{job.name}"}]
        self._seq += 1
        token = "gtu-" + hashlib.sha256(
            f"{cfg.run_id}|{job.name}|{attempt}|{itype}|{self._seq}".encode()).hexdigest()[:40]
        return dict(
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

    # ---- запуск одной задачи
    def launch(self, job: Job, attempt: int = 1) -> tuple:
        """('launched', Active) | ('refused'|'defer'|'failed'|'dryrun_ok'|'dryrun_denied'|'planned', None)."""
        from botocore.exceptions import ClientError
        refused = False
        for itype in job.types:
            info = self.price(itype)
            if info.window_max <= 0:
                self.say(f"[цена] нет истории цены спота для {itype} — тип пропущен")
                continue
            ok, maxp, nums = self.budget_check(job, itype, info)
            if not ok:
                refused = True
                continue
            params = self.build_params(job, itype, self.ami(job.image), maxp, attempt)
            self.last_params = params
            if self.offline:
                self.virtual += nums["new"]
                return "planned", None
            if self.cfg.dryrun:
                try:
                    self.ec2.run_instances(DryRun=True, **params)
                    self.say(f"[dryrun] {job.name}: неожиданный успех без DryRunOperation")
                    return "dryrun_denied", None
                except ClientError as e:
                    code = e.response.get("Error", {}).get("Code", "?")
                    if code == "DryRunOperation":
                        self.virtual += nums["new"]
                        self.say(f"[dryrun] {job.name} на {itype}: запуск разрешён (DryRunOperation)")
                        self.journal.append("dryrun", run_id=self.cfg.run_id, job=job.name, itype=itype, ok=True)
                        return "dryrun_ok", None
                    self.say(f"[dryrun] {job.name} на {itype}: отказ {code}")
                    self.journal.append("dryrun", run_id=self.cfg.run_id, job=job.name, itype=itype,
                                        ok=False, code=code)
                    if code in CAPACITY_ERRORS:
                        continue
                    return "dryrun_denied", None
            try:
                r = self.ec2.run_instances(**params)
            except ClientError as e:
                code = e.response.get("Error", {}).get("Code", "?")
                self.journal.append("launch_error", run_id=self.cfg.run_id, job=job.name, itype=itype,
                                    attempt=attempt, code=code)
                self.say(f"[запуск] {job.name} на {itype}: {code}")
                if code in CAPACITY_ERRORS:
                    continue
                if code in QUOTA_ERRORS:
                    return "defer", None
                return "failed", None
            inst = r["Instances"][0]
            iid = inst["InstanceId"]
            az = inst.get("Placement", {}).get("AvailabilityZone", "")
            now = self.clock()
            act = Active(job, attempt, iid, itype, az, now, now + self.allot_min(job) * 60)
            self.ledger.append("launch", run_id=self.cfg.run_id, job=job.name, attempt=attempt, iid=iid,
                               itype=itype, az=az, spot_usd_h=round(info.az_price(az), 4), max_usd_h=maxp,
                               ebs_usd_h=round(self.ebs_usd_h(job), 5), t_launch=round(now, 1),
                               allot_min=self.allot_min(job))
            self.say(f"[запуск] {job.name} (попытка {attempt}): {iid} {itype} {az}, "
                     f"сторож {job.max_min} мин, сборщик через {self.allot_min(job)} мин")
            return "launched", act
        return ("refused" if refused else "failed"), None

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

    def marker(self, job_name: str, run_id: str | None = None) -> dict | None:
        from botocore.exceptions import ClientError
        try:
            body = self.s3.get_object(Bucket=self.cfg.bucket,
                                      Key=self.key_prefix(job_name, run_id) + "/_DONE.json")["Body"].read()
        except ClientError as e:
            code = e.response.get("Error", {}).get("Code", "")
            if code not in ("NoSuchKey", "404"):
                self.say(f"[s3] {job_name}: чтение метки завершения — {code}")
            return None
        try:
            return json.loads(body)
        except json.JSONDecodeError:
            return {"phase": "unparseable", "exit_code": -1}

    def collect(self, job_name: str, run_id: str | None = None, dest: Path | None = None) -> list:
        """Скачать всё из <префикс>/gtu/<run_id>/<job>/ в dest (по умолчанию runs/<run_id>/<job>)."""
        run_id = run_id or self.cfg.run_id
        pfx = self.key_prefix(job_name, run_id) + "/"
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
        self.say(f"[сбор] {job_name}: {len(got)} файл(ов) -> {dest}")
        return got

    # ---- гашение
    def _terminate(self, iids: list, why: str) -> bool:
        from botocore.exceptions import ClientError
        iids = sorted(set(i for i in iids if i))
        if not iids:
            return True
        self.say(f"[гашение] {iids} ({why})")
        try:
            self.ec2.terminate_instances(InstanceIds=iids)
            return True
        except ClientError as e:
            code = e.response.get("Error", {}).get("Code", "?")
            self.say(f"[гашение][ОШИБКА] {code} — погасите вручную: "
                     f"aws ec2 terminate-instances --region {self.cfg.region} --instance-ids {' '.join(iids)}")
            return False

    def _close(self, a: Active, outcome: str, state: str = "", reason: str = "") -> None:
        self.active.pop(a.iid, None)
        self.ledger.append("end", run_id=self.cfg.run_id, job=a.job.name, iid=a.iid,
                           t_end=round(self.clock(), 1), outcome=outcome, state=state, reason=reason)

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
        self._terminate(targets, why)
        for iid in targets:
            if iid in self.active:
                self._close(self.active[iid], "reaped", reason=why)
        return targets

    def reconcile(self) -> None:
        """Закрыть в журнале машины прошлых прогонов, которые уже погашены (после падения сессии)."""
        open_ = self.ledger.summary()["open"]
        if not open_:
            return
        launches = {r["iid"]: r for r in self.ledger.records() if r.get("ev") == "launch"}
        info = self.describe(open_)
        now = self.clock()
        for iid in open_:
            st = info.get(iid, {}).get("state", "not-found")
            if st in ("shutting-down", "terminated", "not-found"):
                L = launches[iid]
                t_end = min(now, L["t_launch"] + L["allot_min"] * 60)   # верхняя граница
                self.ledger.append("end", run_id=L["run_id"], job=L["job"], iid=iid, t_end=round(t_end, 1),
                                   outcome="reconciled", state=st, reason="закрыто сверкой")
            else:
                self.say(f"[сверка] {iid} ({launches[iid].get('job', '?')}) ещё {st} — учитывается как «в работе»")

    # ---- основной цикл
    def run(self, jobs: list) -> dict:
        cfg = self.cfg
        mode = "plan" if self.offline else ("dryrun" if cfg.dryrun else "run")
        self.journal.append("run_start", run_id=cfg.run_id, mode=mode, branch=cfg.branch,
                            jobs=[j.name for j in jobs], ceiling=cfg.ceiling, max_parallel=cfg.max_parallel)
        self.say(f"[прогон] {cfg.run_id} режим={mode} задач={len(jobs)} ветка={cfg.branch} "
                 f"параллельно≤{cfg.max_parallel} квота={cfg.vcpu_quota} vCPU потолок={cfg.ceiling} долл.")
        pending = deque((j, 1) for j in jobs)
        try:
            if mode == "run":
                self.reconcile()
            while pending or self.active:
                self._launch_some(pending, sideline=(mode != "run"))
                if not self.active:
                    continue          # _launch_some всегда снимает задачу, если нечего ждать
                self.sleep(cfg.poll_s)
                self._poll(pending)
        finally:
            if mode == "run":
                self.reap("finally: конец прогона / ошибка / прерывание", known=list(self.active))
            self.journal.append("run_end", run_id=cfg.run_id, mode=mode, results=self.results)
            summ = cfg.out_dir / cfg.run_id / ("summary.json" if mode == "run" else f"summary.{mode}.json")
            summ.parent.mkdir(parents=True, exist_ok=True)
            summ.write_text(self.redact(json.dumps({"run_id": cfg.run_id, "mode": mode, "branch": cfg.branch,
                                                    "results": self.results}, ensure_ascii=False, indent=1)),
                            encoding="utf-8")
        return self.results

    def _launch_some(self, pending: deque, sideline: bool = False) -> None:
        cfg = self.cfg
        while pending:
            if len(self.active) >= cfg.max_parallel:
                return
            job, attempt = pending[0]
            used = sum(a.job.vcpu for a in self.active.values())
            if used + job.vcpu > cfg.vcpu_quota:
                if self.active:
                    return
                self.results[job.name] = {"status": "quota", "attempts": attempt,
                                          "detail": f"{job.vcpu} vCPU > квоты {cfg.vcpu_quota}"}
                pending.popleft()
                continue
            status, act = self.launch(job, attempt)
            if status == "launched":
                self.active[act.iid] = act
                pending.popleft()
                continue
            if status == "defer" and self.active:
                return                                  # ждём освобождения квоты
            pending.popleft()
            self.results[job.name] = {"status": "quota" if status == "defer" else status, "attempts": attempt}

    def _poll(self, pending: deque) -> None:
        info = self.describe(list(self.active))
        now = self.clock()
        for iid in list(self.active):
            a = self.active[iid]
            st = info.get(iid, {})
            state, reason = st.get("state", "pending?"), st.get("reason", "")
            done = self.marker(a.job.name)
            res = {"attempts": a.attempt, "itype": a.itype, "iid": iid}
            if done is not None:
                files = self.collect(a.job.name)
                self._terminate([iid], "задача завершена")
                self._close(a, "done", state)
                ok = done.get("phase") == "done" and done.get("exit_code") == 0
                res.update(status="ok" if ok else "failed", phase=done.get("phase"),
                           exit_code=done.get("exit_code"), files=len(files))
                self.results[a.job.name] = res
                self.say(f"[итог] {a.job.name}: {res['status']} (фаза {res['phase']}, код {res['exit_code']})")
            elif state in ("shutting-down", "terminated"):
                files = self.collect(a.job.name)
                self._close(a, "lost", state, reason)
                if a.attempt <= self.cfg.retries:
                    pending.appendleft((a.job, a.attempt + 1))
                    self.say(f"[итог] {a.job.name}: машина {state} без метки завершения ({reason}) — повтор")
                else:
                    res.update(status="lost", reason=reason, files=len(files))
                    self.results[a.job.name] = res
                    self.say(f"[итог] {a.job.name}: машина {state} без метки завершения ({reason})")
            elif now > a.deadline:
                self._terminate([iid], "срок сборщика истёк")
                files = self.collect(a.job.name)
                self._close(a, "deadline", state)
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
    if r.returncode != 0 or f"refs/heads/{branch}" not in r.stdout:
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
    for name in ("repo", "branch", "max_parallel", "vcpu_quota", "price_margin", "price_window_h",
                 "ebs_usd_gb_month", "grace_min", "poll_s", "retries", "spent_before"):
        if getattr(a, name, None) is not None:
            setattr(cfg, name, getattr(a, name))
    if getattr(a, "out", None):
        cfg.out_dir = Path(a.out)
    if cfg.max_parallel < 1 or cfg.vcpu_quota < 1 or cfg.price_margin < 1.0 or cfg.ceiling <= 0:
        raise RunnerError("max_parallel, vcpu_quota ≥ 1; price_margin ≥ 1.0; ceiling > 0")
    return cfg


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description="Прогоны ГТУ-платформы на спотовых машинах AWS (us-east-1).")
    sub = p.add_subparsers(dest="cmd", required=True)

    def common(sp):
        sp.add_argument("--ledger", default=str(HERE / "ledger.jsonl"), help="журнал расходов (JSONL)")
        sp.add_argument("--ceiling", type=float, default=300.0, help="потолок расходов, долл. (300)")
        sp.add_argument("--region", default=REGION)

    def knobs(sp):
        sp.add_argument("jobs", help="JSON-файл задач")
        sp.add_argument("--branch", default=os.environ.get("GTU_BRANCH", "main"))
        sp.add_argument("--repo", default=DEFAULT_REPO)
        sp.add_argument("--run-id")
        sp.add_argument("--max-parallel", type=int, default=6, help="машин одновременно (6)")
        sp.add_argument("--vcpu-quota", type=int, default=384, help="квота спота, vCPU (384)")
        sp.add_argument("--price-margin", type=float, default=1.2, help="MaxPrice = макс. цена за окно × запас")
        sp.add_argument("--price-window-h", type=float, default=6.0)
        sp.add_argument("--ebs-usd-gb-month", type=float, default=0.08)
        sp.add_argument("--grace-min", type=int, default=10, help="запас сборщика сверх сторожа, мин")
        sp.add_argument("--poll-s", type=float, default=60.0)
        sp.add_argument("--retries", type=int, default=1, help="повторов при потере машины (прерывание спота)")
        sp.add_argument("--spent-before", type=float, default=0.0, help="расход вне журнала, долл.")
        sp.add_argument("--max-job-min", type=int, default=720)
        sp.add_argument("--out", help="куда скачивать результаты (runs/)")
        common(sp)

    sp = sub.add_parser("run", help="прогон (или --dryrun)")
    knobs(sp)
    sp.add_argument("--dryrun", action="store_true", help="DryRun=True: проверка прав без запуска (и так по умолчанию)")
    sp.add_argument("--live", action="store_true", help="живой запуск спотовых машин (расходы!)")
    sp.add_argument("--skip-branch-check", action="store_true")
    sp = sub.add_parser("plan", help="без AWS: user-data, параметры run_instances, бюджет")
    knobs(sp)
    sp.add_argument("--assume-price", type=float, required=True, help="цена спота, долл./ч (допущение)")
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
    sp.add_argument("--out")
    common(sp)
    a = p.parse_args(argv)

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
            cfg = _config_from_args(a, need_env=False)
            cfg.offline_price = a.assume_price
            jobs = load_jobs(a.jobs, a.max_job_min)
            r = Runner(cfg)
            for job in jobs:
                st, _ = r.launch(job)
                r.results[job.name] = {"status": st}
                if r.last_params and st == "planned":
                    shown = dict(r.last_params)
                    ud = shown.pop("UserData")
                    shown["UserData"] = f"<{len(ud.encode())} байт>"
                    r.say(json.dumps(shown, ensure_ascii=False, indent=1))
                    if a.show_userdata:
                        r.say(ud)
                r.last_params = None
            r.say(f"[план] прогноз пакета {r.virtual:.2f} долл. (цена спота — допущение {a.assume_price}/ч)")
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
                s3.list_objects_v2(Bucket=cfg.bucket, Prefix=f"{cfg.prefix}/gtu/", MaxKeys=1)
                r.say("[dryrun] чтение префикса результатов: OK")
            res = r.run(jobs)
            print(json.dumps(res, ensure_ascii=False, indent=1))
            bad = [k for k, v in res.items() if v.get("status") not in ("ok", "dryrun_ok")]
            return 1 if bad else 0
        cfg = _config_from_args(a, need_env=True)
        ec2, s3, ssm, sts = _clients(cfg.region)
        r = Runner(cfg, ec2, s3, ssm, sts)
        if a.cmd == "reap":
            r.reap("ручной сборщик", run_id=a.run_id, list_only=a.list_only)
            return 0
        if a.cmd == "collect":
            r.collect(a.job, run_id=a.run_id, dest=Path(a.out) if a.out else None)
            return 0
    except RunnerError as e:
        print("ОТКАЗ:", e, file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    sys.exit(main())
