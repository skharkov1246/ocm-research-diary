#!/usr/bin/env python3
"""Шаблоны AWS Batch трека ГТУ -> готовые запросы, план, job.json, manifest.json. AWS не вызывается.

    python3 render.py [--config batch.toml] [--region us-east-1|eu-north-1] [--mode line|point] [--out DIR]

Пишет в <out> (по умолчанию runs/<run_id>/<region>/, в .gitignore):
    compute_environment.json, job_queue.json, job_definition.json, submit_array_job.json — тела запросов
        CreateComputeEnvironment, CreateJobQueue, RegisterJobDefinition, SubmitJob (для aws --cli-input-json);
    s3/inputs/plan.csv, s3/inputs/job.json, s3/manifest.json — содержимое корня прогона в S3;
    commands.sh — команды для владельца (создать, отправить, разобрать); render.py их НЕ выполняет.

Сверка: каждый запрос проверяется по модели API AWS Batch из botocore (то же описание, из которого собран
справочник API: обязательные поля, типы, допустимые значения перечислений) и по ограничениям из текста
документации (см. CHECKS и README). Ошибка — код 1, файлы не пишутся.
"""
from __future__ import annotations

import argparse
import copy
import csv
import datetime as dt
import hashlib
import io
import json
import os
import pathlib
import re
import subprocess
import sys
import tomllib

HERE = pathlib.Path(__file__).resolve().parent
TEMPLATES = {"compute_environment": "CreateComputeEnvironment", "job_queue": "CreateJobQueue",
             "job_definition": "RegisterJobDefinition", "submit_array_job": "SubmitJob"}
REGIONS = ("us-east-1", "eu-north-1")
NAME_RE = re.compile(r"^[A-Za-z0-9_-]{1,128}$")
RUN_ID_RE = re.compile(r"^[A-Za-z0-9._-]{1,48}$")
PLACEHOLDER_RE = re.compile(r"\$\{([A-Z0-9_]+)\}")


# ---------------------------------------------------------------- конфигурация
def load_config(path: pathlib.Path) -> dict:
    with open(path, "rb") as f:
        return tomllib.load(f)


def load_rates(path: pathlib.Path = HERE / "aws_rates.toml") -> dict:
    with open(path, "rb") as f:
        return tomllib.load(f)


def run_prefix(cfg: dict) -> str:
    a = cfg["aws"]
    parts = [a["bucket"], a.get("s3_prefix", "").strip("/"), "gtu", cfg["run"]["run_id"]]
    return "s3://" + "/".join(p for p in parts if p)


def names(cfg: dict) -> dict:
    rid = cfg["run"]["run_id"].replace(".", "-")      # в именах Batch точка недопустима
    return {"CE_NAME": f"gtu-ce-{rid}", "JQ_NAME": f"gtu-jq-{rid}", "JD_NAME": f"gtu-jd-{rid}",
            "JOB_NAME": f"gtu-{rid}"}


# ---------------------------------------------------------------- план
def build_plan(cfg: dict, mode: str) -> list[dict]:
    p = cfg["plan"]
    n_pts = int(p["points_per_line"])
    rows, idx = [], 0
    scale_next = round(p["iters_next"] / p["iters_first"], 4)
    for si, stage in enumerate(p["stages"]):
        for ni, frac in enumerate(p["speed_frac"]):
            line = si * len(p["speed_frac"]) + ni
            for k in range(n_pts):
                pf = p["p_out_frac_first"] + (p["p_out_frac_last"] - p["p_out_frac_first"]) * k / max(n_pts - 1, 1)
                pf = round(pf, 4)
                p_pa = round(pf * p["p_out_design_pa"][si], 1)
                rows.append({
                    "index": idx, "row_id": f"{stage}-n{frac:.2f}-p{k + 1}", "line": line, "stage": stage,
                    "case": p["case_by_stage"][si], "level": p["level"], "plan": p["plan"],
                    "speed_frac": f"{frac:.2f}", "point": k + 1, "p_out_frac": f"{pf:.4f}", "p_out_pa": f"{p_pa:.1f}",
                    "iter_scale": "1.0" if (k == 0 or mode == "point") else f"{scale_next}",
                    "state": f"{frac:.2f},{p_pa:.1f}", "example": "true",
                })
                idx += 1
    return rows


def plan_csv(rows: list[dict]) -> str:
    buf = io.StringIO()
    w = csv.DictWriter(buf, fieldnames=list(rows[0].keys()), lineterminator="\n")
    w.writeheader()
    w.writerows(rows)
    return buf.getvalue()


def array_size(rows: list[dict], mode: str) -> int:
    return len({r["line"] for r in rows}) if mode == "line" else len(rows)


# ---------------------------------------------------------------- подстановка
def variables(cfg: dict, region: str, mode: str, rows: list[dict]) -> dict:
    a, c, j, r = cfg["aws"], cfg["compute"], cfg["job"], cfg["run"]
    nm = names(cfg)
    prefix = run_prefix(cfg)
    v = dict(nm)
    v.update({
        "RUN_ID": r["run_id"], "PROJECT_TAG": r["project_tag"], "REGION": region,
        "MINV_CPUS": int(c["minv_cpus"]), "MAXV_CPUS": int(c["maxv_cpus"]), "INSTANCE_TYPES": list(c["instance_types"]),
        "SUBNETS": list(a["subnets"]), "SECURITY_GROUP_IDS": list(a["security_group_ids"]),
        "INSTANCE_ROLE": a["instance_role"], "SPOT_FLEET_ROLE_ARN": a["spot_fleet_role_arn"],
        "CE_ARN": f"arn:aws:batch:{region}:{a['account_id']}:compute-environment/{nm['CE_NAME']}",
        "RUNNABLE_TIMEOUT_S": int(j["runnable_timeout_s"]),
        "IMAGE_URI": f"{a['account_id']}.dkr.ecr.{region}.amazonaws.com/{a['ecr_repository']}:{a['image_tag']}",
        "JOB_ROLE_ARN": a["job_role_arn"], "VCPU_STR": str(int(j["vcpu"])), "MEMORY_MIB_STR": str(int(j["memory_mib"])),
        "SHARED_MEMORY_MIB": int(j["shared_memory_mib"]), "LOG_STREAM_PREFIX": j["log_stream_prefix"],
        "RETRY_ATTEMPTS": int(j["retry_attempts"]), "ATTEMPT_DURATION_S": int(j["attempt_duration_s"]),
        "JOB_URI": f"{prefix}/inputs/job.json", "ARRAY_SIZE": array_size(rows, mode),
    })
    return v


def substitute(node, v: dict):
    """Строка ровно "${X}" -> значение X своего типа; иначе — подстановка в текст. Неизвестная — KeyError."""
    if isinstance(node, dict):
        return {k: substitute(x, v) for k, x in node.items()}
    if isinstance(node, list):
        return [substitute(x, v) for x in node]
    if isinstance(node, str):
        m = PLACEHOLDER_RE.fullmatch(node)
        if m:
            return copy.deepcopy(v[m.group(1)])
        return PLACEHOLDER_RE.sub(lambda mm: str(v[mm.group(1)]), node)
    return node


def render_templates(v: dict) -> dict:
    out = {}
    for name in TEMPLATES:
        tpl = json.loads((HERE / "templates" / f"{name}.json").read_text(encoding="utf-8"))
        out[name] = substitute(tpl, v)
    return out


# ---------------------------------------------------------------- сверка с моделью API
def batch_model():
    import botocore.session  # noqa: PLC0415 — модель API читается из пакета, без сети и ключей
    return botocore.session.get_session().get_service_model("batch")


def _enum_errors(value, shape, path: str) -> list[str]:
    errs = []
    t = shape.type_name
    if t == "structure" and isinstance(value, dict):
        for k, x in value.items():
            if k in shape.members:
                errs += _enum_errors(x, shape.members[k], f"{path}.{k}")
    elif t == "list" and isinstance(value, list):
        for i, x in enumerate(value):
            errs += _enum_errors(x, shape.member, f"{path}[{i}]")
    elif t == "map" and isinstance(value, dict):
        for k, x in value.items():
            errs += _enum_errors(x, shape.value, f"{path}[{k}]")
    elif t == "string" and getattr(shape, "enum", None):
        if value not in shape.enum:
            errs.append(f"{path}: {value!r} нет среди {shape.enum}")
    return errs


def validate_api(requests: dict, model=None) -> list[str]:
    """ParamValidator botocore (обязательные поля, типы, лишние поля) + допустимые значения перечислений."""
    from botocore.validate import ParamValidator  # noqa: PLC0415
    model = model or batch_model()
    errs = []
    for name, op in TEMPLATES.items():
        shape = model.operation_model(op).input_shape
        rep = ParamValidator().validate(requests[name], shape)
        if rep.has_errors():
            errs.append(f"{name} ({op}): {rep.generate_report()}")
        errs += [f"{name} ({op}): {e}" for e in _enum_errors(requests[name], shape, name)]
    return errs


def validate_limits(req: dict, cfg: dict, rates: dict, region: str) -> list[str]:
    """Ограничения из текста документации AWS Batch (botocore 1.44.3, архив руководства) и решения владельца."""
    try:
        return _limits(req, cfg, rates, region)
    except (KeyError, TypeError, StopIteration) as e:   # нет поля — его уже назвала сверка с моделью API
        return [f"ограничения не проверены: нет поля {e}"]


def _limits(req: dict, cfg: dict, rates: dict, region: str) -> list[str]:
    errs = []
    ce, jq, jd, sj = (req[k] for k in TEMPLATES)
    cr = ce["computeResources"]
    quota, maxinst = cfg["quota"]["vcpu_quota"], cfg["compute"]["max_instance_vcpu"]
    ec2 = rates["regions"][region]["ec2"]
    for n, val in (("computeEnvironmentName", ce["computeEnvironmentName"]), ("jobQueueName", jq["jobQueueName"]),
                   ("jobDefinitionName", jd["jobDefinitionName"]), ("jobName", sj["jobName"])):
        if not NAME_RE.match(val):
            errs.append(f"{n}: до 128 знаков, буквы, цифры, - и _ ({val!r})")
    if not RUN_ID_RE.match(cfg["run"]["run_id"]):
        errs.append("run_id: [A-Za-z0-9._-], до 48 знаков")
    if region not in REGIONS:
        errs.append(f"регион {region}: только {REGIONS} (решение владельца 14.11.2026)")
    if cr["minvCpus"] != 0:
        errs.append("minvCpus должен быть 0 (среда сжимается до нуля)")
    if cr["maxvCpus"] + maxinst > quota:
        errs.append(f"maxvCpus {cr['maxvCpus']} + одна машина {maxinst} > квоты {quota} "
                    "(Batch может превысить maxvCpus не более чем на одну машину)")
    for t in cr["instanceTypes"]:
        if t not in ec2:
            errs.append(f"тип {t} не найден в тарифах {region} (aws_rates.toml)")
        elif ec2[t]["vcpu"] != maxinst:
            errs.append(f"тип {t}: {ec2[t]['vcpu']} vCPU, ожидалось {maxinst} (задание — на всю машину)")
    vcpu = int(next(r["value"] for r in jd["containerProperties"]["resourceRequirements"] if r["type"] == "VCPU"))
    mem = int(next(r["value"] for r in jd["containerProperties"]["resourceRequirements"] if r["type"] == "MEMORY"))
    if mem < 4:
        errs.append("MEMORY: не меньше 4 МиБ")
    for t in cr["instanceTypes"]:
        if t in ec2 and (vcpu > ec2[t]["vcpu"] or mem > 0.9 * ec2[t]["memory_gib"] * 1024):
            errs.append(f"задание {vcpu} vCPU / {mem} МиБ не помещается на {t} (память — с запасом 10 %)")
    rs = jd["retryStrategy"]
    if not 1 <= rs["attempts"] <= 10:
        errs.append("retryStrategy.attempts: от 1 до 10")
    if len(rs.get("evaluateOnExit", [])) > 5:
        errs.append("evaluateOnExit: не больше 5 условий")
    if jd["timeout"]["attemptDurationSeconds"] < 60 or sj["timeout"]["attemptDurationSeconds"] < 60:
        errs.append("attemptDurationSeconds: не меньше 60 с")
    for a in jq.get("jobStateTimeLimitActions", []):
        if not 600 <= a["maxTimeSeconds"] <= 86400:
            errs.append("jobStateTimeLimitActions.maxTimeSeconds: от 600 до 86 400 с")
    if len(jq["computeEnvironmentOrder"]) > 3:
        errs.append("computeEnvironmentOrder: не больше трёх сред")
    if not 2 <= sj["arrayProperties"]["size"] <= 10000:
        errs.append("arrayProperties.size: от 2 до 10 000")
    for e in jd["containerProperties"].get("environment", []):
        if e["name"].startswith("AWS_BATCH"):
            errs.append(f"переменная {e['name']}: префикс AWS_BATCH зарезервирован")
    if cfg["job"]["internal_margin_s"] >= cfg["job"]["attempt_duration_s"]:
        errs.append("internal_margin_s ≥ attempt_duration_s")
    return errs


# ---------------------------------------------------------------- manifest и команды
def git_info(repo: pathlib.Path = HERE) -> dict:
    def g(*args):
        try:
            return subprocess.run(["git", *args], cwd=repo, capture_output=True, text=True, timeout=20).stdout.strip()
        except (OSError, subprocess.SubprocessError):
            return ""
    sha = g("rev-parse", "HEAD") or "unknown"
    return {"sha": sha, "branch": g("rev-parse", "--abbrev-ref", "HEAD"), "dirty": bool(g("status", "--porcelain"))}


def manifest(cfg: dict, region: str, mode: str, req: dict, plan_text: str, rows: list[dict], git: dict) -> dict:
    j = cfg["job"]
    return {
        "schema": "gtu-run-manifest/1",
        "run_id": cfg["run"]["run_id"],
        "created_utc": dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds"),
        "region": region,
        "run_prefix": run_prefix(cfg),
        "code": {"git_sha": git["sha"], "git_branch": git["branch"], "git_dirty": git["dirty"],
                 "image_uri": req["job_definition"]["containerProperties"]["image"],
                 "image_digest": "PLACEHOLDER: sha256 образа после docker push (docker inspect RepoDigests)"},
        "solver": dict(cfg["solver"]),
        "plan": {"file": "inputs/plan.csv", "sha256": hashlib.sha256(plan_text.encode()).hexdigest(),
                 "rows": len(rows), "mode": mode, "array_size": req["submit_array_job"]["arrayProperties"]["size"],
                 "example": True},
        "batch": {"compute_environment": req["compute_environment"]["computeEnvironmentName"],
                  "job_queue": req["job_queue"]["jobQueueName"],
                  "job_definition": req["job_definition"]["jobDefinitionName"],
                  "attempt_duration_s": j["attempt_duration_s"], "retry_attempts": j["retry_attempts"],
                  "vcpu": j["vcpu"], "memory_mib": j["memory_mib"],
                  "instance_types": cfg["compute"]["instance_types"], "maxv_cpus": cfg["compute"]["maxv_cpus"]},
        "layout": {"inputs/": "план, job.json, сетки и конфигурации — читает каждое задание",
                   "outputs/<row_id>/attempt-<n>/": "всё из $GTU_OUT точки; _status.json — последним",
                   "logs/<row_id>/attempt-<n>.log": "вывод команды точки (тот же текст — в CloudWatch Logs)",
                   "manifest.json": "этот файл: прогон, код, образ, решатель, план"},
    }


def commands_sh(cfg: dict, region: str, v: dict) -> str:
    prefix = run_prefix(cfg)
    return f"""#!/usr/bin/env bash
# Команды для владельца (админ-доступ, окружение Coldun). render.py их НЕ выполняет.
# Не запускать без решения владельца: создание среды Batch, отправка заданий и push образа — расходы.
# Регион: {region}. Прогон: {cfg['run']['run_id']}.
set -euo pipefail
R={region}

# 0. Образ (один раз на коммит): сборка в чистой копии коммита, push в ECR
#   aws ecr create-repository --region $R --repository-name {cfg['aws']['ecr_repository']}   # если нет
#   aws ecr get-login-password --region $R | docker login --username AWS --password-stdin {cfg['aws']['account_id']}.dkr.ecr.$R.amazonaws.com
#   docker tag gtu-platform:<sha12> {v['IMAGE_URI']} && docker push {v['IMAGE_URI']}

# 1. Входы и паспорт прогона в S3
aws s3 cp --region $R s3/ {prefix}/ --recursive

# 2. Среда, очередь, определение задания
aws batch create-compute-environment --region $R --cli-input-json file://compute_environment.json
aws batch describe-compute-environments --region $R --compute-environments {v['CE_NAME']} \\
    --query 'computeEnvironments[0].status'          # ждать VALID
aws batch create-job-queue --region $R --cli-input-json file://job_queue.json
aws batch register-job-definition --region $R --cli-input-json file://job_definition.json

# 3. Серия — массив заданий (индекс -> строки плана)
aws batch submit-job --region $R --cli-input-json file://submit_array_job.json

# 4. Ход и итог
#   aws batch describe-jobs --region $R --jobs <jobId>          # arrayProperties.statusSummary
#   aws s3 ls {prefix}/outputs/ --recursive | grep _status.json

# 5. Разбор (после серии; minvCpus=0 — простой среды бесплатен, но разобрать её — порядок)
#   aws batch update-job-queue --region $R --job-queue {v['JQ_NAME']} --state DISABLED
#   aws batch delete-job-queue --region $R --job-queue {v['JQ_NAME']}
#   aws batch update-compute-environment --region $R --compute-environment {v['CE_NAME']} --state DISABLED
#   aws batch delete-compute-environment --region $R --compute-environment {v['CE_NAME']}
#   aws batch deregister-job-definition --region $R --job-definition {v['JD_NAME']}:<ревизия>
"""


def render(cfg: dict, region: str, mode: str, rates: dict, git: dict | None = None) -> tuple[dict, dict, list]:
    rows = build_plan(cfg, mode)
    v = variables(cfg, region, mode, rows)
    req = render_templates(v)
    errs = validate_api(req) + validate_limits(req, cfg, rates, region)
    left = [m for name in TEMPLATES for m in PLACEHOLDER_RE.findall(json.dumps(req[name]))]
    if left:
        errs.append(f"неподставленные поля: {sorted(set(left))}")
    ptext = plan_csv(rows)
    job = {"schema": "gtu-batch-job/1", "run_id": cfg["run"]["run_id"], "run_prefix": run_prefix(cfg),
           "plan": "inputs/plan.csv", "mode": mode, "command": cfg["plan"]["command"], "np": int(cfg["job"]["vcpu"]),
           "time_limit_s": int(cfg["job"]["attempt_duration_s"] - cfg["job"]["internal_margin_s"])}
    files = {f"{k}.json": json.dumps(req[k], ensure_ascii=False, indent=2) + "\n" for k in TEMPLATES}
    files["s3/inputs/plan.csv"] = ptext
    files["s3/inputs/job.json"] = json.dumps(job, ensure_ascii=False, indent=2) + "\n"
    files["s3/manifest.json"] = json.dumps(manifest(cfg, region, mode, req, ptext, rows, git or git_info()),
                                           ensure_ascii=False, indent=2) + "\n"
    files["commands.sh"] = commands_sh(cfg, region, v)
    return req, files, errs


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--config", default=str(HERE / "batch.toml"))
    ap.add_argument("--region", default=None, help="по умолчанию run.region; запасной — run.fallback_region")
    ap.add_argument("--mode", choices=["line", "point"], default=None)
    ap.add_argument("--out", default=None)
    a = ap.parse_args(argv)
    cfg = load_config(pathlib.Path(a.config))
    region = a.region or cfg["run"]["region"]
    if region not in (cfg["run"]["region"], cfg["run"]["fallback_region"]):
        print(f"регион {region}: разрешены {cfg['run']['region']} и запасной {cfg['run']['fallback_region']}")
        return 1
    mode = a.mode or cfg["run"]["mode"]
    req, files, errs = render(cfg, region, mode, load_rates())
    if errs:
        print("ОШИБКИ сверки:\n  " + "\n  ".join(errs))
        return 1
    out = pathlib.Path(a.out or HERE / "runs" / cfg["run"]["run_id"] / region)
    for rel, text in files.items():
        p = out / rel
        p.parent.mkdir(parents=True, exist_ok=True)
        p.write_text(text, encoding="utf-8")
    sj = req["submit_array_job"]
    print(f"[render] {region}, режим {mode}: массив {sj['arrayProperties']['size']} заданий, "
          f"тайм-аут {sj['timeout']['attemptDurationSeconds']} с; сверка с моделью API batch (botocore) — без ошибок")
    print(f"[render] файлы: {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
