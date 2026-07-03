#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
calc/aws_launch.py — boto3-оркестратор AWS-прогонов дневника (EC2 + S3).

Паттерн тот же, что в calc/aws_run.sh (dead-man switch), но запускается ИЗ СЕССИИ:
инстанс поднимается с InstanceInitiatedShutdownBehavior=terminate, в user-data
ставит стек (pyscf), клонирует публичный репозиторий, запускает job-команду,
заливает результаты в S3 и гасит себя. Плановый shutdown в начале user-data —
жёсткий таймаут: инстанс не переживёт зависший расчёт.

Команды:
  python3 calc/aws_launch.py infra                # bucket + IAM-роль (однократно)
  python3 calc/aws_launch.py launch <job>         # запустить джоб из JOBS
  python3 calc/aws_launch.py status               # инстансы + содержимое S3
  python3 calc/aws_launch.py fetch <job>          # забрать результаты в рабочее дерево
  python3 calc/aws_launch.py kill <job|all>       # аварийно погасить

Джобы (JOBS ниже): что считать, чем, на каком типе инстанса, какие файлы забрать.
Честность: сюда попадают только скрипты, уже лежащие в репо, — числа приезжают
result-JSON'ами, руками ничего не вписывается.
"""
import json, os, sys, time

import boto3
from botocore.exceptions import ClientError

REGION = os.environ.get("AWS_DEFAULT_REGION", "us-east-1")
REPO = "https://github.com/skharkov1246/ocm-research-diary.git"
# Существующая инфраструктура аккаунта (создана владельцем): у сессионного
# пользователя alpha-o-agent нет IAM/CreateBucket, но есть ec2:RunInstances,
# iam:PassRole на alpha-o-instance-profile и s3:GetObject на префикс alpha-o/.
# Инстансы пишут результаты через профиль; сессия читает из префикса.
PROFILE = "alpha-o-instance-profile"
BUCKET = "alpha-o-results-097743207937"
PREFIX = "alpha-o"
FALLBACK_AMI = "ami-051bfa33df3949860"   # AMI работающего alpha-o инстанса
TAG = {"Key": "project", "Value": "ocm-research-diary"}

JOBS = {
    # Track #17 (CH4->metanol): реальный IZA-CHA кластер, CASSCF(+NEVPT2) spin ordering
    "cha-sp": dict(
        cmd="CHA_NEVPT2=1 python3 -u calc/fe_zeolite_cha_sp.py",
        instance="r7i.2xlarge", hours=10,
        results=["calc/fe_zeolite_cha_results.json"]),
    # Track #13 (CO2->C2+): полное AVAS-пространство CAS(36e,22o) на эндпойнтах Cu2
    "co2-fullcas": dict(
        cmd="python3 -u routes/co2_to_fuels_fullcas.py",
        instance="r7i.4xlarge", hours=20,
        results=["routes/co2_to_fuels_fullcas_results.json"]),
}


def bucket_name():
    return BUCKET


def ensure_infra():
    """Инфраструктура уже существует (bucket + instance profile). Проверяем только
    читаемость префикса результатов из сессии (NoSuchKey = авторизация есть)."""
    s3 = boto3.client("s3", region_name=REGION)
    try:
        s3.get_object(Bucket=BUCKET, Key=f"{PREFIX}/__probe__")
    except ClientError as e:
        code = e.response["Error"]["Code"]
        if code != "NoSuchKey":
            raise RuntimeError(f"prefix {PREFIX}/ not readable: {code}")
    print(f"infra ok: s3://{BUCKET}/{PREFIX}/ readable, profile {PROFILE}")
    return True


def _ami(ec2):
    try:
        ssm = boto3.client("ssm", region_name=REGION)
        return ssm.get_parameter(
            Name="/aws/service/ami-amazon-linux-latest/al2023-ami-kernel-default-x86_64"
        )["Parameter"]["Value"]
    except Exception:
        return FALLBACK_AMI


def _user_data(job, spec, embed_creds):
    b = bucket_name()
    uploads = "\n".join(
        f"aws s3 cp {f} s3://{b}/{PREFIX}/{job}/$(basename {f}) || true" for f in spec["results"])
    minutes = int(spec["hours"] * 60)
    cred_lines = ""
    if embed_creds:
        c = boto3.Session().get_credentials().get_frozen_credentials()
        cred_lines = (f"export AWS_ACCESS_KEY_ID={c.access_key}\n"
                      f"export AWS_SECRET_ACCESS_KEY={c.secret_key}\n"
                      f"export AWS_DEFAULT_REGION={REGION}\n")
    return f"""#!/bin/bash
set -ux
shutdown -h +{minutes}
{cred_lines}dnf install -y -q git python3-pip || yum install -y -q git python3-pip
python3 -m pip install -q numpy scipy h5py pyscf
cd /root && git clone --depth 1 {REPO} job && cd job
( {spec['cmd']} ) > run.log 2>&1
echo $? > exit_code
aws s3 cp run.log s3://{b}/{PREFIX}/{job}/run.log || true
aws s3 cp exit_code s3://{b}/{PREFIX}/{job}/exit_code || true
{uploads}
shutdown -h now
"""


def launch(job, has_iam):
    spec = JOBS[job]
    ec2 = boto3.client("ec2", region_name=REGION)
    kw = dict(IamInstanceProfile={"Name": PROFILE}) if has_iam else {}
    r = ec2.run_instances(
        ImageId=_ami(ec2), InstanceType=spec["instance"],
        MinCount=1, MaxCount=1,
        InstanceInitiatedShutdownBehavior="terminate",
        UserData=_user_data(job, spec, embed_creds=not has_iam),
        **kw,
        BlockDeviceMappings=[{"DeviceName": "/dev/xvda",
                              "Ebs": {"VolumeSize": 60, "VolumeType": "gp3",
                                      "DeleteOnTermination": True}}],
        TagSpecifications=[{"ResourceType": "instance",
                            "Tags": [TAG, {"Key": "job", "Value": job},
                                     {"Key": "Name", "Value": f"ocm-{job}"}]}])
    iid = r["Instances"][0]["InstanceId"]
    print(f"launched {job}: {iid} ({spec['instance']}, dead-man {spec['hours']}h)")
    return iid


def status():
    ec2 = boto3.client("ec2", region_name=REGION)
    r = ec2.describe_instances(Filters=[
        {"Name": "tag:project", "Values": [TAG["Value"]]},
        {"Name": "instance-state-name",
         "Values": ["pending", "running", "shutting-down", "stopping"]}])
    for res in r["Reservations"]:
        for i in res["Instances"]:
            tags = {t["Key"]: t["Value"] for t in i.get("Tags", [])}
            print(f"  {i['InstanceId']}  {i['State']['Name']:13s}  "
                  f"{i['InstanceType']:12s}  job={tags.get('job','?')}  "
                  f"up since {i.get('LaunchTime')}")
    s3 = boto3.client("s3", region_name=REGION)
    for job, spec in JOBS.items():
        for f in ["exit_code", "run.log"] + [os.path.basename(x) for x in spec["results"]]:
            key = f"{PREFIX}/{job}/{f}"
            try:
                h = s3.head_object(Bucket=BUCKET, Key=key)
                print(f"  s3://{BUCKET}/{key}  {h['ContentLength']}B  {h['LastModified']}")
            except ClientError:
                pass


def fetch(job):
    spec = JOBS[job]
    s3 = boto3.client("s3", region_name=REGION)
    b = bucket_name()
    got = []
    for f in spec["results"] + ["run.log", "exit_code"]:
        key = f"{PREFIX}/{job}/{os.path.basename(f)}"
        dst = f if "/" in f else f"/tmp/{job}_{f}"
        try:
            s3.download_file(b, key, dst)
            got.append(dst)
        except ClientError:
            print(f"  (no s3://{b}/{key} yet)")
    print("fetched:", got)


def kill(which):
    ec2 = boto3.client("ec2", region_name=REGION)
    flt = [{"Name": "tag:project", "Values": [TAG["Value"]]},
           {"Name": "instance-state-name", "Values": ["pending", "running"]}]
    if which != "all":
        flt.append({"Name": "tag:job", "Values": [which]})
    ids = [i["InstanceId"] for r in ec2.describe_instances(Filters=flt)["Reservations"]
           for i in r["Instances"]]
    if ids:
        ec2.terminate_instances(InstanceIds=ids)
    print("terminated:", ids or "nothing running")


if __name__ == "__main__":
    cmd = sys.argv[1] if len(sys.argv) > 1 else "status"
    if cmd == "infra":
        ensure_infra()
    elif cmd == "launch":
        launch(sys.argv[2], has_iam=ensure_infra())
    elif cmd == "fetch":
        fetch(sys.argv[2])
    elif cmd == "kill":
        kill(sys.argv[2] if len(sys.argv) > 2 else "all")
    else:
        status()
