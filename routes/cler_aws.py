#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
routes/cler_aws.py — leak-proof лаунчер: routes/cler_ladder.py (гейт N1-v2,
ClER-лестница на смешанно-оксидных Ru-Pd-Ir центрах — расчётная поддержка живой
анодной программы Кольской ГМК, ванна №304) на c7i.8xlarge в us-west-2.

Защита от утечки денег: terminate-on-shutdown + бортовой `shutdown -P` ПЕРВОЙ
строкой user-data + жнец в finally; DONE-маркер отдельным ключом; IMDSv2
(HttpTokens=required); диск DeleteOnTermination.

Режимы:
  python3 routes/cler_aws.py            # один инстанс, все 5 центров
  DRYRUN=1 python3 routes/cler_aws.py   # только DryRun-проверка прав, без запуска
  CLER_SPLIT=1 python3 routes/cler_aws.py
      # ДВА инстанса параллельно:
      #   #1 (suffix _a) — моноядерные центры ru, ir, pd  (лёгкие, ~107-141 AO)
      #   #2 (suffix _b) — ru_pd, ru_pd_ir                (тяжёлые, ~166-283 AO)
      # каждый пишет СВОЙ cler_ladder_results{suffix}.json — в S3 не сталкиваются
  CLER_ONLY=b python3 routes/cler_aws.py
      # перезапуск ОДНОЙ ветки; результаты подтягиваются из S3 → резюм с места
      # обрыва (тяжёлой ветке _b штатно нужен второй раунд)

Референсы Cl2/H2O/H2 считает КАЖДАЯ ветка сама (научный скрипт всегда их
добавляет) — ветки полностью независимы.

env-ручки: CLER_ITYPE (default c7i.8xlarge), CLER_REGION (us-west-2),
BRANCH (default claude/lukoil-norilsk-optimization-qbrqpk),
CLER_MAX_MIN (default 1080), CLER_STAGE_TIMEOUT (проброс; по умолчанию 5400 для ветки _a и
10800 для тяжёлой _b/одиночного прогона),
CLER_BASIS (def2-svp), CLER_XC ("pbe,pbe").
"""
import os
import sys
import time

import boto3
from botocore.exceptions import ClientError

REGION = os.environ.get("CLER_REGION", "us-west-2")
ACCOUNT = "097743207937"
BUCKET = f"alpha-o-results-{ACCOUNT}"          # us-east-1, cross-region OK
PREFIX = "alpha-o/cler"
PROFILE = "alpha-o-instance-profile"
TAG_PROJECT = "sandbox-agent"
TAG_NAME = "cler"
ITYPE = os.environ.get("CLER_ITYPE", "c7i.8xlarge")   # 32 vCPU
REPO = "https://github.com/skharkov1246/ocm-research-diary"
BRANCH = os.environ.get("BRANCH", "claude/lukoil-norilsk-optimization-qbrqpk")
MAX_MIN = int(os.environ.get("CLER_MAX_MIN", "1080"))
# Стадия ветки _b (би/триядерные кластеры, 166-283 AO) в 1.5 ч не укладывается,
# поэтому по умолчанию даём ей 3 ч; явный CLER_STAGE_TIMEOUT перебивает обе.
STAGE_TIMEOUT_ENV = os.environ.get("CLER_STAGE_TIMEOUT")
STAGE_TIMEOUT_DEFAULT = {"_a": 5400, "_b": 10800, "": 10800}
CLER_BASIS = os.environ.get("CLER_BASIS", "def2-svp")
CLER_XC = os.environ.get("CLER_XC", "pbe,pbe")
SPLIT = os.environ.get("CLER_SPLIT", "") == "1"
ONLY = os.environ.get("CLER_ONLY", "")
DRYRUN = os.environ.get("DRYRUN", "") == "1" or "--dryrun" in sys.argv
JOB_MIN = MAX_MIN - 20          # бортовой watchdog чуть раньше жнеца
SCRIPT_TO_MIN = JOB_MIN - 40    # timeout самого python-запуска
DISK_GB = 60

# Ветки сплита: _a — три моноядерных центра, _b — би/триядерные (дорогие).
LIGHT_CENTERS = "ru,ir,pd"
HEAVY_CENTERS = "ru_pd,ru_pd_ir"
ALL_CENTERS = "ru,ir,pd,ru_pd,ru_pd_ir"
BRANCH_JOBS = {"a": (LIGHT_CENTERS, "_a"), "b": (HEAVY_CENTERS, "_b")}

# ВАЖНО: env-креды нужны ЛАУНЧЕРУ (песочница вызывает EC2/S3 от scoped-юзера);
# на инстанс они не попадают — user-data кредов не содержит, там instance profile.
ec2 = boto3.client("ec2", region_name=REGION)
ssm = boto3.client("ssm", region_name=REGION)
s3 = boto3.client("s3", region_name="us-east-1")


def stage_timeout_for(suffix):
    return int(STAGE_TIMEOUT_ENV) if STAGE_TIMEOUT_ENV \
        else STAGE_TIMEOUT_DEFAULT.get(suffix, 5400)


def userdata(centers_csv, suffix):
    """user-data одного инстанса: watchdog ПЕРВЫМ, клон ветки, pip-зависимости,
    резюм results-JSON из S3, запуск научного скрипта, выгрузка, poweroff."""
    rj = f"cler_ladder_results{suffix}.json"
    log = f"cler{suffix}.log"
    stage_to = stage_timeout_for(suffix)
    return f"""#!/bin/bash
shutdown -P +{JOB_MIN}
exec > /tmp/boot.log 2>&1
set -x
export PATH=$PATH:/usr/local/bin:/root/.local/bin
S3=s3://{BUCKET}/{PREFIX}
export AWS_DEFAULT_REGION=us-east-1
for i in 1 2 3; do timeout 10m dnf install -y python3-pip git gcc gcc-c++ && break; sleep 30; done
for i in 1 2 3; do timeout 10m python3 -m pip install -q awscli && break; sleep 30; done
( while true; do aws s3 cp /tmp/boot.log $S3/boot{suffix}.log >/dev/null 2>&1; sleep 60; done ) &
for i in 1 2 3; do timeout 25m python3 -m pip install -q numpy scipy pyscf pyberny && break; sleep 30; done
for i in 1 2 3; do (cd /root && timeout 10m git clone --depth 1 -b {BRANCH} {REPO} repo) && break; sleep 30; done
cd /root/repo || {{ aws s3 cp /tmp/boot.log $S3/setup_failed{suffix}.log; poweroff; }}
python3 -c "import pyscf" || {{ aws s3 cp /tmp/boot.log $S3/setup_failed{suffix}.log; poweroff; }}
aws s3 cp $S3/{rj} routes/{rj} || true
mkdir -p /root/scratch
export PYSCF_TMPDIR=/root/scratch
export OMP_NUM_THREADS=$(nproc)
export CLER_CENTERS={centers_csv}
export CLER_JSON_SUFFIX={suffix}
export CLER_STAGE_TIMEOUT={stage_to}
export CLER_BASIS={CLER_BASIS}
export CLER_XC={CLER_XC}
sync_up() {{ aws s3 cp routes/{rj} $S3/{rj} >/dev/null 2>&1 || true; aws s3 cp /root/repo/{log} $S3/{log} >/dev/null 2>&1 || true; }}
( while true; do sleep 120; sync_up; done ) &
echo "[aws][cler{suffix}] start centers={centers_csv} basis={CLER_BASIS} xc={CLER_XC} stage_to={stage_to}s" >> /root/repo/{log}
timeout {SCRIPT_TO_MIN}m python3 -u routes/cler_ladder.py >> /root/repo/{log} 2>&1
echo "[aws][cler{suffix}] ALL DONE" >> /root/repo/{log}
sync_up
echo done > /tmp/DONE && aws s3 cp /tmp/DONE $S3/DONE{suffix}
poweroff
"""


def default_subnet():
    vpcs = ec2.describe_vpcs(Filters=[{"Name": "is-default",
                                       "Values": ["true"]}])["Vpcs"]
    if not vpcs:
        sys.exit(f"no default VPC in {REGION}")
    subs = ec2.describe_subnets(Filters=[
        {"Name": "vpc-id", "Values": [vpcs[0]["VpcId"]]}])["Subnets"]
    subs = [s for s in subs if s.get("MapPublicIpOnLaunch")] or subs
    return subs[0]["SubnetId"]


def my_instances():
    r = ec2.describe_instances(Filters=[
        {"Name": "tag:Project", "Values": [TAG_PROJECT]},
        {"Name": "tag:Name", "Values": [TAG_NAME + "*"]},
        {"Name": "instance-state-name",
         "Values": ["pending", "running", "stopping", "stopped"]}])
    return [i["InstanceId"] for res in r["Reservations"] for i in res["Instances"]]


def latest_ami():
    try:
        return ssm.get_parameter(
            Name="/aws/service/ami-amazon-linux-latest/"
                 "al2023-ami-kernel-default-x86_64")["Parameter"]["Value"]
    except ClientError:
        imgs = ec2.describe_images(
            Owners=["amazon"],
            Filters=[{"Name": "name",
                      "Values": ["al2023-ami-2023*-kernel-*-x86_64"]},
                     {"Name": "state", "Values": ["available"]}])["Images"]
        return max(imgs, key=lambda i: i["CreationDate"])["ImageId"]


def launch_spec(ami, subnet, name, centers_csv, suffix):
    return dict(
        ImageId=ami, InstanceType=ITYPE, MinCount=1, MaxCount=1,
        IamInstanceProfile={"Name": PROFILE},
        InstanceInitiatedShutdownBehavior="terminate",
        MetadataOptions={"HttpTokens": "required",
                         "HttpPutResponseHopLimit": 1, "HttpEndpoint": "enabled"},
        NetworkInterfaces=[{"DeviceIndex": 0, "SubnetId": subnet,
                            "AssociatePublicIpAddress": True}],
        BlockDeviceMappings=[{"DeviceName": "/dev/xvda",
                              "Ebs": {"VolumeSize": DISK_GB, "VolumeType": "gp3",
                                      "DeleteOnTermination": True}}],
        TagSpecifications=[
            {"ResourceType": "instance",
             "Tags": [{"Key": "Project", "Value": TAG_PROJECT},
                      {"Key": "Name", "Value": name}]},
            {"ResourceType": "volume",
             "Tags": [{"Key": "Project", "Value": TAG_PROJECT}]}],
        UserData=userdata(centers_csv, suffix))


def main():
    if ONLY in BRANCH_JOBS:                 # перезапуск одной ветки (resume по S3)
        centers_csv, suffix = BRANCH_JOBS[ONLY]
        jobs = [(f"{TAG_NAME}-{ONLY}", centers_csv, suffix)]
    elif ONLY:
        sys.exit(f"CLER_ONLY={ONLY!r}: допустимо только 'a' или 'b'")
    elif SPLIT:
        jobs = [(f"{TAG_NAME}-{b}", *BRANCH_JOBS[b]) for b in ("a", "b")]
    else:
        jobs = [(TAG_NAME, ALL_CENTERS, "")]
    subnet = default_subnet()
    ami = latest_ami()
    print(f"branch={BRANCH} region={REGION} type={ITYPE} subnet={subnet} "
          f"bucket=s3://{BUCKET}/{PREFIX} split={SPLIT} only={ONLY or '-'} "
          f"reap<={MAX_MIN}min", flush=True)

    if DRYRUN:
        for name, centers_csv, suffix in jobs:
            try:
                ec2.run_instances(DryRun=True,
                                  **launch_spec(ami, subnet, name,
                                                centers_csv, suffix))
                print(f"[dryrun] {name}: unexpected success", flush=True)
            except ClientError as e:
                print(f"[dryrun] {name} centers={centers_csv} "
                      f"stage_timeout={stage_timeout_for(suffix)}s: "
                      + ("OK — launch permitted" if "DryRunOperation" in str(e)
                         else f"FAILED: {e}"), flush=True)
        return

    existing = my_instances()
    if existing:
        sys.exit(f"cler instance(s) already exist: {existing}")

    iids = []
    for name, centers_csv, suffix in jobs:
        iid = ec2.run_instances(**launch_spec(ami, subnet, name, centers_csv,
                                              suffix))["Instances"][0]["InstanceId"]
        iids.append(iid)
        print(f"launched {name} (centers={centers_csv}, "
              f"stage_timeout={stage_timeout_for(suffix)}s): {iid}", flush=True)

    done_keys = {f"{PREFIX}/DONE{suffix}": f"cler_ladder_results{suffix}.json"
                 for _n, _c, suffix in jobs}
    pending = set(done_keys)
    t0 = time.time()
    try:
        while time.time() - t0 < MAX_MIN * 60 and pending:
            for dk in sorted(pending):
                try:
                    s3.head_object(Bucket=BUCKET, Key=dk)
                except ClientError:
                    continue
                rj = done_keys[dk]
                try:
                    s3.download_file(BUCKET, f"{PREFIX}/{rj}",
                                     os.path.join("routes", rj))
                    print(f"[ok] {dk} — pulled routes/{rj}", flush=True)
                except ClientError as e:
                    print(f"[ok] {dk}, but pull {rj} failed: {e}", flush=True)
                pending.discard(dk)
            if pending:
                time.sleep(60)
        if pending:
            print(f"[timeout] MAX_MIN reached, pending: {sorted(pending)}",
                  flush=True)
        else:
            print("[ok] all instances DONE — results pulled", flush=True)
    finally:
        leftovers = my_instances()
        if leftovers:
            print(f"[reap] terminating {leftovers}", flush=True)
            try:
                ec2.terminate_instances(InstanceIds=leftovers)
            except ClientError as e:
                print(f"[reap][ERR] {e}", flush=True)


if __name__ == "__main__":
    main()
