#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
routes/lao_aws_smoke.py — leak-proof смок-прогон варианта 2 (ЛАО/ПАО):
lao_cr_trimer.py spins→int на r7i.4xlarge в us-west-2 (свежая квота 116 vCPU).
Смок отвечает: сходится ли CASSCF/NEVPT2 на хромациклах и мультиреференсны ли
они (NOON) — есть ли наша ниша в C6/C8-селективности до постройки TS-стадий.
Защита: terminate-on-shutdown + бортовой shutdown + жнец; DONE-маркер отдельно.

Запуск: python3 routes/lao_aws_smoke.py [--dryrun]
"""
import os
import subprocess
import sys
import time

import boto3
from botocore.exceptions import ClientError

REGION = os.environ.get("LAO_REGION", "us-west-2")
ACCOUNT = "097743207937"
BUCKET = f"alpha-o-results-{ACCOUNT}"          # us-east-1, cross-region OK
PREFIX = "alpha-o/lao-smoke"
PROFILE = "alpha-o-instance-profile"
TAG_PROJECT = "ocm-agent"
TAG_NAME = "lao-smoke"
ITYPE = os.environ.get("LAO_ITYPE", "r7i.4xlarge")
REPO = "https://github.com/skharkov1246/ocm-research-diary"
BRANCH = os.environ.get("LAO_BRANCH") or subprocess.run(
    ["git", "rev-parse", "--abbrev-ref", "HEAD"], capture_output=True,
    text=True).stdout.strip()
MAX_MIN = int(os.environ.get("LAO_MAX_MIN", "420"))
STAGES = os.environ.get("LAO_STAGES", "spins int")
# spins не резюмится файлом → его выход чистим; int резюмится ПО-ТЭГОВО
# (тег без nevpt2-блока пересчитывается) → int.json сохраняем (c5 готов)
_OUT = {"spins": "routes/lao_cr_spins.json", "int": "",
        "bhe": "routes/lao_cr_bhe_result.json",
        "ins": "routes/lao_cr_ins_result.json",
        "desc": "routes/lao_cr_desc_results.json"}
RM_FILES = " ".join(f for f in (_OUT[s] for s in STAGES.split()) if f) or "/dev/null"
JOB_MIN = MAX_MIN - 20
STAGE_TO = max(300, JOB_MIN // 2 - 20)
DISK_GB = 60

for var in ("AWS_ACCESS_KEY_ID", "AWS_SECRET_ACCESS_KEY", "AWS_SESSION_TOKEN"):
    os.environ.pop(var, None)

ec2 = boto3.client("ec2", region_name=REGION)
ssm = boto3.client("ssm", region_name=REGION)
s3 = boto3.client("s3", region_name="us-east-1")

USERDATA = f"""#!/bin/bash
shutdown -P +{JOB_MIN}
exec > /tmp/boot.log 2>&1
set -x
export PATH=$PATH:/usr/local/bin:/root/.local/bin
S3=s3://{BUCKET}/{PREFIX}
export AWS_DEFAULT_REGION=us-east-1
for i in 1 2 3; do timeout 10m dnf install -y python3-pip git gcc gcc-c++ && break; sleep 30; done
for i in 1 2 3; do timeout 10m python3 -m pip install -q awscli && break; sleep 30; done
( while true; do aws s3 cp /tmp/boot.log $S3/boot.log >/dev/null 2>&1; sleep 60; done ) &
for i in 1 2 3; do timeout 25m python3 -m pip install -q numpy scipy pyscf geometric && break; sleep 30; done
for i in 1 2 3; do (cd /root && timeout 10m git clone --depth 1 -b {BRANCH} {REPO} repo) && break; sleep 30; done
cd /root/repo || {{ aws s3 cp /tmp/boot.log $S3/setup_failed.log; poweroff; }}
python3 -c "import pyscf" || {{ aws s3 cp /tmp/boot.log $S3/setup_failed.log; poweroff; }}
rm -f {RM_FILES}
for f in lao_cr_bhe.json lao_cr_ins.json; do aws s3 cp $S3/$f routes/$f || true; done
mkdir -p /root/scratch
export PYSCF_TMPDIR=/root/scratch
export OMP_NUM_THREADS=$(nproc)
sync_up() {{ aws s3 cp routes/ $S3/ --recursive --exclude '*' --include 'lao_cr_*.json' >/dev/null 2>&1 || true; aws s3 cp /root/repo/lao.log $S3/lao.log >/dev/null 2>&1 || true; }}
( while true; do sleep 120; sync_up; done ) &
for st in {STAGES}; do
  echo "[aws][lao] stage $st" >> /root/repo/lao.log
  timeout {STAGE_TO}m python3 -u routes/lao_cr_trimer.py $st >> /root/repo/lao.log 2>&1
  sync_up
done
echo "[aws][lao] ALL DONE" >> /root/repo/lao.log
sync_up
echo done > /tmp/DONE && aws s3 cp /tmp/DONE $S3/DONE
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
        {"Name": "tag:Name", "Values": [TAG_NAME]},
        {"Name": "instance-state-name",
         "Values": ["pending", "running", "stopping", "stopped"]}])
    return [i["InstanceId"] for res in r["Reservations"] for i in res["Instances"]]


def main():
    subnet = default_subnet()
    print(f"branch={BRANCH} region={REGION} type={ITYPE} subnet={subnet} "
          f"bucket=s3://{BUCKET}/{PREFIX} reap<={MAX_MIN}min", flush=True)
    try:
        ami = ssm.get_parameter(
            Name="/aws/service/ami-amazon-linux-latest/"
                 "al2023-ami-kernel-default-x86_64")["Parameter"]["Value"]
    except ClientError:
        imgs = ec2.describe_images(
            Owners=["amazon"],
            Filters=[{"Name": "name",
                      "Values": ["al2023-ami-2023*-kernel-*-x86_64"]},
                     {"Name": "state", "Values": ["available"]}])["Images"]
        ami = max(imgs, key=lambda i: i["CreationDate"])["ImageId"]
    common = dict(
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
                      {"Key": "Name", "Value": TAG_NAME}]},
            {"ResourceType": "volume",
             "Tags": [{"Key": "Project", "Value": TAG_PROJECT}]}],
        UserData=USERDATA)
    if "--dryrun" in sys.argv:
        try:
            ec2.run_instances(DryRun=True, **common)
            print("[dryrun] unexpected success", flush=True)
        except ClientError as e:
            print("[dryrun] OK — launch permitted" if "DryRunOperation" in str(e)
                  else f"[dryrun] FAILED: {e}", flush=True)
        return
    existing = my_instances()
    if existing:
        sys.exit(f"lao-smoke instance(s) already exist: {existing}")
    iid = ec2.run_instances(**common)["Instances"][0]["InstanceId"]
    print(f"launched {iid}", flush=True)
    t0 = time.time()
    try:
        while time.time() - t0 < MAX_MIN * 60:
            try:
                s3.head_object(Bucket=BUCKET, Key=f"{PREFIX}/DONE")
                for f in ("lao_cr_bhe_result.json", "lao_cr_ins_result.json",
                          "lao_cr_desc_results.json"):
                    try:
                        s3.download_file(BUCKET, f"{PREFIX}/{f}", f"routes/{f}")
                    except ClientError:
                        pass
                print("[ok] LAO smoke DONE — results pulled", flush=True)
                return
            except ClientError:
                pass
            time.sleep(60)
        print("[timeout] MAX_MIN reached", flush=True)
    finally:
        iids = my_instances()
        if iids:
            print(f"[reap] terminating {iids}", flush=True)
            try:
                ec2.terminate_instances(InstanceIds=iids)
            except ClientError as e:
                print(f"[reap][ERR] {e}", flush=True)


if __name__ == "__main__":
    main()
