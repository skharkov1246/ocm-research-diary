#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
routes/chemloop_aws.py — запуск вакансионного волкана (chemloop_vacancy.py) на
одном EC2 r7i.2xlarge (8 vCPU): панель из 9 редокс-металлов, dft+corr каждый,
параллельно на ядрах инстанса. Leak-proof: terminate-on-shutdown + бортовой
watchdog. Результаты (chemloop_<M>.json + chemloop_vacancy_results.json) —
в s3://alpha-o-results-097743207937/alpha-o/chemloop/. Прогресс тоже в S3.

Полный EC2-доступ у alpha-o-agent → можем сами терминировать/чистить.
Запуск из песочницы: python3 routes/chemloop_aws.py [--dryrun]
"""
import os
import subprocess
import sys
import time

import boto3
from botocore.exceptions import ClientError

for v in ("AWS_ACCESS_KEY_ID", "AWS_SECRET_ACCESS_KEY", "AWS_SESSION_TOKEN"):
    os.environ.pop(v, None)

REGION = "us-east-1"
ACCOUNT = "097743207937"
BUCKET = f"alpha-o-results-{ACCOUNT}"
PREFIX = "alpha-o/chemloop"
PROFILE = "alpha-o-instance-profile"
SUBNET = "subnet-0b1a363f27ecbbf12"
REPO = "https://github.com/skharkov1246/ocm-research-diary"
BRANCH = subprocess.run(["git", "rev-parse", "--abbrev-ref", "HEAD"],
                        capture_output=True, text=True).stdout.strip()
MAX_MIN = 600
JOB_MIN = MAX_MIN - 20
METALS = "Ti V Cr Mn Fe Co Ni Cu Ce"

ec2 = boto3.client("ec2", region_name=REGION)
ssm = boto3.client("ssm", region_name=REGION)
s3 = boto3.client("s3", region_name=REGION)

USERDATA = f"""#!/bin/bash
shutdown -P +{JOB_MIN}
exec > /tmp/boot.log 2>&1; set -x
export PATH=$PATH:/usr/local/bin:/root/.local/bin
S3=s3://{BUCKET}/{PREFIX}
for i in 1 2 3; do dnf install -y python3-pip git gcc gcc-c++ && break; sleep 20; done
for i in 1 2 3; do python3 -m pip install -q awscli && break; sleep 20; done
( while true; do aws s3 cp /tmp/boot.log $S3/boot.log >/dev/null 2>&1; sleep 60; done ) &
for i in 1 2 3; do python3 -m pip install -q numpy scipy pyscf geometric && break; sleep 20; done
cd /root && git clone --depth 1 -b {BRANCH} {REPO} repo && cd repo || {{ aws s3 cp /tmp/boot.log $S3/setup_failed.log; poweroff; }}
sync_up() {{ aws s3 cp routes/ $S3/ --recursive --exclude '*' --include 'chemloop_*' >/dev/null 2>&1 || true; aws s3 cp /root/repo/cl.log $S3/cl.log >/dev/null 2>&1 || true; }}
( while true; do sleep 120; sync_up; done ) &
# панель: 2 металла параллельно, 4 потока каждый (8 vCPU)
run_metal() {{ M=$1; OMP_NUM_THREADS=4 timeout 120m python3 -u routes/chemloop_vacancy.py dft $M >> cl.log 2>&1; OMP_NUM_THREADS=4 timeout 120m python3 -u routes/chemloop_vacancy.py corr $M >> cl.log 2>&1; echo "[done] $M" >> cl.log; }}
export -f run_metal
echo "{METALS}" | tr ' ' '\\n' | xargs -P 2 -I{{}} bash -c 'run_metal "$@"' _ {{}}
python3 -u routes/chemloop_vacancy.py merge >> cl.log 2>&1
sync_up
echo "[aws] ALL DONE" >> cl.log
sync_up
poweroff
"""


def main():
    print("caller:", boto3.client("sts", region_name=REGION)
          .get_caller_identity()["Arn"], flush=True)
    print(f"branch={BRANCH} bucket=s3://{BUCKET}/{PREFIX} metals={METALS}", flush=True)
    ami = ssm.get_parameter(
        Name="/aws/service/ami-amazon-linux-latest/al2023-ami-kernel-default-x86_64"
    )["Parameter"]["Value"]
    common = dict(
        ImageId=ami, InstanceType="r7i.2xlarge", MinCount=1, MaxCount=1,
        IamInstanceProfile={"Name": PROFILE},
        InstanceInitiatedShutdownBehavior="terminate",
        MetadataOptions={"HttpTokens": "required",
                         "HttpPutResponseHopLimit": 1, "HttpEndpoint": "enabled"},
        NetworkInterfaces=[{"DeviceIndex": 0, "SubnetId": SUBNET,
                            "AssociatePublicIpAddress": True}],
        BlockDeviceMappings=[{"DeviceName": "/dev/xvda",
                              "Ebs": {"VolumeSize": 60, "VolumeType": "gp3",
                                      "DeleteOnTermination": True}}],
        TagSpecifications=[{"ResourceType": "instance",
                            "Tags": [{"Key": "Project", "Value": "ocm-agent"},
                                     {"Key": "Name", "Value": "chemloop"}]},
                           {"ResourceType": "volume",
                            "Tags": [{"Key": "Project", "Value": "ocm-agent"}]}],
        UserData=USERDATA)
    if "--dryrun" in sys.argv:
        try:
            ec2.run_instances(DryRun=True, **common)
        except ClientError as e:
            print("[dryrun] OK" if "DryRunOperation" in str(e) else f"[dryrun] {e}",
                  flush=True)
        return
    ex = ec2.describe_instances(Filters=[
        {"Name": "tag:Name", "Values": ["chemloop"]},
        {"Name": "instance-state-name", "Values": ["pending", "running"]}])
    if any(i for r in ex["Reservations"] for i in r["Instances"]):
        sys.exit("chemloop instance already running")
    iid = ec2.run_instances(**common)["Instances"][0]["InstanceId"]
    print("launched", iid, flush=True)


if __name__ == "__main__":
    main()
