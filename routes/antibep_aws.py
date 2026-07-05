#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
routes/antibep_aws.py — запуск анти-BEP дескриптора A (antibep_carbene.py) на одном
EC2 r7i.2xlarge (8 vCPU): панель металл-free триплетных HAT-медиаторов O/NH/O2/CH2,
каждый в своём timeout (зависание одного не блокирует остальных). Leak-proof:
terminate-on-shutdown + бортовой watchdog + S3-синк. Результаты (antibep_<M>.json +
antibep_carbene_results.json) — в s3://alpha-o-results-097743207937/alpha-o/antibep/.

Запуск из песочницы: python3 routes/antibep_aws.py [--dryrun]
"""
import os
import subprocess
import sys

import boto3
from botocore.exceptions import ClientError

for v in ("AWS_ACCESS_KEY_ID", "AWS_SECRET_ACCESS_KEY", "AWS_SESSION_TOKEN"):
    os.environ.pop(v, None)

REGION = "us-east-1"
ACCOUNT = "097743207937"
BUCKET = f"alpha-o-results-{ACCOUNT}"
PREFIX = "alpha-o/antibep"
PROFILE = "alpha-o-instance-profile"
SUBNET = "subnet-0b1a363f27ecbbf12"
REPO = "https://github.com/skharkov1246/ocm-research-diary"
BRANCH = subprocess.run(["git", "rev-parse", "--abbrev-ref", "HEAD"],
                        capture_output=True, text=True).stdout.strip()
MAX_MIN = 360
JOB_MIN = MAX_MIN - 20
MEDS = "O NH O2 CH2"

ec2 = boto3.client("ec2", region_name=REGION)
ssm = boto3.client("ssm", region_name=REGION)

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
sync_up() {{ aws s3 cp routes/ $S3/ --recursive --exclude '*' --include 'antibep_*.json' >/dev/null 2>&1 || true; aws s3 cp /root/repo/ab.log $S3/ab.log >/dev/null 2>&1 || true; }}
( while true; do sleep 120; sync_up; done ) &
# каждый медиатор в своём timeout: зависание одного не блокирует панель
run_med() {{ Mx=$1; OMP_NUM_THREADS=4 timeout 90m python3 -u routes/antibep_carbene.py run $Mx >> ab.log 2>&1; echo "[done] $Mx" >> ab.log; }}
export -f run_med
echo "{MEDS}" | tr ' ' '\\n' | xargs -P 2 -I{{}} bash -c 'run_med "$@"' _ {{}}
python3 -u routes/antibep_carbene.py merge >> ab.log 2>&1
sync_up
echo "[aws] ALL DONE" >> ab.log
sync_up
poweroff
"""


def main():
    print("caller:", boto3.client("sts", region_name=REGION)
          .get_caller_identity()["Arn"], flush=True)
    print(f"branch={BRANCH} bucket=s3://{BUCKET}/{PREFIX} meds={MEDS}", flush=True)
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
                                     {"Key": "Name", "Value": "antibep"}]},
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
        {"Name": "tag:Name", "Values": ["antibep"]},
        {"Name": "instance-state-name", "Values": ["pending", "running"]}])
    if any(i for r in ex["Reservations"] for i in r["Instances"]):
        sys.exit("antibep instance already running")
    iid = ec2.run_instances(**common)["Instances"][0]["InstanceId"]
    print("launched", iid, flush=True)


if __name__ == "__main__":
    main()
