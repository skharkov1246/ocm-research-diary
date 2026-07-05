#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
routes/halogen_tzvp_aws.py — первый запуск на франкфуртской квоте (eu-central-1,
одобрено 128 vCPU против 16 в us-east-1): рефайн галоген-дескриптора на def2-TZVP
(F/Cl/Br/I, 4 джоба × OMP4 на одном r7i.4xlarge = 16 vCPU). Leak-proof как всегда:
terminate-on-shutdown + watchdog + S3-синк. Результаты halogen_tzvp_<X>.json →
s3://alpha-o-results-097743207937/alpha-o/ffm1/. IAM instance profile глобальный,
S3-бакет глобальный — меняется только регион и сабнет.

Запуск: /opt/qc-venv/bin/python routes/halogen_tzvp_aws.py [--dryrun]
"""
import os
import subprocess
import sys

import boto3
from botocore.exceptions import ClientError

for v in ("AWS_ACCESS_KEY_ID", "AWS_SECRET_ACCESS_KEY", "AWS_SESSION_TOKEN"):
    os.environ.pop(v, None)

REGION = "eu-central-1"                       # Франкфурт: квота 128 vCPU
ACCOUNT = "097743207937"
BUCKET = f"alpha-o-results-{ACCOUNT}"         # бакет в us-east-1, S3 глобален
PREFIX = "alpha-o/ffm1"
PROFILE = "alpha-o-instance-profile"          # IAM глобален
SUBNET = "subnet-057e6f5ca6d3c6d8d"           # default-VPC eu-central-1a
ITYPE = "r7i.4xlarge"                         # 16 vCPU
NAME = "ffm-halogen-tzvp"
REPO = "https://github.com/skharkov1246/ocm-research-diary"
BRANCH = subprocess.run(["git", "rev-parse", "--abbrev-ref", "HEAD"],
                        capture_output=True, text=True).stdout.strip()
MAX_MIN = 480
JOB_MIN = MAX_MIN - 20
XS = "F Cl Br I"

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
sync_up() {{ aws s3 cp routes/ $S3/ --recursive --exclude '*' --include 'halogen_tzvp_*' --include 'halogen_descriptor_tzvp_*' >/dev/null 2>&1 || true; aws s3 cp /root/repo/hx.log $S3/hx.log >/dev/null 2>&1 || true; }}
( while true; do sleep 120; sync_up; done ) &
run_x() {{ Xh=$1; HAL_BASIS=def2-tzvp HAL_SUFFIX=_tzvp OMP_NUM_THREADS=4 timeout 200m python3 -u routes/halogen_descriptor.py run $Xh >> hx.log 2>&1; echo "[done] $Xh" >> hx.log; }}
export -f run_x
echo "{XS}" | tr ' ' '\\n' | xargs -P 4 -I{{}} bash -c 'run_x "$@"' _ {{}}
HAL_SUFFIX=_tzvp python3 -u routes/halogen_descriptor.py merge >> hx.log 2>&1
sync_up
echo "[aws] ALL DONE" >> hx.log
sync_up
poweroff
"""


def main():
    print("caller:", boto3.client("sts", region_name=REGION)
          .get_caller_identity()["Arn"], flush=True)
    print(f"region={REGION} type={ITYPE} branch={BRANCH} "
          f"bucket=s3://{BUCKET}/{PREFIX} X={XS}", flush=True)
    ami = ssm.get_parameter(
        Name="/aws/service/ami-amazon-linux-latest/al2023-ami-kernel-default-x86_64"
    )["Parameter"]["Value"]
    common = dict(
        ImageId=ami, InstanceType=ITYPE, MinCount=1, MaxCount=1,
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
                                     {"Key": "Name", "Value": NAME}]},
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
        {"Name": "tag:Name", "Values": [NAME]},
        {"Name": "instance-state-name", "Values": ["pending", "running"]}])
    if any(i for r in ex["Reservations"] for i in r["Instances"]):
        sys.exit(f"{NAME} instance already running")
    iid = ec2.run_instances(**common)["Instances"][0]["InstanceId"]
    print("launched", iid, flush=True)


if __name__ == "__main__":
    main()
