#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
routes/aws_run.py — УНИВЕРСАЛЬНЫЙ leak-proof лаунчер EC2 под любой регион/тип/
полезную нагрузку. Вместо копипасты одноразовых *_aws.py: полезная нагрузка —
bash-сниппет в payloads/<name>.sh (закоммичен в ветку; инстанс исполняет его из
свежего клона с cwd=/root/repo). Контракт payload'а:
  • пишет прогресс в $LOGF (уже определён), результаты кладёт в routes/;
  • о синке не думает — фоновый sync_up каждые 120с зеркалит routes/ по маске
    $SYNC_GLOB и $LOGF в $S3; после payload'а — финальный синк и ALL DONE.
Всегда: shutdown -P дедлайн + terminate-on-shutdown + IMDSv2 + свой tag Name
(дубль-гард). AMI: SSM, при регион-скоупе — describe-images fallback.

Запуск:
  /opt/qc-venv/bin/python routes/aws_run.py --name ffm-chemloop-ext \\
      --payload payloads/chemloop_ext.sh --prefix alpha-o/ffm2 \\
      [--region eu-central-1] [--itype r7i.4xlarge] [--subnet subnet-...]
      [--sync-glob 'chemloop_*'] [--max-min 600] [--dryrun]
"""
import argparse
import os
import subprocess
import sys

import boto3
from botocore.exceptions import ClientError

for v in ("AWS_ACCESS_KEY_ID", "AWS_SECRET_ACCESS_KEY", "AWS_SESSION_TOKEN"):
    os.environ.pop(v, None)

ACCOUNT = "097743207937"
BUCKET = f"alpha-o-results-{ACCOUNT}"
PROFILE = "alpha-o-instance-profile"
REPO = "https://github.com/skharkov1246/ocm-research-diary"
DEFAULT_SUBNET = {                       # default-VPC сабнеты по регионам
    "us-east-1": "subnet-0b1a363f27ecbbf12",
    "eu-central-1": "subnet-057e6f5ca6d3c6d8d",
}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--name", required=True)
    ap.add_argument("--payload", required=True)
    ap.add_argument("--prefix", required=True)
    ap.add_argument("--region", default="eu-central-1")
    ap.add_argument("--itype", default="r7i.4xlarge")
    ap.add_argument("--subnet", default=None)
    ap.add_argument("--sync-glob", default="*")
    ap.add_argument("--max-min", type=int, default=600)
    ap.add_argument("--dryrun", action="store_true")
    a = ap.parse_args()

    branch = subprocess.run(["git", "rev-parse", "--abbrev-ref", "HEAD"],
                            capture_output=True, text=True).stdout.strip()
    payload = open(a.payload).read()
    subnet = a.subnet or DEFAULT_SUBNET[a.region]
    job_min = a.max_min - 20

    userdata = f"""#!/bin/bash
shutdown -P +{job_min}
exec > /tmp/boot.log 2>&1; set -x
export PATH=$PATH:/usr/local/bin:/root/.local/bin
S3=s3://{BUCKET}/{a.prefix}
for i in 1 2 3; do dnf install -y python3-pip git gcc gcc-c++ && break; sleep 20; done
for i in 1 2 3; do python3 -m pip install -q awscli && break; sleep 20; done
( while true; do aws s3 cp /tmp/boot.log $S3/boot.log >/dev/null 2>&1; sleep 60; done ) &
for i in 1 2 3; do python3 -m pip install -q numpy scipy pyscf geometric && break; sleep 20; done
cd /root && git clone --depth 1 -b {branch} {REPO} repo && cd repo || {{ aws s3 cp /tmp/boot.log $S3/setup_failed.log; poweroff; }}
export LOGF=/root/repo/run.log
sync_up() {{ aws s3 cp routes/ $S3/ --recursive --exclude '*' --include '{a.sync_glob}' >/dev/null 2>&1 || true; aws s3 cp $LOGF $S3/run.log >/dev/null 2>&1 || true; }}
export -f sync_up
( while true; do sleep 120; sync_up; done ) &
# ---------------- payload ({a.payload}) ----------------
{payload}
# --------------------------------------------------------
sync_up
echo "[aws] ALL DONE" >> $LOGF
sync_up
poweroff
"""

    ec2 = boto3.client("ec2", region_name=a.region)
    print("caller:", boto3.client("sts", region_name=a.region)
          .get_caller_identity()["Arn"], flush=True)
    print(f"name={a.name} region={a.region} type={a.itype} branch={branch} "
          f"s3={BUCKET}/{a.prefix} glob={a.sync_glob}", flush=True)
    try:
        ami = boto3.client("ssm", region_name=a.region).get_parameter(
            Name="/aws/service/ami-amazon-linux-latest/al2023-ami-kernel-default-x86_64"
        )["Parameter"]["Value"]
    except ClientError:
        imgs = ec2.describe_images(
            Owners=["amazon"],
            Filters=[{"Name": "name", "Values": ["al2023-ami-2023*-kernel-*-x86_64"]},
                     {"Name": "state", "Values": ["available"]}])["Images"]
        ami = max(imgs, key=lambda i: i["CreationDate"])["ImageId"]
    print("ami:", ami, flush=True)
    common = dict(
        ImageId=ami, InstanceType=a.itype, MinCount=1, MaxCount=1,
        IamInstanceProfile={"Name": PROFILE},
        InstanceInitiatedShutdownBehavior="terminate",
        MetadataOptions={"HttpTokens": "required",
                         "HttpPutResponseHopLimit": 1, "HttpEndpoint": "enabled"},
        NetworkInterfaces=[{"DeviceIndex": 0, "SubnetId": subnet,
                            "AssociatePublicIpAddress": True}],
        BlockDeviceMappings=[{"DeviceName": "/dev/xvda",
                              "Ebs": {"VolumeSize": 60, "VolumeType": "gp3",
                                      "DeleteOnTermination": True}}],
        TagSpecifications=[{"ResourceType": "instance",
                            "Tags": [{"Key": "Project", "Value": "ocm-agent"},
                                     {"Key": "Name", "Value": a.name}]},
                           {"ResourceType": "volume",
                            "Tags": [{"Key": "Project", "Value": "ocm-agent"}]}],
        UserData=userdata)
    if a.dryrun:
        try:
            ec2.run_instances(DryRun=True, **common)
        except ClientError as e:
            print("[dryrun] OK" if "DryRunOperation" in str(e) else f"[dryrun] {e}",
                  flush=True)
        return
    ex = ec2.describe_instances(Filters=[
        {"Name": "tag:Name", "Values": [a.name]},
        {"Name": "instance-state-name", "Values": ["pending", "running"]}])
    if any(i for r in ex["Reservations"] for i in r["Instances"]):
        sys.exit(f"{a.name} already running")
    iid = ec2.run_instances(**common)["Instances"][0]["InstanceId"]
    print("launched", iid, flush=True)


if __name__ == "__main__":
    main()
