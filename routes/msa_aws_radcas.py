#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
routes/msa_aws_radcas.py — leak-proof EC2-прогон MSA radical-CAS фикса
(копилка: NEVPT2-вердикт цепи CH4+SO3 загрязнён кривым CAS радикала CH3SO3• —
барьеры ушли в минус). Гоняем routes/msa_chain.py add→hat→merge с MSA_RADCAS=1
(симметричный CAS по всем трём O). Газофазные молекулы — маленький инстанс
us-east-1 (eu-central занят Этапами 18/19).

Защита от утечек — как в ocm_aws_stage17.py: terminate-on-shutdown, бортовой
`shutdown -P`, внешний жнец по дедлайну. Урок Этапа 18: стейл-JSONы из клона
УДАЛЯЮТСЯ на боте; готовность сигналит отдельный маркер DONE, не файл результата.

Запуск: AWS_DEFAULT_REGION=us-east-1 python3 routes/msa_aws_radcas.py [--dryrun]
"""
import os
import subprocess
import sys
import time

import boto3
from botocore.exceptions import ClientError

REGION = os.environ.get("AWS_DEFAULT_REGION", "us-east-1")
ACCOUNT = "097743207937"
BUCKET = f"alpha-o-results-{ACCOUNT}"
PREFIX = os.environ.get("MSA_PREFIX", "alpha-o/msa-radcas")
PROFILE = "alpha-o-instance-profile"
TAG_PROJECT = "ocm-agent"
TAG_NAME = "msa-radcas"
ITYPE = os.environ.get("MSA_ITYPE", "r7i.4xlarge")      # 16 vCPU — газофаза
SUBNET = {"eu-central-1": "subnet-057e6f5ca6d3c6d8d",
          "us-east-1": "subnet-0b1a363f27ecbbf12"}[REGION]
REPO = "https://github.com/skharkov1246/ocm-research-diary"
BRANCH = os.environ.get("MSA_BRANCH") or subprocess.run(
    ["git", "rev-parse", "--abbrev-ref", "HEAD"], capture_output=True,
    text=True).stdout.strip()
MAX_MIN = int(os.environ.get("MSA_MAX_MIN", "600"))
JOB_MIN = MAX_MIN - 20
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
rm -f routes/msa_add.json routes/msa_hat.json routes/msa_results.json
mkdir -p /root/scratch
export PYSCF_TMPDIR=/root/scratch
export MSA_RADCAS=1
export OMP_NUM_THREADS=$(nproc)
sync_up() {{ aws s3 cp routes/ $S3/ --recursive --exclude '*' --include 'msa_*.json' >/dev/null 2>&1 || true; aws s3 cp /root/repo/msa.log $S3/msa.log >/dev/null 2>&1 || true; }}
( while true; do sleep 120; sync_up; done ) &
for st in add hat merge; do
  echo "[aws][msa] stage $st" >> /root/repo/msa.log
  timeout 480m python3 -u routes/msa_chain.py $st >> /root/repo/msa.log 2>&1
  sync_up
done
echo "[aws][msa] ALL DONE" >> /root/repo/msa.log
sync_up
echo done > /tmp/DONE && aws s3 cp /tmp/DONE $S3/DONE
poweroff
"""


def my_instances():
    r = ec2.describe_instances(Filters=[
        {"Name": "tag:Project", "Values": [TAG_PROJECT]},
        {"Name": "tag:Name", "Values": [TAG_NAME]},
        {"Name": "instance-state-name",
         "Values": ["pending", "running", "stopping", "stopped"]}])
    return [i["InstanceId"] for res in r["Reservations"] for i in res["Instances"]]


def main():
    print(f"branch={BRANCH} region={REGION} type={ITYPE} "
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
        NetworkInterfaces=[{"DeviceIndex": 0, "SubnetId": SUBNET,
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
        sys.exit(f"msa-radcas instance(s) already exist: {existing}")
    iid = ec2.run_instances(**common)["Instances"][0]["InstanceId"]
    print(f"launched {iid}", flush=True)
    t0 = time.time()
    try:
        while time.time() - t0 < MAX_MIN * 60:
            try:
                s3.head_object(Bucket=BUCKET, Key=f"{PREFIX}/DONE")
                for f in ("msa_add.json", "msa_hat.json", "msa_results.json"):
                    try:
                        s3.download_file(BUCKET, f"{PREFIX}/{f}", f"routes/{f}")
                    except ClientError:
                        pass
                print("[ok] MSA radcas DONE — results pulled", flush=True)
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
