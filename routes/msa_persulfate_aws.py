#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
routes/msa_persulfate_aws.py — leak-proof лаунчер: routes/msa_persulfate_echem.py
(гейт T2 — потенциал генерации SO4^•- из персульфата, трек E) на c7i.4xlarge
в us-west-2. Системы мелкие (<=10 атомов, def2-SVP) — 16 vCPU достаточно.
Защита по образцу oer_sulfate_aws.py: terminate-on-shutdown + бортовой
shutdown ПЕРВОЙ строкой + жнец в finally; DONE-маркер отдельно; IMDSv2
(HttpTokens=required). Один инстанс.

Режимы:
  python3 routes/msa_persulfate_aws.py            # запуск
  DRYRUN=1 python3 routes/msa_persulfate_aws.py   # только DryRun-проверка прав

env-ручки: MSAE_ITYPE (default c7i.4xlarge), BRANCH
(default claude/lukoil-norilsk-optimization-qbrqpk), MSAE_MAX_MIN (default 360),
MSAE_STAGE_TIMEOUT (проброс в научный скрипт, default 3600),
MSAE_SPECIES (проброс: подмножество частиц, default все), MSAE_BASIS.
"""
import os
import sys
import time

import boto3
from botocore.exceptions import ClientError

REGION = os.environ.get("MSAE_REGION", "us-west-2")
ACCOUNT = "097743207937"
BUCKET = f"alpha-o-results-{ACCOUNT}"          # us-east-1, cross-region OK
PREFIX = "alpha-o/msa-echem"
PROFILE = "alpha-o-instance-profile"
TAG_PROJECT = "sandbox-agent"
TAG_NAME = "msa-echem"
ITYPE = os.environ.get("MSAE_ITYPE", "c7i.4xlarge")   # 16 vCPU — системы мелкие
REPO = "https://github.com/skharkov1246/ocm-research-diary"
BRANCH = os.environ.get("BRANCH", "claude/lukoil-norilsk-optimization-qbrqpk")
MAX_MIN = int(os.environ.get("MSAE_MAX_MIN", "360"))
STAGE_TIMEOUT = int(os.environ.get("MSAE_STAGE_TIMEOUT", "3600"))
MSAE_SPECIES = os.environ.get("MSAE_SPECIES", "")
MSAE_BASIS = os.environ.get("MSAE_BASIS", "def2-svp")
DRYRUN = os.environ.get("DRYRUN", "") == "1" or "--dryrun" in sys.argv
JOB_MIN = MAX_MIN - 20          # бортовой watchdog чуть раньше жнеца
SCRIPT_TO_MIN = JOB_MIN - 40    # timeout самого python-запуска
DISK_GB = 60

RESULTS_JSON = "msa_persulfate_echem_results.json"
LOG = "msa-echem.log"

# ВАЖНО: env-креды нужны ЛАУНЧЕРУ (песочница вызывает EC2/S3 от scoped-юзера);
# на инстанс они не попадают — user-data кредов не содержит, там instance profile.
ec2 = boto3.client("ec2", region_name=REGION)
ssm = boto3.client("ssm", region_name=REGION)
s3 = boto3.client("s3", region_name="us-east-1")


def userdata():
    """user-data: watchdog ПЕРВЫМ, клон ветки, pip-зависимости, резюм
    results-JSON из S3, запуск научного скрипта, выгрузка, poweroff."""
    return f"""#!/bin/bash
shutdown -P +{JOB_MIN}
exec > /tmp/boot.log 2>&1
set -x
export PATH=$PATH:/usr/local/bin:/root/.local/bin
S3=s3://{BUCKET}/{PREFIX}
export AWS_DEFAULT_REGION=us-east-1
for i in 1 2 3; do timeout 10m dnf install -y python3-pip git gcc gcc-c++ && break; sleep 30; done
for i in 1 2 3; do timeout 10m python3 -m pip install -q awscli && break; sleep 30; done
( while true; do aws s3 cp /tmp/boot.log $S3/boot.log >/dev/null 2>&1; sleep 60; done ) &
for i in 1 2 3; do timeout 25m python3 -m pip install -q numpy scipy pyscf pyberny && break; sleep 30; done
for i in 1 2 3; do (cd /root && timeout 10m git clone --depth 1 -b {BRANCH} {REPO} repo) && break; sleep 30; done
cd /root/repo || {{ aws s3 cp /tmp/boot.log $S3/setup_failed.log; poweroff; }}
python3 -c "import pyscf" || {{ aws s3 cp /tmp/boot.log $S3/setup_failed.log; poweroff; }}
aws s3 cp $S3/{RESULTS_JSON} routes/{RESULTS_JSON} || true
mkdir -p /root/scratch
export PYSCF_TMPDIR=/root/scratch
export OMP_NUM_THREADS=$(nproc)
export MSAE_STAGE_TIMEOUT={STAGE_TIMEOUT}
export MSAE_SPECIES={MSAE_SPECIES}
export MSAE_BASIS={MSAE_BASIS}
sync_up() {{ aws s3 cp routes/{RESULTS_JSON} $S3/{RESULTS_JSON} >/dev/null 2>&1 || true; aws s3 cp /root/repo/{LOG} $S3/{LOG} >/dev/null 2>&1 || true; }}
( while true; do sleep 120; sync_up; done ) &
echo "[aws][msa-echem] start species={MSAE_SPECIES or 'all'}" >> /root/repo/{LOG}
timeout {SCRIPT_TO_MIN}m python3 -u routes/msa_persulfate_echem.py run >> /root/repo/{LOG} 2>&1
echo "[aws][msa-echem] ALL DONE" >> /root/repo/{LOG}
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


def launch_spec(ami, subnet):
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
                      {"Key": "Name", "Value": TAG_NAME}]},
            {"ResourceType": "volume",
             "Tags": [{"Key": "Project", "Value": TAG_PROJECT}]}],
        UserData=userdata())


def main():
    subnet = default_subnet()
    ami = latest_ami()
    print(f"branch={BRANCH} region={REGION} type={ITYPE} subnet={subnet} "
          f"bucket=s3://{BUCKET}/{PREFIX} reap<={MAX_MIN}min "
          f"species={MSAE_SPECIES or 'all'}", flush=True)

    if DRYRUN:
        try:
            ec2.run_instances(DryRun=True, **launch_spec(ami, subnet))
            print(f"[dryrun] {TAG_NAME}: unexpected success", flush=True)
        except ClientError as e:
            print(f"[dryrun] {TAG_NAME}: "
                  + ("OK — launch permitted" if "DryRunOperation" in str(e)
                     else f"FAILED: {e}"), flush=True)
        return

    existing = my_instances()
    if existing:
        sys.exit(f"msa-echem instance(s) already exist: {existing}")

    iid = ec2.run_instances(**launch_spec(ami, subnet))["Instances"][0]["InstanceId"]
    print(f"launched {TAG_NAME}: {iid}", flush=True)

    done_key = f"{PREFIX}/DONE"
    t0 = time.time()
    done = False
    try:
        while time.time() - t0 < MAX_MIN * 60 and not done:
            try:
                s3.head_object(Bucket=BUCKET, Key=done_key)
                done = True
            except ClientError:
                time.sleep(60)
                continue
            try:
                s3.download_file(BUCKET, f"{PREFIX}/{RESULTS_JSON}",
                                 os.path.join("routes", RESULTS_JSON))
                print(f"[ok] DONE — pulled routes/{RESULTS_JSON}", flush=True)
            except ClientError as e:
                print(f"[ok] DONE, but pull {RESULTS_JSON} failed: {e}",
                      flush=True)
        if not done:
            print(f"[timeout] MAX_MIN reached, DONE not seen", flush=True)
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
