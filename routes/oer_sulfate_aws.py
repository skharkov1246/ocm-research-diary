#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
routes/oer_sulfate_aws.py — leak-proof лаунчер: routes/oer_sulfate_ladder.py
(гейт G1 кейса N1, OER на Ni в сульфатной кислой среде) на c7i.8xlarge
в us-west-2. Защита: terminate-on-shutdown + бортовой shutdown ПЕРВОЙ строкой
+ жнец; DONE-маркер отдельно; IMDSv2 (HttpTokens=required).

Режимы:
  python3 routes/oer_sulfate_aws.py            # один инстанс, все species
  DRYRUN=1 python3 routes/oer_sulfate_aws.py   # только DryRun-проверка прав
  OER_SPLIT=1 python3 routes/oer_sulfate_aws.py
      # ДВА инстанса параллельно: #1 — сульфатная ветка (suffix _sulf),
      # #2 — контрольная ветка (suffix _ctrl); каждый пишет СВОЙ
      # oer_sulfate_ladder_results{suffix}.json — в S3 не сталкиваются.

env-ручки: OER_ITYPE (default c7i.8xlarge), BRANCH
(default claude/lukoil-norilsk-optimization-qbrqpk), OER_MAX_MIN (default 720),
OER_STAGE_TIMEOUT (проброс в научный скрипт, default 5400), OER_BASIS, OER_XC.
"""
import os
import sys
import time

import boto3
from botocore.exceptions import ClientError

REGION = os.environ.get("OER_REGION", "us-west-2")
ACCOUNT = "097743207937"
BUCKET = f"alpha-o-results-{ACCOUNT}"          # us-east-1, cross-region OK
PREFIX = "alpha-o/oer-sulfate"
PROFILE = "alpha-o-instance-profile"
TAG_PROJECT = "sandbox-agent"
TAG_NAME = "oer-sulfate"
ITYPE = os.environ.get("OER_ITYPE", "c7i.8xlarge")   # 32 vCPU
REPO = "https://github.com/skharkov1246/ocm-research-diary"
BRANCH = os.environ.get("BRANCH", "claude/lukoil-norilsk-optimization-qbrqpk")
MAX_MIN = int(os.environ.get("OER_MAX_MIN", "720"))
STAGE_TIMEOUT = int(os.environ.get("OER_STAGE_TIMEOUT", "5400"))
OER_BASIS = os.environ.get("OER_BASIS", "def2-svp")
OER_XC = os.environ.get("OER_XC", "pbe,pbe")
SPLIT = os.environ.get("OER_SPLIT", "") == "1"
DRYRUN = os.environ.get("DRYRUN", "") == "1" or "--dryrun" in sys.argv
JOB_MIN = MAX_MIN - 20          # бортовой watchdog чуть раньше жнеца
SCRIPT_TO_MIN = JOB_MIN - 40    # timeout самого python-запуска
DISK_GB = 60

SULF_SPECIES = "sulf_bare,sulf_oh,sulf_o,sulf_ooh,h2o,h2"
CTRL_SPECIES = "ctrl_bare,ctrl_oh,ctrl_o,ctrl_ooh,h2o,h2"
ALL_SPECIES = "sulf_bare,sulf_oh,sulf_o,sulf_ooh,ctrl_bare,ctrl_oh,ctrl_o,ctrl_ooh,h2o,h2"

# ВАЖНО: env-креды нужны ЛАУНЧЕРУ (песочница вызывает EC2/S3 от scoped-юзера);
# на инстанс они не попадают — user-data кредов не содержит, там instance profile.
ec2 = boto3.client("ec2", region_name=REGION)
ssm = boto3.client("ssm", region_name=REGION)
s3 = boto3.client("s3", region_name="us-east-1")


def userdata(species_csv, suffix):
    """user-data одного инстанса: watchdog ПЕРВЫМ, клон ветки, pip-зависимости,
    резюм results-JSON из S3, запуск научного скрипта, выгрузка, poweroff."""
    rj = f"oer_sulfate_ladder_results{suffix}.json"
    log = f"oer-sulfate{suffix}.log"
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
export OER_SPECIES={species_csv}
export OER_JSON_SUFFIX={suffix}
export OER_STAGE_TIMEOUT={STAGE_TIMEOUT}
export OER_BASIS={OER_BASIS}
export OER_XC={OER_XC}
sync_up() {{ aws s3 cp routes/{rj} $S3/{rj} >/dev/null 2>&1 || true; aws s3 cp /root/repo/{log} $S3/{log} >/dev/null 2>&1 || true; }}
( while true; do sleep 120; sync_up; done ) &
echo "[aws][oer-sulfate{suffix}] start species={species_csv}" >> /root/repo/{log}
timeout {SCRIPT_TO_MIN}m python3 -u routes/oer_sulfate_ladder.py >> /root/repo/{log} 2>&1
echo "[aws][oer-sulfate{suffix}] ALL DONE" >> /root/repo/{log}
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


def launch_spec(ami, subnet, name, species_csv, suffix):
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
        UserData=userdata(species_csv, suffix))


def main():
    if SPLIT:
        jobs = [(TAG_NAME + "-sulf", SULF_SPECIES, "_sulf"),
                (TAG_NAME + "-ctrl", CTRL_SPECIES, "_ctrl")]
    else:
        jobs = [(TAG_NAME, ALL_SPECIES, "")]
    subnet = default_subnet()
    ami = latest_ami()
    print(f"branch={BRANCH} region={REGION} type={ITYPE} subnet={subnet} "
          f"bucket=s3://{BUCKET}/{PREFIX} split={SPLIT} reap<={MAX_MIN}min",
          flush=True)

    if DRYRUN:
        for name, species_csv, suffix in jobs:
            try:
                ec2.run_instances(DryRun=True,
                                  **launch_spec(ami, subnet, name,
                                                species_csv, suffix))
                print(f"[dryrun] {name}: unexpected success", flush=True)
            except ClientError as e:
                print(f"[dryrun] {name}: "
                      + ("OK — launch permitted" if "DryRunOperation" in str(e)
                         else f"FAILED: {e}"), flush=True)
        return

    existing = my_instances()
    if existing:
        sys.exit(f"oer-sulfate instance(s) already exist: {existing}")

    iids = []
    for name, species_csv, suffix in jobs:
        iid = ec2.run_instances(**launch_spec(ami, subnet, name, species_csv,
                                              suffix))["Instances"][0]["InstanceId"]
        iids.append(iid)
        print(f"launched {name}: {iid}", flush=True)

    done_keys = {f"{PREFIX}/DONE{suffix}": f"oer_sulfate_ladder_results{suffix}.json"
                 for _n, _s, suffix in jobs}
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
