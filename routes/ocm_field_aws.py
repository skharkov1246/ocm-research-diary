#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
routes/ocm_field_aws.py — leak-proof прогон Этапа 21-F (изобретение I5):
полевая поляризуемость селективности dΔΔE‡/dF на Cr-центре, NEVPT2.
ocm_field_dde.py field→readout на r7i.8xlarge в us-west-2.
Геометрии из гита (Cr-профили Этапа 17/18), скан = чистые синглпойнты;
поточечный чекпойнт ocm_field_scan.json тянется с S3 при буте (урок msa-v7).
Защита: terminate-on-shutdown + бортовой shutdown + жнец; DONE-маркер отдельно.

Запуск: python3 routes/ocm_field_aws.py [--dryrun]
"""
import os
import subprocess
import sys
import time

import boto3
from botocore.exceptions import ClientError

REGION = os.environ.get("FLD_REGION", "us-west-2")
ACCOUNT = "097743207937"
BUCKET = f"alpha-o-results-{ACCOUNT}"          # us-east-1, cross-region OK
PREFIX = "alpha-o/ocm-field"
PROFILE = "alpha-o-instance-profile"
TAG_PROJECT = "ocm-agent"
TAG_NAME = "ocm-field"
ITYPE = os.environ.get("FLD_ITYPE", "r7i.8xlarge")
REPO = "https://github.com/skharkov1246/ocm-research-diary"
BRANCH = os.environ.get("FLD_BRANCH") or subprocess.run(
    ["git", "rev-parse", "--abbrev-ref", "HEAD"], capture_output=True,
    text=True).stdout.strip()
MAX_MIN = int(os.environ.get("FLD_MAX_MIN", "840"))
STAGES = os.environ.get("FLD_STAGES", "field readout")
_OUT = {"field": "", "readout": "routes/ocm_field_dde_results.json",
        "field2": "", "readout2": "routes/ocm_field_dde2_results.json"}
RM_FILES = " ".join(f for f in (_OUT[s] for s in STAGES.split()) if f) or "/dev/null"
JOB_MIN = MAX_MIN - 20
STAGE_TO = JOB_MIN - 30
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
for i in 1 2 3; do timeout 25m python3 -m pip install -q numpy scipy pyscf && break; sleep 30; done
for i in 1 2 3; do (cd /root && timeout 10m git clone --depth 1 -b {BRANCH} {REPO} repo) && break; sleep 30; done
cd /root/repo || {{ aws s3 cp /tmp/boot.log $S3/setup_failed.log; poweroff; }}
python3 -c "import pyscf" || {{ aws s3 cp /tmp/boot.log $S3/setup_failed.log; poweroff; }}
rm -f {RM_FILES}
aws s3 cp $S3/ocm_field_scan.json routes/ocm_field_scan.json || true
aws s3 cp $S3/ocm_field_scan2.json routes/ocm_field_scan2.json || true
mkdir -p /root/scratch
export PYSCF_TMPDIR=/root/scratch
export OMP_NUM_THREADS=$(nproc)
sync_up() {{ aws s3 cp routes/ $S3/ --recursive --exclude '*' --include 'ocm_field_*.json' >/dev/null 2>&1 || true; aws s3 cp /root/repo/fld.log $S3/fld.log >/dev/null 2>&1 || true; }}
( while true; do sleep 120; sync_up; done ) &
for st in {STAGES}; do
  echo "[aws][fld] stage $st" >> /root/repo/fld.log
  timeout {STAGE_TO}m python3 -u routes/ocm_field_dde.py $st >> /root/repo/fld.log 2>&1
  sync_up
done
echo "[aws][fld] ALL DONE" >> /root/repo/fld.log
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
        sys.exit(f"ocm-field instance(s) already exist: {existing}")
    iid = ec2.run_instances(**common)["Instances"][0]["InstanceId"]
    print(f"launched {iid}", flush=True)
    t0 = time.time()
    try:
        while time.time() - t0 < MAX_MIN * 60:
            try:
                s3.head_object(Bucket=BUCKET, Key=f"{PREFIX}/DONE")
                for f in ("ocm_field_scan.json", "ocm_field_scan2.json",
                          "ocm_field_dde_results.json",
                          "ocm_field_dde2_results.json"):
                    try:
                        s3.download_file(BUCKET, f"{PREFIX}/{f}", f"routes/{f}")
                    except ClientError:
                        pass
                print("[ok] field sweep DONE — results pulled", flush=True)
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
