#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
calc/fes_cogef_aws.py — leak-proof boto3 launcher for the Track T3 gate
(calc/fes_cogef.py: COGEF stretch of a terminal Fe–S(thiolate) bond of the
[Fe4S4(SH)4]2- cubane — "does Fe–S break easier or harder than C–C under force?"
the screen-gate for "tailings as a mine in reverse", mechanochemical opening of
the sulfide lattice).

Modeled on calc/fe_ni_s_aws.py. ONE instance (no split). Dead-man layering:
  1. bootstrap schedules `shutdown -P` FIRST (fires even if setup hangs);
     InstanceInitiatedShutdownBehavior=terminate -> a poweroff TERMINATES.
  2. a watchdog `shutdown -P` is scheduled again inside user-data before any
     real work, as the very first command.
  3. sandbox-side reaper: this script terminates the tagged instance on exit.
IMDSv2 required (HttpTokens=required). Results stream to S3 as the job writes.

DEFAULTS
  region  us-west-2
  type    r7i.8xlarge  (32 vCPU, 256 GB; env FES_ITYPE)
  S3      s3://alpha-o-results-097743207937/alpha-o/fes-cogef/
  tags    Project=sandbox-agent  Name=fes-cogef
  branch  claude/lukoil-norilsk-optimization-qbrqpk (env BRANCH)

MODES
  default   one instance: the full COGEF scan + analysis
  DRYRUN=1  validate the launch (ec2 DryRun) without starting anything

ENV
  FES_ITYPE (r7i.8xlarge)  BRANCH  JOB_REGION (us-west-2)  JOB_MAX_MIN (600)
  DRYRUN (0)  and any FES_* knob forwarded to the science job

Run:  python3 calc/fes_cogef_aws.py            (do NOT run from CI; launches EC2)
      DRYRUN=1 python3 calc/fes_cogef_aws.py
"""
import os
import sys
import time

import boto3
from botocore.exceptions import ClientError

REGION = os.environ.get("JOB_REGION", "us-west-2")
ACCOUNT = "097743207937"
BUCKET = f"alpha-o-results-{ACCOUNT}"           # us-east-1 bucket, cross-region OK
PREFIX = "alpha-o/fes-cogef"
PROFILE = "alpha-o-instance-profile"
TAG_PROJECT = "sandbox-agent"
TAG_NAME = "fes-cogef"
ITYPE = os.environ.get("FES_ITYPE", "r7i.8xlarge")
REPO = "https://github.com/skharkov1246/ocm-research-diary"
BRANCH = os.environ.get("BRANCH") or os.environ.get("JOB_BRANCH") \
    or "claude/lukoil-norilsk-optimization-qbrqpk"
MAX_MIN = int(os.environ.get("JOB_MAX_MIN", "600"))
DRYRUN = os.environ.get("DRYRUN", "0") == "1" or "--dryrun" in sys.argv
DISK_GB = 80

# per-instance job budget: leave 20 min of the reaper cap for uploads + terminate
JOB_MIN = MAX_MIN - 20
# forward these science knobs verbatim if the caller set them
FWD = ("FES_STAGE_TIMEOUT", "FES_DMAX", "FES_STEP", "FES_RIGID", "FES_REF",
       "PYSCF_MAX_MEMORY")

# ВАЖНО: env-креды нужны ЛАУНЧЕРУ (песочница вызывает EC2/S3 от scoped-юзера);
# на инстанс они не попадают — user-data кредов не содержит, там instance profile.
ec2 = boto3.client("ec2", region_name=REGION)
ssm = boto3.client("ssm", region_name=REGION)
s3 = boto3.client("s3", region_name="us-east-1")


def _fwd_exports():
    """Shell export lines for forwarded FES_* knobs the caller set."""
    lines = []
    for k in FWD:
        v = os.environ.get(k)
        if v:
            lines.append(f"export {k}={v!r}")
    return "\n".join(lines)


def userdata():
    """user-data: watchdog shutdown FIRST, then clone, install, run the science
    job, stream JSON+log to S3, poweroff."""
    exports = _fwd_exports()
    run_to = JOB_MIN - 15                        # minutes for the python job
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
for i in 1 2 3; do timeout 25m python3 -m pip install -q numpy scipy pyscf geometric && break; sleep 30; done
for i in 1 2 3; do (cd /root && timeout 10m git clone --depth 1 -b {BRANCH} {REPO} repo) && break; sleep 30; done
cd /root/repo || {{ aws s3 cp /tmp/boot.log $S3/setup_failed.log; poweroff; }}
python3 -c "import pyscf" || {{ aws s3 cp /tmp/boot.log $S3/setup_failed.log; poweroff; }}
mkdir -p /root/scratch
export PYSCF_TMPDIR=/root/scratch
export OMP_NUM_THREADS=$(nproc)
export PYSCF_MAX_MEMORY=${{PYSCF_MAX_MEMORY:-$(( $(getconf _PHYS_PAGES) * $(getconf PAGE_SIZE) / 1024 / 1024 * 6 / 10 ))}}
{exports}
# resume: pull any prior checkpoint back down before running
aws s3 cp $S3/fes_cogef_results.json calc/fes_cogef_results.json 2>/dev/null || true
sync_up() {{
  aws s3 cp calc/fes_cogef_results.json $S3/fes_cogef_results.json >/dev/null 2>&1 || true
  aws s3 cp /root/repo/fes-cogef.log $S3/fes-cogef.log >/dev/null 2>&1 || true
}}
( while true; do sleep 120; sync_up; done ) &
echo "[aws][fes-cogef] start" >> /root/repo/fes-cogef.log
timeout {run_to}m python3 -u calc/fes_cogef.py >> /root/repo/fes-cogef.log 2>&1
echo "[aws][fes-cogef] DONE rc=$?" >> /root/repo/fes-cogef.log
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


def my_instances():
    r = ec2.describe_instances(Filters=[
        {"Name": "tag:Project", "Values": [TAG_PROJECT]},
        {"Name": "tag:Name", "Values": [TAG_NAME]},
        {"Name": "instance-state-name",
         "Values": ["pending", "running", "stopping", "stopped"]}])
    return [i["InstanceId"] for res in r["Reservations"] for i in res["Instances"]]


def run_spec(ami, subnet, ud):
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
        UserData=ud)


def main():
    print(f"branch={BRANCH} region={REGION} type={ITYPE} "
          f"bucket=s3://{BUCKET}/{PREFIX} reap<={MAX_MIN}min", flush=True)
    subnet = default_subnet()
    ami = latest_ami()
    spec = run_spec(ami, subnet, userdata())

    if DRYRUN:
        try:
            ec2.run_instances(DryRun=True, **spec)
            print("[dryrun] unexpected success", flush=True)
        except ClientError as e:
            ok = "DryRunOperation" in str(e)
            print("[dryrun] " + ("OK — launch permitted" if ok
                                 else f"FAILED: {e}"), flush=True)
        return

    existing = my_instances()
    if existing:
        sys.exit(f"fes-cogef instance(s) already exist: {existing}")

    iid = ec2.run_instances(**spec)["Instances"][0]["InstanceId"]
    print(f"launched {iid} -> fes_cogef_results.json", flush=True)

    t0 = time.time()
    done = False
    try:
        while time.time() - t0 < MAX_MIN * 60 and not done:
            try:
                s3.head_object(Bucket=BUCKET, Key=f"{PREFIX}/DONE")
                try:
                    s3.download_file(BUCKET, f"{PREFIX}/fes_cogef_results.json",
                                     "calc/fes_cogef_results.json")
                except ClientError:
                    pass
                print("[ok] DONE — pulled calc/fes_cogef_results.json", flush=True)
                done = True
            except ClientError:
                time.sleep(60)
        if not done:
            print("[timeout] MAX_MIN reached without DONE marker", flush=True)
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
