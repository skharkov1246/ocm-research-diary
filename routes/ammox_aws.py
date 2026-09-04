#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
routes/ammox_aws.py — leak-proof boto3 launcher for the Track L7 opener
(routes/ammox_mobi.py: propylene ammoxidation on a Bi-Mo-O Sohio-site proxy —
ddE# = barrier(acrolein over-oxidation) - barrier(allylic H abstraction) at
DFT / CASSCF / NEVPT2 with the mirage-detector gates).

Modeled on calc/fe_ni_s_aws.py + routes/andr_aws.py. Dead-man layering:
  1. bootstrap schedules `shutdown -P` FIRST (fires even if setup hangs);
     InstanceInitiatedShutdownBehavior=terminate -> a poweroff TERMINATES.
  2. the same watchdog is re-scheduled inside user-data before any real work.
  3. sandbox-side reaper: this script terminates the tagged instance(s) on exit.
IMDSv2 required (HttpTokens=required). Results stream to S3 as the job writes.

ONE INSTANCE. The pipeline is a chain (spins -> geom -> polish -> profile ->
merge) over two substrates that must share one spin state and one active space,
so it is NOT split across boxes; the stages themselves are resume-checkpointed,
and a relaunch on the same S3 prefix continues where the previous box died.

DEFAULTS
  region  us-west-2
  type    c7i.8xlarge  (32 vCPU, 64 GB — compute-bound geomeTRIC scan + CASSCF
                        on ~210 basis functions; the job is CPU-, not RAM-bound)
  S3      s3://alpha-o-results-097743207937/alpha-o/ammox/
  tags    Project=sandbox-agent  Name=ammox
  branch  claude/lukoil-norilsk-optimization-qbrqpk (env BRANCH)
  reaper  JOB_MAX_MIN=900 (15 h): 14 correlated scan points is a long chain

ENV
  AMMOX_ITYPE (c7i.8xlarge)  BRANCH  JOB_REGION (us-west-2)  JOB_MAX_MIN (900)
  DRYRUN (0)  and the science knobs, forwarded verbatim:
  AMMOX_STAGE_TIMEOUT  AMMOX_LEVEL  AMMOX_NPTS  AMMOX_SPIN  AMMOX_NDOCC
  AMMOX_NVIR  AMMOX_BASIS  PYSCF_MAX_MEMORY

Run:  python3 routes/ammox_aws.py              (do NOT run from CI; launches EC2)
      DRYRUN=1 python3 routes/ammox_aws.py
"""
import os
import sys
import time

import boto3
from botocore.exceptions import ClientError

REGION = os.environ.get("JOB_REGION", "us-west-2")
ACCOUNT = "097743207937"
BUCKET = f"alpha-o-results-{ACCOUNT}"          # us-east-1 bucket, cross-region OK
PREFIX = "alpha-o/ammox"
PROFILE = "alpha-o-instance-profile"
TAG_PROJECT = "sandbox-agent"
TAG_NAME = "ammox"
ITYPE = os.environ.get("AMMOX_ITYPE", "c7i.8xlarge")
REPO = "https://github.com/skharkov1246/ocm-research-diary"
BRANCH = os.environ.get("BRANCH") or os.environ.get("JOB_BRANCH") \
    or "claude/lukoil-norilsk-optimization-qbrqpk"
MAX_MIN = int(os.environ.get("JOB_MAX_MIN", "900"))
DRYRUN = os.environ.get("DRYRUN", "0") == "1" or "--dryrun" in sys.argv
DISK_GB = 80

# per-instance job budget: leave 20 min of the reaper cap for uploads + terminate
JOB_MIN = MAX_MIN - 20
FWD = ("AMMOX_STAGE_TIMEOUT", "AMMOX_LEVEL", "AMMOX_NPTS", "AMMOX_SPIN",
       "AMMOX_NDOCC", "AMMOX_NVIR", "AMMOX_BASIS", "PYSCF_MAX_MEMORY")

# ВАЖНО: env-креды нужны ЛАУНЧЕРУ (песочница вызывает EC2/S3 от scoped-юзера);
# на инстанс они не попадают — user-data кредов не содержит, там instance profile.
ec2 = boto3.client("ec2", region_name=REGION)
ssm = boto3.client("ssm", region_name=REGION)
s3 = boto3.client("s3", region_name="us-east-1")


def _fwd_exports():
    return "\n".join(f"export {k}={os.environ[k]!r}"
                     for k in FWD if os.environ.get(k))


def userdata():
    """user-data: watchdog shutdown FIRST, then clone, install, pull any prior
    checkpoints back from S3 (resume), run the full stage chain, stream every
    ammox_* JSON + the log to S3 continuously, poweroff."""
    exports = _fwd_exports()
    run_to = JOB_MIN - 15                       # minutes for the python job
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
python3 -c "import pyscf, geometric" || {{ aws s3 cp /tmp/boot.log $S3/setup_failed.log; poweroff; }}
mkdir -p /root/scratch
export PYSCF_TMPDIR=/root/scratch
export OMP_NUM_THREADS=$(nproc)
export PYSCF_MAX_MEMORY=${{PYSCF_MAX_MEMORY:-$(( $(getconf _PHYS_PAGES) * $(getconf PAGE_SIZE) / 1024 / 1024 * 6 / 10 ))}}
{exports}
# resume: pull whatever a previous box already computed for this prefix
aws s3 cp $S3/ routes/ --recursive --exclude '*' --include 'ammox_mobi_*.json' 2>/dev/null || true
sync_up() {{
  aws s3 cp routes/ $S3/ --recursive --exclude '*' --include 'ammox_mobi_*.json' >/dev/null 2>&1 || true
  aws s3 cp /root/repo/ammox.log $S3/ammox.log >/dev/null 2>&1 || true
}}
( while true; do sleep 120; sync_up; done ) &
echo "[aws][ammox] start" >> /root/repo/ammox.log
timeout {run_to}m python3 -u routes/ammox_mobi.py all >> /root/repo/ammox.log 2>&1
echo "[aws][ammox] DONE rc=$?" >> /root/repo/ammox.log
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


PULL = ["ammox_mobi_results.json",
        "ammox_mobi_spins.json",
        "ammox_mobi_c3h6_geom.json", "ammox_mobi_acrolein_geom.json",
        "ammox_mobi_c3h6_profile.json", "ammox_mobi_acrolein_profile.json"]


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
            print("[dryrun] OK — launch permitted" if "DryRunOperation" in str(e)
                  else f"[dryrun] FAILED: {e}", flush=True)
        return

    existing = my_instances()
    if existing:
        sys.exit(f"ammox instance(s) already exist: {existing}")

    iid = ec2.run_instances(**spec)["Instances"][0]["InstanceId"]
    print(f"launched {iid} -> routes/ammox_mobi_results.json", flush=True)

    t0 = time.time()
    try:
        while time.time() - t0 < MAX_MIN * 60:
            try:
                s3.head_object(Bucket=BUCKET, Key=f"{PREFIX}/DONE")
            except ClientError:
                time.sleep(60)
                continue
            for fn in PULL:
                try:
                    s3.download_file(BUCKET, f"{PREFIX}/{fn}", f"routes/{fn}")
                    print(f"[ok] pulled routes/{fn}", flush=True)
                except ClientError:
                    pass
            break
        else:
            print("[timeout] MAX_MIN reached without a DONE marker", flush=True)
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
