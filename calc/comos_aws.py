#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
calc/comos_aws.py — leak-proof boto3 launcher for the Track L6 opener
(calc/comos_hds.py: [M3S9]2- MoS2 edge-motif trimer, unpromoted / Co-promoted /
Ni-promoted — BS-UNO spin structure + sulfur-vacancy HDS descriptor + the
PBE-vs-PBE0 "where DFT lies" block, for hydrotreating catalyst passports).

Modeled on calc/fe_ni_s_aws.py. Dead-man layering:
  1. bootstrap schedules `shutdown -P` FIRST (fires even if setup hangs);
     InstanceInitiatedShutdownBehavior=terminate -> a poweroff TERMINATES.
  2. the same watchdog is re-scheduled inside user-data before any real work.
  3. sandbox-side reaper: this script terminates the tagged instance(s) on exit.
IMDSv2 required (HttpTokens=required). Results stream to S3 as the job writes.

DEFAULTS
  region  us-west-2
  type    r7i.8xlarge  (32 vCPU, 256 GB — many UKS sectors on ~240 basis functions)
  S3      s3://alpha-o-results-097743207937/alpha-o/comos/
  tags    Project=sandbox-agent  Name=comos
  branch  claude/lukoil-norilsk-optimization-qbrqpk (env BRANCH)

MODES
  default          ONE instance: all three compositions in sequence -> results.json
  COMOS_SPLIT=1    THREE instances, one per promoter:
                     #1 COMOS_ONLY=mo -> comos_hds_results_mo.json
                     #2 COMOS_ONLY=co -> comos_hds_results_co.json
                     #3 COMOS_ONLY=ni -> comos_hds_results_ni.json
                   (each box also runs its own PBE0 block on its own composition;
                    the promoter ranking is assembled from the three JSONs)
  DRYRUN=1         validate the launch (ec2 DryRun) without starting anything

ENV
  COMOS_ITYPE (r7i.8xlarge)  BRANCH  JOB_REGION (us-west-2)  JOB_MAX_MIN (600)
  COMOS_SPLIT (0)  DRYRUN (0)  and any COMOS_* science knob (forwarded verbatim)

Run:  python3 calc/comos_aws.py               (do NOT run from CI; launches EC2)
      DRYRUN=1 python3 calc/comos_aws.py
      COMOS_SPLIT=1 python3 calc/comos_aws.py
"""
import os
import sys
import time

import boto3
from botocore.exceptions import ClientError

REGION = os.environ.get("JOB_REGION", "us-west-2")
ACCOUNT = "097743207937"
BUCKET = f"alpha-o-results-{ACCOUNT}"          # us-east-1 bucket, cross-region OK
PREFIX = "alpha-o/comos"
PROFILE = "alpha-o-instance-profile"
TAG_PROJECT = "sandbox-agent"
TAG_NAME = "comos"
ITYPE = os.environ.get("COMOS_ITYPE", "r7i.8xlarge")
REPO = "https://github.com/skharkov1246/ocm-research-diary"
BRANCH = os.environ.get("BRANCH") or os.environ.get("JOB_BRANCH") \
    or "claude/lukoil-norilsk-optimization-qbrqpk"
MAX_MIN = int(os.environ.get("JOB_MAX_MIN", "600"))
SPLIT = os.environ.get("COMOS_SPLIT", "0") == "1"
DRYRUN = os.environ.get("DRYRUN", "0") == "1" or "--dryrun" in sys.argv
DISK_GB = 80

# per-instance job budget: leave 20 min of the reaper cap for uploads + terminate
JOB_MIN = MAX_MIN - 20
# forward these science knobs verbatim if the caller set them
FWD = ("COMOS_STAGE_TIMEOUT", "COMOS_MF", "COMOS_SKIP_STAGE4", "PYSCF_MAX_MEMORY")

# ВАЖНО: env-креды нужны ЛАУНЧЕРУ (песочница вызывает EC2/S3 от scoped-юзера);
# на инстанс они не попадают — user-data кредов не содержит, там instance profile.
ec2 = boto3.client("ec2", region_name=REGION)
ssm = boto3.client("ssm", region_name=REGION)
s3 = boto3.client("s3", region_name="us-east-1")


def _fwd_exports(extra):
    """Shell export lines for forwarded COMOS_* knobs + per-run `extra`."""
    lines = []
    for k in FWD:
        v = os.environ.get(k)
        if v:
            lines.append(f"export {k}={v!r}")
    for k, v in extra.items():
        lines.append(f"export {k}={v!r}")
    return "\n".join(lines)


def userdata(only, json_suffix, total_budget_s):
    """user-data: watchdog shutdown FIRST, then clone, install, run the science
    job, stream JSON+log to S3, poweroff. `only` selects the composition for the
    split mode ('' | 'mo' | 'co' | 'ni'); json_suffix distinguishes S3 output."""
    exports = _fwd_exports({"COMOS_ONLY": only,
                            "COMOS_TOTAL_BUDGET": str(total_budget_s)})
    run_to = JOB_MIN - 15                       # minutes for the python job
    return f"""#!/bin/bash
shutdown -P +{JOB_MIN}
exec > /tmp/boot.log 2>&1
set -x
export PATH=$PATH:/usr/local/bin:/root/.local/bin
S3=s3://{BUCKET}/{PREFIX}
SUF={json_suffix}
export AWS_DEFAULT_REGION=us-east-1
for i in 1 2 3; do timeout 10m dnf install -y python3-pip git gcc gcc-c++ && break; sleep 30; done
for i in 1 2 3; do timeout 10m python3 -m pip install -q awscli && break; sleep 30; done
( while true; do aws s3 cp /tmp/boot.log $S3/boot$SUF.log >/dev/null 2>&1; sleep 60; done ) &
for i in 1 2 3; do timeout 25m python3 -m pip install -q numpy scipy pyscf && break; sleep 30; done
for i in 1 2 3; do (cd /root && timeout 10m git clone --depth 1 -b {BRANCH} {REPO} repo) && break; sleep 30; done
cd /root/repo || {{ aws s3 cp /tmp/boot.log $S3/setup_failed$SUF.log; poweroff; }}
python3 -c "import pyscf" || {{ aws s3 cp /tmp/boot.log $S3/setup_failed$SUF.log; poweroff; }}
mkdir -p /root/scratch
export PYSCF_TMPDIR=/root/scratch
export OMP_NUM_THREADS=$(nproc)
export PYSCF_MAX_MEMORY=${{PYSCF_MAX_MEMORY:-$(( $(getconf _PHYS_PAGES) * $(getconf PAGE_SIZE) / 1024 / 1024 * 6 / 10 ))}}
{exports}
# resume: pull any prior checkpoint for this half back down before running
aws s3 cp $S3/comos_hds_results$SUF.json calc/comos_hds_results.json 2>/dev/null || true
sync_up() {{
  cp calc/comos_hds_results.json calc/comos_hds_results$SUF.json 2>/dev/null || true
  aws s3 cp calc/comos_hds_results$SUF.json $S3/comos_hds_results$SUF.json >/dev/null 2>&1 || true
  aws s3 cp /root/repo/comos$SUF.log $S3/comos$SUF.log >/dev/null 2>&1 || true
  aws s3 cp calc/ $S3/ --recursive --exclude '*' --include 'comos_*.xyz' >/dev/null 2>&1 || true
}}
( while true; do sleep 120; sync_up; done ) &
echo "[aws][comos] start only='{only}' suffix='$SUF'" >> /root/repo/comos$SUF.log
timeout {run_to}m python3 -u calc/comos_hds.py >> /root/repo/comos$SUF.log 2>&1
echo "[aws][comos] DONE rc=$?" >> /root/repo/comos$SUF.log
sync_up
echo done > /tmp/DONE && aws s3 cp /tmp/DONE $S3/DONE$SUF
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


def run_spec(ami, subnet, ud, name):
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
        UserData=ud)


def jobs():
    """(only, json_suffix, done_marker) tuples for the chosen mode."""
    if SPLIT:
        return [(c, f"_{c}", f"DONE_{c}") for c in ("mo", "co", "ni")]
    return [("", "", "DONE")]


def main():
    budget = (JOB_MIN - 15) * 60 - 300     # python budget minus upload headroom
    print(f"branch={BRANCH} region={REGION} type={ITYPE} split={SPLIT} "
          f"bucket=s3://{BUCKET}/{PREFIX} reap<={MAX_MIN}min", flush=True)
    subnet = default_subnet()
    ami = latest_ami()
    specs = [(only, suf, done, run_spec(
                ami, subnet, userdata(only, suf, budget), TAG_NAME))
             for only, suf, done in jobs()]

    if DRYRUN:
        for only, suf, _done, spec in specs:
            try:
                ec2.run_instances(DryRun=True, **spec)
                print(f"[dryrun] only='{only or 'all'}' suffix='{suf}': "
                      "unexpected success", flush=True)
            except ClientError as e:
                ok = "DryRunOperation" in str(e)
                print(f"[dryrun] only='{only or 'all'}' suffix='{suf}': "
                      + ("OK — launch permitted" if ok else f"FAILED: {e}"),
                      flush=True)
        return

    existing = my_instances()
    if existing:
        sys.exit(f"comos instance(s) already exist: {existing}")

    launched = []
    for only, suf, done, spec in specs:
        iid = ec2.run_instances(**spec)["Instances"][0]["InstanceId"]
        launched.append((iid, suf, done))
        print(f"launched {iid}  only='{only or 'all'}' -> "
              f"comos_hds_results{suf}.json", flush=True)

    pending = {done for _iid, _suf, done in launched}
    t0 = time.time()
    try:
        while time.time() - t0 < MAX_MIN * 60 and pending:
            for _iid, suf, done in launched:
                if done not in pending:
                    continue
                try:
                    s3.head_object(Bucket=BUCKET, Key=f"{PREFIX}/{done}")
                    fn = f"comos_hds_results{suf}.json"
                    try:
                        s3.download_file(BUCKET, f"{PREFIX}/{fn}", f"calc/{fn}")
                    except ClientError:
                        pass
                    print(f"[ok] {done} — pulled calc/{fn}", flush=True)
                    pending.discard(done)
                except ClientError:
                    pass
            if pending:
                time.sleep(60)
        if pending:
            print(f"[timeout] MAX_MIN reached, still pending: {pending}",
                  flush=True)
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
