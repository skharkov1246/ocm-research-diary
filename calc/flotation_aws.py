#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
calc/flotation_aws.py — leak-proof boto3 launcher for the Track N9 gate
(calc/flotation_adsorption.py: adsorption of flotation reagents on a
PENTLANDITE vs a PYRRHOTITE motif — the electronic gate of "take the sulfur out
BEFORE the stack", i.e. before the Sulfur Programme has to catch it as SO2).

Modeled on calc/fe_ni_s_aws.py. Dead-man layering:
  1. bootstrap schedules `shutdown -P` FIRST (fires even if setup hangs);
     InstanceInitiatedShutdownBehavior=terminate -> a poweroff TERMINATES.
  2. a watchdog `shutdown -P` is scheduled again inside user-data before any
     real work, as the very first command.
  3. sandbox-side reaper: this script terminates the tagged instance(s) on exit.
IMDSv2 required (HttpTokens=required). Results stream to S3 as the job writes.

DEFAULTS
  region  us-west-2
  type    r7i.8xlarge  (32 vCPU, 256 GB — the complexes reach ~380 basis
          functions, and every point is a broken-symmetry UKS SCF)
  S3      s3://alpha-o-results-097743207937/alpha-o/flotation/
  tags    Project=sandbox-agent  Name=flotation
  branch  claude/lukoil-norilsk-optimization-qbrqpk (env BRANCH)

MODES
  default          one instance: all 9 substrate-adsorbate-site pairs + the
                   PBE0//PBE sign check on both descriptors  (~9 h)
  FLOT_SPLIT=1     two instances (recommended):
                     #1 suffix _coll  xanthate + OH- on all three sites,
                        then the PBE0 sign check of the COLLECTOR descriptor
                     #2 suffix _depr  ethylenediamine (DETA proxy) on all three
                        sites, then the PBE0 sign check of the DEPRESSANT
                        descriptor
                   Each half owns the PBE data its own hybrid check needs, so
                   the two boxes are fully independent; the merged JSONs carry
                   a hybrid-verified sign for BOTH descriptors (the split spec
                   only asked #2 to run a PBE0 push — #1 pushing its own
                   collector pair is the cheap extension that keeps the
                   "look at the sign on a hybrid" rule for both levers).
  DRYRUN=1         validate the launch (ec2 DryRun) without starting anything

ENV
  FLOT_ITYPE (r7i.8xlarge)  BRANCH  JOB_REGION (us-west-2)  JOB_MAX_MIN (600)
  FLOT_SPLIT (0)  DRYRUN (0)  and any FLOT_* science knob forwarded verbatim

Run:  python3 calc/flotation_aws.py            (do NOT run from CI; launches EC2)
      DRYRUN=1 python3 calc/flotation_aws.py
"""
import os
import sys
import time

import boto3
from botocore.exceptions import ClientError

REGION = os.environ.get("JOB_REGION", "us-west-2")
ACCOUNT = "097743207937"
BUCKET = f"alpha-o-results-{ACCOUNT}"          # us-east-1 bucket, cross-region OK
PREFIX = "alpha-o/flotation"
PROFILE = "alpha-o-instance-profile"
TAG_PROJECT = "sandbox-agent"
TAG_NAME = "flotation"
ITYPE = os.environ.get("FLOT_ITYPE", "r7i.8xlarge")
REPO = "https://github.com/skharkov1246/ocm-research-diary"
BRANCH = os.environ.get("BRANCH") or os.environ.get("JOB_BRANCH") \
    or "claude/lukoil-norilsk-optimization-qbrqpk"
MAX_MIN = int(os.environ.get("JOB_MAX_MIN", "600"))
SPLIT = os.environ.get("FLOT_SPLIT", "0") == "1"
DRYRUN = os.environ.get("DRYRUN", "0") == "1" or "--dryrun" in sys.argv
DISK_GB = 80

# per-instance job budget: leave 20 min of the reaper cap for uploads + terminate
JOB_MIN = MAX_MIN - 20
# forward these science knobs verbatim if the caller set them
FWD = ("FLOT_STAGE_TIMEOUT", "FLOT_MODE", "FLOT_SCAN_N", "FLOT_SCAN_DR",
       "FLOT_ADS_RELAX", "PYSCF_MAX_MEMORY")

# ВАЖНО: env-креды нужны ЛАУНЧЕРУ (песочница вызывает EC2/S3 от scoped-юзера);
# на инстанс они не попадают — user-data кредов не содержит, там instance profile.
ec2 = boto3.client("ec2", region_name=REGION)
ssm = boto3.client("ssm", region_name=REGION)
s3 = boto3.client("s3", region_name="us-east-1")


def _fwd_exports(extra):
    """Shell export lines for forwarded FLOT_* knobs + per-run `extra`."""
    lines = []
    for k in FWD:
        v = os.environ.get(k)
        if v:
            lines.append(f"export {k}={v!r}")
    for k, v in extra.items():
        lines.append(f"export {k}={v!r}")
    return "\n".join(lines)


def userdata(only, pbe0_set, json_suffix, total_budget_s):
    """user-data: watchdog shutdown FIRST, then clone, install, run the science
    job, stream JSON+log to S3, poweroff. `only` selects the adsorbate subset
    for the split mode ('' | 'xanthate,oh' | 'en'); json_suffix distinguishes
    the S3 output of the two halves."""
    exports = _fwd_exports({"FLOT_ONLY": only,
                            "FLOT_MF": "pbe",
                            "FLOT_PBE0_SET": pbe0_set,
                            "FLOT_TOTAL_BUDGET": str(total_budget_s)})
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
for i in 1 2 3; do timeout 25m python3 -m pip install -q numpy scipy pyscf geometric && break; sleep 30; done
for i in 1 2 3; do (cd /root && timeout 10m git clone --depth 1 -b {BRANCH} {REPO} repo) && break; sleep 30; done
cd /root/repo || {{ aws s3 cp /tmp/boot.log $S3/setup_failed$SUF.log; poweroff; }}
python3 -c "import pyscf" || {{ aws s3 cp /tmp/boot.log $S3/setup_failed$SUF.log; poweroff; }}
mkdir -p /root/scratch
export PYSCF_TMPDIR=/root/scratch
export OMP_NUM_THREADS=$(nproc)
export PYSCF_MAX_MEMORY=${{PYSCF_MAX_MEMORY:-$(( $(getconf _PHYS_PAGES) * $(getconf PAGE_SIZE) / 1024 / 1024 * 6 / 10 ))}}
{exports}
# resume: pull this half's prior checkpoints (JSON + densities) back down first
aws s3 cp $S3/flotation_adsorption_results$SUF.json calc/flotation_adsorption_results.json 2>/dev/null || true
aws s3 cp --recursive --exclude '*' --include 'flotation_sub_*.npz' --include 'flotation_dm_*.npz' --include 'flotation_geom_*.npz' $S3/npz$SUF calc/ 2>/dev/null || true
sync_up() {{
  cp calc/flotation_adsorption_results.json calc/flotation_adsorption_results$SUF.json 2>/dev/null || true
  aws s3 cp calc/flotation_adsorption_results$SUF.json $S3/flotation_adsorption_results$SUF.json >/dev/null 2>&1 || true
  aws s3 cp /root/repo/flotation$SUF.log $S3/flotation$SUF.log >/dev/null 2>&1 || true
  aws s3 cp calc/flotation_adsorption_geoms.xyz $S3/flotation_adsorption_geoms$SUF.xyz >/dev/null 2>&1 || true
  for f in calc/flotation_sub_*.npz calc/flotation_dm_*.npz calc/flotation_geom_*.npz; do
    [ -f "$f" ] && aws s3 cp "$f" $S3/npz$SUF/$(basename $f) >/dev/null 2>&1 || true
  done
}}
( while true; do sleep 180; sync_up; done ) &
echo "[aws][flotation] start only='{only}' pbe0_set='{pbe0_set}' suffix='$SUF'" >> /root/repo/flotation$SUF.log
timeout {run_to}m python3 -u calc/flotation_adsorption.py >> /root/repo/flotation$SUF.log 2>&1
echo "[aws][flotation] DONE rc=$?" >> /root/repo/flotation$SUF.log
sync_up
for f in calc/flotation_sub_*.xyz; do aws s3 cp "$f" $S3/$(basename $f) >/dev/null 2>&1 || true; done
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
    """(only, pbe0_set, json_suffix, done_marker) tuples for the chosen mode."""
    if SPLIT:
        return [("xanthate,oh", "collector", "_coll", "DONE_coll"),
                ("en", "depressant", "_depr", "DONE_depr")]
    return [("", "both", "", "DONE")]


def main():
    budget = (JOB_MIN - 15) * 60 - 300     # python budget minus upload headroom
    print(f"branch={BRANCH} region={REGION} type={ITYPE} split={SPLIT} "
          f"bucket=s3://{BUCKET}/{PREFIX} reap<={MAX_MIN}min", flush=True)
    subnet = default_subnet()
    ami = latest_ami()
    specs = [(only, suf, done, run_spec(
                ami, subnet, userdata(only, p0, suf, budget), TAG_NAME))
             for only, p0, suf, done in jobs()]

    if DRYRUN:
        for only, suf, _done, spec in specs:
            try:
                ec2.run_instances(DryRun=True, **spec)
                print(f"[dryrun] only='{only or 'all pairs'}' suffix='{suf}': "
                      "unexpected success", flush=True)
            except ClientError as e:
                ok = "DryRunOperation" in str(e)
                print(f"[dryrun] only='{only or 'all pairs'}' suffix='{suf}': "
                      + ("OK — launch permitted" if ok else f"FAILED: {e}"),
                      flush=True)
        return

    existing = my_instances()
    if existing:
        sys.exit(f"flotation instance(s) already exist: {existing}")

    launched = []
    for only, suf, done, spec in specs:
        iid = ec2.run_instances(**spec)["Instances"][0]["InstanceId"]
        launched.append((iid, suf, done))
        print(f"launched {iid}  only='{only or 'all pairs'}' -> "
              f"flotation_adsorption_results{suf}.json", flush=True)

    pending = {done for _iid, _suf, done in launched}
    t0 = time.time()
    try:
        while time.time() - t0 < MAX_MIN * 60 and pending:
            for _iid, suf, done in launched:
                if done not in pending:
                    continue
                try:
                    s3.head_object(Bucket=BUCKET, Key=f"{PREFIX}/{done}")
                    fn = f"flotation_adsorption_results{suf}.json"
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
