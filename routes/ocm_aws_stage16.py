#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
routes/ocm_aws_stage16.py — leak-proof EC2-оркестратор Этапа 16 OCM:
эмбеддед-кластер [O=Mn(–O–Na)(μ-O)₂W(=O)(OH)] + CH₄/C₂H₆, полный пайплайн
spins → geom → polish → profile → refine → merge2 на большом инстансе.

Паттерн — calc/aws_orchestrate.py (alpha-O), та же тройная защита от утечек:
  1. InstanceInitiatedShutdownBehavior=terminate — любой poweroff терминирует;
  2. on-instance watchdog `shutdown -P +N` ставится ПЕРВОЙ строкой user-data;
  3. внешний жнец в песочнице гасит инстанс по тегу к дедлайну (finally).
IMDSv2 required; креды на инстанс не попадают (S3 — через instance profile);
песочница читает прогресс/результаты из s3://<bucket>/alpha-o/ocm16/.

Запуск из песочницы:  python3 routes/ocm_aws_stage16.py [--dryrun]
Env-переопределения: OCM16_ITYPE (r7i.8xlarge), OCM16_MAX_MIN (750),
OCM16_BRANCH (текущая ветка), OCM16_FORCE=1 (игнорировать чужие инстансы нельзя —
жнец фильтрует по Name=ocm-stage16, чужой alpha-o не трогает).
"""
import json
import os
import subprocess
import sys
import time

import boto3
from botocore.exceptions import ClientError

REGION = os.environ.get("AWS_DEFAULT_REGION", "us-east-1")
ACCOUNT = "097743207937"
BUCKET = f"alpha-o-results-{ACCOUNT}"
PREFIX = "alpha-o/ocm16"
PROFILE = "alpha-o-instance-profile"
TAG_PROJECT = "sandbox-agent"          # обязателен: политика terminate по тегу
TAG_NAME = "ocm-stage16"               # жнец фильтрует ТОЛЬКО свои инстансы
ITYPE = os.environ.get("OCM16_ITYPE", "r7i.8xlarge")
REPO = "https://github.com/skharkov1246/ocm-research-diary"
BRANCH = os.environ.get("OCM16_BRANCH") or subprocess.run(
    ["git", "rev-parse", "--abbrev-ref", "HEAD"], capture_output=True,
    text=True).stdout.strip()
MAX_MIN = int(os.environ.get("OCM16_MAX_MIN", "750"))
JOB_MIN = MAX_MIN - 20
JOB_TIMEOUT = JOB_MIN - 20
DISK_GB = 100
FINAL_KEY = f"{PREFIX}/ocm_mnw_emb_final.json"

for var in ("AWS_ACCESS_KEY_ID", "AWS_SECRET_ACCESS_KEY", "AWS_SESSION_TOKEN"):
    os.environ.pop(var, None)          # битые env-креды → ~/.aws/credentials

ec2 = boto3.client("ec2", region_name=REGION)
ssm = boto3.client("ssm", region_name=REGION)
s3 = boto3.client("s3", region_name=REGION)

USERDATA = f"""#!/bin/bash
shutdown -P +{JOB_MIN}
exec > /tmp/boot.log 2>&1
set -x
export PATH=$PATH:/usr/local/bin:/root/.local/bin
S3=s3://{BUCKET}/{PREFIX}
for i in 1 2 3; do timeout 10m dnf install -y python3-pip git gcc gcc-c++ && break; sleep 30; done
for i in 1 2 3; do timeout 10m python3 -m pip install -q awscli && break; sleep 30; done
( while true; do aws s3 cp /tmp/boot.log $S3/boot.log >/dev/null 2>&1; sleep 60; done ) &
for i in 1 2 3; do timeout 25m python3 -m pip install -q numpy scipy pyscf geometric && break; sleep 30; done
for i in 1 2 3; do (cd /root && timeout 10m git clone --depth 1 -b {BRANCH} {REPO} repo) && break; sleep 30; done
cd /root/repo || {{ aws s3 cp /tmp/boot.log $S3/setup_failed.log; poweroff; }}
python3 -c "import pyscf, geometric" || {{ aws s3 cp /tmp/boot.log $S3/setup_failed.log; poweroff; }}
mkdir -p /root/scratch
export PYSCF_TMPDIR=/root/scratch
export PYSCF_MAX_MEMORY=$(( $(getconf _PHYS_PAGES) * $(getconf PAGE_SIZE) / 1024 / 1024 * 6 / 10 ))
export OCM_MODEL=embedded
S3=s3://{BUCKET}/{PREFIX}
sync_up() {{ aws s3 cp routes/ $S3/ --recursive --exclude '*' --include 'ocm_mnw_emb_*' >/dev/null 2>&1 || true; aws s3 cp /root/repo/emb_run.log $S3/emb_run.log >/dev/null 2>&1 || true; }}
( while true; do sleep 180; sync_up; done ) &
run() {{ OMP_NUM_THREADS=$2 OCM_SPIN=${{BEST:-5}} timeout {JOB_TIMEOUT}m python3 -u routes/ocm_mnw_selectivity.py $1 $3 >> emb_run.log 2>&1; }}
echo "[aws] spins ladder" >> emb_run.log
OMP_NUM_THREADS=32 timeout 120m python3 -u routes/ocm_mnw_selectivity.py spins >> emb_run.log 2>&1
BEST=$(python3 -c "import json;print(json.load(open('routes/ocm_mnw_emb_spins.json')).get('ground_2S',3))" 2>/dev/null || echo 3)
echo "[aws] ground 2S=$BEST — full pipeline both substrates" >> emb_run.log
for st in geom polish profile refine; do
  ( run $st 16 ch4 ) &
  P1=$!
  ( run $st 16 c2h6 ) &
  P2=$!
  wait $P1 $P2
  sync_up
  echo "[aws] stage $st done (both substrates)" >> emb_run.log
done
OCM_SPIN=$BEST python3 -u routes/ocm_mnw_selectivity.py merge2 >> emb_run.log 2>&1
sync_up
echo "[aws] ALL DONE" >> emb_run.log
sync_up
poweroff
"""


def my_instances():
    r = ec2.describe_instances(Filters=[
        {"Name": "tag:Project", "Values": [TAG_PROJECT]},
        {"Name": "tag:Name", "Values": [TAG_NAME]},
        {"Name": "instance-state-name",
         "Values": ["pending", "running", "stopping", "stopped"]}])
    return [i["InstanceId"] for res in r["Reservations"] for i in res["Instances"]]


def reap(iids, why):
    iids = sorted(set(i for i in iids if i))
    if not iids:
        return
    print(f"[reap] terminating {iids} ({why})", flush=True)
    try:
        ec2.terminate_instances(InstanceIds=iids)
    except ClientError as e:
        print(f"[reap][ERR] {e} — terminate MANUALLY: "
              f"aws ec2 terminate-instances --instance-ids {' '.join(iids)}",
              flush=True)


def s3_text(key, max_bytes=4000):
    try:
        body = s3.get_object(Bucket=BUCKET, Key=key)["Body"].read()
        return body[-max_bytes:].decode("utf-8", "replace")
    except ClientError:
        return None


def main():
    print("caller:", boto3.client("sts", region_name=REGION)
          .get_caller_identity()["Arn"], flush=True)
    print(f"branch={BRANCH} type={ITYPE} bucket=s3://{BUCKET}/{PREFIX} "
          f"reap<={MAX_MIN}min", flush=True)
    ami = ssm.get_parameter(
        Name="/aws/service/ami-amazon-linux-latest/al2023-ami-kernel-default-x86_64"
    )["Parameter"]["Value"]
    common = dict(
        ImageId=ami, InstanceType=ITYPE, MinCount=1, MaxCount=1,
        IamInstanceProfile={"Name": PROFILE},
        InstanceInitiatedShutdownBehavior="terminate",
        MetadataOptions={"HttpTokens": "required",
                         "HttpPutResponseHopLimit": 1, "HttpEndpoint": "enabled"},
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
        sys.exit(f"ocm-stage16 instance(s) already exist: {existing} — "
                 "terminate or wait first")
    iid = ec2.run_instances(**common)["Instances"][0]["InstanceId"]
    print(f"launched {iid}", flush=True)
    t0 = time.time()
    last_log = ""
    try:
        while time.time() - t0 < MAX_MIN * 60:
            try:
                s3.head_object(Bucket=BUCKET, Key=FINAL_KEY)
                s3.download_file(BUCKET, FINAL_KEY,
                                 "routes/ocm_mnw_emb_final.json")
                for extra in ("ocm_mnw_emb_spins.json",
                              "ocm_mnw_emb_ch4_geom.json",
                              "ocm_mnw_emb_c2h6_geom.json",
                              "ocm_mnw_emb_ch4_profile.json",
                              "ocm_mnw_emb_c2h6_profile.json", "emb_run.log"):
                    try:
                        s3.download_file(BUCKET, f"{PREFIX}/{extra}",
                                         f"routes/{extra}")
                    except ClientError:
                        pass
                v = json.load(open("routes/ocm_mnw_emb_final.json"))
                print("[ok] FINAL:", json.dumps(
                    {k: v.get(k) for k in ("spin_2S", "cas", "ddE_kcal",
                                           "quantum_shift_kcal")},
                    ensure_ascii=False), flush=True)
                return
            except ClientError:
                pass
            tail = s3_text(f"{PREFIX}/emb_run.log", 1500)
            if tail and tail != last_log:
                new = tail[len(last_log):] if tail.startswith(last_log) else tail
                for ln in new.strip().splitlines()[-4:]:
                    print(f"[log] {ln}", flush=True)
                last_log = tail
            try:
                st = ec2.describe_instances(InstanceIds=[iid])["Reservations"][0][
                    "Instances"][0]["State"]["Name"]
            except Exception:
                st = "unknown"
            if st in ("shutting-down", "terminated"):
                print(f"[end] instance {st}; no FINAL in S3 — job died; "
                      f"check s3://{BUCKET}/{PREFIX}/emb_run.log", flush=True)
                return
            print(f"[wait] {int((time.time()-t0)/60)}min, instance={st}",
                  flush=True)
            time.sleep(120)
        print("[deadline] reaper deadline hit", flush=True)
    finally:
        reap([iid] + my_instances(), "done/deadline/error")


if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        print("ERROR:", type(e).__name__, e, file=sys.stderr)
        try:
            reap(my_instances(), "orchestrator exception")
        except Exception:
            pass
        raise
