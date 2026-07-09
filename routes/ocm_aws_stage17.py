#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
routes/ocm_aws_stage17.py — leak-proof EC2-оркестратор Этапа 17 OCM:
АМПЛИФИКАЦИЯ анти-BEP допантом. Этап 16 показал, что голый Mn=O в реальном
Na/W-окружении селективностно-нейтрален (ΔΔE‡≈0), хотя знак квантового сдвига
выжил. Вопрос Этапа 17: даёт ли ДРУГОЙ активный оксо-металл из окна вакансионного
волкана (Cr(IV) d², Fe(IV) d⁴) в ТОМ ЖЕ embedded Na/W-каркасе положительный
(селективный) ΔΔE‡ там, где Mn нейтрален?

Оба металла считаются ПАРАЛЛЕЛЬНО на одном большом инстансе (eu-central-1, где
квота 128 vCPU): 2 металла × 2 субстрата = 4 задачи × (nproc/4) потоков. Каждый
металл: spins → читаем ground_2S → geom/polish/profile/refine (оба субстрата
параллельно) → merge2. Префиксы ocm_mnw_emb_cr_* / _fe_* не пересекаются.

Тройная защита от утечек — как в ocm_aws_stage16.py:
  1. InstanceInitiatedShutdownBehavior=terminate;
  2. бортовой watchdog `shutdown -P +N` первой строкой user-data;
  3. внешний жнец в песочнице гасит инстанс по тегу к дедлайну (finally).

Запуск:  AWS_DEFAULT_REGION=eu-central-1 python3 routes/ocm_aws_stage17.py [--dryrun]
Env: OCM17_ITYPE (r7i.12xlarge), OCM17_MAX_MIN (1100), OCM17_BRANCH (тек. ветка),
     OCM17_METALS ("Cr Fe").
"""
import json
import os
import subprocess
import sys
import time

import boto3
from botocore.exceptions import ClientError

REGION = os.environ.get("AWS_DEFAULT_REGION", "eu-central-1")
ACCOUNT = "097743207937"
BUCKET = f"alpha-o-results-{ACCOUNT}"          # bucket в us-east-1 — cross-region OK
BUCKET_REGION = "us-east-1"
PREFIX = os.environ.get("OCM17_PREFIX", "alpha-o/ocm17")   # отд. префикс для
                                                           # solo-прогона (без клоббера)
PROFILE = "alpha-o-instance-profile"           # IAM instance profile — глобальный
TAG_PROJECT = "ocm-agent"
TAG_NAME = os.environ.get("OCM17_TAG", "ocm-stage17")   # отд. тег → параллельный
                                                        # solo-прогон металла
ITYPE = os.environ.get("OCM17_ITYPE", "r7i.12xlarge")   # 48 vCPU, 384 GB
SUBNET = {"eu-central-1": "subnet-057e6f5ca6d3c6d8d",
          "us-east-1": "subnet-0b1a363f27ecbbf12"}[REGION]
METALS = os.environ.get("OCM17_METALS", "Cr Fe")
GRID = os.environ.get("OCM17_GRID", "")   # "ext" → расширенная сетка (Этап 18)
REPO = "https://github.com/skharkov1246/ocm-research-diary"
BRANCH = os.environ.get("OCM17_BRANCH") or subprocess.run(
    ["git", "rev-parse", "--abbrev-ref", "HEAD"], capture_output=True,
    text=True).stdout.strip()
MAX_MIN = int(os.environ.get("OCM17_MAX_MIN", "1100"))
JOB_MIN = MAX_MIN - 20
JOB_TIMEOUT = JOB_MIN - 20
DISK_GB = 120
FINAL_KEYS = {m: f"{PREFIX}/ocm_mnw_emb_{m.lower()}_final.json"
              for m in METALS.split()}
# глоб синка выводим ИЗ METALS (не хардкод cr/fe): иначе прогон другого металла
# заливал бы закоммиченные чужие JSON и пропускал свои реальные результаты
SYNC_INCLUDES = " ".join(f"--include 'ocm_mnw_emb_{m.lower()}_*'"
                         for m in METALS.split())

for var in ("AWS_ACCESS_KEY_ID", "AWS_SECRET_ACCESS_KEY", "AWS_SESSION_TOKEN"):
    os.environ.pop(var, None)          # битые env-креды → ~/.aws/credentials

ec2 = boto3.client("ec2", region_name=REGION)
ssm = boto3.client("ssm", region_name=REGION)
s3 = boto3.client("s3", region_name=BUCKET_REGION)

# per-metal функция: spins → ground_2S → полный пайплайн обоих субстратов → merge2.
# Оба металла запускаются параллельно (см. хвост user-data); каждый — в subshell,
# чтобы OCM_METAL/OCM_SPIN/PREF были изолированы. Внутри металла оба субстрата
# тоже параллельны. Итого 4 задачи × THREADS потоков.
USERDATA = f"""#!/bin/bash
shutdown -P +{JOB_MIN}
exec > /tmp/boot.log 2>&1
set -x
export PATH=$PATH:/usr/local/bin:/root/.local/bin
S3=s3://{BUCKET}/{PREFIX}
export AWS_DEFAULT_REGION={BUCKET_REGION}
for i in 1 2 3; do timeout 10m dnf install -y python3-pip git gcc gcc-c++ && break; sleep 30; done
for i in 1 2 3; do timeout 10m python3 -m pip install -q awscli && break; sleep 30; done
( while true; do aws s3 cp /tmp/boot.log $S3/boot.log >/dev/null 2>&1; sleep 60; done ) &
for i in 1 2 3; do timeout 25m python3 -m pip install -q numpy scipy pyscf geometric && break; sleep 30; done
for i in 1 2 3; do (cd /root && timeout 10m git clone --depth 1 -b {BRANCH} {REPO} repo) && break; sleep 30; done
cd /root/repo || {{ aws s3 cp /tmp/boot.log $S3/setup_failed.log; poweroff; }}
python3 -c "import pyscf, geometric" || {{ aws s3 cp /tmp/boot.log $S3/setup_failed.log; poweroff; }}
# стейл-финалы из git УДАЛЯЕМ до старта: иначе первый же sync_up заливает их в
# свежий S3-префикс, оркестратор принимает за готовый результат и гасит инстанс
# посреди счёта (случилось на relaunch Этапа 18). merge2 пересоздаст в конце.
for M in {METALS}; do rm -f routes/ocm_mnw_emb_$(echo $M | tr A-Z a-z)_final.json; done
mkdir -p /root/scratch
export PYSCF_TMPDIR=/root/scratch
export OCM_MODEL=embedded
export OCM_GRID={GRID}
NPROC=$(nproc)
THREADS=$(( NPROC / 4 )); [ $THREADS -lt 1 ] && THREADS=1
sync_up() {{ aws s3 cp routes/ $S3/ --recursive --exclude '*' {SYNC_INCLUDES} >/dev/null 2>&1 || true; aws s3 cp /root/repo/stage17.log $S3/stage17.log >/dev/null 2>&1 || true; }}
( while true; do sleep 180; sync_up; done ) &

run_metal() {{
  local M=$1
  (
    export OCM_METAL=$M
    echo "[aws][$M] spins ladder" >> /root/repo/stage17.log
    OMP_NUM_THREADS=$THREADS timeout 120m python3 -u routes/ocm_mnw_selectivity.py spins >> /root/repo/stage17.log 2>&1
    local SP="routes/ocm_mnw_emb_$(echo $M | tr A-Z a-z)_spins.json"
    local BEST=$(python3 -c "import json;print(json.load(open('$SP')).get('ground_2S',2))" 2>/dev/null || echo 2)
    export OCM_SPIN=$BEST
    echo "[aws][$M] ground 2S=$BEST — full pipeline both substrates" >> /root/repo/stage17.log
    local st
    for st in geom polish profile refine; do
      OMP_NUM_THREADS=$THREADS timeout {JOB_TIMEOUT}m python3 -u routes/ocm_mnw_selectivity.py $st ch4 >> /root/repo/stage17.log 2>&1 &
      local A=$!
      OMP_NUM_THREADS=$THREADS timeout {JOB_TIMEOUT}m python3 -u routes/ocm_mnw_selectivity.py $st c2h6 >> /root/repo/stage17.log 2>&1 &
      local B=$!
      wait $A $B
      sync_up
      echo "[aws][$M] stage $st done" >> /root/repo/stage17.log
    done
    python3 -u routes/ocm_mnw_selectivity.py merge2 >> /root/repo/stage17.log 2>&1
    sync_up
    echo "[aws][$M] FINAL done" >> /root/repo/stage17.log
  )
}}

PIDS=""
for M in {METALS}; do run_metal $M & PIDS="$PIDS $!"; done
wait $PIDS
sync_up
echo "[aws] ALL METALS DONE" >> /root/repo/stage17.log
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
              f"aws ec2 --region {REGION} terminate-instances "
              f"--instance-ids {' '.join(iids)}", flush=True)


def s3_text(key, max_bytes=1500):
    try:
        body = s3.get_object(Bucket=BUCKET, Key=key)["Body"].read()
        return body[-max_bytes:].decode("utf-8", "replace")
    except ClientError:
        return None


def pull_metal(m):
    """Скачать FINAL и промежуточные JSON конкретного металла из S3 в routes/."""
    ml = m.lower()
    got = False
    try:
        s3.download_file(BUCKET, FINAL_KEYS[m],
                         f"routes/ocm_mnw_emb_{ml}_final.json")
        got = True
    except ClientError:
        pass
    for sub in ("ch4", "c2h6"):
        for kind in ("geom", "profile"):
            k = f"{PREFIX}/ocm_mnw_emb_{ml}_{sub}_{kind}.json"
            try:
                s3.download_file(BUCKET, k,
                                 f"routes/ocm_mnw_emb_{ml}_{sub}_{kind}.json")
            except ClientError:
                pass
    try:
        s3.download_file(BUCKET, f"{PREFIX}/ocm_mnw_emb_{ml}_spins.json",
                         f"routes/ocm_mnw_emb_{ml}_spins.json")
    except ClientError:
        pass
    return got


def main():
    print("caller:", boto3.client("sts", region_name=REGION)
          .get_caller_identity()["Arn"], flush=True)
    print(f"branch={BRANCH} region={REGION} type={ITYPE} metals={METALS} "
          f"bucket=s3://{BUCKET}/{PREFIX} reap<={MAX_MIN}min", flush=True)
    try:
        ami = ssm.get_parameter(
            Name="/aws/service/ami-amazon-linux-latest/"
                 "al2023-ami-kernel-default-x86_64")["Parameter"]["Value"]
    except ClientError:                     # ssm:GetParameter запрещён в регионе
        imgs = ec2.describe_images(
            Owners=["amazon"],
            Filters=[{"Name": "name",
                      "Values": ["al2023-ami-2023*-kernel-*-x86_64"]},
                     {"Name": "state", "Values": ["available"]}])["Images"]
        ami = max(imgs, key=lambda i: i["CreationDate"])["ImageId"]
    print("ami:", ami, flush=True)
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
        sys.exit(f"ocm-stage17 instance(s) already exist: {existing} — "
                 "terminate or wait first")
    iid = ec2.run_instances(**common)["Instances"][0]["InstanceId"]
    print(f"launched {iid}", flush=True)
    t0 = time.time()
    done = set()
    last_log = ""
    try:
        while time.time() - t0 < MAX_MIN * 60:
            for m in METALS.split():
                if m in done:
                    continue
                try:
                    s3.head_object(Bucket=BUCKET, Key=FINAL_KEYS[m])
                except ClientError:
                    continue
                if pull_metal(m):
                    v = json.load(open(
                        f"routes/ocm_mnw_emb_{m.lower()}_final.json"))
                    print(f"[ok][{m}] FINAL:", json.dumps(
                        {k: v.get(k) for k in ("spin_2S", "cas", "ddE_kcal",
                                               "quantum_shift_kcal")},
                        ensure_ascii=False), flush=True)
                    done.add(m)
            if done == set(METALS.split()):
                print("[ok] all metals FINAL — done", flush=True)
                return
            tail = s3_text(f"{PREFIX}/stage17.log")
            if tail and tail != last_log:
                new = tail[len(last_log):] if tail.startswith(last_log) else tail
                print(new, end="", flush=True)
                last_log = tail
            time.sleep(90)
        print("[timeout] MAX_MIN reached", flush=True)
    finally:
        reap(my_instances(), "deadline/exit")


if __name__ == "__main__":
    main()
