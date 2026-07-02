#!/usr/bin/env python3
"""
Leak-proof EC2 orchestrator — runs FROM the Claude Code sandbox (with the scoped
IAM user's creds in env) to execute the real-CHA alpha-O CASSCF+NEVPT2 job on a
big instance, retrieve the result from S3, and GUARANTEE the instance dies.

Triple cleanup (addresses the adversarial review of aws_run.sh, whose fallback
`shutdown` only STOPPED the instance and didn't catch SIGHUP):
  1. Launch with InstanceInitiatedShutdownBehavior=terminate  -> any poweroff terminates.
  2. On-instance `shutdown -P` watchdog scheduled FIRST in user-data -> terminates
     no matter what the payload does (even if pip/job hangs).
  3. Sandbox-side reaper: this script terminates the instance BY TAG in a finally:
     block after results appear OR a hard deadline -> the reaper is external, so a
     hung instance can't leak. (If the sandbox itself dies, 1+2 still terminate.)
Security: IMDSv2 required (HttpTokens=required, hop-limit 1) so a poisoned pip dep
can't exfil the instance role; the instance uses its INSTANCE PROFILE for S3 (no
AWS keys ever on the box / in user-data). No inbound (default SG), no SSH/key-pair.

Env (from the setup): AWS_* creds, AWS_DEFAULT_REGION, ALPHA_O_S3=s3://bucket/prefix,
ALPHA_O_INSTANCE_PROFILE, ALPHA_O_TAG. Optional: ALPHA_O_ITYPE, ALPHA_O_MAX_MIN,
ALPHA_O_REPO, ALPHA_O_BRANCH, ALPHA_O_DRYRUN=1.
"""
import os, sys, time, json, boto3
from botocore.exceptions import ClientError

REGION = os.environ["AWS_DEFAULT_REGION"]
S3_URI = os.environ["ALPHA_O_S3"].replace("s3://", "")
BUCKET, PREFIX = S3_URI.split("/", 1); PREFIX = PREFIX.rstrip("/")
PROFILE = os.environ["ALPHA_O_INSTANCE_PROFILE"]
TAG = os.environ.get("ALPHA_O_TAG", "sandbox-agent")
ITYPE = os.environ.get("ALPHA_O_ITYPE", "c7i.8xlarge")
REPO = os.environ.get("ALPHA_O_REPO", "https://github.com/skharkov1246/ocm-research-diary")
BRANCH = os.environ.get("ALPHA_O_BRANCH", "claude/affectionate-franklin-jkc8gx")
MAX_MIN = int(os.environ.get("ALPHA_O_MAX_MIN", "300"))     # sandbox reaper deadline
JOB_MIN = MAX_MIN - 20                                        # on-instance job/watchdog budget
KEY = f"{PREFIX}/fe_zeolite_cha_results.json"
LOGKEY = f"{PREFIX}/fe_zeolite_cha.log"

ec2 = boto3.client("ec2", region_name=REGION)
ssm = boto3.client("ssm", region_name=REGION)
s3 = boto3.client("s3", region_name=REGION)

USERDATA = f"""#!/bin/bash
# 2) independent on-instance watchdog FIRST: terminate no matter what the job does
shutdown -P +{JOB_MIN}
set -x
dnf install -y python3-pip git gcc gcc-c++ awscli >/tmp/setup.log 2>&1 || python3 -m pip install -q awscli
python3 -m pip install -q numpy scipy pyscf >>/tmp/setup.log 2>&1
cd /root && git clone --depth 1 -b {BRANCH} {REPO} repo >>/tmp/setup.log 2>&1 && cd repo
export OMP_NUM_THREADS=$(nproc) PYSCF_MAX_MEMORY=$(( $(getconf _PHYS_PAGES) * $(getconf PAGE_SIZE) / 1024 / 1024 * 6 / 10 ))
CHA_NEVPT2=1 timeout {JOB_MIN-10}m python3 -u calc/fe_zeolite_cha_sp.py > calc/fe_zeolite_cha.log 2>&1
aws s3 cp calc/fe_zeolite_cha.log         s3://{BUCKET}/{LOGKEY} || true
aws s3 cp calc/fe_zeolite_cha_results.json s3://{BUCKET}/{KEY}
poweroff   # 1) InstanceInitiatedShutdownBehavior=terminate -> this terminates on success
"""


def find_running():
    r = ec2.describe_instances(Filters=[{"Name": "tag:Project", "Values": [TAG]},
                                        {"Name": "instance-state-name", "Values": ["pending", "running", "stopping", "stopped"]}])
    return [i["InstanceId"] for res in r["Reservations"] for i in res["Instances"]]


def reap(iids, why):
    if not iids:
        return
    print(f"[reap] terminating {iids} ({why})", flush=True)
    try:
        ec2.terminate_instances(InstanceIds=iids)
    except ClientError as e:
        print(f"[reap][ERR] {e} — terminate MANUALLY: aws ec2 terminate-instances --instance-ids {' '.join(iids)}", flush=True)


def main():
    print(f"caller: {boto3.client('sts', region_name=REGION).get_caller_identity()['Arn']}", flush=True)
    ami = ssm.get_parameter(Name="/aws/service/ami-amazon-linux-latest/al2023-ami-kernel-default-x86_64")["Parameter"]["Value"]
    print(f"AMI={ami} type={ITYPE} region={REGION} reap<={MAX_MIN}min", flush=True)
    common = dict(ImageId=ami, InstanceType=ITYPE, MinCount=1, MaxCount=1,
                  IamInstanceProfile={"Name": PROFILE},
                  InstanceInitiatedShutdownBehavior="terminate",
                  MetadataOptions={"HttpTokens": "required", "HttpPutResponseHopLimit": 1, "HttpEndpoint": "enabled"},
                  BlockDeviceMappings=[{"DeviceName": "/dev/xvda", "Ebs": {"VolumeSize": 30, "VolumeType": "gp3", "DeleteOnTermination": True}}],
                  TagSpecifications=[{"ResourceType": "instance", "Tags": [{"Key": "Project", "Value": TAG}, {"Key": "Name", "Value": "alpha-o"}]},
                                     {"ResourceType": "volume", "Tags": [{"Key": "Project", "Value": TAG}]}],
                  UserData=USERDATA)
    if os.environ.get("ALPHA_O_DRYRUN") == "1":
        try:
            ec2.run_instances(DryRun=True, **common)
        except ClientError as e:
            ok = "DryRunOperation" in str(e)
            print(("[dryrun] OK — permissions sufficient to launch" if ok else f"[dryrun] FAILED: {e}"), flush=True)
        return
    iid = ec2.run_instances(**common)["Instances"][0]["InstanceId"]
    print(f"launched {iid}", flush=True)
    t0 = time.time()
    try:
        while time.time() - t0 < MAX_MIN * 60:
            try:
                s3.head_object(Bucket=BUCKET, Key=KEY)
                s3.download_file(BUCKET, KEY, "calc/fe_zeolite_cha_results.json")
                print("[ok] results retrieved from S3", flush=True)
                v = json.load(open("calc/fe_zeolite_cha_results.json"))
                print("VERDICT:", json.dumps({k: v.get(k) for k in ("verdict_e_cas", "verdict_e_nevpt2")}, ensure_ascii=False), flush=True)
                return
            except ClientError:
                pass
            print(f"[wait] {int((time.time()-t0)/60)}min elapsed…", flush=True)
            time.sleep(60)
        print("[deadline] no results by reaper deadline", flush=True)
    finally:
        reap(find_running() or [iid], "job done / deadline / error")


if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        print("ERROR:", type(e).__name__, e, file=sys.stderr)
        reap(find_running(), "orchestrator exception")
        raise
