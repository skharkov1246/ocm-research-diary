#!/usr/bin/env python3
"""Recycle-proof fleet status for the Frankfurt (eu-central-1) diary compute fleet.

Design principle: the sandbox container is ephemeral — it can be recycled at any
time, wiping local scratch, pip packages and even reverting the git checkout.
So this script keeps ZERO durable state locally. It reconstructs everything from
the two things that survive a recycle:
  * S3   — the results/logs each instance uploads (source of truth for outputs)
  * AWS  — describe_instances in eu-central-1 (source of truth for what's alive)

Convention (set by the operator): us-east-1 is the OCM track's region; every
running instance in eu-central-1 belongs to THIS 4-product fleet. So we identify
our jobs by region, not by tags/markers (the minimal IAM has no CreateTags).

Usage:  python3 -m pip install -q boto3 ; python3 calc/fleet_status.py
Env:    AWS_* creds, AWS_DEFAULT_REGION, ALPHA_O_S3=s3://bucket/prefix
"""
import os, json, datetime, boto3
from botocore.exceptions import ClientError

S3URI = os.environ["ALPHA_O_S3"].replace("s3://", "").rstrip("/")
BUCKET, PREFIX = S3URI.split("/", 1)
s3 = boto3.client("s3", region_name=os.environ["AWS_DEFAULT_REGION"])
ec2 = boto3.client("ec2", region_name="eu-central-1")

# product -> (S3 result key basename, repo destination)
RESULTS = {
    "ch4":   ("fe_zeolite_cha_results.json",      "calc/fe_zeolite_cha_results.json"),
    "fe4s4": ("fe4s4_casscf_results.json",        "calc/fe4s4_casscf_results.json"),
    "poly":  ("polyolefins_chiral_results.json",  "routes/polyolefins_chiral_results.json"),
    # co2 was completed earlier by the fleet and is already in the diary.
}

def result_status(name, key_base):
    key = f"{PREFIX}/{key_base}"
    try:
        b = s3.get_object(Bucket=BUCKET, Key=key)["Body"].read()
    except ClientError:
        return "absent", None, 0
    try:
        d = json.loads(b)
        st = d.get("status", "?")
    except Exception:
        d, st = None, "unparseable"
    final = st in ("ok", "done", "complete", "finished")
    return ("FINAL" if final else f"partial({st})"), d, len(b)

def main():
    now = datetime.datetime.now(datetime.timezone.utc)
    print("=== Frankfurt fleet (eu-central-1) instances — all are ours:")
    alive = []
    for r in ec2.describe_instances(Filters=[{"Name": "instance-state-name",
                                              "Values": ["pending", "running", "shutting-down"]}])["Reservations"]:
        for i in r["Instances"]:
            age = (now - i["LaunchTime"]).total_seconds() / 3600
            alive.append(i["InstanceId"])
            print(f"  {i['InstanceId']} {i['InstanceType']} {i['State']['Name']} age={age:.1f}h")
    if not alive:
        print("  (none running)")
    print("=== results in S3:")
    done = []
    for name, (kb, dest) in RESULTS.items():
        tag, d, sz = result_status(name, kb)
        print(f"  {name:6} {tag:16} {sz}B  <- s3://{BUCKET}/{PREFIX}/{kb}")
        if tag == "FINAL":
            done.append((name, dest, kb))
    if done:
        print("=== FINAL results ready to write (download with --pull):")
        for name, dest, kb in done:
            print(f"  {name} -> {dest}")
    return alive, done

if __name__ == "__main__":
    alive, done = main()
    import sys
    if "--pull" in sys.argv:
        for name, dest, kb in done:
            s3.download_file(BUCKET, f"{PREFIX}/{kb}", dest)
            print(f"pulled {name} -> {dest}")
