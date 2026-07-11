#!/usr/bin/env bash
cd /home/user/ocm-research-diary
for att in $(seq 1 12); do
  echo "[loop] attempt $att: r7i.4xlarge (16 vCPU)"
  OCM16_ITYPE=r7i.4xlarge OCM16_MAX_MIN=750 /opt/qc-venv/bin/python -u routes/ocm_aws_stage16.py; rc=$?
  [ $rc -eq 0 ] && { echo "[loop] SUCCESS"; exit 0; }
  echo "[loop] attempt $att: r7i.2xlarge (8 vCPU)"
  OCM16_ITYPE=r7i.2xlarge OCM16_MAX_MIN=900 /opt/qc-venv/bin/python -u routes/ocm_aws_stage16.py; rc=$?
  [ $rc -eq 0 ] && { echo "[loop] SUCCESS"; exit 0; }
  echo "[loop] attempt $att failed (rc=$rc); sleep 180"
  sleep 180
done
echo "[loop] EXHAUSTED"
