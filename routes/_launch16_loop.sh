#!/usr/bin/env bash
cd /home/user/ocm-research-diary
while true; do
  echo "[loop] try r7i.4xlarge (16 vCPU)"
  OCM16_ITYPE=r7i.4xlarge OCM16_MAX_MIN=750 /opt/qc-venv/bin/python -u routes/ocm_aws_stage16.py && break
  echo "[loop] try r7i.2xlarge (8 vCPU)"
  OCM16_ITYPE=r7i.2xlarge OCM16_MAX_MIN=900 /opt/qc-venv/bin/python -u routes/ocm_aws_stage16.py && break
  echo "[loop] both denied; sleep 180"
  sleep 180
done
echo "[loop] launcher finished"
