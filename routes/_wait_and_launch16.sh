#!/usr/bin/env bash
unset AWS_ACCESS_KEY_ID AWS_SECRET_ACCESS_KEY AWS_SESSION_TOKEN
export AWS_DEFAULT_REGION=us-east-1
AWS=/opt/qc-venv/bin/aws
while true; do
  ST=$($AWS ec2 describe-instances --instance-ids i-08c72686014362bb3 \
       --query 'Reservations[0].Instances[0].State.Name' --output text 2>/dev/null || echo gone)
  echo "$(date -u +%H:%M) alpha-o: $ST"
  case "$ST" in
    terminated|gone|shutting-down) break;;
  esac
  sleep 120
done
echo "alpha-o released quota — launching ocm-stage16 (r7i.4xlarge)"
sleep 60
cd /home/user/ocm-research-diary
OCM16_ITYPE=r7i.4xlarge /opt/qc-venv/bin/python -u routes/ocm_aws_stage16.py
