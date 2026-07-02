#!/usr/bin/env bash
# ============================================================================
# AWS launcher — hydrogen track (h2-electrocatalysis), Stages 2b/3 heavy part.
# Same dead-man-switch pattern as calc/aws_run.sh (FeMoco/NiO runs): the EC2
# instance SELF-TERMINATES on exit (success, failure, or timeout).
# Run ON the instance, from the repo root:
#
#   ./routes/h2_aws_run.sh              # dinuclear ladder + quantum probe + VQE
#   TIMEOUT=28800 ./routes/h2_aws_run.sh
#
# Recommended: c7i.4xlarge (16 vCPU / 32 GB) — dinuclear UKS opts + 16-20q VQE.
# NOTE (2026-07-02): the web-session container has NO valid AWS credentials
# (sts InvalidClientTokenId) — this script is prepared for a manually started
# instance, exactly like the 2026-06-07 NiO run. Local 4-core fallback:
# run the same three commands without the wrapper.
# ============================================================================
set -uo pipefail

TIMEOUT="${TIMEOUT:-21600}"   # 6 h hard cap

terminate() {
  echo "[dead-man] terminating instance to avoid leaks…"
  IID="$(curl -s --max-time 3 http://169.254.169.254/latest/meta-data/instance-id || true)"
  REG="$(curl -s --max-time 3 http://169.254.169.254/latest/meta-data/placement/region || true)"
  if [ -n "${IID}" ] && command -v aws >/dev/null 2>&1; then
    aws ec2 terminate-instances --instance-ids "${IID}" --region "${REG}" || sudo shutdown -h now || true
  else
    sudo shutdown -h now || true
  fi
}
trap terminate EXIT

python3 -m pip install -q --upgrade pip
python3 -m pip install -q numpy scipy pyscf pyberny pennylane pennylane-lightning openfermion

run() { echo "[run] $*"; timeout "${TIMEOUT}" "$@" 2>&1 | tee -a routes/h2_aws.log; }

# 1) mononuclear ladder (fills gas references too; resumes if JSON present)
run python3 -u routes/h2_oer_ladder.py
# 2) dinuclear NiOOH-edge ladder
run python3 -u routes/h2_oer_dinuclear.py
# 3) quantum probe + converged UCCSD-VQE on the limiting step (up to 20q here)
H2_VQE=1 H2_VQE_MAXQ=20 run python3 -u routes/h2_oer_quantum.py --source di

echo "[done] results: routes/h2_oer_*_results.json (grab them BEFORE termination)"
# trap fires here -> instance self-terminates
