#!/usr/bin/env bash
# ============================================================================
# AWS launcher for the Fe-zeolite alpha-O cluster validation, with a dead-man
# switch (same pattern the diary uses for FeMoco / NiO / Ti runs): the EC2
# instance SELF-TERMINATES on exit — success, failure, or timeout — so nothing
# leaks. Run this ON the instance, from the repo root.
#
#   ./calc/aws_run.sh                 # idealized model (smoke of the full path)
#   ./calc/aws_run.sh site.xyz 0      # real cluster geometry + total charge
#
# Recommended instance: compute-optimized, >=32 GB RAM (e.g. c7i.4xlarge /
# r7i.2xlarge). CASSCF for the spin states is the core cost; the full-CAS VQE
# statevector is the heavy lift (that is why this is an AWS job, not local).
# ============================================================================
set -uo pipefail

XYZ="${1:-}"
CHARGE="${2:-0}"
TIMEOUT="${TIMEOUT:-21600}"   # 6 h hard cap

terminate() {
  echo "[dead-man] terminating instance to avoid leaks…"
  # prefer a clean self-terminate via the instance's own metadata
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
python3 -m pip install -q numpy scipy pyscf pennylane pennylane-lightning geometric

ARGS=(--full --charge "${CHARGE}")
[ -n "${XYZ}" ] && ARGS+=(--xyz "${XYZ}")

echo "[run] python3 calc/fe_zeolite_alpha_o.py ${ARGS[*]}  (timeout ${TIMEOUT}s)"
timeout "${TIMEOUT}" python3 calc/fe_zeolite_alpha_o.py "${ARGS[@]}" 2>&1 | tee calc/fe_zeolite_alpha_o.log
echo "[done] results in calc/fe_zeolite_alpha_o_results.json"
# trap fires here -> instance self-terminates
