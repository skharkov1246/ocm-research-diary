#!/usr/bin/env bash
# ============================================================================
# AWS launcher for the FAITHFUL real-CHA alpha-O spin-state job (CASSCF+NEVPT2
# on the real IZA CHA framework geometry), with a dead-man switch: the EC2
# instance SELF-TERMINATES on exit — success, failure, or timeout — so nothing
# leaks. Run this ON the instance, from the repo root:
#
#   ./calc/aws_run.sh
#
# Everything needed is committed and pyscf-only — the cluster is already built
# from the real IZA CHA framework and xTB-relaxed (calc/alpha_o_cha_relaxed.xyz),
# so NO ase/zse/tblite on the instance, just pyscf. Curated CAS(16e,11o)=22q.
#
# WHY AWS: on a 4-vCPU box the DF-SCF alone took ~1 h (465 bf) and CASSCF+NEVPT2
# x3 did not finish in hours. Recommended instance: >=16 vCPU, >=64 GB RAM
# (e.g. c7i.8xlarge / r7i.4xlarge). Copy the results off BEFORE it self-terminates.
# ============================================================================
set -uo pipefail
TIMEOUT="${TIMEOUT:-43200}"          # 12 h hard cap

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

export OMP_NUM_THREADS="${OMP_NUM_THREADS:-$(nproc)}"
export PYSCF_MAX_MEMORY="${PYSCF_MAX_MEMORY:-$(( $(getconf _PHYS_PAGES) * $(getconf PAGE_SIZE) / 1024 / 1024 * 8 / 10 ))}"
python3 -m pip install -q --upgrade pip
python3 -m pip install -q numpy scipy pyscf

echo "[run] CHA_NEVPT2=1 python3 calc/fe_zeolite_cha_sp.py  (timeout ${TIMEOUT}s, ${OMP_NUM_THREADS} threads)"
CHA_NEVPT2=1 timeout "${TIMEOUT}" python3 -u calc/fe_zeolite_cha_sp.py 2>&1 | tee calc/fe_zeolite_cha.log
echo "[done] results -> calc/fe_zeolite_cha_results.json"
echo "[IMPORTANT] copy results off the box NOW, e.g.:  aws s3 cp calc/fe_zeolite_cha_results.json s3://YOUR_BUCKET/"
# trap fires here -> instance self-terminates
