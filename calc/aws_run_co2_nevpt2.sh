#!/usr/bin/env bash
# ============================================================================
# calc/aws_run_co2_nevpt2.sh — EC2 user-data / on-instance launcher for the
# co2-to-fuels NEVPT2 job (issue #13): SC-NEVPT2 on top of the fixed-active-
# space CASSCF rehabilitation (routes/co2_to_fuels_nevpt2_aws.py, protocol
# imported 1:1 from routes/co2_to_fuels_casscf.py), reactant/TS x Cu2/CuAl.
#
# Leak-proof pattern of calc/aws_run.sh + calc/aws_orchestrate.py:
#   1) dead-man watchdog `shutdown -P +N` scheduled FIRST — fires even if
#      dnf/pip/git/python hang. Launch the instance with
#      InstanceInitiatedShutdownBehavior=terminate so ANY poweroff TERMINATES
#      (not stops) — see calc/aws_orchestrate.py run_instances(common).
#   2) trap terminate EXIT — explicit IMDSv2 self-terminate as second layer
#      (uploads whatever results exist first).
#   3) hard `timeout` on the payload (TIMEOUT, default 21600 s = 6 h) plus a
#      per-stage timeout inside python (NEVPT2_STAGE_TIMEOUT_S, default 5400 s).
#
# The CASSCF checkpoints (*.chk.npz) are NOT committed — the instance honestly
# recomputes SCF+CASSCF from scratch on the held branch, cross-checks against
# the committed routes/co2_to_fuels_casscf_results.json, then runs NEVPT2.
# Results JSON is atomic-saved after every stage and also pushed to S3 every
# 10 min (progress key), so a killed instance still leaves partial results.
#
# Usage as user-data (root, Amazon Linux 2023 x86_64), env baked in by the
# orchestrator or exported before sourcing:
#   S3_OUT=s3://bucket/prefix    # instance-profile creds; no AWS keys on box
#   TIMEOUT=21600 BRANCH=main    # optional overrides
# Suggested instance: r7i.2xlarge (8 vCPU / 64 GiB) — the job is small
# (105–118 basis functions, CAS <= (14e,11o)); 64 GiB is headroom, not need.
# ============================================================================
set -uo pipefail

TIMEOUT="${TIMEOUT:-21600}"                                   # 6 h payload cap
REPO="${REPO:-https://github.com/skharkov1246/ocm-research-diary}"
BRANCH="${BRANCH:-main}"
S3_OUT="${S3_OUT:-}"
export NEVPT2_STAGE_TIMEOUT_S="${NEVPT2_STAGE_TIMEOUT_S:-5400}"
WATCHDOG_MIN="${WATCHDOG_MIN:-$(( TIMEOUT / 60 + 45 ))}"      # payload+setup+upload

# ---- layer 1: dead-man watchdog FIRST (before anything can hang) -----------
shutdown -P "+${WATCHDOG_MIN}" 2>/dev/null \
  || sudo shutdown -P "+${WATCHDOG_MIN}" 2>/dev/null || true

RESULTS=routes/co2_to_fuels_nevpt2_results.json
LOG=routes/co2_to_fuels_nevpt2_aws.log

upload() {   # best-effort: results, log, restart checkpoints, diagnostics
  [ -n "${S3_OUT}" ] || { echo "[warn] S3_OUT not set — nothing uploaded"; return 0; }
  for f in "${RESULTS}" "${LOG}" \
           routes/co2_to_fuels_nevpt2_Cu2.chk.npz \
           routes/co2_to_fuels_nevpt2_CuAl.chk.npz \
           /tmp/setup.log; do
    if [ -f "${f}" ]; then
      aws s3 cp "${f}" "${S3_OUT%/}/$(basename "${f}")" || true
    fi
  done
  # full user-data trace — the only forensic artifact if setup dies early
  if [ -f /var/log/cloud-init-output.log ]; then
    aws s3 cp /var/log/cloud-init-output.log \
      "${S3_OUT%/}/co2_nevpt2_cloud-init-output.log" || true
  fi
}

terminate() {   # layer 2: explicit self-terminate via IMDSv2; poweroff fallback
  echo "[dead-man] uploading leftovers and terminating instance…"
  upload || true
  TOK="$(curl -sS -X PUT --max-time 5 \
        -H 'X-aws-ec2-metadata-token-ttl-seconds: 60' \
        http://169.254.169.254/latest/api/token || true)"
  IID="$(curl -sS --max-time 5 -H "X-aws-ec2-metadata-token: ${TOK}" \
        http://169.254.169.254/latest/meta-data/instance-id || true)"
  REG="$(curl -sS --max-time 5 -H "X-aws-ec2-metadata-token: ${TOK}" \
        http://169.254.169.254/latest/meta-data/placement/region || true)"
  if [ -n "${IID}" ] && command -v aws >/dev/null 2>&1; then
    aws ec2 terminate-instances --instance-ids "${IID}" --region "${REG}" \
      || shutdown -P now || sudo shutdown -P now || true
  else
    shutdown -P now || sudo shutdown -P now || true
  fi
}
trap terminate EXIT

set -x
export PATH="${PATH}:/usr/local/bin:/root/.local/bin"

# ---- setup: robust installs, every network step under its own timeout ------
timeout 15m dnf install -y python3-pip git gcc gcc-c++ >/tmp/setup.log 2>&1 || true
timeout 25m python3 -m pip install -q --upgrade pip           >>/tmp/setup.log 2>&1
timeout 25m python3 -m pip install -q awscli numpy scipy pyscf >>/tmp/setup.log 2>&1

cd /root || exit 1
timeout 10m git clone --depth 1 -b "${BRANCH}" "${REPO}" co2f-nevpt2 \
  >>/tmp/setup.log 2>&1 || { echo "SETUP-FAILED (clone)"; exit 1; }   # trap fires
cd co2f-nevpt2 || exit 1

mkdir -p /root/scratch
export PYSCF_TMPDIR=/root/scratch
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-$(nproc)}"
# 60% of RAM in MB (pyscf's in-core ceiling); leave room for the 4-RDM peaks
export PYSCF_MAX_MEMORY="${PYSCF_MAX_MEMORY:-$(( $(getconf _PHYS_PAGES) * $(getconf PAGE_SIZE) / 1024 / 1024 * 6 / 10 ))}"

# ---- 10-min progress uploader (partial results visible before the end) -----
if [ -n "${S3_OUT}" ]; then
  ( while true; do
      sleep 600
      [ -f "${RESULTS}" ] && aws s3 cp "${RESULTS}" \
        "${S3_OUT%/}/co2_to_fuels_nevpt2_results.progress.json" \
        >/dev/null 2>&1 || true
    done ) &
  UPLOADER_PID=$!
else
  UPLOADER_PID=""
fi

echo "[run] NEVPT2 job: timeout ${TIMEOUT}s, ${OMP_NUM_THREADS} threads," \
     "PYSCF_MAX_MEMORY=${PYSCF_MAX_MEMORY} MB," \
     "stage timeout ${NEVPT2_STAGE_TIMEOUT_S}s"
timeout "${TIMEOUT}" python3 -u routes/co2_to_fuels_nevpt2_aws.py \
  >"${LOG}" 2>&1
RC=$?
echo "[done] rc=${RC} (124 = wallclock timeout; partial results are still valid)"
tail -n 40 "${LOG}" || true

[ -n "${UPLOADER_PID}" ] && kill "${UPLOADER_PID}" 2>/dev/null || true
upload
# trap fires here -> instance self-terminates (terminate, not stop)
