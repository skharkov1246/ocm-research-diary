#!/usr/bin/env bash
# ============================================================================
# AWS user-data wrapper for the [4Fe-4S] BS->UNO->CASSCF + qubit-Hamiltonian
# job (calc/fe4s4_casscf_aws.py), FeMoco diary Stage 19. Modeled on
# calc/aws_run.sh + the user-data block of calc/aws_orchestrate.py.
#
# Designed to run as EC2 USER-DATA (root, no tty, fresh AL2023 box): it clones
# the repo itself, installs the stack, runs the job, uploads results to S3 and
# self-terminates. It can equally be run by hand on an instance.
#
# Dead-man layering (same philosophy as aws_orchestrate.py):
#   1. `shutdown -P +N` scheduled FIRST — fires even if everything below hangs;
#      with InstanceInitiatedShutdownBehavior=terminate a poweroff TERMINATES.
#   2. terminate() trap on EXIT — explicit terminate via IMDSv2 (the orchestrator
#      launches with HttpTokens=required, so plain IMDSv1 curl would fail),
#      falling back to poweroff.
#   3. The sandbox-side reaper in aws_orchestrate.py (kills by tag).
#
# Env (all optional except S3_OUT if you want the results to survive):
#   TIMEOUT               total seconds for the whole job, default 21600 (6 h)
#   FE4S4_REPO            git repo URL (default: this project)
#   FE4S4_BRANCH          branch to clone, default main
#   S3_OUT                s3://bucket/prefix for results upload
#   FE4S4_STAGE_TIMEOUT   per-stage cap passed to the python job (default 5400)
#   FE4S4_LADDER / FE4S4_MF / FE4S4_LIGAND / FE4S4_MACRO  passed through
# ============================================================================
set -uo pipefail
export PATH="$PATH:/usr/local/bin:/root/.local/bin"

TIMEOUT="${TIMEOUT:-21600}"                                   # 6 h hard cap
REPO="${FE4S4_REPO:-https://github.com/skharkov1246/ocm-research-diary}"
BRANCH="${FE4S4_BRANCH:-main}"
WORK="${FE4S4_WORKDIR:-/root/fe4s4}"

# --- dead-man #1: schedule poweroff FIRST, before anything can hang ---------
shutdown -P +$(( TIMEOUT / 60 + 30 )) 2>/dev/null || \
  sudo shutdown -P +$(( TIMEOUT / 60 + 30 )) 2>/dev/null || true

# --- dead-man #2: terminate on ANY exit (IMDSv2-aware) -----------------------
terminate() {
  echo "[dead-man] terminating instance to avoid leaks…"
  # forensic uploads first — the only trace if setup died early
  if [ -n "${S3_OUT:-}" ] && command -v aws >/dev/null 2>&1; then
    [ -f "${WORK}/setup.log" ] && aws s3 cp "${WORK}/setup.log" \
      "${S3_OUT%/}/fe4s4_setup.log" || true
    [ -f /var/log/cloud-init-output.log ] && aws s3 cp \
      /var/log/cloud-init-output.log \
      "${S3_OUT%/}/fe4s4_cloud-init-output.log" || true
  fi
  TOK="$(curl -sX PUT --max-time 3 \
        -H 'X-aws-ec2-metadata-token-ttl-seconds: 60' \
        http://169.254.169.254/latest/api/token || true)"
  IID="$(curl -s --max-time 3 -H "X-aws-ec2-metadata-token: ${TOK}" \
        http://169.254.169.254/latest/meta-data/instance-id || true)"
  REG="$(curl -s --max-time 3 -H "X-aws-ec2-metadata-token: ${TOK}" \
        http://169.254.169.254/latest/meta-data/placement/region || true)"
  if [ -n "${IID}" ] && command -v aws >/dev/null 2>&1; then
    aws ec2 terminate-instances --instance-ids "${IID}" --region "${REG}" \
      || shutdown -P now || sudo shutdown -P now || true
  else
    shutdown -P now || sudo shutdown -P now || true
  fi
}
trap terminate EXIT

upload() {  # upload <local> <basename>  (best-effort, never fatal)
  [ -n "${S3_OUT:-}" ] && [ -f "$1" ] \
    && aws s3 cp "$1" "${S3_OUT%/}/$2" || true
}

# --- setup: packages + python stack (every network step under timeout) ------
mkdir -p "${WORK}"
timeout 15m dnf install -y python3-pip git gcc gcc-c++ >"${WORK}/setup.log" 2>&1 || true
timeout 30m python3 -m pip install -q --upgrade pip           >>"${WORK}/setup.log" 2>&1 || true
timeout 30m python3 -m pip install -q awscli boto3 numpy scipy pyscf openfermion \
                                                              >>"${WORK}/setup.log" 2>&1
PIP_OK=$?
# pennylane is only the fallback JW backend — best-effort, non-fatal
timeout 15m python3 -m pip install -q pennylane               >>"${WORK}/setup.log" 2>&1 || true
if [ "${PIP_OK}" -ne 0 ]; then
  echo "SETUP-FAILED: pip stack install" >>"${WORK}/setup.log"
  upload "${WORK}/setup.log" fe4s4_setup.log
  exit 1                                   # trap -> terminate
fi

# --- clone the repo (branch parametric, default main) ------------------------
if ! timeout 10m git clone --depth 1 -b "${BRANCH}" "${REPO}" "${WORK}/repo" \
       >>"${WORK}/setup.log" 2>&1; then
  echo "SETUP-FAILED: git clone ${REPO} @ ${BRANCH}" >>"${WORK}/setup.log"
  upload "${WORK}/setup.log" fe4s4_setup.log
  exit 1                                   # trap -> terminate
fi
cd "${WORK}/repo"

# --- environment for the job -------------------------------------------------
mkdir -p /root/scratch
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-$(nproc)}"
export PYSCF_TMPDIR=/root/scratch
export PYSCF_MAX_MEMORY="${PYSCF_MAX_MEMORY:-$(( $(getconf _PHYS_PAGES) * $(getconf PAGE_SIZE) / 1024 / 1024 * 6 / 10 ))}"
export FE4S4_STAGE_TIMEOUT="${FE4S4_STAGE_TIMEOUT:-5400}"
# leave 30 min of the shell cap for uploads + termination
export FE4S4_TOTAL_BUDGET="${FE4S4_TOTAL_BUDGET:-$(( TIMEOUT - 1800 ))}"

# --- run ---------------------------------------------------------------------
echo "[run] fe4s4_casscf_aws.py  (timeout $(( TIMEOUT - 900 ))s, ${OMP_NUM_THREADS} threads, ${PYSCF_MAX_MEMORY} MB)"
timeout "$(( TIMEOUT - 900 ))" python3 -u calc/fe4s4_casscf_aws.py \
  > calc/fe4s4_casscf.log 2>&1
RC=$?
echo "[done] rc=${RC}"
tail -n 40 calc/fe4s4_casscf.log || true

# --- upload results (the python job saves incrementally, so even a timed-out
# --- run leaves a valid, honest partial results JSON) ------------------------
if [ -n "${S3_OUT:-}" ]; then
  echo "[upload] -> ${S3_OUT}"
  upload calc/fe4s4_casscf_results.json  fe4s4_casscf_results.json
  upload calc/fe4s4_casscf.log           fe4s4_casscf.log
  upload calc/fe4s4_casscf_integrals.npz fe4s4_casscf_integrals.npz
  upload calc/fe4s4_cubane_aws.xyz       fe4s4_cubane_aws.xyz
  upload "${WORK}/setup.log"             fe4s4_setup.log
else
  echo "[IMPORTANT] no S3_OUT set — copy calc/fe4s4_casscf_results.json (+log,"
  echo "            +npz) off NOW; the dead-man switch terminates this instance"
  echo "            on exit. Re-run with: S3_OUT=s3://bucket/prefix $0"
fi
# trap fires here -> instance self-terminates
