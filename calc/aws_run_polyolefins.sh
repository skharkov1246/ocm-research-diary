#!/usr/bin/env bash
# ============================================================================
# AWS user-data wrapper for the polyolefins chiral stereo-DDE job
# (routes/polyolefins_chiral_aws.py), polyolefins track / GitHub issue #11.
# Modeled 1:1 on calc/aws_run_fe4s4.sh (the freshest AWS pattern in the repo).
#
# Designed to run as EC2 USER-DATA (root, no tty, fresh AL2023 box): it clones
# the repo itself, installs the stack, runs the job, uploads results to S3 and
# self-terminates. It can equally be run by hand on an instance.
#
# Dead-man layering (same philosophy as aws_orchestrate.py / aws_run_fe4s4.sh):
#   1. `shutdown -P +N` scheduled FIRST — fires even if everything below hangs;
#      with InstanceInitiatedShutdownBehavior=terminate a poweroff TERMINATES.
#   2. terminate() trap on EXIT — explicit terminate via IMDSv2 (launch with
#      HttpTokens=required; plain IMDSv1 curl would fail), poweroff fallback.
#   3. The sandbox-side reaper in aws_orchestrate.py (kills by tag).
#
# Progress upload: a background loop pushes the incremental results JSON and
# the log tail to S3 every PROGRESS_EVERY seconds, so a run can be monitored
# (and salvaged) WITHOUT logging into the box — the python job saves
# atomically after every stage, so the uploaded JSON is always valid.
#
# Env (all optional except S3_OUT if you want the results to survive):
#   TIMEOUT              total seconds for the whole job, default 21600 (6 h)
#   POLY_REPO            git repo URL (default: this project)
#   POLY_BRANCH          branch to clone, default main
#   S3_OUT               s3://bucket/prefix for results upload
#   PROGRESS_EVERY       progress-upload period in seconds, default 600
#   POLY_STAGE_TIMEOUT   per-stage cap passed to the python job (default 4800)
#   POLY_TOTAL_BUDGET    python-side wall budget (default TIMEOUT-1800)
#   POLY_CC_POINTS / POLY_TS_CC / POLY_FACES / POLY_MF / POLY_XC /
#   POLY_OPT_MAXSTEPS / POLY_CAS / POLY_AVAS_THRESHOLD / POLY_DEGEN_TOL
#                        passed through to the python job (see its docstring)
#
# Sizing note (account limit: 16 vCPU total!):
#   c7i.4xlarge (16 vCPU, 32 GB) — uses the WHOLE quota, nothing else can run;
#   c7i.2xlarge (8 vCPU, 16 GB)  — slower; if it doesn't fit TIMEOUT, cut
#   POLY_CC_POINTS to "2.8,2.0" (the DDE pair) or raise TIMEOUT.
# ============================================================================
set -uo pipefail
export PATH="$PATH:/usr/local/bin:/root/.local/bin"

TIMEOUT="${TIMEOUT:-21600}"                                   # 6 h hard cap
REPO="${POLY_REPO:-https://github.com/skharkov1246/ocm-research-diary}"
BRANCH="${POLY_BRANCH:-main}"
WORK="${POLY_WORKDIR:-/root/poly}"
PROGRESS_EVERY="${PROGRESS_EVERY:-600}"

# --- dead-man #1: schedule poweroff FIRST, before anything can hang ---------
shutdown -P +$(( TIMEOUT / 60 + 30 )) 2>/dev/null || \
  sudo shutdown -P +$(( TIMEOUT / 60 + 30 )) 2>/dev/null || true

# --- dead-man #2: terminate on ANY exit (IMDSv2-aware) -----------------------
terminate() {
  echo "[dead-man] terminating instance to avoid leaks…"
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
    && aws s3 cp --only-show-errors "$1" "${S3_OUT%/}/$2" || true
}

# --- setup: packages + python stack (every network step under timeout) ------
mkdir -p "${WORK}"
timeout 15m dnf install -y python3-pip git gcc gcc-c++ >"${WORK}/setup.log" 2>&1 || true
timeout 30m python3 -m pip install -q --upgrade pip           >>"${WORK}/setup.log" 2>&1 || true
timeout 30m python3 -m pip install -q awscli boto3 numpy scipy pyscf geometric \
                                                              >>"${WORK}/setup.log" 2>&1
PIP_OK=$?
if [ "${PIP_OK}" -ne 0 ]; then
  echo "SETUP-FAILED: pip stack install" >>"${WORK}/setup.log"
  upload "${WORK}/setup.log" poly_setup.log
  exit 1                                   # trap -> terminate
fi

# --- clone the repo (branch parametric, default main) ------------------------
if ! timeout 10m git clone --depth 1 -b "${BRANCH}" "${REPO}" "${WORK}/repo" \
       >>"${WORK}/setup.log" 2>&1; then
  echo "SETUP-FAILED: git clone ${REPO} @ ${BRANCH}" >>"${WORK}/setup.log"
  upload "${WORK}/setup.log" poly_setup.log
  exit 1                                   # trap -> terminate
fi
cd "${WORK}/repo"

# --- environment for the job -------------------------------------------------
mkdir -p /root/scratch
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-$(nproc)}"
export PYSCF_TMPDIR=/root/scratch
export PYSCF_MAX_MEMORY="${PYSCF_MAX_MEMORY:-$(( $(getconf _PHYS_PAGES) * $(getconf PAGE_SIZE) / 1024 / 1024 * 6 / 10 ))}"
export POLY_STAGE_TIMEOUT="${POLY_STAGE_TIMEOUT:-4800}"
# leave 30 min of the shell cap for uploads + termination
export POLY_TOTAL_BUDGET="${POLY_TOTAL_BUDGET:-$(( TIMEOUT - 1800 ))}"
# pass-through knobs (exported only if the caller set them)
for V in POLY_CC_POINTS POLY_TS_CC POLY_FACES POLY_MF POLY_XC \
         POLY_OPT_MAXSTEPS POLY_CAS POLY_AVAS_THRESHOLD POLY_DEGEN_TOL; do
  [ -n "$(eval "echo \${${V}:-}")" ] && export "${V?}"
done

RESJSON=routes/polyolefins_chiral_results.json
LOG=routes/polyolefins_chiral.log

# --- preflight smoke: never burn hours on a config typo ----------------------
echo "[smoke] preflight (--smoke --no-scf)"
if ! timeout 5m python3 routes/polyolefins_chiral_aws.py --smoke --no-scf \
       > routes/polyolefins_chiral_smoke.log 2>&1; then
  echo "SMOKE-FAILED — aborting before the heavy run"
  tail -n 30 routes/polyolefins_chiral_smoke.log || true
  upload routes/polyolefins_chiral_smoke.log poly_smoke.log
  upload "${WORK}/setup.log"                poly_setup.log
  exit 1                                   # trap -> terminate
fi

# --- progress uploader: incremental JSON + log tail every PROGRESS_EVERY s ---
if [ -n "${S3_OUT:-}" ]; then
  (
    while sleep "${PROGRESS_EVERY}"; do
      upload "${RESJSON}" polyolefins_chiral_results.json
      tail -n 200 "${LOG}" > /tmp/poly_progress.log 2>/dev/null || true
      upload /tmp/poly_progress.log poly_progress.log
    done
  ) &
  PROG_PID=$!
else
  PROG_PID=""
fi

# --- run ---------------------------------------------------------------------
echo "[run] polyolefins_chiral_aws.py  (timeout $(( TIMEOUT - 900 ))s, ${OMP_NUM_THREADS} threads, ${PYSCF_MAX_MEMORY} MB)"
timeout "$(( TIMEOUT - 900 ))" python3 -u routes/polyolefins_chiral_aws.py \
  > "${LOG}" 2>&1
RC=$?
echo "[done] rc=${RC}"
[ -n "${PROG_PID}" ] && kill "${PROG_PID}" 2>/dev/null || true
tail -n 40 "${LOG}" || true

# --- upload results (the python job saves incrementally + atomically, so even
# --- a timed-out run leaves a valid, honest partial results JSON) ------------
if [ -n "${S3_OUT:-}" ]; then
  echo "[upload] -> ${S3_OUT}"
  upload "${RESJSON}"                            polyolefins_chiral_results.json
  upload "${LOG}"                                polyolefins_chiral.log
  upload routes/polyolefins_chiral_smoke.log     poly_smoke.log
  upload routes/polyolefins_chiral_model.xyz     polyolefins_chiral_model.xyz
  upload routes/polyolefins_chiral_ts_si.xyz     polyolefins_chiral_ts_si.xyz
  upload routes/polyolefins_chiral_ts_re.xyz     polyolefins_chiral_ts_re.xyz
  upload "${WORK}/setup.log"                     poly_setup.log
else
  echo "[IMPORTANT] no S3_OUT set — copy routes/polyolefins_chiral_results.json"
  echo "            (+log, +xyz) off NOW; the dead-man switch terminates this"
  echo "            instance on exit. Re-run with: S3_OUT=s3://bucket/prefix $0"
fi
# trap fires here -> instance self-terminates
