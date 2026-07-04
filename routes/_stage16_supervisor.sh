#!/usr/bin/env bash
# Супервизор Этапа 16 в песочнице: ждёт исход текущего оркестратора; при смерти
# джоба — спасает частичные JSON из S3 в ветку (инстанс-роль write-only, читать
# S3 инстанс не может, поэтому resume едет через git) и перезапускает.
cd /home/user/ocm-research-diary
S3=s3://alpha-o-results-097743207937/alpha-o/ocm16
AWSQ=/usr/local/bin/awsq
CUR_PID=${1:-}

salvage() {
  local got=0
  for f in ocm_mnw_emb_spins.json ocm_mnw_emb_ch4_geom.json ocm_mnw_emb_c2h6_geom.json \
           ocm_mnw_emb_ch4_profile.json ocm_mnw_emb_c2h6_profile.json; do
    $AWSQ s3 cp $S3/$f routes/$f >/dev/null 2>&1 && got=1
  done
  if [ $got -eq 1 ]; then
    python3 -c "
import json,glob
for p in glob.glob('routes/ocm_mnw_emb_*.json'): json.load(open(p))
print('salvaged JSONs valid')" || return 1
    git add routes/ocm_mnw_emb_*.json 2>/dev/null
    git commit -q -m "wip: ocm Stage 16 — salvage partial AWS state for cross-instance resume

Co-Authored-By: Claude Fable 5 <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01YA1aCie7ocDd1QDZSaeQMq" 2>/dev/null && git push -q 2>/dev/null
    echo "[sup] salvaged state committed"
  fi
}

if [ -n "$CUR_PID" ]; then
  while kill -0 $CUR_PID 2>/dev/null; do sleep 60; done
fi
# оркестратор умер; если финал уже в S3 — выходим
$AWSQ s3api head-object --bucket alpha-o-results-097743207937 --key alpha-o/ocm16/ocm_mnw_emb_final.json >/dev/null 2>&1 && { echo "[sup] FINAL present"; exit 0; }
for att in $(seq 1 8); do
  echo "[sup] attempt $att: salvage + relaunch"
  salvage
  OCM16_ITYPE=r7i.4xlarge OCM16_MAX_MIN=750 /opt/qc-venv/bin/python -u routes/ocm_aws_stage16.py; rc=$?
  [ $rc -eq 0 ] && { echo "[sup] SUCCESS"; exit 0; }
  OCM16_ITYPE=r7i.2xlarge OCM16_MAX_MIN=900 /opt/qc-venv/bin/python -u routes/ocm_aws_stage16.py; rc=$?
  [ $rc -eq 0 ] && { echo "[sup] SUCCESS"; exit 0; }
  echo "[sup] attempt $att failed (rc=$rc); sleep 180"; sleep 180
done
echo "[sup] EXHAUSTED"
