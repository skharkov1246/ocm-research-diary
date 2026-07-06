#!/usr/bin/env bash
# Супервизор Этапа 16 в песочнице (переживает рестарты контейнера — просто
# перезапустить без аргументов): пока живёт ocm-stage16-инстанс — ждёт и
# периодически спасает частичные JSON из S3 в ветку; после смерти инстанса —
# если FINAL в S3, выходит; иначе salvage + перезапуск оркестратора.
# Инстанс-роль write-only в S3, читать S3 инстанс не может — resume между
# поколениями едет через git (свежий клон ветки содержит спасённый прогресс).
cd /home/user/ocm-research-diary
S3=s3://alpha-o-results-097743207937/alpha-o/ocm16
AWSQ=/usr/local/bin/awsq

final_present() {
  $AWSQ s3api head-object --bucket alpha-o-results-097743207937 \
    --key alpha-o/ocm16/ocm_mnw_emb_final.json >/dev/null 2>&1
}

my_instance() {
  $AWSQ ec2 describe-instances \
    --filters "Name=tag:Name,Values=ocm-stage16" \
              "Name=instance-state-name,Values=pending,running" \
    --query 'Reservations[0].Instances[0].InstanceId' --output text 2>/dev/null
}

salvage() {
  local got=0
  for f in ocm_mnw_emb_spins.json ocm_mnw_emb_ch4_geom.json ocm_mnw_emb_c2h6_geom.json \
           ocm_mnw_emb_ch4_profile.json ocm_mnw_emb_c2h6_profile.json; do
    $AWSQ s3 cp $S3/$f routes/$f >/dev/null 2>&1 && got=1
  done
  [ $got -eq 0 ] && return 0
  python3 - <<'PY' || return 1
import json, glob
for p in glob.glob('routes/ocm_mnw_emb_*.json'):
    json.load(open(p))
PY
  git add routes/ocm_mnw_emb_*.json 2>/dev/null
  if ! git diff --cached --quiet 2>/dev/null; then
    git commit -q -m "wip: ocm Stage 16 — salvage partial AWS state for cross-instance resume

Co-Authored-By: Claude Fable 5 <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01YA1aCie7ocDd1QDZSaeQMq" && git push -q 2>/dev/null
    echo "[sup] salvaged state committed"
  fi
}

for att in $(seq 1 200); do
  if final_present; then echo "[sup] FINAL present — done"; exit 0; fi
  IID=$(my_instance)
  if [ -n "$IID" ] && [ "$IID" != "None" ]; then
    echo "[sup] instance $IID alive; salvage + wait"
    salvage
    sleep 600
    continue
  fi
  echo "[sup] no instance; attempt $att: salvage + relaunch"
  salvage
  OCM16_ITYPE=r7i.4xlarge OCM16_MAX_MIN=750 /opt/qc-venv/bin/python -u routes/ocm_aws_stage16.py; rc=$?
  if final_present; then echo "[sup] FINAL present — done"; exit 0; fi
  OCM16_ITYPE=r7i.2xlarge OCM16_MAX_MIN=900 /opt/qc-venv/bin/python -u routes/ocm_aws_stage16.py; rc=$?
  if final_present; then echo "[sup] FINAL present — done"; exit 0; fi
  echo "[sup] attempt $att no final (rc=$rc); sleep 180"
  sleep 180
done
echo "[sup] EXHAUSTED"
