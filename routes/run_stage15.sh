#!/usr/bin/env bash
# Супервизор Этапа 15 OCM: доводит недостающие стадии по обоим субстратам.
# Идемпотентен — можно перезапускать после смерти контейнера; geom резюмится
# с последней сохранённой точки, готовые стадии определяются по маркерам в JSON.
set -u
cd "$(dirname "$0")/.."
PY=/opt/qc-venv/bin/python
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-2}

run_substrate() {
  local s=$1
  if ! $PY - "$s" <<'EOF'
import json, sys
try:
    g = json.load(open(f"routes/ocm_mnw_{sys.argv[1]}_geom.json"))
    assert "dft_barrier_kcal" in g and len(g["points"]) == 8
except Exception:
    raise SystemExit(1)
EOF
  then
    $PY -u routes/ocm_mnw_selectivity.py geom "$s" >> "routes/ocm_mnw_${s}_geom.log" 2>&1 || return 1
  fi
  if ! $PY - "$s" <<'EOF'
import json, sys
g = json.load(open(f"routes/ocm_mnw_{sys.argv[1]}_geom.json"))
raise SystemExit(0 if g.get("polish_done") else 1)
EOF
  then
    $PY -u routes/ocm_mnw_selectivity.py polish "$s" >> "routes/ocm_mnw_${s}_geom.log" 2>&1 || return 1
  fi
  if ! $PY - "$s" <<'EOF'
import json, sys
try:
    e = json.load(open(f"routes/ocm_mnw_{sys.argv[1]}_energy.json"))
    assert "nevpt2_barrier_kcal" in e
except Exception:
    raise SystemExit(1)
EOF
  then
    $PY -u routes/ocm_mnw_selectivity.py energy "$s" >> "routes/ocm_mnw_${s}_energy.log" 2>&1 || return 1
  fi
  echo "SUPERVISOR_DONE_$s"
}

run_substrate ch4 &
P1=$!
run_substrate c2h6 &
P2=$!
wait $P1; R1=$?
wait $P2; R2=$?
if [ $R1 -eq 0 ] && [ $R2 -eq 0 ]; then
  $PY -u routes/ocm_mnw_selectivity.py merge && echo SUPERVISOR_ALL_DONE
else
  echo "SUPERVISOR_FAILED ch4=$R1 c2h6=$R2"
fi
