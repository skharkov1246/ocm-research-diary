# payload: металл-free анти-BEP панель I2 (O NH O2 CH2), каждый в своём timeout
run_med() { Mx=$1; OMP_NUM_THREADS=4 timeout 90m python3 -u routes/antibep_carbene.py run $Mx >> $LOGF 2>&1; echo "[done] $Mx" >> $LOGF; }
export -f run_med
echo "O NH O2 CH2" | tr ' ' '\n' | xargs -P 2 -I{} bash -c 'run_med "$@"' _ {}
python3 -u routes/antibep_carbene.py merge >> $LOGF 2>&1
