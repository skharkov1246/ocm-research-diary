# payload: расширение вакансионного волкана — 4d/5d + главная подгруппа
# (Zr Nb Mo W Sn Pb), dft+corr на металл, 3 параллельно × OMP 5 (r7i.4xlarge).
run_m() { M=$1; OMP_NUM_THREADS=5 timeout 150m python3 -u routes/chemloop_vacancy.py dft $M >> $LOGF 2>&1; OMP_NUM_THREADS=5 timeout 150m python3 -u routes/chemloop_vacancy.py corr $M >> $LOGF 2>&1; echo "[done] $M" >> $LOGF; }
export -f run_m
echo "Zr Nb Mo W Sn Pb" | tr ' ' '\n' | xargs -P 3 -I{} bash -c 'run_m "$@"' _ {}
python3 -u routes/chemloop_vacancy.py merge >> $LOGF 2>&1
