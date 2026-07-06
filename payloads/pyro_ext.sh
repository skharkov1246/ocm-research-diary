# payload: волна 2 пиролиз-волкана — 8 новых пар решётки активный×растворитель
# + чистый пере-прогон Cu-Bi с расширенной лестницей спинов (NSPIN=3, суффикс _v2).
run_p() { P=$1; OMP_NUM_THREADS=4 timeout 150m python3 -u routes/pyrolysis_descriptor.py dft $P >> $LOGF 2>&1; OMP_NUM_THREADS=4 timeout 150m python3 -u routes/pyrolysis_descriptor.py corr $P >> $LOGF 2>&1; echo "[done] $P" >> $LOGF; }
export -f run_p
echo "Pt-Sn Pd-In Cu-Sn Fe-Bi Ni-Ga Co-Bi Cu-In Pd-Sn" | tr ' ' '\n' | xargs -P 4 -I{} bash -c 'run_p "$@"' _ {}
# Cu-Bi redo: 3 спина, отдельные файлы pyrolysis_v2_Cu-Bi.json
PYRO_NSPIN=3 PYRO_SUFFIX=_v2 OMP_NUM_THREADS=8 timeout 150m python3 -u routes/pyrolysis_descriptor.py dft Cu-Bi >> $LOGF 2>&1
PYRO_NSPIN=3 PYRO_SUFFIX=_v2 OMP_NUM_THREADS=8 timeout 150m python3 -u routes/pyrolysis_descriptor.py corr Cu-Bi >> $LOGF 2>&1
echo "[done] Cu-Bi-v2" >> $LOGF
python3 -u routes/pyrolysis_descriptor.py merge >> $LOGF 2>&1
