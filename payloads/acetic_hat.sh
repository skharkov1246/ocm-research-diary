# payload: только стадия hat (cc уже готова и лежит в ветке) + merge
OMP_NUM_THREADS=8 timeout 240m python3 -u routes/acetic_radical.py hat >> $LOGF 2>&1
echo "[done] hat" >> $LOGF
python3 -u routes/acetic_radical.py merge >> $LOGF 2>&1
