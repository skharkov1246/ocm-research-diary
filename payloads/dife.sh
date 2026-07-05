# payload: I3 ди-Fe MMO-Q — geom (FM), затем профили AFM(0) и FM(8) параллельно
DIFE_BASIS=def2-svp OMP_NUM_THREADS=16 timeout 200m python3 -u routes/dife_mmo.py geom >> $LOGF 2>&1
echo "[done] geom" >> $LOGF
DIFE_BASIS=def2-svp OMP_NUM_THREADS=16 timeout 240m python3 -u routes/dife_mmo.py profile 0 >> $LOGF 2>&1 &
P1=$!
DIFE_BASIS=def2-svp OMP_NUM_THREADS=16 timeout 240m python3 -u routes/dife_mmo.py profile 8 >> $LOGF 2>&1 &
P2=$!
wait $P1 $P2
echo "[done] profiles" >> $LOGF
python3 -u routes/dife_mmo.py merge >> $LOGF 2>&1
