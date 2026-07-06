# payload: прямая карбоксилация CH4+CO2->CH3COOH — фаза 1 (радикальная цепь)
OMP_NUM_THREADS=8 timeout 300m python3 -u routes/acetic_radical.py cc >> $LOGF 2>&1 &
P1=$!
OMP_NUM_THREADS=8 timeout 300m python3 -u routes/acetic_radical.py hat >> $LOGF 2>&1 &
P2=$!
wait $P1 $P2
echo "[done] cc+hat" >> $LOGF
python3 -u routes/acetic_radical.py merge >> $LOGF 2>&1
