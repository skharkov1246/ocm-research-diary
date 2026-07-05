# payload: MSA-цепь — add (стадии 1+2) и hat (стадия 3) параллельно, потом merge
OMP_NUM_THREADS=8 timeout 200m python3 -u routes/msa_chain.py add >> $LOGF 2>&1 &
P1=$!
OMP_NUM_THREADS=8 timeout 200m python3 -u routes/msa_chain.py hat >> $LOGF 2>&1 &
P2=$!
wait $P1 $P2
echo "[done] add+hat" >> $LOGF
python3 -u routes/msa_chain.py merge >> $LOGF 2>&1
