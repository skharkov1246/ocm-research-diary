# payload: MSA add-стадия с симметричным радикал-CAS (MSA_RADCAS=1) — чистый
# коррелированный вердикт β-scission vs HAT (hat-число уже есть от ffm7)
MSA_RADCAS=1 OMP_NUM_THREADS=16 timeout 340m python3 -u routes/msa_chain.py add >> $LOGF 2>&1
echo "[done] add-radcas" >> $LOGF
