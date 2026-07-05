# payload: I3 ди-Fe MMO-Q — geom (FM 2S=8) + FM-профиль (tractable CAS(12,13)).
# AFM-синглет (2S=0) = CAS(18,15) ~25M детерминантов = DMRG-grade, вне статич. стека
# (задокументировано как потолок; сюда НЕ запускаем — не жжём бокс на неподъёмном).
DIFE_BASIS=def2-svp OMP_NUM_THREADS=16 timeout 200m python3 -u routes/dife_mmo.py geom >> $LOGF 2>&1
echo "[done] geom" >> $LOGF
DIFE_BASIS=def2-svp OMP_NUM_THREADS=32 timeout 300m python3 -u routes/dife_mmo.py profile 8 >> $LOGF 2>&1
echo "[done] profile-FM" >> $LOGF
python3 -u routes/dife_mmo.py merge >> $LOGF 2>&1
