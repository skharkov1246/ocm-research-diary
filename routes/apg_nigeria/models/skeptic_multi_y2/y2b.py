import sys
sys.path.insert(0,'/tmp/claude-0/-home-user-ocm-research-diary/50ead64b-47e2-5bdd-b34b-09f6667b6ca9/scratchpad/skeptic_multi_y2')
import io,contextlib
with contextlib.redirect_stdout(io.StringIO()):
    import y2
run=y2.run; P=('P',)
for Q in (5,15):
  print(Q,'USD-тариф (idx1 lag0), USD-инфл3, простой 1мес, рамп80, WC60', run(Q,P,dep=0.10,idx=1,lag=0,usd_inf=0.03,up=11/12,ramp=0.8,wc=60,debt=(0.6,0.14,7)))
  print(Q,'USD-тариф, без простоя, рамп 100', run(Q,P,dep=0.10,idx=1,lag=0,usd_inf=0.03,wc=60))
  # безубыточный тариф ₦ при B0 (индексация 100% лаг 1)
  lo,hi=200,800
  for _ in range(40):
    m=(lo+hi)/2
    if run(Q,P,dep=0.10,idx=1,usd_inf=0.03,up=11/12,ramp=0.8,wc=60,ov=dict(p_dir=m))['npv15']<0: lo=m
    else: hi=m
  print(Q,'тариф безубыточности NPV15 при B0, ₦/кВт·ч',round(hi))
  lo,hi=200,1200
  for _ in range(40):
    m=(lo+hi)/2
    if run(Q,P,dep=0.10,idx=0.5,usd_inf=0.03,up=11/12,ramp=0.8,wc=60,ov=dict(p_dir=m))['npv15']<0: lo=m
    else: hi=m
  print(Q,'тариф безубыточности NPV15 при B1, ₦/кВт·ч',round(hi))
# эквивалент дизеля: ₦/л -> ₦/кВт·ч топливо при 0.30 л/кВт·ч
for d in (1100,1410,2000,3277): print('дизель',d,'₦/л ->',round(d*0.30),'₦/кВт·ч топливо')
