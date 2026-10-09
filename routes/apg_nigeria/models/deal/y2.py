# Скептик (второй год и деньги) для multi_model.py: энергия-только и полная линия, 5 и 15 MMscf/д.
# Берёт год-1 разбивку из multi_model.site и строит поток с: девальвацией найры, индексацией, USD-инфляцией опекса,
# простоем апстрима, рампом, дебиторкой, падением цены LPG, долгом (DSCR).
import sys,io,contextlib,json
sys.path.insert(0,'/tmp/claude-0/-home-user-ocm-research-diary/50ead64b-47e2-5bdd-b34b-09f6667b6ca9/scratchpad/deal')
with contextlib.redirect_stdout(io.StringIO()):
    import multi_model as M
def run(Q,cfg,sc='base',dep=0.0,idx=0.0,lag=1,usd_inf=0.0,up=1.0,ev=None,ramp=1.0,wc=0,lpg_drift=0.0,delay=0,
        tariff_cut=None,debt=None,years=10,ov=None):
    r=M.site(Q,sc,cfg,verbose=True,ov=ov or {}); C=r['capex']; rows=r['rows']
    cf=[-C]+[0.0]*delay; prev=0; out=[]
    for i,row in enumerate(rows):
        t=i+1; T=t+delay
        fx=(1+dep)**T
        k=max(0,T-lag); pidx=(1+idx*((1+dep)**k-1))/fx      # USD-ценность найрового тарифа
        if tariff_cut and t>=tariff_cut[0]: pidx*=tariff_cut[1]
        u=up*(ev.get(t,1.0) if ev else 1.0)*(ramp if t==1 else 1.0)
        rp=row.get('rev_power',0)*pidx*u
        rc=row.get('rev_cng',0)*pidx*u                         # CNG тоже в найре
        rl=row.get('rev_lpg',0)*((1+lpg_drift)**(t-1))*u        # LPG в $-паритете к импорту/Dangote
        rev=rp+rc+rl
        var=0
        gas=row['gas']*(1+0.02)**(T-1)                          # GSA индексирован по US-индексу ~2%/г [В: индексация; Д: ставка]
        fixed=(row['opex']-row['gas'])*(1+usd_inf)**(T-1)       # прочий опекс (ToP газа не падает при простое)
        e=rev-gas-fixed
        dep_=C/5 if t<=5 else 0
        tx=max(0,(e-dep_)*0.34); ed=min(tx,0.05*C) if t<=5 else 0
        w=wc/365*(rev-prev); prev=rev
        c=e-tx+ed-w+(wc/365*rev if t==years else 0)
        cf.append(c); out.append(dict(t=t,rev=round(rev,2),ebitda=round(e,2),cfads=round(e-tx+ed,2)))
    res=dict(capex=C,y1=out[0]['ebitda'],y2=out[1]['ebitda'],y5=out[4]['ebitda'],npv15=round(M.npv(.15,cf),1),npv20=round(M.npv(.20,cf),1),
             irr=(round(M.irr(cf)*100,1) if M.irr(cf) is not None else None))
    acc=-C;pb=None
    for i,c in enumerate(cf[1:],1):
        if pb is None and acc+c>=0 and c>0: pb=round(i-1+(-acc)/c,1)
        acc+=c
    res['payback']=pb
    if debt:
        sh,rate,ten=debt; D=C*sh; ann=D*rate/(1-(1+rate)**-ten)
        res['debt_service']=round(ann,2); res['dscr_y1_3']=[round(out[j]['cfads']/ann,2) for j in range(3)]
        res['min_dscr']=round(min(o['cfads'] for o in out[:ten])/ann,2)
    return res
P=('P',); F=('P','L','C')
cases=[
 ('R0 репродукция',{}),
 ('N8 найра -8%/г, без индексации',dict(dep=0.08)),
 ('N12 найра -12%/г, без индексации',dict(dep=0.12)),
 ('N16 найра -16%/г (CAGR 2016-26), без индексации',dict(dep=0.16)),
 ('I12 найра -12%/г, индексация 100% лаг 1 г',dict(dep=0.12,idx=1.0)),
 ('H12 найра -12%/г, индексация 50%',dict(dep=0.12,idx=0.5)),
 ('U3 USD-инфляция фикс.опекса 3%/г',dict(usd_inf=0.03)),
 ('O1 простой апстрима 1 мес/г (TNP/Forcados, общины)',dict(up=11/12)),
 ('O3 блокада 3 мес на 2-й год',dict(ev={2:0.75})),
 ('RP рамп 1-го года 70%',dict(ramp=0.7)),
 ('WC дебиторка 75 дней',dict(wc=75)),
 ('D1 задержка 1 г',dict(delay=1)),
 ('T2 ренегоциация тарифа -20% со 2-го года (дизель упал до ₦1 400/л)',dict(tariff_cut=(2,0.8))),
 ('L5 LPG дешевеет 5%/г (Dangote)',dict(lpg_drift=-0.05)),
]
REAL=dict(dep=0.10,idx=0.5,usd_inf=0.03,up=11/12,ramp=0.8,wc=60)
BAD=dict(dep=0.12,idx=0.0,usd_inf=0.03,up=11/12,ev={2:0.75},ramp=0.7,wc=75,delay=1,tariff_cut=(2,0.8),lpg_drift=-0.05)
GOOD=dict(dep=0.10,idx=1.0,usd_inf=0.03,up=11/12,ramp=0.8,wc=60)
cases+=[('B1 «год 2 реалистичный»: найра -10%/г, индексация 50%, USD-инфл 3%, простой 1 мес/г, рамп 80%, WC 60',REAL),
        ('B0 «год 2 с полной индексацией»: как B1, но индексация 100% лаг 1 г',GOOD),
        ('B2 «год 2 плохой»: найра -12% без индексации, +блокада, рамп 70%, WC 75, задержка 1, тариф -20%, LPG -5%/г',BAD)]
out={}
for Q in (5,15):
  for cfg in (P,F):
    for nm,kw in cases:
      k=f'{Q}|{"+".join(cfg)}|{nm}'
      out[k]=run(Q,cfg,debt=(0.6,0.14,7),**kw); print(k,out[k])
# EV база/пессимизм при B1
for Q in (5,15):
  b=run(Q,P,**REAL); p=run(Q,P,sc='pess',**REAL)
  out[f'EV50 {Q}|P|B1']=dict(npv15=round(.5*b['npv15']+.5*p['npv15'],1),npv20=round(.5*b['npv20']+.5*p['npv20'],1),base=b['npv15'],pess=p['npv15'])
  print('EV',Q,out[f'EV50 {Q}|P|B1'])
# порог: какая доля индексации нужна для NPV15=0 при -10%/г (5, P)
for Q in (5,15):
  for d in (0.08,0.10,0.12):
    lo,hi=0,1
    if run(Q,P,dep=d,idx=1,usd_inf=0.03,up=11/12,ramp=0.8,wc=60)['npv15']<0: out[f'idx_need {Q} dep{d}']='>1'; print(Q,d,'>1'); continue
    for _ in range(40):
      m=(lo+hi)/2
      if run(Q,P,dep=d,idx=m,usd_inf=0.03,up=11/12,ramp=0.8,wc=60)['npv15']<0: lo=m
      else: hi=m
    out[f'idx_need {Q} dep{d}']=round(hi,2); print('idx needed',Q,d,round(hi,2))
json.dump(out,open('/tmp/claude-0/-home-user-ocm-research-diary/50ead64b-47e2-5bdd-b34b-09f6667b6ca9/scratchpad/deal/y2_out_copy.json','w'),ensure_ascii=False,indent=1)
