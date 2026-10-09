# Скептик: второй год, деньги, инфляция, простой. Расширение captive.py (база 1 MMscf/д)
import sys,io,contextlib
sys.path.insert(0,'/tmp/claude-0/-home-user-ocm-research-diary/50ead64b-47e2-5bdd-b34b-09f6667b6ca9/scratchpad/power')
with contextlib.redirect_stdout(io.StringIO()):
    import captive as C
FX=1330.0
def model(p0, delay=0, ramp1=1.0, usd_infl=0.0, ngn_dep=0.0, idx=1.0, idx_lag=1,
          outage=None, buyer_loss=None, wc_days=0, top=False, years=10, dscr_debt=None, ngn_cost_share=0.0):
    p=dict(p0); q=1
    mwg=C.MW_PER_MMSCFD; mwn=mwg*(1-C.PARASITIC)
    gen=mwg*1000*p['usd_kw']*p['red']/1e6
    capex=gen+p['cond_fix']+p['cond_var']+p['line_km']*p['usd_km']/1e6+p['subst']+p['dev']
    n=years+delay+1
    cf=[0.0]*n
    # стройка: капекс равными долями на годы 0..delay
    for i in range(delay+1): cf[i]-=capex/(delay+1)
    rows=[]; prev_rev=0
    for t in range(1,years+1):
        T=t+delay            # календарный год от FID
        g=(1-p['decl'])**(t-1)
        LF=p['LF']
        if buyer_loss and t in buyer_loss: LF*=buyer_loss[t]
        u=min(LF,g)*(ramp1 if t==1 else 1)
        av=p['avail']*(outage.get(t,1.0) if outage else 1.0)
        sold=8760*mwn*av*u*(1-p['loss'])/1e3
        fx=FX*(1+ngn_dep)**T
        # тариф в найре: индексация доли idx к найровой девальвации с лагом idx_lag лет
        k=max(0,T-idx_lag); price=p['price_ngn']*(1+idx*((1+ngn_dep)**k-1))
        rev=sold*price/fx*(1-p['baddebt'])
        gv=g if top else u
        gas=C.MMBTU_PER_MMSCFD_Y*gv*(1 if top else av)*p['gas']/1e6
        infl=(1+usd_infl)**T
        fixed=(p['tor']*capex+p['ins']*capex+p['staff']+p['sec'])*infl
        opex=gas+fixed
        e=rev-opex
        tax=max(0,0.34*(e-capex/years))
        wc=wc_days/365*(rev-prev_rev); prev_rev=rev
        c=e-tax-wc+(p['salv']*gen if t==years else 0)+(wc_days/365*rev if t==years else 0)
        cf[T]+=c
        rows.append(dict(t=t,GWh=round(sold,1),rev=round(rev,2),opex=round(opex,2),ebitda=round(e,2),cfads=round(e-tax,2)))
    acc=0;dpb=None
    for i,c in enumerate(cf):
        acc+=c/1.15**i
        if acc>=0 and dpb is None and i>0: dpb=i
    r=dict(capex=round(capex,2),y1=rows[0]['ebitda'],y2=rows[1]['ebitda'],npv15=round(C.npv(.15,cf),2),npv20=round(C.npv(.20,cf),2),
           irr=(round(C.irr(cf)*100,1) if C.irr(cf) is not None else None),disc_pb15_from_FID=dpb)
    if dscr_debt:
        share,rate,ten=dscr_debt
        D=capex*share; ann=D*rate/(1-(1+rate)**-ten)
        r['dscr_y1_y2_y3']=[round(rows[i]['cfads']/ann,2) for i in range(3)]
        r['debt_service']=round(ann,2)
    return r,rows
b=C.S['base']
cases=[
 ('A0 база (репродукция)',{}),
 ('A1 стройка+разрешения 2 г (капекс 3 транша), без задержек выручки после',dict(delay=2)),
 ('A2 первый год 60% нагрузки (Otakikpo)',dict(ramp1=0.6)),
 ('A3 USD-инфляция опекса 3%/г',dict(usd_infl=0.03)),
 ('A4 найра -8%/г, PPA индексирован к курсу с лагом 1 г',dict(ngn_dep=0.08,idx=1.0,idx_lag=1)),
 ('A5 найра -8%/г, индексация 50%',dict(ngn_dep=0.08,idx=0.5)),
 ('A6 найра -8%/г, без индексации',dict(ngn_dep=0.08,idx=0.0)),
 ('A7 блокада/вандализм ЛЭП: 3 мес. простоя на 2-й год',dict(outage={2:0.75})),
 ('A8 ежегодно 1 мес. простоя ЛЭП/общины',dict(outage={t:11/12 for t in range(1,11)})),
 ('A9 уход 1 из 4 покупателей на 2-3 г (на трубный газ)',dict(buyer_loss={2:0.75,3:0.75})),
 ('A10 дебиторка 75 дней',dict(wc_days=75)),
 ('A11 take-or-pay по газу',dict(top=True)),
]
comb=dict(delay=1,ramp1=0.6,usd_infl=0.03,ngn_dep=0.08,idx=1.0,idx_lag=1,outage={2:0.75},buyer_loss={3:0.75},wc_days=75,top=True)
comb_mild=dict(delay=1,ramp1=0.8,usd_infl=0.03,ngn_dep=0.08,idx=1.0,idx_lag=1,outage={t:11/12 for t in range(1,11)},wc_days=60,top=True)
cases+= [('B1 «второй год реалистичный»: стройка 2 г, рамп 80%, инфл 3%, найра -8% с лагом, 1 мес. простоя/г, WC 60, ToP',comb_mild),
         ('B2 «второй год плохой»: + 3 мес. блокада г.2, уход покупателя г.3, рамп 60%, WC 75',comb)]
import json
out={}
for nm,kw in cases:
    r,rows=model(b,dscr_debt=(0.6,0.14,8),**kw)
    out[nm]=r; print(nm,r)
# DSCR в найровом долге (InfraCredit-тип): ставка 20% в найре ~ в долларах при -8%/г: (1.20/1.08)-1=11%
print('naira debt equiv USD rate', round(1.20/1.08-1,3))
for nm,kw in [('A0',{}),('B1',comb_mild),('B2',comb)]:
    r,rows=model(b,dscr_debt=(0.6,0.11,10),**kw); print(nm,'naira-debt-equiv DSCR',r['dscr_y1_y2_y3'],r['debt_service'])
# высокий и низкий сценарии при B1
for k in ['high','low']:
    r,_=model(C.S[k],**comb_mild); print(k,'B1',r)
r,_=model(b); print('discounted payback 15% base (years from FID)',r['disc_pb15_from_FID'])
json.dump(out,open('/tmp/claude-0/-home-user-ocm-research-diary/50ead64b-47e2-5bdd-b34b-09f6667b6ca9/scratchpad/power_y2/y2_out.json','w'),ensure_ascii=False,indent=1)
