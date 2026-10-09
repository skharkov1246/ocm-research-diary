# Скептик CNG: второй год, деньги, инфляция, простой. Расширение cng_model.run
import sys,io,contextlib,math,json
src=open('/tmp/claude-0/-home-user-ocm-research-diary/50ead64b-47e2-5bdd-b34b-09f6667b6ca9/scratchpad/cng/cng_model.py').read()
src=src.split("out={'units'")[0]
g={}; exec(src,g)
S=g['S'];npv=g['npv'];irr=g['irr'];SCM=g['SCM_PER_MMSCF'];MM=g['MMBTU_PER_SCM'];FX=1330.0
def m(Q,p,mode='mult',delay=0,ramp1=1.0,usd_infl=0.0,ngn_dep=0.0,idx=1.0,idx_lag=1,outage=None,buyer_loss=None,
      wc_days=0,top=False,years=10,tax=0.34,debt=None,theft=0.0,carry=True,resid=0.0):
    station=p['st_pkg']*Q**0.7*p['cmult']; inlet=Q*SCM
    sold_d=inlet*p['sell']*p['util']; loads=sold_d/p['usable_scm']
    trip=2*p['km']/p['speed']+p['cyc_extra']
    tractors=math.ceil(loads*trip/p['drive_h']*1.15); ncust=math.ceil(sold_d/p['cust_scm'])
    trailers=math.ceil((loads*12/24+loads*trip/24+ncust)*1.15)
    fleet=trailers*p['trailer']+tractors*p['tractor']+ncust*p['cust_kit']
    capex=station+fleet+p['dev']
    n=years+delay+1; cf=[0.0]*n
    for i in range(delay+1): cf[i]-=capex/(delay+1)
    rows=[];prev=0;loss=0
    for t in range(1,years+1):
        T=t+delay; gg=(1-p['decl'])**(t-1)
        u=p['util']*gg if mode=='mult' else min(p['util'],gg)
        if t==1: u*=ramp1
        if buyer_loss and t in buyer_loss: u*=buyer_loss[t]
        av=outage.get(t,1.0) if outage else 1.0
        sold=inlet*p['sell']*u*av*365
        fx=FX*(1+ngn_dep)**T; k=max(0,T-idx_lag)
        price=p['price']*(1+idx*((1+ngn_dep)**k-1))
        rev=sold*price/fx/1e6*(1-p['bad'])
        infl=(1+usd_infl)**T
        gasv=(gg if top else u*av)
        gas=inlet*gasv*365*MM*p['gas']/1e6
        km=sold/p['usable_scm']*2*p['km']
        truck=km*p['opx_km']/1e6*infl
        fixed=(p['staff']+p['sec']+p['overhead']+tractors*2*6000/1e6+capex*(p['tor']+p['ins']))*infl
        e=rev-gas-truck-fixed-theft*p['trailer']
        dep=capex/7 if t<=7 else 0
        ti=e-dep
        if carry:
            ti2=ti-loss
            if ti2<0: loss=-ti2; tx=0
            else: loss=0; tx=ti2*tax
        else: tx=max(0,ti*tax)
        dwc=wc_days/365*(rev-prev); prev=rev
        c=e-tx-dwc-(tractors*p['tractor']*infl if t==6 else 0)
        if t==years: c+=wc_days/365*rev+resid*trailers*p['trailer']
        cf[T]+=c
        rows.append(dict(t=t,rev=round(rev,2),ebitda=round(e,2),cfads=round(e-tx,2)))
    r=dict(capex=round(capex,2),y1=rows[0]['ebitda'],y2=rows[1]['ebitda'],y3=rows[2]['ebitda'],npv15=round(npv(.15,cf),2),npv20=round(npv(.20,cf),2),irr=(round(irr(cf)*100,1) if irr(cf) is not None else None))
    if debt:
        sh,rate,ten=debt; D=capex*sh; ann=D*rate/(1-(1+rate)**-ten)
        r['dscr_y1_3']=[round(rows[i]['cfads']/ann,2) for i in range(3)]; r['ds']=round(ann,2)
    return r
def be(Q,p,key='npv15',**kw):
    lo,hi=100,4000
    for _ in range(50):
        x=(lo+hi)/2
        if m(Q,dict(p,price=x),**kw)[key]>0: hi=x
        else: lo=x
    return round(hi), round(hi/FX/MM,2)
pc=dict(S['base']);pc.update(km=60,util=0.75,staff=0.55*1.8,sec=0.3*1.6,overhead=0.2*1.5)
cases={'5c_620':(5,dict(pc,price=620)),'5c_750':(5,dict(pc,price=750)),'1_high':(1,S['high'])}
mild=dict(mode='min',delay=1,ramp1=0.8,usd_infl=0.03,ngn_dep=0.08,idx=1.0,idx_lag=1,outage={t:11/12 for t in range(1,11)},wc_days=60,top=True)
bad=dict(mild,ramp1=0.6,outage={t:(0.75 if t==2 else 11/12) for t in range(1,11)},buyer_loss={3:0.75,4:0.75},wc_days=90,theft=1.0)
mild_noidx=dict(mild,idx=0.5)
tests=[('A0 модель как есть',{}),
 ('A1 спрос-лимит min(util,дебит) (за проект)',dict(mode='min')),
 ('A2 стройка+разрешения 1 г',dict(delay=1)),
 ('A3 разгон 1-го года 60%',dict(ramp1=0.6)),
 ('A4 USD-инфляция опекса 3%/г',dict(usd_infl=0.03)),
 ('A5 найра -8%/г, цена индексируется с лагом 1 г',dict(ngn_dep=0.08)),
 ('A6 найра -8%/г, индексация 50%',dict(ngn_dep=0.08,idx=0.5)),
 ('A7 найра -8%/г, без индексации',dict(ngn_dep=0.08,idx=0.0)),
 ('A8 простой оператора 3 мес. на 2-й год (нет DoP)',dict(outage={2:0.75})),
 ('A9 1 мес. простоя каждый год',dict(outage={t:11/12 for t in range(1,11)})),
 ('A10 уход якорного клиента (-25%) на 2-3 г',dict(buyer_loss={2:0.75,3:0.75})),
 ('A11 дебиторка 90 дней',dict(wc_days=90)),
 ('A12 take-or-pay по газу',dict(top=True)),
 ('A13 потеря 1 каскада/г (ДТП/угон)',dict(theft=1.0)),
 ('B1 реалистичный 2-й год',mild),('B1b то же, индексация 50%',mild_noidx),('B2 плохой 2-й год',bad)]
out={}
for nm,(Q,p) in cases.items():
    for lab,kw in tests:
        r=m(Q,p,debt=(0.6,0.14,7),**kw); out[nm+'|'+lab]=r; print(nm,lab,r)
    print()
base5=pc
for lab,kw in [('A0',{}),('A1 min',dict(mode='min')),('B1',mild),('B1b',mild_noidx),('B2',bad)]:
    print('BE 5contract',lab,be(5,base5,**kw),be(5,base5,'npv20',**kw),' | 1_high',be(1,S['high'],**kw))
p3=dict(S['base']);p3.update(km=60,util=0.75,staff=0.55*1.4,sec=0.3*1.3,overhead=0.2*1.25)
for lab,kw in [('A0',{}),('B1',mild),('B2',bad)]: print('BE 3contract',lab,be(3,p3,**kw))
json.dump(out,open('/tmp/claude-0/-home-user-ocm-research-diary/50ead64b-47e2-5bdd-b34b-09f6667b6ca9/scratchpad/cng_y2/y2_out.json','w'),ensure_ascii=False,indent=1)
