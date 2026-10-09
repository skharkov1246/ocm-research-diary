# Скептик (второй год и деньги) для multi_model.py: энергия-только и полная линия, 5 и 15 MMscf/д.
# Берёт год-1 разбивку из multi_model.site и строит поток с: девальвацией найры, индексацией, USD-инфляцией опекса,
# простоем апстрима, рампом, дебиторкой, падением цены LPG, долгом (DSCR).
import sys,io,contextlib,json
sys.path.insert(0,'/tmp/claude-0/-home-user-ocm-research-diary/50ead64b-47e2-5bdd-b34b-09f6667b6ca9/scratchpad/deal')
with contextlib.redirect_stdout(io.StringIO()):
    import multi_model as M
def run(Q,cfg,sc='base',dep=0.0,idx=0.0,lag=1,usd_inf=0.0,up=1.0,ev=None,ramp=1.0,wc=0,lpg_drift=0.0,delay=0,
        tariff_cut=None,debt=None,years=10,ov=None,extra_fee=0.0,gas_extra=0.0,edti=True,depyrs=5):
    r=M.site(Q,sc,cfg,verbose=True,ov=ov or {},extra_rev_per_mscf=extra_fee); C=r['capex']; rows=r['rows']
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
        rf=row.get('rev_operator_fee',0)*u
        rev=rp+rc+rl+rf
        var=0
        gas=row['gas']*(1+0.02)**(T-1)                          # GSA индексирован по US-индексу ~2%/г [В: индексация; Д: ставка]
        fixed=(row['opex']-row['gas']+gas_extra)*(1+usd_inf)**(T-1)       # прочий опекс (ToP газа не падает при простое)
        e=rev-gas-fixed
        dep_=C/depyrs if t<=depyrs else 0
        tx=max(0,(e-dep_)*0.34); ed=(min(tx,0.05*C) if t<=5 else 0) if edti else 0
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
