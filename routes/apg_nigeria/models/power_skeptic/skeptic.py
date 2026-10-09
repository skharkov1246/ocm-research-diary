import sys
sys.path.insert(0,'/tmp/claude-0/-home-user-ocm-research-diary/50ead64b-47e2-5bdd-b34b-09f6667b6ca9/scratchpad/power')
import io,contextlib
with contextlib.redirect_stdout(io.StringIO()):
    import captive as C
FX=1330.0
def run(p0,lag=0,extra_opex=0.0,extra_capex=0.0,gas_avail=1.0,top=False,wc_days=0,lf_div=1.0,fx_dep=0.0,index=1.0,years=10,ramp=1.0):
    p=dict(p0); q=1
    mw_g=C.MW_PER_MMSCFD*q; mw_n=mw_g*(1-C.PARASITIC)
    gen=mw_g*1000*p['usd_kw']*p['red']/1e6
    capex=gen+p['cond_fix']+p['cond_var']*q+p['line_km']*p['usd_km']/1e6+p['subst']+p['dev']+extra_capex
    LF=p['LF']/lf_div
    cf=[0.0]*(years+lag+1)
    # capex spread evenly over construction years (lag years), else at t=0
    if lag==0: cf[0]-=capex
    else:
        for i in range(lag): cf[i]-=capex/lag
    prev_rev=0; e1=None
    for t in range(1,years+1):
        g=gas_avail*(1-p['decl'])**(t-1)
        u=min(LF,g)*(ramp if t==1 else 1)
        sold=8760*mw_n*p['avail']*u*(1-p['loss'])/1e3
        fxt=FX*(1+fx_dep)**(t-1); price_ngn=p['price_ngn']*(1+(index)*(((1+fx_dep)**(t-1))-1))
        rev=sold*price_ngn/fxt*(1-p['baddebt'])
        gasvol=(gas_avail*(1-p['decl'])**(t-1)) if top else u
        gas_c=C.MMBTU_PER_MMSCFD_Y*q*gasvol*(1 if top else p['avail'])*p['gas']/1e6
        opex=gas_c+p['tor']*capex+p['staff']+p['sec']+p['ins']*capex+extra_opex
        e=rev-opex
        if e1 is None: e1=e
        tax=max(0,0.34*(e-capex/years))
        wc=wc_days/365*(rev-prev_rev); prev_rev=rev
        c=e-tax-wc+(p['salv']*gen if t==years else 0)+(wc_days/365*rev if t==years else 0)
        idx=t+lag if lag>0 else t
        if lag>0: idx=lag-1+t+1-1+0  # first operating year right after construction
        cf[t+ (lag-1 if lag>0 else 0)]+=c
    return dict(capex=round(capex,2),e1=round(e1,2),npv15=round(C.npv(0.15,cf),2),npv20=round(C.npv(0.20,cf),2),irr=(round(C.irr(cf)*100,1) if C.irr(cf) is not None else None),pb=round(capex/e1,1) if e1>0 else None)
b=C.S['base']
cases={
 'base_reproduce':dict(),
 'build_12m(capex t0, rev from y2)':dict(lag=1),  # check
}
print('base',run(b))
# lag=1 means capex at t=0, first op year t=1 -> same as base; lag=2: capex split t0,t1, op from t=2
print('build 18-24m (capex split y0/y1, op from y2)',run(b,lag=2))
print('build 12m + 6m ramp 50% in y1',run(b,lag=1,ramp=0.5))
print('missing G&A/community/regulatory +0.3',run(b,extra_opex=0.3))
print('take-or-pay full gas',run(b,top=True))
print('distribution net cluster +1.0 capex',run(b,extra_capex=1.0))
print('WC 60d receivables',run(b,wc_days=60))
print('gas 0.5 of declared',run(b,gas_avail=0.5))
print('gas 0.25 of declared',run(b,gas_avail=0.25))
print('demand /1.3 pre-contract',run(b,lf_div=1.3))
print('naira dep 10%/y, index 50%',run(b,fx_dep=0.10,index=0.5))
print('naira dep 10%/y, index 100% (CPI~FX)',run(b,fx_dep=0.10,index=1.0))
print('COMBINED soft: lag2 + opex+0.3 + top + distnet+1.0 + WC60',run(b,lag=2,extra_opex=0.3,top=True,extra_capex=1.0,wc_days=60))
print('COMBINED + demand/1.3',run(b,lag=2,extra_opex=0.3,top=True,extra_capex=1.0,wc_days=60,lf_div=1.3))
print('COMBINED, 12m build only',run(b,lag=1,ramp=0.5,extra_opex=0.3,top=True,extra_capex=1.0,wc_days=60))
h=C.S['high']; print('HIGH combined lag2',run(h,lag=2,extra_opex=0.3,top=True,extra_capex=0.5,wc_days=60))
# required PPA price for NPV15=0 under combined soft
for pr in range(300,460,10):
    p=dict(b);p['price_ngn']=pr
    r=run(p,lag=2,extra_opex=0.3,top=True,extra_capex=1.0,wc_days=60)
    print(pr,r['npv15'],r['npv20'],r['irr'])
