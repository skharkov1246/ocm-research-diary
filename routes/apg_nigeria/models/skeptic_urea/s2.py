import importlib.util,io,contextlib
P='/tmp/claude-0/-home-user-ocm-research-diary/50ead64b-47e2-5bdd-b34b-09f6667b6ca9/scratchpad/urea_ng/model.py'
spec=importlib.util.spec_from_file_location('m',P); m=importlib.util.module_from_spec(spec)
with contextlib.redirect_stdout(io.StringIO()): spec.loader.exec_module(m)
def irr(cf):
    lo,hi=-0.5,1.0;f=lambda r:m.npv(r,cf)
    if f(lo)*f(hi)>0: return None
    for _ in range(200):
        mid=(lo+hi)/2
        if f(lo)*f(mid)<=0: hi=mid
        else: lo=mid
    return mid
def run(q,sc,price,decl=None,ramp=(1,),wc_m=0,outage=0,ta_every=0,ta_days=0,ta_cost=0,catfill=0,cost_infl=0.0,years=20,tax=0.30,capex_x=1.0,close=True):
    p=dict(m.S[sc]); capex=m.CAPEX[q][sc]*capex_x
    if decl is not None: p['decl']=decl
    g=p['gas']; gj=q*m.MMBTU_HHV_PER_MMSCF*m.GJ_LHV_PER_MMBTU_HHV
    tpd=gj/p['gj_t']; mm=q*m.MMBTU_HHV_PER_MMSCF/tpd
    cf=[-capex/p['build']]*p['build']; cf[-1]-=catfill
    eb=[]; prev_rev=0
    for t in range(1,years+1):
        f=(1-p['decl'])**(t-1); load=f if f>=0.6 else 0
        r=ramp[t-1] if t-1<len(ramp) else 1
        days=p['days']-outage - (ta_days if ta_every and t%ta_every==0 else 0)
        tpy=tpd*days/1e3*load*r
        rev=tpy*price/1e3; var=tpy*(mm*g+p['cons'])/1e3
        fix=(m.STAFF[q][sc]+p['sec']*m.SEC[q]+p['tor']*capex+p['ins']*capex+p['ga'])*(1+cost_infl)**(t-1)
        fix+= ta_cost if ta_every and t%ta_every==0 else 0
        e=rev-var-fix
        if load==0 and close: e=0
        wc=wc_m/12*(rev-prev_rev); prev_rev=rev
        dep=capex/10 if t<=10 else 0
        cf.append(e-max(0,tax*(e-dep))-wc)
        eb.append(round(e,1))
    cf[-1]+=wc_m/12*prev_rev
    r=irr(cf)
    return dict(npv15=round(m.npv(.15,cf),1),npv20=round(m.npv(.20,cf),1),irr=None if r is None else round(r*100,1),y1=eb[0],y2=eb[1],y5=eb[4])
base=dict(ramp=(0.65,0.85),wc_m=2,outage=15,ta_every=4,ta_days=21,ta_cost=3.0,catfill=0.0)
for q,sc in ((15,'opt'),(5,'opt'),(15,'mid')):
  for pr in (480,430,400,325,295):
    a=run(q,sc,pr); b=run(q,sc,pr,**base); c=run(q,sc,pr,decl=0.03,**base)
    print(q,sc,pr,'model',a,'| y2-real',b,'| +decl3%',c)
# cost inflation 3%/yr USD on fixed with flat price
print('15opt 480 infl3',run(15,'opt',480,cost_infl=0.03,**base))
print('15opt 480 capex x1.3',run(15,'opt',480,capex_x=1.3,**base))
# price needed for NPV15=0 with y2-real and decl 3%
def be(q,sc,**kw):
    lo,hi=100,3000
    for _ in range(50):
        mid=(lo+hi)/2
        if run(q,sc,mid,**kw)['npv15']>0: hi=mid
        else: lo=mid
    return round(mid)
for q,sc in ((15,'opt'),(5,'opt'),(1,'opt'),(15,'mid')):
    print('BE',q,sc,be(q,sc),be(q,sc,**base),be(q,sc,decl=0.03,**base))
# debt service 15 opt: 60% debt, 18%, 8 yrs, 2 grace
cap=250; D=0.6*cap; r=0.18; n=8
ann=D*r/(1-(1+r)**-n)
print('annuity',round(ann,1))
for pr in (480,400,325):
    x=run(15,'opt',pr,**base); print(pr,'DSCR y1',round(x['y1']/ann,2),'y2',round(x['y2']/ann,2))
