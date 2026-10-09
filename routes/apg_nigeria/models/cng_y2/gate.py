src=open('/tmp/claude-0/-home-user-ocm-research-diary/50ead64b-47e2-5bdd-b34b-09f6667b6ca9/scratchpad/cng/cng_model.py').read().split("out={'units'")[0]
g={}; exec(src,g); npv=g['npv'];SCM=g['SCM_PER_MMSCF'];MM=g['MMBTU_PER_SCM']
def gate(Q,pr,util,cm,gas,st,se,ov,tor=0.045,ins=0.015,decl=0.06,dev=0.7,years=10,tax=0.34,sell=0.88,fixdep=False,delay=0,ramp=1.0,out=1.0,infl=0.0,mode='mult'):
    capex=2.0*Q**0.7*cm+dev; mm=Q*SCM*sell*365*MM/1e6
    cf=[0.0]*(years+delay+1)
    for i in range(delay+1): cf[i]-=capex/(delay+1)
    for t in range(1,years+1):
        T=t+delay; gg=(1-decl)**(t-1); u=util*gg if mode=='mult' else min(util,gg)
        if t==1: u*=ramp
        u*=out; f=(1+infl)**T
        e=mm*u*pr-Q*SCM*u*365*MM*gas/1e6-(st+se+ov+capex*(tor+ins))*f
        dep=capex/7 if (not fixdep or t<=7) else 0
        cf[T]+=e-max(0,(e-dep)*tax)
    return npv(.15,cf)
def be(**kw):
    lo,hi=0.5,40
    for _ in range(50):
        x=(lo+hi)/2
        if gate(pr=x,**kw)>0: hi=x
        else: lo=x
    return round(hi,2)
k=1.6
base=dict(Q=5,util=0.65,cm=2.6,gas=0.6,st=0.35*k,se=0.25*k,ov=0.12*k)
high=dict(Q=5,util=0.8,cm=2.4,gas=0.25,st=0.30*k,se=0.20*k,ov=0.10*k)
for nm,c in [('base',base),('high',high)]:
    print(nm,'as model',be(**c),'fixdep',be(fixdep=True,**c),'y2 realistic',be(fixdep=True,delay=1,ramp=0.8,out=11/12,infl=0.03,**c),'y2+min',be(fixdep=True,delay=1,ramp=0.8,out=11/12,infl=0.03,mode='min',**c))
