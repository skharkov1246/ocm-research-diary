import sys,io,contextlib
sys.path.insert(0,'/tmp/claude-0/-home-user-ocm-research-diary/50ead64b-47e2-5bdd-b34b-09f6667b6ca9/scratchpad/deal')
with contextlib.redirect_stdout(io.StringIO()):
    import dealcore as D
M=D.M; run=D.run; P=('P',)
DEAL=dict(dep=0.10,idx=1.0,lag=0,usd_inf=0.03,up=11/12,ramp=0.8,wc=60,delay=1)
STRUCT=dict(DEAL,up=1.0,wc=30)
def be(kw,key='npv15',ov=None,**ex):
    lo,hi=150,900
    for _ in range(40):
        m=(lo+hi)/2; o=dict(p_dir=m); o.update(ov or {})
        if run(5,P,ov=o,**ex,**kw)[key]<0: lo=m
        else: hi=m
    return round(hi)
o=dict(bad=0.02)
print('STRUCT BE with EDTI', be(STRUCT,'npv15',o), be(STRUCT,'npv20',o))
print('STRUCT BE no EDTI, dep 5y', be(STRUCT,'npv15',o,edti=False), be(STRUCT,'npv20',o,edti=False))
print('STRUCT BE no EDTI, dep 10y', be(STRUCT,'npv15',o,edti=False,depyrs=10), be(STRUCT,'npv20',o,edti=False,depyrs=10))
print('DEAL BE no EDTI dep10', be(DEAL,'npv15',o,edti=False,depyrs=10), be(DEAL,'npv20',o,edti=False,depyrs=10))
for pr in (330,350):
    for e,dy in ((True,5),(False,10)):
        r=run(5,P,ov=dict(p_dir=pr,bad=0.02),debt=(0.6,0.14,7),edti=e,depyrs=dy,**STRUCT)
        print(pr,'EDTI',e,'dep',dy,{k:r[k] for k in ('y1','y2','npv15','npv20','irr','payback','min_dscr')})
# equity IRR: 60% USD debt 14% 7y, drawn at t0, annuity from year 1 of operation (after delay)
def eq(pr,**kw):
    import dealcore as DD
    # recompute cash flows by calling run internals: replicate quickly
    r=M.site(5,'base',P,verbose=True,ov=dict(p_dir=pr,bad=0.02))
    return r
# simple equity calc using run cfads: reuse run with debt to get cfads list -> need cf; patch: compute via unlevered cf
import math
def eq_irr(pr,edti=True,depyrs=5,share=0.6,rate=0.14,ten=7,kw=STRUCT):
    # rebuild unlevered cf from run by monkeypatch: run returns only summary; replicate run loop
    r=M.site(5,'base',P,verbose=True,ov=dict(p_dir=pr,bad=0.02)); C=r['capex']; rows=r['rows']
    delay=kw['delay']; cf=[-C]+[0.0]*delay; prev=0
    Dd=C*share; ann=Dd*rate/(1-(1+rate)**-ten)
    eqcf=[-C*(1-share)]+[-Dd*rate]*delay   # interest during construction paid by equity
    for i,row in enumerate(rows):
        t=i+1; T=t+delay
        u=kw['up']*(kw['ramp'] if t==1 else 1.0)
        rev=row['rev_power']*u
        gas=row['gas']*(1.02)**(T-1); fixed=(row['opex']-row['gas'])*(1+kw['usd_inf'])**(T-1)
        e=rev-gas-fixed
        d_=C/depyrs if t<=depyrs else 0
        bal=Dd*(1+rate)**(t-1)-ann*((1+rate)**(t-1)-1)/rate if t<=ten else 0
        intr=bal*rate if t<=ten else 0
        tx=max(0,(e-d_-intr)*0.34); ed=(min(tx,0.05*C) if t<=5 else 0) if edti else 0
        w=kw['wc']/365*(rev-prev); prev=rev
        c=e-tx+ed-w+(kw['wc']/365*rev if t==10 else 0)
        eqcf.append(c-(ann if t<=ten else 0))
    return round(M.irr(eqcf)*100,1) if M.irr(eqcf) is not None else None, round(M.npv(0.20,eqcf),1), round(C*(1-share),1)
for pr in (300,330,350):
    print('equity IRR',pr, eq_irr(pr), 'noEDTI dep10', eq_irr(pr,edti=False,depyrs=10))
