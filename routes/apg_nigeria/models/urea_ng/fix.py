import importlib.util,io,contextlib
spec=importlib.util.spec_from_file_location('m','/tmp/claude-0/-home-user-ocm-research-diary/50ead64b-47e2-5bdd-b34b-09f6667b6ca9/scratchpad/urea_ng/model.py')
m=importlib.util.module_from_spec(spec)
with contextlib.redirect_stdout(io.StringIO()): spec.loader.exec_module(m)
def irr(cf):
    lo,hi=-0.5,1.0
    f=lambda r:m.npv(r,cf)
    if f(lo)*f(hi)>0: return None
    for _ in range(200):
        mid=(lo+hi)/2
        if f(lo)*f(mid)<=0: hi=mid
        else: lo=mid
    return mid
def run(q,sc,price=None,decl=None,close=True,years=20,tax=0.30):
    p=dict(m.S[sc]); capex=m.CAPEX[q][sc]
    if decl is not None: p['decl']=decl
    pr=price if price is not None else p['price']; g=p['gas']
    gj_day=q*m.MMBTU_HHV_PER_MMSCF*m.GJ_LHV_PER_MMBTU_HHV
    tpd=gj_day/p['gj_t']; tpy0=tpd*p['days']/1e3; mmbtu_t=q*m.MMBTU_HHV_PER_MMSCF/tpd
    cf=[-capex/p['build']]*p['build']; stop=None; neg_after=0
    for t in range(1,years+1):
        f=(1-p['decl'])**(t-1); load=f if f>=0.6 else 0
        if load==0 and stop is None: stop=t
        tpy=tpy0*load; rev=tpy*pr/1e3; var=tpy*(mmbtu_t*g+p['cons'])/1e3
        fix=m.STAFF[q][sc]+p['sec']*m.SEC[q]+p['tor']*capex+p['ins']*capex+p['ga']
        e=rev-var-fix
        if load==0:
            neg_after+=e
            if close: e=0
        dep=capex/10 if t<=10 else 0
        cf.append(e-max(0,tax*(e-dep)))
    r=irr(cf)
    return dict(npv15=round(m.npv(0.15,cf),1),npv20=round(m.npv(0.20,cf),1),irr=None if r is None else round(r*100,1),stop=stop,fixed_after_stop=round(neg_after,1))
for q in (1,5,15):
    for sc in ('opt','mid','pes'):
        a=run(q,sc,close=False); b=run(q,sc,close=True)
        print(q,sc,'as-is',a,'| closed',b)
