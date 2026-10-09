import sys,json,importlib.util
sys.argv=['x']
spec=importlib.util.spec_from_file_location('m','/tmp/claude-0/-home-user-ocm-research-diary/50ead64b-47e2-5bdd-b34b-09f6667b6ca9/scratchpad/urea_ng/model.py')
import io,contextlib
buf=io.StringIO()
with contextlib.redirect_stdout(buf):
    m=importlib.util.module_from_spec(spec); spec.loader.exec_module(m)
S,CAPEX,STAFF,SEC=m.S,m.CAPEX,m.STAFF,m.SEC
npv,irr=m.npv,m.irr
def run(q,sc,price=None,gas=None,decl=None,netback=0,tax=0.30,years=20,close=True,capex=None):
    p=dict(S[sc]); 
    if decl is not None: p['decl']=decl
    cap=capex if capex is not None else CAPEX[q][sc]
    pr=(price if price is not None else p['price'])-netback; g=gas if gas is not None else p['gas']
    gj=q*1050*1.055*0.9; tpd=gj/p['gj_t']; tpy0=tpd*p['days']/1e3; mmbtu_t=q*1050/tpd
    cf=[-cap/p['build']]*p['build']; e1=None; stop=None; es=[]
    for t in range(1,years+1):
        f=(1-p['decl'])**(t-1); load=f if f>=0.6 else 0
        if load==0 and close:
            stop=stop or t; cf.append(0); continue
        tpy=tpy0*load; rev=tpy*pr/1e3; var=tpy*(mmbtu_t*g+p['cons'])/1e3
        fix=STAFF[q][sc]+p['sec']*SEC[q]+(p['tor']+p['ins'])*cap+p['ga']
        e=rev-var-fix; es.append(e)
        if e1 is None: e1=e
        dep=cap/10 if t<=10 else 0
        cf.append(e-max(0,tax*(e-dep)))
    ir=irr(cf)
    return dict(q=q,sc=sc,cap=cap,price=pr,e1=round(e1,2),stop_year=stop,npv15=round(npv(.15,cf),1),npv20=round(npv(.2,cf),1),irr=None if ir is None else round(ir*100,1),sumE=round(sum(es),1))
for q in (1,5,15):
    for sc in ('opt','mid','pes'):
        a=run(q,sc,close=False); b=run(q,sc)
        print('orig-like',a); print('closed  ',b)
print('--- 15 opt with decline 3%, 5%')
for d in (0.03,0.05): print(run(15,'opt',decl=d))
print('--- netback FOB-> gate: minus $30/40 inland+port+small lot')
for q in (5,15):
    for nb in (30,60): print(run(q,'opt',netback=nb))
print('--- 15 opt at Argus late-Aug-2026 ~$325 FOB')
print(run(15,'opt',price=325)); print(run(5,'opt',price=325))
print('--- tax 30%+4% dev levy')
print(run(15,'opt',tax=0.34))
print('--- capex threshold from key_unknowns: 5MMscf/d capex 120 @450 opt/mid')
print(run(5,'opt',price=450,capex=120)); print(run(5,'mid',price=450,capex=120))
print('--- 1 MMscf/d: $25M capex, $450')
print(run(1,'opt',price=450,capex=25)); print(run(1,'mid',price=450,capex=25))
print('=== capex threshold for NPV15=0 at $450')
def thr(q,sc,price=450,**kw):
    lo,hi=1,1000
    for _ in range(60):
        c=(lo+hi)/2
        if run(q,sc,price=price,capex=c,**kw)['npv15']>0: lo=c
        else: hi=c
    return round(c,1)
for q in (1,5,15):
    for sc in ('opt','mid'):
        c=thr(q,sc); kt={1:dict(opt=14.1,mid=12.2),5:dict(opt=70.6,mid=60.9),15:dict(opt=211.9,mid=182.8)}[q][sc]
        print(q,sc,'capex<=',c,'M  $/t/y=',round(c*1e3/kt))
print('=== combined realism for 15 opt: decl 3%, netback 40')
print(run(15,'opt',decl=0.03,netback=40))
print(run(15,'opt',decl=0.03,netback=40,tax=0.34))
print(run(15,'opt',price=325,decl=0.03))
