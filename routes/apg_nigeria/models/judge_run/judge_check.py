import sys,io,contextlib
sys.path.insert(0,'/tmp/claude-0/-home-user-ocm-research-diary/50ead64b-47e2-5bdd-b34b-09f6667b6ca9/scratchpad/deal')
with contextlib.redirect_stdout(io.StringIO()):
    import dealcore as D
M=D.M; run=D.run; P=('P',)
DEAL=dict(dep=0.10,idx=1.0,lag=0,usd_inf=0.03,up=11/12,ramp=0.8,wc=60,delay=1)
STRUCT=dict(DEAL,up=1.0,wc=30)
ESC=dict(DEAL,wc=30)   # escrow, no DoP
o2=dict(bad=0.02)
def r(pr,kw,ov=None,**ex):
    o=dict(p_dir=pr); o.update(ov or {})
    x=run(5,P,ov=o,debt=(0.6,0.14,7),**ex,**kw); return {k:x.get(k) for k in ('capex','y1','y2','npv15','npv20','irr','payback','min_dscr')}
def be(kw,key,ov=None,**ex):
    lo,hi=150,900
    for _ in range(40):
        m=(lo+hi)/2; o=dict(p_dir=m); o.update(ov or {})
        if run(5,P,ov=o,**ex,**kw)[key]<0: lo=m
        else: hi=m
    return round(hi)
print('377 STRUCT EDTI', r(377,STRUCT,o2))
print('377 STRUCT noEDTI dep10', r(377,STRUCT,o2,edti=False,depyrs=10))
print('377 ESC(noDoP) EDTI', r(377,ESC,o2))
print('377 ESC(noDoP) noEDTI', r(377,ESC,o2,edti=False,depyrs=10))
print('377 DEAL EDTI', r(377,DEAL))
print('275 STRUCT noEDTI', r(275,STRUCT,o2,edti=False,depyrs=10))
print('BE ESC noEDTI', be(ESC,'npv15',o2,edti=False,depyrs=10), be(ESC,'npv20',o2,edti=False,depyrs=10))
print('BE STRUCT noEDTI', be(STRUCT,'npv15',o2,edti=False,depyrs=10), be(STRUCT,'npv20',o2,edti=False,depyrs=10))
# gas price stop at 350 STRUCT EDTI: bisect gas $/MMBtu for NPV20=0 and NPV15=0 (opfee kept 0.15)
for key in ('npv20','npv15'):
    for pr,e,dy in ((350,True,5),(377,False,10)):
        lo,hi=0.0,10.0
        for _ in range(40):
            m=(lo+hi)/2
            x=run(5,P,ov=dict(p_dir=pr,bad=0.02,gas=m),edti=e,depyrs=dy,**STRUCT)[key]
            if x>=0: lo=m
            else: hi=m
        print('gas stop',pr,'EDTI',e,key,'=0 at gas $',round(lo,2),'/MMBtu (+opfee 0.15/Mscf)')
# ToP whole site, earnout 0.30, LF 0.6, block, cut at 350 STRUCT EDTI
base=r(350,STRUCT,o2); print('base350',base)
print('ToP whole site', r(350,STRUCT,o2,gas_extra=0.51))
print('earnout 0.30', r(350,STRUCT,dict(bad=0.02,opfee=0.45)))
print('LF0.6', r(350,STRUCT,dict(bad=0.02,LF=0.6)))
print('cm2.8', r(350,STRUCT,dict(bad=0.02,cm=2.8)))
print('avail0.6 noEDTI', r(377,STRUCT,dict(bad=0.02,avail=0.6),edti=False,depyrs=10))
print('avail0.5 noEDTI 350', r(350,STRUCT,dict(bad=0.02,avail=0.5),edti=False,depyrs=10))
print('delay2 377 noEDTI', r(377,dict(STRUCT,delay=2),o2,edti=False,depyrs=10))
