# v1 10-летняя модель + падение дебита газа (как в captive.py: 10/5/0 %/г) — безубыточный hashprice
import sys; sys.path.insert(0,'../compute')
import io,contextlib
with contextlib.redirect_stdout(io.StringIO()):
    import mining10 as M
from mining_ai import npv
def run(hp0,asic,sc,decl,years=10,refresh=4):
    p=M.S2[sc]; a=M.ASIC[asic]
    mw0=M.MW_GROSS*p['derate']*p['parasit']
    bop=M.MW_GROSS*1000*p['gen_kw']/1e6+p['cond']+p['infra']+p['dev']
    e=a['e']; price=a['usd_th'][sc]
    th=mw0*1e6/e; asic_capex=th*price*p['landed']/1e6
    cf=[-(bop+asic_capex)]
    for y in range(years):
        if y>0 and y%refresh==0:
            e*=0.75; price*=0.8; th=mw0*(1-decl)**y*1e6/e; asic_capex=th*price*p['landed']/1e6; cf[-1]-=asic_capex
        g=(1-decl)**y
        rev=0
        for m in range(12):
            mm=y*12+m; f=1.0
            for h in M.HALVINGS:
                if mm>=h: f*=p['halv']
            rev+=th*min(1,g/((1-decl)**(4*(y//4))))/1000*hp0*(1-p['erosion'])**(mm/12)*f*30.42*p['uptime']/1e6
        ct=bop+asic_capex
        opex=M.MMBTU_D*365*g*p['uptime']*p['gas']/1e6+p['tor']*ct+p['staff']+p['sec']+p['ins']*ct+(p['comm']+p['pool'])*rev+0.03
        eb=rev-opex; dep=asic_capex/refresh+bop/10
        cf.append(eb-max(0,M.TAX*(eb-dep)))
    cf[-1]+=p['salv']*bop
    return npv(0.15,cf)
def be(sc,decl):
    lo,hi=5,600
    for _ in range(60):
        m=(lo+hi)/2
        if run(m,'A_13.5JTH',sc,decl)<0: lo=m
        else: hi=m
    return round(m,1)
for sc,d in [('pess',0.10),('base',0.05),('opt',0.0),('base',0.0),('pess',0.0)]:
    print(sc,d,'BE15',be(sc,d),'NPV15@39.87',round(run(39.87,'A_13.5JTH',sc,d),2))
