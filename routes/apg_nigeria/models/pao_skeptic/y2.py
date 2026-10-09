# Second-year / money stress of scenario A (PAO Nigeria). Base = model.py scenario A
CH4=7300.0; MMBTU=365000*1.05
def run(price=3200,Y=0.35,olig=0.85,capex=12.0,gas=0.25,var=250,logi=150,staff=0.6,sec=0.3,maint=0.035,ins=0.015,
        avail=1.0,ramp=(1.0,1.0),decline=0.0,wc_days=0,tax=0.30,infl=0.0,r=0.15,life=12,build=2,par=0.12,verbose=False):
    t0=CH4*(1-par)*Y*28/32.08*olig
    v=-capex/build-capex/build/(1+r); cf=[]; wc_prev=0
    for y in range(life):
        f=avail*(ramp[y] if y<len(ramp) else 1.0)*(1-decline)**y
        t=t0*f; rev=t*price/1e6
        fixed=(staff+sec+capex*(maint+ins))*(1+infl)**y
        e=rev-MMBTU*f*gas/1e6-t*(var+logi)/1e6-fixed
        wc=rev*wc_days/365; dwc=wc-wc_prev; wc_prev=wc
        dep=capex/10 if y<10 else 0
        tx=max(0,e-dep)*tax
        c=e-tx-dwc
        if y==life-1: c+=wc
        cf.append((y+1,round(t),round(rev,2),round(e,2),round(c,2)))
        v+=c/(1+r)**(build+y)
    if verbose:
        for x in cf[:4]: print('  op-yr',x)
    return round(v,1),cf
def show(name,**k):
    n15,cf=run(r=0.15,**k); n20,_=run(r=0.20,**k)
    print(f"{name}: yr1 EBITDA {cf[0][3]}, yr2 EBITDA {cf[1][3]}, yr5 {cf[4][3]}, NPV15 {n15}, NPV20 {n20}")
show('A as published ($3200)')
show('A + tax 34% (CIT30+levy4)',tax=0.34)
OPS=dict(avail=0.85,ramp=(0.5,0.75),wc_days=90,tax=0.34,infl=0.03)
show('A + ops realism (avail.85, ramp 50/75, WC90d, tax34, infl3%)',**OPS)
show('A + ops + decline10%',decline=0.10,**OPS)
show('A + ops + decline15%',decline=0.15,**OPS)
show('A price EU $2100 (IMARC DE), no ops',price=2100)
show('A price EU $2100 + ops',price=2100,**OPS)
show('A price EU $2100 + ops + decline10',price=2100,decline=0.10,**OPS)
show('hexene-PAO EU $1450 (2100x2200/3200) + ops',price=1450,**OPS)
# breakeven price with ops+decline10 for NPV15=0
for dec in (0.0,0.10):
    lo,hi=1000,8000
    for _ in range(60):
        m=(lo+hi)/2
        if run(price=m,decline=dec,**OPS)[0]>0: hi=m
        else: lo=m
    print('breakeven price NPV15=0, decline',dec,round(m))
# scale 15 MMscf/d at EU price
s=15
show('15 MMscf/d, $3200, ops+decl10', capex=12*s**0.7, staff=0.6*s**0.4, sec=0.3*s**0.3, decline=0.10, **{**OPS})
