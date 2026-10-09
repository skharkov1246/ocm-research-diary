# independent check of standalone mining, 1 MMscf/d, no exit option, 10y
def npv(r,cf): return sum(c/(1+r)**t for t,c in enumerate(cf))
def irr(cf):
    lo,hi=-0.9,3
    if npv(lo,cf)*npv(hi,cf)>0: return None
    for _ in range(200):
        m=(lo+hi)/2
        if npv(lo,cf)*npv(m,cf)<=0: hi=m
        else: lo=m
    return m
def run(hp0,g,h,gen_kw=1700,cond=1.2,infra=2.6,dev=0.7,e=13.5,usd_th=13,landed=1.2,avail=0.88,pit=4.14,
        staff=0.35,sec=0.35,gas=0.5,tor=0.04,ins=0.012,years=10,salv=0.4,r=0.15):
    bop=4.64*gen_kw/1e3+cond+infra+dev
    th=pit*1e6/e; asic=th*usd_th*landed/1e6
    cap=bop+asic; cf=[-cap]
    eff=1.0; ac=asic
    for y in range(1,years+1):
        # monthly hashprice; halvings at month 10 and 58 from start
        rev=0
        for m in range((y-1)*12,y*12):
            hv=(1 if m>=10 else 0)+(1 if m>=58 else 0)
            hp=hp0*(1-g)**(m/12)*h**hv
            rev+=pit*1e3/(e*eff)*hp*avail*365/12/1e6  # PH * $/PH/d
        opex=1050*365*avail*gas/1e6+tor*cap+staff+sec+ins*cap+0.04*rev+0.03
        ebitda=rev-opex
        c=ebitda-max(0,0.34*(ebitda-bop/10-ac/4))
        if y==4: # replace ASIC: same MW, eff x0.75, $/TH x0.8
            eff*=0.75; ac=th/0.75*usd_th*0.8*landed/1e6; c-=ac
        if y==8: eff*=0.75; ac=th/0.75**2*usd_th*0.64*landed/1e6; c-=ac
        if y==years: c+=salv*bop
        cf.append(c)
    return cap,cf
def be(g,h,r=0.15,**k):
    lo,hi=1,2000
    for _ in range(60):
        m=(lo+hi)/2
        if npv(r,run(m,g,h,**k)[1])>0: hi=m
        else: lo=m
    return m
for name,g,h in [('base g10 h0.6',0.10,0.6),('USD-const hashprice',0.0,1.0),('opt g0 h0.7',0.0,0.7)]:
    for hp in (27.7,39.87,50):
        cap,cf=run(hp,g,h); print(name,hp,'cap',round(cap,2),'y1 EBITDA~',round(cf[1],2),'NPV15',round(npv(.15,cf),2),'IRR',irr(cf))
    print(name,'BE15',round(be(g,h),1),'BE15 lean BOP(gen 1000,infra1.5,cond0.8)',round(be(g,h,gen_kw=1000,infra=1.5,cond=0.8),1))
