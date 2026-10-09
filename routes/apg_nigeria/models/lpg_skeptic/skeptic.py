# skeptic re-run of add-on LPG case, 15 MMscf/d; same formulas as lpg_model.run2 + stresses
def npv(c,r): return sum(x/(1+r)**i for i,x in enumerate(c))
def irr(c):
    lo,hi=-0.99,3.0
    if npv(c,lo)*npv(c,hi)>0: return None
    for _ in range(200):
        m=(lo+hi)/2
        if npv(c,lo)*npv(c,m)<=0: hi=m
        else: lo=m
    return m
def run(capex,y=3.0,price=600,LF=0.8,decl=0.08,shrink=1.0,Q=15,infl=0.0,yr2_LF=None,yr2_price=None,price_from3=None,debt=0.0,rate=0.15,tenor=7,label=''):
    staff,sec,misc=0.35,0.25,0.2
    cfs=[-capex*0.4,-capex*0.6]; dep=capex/5; eb=[]; ds=[]
    D=capex*debt; ann=D*rate/(1-(1+rate)**-tenor) if D>0 else 0
    for yr in range(12):
        f=LF*(1-decl)**yr
        if yr==1 and yr2_LF is not None: f=yr2_LF*(1-decl)**yr
        p=price
        if yr==1 and yr2_price is not None: p=yr2_price
        if yr>=2 and price_from3 is not None: p=price_from3
        lpg=y*Q*365*f; cond=0.2*y*Q*365*f
        rev=(lpg*p+cond*400)/1e6
        gas=shrink*(lpg+cond)*47.3/1e6+shrink*0.03*1.05*Q*1000*365*f/1e6
        fixed=(staff+sec+misc+capex*0.06)*(1+infl)**yr
        e=rev-gas-fixed; t=max(0,(e-(dep if yr<5 else 0))*0.34); ed=min(t,0.05*capex) if yr<5 else 0
        cfs.append(e-t+ed); eb.append(e)
        if yr<tenor and ann: ds.append(round((e-t+ed)/ann,2))
    r=irr(cfs)
    print(f"{label:46s} capex {capex:5.1f} EBITDA1 {eb[0]:5.2f} EBITDA2 {eb[1]:5.2f} NPV15 {npv(cfs,.15):6.1f} NPV20 {npv(cfs,.20):6.1f} IRR {'neg' if r is None else round(r*100,1)}"+(f" DSCR y1-7 {ds}" if ds else ''))
# Otakikpo actual: >$60M for LPG 12 MMscf/d + 20 MW (ThisDay 19.03.2025). gas engines $1.0-1.5M/MW [Д] -> LPG share $30-40M
for gen in (1.0,1.5):
    lpg12=60-20*gen; c15=lpg12*(15/12)**0.6
    print('Otakikpo-derived LPG capex @15:',round(lpg12,1),'->',round(c15,1),' ratio to screen 18.2:',round(c15/18.2,2))
print()
run(18.2*2.4,y=4.0,price=650,label='evaluation rich add-on x2.4 (repro)')
run(18.2*1.6,y=4.0,price=650,LF=0.85,decl=0.05,label='evaluation "what it takes" (repro)')
print('-- realistic rich yield (flash 3.4-3.8, NJTD 2.73) --')
for y in (2.73,3.4):
  for cap in (34.3,45.7):
    run(cap,y=y,price=650,label=f'y{y} capex {cap} p650')
    run(cap,y=y,price=650,LF=0.85,decl=0.05,label=f'y{y} capex {cap} p650 LF.85 d5%')
print('-- favourable everything except capex = Otakikpo-derived low --')
run(34.3,y=4.0,price=650,LF=0.85,decl=0.05,label='y4 p650 LF.85 d5% capex34.3')
run(34.3,y=4.0,price=750,LF=0.85,decl=0.05,label='y4 p750 LF.85 d5% capex34.3')
print('-- year-2 stresses on evaluation best case (x1.6, y4, p650, LF.85, d5%) --')
b=dict(capex=29.1,y=4.0,price=650,LF=0.85,decl=0.05)
run(**b,yr2_LF=0.55,label='yr2 outage LF .55 (TNP/FM, Otakikpo ~60%)')
run(**b,price_from3=550,label='price $550 from yr3 (NALPGAM N900-1100 retail)')
run(**b,price_from3=450,label='price $450 from yr3')
run(**b,infl=0.05,label='fixed opex +5%/y USD, price flat')
print('anchor EBITDA per MMBtu: mining27.7',round(2.28e6/(1050*365),2),'power .12/.20',round(2.76e6/(1050*365),2),round(5.69e6/(1050*365),2))
run(**b,shrink=5.95,label='shrinkage at mining EBITDA $5.95/MMBtu')
run(**b,shrink=7.2,label='shrinkage at captive power EBITDA $7.2/MMBtu')
run(**dict(b,decl=0.15),label='decline 15%/y')
run(**dict(b,LF=0.6),label='LF 0.60 all life (Otakikpo trials)')
print('-- financing: 60% debt 15% 7y --')
run(**b,debt=0.6,label='best case, 60% debt')
run(**b,debt=0.6,yr2_LF=0.55,label='best case, 60% debt, yr2 LF .55')
run(34.3,y=3.4,price=650,LF=0.8,decl=0.08,debt=0.6,label='realistic rich y3.4 capex34.3, 60% debt')
