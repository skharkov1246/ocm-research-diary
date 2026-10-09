def npv(r,c): return sum(x/(1+r)**i for i,x in enumerate(c))
def irr(c):
    lo,hi=-0.99,2.0
    if npv(lo,c)*npv(hi,c)>0: return None
    for _ in range(200):
        m=(lo+hi)/2
        if npv(lo,c)*npv(m,c)<=0: hi=m
        else: lo=m
    return m
def run(name,Tn,on,u,P,PS,v,staff,capex,om,ins,logi,gas=0.5,roy=0.0,wc_m=0,build=0,comm=0.0):
    T=Tn*on; rev=T*P/1e6
    op=T*(0.333/u)*PS/1e6+T*v/1e6+staff+capex*(om+ins)+T*logi/1e6+T*0.167*52*gas/1e6+rev*roy
    op*= (1+comm)
    e=rev-op; dep=capex/10
    cf=[-capex] if not build else [-capex/2,-capex/2]
    for y in range(1,13):
        ee=e*(0.5 if y==1 else 1); cf.append(ee-max(0,(ee-(dep if y<=10 else 0))*0.34))
    if wc_m:
        wc=rev*wc_m/12; i0=len(cf)-12; cf[i0]-=wc; cf[-1]+=wc
    r=irr(cf)
    print(f"{name:52s} EBITDA={e:6.2f} $/t={op*1e6/T:5.0f} PB={(capex/e if e>0 else float('nan')):5.1f} NPV15={npv(0.15,cf):7.2f} NPV20={npv(0.20,cf):7.2f} IRR={(r*100 if r is not None else float('nan')):6.1f}")
    return e
# capex with corrosion (wet part 50%, x1.3-1.5) and location 1.1-1.3 that evaluator left out
lo=4.5*2.4*(0.5+0.5*1.3)*1.1; hi=4.5*2.8*(0.5+0.5*1.5)*1.3
print("capex с коррозией и локацией:",round(lo,1),round(hi,1))
print("capex если 4.5 = котировка модулей FOB x5.3-6:",round(4.5*5.3,1),round(4.5*6,1))
# A: no double count: v x1 (sulfur repriced separately) .. x1.7, staff bottom-up
run("A лучший, без дв.счёта, капекс 10.8, P2200",1093,.9,.95,2200,925,230,0.55,10.8,.035,.01,40,0.25)
run("A2 лучший, P2800 (если Comtrade=70%)",1093,.9,.95,2800,925,230,0.55,10.8,.035,.01,40,0.25)
run("A3 лучший, P2800, сера 250",1093,.9,.95,2800,250,230,0.55,10.8,.035,.01,40,0.25)
run("A4 A2 + капекс с коррозией/локацией 13.7",1093,.9,.95,2800,925,230,0.55,lo,.035,.01,40,0.25)
run("B реалистичный: P2200, v x1.7, staff .8, капекс 13.7, роялти 4%, WC 4 мес, стройка 2г, община 2%",1093,.875,.925,2200,945,391,0.8,lo,.04,.015,60,0.6,roy=0.04,wc_m=4,build=2,comm=0.02)
run("C худший: P1700, капекс 20.5, роялти 5%, WC 4м, стройка 2г",1093,.85,.9,1700,965,552,1.0,hi,.05,.02,80,1.0,roy=0.05,wc_m=4,build=2,comm=0.03)
# breakeven capex for NPV15=0 in best case with P 2200/2800
for P in (2200,2800,3000):
    for cap in [x/2 for x in range(2,30)]:
        T=1093*.9; e=T*P/1e6-(T*.3505*925/1e6+T*230/1e6+0.55+cap*.045+T*40/1e6)
        cf=[-cap]+[ (e*(0.5 if y==1 else 1)) - max(0,((e*(0.5 if y==1 else 1))-(cap/10 if y<=10 else 0))*.34) for y in range(1,13)]
        if npv(0.15,cf)<0: print("P",P,"капекс-порог NPV15=0 ≈ $",cap-0.5,"млн"); break
# 8.7 kt
run("8.7кт 0.6-правило x2.4, P2800, S250, v x1",8747,.9,.95,2800,250,230,1.2,37.6,.035,.01,40,0.25)
run("8.7кт 0.6 x2.4 + корр/лок (47.7), P2200, S925",8747,.9,.95,2200,925,230,1.2,4.5*8**0.6*2.4*1.15*1.1,.035,.01,40,0.25)
# what 8.7kt does to price: entry 7.9kt into ~ 40-60kt seaborne trade
# generation with audit
E=4.64*8760*0.85/1e3
for tag,tar,cap,om_kwh,fix in (("лучший",0.25,4.3*2.4,0.02,0.55),("худший",0.15,6.3*2.8,0.03,1.0)):
    eb=E*tar-E*om_kwh-cap*0.04-fix
    print("генерация с аудитом",tag,"капекс",round(cap,1),"EBITDA",round(eb,2),"PB",round(cap/eb,1))
