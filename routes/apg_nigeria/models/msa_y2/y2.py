# Скептик-2: второй год эксплуатации, деньги, инфляция, простой. MSA-скид.
def npv(r,c): return sum(x/(1+r)**i for i,x in enumerate(c))
def irr(c):
    f=lambda r: npv(r,c); lo,hi=-0.9,2.0
    if f(lo)*f(hi)>0: return None
    for _ in range(200):
        m=(lo+hi)/2
        if f(lo)*f(m)<=0: hi=m
        else: lo=m
    return m
def run(name,Tnom,cap,P0,v,staff,om,ins,logi,ramp,S_path,uS=0.95,infl=0.03,pdrift=0.0,
        build=2,wc_days=90,life=12,tax=0.34,sales_cap=None,outage=None,verbose=False):
    cf=[-cap/build]*build   # стройка
    wc_prev=0; rows=[]; cum=0; trough=0
    for y in range(1,life+1):
        on=ramp[y-1] if y-1<len(ramp) else ramp[-1]
        if outage and y in outage: on*=outage[y]
        T=Tnom*on
        if sales_cap: T=min(T,sales_cap)
        P=P0*(1+pdrift)**(y-1); k=(1+infl)**(y-1)
        PS=S_path[y-1] if y-1<len(S_path) else S_path[-1]
        rev=T*P/1e6
        opex=T*(0.333/uS)*PS/1e6 + T*v*k/1e6 + staff*k + cap*om*k + cap*ins + T*logi*k/1e6 + T*0.167/uS*52*0.5/1e6
        e=rev-opex
        wc=rev*wc_days/365; dwc=wc-wc_prev; wc_prev=wc
        taxv=max(0,(e-(cap/10 if y<=10 else 0))*tax)
        c=e-taxv-dwc
        if y==life: c+=wc
        cf.append(c); rows.append((y,round(T),round(rev,2),round(opex,2),round(e,2),round(c,2)))
    cum=0
    for x in cf:
        cum+=x; trough=min(trough,cum)
    ir=irr(cf)
    print(f"\n== {name}: NPV15 {npv(0.15,cf):.2f} NPV20 {npv(0.20,cf):.2f} NPV25 {npv(0.25,cf):.2f} IRR {ir*100 if ir is not None else float('nan'):.1f}% пик финансирования {trough:.2f}")
    for r in rows[:3]: print("  год",r[0],"т",r[1],"выручка",r[2],"опекс",r[3],"EBITDA",r[4],"CF",r[5])
    return cf
ramp_foak=[0.40,0.60,0.75,0.85]       # FOAK [Д]
ramp_eval=[0.45,0.90]                  # как в оценке: 50% в 1-й, далее 0.90
S_crisis=[925]*12
S_revert=[925,500,300]                 # [Д] откат после Ормуза
# 1.1 кт без двойного счёта, лучший край (skeptic A)
run("A оценка-эквивалент (без двойн.счёта, 1.1кт)",1093,10.8,2200,230,0.55,0.035,0.01,40,ramp_eval,S_crisis,build=1,wc_days=0,infl=0)
run("A + стройка 2г, FOAK-рамп, инфл 3%, оборотка 90д, сера откат",1093,10.8,2200,230,0.55,0.035,0.01,40,ramp_foak,S_revert)
run("A + то же + цена 2800 (70%-база Comtrade)",1093,10.8,2800,230,0.55,0.035,0.01,40,ramp_foak,S_revert)
run("A + цена 2800, сера откат, капекс 4.5 (без аудита)",1093,4.5,2800,230,0.55,0.035,0.01,40,ramp_foak,S_revert)
run("A + простой 4 мес во 2-й год (коррозия/перекрытие трубы)",1093,10.8,2200,230,0.55,0.035,0.01,40,ramp_foak,S_revert,outage={2:8/12})
run("A + цена дрейф -3%/г (Китай)",1093,10.8,2200,230,0.55,0.035,0.01,40,ramp_foak,S_revert,pdrift=-0.03)
# 8.7 кт
run("8.7кт капекс x1.3 (оценка, лучший) как в оценке",8747,14.0,2200,230*1.7,1.2,0.035,0.01,40,ramp_eval,S_crisis,build=1,wc_days=0,infl=0)
run("8.7кт капекс x1.3, FOAK, сбыт ≤7.7кт, цена -$200 за объём",8747,14.0,2000,230,1.2,0.035,0.01,40,ramp_foak,S_revert,sales_cap=7700)
run("8.7кт капекс x1.3, сбыт ≤3.7кт",8747,14.0,2200,230,1.2,0.035,0.01,40,ramp_foak,S_revert,sales_cap=3700)
for cap in (37.6,47.0,60.0):
    run(f"8.7кт капекс {cap} (0.6 / бенчмарк BASF), FOAK, сбыт≤7.7кт, $2000",8747,cap,2000,230,1.2,0.035,0.01,40,ramp_foak,S_revert,sales_cap=7700)
run("8.7кт капекс 37.6, цена 2800, сера откат, сбыт 7.7",8747,37.6,2800,230,1.2,0.035,0.01,40,ramp_foak,S_revert,sales_cap=7700)
# BASF бенчмарк: 'higher double-digit million EUR' за прирост до 50 кт
for eur in (50,99):
    for base in (20,):
        print("BASF ~",eur,"MEUR на",base,"кт -> 1.1кт:",round(eur*(1.1/base)**0.6*1.08,1),"M$; 8.7кт:",round(eur*(8.7/base)**0.6*1.08,1),"M$")
