import importlib.util,sys,io,contextlib
with contextlib.redirect_stdout(io.StringIO()):
    spec=importlib.util.spec_from_file_location("m","msa_ng.py"); m=importlib.util.module_from_spec(spec); spec.loader.exec_module(m)
def itemized(name,T_nom,on,u,P,PS,v,staff,capex,om,ins,logi,gas,wc_months=0,build_yrs=0):
    T=T_nom*on; rev=T*P/1e6
    opex=T*(0.333/u)*PS/1e6+T*v/1e6+staff+capex*om+capex*ins+T*logi/1e6+T*0.167/u*52*gas/1e6
    e=rev-opex
    cf=m.cash(capex,e)
    if wc_months: # оборотный капитал на 1-й год, возврат в конце
        wc=rev*wc_months/12; cf[1]-=wc; cf[-1]+=wc
    if build_yrs: cf=[-capex/2,-capex/2]+cf[1:]
    ir=m.irr(cf)
    print(f"{name:45s} T={T:.0f} rev={rev:.2f} opex={opex:.2f} ($/t {opex*1e6/T:.0f}) EBITDA={e:.2f} PB={capex/e if e>0 else float('nan'):.1f} NPV15={m.npv(0.15,cf):.2f} IRR={ir*100 if ir is not None else float('nan'):.1f}")
# без двойного счёта: персонал/ТОиР/страхование статьями, без глобального x1.7-2.4 на фикс
itemized("A без двойного счёта, лучший",1093,0.90,0.95,2200,925,230,0.55,10.8,0.035,0.01,40,0.25)
itemized("B без двойного счёта, худший",1093,0.85,0.90,1700,965,230*1.7,1.0,12.6,0.05,0.02,80,1.0)
itemized("C лучший + цена 2800",1093,0.90,0.95,2800,925,230,0.55,10.8,0.035,0.01,40,0.25)
itemized("D лучший + цена 2800 + сера 250",1093,0.90,0.95,2800,250,230,0.55,10.8,0.035,0.01,40,0.25)
itemized("E D + цена 3000 (если Comtrade за т 70%-раствора)",1093,0.90,0.95,3000,250,230,0.55,10.8,0.035,0.01,40,0.25)
itemized("F A + оборотка 3 мес + стройка 2 г",1093,0.90,0.95,2200,925,230,0.55,10.8,0.035,0.01,40,0.25,3,2)
itemized("G A без аудита капекса (4.5)",1093,0.90,0.95,2200,925,230,0.55,4.5,0.035,0.01,40,0.25)
itemized("H 8.7кт правило 0.6 лучший, без дв.счёта",8747,0.90,0.95,2200,925,230,1.2,37.6,0.035,0.01,40,0.25)
itemized("I 8.7кт 0.6, цена 2800 сера 250",8747,0.90,0.95,2800,250,230,1.2,37.6,0.035,0.01,40,0.25)
itemized("J 8.7кт 0.6 x2.8, худший",8747,0.85,0.90,1700,965,230*1.7,1.8,43.9,0.05,0.02,80,1.0)
# порог капекса для IRR 15% в лучшем 1.1кт
for cap in (2,3,4,5,6):
    itemized(f"K капекс {cap} при A",1093,0.90,0.95,2200,925,230,0.55,cap,0.035,0.01,40,0.25)
# генерация с симметричным аудитом
for lo in (True,False):
    E=4.64*8760*0.85/1e3  # GWh
    tar=0.15 if lo else 0.25; cap=(4.3*2.8 if lo else 6.3*2.4)
    eb=E*tar-E*(0.03 if lo else 0.02)-cap*(0.05 if lo else 0.035)-(1.0 if lo else 0.55)-(0.38 if lo else 0.1)
    print("генерация с аудитом", "худший" if lo else "лучший", "capex",round(cap,1),"EBITDA",round(eb,2), "EBITDA/capex", round(eb/cap*100))
