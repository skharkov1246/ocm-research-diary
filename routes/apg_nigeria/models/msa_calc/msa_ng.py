# MSA-скид Нигерия: скрин корпуса -> поправки аудита NaCN -> Нигерия
import json
CH4_T, S_T = 0.167, 0.333
MMBTU_PER_T_CH4 = 52.0   # HHV ~55.5 MJ/kg -> 52.6 MMBtu/t
def irr(cfs):
    lo,hi=-0.99,1.0
    f=lambda r: sum(c/(1+r)**i for i,c in enumerate(cfs))
    if f(lo)*f(hi)>0: return None
    for _ in range(200):
        m=(lo+hi)/2
        if f(lo)*f(m)<=0: hi=m
        else: lo=m
    return m
def npv(r,cfs): return sum(c/(1+r)**i for i,c in enumerate(cfs))
def cash(capex, ebitda, life=12, tax=0.34, ramp=0.5):
    dep=capex/10; cfs=[-capex]
    for y in range(1,life+1):
        e=ebitda*(ramp if y==1 else 1)
        t=max(0,(e-(dep if y<=10 else 0))*tax)
        cfs.append(e-t)
    return cfs
def case(name, T_nom, onstream, s_util, price, s_price, var_ns, fixed_musd, capex, om_frac, ins_frac, inland, gas_mmbtu):
    T=T_nom*onstream
    rev=T*price/1e6
    sulfur=T*(S_T/s_util)*s_price/1e6
    varns=T*var_ns/1e6
    gas=T*CH4_T/s_util*MMBTU_PER_T_CH4*gas_mmbtu/1e6
    om=capex*om_frac; ins=capex*ins_frac
    logi=T*inland/1e6
    opex=sulfur+varns+gas+fixed_musd+om+ins+logi
    e=rev-opex
    cf=cash(capex,e)
    r=dict(name=name,t_y=round(T),rev=round(rev,2),sulfur=round(sulfur,3),var_ns=round(varns,3),gas=round(gas,4),
           fixed=fixed_musd,om=round(om,3),ins=round(ins,3),logi=round(logi,3),opex=round(opex,2),
           cash_per_t=round(opex*1e6/T),ebitda=round(e,2),capex=capex,
           payback=(round(capex/e,1) if e>0 else None),
           npv15=round(npv(0.15,cf),2),npv20=round(npv(0.20,cf),2),irr=(round(irr(cf)*100,1) if irr(cf) is not None else None))
    return r
R=[]
# 0. Скрин корпуса как есть (msa_cost_model): 1093 т, $2800, сера $150, вар.без серы 230, фикс 400$/т -> 0.437M, капекс 4.5, газ $30/т CH4 ~ $0.58/MMBtu
R.append(case("0 скрин корпуса",1093,1.0,1.0,2800,150,230,0.437,4.5,0,0,0,0.58))
# 1. Только рынок 2026: цена netback и сера, без аудита
R.append(case("1 скрин + цена/сера 2026",1093,1.0,1.0,2200,925,230,0.437,4.5,0,0,0,0.58))
# 2A. Аудит мультипликатором (лучший край): денежная x1.7 на не-серную часть, капекс x2.4, выход 0.9x, сера $925 (сен-2026+фрахт), цена $2200
R.append(case("2A аудит-лучший",1093,0.90,0.95,2200,925,230*1.7,0.437*1.7,4.5*2.4,0.035,0.01,40,0.25))
# 2B. Аудит худший край: x2.4/x2.8, выход 0.85/0.9, сера $965, цена $1700
R.append(case("2B аудит-худший",1093,0.85,0.90,1700,965,230*2.4,0.437*2.4,4.5*2.8,0.05,0.02,80,1.0))
# 3. снизу вверх Нигерия: персонал 16 местн.+экспат+охрана 0.55-1.0M; вар.не-сера x1.7-2.4
R.append(case("3A снизу-вверх лучший",1093,0.90,0.95,2200,925,230*1.7,0.55,4.5*2.4,0.035,0.01,40,0.25))
R.append(case("3B снизу-вверх худший",1093,0.85,0.90,1700,965,230*2.4,1.00,4.5*2.8,0.05,0.02,80,1.0))
# 4. Если сера откатится к $250 (2024-25 уровень+фрахт)
R.append(case("4A аудит-лучший, сера $250",1093,0.90,0.95,2200,250,230*1.7,0.437*1.7,4.5*2.4,0.035,0.01,40,0.25))
R.append(case("4B аудит-худший, сера $250",1093,0.85,0.90,1700,250,230*2.4,0.437*2.4,4.5*2.8,0.05,0.02,80,1.0))
# 4C. скрин-цена $2800 + аудит лучший + сера 925 (если найдётся премиальный покупатель)
R.append(case("4C аудит-лучший, цена $2800",1093,0.90,0.95,2800,925,230*1.7,0.437*1.7,4.5*2.4,0.035,0.01,40,0.25))
# 5. Крупный вариант 8.7 кт (20% конверсии): капекс скрина 4.5*1.3=5.85 -> x2.4-2.8; или правило 0.6: 4.5*8^0.6=15.7 -> x2.4-2.8
for lab,cap_lo,cap_hi in (("5 8.7кт капекс корпуса x1.3",5.85*2.4,5.85*2.8),("6 8.7кт капекс по правилу 0.6",4.5*8**0.6*2.4,4.5*8**0.6*2.8)):
    R.append(case(lab+" лучший",8747,0.90,0.95,2200,925,230*1.7,1.2,round(cap_lo,1),0.035,0.01,40,0.25))
    R.append(case(lab+" худший",8747,0.85,0.90,1700,965,230*2.4,1.8,round(cap_hi,1),0.05,0.02,80,1.0))
for r in R: print(json.dumps(r,ensure_ascii=False))
# газ на скид
for q in (1,5,15):
    print("MMscf/d",q,"доля газа на 1.1кт MSA:",round(182/(7300*q)*100,2),"%; полная конверсия даёт кт MSA:",round(7300*q/0.167/1000,1),"; S кт:",round(7300*q/0.167*0.333/1000,1))
# майнинг на 1 MMscf/d по формуле корпуса
for hp in (27.7,50):
    ph=265.1; g1=ph*hp*365/1e6; tot=g1+g1*0.7; om=tot*0.15; asic=265100*18/1e6; gen=4.64*700/1e3
    net=(tot-om-asic-gen*0.2)/2
    print("hashprice",hp,"gross1",round(g1,2),"EBITDA1",round(g1*0.85,2),"net/yr(корпус)",round(net,2),"capex",round(asic+gen+1,2))
# переработка факела на 1 MMscf/d: штраф 3.5$/Mscf
print("штраф за 1 MMscf/d $M/г:",round(3.5*1000*365/1e6,2))
