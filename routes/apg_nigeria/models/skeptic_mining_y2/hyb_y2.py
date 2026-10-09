# Скептик (год 2, деньги): подвариант «ASIC на избытке captive» с профилем нагрузки, падением дебита, стартом после халвинга
import json
def npv(r,cf): return sum(c/(1+r)**t for t,c in enumerate(cf))
MWN=4.64*0.96   # нетто captive, captive.py
AV=0.92         # готовность станции (captive base) - один раз
ASIC_UP=0.97    # собственная готовность ASIC
LF=0.72
PROFILES={ # доля нетто-мощности, которую берёт клиент, по часам суток; среднее = 0.72
 'flat':[0.72]*24,
 'two_shift':[0.93]*16+[0.30]*8,
 'one_shift_plus_base':[0.98]*12+[0.46]*12,
}
def run(hp0,prof,S=1.15,decl=0.0,halv_month=10,e=13.5,usd_th=13,landed=1.2,gas=0.5,vom=18,erosion=0.10,halv=0.60,pool=0.04,tor=0.04,years=4):
    th=S*1e6/e; capex=th*usd_th*landed/1e6+0.3
    cf=[-capex]; rows=[]
    for y in range(years):
        g=(1-decl)**y   # газ: доступная нетто-мощность падает, клиент в приоритете
        mwh=sum(min(S,max(0.0,MWN*g-MWN*d)) for d in PROFILES[prof])*365*AV*ASIC_UP
        cfac=mwh/(S*8760)
        ph=th/1000*cfac
        rev=sum(ph*hp0*(1-erosion)**((y*12+m)/12)*(halv if y*12+m>=halv_month else 1)*30.42 for m in range(12))/1e6
        opex=mwh*(gas*3.412/0.38+vom)/1e6+tor*capex+pool*rev
        eb=rev-opex; cf.append(eb-max(0,0.34*(eb-capex/4)))
        rows.append(dict(y=y+1,cf=round(cfac,2),rev=round(rev,2),ebitda=round(eb,2)))
    return dict(capex=round(capex,2),rows=rows,npv15=round(npv(.15,cf),2),npv20=round(npv(.20,cf),2))
out={}
for hp in [27.7,39.87,50.0]:
    for prof in PROFILES:
        for decl in [0.0,0.05,0.10]:
            for hm in [10,0]:
                out[f'hp{hp}|{prof}|decl{decl}|halv_m{hm}']=run(hp,prof,decl=decl,halv_month=hm)
json.dump(out,open('/tmp/claude-0/-home-user-ocm-research-diary/50ead64b-47e2-5bdd-b34b-09f6667b6ca9/scratchpad/skeptic_mining_y2/hyb_y2_out.json','w'),indent=1)
for k,v in out.items():
    if ('decl0.0' in k or 'decl0.05' in k): print(k,v['capex'],[(r['cf'],r['ebitda']) for r in v['rows']],v['npv15'],v['npv20'])
