# 10-летняя модель майнинга с заменой парка ASIC каждые 4 года; поиск безубыточного hashprice
import sys, json
sys.path.insert(0,'/tmp/claude-0/-home-user-ocm-research-diary/50ead64b-47e2-5bdd-b34b-09f6667b6ca9/scratchpad/compute')
from mining_ai import S, ASIC, SCALE, MW_GROSS, MMBTU_D, TAX, npv, irr
LEAN=dict(derate=0.95,parasit=0.94,uptime=0.88,gen_kw=700,cond=0.5,infra=1.0,dev=0.3,landed=1.10,
          tor=0.035,staff=0.20,sec=0.15,ins=0.01,comm=0.01,pool=0.02,gas=0.25,erosion=0.20,halv=0.55,salv=0.40)
S2=dict(S); S2['lean_noaudit']=LEAN
HALVINGS=[10,10+48]   # ~04.2028, ~2032
def run(hp0,asic,sc,q=1,years=10,refresh=4,refresh_cost=0.8,eff_gain=0.75):
    p=S2[sc]; a=ASIC[asic]; k=SCALE[q]
    usd_th=a['usd_th'].get(sc,a['usd_th']['opt'] if sc=='lean_noaudit' else None)
    mw_it=MW_GROSS*q*p['derate']*p['parasit']
    bop=MW_GROSS*q*1000*p['gen_kw']*k['gen']/1e6+p['cond']+0.2*(q-1)+p['infra']*q*k['infra']+p['dev']*k['dev']
    e=a['e']; price=usd_th
    th=mw_it*1e6/e; asic_capex=th*price*p['landed']/1e6
    cf=[-(bop+asic_capex)]; eb=[]
    for y in range(years):
        if y>0 and y%refresh==0:
            e*=eff_gain; price*=refresh_cost; th=mw_it*1e6/e; asic_capex=th*price*p['landed']/1e6
            cf[-1]-=asic_capex
        rev=0
        for m in range(12):
            mm=y*12+m; f=1.0
            for h in HALVINGS:
                if mm>=h: f*=p['halv']
            hp=hp0*(1-p['erosion'])**(mm/12)*f
            rev+=th/1000*hp*30.42*p['uptime']/1e6
        capex_tot=bop+asic_capex
        opex=MMBTU_D*q*365*p['uptime']*p['gas']/1e6+p['tor']*capex_tot+p['staff']*k['staff']+p['sec']*k['sec']+p['ins']*capex_tot+(p['comm']+p['pool'])*rev+0.03
        ebitda=rev-opex
        dep=asic_capex/refresh+bop/10
        cf.append(ebitda-max(0,TAX*(ebitda-dep))); eb.append(round(ebitda,2))
    return dict(npv15=round(npv(0.15,cf),2),npv20=round(npv(0.20,cf),2),irr=(round(irr(cf)*100,1) if irr(cf) is not None else None),
                ebitda=eb,capex0=round(-cf[0],2),bop=round(bop,2))
def be(asic,sc,q=1,r=0.15):
    lo,hi=5,400
    for _ in range(60):
        m=(lo+hi)/2
        v=run(m,asic,sc,q)['npv15' if r==0.15 else 'npv20']
        if v<0: lo=m
        else: hi=m
    return round(m,1)
out={}
for sc in ['pess','base','opt','lean_noaudit']:
    for asic in ASIC:
        for hp in [27.7,39.87,50.0]:
            out[f'{sc}_{asic}_hp{hp}']=run(hp,asic,sc)
        out[f'{sc}_{asic}_BE_hp_npv15']=be(asic,sc)
        out[f'{sc}_{asic}_BE_hp_npv20']=be(asic,sc,r=0.20)
for q in [5,15]:
    for sc in ['pess','base','opt']:
        out[f'q{q}_{sc}_A_BE_hp_npv15']=be('A_13.5JTH',sc,q)
        out[f'q{q}_{sc}_A_hp50']=run(50,'A_13.5JTH',sc,q)
        out[f'q{q}_{sc}_A_hp27.7']=run(27.7,'A_13.5JTH',sc,q)
# без халвинга-эрозии (BTC растёт так, что hashprice в $ постоянен) - верхний опцион
for sc in ['base','opt']:
    p=S2[sc]; old=(p['erosion'],p['halv']); p['erosion']=0.0; p['halv']=1.0
    out[f'{sc}_A_flatHP_BE']=be('A_13.5JTH',sc); out[f'{sc}_A_flatHP_hp39.87']=run(39.87,'A_13.5JTH',sc)
    p['erosion'],p['halv']=old
json.dump(out,open('/tmp/claude-0/-home-user-ocm-research-diary/50ead64b-47e2-5bdd-b34b-09f6667b6ca9/scratchpad/compute/mining10_out.json','w'),indent=1)
for k,v in out.items(): print(k,v)
