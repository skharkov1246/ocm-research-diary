# Флаер -> генерация -> майнинг BTC / AI-HPC (Нигерия). Все входы - допущения [Д] с диапазоном, кроме помеченных.
import json, itertools
OUT='/tmp/claude-0/-home-user-ocm-research-diary/50ead64b-47e2-5bdd-b34b-09f6667b6ca9/scratchpad/compute/mining_ai_out.json'
MW_GROSS=4.64          # корпус econ_vs_mining: 1 MMscf/д, 1000 Btu/scf, КПД 0.38
MMBTU_D=1050           # 1.03-1.1 MMBtu/Mscf (контекст C)
def npv(r,cf): return sum(c/(1+r)**t for t,c in enumerate(cf))
def irr(cf):
    lo,hi=-0.99,5.0
    f=lambda r:npv(r,cf)
    if f(lo+1e-6)*f(hi)>0: return None
    for _ in range(200):
        m=(lo+hi)/2
        if f(lo+1e-6)*f(m)<=0: hi=m
        else: lo=m
    return m
S={ # сценарии площадки
 'pess':dict(derate=0.92,parasit=0.92,uptime=0.80,gen_kw=1960,cond=1.8,infra=2.8,dev=1.0,landed=1.25,
             tor=0.05,staff=0.45,sec=0.50,ins=0.015,comm=0.03,pool=0.025,gas=1.0,erosion=0.30,halv=0.50,salv=0.30),
 'base':dict(derate=0.95,parasit=0.94,uptime=0.88,gen_kw=1700,cond=1.2,infra=2.6,dev=0.7,landed=1.20,
             tor=0.04,staff=0.35,sec=0.35,ins=0.012,comm=0.02,pool=0.02,gas=0.5,erosion=0.20,halv=0.55,salv=0.40),
 'opt': dict(derate=0.97,parasit=0.96,uptime=0.93,gen_kw=1500,cond=0.75,infra=2.4,dev=0.5,landed=1.15,
             tor=0.035,staff=0.30,sec=0.25,ins=0.010,comm=0.01,pool=0.015,gas=0.25,erosion=0.10,halv=0.65,salv=0.50),
}
ASIC={'A_13.5JTH':dict(e=13.5,usd_th={'pess':17,'base':13,'opt':9}),     # S21 XP-класс: $8.9-16.7/TH [В]
      'B_8.9JTH': dict(e=8.9, usd_th={'pess':25,'base':20,'opt':15})}    # S23 XP Hyd, цена НЕ НАЙДЕНА [Д]
SCALE={1:dict(gen=1.0,infra=1.0,staff=1.0,sec=1.0,dev=1.0),
       5:dict(gen=0.9,infra=0.9,staff=1.8,sec=1.5,dev=1.5),
       15:dict(gen=0.8,infra=0.8,staff=2.5,sec=2.0,dev=2.0)}
HALVING_MONTH=10       # старт ~06.2027, халвинг ~04.2028 [Д]
YEARS=4                # жизнь парка ASIC в модели
TAX=0.34               # CIT 30% + development levy 4% (NTA 2025) [В]
def mining(hp0,asic,sc,q=1,years=YEARS):
    p=S[sc]; a=ASIC[asic]; k=SCALE[q]
    mw_it=MW_GROSS*q*p['derate']*p['parasit']
    th=mw_it*1e6/a['e']; ph=th/1000
    capex_asic=th*a['usd_th'][sc]*p['landed']/1e6
    capex_gen=MW_GROSS*q*1000*p['gen_kw']*k['gen']/1e6
    capex_cond=p['cond']+0.2*(q-1)
    capex_infra=p['infra']*q*k['infra']
    capex_dev=p['dev']*k['dev']
    bop=capex_gen+capex_cond+capex_infra+capex_dev
    capex=bop+capex_asic
    cf=[-capex]; rows=[]
    for y in range(years):
        rev=0
        for m in range(12):
            mm=y*12+m
            hp=hp0*(1-p['erosion'])**(mm/12)*(p['halv'] if mm>=HALVING_MONTH else 1.0)
            rev+=ph*hp*30.42*p['uptime']/1e6
        gas=MMBTU_D*q*365*p['uptime']*p['gas']/1e6
        opex=gas+p['tor']*capex+p['staff']*k['staff']+p['sec']*k['sec']+p['ins']*capex+(p['comm']+p['pool'])*rev+0.03
        e=rev-opex
        dep=capex_asic/years+bop/10
        tax=max(0,TAX*(e-dep))
        c=e-tax+(p['salv']*bop if y==years-1 else 0)
        cf.append(c); rows.append(dict(y=y+1,rev=round(rev,2),opex=round(opex,2),ebitda=round(e,2)))
    acc=-capex; pb=None
    for t,c in enumerate(cf[1:],1):
        if acc+c>=0 and pb is None: pb=t-1+(-acc)/c
        acc+=c
    gwh=mw_it*8760*p['uptime']/1e3
    return dict(q=q,mw_it=round(mw_it,2),ph=round(ph,1),capex=round(capex,2),capex_asic=round(capex_asic,2),bop=round(bop,2),
                rows=rows,ebitda_y1=rows[0]['ebitda'],ebitda_avg=round(sum(r['ebitda'] for r in rows)/years,2),
                npv15=round(npv(0.15,cf),2),npv20=round(npv(0.20,cf),2),irr=(round(irr(cf)*100,1) if irr(cf) is not None else None),
                payback=(round(pb,1) if pb else None),
                rev_per_mwh_y1=round(rows[0]['rev']*1e6/(gwh*1e3),1),
                fixed_site_cost_per_mwh=round((p['tor']*bop+p['staff']*k['staff']+p['sec']*k['sec']+p['ins']*bop+0.03+MMBTU_D*q*365*p['uptime']*p['gas']/1e6)*1e6/(gwh*1e3),1),
                bop_capital_per_mwh_15=round(bop*1e6*0.199/(gwh*1e3),1))
out={'mining':{}}
for hp,asic,sc in itertools.product([27.7,39.87,50.0],ASIC,S):
    out['mining'][f'hp{hp}_{asic}_{sc}']=mining(hp,asic,sc)
for q in [5,15]:
    for hp in [27.7,39.87,50.0]:
        for sc in S:
            out['mining'][f'q{q}_hp{hp}_B_8.9JTH_{sc}']=mining(hp,'B_8.9JTH',sc,q=q)
# корпусный скрин и аудит на нём
mw=4.64; th=mw*1e6/17.5
corp=dict(capex=9.02,asic=th*18/1e6,gen=mw*700/1e3,infra=1.0)
corp['audit_literal']=round(corp['asic']*5.3+ (corp['gen']+corp['infra'])*2.4,2), round(corp['asic']*6.0+(corp['gen']+corp['infra'])*2.8,2)
corp['audit_justified']=round(corp['asic']*1.15+(corp['gen']+corp['infra'])*2.4,2), round(corp['asic']*1.25+(corp['gen']+corp['infra'])*2.8,2)
# корпусная формула при hp 27.7 (повтор econ_vs_mining)
def corp_net(hp,om=0.15,capex_mult=1.0):
    ph=th/1000; g=ph*hp*365/1e6; tot=g+g*0.7
    asic=th*18/1e6*capex_mult; gen=mw*700/1e3*capex_mult
    return round((tot-tot*om-asic-gen*0.2)/2,2), round(tot/2*(1-om),2)
corp['net_ebitda_hp50']=corp_net(50); corp['net_ebitda_hp27.7']=corp_net(27.7)
corp['net_ebitda_hp50_om_x1.7_2.4']=[corp_net(50,0.15*1.7),corp_net(50,0.15*2.4)]
corp['net_ebitda_hp27.7_om_x2.4_capex_x2.4']=corp_net(27.7,0.36,2.4)
corp['net_ebitda_hp50_om_x2.4_capex_x2.4']=corp_net(50,0.36,2.4)
out['corpus_screen']=corp
# хостинг: продаём э/э чужим ASIC
def hosting(price_kwh,sc):
    r=mining(39.87,'B_8.9JTH',sc); p=S[sc]
    bop=r['bop']; gwh=r['mw_it']*8760*p['uptime']/1e3
    rev=gwh*price_kwh
    opex=MMBTU_D*365*p['uptime']*p['gas']/1e6+p['tor']*bop+p['staff']+p['sec']+p['ins']*bop+0.03
    e=rev-opex; cf=[-bop]+[e*(1-0.0)-max(0,TAX*(e-bop/10))]*10
    return dict(price=price_kwh,bop=bop,rev=round(rev,2),opex=round(opex,2),ebitda=round(e,2),npv15=round(npv(0.15,cf),2),irr=(round(irr(cf)*100,1) if irr(cf) else None))
out['hosting']={f'{pr}_{sc}':hosting(pr,sc) for pr in [0.05,0.07,0.10] for sc in S}
# AI/HPC на 1 MMscf/д
AI={'pess':dict(pue=1.45,price=2.5,util=0.60,gpu_usd=70e3,kw_gpu=2.0,fac_mw=6.0,red=1.25,gen_kw=1960,fiber_km=40,usd_km=30e3,
                staff=4.0,bw=1.5,diesel=1.2,life=4,salv=0.10,sec=1.0),
    'base':dict(pue=1.35,price=3.0,util=0.70,gpu_usd=62e3,kw_gpu=2.0,fac_mw=4.5,red=1.25,gen_kw=1700,fiber_km=20,usd_km=20e3,
                staff=3.0,bw=1.0,diesel=0.8,life=5,salv=0.15,sec=0.7),
    'opt': dict(pue=1.25,price=4.0,util=0.85,gpu_usd=55e3,kw_gpu=2.0,fac_mw=3.0,red=1.15,gen_kw=1500,fiber_km=10,usd_km=15e3,
                staff=2.0,bw=0.5,diesel=0.5,life=5,salv=0.20,sec=0.5)}
def ai(sc):
    p=AI[sc]
    mw_it=MW_GROSS*0.95*0.96/p['pue']
    n=mw_it*1000/p['kw_gpu']
    gpu=n*p['gpu_usd']*1.1/1e6
    fac=mw_it*p['fac_mw']
    gen=MW_GROSS*1000*p['gen_kw']*p['red']/1e6
    other=p['fiber_km']*p['usd_km']/1e6+1.5+1.0  # fiber + diesel/BESS backup + dev
    capex=gpu+fac+gen+other
    rev=n*8760*p['util']*p['price']/1e6
    opex=p['staff']+p['bw']+p['diesel']+p['sec']+0.015*capex+0.04*(fac+gen)+0.3
    e=rev-opex
    cf=[-capex]+[e-max(0,TAX*(e-capex/p['life']))]*p['life']; cf[-1]+=p['salv']*gpu+0.4*(fac+gen)
    return dict(mw_it=round(mw_it,2),gpus=int(n),capex=round(capex,1),gpu_capex=round(gpu,1),rev=round(rev,1),opex=round(opex,1),ebitda=round(e,1),
                energy_share_of_cost=round((0.04*gen+0.3)/(opex+capex/p['life']),3),
                npv15=round(npv(0.15,cf),1),npv20=round(npv(0.20,cf),1),irr=(round(irr(cf)*100,1) if irr(cf) is not None else None),payback=round(capex/e,1) if e>0 else None)
out['ai']={sc:ai(sc) for sc in AI}
json.dump(out,open(OUT,'w'),indent=1,ensure_ascii=False)
for k,v in out['mining'].items():
    print(k,{kk:v[kk] for kk in ['mw_it','ph','capex','capex_asic','bop','ebitda_y1','ebitda_avg','npv15','npv20','irr','payback','rev_per_mwh_y1','fixed_site_cost_per_mwh','bop_capital_per_mwh_15']})
print(json.dumps(out['corpus_screen'],ensure_ascii=False)); print(json.dumps(out['hosting'])); print(json.dumps(out['ai'],ensure_ascii=False))
