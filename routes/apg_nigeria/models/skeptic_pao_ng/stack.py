# Skeptic (Nigeria lens) re-run of PAO scenario A; each correction applied cumulatively.
def npv_cf(cfs,r): return sum(c/(1+r)**t for t,c in enumerate(cfs))
def irr(cfs):
    lo,hi=-0.9,1.5
    if npv_cf(cfs,0)<0: return None
    for _ in range(80):
        m=(lo+hi)/2
        if npv_cf(cfs,m)>0: lo=m
        else: hi=m
    return m
def run(sc=1,Y=0.35,olig=0.85,price=3200,capex=12.0,gas=0.25,var=250,logi=150,staff=0.6,sec=0.3,maint=0.035,ins=0.015,
        ch4=7300,util=1.0,ramp=1.0,tax=0.30,wc=0.0,demand_cap=None,life=12,build=2,dep=10,par=0.12):
    cap=capex*sc**0.7; st=staff*sc**0.4; se=sec*sc**0.3
    t_cap=ch4*sc*(1-par)*Y*28/32.08*olig
    t=t_cap*util
    if demand_cap: t=min(t,demand_cap)
    cfs=[-cap/build]*build; es=[]
    for y in range(life):
        f=ramp if y==0 else 1.0
        tt=t*f; rev=tt*price/1e6
        gasc=365000*1.05*sc*gas/1e6       # take-or-pay: gas paid on full contracted volume
        e=rev-gasc-tt*(var+logi)/1e6-(st+se+cap*(maint+ins))
        es.append(e); d=cap/dep if y<dep else 0
        cf=e-max(0,e-d)*tax
        if y==0: cf-=wc*t*price/1e6
        if y==life-1: cf+=wc*t*price/1e6
        cfs.append(cf)
    i=irr(cfs)
    return dict(t=round(t),e=round(es[1],2),capex=round(cap,1),npv15=round(npv_cf(cfs,.15),1),npv20=round(npv_cf(cfs,.20),1),irr=None if i is None else round(i*100,1))
steps=[('A as evaluated',{}),
 ('CH4 at 60F 6999 t [O]',dict(ch4=6999)),
 ('gas supply util 0.78 (VIIRS 2020-25) x plant 0.92',dict(util=0.78*0.92)),
 ('tax 34% CIT+dev levy [V]',dict(tax=0.34)),
 ('ramp y1 0.6 + WC 20% rev (90-day export repatriation) [D]',dict(ramp=0.6,wc=0.20)),
 ('price $2700 (mid of IMARC US/CN/DE) [V/O]',dict(price=2700)),
 ('price $2100 EU (shipping to Antwerp) [V]',dict(price=2100)),
]
for sc,dc in ((1,None),(5,None),(15,None),(15,12000),(5,12000)):
    cur={}
    print(f'=== scale {sc} demand_cap {dc}')
    for name,ch in steps:
        cur.update(ch); p=dict(cur); p['sc']=sc; p['demand_cap']=dc
        print(f'{name:60s}',run(**p))
# lower-bound yield rule at each scale with realism, price 2700
for sc in (1,5,15):
    for Y,o in ((0.25,0.75),(0.20,0.70)):
        print('Ylow',sc,Y,run(sc=sc,Y=Y,olig=o,ch4=6999,util=0.78*0.92,tax=0.34,ramp=0.6,wc=0.2,price=2700))
# flare gas value: what fraction of revenue is gas at flare vs pipeline
for g in (0.25,1.0,2.18,2.68): print('gas',g,'MUSD/y',round(365000*1.05*g/1e6,2))
