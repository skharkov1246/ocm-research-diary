# Skeptic recheck of PAO-on-flare (Nigeria). All inputs tagged in comments.
import math
# --- mass basis
rho60=16.043/(22.414*288.71/273.15)          # kg/m3 CH4 at 60F,1atm [О]
ch4_t=1e6*0.0283168*rho60/1000*365            # t/y per 1 MMscf/d [О]
print('CH4 t/y per MMscf/d at 60F:',round(ch4_t),'vs corpus 7300 ->',round(ch4_t/7300,3))
mwth=1050*1.05506/86400*1000                  # MW thermal of 1 MMscf/d @1.05 MMBtu/Mscf
print('gas MWth',round(mwth,2),' 12% -> MWe@36%',round(0.12*mwth*0.36,2),' share of gas for 1-2 MWe:',round(1/0.36/mwth,2),round(2/0.36/mwth,2))
X=0.35/0.80; print('per-pass conversion X at Y=.35,S=.8:',round(X,3),' unreacted CH4 MWth',round((1-X)*mwth,2))
# hexene oligomer H2 demand
for n,name in [(3,'C18'),(5,'C30')]:
    mw=84.16*n; print(name,'H2 wt%',round(2.016/mw*100,2))

def npv_cf(cfs,r):
    return sum(c/(1+r)**t for t,c in enumerate(cfs))
def project(Y=0.35,olig=0.85,price=3200,capex=12.0,gas=0.25,var=250,staff=0.6,sec=0.3,logi=150,maint=0.035,ins=0.015,
            ch4=7300,avail=1.0,ramp=1.0,decl=0.0,tax=0.30,wc=0.0,owner=0.0,parasitic=0.12,life=12,build=2,dep=10,r=0.15):
    t0=ch4*(1-parasitic)*Y*28/32.08*olig*avail
    cap=capex*(1+owner)
    cfs=[-cap/build]*build
    gasc=365000*1.05*gas/1e6
    es=[]
    for y in range(life):
        f=(1-decl)**y*(ramp if y==0 else 1.0)
        t=t0*f; rev=t*price/1e6
        e=rev-gasc*f-t*(var+logi)/1e6-(staff+sec+capex*(maint+ins))
        es.append(e)
        d=cap/dep if y<dep else 0
        tx=max(0,e-d)*tax
        cf=e-tx
        if y==0: cf-=wc*rev
        if y==life-1: cf+=wc*t0*price/1e6
        cfs.append(cf)
    return dict(t=round(t0),e1=round(es[1],2),npv15=round(npv_cf(cfs,0.15),1),npv20=round(npv_cf(cfs,0.20),1))
base=dict(Y=0.35,olig=0.85,price=3200,capex=12.0,gas=0.25,var=250,staff=0.6,sec=0.3,logi=150,maint=0.035)
print('A reproduce (no ramp etc.)',project(**base))
steps=[('mass 60F',dict(ch4=ch4_t)),('avail .90',dict(avail=0.90)),('ramp y1 .6',dict(ramp=0.6)),
       ('tax 34% (CIT+levy)',dict(tax=0.34)),('WC 20% rev + owner/startup 10%',dict(wc=0.20,owner=0.10)),
       ('staff 1.2 (5 units, 3 expats)',dict(staff=1.2)),('APG decline 10%/y',dict(decl=0.10)),
       ('olig 0.75',dict(olig=0.75)),('price EU $2100 decene-grade',dict(price=2100)),('capex 22 bottom-up',dict(capex=22.0))]
cur=dict(base)
for name,ch in steps:
    cur.update(ch); print(f'{name:35s}',project(**cur))
# best case kept: decene quality, US price, but realism on plant
best=dict(base); best.update(ch4=ch4_t,avail=0.9,ramp=0.6,tax=0.34,wc=0.2,owner=0.1)
print('A + only mechanical realism (no price/capex/olig/decl/staff change)',project(**best))
# scale with all corrections except price at 5 and 15
for price in (3200,2600,2100):
  for sc in (5,15):
    p=dict(base); p.update(ch4=ch4_t*sc,avail=0.9,ramp=0.6,tax=0.34,wc=0.2,owner=0.1,olig=0.75,decl=0.10,price=price,
       capex=22*sc**0.7,staff=1.2*sc**0.4,sec=0.3*sc**0.3)
    # gas cost scales
    p['gas']=0.25
    r=project(**p); 
    # fix gas scaling: gas cost uses 1 MMscf base; add rest
    extra=365000*1.05*0.25/1e6*(sc-1)
    print('scale',sc,'price',price,'capex',round(22*sc**0.7,1),r,' (gas undercount ~',round(extra,2),'M/y)')
