# Майнинг как гибкая нагрузка на ИЗБЫТКЕ captive-электростанции (BOP оплачен captive-выручкой). Маржинальная экономика.
import json
from mining_ai import npv, irr
def hybrid(hp0, sc):
    P={'pess':dict(spare_mw=0.8,e=13.5,usd_th=17,landed=1.25,gas=1.0,vom=25,erosion=0.20,halv=0.50,uptime=0.80,pool=0.045,tor_asic=0.05),
       'base':dict(spare_mw=1.15,e=13.5,usd_th=13,landed=1.20,gas=0.5,vom=18,erosion=0.10,halv=0.60,uptime=0.88,pool=0.04,tor_asic=0.04),
       'opt': dict(spare_mw=1.4,e=13.5,usd_th=9,landed=1.15,gas=0.25,vom=12,erosion=0.0,halv=0.70,uptime=0.93,pool=0.025,tor_asic=0.035)}[sc]
    th=P['spare_mw']*1e6/P['e']; capex=th*P['usd_th']*P['landed']/1e6+0.3   # +контейнер/щит
    mwh=P['spare_mw']*8760*P['uptime']
    gas_mmbtu=mwh*3.412/0.38
    cf=[-capex]; eb=[]
    for y in range(4):
        rev=sum(th/1000*hp0*(1-P['erosion'])**((y*12+m)/12)*(P['halv'] if y*12+m>=10 else 1)*30.42*P['uptime'] for m in range(12))/1e6
        opex=gas_mmbtu*P['gas']/1e6+mwh*P['vom']/1e6+P['tor_asic']*capex+P['pool']*rev
        e=rev-opex; eb.append(round(e,2)); cf.append(e-max(0,0.34*(e-capex/4)))
    return dict(capex=round(capex,2),mwh=round(mwh),ebitda=eb,npv15=round(npv(0.15,cf),2),npv20=round(npv(0.20,cf),2),irr=(round(irr(cf)*100,1) if irr(cf) is not None else None),
                marg_cost_mwh=round((gas_mmbtu*P['gas']+mwh*P['vom'])/mwh,1))
out={f'{hp}_{sc}':hybrid(hp,sc) for hp in [27.7,39.87,50.0] for sc in ['pess','base','opt']}
json.dump(out,open('hybrid_out.json','w'),indent=1)
for k,v in out.items(): print(k,v)
