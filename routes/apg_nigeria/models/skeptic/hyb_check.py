# Скептик: подвариант "ASIC на избытке captive" с падением дебита (как в captive.py decl) и профилем нагрузки
import sys; sys.path.insert(0,'../compute')
from mining_ai import npv, irr
P={'pess':dict(e=13.5,usd_th=17,landed=1.25,gas=1.0,vom=25,erosion=0.20,halv=0.50,pool=0.045,tor=0.05,decl=0.10,LF=0.60,avail=0.88),
   'base':dict(e=13.5,usd_th=13,landed=1.20,gas=0.5,vom=18,erosion=0.10,halv=0.60,pool=0.04,tor=0.04,decl=0.05,LF=0.72,avail=0.92),
   'opt': dict(e=13.5,usd_th=9,landed=1.15,gas=0.25,vom=12,erosion=0.0,halv=0.70,pool=0.025,tor=0.035,decl=0.0,LF=0.85,avail=0.94)}
MWN=4.64*0.96
def run(hp0,sc,size_mw,capture,asic_up=0.96,decl=None):
    p=P[sc]; d=p['decl'] if decl is None else decl
    th=size_mw*1e6/p['e']; capex=th*p['usd_th']*p['landed']/1e6+0.3
    cf=[-capex]; eb=[]; mw_list=[]
    for y in range(4):
        g=(1-d)**y
        surplus=MWN*p['avail']*max(0,g-p['LF'])          # средний избыток по газу, МВт
        used=min(size_mw*capture, surplus*capture)       # профиль: ASIC берёт долю capture
        mw_list.append(round(used,2))
        mwh=used*8760*asic_up
        rev=sum(used*1e6/p['e']/1000*hp0*(1-p['erosion'])**((y*12+m)/12)*(p['halv'] if y*12+m>=10 else 1)*30.42*asic_up for m in range(12))/1e6
        opex=mwh*(3.412/0.38*p['gas']+p['vom'])/1e6+p['tor']*capex+p['pool']*rev
        e=rev-opex; eb.append(round(e,2)); cf.append(e-max(0,0.34*(e-capex/4)))
    return dict(capex=round(capex,2),mw_used=mw_list,ebitda=eb,npv15=round(npv(0.15,cf),2),irr=(round(irr(cf)*100,1) if irr(cf) is not None else None))
size={'pess':0.8,'base':1.15,'opt':1.4}
for hp in [27.7,39.87,50.0]:
    for sc in ['base','opt']:
        print(hp,sc,'flat,no decl',run(hp,sc,size[sc],1.0,decl=0.0))
        print(hp,sc,'flat,decl   ',run(hp,sc,size[sc],1.0))
        print(hp,sc,'cap0.77,decl',run(hp,sc,size[sc],0.77))
        print(hp,sc,'cap0.5,decl ',run(hp,sc,size[sc],0.5))
print('--- opt sized to consistent surplus 0.63 MW')
for hp in [27.7,39.87,50.0]:
    for cap in [1.0,0.77]:
        print('RS',hp,cap,run(hp,'opt',0.63,cap))
