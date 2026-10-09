import math, importlib.util, sys, io, contextlib
spec=importlib.util.spec_from_file_location('pm','/tmp/claude-0/-home-user-ocm-research-diary/50ead64b-47e2-5bdd-b34b-09f6667b6ca9/scratchpad/pyro/pyro_model.py')
pm=importlib.util.module_from_spec(spec)
# avoid overwriting pyro_out.json: patch json.dump
import json
_d=json.dump; json.dump=lambda *a,**k:None
with contextlib.redirect_stdout(io.StringIO()): spec.loader.exec_module(pm)
json.dump=_d
# equilibrium X vs P (JANAF dfG CH4 1300K 52.26, 1400K 63.27 kJ/mol)
for T in (1338,1373):
    dG=52.26+(T-1300)/100*11.01
    Kp=math.exp(dG*1000/(8.314*T))
    print(T,'Kp',round(Kp,1),[ (P,round(math.sqrt(Kp/(Kp+4*P)),3)) for P in (1,5,10,15,20)])
# basis checks
mol_d=1e6*0.0283168/0.023645
print('t CH4/y per MMscfd',round(mol_d*16.04/1e6*365))
print('MW LHV',round(7300e6/16.04*802.3e3/3.1536e7/1e6,2),'MW el @0.38',round(7300e6/16.04*802.3e3/3.1536e7/1e6*0.38,2))
def r(sc,q=1,power=True,o=None):
    x=pm.run(sc,q,power,over=o); return (x['capex'],x['ebitda_y2'],x['npv15'],x['npv20'],x['irr'])
cases={'opt as is':{},
 'opt X=0.83 (20t melt ~15-20bar)':dict(X=0.83),
 'opt X=0.95 melt 70t':dict(melt_t=70),
 'opt c_net 220 (anchor ceiling)':dict(c_net=220),
 'opt X0.83+c220':dict(X=0.83,c_net=220),
 'opt X0.95 melt70 + c220':dict(melt_t=70,c_net=220),
 'opt X0.83+c220+no wet':dict(X=0.83,c_net=220,wet=1.0),
}
for q in (1,5,15):
    for k,o in cases.items(): print(q,k,r('high',q,True,o))
print('carbon-only opt X0.83 c220',r('high',1,False,dict(X=0.83,c_net=220)))
print('base, Bi 55',r('base',1,True,dict(bi=55)))
print('base FX2000 equiv price 320*1330/2000',r('base',1,True,dict(price_ngn=320*1330/2000)))
print('--- basis fix 7011 t, 4.40 MW')
pm.CH4_T_Y=7011.0; pm.MW_DIRECT=4.40
for q in (1,5,15):
    print(q,'base',r('base',q)); print(q,'opt corrected X0.83 c220 wet1',r('high',q,True,dict(X=0.83,c_net=220,wet=1.0)))
print('carbon-only base',r('base',1,False),'opt corr',r('high',1,False,dict(X=0.83,c_net=220,wet=1.0)))
x=pm.run('base',1); print('base y2',x['y2'],x['balance'])
