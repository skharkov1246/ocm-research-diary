import sys
sys.path.insert(0,'../nacn_ng')
src=open('../nacn_ng/model.py').read().split('print("=== (a)')[0]
exec(src)
def run(Q, netback, var, sold, capex, fixed, ramp_fixed=True, label=''):
    contrib=sold*(netback-var)/1e6
    e=contrib-fixed
    e1=(0.6*contrib-fixed) if ramp_fixed else 0.6*e
    f=[-0.4*capex,-0.6*capex,e1]+[e]*14
    i=irr(f)
    print(f"{label}: Q={Q} sold={sold:.0f} netback={netback:.0f} var={var:.0f} EBITDA={e:.2f} yr1={e1:.2f} pb={capex/e if e>0 else float('inf'):.1f} NPV15={npv(0.15,f):.1f} NPV20={npv(0.20,f):.1f} IRR={None if i is None else round(i*100,1)}")
    return e
r=option_b(10,'opt')
print('base opt 10kt reproduced EBITDA',round(r['ebitda'],2),'capex',round(r['capex'],2),'fixed',round(r['fixed'],2))
run(10,r['netback'],r['var'],r['sold'],r['capex'],r['fixed'],ramp_fixed=False,label='A as-is')
run(10,r['netback'],r['var'],r['sold'],r['capex'],r['fixed'],ramp_fixed=True,label='B ramp w/ fixed')
# price fix: Ghana mine-gate 2161-2253 (corpus) x0.97 - sea 240 - inland 40-80; or BF CIF 2240*0.97+130 -240 - Tema->Ouaga 120-200
nb_gh=[2253*0.97-240-40, 2161*0.97-240-80]
nb_bf=[2240*0.97+130-240-120, 2240*0.97+130-240-200]
print('netback Ghana gate',[round(x) for x in nb_gh],' BF',[round(x) for x in nb_bf])
nb=max(nb_gh+nb_bf)
print('best corrected netback',round(nb))
run(10,nb,r['var'],r['sold'],r['capex'],r['fixed'],label='C price fix')
# favourable fixes: NaOH 617+30, urea export parity 360 delivered
var2=r['var']-(0.857+0.094)*(680-647)-0.929*(400-360)
print('var with cheaper NaOH/urea',round(var2))
run(10,nb,var2,r['sold'],r['capex'],r['fixed'],label='D price fix + cheap reagents')
# volume cap to market: top 11 kt addressable *? sold 9.5 ok at top; center 5.4kt
# central, sold limited to 5.4 kt
import importlib
import os; os.chdir('../nacn_ng'); exec(open('central.py').read().split('for nh3 in')[0]); os.chdir('../skeptic_nacn')
c=option_b_t(10,0.5)
print('central 10kt as is EBITDA',round(c[6],2))
# central with market-limited sales 5.4 kt: recompute
t=0.5;Q=10;s=4
capex=c[0];var=c[1]
fixed_c = c[5]/1.015 - 8750*var/1e6
sold=5400; netc=1718
e=(sold*(netc-var)/1e6 - fixed_c)-0.015*(sold*var/1e6+fixed_c)
print('central 10kt sales capped 5.4kt EBITDA',round(e,2),'fixed',round(fixed_c,2))
# NH3 lever with 0.6 multiplier
e_nh3=5.17+0.6*(8.91-5.17); cap=49.7
f=[-0.4*cap,-0.6*cap,0.6*e_nh3]+[e_nh3]*14
print('NH3 x0.6 EBITDA',round(e_nh3,2),'pb',round(cap/e_nh3,1),'IRR',round(irr(f)*100,1),'NPV15',round(npv(0.15,f),1))
# NH3 with price fix
e_nh3b=e_nh3-(r['netback']-nb)*r['sold']/1e6
f=[-0.4*cap,-0.6*cap,0.6*e_nh3b]+[e_nh3b]*14
print('NH3 x0.6 + price fix EBITDA',round(e_nh3b,2),'pb',round(cap/e_nh3b,1),'IRR',round(irr(f)*100,1),'NPV15',round(npv(0.15,f),1))
# power capex line check at 2.5
print('power capex 2.5kt opt $M', round((0.19*1+0.05*2.5)*2.5*0.9,3))
# what netback needed for IRR=15% at 10kt opt (capex 41.7, fixed 2.74, var 1230)
lo,hi=1500,4000
for _ in range(60):
    m=(lo+hi)/2
    contrib=9500*(m-1230)/1e6; e=contrib-2.74
    f=[-0.4*41.71,-0.6*41.71,0.6*contrib-2.74]+[e]*14
    if npv(0.15,f)>0: hi=m
    else: lo=m
print('netback for NPV15=0 at 10kt opt', round(m), ' CIF-eq', round(m+240))
