import importlib.util,sys
spec=importlib.util.spec_from_file_location('m','model.py')
src=open('model.py').read().split('print("=== (a)')[0]
exec(src)
def L(a,b,t): return a+(b-a)*t
def option_b_t(Q,t, own_nh3=False, price_override=None):
    s=Q/2.5
    base = L(10.22*s**0.62, 11.64*s**0.72*(1.18 if Q>2.5 else 1.0), t)
    solid = base*L(0.35,0.60,t)
    loc = 1.51*s**0.7*L(0.8,1.5,t)
    gascond=L(1.0,2.5,t)*s**0.6
    power=(0.19*s+0.05*Q)*2.5*L(0.9,1.4,t)
    sec=L(0.3,0.8,t); tie=L(0.5,2.0,t)
    nh3cap = (L(8,20,t)*(Q/10)**0.65 if own_nh3 else 0)  # [Д] scaled from 4.4-11.7 at 3.5 t/d -> ~14-16 t/d
    capex=base+solid+loc+gascond+power+sec+tie+nh3cap
    urea = L(0.929*400,1.071*520,t)
    if own_nh3: urea = L(0.511,0.589,t)*L(35,55,t)*0.948*L(0.25,1.0,t)  # GJ/t NH3 -> MMBtu
    alk=L((0.857+0.094)*680,(0.90+0.188)*830,t)
    if own_nh3: alk=L(0.857*680,0.90*830,t)  # no CO2 ballast
    var=urea+alk+L(37*0.25,42*1.0,t)+L(700*0.015,850*0.025,t)+L(10,30,t)+L(20,40,t)+L(50,100,t)+L(40,78,t)+L(22,38,t)+L(50,110,t)
    nloc=L(18,24,t)*s**0.4 + (6 if own_nh3 else 0)
    staff=nloc*L(0.012,0.025,t)+L(0.30,0.80,t)
    fixed=staff+L(0.4,1.0,t)*s**0.3+capex*L(0.035,0.05,t)
    price = price_override if price_override else L(2240*0.97+130,2036*0.92,t)
    netback=price-L(240,500,t)
    sold=Q*1000*L(0.95,0.80,t)
    opex=sold*var/1e6+fixed; opex*=1+0.03*t
    e=sold*netback/1e6-opex
    f=proj(capex,e); i=irr(f)
    return capex,var,opex*1e6/sold,netback,sold*price/1e6,opex,e,npv(0.15,f),npv(0.20,f),i
for nh3 in (False,True):
  for Q in (2.5,5,10):
    for t in (0.0,0.5):
        c,v,cpt,nb,rev,ox,e,n15,n20,i=option_b_t(Q,t,nh3)
        print(f"NH3own={nh3} Q={Q} t={t}: capex {c:.1f} var {v:.0f} cash {cpt:.0f}/t netback {nb:.0f} rev {rev:.1f} opex {ox:.1f} EBITDA {e:.2f} pb {c/e if e>0 else float('inf'):.1f} NPV15 {n15:.1f} NPV20 {n20:.1f} IRR {None if i is None else round(i*100,1)}")
# Ghana-audit replay: what pre-audit screen would have said
print("pre-audit screen in Nigeria: capex 4.2-6.0 at 2.5kt, cash 1228 -> margin at netback 2063-1718:")
for cap in (4.2,6.0):
  for nb in (2063,1718):
    e=2500*(nb-1228)/1e6; print(cap,nb,round(e,2),round(cap/e,1))
