import itertools, json
CH4=7300.0
def mass_y(Y): return Y*28/32.08
MMBTU=365000*1.05
def npv(capex,ebitda,r,life=12,build=2,tax=0.30,dep_years=10):
    # capex spread over build years, then life years ebitda, simple tax on (ebitda-dep)
    v=0
    for t in range(build): v-= capex/build/(1+r)**t
    for t in range(life):
        dep=capex/dep_years if t<dep_years else 0
        tx=max(0,(ebitda-dep))*tax
        v+= (ebitda-tx)/(1+r)**(build+t)
    return v
def irr(capex,ebitda,**k):
    lo,hi=-0.9,2.0
    if npv(capex,ebitda,0.0,**k)<0: return None
    for _ in range(100):
        m=(lo+hi)/2
        if npv(capex,ebitda,m,**k)>0: lo=m
        else: hi=m
    return m
# 1) corpus twin reproduction
def twin(Y=0.35,olig=0.85,price=3200,opex_t=650,fix=0.15,gas=0.219,capex=5.0,opx_mult=1.0,cap_mult=1.0):
    t=CH4*mass_y(Y)*olig
    rev=t*price/1e6
    opex=(t*opex_t/1e6+fix)*opx_mult
    e=rev-gas-opex
    return dict(t=round(t),rev=round(rev,2),opex=round(opex,2),ebitda=round(e,2),capex=capex*cap_mult,pb=round(capex*cap_mult/e,1) if e>0 else None)
print('twin base',twin())
for om,cm in [(1.7,2.4),(2.4,2.8)]:
    print('audit',om,cm,twin(opx_mult=om,cap_mult=cm))
# 2) bottom-up Nigeria
def bottom(Y,olig,price,capex,gas_usd_mmbtu,scale=1.0,var_t=None,staff=None,sec=None,logi=None,maint=None,ins=0.015,parasitic=0.12):
    # parasitic: share of gas burned for power/compression (0.10-0.15)
    t=CH4*(1-parasitic)*mass_y(Y)*olig*scale
    rev=t*price/1e6
    gas=MMBTU*scale*gas_usd_mmbtu/1e6
    var=t*var_t/1e6        # catalysts/MAO/H2/BF3/chemicals/utilities
    lg=t*logi/1e6
    fixed=staff+sec+capex*maint+capex*ins
    e=rev-gas-var-lg-fixed
    return dict(t=round(t),rev=round(rev,2),gas=round(gas,2),var=round(var,2),log=round(lg,2),fixed=round(fixed,2),ebitda=round(e,2),capex=round(capex,1),
                pb=round(capex/e,1) if e>0 else None,npv15=round(npv(capex,e,0.15),1),npv20=round(npv(capex,e,0.20),1),irr=(round(irr(capex,e)*100,1) if e>0 and irr(capex,e) else None))
S={
 'A_corpus_audited (Y35,olig.85,$3200,capex12)':dict(Y=0.35,olig=0.85,price=3200,capex=12.0,gas_usd_mmbtu=0.25,var_t=250,staff=0.6,sec=0.3,logi=150,maint=0.035),
 'A2_corpus_audited_hi (capex14, worse costs)':dict(Y=0.35,olig=0.85,price=3200,capex=14.0,gas_usd_mmbtu=1.0,var_t=400,staff=0.9,sec=0.5,logi=250,maint=0.05),
 'B_hexenePAO_discount (Y35,$2200)':dict(Y=0.35,olig=0.80,price=2200,capex=13.0,gas_usd_mmbtu=0.5,var_t=300,staff=0.75,sec=0.4,logi=200,maint=0.04),
 'C_lowY (Y25,olig.75,$2200)':dict(Y=0.25,olig=0.75,price=2200,capex=13.0,gas_usd_mmbtu=0.5,var_t=300,staff=0.75,sec=0.4,logi=200,maint=0.04),
 'D_state_of_art (Y20,olig.70,$2200)':dict(Y=0.20,olig=0.70,price=2200,capex=13.0,gas_usd_mmbtu=0.5,var_t=300,staff=0.75,sec=0.4,logi=200,maint=0.04),
 'E_sell_1hexene (Y35, LAO 0.88, $1300, capex 9)':dict(Y=0.35,olig=0.88,price=1300,capex=9.0,gas_usd_mmbtu=0.5,var_t=200,staff=0.6,sec=0.4,logi=200,maint=0.04),
 'F_against_self (Y25,olig.75,$2200,capex22 miniGTL-analog)':dict(Y=0.25,olig=0.75,price=2200,capex=22.0,gas_usd_mmbtu=1.0,var_t=400,staff=0.9,sec=0.5,logi=250,maint=0.05),
 'G_best_case_with_decent_capex (Y35,.85,$3200,capex 22)':dict(Y=0.35,olig=0.85,price=3200,capex=22.0,gas_usd_mmbtu=0.5,var_t=300,staff=0.75,sec=0.4,logi=200,maint=0.04),
}
res={}
for k,v in S.items():
    r=bottom(**v); res[k]=r; print(k,r)
# scale 1/5/15 for scenario A and C: capex exponent 0.7, staff ^0.4, sec ^0.3
print('--- scale')
for name in ['A_corpus_audited (Y35,olig.85,$3200,capex12)','C_lowY (Y25,olig.75,$2200)']:
    v=S[name]
    for sc in (1,5,15):
        vv=dict(v); vv['capex']=v['capex']*sc**0.7; vv['staff']=v['staff']*sc**0.4; vv['sec']=v['sec']*sc**0.3; vv['scale']=sc
        r=bottom(**vv); print(name,sc,r)
# mining benchmark
def mining(hp, om=0.15, erosion=0.30, life=2):
    mw=1e6*1000*0.293/24*0.38/1e6
    th=mw*1e6/17.5; ph=th/1000
    g=ph*hp*365/1e6; tot=0; gg=g
    for _ in range(life): tot+=gg; gg*=(1-erosion)
    asic=th*18/1e6; gen=mw*1000*700/1e6; capex=asic+gen+1.0
    ebitda=tot*(1-om)/life
    net=(tot*(1-om)-asic-gen*life/10)/life
    return dict(mw=round(mw,2),ph=round(ph,1),gross1=round(g,2),ebitda_avg=round(ebitda,2),net=round(net,2),capex=round(capex,2))
for hp in (27.7,50): print('mining',hp,mining(hp))
# captive power
for price in (0.10,0.13,0.16):
    mw=4.64; gwh=mw*0.9*8760/1000
    rev=gwh*price; gas=MMBTU*0.25/1e6; gas2=MMBTU*1.0/1e6; om=gwh*0.02
    capex=(3.25+2.0)*1.4
    e_lo=rev-gas2-om-0.4; e_hi=rev-gas-om-0.25
    print('captive',price,round(gwh,1),round(rev,2),round(e_lo,2),round(e_hi,2),round(capex,2),round(npv(capex,e_lo,0.15,life=12),1),round(npv(capex,e_hi,0.20,life=12),1))
print('--- mining NPV 12y, ASIC refresh every 2y, tax 30%')
def mining_npv(hp,r,tax=0.30):
    mw=4.64; th=mw*1e6/17.5; ph=th/1000; asic=th*18/1e6; gen=mw*700/1e3; infra=1.0
    v=-(gen+infra)
    for y in range(12):
        if y%2==0: v-=asic/(1+r)**y
        g=ph*hp*365/1e6*(1-0.30)**(y%2)
        e=g*0.85 - 0.1  # gas $0.25/MMBtu ~0.1
        dep=asic/2+ (gen+infra)/10
        v+=(e-max(0,e-dep)*tax)/(1+r)**(y+0.5)
    return round(v,1)
for hp in (27.7,50):
    print(hp,mining_npv(hp,0.15),mining_npv(hp,0.20))
print('--- captive NPV (no build lag beyond 1y)')
