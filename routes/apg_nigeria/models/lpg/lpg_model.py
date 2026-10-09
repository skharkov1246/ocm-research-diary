# LPG/NGL extraction skid, Nigeria - screen vs NaCN-audit corrected
KMOL_PER_MMSCF = 1e6/379.49*0.45359  # kmol per MMscf (60F,14.696)
MW = dict(C1=16.04,C2=30.07,C3=44.10,iC4=58.12,nC4=58.12,iC5=72.15,nC5=72.15,C6=86.18,C7=100.2)
DENS_LB_GAL = dict(C3=4.23,iC4=4.69,nC4=4.87,iC5=5.20,nC5=5.26,C6=5.53,C7=5.73)
HHV_BTU_SCF = dict(C1=1010,C2=1770,C3=2516,iC4=3252,nC4=3262,iC5=4001,nC5=4009,C6=4756,C7=5503,N2=0,CO2=0,H2O=0)
comps = {
 # Niger Delta lean sample (Modelling FGR, IJOGCE 2022) [В]; pentanes "small" -> 0.1/0.08 [Д]
 'lean_ND': dict(C1=87.92,C2=4.65,C3=0.93,iC4=0.26,nC4=0.29,iC5=0.10,nC5=0.08),
 # Ashton-Jones 1998 mid: C3 2.3 (1.2-3.4), heavier 1.9 split C4 1.2 / C5+ 0.7 [Д]
 'mid_AJ': dict(C1=88.0,C2=4.6,C3=2.3,iC4=0.6,nC4=0.6,iC5=0.35,nC5=0.35),
 # Marginal Oilfield A, Ndunagu et al NJTD 19(1) 2022, Table 1 [П]
 'rich_NJTD': dict(N2=0.10,CO2=0.46,C1=82.72,C2=6.73,C3=5.68,iC4=1.56,nC4=1.60,iC5=0.40,nC5=0.26,C6=0.16,C7=0.20,H2O=0.14),
}
def analyse(c, rC3, rC4):
    k = KMOL_PER_MMSCF
    t = lambda x: c.get(x,0)/100*k*MW[x]/1000
    lpg_full = t('C3')+t('iC4')+t('nC4')
    c5p = sum(t(x) for x in ['iC5','nC5','C6','C7'])
    lpg = t('C3')*rC3 + (t('iC4')+t('nC4'))*rC4
    lbmol = 1e3/379.49
    gpm = sum(c.get(x,0)/100*lbmol*MW[x]/DENS_LB_GAL[x] for x in DENS_LB_GAL)
    hhv = sum(c.get(x,0)/100*HHV_BTU_SCF.get(x,0) for x in c)
    lpg_hhv = sum(c.get(x,0)/100*HHV_BTU_SCF[x] for x in ['C3','iC4','nC4'])
    c5_hhv = sum(c.get(x,0)/100*HHV_BTU_SCF[x] for x in ['iC5','nC5','C6','C7'])
    return dict(GPM_C3plus=round(gpm,2), HHV=round(hhv), lpg_full_t=round(lpg_full,2), lpg_rec_t=round(lpg,2),
                c5p_t=round(c5p,2), energy_share_LPG=round(lpg_hhv/hhv,3), energy_share_C5p=round(c5_hhv/hhv,3))
print("== composition -> t/d per 1 MMscf/d ==")
for n,c in comps.items():
    print(n, 'sum', round(sum(c.values()),2), 'refrig(C3 .70/C4 .95):', analyse(c,.70,.95), '| turboexp(.95/.99):', analyse(c,.95,.99)['lpg_rec_t'])
# NJTD simulated: LPG 163.7 t/d @60 MMscfd
print('NJTD sim t/MMscf', round(163.7/60,2), 'cond', round(33.19/60,2))
fC3 = 5.68*60; sC3 = 4.69*57.15; print('NJTD implied C3 recovery', round(1-sC3/fC3,2), 'iC4', round(1-0.44*57.15/(1.56*60),2), 'nC4', round(1-0.33*57.15/(1.60*60),2))
print('Otakikpo 60/12', 60/12, '| Kokori 74.96e6kg/yr /34.83', round(74.96e3/365/34.83,2))
print('NGFCP R1 600kt/700MMscfd', round(600e3/365/700,2), '| NGFCP R2 170kt/250-300', round(170e3/365/300,2), round(170e3/365/250,2))

import math
def npv(cfs, r): return sum(cf/(1+r)**i for i,cf in enumerate(cfs))
def irr(cfs):
    lo,hi=-0.99,3.0
    if npv(cfs,lo)*npv(cfs,hi)>0: return None
    for _ in range(200):
        m=(lo+hi)/2
        if npv(cfs,lo)*npv(cfs,m)<=0: hi=m
        else: lo=m
    return m
CEPCI = 800/603.1  # 2018 -> 2025 [Д]
def run(Q, case, capex_mult=None, price=None, yld=None, gasp=None, verbose=False):
    S = dict(
      screen=dict(y=4.0, price=750, cond_t=0.8, cond_p=400, LF=0.95, decl=0.03, gas=0.25, staff=0.25, sec=0.10, mnt=0.03, ins=0.01, misc=0.05, cmult=1.0, build=[0.5,0.5]),
      audit =dict(y=2.0, price=600, cond_t=0.4, cond_p=400, LF=0.80, decl=0.08, gas=1.00, staff=0.5 if Q<10 else 0.7, sec=0.3 if Q<10 else 0.5, mnt=0.045, ins=0.015, misc=0.10 if Q<10 else 0.2, cmult=2.6, build=[0.4,0.6]),
    )[case]
    if capex_mult is not None: S['cmult']=capex_mult
    if price is not None: S['price']=price
    if yld is not None: S['y']=yld
    if gasp is not None: S['gas']=gasp
    capex_screen = 12*CEPCI*(Q/12)**0.6
    capex = capex_screen*S['cmult']
    life=12; tax=0.34
    cfs=[-capex*b for b in S['build']]
    dep = capex/5
    rows=[]
    for yr in range(life):
        f = S['LF']*(1-S['decl'])**yr
        lpg_t = S['y']*Q*365*f
        rev = (lpg_t*S['price'] + S['cond_t']*Q*365*f*S['cond_p'])/1e6
        gas = S['gas']*Q*1000*365*f/1e6
        fixed = S['staff']+S['sec']+S['misc']+capex*(S['mnt']+S['ins'])
        opex = gas+fixed
        ebitda = rev-opex
        taxable = ebitda-(dep if yr<5 else 0)
        t = max(0,taxable*tax)
        edti = min(t, 0.05*capex) if yr<5 else 0
        cf = ebitda - t + edti
        cfs.append(cf); rows.append((yr+1,round(lpg_t),round(rev,2),round(opex,2),round(ebitda,2),round(cf,2)))
    # payback from start of operation
    cum=-capex; pb=None
    for i,r in enumerate(rows):
        if cum<0 and cum+r[5]>=0: pb=i+(-cum)/r[5]
        cum+=r[5]
    out=dict(Q=Q,case=case,capex=round(capex,1),capex_screen=round(capex_screen,1),rev1=rows[0][2],opex1=rows[0][3],ebitda1=rows[0][4],
             NPV15=round(npv(cfs,0.15),1),NPV20=round(npv(cfs,0.20),1),IRR=None if irr(cfs) is None else round(irr(cfs)*100,1),payback=None if pb is None else round(pb,1),lpg_t_yr1=rows[0][1])
    if verbose: print(rows)
    return out
print("\n== site economics ($M) ==")
for Q in [1,5,15]:
    for case in ['screen','audit']:
        print(run(Q,case))
    for m in [2.4,2.8]:
        print('  audit capex x',m, run(Q,'audit',capex_mult=m))
print("\n== sensitivities on audit, Q=15 and 5 ==")
for Q in [5,15]:
  for lab,kw in [('price450',dict(price=450)),('price900',dict(price=900)),('yield1.0',dict(yld=1.0)),('yield4',dict(yld=4.0)),('gas0.25',dict(gasp=0.25)),('capex x1.6',dict(capex_mult=1.6)),('y4 p750 x1.6',dict(yld=4,price=750,capex_mult=1.6))]:
    print(Q, lab, run(Q,'audit',**kw))
print(run(15,'audit',verbose=True))
# screen opex x1.7-2.4 check
for Q in [1,5,15]:
    s=run(Q,'screen'); a=run(Q,'audit')
    print('Q',Q,'screen opex',s['opex1'],'x1.7-2.4 ->',round(s['opex1']*1.7,2),round(s['opex1']*2.4,2),'bottom-up audit opex',a['opex1'])

print("\n\n===== V2: add-on (residual gas goes to anchor; LPG pays shrinkage by energy) =====")
MMBTU_PER_T = 47.3   # LPG HHV ~50 MJ/kg [В справочно]
def run2(Q, mode, y=2.0, price=600, cmult=2.6, gas_mmbtu=1.0, LF=0.80, decl=0.08, verbose=False):
    capex_screen = 12*CEPCI*(Q/12)**0.6
    capex = capex_screen*cmult
    staff = 0.5 if Q<10 else 0.7
    sec = 0.3 if Q<10 else 0.5
    misc = 0.1 if Q<10 else 0.2
    if mode=='addon':  # shared site: staff/security partly carried by anchor -> 50%
        staff*=0.5; sec*=0.5
    cfs=[-capex*0.4,-capex*0.6]; dep=capex/5; rows=[]
    for yr in range(12):
        f=LF*(1-decl)**yr
        lpg=y*Q*365*f; cond=0.2*y*Q*365*f
        rev=(lpg*price+cond*400)/1e6
        if mode=='standalone':
            gas=gas_mmbtu*1.05*Q*1000*365*f/1e6   # whole feed, $/MMBtu ~ $/Mscf*1.05
        else:
            gas=gas_mmbtu*(lpg+cond)*MMBTU_PER_T/1e6 + gas_mmbtu*0.03*1.05*Q*1000*365*f/1e6  # shrinkage + 3% own fuel
        opex=gas+staff+sec+misc+capex*(0.045+0.015)
        e=rev-opex; t=max(0,(e-(dep if yr<5 else 0))*0.34); ed=min(t,0.05*capex) if yr<5 else 0
        cfs.append(e-t+ed); rows.append(e)
    cum=-capex; pb=None
    for i,cf in enumerate(cfs[2:]):
        if cum<0 and cum+cf>=0: pb=i+(-cum)/cf
        cum+=cf
    r=irr(cfs)
    return dict(Q=Q,mode=mode,y=y,price=price,cmult=cmult,capex=round(capex,1),rev1=round((y*Q*365*LF*price+0.2*y*Q*365*LF*400)/1e6,2),ebitda1=round(rows[0],2),NPV15=round(npv(cfs,.15),1),NPV20=round(npv(cfs,.20),1),IRR=None if r is None else round(r*100,1),payback=None if pb is None else round(pb,1))
for Q in [1,5,15]:
    for mode in ['standalone','addon']:
        for (y,p) in [(2.0,600),(4.0,650)]:
            for cm in [2.4,2.8]:
                print(run2(Q,mode,y=y,price=p,cmult=cm))
print('-- what it takes (addon, Q=15) --')
for cm in [1.0,1.3,1.6,2.0]:
    for (y,p) in [(2.0,600),(3.0,650),(4.0,650),(4.0,750),(5.0,750)]:
        print(run2(15,'addon',y=y,price=p,cmult=cm,decl=0.05,LF=0.85))
print('-- addon Q=5 --')
for cm in [1.0,1.6,2.4]:
    for (y,p) in [(2.0,600),(4.0,650),(5.0,750)]:
        print(run2(5,'addon',y=y,price=p,cmult=cm,decl=0.05,LF=0.85))

print("\n===== benchmarks per 1 MMscf/d =====")
mw=4.64; ph=265.1; asic_capex=ph*1000*18/1e6; gens=mw*700/1e3
for hp in [27.7,50]:
    g1=ph*hp*365/1e6; g2=g1*0.7; om=(g1+g2)*0.15
    net=(g1+g2-om-asic_capex-gens*0.2)/2
    print('mining hashprice',hp,'gross1',round(g1,2),'EBITDA1',round(g1*0.85,2),'net/yr corpus method',round(net,2),'capex',round(asic_capex+gens+1,2), '$/MMBtu gross', round(ph*hp/1000,2))
for tar in [0.12,0.157,0.20]:
    rev=mw*0.9*8760*tar/1e3*1e3/1e3
    rev=mw*1000*0.9*8760*tar/1e6
    ebit=rev-mw*1000*0.9*8760*0.02/1e6-0.8-0.25*1.05*365/1e3*1000/1e3
    print('captive power tariff',tar,'rev',round(rev,2),'EBITDA~',round(ebit,2),'$/MMBtu gross',round(tar*0.38*293.07,1))
print('LPG $/MMBtu at 600/650/750:',[round(p/MMBTU_PER_T,1) for p in (600,650,750)])
print('LPG revenue per MMscf/d (y=2..4, LF .8, $600-750):', round(2*365*0.8*600/1e6,2), round(4*365*0.8*750/1e6,2))
