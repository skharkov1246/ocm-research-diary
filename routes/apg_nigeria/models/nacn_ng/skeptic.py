src=open('central.py').read().split('for nh3 in')[0]
exec(src)
def run(Q,t,nb_adj=0,wc_days=(30,30,60),decl=0.0,util=None,own=False,label=''):
    s=Q/2.5
    base = L(10.22*s**0.62, 11.64*s**0.72*(1.18 if Q>2.5 else 1.0), t)
    capex=base*(1+L(0.35,0.60,t))+1.51*s**0.7*L(0.8,1.5,t)+L(1.0,2.5,t)*s**0.6+(0.19*s+0.05*Q)*2.5*L(0.9,1.4,t)+L(0.3,0.8,t)+L(0.5,2.0,t)
    if own: capex+=L(8,20,t)*(Q/10)**0.65
    urea=L(0.929*400,1.071*520,t); alk=L((0.857+0.094)*680,(0.90+0.188)*830,t)
    if own: urea=L(0.511,0.589,t)*L(35,55,t)*0.948*L(0.25,1.0,t); alk=L(0.857*680,0.90*830,t)
    var=urea+alk+L(37*0.25,42*1.0,t)+L(700*0.015,850*0.025,t)+L(10,30,t)+L(20,40,t)+L(50,100,t)+L(40,78,t)+L(22,38,t)+L(50,110,t)
    nloc=L(18,24,t)*s**0.4+(6 if own else 0)
    fixed=nloc*L(0.012,0.025,t)+L(0.30,0.80,t)+L(0.4,1.0,t)*s**0.3+capex*L(0.035,0.05,t)
    price=L(2240*0.97+130,2036*0.92,t)
    nb0=price-L(240,500,t)+nb_adj
    u=util if util else L(0.95,0.80,t)
    sold=Q*1000*u
    flows=[-0.4*capex,-0.6*capex]
    es=[]
    for y in range(15):
        nb=nb0-price*(1-(1-decl)**y) if decl else nb0
        k=0.6 if y==0 else 1
        opex=(sold*k*var/1e6+fixed)*(1+0.03*t)
        e=sold*k*nb/1e6-opex; es.append(e); flows.append(e)
    # working capital at full rate: raw 30d, product 30d, receivables 60d
    rev=sold*nb0/1e6; vc=sold*var/1e6
    wc=vc*(wc_days[0]+wc_days[1])/365+rev*wc_days[2]/365
    flows[2]-=wc; flows[-1]+=wc
    i=irr(flows)
    print(f"{label:38s} Q={Q} t={t}: capex {capex:.1f} WC {wc:.1f} EBITDA y2 {es[1]:.2f} y5 {es[4]:.2f} pb {capex/es[1] if es[1]>0 else float('inf'):.1f} NPV15 {npv(0.15,flows):.1f} NPV20 {npv(0.20,flows):.1f} IRR {None if i is None else round(i*100,1)}")
run(10,0,label='модель как есть + оборотка')
run(10,0,nb_adj=-201,label='опт: + фрахт Тема->Уагадугу $201 [П-корп]')
run(10,0,nb_adj=-201,decl=0.03,label='то же + цена -3%/г (тренд CIF Гана)')
run(10,0,nb_adj=-201,util=0.55,label='то же, сбыт 5,5 кт (центр рынка)')
run(10,0.5,util=0.55,label='центр, сбыт 5,5 кт')
run(5,0,nb_adj=-201,label='5 кт опт, фрахт исправлен')
run(10,0,own=True,label='свой NH3 как есть + оборотка')
run(10,0,own=True,nb_adj=-201,label='свой NH3 + фрахт исправлен')
# gas balance own NH3
for t in (0,0.5):
    sold=10000*L(0.95,0.80,t); g1=sold*L(37,42,t); g2=sold*L(0.511,0.589,t)*L(35,55,t)*0.948
    print('газ NaCN+NH3, MMscf/d:', round(g1/1.05/1e3/365,2), '+', round(g2/1.05/1e3/365,2), '=', round((g1+g2)/1.05/1e3/365,2))
