# Мини-GTL (ФТ) на факельном газе Нигерии: экономика 1/5/15 MMscf/д, скрин vs аудит NaCN.
import json
MMBTU_PER_MSCF=1.05
def npv(r,cf): return sum(c/(1+r)**t for t,c in enumerate(cf))
def irr(cf):
    lo,hi=-0.95,2.0
    f=lambda r:npv(r,cf)
    if f(lo)*f(hi)>0: return None
    for _ in range(200):
        m=(lo+hi)/2
        if f(lo)*f(m)<=0: hi=m
        else: lo=m
    return m
# сценарии: Y bpd на MMscf/д, доступность, цена дизеля $/bbl у ворот, нафта=k*дизель,
# капекс $/bpd по номиналу 100 bpd/MMscfd (с поправкой на масштаб), опекс-компоненты
SC={
 'screen': dict(Y=100,avail=0.95,diesel=120,naph_k=0.65,dshare=0.70,capex={1:10,5:50,15:120},
                opex_bbl=10.0,staff=0,sec=0,tor=0,ins=0,cat=0,gas=0.6,levy=0.0,decl=0.0,ramp=1.0,build=1),
 'opt':    dict(Y=95,avail=0.92,diesel=175,naph_k=0.65,dshare=0.70,capex={1:24,5:90,15:200},
                opex_bbl=0,staff=0.6,sec=0.25,tor=0.035,ins=0.010,cat=3.0,gas=0.25,levy=0.02,decl=0.0,ramp=0.8,build=2),
 'base':   dict(Y=80,avail=0.88,diesel=120,naph_k=0.65,dshare=0.70,capex={1:32,5:115,15:270},
                opex_bbl=0,staff=0.8,sec=0.35,tor=0.0425,ins=0.012,cat=4.5,gas=0.6,levy=0.025,decl=0.05,ramp=0.6,build=2),
 'pess':   dict(Y=65,avail=0.80,diesel=90,naph_k=0.60,dshare=0.70,capex={1:50,5:140,15:340},
                opex_bbl=0,staff=1.2,sec=0.5,tor=0.05,ins=0.015,cat=6.0,gas=1.0,levy=0.03,decl=0.10,ramp=0.5,build=3),
}
STAFF_SCALE={1:1.0,5:1.8,15:2.8}
TAX=0.34; LIFE=15
def run(sc,q,diesel=None,capex=None,wax=None,edti=False,tax=True):
    p=dict(SC[sc])
    if diesel is not None: p['diesel']=diesel
    C=capex if capex is not None else p['capex'][q]
    k=STAFF_SCALE[q]
    price=p['dshare']*p['diesel']+(1-p['dshare'])*p['naph_k']*p['diesel']
    if wax:  # вариант без гидрокрекинга: 40% воск, 40% дизель, 20% нафта; +капекс узла воска
        wprice,wcap=wax
        t_per_bbl_wax=0.14
        price=0.4*wprice*t_per_bbl_wax+0.4*p['diesel']+0.2*p['naph_k']*p['diesel']
        C=C+wcap*q**0.6
    cf=[]; nb=p['build']
    for t in range(nb): cf.append(-C/nb)
    rows=[]
    for y in range(1,LIFE+1):
        g=(1-p['decl'])**(y-1)*(p['ramp'] if y==1 else 1.0)
        bbl=p['Y']*q*p['avail']*365*g
        rev=bbl*price/1e6
        gas=MMBTU_PER_MSCF*1000*q*365*p['avail']*g*p['gas']/1e6
        opex=(bbl*(p['opex_bbl']+p['cat'])/1e6 + gas + p['staff']*k + p['sec']*k**0.7
              + (p['tor']+p['ins'])*C + p['levy']*rev)
        e=rev-opex
        dep=C/10 if y<=10 else 0
        tx=max(0,TAX*(e-dep)) if tax else 0
        if edti and y<=5: tx=max(0,tx-0.05*C)
        cf.append(e-tx); rows.append(dict(y=y,bbl=round(bbl),rev=round(rev,2),opex=round(opex,2),ebitda=round(e,2)))
    full=rows[1]  # первый полный год
    acc=0;pb=None
    for t,c in enumerate(cf):
        acc+=c
        if acc>=0 and pb is None and t>0: pb=t
    return dict(capex=round(C,1),price_bbl=round(price,1),bbl_y=full['bbl'],rev=full['rev'],opex=full['opex'],
                ebitda=full['ebitda'],usd_per_bbl_cash=round(full['opex']*1e6/full['bbl'],1),
                simple_pb=(round(C/full['ebitda'],1) if full['ebitda']>0 else None),pb_from_t0=pb,
                npv15=round(npv(0.15,cf),1),npv20=round(npv(0.20,cf),1),irr=(round(irr(cf)*100,1) if irr(cf) is not None else None))
def be_price(sc,q,r=0.15,**kw):
    lo,hi=20,1000
    for _ in range(80):
        m=(lo+hi)/2
        v=npv(r,[0]) # dummy
        res=run(sc,q,diesel=m,**kw)
        if (res['npv15'] if r==0.15 else res['npv20'])<0: lo=m
        else: hi=m
    return round(m,0)
out={}
for sc in SC:
    for q in (1,5,15):
        out[f'{sc}_q{q}']=run(sc,q)
        out[f'{sc}_q{q}_BEdiesel_npv15']=be_price(sc,q)
        out[f'{sc}_q{q}_BEdiesel_npv20']=be_price(sc,q,r=0.20)
# чувствительности на базе
sens={}
for q in (1,5,15):
    sens[f'base_q{q}_diesel200_crisis']=run('base',q,diesel=200)
    sens[f'base_q{q}_screen_capex']=run('base',q,capex=SC['screen']['capex'][q])
    sens[f'base_q{q}_EDTI']=run('base',q,edti=True)
    sens[f'base_q{q}_wax1200']=run('base',q,wax=(1200,3.0))
    sens[f'opt_q{q}_wax1500']=run('opt',q,wax=(1500,2.0))
    sens[f'opt_q{q}_capex_base']=run('opt',q,capex=SC['base']['capex'][q])
out['sens']=sens
json.dump(out,open('/tmp/claude-0/-home-user-ocm-research-diary/50ead64b-47e2-5bdd-b34b-09f6667b6ca9/scratchpad/gtl/gtl_out.json','w'),ensure_ascii=False,indent=1)
for k,v in out.items():
    if k!='sens': print(k,v)
for k,v in sens.items(): print(k,v)
