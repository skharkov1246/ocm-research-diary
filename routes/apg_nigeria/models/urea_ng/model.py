# Мини-аммиак/мочевина на факеле, Нигерия. Все входы — [Д] с диапазоном, см. комментарии.
import json
def npv(r,cf): return sum(c/(1+r)**t for t,c in enumerate(cf))
def irr(cf):
    lo,hi=-0.9,2.0
    if npv(lo,cf)*npv(hi,cf)>0: return None
    for _ in range(200):
        m=(lo+hi)/2
        if npv(lo,cf)*npv(m,cf)<=0: hi=m
        else: lo=m
    return m
MMBTU_HHV_PER_MMSCF=1050.0   # 1.03-1.1 MMBtu/Mscf (контекст C)
# сценарии: opt / mid / pes
S={
 'opt':dict(gj_t=24.0, days=340, capex_k=1.0, price=480, gas=0.25, staff=1.0, sec=0.4, tor=0.035, ins=0.005, ga=0.3, cons=20, decl=0.0, build=2),
 'mid':dict(gj_t=27.0, days=330, capex_k=1.0, price=400, gas=0.60, staff=1.0, sec=0.7, tor=0.042, ins=0.0075,ga=0.4, cons=25, decl=0.03, build=3),
 'pes':dict(gj_t=30.0, days=310, capex_k=1.0, price=300, gas=1.00, staff=1.0, sec=1.0, tor=0.05, ins=0.01, ga=0.5, cons=30, decl=0.08, build=3),
}
# капекс $M по размеру и сценарию ([О] триангуляция, см. derivation)
CAPEX={1:dict(opt=43,mid=60,pes=91),5:dict(opt=120,mid=170,pes=240),15:dict(opt=250,mid=340,pes=450)}
# персонал $M/г: [Д] 1 MMscf/д 40-50 местных+3-4 экспата; растёт медленно с масштабом
STAFF={1:dict(opt=1.5,mid=2.0,pes=2.5),5:dict(opt=2.5,mid=3.0,pes=3.5),15:dict(opt=4.0,mid=5.0,pes=6.0)}
SEC={1:1.0,5:1.4,15:2.0}   # множитель к охране/ОО
GJ_LHV_PER_MMBTU_HHV=1.055*0.90
def run(q,sc,price=None,gas=None,years=20,tax=0.30):
    p=S[sc]; capex=CAPEX[q][sc]
    pr=price if price is not None else p['price']; g=gas if gas is not None else p['gas']
    gj_day=q*MMBTU_HHV_PER_MMSCF*GJ_LHV_PER_MMBTU_HHV
    tpd=gj_day/p['gj_t']; tpy0=tpd*p['days']/1e3   # кт/г
    mmbtu_t=q*MMBTU_HHV_PER_MMSCF/tpd
    cf=[]; # стройка: капекс равными долями
    for i in range(p['build']): cf.append(-capex/p['build'])
    rows=[]
    for t in range(1,years+1):
        f=max((1-p['decl'])**(t-1),0.0)
        load=f if f>=0.6 else 0.0   # ниже 60% турндауна установка стоит
        tpy=tpy0*min(load,1.0)
        rev=tpy*pr/1e3
        var=tpy*(mmbtu_t*g+p['cons'])/1e3
        fix=STAFF[q][sc]+p['sec']*SEC[q]+p['tor']*capex+p['ins']*capex+p['ga']
        e=rev-var-fix
        dep=capex/10 if t<=10 else 0
        cf.append(e-max(0,tax*(e-dep)))
        rows.append(dict(t=t,kt=round(tpy,1),rev=round(rev,2),var=round(var,2),fix=round(fix,2),ebitda=round(e,2)))
    r1=rows[0]
    cash_t=(r1['var']+r1['fix'])/r1['kt']*1e3 if r1['kt'] else None
    ir=irr(cf)
    return dict(q=q,sc=sc,capex=capex,tpd=round(tpd,1),kt_y=round(tpy0,1),mmbtu_t=round(mmbtu_t,1),price=pr,gas=g,
                rev1=r1['rev'],opex1=round(r1['var']+r1['fix'],2),fix1=r1['fix'],ebitda1=r1['ebitda'],cash_usd_t=round(cash_t) if cash_t else None,
                capex_usd_per_tpy=round(capex*1e3/tpy0),payback=(round(capex/r1['ebitda'],1) if r1['ebitda']>0 else None),
                npv15=round(npv(0.15,cf),1),npv20=round(npv(0.20,cf),1),irr=(round(ir*100,1) if ir is not None else None),
                ebitda_sum20=round(sum(x['ebitda'] for x in rows),1))
out={}
for q in (1,5,15):
    for sc in ('opt','mid','pes'):
        out[f'{q}_{sc}']=run(q,sc)
# чувствительность: цена 276 (2024 FOB), 490 (2025 FOB avg), 640000 NGN ex-PH=481
for q in (1,5,15):
    for pr in (276,400,481,600):
        out[f'{q}_mid_price{pr}']=run(q,'mid',price=pr)
    out[f'{q}_mid_gasDBP2.18']=run(q,'mid',gas=2.18)
# цена, нужная для NPV15=0 (mid)
def be(q,sc,r=0.15):
    lo,hi=0,5000
    for _ in range(60):
        m=(lo+hi)/2
        cf=run(q,sc,price=m)['npv15' if r==0.15 else 'npv20']
        if cf>0: hi=m
        else: lo=m
    return round(m)
for q in (1,5,15):
    for sc in ('opt','mid','pes'):
        out[f'breakeven_price_npv15_{q}_{sc}']=be(q,sc)
# аудит-поправка к ряду двойника skid_twin (урея: капекс $6M, опекс $140/т, 8030 т/г, $380/т)
tw=dict(capex=6.0,opex_t=140,tpy=8030,price=380)
aud={}
for k,(ck,ok) in {'low':(2.4,1.7),'high':(2.8,2.4)}.items():
    cap=tw['capex']*ck; opx=tw['opex_t']*ok
    e=tw['tpy']*(tw['price']-opx)/1e6
    aud[k]=dict(capex=round(cap,1),opex_t=round(opx),ebitda=round(e,2),payback=round(cap/e,1))
out['audit_skid_twin']=aud
# майнинг по формуле корпуса при 27.7 и 50
def mining(hp):
    mw=1e6*1000*0.293/24*0.38/1e6; th=mw*1e6/17.5; ph=th/1000
    g1=ph*hp*365/1e6; asic=th*18/1e6; gen=mw*1000*700/1e6; capex=asic+gen+1.0
    gl=g1+g1*0.7; om=gl*0.15; net=(gl-om-asic-gen*0.2)/2
    return dict(hp=hp,mw=round(mw,2),ph=round(ph,1),gross1=round(g1,2),ebitda1=round(g1*0.85,2),capex=round(capex,2),net_per_yr=round(net,2))
out['mining']=[mining(27.7),mining(50)]
json.dump(out,open('/tmp/claude-0/-home-user-ocm-research-diary/50ead64b-47e2-5bdd-b34b-09f6667b6ca9/scratchpad/urea_ng/out.json','w'),ensure_ascii=False,indent=1)
for k,v in out.items(): print(k,v)
