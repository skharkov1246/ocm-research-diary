# Флаер -> газопоршневая генерация -> продажа э/э промпотребителям (Нигерия). Все входы — допущения с диапазоном.
import json
FX=1330.0
MW_PER_MMSCFD=4.64      # корпус econ_vs_mining: 1000 Btu/scf, КПД 0.38
PARASITIC=0.04          # компрессия/подготовка/собств. нужды
MMBTU_PER_MMSCFD_Y=365*1000*1.05

def npv(r,cf): return sum(c/(1+r)**t for t,c in enumerate(cf))
def irr(cf):
    lo,hi=-0.99,3.0
    if npv(lo+1e-9,cf)*npv(hi,cf)>0: return None
    for _ in range(200):
        mid=(lo+hi)/2
        if npv(lo,cf)*npv(mid,cf)<=0: hi=mid
        else: lo=mid
    return mid

S={
 'low':  dict(usd_kw=1960,red=1.25,cond_fix=1.5,cond_var=0.3,line_km=40,usd_km=120e3,subst=0.8,dev=1.0,
              LF=0.60,avail=0.88,price_ngn=250,gas=1.0,tor=0.05,staff=0.45,sec=0.5,ins=0.015,baddebt=0.10,loss=0.04,decl=0.10,salv=0.0,dbd=0.30),
 'base': dict(usd_kw=1700,red=1.25,cond_fix=1.0,cond_var=0.2,line_km=25,usd_km=90e3,subst=0.6,dev=0.7,
              LF=0.72,avail=0.92,price_ngn=320,gas=0.5,tor=0.04,staff=0.35,sec=0.35,ins=0.012,baddebt=0.05,loss=0.03,decl=0.05,salv=0.15,dbd=0.15),
 'high': dict(usd_kw=1500,red=1.20,cond_fix=0.6,cond_var=0.15,line_km=12,usd_km=70e3,subst=0.5,dev=0.5,
              LF=0.85,avail=0.94,price_ngn=400,gas=0.25,tor=0.035,staff=0.30,sec=0.25,ins=0.010,baddebt=0.03,loss=0.025,decl=0.0,salv=0.25,dbd=0.05),
}
def run(q,p0,red=None,direct_share=1.0,disco_ngn=120,scale_kw=1.0,line_mult=1.0,staff_mult=1.0,years=10,gas_override=None):
    p=dict(p0)
    if red: p['red']=red
    mw_g=MW_PER_MMSCFD*q; mw_n=mw_g*(1-PARASITIC)
    gen_capex=mw_g*1000*p['usd_kw']*scale_kw*p['red']/1e6
    cond=p['cond_fix']+p['cond_var']*q
    line=p['line_km']*p['usd_km']*line_mult/1e6+p['subst']*line_mult
    capex=gen_capex+cond+line+p['dev']
    gas=gas_override if gas_override is not None else p['gas']
    price=(direct_share*p['price_ngn']+(1-direct_share)*disco_ngn)/FX  # $/kWh
    # DisCo-часть покупается почти плоско (LF~0.9), прямые клиенты — LF сценария
    LFeff=direct_share*p['LF']+(1-direct_share)*0.9
    cf=[-capex]; rows=[]
    for t in range(1,years+1):
        g=(1-p['decl'])**(t-1)
        sold=8760*mw_n*p['avail']*min(LFeff,g)*(1-p['loss'])/1e3  # GWh
        rev=sold*1e6*((direct_share*p['price_ngn']*(1-p['baddebt'])+(1-direct_share)*disco_ngn*(1-p.get('dbd',0.15)))/FX)/1e6       # $M
        gas_c=MMBTU_PER_MMSCFD_Y*q*min(LFeff,g)*p['avail']*gas/1e6 # платим за сожжённый в двигателях газ
        opex=gas_c+p['tor']*capex+p['staff']*staff_mult+p['sec']*staff_mult+p['ins']*capex
        e=rev-opex
        tax=max(0,0.34*(e-capex/years))
        c=e-tax+(p['salv']*gen_capex if t==years else 0)
        cf.append(c); rows.append((t,round(sold,1),round(rev,2),round(opex,2),round(e,2)))
    e1=rows[0][4]
    pb=None;acc=-capex
    for t,c in enumerate(cf[1:],1):
        acc+=c
        if acc>=0 and pb is None: pb=t-(acc/c)+0 if c else t
    return dict(q=q,capex=round(capex,2),gen=round(gen_capex,2),cond=round(cond,2),line=round(line,2),dev=p['dev'],
                price_usd_kwh=round(price,3),y1_GWh=rows[0][1],y1_rev=rows[0][2],y1_opex=rows[0][3],y1_ebitda=e1,
                simple_payback=round(capex/e1,1) if e1>0 else None, dcf_payback=round(pb,1) if pb else None,
                npv15=round(npv(0.15,cf),2),npv20=round(npv(0.20,cf),2),irr=(round(irr(cf)*100,1) if irr(cf) is not None else None))
out={}
for k,p in S.items():
    out[f'1MMscfd_{k}']=run(1,p)
# 5 и 15 MMscf/д: доля прямых клиентов падает, остальное — DisCo/embedded по ₦100-150
for k,p,ds5,ds15,dn in [('low',S['low'],0.3,0.1,100),('base',S['base'],0.5,0.2,120),('high',S['high'],0.7,0.35,150)]:
    out[f'5MMscfd_{k}']=run(5,p,direct_share=ds5,disco_ngn=dn,red=1.1,scale_kw=0.9,line_mult=2.0,staff_mult=1.8)
    out[f'15MMscfd_{k}']=run(15,p,direct_share=ds15,disco_ngn=dn,red=1.1,scale_kw=0.8,line_mult=4.0,staff_mult=2.5)
# чувствительности на 1 MMscf/д base
b=S['base']
sens={}
for nm,mod in [('price_210_BandA',dict(price_ngn=210)),('price_250',dict(price_ngn=250)),('price_450',dict(price_ngn=450)),
               ('LF_0.5',dict(LF=0.5)),('line_50km',dict(line_km=50)),('line_5km',dict(line_km=5)),
               ('gas_DBP_2.18',dict(gas=2.18)),('gas_2.68_comm',dict(gas=2.68)),('decl_15pct',dict(decl=0.15)),
               ('FX_2000_unindexed',dict(price_ngn=320*1330/2000)),('disco_only_120',dict(price_ngn=120,LF=0.9,baddebt=0.15)),('no_N+1',dict(red=1.0)),('EDTI_like_none',dict())]:
    p=dict(b);p.update(mod);sens[nm]=run(1,p)
out['sens_1MMscfd_base']={k:{kk:v[kk] for kk in ['capex','price_usd_kwh','y1_ebitda','simple_payback','npv15','npv20','irr']} for k,v in sens.items()}
# корпусный «скрин» без поправок: генсет $700/кВт + $1 млн инфра, цена = дизельный паритет ₦915, LF=1, O&M 15% выручки
mw=4.64; capex=mw*700/1e3+1.0; rev=mw*8760*915/FX/1e3; e=rev*0.85
out['screen_no_audit']=dict(capex=round(capex,2),rev=round(rev,2),ebitda=round(e,2),payback_yr=round(capex/e,2))
json.dump(out,open('/tmp/claude-0/-home-user-ocm-research-diary/50ead64b-47e2-5bdd-b34b-09f6667b6ca9/scratchpad/power/captive_out.json','w'),indent=1,ensure_ascii=False)
for k,v in out.items(): print(k,v)
