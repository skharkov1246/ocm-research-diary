# Флаер -> CNG mother station на площадке -> каскады на тягачах -> промпотребители/дочерние АГНКС (Нигерия)
# Все входы помечены: [П] первичный, [В] вторичный, [О] расчёт, [Д] допущение
import json, math
FX=1330.0                      # ₦/$ [Д], NFEM окт-2026 1327-1333 [В]
SCF_PER_SCM=35.31              # [П] определение
MMBTU_PER_MSCF=1.03            # после контроля точки росы [Д 1.02-1.06]
MMBTU_PER_SCM=SCF_PER_SCM/1000*MMBTU_PER_MSCF   # 0.03637 [О]
SCM_PER_MMSCF=1e6/SCF_PER_SCM  # 28 317 [О]
def ngn_scm_to_usd_mmbtu(p): return p/FX/MMBTU_PER_SCM

def npv(r,cf): return sum(c/(1+r)**t for t,c in enumerate(cf))
def irr(cf):
    lo,hi=-0.95,3.0
    if npv(lo,cf)*npv(hi,cf)>0: return None
    for _ in range(200):
        m=(lo+hi)/2
        if npv(lo,cf)*npv(m,cf)<=0: hi=m
        else: lo=m
    return m

S={
 # 'screen' = как выглядел бы промоутерский скрин (пакет оборудования, без поправок аудита)
 'screen': dict(price=520, util=0.90, sell=0.93, gas=0.25, st_pkg=2.0, cmult=1.0, trailer=0.30, tractor=0.10, cust_kit=0.08,
                km=100, speed=45, drive_h=16, cyc_extra=3, usable_scm=7500, cust_scm=6000, opx_km=0.7,
                staff=0.30, sec=0.10, tor=0.03, ins=0.01, bad=0.02, overhead=0.10, decl=0.03, dev=0.3),
 # аудит NaCN применён: капитал x2.4-2.8 к пакету, выход по нижней границе, спрос /1.3-2.7, персонал с экспатом, ТОиР 3.5-5%
 'low':    dict(price=330, util=0.40, sell=0.85, gas=1.00, st_pkg=2.0, cmult=2.8, trailer=0.40, tractor=0.13, cust_kit=0.20,
                km=150, speed=30, drive_h=12, cyc_extra=4, usable_scm=6500, cust_scm=3500, opx_km=1.2,
                staff=0.80, sec=0.50, tor=0.05, ins=0.020, bad=0.08, overhead=0.30, decl=0.10, dev=1.0),
 'base':   dict(price=480, util=0.65, sell=0.88, gas=0.60, st_pkg=2.0, cmult=2.6, trailer=0.33, tractor=0.11, cust_kit=0.15,
                km=100, speed=35, drive_h=14, cyc_extra=3.5, usable_scm=7000, cust_scm=4500, opx_km=0.95,
                staff=0.55, sec=0.30, tor=0.045, ins=0.015, bad=0.05, overhead=0.20, decl=0.06, dev=0.7),
 'high':   dict(price=620, util=0.80, sell=0.90, gas=0.25, st_pkg=2.0, cmult=2.4, trailer=0.30, tractor=0.10, cust_kit=0.12,
                km=60, speed=40, drive_h=16, cyc_extra=3, usable_scm=7500, cust_scm=6000, opx_km=0.8,
                staff=0.40, sec=0.20, tor=0.035, ins=0.012, bad=0.03, overhead=0.15, decl=0.03, dev=0.5),
}
def run(Q, p0, years=10, tax=0.34, **ov):
    p=dict(p0); p.update(ov)
    # ---- станция (capex масштаб Q^0.7) ----
    station=p['st_pkg']*Q**0.7*p['cmult']                       # $M
    # ---- объёмы ----
    inlet_scm=Q*SCM_PER_MMSCF
    sold_scm_d=inlet_scm*p['sell']*p['util']                   # год 1
    loads=sold_scm_d/p['usable_scm']
    trip_h=2*p['km']/p['speed']+p['cyc_extra']                  # рейс туда-обратно + сцепка/досмотр
    tractors=math.ceil(loads*trip_h/p['drive_h']*1.15)          # +15% резерв
    n_cust=math.ceil(sold_scm_d/p['cust_scm'])
    # прицепы: на заправке (заполняются ~12ч каждый в параллели), в пути, по одному у каждого клиента, +15%
    trailers=math.ceil((loads*12/24 + loads*trip_h/24 + n_cust)*1.15)
    fleet=trailers*p['trailer']+tractors*p['tractor']+n_cust*p['cust_kit']
    capex=station+fleet+p['dev']
    rows=[];cf=[-capex]
    for t in range(1,years+1):
        g=(1-p['decl'])**(t-1)
        sold=sold_scm_d*g*365                                    # scm/год
        rev=sold*p['price']/FX/1e6*(1-p['bad'])
        gas=inlet_scm*p['util']*g*365*MMBTU_PER_SCM*p['gas']/1e6   # платим за весь поданный газ (включая собственное топливо)
        km_y=loads*g*2*p['km']*365
        truck=km_y*p['opx_km']/1e6
        drivers=tractors*2*6000/1e6                              # 2 водителя на тягач по $6k/г [Д]
        fixed=p['staff']+p['sec']+p['overhead']+drivers+capex*(p['tor']+p['ins'])
        opex=gas+truck+fixed
        e=rev-opex
        dep=capex/7
        taxv=max(0,(e-dep)*tax)
        repl=(tractors*p['tractor'] if t==6 else 0)              # замена тягачей на 6-й год
        cf.append(e-taxv-repl); rows.append(dict(t=t,rev=round(rev,2),opex=round(opex,2),gas=round(gas,2),truck=round(truck,2),fixed=round(fixed,2),ebitda=round(e,2)))
    e1=rows[0]['ebitda']
    cum=-capex;pb=None
    for t,c in enumerate(cf[1:],1):
        if pb is None and cum+c>=0 and c>0: pb=t-1+(-cum)/c
        cum+=c
    r1=rows[0]; mmbtu1=sold_scm_d*365*MMBTU_PER_SCM/1e6
    return dict(Q=Q,price_ngn_scm=p['price'],price_usd_mmbtu=round(ngn_scm_to_usd_mmbtu(p['price']),2),
        sold_scm_d=round(sold_scm_d),sold_MMBtu_y_M=round(mmbtu1,3),loads_d=round(loads,1),trip_h=round(trip_h,1),
        tractors=tractors,trailers=trailers,customers=n_cust,
        capex=round(capex,2),station=round(station,2),fleet=round(fleet,2),
        rev1=r1['rev'],opex1=r1['opex'],gas1=r1['gas'],truck1=r1['truck'],fixed1=r1['fixed'],ebitda1=e1,
        cash_cost_usd_mmbtu=round(r1['opex']/mmbtu1,2),
        simple_payback=(round(capex/e1,1) if e1>0 else None),dcf_payback=(round(pb,1) if pb else None),
        npv15=round(npv(0.15,cf),2),npv20=round(npv(0.20,cf),2),irr=(round(irr(cf)*100,1) if irr(cf) is not None else None))
out={'units':dict(MMBTU_PER_SCM=round(MMBTU_PER_SCM,5),SCM_PER_MMSCF=round(SCM_PER_MMSCF),
      price_map={n:round(ngn_scm_to_usd_mmbtu(n),2) for n in [200,250,300,320,380,400,450,470,520]})}
for k in ['screen','low','base','high']:
    out[f'1_{k}']=run(1,S[k])
# масштаб: при 5 и 15 MMscf/д нужно больше клиентов, средняя дистанция растёт, утилизация падает
for Q,km_add,util_mult in [(5,60,0.85),(15,150,0.6)]:
    for k in ['low','base','high']:
        p=dict(S[k]); p['km']=p['km']+km_add; p['util']=p['util']*util_mult
        p['staff']*=1.8 if Q==5 else 2.8; p['sec']*=1.6 if Q==5 else 2.5; p['overhead']*=1.5 if Q==5 else 2.2
        out[f'{Q}_{k}']=run(Q,p)
# чувствительности 1 MMscf/д base
b=S['base']; sens={}
for nm,mod in [('km_30',dict(km=30)),('km_200',dict(km=200)),('km_300',dict(km=300)),
               ('price_300',dict(price=300)),('price_400_mobility',dict(price=400)),('price_250_NIPCOlike',dict(price=250)),('price_520_industry_breakeven',dict(price=520)),('price_680_14usd',dict(price=680)),('price_950_half_diesel',dict(price=950)),
               ('util_0.45',dict(util=0.45)),('util_0.85',dict(util=0.85)),('gas_0.25',dict(gas=0.25)),('gas_1.57_GBI',dict(gas=1.57)),('gas_2.68',dict(gas=2.68)),
               ('cmult_1.5_commodity',dict(cmult=1.5)),('cmult_1.0_noaudit',dict(cmult=1.0)),('FX_2000_unindexed',dict(price=480*1330/2000)),
               ('staff_local_only_0.3',dict(staff=0.3)),('decl_15',dict(decl=0.15))]:
    p=dict(b);p.update(mod); r=run(1,p); sens[nm]={k:r[k] for k in ['capex','rev1','ebitda1','cash_cost_usd_mmbtu','simple_payback','npv15','npv20','irr']}
out['sens_1_base']=sens
# distance curve for delivered cash+capital cost $/MMBtu (base, 1 MMscf/d): levelized
def lcod(Q,p,r=0.15,years=10):
    x=run(Q,p,years=years)
    crf=r*(1+r)**years/((1+r)**years-1)
    return round((x['capex']*crf+x['opex1'])/x['sold_MMBtu_y_M'],2)
out['levelized_cost_usd_mmbtu_base_1MMscfd']={km:lcod(1,dict(b,km=km)) for km in [30,60,100,150,200,300,400]}
out['levelized_cost_usd_mmbtu_high_1MMscfd']={km:lcod(1,dict(S['high'],km=km)) for km in [30,60,100,150,200,300,400]}
out['levelized_cost_usd_mmbtu_low_1MMscfd']={km:lcod(1,dict(S['low'],km=km)) for km in [30,60,100,150,200,300,400]}

def gate(Q,price_usd,util,cmult,gas,staff,sec,ovh,tor=0.045,ins=0.015,decl=0.06,dev=0.7,years=10,tax=0.34,sell=0.88):
    station=2.0*Q**0.7*cmult; capex=station+dev
    mm=Q*SCM_PER_MMSCF*sell*util*365*MMBTU_PER_SCM/1e6
    cf=[-capex];e1=None
    for t in range(1,years+1):
        g=(1-decl)**(t-1)
        rev=mm*g*price_usd; gasc=Q*SCM_PER_MMSCF*util*g*365*MMBTU_PER_SCM*gas/1e6
        e=rev-gasc-staff-sec-ovh-capex*(tor+ins)
        if e1 is None: e1=e
        cf.append(e-max(0,(e-capex/7)*tax))
    return dict(Q=Q,gate_price=price_usd,capex=round(capex,2),rev1=round(mm*price_usd,2),ebitda1=round(e1,2),npv15=round(npv(0.15,cf),2),npv20=round(npv(0.20,cf),2),irr=(round(irr(cf)*100,1) if irr(cf) is not None else None),simple_payback=(round(capex/e1,1) if e1>0 else None))
out['gate_sale']={}
for Q in [1,5,15]:
    for nm,pr,ut,cm,ga,st,se,ov in [('low',3.5,0.45,2.8,1.0,0.45,0.35,0.15),('base',4.5,0.65,2.6,0.6,0.35,0.25,0.12),('high',6.0,0.8,2.4,0.25,0.30,0.20,0.10)]:
        k=1 if Q==1 else (1.6 if Q==5 else 2.4)
        out['gate_sale'][f'{Q}_{nm}']=gate(Q,pr,ut,cm,ga,st*k,se*k,ov*k)
json.dump(out,open(__file__.replace('cng_model.py','cng_out.json'),'w'),indent=1,ensure_ascii=False)
for k,v in out.items(): print(k,v)
