# Мини-метанол -> формалин -> КФ-смола, Нигерия. Все входы [Д]/[П]/[В] размечены в отчёте.
import math, json
def npv(r, cfs): return sum(c/(1+r)**i for i,c in enumerate(cfs))
def irr(cfs):
    lo,hi=-0.9,1.0
    if npv(lo,cfs)*npv(hi,cfs)>0: return None
    for _ in range(200):
        m=(lo+hi)/2
        if npv(lo,cfs)*npv(m,cfs)<=0: hi=m
        else: lo=m
    return m
def project(capex, ebitda_full, life=12, ramp=(0.6,0.9), decline=0.0, decl_start=4, tax=0.30, dep_y=10):
    cfs=[-capex/2,-capex/2]
    for y in range(1,life+1):
        f = ramp[y-1] if y<=len(ramp) else 1.0
        if y>=decl_start: f*= (1-decline)**(y-decl_start+1)
        e=ebitda_full*f
        dep=capex/dep_y if y<=dep_y else 0
        t=max(0,(e-dep)*tax)
        cfs.append(e-t)
    return cfs
MMBTU_PER_MSCF=(1.03,1.05,1.10)
DAYS=330
def meoh_t(mmscfd, sc):  # sc 0 low,1 base,2 high
    hr=(38,35,32)[sc]   # MMBtu/t малого SMR [Д]; 38 — нижняя граница выхода
    hhv=(1.03,1.05,1.10)[sc]
    return mmscfd*1000*hhv/hr*DAYS
res={}
# ---------- A. только метанол, 1/5/15 MMscf/d
for q in (1,5,15):
    for sc,name in enumerate(('low','base','high')):
        t=meoh_t(q,sc)
        price=(250,330,420)[sc] if q>1 else (300,375,450)[sc]  # нетбэк у скида, $/т [Д]
        rev=t*price
        gas=q*1000*1.05*DAYS*(1.0,0.6,0.25)[sc]
        staff={1:(1.0,0.8,0.6),5:(1.8,1.5,1.2),15:(2.8,2.3,1.9)}[q][sc]*1e6
        capex_tpy={1:(3500,2500,1500),5:(2200,1700,1200),15:(1750,1300,900)}[q][sc]
        capex=t*capex_tpy if sc!=0 else meoh_t(q,0)*capex_tpy
        maint=capex*(0.05,0.045,0.035)[sc]
        cat=t*(20,15,10)[sc]
        ins_sec=capex*(0.03,0.02,0.015)[sc] + (0.3,0.2,0.1)[sc]*1e6*(1 if q==1 else 1.5 if q==5 else 2)
        opex=gas+staff+maint+cat+ins_sec
        e=rev-opex
        dec=(0.08,0.04,0.0)[sc]
        cf=project(capex,e,decline=dec)
        res[f'A_{q}_{name}']=dict(t_y=round(t),price=price,rev=round(rev/1e6,2),opex=round(opex/1e6,2),
            gas=round(gas/1e6,2),staff=staff/1e6,maint=round(maint/1e6,2),ins_sec=round(ins_sec/1e6,2),
            capex=round(capex/1e6,1),ebitda=round(e/1e6,2),
            payback=round(capex/e,1) if e>0 else None,
            npv15=round(npv(0.15,cf)/1e6,1),npv20=round(npv(0.20,cf)/1e6,1),irr=(round(irr(cf)*100,1) if irr(cf) is not None else None))
# ---------- B. интеграция 1 MMscf/d: метанол -> формалин 37% + КФ-смола 65%
MEOH_PER_FORM=0.45; MEOH_PER_UF=0.27; UREA_PER_UF=0.42
for sc,name in enumerate(('low','base','high')):
    t=meoh_t(1,sc)
    div=(2.7,1.9,1.3)[sc]
    uf=33.7e3/div; form=15.3e3/div
    need=uf*MEOH_PER_UF+form*MEOH_PER_FORM
    if need>t:
        k=t/need; uf*=k; form*=k; need=t
    meoh_sold=t-need
    p_uf=(520,580,620)[sc]; p_form=(350,420,500)[sc]; p_me=(300,375,450)[sc]
    rev=uf*p_uf+form*p_form+meoh_sold*p_me
    urea=uf*UREA_PER_UF*(500,480,420)[sc]+uf*UREA_PER_UF*20
    gas=1000*1.05*DAYS*(1.0,0.6,0.25)[sc]
    capex_me=meoh_t(1,0 if sc==0 else sc)*(3500,2500,1500)[sc]
    capex_dn=(11,8,6)[sc]*1e6   # формалин+смола+баки+склад [Д], уже с поправкой аудита
    capex=capex_me+capex_dn
    staff=((1.0,0.8,0.6)[sc]+(0.6,0.45,0.3)[sc])*1e6
    maint=capex*(0.05,0.045,0.035)[sc]
    util_dn=(uf+form)*(15,12,8)[sc]
    freight=(uf+form)*(60,50,40)[sc]   # до клиентов Лагос/Огун/Ондо [Д]
    ins_sec=capex*(0.03,0.02,0.015)[sc]+(0.3,0.2,0.1)[sc]*1e6
    cat=t*(20,15,10)[sc]
    sga=(0.5,0.4,0.3)[sc]*1e6
    opex=urea+gas+staff+maint+util_dn+freight+ins_sec+cat+sga
    e=rev-opex
    cf=project(capex,e,decline=(0.08,0.04,0.0)[sc])
    res[f'B_1_{name}']=dict(meoh_t=round(t),uf_t=round(uf),form_t=round(form),meoh_sold=round(meoh_sold),
        rev=round(rev/1e6,2),urea=round(urea/1e6,2),freight=round(freight/1e6,2),opex=round(opex/1e6,2),capex=round(capex/1e6,1),
        capex_me=round(capex_me/1e6,1),ebitda=round(e/1e6,2),payback=round(capex/e,1) if e>0 else None,
        npv15=round(npv(0.15,cf)/1e6,1),npv20=round(npv(0.20,cf)/1e6,1),irr=(round(irr(cf)*100,1) if irr(cf) is not None else None))
# ---------- C. контроль: только формалин+смола у клиента, метанол покупной (SuperTech/импорт), без факела
for sc,name in enumerate(('low','base','high')):
    div=(2.7,1.9,1.3)[sc]; uf=33.7e3/div; form=15.3e3/div
    if sc==2: uf=25e3; form=9e3
    meoh=uf*MEOH_PER_UF+form*MEOH_PER_FORM
    p_me_buy=(560,480,420)[sc]
    p_uf=(520,580,620)[sc]; p_form=(350,420,500)[sc]
    rev=uf*p_uf+form*p_form
    urea=uf*UREA_PER_UF*((500,480,420)[sc]+20)
    capex=(11,8,6)[sc]*1e6
    opex=urea+meoh*p_me_buy+(0.6,0.45,0.3)[sc]*1e6+capex*(0.05,0.045,0.035)[sc]+(uf+form)*((15,12,8)[sc]+(25,20,15)[sc])+capex*(0.02,0.015,0.01)[sc]+(0.5,0.4,0.3)[sc]*1e6
    e=rev-opex
    cf=project(capex,e)
    res[f'C_{name}']=dict(uf_t=round(uf),form_t=round(form),meoh_buy_t=round(meoh),rev=round(rev/1e6,2),opex=round(opex/1e6,2),capex=capex/1e6,ebitda=round(e/1e6,2),
        payback=round(capex/e,1) if e>0 else None,npv15=round(npv(0.15,cf)/1e6,1),npv20=round(npv(0.20,cf)/1e6,1),irr=(round(irr(cf)*100,1) if irr(cf) is not None else None))
# ---------- D. путь аудита корпуса: метанол-скид econ_vs_mining (capex 3.0, opex 0.4, EBITDA 1.94)
for k,(cx,ox) in {'audit_low':(2.8,2.4),'audit_high':(2.4,1.7)}.items():
    capex=3.0*cx; gross=2.34; opex=0.4*ox
    res['D_'+k]=dict(capex=capex,opex=round(opex,2),ebitda_screen_audit=round(gross-opex,2))
json.dump(res,open('/tmp/claude-0/-home-user-ocm-research-diary/50ead64b-47e2-5bdd-b34b-09f6667b6ca9/scratchpad/meoh/meoh_ng_results.json','w'),ensure_ascii=False,indent=1)
for k,v in res.items(): print(k,v)
