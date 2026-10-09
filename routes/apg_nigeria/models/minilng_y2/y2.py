# Скептик "2-й год и деньги" для мини-СПГ. Переиспользует формулы minilng_model.py, но с годовыми профилями.
import importlib.util,io,contextlib,json
P='/tmp/claude-0/-home-user-ocm-research-diary/50ead64b-47e2-5bdd-b34b-09f6667b6ca9/scratchpad/minilng/minilng_model.py'
src=open(P).read().split("out={}")[0]          # только функции, без перезаписи out.json
m={}; exec(src,m)
HHV=m['HHV_LNG']; FEED=m['FEED_HHV']
def npv(c,r): return sum(x/(1+r)**i for i,x in enumerate(c))
def irr(c):
    lo,hi=-0.9,3.0
    if npv(c,lo)*npv(c,hi)>0: return None
    for _ in range(200):
        mid=(lo+hi)/2
        if npv(c,lo)*npv(c,mid)<=0: hi=mid
        else: lo=mid
    return mid
def sim(Q=15,C=None,nb=10.0,LF=0.9,decl=0.0,life=15,gas=0.25,kwh=1.1,tax=0.34,
        lf_prof=None,nb_prof=None,cost_infl=0.0,wc_days=0,extra_capex=0.0,out_fixed_share=1.0):
    S=dict(m['CASES']['audit'])
    C=(m['capex'](Q,'audit_mid') if C is None else C)+extra_capex
    lng_d=m['lng_mmbtu_per_day'](Q,kwh,S['ngl_sh'],S['loss'])
    cfs=[-C*0.4,-C*0.6]; dep=C/5; rows=[]; wc_prev=0
    for yr in range(life):
        f=LF*(1-decl)**yr
        if lf_prof and yr in lf_prof: f=lf_prof[yr]
        p=nb_prof.get(yr,nb) if nb_prof else nb
        lng=lng_d*365*f; ngl=S['ngl_t']*Q*365*f
        rev=(lng*p+ngl*S['ngl_p'])/1e6
        g=gas*FEED*1000*Q*365*f/1e6
        infl=(1+cost_infl)**yr
        fixed=(S['staff'](Q)+S['sec'](Q)+S['misc'](Q)+C*(S['mnt']+S['ins']))*infl
        var=g+S['cons']*lng/1e6*infl
        e=rev-fixed-var
        t=max(0,(e-(dep if yr<5 else 0))*tax); ed=min(t,0.05*C) if yr<5 else 0
        wc=rev*wc_days/365; dwc=wc-wc_prev; wc_prev=wc
        cf=e-t+ed-dwc
        if yr==life-1: cf+=wc
        cfs.append(cf); rows.append(dict(yr=yr+1,rev=round(rev,1),opex=round(fixed+var,1),ebitda=round(e,1),cf=round(cf,1)))
    r=irr(cfs)
    return dict(C=round(C,1),NPV15=round(npv(cfs,.15),1),NPV20=round(npv(cfs,.20),1),IRR=None if r is None else round(r*100,1),y1=rows[0],y2=rows[1],rows=rows,cfs=cfs)
out={}
def show(k,d):
    out[k]={x:d[x] for x in ('C','NPV15','NPV20','IRR','y1','y2')}; print(k,out[k])
# 0. воспроизвести "благоприятный" случай
base=sim(); show('fav_nb10_repro',base)
for nb in (8,12): show(f'fav_nb{nb}',sim(nb=nb))
# 1. энергоостров: средн. 10.6 МВт, установлено x1.25 (N+1), $1000-1500/кВт установл. [Д]
mw=1.1*m['lng_mmbtu_per_day'](15,1.1,.05,.03)/HHV/24
pi_lo=mw*1.25*1000/1e3; pi_hi=mw*1.25*1500/1e3
print('avgMW',round(mw,2),'power island $M',round(pi_lo,1),round(pi_hi,1))
show('fav_+power_lo',sim(extra_capex=pi_lo)); show('fav_+power_hi',sim(extra_capex=pi_hi))
# 2. второй год: простой 6 мес. (ремонт компрессора холодильного цикла / блокада общины / shut-in оператора)
show('fav_y2_outage6m',sim(lf_prof={1:0.45}))
show('fav_y2_outage3m',sim(lf_prof={1:0.675}))
# 3. ценовая война с года 2 (Ajaokuta/Greenville): нетбэк $8 и $7 навсегда
show('fav_pricewar_y2_8',sim(nb_prof={i:8.0 for i in range(1,15)}))
show('fav_pricewar_y2_7',sim(nb_prof={i:7.0 for i in range(1,15)}))
# 4. девальвация 30% в год 2, цена в найре индексируется с лагом: -23% в г2, -12% в г3, далее 0
show('fav_deval30_lag',sim(nb_prof={1:7.7,2:8.85}))
show('fav_deval30_noindex',sim(nb_prof={i:7.7 for i in range(1,15)}))
# 5. инфляция USD-затрат (ЗИП, экспат, страховка) 4%/г при цене без индексации
show('fav_costinfl4',sim(cost_infl=0.04))
# 6. оборотный капитал: дебиторка 60 сут у промпотребителя
show('fav_wc60',sim(wc_days=60))
# 7. комбинированный "реалистичный 2-й год": энергоостров середина + простой 3 мес. в г2 + WC 60 сут + нетбэк 9 с г2
comb=sim(extra_capex=(pi_lo+pi_hi)/2,lf_prof={1:0.675},wc_days=60,nb_prof={i:9.0 for i in range(1,15)})
show('fav_combined',comb)
# 8. ФИНАНСИРОВАНИЕ: 60% долга, 2 года grace (стройка), аннуитет 7 лет; ставка USD DFI 12% / naira-эквивалент 25%
def dscr(res,C,rate,share=0.6,n=7):
    D=C*share*(1+rate*0.6)  # капитализация процентов на стройке ~ (упрощ.)
    ann=D*rate/(1-(1+rate)**-n)
    cf=[x for x in res['cfs'][2:]]
    return round(D,1),round(ann,1),[round(cf[i]/ann,2) for i in range(3)]
for lab,res in [('fav',base),('fav_+power_mid',sim(extra_capex=(pi_lo+pi_hi)/2)),('combined',comb)]:
    for rate in (0.12,0.25):
        D,ann,ds=dscr(res,res['C'],rate); out[f'dscr_{lab}_{rate}']=dict(debt=D,annuity=ann,DSCR_y1_3=ds); print('DSCR',lab,rate,D,ann,ds)
# 9. конкурент: Аджаокута, трубный газ $2.68, $5.2 млн/MMscf/д, поезд ~19.4 MMscf/д, LF 0.85, без доп.подготовки
def be(**kw):
    lo,hi=0,40
    for _ in range(60):
        mid=(lo+hi)/2
        if sim(nb=mid,**kw)['NPV15']<0: lo=mid
        else: hi=mid
    return round(hi,2)
aj=be(Q=19.4,C=5.2*19.4,gas=2.68,LF=0.85,life=15)
aj_cash=round(2.68*1.05*19.4*1000/ (m['lng_mmbtu_per_day'](19.4,1.1,.05,.03)) + 0.45 + 0,2)
print('Ajaokuta full-cost netback NPV15=0',aj,'gas+consumables per MMBtu LNG',aj_cash)
out['ajaokuta_fullcost_nb']=aj; out['ajaokuta_gas_cons_per_mmbtu']=aj_cash
out['flare_fav_breakeven_nb']=be(); out['flare_fav_power_breakeven_nb']=be(extra_capex=(pi_lo+pi_hi)/2)
print('flare fav breakeven',out['flare_fav_breakeven_nb'],'with power island',out['flare_fav_power_breakeven_nb'])
json.dump(out,open('/tmp/claude-0/-home-user-ocm-research-diary/50ead64b-47e2-5bdd-b34b-09f6667b6ca9/scratchpad/minilng_y2/y2_out.json','w'),ensure_ascii=False,indent=1)
print('---- extra')
pm=(pi_lo+pi_hi)/2
for cc in ('audit_low','audit_high'):
    C=m['capex'](15,cc)
    for lab,kw in [('audit',dict(LF=0.8,decl=0.08,life=12,gas=1.0,nb=6.0))]:
        a=sim(C=C,extra_capex=(pi_lo if cc=='audit_low' else pi_hi),**kw); print(cc,'+power',a['C'],a['NPV15'],a['NPV20'],a['y1'])
print('fav+power_mid+rampup y1 LF0.6', {k:v for k,v in sim(extra_capex=pm,lf_prof={0:0.6}).items() if k in ('C','NPV15','NPV20','IRR')})
print('fav+power_mid nb12', {k:v for k,v in sim(extra_capex=pm,nb=12).items() if k in ('C','NPV15','NPV20','IRR')})
print('fav+power_mid nb11', {k:v for k,v in sim(extra_capex=pm,nb=11).items() if k in ('C','NPV15','NPV20','IRR')})
r=sim(extra_capex=pm); print('fav+power_mid', {k:r[k] for k in ('C','NPV15','NPV20','IRR')}, 'payback', round(r['C']/r['y1']['cf'],1))
def be20(**kw):
    lo,hi=0,40
    for _ in range(60):
        mid=(lo+hi)/2
        if sim(nb=mid,**kw)['NPV20']<0: lo=mid
        else: hi=mid
    return round(hi,2)
print('NPV20 breakeven fav',be20(),'fav+power',be20(extra_capex=pm))
