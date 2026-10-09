# Скептик: второй год, деньги, простой, инфляция для цепочки метанол->формалин->КФ (B) и контроля C
import json
def npv(r,c): return sum(x/(1+r)**i for i,x in enumerate(c))
def irr(c):
    lo,hi=-0.9,1.0
    if npv(lo,c)*npv(hi,c)>0: return None
    for _ in range(200):
        m=(lo+hi)/2
        if npv(lo,c)*npv(m,c)<=0: hi=m
        else: lo=m
    return m
MEF,MEU,UREA=0.45,0.27,0.42
def run(kind='B',uf=17737,form=8053,p_uf=580,p_form=420,p_me=375,meoh_buy=480,urea_p=480,
        capex_me=24.8e6,capex_dn=8e6,gas_av=1.0,ramp=(0.6,0.9),decl=0.04,wc=0.0,tax=0.30,carry=False,
        dev_shock=None,usd_infl=0.0,debt=None,life=12):
    t_me=9900 if kind=='B' else 0
    capex=(capex_me if kind=='B' else 0)+capex_dn
    cf=[-capex/2,-capex/2]; rows=[]; loss=0; prev_rev=0
    for y in range(1,life+1):
        f=ramp[y-1] if y<=len(ramp) else 1.0
        g=(1-decl)**(y-3) if y>=4 else 1.0
        sales=f   # сбыт продукта
        uf_y=uf*sales; fo_y=form*sales
        need=uf_y*MEU+fo_y*MEF
        own=min(t_me*gas_av*g, t_me) if kind=='B' else 0
        own_used=min(own,need); sold_me=max(0,own-need); bought=max(0,need-own)
        fx=1.0
        if dev_shock and y>=dev_shock[0]:
            # девальвация d, перенос в цену смолы/формалина с лагом 1 год (импортный паритет)
            d=dev_shock[1]; fx=1/(1+d) if y==dev_shock[0] else 1.0
        rev_ngn=uf_y*p_uf+fo_y*p_form
        rev=rev_ngn*fx+sold_me*p_me
        infl=(1+usd_infl)**y
        urea=uf_y*UREA*(urea_p+20)*fx
        gas=(own/9900)*1000*1.05*330*0.6 if kind=='B' else 0
        staff=(1.25e6 if kind=='B' else 0.45e6)
        staff=staff*(0.75*fx+0.25)*infl  # 25% экспат в USD
        maint=capex*0.045*infl
        util=(uf_y+fo_y)*12*fx
        fr=(uf_y+fo_y)*(50 if kind=='B' else 20)*fx
        ins=(capex*0.02+0.2e6 if kind=='B' else capex*0.015)*infl
        cat=own*15
        sga=0.4e6*infl
        buy=bought*meoh_buy
        e=rev-(urea+gas+staff+maint+util+fr+ins+cat+sga+buy)
        dep=capex/10 if y<=10 else 0
        ti=e-dep
        if carry:
            ti2=ti-loss
            if ti2<0: loss=-ti2; tx=0
            else: loss=0; tx=ti2*tax
        else: tx=max(0,ti*tax)
        dwc=wc*(rev-prev_rev); prev_rev=rev
        c=e-tx-dwc
        if y==life: c+=wc*rev
        cf.append(c); rows.append(dict(y=y,rev=round(rev/1e6,2),ebitda=round(e/1e6,2),cfads=round((e-tx)/1e6,2),dwc=round(dwc/1e6,2)))
    r=dict(capex=round(capex/1e6,1),y1=rows[0]['ebitda'],y2=rows[1]['ebitda'],y3=rows[2]['ebitda'],y2_cf=round((rows[1]['cfads']-rows[1]['dwc']),2),
           npv15=round(npv(.15,cf)/1e6,1),npv20=round(npv(.20,cf)/1e6,1),irr=(round(irr(cf)*100,1) if irr(cf) is not None else None))
    # срок возврата по фактическому потоку от FID
    cum=0; pb=None
    for i,x in enumerate(cf):
        cum+=x
        if cum>=0 and pb is None: pb=i
    r['payback_from_FID_y']=pb
    if debt:
        sh,rate,ten,grace=debt; D=capex*sh
        ann=D*rate/(1-(1+rate)**-ten)
        ds=[D*rate if y<=grace else ann for y in range(1,4)]
        r['ds_musd']=round(ann/1e6,2)
        r['dscr_y1_3']=[round(rows[i]['cfads']*1e6/ds[i],2) for i in range(3)]
    return r
out={}
D=(0.6,0.12,8,1)  # 60% долг, 12% USD, 8 лет, 1 год льготы по телу [Д]
out['B_as_assessed']=run('B')
out['B_as_assessed_real']=run('B',tax=0.34,carry=True,wc=0.2,debt=D)
# поправки: цена смолы = импортный паритет порошка в жидком эквиваленте + пошлина/логистика
out['B_p450']=run('B',p_uf=450,tax=0.34,carry=True,wc=0.2,debt=D)
out['B_p450_gas82']=run('B',p_uf=450,tax=0.34,carry=True,wc=0.2,gas_av=0.82,debt=D)
out['B_p450_gas82_ramp47']=run('B',p_uf=450,tax=0.34,carry=True,wc=0.2,gas_av=0.82,ramp=(0.4,0.7),debt=D)
out['B_corr_capex']=run('B',p_uf=450,tax=0.34,carry=True,wc=0.2,gas_av=0.82,ramp=(0.4,0.7),capex_me=24.8e6*1.3,capex_dn=10e6,debt=D)
out['B_corr_capex_8kt']=run('B',uf=8000,p_uf=450,tax=0.34,carry=True,wc=0.2,gas_av=0.82,ramp=(0.4,0.7),capex_me=24.8e6*1.3,capex_dn=9e6,debt=D)
out['B_upside_580_but_dev30_y2']=run('B',tax=0.34,carry=True,wc=0.2,dev_shock=(2,0.30),debt=D)
out['C_as_assessed']=run('C')
out['C_real_580']=run('C',tax=0.34,carry=True,wc=0.2,debt=(0.6,0.12,6,1))
out['C_p450']=run('C',p_uf=450,tax=0.34,carry=True,wc=0.2,debt=(0.6,0.12,6,1))
out['C_p450_ramp47_dn10']=run('C',p_uf=450,tax=0.34,carry=True,wc=0.2,ramp=(0.4,0.7),capex_dn=10e6,debt=(0.6,0.12,6,1))
out['C_p450_8kt']=run('C',uf=8000,p_uf=450,tax=0.34,carry=True,wc=0.2,ramp=(0.4,0.7),capex_dn=9e6)
out['C_p500_8kt']=run('C',uf=8000,p_uf=500,tax=0.34,carry=True,wc=0.2,ramp=(0.4,0.7),capex_dn=9e6)
# B минус C приростно на поправленной базе
b=run('B',p_uf=450,gas_av=0.82); c=run('C',p_uf=450)
out['incr_B_minus_C_p450']=dict(npv15=round(b['npv15']-c['npv15'],1))
# нетто-валютная экспозиция базы (оценка автора: -$2.8 млн при -20%)
rev_ngn=17737*580+8053*420; costs_ngn=17737*0.42*500+ (17737+8053)*(50+12) + 1.25e6*0.75
for d in (0.2,0.3):
    out[f'fx_net_loss_dev{int(d*100)}']=round((rev_ngn-costs_ngn)*(1-1/(1+d))/1e6,2)
# безубыточная цена смолы по NPV15=0 для C (реальная модель)
for p in range(400,900,5):
    if run('C',p_uf=p,tax=0.34,carry=True,wc=0.2,ramp=(0.4,0.7),capex_dn=10e6)['npv15']>=0: out['C_be_price_npv15']=p;break
for p in range(400,1200,10):
    if run('B',p_uf=p,tax=0.34,carry=True,wc=0.2,gas_av=0.82,ramp=(0.4,0.7),capex_me=24.8e6*1.3,capex_dn=10e6)['npv15']>=0: out['B_be_price_npv15']=p;break
json.dump(out,open('/tmp/claude-0/-home-user-ocm-research-diary/50ead64b-47e2-5bdd-b34b-09f6667b6ca9/scratchpad/skeptic_meoh_y2/y2_out.json','w'),ensure_ascii=False,indent=1)
for k,v in out.items(): print(k,v)
