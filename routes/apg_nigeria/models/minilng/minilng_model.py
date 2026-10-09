# Мини-СПГ на факельном газе, Нигерия: скрин vs поправки аудита NaCN. $M, $/MMBtu.
import json
ESC = 800/576.1          # CEPCI 2014 -> 2025 [Д]
HHV_LNG = 52.6           # MMBtu/t (55.5 MJ/kg HHV) [В справочно]
FEED_HHV = 1.05          # MMBtu/Mscf попутного газа [Д 1.03-1.1, база контекста]
def npv(cfs,r): return sum(c/(1+r)**i for i,c in enumerate(cfs))
def irr(cfs):
    lo,hi=-0.99,3.0
    if npv(cfs,lo)*npv(cfs,hi)>0: return None
    for _ in range(200):
        m=(lo+hi)/2
        if npv(cfs,lo)*npv(cfs,m)<=0: hi=m
        else: lo=m
    return m
def capex_parts(Q):
    eq = (1.5607*Q+3.1243)*ESC                 # WB 2016 Fig.31, оборудование [П] -> 2025
    inst = (3.1214*Q+6.2485)*ESC               # WB Fig.32, установленная стоимость [П]
    pretreat = 19*(Q/15)**0.6*ESC              # WB: доп. подготовка попутного газа ~19 MMUSD @15 [П], масштаб 0.6 [Д]
    stor_m3 = 3*15.5*Q/0.43                    # 3 сут хранения, 15.5 т/сут на MMscf/d, 0.43 т/м3
    stor = stor_m3*2000/1e6                    # WB $2000/м3 bullets [П]
    site = 1.0+0.1*Q                           # подключение к факелу, площадка, налив [Д]
    return dict(eq=eq,inst=inst,pretreat=pretreat,stor=stor,site=site)
def capex(Q,case):
    p=capex_parts(Q)
    if case=='screen': return 5.0*Q            # линейно от Аджаокуты ~$5 млн/MMscf/д [В->О]
    if case=='audit_low': return p['eq']*2.4   # пакет оборудования x2.4 (правило аудита)
    if case=='audit_high': return (p['inst']+p['pretreat']+p['stor']+p['site']+0.2*Q)*1.2  # снизу вверх, x1.2 Нигерия [Д]
    if case=='audit_mid': return 0.5*(capex(Q,'audit_low')+capex(Q,'audit_high'))
def lng_mmbtu_per_day(Q, kwh_kg, ngl_share, loss):
    feed = FEED_HHV*1000*Q
    fuel_ratio = kwh_kg/0.38*3.6/(HHV_LNG*1.055)   # МДж топлива на МДж СПГ, ГПУ КПД 0.38
    return (feed*(1-ngl_share))*(1-loss)/(1+fuel_ratio)
CASES = {
 # netback = цена ex-plant, $/MMBtu СПГ
 'screen':     dict(kwh=0.7, ngl_sh=0.05, loss=0.02, ngl_t=2.0, ngl_p=600, LF=0.92, decl=0.00, gas=0.25, nb=10.0,
                    staff=lambda Q:0.25+0.03*Q, sec=lambda Q:0.1, mnt=0.03, ins=0.005, cons=0.20, misc=lambda Q:0.1, build=[0.5,0.5]),
 'audit':      dict(kwh=1.1, ngl_sh=0.05, loss=0.03, ngl_t=1.0, ngl_p=450, LF=0.80, decl=0.08, gas=1.00, nb=6.0,
                    staff=lambda Q:0.5+0.09*Q, sec=lambda Q:0.4+0.04*Q, mnt=0.045, ins=0.015, cons=0.45, misc=lambda Q:0.2+0.02*Q, build=[0.4,0.6]),
}
def run(Q, case='audit', ccase='audit_mid', life=12, tax=0.34, **kw):
    S=dict(CASES[case]); S.update(kw)
    C = capex(Q, ccase if case=='audit' else 'screen') if 'C' not in kw else kw['C']
    lng_d = lng_mmbtu_per_day(Q,S['kwh'],S['ngl_sh'],S['loss'])
    cfs=[-C*b for b in S['build']]; dep=C/5; rows=[]
    for yr in range(life):
        f=S['LF']*(1-S['decl'])**yr
        lng=lng_d*365*f; ngl=S['ngl_t']*Q*365*f
        rev=(lng*S['nb']+ngl*S['ngl_p'])/1e6
        gas=S['gas']*FEED_HHV*1000*Q*365*f/1e6
        fixed=S['staff'](Q)+S['sec'](Q)+S['misc'](Q)+C*(S['mnt']+S['ins'])
        var=gas+S['cons']*lng/1e6
        e=rev-fixed-var
        t=max(0,(e-(dep if yr<5 else 0))*tax); ed=min(t,0.05*C) if yr<5 else 0
        cfs.append(e-t+ed); rows.append(dict(yr=yr+1,lng_t=round(lng/HHV_LNG),rev=round(rev,2),opex=round(fixed+var,2),ebitda=round(e,2),cf=round(e-t+ed,2)))
    cum=-C; pb=None
    for i,r in enumerate(rows):
        if cum<0 and cum+r['cf']>=0: pb=i+(-cum)/r['cf']
        cum+=r['cf']
    r=irr(cfs)
    return dict(Q=Q,case=case,ccase=ccase,capex=round(C,1),lng_tpd_nameplate=round(lng_d/HHV_LNG,1),lng_kt_y1=round(rows[0]['lng_t']/1e3,2),
                rev1=rows[0]['rev'],opex1=rows[0]['opex'],ebitda1=rows[0]['ebitda'],NPV15=round(npv(cfs,.15),1),NPV20=round(npv(cfs,.20),1),
                IRR=None if r is None else round(r*100,1),payback=None if pb is None else round(pb,1))
def breakeven_nb(Q, r=0.15, ccase='audit_mid', **kw):
    lo,hi=0.0,40.0
    for _ in range(60):
        m=(lo+hi)/2
        x=run(Q,'audit',ccase,nb=m,**kw)
        if (x['NPV15'] if r==0.15 else x['NPV20'])<0: lo=m
        else: hi=m
    return round(hi,2)
out={}
print('== capex parts ($M 2025) =='); 
for Q in [1,3,5,10,15]:
    p=capex_parts(Q); print(Q,{k:round(v,1) for k,v in p.items()},'screen',capex(Q,'screen'),'aud_low',round(capex(Q,'audit_low'),1),'aud_high',round(capex(Q,'audit_high'),1),
      'eq x2.8',round(p['eq']*2.8,1), '| $/tpa mid', round(capex(Q,'audit_mid')*1e6/(lng_mmbtu_per_day(Q,1.1,.05,.03)/HHV_LNG*365*0.9)))
    out[f'capex_{Q}']=dict(screen=capex(Q,'screen'),audit_low=round(capex(Q,'audit_low'),1),audit_high=round(capex(Q,'audit_high'),1),audit_mid=round(capex(Q,'audit_mid'),1))
print('LNG t/d per MMscf/d: screen',round(lng_mmbtu_per_day(1,0.7,.05,.02)/HHV_LNG,1),'audit',round(lng_mmbtu_per_day(1,1.1,.05,.03)/HHV_LNG,1))
print('\n== site economics ==')
for Q in [1,5,15]:
    s=run(Q,'screen'); print(s); out[f'screen_{Q}']=s
    for cc in ['audit_low','audit_mid','audit_high']:
        a=run(Q,'audit',cc); print(' ',a); out[f'{cc}_{Q}']=a
    for nb in [4.5,8.0,10.0,12.0]:
        a=run(Q,'audit','audit_mid',nb=nb); print('   nb',nb,a); out[f'audit_mid_{Q}_nb{nb}']=a
    be15=breakeven_nb(Q,0.15); be20=breakeven_nb(Q,0.20)
    be15l=breakeven_nb(Q,0.15,'audit_low'); be15h=breakeven_nb(Q,0.15,'audit_high')
    print('   breakeven netback NPV15=0 mid/low/high capex:',be15,be15l,be15h,' NPV20=0 mid:',be20)
    out[f'breakeven_{Q}']=dict(nb15_mid=be15,nb15_lowcapex=be15l,nb15_highcapex=be15h,nb20_mid=be20)
print('\n== sensitivities Q=5 audit_mid ==')
for lab,kw in [('gas0.25',dict(gas=0.25)),('LF0.5_resource',dict(LF=0.5)),('decl0',dict(decl=0.0)),('decl15',dict(decl=0.15)),('nb6_life20_decl0',dict(life=20,decl=0.0)),
               ('ngl_rich_4t_600',dict(ngl_t=4.0,ngl_p=600,ngl_sh=0.12)),('kwh0.8',dict(kwh=0.8)),('tax_holiday_0',dict(tax=0.0))]:
    for Q in [5,15]:
        a=run(Q,'audit','audit_mid',**kw); print(Q,lab,a); out[f'sens_{Q}_{lab}']=a
# screen opex vs audit opex ratio
for Q in [1,5,15]:
    s=run(Q,'screen'); a=run(Q,'audit','audit_mid')
    # same-capex opex comparison to isolate cash-cost multiplier
    a2=run(Q,'audit','audit_mid',C=s['capex'])
    print('Q',Q,'screen opex',s['opex1'],'audit opex (same capex)',a2['opex1'],'ratio',round(a2['opex1']/s['opex1'],2),'| audit opex own capex',a['opex1'])
    out[f'opex_ratio_{Q}']=round(a2['opex1']/s['opex1'],2)
# cost of liquefaction chain per MMBtu (levelized, 15%, audit_mid) = breakeven netback minus nothing
json.dump(out,open('/tmp/claude-0/-home-user-ocm-research-diary/50ead64b-47e2-5bdd-b34b-09f6667b6ca9/scratchpad/minilng/minilng_out.json','w'),ensure_ascii=False,indent=1)
