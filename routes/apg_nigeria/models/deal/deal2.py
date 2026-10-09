import sys,io,contextlib,json
sys.path.insert(0,'/tmp/claude-0/-home-user-ocm-research-diary/50ead64b-47e2-5bdd-b34b-09f6667b6ca9/scratchpad/deal')
with contextlib.redirect_stdout(io.StringIO()):
    import dealcore as D
M=D.M; run=D.run; P=('P',)
DEAL=dict(dep=0.10,idx=1.0,lag=0,usd_inf=0.03,up=11/12,ramp=0.8,wc=60,delay=1)
STRUCT=dict(DEAL,up=1.0,wc=30)          # DoP оператора + оплата 30 дней через эскроу
def be(kw,key='npv15',ov=None,**ex):
    lo,hi=150,900
    for _ in range(40):
        m=(lo+hi)/2; o=dict(p_dir=m); o.update(ov or {})
        if run(5,P,ov=o,**ex,**kw)[key]<0: lo=m
        else: hi=m
    return round(hi)
out={}
for tag,kw,ov in [('DEAL',DEAL,{}),('STRUCT bad5%',STRUCT,{}),('STRUCT escrow bad2%',STRUCT,dict(bad=0.02)),('DEAL escrow bad2% noDoP',dict(DEAL,wc=30),dict(bad=0.02))]:
    a,b=be(kw,'npv15',ov),be(kw,'npv20',ov); out[tag]=(a,b); print(tag,'BE NPV15',a,'NPV20',b)
    for pr in (275,300,330,350,380):
        o=dict(p_dir=pr); o.update(ov); r=run(5,P,ov=o,debt=(0.6,0.14,7),**kw)
        print('   ₦',pr,{k:r[k] for k in ('y1','y2','npv15','npv20','irr','payback','min_dscr')})
# ёмкость earn-out держателю ($/Mscf на оплачиваемый газ 2,47 MMscf/д) при NPV15=0 и NPV20=0
S2=STRUCT; ov0=dict(bad=0.02)
for pr in (330,350,380):
    for key in ('npv15','npv20'):
        lo,hi=0.0,5.0
        o=dict(p_dir=pr,opfee=0.15); o.update(ov0)
        if run(5,P,ov=o,**S2)[key]<0: print('earnout cap',pr,key,'<0 already'); continue
        for _ in range(40):
            m=(lo+hi)/2; o=dict(p_dir=pr,opfee=0.15+m); o.update(ov0)
            if run(5,P,ov=o,**S2)[key]<0: hi=m
            else: lo=m
        print('earnout cap ₦',pr,key,round(lo,2),'$/Mscf ->', round(lo*2.47*365,2),'$M/yr (на 2,47 MMscf/д)')
# плата оператора за снятие факела x0.6
for fee in (0.6,1.05):
    print('operator fee',fee,'BE NPV15',be(STRUCT,'npv15',dict(bad=0.02),extra_fee=fee),'NPV20',be(STRUCT,'npv20',dict(bad=0.02),extra_fee=fee))
# газ: 0.5 и 0.6 от паспорта при STRUCT ₦350
for av in (0.85,0.6,0.5):
    r=run(5,P,ov=dict(p_dir=350,bad=0.02,avail=av),debt=(0.6,0.14,7),**STRUCT); print('avail',av,{k:r[k] for k in ('y1','y2','y5','npv15','npv20','irr','min_dscr')})
# падение дебита 10%/г вместо 6%
r=run(5,P,ov=dict(p_dir=350,bad=0.02,decl=0.10),debt=(0.6,0.14,7),**STRUCT); print('decl10',{k:r[k] for k in ('y1','y2','y5','npv15','npv20','irr','min_dscr')})
# капитал x2.8 при STRUCT ₦350
r=run(5,P,ov=dict(p_dir=350,bad=0.02,cm=2.8),debt=(0.6,0.14,7),**STRUCT); print('cm2.8',{k:r[k] for k in ('capex','y1','npv15','npv20','irr','min_dscr')})
# задержка 2 г
r=run(5,P,ov=dict(p_dir=350,bad=0.02),debt=(0.6,0.14,7),**dict(STRUCT,delay=2)); print('delay2',{k:r[k] for k in ('y1','npv15','npv20','irr','min_dscr')})
# ренегоциация -20% со 2-го года (дизель падает)
r=run(5,P,ov=dict(p_dir=350,bad=0.02),debt=(0.6,0.14,7),tariff_cut=(2,0.8),**STRUCT); print('cut20',{k:r[k] for k in ('y1','y2','npv15','npv20','irr','min_dscr')})
# блокада 3 мес на 2-й год
r=run(5,P,ov=dict(p_dir=350,bad=0.02),debt=(0.6,0.14,7),**dict(STRUCT,ev={2:0.75})); print('block3m',{k:r[k] for k in ('y1','y2','npv15','npv20','irr','min_dscr')})
# финансирование: debt service
C=29.1; Dd=C*0.6; ann=Dd*0.14/(1-1.14**-7); print('debt',round(Dd,2),'annuity',round(ann,2),'DSRA6m',round(ann/2,2),'IDC1y(half drawn)',round(Dd*0.14*0.5,2))
r=M.site(5,'base',P,verbose=True,ov=dict(p_dir=350,bad=0.02)); print('rev y1 (no frictions) at 350',r['y1'])
