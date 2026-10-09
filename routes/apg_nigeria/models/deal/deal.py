# Узкий вариант «энергия якорю, 5 MMscf/д»: условия сделки. Модель: копия multi_model.py + y2.run (скептик 3, train).
import sys,io,contextlib,json
sys.path.insert(0,'/tmp/claude-0/-home-user-ocm-research-diary/50ead64b-47e2-5bdd-b34b-09f6667b6ca9/scratchpad/deal')
with contextlib.redirect_stdout(io.StringIO()):
    import dealcore as D
M=D.M; run=D.run; P=('P',)
DEAL=dict(dep=0.10,idx=1.0,lag=0,usd_inf=0.03,up=11/12,ramp=0.8,wc=60,delay=1)       # USD-тариф, B1-трения, стройка 1 г
DOP=dict(DEAL,up=1.0)                                                                  # + deliver-or-pay оператора
BAD=dict(DEAL,ev={2:0.75},ramp=0.7,wc=75,delay=2)                                       # + блокада 3 мес, рамп 70, стройка 2 г
NODELAY=dict(dep=0.10,idx=1.0,lag=0)                                                    # USD-тариф, без трений (проверка ₦250/272)
out={}
def show(tag,r): out[tag]=r; print(tag, {k:r.get(k) for k in ('capex','y1','y2','y5','npv15','npv20','irr','payback','min_dscr','dscr_y1_3')})
print('== репродукция: оценка без трений, ₦320, без индексации/константный курс'); show('R0',run(5,P))
for nm,kw in [('NODELAY',NODELAY),('DEAL',DEAL),('DEAL+DoP',DOP),('BAD',BAD)]:
    for pr in (250,275,300,320,350):
        show(f'{nm} ₦{pr}',run(5,P,ov=dict(p_dir=pr),debt=(0.6,0.14,7),**kw))
def be(kw,key='npv15',lo=150,hi=900,**extra):
    for _ in range(40):
        m=(lo+hi)/2
        if run(5,P,ov=dict(p_dir=m,**extra.get('ov',{})),**{k:v for k,v in extra.items() if k!='ov'},**kw)[key]<0: lo=m
        else: hi=m
    return round(hi)
print('== безубыточный тариф, ₦/кВт·ч (USD-эквивалент при ₦1330)')
for nm,kw in [('NODELAY',NODELAY),('DEAL',DEAL),('DEAL+DoP',DOP),('BAD',BAD)]:
    a,b=be(kw,'npv15'),be(kw,'npv20'); out[f'BE {nm}']=(a,b); print(nm,'NPV15=0 при',a,'| NPV20=0 при',b, '| $/кВт·ч', round(a/1330,3), round(b/1330,3))
print('== чувствительности DEAL ₦300')
base=dict(p_dir=300)
for tag,ov,extra in [('cm2.8',dict(cm=2.8),{}),('cm2.4',dict(cm=2.4),{}),('avail0.5',dict(avail=0.5),{}),('avail0.4',dict(avail=0.4),{}),('avail0.25',dict(avail=0.25),{}),
                     ('LF0.6',dict(LF=0.6),{}),('gas DBP2.18',dict(gas=2.18,opfee=0.0),{}),
                     ('earnout +0.25/Mscf',dict(opfee=0.40),{}),('earnout +0.50/Mscf',dict(opfee=0.65),{}),('earnout +1.00/Mscf',dict(opfee=1.15),{}),
                     ('gasToP whole site +0.51',{},dict(gas_extra=0.51)),
                     ('operator fee 1.0 x0.6',{},dict(extra_fee=0.6)),('operator fee 1.75 x0.6',{},dict(extra_fee=1.05))]:
    o=dict(base); o.update(ov)
    show('S ₦300 '+tag,run(5,P,ov=o,debt=(0.6,0.14,7),**extra,**DEAL))
print('== якорь меньше: 6 МВт и 8,5 МВт (станция под спрос)')
for mw in (6.0,8.5,11.0):
    M.DIRECT_MW[5]=mw
    for pr in (300,320,350):
        show(f'anchor{mw} ₦{pr}',run(5,P,ov=dict(p_dir=pr),debt=(0.6,0.14,7),**DEAL))
    out[f'BE anchor{mw}']=be(DEAL); print('  BE NPV15 anchor',mw,out[f'BE anchor{mw}'])
M.DIRECT_MW[5]=11.0
r=M.site(5,'base',P); print('gas use: G_avail',r['G_avail'],'g_pow',r['g_pow'],'share',round(r['g_pow']/r['G_avail'],2),'burn share',round(r['g_pow']*0.72*0.92/r['G_avail'],2),'capex parts',r['capex_parts'],'MW',r['MW_gross'],'GWh',r['y1'].get('GWh'))
json.dump(out,open('deal_out.json','w'),ensure_ascii=False,indent=1,default=str)
