# Скептик: доповерка cng_model.run — капитал на постройку, рабочий капитал, разгон клиентов, амортизация 7 лет (а не 10/7), перенос убытков, остаточная стоимость парка
import sys,io,contextlib,math,json
src=open(sys.argv[1]).read(); g={'__file__':sys.argv[1]}
with contextlib.redirect_stdout(io.StringIO()): exec(src,g)
run,S,npv,irr,conv=g['run'],g['S'],g['npv'],g['irr'],g['ngn_scm_to_usd_mmbtu']
def adj(Q,p,build=1,ramp=(0.5,0.8,1.0),dso=60,resid=0.25,carry=True,dep_fix=True,years=10,tax=0.34,wacc=(0.15,0.20)):
    base=run(Q,p)   # для капекса/опекса 1-го года
    capex=base['capex']; fleet=base['fleet']; station=capex-fleet
    sold1=base['sold_scm_d']; rev_full=base['rev1']
    fixed=base['fixed1']; var_full=base['gas1']+base['truck1']
    cf=[]
    # постройка: станция+девелопмент в t0, парк в конце постройки
    if build>=1:
        cf=[-station]+[0]*(build-1)+[-fleet]
    else: cf=[-capex]
    wc_prev=0; loss=0; e_list=[]
    for t in range(1,years+1):
        r=ramp[t-1] if t-1<len(ramp) else 1.0
        gdec=(1-p['decl'])**(t-1)
        rev=rev_full*r*gdec
        e=rev-var_full*r*gdec-fixed
        e_list.append(e)
        dep=capex/7 if (not dep_fix or t<=7) else 0
        ti=e-dep
        if carry:
            ti2=ti-loss
            if ti2<0: loss=-ti2; taxv=0
            else: loss=0; taxv=ti2*tax
        else: taxv=max(0,ti*tax)
        wc=rev*dso/365; dwc=wc-wc_prev; wc_prev=wc
        repl=base['tractors']*p['tractor'] if t==6 else 0
        c=e-taxv-dwc-repl
        if t==years: c+=wc_prev+resid*fleet
        if len(cf)>build+t: cf[build+t]+=c
        else: cf.append(c)
    return dict(capex=capex,ebitda_y1=round(e_list[0],2),ebitda_peak=round(max(e_list),2),npv15=round(npv(0.15,cf),2),npv20=round(npv(0.20,cf),2),irr=(round(irr(cf)*100,1) if irr(cf) is not None else None))
def be(Q,p,key='npv15',**kw):
    lo,hi=100,4000
    for _ in range(60):
        m=(lo+hi)/2
        if adj(Q,dict(p,price=m),**kw)[key]>0: hi=m
        else: lo=m
    return round(hi),round(conv(hi),2)
pc=dict(S['base']);pc.update(km=60,util=0.75,staff=0.55*1.8,sec=0.3*1.6,overhead=0.2*1.5)
cases={'5_contract':(5,pc),'1_high':(1,S['high']),'1_base':(1,S['base'])}
p5h=dict(S['high']);p5h['km']+=60;p5h['util']*=0.85;p5h['staff']*=1.8;p5h['sec']*=1.6;p5h['overhead']*=1.5
cases['5_high']=(5,p5h)
out={}
for nm,(Q,p) in cases.items():
    # 1) только техфиксы (dep 7 лет, перенос убытков, остаточная 25% парка), без постройки/разгона/WC
    out[nm+'_techfix']=adj(Q,p,build=0,ramp=(),dso=0)
    out[nm+'_full']=adj(Q,p)
    out[nm+'_be15_full']=be(Q,p); out[nm+'_be20_full']=be(Q,p,'npv20')
    out[nm+'_be15_techfix']=be(Q,p,build=0,ramp=(),dso=0)
for pr in [620,750,900]:
    out[f'5_contract_p{pr}_full']=adj(5,dict(pc,price=pr))
    out[f'5_contract_p{pr}_full_build2']=adj(5,dict(pc,price=pr),build=2)
# каждая поправка по отдельности для 5_contract ₦620
for nm,kw in [('only_build1',dict(ramp=(),dso=0)),('only_ramp',dict(build=0,dso=0)),('only_wc60',dict(build=0,ramp=())),('only_wc90',dict(build=0,ramp=(),dso=90)),('none_but_resid',dict(build=0,ramp=(),dso=0))]:
    out['5c620_'+nm]=adj(5,dict(pc,price=620),**kw)
print(json.dumps(out,indent=0,ensure_ascii=False))
p620=dict(pc,price=620)
for nm,kw in [('repro',dict(build=0,ramp=(),dso=0,resid=0,dep_fix=False,carry=False)),('+depfix',dict(build=0,ramp=(),dso=0,resid=0,carry=False)),('+carry',dict(build=0,ramp=(),dso=0,resid=0)),('+resid25',dict(build=0,ramp=(),dso=0))]:
    print(nm,adj(5,p620,**kw))
