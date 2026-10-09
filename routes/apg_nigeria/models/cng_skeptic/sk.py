import sys,io,contextlib,json,math
src=open('../cng/cng_model.py').read()
g={'__file__':'/dev/null/cng_model.py'}
src=src.replace("json.dump(out,open(__file__.replace('cng_model.py','cng_out.json'),'w'),indent=1,ensure_ascii=False)","")
with contextlib.redirect_stdout(io.StringIO()): exec(src,g)
S=g['S'];npv=g['npv'];irr=g['irr'];SCM=g['SCM_PER_MMSCF'];MM=g['MMBTU_PER_SCM'];FX=1330
def run2(Q,p,mode='mult',delay=0,ramp=1.0,years=10,tax=0.34,residual=0.0):
    station=p['st_pkg']*Q**0.7*p['cmult']
    inlet=Q*SCM
    sold_d=inlet*p['sell']*p['util']; loads=sold_d/p['usable_scm']
    trip=2*p['km']/p['speed']+p['cyc_extra']
    tractors=math.ceil(loads*trip/p['drive_h']*1.15); ncust=math.ceil(sold_d/p['cust_scm'])
    trailers=math.ceil((loads*12/24+loads*trip/24+ncust)*1.15)
    capex=station+trailers*p['trailer']+tractors*p['tractor']+ncust*p['cust_kit']+p['dev']
    cf=[-capex]+[0.0]*delay
    e1=None
    for t in range(1,years+1):
        gg=(1-p['decl'])**(t-1)
        f = p['util']*gg if mode=='mult' else min(p['util'],gg)   # доля входа, реально проданная
        if t==1: f*=ramp
        sold=inlet*p['sell']*f*365
        rev=sold*p['price']/FX/1e6*(1-p['bad'])
        gas=inlet*f*365*MM*p['gas']/1e6
        km=(sold/365/p['usable_scm'])*2*p['km']*365
        fixed=p['staff']+p['sec']+p['overhead']+tractors*2*6000/1e6+capex*(p['tor']+p['ins'])
        e=rev-gas-km*p['opx_km']/1e6-fixed
        if e1 is None: e1=e
        c=e-max(0,(e-capex/7)*tax)-(tractors*p['tractor'] if t==6 else 0)
        if t==years: c+=residual*trailers*p['trailer']
        cf.append(c)
    return dict(capex=round(capex,2),ebitda1=round(e1,2),npv15=round(npv(.15,cf),2),npv20=round(npv(.20,cf),2),irr=(round(irr(cf)*100,1) if irr(cf) else None))
def be(Q,p,**kw):
    lo,hi=100,3000
    for _ in range(50):
        m=(lo+hi)/2
        if run2(Q,dict(p,price=m),**kw)['npv15']>0: hi=m
        else: lo=m
    return round(hi)
pc=dict(S['base']);pc.update(km=60,util=0.75,staff=0.55*1.8,sec=0.3*1.6,overhead=0.2*1.5)
cases={'1_base':(1,S['base']),'1_high':(1,S['high']),'5_contract':(5,pc)}
for nm,(Q,p) in cases.items():
    print(nm)
    for lab,kw in [('as_model',{}),('demand_limited_min',dict(mode='min')),('min+delay1+ramp0.6',dict(mode='min',delay=1,ramp=0.6)),('mult+delay1+ramp0.6',dict(delay=1,ramp=0.6)),('min+delay2+ramp0.6',dict(mode='min',delay=2,ramp=0.6)),('min+resid30%',dict(mode='min',residual=0.3))]:
        print(' ',lab,run2(Q,p,**kw),'BE15 N/scm=',be(Q,p,**kw))
# local staff only for CNG (no expat): staff 0.30 at 1 MMscf/d
for nm,(Q,p) in cases.items():
    pp=dict(p); pp['staff']=p['staff']*0.30/0.55
    print('noexpat',nm,run2(Q,pp,mode='min',delay=1,ramp=0.6),'BE',be(Q,pp,mode='min',delay=1,ramp=0.6))
# 5 contract at 620 & 750 under min+delay
for pr in [480,620,750]:
    print('5c',pr,run2(5,dict(pc,price=pr),mode='min',delay=1,ramp=0.6), run2(5,dict(pc,price=pr),mode='min'))
print('----')
for pr in [380,450,520]:
    print('5c',pr,run2(5,dict(pc,price=pr),mode='min',delay=1,ramp=0.6))
# permit-holder margin: gas $1.5; extra conditioning capex: dev +1.0*Q^0.7
for lab,mod in [('gas1.5',dict(gas=1.5)),('cond+',dict(dev=0.7+1.0*5**0.7)),('both',dict(gas=1.5,dev=0.7+1.0*5**0.7))]:
    p=dict(pc);p.update(mod)
    print(lab,run2(5,dict(p,price=620),mode='min',delay=1,ramp=0.6),'BE',be(5,p,mode='min',delay=1,ramp=0.6))
print(620/1330/0.03637, 655/1330/0.03637, 734/1330/0.03637, 790/1330/0.03637)
