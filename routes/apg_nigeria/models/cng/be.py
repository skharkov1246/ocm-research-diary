import sys,io,contextlib,json
src=open(sys.argv[1]).read()
g={'__file__':sys.argv[1]}
with contextlib.redirect_stdout(io.StringIO()): exec(src,g)
run,S,conv=g['run'],g['S'],g['ngn_scm_to_usd_mmbtu']
def be(Q,p,r='npv15'):
    lo,hi=100,3000
    for _ in range(60):
        m=(lo+hi)/2; x=run(Q,dict(p,price=m))
        if x[r]>0: hi=m
        else: lo=m
    return round(hi),round(conv(hi),2)
res={}
for k in ['screen','low','base','high']:
    res[f'1_{k}_NPV15=0']=be(1,S[k]); res[f'1_{k}_NPV20=0']=be(1,S[k],'npv20')
for Q,km_add,um,sm,se,ov in [(5,60,0.85,1.8,1.6,1.5),(15,150,0.6,2.8,2.5,2.2)]:
    for k in ['low','base','high']:
        p=dict(S[k]);p['km']+=km_add;p['util']*=um;p['staff']*=sm;p['sec']*=se;p['overhead']*=ov
        res[f'{Q}_{k}_NPV15=0']=be(Q,p); res[f'{Q}_{k}_NPV20=0']=be(Q,p,'npv20')
# contracted near-hub 5 MMscf/d: base costs, km 60, util 0.75, price 620
p=dict(S['base']);p.update(km=60,util=0.75,staff=0.55*1.8,sec=0.3*1.6,overhead=0.2*1.5)
for pr in [480,620,750]:
    x=run(5,dict(p,price=pr)); res[f'5_contract_km60_u75_p{pr}']={k:x[k] for k in ['capex','rev1','ebitda1','cash_cost_usd_mmbtu','simple_payback','npv15','npv20','irr','loads_d','trailers','tractors','customers']}
res['5_contract_be15']=be(5,p); res['5_contract_be20']=be(5,p,'npv20')
# 3 MMscf/d contracted
p3=dict(S['base']);p3.update(km=60,util=0.75,staff=0.55*1.4,sec=0.3*1.3,overhead=0.2*1.25)
res['3_contract_be15']=be(3,p3); x=run(3,dict(p3,price=620)); res['3_contract_p620']={k:x[k] for k in ['capex','rev1','ebitda1','npv15','npv20','irr','simple_payback']}
print(json.dumps(res,indent=0,ensure_ascii=False))
