import importlib.util,sys,io,contextlib
spec=importlib.util.spec_from_file_location('m','/tmp/claude-0/-home-user-ocm-research-diary/50ead64b-47e2-5bdd-b34b-09f6667b6ca9/scratchpad/urea_ng/model.py')
m=importlib.util.module_from_spec(spec)
with contextlib.redirect_stdout(io.StringIO()): spec.loader.exec_module(m)
# 15 opt stresses
for pr in (480,430,400,380,347):
    r=m.run(15,'opt',price=pr); print('15opt price',pr,r['ebitda1'],r['npv15'],r['irr'])
for d in (0.03,0.05,0.08):
    m.S['opt']['decl']=d
    for pr in (480,430):
        r=m.run(15,'opt',price=pr); print('15opt decl',d,'price',pr,r['ebitda1'],r['npv15'],r['npv20'],r['irr'])
m.S['opt']['decl']=0.0
for gj in (27,):
    m.S['opt']['gj_t']=gj; r=m.run(15,'opt',price=430); print('15opt gj27 p430',r['kt_y'],r['ebitda1'],r['npv15'],r['irr'])
m.S['opt']['gj_t']=24
# 1 opt with netback 430
r=m.run(1,'opt',price=430); print('1opt p430',r['ebitda1'],r['npv15'],r['payback'])
# capex threshold check: 5 MMscf/d, capex 120 at price 450 mid other?
m.CAPEX[5]['mid']=120; r=m.run(5,'mid',price=450); print('5mid capex120 p450',r['ebitda1'],r['npv15'],r['irr'])
m.CAPEX[5]['opt']=120; r=m.run(5,'opt',price=450); print('5opt capex120 p450',r['ebitda1'],r['npv15'],r['irr'])
m.CAPEX[1]['mid']=25; r=m.run(1,'mid',price=450); print('1mid capex25 p450',r['ebitda1'],r['npv15'],r['irr'])
m.CAPEX[1]['opt']=25; r=m.run(1,'opt',price=450); print('1opt capex25 p450',r['ebitda1'],r['npv15'],r['irr'])
# AN check
print('AN skids', 22.3*0.425/(12.2*0.567), 'region', 22+75+33+14.5+2.4, 22.4+118+57+14.5+3.3)
