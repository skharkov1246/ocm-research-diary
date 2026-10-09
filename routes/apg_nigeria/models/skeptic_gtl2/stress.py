# Стресс "второй год и деньги" поверх gtl_model.py (без изменения его формул)
import json
src=open('/tmp/claude-0/-home-user-ocm-research-diary/50ead64b-47e2-5bdd-b34b-09f6667b6ca9/scratchpad/gtl/gtl_model.py').read().split("out={}")[0]
exec(src)
def run2(sc,q,C,diesel=None,Y=None,y2=None,outage=None,esc=0.0,catcap=0.0,life=15,wc=0.0):
    p=dict(SC[sc]); 
    if diesel: p['diesel']=diesel
    if Y: p['Y']=Y
    k=STAFF_SCALE[q]; price=p['dshare']*p['diesel']+(1-p['dshare'])*p['naph_k']*p['diesel']
    cf=[-C/p['build']]*p['build']; e2=None
    for y in range(1,life+1):
        g=(1-p['decl'])**(y-1)*(p['ramp'] if y==1 else 1.0)
        if y2 and y==2: g*=y2
        if outage and y==outage[0]: g*=outage[1]
        bbl=p['Y']*q*p['avail']*365*g; rev=bbl*price/1e6
        gas=MMBTU_PER_MSCF*1000*q*365*p['avail']*g*p['gas']/1e6
        fixed=(p['staff']*k+p['sec']*k**0.7+(p['tor']+p['ins'])*C)*(1+esc)**(y-1)
        opex=bbl*(p['opex_bbl']+p['cat'])/1e6+gas+fixed+p['levy']*rev
        e=rev-opex
        if catcap and y in (5,10): e-=catcap
        if wc and y==1: e-=wc*rev
        dep=C/10 if y<=10 else 0
        tx=max(0,TAX*(e-dep)); cf.append(e-tx)
        if y==2: e2=round(e,2)
    return dict(C=C,ebitda_y2=e2,npv15=round(npv(.15,cf),1),npv20=round(npv(.20,cf),1),irr=(round(irr(cf)*100,1) if irr(cf) is not None else None))
R={}
scr=SC['screen']['capex']
for q in (15,):
  for m in (1.67,2.4,2.8):
    C=scr[q]*m
    R[f'opt_q{q}_x{m}']=run2('opt',q,C)
    R[f'opt_q{q}_x{m}_Y80']=run2('opt',q,C,Y=80)
    R[f'opt_q{q}_x{m}_y2avail0.6']=run2('opt',q,C,y2=0.6/0.92)
    R[f'opt_q{q}_x{m}_outage6m_y4']=run2('opt',q,C,outage=(4,0.5))
    R[f'opt_q{q}_x{m}_fixedesc5']=run2('opt',q,C,esc=0.05)
    R[f'opt_q{q}_x{m}_diesel140']=run2('opt',q,C,diesel=140)
    R[f'opt_q{q}_x{m}_diesel120']=run2('opt',q,C,diesel=120)
    R[f'opt_q{q}_x{m}_ALL']=run2('opt',q,C,Y=80,y2=0.6/0.92,outage=(4,0.5),esc=0.05,diesel=140)
# Exact "opt" as stated (capex 200) with second-year stress only
R['opt_q15_200_y2_0.6_out_esc']=run2('opt',15,200,y2=0.6/0.92,outage=(4,0.5),esc=0.05)
# сколько капитала выдерживает opt q15 (NPV15=0) при Y80 и дизеле 140/175
def becap(**kw):
    lo,hi=1,1000
    for _ in range(60):
        m=(lo+hi)/2
        if run2('opt',15,m,**kw)['npv15']>0: lo=m
        else: hi=m
    return round(m,1)
R['BEcap_opt15_Y80_d175']=becap(Y=80)
R['BEcap_opt15_Y95_d175']=becap()
R['BEcap_opt15_Y80_d140']=becap(Y=80,diesel=140)
R['BEcap_opt15_Y80_d175_npv20']=None
for k,v in R.items(): print(k,v)
json.dump(R,open('stress_out.json','w'),ensure_ascii=False,indent=1)
