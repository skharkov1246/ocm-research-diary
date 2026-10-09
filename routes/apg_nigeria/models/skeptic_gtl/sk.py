import json
src=open('../gtl/gtl_model.py').read().split("out={}")[0]
exec(src)
# 1) IRR "не существует": проверим знак/корни на сетке
def cfs(sc,q,**kw):
    p=dict(SC[sc]); p.update(kw.get('over',{}))
    SC['_t']=p; return p
def cashflows(sc,q,diesel=None,capex=None,pre_decl=False,trk=0.0,gasp=None,top=False,lpg=0.0,extra_capex=0.0,yld=None,fx=0.0):
    p=dict(SC[sc])
    if diesel: p['diesel']=diesel
    if gasp is not None: p['gas']=gasp
    if yld: p['Y']=yld
    C=(capex if capex is not None else p['capex'][q])+extra_capex
    k=STAFF_SCALE[q]
    price=(p['dshare']*p['diesel']+(1-p['dshare'])*p['naph_k']*p['diesel'])*(1-fx)
    cf=[-C/p['build']]*p['build']; e2=None
    g0=(1-p['decl'])**p['build'] if pre_decl else 1.0
    for y in range(1,LIFE+1):
        g=g0*(1-p['decl'])**(y-1)*(p['ramp'] if y==1 else 1.0)
        bbl=p['Y']*q*p['avail']*365*g
        rev=bbl*(price-trk)/1e6 + lpg*q*g*(p['avail'])
        gasvol=MMBTU_PER_MSCF*1000*q*365*(1.0 if top else p['avail']*g)
        opex=(bbl*(p['opex_bbl']+p['cat'])/1e6+gasvol*p['gas']/1e6+p['staff']*k+p['sec']*k**0.7+(p['tor']+p['ins'])*C+p['levy']*rev)
        e=rev-opex; dep=C/10 if y<=10 else 0
        cf.append(e-max(0,TAX*(e-dep)))
        if y==2: e2=e
    return cf,e2,C
def rep(tag,cf,e2,C):
    print(f"{tag:60s} C={C:6.1f} EBITDA2={e2:6.2f} NPV15={npv(.15,cf):7.1f} NPV20={npv(.20,cf):7.1f}")
# IRR base q15
cf,_,_=cashflows('base',15)
import math
for r in (-0.3,-0.2,-0.1,-0.05,0.0,0.05): print('base q15 NPV@',r,round(npv(r,cf),1))
# 2) поправки против GTL (не учтены в базе)
for q in (1,15):
    rep(f'base q{q} as-is',*cashflows('base',q))
    rep(f'base q{q} + падение газа за стройку',*cashflows('base',q,pre_decl=True))
    rep(f'base q{q} + стройка-падение + take-or-pay + доставка $5/bbl + FX 4.7%',*cashflows('base',q,pre_decl=True,top=True,trk=5,fx=0.047))
    rep(f'base q{q} + всё выше + Fourth Schedule $2.18',*cashflows('base',q,pre_decl=True,top=True,trk=5,fx=0.047,gasp=2.18))
# 3) за GTL: LPG из C3+ (−17% энергии на ФТ, +~1000 т/г×$650/т на 1 MMscfd, +$3 млн/MMscfd^0.6 капитал)
for q in (1,15):
    rep(f'base q{q} + LPG ко-продукт',*cashflows('base',q,yld=80*0.83,lpg=0.65,extra_capex=3*q**0.6))
# 4) лучший «справедливый» случай на 15 MMscf/d
for d in (150,175,200):
  for C in (116,162,242,270):
    rep(f'q15 base-yield diesel {d} capex {C}',*cashflows('base',15,diesel=d,capex=C))
for d in (150,175):
    rep(f'q15 opt diesel {d} capex 200',*cashflows('opt',15,diesel=d))
    rep(f'q15 opt diesel {d} capex 200 + стройка-падение(0)',*cashflows('opt',15,diesel=d))
SC['opt5']=dict(SC['opt']); SC['opt5']['decl']=0.05
SC['opt10']=dict(SC['opt']); SC['opt10']['decl']=0.10
rep('q15 opt decl 5%',*cashflows('opt5',15)); rep('q15 opt decl 10%',*cashflows('opt10',15))
rep('q15 opt decl 5% + pre-decl',*cashflows('opt5',15,pre_decl=True))
# 5) энергетика: HHV ФТ-жидкостей
dens=0.765; m=dens*158.987; print('FT HHV MMBtu/bbl',round(m*47.3/1055.06,2),' ceiling eff @100bpd', round(100*m*47.3/1055.06/1050,3), 'base80',round(80*m*47.3/1055.06/1050,3))
