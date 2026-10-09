import sys; sys.argv=['x']
exec(open('gtl_model.py').read().split("out={}")[0])
def be_capex(sc,q,diesel,r='npv15'):
    lo,hi=1,1000
    for _ in range(80):
        m=(lo+hi)/2
        if run(sc,q,diesel=diesel,capex=m)[r]>0: lo=m
        else: hi=m
    return round(m,1)
for q in (1,5,15):
    for d in (90,120,175,200):
        print('BEcapex base q',q,'diesel',d,'->',be_capex('base',q,d),'npv20:',be_capex('base',q,d,'npv20'))
# против себя: ТОиР 2.5%+страх 0.8% на 15 MMscfd, Y=90
SC['base2']=dict(SC['base']); SC['base2'].update(tor=0.025,ins=0.008,Y=90,decl=0.03)
for q in (5,15): print('base_softer q',q,run('base2',q), 'BE diesel', be_price('base2',q))
# аудитный множитель в чистом виде: скрин с капексом x2.4/x2.8 и опексом x1.7/x2.4
for q in (1,5,15):
  for cm,om in ((2.4,1.7),(2.8,2.4)):
    SC['aud']=dict(SC['screen']); SC['aud']['opex_bbl']=10*om; SC['aud']['gas']=0.6
    SC['aud']['capex']={k:v*cm for k,v in SC['screen']['capex'].items()}
    print('screen*audit q',q,cm,om,run('aud',q))
