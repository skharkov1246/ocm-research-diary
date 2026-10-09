import sys
exec(open('gtl_model.py').read().split("out={}")[0])
scr=SC['screen']['capex']
print('--- opt с капиталом по правилу аудита x2.4/x2.8 к скрину')
for q in (1,5,15):
    for m in (2.4,2.8):
        print(q,m,run('opt',q,capex=scr[q]*m))
print('--- opt, капитал аудит + выход по нижней границе Y=80')
SC['optY80']=dict(SC['opt']); SC['optY80']['Y']=80
for q in (15,):
    for m in (2.4,2.8):
        print(q,m,run('optY80',q,capex=scr[q]*m), 'BEdiesel', be_price('optY80',q,capex=scr[q]*m))
print('--- base с капиталом строго x2.4/x2.8')
for q in (1,5,15):
    for m in (2.4,2.8):
        print(q,m,run('base',q,capex=scr[q]*m))
print('--- base при текущей цене ворот Dangote $171/$222')
for q in (15,):
    for d in (171,222):
        print(q,d,run('base',q,diesel=d))
print('--- opt wax 1500 c ценой воска на уровне импортного CIF 1000-1300 и аудит капиталом x2.4')
for w in (1000,1300):
    print(w,run('opt',15,wax=(w,2.0)), run('opt',15,wax=(w,2.0),capex=scr[15]*2.4))
# NDDC 3% бюджета (капекс + опекс): капекс x1.03 и опекс +3%
print('--- opt q15 + NDDC 3% (капекс*1.03, опекс через levy +3% выручки как прокси)')
SC['optN']=dict(SC['opt']); SC['optN']['levy']=SC['opt']['levy']+0.03*15.45/74.95
print(run('optN',15,capex=200*1.03))
# Брент-паритет: нафта смешивается с нефтью -> нафта = 0.6*дизель; и потеря 10% цены на офф-спек плотность
SC['optD']=dict(SC['opt']); 
print('opt q15 diesel 175*0.9 offspec', run('opt',15,diesel=175*0.9))
# 10-летний срок разрешения/ресурса
LIFE=10
print('opt q15 LIFE=10', run('opt',15))
