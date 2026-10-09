import sys
sys.argv=['x','../cng/cng_model.py']
exec(open('adj.py').read().split("print(json.dumps")[0])
# Z метана при 250 бар, 15C ~0.83 [Д]; стальной 40-фт каскад 18-24 м3 водяного объёма [Д]
for wv in (18,22,26): print('steel water',wv,'m3 -> scm',round(wv*250/1.013/0.83))
for nm,ov in [('model',{}),('steel_4500_usable',dict(usable_scm=4500)),('steel_5500_usable',dict(usable_scm=5500)),('composite_7000_0.60M',dict(trailer=0.60)),('composite_7000_0.69M',dict(trailer=0.69))]:
    p=dict(pc,price=620); p.update(ov); r=run(5,p)
    print(nm,{k:r[k] for k in ['capex','fleet','trailers','tractors','loads_d','ebitda1','npv15','irr']}, 'adj_full',adj(5,p), 'be15_full',be(5,dict(pc,**ov)))
