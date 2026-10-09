import io,contextlib
exec(open('pm.py').read())
def show(nm,sc,q=1,o=None,power=True):
    r=run(sc,q,power,over=o)
    y=r['y2']
    print(f"{nm:55s} capex {r['capex']:6.2f} y2EBITDA {y['ebitda']:6.2f} rev_c {y['rev_c']:.2f} rev_e {y['rev_e']:.2f} opex {y['opex']:.2f} NPV15 {r['npv15']:7.2f} NPV20 {r['npv20']:7.2f} IRR {r['irr']} pb {r['simple_payback']}")
    return r
print('--- optimism thermo-consistency (1 MMscf/d)')
show('high as published (X0.95, 20 t melt => ~15 bar)','high')
show('high X0.86 @20t (15 bar, Xeq 0.85-0.88)','high',o=dict(X=0.86))
show('high X0.95 @120t (~2.5 bar)','high',o=dict(melt_t=120))
show('high X0.95 @120t, pyro capex x(120/20)^0.6 on vessel? -> 9.6*1.3','high',o=dict(melt_t=120,pyro_capex=12.5))
for q in (5,15):
    show(f'high q={q} published','high',q)
    show(f'high q={q} X0.86','high',q,dict(X=0.86))
    show(f'high q={q} X0.95 melt120','high',q,dict(melt_t=120))
print('--- year-2 money stresses, optimism')
show('high FX 2000 unindexed (400->266)','high',o=dict(price_ngn=400*1330/2000))
show('high avail 0.6 FOAK','high',o=dict(avail=0.6))
show('high netback 200 (freight 100-150 vs 40-80)','high',o=dict(c_net=200))
show('high X0.86 + FX2000 + netback200','high',o=dict(X=0.86,price_ngn=266,c_net=200))
show('high X0.86 + FX2000 + netback200 + avail0.7','high',o=dict(X=0.86,price_ngn=266,c_net=200,avail=0.7))
print('--- base stresses')
show('base','base')
show('base netback 100','base',o=dict(c_net=100))
show('base avail 0.6','base',o=dict(avail=0.6))
show('base FX2000','base',o=dict(price_ngn=320*1330/2000))
show('base all three','base',o=dict(c_net=100,avail=0.6,price_ngn=213))
show('base CH4 7000 t (phys) not 7300 -> scale by tor', 'base')
