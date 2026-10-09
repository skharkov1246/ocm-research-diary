exec(open('pm.py').read())
def show(nm,sc,q=1,o=None,power=True):
    r=run(sc,q,power,over=o); y=r['y2']
    print(f"{nm:48s} capex {r['capex']:6.2f} y2EBITDA {y['ebitda']:6.2f} NPV15 {r['npv15']:7.2f} NPV20 {r['npv20']:7.2f} IRR {r['irr']} pb {r['simple_payback']} dcfpb {r['dcf_payback_yr']}")
    return r
fix=dict(melt_t=120,c_net=200)
for q in (1,5,15):
    show(f'opt corrected (120t melt, netback 200) q={q}','high',q,fix)
    show(f'opt corr + FX2000 q={q}','high',q,dict(fix,price_ngn=266))
show('opt corr carbon-only','high',1,fix,power=False)
# against myself: genset bypass on raw gas when pyro down -> power availability 0.95 instead of 0.85
# approximate by raising avail only for power: emulate via LF scaling 0.95/0.85
show('base + bypass (power avail 0.95 approx via LF)','base',o=dict(LF=0.72*0.95/0.85))
show('opt corr + bypass','high',o=dict(fix,LF=min(1,0.85*0.95/0.92)))
# DSCR
for nm,cap,e in [('base',23.97,0.29),('opt pub',18.03,5.13),('opt corr+FX',None,None)]:
    if cap: 
        ds=0.6*cap*0.15/(1-1.15**-7); print(nm,'debt service',round(ds,2),'DSCR',round(e/ds,2))
