import sys,io,contextlib,importlib
sys.path.insert(0,'/tmp/claude-0/-home-user-ocm-research-diary/50ead64b-47e2-5bdd-b34b-09f6667b6ca9/scratchpad/skeptic_pyro')
with contextlib.redirect_stdout(io.StringIO()): import m
def show(tag,sc='base',power=True,**o):
    r=m.run(sc,1,power=power,over=o); y=r['y2']
    print(f"{tag:45s} capex {r['capex']:6.2f} MWnet {r['balance']['mw_net']:.2f} revC {y['rev_c']:.2f} revE {y['rev_e']:.2f} opex {y['opex']:.2f} EBITDA {y['ebitda']:6.2f} NPV15 {r['npv15']:7.2f} IRR {r['irr']}")
show('as published')
m.MW_DIRECT=4.40; show('(1) LHV-basis direct genset 4.40 MW')
show('(1)+(2) heat 210', heat=210)
show('(1)+(2)+(3) parasit 2MWh/tC ~ 0.40', heat=210, parasit=0.40)
show('(1)+(2) +(4) bi_loss 2 kg/t', heat=210, bi_loss=2.0)
show('(1)+(2)+(4)+(5) netback 150', heat=210, bi_loss=2.0, c_net=150)
show('(1)+(2)+(4)+(5)+(6) melt 60t', heat=210, bi_loss=2.0, c_net=150, melt_t=60)
show('corrected base', heat=210, bi_loss=2.0, c_net=150, melt_t=60)
show('corrected + pyro x5.3 FOB-rule (21.2)', heat=210, bi_loss=2.0, c_net=150, melt_t=60, pyro_capex=21.2)
show('optimism as published','high')
show('optimism heat 175, MW 4.40','high',heat=175)
show('optimism heat 175 netback 300','high',heat=175,c_net=300)
show('corrected carbon-only',power=False, heat=210, bi_loss=2.0, c_net=150, melt_t=60)
# corrected pessimism
show('pess heat 270','low',heat=270)
