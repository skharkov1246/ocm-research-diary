import importlib.util,sys,io,contextlib
spec=importlib.util.spec_from_file_location('m','/tmp/claude-0/-home-user-ocm-research-diary/50ead64b-47e2-5bdd-b34b-09f6667b6ca9/scratchpad/minilng/minilng_model.py')
m=importlib.util.module_from_spec(spec)
with contextlib.redirect_stdout(io.StringIO()):
    # avoid overwriting their json: patch open? model writes out json; fine (same content)
    spec.loader.exec_module(m)
fav=dict(LF=0.9,decl=0.0,life=15,gas=0.25)
for Q in [5,15]:
    for nb in [8,10,12]:
        print('FAV',Q,nb,m.run(Q,'audit','audit_mid',nb=nb,**fav))
# power island capex: avg MW = 1.1 kWh/kg * t/d
for Q in [5,15]:
    t=m.lng_mmbtu_per_day(Q,1.1,.05,.03)/m.HHV_LNG
    mw=1.1*t*1000/24/1000
    lo=mw*1.1*1360/1e3; hi=mw*1.25*1960/1e3
    print('Q',Q,'avgMW',round(mw,2),'power capex',round(lo,1),round(hi,1))
    Cm=m.capex(Q,'audit_mid'); Ch=m.capex(Q,'audit_high'); Cl=m.capex(Q,'audit_low')
    Ch2=Ch+hi; Cm2=0.5*(Cl+Ch2)
    print(' corrected capex low/high/mid',round(Cl,1),round(Ch2,1),round(Cm2,1))
    for nb in [6,10,12]:
        print('  audit nb',nb,m.run(Q,'audit',nb=nb,C=Cm2))
        print('  fav   nb',nb,m.run(Q,'audit',nb=nb,C=Cm2,**fav))
    # breakeven with corrected capex audit and fav
    for lab,kw in [('audit',{}),('fav',fav)]:
        lo_,hi_=0,50
        for _ in range(60):
            mid=(lo_+hi_)/2
            if m.run(Q,'audit',nb=mid,C=Cm2,**kw)['NPV15']<0: lo_=mid
            else: hi_=mid
        print('  breakeven NPV15',lab,round(hi_,2))
# fav with only one lever relaxed at a time (Q=15, nb10, own mid capex)
for lab,kw in [('LF0.8',dict(LF=0.8)),('decl8',dict(decl=0.08)),('life12',dict(life=12)),('gas1',dict(gas=1.0)),('highcapex',dict(C=m.capex(15,'audit_high')))]:
    f=dict(fav); f.update(kw)
    print('FAV15 nb10 minus',lab,m.run(15,'audit','audit_mid',nb=10,**f))
