import importlib.util, io, contextlib
spec=importlib.util.spec_from_file_location('m','/tmp/claude-0/-home-user-ocm-research-diary/50ead64b-47e2-5bdd-b34b-09f6667b6ca9/scratchpad/lpg/lpg_model.py')
m=importlib.util.module_from_spec(spec)
with contextlib.redirect_stdout(io.StringIO()): spec.loader.exec_module(m)
# Otakikpo-derived capex: >$60M (LPG 12MMscfd + 20MW). power 20MW at 0.9-1.5 $M/MW
for pw in (0.9,1.5):
    lpg12=60-20*pw; c15=lpg12*(15/12)**0.6
    print('power $/MW',pw,'LPG12',lpg12,'scaled15',round(c15,1),'mult vs screen',round(c15/18.2,2))
for mode in ('addon','standalone'):
  for y in (0.71,2.7,3.0,3.5):
    for p in (500,600,650):
      for cm in (1.9,2.4):
        r=m.run2(15,mode,y=y,price=p,cmult=cm)
        print(mode,y,p,cm,'capex',r['capex'],'rev',r['rev1'],'EBITDA',r['ebitda1'],'NPV15',r['NPV15'],'NPV20',r['NPV20'],'IRR',r['IRR'],'pb',r['payback'])
# yield needed for NPV15=0, addon, 15, cm 1.9 and 2.4, price 600/650, LF .8 decl .08 and LF .85 decl .05
for cm in (1.9,2.4):
  for p in (550,600,650):
    for LF,d in ((0.8,0.08),(0.85,0.05)):
      lo,hi=0.5,12
      for _ in range(60):
        mid=(lo+hi)/2
        if m.run2(15,'addon',y=mid,price=p,cmult=cm,LF=LF,decl=d)['NPV15']<0: lo=mid
        else: hi=mid
      print('breakeven y NPV15=0 cm',cm,'p',p,'LF',LF,'decl',d,round(hi,2))
# market: July 2026 NMDPRA
print('import gap july t/y', round((4.4-4.373)*1000*365), 'june gap', round((4.1-3.6)*1000*365))
print('Lagos retail $/t', round(1235e3/1330), round(1776e3/1330))
