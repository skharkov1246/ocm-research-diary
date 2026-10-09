import sys, importlib.util, io, contextlib
spec=importlib.util.spec_from_file_location('m','/tmp/claude-0/-home-user-ocm-research-diary/50ead64b-47e2-5bdd-b34b-09f6667b6ca9/scratchpad/lpg/lpg_model.py')
m=importlib.util.module_from_spec(spec)
with contextlib.redirect_stdout(io.StringIO()): spec.loader.exec_module(m)
# corrected yields from PR flash: rich 3.48 flash-only at -35C/27.5bar, NJTD 2.73 after columns; mid ~0.7
kmol=m.KMOL_PER_MMSCF
def y(z,rec):
    return sum(z[c]/100*kmol*m.MW[c]/1000*rec[c] for c in rec)
rich=dict(C3=5.68,iC4=1.56,nC4=1.60); mid=dict(C3=2.3,iC4=0.6,nC4=0.6)
print('rich -35/27.5 flash', round(y(rich,dict(C3=.56,iC4=.79,nC4=.85)),2),' -40/40', round(y(rich,dict(C3=.71,iC4=.86,nC4=.91)),2))
print('mid  -35/27.5 flash', round(y(mid,dict(C3=.23,iC4=.46,nC4=.58)),2),' -40/40', round(y(mid,dict(C3=.38,iC4=.62,nC4=.72)),2))
for yy in (0.75,2.7,3.0,3.5,4.0):
  for p in (600,650):
    for cm in (1.6,2.4,2.8):
      r=m.run2(15,'addon',y=yy,price=p,cmult=cm)
      print(yy,p,cm,'capex',r['capex'],'rev',r['rev1'],'EBITDA',r['ebitda1'],'NPV15',r['NPV15'],'NPV20',r['NPV20'],'IRR',r['IRR'],'pb',r['payback'])
# shrinkage opportunity cost if anchor gas-limited
for yy in (3.0,4.0):
    e=(yy*1.2)*15*365*0.8*47.3
    print('y',yy,'MMBtu/y removed',round(e), 'lost anchor EBITDA at $6.2 (mining 27.7)',round(e*6.2/1e6,2),'at $7.6-15.6 (captive)',round(e*7.6/1e6,2),round(e*15.6/1e6,2))
