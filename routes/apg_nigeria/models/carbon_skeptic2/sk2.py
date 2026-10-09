import json
S='/tmp/claude-0/-home-user-ocm-research-diary/50ead64b-47e2-5bdd-b34b-09f6667b6ca9/scratchpad'
cng=json.load(open(S+'/cng/cng_out.json'))
def npv(c,r): return sum(x/(1+r)**t for t,x in enumerate(c))
NP={'VCM':2.32,'CORSIA':5.0,'ITMO':13.6}
V100=17.7  # kt/y per MMscf/d fully sold, pure CH4, PE 8% (lower bound)
for Q in (1,5,15):
  print('CNG npv/ebitda host', {k:(cng[f'{Q}_{k}']['ebitda1'],cng[f'{Q}_{k}']['npv15']) for k in ('low','base','high')})
  for cs in ('low','base'):
    frac=cng[f'{Q}_{cs}']['sold_scm_d']/(28321*Q)
    vol=V100*Q*frac
    for p in NP:
      row=[]
      for dem,mm,cm in ((1.3,1.7,2.4),(2.7,2.4,2.8)):
        rev=vol*NP[p]/dem/1000
        mrv=0.08*mm*(1.3 if Q>=5 else 1)
        cap=0.30*cm*(1.5 if Q>=5 else 1)
        e=rev-mrv
        n8=npv([-cap,0]+[e]*8,.15); n3=npv([-cap,0]+[e]*3,.15)
        # EV: issuance prob 0.6 on revenue, costs certain until refusal (MRV only if issued -> keep cost)
        ev8=npv([-cap,0]+[0.6*rev-mrv]*8,.15)
        row.append((round(rev,3),round(e,3),round(n3,2),round(n8,2),round(ev8,2)))
      print(Q,cs,round(frac,2),'kt',round(vol,1),p,row)
# LPG-only route: audit yield 2.0 t/d per MMscf/d, LF .8; screen 4.0, LF .95
for y,lf in ((2.0,0.8),(4.0,0.95)):
  kt=y*365*lf*3.0*0.92/1000; print('LPG kt per MMscf/d',y,lf,round(kt,2))
  for Q in (1,5,15):
    print('  Q',Q,'ITMO best rev',round(kt*Q*13.6/1.3/1000,3),'mrv',round(0.08*1.7*(1.3 if Q>=5 else 1),3))
