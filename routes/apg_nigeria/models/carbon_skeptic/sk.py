def npv(c,r): return sum(x/(1+r)**t for t,x in enumerate(c))
NP={'VCM':2.32,'CORSIA':5.0,'ITMO':13.6}
V0=17.7  # kt/y per 1 MMscf/d fully recovered, pure CH4, PE 8% (lower bound per audit rule)
# CNG sold fraction util*sell from cng_model (low/base/high), scale mult at 5/15
F={'low':0.40*0.85,'base':0.65*0.88,'high':0.80*0.90}
SC={1:1.0,5:0.85,15:0.6}
for Q in (1,5,15):
  for cs in ('low','base','high'):
    vol=V0*Q*F[cs]*SC[Q]
    for p in NP:
      for dem in (1.3,2.7):
        rev=vol*NP[p]/dem/1000
        for mrvm in (1.7,2.4):
          mrv=0.08*mrvm*(1.3 if Q>=5 else 1)
          cap=0.30*(2.4 if mrvm==1.7 else 2.8)*(1.5 if Q>=5 else 1)
          e=rev-mrv
          for yrs in (3,8):
            c=[-cap,0]+[e]*yrs
            if (dem==1.3 and mrvm==1.7) or (dem==2.7 and mrvm==2.4):
              tag='best' if dem==1.3 else 'worst'
              print(Q,cs,p,tag,yrs,'vol_kt',round(vol,1),'rev',round(rev,3),'ebitda',round(e,3),'npv15',round(npv(c,.15),2),'EV',round(0.6*npv(c,.15),2))
# LPG-only route: LPG 2-4 t/d per MMscf/d, ~3.0 t CO2/t LPG, PE 8%
for y in (2,4):
  print('LPG kt/y per MMscf/d', round(y*365*3.0*0.92/1000,2))
# lean ND carbon per mol
print('lean_ND C/mol',0.8792+2*.0465+3*.0093+4*.0055+5*.0018)
