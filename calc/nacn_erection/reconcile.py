# -*- coding: utf-8 -*-
"""Part 6 - reconcile the stick-built network with its own man-hour take-off."""
HRS,TT=60,0.48
MH_CONSTR_STICK=16681.0   # piping+E&I+mech/steel only, from schedule.py
need_mw = MH_CONSTR_STICK/(HRS*TT)
print(f"stick-built construction take-off {MH_CONSTR_STICK:.0f} mh / ({HRS} h x {TT}) "
      f"= {need_mw:.0f} man-weeks REQUIRED for erection+piping+E&I")
NET=[("site establishment / camp / temp utilities", 0, 6,  8, 0),
     ("civil incl. in-situ plinths + 28 d cure",     3,15, 14, 0),
     ("equipment erection + structural steel",      13,23, 15, 1),
     ("piping ~157 Cat-M + ~185 utility welds",     18,34, 26, 1),
     ("E&I 23,100 m cable, 420 home runs",          24,38, 18, 1),
     ("RT backlog, hydro/He tests, punchlist",      34,41,  9, 0),
     ("cold + hot commissioning",                   39,47, 10, 0),
     ("72 h guarantee run",                         47,48, 10, 0)]
print(f"\n  {'activity':45}{'st':>4}{'fin':>5}{'crew':>6}{'man-wk':>8}")
mw_c=0; mw_all=0
for k,s,f,c,isc in NET:
    mw=(f-s)*c; mw_all+=mw
    if isc: mw_c+=mw
    print(f"  {k:45}{s:4}{f:5}{c:6}{mw:8}")
print(f"  construction man-weeks in network = {mw_c} vs required {need_mw:.0f} "
      f"-> {mw_c/need_mw:.0%}  OK")
T=max(f for _,_,f,_,_ in NET)
Wn=T; hist=[0]*Wn
for k,s,f,c,_ in NET:
    for i in range(s,f): hist[i]+=c
print(f"  duration to first NaCN = {T} weeks = {T/4.33:.1f} months")
print(f"  peak direct crew = {max(hist)}, +30% indirect = {round(max(hist)*1.3)}")
print(f"  footprint check: 1200 m2 / {max(hist)} = {1200/max(hist):.0f} m2 per direct man")
print("  crew: "+"".join(f"{x:4}" for x in hist))
print(f"  total man-weeks on site (all disciplines) = {mw_all} = "
      f"{mw_all*HRS:.0f} paid man-hours")

print()
print("="*78); print("22. CORRECTED HEADLINE TABLE"); print("="*78)
MOD={"rfc":8.8,"first":13.2,"bench":17.6,"peak":27,"mh":165*60,"welds":44,"cable":1900}
STK={"rfc":41,"first":48,"bench":52,"peak":round(max(hist)*1.3),"mh":mw_all*60,
     "welds":157,"cable":23100}
rows=[("site weeks to mechanical completion (RFC)",MOD["rfc"],STK["rfc"]),
      ("site weeks to first on-spec NaCN",MOD["first"],STK["first"]),
      ("site weeks incl. 30-d availability run",MOD["bench"],STK["bench"]),
      ("peak headcount on site (direct+indirect)",MOD["peak"],STK["peak"]),
      ("paid man-hours spent inside the mine lease",MOD["mh"],STK["mh"]),
      ("Category M (HCN) welds made in the field",MOD["welds"],STK["welds"]),
      ("instrument cable pulled in the field, m",MOD["cable"],STK["cable"])]
print(f"  {'':46}{'MODULAR':>10}{'STICK':>10}{'ratio':>8}")
for k,a,b in rows:
    print(f"  {k:46}{a:>10,.0f}{b:>10,.0f}{b/a:>7.1f}x")
print(f"\n  END-TO-END (contract award -> first NaCN):")
print(f"    modular    : 34-40 wk China fab+FAT+sea+clearance (PARALLEL with 24 wk")
print(f"                 client civil) + {MOD['first']:.0f} wk site = 47-53 wk")
print(f"    stick-built: 20-26 wk detail eng.+procurement + {STK['first']} wk site "
      f"= 68-74 wk")
print(f"    end-to-end saving = 21 weeks (~29%). The 3.6x is SITE time only - say so.")
