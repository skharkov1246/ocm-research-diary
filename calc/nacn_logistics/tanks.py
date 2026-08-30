# -*- coding: utf-8 -*-
"""Tanks: 3 sourcing options each. (a) assembled ex-China OOG (b) knocked-down ex-China
+ field erection (c) fully local. Field erection is incurred in (b) and is inside k for (c)."""
import sys, json, math
sys.path.insert(0,"/home/user/ocm-research-diary/calc/nacn_logistics")
from packages import CONT, RATE
P=[tuple(x) for x in json.load(open("/home/user/ocm-research-diary/calc/nacn_logistics/pkg.json"))]
D={p[0]:p for p in P}
m=lambda t:D[t][2]; v=lambda t:D[t][3]*D[t][4]*D[t][5]; c=lambda t:D[t][6]*1000.0
K={"T1":1.30,"T2":1.30,"T3":1.25,"T4":1.30,"T5":1.35,"T6":1.30}
BOX=6300.; OOGC=RATE["ocean_40FR_OOG"]+RATE["dest_port_OOG"]+RATE["inland_OOG_tarkwa"]  # 23100
USE_V=CONT["40HC"][3]*0.72; USE_M=24.1; LEV=0.0245
KD   = 0.18     # knocked-down nested volume / assembled envelope volume. SRC=ASSUMPTION
ERECT= 900.0    # $/t Ghana field erection of a KD tank (weld-out, NDT, test). SRC=ASSUMPTION
KDDISC=0.88     # KD shop cost vs assembled shop cost (no shop weld-out)

print(f"{'tag':4} {'m,t':>5} {'V_ass':>7} {'V_kd':>6} | {'(a) assembled OOG':>18} {'(b) KD + erect':>15} {'(c) local':>10}  best")
res={}
for t in K:
    Vk=v(t)*KD
    a = c(t) + OOGC + (c(t)+OOGC)*LEV
    slots_b = max(Vk/USE_V, m(t)/USE_M)
    fb = slots_b*BOX
    b = c(t)*KDDISC + fb + (c(t)*KDDISC+fb)*LEV + m(t)*ERECT
    cc= c(t)*K[t]
    best=min([(a,"a"),(b,"b"),(cc,"c")])
    res[t]=(a,b,cc,best[1])
    print(f"{t:4} {m(t):5.1f} {v(t):7.1f} {Vk:6.1f} | {a:18,.0f} {b:15,.0f} {cc:10,.0f}   {best[1]}"
          f"   (b vs c: {b-cc:+,.0f})")
print()
print(f"  Sum (a) all assembled OOG : ${sum(r[0] for r in res.values()):,.0f}")
print(f"  Sum (b) all KD ex-China   : ${sum(r[1] for r in res.values()):,.0f}")
print(f"  Sum (c) all local Ghana   : ${sum(r[2] for r in res.values()):,.0f}")
print(f"  Sum best-of              : ${sum(min(r[0],r[1],r[2]) for r in res.values()):,.0f}")
print(f"  a->best saving           : ${sum(r[0] for r in res.values())-sum(min(r[0],r[1],r[2]) for r in res.values()):,.0f}")

print("\n  break-even Ghana premium k where local beats KD-ex-China:")
for t in K:
    a,b,cc,_=res[t]
    print(f"    {t:4} k* = {b/c(t):5.2f}   (assumed {K[t]:.2f}) -> {'LOCAL wins' if K[t]<b/c(t) else 'import KD wins'}")

print("\n  NOTE: option (b) already puts the welding, NDT and hydrotest labour in Ghana.")
print("  The real question is not 'local vs import' but 'which HALF of the tank is local':")
print("  plate+rolling ex-China, erection+testing in Ghana is the dominant answer for T1/T2/T5,")
print("  while T3/T4 (big, cheap, non-hazardous) are fully local.")
