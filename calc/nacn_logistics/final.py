# -*- coding: utf-8 -*-
import sys, json, math
sys.path.insert(0,"/home/user/ocm-research-diary/calc/nacn_logistics")
from packages import CONT, RATE
P=[tuple(x) for x in json.load(open("/home/user/ocm-research-diary/calc/nacn_logistics/pkg.json"))]
D={p[0]:p for p in P}
m=lambda t:D[t][2]; v=lambda t:D[t][3]*D[t][4]*D[t][5]; c=lambda t:D[t][6]*1000.0
K={"S1":1.44,"S2":1.44,"S3":1.44,"S4":1.44,"S5":1.44,"T1":1.30,"T2":1.30,"T3":1.25,
   "T4":1.30,"T5":1.35,"T6":1.30,"PP3":1.44,"PP4":1.44,"E3":1.10,"E4":1.48,"E5":1.15}
BOX=6300.; OOGC=RATE["ocean_40FR_OOG"]+RATE["dest_port_OOG"]+RATE["inland_OOG_tarkwa"]
USE_V=CONT["40HC"][3]*0.72; USE_M=24.1; LEV=0.0245

print("=== THE DECISION RULE: value density of shipped volume ===")
for k in (1.15,1.30,1.44,1.60):
    psi=BOX*(1+LEV)/(USE_V*(k-1)); phi=BOX*(1+LEV)/(USE_M*(k-1))
    print(f"  local premium k={k:4.2f}: LOCALISE if FOB value < ${psi:6.0f}/m3 shipped  "
          f"(or < ${phi:6.0f}/t for dense items)")
print("  -> at k=1.44 the line is $267/m3: empty vessels & silos fall below it, "
      "fabricated steel ($550/m3) and cable ($6000/m3) do not.\n")
for t in sorted(K,key=lambda x:c(x)/v(x)):
    print(f"    {t:4} ${c(t)/v(t):7.0f}/m3  ${c(t)/m(t):7.0f}/t   {D[t][1][:46]}")

print("\n=== SCENARIOS ===")
def box_count(mass,vol): return math.ceil(max(vol/USE_V, mass/USE_M))
LOC_OPT=["T1","T3","T4","T5"]
rows=[]
# A baseline
rows.append(("A  naive turnkey ex-China, OOG modules+tanks", 9, 11, 435444, 0))
# C  design-to-gauge + ISO frames + void-fill + urea/1-day scope, import all fabricables
ship=[t for t in D if t not in ("C1","C2","Z2","G2")]
loose=[t for t in ship if not t.startswith("M")]
# void-fill consumed 25.5 t / 72.9 m3 of the dense CN pool
lm=sum(m(t) for t in loose)+8.0-25.5; lv=sum(v(t) for t in loose)+22.0-72.9
nC=box_count(lm,lv)
fobC=sum(c(t) for t in ship)+c("G2")+c("C2")*0.6
frC=6*(BOX+800)+nC*BOX+2804
logC=frC+(fobC+frC)*RATE["insurance_pct"]; logC+= (fobC+logC)*LEV
rows.append((f"C  in-gauge + ISO frames + void-fill, all fab imported ({lm:.0f} t/{lv:.0f} m3 loose)",0,nC,logC,0))
# D  = C + localise ONLY T1,T3,T4,T5
loose2=[t for t in loose if t not in LOC_OPT]
lm2=sum(m(t) for t in loose2)+8.0-25.5; lv2=sum(v(t) for t in loose2)+22.0-72.9
nD=box_count(lm2,lv2)
fobD=fobC-sum(c(t) for t in LOC_OPT)
frD=6*(BOX+800)+nD*BOX+2804
logD=frD+(fobD+frD)*RATE["insurance_pct"]; logD+=(fobD+logD)*LEV
locD=sum(c(t)*K[t] for t in LOC_OPT)
rows.append((f"D  = C + localise T1,T3,T4,T5 only ({lm2:.0f} t/{lv2:.0f} m3 loose)",0,nD,logD,locD-sum(c(t) for t in LOC_OPT)))
# E  = C + blanket localisation of every Ghana-capable item
loose3=[t for t in loose if t not in K]
lm3=sum(m(t) for t in loose3)+8.0-25.5; lv3=sum(v(t) for t in loose3)+22.0-72.9
nE=box_count(max(lm3,0),max(lv3,0))
fobE=fobC-sum(c(t) for t in K)
frE=6*(BOX+800)+nE*BOX+2804
logE=frE+(fobE+frE)*RATE["insurance_pct"]; logE+=(fobE+logE)*LEV
locE=sum(c(t)*(K[t]-1) for t in K)
rows.append((f"E  = C + blanket localisation ({lm3:.0f} t/{lv3:.0f} m3 loose)",0,nE,logE,locE))

print(f"{'scenario':66} {'OOG':>4} {'box':>4} {'units':>6} {'logistics$':>11} {'loc.premium$':>13} {'TOTAL$':>11}")
base=None
for n,o,b,lg,pr in rows:
    u=o+b+(0 if n.startswith("A") else 6)
    if n.startswith("A"): u=o+b
    tot=lg+pr
    if base is None: base=tot
    print(f"{n:66} {o:>4} {b:>4} {u:>6} {lg:>11,.0f} {pr:>13,.0f} {tot:>11,.0f}   "
          f"{'' if base==tot else f'-{base-tot:,.0f} ({100*(base-tot)/base:.0f}%)'}")

print("\n=== SHIPPED VOLUME / MASS ===")
tot_v=sum(v(t) for t in D); tot_m=sum(m(t) for t in D)
print(f"  A: 9 OOG flat racks + 11 x 40'HC = 20 units;  ~{sum(m(t) for t in D if t!='C1'):.0f} t / {sum(v(t) for t in D):.0f} m3 crossing the ocean")
print(f"  D: 6 ISO-frame modules + {nD} x 40'HC = {6+nD} units;  {113+25.5:.0f} t in modules + {lm2:.0f} t loose")
print(f"     -> transport places  20 -> {6+nD}  ({100*(1-(6+nD)/20):.0f}% fewer)")
print(f"     -> heavy-truck movements Tema->Tarkwa  20 -> {6+nD} (abnormal-load convoys 9 -> 0)")
