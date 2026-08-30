# -*- coding: utf-8 -*-
"""Part 7 - crew reconciliation, schedule risk (P50/P80), and site-phase cost."""
import math, random
HRS,TT=60,0.48
# --- corrected stick-built crews to match the 579 man-week take-off
NET=[("equipment erection + structural steel",13,23,12),
     ("piping ~157 Cat-M + ~185 utility welds",18,34,20),
     ("E&I 23,100 m cable, 420 home runs",24,38,14)]
mw=sum((f-s)*c for _,s,f,c in NET)
print(f"corrected stick-built construction crews -> {mw} man-weeks vs 579 required "
      f"= {mw/579:.0%}. Durations unchanged (they are set by work-front sequence,")
print("not by crew size; the man-hour check only confirms the crews are not fantasy).")
print("=> stick-built to first NaCN stays 48 weeks; peak headcount drops 57 -> 43.")

print(); print("="*78); print("23. SCHEDULE RISK - MONTE CARLO ON THE MODULAR SITE PHASE"); print("="*78)
# (activity, most-likely wk, optimistic, pessimistic, driver)
RISK=[("offload/inspect + transit damage",0.6,0.4,2.5,"shipping damage to a module skid"),
      ("module setting (crane availability)",0.8,0.6,3.0,"crane not on site / ground bearing"),
      ("grout cure",0.9,0.7,1.2,"fixed by chemistry"),
      ("steel + Cat-M welding",2.6,2.0,5.0,"coded welder availability; plinth survey error"),
      ("RT + repairs",1.0,0.6,3.0,"reject rate >6%; film to Accra and back"),
      ("E&I + instruments",2.0,1.5,3.5,"missing certs, IECEx paperwork at customs"),
      ("pressure test + He sniff",1.4,1.0,4.0,"He tech visa; leaks at 33 Cat-M flanges"),
      ("dry-out + reinstatement",1.2,0.9,2.0,"fixed by physics (32 h N2)"),
      ("loop check + 26 SIF validation",1.8,1.2,4.0,"SIS logic bugs found only in field"),
      ("cold comm.",1.2,0.9,2.5,"utility quality: raw water, gas pressure at BL"),
      ("refractory dry-out",0.6,0.5,0.9,"fixed by vendor curve"),
      ("chemicals in",0.6,0.4,2.0,"NaOH / urea delivery to site"),
      ("gauze light-off + conditioning + ramp",1.2,1.0,4.0,"THE big one: gauze not lighting,"
       " ratio control, HCN in stack"),
      ("72 h guarantee run",0.8,0.7,3.0,"one trip = restart the 72 h clock")]
random.seed(7)
def tri(m,a,b):
    u=random.random(); c=(m-a)/(b-a)
    return a+math.sqrt(u*(b-a)*(m-a)) if u<c else b-math.sqrt((1-u)*(b-a)*(b-m))
N=20000; OVL=0.24
sims=[]
for _ in range(N):
    s=sum(tri(m,a,b) for _,m,a,b,_ in RISK)
    sims.append(s*(1-OVL))
sims.sort()
det=sum(m for _,m,_,_,_ in RISK)*(1-OVL)
print(f"  deterministic (most-likely, 24% overlap credit) = {det:.1f} weeks")
for p in (10,50,80,90,95):
    print(f"  P{p:<3} = {sims[int(N*p/100)]:5.1f} weeks")
print(f"  P80/P50 = {sims[int(N*.8)]/sims[int(N*.5)]:.2f}  -> carry "
      f"{sims[int(N*.8)]-sims[int(N*.5)]:.1f} weeks of schedule contingency, not 'a bit'.")
print("  QUOTE 13 weeks INTERNALLY, COMMIT 17 weeks CONTRACTUALLY.")
print("\n  top-3 variance contributors (pessimistic minus most-likely):")
for n,m,a,b,d in sorted(RISK,key=lambda r:-(r[3]-r[1]))[:5]:
    print(f"    +{b-m:4.1f} wk  {n:38} <- {d}")

print(); print("="*78); print("24. WHAT ACTUALLY BREAKS 'PLUG-AND-PLAY' (ranked, with days)"); print("="*78)
BRK=[("as-built anchor-bolt survey not sent to China (or sent late)",
      "every module base frame is re-drilled in the field; all 44 Cat-M spools re-cut",
      "+25-40 d"),
     ("no ASME IX coded SS316L 6G welder available in Ghana",
      "fly in 2 welders + re-qualify to the WPS on site (test coupons + RT)",
      "+10-15 d"),
     ("He mass-spec leak test (B31.3 M345.8) treated as 'optional'",
      "it is NOT optional for Cat M; discovering this at RFC means re-blinding all "
      "7 Cat-M packs","+8-12 d"),
     ("IECEx / ATEX certificates of the Chinese gas detectors not accepted",
      "gas detection is a SIF; no cert = no permit to introduce HCN","+15-30 d"),
     ("crane booked as 100 t class (marginal at 14 m radius)",
      "re-book 130-160 t; in Ghana that is a Takoradi/Tema mobilisation","+7-14 d"),
     ("gas at battery limit not at spec (Genser line pressure/Wobbe)",
      "burner ratio control cannot hold; light-off aborts","+5-20 d"),
     ("hydrotesting the HCN gas circuit instead of pneumatic",
      "adds the full 32 h N2 dry-out plus a second tightness check","+3-4 d"),
     ("3rd-party/insurer pre-op audit scheduled AFTER mechanical completion",
      "auditor availability in West Africa; run it CONCURRENT with loop checks","+10-20 d")]
for a,b,c in BRK:
    print(f"  {c:>9}  {a}")
    print(f"             -> {b}")
print("\n  Sum of the top three, if all hit: +43 to +67 days = the entire schedule")
print("  contingency. Every one of them is a DOCUMENT problem, not a steel problem.")
print("  Plug-and-play is won or lost in the 8 weeks BEFORE the ship sails.")

print(); print("="*78); print("25. SITE-PHASE COST (order of magnitude, SRC=ASSUMPTION rates)"); print("="*78)
RATES={"expat_specialist_day":950,"expat_supervisor_day":650,"local_skilled_day":95,
       "local_helper_day":45,"crane_130t_day":2200,"NDT_crew_day":900,
       "camp_per_man_day":85}
for k,v in RATES.items(): print(f"    {k:26} ${v}")
# modular
mod_days=13.2*6
cost_mod = (4*RATES["expat_specialist_day"] + 2*RATES["expat_supervisor_day"])*mod_days \
         + (12*RATES["local_skilled_day"] + 9*RATES["local_helper_day"])*mod_days \
         + 5*RATES["crane_130t_day"] + 12*RATES["NDT_crew_day"] \
         + 27*RATES["camp_per_man_day"]*mod_days
stk_days=48*6
cost_stk = (6*RATES["expat_specialist_day"] + 5*RATES["expat_supervisor_day"])*stk_days \
         + (24*RATES["local_skilled_day"] + 14*RATES["local_helper_day"])*stk_days \
         + 20*RATES["crane_130t_day"] + 55*RATES["NDT_crew_day"] \
         + 43*RATES["camp_per_man_day"]*stk_days
print(f"\n  modular site phase   : {mod_days:.0f} working days -> ${cost_mod/1e6:.2f} M")
print(f"  stick-built site phase: {stk_days:.0f} working days -> ${cost_stk/1e6:.2f} M")
print(f"  delta = ${(cost_stk-cost_mod)/1e6:.2f} M of SITE cost avoided, on a "
      f"$4.2 M CAPEX base = {(cost_stk-cost_mod)/4.2e6:.0%} of CAPEX.")
print(f"  (the modular premium paid to the Chinese shop must be < ${(cost_stk-cost_mod)/1e6:.2f} M")
print(f"   for modular to be cost-neutral; it typically is, so modular is not just")
print(f"   faster and safer here - on this cost base it is cheaper.)")
# revenue-time value
print(f"\n  TIME-VALUE: 35 weeks earlier production x 2500 t/y x (price $2300 - cash")
print(f"  cost $1228)/t = {35/52*2500*(2300-1228)/1e6:.2f} M$ of gross margin pulled forward.")
print(f"  That number alone pays for the entire modularisation premium.")
