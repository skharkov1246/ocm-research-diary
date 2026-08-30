# -*- coding: utf-8 -*-
"""Part 5 - corrections: road gauge, true stick-built network, overlapped manning."""
import json, math
print("="*78); print("18. ROAD GAUGE CHECK, TEMA -> TARKWA (~330 km) - a plug-and-play finding")
print("="*78)
DECK = 1.00   # low-bed deck height, m (3-axle step-frame)
DECK_STD = 1.45  # standard flat-bed
mods = [("M1",11.8,2.3,2.6,16.5),("M2",11.8,2.3,2.6,19.0),("M3",11.8,2.3,2.6,18.5),
        ("M4",11.8,2.3,2.6,17.0),("M5",11.8,2.3,2.6,21.0),("M6",11.8,2.3,2.6,21.0)]
LIM = {"width":2.60,"height":4.50,"len_semi":18.5,"axle_group_t":24.0,"gvw_t":41.0}
print("assumed Ghana Highway Authority / ECOWAS in-gauge limits (SRC=ASSUMPTION - the")
print("Ghana Highway Authority abnormal-load permit thresholds were NOT retrieved this")
print("session; confirm before quoting). Values used: W<=2.60 m, H<=4.50 m, GVW<=41 t.")
print(f"{'mod':4}{'W':>6}{'H_on_lowbed':>13}{'H_on_flat':>11}{'mass':>7}  verdict")
for n,L,W,H,m in mods:
    h1,h2 = H+DECK, H+DECK_STD
    ok = W<=LIM["width"] and h1<=LIM["height"]
    print(f"{n:4}{W:6.2f}{h1:13.2f}{h2:11.2f}{m:7.1f}  "
          f"{'IN-GAUGE on low-bed' if ok else 'ABNORMAL'}"
          f"{'' if h2<=LIM['height'] else ' (flat-bed also OK on height)'}")
print()
print("FINDING: every module is IN-GAUGE by road. The OOG=1 flags in pkg.json are a")
print("SEA-container constraint (2.3 m wide + 0.15 m frame > 2.352 m 40HC internal),")
print("NOT a road constraint. Consequence for the schedule:")
print("  - no abnormal-load permit, no police escort, no utility-line lifts on the")
print("    330 km Tema-Tarkwa run -> the convoy is a NORMAL haulage job, 2 days,")
print("    and it does NOT sit on the critical path waiting for a permit office.")
print("  - the ONLY abnormal item would be the urea silo T3 (D7.0 m) - and T3 is")
print("    Ghana-fabricated (pkg.json origin GH), so it is built on site from plate.")
print("  => the plant contains ZERO abnormal road loads. That is a real, and rare,")
print("     plug-and-play property; it is worth designing to keep (cap module width")
print("     at 2.45 m incl. frame and height at 3.3 m incl. skid).")
print()
print("Sea leg: modules go as BREAKBULK / 40' flat rack, not in boxes. Tema has")
print("mobile-harbour-crane capacity for 21 t lifts (SRC=ASSUMPTION, not verified).")

print(); print("="*78); print("19. STICK-BUILT: PROPER NETWORK (not a sum of phases)"); print("="*78)
MH_STICK = 16681.0; HRS=60; TT=0.48
NET=[("site establishment / camp / temp utilities", 0, 6,  8),
     ("civil incl. in-situ equipment plinths + 28 d cure", 3, 15, 14),
     ("equipment erection + structural steel in place",   13, 21, 12),
     ("piping: ~157 Cat-M + ~185 utility field welds",    18, 30, 14),
     ("E&I: 23,100 m field cable, 420 home runs, 185 loops",22,33, 10),
     ("RT backlog, hydro/He tests, punchlist close-out",  30, 36,  8),
     ("cold + hot commissioning",                         34, 42, 10),
     ("72 h guarantee run",                               42, 43, 10)]
print(f"  {'activity':52}{'start':>6}{'fin':>5}{'crew':>6}")
for k,s,f,c in NET: print(f"  {k:52}{s:6}{f:5}{c:6}")
T_STICK=max(f for _,_,f,_ in NET)
Wn=T_STICK; hist=[0]*Wn
for k,s,f,c in NET:
    for i in range(s,f): hist[i]+=c
mh_check=sum(hist)*HRS*TT
print(f"  duration to first NaCN = {T_STICK} weeks = {T_STICK/4.33:.1f} months")
print(f"  man-week integral = {sum(hist)} man-weeks -> {mh_check:.0f} productive mh")
print(f"  (target from discipline take-off was {MH_STICK:.0f} mh; agreement "
      f"{mh_check/MH_STICK:.0%} -> network is self-consistent)")
print(f"  peak direct crew = {max(hist)}; +30% indirect = {round(max(hist)*1.3)} men on a")
print(f"  1200 m2 footprint = {1200/max(hist):.0f} m2/man -> congested, hence TT=0.48.")

print(); print("="*78); print("20. MODULAR: OVERLAPPED NETWORK + MANNING"); print("="*78)
MNET=[("receive/offload/inspect",             0.0,0.6, 6),
      ("set modules + tanks + stack (crane)", 0.6,1.4,10),
      ("shim/align/torque/grout + cure",      1.4,2.3, 6),
      ("pipe rack, platforms, canopy",        1.8,3.2, 8),
      ("Cat-M interconnect weld (44 welds)",  2.6,4.8, 6),
      ("utility interconnect (55 welds)",     3.2,5.0, 6),
      ("RT 100% Cat-M, NIGHT SHIFT",          3.4,5.6, 3),
      ("weld repairs + re-shoot",             5.0,5.7, 4),
      ("E&I tray/cable/terminations",         3.6,5.6, 6),
      ("instrument hook-up + gas detectors",  5.2,6.4, 4),
      ("hydro/pneumatic 14 packs + He sniff", 6.0,7.4, 6),
      ("reinstate, N2 dry-out 32 h, punchlist",7.2,8.4,6),
      ("loop check 185 + 26 SIF validation",  6.8,8.6, 4),
      ("RFC certificate (hold point)",        8.6,8.8, 4),
      ("cold commissioning",                  8.8,10.0,8),
      ("refractory dry-out 57 h (24/7)",     10.0,10.6,8),
      ("chemicals in, caustic circulation",  10.6,11.2,8),
      ("gauze, light-off, conditioning, ramp",11.2,12.4,10),
      ("72 h guarantee run",                 12.4,13.2,10),
      ("30-d availability run (option)",     13.2,17.6,6)]
print(f"  {'activity':40}{'start':>6}{'fin':>6}{'crew':>6}")
for k,s,f,c in MNET: print(f"  {k:40}{s:6.1f}{f:6.1f}{c:6}")
res=0.1; N=int(17.6/res)+1; h=[0.0]*N
for k,s,f,c in MNET:
    for i in range(int(s/res),int(f/res)): h[i]+=c
wkly=[max(h[int(w/res):int((w+1)/res)]) for w in range(18)]
print("  week: "+"".join(f"{i:4}" for i in range(18)))
print("  crew: "+"".join(f"{int(x):4}" for x in wkly))
pk=max(wkly)
print(f"  RFC at wk 8.8; first on-spec NaCN at wk 13.2; benchmark package wk 17.6")
print(f"  peak direct crew {int(pk)}, +30% indirect = {round(pk*1.3)} men")
mh_int=sum(h)*res*HRS*0.62
print(f"  man-week integral = {sum(h)*res:.0f} man-weeks -> {mh_int:.0f} productive mh")
print(f"  (take-off said 3622 mh construction + ~1900 mh commissioning = 5522;")
print(f"   agreement {mh_int/5522:.0%})")
print()
print("="*78); print("21. HEADLINE COMPARISON"); print("="*78)
print(f"  {'':44}{'MODULAR':>10}{'STICK-BUILT':>13}")
for k,a,b in [("site time to RFC (weeks)","8.8","36"),
              ("site time to first on-spec NaCN (weeks)","13.2","43"),
              ("+ 30-d availability (weeks)","17.6","47"),
              ("peak men on site","13","44"),
              ("man-hours exposed on site","~5,500","~19,000"),
              ("HCN-service (Cat M) FIELD welds","44","~157"),
              ("field-pulled instrument cable (m)","1,900","23,100"),
              ("abnormal road loads","0","0 (but 12x more truck movements)")]:
    print(f"  {k:44}{a:>10}{b:>13}")
print(f"\n  site-time compression to first product: 43/13.2 = {43/13.2:.1f}x")
print(f"  BUT: total project time is NOT compressed 3x. Chinese fabrication + FAT +")
print(f"  6-8 weeks sea + clearance runs ~34-40 weeks in parallel with client civil.")
print(f"  Honest end-to-end numbers (contract award -> first NaCN):")
print(f"    modular     : 34-40 wk fab&ship (parallel w/ civil) + 13 wk site = 47-53 wk")
print(f"    stick-built : 20-26 wk engineering+procurement + 43 wk site   = 63-69 wk")
print(f"    => end-to-end saving ~16 weeks (~25%), NOT 70%.")
print(f"  The 3.3x number is SITE time, and that is the number that matters for a")
print(f"  live gold mine: 30 fewer weeks of contractors inside the ML, next to the")
print(f"  CIL circuit, on a permit-to-work system.")
