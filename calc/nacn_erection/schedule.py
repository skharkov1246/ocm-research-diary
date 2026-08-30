# -*- coding: utf-8 -*-
"""Part 4 - week-by-week schedule, crewing, and stick-built comparison."""
import json, math
E=json.load(open("/home/user/ocm-research-diary/calc/nacn_erection/erect_out.json"))
C=json.load(open("/home/user/ocm-research-diary/calc/nacn_erection/comm_out.json"))

HRS_WK = 6*10          # 6-day, 10-h roster (mine-site norm)
TT_MOD = 0.62          # tool-time factor, modular job (small crew, pre-fab, low scaffold)
TT_STK = 0.48          # tool-time factor, stick-built in West Africa
print("="*78)
print("13. MANHOUR -> DURATION CONVERSION")
print("="*78)
print(f"roster {HRS_WK} paid h/person/week; tool-time factor modular {TT_MOD}, stick-built {TT_STK}")
print(f"=> effective productive h/person/week: modular {HRS_WK*TT_MOD:.0f}, stick {HRS_WK*TT_STK:.0f}")

# ---- activities: (name, weeks_start, weeks_dur, crew, note) built from mh above
def wk(mh, crew, tt=TT_MOD): return mh/(crew*HRS_WK*tt)

print()
print("="*78); print("14. SCHEDULE - MODULAR (t=0 is the day the last module lands on site)"); print("="*78)
NEG = [
 ("-24","geotech, site survey, permit-to-construct in hand",            "client"),
 ("-20","earthworks, laterite platform, haul road upgrade Tema->Tarkwa","client civil, 8"),
 ("-16","piling/blinding, rebar, formwork, C1 193 m3 + C2 35 m3 pours", "client civil, 12"),
 ("-10","HDPE liner + bund + sump (ICMI Std 4), anchor bolts to template","client civil, 10"),
 (" -8","28-d cure complete; AS-BUILT anchor-bolt survey -> to China",  "surveyor, 2"),
 (" -6","BL tie-ins stubbed & valved: gas, HV, CW, product to leach",   "client, 6"),
 (" -4","modules leave Tema port; abnormal-load permits, 330 km convoy","logistics, 4"),
]
for a,b,c in NEG: print(f"  wk {a:>3}  {b:58} [{c}]")
print("  ^ ALL of the above is client scope and OFF the vendor critical path.")
print("    The single hard interlock: as-built anchor-bolt survey must reach the")
print("    Chinese shop >=8 weeks before shipment, or module base frames will not")
print("    match the plinths and every one of the 44 Cat-M spools is a field re-cut.")

ACT = [
 # name, dur_wk, crew, discipline note, depends
 ("Receive, offload, inspect, damage survey, lay-down",      0.6, 6,  "riggers+QC"),
 ("Set 6 modules + tanks/stack (crane on site 5 d)",         0.8, 10, "crane crew + surveyor"),
 ("Shim, align, anchor-bolt torque, grout + 3 d cure",       0.9, 6,  "millwrights"),
 ("Pipe rack S1 / platforms S2 / frames S3 erect",           1.2, 8,  "steel erectors"),
 ("Cat-M interconnect: fit-up + weld (44 welds, 138 DI)",    2.2, 6,  "2 coded welders + 2 fitters + 2 helpers"),
 ("Utility interconnect + BL tie-ins (55 welds)",            1.8, 6,  "pipefitters"),
 ("100% RT of Cat-M welds (night shift) + film review",      1.4, 3,  "NDT 2 + RT level II"),
 ("Weld repairs (6% reject) + re-shoot",                     0.7, 4,  "welders"),
 ("E&I: tray, 1900 m trunk cable, 760 cores terminated",     1.8, 6,  "electricians"),
 ("Instrument hook-up, gas-detector siting + bump test",     1.0, 4,  "instrument techs"),
 ("Hydro/pneumatic test, 14 packs (2/day) + He sniff Cat-M", 1.4, 6,  "test crew + client witness"),
 ("Reinstatement, N2 dry-out 32 h, blind removal, punchlist",1.0, 6,  "mech"),
 ("Loop check 185 loops + 26 SIF end-to-end (IEC 61511 cl.15)",1.2, 4, "DCS/SIS engineers"),
 ("MECHANICAL COMPLETION / RFC certificate",                 0.2, 4,  "client + vendor + 3rd party"),
 ("Cold commissioning: water/air runs, pumps, blower, deluge",1.2, 8,  "commissioning + ops trainees"),
 ("Refractory dry-out 57 h continuous + oxidiser on fuel gas",0.5, 8, "24/7, 4 per shift"),
 ("Chemicals in: NaOH, urea; caustic circulation, pH>11.5",  0.6, 8,  "24/7"),
 ("Gauze install, light-off, conditioning 48 h, ramp to 100%",1.1, 10, "24/7 + vendor process eng"),
 ("72 h guarantee run + emissions/HCN monitoring",           0.7, 10, "24/7 + client witness"),
 ("30-d availability run (SAFETY-BENCHMARK option, overlaps ops)",4.3,6,"client ops + 2 vendor"),
]
t=0.0; rows=[]
print()
print(f"  {'wk':>5} {'act':62}{'dur':>5}{'crew':>5}")
for n,d,c,note in ACT:
    rows.append((t,d,c,n,note))
    print(f"  {t:5.1f} {n:62}{d:5.1f}{c:5}")
    t+=d
T_MC   = sum(a[1] for a in ACT[:14])
T_RFSU = sum(a[1] for a in ACT[:19])
T_ALL  = t
print("-"*78)
print(f"  Mechanical completion (RFC)      : {T_MC:5.1f} weeks after last module lands")
print(f"  Ready for start-up -> 72h guarantee passed : {T_RFSU:5.1f} weeks")
print(f"  incl. 30-day availability run    : {T_ALL:5.1f} weeks")

# overlap credit: several activities run in parallel on a small site
OVL = 0.24
print(f"\n  The list above is SEQUENTIAL. Real overlap (steel//E&I tray, RT at night,")
print(f"  utility welding // Cat-M welding, cold comm. of utilities // loop checks):")
print(f"  credit {OVL:.0%} -> RFC {T_MC*(1-OVL):.1f} wk, first NaCN {T_RFSU*(1-OVL):.1f} wk,")
print(f"  guarantee+availability {T_ALL*(1-OVL):.1f} wk.")
T_MC_R, T_RF_R, T_ALL_R = T_MC*(1-OVL), T_RFSU*(1-OVL), T_ALL*(1-OVL)
print(f"\n  HEADLINE: last module on the plinth -> first on-spec NaCN in "
      f"{T_RF_R:.0f} weeks (~{T_RF_R/4.33:.1f} months).")
print(f"           full benchmark package (30-d availability) {T_ALL_R:.0f} weeks.")

# ---- crew histogram
print()
print("="*78); print("15. MANNING PROFILE"); print("="*78)
W=int(math.ceil(t)); hist=[0]*W
for st,d,c,n,note in rows:
    for i in range(int(st), int(math.ceil(st+d))):
        if i<W: hist[i]=max(hist[i],c)
peak=max(hist)
print("  week: " + "".join(f"{i:3}" for i in range(W)))
print("  crew: " + "".join(f"{h:3}" for h in hist))
print(f"  peak site crew (modular) = {peak} incl. supervision -> add 30% indirect "
      f"(HSE, QA/QC, storeman, driver, medic) = {round(peak*1.3)}")
print("  work-front limit check: process footprint ~40x30 m = 1200 m2; at 30 m2/man")
print(f"  the physical limit is ~{1200//30} men. Peak {round(peak*1.3)} is inside it. "
      "A stick-built job")
print("  would need 60-90 men on the same 1200 m2 -> congestion, which is exactly")
print("  where the stick-built productivity factor 0.48 comes from.")

# ---- KEY ROLES
print()
print("  ROLES THAT MUST BE PRESENT (and are the actual scarce resource):")
for r in [
 "2 x ASME IX coded welders, SS316L GTAW 6G, qualified for the WPS in use",
 "1 x ASNT/ISO 9712 Level II RT interpreter (film review is a HOLD POINT)",
 "1 x He mass-spec leak-test technician (M345.8 - rare in Ghana, likely fly-in)",
 "1 x rigging engineer / appointed person for the 30 t lift plan",
 "1 x SIS engineer for 26 SIF validations (IEC 61511 cl.15 - cannot be a fitter)",
 "1 x Andrussow process engineer from the licensor for gauze light-off",
 "1 x HCN-qualified medic + on-site antidote kit from the first HCN-service test",
]: print("    - "+r)

# ---- STICK-BUILT COMPARISON
print()
print("="*78); print("16. STICK-BUILT COMPARISON (same plant, built in place at Tarkwa)"); print("="*78)
# what modularisation moves to the shop
SHOP_FRAC = {"piping":0.72,"equipment set/align":0.55,"E&I":0.82,"structural":0.60}
mh_field_mod = E["MH_CONSTR"]
print(f"modular FIELD man-hours (calculated in erect.py) = {mh_field_mod:.0f} mh")
# reconstruct total install hours = field + what the shop did
mh_pipe_tot   = E["MH_PIPE"]/(1-SHOP_FRAC["piping"])
mh_ei_tot     = E["mh_ei_tot"]/(1-SHOP_FRAC["E&I"])
mh_oth_tot    = E["mh_oth"]/(1-SHOP_FRAC["structural"])
mh_install_tot= mh_pipe_tot+mh_ei_tot+mh_oth_tot
print(f"  piping   : field {E['MH_PIPE']:.0f} mh is {1-SHOP_FRAC['piping']:.0%} of total "
      f"-> total install {mh_pipe_tot:.0f} mh")
print(f"  E&I      : field {E['mh_ei_tot']:.0f} mh is {1-SHOP_FRAC['E&I']:.0%} of total "
      f"-> total install {mh_ei_tot:.0f} mh")
print(f"  mech/steel: field {E['mh_oth']:.0f} mh is {1-SHOP_FRAC['structural']:.0%} of total "
      f"-> total install {mh_oth_tot:.0f} mh")
print(f"  TOTAL installation content of this plant ~ {mh_install_tot:.0f} mh")
PF = 1.0/ (TT_STK/TT_MOD)
mh_stick = mh_install_tot*(TT_MOD/TT_STK)
print(f"\nstick-built does ALL of it in the field at the field tool-time {TT_STK}:")
print(f"  equivalent field man-hours = {mh_install_tot:.0f} x ({TT_MOD}/{TT_STK}) "
      f"= {mh_stick:.0f} mh")
print(f"  ratio field(stick)/field(modular) = {mh_stick/mh_field_mod:.1f}x")
for crew in (25,40,55):
    d = mh_stick/(crew*HRS_WK*TT_STK)
    print(f"   at {crew} direct men: {d:5.1f} weeks of construction alone")
crew_s=40; d_constr=mh_stick/(crew_s*HRS_WK*TT_STK)
STICK = [("site establishment, camp, temp power/water, laydown",6),
         ("civil (same 228 m3 + more: equipment plinths in situ)",10),
         ("structural steel + equipment erection in place",       8),
         ("piping (all 99+ field welds become ~340 field welds)", d_constr*0.42),
         ("E&I (23,100 m of field-pulled cable, 420 home runs)",  d_constr*0.30),
         ("testing, RT backlog, punchlist",                        5),
         ("cold+hot commissioning (same physics, longer punchlist)",7),
         ("guarantee run",                                         1)]
ts=sum(v for _,v in STICK)
print(f"\n  stick-built phase build-up (40 direct men):")
for k,v in STICK: print(f"    {k:56}{v:6.1f} wk")
print(f"    {'TOTAL, first NaCN':56}{ts:6.1f} wk = {ts/4.33:.1f} months")
print(f"\n  MODULAR   : {T_RF_R:5.1f} wk on site (+ {24} wk of client civil done in parallel")
print(f"              with Chinese fabrication - NOT additive)")
print(f"  STICK-BUILT:{ts:5.1f} wk on site")
print(f"  SITE-TIME COMPRESSION = {ts/T_RF_R:.1f}x  ({ts-T_RF_R:.0f} weeks removed from site)")
print(f"  MAN-WEEKS ON SITE: modular {mh_field_mod/HRS_WK:.0f} vs stick {mh_stick/HRS_WK:.0f} "
      f"-> {1-mh_field_mod/mh_stick:.0%} fewer man-hours exposed to the hazard/climate/site")

print()
print("="*78); print("17. WHAT THIS BUYS IN SAFETY TERMS (the client's actual priority)"); print("="*78)
exp_mod=mh_field_mod; exp_stk=mh_stick
print(f"  exposure hours on a live mine site: {exp_mod:.0f} vs {exp_stk:.0f} mh "
      f"({exp_stk/exp_mod:.1f}x)")
print(f"  at a TRIFR of 3.0 per 200,000 h (typical W-African contractor mining norm,")
print(f"  SRC=ASSUMPTION, not verified): expected recordables "
      f"{exp_mod*3.0/2e5:.2f} vs {exp_stk*3.0/2e5:.2f}")
print(f"  field welds in HCN service: {44} (modular) vs ~{int(44/(1-SHOP_FRAC['piping'])):d} "
      f"(stick-built)")
print(f"  every field weld avoided = one fewer lethal-service leak path made under a")
print(f"  tarpaulin in 32 C / 85% RH instead of in a climate-controlled shop with an")
print(f"  automatic orbital welder and a permanent RT bay.")
print(f"  THIS is the safety argument for modular, and it is quantitative:")
print(f"    {int(44/(1-SHOP_FRAC['piping']))-44} Cat-M field joints eliminated.")
json.dump({"T_MC":T_MC_R,"T_RFSU":T_RF_R,"T_ALL":T_ALL_R,"t_stick":ts,
  "mh_stick":mh_stick,"mh_mod":mh_field_mod,"peak":round(peak*1.3)},
  open("/home/user/ocm-research-diary/calc/nacn_erection/sched_out.json","w"),indent=1)
