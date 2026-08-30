# -*- coding: utf-8 -*-
"""Scenario A (naive turnkey ex-China) vs Scenario B (logistics-minimised)."""
import sys, json, math
sys.path.insert(0,"/home/user/ocm-research-diary/calc/nacn_logistics")
from packages import CONT, RATE, GH, CN

P=[tuple(x) for x in json.load(open("/home/user/ocm-research-diary/calc/nacn_logistics/pkg.json"))]
D={p[0]:p for p in P}
def m(t): return D[t][2]
def v(t): return D[t][3]*D[t][4]*D[t][5]
def c(t): return D[t][6]

# --- local-price ratio k = Ghana delivered-to-site / China FOB (same scope) ---
K = {  # SRC=ASSUMPTION, needs 3 quotes each
 "S1":1.44,"S2":1.44,"S3":1.44,"S4":1.44,"S5":1.44,          # structural steel
 "T1":1.30,"T2":1.30,"T3":1.25,"T4":1.30,"T5":1.35,"T6":1.30, # tanks
 "PP3":1.44,"PP4":1.44,                                        # low-class pipe
 "E3":1.10,"E4":1.48,"E5":1.15,                                # cable / tray / genset
 "C1":1.00,"C2":1.00,"Z2":1.00,
}
NEVER_SHIPPED = {"C1"}          # concrete/rebar/roads: local in every scenario
C2_SHIP_T, C2_SHIP_V = 8.0, 22.0   # HDPE liner rolls + geotextile only

ROAD_CAP_T   = 28.0   # self-imposed road cap Tarkwa (ISO limit 30.48 t) - VERIFY with GHA
BOX_PAYLOAD  = CONT["40HC"][4]
BOX_VOL      = CONT["40HC"][3]
STOW_LOOSE   = 0.72   # crated general cargo stowage efficiency in a box
STOW_VOID    = 0.34   # fraction of a module's free volume actually usable for void-fill crates
CRATE_DENS   = 620.0  # kg/m3 of dense crated goods (valves, spools, spares)

def levies(cif, exempt):
    """Ghana import charges. SRC=UNVERIFIED this session; structure per ECOWAS CET + GRA levies."""
    if exempt:  # Minerals and Mining Act 2006 (Act 703) s.30 'mining list' - UNCONFIRMED applicability
        return cif*(0.005+0.002+0.0075+0.01)      # ECOWAS 0.5 + AU 0.2 + EXIM 0.75 + inspection 1.0
    return cif*(0.05+0.005+0.002+0.0075+0.01)     # + 5% CET duty band for plant

def report(name, rows, total):
    print(f"\n--- {name} ---")
    for r in rows: print("   "+r)
    print(f"   {'TOTAL LOGISTICS':<48} ${total:>10,.0f}")

# =====================================================================
# SCENARIO A : naive turnkey. Vendor ships everything; modules built to
# his default 12.2 x 2.9 x 3.2 m envelope; tanks shipped assembled.
# =====================================================================
def scenario_A():
    rows=[]; cost=0.0
    # 6 modules, over-width (2.9>2.438) and over-height (3.2>2.591) -> OOG flat rack each
    n_oog_mod = 6
    # assembled tanks T1 (D6.0), T3 (D7.2), T4 (D8.0) -> OOG; T2,T5,T6 in-gauge
    oog_tanks = ["T1","T3","T4"]
    n_oog = n_oog_mod + len(oog_tanks)
    c_oog = n_oog*(RATE["ocean_40FR_OOG"]+RATE["dest_port_OOG"]+RATE["inland_OOG_tarkwa"])
    rows.append(f"{'OOG flat racks (6 modules + 3 tanks)':<48} {n_oog:>3} x $23,100 = ${c_oog:>9,.0f}")
    cost += c_oog
    # everything else loose in boxes
    loose = [t for t in D if t not in NEVER_SHIPPED and t not in oog_tanks and not t.startswith("M")]
    lv = sum(v(t) for t in loose) + C2_SHIP_V
    lm = sum(m(t) for t in loose) + C2_SHIP_T
    # C1 excluded; C2 partial -> subtract full C2 then add shipped part
    lv -= v("C2"); lm -= m("C2")
    by_vol  = lv/(BOX_VOL*STOW_LOOSE)
    by_mass = lm/min(BOX_PAYLOAD, ROAD_CAP_T-3.9)
    n_box = math.ceil(max(by_vol, by_mass))
    c_box = n_box*(RATE["ocean_40HC"]+RATE["dest_port_40"]+RATE["inland_40_tarkwa"])
    rows.append(f"{'loose cargo: '+f'{lm:.0f} t / {lv:.0f} m3':<48} {n_box:>3} x  $6,300 = ${c_box:>9,.0f}")
    rows.append(f"{'   (volume-driven '+f'{by_vol:.1f}'+' vs mass-driven '+f'{by_mass:.1f}'+' boxes)':<48}")
    cost += c_box
    fob = sum(c(t) for t in D if t not in NEVER_SHIPPED)*1000 - c("C2")*1000*0.6
    ins = (fob+cost)*RATE["insurance_pct"]; cost += ins
    rows.append(f"{'marine insurance 1.1% of CIF':<48}              ${ins:>9,.0f}")
    for ex in (True,False):
        lv_ = levies(fob+cost, ex)
        rows.append(f"{'import levies (mining-list exempt='+str(ex)+')':<48}              ${lv_:>9,.0f}")
    cost_ex  = cost + levies(fob+cost, True)
    cost_nex = cost + levies(fob+cost, False)
    return rows, cost_ex, cost_nex, fob, n_oog, n_box

rowsA, A_ex, A_nex, fobA, noogA, nboxA = scenario_A()
report("SCENARIO A  naive turnkey ex-China", rowsA, A_ex)
print(f"   {'(non-exempt variant)':<48} ${A_nex:>10,.0f}")
print(f"   {'FOB value shipped':<48} ${fobA:>10,.0f}")

# =====================================================================
# SCENARIO B : logistics-minimised
#   B1 modules re-designed to 40'HC ISO gauge, built ON ISO corner-casting
#      frames (the frame IS the container: no box hire, no tare, no stuffing)
#   B2 everything Ghana-fabricable procured in Ghana
#   B3 module voids filled with dense China-origin cargo
#   B4 zero out-of-gauge pieces
# =====================================================================
MOD_EXT = (12.192, 2.438, 2.896)          # 40'HC ISO envelope
MOD_INT_VOL = 11.90*2.30*2.75             # usable inside the frame, m3
DISP = {"M1":0.42,"M2":0.45,"M3":0.45,"M4":0.50,"M5":0.85,"M6":0.48}
CN_ITEMS = [t for t in D if D[t][7]=="CN"]
GH_ITEMS = [t for t in D if D[t][7]=="GH"]
VOIDFILL_POOL = ["PP2","E1","PP1","E2","Z1","G1"]   # densest first

print("\n--- SCENARIO B  logistics-minimised ---")
print(f"   module usable internal volume {MOD_INT_VOL:.1f} m3 (40'HC ISO frame)")
remaining = {t:[m(t), v(t)] for t in VOIDFILL_POOL}
filled_t = filled_v = 0.0
for mod in ["M1","M2","M3","M4","M5","M6"]:
    free_v = (1-DISP[mod])*MOD_INT_VOL*STOW_VOID
    head_t = ROAD_CAP_T - m(mod)
    got_t=got_v=0.0
    for t in VOIDFILL_POOL:
        if remaining[t][0]<=1e-6: continue
        take_v = min(free_v-got_v, remaining[t][1])
        if take_v<=0.05: continue
        dens = remaining[t][0]/remaining[t][1]
        take_m = min(take_v*dens, head_t-got_t)
        take_v = take_m/dens
        if take_m<=0.02: continue
        remaining[t][0]-=take_m; remaining[t][1]-=take_v
        got_t+=take_m; got_v+=take_v
        if got_v>=free_v-0.05 or got_t>=head_t-0.05: break
    filled_t+=got_t; filled_v+=got_v
    print(f"   {mod}: equip {m(mod):5.1f} t, free void {free_v:5.1f} m3, "
          f"road headroom {head_t:5.1f} t -> void-fill {got_t:5.1f} t / {got_v:5.1f} m3 "
          f"-> gross {m(mod)+got_t:5.1f} t")
print(f"   void-filled total: {filled_t:.1f} t / {filled_v:.1f} m3 carried at ZERO extra freight")

rest_t = sum(remaining[t][0] for t in VOIDFILL_POOL) + m("G2")
rest_v = sum(remaining[t][1] for t in VOIDFILL_POOL) + v("G2")
rest_t -= m("G2"); rest_v -= v("G2")     # G2 = Pt gauze -> AIR freight, separate
n_box_B = math.ceil(max(rest_v/(BOX_VOL*STOW_LOOSE), rest_t/min(BOX_PAYLOAD, ROAD_CAP_T-3.9)))
print(f"   residual China cargo {rest_t:.1f} t / {rest_v:.1f} m3 -> {n_box_B} x 40'HC")

rowsB=[]; costB=0.0
c_mod = 6*(RATE["ocean_40HC"]+RATE["dest_port_40"]+RATE["inland_40_tarkwa"]) + 6*800
rowsB.append(f"{'6 ISO-frame modules (slot rate, no box hire)':<48}   6 x  $7,100 = ${c_mod:>9,.0f}")
costB += c_mod
c_bx = n_box_B*(RATE["ocean_40HC"]+RATE["dest_port_40"]+RATE["inland_40_tarkwa"])
rowsB.append(f"{'residual crated China cargo':<48} {n_box_B:>3} x  $6,300 = ${c_bx:>9,.0f}")
costB += c_bx
c_air = 3.2*95 + 2500      # 3.2 kg PGM air freight + secure handling/escort. SRC=ASSUMPTION
rowsB.append(f"{'Pt/Rh gauze by air (3.2 kg), secure':<48}              ${c_air:>9,.0f}")
costB += c_air
rowsB.append(f"{'OOG flat racks':<48}   0            ${0:>9,.0f}")
fobB = sum(c(t) for t in CN_ITEMS)*1000
insB = (fobB+costB)*RATE["insurance_pct"]; costB += insB
rowsB.append(f"{'marine insurance 1.1% of CIF':<48}              ${insB:>9,.0f}")
B_ex  = costB + levies(fobB+costB, True)
B_nex = costB + levies(fobB+costB, False)
rowsB.append(f"{'import levies (exempt / non-exempt)':<48}    ${levies(fobB+costB,True):>7,.0f} / ${levies(fobB+costB,False):,.0f}")
report("SCENARIO B totals", rowsB, B_ex)
print(f"   {'(non-exempt variant)':<48} ${B_nex:>10,.0f}")
print(f"   {'FOB value shipped from China':<48} ${fobB:>10,.0f}")

# --- local procurement cost of the Ghana scope ---
gh_cn_fob = sum(c(t) for t in GH_ITEMS if t not in NEVER_SHIPPED)*1000
gh_local  = sum(c(t)*1000*K[t] for t in GH_ITEMS if t not in NEVER_SHIPPED)
print(f"\n   Ghana-scope ex-works: China FOB ${gh_cn_fob:,.0f}  ->  Ghana delivered ${gh_local:,.0f} "
      f"(+${gh_local-gh_cn_fob:,.0f}, +{100*(gh_local/gh_cn_fob-1):.1f}%)")

# =====================================================================
# TOTAL DELIVERED-TO-SITE COMPARISON
# =====================================================================
civil = c("C1")*1000
A_tot_ex  = fobA + A_ex  + civil
A_tot_nex = fobA + A_nex + civil
B_tot_ex  = fobB + B_ex  + gh_local + civil
B_tot_nex = fobB + B_nex + gh_local + civil
print("\n=== DELIVERED-TO-SITE TOTAL (materials+equipment+freight+levies+civil, no erection) ===")
print(f"   A exempt      ${A_tot_ex:>12,.0f}      A non-exempt  ${A_tot_nex:>12,.0f}")
print(f"   B exempt      ${B_tot_ex:>12,.0f}      B non-exempt  ${B_tot_nex:>12,.0f}")
print(f"   SAVING        ${A_tot_ex-B_tot_ex:>12,.0f} ({100*(A_tot_ex-B_tot_ex)/A_tot_ex:5.2f}%)"
      f"   |  ${A_tot_nex-B_tot_nex:>10,.0f} ({100*(A_tot_nex-B_tot_nex)/A_tot_nex:5.2f}%)")
print(f"   logistics line only:  A ${A_ex:,.0f} -> B ${B_ex:,.0f}  "
      f"(-{100*(1-B_ex/A_ex):.1f}%) exempt;  A ${A_nex:,.0f} -> B ${B_nex:,.0f} (-{100*(1-B_nex/A_nex):.1f}%) non-exempt")
print(f"   shipping units: A = {noogA} OOG + {nboxA} box = {noogA+nboxA};  B = 6 ISO-frame + {n_box_B} box = {6+n_box_B}")
