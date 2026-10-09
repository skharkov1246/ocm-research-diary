import json
MOLV=23.645e-3  # m3/mol at 60F,1atm
mol_per_MMscf=28316.85/MOLV  # mol
comps={
 "pure_CH4":{"CH4":1.0},
 "wet_APG":{"CH4":0.80,"C2H6":0.10,"C3H8":0.05,"C4H10":0.03,"CO2":0.02},
}
Cn={"CH4":1,"C2H6":2,"C3H8":3,"C4H10":4,"CO2":1}
HHV={"CH4":1010,"C2H6":1770,"C3H8":2516,"C4H10":3262,"CO2":0} # Btu/scf
out={}
for k,c in comps.items():
    nC=sum(c[s]*Cn[s] for s in c)
    co2_t=mol_per_MMscf*nC*44.01e-6       # t CO2 per MMscf if fully burned (incl. inert CO2)
    ch4_t=mol_per_MMscf*c.get("CH4",0)*16.04e-6
    hhv=sum(c[s]*HHV[s] for s in c)
    yr=365
    r={"t_CO2_per_MMscf":round(co2_t,1),"t_CH4_per_MMscf":round(ch4_t,2),"Btu_scf":round(hhv),
       "CH4_t_y":round(ch4_t*yr),"CO2_full_t_y":round(co2_t*yr)}
    # baseline flare emissions (t CO2e/y) per 1 MMscf/d
    for gwp_name,gwp in [("AR5",28),("AR6",29.8)]:
        for dre in [0.98,0.95,0.911,0.0]:
            # combustion efficiency applied to all C; uncombusted carbon assumed released as CH4-equivalent for CH4 share only (conservative: non-CH4 slip ignored)
            co2=co2_t*yr*dre if dre>0 else (co2_t-ch4_t*44.01/16.04)*yr*0  # venting: no combustion CO2
            slip=ch4_t*yr*(1-dre)
            r[f"BE_{gwp_name}_DRE{dre}"]=round(co2+slip*gwp)
    out[k]=r
# creditable volumes per route (wet APG basis), per 1 MMscf/d, before availability
w=out["wet_APG"]; p=out["pure_CH4"]
routes={}
# R1 AM0009: BE = CO2 of combustion of recovered gas (displaces NG); PE = 3-8% (compression/processing/transport)
for nm,base in [("pure",p["CO2_full_t_y"]),("wet",w["CO2_full_t_y"])]:
    routes[f"AM0009_{nm}"]=[round(base*(1-0.08)),round(base*(1-0.03))]
# R2 fuel switch captive power: GWh from captive_out (low/base/high 19.8/25.1/30.4)
for g in [19.8,25.1,30.4]:
    routes[f"fuelswitch_conservative_{g}GWh"]=[round(g*1000*(0.70-0.60)),round(g*1000*(0.80-0.50))]
    routes[f"fuelswitch_consequential_{g}GWh"]=[round(g*1000*0.70),round(g*1000*0.80)]
# R3 mining: slip difference genset 99.5% vs flare baseline
for dre in [0.98,0.95,0.911]:
    routes[f"mining_slip_DRE{dre}"]=round(p["CH4_t_y"]*0.80/1.0*(dre and (1-dre-0.005))*28) # wet gas CH4 share 0.8
out["routes_t_y"]=routes
json.dump(out,open("/tmp/claude-0/-home-user-ocm-research-diary/50ead64b-47e2-5bdd-b34b-09f6667b6ca9/scratchpad/carbon/cc_out.json","w"),indent=1)
print(json.dumps(out,indent=1))
