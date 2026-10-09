import json
cap=json.load(open("/tmp/claude-0/-home-user-ocm-research-diary/50ead64b-47e2-5bdd-b34b-09f6667b6ca9/scratchpad/power/captive_out.json"))
def npv(cfs,r): return sum(cf/(1+r)**t for t,cf in enumerate(cfs,start=0))
# price scenarios gross $/t and deduction stack
SC={
 "low_VCM":   dict(P=4.0,  ca=0.0, ndc=0.0, broker=0.12, share=0.25, reg=0.20),
 "mid_CORSIA":dict(P=10.0, ca=0.08,ndc=0.05,broker=0.10, share=0.25, reg=0.20),
 "high_ITMO": dict(P=20.0, ca=0.08,ndc=0.05,broker=0.08, share=0.10, reg=0.20),
}
def netp(s): return s["P"]*(1-s["ca"]-s["ndc"]-s["broker"]-s["share"])-s["reg"]
# per-site fixed costs (pre-audit screen) and audit multipliers
DEV_capex=(0.15,0.30)   # $M: PDD, validation, NCCC fees $4k, monitoring meters/GC
MRV_opex=(0.04,0.08)    # $M/y verification, MRV staff, registry account
AUD_CAPEX=(2.4,2.8); AUD_OPEX=(1.7,2.4); DEMAND=(1.3,2.7); P_SUCCESS=0.6
# creditable volume kt/y per 1 MMscf/d, route A (AM0009: CNG/LPG/pipeline), lower bound applies availability
VOL_A=(17.7*0.85, 24.07*0.95)
res={"netp":{k:round(netp(v),2) for k,v in SC.items()}}
def site(q,vol_kt,scn,audit,years_cred,delay=2,r=0.15,scale_fixed=True):
    s=SC[scn]; n=netp(s)
    rev=vol_kt*1000*q*n/1e6  # $M/y
    dcap=DEV_capex[0 if not audit else 1]*(1 if q<5 else 1.5)
    mrv=MRV_opex[0 if not audit else 1]*(1 if q<5 else 1.3)
    if audit:
        dcap*=AUD_CAPEX[0]; mrv*=AUD_OPEX[0]; rev/=DEMAND[0]
    ebitda=rev-mrv
    cfs=[-dcap]+[0]*(delay-1)+[ebitda]*years_cred
    return dict(rev=round(rev,3),mrv=round(mrv,3),ebitda=round(ebitda,3),dev_capex=round(dcap,3),
                npv15=round(npv(cfs,0.15),3),npv20=round(npv(cfs,0.20),3),
                ev_npv15=round(npv(cfs,0.15)*(P_SUCCESS if audit else 1),3))
out={}
for q in (1,5,15):
    for scn in SC:
        for vol_name,vol in (("lo",VOL_A[0]),("hi",VOL_A[1])):
            for yrs in (3,8):
                out[f"A_q{q}_{scn}_{vol_name}_{yrs}y_pre"]=site(q,vol,scn,False,yrs)
                out[f"A_q{q}_{scn}_{vol_name}_{yrs}y_aud"]=site(q,vol,scn,True,yrs)
# corpus kicker reproduction
out["corpus_kicker"]={"vol_kt":20,"P":[30,80],"M$":[0.6,1.6]}
# fuel-switch captive power (route B), conservative t/MWh 0.1-0.3, base GWh
for q,key in ((1,"1MMscfd_base"),(5,"5MMscfd_base"),(15,"15MMscfd_base")):
    c=cap[key]; gwh=c["y1_GWh"]
    for scn in SC:
        n=netp(SC[scn])
        for ef in (0.1,0.3):
            vol=gwh*1000*ef  # t
            rev=vol*n/1e6
            mrv=MRV_opex[1]*AUD_OPEX[0]*(1 if q<5 else 1.3)
            add=rev/DEMAND[0]-mrv
            pb0=c["capex"]/c["y1_ebitda"]
            pb1=c["capex"]/(c["y1_ebitda"]+max(add,-1e9))
            out[f"B_q{q}_{scn}_ef{ef}"]=dict(GWh=gwh,kt=round(vol/1000,1),add_ebitda=round(add,3),
                 base_ebitda=c["y1_ebitda"],capex=c["capex"],payback0=round(pb0,2),payback1=round(pb1,2),
                 share_of_ebitda_pct=round(100*add/c["y1_ebitda"],1))
# penalty shadow: $3.50/Mscf
out["penalty_per_MMscfd_M$"]=round(3.5*365000/1e6,3)
out["penalty_per_tCO2"]=[round(3.5*1000/68.0,1),round(3.5*1000/52.7,1)]
out["nuprc_collect_rate"]=dict(collected_M=round(521.87e9/1330/1e6,0),theor_WB=round(639*365000*3.5/1e6,0),theor_NUPRC=round(203.97e6*3.5/1e6,0))
json.dump(out,open("/tmp/claude-0/-home-user-ocm-research-diary/50ead64b-47e2-5bdd-b34b-09f6667b6ca9/scratchpad/carbon/econ_out.json","w"),indent=1)
print(json.dumps(res))
for k,v in out.items():
    if k.startswith("A_q") and ("_8y_" in k or "_3y_" in k): print(k,v)
for k,v in out.items():
    if not k.startswith("A_q"): print(k,v)
