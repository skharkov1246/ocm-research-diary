import sys, json
sys.path.insert(0,'/tmp/claude-0/-home-user-ocm-research-diary/50ead64b-47e2-5bdd-b34b-09f6667b6ca9/scratchpad/minilng')
import importlib.util
spec=importlib.util.spec_from_file_location('m','/tmp/claude-0/-home-user-ocm-research-diary/50ead64b-47e2-5bdd-b34b-09f6667b6ca9/scratchpad/minilng/minilng_model.py')
# re-implement pre-tax levelized cost to keep both chains on same basis
ESC=800/576.1; HHV_LNG=52.6
def crf(r,n): return r*(1+r)**n/((1+r)**n-1)
def lng_plant(Q, C, LF=0.85, gas=0.5, kwh=1.1, r=0.15, n=15):
    feed=1.05*1000*Q
    fuel_ratio=kwh/0.38*3.6/(HHV_LNG*1.055)
    lng=feed*0.95*0.97/(1+fuel_ratio)*365*LF          # MMBtu/yr
    ngl=1.0*Q*365*LF*450/1e6
    fixed=(0.5+0.09*Q)+(0.4+0.04*Q)+(0.2+0.02*Q)+C*0.06
    var=gas*feed*365*LF/1e6+0.45*lng/1e6
    return (C*crf(r,n)+fixed+var-ngl)*1e6/lng
def cng_mother(Q, C, LF=0.85, gas=0.5, r=0.15, n=15):
    feed=1.05*1000*Q; out=feed*0.94*365*LF            # 6% на привод компрессоров [Д]
    fixed=(0.4+0.05*Q)+(0.3+0.03*Q)+0.1+C*0.06
    var=gas*feed*365*LF/1e6+0.15*out/1e6               # расходники/ремонт компрессоров $0.15/MMBtu [Д]
    return (C*crf(r,n)+fixed+var)*1e6/out
res={}
for Q in [1,5,15]:
    for lab,C in [('low',None),('high',None)]:
        pass
for Q in [1,5,15]:
    C_lng_lo=(1.5607*Q+3.1243)*ESC*2.4
    C_lng_hi=((3.1214*Q+6.2485)*ESC+19*(Q/15)**0.6*ESC+3*15.5*Q/0.43*2000/1e6+1.0+0.1*Q+0.2*Q)*1.2
    C_cng_lo=(1.2+0.8*Q)*2.4/1.0 if False else (1.0+0.6*Q)*2.4   # $1.2 млн/станция (Nairametrics 2023) ~ 1-2 MMscf/d [В]; x2.4-2.8 аудит
    C_cng_hi=(1.0+0.6*Q)*2.8+0.5
    L_lo=lng_plant(Q,C_lng_lo); L_hi=lng_plant(Q,C_lng_hi)
    K_lo=cng_mother(Q,C_cng_lo); K_hi=cng_mother(Q,C_cng_hi)
    res[Q]=dict(capex_lng=(round(C_lng_lo,1),round(C_lng_hi,1)),capex_cng=(round(C_cng_lo,1),round(C_cng_hi,1)),
                lng_plant_usd_mmbtu=(round(L_lo,2),round(L_hi,2)),cng_station_usd_mmbtu=(round(K_lo,2),round(K_hi,2)))
    # distance breakeven: LNG: L + 2Dc/894 + R_l ; CNG: K + 2Dc/P_c + R_c
    for c in [1.5,2.5]:
        for Pc in [496,677]:
            for (L,K,tag) in [(L_lo,K_lo,'lowcapex'),(L_hi,K_hi,'highcapex')]:
                Rl,Rc=1.7,0.75
                prem=(L+Rl)-(K+Rc)
                slope=2*c*(1/Pc-1/894)
                D=prem/slope
                res[Q][f'D*_c{c}_P{Pc}_{tag}']=round(D)
print(json.dumps(res,indent=1,ensure_ascii=False))
# delivered cost at 500 km for Q=5 mid
for Q in [1,5,15]:
    r=res[Q]; Lm=sum(r['lng_plant_usd_mmbtu'])/2; Km=sum(r['cng_station_usd_mmbtu'])/2
    for D in [100,250,500,800]:
        tl=2*D*2.0/894; tc=2*D*2.0/586
        print(Q,D,'LNG delivered',round(Lm+tl+1.7,2),'CNG delivered',round(Km+tc+0.75,2))
json.dump(res,open('/tmp/claude-0/-home-user-ocm-research-diary/50ead64b-47e2-5bdd-b34b-09f6667b6ca9/scratchpad/minilng/cng_vs_lng_out.json','w'),indent=1)
