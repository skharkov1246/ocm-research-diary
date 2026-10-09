import sys
sys.path.insert(0,'/tmp/claude-0/-home-user-ocm-research-diary/50ead64b-47e2-5bdd-b34b-09f6667b6ca9/scratchpad/meoh')
def npv(r,c): return sum(x/(1+r)**i for i,x in enumerate(c))
def irr(c):
    lo,hi=-0.9,1.0
    if npv(lo,c)*npv(hi,c)>0: return None
    for _ in range(200):
        m=(lo+hi)/2
        if npv(lo,c)*npv(m,c)<=0: hi=m
        else: lo=m
    return m
def proj(capex,e,wc=0,tax=0.30,life=12,ramp=(0.6,0.9),dec=0.04):
    cf=[-capex/2,-capex/2-wc]
    for y in range(1,life+1):
        f=ramp[y-1] if y<=2 else 1.0
        if y>=4: f*=(1-dec)**(y-3)
        ee=e*f; dep=capex/10 if y<=10 else 0
        cf.append(ee-max(0,(ee-dep)*tax))
    cf[-1]+=wc
    return cf
# HCHO capacity check
uf=17737; form=8053
hcho=uf*0.231+form*0.37
print('HCHO t/y',round(hcho),'t/d',round(hcho/330,1),'formalin37-eq t/d',round(hcho/0.37/330,1))
# powder->liquid equivalent price
for p in (552,):
    print('liquid65 eq of powder95 CIF',round(p*0.65/0.95),'; of powder at 100% solids',round(p*0.65))
# chain B base re-run with price & demand variants
def chainB(uf,form,p_uf,p_form=420,p_me=375,t=9900,capex=32.8e6,urea_p=480,wc_frac=0.0,tax=0.30,meoh_cost_buy=None):
    need=uf*0.27+form*0.45
    sold=max(0,t-need) if meoh_cost_buy is None else 0
    rev=uf*p_uf+form*p_form+sold*p_me
    urea=uf*0.42*(urea_p+20)
    gas=1000*1.05*330*0.6 if meoh_cost_buy is None else 0
    staff=1.25e6 if meoh_cost_buy is None else 0.45e6
    maint=capex*0.045
    util=(uf+form)*12; fr=(uf+form)*50 if meoh_cost_buy is None else (uf+form)*20
    ins=capex*0.02+(0.2e6 if meoh_cost_buy is None else 0) if meoh_cost_buy is None else capex*0.015
    cat=t*15 if meoh_cost_buy is None else 0
    buy=need*meoh_cost_buy if meoh_cost_buy else 0
    opex=urea+gas+staff+maint+util+fr+ins+cat+0.4e6+buy
    e=rev-opex
    wc=wc_frac*rev
    cf=proj(capex,e,wc=wc,tax=tax)
    return dict(rev=round(rev/1e6,2),opex=round(opex/1e6,2),ebitda=round(e/1e6,2),payback=round(capex/e,1) if e>0 else None,npv15=round(npv(0.15,cf)/1e6,1),npv20=round(npv(0.20,cf)/1e6,1),irr=(round(irr(cf)*100,1) if irr(cf) is not None else None))
print('B base as-is           ',chainB(17737,8053,580))
print('B base, UF liquid $380 ',chainB(17737,8053,380))
print('B base, UF $380, WC 20%, tax 34%',chainB(17737,8053,380,wc_frac=0.2,tax=0.34))
print('B board-demand UF 8kt $580',chainB(8000,8053,580))
print('B board-demand UF 8kt $380',chainB(8000,8053,380))
print('B UF 8kt $450, form $350',chainB(8000,8053,450,p_form=350))
print('B base $380 urea $380 (export parity)',chainB(17737,8053,380,urea_p=380))
print('C control as-is ',chainB(17737,8053,580,capex=8e6,meoh_cost_buy=480))
print('C control $380  ',chainB(17737,8053,380,capex=8e6,meoh_cost_buy=480))
print('C control $450  ',chainB(17737,8053,450,capex=8e6,meoh_cost_buy=480))
print('C control 8kt $450 ',chainB(8000,8053,450,capex=8e6,meoh_cost_buy=480))
# breakeven UF price for chain B base NPV15=0
for p in range(380,800,10):
    r=chainB(17737,8053,p)
    if r['npv15']>=0: print('B NPV15=0 at UF price ~',p); break
for p in range(300,800,5):
    r=chainB(17737,8053,p,capex=8e6,meoh_cost_buy=480)
    if r['npv15']>=0: print('C NPV15=0 at UF price ~',p); break
# raw material cost floor of liquid UF in Nigeria
print('raw mat/t UF: urea',0.42*500,'+ MeOH',round(0.27*480),'=',round(0.42*500+0.27*480))
for p in range(580,1500,10):
    if chainB(17737,8053,p)['npv15']>=0: print('B NPV15=0 at UF',p); break
