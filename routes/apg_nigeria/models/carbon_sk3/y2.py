def npv(c,r): return sum(x/(1+r)**t for t,x in enumerate(c))
def netp(P,ca,ndc,br,sh,reg=0.2): return P*(1-ca-ndc-br-sh)-reg
NP={'VCM4':netp(4,0,0,.12,.25),'VCM6.8':netp(6.8,0,0,.10,.25),'C10':netp(10,.08,.05,.10,.25),'I20':netp(20,.08,.05,.08,.10)}
print({k:round(v,2) for k,v in NP.items()})
V0=17.7
FR={'CNG_low':0.40*0.85,'CNG_base':0.65*0.88,'pipe':0.85}
SCL={1:1.0,5:0.85,15:0.6}
def run(Q,route,price,dem,mrvm,capm,decl,yrs,delay,infl,r=0.15):
    f=FR[route]*(SCL[Q] if route!='pipe' else 1)
    v1=V0*Q*f  # kt/y year1
    cap=0.30*capm*(1.5 if Q>=5 else 1)
    mrv0=0.08*mrvm*(1.3 if Q>=5 else 1)
    cfs=[-cap]; cum=-cap; rows=[]
    # crediting vintages 1..yrs ; vintage t issued at t+delay-1 ; MRV paid each crediting year
    T=yrs+delay
    cash=[0.0]*(T+1); cash[0]=-cap
    for t in range(1,yrs+1):
        cash[t]-=mrv0*(1+infl)**t
        vol=v1*(1-decl)**(t-1)
        cash[t+delay-1]+=vol*NP[price]/dem/1000
    return round(v1,1),round(cap,2),round(mrv0,3),[round(x,3) for x in cash],round(npv(cash,r),2),round(npv(cash,0.20),2)
for Q in (1,5,15):
  for route in ('CNG_low','CNG_base','pipe'):
    for price in ('VCM4','C10','I20'):
      b=run(Q,route,price,1.3,1.7,2.4,0.06,4,2,0.03)
      w=run(Q,route,price,2.7,2.4,2.8,0.12,4,3,0.06)
      b8=run(Q,route,price,1.3,1.7,2.4,0.06,8,2,0.03)
      print(Q,route,price,'v1',b[0],'| best4y npv15',b[4],'cum y2',round(sum(b[3][:3]),2),'| best8y npv15',b8[4],'| worst4y npv15',w[4])
# EV correction example
for q,n,c in ((15,1.617,1.08),(15,3.378,1.08),(5,1.301,1.08),(1,-1.146,0.72)):
    sunk=0.5*c  # validation+PDD spent before LoA, half of dev capex
    print('EV',q,n,round(0.6*n,2),round(0.6*n-0.4*sunk,2))
