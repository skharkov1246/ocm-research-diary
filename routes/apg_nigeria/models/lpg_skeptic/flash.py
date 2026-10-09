# Peng-Robinson isothermal flash, kij=0 (crude), composition (a) NJTD Table 1
import math
R=8.314
comp={'N2':(0.10,126.2,33.9e5,0.037),'CO2':(0.46,304.1,73.8e5,0.225),'C1':(82.72,190.6,46.0e5,0.011),
'C2':(6.73,305.3,48.7e5,0.099),'C3':(5.68,369.8,42.5e5,0.152),'iC4':(1.56,408.1,36.5e5,0.181),
'nC4':(1.60,425.1,38.0e5,0.200),'iC5':(0.40,460.4,33.8e5,0.227),'nC5':(0.26,469.7,33.7e5,0.251),
'C6':(0.16,507.6,30.2e5,0.301),'C7':(0.20,540.2,27.4e5,0.350)}
names=list(comp); z=[comp[n][0] for n in names]; s=sum(z); z=[x/s for x in z]
def cubic(a2,a1,a0):
    p=a1-a2*a2/3; q=2*a2**3/27-a2*a1/3+a0; D=(q/2)**2+(p/3)**3
    if D>0:
        u=(-q/2+math.sqrt(D)); v=(-q/2-math.sqrt(D))
        cr=lambda t: math.copysign(abs(t)**(1/3),t)
        return [cr(u)+cr(v)-a2/3]
    r=math.sqrt(-(p/3)**3); ph=math.acos(max(-1,min(1,-q/(2*r))))
    m=2*math.sqrt(-p/3)
    return [m*math.cos((ph+2*math.pi*k)/3)-a2/3 for k in range(3)]
def params(T):
    a=[];b=[]
    for n in names:
        _,Tc,Pc,w=comp[n]; m=0.37464+1.54226*w-0.26992*w*w
        al=(1+m*(1-math.sqrt(T/Tc)))**2
        a.append(0.45724*R*R*Tc*Tc/Pc*al); b.append(0.07780*R*Tc/Pc)
    return a,b
def phis(x,T,P,phase):
    a,b=params(T); n=len(x)
    am=sum(x[i]*x[j]*math.sqrt(a[i]*a[j]) for i in range(n) for j in range(n)); bm=sum(x[i]*b[i] for i in range(n))
    A=am*P/(R*T)**2; B=bm*P/(R*T)
    # cubic Z^3-(1-B)Z^2+(A-3B^2-2B)Z-(AB-B^2-B^3)=0
    c=[1,-(1-B),A-3*B*B-2*B,-(A*B-B*B-B**3)]
    r=[x for x in cubic(c[1],c[2],c[3]) if x>B]
    Z=max(r) if phase=='v' else min(r)
    out=[]
    for i in range(n):
        sa=sum(x[j]*math.sqrt(a[i]*a[j]) for j in range(n))
        lnphi=b[i]/bm*(Z-1)-math.log(Z-B)-A/(2*math.sqrt(2)*B)*(2*sa/am-b[i]/bm)*math.log((Z+(1+math.sqrt(2))*B)/(Z+(1-math.sqrt(2))*B))
        out.append(math.exp(lnphi))
    return out
def flash(T,P):
    K=[comp[n][2]/P*math.exp(5.373*(1+comp[n][3])*(1-comp[n][1]/T)) for n in names]
    for it in range(300):
        f=lambda V: sum(z[i]*(K[i]-1)/(1+V*(K[i]-1)) for i in range(len(z)))
        lo,hi=1e-9,1-1e-9
        if f(lo)<0: V=0.0
        elif f(hi)>0: V=1.0
        else:
            for _ in range(100):
                m=(lo+hi)/2
                if f(m)>0: lo=m
                else: hi=m
            V=(lo+hi)/2
        x=[z[i]/(1+V*(K[i]-1)) for i in range(len(z))]; y=[K[i]*x[i] for i in range(len(z))]
        sx=sum(x); sy=sum(y); x=[v/sx for v in x]; y=[v/sy for v in y]
        pl=phis(x,T,P,'l'); pv=phis(y,T,P,'v')
        Kn=[pl[i]/pv[i] for i in range(len(z))]
        if max(abs(Kn[i]/K[i]-1) for i in range(len(z)))<1e-8: K=Kn;break
        K=Kn
    L=1-V
    rec={n:L*x[i]/z[i] for i,n in enumerate(names)}
    return V,rec
for T_C,Pbar in [(-35,27.5),(-35,40),(-40,27.5),(-30,27.5),(-35,20),(-90,27.5)]:
    V,rec=flash(T_C+273.15,Pbar*1e5)
    print(T_C,'C',Pbar,'bar  Lfrac',round(1-V,4),' rec C2 %.2f C3 %.2f iC4 %.2f nC4 %.2f C5 %.2f'%(rec['C2'],rec['C3'],rec['iC4'],rec['nC4'],rec['iC5']))
