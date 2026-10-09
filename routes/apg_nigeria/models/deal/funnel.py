import csv, math, sys
from collections import defaultdict
rows=list(csv.DictReader(open(sys.argv[1])))
hubs={'PH_TransAmadi':(4.815,7.035),'Onne':(4.71,7.15),'Eleme':(4.79,7.12),'Aba':(5.12,7.37),'Warri':(5.52,5.75),
      'Ughelli':(5.49,5.99),'Sapele':(5.89,5.68),'Owerri':(5.48,7.03),'Benin':(6.34,5.63),'Yenagoa':(4.93,6.27)}
def hav(a,b,c,e):
    R=6371;p=math.radians
    dl=p(c-a);dn=p(e-b)
    h=math.sin(dl/2)**2+math.cos(p(a))*math.cos(p(c))*math.sin(dn/2)**2
    return 2*R*math.asin(math.sqrt(h))
hist=defaultdict(dict)
for r in rows:
    try: hist[r['id']][int(r['Year'])]=hist[r['id']].get(int(r['Year']),0)+float(r['MMSCFD'])
    except: pass
d=[r for r in rows if r['Year']=='2025' and float(r['BCM'])>0]
for r in d:
    la,lo=float(r['lat']),float(r['lon'])
    ds={h:hav(a,b,la,lo) for h,(a,b) in hubs.items()}
    r['hub']=min(ds,key=ds.get); r['km']=ds[r['hub']]; r['q']=float(r['MMSCFD'])
print('all 2025:',len(d),round(sum(r['q'] for r in d),1))
on=[r for r in d if r['Location']=='ONSHORE']
print('onshore:',len(on),round(sum(r['q'] for r in on),1))
g5=[r for r in on if r['q']>=5]
print('onshore >=5:',len(g5),round(sum(r['q'] for r in g5),1))
for R in [15,20,25,30,40]:
    y=[r for r in g5 if r['km']<=R]; print(f' within {R} km: n={len(y)} vol={sum(r["q"] for r in y):.1f}')
g5.sort(key=lambda r:r['km'])
print('\nname | operator | 2025 | hub | km | min23-25 | mean19-25 | 2019..2025')
for r in g5:
    h=hist[r['id']]
    s=' '.join(f'{h.get(y,0):.1f}' for y in range(2019,2026))
    mn=min(h.get(y,0) for y in (2023,2024,2025)); mean=sum(h.get(y,0) for y in range(2019,2026))/7
    print(f"{r['Field name']} | {r['Operator']} | {r['q']:.1f} | {r['hub']} | {r['km']:.1f} | {mn:.1f} | {mean:.1f} | {s} | {r['lat'][:6]},{r['lon'][:6]}")
# also 3-5 within 15 km
print('\nonshore 3-5 within 20km:')
for r in sorted([r for r in on if 3<=r['q']<5 and r['km']<=20],key=lambda r:r['km']):
    h=hist[r['id']]; s=' '.join(f'{h.get(y,0):.1f}' for y in range(2019,2026))
    print(f"{r['Field name']} | {r['Operator']} | {r['q']:.1f} | {r['hub']} | {r['km']:.1f} | {s}")
