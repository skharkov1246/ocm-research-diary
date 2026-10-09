# Газ для станции под якорь: эмпирика VIIRS (World Bank, 2012-2025) для наземных факелов >=5 MMscf/д.
# Вопрос: какую долю времени станция, рассчитанная на долю f от дебита года 0, получает свой газ,
# и каково годовое падение дебита. Метка результатов: [О] по ../../wb/nigeria_flares.csv [П: World Bank VIIRS].
import csv, math, os, statistics as st, json
HERE=os.path.dirname(os.path.abspath(__file__))
CSV=os.path.join(HERE,'..','..','wb','nigeria_flares.csv')
rows=list(csv.DictReader(open(CSV)))
H={}; meta={}
for r in rows:
    try: q=float(r['MMSCFD']); y=int(r['Year'])
    except: continue
    H.setdefault(r['id'],{}); H[r['id']][y]=H[r['id']].get(y,0)+q
    meta[r['id']]=(r['Field name'],r['Operator'],r['Location'])
def util(f, y0s=range(2014,2020), horizon=6, qmin=5.0):
    """средняя по факелам доля min(q_t, f*q0)/(f*q0) за горизонт после года y0 (станция на f от дебита года 0)"""
    us=[]
    for i,h in H.items():
        if meta[i][2]!='ONSHORE': continue
        for y0 in y0s:
            q0=h.get(y0,0)
            if q0<qmin: continue
            cap=f*q0; ys=[h.get(y0+k,0) for k in range(1,horizon+1) if (y0+k)<=2025]
            if len(ys)<horizon: continue
            us.append(sum(min(q,cap) for q in ys)/(cap*len(ys)))
    return round(st.mean(us),3), round(st.median(us),3), len(us)
def decline(y0s=range(2012,2020), horizon=6, qmin=5.0):
    """медиана логарифмического наклона дебита за horizon лет (всех наземных >=qmin в году 0)"""
    sl=[]
    for i,h in H.items():
        if meta[i][2]!='ONSHORE': continue
        for y0 in y0s:
            q0=h.get(y0,0)
            if q0<qmin: continue
            ys=[h.get(y0+k,0) for k in range(0,horizon+1)]
            if min(ys)<=0.05: sl.append(None); continue   # факел погас/ушёл в трубу
            n=len(ys); xs=list(range(n)); lx=[math.log(v) for v in ys]
            mx=sum(xs)/n; ml=sum(lx)/n
            b=sum((x-mx)*(l-ml) for x,l in zip(xs,lx))/sum((x-mx)**2 for x in xs)
            sl.append(b)
    ok=[s for s in sl if s is not None]
    return dict(n=len(sl), gone_share=round(1-len(ok)/len(sl),3), median_annual=round(math.exp(st.median(ok))-1,3),
                mean_annual=round(math.exp(st.mean(ok))-1,3), p25=round(math.exp(sorted(ok)[len(ok)//4])-1,3))
out={}
for f in (1.0,0.75,0.52,0.39,0.17):
    out[f'util_f{f}']=util(f)
out['util_f1.0_2020_25']=util(1.0,y0s=[2019],horizon=6)
out['decline_6y']=decline()
out['decline_3y_2019_22']=decline(y0s=range(2019,2023),horizon=3)
for name in ('Ughelli East','Oredo','Sapele','Oben'):
    i=[k for k,v in meta.items() if v[0]==name][0]
    out[name]={y:round(H[i].get(y,0),1) for y in range(2012,2026)}
for k,v in out.items(): print(k,v)
json.dump(out,open(os.path.join(HERE,'viirs_supply_out.json'),'w'),ensure_ascii=False,indent=1)

# --- условно на фильтр воронки: факел стабилен (минимум трёх лет до y0 включительно >= 5 MMscf/д) ---
def util_stable(f, y0s=range(2014,2020), horizon=6, qmin=5.0, ref='min3'):
    us=[]; gone=0
    for i,h in H.items():
        if meta[i][2]!='ONSHORE': continue
        for y0 in y0s:
            hist=[h.get(y0-k,0) for k in range(3)]
            if min(hist)<qmin: continue
            q0=min(hist) if ref=='min3' else h.get(y0,0)
            cap=f*q0; ys=[h.get(y0+k,0) for k in range(1,horizon+1)]
            if any(q<0.05 for q in ys): gone+=1
            us.append(sum(min(q,cap) for q in ys)/(cap*len(ys)))
    return dict(mean=round(st.mean(us),3), median=round(st.median(us),3), p25=round(sorted(us)[len(us)//4],3), n=len(us), gone=gone)
out2={}
for f in (1.0,0.52,0.39,0.17):
    out2[f'stable_min3_f{f}']=util_stable(f)
    out2[f'stable_y0_f{f}']=util_stable(f,ref='y0')
for k,v in out2.items(): print(k,v)
out.update(out2)
json.dump(out,open(os.path.join(HERE,'viirs_supply_out.json'),'w'),ensure_ascii=False,indent=1)
