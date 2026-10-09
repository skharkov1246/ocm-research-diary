# Скептик BOOT: что на самом деле означает «загрузка 0.78» и не учтено ли падение дебита дважды.
# Эмпирика по VIIRS (wb/nigeria_flares.csv): площадка размером по году y, реализованная загрузка в годы y+1..y+k
# = mean(min(Q_{y+k}, cap)/Q_y). Сравниваем с моделью boot.py: lf*(1-decl)^(t-1) (без готовности).
import csv, collections, os, statistics as st
HERE = os.path.dirname(os.path.abspath(__file__))
P = os.path.join(HERE, '..', '..', 'wb', 'nigeria_flares.csv')
s = collections.defaultdict(dict)
for r in csv.DictReader(open(P)):
    s[r['id']][int(r['Year'])] = float(r['MMSCFD'])
def emp(lo, hi, cap=1.0, H=(1, 2, 3, 4, 5, 6)):
    out = {}
    for h in H:
        v = []
        for k, d in s.items():
            for y, q in d.items():
                if lo <= q <= hi and (y+h) <= 2025:
                    v.append(min(d.get(y+h, 0.0), cap*q)/q)
        out[h] = (round(st.mean(v), 3), len(v))
    return out
for lo, hi in ((0.5, 2), (3, 8), (8, 30), (3, 30)):
    print(f'Q {lo}-{hi} MMscf/d, cap 1.0xQ :', emp(lo, hi))
    print(f'Q {lo}-{hi} MMscf/d, cap 1.5xQ :', emp(lo, hi, 1.5))
print('model boot.py base lf*(1-0.08)^(h-1):', {h: round(0.78*0.92**(h-1), 3) for h in range(1, 7)})
print('model mean 10y lf*decl:', round(sum(0.78*0.92**(t-1) for t in range(1, 11))/10, 3),
      ' x up 0.92:', round(0.92*sum(0.78*0.92**(t-1) for t in range(1, 11))/10, 3))
# Именованные площадки: ряд 2019-2025
names = ['Ughelli East', 'Oredo', 'Odidi', 'Obiafu', 'Kwale', 'Oben', 'Afam', 'Opuama', 'Saghara']
for r in csv.DictReader(open(P)):
    pass
rows = list(csv.DictReader(open(P)))
for n in names:
    ids = sorted({r['id'] for r in rows if n.lower() in r['Field name'].lower()})
    for i in ids:
        d = s[i]; print(n, i, {y: round(d[y], 1) for y in sorted(d) if y >= 2019})
