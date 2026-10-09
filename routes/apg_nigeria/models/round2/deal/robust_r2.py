# Раунд 2: устойчивость рекомендованного порога ворот 0 (₦495, DoP-газ+эскроу, EDTI не подтверждён,
# плановая стройка 2 г, газ на стопе $1.5) и судейских ₦316/350/377 на исправленной модели.
import os, sys, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from deal_r2 import *
S = STRUCTS; TL = TL_NEW; OPF = 0.15/1.05; OUT = {}
def kw_e(e): return dict(edti=e, depyrs=5 if e else 10)
B = {'1y': (dict(TL, delay=1), {}), '2y': (dict(TL, delay=2), dict(capex_split=[0.5, 0.5])), '2y_slip': (dict(TL, delay=2), {})}
for T in (316, 350, 377, 411, 440, 495, 530):
    for sn in ('25km', 'oredo30'):
        for b, (tl, ex) in B.items():
            for e in (True, False):
                for g in (0.60+OPF, 1.5):
                    x = run(dict(site(sn), p_dir=T, gas=g-OPF), S['dopgas'], tl, debt=(0.6, 0.14, 7), **ex, **kw_e(e))
                    OUT[f'T{T}|{sn}|{b}|{"E" if e else "N"}|gas{g:.2f}'] = dict(npv15=x['npv15'], npv20=x['npv20'], irr=x['irr'], min_dscr=x['min_dscr'])
for k, v in OUT.items():
    if k.startswith(('T495|', 'T350|25km', 'T316|25km', 'T377|25km', 'T440|')) and 'gas0.74' in k or k.startswith('T495|'):
        print(k, v)
json.dump(OUT, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'robust_r2_out.json'), 'w'), ensure_ascii=False, indent=1)

# потолок капитала (FEED/EPC) при подписанном тарифе: множитель cm и $ млн, при котором NPV20=0
caps = {}
for T in (440, 460, 495):
    for b, (tl, ex) in B.items():
        for e in (True, False):
            for g in (0.60+OPF, 1.5):
                lo, hi = 1.0, 6.0
                for _ in range(40):
                    m = (lo+hi)/2
                    if run(dict(site('25km'), p_dir=T, gas=g-OPF, cm=m), S['dopgas'], tl, **ex, **kw_e(e))['npv20'] >= 0: lo = m
                    else: hi = m
                caps[f'T{T}|25km|{b}|{"E" if e else "N"}|gas{g:.2f}'] = dict(cm=round(lo, 2), capex=round(station(dict(site('25km'), cm=lo))['C'], 1))
for k, v in caps.items():
    if '2y_slip' not in k: print('cap', k, v)
OUT['capex_cap'] = caps
json.dump(OUT, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'robust_r2_out.json'), 'w'), ensure_ascii=False, indent=1)
