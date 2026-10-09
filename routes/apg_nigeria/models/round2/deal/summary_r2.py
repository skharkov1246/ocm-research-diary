# Раунд 2: сводные таблицы для решения совета (порог LOI, газ, замер, якорь, часы сети) по трём графикам стройки:
#   1y      -- стройка 1 год, капитал в t0 (как матрица судьи);
#   2y      -- плановая стройка 2 года, капитал 50/50 по годам (типичная задержка 12-36 мес. в Нигерии [В]);
#   2y_slip -- стройка 1 год + срыв на 1 год, весь капитал в t0 (как «delay2» судьи).
# Запуск после matrix_r2.py не обязателен: всё считается заново. Вывод: summary_r2_out.json + печать.
import os, sys, json, math
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from deal_r2 import *

HERE = os.path.dirname(os.path.abspath(__file__))
S = STRUCTS; TL = TL_NEW; OUT = {}
def kw_e(e): return dict(edti=e, depyrs=5 if e else 10)
def ceil5(x): return int(math.ceil(x/5.0)*5)
def USD(n): return round(n/FX, 3)
OPF = 0.15/1.05
def g_allin(a): return a-OPF
STOP = 1.5
BUILDS = {'1y': (dict(TL, delay=1), {}), '2y': (dict(TL, delay=2), dict(capex_split=[0.5, 0.5])),
          '2y_slip': (dict(TL, delay=2), {})}
STR = ('dopgas', 'esc', 'none', 'dopmargin')
SITES = ('25km', 'oredo30')

# ---- 1. матрица BE15/BE20 при базовом газе и порог LOI при газе на стопе
MAT = {}
for sn in SITES:
    for mw in (11.0, 8.5):
        for b, (tl, ex) in BUILDS.items():
            for s in STR:
                for e in (True, False):
                    p = dict(site(sn), anchor=mw)
                    a = be(p, S[s], tl, 'npv15', **ex, **kw_e(e)); c = be(p, S[s], tl, 'npv20', **ex, **kw_e(e))
                    t = ceil5(be(dict(p, gas=g_allin(STOP)), S[s], tl, 'npv20', **ex, **kw_e(e)))
                    tA = ceil5(be(dict(p, gas=g_allin(STOP)), S[s], tl, 'npv20', earnout=(0.30, 2), **ex, **kw_e(e)))
                    MAT[f'{sn}|{mw}|{b}|{s}|{"E" if e else "N"}'] = dict(be15=round(a), be20=round(c), T_LOI=t, T_LOI_A=tA)
OUT['matrix'] = MAT
print('BE15/BE20 при газе $0.74 всё включено | T_LOI путь Б / путь А (газ на стопе $1.5, округл. вверх до ₦5)')
for sn in SITES:
    for mw in (11.0, 8.5):
        print(f'== {sn}, {mw} МВт   [EDTI: 1г | 2г(50/50) | 2г срыв]   [без EDTI: 1г | 2г(50/50) | 2г срыв]')
        for s in STR:
            cells = []
            for e in ('E', 'N'):
                for b in BUILDS:
                    m = MAT[f'{sn}|{mw}|{b}|{s}|{e}']
                    cells.append(f"{m['be15']}/{m['be20']}→{m['T_LOI']}")
            print(f'  {SLABEL[s][:44]:44s}', ' | '.join(cells))

# ---- 2. газ: точка NPV20=0 по цене газа при подписанном тарифе T_LOI+dT
gp = {}
for b, (tl, ex) in BUILDS.items():
    for e in (True, False):
        T0 = MAT[f'25km|11.0|{b}|dopgas|{"E" if e else "N"}']['T_LOI']
        for dT in (0, 25, 50):
            lo, hi = 0.0, 15.0
            for _ in range(40):
                m = (lo+hi)/2
                if run(dict(NEW, p_dir=T0+dT, gas=g_allin(m)), S['dopgas'], tl, **ex, **kw_e(e))['npv20'] >= 0: lo = m
                else: hi = m
            gp[f'{b}|{"E" if e else "N"}|T={T0+dT}'] = round(lo, 2)
OUT['gas_point'] = gp
print('\nцена газа «всё включено» $/MMBtu, при которой NPV20=0 (25 км, 11 МВт, DoP-газ):', gp)

# ---- 3. порог замера газа (год 1, MMscf/д) при T_LOI, газ 0.74 и на стопе, падение 8% и 12%
meas = {}
for sn in ('25km',):
    for mw in (11.0, 8.5):
        for b, (tl, ex) in BUILDS.items():
            if b == '2y_slip': continue
            for e in (True, False):
                T0 = MAT[f'{sn}|{mw}|{b}|dopgas|{"E" if e else "N"}']['T_LOI']
                for allin in (0.60+OPF, STOP):
                    for dcl in (0.08, 0.12):
                        lo, hi = 0.5, 12.0
                        for _ in range(40):
                            m = (lo+hi)/2
                            p = dict(NEW, anchor=mw, Q=m, avail=1.0, decl=dcl, gas=g_allin(allin), p_dir=T0)
                            if run(p, S['dopgas'], tl, **ex, **kw_e(e))['npv20'] >= 0: hi = m
                            else: lo = m
                        meas[f'{mw}|{b}|{"E" if e else "N"}|T={T0}|gas{allin:.2f}|decl{dcl}'] = round(hi, 2)
OUT['gas_measure_min'] = meas
print('\nминимальный доступный газ года 1, MMscf/д (DCQ 11 МВт = %.2f; 8.5 МВт = %.2f):' % (
    station(NEW)['g_pow'], station(dict(NEW, anchor=8.5))['g_pow']))
for k, v in meas.items(): print('  ', k, v)

# ---- 4. порог якоря
LEGAL = 6.0/0.72
MWS = (5.0, 6.0, 7.0, 8.0, 8.5, 9.0, 10.0, 11.0, 12.0, 13.0, 15.0)
curve = {}
for sn in SITES:
    for b, (tl, ex) in BUILDS.items():
        for s in STR:
            for e in (True, False):
                curve[f'{sn}|{b}|{s}|{"E" if e else "N"}'] = [
                    (mw, round(be(dict(site(sn), anchor=mw, gas=g_allin(STOP)), S[s], tl, 'npv20', **ex, **kw_e(e)))) for mw in MWS]
OUT['anchor_curve_T20_gasstop'] = curve
def mw_for(xs, ceiling):
    if xs[-1][1] > ceiling: return None
    if xs[0][1] <= ceiling: return xs[0][0]
    for (m1, b1), (m2, b2) in zip(xs, xs[1:]):
        if b1 > ceiling >= b2: return round(m1+(b1-ceiling)*(m2-m1)/(b1-b2), 1)
ALTS = {'дизель 530': 530, 'дизель 650': 650, 'дизель 770': 770,
        'сеть 8 ч + дизель 650': (8*209.5+16*650)/24, 'сеть 12 ч + дизель 650': (12*209.5+12*650)/24}
thr = {}
for k, xs in curve.items():
    r = {}
    for nm, alt in ALTS.items():
        m = mw_for(xs, 0.85*alt)
        r[nm] = None if m is None else round(max(m, LEGAL), 1)
    thr[k] = r
OUT['anchor_threshold'] = thr; OUT['legal_floor'] = round(LEGAL, 2)
print('\nпорог якоря, МВт нетто = max(%.1f юридический, экономический при цене <= 0.85 x альтернатива покупателя); None = не проходит и при 15 МВт' % LEGAL)
for k, v in thr.items():
    if '|dopmargin|' in k: continue
    print('  ', k, v)
for k in ('25km|1y|dopgas|E', '25km|2y|dopgas|N', 'oredo30|2y|dopgas|N', 'oredo30|1y|dopgas|E'):
    print('   кривая', k, curve[k])

# ---- 5. часы сети
def hmax(T, Pd, Pg=209.5, s=0.0):
    return round(max(0.0, min(24.0, (Pd-T/(1-s))*24/(Pd-Pg))), 1)
grid = {}
for sn in SITES:
    for mw in (11.0, 8.5):
        for b in BUILDS:
            for s in ('dopgas', 'esc', 'none', 'dopmargin'):
                for e in ('E', 'N'):
                    T = MAT[f'{sn}|{mw}|{b}|{s}|{e}']['T_LOI']
                    grid[f'{sn}|{mw}|{b}|{s}|{e}'] = dict(T=T, **{f'Pd{Pd}_s{int(sv*100)}': hmax(T, Pd, s=sv) for Pd in (530, 650, 770) for sv in (0.0, 0.15)},
                                                          **{f'Pd{Pd}_s15_Pg225': hmax(T, Pd, 225.0, 0.15) for Pd in (650,)})
OUT['grid_hours'] = grid
print('\nмакс. часов сети в сутки (Band A ₦209.5; s = требуемая экономия покупателя)')
for k, v in grid.items():
    if '|dopgas|' in k or ('|none|' in k and '|11.0|' in k): print('  ', k, v)

# ---- 6. ценность DoP-маржи и EDTI в NPV15 при T_LOI базовой структуры
val = {}
for b, (tl, ex) in BUILDS.items():
    for e in (True, False):
        T = MAT[f'25km|11.0|{b}|dopgas|{"E" if e else "N"}']['T_LOI']
        x_g = run(dict(NEW, p_dir=T), S['dopgas'], tl, **ex, **kw_e(e))
        x_m = run(dict(NEW, p_dir=T), S['dopmargin'], tl, **ex, **kw_e(e))
        x_n = run(dict(NEW, p_dir=T), S['dopgas'], tl, **ex, **kw_e(not e))
        val[f'{b}|{"E" if e else "N"}|T={T}'] = dict(npv15=x_g['npv15'], dop_margin_vs_gas_npv15=round(x_m['npv15']-x_g['npv15'], 2),
                                                     edti_toggle_npv15=round(x_g['npv15']-x_n['npv15'], 2))
OUT['values'] = val
print('\nценность (NPV15, $ млн) при T_LOI: DoP-маржа против DoP-газа; EDTI вкл/выкл:', val)
json.dump(OUT, open(os.path.join(HERE, 'summary_r2_out.json'), 'w'), ensure_ascii=False, indent=1)
print('\nзаписано', os.path.join(HERE, 'summary_r2_out.json'))
