# Раунд 2: матрица безубыточного тарифа, правила ворот, порог якоря, условие по часам сети.
# Запуск: python3 matrix_r2.py  -> печать + matrix_r2_out.json. Все числа [О] на модели deal_r2.py,
# входы с метками [П]/[В]/[Д] -- в deal_r2.py и в комментариях ниже.
import os, sys, json, math
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from deal_r2 import *

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = {}
S = STRUCTS; TL = TL_NEW
USD = lambda n: round(n/FX, 3)
def kw_e(e): return dict(edti=e, depyrs=5 if e else 10)
def ceil5(x): return int(math.ceil(x/5.0)*5)
GAS_OPFEE_MMBTU = 0.15/1.05            # плата оператора 0.15 $/Mscf в пересчёте на MMBtu
def gas_from_allin(a): return a-GAS_OPFEE_MMBTU   # «всё включено» $/MMBtu -> цена газа без платы оператора
ALLIN_BASE = 0.60+GAS_OPFEE_MMBTU      # 0.743 $/MMBtu всё включено в базе модели
GAS_STOP_ALLIN = 1.5                   # стоп по газу «всё включено» [правило судьи, Д]

# ------------------------------------------------------------------ 1. Каскад исправлений
print('=== 1. Каскад исправлений (площадка 5 MMscf/д, 25 км, якорь 11 МВт, стройка 1 г), BE ₦/кВт·ч NPV15=0 / NPV20=0')
steps = [
    ('R1 судья (4.64; доступно 0.85, падение 6%)', OLD, TL_OLD),
    ('+ 4.40 МВт на MMscf/д', dict(OLD, mw_per=4.40), TL_OLD),
    ('+ газ VIIRS: доступно 0.78, падение 8%/г', NEW, TL_OLD),
    ('+ курсовой лаг: квартальный пересчёт + дебиторка', NEW, TL_NEW),
]
casc = {}
for nm, p, tl in steps:
    row = {}
    for s in ('dopmargin', 'dopgas', 'esc', 'none'):
        for e in (True, False):
            row[f'{s}|{"EDTI" if e else "noEDTI"}'] = (round(be(p, S[s], tl, 'npv15', **kw_e(e))), round(be(p, S[s], tl, 'npv20', **kw_e(e))))
    casc[nm] = row
    print(f'{nm:52s}', ' '.join(f'{k}:{v[0]}/{v[1]}' for k, v in row.items()))
OUT['cascade'] = casc
st_old = station(OLD); st_new = station(NEW)
OUT['gas_volumes'] = dict(DCQ_old=round(st_old['g_pow'], 3), DCQ_new=round(st_new['g_pow'], 3),
                          burned_new_y2=round(run(dict(NEW, p_dir=400), S['dopgas'], TL)['rows'][1]['burned'], 3),
                          kwh_per_mscf=round(st_new['kwh_per_mscf'], 1), capex=round(st_new['C'], 2),
                          capex_parts={k: round(v, 2) for k, v in st_new['capex'].items()})
print('газ: DCQ 11 МВт', OUT['gas_volumes'])

# ------------------------------------------------------------------ 2. Основная матрица
print('\n=== 2. Матрица BE ₦/кВт·ч NPV15=0 / NPV20=0 (газ база $0.74/MMBtu всё включено)')
MAT = {}
for site_nm in ('25km', 'oredo30'):
    for delay in (1, 2):
        for anchor in (11.0, 8.5):
            for s in ('dopgas', 'esc', 'none', 'dopmargin'):
                for e in (True, False):
                    p = dict(site(site_nm), anchor=anchor); tl = dict(TL, delay=delay)
                    a = be(p, S[s], tl, 'npv15', **kw_e(e)); b = be(p, S[s], tl, 'npv20', **kw_e(e))
                    MAT[f'{site_nm}|{delay}y|{anchor}MW|{s}|{"EDTI" if e else "noEDTI"}'] = (round(a), round(b))
for site_nm in ('25km', 'oredo30'):
    for anchor in (11.0, 8.5):
        print(f'-- {site_nm}, якорь {anchor} МВт: структура | EDTI 1г | EDTI 2г | без EDTI 1г | без EDTI 2г')
        for s in ('dopgas', 'esc', 'none', 'dopmargin'):
            cells = [MAT[f'{site_nm}|{d}y|{anchor}MW|{s}|{e}'] for e in ('EDTI', 'noEDTI') for d in (1, 2)]
            print(f'   {SLABEL[s]:48s}', ' | '.join(f'{a}/{b} (${USD(a):.3f}/${USD(b):.3f})' for a, b in cells))
OUT['matrix'] = MAT
# капитал станции по площадкам
OUT['capex'] = {f'{sn}|{a}MW': round(station(dict(site(sn), anchor=a))['C'], 2) for sn in ('25km', 'oredo30') for a in (11.0, 8.5)}
print('капитал $ млн:', OUT['capex'])

# ------------------------------------------------------------------ 3. Таблица NPV/IRR/DSCR для базовой структуры
print('\n=== 3. База р.2 (DoP-газ + эскроу), 25 км, 11 МВт, стройка 1 г: NPV15; NPV20; IRR; мин. DSCR (60%/14%/7 лет)')
tab = {}
for pr in (316, 350, 377, 411, 450, 500):
    for e in (True, False):
        x = run(dict(NEW, p_dir=pr), S['dopgas'], TL, debt=(0.6, 0.14, 7), **kw_e(e))
        tab[f'{pr}|{"EDTI" if e else "noEDTI"}'] = dict(npv15=x['npv15'], npv20=x['npv20'], irr=x['irr'], min_dscr=x['min_dscr'], y1=x['y1'], y2=x['y2'])
    print(f'₦{pr} (${USD(pr):.3f})', tab[f'{pr}|EDTI'], '| без EDTI', tab[f'{pr}|noEDTI'])
OUT['npv_table_base'] = tab

# ------------------------------------------------------------------ 4. Цена газа: BE NPV20 и стоп
print('\n=== 4. Газ. BE NPV20 (₦) при цене газа «всё включено» $/MMBtu; база р.2, 25 км, 11 МВт, 1 г')
gas_tab = {}
for s in ('dopgas', 'esc', 'none'):
    for e in (True, False):
        for d in (1, 2):
            for allin in (ALLIN_BASE, 1.0, 1.5, 2.18+GAS_OPFEE_MMBTU):
                p = dict(NEW, gas=gas_from_allin(allin))
                gas_tab[f'{s}|{"EDTI" if e else "noEDTI"}|{d}y|{allin:.2f}'] = round(be(p, S[s], dict(TL, delay=d), 'npv20', **kw_e(e)))
for k, v in gas_tab.items():
    if '|1y|' in k: print('  ', k, v)
OUT['be20_vs_gas'] = gas_tab

def gas_point(T, p, s, tl, e, key='npv20'):
    """цена газа «всё включено» $/MMBtu, при которой NPV=0 при тарифе T"""
    lo, hi = 0.0, 15.0
    for _ in range(45):
        m = (lo+hi)/2
        if run(dict(p, p_dir=T, gas=gas_from_allin(m)), s, tl, **kw_e(e))[key] >= 0: lo = m
        else: hi = m
    return lo

# ------------------------------------------------------------------ 5. Правила: порог LOI, стоп по газу, пол дизельной оговорки
print('\n=== 5. Правила по структурам: порог LOI = NPV20=0 при газе-стопе $1.5 всё включено (округл. вверх до ₦5)')
RULES = {}
for site_nm in ('25km', 'oredo30'):
    for anchor in (11.0, 8.5):
        for d in (1, 2):
            for s in ('dopgas', 'esc', 'none', 'dopmargin'):
                for e in (True, False):
                    p = dict(site(site_nm), anchor=anchor); tl = dict(TL, delay=d)
                    T_stop = be(dict(p, gas=gas_from_allin(GAS_STOP_ALLIN)), S[s], tl, 'npv20', **kw_e(e))
                    T = ceil5(T_stop)
                    gp = gas_point(T, p, S[s], tl, e)
                    # пол дизельной оговорки: тариф со 2-го года, при котором NPV20=0 (газ на стопе)
                    lo, hi = 0.3, 1.0
                    for _ in range(40):
                        m = (lo+hi)/2
                        if run(dict(p, p_dir=T, gas=gas_from_allin(GAS_STOP_ALLIN)), S[s], tl, tariff_cut=(2, m), **kw_e(e))['npv20'] >= 0: hi = m
                        else: lo = m
                    floor = T*hi
                    x = run(dict(p, p_dir=T, gas=gas_from_allin(GAS_STOP_ALLIN)), S[s], tl, debt=(0.6, 0.14, 7), **kw_e(e))
                    RULES[f'{site_nm}|{anchor}MW|{d}y|{s}|{"EDTI" if e else "noEDTI"}'] = dict(
                        T_LOI=T, T_LOI_usd=USD(T), gas_point_npv20=round(gp, 2), diesel_floor=round(floor), floor_usd=USD(floor),
                        min_dscr=x['min_dscr'], npv15=x['npv15'], irr=x['irr'])
for k, v in RULES.items():
    if '|11.0MW|' in k: print('  ', k, v)
OUT['rules'] = RULES

# ------------------------------------------------------------------ 6. Earn-out держателю: на реально сожжённом газе
print('\n=== 6. Earn-out $0.30/Mscf реально сожжённого газа со 2-го года (путь А)')
eo = {}
for s in ('dopgas',):
    for e in (True, False):
        p = NEW
        T = RULES[f'25km|11.0MW|1y|{s}|{"EDTI" if e else "noEDTI"}']['T_LOI']
        base_ = run(dict(p, p_dir=T), S[s], TL, **kw_e(e))
        w = run(dict(p, p_dir=T), S[s], TL, earnout=(0.30, 2), **kw_e(e))
        old_basis = run(dict(p, p_dir=T, opfee=0.15+0.30), S[s], TL, **kw_e(e))   # как в раунде 1: на DCQ, все годы
        be_w = be(p, S[s], TL, 'npv20', earnout=(0.30, 2), **kw_e(e)); be_0 = be(p, S[s], TL, 'npv20', **kw_e(e))
        # потолок earn-out при T_LOI+25 и +50 (газ база)
        caps = {}
        for dT in (0, 25, 50):
            lo, hi = 0.0, 10.0
            if run(dict(p, p_dir=T+dT), S[s], TL, **kw_e(e))['npv20'] < 0: caps[dT] = 0.0; continue
            for _ in range(40):
                m = (lo+hi)/2
                if run(dict(p, p_dir=T+dT), S[s], TL, earnout=(m, 2), **kw_e(e))['npv20'] >= 0: lo = m
                else: hi = m
            caps[dT] = round(lo, 2)
        yr = [round(r['eo'], 3) for r in w['rows']]
        eo[f'{s}|{"EDTI" if e else "noEDTI"}'] = dict(T=T, burned_mmscfd_y2=round(w['rows'][1]['burned'], 3), eo_musd_by_year=yr,
            d_npv15=round(w['npv15']-base_['npv15'], 2), d_npv20=round(w['npv20']-base_['npv20'], 2),
            d_npv15_old_basis=round(old_basis['npv15']-base_['npv15'], 2), be20_shift=round(be_w-be_0, 1),
            cap_usd_per_mscf_at_T_plus={k: v for k, v in caps.items()})
        print('  ', s, 'EDTI' if e else 'noEDTI', eo[f'{s}|{"EDTI" if e else "noEDTI"}'])
OUT['earnout'] = eo

# ------------------------------------------------------------------ 7. Порог замера газа
print('\n=== 7. Минимальный доступный газ года 1 (MMscf/д) для NPV20>=0; 25 км, 11 МВт, 1 г, DoP-газ')
meas = {}
for e in (True, False):
    T = RULES[f'25km|11.0MW|1y|dopgas|{"EDTI" if e else "noEDTI"}']['T_LOI']
    for dT in (0, 25, 50):
        for allin in (ALLIN_BASE, GAS_STOP_ALLIN):
            for decl in (0.08, 0.12):
                lo, hi = 0.5, 10.0
                for _ in range(40):
                    m = (lo+hi)/2
                    p = dict(NEW, Q=m, avail=1.0, decl=decl, gas=gas_from_allin(allin))
                    if run(dict(p, p_dir=T+dT), S['dopgas'], TL, **kw_e(e))['npv20'] >= 0: hi = m
                    else: lo = m
                meas[f'{"EDTI" if e else "noEDTI"}|T={T+dT}|gas{allin:.2f}|decl{decl}'] = round(hi, 2)
for k, v in meas.items(): print('  ', k, v)
OUT['gas_measure_min'] = meas

# ------------------------------------------------------------------ 8. Длина трассы
print('\n=== 8. BE NPV20 (₦) от длины трассы ЛЭП; DoP-газ, 11 МВт, газ не ограничивает (Oredo-подобный факел)')
line = {}
for km in (10, 15, 20, 25, 30, 33, 35, 40):
    for e in (True, False):
        for d in (1, 2):
            p = dict(site('oredo30'), line_km=float(km), loss=0.03*km/25)
            line[f'{km}|{"EDTI" if e else "noEDTI"}|{d}y'] = round(be(p, S['dopgas'], dict(TL, delay=d), 'npv20', **kw_e(e)))
    print(f'  {km} км:', {k.split("|",1)[1]: v for k, v in line.items() if k.startswith(f"{km}|")})
OUT['be20_vs_line_km'] = line

# ------------------------------------------------------------------ 9. Порог якоря
print('\n=== 9. BE NPV20 (₦) от размера якоря, МВт нетто; 25 км и Oredo 30 км, стройка 1 и 2 г')
anc = {}
MWS = (5.0, 6.0, 7.0, 8.0, 8.5, 9.0, 10.0, 11.0, 12.0, 13.0)
for site_nm in ('25km', 'oredo30'):
    for s in ('dopgas', 'esc', 'none', 'dopmargin'):
        for e in (True, False):
            for d in (1, 2):
                for mw in MWS:
                    p = dict(site(site_nm), anchor=mw)
                    anc[f'{site_nm}|{s}|{"EDTI" if e else "noEDTI"}|{d}y|{mw}'] = round(be(dict(p, gas=gas_from_allin(GAS_STOP_ALLIN)), S[s], dict(TL, delay=d), 'npv20', **kw_e(e)))
OUT['be20_vs_anchor_gasstop'] = anc
# потолок цены покупателя на дизеле: (1-s)*Pd, s -- требуемая экономия [Д 10-20%, центр 15%]
LEGAL_MW = 6.0/0.72   # eligible customer >=6 МВт·ч/ч в среднем за 90 дней [В: Mondaq/NERC-R-001-2024] при LF 0.72 [Д]
def min_anchor(site_nm, s, e, d, ceiling):
    xs = [(mw, anc[f'{site_nm}|{s}|{e}|{d}y|{mw}']) for mw in MWS]
    if xs[-1][1] > ceiling: return None
    for (m1, b1), (m2, b2) in zip(xs, xs[1:]):
        if b1 > ceiling >= b2: return round(m1+(b1-ceiling)*(m2-m1)/(b1-b2), 1)
    return xs[0][0] if xs[0][1] <= ceiling else None
thr = {}
for site_nm in ('25km', 'oredo30'):
    for s in ('dopgas', 'esc', 'none', 'dopmargin'):
        for e in ('EDTI', 'noEDTI'):
            for d in (1, 2):
                r = {}
                for Pd in (530, 650, 770):
                    c = 0.85*Pd
                    m = min_anchor(site_nm, s, e, d, c)
                    r[f'Pd{Pd}'] = None if m is None else max(m, round(LEGAL_MW, 1))
                    r[f'Pd{Pd}_econ'] = m
                thr[f'{site_nm}|{s}|{e}|{d}y'] = r
for k, v in thr.items():
    if k.startswith('25km') or '|dopgas|' in k: print('  ', k, v)
OUT['anchor_threshold'] = thr
OUT['legal_floor_MW'] = round(LEGAL_MW, 2)
print('  юридический пол якоря (6 МВт·ч/ч / LF 0.72):', round(LEGAL_MW, 2), 'МВт')
for site_nm in ('25km',):
    for s in ('dopgas',):
        for e in ('EDTI', 'noEDTI'):
            print('  кривая', site_nm, s, e, '1г:', [(mw, anc[f'{site_nm}|{s}|{e}|1y|{mw}']) for mw in MWS])
            print('  кривая', site_nm, s, e, '2г:', [(mw, anc[f'{site_nm}|{s}|{e}|2y|{mw}']) for mw in MWS])

# ------------------------------------------------------------------ 10. Часы сети: смешанная цена покупателя
print('\n=== 10. Максимум часов сети в сутки h, при котором покупатель проходит: T <= (1-s)*[(h*Pg+(24-h)*Pd)/24]')
def hmax(T, Pd, Pg=209.5, s=0.0):
    v = (Pd-T/(1-s))*24/(Pd-Pg)
    return round(max(0.0, min(24.0, v)), 1)
grid = {}
for key in ('25km|11.0MW|1y', '25km|11.0MW|2y', 'oredo30|11.0MW|1y', 'oredo30|11.0MW|2y', '25km|8.5MW|1y', 'oredo30|8.5MW|1y'):
    for s in ('dopgas', 'esc', 'none'):
        for e in ('EDTI', 'noEDTI'):
            T = RULES[f'{key}|{s}|{e}']['T_LOI']
            g = {}
            for Pd in (530, 650, 770):
                for sv in (0.0, 0.15):
                    for Pg in (209.5, 225.0):
                        g[f'Pd{Pd}|s{int(sv*100)}|Pg{Pg}'] = hmax(T, Pd, Pg, sv)
            grid[f'{key}|{s}|{e}'] = dict(T=T, h=g)
for k, v in grid.items():
    if '|11.0MW|1y' in k or '|11.0MW|2y' in k:
        h = v['h']; print(f'  {k:32s} T=₦{v["T"]}: h_max (s=0) Pd530/650/770 =', h['Pd530|s0|Pg209.5'], h['Pd650|s0|Pg209.5'], h['Pd770|s0|Pg209.5'],
                          '| s=15%:', h['Pd530|s15|Pg209.5'], h['Pd650|s15|Pg209.5'], h['Pd770|s15|Pg209.5'])
OUT['grid_hours'] = grid
# Band A 20 ч + дизель 4 ч
OUT['blend_bandA20'] = {Pd: round((20*209.5+4*Pd)/24) for Pd in (530, 650, 770)}
print('  Band A 20 ч + 4 ч дизеля, ₦/кВт·ч:', OUT['blend_bandA20'])

# ------------------------------------------------------------------ 11. Стрессы на пороге LOI базовой структуры
print('\n=== 11. Стрессы при T_LOI (DoP-газ, 25 км, 11 МВт, 1 г, газ на стопе $1.5): NPV20, $ млн')
stress = {}
for e in ('EDTI', 'noEDTI'):
    T = RULES[f'25km|11.0MW|1y|dopgas|{e}']['T_LOI']; kw = kw_e(e == 'EDTI')
    p = dict(NEW, p_dir=T, gas=gas_from_allin(GAS_STOP_ALLIN))
    cases = {
        'база (=0 по построению)': (p, TL, {}),
        'спрос /1.3 до контракта (нет ToP 70%)': (dict(p, LF=0.72/1.3), TL, {}),
        'спрос /2.7 до контракта': (dict(p, LF=0.72/2.7), TL, {}),
        'LF 0.60': (dict(p, LF=0.60), TL, {}),
        'капитал x2.8': (dict(p, cm=2.8), TL, {}),
        'капитал x2.4 (внешние $1 500-1 667/кВт)': (dict(p, cm=2.4), TL, {}),
        'падение дебита 12%/г': (dict(p, decl=0.12), TL, {}),
        'VIIRS p25: доступно 0.6': (dict(p, avail=0.6), TL, {}),
        'блокада 3 мес во 2-й год': (p, TL, dict(ev={2: 0.75})),
        'дизельная оговорка -20% со 2-го года': (p, TL, dict(tariff_cut=(2, 0.8))),
        'девальвация 16%/г (курсовой лаг)': (p, dict(TL, dep=0.16), {}),
        'стройка 2 г, капитал 50/50 по годам': (p, dict(TL, delay=2), dict(capex_split=[0.5, 0.5])),
        'стройка 2 г, весь капитал в t0 (срыв графика)': (p, dict(TL, delay=2), {}),
        'тариф в найре без индексации, девальвация 8%': (p, dict(TL, idx=0.0, dep=0.08), {}),
        'тариф в найре без индексации, девальвация 12%': (p, dict(TL, idx=0.0, dep=0.12), {}),
        'тариф в найре без индексации, девальвация 16%': (p, dict(TL, idx=0.0, dep=0.16), {}),
        'earn-out $0.30/Mscf сожжённого газа со 2-го года': (p, TL, dict(earnout=(0.30, 2))),
    }
    for nm, (pp, tl, ex) in cases.items():
        x = run(pp, S['dopgas'], tl, **kw, **ex)
        stress[f'{e}|{nm}'] = dict(npv15=x['npv15'], npv20=x['npv20'], irr=x['irr'])
        print(f'  {e:6s} T=₦{T} {nm:52s} NPV15 {x["npv15"]:+.1f}  NPV20 {x["npv20"]:+.1f}  IRR {x["irr"]}')
OUT['stress_at_T_LOI'] = stress
# BE в найре без индексации (тариф фиксирован в найре)
naira = {}
for dep in (0.08, 0.12, 0.16):
    for e in ('EDTI', 'noEDTI'):
        naira[f'dep{dep}|{e}'] = round(be(dict(NEW, gas=gas_from_allin(GAS_STOP_ALLIN)), S['dopgas'], dict(TL, idx=0.0, dep=dep), 'npv20', hi=5000, **kw_e(e == 'EDTI')))
print('  BE NPV20 при тарифе в найре без индексации, ₦/кВт·ч года 1:', naira)
OUT['be20_naira_unindexed'] = naira

# ------------------------------------------------------------------ 12. Поправки аудита NaCN: проверка кратностей
scr = dict(NEW, cm=1.0, red=1.0, dev=0.4, line_km=12.0, usd_km=70e3, subst=0.5,
           staff_core=0.30, sec=0.15, comm=0.05, ovh=0.10, tor=0.03, ins=0.01)
st_s = station(scr); st_b = station(NEW)
fs = st_s['rows'][0]['fixed']; fb = st_b['rows'][0]['fixed']
OUT['audit_check'] = dict(capex_screen=round(st_s['C'], 2), capex_base=round(st_b['C'], 2), capex_ratio=round(st_b['C']/st_s['C'], 2),
                          power_usd_per_kw_base=round(st_b['capex']['power']*1e6/(st_b['mw_gross']*1000)),
                          fixed_cash_screen=round(fs, 2), fixed_cash_base=round(fb, 2), fixed_ratio=round(fb/fs, 2),
                          cash_cost_ratio_incl_gas=round((fb+st_b['rows'][0]['gas'])/(fs+st_s['rows'][0]['gas']), 2))
print('\n=== 12. Проверка поправок аудита:', OUT['audit_check'])

json.dump(OUT, open(os.path.join(HERE, 'matrix_r2_out.json'), 'w'), ensure_ascii=False, indent=1)
print('\nзаписано', os.path.join(HERE, 'matrix_r2_out.json'))

# ------------------------------------------------------------------ 13. Дополнения
print('\n=== 13а. Стройка 2 г при капитале 50/50 по годам (плановая стройка, не срыв): BE NPV15/NPV20')
spl = {}
for site_nm in ('25km', 'oredo30'):
    for anchor in (11.0, 8.5):
        for s in ('dopgas', 'esc', 'none'):
            for e in (True, False):
                p = dict(site(site_nm), anchor=anchor); tl = dict(TL, delay=2)
                spl[f'{site_nm}|{anchor}MW|{s}|{"EDTI" if e else "noEDTI"}'] = (
                    round(be(p, S[s], tl, 'npv15', capex_split=[0.5, 0.5], **kw_e(e))),
                    round(be(p, S[s], tl, 'npv20', capex_split=[0.5, 0.5], **kw_e(e))),
                    ceil5(be(dict(p, gas=gas_from_allin(GAS_STOP_ALLIN)), S[s], tl, 'npv20', capex_split=[0.5, 0.5], **kw_e(e))))
for k, v in spl.items(): print('  ', k, 'BE15/BE20/T_LOI', v)
OUT['matrix_2y_split'] = spl

print('\n=== 13б. Порог пути А (СП с держателем NGFCP, earn-out $0.30/Mscf сожжённого газа со 2-го года), газ на стопе')
pa = {}
for site_nm in ('25km', 'oredo30'):
    for anchor in (11.0, 8.5):
        for d in (1, 2):
            for e in (True, False):
                p = dict(site(site_nm), anchor=anchor, gas=gas_from_allin(GAS_STOP_ALLIN))
                pa[f'{site_nm}|{anchor}MW|{d}y|dopgas|{"EDTI" if e else "noEDTI"}'] = ceil5(be(p, S['dopgas'], dict(TL, delay=d), 'npv20', earnout=(0.30, 2), **kw_e(e)))
for k, v in pa.items(): print('  ', k, 'T_LOI путь А =', v, f'(${USD(v):.3f})', '| путь Б =', RULES[k]['T_LOI'])
OUT['T_LOI_pathA'] = pa

print('\n=== 13в. Потолок earn-out $/Mscf сожжённого газа (NPV20=0) при тарифе T_LOI+dT и цене газа «всё включено»')
cap = {}
for e in (True, False):
    T0 = RULES[f'25km|11.0MW|1y|dopgas|{"EDTI" if e else "noEDTI"}']['T_LOI']
    for dT in (0, 25, 50):
        for allin in (ALLIN_BASE, 1.0, GAS_STOP_ALLIN):
            p = dict(NEW, p_dir=T0+dT, gas=gas_from_allin(allin))
            if run(p, S['dopgas'], TL, **kw_e(e))['npv20'] < 0: cap[f'{"EDTI" if e else "noEDTI"}|T={T0+dT}|gas{allin:.2f}'] = 0.0; continue
            lo, hi = 0.0, 15.0
            for _ in range(40):
                m_ = (lo+hi)/2
                if run(p, S['dopgas'], TL, earnout=(m_, 2), **kw_e(e))['npv20'] >= 0: lo = m_
                else: hi = m_
            cap[f'{"EDTI" if e else "noEDTI"}|T={T0+dT}|gas{allin:.2f}'] = round(lo, 2)
for k, v in cap.items(): print('  ', k, v)
OUT['earnout_cap'] = cap

print('\n=== 13г. ToP газа: 100% DCQ (база) против 70% DCQ (back-to-back с ToP 70% PPA): BE NPV20')
tp = {}
for e in (True, False):
    for allin in (ALLIN_BASE, GAS_STOP_ALLIN):
        p = dict(NEW, gas=gas_from_allin(allin))
        a = be(p, S['dopgas'], TL, 'npv20', **kw_e(e)); b = be(p, S['dopgas'], TL, 'npv20', top=0.70, **kw_e(e))
        tp[f'{"EDTI" if e else "noEDTI"}|gas{allin:.2f}'] = (round(a), round(b))
print('  ', tp); OUT['top_share'] = tp

print('\n=== 13д. Экономия покупателя при T_LOI (11 МВт, 56.7 ГВт·ч/г с простоем 1 мес; 61.9 без простоя), $ млн/г')
gwh = run(dict(NEW, p_dir=400), S['dopgas'], TL)['rows'][1]
GWH_FULL = station(NEW)['rows'][0]['gwh']; GWH_UP = GWH_FULL*11/12
sav = {}
for e in ('EDTI', 'noEDTI'):
    T = RULES[f'25km|11.0MW|1y|dopgas|{e}']['T_LOI']
    for nm, alt in (('дизель 530', 530), ('дизель 650', 650), ('дизель 770', 770),
                    ('сеть 12 ч + дизель 650', (12*209.5+12*650)/24), ('Band A 20 ч + дизель 650', (20*209.5+4*650)/24)):
        sav[f'{e}|T={T}|{nm}'] = dict(alt=round(alt), saving_ngn_kwh=round(alt-T), saving_musd=round((alt-T)*GWH_UP/FX, 2))
for k, v in sav.items(): print('  ', k, v)
OUT['buyer_saving'] = sav
OUT['gwh'] = dict(full=round(GWH_FULL, 1), with_outage=round(GWH_UP, 1))

print('\n=== 13е. Правила для якоря 8.5 МВт')
for k, v in RULES.items():
    if '|8.5MW|' in k and ('dopgas' in k or 'esc' in k): print('  ', k, v)

json.dump(OUT, open(os.path.join(HERE, 'matrix_r2_out.json'), 'w'), ensure_ascii=False, indent=1)
print('\nзаписано (с дополнениями)', os.path.join(HERE, 'matrix_r2_out.json'))
