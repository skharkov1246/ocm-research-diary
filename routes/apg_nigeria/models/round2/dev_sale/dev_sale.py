# Раунд 2. Девелоперская модель: довести площадку до готовности и продать стратегу.
# Запуск: python3 dev_sale.py  -> печать + dev_sale_out.json
# Все расчётные числа [О]. Входы помечены: [П] первичный, [В] вторичный, [Д] допущение (диапазон).
# Экономика станции, которую покупает стратег, берётся из ../deal/deal_r2.py (исправленная модель раунда 2:
# 4.40 МВт на MMscf/д, загрузка факела 0.78 по VIIRS, падение 8%/г, курсовой лаг, стройка 1-2 г, налог 34%).
import os, sys, json, math, csv
from collections import defaultdict
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'deal'))
from deal_r2 import NEW, STRUCTS, TL_NEW, run, be, site, station, FX, npv

OUT = {}
S = STRUCTS; TL = TL_NEW
def kw_e(e): return dict(edti=e, depyrs=5 if e else 10)
GAS_OPFEE_MMBTU = 0.15/1.05
def gas_from_allin(a): return a-GAS_OPFEE_MMBTU

# ---------------------------------------------------------------- 0. Именованные площадки (World Bank VIIRS)
CSV = os.path.join(HERE, '..', '..', 'wb', 'nigeria_flares.csv')
hubs = {'Warri': (5.52, 5.75), 'Ughelli': (5.49, 5.99), 'Sapele': (5.89, 5.68), 'Benin': (6.34, 5.63)}
def hav(a, b, c, e):
    R = 6371; p = math.radians
    h = math.sin(p(c-a)/2)**2+math.cos(p(a))*math.cos(p(c))*math.sin(p(e-b)/2)**2
    return 2*R*math.asin(math.sqrt(h))
hist = defaultdict(dict); meta = {}
for r in csv.DictReader(open(CSV)):
    try: hist[r['id']][int(r['Year'])] = hist[r['id']].get(int(r['Year']), 0)+float(r['MMSCFD'])
    except Exception: continue
    meta[r['id']] = r
named = []
for fid, h in hist.items():
    r = meta[fid]
    if r['Location'] != 'ONSHORE': continue
    la, lo = float(r['lat']), float(r['lon'])
    ds = {k: hav(a, b, la, lo) for k, (a, b) in hubs.items()}
    hub = min(ds, key=ds.get)
    q25 = h.get(2025, 0); mn = min(h.get(y, 0) for y in (2023, 2024, 2025))
    named.append(dict(name=r['Field name'], op=r['Operator'], q2025=round(q25, 2), min23_25=round(mn, 2),
                      hub=hub, km=round(ds[hub], 1)))
# класс 1 MMscf/д: стабильные 0.8-1.6 (минимум 2023-25 >= 0.7) до 20 км от промузлов Delta/Edo
c1 = sorted([n for n in named if 0.8 <= n['q2025'] <= 1.6 and n['min23_25'] >= 0.7 and n['km'] <= 20], key=lambda n: n['km'])
c5 = [n for n in named if n['name'] in ('Ughelli East',)]
c15 = [n for n in named if n['name'] in ('Oredo', 'Sapele', 'Oben')]
OUT['named'] = dict(class1=c1[:8], class5=c5, class15=c15)
print('=== 0. Именованные площадки [О, VIIRS World Bank 2012-2025]')
for k, v in OUT['named'].items():
    for n in v: print(f'  {k}: {n}')

# ---------------------------------------------------------------- 1. Стоимость довести до готовности
# Бюджеты ворот 0/1/2 -- из решения судьи (сумма статей) [Д по статьям]; сверху -- то, что судья вынес за скобки.
GATES = dict(
    g0=(0.047, 0.112),        # реестр, юрист, данные NGFCP, котировки, MAIN [Д; данные NGFCP $2k+$1k/площадка П]
    g1=(0.110, 0.270),        # term sheet СП/опцион или GCA, налог. заключение, расходомер+пробы, трасса, община [Д]
    g2=(0.272, 0.442),        # 6 мес. замера, PPA, FEED с твёрдыми котировками, кредитор, лицензия NMDPRA [Д]
)
TEAM = dict(                  # команда 2 чел. (коммерция+юрист/инженер), поездки, охрана полевых работ [Д]
    g0=(0.030, 0.080),        # 2 мес.: ~$15-40k/мес [Д]
    g1=(0.060, 0.160),        # 4 мес.
    g2=(0.090, 0.240),        # 6 мес.
)
SALE_COST = (0.04, 0.08)      # юристы/консультант/дата-рум по сделке продажи, доля цены, мин. $30k [Д]
def stage_cost(g, k):         # k=0 низ, 1 верх
    return GATES[g][k]+TEAM[g][k]
dev_cost = {lvl: (round(sum(stage_cost(g, 0) for g in gs), 3), round(sum(stage_cost(g, 1) for g in gs), 3))
            for lvl, gs in (('R1', ('g0', 'g1')), ('R2', ('g0', 'g1', 'g2')))}
OUT['dev_cost_musd'] = dev_cost
print('\n=== 1. Стоимость готовности, $ млн [О: сумма статей]:', dev_cost)
# 1 MMscf/д: те же фиксированные статьи (юрист, реестр, расходомер), FEED дешевле; берём R1 как у 5 [Д].

# ---------------------------------------------------------------- 2. Сколько стоит готовая площадка стратегу
# Покупатель-стратег энергии (Genesis, Elektron, Globeleq/CPGNL, BlueCore/Axxela -- тип, не подтверждённый интерес):
#   - капитал x2.2 вместо x2.6 (свой парк и закупки) [Д 2.1-2.6; внешние точки x2.1-2.4 по генсетам, В];
#   - персонал без экспата: множитель 1.4 вместо 1.8 [Д 1.3-1.8];
#   - ставка 15% [Д 12-18%];
#   - остаток девелопмента внутри капитала у покупателя 0.3 вместо 0.8 $ млн (FEED/лицензии сделаны нами) [Д].
BUY = dict(cm=2.2, staff_mult=1.4, dev=0.3)
OWN = dict()   # наша стройка = модель раунда 2 как есть (cm 2.6, staff_mult 1.8, dev 0.8)
SITES = {
    '5_UghelliEast_25km': dict(site('25km')),                     # Q=5, трасса 25 км
    '15_Oredo_30km_Q8.5': dict(site('oredo30')),                  # Oredo: минимум 2023-25 8.5, трасса 30 км
    '15_generic_25km': dict(site('25km'), Q=15.0),                # условный факел 15, якорь тот же 11 МВт
}
TARIFFS = [350, 377, 405, 440, 470]   # ₦/кВт·ч в долларовом эквиваленте при ₦1 330/$
def val(p, T, e, delay, extra):
    pp = dict(p, p_dir=T, gas=gas_from_allin(1.0), **extra)    # газ $1.0/MMBtu всё включено (цель GSA) [Д]
    tl = dict(TL, delay=delay)
    kw = kw_e(e)
    if delay == 2: kw['capex_split'] = [0.5, 0.5]
    return run(pp, S['dopgas'], tl, **kw)
buyer = {}
print('\n=== 2. NPV станции для покупателя-стратега и для нас ($ млн), DoP возвращает газ + эскроу, газ $1.0 всё вкл.')
print('   site | T | EDTI | стройка | capex наш/стратег | NPV15 наш | NPV15 стратег | NPV20 стратег | IRR стратег')
for sn, p in SITES.items():
    for T in TARIFFS:
        for e in (True, False):
            for d in (1, 2):
                o = val(p, T, e, d, OWN); b = val(p, T, e, d, BUY)
                buyer[f'{sn}|{T}|{"E" if e else "N"}|{d}y'] = dict(capex_own=o['capex'], capex_buy=b['capex'],
                    npv15_own=o['npv15'], npv20_own=o['npv20'], npv15_buy=b['npv15'], npv20_buy=b['npv20'],
                    irr_buy=b['irr'], irr_own=o['irr'], ebitda_y1_own=o['y1'], ebitda_y2_own=o['y2'])
                if e and d in (1, 2) and T in (350, 405, 440):
                    print(f'   {sn} | {T} | {"E" if e else "N"} | {d}y | {o["capex"]}/{b["capex"]} | {o["npv15"]} | {b["npv15"]} | {b["npv20"]} | {b["irr"]}')
OUT['buyer_value'] = buyer

# ---------------------------------------------------------------- 3. Премия за девелопмент
# Ориентиры рынка:
#   - «типичная плата за девелопмент 3-5% капзатрат», выплата на финансовом закрытии [В: Trinity LLP, практикующий юрист, не обзор рынка];
#   - RTB-проект малой распределённой генерации, Бразилия: ~EUR100 тыс./МВт против ~EUR1.3 млн/МВт действующего [В: Mercom, Enerside->IVI, слабый аналог];
#   - <30% ВИЭ-проектов в Африке доходят до финансового закрытия [В: Castalia, Engineering News 2016];
#   - Ardova->Powergas: цена не раскрыта [В]; Globeleq->CPGNL 74% (12 газ. станций, 58 МВт): цена не раскрыта [В];
#     BlueCore->Axxela: долговое финансирование покупки $285 млн [В], цена/МВт не выделяется.
#   => прецедента «готовая факельная площадка продана в Нигерии» НЕ НАЙДЕНО.
# Модель премии: P = min(alpha * NPV15_стратега, cap * capex_стратега), не ниже 0.
#   alpha -- доля стоимости, которую забирает девелопер [Д 0.2-0.5; R1 (до PPA) 0.1-0.25];
#   cap   -- 3-5% капзатрат (R2, подписанный PPA) [В/Д]; для R1 cap 1-2% [Д].
SC = {   # пессимизм / база / оптимизм
    'pess': dict(alpha_R1=0.10, alpha_R2=0.20, cap_R1=0.010, cap_R2=0.03, p0=0.10, p1=0.40, p2=0.35,
                 ps_R1=0.15, ps_R2=0.30, pfc=0.50, upfront=0.30, cost_k=1, sale_cost=SALE_COST[1], bypass=0.6),
    'base': dict(alpha_R1=0.175, alpha_R2=0.35, cap_R1=0.015, cap_R2=0.04, p0=0.20, p1=0.55, p2=0.50,
                 ps_R1=0.30, ps_R2=0.50, pfc=0.65, upfront=0.40, cost_k=0.5, sale_cost=sum(SALE_COST)/2, bypass=0.45),
    'opt':  dict(alpha_R1=0.25, alpha_R2=0.50, cap_R1=0.020, cap_R2=0.05, p0=0.30, p1=0.70, p2=0.65,
                 ps_R1=0.45, ps_R2=0.70, pfc=0.80, upfront=0.50, cost_k=0, sale_cost=SALE_COST[0], bypass=0.3),
}
# p0 -- LOI якоря >=8.3-10 МВт на дизеле у названного факела по цене >= порога [Д; в данных 0 подписанных покупателей
#       у 9 факелов, пул без трубы 12-70 МВт на весь юго-юг/юго-восток -- поэтому 0.1-0.3];
# p1 -- опцион/СП с держателем NGFCP или GCA с NEPL + согласие NUPRC [Д];
# p2 -- 6 мес. замера >= 3.0 MMscf/д и подписанный PPA [Д];
# ps -- стратег подписывает покупку при готовности [Д]; pfc -- финансовое закрытие после покупки [Д; Castalia <30% от раннего этапа].
# Сроки (лет от 10.2026): R1 готов ~0.5, сделка ~1.0, FC ~2.0; R2 готов ~1.0, сделка ~1.3, FC ~2.0 [Д].
TIMES = dict(R1=dict(c=[0.1, 0.35], sale=1.0, fc=2.0), R2=dict(c=[0.1, 0.35, 0.75], sale=1.3, fc=2.0))
CAPDEV = 0.15        # ставка дисконтирования девелопера для отчёта NPV15; NPV20 -- то же при 20%
NAIRA_DEP = 0.12     # если отложенная часть цены в найре: девальвация 8-16%/г, центр 12% [Д]

def premium(lvl, sc, bv):
    a = sc[f'alpha_{lvl}']; cap = sc[f'cap_{lvl}']
    return max(0.0, min(a*bv['npv15_buy'], cap*bv['capex_buy']))

def dev_ev(lvl, sc, bv, r, k_cost, naira_deferred=False, bypass=False):
    """Ожидаемая стоимость девелопера, $ млн. k_cost: 0 = нижняя граница затрат, 1 = верхняя."""
    gs = ['g0', 'g1'] if lvl == 'R1' else ['g0', 'g1', 'g2']
    probs = [1.0, sc['p0']] + ([sc['p0']*sc['p1']] if lvl == 'R2' else [])   # вероятность дойти до трат стадии
    t = TIMES[lvl]
    ev = 0.0
    for g, pr, tt in zip(gs, probs, t['c']):
        c = stage_cost(g, 0)*(1-k_cost)+stage_cost(g, 1)*k_cost
        ev -= pr*c/(1+r)**tt
    p_ready = sc['p0']*sc['p1']*(sc['p2'] if lvl == 'R2' else 1.0)
    P = premium(lvl, sc, bv)
    if bypass: P *= (1-sc['bypass'])              # данные/LOI раскрыты до эксклюзива: стратег идёт к якорю/держателю напрямую
    basis = sum(stage_cost(g, 0)*(1-k_cost)+stage_cost(g, 1)*k_cost for g in gs)
    up = P*sc['upfront']; df = P*(1-sc['upfront'])
    tx_up = max(0.0, (up-basis)*0.34)              # налог 34% на прибыль от продажи; убытки без щита [Д]
    sale_c = max(0.03, P*sc['sale_cost']) if P > 0 else 0.0
    fc_mult = (1-NAIRA_DEP)**(t['fc']-t['sale']) if naira_deferred else 1.0
    tx_df = max(0.0, (df*fc_mult - max(0.0, basis-up))*0.34)
    ps = sc[f'ps_{lvl}']
    ev += p_ready*ps*((up-tx_up-sale_c)/(1+r)**t['sale'] + sc['pfc']*(df*fc_mult-tx_df)/(1+r)**t['fc'])
    return ev, P, p_ready*ps, p_ready*ps*sc['pfc']

res = {}
print('\n=== 3. Девелопер: премия и ожидаемая стоимость ($ млн), тариф LOI = T, EDTI подтверждён, стройка у покупателя 1 г')
for sn in SITES:
    for T in (350, 405, 440):
        bv = buyer[f'{sn}|{T}|E|1y']
        for lvl in ('R1', 'R2'):
            for scn, sc in SC.items():
                kc = sc['cost_k']
                ev15, P, psale, pfc = dev_ev(lvl, sc, bv, 0.15, kc)
                ev20, _, _, _ = dev_ev(lvl, sc, bv, 0.20, kc)
                evn, _, _, _ = dev_ev(lvl, sc, bv, 0.15, kc, naira_deferred=True)
                evb, _, _, _ = dev_ev(lvl, sc, bv, 0.15, kc, bypass=True)
                # аукцион 3-4 стратегов без потолка 3-5%: премия = alpha_R2 x NPV15 стратега [Д, верхняя оценка]
                sa = dict(sc, cap_R1=1.0, cap_R2=1.0)
                eva, Pa, _, _ = dev_ev(lvl, sa, bv, 0.15, kc)
                key = f'{sn}|T{T}|{lvl}|{scn}'
                res[key] = dict(premium=round(P, 3), p_sale=round(psale, 3), p_fc=round(pfc, 3),
                                ev15=round(ev15, 3), ev20=round(ev20, 3), ev15_naira_deferred=round(evn, 3),
                                ev15_bypass_open_data=round(evb, 3), npv15_buy=bv['npv15_buy'], capex_buy=bv['capex_buy'],
                                premium_auction=round(Pa, 2), ev15_auction=round(eva, 3))
                if scn == 'base' or T == 405:
                    print(f'  {key:42s} premium {P:5.2f}  P(sale) {psale:.3f}  P(FC) {pfc:.3f}  EV15 {ev15:+.3f}  EV20 {ev20:+.3f}  '
                          f'найра {evn:+.3f}  обход {evb:+.3f}  аукцион P {Pa:.2f} EV15 {eva:+.3f}  (NPV15 стратега {bv["npv15_buy"]})')
OUT['dev_ev'] = res

# ---------------------------------------------------------------- 4. 1 MMscf/д: есть ли что продавать
# Энергия: 1 x 0.78 x 4.40 x 0.96 = 3.3 МВт нетто < юридического порога 8.3 МВт якоря (раунд 2) -> законного пути нет.
mw1 = 1.0*0.78*4.40*0.96
# Газ в трубу (round2/gas_to_pipe): Q1 при 1-8 км NPV15 от -6.8 до -38 во всех сценариях -> стратегу стоимость 0.
# CNG (cng_skeptic, cng_y2): 1 MMscf/д NPV15 base -10.4, high -2.6 -> 0.
g2p = json.load(open(os.path.join(HERE, '..', 'gas_to_pipe', 'g2p_out.json')))
OUT['q1'] = dict(mw_net=round(mw1, 2), legal_floor_mw=8.33, power_value=0.0,
                 note='газ в трубу и CNG на 1 MMscf/д отрицательны во всех сценариях -> премия 0, EV = -затраты ворот 0')
for sc in SC.values():
    pass
ev_q1 = {k: round(-(stage_cost('g0', 0)*(1-kk)+stage_cost('g0', 1)*kk), 3) for k, kk in (('opt', 0), ('base', 0.5), ('pess', 1))}
OUT['q1']['ev_if_started'] = ev_q1
print('\n=== 4. 1 MMscf/д: МВт нетто', round(mw1, 2), '< 8.33; EV если начать ворота 0:', ev_q1)

# ---------------------------------------------------------------- 5. Пороги: при чём EV девелопера = 0
def solve(f, lo, hi, it=50):
    if f(lo)*f(hi) > 0: return None
    for _ in range(it):
        m = (lo+hi)/2
        if f(lo)*f(m) <= 0: hi = m
        else: lo = m
    return (lo+hi)/2
th = {}
for sn in ('5_UghelliEast_25km', '15_Oredo_30km_Q8.5'):
    for lvl in ('R1', 'R2'):
        bv = buyer[f'{sn}|405|E|1y']
        sc = dict(SC['base'])
        # (а) минимальная премия для EV15=0
        def fP(P):
            s2 = dict(sc, **{f'alpha_{lvl}': 1.0, f'cap_{lvl}': P/bv['capex_buy']})
            return dev_ev(lvl, s2, dict(bv, npv15_buy=1e9), 0.15, 0.5)[0]
        Pmin = solve(fP, 0.0, 200.0)
        # (б) минимальная p0 при базовой премии
        fp0 = lambda x: dev_ev(lvl, dict(sc, p0=x), bv, 0.15, 0.5)[0]
        p0min = solve(fp0, 0.001, 1.0)
        # (б2) требуемая общая вероятность «готово и продано и закрыто» при базовой премии:
        #      P* = PV затрат / PV чистой премии (затраты все стадии как если бы шли до конца)
        gs = ['g0', 'g1'] if lvl == 'R1' else ['g0', 'g1', 'g2']
        pvc = sum(stage_cost(g, 0)*0.5+stage_cost(g, 1)*0.5 for g in gs)
        Pb = premium(lvl, sc, bv); pvp = Pb*0.66*0.9/1.15**1.5
        pstar = pvc/pvp if pvp > 0 else None
        # (в) минимальный тариф LOI, при котором EV20=0 (базовые вероятности)
        def fT(T):
            b2 = val(SITES[sn], T, True, 1, BUY)
            bvv = dict(npv15_buy=b2['npv15'], capex_buy=b2['capex'])
            return dev_ev(lvl, sc, bvv, 0.20, 0.5)[0]
        Tmin = solve(fT, 250.0, 700.0)
        # (г) порог по NPV стратега: тариф, при котором NPV15 стратега = 0
        Tb0 = be(dict(SITES[sn], gas=gas_from_allin(1.0), **BUY), S['dopgas'], TL, 'npv15', **kw_e(True))
        th[f'{sn}|{lvl}'] = dict(premium_min_musd=None if Pmin is None else round(Pmin, 3),
                                 p0_min=None if p0min is None else round(p0min, 3),
                                 T_min_EV20=None if Tmin is None else round(Tmin),
                                 T_buyer_npv15_0=round(Tb0), premium_base=round(Pb, 2),
                                 p_success_required_if_all_costs=None if pstar is None else round(pstar, 2),
                                 p_success_base=round(sc['p0']*sc['p1']*(sc['p2'] if lvl == 'R2' else 1)*sc[f'ps_{lvl}']*sc['pfc'], 3))
        print(f'  пороги {sn}|{lvl}: {th[f"{sn}|{lvl}"]}')
OUT['thresholds'] = th

# ---------------------------------------------------------------- 6. Сравнение с «энергией якорю» (своя стройка)
# Своя стройка: EV = P(ворота 0-2) * NPV станции - затраты ворот; капитал под риском $29-34 млн.
cmp = {}
for sn in ('5_UghelliEast_25km', '15_Oredo_30km_Q8.5'):
    for T in (350, 405, 440):
        for d in (1, 2):
            o = buyer[f'{sn}|{T}|E|{d}y']
            for scn, sc in SC.items():
                k = sc['cost_k']
                pr = sc['p0']*sc['p1']*sc['p2']*sc['pfc']      # свой FID тоже требует кредитора [Д тот же pfc]
                def gc(r):
                    probs = [1.0, sc['p0'], sc['p0']*sc['p1']]
                    return sum(pp*(stage_cost(g, 0)*(1-k)+stage_cost(g, 1)*k)/(1+r)**tt
                               for g, pp, tt in zip(('g0', 'g1', 'g2'), probs, TIMES['R2']['c']))
                ev_build15 = pr*o['npv15_own']/(1.15**1.0) - gc(0.15)   # FID ~10.2027 (t=1) [Д]
                ev_build20 = pr*o['npv20_own']/(1.20**1.0) - gc(0.20)
                dv = res.get(f'{sn}|T{T}|R2|{scn}')
                cmp[f'{sn}|T{T}|{d}y|{scn}'] = dict(npv15_build=o['npv15_own'], npv20_build=o['npv20_own'],
                    ev15_build=round(ev_build15, 2), ev20_build=round(ev_build20, 2),
                    ev15_dev_R2=None if dv is None else dv['ev15'], ev15_dev_R1=res.get(f'{sn}|T{T}|R1|{scn}', {}).get('ev15'),
                    capex_build=o['capex_own'], ebitda_y1_own=o['ebitda_y1_own'], ebitda_y2_own=o['ebitda_y2_own'], irr_own=o['irr_own'])
OUT['vs_build'] = cmp
print('\n=== 6. Своя стройка vs продажа на R2 (база, стройка 1 г, EDTI):')
for k, v in cmp.items():
    if '|1y|base' in k or '|2y|base' in k: print('  ', k, v)

# ---------------------------------------------------------------- 7. Вариант «девелопмент как услуга» держателю NGFCP
# Держатель платит за работу ворот 0-1 (cost-plus 15-30% [Д]); риск неплатежа держателя 20-50% [Д: ~22 из 42 отстают, первого газа нет].
svc = {}
for nm, (m, nonpay, kc) in dict(pess=(0.15, 0.5, 1), base=(0.22, 0.35, 0.5), opt=(0.30, 0.2, 0)).items():
    c = sum(stage_cost(g, 0)*(1-kc)+stage_cost(g, 1)*kc for g in ('g0', 'g1'))
    fee = c*(1+m)
    ev = fee*(1-nonpay) - c
    svc[nm] = dict(cost=round(c, 3), fee=round(fee, 3), ev_pre_tax=round(ev, 3), ev_after_tax=round(ev*(0.66 if ev > 0 else 1), 3))
OUT['service_fee'] = svc
print('\n=== 7. Девелопмент как услуга держателю (до налога/после):', svc)

# ---------------------------------------------------------------- 8. Продажа как страховка затрат: условно на «готово»
# (а) все вероятности = 1: покрывает ли премия 3-5% капзатрат сами затраты готовности?
# (б) условно на R2 достигнут: чистые поступления от продажи против уже потраченного.
cond = {}
for sn in ('5_UghelliEast_25km', '15_Oredo_30km_Q8.5'):
    bv = buyer[f'{sn}|405|E|1y']
    for scn, sc in SC.items():
        for lvl in ('R1', 'R2'):
            s1 = dict(sc, p0=1.0, p1=1.0, p2=1.0, ps_R1=1.0, ps_R2=1.0, pfc=1.0)
            ev_cert, P, _, _ = dev_ev(lvl, s1, bv, 0.15, sc['cost_k'])
            gs = ['g0', 'g1'] if lvl == 'R1' else ['g0', 'g1', 'g2']
            sunk = sum(stage_cost(g, 0)*(1-sc['cost_k'])+stage_cost(g, 1)*sc['cost_k'] for g in gs)
            t = TIMES[lvl]; ready_t = 0.5 if lvl == 'R1' else 1.0
            sale_c = max(0.03, P*sc['sale_cost'])
            up = P*sc['upfront']; df = P*(1-sc['upfront'])
            net_at_ready = (sc[f'ps_{lvl}']*((up-sale_c)/(1.15)**(t['sale']-ready_t) + sc['pfc']*df/(1.15)**(t['fc']-ready_t)))
            cond[f'{sn}|{lvl}|{scn}'] = dict(premium=round(P, 2), sunk=round(sunk, 3), ev15_all_prob_1=round(ev_cert, 3),
                                             expected_net_proceeds_given_ready=round(net_at_ready, 3),
                                             recovery_share=round(net_at_ready/sunk, 2))
OUT['conditional'] = cond
print('\n=== 8. Условно на готовность (T=₦405, EDTI, 1 г):')
for k, v in cond.items(): print('  ', k, v)

json.dump(OUT, open(os.path.join(HERE, 'dev_sale_out.json'), 'w'), ensure_ascii=False, indent=1)
print('\nOK ->', os.path.join(HERE, 'dev_sale_out.json'))
