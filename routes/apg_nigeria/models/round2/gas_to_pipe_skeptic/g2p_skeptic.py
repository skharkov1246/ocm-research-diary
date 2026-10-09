# Скептик к g2p.py («газ в трубу»). Проверяем, ломается ли вердикт «не делать», если снять
# методические перекосы модели в сторону минуса и добавить пропущенные плюсы. Метки: [О] расчёт, [Д] допущение.
# Запуск: python3 g2p_skeptic.py -> g2p_skeptic_out.json рядом.
import json, os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'gas_to_pipe'))
import g2p
HERE = os.path.dirname(os.path.abspath(__file__))
TAX = g2p.TAX

def run(Q, km, sc='base', ov=None, contract=True, fix=(), years=10):
    """Копия g2p.run с переключаемыми поправками:
       'up'   - не умножать на работу апстрима 0.92: среднегодовой VIIRS уже включает простои промысла [О]
       'own'  - готовность станции 0.98 при резерве 3x50% (N+1) вместо 0.95 [Д 0.97-0.99]
       'tor'  - ТОиР 4.25% только на смонтированный пакет (pkg x cm); труба 1%/г; страховка 1% от капитала [Д]
       'nol'  - перенос убытков для налога (NTA: убытки переносятся без срока) [Д]
       'size' - компрессоры под среднюю загрузку: Q*lf*1.5 вместо Q*1.5 [Д]"""
    p = dict(g2p.SC[sc]); p.update(ov or {})
    if 'up' in fix: p['up'] = 1.0
    if 'own' in fix: p['own'] = 0.98
    Qcap = Q*p['lf'] if 'size' in fix else Q
    cx = g2p.capex(Qcap, km, p)
    # труба не зависит от размера компрессоров: пересчитать по полному Q (диаметр)
    if 'size' in fix:
        cx_full = g2p.capex(Q, km, p); dp = cx_full['pipe'] - cx['pipe']
        cx['total'] = round(cx['total'] + dp*(1+0.14*p['build']/2), 2); cx['pipe'] = cx_full['pipe']
    C = cx['total']
    dem = 1.0 if contract else p['demand']
    ss = g2p.staff_scale(Q)
    if 'tor' in fix:
        main = cx['pkg']*p['cm']*p['tor'] + cx['pipe']*0.01 + C*0.01
    else:
        main = C*(p['tor']+p['ins'])
    fixed_bu = (p['staff']+p['ovh'])*ss + p['sec'] + 0.01*cx['pipe_len'] + p['comm'] + main
    fixed_screen = 0.20*ss + 0.10 + cx['screen']*(0.02+0.005)
    fixed = max(fixed_bu, fixed_screen*p['opx_screen_mult'])
    b = p['build']; cf = [-C/b]*b; prev = 0.0; loss = 0.0; e2 = None; e1 = None
    for t in range(1, years+1):
        Gin = Q*p['lf']*(1-p['decl'])**(t-1)*p['up']*p['own']*(p['ramp'] if t == 1 else 1.0)
        take = Gin*dem; sales = take*(1-g2p.SHRINK-0.01)
        price = (p['price']*(1+p['pdrift'])**(t-1) - p['transport'])*p['hv']
        rev = sales*1000*365*price/1e6*p['coll']*(1-p['fxloss']*p['dso']/365)
        opex = fixed*1.03**(t-1) + (Gin if p['top_all'] else take)*1000*365*p['gas_in']/1e6
        e = rev - opex
        if t == 1: e1 = e
        if t == 2: e2 = e
        ti = e - C/10
        if 'nol' in fix:
            if ti < 0: loss += -ti; tx = 0.0
            else:
                use = min(loss, ti); loss -= use; tx = (ti-use)*TAX
        else: tx = max(0.0, ti*TAX)
        wc = p['dso']/365*(rev-prev); prev = rev
        cf.append(e - tx - wc + (p['dso']/365*rev if t == years else 0))
    i = g2p.irr(cf)
    return dict(capex=C, fixed=round(fixed, 2), e1=round(e1, 2), e2=round(e2, 2), npv15=round(g2p.npv(.15, cf), 2),
                npv20=round(g2p.npv(.20, cf), 2), irr=None if i is None else round(i*100, 1))

def solve(fn, lo, hi):
    return g2p.solve(fn, lo, hi)

SITES = {'Ughelli East 6.7/3.7 км': (6.7, 3.7), 'Afam Umuosi 5.6/0.9 км': (5.6, 0.9), 'Oben 8.7/2 км': (8.7, 2.0),
         'типовой 5/4 км': (5.0, 4.0), 'типовой 15/4 км': (15.0, 4.0),
         # кластер Ughelli: East 6.7 + West 2.9 + Utorogu 4.3 (VIIRS 2025) = 13.9; сборные линии ~10 км + отвод 3.7 [Д]
         'кластер Ughelli 13.9/14 км': (13.9, 14.0)}
ALL = ('up', 'own', 'tor', 'nol', 'size')
out = {}
for nm, (Q, km) in SITES.items():
    r = {}
    r['модель как есть'] = run(Q, km)
    for f in ALL: r[f'+{f}'] = run(Q, km, fix=(f,))
    r['все поправки'] = run(Q, km, fix=ALL)
    # плюсы, которых в модели нет
    r['все + газ по полу NGFCP $0.25'] = run(Q, km, ov=dict(gas_in=0.25), fix=ALL)
    r['все + $0.25 + эскроу (coll 1, DSO 45, fx 0)'] = run(Q, km, ov=dict(gas_in=0.25, coll=1.0, dso=45, fxloss=0.0), fix=ALL)
    # доставка на забор ТЭС: покупатель экономит транспорт NGC; размер тарифа НЕ НАЙДЕН — премия как чувствительность
    for prem in (0.30, 0.80):
        r[f'все + $0.25 + эскроу + премия за доставку {prem}'] = run(Q, km, ov=dict(gas_in=0.25, coll=1.0, dso=45, fxloss=0.0,
                                                                                transport=-prem), fix=ALL)
    # капитал по внешнему ориентиру Seplat $0.22-0.36 млн на MMscf/д [В] x2.6 (правило аудита) вместо GPSA x $1250/л.с. x1.5
    hpq = g2p.hp_per_mmscfd(2.0, 50); usd_hp_seplat = 0.36e6/(hpq*1.5)/2.6*2.6   # $/л.с. так, чтобы pkg_compr = 0.36 млн/MMscf/д
    r['все + капитал компрессии по Seplat 0.36 млн/MMscf/д (x1 к ориентиру)'] = run(Q, km, ov=dict(usd_hp=0.36e6/(hpq*1.5), cm=1.0), fix=ALL)
    r['максимум: всё выше + Seplat-капитал + премия 0.80'] = run(Q, km, ov=dict(usd_hp=0.36e6/(hpq*1.5), cm=1.0, gas_in=0.25, coll=1.0,
                                                                         dso=45, fxloss=0.0, transport=-0.80), fix=ALL)
    r['порог цены NPV20=0, все поправки, $/MMBtu'] = solve(lambda x: run(Q, km, ov=dict(price=x), fix=ALL)['npv20'], 0.5, 40)
    r['порог цены NPV20=0, все + $0.25 + эскроу'] = solve(lambda x: run(Q, km, ov=dict(price=x, gas_in=0.25, coll=1.0, dso=45, fxloss=0.0), fix=ALL)['npv20'], 0.5, 40)
    out[nm] = r
# оператор: сверка с «выявленным предпочтением» — Ughelli East горит 6.3-9.1 MMscf/д с 2020 г. [О: VIIRS]
op = {}
for pen in (0.0, 0.5, 1.0, 1.75):
    op[f'pen{pen}'] = g2p.run(6.7, 3.7, 'operator', ov=dict(gas_in=-pen/(1-TAX)), contract=True)
    op[f'pen{pen}'] = {k: op[f'pen{pen}'][k] for k in ('npv15', 'npv20', 'irr')}
op['порог штрафа для оператора NPV15=0'] = solve(lambda x: g2p.run(6.7, 3.7, 'operator', ov=dict(gas_in=-x/(1-TAX)), contract=True)['npv15'], 0, 5)
out['оператор Ughelli East'] = op
# разложение маржи на Mscf входа, база
p = g2p.SC['base']
m = 2.18*1.05*(1-0.065)*0.85*(1-0.12*120/365) - 0.78
out['маржа база $/Mscf входа'] = round(m, 3)
out['маржа: $0.25 + эскроу'] = round(2.18*1.05*(1-0.065) - 0.25, 3)
out['маржа: $0.25 + эскроу + премия 0.80'] = round(2.98*1.05*(1-0.065) - 0.25, 3)
json.dump(out, open(os.path.join(HERE, 'g2p_skeptic_out.json'), 'w'), indent=1, ensure_ascii=False)
for k, v in out.items():
    print('==', k)
    if isinstance(v, dict):
        for kk, vv in v.items(): print('  ', kk, vv)
    else: print('  ', v)
