# Скептик BOOT, оптика «второй год и деньги»: банкуемость, валюта, что сломается, что сделает оператор/конкурент.
# Своя (независимая) реализация денежного потока; сверка с ../compression_boot/boot.py на базовых точках.
# Метки: [П] первичный, [В] вторичный, [О] расчёт, [Д] допущение с диапазоном.
# Запуск: python3 sk_money.py -> sk_money_out.json и печать.
import sys, os, json
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'compression_boot'))
import boot as B          # берём только capex() и SC/CONF (физика и аудит-поправки), денежный поток считаем сами
out = {}
def pr(k, v): out[k] = v; print(k, v)
def npv(r, cf): return sum(c/(1+r)**t for t, c in enumerate(cf))
def solve(fn, lo, hi, it=60):
    flo = fn(lo); fhi = fn(hi)
    if (flo > 0) == (fhi > 0): return None
    for _ in range(it):
        m = (lo+hi)/2; fm = fn(m)
        if (fm > 0) == (flo > 0): lo, flo = m, fm
        else: hi = m
    return round((lo+hi)/2, 3)

def cash(Q, sc='base', fee=1.375, conf='boost', ov=None, usd_share=1.0, deval=0.0, cost_ngn=0.0,
         fee_after=None, reneg_year=None, stop_year=None, term_pay=0.0, dso=None, coll=None, years=10):
    """Годовой денежный поток подрядчика, $ млн. usd_share — доля тарифа в долларах, остальное в найре,
    фиксированной при подписании; deval — девальвация найры в год; cost_ngn — доля постоянного опекса в найре,
    индексируемой местной инфляцией (реальная девальвация = 0 -> в $ неизменна; здесь считаем, что найровый опекс
    растёт с инфляцией ~= девальвации, т.е. в долларах не дешевеет — консервативно [Д]).
    reneg_year/fee_after — с этого года оператор пересматривает тариф до fee_after (hold-up).
    stop_year — оператор прекращает платить с этого года (смена владельца, свой компрессор), term_pay — выплата
    при расторжении, доля остаточной балансовой стоимости."""
    p = dict(B.SC[sc]); p.update(ov or {})
    if dso is not None: p['dso'] = dso
    if coll is not None: p['coll'] = coll
    cx = B.capex(Q, conf, p); C = cx['total']
    ss = B.staff_scale(Q)
    fixed = max((p['staff']+p['ovh'])*ss + p['sec'] + C*(p['tor']+p['ins']),
                (0.20*ss + 0.10 + cx['screen']*0.025)*p['mult'])
    b = p['build']; cf = [-C/b]*b; prev = 0.0; rows = []
    for t in range(1, years+1):
        G = Q*p['lf']*(1-p['decl'])**(t-1)*p['up']*(p['ramp'] if t == 1 else 1.0)
        f = fee*(1+p['fee_drift'])**(t-1)
        if reneg_year and t >= reneg_year: f = min(f, fee_after)
        f_usd = f*usd_share + f*(1-usd_share)/(1+deval)**(t-1)
        rev = G*365*f_usd/1000*p['coll']*(1-B.NCDMB)
        if stop_year and t >= stop_year: rev = 0.0
        opex = fixed*1.03**(t-1) + G*365*p['var']/1000
        if stop_year and t >= stop_year: opex = 0.0          # станцию останавливаем/консервируем
        e = rev - opex
        tx = max(0.0, (e - C/10)*B.TAX)
        wc = p['dso']/365*(rev-prev); prev = rev
        c = e - tx - wc
        if stop_year and t == stop_year:
            c += term_pay*C*(1-(t-1)/10)                      # компенсация доли остаточной стоимости
        if t == years: c += p['dso']/365*rev
        cf.append(c); rows.append(dict(t=t, G=round(G, 2), rev=round(rev, 2), opex=round(opex, 2), ebitda=round(e, 2),
                                       tax=round(tx, 2), wc=round(wc, 2), cfads=round(e-tx, 2)))
        if stop_year and t >= stop_year: break
    return dict(C=round(C, 2), fixed=round(fixed, 2), cf=cf, rows=rows, npv15=round(npv(.15, cf), 2),
                npv20=round(npv(.20, cf), 2))

# 0. Сверка с boot.py [О]
for Q, fee in ((5, 1.375), (15, 1.375), (15, 1.75)):
    a = cash(Q, 'base', fee); b_ = B.run(Q, 'boost', 'base', fee)
    pr(f'0_check_Q{Q}_fee{fee}', dict(mine=(a['C'], a['rows'][0]['ebitda'], a['rows'][1]['ebitda'], a['npv15'], a['npv20']),
                                      boot=(b_['capex']['total'], b_['y1']['ebitda'], b_['y2']['ebitda'], b_['npv15'], b_['npv20'])))
thr = {}
for Q in (5, 15):
    for sc in ('opt', 'base'):
        thr[(Q, sc)] = solve(lambda x: cash(Q, sc, x)['npv20'], 0.01, 40)
pr('0_thr_npv20', {f'Q{k[0]}_{k[1]}': v for k, v in thr.items()})

# 1. Банкуемость: долг 60% капитала, 12% USD [Д: 10-14%], 7 лет равными аннуитетами после 1 г. стройки.
#    DSCR года 2 и минимальный DSCR при пороговом тарифе (NPV20=0) и при «рыночном» $1.75. Банк хочет >=1.3-1.4 [Д].
def dscr(Q, sc, fee, gear=0.6, r=0.12, n=7, **kw):
    a = cash(Q, sc, fee, **kw); D = a['C']*gear
    ds = D*r/(1-(1+r)**-n)
    cfads = [row['cfads'] - row['wc'] for row in a['rows']]
    return dict(debt=round(D, 2), debt_service=round(ds, 2), dscr_y1=round(cfads[0]/ds, 2), dscr_y2=round(cfads[1]/ds, 2),
                dscr_min_y1_7=round(min(cfads[:n])/ds, 2), dscr_y7=round(cfads[n-1]/ds, 2))
for Q in (5, 15):
    for sc in ('opt', 'base'):
        pr(f'1_DSCR_Q{Q}_{sc}_fee1.75', dscr(Q, sc, 1.75))
        pr(f'1_DSCR_Q{Q}_{sc}_at_thr{thr[(Q, sc)]}', dscr(Q, sc, thr[(Q, sc)]))
# тариф, при котором минимальный DSCR лет 1-7 = 1.3 (кредитор, а не NPV, задаёт цену)
for Q in (5, 15):
    for sc in ('opt', 'base'):
        pr(f'1_fee_for_minDSCR1.3_Q{Q}_{sc}', solve(lambda x: dscr(Q, sc, x)['dscr_min_y1_7']-1.3, 0.01, 60))
        pr(f'1_fee_for_DSCR1.3_y2on_Q{Q}_{sc}', solve(lambda x: min(dscr(Q, sc, x)['dscr_y2'], dscr(Q, sc, x)['dscr_y7'])-1.3, 0.01, 60))

# 2. Валюта. boot.py: тариф 100% $ («главный плюс»). Для подрядчика с регистрацией NCDMB и Nigerian Content Plan
#    долларовая доля контракта НЕ НАЙДЕНА; [Д] 50-100% в $. Найровая часть фиксирована, девальвация 8/16%/г.
for Q in (15,):
    for us in (1.0, 0.7, 0.5):
        for dv in (0.08, 0.16):
            a = cash(Q, 'base', thr[(Q, 'base')], usd_share=us, deval=dv)
            pr(f'2_fx_Q{Q}_base_fee{thr[(Q, "base")]}_usd{us}_dev{dv}', dict(npv20=a['npv20'], npv15=a['npv15']))
        pr(f'2_fee_npv20_Q{Q}_base_usd{us}_dev0.16', solve(lambda x: cash(Q, 'base', x, usd_share=us, deval=0.16)['npv20'], 0.01, 60))
        pr(f'2_fee_npv20_Q{Q}_opt_usd{us}_dev0.16', solve(lambda x: cash(Q, 'opt', x, usd_share=us, deval=0.16)['npv20'], 0.01, 60))

# 3. Hold-up: станция специфична и не перевозится дёшево (монтаж ×2.6 — невозвратная часть). После года N оператор
#    требует снизить тариф до своей себестоимости «сделать сам» (DIY, по симметричному аудиту $2.19 для 15 база,
#    sk_boot.py) или до аренды с обвязкой ($0.44-0.74 + $0.56). Какой стартовый тариф тогда нужен для NPV20=0?
for Q, sc, diy in ((15, 'base', 2.19), (15, 'opt', 1.11), (5, 'base', 3.12)):
    for ry in (3, 4, 6):
        pr(f'3_holdup_Q{Q}_{sc}_reneg_y{ry}_to{diy}', solve(lambda x: cash(Q, sc, x, reneg_year=ry, fee_after=diy)['npv20'], 0.01, 80))
    pr(f'3_holdup_Q{Q}_{sc}_reneg_y3_to1.30(rent+BoP)', solve(lambda x: cash(Q, sc, x, reneg_year=3, fee_after=1.30)['npv20'], 0.01, 80))

# 4. Остановка платежей (смена владельца, отзыв лицензии, свой компрессор оператора): с года 3/5, без компенсации и
#    с компенсацией 100% остаточной стоимости. NPV при пороговом тарифе.
for Q, sc in ((15, 'base'), (15, 'opt')):
    f0 = thr[(Q, sc)]
    for sy in (3, 5):
        for tp in (0.0, 1.0):
            a = cash(Q, sc, f0, stop_year=sy, term_pay=tp)
            pr(f'4_stop_Q{Q}_{sc}_fee{f0}_from_y{sy}_comp{tp}', dict(npv15=a['npv15'], npv20=a['npv20']))
# вероятностно: риск смены владельца/остановки p в год (Renaissance/SPDC, Oando/NAOC, Seplat/MPNU — 3 сделки за 2024-25) [Д 5-10%/г]
def exp_npv(Q, sc, fee, h, tp=0.0, r=.20):
    tot = 0; surv = 1.0
    for sy in range(2, 11):
        pstop = surv*h; tot += pstop*npv(r, cash(Q, sc, fee, stop_year=sy, term_pay=tp)['cf']); surv -= pstop
    tot += surv*npv(r, cash(Q, sc, fee)['cf'])
    return tot
for h in (0.05, 0.10):
    for tp in (0.0, 1.0):
        pr(f'4_fee_npv20_Q15_base_hazard{h}_comp{tp}', solve(lambda x: exp_npv(15, 'base', x, h, tp), 0.01, 60))
        pr(f'4_fee_npv20_Q15_opt_hazard{h}_comp{tp}', solve(lambda x: exp_npv(15, 'opt', x, h, tp), 0.01, 60))

# 5. Дебиторка NEPL: история cash call ~$8.5 млрд к 2016 г. [В]. DSO 180/365, сбор 0.85/0.8 — порог.
for d, c in ((90, 0.92), (180, 0.85), (365, 0.80)):
    pr(f'5_fee_npv20_Q15_base_dso{d}_coll{c}', solve(lambda x: cash(15, 'base', x, dso=d, coll=c)['npv20'], 0.01, 60))
    pr(f'5_fee_npv20_Q15_opt_dso{d}_coll{c}', solve(lambda x: cash(15, 'opt', x, dso=d, coll=c)['npv20'], 0.01, 60))

# 6. Конкурент: Enerflex-подобный (капитал ×1.6, без экспата) — порог $2.74 (boot.py extra). Если он ставит
#    свою цену NPV20=0 — что получаем мы по нашей базе при этой цене. И сколько маржи у оператора остаётся в
#    сравнении со штрафом+газом: WTP = P/(1-t) + нетбэк (0-1.0) [Д].
for f in (2.74, 2.55, 2.19):
    a = cash(15, 'base', f); o = cash(15, 'opt', f)
    pr(f'6_us_at_competitor_fee{f}', dict(base_npv20=a['npv20'], base_npv15=a['npv15'], opt_npv20=o['npv20'], opt_npv15=o['npv15'],
                                         base_E1=a['rows'][0]['ebitda'], base_E2=a['rows'][1]['ebitda']))
pr('6_operator_WTP_t0.34', {f'P{P}_N{N}': round(P/(1-0.34)+N, 2) for P in (1.6, 1.9) for N in (0.0, 0.5, 1.0)})

# 7. Второй год по деньгам (база 15, тариф по порогу и по $1.75): EBITDA, налог, оборотный капитал, CFADS
for f in (1.75, thr[(15, 'base')]):
    a = cash(15, 'base', f)
    pr(f'7_years1_3_Q15_base_fee{f}', a['rows'][:3])

json.dump(out, open(os.path.join(HERE, 'sk_money_out.json'), 'w'), indent=1, ensure_ascii=False, default=str)
