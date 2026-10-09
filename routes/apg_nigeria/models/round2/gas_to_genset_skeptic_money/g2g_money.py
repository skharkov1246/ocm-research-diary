# Скептик раунда 2, оптика «второй год и деньги» к варианту «газ генераторам покупателя» (g2g.py).
# Модель g2g не меняем: импортируем design/константы и повторяем run() с крючками по годам
# (сбор, цена, доступность газа, спрос, штраф за недопоставку, простой «по времени»).
# Метки: [П] первичный, [В] вторичный, [О] расчёт, [Д] допущение с диапазоном.
import importlib.util, os, json
HERE = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location('g', os.path.join(HERE, '..', 'gas_to_genset', 'g2g.py'))
g = importlib.util.module_from_spec(spec); spec.loader.exec_module(g)

U8 = (6.3, g.D_ENG, 1, 6.5*1.3)          # Ughelli East -> Ughelli, 8.45 км трассы, якорь на двигателях
K26 = (5, g.D_ENG, 1, 26)                # 5 MMscf/д, 26 км, якорь
BIG = (15, 6.0, 8, 26)                   # большой кластер 6 MMscf/д
LEAN = dict(price=8.0, staff=0.25, sec=0.20, surv=0.005, comm=0.08, ovh=0.10, scr=1.7, ng=0.6, cm=2.4)

def run2(Q, dem, nb, km, sc='base', ov=None, edti=False, coll_t=None, pmult_t=None, avail_t=None, dem_t=None,
         time_outage=False, short_pen=0.0, buyer_top=0.0, years=g.YEARS):
    """Копия g.run с крючками. time_outage=True: простой добычи (1-up) режет продажи по времени
    (когда месторождение стоит, газа нет вовсе), а не только средний объём факела.
    short_pen: штраф покупателю за недопоставленный DCQ, доля цены [Д]. buyer_top: take-or-pay покупателя."""
    p = dict(g.SC[sc]); p.update(ov or {})
    d = g.design(Q, dem, nb, km, p); C = d['total']
    sc_mult = 1.0 if d['Qdes'] <= 2 else (1.2 if d['Qdes'] <= 7 else 1.5)
    fixed_bu = ((p['staff']+p['ovh'])*sc_mult + p['sec'] + p['surv']*km + p['comm'] + d['fee_yr']
                + d['eq']*p['tor_eq'] + d['pipe']*p['tor_pipe'] + C*p['ins'])
    fixed_scr = 0.20*sc_mult + 0.10 + d['screen']*(0.02+0.005)
    fixed = max(fixed_bu, fixed_scr*p['scr'])
    b = p['build']; cf = [-C/b]*b
    fxh = 1 - p['dev_fx']*(p['lag_q']/2 + p['dso']/365)
    prev = 0.0; rows = []
    for t in range(1, years+1):
        T = t + b
        up = p['up']
        avail = Q*p['lf']*(1-p['decl'])**T*(1.0 if time_outage else up)*(avail_t(t) if avail_t else 1.0)
        demand = dem*(dem_t(t) if dem_t else 1.0)
        ramp = p['ramp'] if t == 1 else 1.0
        sales = min(avail*(1-p['shrink']), demand)*ramp*(up if time_outage else 1.0)
        short = max(0.0, demand*ramp - sales)
        take = sales/(1-p['shrink'])
        paid = max(take, p['top']*min(avail, demand/(1-p['shrink'])))
        price = p['price']*(1+p['pdrift'])**(t-1)*fxh*(pmult_t(t) if pmult_t else 1.0)
        billed = max(sales, buyer_top*dem*ramp*(up if time_outage else 1.0))   # ToP покупателя на законтрактованный DCQ
        rev_inv = billed*1000*365*g.HHV*price/1e6
        coll = p['coll']*(coll_t(t) if coll_t else 1.0)
        rev = rev_inv*coll
        pen = short*1000*365*g.HHV*price*short_pen/1e6
        gas = paid*1000*365*(p['gas_in']*1.02**(T-1))/1e6 + sales*1000*365*g.LIC_MSCF/1e6 + g.MDGIF*rev_inv
        fx_t = fixed*1.03**(T-1)
        e = rev - gas - fx_t - pen
        dep = C/10
        tx = max(0.0, (e-dep)*g.TAX)
        ed = (min(tx, 0.05*C) if t <= 5 else 0.0) if edti else 0.0
        wc = p['dso']/365*(rev-prev); prev = rev
        c = e - tx + ed - wc + (p['dso']/365*rev if t == years else 0)
        cf.append(c)
        rows.append(dict(t=t, avail=round(avail, 2), sales=round(sales, 2), rev=round(rev, 2), ebitda=round(e, 2),
                         cfads=round(e - tx, 2)))
    return dict(cap=round(C, 2), fixed=round(fixed, 2), e2=rows[1]['ebitda'], rows=rows, cf=cf,
                npv15=round(g.npv(.15, cf), 2), npv20=round(g.npv(.20, cf), 2),
                irr=None if g.irr(cf) is None else round(g.irr(cf)*100, 1))

def be(args, key='npv20', **kw):
    ov = kw.pop('ov', {})
    return g.solve(lambda x: run2(*args, ov=dict(ov, price=x), **kw)[key], 0.5, 120)

def brief(r): return dict(cap=r['cap'], e2=r['e2'], npv15=r['npv15'], npv20=r['npv20'], irr=r['irr'])

out = {}
# 0. Воспроизведение (run2 без крючков обязан совпасть с g.run)
for nm, a in (('U8', U8), ('K26', K26), ('BIG', BIG)):
    r0 = g.run(*a, 'base'); r = run2(*a)
    out[f'repro_{nm}'] = dict(g=dict(npv15=r0['npv15'], npv20=r0['npv20'], e2=r0['y2']['ebitda']), run2=brief(r),
                              be15=be(a, 'npv15'), be20=be(a))

# 1. Простой добычи по времени: up=0.92 в модели режет только средний объём факела; при факеле >> спроса
#    (U8: 4.5 против 1.77 MMscf/д на входе) простой никогда не связывает -> продажи завышены на 8% [О].
for nm, a in (('U8', U8), ('K26', K26), ('BIG', BIG)):
    r = run2(*a, time_outage=True)
    out[f'outage_{nm}'] = dict(brief(r), be15=be(a, 'npv15', time_outage=True), be20=be(a, time_outage=True))
#    + штраф покупателю за недопоставку 50% цены на недопоставленный объём [Д 25-100%]: покупатель в простой
#    сидит на дизеле $36.7-41.6/MMBtu [О g2g], поэтому без штрафа он не подпишет 10-летний GSA.
for nm, a in (('U8', U8), ('K26', K26)):
    r = run2(*a, time_outage=True, short_pen=0.5)
    out[f'outage_pen50_{nm}'] = dict(brief(r), be20=be(a, time_outage=True, short_pen=0.5))

# 2. Валютный шок второго года. Факт: ₦460 -> ~₦1500/$ в 2023-24 [В, nigeria_facts: долг ЦБ ~$7 млрд, 09.2023];
#    в 2026 ₦1446 -> ₦1344 [В]. Модель берёт ровную девальвацию 12%/г и сбор 95% все 10 лет.
#    Сценарий S1: шок в год 2 (х1.4 ₦/$ за год) -> покупатель, чья выручка в найре, платит хуже:
#    сбор 0.65 в годы 2-3, затем 0.90 [Д]. S2: + пересмотр цены вниз на 25% с года 3 (покупатель грозит уйти на CNG
#    или дизель) [Д 15-35%]. S3: потеря 1 года выручки (неплатёж, арбитраж) в год 3 [Д].
S1 = dict(coll_t=lambda t: 0.65/0.95 if t in (2, 3) else (0.90/0.95 if t > 3 else 1.0))
S2 = dict(S1, pmult_t=lambda t: 0.75 if t >= 3 else 1.0)
S3 = dict(coll_t=lambda t: 0.0 if t == 3 else 1.0)
for nm, a in (('U8', U8), ('K26', K26), ('BIG', BIG)):
    for sn, s in (('S1_shock_coll', S1), ('S2_shock_coll_reprice25', S2), ('S3_lost_year3', S3)):
        r = run2(*a, time_outage=True, **s)
        out[f'fx_{sn}_{nm}'] = dict(brief(r), be20=be(a, time_outage=True, **s))

# 3. Оператор забирает газ (NEPL/Seplat используют сами, отзыв/пересмотр разрешения NGFCP, DGDO): факела нет с года k.
#    GSA NGFCP: «оператор не отвечает за недопоставку по объёму и качеству» [П, презентация NUPRC в nigeria_facts].
for k in (3, 5):
    for nm, a in (('U8', U8), ('K26', K26)):
        r = run2(*a, time_outage=True, avail_t=lambda t, k=k: 0.0 if t >= k else 1.0)
        out[f'operator_takes_from_y{k}_{nm}'] = brief(r)
# ожидаемое значение при вероятности 25% / 40% отъёма с года 4 [Д]
for nm, a in (('U8', U8), ('K26', K26)):
    r_ok = run2(*a, time_outage=True); r_cut = run2(*a, time_outage=True, avail_t=lambda t: 0.0 if t >= 4 else 1.0)
    for pr in (0.25, 0.40):
        out[f'operator_cut_y4_p{pr}_{nm}'] = dict(npv15=round((1-pr)*r_ok['npv15']+pr*r_cut['npv15'], 2),
                                                  npv20=round((1-pr)*r_ok['npv20']+pr*r_cut['npv20'], 2))

# 4. Спрос покупателя: без его take-or-pay. Загрузка промышленности Нигерии ~55-60% [Д]; сеть/солнце/сокращение
#    выпуска: отбор -25% с года 3 [Д]. И зеркально: ToP покупателя 70% как защита.
D25 = dict(dem_t=lambda t: 0.75 if t >= 3 else 1.0)
for nm, a in (('U8', U8), ('K26', K26)):
    r = run2(*a, time_outage=True, **D25); out[f'buyer_minus25_{nm}'] = dict(brief(r), be20=be(a, time_outage=True, **D25))
    r = run2(*a, time_outage=True, buyer_top=0.7, **D25); out[f'buyer_minus25_top70_{nm}'] = brief(r)

# 5. Сводный «второй год и деньги»: простой по времени + штраф 50% + S1 (шок и сбор) + отъём с вероятностью 25%
def combined(a, sc='base', ov=None, edti=False):
    kw = dict(time_outage=True, short_pen=0.5, **S1)
    ok = lambda x: run2(*a, sc, ov=dict(ov or {}, price=x), edti=edti, **kw)
    cut = lambda x: run2(*a, sc, ov=dict(ov or {}, price=x), edti=edti, avail_t=lambda t: 0.0 if t >= 4 else 1.0, **kw)
    f = lambda x, key: 0.75*ok(x)[key] + 0.25*cut(x)[key]
    p0 = dict(g.SC[sc]); p0.update(ov or {})
    return dict(npv15_at_own_price=round(f(p0['price'], 'npv15'), 2), npv20_at_own_price=round(f(p0['price'], 'npv20'), 2),
                be15=g.solve(lambda x: f(x, 'npv15'), 0.5, 120), be20=g.solve(lambda x: f(x, 'npv20'), 0.5, 120))
for nm, a in (('U8', U8), ('K26', K26), ('BIG', BIG)):
    out[f'combined_base_{nm}'] = combined(a)
    out[f'combined_opt_{nm}'] = combined(a, 'opt')
out['combined_lean_EDTI_U8'] = combined(U8, ov=LEAN, edti=True)
# сводный стресс на потолке рынка (CNG + приём у покупателя $8.5-12.25; СПГ $5.2-6.4)
for nm, a in (('U8', U8), ('K26', K26), ('BIG', BIG)):
    for pr in (6.4, 9.0, 12.25):
        out[f'combined_base_{nm}_p{pr}'] = {k: v for k, v in combined(a, ov=dict(price=pr)).items() if k.startswith('npv')}
        out[f'combined_opt_{nm}_p{pr}'] = {k: v for k, v in combined(a, 'opt', ov=dict(price=pr)).items() if k.startswith('npv')}

# 6. Банкуемость: долговая ёмкость по CFADS лет 2-8 при DSCR 1.35, долларовая ставка 12-14% [Д], срок 7 лет.
def debt_capacity(a, price, rate=0.13, dscr=1.35, **kw):
    r = run2(*a, ov=dict(price=price), time_outage=True, **kw)
    cfads = [max(0.0, x['cfads']) for x in r['rows'][1:8]]
    return round(sum(c/dscr/(1+rate)**(i+1) for i, c in enumerate(cfads)), 2), r['cap']
for nm, a in (('U8', U8), ('K26', K26), ('BIG', BIG)):
    for pr in (6.0, 8.0, 10.75):
        dc, cap = debt_capacity(a, pr)
        out[f'debt_{nm}_p{pr}'] = dict(debt_musd=dc, capex=cap, gearing=round(dc/cap, 2))

# 7. Паритет покупателя с нашей энергией при ставке покупателя: модель берёт CRF 15% (долларовые деньги).
#    Покупатель с выручкой в найре: заём в найре ~25-35% [Д] или $ 15-20% [Д].
for r_b in (0.15, 0.20, 0.25, 0.30):
    for T in (350, 383):
        l0, hr = g.buyer(0.0, crf_r=r_b)
        out[f'buyer_parity_{T}_crf{int(r_b*100)}'] = round((T/g.FX - l0)/hr, 2)
# Потолок рынка с затратами покупателя на приём «виртуальной трубы» (хранилище, регазификация/редуцирование):
#    +$0.5-1.5/MMBtu [Д] к доставленной цене. СПГ Greenville ₦230/scm -> $4.7-4.9 [В]; CNG ₦380-520/scm -> $8-10.75 [В].
out['ceiling'] = dict(lng=[round(4.7+0.5, 2), round(4.9+1.5, 2)], cng=[round(8.0+0.5, 2), round(10.75+1.5, 2)])
for nm, a in (('U8', U8), ('K26', K26)):
    for pr in (6.4, 9.0, 12.25):
        r = run2(*a, ov=dict(price=pr), time_outage=True)
        out[f'at_ceiling_{nm}_p{pr}'] = brief(r)

# 8. Оператор делает сам: тот же объект без экспата (персонал и охрана площадки уже есть), газ «бесплатно»,
#    плюс избегаемый штраф за сжигание $2.0/Mscf (≥10 тыс. барр./сут) [П, Flare Regs 2018] -- gas_in < 0.
#    Сравниваем его безубыточную цену с нашей: кто кого перебьёт на переговорах с тем же покупателем.
OPS = dict(staff=0.10, ovh=0.05, sec=0.10, surv=0.006, comm=0.05, scr=0.0, top=0.0, dev=0.4)
for nm, a in (('U8', U8), ('K26', K26)):
    for fine in (0.0, 0.5, 2.0):
        out[f'operator_self_{nm}_fine{fine}'] = be(a, ov=dict(OPS, gas_in=-fine), time_outage=True)
    out[f'us_{nm}_be20_time_outage'] = be(a, time_outage=True)

# 9. Валюта у нас: если цена фиксирована в найре (покупатель настоял), девальвация 12% и шок S1
for nm, a in (('U8', U8), ('K26', K26)):
    r = g.run(*a, 'base', naira_price=True); out[f'naira_fixed_{nm}'] = dict(npv15=r['npv15'], npv20=r['npv20'])

# 10. Системный капитал на 61.9 ГВт·ч/г (наш + двигатели покупателя $1 500-1 667/кВт [В]) против «энергии якорю» $29.15 млн [О]
out['system_capital'] = dict(U8=[round(out['repro_U8']['run2']['cap']+11.5*x/1e3, 1) for x in (1500, 1667)],
                             K26=[round(out['repro_K26']['run2']['cap']+11.5*x/1e3, 1) for x in (1500, 1667)],
                             power=29.15)

json.dump(out, open(os.path.join(HERE, 'g2g_money_out.json'), 'w'), indent=1, ensure_ascii=False)
txt = '\n'.join(f'{k}: {v}' for k, v in out.items()); print(txt)
open(os.path.join(HERE, 'g2g_money_stdout.txt'), 'w').write(txt)
