# Скептик к power_plus_pipe/bundle.py («энергия якорю + газ в трубу»). Метки: [П] [В] [О] [Д].
# Проверяем: (1) держится ли «не делать», если снять перекосы в минус (двойной учёт простоя VIIRS в трубе,
# ТОиР 4.25% на линейную часть трубы, газ по полу NGFCP, эскроу у ТЭС, транспортная премия за доставку до забора ТЭС);
# (2) что модель пропускает в плюс связке: почасовой профиль нагрузки якоря и остановы станции против
# мощности компрессии 1.25*Qd; усадка на точке росы по УВ; (3) альтернативу «врезка в действующий газовый завод
# Ughelli East (NPDC/NEPL, ~90 MMscf/д, питает Delta/Transcorp Ughelli) [В: Wikipedia, List of natural gas processing
# plants in Nigeria]» вместо собственного отвода 4.8 км; (4) сравнение «оператор в 9-14 раз дешевле».
# Запуск: python3 sk_bundle.py -> sk_bundle_out.json рядом.
import os, sys, json
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'power_plus_pipe'))
import bundle as B
G = B.G
OUT = {}
SITES = {k: B.SITES[k] for k in ('Ughelli East 2025 (6.7) -> Transcorp Ughelli', 'Oredo min23-25 (8.5) -> Sapele PS',
                                 'Oredo 2025 (15.5) -> Sapele PS', 'типовой 5 MMscf/д, отвод 4 км', 'типовой 15 MMscf/д, отвод 8 км')}
T = 411; E = True
orig_fixed = B.pipe_fixed
def fixed_tor_split(Q_site, cx, pp, shared):
    # ТОиР 4.25% только на смонтированный пакет (pkg x cm), линейная часть 1%/г, страховка 1% от капитала [Д, как g2p_skeptic]
    main = cx['pkg']*pp['cm']*pp['tor'] + cx['pipe']*0.01 + cx['C']*0.01
    if shared:
        bu = B.SHARE['staff'] + B.SHARE['ovh'] + B.SHARE['sec'] + B.SHARE['comm'] + 0.01*cx['L'] + main
        scr = 0.08 + cx['screen']*0.025
    else:
        ss = G.staff_scale(Q_site)
        bu = (pp['staff']+pp['ovh'])*ss + pp['sec'] + 0.01*cx['L'] + pp['comm'] + main
        scr = 0.20*ss + 0.10 + cx['screen']*0.025
    return max(bu, scr*pp['opx_screen_mult'])

def run(Q, lk, pk, sc='base', fixes=(), **kw):
    B.UP_PROD = 1.0 if 'up' in fixes else 11/12
    B.pipe_fixed = fixed_tor_split if 'tor' in fixes else orig_fixed
    try:
        return B.model(Q, lk, pk, sc, T, E, kw.pop('mode', 'bundle'), **kw)
    finally:
        B.UP_PROD = 11/12; B.pipe_fixed = orig_fixed

ESC = dict(coll=0.98, dso=45)
STACKS = {
  'база (как в оценке)': dict(),
  'F1: без двойного простоя + ТОиР трубы 1%': dict(fixes=('up', 'tor')),
  'F2: F1 + газ пол NGFCP $0.40/Mscf + эскроу': dict(fixes=('up', 'tor'), gas_allin=0.40/1.05, pov=dict(ESC)),
  'F3: F2 + премия за доставку $0.30/MMBtu [Д]': dict(fixes=('up', 'tor'), gas_allin=0.40/1.05, pov=dict(ESC, transport=-0.30)),
  'F4: F2 + премия за доставку $0.80/MMBtu [Д]': dict(fixes=('up', 'tor'), gas_allin=0.40/1.05, pov=dict(ESC, transport=-0.80)),
  'F5: F4 + цена $2.68 (коммерч., не ТЭС) [Д]': dict(fixes=('up', 'tor'), gas_allin=0.40/1.05, price=2.68, pov=dict(ESC, transport=-0.80)),
  'F6: F5 + сильная синергия': dict(fixes=('up', 'tor'), gas_allin=0.40/1.05, price=2.68, pov=dict(ESC, transport=-0.80),
                                    shr=dict(staff=0.12, sec=0.0, dev=0.2, corridor_k=0.8, dehy_k=0.6)),
}
OUT['stacks'] = {}
for nm, (Q, lk, pk) in SITES.items():
    for sc in ('base', 'opt'):
        p0 = run(Q, lk, pk, sc, mode='power')
        for sn, kw in STACKS.items():
            kw = dict(kw); fx = kw.pop('fixes', ())
            r = run(Q, lk, pk, sc, fixes=fx, **kw)
            OUT['stacks'][f'{nm}|{sc}|{sn}'] = dict(capex=r['capex'], cap_pipe=r['capex_pipe'], y2_pipe=r['y2_pipe'],
                npv15=r['npv15'], inc15=round(r['npv15']-p0['npv15'], 2), inc20=round(r['npv20']-p0['npv20'], 2), irr=r['irr'])

# (2) почасовой профиль: якорь работает на 1.0 долю f времени и на L остальное (средняя нагрузка 0.72 [Д deal_r2]),
# станция недоступна 8% времени (готовность 0.92). Труба берёт min(cap, avail - burn_h). Сколько газа, который модель
# на годовых средних отправляет в трубу, на деле не помещается в компрессию 1.25*Qd. [О]
def clip_share(Q, lk, pk, sc='base', f=0.6, av=0.92, LF=0.72):
    r = run(Q, lk, pk, sc); row = r['rows'][1]
    cap = r['Qd_pipe']*1.25
    avail = row['avail']; burned_mean = row['burned']
    # burn при полной нагрузке, нормированный так, чтобы среднее совпало со средним модели
    L = (LF - f)/(1 - f)
    full = burned_mean/(av*LF)
    states = [(av*f, full*1.0), (av*(1-f), full*L), (1-av, 0.0)]
    take_h = sum(w*min(cap, max(0.0, avail - b)) for w, b in states)
    take_avg = min(cap, max(0.0, avail - burned_mean))
    return dict(cap=round(cap, 2), avail=avail, burned=burned_mean, L=round(L, 2),
                pipe_take_avg_model=round(take_avg, 3), pipe_take_hourly=round(take_h, 3),
                lost_share=round(1 - take_h/take_avg, 3) if take_avg > 0 else None)
OUT['hourly'] = {f'{nm}|f={f}': clip_share(*SITES[nm], f=f) for nm in SITES for f in (0.5, 0.6, 0.72)}

# (2б) усадка на точке росы по УВ (C3+ уходит из потока в трубу) +3/+8% [Д; состав газа НЕ НАЙДЕН]
OUT['hcdp'] = {}
for nm, (Q, lk, pk) in SITES.items():
    p0 = run(Q, lk, pk, 'base', mode='power')
    for add in (0.03, 0.08):
        B.SHRINK = G.SHRINK + 0.01 + add
        r = run(Q, lk, pk, 'base')
        B.SHRINK = G.SHRINK + 0.01
        OUT['hcdp'][f'{nm}|+{add}'] = dict(inc15=round(r['npv15']-p0['npv15'], 2))

# (3) врезка в действующий газовый завод Ughelli East (NPDC/NEPL, ~90 MMscf/д, питает Delta PS [В]):
# отвод 0.5 км по прямой, нагнетание 25 бар (вход завода, НЕ НАЙДЕНО), покупатель — оператор; цена сырого газа
# на входе завода НЕ НАЙДЕНА: диапазон $1.0-2.18/MMBtu [Д]. Сбор у NEPL как у ТЭС (0.85) и как эскроу (0.98).
OUT['tie_in'] = {}
Q, lk, _ = SITES['Ughelli East 2025 (6.7) -> Transcorp Ughelli']
for sc in ('base', 'opt'):
    p0 = run(Q, lk, 3.7, sc, mode='power')
    for pr in (1.0, 1.5, 2.18):
        for cn, cov in (('сбор как в базе', {}), ('эскроу', ESC)):
            r = run(Q, lk, 0.5, sc, price=pr, pov=dict(cov, p_dis=25))
            OUT['tie_in'][f'{sc}|${pr}|{cn}'] = dict(cap_pipe=r['capex_pipe'], hp=r['hp'], y2_pipe=r['y2_pipe'],
                                                    inc15=round(r['npv15']-p0['npv15'], 2), inc20=round(r['npv20']-p0['npv20'], 2))
    # и с F2-поправками
    for pr in (1.5, 2.18):
        r = run(Q, lk, 0.5, sc, fixes=('up', 'tor'), price=pr, gas_allin=0.40/1.05, pov=dict(ESC, p_dis=25))
        OUT['tie_in'][f'{sc}|${pr}|F2'] = dict(cap_pipe=r['capex_pipe'], y2_pipe=r['y2_pipe'],
                                               inc15=round(r['npv15']-p0['npv15'], 2), inc20=round(r['npv20']-p0['npv20'], 2))

# (4) «оператор в 9-14 раз дешевле»: сравнение одной и той же величины [О]
r = run(*SITES['Ughelli East 2025 (6.7) -> Transcorp Ughelli'])
pp = G.SC['base']; Qd = r['Qd_pipe']
compr_raw = r['hp']*pp['usd_hp']/1e6
OUT['operator_ratio'] = dict(Qd=Qd, hp=r['hp'], compr_pkg_raw_musd=round(compr_raw, 2),
    per_mmscfd_raw=round(compr_raw/Qd, 3), per_mmscfd_x2_6=round(compr_raw*2.6/Qd, 3),
    whole_block_per_mmscfd=round(r['capex_pipe']/Qd, 2),
    ratio_like_for_like_raw=[round(compr_raw/Qd/x, 1) for x in (0.36, 0.22)],
    ratio_like_for_like_x2_6=[round(compr_raw*2.6/Qd/x, 1) for x in (0.36, 0.22)],
    ratio_in_assessment=[round(r['capex_pipe']/Qd/x, 1) for x in (0.36, 0.22)])

# (5) рабочие часы: доля факела, реально снятая (газ станции + труба) против доступного по VIIRS (без 0.78) — справка
json.dump(OUT, open(os.path.join(HERE, 'sk_bundle_out.json'), 'w'), indent=1, ensure_ascii=False)
for sec, d in OUT.items():
    print('==', sec)
    if isinstance(d, dict):
        for k, v in d.items(): print('  ', k, v)

# (6) оптимизм + «стек A» оценки (цена $2.68 = уже opt, газ по полу NGFCP, эскроу) БЕЗ моих поправок up/tor:
# оценка утверждает, что для оптимизма нужна цена от $3.3/MMBtu (выше законного верха) или плата оператора от $0.6/Mscf.
OUT['opt_stackA'] = {}
for nm, (Q, lk, pk) in SITES.items():
    for TT in (405, 411):
        p0 = B.model(Q, lk, pk, 'opt', TT, True, 'power')
        for lab, kw in (('стек A', dict(price=2.68, gas_allin=0.40/1.05, pov=dict(ESC))),
                        ('только газ по полу', dict(gas_allin=0.40/1.05))):
            r = B.model(Q, lk, pk, 'opt', TT, True, 'bundle', **kw)
            OUT['opt_stackA'][f'{nm}|T{TT}|{lab}'] = dict(inc15=round(r['npv15']-p0['npv15'], 2), inc20=round(r['npv20']-p0['npv20'], 2))
json.dump(OUT, open(os.path.join(HERE, 'sk_bundle_out.json'), 'w'), indent=1, ensure_ascii=False)
print('== opt_stackA')
for k, v in OUT['opt_stackA'].items(): print('  ', k, v)

# (7) ОШИБКА СРАВНЕНИЯ: в bundle._model gas_allin меняет p['gas'] и для энергоблока (D.station), а приращение
# в stack.py / thresholds['gsa_allin'] / sens считается против power БЕЗ этой цены. Часть «приращения трубы» —
# это удешевление газа станции. Пересчёт: power с той же gas_allin. [О]
OUT['gas_fix'] = {}
for nm, (Q, lk, pk) in SITES.items():
    for sc in ('base', 'opt'):
        for lab, kw in (('стек A', dict(price=2.68, gas_allin=0.40/1.05, pov=dict(ESC))),
                        ('газ по полу', dict(gas_allin=0.40/1.05)),
                        ('F2 (up+tor+пол+эскроу)', dict(fixes=('up', 'tor'), gas_allin=0.40/1.05, pov=dict(ESC))),
                        ('F4 (F2 + премия доставки $0.80)', dict(fixes=('up', 'tor'), gas_allin=0.40/1.05, pov=dict(ESC, transport=-0.80))),
                        ('F6 (F4 + $2.68 + сильная синергия)', dict(fixes=('up', 'tor'), gas_allin=0.40/1.05, price=2.68,
                              pov=dict(ESC, transport=-0.80), shr=dict(staff=0.12, sec=0.0, dev=0.2, corridor_k=0.8, dehy_k=0.6)))):
            kw = dict(kw); fx = kw.pop('fixes', ())
            p_old = run(Q, lk, pk, sc, mode='power')
            p_new = run(Q, lk, pk, sc, mode='power', gas_allin=kw['gas_allin'])
            r = run(Q, lk, pk, sc, fixes=fx, **kw)
            OUT['gas_fix'][f'{nm}|{sc}|{lab}'] = dict(inc15_as_in_assessment=round(r['npv15']-p_old['npv15'], 2),
                inc15_correct=round(r['npv15']-p_new['npv15'], 2), inc20_correct=round(r['npv20']-p_new['npv20'], 2),
                power_gain_from_gas15=round(p_new['npv15']-p_old['npv15'], 2))
    # порог платы оператора поверх стека A с правильной базой (NPV20 приращения = 0), base, T=405
    kwA = dict(price=2.68, gas_allin=0.40/1.05, pov=dict(ESC))
    p_new = B.model(Q, lk, pk, 'base', 405, True, 'power', gas_allin=0.40/1.05)
    OUT['gas_fix'][f'{nm}|base|порог платы оператора поверх стека A, верно'] = B.solve(
        lambda x: B.model(Q, lk, pk, 'base', 405, True, 'bundle', fee_op=x, **kwA)['npv20'] - p_new['npv20'], 0, 15)
json.dump(OUT, open(os.path.join(HERE, 'sk_bundle_out.json'), 'w'), indent=1, ensure_ascii=False)
print('== gas_fix')
for k, v in OUT['gas_fix'].items(): print('  ', k, v)
