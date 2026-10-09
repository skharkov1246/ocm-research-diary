# Скептик связки «энергия якорю + газ в трубу» через оптику «второй год и деньги»: банкуемость, валюта,
# что сломается, что сделают оператор и конкурент. Модель ../power_plus_pipe/bundle.py используется без изменений
# логики; добавлены только флаги через патч исходника (раздельный налог при ring-fence, обрыв потока трубы).
# Метки: [П] первичный, [В] вторичный, [О] расчёт, [Д] допущение.
# Запуск: python3 sk_money.py -> sk_money_out.json рядом.
import os, sys, json, inspect, textwrap
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'power_plus_pipe'))
import bundle as B
D = B.D
npv, irr = B.npv, B.irr

# --- патч: флаги CONSOL (консолидированный налог в связке) и STRAND (доля потока трубы по годам) ---
src = textwrap.dedent(inspect.getsource(B._model))
src = src.replace("        if shared:   # консолидированный налог", "        if shared and CONSOL:   # консолидированный налог")
src = src.replace("            avail = base_flow[i]\n", "            avail = base_flow[i]*STRAND(t)\n")
assert 'CONSOL' in src and 'STRAND(t)' in src
B.CONSOL = True; B.STRAND = lambda t: 1.0
exec(compile(src, 'patched', 'exec'), B.__dict__)

SITES = {k: B.SITES[k] for k in ('Ughelli East 2025 (6.7) -> Transcorp Ughelli', 'Oredo min23-25 (8.5) -> Sapele PS',
                                 'Oredo 2025 (15.5) -> Sapele PS', 'типовой 5 MMscf/д, отвод 4 км', 'типовой 15 MMscf/д, отвод 8 км')}
T, E = 411, True
OUT = {}

def m(nm, mode, sc='base', **kw):
    Q, lk, pk = B.SITES[nm]
    return B.model(Q, lk, pk, sc, kw.pop('T', T), kw.pop('edti', E), mode, **kw)

# 0. регрессия ключевых чисел оценки
OUT['regress'] = {}
for nm in SITES:
    p, b = m(nm, 'power'), m(nm, 'bundle')
    OUT['regress'][nm] = dict(power=dict(capex=p['capex'], y2=p['y2'], npv15=p['npv15'], npv20=p['npv20'], irr=p['irr']),
                              bundle=dict(capex=b['capex'], y2=b['y2'], y2_pipe=b['y2_pipe'], npv15=b['npv15'], npv20=b['npv20'],
                                          irr=b['irr'], pb=b['payback_from_t0']))

# 1. Банкуемость. Долг как в matrix_r2: 60% / 14% / 7 лет [Д, как deal]. CFADS = денежный поток проекта после налога
#    (cf модели, без налогового щита процентов — одинаково для обоих режимов). Стройка 1 г: долг выдан в t=0, платежи с t=2.
def debt_metrics(res, delay=1, rate=0.14, ten=7, share=0.60, dscr_target=1.40):
    cf = res['cf']; C = res['capex']; cfads = cf[delay+1:delay+1+ten]
    ann = lambda Dm: Dm*rate/(1-(1+rate)**-ten)
    D60 = C*share
    min_dscr60 = min(cfads)/ann(D60)
    # размер долга по минимальному DSCR (кредитор считает по худшему году) [Д 1.30-1.50; Нигерия, USD]
    Dsz = min(min(cfads)/dscr_target*(1-(1+rate)**-ten)/rate, C*0.70)
    Dsz = max(Dsz, 0.0)
    eq = C - Dsz
    lev = [-eq] + [0.0]*delay + [c - (ann(Dsz) if k < ten else 0.0) for k, c in enumerate(cf[delay+1:])]
    r = irr(lev)
    return dict(capex=C, min_dscr_at60=round(min_dscr60, 2), debt_sized=round(Dsz, 1), equity=round(eq, 1),
                eq_npv20=round(npv(.20, lev), 2), eq_npv25=round(npv(.25, lev), 2), eq_irr=None if r is None else round(r*100, 1),
                y2_cfads=round(cf[delay+1+1], 2))
OUT['bank'] = {}
for nm in SITES:
    for sc in ('base', 'opt'):
        p, b = m(nm, 'power', sc), m(nm, 'bundle', sc)
        dl = B.SCN[sc]['delay']
        OUT['bank'][f'{nm}|{sc}'] = dict(power=debt_metrics(p, dl), bundle=debt_metrics(b, dl))

# 2. Ring-fence: кредитор энергоблока требует отдельную SPV/водопад для трубы -> налог раздельный, общие затраты остаются
OUT['ringfence'] = {}
for nm in SITES:
    for sc in ('base', 'opt'):
        p = m(nm, 'power', sc); b = m(nm, 'bundle', sc)
        B.CONSOL = False
        try: r = m(nm, 'bundle', sc)
        finally: B.CONSOL = True
        OUT['ringfence'][f'{nm}|{sc}'] = dict(bundle_consol=b['npv15'], bundle_ringfenced=r['npv15'],
                                              tax_consol_value=round(b['npv15']-r['npv15'], 2),
                                              inc15_ringfenced=round(r['npv15']-p['npv15'], 2), inc20_ringfenced=round(r['npv20']-p['npv20'], 2))

# 3. ToP «хвост» — артефакт GSA «только энергии»: при ToP на 70% DCQ или pay-as-burned у станции 0.81 MMscf/д не оплачиваются
#    впустую, «только энергия» дорожает, а синергия связки по ToP исчезает. D.run(top=...) — штатный параметр deal_r2.
OUT['top'] = {}
for nm in SITES:
    Q, lk, pk = B.SITES[nm]
    pp = B.G.SC['base']
    pr = dict(D.NEW, Q=Q, line_km=lk, loss=0.03*lk/25, avail=pp['lf'], decl=pp['decl'], cm=2.6, p_dir=T)
    tl = dict(D.TL_NEW, dep=0.12, delay=1); s = D.STRUCTS['dopgas']
    ref = D.run(pr, s, tl, edti=E)
    b = m(nm, 'bundle')
    row = dict(power_bundle_py=m(nm, 'power')['npv15'], power_deal_r2=ref['npv15'], bundle=b['npv15'])
    for top in (0.8, 0.7, 0.0):
        x = D.run(pr, s, tl, edti=E, top=top)
        row[f'power_top{top}'] = x['npv15']; row[f'inc15_top{top}'] = round(b['npv15']-x['npv15'], 2)
    OUT['top'][nm] = row

# 4. Оператор/конкурент «во второй год»: поток трубы обрывается с операционного года k (самостройка NEPL,
#    перевод газа в свой газовый завод, зачёт в DGDO, ТЭС уходит на газ NEPL). Энергоблок не трогаем.
OUT['strand'] = {}
for nm in SITES:
    p = m(nm, 'power')
    for k in (3, 4, 6):
        B.STRAND = (lambda kk: (lambda t: 1.0 if t < kk else 0.0))(k)
        try: r = m(nm, 'bundle')
        finally: B.STRAND = lambda t: 1.0
        OUT['strand'][f'{nm}|обрыв с года {k}'] = dict(npv15=r['npv15'], inc15=round(r['npv15']-p['npv15'], 2), inc20=round(r['npv20']-p['npv20'], 2))

# 5. «Оптимизм» с законной ценой для энергетики: s.167(5) — ТЭС платит DBP $2.18, а не коммерческие $2.68 [П]
OUT['opt_legal'] = {}
for nm in SITES:
    p = m(nm, 'power', 'opt'); b = m(nm, 'bundle', 'opt'); bl = m(nm, 'bundle', 'opt', price=2.18)
    OUT['opt_legal'][nm] = dict(power=p['npv15'], bundle_2_68=b['npv15'], bundle_2_18=bl['npv15'],
                                inc15_2_68=round(b['npv15']-p['npv15'], 2), inc15_2_18=round(bl['npv15']-p['npv15'], 2),
                                inc20_2_18=round(bl['npv20']-p['npv20'], 2))

# 6. Валюта: GSA ТЭС номинирован в $, платится в найре по курсу ЦБ [П: условия NGFCP/GSA; nigeria_facts]; дебиторка gencos
#    по Transcorp I кв.2025: ~305 дней, сбор ~48% [О из В, power_plus_pipe_ng]. Девальвация 8/12/16% [Д] на этой дебиторке.
OUT['fx'] = {}
for nm in ('Ughelli East 2025 (6.7) -> Transcorp Ughelli', 'типовой 15 MMscf/д, отвод 8 км'):
    p = m(nm, 'power')
    for dep in (0.08, 0.12, 0.16):
        for coll, dso in ((0.85, 120), (0.85, 305), (0.50, 305)):
            saved = dict(B.SCN['base']); B.SCN['base']['dep'] = dep
            try:
                p2 = m(nm, 'power'); r = m(nm, 'bundle', pov=dict(coll=coll, dso=dso))
            finally: B.SCN['base'].clear(); B.SCN['base'].update(saved)
            OUT['fx'][f'{nm}|dep{dep}|coll{coll}|dso{dso}'] = dict(power=p2['npv15'], bundle=r['npv15'], inc15=round(r['npv15']-p2['npv15'], 2),
                                                                   y2_pipe=r['y2_pipe'])

# 7. Второй год: EBITDA трубы и её доля в обслуживании долга (Ughelli, база)
nm = 'Ughelli East 2025 (6.7) -> Transcorp Ughelli'
b = m(nm, 'bundle')
OUT['y2_rows'] = b['rows'][:4] + [b['rows'][-1]]

# 8. «В 9-14 раз дешевле»: сравнение в одних единицах — только компрессорный пакет 3x50%
Gs = B.G
hpq = Gs.hp_per_mmscfd(2.0, 50)
pkg = hpq*1.5*1250/1e6
OUT['op_cost_cmp'] = dict(pkg_per_mmscfd=round(pkg, 2), pkg_x2_6=round(pkg*2.6, 2), seplat=(0.22, 0.36),
                          ratio_pkg=(round(pkg/0.36, 1), round(pkg/0.22, 1)), ratio_x2_6=(round(pkg*2.6/0.36, 1), round(pkg*2.6/0.22, 1)))

json.dump(OUT, open(os.path.join(HERE, 'sk_money_out.json'), 'w'), indent=1, ensure_ascii=False)
for k, v in OUT.items():
    print('\n===', k)
    if isinstance(v, dict):
        for kk, vv in v.items(): print(' ', kk, vv)
    else: print(' ', v)
