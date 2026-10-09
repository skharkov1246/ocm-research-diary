# Раунд 2. Вариант «компрессия / обратная закачка / газлифт как услуга с оплатой оператором (BOOT)».
# Мы (подрядчик) строим и эксплуатируем компрессию у факела, оператор (NEPL/Seplat/Oando/Renaissance)
# платит тариф $/Mscf в долларах. Газ остаётся у оператора: он идёт (а) на вход его газового завода/в экспорт,
# (б) в нагнетательную скважину оператора, (в) в систему газлифта.
# Метки: [П] первичный, [В] вторичный, [О] расчёт, [Д] допущение с диапазоном.
# Поправки аудита NaCN: капитал x2.4/2.6/2.8 к пакету; денежный опекс >= скрин x1.7/2.0/2.4; персонал с экспатом;
# ТОиР 3.5-5%; выход по нижней границе; спрос /1.3-2.7 до контракта. Второй год: стройка, разгон, загрузка 0.78,
# падение ~8%/г, оборотный капитал, налог 34%, девальвация 8-16%/г при тарифе в найре (вариант ngn).
# Запуск: python3 boot.py -> boot_out.json и печать рядом.
import json, os
HERE = os.path.dirname(os.path.abspath(__file__))
TAX = 0.34          # CIT 30% + development levy 4% (NTA 2025) [В]
NCDMB = 0.01        # 1% стоимости upstream-контракта в NCDF (NOGICD Act s.104(2)) [В]
YEARS = 10          # срок BOOT, затем передача оператору за $0 [Д]

def npv(r, cf): return sum(c/(1+r)**t for t, c in enumerate(cf))
def irr(cf):
    lo, hi = -0.99, 3.0
    if npv(lo, cf)*npv(hi, cf) > 0: return None
    for _ in range(200):
        m = (lo+hi)/2
        if npv(lo, cf)*npv(m, cf) <= 0: hi = m
        else: lo = m
    return m

def hp_per_mmscfd(p_suc, p_dis, stages):
    # GPSA: BHP = 22 * R * N * F * Q (Q в MMscf/д), F по числу ступеней [В: GPSA Data Book, гл.13] -> [О]
    R = (p_dis/p_suc)**(1/stages); F = {1: 1.0, 2: 1.08, 3: 1.10, 4: 1.12, 5: 1.13}[stages]
    return 22*R*stages*F

# Конфигурации: давление нагнетания (бар а) [Д], число ступеней, $/л.с. пакета [Д], отвод (км) [Д]
CONF = {
 'boost':   dict(p_dis=(30, 50, 70),    stages=3, usd_hp=(1000, 1250, 1500), km=(0.5, 1.0, 2.0)),   # на вход завода / в экспорт
 'reinj':   dict(p_dis=(150, 200, 250), stages=4, usd_hp=(1400, 1800, 2300), km=(1.0, 2.0, 3.0)),   # в скважину оператора
 'gaslift': dict(p_dis=(70, 90, 110),   stages=4, usd_hp=(1200, 1500, 1800), km=(1.0, 2.0, 3.0)),   # в газлифтный коллектор
}
SC = {
 'opt':  dict(i=0, lf=0.85, decl=0.05, up=0.97, p_suc=3.0, cm=2.4, pipe_km=0.5, dev=0.4, build=1, ramp=0.8,
              fee_drift=0.02, coll=0.98, dso=45, staff=0.45, sec=0.25, ovh=0.15, tor=0.035, ins=0.010, var=0.02,
              mult=1.7, demand=1/1.3),
 'base': dict(i=1, lf=0.78, decl=0.08, up=0.92, p_suc=2.0, cm=2.6, pipe_km=0.8, dev=0.6, build=1, ramp=0.7,
              fee_drift=0.0, coll=0.92, dso=90, staff=0.55, sec=0.35, ovh=0.20, tor=0.0425, ins=0.015, var=0.03,
              mult=2.0, demand=1/2.0),
 'pess': dict(i=2, lf=0.60, decl=0.12, up=0.85, p_suc=1.5, cm=2.8, pipe_km=1.2, dev=0.8, build=2, ramp=0.5,
              fee_drift=0.0, coll=0.80, dso=180, staff=0.65, sec=0.50, ovh=0.30, tor=0.050, ins=0.020, var=0.04,
              mult=2.4, demand=1/2.7),
}
# Оператор делает сам: площадка, персонал, охрана уже есть; капитал x1.3 к пакету (ориентир Seplat $10.8 млн
# на компрессию Amukpe/Oben/Sapele под снижение факела ~30+20 MMscf/д [В] -> $0.22 млн на MMscf/д [О]).
OPER = dict(cm=1.3, dev=0.15, staff=0.10, sec=0.05, ovh=0.05, tor=0.035, ins=0.01, mult=1.0, coll=1.0, dso=0, demand=1.0)

def staff_scale(Q): return 1.0 if Q <= 1 else (1.3 if Q <= 7 else 1.6)

def capex(Q, conf, p, km=None):
    c = CONF[conf]; i = p['i']
    hpq = hp_per_mmscfd(p['p_suc'], c['p_dis'][i], c['stages'])
    hp = hpq*Q*1.5                                       # 3x50% (N+1), размер под VIIRS 2025 [Д]
    usd_hp = p.get('usd_hp', c['usd_hp'][i])
    pkg = hp*usd_hp/1e6 + 0.25*Q**0.6 + 0.30 + 0.03*Q    # компрессоры + сепарация/осушка + учёт [Д]
    L = c['km'][i] if km is None else km
    pipe = L*p['pipe_km']*(1.3 if conf == 'reinj' else 1.0)   # HP-линия дороже [Д]
    inst = pkg*p['cm'] + pipe + p['dev']
    idc = inst*0.14*p['build']/2                          # проценты на стройке 14% [Д]
    screen = pkg + L*0.4 + 0.2
    return dict(hp_per=round(hpq), hp=round(hp), pkg=round(pkg, 2), pipe=round(pipe, 2), total=round(inst+idc, 2),
                screen=round(screen, 2))

def run(Q, conf='boost', sc='base', fee=1.375, contract=True, ov=None, km=None, ngn=None, top=None, years=YEARS,
        operator=False, rate_tax=TAX):
    p = dict(SC[sc]); p.update(OPER if operator else {}); p.update(ov or {})
    cx = capex(Q, conf, p, km); C = cx['total']
    dem = 1.0 if (contract or operator) else p['demand']
    ss = staff_scale(Q)
    fixed_bu = (p['staff']+p['ovh'])*ss + p['sec'] + C*(p['tor']+p['ins'])
    fixed_screen = 0.20*ss + 0.10 + cx['screen']*0.025
    fixed = max(fixed_bu, fixed_screen*p['mult'])
    b = p['build']; cf = [-C/b]*b; prev = 0.0; rows = []
    G1 = None
    for t in range(1, years+1):
        Gin = Q*p['lf']*(1-p['decl'])**(t-1)*p['up']*(p['ramp'] if t == 1 else 1.0)*dem   # MMscf/д
        if G1 is None: G1 = Q*p['lf']*p['up']*dem
        billed = max(Gin, top*Q*(1-p['decl'])**(t-1)) if top else Gin   # ToP: доля объёма VIIRS 2025 с падением дебита [Д]
        f = fee*(1+p['fee_drift'])**(t-1)
        if ngn is not None: f = f/(1+ngn)**(t-1)            # тариф в найре, фиксирован при подписании [Д]
        rev = billed*365*f/1000*p['coll']*(1-NCDMB) if not operator else 0.0   # $ млн
        opex = fixed*1.03**(t-1) + Gin*365*p['var']/1000   # var $/Mscf -> $ млн
        e = rev - opex
        tx = max(0.0, (e - C/10)*rate_tax)
        wc = p['dso']/365*(rev-prev); prev = rev
        cf.append(e - tx - wc + (p['dso']/365*rev if t == years else 0))
        rows.append(dict(t=t, G=round(Gin, 2), rev=round(rev, 2), opex=round(opex, 2), ebitda=round(e, 2)))
    i = irr(cf); acc = 0; pb = None
    for k, c in enumerate(cf):
        if pb is None and k > 0 and acc < 0 and acc+c >= 0: pb = round(k-1+(-acc)/c-(b-1), 1)
        acc += c
    return dict(Q=Q, conf=conf, sc=sc, fee=fee, capex=cx, fixed=round(fixed, 2), x_screen=round(fixed/fixed_screen, 2),
                y1=rows[0], y2=rows[1], npv15=round(npv(.15, cf), 2), npv20=round(npv(.20, cf), 2),
                irr=None if i is None else round(i*100, 1), payback=pb, cf=[round(x, 2) for x in cf])

def solve(fn, lo, hi, it=70):
    flo, fhi = fn(lo), fn(hi)
    if (flo > 0) == (fhi > 0): return None
    for _ in range(it):
        m = (lo+hi)/2; fm = fn(m)
        if (fm > 0) == (flo > 0): lo, flo = m, fm
        else: hi = m
    return round((lo+hi)/2, 3)

def levelized_operator(Q, conf, sc, r):
    """Себестоимость $/Mscf для оператора, делающего сам (до налога: и тариф, и свой капитал вычитаемы)."""
    o = run(Q, conf, sc, operator=True, rate_tax=0.0)
    p = dict(SC[sc])
    # объёмы по годам
    G = []
    for t in range(1, YEARS+1):
        G.append(Q*p['lf']*(1-p['decl'])**(t-1)*p['up']*(p['ramp'] if t == 1 else 1.0))
    pvv = sum(g*365/1000/(1+r)**(t+p['build']-1) for t, g in enumerate(G, 1))   # млн Mscf
    pvc = -npv(r, o['cf'])
    return round(pvc/pvv, 3)

def contract_compression_benchmark(Q, conf, sc, prem):
    """Рыночная аренда компрессии: USAC 2025 $21.38/л.с.-мес [П] x международная надбавка [Д], $/Mscf."""
    p = SC[sc]; cx = capex(Q, conf, p)
    hp_run = cx['hp']/1.5*1.0                     # оплачиваемая мощность = рабочая (без резерва) [Д]
    G_avg = sum(Q*p['lf']*(1-p['decl'])**(t-1)*p['up'] for t in range(1, YEARS+1))/YEARS
    return round(hp_run*21.38*prem*12/(G_avg*365*1000), 3)

NAMED = {  # VIIRS 2025, минимум 2023-25 [О: wb/nigeria_flares.csv]; конфигурация и отвод по наличию газового узла [Д]
 'Ughelli East (NEPL)':      (6.7, 6.3, 'boost', 1.0),
 'Oredo (NEPL, IGHF)':       (15.5, 8.5, 'boost', 1.0),
 'Opuama (NEPL, болото)':    (11.6, 9.4, 'reinj', 2.0),
 'Odidi (NEPL)':             (20.6, 5.0, 'boost', 1.0),
 'Obiafu-Obrikom (Oando)':   (18.0, 10.4, 'boost', 1.0),
 'Kwale (Oando)':            (10.0, 9.9, 'boost', 1.0),
 'Saghara (Renaissance)':    (11.7, 4.3, 'reinj', 2.0),
 'Oben (Seplat)':            (8.7, 8.7, 'boost', 0.5),
 'Afam Umuosi (NEPL)':       (5.6, 2.9, 'boost', 1.0),
}
FEES = (1.0, 1.375, 1.75)

if __name__ == '__main__':
    out = dict(generic={}, named={}, thresholds={}, operator={}, benchmark={}, sens={}, wtp={})
    for conf in CONF:
        for Q in (1, 5, 15):
            for sc in SC:
                for fee in FEES:
                    for con in (True, False):
                        r = run(Q, conf, sc, fee, contract=con); r.pop('cf')
                        out['generic'][f'{conf}|Q{Q}|{sc}|fee{fee}|{"contract" if con else "pre"}'] = r
            for sc in SC:
                out['thresholds'][f'{conf}|Q{Q}|{sc}'] = dict(
                    fee_npv20_contract=solve(lambda x: run(Q, conf, sc, x)['npv20'], 0.01, 30),
                    fee_npv15_contract=solve(lambda x: run(Q, conf, sc, x)['npv15'], 0.01, 30),
                    fee_npv20_pre=solve(lambda x: run(Q, conf, sc, x, contract=False)['npv20'], 0.01, 60),
                    Q_npv20_fee1375=solve(lambda x: run(x, conf, sc, 1.375)['npv20'], 0.1, 300),
                    Q_npv20_fee175=solve(lambda x: run(x, conf, sc, 1.75)['npv20'], 0.1, 300),
                    cm_npv20_fee175=solve(lambda x: run(Q, conf, sc, 1.75, ov=dict(cm=x))['npv20'], 0.5, 8))
                out['operator'][f'{conf}|Q{Q}|{sc}'] = dict(
                    lcoc15=levelized_operator(Q, conf, sc, .15), lcoc20=levelized_operator(Q, conf, sc, .20),
                    capex_oper=run(Q, conf, sc, operator=True)['capex']['total'])
                out['benchmark'][f'{conf}|Q{Q}|{sc}'] = {f'prem{k}': contract_compression_benchmark(Q, conf, sc, k)
                                                          for k in (1.0, 1.5, 2.5)}
    for name, (q, qmin, conf, km) in NAMED.items():
        for qn, qq in (('q2025', q), ('qmin', qmin)):
            for sc in SC:
                for fee in FEES:
                    r = run(qq, conf, sc, fee, km=km, ov=None) if qn == 'q2025' else \
                        run(q, conf, sc, fee, km=km, ov=dict(lf=SC[sc]['lf']*qmin/q))   # размер по 2025, газ по минимуму
                    r.pop('cf'); out['named'][f'{name}|{qn}|{sc}|fee{fee}'] = r
            out['thresholds'][f'NAMED {name}|{qn}|base'] = dict(fee_npv20_contract=solve(
                (lambda x: run(qq, conf, 'base', x, km=km)['npv20']) if qn == 'q2025' else
                (lambda x: run(q, conf, 'base', x, km=km, ov=dict(lf=0.78*qmin/q))['npv20']), 0.01, 30))
    # чувствительности: 5 MMscf/д, boost, base, тариф 1.375, контракт
    for nm, kw in (('base', {}), ('ngn_8', dict(ngn=0.08)), ('ngn_16', dict(ngn=0.16)),
                   ('ToP70', dict(top=0.7)), ('coll_0.6_dso240', dict(ov=dict(coll=0.6, dso=240))),
                   ('build_2y', dict(ov=dict(build=2, ramp=0.5))), ('decl_15', dict(ov=dict(decl=0.15))),
                   ('cm_1.6', dict(ov=dict(cm=1.6))), ('no_expat', dict(ov=dict(staff=0.25, ovh=0.10))),
                   ('lf_0.6', dict(ov=dict(lf=0.6)))):
        for conf in ('boost', 'reinj'):
            r = run(5, conf, 'base', 1.375, **kw)
            out['sens'][f'{conf}|{nm}'] = dict(capex=r['capex']['total'], fixed=r['fixed'], e1=r['y1']['ebitda'],
                                               e2=r['y2']['ebitda'], npv15=r['npv15'], npv20=r['npv20'], irr=r['irr'])
    # готовность оператора платить, $/Mscf до налога: штраф не вычитается (PIA s.104(3) [П]), тариф вычитается.
    # Налог оператора: CIT 30% + HT 30% (PIA, конвертированные PML) ~51% [В/О]; PPT 85% (неконвертированные JV) [В].
    for pen, lab in ((1.6, 'fact_lo'), (1.9, 'fact_hi'), (3.5, 'statute')):
        out['wtp'][lab] = {f't{t}': round(pen/(1-t), 2) for t in (0.34, 0.51, 0.85)}
    # «подрядчик уровня Enerflex»: капитал x1.6, без экспата, база; и ориентир Seplat ($0.22 млн/MMscf/д) как x~0.6
    out['extra'] = {}
    for Q in (5, 15):
        for nm, ov in (('cm1.6_noexpat', dict(cm=1.6, staff=0.25, ovh=0.10)),
                       ('cm1.6_noexpat_ToP70', dict(cm=1.6, staff=0.25, ovh=0.10)),
                       ('seplat_like_cm0.6_noexpat', dict(cm=0.6, staff=0.25, ovh=0.10, dev=0.2))):
            kw = dict(top=0.7) if 'ToP' in nm else {}
            r = run(Q, 'boost', 'base', 1.75, ov=ov, **kw)
            out['extra'][f'boost|Q{Q}|{nm}'] = dict(capex=r['capex']['total'], e1=r['y1']['ebitda'], e2=r['y2']['ebitda'],
                npv15=r['npv15'], npv20=r['npv20'], irr=r['irr'],
                fee_npv20=solve(lambda x: run(Q, 'boost', 'base', x, ov=ov, **kw)['npv20'], 0.01, 30))
    # лучший случай: 15 MMscf/д, opt, тариф 1.75, контракт + ToP 70%
    r = run(15, 'boost', 'opt', 1.75, top=0.7); r.pop('cf'); out['extra']['best_Q15_opt_fee1.75_ToP70'] = r
    json.dump(out, open(os.path.join(HERE, 'boot_out.json'), 'w'), indent=1, ensure_ascii=False)

    def ln(k, r): print(f"{k:52s} cap {r['capex']['total']:6.2f} (pkg {r['capex']['pkg']:.2f} hp {r['capex']['hp']}) F {r['fixed']:.2f}(x{r['x_screen']}) "
                        f"E1 {r['y1']['ebitda']:6.2f} E2 {r['y2']['ebitda']:6.2f} NPV15 {r['npv15']:6.2f} NPV20 {r['npv20']:6.2f} IRR {r['irr']} PB {r['payback']}")
    print('== generic (contract)')
    for k, r in out['generic'].items():
        if 'contract' in k: ln(k, r)
    print('== generic pre-contract, base')
    for k, r in out['generic'].items():
        if '|pre' in k and '|base|' in k: ln(k, r)
    print('== named'); [ln(k, r) for k, r in out['named'].items()]
    for s in ('thresholds', 'operator', 'benchmark', 'sens', 'wtp', 'extra'):
        print('==', s); [print(' ', k, v) for k, v in out[s].items()]
