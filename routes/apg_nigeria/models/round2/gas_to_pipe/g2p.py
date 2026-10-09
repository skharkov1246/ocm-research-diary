# Раунд 2. Вариант «газ в трубу»: компрессия + осушка/точка росы + отвод до газовой ТЭС или в сеть.
# Метки: [П] первичный, [В] вторичный, [О] расчёт, [Д] допущение с диапазоном.
# Поправки аудита NaCN: капитал x2.4/2.6/2.8 к пакету оборудования; денежный опекс >= скрин x1.7/2.0/2.4;
# персонал с экспатом; ТОиР 3.5-5%; выход по нижней границе; спрос /1.3-2.7 до контракта.
# Второй год: стройка, разгон, загрузка факела 0.78 (VIIRS), падение дебита ~8%/г, оборотный капитал,
# девальвация на дебиторке (оплата в найре), налог 34%.
# Запуск: python3 g2p.py  -> g2p_out.json рядом.
import json, math, os, csv, collections
HERE = os.path.dirname(os.path.abspath(__file__))
TAX = 0.34                                   # CIT 30% + development levy 4% (NTA 2025) [В]
YEARS = 10                                   # лет эксплуатации [Д]

def npv(r, cf): return sum(c/(1+r)**t for t, c in enumerate(cf))
def irr(cf):
    lo, hi = -0.99, 3.0
    if npv(lo, cf)*npv(hi, cf) > 0: return None
    for _ in range(200):
        m = (lo+hi)/2
        if npv(lo, cf)*npv(m, cf) <= 0: hi = m
        else: lo = m
    return m

def hp_per_mmscfd(p_suc, p_dis, stages=3):
    # GPSA-эмпирика: BHP = 22 * R * N * F * Q, F=1.10 при 3 ступенях [В: GPSA Data Book, гл.13] -> [О]
    R = (p_dis/p_suc)**(1/stages); F = {1: 1.0, 2: 1.08, 3: 1.10, 4: 1.12}[stages]
    return 22*R*stages*F

SC = {
 # opt / base / pess
 'opt':  dict(lf=0.85, decl=0.05, up=0.97, own=0.97, p_suc=3.0, p_dis=40, usd_hp=1000, cm=2.4, dehy=0.30, meter=0.30,
              pipe_km=0.5, route=1.2, dev=0.6, build=1, ramp=0.8,
              gas_in=0.25, hv=1.10, price=2.68, transport=0.0, pdrift=0.02, coll=0.97, dso=45, fxloss=0.0,
              staff=0.45, sec=0.25, comm=0.10, tor=0.035, ins=0.010, ovh=0.15, opx_screen_mult=1.7, demand=1/1.3, top_all=False),
 'base': dict(lf=0.78, decl=0.08, up=0.92, own=0.95, p_suc=2.0, p_dis=50, usd_hp=1250, cm=2.6, dehy=0.35, meter=0.40,
              pipe_km=0.8, route=1.3, dev=0.8, build=1, ramp=0.7,
              gas_in=0.60*1.05+0.15, hv=1.05, price=2.18, transport=0.0, pdrift=0.0, coll=0.85, dso=120, fxloss=0.12,
              staff=0.55, sec=0.35, comm=0.15, tor=0.0425, ins=0.015, ovh=0.20, opx_screen_mult=2.0, demand=1/2.0, top_all=False),
 'pess': dict(lf=0.60, decl=0.12, up=0.85, own=0.92, p_suc=1.5, p_dis=70, usd_hp=1500, cm=2.8, dehy=0.45, meter=0.50,
              pipe_km=1.2, route=1.5, dev=1.0, build=2, ramp=0.5,
              gas_in=1.00*1.02+0.30, hv=1.02, price=2.18, transport=0.30, pdrift=-0.01, coll=0.60, dso=240, fxloss=0.16,
              staff=0.65, sec=0.50, comm=0.30, tor=0.05, ins=0.020, ovh=0.30, opx_screen_mult=2.4, demand=1/2.7, top_all=True),
}
# Оператор месторождения (NEPL/Seplat) делает то же сам: персонал, охрана, площадка и газовый завод уже есть.
# Капитал x1.0-1.5 к пакету (ориентир Seplat: $10.8 млн на компрессию Amukpe/Oben/Sapele под снижение
# факела ~30 MMscf/д в 2023 + 20 в 2024 [В] -> $0.22-0.36 млн на MMscf/д [О]); газ свой, штраф избегается.
SC['operator'] = dict(SC['base'], usd_hp=1250, cm=1.3, dev=0.2, staff=0.10, sec=0.05, comm=0.05, ovh=0.05, tor=0.035, ins=0.01,
                      opx_screen_mult=1.0, gas_in=0.0, demand=1.0)
SHRINK = 0.055   # топливо газопоршневых приводов компрессоров ~8000 Btu/hp-ч [Д] -> ~5-6% входящего газа [О]; +1% потери

def staff_scale(Q):
    return 1.0 if Q <= 1 else (1.3 if Q <= 7 else 1.6)

def capex(Q, km_straight, p, ov=None):
    p = dict(p); p.update(ov or {})
    hpq = hp_per_mmscfd(p['p_suc'], p['p_dis'])
    hp_inst = hpq*Q*1.5                                  # 3x50% (N+1), размер под замер 2025 [Д]
    pkg = dict(compr=hp_inst*p['usd_hp']/1e6, dehy=p['dehy']*Q**0.6,          # TEG + точка росы по УВ, пакет [Д]
               meter=p['meter']+0.03*Q)
    pkg_sum = sum(pkg.values())
    pipe_len = km_straight*p['route']
    pipe = pipe_len*p['pipe_km']*(1 + 0.15*(Q > 7))       # 8"→10" на 15 MMscf/д [Д]
    inst = pkg_sum*p['cm'] + pipe + p['dev'] + 0.03       # + лицензия NMDPRA ~$22-42k [В]
    idc = inst*0.14*p['build']/2                          # проценты на стройке 14% [Д]
    screen = pkg_sum + pipe_len*0.4 + 0.3                 # промоутерский скрин: пакет x1, труба $0.4M/км [Д]
    return dict(hp_per=round(hpq), hp=round(hp_inst), pkg=round(pkg_sum, 2), pipe_len=round(pipe_len, 1), pipe=round(pipe, 2),
                total=round(inst+idc, 2), screen=round(screen, 2))

def run(Q, km_straight, sc='base', ov=None, contract=False, edti=False, years=YEARS):
    p = dict(SC[sc]); p.update(ov or {})
    cx = capex(Q, km_straight, p); C = cx['total']
    dem = 1.0 if contract else p['demand']
    hpq = cx['hp_per']
    # денежный опекс: снизу-вверх и проверка против скрина x1.7-2.4
    ss = staff_scale(Q)
    fixed_bu = (p['staff']+p['ovh'])*ss + p['sec'] + 0.01*cx['pipe_len'] + p['comm'] + C*(p['tor']+p['ins'])
    fixed_screen = 0.20*ss + 0.10 + cx['screen']*(0.02+0.005)       # скрин: ТОиР 2%, страховка 0.5% [Д]
    fixed = max(fixed_bu, fixed_screen*p['opx_screen_mult'])
    b = p['build']
    cf = [-C/b]*b                       # стройка: капитал равными долями в t=0..b-1
    prev_rev = 0.0; rows = []
    for t in range(1, years+1):
        Gin = Q*p['lf']*(1-p['decl'])**(t-1)*p['up']*p['own']*(p['ramp'] if t == 1 else 1.0)   # MMscf/д на входе
        take = Gin*dem
        sales = take*(1-SHRINK-0.01)                                                           # MMscf/д товарного газа
        price = (p['price']*(1+p['pdrift'])**(t-1) - p['transport'])*p['hv']                    # $/Mscf
        rev_inv = sales*1000*365*price/1e6
        rev = rev_inv*p['coll']*(1 - p['fxloss']*p['dso']/365)
        gas_paid = (Gin if p['top_all'] else take)*1000*365*p['gas_in']/1e6
        opex = fixed*(1.03)**(t-1) + gas_paid
        e = rev - opex
        dep = C/10
        tx = max(0.0, (e-dep)*TAX)
        ed = (min(tx, 0.05*C) if t <= 5 else 0.0) if edti else 0.0
        wc = p['dso']/365*(rev-prev_rev); prev_rev = rev
        c = e - tx + ed - wc + (p['dso']/365*rev if t == years else 0)
        cf.append(c)
        rows.append(dict(t=t, Gin=round(Gin, 2), sales=round(sales, 2), rev=round(rev, 2), opex=round(opex, 2), ebitda=round(e, 2)))
    i = irr(cf)
    acc = 0; pb = None
    for k, c in enumerate(cf):
        if pb is None and acc+c >= 0 and k > 0 and acc < 0: pb = round(k-1+(-acc)/c - (b-1), 1)
        acc += c
    return dict(Q=Q, km=km_straight, sc=sc, contract=contract, capex=cx, fixed_opex=round(fixed, 2),
                opex_vs_screen=round(fixed/fixed_screen, 2), y1=rows[0], y2=rows[1], y5=rows[4],
                npv15=round(npv(.15, cf), 2), npv20=round(npv(.20, cf), 2), irr=None if i is None else round(i*100, 1), payback=pb)

def solve(fn, lo, hi, it=80):
    flo, fhi = fn(lo), fn(hi)
    if (flo > 0) == (fhi > 0): return None            # корня в диапазоне нет
    for _ in range(it):
        m = (lo+hi)/2; fm = fn(m)
        if (fm > 0) == (flo > 0): lo, flo = m, fm
        else: hi = m
    return round((lo+hi)/2, 2)

NAMED = {   # (VIIRS 2025, минимум 2023-25, км по прямой до ТЭС) [О: wb/nigeria_flares.csv; координаты ТЭС — GEM wiki [В]]
  'Ughelli East (NEPL) -> Transcorp Ughelli': (6.7, 6.3, 3.7),
  'Oredo (NEPL) -> Sapele PS':                 (15.5, 8.5, 16.4),
  'Afam Umuosi (NEPL) -> Afam VI/IV-V':        (5.6, 2.9, 0.9),
  'Sapele (Seplat) -> Sapele PS':              (6.0, 6.0, 7.4),
  'Oben (Seplat) -> Oben GP/ELPS [Д 2 км]':    (8.7, 8.7, 2.0),
  'Utorogu (NEPL) -> Transcorp Ughelli':       (4.3, 4.3, 10.9),
}

if __name__ == '__main__':
    out = {'generic': {}, 'named': {}, 'thresholds': {}, 'sens': {}, 'operator': {}}
    for Q in (1, 5, 15):
        for km in (1, 4, 8, 16):
            for sc in ('opt', 'base', 'pess'):
                for con in (False, True):
                    out['generic'][f'Q{Q}_{km}km_{sc}_{"GSA" if con else "pre"}'] = run(Q, km, sc, contract=con)
    for name, (q25, qmin, km) in NAMED.items():
        for qn, q in (('q2025', q25), ('qmin', qmin)):
            for sc in ('opt', 'base', 'pess'):
                for con in (False, True):
                    out['named'][f'{name}|{qn}|{sc}|{"GSA" if con else "pre"}'] = run(q, km, sc, contract=con)
    # пороги NPV20 = 0, база
    for label, Q, km in (('Q1_4km', 1, 4), ('Q5_4km', 5, 4), ('Q5_8km', 5, 8), ('Q15_8km', 15, 8),
                         ('UghelliEast', 6.7, 3.7), ('AfamUmuosi', 5.6, 0.9), ('Oredo_qmin', 8.5, 16.4)):
        for sc in ('base', 'opt'):
            for con in (True, False):
                k = f'{label}_{sc}_{"GSA" if con else "pre"}'
                out['thresholds'][k] = dict(
                  price_usd_mmbtu=solve(lambda x: run(Q, km, sc, ov=dict(price=x), contract=con)['npv20'], 0.5, 40),
                  fee_from_operator_usd_mscf=(lambda g: None if g is None else round(-g, 2))(
                      solve(lambda x: run(Q, km, sc, ov=dict(gas_in=x), contract=con)['npv20'], -15, 3)),
                  Q_mmscfd=solve(lambda x: run(x, km, sc, contract=con)['npv20'], 0.2, 300),
                  collection=solve(lambda x: run(Q, km, sc, ov=dict(coll=x), contract=con)['npv20'], 0.05, 1.0))
    # оператор делает сам: газ свой, штраф избегается. Избегаемый штраф не вычитается из базы (PIA s.104(3) [П]),
    # поэтому до налога он эквивалентен штраф/(1-0.34). Ставки: закон $3.50/Mscf [В], фактический сбор $1.6-1.9 [О фаз 1-2].
    for name, (q25, qmin, km) in NAMED.items():
        for pen in (0.0, 1.75, 3.50):
            out['operator'][f'{name}|pen{pen}'] = run(q25, km, 'operator', ov=dict(gas_in=-pen/(1-TAX)), contract=True)
    # грубая оценка критика: 4.25 MMscf/д, маржа 1.43/1.93 $/Mscf, 10 лет, 15%, без опекса
    out['critic_repro_PV15'] = {k: round(4.25*365*1000*v/1e6*(1-1.15**-10)/0.15, 1) for k, v in (('DBP', 1.43), ('comm', 1.93))}
    for nm, ov, kw in (('base', {}, {}), ('EDTI', {}, dict(edti=True)), ('coll_0.6', dict(coll=0.6), {}),
                       ('escrow_coll1_dso45_nofx', dict(coll=1.0, fxloss=0, dso=45), {}),
                       ('price_2.68', dict(price=2.68), {}), ('pipe_x2', dict(pipe_km=1.6), {}),
                       ('cm_1.0', dict(cm=1.0), {}), ('build_2y', dict(build=2, ramp=0.5), {}), ('decl_15', dict(decl=0.15), {}),
                       ('opex_screen_x1.7_only', dict(staff=0.2, sec=0.1, comm=0.05, ovh=0.05, tor=0.035, ins=0.01, opx_screen_mult=1.7), {}),
                       ('fee_from_operator_1.0', dict(gas_in=-1.0), {}), ('fee_from_operator_1.75', dict(gas_in=-1.75), {})):
        r = run(5, 4, 'base', ov=ov, contract=True, **kw)
        out['sens'][nm] = dict(capex=r['capex']['total'], fixed=r['fixed_opex'], ebitda_y2=r['y2']['ebitda'], npv15=r['npv15'], npv20=r['npv20'], irr=r['irr'])
    json.dump(out, open(os.path.join(HERE, 'g2p_out.json'), 'w'), indent=1, ensure_ascii=False)

    def line(k, r): print(f"{k:58s} cap {r['capex']['total']:6.1f} (pkg {r['capex']['pkg']:.2f} pipe {r['capex']['pipe']:.1f} hp {r['capex']['hp']}) "
                          f"F {r['fixed_opex']:.2f}(x{r['opex_vs_screen']}) E1 {r['y1']['ebitda']:6.2f} E2 {r['y2']['ebitda']:6.2f} "
                          f"NPV15 {r['npv15']:6.1f} NPV20 {r['npv20']:6.1f} IRR {r['irr']} PB {r['payback']}")
    for sec in ('generic', 'named', 'operator'):
        print('==', sec)
        for k, r in out[sec].items():
            if sec != 'generic' or ('_4km_' in k or '_1km_' in k or '_8km_' in k): line(k, r)
    print('== thresholds'); [print(k, v) for k, v in out['thresholds'].items()]
    print('== critic', out['critic_repro_PV15'])
    print('== sens Q5 4km base GSA'); [print(k, v) for k, v in out['sens'].items()]
